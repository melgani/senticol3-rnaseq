# RNA-seq analysis with limma/voom
# ============================================================

# Libraries
library(DESeq2)
library(limma)
library(edgeR)
library(FactoMineR)
library(factoextra)
library(pheatmap)
library(RColorBrewer)
library(corrplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(MCPcounter)
library(mMCPcounter)
library(devtools)


# ============================================================
# 1. Data loading
# ============================================================

# Define path to data
countsdatapath <- "path/to/data"  # TODO: update path

# Raw count tables
D1308 <- as.matrix(read.csv(file.path(countsdatapath, "tablecounts_raw_D1308.csv"), row.names = 1))
D1315 <- as.matrix(read.csv(file.path(countsdatapath, "tablecounts_raw_D1315.csv"), row.names = 1))
d <- cbind(D1308, D1315)

# TPM normalized count tables
D1308.tpm <- as.matrix(read.csv(file.path(countsdatapath, "tablecounts_tpm_D1308.csv"), row.names = 1))
D1315.tpm <- as.matrix(read.csv(file.path(countsdatapath, "tablecounts_tpm_D1315.csv"), row.names = 1))
d.tpm <- cbind(D1308.tpm, D1315.tpm)

# Sample plan
splan <- read.csv(file.path(countsdatapath, "sampledata.csv"), row.names = 1, header = TRUE)

# Update column names
colnames(d)     <- as.character(splan[colnames(d), "sname"])
colnames(d.tpm) <- as.character(splan[colnames(d.tpm), "sname"])


# ============================================================
# 2. Restrict to coding genes
# ============================================================

d.annot <- read.csv("tableannot.csv", row.names = 1)
rownames(d.annot) <- d.annot[, 1]
coding.annot <- read.table("gencode.v19.annotation_proteinCoding_gene.bed")
coding.ids   <- d.annot$gene_id[which(d.annot$gene_name %in% coding.annot[, 4])]
d.coding     <- d[coding.ids, ]
d.tpm.coding <- d.tpm[coding.ids, ]


# ============================================================
# 3. Exploratory data mining
# ============================================================

dim(d.coding)
colSums(d.coding)

# Genes with zero counts across all samples
nbgenes_at_zeros <- length(which(rowSums(d.coding) == 0))

# Genes expressed (>1 count) per sample
number_expressed <- function(x, mincounts = 1) length(which(x > mincounts))
nbgenes_per_sample <- apply(d.coding, 2, number_expressed)

# Distribution of raw counts (log2)
boxplot(log2(1 + d),
        las = 2, ylab = "raw counts (log2)", col = "gray50",
        pch = 16, cex = 0.3, cex.axis = 0.70, cex.lab = 1)

# Genes expressed (TPM >= 1) per sample
nb_expressed_genes <- apply(d.tpm.coding, 2, function(x) length(which(x >= 1)))
barplot(nb_expressed_genes,
        las = 2, ylab = "expressed genes (TPM>=1)", col = "gray50",
        pch = 16, cex = 0.55, cex.axis = 0.6, cex.lab = 1, ylim = c(0, 20000))


# ============================================================
# 4. Exploratory analysis
# ============================================================

# --- Variance stabilization with DESeq2 ---
dds <- DESeqDataSetFromMatrix(
  countData = d.coding,
  colData   = DataFrame(condition = splan$Histological_type_to_keep),
  design    = ~condition
)
dds <- estimateSizeFactors(dds)
dds <- dds[rowSums(counts(dds)) > 0, ]
rld <- rlog(dds, blind = TRUE)

# Comparison of normalization methods
par(mfrow = c(1, 3))
plot(counts(dds, normalized = TRUE)[, 3:4],
     pch = 16, cex = 0.3, xlim = c(0, 20e3), ylim = c(0, 20e3), main = "normalized counts")
plot(log2(counts(dds, normalized = TRUE)[, 3:4] + 1),
     pch = 16, cex = 0.3, main = "log2 normalized counts")
plot(assay(rld)[, 3:4],
     pch = 16, cex = 0.3, main = "rlog normalized counts")


# --- PCA ---
d.rlog <- assay(rld)
gvar         <- apply(d.rlog, 1, var)
mostvargenes <- order(gvar, decreasing = TRUE)[1:1000]
res_pca      <- PCA(t(d.rlog[mostvargenes, ]), ncp = 3, graph = FALSE)

fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 65))
fviz_pca_ind(res_pca, axes = c(1, 2), col.ind = splan$Event, label = "none",
             legend.title = "Group", mean.point = FALSE, pointshape = 19,
             labelsize = 4, pointsize = 1.5, addEllipses = FALSE, ggtheme = theme_gray())
fviz_pca_ind(res_pca, axes = c(2, 3), col.ind = splan$LVSI_to_keep, label = "none",
             legend.title = "Group", mean.point = FALSE, pointshape = 19,
             labelsize = 4, pointsize = 1.5, addEllipses = FALSE, ggtheme = theme_gray())

# Top contributors to PC1
best.contrib <- names(sort(res_pca$var$contrib[, "Dim.1"], decreasing = TRUE))
par(mfrow = c(1, 3))
for (i in 1:3) {
  barplot(d.tpm[best.contrib[i], ], col = "gray50", border = "white",
          las = 2, main = best.contrib[i])
}

# Association between PCA and clinical variables
d.rlog.mostvargenes <- as.data.frame(t(d.rlog[mostvargenes, ]))
d.rlog.splan        <- cbind(d.rlog.mostvargenes, splan)
res.pca.quali <- PCA(d.rlog.splan, scale.unit = TRUE, ncp = 3,
                     quali.sup = c(1001:1006), graph = FALSE)
v.test <- res.pca.quali$quali.sup$v.test

corrplot(res.pca.quali$quali.sup$v.test[267:270, ], is.corr = FALSE,
         tl.col = "black", cl.ratio = 0.60, method = "circle",
         col = COL2("RdBu", 6), col.lim = c(-3, 3))


# --- Hierarchical clustering ---
sampleDist <- dist(t(d.rlog), method = "euclidean")
hc <- hclust(sampleDist, method = "ward.D2")
plot(hc, cex = 0.7, hang = -1)
rect.hclust(hc, k = 2)
groupes.cah <- cutree(hc, k = 2)
print(sort(groupes.cah))

# Heatmap of the 100 most variable genes
mostvargenes100 <- order(gvar, decreasing = TRUE)[1:100]
annot <- data.frame(splan$Event, splan$Histological_type_to_keep, row.names = splan$sname)
pheatmap(d.rlog[mostvargenes100, ],
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method        = "ward.D2",
         color                    = colorRampPalette(c("cornflowerblue", "white", "coral3"))(50),
         annotation_col           = annot,
         show_rownames            = FALSE,
         scale                    = "row")


# ============================================================
# 4b. Exploratory analysis on a sub-list of DEGs
# ============================================================

d.rlog_1275 <- as.matrix(read.csv(file.path(countsdatapath, "d.rlog_39.csv"), row.names = 1))
colnames(d.rlog_1275) <- as.character(splan[colnames(d.rlog_1275), "sname"])

res_pca.1275 <- PCA(t(d.rlog_1275), ncp = 4, graph = FALSE)
fviz_eig(res_pca.1275, addlabels = TRUE, ylim = c(0, 100))
fviz_pca_ind(res_pca.1275, axes = c(1, 2),
             col.ind = splan$Histological_type_to_keep, label = "none",
             legend.title = "Group", mean.point = FALSE, pointshape = 19,
             labelsize = 4, pointsize = 2, addEllipses = FALSE, ggtheme = theme_gray())
fviz_pca_ind(res_pca.1275, axes = c(2, 3),
             col.ind = splan$Histological_type_to_keep, label = "none",
             legend.title = "Group", mean.point = FALSE, pointshape = 19,
             labelsize = 3, pointsize = 1.5, addEllipses = FALSE, ggtheme = theme_gray())

# Hierarchical clustering on d.rlog_1275
sampleDist    <- dist(t(d.rlog_1275), method = "euclidean")
hc            <- hclust(sampleDist, method = "ward.D2")
plot(hc, cex = 0.7, hang = -1)
rect.hclust(hc, k = 2)
groupes.cah <- cutree(hc, k = 2)
print(sort(groupes.cah))

gvar_1275       <- apply(d.rlog_1275, 1, var)
mostvargenes100 <- order(gvar_1275, decreasing = TRUE)[1:100]
annot <- data.frame(splan$Event, splan$Histological_type_to_keep, row.names = splan$sname)
pheatmap(d.rlog_1275[mostvargenes100, ],
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method        = "ward.D2",
         color                    = colorRampPalette(c("cornflowerblue", "white", "coral3"))(50),
         annotation_col           = annot,
         show_rownames            = FALSE,
         scale                    = "row")


# ============================================================
# 5. Differential analysis (limma/voom)
# ============================================================

# --- Gene filtering: expressed genes (TPM >= 1 in at least one sample) ---
isexpr <- which(apply(d.tpm.coding, 1, function(x) length(which(x >= 1))) >= 1)
d.f    <- d.coding[isexpr, ]

par(mfrow = c(1, 2))
hist(log2(0.1 + rowSums(d.tpm.coding)),
     breaks = 100, xlab = "log2(sum of counts per gene)",
     ylab = "Number of genes", las = 1,
     main = "Gene expression across all samples (TPM)",
     col = "#BBBBFF", ylim = c(0, 3000))
hist(log2(0.1 + rowSums(d.tpm.coding[isexpr, ])),
     breaks = 100, add = TRUE, col = "red", border = "white")


# --- TMM normalization ---
y      <- DGEList(counts = d.f)
y      <- calcNormFactors(y, method = "TMM")
d.norm <- cpm(y, log = TRUE)


# --- Voom + limma ---
design   <- model.matrix(~0 + LVSI_to_keep + Histological_type_to_keep + Event + FIGO_to_keep,
                         data = splan)
v        <- voom(y, design, plot = FALSE)
fit      <- lmFit(v, design)
contrast <- makeContrasts(LVSI_to_keepYes - LVSI_to_keepNo, levels = design)
fit2     <- eBayes(contrasts.fit(fit, contrast))

# Results
res      <- topTable(fit2, number = 1e6, adjust.method = "BH")
hist(res$P.Value, main = "P-value histogram", col = "grey50", border = "white")

# DEGs (FDR < 0.05 and |logFC| > 1)
idx.sign <- which(res$adj.P.Val < 0.05 & abs(res$logFC) > 1)
deg      <- rownames(res[idx.sign, ])


# --- Visualization ---
mygene <- rownames(res)[2]
barplot(d.norm[mygene, ], las = 2, main = mygene, ylab = "Normalized counts",
        col = "grey50", border = "white")

# Heatmap of top DEGs
idx.sub  <- which(splan$Group == "T" | splan$Group == "N")
data.sub <- d.norm[deg[1:5], idx.sub]
pheatmap(data.sub, cutree_cols = 2, show_rownames = FALSE)

# MA plot and volcano plot
plot(res$AveExpr, res$logFC,
     xlab = "A - Mean Expression", ylab = "M - logFC",
     col = ifelse(res$adj.P.Val < 0.05, "red", "black"), pch = 16, cex = 0.5)
volcanoplot(fit2, highlight = 100)


# --- Gene Ontology enrichment ---
entrez.sign <- bitr(deg, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe    <- bitr(rownames(res), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego.BP      <- enrichGO(gene = entrez.sign$ENTREZID, universe = universe$ENTREZID,
                        OrgDb = org.Hs.eg.db, ont = "BP")
dotplot(dropGO(ego.BP, level = c(1:3)), showCategory = 20)


# ============================================================
# 6. Immune deconvolution
# ============================================================

# Install MCPcounter (run once)
# install.packages(c("devtools", "curl"))
# install_github("ebecht/MCPcounter", ref = "master", subdir = "Source", force = TRUE)

# Load TPM data
d.tpm <- as.matrix(read.csv(file.path(countsdatapath, "d.tpm.csv"), row.names = 1))

# Human MCPcounter
res_MCPcounter <- MCPcounter.estimate(
  d.tpm,
  featuresType = "ENSEMBL_ID",
  probesets = read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),
                         sep = "\t", stringsAsFactors = FALSE, colClasses = "character"),
  genes     = read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE,
                         colClasses = "character", check.names = FALSE)
)

pheatmap(res_MCPcounter,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method        = "ward.D2",
         color                    = colorRampPalette(c("cornflowerblue", "white", "coral3"))(50),
         show_rownames            = TRUE,
         scale                    = "row")


# ============================================================
# 7. Custom gene signature scoring
# ============================================================

d.tpm.coding <- as.matrix(read.csv(file.path(countsdatapath, "d.tpm.coding.csv"), row.names = 1))
colnames(d.tpm.coding) <- as.character(splan[colnames(d.tpm.coding), "sname"])

# Load signature gene lists
signCAF <- read.csv(file.path(countsdatapath, "Liste genes.csv"), row.names = 1)
signCAF <- as.list(signCAF)

# Compute mean expression per signature
Sig <- lapply(signCAF, function(g) colMeans(d.tpm.coding[intersect(g, rownames(d.tpm.coding)), ]))
Sig_CAF <- data.frame(
  fibroblasts  = Sig$Fibroblasts,
  endotheliale = Sig$Endothelial.cells,
  tcells       = Sig$T.cells
)

pheatmap(t(Sig_CAF),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method        = "ward.D2",
         show_rownames            = TRUE,
         scale                    = "row",
         angle_col                = 90)
