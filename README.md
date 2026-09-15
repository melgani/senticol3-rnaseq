SentiCOL3 — RNA-seq analysis

This repository contains the R scripts used for the RNA-seq analysis performed as part of the SentiCOL3 study.

The analysis implements a reproducible workflow for transcriptomic data analysis, including differential expression analysis using the limma-voom framework.

Repository contents
senticol3-rnaseq/
├── README.md
├── RNAseq_Senticol_limmavoom_clean.R
├── LICENSE
└── CITATION.cff
Analysis

The main analysis script is:

RNAseq_Senticol_limmavoom_clean.R

The script contains the statistical workflow used for RNA-seq differential expression analysis.

The analysis is based on R and the Bioconductor ecosystem, including the limma package and the voom methodology.

Requirements

The analysis requires:

R
Bioconductor
limma
the additional R/Bioconductor packages required by the analysis script

The versions of R and the packages should be recorded to ensure reproducibility.

Usage

Clone the repository:

git clone https://github.com/melgani/senticol3-rnaseq.git
cd senticol3-rnaseq

Open the main R script:

RNAseq_Senticol_limmavoom_clean.R

Before running the analysis, adapt the input file paths and experimental metadata to the local data structure, if required.

Data availability

The raw and/or processed sequencing data are not included in this repository.

When applicable, sequencing data are available through the corresponding public repository under the accession number reported in the associated publication.

Reproducibility

This repository is intended to support the reproducibility of the analyses reported in the associated research article.

For reproducibility, users should record the versions of:

R
Bioconductor
limma
all other R packages used by the workflow

A versioned release of this repository should preferably be used when citing the analysis.

Citation

If you use this code or reproduce analyses from this repository, please cite the associated research article and this repository.

Contact

For questions regarding the analysis or the repository, please open a GitHub issue or contact the corresponding author of the associated publication.
