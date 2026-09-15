SENTICOL III — RNA-seq analysis

This repository contains the R scripts used for the RNA-seq analysis performed as part of the SENTICOL III GINECO trial.

The analyses contribute to the study:

Molecular alterations and relapse risk in early-staged cervical cancer from the SENTICOL III GINECO trial

The repository is intended to facilitate transparency and reproducibility of the computational analyses reported in the associated research work.

Study

SENTICOL III is a multicenter clinical trial investigating molecular and clinical factors associated with outcomes in patients with early-stage cervical cancer.

The present analysis focuses on transcriptomic data and differential gene expression analysis.

Analysis workflow

The main analysis was performed in R using the limma-voom framework.

The workflow includes the processing and statistical analysis of RNA-seq data and the identification of differentially expressed genes between the relevant experimental groups.

Repository structure
senticol3-rnaseq/
│
├── README.md
├── CITATION.cff
├── LICENSE
│
├── scripts/
│   └── RNAseq_Sententicol_limmavoom_clean.R
│
├── data/
│   └── README.md
│
└── results/
    └── README.md

Note: The repository does not contain identifiable patient-level data or raw sequencing files.

Main script

The main analysis script is:

scripts/RNAseq_Senticol_limmavoom_clean.R

It implements the RNA-seq differential expression analysis using the limma and voom methodology.

Requirements

The analysis requires:

R
Bioconductor
limma
edgeR
and the additional R packages required by the analysis script

The exact package versions should be recorded to facilitate reproducibility.

Data availability

Raw sequencing data and/or patient-level clinical data are not distributed in this repository.

Where applicable, the sequencing data are deposited in an appropriate public repository under the accession number reported in the associated publication.

Because the study involves human samples, access to clinical and molecular data may be subject to ethical, legal, and data-protection restrictions.

Reproducibility

This repository provides the computational code used for the RNA-seq analyses associated with the study.

For reproducibility, users should use the versioned release of this repository corresponding to the publication.

The software environment should include the R and Bioconductor versions and package versions used to perform the analysis.

Citation

If you use the code or reproduce analyses from this repository, please cite the associated research article.

Please also cite the specific version of this repository used for the analysis.

Citation information is provided in CITATION.cff.

Authors

Maryame El Gani, Sabrina Ibadioune, Abderaouf Hamza, Zakhia El Beaino, Sophie Vacher, Anne Schnitzler, Emmanuelle Jeannot, Julien Masliah-Planchon, Vincent Cockenpot, Alexandre Degnieau, Gwenaël Ferron, François Golfier, Eric Lambaudie, Fabrice Narducci, Cécile Loaec, Jennifer Uzan, Frederic Marchal, Anne-Sophie Bats, Martin Koskas, Virginie Fourchotte, Estelle Wafo, Nicolas Bourdel, Raffaèle Fauvet, Marie Plante, Patrice Mathevet, Maud Kamal, Fabrice Lecuru, Ivan Bieche.

License

The source code is distributed under the license specified in LICENSE.

Contact

For questions regarding the computational analysis, please open an issue in this repository or contact the corresponding author of the associated publication.
