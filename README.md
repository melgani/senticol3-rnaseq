# SENTICOL III — RNA-seq analysis

This repository contains the R scripts used for the RNA-seq analysis performed as part of the **SENTICOL III GINECO trial**.

The analysis is associated with the study:

> **Molecular alterations and relapse risk in early-staged cervical cancer from the SENTICOL III GINECO trial**

The purpose of this repository is to provide the computational code used for the transcriptomic analyses and to facilitate transparency and reproducibility of the research.

## Study

SENTICOL III is a multicenter GINECO clinical trial investigating molecular and clinical factors associated with relapse risk and outcomes in patients with early-stage cervical cancer.

This repository specifically contains the computational workflow used for the RNA-seq analysis.

## RNA-seq analysis

RNA-seq data were analyzed using R and the Bioconductor ecosystem.

Differential gene expression analysis was performed using the **limma-voom** framework.

The main analysis script is:

```text
RNAseq_Senticol_limmavoom_clean.R
```

## Repository contents

```text
senticol3-rnaseq/
│
├── CITATION.CFF
├── README.md
└── RNAseq_Senticol_limmavoom_clean.R
```

Additional files and documentation may be added as the repository is updated.

## Requirements

The analysis requires:

* R
* Bioconductor
* limma
* voom
* Other R/Bioconductor packages required by the analysis script

The analysis should preferably be performed using the same or compatible versions of R and the required packages to ensure reproducibility.

## Usage

Clone the repository:

```bash
git clone https://github.com/melgani/senticol3-rnaseq.git
```

Then move to the repository directory:

```bash
cd senticol3-rnaseq
```

The main analysis script can be opened and executed in R or RStudio:

```text
RNAseq_Senticol_limmavoom_clean.R
```

Before running the script, the input data paths and other parameters may need to be adapted to the local environment.

## Data availability

Patient-level clinical data and raw sequencing data are **not included in this repository**.

The data are subject to the applicable ethical, legal, and data-protection requirements governing the SENTICOL III GINECO trial.

Where applicable, information regarding access to the underlying molecular or clinical datasets will be provided in the associated publication.

## Reproducibility

This repository provides the R code used for the RNA-seq analyses associated with the study.

For reproducibility, users should record the versions of:

* R
* Bioconductor
* limma
* voom
* Other R/Bioconductor packages used in the analysis

The repository may be updated over time. When available, users are encouraged to use the specific version or release of the repository corresponding to the published analysis.

## Citation

If you use this code or reproduce analyses from this repository, please cite the associated research article:

> **El Gani M, Ibadioune S, Hamza A, et al. Molecular alterations and relapse risk in early-staged cervical cancer from the SENTICOL III GINECO trial.**

Please also cite the specific version of this repository used for your analysis.

A `CITATION.cff` file will be provided to facilitate citation of the repository.

## Authors

**Maryame El Gani**, **Sabrina Ibadioune**, **Abderaouf Hamza**, **Zakhia El Beaino**, **Sophie Vacher**, **Anne Schnitzler**, **Emmanuelle Jeannot**, **Julien Masliah-Planchon**, **Vincent Cockenpot**, **Alexandre Degnieau**, **Gwenaël Ferron**, **François Golfier**, **Eric Lambaudie**, **Fabrice Narducci**, **Cécile Loaec**, **Jennifer Uzan**, **Frederic Marchal**, **Anne-Sophie Bats**, **Martin Koskas**, **Virginie Fourchotte**, **Estelle Wafo**, **Nicolas Bourdel**, **Raffaèle Fauvet**, **Marie Plante**, **Patrice Mathevet**, **Maud Kamal**, **Fabrice Lecuru**, **Ivan Bieche**.

## License

The source code is made available under the license specified in the `LICENSE` file.

## Contact

For questions regarding the computational analysis or this repository, please open an issue on GitHub or contact the authors of the associated publication.

## Repository

https://github.com/melgani/senticol3-rnaseq
