# TUBOMIC — Airway Microbiota Analysis

**Study:** TUBOMIC — Influence of an endotracheal tube with subglottic suction and continuous cuff pressure measurement on the pulmonary microbiome

**Registration:** DRKS00029176

**Institution:** Surgical Intensive Care Unit, Heidelberg University Hospital

---

## Overview

This repository contains all R analysis scripts for the TUBOMIC pilot randomized controlled trial, comparing the Venner PneuX® endotracheal tube system (VT) with a standard endotracheal tube (ST) in patients with acute respiratory failure. The primary outcome was the change in airway microbiota beta diversity (Morisita–Horn distances) between tube-associated, upper-airway, and lower-airway communities from intubation (T1) to after four days of ventilation (T2).

---

## Repository Structure

```
tubomic/
├── README.md
├── data/
│   ├── ASV_table.csv               # DADA2 output (ASV × sample counts)
│   ├── taxonomy.txt                # SILVA 138.2 taxonomy assignments
│   ├── sample_info.csv             # Sample metadata
│   └── ps_clean_bio_filtered.rds   # Cleaned phyloseq object (after running 00_setup.R)
├── scripts/
│   ├── 00_setup.R
│   ├── 01_prepare_metadata.R
│   ├── 02_figure2_clinical_scores.R
│   ├── 03_figure3_alpha_diversity.R
│   ├── 04_figure4_beta_network.R
│   ├── 05_figure5_clustering.R
│   ├── 06_supplementary_S1_S2.R
│   ├── 07_supplementary_S5_S6.R
│   ├── 08_supplementary_S3_S7_S8.R
│   └── 09_reviewer_reanalysis_beta.R
└── output/
    └── figures/
```

---

## Run Order

Scripts must be run **in numerical order**. Each script depends on outputs from previous scripts.

```
00_setup.R                    ← Run first. Produces ps_clean_bio_filtered.rds
01_prepare_metadata.R         ← Produces metafile_master.csv, hclust_ward2.rds
02_figure2_clinical_scores.R
03_figure3_alpha_diversity.R
04_figure4_beta_network.R     ← Produces niche_network_full.rds, niche_network_all.rds
05_figure5_clustering.R
06_supplementary_S1_S2.R
07_supplementary_S5_S6.R
08_supplementary_S3_S7_S8.R  ← Requires outputs from 04
09_reviewer_reanalysis_beta.R ← Requires outputs from 04
```

---

## Getting Started

1. Clone this repository
2. Open R and set your working directory to the repository root:
   ```r
   setwd("path/to/tubomic")
   ```
3. Install required packages (see below)
4. Run scripts in order starting from `00_setup.R`

---

## Required R Packages

```r
install.packages(c(
  "tidyverse", "magrittr", "vegan", "lme4", "lmerTest",
  "rcompanion", "effsize", "RColorBrewer", "patchwork",
  "ggraph", "igraph", "ggrepel", "ggdist", "dendextend",
  "ggdendro", "colorspace", "cluster", "fpc", "svglite",
  "flextable", "officer", "conflicted"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("phyloseq", "decontam"))

# GitHub packages
remotes::install_github("ChiLiubio/microeco")
remotes::install_github("jbisanz/file2meco")
```

---

## Figure Index

| Figure | Description | Script |
|--------|-------------|--------|
| Figure 2 | Clinical severity scores (APACHE II, SAPS II, SOFA) | 02 |
| Figure 3 | Alpha diversity: Shannon index and relative dominance | 03 |
| Figure 4 | Beta diversity network graphs | 04 |
| Figure 5 | Hierarchical clustering heatmap (Respirotypes) | 05 |
| Figure S1 | QC: read depth and mock community validation | 06 |
| Figure S2 | Sample availability across patients and niches | 06 |
| Figure S3 | Baseline (T1) niche-pair distance comparison ST vs VT | 08 |
| Figure S4 | Between-group delta distance comparison (reviewer analysis) | 09 |
| Figure S5 | Respirotype cluster transitions — Standard group | 07 |
| Figure S6 | Respirotype cluster transitions — Venner group | 07 |
| Figure S7 | Paired niche distances T1 vs T2 — Standard group | 08 |
| Figure S8 | Paired niche distances T1 vs T2 — Venner group | 08 |

---

## Data Availability

Raw sequencing data are deposited at NCBI SRA under accession number [XXXXX].

---

## Citation

[Citation will be added upon publication]

---

## License

MIT License. See LICENSE for details.
