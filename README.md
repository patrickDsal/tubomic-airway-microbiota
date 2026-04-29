# Tubomic — Randomised Pilot Study of Airway Microbiota

**Standard vs. Venner PneuX® Endotracheal Tube: Impact on Airway Microbiota**

Reproducible analysis code for the manuscript:

> *Airway microbial compartmentalisation under mechanical ventilation: a randomized pilot study comparing
>a subglottic-suction endotracheal tube (Venner PneuX) with standard intubation*  
> Anna Pilkowski et al., Heidelberg University, Medical Faculty Heidelberg, Department of Anesthesiology

---

## Study overview

Ventilator-associated pneumonia (VAP) is partly driven by microaspiration along 
the endotracheal tube cuff. Because microbial communities differ between airway 
compartments, changes in their similarity over time may serve as an indirect 
readout of microaspiration.

This prospective, randomised, single-centre pilot study (DRKS00029176) compared 
the Venner PneuX® tube (VT; n = 25), combining cuff-pressure monitoring, 
subglottic suction, and a biofilm-resistant coating, against a standard 
endotracheal tube (ST; n = 22) in adults with acute respiratory failure.

Microbial communities from five airway niches (throat, tracheal secretions, 
right upper, right lower, and left lower lung lobes) were sampled by 16S rRNA 
sequencing at intubation (T1), after four days of ventilation (T2), and at the 
tube tip at extubation (T3). Twenty-one patients had complete microbiota datasets 
(9 ST, 12 VT).

**Key finding:** Standard intubation produced progressive microbial convergence 
between airway compartments (e.g. TS–Tube Morisita–Horn distance 0.55 → 0.30, 
p = 0.001); the VT system did not. The findings provide biological plausibility 
for previously reported VAP reductions with the VT system.

---

## Repository structure

```
tubomic/
├── README.md
├── data/
│   ├── raw/
│   │   ├── ASV_table.csv               # DADA2 output (ASV × sample counts)
│   │   ├── taxonomy.txt                # SILVA 138.2 taxonomy assignments
│   │   └── sample_info.csv             # Sample metadata
│   └── processed/
│       └── ps_clean_bio_filtered.rds   # Cleaned phyloseq object (ready to use)
├── scripts/
│   ├── 00_setup.R                      # Data loading, decontam, phyloseq → run first
│   ├── 01_figure2_clinical_scores.R    # Figure 2: APACHE II, SAPS II, SOFA
│   ├── 02_figure3_alpha_diversity.R    # Figure 3: Shannon, Relative Dominance
│   ├── 03_figure4_beta_network.R       # Figure 4: Morisita-Horn network graphs
│   ├── 04_figure5_clustering.R         # Figure 5: Respirotype heatmap
│   ├── 05_supplementary_S1_S2.R        # S1: QC plots; S2: Sample availability
│   ├── 05_supplementary_S3_S4.R        # S3/S4: Per-patient cluster transitions
│   └── 05_supplementary_S5_S6_S7.R    # S5/S6: Beta diversity; S7: ST vs VT at T1
└── output/
    └── figures/                        # All figure outputs written here
```

---

## How to reproduce the figures

### Prerequisites

R ≥ 4.4.2. Install required packages by running:

```r
install.packages(c(
  "tidyverse", "magrittr", "phyloseq", "vegan", "microeco", "file2meco",
  "microbiome", "decontam", "lmerTest", "patchwork", "ggtext", "svglite",
  "ggdist", "ggforce", "ggrepel", "ggraph", "igraph", "ggdendro", "dendextend",
  "RColorBrewer", "colorspace", "scales", "rcompanion", "effsize",
  "rstatix", "cluster", "fpc", "conflicted"
))
```

> **Note:** `microeco` and `file2meco` are available on CRAN.
> `phyloseq` is available via Bioconductor:
> ```r
> BiocManager::install("phyloseq")
> ```

### Option A — Start from raw data (full pipeline)

Run scripts in order from the `scripts/` directory:

```r
source("scripts/00_setup.R")               # ~5 min; produces ps_clean_bio_filtered.rds
source("scripts/01_figure2_clinical_scores.R")
source("scripts/02_figure3_alpha_diversity.R")
source("scripts/03_figure4_beta_network.R")
source("scripts/04_figure5_clustering.R")  # ~10 min; produces metafile_cluster.csv
source("scripts/05_supplementary_S1_S2.R")
source("scripts/05_supplementary_S3_S4.R")
source("scripts/05_supplementary_S5_S6_S7.R")
```

### Option B — Start from the cleaned phyloseq object

If you want to skip preprocessing and jump straight to a specific figure,
load the cleaned phyloseq object directly:

```r
ps <- readRDS("data/processed/ps_clean_bio_filtered.rds")
```

Then run any figure script individually. Each script documents its own inputs
at the top.

---

## Data

| File | Description |
|------|-------------|
| `data/raw/ASV_table.csv` | Raw ASV count matrix (ASVs × samples), DADA2 output |
| `data/raw/taxonomy.txt` | Taxonomic assignments against SILVA 138.2 |
| `data/raw/sample_info.csv` | Sample metadata (patient ID, timepoint, niche, method, DNA concentration) |
| `data/processed/ps_clean_bio_filtered.rds` | Cleaned phyloseq object after decontamination and read-depth filtering (≥ 2,300 reads) |

**Patient privacy:** `sample_info.csv` contains only anonymised patient codes
(e.g., `P01`). No direct identifiers (names, dates of birth, admission dates)
are included.

---

## Bioinformatic processing summary

| Step | Tool / parameter |
|------|-----------------|
| Primer trimming | Cutadapt 4.4, linked adapters, min length 200 bp |
| Quality filtering | DADA2 1.34: truncLen = (220, 200), maxEE = (2, 2) |
| Denoising | DADA2, pool = FALSE (per-sample) |
| Chimera removal | removeBimeraDenovo, method = "consensus" |
| Taxonomy | SILVA 138.2, assignTaxonomy(), minBoot = 50 |
| Decontamination | decontam 1.26, method = "combined", threshold = 0.1 |
| Read-depth filter | ≥ 2,300 reads per sample |
| Diversity metric | Morisita-Horn (vegan::vegdist, method = "horn") |
| Clustering | Ward.D2 hierarchical, k = 30, silhouette-optimised |

---

## Figure index

| Figure | Script | Description |
|--------|--------|-------------|
| Figure 2 | `01_figure2_clinical_scores.R` | APACHE II, SAPS II, SOFA over T1–T3 |
| Figure 3 | `02_figure3_alpha_diversity.R` | Shannon diversity and relative dominance |
| Figure 4 | `03_figure4_beta_network.R` | Morisita-Horn niche similarity networks |
| Figure 5 | `04_figure5_clustering.R` | Respirotype heatmap (k = 30 clusters) |
| Figure S1 | `05_supplementary_S1_S2.R` | Read-depth QC and mock community validation |
| Figure S2 | `05_supplementary_S1_S2.R` | Sample availability across patients and niches |
| Figure S3 | `05_supplementary_S3_S4.R` | Cluster transitions T1 → T2, Standard group |
| Figure S4 | `05_supplementary_S3_S4.R` | Cluster transitions T1 → T2, Venner group |
| Figure S5 | `05_supplementary_S5_S6_S7.R` | Paired niche distances T1 vs T2, Standard |
| Figure S6 | `05_supplementary_S5_S6_S7.R` | Paired niche distances T1 vs T2, Venner |
| Figure S7 | `05_supplementary_S5_S6_S7.R` | Standard vs Venner niche distances at T1 |

---

## Session info

```
R version 4.4.2
Platform: x86_64-pc-linux-gnu
```

Full session info can be reproduced by running `sessionInfo()` after sourcing
any figure script.

---

## Contact

Patrick Schaal  
Institute of Medical Microbiology, University of Luebeck, Germany

---

## License

Code: MIT License  
Data: available for research use; please cite the manuscript if used.
