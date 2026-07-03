# =============================================================================
# 01_prepare_metadata.R
# Master metadata preparation — run SECOND after 00_setup.R
#
# This script consolidates all metadata enrichment steps into one place so
# that every downstream figure script can simply load data/metafile_master.csv
# without depending on another figure script having run first.
#
# Steps performed:
#   1. Load cleaned phyloseq object (from 00_setup.R)
#   2. Identify patients with complete sample sets (FullSet flag)
#   3. Compute alpha diversity metrics (Shannon, Relative Dominance, etc.)
#   4. Compute Morisita-Horn distances and hierarchical clustering (k = 30)
#   5. Assign cluster names
#   6. Save master metadata file
#
# Input:  data/ps_clean_bio_filtered.rds   (from 00_setup.R)
#         data/otu_table_decontam.csv
#         data/taxonomy.txt
#         data/sample_data_decontam.csv
#
# Output: data/metafile_master.csv         ← single source of truth for all
#                                             downstream scripts
#         data/hclust_ward2.rds            ← saved hclust object (reused by
#                                             05_figure5_clustering.R)
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Libraries
# -----------------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(phyloseq)
library(microeco)
library(file2meco)
library(microbiome)
library(vegan)
library(cluster)
library(fpc)
library(RColorBrewer)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

source("scripts/00_setup.R")   # loads theme_pub(), data_dir

data_dir <- "data/"


# =============================================================================
# STEP 1 — Load cleaned phyloseq object
# =============================================================================

cat("\n--- Loading cleaned phyloseq object ---\n")

ps <- readRDS(file.path(data_dir, "ps_clean_bio_filtered.rds"))

cat("Samples:", nsamples(ps), "\n")
cat("ASVs:   ", ntaxa(ps),    "\n")


# =============================================================================
# STEP 2 — Identify patients with complete sample sets
#
# A complete set = T1 and T2 samples for all 5 niches (T, TS, URL, LRL, LLL)
#                 PLUS a T3 Tube sample
# =============================================================================

cat("\n--- Computing FullSet flag ---\n")

REQUIRED_NICHE_TIMES <- c(
  "T1_LLL", "T2_LLL",
  "T1_LRL", "T2_LRL",
  "T1_URL", "T2_URL",
  "T1_TS",  "T2_TS",
  "T1_T",   "T2_T",
  "T3_Tube"
)

sample_df <- as(sample_data(ps), "data.frame") %>%
  rownames_to_column("Sample_ID")

# Build Niche_Time label if not already present
if (!"Niche_Time" %in% colnames(sample_df)) {
  sample_df <- sample_df %>%
    mutate(Niche_Time = paste(Timepoint, Niche_short, sep = "_"))
}

fullset_flag <- sample_df %>%
  filter(Environment == "Human") %>%
  group_by(Patient_ID) %>%
  summarise(
    FullSet = ifelse(
      all(REQUIRED_NICHE_TIMES %in% Niche_Time), "Yes", "No"
    ),
    .groups = "drop"
  )

sample_df <- sample_df %>%
  left_join(fullset_flag, by = "Patient_ID")

n_full <- sum(fullset_flag$FullSet == "Yes", na.rm = TRUE)
cat("Patients with complete datasets:", n_full, "\n")
cat("  Standard:", sum(
  sample_df$FullSet == "Yes" & sample_df$Method == "Standard",
  na.rm = TRUE
) / 11, "\n")   # 11 samples per patient
cat("  Venner:  ", sum(
  sample_df$FullSet == "Yes" & sample_df$Method == "Venner",
  na.rm = TRUE
) / 11, "\n")


# =============================================================================
# STEP 3 — Alpha diversity metrics
# =============================================================================

cat("\n--- Computing alpha diversity ---\n")

# Rebuild taxonomy for unique genus labels (needed for clustering step too)
taxonomy_raw <- read.table(
  file.path(data_dir, "taxonomy.txt"),
  check.names = FALSE, header = TRUE,
  dec = ".", sep = "\t", row.names = 1,
  comment.char = ""
) %>% tidy_taxonomy()

otu_table_raw <- read.table(
  file.path(data_dir, "otu_table_decontam.csv"),
  check.names = FALSE, header = TRUE,
  dec = ".", sep = ",", row.names = 1,
  comment.char = ""
)
otu_table_raw[, c("taxonomy", "lineage")] <-
  lapply(c("taxonomy", "lineage"),
         function(col) if (col %in% colnames(otu_table_raw)) NULL)
otu_table_raw <- otu_table_raw[, !sapply(otu_table_raw, is.character)]

# Use the sample_df (with FullSet flag) as the sample table
rownames(sample_df) <- sample_df$Sample_ID

dataset <- microtable$new(
  sample_table = sample_df,
  otu_table    = otu_table_raw,
  tax_table    = taxonomy_raw
)
dataset$tidy_dataset()
pseq <- meco2phyloseq(dataset)

# Standard phyloseq alpha diversity
alpha_div <- estimate_richness(
  pseq,
  measures = c("Shannon", "InvSimpson", "Chao1", "Observed")
)

# Relative dominance: proportion of most abundant taxon
otu_mat <- as.matrix(otu_table(pseq))
alpha_div$Relative_Dominance <-
  apply(otu_mat, 2, max) / colSums(otu_mat)

# Evenness
alpha_div$Evenness <- alpha_div$Shannon / log(alpha_div$Observed)

cat("Alpha diversity computed for", nrow(alpha_div), "samples\n")


# =============================================================================
# STEP 4 — Hierarchical clustering (Morisita-Horn, ward.D2, k = 30)
# =============================================================================

cat("\n--- Running hierarchical clustering ---\n")

# Rename ASVs by unique genus label (for interpretable cluster names)
taxonomy_df       <- as.data.frame(tax_table(pseq))
taxonomy_df$Genus <- ifelse(
  is.na(taxonomy_df$Genus) | taxonomy_df$Genus == "g__",
  "g__Notassigned",
  taxonomy_df$Genus
)
taxonomy_df$Genus <- ave(
  taxonomy_df$Genus, taxonomy_df$Genus,
  FUN = function(x) paste0(x, seq_along(x))
)
# Strip g__ prefix for readability
taxonomy_df$Genus <- sub("^g__", "", taxonomy_df$Genus)
tax_table(pseq) <- as.matrix(taxonomy_df)

pseq_RA <- transform_sample_counts(pseq, function(x) x / sum(x))

otu_np <- as(otu_table(pseq_RA), "matrix")
if (taxa_are_rows(pseq_RA)) otu_np <- t(otu_np)
otu_np <- as.data.frame(otu_np)

ord_np <- vegdist(otu_np, method = "horn")

# Ward.D2 selected based on cophenetic correlation (see 05_figure5_clustering.R)
set.seed(123)
hc_ward2 <- hclust(ord_np, method = "ward.D2")

# Save for reuse in 05_figure5_clustering.R
saveRDS(hc_ward2, file.path(data_dir, "hclust_ward2.rds"))
cat("Hierarchical clustering complete — hclust object saved\n")

# Cut at k = 30
K_FINAL  <- 30
clusters <- cutree(hc_ward2, k = K_FINAL)

# Map cluster numbers to names (curated from top genera per cluster)
cluster_map <- c(
  "1"  = "Nei1",        "2"  = "Sten1/Sten4", "3"  = "Por2/Str5",
  "4"  = "Mix",         "5"  = "Kle1",         "6"  = "Str1/Por1",
  "7"  = "Myc1",        "8"  = "Sta1",          "9"  = "Bra1",
  "10" = "Ent1",        "11" = "Esch1",         "12" = "Lac1/Esch1",
  "13" = "Bac1/Esch1",  "14" = "Str2",          "15" = "Hae3/Seg1",
  "16" = "Rot4/Myc1",   "17" = "Cor1/Sta1",     "18" = "Esch1/Mix",
  "19" = "Str3",        "20" = "Pse1",           "21" = "Myc4/Kle1",
  "22" = "Pre1/Hae2",   "23" = "Myc2/Urea1",    "24" = "Mix2",
  "25" = "Nei2/Esch2",  "26" = "Kle2/Kle1",     "27" = "Ser1",
  "28" = "Pre6/Vei1",   "29" = "Hae1/Hae3",     "30" = "Str4"
)

# Assign clusters to samples in the same order as the distance matrix
sample_cluster <- data.frame(
  Sample_ID    = hc_ward2$labels,
  Cluster      = clusters[match(hc_ward2$labels,
                                names(clusters))],
  stringsAsFactors = FALSE
) %>%
  mutate(Cluster_name = cluster_map[as.character(Cluster)])

cat("Cluster sizes:\n")
print(table(sample_cluster$Cluster_name))


# =============================================================================
# STEP 5 — Assemble master metadata file
# =============================================================================

cat("\n--- Assembling master metadata file ---\n")

# Base: sample_df already has FullSet flag
master <- sample_df %>%
  # Join alpha diversity
  left_join(
    alpha_div %>% rownames_to_column("Sample_ID"),
    by = "Sample_ID"
  ) %>%
  # Join cluster assignments
  left_join(
    sample_cluster,
    by = "Sample_ID"
  )

# Convenience columns used by multiple downstream scripts
master <- master %>%
  mutate(
    # Short identifier used in some plots
    Identifier = Sample_ID,
    # Combined niche × timepoint label
    Niche_Time = if_else(
      is.na(Niche_Time),
      paste(Timepoint, Niche_short, sep = "_"),
      Niche_Time
    )
  )

master <- master %>% column_to_rownames("Sample_ID")
write.csv(master, file.path(data_dir, "metafile_master.csv"), row.names = TRUE)

cat("Master metadata saved: data/metafile_master.csv\n")
cat("Rows:", nrow(master), "\n")
cat("Columns:", ncol(master), "\n")

# Also save FullSet-only subset for scripts that need it

master_fullset <- master %>% filter(FullSet == "Yes", Environment == "Human")
write.csv(master_fullset, file.path(data_dir, "metafile_master_FullSet.csv"), row.names = TRUE)

cat("Full-set subset saved: data/metafile_master_FullSet.csv\n")
cat("Rows:", nrow(master_fullset), "\n")

# Print quick summary
cat("\n--- Quick summary ---\n")
master %>%
  filter(Environment == "Human") %>%
  group_by(Method, FullSet) %>%
  summarise(n_samples  = n(),
            n_patients = n_distinct(Patient_ID),
            .groups    = "drop") %>%
  print()

cat("\n01_prepare_metadata.R complete.\n")
cat("Suggested run order:\n")
cat("  00_setup.R\n")
cat("  01_prepare_metadata.R   ← this script\n")
cat("  02_figure2_clinical_scores.R\n")
cat("  03_figure3_alpha_diversity.R\n")
cat("  04_figure4_beta_network.R\n")
cat("  05_figure5_clustering.R\n")
cat("  06_supplementary_S1_S2.R\n")
cat("  07_supplementary_S3_S4.R\n")
cat("  08_supplementary_S5_S6_S7.R\n")
cat("  09_reviewer_reanalysis_beta.R\n")
