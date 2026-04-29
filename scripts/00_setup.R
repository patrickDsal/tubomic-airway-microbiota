# =============================================================================
# 00_setup.R
# Data loading, quality control, decontamination, and phyloseq construction
#
# Study: ST vs VT intubation — airway microbiota pilot RCT
# Output: ps_clean_bio_filtered.rds  (clean phyloseq object used by all
#         downstream figure scripts)
#
# Run this script FIRST before any figure script.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Libraries
# -----------------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(phyloseq)
library(microeco)
library(file2meco)
library(decontam)
library(ggplot2)
library(ggdist)
library(patchwork)
library(scales)

# Resolve common namespace conflicts
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("lag",    "dplyr")


# -----------------------------------------------------------------------------
# 2. Shared ggplot theme (sourced by all figure scripts)
# -----------------------------------------------------------------------------

theme_pub <- function(base_size = 18) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.text.x      = element_text(size = rel(0.9)),
      axis.text.y      = element_text(size = rel(0.9)),
      axis.title.x     = element_text(size = rel(1.0), face = "bold"),
      axis.title.y     = element_text(size = rel(1.0), face = "bold"),
      legend.text      = element_text(size = rel(0.9)),
      legend.title     = element_text(size = rel(1.0)),
      plot.title       = element_text(hjust = 0.5, size = rel(1.3), face = "bold"),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.25)
    )
}

# =============================================================================
# Resolve namespace conflicts — declared once here, applies to all scripts
# =============================================================================
library(conflicted)

conflicts_prefer(
  # dplyr wins for data manipulation
  dplyr::filter,
  dplyr::select,
  dplyr::lag,
  dplyr::rename,
  
  # purrr wins for functional programming
  purrr::set_names,
  purrr::map,
  
  # lmerTest::lmer gives p-values (wraps lme4), always prefer it
  lmerTest::lmer,
  
  # rstatix for stats helpers
  rstatix::filter
)

# -----------------------------------------------------------------------------
# 3. Load raw input files
# -----------------------------------------------------------------------------
setwd("")
# Adjust this path to wherever your input files live
data_dir <- "data/"

sample_info    <- read.csv(file.path(data_dir, "sample_info.csv"),
                           row.names = NULL, header = TRUE,
                           stringsAsFactors = FALSE)

asv_table      <- read.csv(file.path(data_dir, "ASV_table.csv"),
                           row.names = 1, header = TRUE,
                           stringsAsFactors = FALSE)

taxonomy_table <- read.table(file.path(data_dir, "taxonomy.txt"),
                             check.names = FALSE, header = TRUE,
                             dec = ".", sep = "\t",
                             row.names = 1, comment.char = "") %>%
  tidy_taxonomy()


# -----------------------------------------------------------------------------
# 4. Sanity check: match sample IDs between sample_info and ASV table
# -----------------------------------------------------------------------------

sample_info_ids  <- sample_info$Sample_ID
asv_colnames     <- colnames(asv_table)

cat("--- Sample ID matching ---\n")
cat("Matching IDs:        ", length(intersect(sample_info_ids, asv_colnames)), "\n")
cat("In sample_info only: ", length(setdiff(sample_info_ids, asv_colnames)), "\n")
cat("In ASV table only:   ", length(setdiff(asv_colnames, sample_info_ids)),  "\n")


# -----------------------------------------------------------------------------
# 5. Build microeco object and initial filtering
# -----------------------------------------------------------------------------

rownames(sample_info) <- sample_info[[1]]
sample_info           <- sample_info[, -1]

dataset <- microtable$new(
  sample_table = sample_info,
  otu_table    = asv_table,
  tax_table    = taxonomy_table
)
dataset$tidy_dataset()

# Keep only Bacteria; remove plastid and environmental contaminants
dataset$tax_table %<>% base::subset(Kingdom == "k__Bacteria")
dataset$filter_pollution(
  taxa = c("mitochondria", "Chloroplast", "Cyanobacteria",
           "Enhydrobacter", "Cutibacterium")
)
dataset$tidy_dataset()

cat("\n--- After initial filtering ---\n")
cat("Samples:", ncol(dataset$otu_table), "\n")
cat("ASVs:   ", nrow(dataset$otu_table), "\n")
cat("Reads:  ", sum(as.matrix(dataset$otu_table)), "\n")


# -----------------------------------------------------------------------------
# 6. Decontamination with decontam (combined frequency + prevalence method)
# -----------------------------------------------------------------------------

# Subset to human samples + controls (exclude mock community positive controls
# from decontam; they would confound the frequency model)
dataset_for_decontam <- clone(dataset)
dataset_for_decontam$sample_table <- subset(
  dataset_for_decontam$sample_table,
  Environment %in% c("Human", "negative_control", "positive_control")
)
dataset_for_decontam$tidy_dataset()

# Convert to phyloseq for decontam
ps_raw <- meco2phyloseq(dataset_for_decontam)

# Only samples with a measured DNA concentration can be used in the combined model
ps_with_conc <- prune_samples(!is.na(sample_data(ps_raw)$concentration), ps_raw)
is_neg        <- sample_data(ps_with_conc)$Environment == "negative_control"

contam_result <- isContaminant(
  ps_with_conc,
  method = "combined",
  neg    = is_neg,
  conc   = sample_data(ps_with_conc)$concentration,
  threshold = 0.1
)

contaminant_asvs <- taxa_names(ps_raw)[contam_result$contaminant]
cat("\n--- Decontamination ---\n")
cat("Contaminant ASVs identified:", length(contaminant_asvs), "\n")

# Save contaminant list for transparency
contam_tax <- as.data.frame(tax_table(ps_raw)[contaminant_asvs, ])
contam_tax$ASV <- contaminant_asvs
write.csv(contam_tax, file.path(data_dir, "Contaminant_ASVs.csv"), row.names = FALSE)

# Remove contaminants and controls from the phyloseq object
ps_clean <- prune_taxa(!taxa_names(ps_raw) %in% contaminant_asvs, ps_raw)
ps_clean <- prune_samples(sample_data(ps_clean)$Environment == "Human", ps_clean)

cat("Samples after removing controls:", nsamples(ps_clean), "\n")
cat("ASVs after decontam:            ", ntaxa(ps_clean),    "\n")


# -----------------------------------------------------------------------------
# 7. Read-depth QC and low-depth sample removal (threshold: 2 300 reads)
# -----------------------------------------------------------------------------

READ_THRESHOLD <- 2300

reads_per_sample <- data.frame(
  Sample   = names(sample_sums(ps_clean)),
  Reads    = sample_sums(ps_clean)
)
reads_per_sample$log_reads <- log10(reads_per_sample$Reads + 1)

ps_clean_filtered <- prune_samples(sample_sums(ps_clean) >= READ_THRESHOLD, ps_clean)

cat("\n--- After read-depth filtering (threshold:", READ_THRESHOLD, "reads) ---\n")
cat("Samples:", nsamples(ps_clean_filtered), "\n")
cat("ASVs:   ", ntaxa(ps_clean_filtered),    "\n")
cat("Reads:  ", sum(otu_table(ps_clean_filtered)), "\n")


# -----------------------------------------------------------------------------
# 8. Save cleaned phyloseq object and supporting tables
# -----------------------------------------------------------------------------

saveRDS(ps_clean_filtered, file.path(data_dir, "ps_clean_bio_filtered.rds"))

write.csv(as.data.frame(otu_table(ps_clean_filtered)),
          file.path(data_dir, "otu_table_decontam.csv"))

write.csv(as.data.frame(tax_table(ps_clean_filtered)),
          file.path(data_dir, "tax_table_decontam.csv"))

sample_data_df <- as(sample_data(ps_clean_filtered), "data.frame")
sample_data_df <- merge(sample_data_df,
                        reads_per_sample[, c("Sample", "Reads")],
                        by.x = "row.names", by.y = "Sample",
                        all.x = TRUE)
write.csv(sample_data_df,
          file.path(data_dir, "sample_data_decontam.csv"),
          row.names = FALSE)

cat("\nSetup complete. Cleaned phyloseq saved to data/ps_clean_bio_filtered.rds\n")
