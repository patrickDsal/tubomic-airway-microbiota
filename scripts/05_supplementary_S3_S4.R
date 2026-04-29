# =============================================================================
# 05_supplementary_S3_S4.R
# Figure S3: Per-patient cluster transitions T1 → T2 (Standard group)
# Figure S4: Per-patient cluster transitions T1 → T2 (Venner group)
#
# Requires: 04_figure5_clustering.R to have been run first
#           (produces data/metafile_cluster.csv)
# Input:    data/metafile_cluster.csv
#           data/otu_table_decontam.csv
#           data/taxonomy.txt
# Output:   output/figures/FigS3_Cluster_Transitions_Standard.svg/.pdf
#           output/figures/FigS4_Cluster_Transitions_Venner.svg/.pdf
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Libraries
# -----------------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(phyloseq)
library(microeco)
library(file2meco)
library(microbiome)    # meta()
library(ggforce)       # geom_parallel_sets()
library(vegan)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

source("00_setup.R")   # loads theme_pub(), data_dir

data_dir   <- "data/"
output_dir <- "output/figures/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 2. Colour palettes (must match Figure 5)
# -----------------------------------------------------------------------------

niche_colors <- c(
  Tube = "#6B8E23", T   = "#F4A300", TS  = "#E56B6F",
  URL  = "#4682B4", LRL = "#5F9EA0", LLL = "#7B68EE"
)

# Niche order used throughout (Tube excluded — not present at T1/T2)
niche_order <- c("T", "TS", "URL", "LLL", "LRL")


# -----------------------------------------------------------------------------
# 3. Load data and build phyloseq
# -----------------------------------------------------------------------------

sample_info   <- read.csv(file.path(data_dir, "metafile_cluster.csv"),
                          header = TRUE, row.names = 1,
                          stringsAsFactors = FALSE)

otu_table_raw <- read.table(file.path(data_dir, "otu_table_decontam.csv"),
                            check.names = FALSE, header = TRUE,
                            dec = ".", sep = ",", row.names = 1,
                            comment.char = "")
otu_table_raw <- otu_table_raw[, !colnames(otu_table_raw) %in% c("taxonomy", "lineage")]

taxonomy_raw  <- read.table(file.path(data_dir, "taxonomy.txt"),
                            check.names = FALSE, header = TRUE,
                            dec = ".", sep = "\t", row.names = 1,
                            comment.char = "") %>%
  tidy_taxonomy()

dataset <- microtable$new(sample_table = sample_info,
                          otu_table    = otu_table_raw,
                          tax_table    = taxonomy_raw)
dataset$tidy_dataset()

pseq    <- meco2phyloseq(dataset)
pseq_RA <- transform_sample_counts(pseq, function(x) x / sum(x))


# -----------------------------------------------------------------------------
# 4. Build T1 vs T2 cluster comparison table
#    (restricted to patients with complete sample sets)
# -----------------------------------------------------------------------------

meta_df <- meta(pseq_RA) %>%
  filter(Timepoint %in% c("T1", "T2"), FullSet == "Yes")

# One row per patient × niche, with T1 and T2 cluster side by side
comparison_table <- meta_df %>%
  select(Patient_ID, Niche_short, Timepoint, Cluster_name, Method) %>%
  pivot_wider(names_from  = Timepoint,
              values_from = Cluster_name,
              names_prefix = "Cluster_") %>%
  filter(!is.na(Cluster_T1), !is.na(Cluster_T2), !is.na(Niche_short)) %>%
  mutate(Niche_short = factor(Niche_short, levels = niche_order))


# -----------------------------------------------------------------------------
# 5. Helper: build one parallel-sets panel per method
# -----------------------------------------------------------------------------

build_parallel_sets <- function(data, method_label) {

  # Reshape to long format required by geom_parallel_sets
  long_data <- data %>%
    pivot_longer(cols      = c(Cluster_T1, Cluster_T2),
                 names_to  = "Timepoint",
                 values_to = "Cluster") %>%
    mutate(Timepoint = factor(Timepoint,
                              levels = c("Cluster_T1", "Cluster_T2"),
                              labels = c("T1", "T2")))

  ggplot(long_data,
         aes(x = Timepoint, id = Niche_short,
             split = Cluster, value = 1)) +
    geom_parallel_sets(
      aes(fill = Niche_short),
      alpha      = 0.9,
      axis.width = 0.5
    ) +
    geom_parallel_sets_axes(
      axis.width = 0.5,
      fill       = "gray98",
      color      = "gray90"
    ) +
    geom_parallel_sets_labels(
      aes(label = Cluster),
      angle    = 0,
      size     = 4,
      fontface = "bold",
      color    = "black",
      hjust    = 0.5,
      vjust    = 0.5
    ) +
    scale_fill_manual(values = niche_colors, name = "Niche") +
    scale_x_discrete(expand = c(0.01, 0.01)) +
    facet_wrap(~ Patient_ID, shrink = TRUE) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y   = element_blank(),
      axis.title.y  = element_blank(),
      axis.ticks.y  = element_blank(),
      panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.7),
      panel.spacing = unit(0.1, "lines"),
      plot.margin   = margin(2, 2, 2, 2),
      legend.position = "bottom"
    )
}


# -----------------------------------------------------------------------------
# 6. Figure S3 — Standard group
# -----------------------------------------------------------------------------

standard_data <- comparison_table %>% filter(Method == "Standard")

fig_s3 <- build_parallel_sets(standard_data, "Standard")
fig_s3

ggsave(file.path(output_dir, "FigS3_Cluster_Transitions_Standard.svg"),
       fig_s3, width = 26, height = 18, units = "cm")
ggsave(file.path(output_dir, "FigS3_Cluster_Transitions_Standard.pdf"),
       fig_s3, width = 26, height = 18, units = "cm", device = cairo_pdf)


# -----------------------------------------------------------------------------
# 7. Figure S4 — Venner group
# -----------------------------------------------------------------------------

venner_data <- comparison_table %>% filter(Method == "Venner")

fig_s4 <- build_parallel_sets(venner_data, "Venner")
fig_s4

ggsave(file.path(output_dir, "FigS4_Cluster_Transitions_Venner.svg"),
       fig_s4, width = 26, height = 18, units = "cm")
ggsave(file.path(output_dir, "FigS4_Cluster_Transitions_Venner.pdf"),
       fig_s4, width = 26, height = 18, units = "cm", device = cairo_pdf)

cat("\nFigures S3 and S4 saved to", output_dir, "\n")
