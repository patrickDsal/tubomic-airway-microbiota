# =============================================================================
# 05_figure5_clustering.R
# Figure 5: Hierarchical clustering heatmap with dendrogram and annotation bars
#
# The cluster assignments and hclust object are pre-computed in
# 01_prepare_metadata.R. This script handles PLOTTING ONLY.
#
# Requires: 00_setup.R and 01_prepare_metadata.R to have been run first
# Input:    data/metafile_master.csv        ← cluster assignments already here
#           data/hclust_ward2.rds           ← pre-computed hclust object
#           data/otu_table_decontam.csv
#           data/taxonomy.txt
#
# Output:   output/figures/Fig5_Clustering_Heatmap.svg/.pdf
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
library(dendextend)
library(ggdendro)
library(patchwork)
library(RColorBrewer)
library(colorspace)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

source("scripts/00_setup.R")

data_dir   <- "data/"
output_dir <- "output/figures/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)


# -----------------------------------------------------------------------------
# 2. Load master metadata and pre-computed hclust object
# -----------------------------------------------------------------------------

master <- read.csv(file.path(data_dir, "metafile_master.csv"),
                   header = TRUE, row.names = 1,
                   stringsAsFactors = FALSE)

hc_ward2 <- readRDS(file.path(data_dir, "hclust_ward2.rds"))

cat("Loaded hclust object with", length(hc_ward2$labels), "samples\n")


# -----------------------------------------------------------------------------
# 3. Rebuild phyloseq (relative abundance) for OTU heatmap values
# -----------------------------------------------------------------------------

otu_table_raw <- read.table(file.path(data_dir, "otu_table_decontam.csv"),
                            check.names = FALSE, header = TRUE,
                            dec = ".", sep = ",", row.names = 1,
                            comment.char = "")
otu_table_raw[, c("taxonomy", "lineage")] <- NULL
otu_table_raw <- otu_table_raw[, !sapply(otu_table_raw, is.character)]

taxonomy_raw  <- read.table(file.path(data_dir, "taxonomy.txt"),
                            check.names = FALSE, header = TRUE,
                            dec = ".", sep = "\t", row.names = 1,
                            comment.char = "") %>%
  tidy_taxonomy()

# Use only Human samples that appear in the master metadata
sample_meta <- master %>%
  filter(Environment == "Human") %>%
  column_to_rownames("Sample_ID")

dataset <- microtable$new(sample_table = sample_meta,
                          otu_table    = otu_table_raw,
                          tax_table    = taxonomy_raw)
dataset$tidy_dataset()
pseq    <- meco2phyloseq(dataset)
pseq_RA <- transform_sample_counts(pseq, function(x) x / sum(x))

# Apply same unique-genus renaming as in 01_prepare_metadata.R
taxonomy_df       <- as.data.frame(tax_table(pseq_RA))
taxonomy_df$Genus <- ifelse(
  is.na(taxonomy_df$Genus) | taxonomy_df$Genus == "g__",
  "g__Notassigned", taxonomy_df$Genus
)
taxonomy_df$Genus <- ave(
  taxonomy_df$Genus, taxonomy_df$Genus,
  FUN = function(x) paste0(x, seq_along(x))
)
taxonomy_df$Genus <- sub("^g__", "", taxonomy_df$Genus)
tax_table(pseq_RA) <- as.matrix(taxonomy_df)


# -----------------------------------------------------------------------------
# 4. Colour palettes
# -----------------------------------------------------------------------------

base_colors <- c(
  Nei = "#89746A", Sten = "#7A5A3C", Por = "#4B7F4D",  Mix  = "#A44376",
  Kle = "#DBCB9A", Str  = "#D0B84E", Myc = "#4A2C4D",  Sta  = "#C24C5B",
  Bra = "#B26A4E", Ent  = "#5D8D8D", Esch = "#6E8DA0", Bac  = "#8F6B40",
  Hae = "#7B9E5F", Rot  = "#8B5A2B", Pse = "#9E2A2F",  Ser  = "#3E5B5B",
  Pre = "#6E3B58"
)

cluster_colors <- c(
  "Nei1"         = unname(base_colors["Nei"]),
  "Sten1/Sten4"  = darken(unname(base_colors["Sten"]), 0.1),
  "Por2/Str5"    = darken(unname(base_colors["Str"]),  0.1),
  "Mix"          = unname(base_colors["Mix"]),
  "Kle1"         = lighten(unname(base_colors["Kle"]), 0.1),
  "Str1/Por1"    = unname(base_colors["Str"]),
  "Myc1"         = unname(base_colors["Myc"]),
  "Sta1"         = unname(base_colors["Sta"]),
  "Bra1"         = unname(base_colors["Bra"]),
  "Ent1"         = unname(base_colors["Ent"]),
  "Esch1"        = unname(base_colors["Esch"]),
  "Lac1/Esch1"   = lighten(unname(base_colors["Esch"]), 0.3),
  "Bac1/Esch1"   = lighten(unname(base_colors["Esch"]), 0.4),
  "Str2"         = lighten(unname(base_colors["Str"]),  0.2),
  "Hae3/Seg1"    = darken(unname(base_colors["Hae"]),   0.1),
  "Rot4/Myc1"    = unname(base_colors["Rot"]),
  "Cor1/Sta1"    = lighten(unname(base_colors["Sta"]),  0.2),
  "Esch1/Mix"    = lighten(unname(base_colors["Esch"]), 0.2),
  "Str3"         = lighten(unname(base_colors["Str"]),  0.5),
  "Pse1"         = unname(base_colors["Pse"]),
  "Myc4/Kle1"   = darken(unname(base_colors["Myc"]),   0.2),
  "Pre1/Hae2"    = lighten(unname(base_colors["Pre"]),  0.1),
  "Myc2/Urea1"   = darken(unname(base_colors["Myc"]),   0.4),
  "Mix2"         = darken(unname(base_colors["Mix"]),    0.2),
  "Nei2/Esch2"   = lighten(unname(base_colors["Esch"]), 0.2),
  "Kle2/Kle1"    = lighten(unname(base_colors["Kle"]),  0.3),
  "Ser1"         = unname(base_colors["Ser"]),
  "Pre6/Vei1"    = unname(base_colors["Pre"]),
  "Hae1/Hae3"    = darken(unname(base_colors["Hae"]),   0.2),
  "Str4"         = lighten(unname(base_colors["Str"]),   0.7)
)

niche_colors <- c(
  Tube = "#6B8E23", T   = "#F4A300", TS  = "#E56B6F",
  URL  = "#4682B4", LRL = "#5F9EA0", LLL = "#7B68EE"
)

method_colors <- c(Standard = "#f0c571", Venner = "#f0746e")

time_colors <- c(T1 = "#3B6A4D", T2 = "#6A4E23", T3 = "#B19C4E")


# -----------------------------------------------------------------------------
# 5. Patient colours (ordered numerically)
# -----------------------------------------------------------------------------

meta_hm <- master %>%
  filter(Environment == "Human") %>%
  filter(Sample_ID %in% hc_ward2$labels)

# Order by numeric patient ID
patient_order <- unique(meta_hm$Patient_ID[
  order(as.integer(stringr::str_extract(meta_hm$Patient_ID, "\\d+")))
])

patient_colors <- colorRampPalette(brewer.pal(12, "Set3"))(
  length(patient_order)
) %>% setNames(patient_order)


# -----------------------------------------------------------------------------
# 6. Panel A — Dendrogram
# -----------------------------------------------------------------------------

build_dendrogram_panel <- function(hc, hang_height = 0.05,
                                   seg_size = 0.15, height_power = 0.5) {
  transform_y <- function(y) y^height_power

  dend_bl <- hc %>% as.dendrogram() %>%
    hang.dendrogram(hang_height = hang_height) %>% dendro_data()
  dend_gr <- hc %>% as.dendrogram() %>% dendro_data()

  dend_bl$segments$y    <- transform_y(dend_bl$segments$y)
  dend_bl$segments$yend <- transform_y(dend_bl$segments$yend)
  dend_gr$segments$y    <- transform_y(dend_gr$segments$y)
  dend_gr$segments$yend <- transform_y(dend_gr$segments$yend)

  n_samples <- length(hc$order)

  p <- ggplot(segment(dend_gr)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "white", linewidth = seg_size) +
    geom_segment(data = segment(dend_bl),
                 aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = seg_size) +
    scale_x_continuous(expand = rep(1 / n_samples / 2, 2)) +
    scale_y_continuous(expand = c(0, 0.02)) +
    theme_void()

  list(plot = p, hc_order = dend_gr$labels$label)
}

dendro <- build_dendrogram_panel(hc_ward2, seg_size = 0.15, height_power = 0.5)


# -----------------------------------------------------------------------------
# 7. Prepare OTU matrix ordered by dendrogram
# -----------------------------------------------------------------------------

otu_hm      <- as(otu_table(pseq_RA), "matrix")
tax_df      <- as.data.frame(tax_table(pseq_RA))
genus_lookup <- setNames(tax_df$Genus, rownames(tax_df))
rownames(otu_hm) <- genus_lookup[rownames(otu_hm)]

force_include <- c("Rothia4", "Segatella1", "Mycoplasma4",
                   "Mycoplasma2", "Ureaplasma1")
top30_taxa    <- names(sort(rowMeans(otu_hm), decreasing = TRUE))[1:30]
selected_taxa <- unique(c(top30_taxa, force_include))

# Order columns by dendrogram; keep only samples present in both
shared_samples <- intersect(hc_ward2$labels[hc_ward2$order],
                            colnames(otu_hm))
otu_hm_ord     <- otu_hm[selected_taxa, shared_samples]
meta_hm_ord    <- meta_hm[match(shared_samples, meta_hm$Sample_ID), ]
meta_hm_ord$Sample_ID <- factor(meta_hm_ord$Sample_ID,
                                 levels = meta_hm_ord$Sample_ID)


# -----------------------------------------------------------------------------
# 8. Panel G — Heatmap
# -----------------------------------------------------------------------------

hm_data <- otu_hm_ord %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("Taxon") %>%
  pivot_longer(-Taxon, names_to = "Sample_ID", values_to = "RA") %>%
  mutate(
    RA        = if_else(RA == 0, NA_real_, RA),
    Taxon     = fct_inorder(Taxon) %>% fct_rev(),
    Sample_ID = fct_inorder(Sample_ID)
  )

hm <- ggplot(hm_data, aes(x = Sample_ID, y = Taxon, fill = RA)) +
  geom_tile() +
  scale_fill_gradientn(
    name   = "Relative abundance",
    colors = c("grey90", "#f67979", "#c22d2d", "#942222", "#540808"),
    na.value = "white"
  ) +
  theme_grey(base_size = 10) +
  theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.text.y       = element_text(face = "italic"),
    legend.position   = "bottom",
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 6),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width  = unit(0.5, "cm"),
    legend.spacing.x  = unit(0.1, "cm"),
    plot.background   = element_blank()
  ) +
  labs(x = "Samples", y = "Taxa (Genus)")


# -----------------------------------------------------------------------------
# 9. Annotation bar builder
# -----------------------------------------------------------------------------

make_ann_bar <- function(data, y_var, fill_var, color_scale,
                         legend_title = fill_var, show_legend = TRUE) {
  ggplot(data, aes(x = Sample_ID,
                   y = !!sym(y_var),
                   fill = !!sym(fill_var))) +
    theme_void() +
    geom_tile() +
    color_scale +
    theme(
      axis.text.y       = element_text(size = 8),
      legend.position   = if (show_legend) "bottom" else "none",
      legend.title      = element_text(size = 8),
      legend.text       = element_text(size = 6),
      legend.key.height = unit(0.2, "cm"),
      legend.key.width  = unit(0.5, "cm"),
      legend.spacing.x  = unit(0.1, "cm")
    ) +
    labs(fill = legend_title)
}

hm_ann_patient  <- make_ann_bar(meta_hm_ord, "Patient_ID",   "Patient_ID",
                                scale_fill_manual(values = patient_colors),
                                show_legend = FALSE)

hm_ann_time     <- make_ann_bar(meta_hm_ord, "Timepoint",    "Timepoint",
                                scale_fill_manual(values = time_colors))

hm_ann_niche    <- make_ann_bar(meta_hm_ord, "Niche_short",  "Niche_short",
                                scale_fill_manual(values = niche_colors),
                                legend_title = "Niche")

hm_ann_method   <- make_ann_bar(meta_hm_ord, "Method",       "Method",
                                scale_fill_manual(values = method_colors))

hm_ann_cluster  <- make_ann_bar(meta_hm_ord, "Cluster_name", "Cluster_name",
                                scale_fill_manual(values = cluster_colors),
                                legend_title = "Cluster")


# -----------------------------------------------------------------------------
# 10. Assemble Figure 5
# -----------------------------------------------------------------------------

layout <- c(
  area(1,  1, 4,  1),   # dendrogram
  area(5,  1, 5,  1),   # patient
  area(6,  1, 6,  1),   # timepoint
  area(7,  1, 7,  1),   # niche
  area(8,  1, 8,  1),   # method
  area(9,  1, 9,  1),   # cluster
  area(10, 1, 36, 1)    # heatmap
)

fig5 <- dendro$plot +
  hm_ann_patient + hm_ann_time + hm_ann_niche +
  hm_ann_method  + hm_ann_cluster +
  hm +
  plot_layout(design = layout, guides = "collect") &
  theme(legend.position = "bottom")

fig5

ggsave(file.path(output_dir, "Fig5_Clustering_Heatmap.svg"),
       fig5, width = 28, height = 24, units = "cm")
ggsave(file.path(output_dir, "Fig5_Clustering_Heatmap.pdf"),
       fig5, width = 28, height = 24, units = "cm", device = cairo_pdf)

cat("\nFigure 5 saved to", output_dir, "\n")
