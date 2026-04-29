# =============================================================================
# 04_figure5_clustering.R
# Figure 5: Hierarchical clustering heatmap with dendrogram and annotation bars
#           (Respirotype analysis)
#
# Requires: 00_setup.R to have been run first
# Input:    data/ps_clean_bio_filtered.rds
# Output:   output/figures/Fig5_Clustering_Heatmap.svg/.pdf
#           data/metafile_cluster.csv       (sample-level cluster assignments)
#           data/top5_genera_summary.csv    (cluster characterisation table)
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Libraries
# -----------------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(phyloseq)
library(microeco)
library(file2meco)
library(vegan)
library(dendextend)
library(ggdendro)
library(patchwork)
library(RColorBrewer)
library(colorspace)   # lighten() / darken()
library(microbiome)   # meta()
library(conflicted)

conflict_prefer("filter",  "dplyr")
conflict_prefer("select",  "dplyr")

source("00_setup.R")   # loads theme_pub(), data_dir

data_dir   <- "data/"
output_dir <- "output/figures/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)


# -----------------------------------------------------------------------------
# 2. Load data and rebuild microeco / phyloseq objects
# -----------------------------------------------------------------------------

sample_info    <- read.csv(file.path(data_dir, "sample_data_decontam.csv"),
                           header = TRUE, row.names = 1,
                           stringsAsFactors = FALSE)

otu_table_raw  <- read.table(file.path(data_dir, "otu_table_decontam.csv"),
                             check.names = FALSE, header = TRUE,
                             dec = ".", sep = ",", row.names = 1,
                             comment.char = "")
# Drop non-count columns if present
otu_table_raw[, c("taxonomy", "lineage")] <-
  lapply(c("taxonomy", "lineage"),
         function(col) if (col %in% colnames(otu_table_raw)) NULL)
otu_table_raw <- otu_table_raw[, !sapply(otu_table_raw, is.character)]

taxonomy_raw   <- read.table(file.path(data_dir, "taxonomy.txt"),
                             check.names = FALSE, header = TRUE,
                             dec = ".", sep = "\t", row.names = 1,
                             comment.char = "") %>%
  tidy_taxonomy()

dataset <- microtable$new(sample_table = sample_info,
                          otu_table    = otu_table_raw,
                          tax_table    = taxonomy_raw)
dataset$tidy_dataset()
pseq <- meco2phyloseq(dataset)


# -----------------------------------------------------------------------------
# 3. Rename ASVs by unique genus labels
#    (ensures each ASV row has a readable, non-duplicated genus name)
# -----------------------------------------------------------------------------

taxonomy_df       <- as.data.frame(tax_table(pseq))

# Replace missing/unassigned genus with a placeholder
taxonomy_df$Genus <- ifelse(
  is.na(taxonomy_df$Genus) | taxonomy_df$Genus == "g__",
  "g__Notassigned",
  taxonomy_df$Genus
)

# Make genus names unique by appending a sequence number within each genus
taxonomy_df$Genus <- ave(
  taxonomy_df$Genus, taxonomy_df$Genus,
  FUN = function(x) paste0(x, seq_along(x))
)
tax_table(pseq) <- as.matrix(taxonomy_df)

# Strip the "g__" prefix for display
taxonomy_df$Genus <- sub("^g__", "", taxonomy_df$Genus)
tax_table(pseq)   <- as.matrix(taxonomy_df)


# -----------------------------------------------------------------------------
# 4. Compute Morisita-Horn distance matrix and evaluate linkage methods
# -----------------------------------------------------------------------------

pseq_RA  <- transform_sample_counts(pseq, function(x) x / sum(x))
otu_np   <- as(otu_table(pseq_RA), "matrix")
if (taxa_are_rows(pseq_RA)) otu_np <- t(otu_np)
otu_np   <- as.data.frame(otu_np)
ord_np   <- vegdist(otu_np, method = "horn")

# Compare cophenetic correlations across linkage methods
linkage_methods <- c("ward.D", "ward.D2", "single", "complete",
                     "average", "mcquitty", "median", "centroid")

coph_cors <- sapply(linkage_methods, function(m) {
  hc  <- hclust(ord_np, method = m)
  cor(ord_np, cophenetic(hc), method = "spearman")
})
cat("\n--- Cophenetic correlations ---\n")
print(round(coph_cors, 3))
# ward.D2 selected as final linkage method

hc.ward2 <- hclust(ord_np, method = "ward.D2")


# -----------------------------------------------------------------------------
# 5. Select optimal k with Calinski-Harabasz and silhouette indices
# -----------------------------------------------------------------------------

compute_cluster_indices <- function(dist, k_range, method = "ward.D2") {
  results <- lapply(k_range, function(k) {
    labels <- cutree(hclust(dist, method = method), k)
    stats  <- fpc::cluster.stats(dist, labels)
    data.frame(k              = k,
               `Calinski-Harabasz` = stats$ch,
               Silhouette     = stats$avg.silwidth)
  })
  bind_rows(results)
}

k_range      <- 2:50
clust_indices <- compute_cluster_indices(ord_np, k_range)

# Plot indices (diagnostic — not saved to manuscript figures)
clust_diag <- clust_indices %>%
  pivot_longer(-k, names_to = "Index", values_to = "Value") %>%
  ggplot(aes(x = k, y = Value)) +
  geom_point() + geom_line() +
  facet_wrap(~ Index, scales = "free_y") +
  scale_x_continuous(breaks = seq(2, 50, 2)) +
  labs(title = "Cluster index diagnostics",
       x = "Number of clusters (k)", y = "Metric value") +
  theme_pub()

ggsave(file.path(output_dir, "Fig5_cluster_index_diagnostics.svg"),
       clust_diag, width = 18, height = 10, units = "cm")


# -----------------------------------------------------------------------------
# 6. Assign k = 30 clusters and name them by dominant genus
# -----------------------------------------------------------------------------

K_FINAL <- 30

clusters <- cutree(hc.ward2, k = K_FINAL)
sample_data(pseq_RA)$Cluster <- as.factor(
  clusters[match(hc.ward2$labels, sample_names(pseq_RA))]
)

# Manually curated cluster names (dominant genus / top-2 genera per cluster)
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

sample_data_df <- meta(pseq_RA)
sample_data_df$Cluster_name <- cluster_map[as.character(sample_data_df$Cluster)]
sample_data(pseq_RA)        <- sample_data(sample_data_df)

write.csv(sample_data_df,
          file.path(data_dir, "metafile_cluster.csv"),
          row.names = TRUE)


# -----------------------------------------------------------------------------
# 7. Cluster summary table (top 5 genera per cluster)
# -----------------------------------------------------------------------------

otu_long <- otu_np %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "Abundance")

# Top 5 genera per cluster
top5_genera <- otu_long %>%
  left_join(sample_data_df %>% rownames_to_column("sample_id") %>%
              select(sample_id, Cluster),
            by = "sample_id") %>%
  group_by(Cluster, ASV) %>%
  summarise(mean_abundance = mean(Abundance),
            sd_abundance   = sd(Abundance), .groups = "drop") %>%
  group_by(Cluster) %>%
  slice_max(order_by = mean_abundance, n = 5) %>%
  left_join(as.data.frame(tax_table(pseq_RA)) %>%
              rownames_to_column("ASV") %>% select(ASV, Genus),
            by = "ASV") %>%
  mutate(mean_abundance_percent = mean_abundance * 100)

# Niche and timepoint breakdowns
niche_counts <- sample_data_df %>%
  group_by(Cluster, Niche_short) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Niche_short, values_from = n, values_fill = 0)

timepoint_counts <- sample_data_df %>%
  group_by(Cluster, Timepoint) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Timepoint, values_from = n, values_fill = 0)

cluster_metadata <- sample_data_df %>%
  group_by(Cluster) %>%
  summarise(n_samples   = n(),
            n_patients  = n_distinct(Patient_ID),
            n_venner    = sum(Method == "Venner",   na.rm = TRUE),
            n_standard  = sum(Method == "Standard", na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Cluster_name = cluster_map[as.character(Cluster)])

top5_summary <- top5_genera %>%
  left_join(cluster_metadata,   by = "Cluster") %>%
  left_join(niche_counts,       by = "Cluster") %>%
  left_join(timepoint_counts,   by = "Cluster")

write.csv(top5_summary,
          file.path(data_dir, "top5_genera_summary.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 8. Colour palettes
# -----------------------------------------------------------------------------

# Base colours keyed by genus prefix
base_colors <- c(
  Nei  = "#89746A", Sten = "#7A5A3C", Por  = "#4B7F4D", Mix  = "#A44376",
  Kle  = "#DBCB9A", Str  = "#D0B84E", Myc  = "#4A2C4D", Sta  = "#C24C5B",
  Bra  = "#B26A4E", Ent  = "#5D8D8D", Esch = "#6E8DA0", Bac  = "#8F6B40",
  Hae  = "#7B9E5F", Rot  = "#8B5A2B", Pse  = "#9E2A2F", Ser  = "#3E5B5B",
  Pre  = "#6E3B58"
)

cluster_colors <- c(
  "Nei1"         = base_colors["Nei"],
  "Sten1/Sten4"  = darken(base_colors["Sten"], 0.1),
  "Por2/Str5"    = darken(base_colors["Str"],  0.1),
  "Mix"          = base_colors["Mix"],
  "Kle1"         = lighten(base_colors["Kle"], 0.1),
  "Str1/Por1"    = base_colors["Str"],
  "Myc1"         = base_colors["Myc"],
  "Sta1"         = base_colors["Sta"],
  "Bra1"         = base_colors["Bra"],
  "Ent1"         = base_colors["Ent"],
  "Esch1"        = base_colors["Esch"],
  "Lac1/Esch1"   = lighten(base_colors["Esch"], 0.3),
  "Bac1/Esch1"   = lighten(base_colors["Esch"], 0.4),
  "Str2"         = lighten(base_colors["Str"],  0.2),
  "Hae3/Seg1"    = darken(base_colors["Hae"],   0.1),
  "Rot4/Myc1"    = base_colors["Rot"],
  "Cor1/Sta1"    = lighten(base_colors["Sta"],  0.2),
  "Esch1/Mix"    = lighten(base_colors["Esch"], 0.2),
  "Str3"         = lighten(base_colors["Str"],  0.5),
  "Pse1"         = base_colors["Pse"],
  "Myc4/Kle1"    = darken(base_colors["Myc"],   0.2),
  "Pre1/Hae2"    = lighten(base_colors["Pre"],  0.1),
  "Myc2/Urea1"   = darken(base_colors["Myc"],   0.4),
  "Mix2"         = darken(base_colors["Mix"],    0.2),
  "Nei2/Esch2"   = lighten(base_colors["Esch"], 0.2),
  "Kle2/Kle1"    = lighten(base_colors["Kle"],  0.3),
  "Ser1"         = base_colors["Ser"],
  "Pre6/Vei1"    = base_colors["Pre"],
  "Hae1/Hae3"    = darken(base_colors["Hae"],   0.2),
  "Str4"         = lighten(base_colors["Str"],   0.7)
)
cluster_colors <- unname(cluster_colors) |>
  setNames(names(cluster_colors))   # strip named-vector nesting from base_colors

niche_colors <- c(
  Tube = "#6B8E23", T = "#F4A300", TS  = "#E56B6F",
  URL  = "#4682B4", LRL = "#5F9EA0", LLL = "#7B68EE"
)

method_colors <- c(Standard = "#f0c571", Venner = "#f0746e")

time_colors <- c(T1 = "#3B6A4D", T2 = "#6A4E23", T3 = "#B19C4E")

patient_colors <- colorRampPalette(brewer.pal(12, "Set3"))(
  length(unique(sample_data_df$Patient_ID))
) |> setNames(
  unique(sample_data_df$Patient_ID[
    order(as.integer(str_extract(sample_data_df$Patient_ID, "\\d+")))
  ])
)


# -----------------------------------------------------------------------------
# 9. Panel A — Dendrogram
# -----------------------------------------------------------------------------

build_dendrogram_panel <- function(hc, hang_height = 0.05,
                                   seg_size = 0.15, height_power = 0.5) {
  dend_bl <- hc %>% as.dendrogram() %>%
    hang.dendrogram(hang_height = hang_height) %>%
    dendro_data()

  dend_gr <- hc %>% as.dendrogram() %>% dendro_data()

  # Power-transform heights to compress top and expand bottom
  transform_y <- function(y) y^height_power
  for (d in list(dend_bl, dend_gr)) {
    d$segments$y    <- transform_y(d$segments$y)
    d$segments$yend <- transform_y(d$segments$yend)
  }
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

dendro <- build_dendrogram_panel(hc.ward2, seg_size = 0.15, height_power = 0.5)


# -----------------------------------------------------------------------------
# 10. Prepare OTU matrix ordered by dendrogram
# -----------------------------------------------------------------------------

otu_hm      <- as(otu_table(pseq_RA), "matrix")
tax_df      <- as.data.frame(tax_table(pseq_RA))
genus_lookup <- setNames(tax_df$Genus, rownames(tax_df))
rownames(otu_hm) <- genus_lookup[rownames(otu_hm)]

# Top 30 taxa by mean abundance, plus a few manually forced inclusions
force_include <- c("Rothia4", "Segatella1", "Mycoplasma4",
                   "Mycoplasma2", "Ureaplasma1")
top30_taxa    <- names(sort(rowMeans(otu_hm), decreasing = TRUE))[1:30]
selected_taxa <- unique(c(top30_taxa, force_include))

# Order samples by dendrogram
otu_hm_ord  <- otu_hm[selected_taxa, hc.ward2$order]
meta_hm_ord <- meta(pseq_RA)[hc.ward2$order, ]
meta_hm_ord$sample_id <- factor(rownames(meta_hm_ord),
                                 levels = rownames(meta_hm_ord))

stopifnot(all(colnames(otu_hm_ord) == hc.ward2$labels[hc.ward2$order]))


# -----------------------------------------------------------------------------
# 11. Panel G — Heatmap
# -----------------------------------------------------------------------------

hm_data <- otu_hm_ord %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("Taxon") %>%
  pivot_longer(-Taxon, names_to = "sample_id", values_to = "RA") %>%
  mutate(
    RA        = if_else(RA == 0, NA_real_, RA),
    Taxon     = fct_inorder(Taxon) %>% fct_rev(),
    sample_id = fct_inorder(sample_id)
  )

hm <- ggplot(hm_data, aes(x = sample_id, y = Taxon, fill = RA)) +
  geom_tile() +
  scale_fill_gradientn(
    name   = "Relative abundance",
    colors = c("grey90", "#f67979", "#c22d2d", "#942222", "#540808"),
    na.value = "white"
  ) +
  theme_grey(base_size = 10) +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(face = "italic"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 6),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width  = unit(0.5, "cm"),
    legend.spacing.x  = unit(0.1, "cm"),
    plot.background   = element_blank()
  ) +
  labs(x = "Samples", y = "Taxa (Genus)")


# -----------------------------------------------------------------------------
# 12. Annotation bar builder (shared helper)
# -----------------------------------------------------------------------------

make_ann_bar <- function(data, y_var, fill_var, color_scale,
                         legend_title = fill_var, show_legend = TRUE) {
  ggplot(data, aes(x = sample_id, y = !!sym(y_var), fill = !!sym(fill_var))) +
    theme_void() +
    geom_tile() +
    color_scale +
    theme(
      axis.text.y      = element_text(size = 8),
      legend.position  = if (show_legend) "bottom" else "none",
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 6),
      legend.key.height = unit(0.2, "cm"),
      legend.key.width  = unit(0.5, "cm"),
      legend.spacing.x  = unit(0.1, "cm")
    ) +
    labs(fill = legend_title)
}

hm_ann1 <- make_ann_bar(meta_hm_ord, "Patient_ID",   "Patient_ID",
                         scale_fill_manual(values = patient_colors),
                         show_legend = FALSE)

hm_ann2 <- make_ann_bar(meta_hm_ord, "Timepoint",    "Timepoint",
                         scale_fill_manual(values = time_colors))

hm_ann3 <- make_ann_bar(meta_hm_ord, "Niche_short",  "Niche_short",
                         scale_fill_manual(values = niche_colors),
                         legend_title = "Niche")

hm_ann4 <- make_ann_bar(meta_hm_ord, "Method",       "Method",
                         scale_fill_manual(values = method_colors))

hm_ann5 <- make_ann_bar(meta_hm_ord, "Cluster_name", "Cluster_name",
                         scale_fill_manual(values = cluster_colors),
                         legend_title = "Cluster")


# -----------------------------------------------------------------------------
# 13. Assemble Figure 5
# -----------------------------------------------------------------------------

# Layout: dendrogram (A) on top, annotation rows (B-F), heatmap (G) at bottom
# Relative heights match the original figure proportions
layout <- c(
  area(1,  1, 4,  1),   # A  dendrogram
  area(5,  1, 5,  1),   # B  patient
  area(6,  1, 6,  1),   # C  timepoint
  area(7,  1, 7,  1),   # D  niche
  area(8,  1, 8,  1),   # E  method
  area(9,  1, 9,  1),   # F  cluster
  area(10, 1, 36, 1)    # G  heatmap
)

fig5 <- dendro$plot +
  hm_ann1 + hm_ann2 + hm_ann3 + hm_ann4 + hm_ann5 +
  hm +
  plot_layout(design = layout, guides = "collect") &
  theme(legend.position = "bottom")

fig5

ggsave(file.path(output_dir, "Fig5_Clustering_Heatmap.svg"),
       fig5, width = 28, height = 24, units = "cm")
ggsave(file.path(output_dir, "Fig5_Clustering_Heatmap.pdf"),
       fig5, width = 28, height = 24, units = "cm", device = cairo_pdf)

cat("\nFigure 5 saved to", output_dir, "\n")
