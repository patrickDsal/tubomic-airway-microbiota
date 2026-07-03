# =============================================================================
# 04_figure4_beta_network.R
# Figure 4: Beta Diversity Similarity Graph
#
# Two datasets:
#   Full set (metafile_master_FullSet.csv)  → network visualisation only
#   All samples (metafile_master.csv)       → statistics and asterisks
#
# Requires: 00_setup.R and 01_prepare_metadata.R to have been run first
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
library(igraph)
library(ggraph)
library(ggrepel)
library(patchwork)
library(rcompanion)
library(effsize)
library(RColorBrewer)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

source("scripts/00_setup.R")

data_dir   <- "data/"
output_dir <- "output/figures/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 2. Colour palettes
# -----------------------------------------------------------------------------

niche_colors <- c(
  Tube = "#6B8E23", T   = "#F4A300", TS  = "#E56B6F",
  URL  = "#4682B4", LRL = "#5F9EA0", LLL = "#7B68EE"
)

method_colors <- c(Standard = "#f0c571", Venner = "#f0746e")


# -----------------------------------------------------------------------------
# 3. Shared OTU and taxonomy tables
# -----------------------------------------------------------------------------

otu_table_1 <- read.table(file.path(data_dir, "otu_table_decontam.csv"),
                          check.names = FALSE, header = TRUE,
                          dec = ".", sep = ",", row.names = 1,
                          comment.char = "")
otu_table_1[, "taxonomy"] <- NULL
otu_table_1[, "lineage"]  <- NULL

taxonomy_table_1 <- read.table(file.path(data_dir, "taxonomy.txt"),
                               check.names = FALSE, header = TRUE,
                               dec = ".", sep = "\t", row.names = 1,
                               comment.char = "") %>%
  tidy_taxonomy()


# -----------------------------------------------------------------------------
# 4. Helper: build niche_network_df from a metadata CSV
#    Faithful to original script — no Niche1 != Niche2 filter,
#    no upper-triangle deduplication
# -----------------------------------------------------------------------------

build_niche_network <- function(meta_csv) {

  sample_info <- read.csv(file.path(data_dir, meta_csv),
                          header = TRUE, row.names = 1,
                          stringsAsFactors = FALSE)
  rownames(sample_info) <- sample_info$ID

  dataset <- microtable$new(sample_table = sample_info,
                            otu_table    = otu_table_1,
                            tax_table    = taxonomy_table_1)
  dataset$tidy_dataset()

  pseq    <- meco2phyloseq(dataset)
  pseq_RA <- transform_sample_counts(pseq, function(x) x / sum(x))

  meta_np <- meta(pseq_RA)
  meta_np <- rownames_to_column(meta_np, var = "Sample_ID")

  otu_np <- as(otu_table(pseq_RA), "matrix")
  if (taxa_are_rows(pseq_RA)) otu_np <- t(otu_np)
  otu_np <- as.data.frame(otu_np)

  ord_np <- vegdist(otu_np, method = "horn")

  distance_df <- as.data.frame(as.table(as.matrix(ord_np))) %>%
    rename(Sample1 = Var1, Sample2 = Var2, Distance = Freq) %>%
    filter(Sample1 != Sample2)

  distance_df <- distance_df %>%
    left_join(meta_np %>% select(Sample_ID, Patient_ID, Niche_short,
                                 Timepoint, Method),
              by = c("Sample1" = "Sample_ID")) %>%
    rename(Patient1 = Patient_ID, Niche1 = Niche_short,
           Timepoint1 = Timepoint, Method1 = Method) %>%
    left_join(meta_np %>% select(Sample_ID, Patient_ID, Niche_short,
                                 Timepoint, Method),
              by = c("Sample2" = "Sample_ID")) %>%
    rename(Patient2 = Patient_ID, Niche2 = Niche_short,
           Timepoint2 = Timepoint, Method2 = Method)

  # Exactly as in original script
  distance_df %>%
    filter(Patient1 == Patient2,
           Niche1   != "control",
           Niche2   != "control")
}


# -----------------------------------------------------------------------------
# 5. Build the two distance matrices
# -----------------------------------------------------------------------------

cat("Building distance matrix — full set (network plots)...\n")
niche_network_full <- build_niche_network("metafile_master_FullSet.csv")
saveRDS(niche_network_full, file.path(data_dir, "niche_network_full.rds"))

cat("Building distance matrix — all available samples (statistics)...\n")
niche_network_all <- build_niche_network("metafile_master.csv")
saveRDS(niche_network_all, file.path(data_dir, "niche_network_all.rds"))


# -----------------------------------------------------------------------------
# 6. Network visualisation (full set)
#    Mean distances per method × timepoint, deduplicated for plotting only
# -----------------------------------------------------------------------------

get_mean_distances <- function(df, method_sel, tp_sel) {
  df %>%
    filter(Method1 == method_sel, Method2 == method_sel) %>%
    filter(
      (Timepoint1 == tp_sel & Timepoint2 == tp_sel) |
      (Timepoint1 == "T3"   & Timepoint2 == tp_sel) |
      (Timepoint1 == tp_sel & Timepoint2 == "T3")
    ) %>%
    filter(Niche1 != Niche2) %>%
    mutate(sorted_niches = pmap_chr(list(Niche1, Niche2),
                                    ~ paste(sort(c(..1, ..2)), collapse = "-"))) %>%
    group_by(sorted_niches) %>%
    summarise(mean_dist = mean(Distance, na.rm = TRUE),
              sd_dist   = sd(Distance,   na.rm = TRUE),
              .groups   = "drop") %>%
    separate(sorted_niches, into = c("NicheA", "NicheB"), sep = "-")
}

dist_ST_T1 <- get_mean_distances(niche_network_full, "Standard", "T1")
dist_ST_T2 <- get_mean_distances(niche_network_full, "Standard", "T2")
dist_VT_T1 <- get_mean_distances(niche_network_full, "Venner",   "T1")
dist_VT_T2 <- get_mean_distances(niche_network_full, "Venner",   "T2")


# -----------------------------------------------------------------------------
# 7. Network plot helper
# -----------------------------------------------------------------------------

scale_distance_inverse <- function(dist, new_min = 0.2, new_max = 1.5) {
  ((1 - dist) / (1 - 0)) * (new_max - new_min) + new_min
}

make_network_plot <- function(dist_summary, title_text, edge_color = "#888888") {

  g <- graph_from_data_frame(
    dist_summary %>% rename(from = NicheA, to = NicheB, weight = mean_dist),
    directed = FALSE
  )

  set.seed(123)
  ggraph(g, layout = "graphopt") +
    geom_edge_link(
      aes(edge_width = scale_distance_inverse(weight),
          edge_alpha = scale_distance_inverse(weight)),
      color       = edge_color,
      show.legend = FALSE
    ) +
    geom_node_point(aes(color = name), size = 8, stroke = 1.5) +
    geom_node_text(aes(label = name), repel = TRUE,
                   size = 5, fontface = "bold") +
    scale_color_manual(values = niche_colors) +
    theme_void() +
    labs(title = title_text) +
    theme(
      plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "none",
      plot.margin     = margin(10, 10, 10, 10)
    )
}

pA <- make_network_plot(dist_ST_T1, "A  Standard — T1",
                        edge_color = method_colors["Standard"])
pB <- make_network_plot(dist_ST_T2, "B  Standard — T2",
                        edge_color = method_colors["Standard"])
pC <- make_network_plot(dist_VT_T1, "C  Venner — T1",
                        edge_color = method_colors["Venner"])
pD <- make_network_plot(dist_VT_T2, "D  Venner — T2",
                        edge_color = method_colors["Venner"])


# -----------------------------------------------------------------------------
# 8. Statistics (all available samples)
#    Faithful to original script — no distinct(), no arrange()
# -----------------------------------------------------------------------------

run_beta_stats <- function(niche_network_df, method_sel) {

  niche_network_method <- niche_network_df %>%
    filter(Method1 == method_sel, Method2 == method_sel) %>%
    mutate(
      Niche_Pair = ifelse(Niche1 < Niche2,
                          paste(Niche1, Niche2, sep = "_"),
                          paste(Niche2, Niche1, sep = "_")),
      Time_Comparison = case_when(
        (Timepoint1 == "T1" & Timepoint2 == "T1") ~ "T1",
        (Timepoint1 == "T2" & Timepoint2 == "T2") ~ "T2",
        (Timepoint1 == "T1" & Timepoint2 == "T2") |
        (Timepoint1 == "T2" & Timepoint2 == "T1") ~ "T1_T2",
        (Timepoint1 == "T3" & Timepoint2 == "T1") |
        (Timepoint1 == "T1" & Timepoint2 == "T3") ~ "T1",
        (Timepoint1 == "T3" & Timepoint2 == "T2") |
        (Timepoint1 == "T2" & Timepoint2 == "T3") ~ "T2",
        TRUE ~ NA_character_
      )
    )

  filtered_data <- niche_network_method %>%
    filter(Time_Comparison %in% c("T1", "T2"))

  filtered_data <- filtered_data %>%
    group_by(Patient1, Niche_Pair) %>%
    filter(all(c("T1", "T2") %in% Time_Comparison)) %>%
    ungroup()

  # Summary statistics
  summary_stats <- filtered_data %>%
    group_by(Niche_Pair) %>%
    filter(all(c("T1", "T2") %in% Time_Comparison)) %>%
    summarise(
      n_pairs = n_distinct(Patient1),
      mean_T1 = mean(Distance[Time_Comparison == "T1"], na.rm = TRUE),
      sd_T1   = sd(Distance[Time_Comparison == "T1"],   na.rm = TRUE),
      mean_T2 = mean(Distance[Time_Comparison == "T2"], na.rm = TRUE),
      sd_T2   = sd(Distance[Time_Comparison == "T2"],   na.rm = TRUE),
      .groups = "drop"
    )

  # Wilcoxon + effect sizes — exactly as in original script
  stats_full <- filtered_data %>%
    group_by(Niche_Pair) %>%
    filter(all(c("T1", "T2") %in% Time_Comparison)) %>%
    summarise(
      p_value = wilcox.test(
        Distance[Time_Comparison == "T1"],
        Distance[Time_Comparison == "T2"],
        paired = TRUE
      )$p.value,
      cohen_d = cohen.d(
        Distance[Time_Comparison == "T1"],
        Distance[Time_Comparison == "T2"],
        paired = TRUE
      )$estimate,
      rank_biserial = wilcoxonR(
        Distance[Time_Comparison == "T1"],
        Distance[Time_Comparison == "T2"],
        paired = TRUE
      ),
      .groups = "drop"
    )

  summary_stats %>%
    left_join(stats_full, by = "Niche_Pair") %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           sig   = case_when(
             p_adj < 0.001 ~ "***",
             p_adj < 0.01  ~ "**",
             p_adj < 0.05  ~ "*",
             TRUE          ~ "ns"
           ))
}

stats_ST <- run_beta_stats(niche_network_all, "Standard")
stats_VT <- run_beta_stats(niche_network_all, "Venner")

write.csv(stats_ST,
          file.path(data_dir, "Standard_Beta_Stats_Wilcox.csv"),
          row.names = FALSE)
write.csv(stats_VT,
          file.path(data_dir, "Venner_Beta_Stats_Wilcox.csv"),
          row.names = FALSE)

cat("\n=== Standard group ===\n"); print(stats_ST)
cat("\n=== Venner group ===\n");   print(stats_VT)


# -----------------------------------------------------------------------------
# 9. Assemble and save Figure 4
# -----------------------------------------------------------------------------

fig4 <- (pA | pB) / (pC | pD)

ggsave(file.path(output_dir, "Fig4_Beta_Network.svg"),
       fig4, width = 28, height = 24, units = "cm")
ggsave(file.path(output_dir, "Fig4_Beta_Network.pdf"),
       fig4, width = 28, height = 24, units = "cm", device = cairo_pdf)

cat("\nFigure 4 saved to", output_dir, "\n")
