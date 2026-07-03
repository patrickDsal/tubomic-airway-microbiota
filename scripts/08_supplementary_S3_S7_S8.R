# =============================================================================
# 08_supplementary_S5_S6_S7.R
# Figure S5: Paired niche-niche distances T1 → T2, Standard group
# Figure S6: Paired niche-niche distances T1 → T2, Venner group
# Figure S7: ST vs VT baseline (T1) comparison of niche-niche distances
#
# Requires: 04_figure4_beta_network.R to have been run first
#           (produces niche_network_full.rds and niche_network_all.rds)
#
# S5/S6: all available samples, one observation per patient per niche pair
#         per timepoint (deduplicated)
# S7:    all available samples, baseline T1 comparison ST vs VT
#
# Input:    data/niche_network_full.rds
#           data/niche_network_all.rds
# Output:   output/figures/FigS7_Beta_Standard_Paired.svg/.pdf
#           output/figures/FigS8_Beta_Venner_Paired.svg/.pdf
#           output/figures/FigS3_Beta_Baseline_Comparison.svg/.pdf
#           data/Standard_Beta_Stats_Wilcox_ALL.csv
#           data/Venner_Beta_Stats_Wilcox_ALL.csv
#           data/Baseline_Beta_Stats_ST_vs_VT.csv
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
library(rcompanion)
library(effsize)
library(RColorBrewer)
library(patchwork)
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
# 3. Load pre-computed distance data from script 04
# -----------------------------------------------------------------------------

niche_network_full <- readRDS(file.path(data_dir, "niche_network_full.rds"))
niche_network_all  <- readRDS(file.path(data_dir, "niche_network_all.rds"))


# -----------------------------------------------------------------------------
# 4. Core analysis function
#    Uses all available samples, deduplicated to one observation per
#    patient × niche pair × timepoint — statistically correct pairing
# -----------------------------------------------------------------------------

run_beta_analysis <- function(niche_net, method_sel) {

  niche_network_method <- niche_net %>%
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

  # Wilcoxon + effect sizes — identical to script 04
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

  list(
    filtered_data = filtered_data,
    stats         = summary_stats %>%
      left_join(stats_full, by = "Niche_Pair") %>%
      mutate(
        p_adj = p.adjust(p_value, method = "BH"),
        sig   = case_when(
          p_adj < 0.001 ~ "***",
          p_adj < 0.01  ~ "**",
          p_adj < 0.05  ~ "*",
          TRUE          ~ "ns"
        )
      )
  )
}


# -----------------------------------------------------------------------------
# 5. Run analysis
#    S5/S6 plots use full set for visualisation
#    Stats exported use all available samples
# -----------------------------------------------------------------------------

res_ST_full <- run_beta_analysis(niche_network_full, "Standard")
res_VT_full <- run_beta_analysis(niche_network_full, "Venner")
res_ST_all  <- run_beta_analysis(niche_network_all,  "Standard")
res_VT_all  <- run_beta_analysis(niche_network_all,  "Venner")

write.csv(res_ST_all$stats,
          file.path(data_dir, "Standard_Beta_Stats_Wilcox_ALL.csv"),
          row.names = FALSE)
write.csv(res_VT_all$stats,
          file.path(data_dir, "Venner_Beta_Stats_Wilcox_ALL.csv"),
          row.names = FALSE)

cat("\n=== Standard — all available ===\n"); print(res_ST_all$stats)
cat("\n=== Venner — all available ===\n");   print(res_VT_all$stats)


# -----------------------------------------------------------------------------
# 6. Paired boxplot helper
# -----------------------------------------------------------------------------

build_paired_plot <- function(filtered_data, stats_df, plot_title) {

  pt_colors <- colorRampPalette(brewer.pal(12, "Set3"))(
    length(unique(filtered_data$Patient1))
  ) %>% setNames(sort(unique(filtered_data$Patient1)))

  p_values <- filtered_data %>%
    group_by(Niche_Pair) %>%
    summarise(max_distance = max(Distance, na.rm = TRUE), .groups = "drop") %>%
    left_join(stats_df %>% select(Niche_Pair, p_value),
              by = "Niche_Pair") %>%
    mutate(
      y_position = max_distance + 0.2,
      label      = sprintf("p = %.3f", p_value)
    )

  ggplot(filtered_data, aes(x = Time_Comparison, y = Distance)) +
    geom_boxplot(alpha = 0.9, outlier.shape = NA,
                 position = position_dodge(width = 0.3)) +
    geom_point(aes(group = Patient1, color = Patient1),
               position = position_dodge(width = 0.3),
               size = 2.5, alpha = 1) +
    geom_line(aes(group = Patient1, color = Patient1),
              position = position_dodge(width = 0.3),
              alpha = 0.5, linewidth = 1, linetype = "dashed") +
    geom_text(data = p_values,
              aes(x = 1.5, y = y_position, label = label),
              inherit.aes = FALSE, size = 4) +
    facet_wrap(~ Niche_Pair, scales = "free_y") +
    scale_color_manual(values = pt_colors) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                       limits = c(0, 1.3)) +
    labs(title = plot_title,
         x     = "Timepoint",
         y     = "Distance (Morisita-Horn)") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position  = "none",
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text       = element_text(face = "bold"),
      panel.border     = element_rect(color = "black", fill = NA,
                                      linewidth = 1.2),
      panel.grid.major = element_line(color = "gray80", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 16, hjust = 0.5)
    )
}


# -----------------------------------------------------------------------------
# 7. Figure S5 — Standard (full set)
# -----------------------------------------------------------------------------

figS5 <- build_paired_plot(res_ST_full$filtered_data,
                           res_ST_full$stats,
                           "Standard — paired niche distances T1 vs T2")
figS5

ggsave(file.path(output_dir, "FigS7_Beta_Standard_Paired.svg"),
       figS5, width = 22, height = 18, units = "cm")
ggsave(file.path(output_dir, "FigS7_Beta_Standard_Paired.pdf"),
       figS5, width = 22, height = 18, units = "cm", device = cairo_pdf)


# -----------------------------------------------------------------------------
# 8. Figure S6 — Venner (full set)
# -----------------------------------------------------------------------------

figS6 <- build_paired_plot(res_VT_full$filtered_data,
                           res_VT_full$stats,
                           "Venner — paired niche distances T1 vs T2")
figS6

ggsave(file.path(output_dir, "FigS8_Beta_Venner_Paired.svg"),
       figS6, width = 22, height = 18, units = "cm")
ggsave(file.path(output_dir, "FigS8_Beta_Venner_Paired.pdf"),
       figS6, width = 22, height = 18, units = "cm", device = cairo_pdf)


# -----------------------------------------------------------------------------
# 9. Figure S7 — ST vs VT at baseline T1 (all available samples)
#    Exactly as in original script
# -----------------------------------------------------------------------------

niche_network_filtered <- niche_network_all %>%
  filter(Timepoint1 == "T1", Timepoint2 == "T1",
         Niche1 != Niche2)

niche_network_long <- niche_network_filtered %>%
  pivot_longer(cols      = c(Sample1, Sample2),
               names_to  = "Sample_Type",
               values_to = "Sample_ID") %>%
  mutate(
    Timepoint = "T1",
    Method    = ifelse(Sample_Type == "Sample1", Method1, Method2),
    Niche     = ifelse(Sample_Type == "Sample1", Niche1,  Niche2),
    Patient   = ifelse(Sample_Type == "Sample1", Patient1, Patient2)
  ) %>%
  mutate(Niche_Pair = ifelse(Niche1 < Niche2,
                             paste(Niche1, Niche2, sep = "_"),
                             paste(Niche2, Niche1, sep = "_"))) %>%
  distinct(Patient, Niche_Pair, Timepoint, Method, .keep_all = TRUE)

p_values_s7 <- niche_network_long %>%
  group_by(Niche_Pair) %>%
  filter(all(c("Standard", "Venner") %in% Method)) %>%
  summarise(
    p_value = wilcox.test(
      Distance[Method == "Standard"],
      Distance[Method == "Venner"],
      paired = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj   = p.adjust(p_value, method = "BH"),
    p_label = formatC(p_value, format = "f", digits = 3)
  )

write.csv(p_values_s7,
          file.path(data_dir, "Baseline_Beta_Stats_ST_vs_VT.csv"),
          row.names = FALSE)

max_distance_s7   <- max(niche_network_long$Distance, na.rm = TRUE)
y_position_s7     <- max_distance_s7 + 0.1
n_patients_s7     <- length(unique(niche_network_long$Patient))
patient_colors_s7 <- colorRampPalette(brewer.pal(12, "Set3"))(n_patients_s7)

figS7 <- ggplot(niche_network_long,
                aes(x = Method, y = Distance, fill = Method)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA,
               position = position_dodge(width = 0.3)) +
  geom_point(aes(group = Patient, color = Patient),
             position = position_dodge(width = 0.3), size = 2) +
  geom_line(aes(group = Patient, color = Patient),
            position = position_dodge(width = 0.3),
            alpha = 0.5, linewidth = 1) +
  geom_text(data = p_values_s7,
            aes(x = 1.5, y = y_position_s7,
                label = paste0("p = ", p_label)),
            inherit.aes = FALSE, size = 4) +
  facet_wrap(~ Niche_Pair) +
  scale_fill_manual(values  = method_colors) +
  scale_color_manual(values = patient_colors_s7) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                     limits = c(0, 1.15)) +
  labs(title = "Baseline (T1) — Standard vs. Venner",
       x     = NULL,
       y     = "Distance (Morisita-Horn)") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text       = element_text(face = "bold"),
    panel.border     = element_rect(color = "black", fill = NA,
                                    linewidth = 1.2),
    panel.grid.major = element_line(color = "gray80", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 16, hjust = 0.5)
  )

figS7

ggsave(file.path(output_dir, "FigS3_Beta_Baseline_Comparison.svg"),
       figS7, width = 22, height = 18, units = "cm")
ggsave(file.path(output_dir, "FigS3_Beta_Baseline_Comparison.pdf"),
       figS7, width = 22, height = 18, units = "cm", device = cairo_pdf)

cat("\nFigures S5, S6, S7 saved to", output_dir, "\n")
