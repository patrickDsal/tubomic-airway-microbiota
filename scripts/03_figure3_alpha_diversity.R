# =============================================================================
# 03_figure3_alpha_diversity.R
# Figure 3: Alpha diversity (Shannon + Relative Dominance) across niches
#           and timepoints, with paired Wilcoxon tests
#
# Requires: 00_setup.R and 01_prepare_metadata.R to have been run first
# Input:    data/metafile_master.csv        ← single source of truth
#           data/otu_table_decontam.csv
#           data/taxonomy.txt
#
# Output:   output/figures/Fig3A_Shannon.svg/.pdf
#           output/figures/Fig3B_RelativeDominance.svg/.pdf
#           data/stats_Shannon_paired_wilcoxon.csv
#           data/stats_RelativeDominance_paired_wilcoxon.csv
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
library(rcompanion)
library(effsize)
library(RColorBrewer)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

source("scripts/00_setup.R")   # loads theme_pub(), data_dir

data_dir   <- "data/"
output_dir <- "output/figures/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 2. Load master metadata
#    Alpha diversity metrics are already computed in 01_prepare_metadata.R
# -----------------------------------------------------------------------------

master <- read.csv(file.path(data_dir, "metafile_master.csv"),
                   header = TRUE, row.names = 1,
                   stringsAsFactors = FALSE)

# Figure 3 uses only patients with complete sample sets, T1 and T2 only,
# with matched T1/T2 pairs per niche
analysis_df <- master %>%
  filter(
    Environment  == "Human",
    FullSet      == "Yes",
    Timepoint    %in% c("T1", "T2"),
    Niche_short  != "Tube"
  ) %>%
  group_by(Patient_ID, Niche_short) %>%
  filter(all(c("T1", "T2") %in% Timepoint)) %>%
  ungroup() %>%
  mutate(
    Niche_short = factor(Niche_short,
                         levels = c("T", "TS", "URL", "LRL", "LLL")),
    Method      = factor(Method, levels = c("Standard", "Venner")),
    Patient_ID  = factor(Patient_ID)
  )

cat("Patients in analysis:", n_distinct(analysis_df$Patient_ID), "\n")
cat("  Standard:", n_distinct(analysis_df$Patient_ID[analysis_df$Method == "Standard"]), "\n")
cat("  Venner:  ", n_distinct(analysis_df$Patient_ID[analysis_df$Method == "Venner"]),   "\n")

# Patient colour palette
patient_colors <- colorRampPalette(brewer.pal(12, "Set3"))(
  nlevels(analysis_df$Patient_ID)
) %>% setNames(levels(analysis_df$Patient_ID))


# -----------------------------------------------------------------------------
# 3. Paired Wilcoxon test helper
#    Runs per Niche × Method, returns p-value, effect size, and n_paired
# -----------------------------------------------------------------------------

run_paired_wilcoxon <- function(df, metric) {

  df %>%
    group_by(Niche_short, Method) %>%
    group_modify(~ {
      x <- .x[[metric]][.x$Timepoint == "T1"]
      y <- .x[[metric]][.x$Timepoint == "T2"]

      ok <- complete.cases(x, y)
      x  <- x[ok]; y <- y[ok]
      n  <- length(x)

      if (n < 3) {
        return(data.frame(
          n_paired      = n,
          p_value_raw   = NA_real_,
          rank_biserial = NA_real_,
          cohen_d       = NA_real_
        ))
      }

      data.frame(
        n_paired      = n,
        p_value_raw   = wilcox.test(x, y, paired = TRUE)$p.value,
        rank_biserial = wilcoxonPairedR(x, y),
        cohen_d       = cohen.d(x, y, paired = TRUE)$estimate
      )
    }) %>%
    ungroup() %>%
    mutate(p_value_adj = p.adjust(p_value_raw, method = "fdr"))
}


# -----------------------------------------------------------------------------
# 4. Paired panel plot helper
# -----------------------------------------------------------------------------

build_alpha_panel <- function(df, metric, y_label, stats_df) {

  max_val     <- max(df[[metric]], na.rm = TRUE)
  y_pos_label <- max_val + (max_val * 0.12)
  y_upper     <- max_val + (max_val * 0.22)

  label_df <- stats_df %>%
    mutate(
      label      = sprintf("p = %.3f (n = %d)", p_value_raw, n_paired),
      y_position = y_pos_label
    )

  ggplot(df, aes(x = Timepoint, y = .data[[metric]])) +
    geom_boxplot(
      alpha         = 0.9,
      outlier.shape = NA,
      position      = position_dodge(width = 0.3)
    ) +
    geom_point(
      aes(group = Patient_ID, color = Patient_ID),
      position  = position_dodge(width = 0.3),
      size      = 2.5,
      alpha     = 1
    ) +
    geom_line(
      aes(group = Patient_ID, color = Patient_ID),
      position  = position_dodge(width = 0.3),
      alpha     = 0.5,
      linewidth = 1,
      linetype  = "dashed"
    ) +
    geom_text(
      data        = label_df,
      aes(x = 1.5, y = y_position, label = label),
      inherit.aes = FALSE,
      size        = 4
    ) +
    facet_grid(Niche_short ~ Method, scales = "free_y") +
    scale_color_manual(values = patient_colors) +
    scale_y_continuous(
      limits = c(0, y_upper),
      breaks = pretty(c(0, max_val), n = 4)
    ) +
    labs(x = "Timepoint", y = y_label) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position  = "none",
      axis.text.x      = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text       = element_text(face = "bold"),
      panel.border     = element_rect(color = "black", fill = NA,
                                      linewidth = 1.2),
      panel.grid.major = element_line(color = "gray80", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
}


# -----------------------------------------------------------------------------
# 5. Figure 3A — Shannon diversity
# -----------------------------------------------------------------------------

stats_shannon <- run_paired_wilcoxon(analysis_df, "Shannon")

write.csv(stats_shannon,
          file.path(data_dir, "stats_Shannon_paired_wilcoxon.csv"),
          row.names = FALSE)

cat("\n=== Shannon paired Wilcoxon results ===\n")
print(stats_shannon)

fig3a <- build_alpha_panel(analysis_df, "Shannon",
                           "Shannon Index", stats_shannon)
fig3a

ggsave(file.path(output_dir, "Fig3A_Shannon.svg"),
       fig3a, width = 14, height = 18, units = "cm")
ggsave(file.path(output_dir, "Fig3A_Shannon.pdf"),
       fig3a, width = 14, height = 18, units = "cm", device = cairo_pdf)


# -----------------------------------------------------------------------------
# 6. Figure 3B — Relative Dominance
# -----------------------------------------------------------------------------

stats_dominance <- run_paired_wilcoxon(analysis_df, "Relative_Dominance")

write.csv(stats_dominance,
          file.path(data_dir, "stats_RelativeDominance_paired_wilcoxon.csv"),
          row.names = FALSE)

cat("\n=== Relative Dominance paired Wilcoxon results ===\n")
print(stats_dominance)

fig3b <- build_alpha_panel(analysis_df, "Relative_Dominance",
                           "Relative Dominance", stats_dominance)
fig3b

ggsave(file.path(output_dir, "Fig3B_RelativeDominance.svg"),
       fig3b, width = 14, height = 18, units = "cm")
ggsave(file.path(output_dir, "Fig3B_RelativeDominance.pdf"),
       fig3b, width = 14, height = 18, units = "cm", device = cairo_pdf)

cat("\nFigure 3 panels saved to", output_dir, "\n")
