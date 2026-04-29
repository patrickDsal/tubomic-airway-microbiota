# =============================================================================
# 01_figure2_clinical_scores.R
# Figure 2: Clinical severity scores over time (APACHE II, SAPS II, SOFA)
#
# This script uses CLINICAL data only — it does NOT require the phyloseq object.
# Input:  data/TUBOMIC_Scores_fullset.csv
#         Columns: Patient_ID, Group (ST/VT),
#                  APACHEScoreT1/T2/T3, SAPS2ScoreT1/T2/T3, SOFAScoreT1/T2/T3
# Output: output/figures/Fig2_Clinical_Scores.svg/.pdf
#         output/figures/Fig2A_APACHE.svg
#         output/figures/Fig2B_SAPS.svg
#         output/figures/Fig2C_SOFA.svg
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Libraries
# -----------------------------------------------------------------------------

library(tidyverse)
library(lmerTest)   # LMM with Satterthwaite p-values
library(patchwork)
library(ggtext)     # element_markdown() for italic p in subtitles
library(svglite)


# -----------------------------------------------------------------------------
# 2. Paths
# -----------------------------------------------------------------------------

data_dir   <- "data/"
output_dir <- "output/figures/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 3. Appearance settings
#    Edit these to adjust the look without touching the plot code
# -----------------------------------------------------------------------------

color_ST     <- "#f0c571"   # Standard endotracheal tube
color_VT     <- "#f0746e"   # Venner PneuX tube
point_size   <- 2.2
line_alpha   <- 0.15        # individual patient lines (semi-transparent)
mean_size    <- 3.5         # group mean point size
jitter_width <- 0.08        # horizontal separation between ST and VT
base_font    <- 11


# -----------------------------------------------------------------------------
# 4. Load and reshape data
# -----------------------------------------------------------------------------

raw <- read.csv(file.path(data_dir, "TUBOMIC_Scores_fullset.csv"),
                stringsAsFactors = FALSE)

# Pivot to long format: one row per patient × timepoint × score
long <- raw %>%
  pivot_longer(
    cols         = -c(Patient_ID, Group),
    names_to     = c("Score", "Timepoint"),
    names_pattern = "(APACHE|SAPS2|SOFA)Score(T[123])",
    values_to    = "Value"
  ) %>%
  mutate(
    Score      = factor(Score,     levels = c("APACHE", "SAPS2", "SOFA")),
    Timepoint  = factor(Timepoint, levels = c("T1", "T2", "T3")),
    Group      = factor(Group,     levels = c("ST", "VT")),
    Patient_ID = factor(Patient_ID)
  )


# -----------------------------------------------------------------------------
# 5. Linear mixed models
#    Model: Value ~ Timepoint * Group + (1 | Patient_ID)
#    Extracts p(time), p(group), p(group × time) for each score
# -----------------------------------------------------------------------------

run_lmm <- function(df) {
  m   <- lmer(Value ~ Timepoint * Group + (1 | Patient_ID), data = df)
  aov <- anova(m)
  list(
    p_time  = aov["Timepoint",       "Pr(>F)"],
    p_group = aov["Group",           "Pr(>F)"],
    p_inter = aov["Timepoint:Group", "Pr(>F)"]
  )
}

conflicts_prefer(purrr::set_names)

stats_list <- long %>%
  group_by(Score) %>%
  group_split() %>%
  set_names(levels(long$Score)) %>%
  map(run_lmm)

# Print to console for quick inspection
fmt_p <- function(p) {
  if (is.na(p))   return("NA")
  if (p < 0.001)  return("< 0.001")
  sprintf("= %.3f", p)
}

cat("\n=== LMM results ===\n")
for (s in names(stats_list)) {
  cat(sprintf("%-8s  p(time) = %-8s  p(group) = %-8s  p(group*time) = %s\n",
              s,
              fmt_p(stats_list[[s]]$p_time),
              fmt_p(stats_list[[s]]$p_group),
              fmt_p(stats_list[[s]]$p_inter)))
}


# -----------------------------------------------------------------------------
# 6. Plot helper
# -----------------------------------------------------------------------------

make_subtitle <- function(stats) {
  paste0(
    "<i>p</i> (time) ",          fmt_p(stats$p_time),  "<br>",
    "<i>p</i> (group) ",         fmt_p(stats$p_group), "<br>",
    "<i>p</i> (group × time) ",  fmt_p(stats$p_inter)
  )
}

make_panel <- function(score_name, y_label, title_text, stats) {

  dat <- long %>% filter(Score == score_name)

  # Per-group summary (mean ± SE)
  smry <- dat %>%
    group_by(Group, Timepoint) %>%
    summarise(mean = mean(Value, na.rm = TRUE),
              se   = sd(Value,   na.rm = TRUE) / sqrt(sum(!is.na(Value))),
              .groups = "drop")

  # Horizontal dodge: ST left, VT right of each timepoint tick
  dat  <- dat  %>%
    mutate(x = as.numeric(Timepoint) +
             ifelse(Group == "ST", -jitter_width, jitter_width))
  smry <- smry %>%
    mutate(x = as.numeric(Timepoint) +
             ifelse(Group == "ST", -jitter_width, jitter_width))

  ggplot() +
    # Individual patient trajectories (thin, semi-transparent)
    geom_line(data  = dat,
              aes(x = x, y = Value, group = Patient_ID, color = Group),
              alpha = line_alpha, linewidth = 0.4) +
    geom_point(data = dat,
               aes(x = x, y = Value, color = Group),
               alpha = line_alpha + 0.15, size = point_size - 0.8) +
    # Group means ± SE (bold, filled)
    geom_errorbar(data = smry,
                  aes(x = x, ymin = mean - se, ymax = mean + se,
                      color = Group),
                  width = 0.06, linewidth = 0.7) +
    geom_line(data = smry,
              aes(x = x, y = mean, color = Group, group = Group),
              linewidth = 1.0) +
    geom_point(data = smry,
               aes(x = x, y = mean, fill = Group),
               shape = 21, color = "black",
               size = mean_size, stroke = 0.6) +
    scale_x_continuous(breaks = 1:3, labels = c("T1", "T2", "T3"),
                       limits = c(0.7, 3.45)) +
    scale_color_manual(values = c(ST = color_ST, VT = color_VT),
                       labels = c("Standard", "Venner PneuX")) +
    scale_fill_manual (values = c(ST = color_ST, VT = color_VT),
                       labels = c("Standard", "Venner PneuX")) +
    labs(title    = title_text,
         subtitle = make_subtitle(stats),
         x        = NULL,
         y        = y_label,
         color    = NULL,
         fill     = NULL) +
    theme_classic(base_size = base_font) +
    theme(
      plot.title    = element_text(face = "bold", size = base_font + 1,
                                   hjust = 0, margin = margin(b = 6)),
      plot.subtitle = element_markdown(size = base_font - 2,
                                       lineheight = 1.2,
                                       margin = margin(b = 6)),
      axis.title.y  = element_text(size = base_font, margin = margin(r = 6)),
      axis.text     = element_text(color = "black"),
      legend.position = "bottom",
      legend.key.size = unit(0.6, "lines"),
      plot.margin   = margin(8, 10, 5, 8)
    )
}


# -----------------------------------------------------------------------------
# 7. Build panels
# -----------------------------------------------------------------------------

p_apache <- make_panel("APACHE", "APACHE II score", "A  APACHE II",
                       stats_list$APACHE)
p_saps   <- make_panel("SAPS2",  "SAPS II score",   "B  SAPS II",
                       stats_list$SAPS2)
p_sofa   <- make_panel("SOFA",   "SOFA score",      "C  SOFA",
                       stats_list$SOFA)


# -----------------------------------------------------------------------------
# 8. Save individual panels and combined figure
# -----------------------------------------------------------------------------

ggsave(file.path(output_dir, "Fig2A_APACHE.svg"),
       p_apache, width = 3.5, height = 3.6, units = "in", device = svglite)
ggsave(file.path(output_dir, "Fig2B_SAPS.svg"),
       p_saps,   width = 3.5, height = 3.6, units = "in", device = svglite)
ggsave(file.path(output_dir, "Fig2C_SOFA.svg"),
       p_sofa,   width = 3.5, height = 3.6, units = "in", device = svglite)

fig2 <- (p_apache | p_saps | p_sofa) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "Fig2_Clinical_Scores.svg"),
       fig2, width = 10.5, height = 3.8, units = "in", device = svglite)
ggsave(file.path(output_dir, "Fig2_Clinical_Scores.pdf"),
       fig2, width = 10.5, height = 3.8, units = "in")

cat("\nFigure 2 saved to", output_dir, "\n")
