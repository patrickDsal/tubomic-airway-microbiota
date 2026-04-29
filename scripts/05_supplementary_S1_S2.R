# =============================================================================
# 05_supplementary_S1_S2.R
# Figure S1: Read depth QC (raincloud) + Mock community validation
# Figure S2: Sample availability across patients, timepoints, and niches
#
# Requires: 00_setup.R to have been run first
# Input:    data/ps_clean_bio_filtered.rds
#           data/sample_data_decontam.csv
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Libraries and shared theme
# -----------------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(phyloseq)
library(ggplot2)
library(ggdist)
library(patchwork)
library(scales)
library(conflicted)

conflict_prefer("filter", "dplyr")

source("00_setup.R")   # loads theme_pub() and data_dir

data_dir   <- "data/"
output_dir <- "output/figures/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# FIGURE S1A — Raincloud plot: reads per sample (log-scale)
# =============================================================================

ps <- readRDS(file.path(data_dir, "ps_clean_bio_filtered.rds"))

reads_per_sample <- data.frame(
  Sample    = names(sample_sums(ps)),
  Reads     = sample_sums(ps),
  log_reads = log10(sample_sums(ps) + 1)
)

fig_s1a <- ggplot(reads_per_sample, aes(x = log_reads)) +
  # Half-eye (density) above the axis
  ggdist::stat_halfeye(
    adjust       = 0.5,
    width        = 0.6,
    .width       = 0,
    justification = -0.3,
    point_colour = NA,
    fill         = "#B84B14"
  ) +
  # Compact boxplot
  geom_boxplot(
    width         = 0.25,
    outlier.shape = NA,
    color         = "black",
    fill          = "#B84B14",
    linewidth     = 1
  ) +
  # Raw data points (jittered along y only)
  geom_jitter(
    aes(y = 0),
    width  = 0.05,
    height = 0.1,
    size   = 1.5,
    alpha  = 0.5,
    color  = "#B84B14"
  ) +
  # Minimum read threshold
  geom_vline(
    xintercept = log10(2300),
    linetype   = "dashed",
    color      = "#B84B14",
    linewidth  = 1
  ) +
  scale_x_continuous(
    breaks = log10(c(1000, 2000, 5000, 10000, 50000, 100000)),
    labels = c("1k", "2k", "5k", "10k", "50k", "100k")
  ) +
  labs(x = "Reads per Sample", y = NULL) +
  theme_pub()

fig_s1a


# =============================================================================
# FIGURE S1C — Mock community: observed vs expected relative abundance
# =============================================================================

# ── Colour palette (consistent across all mock panels) ──────────────────────
mock_colors <- c(
  "Listeria"            = "#D73027",
  "Pseudomonas"         = "#4575B4",
  "Bacillus"            = "#66C2A5",
  "Escherichia"         = "#8C6BB1",
  "Salmonella"          = "#FDAE61",
  "Limosilactobacillus" = "#A6D96A",
  "Enterococcus"        = "#9E9AC8",
  "Staphylococcus"      = "#F47D4A"
)

italic_labels <- function(x) parse(text = paste0("italic('", x, "')"))

# ── Expected community composition ──────────────────────────────────────────
expected <- data.frame(
  SampleID  = "Expected",
  Genus     = names(mock_colors),
  Abundance = c(0.141, 0.042, 0.174, 0.101, 0.104, 0.184, 0.099, 0.155)
)

# ── Observed relative abundances from the positive control samples ───────────
# Pull mock community samples from the full (unfiltered) phyloseq
# NOTE: ps_raw is built inside 00_setup.R; if re-running standalone,
#       reload the raw phyloseq before this step.
ps_mock <- prune_samples(
  sample_data(ps_raw)$Environment == "positive_control", ps_raw
)
ps_mock_rel <- transform_sample_counts(ps_mock, function(x) x / sum(x))

# Melt to long format and keep only the 8 mock genera
abund_long <- psmelt(ps_mock_rel) %>%
  rename(Genus = Genus) %>%                 # already named Genus by psmelt
  filter(Genus %in% names(mock_colors)) %>%
  group_by(Sample, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  rename(SampleID = Sample)

# Combine observed + expected for the stacked bar
sample_order  <- c(sort(unique(abund_long$SampleID)), "Expected")
plot_data_bar <- bind_rows(abund_long, expected) %>%
  mutate(SampleID = factor(SampleID, levels = sample_order))

fig_s1c <- ggplot(plot_data_bar, aes(x = SampleID, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.75, colour = NA) +
  geom_vline(
    xintercept = length(sample_order) - 0.5,
    linewidth  = 0.4, linetype = "dashed", colour = "grey50"
  ) +
  scale_fill_manual(
    values = mock_colors,
    labels = italic_labels,
    name   = NULL
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02)),
    limits = c(0, 1.01)
  ) +
  labs(x = NULL, y = "Relative abundance") +
  theme_pub() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right"
  )

fig_s1c


# =============================================================================
# Assemble Figure S1 (A + C side by side)
# =============================================================================

fig_s1 <- fig_s1a + fig_s1c +
  plot_layout(ncol = 2, widths = c(1, 1.4)) +
  plot_annotation(tag_levels = "A")

fig_s1

ggsave(
  file.path(output_dir, "FigS1_QC.svg"),
  fig_s1, width = 30, height = 12, units = "cm"
)
ggsave(
  file.path(output_dir, "FigS1_QC.pdf"),
  fig_s1, width = 30, height = 12, units = "cm", device = cairo_pdf
)


# =============================================================================
# FIGURE S2 — Sample availability across patients, timepoints, and niches
# =============================================================================

sample_df <- read.csv(file.path(data_dir, "sample_data_decontam.csv"),
                      stringsAsFactors = FALSE)

# Define the expected niche × timepoint order for the x-axis
niche_order <- c(
  "T1_LLL", "T2_LLL",
  "T1_LRL", "T2_LRL",
  "T1_URL", "T2_URL",
  "T1_TS",  "T2_TS",
  "T1_T",   "T2_T",
  "T3_Tube"
)

# Keep only human samples and flag completeness per patient
plot_df <- sample_df %>%
  filter(Environment == "Human") %>%
  mutate(
    Patient_Num = as.numeric(str_extract(Patient_ID, "\\d+")),
    Niche_Time  = factor(Niche_Time, levels = niche_order)
  ) %>%
  arrange(Patient_Num) %>%
  mutate(Patient_ID = factor(Patient_ID, levels = rev(unique(Patient_ID))))

# Determine which patients have a complete sample set
full_set_status <- plot_df %>%
  group_by(Patient_ID) %>%
  summarise(
    Full_Set = ifelse(all(niche_order %in% Niche_Time), "Yes", "No"),
    .groups  = "drop"
  )

plot_df <- plot_df %>%
  left_join(full_set_status, by = "Patient_ID")

fig_s2 <- ggplot(plot_df,
                 aes(x = Niche_Time, y = Patient_ID, colour = Method)) +
  geom_point(size = 4, shape = 16, alpha = 0.8) +
  scale_colour_manual(values = c("Standard" = "#f0c571", "Venner" = "#f0746e")) +
  facet_wrap(~ Full_Set, scales = "free_y", ncol = 2,
             labeller = labeller(Full_Set = c("Yes" = "Yes", "No" = "No"))) +
  labs(
    x      = "Sample (Niche × Timepoint)",
    y      = "Patient",
    colour = "Intubation method"
  ) +
  theme_pub(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position  = "top"
  )

fig_s2

ggsave(
  file.path(output_dir, "FigS2_SampleAvailability.svg"),
  fig_s2, width = 24, height = 18, units = "cm"
)
ggsave(
  file.path(output_dir, "FigS2_SampleAvailability.pdf"),
  fig_s2, width = 24, height = 18, units = "cm", device = cairo_pdf
)

cat("\nFigures S1 and S2 saved to", output_dir, "\n")
