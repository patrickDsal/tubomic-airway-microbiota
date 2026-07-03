# =============================================================================
# 09_reviewer_reanalysis_beta.R
# Reviewer Comment 1 — Formal between-group comparison of T2−T1 changes
#
# The reviewer asked for a direct between-group test showing whether the
# temporal changes in niche-niche distances differed between Standard and
# Venner groups.
#
# Approach:
#   (a) Compute T2−T1 delta per patient per niche pair
#   (b) Compare deltas between ST and VT with unpaired Wilcoxon rank-sum
#   (c) Linear mixed-effects model: Distance ~ Time * Method + (1|Patient)
#       — p-value for the Time:Method interaction is the formal test
#
# Requires: 04_figure4_beta_network.R to have been run first
#           (produces niche_network_full.rds and niche_network_all.rds)
#
# Output:   data/delta_distances_FullSet.csv
#           data/delta_distances_ALL.csv
#           data/delta_wilcox_results.csv
#           data/lme_groupXtime_results.csv
#           data/rebuttal_beta_results_combined.csv
#           output/figures/FigS4_Delta_Comparison.svg/.pdf
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Libraries
# -----------------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(lme4)
library(lmerTest)
library(rcompanion)
library(effsize)
library(RColorBrewer)
library(patchwork)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflicts_prefer(lmerTest::lmer)

source("scripts/00_setup.R")

data_dir   <- "data/"
output_dir <- "output/figures/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

method_colors <- c(Standard = "#f0c571", Venner = "#f0746e")


# -----------------------------------------------------------------------------
# 2. Load pre-computed distance data from script 04
# -----------------------------------------------------------------------------

niche_network_full <- readRDS(file.path(data_dir, "niche_network_full.rds"))
niche_network_all  <- readRDS(file.path(data_dir, "niche_network_all.rds"))


# -----------------------------------------------------------------------------
# 3. Helper: prepare filtered long data
#    Same logic as script 04 run_beta_stats — faithful to original
# -----------------------------------------------------------------------------

prepare_long <- function(niche_net, method_sel) {

  niche_net %>%
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
    ) %>%
    filter(Time_Comparison %in% c("T1", "T2")) %>%
    group_by(Patient1, Niche_Pair) %>%
    filter(all(c("T1", "T2") %in% Time_Comparison)) %>%
    ungroup() %>%
    mutate(Method = method_sel)
}


# -----------------------------------------------------------------------------
# 4. Compute T2−T1 deltas per patient per niche pair
#    Pivot wide so each row is one patient — pairing guaranteed
# -----------------------------------------------------------------------------

compute_deltas <- function(niche_net) {

  # Combine Standard and Venner
  bind_rows(
    prepare_long(niche_net, "Standard"),
    prepare_long(niche_net, "Venner")
  ) %>%
    # Average duplicates within patient × niche pair × timepoint × method
    group_by(Patient1, Niche_Pair, Method, Time_Comparison) %>%
    summarise(Distance = mean(Distance, na.rm = TRUE), .groups = "drop") %>%
    # Pivot wide — one row per patient per niche pair
    pivot_wider(names_from  = Time_Comparison,
                values_from = Distance) %>%
    filter(!is.na(T1), !is.na(T2)) %>%
    mutate(Delta = T2 - T1)
}

deltas_full <- compute_deltas(niche_network_full)
deltas_all  <- compute_deltas(niche_network_all)

write.csv(deltas_full,
          file.path(data_dir, "delta_distances_FullSet.csv"),
          row.names = FALSE)
write.csv(deltas_all,
          file.path(data_dir, "delta_distances_ALL.csv"),
          row.names = FALSE)

cat("\n--- Delta counts per method per niche pair (full set) ---\n")
print(deltas_full %>%
        count(Method, Niche_Pair) %>%
        pivot_wider(names_from = Method, values_from = n))


# =============================================================================
# PART A — WILCOXON RANK-SUM ON DELTAS (ST vs VT per niche pair)
# =============================================================================

run_delta_wilcox <- function(deltas, label = "") {

  deltas %>%
    group_by(Niche_Pair) %>%
    filter(n_distinct(Method) == 2) %>%
    summarise(
      n_ST          = sum(Method == "Standard"),
      n_VT          = sum(Method == "Venner"),
      mean_delta_ST = mean(Delta[Method == "Standard"], na.rm = TRUE),
      sd_delta_ST   = sd(Delta[Method == "Standard"],   na.rm = TRUE),
      mean_delta_VT = mean(Delta[Method == "Venner"],   na.rm = TRUE),
      sd_delta_VT   = sd(Delta[Method == "Venner"],     na.rm = TRUE),
      p_value       = wilcox.test(
        Delta[Method == "Standard"],
        Delta[Method == "Venner"],
        exact = FALSE
      )$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      p_adj      = p.adjust(p_value, method = "BH"),
      sig        = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      ),
      sample_set = label
    )
}

delta_wilcox_full <- run_delta_wilcox(deltas_full, "FullSet")
delta_wilcox_all  <- run_delta_wilcox(deltas_all,  "AllAvailable")

write.csv(bind_rows(delta_wilcox_full, delta_wilcox_all),
          file.path(data_dir, "delta_wilcox_results.csv"),
          row.names = FALSE)

cat("\n=== Wilcoxon on deltas — full set ===\n")
print(delta_wilcox_full %>%
        select(Niche_Pair, n_ST, n_VT,
               mean_delta_ST, mean_delta_VT, p_value, p_adj, sig))

cat("\n=== Wilcoxon on deltas — all available ===\n")
print(delta_wilcox_all %>%
        select(Niche_Pair, n_ST, n_VT,
               mean_delta_ST, mean_delta_VT, p_value, p_adj, sig))


# =============================================================================
# PART B — LINEAR MIXED-EFFECTS MODEL: Distance ~ Time * Method + (1|Patient)
# =============================================================================

prepare_lme_data <- function(niche_net) {

  bind_rows(
    prepare_long(niche_net, "Standard"),
    prepare_long(niche_net, "Venner")
  ) %>%
    group_by(Patient1, Niche_Pair, Method, Time_Comparison) %>%
    summarise(Distance = mean(Distance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Patient1, Niche_Pair) %>%
    filter(all(c("T1", "T2") %in% Time_Comparison)) %>%
    ungroup() %>%
    mutate(
      Method         = factor(Method, levels = c("Standard", "Venner")),
      Time_Comparison = factor(Time_Comparison, levels = c("T1", "T2")),
      Patient1       = factor(Patient1)
    )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

run_lme_per_pair <- function(data, label = "") {

  map_dfr(unique(data$Niche_Pair), function(np) {

    df <- data %>% filter(Niche_Pair == np)

    if (nlevels(droplevels(df$Method)) < 2 ||
        nlevels(droplevels(df$Time_Comparison)) < 2 ||
        nrow(df) < 8) {
      return(data.frame(
        Niche_Pair       = np,
        n_obs            = nrow(df),
        n_patients       = n_distinct(df$Patient1),
        beta_interaction = NA_real_,
        p_interaction    = NA_real_,
        sample_set       = label
      ))
    }

    m <- tryCatch(
      lmer(Distance ~ Time_Comparison * Method + (1 | Patient1),
           data = df, REML = FALSE),
      error   = function(e) NULL,
      warning = function(w) suppressWarnings(
        lmer(Distance ~ Time_Comparison * Method + (1 | Patient1),
             data = df, REML = FALSE)
      )
    )

    if (is.null(m)) {
      return(data.frame(
        Niche_Pair       = np,
        n_obs            = nrow(df),
        n_patients       = n_distinct(df$Patient1),
        beta_interaction = NA_real_,
        p_interaction    = NA_real_,
        sample_set       = label
      ))
    }

    coef_tbl <- as.data.frame(summary(m)$coefficients)
    int_row  <- grep("Time_Comparison.*Method|Method.*Time_Comparison",
                     rownames(coef_tbl), value = TRUE)

    data.frame(
      Niche_Pair       = np,
      n_obs            = nrow(df),
      n_patients       = n_distinct(df$Patient1),
      beta_interaction = if (length(int_row)) coef_tbl[int_row, "Estimate"] else NA_real_,
      p_interaction    = if (length(int_row)) coef_tbl[int_row, "Pr(>|t|)"] else NA_real_,
      sample_set       = label
    )
  }) %>%
    mutate(
      p_adj = p.adjust(p_interaction, method = "BH"),
      sig   = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    )
}

lme_data_full <- prepare_lme_data(niche_network_full)
lme_data_all  <- prepare_lme_data(niche_network_all)

lme_full <- run_lme_per_pair(lme_data_full, "FullSet")
lme_all  <- run_lme_per_pair(lme_data_all,  "AllAvailable")

write.csv(bind_rows(lme_full, lme_all),
          file.path(data_dir, "lme_groupXtime_results.csv"),
          row.names = FALSE)

cat("\n=== LME group × time interaction — full set ===\n")
print(lme_full %>% select(Niche_Pair, n_obs, n_patients,
                          beta_interaction, p_interaction, p_adj, sig))

cat("\n=== LME group × time interaction — all available ===\n")
print(lme_all %>% select(Niche_Pair, n_obs, n_patients,
                         beta_interaction, p_interaction, p_adj, sig))


# =============================================================================
# PART C — FIGURE FOR REBUTTAL
# =============================================================================

label_delta <- lme_full %>%
  select(Niche_Pair, p_interaction) %>%
  left_join(
    deltas_full %>%
      group_by(Niche_Pair) %>%
      summarise(y_pos = max(abs(Delta), na.rm = TRUE) + 0.1,
                .groups = "drop"),
    by = "Niche_Pair"
  ) %>%
  mutate(label = sprintf("p = %.3f", p_interaction))

pt_colors <- colorRampPalette(brewer.pal(12, "Set3"))(
  n_distinct(deltas_full$Patient1)
) %>% setNames(sort(unique(deltas_full$Patient1)))

figR1 <- ggplot(deltas_full, aes(x = Method, y = Delta, fill = Method)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = Patient1), width = 0.15,
              size = 3, alpha = 0.85) +
  geom_text(data = label_delta,
            aes(x = 1.5, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  facet_wrap(~ Niche_Pair, scales = "free_y", ncol = 4) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  scale_fill_manual(values  = method_colors,
                    labels  = c("Standard" = "Standard", "Venner" = "Venner")) +
  scale_color_manual(values = pt_colors) +
  scale_x_discrete(labels = c("Standard" = "Standard", "Venner" = "Venner")) +
  labs(
    x     = NULL,
    y     = "Delta Morisita-Horn distance (T2 - T1)",
    title = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(face = "bold", size = 14),
    axis.text.y      = element_text(size = 12),
    axis.title.y     = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text       = element_text(face = "bold", size = 13),
    panel.border     = element_rect(color = "black", fill = NA,
                                    linewidth = 1.1),
    panel.grid.major = element_line(color = "grey85", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5)
  )

figR1

ggsave(file.path(output_dir, "FigS4_Delta_Comparison.svg"),
       figR1, width = 28, height = 18, units = "cm")
ggsave(file.path(output_dir, "FigS4_Delta_Comparison.pdf"),
       figR1, width = 28, height = 18, units = "cm", device = cairo_pdf)


# =============================================================================
# PART D — COMBINED SUMMARY TABLE FOR REBUTTAL LETTER
# =============================================================================

summary_rebuttal <- delta_wilcox_full %>%
  select(Niche_Pair, n_ST, n_VT,
         mean_delta_ST, sd_delta_ST,
         mean_delta_VT, sd_delta_VT,
         p_wilcox_full   = p_value,
         q_wilcox_full   = p_adj,
         sig_wilcox_full = sig) %>%
  left_join(
    delta_wilcox_all %>% select(Niche_Pair,
                                p_wilcox_all   = p_value,
                                q_wilcox_all   = p_adj,
                                sig_wilcox_all = sig),
    by = "Niche_Pair"
  ) %>%
  left_join(
    lme_full %>% select(Niche_Pair,
                        beta_lme_full = beta_interaction,
                        p_lme_full    = p_interaction,
                        q_lme_full    = p_adj,
                        sig_lme_full  = sig),
    by = "Niche_Pair"
  ) %>%
  left_join(
    lme_all %>% select(Niche_Pair,
                       beta_lme_all = beta_interaction,
                       p_lme_all    = p_interaction,
                       q_lme_all    = p_adj,
                       sig_lme_all  = sig),
    by = "Niche_Pair"
  )

write.csv(summary_rebuttal,
          file.path(data_dir, "rebuttal_beta_results_combined.csv"),
          row.names = FALSE)

cat("\n=== COMBINED REBUTTAL RESULTS ===\n")
print(summary_rebuttal)

cat("\nScript 09 complete.\n")

# Create clean table
library(flextable)
library(officer)

table_s1 <- summary_rebuttal %>%
  mutate(
    `Niche Pair`           = Niche_Pair,
    `n ST`                 = n_ST,
    `n VT`                 = n_VT,
    `Delta ST (mean ± SD)` = sprintf("%.2f ± %.2f", mean_delta_ST, sd_delta_ST),
    `Delta VT (mean ± SD)` = sprintf("%.2f ± %.2f", mean_delta_VT, sd_delta_VT),
    `p`                    = sprintf("%.3f", p_lme_full),
    `q`                    = sprintf("%.3f", q_lme_full),
    `Sig`                  = sig_lme_full
  ) %>%
  select(`Niche Pair`, `n ST`, `n VT`,
         `Delta ST (mean ± SD)`, `Delta VT (mean ± SD)`,
         `p`, `q`, `Sig`)

ft <- flextable(table_s1) %>%
  bold(part = "header") %>%
  border_remove() %>%
  hline_top(part = "header", 
            border = fp_border(width = 1.5)) %>%
  hline_bottom(part = "header", 
               border = fp_border(width = 1)) %>%
  hline_bottom(part = "body", 
               border = fp_border(width = 1.5)) %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Calibri", part = "all") %>%
  autofit()

doc <- read_docx() %>%
  body_add_par("Table S1: Between-group comparison of T2−T1 changes in niche-pair Morisita–Horn distances. Negative delta values indicate convergence; positive values indicate divergence. P-values from the Timepoint × Method interaction term (Distance ~ Timepoint × Method + (1|Patient)) with Benjamini–Hochberg FDR correction. T = Throat, TS = Trachea, URL = Upper Right Lobe, LLL = Lower Left Lobe, LRL = Lower Right Lobe.", 
               style = "Normal") %>%
  body_add_flextable(ft)

print(doc, target = file.path(output_dir, "TableS1_Delta_Summary.docx"))
