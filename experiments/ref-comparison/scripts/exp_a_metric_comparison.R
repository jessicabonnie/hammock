#!/usr/bin/env Rscript
# =============================================================================
# scripts/exp_a_metric_comparison.R
#
# Sanity-check the new Mode D columns (containment + cosketch) on the
# cross-reference robustness question at one (k, w) cell. For each of the
# 12 columns (6 minimizer-only + 6 with-ends), compute:
#   - median(cross_ref) − median(diff_tissue)
#   - Wilcoxon p (one-sided, cross_ref > diff_tissue)
# and plot effect size + −log10(p), grouped by sketch flavor.
#
# Usage:
#   Rscript scripts/exp_a_metric_comparison.R <csv> <metadata_tsv> <out_png> <out_tsv> [peak_type_label]
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(tibble)
  library(purrr); library(ggplot2); library(cowplot); library(Cairo)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4)
csv_file <- args[1]; meta_file <- args[2]
out_png  <- args[3]; out_tsv   <- args[4]
peak_lbl <- if (length(args) >= 5) args[5] else "broad"

metric_groups <- tribble(
  ~metric,                          ~flavor,           ~kind,
  "jaccard_similarity",             "minimizer",       "jaccard",
  "containment_AB",                 "minimizer",       "containment_AB",
  "containment_BA",                 "minimizer",       "containment_BA",
  "cosketch_geom",                  "minimizer",       "cosketch_geom",
  "cosketch_arith",                 "minimizer",       "cosketch_arith",
  "cosketch_max",                   "minimizer",       "cosketch_max",
  "jaccard_similarity_with_ends",   "minimizer+ends",  "jaccard",
  "containment_AB_with_ends",       "minimizer+ends",  "containment_AB",
  "containment_BA_with_ends",       "minimizer+ends",  "containment_BA",
  "cosketch_geom_with_ends",        "minimizer+ends",  "cosketch_geom",
  "cosketch_arith_with_ends",       "minimizer+ends",  "cosketch_arith",
  "cosketch_max_with_ends",         "minimizer+ends",  "cosketch_max"
)

mat  <- read_csv(csv_file,  show_col_types = FALSE)
meta <- read_tsv(meta_file, show_col_types = FALSE)

pairs <- mat %>%
  mutate(sample_a = str_remove(basename(file1), "\\.fa$"),
         ref_a    = basename(dirname(file1)),
         sample_b = str_remove(basename(file2), "\\.fa$"),
         ref_b    = basename(dirname(file2))) %>%
  filter(file1 != file2) %>%
  left_join(meta %>% select(sample_id, tissue) %>% distinct(),
            by = c("sample_a" = "sample_id")) %>%
  rename(tissue_a = tissue) %>%
  left_join(meta %>% select(sample_id, tissue) %>% distinct(),
            by = c("sample_b" = "sample_id")) %>%
  rename(tissue_b = tissue) %>%
  mutate(pair_type = case_when(
    tissue_a == tissue_b & ref_a != ref_b ~ "cross_ref",
    tissue_a != tissue_b                  ~ "diff_tissue",
    TRUE                                  ~ "other"))

compute_row <- function(metric, flavor, kind) {
  xref <- pairs %>% filter(pair_type == "cross_ref")   %>% pull(.data[[metric]])
  ddif <- pairs %>% filter(pair_type == "diff_tissue") %>% pull(.data[[metric]])
  test <- tryCatch(
    wilcox.test(xref, ddif, alternative = "greater"),
    error = function(e) list(p.value = NA_real_))
  tibble(metric        = metric,
         flavor        = flavor,
         kind          = kind,
         median_xref   = median(xref),
         median_diff   = median(ddif),
         delta         = median(xref) - median(ddif),
         wilcoxon_p    = test$p.value,
         neg_log10_p   = -log10(test$p.value))
}
stats <- purrr::pmap_dfr(metric_groups, compute_row)

write_tsv(stats, out_tsv)

# Plot 1: delta (effect size)
kind_order <- c("jaccard", "containment_AB", "containment_BA",
                "cosketch_geom", "cosketch_arith", "cosketch_max")
stats <- stats %>%
  mutate(kind   = factor(kind,   levels = kind_order),
         flavor = factor(flavor, levels = c("minimizer", "minimizer+ends")))

p_delta <- ggplot(stats, aes(x = kind, y = delta, fill = flavor)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", delta)),
            position = position_dodge(0.8), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c(minimizer = "#1f77b4", `minimizer+ends` = "#ff7f0e")) +
  labs(title = "Effect size: median(cross-ref) − median(diff-tissue)",
       x = NULL, y = "Δ median",
       fill = "Sketch") +
  theme_cowplot(11) +
  theme(axis.text.x  = element_text(angle = 25, hjust = 1, size = 9),
        plot.title   = element_text(size = 11, face = "bold"),
        legend.position = "top") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)))

p_p <- ggplot(stats, aes(x = kind, y = neg_log10_p, fill = flavor)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = formatC(wilcoxon_p, format = "e", digits = 1)),
            position = position_dodge(0.8), vjust = -0.4, size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c(minimizer = "#1f77b4", `minimizer+ends` = "#ff7f0e"),
                    guide = "none") +
  labs(title = "Discriminative power: −log₁₀(Wilcoxon p)",
       x = NULL, y = "−log₁₀(p)") +
  theme_cowplot(11) +
  theme(axis.text.x  = element_text(angle = 25, hjust = 1, size = 9),
        plot.title   = element_text(size = 11, face = "bold")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)))

title <- ggdraw() + draw_label(
  sprintf("Experiment A: discrimination across all Mode D metrics (k=10, w=10, %s)", peak_lbl),
  fontface = "bold", x = 0.02, hjust = 0, size = 13)
subtitle <- ggdraw() + draw_label(
  paste0("Cross-ref-of-same-tissue (n=", sum(pairs$pair_type=="cross_ref"),
         ") vs different-tissue (n=", sum(pairs$pair_type=="diff_tissue"), ") pairs."),
  x = 0.02, hjust = 0, size = 10, color = "grey30")

body <- plot_grid(p_delta, p_p, ncol = 1, align = "v", labels = c("A", "B"))
out  <- plot_grid(title, subtitle, body, ncol = 1,
                  rel_heights = c(0.04, 0.03, 1))

CairoPNG(out_png, width = 11, height = 8, units = "in", res = 300)
print(out)
invisible(dev.off())
message("Wrote ", out_png, " and ", out_tsv)
