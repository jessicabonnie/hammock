#!/usr/bin/env Rscript
# =============================================================================
# scripts/exp_a_metric_comparison.R
#
# Sanity-check the new Mode D columns (containment + cosketch) on the
# cross-reference robustness question at one (k, w) cell. For each similarity
# column present in the input, compute:
#   - median(cross_ref) − median(diff_tissue)
#   - Wilcoxon p (one-sided, cross_ref > diff_tissue)
# and plot effect size + −log10(p), grouped by sketch flavor.
#
# Schema note: hammock < 0.5.0 emitted 12 such columns (6 minimizer-only +
# 6 with-ends); 0.5.0 added jaccard_similarity_ie and its _with_ends twin, for
# 14; 0.6.0 REMOVED all 7 `_with_ends` columns (CLAUDE.md divergence #8),
# leaving 7. The table below lists all 14 and the script keeps only those
# actually present, so it runs unchanged on CSVs of any vintage -- rather than
# hard-erroring inside pull(.data[[metric]]), which is what a missing column
# used to do. It prints which ones it dropped. On a v0.6.0 CSV the
# "minimizer+ends" flavor simply drops out of the figure.
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

mat  <- read_csv(csv_file,  show_col_types = FALSE)

# Resolve the live register-equality column name before building the
# tribble, so both the read (via metric_groups$metric, below) and the
# re-emitted `metric` value in out_tsv use the resolved name. Prefer the
# post-rename `reg_eq_similarity`; fall back to the legacy
# `jaccard_similarity` for archived CSVs written before
# docs/seed-jaccard-reg-eq-rename.md's Step 1 -- logged only when the
# fallback actually fires (i.e. `reg_eq_similarity` is genuinely absent),
# kept as its own message distinct from the generic "skipping" note below so
# the two don't read as duplicates of each other.
reg_eq_col <- if ("reg_eq_similarity" %in% names(mat)) "reg_eq_similarity" else "jaccard_similarity"
if (!("reg_eq_similarity" %in% names(mat))) {
  message("exp_a_metric_comparison.R: 'reg_eq_similarity' column not found ",
          "in ", basename(csv_file), "; falling back to legacy ",
          "'jaccard_similarity' column name (pre-rename hammock output).")
}

metric_groups <- tribble(
  ~metric,                          ~flavor,           ~kind,
  reg_eq_col,                       "minimizer",       "jaccard",
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
  "cosketch_max_with_ends",         "minimizer+ends",  "cosketch_max",
  # hammock >= 0.5.0 only; absent from CSVs written before 2026-07-31.
  "jaccard_similarity_ie",          "minimizer",       "jaccard_ie",
  "jaccard_similarity_ie_with_ends","minimizer+ends",  "jaccard_ie"
)

missing_metrics <- setdiff(metric_groups$metric, names(mat))
if (length(missing_metrics) > 0) {
  message("note: ", length(missing_metrics), " metric column(s) absent from ",
          basename(csv_file), " -- skipping: ",
          paste(missing_metrics, collapse = ", "))
  metric_groups <- metric_groups %>% filter(!metric %in% missing_metrics)
}
if (nrow(metric_groups) == 0) stop("no known similarity columns in ", csv_file)
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
# Must list every `kind` in metric_groups: factor() maps an unlisted level to
# NA and ggplot then drops those bars with only a warning.
kind_order <- c("jaccard", "jaccard_ie", "containment_AB", "containment_BA",
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
