#!/usr/bin/env Rscript
# =============================================================================
# scripts/exp_a_validate_plot.R
# Experiment A: Cross-Reference Robustness Validation
#
# Input:  similarity_matrix.tsv  (all-vs-all pairwise Jaccard similarities)
#         exp_a_metadata.tsv     (sample_id, tissue, ref, mark)
# Output: cross_ref_validation.png  (boxplot + individual comparisons)
#         cross_ref_stats.tsv       (Wilcoxon test results)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggbeeswarm)
  library(cowplot)
})

# ── Snakemake I/O ─────────────────────────────────────────────────────────
mat_file  <- snakemake@input[["matrix"]]
meta_file <- snakemake@input[["metadata"]]
out_plot  <- snakemake@output[["plot"]]
out_stats <- snakemake@output[["stats"]]

# ── Load data ─────────────────────────────────────────────────────────────
mat  <- read_tsv(mat_file,  show_col_types = FALSE)
meta <- read_tsv(meta_file, show_col_types = FALSE)

# Convert wide matrix to long pairs
mat_long <- mat %>%
  pivot_longer(-sample_a, names_to = "sample_b", values_to = "similarity") %>%
  filter(sample_a != sample_b)   # remove self-comparisons

# Annotate pairs with metadata
pairs <- mat_long %>%
  left_join(meta, by = c("sample_a" = "sample_id")) %>%
  rename(tissue_a = tissue, ref_a = ref) %>%
  left_join(meta, by = c("sample_b" = "sample_id")) %>%
  rename(tissue_b = tissue, ref_b = ref)

# Classify each pair
pairs <- pairs %>%
  mutate(pair_type = case_when(
    tissue_a == tissue_b & ref_a != ref_b ~ "Same tissue, cross-reference",
    tissue_a == tissue_b & ref_a == ref_b ~ "Same tissue, same reference",
    tissue_a != tissue_b                  ~ "Different tissue",
    TRUE                                  ~ "Other"
  ))

# ── Statistical test ──────────────────────────────────────────────────────
cross_ref    <- pairs %>% filter(pair_type == "Same tissue, cross-reference") %>% pull(similarity)
diff_tissue  <- pairs %>% filter(pair_type == "Different tissue")             %>% pull(similarity)

test_result <- wilcox.test(cross_ref, diff_tissue, alternative = "greater")

stats_df <- tibble(
  comparison            = "Same-tissue cross-reference vs Different-tissue",
  n_cross_ref           = length(cross_ref),
  n_diff_tissue         = length(diff_tissue),
  median_cross_ref      = median(cross_ref),
  median_diff_tissue    = median(diff_tissue),
  wilcoxon_p            = test_result$p.value,
  wilcoxon_statistic    = test_result$statistic,
)
write_tsv(stats_df, out_stats)

# ── Plot ──────────────────────────────────────────────────────────────────
pair_levels <- c(
  "Same tissue, same reference",
  "Same tissue, cross-reference",
  "Different tissue"
)

plot_data <- pairs %>%
  filter(pair_type %in% pair_levels) %>%
  mutate(pair_type = factor(pair_type, levels = pair_levels))

p <- ggplot(plot_data, aes(x = pair_type, y = similarity, fill = pair_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_beeswarm(size = 1.5, alpha = 0.6, cex = 1.5) +
  scale_fill_manual(values = c("#2166ac", "#f4a582", "#d1e5f0")) +
  labs(
    title    = "Experiment A: Cross-Reference Robustness",
    subtitle = sprintf("Same-tissue cross-ref vs different-tissue: Wilcoxon p = %.2e",
                       test_result$p.value),
    x        = "Comparison type",
    y        = "Minimizer sketch similarity (Jaccard)",
    fill     = NULL
  ) +
  theme_cowplot(12) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 15, hjust = 1)
  )

ggsave(out_plot, p, width = 8, height = 6, dpi = 300)
message("Done. Plot saved to ", out_plot)

