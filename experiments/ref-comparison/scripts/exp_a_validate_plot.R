#!/usr/bin/env Rscript
# =============================================================================
# scripts/exp_a_validate_plot.R
# Experiment A: Cross-Reference Robustness Validation
#
# Input:  hammock similarity CSV (all-vs-all pairwise Jaccard from --full-paths)
#         exp_a_metadata.tsv     (sample_id, tissue, ref, mark)
# Output: cross_ref_validation.png  (boxplot + individual comparisons)
#         cross_ref_stats.tsv       (Wilcoxon test results)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(Cairo)
})

# Wrapper so ggsave's width/height/dpi (inches) reach Cairo::CairoPNG correctly
cairo_png_in <- function(filename, width, height, ...) {
  Cairo::CairoPNG(filename, width = width, height = height,
                  units = "in", res = 300, ...)
}

# ── Snakemake I/O ─────────────────────────────────────────────────────────
mat_file  <- snakemake@input[["matrix"]]
meta_file <- snakemake@input[["metadata"]]
out_plot  <- snakemake@output[["plot"]]
out_stats <- snakemake@output[["stats"]]
# Snakemake-injected literal, now only a legacy default -- see the
# reg_eq_similarity-preferred resolution below, after the data is loaded
# (docs/seed-jaccard-reg-eq-rename.md Step 2, round 2 correction:
# workflow/Snakefile:180's literal is left untouched and unused for the
# preferred path).
snakemake_sim_col <- snakemake@params[["sim_col"]]

# ── Load data ─────────────────────────────────────────────────────────────
mat  <- read_csv(mat_file, show_col_types = FALSE)
meta <- read_tsv(meta_file, show_col_types = FALSE)

# Prefer the post-rename `reg_eq_similarity` column if the loaded matrix has
# it; otherwise fall back to the Snakemake-passed `sim_col` param (which for
# pre-rename runs/configs is "jaccard_similarity"). Logged only when the
# fallback actually fires, so a fresh post-rename run stays silent.
if ("reg_eq_similarity" %in% names(mat)) {
  sim_col <- "reg_eq_similarity"
} else {
  sim_col <- snakemake_sim_col
  message("exp_a_validate_plot.R: 'reg_eq_similarity' column not found in ",
          mat_file, "; falling back to Snakemake-passed sim_col param (",
          snakemake_sim_col, ").")
}

# Hammock output is long-form with columns like:
#   file1, file2, jaccard_similarity, jaccard_similarity_with_ends, ...
# Extract sample_id and ref from file paths:
#   .../fastas/{peak_type}/{ref}/{sample}.fa → sample, ref
extract_info <- function(path) {
  parts <- str_split(basename(path), "\\.")[[1]]
  sample <- parts[1]
  # ref is the parent directory of the file
  ref <- basename(dirname(path))
  tibble(sample_id = sample, ref = ref)
}

pairs <- mat %>%
  rename(similarity = !!sym(sim_col)) %>%
  mutate(
    sample_a = str_remove(basename(file1), "\\.fa$"),
    ref_a    = basename(dirname(file1)),
    sample_b = str_remove(basename(file2), "\\.fa$"),
    ref_b    = basename(dirname(file2)),
  ) %>%
  filter(file1 != file2)  # remove self-comparisons

# Annotate with tissue from metadata
pairs <- pairs %>%
  left_join(meta %>% select(sample_id, tissue) %>% distinct(),
            by = c("sample_a" = "sample_id")) %>%
  rename(tissue_a = tissue) %>%
  left_join(meta %>% select(sample_id, tissue) %>% distinct(),
            by = c("sample_b" = "sample_id")) %>%
  rename(tissue_b = tissue)

# Classify each pair
pairs <- pairs %>%
  mutate(pair_type = case_when(
    tissue_a == tissue_b & ref_a != ref_b ~ "Same tissue, cross-reference",
    tissue_a == tissue_b & ref_a == ref_b ~ "Same tissue, same reference",
    tissue_a != tissue_b                  ~ "Different tissue",
    TRUE                                  ~ "Other"
  ))

# ── Statistical test ──────────────────────────────────────────────────────
cross_ref   <- pairs %>% filter(pair_type == "Same tissue, cross-reference") %>% pull(similarity)
diff_tissue <- pairs %>% filter(pair_type == "Different tissue")             %>% pull(similarity)

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

# ── Plot 1 of 2: comparison-type boxplot ──────────────────────────────────
pair_levels <- c(
  "Same tissue, same reference",
  "Same tissue, cross-reference",
  "Different tissue"
)

plot_data <- pairs %>%
  filter(pair_type %in% pair_levels) %>%
  mutate(pair_type = factor(pair_type, levels = pair_levels))

p_box <- ggplot(plot_data, aes(x = pair_type, y = similarity, fill = pair_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = c("#2166ac", "#f4a582", "#d1e5f0")) +
  labs(
    subtitle = sprintf("Wilcoxon p = %.2e", test_result$p.value),
    x        = "Comparison type",
    y        = "Minimizer-sketch Jaccard"
  ) +
  theme_cowplot(11) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 20, hjust = 1, size = 10),
    axis.title      = element_text(size = 11)
  )

# ── Plot 2 of 2: 9x9 sample-by-sample similarity heatmap ──────────────────
# Order rows/cols by tissue first, then by reference, so tissue-blocks are visible.
ref_order <- c("GRCh37", "GRCh38", "CHM13")     # adjust if metadata adds more
tissue_order <- meta %>%
  arrange(tissue) %>% pull(tissue) %>% unique()

label_df <- meta %>%
  select(sample_id, tissue) %>% distinct() %>%
  tidyr::expand_grid(ref = ref_order) %>%
  mutate(
    key   = paste(sample_id, ref, sep = "__"),
    label = paste(tissue, ref, sep = " · "),
    tissue = factor(tissue, levels = tissue_order)
  ) %>%
  arrange(tissue, factor(ref, levels = ref_order))

ordered_keys   <- label_df$key
ordered_labels <- label_df$label

# Long → matrix-like, including the diagonal (filtered out earlier as self-self)
mat_long <- pairs %>%
  mutate(key_a = paste(sample_a, ref_a, sep = "__"),
         key_b = paste(sample_b, ref_b, sep = "__")) %>%
  select(key_a, key_b, similarity)

# Add a self-self row per key (similarity = 1.0)
mat_long <- bind_rows(
  mat_long,
  tibble(key_a = ordered_keys, key_b = ordered_keys, similarity = 1.0)
)

heatmap_df <- mat_long %>%
  filter(key_a %in% ordered_keys, key_b %in% ordered_keys) %>%
  mutate(
    key_a = factor(key_a, levels = ordered_keys, labels = ordered_labels),
    key_b = factor(key_b, levels = ordered_keys, labels = ordered_labels)
  )

p_heat <- ggplot(heatmap_df, aes(x = key_b, y = key_a, fill = similarity)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1)) +
  scale_y_discrete(limits = rev) +
  coord_fixed() +
  labs(
    subtitle = "Pairwise similarity (rows/cols ordered by tissue → reference)",
    x = "Sample × Reference",
    y = "Sample × Reference",
    fill = "Jaccard"
  ) +
  theme_cowplot(11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y  = element_text(size = 9),
    axis.title.x = element_text(size = 11, margin = margin(t = 8)),
    axis.title.y = element_text(size = 11, margin = margin(r = 8)),
    axis.line    = element_line(color = "grey40"),
    axis.ticks   = element_line(color = "grey40"),
    panel.grid   = element_blank()
  )

# ── Combined output ───────────────────────────────────────────────────────
combined <- cowplot::plot_grid(
  p_box, p_heat,
  ncol = 2, rel_widths = c(1, 1.4),
  labels = c("A", "B")
)
title <- cowplot::ggdraw() +
  cowplot::draw_label("Experiment A: Cross-Reference Robustness",
                      fontface = "bold", x = 0.02, hjust = 0)
combined <- cowplot::plot_grid(title, combined, ncol = 1, rel_heights = c(0.06, 1))

ggsave(out_plot, combined, width = 14, height = 6, dpi = 300, device = cairo_png_in)
message("Done. Plot saved to ", out_plot)
