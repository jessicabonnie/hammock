#!/usr/bin/env Rscript
# =============================================================================
# scripts/exp_a_sweep_summary.R
# Experiment A — across-(k,w) sweep summary.
#
# Reads every per-(k,w) cross_ref_stats.tsv for ONE peak_type and produces
# a (k × w) heatmap of effect size = median(same-tissue cross-ref) - median(different-tissue).
# Bigger = better separation between within-tissue cross-ref similarity and
# different-tissue similarity, i.e., more reference-robust at that (k,w).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(Cairo)
})

cairo_png_in <- function(filename, width, height, ...) {
  Cairo::CairoPNG(filename, width = width, height = height,
                  units = "in", res = 300, ...)
}

# ── Snakemake I/O ─────────────────────────────────────────────────────────
stats_files <- snakemake@input[["stats"]]
out_plot    <- snakemake@output[["plot"]]
peak_type   <- snakemake@wildcards[["peak_type"]]

# ── Read stats from each (k, w) ───────────────────────────────────────────
# Path pattern: .../exp_a/<peak_type>/k<k>_w<w>/cross_ref_stats.tsv
all_stats <- tibble(file = stats_files) %>%
  mutate(
    k = as.integer(str_extract(file, "(?<=/k)\\d+(?=_w)")),
    w = as.integer(str_extract(file, "(?<=_w)\\d+(?=/)"))
  ) %>%
  rowwise() %>%
  mutate(stats = list(read_tsv(file, show_col_types = FALSE))) %>%
  ungroup() %>%
  unnest(stats) %>%
  mutate(
    effect = median_cross_ref - median_diff_tissue,
    log10p = -log10(pmax(wilcoxon_p, 1e-300))
  )

# ── Heatmap ───────────────────────────────────────────────────────────────
k_levels <- sort(unique(all_stats$k))
w_levels <- sort(unique(all_stats$w))

p <- ggplot(all_stats, aes(x = factor(k, levels = k_levels),
                           y = factor(w, levels = w_levels),
                           fill = effect)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f\np=%.0e", effect, wilcoxon_p)),
            size = 2.6, lineheight = 0.85) +
  scale_fill_gradient2(
    low      = "#2166ac",
    mid      = "white",
    high     = "#b2182b",
    midpoint = 0,
    name     = "median(cross-ref)\n − median(diff-tissue)"
  ) +
  coord_fixed() +
  labs(
    title    = sprintf("Experiment A — %s peaks: parameter sweep", peak_type),
    subtitle = "Each cell: effect size and Wilcoxon p for that (k, w). Higher (red) = better cross-ref robustness.",
    x        = "k-mer size",
    y        = "window size"
  ) +
  theme_cowplot(11) +
  theme(panel.grid = element_blank())

ggsave(out_plot, p, width = 9, height = 7, dpi = 300, device = cairo_png_in)
message("Done. Sweep summary saved to ", out_plot)
