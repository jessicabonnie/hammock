#!/usr/bin/env Rscript
# analyze_part2_synthetic.R — Part 2 of modeD_flanking.
#
# Joins:
#   results/part2_sweep.csv           (Mode D both columns × pairs × k×w×p)
#   results/part2_synthetic_truth.tsv (exact canonical k-mer Jaccards per k)
#
# Computes error of each column vs the exact k-mer Jaccard at the matching k,
# the flanking fraction φ, and emits:
#   results/part2_summary.csv
#   figures/synthetic_delta_vs_phi.png
#   figures/synthetic_phase_diagram.png
#
# Usage:
#   ml gcc/9.3.0 r/4.3.0 libjpeg/9c
#   Rscript analyze_part2_synthetic.R

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(scales); library(Cairo)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) >= 1) {
  dirname(normalizePath(sub("--file=", "", file_arg[1])))
} else { getwd() }

results_dir <- file.path(script_dir, "results")
figures_dir <- file.path(script_dir, "figures")
sweep_path  <- file.path(results_dir, "part2_sweep.csv")
truth_path  <- file.path(results_dir, "part2_synthetic_truth.tsv")

stopifnot(file.exists(sweep_path), file.exists(truth_path))
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(path, plot, width = 8, height = 5, dpi = 150) {
  CairoPNG(filename = path, width = width, height = height,
           units = "in", res = dpi)
  on.exit(dev.off())
  print(plot)
}

# ── inputs ──────────────────────────────────────────────────────────────────
sweep <- read_csv(sweep_path, show_col_types = FALSE)
truth <- read_tsv(truth_path, show_col_types = FALSE) %>%
  select(label, k, j_truth = jaccard,
         n_kmers_A, n_kmers_B, n_intersection, n_union)

# Manifest carries realised total_len per file (lognormal makes nominal
# mean_len × n_intervals a poor proxy at small n_intervals or large σ).
manifest_path <- file.path(script_dir, "data", "synthetic", "manifest.tsv")
stopifnot(file.exists(manifest_path))
manifest <- read_tsv(manifest_path, show_col_types = FALSE) %>%
  select(label, total_len_A, total_len_B)

# ── join + per-pair flanking budget ─────────────────────────────────────────
# φ_pair = sqrt(φ_A · φ_B) where φ_X = 2(k-1) · n_intervals / (total_len_X / w)
pairs <- sweep %>%
  inner_join(truth,    by = c("label", "k")) %>%
  inner_join(manifest, by = "label") %>%
  mutate(
    phi_A         = (2 * (k - 1) * n_intervals) / pmax(1, total_len_A / w),
    phi_B         = (2 * (k - 1) * n_intervals) / pmax(1, total_len_B / w),
    phi           = sqrt(phi_A * phi_B),
    err_no_ends   = j_no_ends   - j_truth,
    err_with_ends = j_with_ends - j_truth,
    abs_no        = abs(err_no_ends),
    abs_with      = abs(err_with_ends),
    delta_abs     = abs_with - abs_no,           # >0 ⇒ no_ends better
    winner        = if_else(delta_abs > 0, "no_ends", "with_ends")
  )

write_csv(pairs, file.path(results_dir, "part2_pairs.csv"))
cat("Wrote part2_pairs.csv (", nrow(pairs), "rows)\n")

# ── aggregated summary per (k, w, p, mutation, dist, n_intervals, mean_len) ─
agg <- pairs %>%
  group_by(k, w, p, mutation, dist, n_intervals, mean_len) %>%
  summarise(
    n            = n(),
    j_truth      = mean(j_truth),
    j_no_ends    = mean(j_no_ends),
    j_with_ends  = mean(j_with_ends),
    mae_no_ends  = mean(abs_no),
    mae_with_ends= mean(abs_with),
    delta_mae    = mae_with_ends - mae_no_ends,
    mean_phi     = mean(phi),
    .groups = "drop"
  )
write_csv(agg, file.path(results_dir, "part2_summary.csv"))
cat("Wrote part2_summary.csv (", nrow(agg), "rows)\n")

# ── plot 1: delta_abs vs φ, coloured by mutation, faceted by k ──────────────
p1 <- ggplot(pairs,
             aes(x = phi, y = delta_abs,
                 colour = factor(mutation))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.5, size = 1.3) +
  facet_grid(k ~ p, labeller = label_both) +
  scale_x_log10() +
  labs(x = expression("flanking fraction "*phi*" (log)"),
       y = "|err with_ends| - |err no_ends|",
       title = "Synthetic Mode D: sign-of-flanking-gain vs φ",
       subtitle = "y < 0 ⇒ with_ends is better; y > 0 ⇒ no_ends is better",
       colour = "mutation") +
  theme_minimal(base_size = 10)
save_png(file.path(figures_dir, "synthetic_delta_vs_phi.png"), p1,
         width = 10, height = 7)

# ── plot 2: phase diagram — fraction-of-pairs-won-by-with_ends in (φ, μ) ───
phase <- pairs %>%
  mutate(phi_bin = cut(phi, breaks = quantile(phi, probs = seq(0, 1, 0.1),
                                              na.rm = TRUE),
                       include.lowest = TRUE)) %>%
  group_by(phi_bin, mutation) %>%
  summarise(frac_with_ends = mean(winner == "with_ends"),
            n_cells = n(),
            .groups = "drop")

p2 <- ggplot(phase,
             aes(x = phi_bin, y = factor(mutation),
                 fill = frac_with_ends)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", frac_with_ends)), size = 3) +
  scale_fill_gradient2(midpoint = 0.5, low = "#d7191c",
                       mid = "#ffffbf", high = "#2c7bb6", limits = c(0, 1)) +
  labs(x = "φ bin (deciles)", y = "mutation rate",
       fill = "frac. wins by\nwith_ends",
       title = "Phase diagram: fraction of (k,w,p) cells where with_ends beats no_ends") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_png(file.path(figures_dir, "synthetic_phase_diagram.png"), p2,
         width = 10, height = 5)

cat("Done.\n")
