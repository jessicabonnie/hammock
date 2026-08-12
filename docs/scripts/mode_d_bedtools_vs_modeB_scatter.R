#!/usr/bin/env Rscript
# Paper Fig 3 — sequence mode: agreement with bedtools vs agreement with
# interval mode, per (k, w, p) config.
#
# One point per config: x = Pearson r of sequence mode vs bedtools,
# y = Pearson r of sequence mode vs interval mode. Points on the y=x
# diagonal mean interval mode and bedtools are interchangeable references —
# once sequence mode is validated against bedtools it is validated against
# any reasonable interval-Jaccard estimator.
#
# (Sequence mode = orig "Mode D"; interval mode = orig "Mode B". The data
# file still encodes the interval-mode reference as "mode_B".)
#
# Data (copied into docs/data/ from
#   experiments/maurano_dhs_validation/results/mode_d_summary.csv):
#   mode_d_summary.csv
#
# Usage:
#   ml r/4.3.0 libjpeg/9c && Rscript docs/scripts/mode_d_bedtools_vs_modeB_scatter.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(Cairo)
})

script_dir <- dirname(normalizePath(sub("--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
data_dir    <- file.path(script_dir, "..", "data")
figures_dir <- file.path(script_dir, "..", "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read_csv(file.path(data_dir, "mode_d_summary.csv"),
                show_col_types = FALSE)

# `column` holds a data value naming which metric a summary row is about, not
# a CSV header, so the reg_eq_similarity/jaccard_similarity fallback here is a
# value comparison against unique(raw$column), not a names(df) presence check.
# See docs/seed-jaccard-reg-eq-rename.md Step 2.
reg_eq_value <- if ("reg_eq_similarity" %in% unique(raw$column)) {
  "reg_eq_similarity"
} else {
  message(
    "reg_eq_similarity not found in 'column' values; falling back to jaccard_similarity"
  )
  "jaccard_similarity"
}

refs_wide <- raw |>
  # k=5 has only a single surviving (w=20) point at p=24 — not a real sweep,
  # so drop it rather than imply coverage it doesn't have.
  filter(column == reg_eq_value, !is.na(pearson), precision == 24,
         k >= 8) |>
  select(k, w, precision, reference, pearson) |>
  pivot_wider(names_from = reference, values_from = pearson) |>
  filter(!is.na(bedtools), !is.na(mode_B)) |>
  mutate(w = factor(w, levels = sort(unique(w))))

p_refs <- ggplot(refs_wide,
                 aes(x = bedtools, y = mode_B,
                     colour = factor(k), size = w)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed") +
  geom_point(alpha = 0.5) +
  scale_size_manual(
    values = setNames(seq(1.5, 7, length.out = nlevels(refs_wide$w)),
                      levels(refs_wide$w)),
    name = "minimizer window w") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Pearson r vs bedtools", y = "Pearson r vs interval mode (p=21)",
       colour = "k-mer size",
       title = "Agreement between Sequence and Interval modes",
       caption = paste(
         "Each point is one sequence-mode (k, w) config at HLL precision p=24:",
         "its Pearson r vs bedtools (x) against its Pearson r vs interval mode (y).",
         "Points on y=x ⇒ interval mode and bedtools are interchangeable as references.",
         sep = "\n")) +
  theme_minimal(base_size = 11) +
  theme(plot.caption = element_text(hjust = 0, size = 9, colour = "gray30"))

CairoPNG(file.path(figures_dir, "mode_d_bedtools_vs_modeB_scatter.png"),
         width = 1100, height = 900, res = 150)
print(p_refs); dev.off()
cat("wrote mode_d_bedtools_vs_modeB_scatter.png\n")
