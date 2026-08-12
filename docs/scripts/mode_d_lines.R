#!/usr/bin/env Rscript
# Paper Fig 4 — sequence mode Pearson r and ARI vs window size, at p = 24.
#
# Two-panel line plot, no_ends column (jaccard_similarity), one line per k.
# Left: Pearson r vs bedtools. Right: ARI vs the 10 fetal-tissue labels.
# Carries the "Pearson is a ridge, ARI is a single peak" contrast.
#
# Colours: Dark2 qualitative palette, matching the other paper figures
# (Fig 1 / Fig S5 use the first three Dark2 colours).
#
# Data (copied into docs/data/ from
#   experiments/maurano_dhs_validation/results/mode_d_summary.csv):
#   mode_d_summary.csv
#
# Usage:
#   ml r/4.3.0 libjpeg/9c && Rscript docs/scripts/mode_d_lines.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(cowplot)
  library(Cairo)
})

script_dir <- dirname(normalizePath(sub("--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
data_dir    <- file.path(script_dir, "..", "data")
figures_dir <- file.path(script_dir, "..", "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Dark2, ordered to the five k levels — same family as Fig 1 / Fig S5.
k_levels <- c("8", "10", "15", "20", "25")
k_pal <- setNames(
  c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"),
  k_levels)

theme_pub <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title   = element_text(face = "bold", size = 13),
        legend.title = element_text(face = "bold"))

raw <- read_csv(file.path(data_dir, "mode_d_summary.csv"),
                show_col_types = FALSE)

d24 <- raw |>
  filter(column == "jaccard_similarity",
         reference == "bedtools",
         precision == 24,
         k %in% as.integer(k_levels)) |>   # drop k=5 (lone w=20 point, not a sweep)
  mutate(k = factor(as.character(k), levels = k_levels))

panel <- function(yvar, ylab, title) {
  ggplot(d24, aes(w, .data[[yvar]], colour = k, group = k)) +
    geom_line(linewidth = 0.6) + geom_point(size = 2) +
    scale_x_log10() + ylim(0, 1) +
    scale_colour_manual(values = k_pal, name = "k") +
    labs(x = "window size w", y = ylab, title = title) +
    theme_pub
}

p_left  <- panel("pearson", "Pearson r vs bedtools",   "Pearson r vs bedtools")
p_right <- panel("ari",     "ARI vs tissue labels",     "ARI vs tissue labels")

shared_legend <- cowplot::get_legend(p_left + theme(legend.position = "right"))
panels <- cowplot::plot_grid(p_left  + theme(legend.position = "none"),
                             p_right + theme(legend.position = "none"),
                             ncol = 2)
panels_with_legend <- cowplot::plot_grid(panels, shared_legend,
                                         ncol = 2, rel_widths = c(1, 0.08))
title <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Sequence mode — Pearson r and ARI vs window size, at p = 24",
    fontface = "bold", size = 14, hjust = 0.5, x = 0.5, y = 0.5)
p_lines <- cowplot::plot_grid(title, panels_with_legend,
                              ncol = 1, rel_heights = c(0.1, 1))

CairoPNG(file.path(figures_dir, "mode_d_lines_p24.png"),
         width = 1400, height = 600, res = 150)
print(p_lines); dev.off()
cat("wrote mode_d_lines_p24.png\n")
