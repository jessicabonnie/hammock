#!/usr/bin/env Rscript
# Paper Fig 6 — sequence mode: numerical agreement (Pearson r vs bedtools)
# against biological-clustering recovery (ARI vs tissue labels), per config.
#
# The headline point: the best-numerical config and the best-clustering config
# are NOT the same cell. One point per sequence-mode (k, w) config at HLL
# precision p=24; circled points are the best of Pearson r, ARI, and MAE.
#
# Aesthetics match Fig 3 (mode_d_bedtools_vs_modeB_scatter): colour = k-mer
# size, point size = minimizer window w (discrete so small windows stay
# legible). k=5 is dropped (lone w=20 point at p=24, not a real sweep).
#
# (Sequence mode = orig "Mode D".)
#
# Data (copied into docs/data/ from
#   experiments/maurano_dhs_validation/results/mode_d_summary.csv):
#   mode_d_summary.csv
#
# Usage:
#   ml r/4.3.0 libjpeg/9c && Rscript docs/scripts/mode_d_metric_tradeoff.R

suppressPackageStartupMessages({
  library(dplyr)
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

tr <- raw |>
  filter(column == "jaccard_similarity", reference == "bedtools",
         precision == 24, !is.na(pearson), !is.na(ari),
         k >= 8) |>
  mutate(w = factor(w, levels = sort(unique(w))),
         label = sprintf("k=%d w=%s", k, as.character(w)))

best_pts <- bind_rows(
  tr |> arrange(desc(pearson)) |> slice(1),
  tr |> arrange(desc(ari))     |> slice(1),
  tr |> arrange(mae)           |> slice(1))

p_tr <- ggplot(tr, aes(x = pearson, y = ari, colour = factor(k), size = w)) +
  geom_point(alpha = 0.5) +
  geom_point(data = best_pts, aes(x = pearson, y = ari),
             size = 4, shape = 21, fill = NA, colour = "black",
             stroke = 1.2, inherit.aes = FALSE) +
  # best-Pearson and best-MAE sit on top of each other near (1, 0.69), at the
  # right panel edge; separate their labels vertically and right-align them
  # (nudged left) so neither clips. Order: Pearson, ARI, MAE.
  geom_text(data = best_pts, aes(x = pearson, y = ari, label = label),
            nudge_x = c(-0.015, 0, -0.015), nudge_y = c(0.05, 0.035, -0.05),
            hjust = c(1, 0.5, 1), size = 3, colour = "black",
            inherit.aes = FALSE) +
  scale_size_manual(
    values = setNames(seq(1.5, 7, length.out = nlevels(tr$w)), levels(tr$w)),
    name = "minimizer window w") +
  labs(x = "Pearson r vs bedtools (numerical agreement)",
       y = "ARI vs tissue labels (clustering recovery)",
       colour = "k-mer size",
       title = "Sequence Mode: numerical agreement vs biological-clustering recovery",
       caption = paste(
         "Each point is one sequence-mode (k, w) config at HLL precision p=24.",
         "Circled points are the best config for Pearson r, ARI, and MAE —",
         "the best-numerical and best-clustering configs are different cells.",
         sep = "\n")) +
  theme_minimal(base_size = 11) +
  theme(plot.caption = element_text(hjust = 0, size = 9, colour = "gray30"))

CairoPNG(file.path(figures_dir, "mode_d_metric_tradeoff.png"),
         width = 1200, height = 900, res = 150)
print(p_tr); dev.off()
cat("wrote mode_d_metric_tradeoff.png\n")
