#!/usr/bin/env Rscript
# make_metric_plots.R
# Generates three Mode D Pearson-vs-ARI summary figures from
# results/mode_d_summary.csv. Spec: ../PLOT_GENERATION.md
#
# Outputs (to ../figures/):
#   - mode_d_lines_p24.png        (Option B; replaces Fig 4)
#   - mode_d_max_lines_p24.png    (Option C; candidate)
#   - mode_d_violins_by_k.png     (replaces Fig 7)
#
# Run from repo root:
#   ml r/4.3.0
#   Rscript experiments/maurano_dhs_validation/scripts/make_metric_plots.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(cowplot)
  library(Cairo)
})

here <- "experiments/maurano_dhs_validation"
csv_path <- file.path(here, "results", "mode_d_summary.csv")
fig_dir  <- file.path(here, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

raw <- read_csv(csv_path, show_col_types = FALSE)

# ------------------------------------------------------------------
# Plot 1 — Option B: two-panel lines at p = 24
# ------------------------------------------------------------------
d24 <- raw |>
  filter(column == "jaccard_similarity",
         reference == "bedtools",
         precision == 24) |>
  mutate(k = factor(k))

p_left <- d24 |>
  ggplot(aes(w, pearson, colour = k, group = k)) +
  geom_line() + geom_point() +
  scale_x_log10() + ylim(0, 1) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = "window size w", y = "Pearson r vs bedtools",
       title = "Pearson r vs bedtools")

p_right <- d24 |>
  ggplot(aes(w, ari, colour = k, group = k)) +
  geom_line() + geom_point() +
  scale_x_log10() + ylim(0, 1) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = "window size w", y = "ARI vs tissue labels",
       title = "ARI vs tissue labels")

# cowplot doesn't have plot_annotation; put title/subtitle on a top strip
# and share the k legend by stripping per-panel legends then re-attaching.
shared_legend <- cowplot::get_legend(p_left + theme(legend.position = "right"))
p_left_nl  <- p_left  + theme(legend.position = "none")
p_right_nl <- p_right + theme(legend.position = "none")
panels <- cowplot::plot_grid(p_left_nl, p_right_nl, ncol = 2)
panels_with_legend <- cowplot::plot_grid(panels, shared_legend,
                                         ncol = 2, rel_widths = c(1, 0.08))
title <- cowplot::ggdraw() +
  cowplot::draw_label("Mode D — Pearson r and ARI vs window size, at p = 24",
                       fontface = "bold", size = 14, hjust = 0.5, x = 0.5, y = 0.7) +
  cowplot::draw_label("no_ends column only; lines connect points within the same k",
                       size = 11, hjust = 0.5, x = 0.5, y = 0.25)
p_lines <- cowplot::plot_grid(title, panels_with_legend,
                              ncol = 1, rel_heights = c(0.12, 1))

CairoPNG(file.path(fig_dir, "mode_d_lines_p24.png"),
         width = 1400, height = 600, res = 150)
print(p_lines); dev.off()

# ------------------------------------------------------------------
# Plot 2 — Option C: best-per-k at p = 24
# ------------------------------------------------------------------
best <- d24 |>
  mutate(k = as.integer(as.character(k))) |>
  group_by(k) |>
  summarise(pearson_max = max(pearson, na.rm = TRUE),
            pearson_w   = w[which.max(pearson)],
            ari_max     = max(ari, na.rm = TRUE),
            ari_w       = w[which.max(ari)],
            .groups = "drop") |>
  pivot_longer(c(pearson_max, ari_max),
               names_to = "metric", values_to = "value") |>
  mutate(best_w = ifelse(metric == "pearson_max", pearson_w, ari_w))

p_max <- best |>
  ggplot(aes(k, value, colour = metric, group = metric)) +
  geom_line() + geom_point(size = 3) +
  geom_text(aes(label = paste0("w=", best_w)),
            vjust = -1, size = 3, show.legend = FALSE) +
  ylim(0, 1) +
  scale_x_continuous(breaks = c(8, 10, 15, 20, 25)) +
  scale_colour_manual(values = c(pearson_max = "#1f77b4",
                                  ari_max     = "#d62728"),
                      labels = c(pearson_max = "max Pearson r",
                                  ari_max     = "max ARI")) +
  labs(x = "k-mer size k", y = "metric value", colour = NULL,
       title = "Mode D — best Pearson and best ARI per k, at p = 24",
       subtitle = "no_ends only; each point is the max over w at that (k, p); annotation shows the winning w")

CairoPNG(file.path(fig_dir, "mode_d_max_lines_p24.png"),
         width = 1100, height = 700, res = 150)
print(p_max); dev.off()

# ------------------------------------------------------------------
# Plot 3 — Violins per k across full (w, p) sweep
# ------------------------------------------------------------------
d_long <- raw |>
  filter(column == "jaccard_similarity", reference == "bedtools") |>
  select(k, w, precision, pearson, ari) |>
  pivot_longer(c(pearson, ari), names_to = "metric", values_to = "value") |>
  mutate(k = factor(k),
         metric = recode(metric,
                         pearson = "Pearson r vs bedtools",
                         ari     = "ARI vs tissue labels"))

p_violin <- d_long |>
  ggplot(aes(k, value, fill = metric)) +
  geom_violin(position = position_dodge(width = 0.85),
              alpha = 0.6, scale = "width") +
  geom_boxplot(position = position_dodge(width = 0.85),
               width = 0.12, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(aes(colour = metric),
              position = position_jitterdodge(jitter.width = 0.12,
                                              dodge.width = 0.85),
              size = 0.6, alpha = 0.4, show.legend = FALSE) +
  ylim(0, 1) +
  scale_fill_manual(values = c("Pearson r vs bedtools" = "#1f77b4",
                                "ARI vs tissue labels" = "#d62728")) +
  scale_colour_manual(values = c("Pearson r vs bedtools" = "#1f77b4",
                                  "ARI vs tissue labels" = "#d62728")) +
  labs(x = "k-mer size k", y = "metric value", fill = NULL,
       title = "Mode D metric distributions per k, across the full (w, p) sweep",
       subtitle = paste0(
         "no_ends only; Pearson migrates upward with k.\n",
         "k=10 ARI is tall/multimodal (sub-plateaus at ~0.09, 0.55, 0.69, 0.91); ",
         "k>=15 collapses to ARI~0.69 for every (w, p)."))

CairoPNG(file.path(fig_dir, "mode_d_violins_by_k.png"),
         width = 1300, height = 700, res = 150)
print(p_violin); dev.off()

# ------------------------------------------------------------------
# Sanity log
# ------------------------------------------------------------------
cat("=== sanity checks ===\n")
cat(sprintf("max Pearson at p=24:  %.4f  (expect ~0.9996)\n",
            max(d24$pearson, na.rm = TRUE)))
cat(sprintf("max ARI     at p=24:  %.4f  (expect ~0.9102)\n",
            max(d24$ari, na.rm = TRUE)))
cat("wrote:\n")
cat("  ", file.path(fig_dir, "mode_d_lines_p24.png"), "\n")
cat("  ", file.path(fig_dir, "mode_d_max_lines_p24.png"), "\n")
cat("  ", file.path(fig_dir, "mode_d_violins_by_k.png"), "\n")
