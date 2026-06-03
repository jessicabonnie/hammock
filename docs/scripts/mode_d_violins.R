#!/usr/bin/env Rscript
# Paper Fig 7 — sequence mode Pearson r and ARI distributions per k,
# across the full (w, p) sweep (no_ends column).
#
# Two violins per k: Pearson r vs bedtools, and ARI vs the 10 fetal-tissue
# labels. Pearson migrates upward with k; k=10 ARI is tall/multimodal while
# k>=15 collapses to a single ARI ~0.69 partition.
#
# Colours: Dark2 qualitative palette, matching the other paper figures
# (Fig 1 / Fig S5 / Fig 4).
#
# Data (copied into docs/data/ from
#   experiments/maurano_dhs_validation/results/mode_d_summary.csv):
#   mode_d_summary.csv
#
# Usage:
#   ml r/4.3.0 libjpeg/9c && Rscript docs/scripts/mode_d_violins.R

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

# Dark2 — same family as Fig 1 / Fig S5 / Fig 4.
metric_pal <- c("Pearson r vs bedtools" = "#1b9e77",
                "ARI vs tissue labels"  = "#d95f02")

theme_pub <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, colour = "gray30"),
        legend.position = "top",
        legend.title  = element_text(face = "bold"))

raw <- read_csv(file.path(data_dir, "mode_d_summary.csv"),
                show_col_types = FALSE)

d_long <- raw |>
  filter(column == "jaccard_similarity", reference == "bedtools") |>
  select(k, w, precision, pearson, ari) |>
  pivot_longer(c(pearson, ari), names_to = "metric", values_to = "value") |>
  mutate(k = factor(k),
         metric = recode(metric,
                         pearson = "Pearson r vs bedtools",
                         ari     = "ARI vs tissue labels"),
         metric = factor(metric, levels = names(metric_pal)))

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
  scale_fill_manual(values = metric_pal, name = NULL) +
  scale_colour_manual(values = metric_pal, name = NULL) +
  labs(x = "k-mer size k", y = "metric value",
       title = "Sequence mode metric distributions per k, across the full (w, p) sweep",
       subtitle = paste0(
         "no_ends only; Pearson migrates upward with k.\n",
         "k=10 ARI is tall/multimodal (sub-plateaus at ~0.09, 0.55, 0.69, 0.91); ",
         "k>=15 collapses to ARI~0.69 for every (w, p).")) +
  theme_pub

CairoPNG(file.path(figures_dir, "mode_d_violins_by_k.png"),
         width = 1300, height = 700, res = 150)
print(p_violin); dev.off()
cat("wrote mode_d_violins_by_k.png\n")
