#!/usr/bin/env Rscript
# =============================================================================
# scripts/exp_a_dendrogram.R
#
# Hierarchical clustering of the 9 Exp A sample x reference combinations.
# If cross-reference robustness holds, the 3 references of each tissue should
# cluster together (tissue identity dominates reference choice).
#
# Usage:
#   Rscript scripts/exp_a_dendrogram.R \
#     <broad_csv> <narrow_csv> <metadata_tsv> <out_png> [kw_label]
#
# CSV is hammock Mode D long-form (file1,file2,...,jaccard_similarity_with_ends).
# kw_label appears in the figure title (e.g. "k=15, w=15"); defaults to "k=10, w=10".
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(tibble); library(ggplot2); library(ggdendro); library(cowplot)
  library(Cairo)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4)
broad_csv <- args[1]; narrow_csv <- args[2]
meta_file <- args[3]; out_png   <- args[4]
kw_label  <- if (length(args) >= 5) args[5] else "k=10, w=10"

meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  mutate(key = paste(sample_id, ref, sep = "__"))

# Tissue palette — fixed across panels for readable comparison.
tissue_pal <- c(heart = "#c1272d", liver = "#7f3f00", lung = "#4575b4")

build_dendro <- function(csv_file, peak_type) {
  mat <- read_csv(csv_file, show_col_types = FALSE) %>%
    mutate(sample_a = str_remove(basename(file1), "\\.fa$"),
           ref_a    = basename(dirname(file1)),
           sample_b = str_remove(basename(file2), "\\.fa$"),
           ref_b    = basename(dirname(file2)),
           key_a    = paste(sample_a, ref_a, sep = "__"),
           key_b    = paste(sample_b, ref_b, sep = "__"))

  keys <- sort(unique(c(mat$key_a, mat$key_b)))
  n    <- length(keys)
  sim  <- matrix(NA_real_, n, n, dimnames = list(keys, keys))
  for (i in seq_len(nrow(mat))) {
    sim[mat$key_a[i], mat$key_b[i]] <- mat$jaccard_similarity_with_ends[i]
  }
  # Symmetrize defensively, then derive distance.
  sim[is.na(sim)] <- t(sim)[is.na(sim)]
  diag(sim) <- 1.0
  dist_mat <- as.dist(1 - sim)

  hc   <- hclust(dist_mat, method = "average")
  dd   <- dendro_data(as.dendrogram(hc), type = "rectangle")

  # Render leaves on the right, root on the left, without coord_flip:
  # plot-x = original height (distance), plot-y = original x (leaf position).
  s    <- segment(dd)
  segs <- data.frame(x    = s$y,    y    = s$x,
                     xend = s$yend, yend = s$xend)

  lbls <- label(dd) %>%
    left_join(meta %>% select(key, tissue, ref), by = c("label" = "key")) %>%
    mutate(pretty = paste(tissue, ref, sep = " · "),
           plot_x = 0, plot_y = x)

  x_max <- max(segs$x, segs$xend)
  pad_x <- x_max * 0.45  # room for labels on the right

  ggplot() +
    geom_segment(data = segs,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.5) +
    geom_label(data = lbls,
               aes(x = plot_x, y = plot_y, label = pretty, color = tissue),
               hjust = 0, nudge_x = -x_max * 0.03, size = 3.4,
               fill = "white", label.size = 0, label.padding = unit(0.1, "lines")) +
    scale_color_manual(values = tissue_pal, guide = "none") +
    scale_x_reverse(limits = c(x_max * 1.05, -pad_x),
                    breaks = scales::pretty_breaks(4)(c(0, x_max)),
                    expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = 0.08)) +
    labs(title = sprintf("%s peaks", peak_type),
         x = "1 − Jaccard (minimizer + ends)", y = NULL) +
    theme_cowplot(11) +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank(),
          plot.title   = element_text(face = "bold", size = 12),
          plot.margin  = margin(8, 12, 8, 8))
}

p_broad  <- build_dendro(broad_csv,  "Broad")
p_narrow <- build_dendro(narrow_csv, "Narrow")

title <- ggdraw() +
  draw_label(
    sprintf("Experiment A: tissue clusters dominate reference choice (%s)", kw_label),
    fontface = "bold", x = 0.02, hjust = 0, size = 13)

subtitle <- ggdraw() +
  draw_label(
    paste("Hierarchical clustering (UPGMA on 1 − Jaccard).",
          "Tissue color = heart (red), liver (brown), lung (blue)."),
    x = 0.02, hjust = 0, size = 10, color = "grey30")

body <- plot_grid(p_broad, p_narrow, ncol = 2)
out  <- plot_grid(title, subtitle, body, ncol = 1,
                  rel_heights = c(0.06, 0.04, 1))

# 11" x 5" at 300 dpi — fits next to existing figures.
Cairo::CairoPNG(out_png, width = 11, height = 5, units = "in", res = 300)
print(out)
invisible(dev.off())
message("Wrote ", out_png)
