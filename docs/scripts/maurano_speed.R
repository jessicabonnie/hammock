#!/usr/bin/env Rscript
# Paper figures for the Maurano speed story (subB_mixed_stride experiment).
#
# Produces two figures:
#   1. maurano_speed_bars.png  (MAIN, Fig 1) — grouped wall-time bars:
#      bedtools vs hammock at no-subsample / subB=0.1 / subB=0.01
#      (mixed-stride), each hammock bar annotated with speedup over
#      bedtools and the per-pair Jaccard change vs hammock's OWN
#      no-subsample output (not vs bedtools).
#   2. maurano_subB_pareto_scatter.png  (SUPPLEMENT) — de-zigzagged
#      replacement for the old headline pareto: all three subB methods,
#      points only (no subB-order connecting line), so the doubling-back
#      is gone. Shows that only mixed-stride converts subsampling into
#      speed.
#
# Data (copied into docs/data/ from experiments/subB_mixed_stride/results/):
#   maurano_subB_summary.csv  — per (method, subB) wall_median + mae
#   maurano_bedtools.csv      — bedtools per-pair runs (for reference wall)
#
# Usage:
#   ml r/4.3.0 && Rscript docs/scripts/maurano_speed.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(Cairo)
})

script_dir <- dirname(normalizePath(sub("--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
data_dir    <- file.path(script_dir, "..", "data")
figures_dir <- file.path(script_dir, "..", "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

wrapcap <- function(s, width = 110) paste(strwrap(s, width = width), collapse = "\n")

save_png <- function(path, plot, width, height, dpi = 200) {
  CairoPNG(filename = path, width = width, height = height,
           units = "in", res = dpi)
  on.exit(dev.off())
  print(plot)
}

method_levels <- c("hash-threshold", "mixed-stride", "single-hash")
method_pal <- c(
  "hash-threshold" = "#1b9e77",
  "mixed-stride"   = "#d95f02",
  "single-hash"    = "#7570b3"
)
bedtools_col <- "#333333"

theme_pub <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 15),
        plot.subtitle    = element_text(size = 11, color = "gray30"),
        plot.caption     = element_text(size = 9,  color = "gray40", hjust = 0),
        legend.position  = "top",
        legend.title     = element_text(size = 11, face = "bold"),
        legend.text      = element_text(size = 10))

# ---------- load ----------
summ <- read_csv(file.path(data_dir, "maurano_subB_summary.csv"),
                 show_col_types = FALSE)
summ$method <- factor(summ$method, levels = method_levels)

bt <- read_csv(file.path(data_dir, "maurano_bedtools.csv"),
               show_col_types = FALSE)
bt_wall <- bt %>%
  distinct(rep, run_id, wall_time) %>%
  summarise(median = median(wall_time)) %>% pull(median)
cat(sprintf("bedtools median wall: %.3fs\n", bt_wall))

summ <- summ %>% mutate(speedup_bt = bt_wall / wall_median)

# ========== 1. MAIN — grouped wall-time bars (mixed-stride path) ==========
ms <- summ %>% filter(method == "mixed-stride")
ms_nosub <- ms %>% filter(subB == 1.0) %>% pull(wall_median)

bars <- tibble(
  label = c("bedtools", "hammock\n(no subsample)",
            "hammock\nsubB = 0.1", "hammock\nsubB = 0.01"),
  tool  = c("bedtools", "hammock", "hammock", "hammock"),
  wall  = c(bt_wall,
            ms %>% filter(subB == 1.0)  %>% pull(wall_median),
            ms %>% filter(subB == 0.1)  %>% pull(wall_median),
            ms %>% filter(subB == 0.01) %>% pull(wall_median)),
  mae   = c(NA,
            0,
            ms %>% filter(subB == 0.1)  %>% pull(mae),
            ms %>% filter(subB == 0.01) %>% pull(mae))
)
bars <- bars %>%
  mutate(order = row_number(),
         label = factor(label, levels = label[order(order)]),
         speedup = bt_wall / wall,
         # annotation above each bar
         tag = ifelse(tool == "bedtools",
                      sprintf("%.1f s\n1.00× (ref)", wall),
                      sprintf("%.1f s\n%.2f× faster\nΔJ = %s vs no-sub",
                              wall, speedup,
                              ifelse(mae == 0, "0",
                                     formatC(mae, format = "e", digits = 0)))))

p_bars <- ggplot(bars, aes(x = label, y = wall, fill = tool)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = tag), vjust = -0.25, size = 3.5,
            lineheight = 0.95, fontface = "plain", color = "gray15") +
  scale_fill_manual(values = c("bedtools" = bedtools_col,
                               "hammock"  = method_pal[["mixed-stride"]]),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 14),
                     breaks = seq(0, 14, 2),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "hammock is faster than bedtools at every subsampling level",
       subtitle = sprintf(
         "Maurano fetal-tissue DHS: 20 BEDs, 190 pairs, interval mode, p=18, 8 threads. Bedtools wall = %.2f s (8-way GNU parallel).",
         bt_wall),
       x = NULL, y = "wall time (s)  —  lower is faster",
       caption = wrapcap("Even with no subsampling the HLL path beats bedtools (1.16×); mixed-stride subB extends the lead to ~3×. ΔJ is the mean per-pair Jaccard change vs hammock's own no-subsample output (the speed knob is near-free), not the gap to bedtools.")) +
  theme_pub +
  theme(axis.text.x = element_text(size = 11))
save_png(file.path(figures_dir, "maurano_speed_bars.png"),
         p_bars, width = 9.5, height = 6)

# ========== 2. SUPPLEMENT — de-zigzagged subB-method scatter ==========
# Same axes as the retired pareto, but points only (no subB-order path),
# so the structural doubling-back is gone. subB=1.0 dropped (MAE=0 -> -Inf
# on log x). Label subB on the mixed-stride points only to keep it legible.
sc <- summ %>% filter(subB < 1.0)
ms_lab <- sc %>% filter(method == "mixed-stride")

p_sc <- ggplot(sc, aes(x = mae, y = speedup_bt, color = method)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = bedtools_col,
             linewidth = 0.35) +
  annotate("text", x = max(sc$mae), y = 1.02, label = "bedtools  →  1.00×",
           hjust = 1, vjust = 0, color = bedtools_col, size = 3.6,
           fontface = "italic") +
  geom_point(size = 3.8) +
  geom_text(data = ms_lab, aes(label = paste0("subB=", subB)),
            color = "gray20", size = 3.0, fontface = "bold",
            nudge_y = 0.12, show.legend = FALSE) +
  scale_color_manual(values = method_pal, name = NULL) +
  scale_x_log10(labels = label_number(accuracy = 0.0001)) +
  scale_y_continuous(breaks = seq(0.5, 3.5, 0.5), limits = c(0.55, 3.4)) +
  labs(title = "Only mixed-stride converts subsampling into speed",
       subtitle = sprintf(
         "Maurano DHS (20 BEDs, 190 pairs, interval mode, p=18, 8 threads). Bedtools wall = %.2f s.\nPoints = one --subB level per method (0.5, 0.25, 0.1, 0.05, 0.01); top-left is best.",
         bt_wall),
       x = "Mean absolute per-pair Jaccard error vs no-subsample (log scale)",
       y = "Speedup vs bedtools",
       caption = wrapcap("Replaces the earlier connected-line version, whose lines doubled back because they joined points in subB order (not an axis). Mixed-stride reaches ~2× at subB=0.1 and ~3× at subB=0.01; hash-threshold and single-hash stay near 1× — their per-position gate hash masks the savings from fewer HLL ingests.")) +
  theme_pub
save_png(file.path(figures_dir, "maurano_subB_pareto_scatter.png"),
         p_sc, width = 10, height = 6.2)

cat("wrote maurano_speed_bars.png + maurano_subB_pareto_scatter.png\n")
