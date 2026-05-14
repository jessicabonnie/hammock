#!/usr/bin/env Rscript
# Publication-style headline figures for the subB-method experiment.
#
# Three plots:
#   1. headline_maurano_pareto.png  — speedup-vs-accuracy tradeoff on real
#      DHS data, with bedtools reference. The money plot.
#   2. headline_hammock_vs_bedtools.png — wall-time bars on Maurano showing
#      that even hammock at subB=1.0 beats bedtools.
#   3. headline_synthetic_scaling.png — peak mixed-stride speedup across
#      synthetic size classes; shows how the win compounds with corpus size.
#
# Run after analyze.R has been executed at least once.
#
# Usage:
#   ml r/4.3.0 && Rscript headline_figures.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(Cairo)
})

script_dir <- dirname(normalizePath(sub("--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
results_dir <- file.path(script_dir, "results")
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(path, plot, width, height, dpi = 200) {
  CairoPNG(filename = path, width = width, height = height,
           units = "in", res = dpi)
  on.exit(dev.off())
  print(plot)
}

pick_latest <- function(pattern) {
  files <- sort(list.files(results_dir, pattern = pattern,
                           full.names = TRUE), decreasing = TRUE)
  if (length(files) == 0) return(NULL)
  files[1]
}

method_levels <- c("hash-threshold", "mixed-stride", "single-hash")
method_pal <- c(
  "hash-threshold" = "#1b9e77",
  "mixed-stride"   = "#d95f02",
  "single-hash"    = "#7570b3"
)
bedtools_col <- "#333333"

# ---------- load data ----------
maur_csv <- pick_latest("^sweep_maurano_[0-9].*\\.csv$")
syn_csv  <- pick_latest("^sweep_synthetic_.*\\.csv$")
if (is.null(syn_csv)) syn_csv <- pick_latest("^sweep_[0-9].*\\.csv$")
bt_csv   <- pick_latest("^sweep_maurano_bedtools_.*\\.csv$")
if (is.null(maur_csv)) stop("missing maurano sweep CSV")
if (is.null(bt_csv))   stop("missing bedtools maurano CSV")
if (is.null(syn_csv))  stop("missing synthetic sweep CSV")
cat("maurano: ", maur_csv, "\nbedtools:", bt_csv, "\nsynthetic:", syn_csv, "\n", sep = "")

maur <- read_csv(maur_csv, show_col_types = FALSE)
syn  <- read_csv(syn_csv,  show_col_types = FALSE)
bt   <- read_csv(bt_csv,   show_col_types = FALSE)

bt_wall <- bt %>%
  distinct(rep, run_id, wall_time) %>%
  summarise(median = median(wall_time)) %>% pull(median)
cat(sprintf("bedtools median wall: %.3fs\n", bt_wall))

# ---------- summaries ----------
summarise_sweep <- function(df) {
  df$method <- factor(df$method, levels = method_levels)
  gt <- df %>% filter(subB == 1.0) %>%
    group_by(size_class, file_a, file_b) %>%
    summarise(j_truth = mean(jaccard), .groups = "drop")
  acc <- df %>%
    inner_join(gt, by = c("size_class", "file_a", "file_b")) %>%
    mutate(abs_err = abs(jaccard - j_truth)) %>%
    group_by(method, size_class, subB) %>%
    summarise(mae = mean(abs_err),
              max_err = max(abs_err),
              .groups = "drop")
  runs <- df %>% distinct(method, size_class, subB, rep, run_id, wall_time)
  speed <- runs %>%
    group_by(method, size_class, subB) %>%
    summarise(wall_median = median(wall_time), .groups = "drop")
  speed %>% left_join(acc, by = c("method", "size_class", "subB"))
}

maur_sum <- summarise_sweep(maur) %>%
  mutate(speedup_bt = bt_wall / wall_median)
syn_sum  <- summarise_sweep(syn)

# baseline (no-subsample wall, averaged across methods) per size_class
syn_baseline <- syn_sum %>% filter(subB == 1.0) %>%
  group_by(size_class) %>%
  summarise(wall_nosub = median(wall_median), .groups = "drop")
syn_sum <- syn_sum %>%
  left_join(syn_baseline, by = "size_class") %>%
  mutate(speedup_nosub = wall_nosub / wall_median)

theme_pub <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 15),
        plot.subtitle    = element_text(size = 11, color = "gray30"),
        plot.caption     = element_text(size = 9,  color = "gray40", hjust = 0),
        legend.position  = "top",
        legend.box       = "horizontal",
        legend.title     = element_text(size = 11, face = "bold"),
        legend.text      = element_text(size = 10))

# ========== 1. Maurano Pareto: speedup vs accuracy ==========
# Drop subB=1.0 (MAE=0, would compress the x-axis). Annotate the rec config.
pareto <- maur_sum %>% filter(subB < 1.0)
rec <- pareto %>% filter(method == "mixed-stride", subB == 0.1)

# Only label the mixed-stride curve so the cloud is readable
ms_labels <- pareto %>% filter(method == "mixed-stride")

p1 <- ggplot(pareto, aes(x = mae, y = speedup_bt, color = method)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = bedtools_col,
             linewidth = 0.35) +
  annotate("text", x = max(pareto$mae), y = 1.02,
           label = "bedtools  →  1.00×",
           hjust = 1, vjust = 0, color = bedtools_col, size = 3.6,
           fontface = "italic") +
  geom_path(aes(group = method), linewidth = 0.5, alpha = 0.9) +
  geom_point(size = 3.6) +
  geom_text(data = ms_labels, aes(label = paste0("subB=", subB)),
            color = "gray20", size = 3.0, fontface = "bold",
            nudge_y = 0.13, show.legend = FALSE) +
  geom_point(data = rec, aes(x = mae, y = speedup_bt),
             shape = 21, size = 8, stroke = 0.8,
             color = method_pal[["mixed-stride"]], fill = NA) +
  annotate("text", x = rec$mae * 1.06, y = rec$speedup_bt - 0.13,
           label = "recommended",
           hjust = 0, vjust = 1, fontface = "italic", size = 3.6,
           color = method_pal[["mixed-stride"]]) +
  scale_color_manual(values = method_pal, name = NULL) +
  scale_x_log10(labels = label_number(accuracy = 0.0001)) +
  scale_y_continuous(breaks = seq(0.5, 3.5, 0.5),
                     limits = c(0.55, 3.4)) +
  labs(title = "Hammock Mode B beats bedtools on real DHS data",
       subtitle = sprintf(
         "20 Maurano fetal-tissue DHS BEDs, 190 pairs, Mode B p=18, 8 threads. Bedtools wall = %.2fs.\nEach line sweeps subB ∈ {0.5, 0.25, 0.1, 0.05, 0.01}; top-left is better.",
         bt_wall),
       x = "Mean absolute Jaccard error vs subB=1.0 (log scale)",
       y = "Speedup vs bedtools",
       caption = "Mixed-stride: ~2× faster than bedtools at MAE < 0.001, up to 3× at MAE < 0.002. The other two methods barely move.") +
  theme_pub
save_png(file.path(figures_dir, "headline_maurano_pareto.png"),
         p1, width = 10, height = 6.2)

# ========== 2. Wall time vs subB lines, with bedtools reference ==========
# Bedtools is a single fixed cost; only mixed-stride's wall time drops with
# subB. Line plot makes that contrast immediate.
ms_rec <- maur_sum %>% filter(method == "mixed-stride", subB == 0.1)
xs <- sort(unique(maur_sum$subB))

p2 <- ggplot(maur_sum, aes(x = subB, y = wall_median, color = method)) +
  geom_hline(yintercept = bt_wall, linetype = "dashed",
             color = bedtools_col, linewidth = 0.4) +
  annotate("text", x = min(xs), y = bt_wall,
           label = sprintf("bedtools (%.1fs)", bt_wall),
           hjust = 0, vjust = -0.6, color = bedtools_col,
           fontface = "italic", size = 3.8) +
  geom_line(linewidth = 0.5, alpha = 0.95) +
  geom_point(size = 3.4) +
  geom_point(data = ms_rec,
             shape = 21, size = 8, stroke = 0.8,
             color = method_pal[["mixed-stride"]], fill = NA) +
  annotate("text", x = ms_rec$subB * 0.85, y = ms_rec$wall_median,
           label = sprintf("mixed-stride @ subB=0.1\n%.1fs (%.2f× vs bedtools)",
                           ms_rec$wall_median, bt_wall / ms_rec$wall_median),
           hjust = 1, vjust = 0.5, lineheight = 1.05,
           fontface = "italic", size = 3.6,
           color = method_pal[["mixed-stride"]]) +
  scale_color_manual(values = method_pal, name = NULL) +
  scale_x_log10(breaks = xs,
                labels = function(x) format(x, drop0trailing = TRUE)) +
  scale_y_continuous(limits = c(0, NA),
                     breaks = seq(0, 16, 2),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Bedtools is fixed cost; only mixed-stride drops with subB",
       subtitle = "Wall time on the Maurano corpus (20 BEDs, 190 pairs, 8 threads). Sort time not charged to bedtools.",
       x = "subB (log scale; 1.0 = no subsampling)",
       y = "wall time (s)",
       caption = "Hash-threshold and single-hash barely move — their per-point gate cost masks the savings from fewer HLL ingests. Mixed-stride's chr-keyed deterministic stride scales naturally with subB.") +
  theme_pub
save_png(file.path(figures_dir, "headline_hammock_vs_bedtools.png"),
         p2, width = 10, height = 6)

# ========== 3. Mixed-stride scaling across synthetic sizes ==========
# Show how speedup grows with corpus size at fixed subB=0.1.
scale_df <- syn_sum %>%
  filter(subB == 0.1) %>%
  mutate(size_class = factor(size_class, levels = c("10k", "100k", "1M")))

p3 <- ggplot(scale_df, aes(x = size_class, y = speedup_nosub,
                            fill = method)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  geom_text(aes(label = sprintf("%.2f×", speedup_nosub)),
            position = position_dodge(width = 0.78),
            vjust = -0.4, size = 3.3, fontface = "bold") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = method_pal, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Mixed-stride is the only method whose speedup tracks subB",
       subtitle = "Speedup vs hammock at subB=1.0, synthetic random BEDs, subB=0.1, Mode B p=18, 8 threads.",
       x = "intervals per BED (6 files per class, all-vs-all)",
       y = "speedup vs no subsample",
       caption = "Hash-threshold and single-hash plateau at ~1.1× regardless of corpus size; mixed-stride scales from 3× to 4× as files grow.") +
  theme_pub
save_png(file.path(figures_dir, "headline_synthetic_scaling.png"),
         p3, width = 9, height = 5)

cat("\nWrote headline_* figures into", figures_dir, "\n")
