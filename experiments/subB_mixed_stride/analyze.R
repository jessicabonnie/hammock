#!/usr/bin/env Rscript
# Analyze the subB / mixed-stride sweep.
#
# Inputs:  results/sweep_<stamp>.csv (long-form; one row per pair × subB × rep)
# Outputs: figures/wall_time_vs_subB.png, figures/jaccard_error_vs_subB.png,
#          figures/summary_grid.png, results/summary_<stamp>.csv
#
# Ground truth: hammock Mode B with subB = 1.0 (no point subsampling), same
# mixed-stride config, same pair. Per-pair accuracy = |J(subB=x) - J(subB=1.0)|.
#
# Usage:
#   ml r/4.3.0
#   Rscript analyze.R [path/to/sweep.csv]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(Cairo)
})

# Headless node: base R has no png()/cairo support; route ggsave through Cairo::CairoPNG
save_png <- function(path, plot, width = 7, height = 4.5, dpi = 150) {
  CairoPNG(filename = path, width = width, height = height,
           units = "in", res = dpi)
  on.exit(dev.off())
  print(plot)
}

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(sub("--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
results_dir <- file.path(script_dir, "results")
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

if (length(args) >= 1) {
  in_csv <- args[1]
} else {
  csvs <- sort(list.files(results_dir, pattern = "^sweep_.*\\.csv$",
                          full.names = TRUE), decreasing = TRUE)
  if (length(csvs) == 0) stop("No sweep CSV found in ", results_dir)
  in_csv <- csvs[1]
}
cat("Reading:", in_csv, "\n")
df <- read_csv(in_csv, show_col_types = FALSE)

stamp <- sub("^sweep_(.*)\\.csv$", "\\1", basename(in_csv))

size_levels <- c("10k", "100k", "1M")
df$size_class <- factor(df$size_class, levels = size_levels)

# ---- accuracy: per-pair error vs subB=1.0 ground truth ----
gt <- df %>%
  filter(subB == 1.0) %>%
  group_by(size_class, file_a, file_b) %>%
  summarise(j_truth = mean(jaccard), .groups = "drop")

acc <- df %>%
  inner_join(gt, by = c("size_class", "file_a", "file_b")) %>%
  mutate(abs_err = abs(jaccard - j_truth),
         signed_err = jaccard - j_truth)

acc_summary <- acc %>%
  group_by(size_class, subB) %>%
  summarise(
    n_pairs    = n(),
    mae        = mean(abs_err),
    max_err    = max(abs_err),
    median_err = median(abs_err),
    mean_j     = mean(jaccard),
    mean_truth = mean(j_truth),
    .groups = "drop"
  )

# ---- timing: one row per run (collapse the 15 per-pair duplicates) ----
runs <- df %>%
  distinct(size_class, subB, rep, run_id, wall_time, cpu_time, max_rss_mb,
           sketch_creation_time, comparison_time)

speed_summary <- runs %>%
  group_by(size_class, subB) %>%
  summarise(
    n_reps          = n(),
    wall_median     = median(wall_time),
    wall_mean       = mean(wall_time),
    wall_sd         = if (n() > 1) sd(wall_time) else 0,
    cpu_median      = median(cpu_time),
    rss_median_mb   = median(max_rss_mb),
    sketch_median_s = median(sketch_creation_time, na.rm = TRUE),
    pair_median_s   = median(comparison_time,   na.rm = TRUE),
    .groups = "drop"
  )

# speedup vs subB=1.0 baseline
baseline <- speed_summary %>%
  filter(subB == 1.0) %>%
  select(size_class, wall_baseline = wall_median)
speed_summary <- speed_summary %>%
  left_join(baseline, by = "size_class") %>%
  mutate(speedup = wall_baseline / wall_median)

# write joined summary
summary_path <- file.path(results_dir, paste0("summary_", stamp, ".csv"))
joined <- speed_summary %>% left_join(acc_summary, by = c("size_class", "subB"))
write_csv(joined, summary_path)
cat("Wrote:", summary_path, "\n")
cat("\n--- Summary ---\n")
print(as.data.frame(joined %>% select(size_class, subB, wall_median, speedup,
                                       mae, max_err, n_pairs)),
      row.names = FALSE, digits = 4)

# ---- plots ----
theme_set(theme_bw(base_size = 12))
size_pal <- c("10k" = "#1b9e77", "100k" = "#d95f02", "1M" = "#7570b3")

p_speed <- ggplot(runs, aes(x = subB, y = wall_time, color = size_class)) +
  geom_jitter(width = 0.02, height = 0, alpha = 0.45, size = 1.6) +
  stat_summary(fun = median, geom = "line", linewidth = 0.9) +
  stat_summary(fun = median, geom = "point", size = 2.4) +
  scale_x_log10(breaks = sort(unique(runs$subB)),
                labels = function(x) format(x, drop0trailing = TRUE)) +
  scale_y_log10() +
  scale_color_manual(values = size_pal) +
  labs(title = "Mode B wall time vs subB (mixed-stride on)",
       subtitle = "Median across replicates; dots are individual runs",
       x = "subB (log scale)", y = "wall time (s, log)",
       color = "size class")
save_png(file.path(figures_dir, "wall_time_vs_subB.png"), p_speed,
       width = 7, height = 4.5, dpi = 150)

p_speedup <- ggplot(speed_summary, aes(x = subB, y = speedup, color = size_class)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
  geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
  scale_x_log10(breaks = sort(unique(speed_summary$subB)),
                labels = function(x) format(x, drop0trailing = TRUE)) +
  scale_color_manual(values = size_pal) +
  labs(title = "Speedup vs subB=1.0 baseline",
       x = "subB (log scale)", y = "wall_time(subB=1.0) / wall_time(subB)",
       color = "size class")
save_png(file.path(figures_dir, "speedup_vs_subB.png"), p_speedup,
       width = 7, height = 4.5, dpi = 150)

p_acc <- ggplot(acc %>% filter(subB < 1.0),
                aes(x = factor(subB), y = abs_err, fill = size_class)) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.4, alpha = 0.6,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = size_pal) +
  labs(title = "Jaccard absolute error vs subB=1.0 ground truth",
       subtitle = "Per pair, across replicates",
       x = "subB", y = "|J(subB) - J(subB=1.0)|",
       fill = "size class")
save_png(file.path(figures_dir, "jaccard_error_vs_subB.png"), p_acc,
       width = 7, height = 4.5, dpi = 150)

p_mae <- ggplot(acc_summary %>% filter(subB < 1.0),
                aes(x = subB, y = mae, color = size_class)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
  scale_x_log10(breaks = sort(unique(acc_summary$subB[acc_summary$subB < 1])),
                labels = function(x) format(x, drop0trailing = TRUE)) +
  scale_color_manual(values = size_pal) +
  labs(title = "Mean absolute Jaccard error vs subB",
       x = "subB (log scale)", y = "MAE vs subB=1.0",
       color = "size class")
save_png(file.path(figures_dir, "mae_vs_subB.png"), p_mae,
       width = 7, height = 4.5, dpi = 150)

cat("\nFigures written to:", figures_dir, "\n")
