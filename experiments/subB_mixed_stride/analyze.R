#!/usr/bin/env Rscript
# Analyze the subB / subB-method sweep on one or more corpora.
#
# Inputs:  results/sweep_<corpus>_<stamp>.csv (one row per pair × method ×
#          subB × rep). Legacy sweep_<stamp>.csv files without a corpus tag
#          are treated as synthetic.
# Outputs: figures/<corpus>_*.png and results/summary_<corpus>_<stamp>.csv
#
# Ground truth: hammock Mode B at subB = 1.0 (no point subsampling), averaged
# across methods. Self-consistent — isolates subsample-induced drift.
#
# Usage:
#   ml r/4.3.0
#   Rscript analyze.R                       # process latest CSV per corpus
#   Rscript analyze.R results/sweep_X.csv   # process this CSV only

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(Cairo)
})

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

method_levels <- c("hash-threshold", "mixed-stride", "single-hash")
method_pal <- c(
  "hash-threshold" = "#1b9e77",
  "mixed-stride"   = "#d95f02",
  "single-hash"    = "#7570b3"
)

# ---------- corpus dispatch ----------
# If the user passed an explicit CSV, use it; otherwise pick the latest per
# corpus (sweep_synthetic_*.csv, sweep_maurano_*.csv, and legacy sweep_*.csv
# which we treat as synthetic).
pick_latest <- function(pattern) {
  files <- sort(list.files(results_dir, pattern = pattern,
                           full.names = TRUE), decreasing = TRUE)
  if (length(files) == 0) return(NULL)
  files[1]
}

if (length(args) >= 1) {
  job_list <- list(args[1])
} else {
  jobs <- list()
  # The [0-9] is load-bearing: pick_latest sorts decreasing, so a bare
  # `^sweep_synthetic_.*` would adopt any sweep_synthetic_<word>_* file (e.g. an
  # _ie or _bedtools side-run) as *the* synthetic sweep, since 'i'/'b' > '2'.
  # Mirrors the maurano pattern below, which already required a digit.
  syn <- pick_latest("^sweep_synthetic_[0-9].*\\.csv$")
  if (is.null(syn)) syn <- pick_latest("^sweep_[0-9].*\\.csv$")   # legacy
  if (!is.null(syn)) jobs <- c(jobs, syn)
  maur <- pick_latest("^sweep_maurano_[0-9].*\\.csv$")
  if (!is.null(maur)) jobs <- c(jobs, maur)
  if (length(jobs) == 0) stop("No sweep CSV found in ", results_dir)
  job_list <- jobs
}

# Bedtools reference (optional). One median wall time per corpus, used to draw
# a horizontal reference line on the speed/speedup plots.
load_bedtools_ref <- function(corpus) {
  pattern <- sprintf("^sweep_%s_bedtools_.*\\.csv$", corpus)
  path <- pick_latest(pattern)
  if (is.null(path)) return(NULL)
  bt <- read_csv(path, show_col_types = FALSE)
  bt %>%
    distinct(rep, run_id, wall_time) %>%
    summarise(wall_median = median(wall_time),
              wall_min    = min(wall_time),
              wall_max    = max(wall_time),
              n_reps      = n()) %>%
    as.list()
}

analyze_one <- function(in_csv) {
  cat("\n========== ", in_csv, " ==========\n", sep = "")
  df <- read_csv(in_csv, show_col_types = FALSE)

  # Back-fill corpus column for legacy CSVs that pre-date the dimension.
  if (!"corpus" %in% names(df)) df$corpus <- "synthetic"
  corpus <- unique(df$corpus)
  if (length(corpus) != 1) {
    cat("WARN: multiple corpora in one CSV; keeping first =", corpus[1], "\n")
    corpus <- corpus[1]
    df <- df %>% filter(corpus == !!corpus)
  }
  stamp <- sub("^sweep_.*?_?([0-9]{8}_[0-9]{6})\\.csv$", "\\1", basename(in_csv))

  size_levels <- unique(as.character(df$size_class))
  # Preserve a sensible order for synthetic; pass through whatever else
  if (all(c("10k", "100k", "1M") %in% size_levels)) {
    size_levels <- c("10k", "100k", "1M")
  }
  df$size_class <- factor(df$size_class, levels = size_levels)
  df$method     <- factor(df$method,     levels = method_levels)

  # ---- ground truth at subB=1.0, averaged across methods ----
  gt <- df %>%
    filter(subB == 1.0) %>%
    group_by(size_class, file_a, file_b) %>%
    summarise(j_truth = mean(jaccard), .groups = "drop")

  acc <- df %>%
    inner_join(gt, by = c("size_class", "file_a", "file_b")) %>%
    mutate(abs_err    = abs(jaccard - j_truth),
           signed_err = jaccard - j_truth)

  acc_summary <- acc %>%
    group_by(method, size_class, subB) %>%
    summarise(
      n_pairs    = n(),
      mae        = mean(abs_err),
      max_err    = max(abs_err),
      median_err = median(abs_err),
      .groups = "drop"
    )

  # ---- timing: collapse the per-pair duplicates ----
  runs <- df %>%
    distinct(method, size_class, subB, rep, run_id, wall_time, cpu_time,
             max_rss_mb, sketch_creation_time, comparison_time)

  speed_summary <- runs %>%
    group_by(method, size_class, subB) %>%
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

  baseline_nosub <- runs %>%
    filter(subB == 1.0) %>%
    group_by(size_class) %>%
    summarise(wall_nosub = median(wall_time), .groups = "drop")

  baseline_within <- speed_summary %>%
    filter(subB == 1.0) %>%
    select(method, size_class, wall_within = wall_median)

  speed_summary <- speed_summary %>%
    left_join(baseline_nosub,    by = "size_class") %>%
    left_join(baseline_within,   by = c("method", "size_class")) %>%
    mutate(
      speedup_vs_nosub  = wall_nosub  / wall_median,
      speedup_vs_within = wall_within / wall_median
    )

  joined <- speed_summary %>%
    left_join(acc_summary, by = c("method", "size_class", "subB"))
  summary_path <- file.path(results_dir,
                            paste0("summary_", corpus, "_", stamp, ".csv"))
  write_csv(joined, summary_path)

  # Optional bedtools reference for this corpus
  bt <- load_bedtools_ref(corpus)
  if (!is.null(bt)) {
    cat(sprintf("bedtools ref (%s): median %.2fs across %d reps\n",
                corpus, bt$wall_median, bt$n_reps))
  }
  cat("Wrote:", summary_path, "\n")
  cat("\n--- Summary (", corpus, ") ---\n", sep = "")
  print(as.data.frame(joined %>%
                       select(method, size_class, subB,
                              wall_median, speedup_vs_nosub,
                              mae, max_err)),
        row.names = FALSE, digits = 4)

  # ---- plots, prefixed by corpus ----
  prefix <- paste0(corpus, "_")
  # Synthetic has multi-class facets; maurano has one. Choose a width that
  # works for either.
  fw <- if (length(size_levels) > 1) 10 else 6.5
  fh <- 4.2

  p_speed <- ggplot(runs, aes(x = subB, y = wall_time, color = method)) +
    geom_jitter(width = 0.04, height = 0, alpha = 0.35, size = 1.2) +
    stat_summary(fun = median, geom = "line", linewidth = 0.9) +
    stat_summary(fun = median, geom = "point", size = 2.2) +
    scale_x_log10(breaks = sort(unique(runs$subB)),
                  labels = function(x) format(x, drop0trailing = TRUE)) +
    scale_y_log10() +
    scale_color_manual(values = method_pal) +
    facet_wrap(~ size_class, scales = "free_y", labeller = label_both) +
    labs(title = paste0("Mode B wall time vs subB (", corpus, ")"),
         subtitle = "Median across replicates; dots are individual runs",
         x = "subB (log scale)", y = "wall time (s, log)",
         color = "subB-method")
  if (!is.null(bt)) {
    p_speed <- p_speed +
      geom_hline(yintercept = bt$wall_median, linetype = "dashed",
                 color = "black") +
      annotate("text", x = min(runs$subB), y = bt$wall_median,
               label = sprintf("bedtools (%.1fs)", bt$wall_median),
               hjust = 0, vjust = -0.4, size = 3, color = "black")
  }
  save_png(file.path(figures_dir, paste0(prefix, "wall_time_vs_subB.png")),
           p_speed, width = fw, height = fh)

  p_sp_nosub <- ggplot(speed_summary,
                       aes(x = subB, y = speedup_vs_nosub, color = method)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
    geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
    scale_x_log10(breaks = sort(unique(speed_summary$subB)),
                  labels = function(x) format(x, drop0trailing = TRUE)) +
    scale_color_manual(values = method_pal) +
    facet_wrap(~ size_class, labeller = label_both) +
    labs(title = paste0("Speedup vs Mode B without subsample (", corpus, ")"),
         subtitle = "Baseline = median wall time across methods at subB=1.0",
         x = "subB (log scale)",
         y = "wall_time(no subsample) / wall_time(method, subB)",
         color = "subB-method")
  if (!is.null(bt)) {
    # bedtools speedup = wall_nosub / bedtools_wall, per size_class
    bt_lines <- baseline_nosub %>%
      mutate(bt_speedup = wall_nosub / bt$wall_median)
    p_sp_nosub <- p_sp_nosub +
      geom_hline(data = bt_lines, aes(yintercept = bt_speedup),
                 linetype = "dashed", color = "black") +
      geom_text(data = bt_lines,
                aes(x = min(speed_summary$subB), y = bt_speedup,
                    label = sprintf("bedtools (%.2f×)", bt_speedup)),
                hjust = 0, vjust = -0.4, size = 3, color = "black",
                inherit.aes = FALSE)
  }
  save_png(file.path(figures_dir, paste0(prefix, "speedup_vs_nosub.png")),
           p_sp_nosub, width = fw, height = fh)

  p_sp_within <- ggplot(speed_summary,
                        aes(x = subB, y = speedup_vs_within, color = method)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
    geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
    scale_x_log10(breaks = sort(unique(speed_summary$subB)),
                  labels = function(x) format(x, drop0trailing = TRUE)) +
    scale_color_manual(values = method_pal) +
    facet_wrap(~ size_class, labeller = label_both) +
    labs(title = paste0("Speedup vs each method's own subB=1.0 baseline (", corpus, ")"),
         x = "subB (log scale)",
         y = "wall_time(method, 1.0) / wall_time(method, subB)",
         color = "subB-method")
  save_png(file.path(figures_dir, paste0(prefix, "speedup_within_method.png")),
           p_sp_within, width = fw, height = fh)

  p_mae <- ggplot(acc_summary %>% filter(subB < 1.0),
                  aes(x = subB, y = mae, color = method)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
    scale_x_log10(breaks = sort(unique(acc_summary$subB[acc_summary$subB < 1])),
                  labels = function(x) format(x, drop0trailing = TRUE)) +
    scale_color_manual(values = method_pal) +
    facet_wrap(~ size_class, scales = "free_y", labeller = label_both) +
    labs(title = paste0("Mean absolute Jaccard error vs subB (", corpus, ")"),
         subtitle = "Ground truth = subB=1.0, averaged across methods",
         x = "subB (log scale)", y = "MAE",
         color = "subB-method")
  save_png(file.path(figures_dir, paste0(prefix, "mae_vs_subB.png")),
           p_mae, width = fw, height = fh)

  p_acc <- ggplot(acc %>% filter(subB < 1.0),
                  aes(x = factor(subB), y = abs_err, fill = method)) +
    geom_boxplot(outlier.size = 0.5, linewidth = 0.3, alpha = 0.6,
                 position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = method_pal) +
    facet_wrap(~ size_class, scales = "free_y", labeller = label_both) +
    labs(title = paste0("Jaccard absolute error vs subB=1.0 ground truth (", corpus, ")"),
         subtitle = "Per pair, across replicates",
         x = "subB", y = "|J(method, subB) - J(no subsample)|",
         fill = "subB-method")
  save_png(file.path(figures_dir, paste0(prefix, "jaccard_error_vs_subB.png")),
           p_acc, width = fw, height = 4.5)

  cat("Figures written to:", figures_dir, "(prefix:", prefix, ")\n")
}

for (csv in job_list) analyze_one(csv)
