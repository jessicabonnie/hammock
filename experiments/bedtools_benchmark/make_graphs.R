#!/usr/bin/env Rscript
# Generate paper-quality summary graphs from the bedtools_benchmark CSVs.
# R port of make_graphs.py — same inputs, same output filenames, ggplot2-based.
#
# Usage:
#   Rscript make_graphs.R --files-csv <path>                     # files-axis plots
#   Rscript make_graphs.R --pairs-csv <path>                     # jaccard plots
#   Rscript make_graphs.R --files-csv <p1> --pairs-csv <p2> [--out-dir <d>]
#
# Requires r/4.3.0 module (`ml r/4.3.0`) and: ggplot2, scales, dplyr, readr, tidyr.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(Cairo)  # headless PNG device (this R build lacks native png/cairo)
})

# Headless PNG: use Cairo's device through ggsave.
ggsave_png <- function(filename, plot, width, height, dpi = 200) {
  CairoPNG(filename, width = width, height = height, units = "in",
           res = dpi, bg = "white")
  on.exit(dev.off(), add = TRUE)
  print(plot)
}

SCRIPT_DIR  <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=^--file=).*", args, perl = TRUE))
  if (length(m) == 0) getwd() else dirname(normalizePath(m[1]))
}, error = function(e) getwd())
FIGURES_DIR <- file.path(SCRIPT_DIR, "figures")

# ---------- arg parsing ----------

parse_args <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  out  <- list(files_csv = NULL, pairs_csv = NULL, out_dir = FIGURES_DIR)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--files-csv") { out$files_csv <- argv[i + 1]; i <- i + 2; next }
    if (a == "--pairs-csv") { out$pairs_csv <- argv[i + 1]; i <- i + 2; next }
    if (a == "--out-dir")   { out$out_dir   <- argv[i + 1]; i <- i + 2; next }
    if (a %in% c("-h", "--help")) {
      cat("Usage: Rscript make_graphs.R --files-csv <path> [--pairs-csv <path>] [--out-dir <dir>]\n")
      quit(save = "no", status = 0)
    }
    stop(sprintf("Unknown argument: %s", a))
    i <- i + 1
  }
  if (is.null(out$files_csv) && is.null(out$pairs_csv)) {
    stop("supply at least one of --files-csv or --pairs-csv")
  }
  out
}

# ---------- files-axis plots ----------

load_files_csv <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  bt <- df %>% filter(tool == "bedtools") %>%
    select(num_files, num_threads, precision,
           bt_wall = mean_wall_time, bt_wall_err = std_wall_time,
           sort_wall = mean_sort_time)
  hm <- df %>% filter(tool == "hammock_cpp_B") %>%
    select(num_files,
           hm_wall = mean_wall_time, hm_wall_err = std_wall_time,
           hm_sketch = mean_sketch_creation_time,
           hm_compare = mean_comparison_time)
  inner_join(bt, hm, by = "num_files") %>% arrange(num_files)
}

# Persistent crossover: largest i where bt[i] < hm[i] and bt[i+1] > hm[i+1];
# return the first N at which hammock wins persistently (the "+1" point).
find_persistent_crossover <- function(d) {
  idx <- NA_integer_
  for (i in seq_len(nrow(d) - 1)) {
    if (d$bt_wall[i] < d$hm_wall[i] && d$bt_wall[i + 1] > d$hm_wall[i + 1]) {
      idx <- i + 1L
    }
  }
  idx
}

plot_files_sketch_compare_split <- function(files_csv_path, out_path) {
  d <- load_files_csv(files_csv_path)
  t <- d$num_threads[1]; p <- d$precision[1]

  # Compare values can be 0 at small N; clip for log axis.
  d$hm_compare_plot <- pmax(d$hm_compare, 1e-4)

  series <- bind_rows(
    d %>% transmute(num_files, y = bt_wall,        yerr = bt_wall_err,
                    series = sprintf("bedtools wall (parallel-wrapped, t=%d)", t)),
    d %>% transmute(num_files, y = hm_sketch,      yerr = hm_wall_err,
                    series = sprintf("hammock sketching (t=%d, p=%d) — dominates hammock wall", t, p)),
    d %>% transmute(num_files, y = hm_compare_plot, yerr = NA_real_,
                    series = "hammock pairwise compare only"),
    d %>% transmute(num_files, y = sort_wall,       yerr = NA_real_,
                    series = "bedtools pre-sort (excluded from bedtools wall)")
  )

  series_levels <- c(
    sprintf("bedtools wall (parallel-wrapped, t=%d)", t),
    sprintf("hammock sketching (t=%d, p=%d) — dominates hammock wall", t, p),
    "hammock pairwise compare only",
    "bedtools pre-sort (excluded from bedtools wall)"
  )
  series$series <- factor(series$series, levels = series_levels)

  color_map <- setNames(c("#1f77b4", "#ff7f0e", "#d62728", "#7f7f7f"), series_levels)
  lty_map   <- setNames(c("solid",   "solid",   "solid",   "dotted"),  series_levels)
  shape_map <- setNames(c(16,        15,        18,        4),         series_levels)

  cross_idx <- find_persistent_crossover(d)

  g <- ggplot(series, aes(x = num_files, y = y,
                          color = series, linetype = series, shape = series))
  if (!is.na(cross_idx)) {
    g <- g + geom_vline(xintercept = d$num_files[cross_idx],
                        color = "grey", linetype = "dashed", alpha = 0.4)
  }
  g <- g +
    geom_errorbar(data = subset(series, !is.na(yerr)),
                  aes(ymin = pmax(y - yerr, 1e-6), ymax = y + yerr),
                  width = 0.05, alpha = 0.7, show.legend = FALSE) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.4) +
    scale_color_manual(values = color_map, name = NULL) +
    scale_linetype_manual(values = lty_map, name = NULL) +
    scale_shape_manual(values = shape_map, name = NULL) +
    scale_x_continuous(
      trans = log2_trans(),
      breaks = d$num_files,
      labels = d$num_files,
      sec.axis = sec_axis(~ . ^ 2, name = "Pairwise compares (N²)",
                          breaks = d$num_files ^ 2,
                          labels = function(x) format(x, big.mark = ",", scientific = FALSE))
    ) +
    scale_y_continuous(trans = log10_trans(), labels = label_log()) +
    labs(
      x = "Number of files (N) — pairwise compares = N²",
      y = "Wall time (s)",
      title = "Hammock vs bedtools: where the time goes",
      subtitle = sprintf("p=%d, t=%d, 10 000 intervals/file, 3 runs/config", p, t)
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
      legend.key = element_rect(fill = NA),
      panel.grid.minor = element_line(linewidth = 0.2, color = "grey90"),
      plot.title = element_text(face = "bold")
    )

  if (!is.na(cross_idx)) {
    ncross <- d$num_files[cross_idx]
    g <- g + annotate("text",
                      x = ncross * 0.55, y = d$bt_wall[cross_idx] * 4.5,
                      label = sprintf("persistent crossover\nat N=%d", ncross),
                      color = "grey40", size = 3.2, hjust = 1)
  }

  # Annotate the compare-vs-sketch ratio at the rightmost N.
  n_last <- tail(d$num_files, 1)
  cmp_last <- tail(d$hm_compare, 1)
  sk_last  <- tail(d$hm_sketch, 1)
  cmp_plot <- tail(d$hm_compare_plot, 1)
  ratio <- if (sk_last > 0) cmp_last / sk_last else 0
  g <- g + annotate("text",
                    x = n_last * 0.6, y = cmp_plot * 30,
                    label = sprintf("compare = %.3fs\nvs sketch %.1fs\n(%.2f%% of sketch)",
                                    cmp_last, sk_last, ratio * 100),
                    color = "#d62728", size = 3.0, hjust = 1)

  ggsave_png(out_path, g, width = 9, height = 6, dpi = 200)
  cat(sprintf("Wrote %s\n", out_path))
}

plot_cost_per_pair <- function(files_csv_path, out_path) {
  d <- load_files_csv(files_csv_path)
  t <- d$num_threads[1]; p <- d$precision[1]
  d <- d %>% mutate(
    N2     = num_files ^ 2,
    bt_pp  = bt_wall / N2,
    hm_pp  = hm_wall / N2,
    sk_pp  = hm_sketch / N2,
    cmp_pp = pmax(hm_compare / N2, 1e-9)
  )

  series <- bind_rows(
    d %>% transmute(num_files, y = bt_pp,
                    series = sprintf("bedtools (parallel-wrapped, t=%d)", t)),
    d %>% transmute(num_files, y = hm_pp,
                    series = sprintf("hammock total (t=%d, p=%d)", t, p)),
    d %>% transmute(num_files, y = sk_pp,
                    series = "hammock sketching"),
    d %>% transmute(num_files, y = cmp_pp,
                    series = "hammock pairwise compare")
  )
  series_levels <- c(
    sprintf("bedtools (parallel-wrapped, t=%d)", t),
    sprintf("hammock total (t=%d, p=%d)", t, p),
    "hammock sketching",
    "hammock pairwise compare"
  )
  series$series <- factor(series$series, levels = series_levels)

  color_map <- setNames(c("#1f77b4", "#ff7f0e", "#ff7f0e", "#d62728"), series_levels)
  lty_map   <- setNames(c("solid",   "solid",   "dashed",  "solid"),   series_levels)
  shape_map <- setNames(c(16,        15,        17,        18),        series_levels)
  alpha_map <- setNames(c(1, 1, 0.55, 1), series_levels)

  g <- ggplot(series, aes(x = num_files, y = y,
                          color = series, linetype = series, shape = series, alpha = series)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.4) +
    scale_color_manual(values = color_map, name = NULL) +
    scale_linetype_manual(values = lty_map, name = NULL) +
    scale_shape_manual(values = shape_map, name = NULL) +
    scale_alpha_manual(values = alpha_map, name = NULL) +
    scale_x_continuous(trans = log2_trans(),
                       breaks = d$num_files, labels = d$num_files) +
    scale_y_continuous(trans = log10_trans(), labels = label_log()) +
    labs(
      x = "Number of files (N)",
      y = "Wall time per pairwise compare = wall / N²  (s)",
      title = "Cost per pair: bedtools is flat, hammock amortizes",
      subtitle = sprintf("p=%d, t=%d, 10 000 intervals/file, 3 runs/config", p, t)
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
      panel.grid.minor = element_line(linewidth = 0.2, color = "grey90"),
      plot.title = element_text(face = "bold")
    )

  ggsave_png(out_path, g, width = 9, height = 6, dpi = 200)
  cat(sprintf("Wrote %s\n", out_path))
}

# ---------- jaccard pair plots ----------

load_pairs_csv <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      bedtools_jaccard = suppressWarnings(as.numeric(bedtools_jaccard)),
      hammock_jaccard  = suppressWarnings(as.numeric(hammock_jaccard)),
      precision        = as.integer(precision)
    ) %>%
    filter(!is.na(bedtools_jaccard), !is.na(hammock_jaccard))
}

plot_jaccard_scatter <- function(pairs_csv_path, out_path) {
  d <- load_pairs_csv(pairs_csv_path)
  if (nrow(d) == 0) { warning("No pairs found in ", pairs_csv_path); return(invisible(NULL)) }

  counts <- d %>% group_by(precision) %>% summarise(n = n(), .groups = "drop")
  d <- d %>% left_join(counts, by = "precision") %>%
    mutate(label = sprintf("p=%d  (n=%d)", precision, n))
  legend_order <- counts %>% arrange(precision) %>%
    mutate(label = sprintf("p=%d  (n=%d)", precision, n)) %>% pull(label)
  d$label <- factor(d$label, levels = legend_order)

  xpad <- max(0.005, diff(range(d$bedtools_jaccard)) * 0.15)
  ypad <- max(0.005, diff(range(d$hammock_jaccard)) * 0.15)
  xlim <- range(d$bedtools_jaccard) + c(-xpad, xpad)
  ylim <- range(d$hammock_jaccard)  + c(-ypad, ypad)
  lo <- min(xlim[1], ylim[1]); hi <- max(xlim[2], ylim[2])

  bt_med <- median(d$bedtools_jaccard)
  hm_med <- median(d$hammock_jaccard)

  g <- ggplot(d, aes(x = bedtools_jaccard, y = hammock_jaccard, color = label)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
    geom_point(size = 1.0, alpha = 0.35) +
    annotate("segment", x = bt_med, xend = bt_med, y = bt_med, yend = hm_med,
             arrow = arrow(length = unit(0.18, "cm")), color = "black", linewidth = 0.7) +
    annotate("text",
             x = bt_med + (hi - bt_med) * 0.10,
             y = (bt_med + hm_med) / 2,
             label = sprintf("median gap = %.3f\n(bedtools med = %.3f,\n hammock med = %.3f)",
                             hm_med - bt_med, bt_med, hm_med),
             hjust = 0, size = 3.0) +
    scale_color_viridis_d(name = NULL, end = 0.95,
                          guide = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
    coord_fixed(xlim = c(lo, hi), ylim = c(lo, hi)) +
    labs(
      x = "bedtools jaccard  (bp set-Jaccard, ground truth)",
      y = "hammock jaccard  (register-equality Jaccard)",
      title = "Definitional gap (zoomed to data cluster)",
      subtitle = "Points above the y=x diagonal = hammock overestimates set-Jaccard"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
      plot.title = element_text(face = "bold")
    )

  # Diagonal y=x label sits along the line outside the legend
  g <- g + annotate("text", x = lo + (hi - lo) * 0.08, y = lo + (hi - lo) * 0.04,
                    label = "y = x  (set-Jaccard equality)", angle = 45,
                    color = "black", alpha = 0.7, size = 3.0)

  ggsave_png(out_path, g, width = 8.5, height = 8, dpi = 200)
  cat(sprintf("Wrote %s\n", out_path))
}

plot_jaccard_delta <- function(pairs_csv_path, out_path) {
  d <- load_pairs_csv(pairs_csv_path)
  if (nrow(d) == 0) { warning("No pairs found in ", pairs_csv_path); return(invisible(NULL)) }
  d <- d %>% mutate(delta = hammock_jaccard - bedtools_jaccard)

  counts <- d %>% group_by(precision) %>% summarise(n = n(), .groups = "drop")
  d <- d %>% left_join(counts, by = "precision") %>%
    mutate(label = sprintf("p=%d  (n=%d)", precision, n))
  legend_order <- counts %>% arrange(precision) %>%
    mutate(label = sprintf("p=%d  (n=%d)", precision, n)) %>% pull(label)
  d$label <- factor(d$label, levels = legend_order)

  bins <- seq(min(d$bedtools_jaccard), max(d$bedtools_jaccard), length.out = 12)
  binned <- d %>%
    mutate(bin = cut(bedtools_jaccard, breaks = bins, include.lowest = TRUE)) %>%
    group_by(precision, label, bin) %>%
    summarise(
      n = n(),
      x = mean(c(bins[as.integer(bin)], bins[as.integer(bin) + 1])),
      mean_delta = mean(delta),
      .groups = "drop"
    ) %>% filter(n >= 3) %>% arrange(precision, x)

  med_gap <- median(d$delta)

  g <- ggplot(d, aes(x = bedtools_jaccard, y = delta, color = label)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8, alpha = 0.7) +
    geom_hline(yintercept = med_gap, color = "red", linewidth = 0.6,
               linetype = "dotted", alpha = 0.7) +
    geom_point(size = 0.9, alpha = 0.25) +
    geom_line(data = binned, aes(x = x, y = mean_delta, color = label, group = label),
              linewidth = 1.0) +
    annotate("text", x = max(d$bedtools_jaccard), y = med_gap,
             label = sprintf("  median gap = %.3f", med_gap),
             color = "red", hjust = 1, vjust = -0.5, size = 3.0) +
    scale_color_viridis_d(name = NULL, end = 0.95,
                          guide = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
    labs(
      x = "bedtools jaccard  (bp set-Jaccard, ground truth)",
      y = "hammock − bedtools  (register-equality − set-Jaccard)",
      title = "Definitional gap, isolated: hammock systematically overestimates set Jaccard",
      subtitle = "Lines = per-precision running mean; gap is roughly constant across p (not HLL noise)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
      plot.title = element_text(face = "bold")
    )

  ggsave_png(out_path, g, width = 10, height = 6, dpi = 200)
  cat(sprintf("Wrote %s\n", out_path))
}

# ---------- main ----------

main <- function() {
  opts <- parse_args()
  dir.create(opts$out_dir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(opts$files_csv)) {
    stem <- tools::file_path_sans_ext(basename(opts$files_csv))
    plot_files_sketch_compare_split(
      opts$files_csv,
      file.path(opts$out_dir, paste0(stem, "_sketch_compare_split.png"))
    )
    plot_cost_per_pair(
      opts$files_csv,
      file.path(opts$out_dir, paste0(stem, "_cost_per_pair.png"))
    )
  }

  if (!is.null(opts$pairs_csv)) {
    stem <- tools::file_path_sans_ext(basename(opts$pairs_csv))
    plot_jaccard_scatter(
      opts$pairs_csv,
      file.path(opts$out_dir, paste0(stem, "_jaccard_scatter.png"))
    )
    plot_jaccard_delta(
      opts$pairs_csv,
      file.path(opts$out_dir, paste0(stem, "_jaccard_delta.png"))
    )
  }
}

main()
