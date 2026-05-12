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
  # Backwards compat: pre-subB CSVs have no sub_b column.
  if (!"sub_b" %in% names(df)) df$sub_b <- ifelse(df$tool == "bedtools", NA_real_, 1.0)
  # In the new format the sub_b cell for bedtools is empty; coerce numeric.
  df$sub_b <- suppressWarnings(as.numeric(df$sub_b))

  bt <- df %>% filter(tool == "bedtools") %>%
    select(num_files, num_threads, precision,
           bt_wall = mean_wall_time, bt_wall_err = std_wall_time,
           sort_wall = mean_sort_time)
  hm <- df %>% filter(grepl("^hammock_cpp_B", tool)) %>%
    mutate(sub_b = ifelse(is.na(sub_b), 1.0, sub_b)) %>%
    select(num_files, sub_b,
           hm_wall = mean_wall_time, hm_wall_err = std_wall_time,
           hm_sketch = mean_sketch_creation_time,
           hm_compare = mean_comparison_time)
  list(bt = arrange(bt, num_files),
       hm = arrange(hm, num_files, desc(sub_b)),
       sub_bs = sort(unique(hm$sub_b), decreasing = TRUE))
}

# Persistent crossover: largest i where bt[i] < hm[i] and bt[i+1] > hm[i+1];
# return the first N at which hammock wins persistently (the "+1" point).
# Operates on per-subB joined frame.
find_persistent_crossover <- function(d) {
  idx <- NA_integer_
  for (i in seq_len(nrow(d) - 1)) {
    if (d$bt_wall[i] < d$hm_wall[i] && d$bt_wall[i + 1] > d$hm_wall[i + 1]) {
      idx <- i + 1L
    }
  }
  idx
}

subb_label <- function(s) ifelse(abs(s - 1.0) < 1e-9, "subB=1.0",
                                 sprintf("subB=%g", s))

plot_files_sketch_compare_split <- function(files_csv_path, out_path) {
  loaded <- load_files_csv(files_csv_path)
  bt <- loaded$bt; hm <- loaded$hm; sub_bs <- loaded$sub_bs
  t <- bt$num_threads[1]; p <- bt$precision[1]

  # One hammock-total/sketching line per subB. Compare is sub_b-agnostic
  # in practice (N² register comparisons regardless of sketch contents),
  # so we plot just the subB=1.0 compare line (or the first one available).
  hm <- hm %>% mutate(subB_label = subb_label(sub_b))
  primary <- if (1.0 %in% sub_bs) 1.0 else sub_bs[1]
  hm_primary <- hm %>% filter(abs(sub_b - primary) < 1e-9)

  bedtools_series <- sprintf("bedtools wall (parallel-wrapped, t=%d)", t)
  sort_series     <- "bedtools pre-sort (excluded from bedtools wall)"
  compare_series  <- "hammock pairwise compare (≈ subB-agnostic)"

  sketch_series_label <- function(s) {
    sprintf("hammock sketching (t=%d, p=%d, %s)", t, p, subb_label(s))
  }

  # Build long frame for plotting.
  bt_rows   <- bt %>% transmute(num_files, y = bt_wall, yerr = bt_wall_err,
                                series = bedtools_series)
  sort_rows <- bt %>% transmute(num_files, y = sort_wall, yerr = NA_real_,
                                series = sort_series)
  sketch_rows <- hm %>% transmute(
    num_files, y = hm_sketch, yerr = hm_wall_err,
    series = sketch_series_label(sub_b)
  )
  compare_rows <- hm_primary %>% transmute(
    num_files, y = pmax(hm_compare, 1e-4), yerr = NA_real_,
    series = compare_series
  )
  series_df <- bind_rows(bt_rows, sketch_rows, compare_rows, sort_rows)

  sketch_levels <- vapply(sub_bs, sketch_series_label, character(1))
  series_levels <- c(bedtools_series, sketch_levels, compare_series, sort_series)
  series_df$series <- factor(series_df$series, levels = series_levels)

  # Colors: bedtools blue, sketch subBs from an orange→brown ramp, compare red, sort grey.
  if (length(sub_bs) == 1) {
    sketch_colors <- "#ff7f0e"
  } else {
    # Distinct shades per subB
    palette <- c("#ff7f0e", "#8c564b", "#9467bd", "#17becf")
    sketch_colors <- palette[seq_along(sub_bs)]
  }
  color_map <- setNames(c("#1f77b4", sketch_colors, "#d62728", "#7f7f7f"), series_levels)
  lty_map <- setNames(c("solid", rep("solid", length(sub_bs)), "solid", "dotted"), series_levels)
  shape_map <- setNames(c(16, rep(15, length(sub_bs)), 18, 4), series_levels)

  # Crossover annotation uses subB=1.0 (the default) only.
  d_cross <- inner_join(bt, hm_primary, by = "num_files") %>% arrange(num_files)
  cross_idx <- find_persistent_crossover(d_cross)

  g <- ggplot(series_df, aes(x = num_files, y = y,
                             color = series, linetype = series, shape = series))
  if (!is.na(cross_idx)) {
    g <- g + geom_vline(xintercept = d_cross$num_files[cross_idx],
                        color = "grey", linetype = "dashed", alpha = 0.4)
  }
  g <- g +
    geom_errorbar(data = subset(series_df, !is.na(yerr)),
                  aes(ymin = pmax(y - yerr, 1e-6), ymax = y + yerr),
                  width = 0.05, alpha = 0.7, show.legend = FALSE) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.4) +
    scale_color_manual(values = color_map, name = NULL) +
    scale_linetype_manual(values = lty_map, name = NULL) +
    scale_shape_manual(values = shape_map, name = NULL) +
    scale_x_continuous(
      trans = log2_trans(),
      breaks = bt$num_files,
      labels = bt$num_files,
      sec.axis = sec_axis(~ . ^ 2, name = "Pairwise compares (N²)",
                          breaks = bt$num_files ^ 2,
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
    ncross <- d_cross$num_files[cross_idx]
    g <- g + annotate("text",
                      x = ncross * 0.55, y = d_cross$bt_wall[cross_idx] * 4.5,
                      label = sprintf("persistent crossover\n(subB=1.0) at N=%d", ncross),
                      color = "grey40", size = 3.2, hjust = 1)
  }

  ggsave_png(out_path, g, width = 9.5, height = 6, dpi = 200)
  cat(sprintf("Wrote %s\n", out_path))
}

plot_cost_per_pair <- function(files_csv_path, out_path) {
  loaded <- load_files_csv(files_csv_path)
  bt <- loaded$bt; hm <- loaded$hm; sub_bs <- loaded$sub_bs
  t <- bt$num_threads[1]; p <- bt$precision[1]

  bt2 <- bt %>% mutate(N2 = num_files ^ 2, bt_pp = bt_wall / N2)
  hm2 <- hm %>% mutate(N2 = num_files ^ 2,
                       hm_pp = hm_wall / N2,
                       sk_pp = hm_sketch / N2,
                       cmp_pp = pmax(hm_compare / N2, 1e-9))

  bedtools_series <- sprintf("bedtools (parallel-wrapped, t=%d)", t)
  total_label  <- function(s) sprintf("hammock total (t=%d, p=%d, %s)", t, p, subb_label(s))
  sketch_label <- function(s) sprintf("hammock sketching (%s)", subb_label(s))
  compare_series <- "hammock pairwise compare"

  bt_rows <- bt2 %>% transmute(num_files, y = bt_pp, series = bedtools_series)
  total_rows  <- hm2 %>% transmute(num_files, y = hm_pp, series = total_label(sub_b))
  sketch_rows <- hm2 %>% transmute(num_files, y = sk_pp, series = sketch_label(sub_b))
  primary <- if (1.0 %in% sub_bs) 1.0 else sub_bs[1]
  cmp_rows <- hm2 %>% filter(abs(sub_b - primary) < 1e-9) %>%
    transmute(num_files, y = cmp_pp, series = compare_series)

  total_levels <- vapply(sub_bs, total_label, character(1))
  sketch_levels <- vapply(sub_bs, sketch_label, character(1))
  series_levels <- c(bedtools_series, total_levels, sketch_levels, compare_series)
  series_df <- bind_rows(bt_rows, total_rows, sketch_rows, cmp_rows)
  series_df$series <- factor(series_df$series, levels = series_levels)

  total_palette <- if (length(sub_bs) == 1) "#ff7f0e" else c("#ff7f0e", "#8c564b", "#9467bd", "#17becf")[seq_along(sub_bs)]
  color_map <- setNames(c("#1f77b4", total_palette, total_palette, "#d62728"), series_levels)
  lty_map <- setNames(c("solid", rep("solid", length(sub_bs)), rep("dashed", length(sub_bs)), "solid"),
                      series_levels)
  shape_map <- setNames(c(16, rep(15, length(sub_bs)), rep(17, length(sub_bs)), 18), series_levels)
  alpha_map <- setNames(c(1, rep(1, length(sub_bs)), rep(0.55, length(sub_bs)), 1), series_levels)

  g <- ggplot(series_df, aes(x = num_files, y = y,
                             color = series, linetype = series, shape = series, alpha = series)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.4) +
    scale_color_manual(values = color_map, name = NULL) +
    scale_linetype_manual(values = lty_map, name = NULL) +
    scale_shape_manual(values = shape_map, name = NULL) +
    scale_alpha_manual(values = alpha_map, name = NULL) +
    scale_x_continuous(trans = log2_trans(),
                       breaks = bt$num_files, labels = bt$num_files) +
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

  ggsave_png(out_path, g, width = 9.5, height = 6, dpi = 200)
  cat(sprintf("Wrote %s\n", out_path))
}

# ---------- jaccard pair plots ----------

load_pairs_csv <- function(path) {
  d <- read_csv(path, show_col_types = FALSE)
  if (!"sub_b" %in% names(d)) d$sub_b <- 1.0
  d %>%
    mutate(
      bedtools_jaccard = suppressWarnings(as.numeric(bedtools_jaccard)),
      hammock_jaccard  = suppressWarnings(as.numeric(hammock_jaccard)),
      precision        = as.integer(precision),
      sub_b            = suppressWarnings(as.numeric(sub_b))
    ) %>%
    mutate(sub_b = ifelse(is.na(sub_b), 1.0, sub_b)) %>%
    filter(!is.na(bedtools_jaccard), !is.na(hammock_jaccard))
}

plot_jaccard_scatter <- function(pairs_csv_path, out_path) {
  d <- load_pairs_csv(pairs_csv_path)
  if (nrow(d) == 0) { warning("No pairs found in ", pairs_csv_path); return(invisible(NULL)) }

  sub_bs <- sort(unique(d$sub_b), decreasing = TRUE)
  d$subb_facet <- factor(subb_label(d$sub_b),
                         levels = subb_label(sub_bs))

  counts <- d %>% group_by(precision, sub_b) %>% summarise(n = n(), .groups = "drop")
  d <- d %>% left_join(counts, by = c("precision", "sub_b")) %>%
    mutate(label = sprintf("p=%d  (n=%d)", precision, n))
  legend_order <- counts %>% filter(abs(sub_b - sub_bs[1]) < 1e-9) %>% arrange(precision) %>%
    mutate(label = sprintf("p=%d  (n=%d)", precision, n)) %>% pull(label)
  d$label <- factor(d$label, levels = legend_order)

  xpad <- max(0.005, diff(range(d$bedtools_jaccard)) * 0.15)
  ypad <- max(0.005, diff(range(d$hammock_jaccard)) * 0.15)
  xlim <- range(d$bedtools_jaccard) + c(-xpad, xpad)
  ylim <- range(d$hammock_jaccard)  + c(-ypad, ypad)
  lo <- min(xlim[1], ylim[1]); hi <- max(xlim[2], ylim[2])

  med_df <- d %>% group_by(subb_facet) %>%
    summarise(bt_med = median(bedtools_jaccard),
              hm_med = median(hammock_jaccard), .groups = "drop") %>%
    mutate(label = sprintf("median gap = %.3f\n(bedtools med = %.3f,\n hammock med = %.3f)",
                           hm_med - bt_med, bt_med, hm_med))

  g <- ggplot(d, aes(x = bedtools_jaccard, y = hammock_jaccard, color = label)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
    geom_point(size = 1.0, alpha = 0.35) +
    geom_segment(data = med_df,
                 aes(x = bt_med, xend = bt_med, y = bt_med, yend = hm_med),
                 arrow = arrow(length = unit(0.18, "cm")), color = "black",
                 linewidth = 0.7, inherit.aes = FALSE) +
    geom_text(data = med_df,
              aes(x = bt_med + (hi - bt_med) * 0.10,
                  y = (bt_med + hm_med) / 2, label = label),
              color = "black", hjust = 0, size = 3.0,
              inherit.aes = FALSE) +
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
  if (length(sub_bs) > 1) g <- g + facet_wrap(~ subb_facet)

  width <- if (length(sub_bs) > 1) 8.5 * length(sub_bs) else 8.5
  ggsave_png(out_path, g, width = width, height = 8, dpi = 200)
  cat(sprintf("Wrote %s\n", out_path))
}

plot_jaccard_delta <- function(pairs_csv_path, out_path) {
  d <- load_pairs_csv(pairs_csv_path)
  if (nrow(d) == 0) { warning("No pairs found in ", pairs_csv_path); return(invisible(NULL)) }
  d <- d %>% mutate(delta = hammock_jaccard - bedtools_jaccard)

  sub_bs <- sort(unique(d$sub_b), decreasing = TRUE)
  d$subb_facet <- factor(subb_label(d$sub_b),
                         levels = subb_label(sub_bs))

  counts <- d %>% group_by(precision, sub_b) %>% summarise(n = n(), .groups = "drop")
  d <- d %>% left_join(counts, by = c("precision", "sub_b")) %>%
    mutate(label = sprintf("p=%d  (n=%d)", precision, n))
  legend_order <- counts %>% filter(abs(sub_b - sub_bs[1]) < 1e-9) %>% arrange(precision) %>%
    mutate(label = sprintf("p=%d  (n=%d)", precision, n)) %>% pull(label)
  d$label <- factor(d$label, levels = legend_order)

  bins <- seq(min(d$bedtools_jaccard), max(d$bedtools_jaccard), length.out = 12)
  binned <- d %>%
    mutate(bin = cut(bedtools_jaccard, breaks = bins, include.lowest = TRUE)) %>%
    group_by(precision, sub_b, subb_facet, label, bin) %>%
    summarise(
      n = n(),
      x = mean(c(bins[as.integer(bin)], bins[as.integer(bin) + 1])),
      mean_delta = mean(delta),
      .groups = "drop"
    ) %>% filter(n >= 3) %>% arrange(precision, x)

  med_df <- d %>% group_by(subb_facet) %>%
    summarise(med_gap = median(delta), .groups = "drop")

  g <- ggplot(d, aes(x = bedtools_jaccard, y = delta, color = label)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8, alpha = 0.7) +
    geom_hline(data = med_df, aes(yintercept = med_gap),
               color = "red", linewidth = 0.6, linetype = "dotted",
               alpha = 0.7) +
    geom_point(size = 0.9, alpha = 0.25) +
    geom_line(data = binned, aes(x = x, y = mean_delta, color = label, group = label),
              linewidth = 1.0) +
    geom_text(data = med_df,
              aes(x = max(d$bedtools_jaccard), y = med_gap,
                  label = sprintf("  median gap = %.3f", med_gap)),
              color = "red", hjust = 1, vjust = -0.5, size = 3.0,
              inherit.aes = FALSE) +
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
  if (length(sub_bs) > 1) g <- g + facet_wrap(~ subb_facet)

  width <- if (length(sub_bs) > 1) 10 * length(sub_bs) else 10
  ggsave_png(out_path, g, width = width, height = 6, dpi = 200)
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
