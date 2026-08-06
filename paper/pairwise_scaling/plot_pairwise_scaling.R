# Figure 3 — Pairwise scaling of hammock versus BEDTools
# Creates a publication-ready two-panel PNG using CairoPNG.
#
# Panel A was rebuilt in Aug 2026; see docs/figure3-panel-a-rebuild.md for the
# evidence. Three things there are load-bearing and easy to undo by accident:
#
#   1. The secondary axis counts N^2. The harness builds two DISJOINT sets of N
#      files and both tools run the full cross product, so N=512 is 262,144
#      pairs, not 130,816. N(N-1)/2 is Panel B's denominator, not Panel A's.
#   2. There is no pmax() floor on comparison_time, and its absence is checked.
#      hammock-cpp reported integer milliseconds before v0.7.0, so a zero here
#      means the CSV predates the microsecond timers -- fail, do not floor.
#   3. Breaks are pinned to decades. Left to itself log_breaks() picks breaks
#      100x apart over this range and the axis silently reads two decades per
#      gridline.

required_packages <- c("dplyr", "readr", "ggplot2", "scales", "patchwork", "Cairo")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(Cairo)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
  stop("Could not determine the script path. Run with Rscript.", call. = FALSE)
}
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

synthetic_csv <- file.path(
  repo_root, "docs", "data", "cpp_vs_bedtools_t16_20260804_172242.csv"
)
maurano_summary_csv <- file.path(
  repo_root, "docs", "data", "maurano_subB_summary.csv"
)
maurano_bedtools_csv <- file.path(
  repo_root, "docs", "data", "maurano_bedtools.csv"
)

argv <- commandArgs(trailingOnly = TRUE)
out_png <- if (length(argv) >= 1) {
  normalizePath(argv[1], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "pairwise_scaling.png")
}
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

for (path in c(synthetic_csv, maurano_summary_csv, maurano_bedtools_csv)) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
}

COL_BEDTOOLS <- "#46515C"
COL_HAMMOCK <- "#007C83"
COL_COMPARE <- "#D28B35"
COL_COMPARE_IE <- "#9B4D9B"
COL_GRID <- "#D9DEE3"
COL_TEXT <- "#20262D"
base_family <- "sans"

theme_paper <- function(base_size = 10.5) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        face = "bold", size = rel(1.06), color = COL_TEXT,
        lineheight = 1.05, margin = margin(b = 6)
      ),
      plot.subtitle = element_blank(),
      axis.title = element_text(color = COL_TEXT),
      axis.text = element_text(color = COL_TEXT),
      axis.line = element_line(color = "#6B747D", linewidth = 0.35),
      axis.ticks = element_line(color = "#6B747D", linewidth = 0.35),
      panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.8), color = COL_TEXT),
      legend.key.width = grid::unit(1.1, "lines"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

# Panel A: synthetic scaling ----------------------------------------------------
synthetic_raw <- read_csv(synthetic_csv, show_col_types = FALSE)

synthetic_bt <- synthetic_raw %>%
  filter(tool == "bedtools") %>%
  transmute(
    num_files,
    wall_time = mean_wall_time,
    wall_sd = std_wall_time,
    threads = num_threads,
    precision
  )

synthetic_hm <- synthetic_raw %>%
  mutate(sub_b = suppressWarnings(as.numeric(sub_b))) %>%
  filter(grepl("^hammock_cpp_B", tool), abs(sub_b - 1) < 1e-9) %>%
  transmute(
    num_files,
    wall_time = mean_wall_time,
    wall_sd = std_wall_time,
    comparison_time = mean_comparison_time,
    threads = num_threads,
    precision
  )

# A different row, not a different column of synthetic_hm, so it cannot be
# transmuted out of it. left_join, not inner_join: against a CSV predating the
# --metrics-arm option this must yield NA and drop the series, not empty
# Panel A and trip the stop() below with a misleading message.
synthetic_ie <- synthetic_raw %>%
  filter(tool == "hammock_ie_B") %>%
  transmute(
    num_files,
    threads = num_threads,
    precision,
    comparison_time_ie = mean_comparison_time
  )

synthetic <- inner_join(
  synthetic_bt,
  synthetic_hm,
  by = c("num_files", "threads", "precision"),
  suffix = c("_bedtools", "_hammock")
) %>%
  left_join(synthetic_ie, by = c("num_files", "threads", "precision")) %>%
  arrange(num_files) %>%
  mutate(
    # Full cross product of two disjoint N-file sets, not within-set all-pairs.
    cross_pairs = num_files * num_files,
    speedup = wall_time_bedtools / wall_time_hammock
  )

if (nrow(synthetic) == 0) {
  stop("No matching unsubsampled synthetic benchmark rows were found.", call. = FALSE)
}

# Without the pmax() floor a zero lands at -Inf on the log axis and the point
# vanishes silently. A zero here means the CSV predates the microsecond timers
# (v0.7.0), where the comparison phase read 0.000000 for N <= 16.
if (any(synthetic$comparison_time <= 0, na.rm = TRUE)) {
  stop(
    "comparison_time <= 0: this CSV predates the microsecond timers in ",
    "hammock-cpp 0.7.0. Re-run the benchmark rather than restoring a floor.",
    call. = FALSE
  )
}
if (!"hammock_ie_B" %in% synthetic_raw$tool) {
  message(
    "note: no hammock_ie_B rows in ", basename(synthetic_csv),
    " -- the +IE series will be absent. Re-run with --metrics-arm."
  )
}

synthetic_long <- bind_rows(
  synthetic %>% transmute(
    num_files, value = wall_time_bedtools, error = wall_sd_bedtools,
    series = "BEDTools total"
  ),
  synthetic %>% transmute(
    num_files, value = wall_time_hammock, error = wall_sd_hammock,
    series = "hammock total"
  ),
  synthetic %>% transmute(
    num_files, value = comparison_time, error = NA_real_,
    series = "hammock sketch comparison"
  ),
  synthetic %>%
    filter(!is.na(comparison_time_ie)) %>%
    transmute(
      num_files, value = comparison_time_ie, error = NA_real_,
      series = "hammock sketch comparison (+IE)"
    )
) %>%
  mutate(
    # A series absent from levels becomes NA and ggplot drops it with only a
    # warning, so the IE level is declared whether or not any row carries it.
    series = factor(
      series,
      levels = c("BEDTools total", "hammock total",
                 "hammock sketch comparison", "hammock sketch comparison (+IE)")
    )
  )

n_breaks <- synthetic$num_files
pair_labels <- format(
  synthetic$cross_pairs,
  scientific = FALSE,
  big.mark = ",",
  trim = TRUE
)
largest <- synthetic %>% slice_max(num_files, n = 1, with_ties = FALSE)

panel_a <- ggplot(
  synthetic_long,
  aes(x = num_files, y = value, color = series, linetype = series, shape = series)
) +
  geom_errorbar(
    data = filter(synthetic_long, !is.na(error)),
    aes(ymin = pmax(value - error, 1e-4), ymax = value + error),
    width = 0.05,
    linewidth = 0.35,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.25, stroke = 0.4) +
  annotate(
    "segment",
    x = largest$num_files,
    xend = largest$num_files,
    y = largest$wall_time_hammock * 1.18,
    yend = largest$wall_time_bedtools / 1.18,
    linewidth = 0.45,
    color = "#69737D",
    arrow = grid::arrow(ends = "both", length = grid::unit(0.08, "in"))
  ) +
  annotate(
    "label",
    x = largest$num_files / 1.18,
    y = sqrt(largest$wall_time_hammock * largest$wall_time_bedtools),
    label = sprintf("%.1f× faster", largest$speedup),
    size = 3,
    linewidth = 0,
    fill = alpha("white", 0.9),
    color = COL_TEXT,
    hjust = 1
  ) +
  scale_color_manual(values = c(
    "BEDTools total" = COL_BEDTOOLS,
    "hammock total" = COL_HAMMOCK,
    "hammock sketch comparison" = COL_COMPARE,
    "hammock sketch comparison (+IE)" = COL_COMPARE_IE
  ), drop = FALSE) +
  scale_linetype_manual(values = c(
    "BEDTools total" = "solid",
    "hammock total" = "solid",
    "hammock sketch comparison" = "22",
    "hammock sketch comparison (+IE)" = "42"
  ), drop = FALSE) +
  scale_shape_manual(values = c(
    "BEDTools total" = 16,
    "hammock total" = 17,
    "hammock sketch comparison" = 15,
    "hammock sketch comparison (+IE)" = 18
  ), drop = FALSE) +
  scale_x_continuous(
    trans = log2_trans(),
    breaks = n_breaks,
    labels = n_breaks,
    sec.axis = sec_axis(
      ~ ., breaks = n_breaks, labels = pair_labels,
      name = "File pairs compared, N²"
    )
  ) +
  scale_y_log10(
    # Breaks are pinned to decades. The data spans 6.6 decades, over which
    # scale_y_log10()'s default log_breaks() picks breaks 100x apart
    # (1e-4, 0.01, 1, 100) -- so consecutive gridlines are TWO decades and the
    # axis appears to jump 10 -> 1,000. On a log axis that reads as a decade
    # step unless you check the numbers.
    breaks = 10^(-4:3),
    minor_breaks = NULL,
    # Not label_number(accuracy = 0.1) either: the comparison series reaches
    # 2e-4 s, and a fixed 0.1 accuracy renders every sub-decisecond break as
    # "0.0". Three significant figures in fixed notation labels 0.0001 and
    # 1,000 correctly with one rule.
    labels = function(x) formatC(x, format = "fg", digits = 3, big.mark = ","),
    expand = expansion(mult = c(0.06, 0.15))
  ) +
  labs(
    title = "A  Sketch reuse increases the advantage as collections grow",
    x = "Number of BED files (N)",
    y = "Wall time (seconds, log scale)"
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_paper() +
  theme(
    axis.text.x.top = element_text(size = 7.7, color = "#59636D"),
    axis.title.x.top = element_text(size = 8.5, color = "#59636D"),
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.margin = margin(b = 2)
  )

# Panel B: Maurano real-data benchmark ----------------------------------------
maurano_summary <- read_csv(maurano_summary_csv, show_col_types = FALSE)
maurano_bedtools <- read_csv(maurano_bedtools_csv, show_col_types = FALSE)

required_summary <- c("method", "subB", "wall_median", "mae")
missing_summary <- setdiff(required_summary, names(maurano_summary))
if (length(missing_summary) > 0) {
  stop(
    "Maurano summary lacks columns: ",
    paste(missing_summary, collapse = ", "),
    call. = FALSE
  )
}
if (!all(c("rep", "run_id", "wall_time") %in% names(maurano_bedtools))) {
  stop("Maurano BEDTools input lacks rep, run_id, or wall_time.", call. = FALSE)
}

bt_runs <- maurano_bedtools %>% distinct(rep, run_id, wall_time)
bt_wall <- median(bt_runs$wall_time, na.rm = TRUE)

mixed_stride <- maurano_summary %>%
  filter(method == "mixed-stride", subB %in% c(1, 0.1, 0.01)) %>%
  arrange(match(subB, c(1, 0.1, 0.01)))

if (nrow(mixed_stride) != 3) {
  stop("Expected mixed-stride rows for subB = 1, 0.1, and 0.01.", call. = FALSE)
}

bars <- bind_rows(
  tibble(
    condition = "BEDTools",
    tool = "BEDTools",
    wall = bt_wall,
    mae = NA_real_
  ),
  mixed_stride %>% transmute(
    condition = case_when(
      subB == 1 ~ "hammock\nno subsampling",
      subB == 0.1 ~ "hammock\nsubB = 0.1",
      subB == 0.01 ~ "hammock\nsubB = 0.01"
    ),
    tool = "hammock",
    wall = wall_median,
    mae
  )
) %>%
  mutate(
    condition = factor(condition, levels = condition),
    tool = factor(tool, levels = c("BEDTools", "hammock")),
    speedup = bt_wall / wall,
    # "mean |ΔJ|", not "ΔJ": mae is mean(abs(j - j_truth)) over the 950
    # pair-by-replicate comparisons, so it carries magnitude but no direction.
    # And subB = 1.0 IS the baseline it is measured against, so its zero is true
    # by construction -- printing a bare "0" beside "1.16x faster" reads as
    # "agrees exactly with BEDTools", which is false by ~0.16 here (the
    # register-equality chance floor, see CLAUDE.md divergence #2).
    label = case_when(
      tool == "BEDTools" ~ sprintf("%.1f s\nspeed reference", wall),
      mae == 0 ~ sprintf("%.1f s\n%.2f× faster\nΔJ baseline", wall, speedup),
      TRUE ~ sprintf(
        "%.1f s\n%.2f× faster\nmean |ΔJ| = %s",
        wall, speedup, formatC(mae, format = "e", digits = 1)
      )
    )
  )

panel_b <- ggplot(bars, aes(x = condition, y = wall, fill = tool)) +
  geom_col(width = 0.68) +
  geom_text(
    aes(label = label),
    vjust = -0.28,
    size = 3.05,
    lineheight = 0.95,
    color = COL_TEXT
  ) +
  scale_fill_manual(values = c(
    "BEDTools" = COL_BEDTOOLS,
    "hammock" = COL_HAMMOCK
  )) +
  scale_y_continuous(
    labels = label_number(accuracy = 1),
    expand = expansion(mult = c(0, 0.27))
  ) +
  labs(
    title = "B  Subsampling further reduces runtime",
    x = NULL,
    y = "Wall time (seconds)"
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(size = 8.8, lineheight = 0.95),
    legend.position = "none",
    plot.title = element_text(size = 10.5, lineheight = 1.05)
  )

figure <- panel_a + panel_b +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(
    title = "Hammock expands feasible all-pairs comparison as interval collections grow",
    subtitle = paste0(
      "Reusable sketches reduce repeated full-file processing; optional subsampling ",
      "further lowers sketch-construction time with little change in estimated similarity."
    ),
    theme = theme(
      plot.title = element_text(
        family = base_family, face = "bold", size = 15.5, color = COL_TEXT,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        family = base_family, size = 10.3, color = "#56616C",
        margin = margin(b = 8)
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

CairoPNG(
  filename = out_png,
  width = 14.2,
  height = 6.5,
  units = "in",
  res = 300,
  bg = "white"
)
print(figure)
dev.off()

message("Wrote: ", out_png)
