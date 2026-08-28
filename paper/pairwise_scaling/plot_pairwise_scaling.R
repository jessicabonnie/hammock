# Figure 3 — Pairwise scaling of hammock versus BEDTools
# Creates a publication-ready two-panel PNG using CairoPNG.
#
# Main-text version: both panels report the +IE (jaccard_similarity_ie)
# column exclusively, in place of register-equality (jaccard_similarity).
# Simplified 2026-08-10 from a dual-metric design that plotted both hammock
# variants side by side; register-equality carries a chance-agreement floor
# and is not comparable to BEDTools on the same scale (CLAUDE.md divergence
# #2), so showing it next to a bedtools-calibrated column in the main text
# invites exactly that miscomparison. The full register-equality-vs-+IE
# picture (accuracy for BOTH now measured against exact BEDTools, so they
# ARE comparable) lives in plot_pairwise_scaling_supplement.R, for the
# supplement or advisor review, not this figure.
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

# p=18 (the CLI default), job 29652408, node c595, one exclusive allocation,
# 2026-08-09. Replaces cpp_vs_bedtools_t16_20260804_172242.csv, which was p=14 --
# a precision used nowhere else in the paper -- and whose bedtools leg ran under
# the three-process-per-pair bedtools.sh with its achieved parallelism
# unrecorded. Every bedtools row here carries mean_bedtools_parallel_eff; it
# reads 0.10 at N >= 64, i.e. the baseline converted ~1.6 of its 16 cores into
# throughput. Quote that alongside any speedup taken from this figure.
argv <- commandArgs(trailingOnly = TRUE)
out_png <- if (length(argv) >= 1) {
  normalizePath(argv[1], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "pairwise_scaling.png")
}
synthetic_csv <- if (length(argv) >= 2) {
  normalizePath(argv[2], mustWork = TRUE)
} else {
  file.path(repo_root, "docs", "data", "cpp_vs_bedtools_t16_p18.csv")
}
maurano_ie_summary_csv <- if (length(argv) >= 3) {
  normalizePath(argv[3], mustWork = TRUE)
} else {
  file.path(repo_root, "docs", "data", "maurano_subB_ie_summary.csv")
}
maurano_bedtools_csv <- if (length(argv) >= 4) {
  normalizePath(argv[4], mustWork = TRUE)
} else {
  file.path(repo_root, "docs", "data", "maurano_bedtools.csv")
}
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

for (path in c(synthetic_csv, maurano_bedtools_csv, maurano_ie_summary_csv)) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
}

COL_BEDTOOLS <- "#46515C"
COL_HAMMOCK <- "#007C83"
COL_COMPARE <- "#D28B35"
COL_COMPARE_IE <- "#9B4D9B"
# Distinct from COL_COMPARE_IE on purpose: 4 independent reviewers (2026-08-09)
# flagged that "hammock total (+IE)" reusing COL_COMPARE_IE made it collide,
# visually and in the legend, with the unrelated "hammock sketch comparison
# (+IE)" series, on top of already nearly-overlapping COL_HAMMOCK in the line
# itself. This color is not reused anywhere else in the panel.
COL_TOTAL_IE <- "#C2185B"
COL_GRID <- "#D9DEE3"
COL_TEXT <- "#20262D"
base_family <- "sans"

theme_paper <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        face = "bold", size = rel(1.06), color = COL_TEXT,
        lineheight = 1.05, margin = margin(b = 6)
      ),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 12, color = COL_TEXT),
      axis.text = element_text(size = 12, color = COL_TEXT),
      axis.line = element_line(color = "#6B747D", linewidth = 0.35),
      axis.ticks = element_line(color = "#6B747D", linewidth = 0.35),
      panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12, color = COL_TEXT),
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

# Uses the +IE (jaccard_similarity_ie) arm exclusively, in place of
# register-equality -- see the file header. wall_time/comparison_time here
# are the +IE arm's, i.e. the full metrics block (--metrics), not the
# reduced-column arm's timing (--register-equality, was --no-metrics) a
# pre-2026-08-10 reader of this file would expect.
if (!"hammock_ie_B" %in% synthetic_raw$tool) {
  stop(
    "No hammock_ie_B rows in ", basename(synthetic_csv), " -- Panel A now ",
    "requires the +IE arm (re-run the benchmark with --metrics-arm).",
    call. = FALSE
  )
}
synthetic_hm <- synthetic_raw %>%
  filter(tool == "hammock_ie_B") %>%
  transmute(
    num_files,
    wall_time = mean_wall_time,
    wall_sd = std_wall_time,
    comparison_time = mean_comparison_time,
    threads = num_threads,
    precision
  )

synthetic <- inner_join(
  synthetic_bt,
  synthetic_hm,
  by = c("num_files", "threads", "precision"),
  suffix = c("_bedtools", "_hammock")
) %>%
  arrange(num_files) %>%
  mutate(
    # Full cross product of two disjoint N-file sets, not within-set all-pairs.
    cross_pairs = num_files * num_files,
    speedup = wall_time_bedtools / wall_time_hammock
  )

if (nrow(synthetic) == 0) {
  stop("No matching +IE synthetic benchmark rows were found.", call. = FALSE)
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
  )
) %>%
  mutate(
    series = factor(
      series,
      levels = c("BEDTools total", "hammock total", "hammock sketch comparison")
    )
  )

n_breaks <- synthetic$num_files
pair_power_labels <- parse(text = paste0("2^", 2 * seq_along(n_breaks)))
largest <- synthetic %>% slice_max(num_files, n = 1, with_ties = FALSE)

speedup_label <- sprintf("%.1f× faster", largest$speedup)

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
    label = speedup_label,
    size = 3,
    linewidth = 0,
    fill = alpha("white", 0.9),
    color = COL_TEXT,
    hjust = 1,
    lineheight = 0.95
  ) +
  scale_color_manual(values = c(
    "BEDTools total" = COL_BEDTOOLS,
    "hammock total" = COL_HAMMOCK,
    "hammock sketch comparison" = COL_COMPARE
  ), drop = FALSE) +
  scale_linetype_manual(values = c(
    "BEDTools total" = "solid",
    "hammock total" = "solid",
    "hammock sketch comparison" = "22"
  ), drop = FALSE) +
  scale_shape_manual(values = c(
    "BEDTools total" = 16,
    "hammock total" = 17,
    "hammock sketch comparison" = 15
  ), drop = FALSE) +
  scale_x_continuous(
    trans = log2_trans(),
    breaks = n_breaks,
    labels = n_breaks,
    sec.axis = sec_axis(
      ~ ., breaks = n_breaks, labels = pair_power_labels,
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
    title = "A",
    x = "Number of BED files (N)",
    y = "Wall time (seconds, log scale)"
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE),
    linetype = guide_legend(nrow = 1, byrow = TRUE),
    shape = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme_paper() +
  theme(
    axis.text.x.top = element_text(size = 12, color = COL_TEXT),
    axis.title.x.top = element_text(size = 12, color = COL_TEXT),
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    legend.margin = margin(t = 5),
    plot.title = element_text(size = 16, lineheight = 1.05, margin = margin(b = 6)),
    plot.margin = margin(28, 8, 6, 6)
  )

# Panel B: Maurano real-data benchmark ----------------------------------------
maurano_bedtools <- read_csv(maurano_bedtools_csv, show_col_types = FALSE)
maurano_ie_summary <- read_csv(maurano_ie_summary_csv, show_col_types = FALSE)

required_ie <- c("subB", "wall_median", "mae_ie_vs_bedtools")
missing_ie <- setdiff(required_ie, names(maurano_ie_summary))
if (length(missing_ie) > 0) {
  stop(
    "Maurano +IE summary lacks columns: ",
    paste(missing_ie, collapse = ", "),
    call. = FALSE
  )
}
if (!all(c("rep", "run_id", "wall_time") %in% names(maurano_bedtools))) {
  stop("Maurano BEDTools input lacks rep, run_id, or wall_time.", call. = FALSE)
}

bt_runs <- maurano_bedtools %>% distinct(rep, run_id, wall_time)
bt_wall <- median(bt_runs$wall_time, na.rm = TRUE)

required_subb <- c(1, 0.1, 0.01)
if (!all(required_subb %in% maurano_ie_summary$subB)) {
  stop("Expected +IE summary to include subB = 1, 0.1, and 0.01.", call. = FALSE)
}

condition_of <- function(subb) case_when(
  subb == 1 ~ "no\nsubsampling",
  TRUE ~ paste0(format(subb, scientific = FALSE, trim = TRUE), "\nsubsample")
)

subb_levels <- maurano_ie_summary %>%
  arrange(desc(subB)) %>%
  pull(subB) %>%
  condition_of()

# One bar per condition: BEDTools, then hammock at the +IE arm (full metrics
# block) for each subB level -- in place of register-equality, same
# rationale as Panel A. Wall time here is the +IE arm's, costlier than the
# register-equality (--register-equality, was --no-metrics) timing this
# panel plotted before 2026-08-10; see plot_pairwise_scaling_supplement.R
# for that comparison.
bars <- bind_rows(
  tibble(
    condition = "BEDTools",
    tool = "BEDTools",
    wall = bt_wall,
    mae = NA_real_
  ),
  maurano_ie_summary %>% transmute(
    condition = condition_of(subB),
    tool = "hammock (+IE)",
    wall = wall_median,
    mae = mae_ie_vs_bedtools
  )
) %>%
  mutate(
    condition = factor(
      condition,
      levels = c("BEDTools", subb_levels)
    ),
    tool = factor(tool, levels = c("BEDTools", "hammock (+IE)")),
    speedup = bt_wall / wall,
    # Never hardcode "faster" -- word it from the sign: the corrected
    # BEDTools baseline makes the unsubsampled bar a genuine slowdown, and
    # "0.90x faster" would read as a speedup to anyone skimming.
    ratio_txt = ifelse(speedup >= 1,
                       sprintf("%.2f× faster", speedup),
                       sprintf("%.2f× slower", 1 / speedup)),
    # Keep the speed comparison on its own line so the larger annotation text
    # remains centered over the rightmost bar without clipping the panel edge.
    label = case_when(
      tool == "BEDTools" ~ sprintf("%.1f s\nspeed reference", wall),
      TRUE ~ sprintf(
        "%.1f s\n(%s)\nMAE %s",
        wall, ratio_txt, formatC(mae, format = "e", digits = 1)
      )
    )
  )

expanded_panel_b <- nrow(maurano_ie_summary) > 3
panel_b <- ggplot(bars, aes(x = condition, y = wall, fill = tool)) +
  geom_col(width = 0.68) +
  geom_text(
    aes(label = label),
    vjust = -0.22,
    size = if (expanded_panel_b) 2.7 else 4.1,
    lineheight = 0.95,
    color = COL_TEXT
  ) +
  scale_fill_manual(values = c(
    "BEDTools" = COL_BEDTOOLS,
    "hammock (+IE)" = COL_HAMMOCK
  )) +
  scale_y_continuous(
    labels = label_number(accuracy = 1),
    expand = expansion(mult = c(0, 0.32))
  ) +
  labs(
    title = "B",
    x = NULL,
    # Deliberately not "per pairwise comparison": each bar is the median
    # TOTAL wall time to sketch all 20 files and run all 400 pairwise
    # comparisons in one benchmark invocation, not divided by pair count.
    y = "Wall time, 20-file corpus (s)"
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(
      size = if (expanded_panel_b) 8.5 else 12,
      lineheight = 0.95
    ),
    legend.position = "none",
    plot.title = element_text(size = 16, lineheight = 1.05, margin = margin(b = 6)),
    plot.margin = margin(28, 8, 6, 6)
  )

figure <- panel_a + panel_b +
  plot_layout(widths = c(1.35, 1))

CairoPNG(
  filename = out_png,
  width = if (expanded_panel_b) 17.5 else 14.2,
  height = 6.5,
  units = "in",
  res = 300,
  bg = "white"
)
print(figure)
dev.off()

message("Wrote: ", out_png)
