#!/usr/bin/env Rscript

# Figure 3 — Pairwise scaling of hammock versus BEDTools
#
# Creates a publication-ready two-panel PNG without requiring X11 or Cairo.
#
# Usage:
#   Rscript paper/pairwise_scaling/plot_pairwise_scaling.R
#   Rscript paper/pairwise_scaling/plot_pairwise_scaling.R path/to/output.png

required_packages <- c(
  "dplyr", "readr", "ggplot2", "scales", "patchwork", "ragg"
)
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
})

# Resolve repository-relative paths.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
  stop("Could not determine the script path. Run with Rscript.", call. = FALSE)
}
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

synthetic_csv <- file.path(
  repo_root, "docs", "data", "cpp_vs_bedtools_t16_20260512_160412.csv"
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

# Shared publication palette.
COL_BEDTOOLS <- "#46515C"
COL_HAMMOCK <- "#007C83"
COL_COMPARE <- "#D28B35"
COL_HAMMOCK_LIGHT <- "#6BB7B5"
COL_GRID <- "#D9DEE3"
COL_TEXT <- "#20262D"
base_family <- "sans"

theme_paper <- function(base_size = 10.5) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        face = "bold", size = rel(1.08), color = COL_TEXT,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = rel(0.88), color = "#56616C", lineheight = 1.05,
        margin = margin(b = 8)
      ),
      axis.title = element_text(color = COL_TEXT),
      axis.text = element_text(color = COL_TEXT),
      axis.line = element_line(color = "#6B747D", linewidth = 0.35),
      axis.ticks = element_line(color = "#6B747D", linewidth = 0.35),
      panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.box.just = "left",
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.82), color = COL_TEXT),
      legend.key.width = grid::unit(1.2, "lines"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

# -----------------------------------------------------------------------------
# Panel A: synthetic scaling
# -----------------------------------------------------------------------------
synthetic_raw <- read_csv(synthetic_csv, show_col_types = FALSE)

required_synthetic <- c(
  "num_files", "num_threads", "precision", "sub_b", "tool",
  "mean_wall_time", "std_wall_time", "mean_comparison_time"
)
missing_synthetic <- setdiff(required_synthetic, names(synthetic_raw))
if (length(missing_synthetic) > 0) {
  stop(
    "Synthetic input lacks columns: ",
    paste(missing_synthetic, collapse = ", "), call. = FALSE
  )
}

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

synthetic <- inner_join(
  synthetic_bt,
  synthetic_hm,
  by = c("num_files", "threads", "precision"),
  suffix = c("_bedtools", "_hammock")
) %>%
  arrange(num_files) %>%
  mutate(
    unique_pairs = num_files * (num_files - 1) / 2,
    speedup = wall_time_bedtools / wall_time_hammock
  )

if (nrow(synthetic) == 0) {
  stop("No matching unsubsampled synthetic benchmark rows were found.", call. = FALSE)
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
    num_files, value = pmax(comparison_time, 1e-4), error = NA_real_,
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
pair_labels <- format(
  synthetic$unique_pairs,
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
    "hammock sketch comparison" = COL_COMPARE
  )) +
  scale_linetype_manual(values = c(
    "BEDTools total" = "solid",
    "hammock total" = "solid",
    "hammock sketch comparison" = "22"
  )) +
  scale_shape_manual(values = c(
    "BEDTools total" = 16,
    "hammock total" = 17,
    "hammock sketch comparison" = 15
  )) +
  scale_x_continuous(
    trans = log2_trans(),
    breaks = n_breaks,
    labels = n_breaks,
    sec.axis = sec_axis(
      ~ ., breaks = n_breaks, labels = pair_labels,
      name = "Unique file pairs, N(N−1)/2"
    )
  ) +
  scale_y_log10(
    labels = label_number(accuracy = 0.1, big.mark = ","),
    expand = expansion(mult = c(0.06, 0.15))
  ) +
  labs(
    title = "A  Sketch reuse increases the advantage as collections grow",
    subtitle = sprintf(
      "Synthetic BED collections; 10,000 intervals per file; p = %d; %d threads; three runs per configuration",
      synthetic$precision[1], synthetic$threads[1]
    ),
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
    axis.text.x.top = element_text(size = 7.7, color = "#59636D"),
    axis.title.x.top = element_text(size = 8.5, color = "#59636D"),
    legend.margin = margin(b = 2)
  )

# -----------------------------------------------------------------------------
# Panel B: Maurano real-data benchmark
# -----------------------------------------------------------------------------
maurano_summary <- read_csv(maurano_summary_csv, show_col_types = FALSE)
maurano_bedtools <- read_csv(maurano_bedtools_csv, show_col_types = FALSE)

required_summary <- c("method", "subB", "wall_median", "mae")
missing_summary <- setdiff(required_summary, names(maurano_summary))
if (length(missing_summary) > 0) {
  stop(
    "Maurano summary lacks columns: ",
    paste(missing_summary, collapse = ", "), call. = FALSE
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
    speedup = bt_wall / wall,
    label = if_else(
      tool == "BEDTools",
      sprintf("%.1f s\nreference", wall),
      sprintf(
        "%.1f s\n%.2f× faster\nΔJ = %s",
        wall,
        speedup,
        if_else(mae == 0, "0", formatC(mae, format = "e", digits = 1))
      )
    )
  )

panel_b <- ggplot(bars, aes(x = condition, y = wall, fill = condition)) +
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
    "hammock\nno subsampling" = COL_HAMMOCK,
    "hammock\nsubB = 0.1" = COL_HAMMOCK_LIGHT,
    "hammock\nsubB = 0.01" = COL_COMPARE
  )) +
  scale_y_continuous(
    labels = label_number(accuracy = 1),
    expand = expansion(mult = c(0, 0.27))
  ) +
  labs(
    title = "B  Real-data speed can be increased with controlled approximation",
    subtitle = paste0(
      "Maurano fetal-tissue DHS; 20 BED files; 190 unique pairs; ",
      "interval mode; p = 18; 8 threads"
    ),
    x = NULL,
    y = "Wall time (seconds)"
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(size = 8.8, lineheight = 0.95),
    legend.position = "none"
  )

# Assemble and save through ragg, which is fully headless and does not use X11.
figure <- panel_a + panel_b +
  plot_layout(widths = c(1.35, 1), guides = "collect") +
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
  ) &
  theme(legend.position = "top")

ragg::agg_png(
  filename = out_png,
  width = 14.2,
  height = 6.6,
  units = "in",
  res = 300,
  background = "white"
)
print(figure)
dev.off()

message("Wrote: ", out_png)
