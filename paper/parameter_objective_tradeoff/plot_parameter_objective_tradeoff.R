#!/usr/bin/env Rscript

# Figure 7 — Sequence-mode parameter objective trade-off
#
# Each point is one (k, w) configuration. The x-axis measures numerical
# agreement with exact BEDTools Jaccard (Pearson r), while the y-axis measures
# recovery of annotated tissue structure (adjusted Rand index). This makes the
# central result explicit: parameter settings that best reproduce coordinate-
# based similarity are not necessarily those that best preserve biological
# organization.
#
# Data source: docs/data/mode_d_summary.csv, supplemented from the experiment
# summary when the staged paper copy is missing a current metric.
# Filters: p = 24, BEDTools reference; one figure per Jaccard estimator.
#
# Usage:
#   Rscript paper/parameter_objective_tradeoff/plot_parameter_objective_tradeoff.R \
#     [register_equality.png] [inclusion_exclusion.png]

required_packages <- c("dplyr", "readr", "ggplot2", "scales", "Cairo")
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
  library(Cairo)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
  stop("Could not determine the script path. Run with Rscript.", call. = FALSE)
}
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

summary_csv <- file.path(repo_root, "docs", "data", "mode_d_summary.csv")
experiment_summary_csv <- file.path(
  repo_root, "experiments", "maurano_dhs_validation",
  "results", "mode_d_summary.csv"
)
argv <- commandArgs(trailingOnly = TRUE)
out_png <- if (length(argv) >= 1) {
  normalizePath(argv[1], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "parameter_objective_tradeoff.png")
}
out_ie_png <- if (length(argv) >= 2) {
  normalizePath(argv[2], mustWork = FALSE)
} else {
  sub("\\.png$", "_ie.png", out_png)
}
for (path in c(out_png, out_ie_png)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
if (!file.exists(summary_csv)) stop("Input file not found: ", summary_csv, call. = FALSE)

PRECISION <- 24
SIM_COLUMNS <- c("jaccard_similarity", "jaccard_similarity_ie")
REFERENCE <- "bedtools"
FIG6_K <- 10
FIG6_W <- 30

COL_TEXT <- "#20262D"
COL_GRID <- "#D9DEE3"
COL_HIGHLIGHT <- "#D94F2B"
base_family <- "sans"

# Sequential blue palette with larger lightness steps between discrete k values.
# Extend if the sweep gains additional k values.
K_COLORS <- c("#B3DDF2", "#6DB7E3", "#2F8FD3", "#1764A0", "#08345E")

raw <- read_csv(summary_csv, show_col_types = FALSE)
required_cols <- c("precision", "k", "w", "column", "reference", "pearson", "ari")
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop(
    basename(summary_csv), " lacks columns: ",
    paste(missing_cols, collapse = ", "), call. = FALSE
  )
}

missing_metrics <- setdiff(SIM_COLUMNS, unique(raw$column))
if (length(missing_metrics) > 0) {
  if (!file.exists(experiment_summary_csv)) {
    stop(
      basename(summary_csv), " lacks metrics: ",
      paste(missing_metrics, collapse = ", "),
      "; experiment summary not found at ", experiment_summary_csv,
      call. = FALSE
    )
  }
  experiment_raw <- read_csv(experiment_summary_csv, show_col_types = FALSE)
  still_missing <- setdiff(missing_metrics, unique(experiment_raw$column))
  if (length(still_missing) > 0) {
    stop(
      "No summary rows available for metrics: ",
      paste(still_missing, collapse = ", "), call. = FALSE
    )
  }
  message(
    "Supplementing ", basename(summary_csv), " from experiment summary for: ",
    paste(missing_metrics, collapse = ", ")
  )
  raw <- bind_rows(
    raw,
    experiment_raw %>% filter(column %in% missing_metrics)
  )
}

render_tradeoff <- function(sim_column, output_path) {
sweep <- raw %>%
  filter(
    precision == PRECISION,
    column == sim_column,
    reference == REFERENCE,
    k >= 8
  ) %>%
  mutate(
    pearson = as.numeric(pearson),
    ari = as.numeric(ari)
  ) %>%
  filter(!is.na(pearson), !is.na(ari))

if (nrow(sweep) == 0) {
  stop("No usable rows for p = ", PRECISION, ", column ", sim_column, call. = FALSE)
}

k_levels <- sort(unique(sweep$k))
if (length(k_levels) > length(K_COLORS)) {
  stop("More k values than palette steps; extend K_COLORS.", call. = FALSE)
}
k_colors <- setNames(K_COLORS[seq_along(k_levels)], as.character(k_levels))

sweep <- sweep %>%
  mutate(k_label = factor(as.character(k), levels = as.character(k_levels)))

best_numeric <- sweep %>%
  slice_max(pearson, n = 1, with_ties = FALSE)

best_biological <- sweep %>%
  arrange(desc(ari), desc(pearson)) %>%
  slice(1)

fig6 <- sweep %>% filter(k == FIG6_K, w == FIG6_W)
if (nrow(fig6) != 1) {
  stop("Expected exactly one k = ", FIG6_K, ", w = ", FIG6_W, " row.", call. = FALSE)
}

message(sprintf(
  "Best numerical agreement: k = %d, w = %d -> r = %.4f, ARI = %.3f",
  best_numeric$k, best_numeric$w, best_numeric$pearson, best_numeric$ari
))
message(sprintf(
  "Best biological resolution: k = %d, w = %d -> r = %.4f, ARI = %.3f",
  best_biological$k, best_biological$w, best_biological$pearson, best_biological$ari
))
message(sprintf(
  "Figure 6 setting: k = %d, w = %d -> r = %.4f, ARI = %.3f",
  fig6$k, fig6$w, fig6$pearson, fig6$ari
))

# If Figure 6 is itself the biological optimum, do not double-label it.
fig6_is_bio_optimum <- (
  fig6$k == best_biological$k && fig6$w == best_biological$w
)

label_numeric <- sprintf(
  "Best numerical agreement\nk = %d, w = %d",
  best_numeric$k, best_numeric$w
)

label_bio <- sprintf(
  "Best tissue group recovery\nk = %d, w = %d",
  best_biological$k, best_biological$w
)

# Fixed positions keep the callouts separate. The numerical label sits
# below-left of its point so its short leader avoids the high-agreement row.
numeric_label_x <- 0.855
numeric_label_y <- 0.585
bio_label_x <- 0.785
bio_label_y <- 0.940

x_axis_label <- if (identical(sim_column, "jaccard_similarity_ie")) {
  "Agreement with exact BEDTools Jaccard (Pearson r)\nInclusion–exclusion estimator"
} else {
  "Agreement with exact BEDTools Jaccard (Pearson r)"
}

p <- ggplot(sweep, aes(x = pearson, y = ari)) +
  geom_hline(
    yintercept = 0,
    color = "#6B747D",
    linewidth = 0.65
  ) +
  annotate(
    "text",
    x = 0.995, y = 0.015,
    label = "ARI = 0 (chance agreement)",
    hjust = 1, vjust = 0,
    size = 2.8,
    color = "#59636D",
    family = base_family
  ) +
  geom_point(
    aes(color = k_label, size = w),
    alpha = 0.78,
    stroke = 0.25
  ) +
  geom_point(
    data = best_numeric,
    shape = 21, size = 6.1, stroke = 1.15,
    fill = NA, color = COL_HIGHLIGHT,
    inherit.aes = FALSE,
    aes(x = pearson, y = ari)
  ) +
  geom_point(
    data = best_biological,
    shape = 21, size = 6.1, stroke = 1.15,
    fill = NA, color = COL_HIGHLIGHT,
    inherit.aes = FALSE,
    aes(x = pearson, y = ari)
  ) +
  {if (!fig6_is_bio_optimum) geom_point(
    data = fig6,
    shape = 21, size = 6.1, stroke = 1.15,
    fill = NA, color = COL_HIGHLIGHT,
    inherit.aes = FALSE,
    aes(x = pearson, y = ari)
  )} +
  annotate(
    "segment",
    x = numeric_label_x + 0.115, y = numeric_label_y,
    xend = best_numeric$pearson - 0.006, yend = best_numeric$ari,
    linewidth = 0.55, color = COL_HIGHLIGHT
  ) +
  annotate(
    "segment",
    x = bio_label_x + 0.115, y = bio_label_y,
    xend = best_biological$pearson - 0.006, yend = best_biological$ari,
    linewidth = 0.55, color = COL_HIGHLIGHT
  ) +
  annotate(
    "label",
    x = numeric_label_x,
    y = numeric_label_y,
    label = label_numeric,
    hjust = 0, vjust = 0.5,
    size = 3.05, lineheight = 1.12,
    label.size = 0.25,
    color = COL_TEXT,
    fill = "white",
    family = base_family
  ) +
  annotate(
    "label",
    x = bio_label_x,
    y = bio_label_y,
    label = label_bio,
    hjust = 0, vjust = 0.5,
    size = 3.05, lineheight = 1.12,
    label.size = 0.25,
    color = COL_TEXT,
    fill = "white",
    family = base_family
  ) +
  {if (!fig6_is_bio_optimum) annotate(
    "label",
    x = fig6$pearson,
    y = fig6$ari,
    label = sprintf(
      "Figure 6 setting\nk = %d, w = %d\nr = %.4f; ARI = %.3f",
      fig6$k, fig6$w, fig6$pearson, fig6$ari
    ),
    hjust = -0.05, vjust = -0.15,
    size = 3.0, lineheight = 1.08,
    label.size = 0.25,
    color = COL_TEXT,
    fill = alpha("white", 0.92),
    family = base_family
  )} +
  scale_color_manual(values = k_colors, name = "k-mer size (k)") +
  # Explicit positive breaks: pretty_breaks([8,500]) includes 0, and
  # log10(0) = -Inf which NaNs the size legend (points themselves were fine).
  scale_size_continuous(
    trans = "log10",
    range = c(2.2, 5.4),
    breaks = c(10, 30, 100, 300),
    name = "Window size (w)"
  ) +
  scale_x_continuous(
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.05, 0.11))
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.06, 0.13))
  ) +
  labs(
    title = "Numerical agreement and tissue recovery favor different settings",
    subtitle = "Each point is one (k, w) configuration; color shows k and size shows w",
    x = x_axis_label,
    y = "Tissue recovery (adjusted Rand index)"
  ) +
  theme_classic(base_size = 11, base_family = base_family) +
  theme(
    plot.title = element_text(
      face = "bold", size = 14.2, color = COL_TEXT,
      lineheight = 1.05, margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 10.2, color = "#4C5661", margin = margin(b = 10)
    ),
    axis.title = element_text(color = COL_TEXT),
    axis.text = element_text(color = COL_TEXT),
    axis.line = element_line(color = "#6B747D", linewidth = 0.35),
    axis.ticks = element_line(color = "#6B747D", linewidth = 0.35),
    panel.grid.major = element_line(color = COL_GRID, linewidth = 0.32),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 9.4, color = COL_TEXT),
    legend.text = element_text(size = 9.0, color = COL_TEXT),
    plot.margin = margin(12, 18, 12, 12)
  )

CairoPNG(
  filename = output_path,
  width = 8.6,
  height = 6.3,
  units = "in",
  res = 300,
  bg = "white"
)
print(p)
dev.off()

message("Wrote ", sim_column, ": ", output_path)
}

render_tradeoff(SIM_COLUMNS[[1]], out_png)
render_tradeoff(SIM_COLUMNS[[2]], out_ie_png)
