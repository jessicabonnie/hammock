#!/usr/bin/env Rscript

# Figure 7 — sequence-mode parameter objective trade-off
#
# Show the numerical-agreement and biological-recovery objectives for the
# inclusion–exclusion Jaccard estimate at a selected precision (default p=24).
#
# Each point is one (k, w) configuration:
#   x = Pearson r with exact BEDTools Jaccard
#   y = adjusted Rand index (ARI) for the 10-tissue clustering
#   color = k-mer size k
#   size = minimizer window w
#
# Usage:
#   Rscript paper/parameter_objective_tradeoff/plot_parameter_objective_tradeoff_estimators.R [output.png] [precision]

required_packages <- c("dplyr", "readr", "ggplot2", "scales", "Cairo")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
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
if (length(script_arg) != 1) stop("Run this script with Rscript.", call. = FALSE)
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

summary_csv <- file.path(repo_root, "docs", "data", "mode_d_summary.csv")
experiment_summary_csv <- file.path(
  repo_root, "experiments", "maurano_dhs_validation", "results", "mode_d_summary.csv"
)
argv <- commandArgs(trailingOnly = TRUE)
out_png <- if (length(argv) >= 1) {
  normalizePath(argv[1], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "parameter_objective_tradeoff_estimators.png")
}
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

PRECISION <- if (length(argv) >= 2) as.integer(argv[2]) else 24L
if (is.na(PRECISION)) stop("Precision must be an integer.", call. = FALSE)
REFERENCE <- "bedtools"
FIG6_K <- 10
FIG6_W <- 30

COL_TEXT <- "#20262D"
COL_GRID <- "#D9DEE3"
COL_NUMERIC <- "#E69F00"
COL_BIO <- "#D81B60"
K_COLORS <- c("#B3DDF2", "#6DB7E3", "#2F8FD3", "#1764A0", "#08345E")
base_family <- "sans"

if (!file.exists(summary_csv)) stop("Input file not found: ", summary_csv, call. = FALSE)
raw <- read_csv(summary_csv, show_col_types = FALSE)
required_cols <- c("precision", "k", "w", "column", "reference", "pearson", "ari")
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop(basename(summary_csv), " lacks columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

SIM_COLUMN <- "jaccard_similarity_ie"
if (!SIM_COLUMN %in% unique(raw$column)) {
  if (!file.exists(experiment_summary_csv)) {
    stop("No rows available for ", SIM_COLUMN, call. = FALSE)
  }
  experiment_raw <- read_csv(experiment_summary_csv, show_col_types = FALSE)
  if (!SIM_COLUMN %in% unique(experiment_raw$column)) {
    stop("No rows available for ", SIM_COLUMN, call. = FALSE)
  }
  message("Supplementing staged summary from experiment results for ", SIM_COLUMN)
  raw <- bind_rows(raw, experiment_raw %>% filter(column == SIM_COLUMN))
}

sweep <- raw %>%
  filter(
    precision == PRECISION,
    reference == REFERENCE,
    column == SIM_COLUMN,
    k >= 8
  ) %>%
  mutate(
    pearson = as.numeric(pearson),
    ari = as.numeric(ari),
    k = as.integer(k),
    w = as.numeric(w)
  ) %>%
  filter(!is.na(pearson), !is.na(ari))

k_levels <- sort(unique(sweep$k))
if (length(k_levels) > length(K_COLORS)) {
  stop("More k values than palette steps; extend K_COLORS.", call. = FALSE)
}
k_colors <- setNames(K_COLORS[seq_along(k_levels)], as.character(k_levels))
sweep <- sweep %>% mutate(k_label = factor(as.character(k), levels = as.character(k_levels)))

best_numeric <- sweep %>%
  arrange(desc(pearson), desc(ari), k, w) %>%
  slice(1)

fig6_points <- sweep %>% filter(k == FIG6_K, w == FIG6_W)
if (nrow(fig6_points) != 1) {
  stop("Expected one k=10, w=30 inclusion–exclusion cell.", call. = FALSE)
}
if (fig6_points$ari < max(sweep$ari) - 1e-12) {
  warning("k=10, w=30 is not a maximum-ARI inclusion–exclusion cell.")
}

message(sprintf(
  "Numerical optimum: k=%d w=%g r=%.6f ARI=%.3f",
  best_numeric$k, best_numeric$w, best_numeric$pearson, best_numeric$ari
))
message(sprintf(
  "Biological optimum shown: k=%d w=%g r=%.6f ARI=%.3f",
  fig6_points$k, fig6_points$w, fig6_points$pearson, fig6_points$ari
))

p <- ggplot(sweep, aes(x = pearson, y = ari)) +
  geom_hline(yintercept = 0, linewidth = 0.55, color = "#707981") +
  geom_point(aes(color = k_label, size = w), alpha = 0.76, stroke = 0.25) +
  # The two highlighted configurations are the only non-dominated points in
  # this objective space; connect them to show the empirical Pareto trade-off.
  annotate(
    "segment",
    x = fig6_points$pearson, y = fig6_points$ari,
    xend = best_numeric$pearson, yend = best_numeric$ari,
    linewidth = 0.8, linetype = "dashed", color = "#6B747D",
    alpha = 0.72
  ) +
  # Numerical optimum: saturated-orange outline on the underlying point.
  geom_point(
    data = best_numeric,
    aes(x = pearson, y = ari, size = w),
    inherit.aes = FALSE,
    shape = 21, fill = NA, stroke = 1.35, color = COL_NUMERIC,
    show.legend = FALSE
  ) +
  # Manuscript biological optimum: magenta outline on the underlying point.
  geom_point(
    data = fig6_points,
    aes(x = pearson, y = ari, size = w),
    inherit.aes = FALSE,
    shape = 21, fill = NA, stroke = 1.35, color = COL_BIO,
    show.legend = FALSE
  ) +
  scale_color_manual(values = k_colors, name = "k-mer size (k)") +
  scale_fill_manual(values = k_colors, guide = "none") +
  # Explicit positive breaks: pretty_breaks([8,500]) includes 0, and
  # log10(0) = -Inf which NaNs the size legend (points themselves were fine).
  scale_size_continuous(
    trans = "log10", range = c(2.2, 5.2),
    breaks = c(10, 30, 100, 300), name = "Window size (w)"
  ) +
  scale_x_continuous(
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.06, 0.08))
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.06, 0.10))
  ) +
  labs(
    x = "Agreement with exact BEDTools Jaccard (Pearson r)",
    y = "Tissue recovery (adjusted Rand index)"
  ) +
  theme_classic(base_size = 11, base_family = base_family) +
  theme(
    axis.title = element_text(face = "bold", color = COL_TEXT),
    axis.text = element_text(color = COL_TEXT),
    axis.line = element_line(color = "#6B747D", linewidth = 0.35),
    axis.ticks = element_line(color = "#6B747D", linewidth = 0.35),
    panel.grid.major = element_line(color = COL_GRID, linewidth = 0.30),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 9.3, face = "bold", color = COL_TEXT),
    legend.text = element_text(size = 8.9, color = COL_TEXT),
    plot.margin = margin(12, 18, 12, 12)
  )

CairoPNG(
  filename = out_png,
  width = 8.2,
  height = 6.7,
  units = "in",
  res = 300,
  bg = "white"
)
print(p)
dev.off()
message("Wrote: ", out_png)
