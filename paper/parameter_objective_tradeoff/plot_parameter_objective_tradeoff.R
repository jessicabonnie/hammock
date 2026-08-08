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
# Data source: docs/data/mode_d_summary.csv
# Default filter: p = 24, minimizer-only jaccard_similarity, BEDTools reference.
#
# Usage:
#   Rscript paper/parameter_objective_tradeoff/plot_parameter_objective_tradeoff.R [output.png]

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
argv <- commandArgs(trailingOnly = TRUE)
out_png <- if (length(argv) >= 1) {
  normalizePath(argv[1], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "parameter_objective_tradeoff.png")
}
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
if (!file.exists(summary_csv)) stop("Input file not found: ", summary_csv, call. = FALSE)

PRECISION <- 24
SIM_COLUMN <- "jaccard_similarity"
REFERENCE <- "bedtools"
FIG6_K <- 10
FIG6_W <- 30

COL_TEXT <- "#20262D"
COL_GRID <- "#D9DEE3"
COL_FRONTIER <- "#5D6670"
COL_HIGHLIGHT <- "#20262D"
base_family <- "sans"

# Centralized palette for k. Extend if the sweep gains additional k values.
K_COLORS <- c("#86B6EF", "#5598E7", "#2A78D6", "#184F95", "#0D366B")

raw <- read_csv(summary_csv, show_col_types = FALSE)
required_cols <- c("precision", "k", "w", "column", "reference", "pearson", "ari")
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop(
    basename(summary_csv), " lacks columns: ",
    paste(missing_cols, collapse = ", "), call. = FALSE
  )
}

sweep <- raw %>%
  filter(
    precision == PRECISION,
    column == SIM_COLUMN,
    reference == REFERENCE,
    k >= 8
  ) %>%
  mutate(
    pearson = as.numeric(pearson),
    ari = as.numeric(ari)
  ) %>%
  filter(!is.na(pearson), !is.na(ari))

if (nrow(sweep) == 0) {
  stop("No usable rows for p = ", PRECISION, ", column ", SIM_COLUMN, call. = FALSE)
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

# A point is Pareto-optimal when no other configuration is at least as good on
# both objectives and strictly better on one of them.
is_pareto <- function(i, df) {
  dominated <- (
    df$pearson >= df$pearson[i] &
    df$ari >= df$ari[i] &
    (df$pearson > df$pearson[i] | df$ari > df$ari[i])
  )
  !any(dominated, na.rm = TRUE)
}
pareto <- sweep[vapply(seq_len(nrow(sweep)), is_pareto, logical(1), df = sweep), ] %>%
  arrange(pearson, ari)

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
  "Best numerical agreement\nk = %d, w = %d\nr = %.4f; ARI = %.3f",
  best_numeric$k, best_numeric$w,
  best_numeric$pearson, best_numeric$ari
)

label_bio <- if (fig6_is_bio_optimum) {
  sprintf(
    "Best tissue recovery / Figure 6\nk = %d, w = %d\nr = %.4f; ARI = %.3f",
    best_biological$k, best_biological$w,
    best_biological$pearson, best_biological$ari
  )
} else {
  sprintf(
    "Best tissue recovery\nk = %d, w = %d\nr = %.4f; ARI = %.3f",
    best_biological$k, best_biological$w,
    best_biological$pearson, best_biological$ari
  )
}

p <- ggplot(sweep, aes(x = pearson, y = ari)) +
  geom_path(
    data = pareto,
    aes(x = pearson, y = ari, group = 1),
    inherit.aes = FALSE,
    color = COL_FRONTIER,
    linewidth = 0.65,
    linetype = "22"
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
    "label",
    x = best_numeric$pearson,
    y = best_numeric$ari,
    label = label_numeric,
    hjust = 1.05, vjust = -0.25,
    size = 3.05, lineheight = 1.08,
    label.size = 0.25,
    color = COL_TEXT,
    fill = alpha("white", 0.92),
    family = base_family
  ) +
  annotate(
    "label",
    x = best_biological$pearson,
    y = best_biological$ari,
    label = label_bio,
    hjust = -0.05, vjust = 1.15,
    size = 3.05, lineheight = 1.08,
    label.size = 0.25,
    color = COL_TEXT,
    fill = alpha("white", 0.92),
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
  scale_size_continuous(
    trans = "log10",
    range = c(2.2, 5.4),
    breaks = pretty_breaks(n = 4),
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
    title = "Numerical agreement and biological resolution favor different sequence-mode settings",
    subtitle = "Each point is one (k, w) configuration; dashed line connects Pareto-optimal settings",
    x = "Agreement with exact BEDTools Jaccard (Pearson r)",
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
  filename = out_png,
  width = 8.6,
  height = 6.3,
  units = "in",
  res = 300,
  bg = "white"
)
print(p)
dev.off()

message("Wrote: ", out_png)
