#!/usr/bin/env Rscript

# Figure 7 — Sequence-mode parameter response
#
# Two paired line panels over the same ordered minimizer-window axis:
#   A  Pearson r against exact BEDTools Jaccard   (numerical agreement)
#   B  Adjusted Rand index against tissue labels  (biological resolution)
#
# Both panels are read from docs/data/mode_d_summary.csv at HyperLogLog
# precision p = 24, similarity column `jaccard_similarity` (minimizer only),
# reference `bedtools`. Reported values are computed from that file.
#
# Usage:
#   Rscript paper/parameter_response/plot_parameter_response.R [output.png]

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

summary_csv <- file.path(repo_root, "docs", "data", "mode_d_summary.csv")
argv <- commandArgs(trailingOnly = TRUE)
out_png <- if (length(argv) >= 1) {
  normalizePath(argv[1], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "parameter_response.png")
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
COL_MARK <- "#20262D"
base_family <- "sans"

# These values are intentionally centralized so a future paper-wide palette can
# replace them without touching the plotting logic.
K_COLORS <- c("#86B6EF", "#5598E7", "#2A78D6", "#184F95", "#0D366B")
K_SHAPES <- c(21, 24, 22, 23, 25)
K_LINETYPES <- c("solid", "22", "42", "13", "longdash")

theme_paper <- function(base_size = 10.5) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        face = "bold", size = rel(1.04), color = COL_TEXT,
        lineheight = 1.05, margin = margin(b = 7)
      ),
      axis.title = element_text(color = COL_TEXT),
      axis.text = element_text(color = COL_TEXT),
      axis.line = element_line(color = "#6B747D", linewidth = 0.35),
      axis.ticks = element_line(color = "#6B747D", linewidth = 0.35),
      panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = rel(0.9), color = COL_TEXT),
      legend.text = element_text(size = rel(0.85), color = COL_TEXT),
      legend.key.width = grid::unit(1.7, "lines"),
      plot.margin = margin(10, 12, 10, 10)
    )
}

raw <- read_csv(summary_csv, show_col_types = FALSE)
required_cols <- c("precision", "k", "w", "column", "reference", "pearson", "mae", "ari")
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
    reference == REFERENCE
  ) %>%
  mutate(across(c(pearson, mae, ari), as.numeric))

if (nrow(sweep) == 0) {
  stop("No rows for p = ", PRECISION, ", column ", SIM_COLUMN, call. = FALSE)
}

# k = 5 has too few usable windows to form a response curve.
dropped <- sweep %>% filter(k < 8)
sweep <- sweep %>% filter(k >= 8)
if (nrow(dropped) > 0) {
  message(sprintf(
    "Excluded %d row(s) with k < 8: %s",
    nrow(dropped),
    paste(sprintf("k=%d,w=%d", dropped$k, dropped$w), collapse = "; ")
  ))
}

w_levels <- sort(unique(sweep$w))
k_levels <- sort(unique(sweep$k))
if (length(k_levels) > length(K_COLORS)) {
  stop("More k values than palette steps; extend K_COLORS/K_SHAPES/K_LINETYPES.", call. = FALSE)
}

k_colors <- setNames(K_COLORS[seq_along(k_levels)], as.character(k_levels))
k_shapes <- setNames(K_SHAPES[seq_along(k_levels)], as.character(k_levels))
k_linetypes <- setNames(K_LINETYPES[seq_along(k_levels)], as.character(k_levels))

sweep <- sweep %>%
  mutate(
    x = match(w, w_levels),
    k_label = factor(as.character(k), levels = as.character(k_levels))
  )

guide_x <- match(FIG6_W, w_levels)
if (is.na(guide_x)) stop("w = ", FIG6_W, " is absent from the sweep.", call. = FALSE)

best_numeric <- sweep %>% filter(!is.na(pearson)) %>% slice_max(pearson, n = 1, with_ties = FALSE)
best_mae <- sweep %>% filter(!is.na(mae)) %>% slice_min(mae, n = 1, with_ties = FALSE)
fig6 <- sweep %>% filter(k == FIG6_K, w == FIG6_W)
if (nrow(fig6) != 1) {
  stop("Expected exactly one k = ", FIG6_K, ", w = ", FIG6_W, " row.", call. = FALSE)
}

message(sprintf(
  "Best Pearson r:  k = %d, w = %d -> r = %.4f, MAE = %.4f, ARI = %.3f",
  best_numeric$k, best_numeric$w, best_numeric$pearson, best_numeric$mae, best_numeric$ari
))
message(sprintf(
  "Lowest MAE:      k = %d, w = %d -> r = %.4f, MAE = %.4f, ARI = %.3f",
  best_mae$k, best_mae$w, best_mae$pearson, best_mae$mae, best_mae$ari
))
message(sprintf(
  "Figure 6 setting: k = %d, w = %d -> r = %.4f, MAE = %.4f, ARI = %.3f",
  fig6$k, fig6$w, fig6$pearson, fig6$mae, fig6$ari
))

annot_numeric <- sprintf(
  "Best numerical agreement\nk = %d, w = %d\nr = %.4f; ARI = %.3f",
  best_numeric$k, best_numeric$w, best_numeric$pearson, best_numeric$ari
)
annot_fig6 <- sprintf(
  "Figure 6 setting\nk = %d, w = %d\nr = %.4f; ARI = %.3f",
  fig6$k, fig6$w, fig6$pearson, fig6$ari
)

high_k <- sweep %>% filter(k >= 15, !is.na(ari))
high_k_ari <- unique(round(high_k$ari, 12))
high_k_label <- if (length(high_k_ari) == 1) {
  sprintf("k = 15, 20, and 25 overlap\nARI = %.3f", high_k_ari)
} else {
  NULL
}

series_layers <- function(df, yvar) {
  list(
    annotate(
      "segment",
      x = guide_x, xend = guide_x, y = -Inf, yend = Inf,
      color = "#C3CAD1", linewidth = 0.4, linetype = "22"
    ),
    geom_line(
      data = df,
      aes(
        x = x, y = .data[[yvar]], color = k_label,
        linetype = k_label, group = k_label
      ),
      linewidth = 0.85, na.rm = TRUE
    ),
    geom_point(
      data = df,
      aes(x = x, y = .data[[yvar]], fill = k_label, shape = k_label),
      color = "white", stroke = 0.3, size = 2.4, na.rm = TRUE
    ),
    scale_color_manual(values = k_colors, name = "k-mer size"),
    scale_fill_manual(values = k_colors, name = "k-mer size"),
    scale_shape_manual(values = k_shapes, name = "k-mer size"),
    scale_linetype_manual(values = k_linetypes, name = "k-mer size"),
    scale_x_continuous(
      breaks = seq_along(w_levels),
      labels = w_levels,
      limits = c(0.55, length(w_levels) + 0.45),
      expand = c(0, 0)
    )
  )
}

panel_a <- ggplot() +
  series_layers(sweep, "pearson") +
  geom_point(
    data = best_numeric, aes(x = x, y = pearson),
    shape = 21, size = 4.8, stroke = 1.0, fill = NA, color = COL_MARK
  ) +
  annotate(
    "segment",
    x = best_numeric$x - 0.05, xend = best_numeric$x - 0.30,
    y = best_numeric$pearson + 0.005, yend = best_numeric$pearson + 0.033,
    color = COL_MARK, linewidth = 0.35
  ) +
  annotate(
    "text",
    x = best_numeric$x - 0.36, y = best_numeric$pearson + 0.038,
    label = annot_numeric, hjust = 1, vjust = 0,
    size = 2.85, lineheight = 1.12, color = COL_TEXT, family = base_family
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.1),
    limits = c(0.55, 1.12),
    breaks = seq(0.6, 1.0, by = 0.1),
    expand = c(0, 0)
  ) +
  labs(
    title = "A  Agreement with exact interval similarity",
    x = "Minimizer window size w",
    y = "Pearson r vs BEDTools Jaccard"
  ) +
  theme_paper()

panel_b <- ggplot() +
  series_layers(sweep, "ari") +
  geom_point(
    data = fig6, aes(x = x, y = ari),
    shape = 21, size = 6.2, stroke = 1.2, fill = NA, color = COL_MARK
  ) +
  annotate(
    "segment",
    x = fig6$x + 0.10, xend = fig6$x + 0.34,
    y = fig6$ari + 0.02, yend = fig6$ari + 0.055,
    color = COL_MARK, linewidth = 0.35
  ) +
  annotate(
    "text",
    x = fig6$x + 0.40, y = fig6$ari + 0.058,
    label = annot_fig6, hjust = 0, vjust = 0,
    size = 2.85, lineheight = 1.12, color = COL_TEXT, family = base_family
  ) +
  {if (!is.null(high_k_label)) annotate(
    "text",
    x = length(w_levels) - 0.15, y = high_k_ari + 0.045,
    label = high_k_label, hjust = 1, vjust = 0,
    size = 2.75, lineheight = 1.1, color = "#46515C", family = base_family
  )} +
  scale_y_continuous(
    labels = label_number(accuracy = 0.1),
    limits = c(-0.09, 1.12),
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)
  ) +
  labs(
    title = "B  Recovery of tissue organization",
    x = "Minimizer window size w",
    y = "Adjusted Rand index vs tissue labels"
  ) +
  theme_paper()

figure <- panel_a + panel_b +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Numerical agreement and biological resolution favor different sequence-mode settings",
    theme = theme(
      plot.title = element_text(
        family = base_family, face = "bold", size = 14.5,
        color = COL_TEXT, margin = margin(b = 10)
      ),
      plot.margin = margin(12, 16, 10, 12)
    )
  ) &
  guides(
    color = guide_legend(nrow = 1, order = 1),
    fill = guide_legend(nrow = 1, order = 1),
    shape = guide_legend(nrow = 1, order = 1),
    linetype = guide_legend(nrow = 1, order = 1)
  ) &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.margin = margin(t = 2)
  )

CairoPNG(
  filename = out_png,
  width = 11.6,
  height = 5.6,
  units = "in",
  res = 300,
  bg = "white"
)
print(figure)
dev.off()

message("Wrote: ", out_png)
