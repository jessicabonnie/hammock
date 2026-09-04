#!/usr/bin/env Rscript

# Figure 7 — sequence-mode parameter objective trade-off
#
# Show the numerical-agreement and biological-recovery objectives for the
# production rehash-selector64 sequence-HLL estimator. Each cell is summarized
# by the median across the eight prespecified, independent rehash seeds.
#
# Each point is one (k, w) configuration:
#   x = Pearson r with exact BEDTools coordinate Jaccard
#   y = adjusted Rand index (ARI) for the 10-tissue clustering
#   color = minimizer window w
#   shape = k-mer size k
#
# Usage:
#   Rscript paper/parameter_objective_tradeoff/plot_parameter_objective_tradeoff_estimators.R [output.png] [precision] [classes] [x_min]

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

summary_csv <- file.path(
  repo_root, "experiments", "maurano_dhs_rehash_validation", "results",
  "summaries", "rehash_run_scores.csv"
)
argv <- commandArgs(trailingOnly = TRUE)
out_png <- if (length(argv) >= 1) {
  normalizePath(argv[1], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "parameter_objective_tradeoff_estimators_10class.png")
}
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

PRECISION <- if (length(argv) >= 2) as.integer(argv[2]) else 18L
if (is.na(PRECISION)) stop("Precision must be an integer.", call. = FALSE)
CLASSES <- if (length(argv) >= 3) as.integer(argv[3]) else 10L
if (!CLASSES %in% c(8L, 10L)) stop("Classes must be 8 or 10.", call. = FALSE)
X_MIN <- if (length(argv) >= 4) as.numeric(argv[4]) else NA_real_
if (!is.na(X_MIN) && (X_MIN < 0 || X_MIN >= 1)) {
  stop("x_min must be at least 0 and less than 1.", call. = FALSE)
}
FIG6_K <- 10
FIG6_W <- 30
HASH_SCHEME <- "rehash-selector64"
PRESPECIFIED_SEEDS <- c(1L, 2L, 3L, 17L, 42L, 43L, 99L, 31337L)
EXPECTED_CELLS <- 37L

COL_TEXT <- "#20262D"
COL_GRID <- "#D9DEE3"
COL_ORANGE <- "#FF8C00"
COL_MAGENTA <- "#D81B60"
K_COLORS <- c("#B3DDF2", "#6DB7E3", "#2F8FD3", "#1764A0", "#08345E")
K_SHAPES <- c(12, 7, 4, 3, 8)
base_family <- "sans"

if (!file.exists(summary_csv)) stop("Input file not found: ", summary_csv, call. = FALSE)
raw <- read_csv(summary_csv, show_col_types = FALSE)
required_cols <- c(
  "hash_scheme", "precision", "seed", "k", "w", "bedtools_pearson",
  "ari_10class", "ari_8class"
)
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop(basename(summary_csv), " lacks columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

sweep_runs <- raw %>%
  filter(
    hash_scheme == HASH_SCHEME,
    precision == PRECISION,
    seed %in% PRESPECIFIED_SEEDS,
    k >= 8
  ) %>%
  mutate(
    bedtools_pearson = as.numeric(bedtools_pearson),
    ari_10class = as.numeric(ari_10class),
    ari_8class = as.numeric(ari_8class),
    seed = as.integer(seed),
    k = as.integer(k),
    w = as.numeric(w)
  ) %>%
  filter(
    !is.na(bedtools_pearson), !is.na(ari_10class), !is.na(ari_8class)
  )

duplicate_runs <- sweep_runs %>%
  count(k, w, seed, name = "copies") %>%
  filter(copies != 1L)
if (nrow(duplicate_runs) > 0L) {
  stop("Duplicate rehash runs found for one or more (k, w, seed) keys.", call. = FALSE)
}

cell_counts <- sweep_runs %>%
  group_by(k, w) %>%
  summarise(seeds_observed = n_distinct(seed), .groups = "drop")
if (nrow(cell_counts) != EXPECTED_CELLS ||
    any(cell_counts$seeds_observed != length(PRESPECIFIED_SEEDS))) {
  stop(
    sprintf(
      "Expected %d cells with %d prespecified rehash seeds each at p=%d.",
      EXPECTED_CELLS, length(PRESPECIFIED_SEEDS), PRECISION
    ),
    call. = FALSE
  )
}

sweep <- sweep_runs %>%
  group_by(k, w) %>%
  summarise(
    pearson = median(bedtools_pearson),
    ari_10class = median(ari_10class),
    ari_8class = median(ari_8class),
    .groups = "drop"
  ) %>%
  mutate(ari = ari_10class)

if (CLASSES == 8L) {
  sweep <- sweep %>%
    select(-ari) %>%
    rename(ari = ari_8class)
}

k_levels <- sort(unique(sweep$k))
if (length(k_levels) > length(K_COLORS) || length(k_levels) > length(K_SHAPES)) {
  stop("More k values than available colors or point shapes.", call. = FALSE)
}
k_shapes <- setNames(K_SHAPES[seq_along(k_levels)], as.character(k_levels))
sweep <- sweep %>%
  mutate(k_label = factor(as.character(k), levels = as.character(k_levels)))
window_limits <- range(sweep$w)

best_numeric <- sweep %>%
  arrange(desc(pearson), desc(ari), k, w) %>%
  slice(1)

fig6_points <- sweep %>% filter(k == FIG6_K, w == FIG6_W)
if (nrow(fig6_points) != 1) {
  stop("Expected one k=10, w=30 inclusion–exclusion cell.", call. = FALSE)
}
if (CLASSES == 10L && fig6_points$ari < max(sweep$ari) - 1e-12) {
  warning("k=10, w=30 is not a maximum-ARI inclusion–exclusion cell.")
}
bio_points <- if (CLASSES == 8L) {
  sweep %>% filter(abs(ari - max(ari)) <= 1e-12)
} else fig6_points
plot_sweep <- if (is.na(X_MIN)) sweep else
  sweep %>% filter(abs(ari - max(ari)) <= 1e-12)
highlight_numeric <- if (is.na(X_MIN)) best_numeric else
  best_numeric %>% filter(abs(ari - max(sweep$ari)) <= 1e-12)

message(sprintf(
  "Numerical optimum: k=%d w=%g r=%.6f ARI=%.3f",
  best_numeric$k, best_numeric$w, best_numeric$pearson, best_numeric$ari
))
message(sprintf("Biological optimum contains %d cell(s); maximum ARI=%.6f",
                nrow(bio_points), max(sweep$ari)))

p <- ggplot(plot_sweep, aes(x = pearson, y = ari)) +
  geom_hline(yintercept = 0, linewidth = 0.55, color = "#707981")

if (CLASSES == 10L) {
  p <- p +
    geom_point(
      aes(color = w, shape = k_label),
      size = 3.8, alpha = 0.82, stroke = 0.4
    ) +
    # Color the numerical and biological optima using their native k shapes.
    geom_point(
      data = best_numeric,
      aes(x = pearson, y = ari, shape = k_label),
      inherit.aes = FALSE,
      size = 3.8, stroke = 0.8, color = COL_MAGENTA,
      show.legend = FALSE
    ) +
    geom_point(
      data = fig6_points,
      aes(x = pearson, y = ari, shape = k_label),
      inherit.aes = FALSE,
      size = 3.8, stroke = 0.8, color = COL_ORANGE,
      show.legend = FALSE
    ) +
    scale_color_gradientn(
      colors = K_COLORS, trans = "log10",
      limits = window_limits, breaks = c(10, 30, 100, 300),
      name = "Window size (w)"
    ) +
    scale_shape_manual(
      values = k_shapes, limits = names(k_shapes), drop = FALSE,
      name = "k-mer size (k)"
    ) +
    guides(
      color = guide_colorbar(order = 1),
      shape = guide_legend(
        order = 2, override.aes = list(color = "#333333", alpha = 1)
      )
    )
} else {
  p <- p +
    geom_point(
      aes(color = w, shape = k_label),
      size = 3.8, alpha = 0.82, stroke = 0.4
    ) +
    # Color every tied biological optimum bright orange.
    geom_point(
      data = if (is.na(X_MIN)) bio_points else bio_points[0, ],
      aes(x = pearson, y = ari, shape = k_label),
      inherit.aes = FALSE,
      size = 3.8, stroke = 0.8, color = COL_ORANGE,
      show.legend = FALSE
    ) +
    # Color the numerical optimum's own k-mer-size shape magenta.
    geom_point(
      data = highlight_numeric,
      aes(x = pearson, y = ari, shape = k_label),
      inherit.aes = FALSE,
      size = 3.8, stroke = 0.8, color = COL_MAGENTA,
      show.legend = FALSE
    ) +
    scale_color_gradientn(
      colors = K_COLORS, trans = "log10",
      limits = window_limits, breaks = c(10, 30, 100, 300),
      name = "Window size (w)"
    ) +
    scale_shape_manual(
      values = k_shapes, limits = names(k_shapes), drop = FALSE,
      name = "k-mer size (k)"
    ) +
    guides(
      color = guide_colorbar(order = 1),
      shape = guide_legend(
        order = 2, override.aes = list(color = "#333333", alpha = 1)
      )
    )
}

p <- p +
  scale_x_continuous(
    labels = label_number(accuracy = 0.01),
    breaks = if (is.na(X_MIN)) waiver() else seq(X_MIN, 1.00, length.out = 6),
    expand = if (is.na(X_MIN)) expansion(mult = c(0.06, 0.08)) else
      expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.06, 0.10))
  ) +
  labs(
    x = if (CLASSES == 8L) "Median agreement with exact BEDTools Jaccard (Pearson r)" else
      "Agreement with exact BEDTools Jaccard (Pearson r)",
    y = if (CLASSES == 8L) "Median eight-class tissue recovery (adjusted Rand index)" else
      "Tissue recovery (adjusted Rand index)"
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

if (!is.na(X_MIN)) {
  p <- p + coord_cartesian(xlim = c(X_MIN, 1.00), clip = "off")
}

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
