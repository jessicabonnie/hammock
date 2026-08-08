#!/usr/bin/env Rscript

# Figure 7 — paired sequence-mode estimator objective trade-off
#
# Purpose:
#   Compare the two HLL similarity estimators on the same p=24 sequence-mode
#   parameter sweep. The figure is intentionally faceted rather than rendered
#   as two independent plots so the reader can compare the numerical-calibration
#   shift directly while seeing that the useful biological optimum is stable.
#
# Each point is one (k, w) configuration:
#   x = Pearson r with exact BEDTools Jaccard
#   y = adjusted Rand index (ARI) for the 10-tissue clustering
#   color = k-mer size k
#   size = minimizer window w
#
# Panel A: register equality (`jaccard_similarity`)
# Panel B: inclusion–exclusion Jaccard (`jaccard_similarity_ie`)
#
# The script leaves plot_parameter_objective_tradeoff.R unchanged.
#
# Usage:
#   Rscript paper/parameter_objective_tradeoff/plot_parameter_objective_tradeoff_estimators.R [output.png]

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
  file.path(
    repo_root, "paper", "figures",
    "parameter_objective_tradeoff_estimators.png"
  )
}
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

PRECISION <- 24
REFERENCE <- "bedtools"
SIM_COLUMNS <- c("jaccard_similarity", "jaccard_similarity_ie")
FIG6_K <- 10
FIG6_W <- 30
ARI_TOL <- 1e-12

COL_TEXT <- "#20262D"
COL_GRID <- "#D9DEE3"
COL_NUMERIC <- "#D94F2B"
COL_BIO <- "#20262D"
COL_DISCORD <- "#7A828A"
base_family <- "sans"

K_COLORS <- c("#B3DDF2", "#6DB7E3", "#2F8FD3", "#1764A0", "#08345E")

if (!file.exists(summary_csv)) {
  stop("Input file not found: ", summary_csv, call. = FALSE)
}

raw <- read_csv(summary_csv, show_col_types = FALSE)
required_cols <- c("precision", "k", "w", "column", "reference", "pearson", "ari")
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop(
    basename(summary_csv), " lacks columns: ",
    paste(missing_cols, collapse = ", "), call. = FALSE
  )
}

# The staged paper copy may lag the experiment summary for a newly added metric.
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
  raw <- bind_rows(
    raw,
    experiment_raw %>% filter(column %in% missing_metrics)
  )
}

panel_levels <- c(
  "A  Register equality",
  "B  Inclusion–exclusion Jaccard"
)
panel_lookup <- c(
  jaccard_similarity = panel_levels[[1]],
  jaccard_similarity_ie = panel_levels[[2]]
)

sweep <- raw %>%
  filter(
    precision == PRECISION,
    reference == REFERENCE,
    column %in% SIM_COLUMNS,
    k >= 8
  ) %>%
  mutate(
    pearson = as.numeric(pearson),
    ari = as.numeric(ari),
    k = as.integer(k),
    w = as.numeric(w),
    estimator = factor(unname(panel_lookup[column]), levels = panel_levels)
  ) %>%
  filter(!is.na(pearson), !is.na(ari))

if (nrow(sweep) == 0) {
  stop("No usable p=24 estimator rows found.", call. = FALSE)
}

counts <- sweep %>% count(column, name = "n")
if (nrow(counts) != 2 || length(unique(counts$n)) != 1) {
  stop(
    "The two estimator panels do not contain the same number of p=24 cells.",
    call. = FALSE
  )
}

k_levels <- sort(unique(sweep$k))
if (length(k_levels) > length(K_COLORS)) {
  stop("More k values than palette steps; extend K_COLORS.", call. = FALSE)
}
k_colors <- setNames(K_COLORS[seq_along(k_levels)], as.character(k_levels))

sweep <- sweep %>%
  mutate(k_label = factor(as.character(k), levels = as.character(k_levels)))

# Match configurations across estimators to quantify biological agreement at p=24.
ari_compare <- sweep %>%
  select(column, k, w, ari) %>%
  tidyr::pivot_wider(names_from = column, values_from = ari) %>%
  filter(!is.na(jaccard_similarity), !is.na(jaccard_similarity_ie)) %>%
  mutate(
    delta_ari = jaccard_similarity_ie - jaccard_similarity,
    differs = abs(delta_ari) > ARI_TOL
  )

# Avoid adding tidyr as a required package: recreate the comparison if tidyr is
# not installed using base merge. This branch is normally skipped when tidyr is
# available through the plotting environment.
if (!requireNamespace("tidyr", quietly = TRUE)) {
  re <- sweep %>%
    filter(column == "jaccard_similarity") %>%
    select(k, w, ari_re = ari)
  ie <- sweep %>%
    filter(column == "jaccard_similarity_ie") %>%
    select(k, w, ari_ie = ari)
  ari_compare <- merge(re, ie, by = c("k", "w")) %>%
    mutate(
      delta_ari = ari_ie - ari_re,
      differs = abs(delta_ari) > ARI_TOL
    )
}

n_cells <- nrow(ari_compare)
n_diff <- sum(ari_compare$differs)
max_delta <- if (n_cells > 0) max(abs(ari_compare$delta_ari)) else NA_real_

discordant_keys <- ari_compare %>% filter(differs) %>% select(k, w)
discordant_points <- if (nrow(discordant_keys) > 0) {
  inner_join(sweep, discordant_keys, by = c("k", "w"))
} else {
  sweep[0, ]
}

best_numeric <- sweep %>%
  group_by(column, estimator) %>%
  arrange(desc(pearson), desc(ari), k, w, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

max_ari <- sweep %>%
  group_by(column, estimator) %>%
  summarise(max_ari = max(ari), .groups = "drop")

best_biological <- sweep %>%
  inner_join(max_ari, by = c("column", "estimator")) %>%
  filter(abs(ari - max_ari) <= ARI_TOL) %>%
  arrange(column, k, w)

fig6_points <- sweep %>%
  filter(k == FIG6_K, w == FIG6_W)
if (nrow(fig6_points) != 2) {
  stop("Expected the k=10, w=30 cell under both estimators.", call. = FALSE)
}

# Verify the manuscript cell is on the biological optimum for both estimators.
fig6_optimum_check <- fig6_points %>%
  inner_join(max_ari, by = c("column", "estimator")) %>%
  mutate(is_optimum = abs(ari - max_ari) <= ARI_TOL)
if (!all(fig6_optimum_check$is_optimum)) {
  warning("k=10, w=30 is not the maximum-ARI cell under both estimators.")
}

bio_ari <- mean(fig6_points$ari)

# Labels are placed at stable panel-relative positions so the matched facets
# read as one figure rather than two unrelated plots.
panel_ranges <- sweep %>%
  group_by(column, estimator) %>%
  summarise(
    xmin = min(pearson), xmax = max(pearson),
    ymin = min(ari), ymax = max(ari),
    .groups = "drop"
  ) %>%
  mutate(
    xspan = pmax(xmax - xmin, 0.05),
    yspan = pmax(ymax - ymin, 0.10)
  )

numeric_labels <- best_numeric %>%
  inner_join(panel_ranges, by = c("column", "estimator")) %>%
  mutate(
    label_x = pmax(xmin + 0.04 * xspan, pearson - 0.34 * xspan),
    label_y = pmax(ymin + 0.08 * yspan, ari - 0.13 * yspan),
    label = sprintf("Best numerical agreement\nk=%d, w=%g\nr=%.4f", k, w, pearson)
  )

bio_labels <- fig6_points %>%
  inner_join(panel_ranges, by = c("column", "estimator")) %>%
  mutate(
    label_x = pmin(xmax - 0.02 * xspan, pearson + 0.06 * xspan),
    label_y = pmin(ymax - 0.02 * yspan, ari + 0.045 * yspan),
    label = sprintf("Biological optimum\nk=%d, w=%g\nARI=%.3f", k, w, ari)
  )

agreement_note <- sprintf(
  "At p=24, ARI differs between estimators at %d of %d parameter cells; both peak at k=%d, w=%d (ARI=%.3f).",
  n_diff, n_cells, FIG6_K, FIG6_W, bio_ari
)
if (is.finite(max_delta) && n_diff > 0) {
  agreement_note <- paste0(
    agreement_note,
    sprintf(" Maximum |ΔARI| at p=24 = %.3f.", max_delta)
  )
}

message(agreement_note)
for (i in seq_len(nrow(best_numeric))) {
  message(sprintf(
    "%s numerical optimum: k=%d w=%g r=%.6f ARI=%.3f",
    as.character(best_numeric$estimator[i]),
    best_numeric$k[i], best_numeric$w[i],
    best_numeric$pearson[i], best_numeric$ari[i]
  ))
}

p <- ggplot(sweep, aes(x = pearson, y = ari)) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.55,
    color = "#707981"
  ) +
  geom_point(
    aes(color = k_label, size = w),
    alpha = 0.76,
    stroke = 0.25
  ) +
  # The rare p=24 estimator-discordant configuration(s) receive a quiet grey
  # outline so the near-invariance claim remains visually honest.
  geom_point(
    data = discordant_points,
    aes(x = pearson, y = ari),
    inherit.aes = FALSE,
    shape = 21,
    size = 6.0,
    stroke = 0.8,
    fill = NA,
    color = COL_DISCORD
  ) +
  # Estimator-specific numerical optimum.
  geom_point(
    data = best_numeric,
    aes(x = pearson, y = ari),
    inherit.aes = FALSE,
    shape = 21,
    size = 6.5,
    stroke = 1.25,
    fill = NA,
    color = COL_NUMERIC
  ) +
  # Shared manuscript biological optimum.
  geom_point(
    data = fig6_points,
    aes(x = pearson, y = ari),
    inherit.aes = FALSE,
    shape = 23,
    size = 5.8,
    stroke = 1.1,
    fill = "white",
    color = COL_BIO
  ) +
  geom_segment(
    data = numeric_labels,
    aes(x = label_x, y = label_y, xend = pearson, yend = ari),
    inherit.aes = FALSE,
    linewidth = 0.5,
    color = COL_NUMERIC
  ) +
  geom_label(
    data = numeric_labels,
    aes(x = label_x, y = label_y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 2.9,
    lineheight = 1.08,
    label.size = 0.22,
    color = COL_TEXT,
    fill = alpha("white", 0.94),
    family = base_family
  ) +
  geom_segment(
    data = bio_labels,
    aes(x = label_x, y = label_y, xend = pearson, yend = ari),
    inherit.aes = FALSE,
    linewidth = 0.5,
    color = COL_BIO
  ) +
  geom_label(
    data = bio_labels,
    aes(x = label_x, y = label_y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 2.9,
    lineheight = 1.08,
    label.size = 0.22,
    color = COL_TEXT,
    fill = alpha("white", 0.94),
    family = base_family
  ) +
  facet_wrap(~ estimator, nrow = 1, scales = "fixed") +
  scale_color_manual(values = k_colors, name = "k-mer size (k)") +
  scale_size_continuous(
    trans = "log10",
    range = c(2.2, 5.2),
    breaks = pretty_breaks(n = 4),
    name = "Window size (w)"
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
    title = "Sequence parameters set tissue recovery; estimator choice shifts calibration",
    subtitle = agreement_note,
    x = "Agreement with exact BEDTools Jaccard (Pearson r)",
    y = "Tissue recovery (adjusted Rand index)",
    caption = paste(
      "Orange ring: estimator-specific numerical optimum.",
      "Diamond: k=10, w=30 biological optimum used in Figure 6.",
      "Grey ring: p=24 cell whose ARI differs between estimators."
    )
  ) +
  theme_classic(base_size = 11, base_family = base_family) +
  theme(
    plot.title = element_text(
      face = "bold", size = 14.2, color = COL_TEXT,
      lineheight = 1.05, margin = margin(b = 4)
    ),
    plot.subtitle = element_text(
      size = 9.7, color = "#4C5661", margin = margin(b = 10)
    ),
    plot.caption = element_text(
      size = 8.6, color = "#59636D", hjust = 0, margin = margin(t = 8)
    ),
    strip.text = element_text(face = "bold", size = 10.5, color = COL_TEXT),
    strip.background = element_blank(),
    axis.title = element_text(face = "bold", color = COL_TEXT),
    axis.text = element_text(color = COL_TEXT),
    axis.line = element_line(color = "#6B747D", linewidth = 0.35),
    axis.ticks = element_line(color = "#6B747D", linewidth = 0.35),
    panel.grid.major = element_line(color = COL_GRID, linewidth = 0.30),
    panel.grid.minor = element_blank(),
    panel.spacing.x = grid::unit(1.2, "lines"),
    legend.position = "right",
    legend.title = element_text(size = 9.3, face = "bold", color = COL_TEXT),
    legend.text = element_text(size = 8.9, color = COL_TEXT),
    plot.margin = margin(12, 18, 12, 12)
  )

CairoPNG(
  filename = out_png,
  width = 12.8,
  height = 6.7,
  units = "in",
  res = 300,
  bg = "white"
)
print(p)
dev.off()

message("Wrote: ", out_png)
