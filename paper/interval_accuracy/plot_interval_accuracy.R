#!/usr/bin/env Rscript

# Figure 4 — Interval-mode similarity versus BEDTools Jaccard
# Creates a publication-ready two-panel PNG using CairoPNG.

required_packages <- c(
  "dplyr", "readr", "ggplot2", "scales", "patchwork", "Cairo"
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
  library(Cairo)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
  stop("Could not determine the script path. Run with Rscript.", call. = FALSE)
}
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

data_dir <- file.path(repo_root, "docs", "data")
bedtools_tsv <- file.path(data_dir, "maurano_bedtools_ref.tsv")
PRECISIONS <- c(18, 21, 23)
REFERENCE_PRECISION <- 21
hammock_csvs <- setNames(
  file.path(data_dir, sprintf("hammock_hll_p%d_jaccB.csv", PRECISIONS)),
  PRECISIONS
)

argv <- commandArgs(trailingOnly = TRUE)
out_png <- if (length(argv) >= 1) {
  normalizePath(argv[1], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "interval_accuracy.png")
}
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

for (path in c(bedtools_tsv, hammock_csvs)) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
}

COL_BEDTOOLS <- "#46515C"
COL_HAMMOCK <- "#007C83"
COL_COMPARE <- "#D28B35"
COL_GRID <- "#D9DEE3"
COL_TEXT <- "#20262D"
base_family <- "sans"

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
      panel.grid.major = element_line(color = COL_GRID, linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.justification = "center",
      legend.box.just = "center",
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.82), color = COL_TEXT),
      legend.key.width = grid::unit(1.25, "lines"),
      plot.margin = margin(10, 12, 10, 10)
    )
}

unordered_pairs <- function(df, file1, file2) {
  df %>%
    mutate(
      .a = pmin({{ file1 }}, {{ file2 }}),
      .b = pmax({{ file1 }}, {{ file2 }}),
      is_self = .a == .b
    ) %>%
    distinct(.a, .b, .keep_all = TRUE)
}

bedtools_raw <- read_tsv(bedtools_tsv, show_col_types = FALSE)
required_bedtools <- c("file1", "file2", "jaccard")
missing_bedtools <- setdiff(required_bedtools, names(bedtools_raw))
if (length(missing_bedtools) > 0) {
  stop(
    "BEDTools reference lacks columns: ",
    paste(missing_bedtools, collapse = ", "), call. = FALSE
  )
}

bedtools_pairs <- bedtools_raw %>%
  unordered_pairs(file1, file2) %>%
  transmute(.a, .b, is_self, bedtools_jaccard = jaccard)

read_hammock <- function(path, precision) {
  raw <- read_csv(path, show_col_types = FALSE)
  required_hammock <- c("file1", "file2", "jaccard_similarity")
  missing_hammock <- setdiff(required_hammock, names(raw))
  if (length(missing_hammock) > 0) {
    stop(
      basename(path), " lacks columns: ",
      paste(missing_hammock, collapse = ", "), call. = FALSE
    )
  }
  raw %>%
    unordered_pairs(file1, file2) %>%
    transmute(.a, .b, precision = precision, hammock_jaccard = jaccard_similarity)
}

hammock_pairs <- bind_rows(
  lapply(names(hammock_csvs), function(p) {
    read_hammock(hammock_csvs[[p]], as.integer(p))
  })
)

paired <- inner_join(hammock_pairs, bedtools_pairs, by = c(".a", ".b")) %>%
  mutate(gap = hammock_jaccard - bedtools_jaccard)

if (nrow(paired) == 0) {
  stop("No pairs joined between hammock and BEDTools inputs.", call. = FALSE)
}

cross <- paired %>%
  filter(!is_self) %>%
  mutate(
    precision_label = factor(
      sprintf("p = %d", precision),
      levels = sprintf("p = %d", PRECISIONS)
    )
  )

n_pairs_per_precision <- cross %>% count(precision, name = "n")
if (length(unique(n_pairs_per_precision$n)) != 1) {
  stop("Precisions cover different pair sets; refusing to plot.", call. = FALSE)
}

stats <- cross %>%
  group_by(precision) %>%
  summarise(
    n = n(),
    pearson = cor(hammock_jaccard, bedtools_jaccard, method = "pearson"),
    spearman = cor(hammock_jaccard, bedtools_jaccard, method = "spearman"),
    kendall = cor(hammock_jaccard, bedtools_jaccard, method = "kendall"),
    mae = mean(abs(gap)),
    .groups = "drop"
  )

message("Per-precision agreement; self-comparisons excluded:")
print(as.data.frame(stats), digits = 5)

ref_stats <- stats %>% filter(precision == REFERENCE_PRECISION)
if (nrow(ref_stats) != 1) {
  stop("No statistics for p = ", REFERENCE_PRECISION, call. = FALSE)
}
ref_points <- cross %>% filter(precision == REFERENCE_PRECISION)

annotation_a <- sprintf(
  "Spearman ρ = %.4f\nKendall τ = %.4f\nPearson r = %.4f\nn = %d",
  ref_stats$spearman,
  ref_stats$kendall,
  ref_stats$pearson,
  ref_stats$n
)

# Crop Panel A to the observed off-diagonal region, with equal axis limits so
# the y = x line remains interpretable without allowing self-comparisons at
# (1, 1) to compress the data cloud.
combined_range <- range(
  c(ref_points$bedtools_jaccard, ref_points$hammock_jaccard),
  finite = TRUE
)
pad <- diff(combined_range) * 0.08
if (!is.finite(pad) || pad <= 0) pad <- 0.02
plot_limits <- c(
  max(0, combined_range[1] - pad),
  min(1, combined_range[2] + pad)
)

panel_a <- ggplot(ref_points, aes(x = bedtools_jaccard, y = hammock_jaccard)) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "22", linewidth = 0.5, color = "#8A939C"
  ) +
  geom_smooth(
    method = "loess", formula = y ~ x, span = 0.95,
    se = TRUE, level = 0.95,
    color = COL_COMPARE, fill = alpha(COL_COMPARE, 0.17), linewidth = 0.85
  ) +
  geom_point(
    shape = 21, size = 2.15, stroke = 0.3,
    fill = alpha(COL_HAMMOCK, 0.58), color = "white"
  ) +
  annotate(
    "label",
    x = plot_limits[1] + 0.025 * diff(plot_limits),
    y = plot_limits[2] - 0.025 * diff(plot_limits),
    label = annotation_a,
    hjust = 0, vjust = 1,
    size = 2.85, lineheight = 1.04,
    linewidth = 0, fill = alpha("white", 0.88), color = COL_TEXT
  ) +
  coord_equal(xlim = plot_limits, ylim = plot_limits, expand = FALSE) +
  scale_x_continuous(labels = label_number(accuracy = 0.05)) +
  scale_y_continuous(labels = label_number(accuracy = 0.05)) +
  labs(
    title = "A  Pairwise similarity is strongly preserved",
    x = "BEDTools Jaccard",
    y = "Hammock interval-mode Jaccard"
  ) +
  theme_paper() +
  theme(legend.position = "none")

# Precision is encoded by line type rather than color so that the panel does
# not visually exaggerate the small differences among precision settings.
panel_b <- ggplot(cross, aes(x = bedtools_jaccard, y = gap)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "#8A939C") +
  geom_point(color = COL_BEDTOOLS, size = 1.35, alpha = 0.18) +
  geom_smooth(
    aes(linetype = precision_label, group = precision_label),
    method = "loess", formula = y ~ x, span = 0.95,
    se = FALSE, color = COL_HAMMOCK, linewidth = 0.85
  ) +
  scale_linetype_manual(values = c("solid", "22", "42")) +
  scale_x_continuous(labels = label_number(accuracy = 0.05)) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.05),
    limits = c(0, NA),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  labs(
    title = "B  The numerical gap decreases with overlap",
    x = "BEDTools Jaccard",
    y = "Hammock − BEDTools"
  ) +
  guides(
    linetype = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(color = COL_HAMMOCK, linewidth = 0.9)
    )
  ) +
  theme_paper() +
  theme(
    legend.position = "top",
    legend.margin = margin(b = 4)
  )

figure <- panel_a + panel_b +
  plot_layout(widths = c(1, 1.05)) +
  plot_annotation(
    title = "Interval sketches preserve pairwise similarity structure",
    theme = theme(
      plot.title = element_text(
        family = base_family, face = "bold", size = 15,
        color = COL_TEXT, margin = margin(b = 10)
      ),
      plot.margin = margin(12, 16, 12, 12)
    )
  )

CairoPNG(
  filename = out_png,
  width = 11.6,
  height = 5.9,
  units = "in",
  res = 300,
  bg = "white"
)
print(figure)
dev.off()

message("Wrote: ", out_png)
