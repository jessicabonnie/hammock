#!/usr/bin/env Rscript

# Figure 4 — Interval-mode inclusion-exclusion Jaccard versus BEDTools
#
# Produces two figures:
#   1. paper/figures/interval_accuracy.png
#      Main-text, single-panel comparison of hammock inclusion-exclusion
#      Jaccard with exact BEDTools Jaccard.
#   2. paper/figures/interval_accuracy_bothmetrics.png
#      Two-panel comparison retaining both inclusion-exclusion and the legacy
#      register-equality statistic for possible supplementary use.
#
# The recoverable caption for the two-metric figure is also written to
# paper/interval_accuracy/interval_accuracy_bothmetrics_caption.txt.

required_packages <- c(
  "dplyr", "readr", "tidyr", "ggplot2", "scales", "patchwork", "Cairo"
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
  library(tidyr)
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
main_png <- if (length(argv) >= 1) {
  normalizePath(argv[1], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "interval_accuracy.png")
}
both_png <- if (length(argv) >= 2) {
  normalizePath(argv[2], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "interval_accuracy_bothmetrics.png")
}
caption_txt <- file.path(
  script_dir, "interval_accuracy_bothmetrics_caption.txt"
)

dir.create(dirname(main_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(both_png), recursive = TRUE, showWarnings = FALSE)

for (path in c(bedtools_tsv, hammock_csvs)) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
}

COL_BEDTOOLS <- "#46515C"
COL_IE <- "#007C83"
COL_RE <- "#D28B35"
COL_GRID <- "#D9DEE3"
COL_TEXT <- "#20262D"
base_family <- "sans"

EST_RE <- "Register-equality (jaccard_similarity)"
EST_IE <- "Inclusion–exclusion (jaccard_similarity_ie)"
EST_LEVELS <- c(EST_IE, EST_RE)
EST_COLORS <- setNames(c(COL_IE, COL_RE), EST_LEVELS)

theme_paper <- function(base_size = 10.5) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        face = "bold", size = rel(1.04), color = COL_TEXT,
        lineheight = 1.05, margin = margin(b = 7)
      ),
      plot.subtitle = element_text(
        size = rel(0.88), color = "#56616B", margin = margin(b = 8)
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
    paste(missing_bedtools, collapse = ", "),
    call. = FALSE
  )
}

bedtools_pairs <- bedtools_raw %>%
  unordered_pairs(file1, file2) %>%
  transmute(.a, .b, is_self, bedtools_jaccard = jaccard)

read_hammock <- function(path, precision) {
  raw <- read_csv(path, show_col_types = FALSE)
  required_hammock <- c(
    "file1", "file2", "jaccard_similarity", "jaccard_similarity_ie"
  )
  missing_hammock <- setdiff(required_hammock, names(raw))
  if (length(missing_hammock) > 0) {
    stop(
      basename(path), " lacks columns: ",
      paste(missing_hammock, collapse = ", "),
      "\nRegenerate these outputs with a current hammock version.",
      call. = FALSE
    )
  }

  raw %>%
    unordered_pairs(file1, file2) %>%
    transmute(
      .a, .b, is_self,
      precision = as.integer(precision),
      `Register-equality (jaccard_similarity)` = jaccard_similarity,
      `Inclusion–exclusion (jaccard_similarity_ie)` = jaccard_similarity_ie
    ) %>%
    pivot_longer(
      cols = all_of(c(EST_RE, EST_IE)),
      names_to = "estimator",
      values_to = "hammock_jaccard"
    )
}

hammock_pairs <- bind_rows(
  lapply(names(hammock_csvs), function(p) {
    read_hammock(hammock_csvs[[p]], as.integer(p))
  })
)

cross <- inner_join(
  hammock_pairs,
  bedtools_pairs,
  by = c(".a", ".b", "is_self")
) %>%
  filter(!is_self) %>%
  mutate(
    estimator = factor(estimator, levels = EST_LEVELS),
    precision_label = factor(
      sprintf("p = %d", precision),
      levels = sprintf("p = %d", PRECISIONS)
    ),
    gap = hammock_jaccard - bedtools_jaccard,
    abs_gap = abs(gap)
  )

if (nrow(cross) == 0) {
  stop("No off-diagonal pairs joined between hammock and BEDTools.", call. = FALSE)
}

coverage <- cross %>% count(precision, estimator, name = "n")
if (length(unique(coverage$n)) != 1) {
  stop("Precisions or estimators cover different pair sets.", call. = FALSE)
}

stats <- cross %>%
  group_by(precision, estimator) %>%
  summarise(
    n = n(),
    pearson = cor(hammock_jaccard, bedtools_jaccard, method = "pearson"),
    spearman = cor(hammock_jaccard, bedtools_jaccard, method = "spearman"),
    kendall = cor(hammock_jaccard, bedtools_jaccard, method = "kendall"),
    mae = mean(abs_gap),
    rmse = sqrt(mean(gap^2)),
    .groups = "drop"
  )

stats_csv <- file.path(script_dir, "interval_accuracy_stats.csv")
write_csv(stats, stats_csv)
message("Per-precision agreement; self-comparisons excluded:")
print(as.data.frame(stats), digits = 5)
message("Wrote: ", stats_csv)

format_stats <- function(row, include_name = FALSE) {
  prefix <- if (include_name) paste0(as.character(row$estimator), "\n") else ""
  paste0(
    prefix,
    sprintf("Pearson r = %.5f\n", row$pearson),
    sprintf("Spearman ρ = %.5f\n", row$spearman),
    sprintf("Kendall τ = %.4f\n", row$kendall),
    sprintf("MAE = %.5f\n", row$mae),
    sprintf("n = %d pairs", row$n)
  )
}

# ---------------------------------------------------------------------------
# Main-text Figure 4: inclusion-exclusion only
# ---------------------------------------------------------------------------
main_points <- cross %>%
  filter(precision == REFERENCE_PRECISION, estimator == EST_IE)
main_stats <- stats %>%
  filter(precision == REFERENCE_PRECISION, estimator == EST_IE)
if (nrow(main_stats) != 1) {
  stop("Expected one inclusion-exclusion statistics row at p = ",
       REFERENCE_PRECISION, call. = FALSE)
}

main_range <- range(
  c(main_points$bedtools_jaccard, main_points$hammock_jaccard),
  finite = TRUE
)
main_pad <- diff(main_range) * 0.08
if (!is.finite(main_pad) || main_pad <= 0) main_pad <- 0.02
main_limits <- c(
  max(0, main_range[1] - main_pad),
  min(1, main_range[2] + main_pad)
)

main_figure <- ggplot(
  main_points,
  aes(x = bedtools_jaccard, y = hammock_jaccard)
) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "22", linewidth = 0.55, color = "#8A939C"
  ) +
  geom_point(
    shape = 21, size = 2.35, stroke = 0.35,
    color = "white", fill = COL_IE, alpha = 0.70
  ) +
  geom_smooth(
    method = "loess", formula = y ~ x, span = 0.95,
    se = FALSE, linewidth = 0.9, color = COL_IE
  ) +
  annotate(
    "label",
    x = main_limits[1] + 0.025 * diff(main_limits),
    y = main_limits[2] - 0.025 * diff(main_limits),
    label = format_stats(main_stats),
    hjust = 0, vjust = 1,
    size = 3.05, lineheight = 1.08,
    linewidth = 0, fill = alpha("white", 0.90), color = COL_TEXT
  ) +
  coord_equal(xlim = main_limits, ylim = main_limits, expand = FALSE) +
  scale_x_continuous(labels = label_number(accuracy = 0.02)) +
  scale_y_continuous(labels = label_number(accuracy = 0.02)) +
  labs(
    title = "Hammock inclusion–exclusion Jaccard reproduces exact interval overlap",
    subtitle = sprintf(
      "Maurano fetal DNase hypersensitivity data; HLL precision p = %d",
      REFERENCE_PRECISION
    ),
    x = "BEDTools exact base-pair Jaccard",
    y = "Hammock inclusion–exclusion Jaccard"
  ) +
  theme_paper(base_size = 11.5) +
  theme(legend.position = "none")

CairoPNG(
  filename = main_png,
  width = 6.7,
  height = 6.2,
  units = "in",
  res = 300,
  bg = "white"
)
print(main_figure)
dev.off()
message("Wrote: ", main_png)

# ---------------------------------------------------------------------------
# Recoverable two-metric figure for possible supplementary use
# ---------------------------------------------------------------------------
both_points <- cross %>% filter(precision == REFERENCE_PRECISION)
both_stats <- stats %>% filter(precision == REFERENCE_PRECISION)

both_range <- range(
  c(both_points$bedtools_jaccard, both_points$hammock_jaccard),
  finite = TRUE
)
both_pad <- diff(both_range) * 0.08
if (!is.finite(both_pad) || both_pad <= 0) both_pad <- 0.02
both_limits <- c(
  max(0, both_range[1] - both_pad),
  min(1, both_range[2] + both_pad)
)

annotation_both <- paste(
  vapply(
    EST_LEVELS,
    function(e) {
      format_stats(both_stats %>% filter(estimator == e), include_name = TRUE)
    },
    character(1)
  ),
  collapse = "\n\n"
)

panel_a <- ggplot(
  both_points,
  aes(
    x = bedtools_jaccard,
    y = hammock_jaccard,
    color = estimator,
    fill = estimator
  )
) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "22", linewidth = 0.5, color = "#8A939C"
  ) +
  geom_point(
    shape = 21, size = 2.15, stroke = 0.3,
    color = "white", alpha = 0.62
  ) +
  geom_smooth(
    method = "loess", formula = y ~ x, span = 0.95,
    se = FALSE, linewidth = 0.85
  ) +
  annotate(
    "label",
    x = both_limits[1] + 0.025 * diff(both_limits),
    y = both_limits[2] - 0.025 * diff(both_limits),
    label = annotation_both,
    hjust = 0, vjust = 1,
    size = 2.35, lineheight = 1.03,
    linewidth = 0, fill = alpha("white", 0.90), color = COL_TEXT
  ) +
  scale_color_manual(values = EST_COLORS, drop = FALSE) +
  scale_fill_manual(values = EST_COLORS, drop = FALSE) +
  coord_equal(xlim = both_limits, ylim = both_limits, expand = FALSE) +
  scale_x_continuous(labels = label_number(accuracy = 0.05)) +
  scale_y_continuous(labels = label_number(accuracy = 0.05)) +
  labs(
    title = sprintf("A  Metric behavior at p = %d", REFERENCE_PRECISION),
    x = "BEDTools exact base-pair Jaccard",
    y = "Hammock interval-mode estimate"
  ) +
  guides(
    color = guide_legend(
      nrow = 2,
      override.aes = list(alpha = 1, size = 2.6)
    ),
    fill = "none"
  ) +
  theme_paper() +
  theme(legend.position = "top", legend.margin = margin(b = 2))

panel_b <- ggplot(
  cross,
  aes(x = bedtools_jaccard, y = abs_gap)
) +
  geom_point(aes(color = estimator), size = 1.2, alpha = 0.16) +
  geom_smooth(
    aes(
      linetype = precision_label,
      group = interaction(precision_label, estimator),
      color = estimator
    ),
    method = "loess", formula = y ~ x, span = 0.95,
    se = FALSE, linewidth = 0.85
  ) +
  scale_linetype_manual(values = c("solid", "22", "42")) +
  scale_color_manual(values = EST_COLORS, drop = FALSE, guide = "none") +
  scale_x_continuous(labels = label_number(accuracy = 0.05)) +
  scale_y_log10(
    labels = label_log(digits = 2),
    breaks = 10^(-6:0),
    expand = expansion(mult = c(0.05, 0.10))
  ) +
  annotation_logticks(
    sides = "l", size = 0.3, color = "#8A939C",
    short = unit(0.04, "cm"), mid = unit(0.07, "cm"),
    long = unit(0.11, "cm")
  ) +
  labs(
    title = "B  Absolute deviation across HLL precisions",
    x = "BEDTools exact base-pair Jaccard",
    y = "|Hammock − BEDTools|  (log scale)"
  ) +
  guides(
    linetype = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(color = COL_TEXT, linewidth = 0.9)
    )
  ) +
  theme_paper() +
  theme(legend.position = "top", legend.margin = margin(b = 4))

both_figure <- panel_a + panel_b +
  plot_layout(widths = c(1, 1.05)) +
  plot_annotation(
    title = paste(
      "Inclusion–exclusion estimates BEDTools set Jaccard, while",
      "register equality is a distinct compatibility statistic"
    ),
    theme = theme(
      plot.title = element_text(
        family = base_family, face = "bold", size = 13.2,
        color = COL_TEXT, margin = margin(b = 10)
      ),
      plot.margin = margin(12, 16, 12, 12)
    )
  )

CairoPNG(
  filename = both_png,
  width = 11.6,
  height = 5.9,
  units = "in",
  res = 300,
  bg = "white"
)
print(both_figure)
dev.off()
message("Wrote: ", both_png)

both_caption <- paste0(
  "Figure S?. Inclusion–exclusion and register-equality statistics have ",
  "different relationships to exact interval Jaccard. Pairwise comparisons ",
  "were performed on 20 Maurano fetal-tissue DNase hypersensitivity BED files, ",
  "yielding 190 unique off-diagonal pairs; self-comparisons were excluded. ",
  "BEDTools reports exact set Jaccard over covered base pairs. Hammock interval ",
  "mode reports both an inclusion–exclusion estimate of the same set Jaccard, ",
  "computed from HyperLogLog cardinality estimates of A, B, and their union, ",
  "and a register-equality statistic retained for compatibility with the ",
  "original hammock implementation. (A) At HyperLogLog precision p = ",
  REFERENCE_PRECISION,
  ", inclusion–exclusion estimates lie near the identity line with BEDTools, ",
  "whereas register equality is shifted upward because equal register values ",
  "can arise without equality of the underlying sets. The inset reports ",
  "Pearson correlation, Spearman correlation, Kendall rank correlation, mean ",
  "absolute error, and the number of pairs for each statistic. (B) Absolute ",
  "deviation from BEDTools is shown across precisions p = 18, 21, and 23 on a ",
  "logarithmic scale. Inclusion–exclusion estimates the same target quantity ",
  "as BEDTools and remain close to zero deviation, whereas register equality ",
  "retains a substantially larger definitional gap. This comparison is ",
  "provided to distinguish the recommended set-Jaccard estimate from the ",
  "legacy compatibility statistic; the main-text analysis uses ",
  "jaccard_similarity_ie."
)
writeLines(both_caption, caption_txt, useBytes = TRUE)
message("Wrote: ", caption_txt)
