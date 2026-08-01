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

# The two estimators plotted side by side. `jaccard_similarity` is
# register-equality: the fraction of active registers whose values agree. It is
# not set Jaccard -- registers tie by chance, which puts a floor under it -- so
# its offset from y = x is definitional, not error. `jaccard_similarity_ie`
# (hammock >= 0.5.0) is the inclusion-exclusion estimate |A|+|B|-|A∪B| over the
# union, which does estimate the same quantity BEDTools reports.
EST_RE <- "Register-equality (jaccard_similarity)"
EST_IE <- "Inclusion\u2013exclusion (jaccard_similarity_ie)"
EST_LEVELS <- c(EST_RE, EST_IE)

read_hammock <- function(path, precision) {
  raw <- read_csv(path, show_col_types = FALSE)
  required_hammock <- c("file1", "file2",
                        "jaccard_similarity", "jaccard_similarity_ie")
  missing_hammock <- setdiff(required_hammock, names(raw))
  if (length(missing_hammock) > 0) {
    stop(
      basename(path), " lacks columns: ",
      paste(missing_hammock, collapse = ", "),
      "\n  CSVs written before 2026-05-14 carry the placeholder `containment` ",
      "column instead of the metric block and cannot supply the IE estimator; ",
      "regenerate them with hammock >= 0.5.0.", call. = FALSE
    )
  }
  base <- raw %>% unordered_pairs(file1, file2)
  bind_rows(
    base %>% transmute(.a, .b, precision = precision,
                       estimator = EST_RE, hammock_jaccard = jaccard_similarity),
    base %>% transmute(.a, .b, precision = precision,
                       estimator = EST_IE, hammock_jaccard = jaccard_similarity_ie)
  )
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
    ),
    estimator = factor(estimator, levels = EST_LEVELS)
  )

n_pairs_per_precision <- cross %>% count(precision, estimator, name = "n")
if (length(unique(n_pairs_per_precision$n)) != 1) {
  stop("Precisions cover different pair sets; refusing to plot.", call. = FALSE)
}

stats <- cross %>%
  group_by(precision, estimator) %>%
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
if (nrow(ref_stats) != length(EST_LEVELS)) {
  stop("No statistics for p = ", REFERENCE_PRECISION, call. = FALSE)
}
ref_points <- cross %>% filter(precision == REFERENCE_PRECISION)

fmt_stats <- function(row) {
  sprintf(
    "%s\n  Pearson r = %.5f   Kendall \u03c4 = %.4f\n  MAE = %.4f",
    sub(" \\(.*", "", row$estimator), row$pearson, row$kendall, row$mae
  )
}
annotation_a <- paste(
  vapply(EST_LEVELS,
         function(e) fmt_stats(ref_stats %>% filter(estimator == e)),
         character(1)),
  collapse = "\n"
)
annotation_a <- paste0(annotation_a, sprintf("\nn = %d pairs", ref_stats$n[1]))

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

EST_COLORS <- setNames(c(COL_HAMMOCK, COL_COMPARE), EST_LEVELS)

panel_a <- ggplot(ref_points,
                  aes(x = bedtools_jaccard, y = hammock_jaccard,
                      color = estimator, fill = estimator)) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "22", linewidth = 0.5, color = "#8A939C"
  ) +
  geom_point(shape = 21, size = 2.15, stroke = 0.3, color = "white",
             alpha = 0.62) +
  geom_smooth(
    method = "loess", formula = y ~ x, span = 0.95,
    se = FALSE, linewidth = 0.85
  ) +
  annotate(
    "label",
    x = plot_limits[1] + 0.025 * diff(plot_limits),
    y = plot_limits[2] - 0.025 * diff(plot_limits),
    label = annotation_a,
    hjust = 0, vjust = 1,
    size = 2.5, lineheight = 1.06,
    linewidth = 0, fill = alpha("white", 0.88), color = COL_TEXT
  ) +
  scale_color_manual(values = EST_COLORS, drop = FALSE) +
  scale_fill_manual(values = EST_COLORS, drop = FALSE) +
  coord_equal(xlim = plot_limits, ylim = plot_limits, expand = FALSE) +
  scale_x_continuous(labels = label_number(accuracy = 0.05)) +
  scale_y_continuous(labels = label_number(accuracy = 0.05)) +
  labs(
    title = sprintf("A  Inclusion\u2013exclusion lands on y = x (p = %d)",
                    REFERENCE_PRECISION),
    x = "BEDTools Jaccard",
    y = "Hammock interval-mode estimate"
  ) +
  guides(color = guide_legend(nrow = 2, override.aes = list(alpha = 1, size = 2.6)),
         fill = "none") +
  theme_paper() +
  theme(legend.position = "top", legend.margin = margin(b = 2))

# Precision is encoded by line type rather than color so that the panel does
# not visually exaggerate the small differences among precision settings.
# Absolute gap on a log axis. A signed linear axis cannot show this panel's
# point: the register-equality gap (~0.14) and the IE gap (~2e-4) differ by
# nearly three orders of magnitude, so on a linear scale the three IE
# precision curves collapse onto zero and the panel silently fails to display
# the convergence its title claims. No gap is exactly zero (min 1.6e-6 across
# all three precisions), so the log transform drops nothing.
cross_abs <- cross %>% mutate(abs_gap = abs(gap))

panel_b <- ggplot(cross_abs, aes(x = bedtools_jaccard, y = abs_gap)) +
  geom_point(aes(color = estimator), size = 1.2, alpha = 0.16) +
  geom_smooth(
    aes(linetype = precision_label, group = interaction(precision_label, estimator),
        color = estimator),
    method = "loess", formula = y ~ x, span = 0.95,
    se = FALSE, linewidth = 0.85
  ) +
  scale_linetype_manual(values = c("solid", "22", "42")) +
  scale_color_manual(values = EST_COLORS, drop = FALSE, guide = "none") +
  scale_x_continuous(labels = label_number(accuracy = 0.05)) +
  scale_y_log10(
    labels = label_log(digits = 2),
    breaks = 10^(-5:0),
    expand = expansion(mult = c(0.05, 0.10))
  ) +
  annotation_logticks(sides = "l", size = 0.3, color = "#8A939C",
                      short = unit(0.04, "cm"), mid = unit(0.07, "cm"),
                      long = unit(0.11, "cm")) +
  labs(
    title = "B  The register-equality gap is a floor, not an error",
    x = "BEDTools Jaccard",
    y = "|Hammock − BEDTools|  (log scale)"
  ) +
  guides(
    linetype = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(color = COL_TEXT, linewidth = 0.9)
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
    title = paste("Inclusion\u2013exclusion reproduces BEDTools Jaccard;",
                  "register-equality only its ordering"),
    theme = theme(
      plot.title = element_text(
        family = base_family, face = "bold", size = 13.5,
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
