#!/usr/bin/env Rscript

# Figure 4 — Interval-mode similarity versus BEDTools Jaccard
#
# Per-pair agreement between hammock's interval-mode estimate and the exact
# BEDTools bp set-Jaccard on the Maurano fetal-tissue DHS collection
# (20 BED files; 190 unique off-diagonal pairs).
#
#   Panel A: scatter with the y = x diagonal and a smoothed curve showing the
#            empirical mapping f(J), plus Pearson r and Spearman rho.
#   Panel B: the gap (hammock - BEDTools) against BEDTools Jaccard, across all
#            three precisions, showing that it shrinks as overlap grows and is
#            essentially precision-independent.
#
# The two estimators are not defined over the same universe: BEDTools reports
# a set-Jaccard over covered base pairs, while hammock reports a
# register-equality Jaccard over the sketched positions. Chance register ties
# put a positive floor under the latter, so hammock returns f(J) > J with f
# monotone and f(1) = 1. The claim this figure supports is preserved
# *structure* (r, rho), not numerical identity.
#
# Framing follows docs/jaccard-definitional-gap.md, which prescribes the
# scatter + y = x + smoothed-mapping presentation and explicitly warns against
# subtracting a constant offset: the gap is not constant, it goes to 0 at
# J = 1. The self-comparisons are plotted (but excluded from the statistics)
# precisely because they anchor f(1) = 1.
#
# Data (copied into docs/data/ from experiments/maurano_dhs_validation/):
#   hammock_hll_p{18,21,23}_jaccB.csv   -- results/raw_abc/
#   maurano_bedtools_ref.tsv            -- data/
#
# Usage:
#   ml r/4.3.0 libjpeg/9c
#   Rscript paper/interval_accuracy/plot_interval_accuracy.R
#   Rscript paper/interval_accuracy/plot_interval_accuracy.R path/to/output.png

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

# Resolve repository-relative paths.
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

# Shared publication palette (matches paper/pairwise_scaling/).
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
      panel.grid.major = element_line(color = COL_GRID, linewidth = 0.35),
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
# Load and join per-pair estimates
# -----------------------------------------------------------------------------

# Both inputs list every ordered pair including self-comparisons. Key each row
# on the unordered pair so the join is direction-independent, then keep one row
# per unordered pair. Self-comparisons are retained but flagged: they are 1.0
# by construction in both tools, so they must not enter the correlations, but
# they are the one place where the two estimators provably agree.
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

cross <- paired %>% filter(!is_self)
self_pairs <- paired %>% filter(is_self)

n_pairs_per_precision <- cross %>% count(precision, name = "n")
if (length(unique(n_pairs_per_precision$n)) != 1) {
  stop("Precisions cover different pair sets; refusing to plot.", call. = FALSE)
}
n_pairs <- n_pairs_per_precision$n[1]

# Sanity check the anchor: both tools must report 1.0 for a file against
# itself, or the "f(1) = 1" reading of the figure is wrong.
if (nrow(self_pairs) > 0) {
  max_self_dev <- max(abs(c(
    self_pairs$bedtools_jaccard - 1, self_pairs$hammock_jaccard - 1
  )))
  if (max_self_dev > 1e-6) {
    stop("Self-comparisons are not 1.0 (max deviation ",
         signif(max_self_dev, 3), "); check the inputs.", call. = FALSE)
  }
}

stats <- cross %>%
  group_by(precision) %>%
  summarise(
    n = n(),
    pearson = cor(hammock_jaccard, bedtools_jaccard, method = "pearson"),
    spearman = cor(hammock_jaccard, bedtools_jaccard, method = "spearman"),
    mae = mean(abs(gap)),
    gap_at_low_J = gap[which.min(bedtools_jaccard)],
    gap_at_high_J = gap[which.max(bedtools_jaccard)],
    .groups = "drop"
  )

message("Per-precision agreement (", n_pairs,
        " unique off-diagonal pairs; self-comparisons excluded):")
print(as.data.frame(stats), digits = 5)

ref_stats <- stats %>% filter(precision == REFERENCE_PRECISION)
if (nrow(ref_stats) != 1) {
  stop("No statistics for the reference precision p = ", REFERENCE_PRECISION,
       call. = FALSE)
}
ref_points <- cross %>% filter(precision == REFERENCE_PRECISION)
ref_self <- self_pairs %>% filter(precision == REFERENCE_PRECISION)

precision_labels <- sprintf("p = %d", PRECISIONS)
precision_colors <- setNames(
  c(COL_BEDTOOLS, COL_HAMMOCK, COL_COMPARE), precision_labels
)
cross <- cross %>%
  mutate(precision_label = factor(sprintf("p = %d", precision),
                                  levels = precision_labels))

# -----------------------------------------------------------------------------
# Panel A: per-pair scatter and the empirical mapping f(J)
# -----------------------------------------------------------------------------

# The full unit square, so the y = x diagonal and the f(1) = 1 anchor are both
# in frame. The observed pairs occupy only part of it; that is the honest
# picture of where f(J) has actually been characterised.
annotation_a <- sprintf(
  "Pearson r = %.4f\nSpearman ρ = %.4f\nn = %d pairs",
  ref_stats$pearson, ref_stats$spearman, ref_stats$n
)

observed_band <- sprintf(
  "observed range\nJ = %.2f – %.2f",
  min(ref_points$bedtools_jaccard), max(ref_points$bedtools_jaccard)
)

panel_a <- ggplot(ref_points, aes(x = bedtools_jaccard, y = hammock_jaccard)) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "22", linewidth = 0.45, color = "#8A939C"
  ) +
  annotate(
    "text",
    x = 0.985, y = 0.985,
    label = "y = x",
    hjust = 1, vjust = 1.6, size = 3, color = "#8A939C", angle = 45
  ) +
  geom_smooth(
    method = "loess", formula = y ~ x, span = 0.95,
    se = TRUE, level = 0.95,
    color = COL_COMPARE, fill = alpha(COL_COMPARE, 0.18), linewidth = 0.8
  ) +
  geom_point(
    shape = 21, size = 2.1, stroke = 0.3,
    fill = alpha(COL_HAMMOCK, 0.55), color = "white"
  ) +
  geom_point(
    data = ref_self,
    aes(x = bedtools_jaccard, y = hammock_jaccard),
    shape = 21, size = 2.6, stroke = 0.7,
    fill = "white", color = COL_TEXT
  ) +
  annotate(
    "text",
    x = 1, y = 1,
    label = "self-comparisons\n(J = 1 in both tools)",
    hjust = 1.06, vjust = 0.9, size = 2.85, lineheight = 1.05, color = COL_TEXT
  ) +
  annotate(
    "label",
    x = 0.02, y = 0.98,
    label = annotation_a,
    hjust = 0, vjust = 1,
    size = 3.05, lineheight = 1.05,
    label.size = 0, fill = alpha("white", 0.85), color = COL_TEXT
  ) +
  annotate(
    "text",
    x = mean(range(ref_points$bedtools_jaccard)), y = 0.045,
    label = observed_band,
    hjust = 0.5, vjust = 0, size = 2.85, lineheight = 1.05, color = "#56616C"
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "A  Monotone, but not the identity",
    subtitle = sprintf(
      "Maurano fetal-tissue DHS; 20 BED files;\ninterval mode; p = %d; loess fit with 95%% CI",
      REFERENCE_PRECISION
    ),
    x = "BEDTools Jaccard (exact, bp set)",
    y = "hammock interval-mode Jaccard"
  ) +
  theme_paper() +
  theme(legend.position = "none")

# -----------------------------------------------------------------------------
# Panel B: how the gap varies with overlap and precision
# -----------------------------------------------------------------------------
panel_b <- ggplot(
  cross,
  aes(x = bedtools_jaccard, y = gap, color = precision_label)
) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "#8A939C") +
  geom_point(size = 1.7, alpha = 0.5) +
  geom_smooth(
    method = "loess", formula = y ~ x, span = 0.95, se = FALSE,
    color = COL_TEXT, linewidth = 0.7
  ) +
  annotate(
    "text",
    x = max(cross$bedtools_jaccard), y = 0,
    label = "no gap",
    hjust = 1, vjust = -0.7, size = 2.85, color = "#56616C"
  ) +
  scale_color_manual(values = precision_colors) +
  scale_x_continuous(labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.05),
    limits = c(0, NA), expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(
    title = "B  The gap closes as overlap grows",
    subtitle = paste0(
      "Chance register ties inflate low-overlap pairs most;\n",
      "the gap is nearly independent of sketch precision"
    ),
    x = "BEDTools Jaccard (exact, bp set)",
    y = "hammock − BEDTools"
  ) +
  guides(color = guide_legend(
    nrow = 1, byrow = TRUE, override.aes = list(alpha = 1, size = 2.4)
  )) +
  theme_paper() +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.box.just = "center",
    legend.margin = margin(b = 2)
  )

# Assemble and save with the same CairoPNG pattern used by docs/scripts/.
figure <- panel_a + panel_b +
  plot_layout(widths = c(1, 1.05)) +
  plot_annotation(
    title = "Interval sketches preserve exact-overlap similarity structure",
    subtitle = paste0(
      "Rank and linear agreement with BEDTools are near-perfect and stable ",
      "across precision; the two tools estimate different quantities, so the ",
      "mapping is monotone rather than the identity."
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
  )

CairoPNG(
  filename = out_png,
  width = 11.0,
  height = 5.8,
  units = "in",
  res = 300,
  bg = "white"
)
print(figure)
dev.off()

message("Wrote: ", out_png)
