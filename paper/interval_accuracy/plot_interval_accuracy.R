#!/usr/bin/env Rscript

# Figure 4 — Interval-mode accuracy and precision/runtime frontier
#
# Produces two figures:
#   1. paper/figures/interval_accuracy.png
#      Main-text, two-panel figure: inclusion-exclusion Jaccard versus exact
#      BEDTools Jaccard (A), and the precision/runtime frontier (B).
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
precision_frontier_csv <- file.path(
  data_dir, "sweep_precision_maurano_p18_t16.csv"
)
PRECISIONS <- c(18, 21, 23)
# "_full" tag: this script reads both reg_eq_similarity (register equality)
# and jaccard_similarity_ie, which only the full metrics block (--metrics)
# emits together (python/hammock/outprefix.py; the file was renamed to match
# in the metrics-column restructure, docs/seed-metrics-column-restructure.md
# -- same content, tag added).
hammock_csvs <- setNames(
  file.path(data_dir, sprintf("hammock_hll_p%d_jaccB_full.csv", PRECISIONS)),
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
REFERENCE_PRECISION <- if (length(argv) >= 3) as.integer(argv[3]) else 21L
if (is.na(REFERENCE_PRECISION) || !REFERENCE_PRECISION %in% PRECISIONS) {
  stop("Reference precision must be one of: ", paste(PRECISIONS, collapse = ", "),
       call. = FALSE)
}
caption_txt <- if (REFERENCE_PRECISION == 21L) {
  file.path(script_dir, "interval_accuracy_bothmetrics_caption.txt")
} else {
  file.path(
    script_dir,
    sprintf("interval_accuracy_bothmetrics_caption_p%d.txt", REFERENCE_PRECISION)
  )
}

dir.create(dirname(main_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(both_png), recursive = TRUE, showWarnings = FALSE)

for (path in c(bedtools_tsv, hammock_csvs, precision_frontier_csv)) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
}

COL_BEDTOOLS <- "#46515C"
COL_IE <- "#007C83"
COL_RE <- "#D28B35"
COL_GRID <- "#D9DEE3"
COL_TEXT <- "#20262D"
base_family <- "sans"

EST_RE <- "Register-equality (reg_eq_similarity)"
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
      # Anchor title/subtitle to the plot edge, not the panel edge; the panel
      # edge sits ~1 in right of it once the y-axis label and ticks are drawn,
      # and long titles run off the device from there.
      plot.title.position = "plot",
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

# Prefer the renamed register-equality column; archived pre-Step-1 CSVs only
# have jaccard_similarity, so fall back to it and log that we did (once per
# script run, not once per input file) -- see
# docs/seed-jaccard-reg-eq-rename.md Step 2.
.sim_col_fallback_logged <- FALSE
resolve_sim_col <- function(df, context = NULL) {
  if ("reg_eq_similarity" %in% names(df)) return("reg_eq_similarity")
  if (!.sim_col_fallback_logged) {
    message(
      "reg_eq_similarity not found",
      if (!is.null(context)) paste0(" (", context, ")") else "",
      "; falling back to jaccard_similarity"
    )
    .sim_col_fallback_logged <<- TRUE
  }
  "jaccard_similarity"
}

read_hammock <- function(path, precision) {
  raw <- read_csv(path, show_col_types = FALSE)
  sim_col <- resolve_sim_col(raw, basename(path))
  required_hammock <- c(
    "file1", "file2", sim_col, "jaccard_similarity_ie"
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
      !!EST_RE := .data[[sim_col]],
      !!EST_IE := jaccard_similarity_ie
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

message("Per-precision agreement; self-comparisons excluded:")
print(as.data.frame(stats), digits = 5)

# ---------------------------------------------------------------------------
# Provenance for the numbers quoted in docs/paper_outline.md and the
# Supplementary Note S4 in paper/outline.md.
#
# Everything in this section is reporting, not plotting: it exists so that the
# discordance count, the largest inverted gap, the chance-agreement floor and
# the residual-vs-size-ratio correlation have a committed generator instead of
# living only in prose. It was dropped by accident in commit d00b611 (a
# figure-splitting rewrite) and restored on 2026-08-08 -- the paper cites
# 439/17,955 (2.45%), the 48 IE inversions (0.27%), the ~0.14 offset and the
# DeltaJ <= 0.025 bound, none of which any other file in the repo can produce.
# Do not delete it again without re-pointing those citations.
#
# Conventions are pinned here because they move the reported digits:
#
#   * pairs are unordered and self-comparisons are dropped, so n = 190 pairs
#     and C(190, 2) = 17,955 comparisons;
#   * discordance is counted over comparison *pairs*, not rows;
#   * tau-a and tau-b are both reported. They coincide when no two pairs tie,
#     which holds for BEDTools truth here but need not hold for a censored
#     estimator -- an exact 0.0 in the IE column is the `>= 0` clamp in
#     HLLSketch::intersection_size firing, and clamped rows tie with each other.
# ---------------------------------------------------------------------------

rank_agreement <- function(truth, est) {
  dt <- outer(truth, truth, "-")
  de <- outer(est, est, "-")
  upper <- upper.tri(dt)
  dt <- dt[upper]
  de <- de[upper]
  s <- sign(dt) * sign(de)
  n_comparisons <- length(s)
  n_disc <- sum(s < 0)
  n_tied <- sum(s == 0)
  # Largest true gap that the estimator puts the wrong way round.
  worst <- if (n_disc > 0) max(abs(dt[s < 0])) else 0
  tibble(
    comparisons = n_comparisons,
    discordant = n_disc,
    discordant_pct = 100 * n_disc / n_comparisons,
    tied = n_tied,
    tau_a = (n_comparisons - 2 * n_disc - n_tied) / n_comparisons,
    tau_b = cor(truth, est, method = "kendall"),
    max_inverted_gap = worst
  )
}

rank_stats <- cross %>%
  group_by(precision, estimator) %>%
  group_modify(~ rank_agreement(.x$bedtools_jaccard, .x$hammock_jaccard)) %>%
  ungroup() %>%
  left_join(
    cross %>%
      group_by(precision, estimator) %>%
      summarise(clamped_at_zero = sum(hammock_jaccard == 0), .groups = "drop"),
    by = c("precision", "estimator")
  )

message("\nRank agreement against BEDTools (unordered, off-diagonal):")
print(as.data.frame(rank_stats), digits = 5)

# --- chance-agreement floor c, three ways -----------------------------------
# The register-equality column is approximately c + (1 - c) * J. The three fits
# disagree in the second digit and the document quotes all three, so compute
# them rather than restate them. `constrained` is the one-parameter fit that
# respects the (1, 1) fixed point; `pointwise` inverts the relation per pair and
# averages, which is what an over-dispersed reading of the same model gives.
floor_fits <- function(j_bt, j_re) {
  ols <- coef(lm(j_re ~ j_bt))
  num <- sum((1 - j_bt) * (j_re - j_bt))
  den <- sum((1 - j_bt)^2)
  tibble(
    c_ols_intercept = unname(ols[1]),
    c_constrained = num / den,
    c_pointwise = mean((j_re - j_bt) / (1 - j_bt))
  )
}

floor_stats <- cross %>%
  filter(estimator == EST_RE) %>%
  group_by(precision) %>%
  group_modify(~ floor_fits(.x$bedtools_jaccard, .x$hammock_jaccard)) %>%
  ungroup()

message("\nChance-agreement floor c for the register-equality column:")
print(as.data.frame(floor_stats), digits = 5)

# --- residual versus cardinality ratio --------------------------------------
# c depends on |A|/|B|, so the register-equality residual should track the size
# imbalance of the pair. The ratio needs no sketching: the BEDTools reference's
# self-comparison rows give each file's exact covered bp (union of a file with
# itself is its own size).
file_sizes <- bedtools_raw %>%
  filter(file1 == file2) %>%
  transmute(file = file1, size_bp = as.numeric(union))

ratio_stats <- cross %>%
  filter(estimator == EST_RE) %>%
  left_join(file_sizes, by = c(".a" = "file")) %>%
  rename(size_a = size_bp) %>%
  left_join(file_sizes, by = c(".b" = "file")) %>%
  rename(size_b = size_bp) %>%
  filter(!is.na(size_a), !is.na(size_b)) %>%
  mutate(
    log_ratio = abs(log(size_a / size_b)),
    # Three residual definitions, all three reported because the prose in
    # docs/ quotes values from more than one and does not say which is which.
    # `raw` is the bare gap and is dominated by the floor itself, so it mostly
    # reports the floor's dependence on J; `constrained` and `ols` remove a
    # fitted c + (1 - c) J model first and differ only in how c was fitted.
    resid_raw = hammock_jaccard - bedtools_jaccard,
    resid_constrained = hammock_jaccard -
      (bedtools_jaccard + (1 - bedtools_jaccard) *
         floor_stats$c_constrained[match(precision, floor_stats$precision)]),
    resid_ols = hammock_jaccard -
      (bedtools_jaccard + (1 - bedtools_jaccard) *
         floor_stats$c_ols_intercept[match(precision, floor_stats$precision)]),
    resid_pointwise = hammock_jaccard -
      (bedtools_jaccard + (1 - bedtools_jaccard) *
         floor_stats$c_pointwise[match(precision, floor_stats$precision)])
  ) %>%
  group_by(precision) %>%
  summarise(
    n = n(),
    max_log_ratio = max(log_ratio),
    max_size_ratio = exp(max(log_ratio)),
    cor_raw_logratio = cor(resid_raw, log_ratio),
    cor_constrained_logratio = cor(resid_constrained, log_ratio),
    cor_ols_logratio = cor(resid_ols, log_ratio),
    cor_pointwise_logratio = cor(resid_pointwise, log_ratio),
    # The free-slope two-parameter fit, i.e. residuals of lm(j_re ~ j_bt).
    cor_freeslope_logratio = cor(resid(lm(hammock_jaccard ~ bedtools_jaccard)),
                                 log_ratio),
    .groups = "drop"
  )

message("\nRegister-equality residual versus |log(|A|/|B|)|:")
print(as.data.frame(ratio_stats), digits = 5)

# --- leave-one-file-out jackknife -------------------------------------------
# The 17,955 comparisons come from 190 pairs over 20 files, each file appearing
# in 19 pairs, so they are far from independent and a binomial SE on the
# comparison count understates the uncertainty by roughly an order of magnitude.
# Drop one file at a time and recompute. This is the evidence behind
# docs/estimator-analysis-findings.md section 4.
jackknife <- function(df) {
  files <- sort(unique(c(df$.a, df$.b)))
  reps <- lapply(files, function(f) {
    sub <- df %>% filter(.a != f, .b != f)
    rank_agreement(sub$bedtools_jaccard, sub$hammock_jaccard)
  })
  reps <- bind_rows(reps)
  n <- length(files)
  # Jackknife SE: sqrt((n-1)/n * sum (theta_i - mean)^2).
  se <- function(x) sqrt((n - 1) / n * sum((x - mean(x))^2))
  tibble(
    n_files = n,
    tau_b_se = se(reps$tau_b),
    tau_b_min = min(reps$tau_b), tau_b_max = max(reps$tau_b),
    discordant_pct_se = se(reps$discordant_pct),
    discordant_pct_min = min(reps$discordant_pct),
    discordant_pct_max = max(reps$discordant_pct)
  )
}

jack_stats <- cross %>%
  filter(precision == REFERENCE_PRECISION) %>%
  group_by(estimator) %>%
  group_modify(~ jackknife(.x)) %>%
  ungroup()

message(sprintf(
  "\nLeave-one-file-out jackknife at p = %d (20 files, each in 19 pairs):",
  REFERENCE_PRECISION))
print(as.data.frame(jack_stats), digits = 5)

# `stats` itself stays as the plotting frame; only the CSV carries the joins.
stats_csv <- file.path(script_dir, "interval_accuracy_stats.csv")
write_csv(
  stats %>%
    left_join(rank_stats, by = c("precision", "estimator")) %>%
    left_join(floor_stats, by = "precision") %>%
    left_join(ratio_stats %>% select(-n), by = "precision") %>%
    left_join(jack_stats, by = "estimator"),
  stats_csv
)
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
# Keep the pad small and let the scales expand: with expand = FALSE the panel
# edge lands on the data limit, and the first break (0.10 here, ~0.5% of the
# range inside) draws a light gridline a few pixels from the dark axis line,
# which reads as a doubled axis.
main_pad <- diff(main_range) * 0.02
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
  coord_equal(xlim = main_limits, ylim = main_limits) +
  scale_x_continuous(
    labels = label_number(accuracy = 0.02),
    expand = expansion(mult = 0.04)
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.02),
    expand = expansion(mult = 0.04)
  ) +
  labs(
    x = "BEDTools exact base-pair Jaccard",
    y = "hammock Jaccard estimate"
  ) +
  theme_paper(base_size = 11.5) +
  theme(legend.position = "none")

# ---------------------------------------------------------------------------
# Panel B: precision/runtime frontier
# ---------------------------------------------------------------------------

read_precision_frontier <- function(path) {
  raw <- read_csv(path, show_col_types = FALSE)
  required <- c(
    "tool", "precision", "threads", "wall_time", "jaccard_n_pairs",
    "jaccard_ie_mae_vs_bt"
  )
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0) {
    stop(
      basename(path), " lacks columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  bedtools_rows <- raw %>% filter(tool == "bedtools")
  if (nrow(bedtools_rows) == 0) {
    stop("No BEDTools rows in ", basename(path), call. = FALSE)
  }
  bedtools_wall <- median(bedtools_rows$wall_time)

  frontier <- raw %>%
    filter(tool != "bedtools", !is.na(precision)) %>%
    group_by(precision) %>%
    summarise(
      wall = median(wall_time),
      wall_min = min(wall_time),
      wall_max = max(wall_time),
      mae_ie = median(jaccard_ie_mae_vs_bt),
      n_pairs = median(jaccard_n_pairs),
      threads = first(threads),
      .groups = "drop"
    ) %>%
    mutate(
      speedup = bedtools_wall / wall,
      speedup_low = bedtools_wall / wall_max,
      speedup_high = bedtools_wall / wall_min
    )

  pair_counts <- unique(frontier$n_pairs)
  if (length(pair_counts) != 1 || pair_counts != 190) {
    stop(
      "Expected 190 unique off-diagonal pairs in ", basename(path),
      "; got ", paste(pair_counts, collapse = ", "), call. = FALSE
    )
  }
  frontier
}

frontier <- read_precision_frontier(precision_frontier_csv)
frontier_labels <- frontier %>% mutate(label = paste0("p=", precision))
default_precision <- frontier %>% filter(precision == 18)

x_low <- min(frontier$mae_ie) / 1.45
x_high <- max(frontier$mae_ie) * 1.75
y_low <- min(frontier$speedup) / 1.30
y_high <- max(frontier$speedup) * 1.12

frontier_figure <- ggplot(frontier, aes(mae_ie, speedup)) +
  geom_linerange(
    aes(ymin = speedup_low, ymax = speedup_high),
    color = COL_IE, linewidth = 0.5, alpha = 0.8
  ) +
  geom_path(color = COL_IE, linewidth = 0.6, alpha = 0.8) +
  geom_point(color = COL_IE, size = 2.2) +
  geom_point(
    data = default_precision, shape = 21, size = 5, stroke = 0.9,
    color = "#B8420F", fill = NA
  ) +
  geom_text(
    data = frontier_labels, aes(label = label), size = 3.2,
    hjust = -0.35, vjust = -0.55, color = COL_TEXT
  ) +
  scale_x_log10(
    breaks = breaks_log(n = 6), labels = label_scientific(digits = 2)
  ) +
  scale_y_log10(
    breaks = breaks_log(n = 6),
    labels = label_number(accuracy = 0.01, drop0trailing = TRUE)
  ) +
  labs(
    x = expression(
      "Mean absolute error of " * italic(J)[IE] * " vs exact BEDTools  (log)"
    ),
    y = "Speedup vs BEDTools  (log)"
  ) +
  coord_cartesian(
    xlim = c(x_low, x_high), ylim = c(y_low, y_high), expand = FALSE
  ) +
  theme_paper(base_size = 11.5) +
  theme(legend.position = "none")

combined_figure <- main_figure + frontier_figure +
  plot_layout(widths = c(1, 1.12)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 15, color = COL_TEXT)
    )
  )

CairoPNG(
  filename = main_png,
  width = 13.2,
  height = 6.2,
  units = "in",
  res = 300,
  bg = "white"
)
print(combined_figure)
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
both_pad <- diff(both_range) * 0.02  # see the note on main_pad above
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
  coord_equal(xlim = both_limits, ylim = both_limits) +
  scale_x_continuous(
    labels = label_number(accuracy = 0.05),
    expand = expansion(mult = 0.04)
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.05),
    expand = expansion(mult = 0.04)
  ) +
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
      "Inclusion–exclusion estimates BEDTools set Jaccard,",
      "while register equality is a distinct compatibility statistic",
      sep = "\n"
    ),
    theme = theme(
      plot.title = element_text(
        family = base_family, face = "bold", size = 13.2,
        color = COL_TEXT, lineheight = 1.15, margin = margin(b = 10)
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
