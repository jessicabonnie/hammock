#!/usr/bin/env Rscript

# Figure 5 — Sequence sketches group samples by tissue rather than by reference
# Creates a publication-ready two-panel PNG using CairoPNG.

required_packages <- c(
  "dplyr", "readr", "stringr", "ggplot2", "scales", "ggdendro", "patchwork",
  "Cairo"
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
  library(stringr)
  library(ggplot2)
  library(scales)
  library(ggdendro)
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
peaks_csv <- file.path(data_dir, "exp_a_broad_k10_w10.csv")
meta_tsv <- file.path(data_dir, "exp_a_metadata.tsv")

SIM_COL <- "jaccard_similarity"
KW_LABEL <- "k = 10, w = 10"
PEAK_LABEL <- "broad peaks"

argv <- commandArgs(trailingOnly = TRUE)
out_png <- if (length(argv) >= 1) {
  normalizePath(argv[1], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "cross_reference_identity.png")
}
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

for (path in c(peaks_csv, meta_tsv)) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
}

COL_BEDTOOLS <- "#46515C"
COL_HAMMOCK <- "#007C83"
COL_COMPARE <- "#D28B35"
COL_GRID <- "#D9DEE3"
COL_TEXT <- "#20262D"
base_family <- "sans"
TISSUE_PAL <- c(heart = "#c1272d", liver = "#7f3f00", lung = "#4575b4")

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

meta <- read_tsv(meta_tsv, show_col_types = FALSE) %>%
  mutate(key = paste(sample_id, ref, sep = "__"))

raw <- read_csv(peaks_csv, show_col_types = FALSE)
if (!SIM_COL %in% names(raw)) {
  stop(
    "Similarity column '", SIM_COL, "' not present in ", basename(peaks_csv),
    call. = FALSE
  )
}

pairs <- raw %>%
  mutate(
    key_a = paste(
      str_remove(basename(file1), "\\.fa$"),
      basename(dirname(file1)),
      sep = "__"
    ),
    key_b = paste(
      str_remove(basename(file2), "\\.fa$"),
      basename(dirname(file2)),
      sep = "__"
    ),
    similarity = .data[[SIM_COL]]
  ) %>%
  select(key_a, key_b, similarity)

keys <- sort(unique(c(pairs$key_a, pairs$key_b)))
missing_keys <- setdiff(keys, meta$key)
if (length(missing_keys) > 0) {
  stop(
    "Metadata is missing entries for: ", paste(missing_keys, collapse = ", "),
    call. = FALSE
  )
}

# Panel A: UPGMA over 1 - Jaccard ---------------------------------------------
n <- length(keys)
sim <- matrix(NA_real_, n, n, dimnames = list(keys, keys))
for (i in seq_len(nrow(pairs))) {
  sim[pairs$key_a[i], pairs$key_b[i]] <- pairs$similarity[i]
}
sim[is.na(sim)] <- t(sim)[is.na(sim)]
if (anyNA(sim)) stop("Similarity matrix is incomplete.", call. = FALSE)
diag(sim) <- 1

hc <- hclust(as.dist(1 - sim), method = "average")
dd <- dendro_data(as.dendrogram(hc), type = "rectangle")
s <- segment(dd)
segs <- data.frame(x = s$y, y = s$x, xend = s$yend, yend = s$xend)

leaf_labels <- label(dd) %>%
  left_join(meta %>% select(key, tissue, ref), by = c("label" = "key")) %>%
  mutate(pretty = paste(tissue, ref, sep = " · "), plot_y = x)

x_max <- max(segs$x, segs$xend)

panel_a <- ggplot() +
  geom_segment(
    data = segs,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5,
    color = COL_TEXT
  ) +
  geom_text(
    data = leaf_labels,
    aes(x = 0, y = plot_y, label = pretty, color = tissue),
    hjust = 0,
    nudge_x = x_max * 0.03,
    size = 3.2
  ) +
  scale_color_manual(values = TISSUE_PAL, guide = "none") +
  scale_x_reverse(
    limits = c(x_max * 1.05, -x_max * 0.62),
    labels = label_number(accuracy = 0.01)
  ) +
  labs(
    title = sprintf(
      "A  Tissue identity dominates reference choice\n(%s; %s)",
      KW_LABEL,
      PEAK_LABEL
    ),
    x = "1 − Jaccard (UPGMA)",
    y = NULL
  ) +
  theme_paper() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )

# Panel B: three comparison classes over unique unordered pairs ---------------
# The 36 unique pairs split into three classes rather than two. Splitting the
# old "different tissue, any reference" class separates the two situations it
# conflated, and makes the tissue-vs-reference claim a balanced 9-vs-9 contrast
# (same tissue with the reference changed vs. same reference with the tissue
# changed) instead of a 9-vs-27 comparison.
LAB_ST_DR <- "Same tissue\nDifferent reference"
LAB_DT_SR <- "Different tissue\nSame reference"
LAB_DT_DR <- "Different tissue\nDifferent reference"
CLASS_LEVELS <- c(LAB_ST_DR, LAB_DT_SR, LAB_DT_DR)

classified <- pairs %>%
  filter(key_a != key_b) %>%
  mutate(
    .a = pmin(key_a, key_b),
    .b = pmax(key_a, key_b)
  ) %>%
  distinct(.a, .b, .keep_all = TRUE) %>%
  left_join(
    meta %>% select(key, tissue_a = tissue, ref_a = ref),
    by = c(".a" = "key")
  ) %>%
  left_join(
    meta %>% select(key, tissue_b = tissue, ref_b = ref),
    by = c(".b" = "key")
  ) %>%
  mutate(
    same_tissue = tissue_a == tissue_b,
    same_ref = ref_a == ref_b,
    pair_class = factor(
      case_when(
        same_tissue & !same_ref ~ LAB_ST_DR,
        !same_tissue & same_ref ~ LAB_DT_SR,
        TRUE ~ LAB_DT_DR
      ),
      levels = CLASS_LEVELS
    ),
    x_pos = as.numeric(pair_class)
  )

expected_pairs <- n * (n - 1) / 2
if (nrow(classified) != expected_pairs) {
  stop(
    "Expected ", expected_pairs, " unique pairs, got ", nrow(classified),
    call. = FALSE
  )
}
# Each tissue appears exactly once per reference, so no pair can share both.
if (any(classified$same_tissue & classified$same_ref)) {
  stop("Design is not one sample per tissue per reference.", call. = FALSE)
}

by_class <- split(classified$similarity, classified$pair_class)
st_dr <- by_class[[LAB_ST_DR]]
dt_sr <- by_class[[LAB_DT_SR]]
dt_dr <- by_class[[LAB_DT_DR]]

auc_gt <- function(x, y) {
  mean(outer(x, y, ">") + 0.5 * outer(x, y, "=="))
}

# Exact constrained label permutation -----------------------------------------
# The nine peak sets form a 3 x 3 tissue-by-reference design, so the 36 pairs
# are not 36 independent observations and a rank-sum test on them is not valid.
# Instead permute the tissue labels *within each reference*, which preserves the
# reference composition (every reference still contributes one sample per label)
# and asks only whether tissue identity is related to sequence similarity.
#
# There are (3!)^3 = 216 such labelings. Relabeling the three tissue names
# globally leaves every comparison class unchanged, so those 216 labelings
# induce only 216 / 3! = 36 distinct classifications and the smallest attainable
# p-value is 6 / 216 = 1 / 36 ~ 0.028.
#
# Note that the "different tissue, same reference" class is invariant: within a
# reference the three permuted labels are always distinct. The null therefore
# varies only through which cross-reference pairs count as same-tissue, which is
# exactly the comparison of interest.
all_perms <- function(x) {
  if (length(x) <= 1) return(list(x))
  do.call(c, lapply(seq_along(x), function(i) {
    lapply(all_perms(x[-i]), function(rest) c(x[i], rest))
  }))
}

refs <- sort(unique(meta$ref))
ref_units <- lapply(refs, function(r) {
  meta %>% filter(ref == r) %>% arrange(sample_id) %>% select(key, tissue)
})
n_per_ref <- unique(vapply(ref_units, nrow, integer(1)))
if (length(n_per_ref) != 1) {
  stop("References contribute differing numbers of samples.", call. = FALSE)
}

within_perms <- all_perms(seq_len(n_per_ref))
perm_grid <- as.matrix(
  expand.grid(rep(list(seq_along(within_perms)), length(refs)))
)

# AUC is the test statistic: it is a rank statistic of the full class ordering,
# so it does not hinge on a single order statistic the way a median difference
# does. The median difference is reported as a descriptive effect size only.
stat_from_labels <- function(labels) {
  same <- labels[classified$.a] == labels[classified$.b]
  g_st_dr <- classified$similarity[same & !classified$same_ref]
  g_dt_sr <- classified$similarity[!same & classified$same_ref]
  c(
    auc = auc_gt(g_st_dr, g_dt_sr),
    delta = median(g_st_dr) - median(g_dt_sr)
  )
}

labels_for <- function(perm_row) {
  out <- character(0)
  for (i in seq_along(ref_units)) {
    unit <- ref_units[[i]]
    out[unit$key] <- unit$tissue[within_perms[[perm_row[i]]]]
  }
  out
}

observed_labels <- setNames(meta$tissue, meta$key)
observed <- stat_from_labels(observed_labels)
null_stats <- t(apply(perm_grid, 1, function(r) stat_from_labels(labels_for(r))))

n_perms <- nrow(null_stats)
perm_p_auc <- mean(null_stats[, "auc"] >= observed[["auc"]])
perm_p_delta <- mean(null_stats[, "delta"] >= observed[["delta"]])
min_attainable_p <- (factorial(n_per_ref)) / n_perms

auc <- observed[["auc"]]
delta <- observed[["delta"]]

stats <- tibble::tibble(
  metric = SIM_COL,
  n_st_dr = length(st_dr),
  n_dt_sr = length(dt_sr),
  n_dt_dr = length(dt_dr),
  median_st_dr = median(st_dr),
  median_dt_sr = median(dt_sr),
  median_dt_dr = median(dt_dr),
  delta_balanced = delta,
  auc_balanced = auc,
  n_permutations = n_perms,
  perm_p_auc = perm_p_auc,
  perm_p_delta = perm_p_delta,
  min_attainable_p = min_attainable_p
)
message("Comparison classes over ", nrow(classified), " unique pairs:")
print(as.data.frame(stats), digits = 5)
message(
  "Balanced contrast: ", LAB_ST_DR %>% str_replace("\n", ", "), " (n = ",
  length(st_dr), ") vs ", LAB_DT_SR %>% str_replace("\n", ", "), " (n = ",
  length(dt_sr), ")"
)
message(
  "Exact constrained permutation over ", n_perms, " labelings; ",
  "largest AUC under the null = ",
  signif(max(null_stats[null_stats[, "auc"] < auc, "auc"]), 4),
  "; largest median difference under the null = ",
  signif(max(null_stats[null_stats[, "delta"] < delta, "delta"]), 4)
)

annotation_b <- sprintf(
  "AUC = %.2f\nΔ median = %+.3f\nExact permutation p = %.3f",
  auc,
  delta,
  perm_p_auc
)

other_max <- max(c(dt_sr, dt_dr))
separation <- min(st_dr) - other_max
gap_layer <- if (separation > 0) {
  list(
    annotate(
      "rect",
      xmin = 0.4,
      xmax = length(CLASS_LEVELS) + 0.6,
      ymin = other_max,
      ymax = min(st_dr),
      fill = COL_COMPARE,
      alpha = 0.10
    ),
    annotate(
      "text",
      x = length(CLASS_LEVELS) + 0.55,
      y = mean(c(other_max, min(st_dr))),
      label = "no overlap",
      hjust = 1,
      vjust = 0.5,
      size = 2.9,
      color = COL_COMPARE
    )
  )
} else {
  list()
}

# Reserve space above the observations for the statistics box so it cannot
# obscure points in either comparison class.
y_range <- range(classified$similarity, finite = TRUE)
y_pad <- max(diff(y_range) * 0.32, 0.015)
y_upper <- y_range[2] + y_pad

set.seed(1)
panel_b <- ggplot(classified, aes(x = x_pos, y = similarity)) +
  gap_layer +
  geom_point(
    aes(color = pair_class),
    position = position_jitter(width = 0.14, height = 0),
    size = 2.1,
    alpha = 0.75
  ) +
  stat_summary(
    aes(group = pair_class),
    fun = median,
    geom = "crossbar",
    width = 0.42,
    linewidth = 0.4,
    color = COL_TEXT,
    fatten = 0
  ) +
  # Bracket marking the two classes the permutation test contrasts.
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = y_upper - y_pad * 0.86,
    yend = y_upper - y_pad * 0.86,
    linewidth = 0.35,
    color = COL_TEXT
  ) +
  annotate(
    "label",
    x = 1.5,
    y = y_upper,
    label = annotation_b,
    hjust = 0.5,
    vjust = 1,
    size = 2.85,
    lineheight = 1.04,
    linewidth = 0,
    fill = alpha("white", 0.92),
    color = COL_TEXT
  ) +
  scale_color_manual(
    values = setNames(
      c(COL_HAMMOCK, COL_BEDTOOLS, "#98A2AC"),
      CLASS_LEVELS
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(classified$pair_class)),
    labels = levels(classified$pair_class)
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.01),
    limits = c(y_range[1], y_upper),
    expand = expansion(mult = c(0.03, 0.02))
  ) +
  coord_cartesian(xlim = c(0.4, length(CLASS_LEVELS) + 0.6)) +
  labs(
    title = sprintf(
      "B  Same-tissue pairs rank above every different-tissue pair\n(n = %d, %d, and %d unique pairs)",
      length(st_dr),
      length(dt_sr),
      length(dt_dr)
    ),
    x = NULL,
    y = "Sequence-mode Jaccard"
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(size = 8.4, lineheight = 0.95),
    legend.position = "none"
  )

figure <- panel_a + panel_b +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Sequence sketches group samples by tissue, not by reference genome",
    theme = theme(
      plot.title = element_text(
        family = base_family,
        face = "bold",
        size = 15,
        color = COL_TEXT,
        margin = margin(b = 10)
      ),
      plot.margin = margin(12, 16, 12, 12)
    )
  )

CairoPNG(
  filename = out_png,
  width = 12.4,
  height = 5.9,
  units = "in",
  res = 300,
  bg = "white"
)
print(figure)
dev.off()

message("Wrote: ", out_png)
