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

# Panel B: comparison classes over unique unordered pairs ---------------------
XREF_LABEL <- "Same tissue,\ndifferent reference"
DIFF_LABEL <- "Different tissue,\nany reference"

classified <- pairs %>%
  filter(key_a != key_b) %>%
  mutate(
    .a = pmin(key_a, key_b),
    .b = pmax(key_a, key_b)
  ) %>%
  distinct(.a, .b, .keep_all = TRUE) %>%
  left_join(meta %>% select(key, tissue_a = tissue), by = c(".a" = "key")) %>%
  left_join(meta %>% select(key, tissue_b = tissue), by = c(".b" = "key")) %>%
  mutate(
    pair_class = factor(
      if_else(tissue_a == tissue_b, XREF_LABEL, DIFF_LABEL),
      levels = c(XREF_LABEL, DIFF_LABEL)
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

xref <- classified %>%
  filter(pair_class == XREF_LABEL) %>%
  pull(similarity)
diff_tissue <- classified %>%
  filter(pair_class == DIFF_LABEL) %>%
  pull(similarity)

auc <- mean(
  outer(xref, diff_tissue, ">") + 0.5 * outer(xref, diff_tissue, "==")
)
test_result <- wilcox.test(xref, diff_tissue, alternative = "greater")

stats <- tibble::tibble(
  metric = SIM_COL,
  n_xref = length(xref),
  n_diff = length(diff_tissue),
  median_xref = median(xref),
  median_diff = median(diff_tissue),
  delta = median(xref) - median(diff_tissue),
  auc = auc,
  wilcoxon_p = test_result$p.value
)
message("Comparison classes over ", nrow(classified), " unique pairs:")
print(as.data.frame(stats), digits = 5)

annotation_b <- sprintf(
  "AUC = %.2f\nMann–Whitney p = %s\nΔ median = %+.3f",
  auc,
  format(signif(test_result$p.value, 2), scientific = TRUE),
  stats$delta
)

separation <- min(xref) - max(diff_tissue)
gap_layer <- if (separation > 0) {
  list(
    annotate(
      "rect",
      xmin = 0.4,
      xmax = 2.6,
      ymin = max(diff_tissue),
      ymax = min(xref),
      fill = COL_COMPARE,
      alpha = 0.10
    ),
    annotate(
      "text",
      x = 2.55,
      y = mean(c(max(diff_tissue), min(xref))),
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
      c(COL_HAMMOCK, COL_BEDTOOLS),
      c(XREF_LABEL, DIFF_LABEL)
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
  coord_cartesian(xlim = c(0.4, 2.6)) +
  labs(
    title = sprintf(
      "B  Same-tissue pairs rank above different-tissue pairs\n(n = %d and %d unique pairs)",
      length(xref),
      length(diff_tissue)
    ),
    x = NULL,
    y = "Sequence-mode Jaccard"
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(size = 9, lineheight = 0.95),
    legend.position = "none"
  )

figure <- panel_a + panel_b +
  plot_layout(widths = c(1.25, 1)) +
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
  width = 11.6,
  height = 5.9,
  units = "in",
  res = 300,
  bg = "white"
)
print(figure)
dev.off()

message("Wrote: ", out_png)
