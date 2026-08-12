#!/usr/bin/env Rscript

# Figure 5 — Sequence sketches group samples by tissue rather than by reference
# Single-panel dendrogram rendered to PNG using CairoPNG.

required_packages <- c(
  "dplyr", "readr", "stringr", "ggplot2", "scales", "ggdendro", "png", "Cairo"
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
  library(png)
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

argv <- commandArgs(trailingOnly = TRUE)
out_png <- if (length(argv) >= 1) {
  normalizePath(argv[1], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "cross_reference_identity.png")
}
# argv[2] selects the similarity column explicitly. `jaccard_similarity_ie` is
# accepted even though this CSV predates that column: it carries
# containment_AB and containment_BA, from which it is exactly recoverable
# (see below). When not given, resolved below (after raw is loaded) preferring
# reg_eq_similarity and falling back to jaccard_similarity for archived
# pre-Step-1 CSVs -- see docs/seed-jaccard-reg-eq-rename.md Step 2.
SIM_COL_ARG <- if (length(argv) >= 2 && nzchar(argv[2])) argv[2] else NA_character_
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

for (path in c(peaks_csv, meta_tsv)) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
}

COL_TEXT <- "#20262D"
base_family <- "sans"
# Okabe-Ito hues (colour-vision-deficiency safe), darkened where needed so the
# leaf labels clear ~4.5:1 contrast on white as body text. These same values
# recolor the organ outlines, so tissue labels and icons have one color source.
TISSUE_PAL <- c(heart = "#B8500B", liver = "#00745B", lung = "#0072B2")

organ_dir <- file.path(repo_root, "paper", "organs")
ORGAN_FILES <- c(
  heart = file.path(organ_dir, "heart_outline.png"),
  liver = file.path(organ_dir, "liver_outline.png"),
  lung = file.path(organ_dir, "lung_outline.png"),
  stomach = file.path(organ_dir, "stomach_outline.png"),
  kidney = file.path(organ_dir, "kidney_outline.png")
)

meta <- read_tsv(meta_tsv, show_col_types = FALSE) %>%
  mutate(key = paste(sample_id, ref, sep = "__"))

jaccard_ie_from_containments <- function(c_ab, c_ba) {
  # Set-Jaccard from the two directional containments:
  #   1 / (1/C_AB + 1/C_BA - 1) == |A n B| / (|A| + |B| - |A n B|).
  # Mirrors python/hammock/runner.py::_jaccard_ie_from_containments, clamp
  # included, and matches experiments/maurano_dhs_validation/analyze.R:166-171.
  cab <- pmin(as.numeric(c_ab), 1)
  cba <- pmin(as.numeric(c_ba), 1)
  ifelse(cab <= 0 | cba <= 0, 0, 1 / (1 / cab + 1 / cba - 1))
}

# Oracle values from tests/test_jaccard_ie.py, which pins the canonical helper.
local({
  probe <- jaccard_ie_from_containments(
    c(0.5, 1.0, 1.0000000000050957, 0.0, 0.5),
    c(0.5, 1.0, 1.0,                0.5, 0.0)
  )
  stopifnot(
    isTRUE(all.equal(probe[1], 1 / 3)),
    probe[2] == 1, probe[3] == 1, probe[4] == 0, probe[5] == 0
  )
})

raw <- read_csv(peaks_csv, show_col_types = FALSE)
if (!("jaccard_similarity_ie" %in% names(raw)) &&
    all(c("containment_AB", "containment_BA") %in% names(raw))) {
  raw[["jaccard_similarity_ie"]] <-
    jaccard_ie_from_containments(raw$containment_AB, raw$containment_BA)
}
resolve_sim_col <- function(explicit, df) {
  if (!is.na(explicit)) return(explicit)
  if ("reg_eq_similarity" %in% names(df)) return("reg_eq_similarity")
  message(
    "reg_eq_similarity not found in input CSV; falling back to jaccard_similarity"
  )
  "jaccard_similarity"
}
SIM_COL <- resolve_sim_col(SIM_COL_ARG, raw)
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
  mutate(pretty = paste(str_to_title(tissue), ref, sep = " · "), plot_y = x)

x_max <- max(segs$x, segs$xend)

# Panel width = the tree (left of zero) plus a gutter that holds the leaf
# labels and the organ glyphs. Ticks and the axis line stay inside the tree
# half: the gutter is not distance, so it must not pick up (negative) breaks
# or carry the axis rule.
TREE_PAD <- 1.04
# Wider label gutter so the larger organ glyphs sit clear of the leaf text.
LABEL_PAD <- 0.78
tree_frac <- TREE_PAD / (TREE_PAD + LABEL_PAD)
ICON_X <- -x_max * 0.52

x_breaks <- pretty(c(0, x_max), n = 4)
x_breaks <- x_breaks[x_breaks >= 0 & x_breaks <= x_max]

# --- Organ glyphs ----------------------------------------------------------
# Publication-ready transparent PNG outlines live in paper/organs/. The files
# themselves are monochrome masks; recoloring happens here from TISSUE_PAL so
# changing a tissue color automatically updates both its labels and its icon.
# Stomach and kidney assets are included for reuse by later tissue figures even
# though this three-tissue cross-reference panel uses heart, liver, and lung.
ICON_SIZE <- grid::unit(0.72, "in")
# Dilate the alpha mask by this many pixels (in source-image space) so thin
# outline strokes read as bold at figure size without editing the PNGs.
ICON_STROKE_PX <- 5L

# Disk max-filter on a 2-D alpha matrix. Thickens anti-aliased outline strokes
# before recoloring.
dilate_alpha <- function(mat, radius) {
  if (radius <= 0L) {
    return(mat)
  }
  nr <- nrow(mat)
  nc <- ncol(mat)
  pad <- matrix(0, nr + 2L * radius, nc + 2L * radius)
  pad[(radius + 1L):(radius + nr), (radius + 1L):(radius + nc)] <- mat
  out <- matrix(0, nr, nc)
  for (dy in -radius:radius) {
    for (dx in -radius:radius) {
      if (dx * dx + dy * dy > radius * radius) {
        next
      }
      ys <- (radius + 1L + dy):(radius + nr + dy)
      xs <- (radius + 1L + dx):(radius + nc + dx)
      out <- pmax(out, pad[ys, xs])
    }
  }
  out
}

recolor_icon <- function(path, color) {
  if (!file.exists(path)) {
    stop("Organ icon not found: ", path, call. = FALSE)
  }

  img <- png::readPNG(path)
  if (length(dim(img)) != 3 || dim(img)[3] < 4) {
    stop("Organ icon must be an RGBA PNG with transparency: ", path, call. = FALSE)
  }

  alpha <- dilate_alpha(img[, , 4], ICON_STROKE_PX)
  rgb <- grDevices::col2rgb(color) / 255
  img[, , 1] <- rgb[1]
  img[, , 2] <- rgb[2]
  img[, , 3] <- rgb[3]
  img[, , 4] <- alpha
  img
}

organ_grob <- function(tissue, color) {
  path <- unname(ORGAN_FILES[tissue])
  if (length(path) != 1 || is.na(path)) {
    stop("No organ PNG defined for tissue: ", tissue, call. = FALSE)
  }

  # The outer grobTree absorbs the viewport annotation_custom() assigns to the
  # top-level grob; the raster child retains a fixed physical size and is not
  # distorted by the dendrogram panel aspect ratio.
  grid::grobTree(
    grid::rasterGrob(
      recolor_icon(path, color),
      width = ICON_SIZE,
      height = ICON_SIZE,
      interpolate = TRUE
    )
  )
}

organ_positions <- leaf_labels %>%
  group_by(tissue) %>%
  summarise(y = mean(plot_y), .groups = "drop")

organ_layers <- lapply(seq_len(nrow(organ_positions)), function(i) {
  tissue <- organ_positions$tissue[i]
  annotation_custom(
    organ_grob(tissue, unname(TISSUE_PAL[tissue])),
    # Extents are in data space: the coord maps them through scale_x_reverse,
    # so ICON_X goes in as-is (negating it lands the glyph inside the tree).
    xmin = ICON_X, xmax = ICON_X,
    ymin = organ_positions$y[i], ymax = organ_positions$y[i]
  )
})

panel <- ggplot() +
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
    size = 3.2,
    fontface = "bold"
  ) +
  organ_layers +
  annotate(
    "segment",
    x = x_max * TREE_PAD, xend = 0, y = -Inf, yend = -Inf,
    color = "#6B747D", linewidth = 0.4
  ) +
  scale_color_manual(values = TISSUE_PAL, guide = "none") +
  scale_x_reverse(
    limits = c(x_max * TREE_PAD, -x_max * LABEL_PAD),
    breaks = x_breaks,
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Sequence sketches group samples by tissue across references",
    # Named only on the non-default column, so the manuscript figure is
    # unchanged. It goes in the subtitle rather than the title because the
    # title is already near the panel width and appending to it clips.
    subtitle = if (is.na(SIM_COL_ARG)) NULL else {
      sprintf("similarity = %s", SIM_COL)
    },
    x = "1 − Jaccard (UPGMA)",
    y = NULL
  ) +
  theme_classic(base_size = 11, base_family = base_family) +
  theme(
    plot.title = element_text(
      size = 12, face = "bold", color = COL_TEXT, margin = margin(b = 8)
    ),
    axis.title = element_text(color = COL_TEXT),
    # Centre the x title under the tree, not under the label gutter.
    axis.title.x = element_text(
      color = COL_TEXT, hjust = tree_frac / 2, margin = margin(t = 4)
    ),
    axis.text = element_text(color = COL_TEXT),
    axis.ticks = element_line(color = "#6B747D", linewidth = 0.4),
    # The x rule is drawn as an annotate() segment so it stops at zero.
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(12, 14, 10, 12)
  )

CairoPNG(
  filename = out_png,
  width = 7.0,
  height = 3.9,
  units = "in",
  res = 300,
  bg = "white"
)
print(panel)
dev.off()

message("Wrote: ", out_png)
