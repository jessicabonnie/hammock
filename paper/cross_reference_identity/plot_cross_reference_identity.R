#!/usr/bin/env Rscript

# Figure 5 — Sequence sketches group samples by tissue rather than by reference
# Single-panel dendrogram rendered to PNG using CairoPNG.

required_packages <- c(
  "dplyr", "readr", "stringr", "ggplot2", "scales", "ggdendro", "Cairo"
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

COL_TEXT <- "#20262D"
base_family <- "sans"
# Okabe-Ito hues (colour-vision-deficiency safe), darkened where needed so the
# leaf labels clear ~4.5:1 contrast on white as body text.
TISSUE_PAL <- c(heart = "#B8500B", liver = "#00745B", lung = "#0072B2")

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
LABEL_PAD <- 0.56
tree_frac <- TREE_PAD / (TREE_PAD + LABEL_PAD)
ICON_X <- -x_max * 0.39

x_breaks <- pretty(c(0, x_max), n = 4)
x_breaks <- x_breaks[x_breaks >= 0 & x_breaks <= x_max]

# --- Organ glyphs ----------------------------------------------------------
# Hand-drawn silhouettes, defined in a unit box and rendered at a fixed
# absolute size so the panel aspect never distorts them. Each is stroked in
# its tissue colour, matching the leaf labels it sits beside.
ICON_SIZE <- grid::unit(0.44, "in")

outline <- function(x, y, gp) {
  grid::xsplineGrob(x = x, y = y, shape = -0.75, open = FALSE, gp = gp)
}

stroke <- function(x, y, gp) {
  grid::xsplineGrob(x = x, y = y, shape = -0.75, open = TRUE, gp = gp)
}

ORGAN_SHAPES <- list(
  # Anatomical heart: tilted ventricular mass with the apex low-left, an aortic
  # arch over the base, vena cava and pulmonary stubs, interventricular groove.
  heart = function(gp) {
    body_x <- c(0.30, 0.16, 0.14, 0.26, 0.48, 0.70, 0.82, 0.80, 0.62, 0.44)
    body_y <- c(0.06, 0.28, 0.52, 0.68, 0.72, 0.68, 0.52, 0.32, 0.12, 0.04)
    grid::gList(
      outline(body_x, body_y, gp),
      # aortic arch
      stroke(c(0.40, 0.39, 0.50, 0.61, 0.62), c(0.70, 0.88, 0.96, 0.87, 0.72), gp),
      # superior vena cava
      stroke(c(0.71, 0.76, 0.74), c(0.66, 0.80, 0.93), gp),
      # pulmonary trunk
      stroke(c(0.31, 0.26, 0.28), c(0.70, 0.82, 0.93), gp),
      # interventricular groove
      stroke(c(0.47, 0.40, 0.31), c(0.71, 0.42, 0.10), gp)
    )
  },
  liver = function(gp) {
    lx <- c(0.05, 0.13, 0.32, 0.60, 0.85, 0.95, 0.90, 0.70, 0.44, 0.20)
    ly <- c(0.36, 0.60, 0.75, 0.79, 0.69, 0.51, 0.33, 0.21, 0.22, 0.28)
    grid::gList(
      outline(lx, ly, gp),
      # falciform ligament — the notch that reads as "liver"
      grid::linesGrob(c(0.56, 0.52), c(0.78, 0.46), gp = gp)
    )
  },
  lung = function(gp) {
    lobe_x <- c(0.56, 0.72, 0.88, 0.92, 0.82, 0.66, 0.58, 0.55)
    lobe_y <- c(0.70, 0.72, 0.55, 0.32, 0.10, 0.08, 0.26, 0.48)
    grid::gList(
      outline(lobe_x, lobe_y, gp),
      outline(1 - lobe_x, lobe_y, gp),
      grid::linesGrob(c(0.5, 0.5), c(0.96, 0.80), gp = gp),          # trachea
      grid::linesGrob(c(0.33, 0.5, 0.67), c(0.70, 0.80, 0.70), gp = gp)  # bronchi
    )
  }
)

organ_grob <- function(tissue, color) {
  if (is.null(ORGAN_SHAPES[[tissue]])) {
    stop("No organ glyph defined for tissue: ", tissue, call. = FALSE)
  }
  gp <- grid::gpar(
    col = color, fill = NA, lwd = 1.3, lineend = "round", linejoin = "round"
  )
  # The outer grobTree absorbs the viewport annotation_custom() forces onto the
  # top-level grob; the inner gTree keeps the glyph at its absolute size.
  grid::grobTree(grid::gTree(
    children = ORGAN_SHAPES[[tissue]](gp),
    vp = grid::viewport(width = ICON_SIZE, height = ICON_SIZE)
  ))
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
