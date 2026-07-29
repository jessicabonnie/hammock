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
COL_GRID <- "#D9DEE3"
base_family <- "sans"
TISSUE_PAL <- c(heart = "#c1272d", liver = "#7f3f00", lung = "#4575b4")

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
  mutate(pretty = paste(tissue, ref, sep = " · "), plot_y = x)

x_max <- max(segs$x, segs$xend)

panel <- ggplot() +
  geom_segment(
    data = segs,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.6,
    color = COL_TEXT
  ) +
  geom_text(
    data = leaf_labels,
    aes(x = 0, y = plot_y, label = pretty, color = tissue),
    hjust = 0,
    nudge_x = x_max * 0.03,
    size = 4.0
  ) +
  scale_color_manual(values = TISSUE_PAL, guide = "none") +
  scale_x_reverse(
    limits = c(x_max * 1.05, -x_max * 0.72),
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Sequence sketches group samples by tissue across references",
    x = "1 − Jaccard (UPGMA)",
    y = NULL
  ) +
  theme_classic(base_size = 12, base_family = base_family) +
  theme(
    plot.title = element_text(
      face = "bold", color = COL_TEXT, margin = margin(b = 10)
    ),
    axis.title = element_text(color = COL_TEXT),
    axis.text = element_text(color = COL_TEXT),
    axis.line = element_line(color = "#6B747D", linewidth = 0.4),
    axis.ticks = element_line(color = "#6B747D", linewidth = 0.4),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.x = element_line(color = COL_GRID, linewidth = 0.35),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(14, 18, 14, 14)
  )

CairoPNG(
  filename = out_png,
  width = 8.8,
  height = 6.4,
  units = "in",
  res = 300,
  bg = "white"
)
print(panel)
dev.off()

message("Wrote: ", out_png)
