#!/usr/bin/env Rscript

# Figure S5c — pairwise sample-distance heatmap, companion to Figure 6
#
# Figure 6's dendrogram only shows the true pairwise distance between two
# leaves when they are direct siblings ("cherries"); between any other pair
# of leaves, the height where they visually join is the *cophenetic*
# distance (height of their nearest common ancestor), which can distort the
# real pairwise distance substantially. This figure plots the full 1-J
# pairwise distance matrix directly instead, ordered to match Figure 6's
# leaf order, so every one of the 190 unique pairwise distances is visible
# at once rather than only the 19 that appear as tree merges.
#
# Default analysis matches Figure 6 exactly:
#   similarity = jaccard_similarity_ie (inclusion-exclusion Jaccard)
#   k = 10, w = 30, p = 18; linkage = average. The heatmap is ordered by the
#   dendrogram leaves and outlines the k=10 cut used for the reported ARI.
#
# Usage:
#   Rscript paper/sequence_tissue_clustering/plot_distance_heatmap.R
#
# Optional overrides:
#   Rscript ... <similarity_csv> <tissue_key.tsv> <output.png> [similarity_column] [panel_label] [cluster_strip] [legend] [precision_label] [black_accession_labels] [clusters]

required_packages <- c("dplyr", "readr", "ggplot2", "scales", "Cairo")
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
  library(Cairo)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
  stop("Could not determine the script path. Run with Rscript.", call. = FALSE)
}
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

K <- 10
W <- 30
P <- 18

experiment_dir <- file.path(repo_root, "experiments", "maurano_dhs_validation")
default_csv <- file.path(
  experiment_dir, "results", "raw_d",
  sprintf("hammock_mnmzr_p%d_jaccD_k%d_w%d.csv", P, K, W)
)
default_key <- file.path(experiment_dir, "data", "maurano_filenames_key.tsv")
default_output <- file.path(repo_root, "paper", "figures",
                            "sequence_tissue_distance_heatmap.png")

argv <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(argv) >= 1) argv[1] else default_csv
key_tsv <- if (length(argv) >= 2) argv[2] else default_key
out_png <- if (length(argv) >= 3) argv[3] else default_output
sim_col_arg <- if (length(argv) >= 4 && nzchar(argv[4])) argv[4] else NA_character_
panel_label <- if (length(argv) >= 5 && nzchar(argv[5])) argv[5] else "B"
show_cluster_strip <- if (length(argv) >= 6 && nzchar(argv[6])) {
  tolower(argv[6]) %in% c("1", "true", "yes", "cluster-strip")
} else FALSE
show_legend <- if (length(argv) >= 7 && nzchar(argv[7])) {
  tolower(argv[7]) %in% c("1", "true", "yes", "legend")
} else TRUE
precision_label <- if (length(argv) >= 8 && nzchar(argv[8])) {
  precision_value <- suppressWarnings(as.integer(argv[8]))
  if (is.na(precision_value)) argv[8] else bquote(italic(p) == .(precision_value))
} else NULL
black_accession_labels <- if (length(argv) >= 9 && nzchar(argv[9])) {
  tolower(argv[9]) %in% c("1", "true", "yes", "black")
} else FALSE
cluster_count <- if (length(argv) >= 10 && nzchar(argv[10])) {
  suppressWarnings(as.integer(argv[10]))
} else NA_integer_
if (!is.na(cluster_count) && cluster_count < 2) {
  stop("clusters must be an integer of at least 2.", call. = FALSE)
}

dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
for (path in c(input_csv, key_tsv)) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
}

strip_ext <- function(x) sub("\\.(bed|fa|fasta|fna)$", "", basename(x), ignore.case = TRUE)
short_label <- function(stem) sub("^[^-]*-", "", sub("\\..*$", "", stem))

GROUP_DISPLAY <- c(
  fBrain = "Brain",
  fHeart = "Heart",
  fIntestine_Sm = "Small\nintestine",
  fKidney_renal_cortex_L = "Kidney",
  fLung = "Lung",
  fMuscle_arm = "Arm",
  fMuscle_back = "Back",
  fMuscle_leg = "Leg",
  fMuscle = "Muscle",
  fSkin_fibro_bicep_R = "Skin",
  fStomach = "Stomach"
)

BRACKET_GROUP <- c(
  fMuscle_arm = "Muscle",
  fMuscle_back = "Muscle",
  fMuscle_leg = "Muscle"
)

display_group <- function(x) {
  out <- unname(GROUP_DISPLAY[x])
  fallback <- gsub("_", " ", sub("^f(?=[A-Z])", "", x, perl = TRUE))
  fallback <- paste0(toupper(substr(fallback, 1, 1)), substr(fallback, 2, nchar(fallback)))
  ifelse(is.na(out), fallback, out)
}

jaccard_ie_from_containments <- function(c_ab, c_ba) {
  # Mirrors plot_sequence_tissue_clustering.R / runner._jaccard_ie_from_containments.
  cab <- pmin(as.numeric(c_ab), 1)
  cba <- pmin(as.numeric(c_ba), 1)
  ifelse(cab <= 0 | cba <= 0, 0, 1 / (1 / cab + 1 / cba - 1))
}
add_jaccard_ie <- function(df) {
  if (!("jaccard_similarity_ie" %in% names(df)) &&
      all(c("containment_AB", "containment_BA") %in% names(df))) {
    df[["jaccard_similarity_ie"]] <-
      jaccard_ie_from_containments(df$containment_AB, df$containment_BA)
  }
  df
}

key <- read_tsv(key_tsv, show_col_types = FALSE) %>%
  transmute(stem = strip_ext(File), tissue = Biosample_term_name)
if (identical(cluster_count, 8L)) {
  key <- key %>% mutate(tissue = if_else(grepl("^fMuscle_", tissue), "fMuscle", tissue))
}
tissue_by_stem <- setNames(key$tissue, key$stem)

raw <- add_jaccard_ie(read_csv(input_csv, show_col_types = FALSE))

resolve_sim_col <- function(explicit, df) {
  if (!is.na(explicit)) return(explicit)
  "jaccard_similarity_ie"
}
sim_col <- resolve_sim_col(sim_col_arg, raw)

required_cols <- c("file1", "file2", sim_col)
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("Input lacks columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

stems <- sort(unique(c(strip_ext(raw$file1), strip_ext(raw$file2))))
missing_stems <- setdiff(stems, key$stem)
if (length(missing_stems) > 0) {
  stop("Tissue key is missing: ", paste(missing_stems, collapse = ", "), call. = FALSE)
}

similarity_matrix <- function(df, column) {
  pairs <- df %>%
    transmute(stem1 = strip_ext(file1), stem2 = strip_ext(file2),
              similarity = .data[[column]])
  m <- matrix(NA_real_, length(stems), length(stems), dimnames = list(stems, stems))
  for (i in seq_len(nrow(pairs))) {
    m[pairs$stem1[i], pairs$stem2[i]] <- pairs$similarity[i]
  }
  m[is.na(m)] <- t(m)[is.na(m)]
  if (anyNA(m)) stop("Similarity matrix is incomplete.", call. = FALSE)
  diag(m) <- 1
  m[] <- pmin(pmax(m, 0), 1)
  m
}

mat <- similarity_matrix(raw, sim_col)
# Leaf order only -- matches Figure 6's average-linkage dendrogram exactly,
# so the two figures are directly comparable side by side. No cut is taken.
hc <- hclust(as.dist(1 - mat), method = "average")
ord <- hc$labels[hc$order]

# The ten-cluster cut used for the reported ARI. Hierarchical-clustering
# clusters are contiguous in leaf order, so each becomes a square on the
# heatmap diagonal. Outlining those squares exposes tissue groups that the cut
# splits even when their samples sit next to one another in the dendrogram.
n_tissues <- if (is.na(cluster_count)) length(unique(tissue_by_stem[stems])) else cluster_count
predicted <- cutree(hc, k = n_tissues)
predicted_ordered <- unname(predicted[ord])
predicted_runs <- split(
  seq_along(predicted_ordered),
  cumsum(c(1, predicted_ordered[-1] != predicted_ordered[-length(predicted_ordered)]))
)
predicted_boxes <- data.frame(
  xmin = vapply(predicted_runs, min, numeric(1)) - 0.5,
  xmax = vapply(predicted_runs, max, numeric(1)) + 0.5,
  cluster = vapply(predicted_runs, function(run) predicted_ordered[run[1]], numeric(1)),
  stringsAsFactors = FALSE
)
predicted_boxes$ymin <- length(ord) + 1 - predicted_boxes$xmax
predicted_boxes$ymax <- length(ord) + 1 - predicted_boxes$xmin

# A narrow dark-green strip makes the inferred ten-cluster cut visible next to
# the existing annotated-tissue strip. Small gaps distinguish adjacent inferred
# clusters without assigning biologically suggestive colours to arbitrary
# cluster IDs.
predicted_annotations <- predicted_boxes %>%
  transmute(
    start = xmin + 0.10,
    end = xmax - 0.10,
    center = (xmin + xmax) / 2,
    cluster = cluster
  )

tissue_levels <- unique(tissue_by_stem[stems])
pal <- setNames(hue_pal()(length(tissue_levels)), tissue_levels)

dist_long <- expand.grid(row = ord, col = ord, stringsAsFactors = FALSE) %>%
  mutate(
    distance = 1 - mat[cbind(row, col)],
    row_lab = factor(short_label(row), levels = short_label(ord)),
    col_lab = factor(short_label(col), levels = rev(short_label(ord)))
  )

row_colors <- if (black_accession_labels) {
  rep("black", length(ord))
} else {
  pal[tissue_by_stem[ord]]
}
col_colors <- rev(row_colors)

# Build annotations from contiguous tissue runs in the computed leaf order.
# Their locations therefore remain correct if a different input or similarity
# column rearranges the leaves.
ordered_tissues <- unname(tissue_by_stem[ord])
group_runs <- split(
  seq_along(ordered_tissues),
  cumsum(c(1, ordered_tissues[-1] != ordered_tissues[-length(ordered_tissues)]))
)
group_tissues <- vapply(
  group_runs, function(run) ordered_tissues[run[1]], character(1)
)
group_annotations <- data.frame(
  start = vapply(group_runs, min, numeric(1)) - 0.45,
  end = vapply(group_runs, max, numeric(1)) + 0.45,
  center = vapply(group_runs, function(run) mean(range(run)), numeric(1)),
  label = display_group(group_tissues),
  colour = unname(pal[group_tissues]),
  bracket = unname(BRACKET_GROUP[group_tissues]),
  stringsAsFactors = FALSE
)

n_samples <- length(ord)
# Tiles end at n + 0.5. Leave a visible white gutter before the inferred-cut
# strip, then place the annotated-tissue strip and labels farther outward.
top_cluster_y <- n_samples + 0.68
top_group_y <- n_samples + 1.03
top_label_y <- n_samples + 1.48
right_cluster_x <- n_samples + 0.68
right_group_x <- n_samples + 1.03
right_label_x <- n_samples + 1.48
top_bracket_y <- n_samples + 2.62
top_bracket_tick_y <- n_samples + 2.49
top_bracket_label_y <- n_samples + 3.02
right_bracket_x <- n_samples + 2.62
right_bracket_tick_x <- n_samples + 2.49
right_bracket_label_x <- n_samples + 3.02

# Apply one Muscle bracket to each contiguous run of muscle subtypes. At p=12,
# Arm is displaced from Back/Leg and therefore receives its own short bracket,
# matching the corresponding dendrogram panel.
muscle_members <- which(group_annotations$bracket == "Muscle")
muscle_spans <- if (length(muscle_members) > 0) {
  muscle_runs <- split(
    muscle_members,
    cumsum(c(1, diff(muscle_members) != 1))
  )
  bind_rows(lapply(muscle_runs, function(members) {
    start <- min(group_annotations$start[members])
    end <- max(group_annotations$end[members])
    data.frame(start = start, end = end, center = mean(c(start, end)))
  }))
} else NULL

p <- ggplot(dist_long, aes(x = row_lab, y = col_lab, fill = distance)) +
  geom_tile() +
  geom_rect(
    data = predicted_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, fill = NA, colour = "#7FFF00", linewidth = 0.8
  ) +
  geom_segment(
    data = group_annotations,
    aes(x = start, xend = end, y = top_group_y, yend = top_group_y,
        colour = colour),
    inherit.aes = FALSE, linewidth = 1.8, show.legend = FALSE
  ) +
  geom_segment(
    data = group_annotations,
    aes(x = right_group_x, xend = right_group_x,
        y = n_samples + 1 - end, yend = n_samples + 1 - start,
        colour = colour),
    inherit.aes = FALSE, linewidth = 1.8, show.legend = FALSE
  ) +
  geom_text(
    data = group_annotations,
    aes(x = center, y = top_label_y, label = label),
    inherit.aes = FALSE, angle = 45, hjust = 0, vjust = 0.5,
    size = 3.2, lineheight = 0.9
  ) +
  geom_text(
    data = group_annotations,
    aes(x = right_label_x, y = n_samples + 1 - center, label = label),
    inherit.aes = FALSE, hjust = 0, size = 3.2, lineheight = 0.9
  ) +
  scale_fill_viridis_c(name = "1 − Jaccard\n(distance)", option = "C", direction = -1) +
  scale_colour_identity() +
  coord_cartesian(
    xlim = c(0.5, n_samples + 3.25),
    ylim = c(0.5, n_samples + 3.25),
    clip = "off",
    expand = FALSE
  ) +
  labs(
    x = NULL, y = NULL, title = NULL, subtitle = precision_label,
    tag = panel_label
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = row_colors, size = 8),
    axis.text.y = element_text(colour = col_colors, size = 8),
    panel.grid = element_blank(),
    aspect.ratio = 1,
    legend.position = if (show_legend) "right" else "none",
    plot.subtitle = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1),
    plot.margin = margin(t = 50, r = 95, b = 10, l = 10)
  )

if (show_cluster_strip) {
  p <- p +
    geom_segment(
      data = predicted_annotations,
      aes(x = start, xend = end, y = top_cluster_y, yend = top_cluster_y),
      inherit.aes = FALSE, colour = "#166534", linewidth = 1.8,
      show.legend = FALSE
    ) +
    geom_segment(
      data = predicted_annotations,
      aes(x = right_cluster_x, xend = right_cluster_x,
          y = n_samples + 1 - end, yend = n_samples + 1 - start),
      inherit.aes = FALSE, colour = "#166534", linewidth = 1.8,
      show.legend = FALSE
    )
}

if (!is.null(muscle_spans)) {
  for (i in seq_len(nrow(muscle_spans))) {
    muscle_span <- unlist(muscle_spans[i, c("start", "end")])
    muscle_mid <- muscle_spans$center[i]
    p <- p +
      annotate(
        "segment", x = muscle_span[1], xend = muscle_span[2],
        y = top_bracket_y, yend = top_bracket_y, linewidth = 0.4
      ) +
      annotate(
        "segment", x = c(muscle_span[1], muscle_span[2]),
        xend = c(muscle_span[1], muscle_span[2]),
        y = top_bracket_y, yend = top_bracket_tick_y, linewidth = 0.4
      ) +
      annotate(
        "text", x = muscle_mid, y = top_bracket_label_y,
        label = "Muscle", size = 3.2
      ) +
      annotate(
        "segment", x = right_bracket_x, xend = right_bracket_x,
        y = n_samples + 1 - muscle_span[2],
        yend = n_samples + 1 - muscle_span[1], linewidth = 0.4
      ) +
      annotate(
        "segment", x = right_bracket_x, xend = right_bracket_tick_x,
        y = c(n_samples + 1 - muscle_span[2], n_samples + 1 - muscle_span[1]),
        yend = c(n_samples + 1 - muscle_span[2], n_samples + 1 - muscle_span[1]),
        linewidth = 0.4
      ) +
      annotate(
        "text", x = right_bracket_label_x,
        y = n_samples + 1 - muscle_mid,
        label = "Muscle", angle = 270, size = 3.2
      )
  }
}

CairoPNG(out_png, width = 9, height = 8, units = "in", res = 300, bg = "white")
print(p)
dev.off()
message("Wrote: ", out_png)
