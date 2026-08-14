#!/usr/bin/env Rscript

# Plot the five shared Roadmap/BLUEPRINT histone marks with the largest number
# of cross-reference comparisons blocked by the hg19/hg38 mismatch.
#
# For each marker, show three side-by-side bars:
# 1. within-Roadmap comparisons on hg19: choose(R, 2)
# 2. within-BLUEPRINT comparisons on hg38: choose(B, 2)
# 3. blocked Roadmap x BLUEPRINT comparisons: R * B

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(Cairo)
})

input_path <- "paper/database_counts/results/roadmap_blueprint_top5_markers.tsv"
output_dir <- "paper/figures"
output_path <- file.path(
  output_dir,
  "roadmap_blueprint_top5_marker_pairwise_comparisons.png"
)

if (!file.exists(input_path)) {
  stop(
    paste0(
      "Input file not found: ", input_path, "\n",
      "Run first:\n",
      "python paper/database_counts/scripts/count_filer_roadmap_blueprint_top_markers.py"
    )
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(
  input_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "marker",
  "roadmap_hg19_bed_count",
  "blueprint_hg38_bed_count",
  "blocked_cross_reference_pairs",
  "blocked_percent"
)

missing_columns <- setdiff(required_columns, names(counts))
if (length(missing_columns) > 0) {
  stop(
    paste(
      "Input is missing required columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

if (nrow(counts) == 0) {
  stop("Input table contains no markers")
}

# Recompute the three comparison classes directly from the file counts so that
# the plotting script remains transparent and independently checkable.
counts$roadmap_within_pairs <- choose(counts$roadmap_hg19_bed_count, 2)
counts$blueprint_within_pairs <- choose(counts$blueprint_hg38_bed_count, 2)
counts$blocked_cross_reference_pairs <-
  counts$roadmap_hg19_bed_count * counts$blueprint_hg38_bed_count

# Rank markers by the blocked cross-reference workload.
counts <- counts[
  order(counts$blocked_cross_reference_pairs, decreasing = TRUE),
]
counts$marker <- factor(counts$marker, levels = counts$marker)

plot_data <- rbind(
  data.frame(
    marker = counts$marker,
    comparison_class = "Roadmap hg19: within-resource pairs",
    pair_count = counts$roadmap_within_pairs,
    stringsAsFactors = FALSE
  ),
  data.frame(
    marker = counts$marker,
    comparison_class = "BLUEPRINT hg38: within-resource pairs",
    pair_count = counts$blueprint_within_pairs,
    stringsAsFactors = FALSE
  ),
  data.frame(
    marker = counts$marker,
    comparison_class = "Blocked Roadmap hg19 x BLUEPRINT hg38 pairs",
    pair_count = counts$blocked_cross_reference_pairs,
    stringsAsFactors = FALSE
  )
)

plot_data$comparison_class <- factor(
  plot_data$comparison_class,
  levels = c(
    "Roadmap hg19: within-resource pairs",
    "BLUEPRINT hg38: within-resource pairs",
    "Blocked Roadmap hg19 x BLUEPRINT hg38 pairs"
  )
)

# A three-category dodge with width 0.82 places the third bar one-third of the
# dodge width to the right of the marker center. Set that x position explicitly
# because an annotation layer containing only the red category would otherwise
# be centered by position_dodge().
dodge_width <- 0.82
third_bar_offset <- dodge_width / 3

blocked_labels <- data.frame(
  marker_x = seq_len(nrow(counts)) + third_bar_offset,
  pair_count = counts$blocked_cross_reference_pairs,
  label = paste0(
    comma(counts$blocked_cross_reference_pairs),
    "\n(",
    number(counts$blocked_percent, accuracy = 0.1),
    "% of all pairs)"
  ),
  stringsAsFactors = FALSE
)

within_labels <- plot_data[
  plot_data$comparison_class !=
    "Blocked Roadmap hg19 x BLUEPRINT hg38 pairs",
]
within_labels$marker_x <-
  match(within_labels$marker, levels(counts$marker)) +
  ifelse(
    within_labels$comparison_class ==
      "Roadmap hg19: within-resource pairs",
    -third_bar_offset,
    0
  )

file_labels <- data.frame(
  marker = counts$marker,
  label = paste0(
    "Roadmap n = ", comma(counts$roadmap_hg19_bed_count),
    "\nBLUEPRINT n = ", comma(counts$blueprint_hg38_bed_count)
  ),
  stringsAsFactors = FALSE
)

max_pairs <- max(plot_data$pair_count)
lower_room <- max_pairs * 0.16
upper_room <- max_pairs * 0.14

pairwise_plot <- ggplot(
  plot_data,
  aes(x = marker, y = pair_count, fill = comparison_class)
) +
  geom_col(
    position = position_dodge(width = dodge_width),
    width = 0.72
  ) +
  geom_text(
    data = within_labels,
    aes(x = marker_x, y = pair_count, label = comma(pair_count)),
    inherit.aes = FALSE,
    vjust = -0.35,
    size = 3.5,
    show.legend = FALSE
  ) +
  geom_text(
    data = blocked_labels,
    aes(x = marker_x, y = pair_count, label = label),
    inherit.aes = FALSE,
    vjust = -0.25,
    size = 3.5,
    lineheight = 0.95
  ) +
  geom_text(
    data = file_labels,
    aes(x = marker, y = -lower_room * 0.32, label = label),
    inherit.aes = FALSE,
    size = 3.3,
    lineheight = 0.95,
    vjust = 1
  ) +
  scale_y_continuous(
    labels = label_number(scale_cut = cut_short_scale()),
    limits = c(-lower_room, max_pairs + upper_room),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(
    values = c(
      "Roadmap hg19: within-resource pairs" = "#88B04B",
      "BLUEPRINT hg38: within-resource pairs" = "#6FA8DC",
      "Blocked Roadmap hg19 x BLUEPRINT hg38 pairs" = "#D9827C"
    )
  ) +
  labs(
    title = "Reference mismatch blocks histone-mark ChIP-seq peak comparisons",
    x = NULL,
    y = "Number of pairwise ChIP-seq peak BED-file comparisons",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold", size = 12.5),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 15, lineheight = 1.05),
    plot.margin = margin(16, 18, 16, 18)
  )

print(pairwise_plot)

CairoPNG(
  filename = output_path,
  width = 10.8,
  height = 7.2,
  units = "in",
  res = 300,
  bg = "white"
)
print(pairwise_plot)
dev.off()

message("Wrote: ", output_path)
