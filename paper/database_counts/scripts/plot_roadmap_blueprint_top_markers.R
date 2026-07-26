#!/usr/bin/env Rscript

# Plot the five shared Roadmap/BLUEPRINT histone marks with the largest number
# of cross-reference comparisons blocked by the hg19/hg38 mismatch.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

input_path <- "paper/database_counts/results/roadmap_blueprint_top5_markers.tsv"
output_dir <- "paper/database_counts/results"
output_path <- file.path(output_dir, "roadmap_blueprint_top5_marker_pairwise_comparisons.png")

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
  "directly_comparable_pairs",
  "blocked_cross_reference_pairs",
  "all_possible_pairs",
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

# Preserve the Python ranking, with the highest blocked workload at the top.
counts <- counts[order(counts$blocked_cross_reference_pairs, decreasing = TRUE), ]
counts$marker <- factor(counts$marker, levels = rev(counts$marker))

plot_data <- rbind(
  data.frame(
    marker = counts$marker,
    comparison_type = "Directly comparable within reference",
    pair_count = counts$directly_comparable_pairs,
    stringsAsFactors = FALSE
  ),
  data.frame(
    marker = counts$marker,
    comparison_type = "Blocked across hg19 and hg38",
    pair_count = counts$blocked_cross_reference_pairs,
    stringsAsFactors = FALSE
  )
)

plot_data$comparison_type <- factor(
  plot_data$comparison_type,
  levels = c(
    "Directly comparable within reference",
    "Blocked across hg19 and hg38"
  )
)

annotation_data <- counts
annotation_data$file_label <- paste0(
  comma(annotation_data$roadmap_hg19_bed_count), " Roadmap hg19 + ",
  comma(annotation_data$blueprint_hg38_bed_count), " BLUEPRINT hg38"
)
annotation_data$blocked_label <- paste0(
  comma(annotation_data$blocked_cross_reference_pairs),
  " blocked (",
  number(annotation_data$blocked_percent, accuracy = 0.1),
  "%)"
)

max_pairs <- max(counts$all_possible_pairs)
label_padding <- max_pairs * 0.025

pairwise_plot <- ggplot(
  plot_data,
  aes(x = pair_count, y = marker, fill = comparison_type)
) +
  geom_col(width = 0.68) +
  geom_text(
    data = annotation_data,
    aes(
      x = all_possible_pairs + label_padding,
      y = marker,
      label = blocked_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.7
  ) +
  geom_text(
    data = annotation_data,
    aes(
      x = 0,
      y = marker,
      label = file_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    nudge_y = 0.32,
    size = 3.2
  ) +
  scale_x_continuous(
    labels = label_number(scale_cut = cut_short_scale()),
    expand = expansion(mult = c(0, 0.26))
  ) +
  scale_fill_manual(
    values = c(
      "Directly comparable within reference" = "#9CC4A1",
      "Blocked across hg19 and hg38" = "#D9827C"
    )
  ) +
  labs(
    title = "Reference mismatch blocks comparisons across major histone marks",
    subtitle = paste0(
      "Top five shared Roadmap and BLUEPRINT marks ranked by blocked comparisons; ",
      "each full bar represents all possible file pairs"
    ),
    x = "Possible pairwise BED-file comparisons",
    y = NULL,
    fill = NULL,
    caption = paste0(
      "Roadmap peak BEDs use hg19; BLUEPRINT peak BEDs use hg38. ",
      "Blocked pairs are Roadmap × BLUEPRINT comparisons that cannot be evaluated ",
      "directly with coordinate overlap on the native files."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 15, lineheight = 1.05),
    plot.subtitle = element_text(size = 10.5, lineheight = 1.1, margin = margin(b = 12)),
    plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.1),
    plot.margin = margin(16, 28, 16, 16)
  )

print(pairwise_plot)

ggsave(
  filename = output_path,
  plot = pairwise_plot,
  width = 10.5,
  height = 6.7,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Wrote: ", output_path)
