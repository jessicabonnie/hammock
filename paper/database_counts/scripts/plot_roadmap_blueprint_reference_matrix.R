#!/usr/bin/env Rscript

# Roadmap–BLUEPRINT H3K27ac reference-mismatch comparison matrix
#
# Roadmap Epigenomics: 462 native hg19 BED files
# BLUEPRINT:           173 native hg38 BED files
#
# Diagonal blocks represent comparisons possible using native BED coordinates.
# Off-diagonal blocks represent comparisons blocked by the hg19/hg38 mismatch.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

n_roadmap <- 462
n_blueprint <- 173
n_total <- n_roadmap + n_blueprint

roadmap_pairs <- choose(n_roadmap, 2)
blueprint_pairs <- choose(n_blueprint, 2)
within_reference_pairs <- roadmap_pairs + blueprint_pairs
blocked_pairs <- n_roadmap * n_blueprint
all_pairs <- choose(n_total, 2)

blocked_percent <- 100 * blocked_pairs / all_pairs
comparable_percent <- 100 * within_reference_pairs / all_pairs

stopifnot(blocked_pairs == 79926)
stopifnot(all_pairs == 201295)

blocks <- data.frame(
  block = c(
    "Roadmap × Roadmap",
    "BLUEPRINT × BLUEPRINT",
    "Roadmap × BLUEPRINT",
    "BLUEPRINT × Roadmap"
  ),
  comparison_type = c(
    "Directly comparable",
    "Directly comparable",
    "Reference mismatch",
    "Reference mismatch"
  ),
  xmin = c(0, n_roadmap, 0, n_roadmap),
  xmax = c(n_roadmap, n_total, n_roadmap, n_total),
  ymin = c(0, n_roadmap, n_roadmap, 0),
  ymax = c(n_roadmap, n_total, n_total, n_roadmap)
)

# Keep labels deliberately compact so that they fit within the smaller blocks.
block_labels <- transform(
  blocks,
  x = (xmin + xmax) / 2,
  y = (ymin + ymax) / 2,
  label = c(
    paste0("Roadmap\nhg19\n", comma(roadmap_pairs), " pairs"),
    paste0("BLUEPRINT\nhg38\n", comma(blueprint_pairs), " pairs"),
    paste0("Reference mismatch\n", comma(blocked_pairs), " blocked"),
    paste0("Reference mismatch\n", comma(blocked_pairs), " blocked")
  ),
  text_size = c(4.1, 3.1, 3.2, 3.2)
)

axis_breaks <- c(
  n_roadmap / 2,
  n_roadmap + n_blueprint / 2
)

axis_labels <- c(
  paste0("Roadmap\nhg19 (n = ", comma(n_roadmap), ")"),
  paste0("BLUEPRINT\nhg38 (n = ", comma(n_blueprint), ")")
)

comparison_matrix <- ggplot() +
  geom_rect(
    data = blocks,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = comparison_type
    ),
    linewidth = 0.8
  ) +
  geom_vline(xintercept = n_roadmap, linewidth = 1) +
  geom_hline(yintercept = n_roadmap, linewidth = 1) +
  geom_text(
    data = block_labels,
    aes(x = x, y = y, label = label, size = text_size),
    lineheight = 1.05,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_size_identity() +
  scale_x_continuous(
    breaks = axis_breaks,
    labels = axis_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = axis_breaks,
    labels = axis_labels,
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "Directly comparable" = "#D9EAD3",
      "Reference mismatch" = "#F4CCCC"
    )
  ) +
  coord_fixed(clip = "off") +
  labs(
    title = "Reference mismatch blocks cross-resource\nH3K27ac comparisons",
    subtitle = paste0(
      comma(n_roadmap), " Roadmap hg19 files and ",
      comma(n_blueprint), " BLUEPRINT hg38 files"
    ),
    x = NULL,
    y = NULL,
    fill = NULL,
    caption = paste0(
      comma(blocked_pairs), " of ", comma(all_pairs),
      " unique pairs (", number(blocked_percent, accuracy = 0.1),
      "%) cross the hg19–hg38 boundary."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9, margin = margin(t = 7)),
    axis.text.y = element_text(size = 9, margin = margin(r = 7)),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 15, lineheight = 1.05),
    plot.subtitle = element_text(size = 10.5, margin = margin(b = 10)),
    plot.caption = element_text(hjust = 0, size = 9, margin = margin(t = 10)),
    plot.margin = margin(16, 24, 16, 18)
  )

print(comparison_matrix)

output_dir <- "paper/database_counts/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# PNG only for now. A wider canvas gives the title and axis labels room.
ggsave(
  filename = file.path(output_dir, "roadmap_blueprint_reference_matrix.png"),
  plot = comparison_matrix,
  width = 9,
  height = 7.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

summary_table <- data.frame(
  category = c(
    "Roadmap within-reference pairs",
    "BLUEPRINT within-reference pairs",
    "Total directly comparable pairs",
    "Cross-reference blocked pairs",
    "All possible unique pairs"
  ),
  pair_count = c(
    roadmap_pairs,
    blueprint_pairs,
    within_reference_pairs,
    blocked_pairs,
    all_pairs
  ),
  percent_of_all_pairs = c(
    100 * roadmap_pairs / all_pairs,
    100 * blueprint_pairs / all_pairs,
    comparable_percent,
    blocked_percent,
    100
  )
)

write.table(
  summary_table,
  file = file.path(output_dir, "roadmap_blueprint_reference_matrix_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(summary_table)
