#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Cairo)
  library(ggplot2)
  library(patchwork)
  library(png)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
  stop("Could not determine script path. Run with Rscript.", call. = FALSE)
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
repo_root <- normalizePath(file.path(script_dir, "..", ".."))

argv <- commandArgs(trailingOnly = TRUE)
panel_a_png <- if (length(argv) >= 1) argv[1] else file.path(
  repo_root, "paper", "figures", "cross_reference_identity_panel_a.png"
)
panel_b_png <- if (length(argv) >= 2) argv[2] else file.path(
  repo_root, "paper", "figures", "cross_reference_parameter_plateau.png"
)
output_png <- if (length(argv) >= 3) argv[3] else file.path(
  repo_root, "paper", "figures", "cross_reference_identity.png"
)

for (path in c(panel_a_png, panel_b_png)) {
  if (!file.exists(path)) stop("Panel image not found: ", path, call. = FALSE)
}

raster_panel <- function(path) {
  patchwork::wrap_elements(
    full = grid::rasterGrob(png::readPNG(path), interpolate = TRUE)
  )
}

combined <- raster_panel(panel_a_png) / raster_panel(panel_b_png) +
  plot_layout(heights = c(1, 1.03)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "bold", color = "#20262D"),
    plot.tag.position = c(0.006, 0.985),
    plot.margin = margin(8, 8, 8, 8)
  )

dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
CairoPNG(output_png, width = 12.5, height = 8.5, units = "in", res = 300,
         bg = "white")
print(combined)
dev.off()
message("Wrote: ", output_png)
