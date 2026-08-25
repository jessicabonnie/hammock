#!/usr/bin/env Rscript

# Place p=12, p=18, and p=24 distance heatmaps side by side as a visual
# alternative to the stacked precision dendrograms. Each component heatmap is
# ordered by its own average-linkage hierarchy, matching the corresponding
# dendrogram panel rather than imposing a shared display order.

required_packages <- c("png", "grid", "Cairo", "viridisLite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

argv <- commandArgs(trailingOnly = TRUE)
precisions <- c(12, 18, 24)
panel_labels <- c("A", "B", "C")
generated_inputs <- length(argv) < 3
if (generated_inputs) {
  plot_script <- file.path(script_dir, "plot_distance_heatmap.R")
  data_dir <- file.path(script_dir, "data")
  key_tsv <- file.path(data_dir, "maurano_filenames_key.tsv")
  inputs <- vapply(
    precisions,
    function(p) tempfile(
      pattern = paste0("sequence_tissue_distance_heatmap_p", p, "_"),
      fileext = ".png"
    ),
    character(1)
  )
  on.exit(unlink(inputs), add = TRUE)
  for (i in seq_along(precisions)) {
    input_csv <- file.path(
      data_dir,
      sprintf("p%d_seed00000_k10_w30.csv", precisions[i])
    )
    plot_args <- c(
      plot_script, input_csv, key_tsv, inputs[i],
      "jaccard_similarity_ie", panel_labels[i],
      "true", "false", as.character(precisions[i]), "true"
    )
    status <- system2(
      "Rscript",
      vapply(plot_args, shQuote, character(1))
    )
    if (status != 0) {
      stop("Failed to generate heatmap panel for p=", precisions[i],
           call. = FALSE)
    }
  }
} else {
  inputs <- argv[1:3]
}
output <- if (length(argv) >= 4) argv[4] else file.path(
  repo_root, "paper", "figures", "sequence_tissue_distance_heatmap_precision.png"
)

missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs) > 0) {
  stop("Input file(s) not found: ", paste(missing_inputs, collapse = ", "),
       call. = FALSE)
}

panels <- lapply(inputs, png::readPNG)

# Remove empty outer columns before composition. Keeping a small pad preserves
# panel letters and labels while avoiding the doubled white margin between A
# and B.
trim_horizontal <- function(img, pad = 18L) {
  rgb <- img[, , seq_len(min(3, dim(img)[3])), drop = FALSE]
  nonwhite <- apply(rgb < 0.995, 2, any)
  occupied <- which(nonwhite)
  lo <- max(1L, min(occupied) - pad)
  hi <- min(dim(img)[2], max(occupied) + pad)
  img[, lo:hi, , drop = FALSE]
}
panels <- lapply(panels, trim_horizontal)
panel_heights <- vapply(panels, function(x) dim(x)[1], integer(1))
if (length(unique(panel_heights)) != 1) {
  stop("Component panels must have identical pixel heights.", call. = FALSE)
}

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
panel_widths <- vapply(panels, function(x) dim(x)[2], integer(1))
panel_height <- panel_heights[1]

Cairo::CairoPNG(output, width = sum(panel_widths[1:2]),
                height = panel_height * 2,
                bg = "white")
grid::grid.newpage()
grid::pushViewport(grid::viewport(
  layout = grid::grid.layout(
    2, 2,
    widths = grid::unit(panel_widths[1:2], "null")
  )
))
grid::grid.raster(
  panels[[1]], interpolate = FALSE,
  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1)
)
grid::grid.raster(
  panels[[2]], interpolate = FALSE,
  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2)
)
# Panel C occupies one panel-width centered beneath A and B. Its shared legend
# sits in the same lower-row canvas, immediately to its right.
grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
c_width <- panel_widths[3] / sum(panel_widths[1:2])
grid::grid.raster(
  panels[[3]], interpolate = FALSE,
  x = 0.5, y = 0.5,
  width = c_width, height = 1
)

# One shared continuous legend for all three panels. The component plots are
# rendered without legends so their heatmap matrices retain identical sizes.
legend_vp <- grid::viewport(
  x = 0.86, y = 0.50, width = 0.18, height = 0.88
)
grid::pushViewport(legend_vp)
legend_colours <- viridisLite::viridis(256, option = "C", direction = 1)
grid::grid.raster(
  matrix(legend_colours, ncol = 1),
  x = 0.28, y = 0.45, width = 0.22, height = 0.76,
  interpolate = TRUE
)
grid::grid.text(
  "1 − Jaccard\n(distance)",
  x = 0.28, y = 0.88, just = c("center", "bottom"),
  gp = grid::gpar(fontsize = 30)
)
legend_max <- 0.18
tick_values <- c(0, 0.05, 0.10, 0.15)
tick_y <- 0.07 + 0.76 * tick_values / legend_max
for (i in seq_along(tick_values)) {
  grid::grid.segments(
    x0 = 0.39, x1 = 0.45, y0 = tick_y[i], y1 = tick_y[i],
    gp = grid::gpar(col = "black", lwd = 0.8)
  )
  grid::grid.text(
    sprintf("%.2f", tick_values[i]),
    x = 0.50, y = tick_y[i], just = "left",
    gp = grid::gpar(fontsize = 30)
  )
}
grid::popViewport()
grid::popViewport()
grDevices::dev.off()

message("Wrote: ", output)
