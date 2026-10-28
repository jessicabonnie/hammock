#!/usr/bin/env Rscript

# Stack the p=12, p=18, and p=24 fixed-seed dendrogram panels into a single
# supplementary figure. By default, the component panels are generated
# in temporary files by plot_sequence_tissue_clustering.R with panel labels A,
# B, and C; three explicit panel paths may instead be supplied as arguments.

required_packages <- c("png")
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
cluster_count <- if (length(argv) == 2) as.integer(argv[2]) else 10L
if (length(argv) < 3) {
  plot_script <- file.path(script_dir, "plot_sequence_tissue_clustering.R")
  data_dir <- file.path(script_dir, "data")
  key_tsv <- file.path(data_dir, "maurano_filenames_key.tsv")
  panel_pdfs <- vapply(
    precisions,
    function(p) tempfile(
      pattern = paste0("sequence_tissue_clustering_p", p, "_"),
      fileext = ".pdf"
    ),
    character(1)
  )
  panel_png_bases <- sub("\\.pdf$", "", panel_pdfs)
  inputs <- paste0(panel_png_bases, ".png")
  on.exit(unlink(c(panel_pdfs, inputs)), add = TRUE)
  for (i in seq_along(precisions)) {
    input_csv <- file.path(
      data_dir,
      sprintf("p%d_seed00042_k10_w30.csv", precisions[i])
    )
    plot_args <- c(
      plot_script, input_csv, key_tsv, panel_pdfs[i],
      "jaccard_similarity_ie", panel_labels[i], "0.05", "4.5", "true",
      as.character(precisions[i]), as.character(cluster_count), "", "none",
      "tissue"
    )
    status <- system2(
      "Rscript",
      vapply(plot_args, shQuote, character(1))
    )
    if (status != 0) {
      stop("Failed to generate dendrogram panel for p=", precisions[i],
           call. = FALSE)
    }
    status <- system2(
      "pdftoppm",
      c("-png", "-r", "300", "-singlefile", shQuote(panel_pdfs[i]),
        shQuote(panel_png_bases[i]))
    )
    if (status != 0 || !file.exists(inputs[i])) {
      stop("Failed to rasterize dendrogram panel for p=", precisions[i],
           call. = FALSE)
    }
  }
} else {
  inputs <- argv[1:3]
}
default_output <- file.path(
  repo_root, "paper", "figures", "sequence_tissue_clustering_precision.png"
)
output <- if (length(argv) >= 4) {
  argv[4]
} else if (length(argv) < 3 && length(argv) >= 1) {
  argv[1]
} else {
  default_output
}

missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs) > 0) {
  stop("Input file(s) not found: ", paste(missing_inputs, collapse = ", "),
       call. = FALSE)
}

panels <- lapply(inputs, png::readPNG)
panel_dims <- vapply(panels, function(x) paste(dim(x)[1:2], collapse = "x"), character(1))
if (length(unique(panel_dims)) != 1) {
  stop("Component panels must have identical pixel dimensions.", call. = FALSE)
}

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
# Trim only outer white rows. Retaining a small fixed pad keeps panel letters
# clear while removing the repeated device margins between stacked rows.
trim_vertical <- function(img, pad = 18L) {
  rgb <- img[, , seq_len(min(3, dim(img)[3])), drop = FALSE]
  nonwhite <- apply(rgb < 0.995, 1, any)
  occupied <- which(nonwhite)
  lo <- max(1L, min(occupied) - pad)
  hi <- min(dim(img)[1], max(occupied) + pad)
  img[lo:hi, , , drop = FALSE]
}
panels <- lapply(panels, trim_vertical)
row_heights <- vapply(panels, function(x) dim(x)[1], integer(1))

combined <- array(1, dim = c(sum(row_heights), dim(panels[[1]])[2], 3))
row_start <- 1L
for (panel in panels) {
  row_end <- row_start + dim(panel)[1] - 1L
  combined[row_start:row_end, , ] <- panel[, , seq_len(3), drop = FALSE]
  row_start <- row_end + 1L
}
png::writePNG(combined, output)

message("Wrote: ", output)
