#!/usr/bin/env Rscript

# Stack two dendrograms computed from the same representative k=20, w=20,
# p=18, seed-12 matrix. Panel A shows an eight-cluster cut and panel B a
# ten-cluster cut. Component styling matches the three-precision dendrogram:
# average linkage on 1-J, inferred-cluster strip, compact vertical stacking,
# and a merge-preserving axis break. The 0.35 break is the largest 0.05 step
# below the first merge (0.3671486), so no observed merge is compressed.

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
output <- if (length(argv) >= 1) argv[1] else file.path(
  repo_root, "paper", "figures", "sequence_tissue_clustering_k20_w20_cluster_cuts.png"
)

plot_script <- file.path(script_dir, "plot_sequence_tissue_clustering.R")
data_dir <- file.path(script_dir, "data")
input_csv <- file.path(data_dir, "p18_seed00012_k20_w20.csv")
key_tsv <- file.path(data_dir, "maurano_filenames_key.tsv")
for (path in c(plot_script, input_csv, key_tsv)) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
}

cuts <- c(8L, 10L)
panel_labels <- c("A", "B")
panel_titles <- c("Eight-cluster cut", "Ten-cluster cut")
panel_pdfs <- vapply(
  cuts,
  function(cut) tempfile(pattern = paste0("k20_w20_cut", cut, "_"), fileext = ".pdf"),
  character(1)
)
panel_png_bases <- vapply(
  cuts,
  function(cut) tempfile(pattern = paste0("k20_w20_cut", cut, "_")),
  character(1)
)
panel_pngs <- paste0(panel_png_bases, ".png")
on.exit(unlink(c(panel_pdfs, panel_pngs)), add = TRUE)

for (i in seq_along(cuts)) {
  highlight_splits_from <- if (cuts[i] == 10L) "8" else ""
  plot_args <- c(
    plot_script, input_csv, key_tsv, panel_pdfs[i],
    "jaccard_similarity_ie", panel_labels[i], "0.35", "4.5", "true",
    "18", as.character(cuts[i]), panel_titles[i], "none", "tissue",
    highlight_splits_from
  )
  status <- system2("Rscript", vapply(plot_args, shQuote, character(1)))
  if (status != 0) {
    stop("Failed to generate dendrogram panel for cut=", cuts[i], call. = FALSE)
  }
  status <- system2(
    "pdftoppm",
    c("-png", "-r", "300", "-singlefile", shQuote(panel_pdfs[i]),
      shQuote(panel_png_bases[i]))
  )
  if (status != 0 || !file.exists(panel_pngs[i])) {
    stop("Failed to rasterize dendrogram panel for cut=", cuts[i], call. = FALSE)
  }
}

panels <- lapply(panel_pngs, function(path) {
  image <- png::readPNG(path)
  image[, , seq_len(min(3, dim(image)[3])), drop = FALSE]
})
panel_dims <- vapply(panels, function(x) paste(dim(x)[1:2], collapse = "x"), character(1))
if (length(unique(panel_dims)) != 1) {
  stop("Component panels must have identical pixel dimensions.", call. = FALSE)
}

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

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
combined <- array(1, dim = c(sum(row_heights), dim(panels[[1]])[2], 3))
row_start <- 1L
for (panel in panels) {
  row_end <- row_start + dim(panel)[1] - 1L
  combined[row_start:row_end, , ] <- panel
  row_start <- row_end + 1L
}
png::writePNG(combined, output)
message("Wrote: ", output)
