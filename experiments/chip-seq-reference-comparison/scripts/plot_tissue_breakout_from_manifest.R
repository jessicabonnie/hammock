#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: plot_tissue_breakout_from_manifest.R <manifest_tsv> <output_png>")
}

manifest_path <- args[[1]]
output_png <- args[[2]]

tab <- read.delim(manifest_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

if (!("organism" %in% names(tab)) || !("tissue_type" %in% names(tab))) {
  stop("Manifest must contain organism and tissue_type columns")
}

assembly <- if ("assembly" %in% names(tab)) tab$assembly else rep("", nrow(tab))
reference_build <- if ("reference_build" %in% names(tab)) tab$reference_build else rep("", nrow(tab))

assembly <- trimws(ifelse(is.na(assembly), "", assembly))
reference_build <- trimws(ifelse(is.na(reference_build), "", reference_build))
organism <- trimws(ifelse(is.na(tab$organism), "unknown", tab$organism))
tissue <- trimws(ifelse(is.na(tab$tissue_type) | tab$tissue_type == "", "unknown", tab$tissue_type))

effective_build <- ifelse(
  assembly != "",
  assembly,
  ifelse(
    reference_build != "",
    reference_build,
    ifelse(
      grepl("Mus musculus", organism, fixed = TRUE),
      "mm10",
      ifelse(grepl("Homo sapiens", organism, fixed = TRUE), "GRCh38", "unknown")
    )
  )
)

df <- data.frame(
  tissue_type = tissue,
  organism = organism,
  reference_build = effective_build,
  stringsAsFactors = FALSE
)

agg <- aggregate(
  x = list(count = rep(1, nrow(df))),
  by = list(
    tissue_type = df$tissue_type,
    organism = df$organism,
    reference_build = df$reference_build
  ),
  FUN = sum
)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required in the loaded R module")
}

p <- ggplot2::ggplot(
  agg,
  ggplot2::aes(
    x = tissue_type,
    y = count,
    fill = reference_build,
    group = organism
  )
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge2(width = 0.9, preserve = "single"),
    width = 0.82,
    color = "#333333",
    linewidth = 0.2
  ) +
  ggplot2::labs(
    title = "Tissue types by species with build composition",
    subtitle = "Two bars per tissue (organisms), colored by reference build",
    x = "Tissue type",
    y = "BED files",
    fill = "Reference build"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold")
  )

png_device <- function(filename, width, height, ...) {
  grDevices::png(
    filename = filename,
    width = width,
    height = height,
    type = "cairo",
    units = "in",
    res = 180,
    ...
  )
}

ggplot2::ggsave(
  filename = output_png,
  plot = p,
  width = 13.5,
  height = 5.0,
  device = png_device
)

message("Wrote ", output_png)
