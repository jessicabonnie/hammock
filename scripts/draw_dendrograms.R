#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(stats)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# Resolve scripts directory reliably
args_full <- commandArgs(trailingOnly = FALSE)
script_file <- sub('^--file=', '', args_full[grep('--file=', args_full)])
if (length(script_file) > 0) {
  scripts_dir <- dirname(normalizePath(script_file))
} else {
  scripts_dir <- dirname(normalizePath('scripts/draw_dendrograms.R'))
}
source(file.path(scripts_dir, 'utils.R'))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: draw_dendrograms.R <input_matrix> <accession_key> [output_pdf] [linkages=average,complete,single] [include_life_stage=TRUE]")
}

input_matrix <- args[1]
accession_key <- args[2]
output_pdf <- if (length(args) >= 3) args[3] else paste0(basename(input_matrix), "_dend.pdf")
linkages <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else c("average", "complete", "single")
include_life_stage <- if (length(args) >= 5) as.logical(args[5]) else TRUE

fmt <- detect_file_format(input_matrix)
sim_df <- if (fmt == 'hammock') parse_hammock_format(input_matrix) else parse_bedtools_format(input_matrix)
lab_map <- load_accession_key(accession_key, include_life_stage = include_life_stage)
aligned <- align_matrix_and_labels(sim_df, lab_map, quiet = FALSE)
sim_sub <- aligned[[1]]
labels <- aligned[[2]]

# Check if we have Life_stage information for enhanced visualization
has_life_stage <- has_life_stage_info(labels)
if (has_life_stage) {
  message(sprintf("Using enhanced labels with Life_stage information (%d unique labels)", length(unique(labels))))
} else {
  message(sprintf("Using standard Biosample labels (%d unique labels)", length(unique(labels))))
}

if (nrow(sim_sub) < 2) stop('Need at least 2 samples')

distance_mat <- 1 - as.matrix(sim_sub)
diag(distance_mat) <- 0

pdf(output_pdf, width = 12, height = 8)
on.exit(dev.off(), add = TRUE)

sample_names <- rownames(distance_mat)
lab_vec <- labels[sample_names]
lab_vec[is.na(lab_vec)] <- 'unknown'
uniq_labs <- sort(unique(lab_vec))

# Create color palette - enhanced for Life_stage if available
if (has_life_stage) {
  lab_to_col <- create_color_palette_with_life_stage(uniq_labs)
} else {
  # Traditional color palette
  palette_fun <- if (requireNamespace('scales', quietly = TRUE)) scales::hue_pal() else function(n) grDevices::rainbow(n, s = 0.7, v = 0.9)
  palette <- palette_fun(max(3, length(uniq_labs)))
  lab_to_col <- setNames(palette[seq_along(uniq_labs)], uniq_labs)
}

for (method in linkages) {
  eff_method <- if (tolower(method) == 'ward') 'complete' else method
  hc <- hclust(as.dist(distance_mat), method = eff_method)

  # Convert to dendrogram and color branches/labels (base dendrogram operations)
  dnd <- as.dendrogram(hc)
  label_of_leaf <- lab_vec

  colorize <- function(node) {
    if (is.leaf(node)) {
      lab <- attr(node, 'label')
      col <- lab_to_col[[ label_of_leaf[[lab]] ]] %||% '#888888'
      # color edge leading to leaf
      ep <- attr(node, 'edgePar'); if (is.null(ep)) ep <- list()
      ep$col <- col
      attr(node, 'edgePar') <- ep
      # color label
      np <- attr(node, 'nodePar'); if (is.null(np)) np <- list()
      np$lab.col <- col
      attr(node, 'nodePar') <- np
      return(node)
    } else {
      # recurse on children
      for (i in seq_along(node)) node[[i]] <- colorize(node[[i]])
      leaves <- labels(node)
      labs <- unique(label_of_leaf[leaves])
      col <- if (length(labs) == 1) lab_to_col[[ labs[[1]] ]] %||% '#888888' else '#888888'
      ep <- attr(node, 'edgePar'); if (is.null(ep)) ep <- list()
      ep$col <- col
      attr(node, 'edgePar') <- ep
      return(node)
    }
  }

  dnd <- colorize(dnd)
  
  # Create title with Life_stage information if available
  title_suffix <- ifelse(tolower(method)=='ward', " [ward uses 'complete' on precomputed distances]", "")
  life_stage_suffix <- if (has_life_stage) " with Life Stage" else ""
  main_title <- sprintf('Hierarchical Clustering Dendrogram (%s)%s%s', method, life_stage_suffix, title_suffix)
  
  plot(dnd, main = main_title)
  
  # Create legend - handle long labels for Life_stage enhanced labels
  if (has_life_stage && length(uniq_labs) > 15) {
    # For many enhanced labels, create a more compact legend
    legend('topright', legend = paste("Enhanced labels with Life_stage"), 
           col = "black", pch = NA, cex = 0.7, bty = 'n')
    legend('topright', legend = sprintf("(%d unique combinations)", length(uniq_labs)), 
           col = "black", pch = NA, cex = 0.6, bty = 'n', inset = c(0, 0.05))
  } else {
    # Standard legend for manageable number of labels
    legend('topright', legend = uniq_labs, col = lab_to_col[uniq_labs], pch = 15, 
           cex = if (has_life_stage) 0.6 else 0.8, bty = 'n')
  }
}


