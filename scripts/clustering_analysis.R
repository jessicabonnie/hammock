#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(stats)
  library(utils)
  library(cluster)  # silhouette
  suppressWarnings({
    library(aricode)  # NMI
  })
})

# Reuse shared helpers from scripts/utils.R with robust path resolution
try({
  # Try to detect script directory when run via Rscript
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  script_dir <- if (length(file_arg) > 0) dirname(sub("^--file=", "", file_arg[1])) else getwd()
  candidates <- c(
    file.path(script_dir, "utils.R"),
    file.path(script_dir, "scripts", "utils.R"),
    "scripts/utils.R",
    "utils.R"
  )
  for (p in candidates) {
    if (file.exists(p)) { source(p, local = TRUE); break }
  }
}, silent = TRUE)

## All helper functions are provided by scripts/utils.R (sourced above)

# --- Core analysis ---

perform_hclust_labels <- function(sim_matrix, n_clusters, linkage) {
  # Build a proper distance matrix from similarity
  sim_matrix <- as.matrix(sim_matrix)
  storage.mode(sim_matrix) <- "numeric"
  dist_mat <- 1 - sim_matrix
  dist_mat[dist_mat < 0] <- 0
  diag(dist_mat) <- 0
  d <- stats::as.dist(dist_mat)
  # Map method to Ward.D2 for parity with sklearn
  method <- switch(tolower(linkage),
                   ward = 'ward.D2',
                   complete = 'complete',
                   average = 'average',
                   single = 'single',
                   linkage)
  # Prefer fastcluster if available; fall back to stats::hclust
  if (requireNamespace('fastcluster', quietly = TRUE)) {
    fit <- fastcluster::hclust(d, method = method)
  } else {
    fit <- stats::hclust(d, method = method)
  }
  stats::cutree(fit, k = n_clusters)
}

evaluate_nmi_long_table <- function(sim_matrix, true_labels_map, n_clusters_range, linkage_methods) {
  # Ensure we have a numeric square matrix with names
  sim_matrix <- as.matrix(sim_matrix)
  storage.mode(sim_matrix) <- "numeric"
  files <- rownames(sim_matrix)
  # Align true labels to matrix order and sanitize
  true_labels <- ifelse(files %in% names(true_labels_map), true_labels_map[files], 'unknown')
  true_labels <- as.character(true_labels)
  true_labels[is.na(true_labels) | !nzchar(true_labels)] <- 'unknown'
  # Convert similarity to clipped distance once
  dist_sq <- 1 - sim_matrix
  dist_sq[dist_sq < 0] <- 0
  diag(dist_sq) <- 0
  d <- stats::as.dist(dist_sq)
  # Debug prints: matrix and labels sanity
  message(sprintf("[R] sim_matrix: %dx%d, symmetric=%s", nrow(sim_matrix), ncol(sim_matrix), isTRUE(all.equal(sim_matrix, t(sim_matrix)))))
  message(sprintf("[R] labels: %d samples, unique=%d", length(true_labels), length(unique(true_labels))))

  rows <- list()
  for (meth in linkage_methods) {
    message(sprintf("[R] linkage=%s", meth))
    # Build one dendrogram per linkage, then cut for all k
    method_mapped <- switch(tolower(meth),
                            ward = 'ward.D2',
                            complete = 'complete',
                            average = 'average',
                            single = 'single',
                            meth)
    fit <- tryCatch({
      if (requireNamespace('fastcluster', quietly = TRUE)) fastcluster::hclust(d, method = method_mapped)
      else stats::hclust(d, method = method_mapped)
    }, error = function(e) NULL)
    if (is.null(fit)) {
      for (k in n_clusters_range) {
        rows[[length(rows) + 1]] <- data.frame(n_clusters = as.integer(k), linkage = as.character(meth), nmi = NA_real_, silhouette = NA_real_, stringsAsFactors = FALSE)
      }
      next
    }
    for (k in n_clusters_range) {
      labs <- tryCatch(stats::cutree(fit, k = k), error = function(e) NULL)
      nmi <- NA_real_; sil <- NA_real_
      if (!is.null(labs)) {
        labs <- as.integer(labs)
        if (length(labs) == length(true_labels)) {
          nmi <- tryCatch({ aricode::NMI(true_labels, labs, variant = "sum") }, error = function(e) NA_real_)
          if (length(unique(labs)) >= 2) {
            sil <- tryCatch({
              sw <- cluster::silhouette(labs, d)
              mean(sw[, 3])
            }, error = function(e) NA_real_)
          } else {
            sil <- NA_real_
          }
          if (k %in% c(min(n_clusters_range), max(n_clusters_range), 10, 20)) {
            message(sprintf("[R] k=%d: NMI=%s, sil=%s",
                            k,
                            ifelse(is.na(nmi), "NA", sprintf("%.4f", nmi)),
                            ifelse(is.na(sil), "NA", sprintf("%.4f", sil))))
          }
        }
      }
      rows[[length(rows) + 1]] <- data.frame(
        n_clusters = as.integer(k),
        linkage = as.character(meth),
        nmi = as.numeric(nmi),
        silhouette = as.numeric(sil),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

infer_label_from_basename <- function(basename) {
  token <- strsplit(basename, '-', fixed = TRUE)[[1]][1]
  if (startsWith(token, 'f') && nchar(token) > 1) token <- substring(token, 2)
  token
}

# Public function: returns the long-table and writes it to file
clustering_analysis <- function(hammock_output, accession_key,
                                clusters = 2:30,
                                linkage_methods = c('average', 'complete'),
                                out = NULL,
                                include_life_stage = TRUE) {
  file_format <- detect_file_format(hammock_output)
  if (file_format == 'hammock') {
    sim_df <- parse_hammock_format(hammock_output)
  } else {
    sim_df <- parse_bedtools_format(hammock_output)
  }
  sim_matrix <- as.matrix(sim_df)
  labels_map <- tryCatch(load_accession_key(accession_key, include_life_stage = include_life_stage), error = function(e) NULL)
  files_in_matrix <- rownames(sim_matrix)
  has_overlap <- !is.null(labels_map) && any(files_in_matrix %in% names(labels_map))
  
  if (has_overlap) {
    aligned_result <- align_matrix_and_labels(sim_matrix, labels_map, quiet = FALSE)
    sim_matrix <- aligned_result[[1]]
    effective_labels <- aligned_result[[2]]
    
    # Report on Life_stage usage
    if (has_life_stage_info(effective_labels)) {
      message(sprintf("[R] Using enhanced labels with Life_stage information (%d unique labels)", 
                     length(unique(effective_labels))))
    } else {
      message(sprintf("[R] Using standard Biosample labels (%d unique labels)", 
                     length(unique(effective_labels))))
    }
  } else {
    # Fallback to inference
    eff <- setNames(vapply(files_in_matrix, infer_label_from_basename, character(1)), files_in_matrix)
    attr(eff, "has_life_stage") <- FALSE
    effective_labels <- eff
    message("[R] No overlap with accession key; using inferred labels from basenames")
  }

  long_table <- evaluate_nmi_long_table(sim_matrix, effective_labels, clusters, linkage_methods)

  # Attach filename-derived parameters
  params <- if (file_format == 'hammock') extract_params_from_filename(hammock_output) else
    list(mode = NA, klen = NA, window = NA, precision = NA, sketch = NA, subA = NA, subB = NA, expA = NA)
  if (isTRUE(toupper(params$mode) == 'C') && (is.na(params$expA) || !is.finite(params$expA))) {
    params$expA <- detect_hammock_expA(hammock_output)
  }

  # Insert parameter columns at the front in the specified order
  long_table$sketch <- params$sketch
  long_table$mode <- params$mode
  long_table$precision <- params$precision
  long_table$expA <- params$expA
  long_table$subA <- params$subA
  long_table$subB <- params$subB
  long_table$window <- params$window
  long_table$klen <- params$klen
  # Reorder columns
  long_table <- long_table[, c('sketch','mode','precision','expA','subA','subB','window','klen',
                               'n_clusters','linkage','nmi','silhouette')]

  # Write output to file (default path if none provided)
  if (is.null(out) || !nzchar(out)) {
    in_path <- normalizePath(hammock_output, mustWork = FALSE)
    stem <- tools::file_path_sans_ext(basename(in_path))
    out <- file.path(dirname(in_path), paste0(stem, '_cluster.tsv'))
  }
  utils::write.table(long_table[order(long_table$n_clusters, long_table$linkage), ], file = out,
                     sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
  long_table
}


# --- CLI entry point ---
if (identical(environment(), globalenv()) && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop('Usage: clustering_analysis.R <hammock_or_bedtools_output> <accession_key.tsv> [out.tsv] [include_life_stage=TRUE]')
  }
  input_file <- args[[1]]
  acc_key <- args[[2]]
  out_file <- if (length(args) >= 3) args[[3]] else NULL
  include_life_stage <- if (length(args) >= 4) as.logical(args[[4]]) else TRUE
  
  res <- clustering_analysis(input_file, acc_key, out = out_file, include_life_stage = include_life_stage)
  # Also echo the output path for convenience
  message(sprintf('Long-table NMI results written to: %s', if (is.null(out_file) || !nzchar(out_file)) {
    in_path <- normalizePath(input_file, mustWork = FALSE)
    stem <- tools::file_path_sans_ext(basename(in_path))
    file.path(dirname(in_path), paste0(stem, '_cluster.tsv'))
  } else out_file))
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a) && nchar(as.character(a)) > 0) a else b


