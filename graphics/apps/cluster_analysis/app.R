#.libPaths(c("/home/jbonnie1/R/x86_64-pc-linux-gnu-library/4.3", .libPaths()))
library(shiny)
library(ggplot2)
library(shinyFiles)
library(plotly)
library(cluster)  # For silhouette calculations

# Fix RStudio graphics device recursion issue
if (exists(".rs.restartR")) {
  # Running in RStudio - handle graphics device properly
  options(device = function(...) {
    if (capabilities("png")) {
      png(...)
    } else {
      pdf(...)
    }
  })
  # Reset graphics state to prevent recursion
  while (dev.cur() > 1) dev.off()
}

# Source utils.R at app startup to ensure detect_file_format is available
app_dir <- getwd()
base_dirs <- unique(c(app_dir, dirname(app_dir), dirname(dirname(app_dir)), dirname(dirname(dirname(app_dir)))))
utils_sourced <- FALSE
for (u in file.path(base_dirs, "scripts", "utils.R")) {
  if (file.exists(u)) { 
    try({ source(u); utils_sourced <- TRUE }, silent = TRUE) 
    if (utils_sourced) break 
  }
}
if (!utils_sourced) stop("Cannot find scripts/utils.R - please run from correct directory")

# Normal RStudio graphics (removed server workarounds)

# Helper to extract biosample and lifestage metadata from accession key
extract_accession_metadata <- function(key_path) {
  if (is.null(key_path) || !nzchar(key_path) || !file.exists(key_path)) {
    return(list(biosamples = character(0), lifestages = character(0), full_data = NULL))
  }
  
  df <- tryCatch({
    if (requireNamespace("data.table", quietly = TRUE)) data.table::fread(key_path, sep = "\t", header = TRUE, showProgress = FALSE)
    else utils::read.table(key_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(df) || !nrow(df)) {
    return(list(biosamples = character(0), lifestages = character(0), full_data = NULL))
  }
  
  # Extract unique biosamples and lifestages
  biosamples <- character(0)
  lifestages <- character(0)
  
  if ("Biosample_term_name" %in% names(df)) {
    biosamples <- sort(unique(as.character(df$Biosample_term_name)))
    biosamples <- biosamples[nzchar(biosamples) & !is.na(biosamples)]
  }
  
  if ("Life_stage" %in% names(df)) {
    lifestages <- sort(unique(as.character(df$Life_stage)))
    lifestages <- lifestages[nzchar(lifestages) & !is.na(lifestages)]
  }
  
  # Add file basename for easier matching
  if ("File" %in% names(df)) {
    df$file_base <- vapply(df$File, function(x) tools::file_path_sans_ext(basename(as.character(x))), character(1))
  }
  
  return(list(biosamples = biosamples, lifestages = lifestages, full_data = df))
}

# Helper to save filtered clustering results to RDS file
save_filtered_results_to_rds <- function(source_file, filter_key, filtered_results, dir_path) {
  candidate_path <- file.path(dir_path, source_file)
  if (!file.exists(candidate_path)) {
    hits <- Sys.glob(file.path(dir_path, "**", source_file))
    if (length(hits) >= 1) candidate_path <- hits[1]
  }
  
  if (!file.exists(candidate_path)) return(FALSE)
  
  cache_path <- file.path(dirname(candidate_path), paste0(basename(candidate_path), ".prep.rds"))
  
  # Load existing cache
  cached <- NULL
  if (file.exists(cache_path)) {
    tryCatch({
      cached <- readRDS(cache_path)
    }, error = function(e) {
      cached <<- NULL
    })
  }
  
  if (is.null(cached)) return(FALSE)
  
  # Initialize filtered_results if not present
  if (is.null(cached$filtered_results)) {
    cached$filtered_results <- list()
  }
  
  # Add the new filtered results
  cached$filtered_results[[filter_key]] <- filtered_results
  cached$cache_version <- "v3"  # Update version
  
  # Save back to file
  tryCatch({
    saveRDS(cached, cache_path)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

# Helper to load filtered clustering results from RDS file
load_filtered_results_from_rds <- function(source_file, filter_key, dir_path) {
  candidate_path <- file.path(dir_path, source_file)
  if (!file.exists(candidate_path)) {
    hits <- Sys.glob(file.path(dir_path, "**", source_file))
    if (length(hits) >= 1) candidate_path <- hits[1]
  }
  
  if (!file.exists(candidate_path)) return(NULL)
  
  cache_path <- file.path(dirname(candidate_path), paste0(basename(candidate_path), ".prep.rds"))
  
  if (!file.exists(cache_path)) return(NULL)
  
  # Load cache
  cached <- NULL
  tryCatch({
    cached <- readRDS(cache_path)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(cached) || is.null(cached$filtered_results)) return(NULL)
  
  # Return the specific filtered results if they exist
  if (filter_key %in% names(cached$filtered_results)) {
    return(cached$filtered_results[[filter_key]])
  }
  
  return(NULL)
}

# Helper to filter similarity matrix based on accession key criteria
filter_similarity_matrix <- function(sim_matrix, accession_meta, selected_biosamples, selected_lifestage) {
  if (is.null(sim_matrix) || is.null(accession_meta) || is.null(accession_meta$full_data)) {
    return(sim_matrix)
  }
  
  # If no filtering, return original matrix
  if ((is.null(selected_biosamples) || length(selected_biosamples) == 0) && 
      (is.null(selected_lifestage) || selected_lifestage == "ALL")) {
    return(sim_matrix)
  }
  
  # Add debugging for DNase-hypersensitivity dataset
  cat(sprintf("DEBUG: Filtering matrix with %d samples\n", nrow(sim_matrix)))
  cat(sprintf("DEBUG: Selected biosamples: %s\n", paste(selected_biosamples, collapse = ", ")))
  cat(sprintf("DEBUG: Selected lifestage: %s\n", selected_lifestage))
  
  acc_df <- accession_meta$full_data
  if (!"file_base" %in% names(acc_df)) {
    return(sim_matrix)
  }
  
  # Filter accession data by selections
  filtered_acc <- acc_df
  if (!is.null(selected_biosamples) && length(selected_biosamples) > 0) {
    if ("Biosample_term_name" %in% names(filtered_acc)) {
      filtered_acc <- filtered_acc[filtered_acc$Biosample_term_name %in% selected_biosamples, , drop = FALSE]
    }
  }
  
  if (!is.null(selected_lifestage) && selected_lifestage != "ALL" && nzchar(selected_lifestage)) {
    if ("Life_stage" %in% names(filtered_acc)) {
      filtered_acc <- filtered_acc[filtered_acc$Life_stage == selected_lifestage, , drop = FALSE]
    }
  }
  
  # Get the samples that should be included
  included_samples <- unique(filtered_acc$file_base)
  cat(sprintf("DEBUG: Included samples from accession key: %d\n", length(included_samples)))
  
  # Filter similarity matrix to only include these samples
  sample_names <- rownames(sim_matrix)
  valid_samples <- intersect(sample_names, included_samples)
  cat(sprintf("DEBUG: Valid samples after intersection: %d\n", length(valid_samples)))
  cat(sprintf("DEBUG: Sample names in matrix: %s\n", paste(head(sample_names, 3), collapse = ", ")))
  cat(sprintf("DEBUG: Included samples from key: %s\n", paste(head(included_samples, 3), collapse = ", ")))
  
  if (length(valid_samples) < 2) {
    # Not enough samples for clustering
    return(NULL)
  }
  
  # Subset the matrix
  filtered_matrix <- sim_matrix[valid_samples, valid_samples, drop = FALSE]
  return(filtered_matrix)
}

# Recompute clustering analysis on filtered similarity matrix
recompute_clustering_on_filtered_matrix <- function(filtered_sim_matrix, labels_map, clusters_vec = 2:30, linkage_methods = c("average", "complete")) {
  if (is.null(filtered_sim_matrix) || nrow(filtered_sim_matrix) < 2) {
    cat("Warning: Filtered matrix is NULL or has < 2 rows\n")
    return(NULL)
  }
  
  cat(sprintf("Recomputing clustering on %dx%d filtered matrix\n", nrow(filtered_sim_matrix), ncol(filtered_sim_matrix)))
  cat(sprintf("Sample names: %s\n", paste(rownames(filtered_sim_matrix)[1:min(5, nrow(filtered_sim_matrix))], collapse = ", ")))
  cat(sprintf("Similarity range: %.3f to %.3f\n", min(filtered_sim_matrix, na.rm = TRUE), max(filtered_sim_matrix, na.rm = TRUE)))
  
  # Check if labels_map has entries for our samples
  sample_names <- rownames(filtered_sim_matrix)
  matched_labels <- sum(sample_names %in% names(labels_map))
  cat(sprintf("Labels matched: %d/%d samples\n", matched_labels, length(sample_names)))
  
  # Check the similarity matrix before passing to evaluate_nmi_long_table
  sim_df <- as.data.frame(filtered_sim_matrix)
  cat(sprintf("Similarity matrix diagnostics:\n"))
  cat(sprintf("  Dimensions: %dx%d\n", nrow(sim_df), ncol(sim_df)))
  cat(sprintf("  Value range: %.6f to %.6f\n", min(sim_df, na.rm = TRUE), max(sim_df, na.rm = TRUE)))
  cat(sprintf("  Any NA values: %s\n", any(is.na(sim_df))))
  cat(sprintf("  Any infinite values: %s\n", any(is.infinite(as.matrix(sim_df)))))
  cat(sprintf("  Any values > 1: %s\n", any(sim_df > 1, na.rm = TRUE)))
  cat(sprintf("  Any values < 0: %s\n", any(sim_df < 0, na.rm = TRUE)))
  cat(sprintf("  Diagonal values: %s\n", paste(round(diag(as.matrix(sim_df))[1:min(3, nrow(sim_df))], 3), collapse = ", ")))
  
  # Check for clustering-problematic conditions
  sim_matrix_vals <- as.matrix(sim_df)
  variance_check <- var(sim_matrix_vals[upper.tri(sim_matrix_vals)], na.rm = TRUE)
  cat(sprintf("  Similarity variance (upper triangle): %.6f\n", variance_check))
  
  # Check label diversity
  sample_labels <- labels_map[sample_names]
  sample_labels[is.na(sample_labels)] <- "unknown"
  unique_labels <- length(unique(sample_labels))
  cat(sprintf("  Unique biological labels: %d\n", unique_labels))
  cat(sprintf("  Label distribution: %s\n", paste(names(table(sample_labels)), collapse = ", ")))
  
  # Warn about potential issues
  if (nrow(sim_df) < 4) {
    cat("WARNING: Very small matrix (< 4 samples) - silhouette may be unreliable\n")
  }
  if (variance_check < 1e-6) {
    cat("WARNING: Very low similarity variance - clustering may not be meaningful\n")
  }
  if (unique_labels < 2) {
    cat("WARNING: Only one unique biological label - NMI will be undefined\n")
  }
  
  # Use the evaluate_nmi_long_table function from clustering_analysis.R
  tryCatch({
    result <- evaluate_nmi_long_table(sim_df, labels_map, clusters_vec, linkage_methods)
    
    # Check the result
    if (!is.null(result) && nrow(result) > 0) {
      na_silhouette <- sum(is.na(result$silhouette))
      cat(sprintf("Clustering result: %d rows, %d NA silhouette scores\n", nrow(result), na_silhouette))
      if (na_silhouette > 0) {
        cat("Sample of NA silhouette rows:\n")
        na_rows <- result[is.na(result$silhouette), ]
        print(head(na_rows[, c("n_clusters", "linkage", "nmi", "silhouette")], 3))
      }
    }
    
    return(result)
  }, error = function(e) {
    cat(sprintf("Error in recompute_clustering: %s\n", conditionMessage(e)))
    return(NULL)
  })
}

# Run Python or R clustering on a directory of hammock outputs to produce long-tables in-memory
run_clustering_on_dir <- function(dir_path, acc_key, files = NULL, use_cache = FALSE, progress_cb = NULL) {
  diag <- character(0)
  if (is.null(files)) files <- list.files(dir_path, pattern = "\\.(csv|tsv|txt)$", ignore.case = TRUE, full.names = TRUE)
  diag <- c(diag, sprintf("dir: %s", dir_path))
  diag <- c(diag, sprintf("accession_key: %s", acc_key))
  diag <- c(diag, sprintf("files_scanned: %d", length(files)))
  if (length(files) == 0) return(list(df = NULL, diag = c(diag, "No files found in directory.")))
  if (is.null(acc_key) || !nzchar(acc_key) || !file.exists(acc_key)) return(list(df = NULL, diag = c(diag, "Accession key missing or not found.")))
  # Load accession key once for overlap diagnostics
  labels_map_diag <- NULL
  tryCatch({ labels_map_diag <- load_accession_key(acc_key) }, error = function(e) {
    diag <- c(diag, sprintf("accession_key_error: %s", conditionMessage(e)))
    labels_map_diag <- NULL
  })
  # Prefer R implementation for in-process results
  res_list <- list()
  # Source once: try multiple candidate locations relative to the app directory
  app_dir <- getwd()
  base_dirs <- unique(c(app_dir, dirname(app_dir), dirname(dirname(app_dir)), dirname(dirname(dirname(app_dir)))))
  cand_paths <- file.path(base_dirs, "scripts", "clustering_analysis.R")
  sourced_ok <- FALSE
  for (p in cand_paths) {
    if (file.exists(p)) { try({ source(p); sourced_ok <- TRUE; diag <- c(diag, sprintf("sourced: %s", p)) }, silent = TRUE); if (sourced_ok) break }
  }
  if (!sourced_ok && !exists("clustering_analysis")) return(list(df = NULL, diag = c(diag, "Failed to source scripts/clustering_analysis.R")))
  # Also ensure utils.R is sourced so helper functions are available when this file is sourced (not run via Rscript)
  utils_ok <- exists("detect_file_format")
  if (!utils_ok) {
    utils_cands <- file.path(base_dirs, "scripts", "utils.R")
    for (u in utils_cands) {
      if (file.exists(u)) { try({ source(u); utils_ok <- TRUE; diag <- c(diag, sprintf("sourced: %s", u)) }, silent = TRUE); if (utils_ok) break }
    }
  }
  # Normalizer to ensure uniform columns/order
  required_cols <- c('sketch','mode','precision','expA','subA','subB','window','klen',
                     'n_clusters','linkage','nmi','silhouette','source_file')
  total_files <- length(files)
  file_index <- 0
  for (f in files) {
    file_index <- file_index + 1
    diag <- c(diag, sprintf("scan: %d/%d %s", file_index, total_files, basename(f)))
    if (is.function(progress_cb)) try(progress_cb(0, sprintf("Scanning %d/%d: %s", file_index, total_files, basename(f))), silent = TRUE)
    # Decide if hammock or bedtools based on header; let the R function decide
    out_tmp <- tempfile(fileext = ".tsv")
    tbl <- NULL
    err <- NULL
    # Determine cluster upper bound up to 100 based on matrix size
    clusters_vec <- 2:30
    # Try to compute size
    fmt_f <- NULL
    n_files <- NA_integer_
    try({
      fmt_f <- detect_file_format(f)
      sim_df_quick <- if (fmt_f == 'hammock') parse_hammock_format(f) else parse_bedtools_format(f)
      n_files <- nrow(sim_df_quick)
      upper <- max(2, min(100, n_files))
      clusters_vec <- 2:upper
    }, silent = TRUE)
    # Quick meta parse and report expA
    meta_quick <- tryCatch(extract_params_from_filename(f), error = function(e) NULL)
    if (!is.null(meta_quick)) {
      exp_str <- if (!is.null(meta_quick$expA) && is.finite(meta_quick$expA)) sprintf("%.2f", as.numeric(meta_quick$expA)) else "NA"
      diag <- c(diag, sprintf("expA_parse: %s => %s", basename(f), exp_str))
    }
    # Overlap diagnostic (how many filenames match the accession key)
    if (!is.null(labels_map_diag) && !is.null(sim_df_quick)) {
      sample_names <- rownames(sim_df_quick)
      ov <- sum(sample_names %in% names(labels_map_diag))
      diag <- c(diag, sprintf("overlap: %s => %d/%d", basename(f), ov, length(sample_names)))
    }
    # Cache path (same directory as input file)
    base_name <- basename(f)
    cache_path <- file.path(dirname(f), paste0(base_name, ".prep.rds"))
    used_cache <- FALSE
    if (isTRUE(use_cache) && file.exists(cache_path)) {
      used_cache <- TRUE
      diag <- c(diag, sprintf("cache_load: %s", cache_path))
      cached <- NULL
      tryCatch({ 
        cached <- readRDS(cache_path) 
      }, error = function(e) { 
        err <<- conditionMessage(e) 
      })
      if (!is.null(cached) && is.list(cached) && (!is.null(cached$long_table_raw) || !is.null(cached$similarity))) {
        # If long-table cached with same key, reuse directly
        key_same <- FALSE
        if (!is.null(cached$key_path) && nzchar(cached$key_path) && file.exists(cached$key_path)) {
          key_same <- tryCatch({
            normalizePath(cached$key_path) == normalizePath(acc_key)
          }, error = function(e) FALSE)
        }
        if (key_same && !is.null(cached$long_table_raw) && is.data.frame(cached$long_table_raw)) {
          lt <- cached$long_table_raw
        } else if (!is.null(cached$similarity)) {
          # Recompute from cached similarity
          key_to_use <- if (!is.null(cached$key_path) && nzchar(cached$key_path) && file.exists(cached$key_path)) cached$key_path else acc_key
          if (!identical(key_to_use, acc_key)) diag <- c(diag, sprintf("cache_key: %s", key_to_use))
          labels_map <- load_accession_key(key_to_use)
          lt <- evaluate_nmi_long_table(as.data.frame(cached$similarity), labels_map, clusters_vec, c("average","complete"))
        } else {
          lt <- NULL
        }
        if (!is.null(lt)) {
          # Recompute meta from filename to ensure latest parsing (avoids stale cached meta)
          meta <- extract_params_from_filename(f)
          if (!is.null(meta$mode) && toupper(meta$mode) == 'C') {
            # Try to detect missing parameters from file header if not present in name
            if (is.null(meta$expA) || is.na(meta$expA) || !is.finite(meta$expA)) {
              meta$expA <- tryCatch(detect_hammock_expA(f), error = function(e) NA_real_)
            }
            if (is.null(meta$subA) || is.na(meta$subA) || !is.finite(meta$subA)) {
              meta$subA <- tryCatch(detect_hammock_subA(f), error = function(e) NA_real_)
            }
            if (is.null(meta$subB) || is.na(meta$subB) || !is.finite(meta$subB)) {
              meta$subB <- tryCatch(detect_hammock_subB(f), error = function(e) NA_real_)
            }
          }
          # Safely add expected columns; fill missing with NA
          lt$sketch <- meta$sketch
          lt$mode <- meta$mode
          lt$precision <- meta$precision
          lt$expA <- meta$expA
          lt$subA <- meta$subA
          lt$subB <- meta$subB
          lt$window <- meta$window
          lt$klen <- meta$klen
          lt$source_file <- basename(f)
          for (mc in setdiff(required_cols, names(lt))) lt[[mc]] <- NA
          lt <- lt[, required_cols, drop = FALSE]
          res_list[[length(res_list) + 1]] <- lt
          next
        }
      }
    }
    diag <- c(diag, sprintf("cluster: %d/%d %s", file_index, total_files, basename(f)))
    if (is.function(progress_cb)) try(progress_cb(0, sprintf("Clustering %d/%d: %s", file_index, total_files, basename(f))), silent = TRUE)
    tryCatch({ 
      tbl <- clustering_analysis(f, acc_key, clusters = clusters_vec, out = out_tmp) 
    }, error = function(e) { 
      err <<- conditionMessage(e) 
    })
    if (!is.null(tbl) && is.data.frame(tbl)) {
      # Add source and normalize
      tbl$source_file <- basename(f)
      # Fill missing required cols with NA
      missing_cols <- setdiff(required_cols, names(tbl))
      for (mc in missing_cols) tbl[[mc]] <- NA
      # Drop extra columns
      keep <- intersect(required_cols, names(tbl))
      tbl <- tbl[, keep]
      # Reorder
      tbl <- tbl[, required_cols]
      # Standardize types for key fields
      if ("mode" %in% names(tbl)) {
        tbl$mode <- toupper(substr(as.character(tbl$mode), 1, 1))
        tbl$mode[is.na(tbl$mode) | !(tbl$mode %in% c("A","B","C","D"))] <- "BEDTOOLS"
      }
      if ("linkage" %in% names(tbl)) tbl$linkage <- tolower(as.character(tbl$linkage))
      res_list[[length(res_list) + 1]] <- tbl
       diag <- c(diag, sprintf("ok: %s", basename(f)))
       if (is.function(progress_cb)) {
         try(progress_cb(1/total_files, sprintf("Done %d/%d: %s", file_index, total_files, basename(f))), silent = TRUE)
       }
      # Save cache if requested and not loaded
      if (!used_cache && isTRUE(use_cache)) {
        # Build similarity for cache
        sim_df_for_cache <- NULL
        fmt_tmp <- NA_character_
        try({
          fmt_tmp <- detect_file_format(f)
          sim_df_for_cache <- if (fmt_tmp == 'hammock') parse_hammock_format(f) else parse_bedtools_format(f)
        }, silent = TRUE)
        if (!is.null(sim_df_for_cache)) {
          meta <- extract_params_from_filename(f)
          # Extract accession metadata and sample names for more efficient filtering
          acc_meta <- extract_accession_metadata(acc_key)
          sample_names <- rownames(sim_df_for_cache)
          
          # Create sample metadata lookup for this file
          sample_metadata <- NULL
          if (!is.null(acc_meta$full_data) && "file_base" %in% names(acc_meta$full_data)) {
            acc_df <- acc_meta$full_data
            sample_meta_list <- list()
            for (sname in sample_names) {
              matches <- acc_df[acc_df$file_base == sname, , drop = FALSE]
              if (nrow(matches) > 0) {
                sample_meta_list[[sname]] <- list(
                  biosample = if ("Biosample_term_name" %in% names(matches)) as.character(matches$Biosample_term_name[1]) else "unknown",
                  lifestage = if ("Life_stage" %in% names(matches)) as.character(matches$Life_stage[1]) else "unknown",
                  organism = if ("Organism" %in% names(matches)) as.character(matches$Organism[1]) else "unknown"
                )
              } else {
                sample_meta_list[[sname]] <- list(biosample = "unknown", lifestage = "unknown", organism = "unknown")
              }
            }
            sample_metadata <- sample_meta_list
          }
          
          obj <- list(
            similarity = as.matrix(sim_df_for_cache), 
            meta = meta, 
            source_file = basename(f), 
            key_path = acc_key, 
            long_table_raw = tbl,
            accession_metadata = acc_meta,
            sample_metadata = sample_metadata,
            filtered_results = list(),  # Store filtered clustering results by filter key
            cache_version = "v3"  # Version to track cache format changes
          )
          try({ saveRDS(obj, cache_path); diag <- c(diag, sprintf("cache_save: %s", cache_path)) }, silent = TRUE)
        }
      }
    }
    else {
      diag <- c(diag, sprintf("fail: %s :: %s", basename(f), ifelse(is.null(err), "unknown error", err)))
      if (is.function(progress_cb)) {
        try(progress_cb(1/total_files, sprintf("Failed %d/%d: %s", file_index, total_files, basename(f))), silent = TRUE)
      }
    }
  }
  if (length(res_list) == 0) return(list(df = NULL, diag = c(diag, "No tables produced from any files.")))
  # Bind safely (all have same columns now)
  # Ensure all have identical columns before binding
  cols_ref <- required_cols
  res_list <- lapply(res_list, function(x) { x[, cols_ref, drop = FALSE] })
  df <- do.call(rbind, res_list)
  # Coerce numeric where appropriate
  num_cols <- c("precision", "expA", "subA", "subB", "window", "klen", "n_clusters", "nmi", "silhouette")
  for (col in intersect(names(df), num_cols)) {
    suppressWarnings(df[[col]] <- as.numeric(df[[col]]))
  }
  

  
  # Return both the data frame and any loaded cache data for efficiency
  loaded_cache <- list()
  if (isTRUE(use_cache)) {
    # Collect any cache data that was loaded during processing
    for (f in files) {
      base_name <- basename(f)
      cache_path <- file.path(dirname(f), paste0(base_name, ".prep.rds"))
      if (file.exists(cache_path)) {
        tryCatch({
          cached <- readRDS(cache_path)
          if (!is.null(cached$similarity)) {
            loaded_cache[[base_name]] <- cached$similarity
          }
        }, error = function(e) {
          # Silently ignore cache read errors
        })
      }
    }
  }
  
  list(df = df, diag = c(diag, sprintf("rows: %d", nrow(df))), loaded_cache = loaded_cache)
}

ui <- fluidPage(
  titlePanel("Cluster Analysis Explorer (on hammock outputs)"),
  tabsetPanel(
    tabPanel(
      "Inputs",
      # Filter summary text box at the top
      fluidRow(
        column(12,
          div(
            style = "background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; padding: 10px; margin-bottom: 15px;",
            h5("Active Filters Summary", style = "margin-top: 0; margin-bottom: 8px; color: #495057; font-weight: bold;"),
            verbatimTextOutput("filterSummaryInputs", placeholder = TRUE)
          )
        )
      ),
      fluidRow(
        column(4,
          h4("File Selection"),
          shinyDirButton("dir", "Select directory of hammock outputs", "Select"),
          shinyFilesButton("acc", "Select accession key TSV", "Select", multiple = FALSE),
          textInput("dirText", "Or paste outputs directory", value = ""),
          textInput("accText", "Or paste accession key TSV", value = ""),
          br(),
          checkboxInput("useCache", "Use cached RDS if available (same folder as inputs)", value = TRUE)
        ),
        column(4,
          h4("Life Stage Filtering"),
          radioButtons("lifestageSelect", "Life Stage:", 
                      choices = c("ALL" = "ALL"), 
                      selected = "ALL",
                      inline = TRUE)
        ),
        column(4,
          h4("Biosample Filtering"),
          radioButtons("biosampleMode", "Biosample Selection:", 
                      choices = c("All biosamples" = "ALL", "Custom subset" = "SUBSET"), 
                      selected = "ALL",
                      inline = TRUE),
          conditionalPanel(
            condition = "input.biosampleMode == 'REGULAR' || input.biosampleMode == 'SUBSET'",
            div(
              style = "border: 1px solid #ddd; padding: 10px; border-radius: 5px; background-color: #f9f9f9;",
              div(
                style = "column-count: 2; column-gap: 10px;",
                checkboxGroupInput("biosampleCheckboxes", "Select Biosamples:", 
                                 choices = character(0), 
                                 selected = character(0))
              ),
              br(),
              fluidRow(
                column(8,
                  textInput("selectedBiosamples", "Selected Biosamples:", 
                           value = "", 
                           placeholder = "No biosamples selected",
                           width = "100%")
                ),
                column(4,
                  actionButton("applyBiosampleFilter", "Apply Filter", 
                              class = "btn-primary", 
                              style = "margin-top: 25px; width: 100%;")
                )
              )
            )
          )
        )
      ),
      tags$hr(),
      fluidRow(
        column(4,
          tags$h4("Data Summary by Biosample & Life Stage"),
          tags$p("Shows the distribution of samples across biological conditions.", 
                 style = "font-size: 12px; color: #666; margin-bottom: 10px;"),
          tableOutput("biologicalSummaryTable")
        ),
        column(4,
          tags$h4("Filtered Data Summary"),
          tags$p("Shows the distribution of samples based on current biological filter selections.", 
                 style = "font-size: 12px; color: #666; margin-bottom: 10px;"),
          tableOutput("filteredDataSummaryTable")
        ),
        column(4,
          tags$h4("Top NMI Results per File"),
          fluidRow(
            column(8,
              textInput("downloadFilename", "Filename:", 
                       value = "top_nmi_results", 
                       placeholder = "Enter filename (without extension)",
                       width = "100%")
            ),
            column(4,
              br(),
              downloadButton("downloadSummaryTable", "Download CSV", 
                           class = "btn-primary", 
                           style = "width: 100%;")
            )
          ),
          br(),
          tableOutput("summaryTable")
        )
      ),
      tags$hr(),
      fluidRow(
        column(6,
      tags$div(
        style = "background-color:#FFF8C6;border:1px solid #E0C97F;padding:10px;border-radius:6px;",
        tags$h4("Diagnostics", style = "margin-top:0;font-weight:bold;color:#8A6D3B;"),
        verbatimTextOutput("diagText")
          )
        ),
        column(6,
          tags$div(
            style = "background-color:#E8F5E8;border:1px solid #4CAF50;padding:10px;border-radius:6px;",
            tags$h4("Cache Status", style = "margin-top:0;font-weight:bold;color:#2E7D32;"),
            verbatimTextOutput("cacheStatus")
          )
        )
      )
    ),
    tabPanel(
      "Results",
      # Filter summary text box at the top
      fluidRow(
        column(12,
          div(
            style = "background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; padding: 10px; margin-bottom: 15px;",
            h5("Active Filters Summary", style = "margin-top: 0; margin-bottom: 8px; color: #495057; font-weight: bold;"),
            verbatimTextOutput("filterSummaryResults", placeholder = TRUE)
          )
        )
      ),
      fluidRow(
        column(6,
          radioButtons("modeSelect", "Mode", choices = c("A", "B", "C", "D", "BEDTOOLS"), inline = TRUE),
          conditionalPanel(
            condition = "input.modeSelect != 'BEDTOOLS'",
            selectInput("precision", "Precision (p)", choices = character(0))
          ),
          selectInput("linkageMethod", "Clustering approach (linkage)", choices = character(0)),
          conditionalPanel(
            condition = "input.modeSelect == 'C'",
            selectInput("expA", "expA", choices = character(0), selected = "0"),
            selectInput("subB", "subB", choices = character(0), selected = "1")
          ),
          conditionalPanel(
            condition = "input.modeSelect == 'D'",
            selectInput("klen", "k-mer size (k)", choices = character(0)),
            selectInput("window", "window size (w)", choices = character(0))
          ),
          conditionalPanel(
            condition = "input.modeSelect == 'BEDTOOLS'",
            div(
              style = "background-color: #e8f4fd; border: 1px solid #bee5eb; border-radius: 5px; padding: 15px; margin-top: 10px;",
              h5("BEDTOOLS Mode", style = "color: #0c5460; margin-top: 0;"),
              p("Using bedtools pairwise Jaccard similarity data. No additional parameters needed.", 
                style = "color: #0c5460; margin-bottom: 0; font-size: 14px;")
            )
          )
        ),
        column(6,
          h4("Dendrogram controls"),
          selectInput("sourceFile", "Input file", choices = character(0)),
          selectInput("dendLinkage", "Dendrogram linkage", choices = c("average","complete","single","ward"))
        )
      ),
      tags$hr(),
      verbatimTextOutput("statusText"),
      tags$hr(),
      fluidRow(
        column(6, plotlyOutput("nmiPlot", height = "520px")),
        column(6, plotOutput("dendPlot", height = "520px"))
      ),
      tags$hr()
    )
    ,
    tabPanel(
      "Heatmaps",
      # Filter summary text box at the top
      fluidRow(
        column(12,
          div(
            style = "background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; padding: 10px; margin-bottom: 15px;",
            h5("Active Filters Summary", style = "margin-top: 0; margin-bottom: 8px; color: #495057; font-weight: bold;"),
            verbatimTextOutput("filterSummaryHeatmaps", placeholder = TRUE)
          )
        )
      ),
      fluidRow(
        column(6,
          selectInput("heatmapFile1", "Heatmap A input file", choices = character(0)),
          plotOutput("heatmap1", height = "600px")
        ),
        column(6,
          selectInput("heatmapFile2", "Heatmap B input file", choices = character(0)),
          plotOutput("heatmap2", height = "600px")
        )
      ),
      tags$hr(),
      tags$div(
        tags$h4("Heatmap details"),
        verbatimTextOutput("heatmapInfo")
      )
    ),
    tabPanel(
      "Dendrograms",
      # Filter summary text box at the top
      fluidRow(
        column(12,
          div(
            style = "background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; padding: 10px; margin-bottom: 15px;",
            h5("Active Filters Summary", style = "margin-top: 0; margin-bottom: 8px; color: #495057; font-weight: bold;"),
            verbatimTextOutput("filterSummaryDendrograms", placeholder = TRUE)
          )
        )
      ),
      fluidRow(
        column(6,
          h4("Dendrogram A"),
          selectInput("dendFile1", "Input file", choices = character(0)),
          selectInput("dendLinkage1", "Clustering approach", 
                     choices = c("average", "complete", "single", "ward"), 
                     selected = "average"),
          plotOutput("dendPlotA", height = "520px")
        ),
        column(6,
          h4("Dendrogram B"),
          selectInput("dendFile2", "Input file", choices = character(0)),
          selectInput("dendLinkage2", "Clustering approach", 
                     choices = c("average", "complete", "single", "ward"), 
                     selected = "complete"),
          plotOutput("dendPlotB", height = "520px")
        )
      )
    ),
    tabPanel(
      "PCA",
      # Filter summary text box at the top
      fluidRow(
        column(12,
          div(
            style = "background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; padding: 10px; margin-bottom: 15px;",
            h5("Active Filters Summary", style = "margin-top: 0; margin-bottom: 8px; color: #495057; font-weight: bold;"),
            verbatimTextOutput("filterSummaryPCA", placeholder = TRUE)
          )
        )
      ),
      fluidRow(
        column(6,
          h4("PCA Plot 1"),
          selectInput("pcaFile1", "Select File:", 
                     choices = character(0), 
                     selected = character(0),
                     width = "100%"),
          plotlyOutput("pcaPlot1", height = "500px")
        ),
        column(6,
          h4("PCA Plot 2"),
          selectInput("pcaFile2", "Select File:", 
                     choices = character(0), 
                     selected = character(0),
                     width = "100%"),
          plotlyOutput("pcaPlot2", height = "500px")
        )
      ),
      br(),
      fluidRow(
        column(12,
          h4("PCA Information"),
          verbatimTextOutput("pcaInfo")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  data_all <- reactiveVal(NULL)
  accession_metadata <- reactiveVal(NULL)
  cached_similarity_matrices <- reactiveVal(NULL)  # Store similarity matrices for recalculation
  clustering_results_cache <- reactiveVal(list())  # Cache computed clustering results by filter combination
  approved_biosamples <- reactiveVal(NULL)  # Store approved biosample selection
  
  # Default diagnostics panel content
  output$diagText <- renderText("(no diagnostics yet)")
  
  # Cache status display
  output$cacheStatus <- renderText({
    # Force reactive invalidation by accessing the reactive values
    sim_cache <- cached_similarity_matrices()
    results_cache <- clustering_results_cache()
    # Also access accession metadata to trigger updates during loading
    acc_meta <- accession_metadata()
    
    if (is.null(sim_cache) && is.null(results_cache)) {
      return("No data loaded")
    }
    
    sim_files <- if (is.null(sim_cache)) 0 else length(sim_cache)
    session_cached_combinations <- if (is.null(results_cache)) 0 else length(results_cache)
    
    # Simple loading indicator - only show if we have some files but not many
    loading_status <- ""
    if (sim_files > 0 && sim_files < 10) {
      loading_status <- sprintf("\n🔄 Loading: %d files cached", sim_files)
    }
    
    # Check RDS files for persistent cache info
    dir_path <- effective_dir()
    rds_cache_info <- ""
    if (!is.na(dir_path) && dir.exists(dir_path) && !is.null(sim_cache)) {
      total_rds_combinations <- 0
      rds_details <- c()
      
      for (src_file in names(sim_cache)) {
        candidate_path <- file.path(dir_path, src_file)
        if (!file.exists(candidate_path)) {
          hits <- Sys.glob(file.path(dir_path, "**", src_file))
          if (length(hits) >= 1) candidate_path <- hits[1]
        }
        
        if (file.exists(candidate_path)) {
          cache_path <- file.path(dirname(candidate_path), paste0(basename(candidate_path), ".prep.rds"))
          
          if (file.exists(cache_path)) {
            cached <- tryCatch({
              readRDS(cache_path)
            }, error = function(e) NULL)
            
            if (!is.null(cached) && !is.null(cached$filtered_results)) {
              file_combinations <- length(cached$filtered_results)
              total_rds_combinations <- total_rds_combinations + file_combinations
              if (file_combinations > 0) {
                rds_details <- c(rds_details, sprintf("  %s: %d combinations", src_file, file_combinations))
              }
            }
          }
        }
      }
      
      rds_cache_info <- paste(
        sprintf("RDS persistent cache: %d combinations across all files", total_rds_combinations),
        if (length(rds_details) > 0) paste(rds_details, collapse = "\n") else "  No RDS cached combinations yet",
        sep = "\n"
      )
    }
    
    paste(
      sprintf("✓ Similarity matrices cached: %d files", sim_files),
      loading_status,
      sprintf("✓ Session cache: %d filter combinations", session_cached_combinations),
      if (session_cached_combinations > 0) {
        paste("Session combinations:", paste(names(results_cache), collapse = ", "))
      } else {
        "No session combinations yet"
      },
      rds_cache_info,
      if (isTRUE(input$useCache)) {
        "\n💾 Cache enabled - faster loading"
      } else {
        "\n⚠️ Cache disabled - slower loading"
      },
      sep = "\n"
    )
  })
  
  # Helper to create cache key for filter combination
  create_filter_cache_key <- function(selected_biosamples, selected_lifestage) {
    biosamples_key <- if (is.null(selected_biosamples) || length(selected_biosamples) == 0) {
      "ALL_BIOSAMPLES"
    } else {
      paste(sort(selected_biosamples), collapse = "|")
    }
    
    lifestage_key <- if (is.null(selected_lifestage) || selected_lifestage == "ALL") {
      "ALL_LIFESTAGES"
    } else {
      selected_lifestage
    }
    
    paste(biosamples_key, lifestage_key, sep = ":::")
  }
  
  # Observer to update the selected biosamples text field when checkboxes change
  observeEvent(input$biosampleCheckboxes, {
    if (!is.null(input$biosampleCheckboxes) && length(input$biosampleCheckboxes) > 0) {
      updateTextInput(session, "selectedBiosamples", 
                     value = paste(input$biosampleCheckboxes, collapse = ", "))
    } else {
      updateTextInput(session, "selectedBiosamples", 
                     value = "",
                     placeholder = "No biosamples selected")
    }
  }, ignoreInit = TRUE)
  
  # Observer to handle apply button click
  observeEvent(input$applyBiosampleFilter, {
    if (input$biosampleMode == "SUBSET") {
      if (!is.null(input$biosampleCheckboxes) && length(input$biosampleCheckboxes) > 0) {
        approved_biosamples(input$biosampleCheckboxes)
        showNotification(
          paste("Applied filter for", length(input$biosampleCheckboxes), "biosamples:", 
                paste(input$biosampleCheckboxes, collapse = ", ")),
          type = "message",
          duration = 3
        )
      } else {
        showNotification("Please select at least one biosample before applying the filter.", 
                        type = "warning", duration = 3)
      }
    }
  })
  
  # Reactive to determine biosample mode choices based on available data
  biosample_mode_choices <- reactive({
    acc_meta <- accession_metadata()
    if (!is.null(acc_meta) && !is.null(acc_meta$biosamples)) {
      available_biosamples <- acc_meta$biosamples
      
      # Define the regular subset tissues
      regular_tissues <- c("kidney", "lung", "heart", "brain", "large intestine", "liver", "stomach")
      
      # Check if any regular tissues are present (exact matches only, not prefixed)
      has_regular_tissues <- FALSE
      for (tissue in regular_tissues) {
        # Use word boundaries to ensure exact matches, not prefixed matches
        pattern <- paste0("\\b", tissue, "\\b")
        if (any(grepl(pattern, available_biosamples, ignore.case = TRUE))) {
          has_regular_tissues <- TRUE
          break
        }
      }
      
      if (has_regular_tissues) {
        c("All biosamples" = "ALL", "Regular subset" = "REGULAR", "Custom subset" = "SUBSET")
      } else {
        c("All biosamples" = "ALL", "Custom subset" = "SUBSET")
      }
    } else {
      c("All biosamples" = "ALL", "Custom subset" = "SUBSET")
    }
  })
  
  # Observer to update biosample mode choices based on available data
  observeEvent(biosample_mode_choices(), {
    choices <- biosample_mode_choices()
    cat("Biosample mode choices updated:", paste(names(choices), collapse = ", "), "\n")
    
    # Determine the new selection
    current_mode <- input$biosampleMode
    new_selection <- if (current_mode == "REGULAR" && !("REGULAR" %in% choices)) "ALL" else current_mode
    
    cat("Current mode:", current_mode, "New selection:", new_selection, "\n")
    
    # Update the radio buttons
    updateRadioButtons(session, "biosampleMode", 
                      choices = choices,
                      selected = new_selection)
    
    # Clear approved biosamples if switching from REGULAR
    if (current_mode == "REGULAR" && new_selection != "REGULAR") {
      approved_biosamples(NULL)
    }
  })
  
  # Observer to reset approved biosamples when switching back to "ALL"
  observeEvent(list(input$biosampleMode, accession_metadata()), {
    # Safety check: if mode is REGULAR but no regular tissues are available, force switch to ALL
    if (input$biosampleMode == "REGULAR") {
      acc_meta <- accession_metadata()
      if (!is.null(acc_meta) && !is.null(acc_meta$biosamples)) {
        available_biosamples <- acc_meta$biosamples
        regular_tissues <- c("kidney", "lung", "heart", "brain", "large intestine", "liver", "stomach")
        
        has_regular_tissues <- FALSE
        for (tissue in regular_tissues) {
          # Use word boundaries to ensure exact matches, not prefixed matches
          pattern <- paste0("\\b", tissue, "\\b")
          if (any(grepl(pattern, available_biosamples, ignore.case = TRUE))) {
            has_regular_tissues <- TRUE
            break
          }
        }
        
        # If no regular tissues are available, switch to ALL mode immediately
        if (!has_regular_tissues) {
          cat("Safety check: REGULAR mode selected but no regular tissues available. Switching to ALL.\n")
          updateRadioButtons(session, "biosampleMode", 
                            choices = c("All biosamples" = "ALL", "Custom subset" = "SUBSET"),
                            selected = "ALL")
          approved_biosamples(NULL)
          return()
        }
      }
    }
    
    if (input$biosampleMode == "ALL") {
      approved_biosamples(NULL)
    } else if (input$biosampleMode == "REGULAR") {
      # Define the regular subset tissues
      regular_tissues <- c("kidney", "lung", "heart", "brain", "large intestine", "liver", "stomach")
      
      # Get current available biosamples
      acc_meta <- accession_metadata()
      cat("REGULAR mode selected. Accession metadata:", if(is.null(acc_meta)) "NULL" else "available", "\n")
      
      if (!is.null(acc_meta) && !is.null(acc_meta$biosamples)) {
        available_biosamples <- acc_meta$biosamples
        cat("Available biosamples count:", length(available_biosamples), "\n")
        
        # Find matches for regular tissues (case-insensitive, partial matching)
        selected_biosamples <- c()
        for (tissue in regular_tissues) {
          matches <- available_biosamples[grepl(tissue, available_biosamples, ignore.case = TRUE)]
          cat("Tissue '", tissue, "' matches:", length(matches), "biosamples\n")
          selected_biosamples <- c(selected_biosamples, matches)
        }
        
        # Remove duplicates and update the checkbox selection
        selected_biosamples <- unique(selected_biosamples)
        cat("Total unique selected biosamples:", length(selected_biosamples), "\n")
        
        if (length(selected_biosamples) > 0) {
          updateCheckboxGroupInput(session, "biosampleCheckboxes", selected = selected_biosamples)
          
          # Show notification about auto-selection
          showNotification(
            paste("Auto-selected", length(selected_biosamples), "regular tissues:", 
                  paste(head(selected_biosamples, 3), collapse = ", "), 
                  if (length(selected_biosamples) > 3) "..." else ""),
            type = "message",
            duration = 3
          )
        } else {
          showNotification(
            "No regular tissues found in available biosamples. Please use Custom subset instead.",
            type = "warning",
            duration = 5
          )
        }
      } else {
        cat("No accession metadata or biosamples available yet\n")
        showNotification(
          "Accession metadata not loaded yet. Please wait for data to load, then try selecting Regular subset again.",
          type = "warning",
          duration = 5
        )
      }
    }
  })
  
  # Helper to draw a dendrogram from a similarity matrix (for BEDTOOLS mode)
  draw_dendrogram_from_matrix <- function(sim_matrix, labels_map, linkage, mode_label = "BEDTOOLS") {
    if (is.null(sim_matrix) || nrow(sim_matrix) < 2) {
      plot.new(); title("Not enough samples for dendrogram"); return(invisible())
    }
    
    # Convert similarity to distance
    dist_matrix <- 1 - sim_matrix
    diag(dist_matrix) <- 0
    
    # Create distance object
    d <- as.dist(dist_matrix)
    
    # Perform hierarchical clustering
    hc <- hclust(d, method = linkage)
    
    # Create dendrogram
    dend <- as.dendrogram(hc)
    
    # Color the dendrogram if we have labels
    if (!is.null(labels_map)) {
      # Get labels for samples in the matrix
      sample_names <- rownames(sim_matrix)
      sample_labels <- sapply(sample_names, function(name) {
        if (name %in% names(labels_map)) {
          labels_map[[name]]
        } else {
          "unknown"
        }
      })
      
      # Create color palette
      unique_labels <- unique(sample_labels)
      palette_fun <- if (requireNamespace('scales', quietly = TRUE)) scales::hue_pal() else function(n) grDevices::rainbow(n, s = 0.7, v = 0.9)
      palette <- palette_fun(max(3, length(unique_labels)))
      lab_to_col <- setNames(palette[seq_along(unique_labels)], unique_labels)
      
      # Colorize the dendrogram
      colorize <- function(node) {
        if (is.leaf(node)) {
          lab <- attr(node, 'label')
          col <- lab_to_col[[ sample_labels[[lab]] ]] %||% '#888888'
          ep <- attr(node, 'edgePar'); if (is.null(ep)) ep <- list(); ep$col <- col; attr(node, 'edgePar') <- ep
          np <- attr(node, 'nodePar'); if (is.null(np)) np <- list(); np$lab.col <- col; attr(node, 'nodePar') <- np
          return(node)
        } else {
          for (i in seq_along(node)) node[[i]] <- colorize(node[[i]])
          leaves <- labels(node)
          labs <- unique(sample_labels[leaves])
          col <- if (length(labs) == 1) lab_to_col[[ labs[[1]] ]] %||% '#888888' else '#888888'
          ep <- attr(node, 'edgePar'); if (is.null(ep)) ep <- list(); ep$col <- col; attr(node, 'edgePar') <- ep
          return(node)
        }
      }
      dend <- colorize(dend)
    }
    
    # Plot the dendrogram
    plot(dend, main = paste("Hierarchical Clustering -", mode_label), 
         xlab = "Samples", ylab = "Distance")
  }
  
  # Helper to compute and draw a colored dendrogram for a given file (with biological filtering support)
  draw_dendrogram_for <- function(input_file, accession_key, linkage) {
    # Check if we should use filtered similarity matrix
    src_file <- basename(input_file)
    sim_matrices <- cached_similarity_matrices()
    acc_meta <- accession_metadata()
    
    # Get effective biosample selection
    effective_biosamples <- if (input$biosampleMode == "ALL") {
      acc_meta$biosamples
    } else if (input$biosampleMode == "REGULAR") {
      # For REGULAR mode, use the selected checkboxes directly
      # Add NULL check to prevent errors when no regular tissues are available
      if (!is.null(input$biosampleCheckboxes) && length(input$biosampleCheckboxes) > 0) {
        input$biosampleCheckboxes
      } else {
        # Fallback to all biosamples if no checkboxes are selected
        acc_meta$biosamples
      }
    } else {
      # For SUBSET mode, use approved_biosamples
      approved_biosamples()
    }
    
    has_biosample_filter <- FALSE
    if ((input$biosampleMode == "SUBSET" || input$biosampleMode == "REGULAR") && 
        !is.null(effective_biosamples) && 
        length(effective_biosamples) > 0 && 
        !is.null(acc_meta) && 
        !is.null(acc_meta$biosamples)) {
      has_biosample_filter <- !setequal(effective_biosamples, acc_meta$biosamples)
    }
    has_lifestage_filter <- !is.null(input$lifestageSelect) && input$lifestageSelect != "ALL"
    
    sim_df <- NULL
    
    # Try to use filtered similarity matrix if biological filtering is active
    if ((has_biosample_filter || has_lifestage_filter) && !is.null(sim_matrices) && src_file %in% names(sim_matrices)) {
      original_matrix <- sim_matrices[[src_file]]
      if (!is.null(original_matrix)) {
        filtered_matrix <- filter_similarity_matrix(original_matrix, acc_meta, effective_biosamples, input$lifestageSelect)
        if (!is.null(filtered_matrix)) {
          sim_df <- as.data.frame(filtered_matrix)
        }
      }
    }
    
    # Fallback to reading the file directly if no filtered matrix available
    if (is.null(sim_df)) {
    fmt <- detect_file_format(input_file)
    sim_df <- if (fmt == 'hammock') parse_hammock_format(input_file) else parse_bedtools_format(input_file)
    }
    
    lab_map <- load_accession_key(accession_key)
    aligned <- align_matrix_and_labels(sim_df, lab_map, quiet = TRUE)
    sim_sub <- aligned[[1]]
    labels <- aligned[[2]]
    if (nrow(sim_sub) < 2) { plot.new(); title("Not enough samples for dendrogram"); return(invisible()) }
    distance_mat <- 1 - as.matrix(sim_sub)
    diag(distance_mat) <- 0
    eff_method <- if (tolower(linkage) == 'ward') 'complete' else linkage
    hc <- hclust(as.dist(distance_mat), method = eff_method)
    # Build colored dendrogram (base)
    dnd <- as.dendrogram(hc)
    sample_names <- rownames(distance_mat)
    lab_vec <- labels[sample_names]; lab_vec[is.na(lab_vec)] <- 'unknown'
    uniq_labs <- sort(unique(lab_vec))
    palette_fun <- if (requireNamespace('scales', quietly = TRUE)) scales::hue_pal() else function(n) grDevices::rainbow(n, s = 0.7, v = 0.9)
    palette <- palette_fun(max(3, length(uniq_labs)))
    lab_to_col <- setNames(palette[seq_along(uniq_labs)], uniq_labs)
    colorize <- function(node) {
      if (is.leaf(node)) {
        lab <- attr(node, 'label')
        col <- lab_to_col[[ lab_vec[[lab]] ]] %||% '#888888'
        ep <- attr(node, 'edgePar'); if (is.null(ep)) ep <- list(); ep$col <- col; attr(node, 'edgePar') <- ep
        np <- attr(node, 'nodePar'); if (is.null(np)) np <- list(); np$lab.col <- col; attr(node, 'nodePar') <- np
        return(node)
      } else {
        for (i in seq_along(node)) node[[i]] <- colorize(node[[i]])
        leaves <- labels(node)
        labs <- unique(lab_vec[leaves])
        col <- if (length(labs) == 1) lab_to_col[[ labs[[1]] ]] %||% '#888888' else '#888888'
        ep <- attr(node, 'edgePar'); if (is.null(ep)) ep <- list(); ep$col <- col; attr(node, 'edgePar') <- ep
        return(node)
      }
    }
    dnd <- colorize(dnd)
    
    # Add filtering info to title if active
    filter_info <- ""
    if (has_biosample_filter || has_lifestage_filter) {
      filter_parts <- c()
      if (has_biosample_filter) {
        filter_parts <- c(filter_parts, sprintf("%d biosamples", length(effective_biosamples)))
      }
      if (has_lifestage_filter) {
        filter_parts <- c(filter_parts, sprintf("lifestage: %s", input$lifestageSelect))
      }
      filter_info <- sprintf(" [Filtered: %s]", paste(filter_parts, collapse = ", "))
    }
    
    plot(dnd, main = sprintf('Hierarchical Clustering (%s)%s%s', linkage, ifelse(tolower(linkage)=='ward', " [ward uses 'complete' on precomputed distances]", ""), filter_info))
    legend('topright', legend = uniq_labs, col = lab_to_col[uniq_labs], pch = 15, cex = 0.8, bty = 'n')
  }

  # Helper to filter data based on biosample and lifestage selections using cached metadata when available
  filter_by_biological_metadata <- function(data, selected_biosamples, selected_lifestage, accession_meta) {
    if (is.null(data) || nrow(data) == 0) {
      return(data)
    }
    
    # If no filtering selections made, return all data
    if ((is.null(selected_biosamples) || length(selected_biosamples) == 0) && 
        (is.null(selected_lifestage) || selected_lifestage == "ALL")) {
      return(data)
    }
    
    valid_source_files <- character(0)
    
    for (src_file in unique(data$source_file)) {
      # First try to use cached metadata for fast filtering
      cache_path <- NULL
      cached_sample_metadata <- NULL
      
      # Try to find and load cache for this source file
      dir_path <- effective_dir()
      candidate_path <- file.path(dir_path, src_file)
      if (!file.exists(candidate_path)) {
        hits <- Sys.glob(file.path(dir_path, "**", src_file))
        if (length(hits) >= 1) candidate_path <- hits[1]
      }
      
      if (file.exists(candidate_path)) {
        cache_path <- file.path(dirname(candidate_path), paste0(basename(candidate_path), ".prep.rds"))
        
        if (file.exists(cache_path)) {
          tryCatch({
            cached <- readRDS(cache_path)
            # Support both v2 and v3 cache versions for sample metadata
            if (!is.null(cached$sample_metadata) && 
                (!is.null(cached$cache_version) && cached$cache_version %in% c("v2", "v3"))) {
              cached_sample_metadata <- cached$sample_metadata
            }
          }, error = function(e) {
            cached_sample_metadata <<- NULL
          })
        }
      }
      
      # Use cached metadata if available, otherwise fall back to reading the file
      if (!is.null(cached_sample_metadata)) {
        # Fast path: use cached sample metadata
        file_should_be_included <- FALSE
        
        for (sample_name in names(cached_sample_metadata)) {
          sample_meta <- cached_sample_metadata[[sample_name]]
          
          # Check biosample filter
          biosample_match <- TRUE
          if (!is.null(selected_biosamples) && length(selected_biosamples) > 0) {
            biosample_match <- sample_meta$biosample %in% selected_biosamples
          }
          
          # Check lifestage filter
          lifestage_match <- TRUE
          if (!is.null(selected_lifestage) && selected_lifestage != "ALL" && nzchar(selected_lifestage)) {
            lifestage_match <- sample_meta$lifestage == selected_lifestage
          }
          
          if (biosample_match && lifestage_match) {
            file_should_be_included <- TRUE
            break
          }
        }
        
        if (file_should_be_included) {
          valid_source_files <- c(valid_source_files, src_file)
        }
        
      } else {
        # Slow path: read file and check against accession key
        if (is.null(accession_meta) || is.null(accession_meta$full_data)) {
          # If no accession metadata, include all files
          valid_source_files <- c(valid_source_files, src_file)
          next
        }
        
        acc_df <- accession_meta$full_data
        if (is.null(acc_df) || !"file_base" %in% names(acc_df)) {
          valid_source_files <- c(valid_source_files, src_file)
          next
        }
        
        # Filter accession data by selections
        filtered_acc <- acc_df
        if (!is.null(selected_biosamples) && length(selected_biosamples) > 0) {
          if ("Biosample_term_name" %in% names(filtered_acc)) {
            filtered_acc <- filtered_acc[filtered_acc$Biosample_term_name %in% selected_biosamples, , drop = FALSE]
          }
        }
        
        if (!is.null(selected_lifestage) && selected_lifestage != "ALL" && nzchar(selected_lifestage)) {
          if ("Life_stage" %in% names(filtered_acc)) {
            filtered_acc <- filtered_acc[filtered_acc$Life_stage == selected_lifestage, , drop = FALSE]
          }
        }
        
        included_files <- unique(filtered_acc$file_base)
        
        # Check if this source file contains any of the included samples
        if (file.exists(candidate_path)) {
          tryCatch({
            fmt <- detect_file_format(candidate_path)
            sim_df <- if (fmt == 'hammock') parse_hammock_format(candidate_path) else parse_bedtools_format(candidate_path)
            sample_names <- rownames(sim_df)
            
            overlap <- intersect(sample_names, included_files)
            if (length(overlap) > 0) {
              valid_source_files <- c(valid_source_files, src_file)
            }
          }, error = function(e) {
            # If we can't read the file, include it by default
            valid_source_files <<- c(valid_source_files, src_file)
          })
        } else {
          # If we can't find the file, include it by default
          valid_source_files <- c(valid_source_files, src_file)
        }
      }
    }
    
    # Filter the data to only include valid source files
    if (length(valid_source_files) > 0) {
      data <- data[data$source_file %in% valid_source_files, , drop = FALSE]
    } else {
      # If no files match, return empty data frame with same structure
      data <- data[FALSE, , drop = FALSE]
    }
    
    return(data)
  }

  # Helper to read accession key and return labels for tissue and organism
  get_accession_labels_df <- function(key_path, sample_names) {
    if (is.null(key_path) || !nzchar(key_path) || !file.exists(key_path)) {
      return(data.frame(file_base = sample_names, tissue = rep("unknown", length(sample_names)), organism = rep("unknown", length(sample_names)), stringsAsFactors = FALSE))
    }
    df <- tryCatch({
      if (requireNamespace("data.table", quietly = TRUE)) data.table::fread(key_path, sep = "\t", header = TRUE, showProgress = FALSE)
      else utils::read.table(key_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    }, error = function(e) NULL)
    if (is.null(df) || !nrow(df)) return(data.frame(file_base = sample_names, tissue = rep("unknown", length(sample_names)), organism = rep("unknown", length(sample_names)), stringsAsFactors = FALSE))
    file_col <- if ("File" %in% names(df)) "File" else if ("file" %in% names(df)) "file" else NULL
    tissue_col <- if ("Biosample_term_name" %in% names(df)) "Biosample_term_name" else if ("tissue" %in% names(df)) "tissue" else if ("Tissue" %in% names(df)) "Tissue" else NULL
    org_col <- if ("Organism" %in% names(df)) "Organism" else if ("organism" %in% names(df)) "organism" else NULL
    if (is.null(file_col)) return(data.frame(file_base = sample_names, tissue = rep("unknown", length(sample_names)), organism = rep("unknown", length(sample_names)), stringsAsFactors = FALSE))
    file_base <- vapply(df[[file_col]], function(x) tools::file_path_sans_ext(basename(as.character(x))), character(1))
    tissue <- if (!is.null(tissue_col)) as.character(df[[tissue_col]]) else rep("unknown", length(file_base))
    organism <- if (!is.null(org_col)) as.character(df[[org_col]]) else rep("unknown", length(file_base))
    lab_df <- data.frame(file_base = file_base, tissue = tissue, organism = organism, stringsAsFactors = FALSE)
    # Keep only relevant sample_names and fill missing with 'unknown'
    out <- data.frame(file_base = sample_names, stringsAsFactors = FALSE)
    out <- merge(out, lab_df, by = "file_base", all.x = TRUE, sort = FALSE)
    out$tissue[is.na(out$tissue) | !nzchar(out$tissue)] <- "unknown"
    out$organism[is.na(out$organism) | !nzchar(out$organism)] <- "unknown"
    out
  }

  # Helper to draw a similarity heatmap for a given file with side color annotations (with biological filtering support)
  draw_heatmap_for <- function(input_file, accession_key) {
    # Check if we should use filtered similarity matrix
    src_file <- basename(input_file)
    sim_matrices <- cached_similarity_matrices()
    acc_meta <- accession_metadata()
    
    # Get effective biosample selection
    effective_biosamples <- if (input$biosampleMode == "ALL") {
      acc_meta$biosamples
    } else if (input$biosampleMode == "REGULAR") {
      # For REGULAR mode, use the selected checkboxes directly
      # Add NULL check to prevent errors when no regular tissues are available
      if (!is.null(input$biosampleCheckboxes) && length(input$biosampleCheckboxes) > 0) {
        input$biosampleCheckboxes
      } else {
        # Fallback to all biosamples if no checkboxes are selected
        acc_meta$biosamples
      }
    } else {
      # For SUBSET mode, use approved_biosamples
      approved_biosamples()
    }
    
    has_biosample_filter <- FALSE
    if ((input$biosampleMode == "SUBSET" || input$biosampleMode == "REGULAR") && 
        !is.null(effective_biosamples) && 
        length(effective_biosamples) > 0 && 
        !is.null(acc_meta) && 
        !is.null(acc_meta$biosamples)) {
      has_biosample_filter <- !setequal(effective_biosamples, acc_meta$biosamples)
    }
    has_lifestage_filter <- !is.null(input$lifestageSelect) && input$lifestageSelect != "ALL"
    
    mat <- NULL
    
    # Try to use filtered similarity matrix if biological filtering is active
    if ((has_biosample_filter || has_lifestage_filter) && !is.null(sim_matrices) && src_file %in% names(sim_matrices)) {
      original_matrix <- sim_matrices[[src_file]]
      if (!is.null(original_matrix)) {
        filtered_matrix <- filter_similarity_matrix(original_matrix, acc_meta, effective_biosamples, input$lifestageSelect)
        if (!is.null(filtered_matrix)) {
          mat <- filtered_matrix
        }
      }
    }
    
    # Fallback to reading the file directly if no filtered matrix available
    if (is.null(mat)) {
    fmt <- detect_file_format(input_file)
    sim_df <- if (fmt == 'hammock') parse_hammock_format(input_file) else parse_bedtools_format(input_file)
    mat <- as.matrix(sim_df)
    }
    
    if (nrow(mat) < 2 || ncol(mat) < 2) { plot.new(); title("Not enough samples for heatmap"); return(invisible()) }
    # Labels for side colors
    samples <- rownames(mat)
    labs_df <- get_accession_labels_df(accession_key, samples)
    # Convert to long form
    df_long <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    names(df_long) <- c("row", "col", "val")
    df_long$row <- as.character(df_long$row)
    df_long$col <- as.character(df_long$col)
    df_long$val <- as.numeric(df_long$val)
    # Extend x-axis with side columns
    cols <- colnames(mat)
    col_levels <- c("_tissue", cols, "_organism")
    df_long$col <- factor(df_long$col, levels = col_levels)
    df_long$row <- factor(df_long$row, levels = rownames(mat))
    # Side annotation points
    df_tissue <- data.frame(col = factor(rep("_tissue", length(samples)), levels = col_levels), row = factor(samples, levels = rownames(mat)), label = labs_df$tissue, stringsAsFactors = FALSE)
    df_org <- data.frame(col = factor(rep("_organism", length(samples)), levels = col_levels), row = factor(samples, levels = rownames(mat)), label = labs_df$organism, stringsAsFactors = FALSE)
    # Palettes for side bars and label coloring
    uniq_labels <- sort(unique(c(as.character(df_tissue$label), as.character(df_org$label))))
    palette_fun <- if (requireNamespace('scales', quietly = TRUE)) scales::hue_pal() else function(n) grDevices::rainbow(n, s = 0.7, v = 0.9)
    col_map <- setNames(palette_fun(max(3, length(uniq_labels))), uniq_labels)
    # Colored axis labels: x by tissue, y by organism (using ggtext if available)
    samples_vec <- as.character(samples)
    tissue_for <- setNames(as.character(labs_df$tissue), labs_df$file_base)
    org_for <- setNames(as.character(labs_df$organism), labs_df$file_base)
    color_for_tissue <- function(s) col_map[[ tissue_for[[s]] %||% "unknown" ]] %||% "#000000"
    color_for_org <- function(s) col_map[[ org_for[[s]] %||% "unknown" ]] %||% "#000000"
    use_md <- requireNamespace('ggtext', quietly = TRUE)
    color_span <- function(label, color) if (use_md) paste0("<span style='color:", color, "'>", label, "</span>") else label
    x_labels_named <- setNames(c("", vapply(samples_vec, function(s) color_span(s, color_for_tissue(s)), character(1)), ""), col_levels)
    y_labels_named <- setNames(vapply(samples_vec, function(s) color_span(s, color_for_org(s)), character(1)), rownames(mat))
    # Plot
    p <- ggplot(df_long, aes(x = col, y = row, fill = val)) +
      geom_raster() +
      scale_fill_viridis_c(option = "viridis", limits = c(0, 1), oob = scales::squish, na.value = "grey95") +
      geom_point(data = df_tissue, aes(x = col, y = row, color = label), inherit.aes = FALSE, shape = 15, size = 2.8) +
      geom_point(data = df_org, aes(x = col, y = row, color = label), inherit.aes = FALSE, shape = 15, size = 2.8) +
      scale_color_manual(values = col_map, drop = FALSE) +
      scale_x_discrete(labels = x_labels_named) +
      scale_y_discrete(labels = y_labels_named) +
      coord_fixed() +
      labs(x = NULL, y = NULL, fill = "Similarity", color = "Labels", title = paste0(basename(input_file), 
           if (has_biosample_filter || has_lifestage_filter) {
             filter_parts <- c()
             if (has_biosample_filter) filter_parts <- c(filter_parts, sprintf("%d biosamples", length(effective_biosamples)))
             if (has_lifestage_filter) filter_parts <- c(filter_parts, sprintf("lifestage: %s", input$lifestageSelect))
             sprintf(" [Filtered: %s]", paste(filter_parts, collapse = ", "))
           } else "")) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = if (use_md) ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1, size = 7) else element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        axis.text.y = if (use_md) ggtext::element_markdown(size = 7) else element_text(size = 7),
        panel.grid = element_blank()
      )
    print(p)
  }
  # Build more precise roots: project root and experiments dir if available
  app_wd <- getwd()
  project_root <- tryCatch(normalizePath(file.path(app_wd, "..", "..", "..")), error = function(e) getwd())
  experiments_dir <- file.path(project_root, "experiments")
  vols <- c(
    if (dir.exists(project_root)) c(Project = project_root) else NULL,
    if (dir.exists(experiments_dir)) c(Experiments = experiments_dir) else NULL,
    c(Home = path.expand("~"), Root = "/", Tmp = tempdir())
  )
  volumes <- unlist(vols)
  shinyDirChoose(input, "dir", roots = volumes, session = session)
  shinyFileChoose(input, "acc", roots = volumes, session = session,
                  defaultRoot = if ("Project" %in% names(volumes)) "Project" else names(volumes)[1],
                  defaultPath = ".")
  # No separate cache file picker when using same-folder caches
  dir_path <- reactive({
    sel <- input$dir
    if (is.null(sel) || length(sel) == 0) return(NA_character_)
    tryCatch({
      p <- as.character(shinyFiles::parseDirPath(volumes, sel))
      if (length(p) == 1 && nzchar(p)) p else NA_character_
    }, error = function(e) NA_character_)
  })
  acc_path <- reactive({
    sel <- input$acc
    if (is.null(sel) || length(sel) == 0) return(NA_character_)
    tryCatch({
      p <- as.character(shinyFiles::parseFilePaths(volumes, sel)$datapath)
      if (length(p) >= 1) p[[1]] else NA_character_
    }, error = function(e) NA_character_)
  })
  # Sync file browser selections into text boxes
  observeEvent(dir_path(), {
    p <- dir_path()
    if (!is.null(p) && !is.na(p) && nzchar(p)) {
      updateTextInput(session, "dirText", value = p)
    }
  }, ignoreInit = TRUE)
  observeEvent(acc_path(), {
    p <- acc_path()
    if (!is.null(p) && !is.na(p) && nzchar(p)) {
      updateTextInput(session, "accText", value = p)
    }
  }, ignoreInit = TRUE)

  # Text overrides: if user pastes a path and it exists, prefer it
  effective_dir <- reactive({
    txt <- trimws(input$dirText %||% "")
    if (nzchar(txt) && dir.exists(txt)) return(normalizePath(txt))
    dir_path()
  })
  effective_acc <- reactive({
    txt <- trimws(input$accText %||% "")
    if (nzchar(txt) && file.exists(txt)) return(normalizePath(txt))
    acc_path()
  })

  silhouette_limits <- reactive({
    df <- data_all()
    if (is.null(df) || !("silhouette" %in% names(df))) return(c(0, 1))
    vals <- df$silhouette
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) return(c(0, 1))
    range(vals, na.rm = TRUE)
  })

  observeEvent(list(effective_dir(), effective_acc()), {
    path <- effective_dir()
    ak <- effective_acc()
    if (is.na(ak) || !file.exists(ak)) {
      data_all(NULL)
      output$statusText <- renderText("Please select an accession key TSV.")
      return()
    }
    if (is.na(path) || !dir.exists(path)) {
      data_all(NULL)
      output$statusText <- renderText("Select a directory of hammock or bedtools output files.")
      return()
    }
    
    # Show initial loading message
    output$statusText <- renderText("Scanning directory for valid files...")
    
    # Pre-scan for hammock or bedtools outputs
    all_files <- list.files(path, pattern = "\\.(csv|tsv|txt)$", ignore.case = TRUE, full.names = TRUE)
    is_valid_format <- function(fp) {
      # Use the detect_file_format function from utils.R
      tryCatch({
        fmt <- detect_file_format(fp)
        fmt %in% c('hammock', 'bedtools')
      }, error = function(e) FALSE)
    }
    valid_files <- Filter(is_valid_format, all_files)
    if (length(valid_files) == 0) {
      data_all(NULL)
      output$statusText <- renderText("No hammock or bedtools outputs detected in the selected directory.")
      output$diagText <- renderText(sprintf("Scanned %d files in %s", length(all_files), path))
      return()
    }
    
    # Extract accession metadata early for immediate UI updates
    output$statusText <- renderText("Loading accession metadata...")
    acc_meta <- extract_accession_metadata(ak)
    accession_metadata(acc_meta)
    
    # Update UI immediately with basic information
    output$statusText <- renderText(sprintf("Found %d valid files. Loading data...", length(valid_files)))
    
    # Update biosample and lifestage choices immediately
    updateCheckboxGroupInput(session, "biosampleCheckboxes", 
                           choices = setNames(acc_meta$biosamples, acc_meta$biosamples),
                           selected = character(0))
    
    lifestage_choices <- c("ALL" = "ALL")
    if (length(acc_meta$lifestages) > 0) {
      lifestage_choices <- c(lifestage_choices, setNames(acc_meta$lifestages, acc_meta$lifestages))
    }
    updateRadioButtons(session, "lifestageSelect", choices = lifestage_choices, selected = "ALL")
    
    # Start loading data in background
    output$statusText <- renderText(sprintf("Running clustering on %d files... (key: %s)", length(valid_files), ak))
    res <- NULL
    err_top <- NULL
    
    # Simpler direct call (avoid nested progress scoping issues)
    tryCatch({
      res <- run_clustering_on_dir(
        dir_path = path,
        acc_key = ak,
        files = valid_files,
        use_cache = isTRUE(input$useCache)
      )
    }, error = function(e) { 
      err_top <<- conditionMessage(e); 
      res <<- NULL 
    })
    
    if (is.null(res)) {
      data_all(NULL)
      output$statusText <- renderText("Clustering failed to return a result. See diagnostics.")
      
      # Enhanced error diagnostics
      error_debug_info <- c(
        paste("=== ERROR DIAGNOSTICS ==="),
        paste("Error message:", if (!is.null(err_top)) err_top else "(no result from run)"),
        paste("=== BIOSAMPLE MODE DEBUG ==="),
        paste("Current biosample mode:", input$biosampleMode),
        paste("Biosample checkboxes:", if(is.null(input$biosampleCheckboxes)) "NULL" else paste(input$biosampleCheckboxes, collapse = ", ")),
        paste("Approved biosamples:", if(is.null(approved_biosamples())) "NULL" else paste(approved_biosamples(), collapse = ", ")),
        paste("=== ACCESSION METADATA DEBUG ==="),
        paste("Accession metadata loaded:", if(is.null(accession_metadata())) "NO" else "YES"),
        if(!is.null(accession_metadata())) {
          acc_meta <- accession_metadata()
          c(
            paste("Available biosamples:", paste(acc_meta$biosamples, collapse = ", ")),
            paste("Available lifestages:", paste(acc_meta$lifestages, collapse = ", ")),
            paste("=== REGULAR TISSUE CHECK ==="),
            paste("Regular tissues to check: kidney, lung, heart, brain, large intestine, liver, stomach"),
            paste("Has regular tissues:", any(sapply(c("kidney", "lung", "heart", "brain", "large intestine", "liver", "stomach"), 
                                                    function(t) any(grepl(paste0("\\b", t, "\\b"), acc_meta$biosamples, ignore.case = TRUE)))))
          )
        } else {
          "No accession metadata available"
        }
      )
      
      output$diagText <- renderText(paste(error_debug_info, collapse = "\n"))
      return()
    }
    if (is.null(res$df)) {
      data_all(NULL)
      output$statusText <- renderText("No valid hammock outputs were processed. Check files and key.")
      output$diagText <- renderText(paste(res$diag, collapse = "\n"))
      return()
    }
    
    # Set data immediately
    data_all(res$df)
    
    # Enhanced diagnostics with debugging information
    debug_info <- c(
      paste("=== LOADING DIAGNOSTICS ==="),
      paste("Data loaded successfully:", length(unique(res$df$source_file)), "files,", nrow(res$df), "rows"),
      paste("Accession key:", ak),
      paste("Cache enabled:", isTRUE(input$useCache)),
      paste("=== BIOSAMPLE MODE DEBUG ==="),
      paste("Current biosample mode:", input$biosampleMode),
      paste("Biosample checkboxes:", if(is.null(input$biosampleCheckboxes)) "NULL" else paste(input$biosampleCheckboxes, collapse = ", ")),
      paste("Approved biosamples:", if(is.null(approved_biosamples())) "NULL" else paste(approved_biosamples(), collapse = ", ")),
      paste("=== ACCESSION METADATA DEBUG ==="),
      paste("Accession metadata loaded:", if(is.null(acc_meta)) "NO" else "YES"),
      if(!is.null(acc_meta)) {
        c(
          paste("Available biosamples:", paste(acc_meta$biosamples, collapse = ", ")),
          paste("Available lifestages:", paste(acc_meta$lifestages, collapse = ", ")),
          paste("=== REGULAR TISSUE CHECK ==="),
          paste("Regular tissues to check: kidney, lung, heart, brain, large intestine, liver, stomach"),
                      paste("Has regular tissues:", any(sapply(c("kidney", "lung", "heart", "brain", "large intestine", "liver", "stomach"), 
                                                    function(t) any(grepl(paste0("\\b", t, "\\b"), acc_meta$biosamples, ignore.case = TRUE))))),
          paste("=== EFFECTIVE BIOSAMPLES DEBUG ==="),
          paste("Effective biosamples (ALL mode):", paste(acc_meta$biosamples, collapse = ", ")),
          if(input$biosampleMode == "REGULAR") {
            paste("Effective biosamples (REGULAR mode):", 
                  if(is.null(input$biosampleCheckboxes) || length(input$biosampleCheckboxes) == 0) "NULL/EMPTY" 
                  else paste(input$biosampleCheckboxes, collapse = ", "))
          } else NULL,
          if(input$biosampleMode == "SUBSET") {
            paste("Effective biosamples (SUBSET mode):", 
                  if(is.null(approved_biosamples()) || length(approved_biosamples()) == 0) "NULL/EMPTY" 
                  else paste(approved_biosamples(), collapse = ", "))
          } else NULL
        )
      } else {
        "No accession metadata available"
      },
      paste("=== ORIGINAL DIAGNOSTICS ==="),
      res$diag
    )
    
    output$diagText <- renderText(paste(debug_info, collapse = "\n"))
    
    # Extract and store accession metadata for filtering
    cat(sprintf("Extracting accession metadata from: %s\n", ak))
    acc_meta <- extract_accession_metadata(ak)
    cat(sprintf("Found %d biosamples and %d lifestages\n", length(acc_meta$biosamples), length(acc_meta$lifestages)))
    accession_metadata(acc_meta)
    
    # Cache similarity matrices for each source file for real-time recalculation
    # Note: When cache is enabled, run_clustering_on_dir already loads the data efficiently
    # We still need to load similarity matrices for filtering, but we can optimize this
    output$statusText <- renderText("Loading similarity matrices for filtering...")
    similarity_cache <- list()
    source_files <- unique(res$df$source_file)
    cat(sprintf("Loading similarity matrices for %d source files\n", length(source_files)))
    
    # Load similarity matrices efficiently
    for (i in seq_along(source_files)) {
      src_file <- source_files[i]
      
      if (i %% 10 == 0) {  # Update progress every 10 files for better performance
        output$statusText <- renderText(sprintf("Loading similarity matrices... (%d/%d)", i, length(source_files)))
      }
      
      candidate_path <- file.path(path, src_file)
      if (!file.exists(candidate_path)) {
        hits <- Sys.glob(file.path(path, "**", src_file))
        if (length(hits) >= 1) candidate_path <- hits[1]
      }
      
      if (file.exists(candidate_path)) {
        # Try to load from cache first (if cache is enabled)
        cache_path <- file.path(dirname(candidate_path), paste0(basename(candidate_path), ".prep.rds"))
        sim_matrix <- NULL
        
        if (isTRUE(input$useCache) && file.exists(cache_path)) {
          tryCatch({
            cached <- readRDS(cache_path)
            if (!is.null(cached$similarity)) {
              sim_matrix <- cached$similarity
            }
          }, error = function(e) {
            cat(sprintf("Warning: Could not load cache file %s: %s\n", cache_path, conditionMessage(e)))
            sim_matrix <<- NULL
          })
        }
        
        # If not in cache, read directly
        if (is.null(sim_matrix)) {
          tryCatch({
            fmt <- detect_file_format(candidate_path)
            sim_df <- if (fmt == 'hammock') parse_hammock_format(candidate_path) else parse_bedtools_format(candidate_path)
            sim_matrix <- as.matrix(sim_df)
          }, error = function(e) {
            sim_matrix <<- NULL
          })
        }
        
        if (!is.null(sim_matrix)) {
          similarity_cache[[src_file]] <- sim_matrix
        }
      }
    }
    
    cached_similarity_matrices(similarity_cache)
    
    # Clear clustering results cache when new data is loaded
    clustering_results_cache(list())
    
    # Final UI updates
    output$statusText <- renderText("Initializing UI controls...")
    
    # Update biosample and lifestage choices (already done earlier, but ensure they're current)
    cat(sprintf("Updating UI with %d biosamples: %s\n", length(acc_meta$biosamples), paste(acc_meta$biosamples, collapse = ", ")))
    cat(sprintf("Updating UI with %d lifestages: %s\n", length(acc_meta$lifestages), paste(acc_meta$lifestages, collapse = ", ")))
    
    # Initialize controls
    modes <- sort(unique(na.omit(res$df$mode)))
    if (length(modes)) updateRadioButtons(session, "modeSelect", selected = modes[1])
    linkages <- sort(unique(na.omit(res$df$linkage)))
    if (!length(linkages)) linkages <- c("average","complete","single","ward")
    updateSelectInput(session, "linkageMethod", choices = linkages, selected = if (length(linkages)) linkages[1] else NULL)
    # Initialize precision based on first mode
    if (length(modes) > 0) {
      mode_precisions <- sort(unique(na.omit(res$df$precision[res$df$mode == modes[1]])))
      updateSelectInput(session, "precision", choices = mode_precisions, selected = if (length(mode_precisions)) mode_precisions[1] else NULL)
    }
    # Collect expA across files by re-parsing filenames if column missing or NA
    exp_vals <- suppressWarnings(as.numeric(na.omit(res$df$expA)))
    if (length(exp_vals) == 0 || all(is.na(exp_vals))) {
      # Fall back: parse from source_file names
      exp_from_name <- lapply(unique(res$df$source_file), function(sf) {
        p <- tryCatch(extract_params_from_filename(sf), error = function(e) NULL)
        if (!is.null(p) && is.list(p)) p$expA else NA_real_
      })
      exp_vals <- suppressWarnings(as.numeric(na.omit(unlist(exp_from_name))))
    }
    exp_vals <- sort(unique(exp_vals))
    exp_labels <- if (length(exp_vals)) sprintf("%.2f", exp_vals) else character(0)
    exp_choices <- if (length(exp_vals)) stats::setNames(as.character(exp_vals), exp_labels) else character(0)
    updateSelectInput(session, "expA", choices = exp_choices, selected = if (length(exp_vals)) as.character(exp_vals[1]) else "0")
    
    # Collect subB across files by re-parsing filenames if column missing or NA
    subB_vals <- suppressWarnings(as.numeric(na.omit(res$df$subB)))
    if (length(subB_vals) == 0 || all(is.na(subB_vals))) {
      # Fall back: parse from source_file names
      subB_from_name <- lapply(unique(res$df$source_file), function(sf) {
        p <- tryCatch(extract_params_from_filename(sf), error = function(e) NULL)
        if (!is.null(p) && is.list(p)) p$subB else NA_real_
      })
      subB_vals <- suppressWarnings(as.numeric(na.omit(unlist(subB_from_name))))
    }
    subB_vals <- sort(unique(subB_vals))
    subB_labels <- if (length(subB_vals)) sprintf("%.2f", subB_vals) else character(0)
    subB_choices <- if (length(subB_vals)) stats::setNames(as.character(subB_vals), subB_labels) else character(0)
    updateSelectInput(session, "subB", choices = subB_choices, selected = if (length(subB_vals)) as.character(subB_vals[1]) else "1")
    k_choices <- sort(unique(na.omit(res$df$klen)))
    w_choices <- sort(unique(na.omit(res$df$window)))
    updateSelectInput(session, "klen", choices = k_choices, selected = if (length(k_choices)) k_choices[1] else NULL)
    updateSelectInput(session, "window", choices = w_choices, selected = if (length(w_choices)) w_choices[1] else NULL)
    srcs <- sort(unique(res$df$source_file))
    cat("DEBUG: All source files:", paste(srcs, collapse = ", "), "\n")
    updateSelectInput(session, "sourceFile", choices = srcs, selected = if (length(srcs)) srcs[1] else NULL)
    updateSelectInput(session, "heatmapFile1", choices = srcs, selected = if (length(srcs)) srcs[1] else NULL)
    updateSelectInput(session, "heatmapFile2", choices = srcs, selected = if (length(srcs) >= 2) srcs[2] else if (length(srcs)) srcs[1] else NULL)
    updateSelectInput(session, "dendFile1", choices = srcs, selected = if (length(srcs)) srcs[1] else NULL)
    updateSelectInput(session, "dendFile2", choices = srcs, selected = if (length(srcs) >= 2) srcs[2] else if (length(srcs)) srcs[1] else NULL)
    
    # Update PCA file selectors
    updateSelectInput(session, "pcaFile1", choices = srcs, selected = if (length(srcs)) srcs[1] else NULL)
    updateSelectInput(session, "pcaFile2", choices = srcs, selected = if (length(srcs) >= 2) srcs[2] else if (length(srcs)) srcs[1] else NULL)
    
    # Update bedtools file selector and auto-select for BEDTOOLS mode
    bedtools_files <- srcs[grepl("bedtools", srcs, ignore.case = TRUE)]
    cat("DEBUG: Found bedtools files:", paste(bedtools_files, collapse = ", "), "\n")
    updateSelectInput(session, "bedtoolsFile", choices = bedtools_files, 
                     selected = if (length(bedtools_files)) bedtools_files[1] else NULL)
    
    # Auto-select bedtools file in sourceFile when BEDTOOLS mode is selected
    if (!is.null(input$modeSelect) && input$modeSelect == "BEDTOOLS" && length(bedtools_files) > 0) {
      updateSelectInput(session, "sourceFile", selected = bedtools_files[1])
    }
    
    # Keep diagnostics minimal (overlap/expA from run_clustering_on_dir)
    output$diagText <- renderText(paste(res$diag, collapse = "\n"))
    output$statusText <- renderText(sprintf("✓ Data loaded successfully! %d files (%d rows) using key: %s", length(unique(res$df$source_file)), nrow(res$df), ak))
  })

  # Reactive that recomputes clustering on filtered similarity matrices with caching
  filtered_clustering_data <- reactive({
    original_data <- data_all()
    if (is.null(original_data)) return(NULL)
    
    acc_meta <- accession_metadata()
    sim_matrices <- cached_similarity_matrices()
    
    # Get effective biosample selection based on mode
    effective_biosamples <- if (input$biosampleMode == "ALL") {
      acc_meta$biosamples
    } else if (input$biosampleMode == "REGULAR") {
      # For REGULAR mode, use the selected checkboxes directly
      # Add NULL check to prevent errors when no regular tissues are available
      if (!is.null(input$biosampleCheckboxes) && length(input$biosampleCheckboxes) > 0) {
        input$biosampleCheckboxes
      } else {
        # Fallback to all biosamples if no checkboxes are selected
        acc_meta$biosamples
      }
    } else {
      # For SUBSET mode, use approved_biosamples
      approved_biosamples()
    }
    
    # Check if we have biological filtering active
    has_biosample_filter <- FALSE
    if ((input$biosampleMode == "SUBSET" || input$biosampleMode == "REGULAR") && 
        !is.null(effective_biosamples) && 
        length(effective_biosamples) > 0 && 
        !is.null(acc_meta) && 
        !is.null(acc_meta$biosamples)) {
      has_biosample_filter <- !setequal(effective_biosamples, acc_meta$biosamples)
    }
    has_lifestage_filter <- !is.null(input$lifestageSelect) && input$lifestageSelect != "ALL"
    
    # Special case: Check if filtering would produce the same matrix (e.g., all samples are fetal)
    if (has_lifestage_filter && !is.null(acc_meta$full_data) && "Life_stage" %in% names(acc_meta$full_data)) {
      acc_df <- acc_meta$full_data
      selected_lifestage <- input$lifestageSelect
      
      # Check if all samples in all matrices are already of the selected lifestage
      all_matrices_same_lifestage <- TRUE
      for (src_file in names(sim_matrices)) {
        sim_matrix <- sim_matrices[[src_file]]
        if (!is.null(sim_matrix)) {
          sample_names <- rownames(sim_matrix)
          lifestage_samples <- acc_df[acc_df$Life_stage == selected_lifestage, "file_base", drop = FALSE]
          all_samples_same_lifestage <- all(sample_names %in% lifestage_samples$file_base)
          if (!all_samples_same_lifestage) {
            all_matrices_same_lifestage <- FALSE
            break
          }
        }
      }
      
      if (all_matrices_same_lifestage) {
        cat(sprintf("DEBUG: All samples are already %s, returning original data\n", selected_lifestage))
        return(original_data)
      }
    }
    
      # Special case: Check if filtering would produce the same matrix (e.g., all samples are fetal)
  if (has_lifestage_filter && !is.null(acc_meta$full_data) && "Life_stage" %in% names(acc_meta$full_data)) {
    acc_df <- acc_meta$full_data
    selected_lifestage <- input$lifestageSelect
    
    # Check if all samples in all matrices are already of the selected lifestage
    all_matrices_same_lifestage <- TRUE
    for (src_file in names(sim_matrices)) {
      sim_matrix <- sim_matrices[[src_file]]
      if (!is.null(sim_matrix)) {
        sample_names <- rownames(sim_matrix)
        lifestage_samples <- acc_df[acc_df$Life_stage == selected_lifestage, "file_base", drop = FALSE]
        all_samples_same_lifestage <- all(sample_names %in% lifestage_samples$file_base)
        if (!all_samples_same_lifestage) {
          all_matrices_same_lifestage <- FALSE
          break
        }
      }
    }
    
    if (all_matrices_same_lifestage) {
      cat(sprintf("DEBUG: All samples are already %s, returning original data\n", selected_lifestage))
      return(original_data)
    }
  }
  
  # If no biological filtering is active, return original data
  if (!has_biosample_filter && !has_lifestage_filter) {
    return(original_data)
  }
    
    # Create cache key for this filter combination
    cache_key <- create_filter_cache_key(effective_biosamples, input$lifestageSelect)
    
    # Check session cache first (fastest)
    results_cache <- clustering_results_cache()
    if (!is.null(results_cache) && cache_key %in% names(results_cache)) {
      cat(sprintf("Using session-cached clustering results for filter: %s\n", cache_key))
      return(results_cache[[cache_key]])
    }
    
    # Check RDS files for pre-computed filtered results (persistent cache)
    # Cache effective_dir() call to avoid repeated computation
    dir_path <- effective_dir()
    if (!is.na(dir_path) && dir.exists(dir_path)) {
      all_rds_results <- list()
      has_rds_cache <- FALSE
      
      for (src_file in names(sim_matrices)) {
        rds_result <- load_filtered_results_from_rds(src_file, cache_key, dir_path)
        if (!is.null(rds_result)) {
          all_rds_results[[src_file]] <- rds_result
          has_rds_cache <- TRUE
        }
      }
      
      if (has_rds_cache && length(all_rds_results) == length(sim_matrices)) {
        # We have complete RDS cache for all files
        cat(sprintf("Using RDS-cached clustering results for filter: %s\n", cache_key))
        
        # Combine RDS results
        combined_data <- do.call(rbind, all_rds_results)
        
        # Store in session cache for even faster future access
        updated_cache <- results_cache
        if (is.null(updated_cache)) updated_cache <- list()
        updated_cache[[cache_key]] <- combined_data
        clustering_results_cache(updated_cache)
        
        return(combined_data)
      }
    }
    
    cat(sprintf("Computing clustering results for filter: %s\n", cache_key))
    
    # If biological filtering is active, recompute clustering
    if (is.null(acc_meta) || is.null(sim_matrices)) {
      return(original_data)  # Fallback to original data
    }
    
    # Load labels map for clustering
    # Cache effective_acc() call to avoid repeated computation
    ak <- effective_acc()
    if (is.na(ak) || !file.exists(ak)) {
      return(original_data)
    }
    
    labels_map <- tryCatch({
      load_accession_key(ak)
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(labels_map)) {
      return(original_data)
    }
    
    # Recompute clustering for each source file
    recomputed_results <- list()
    
    for (src_file in names(sim_matrices)) {
      sim_matrix <- sim_matrices[[src_file]]
      if (is.null(sim_matrix)) next
      
      # Filter the similarity matrix
      filtered_matrix <- filter_similarity_matrix(sim_matrix, acc_meta, effective_biosamples, input$lifestageSelect)
      if (is.null(filtered_matrix)) next
      
      # Recompute clustering
      new_clustering <- recompute_clustering_on_filtered_matrix(filtered_matrix, labels_map)
      if (is.null(new_clustering)) next
      
      # Add metadata from original data
      original_rows <- original_data[original_data$source_file == src_file, , drop = FALSE]
      if (nrow(original_rows) > 0) {
        # Get the first row to extract metadata
        meta_row <- original_rows[1, , drop = FALSE]
        
        # Add metadata columns to new clustering results
        new_clustering$sketch <- meta_row$sketch
        new_clustering$mode <- meta_row$mode
        new_clustering$precision <- meta_row$precision
        new_clustering$expA <- meta_row$expA
        new_clustering$subA <- meta_row$subA
        new_clustering$subB <- meta_row$subB
        new_clustering$window <- meta_row$window
        new_clustering$klen <- meta_row$klen
        new_clustering$source_file <- src_file
        
        recomputed_results[[length(recomputed_results) + 1]] <- new_clustering
      }
    }
    
    if (length(recomputed_results) == 0) {
      combined_data <- data.frame()  # Return empty data frame if no valid results
    } else {
      # Combine all recomputed results
      combined_data <- do.call(rbind, recomputed_results)
      
      # Ensure all required columns exist
      required_cols <- c('sketch','mode','precision','expA','subA','subB','window','klen',
                        'n_clusters','linkage','nmi','silhouette','source_file')
      for (col in setdiff(required_cols, names(combined_data))) {
        combined_data[[col]] <- NA
      }
      
      # Reorder columns
      combined_data <- combined_data[, required_cols, drop = FALSE]
    }
    
    # Cache the results for this filter combination (session cache)
    updated_cache <- results_cache
    if (is.null(updated_cache)) updated_cache <- list()
    updated_cache[[cache_key]] <- combined_data
    clustering_results_cache(updated_cache)
    
    # Also save individual file results to RDS files (persistent cache)
    dir_path <- effective_dir()
    if (!is.na(dir_path) && dir.exists(dir_path)) {
      for (i in seq_along(recomputed_results)) {
        result <- recomputed_results[[i]]
        src_file <- result$source_file[1]  # Get source file name
        
        # Save this file's filtered results to its RDS file
        save_success <- save_filtered_results_to_rds(src_file, cache_key, result, dir_path)
        if (save_success) {
          cat(sprintf("Saved filtered results to RDS for %s, filter: %s\n", src_file, cache_key))
        }
      }
    }
    
    return(combined_data)
  })

  # Optimized: Base filtered data with biological filtering only
  base_filtered_data <- reactive({
    # Use recomputed clustering data when biological filters are active
    df <- filtered_clustering_data()
    if (is.null(df)) return(NULL)
    req_cols <- c("mode", "linkage", "precision", "n_clusters", "nmi", "silhouette")
    if (!all(req_cols %in% names(df))) {
      return(NULL)
    }
    d <- df
    
    # Only apply biological filtering, no technical filtering
    # This ensures the NMI table shows all technical combinations for the selected biological subset
    if (nrow(d) == 0) {
      return(NULL)
    }
    d
  })
  
  # Data with only biological filtering (for Inputs tab NMI table)
  biologically_filtered_data <- reactive({
    base_filtered_data()
  })
  
  filtered_data <- reactive({
    # Handle BEDTOOLS mode specially
    if (!is.null(input$modeSelect) && input$modeSelect == "BEDTOOLS") {
      cat("DEBUG: BEDTOOLS mode detected\n")
      # For BEDTOOLS mode, we need to create clustering data from the selected bedtools file
      bedtools_file <- input$bedtoolsFile
      if (is.null(bedtools_file) || !nzchar(bedtools_file)) {
        # Fallback: try to use sourceFile if it's a bedtools file
        if (!is.null(input$sourceFile) && nzchar(input$sourceFile) && grepl("bedtools", input$sourceFile, ignore.case = TRUE)) {
          bedtools_file <- input$sourceFile
          cat("DEBUG: Using sourceFile as bedtools file:", bedtools_file, "\n")
        } else {
          cat("DEBUG: No bedtools file selected\n")
          return(NULL)
        }
      } else {
        cat("DEBUG: Selected bedtools file:", bedtools_file, "\n")
      }
      
      # Get the bedtools file path
      dir_path <- effective_dir()
      if (is.na(dir_path) || !dir.exists(dir_path)) {
        return(NULL)
      }
      
              # Find the bedtools file
        candidate_path <- file.path(dir_path, bedtools_file)
        if (!file.exists(candidate_path)) {
          hits <- Sys.glob(file.path(dir_path, "**", bedtools_file))
          if (length(hits) >= 1) candidate_path <- hits[1]
        }
      
      if (!file.exists(candidate_path)) {
        return(NULL)
      }
      
      # Load the bedtools file and create clustering data
      tryCatch({
        # Parse the bedtools file
        sim_matrix <- parse_bedtools_format(candidate_path)
        
        # Get accession key for labels
        ak <- effective_acc()
        labels_map <- NULL
        if (!is.na(ak) && file.exists(ak)) {
          labels_map <- tryCatch(load_accession_key(ak), error = function(e) NULL)
        }
        
                # Create clustering data for different numbers of clusters and selected linkage method
        clusters_vec <- 2:30
        linkage_method <- input$linkageMethod %||% "average"
        linkage_methods <- c(linkage_method)
        
        # Use the evaluate_nmi_long_table function
        long_table <- evaluate_nmi_long_table(sim_matrix, labels_map, clusters_vec, linkage_methods)
      
      cat("DEBUG: BEDTOOLS long_table created with", nrow(long_table), "rows\n")
      cat("DEBUG: Columns:", paste(names(long_table), collapse = ", "), "\n")
      
              # Add metadata columns
        long_table$source_file <- bedtools_file
      long_table$mode <- "BEDTOOLS"
      long_table$precision <- NA
      long_table$expA <- NA
      long_table$subA <- NA
      long_table$subB <- NA
      long_table$window <- NA
      long_table$klen <- NA
      long_table$sketch <- NA
      
      cat("DEBUG: Returning BEDTOOLS data with", nrow(long_table), "rows\n")
      return(long_table)
        
      }, error = function(e) {
        cat(sprintf("Error processing bedtools file %s: %s\n", candidate_path, conditionMessage(e)))
        return(NULL)
      })
    }
    
    # Regular mode handling
    # Use base filtered data and apply technical filtering
    d <- base_filtered_data()
    if (is.null(d)) return(NULL)
    
    # Apply technical parameter filtering
    if (!is.null(input$modeSelect) && nzchar(input$modeSelect)) {
      cat("DEBUG: Filtering by mode:", input$modeSelect, "\n")
      cat("DEBUG: Available modes:", unique(d$mode), "\n")
      d <- d[d$mode == input$modeSelect, , drop = FALSE]
      cat("DEBUG: After mode filtering:", nrow(d), "rows\n")
    }
    if (!is.null(input$linkageMethod) && nzchar(input$linkageMethod)) d <- d[d$linkage == input$linkageMethod, , drop = FALSE]
    if (!is.null(input$precision) && nzchar(input$precision)) {
      suppressWarnings(pv <- as.numeric(input$precision)); if (!is.na(pv)) d <- d[!is.na(d$precision) & d$precision == pv, , drop = FALSE]
    }
    if (identical(input$modeSelect, "C") && "expA" %in% names(d) && !is.null(input$expA) && nzchar(input$expA)) {
      suppressWarnings(v <- as.numeric(input$expA))
      if (!is.na(v)) d <- d[!is.na(d$expA) & abs(d$expA - v) < 1e-9, , drop = FALSE]
    }
    if (identical(input$modeSelect, "C") && "subB" %in% names(d) && !is.null(input$subB) && nzchar(input$subB)) {
      suppressWarnings(v <- as.numeric(input$subB))
      if (!is.na(v)) d <- d[!is.na(d$subB) & abs(d$subB - v) < 1e-9, , drop = FALSE]
    }
    if (identical(input$modeSelect, "D")) {
      if ("klen" %in% names(d) && !is.null(input$klen) && nzchar(input$klen)) { suppressWarnings(kv <- as.numeric(input$klen)); if (!is.na(kv)) d <- d[!is.na(d$klen) & d$klen == kv, , drop = FALSE] }
      if ("window" %in% names(d) && !is.null(input$window) && nzchar(input$window)) { suppressWarnings(wv <- as.numeric(input$window)); if (!is.na(wv)) d <- d[!is.na(d$window) & d$window == wv, , drop = FALSE] }
    }
    if (nrow(d) == 0) return(NULL)
    d
  })

  output$nmiPlot <- renderPlotly({
    # Safe graphics device handling for RStudio
    tryCatch({
      # Ensure we have a clean graphics state
      if (exists(".rs.restartR") && dev.cur() > 1) {
        while (dev.cur() > 1) try(dev.off(), silent = TRUE)
      }
      
      d <- filtered_data()
      cat("DEBUG: NMI plot - filtered_data returned", if(is.null(d)) "NULL" else paste(nrow(d), "rows"), "\n")
      if (!is.null(d)) {
        cat("DEBUG: NMI plot - columns:", paste(names(d), collapse = ", "), "\n")
        cat("DEBUG: NMI plot - mode:", unique(d$mode), "\n")
      }
      validate(need(!is.null(d), "No data to plot with the current filters."))
      cat("DEBUG: NMI plot - before filtering:", nrow(d), "rows\n")
      d <- d[!is.na(d$n_clusters) & !is.na(d$nmi), , drop = FALSE]
      cat("DEBUG: NMI plot - after filtering:", nrow(d), "rows\n")
      validate(need(nrow(d) > 0, "No rows with NMI and cluster numbers available."))
      limits <- silhouette_limits()
      
      # Subtitle summarizing parameters
      subtitle <- paste(
        c(
          if (!is.null(input$modeSelect) && nzchar(input$modeSelect)) paste0("mode ", input$modeSelect) else NULL,
          if (!is.null(input$linkageMethod) && nzchar(input$linkageMethod)) paste0("linkage ", input$linkageMethod) else NULL,
          if (!is.null(input$precision) && nzchar(input$precision)) paste0("p=", input$precision) else NULL,
          if (identical(input$modeSelect, "C") && !is.null(input$expA) && nzchar(input$expA)) paste0("expA=", input$expA) else NULL,
          if (identical(input$modeSelect, "C") && !is.null(input$subB) && nzchar(input$subB)) paste0("subB=", input$subB) else NULL,
          if (identical(input$modeSelect, "D") && !is.null(input$klen) && nzchar(input$klen)) paste0("k=", input$klen) else NULL,
          if (identical(input$modeSelect, "D") && !is.null(input$window) && nzchar(input$window)) paste0("w=", input$window) else NULL
        ), collapse = ", "
      )
      
      p <- ggplot(d, aes(x = n_clusters, y = nmi, color = silhouette, text = paste0(
          "clusters: ", n_clusters,
          "<br>NMI: ", sprintf("%.4f", nmi),
          "<br>silhouette: ", sprintf("%.4f", silhouette),
          "<br>linkage: ", linkage
        ))) +
        geom_point(size = 2.5, alpha = 0.9) +
        scale_color_viridis_c(option = "viridis", na.value = "grey80", limits = limits, oob = scales::squish) +
        labs(x = "Number of clusters", y = "NMI", color = "Silhouette",
             title = "NMI vs #Clusters",
             subtitle = subtitle) +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor = element_blank())
      
      ggplotly(p, tooltip = "text") %>% layout(hovermode = "closest")
      
    }, error = function(e) {
      # If graphics fail, show a simple plotly with error message
      if (grepl("evaluation nested too deeply|RStudioGD", conditionMessage(e))) {
        plotly::plot_ly() %>% 
          plotly::add_text(x = 0.5, y = 0.5, text = "Graphics device error - restart RStudio session") %>%
          plotly::layout(title = "Graphics Error", showlegend = FALSE)
      } else {
        stop(e)  # Re-throw other errors
      }
    })
  })

  # Biological summary table: counts by biosample and lifestage
  output$biologicalSummaryTable <- renderTable({
    acc_meta <- accession_metadata()
    sim_cache <- cached_similarity_matrices()
    
    if (is.null(acc_meta) || is.null(acc_meta$full_data) || is.null(sim_cache)) {
      return(data.frame(Message = "Load data to view biological summary"))
    }
    
    acc_df <- acc_meta$full_data
    if (!"file_base" %in% names(acc_df)) {
      return(data.frame(Message = "No file mapping available"))
    }
    
    if (length(sim_cache) == 0) {
      return(data.frame(Message = "No similarity matrices cached"))
    }
    
    # Simple approach: just count the total samples across all files
    all_samples <- c()
    for (src_file in names(sim_cache)) {
      sim_matrix <- sim_cache[[src_file]]
      if (!is.null(sim_matrix)) {
        all_samples <- c(all_samples, rownames(sim_matrix))
      }
    }
    
    if (length(all_samples) == 0) {
      return(data.frame(Message = "No samples found in similarity matrices"))
    }
    
    # Get biological metadata for all samples
    sample_meta_list <- list()
    for (sample in all_samples) {
      meta_row <- acc_df[acc_df$file_base == sample, , drop = FALSE]
      if (nrow(meta_row) > 0) {
        biosample <- if ("Biosample_term_name" %in% names(meta_row)) {
          meta_row$Biosample_term_name[1]
        } else {
          "unknown"
        }
        
        lifestage <- if ("Life_stage" %in% names(meta_row)) {
          meta_row$Life_stage[1]
        } else {
          "unknown"
        }
        
        sample_meta_list[[sample]] <- data.frame(
          sample = sample,
          biosample = biosample,
          lifestage = lifestage,
          stringsAsFactors = FALSE
        )
      }
    }
    
    if (length(sample_meta_list) == 0) {
      return(data.frame(Message = "No biological metadata found"))
    }
    
    # Combine all metadata
    all_meta <- do.call(rbind, sample_meta_list)
    
    # Count by biosample and lifestage
    bio_counts <- table(all_meta$biosample, all_meta$lifestage, useNA = "ifany")
    
    # Convert to data frame
    result_list <- list()
    for (i in 1:nrow(bio_counts)) {
      for (j in 1:ncol(bio_counts)) {
        if (bio_counts[i, j] > 0) {
          biosample <- rownames(bio_counts)[i]
          lifestage <- colnames(bio_counts)[j]
          count <- bio_counts[i, j]
          
          key <- paste(biosample, lifestage, sep = "_")
          result_list[[key]] <- data.frame(
            Biosample = biosample,
            Life_Stage = lifestage,
            Count = count,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    if (length(result_list) == 0) {
      return(data.frame(Message = "No biological combinations found"))
    }
    
    result <- do.call(rbind, result_list)
    rownames(result) <- NULL
    
    # Sort by biosample then lifestage
    result <- result[order(result$Biosample, result$Life_Stage), ]
    
    # Add totals row
    totals_row <- data.frame(
      Biosample = "TOTAL",
      Life_Stage = "",
      Count = sum(result$Count),
      stringsAsFactors = FALSE
    )
    
    result <- rbind(result, totals_row)
    
    # Clean up column names for display
    names(result) <- c("Biosample", "Life Stage", "Samples")
    
    return(result)
  }, rownames = FALSE)

  # Filtered data summary table: shows data based on current biological filter selections
  output$filteredDataSummaryTable <- renderTable({
    acc_meta <- accession_metadata()
    sim_cache <- cached_similarity_matrices()
    
    if (is.null(acc_meta) || is.null(acc_meta$full_data) || is.null(sim_cache)) {
      return(data.frame(Message = "Load data to view filtered summary"))
    }
    
    acc_df <- acc_meta$full_data
    if (!"file_base" %in% names(acc_df)) {
      return(data.frame(Message = "No file mapping available"))
    }
    
    if (length(sim_cache) == 0) {
      return(data.frame(Message = "No similarity matrices cached"))
    }
    
    # Get effective biosample selection based on mode
    effective_biosamples <- if (input$biosampleMode == "ALL") {
      acc_meta$biosamples
    } else if (input$biosampleMode == "REGULAR") {
      # For REGULAR mode, use the selected checkboxes directly
      # Add NULL check to prevent errors when no regular tissues are available
      if (!is.null(input$biosampleCheckboxes) && length(input$biosampleCheckboxes) > 0) {
        input$biosampleCheckboxes
      } else {
        # Fallback to all biosamples if no checkboxes are selected
        acc_meta$biosamples
      }
    } else {
      # For SUBSET mode, use approved_biosamples
      approved_biosamples()
    }
    
    # Get effective lifestage selection
    effective_lifestage <- if (is.null(input$lifestageSelect) || input$lifestageSelect == "ALL") {
      acc_meta$lifestages
    } else {
      input$lifestageSelect
    }
    
    # Filter samples based on biological selections
    filtered_samples <- c()
    for (src_file in names(sim_cache)) {
      sim_matrix <- sim_cache[[src_file]]
      if (!is.null(sim_matrix)) {
        file_samples <- rownames(sim_matrix)
        
        # Apply biological filtering to samples
        for (sample in file_samples) {
          meta_row <- acc_df[acc_df$file_base == sample, , drop = FALSE]
          if (nrow(meta_row) > 0) {
            biosample <- if ("Biosample_term_name" %in% names(meta_row)) {
              meta_row$Biosample_term_name[1]
            } else {
              "unknown"
            }
            
            lifestage <- if ("Life_stage" %in% names(meta_row)) {
              meta_row$Life_stage[1]
            } else {
              "unknown"
            }
            
            # Check if sample passes biological filters
            biosample_match <- biosample %in% effective_biosamples
            lifestage_match <- lifestage %in% effective_lifestage
            
            if (biosample_match && lifestage_match) {
              filtered_samples <- c(filtered_samples, sample)
            }
          }
        }
      }
    }
    
    if (length(filtered_samples) == 0) {
      return(data.frame(Message = "No samples match current biological filter selections"))
    }
    
    # Get biological metadata for filtered samples
    sample_meta_list <- list()
    for (sample in filtered_samples) {
      meta_row <- acc_df[acc_df$file_base == sample, , drop = FALSE]
      if (nrow(meta_row) > 0) {
        biosample <- if ("Biosample_term_name" %in% names(meta_row)) {
          meta_row$Biosample_term_name[1]
        } else {
          "unknown"
        }
        
        lifestage <- if ("Life_stage" %in% names(meta_row)) {
          meta_row$Life_stage[1]
        } else {
          "unknown"
        }
        
        sample_meta_list[[sample]] <- data.frame(
          sample = sample,
          biosample = biosample,
          lifestage = lifestage,
          stringsAsFactors = FALSE
        )
      }
    }
    
    if (length(sample_meta_list) == 0) {
      return(data.frame(Message = "No biological metadata found for filtered samples"))
    }
    
    # Combine all metadata
    all_meta <- do.call(rbind, sample_meta_list)
    
    # Count by biosample and lifestage
    bio_counts <- table(all_meta$biosample, all_meta$lifestage, useNA = "ifany")
    
    # Convert to data frame
    result_list <- list()
    for (i in 1:nrow(bio_counts)) {
      for (j in 1:ncol(bio_counts)) {
        if (bio_counts[i, j] > 0) {
          biosample <- rownames(bio_counts)[i]
          lifestage <- colnames(bio_counts)[j]
          count <- bio_counts[i, j]
          
          key <- paste(biosample, lifestage, sep = "_")
          result_list[[key]] <- data.frame(
            Biosample = biosample,
            Life_Stage = lifestage,
            Count = count,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    if (length(result_list) == 0) {
      return(data.frame(Message = "No biological combinations found in filtered data"))
    }
    
    result <- do.call(rbind, result_list)
    rownames(result) <- NULL
    
    # Sort by biosample then lifestage
    result <- result[order(result$Biosample, result$Life_Stage), ]
    
    # Add totals row
    totals_row <- data.frame(
      Biosample = "TOTAL",
      Life_Stage = "",
      Count = sum(result$Count),
      stringsAsFactors = FALSE
    )
    
    result <- rbind(result, totals_row)
    
    # Clean up column names for display
    names(result) <- c("Biosample", "Life Stage", "Samples")
    
    return(result)
  }, rownames = FALSE)

  # Filter summary outputs for all tabs
  output$filterSummaryInputs <- renderText({
    # For Inputs tab, only show biological filters and data summary
    summary_parts <- c()
    
    # Biological filters
    bio_summary <- biological_filter_summary()
    summary_parts <- c(summary_parts, bio_summary)
    
    # Data summary
    data_summary_text <- data_summary()
    if (nzchar(data_summary_text)) {
      summary_parts <- c(summary_parts, data_summary_text)
    }
    
    if (length(summary_parts) == 0) {
      return("No biological filters applied")
    }
    
    paste(summary_parts, collapse = " | ")
  })
  
  output$filterSummaryResults <- renderText({
    # Results tab: show all filters (technical + biological + data)
    filter_summary()
  })
  
  output$filterSummaryHeatmaps <- renderText({
    # Heatmaps tab: show only biological filters and data summary
    summary_parts <- c()
    
    # Biological filters
    bio_summary <- biological_filter_summary()
    summary_parts <- c(summary_parts, bio_summary)
    
    # Data summary
    data_summary_text <- data_summary()
    if (nzchar(data_summary_text)) {
      summary_parts <- c(summary_parts, data_summary_text)
    }
    
    if (length(summary_parts) == 0) {
      return("No biological filters applied")
    }
    
    paste(summary_parts, collapse = " | ")
  })
  
  output$filterSummaryDendrograms <- renderText({
    # Dendrograms tab: show only biological filters and data summary
    summary_parts <- c()
    
    # Biological filters
    bio_summary <- biological_filter_summary()
    summary_parts <- c(summary_parts, bio_summary)
    
    # Data summary
    data_summary_text <- data_summary()
    if (nzchar(data_summary_text)) {
      summary_parts <- c(summary_parts, data_summary_text)
    }
    
    if (length(summary_parts) == 0) {
      return("No biological filters applied")
    }
    
    paste(summary_parts, collapse = " | ")
  })
  
  output$filterSummaryPCA <- renderText({
    # PCA tab: show only biological filters and data summary
    summary_parts <- c()
    
    # Biological filters
    bio_summary <- biological_filter_summary()
    summary_parts <- c(summary_parts, bio_summary)
    
    # Data summary
    data_summary_text <- data_summary()
    if (nzchar(data_summary_text)) {
      summary_parts <- c(summary_parts, data_summary_text)
    }
    
    if (length(summary_parts) == 0) {
      return("No biological filters applied")
    }
    
    paste(summary_parts, collapse = " | ")
  })
  
  # Summary table: for each (source_file, linkage), show the row with the highest NMI
  # Uses biologically_filtered_data() to reflect only biological filtering (biosample/lifestage selection)
  output$summaryTable <- renderTable({
    tryCatch({
      df <- biologically_filtered_data()
      validate(need(!is.null(df), "Load data to view summary."))
      req_cols <- c("source_file", "linkage", "nmi", "n_clusters", "silhouette")
      if (!all(req_cols %in% names(df))) {
        return(data.frame())
      }
      # Ensure numeric types for metrics
      suppressWarnings({
        df$nmi <- as.numeric(df$nmi)
        df$n_clusters <- as.numeric(df$n_clusters)
        df$silhouette <- as.numeric(df$silhouette)
      })
      df <- df[is.finite(df$nmi), , drop = FALSE]
      if (!nrow(df)) return(data.frame())
      # group by file and linkage
      df$key <- paste(df$source_file, df$linkage, sep = "\t")
      parts <- split(df, df$key)
      pick_top <- function(d) {
        d <- d[order(-d$nmi, -d$silhouette, d$n_clusters), , drop = FALSE]
        d[1, c("source_file", "linkage", "nmi", "n_clusters", "silhouette"), drop = FALSE]
      }
      best_list <- lapply(parts, pick_top)
      best <- do.call(rbind, best_list)
      best <- best[order(-best$nmi, -best$silhouette), , drop = FALSE]
      rownames(best) <- NULL
      names(best) <- c("File", "Linkage", "NMI", "Clusters", "Silhouette")
      best
    }, error = function(e) {
      data.frame(Error = paste("Error:", conditionMessage(e)))
    })
  }, rownames = FALSE, digits = 4)

  # Download handler for summary table
  output$downloadSummaryTable <- downloadHandler(
    filename = function() {
      # Get the filename from input, default to "top_nmi_results" if empty
      filename <- input$downloadFilename
      if (is.null(filename) || nchar(trimws(filename)) == 0) {
        filename <- "top_nmi_results"
      }
      # Ensure filename doesn't contain invalid characters
      filename <- gsub("[^A-Za-z0-9_-]", "_", filename)
      paste0(filename, ".csv")
    },
    content = function(file) {
      tryCatch({
        # Get the same data that's displayed in the table
        df <- biologically_filtered_data()
        if (is.null(df)) {
          # Create empty data frame with proper structure if no data
          empty_df <- data.frame(
            File = character(0),
            Linkage = character(0),
            NMI = numeric(0),
            Clusters = numeric(0),
            Silhouette = numeric(0)
          )
          write.csv(empty_df, file, row.names = FALSE)
          return()
        }
        
        req_cols <- c("source_file", "linkage", "nmi", "n_clusters", "silhouette")
        if (!all(req_cols %in% names(df))) {
          # Create empty data frame with proper structure if columns missing
          empty_df <- data.frame(
            File = character(0),
            Linkage = character(0),
            NMI = numeric(0),
            Clusters = numeric(0),
            Silhouette = numeric(0)
          )
          write.csv(empty_df, file, row.names = FALSE)
          return()
        }
        
        # Process data the same way as the table
        suppressWarnings({
          df$nmi <- as.numeric(df$nmi)
          df$n_clusters <- as.numeric(df$n_clusters)
          df$silhouette <- as.numeric(df$silhouette)
        })
        df <- df[is.finite(df$nmi), , drop = FALSE]
        
        if (nrow(df) == 0) {
          # Create empty data frame with proper structure if no valid data
          empty_df <- data.frame(
            File = character(0),
            Linkage = character(0),
            NMI = numeric(0),
            Clusters = numeric(0),
            Silhouette = numeric(0)
          )
          write.csv(empty_df, file, row.names = FALSE)
          return()
        }
        
        # Group by file and linkage, pick top NMI for each
        df$key <- paste(df$source_file, df$linkage, sep = "\t")
        parts <- split(df, df$key)
        pick_top <- function(d) {
          d <- d[order(-d$nmi, -d$silhouette, d$n_clusters), , drop = FALSE]
          d[1, c("source_file", "linkage", "nmi", "n_clusters", "silhouette"), drop = FALSE]
        }
        best_list <- lapply(parts, pick_top)
        best <- do.call(rbind, best_list)
        best <- best[order(-best$nmi, -best$silhouette), , drop = FALSE]
        rownames(best) <- NULL
        names(best) <- c("File", "Linkage", "NMI", "Clusters", "Silhouette")
        
        # Write to CSV
        write.csv(best, file, row.names = FALSE)
      }, error = function(e) {
        # Create error data frame if something goes wrong
        error_df <- data.frame(
          Error = paste("Error generating download:", conditionMessage(e))
        )
        write.csv(error_df, file, row.names = FALSE)
      })
    }
  )

  # Optimized: Shared biological filter logic
  biological_filter_summary <- reactive({
    bio_filters <- c()
    
    # Biosample filtering
    if (!is.null(input$biosampleMode) && input$biosampleMode != "ALL") {
      approved_bs <- approved_biosamples()
      if (!is.null(approved_bs) && length(approved_bs) > 0) {
        bio_filters <- c(bio_filters, paste("Biosamples:", paste(approved_bs, collapse = ", ")))
      } else {
        bio_filters <- c(bio_filters, "Biosamples: Custom subset")
      }
    } else {
      bio_filters <- c(bio_filters, "Biosamples: All")
    }
    
    # Lifestage filtering
    if (!is.null(input$lifestageSelect) && input$lifestageSelect != "ALL") {
      bio_filters <- c(bio_filters, paste("Life Stage:", input$lifestageSelect))
    } else {
      bio_filters <- c(bio_filters, "Life Stage: All")
    }
    
    if (length(bio_filters) == 0) {
      return("Biological: All")
    }
    
    paste("Biological:", paste(bio_filters, collapse = ", "))
  })
  
  # Optimized: Shared technical filter logic
  technical_filter_summary <- reactive({
    tech_filters <- c()
    if (!is.null(input$modeSelect) && nzchar(input$modeSelect)) {
      tech_filters <- c(tech_filters, paste("Mode:", input$modeSelect))
    }
    if (!is.null(input$linkageMethod) && nzchar(input$linkageMethod)) {
      tech_filters <- c(tech_filters, paste("Linkage:", input$linkageMethod))
    }
    if (!is.null(input$precision) && nzchar(input$precision)) {
      tech_filters <- c(tech_filters, paste("Precision:", input$precision))
    }
    if (identical(input$modeSelect, "C")) {
      if (!is.null(input$expA) && nzchar(input$expA)) {
        tech_filters <- c(tech_filters, paste("expA:", input$expA))
      }
      if (!is.null(input$subB) && nzchar(input$subB)) {
        tech_filters <- c(tech_filters, paste("subB:", input$subB))
      }
    }
    if (identical(input$modeSelect, "D")) {
      if (!is.null(input$klen) && nzchar(input$klen)) {
        tech_filters <- c(tech_filters, paste("k:", input$klen))
      }
      if (!is.null(input$window) && nzchar(input$window)) {
        tech_filters <- c(tech_filters, paste("w:", input$window))
      }
    }
    
    if (length(tech_filters) == 0) {
      return("")
    }
    
    paste("Technical:", paste(tech_filters, collapse = ", "))
  })
  
  # Optimized: Shared data summary logic
  data_summary <- reactive({
    df <- base_filtered_data()
    if (!is.null(df) && nrow(df) > 0) {
      unique_files <- length(unique(df$source_file))
      total_rows <- nrow(df)
      paste("Data:", unique_files, "files,", total_rows, "rows")
    } else {
      ""
    }
  })
  
  # Reactive function to generate filter summary text
  filter_summary <- reactive({
    # Handle BEDTOOLS mode
    if (!is.null(input$modeSelect) && input$modeSelect == "BEDTOOLS") {
      summary_parts <- c("Mode: BEDTOOLS")
      
      if (!is.null(input$bedtoolsFile) && nzchar(input$bedtoolsFile)) {
        summary_parts <- c(summary_parts, paste("File:", input$bedtoolsFile))
      }
      
      # Biological filters (still apply to bedtools data)
      bio_summary <- biological_filter_summary()
      summary_parts <- c(summary_parts, bio_summary)
      
      # Data summary
      data_summary_text <- data_summary()
      if (nzchar(data_summary_text)) {
        summary_parts <- c(summary_parts, data_summary_text)
      }
      
      return(paste(summary_parts, collapse = " | "))
    }
    
    # Regular mode
    summary_parts <- c()
    
    # Technical filters
    tech_summary <- technical_filter_summary()
    if (nzchar(tech_summary)) {
      summary_parts <- c(summary_parts, tech_summary)
    }
    
    # Biological filters
    bio_summary <- biological_filter_summary()
    summary_parts <- c(summary_parts, bio_summary)
    
    # Data summary
    data_summary_text <- data_summary()
    if (nzchar(data_summary_text)) {
      summary_parts <- c(summary_parts, data_summary_text)
    }
    
    if (length(summary_parts) == 0) {
      return("No filters applied")
    }
    
    paste(summary_parts, collapse = " | ")
  })
  
  # Auto-select bedtools file when BEDTOOLS mode is selected
  observeEvent(input$modeSelect, {
    if (!is.null(input$modeSelect) && input$modeSelect == "BEDTOOLS") {
      # Get available bedtools files
      df <- data_all()
      if (!is.null(df) && nrow(df) > 0) {
        srcs <- sort(unique(df$source_file))
        bedtools_files <- srcs[grepl("bedtools", srcs, ignore.case = TRUE)]
        if (length(bedtools_files) > 0) {
          updateSelectInput(session, "sourceFile", selected = bedtools_files[1])
        }
      }
    }
  })
  
  # Update precision choices when mode changes
  observeEvent(input$modeSelect, {
    if (!is.null(input$modeSelect) && input$modeSelect != "BEDTOOLS") {
      df <- data_all()
      if (!is.null(df) && nrow(df) > 0) {
        # Get precision values for the selected mode
        mode_precisions <- sort(unique(na.omit(df$precision[df$mode == input$modeSelect])))
        if (length(mode_precisions) > 0) {
          updateSelectInput(session, "precision", choices = mode_precisions, selected = mode_precisions[1])
        } else {
          updateSelectInput(session, "precision", choices = character(0))
        }
      }
    }
  })
  
  # Auto-select a matching file for the dendrogram when NMI parameters change
  observeEvent(list(input$modeSelect, input$linkageMethod, input$precision,
                    input$expA, input$subB, input$klen, input$window, filtered_data()), {
    d <- filtered_data()
    if (is.null(d) || !("source_file" %in% names(d))) return()
    srcs <- sort(unique(na.omit(as.character(d$source_file))))
    if (!length(srcs)) return()
    current <- input$sourceFile %||% ""
    if (nzchar(current) && current %in% srcs) return()
    # Prefer file with highest mean silhouette, then highest mean NMI, otherwise first
    safe_mean <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      m <- tryCatch(mean(x, na.rm = TRUE), error = function(e) NA_real_)
      if (is.infinite(m)) NA_real_ else m
    }
    sil_means <- tryCatch(tapply(d$silhouette, d$source_file, safe_mean), error = function(e) NULL)
    best <- NULL
    if (!is.null(sil_means) && length(sil_means)) {
      if (!all(is.na(sil_means))) best <- names(sil_means)[which.max(sil_means)]
    }
    if (is.null(best)) {
      nmi_means <- tryCatch(tapply(d$nmi, d$source_file, safe_mean), error = function(e) NULL)
      if (!is.null(nmi_means) && length(nmi_means) && !all(is.na(nmi_means))) {
        best <- names(nmi_means)[which.max(nmi_means)]
      }
    }
    if (is.null(best) || !nzchar(best)) best <- srcs[1]
    if (!identical(best, current)) updateSelectInput(session, "sourceFile", selected = best)
  }, ignoreInit = TRUE)

  output$dendPlot <- renderPlot({
    sf <- input$sourceFile
    lk <- input$dendLinkage %||% "average"
    validate(need(!is.null(sf) && nzchar(sf), "Select an input file for dendrogram."))
    
    # find matching original file path by name from current directory
    base_dir <- effective_dir()
    candidate <- file.path(base_dir, sf)
    if (!file.exists(candidate)) {
      # fallback: search recursively under base_dir
      hits <- Sys.glob(file.path(base_dir, "**", sf))
      if (length(hits) >= 1) candidate <- hits[1]
    }
    if (!file.exists(candidate)) {
      plot.new(); title(sprintf("Cannot find file: %s", sf)); return(invisible())
    }
    
    # Handle BEDTOOLS mode specially
    if (!is.null(input$modeSelect) && input$modeSelect == "BEDTOOLS") {
      cat("DEBUG: Dendrogram - BEDTOOLS mode detected\n")
      # For BEDTOOLS mode, we need to parse the file and create dendrogram directly
      tryCatch({
        # Parse the bedtools file
        sim_matrix <- parse_bedtools_format(candidate)
        
        # Get accession key for labels
        ak <- effective_acc()
        labels_map <- NULL
        if (!is.na(ak) && file.exists(ak)) {
          labels_map <- tryCatch(load_accession_key(ak), error = function(e) NULL)
        }
        
        # Create dendrogram using the same logic as draw_dendrogram_for but for bedtools
        draw_dendrogram_from_matrix(sim_matrix, labels_map, lk, "BEDTOOLS")
        
      }, error = function(e) {
        plot.new(); title(sprintf("Error processing bedtools file: %s", conditionMessage(e)))
      })
    } else {
      # Regular mode - use existing function
      draw_dendrogram_for(candidate, effective_acc(), lk)
    }
  })

  # Paired dendrograms tab
  output$dendPlotA <- renderPlot({
    df <- filtered_data()
    validate(need(!is.null(df), "Load data to view dendrograms."))
    sf <- input$dendFile1
    lk <- input$dendLinkage1 %||% "average"
    validate(need(!is.null(sf) && nzchar(sf), "Select an input file for Dendrogram A."))
    base_dir <- effective_dir()
    candidate <- file.path(base_dir, sf)
    if (!file.exists(candidate)) {
      hits <- Sys.glob(file.path(base_dir, "**", sf))
      if (length(hits) >= 1) candidate <- hits[1]
    }
    if (!file.exists(candidate)) { plot.new(); title(sprintf("Cannot find file: %s", sf)); return(invisible()) }
    draw_dendrogram_for(candidate, effective_acc(), lk)
  })

  output$dendPlotB <- renderPlot({
    df <- filtered_data()
    validate(need(!is.null(df), "Load data to view dendrograms."))
    sf <- input$dendFile2
    lk <- input$dendLinkage2 %||% "complete"
    validate(need(!is.null(sf) && nzchar(sf), "Select an input file for Dendrogram B."))
    base_dir <- effective_dir()
    candidate <- file.path(base_dir, sf)
    if (!file.exists(candidate)) {
      hits <- Sys.glob(file.path(base_dir, "**", sf))
      if (length(hits) >= 1) candidate <- hits[1]
    }
    if (!file.exists(candidate)) { plot.new(); title(sprintf("Cannot find file: %s", sf)); return(invisible()) }
    draw_dendrogram_for(candidate, effective_acc(), lk)
  })

  # Heatmap info panel
  output$heatmapInfo <- renderText({
    df <- filtered_data()
    if (is.null(df)) return("Load data to view heatmap details.")
    files <- c(input$heatmapFile1, input$heatmapFile2)
    files <- files[nzchar(files)]
    if (!length(files)) return("Select files to view heatmap details.")
    # Describe formats and approach
    base_dir <- effective_dir()
    fmt_for <- function(sf) {
      candidate <- file.path(base_dir, sf)
      if (!file.exists(candidate)) {
        hits <- Sys.glob(file.path(base_dir, "**", sf))
        if (length(hits) >= 1) candidate <- hits[1]
      }
      if (file.exists(candidate)) {
        tryCatch({
          f <- detect_file_format(candidate)
          paste0(sf, ": ", if (f == 'hammock') "hammock similarity matrix (jaccard_similarity[[_with_ends]])" else if (f == 'bedtools') "bedtools jaccard matrix" else f)
        }, error = function(e) paste0(sf, ": unknown format"))
      } else paste0(sf, ": not found")
    }
    fmts <- vapply(files, fmt_for, character(1))
    paste(c(
      "Heatmap rendering: ggplot2 geom_raster with viridis scale (0..1 similarity).",
      "Side color bars: tissue (left) and organism (right) derived from accession key.",
      paste0("Selected files:", if (length(fmts)) paste0("\n - ", fmts, collapse = "") else " none"),
      "Dendrogram linkage for clustering visuals: follows the selected linkage control."
    ), collapse = "\n")
  })

  output$heatmap1 <- renderPlot({
    df <- filtered_data()
    validate(need(!is.null(df), "Load data to view heatmaps."))
    sf <- input$heatmapFile1
    validate(need(!is.null(sf) && nzchar(sf), "Select an input file for Heatmap A."))
    base_dir <- effective_dir()
    candidate <- file.path(base_dir, sf)
    if (!file.exists(candidate)) {
      hits <- Sys.glob(file.path(base_dir, "**", sf))
      if (length(hits) >= 1) candidate <- hits[1]
    }
    if (!file.exists(candidate)) { plot.new(); title(sprintf("Cannot find file: %s", sf)); return(invisible()) }
    draw_heatmap_for(candidate, effective_acc())
  })

  output$heatmap2 <- renderPlot({
    df <- filtered_data()
    validate(need(!is.null(df), "Load data to view heatmaps."))
    sf <- input$heatmapFile2
    validate(need(!is.null(sf) && nzchar(sf), "Select an input file for Heatmap B."))
    base_dir <- effective_dir()
    candidate <- file.path(base_dir, sf)
    if (!file.exists(candidate)) {
      hits <- Sys.glob(file.path(base_dir, "**", sf))
      if (length(hits) >= 1) candidate <- hits[1]
    }
    if (!file.exists(candidate)) { plot.new(); title(sprintf("Cannot find file: %s", sf)); return(invisible()) }
    draw_heatmap_for(candidate, effective_acc())
  })
  
  # PCA computation helper function
  compute_pca_from_similarity <- function(sim_matrix, accession_meta) {
    if (is.null(sim_matrix) || nrow(sim_matrix) < 3) {
      return(NULL)
    }
    
    # Convert similarity to distance matrix
    # Similarity ranges from 0-1, so distance = 1 - similarity
    dist_matrix <- 1 - sim_matrix
    
    # Ensure it's a proper distance matrix (symmetric, zero diagonal)
    diag(dist_matrix) <- 0
    dist_matrix <- (dist_matrix + t(dist_matrix)) / 2
    
    # Convert to distance object and perform MDS (equivalent to PCA on distance data)
    tryCatch({
      # Use classical multidimensional scaling
      mds_result <- cmdscale(as.dist(dist_matrix), k = min(10, nrow(dist_matrix) - 1), eig = TRUE)
      
      if (is.null(mds_result) || is.null(mds_result$points)) {
        return(NULL)
      }
      
      # Create data frame with PC coordinates
      pca_data <- data.frame(
        Sample = rownames(sim_matrix),
        PC1 = mds_result$points[, 1],
        PC2 = mds_result$points[, 2],
        stringsAsFactors = FALSE
      )
      
      # Add more PCs if available
      if (ncol(mds_result$points) >= 3) {
        pca_data$PC3 <- mds_result$points[, 3]
      }
      
      # Add biological metadata if available
      if (!is.null(accession_meta) && !is.null(accession_meta$full_data)) {
        acc_df <- accession_meta$full_data
        
        for (i in 1:nrow(pca_data)) {
          sample_name <- pca_data$Sample[i]
          meta_row <- acc_df[acc_df$file_base == sample_name, , drop = FALSE]
          
          if (nrow(meta_row) > 0) {
            pca_data$Biosample[i] <- if ("Biosample_term_name" %in% names(meta_row)) {
              meta_row$Biosample_term_name[1]
            } else {
              "unknown"
            }
            
            pca_data$Life_Stage[i] <- if ("Life_stage" %in% names(meta_row)) {
              meta_row$Life_stage[1]
            } else {
              "unknown"
            }
            
            pca_data$Organism[i] <- if ("Organism" %in% names(meta_row)) {
              meta_row$Organism[1]
            } else {
              "unknown"
            }
          } else {
            pca_data$Biosample[i] <- "unknown"
            pca_data$Life_Stage[i] <- "unknown"
            pca_data$Organism[i] <- "unknown"
          }
        }
      } else {
        pca_data$Biosample <- "unknown"
        pca_data$Life_Stage <- "unknown"
        pca_data$Organism <- "unknown"
      }
      
      # Calculate variance explained
      eigenvalues <- mds_result$eig
      if (!is.null(eigenvalues)) {
        total_var <- sum(abs(eigenvalues))
        var_explained <- abs(eigenvalues) / total_var * 100
      } else {
        var_explained <- rep(NA, ncol(mds_result$points))
      }
      
      return(list(
        data = pca_data,
        var_explained = var_explained,
        method = "Classical MDS on similarity matrix"
      ))
      
    }, error = function(e) {
      cat("PCA computation error:", e$message, "\n")
      return(NULL)
    })
  }
  
  # PCA plot rendering function
  render_pca_plot <- function(file_name) {
    if (is.null(file_name) || file_name == "") {
      return(NULL)
    }
    
    # Get filtered similarity matrix
    sim_cache <- cached_similarity_matrices()
    acc_meta <- accession_metadata()
    
    if (is.null(sim_cache) || !file_name %in% names(sim_cache)) {
      return(NULL)
    }
    
    sim_matrix <- sim_cache[[file_name]]
    
    # Apply biological filtering if active
    biosample_mode <- input$biosampleMode
    approved_bs <- approved_biosamples()
    lifestage_sel <- input$lifestageSelect
    
    if (!is.null(biosample_mode) && biosample_mode == "SUBSET" && !is.null(approved_bs) && length(approved_bs) > 0) {
      # Filter the similarity matrix
      sim_matrix <- filter_similarity_matrix(sim_matrix, acc_meta, approved_bs, lifestage_sel)
    } else if (!is.null(lifestage_sel) && lifestage_sel != "ALL") {
      # Filter by lifestage only
      sim_matrix <- filter_similarity_matrix(sim_matrix, acc_meta, NULL, lifestage_sel)
    }
    
    if (is.null(sim_matrix) || nrow(sim_matrix) < 3) {
      return(NULL)
    }
    
    # Compute PCA
    pca_result <- compute_pca_from_similarity(sim_matrix, acc_meta)
    
    if (is.null(pca_result)) {
      return(NULL)
    }
    
    pca_data <- pca_result$data
    var_exp <- pca_result$var_explained
    
    # Create plotly plot
    p <- plot_ly(
      data = pca_data,
      x = ~PC1, 
      y = ~PC2,
      color = ~Biosample,
      symbol = ~Life_Stage,
      text = ~paste("Sample:", Sample, "<br>Biosample:", Biosample, "<br>Life Stage:", Life_Stage),
      hovertemplate = "%{text}<extra></extra>",
      type = "scatter",
      mode = "markers",
      marker = list(size = 8, line = list(width = 1, color = "black"))
    ) %>%
      layout(
        title = paste("PCA:", file_name, 
                     if (!is.null(var_exp) && length(var_exp) >= 2) {
                       paste0(" (PC1: ", round(var_exp[1], 1), "%, PC2: ", round(var_exp[2], 1), "%)")
                     } else {
                       ""
                     }),
        xaxis = list(title = paste("PC1", if (!is.null(var_exp)) paste0(" (", round(var_exp[1], 1), "%)") else "")),
        yaxis = list(title = paste("PC2", if (!is.null(var_exp) && length(var_exp) >= 2) paste0(" (", round(var_exp[2], 1), "%)") else "")),
        showlegend = TRUE,
        hovermode = "closest"
      )
    
    return(p)
  }
  
  # PCA plot outputs
  output$pcaPlot1 <- renderPlotly({
    render_pca_plot(input$pcaFile1)
  })
  
  output$pcaPlot2 <- renderPlotly({
    render_pca_plot(input$pcaFile2)
  })
  
  # PCA information output
  output$pcaInfo <- renderText({
    file1 <- input$pcaFile1
    file2 <- input$pcaFile2
    
    sim_cache <- cached_similarity_matrices()
    acc_meta <- accession_metadata()
    
    if (is.null(sim_cache) || is.null(acc_meta)) {
      return("Load data to view PCA information")
    }
    
    info_lines <- c()
    
    # Information about filtering
    biosample_mode <- input$biosampleMode
    approved_bs <- approved_biosamples()
    lifestage_sel <- input$lifestageSelect
    
    if (!is.null(biosample_mode) && biosample_mode == "SUBSET" && !is.null(approved_bs) && length(approved_bs) > 0) {
      info_lines <- c(info_lines, paste("Biosample filter active:", paste(approved_bs, collapse = ", ")))
    } else {
      info_lines <- c(info_lines, "Biosample filter: All biosamples")
    }
    
    if (!is.null(lifestage_sel) && lifestage_sel != "ALL") {
      info_lines <- c(info_lines, paste("Life stage filter:", lifestage_sel))
    } else {
      info_lines <- c(info_lines, "Life stage filter: All stages")
    }
    
    info_lines <- c(info_lines, "")
    
    # Information about each selected file
    for (file_name in c(file1, file2)) {
      if (!is.null(file_name) && file_name != "" && file_name %in% names(sim_cache)) {
        sim_matrix <- sim_cache[[file_name]]
        
        # Apply filtering
        if (!is.null(biosample_mode) && biosample_mode == "SUBSET" && !is.null(approved_bs) && length(approved_bs) > 0) {
          sim_matrix <- filter_similarity_matrix(sim_matrix, acc_meta, approved_bs, lifestage_sel)
        } else if (!is.null(lifestage_sel) && lifestage_sel != "ALL") {
          sim_matrix <- filter_similarity_matrix(sim_matrix, acc_meta, NULL, lifestage_sel)
        }
        
        if (!is.null(sim_matrix)) {
          info_lines <- c(info_lines, paste(file_name, ": ", nrow(sim_matrix), " samples"))
          
          # Count biological categories
          if (!is.null(acc_meta$full_data)) {
            sample_names <- rownames(sim_matrix)
            acc_df <- acc_meta$full_data
            
            bio_counts <- table(sapply(sample_names, function(s) {
              meta_row <- acc_df[acc_df$file_base == s, , drop = FALSE]
              if (nrow(meta_row) > 0 && "Biosample_term_name" %in% names(meta_row)) {
                meta_row$Biosample_term_name[1]
              } else {
                "unknown"
              }
            }))
            
            info_lines <- c(info_lines, paste("  Biosamples:", paste(names(bio_counts), " (", bio_counts, ")", sep = "", collapse = ", ")))
          }
        } else {
          info_lines <- c(info_lines, paste(file_name, ": No samples after filtering"))
        }
      }
    }
    
    info_lines <- c(info_lines, "", "Method: Classical multidimensional scaling (MDS) on similarity matrices")
    info_lines <- c(info_lines, "Distance = 1 - Similarity")
    
    return(paste(info_lines, collapse = "\n"))
  })
}

`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  # If vector, consider any non-empty element a hit
  if (length(a) > 1) {
    nz <- nzchar(as.character(a))
    if (any(nz)) return(a[which(nz)[1]])
    return(b)
  }
  if (is.na(a)) return(b)
  if (is.character(a) || is.numeric(a)) {
    return(if (nzchar(as.character(a))) a else b)
  }
  b
}

shinyApp(ui, server)
