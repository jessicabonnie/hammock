#!/usr/bin/env Rscript

# Utility functions for R scripts working with hammock/bedtools outputs

suppressWarnings({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    message("Note: data.table not installed; falling back to base readers.")
  }
})

get_basename <- function(filename) {
  name <- as.character(filename)
  if (grepl("\\.bed\\.gz$", name)) {
    sub("\\.bed\\.gz$", "", name)
  } else if (grepl("\\.fa\\.gz$", name)) {
    sub("\\.fa\\.gz$", "", name)
  } else {
    tools::file_path_sans_ext(basename(name))
  }
}

detect_file_format <- function(filepath) {
  # Check filename first - if "bedtools" appears in the filename, treat as bedtools format
  filename <- basename(filepath)
  if (grepl("bedtools", filename, ignore.case = TRUE)) {
    return('bedtools')
  }
  
  con <- file(filepath, "r"); on.exit(close(con))
  first_line <- readLines(con, n = 1, warn = FALSE)
  if (grepl('^file1\\tfile2\\tintersection\\tunion\\tjaccard', first_line) ||
      grepl('^file1 file2 intersection union jaccard', first_line)) {
    return('bedtools')
  }
  if (grepl('^file1,file2', first_line) && (grepl('jaccard_similarity', first_line) || grepl('jaccard_similarity_with_ends', first_line))) {
    return('hammock')
  }
  # Fallback: inspect columns
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- tryCatch(data.table::fread(filepath, nrows = 5, showProgress = FALSE), error = function(e) NULL)
  } else {
    dt <- tryCatch(utils::read.csv(filepath, nrows = 5, header = TRUE), error = function(e) NULL)
  }
  if (is.null(dt)) stop("Cannot determine file format for ", filepath)
  cols <- names(dt)
  if (all(c('jaccard', 'intersection') %in% cols)) return('bedtools')
  if (any(c('jaccard_similarity', 'jaccard_similarity_with_ends') %in% cols)) return('hammock')
  stop("Cannot determine file format for ", filepath)
}

parse_hammock_format <- function(filepath) {
  df <- if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(filepath, showProgress = FALSE)
  } else {
    utils::read.csv(filepath, stringsAsFactors = FALSE)
  }
  jcol <- if ('jaccard_similarity_with_ends' %in% names(df)) 'jaccard_similarity_with_ends' else if ('jaccard_similarity' %in% names(df)) 'jaccard_similarity' else stop('No jaccard column found')
  df$file1_base <- vapply(df$file1, get_basename, character(1))
  df$file2_base <- vapply(df$file2, get_basename, character(1))
  all_files <- sort(unique(c(df$file1_base, df$file2_base)))
  n <- length(all_files)
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(all_files, all_files))
  idx <- setNames(seq_along(all_files), all_files)
  for (i in seq_len(nrow(df))) {
    a <- df$file1_base[i]; b <- df$file2_base[i]; s <- as.numeric(df[[jcol]][i])
    ia <- idx[[a]]; ib <- idx[[b]]
    mat[ia, ib] <- s; mat[ib, ia] <- s
  }
  diag(mat) <- 1
  mat
}

parse_bedtools_format <- function(filepath) {
  # Always use read.table for bedtools files as it handles whitespace-separated files better
  # than fread when headers have multiple spaces between column names
  df <- utils::read.table(filepath, header = TRUE, stringsAsFactors = FALSE)
  if (!all(c('file1', 'file2', 'jaccard') %in% names(df))) stop('Bedtools file missing required columns: file1, file2, jaccard')
  df$file1_base <- vapply(df$file1, get_basename, character(1))
  df$file2_base <- vapply(df$file2, get_basename, character(1))
  all_files <- sort(unique(c(df$file1_base, df$file2_base)))
  n <- length(all_files)
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(all_files, all_files))
  idx <- setNames(seq_along(all_files), all_files)
  for (i in seq_len(nrow(df))) {
    a <- df$file1_base[i]; b <- df$file2_base[i]; s <- as.numeric(df$jaccard[i])
    ia <- idx[[a]]; ib <- idx[[b]]
    mat[ia, ib] <- s; mat[ib, ia] <- s
  }
  diag(mat) <- 1
  mat
}

load_accession_key <- function(filepath, include_life_stage = TRUE) {
  df <- if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(filepath, sep = "\t", header = TRUE, showProgress = FALSE)
  } else {
    utils::read.table(filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  }
  if (!all(c('File', 'Biosample_term_name') %in% names(df))) stop('Accession key missing required columns: File, Biosample_term_name')
  
  base <- vapply(df$File, function(x) tools::file_path_sans_ext(basename(x)), character(1))
  
  # Check if Life_stage column exists and include_life_stage is TRUE
  has_life_stage <- 'Life_stage' %in% names(df) && include_life_stage
  
  if (has_life_stage) {
    # Create enhanced labels combining Biosample + Life_stage
    enhanced_labels <- create_enhanced_labels(df$Biosample_term_name, df$Life_stage)
    result <- setNames(enhanced_labels, base)
    attr(result, "has_life_stage") <- TRUE
    attr(result, "raw_data") <- df
    return(result)
  } else {
    # Return traditional Biosample labels only
    result <- setNames(df$Biosample_term_name, base)
    attr(result, "has_life_stage") <- FALSE
    if (has_life_stage) {
      attr(result, "raw_data") <- df
    }
    return(result)
  }
}

create_enhanced_labels <- function(biosample_terms, life_stages) {
  # Clean up life stage values (handle NAs, empty strings)
  clean_life_stages <- ifelse(is.na(life_stages) | life_stages == "" | is.null(life_stages), 
                             "unknown", as.character(life_stages))
  
  # Create enhanced labels: "biosample (life_stage)"
  enhanced <- paste0(biosample_terms, " (", clean_life_stages, ")")
  
  return(enhanced)
}

get_biosample_from_enhanced_label <- function(enhanced_label) {
  # Extract biosample term from enhanced label "biosample (life_stage)"
  gsub("\\s*\\([^)]*\\)$", "", enhanced_label)
}

get_life_stage_from_enhanced_label <- function(enhanced_label) {
  # Extract life stage from enhanced label "biosample (life_stage)"
  matches <- regmatches(enhanced_label, regexpr("\\(([^)]*)\\)$", enhanced_label))
  if (length(matches) > 0) {
    gsub("[()]", "", matches)
  } else {
    "unknown"
  }
}

has_life_stage_info <- function(label_map) {
  # Check if the label map includes life stage information
  # Handle NULL values properly to avoid "missing value where TRUE/FALSE needed" error
  has_attr <- !is.null(attr(label_map, "has_life_stage"))
  if (!has_attr) return(FALSE)
  attr_value <- attr(label_map, "has_life_stage")
  if (is.null(attr_value) || is.na(attr_value)) return(FALSE)
  return(isTRUE(attr_value))
}

create_color_palette_with_life_stage <- function(enhanced_labels, base_palette_fun = NULL) {
  # Create a color palette that groups by biosample but varies by life stage
  if (is.null(base_palette_fun)) {
    base_palette_fun <- if (requireNamespace('scales', quietly = TRUE)) {
      scales::hue_pal()
    } else {
      function(n) grDevices::rainbow(n, s = 0.7, v = 0.9)
    }
  }
  
  # Extract unique biosamples and life stages
  biosamples <- unique(get_biosample_from_enhanced_label(enhanced_labels))
  life_stages <- unique(get_life_stage_from_enhanced_label(enhanced_labels))
  
  # Create base colors for biosamples
  base_colors <- base_palette_fun(length(biosamples))
  names(base_colors) <- biosamples
  
  # Create variations for life stages using different saturations/brightness
  life_stage_factors <- c("adult" = 1.0, "embryonic" = 0.7, "postnatal" = 0.85, 
                         "child" = 0.9, "unknown" = 0.6)
  
  enhanced_colors <- character(length(enhanced_labels))
  names(enhanced_colors) <- enhanced_labels
  
  for (label in enhanced_labels) {
    biosample <- get_biosample_from_enhanced_label(label)
    life_stage <- get_life_stage_from_enhanced_label(label)
    
    base_color <- base_colors[biosample]
    factor <- life_stage_factors[life_stage]
    if (is.na(factor)) factor <- 0.8
    
    # Adjust color brightness/saturation based on life stage
    rgb_vals <- col2rgb(base_color) / 255
    adjusted_vals <- rgb_vals * factor + (1 - factor) * 0.3  # blend with gray
    enhanced_colors[label] <- rgb(adjusted_vals[1], adjusted_vals[2], adjusted_vals[3])
  }
  
  return(enhanced_colors)
}

infer_label_from_basename <- function(basename) {
  token <- strsplit(basename, '-', fixed = TRUE)[[1]][1]
  if (startsWith(token, 'f') && nchar(token) > 1) token <- substring(token, 2)
  token
}

align_matrix_and_labels <- function(sim_df, accession_labels, quiet = FALSE) {
  files_in_matrix <- rownames(sim_df)
  common <- intersect(files_in_matrix, names(accession_labels))
  if (length(common) == 0) {
    if (!quiet) message('Warning: No overlap with accession key; inferring labels from basenames.')
    inferred <- setNames(vapply(files_in_matrix, infer_label_from_basename, character(1)), files_in_matrix)
    # Preserve attributes for consistency
    attr(inferred, "has_life_stage") <- FALSE
    list(sim_df, inferred)
  } else {
    sim_sub <- sim_df[common, common, drop = FALSE]
    # Subset the labels but preserve attributes
    labels_sub <- accession_labels[common]
    attributes_to_preserve <- attributes(accession_labels)
    attributes_to_preserve$names <- names(labels_sub)
    attributes(labels_sub) <- attributes_to_preserve
    list(sim_sub, labels_sub)
  }
}


# --- Additional helpers for parameter extraction / detection ---

detect_hammock_expA <- function(filepath) {
  head <- tryCatch(utils::read.csv(filepath, nrows = 10), error = function(e) NULL)
  if (is.null(head)) return(NA_real_)
  for (col in names(head)) {
    low <- tolower(col)
    if (low %in% c('expa', 'exp')) {
      vals <- unique(na.omit(head[[col]]))
      if (length(vals) > 0) {
        v <- suppressWarnings(as.numeric(vals[1]))
        if (!is.na(v)) return(v)
      }
    }
  }
  NA_real_
}

detect_hammock_subA <- function(filepath) {
  head <- tryCatch(utils::read.csv(filepath, nrows = 10), error = function(e) NULL)
  if (is.null(head)) return(NA_real_)
  for (col in names(head)) {
    low <- tolower(col)
    if (low %in% c('suba', 'sub_a')) {
      vals <- unique(na.omit(head[[col]]))
      if (length(vals) > 0) {
        v <- suppressWarnings(as.numeric(vals[1]))
        if (!is.na(v)) return(v)
      }
    }
  }
  NA_real_
}

detect_hammock_subB <- function(filepath) {
  head <- tryCatch(utils::read.csv(filepath, nrows = 10), error = function(e) NULL)
  if (is.null(head)) return(NA_real_)
  for (col in names(head)) {
    low <- tolower(col)
    if (low %in% c('subb', 'sub_b')) {
      vals <- unique(na.omit(head[[col]]))
      if (length(vals) > 0) {
        v <- suppressWarnings(as.numeric(vals[1]))
        if (!is.na(v)) return(v)
      }
    }
  }
  NA_real_
}

extract_params_from_filename <- function(filename) {
  name <- basename(filename)
  # Mode token
  mode <- {
    m <- regexpr('jacc([A-Da-d])', name, perl = TRUE)
    if (m > 0) toupper(sub('jacc', '', regmatches(name, m))) else NA_character_
  }
  # klen, window, precision
  klen <- {
    m <- regexpr('k(\\d+)', name, perl = TRUE)
    if (m > 0) as.numeric(sub('k', '', regmatches(name, m))) else NA_real_
  }
  window <- {
    m <- regexpr('w(\\d+)', name, perl = TRUE)
    if (m > 0) as.numeric(sub('w', '', regmatches(name, m))) else NA_real_
  }
  precision <- {
    m <- regexpr('p(\\d+)', name, perl = TRUE)
    if (m > 0) as.numeric(sub('p', '', regmatches(name, m))) else NA_real_
  }
  lowered <- tolower(name)
  sketch <- if (grepl('(^|[_.-])(hll|hyperloglog)([_.-]|$)', lowered)) 'hyperloglog'
            else if (grepl('(^|[_.-])(mnmzr|minimizer)([_.-]|$)', lowered)) 'minimizer'
            else if (grepl('(^|[_.-])(kmv|minhash)([_.-]|$)', lowered)) 'kmv' else NA_character_
  # Robust regex extraction for expA, subA, subB (preserve decimals; don't split on '.')
  stem <- tools::file_path_sans_ext(name)
  subA <- NA_real_; subB <- NA_real_; expA <- NA_real_
  # expA like expA2.50 (case-insensitive)
  m_exp <- regexec("(?i)expA([0-9]+(?:\\.[0-9]+)?)", stem, perl = TRUE)
  r_exp <- regmatches(stem, m_exp)
  if (length(r_exp) > 0 && length(r_exp[[1]]) > 1) {
    expA <- suppressWarnings(as.numeric(r_exp[[1]][2]))
  }
  # subA like A0.50
  m_a <- regexec("(^|[^A-Za-z0-9])A([0-9]+(?:\\.[0-9]+)?)", stem, perl = TRUE)
  r_a <- regmatches(stem, m_a)
  if (length(r_a) > 0 && length(r_a[[1]]) > 2) {
    subA <- suppressWarnings(as.numeric(r_a[[1]][3]))
  }
  # subB like B1.25
  m_b <- regexec("(^|[^A-Za-z0-9])B([0-9]+(?:\\.[0-9]+)?)", stem, perl = TRUE)
  r_b <- regmatches(stem, m_b)
  if (length(r_b) > 0 && length(r_b[[1]]) > 2) {
    subB <- suppressWarnings(as.numeric(r_b[[1]][3]))
  }

  
  list(mode = ifelse(is.na(mode), NA, mode),
       klen = ifelse(is.na(klen), NA, klen),
       window = ifelse(is.na(window), NA, window),
       precision = ifelse(is.na(precision), NA, precision),
       sketch = ifelse(is.na(sketch), NA, sketch),
       subA = ifelse(is.na(subA), NA, subA),
       subB = ifelse(is.na(subB), NA, subB),
       expA = ifelse(is.na(expA), NA, expA))
}

