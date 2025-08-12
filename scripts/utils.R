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
  df <- if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(filepath, sep = "\t", header = TRUE, showProgress = FALSE)
  } else {
    utils::read.table(filepath, header = TRUE)
  }
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

load_accession_key <- function(filepath) {
  df <- if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(filepath, sep = "\t", header = TRUE, showProgress = FALSE)
  } else {
    utils::read.table(filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  }
  if (!all(c('File', 'Biosample_term_name') %in% names(df))) stop('Accession key missing required columns: File, Biosample_term_name')
  base <- vapply(df$File, function(x) tools::file_path_sans_ext(basename(x)), character(1))
  setNames(df$Biosample_term_name, base)
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
    list(sim_df, inferred)
  } else {
    sim_sub <- sim_df[common, common, drop = FALSE]
    list(sim_sub, accession_labels)
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
  # Tokens for subA, subB, expA
  stem <- tools::file_path_sans_ext(name)
  tokens <- unlist(strsplit(stem, "[_.-]+"))
  subA <- NA_real_; subB <- NA_real_; expA <- NA_real_
  for (tok in tokens) {
    if (grepl('^expA[0-9.]+$', tok, ignore.case = TRUE)) {
      expA <- suppressWarnings(as.numeric(sub('(?i)expA', '', tok, perl = TRUE)))
    } else if (grepl('^A[0-9.]+$', tok)) {
      subA <- suppressWarnings(as.numeric(sub('A', '', tok)))
    } else if (grepl('^B[0-9.]+$', tok)) {
      subB <- suppressWarnings(as.numeric(sub('B', '', tok)))
    }
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

