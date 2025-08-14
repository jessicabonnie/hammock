library(shiny)
library(ggplot2)
library(shinyFiles)
library(plotly)

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
      tryCatch({ cached <- readRDS(cache_path) }, error = function(e) { err <<- conditionMessage(e) })
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
          if (!is.null(meta$mode) && toupper(meta$mode) == 'C' && (is.null(meta$expA) || is.na(meta$expA) || !is.finite(meta$expA))) {
            # Try to detect expA from file header if not present in name
            meta$expA <- tryCatch(detect_hammock_expA(f), error = function(e) NA_real_)
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
    tryCatch({ tbl <- clustering_analysis(f, acc_key, clusters = clusters_vec, out = out_tmp) }, error = function(e) { err <<- conditionMessage(e) })
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
          obj <- list(similarity = as.matrix(sim_df_for_cache), meta = meta, source_file = basename(f), key_path = acc_key, long_table_raw = tbl)
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
  list(df = df, diag = c(diag, sprintf("rows: %d", nrow(df))))
}

ui <- fluidPage(
  titlePanel("Cluster Analysis Explorer (on hammock outputs)"),
  tabsetPanel(
    tabPanel(
      "Inputs",
      fluidRow(
        column(6,
          shinyDirButton("dir", "Select directory of hammock outputs", "Select"),
          shinyFilesButton("acc", "Select accession key TSV", "Select", multiple = FALSE),
          textInput("dirText", "Or paste outputs directory", value = ""),
          textInput("accText", "Or paste accession key TSV", value = "")
        ),
        column(6,
          h4("Cache options"),
          checkboxInput("useCache", "Use cached RDS if available (same folder as inputs)", value = TRUE)
        )
      ),
      tags$hr(),
      tags$div(
        style = "background-color:#FFF8C6;border:1px solid #E0C97F;padding:10px;border-radius:6px;",
        tags$h4("Diagnostics", style = "margin-top:0;font-weight:bold;color:#8A6D3B;"),
        verbatimTextOutput("diagText")
      )
    ),
    tabPanel(
      "Results",
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
            selectInput("expA", "expA", choices = character(0))
          ),
          conditionalPanel(
            condition = "input.modeSelect == 'D'",
            selectInput("klen", "k-mer size (k)", choices = character(0)),
            selectInput("window", "window size (w)", choices = character(0))
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
      fluidRow(
        column(6,
          selectInput("heatmapFile1", "Heatmap A input file", choices = character(0)),
          plotOutput("heatmap1", height = "600px")
        ),
        column(6,
          selectInput("heatmapFile2", "Heatmap B input file", choices = character(0)),
          plotOutput("heatmap2", height = "600px")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  data_all <- reactiveVal(NULL)
  # Default diagnostics panel content
  output$diagText <- renderText("(no diagnostics yet)")
  # Helper to compute and draw a colored dendrogram for a given file
  draw_dendrogram_for <- function(input_file, accession_key, linkage) {
    fmt <- detect_file_format(input_file)
    sim_df <- if (fmt == 'hammock') parse_hammock_format(input_file) else parse_bedtools_format(input_file)
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
    plot(dnd, main = sprintf('Hierarchical Clustering (%s)%s', linkage, ifelse(tolower(linkage)=='ward', " [ward uses 'complete' on precomputed distances]", "")))
    legend('topright', legend = uniq_labs, col = lab_to_col[uniq_labs], pch = 15, cex = 0.8, bty = 'n')
  }

  # Helper to draw a similarity heatmap for a given file
  draw_heatmap_for <- function(input_file) {
    fmt <- detect_file_format(input_file)
    sim_df <- if (fmt == 'hammock') parse_hammock_format(input_file) else parse_bedtools_format(input_file)
    mat <- as.matrix(sim_df)
    if (nrow(mat) < 2 || ncol(mat) < 2) { plot.new(); title("Not enough samples for heatmap"); return(invisible()) }
    # Convert to long form without extra deps
    df_long <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    names(df_long) <- c("row", "col", "val")
    df_long$row <- as.character(df_long$row)
    df_long$col <- as.character(df_long$col)
    df_long$val <- as.numeric(df_long$val)
    p <- ggplot(df_long, aes(x = col, y = row, fill = val)) +
      geom_raster() +
      scale_fill_viridis_c(option = "viridis", limits = c(0, 1), oob = scales::squish, na.value = "grey95") +
      coord_fixed() +
      labs(x = NULL, y = NULL, fill = "Similarity", title = basename(input_file)) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
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
      output$statusText <- renderText("Select a directory of hammock output files.")
      return()
    }
    # Pre-scan for hammock-like outputs
    all_files <- list.files(path, pattern = "\\.(csv|tsv|txt)$", ignore.case = TRUE, full.names = TRUE)
    is_hammock <- function(fp) {
      first <- tryCatch(readLines(fp, n = 1, warn = FALSE), error = function(e) "")
      if (startsWith(first, "file1,file2")) {
        grepl("jaccard_similarity", first, fixed = TRUE) || grepl("jaccard_similarity_with_ends", first, fixed = TRUE)
      } else {
        # Fallback: sniff header with read.csv
        hdr <- tryCatch(utils::read.csv(fp, nrows = 1, check.names = FALSE), error = function(e) NULL)
        if (is.null(hdr)) return(FALSE)
        any(c("jaccard_similarity", "jaccard_similarity_with_ends") %in% names(hdr))
      }
    }
    hammock_files <- Filter(is_hammock, all_files)
    if (length(hammock_files) == 0) {
      data_all(NULL)
      output$statusText <- renderText("No hammock outputs detected in the selected directory.")
      output$diagText <- renderText(sprintf("Scanned %d files in %s", length(all_files), path))
      return()
    }
    output$statusText <- renderText(sprintf("Running clustering on %d files... (key: %s)", length(hammock_files), ak))
    res <- NULL
    err_top <- NULL
    # Simpler direct call (avoid nested progress scoping issues)
    tryCatch({
      res <- run_clustering_on_dir(
        dir_path = path,
        acc_key = ak,
        files = hammock_files,
        use_cache = isTRUE(input$useCache)
      )
    }, error = function(e) { err_top <<- conditionMessage(e); res <<- NULL })
    if (is.null(res)) {
      data_all(NULL)
      output$statusText <- renderText("Clustering failed to return a result. See diagnostics.")
      output$diagText <- renderText(if (!is.null(err_top)) paste("error:", err_top) else "(no result from run)")
      return()
    }
    if (is.null(res$df)) {
      data_all(NULL)
      output$statusText <- renderText("No valid hammock outputs were processed. Check files and key.")
      output$diagText <- renderText(paste(res$diag, collapse = "\n"))
      return()
    }
    data_all(res$df)
    output$diagText <- renderText(paste(res$diag, collapse = "\n"))
    # Initialize controls
    modes <- sort(unique(na.omit(res$df$mode)))
    if (length(modes)) updateRadioButtons(session, "modeSelect", selected = modes[1])
    linkages <- sort(unique(na.omit(res$df$linkage)))
    if (!length(linkages)) linkages <- c("average","complete","single","ward")
    updateSelectInput(session, "linkageMethod", choices = linkages, selected = if (length(linkages)) linkages[1] else NULL)
    precisions <- sort(unique(na.omit(res$df$precision)))
    updateSelectInput(session, "precision", choices = precisions, selected = if (length(precisions)) precisions[1] else NULL)
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
    updateSelectInput(session, "expA", choices = exp_choices, selected = if (length(exp_vals)) as.character(exp_vals[1]) else NULL)
    k_choices <- sort(unique(na.omit(res$df$klen)))
    w_choices <- sort(unique(na.omit(res$df$window)))
    updateSelectInput(session, "klen", choices = k_choices, selected = if (length(k_choices)) k_choices[1] else NULL)
    updateSelectInput(session, "window", choices = w_choices, selected = if (length(w_choices)) w_choices[1] else NULL)
    srcs <- sort(unique(res$df$source_file))
    updateSelectInput(session, "sourceFile", choices = srcs, selected = if (length(srcs)) srcs[1] else NULL)
    updateSelectInput(session, "heatmapFile1", choices = srcs, selected = if (length(srcs)) srcs[1] else NULL)
    updateSelectInput(session, "heatmapFile2", choices = srcs, selected = if (length(srcs) >= 2) srcs[2] else if (length(srcs)) srcs[1] else NULL)
    # Keep diagnostics minimal (overlap/expA from run_clustering_on_dir)
    output$diagText <- renderText(paste(res$diag, collapse = "\n"))
    output$statusText <- renderText(sprintf("Computed long-tables from %d files (%d rows) using key: %s", length(unique(res$df$source_file)), nrow(res$df), ak))
  })

  filtered_data <- reactive({
    df <- data_all()
    if (is.null(df)) return(NULL)
    req_cols <- c("mode", "linkage", "precision", "n_clusters", "nmi", "silhouette")
    if (!all(req_cols %in% names(df))) return(NULL)
    d <- df
    if (!is.null(input$modeSelect) && nzchar(input$modeSelect)) d <- d[d$mode == input$modeSelect, , drop = FALSE]
    if (!is.null(input$linkageMethod) && nzchar(input$linkageMethod)) d <- d[d$linkage == input$linkageMethod, , drop = FALSE]
    if (!is.null(input$precision) && nzchar(input$precision)) {
      suppressWarnings(pv <- as.numeric(input$precision)); if (!is.na(pv)) d <- d[!is.na(d$precision) & d$precision == pv, , drop = FALSE]
    }
    if (identical(input$modeSelect, "C") && "expA" %in% names(d) && !is.null(input$expA) && nzchar(input$expA)) {
      suppressWarnings(v <- as.numeric(input$expA))
      if (!is.na(v)) d <- d[!is.na(d$expA) & abs(d$expA - v) < 1e-9, , drop = FALSE]
    }
    if (identical(input$modeSelect, "D")) {
      if ("klen" %in% names(d) && !is.null(input$klen) && nzchar(input$klen)) { suppressWarnings(kv <- as.numeric(input$klen)); if (!is.na(kv)) d <- d[!is.na(d$klen) & d$klen == kv, , drop = FALSE] }
      if ("window" %in% names(d) && !is.null(input$window) && nzchar(input$window)) { suppressWarnings(wv <- as.numeric(input$window)); if (!is.na(wv)) d <- d[!is.na(d$window) & d$window == wv, , drop = FALSE] }
    }
    if (nrow(d) == 0) return(NULL)
    d
  })

  output$nmiPlot <- renderPlotly({
    d <- filtered_data()
    validate(need(!is.null(d), "No data to plot with the current filters."))
    d <- d[!is.na(d$n_clusters) & !is.na(d$nmi), , drop = FALSE]
    validate(need(nrow(d) > 0, "No rows with NMI and cluster numbers available."))
    limits <- silhouette_limits()
    p <- ggplot(d, aes(x = n_clusters, y = nmi, color = silhouette, text = paste0(
        "clusters: ", n_clusters,
        "<br>NMI: ", sprintf("%.4f", nmi),
        "<br>silhouette: ", sprintf("%.4f", silhouette),
        "<br>linkage: ", linkage
      ))) +
      geom_point(size = 2.5, alpha = 0.9) +
      scale_color_viridis_c(option = "viridis", na.value = "grey80", limits = limits, oob = scales::squish) +
      labs(x = "Number of clusters", y = "NMI", color = "Silhouette",
           title = sprintf("NMI vs #Clusters (mode %s, linkage: %s)", input$modeSelect %||% "", input$linkageMethod %||% "")) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
    ggplotly(p, tooltip = "text") %>% layout(hovermode = "closest")
  })

  output$dendPlot <- renderPlot({
    df <- data_all()
    validate(need(!is.null(df), "Load data to view dendrogram."))
    sf <- input$sourceFile
    lk <- input$dendLinkage %||% "average"
    validate(need(!is.null(sf) && nzchar(sf), "Select an input file for dendrogram."))
    # find matching original file path by name from current directory
    # Try both pasted dir and picker dir
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
    draw_dendrogram_for(candidate, effective_acc(), lk)
  })

  output$heatmap1 <- renderPlot({
    df <- data_all()
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
    draw_heatmap_for(candidate)
  })

  output$heatmap2 <- renderPlot({
    df <- data_all()
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
    draw_heatmap_for(candidate)
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
