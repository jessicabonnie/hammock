library(shiny)
library(ggplot2)
library(shinyFiles)

# Run Python or R clustering on a directory of hammock outputs to produce long-tables in-memory
run_clustering_on_dir <- function(dir_path, acc_key, files = NULL) {
  diag <- character(0)
  if (is.null(files)) files <- list.files(dir_path, pattern = "\\.(csv|tsv|txt)$", ignore.case = TRUE, full.names = TRUE)
  diag <- c(diag, sprintf("dir: %s", dir_path))
  diag <- c(diag, sprintf("accession_key: %s", acc_key))
  diag <- c(diag, sprintf("files_scanned: %d", length(files)))
  if (length(files) == 0) return(list(df = NULL, diag = c(diag, "No files found in directory.")))
  if (is.null(acc_key) || !nzchar(acc_key) || !file.exists(acc_key)) return(list(df = NULL, diag = c(diag, "Accession key missing or not found.")))
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
  for (f in files) {
    # Decide if hammock or bedtools based on header; let the R function decide
    out_tmp <- tempfile(fileext = ".tsv")
    tbl <- NULL
    err <- NULL
    tryCatch({ tbl <- clustering_analysis(f, acc_key, out = out_tmp) }, error = function(e) { err <<- conditionMessage(e) })
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
    }
    else {
      diag <- c(diag, sprintf("fail: %s :: %s", basename(f), ifelse(is.null(err), "unknown error", err)))
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
  # Top row: three columns of inputs
  fluidRow(
    column(4,
      shinyDirButton("dir", "Select directory of hammock outputs", "Select"),
      shinyFilesButton("acc", "Select accession key TSV", "Select", multiple = FALSE),
      textInput("dirText", "Or paste outputs directory", value = ""),
      textInput("accText", "Or paste accession key TSV", value = "")
    ),
    column(4,
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
    column(4,
      h4("Dendrogram controls"),
      selectInput("sourceFile", "Input file", choices = character(0)),
      selectInput("dendLinkage", "Dendrogram linkage", choices = c("average","complete","single","ward"))
    )
  ),
  tags$hr(),
  verbatimTextOutput("statusText"),
  # Second row: two plots side-by-side
  fluidRow(
    column(6, plotOutput("nmiPlot", height = "520px")),
    column(6, plotOutput("dendPlot", height = "520px"))
  ),
  tags$hr(),
  tags$div(
    style = "background-color:#FFF8C6;border:1px solid #E0C97F;padding:10px;border-radius:6px;",
    tags$h4("Diagnostics", style = "margin-top:0;font-weight:bold;color:#8A6D3B;"),
    verbatimTextOutput("diagText")
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
  shinyFileChoose(input, "acc", roots = volumes, session = session)
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
    res <- run_clustering_on_dir(path, ak, files = hammock_files)
    if (is.null(res$df)) {
      data_all(NULL)
      output$statusText <- renderText("No valid hammock outputs were processed. Check files and key.")
      output$diagText <- renderText(paste(res$diag, collapse = "\n"))
      return()
    }
    data_all(res$df)
    output$diagText <- renderText(paste(res$diag, collapse = "\n"))
    # Initialize controls
    modes <- sort(unique(na.omit(res$df$mode))); if (length(modes)) updateRadioButtons(session, "modeSelect", selected = modes[1])
    linkages <- sort(unique(na.omit(res$df$linkage))); updateSelectInput(session, "linkageMethod", choices = linkages, selected = if (length(linkages)) linkages[1] else character(0))
    precisions <- sort(unique(na.omit(res$df$precision))); updateSelectInput(session, "precision", choices = precisions, selected = if (length(precisions)) precisions[1] else character(0))
    exp_vals <- sort(unique(na.omit(res$df$expA)))
    updateSelectInput(session, "expA", choices = exp_vals, selected = if (length(exp_vals)) exp_vals[1] else character(0))
    k_choices <- sort(unique(na.omit(res$df$klen)))
    w_choices <- sort(unique(na.omit(res$df$window)))
    updateSelectInput(session, "klen", choices = k_choices, selected = if (length(k_choices)) k_choices[1] else character(0))
    updateSelectInput(session, "window", choices = w_choices, selected = if (length(w_choices)) w_choices[1] else character(0))
    updateSelectInput(session, "sourceFile", choices = sort(unique(res$df$source_file)), selected = sort(unique(res$df$source_file))[1])
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
    if (identical(input$modeSelect, "C") && "expA" %in% names(d) && length(input$expA) == 2) {
      rng <- sort(as.numeric(input$expA)); d <- d[!is.na(d$expA) & d$expA >= rng[1] & d$expA <= rng[2], , drop = FALSE]
    }
    if (identical(input$modeSelect, "D")) {
      if ("klen" %in% names(d) && !is.null(input$klen) && nzchar(input$klen)) { suppressWarnings(kv <- as.numeric(input$klen)); if (!is.na(kv)) d <- d[!is.na(d$klen) & d$klen == kv, , drop = FALSE] }
      if ("window" %in% names(d) && !is.null(input$window) && nzchar(input$window)) { suppressWarnings(wv <- as.numeric(input$window)); if (!is.na(wv)) d <- d[!is.na(d$window) & d$window == wv, , drop = FALSE] }
    }
    if (nrow(d) == 0) return(NULL)
    d
  })

  output$nmiPlot <- renderPlot({
    d <- filtered_data()
    validate(need(!is.null(d), "No data to plot with the current filters."))
    d <- d[!is.na(d$n_clusters) & !is.na(d$nmi), , drop = FALSE]
    validate(need(nrow(d) > 0, "No rows with NMI and cluster numbers available."))
    limits <- silhouette_limits()
    ggplot(d, aes(x = n_clusters, y = nmi, color = silhouette)) +
      geom_point(size = 2.5, alpha = 0.9) +
      scale_color_viridis_c(option = "viridis", na.value = "grey80", limits = limits, oob = scales::squish) +
      labs(x = "Number of clusters", y = "NMI", color = "Silhouette",
           title = sprintf("NMI vs #Clusters (mode %s, linkage: %s)", input$modeSelect %||% "", input$linkageMethod %||% "")) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
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
}

`%||%` <- function(a, b) if (!is.null(a) && nchar(a) > 0) a else b

shinyApp(ui, server)
