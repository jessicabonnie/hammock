#!/usr/bin/env Rscript
# Analyze maurano_dhs_validation sweep outputs vs two references.
#
# Inputs:
#   data/maurano_bedtools_ref.tsv        -- bedtools pairwise Jaccard
#   data/maurano_filenames_key.tsv       -- tissue metadata
#   results/raw_abc/*.csv                -- A/B/C configs
#   results/raw_d/*.csv                  -- Mode D (k, w, p) configs
#
# Outputs (figures/ and results/):
#   figures/abc_pearson_vs_bedtools.png       -- A/B/C accuracy bar
#   figures/mode_d_pearson_heatmap.png        -- Pearson r vs bedtools
#   figures/mode_d_mae_heatmap.png            -- MAE vs bedtools
#   figures/mode_d_pearson_spearman_heatmap.png -- Pearson | Spearman, side-by-side
#   figures/mode_d_vs_bedtools_vs_modeB.png   -- vs bedtools | vs Mode B, side-by-side
#   figures/mode_d_clustering_ari.png         -- ARI(tissue labels) heatmap
#   figures/mode_d_clustering_nmi.png         -- NMI(tissue labels) heatmap
#   figures/mode_d_best_scatter.png           -- best D vs bedtools | vs Mode B
#   figures/bedtools_dendrogram.png           -- bedtools reference clustering
#   figures/mode_d_best_dendrogram.png        -- best Mode D clustering
#   results/abc_summary.csv                   -- per-config Mode A/B/C accuracy
#   results/mode_d_summary.csv                -- per-config Mode D accuracy
#       columns: k, w, precision, column, reference, n, pearson, spearman,
#                mae, rmse, frob, max_err, ari, nmi, path
#       (ari, nmi computed once per (config, column); replicated across references)
#
# Usage:
#   ml r/4.3.0
#   Rscript analyze.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(Cairo)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) >= 1) {
  dirname(normalizePath(sub("--file=", "", file_arg[1])))
} else {
  getwd()
}

data_dir    <- file.path(script_dir, "data")
results_dir <- file.path(script_dir, "results")
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(path, plot, width = 7, height = 4.5, dpi = 150) {
  CairoPNG(filename = path, width = width, height = height,
           units = "in", res = dpi)
  on.exit(dev.off())
  print(plot)
}

# Strip recognised genomic extensions (.bed, .fa, .fasta, .fna). Mode D
# CSVs key on the FASTA basename while the bedtools reference keys on
# the BED basename; we normalise both sides to a common stem before joining.
strip_ext <- function(x) {
  sub("\\.(bed|fa|fasta|fna)$", "", x, ignore.case = TRUE)
}

# ── reg_eq_similarity / jaccard_similarity fallback (rename in progress) ────
# jaccard_similarity was renamed reg_eq_similarity
# (docs/seed-jaccard-reg-eq-rename.md Step 1) -- it's the register-equality
# statistic, not set Jaccard. This script both reads hammock output and
# re-emits the value into several downstream CSVs, so every read site below
# prefers the new name and falls back to the old one per-file for archived
# CSVs Step 1 didn't rewrite. jaccard_similarity_ie is a different, unrelated
# column and is never touched by this fallback.
REG_EQ_COL        <- "reg_eq_similarity"
REG_EQ_LEGACY_COL <- "jaccard_similarity"
REG_EQ_VALUES     <- c(REG_EQ_COL, REG_EQ_LEGACY_COL)

# Log the very first fallback per run, not per row/file -- a bare presence
# check alone can't distinguish "this is a legitimately old archived file"
# from "a fresh run silently regressed and is still writing the old name"
# (e.g. a partial merge that reverted runner.py's rename but not
# hammock_cli.cpp's).
.reg_eq_fallback_logged <- FALSE
log_reg_eq_fallback <- function(context) {
  if (!.reg_eq_fallback_logged) {
    message("Note: '", REG_EQ_COL, "' not found, falling back to legacy '",
            REG_EQ_LEGACY_COL, "' (first seen at: ", context, "). Expected ",
            "for archived pre-rename files; unexpected for a freshly ",
            "regenerated one.")
    .reg_eq_fallback_logged <<- TRUE
  }
}

# ── matrix + clustering helpers (used both for figures and for ARI/NMI) ─────
to_matrix <- function(df) {
  # df has stem1/stem2/j_hat (extension stripped). Symmetrize, fill diag=1.
  files <- sort(unique(c(df$stem1, df$stem2)))
  m <- matrix(NA_real_, nrow = length(files), ncol = length(files),
              dimnames = list(files, files))
  for (i in seq_len(nrow(df))) {
    m[df$stem1[i], df$stem2[i]] <- df$j_hat[i]
  }
  for (i in seq_along(files)) for (j in seq_along(files)) {
    if (is.na(m[i, j]) && !is.na(m[j, i])) m[i, j] <- m[j, i]
  }
  diag(m) <- 1
  m
}

# ── ARI / NMI implementations (no external package deps) ───────────────────
adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  c_combs <- function(x) sum(choose(x, 2))
  idx     <- c_combs(as.vector(tab))
  exp_idx <- c_combs(rowSums(tab)) * c_combs(colSums(tab)) / choose(n, 2)
  max_idx <- (c_combs(rowSums(tab)) + c_combs(colSums(tab))) / 2
  if (isTRUE(all.equal(max_idx, exp_idx))) return(0)
  (idx - exp_idx) / (max_idx - exp_idx)
}

nmi <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  pi_ <- rowSums(tab) / n
  pj_ <- colSums(tab) / n
  pij <- tab / n
  H_a <- -sum(pi_[pi_ > 0] * log(pi_[pi_ > 0]))
  H_b <- -sum(pj_[pj_ > 0] * log(pj_[pj_ > 0]))
  if (H_a == 0 || H_b == 0) return(0)
  outer_p <- outer(pi_, pj_)
  mask    <- (pij > 0) & (outer_p > 0)
  I       <- sum(pij[mask] * log(pij[mask] / outer_p[mask]))
  2 * I / (H_a + H_b)
}

# ---------- bedtools reference ----------
ref_path <- file.path(data_dir, "maurano_bedtools_ref.tsv")
if (!file.exists(ref_path)) {
  stop("Missing ", ref_path, " -- run ./prepare_data.sh first.")
}
ref <- read_tsv(ref_path, show_col_types = FALSE) %>%
  transmute(stem1 = strip_ext(file1),
            stem2 = strip_ext(file2),
            j_truth = as.numeric(jaccard))

# ---------- tissue metadata ----------
key_path <- file.path(data_dir, "maurano_filenames_key.tsv")
key <- if (file.exists(key_path)) {
  read_tsv(key_path, show_col_types = FALSE) %>%
    transmute(file = File, tissue = Biosample_term_name)
} else {
  tibble(file = character(), tissue = character())
}

# stem -> tissue label, used for ARI/NMI.
tissue_for <- function(stems) {
  if (nrow(key) == 0) return(setNames(rep("?", length(stems)), stems))
  setNames(key$tissue[match(stems, strip_ext(key$file))], stems)
}
tissue_label_vec <- tissue_for(unique(c(ref$stem1, ref$stem2)))
n_tissues <- length(unique(tissue_label_vec))

# ---------- second reference: Mode B at p=21 ----------
# Mode B is r=0.998 vs bedtools, so it's effectively the interval-mode
# reading of the same biology. Asking whether Mode D agrees with Mode B
# is a sharper question than asking whether it agrees with bedtools.
modeB_path <- file.path(results_dir, "raw_abc", "hammock_hll_p21_jaccB_full.csv")
modeB_ref <- if (file.exists(modeB_path)) {
  modeB_raw <- read_csv(modeB_path, show_col_types = FALSE)
  modeB_col <- if (REG_EQ_COL %in% names(modeB_raw)) {
    REG_EQ_COL
  } else {
    log_reg_eq_fallback(modeB_path)
    REG_EQ_LEGACY_COL
  }
  modeB_raw %>%
    transmute(stem1 = strip_ext(file1), stem2 = strip_ext(file2),
              j_truth = as.numeric(.data[[modeB_col]]))
} else {
  warning("Mode B p=21 reference not found at ", modeB_path)
  NULL
}

# ---------- helpers ----------
read_hammock_csv <- function(path, jcol = NULL) {
  # All hammock CSVs have file1, file2 + a similarity column. Mode D used to
  # emit a second `jaccard_similarity_with_ends`; that column was removed in
  # hammock v0.6.0 (CLAUDE.md divergence #8) after it was shown to track
  # bedtools less cleanly on this very corpus. Archived CSVs still have it.
  df <- read_csv(path, show_col_types = FALSE)
  if (is.null(jcol)) jcol <- REG_EQ_COL
  # reg_eq_similarity was renamed from jaccard_similarity
  # (docs/seed-jaccard-reg-eq-rename.md Step 1); CSVs written before that
  # rename still carry the old header. Prefer the new name, fall back to the
  # old one per-file, logging the fallback once per run.
  if (identical(jcol, REG_EQ_COL) && !(jcol %in% names(df)) &&
      REG_EQ_LEGACY_COL %in% names(df)) {
    log_reg_eq_fallback(path)
    jcol <- REG_EQ_LEGACY_COL
  }
  # The CSVs under raw_d/ predate the jaccard_similarity_ie column (v0.4.0) but
  # already carry the containments it is derived from, so recover it here rather
  # than resweeping 235 configurations. Mirrors
  # runner._jaccard_ie_from_containments, clamp included: containments can
  # exceed 1 by an ulp, and clamping is what makes denom >= 1.
  if (identical(jcol, "jaccard_similarity_ie") && !(jcol %in% names(df)) &&
      all(c("containment_AB", "containment_BA") %in% names(df))) {
    cab <- pmin(as.numeric(df$containment_AB), 1)
    cba <- pmin(as.numeric(df$containment_BA), 1)
    df[[jcol]] <- ifelse(cab <= 0 | cba <= 0, 0, 1 / (1 / cab + 1 / cba - 1))
  }
  if (!(jcol %in% names(df))) {
    stop("column ", jcol, " not in ", path)
  }
  df %>% transmute(stem1 = strip_ext(file1),
                   stem2 = strip_ext(file2),
                   j_hat = .data[[jcol]])
}

compute_metrics_vs_ref <- function(j_hat_df, ref_df, ref_label) {
  joined <- j_hat_df %>% inner_join(ref_df, by = c("stem1", "stem2"))
  if (nrow(joined) == 0) {
    return(tibble(reference = ref_label, n = 0L,
                  pearson = NA_real_, spearman = NA_real_,
                  mae = NA_real_, rmse = NA_real_,
                  frob = NA_real_, max_err = NA_real_))
  }
  diffs <- joined$j_hat - joined$j_truth
  sd_ok <- sd(joined$j_hat) > 0 && sd(joined$j_truth) > 0
  tibble(
    reference = ref_label,
    n        = nrow(joined),
    pearson  = if (sd_ok) cor(joined$j_hat, joined$j_truth, method = "pearson")  else NA_real_,
    spearman = if (sd_ok) cor(joined$j_hat, joined$j_truth, method = "spearman") else NA_real_,
    mae      = mean(abs(diffs)),
    rmse     = sqrt(mean(diffs^2)),
    frob     = sqrt(sum(diffs^2)),
    max_err  = max(abs(diffs))
  )
}

# Cluster-agreement vs true tissue labels. Reference-independent (the
# clustering is built from j_hat alone), so the same (ari, nmi) gets
# joined onto every reference row for this (config, column).
compute_cluster_metrics <- function(j_hat_df) {
  mat <- tryCatch(to_matrix(j_hat_df), error = function(e) NULL)
  if (is.null(mat) || any(is.na(mat))) return(tibble(ari = NA_real_, nmi = NA_real_))
  # Clamp to [0,1] (HLL noise can drift slightly outside) without losing
  # matrix dimensions — pmin/pmax flatten the matrix to a vector.
  mat[mat < 0] <- 0
  mat[mat > 1] <- 1
  d  <- as.dist(1 - mat)
  hc <- hclust(d, method = "average")
  pred <- cutree(hc, k = n_tissues)
  true <- tissue_label_vec[rownames(mat)]
  tibble(ari = adjusted_rand(true, pred),
         nmi = nmi(true, pred))
}

# Parse hammock filenames produced by python/hammock/outprefix.py.
# Trailing "_full" is the metrics-shape tag (always present now that these
# runs pass --metrics; python/hammock/outprefix.py,
# docs/seed-metrics-column-restructure.md). Examples (stem = basename minus
# .csv):
#   hammock_hll_p21_jaccA_full.csv
#   hammock_hll_p21_jaccB_B0.10_full.csv
#   hammock_hll_p21_jaccC_expA1.50_full.csv
#   hammock_hll_p21_jaccC_B0.25_full.csv
#   hammock_mnmzr_p24_jaccD_k10_w20_full.csv
parse_abc_name <- function(stem) {
  m <- regmatches(stem, regexec(
    "^hammock_hll_p(\\d+)_jacc([ABC])(?:_expA(\\d+\\.\\d+))?(?:_B(\\d+\\.\\d+))?_full$",
    stem, perl = TRUE))[[1]]
  if (length(m) == 0) return(NULL)
  tibble(precision = as.integer(m[2]),
         mode = m[3],
         expA = if (nzchar(m[4])) as.numeric(m[4]) else 0,
         subB = if (nzchar(m[5])) as.numeric(m[5]) else 1)
}

parse_d_name <- function(stem) {
  m <- regmatches(stem, regexec(
    "^hammock_mnmzr_p(\\d+)_jaccD_k(\\d+)_w(\\d+)_full$", stem, perl = TRUE))[[1]]
  if (length(m) == 0) return(NULL)
  tibble(precision = as.integer(m[2]), k = as.integer(m[3]),
         w = as.integer(m[4]))
}

scan_dir <- function(dir, parser, jcols = REG_EQ_COL,
                     refs = list(bedtools = NULL),
                     do_clustering = FALSE) {
  # `refs` is a named list of reference dataframes (each with stem1, stem2,
  # j_truth). Output gets one row per (file × jcol × reference). If
  # do_clustering, also attach (ari, nmi) replicated across reference rows.
  files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  out <- list()
  for (f in files) {
    stem <- tools::file_path_sans_ext(basename(f))
    meta <- parser(stem)
    if (is.null(meta)) next
    cat("  reading", basename(f), "\n")
    for (jc in jcols) {
      j_hat <- read_hammock_csv(f, jcol = jc)
      cm <- if (do_clustering) compute_cluster_metrics(j_hat)
            else tibble(ari = NA_real_, nmi = NA_real_)
      for (ref_name in names(refs)) {
        ref_df <- refs[[ref_name]]
        if (is.null(ref_df)) next
        m <- compute_metrics_vs_ref(j_hat, ref_df, ref_name)
        out[[length(out) + 1L]] <- bind_cols(meta, tibble(column = jc),
                                             m, cm, tibble(path = f))
      }
    }
  }
  if (length(out) == 0) return(tibble())
  bind_rows(out)
}

# ---------- scan ABC ----------
abc_dir <- file.path(results_dir, "raw_abc")
cat("Scanning", abc_dir, "...\n")
abc <- scan_dir(abc_dir, parse_abc_name,
                refs = list(bedtools = ref),
                do_clustering = FALSE)

if (nrow(abc) > 0) {
  abc_out <- abc %>% arrange(mode, precision, expA, subB)
  write_csv(abc_out, file.path(results_dir, "abc_summary.csv"))
  cat("Wrote results/abc_summary.csv (", nrow(abc_out), "rows)\n")

  # Bar plot: Pearson r per A/B/C config. Identify configs by label.
  abc_plot <- abc_out %>%
    mutate(label = case_when(
      mode == "A" ~ sprintf("A p=%d", precision),
      mode == "B" ~ sprintf("B p=%d", precision),
      mode == "C" & expA > 0 ~ sprintf("C p=%d expA=%.2f", precision, expA),
      mode == "C" & subB < 1 ~ sprintf("C p=%d subB=%.2f", precision, subB),
      TRUE ~ sprintf("C p=%d", precision)
    ))
  p1 <- ggplot(abc_plot, aes(x = reorder(label, pearson), y = pearson, fill = mode)) +
    geom_col() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    coord_flip() +
    labs(x = NULL, y = "Pearson r vs bedtools Jaccard",
         title = "Modes A/B/C accuracy on Maurano DHS",
         fill = "Mode") +
    theme_minimal(base_size = 11)
  save_png(file.path(figures_dir, "abc_pearson_vs_bedtools.png"),
           p1, width = 8, height = 0.25 * nrow(abc_plot) + 2)

  # Mode C subB and expA interpolation plots (|C-A| vs |C-B|) live in
  # mode_c_interpolation.R — that script reads the same raw CSVs and
  # writes mode_c_subB_interpolation_agg.png / mode_c_expA_interpolation_agg.png.
} else {
  cat("No A/B/C outputs found in", abc_dir, "-- skipping.\n")
}

# ---------- scan Mode D (all current metrics × references + clustering) ────
d_dir <- file.path(results_dir, "raw_d")
cat("Scanning", d_dir, "...\n")
d_refs <- list(bedtools = ref)
if (!is.null(modeB_ref)) d_refs$mode_B <- modeB_ref
# Keep this canonical order aligned with runner._write_mode_d_csv. Discover
# availability from the files instead of maintaining a separate "Jaccard-only"
# list: current Mode D output carries this entire minimizer block. The removed
# `_with_ends` family is intentionally absent even when archived CSVs contain it.
d_metric_order <- c(
  REG_EQ_COL, "jaccard_similarity_ie",
  "containment_AB", "containment_BA",
  "cosketch_geom", "cosketch_arith", "cosketch_max"
)
d_files <- list.files(d_dir, pattern = "\\.csv$", full.names = TRUE)
d_headers <- lapply(d_files, function(p) {
  names(read_csv(p, n_max = 0, show_col_types = FALSE))
})
# reg_eq_similarity was renamed from jaccard_similarity
# (docs/seed-jaccard-reg-eq-rename.md Step 1). Accept either name per-file,
# mapping both to the one canonical REG_EQ_COL label below -- the old
# `metric %in% header` check alone required the literal old name in *every*
# scanned file's header, which would silently drop the metric from discovery
# the moment raw_d/ contains a mix of pre- and post-rename CSVs.
metric_available <- function(metric, header) {
  metric %in% header ||
    (identical(metric, REG_EQ_COL) && REG_EQ_LEGACY_COL %in% header) ||
    (identical(metric, "jaccard_similarity_ie") &&
       all(c("containment_AB", "containment_BA") %in% header))
}
if (any(vapply(d_headers, function(h)
    !(REG_EQ_COL %in% h) && (REG_EQ_LEGACY_COL %in% h), logical(1)))) {
  log_reg_eq_fallback("Mode D discovery gate (raw_d/*.csv headers)")
}
d_metrics <- Filter(
  function(metric) {
    length(d_headers) > 0 &&
      all(vapply(d_headers, function(header) metric_available(metric, header),
                 logical(1)))
  },
  d_metric_order
)
cat("  Mode D metrics:", paste(d_metrics, collapse = ", "), "\n")
d <- scan_dir(d_dir, parse_d_name,
              jcols = d_metrics,
              refs = d_refs,
              do_clustering = TRUE)

best_d_path <- NULL
best_d_col  <- NULL
if (nrow(d) > 0) {
  d_out <- d %>% arrange(column, precision, k, w)
  write_csv(d_out, file.path(results_dir, "mode_d_summary.csv"))
  cat("Wrote results/mode_d_summary.csv (", nrow(d_out), "rows)\n")

  # Both the canonical REG_EQ_COL ("reg_eq_similarity") and the pre-rename
  # REG_EQ_LEGACY_COL ("jaccard_similarity") map to the same output label --
  # d_out$column can carry either depending on which name a given raw_d/
  # file's header used (see read_hammock_csv()'s fallback above) -- so both
  # keys must resolve identically for downstream consumers (e.g. the
  # pivot_wider(names_from = col_short, ...) below, which bakes this exact
  # string into hardcoded column names like `ari_register_equality`).
  short_col <- function(x) dplyr::recode(
    x,
    reg_eq_similarity = "register_equality",
    jaccard_similarity = "register_equality",
    jaccard_similarity_ie = "ie",
    .default = x
  )
  d_out$col_short <- short_col(d_out$column)

  # Keep the existing published figures on their declared primary estimator.
  # The complete discovered metric block is available in mode_d_summary.csv;
  # adding summary metrics must not silently redraw every historical heatmap.
  d_plot     <- d_out %>% filter(column %in% REG_EQ_VALUES)
  d_bedtools <- d_plot %>% filter(reference == "bedtools")
  d_modeB    <- d_plot %>% filter(reference == "mode_B")
  n_p <- length(unique(d_plot$precision))

  # ── does the estimator choice change anything? (see
  #    docs/estimator-analysis-findings.md section 9) ────────────────────────
  ie_cmp <- d_out %>%
    filter(reference == "bedtools",
           col_short %in% c("register_equality", "ie")) %>%
    select(precision, k, w, col_short, pearson, mae, ari) %>%
    tidyr::pivot_wider(names_from = col_short,
                       values_from = c(pearson, mae, ari))
  if (nrow(ie_cmp) > 0 && "ari_ie" %in% names(ie_cmp)) {
    ie_arm_summary <- ie_cmp %>%
      group_by(precision) %>%
      summarise(n = n(),
                ari_differs = sum(!is.na(ari_register_equality) &
                                    !is.na(ari_ie) &
                                    abs(ari_register_equality - ari_ie) > 1e-9),
                max_abs_ari_delta = max(
                  abs(ari_register_equality - ari_ie), na.rm = TRUE),
                median_mae_register_equality =
                  median(mae_register_equality, na.rm = TRUE),
                median_mae_ie = median(mae_ie, na.rm = TRUE),
                .groups = "drop")
    cat("\nEstimator arm comparison (Mode D, reference = bedtools):\n")
    print(ie_arm_summary, n = 50)
    # Persisted so the supplement can cite a file rather than a console dump;
    # this table is the evidence for "the two columns agree where the figures
    # make their claims but not across the whole sweep".
    write_csv(ie_arm_summary,
              file.path(results_dir, "estimator_arm_comparison.csv"))
    cat("Wrote results/estimator_arm_comparison.csv\n")

    # Per-cell, not just the per-precision aggregate above. The aggregate can
    # only answer "how many cells differ"; this answers "which ones, and by how
    # much", which is the form the question actually gets asked in. Note the
    # answer is NOT "ARI is identical at every (k, w)" --- it is identical at
    # the cells the published figures use, and differs at 48 of 235 overall.
    ie_by_config <- ie_cmp %>%
      transmute(precision, k, w,
                ari_register_equality, ari_ie,
                ari_delta = ari_ie - ari_register_equality,
                ari_identical = !is.na(ari_register_equality) & !is.na(ari_ie) &
                  abs(ari_ie - ari_register_equality) <= 1e-9,
                mae_register_equality, mae_ie,
                pearson_register_equality, pearson_ie) %>%
      arrange(desc(abs(ari_delta)), precision, k, w)
    write_csv(ie_by_config,
              file.path(results_dir, "estimator_ari_by_config.csv"))
    cat(sprintf(
      "Wrote results/estimator_ari_by_config.csv (%d cells, %d with identical ARI)\n",
      nrow(ie_by_config), sum(ie_by_config$ari_identical, na.rm = TRUE)))
  }

  # ── existing Pearson + MAE vs bedtools (filtered now) ────────────────────
  p4 <- ggplot(d_bedtools, aes(x = factor(w), y = factor(k), fill = pearson)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Pearson r", limits = c(0, 1.0),
                         oob = squish) +
    facet_grid(col_short ~ precision, labeller = label_both) +
    labs(x = "w", y = "k", title = "Mode D Pearson r vs bedtools") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_png(file.path(figures_dir, "mode_d_pearson_heatmap.png"), p4,
           width = max(10, 1.6 * n_p), height = 6)

  p5 <- ggplot(d_bedtools, aes(x = factor(w), y = factor(k), fill = mae)) +
    geom_tile() +
    scale_fill_viridis_c(name = "MAE", option = "magma", direction = -1,
                         limits = c(0, 0.8), oob = squish) +
    facet_grid(col_short ~ precision, labeller = label_both) +
    labs(x = "w", y = "k", title = "Mode D MAE vs bedtools") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_png(file.path(figures_dir, "mode_d_mae_heatmap.png"), p5,
           width = max(10, 1.6 * n_p), height = 6)

  # ── Pearson vs Spearman per config (scatter, one panel) ─────────────────
  ps <- d_bedtools %>%
    filter(column %in% REG_EQ_VALUES, !is.na(pearson), !is.na(spearman))
  p_ps <- ggplot(ps, aes(x = pearson, y = spearman,
                         colour = factor(precision))) +
    geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed") +
    geom_point(alpha = 0.8, size = 2) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Pearson r vs bedtools", y = "Spearman rho vs bedtools",
         colour = "precision",
         title = "Mode D: Pearson vs Spearman per config (no_ends column)",
         subtitle = "Points on y=x ⇒ relationship is monotone; gap below y=x ⇒ residual rank shuffling") +
    theme_minimal(base_size = 11)
  save_png(file.path(figures_dir, "mode_d_pearson_vs_spearman_scatter.png"),
           p_ps, width = 7, height = 6)

  # ── r vs bedtools  vs  r vs Mode B per config (scatter, one panel) ──────
  if (!is.null(modeB_ref)) {
    refs_wide <- d_out %>%
      filter(column %in% REG_EQ_VALUES, !is.na(pearson)) %>%
      select(k, w, precision, reference, pearson) %>%
      pivot_wider(names_from = reference, values_from = pearson) %>%
      filter(!is.na(bedtools), !is.na(mode_B))
    p_refs <- ggplot(refs_wide,
                     aes(x = bedtools, y = mode_B,
                         colour = factor(precision))) +
      geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed") +
      geom_point(alpha = 0.8, size = 2) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = "Pearson r vs bedtools", y = "Pearson r vs interval mode (p=21)",
           colour = "precision",
           title = "Sequence mode: agreement with bedtools vs agreement with interval mode (per config)",
           subtitle = "Points on y=x ⇒ interval mode and bedtools are interchangeable as references") +
      theme_minimal(base_size = 11)
    save_png(file.path(figures_dir, "mode_d_bedtools_vs_modeB_scatter.png"),
             p_refs, width = 7, height = 6)
  }

  # ── Metric tradeoff: Pearson r vs ARI per config (the headline tradeoff) ─
  tr <- d_bedtools %>%
    filter(column %in% REG_EQ_VALUES, !is.na(pearson), !is.na(ari)) %>%
    mutate(label = sprintf("k=%d w=%d p=%d", k, w, precision))
  best_pear <- tr %>% arrange(desc(pearson)) %>% slice(1)
  best_ari  <- tr %>% arrange(desc(ari))     %>% slice(1)
  best_mae  <- tr %>% arrange(mae)           %>% slice(1)
  best_pts <- bind_rows(best_pear, best_ari, best_mae)
  p_tr <- ggplot(tr, aes(x = pearson, y = ari, colour = factor(precision))) +
    geom_point(alpha = 0.7, size = 2) +
    geom_point(data = best_pts,
               size = 4, shape = 21, fill = NA, colour = "black", stroke = 1.2) +
    geom_text(data = best_pts, aes(label = label),
              nudge_y = 0.03, size = 3, inherit.aes = TRUE, show.legend = FALSE) +
    labs(x = "Pearson r vs bedtools (numerical agreement)",
         y = "ARI vs tissue labels (clustering recovery)",
         colour = "precision",
         title = "Mode D: numerical agreement vs biological-clustering recovery",
         subtitle = "Circled = best of each metric. Best-numerical and best-clustering are different configs.") +
    theme_minimal(base_size = 11)
  save_png(file.path(figures_dir, "mode_d_metric_tradeoff.png"),
           p_tr, width = 8, height = 6)

  # ── ARI / NMI heatmaps (reference-independent) ──────────────────────────
  d_cluster <- d_plot %>% filter(reference == "bedtools")  # drop duplicate rows
  p_ari <- ggplot(d_cluster, aes(x = factor(w), y = factor(k), fill = ari)) +
    geom_tile() +
    scale_fill_viridis_c(name = "ARI", option = "C",
                         limits = c(0, 1), oob = squish) +
    facet_grid(col_short ~ precision, labeller = label_both) +
    labs(x = "w", y = "k",
         title = "Mode D clustering Adjusted Rand Index vs tissue labels",
         subtitle = sprintf("cutree(k=%d tissues), avg linkage on (1-Jaccard)",
                            n_tissues)) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_png(file.path(figures_dir, "mode_d_clustering_ari.png"),
           p_ari, width = max(10, 1.6 * n_p), height = 6)

  p_nmi <- ggplot(d_cluster, aes(x = factor(w), y = factor(k), fill = nmi)) +
    geom_tile() +
    scale_fill_viridis_c(name = "NMI", option = "C",
                         limits = c(0, 1), oob = squish) +
    facet_grid(col_short ~ precision, labeller = label_both) +
    labs(x = "w", y = "k",
         title = "Mode D clustering Normalized Mutual Information vs tissue labels",
         subtitle = sprintf("cutree(k=%d tissues), avg linkage on (1-Jaccard)",
                            n_tissues)) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_png(file.path(figures_dir, "mode_d_clustering_nmi.png"),
           p_nmi, width = max(10, 1.6 * n_p), height = 6)

  # ── Best configs per criterion (reference=bedtools, jaccard_similarity) ──
  pick_best <- function(df, col_name, descending = TRUE) {
    z <- df %>% filter(!is.na(.data[[col_name]]))
    if (descending) z <- z %>% arrange(desc(.data[[col_name]]))
    else            z <- z %>% arrange(.data[[col_name]])
    z %>% slice(1)
  }
  b_pool <- d_bedtools %>% filter(column %in% REG_EQ_VALUES)
  cat("\nBest Mode D configs (column=jaccard_similarity, ref=bedtools):\n")
  for (m in c("pearson", "spearman", "ari", "nmi")) {
    b <- pick_best(b_pool, m, descending = TRUE)
    cat(sprintf("  by %-9s : k=%2d w=%3d p=%2d  -> %s=%.4f  (r_pearson=%.3f, ari=%.3f)\n",
                m, b$k, b$w, b$precision, m, b[[m]],
                b$pearson, b$ari))
  }
  cat("  by MAE     : ")
  b <- pick_best(b_pool, "mae", descending = FALSE)
  cat(sprintf("k=%2d w=%3d p=%2d  -> mae=%.4f\n",
              b$k, b$w, b$precision, b$mae))

  # ── Per-column optimum, read-only ────────────────────────────────────────
  # The selection above drives every published Mode D figure and deliberately
  # stays on jaccard_similarity. This block answers the separate question of
  # which cell each estimator would pick, and writes it to its own file --- it
  # must not feed the scatter, dendrogram or contingency outputs below, which
  # write to fixed filenames (mode_d_best_dendrogram.png is Fig 5a in
  # docs/paper_outline.md).
  #
  # Sourced from d_out, not b_pool: b_pool descends from d_plot, which is
  # already filtered to jaccard_similarity, so re-filtering it is a no-op.
  #
  # Ties are recorded rather than broken. pick_best resolves them via dplyr's
  # stable arrange over a frame sorted by (column, precision, k, w), so under
  # jaccard_similarity_ie the ARI optimum lands on p=12 k=10 w=20 --- a
  # sort artifact, not a preference. That tie is 9-way across precisions, and
  # the jaccard_similarity one is already 8-way, so an arbitrary tie-break over
  # precision is not new to the IE arm.
  maximise <- c(pearson = TRUE, spearman = TRUE, ari = TRUE, nmi = TRUE,
                mae = FALSE)
  best_by_column <- bind_rows(lapply(sort(unique(d_out$column)), function(cl) {
    pool <- d_out %>% filter(reference == "bedtools", column == cl)
    bind_rows(lapply(names(maximise), function(obj) {
      z <- pool %>% filter(!is.na(.data[[obj]]))
      if (nrow(z) == 0) return(NULL)
      target <- if (maximise[[obj]]) max(z[[obj]]) else min(z[[obj]])
      ties <- z %>% filter(abs(.data[[obj]] - target) < 1e-9)
      sel <- ties %>% slice(1)  # == pick_best's stable-sort winner
      data.frame(
        column = cl, objective = obj, value = target,
        k = sel$k, w = sel$w, precision = sel$precision,
        ari = sel$ari, nmi = sel$nmi, pearson = sel$pearson, mae = sel$mae,
        n_tied = nrow(ties),
        tie_set = paste(sprintf("p%d/k%d/w%d", ties$precision, ties$k, ties$w),
                        collapse = " "),
        stringsAsFactors = FALSE
      )
    }))
  }))
  write_csv(best_by_column, file.path(results_dir, "best_config_by_column.csv"))
  cat("Wrote results/best_config_by_column.csv (", nrow(best_by_column), "rows)\n")
  cat("\nPer-column optimum (ref=bedtools; ties recorded, not broken):\n")
  for (i in seq_len(nrow(best_by_column))) {
    r <- best_by_column[i, ]
    cat(sprintf("  %-28s by %-8s : k=%2d w=%3d p=%2d -> %.4f  (%d tied)\n",
                r$column, r$objective, r$k, r$w, r$precision, r$value, r$n_tied))
  }

  # The "best" config used for the scatter + dendrogram below is now the
  # ARI-best (the question the user actually asks: do clusters match
  # the biology?). Falls back to Pearson if ARI is degenerate.
  best <- pick_best(b_pool, "ari", descending = TRUE)
  if (is.na(best$ari) || best$ari == 0) best <- pick_best(b_pool, "pearson", descending = TRUE)
  cat(sprintf("\nUsing config k=%d w=%d p=%d for scatter+dendrogram (ARI=%.3f)\n",
              best$k, best$w, best$precision, best$ari))
  best_d_path <- best$path
  best_d_col  <- best$column

  # ── Best scatter, side-by-side: vs bedtools | vs Mode B ──────────────────
  best_hat <- read_hammock_csv(best$path, jcol = best$column)
  best_bed <- best_hat %>% inner_join(ref,       by = c("stem1", "stem2")) %>%
    mutate(reference = "vs bedtools")
  best_mB  <- if (!is.null(modeB_ref))
    best_hat %>% inner_join(modeB_ref, by = c("stem1", "stem2")) %>%
      mutate(reference = "vs Mode B (p=21)")
  else NULL
  best_long <- bind_rows(best_bed, best_mB)
  best_long$reference <- factor(best_long$reference,
                                levels = c("vs bedtools", "vs Mode B (p=21)"))
  p6 <- ggplot(best_long, aes(x = j_truth, y = j_hat)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey50") +
    geom_point(alpha = 0.6, size = 1.6) +
    scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
    facet_wrap(~ reference) +
    labs(x = "reference Jaccard",
         y = sprintf("Mode D %s", best$column),
         title = sprintf("Mode D k=%d w=%d p=%d (ARI=%.3f, r_bedtools=%.3f)",
                         best$k, best$w, best$precision,
                         best$ari, best$pearson)) +
    theme_minimal(base_size = 11)
  save_png(file.path(figures_dir, "mode_d_best_scatter.png"), p6,
           width = 9, height = 4.5)
} else {
  cat("No Mode D outputs found in", d_dir, "-- skipping.\n")
}

# ---------- dendrograms (tissue-colored labels) ----------
# to_matrix() and tissue_for() are defined near the top of this file.

short_label <- function(stem) {
  # Maurano stems look like "fBrain-DS14718.hotspot..." or
  # "fIntestine_Sm-DS16712.hg19.hotspot...". The organ-and-accession is
  # the part before the first "."; everything after is process metadata
  # (hg19, hotspot, twopass, fdr0.05, merge) we don't want on the leaf.
  sub("\\..*$", "", stem)
}

dendro_png <- function(out_path, mat, title, draw_clusters = FALSE,
                       k_clusters = NULL) {
  d <- as.dist(1 - mat)
  hc <- hclust(d, method = "average")
  tissues <- tissue_for(rownames(mat))
  pal <- setNames(scales::hue_pal()(length(unique(tissues))),
                  unique(tissues))
  # cols is parallel to hc$labels (input order). For visual placement we
  # index by hc$order to convert input position → visual position.
  cols <- pal[tissues[hc$labels]]
  ord  <- hc$order
  hc$labels <- short_label(hc$labels)
  CairoPNG(filename = out_path, width = 11, height = 6,
           units = "in", res = 150)
  on.exit(dev.off())
  op <- par(mar = c(6, 4, 3, 10))   # bigger right margin for the legend
  on.exit(par(op), add = TRUE)
  plot(hc, hang = -1, labels = FALSE,
       main = title, xlab = "", sub = "", cex = 0.85)
  # Coloured leaf labels in visual order.
  mtext(text = hc$labels[ord], side = 1, at = seq_along(ord),
        col = cols[ord], las = 2, line = 0.5, cex = 0.8)
  if (draw_clusters) {
    rect.hclust(hc, k = if (is.null(k_clusters)) n_tissues else k_clusters,
                border = "steelblue")
  }
  # Tissue → colour legend in the right margin (xpd = NA so it can sit
  # outside the plot region).
  legend(x = par("usr")[2] + diff(par("usr")[1:2]) * 0.02,
         y = par("usr")[4],
         legend = names(pal), fill = pal,
         border = NA, bty = "n", cex = 0.75, xpd = NA,
         title = "tissue", title.adj = 0)
}

ref_mat <- to_matrix(ref %>% transmute(stem1, stem2, j_hat = j_truth))
dendro_png(file.path(figures_dir, "bedtools_dendrogram.png"), ref_mat,
           "bedtools pairwise Jaccard (reference)",
           draw_clusters = TRUE)

if (!is.null(best_d_path)) {
  d_mat <- to_matrix(read_hammock_csv(best_d_path, jcol = best_d_col))

  # Recompute the predicted clusters at the ARI-best config so we can
  # record them and show them on the dendrogram.
  hc_best <- hclust(as.dist(1 - {
    m <- d_mat; m[m < 0] <- 0; m[m > 1] <- 1; m
  }), method = "average")
  pred <- cutree(hc_best, k = n_tissues)
  true <- tissue_label_vec[names(pred)]
  best_ari_val <- adjusted_rand(true, pred)
  best_nmi_val <- nmi(true, pred)

  # Per-sample assignment table.
  assign_df <- tibble(
    stem               = names(pred),
    short_label        = short_label(names(pred)),
    true_tissue        = unname(true),
    predicted_cluster  = unname(pred)
  ) %>% arrange(predicted_cluster, true_tissue)
  write_csv(assign_df, file.path(results_dir, "best_cluster_assignment.csv"))
  cat("Wrote results/best_cluster_assignment.csv\n")

  # Contingency table (predicted cluster × true tissue).
  cont <- as.data.frame.matrix(table(predicted_cluster = pred, true_tissue = true))
  cont_out <- cbind(predicted_cluster = rownames(cont), cont)
  write_csv(as_tibble(cont_out),
            file.path(results_dir, "best_cluster_contingency.csv"))
  cat("Wrote results/best_cluster_contingency.csv\n")

  # Confusion matrix as a meaningful heatmap (predicted × true).
  cont_long <- as.data.frame(table(pred = pred, tissue = true)) %>%
    mutate(pred = factor(pred), tissue = factor(tissue))
  p_conf <- ggplot(cont_long, aes(x = tissue, y = pred, fill = Freq)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(Freq == 0, "", Freq)), size = 3.5) +
    scale_fill_gradient(low = "#f7f7f7", high = "#2c7bb6") +
    labs(x = "true tissue", y = "predicted cluster",
         fill = "n samples",
         title = sprintf("Mode D cluster assignment (k=%d w=%d p=%d, %s)",
                         best$k, best$w, best$precision, best_d_col),
         subtitle = sprintf("ARI = %.3f, NMI = %.3f against tissue labels",
                            best_ari_val, best_nmi_val)) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_png(file.path(figures_dir, "mode_d_cluster_confusion.png"),
           p_conf, width = 8, height = 5)

  # Dendrogram with predicted-cluster boxes + ARI/NMI in title.
  dendro_png(
    file.path(figures_dir, "mode_d_best_dendrogram.png"),
    d_mat,
    sprintf("Mode D k=%d w=%d p=%d (%s)  —  ARI = %.3f, NMI = %.3f",
            best$k, best$w, best$precision,
            sub(paste0(REG_EQ_COL, "|", REG_EQ_LEGACY_COL), "no_ends", best_d_col),
            best_ari_val, best_nmi_val),
    draw_clusters = TRUE
  )
}

cat("Done.\n")
