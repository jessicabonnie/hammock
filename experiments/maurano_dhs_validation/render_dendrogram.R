#!/usr/bin/env Rscript
# Render a tissue-coloured dendrogram from a single Mode D CSV.
#
# Usage:
#   Rscript render_dendrogram.R <input_csv_or_tsv> <out_png> \
#       [jaccard_column] [title]
#
# Default jaccard column is `jaccard_similarity` (no-ends). Default title
# is "Mode D dendrogram (<stem>, <jaccard_column>)". Pass a 4th arg to
# override the title (e.g. "bedtools pairwise Jaccard (reference)").
# Input is read as TSV if its extension is .tsv, else CSV.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(scales)
  library(Cairo)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: render_dendrogram.R <input_csv_or_tsv> <out_png> ",
       "[jaccard_column] [title]")
}
csv_path  <- args[1]
out_png   <- args[2]
jcol      <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "jaccard_similarity"
title_arg <- if (length(args) >= 4 && nzchar(args[4])) args[4] else NULL

script_dir <- dirname(normalizePath(
  sub("--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE),
                          value = TRUE)[1])))
data_dir <- file.path(script_dir, "data")

strip_ext <- function(x) sub("\\.(bed|fa|fasta|fna)$", "", x, ignore.case = TRUE)

key_path <- file.path(data_dir, "maurano_filenames_key.tsv")
key <- if (file.exists(key_path)) {
  read_tsv(key_path, show_col_types = FALSE) %>%
    transmute(file = File, tissue = Biosample_term_name)
} else {
  tibble(file = character(), tissue = character())
}

tissue_for <- function(stems) {
  if (nrow(key) == 0) return(setNames(rep("?", length(stems)), stems))
  setNames(key$tissue[match(stems, strip_ext(key$file))], stems)
}

# ARI / NMI implementations (ported from analyze.R — no external deps).
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

to_matrix <- function(df) {
  files <- sort(unique(c(df$stem1, df$stem2)))
  m <- matrix(NA_real_, nrow = length(files), ncol = length(files),
              dimnames = list(files, files))
  for (i in seq_len(nrow(df))) m[df$stem1[i], df$stem2[i]] <- df$j_hat[i]
  for (i in seq_along(files)) for (j in seq_along(files)) {
    if (is.na(m[i, j]) && !is.na(m[j, i])) m[i, j] <- m[j, i]
  }
  diag(m) <- 1
  m[m < 0] <- 0; m[m > 1] <- 1
  m
}

# Leaf label = accession only. Stems look like:
#   "fBrain-DS14718.hotspot.twopass.fdr0.05.merge"
#   "fIntestine_Sm-DS16712.hg19.hotspot.twopass.fdr0.05.merge"
# Pull "DS#####" — the substring between the first "-" and the next ".".
short_label <- function(stem) sub("^[^-]*-([^.]+)\\..*$", "\\1", stem)

# Tissue display name.
#   - Strip the leading "f" (Maurano fetal prefix): fBrain -> Brain.
#   - Muscle_arm / Muscle_back / Muscle_leg -> "Muscle (Arm)" etc.
#   - Other underscored names (Intestine_Sm, Kidney_renal_cortex_L,
#     Skin_fibro_bicep_R, ...) collapse to the first underscore-delimited
#     token.
display_tissue <- function(t) {
  t <- sub("^f", "", t)
  vapply(t, function(x) {
    if (grepl("^Muscle_", x)) {
      sub_part <- sub("^Muscle_", "", x)
      sub_part <- paste0(toupper(substr(sub_part, 1, 1)),
                         substr(sub_part, 2, nchar(sub_part)))
      paste0("Muscle (", sub_part, ")")
    } else {
      sub("_.*$", "", x)
    }
  }, character(1), USE.NAMES = FALSE)
}

dendro_png <- function(out_path, mat, title) {
  d <- as.dist(1 - mat)
  hc <- hclust(d, method = "average")
  tissues <- tissue_for(rownames(mat))
  uniq <- unique(tissues)
  pal <- setNames(scales::hue_pal()(length(uniq)), uniq)
  cols <- pal[tissues[hc$labels]]
  ord  <- hc$order
  # Snapshot tissue-per-visual-position BEFORE renaming hc$labels — the
  # `tissues` dict is keyed on the full original labels.
  tissue_at_pos <- tissues[hc$labels[ord]]
  # Cluster-vs-tissue agreement at k = n_tissues, the same cut analyze.R
  # uses for the mode_d_summary.csv ARI/NMI columns.
  pred <- cutree(hc, k = length(uniq))
  true <- tissues[hc$labels]
  ari_val <- adjusted_rand(true, pred)
  nmi_val <- nmi(true, pred)
  hc$labels <- short_label(hc$labels)

  # Per-tissue contiguous-run boxes. For each tissue, find its visual
  # positions, split into maximal runs of consecutive positions, draw
  # one box per run sized to the MRCA merge height of the contained
  # leaves (cophenetic distance — height at which all leaves in the run
  # first merge into a single cluster). Single-leaf runs get a small
  # fixed box height so they remain visible.
  coph        <- as.matrix(stats::cophenetic(hc))
  max_h       <- max(hc$height)
  singleton_h <- max_h * 0.05

  CairoPNG(filename = out_path, width = 11, height = 6,
           units = "in", res = 150)
  on.exit(dev.off())
  op <- par(mar = c(6, 4, 3, 10))
  on.exit(par(op), add = TRUE)
  plot(hc, hang = -1, labels = FALSE,
       main = title, xlab = "", sub = "", cex = 0.85)
  subtitle <- sprintf("ARI = %.3f   NMI = %.3f   (cut at k = %d tissues)",
                      ari_val, nmi_val, length(uniq))
  mtext(subtitle, side = 3, line = 0.1, cex = 0.9)
  mtext(text = hc$labels[ord], side = 1, at = seq_along(ord),
        col = cols[ord], las = 2, line = 0.5, cex = 0.8)
  for (t in unique(tissue_at_pos)) {
    pos <- which(tissue_at_pos == t)
    if (length(pos) == 0) next
    run_id <- cumsum(c(1, diff(pos) != 1))
    for (run in split(pos, run_id)) {
      leaves <- ord[run]
      box_top <- if (length(leaves) >= 2) max(coph[leaves, leaves]) * 1.02
                 else singleton_h
      rect(xleft = min(run) - 0.4, ybottom = par("usr")[3],
           xright = max(run) + 0.4, ytop = box_top,
           border = pal[[t]], lwd = 2, xpd = NA)
    }
  }
  legend(x = par("usr")[2] + diff(par("usr")[1:2]) * 0.02,
         y = par("usr")[4],
         legend = display_tissue(names(pal)), fill = pal,
         border = NA, bty = "n", cex = 0.75, xpd = NA,
         title = "Tissue", title.adj = 0)
}

df <- if (grepl("\\.tsv$", csv_path, ignore.case = TRUE)) {
  read_tsv(csv_path, show_col_types = FALSE)
} else {
  read_csv(csv_path, show_col_types = FALSE)
}
if (!(jcol %in% names(df))) {
  stop("column ", jcol, " not in ", csv_path,
       " (have: ", paste(names(df), collapse = ", "), ")")
}
j_hat <- df %>% transmute(stem1 = strip_ext(file1),
                          stem2 = strip_ext(file2),
                          j_hat = .data[[jcol]])
mat <- to_matrix(j_hat)
stem <- tools::file_path_sans_ext(basename(csv_path))
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
title <- if (!is.null(title_arg)) {
  title_arg
} else {
  paste0("Mode D dendrogram (", stem, ", ", jcol, ")")
}
dendro_png(out_png, mat, title)
cat("Wrote ", out_png, "\n", sep = "")
