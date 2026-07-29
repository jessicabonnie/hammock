#!/usr/bin/env Rscript

# Figure 6 — Sequence-mode recovery of fetal-tissue organization
#
# This is the paper-local version of the tissue-coloured, cluster-boxed
# dendrogram produced in experiments/maurano_dhs_validation/analyze.R.
# It intentionally preserves the original visual logic while isolating only
# the code and inputs needed to reproduce the manuscript figure.
#
# Default analysis:
#   similarity = jaccard_similarity (minimizer-only / no ends)
#   k = 10, w = 30, p = 24
#   linkage = average
#   number of displayed clusters = number of annotated tissue labels (10)
#
# Usage:
#   Rscript paper/sequence_tissue_clustering/plot_sequence_tissue_clustering.R
#
# Optional overrides:
#   Rscript ... <similarity_csv> <tissue_key.tsv> <output.png> [similarity_column]

required_packages <- c("dplyr", "readr", "scales", "Cairo")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(scales)
  library(Cairo)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
  stop("Could not determine the script path. Run with Rscript.", call. = FALSE)
}
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

K <- 10
W <- 30
P <- 24
DEFAULT_SIM_COL <- "jaccard_similarity"

experiment_dir <- file.path(repo_root, "experiments", "maurano_dhs_validation")
default_csv <- file.path(
  experiment_dir, "results", "raw_d",
  sprintf("hammock_mnmzr_p%d_jaccD_k%d_w%d.csv", P, K, W)
)
default_key <- file.path(experiment_dir, "data", "maurano_filenames_key.tsv")
default_output <- file.path(repo_root, "paper", "figures", "sequence_tissue_clustering.png")

argv <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(argv) >= 1) argv[1] else default_csv
key_tsv <- if (length(argv) >= 2) argv[2] else default_key
out_png <- if (length(argv) >= 3) argv[3] else default_output
sim_col <- if (length(argv) >= 4 && nzchar(argv[4])) argv[4] else DEFAULT_SIM_COL

dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
for (path in c(input_csv, key_tsv)) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
}

strip_ext <- function(x) {
  sub("\\.(bed|fa|fasta|fna)$", "", basename(x), ignore.case = TRUE)
}

short_label <- function(stem) {
  # Leaves are labelled by accession only; tissue identity is carried by the
  # label colour and the legend.
  sub("^[^-]*-", "", sub("\\..*$", "", stem))
}

# Readable legend text for the ENCODE Biosample_term_name values. Anything not
# listed falls back to a generic cleanup (drop the fetal "f", underscores to
# spaces) so a new tissue never breaks the figure.
TISSUE_DISPLAY <- c(
  fBrain = "Brain",
  fHeart = "Heart",
  fIntestine_Sm = "Small intestine",
  fKidney_renal_cortex_L = "Kidney (renal cortex)",
  fLung = "Lung",
  fMuscle_arm = "Muscle (arm)",
  fMuscle_back = "Muscle (back)",
  fMuscle_leg = "Muscle (leg)",
  fSkin_fibro_bicep_R = "Skin fibroblast (bicep)",
  fStomach = "Stomach"
)

display_tissue <- function(x) {
  out <- unname(TISSUE_DISPLAY[x])
  fallback <- gsub("_", " ", sub("^f(?=[A-Z])", "", x, perl = TRUE))
  fallback <- paste0(toupper(substr(fallback, 1, 1)), substr(fallback, 2, nchar(fallback)))
  ifelse(is.na(out), fallback, out)
}

# Coarse organ grouping used for the cluster boxes: subtypes of the same organ
# (fMuscle_arm/_back/_leg) are boxed as one group.
tissue_group <- function(x) sub("_.*$", "", sub("^f(?=[A-Z])", "", x, perl = TRUE))

adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  choose2 <- function(x) sum(choose(x, 2))
  index <- choose2(as.vector(tab))
  expected <- choose2(rowSums(tab)) * choose2(colSums(tab)) / choose(n, 2)
  maximum <- (choose2(rowSums(tab)) + choose2(colSums(tab))) / 2
  if (isTRUE(all.equal(maximum, expected))) return(0)
  (index - expected) / (maximum - expected)
}

normalized_mi <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  pi <- rowSums(tab) / n
  pj <- colSums(tab) / n
  pij <- tab / n
  h_a <- -sum(pi[pi > 0] * log(pi[pi > 0]))
  h_b <- -sum(pj[pj > 0] * log(pj[pj > 0]))
  if (h_a == 0 || h_b == 0) return(0)
  outer_p <- outer(pi, pj)
  mask <- pij > 0 & outer_p > 0
  mutual_info <- sum(pij[mask] * log(pij[mask] / outer_p[mask]))
  2 * mutual_info / (h_a + h_b)
}

key <- read_tsv(key_tsv, show_col_types = FALSE) %>%
  transmute(stem = strip_ext(File), tissue = Biosample_term_name)

raw <- read_csv(input_csv, show_col_types = FALSE)
required_cols <- c("file1", "file2", sim_col)
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("Input lacks columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

pairs <- raw %>%
  transmute(
    stem1 = strip_ext(file1),
    stem2 = strip_ext(file2),
    similarity = .data[[sim_col]]
  )

stems <- sort(unique(c(pairs$stem1, pairs$stem2)))
missing_stems <- setdiff(stems, key$stem)
if (length(missing_stems) > 0) {
  stop("Tissue key is missing: ", paste(missing_stems, collapse = ", "), call. = FALSE)
}

mat <- matrix(NA_real_, length(stems), length(stems), dimnames = list(stems, stems))
for (i in seq_len(nrow(pairs))) {
  mat[pairs$stem1[i], pairs$stem2[i]] <- pairs$similarity[i]
}
mat[is.na(mat)] <- t(mat)[is.na(mat)]
if (anyNA(mat)) stop("Similarity matrix is incomplete.", call. = FALSE)
diag(mat) <- 1
mat[] <- pmin(pmax(mat, 0), 1)

hc <- hclust(as.dist(1 - mat), method = "average")
tissue_by_stem <- setNames(key$tissue, key$stem)
true_tissue <- tissue_by_stem[hc$labels]
n_tissues <- length(unique(true_tissue))
predicted <- cutree(hc, k = n_tissues)
ari <- adjusted_rand(true_tissue, predicted)
nmi <- normalized_mi(true_tissue, predicted)

# Match the original experiment figure: one colour per tissue, tissue-coloured
# leaf labels, and blue rectangles outlining contiguous organ-level clades.
tissues <- tissue_by_stem[hc$labels]
tissue_levels <- unique(tissues)
tissue_palette <- setNames(hue_pal()(length(tissue_levels)), tissue_levels)
label_colors <- tissue_palette[tissues]
ord <- hc$order

# Boxes mark organ-level clades rather than the k = n_tissues cut. The cut can
# split a genuine clade (the two fBrain samples land in separate singleton
# clusters at k = 10) and it separates the muscle subtypes, so each contiguous
# run of one organ group gets a single box drawn up to that run's MRCA height.
coph <- as.matrix(cophenetic(hc))
group_by_leaf <- tissue_group(tissues)
ordered_stems <- hc$labels[ord]
ordered_groups <- group_by_leaf[ord]
group_runs <- split(
  seq_along(ordered_groups),
  cumsum(c(1, diff(as.integer(factor(ordered_groups, levels = unique(ordered_groups)))) != 0))
)
max_height <- max(hc$height)

box_tops <- vapply(group_runs, function(run) {
  if (length(run) < 2) return(max_height * 0.05)
  inside <- ordered_stems[run]
  mrca <- max(coph[inside, inside])
  outside <- setdiff(ordered_stems, inside)
  # A non-monophyletic run would be drawn as a box enclosing foreign leaves.
  if (length(outside) > 0 && mrca >= min(coph[inside, outside])) {
    warning(sprintf(
      "Group '%s' is not monophyletic; its box spans leaves from other groups.",
      ordered_groups[run[1]]
    ))
  }
  mrca * 1.03
}, numeric(1))

hc$labels <- short_label(hc$labels)

CairoPNG(
  filename = out_png,
  width = 11,
  height = 6,
  units = "in",
  res = 300,
  bg = "white"
)
op <- par(mar = c(6, 4, 3, 10), family = "sans")
on.exit({ par(op); dev.off() }, add = TRUE)

plot(
  hc,
  hang = -1,
  labels = FALSE,
  main = "Sequence sketches recover fetal-tissue organization",
  xlab = "",
  sub = "",
  ylab = "1 − Jaccard",
  cex = 0.85
)

mtext(
  text = hc$labels[ord],
  side = 1,
  at = seq_along(ord),
  col = label_colors[ord],
  las = 2,
  line = 0.5,
  cex = 0.8
)

for (i in seq_along(group_runs)) {
  run <- group_runs[[i]]
  rect(
    xleft = min(run) - 0.45,
    ybottom = par("usr")[3],
    xright = max(run) + 0.45,
    ytop = box_tops[i],
    border = "steelblue"
  )
}

legend(
  x = par("usr")[2] + diff(par("usr")[1:2]) * 0.02,
  y = par("usr")[4],
  legend = display_tissue(names(tissue_palette)),
  fill = tissue_palette,
  border = NA,
  bty = "n",
  cex = 0.75,
  xpd = NA,
  title = "Fetal tissue",
  title.adj = 0
)

message(sprintf("ARI = %.3f; NMI = %.3f", ari, nmi))
message("Wrote: ", out_png)
