#!/usr/bin/env Rscript

# Figure 6 — Sequence-mode recovery of fetal-tissue organization
# Adapts the Maurano DHS validation dendrogram for the manuscript.

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

# ARI-optimal sequence-mode configuration from the Maurano sweep.
K <- 10
W <- 30
P <- 12
SIM_COL <- "jaccard_similarity"
EXPECTED_ARI <- 0.910
EXPECTED_NMI <- 0.961

# -----------------------------------------------------------------------------
# Paper palette interface
# -----------------------------------------------------------------------------
# Replace these defaults when a paper-wide palette is finalized. Every color
# used by this figure is obtained from PAPER_PALETTE or TISSUE_COLORS below.
# A future shared palette file can define the same objects and be sourced here.
PAPER_PALETTE <- list(
  text = "#20262D",
  branch = "#46515C",
  axis = "#6B747D",
  cluster_outline = "#8A939C",
  background = "#FFFFFF"
)

TISSUE_COLORS <- c(
  fBrain = "#4477AA",
  fHeart = "#CC6677",
  fIntestine_Sm = "#DDCC77",
  fKidney = "#117733",
  fLung = "#88CCEE",
  fMuscle_arm = "#AA4499",
  fMuscle_back = "#882255",
  fMuscle_leg = "#CC99BB",
  fSkin = "#44AA99",
  fStomach = "#999933"
)

shared_palette <- file.path(repo_root, "paper", "palette.R")
if (file.exists(shared_palette)) {
  palette_env <- new.env(parent = emptyenv())
  sys.source(shared_palette, envir = palette_env)
  if (exists("PAPER_PALETTE", envir = palette_env, inherits = FALSE)) {
    PAPER_PALETTE <- get("PAPER_PALETTE", envir = palette_env)
  }
  if (exists("TISSUE_COLORS", envir = palette_env, inherits = FALSE)) {
    TISSUE_COLORS <- get("TISSUE_COLORS", envir = palette_env)
  }
}

experiment_dir <- file.path(repo_root, "experiments", "maurano_dhs_validation")
default_csv <- file.path(
  experiment_dir, "results", "raw_d",
  sprintf("hammock_mnmzr_p%d_jaccD_k%d_w%d.csv", P, K, W)
)
default_key <- file.path(experiment_dir, "data", "maurano_filenames_key.tsv")

argv <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(argv) >= 1) normalizePath(argv[1], mustWork = FALSE) else default_csv
key_tsv <- if (length(argv) >= 2) normalizePath(argv[2], mustWork = FALSE) else default_key
out_png <- if (length(argv) >= 3) {
  normalizePath(argv[3], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "sequence_tissue_clustering.png")
}
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

for (path in c(input_csv, key_tsv)) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
}

strip_ext <- function(x) sub("\\.(bed|fa|fasta|fna)$", "", basename(x), ignore.case = TRUE)
short_label <- function(stem) sub("^[^-]*-([^.]+)\\..*$", "\\1", stem)

display_tissue <- function(x) {
  x <- sub("^f", "", x)
  vapply(x, function(value) {
    if (grepl("^Muscle_", value)) {
      part <- sub("^Muscle_", "", value)
      part <- paste0(toupper(substr(part, 1, 1)), substr(part, 2, nchar(part)))
      paste0("Muscle (", part, ")")
    } else if (grepl("^Intestine_Sm", value)) {
      "Small intestine"
    } else {
      sub("_.*$", "", value)
    }
  }, character(1), USE.NAMES = FALSE)
}

adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
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
  pi <- rowSums(tab) / n
  pj <- colSums(tab) / n
  pij <- tab / n
  h_a <- -sum(pi[pi > 0] * log(pi[pi > 0]))
  h_b <- -sum(pj[pj > 0] * log(pj[pj > 0]))
  outer_p <- outer(pi, pj)
  mask <- pij > 0 & outer_p > 0
  mutual_info <- sum(pij[mask] * log(pij[mask] / outer_p[mask]))
  2 * mutual_info / (h_a + h_b)
}

key <- read_tsv(key_tsv, show_col_types = FALSE) %>%
  transmute(stem = strip_ext(File), tissue = Biosample_term_name)

raw <- read_csv(input_csv, show_col_types = FALSE)
required_cols <- c("file1", "file2", SIM_COL)
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("Input lacks columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

pairs <- raw %>%
  transmute(
    stem1 = strip_ext(file1),
    stem2 = strip_ext(file2),
    similarity = .data[[SIM_COL]]
  )

stems <- sort(unique(c(pairs$stem1, pairs$stem2)))
missing_stems <- setdiff(stems, key$stem)
if (length(missing_stems) > 0) {
  stop("Tissue key is missing: ", paste(missing_stems, collapse = ", "), call. = FALSE)
}

mat <- matrix(NA_real_, length(stems), length(stems), dimnames = list(stems, stems))
for (i in seq_len(nrow(pairs))) mat[pairs$stem1[i], pairs$stem2[i]] <- pairs$similarity[i]
mat[is.na(mat)] <- t(mat)[is.na(mat)]
if (anyNA(mat)) stop("Similarity matrix is incomplete.", call. = FALSE)
diag(mat) <- 1
mat[] <- pmin(pmax(mat, 0), 1)

hc <- hclust(as.dist(1 - mat), method = "average")
tissue_by_stem <- setNames(key$tissue, key$stem)
true_tissue <- tissue_by_stem[hc$labels]
n_clusters <- length(unique(true_tissue))
predicted <- cutree(hc, k = n_clusters)
ari <- adjusted_rand(true_tissue, predicted)
nmi <- normalized_mi(true_tissue, predicted)

if (abs(ari - EXPECTED_ARI) > 0.005 || abs(nmi - EXPECTED_NMI) > 0.005) {
  warning(sprintf(
    "Observed ARI/NMI (%.3f/%.3f) differ from expected manuscript values (%.3f/%.3f).",
    ari, nmi, EXPECTED_ARI, EXPECTED_NMI
  ))
}

ordered_stems <- hc$labels[hc$order]
ordered_tissues <- tissue_by_stem[ordered_stems]
ordered_clusters <- predicted[ordered_stems]
unique_tissues <- unique(key$tissue)

missing_colors <- setdiff(unique_tissues, names(TISSUE_COLORS))
if (length(missing_colors) > 0) {
  fallback <- hue_pal(l = 58, c = 90)(length(missing_colors))
  TISSUE_COLORS <- c(TISSUE_COLORS, setNames(fallback, missing_colors))
}
label_colors <- TISSUE_COLORS[ordered_tissues]
hc$labels <- short_label(hc$labels)

# Build one subtle outline per contiguous predicted-cluster run. The outline
# rises to the MRCA height of the leaves in the run, matching the dendrogram cut
# without competing with the tissue colors.
coph <- as.matrix(cophenetic(hc))
max_h <- max(hc$height)
singleton_h <- max_h * 0.045
run_id <- cumsum(c(1, diff(ordered_clusters) != 0))
cluster_runs <- split(seq_along(ordered_clusters), run_id)

CairoPNG(
  filename = out_png,
  width = 10.5,
  height = 6.8,
  units = "in",
  res = 300,
  bg = PAPER_PALETTE$background
)
op <- par(mar = c(7.2, 4.5, 3.1, 9.5), family = "sans")
on.exit({ par(op); dev.off() }, add = TRUE)

plot(
  hc,
  hang = -1,
  labels = FALSE,
  main = "Sequence sketches recover fetal-tissue organization",
  xlab = "",
  sub = "",
  ylab = "1 − Jaccard",
  axes = TRUE,
  cex.axis = 0.9,
  col = PAPER_PALETTE$branch
)

for (run in cluster_runs) {
  run_stems <- ordered_stems[run]
  box_top <- if (length(run) >= 2) {
    max(coph[run_stems, run_stems, drop = FALSE]) * 1.02
  } else {
    singleton_h
  }
  rect(
    xleft = min(run) - 0.42,
    ybottom = par("usr")[3],
    xright = max(run) + 0.42,
    ytop = box_top,
    border = adjustcolor(PAPER_PALETTE$cluster_outline, alpha.f = 0.72),
    lwd = 1.15,
    xpd = NA
  )
}

mtext(
  hc$labels[hc$order],
  side = 1,
  at = seq_along(hc$order),
  col = label_colors,
  las = 2,
  line = 0.45,
  cex = 0.80
)

legend(
  x = par("usr")[2] + diff(par("usr")[1:2]) * 0.025,
  y = par("usr")[4],
  legend = display_tissue(unique_tissues),
  fill = TISSUE_COLORS[unique_tissues],
  border = NA,
  bty = "n",
  cex = 0.72,
  xpd = NA,
  title = "Tissue",
  title.adj = 0
)

message(sprintf("ARI = %.3f; NMI = %.3f", ari, nmi))
message("Wrote: ", out_png)
