#!/usr/bin/env Rscript

# Figure 6 — Sequence-mode recovery of fetal-tissue organization
#
# This is the paper-local version of the tissue-coloured, cluster-boxed
# dendrogram produced in experiments/maurano_dhs_validation/analyze.R.
# It intentionally preserves the original visual logic while isolating only
# the code and inputs needed to reproduce the manuscript figure.
#
# Default analysis:
#   similarity = jaccard_similarity_ie (inclusion-exclusion Jaccard)
#   k = 10, w = 30, p = 18
#   linkage = average
#   number of displayed clusters = number of annotated tissue labels (10)
#
# Usage:
#   Rscript paper/sequence_tissue_clustering/plot_sequence_tissue_clustering.R
#
# Optional overrides:
#   Rscript ... <similarity_csv> <tissue_key.tsv> <output.png> [similarity_column] [panel_label]
#
# `jaccard_similarity_ie` is accepted as the similarity column even though the
# archived sweep predates it: the CSVs carry containment_AB/containment_BA, from
# which it is exactly recoverable (see jaccard_ie_from_containments below). Every
# run also writes estimator_agreement_stats.csv next to this script, scoring both
# columns on the same clustering and recording whether they induce the same
# partition.

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
P <- 18

experiment_dir <- file.path(repo_root, "experiments", "maurano_dhs_validation")
default_csv <- file.path(
  experiment_dir, "results", "raw_d",
  # This archived extended-sweep output predates the current _full filename
  # tag but contains the directional containments needed to reconstruct IE.
  sprintf("hammock_mnmzr_p%d_jaccD_k%d_w%d.csv", P, K, W)
)
default_key <- file.path(experiment_dir, "data", "maurano_filenames_key.tsv")
default_output <- file.path(repo_root, "paper", "figures", "sequence_tissue_clustering.png")

argv <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(argv) >= 1) argv[1] else default_csv
key_tsv <- if (length(argv) >= 2) argv[2] else default_key
out_png <- if (length(argv) >= 3) argv[3] else default_output
# Resolved below after `raw` is loaded. An explicit CLI override always wins;
# otherwise the figure uses inclusion-exclusion Jaccard.
sim_col_arg <- if (length(argv) >= 4 && nzchar(argv[4])) argv[4] else NA_character_
panel_label <- if (length(argv) >= 5 && nzchar(argv[5])) argv[5] else "A"

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

# Name written under each box. Short by design — a box can be a single leaf, and
# the full legend text ("Kidney (renal cortex)") would overrun its neighbours.
GROUP_DISPLAY <- c(
  fBrain = "Brain",
  fHeart = "Heart",
  fIntestine_Sm = "Small\nintestine",
  fKidney_renal_cortex_L = "Kidney",
  fLung = "Lung",
  fMuscle_arm = "Arm",
  fMuscle_back = "Back",
  fMuscle_leg = "Leg",
  fSkin_fibro_bicep_R = "Skin",
  fStomach = "Stomach"
)

# Tissues that share an organ carry only their body part above, and the organ is
# named once in a bracket spanning their boxes.
BRACKET_GROUP <- c(
  fMuscle_arm = "Muscle",
  fMuscle_back = "Muscle",
  fMuscle_leg = "Muscle"
)

display_group <- function(x) {
  out <- unname(GROUP_DISPLAY[x])
  ifelse(is.na(out), display_tissue(x), out)
}

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

jaccard_ie_from_containments <- function(c_ab, c_ba) {
  # Set-Jaccard from the two directional containments:
  #   1 / (1/C_AB + 1/C_BA - 1) == |A n B| / (|A| + |B| - |A n B|).
  # Mirrors python/hammock/runner.py::_jaccard_ie_from_containments, clamp
  # included, and matches experiments/maurano_dhs_validation/analyze.R:166-171.
  # The clamp absorbs Ertl noise that can push a containment a few ulp past 1
  # and is what guarantees denom >= 1; it does not fire on this corpus (no
  # containment in the 235 archived raw_d CSVs exceeds 1.0), but the extreme
  # size-ratio regime it guards against is real elsewhere. A zero containment
  # means the intersection estimate was zero -- genuinely empty, or clamped
  # from a negative -- and is scored 0.0.
  cab <- pmin(as.numeric(c_ab), 1)
  cba <- pmin(as.numeric(c_ba), 1)
  ifelse(cab <= 0 | cba <= 0, 0, 1 / (1 / cab + 1 / cba - 1))
}

# Oracle values lifted from tests/test_jaccard_ie.py, which pins the canonical
# Python helper. Cheap, and it localises a clamp bug that the end-to-end ARI
# check downstream would only show as a wrong number.
local({
  probe <- jaccard_ie_from_containments(
    c(0.5, 1.0, 1.0000000000050957, 0.0, 0.5),
    c(0.5, 1.0, 1.0,                0.5, 0.0)
  )
  stopifnot(
    isTRUE(all.equal(probe[1], 1 / 3)),
    probe[2] == 1, probe[3] == 1, probe[4] == 0, probe[5] == 0
  )
})

# The archived Mode D sweep predates the jaccard_similarity_ie column (v0.5.0)
# but already carries the containments it is derived from, so recover it here
# rather than resweeping. Native _ie is passed through untouched when present.
add_jaccard_ie <- function(df) {
  if (!("jaccard_similarity_ie" %in% names(df)) &&
      all(c("containment_AB", "containment_BA") %in% names(df))) {
    df[["jaccard_similarity_ie"]] <-
      jaccard_ie_from_containments(df$containment_AB, df$containment_BA)
  }
  df
}

key <- read_tsv(key_tsv, show_col_types = FALSE) %>%
  transmute(stem = strip_ext(File), tissue = Biosample_term_name)

raw <- add_jaccard_ie(read_csv(input_csv, show_col_types = FALSE))

resolve_sim_col <- function(explicit, df) {
  if (!is.na(explicit)) return(explicit)
  "jaccard_similarity_ie"
}

# Resolve register equality separately for the estimator-agreement diagnostic.
# Archived pre-rename CSVs use jaccard_similarity for this legacy statistic.
resolve_reg_eq_col <- function(df) {
  if ("reg_eq_similarity" %in% names(df)) return("reg_eq_similarity")
  message(
    "reg_eq_similarity not found in input CSV; falling back to jaccard_similarity"
  )
  "jaccard_similarity"
}
reg_eq_col <- resolve_reg_eq_col(raw)
sim_col <- resolve_sim_col(sim_col_arg, raw)

required_cols <- c("file1", "file2", sim_col)
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("Input lacks columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

stems <- sort(unique(c(strip_ext(raw$file1), strip_ext(raw$file2))))
missing_stems <- setdiff(stems, key$stem)
if (length(missing_stems) > 0) {
  stop("Tissue key is missing: ", paste(missing_stems, collapse = ", "), call. = FALSE)
}
tissue_by_stem <- setNames(key$tissue, key$stem)

similarity_matrix <- function(df, column) {
  pairs <- df %>%
    transmute(
      stem1 = strip_ext(file1),
      stem2 = strip_ext(file2),
      similarity = .data[[column]]
    )
  m <- matrix(NA_real_, length(stems), length(stems), dimnames = list(stems, stems))
  for (i in seq_len(nrow(pairs))) {
    m[pairs$stem1[i], pairs$stem2[i]] <- pairs$similarity[i]
  }
  m[is.na(m)] <- t(m)[is.na(m)]
  if (anyNA(m)) stop("Similarity matrix is incomplete.", call. = FALSE)
  diag(m) <- 1
  m[] <- pmin(pmax(m, 0), 1)
  m
}

mat <- similarity_matrix(raw, sim_col)
hc <- hclust(as.dist(1 - mat), method = "average")
true_tissue <- tissue_by_stem[hc$labels]
n_tissues <- length(unique(true_tissue))
predicted <- cutree(hc, k = n_tissues)
ari <- adjusted_rand(true_tissue, predicted)
nmi <- normalized_mi(true_tissue, predicted)

# ---- estimator agreement ---------------------------------------------------
# Figure 6 is drawn on inclusion-exclusion Jaccard by default. Record whether
# its k = n_tissues partition differs from the archived register-equality
# column so the local estimator agreement remains explicit.
# cutree's cluster ids are arbitrary, so partitions are compared as sets of
# member sets, not elementwise. Deliberately keyed on reg_eq_col here, not
# sim_col: this comparison must stay register-equality-vs-IE regardless of
# which column the S5a invocation asks the dendrogram itself to be drawn on
# -- see reg_eq_col's own comment above.
partition_signature <- function(p) {
  paste(sort(vapply(split(names(p), p),
                    function(g) paste(sort(g), collapse = ","),
                    character(1))),
        collapse = " | ")
}

agreement <- bind_rows(lapply(
  intersect(c(reg_eq_col, "jaccard_similarity_ie"), names(raw)),
  function(cl) {
    h <- hclust(as.dist(1 - similarity_matrix(raw, cl)), method = "average")
    p <- cutree(h, k = n_tissues)
    truth <- tissue_by_stem[h$labels]
    data.frame(
      column = cl,
      ari = adjusted_rand(truth, p),
      nmi = normalized_mi(truth, p),
      signature = partition_signature(p),
      stringsAsFactors = FALSE
    )
  }
))

ref_signature <- agreement$signature[agreement$column == reg_eq_col]
agreement$partition_identical <- if (length(ref_signature) == 1) {
  agreement$signature == ref_signature
} else {
  NA
}
agreement$signature <- NULL
agreement <- cbind(
  data.frame(input = basename(input_csv), n_clusters = n_tissues,
             stringsAsFactors = FALSE),
  agreement
)

stats_csv <- file.path(script_dir, "estimator_agreement_stats.csv")
write_csv(agreement, stats_csv)
message("Wrote: ", stats_csv)
print(agreement, row.names = FALSE)

# Match the original experiment figure: one colour per tissue, tissue-coloured
# leaf labels, and rectangles outlining contiguous same-tissue label runs.
tissues <- tissue_by_stem[hc$labels]
tissue_levels <- unique(tissues)
tissue_palette <- setNames(hue_pal()(length(tissue_levels)), tissue_levels)
label_colors <- tissue_palette[tissues]
ord <- hc$order

# Boxes group the leaf labels, not the tree itself: each contiguous run of one
# tissue in the leaf order gets a box drawn around its labels in the bottom
# margin. Grouping is by annotated tissue rather than by the k = n_tissues cut,
# which can split a genuine clade (the two fBrain samples land in separate
# singleton clusters at k = 10).
coph <- as.matrix(cophenetic(hc))
ordered_stems <- hc$labels[ord]
ordered_groups <- tissues[ord]
group_runs <- split(
  seq_along(ordered_groups),
  cumsum(c(1, diff(as.integer(factor(ordered_groups, levels = unique(ordered_groups)))) != 0))
)

# Diagnostic only — a boxed label run that is not monophyletic means the tree
# interleaves that tissue with others even though its labels sit side by side.
for (run in group_runs) {
  if (length(run) < 2) next
  inside <- ordered_stems[run]
  outside <- setdiff(ordered_stems, inside)
  if (length(outside) > 0 && max(coph[inside, inside]) >= min(coph[inside, outside])) {
    warning(sprintf(
      "Group '%s' is not monophyletic; its labels are boxed but its clade is not.",
      ordered_groups[run[1]]
    ))
  }
}

hc$labels <- short_label(hc$labels)

CairoPNG(
  filename = out_png,
  width = 11,
  height = 6,
  units = "in",
  res = 300,
  bg = "white"
)
op <- par(mar = c(7, 4, 2, 2), family = "sans")
on.exit({ par(op); dev.off() }, add = TRUE)

plot(
  hc,
  hang = -1,
  labels = FALSE,
  main = "",
  xlab = "",
  sub = "",
  ylab = "1 − Jaccard",
  cex = 0.85
)
mtext(panel_label, side = 3, adj = 0, line = 0.15, cex = 16 / 12, font = 2)

LABEL_CEX <- 0.8
GROUP_CEX <- 0.7

# Everything below the tree is placed relative to the leaf tips (y = 0) rather
# than to the plot region edge: par("usr")[3] sits ~4% of the height below the
# tips, and hanging the labels off it left a wide empty band. Vertical distances
# are specified in inches and converted, so they hold at any device size.
usr <- par("usr")
user_y <- function(inches) inches * diff(usr[3:4]) / par("pin")[2]

gap <- user_y(0.10)  # tree tips to the top of the label block
pad <- user_y(0.04)  # box padding around the labels
label_top <- -gap
# Labels are rotated, so their vertical extent is the string *width*.
label_bottom <- label_top -
  user_y(max(strwidth(hc$labels, units = "inches", cex = LABEL_CEX)))

text(
  x = seq_along(ord),
  y = label_top,
  labels = hc$labels[ord],
  col = label_colors[ord],
  srt = 90,
  adj = c(1, 0.5),
  cex = LABEL_CEX,
  xpd = NA
)

box_top <- label_top + pad
box_bottom <- label_bottom - pad

# Each box takes the colour of the tissue it encompasses; a run is one tissue by
# construction, so the first leaf's colour is the run's colour.
group_colors <- label_colors[ord][vapply(group_runs, function(run) run[1], integer(1))]

for (i in seq_along(group_runs)) {
  run <- group_runs[[i]]
  rect(
    xleft = min(run) - 0.45,
    ybottom = box_bottom,
    xright = max(run) + 0.45,
    ytop = box_top,
    border = group_colors[i],
    xpd = NA
  )
}

# Tissue name under each box. The names are slanted because several boxes are a
# single leaf wide and a horizontal name would overrun its neighbours.
group_tissues <- vapply(group_runs, function(run) ordered_groups[run[1]], character(1))
group_labels <- display_group(group_tissues)
group_centers <- vapply(group_runs, function(run) mean(range(run)), numeric(1))

# Apply the organ bracket only within contiguous runs of muscle subtypes. A
# separated subtype (as for Arm at p=12) receives a complete two-line label
# instead of being connected by a bracket across unrelated tissues.
muscle_groups <- which(group_tissues %in% names(BRACKET_GROUP))
group_brackets <- rep(NA_character_, length(group_tissues))
if (length(muscle_groups) > 0) {
  muscle_runs <- split(
    muscle_groups,
    cumsum(c(1, diff(muscle_groups) != 1))
  )
  for (members in muscle_runs) {
    if (length(members) >= 2) {
      group_brackets[members] <- "Muscle"
    } else {
      subtype <- sub("^fMuscle_", "", group_tissues[members])
      group_labels[members] <- paste0(
        toupper(substr(subtype, 1, 1)), substr(subtype, 2, nchar(subtype)),
        "\nmuscle"
      )
    }
  }
}

user_x <- function(inches) inches * diff(usr[1:2]) / par("pin")[1]
line_h <- par("cin")[2] * GROUP_CEX * 1.15  # inches per line of a name
name_parts <- strsplit(group_labels, "\n", fixed = TRUE)
name_widths <- lapply(name_parts, strwidth, units = "inches", cex = GROUP_CEX)
name_lines <- lengths(name_parts)

# Slanted names run parallel to each other, so a name overlaps its right-hand
# neighbour unless the horizontal gap between them, projected onto the
# perpendicular, clears its block: gap * sin(angle) > lines * line_h. A two-line
# name beside a single-leaf box is the binding case, so the angle is derived
# rather than fixed, and only steepens past 45 when it has to.
anchor_gaps <- diff(group_centers) * par("pin")[1] / diff(usr[1:2])
needed_sine <- if (length(anchor_gaps) > 0) {
  max(tail(name_lines, -1) * line_h / anchor_gaps) * 1.15
} else 0
GROUP_ANGLE <- min(70, max(45, asin(min(1, needed_sine)) * 180 / pi))
slant <- GROUP_ANGLE * pi / 180

# Every name ends at the same point relative to its box — just under the box
# centre — and hangs down the slant from there. Drawn one line at a time:
# text() right-justifies the lines of a multi-line string, which would push the
# short "Small" up the slant into the leaf labels. Left-aligning the lines
# instead leaves the short one hanging away from them.
#
# Multi-line names are stacked so their *last* line sits where a one-line name
# would, with earlier lines above it. Stacking downwards instead would leave
# "intestine" a full line further from the leaves than every other name.
name_anchor <- box_bottom - pad

for (i in seq_along(group_labels)) {
  widths <- name_widths[[i]]
  for (j in seq_along(name_parts[[i]])) {
    along <- widths[j] - max(widths)  # short lines start back down the slant
    perp <- (name_lines[i] - j - 0.5) * line_h  # last line at the box, rest above
    text(
      x = group_centers[i] + user_x(along * cos(slant) - perp * sin(slant)),
      y = name_anchor + user_y(along * sin(slant) + perp * cos(slant)),
      labels = name_parts[[i]][j],
      srt = GROUP_ANGLE,
      adj = c(1, 0.5),
      cex = GROUP_CEX,
      xpd = NA
    )
  }
}

# Organ bracket. Ends tick up towards the names it spans; the organ word sits in
# a gap left in the middle of the horizontal line.
for (organ in unique(group_brackets[!is.na(group_brackets)])) {
  members <- which(group_brackets %in% organ)
  leaves <- unlist(group_runs[members])
  # Clear the deepest name in the span: the lower-left corner of its last line.
  depth <- max(vapply(members, function(i) {
    max(name_widths[[i]]) * sin(slant) + line_h * cos(slant)
  }, numeric(1)))
  y_line <- name_anchor - user_y(depth + 0.08)
  tick <- user_y(0.06)
  x_left <- min(leaves) - 0.45
  x_right <- max(leaves) + 0.45
  x_mid <- (x_left + x_right) / 2
  gap_half <- user_x(strwidth(organ, units = "inches", cex = GROUP_CEX) / 2 + 0.05)
  segments(
    x0 = c(x_left, x_right, x_left, x_mid + gap_half),
    y0 = c(y_line, y_line, y_line, y_line),
    x1 = c(x_left, x_right, x_mid - gap_half, x_right),
    y1 = c(y_line + tick, y_line + tick, y_line, y_line),
    xpd = NA
  )
  text(x_mid, y_line, organ, adj = c(0.5, 0.5), cex = GROUP_CEX, xpd = NA)
}

message(sprintf("ARI = %.3f; NMI = %.3f", ari, nmi))
message("Wrote: ", out_png)
