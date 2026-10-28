#!/usr/bin/env Rscript
# Figure 2: hammock workflow diagram, drawn with grid so text stays vector and
# prints at 7--9 pt at a 6.5 in text width. Writes SVG (editable), PDF, and
# PNG versions of the same drawing.
#
# Usage: Rscript plot_workflow_diagram.R [output_prefix]

suppressPackageStartupMessages({
  library(grid)
  library(Cairo)
})

argv <- commandArgs(trailingOnly = TRUE)
out_prefix <- if (length(argv) >= 1) argv[1] else "hammock_workflow"

W <- 7.0
H <- 5.10

SANS <- "Nimbus Sans"
MONO <- "Nimbus Mono PS"
COL_TEXT <- "#20262D"
COL_MUTED <- "#5A6570"
COL_TRACK <- "#8A949E"
COL_INT <- "#2C6DB0"       # interval mode
COL_INT_FILL <- "#E8F0F9"
COL_SEQ <- "#2A7F5C"       # sequence mode
COL_SEQ_FILL <- "#E7F3EC"
COL_BAR <- "#B03A2E"

FS_TAG <- 11
FS_TITLE <- 9.5
FS_STEP <- 8.0
FS_BODY <- 7.4

# --- helpers ---------------------------------------------------------------
# All coordinates are inches from the bottom-left of the canvas.
txt <- function(label, x, y, size = FS_BODY, col = COL_TEXT, face = "plain",
                hjust = 0.5, vjust = 0.5, family = SANS, rot = 0) {
  grid.text(label, x = unit(x, "in"), y = unit(y, "in"),
            hjust = hjust, vjust = vjust, rot = rot,
            gp = gpar(fontsize = size, col = col, fontface = face,
                      fontfamily = family))
}

box <- function(x0, y0, x1, y1, col = COL_TEXT, fill = NA, lwd = 0.8,
                r = 0.06, lty = 1) {
  grid.roundrect(x = unit(x0, "in"), y = unit(y0, "in"),
                 width = unit(x1 - x0, "in"), height = unit(y1 - y0, "in"),
                 just = c("left", "bottom"), r = unit(r, "in"),
                 gp = gpar(col = col, fill = fill, lwd = lwd, lty = lty))
}

rect <- function(x0, y0, x1, y1, col = NA, fill = "grey", lwd = 0.5) {
  grid.rect(x = unit(x0, "in"), y = unit(y0, "in"),
            width = unit(x1 - x0, "in"), height = unit(y1 - y0, "in"),
            just = c("left", "bottom"),
            gp = gpar(col = col, fill = fill, lwd = lwd))
}

seg <- function(x0, y0, x1, y1, col = COL_TRACK, lwd = 0.7, lty = 1) {
  grid.segments(unit(x0, "in"), unit(y0, "in"), unit(x1, "in"),
                unit(y1, "in"), gp = gpar(col = col, lwd = lwd, lty = lty))
}

arrow_right <- function(x0, x1, y, col = COL_TEXT) {
  grid.segments(unit(x0, "in"), unit(y, "in"), unit(x1, "in"), unit(y, "in"),
                arrow = arrow(length = unit(0.06, "in"), type = "closed"),
                gp = gpar(col = col, fill = col, lwd = 1.0))
}

step_header <- function(n, label, x, y, col) {
  grid.circle(x = unit(x, "in"), y = unit(y, "in"), r = unit(0.075, "in"),
              gp = gpar(col = NA, fill = col))
  txt(n, x, y - 0.003, size = 6.8, col = "white", face = "bold")
  lines <- strsplit(label, "\n")[[1]]
  for (i in seq_along(lines)) {
    txt(lines[i], x + 0.12, y - (i - 1) * 0.13, size = FS_STEP, face = "bold",
        hjust = 0)
  }
}

doc_icon <- function(x, y, w = 0.26, h = 0.32, label = "BED", col = COL_INT) {
  fold <- 0.07
  grid.polygon(
    x = unit(c(x, x + w - fold, x + w, x + w, x), "in"),
    y = unit(c(y + h, y + h, y + h - fold, y, y), "in"),
    gp = gpar(col = COL_TEXT, fill = "white", lwd = 0.7)
  )
  rect(x, y + 0.03, x + w, y + 0.13, fill = col)
  txt(label, x + w / 2, y + 0.08, size = 5.8, col = "white", face = "bold")
}

# Interval track with blocks at the given [start, end] pairs (in inches).
track <- function(x0, x1, y, blocks, col = COL_INT, h = 0.07) {
  seg(x0, y, x1, y)
  for (b in blocks) rect(b[1], y - h / 2, b[2], y + h / 2, fill = col)
}

# Covered positions drawn as per-base ticks inside each block.
tick_track <- function(x0, x1, y, blocks, col = COL_INT, step = 0.025,
                       h = 0.08, keep_every = 1, grid_x0 = NULL,
                       faint = "#DCE5EF") {
  seg(x0, y, x1, y)
  for (b in blocks) {
    # Snap ticks to a global grid so both files share retained positions.
    origin <- if (is.null(grid_x0)) x0 else grid_x0
    idx <- seq(ceiling((b[1] - origin) / step), floor((b[2] - origin) / step))
    for (ii in idx) {
      xi <- origin + ii * step
      kept <- is.null(grid_x0) || ii %% keep_every == 0
      seg(xi, y - h / 2, xi, y + h / 2, col = if (kept) col else faint,
          lwd = if (kept) 1.3 else 0.4)
    }
  }
}

registers <- function(x0, y, n = 10, w = 0.085, h = 0.1, col = COL_INT,
                      seed = 1) {
  set.seed(seed)
  shades <- colorRampPalette(c("white", col))(6)
  vals <- sample(2:6, n, replace = TRUE)
  for (i in seq_len(n)) {
    rect(x0 + (i - 1) * w, y, x0 + i * w, y + h, col = COL_TEXT,
         fill = shades[vals[i]], lwd = 0.4)
  }
}

heatmap <- function(x0, y0, size, mat, col) {
  n <- nrow(mat)
  cell <- size / n
  pal <- colorRampPalette(c("white", col))(101)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    rect(x0 + (j - 1) * cell, y0 + (n - i) * cell,
         x0 + j * cell, y0 + (n - i + 1) * cell,
         col = "white", fill = pal[round(mat[i, j] * 100) + 1], lwd = 0.3)
  }
  rect(x0, y0, x0 + size, y0 + size, col = COL_TEXT, fill = NA, lwd = 0.5)
}

panel <- function(x0, y0, x1, y1, tag, title, col, fill) {
  box(x0, y0, x1, y1, col = col, fill = NA, lwd = 1.0)
  txt(tag, x0 + 0.08, y1 - 0.13, size = FS_TAG, face = "bold", hjust = 0)
  txt(title, x0 + 0.28, y1 - 0.13, size = FS_TITLE, face = "bold",
      hjust = 0, col = col)
}

footer <- function(x0, x1, y, question, note, col, fill) {
  box(x0, y - 0.12, x1, y + 0.12, col = col, fill = fill, lwd = 0.6, r = 0.04)
  txt(question, x0 + 0.12, y, size = FS_BODY, hjust = 0)
  txt(note, x1 - 0.12, y, size = FS_BODY, hjust = 1, col = COL_MUTED,
      face = "italic")
}

# --- drawing ---------------------------------------------------------------
draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(col = NA, fill = "white"))

  # Step columns shared by panels B and C.
  cols <- c(0.18, 1.55, 2.92, 4.29, 5.66)
  colw <- 1.18
  arrow_at <- function(i, y, col) {
    arrow_right(cols[i] + colw + 0.02, cols[i + 1] - 0.05, y, col = COL_TEXT)
  }

  # ---- Panel A: inputs and barriers --------------------------------------
  pa <- c(0.05, 3.98, 6.95, 5.05)
  panel(pa[1], pa[2], pa[3], pa[4], "A",
        "Inputs and barriers to systematic comparison", COL_TEXT, NA)
  sub_y <- 4.72

  # (i) interval files defined on different assemblies
  refs <- c("GRCh37", "GRCh38", "CHM13")
  for (i in seq_along(refs)) {
    x <- 0.22 + (i - 1) * 0.82
    doc_icon(x, 4.17, w = 0.24, h = 0.30)
    txt(refs[i], x + 0.29, 4.32, size = FS_BODY, hjust = 0, col = COL_INT,
        face = "bold")
  }
  txt("Interval files on different assemblies",
      0.22, sub_y, size = FS_BODY, hjust = 0, col = COL_MUTED)

  # (ii) quadratic growth in comparisons
  gx <- 2.85
  txt(expression("Comparing" ~ italic(N) ~ "files with" ~ italic(N) ~ "files"),
      gx, sub_y, size = FS_BODY, hjust = 0, col = COL_MUTED)
  heatmap(gx, 4.12, 0.42,
          matrix(c(1, .3, .2, .1, .3, 1, .4, .2, .2, .4, 1, .3,
                   .1, .2, .3, 1), 4), col = "#7F8B96")
  txt(expression(italic(N)^2 ~ "pairwise comparisons"), gx + 0.52, 4.33,
      size = FS_BODY, hjust = 0)

  # (iii) the same locus has different coordinates on each assembly
  rx <- 4.55
  tx0 <- 5.38
  tx1 <- 6.70
  y37 <- 4.46
  y38 <- 4.20
  txt("The same locus has different coordinates", rx, sub_y,
      size = FS_BODY, hjust = 0, col = COL_MUTED)
  txt("GRCh37 chr1", rx, y37, size = 7, hjust = 0, col = COL_INT,
      face = "bold")
  txt("GRCh38 chr1", rx, y38, size = 7, hjust = 0, col = COL_INT,
      face = "bold")
  ticks <- c(5.50, 6.05, 6.60)
  tick_lab <- c("1.20 Mb", "1.25 Mb", "1.30 Mb")
  seg(tx0, y37, tx1, y37)
  seg(tx0, y38, tx1, y38)
  for (i in seq_along(ticks)) {
    seg(ticks[i], y37, ticks[i], y37 + 0.04)
    txt(tick_lab[i], ticks[i], y37 + 0.10, size = 7, col = COL_MUTED)
    seg(ticks[i], y38, ticks[i], y38 - 0.04)
    txt(tick_lab[i], ticks[i], y38 - 0.10, size = 7, col = COL_MUTED)
  }
  h <- 0.07
  rect(5.62, y37 - h / 2, 5.82, y37 + h / 2, fill = COL_INT)
  rect(6.20, y38 - h / 2, 6.40, y38 + h / 2, fill = COL_INT)
  grid.segments(unit(5.74, "in"), unit(y37 - 0.05, "in"),
                unit(6.26, "in"), unit(y38 + 0.05, "in"),
                arrow = arrow(length = unit(0.05, "in"), type = "closed",
                              ends = "both"),
                gp = gpar(col = COL_BAR, fill = COL_BAR, lwd = 0.9, lty = 2))

  # The two modes are alternatives, not consecutive steps.
  txt(expression(italic("hammock") ~ "compares files in one of two alternative modes"),
      3.5, 3.84, size = FS_STEP, face = "bold")

  # ---- Panel B: interval mode --------------------------------------------
  pb <- c(0.05, 1.98, 6.95, 3.70)
  panel(pb[1], pb[2], pb[3], pb[4], "B",
        "Interval mode: coordinate similarity within a shared reference",
        COL_INT, COL_INT_FILL)
  hy <- 3.38
  step_header("1", "BED\nintervals", cols[1], hy, COL_INT)
  step_header("2", "Covered\npositions", cols[2], hy, COL_INT)
  step_header("3", "Subsampling\n(optional)", cols[3], hy, COL_INT)
  step_header("4", "HyperLogLog\nsketch", cols[4], hy, COL_INT)
  step_header("5", "All-pairs\ncomparison", cols[5], hy, COL_INT)

  ya <- 2.94
  yb <- 2.72
  mid <- (ya + yb) / 2
  blocksA <- function(x) list(c(x + 0.10, x + 0.42), c(x + 0.62, x + 0.92))
  blocksB <- function(x) list(c(x + 0.25, x + 0.55), c(x + 0.72, x + 1.02))

  # step 1
  x <- cols[1] + 0.12
  txt("A", cols[1] + 0.02, ya, size = 6.8, col = COL_MUTED)
  txt("B", cols[1] + 0.02, yb, size = 6.8, col = COL_MUTED)
  track(x, x + 1.04, ya, blocksA(x))
  track(x, x + 1.04, yb, blocksB(x))
  txt("one file per sample", cols[1] + 0.59, 2.50, col = COL_MUTED)
  arrow_at(1, mid, COL_INT)

  # step 2
  x <- cols[2] + 0.06
  tick_track(x, x + 1.04, ya, blocksA(x))
  tick_track(x, x + 1.04, yb, blocksB(x))
  txt("each covered base", cols[2] + 0.59, 2.55, col = COL_MUTED)
  txt("is one element", cols[2] + 0.59, 2.42, col = COL_MUTED)
  arrow_at(2, mid, COL_INT)

  # step 3: shared grid; every fourth position retained
  x <- cols[3] + 0.06
  tick_track(x, x + 1.04, ya, blocksA(x), keep_every = 4, grid_x0 = x)
  tick_track(x, x + 1.04, yb, blocksB(x), keep_every = 4, grid_x0 = x)
  for (gxi in seq(x, x + 1.04, by = 0.1)) {
    seg(gxi, yb - 0.08, gxi, ya + 0.08, col = "#C9D3DC", lwd = 0.4, lty = 3)
  }
  txt("same coordinate grid", cols[3] + 0.59, 2.55, col = COL_MUTED)
  txt("in every file", cols[3] + 0.59, 2.42, col = COL_MUTED)
  arrow_at(3, mid, COL_INT)

  # step 4
  x <- cols[4] + 0.10
  txt("A", cols[4] + 0.02, ya, size = 6.8, col = COL_MUTED)
  txt("B", cols[4] + 0.02, yb, size = 6.8, col = COL_MUTED)
  registers(x, ya - 0.05, seed = 3)
  registers(x, yb - 0.05, seed = 4)
  txt("fixed size, built once", cols[4] + 0.59, 2.55, col = COL_MUTED)
  txt("per file and reused", cols[4] + 0.59, 2.42, col = COL_MUTED)
  arrow_at(4, mid, COL_INT)

  # step 5
  set.seed(7)
  m <- matrix(runif(36, 0.05, 0.45), 6)
  m <- (m + t(m)) / 2
  diag(m) <- 1
  heatmap(cols[5] + 0.33, 2.62, 0.48, m, COL_INT)
  txt("Jaccard estimate", cols[5] + 0.57, 2.50, col = COL_MUTED)
  txt("for every pair", cols[5] + 0.57, 2.37, col = COL_MUTED)

  footer(pb[1] + 0.1, pb[3] - 0.1, 2.14,
         "Question answered: do these files cover similar genomic positions?",
         "Requires a shared reference", COL_INT, COL_INT_FILL)

  # ---- Panel C: sequence mode --------------------------------------------
  pc <- c(0.05, 0.05, 6.95, 1.90)
  panel(pc[1], pc[2], pc[3], pc[4], "C",
        "Sequence mode: sequence similarity across references",
        COL_SEQ, COL_SEQ_FILL)
  hy <- 1.58
  step_header("1", "BED + its\nreference FASTA", cols[1], hy, COL_SEQ)
  step_header("2", "Interval\nsequences", cols[2], hy, COL_SEQ)
  step_header("3", "Minimizer\nselection", cols[3], hy, COL_SEQ)
  step_header("4", "HyperLogLog\nsketch", cols[4], hy, COL_SEQ)
  step_header("5", "Cross-reference\ncomparison", cols[5], hy, COL_SEQ)
  mid <- 1.02

  # step 1
  doc_icon(cols[1] + 0.06, 0.98, w = 0.24, h = 0.30, label = "BED",
           col = COL_INT)
  txt("+", cols[1] + 0.42, 1.13, size = 11, face = "bold")
  doc_icon(cols[1] + 0.54, 0.98, w = 0.30, h = 0.30, label = "FASTA",
           col = COL_SEQ)
  txt("each list with its", cols[1] + 0.48, 0.80, col = COL_MUTED)
  txt("reference assembly", cols[1] + 0.48, 0.67, col = COL_MUTED)
  arrow_at(1, mid + 0.11, COL_SEQ)

  # step 2
  box(cols[2] + 0.04, 0.86, cols[2] + 1.16, 1.32, col = COL_TEXT,
      fill = "white", lwd = 0.6, r = 0.03)
  fa <- c(">chr1:1200-1290", "ACCTGAGTCGATCCGA…",
          ">chr1:5400-5480", "TGCTTGCAAGGCTTAC…")
  for (i in seq_along(fa)) {
    txt(fa[i], cols[2] + 0.09, 1.26 - (i - 1) * 0.11, size = 6.6,
        family = MONO, hjust = 0,
        col = if (startsWith(fa[i], ">")) COL_MUTED else COL_TEXT)
  }
  txt("one record", cols[2] + 0.60, 0.74, col = COL_MUTED)
  txt("per interval", cols[2] + 0.60, 0.62, col = COL_MUTED)
  arrow_at(2, mid + 0.11, COL_SEQ)

  # step 3: one window of w k-mers; the smallest-hash k-mer is retained
  sx <- cols[3] + 0.07
  seq_str <- strsplit("ACCTGAGTCGATCC", "")[[1]]
  cw <- 0.075
  for (i in seq_along(seq_str)) {
    txt(seq_str[i], sx + (i - 0.5) * cw, 1.27, size = 7, family = MONO)
  }
  k <- 5
  nk <- 6
  for (j in seq_len(nk)) {
    yk <- 1.18 - (j - 1) * 0.048
    sel <- j == 4
    rect(sx + (j - 1) * cw + 0.004, yk - 0.018, sx + (j - 1 + k) * cw - 0.004,
         yk + 0.018, fill = if (sel) COL_SEQ else "#C7D9CF")
  }
  wx0 <- sx
  wx1 <- sx + (nk - 1 + k) * cw
  seg(wx0, 0.89, wx1, 0.89, col = COL_TEXT, lwd = 0.6)
  seg(wx0, 0.89, wx0, 0.92, col = COL_TEXT, lwd = 0.6)
  seg(wx1, 0.89, wx1, 0.92, col = COL_TEXT, lwd = 0.6)
  txt(expression("window of" ~ italic(w) ~ italic(k) * "-mers;"),
      (wx0 + wx1) / 2, 0.74, col = COL_MUTED)
  txt("smallest hash retained", (wx0 + wx1) / 2, 0.62, col = COL_SEQ)
  arrow_at(3, mid + 0.11, COL_SEQ)

  # step 4
  x <- cols[4] + 0.10
  registers(x, 1.10, col = COL_SEQ, seed = 5)
  registers(x, 0.93, col = COL_SEQ, seed = 6)
  txt("fixed size, built once", cols[4] + 0.59, 0.74, col = COL_MUTED)
  txt("per file and reused", cols[4] + 0.59, 0.62, col = COL_MUTED)
  arrow_at(4, mid + 0.11, COL_SEQ)

  # step 5: block-structured matrix; files group by content, not assembly
  b <- matrix(0.25, 6, 6)
  b[1:3, 1:3] <- 0.85
  b[4:6, 4:6] <- 0.85
  diag(b) <- 1
  heatmap(cols[5] + 0.33, 0.86, 0.48, b, COL_SEQ)
  txt("files group by content,", cols[5] + 0.57, 0.74, col = COL_MUTED)
  txt("not by assembly", cols[5] + 0.57, 0.62, col = COL_MUTED)

  footer(pc[1] + 0.1, pc[3] - 0.1, 0.30,
         "Question answered: do these files contain similar underlying sequence?",
         "No shared coordinates required", COL_SEQ, COL_SEQ_FILL)

  # "or" badge between the alternative mode panels
  box(3.28, 1.85, 3.72, 2.03, col = COL_TEXT, fill = "white", lwd = 0.8,
      r = 0.09)
  txt("or", 3.5, 1.94, size = FS_STEP, face = "bold")
}

CairoSVG(paste0(out_prefix, ".svg"), width = W, height = H, bg = "white")
draw_figure()
invisible(dev.off())

CairoPDF(paste0(out_prefix, ".pdf"), width = W, height = H, bg = "white")
draw_figure()
invisible(dev.off())

CairoPNG(paste0(out_prefix, ".png"), width = W, height = H, units = "in",
         res = 300, bg = "white")
draw_figure()
invisible(dev.off())

message("Wrote: ", out_prefix, ".{svg,pdf,png}")
