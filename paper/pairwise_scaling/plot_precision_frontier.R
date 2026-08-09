# Precision frontier -- accuracy bought against speed paid, one point per p.
#
# PLOTTED DIRECTLY, NOT ON TWIN AXES. The obvious rendering is wall time and MAE
# against p on left and right y-axes, but twin log axes have no canonical
# registration: sliding one relative to the other moves the apparent crossing
# anywhere you like, so "the curves cross at p=18" would be a statement about the
# axis limits, not the data. Putting accuracy on x and speed on y removes the
# free parameter. Iso-accuracy is then a vertical line, and "slower than
# bedtools" is the region below y=1, which is a fact about the data alone.
#
# WHAT THIS PANEL SAYS, AND THE TRAP IN IT. At N=20 hammock is at best ~1.1x
# bedtools, and at p>=18 it is slower. That is not in tension with Figure 3A's
# 27.6x at N=512 -- it is the same crossover seen from the other axis. Maurano is
# 20 files, well below the N~64 crossover, so this panel is measuring the regime
# where sketching cost has not yet been amortised over enough pairs. Say that in
# the caption; a reader who takes the 27.6x and this panel as competing claims
# has been misled by the figure, not by the tools.
#
# The x-axis is MAE of jaccard_similarity_ie, NOT jaccard_similarity. The
# register-equality column carries a chance-agreement floor and sits at ~0.138
# MAE at every precision -- it does not improve with p at all, so plotting it
# here would draw a vertical line and imply precision buys nothing. Both are
# printed to stdout so the contrast is visible; see CLAUDE.md divergence #2.

required_packages <- c("dplyr", "readr", "ggplot2", "scales", "Cairo")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2)
  library(scales); library(Cairo)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Run with Rscript.", call. = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

argv <- commandArgs(trailingOnly = TRUE)
in_csv <- if (length(argv) >= 1) argv[1] else
  file.path(repo_root, "docs", "data", "sweep_precision_maurano_p18_t16.csv")
alt_csv <- if (length(argv) >= 2) argv[2] else
  file.path(repo_root, "docs", "data", "sweep_precision_maurano_p18_t8.csv")
out_png <- if (length(argv) >= 3) argv[3] else
  file.path(repo_root, "paper", "figures", "precision_frontier.png")

if (!file.exists(in_csv)) stop("Input not found: ", in_csv, call. = FALSE)

# --- reduce to medians, and check the run is what it claims to be -----------
load_one <- function(path) {
  raw <- read_csv(path, show_col_types = FALSE)
  for (need in c("tool", "precision", "wall_time", "jaccard_ie_mae_vs_bt",
                 "jaccard_mae_vs_bt", "jaccard_n_pairs")) {
    if (!need %in% names(raw))
      stop("Column '", need, "' missing from ", basename(path),
           " -- this CSV predates the accuracy/provenance columns.", call. = FALSE)
  }

  bt_rows <- raw %>% filter(tool == "bedtools")
  if (nrow(bt_rows) == 0)
    stop("No bedtools rows in ", basename(path), "; there is no denominator ",
         "for a speedup and this panel cannot be drawn.", call. = FALSE)
  bt <- median(bt_rows$wall_time)

  hm <- raw %>%
    filter(tool != "bedtools", !is.na(precision)) %>%
    group_by(precision) %>%
    summarise(wall = median(wall_time),
              lo = min(wall_time), hi = max(wall_time),
              mae_ie = median(jaccard_ie_mae_vs_bt),
              mae_reg = median(jaccard_mae_vs_bt),
              n_pairs = median(jaccard_n_pairs),
              threads = dplyr::first(threads),
              .groups = "drop") %>%
    mutate(speedup = bt / wall, lo_s = bt / hi, hi_s = bt / lo,
           bt_wall = bt)

  # Self-pairs must be excluded or the MAE is diluted by 20 exact-1.0 rows.
  np <- unique(hm$n_pairs)
  if (length(np) != 1)
    stop("Pair count varies across precisions in ", basename(path),
         " (", paste(np, collapse = ", "), "); the accuracy series is not ",
         "comparable across p.", call. = FALSE)
  if (np != 380)
    stop("Expected 380 ordered pairs (400 minus 20 self-pairs), got ", np,
         " in ", basename(path), ". If this is 400 the self-pair filter did ",
         "not fire and every MAE here is ~5% low.", call. = FALSE)
  hm
}

d <- load_one(in_csv)
threads_used <- unique(d$threads)

cat(sprintf("\n=== precision frontier, %d files, %d ordered pairs, t=%s ===\n",
            20L, as.integer(unique(d$n_pairs)), paste(threads_used, collapse = "/")))
cat(sprintf("bedtools reference: %.2f s\n\n", d$bt_wall[1]))
d %>% transmute(p = precision, wall = round(wall, 2), speedup = round(speedup, 2),
                mae_ie = signif(mae_ie, 4), mae_regeq = signif(mae_reg, 4)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# The register-equality column is flat in p -- state it numerically rather than
# letting a reader assume the two estimators behave alike.
reg_spread <- max(d$mae_reg) / min(d$mae_reg)
ie_spread  <- max(d$mae_ie)  / min(d$mae_ie)
cat(sprintf(
  "\njaccard_similarity_ie MAE improves %.0fx from p=%d to p=%d.\n",
  ie_spread, min(d$precision), max(d$precision)))
cat(sprintf(
  "jaccard_similarity (register-equality) MAE varies only %.2fx over the same\n",
  reg_spread))
cat("range and never falls below 0.13 -- it is dominated by the chance-agreement\n")
cat("floor, not by sampling noise, so precision cannot fix it. Hence the x-axis.\n\n")

# --- consistency check against the other thread setting ---------------------
# Same job, same node, so this is a real check and not a cross-allocation guess.
if (file.exists(alt_csv)) {
  a <- load_one(alt_csv)
  j <- inner_join(d %>% select(precision, mae_ie),
                  a %>% select(precision, mae_ie), by = "precision",
                  suffix = c("_main", "_alt"))
  worst <- max(abs(j$mae_ie_main - j$mae_ie_alt))
  cat(sprintf("Accuracy cross-check vs %s: max |dMAE| = %.3g",
              basename(alt_csv), worst))
  if (worst > 1e-12) {
    cat("  <-- UNEXPECTED\n")
    cat("*** Accuracy must be identical across thread counts: the sketch is\n",
        "*** deterministic and threads only change how fast it is built. A\n",
        "*** nonzero difference means a race or a seed leak, not noise.\n", sep = "")
  } else {
    cat(" (bit-identical, as it must be)\n")
  }
  cat(sprintf("Speed at t=%s vs t=%s, p=18: %.2fx vs %.2fx\n\n",
              paste(unique(d$threads), collapse = ""),
              paste(unique(a$threads), collapse = ""),
              d$speedup[d$precision == 18], a$speedup[a$precision == 18]))
}

# --- plot -------------------------------------------------------------------
lab <- d %>% mutate(txt = paste0("p=", precision))
default_p <- d %>% filter(precision == 18)

xlo <- min(d$mae_ie) / 1.45; xhi <- max(d$mae_ie) * 1.75
ylo <- min(d$speedup) / 1.30; yhi <- max(d$speedup) * 1.12

fig <- ggplot(d, aes(mae_ie, speedup)) +
  annotate("rect", xmin = xlo / 10, xmax = xhi * 10,
           ymin = ylo / 10, ymax = 1, fill = "grey80", alpha = 0.5) +
  annotate("text", x = xhi / 1.05, y = 1, vjust = 1.8, hjust = 1,
           label = "slower than BEDTools", size = 3.1, colour = "grey25") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40",
             linewidth = 0.4) +
  geom_linerange(aes(ymin = lo_s, ymax = hi_s), colour = "#1b7f7f",
                 linewidth = 0.5, alpha = 0.8) +
  geom_path(colour = "#1b7f7f", linewidth = 0.6, alpha = 0.8) +
  geom_point(colour = "#1b7f7f", size = 2.2) +
  geom_point(data = default_p, shape = 21, size = 5, stroke = 0.9,
             colour = "#b8420f", fill = NA) +
  geom_text(data = lab, aes(label = txt), size = 3.2, hjust = -0.35,
            vjust = -0.55, colour = "grey20") +
  annotate("text", x = default_p$mae_ie, y = default_p$speedup,
           label = "  CLI default", hjust = 0, vjust = 2.6, size = 3.1,
           colour = "#b8420f") +
  scale_x_log10(breaks = breaks_log(n = 6),
                labels = label_scientific(digits = 2)) +
  scale_y_log10(breaks = breaks_log(n = 6),
                labels = label_number(accuracy = 0.01, drop0trailing = TRUE)) +
  labs(
    x = expression("Mean absolute error of " * italic(J)[IE] * " vs exact BEDTools  (log)"),
    y = "Speedup vs BEDTools  (log)",
    title = "Precision frontier: what accuracy costs",
    subtitle = sprintf(
      paste0("Up is faster; left is more accurate. The frontier runs from cheap-and-rough (top right) to precise-and-slow (bottom left).\n",
             "20 Maurano DHS files, %d ordered pairs, t=%s, subB=1.0 -- N=20 is below the N~64 crossover, so nothing here is far above 1x."),
      as.integer(unique(d$n_pairs)), paste(threads_used, collapse = "/"))) +
  coord_cartesian(xlim = c(xlo, xhi), ylim = c(ylo, yhi), expand = FALSE) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8.2, colour = "grey30", lineheight = 1.25))

dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
CairoPNG(out_png, width = 1500, height = 1150, res = 200)
print(fig)
invisible(dev.off())
cat("Wrote:", out_png, "\n")
