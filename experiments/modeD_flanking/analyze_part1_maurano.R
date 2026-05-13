#!/usr/bin/env Rscript
# analyze_part1_maurano.R — Part 1 of modeD_flanking.
#
# Re-reads every Mode D CSV in maurano_dhs_validation/results/raw_d/ and
# compares both Jaccard columns (with vs without ends) to the bedtools
# pairwise Jaccard reference. Produces:
#
#   results/part1_flanking_vs_bedtools.csv
#     one row per (k, w, p, file_pair) with columns:
#       j_truth, j_no_ends, j_with_ends, err_no_ends, err_with_ends,
#       delta_err, phi_pair
#
#   results/part1_summary.csv
#     aggregated to one row per (k, w, p): mean Pearson r per column,
#     mean MAE per column, mean delta_err, mean phi
#
#   figures/maurano_delta_r_vs_w.png   — Pearson r of each column vs w, by k
#   figures/maurano_delta_mae_vs_phi.png — Δmae as a function of φ
#
# No new compute on hammock side. Run after maurano_dhs_validation has
# finished its sweep.
#
# Usage:
#   ml gcc/9.3.0 r/4.3.0 libjpeg/9c
#   Rscript analyze_part1_maurano.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
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

results_dir <- file.path(script_dir, "results")
figures_dir <- file.path(script_dir, "figures")
maurano    <- file.path(script_dir, "data", "maurano_link")
raw_d_dir  <- file.path(maurano, "results", "raw_d")
ref_path   <- file.path(maurano, "data", "maurano_bedtools_ref.tsv")
fastas_dir <- file.path(maurano, "data", "fastas")

stopifnot(dir.exists(raw_d_dir), file.exists(ref_path))
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(path, plot, width = 8, height = 5, dpi = 150) {
  CairoPNG(filename = path, width = width, height = height,
           units = "in", res = dpi)
  on.exit(dev.off())
  print(plot)
}

strip_ext <- function(x) sub("\\.(bed|fa|fasta|fna)$", "", x, ignore.case = TRUE)

# ── bedtools reference ──────────────────────────────────────────────────────
ref <- read_tsv(ref_path, show_col_types = FALSE) %>%
  transmute(stem1 = strip_ext(file1), stem2 = strip_ext(file2),
            j_truth = as.numeric(jaccard))

# ── per-file flanking budget φ_file = 2(k-1) * n_intervals / (total_len / w) ─
# φ_pair is the geometric mean of φ_file for the two files in the pair.
fa_stats <- function(fa_path) {
  # Quick scan: count records (>) and total residues. We don't need
  # exact accuracy; just length distribution scale.
  records <- 0L
  total_len <- 0L
  con <- file(fa_path, "r")
  on.exit(close(con))
  while (length(line <- readLines(con, n = 1L, warn = FALSE)) > 0) {
    if (startsWith(line, ">")) {
      records <- records + 1L
    } else {
      total_len <- total_len + nchar(line)
    }
  }
  list(records = records, total_len = total_len)
}

cat("Cataloguing FASTA length/interval stats from", fastas_dir, "...\n")
fa_files <- list.files(fastas_dir, pattern = "\\.fa$", full.names = TRUE)
fa_meta <- bind_rows(lapply(fa_files, function(f) {
  s <- fa_stats(f)
  tibble(stem = strip_ext(basename(f)),
         n_intervals = s$records,
         total_len   = s$total_len)
}))

phi_pair <- function(stem1, stem2, k, w, meta = fa_meta) {
  # Boundary k-mers per file ≈ 2 (k-1) n_intervals.
  # Interior minimizers per file ≈ total_len / w.
  # φ_file = boundary / interior, φ_pair = sqrt(φ1 * φ2).
  f1 <- meta[match(stem1, meta$stem), ]
  f2 <- meta[match(stem2, meta$stem), ]
  num1 <- 2 * (k - 1) * f1$n_intervals
  num2 <- 2 * (k - 1) * f2$n_intervals
  den1 <- pmax(1, f1$total_len / w)
  den2 <- pmax(1, f2$total_len / w)
  sqrt((num1 / den1) * (num2 / den2))
}

parse_d_name <- function(stem) {
  m <- regmatches(stem, regexec(
    "^hammock_mnmzr_p(\\d+)_jaccD_k(\\d+)_w(\\d+)$", stem, perl = TRUE))[[1]]
  if (length(m) == 0) return(NULL)
  tibble(precision = as.integer(m[2]), k = as.integer(m[3]), w = as.integer(m[4]))
}

# ── scan all Mode D CSVs ────────────────────────────────────────────────────
cat("Scanning", raw_d_dir, "...\n")
csvs <- list.files(raw_d_dir, pattern = "\\.csv$", full.names = TRUE)
cat("Found", length(csvs), "CSVs\n")

pair_rows <- list()
for (f in csvs) {
  meta <- parse_d_name(tools::file_path_sans_ext(basename(f)))
  if (is.null(meta)) next
  df <- read_csv(f, show_col_types = FALSE) %>%
    transmute(stem1 = strip_ext(file1), stem2 = strip_ext(file2),
              j_no_ends   = jaccard_similarity,
              j_with_ends = jaccard_similarity_with_ends)
  j <- inner_join(df, ref, by = c("stem1", "stem2"))
  if (nrow(j) == 0) next
  j$k <- meta$k; j$w <- meta$w; j$precision <- meta$precision
  j$phi_pair <- phi_pair(j$stem1, j$stem2, meta$k, meta$w)
  pair_rows[[length(pair_rows) + 1L]] <- j
}
pairs <- bind_rows(pair_rows) %>%
  mutate(err_no_ends   = j_no_ends   - j_truth,
         err_with_ends = j_with_ends - j_truth,
         delta_err     = abs(err_with_ends) - abs(err_no_ends))

write_csv(pairs, file.path(results_dir, "part1_flanking_vs_bedtools.csv"))
cat("Wrote part1_flanking_vs_bedtools.csv (", nrow(pairs), "rows)\n")

# ── aggregated summary per (k, w, p) ────────────────────────────────────────
agg <- pairs %>%
  group_by(k, w, precision) %>%
  summarise(
    n             = n(),
    r_no_ends     = cor(j_no_ends,   j_truth),
    r_with_ends   = cor(j_with_ends, j_truth),
    mae_no_ends   = mean(abs(err_no_ends)),
    mae_with_ends = mean(abs(err_with_ends)),
    mean_phi      = mean(phi_pair),
    mean_delta    = mean(delta_err),
    .groups = "drop"
  ) %>%
  mutate(delta_r   = r_no_ends   - r_with_ends,
         delta_mae = mae_with_ends - mae_no_ends)
write_csv(agg, file.path(results_dir, "part1_summary.csv"))
cat("Wrote part1_summary.csv (", nrow(agg), "rows)\n")

# ── plots ───────────────────────────────────────────────────────────────────
# 1. r vs w, per k, faceted by precision, lines coloured by column.
long_r <- agg %>%
  select(k, w, precision, r_no_ends, r_with_ends) %>%
  pivot_longer(c(r_no_ends, r_with_ends),
               names_to = "column", values_to = "pearson") %>%
  mutate(column = sub("^r_", "", column))
p1 <- ggplot(long_r,
             aes(x = w, y = pearson, colour = column,
                 group = interaction(k, column))) +
  geom_line() + geom_point(size = 1.5) +
  facet_grid(k ~ precision, labeller = label_both) +
  scale_x_log10() +
  labs(x = "w (log)", y = "Pearson r vs bedtools",
       title = "Mode D: with_ends vs no_ends across (k, w, p) — Maurano corpus",
       colour = "column") +
  theme_minimal(base_size = 10)
save_png(file.path(figures_dir, "maurano_delta_r_vs_w.png"), p1,
         width = 10, height = 8)

# 2. mean_delta_err vs phi (one point per (k, w, p) cell).
p2 <- ggplot(agg, aes(x = mean_phi, y = mean_delta,
                      colour = factor(precision))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.8) +
  scale_x_log10() +
  labs(x = expression("flanking fraction "*phi*" (log)"),
       y = "mean(|err with_ends|) - mean(|err no_ends|)",
       title = "Sign-of-flanking-gain vs φ on the Maurano corpus",
       subtitle = "y < 0 ⇒ with_ends is better; y > 0 ⇒ no_ends is better",
       colour = "precision") +
  theme_minimal(base_size = 11)
save_png(file.path(figures_dir, "maurano_delta_mae_vs_phi.png"), p2,
         width = 8, height = 5)

cat("Done.\n")
