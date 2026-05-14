#!/usr/bin/env Rscript
# Mode C interpolation between Mode A and Mode B, aggregate view.
#
# For every Mode C CSV at p=21, compute per-pair distance from each
# endpoint:
#   d_A = | j_C - j_A |
#   d_B = | j_C - j_B |
# Then aggregate across pairs and overlay the two curves vs the swept
# parameter. Two figures:
#
#   figures/mode_c_subB_interpolation_agg.png   subB sweep at expA=0
#   figures/mode_c_expA_interpolation_agg.png   expA sweep at subB=1.0
#
# Crossover is the swept value at which Mode C is equidistant from A
# and B; left of it Mode C is A-like, right of it B-like.
#
# Usage:
#   ml gcc/9.3.0 r/4.3.0 libjpeg/9c
#   Rscript mode_c_interpolation.R

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(scales); library(Cairo)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) >= 1) {
  dirname(normalizePath(sub("--file=", "", file_arg[1])))
} else getwd()
abc_dir <- file.path(script_dir, "results", "raw_abc")
fig_dir <- file.path(script_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(path, plot, w = 8, h = 5, dpi = 150) {
  CairoPNG(filename = path, width = w, height = h, units = "in", res = dpi)
  on.exit(dev.off()); print(plot)
}

P <- 21

read_jacc <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    select(file1, file2, jaccard_similarity)
}

mode_a <- read_jacc(file.path(abc_dir, sprintf("hammock_hll_p%d_jaccA.csv", P))) %>%
  rename(j_A = jaccard_similarity)
mode_b <- read_jacc(file.path(abc_dir, sprintf("hammock_hll_p%d_jaccB.csv", P))) %>%
  rename(j_B = jaccard_similarity)

# Helper: build the d_A / d_B aggregate from a list of (param_value, csv) pairs.
build_agg <- function(rows) {
  df <- bind_rows(lapply(rows, function(r) {
    read_jacc(r$path) %>% mutate(param = r$value)
  })) %>% rename(j_C = jaccard_similarity) %>%
    filter(file1 < file2) %>%
    inner_join(mode_a, by = c("file1", "file2")) %>%
    inner_join(mode_b, by = c("file1", "file2")) %>%
    mutate(d_A = abs(j_C - j_A), d_B = abs(j_C - j_B))
  df %>%
    group_by(param) %>%
    summarise(
      n_pairs = n(),
      mean_dA = mean(d_A), q25_dA = quantile(d_A, .25), q75_dA = quantile(d_A, .75),
      mean_dB = mean(d_B), q25_dB = quantile(d_B, .25), q75_dB = quantile(d_B, .75),
      .groups = "drop"
    )
}

plot_agg <- function(agg, x_label, log_x, title_suffix, out_file) {
  agg_long <- bind_rows(
    agg %>% transmute(param, ref = "|Mode C - Mode A|",
                      mean = mean_dA, q25 = q25_dA, q75 = q75_dA),
    agg %>% transmute(param, ref = "|Mode C - Mode B|",
                      mean = mean_dB, q25 = q25_dB, q75 = q75_dB)
  )
  p <- ggplot(agg_long, aes(x = param, y = mean, colour = ref, fill = ref)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = c("|Mode C - Mode A|" = "darkred",
                                   "|Mode C - Mode B|" = "darkgreen")) +
    scale_fill_manual(values   = c("|Mode C - Mode A|" = "darkred",
                                   "|Mode C - Mode B|" = "darkgreen")) +
    labs(x = x_label, y = "mean |delta Jaccard| across pairs",
         colour = NULL, fill = NULL,
         title = sprintf("Mode C interpolation between A and B (p=%d, %s)",
                         P, title_suffix),
         subtitle = sprintf("aggregated over %d non-self pairs; ribbons = IQR",
                            agg$n_pairs[1])) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
  if (log_x) p <- p + scale_x_log10()
  save_png(file.path(fig_dir, out_file), p, w = 8, h = 5)
  cat("Wrote", out_file, "\n")
  print(agg)
}

# ── subB sweep (expA = 0) ───────────────────────────────────────────────────
subB_files <- list.files(abc_dir,
  pattern = sprintf("^hammock_hll_p%d_jaccC_B[0-9.]+\\.csv$", P),
  full.names = TRUE)
subB_rows <- lapply(subB_files, function(p) {
  list(value = as.numeric(sub(".*_B([0-9.]+)\\.csv$", "\\1", p)),
       path  = p)
})
agg_subB <- build_agg(subB_rows)
plot_agg(agg_subB, x_label = "subB (log)", log_x = TRUE,
         title_suffix = "expA = 0",
         out_file = "mode_c_subB_interpolation_agg.png")

# ── expA sweep (subB = 1.0) ─────────────────────────────────────────────────
expA_files <- list.files(abc_dir,
  pattern = sprintf("^hammock_hll_p%d_jaccC_expA[0-9.]+\\.csv$", P),
  full.names = TRUE)
expA_rows <- lapply(expA_files, function(p) {
  list(value = as.numeric(sub(".*_expA([0-9.]+)\\.csv$", "\\1", p)),
       path  = p)
})
if (length(expA_rows) > 0) {
  agg_expA <- build_agg(expA_rows)
  plot_agg(agg_expA, x_label = "expA", log_x = FALSE,
           title_suffix = "subB = 1.0",
           out_file = "mode_c_expA_interpolation_agg.png")
} else {
  cat("No expA CSVs found.\n")
}
