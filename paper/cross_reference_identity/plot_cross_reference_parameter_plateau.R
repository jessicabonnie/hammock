#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Cairo)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
  stop("Could not determine script path. Run with Rscript.", call. = FALSE)
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
repo_root <- normalizePath(file.path(script_dir, "..", ".."))

argv <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(argv) >= 1) normalizePath(argv[1]) else file.path(
  repo_root, "docs", "data", "exp_a_estimator_delta_expanded.csv"
)
output_png <- if (length(argv) >= 2) normalizePath(argv[2], mustWork = FALSE) else file.path(
  repo_root, "paper", "figures", "cross_reference_parameter_plateau.png"
)

required <- c(
  "peak_type", "k", "w", "median_xref_ie", "median_diff_ie"
)
raw <- read_csv(input_csv, show_col_types = FALSE)
missing <- setdiff(required, names(raw))
if (length(missing) > 0) {
  stop("Input lacks columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

text_color <- "#263238"

d <- raw %>%
  filter(peak_type == "broad") %>%
  transmute(
    k = as.numeric(k),
    w = as.numeric(w),
    same_tissue = median_xref_ie,
    different_tissue = median_diff_ie,
    optimum = k == 20 & w == 20
  ) %>%
  ungroup()

if (nrow(d) != 63 || n_distinct(d$k) != 10) {
  stop("Expected 63 broad-peak cells across 10 k values; found ",
       nrow(d), " cells across ", n_distinct(d$k), " k values", call. = FALSE)
}

# Retain the complete sweep in the source table, but display a regular subset:
# w=k, the available rounded-up multiple of 10, and available multiples of 25.
# Because the axis is discrete, every tick remains an actually tested value.
d <- d %>%
  filter(k < 50) %>%
  filter(
    w == k |
      w == ceiling(k / 10) * 10 |
      (w %% 25 == 0 & w >= ceiling(k / 25) * 25)
  )

points <- d %>%
  pivot_longer(
    cols = c(same_tissue, different_tissue),
    names_to = "comparison",
    values_to = "median_similarity"
  ) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c("same_tissue", "different_tissue"),
      labels = c("Same tissue across references", "Different tissues")
    )
  )

W_LEVELS <- sort(unique(d$w))

p <- ggplot(d, aes(x = factor(w, levels = W_LEVELS))) +
  # Train every free-x facet on its complete ordered w set before the optimum
  # and non-optimum rows are split across drawing layers. Otherwise k=20 is
  # first trained without its highlighted w=20 row and appends that tick last.
  geom_blank(aes(y = same_tissue)) +
  geom_segment(
    data = d %>% filter(!optimum),
    aes(
      y = different_tissue,
      yend = same_tissue,
      xend = factor(w, levels = W_LEVELS)
    ),
    color = "#AEB7BF", linewidth = 0.35, alpha = 0.78
  ) +
  geom_segment(
    data = d %>% filter(optimum),
    aes(
      y = different_tissue,
      yend = same_tissue,
      xend = factor(w, levels = W_LEVELS)
    ),
    color = "#D81B60", linewidth = 1.0, alpha = 1
  ) +
  geom_point(
    data = points,
    aes(y = median_similarity, color = comparison),
    size = 3.0, alpha = 1
  ) +
  facet_grid(. ~ k, scales = "free_x", space = "free_x", labeller = label_both) +
  scale_color_manual(
    values = c(
      "Same tissue across references" = "#00A6B2",
      "Different tissues" = "#65727E"
    ),
    name = NULL
  ) +
  scale_y_continuous(limits = c(0.30, 1.02), breaks = seq(0.4, 1.0, 0.2)) +
  labs(
    x = "Minimizer window size (w)",
    y = "Median Jaccard estimate"
  ) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title = element_text(size = 13, color = text_color),
    axis.text = element_text(size = 12, color = text_color),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = grid::unit(0.18, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold", color = text_color),
    legend.position = "bottom",
    legend.text = element_text(size = 11.5, color = text_color),
    plot.margin = margin(12, 16, 12, 12)
  )

dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
CairoPNG(output_png, width = 12.5, height = 4.6, units = "in", res = 300, bg = "white")
print(p)
dev.off()
message("Wrote: ", output_png)
