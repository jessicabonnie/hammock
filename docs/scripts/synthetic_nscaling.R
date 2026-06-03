#!/usr/bin/env Rscript
# Paper Fig 2 — synthetic N-scaling: hammock vs bedtools wall time.
# Source CSV: docs/data/cpp_vs_bedtools_t16_20260512_160412.csv
# Output:     docs/figures/synthetic_nscaling.png
#
# Usage (from any cwd):
#   ml gcc/9.3.0 r/4.3.0 libjpeg/9c
#   Rscript docs/scripts/synthetic_nscaling.R
#
# Optional: pass an alternate CSV / output path:
#   Rscript docs/scripts/synthetic_nscaling.R <files-csv> <out-png>

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(readr)
  library(Cairo)
})

# ---------- resolve script-relative paths ----------

SCRIPT_DIR <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=^--file=).*", args, perl = TRUE))
  if (length(m) == 0) getwd() else dirname(normalizePath(m[1]))
}, error = function(e) getwd())
DOCS_DIR <- dirname(SCRIPT_DIR)

argv <- commandArgs(trailingOnly = TRUE)
files_csv <- if (length(argv) >= 1) {
  argv[1]
} else {
  file.path(DOCS_DIR, "data", "cpp_vs_bedtools_t16_20260512_160412.csv")
}
out_png <- if (length(argv) >= 2) {
  argv[2]
} else {
  file.path(DOCS_DIR, "figures", "synthetic_nscaling.png")
}

# ---------- load ----------

df <- read_csv(files_csv, show_col_types = FALSE)
if (!"sub_b" %in% names(df)) df$sub_b <- ifelse(df$tool == "bedtools", NA_real_, 1.0)
df$sub_b <- suppressWarnings(as.numeric(df$sub_b))

bt <- df %>% filter(tool == "bedtools") %>%
  select(num_files, num_threads, precision,
         bt_wall = mean_wall_time, bt_wall_err = std_wall_time,
         sort_wall = mean_sort_time) %>%
  arrange(num_files)
hm <- df %>% filter(grepl("^hammock_cpp_B", tool)) %>%
  mutate(sub_b = ifelse(is.na(sub_b), 1.0, sub_b)) %>%
  filter(abs(sub_b - 1.0) < 1e-9) %>%
  select(num_files,
         hm_wall_err = std_wall_time,
         hm_sketch = mean_sketch_creation_time,
         hm_compare = mean_comparison_time) %>%
  arrange(num_files)

t <- bt$num_threads[1]; p <- bt$precision[1]

# ---------- plot ----------

bedtools_series <- sprintf("bedtools wall (parallel-wrapped, t=%d)", t)
sort_series     <- "bedtools pre-sort (excluded from bedtools wall)"
compare_series  <- "hammock pairwise compare"
sketch_series   <- sprintf("hammock sketching (t=%d, p=%d)", t, p)

bt_rows   <- bt %>% transmute(num_files, y = bt_wall, yerr = bt_wall_err,
                              series = bedtools_series)
sort_rows <- bt %>% transmute(num_files, y = sort_wall, yerr = NA_real_,
                              series = sort_series)
sketch_rows  <- hm %>% transmute(num_files, y = hm_sketch, yerr = hm_wall_err,
                                 series = sketch_series)
compare_rows <- hm %>% transmute(num_files, y = pmax(hm_compare, 1e-4),
                                 yerr = NA_real_, series = compare_series)
series_df <- bind_rows(bt_rows, sketch_rows, compare_rows, sort_rows)

series_levels <- c(bedtools_series, sketch_series, compare_series, sort_series)
series_df$series <- factor(series_df$series, levels = series_levels)

color_map <- setNames(c("#1f77b4", "#ff7f0e", "#d62728", "#7f7f7f"), series_levels)
lty_map   <- setNames(c("solid", "solid", "solid", "dotted"), series_levels)
shape_map <- setNames(c(16, 15, 18, 4), series_levels)

g <- ggplot(series_df, aes(x = num_files, y = y,
                           color = series, linetype = series, shape = series)) +
  geom_errorbar(data = subset(series_df, !is.na(yerr)),
                aes(ymin = pmax(y - yerr, 1e-6), ymax = y + yerr),
                width = 0.05, alpha = 0.7, show.legend = FALSE) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  scale_color_manual(values = color_map, name = NULL) +
  scale_linetype_manual(values = lty_map, name = NULL) +
  scale_shape_manual(values = shape_map, name = NULL) +
  scale_x_continuous(
    trans = log2_trans(),
    breaks = bt$num_files,
    labels = bt$num_files,
    sec.axis = sec_axis(~ . ^ 2, name = "Pairwise compares (N²)",
                        breaks = bt$num_files ^ 2,
                        labels = function(x) format(x, big.mark = ",", scientific = FALSE))
  ) +
  scale_y_continuous(trans = log10_trans(), labels = label_log()) +
  labs(
    x = "Number of files (N) — pairwise compares = N²",
    y = "Wall time (s)",
    title = "Hammock vs bedtools: where the time goes",
    subtitle = sprintf("p=%d, t=%d, 10 000 intervals/file, 3 runs/config", p, t)
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
    legend.key = element_rect(fill = NA),
    panel.grid.minor = element_line(linewidth = 0.2, color = "grey90"),
    plot.title = element_text(face = "bold")
  )

dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
CairoPNG(out_png, width = 9.5, height = 6, units = "in", res = 200, bg = "white")
print(g); dev.off()
cat(sprintf("Wrote %s\n", out_png))
