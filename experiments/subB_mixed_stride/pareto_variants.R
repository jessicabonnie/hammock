#!/usr/bin/env Rscript
# Four variants of the headline_maurano_pareto plot, dumped into a single
# multi-page PDF for side-by-side comparison.
#
# Variant A: original — geom_path in subB-sweep order (current behaviour)
# Variant B: geom_path re-sorted by MAE so the line traces the Pareto frontier
# Variant C: no path at all — just labelled points
# Variant D: variant C + per-replicate IQR error bars
#
# Output: figures/pareto_variants.pdf (4 pages)
#
# Usage:
#   ml gcc/9.3.0 r/4.3.0 libjpeg/9c
#   Rscript pareto_variants.R

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(scales); library(Cairo)
})

script_dir <- dirname(normalizePath(sub("--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
results_dir <- file.path(script_dir, "results")
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

pick_latest <- function(pattern) {
  files <- sort(list.files(results_dir, pattern = pattern,
                           full.names = TRUE), decreasing = TRUE)
  if (length(files) == 0) NULL else files[1]
}

method_levels <- c("hash-threshold", "mixed-stride", "single-hash")
method_pal <- c("hash-threshold" = "#1b9e77",
                "mixed-stride"   = "#d95f02",
                "single-hash"    = "#7570b3")
bedtools_col <- "#333333"

# ---------- data ----------
maur_csv <- pick_latest("^sweep_maurano_[0-9].*\\.csv$")
bt_csv   <- pick_latest("^sweep_maurano_bedtools_.*\\.csv$")
if (is.null(maur_csv)) stop("missing sweep_maurano_*.csv")
if (is.null(bt_csv))   stop("missing sweep_maurano_bedtools_*.csv")

maur <- read_csv(maur_csv, show_col_types = FALSE) %>%
  mutate(method = factor(method, levels = method_levels))
bt   <- read_csv(bt_csv,   show_col_types = FALSE)

bt_wall <- bt %>% distinct(rep, run_id, wall_time) %>%
  summarise(m = median(wall_time)) %>% pull(m)

# ---------- per-(method, subB, rep) MAE + wall_time ----------
gt <- maur %>% filter(subB == 1.0) %>%
  group_by(size_class, file_a, file_b) %>%
  summarise(j_truth = mean(jaccard), .groups = "drop")

per_rep <- maur %>%
  inner_join(gt, by = c("size_class", "file_a", "file_b")) %>%
  mutate(abs_err = abs(jaccard - j_truth)) %>%
  group_by(method, subB, rep, run_id) %>%
  summarise(mae = mean(abs_err),
            wall_time = first(wall_time),
            .groups = "drop") %>%
  mutate(speedup = bt_wall / wall_time)

# Cell summary (median + IQR across replicates).
cell <- per_rep %>%
  group_by(method, subB) %>%
  summarise(
    mae_med  = median(mae),
    mae_lo   = quantile(mae, 0.25),
    mae_hi   = quantile(mae, 0.75),
    sp_med   = median(speedup),
    sp_lo    = quantile(speedup, 0.25),
    sp_hi    = quantile(speedup, 0.75),
    n_rep    = n(),
    .groups = "drop"
  ) %>% filter(subB < 1.0)

cat("Per-cell summary (median MAE / median speedup, IQR):\n")
print(cell %>% arrange(method, subB))

rec <- cell %>% filter(method == "mixed-stride", subB == 0.1)
ms_labels <- cell %>% filter(method == "mixed-stride")

theme_pub <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 14),
        plot.subtitle    = element_text(size = 10.5, color = "gray30"),
        plot.caption     = element_text(size = 9, color = "gray40", hjust = 0),
        legend.position  = "top")

axis_layers <- list(
  geom_hline(yintercept = 1, linetype = "dashed", color = bedtools_col,
             linewidth = 0.35),
  annotate("text", x = max(cell$mae_med), y = 1.02,
           label = "bedtools  →  1.00×", hjust = 1, vjust = 0,
           color = bedtools_col, size = 3.4, fontface = "italic"),
  scale_color_manual(values = method_pal, name = NULL),
  scale_fill_manual(values = method_pal, name = NULL),
  scale_x_log10(labels = label_number(accuracy = 0.0001)),
  scale_y_continuous(breaks = seq(0.5, 3.5, 0.5),
                     limits = c(0.55, 3.4)),
  labs(x = "Mean absolute Jaccard error vs subB=1.0 (log scale)",
       y = "Speedup vs bedtools"),
  theme_pub
)

label_layer <- list(
  geom_text(data = ms_labels, aes(label = paste0("subB=", subB)),
            color = "gray20", size = 3.0, fontface = "bold",
            nudge_y = 0.13, show.legend = FALSE),
  geom_point(data = rec, aes(x = mae_med, y = sp_med),
             shape = 21, size = 7, stroke = 0.8,
             color = method_pal[["mixed-stride"]], fill = NA),
  annotate("text", x = rec$mae_med * 1.06, y = rec$sp_med - 0.13,
           label = "recommended", hjust = 0, vjust = 1,
           fontface = "italic", size = 3.4,
           color = method_pal[["mixed-stride"]])
)

# ── Variant A: original (geom_path in subB sweep order) ────────────────────
cell_A <- cell %>% arrange(method, desc(subB))   # subB descending = sweep order
pA <- ggplot(cell_A, aes(x = mae_med, y = sp_med, color = method)) +
  geom_path(aes(group = method), linewidth = 0.5, alpha = 0.9) +
  geom_point(size = 3.4) +
  label_layer + axis_layers +
  labs(title = "Variant A: path in subB-sweep order (original)",
       subtitle = "Lines connect points in subB sweep order (0.5 → 0.01); they zigzag where MAE isn't monotone in subB.")

# ── Variant B: geom_path re-sorted by MAE ──────────────────────────────────
cell_B <- cell %>% arrange(method, mae_med)
pB <- ggplot(cell_B, aes(x = mae_med, y = sp_med, color = method)) +
  geom_path(aes(group = method), linewidth = 0.5, alpha = 0.9) +
  geom_point(size = 3.4) +
  label_layer + axis_layers +
  labs(title = "Variant B: path re-sorted by MAE",
       subtitle = "Each line is monotone in MAE — but you lose the 'as subB decreases' reading.")

# ── Variant C: no path, just points + labels ───────────────────────────────
pC <- ggplot(cell, aes(x = mae_med, y = sp_med, color = method)) +
  geom_point(size = 3.4) +
  label_layer + axis_layers +
  labs(title = "Variant C: points only, no connecting path",
       subtitle = "No false 'curve' implication. Each marker is one (method × subB) cell.")

# ── Variant D: points + per-rep IQR error bars ─────────────────────────────
pD <- ggplot(cell,
             aes(x = mae_med, y = sp_med,
                 color = method, fill = method)) +
  geom_errorbar(aes(ymin = sp_lo, ymax = sp_hi),
                width = 0, linewidth = 0.4, alpha = 0.7) +
  geom_errorbarh(aes(xmin = mae_lo, xmax = mae_hi),
                 height = 0, linewidth = 0.4, alpha = 0.7) +
  geom_point(size = 3.4) +
  label_layer + axis_layers +
  labs(title = "Variant D: points + per-replicate IQR error bars",
       subtitle = sprintf("Bars are the 25–75th percentile across %d replicates per cell. Width tiny ⇒ zigzag is structural, not noise.",
                          unique(cell$n_rep)[1]))

# ---------- multi-page PDF ----------
out_pdf <- file.path(figures_dir, "pareto_variants.pdf")
CairoPDF(out_pdf, width = 10, height = 6.2, onefile = TRUE)
for (p in list(pA, pB, pC, pD)) print(p)
dev.off()
cat("Wrote", out_pdf, "\n")
