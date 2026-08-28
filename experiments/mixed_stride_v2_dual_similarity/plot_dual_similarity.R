#!/usr/bin/env Rscript

required_packages <- c("dplyr", "readr", "ggplot2", "scales", "Cairo")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(Cairo)
})

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 3) {
  stop("usage: plot_dual_similarity.R DUAL_CSV BEDTOOLS_CSV OUT_PNG", call. = FALSE)
}

dual <- read_csv(argv[1], show_col_types = FALSE)
bedtools <- read_csv(argv[2], show_col_types = FALSE)
out_png <- argv[3]
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

rates <- c(1, 0.1, 0.01, 0.001)
required <- c("similarity", "subB", "wall_median", "mae_vs_bedtools")
if (!all(required %in% names(dual))) {
  stop("dual input lacks required columns", call. = FALSE)
}
if (nrow(dual) != 8 || !all(rates %in% dual$subB)) {
  stop("expected two similarity rows at each of four subB rates", call. = FALSE)
}
if (!all(c("rep", "run_id", "wall_time") %in% names(bedtools))) {
  stop("BEDTools input lacks timing columns", call. = FALSE)
}

bt_wall <- bedtools %>% distinct(rep, run_id, wall_time) %>% pull(wall_time) %>% median()
rate_label <- c("1" = "no\nsubsampling", "0.1" = "0.1\nsubsample",
                "0.01" = "0.01\nsubsample", "0.001" = "0.001\nsubsample")

bars <- dual %>%
  mutate(
    condition = factor(rate_label[as.character(subB)], levels = unname(rate_label)),
    similarity = factor(
      similarity,
      levels = c("register-equality", "inclusion-exclusion"),
      labels = c("register-equality", "inclusion-exclusion (+IE)")
    ),
    label = sprintf("%.2f s\nMAE %.1e", wall_median, mae_vs_bedtools)
  )

COL_RE <- "#D28B35"
COL_IE <- "#007C83"
COL_BT <- "#46515C"
COL_GRID <- "#D9DEE3"
COL_TEXT <- "#20262D"
dodge <- position_dodge2(width = 0.9, padding = 0.18, preserve = "single")

plot <- ggplot(bars, aes(x = condition, y = wall_median, fill = similarity)) +
  geom_hline(yintercept = bt_wall, color = COL_BT, linewidth = 0.65,
             linetype = "22") +
  geom_col(position = dodge, width = 0.82) +
  geom_text(aes(label = label), position = dodge, vjust = -0.22,
            size = 3.65, lineheight = 0.95, color = COL_TEXT,
            show.legend = FALSE) +
  annotate("label", x = 4.48, y = bt_wall,
           label = sprintf("BEDTools %.2f s", bt_wall), hjust = 1,
           vjust = -0.35, size = 3.7, color = COL_BT,
           fill = alpha("white", 0.9), linewidth = 0) +
  scale_fill_manual(values = c("register-equality" = COL_RE,
                               "inclusion-exclusion (+IE)" = COL_IE)) +
  scale_y_continuous(labels = label_number(accuracy = 1),
                     expand = expansion(mult = c(0, 0.24))) +
  scale_x_discrete(expand = expansion(add = c(0.55, 0.7))) +
  labs(x = NULL, y = "Wall time, 20-file corpus (s)", fill = NULL) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(
    text = element_text(color = COL_TEXT),
    axis.text = element_text(color = COL_TEXT),
    axis.text.x = element_text(size = 13, lineheight = 0.95),
    axis.title.y = element_text(size = 14),
    axis.line = element_line(color = "#6B747D", linewidth = 0.4),
    axis.ticks = element_line(color = "#6B747D", linewidth = 0.4),
    panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.4),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 13),
    legend.margin = margin(t = 8),
    plot.margin = margin(28, 18, 8, 10)
  )

CairoPNG(filename = out_png, width = 12.5, height = 6.5,
         units = "in", res = 300, bg = "white")
print(plot)
dev.off()
message("Wrote: ", out_png)
