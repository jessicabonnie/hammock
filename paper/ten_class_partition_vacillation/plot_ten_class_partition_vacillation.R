#!/usr/bin/env Rscript

required_packages <- c("readr", "ggplot2", "Cairo")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "\nOn Rockfish, run: ml libjpeg/9c r/4.3.0",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(Cairo)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Run this script with Rscript.", call. = FALSE)
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

argv <- commandArgs(trailingOnly = TRUE)
out_png <- if (length(argv) >= 1) {
  normalizePath(argv[1], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "ten_class_partition_vacillation.png")
}
source_tsv <- file.path(script_dir, "ten_class_partition_vacillation.tsv")

dat <- read_tsv(source_tsv, show_col_types = FALSE)
dat$left_label <- gsub("\\\\n", "\n", dat$left_label)
dat$right_label <- gsub("\\\\n", "\n", dat$right_label)

family_fill <- c("Heart" = "#D9E6F3", "Muscle" = "#FDE0BE")
family_border <- c("Heart" = "#4C78A8", "Muscle" = "#F58518")

left_boxes <- data.frame(
  x = 0.235,
  xmin = 0.05,
  xmax = 0.42,
  y = dat$left_y,
  ymin = dat$left_y - 0.07,
  ymax = dat$left_y + 0.07,
  label = dat$left_label,
  family = dat$family,
  stringsAsFactors = FALSE
)
right_boxes <- data.frame(
  x = 0.765,
  xmin = 0.58,
  xmax = 0.95,
  y = dat$right_y,
  ymin = c(0.69, 0.27, 0.27),
  ymax = c(0.83, 0.63, 0.63),
  label = dat$right_label,
  family = dat$family,
  stringsAsFactors = FALSE
)
right_boxes <- unique(right_boxes)

arrows <- data.frame(
  x = 0.44,
  xend = 0.56,
  y = dat$left_y,
  yend = dat$right_y,
  stringsAsFactors = FALSE
)

p <- ggplot() +
  geom_rect(
    data = rbind(left_boxes, right_boxes),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = family, color = family),
    linewidth = 0.45
  ) +
  geom_segment(
    data = arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.08, "in")),
    linewidth = 0.35,
    color = "grey35"
  ) +
  geom_text(
    data = rbind(left_boxes, right_boxes),
    aes(x = x, y = y, label = label),
    size = 2.7,
    lineheight = 0.9
  ) +
  annotate("text", x = 0.235, y = 0.93, label = "Selected-feature partition",
           fontface = "bold", size = 3.2) +
  annotate("text", x = 0.765, y = 0.93, label = "Alternate p=18 seed partition",
           fontface = "bold", size = 3.2) +
  scale_fill_manual(values = family_fill, guide = "none") +
  scale_color_manual(values = family_border, guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.17, 0.98), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(10, 12, 10, 12))

dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
CairoPNG(out_png, width = 7.2, height = 3.8, units = "in", res = 300,
         bg = "white")
print(p)
dev.off()
message("Wrote: ", out_png)
