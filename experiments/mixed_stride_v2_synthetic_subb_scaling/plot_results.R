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
if (length(argv) != 2) stop("usage: plot_results.R SUMMARY_CSV OUT_PNG", call. = FALSE)
data <- read_csv(argv[1], show_col_types = FALSE)
out_png <- argv[2]
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

rates <- c(1, 0.1, 0.01, 0.001)
rate_labels <- c("1", "0.1", "0.01", "0.001")
hm <- data %>%
  filter(tool == "hammock") %>%
  mutate(
    similarity = factor(similarity,
                        levels = c("register-equality", "inclusion-exclusion"),
                        labels = c("register-equality", "inclusion-exclusion (+IE)")),
    subB = factor(subB, levels = rates, labels = rate_labels)
  )
bt <- data %>% filter(tool == "bedtools")

hm_long <- bind_rows(
  hm %>% transmute(num_files, similarity, subB, measure = "Wall time (s)",
                   value = median_wall_time),
  hm %>% transmute(num_files, similarity, subB, measure = "Peak RSS (MiB)",
                   value = median_max_rss_mb),
  hm %>% transmute(num_files, similarity, subB, measure = "MAE vs BEDTools",
                   value = mean_mae_vs_bedtools)
) %>% mutate(measure = factor(measure,
                              levels = c("Wall time (s)", "Peak RSS (MiB)",
                                         "MAE vs BEDTools")))

bt_long <- bind_rows(
  bt %>% transmute(num_files, measure = "Wall time (s)", value = median_wall_time),
  bt %>% transmute(num_files, measure = "Peak RSS (MiB)", value = median_max_rss_mb)
) %>% mutate(measure = factor(measure, levels = levels(hm_long$measure)))

colors <- c("1" = "#007C83", "0.1" = "#D28B35",
            "0.01" = "#9B4D9B", "0.001" = "#C2185B")

plot <- ggplot(hm_long, aes(x = num_files, y = value, color = subB,
                            group = subB)) +
  geom_line(data = bt_long, aes(x = num_files, y = value),
            inherit.aes = FALSE, color = "#46515C", linewidth = 0.8,
            linetype = "22") +
  geom_point(data = bt_long, aes(x = num_files, y = value),
             inherit.aes = FALSE, color = "#46515C", size = 2.1) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.2) +
  scale_color_manual(values = colors) +
  scale_x_continuous(trans = log2_trans(), breaks = 2^(1:9), labels = 2^(1:9)) +
  scale_y_log10(labels = label_number(), expand = expansion(mult = c(0.08, 0.12))) +
  facet_grid(measure ~ similarity, scales = "free_y") +
  labs(x = "Number of BED files per set (N)", y = NULL,
       color = "subB",
       caption = "Dashed gray: exact BEDTools baseline (time and RSS only). p=18; 16 threads; 10,000 intervals/file; N² cross-product pairs.") +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    text = element_text(color = "#20262D"),
    axis.text = element_text(color = "#20262D"),
    panel.grid.major.y = element_line(color = "#D9DEE3", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#EEF1F3", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, size = 9),
    plot.margin = margin(8, 12, 8, 8)
  )

CairoPNG(filename = out_png, width = 13.5, height = 11,
         units = "in", res = 300, bg = "white")
print(plot)
dev.off()
message("Wrote: ", out_png)
