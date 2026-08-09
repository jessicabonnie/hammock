# Supplementary Figure — How hammock and BEDTools use the cores they are given
#
# Companion to Figure 3 (plot_pairwise_scaling.R, same directory). Figure 3
# quotes a speedup at a fixed thread count; that number is only interpretable if
# the reader knows what each tool did with its cores. They differ structurally:
#
#   * hammock parallelises INSIDE one process (OpenMP over a shared sketch
#     pool) and reads each input file once.
#   * bedtools has no batch mode, so a pairwise workflow is N^2 independent
#     process launches under a GNU `parallel` wrapper, re-reading each file N
#     times. On the benchmark nodes process creation caps near 123 exec/s and
#     does not scale with cores, so this saturates early.
#
# The saturation is a property of the workflow, not a defect in bedtools, and
# not something the benchmark harness can fix. It is also not specific to
# bedtools: `md5sum` on node-local files measures 0.46x at 16-way on the same
# nodes, and `xargs -P16` hits the same ceiling as GNU parallel. Reporting it is
# what keeps "N x faster" honest -- without it, a speedup measured against a
# wrapper that failed to parallelise looks identical to one measured against a
# baseline that used its cores.
#
# THREE THINGS THAT ARE LOAD-BEARING HERE, mirroring the guardrails in
# plot_pairwise_scaling.R:
#
#   1. Efficiency is computed WITHIN a run, as speedup(t)/t against that same
#      run's own t=1 point. Absolutes on this cluster are not portable between
#      runs (see docs/seed-benchmark-methodology.md), so a t=1 reference taken
#      from a different job would silently import that error. The threads axis
#      is self-calibrating precisely because t=1 is one of its own points.
#   2. The ideal line is drawn, and the y axis is NOT clipped to it. Efficiency
#      above 1 is a real measurement (cache effects at low t), and clipping
#      would hide a run that had gone wrong.
#   3. Breaks are pinned to the measured thread counts. Left to itself the log2
#      axis invents intermediate breaks at thread counts that were never run.

required_packages <- c("dplyr", "readr", "ggplot2", "scales", "patchwork", "Cairo")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2)
  library(scales); library(patchwork); library(Cairo)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Run with Rscript.", call. = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

argv <- commandArgs(trailingOnly = TRUE)
threads_csv <- if (length(argv) >= 1) argv[1] else
  file.path(repo_root, "docs", "data", "sweep_threads_p18.csv")
out_png <- if (length(argv) >= 2) argv[2] else
  file.path(repo_root, "paper", "figures", "threading_supplement.png")

if (!file.exists(threads_csv)) {
  stop("Input not found: ", threads_csv,
       "\nGenerate it with experiments/bedtools_benchmark/sbatch_fig3_threading.sh",
       " and copy the resulting sweep_threads_*.csv here.", call. = FALSE)
}

raw <- read_csv(threads_csv, show_col_types = FALSE)

# One row per (tool, threads): median over replicates. Median, not mean, because
# a single slow replicate on a shared filesystem is the common failure mode here
# and it should not move the curve.
dat <- raw %>%
  filter(!is.na(wall_time), wall_time > 0) %>%
  mutate(tool_label = case_when(
    tool == "bedtools" ~ "BEDTools (GNU parallel, N² processes)",
    grepl("^hammock_cpp_B", tool) ~ "hammock (OpenMP, one process)",
    TRUE ~ tool)) %>%
  filter(tool_label != tool) %>%          # drop anything unrecognised, loudly below
  group_by(tool_label, threads) %>%
  summarise(wall = median(wall_time),
            wall_lo = min(wall_time), wall_hi = max(wall_time),
            cpu = median(cpu_time), n = n(), .groups = "drop")

if (nrow(dat) == 0) stop("No usable rows in ", threads_csv, call. = FALSE)

# Within-run t=1 reference. Fail rather than silently rebase on the smallest
# thread count present: without a real t=1 point, "efficiency" is not defined
# and a plot drawn anyway would be quietly meaningless.
base <- dat %>% filter(threads == 1) %>% select(tool_label, wall1 = wall)
missing_base <- setdiff(unique(dat$tool_label), base$tool_label)
if (length(missing_base) > 0) {
  stop("No t=1 row for: ", paste(missing_base, collapse = ", "),
       ". Efficiency is defined against this run's own single-thread point; ",
       "rerun the sweep with 1 in --thread-list.", call. = FALSE)
}

dat <- dat %>%
  left_join(base, by = "tool_label") %>%
  mutate(speedup = wall1 / wall, efficiency = speedup / threads)

thread_breaks <- sort(unique(dat$threads))
pal <- c("hammock (OpenMP, one process)" = "#ff7f0e",
         "BEDTools (GNU parallel, N² processes)" = "#000000")

common <- list(
  scale_x_continuous(trans = "log2", breaks = thread_breaks,
                     labels = thread_breaks),
  scale_colour_manual(values = pal, name = NULL),
  scale_fill_manual(values = pal, guide = "none"),
  labs(x = "Threads / concurrent jobs"),
  theme_bw(base_size = 11),
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())
)

# (A) wall time. min-max band across replicates, not a CI -- 3 replicates cannot
# support one, and drawing an SE ribbon would overstate what was measured.
pA <- ggplot(dat, aes(threads, wall, colour = tool_label, fill = tool_label)) +
  geom_ribbon(aes(ymin = wall_lo, ymax = wall_hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  # breaks_log, NOT trans_breaks("log10"). The latter emits fractional exponents
  # whenever the data span less than a decade -- the smoke-test sweep labelled
  # its axis 10^-0.15, 10^-0.2 -- which is unreadable and is the same class of
  # defect plot_pairwise_scaling.R pins decade breaks to avoid. breaks_log picks
  # 1-2-5 style values that stay meaningful at any span.
  scale_y_log10(breaks = breaks_log(n = 6),
                labels = label_number(drop0trailing = TRUE)) +
  common + labs(y = "Wall time (s)", title = "A. Wall time")

# (B) speedup vs this run's own t=1, with the ideal line.
pB <- ggplot(dat, aes(threads, speedup, colour = tool_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey55", linewidth = 0.4) +
  annotate("text", x = max(thread_breaks) * 0.62, y = max(thread_breaks) * 0.92,
           label = "ideal", colour = "grey40", size = 3, hjust = 0) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  scale_y_continuous(trans = "log2", breaks = thread_breaks) +
  common + labs(y = expression("Speedup vs " * italic(t) * " = 1"),
                title = "B. Parallel speedup")

# (C) efficiency. The headline: hammock converts cores into throughput; the
# bedtools workflow stops doing so almost immediately.
pC <- ggplot(dat, aes(threads, efficiency, colour = tool_label)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55",
             linewidth = 0.4) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, NA)) +
  common + labs(y = "Parallel efficiency (speedup / threads)",
                title = "C. Efficiency")

fig <- (pA | pB | pC) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
CairoPNG(out_png, width = 2400, height = 900, res = 200)
print(fig)
invisible(dev.off())
cat("Wrote:", out_png, "\n\n")

# Print the table too. The figure shows the shape; a caption needs the numbers,
# and this is what should be quoted rather than re-read off the axes.
dat %>%
  arrange(tool_label, threads) %>%
  mutate(wall = round(wall, 2), speedup = round(speedup, 2),
         efficiency = round(efficiency, 3),
         cpu_over_wall = round(cpu / wall, 1)) %>%
  select(tool_label, threads, n, wall, speedup, efficiency, cpu_over_wall) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
