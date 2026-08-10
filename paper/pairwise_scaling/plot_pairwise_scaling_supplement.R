# Figure 3 SUPPLEMENT — full register-equality vs +IE comparison, for the
# supplement / advisor review, not the main-text Figure 3. Forked 2026-08-10
# from plot_pairwise_scaling.R at the point where that script still carried
# both hammock variants (register-equality and +IE) in both panels; the
# main-text script has since been simplified to show +IE only, in place of
# register-equality, per CLAUDE.md divergence #2 (register-equality carries a
# chance-agreement floor and is not on bedtools' scale).
#
# This script exists to show the full picture: how much worse
# register-equality actually is when it is judged on the SAME footing as +IE
# -- MAE against exact BEDTools truth, not drift from hammock's own
# unsubsampled baseline (which is what the pre-2026-08-10 dual-bar Panel B
# reported and is a much smaller, more flattering number for a reason that
# has nothing to do with accuracy: subB=1.0 *is* that baseline, so its
# "drift" is a tautological zero). Register-equality's vs-bedtools MAE here
# is ~0.137 -- about 100x the +IE MAE (~0.0012-0.0016) -- which is the
# chance-agreement floor c (~0.17 per CLAUDE.md) dominating the estimate, not
# a subsampling effect. See docs/data/maurano_subB_re_vs_bedtools.csv for the
# source numbers and their provenance.
#
# Panel A was rebuilt in Aug 2026; see docs/figure3-panel-a-rebuild.md for the
# evidence. Three things there are load-bearing and easy to undo by accident:
#
#   1. The secondary axis counts N^2. The harness builds two DISJOINT sets of N
#      files and both tools run the full cross product, so N=512 is 262,144
#      pairs, not 130,816. N(N-1)/2 is Panel B's denominator, not Panel A's.
#   2. There is no pmax() floor on comparison_time, and its absence is checked.
#      hammock-cpp reported integer milliseconds before v0.7.0, so a zero here
#      means the CSV predates the microsecond timers -- fail, do not floor.
#   3. Breaks are pinned to decades. Left to itself log_breaks() picks breaks
#      100x apart over this range and the axis silently reads two decades per
#      gridline.

required_packages <- c("dplyr", "readr", "ggplot2", "scales", "patchwork", "Cairo")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(Cairo)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
  stop("Could not determine the script path. Run with Rscript.", call. = FALSE)
}
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

# p=18 (the CLI default), job 29652408, node c595, one exclusive allocation,
# 2026-08-09. Replaces cpp_vs_bedtools_t16_20260804_172242.csv, which was p=14 --
# a precision used nowhere else in the paper -- and whose bedtools leg ran under
# the three-process-per-pair bedtools.sh with its achieved parallelism
# unrecorded. Every bedtools row here carries mean_bedtools_parallel_eff; it
# reads 0.10 at N >= 64, i.e. the baseline converted ~1.6 of its 16 cores into
# throughput. Quote that alongside any speedup taken from this figure.
synthetic_csv <- file.path(
  repo_root, "docs", "data", "cpp_vs_bedtools_t16_p18.csv"
)
maurano_ie_summary_csv <- file.path(
  repo_root, "docs", "data", "maurano_subB_ie_summary.csv"
)
maurano_summary_csv <- file.path(
  repo_root, "docs", "data", "maurano_subB_summary.csv"
)
maurano_bedtools_csv <- file.path(
  repo_root, "docs", "data", "maurano_bedtools.csv"
)
# Register-equality MAE against exact BEDTools -- same yardstick as the +IE
# MAE already in maurano_ie_summary_csv, computed from the same underlying
# grid (experiments/subB_mixed_stride/results/sweep_maurano_ie_20260809_200658.csv,
# mixed-stride, p=18, t=8, rep=0 -- reps are byte-identical, see RESULTS.md
# "Harness validation") against docs/data/maurano_bedtools_ref.tsv.
maurano_re_bt_csv <- file.path(
  repo_root, "docs", "data", "maurano_subB_re_vs_bedtools.csv"
)

argv <- commandArgs(trailingOnly = TRUE)
out_png <- if (length(argv) >= 1) {
  normalizePath(argv[1], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "pairwise_scaling_supplement.png")
}
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

for (path in c(synthetic_csv, maurano_summary_csv, maurano_bedtools_csv,
               maurano_ie_summary_csv, maurano_re_bt_csv)) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
}

COL_BEDTOOLS <- "#46515C"
COL_HAMMOCK <- "#007C83"
COL_COMPARE <- "#D28B35"
COL_COMPARE_IE <- "#9B4D9B"
# Distinct from COL_COMPARE_IE on purpose: 4 independent reviewers (2026-08-09)
# flagged that "hammock total (+IE)" reusing COL_COMPARE_IE made it collide,
# visually and in the legend, with the unrelated "hammock sketch comparison
# (+IE)" series, on top of already nearly-overlapping COL_HAMMOCK in the line
# itself. This color is not reused anywhere else in the panel.
COL_TOTAL_IE <- "#C2185B"
COL_GRID <- "#D9DEE3"
COL_TEXT <- "#20262D"
base_family <- "sans"

theme_paper <- function(base_size = 10.5) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        face = "bold", size = rel(1.06), color = COL_TEXT,
        lineheight = 1.05, margin = margin(b = 6)
      ),
      plot.subtitle = element_blank(),
      axis.title = element_text(color = COL_TEXT),
      axis.text = element_text(color = COL_TEXT),
      axis.line = element_line(color = "#6B747D", linewidth = 0.35),
      axis.ticks = element_line(color = "#6B747D", linewidth = 0.35),
      panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.8), color = COL_TEXT),
      legend.key.width = grid::unit(1.1, "lines"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

# Panel A: synthetic scaling ----------------------------------------------------
synthetic_raw <- read_csv(synthetic_csv, show_col_types = FALSE)

synthetic_bt <- synthetic_raw %>%
  filter(tool == "bedtools") %>%
  transmute(
    num_files,
    wall_time = mean_wall_time,
    wall_sd = std_wall_time,
    threads = num_threads,
    precision
  )

synthetic_hm <- synthetic_raw %>%
  mutate(sub_b = suppressWarnings(as.numeric(sub_b))) %>%
  filter(grepl("^hammock_cpp_B", tool), abs(sub_b - 1) < 1e-9) %>%
  transmute(
    num_files,
    wall_time = mean_wall_time,
    wall_sd = std_wall_time,
    comparison_time = mean_comparison_time,
    threads = num_threads,
    precision
  )

# A different row, not a different column of synthetic_hm, so it cannot be
# transmuted out of it. left_join, not inner_join: against a CSV predating the
# --metrics-arm option this must yield NA and drop the series, not empty
# Panel A and trip the stop() below with a misleading message.
synthetic_ie <- synthetic_raw %>%
  filter(tool == "hammock_ie_B") %>%
  transmute(
    num_files,
    threads = num_threads,
    precision,
    comparison_time_ie = mean_comparison_time,
    # Total wall time for the +IE arm, so the plotted curve shows the cost a
    # reader actually pays end to end if they use the recommended
    # jaccard_similarity_ie column, not just its comparison-phase slice.
    # Added 2026-08-09, validated against the corrected Panel A rerun (job
    # 29671317, single-node, small-N cells fixed) 2026-08-10 -- no longer
    # provisional.
    wall_time_ie = mean_wall_time,
    wall_sd_ie = std_wall_time
  )

# Guards the left_join below: today hammock_ie_B is always written at exactly
# one sub_b (1.0) per (num_files, threads, precision), so the join can't fan
# out -- but nothing enforces that upstream, and a silent fan-out here would
# quietly multiply every OTHER series in `synthetic`, not just the IE one.
ie_key <- c("num_files", "threads", "precision")
if (anyDuplicated(synthetic_ie[ie_key]) != 0) {
  stop(
    "synthetic_ie has duplicate (num_files, threads, precision) keys -- a ",
    "left_join below would fan out every series in `synthetic`, not just ",
    "the +IE one. Check for a second hammock_ie_B row (e.g. a stray sub_b) ",
    "in ", basename(synthetic_csv), ".",
    call. = FALSE
  )
}

synthetic <- inner_join(
  synthetic_bt,
  synthetic_hm,
  by = c("num_files", "threads", "precision"),
  suffix = c("_bedtools", "_hammock")
) %>%
  left_join(synthetic_ie, by = c("num_files", "threads", "precision")) %>%
  arrange(num_files) %>%
  mutate(
    # Full cross product of two disjoint N-file sets, not within-set all-pairs.
    cross_pairs = num_files * num_files,
    speedup = wall_time_bedtools / wall_time_hammock,
    # NA, not an error, where wall_time_ie is NA (CSV predates --metrics-arm,
    # or this num_files has no IE row) -- the annotation below already
    # handles an NA speedup_ie by omitting the second number.
    speedup_ie = wall_time_bedtools / wall_time_ie
  )

if (nrow(synthetic) == 0) {
  stop("No matching unsubsampled synthetic benchmark rows were found.", call. = FALSE)
}

# Without the pmax() floor a zero lands at -Inf on the log axis and the point
# vanishes silently. A zero here means the CSV predates the microsecond timers
# (v0.7.0), where the comparison phase read 0.000000 for N <= 16.
if (any(synthetic$comparison_time <= 0, na.rm = TRUE)) {
  stop(
    "comparison_time <= 0: this CSV predates the microsecond timers in ",
    "hammock-cpp 0.7.0. Re-run the benchmark rather than restoring a floor.",
    call. = FALSE
  )
}
if (!"hammock_ie_B" %in% synthetic_raw$tool) {
  message(
    "note: no hammock_ie_B rows in ", basename(synthetic_csv),
    " -- the +IE series will be absent. Re-run with --metrics-arm."
  )
}

synthetic_long <- bind_rows(
  synthetic %>% transmute(
    num_files, value = wall_time_bedtools, error = wall_sd_bedtools,
    series = "BEDTools total"
  ),
  synthetic %>% transmute(
    num_files, value = wall_time_hammock, error = wall_sd_hammock,
    series = "hammock total"
  ),
  synthetic %>% transmute(
    num_files, value = comparison_time, error = NA_real_,
    series = "hammock sketch comparison"
  ),
  synthetic %>%
    filter(!is.na(comparison_time_ie)) %>%
    transmute(
      num_files, value = comparison_time_ie, error = NA_real_,
      series = "hammock sketch comparison (+IE)"
    ),
  # See the wall_time_ie comment above -- total wall time for the +IE arm,
  # not just its comparison-phase slice.
  synthetic %>%
    filter(!is.na(wall_time_ie)) %>%
    transmute(
      num_files, value = wall_time_ie, error = wall_sd_ie,
      series = "hammock total (+IE)"
    )
) %>%
  mutate(
    # A series absent from levels becomes NA and ggplot drops it with only a
    # warning, so the IE level is declared whether or not any row carries it.
    series = factor(
      series,
      levels = c("BEDTools total", "hammock total", "hammock total (+IE)",
                 "hammock sketch comparison", "hammock sketch comparison (+IE)")
    )
  )

n_breaks <- synthetic$num_files
pair_labels <- format(
  synthetic$cross_pairs,
  scientific = FALSE,
  big.mark = ",",
  trim = TRUE
)
largest <- synthetic %>% slice_max(num_files, n = 1, with_ties = FALSE)

# Both numbers in one label, not one number next to an arrow that visually
# terminates near BOTH the "hammock total" and "hammock total (+IE)" markers
# (they sit ~9% apart at N=512, nothing on a log axis spanning 6.6 decades --
# 2 independent reviewers flagged the single-number label as readable as
# belonging to either line). Stating both numbers in text removes the
# ambiguity without having to distort the arrow's real endpoints.
speedup_label <- if (!is.na(largest$speedup_ie)) {
  sprintf("%.1f× faster\n(%.1f× at +IE, the recommended column)",
          largest$speedup, largest$speedup_ie)
} else {
  sprintf("%.1f× faster", largest$speedup)
}

panel_a <- ggplot(
  synthetic_long,
  aes(x = num_files, y = value, color = series, linetype = series, shape = series)
) +
  geom_errorbar(
    data = filter(synthetic_long, !is.na(error)),
    aes(ymin = pmax(value - error, 1e-4), ymax = value + error),
    width = 0.05,
    linewidth = 0.35,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.25, stroke = 0.4) +
  annotate(
    "segment",
    x = largest$num_files,
    xend = largest$num_files,
    y = largest$wall_time_hammock * 1.18,
    yend = largest$wall_time_bedtools / 1.18,
    linewidth = 0.45,
    color = "#69737D",
    arrow = grid::arrow(ends = "both", length = grid::unit(0.08, "in"))
  ) +
  annotate(
    "label",
    x = largest$num_files / 1.18,
    y = sqrt(largest$wall_time_hammock * largest$wall_time_bedtools),
    label = speedup_label,
    size = 3,
    linewidth = 0,
    fill = alpha("white", 0.9),
    color = COL_TEXT,
    hjust = 1,
    lineheight = 0.95
  ) +
  scale_color_manual(values = c(
    "BEDTools total" = COL_BEDTOOLS,
    "hammock total" = COL_HAMMOCK,
    "hammock total (+IE)" = COL_TOTAL_IE,
    "hammock sketch comparison" = COL_COMPARE,
    "hammock sketch comparison (+IE)" = COL_COMPARE_IE
  ), drop = FALSE) +
  scale_linetype_manual(values = c(
    "BEDTools total" = "solid",
    "hammock total" = "solid",
    "hammock total (+IE)" = "dashed",
    "hammock sketch comparison" = "22",
    "hammock sketch comparison (+IE)" = "42"
  ), drop = FALSE) +
  scale_shape_manual(values = c(
    "BEDTools total" = 16,
    "hammock total" = 17,
    # Not 17 (shared with "hammock total", which this line sits almost on top
    # of at this panel's scale) -- 8 (asterisk) stays legible where the two
    # markers coincide, per the same reviewer note as COL_TOTAL_IE above.
    "hammock total (+IE)" = 8,
    "hammock sketch comparison" = 15,
    "hammock sketch comparison (+IE)" = 18
  ), drop = FALSE) +
  scale_x_continuous(
    trans = log2_trans(),
    breaks = n_breaks,
    labels = n_breaks,
    sec.axis = sec_axis(
      ~ ., breaks = n_breaks, labels = pair_labels,
      name = "File pairs compared, N²"
    )
  ) +
  scale_y_log10(
    # Breaks are pinned to decades. The data spans 6.6 decades, over which
    # scale_y_log10()'s default log_breaks() picks breaks 100x apart
    # (1e-4, 0.01, 1, 100) -- so consecutive gridlines are TWO decades and the
    # axis appears to jump 10 -> 1,000. On a log axis that reads as a decade
    # step unless you check the numbers.
    breaks = 10^(-4:3),
    minor_breaks = NULL,
    # Not label_number(accuracy = 0.1) either: the comparison series reaches
    # 2e-4 s, and a fixed 0.1 accuracy renders every sub-decisecond break as
    # "0.0". Three significant figures in fixed notation labels 0.0001 and
    # 1,000 correctly with one rule.
    labels = function(x) formatC(x, format = "fg", digits = 3, big.mark = ","),
    expand = expansion(mult = c(0.06, 0.15))
  ) +
  labs(
    title = "A  Sketch reuse increases the advantage as collections grow",
    x = "Number of BED files (N)",
    y = "Wall time (seconds, log scale)"
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_paper() +
  theme(
    axis.text.x.top = element_text(size = 7.7, color = "#59636D"),
    axis.title.x.top = element_text(size = 8.5, color = "#59636D"),
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.margin = margin(b = 2)
  )

# Panel B: Maurano real-data benchmark ----------------------------------------
maurano_summary <- read_csv(maurano_summary_csv, show_col_types = FALSE)
maurano_bedtools <- read_csv(maurano_bedtools_csv, show_col_types = FALSE)
maurano_ie_summary <- read_csv(maurano_ie_summary_csv, show_col_types = FALSE)
maurano_re_bt <- read_csv(maurano_re_bt_csv, show_col_types = FALSE)

required_summary <- c("method", "subB", "wall_median", "mae")
missing_summary <- setdiff(required_summary, names(maurano_summary))
if (length(missing_summary) > 0) {
  stop(
    "Maurano summary lacks columns: ",
    paste(missing_summary, collapse = ", "),
    call. = FALSE
  )
}
if (!all(c("rep", "run_id", "wall_time") %in% names(maurano_bedtools))) {
  stop("Maurano BEDTools input lacks rep, run_id, or wall_time.", call. = FALSE)
}
required_ie <- c("subB", "wall_median", "mae_ie_vs_bedtools")
missing_ie <- setdiff(required_ie, names(maurano_ie_summary))
if (length(missing_ie) > 0) {
  stop(
    "Maurano +IE summary lacks columns: ",
    paste(missing_ie, collapse = ", "),
    call. = FALSE
  )
}

bt_runs <- maurano_bedtools %>% distinct(rep, run_id, wall_time)
bt_wall <- median(bt_runs$wall_time, na.rm = TRUE)

mixed_stride <- maurano_summary %>%
  filter(method == "mixed-stride", subB %in% c(1, 0.1, 0.01)) %>%
  arrange(match(subB, c(1, 0.1, 0.01)))

if (nrow(mixed_stride) != 3) {
  stop("Expected mixed-stride rows for subB = 1, 0.1, and 0.01.", call. = FALSE)
}
if (nrow(maurano_ie_summary) != 3) {
  stop("Expected +IE summary rows for subB = 1, 0.1, and 0.01.", call. = FALSE)
}
if (nrow(maurano_re_bt) != 3) {
  stop("Expected register-equality vs-BEDTools rows for subB = 1, 0.1, and 0.01.",
       call. = FALSE)
}

mixed_stride <- mixed_stride %>%
  left_join(
    maurano_re_bt %>% transmute(subB = as.numeric(subB), mae_re_vs_bedtools),
    by = "subB"
  )
if (anyNA(mixed_stride$mae_re_vs_bedtools)) {
  stop(
    "maurano_subB_re_vs_bedtools.csv is missing a subB level present in ",
    basename(maurano_summary_csv), ".",
    call. = FALSE
  )
}

condition_of <- function(subb) case_when(
  subb == 1 ~ "no\nsubsampling",
  subb == 0.1 ~ "subB = 0.1",
  subb == 0.01 ~ "subB = 0.01"
)

# Two bars per condition, paired tightly (see position_dodge2's padding
# below): the existing register-equality (--no-metrics) timing, and the
# +IE (full metrics block) timing added 2026-08-10. They come from DIFFERENT
# hammock-cpp invocations at the same subB -- the +IE bars cost more because
# they compute the full containment/co-sketch block, not because of anything
# to do with subsampling itself (same story as Figure 3A's two total-wall
# curves). Both accuracy labels are now on the SAME scale, both vs exact
# BEDTools (changed 2026-08-10 -- see the file header). Register-equality's
# MAE here is ~100x the +IE MAE at every subB level; that gap is the
# chance-agreement floor (CLAUDE.md divergence #2), not a subsampling
# effect -- it barely moves across subB=1.0/0.1/0.01 (0.1379/0.1376/0.1370)
# while +IE's MAE does (0.00115/0.00156/0.00139). This is the whole point of
# this supplement panel: the main-text figure drops register-equality
# entirely rather than let a reader compare it to +IE at face value.
bars <- bind_rows(
  tibble(
    condition = "BEDTools",
    variant = "BEDTools",
    wall = bt_wall,
    label_extra = NA_character_
  ),
  mixed_stride %>% transmute(
    condition = condition_of(subB),
    variant = "register-equality",
    wall = wall_median,
    label_extra = sprintf(
      "MAE %s", formatC(mae_re_vs_bedtools, format = "e", digits = 1)
    )
  ),
  maurano_ie_summary %>% transmute(
    condition = condition_of(subB),
    variant = "+IE",
    wall = wall_median,
    label_extra = sprintf(
      "MAE %s", formatC(mae_ie_vs_bedtools, format = "e", digits = 1)
    )
  )
) %>%
  mutate(
    condition = factor(
      condition,
      levels = c("BEDTools", "no\nsubsampling", "subB = 0.1", "subB = 0.01")
    ),
    variant = factor(variant, levels = c("BEDTools", "register-equality", "+IE")),
    speedup = bt_wall / wall,
    # Never hardcode "faster" -- word it from the sign (see the original note
    # this replaced: the corrected BEDTools baseline flips the unsubsampled
    # bar to a genuine slowdown, and "0.90x faster" reads as a speedup to
    # anyone skimming).
    ratio_txt = ifelse(speedup >= 1,
                       sprintf("%.2f× faster", speedup),
                       sprintf("%.2f× slower", 1 / speedup)),
    # Two lines, not three: wall time + speedup share a line so the label's
    # vertical footprint fits inside a paired bar's height gap (as small as
    # 0.5s here) without colliding with the neighbor's label.
    label = case_when(
      variant == "BEDTools" ~ sprintf("%.1f s\nspeed reference", wall),
      TRUE ~ sprintf("%.1f s (%s)\n%s", wall, ratio_txt, label_extra)
    ),
    # Leans each paired bar's label AWAY from its partner instead of
    # centering it (the default) -- 2 independent reviewers (2026-08-10)
    # found the centered version bleeding text onto the neighboring bar at
    # this tight a padding. register-equality sits on the left of its pair,
    # so hjust=1 (right-aligned at its own bar's center) makes the text hang
    # leftward, into the gap toward the previous category; +IE sits on the
    # right, so hjust=0 hangs it rightward. Standalone BEDTools stays
    # centered (0.5).
    # Lean each paired bar's label AWAY from its partner instead of
    # centering it (the default) -- 2 independent reviewers (2026-08-10)
    # found the centered version bleeding text onto the neighboring bar at
    # this tight a padding. register-equality sits on the left of its pair,
    # so hjust=1 (right-aligned at its own bar's center) makes the text hang
    # leftward, into the gap toward the previous category; +IE sits on the
    # right, so hjust=0 hangs it rightward. Standalone BEDTools stays
    # centered (0.5). (A same-height-for-both-bars-in-a-pair variant was
    # tried and reverted -- it fixed the vertical overlap but broke the
    # y-axis's decimal breaks and pushed the rightmost pair's text off the
    # panel edge; per-own-bar height plus the wider padding/margin below
    # is the more robust fix.)
    label_hjust = case_when(
      variant == "BEDTools" ~ 0.5,
      variant == "register-equality" ~ 1.0,
      variant == "+IE" ~ 0.0
    )
  )

panel_b <- ggplot(bars, aes(x = condition, y = wall, fill = variant)) +
  geom_col(
    # padding is the gap WITHIN a pair, as a fraction of each bar's slot.
    # Widened stepwise (0.06 -> 0.11 -> 0.16) after each smaller value left
    # the hjust-leaned labels below still clipping the last pair or bleeding
    # onto the neighbor -- still "a little space" relative to the bar width,
    # not a wide gap. preserve="single" keeps every bar the same width
    # whether it's alone (BEDTools) or paired (the three hammock groups),
    # rather than stretching the lone BEDTools bar to fill its slot.
    position = position_dodge2(width = 0.82, padding = 0.16, preserve = "single")
  ) +
  geom_text(
    aes(label = label, hjust = label_hjust),
    position = position_dodge2(width = 0.82, padding = 0.16, preserve = "single"),
    vjust = -0.18,
    size = 1.95,
    lineheight = 0.92,
    color = COL_TEXT,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(
    "BEDTools" = COL_BEDTOOLS,
    "register-equality" = COL_HAMMOCK,
    "+IE" = COL_TOTAL_IE
  )) +
  scale_x_discrete(
    # Extra room past the last category specifically -- the rightmost
    # pair's +IE label leans further right (hjust=0) than default discrete
    # expansion (0.6 categories) left room for, and was clipping the panel
    # edge.
    expand = expansion(add = c(0.6, 1.0))
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 1),
    breaks = scales::breaks_pretty(n = 5),
    expand = expansion(mult = c(0, 0.30))
  ) +
  labs(
    title = "B  Subsampling further reduces runtime",
    x = NULL,
    # Deliberately not "per pairwise comparison": each bar is the median
    # TOTAL wall time to sketch all 20 files and run all 400 pairwise
    # comparisons in one benchmark invocation, not divided by pair count.
    y = "Wall time, 20-file corpus (s)"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_paper() +
  theme(
    axis.text.x = element_text(size = 8.8, lineheight = 0.95),
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.margin = margin(b = 2),
    legend.title = element_blank(),
    plot.title = element_text(size = 10.5, lineheight = 1.05)
  )

figure <- panel_a + panel_b +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(
    title = paste0(
      "Hammock expands feasible all-pairs comparison as interval collections grow  ",
      "(Supplement)"
    ),
    subtitle = paste0(
      "Panel B accuracy for BOTH hammock variants is MAE against exact BEDTools -- ",
      "register-equality's ~0.137 floor is chance agreement, not a subsampling effect."
    ),
    theme = theme(
      plot.title = element_text(
        family = base_family, face = "bold", size = 15.5, color = COL_TEXT,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        family = base_family, size = 10.3, color = "#56616C",
        margin = margin(b = 8)
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

CairoPNG(
  filename = out_png,
  width = 14.2,
  height = 6.5,
  units = "in",
  res = 300,
  bg = "white"
)
print(figure)
dev.off()

message("Wrote: ", out_png)
