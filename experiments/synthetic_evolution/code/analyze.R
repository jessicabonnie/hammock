#!/usr/bin/env Rscript
# =============================================================================
# code/analyze.R — synthetic_evolution analysis
# =============================================================================
#
# Loads hammock CSVs from results/hammock/, parses (rate, type, mode, expA,
# subB, generation) out of the file paths, and emits:
#   - figures/jaccard_by_generation.png       (headline A/B/C grid)
#   - figures/sensitivity_heatmap.png         (Jaccard@gen=50, A vs B vs C)
#   - figures/mode_c_expA_sweep.png           (Mode C expA effect, per
#                                              (rate, type), Jaccard vs gen)
#   - figures/mode_overlay_by_type.png        (per evolution type: A, B, and
#                                              Mode C across expA values
#                                              overlaid vs generation)
#   - figures/mode_c_expA_interpolation.png   (X=expA, Y=Jaccard@gen=50,
#                                              with Mode A and Mode B
#                                              reference lines)
#   - results/summary.csv                     (one row per
#                                              (rate, type, mode, expA, subB))
#
# Hammock CSV filename convention (set by Makefile):
#   synth_r{RATE}_{TYPE}_p{PRECISION}_jacc{MODE}[_B{SUBB}].csv
# Inside each CSV, bed2 paths end in:
#   synthetic_g{G}_r{RATE}_{TYPE}_{GEN3}.bed
#
# Usage:
#   Rscript code/analyze.R --hammock-dir results/hammock \
#                          --out-figures figures \
#                          --out-summary results/summary.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

# Cairo is required on Rockfish r/4.3.0 (no X11; ggsave default device fails).
# Fall back gracefully if Cairo isn't installed (e.g. dev box without it).
have_cairo <- requireNamespace("Cairo", quietly = TRUE)
if (have_cairo) library(Cairo)

# Minimal flag parser — kwarg-only: --foo bar. No positionals expected.
parse_flags <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- defaults
  i <- 1
  while (i <= length(args)) {
    flag <- sub("^--", "", args[[i]])
    if (i + 1 > length(args)) stop(sprintf("Flag --%s requires a value", flag))
    val <- args[[i + 1]]
    if (!flag %in% names(defaults))
      stop(sprintf("Unknown flag --%s", flag))
    if (is.integer(defaults[[flag]])) val <- as.integer(val)
    else if (is.numeric(defaults[[flag]])) val <- as.numeric(val)
    opt[[flag]] <- val
    i <- i + 2
  }
  opt
}

opt <- parse_flags(list(
  `hammock-dir` = "results/hammock",
  `out-figures` = "figures",
  `out-summary` = "results/summary.csv",
  `gen-snapshot` = 50L
))

dir.create(opt[["out-figures"]], showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt[["out-summary"]]), showWarnings = FALSE, recursive = TRUE)

save_png <- function(p, path, width = 10, height = 7, dpi = 120) {
  if (have_cairo) {
    Cairo::CairoPNG(path, width = width, height = height, units = "in",
                    res = dpi)
    print(p); dev.off()
  } else {
    ggsave(path, p, width = width, height = height, dpi = dpi)
  }
  message(sprintf("  wrote %s", path))
}

# ── Load all CSVs ────────────────────────────────────────────────────────────
csv_files <- list.files(opt[["hammock-dir"]], pattern = "^synth_r.*\\.csv$",
                        full.names = TRUE)
if (length(csv_files) == 0) {
  stop(sprintf("No CSVs in %s — did `make hammock` run?", opt[["hammock-dir"]]))
}
message(sprintf("Loading %d CSV file(s) from %s", length(csv_files),
                opt[["hammock-dir"]]))

# Parse metadata out of the filename. Filename (set by hammock outprefix.py):
#   synth_r{RATE}_{TYPE}_hll_p{PRECISION}_jacc{MODE}[_expA{E}][_B{SUBB}].csv
# - subB suffix omitted by hammock when subB == 1.0 → treat as subB=1.0.
# - expA suffix omitted when expA == 0 → treat as expA=0.0.
parse_meta <- function(path) {
  bn <- basename(path)
  m <- str_match(bn,
    "^synth_r([0-9.]+)_([ABC])_hll_p(\\d+)_jacc([ABC])(?:_expA([0-9.]+))?(?:_B([0-9.]+))?\\.csv$")
  if (any(is.na(m[1, 1])))
    stop(sprintf("Cannot parse hammock CSV name: %s", bn))
  mode <- m[1, 5]
  expA <- if (is.na(m[1, 6])) 0.0 else as.numeric(m[1, 6])
  subB <- if (mode == "C") {
    if (is.na(m[1, 7])) 1.0 else as.numeric(m[1, 7])
  } else NA_real_
  tibble(
    file        = path,
    rate        = as.numeric(m[1, 2]),
    evol_type   = m[1, 3],
    precision   = as.integer(m[1, 4]),
    mode        = mode,
    expA        = expA,
    subB        = subB
  )
}
meta <- map_dfr(csv_files, parse_meta)

# Read each CSV; tolerate either old column naming (jaccard, bed1/bed2) or new
# (jaccard_similarity, file1/file2).
load_csv <- function(path) {
  df <- suppressMessages(read_csv(path, show_col_types = FALSE))
  cols <- names(df)
  bed2col <- intersect(c("bed2", "file2"), cols)[1]
  jaccol  <- intersect(c("jaccard_similarity", "jaccard"), cols)[1]
  if (is.na(bed2col) || is.na(jaccol))
    stop(sprintf("Unexpected columns in %s: %s", basename(path),
                 paste(cols, collapse = ", ")))
  tibble(
    bed2    = basename(df[[bed2col]]),
    jaccard = df[[jaccol]]
  )
}

rows <- meta %>%
  mutate(data = map(file, load_csv)) %>%
  select(-file) %>%
  unnest(data) %>%
  mutate(
    generation = as.integer(str_extract(bed2, "(?<=_)\\d{3}(?=\\.bed$)"))
  ) %>%
  filter(!is.na(generation))

if (nrow(rows) == 0) stop("No rows parsed — check CSV bed2 column format.")

# Friendly names for evolution types (used everywhere downstream).
# A = Indel    : intervals appear or disappear
# B = Jitter   : interval positions drift by small amounts (±100bp)
# C = Combined : 50/50 mix of Indel and Jitter
EVOL_LABEL_MAP <- c(A = "Indel", B = "Jitter", C = "Combined")
EVOL_LABEL_LEVELS <- c("Indel", "Jitter", "Combined")
rows <- rows %>%
  mutate(evol_label = factor(EVOL_LABEL_MAP[evol_type],
                             levels = EVOL_LABEL_LEVELS))

# For the headline plot we show Mode A, Mode B, and the "default" Mode C
# (expA=0 at the experiment-wide subB pin, typically 0.25).
default_subB <- rows %>% filter(mode == "C", expA == 0.0) %>%
  pull(subB) %>% unique()
if (length(default_subB) != 1) {
  warning(sprintf("Multiple subB values for default Mode C: %s",
                  paste(default_subB, collapse = ", ")))
  default_subB <- default_subB[1]
}
headline <- rows %>%
  filter(mode %in% c("A", "B") |
         (mode == "C" & expA == 0.0 & subB == default_subB))

# ── Figure 1: Jaccard vs generation ──────────────────────────────────────────
message("Building Figure 1: Jaccard vs generation")
fig1 <- ggplot(headline, aes(x = generation, y = jaccard, color = mode)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  facet_grid(rows = vars(evol_label),
             cols = vars(sprintf("rate = %g", rate)),
             labeller = label_value,
             scales = "free_y") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c(A = "#d95f02", B = "#1b9e77", C = "#7570b3")) +
  labs(
    title    = "Hammock Jaccard vs. evolution generation",
    subtitle = sprintf("100k synthetic intervals; HLL p=%d; rows = Indel / Jitter / Combined evolution",
                       headline$precision[1]),
    x        = "Generation",
    y        = "Jaccard similarity to pre-evolution BED",
    color    = "Hammock mode"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

save_png(fig1, file.path(opt[["out-figures"]], "jaccard_by_generation.png"),
         width = 12, height = 8)

# ── Figure 2: sensitivity heatmap (Jaccard at snapshot generation) ──────────
message(sprintf("Building Figure 2: Jaccard at gen=%d", opt[["gen-snapshot"]]))
snap <- headline %>%
  filter(generation == opt[["gen-snapshot"]]) %>%
  group_by(rate, evol_label, mode) %>%
  summarise(jaccard = mean(jaccard), .groups = "drop")

# A lower Jaccard at fixed generation = mode is more *sensitive* to the change.
# Annotate the most-sensitive (lowest-Jaccard) mode per (rate, evol_label).
best <- snap %>%
  group_by(rate, evol_label) %>%
  slice_min(jaccard, n = 1) %>%
  ungroup()

fig2 <- ggplot(snap, aes(x = mode, y = factor(rate), fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 3.2, color = "black") +
  geom_point(data = best,
             shape = 8, color = "yellow", fill = "yellow",
             size = 4, stroke = 1.2) +
  facet_wrap(~ evol_label, nrow = 1) +
  scale_fill_gradient(low = "#fff7ec", high = "#7f0000", limits = c(0, 1)) +
  labs(
    title    = sprintf("Jaccard at generation %d (lower => more sensitive)",
                       opt[["gen-snapshot"]]),
    subtitle = "* marks the most-sensitive mode per (rate, evolution type)",
    x        = "Hammock mode",
    y        = "Per-generation change rate",
    fill     = "Jaccard"
  ) +
  theme_bw(base_size = 11)

save_png(fig2, file.path(opt[["out-figures"]], "sensitivity_heatmap.png"),
         width = 10, height = 4.5)

# ── Figure 3: Mode C expA sweep (subB held at the experiment pin) ───────────
message("Building Figure 3: Mode C expA sweep")
expa_df <- rows %>% filter(mode == "C", subB == default_subB)
if (nrow(expa_df) > 0 && length(unique(expa_df$expA)) > 1) {
  fig3 <- ggplot(expa_df, aes(x = generation, y = jaccard,
                              color = factor(sprintf("%.2f", expA)))) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.2) +
    facet_grid(rows = vars(evol_label),
               cols = vars(sprintf("rate = %g", rate)),
               scales = "free_y") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title    = "Mode C: effect of --expA on Jaccard decay",
      subtitle = sprintf("HLL p=%d; subB pinned at %.2f; expA=0 is the default",
                         expa_df$precision[1], default_subB),
      x        = "Generation",
      y        = "Jaccard similarity to pre-evolution BED",
      color    = "expA"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  save_png(fig3, file.path(opt[["out-figures"]], "mode_c_expA_sweep.png"),
           width = 12, height = 8)
} else {
  message("  skipping (no Mode C expA sweep data found)")
}

# ── Figure 4: 3 graphs per evolution type — one PNG per rate ────────────────
# Each PNG: 1 row × 3 cols (Type A | Type B | Type C). Each panel: X = gen,
# Y = jaccard, lines = Mode A, Mode B, Mode C default, Mode C expA = {0.5,
# 1, 2, 3}. The literal answer to "three graphs for three evolution types".
message("Building Figure 4 series: per-type comparison, one PNG per rate")
overlay <- rows %>%
  filter(mode %in% c("A", "B") |
         (mode == "C" & subB == default_subB)) %>%
  mutate(series = case_when(
    mode == "A" ~ "Mode A",
    mode == "B" ~ "Mode B",
    mode == "C" & expA == 0 ~ sprintf("Mode C (expA=0)"),
    mode == "C"             ~ sprintf("Mode C (expA=%.1f)", expA)
  ))
series_levels <- c(
  "Mode A", "Mode B",
  "Mode C (expA=0)",
  sprintf("Mode C (expA=%.1f)",
          sort(unique(overlay$expA[overlay$mode == "C" & overlay$expA > 0])))
)
overlay$series <- factor(overlay$series, levels = series_levels)

n_c_variants <- length(series_levels) - 2L
series_colors <- c(
  "Mode A" = "#d95f02",
  "Mode B" = "#1b9e77",
  setNames(colorRampPalette(c("#bcbddc", "#3f007d"))(n_c_variants),
           series_levels[-c(1, 2)])
)
series_ltypes <- c(
  "Mode A" = "solid",
  "Mode B" = "solid",
  setNames(rep("dashed", n_c_variants), series_levels[-c(1, 2)])
)

for (r in sort(unique(overlay$rate))) {
  sub <- overlay %>% filter(rate == r)
  fig <- ggplot(sub, aes(x = generation, y = jaccard,
                         color = series, linetype = series)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.2) +
    facet_wrap(~ evol_label, nrow = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = series_colors) +
    scale_linetype_manual(values = series_ltypes) +
    labs(
      title    = sprintf("Hammock Jaccard vs generation, by evolution type (rate = %g)", r),
      subtitle = sprintf("HLL p=%d; subB pinned at %.2f for all Mode C runs",
                         sub$precision[1], default_subB),
      x        = "Generation",
      y        = "Jaccard similarity to pre-evolution BED",
      color    = NULL, linetype = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9))
  path <- file.path(opt[["out-figures"]],
                    sprintf("mode_per_type_r%g.png", r))
  save_png(fig, path, width = 14, height = 5.5)
}

# ── Figure 5: expA-as-interpolator ──────────────────────────────────────────
# For each rate × evol_type, plot Jaccard@gen-snapshot vs expA. Add
# horizontal reference lines for Mode A and Mode B Jaccards at the same
# generation. Visually shows expA sliding Mode C between the two extremes.
message("Building Figure 5: expA interpolation")
modec_snap <- rows %>%
  filter(mode == "C", subB == default_subB,
         generation == opt[["gen-snapshot"]]) %>%
  group_by(rate, evol_label, expA) %>%
  summarise(jaccard = mean(jaccard), .groups = "drop")

ab_ref <- rows %>%
  filter(mode %in% c("A", "B"),
         generation == opt[["gen-snapshot"]]) %>%
  group_by(rate, evol_label, mode) %>%
  summarise(jaccard = mean(jaccard), .groups = "drop")

if (nrow(modec_snap) > 0) {
  fig5 <- ggplot(modec_snap, aes(x = expA, y = jaccard,
                                 color = evol_label)) +
    geom_hline(data = ab_ref,
               aes(yintercept = jaccard, color = evol_label, linetype = mode),
               linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    facet_wrap(~ sprintf("rate = %g", rate), nrow = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_brewer(palette = "Set1") +
    scale_linetype_manual(values = c("A" = "dashed", "B" = "dotted")) +
    labs(
      title    = sprintf("Mode C as interpolator: Jaccard at gen %d vs expA",
                         opt[["gen-snapshot"]]),
      subtitle = sprintf(
        "HLL p=%d; subB pinned at %.2f. Dashed lines = Mode A reference; dotted = Mode B reference.",
        rows$precision[1], default_subB),
      x        = "expA (0 = no interval boost; higher ⇒ more like Mode A)",
      y        = sprintf("Jaccard (vs. pre-evolution BED) at generation %d",
                         opt[["gen-snapshot"]]),
      color    = "Evolution",
      linetype = "Reference mode"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  save_png(fig5,
           file.path(opt[["out-figures"]], "mode_c_expA_interpolation.png"),
           width = 13, height = 6)
} else {
  message("  skipping (no Mode C data at snapshot generation)")
}

# ── Summary CSV ──────────────────────────────────────────────────────────────
message("Writing summary CSV")
summary_df <- rows %>%
  filter(generation == opt[["gen-snapshot"]]) %>%
  arrange(rate, evol_type, mode, expA, subB) %>%
  select(rate, evol_type, evol_label,
         mode, expA, subB, generation, jaccard)
write_csv(summary_df, opt[["out-summary"]])
message(sprintf("  wrote %s (%d rows)", opt[["out-summary"]], nrow(summary_df)))

message("Done.")
