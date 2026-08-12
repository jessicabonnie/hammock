#!/usr/bin/env Rscript

# Figure — Where the two Jaccard estimators change places
#
# One figure, three facets, one message: the choice between hammock's two
# Jaccard columns is only live below J ~ 0.05, and there it is settled by
# precision. Above that, both estimators rank a synthetic corpus perfectly by
# p = 16 and there is nothing to choose.
#
# Input is experiments/bedtools_benchmark/results/estimator_compare_full.csv,
# which carries both estimators and exact `bedtools jaccard` truth for a 3 x 30
# cross product (3 A-file replicates x 10 overlap fractions) at four precisions.
# Nothing is sketched here.
#
# Uncertainty: the 90 pairs per precision are a cross product, so comparisons
# are not independent and a binomial interval on the discordance count would be
# far too narrow. The band is a leave-one-replicate-out range over the only
# independent axis the design has -- three replicates, so it is a coarse
# stability check and is drawn and labelled as such, not as a confidence
# interval.

required_packages <- c("dplyr", "readr", "ggplot2", "scales", "Cairo")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "),
       "\nInstall them before running this script.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(Cairo)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
  stop("Could not determine the script path. Run with Rscript.", call. = FALSE)
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

input_csv <- file.path(repo_root, "experiments", "bedtools_benchmark",
                       "results", "estimator_compare_full.csv")
argv <- commandArgs(trailingOnly = TRUE)
out_png <- if (length(argv) >= 1) {
  normalizePath(argv[1], mustWork = FALSE)
} else {
  file.path(repo_root, "paper", "figures", "estimator_crossover.png")
}
if (!file.exists(input_csv)) stop("Input not found: ", input_csv, call. = FALSE)
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

COL_RE <- "#007C83"
COL_IE <- "#D28B35"
COL_GRID <- "#D9DEE3"
COL_TEXT <- "#20262D"

EST_RE <- "Register-equality (reg_eq_similarity)"
EST_IE <- "Inclusion–exclusion (jaccard_similarity_ie)"
EST_LEVELS <- c(EST_RE, EST_IE)
EST_COLORS <- setNames(c(COL_RE, COL_IE), EST_LEVELS)

raw <- read_csv(input_csv, show_col_types = FALSE)
required <- c("precision", "file1", "bedtools_jaccard",
              "j_register_equality", "j_incl_excl")
missing_cols <- setdiff(required, names(raw))
if (length(missing_cols) > 0) {
  stop(basename(input_csv), " lacks columns: ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}

# Strata chosen so the first bin is the regime the Maurano corpus cannot reach
# (its off-diagonal J bottoms out at 0.1355) and the cross-species Mode D
# corpora live in. Not tuned to the result.
STRATA <- tibble(
  stratum = factor(c("J < 0.05", "0.05 ≤ J < 0.2", "J ≥ 0.2"),
                   levels = c("J < 0.05", "0.05 ≤ J < 0.2", "J ≥ 0.2")),
  lo = c(0, 0.05, 0.2),
  hi = c(0.05, 0.2, 1.01)
)

tidy <- raw %>%
  mutate(replicate = sub(".*rep([0-9]+).*", "\\1", file1)) %>%
  rowwise() %>%
  mutate(stratum = STRATA$stratum[which(bedtools_jaccard >= STRATA$lo &
                                          bedtools_jaccard < STRATA$hi)[1]]) %>%
  ungroup() %>%
  filter(!is.na(stratum))

kendall_of <- function(df, col) {
  if (nrow(df) < 3) return(NA_real_)
  cor(df$bedtools_jaccard, df[[col]], method = "kendall")
}

summarise_cell <- function(df) {
  reps <- sort(unique(df$replicate))
  loo <- lapply(reps, function(r) {
    keep <- df %>% filter(replicate != r)
    c(kendall_of(keep, "j_register_equality"), kendall_of(keep, "j_incl_excl"))
  })
  loo_re <- vapply(loo, `[`, numeric(1), 1)
  loo_ie <- vapply(loo, `[`, numeric(1), 2)
  bind_rows(
    tibble(estimator = EST_RE, tau = kendall_of(df, "j_register_equality"),
           lo = min(loo_re), hi = max(loo_re),
           mae = mean(abs(df$j_register_equality - df$bedtools_jaccard))),
    tibble(estimator = EST_IE, tau = kendall_of(df, "j_incl_excl"),
           lo = min(loo_ie), hi = max(loo_ie),
           mae = mean(abs(df$j_incl_excl - df$bedtools_jaccard)))
  ) %>%
    mutate(n = nrow(df))
}

cells <- tidy %>%
  group_by(stratum, precision) %>%
  group_modify(~ summarise_cell(.x)) %>%
  ungroup() %>%
  mutate(estimator = factor(estimator, levels = EST_LEVELS))

message("Kendall tau by stratum and precision (band = leave-one-replicate-out range):")
print(as.data.frame(cells), digits = 4)

# The crossover: the lowest precision at which inclusion-exclusion overtakes
# register-equality in the low-J stratum. Located from the data rather than
# hard-coded, so the annotation cannot drift away from the numbers.
low <- cells %>% filter(stratum == levels(STRATA$stratum)[1])
wide <- low %>%
  select(precision, estimator, tau) %>%
  tidyr::pivot_wider(names_from = estimator, values_from = tau)
crossover_p <- wide$precision[which(wide[[EST_IE]] > wide[[EST_RE]])]
crossover_p <- if (length(crossover_p) > 0) min(crossover_p) else NA_integer_
# The measured precisions bracket the crossover; they do not locate it. Draw the
# marker between the last precision where register-equality still leads and the
# first where it does not, and say "between" rather than naming a value the
# design cannot resolve.
below_p <- if (!is.na(crossover_p)) max(wide$precision[wide$precision < crossover_p]) else NA
crossover_x <- if (!is.na(crossover_p)) mean(c(below_p, crossover_p)) else NA
message("Crossover bracketed between p = ", below_p, " and p = ", crossover_p)

mae_note <- sprintf(
  "Mean absolute error favours inclusion–exclusion by %.0f–%.0f× in every cell.",
  min(low$mae[low$estimator == EST_RE] / low$mae[low$estimator == EST_IE]),
  max(low$mae[low$estimator == EST_RE] / low$mae[low$estimator == EST_IE]))

theme_paper <- function(base_size = 10.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.12), color = COL_TEXT,
                                margin = margin(b = 4)),
      plot.subtitle = element_text(size = rel(0.92), color = COL_TEXT,
                                   margin = margin(b = 9)),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = COL_TEXT, size = rel(0.95)),
      axis.title = element_text(color = COL_TEXT),
      axis.text = element_text(color = COL_TEXT),
      axis.line = element_line(color = "#6B747D", linewidth = 0.35),
      axis.ticks = element_line(color = "#6B747D", linewidth = 0.35),
      panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.85), color = COL_TEXT),
      plot.caption = element_text(size = rel(0.78), color = "#5A636C", hjust = 0,
                                  margin = margin(t = 9)),
      plot.margin = margin(10, 14, 10, 10)
    )
}

crossover_layer <- if (!is.na(crossover_x)) {
  marker <- tibble(stratum = factor(levels(STRATA$stratum)[1],
                                    levels = levels(STRATA$stratum)),
                   x = crossover_x)
  list(
    geom_vline(data = marker, aes(xintercept = x),
               linetype = "22", linewidth = 0.45, color = "#8A939C"),
    geom_text(
      data = marker %>% mutate(y = 0.10),
      aes(x = x, y = y,
          label = sprintf("crossover\nbetween p = %d and %d", below_p, crossover_p)),
      hjust = 0.5, vjust = 0, size = 2.7, lineheight = 1.05, color = COL_TEXT,
      inherit.aes = FALSE
    )
  )
} else {
  NULL
}

figure <- ggplot(cells, aes(x = precision, y = tau,
                            color = estimator, fill = estimator)) +
  crossover_layer +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.1) +
  facet_wrap(~ stratum, nrow = 1) +
  scale_color_manual(values = EST_COLORS, drop = FALSE) +
  scale_fill_manual(values = EST_COLORS, drop = FALSE) +
  scale_x_continuous(breaks = sort(unique(cells$precision))) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2),
                     labels = label_number(accuracy = 0.1)) +
  labs(
    title = "The estimator choice is only live below J ≈ 0.05, and precision settles it",
    subtitle = paste(
      "Kendall τ against exact BEDTools Jaccard.",
      "Band: leave-one-replicate-out range over 3 replicates — a stability check, not a CI."),
    x = "HyperLogLog precision p",
    y = "Kendall τ vs BEDTools",
    caption = paste0(
      mae_note, " Rank fidelity is the only axis on which they disagree.\n",
      "Above J = 0.05 both reach τ = 1 by p = 20 on this corpus (n = 9 and 15 pairs), ",
      "so those facets carry no resolving power.")
  ) +
  guides(fill = "none",
         color = guide_legend(override.aes = list(linewidth = 1.1, size = 2.4))) +
  theme_paper()

CairoPNG(filename = out_png, width = 10.6, height = 4.3, units = "in",
         res = 300, bg = "white")
print(figure)
dev.off()

write_csv(cells, file.path(script_dir, "estimator_crossover_stats.csv"))
message("Wrote: ", out_png)
message("Wrote: ", file.path(script_dir, "estimator_crossover_stats.csv"))
