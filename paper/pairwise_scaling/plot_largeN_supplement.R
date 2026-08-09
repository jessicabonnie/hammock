# Supplementary Figure — hammock at catalog scale (N up to 2048)
#
# A SEPARATE figure, not a panel of Figure 3. Figure 3 compares two measured
# curves; this one shows a measured hammock curve against a PROJECTED bedtools
# curve, and mixing a measurement with an extrapolation inside the same panel as
# two measurements is exactly how a projection gets quoted later as if it had
# been run.
#
# WHY PROJECT AT ALL. bedtools is Theta(N^2) with a large constant -- measured
# 3.97x per doubling, 495 s per replicate at N=256 -- so a single replicate costs
# ~2.2 h at N=1024 and ~8.6 h at N=2048. Three replicates at N=2048 is over a day
# of node time to extend a comparison that is already decided by N=256 (14.75x).
# hammock is near-linear in the same regime because at p=18 it is
# sketching-dominated, so measuring IT at catalog scale is minutes, not days.
# N~2048 is a real corpus, not a round number: the ChIP-Atlas manifest holds
# 2,206 verified hg38 CTCF files.
#
# THREE THINGS THAT ARE LOAD-BEARING:
#
#   1. The projected segment is drawn dashed, in a lighter shade, labelled in
#      the legend as "projected", and annotated on the panel. It must never be
#      mistakable for a measurement at a glance.
#   2. The fitted exponent is PRINTED, not assumed to be 2. If bedtools' fit
#      comes out far from 2 the projection is not trustworthy and the caption
#      should say so rather than quoting a speedup.
#   3. N=512 is measured by BOTH jobs and the script reports the disagreement.
#      Two separate allocations on this cluster are not guaranteed comparable
#      (see docs/seed-benchmark-methodology.md); the overlap point is the only
#      thing that makes the join checkable, so a large gap there invalidates the
#      combined figure and the script says so out loud.

required_packages <- c("dplyr", "readr", "ggplot2", "scales", "patchwork", "Cairo")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
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
panelA_csv <- if (length(argv) >= 1) argv[1] else
  file.path(repo_root, "docs", "data", "cpp_vs_bedtools_t16_p18.csv")
largeN_csv <- if (length(argv) >= 2) argv[2] else
  file.path(repo_root, "docs", "data", "cpp_vs_bedtools_t16_p18_largeN.csv")
out_png <- if (length(argv) >= 3) argv[3] else
  file.path(repo_root, "paper", "figures", "largeN_supplement.png")

for (f in c(panelA_csv, largeN_csv)) {
  if (!file.exists(f)) stop("Input not found: ", f, call. = FALSE)
}

read_tool <- function(path, want) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(tool == want, !is.na(mean_wall_time)) %>%
    transmute(num_files, wall = mean_wall_time,
              lo = min_wall_time, hi = max_wall_time)
}

bt   <- read_tool(panelA_csv, "bedtools")
hm_a <- read_tool(panelA_csv, "hammock_cpp_B")
hm_b <- read_tool(largeN_csv, "hammock_cpp_B")

if (nrow(bt) == 0 || nrow(hm_a) == 0 || nrow(hm_b) == 0)
  stop("One of the three series is empty; check the two input CSVs.", call. = FALSE)

# --- the overlap check, before anything is plotted -------------------------
overlap <- inner_join(hm_a, hm_b, by = "num_files", suffix = c("_A", "_B"))
if (nrow(overlap) == 0) {
  stop("The two jobs share no N. Without an overlap point the join between a ",
       "Panel A allocation and a large-N allocation cannot be checked, and on ",
       "this cluster two allocations are not guaranteed comparable. Re-run the ",
       "large-N job with an N that Panel A also measured.", call. = FALSE)
}
overlap <- overlap %>% mutate(ratio = wall_B / wall_A)
cat("Overlap check (hammock measured by both jobs):\n")
overlap %>%
  transmute(N = num_files, panelA = round(wall_A, 2), largeN = round(wall_B, 2),
            ratio = round(ratio, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
worst <- max(abs(log(overlap$ratio)))
if (worst > log(1.25)) {
  cat("\n*** WARNING: the two jobs disagree by more than 25% on a shared N.\n",
      "*** Treat the combined figure as unreliable and say so, or re-run both\n",
      "*** in one allocation. Not stopping, because seeing the plot is how you\n",
      "*** diagnose this -- but do not publish it in this state.\n\n", sep = "")
} else {
  cat(sprintf("\nMax disagreement %.1f%% -- within the noise this cluster shows; join OK.\n\n",
              100 * (exp(worst) - 1)))
}

# --- fit bedtools' scaling on the MEASURED points only ---------------------
# Fitted on N >= 32: below that the run is dominated by fixed process/module
# startup that does not scale with pair count, and including it drags the
# exponent down and would make the projection optimistic in our favour.
bt_fit_data <- bt %>% filter(num_files >= 32)
fit <- lm(log(wall) ~ log(num_files), data = bt_fit_data)
expo <- unname(coef(fit)[2])
cat(sprintf("bedtools fitted exponent on N >= 32: %.3f (theory 2.0), R^2 = %.4f\n",
            expo, summary(fit)$r.squared))
if (abs(expo - 2) > 0.25)
  cat("*** WARNING: exponent is far from 2; the projection below is not trustworthy.\n")

proj_N <- setdiff(hm_b$num_files, bt$num_files)
bt_proj <- tibble(num_files = proj_N,
                  wall = exp(predict(fit, newdata = data.frame(num_files = proj_N))))
cat("\nProjected bedtools (NOT measured):\n")
bt_proj %>% mutate(hours = round(wall / 3600, 2), wall = round(wall)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# --- assemble ---------------------------------------------------------------
hammock_all <- bind_rows(hm_a %>% mutate(src = "Panel A job"),
                         hm_b %>% mutate(src = "large-N job")) %>%
  group_by(num_files) %>%                       # average the overlap point
  summarise(wall = mean(wall), lo = min(lo), hi = max(hi), .groups = "drop")

series <- bind_rows(
  bt      %>% mutate(series = "BEDTools (measured)",  kind = "measured"),
  bt_proj %>% mutate(lo = NA, hi = NA,
                     series = "BEDTools (projected)", kind = "projected"),
  hammock_all %>% mutate(series = "hammock (measured)", kind = "measured"))

# Connect the measured bedtools curve to its projection so the eye follows one
# line, but keep the projected part visually distinct.
bt_link <- bind_rows(bt %>% filter(num_files == max(num_files)) %>%
                       transmute(num_files, wall),
                     bt_proj) %>%
  mutate(series = "BEDTools (projected)", kind = "projected")

pal <- c("BEDTools (measured)"  = "#000000",
         "BEDTools (projected)" = "#8c8c8c",
         "hammock (measured)"   = "#ff7f0e")
Ns <- sort(unique(series$num_files))

pA <- ggplot() +
  geom_ribbon(data = series %>% filter(kind == "measured", !is.na(lo)),
              aes(num_files, ymin = lo, ymax = hi, fill = series),
              alpha = 0.15, colour = NA) +
  geom_line(data = series %>% filter(series != "BEDTools (projected)"),
            aes(num_files, wall, colour = series), linewidth = 0.7) +
  geom_line(data = bt_link, aes(num_files, wall, colour = series),
            linewidth = 0.7, linetype = "22") +
  geom_point(data = series %>% filter(kind == "measured"),
             aes(num_files, wall, colour = series), size = 1.8) +
  geom_point(data = bt_proj %>% mutate(series = "BEDTools (projected)"),
             aes(num_files, wall, colour = series), size = 1.8, shape = 1) +
  annotate("text", x = max(Ns), y = max(bt_proj$wall) * 1.6,
           label = "projected, not measured", hjust = 1, size = 3, colour = "grey35") +
  scale_x_continuous(trans = "log2", breaks = Ns) +
  scale_y_log10(breaks = breaks_log(n = 7),
                labels = label_number(drop0trailing = TRUE)) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(x = "Files per side (N)", y = "Wall time (s)",
       title = "A. Cost at catalog scale") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

speed <- series %>%
  filter(series != "hammock (measured)") %>%
  select(num_files, bt_wall = wall, kind) %>%
  inner_join(hammock_all %>% select(num_files, hm_wall = wall), by = "num_files") %>%
  mutate(speedup = bt_wall / hm_wall)

pB <- ggplot(speed, aes(num_files, speedup)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_line(data = speed %>% filter(kind == "measured"), linewidth = 0.7,
            colour = "#ff7f0e") +
  geom_line(data = speed %>% filter(num_files >= max(bt$num_files)),
            linewidth = 0.7, colour = "grey50", linetype = "22") +
  geom_point(aes(shape = kind), size = 2, colour = "#ff7f0e") +
  scale_shape_manual(values = c(measured = 16, projected = 1), name = NULL,
                     labels = c(measured = "measured", projected = "projected")) +
  scale_x_continuous(trans = "log2", breaks = Ns) +
  scale_y_log10(breaks = breaks_log(n = 7),
                labels = label_number(drop0trailing = TRUE)) +
  labs(x = "Files per side (N)", y = "hammock speedup over BEDTools",
       title = "B. Speedup, measured and projected") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

fig <- pA | pB
dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
CairoPNG(out_png, width = 2000, height = 900, res = 200)
print(fig)
invisible(dev.off())
cat("\nWrote:", out_png, "\n\n")

speed %>%
  transmute(N = num_files, bedtools = round(bt_wall, 1),
            hammock = round(hm_wall, 2), speedup = round(speedup, 1), kind) %>%
  as.data.frame() %>% print(row.names = FALSE)
