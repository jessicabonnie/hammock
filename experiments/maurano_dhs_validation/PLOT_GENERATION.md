# Plot generation: Mode D Pearson vs ARI sweep summaries

You're being handed a self-contained plotting task. The paper outline
(`docs/paper_outline.md`) has placeholder captions for **Fig 4** and
**Fig 7** that need to be filled with three candidate figures (one is
the committed choice; two are candidates the project lead wants to
compare side-by-side before final selection).

## Context

The Mode D sweep on Maurano DHS produces, for every `(precision p, k,
w, column-flavor)` config, a Pearson r vs bedtools and an ARI vs the
10 fetal-tissue labels. The headline scientific point is that the
**(k, w, p) optimum for matching bedtools is in a different part of
the sweep than the optimum for recovering tissue structure** —
Pearson is a high-k / high-w / high-p ridge; ARI is a single
(k=10, w=30) peak that lights up starting at p = 12.

The current Fig 4 (`mode_d_pearson_heatmap.png`) and Fig 7
(`mode_d_clustering_ari.png`) are both 9-precision × 2-column-flavor
grids of (k, w) heatmaps. Project lead found these hard to read and
asked for line / violin alternatives that show the contrast more
directly.

## Input data

CSV: `experiments/maurano_dhs_validation/results/mode_d_summary.csv`
(symlink into `/vast/blangme2/jbonnie/hammock_claude_experiments/`).

Columns:
```
precision, k, w, column, reference, n, pearson, spearman, mae, rmse,
frob, max_err, ari, nmi, path
```

For every plot below:
- **Filter to `column == "jaccard_similarity"`** (no_ends only — the
  with_ends flavor is dominated for both metrics on Maurano and lives
  in §4.5)
- **Filter to `reference == "bedtools"`** (the per-config Pearson /
  Spearman / MAE in the CSV are against either `bedtools` or `Mode B`;
  use bedtools)

After filtering, expect roughly one row per `(precision, k, w)` cell.
There are 9 precisions (10, 12, 14, 16, 18, 20, 22, 23, 24) and
varying k (8, 10, 15, 20, 25) × w (8, 10, 20, 30, 50, 100, 200,
300, 500) cells.

## Environment

```bash
ml r/4.3.0
# (the claude-ref-comparison conda env does not bundle R)
```

**Use `Cairo::CairoPNG()` for PNG output** — system R 4.3.0 has no
native png/cairo capability and `ggsave(..., device="png")` fails
with "X11 is not available". `library(Cairo)` then wrap each save:

```r
library(Cairo)
CairoPNG(out_path, width = 1200, height = 700, res = 150)
print(p)
dev.off()
```

## Outputs

Save all three to `experiments/maurano_dhs_validation/figures/`:

1. `mode_d_lines_p24.png` — **Option B** (committed: replaces Fig 4
   in the paper outline)
2. `mode_d_max_lines_p24.png` — **Option C** (candidate: project
   lead wants to evaluate vs Option B)
3. `mode_d_violins_by_k.png` — **violin plot** (committed: replaces
   Fig 7 in the paper outline)

Save the generating script as
`experiments/maurano_dhs_validation/scripts/make_metric_plots.R` so
the figures are regenerable.

## Plot 1 — Option B (`mode_d_lines_p24.png`)

Two-panel line plot side-by-side, both at p = 24.

- **Left panel:** Pearson r vs w; one colored line per k.
- **Right panel:** ARI vs w; same color/legend.
- Shared y-axis on `[0, 1]`.
- x-axis: `w` on a log-2 scale (8, 10, 20, 30, 50, 100, 200, 300, 500
  spans roughly 6 doublings — log scale keeps the points evenly spaced).
- Color: `k` as a discrete categorical with 5 levels (8, 10, 15, 20,
  25), viridis or any sequential palette.
- Title: `Mode D — Pearson r and ARI vs window size, at p = 24`.
- Subtitle: `no_ends column only; lines connect points within the same k.`

```r
library(tidyverse); library(Cairo)

d <- read_csv("experiments/maurano_dhs_validation/results/mode_d_summary.csv") |>
  filter(column == "jaccard_similarity",
         reference == "bedtools",
         precision == 24) |>
  mutate(k = factor(k))

p_left <- d |> ggplot(aes(w, pearson, colour = k, group = k)) +
  geom_line() + geom_point() +
  scale_x_log10() + ylim(0, 1) +
  labs(x = "window size w", y = "Pearson r vs bedtools",
       title = "Pearson r vs bedtools")

p_right <- d |> ggplot(aes(w, ari, colour = k, group = k)) +
  geom_line() + geom_point() +
  scale_x_log10() + ylim(0, 1) +
  labs(x = "window size w", y = "ARI vs tissue labels",
       title = "ARI vs tissue labels")

p <- patchwork::wrap_plots(p_left, p_right, ncol = 2, guides = "collect") +
     patchwork::plot_annotation(
       title = "Mode D — Pearson r and ARI vs window size, at p = 24",
       subtitle = "no_ends column only; lines connect points within the same k"
     )

CairoPNG("experiments/maurano_dhs_validation/figures/mode_d_lines_p24.png",
         width = 1400, height = 600, res = 150)
print(p); dev.off()
```

## Plot 2 — Option C (`mode_d_max_lines_p24.png`)

Single panel; **for each k**, take the max-over-w Pearson and
max-over-w ARI at `p = 24`, then plot both as a function of k.

- x-axis: `k` (8, 10, 15, 20, 25), numeric (not factor — meaningful
  spacing matters here).
- Two lines:
  - "best Pearson at this k" (max over w)
  - "best ARI at this k" (max over w)
- Shared y-axis `[0, 1]`.
- Distinct colors for the two lines.
- Annotate each point with the `w` that achieved the maximum, so the
  reader sees that Pearson-best lives at high w (100+) while ARI-best
  lives at w = 30.
- Title: `Mode D — best Pearson and best ARI per k, at p = 24`.
- Subtitle: `no_ends only; each point is the (w) maximum at that (k, p).`

```r
d24 <- read_csv("experiments/maurano_dhs_validation/results/mode_d_summary.csv") |>
  filter(column == "jaccard_similarity",
         reference == "bedtools",
         precision == 24)

best <- d24 |>
  group_by(k) |>
  summarise(pearson_max = max(pearson),
            pearson_w   = w[which.max(pearson)],
            ari_max     = max(ari),
            ari_w       = w[which.max(ari)],
            .groups = "drop") |>
  pivot_longer(c(pearson_max, ari_max),
               names_to = "metric", values_to = "value") |>
  mutate(best_w = ifelse(metric == "pearson_max", pearson_w, ari_w))

p <- best |> ggplot(aes(k, value, colour = metric, group = metric)) +
  geom_line() + geom_point(size = 3) +
  geom_text(aes(label = paste0("w=", best_w)),
            vjust = -1, size = 3, show.legend = FALSE) +
  ylim(0, 1) +
  scale_colour_manual(values = c(pearson_max = "#1f77b4",
                                  ari_max     = "#d62728"),
                      labels = c(pearson_max = "max Pearson r",
                                  ari_max     = "max ARI")) +
  labs(x = "k-mer size k", y = "metric value",
       title = "Mode D — best Pearson and best ARI per k, at p = 24",
       subtitle = "no_ends only; each point is the max over w at that (k, p); annotation shows the winning w")

CairoPNG("experiments/maurano_dhs_validation/figures/mode_d_max_lines_p24.png",
         width = 1100, height = 700, res = 150)
print(p); dev.off()
```

## Plot 3 — Violin (`mode_d_violins_by_k.png`)

Per-k violin plot, two metrics, **across the full sweep of (w, p)**
at each k (not restricted to p = 24).

- x-axis: `k` (8, 10, 15, 20, 25) as a factor.
- y-axis: metric value on `[0, 1]`.
- Two violins per k, side-by-side via `position_dodge`:
  - Pearson r vs bedtools (one color)
  - ARI vs tissue labels (a different color)
- Each violin's distribution is over every (w, p) combination at that
  k, no_ends column, bedtools reference. Expect roughly 9 w × 9 p =
  ~81 points per violin (fewer where high-k × low-p cells were
  dropped upstream).
- Show interior box (`geom_boxplot(width=0.05)`) and individual points
  (`geom_jitter`) inside each violin so the bimodality at k = 10 ARI
  is visible.
- Title: `Mode D metric distributions per k, across the full (w, p) sweep`.
- Subtitle: `no_ends only; Pearson migrates upward with k, ARI is bimodal at k=10 only.`

```r
d_long <- read_csv("experiments/maurano_dhs_validation/results/mode_d_summary.csv") |>
  filter(column == "jaccard_similarity", reference == "bedtools") |>
  select(k, w, precision, pearson, ari) |>
  pivot_longer(c(pearson, ari), names_to = "metric", values_to = "value") |>
  mutate(k = factor(k),
         metric = recode(metric,
                         pearson = "Pearson r vs bedtools",
                         ari     = "ARI vs tissue labels"))

p <- d_long |>
  ggplot(aes(k, value, fill = metric)) +
  geom_violin(position = position_dodge(width = 0.85),
              alpha = 0.6, scale = "width") +
  geom_boxplot(position = position_dodge(width = 0.85),
               width = 0.12, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(aes(colour = metric),
              position = position_jitterdodge(jitter.width = 0.12,
                                              dodge.width = 0.85),
              size = 0.6, alpha = 0.4, show.legend = FALSE) +
  ylim(0, 1) +
  scale_fill_manual(values = c("Pearson r vs bedtools" = "#1f77b4",
                                "ARI vs tissue labels" = "#d62728")) +
  scale_colour_manual(values = c("Pearson r vs bedtools" = "#1f77b4",
                                  "ARI vs tissue labels" = "#d62728")) +
  labs(x = "k-mer size k", y = "metric value", fill = NULL,
       title = "Mode D metric distributions per k, across the full (w, p) sweep",
       subtitle = "no_ends only; Pearson migrates upward with k; ARI is bimodal at k=10 only")

CairoPNG("experiments/maurano_dhs_validation/figures/mode_d_violins_by_k.png",
         width = 1300, height = 700, res = 150)
print(p); dev.off()
```

## After generation

1. Inspect each PNG. The violin plot is the most novel — confirm that
   the ARI violin at k = 10 visibly stretches to ~0.91 with a long
   tail / bimodal shape; if it looks unimodal-low you've probably
   filtered out the optimum cell by mistake.
2. Update `docs/paper_outline.md`:
   - Remove the `<!-- ... -->` HTML comment around the Fig 4 image
     link so the new `mode_d_lines_p24.png` renders.
   - Same for Fig 7's `mode_d_violins_by_k.png`.
   - Drop the "(planned — pending generation)" tag from both captions
     and trim the "generation spec" pointer line.
3. `mode_d_max_lines_p24.png` (Option C) is **not** wired into the
   outline yet — leave it on disk for the project lead to evaluate.
   They will decide whether to swap it in over Option B as Fig 4.
4. Save the R script as `experiments/maurano_dhs_validation/scripts/make_metric_plots.R`
   so all three figures are regenerable in one shot.
5. Commit:
   - `figures/mode_d_lines_p24.png`
   - `figures/mode_d_max_lines_p24.png`
   - `figures/mode_d_violins_by_k.png`
   - `scripts/make_metric_plots.R`
   - `docs/paper_outline.md` (the placeholder unwiring)

## Sanity checks

- At p = 24, max Pearson should be **0.9996** (at k = 20, w = 100
  approximately). If your Option B left panel doesn't reach the top
  of the y-axis somewhere on the k = 20 line, you're filtering wrong.
- At p = 24, max ARI should be **0.910** (at k = 10, w = 30). If your
  Option B right panel doesn't peak there, same diagnosis.
- The ARI plateau across precision is `0.9101796407185628` for every
  p ∈ {12, 14, 16, 18, 20, 22, 23, 24} at (k=10, w=30) — confirm in
  the violin's upper lobe.
- **k ≥ 15 ARI degeneracy:** at k ∈ {15, 20, 25} every (w, p) cell
  in the sweep gives ARI = `0.6931...` exactly (because the high-k
  clustering hits a fixed basin). In the violin plot these three k
  values should collapse to a flat horizontal stripe at 0.69 — *not*
  a violin shape. If you see any spread at all in the k = 15 / 20 /
  25 ARI violins you're either picking up a different `column` or
  `reference` than the filter intends.
- The k = 10 ARI distribution spans roughly the full [0, 0.91] range
  with multiple sub-plateaus (around 0.09, 0.55, 0.69, 0.91); the
  violin will be tall and multimodal, not a smooth single bump.

## Notes

- Don't include `with_ends` data in any of these three figures. The
  `jaccard_similarity_with_ends` flavor is dominated on Maurano (Fig
  S4 / §4.5 already cover it). Keeping the new figures `no_ends`-only
  keeps them legible.
- The script should be deterministic and read directly from
  `results/mode_d_summary.csv` (a symlink into vast storage). No SLURM
  job needed — this runs in a few seconds on the head node.
- The CairoPNG note above is load-bearing — `ggsave()` defaults to a
  PNG device that requires X11; on this cluster's R 4.3.0 module
  that's missing.
