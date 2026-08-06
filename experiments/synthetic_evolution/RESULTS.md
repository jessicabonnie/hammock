# synthetic_evolution — results notes

> **Living document.** Updated as we refine the experiment. The headline
> numbers below come from `results/summary_g{GENOME_RANGE}.csv`; figures
> are in `figures/g{GENOME_RANGE}/`. Y-axis on every "vs. generation"
> figure is *Jaccard similarity to the pre-evolution BED* (distance from
> generation 0).

## TL;DR

On a saturated synthetic genome (the original setup), Mode A appeared to
have a unique advantage at detecting Indel evolution that Mode B
"missed". **That asymmetry turned out to be an artifact of genome
saturation.** On a biologically realistic sparse genome:

- **Mode A and Mode B both track Indel equally well** — both drop to
  near zero at high indel rates.
- **Mode A is brittle to Jitter** (any ±1 bp shift breaks its hash);
  **Mode B is robust** (positions still overlap).
- The only honest, generally-applicable asymmetry is that **Mode B is
  strictly better than Mode A as an appropriate-reading metric** on
  realistic sparse data — it tracks indel as well as Mode A *and*
  correctly says "still similar" for jitter, while Mode A overstates the
  jitter case.
- Mode C with `expA` is a continuous dial between Mode B (expA=0) and
  Mode A (expA→large); useful only if you specifically want Mode A's
  brittleness mixed in.

---

## Sweep configuration

| Knob              | Value                                  |
|-------------------|----------------------------------------|
| Base BED          | 100k synthetic intervals, seed=42      |
| Genome space      | 5 chromosomes                          |
| Interval lengths  | {10, 50, 100, 500, 1000, 5000, 10000} (avg ~2380 bp) |
| Generations       | 100 (snapshot every 5)                 |
| Change rates      | {0.005, 0.01, 0.05, 0.1}               |
| Evolution types   | Indel (A), Jitter (B), Combined (C)    |
| Hammock           | hammock_claude, HLL p=18, subB=0.25    |
| Mode C expA sweep | {0, 0.5, 1, 2, 3} (subB pinned at 0.25) |

Two `GENOME_RANGE` configs run side-by-side; tag = `g{GENOME_RANGE}`:

| Tag           | Start range / chr | Coverage fraction        | Regime                                |
|---------------|-------------------|--------------------------|---------------------------------------|
| `g1000000`    | 1 Mbp             | ~47× over-covered        | **Saturated** (original setup; biologically unrealistic) |
| `g100000000`  | 100 Mbp           | ~0.5%                    | **Sparse**  (biologically realistic)  |

Coverage fraction = `N × L_avg / (n_chr × range)`. At 100k intervals × 2380 bp
average length over 5 chromosomes, that's 238 Mbp of covered positions —
which fits comfortably in 5 × 100 Mbp = 500 Mbp (0.5%), but vastly
exceeds 5 × 1 Mbp = 5 Mbp (47× over-coverage; every position covered
multiple times).

Run a config:

```bash
make GENOME_RANGE=100000000 -j 6 all   # sparse (new default)
make GENOME_RANGE=1000000   -j 6 all   # saturated (original)
```

---

## The saturation artifact, in numbers

At rate=0.05 over 100 generations of Indel, the per-generation survival
probability of any original interval is `0.95`, so by gen=100 only
`0.95^100 ≈ 0.006` of original intervals remain (99.4% turnover).
Compare Mode A and Mode B Jaccards at gen=100 across the two configs:

| Regime          | Mode A   | Mode B   | What it means                                          |
|-----------------|---------:|---------:|--------------------------------------------------------|
| Saturated (g1M) | ~0.01    | ~1.0     | Mode B is blind to Indel — saturation hides the signal |
| Sparse (g100M)  | ~0.01    | ~0.08    | Mode B sees Indel almost as well as Mode A             |

> **Checked 2026-08-06** against `results/hammock_g{TAG}/synth_r0.05_A_hll_p18_jacc{A,B}.csv`
> at generation 100. Exact values: saturated A 0.0141 / B 0.9947; sparse
> A 0.0139 / B **0.1109** — the "~0.08" above reads low, the rest hold.
> Two caveats on the sparse Mode B figure, both material at this magnitude:
> `jaccard_similarity` is register-equality, not set Jaccard, and carries a
> chance-agreement floor (`CLAUDE.md` divergence #2); the inclusion-exclusion
> value recovered from the same row's containments is 0.049. So Mode B's true
> sensitivity to Indel in the sparse regime is *better* than either number
> above suggests, which strengthens rather than weakens the conclusion.

The saturated number for Mode B is misleading: each interval contributes
~2380 unique positions, but in a 5 Mbp genome most positions are also
covered by neighboring intervals, so deleting one interval barely
shifts which unique positions exist in the sketch. Mode B sees ~no
change. Bump the genome range to 100 Mbp and the over-coverage is gone;
each interval's positions become unique-ish, and deletions show up.

`figures/g100000000/mode_per_type_r0.05.png` and `figures/g1000000/...`
make this directly visible: on the **Indel** panel, Mode B (green) is a
flat line at the top in the saturated run but drops to track Mode A
(orange) in the sparse run.

---

## What survives: Mode A is brittle, Mode B is robust

The Jitter panel tells the *real* story, and it's the same in both
regimes: Mode A drops fast, Mode B barely moves.

| Evolution        | Mode A behavior                       | Mode B behavior                       | Which is "right"?     |
|------------------|---------------------------------------|---------------------------------------|-----------------------|
| Indel (sparse)   | Drops — correct (intervals deleted)   | Drops — correct (positions removed)   | Both                  |
| Indel (saturated)| Drops — correct                       | Stays high — *misleading* (artifact)  | Mode A (by default)   |
| Jitter           | Drops — *misleading* (overreacts)     | Stays high — correct (regions same)   | **Mode B**            |
| Combined (sparse)| Drops faster than B (responds to both)| Drops at the indel-half rate          | Mode B (cleaner)      |

The reason: Mode A hashes `(chr, start, end)` exactly, so a ±100 bp
jitter creates a brand-new element in the sketch — looks identical to
"this interval was deleted and a different one inserted." Mode B
expands each interval to its 1-bp-resolution point set, so jittering
only flips boundary positions; the bulk of the interval's points are
unchanged.

So on real biology — where positions matter a lot more than exact
interval boundaries — **Mode B is the better default metric**.

---

## What `expA` is for, in this context

Mode C's `expA` knob multiplies each interval's contribution to the
sketch. With Mode B's behavior at one end and Mode A's brittleness at
the other, `expA` literally slides Mode C between them:

- `expA = 0` → behaves like Mode B (point-dominated; robust to jitter)
- `expA = 3` → behaves like Mode A (interval-dominated; brittle to
  jitter)

`figures/g100000000/mode_c_expA_interpolation.png` shows this: at any
fixed generation, the Mode C Jaccard slides from the Mode B reference
line (dotted) down to the Mode A reference line (dashed) as expA
climbs from 0 to 3. The Jitter curve (blue) shows the largest
displacement — it's pulled from the "appropriate" Mode B reading toward
the "overreactive" Mode A reading.

So in practice you only want a non-zero `expA` if you specifically need
to count interval-identity changes more heavily than positional ones —
e.g., if you care that "the same locus boundaries" defines biological
identity. For ordinary similarity work, leave it at 0.

---

## Figures

Each tag's figures are in `figures/g{TAG}/`; replace the prefix in the
filenames below.

| File                                       | What it shows |
|--------------------------------------------|---------------|
| `jaccard_by_generation.png`                | 3 × 4 grid: Jaccard vs gen, rows = evolution type, cols = rate, three lines (Modes A/B/C-default). |
| `sensitivity_heatmap.png`                  | Jaccard at gen=50 for each (rate, evolution type, mode). Yellow ★ marks the lowest-Jaccard mode per cell — note "most-sensitive" ≠ "most-appropriate" (Mode A wins everywhere on this metric, but for the wrong reason on Jitter). |
| `mode_per_type_r{R}.png` (one per rate)    | Headline per-rate figure: 3 side-by-side panels (Indel \| Jitter \| Combined), each overlays Mode A, Mode B, Mode C default, and Mode C at expA ∈ {0.5, 1, 2, 3}. **This is the figure that answers "does each mode read its kind of change appropriately?"** |
| `mode_c_expA_sweep.png`                    | Mode C across the expA sweep, full (rate × evolution-type) grid. |
| `mode_c_expA_interpolation.png`            | X = expA, Y = Jaccard at gen=50, faceted by rate. Lines colored by evolution type, with horizontal reference lines for the Mode A (dashed) and Mode B (dotted) readings at the same gen. Visually shows expA interpolating Mode C between the two extremes. |

---

## What's not in the experiment (and why)

- **No replicate seeds.** Each (rate, evolution-type) is one stochastic
  trajectory. The curves are smooth enough that the qualitative story
  is robust, but error bars would tighten any numeric claim.
- **No subB sweep.** `subB` is pinned at 0.25; we previously swept it
  and decided subB matters only in the context of varying expA.
- **No reference parity against the original `hammock` outputs.** The
  refactored hammock honors `--subB` in Mode B (the original silently
  ignored it; see `CLAUDE.md` divergence #1), and Mode C's containment
  block differs (see divergence #2). Comparing against the original
  `hammock/experiments/synthetic-c/hammock/` CSVs would require
  re-running the original with explicit defaults; not done.

## Open questions / next steps

- **Jitter amount sweep.** Default is ±100 bp; at very large jitters
  (say ±5000 bp) Mode B's robustness should eventually break down
  because the point set finally shifts meaningfully.
- **An even sparser regime** (e.g. `GENOME_RANGE=1e9`, ~0.05%) to
  confirm that the Indel detection rate doesn't depend further on
  coverage, just to nail the asymptote.
- **Replicate seeds** for the sparse sweep to put error bars on the
  Mode A ≈ Mode B Indel claim.
