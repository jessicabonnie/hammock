# synthetic_evolution

Recreates the original `synthetic-c` experiment from
`../hammock/experiments/synthetic-c` against the refactored hammock in this
repo. Tests the hypothesis that:

| Evolution type | Mechanism | Predicted best detector |
|----------------|-----------|-------------------------|
| **Indel** (A)    | Intervals appear/disappear (set membership churns) | Mode A (interval-level)    |
| **Jitter** (B)   | Interval positions drift ±100 bp (membership preserved) | Mode B (point-level)       |
| **Combined** (C) | 50/50 mix of Indel and Jitter changes               | Mode C (interval + points) |

## Layout

```
experiments/synthetic_evolution/
├── Makefile          # drives data → evolve → hammock → analyze
├── code/
│   ├── createBED.py  # base BED generator (ported from original)
│   ├── evolve.py     # generation-stepped evolution simulator
│   └── analyze.R     # CSVs → figures + summary
├── data/      → /vast/…/synthetic_evolution/data
│   ├── synthetic_100000.bed         # base BED (seed-fixed)
│   ├── synthetic_base.txt           # 1-line path file (bed1 for hammock)
│   └── evolved/                     # 20 snapshot BEDs per (rate, type)
│       ├── synthetic_g100_r0.05_A_005.bed … _100.bed
│       └── list_r0.05_A.txt         # path file (bed2 for hammock)
├── results/   → /vast/…/synthetic_evolution/results
│   ├── hammock/      # per-(rate, type, mode[, subB]) CSVs
│   └── summary.csv   # Jaccard at gen=50, one row per (rate, type, mode, subB)
├── logs/      → /vast/…/synthetic_evolution/logs
└── figures/                         # PNGs (committed)
    ├── jaccard_by_generation.png    # headline grid
    ├── sensitivity_heatmap.png      # Jaccard at gen=50 by mode × type × rate
    └── mode_c_subB_sweep.png        # Mode C subB effect
```

Bulk artifacts (`data/`, `results/`, `logs/`) are symlinks to
`/vast/blangme2/jbonnie/hammock_claude_experiments/synthetic_evolution/`;
only `code/`, `Makefile`, `README.md`, and `figures/` live in git.

## Run

```bash
# Full sweep (≈10 min single-node, ≈3 min with `make -j 6`).
cd experiments/synthetic_evolution
make all

# Smoke-test (one rate, one type, fast).
make RATES=0.05 TYPES=A all
```

Knobs (override on the make command line):

| Variable      | Default                       | Note                      |
|---------------|-------------------------------|---------------------------|
| `RATES`       | `0.005 0.01 0.05 0.1`         | Per-generation change rate |
| `TYPES`       | `A B C`                       | Evolution types            |
| `C_SUBB`      | `0.25 0.50 0.75 1.00`         | Mode C subsampling sweep   |
| `N_INTERVALS` | `100000`                      | Base BED size              |
| `GENERATIONS` | `100`                         | Generations to evolve      |
| `STEP`        | `5`                           | Snapshot every N gens      |
| `PRECISION`   | `24`                          | HLL precision              |
| `THREADS`     | `8`                           | Per-hammock-invocation     |
| `SEED`        | `42`                          | Fixed for reproducibility  |
| `HAMMOCK`     | `…/claude-ref-comparison/bin/hammock` | hammock_claude binary |

## What the analysis tests

For each (rate, type, mode), hammock computes Jaccard(base, evolved_at_gen_t)
for t ∈ {5, 10, …, 100}. The **rate of Jaccard decay** with generation tells
us how sensitive that mode is to that flavor of change. The headline figure
plots all 12 (rate × type) panels with all three modes overlaid; the
sensitivity heatmap collapses to one number per cell (Jaccard at gen 50,
lower ⇒ more sensitive).

Mode C in `hammock_claude` defaults to `--subB-method=mixed-stride`, which is
**not** parity-comparable to the original hammock's hash-threshold default.
See `CLAUDE.md` divergence #3. Pass `--subB-method=hash-threshold` to
hammock for original parity if needed.

## Divergences from the original

- 100k intervals (not 10k as the README in `../hammock/experiments/synthetic-c`
  inconsistently claims — the actual filenames there encode 100k).
- Focused rate sweep `{0.005, 0.01, 0.05, 0.1}`; the higher rates in the
  original (`0.5, 0.8, 1.0`) saturate Jaccard near 0 within a few
  generations and don't add resolution.
- `--expA` removed (it was a placeholder-column multiplier and is gone from
  hammock_claude; see `CLAUDE.md`).
- Mode C subB sweep is run at **all** rates, not just at the saturated `r=1.0`
  as in the original.
- HLL containment/cosketch columns are emitted by hammock_claude but ignored
  by this analysis — the headline plot uses `jaccard_similarity` only.
