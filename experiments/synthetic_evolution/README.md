# synthetic_evolution

Recreates the original `synthetic-c` experiment from
`../hammock/experiments/synthetic-c` against the refactored hammock in this
repo. Tests the hypothesis that:

| Evolution type | Mechanism | Predicted best detector |
|----------------|-----------|-------------------------|
| **Indel** (A)    | Intervals appear/disappear (set membership churns) | Mode A (interval-level)    |
| **Jitter** (B)   | Interval positions drift ±100 bp (membership preserved) | Mode B (point-level)       |
| **Combined** (C) | 50/50 mix of Indel and Jitter changes               | Mode C (interval + points) |

**The hypothesis as stated did not survive.** `RESULTS.md` is the findings
record and supersedes this table: the Mode A / Mode B asymmetry on Indel was
an artifact of genome saturation, and on a sparse genome both modes track
Indel equally. Read `RESULTS.md` before citing anything from here.

## Layout

```
experiments/synthetic_evolution/
├── Makefile          # drives data → evolve → hammock → analyze
├── RESULTS.md        # findings record (read this, not the table above)
├── code/
│   ├── createBED.py  # base BED generator (ported from original)
│   ├── evolve.py     # generation-stepped evolution simulator
│   └── analyze.R     # CSVs → figures + summary
├── data/      → /vast/…/synthetic_evolution/data
│   ├── synthetic_n100000_g{RANGE}.bed   # base BED (seed-fixed)
│   ├── synthetic_base_g{RANGE}.txt      # 1-line path file (bed1 for hammock)
│   └── evolved_g{RANGE}/                # 20 snapshot BEDs per (rate, type)
│       ├── synthetic_g100_r0.05_A_005.bed … _100.bed
│       └── list_r0.05_A.txt             # path file (bed2 for hammock)
├── results/   → /vast/…/synthetic_evolution/results
│   ├── hammock_g{RANGE}/    # per-(rate, type, mode[, expA][, subB]) CSVs
│   └── summary_g{RANGE}.csv # Jaccard at gen=50, one row per (rate, type, mode, expA, subB)
├── logs/      → /vast/…/synthetic_evolution/logs
└── figures/g{RANGE}/                # PNGs (committed)
    ├── jaccard_by_generation.png    # headline grid
    ├── sensitivity_heatmap.png      # Jaccard at gen=50 by mode × type × rate
    ├── mode_per_type_r{RATE}.png    # per-rate 3-panel headline (one per rate)
    ├── mode_c_expA_sweep.png        # Mode C across the expA sweep
    └── mode_c_expA_interpolation.png # expA interpolating Mode C between A and B
```

Every artifact is tagged `g{GENOME_RANGE}` so the two regimes coexist; both
`g1000000` (saturated) and `g100000000` (sparse) are on disk. See
`RESULTS.md` for what the two tags mean.

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

| Variable       | Default                       | Note                      |
|----------------|-------------------------------|---------------------------|
| `RATES`        | `0.005 0.01 0.05 0.1`         | Per-generation change rate |
| `TYPES`        | `A B C`                       | Evolution types            |
| `SUBB`         | `0.25`                        | Mode C subB — pinned, not swept |
| `C_EXPA_POS`   | `0.50 1.00 2.00 3.00`         | Mode C expA sweep (expA=0 run separately) |
| `GENOME_RANGE` | `100000000`                   | Max start per chr; sets the `g{RANGE}` tag |
| `N_INTERVALS`  | `100000`                      | Base BED size              |
| `GENERATIONS`  | `100`                         | Generations to evolve      |
| `STEP`         | `5`                           | Snapshot every N gens      |
| `PRECISION`    | `18`                          | HLL precision              |
| `THREADS`      | `8`                           | Per-hammock-invocation     |
| `SEED`         | `42`                          | Fixed for reproducibility  |
| `HAMMOCK`      | `…/claude-ref-comparison/bin/hammock` | hammock_claude binary |

`make all` needs `Rscript` on PATH for the analysis step (`ml r/4.3.0` on the
cluster). `make help` prints the current values; `make list-targets` prints
every CSV the sweep will produce.

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
- ~~`--expA` removed (it was a placeholder-column multiplier and is gone from
  hammock_claude; see `CLAUDE.md`).~~ **Corrected 2026-08-06: wrong.** `--expA`
  is alive and this experiment sweeps it (`C_EXPA_POS`; see `RESULTS.md` "What
  `expA` is for"). What v0.4-era hammock dropped was the `** expA` *exponent
  applied to the containment column*, not the flag — `expA` still changes
  interval multiplicity upstream of the sketch (`CLAUDE.md` divergence #2).
- Mode C subB is **pinned at 0.25**, not swept (it was swept in an earlier
  pass; the sweep axis is now `expA`). The original swept subB only at the
  saturated `r=1.0`.
- HLL containment/cosketch columns are emitted by hammock_claude but ignored
  by this analysis — the headline plot uses `jaccard_similarity` only.

  **Caveat added 2026-08-06.** `jaccard_similarity` is the register-equality
  estimator, which carries a chance-agreement floor and is *not* set Jaccard
  (`CLAUDE.md` divergence #2). It inflates low values, which is exactly the
  regime the headline decay curves end in. Worked example from
  `results/hammock_g100000000/synth_r0.05_A_hll_p18_jaccB.csv`, gen=100:
  `jaccard_similarity = 0.111`, but the inclusion-exclusion value recovered
  from the same row's containments (`J = 1/(1/C_AB + 1/C_BA − 1)`, C_AB
  0.0522 / C_BA 0.4101) is **0.049**. The qualitative decay story is
  unaffected; any absolute low-J number quoted off these CSVs is not.
  These CSVs predate the `jaccard_similarity_ie` column, so it must be
  derived from the containments rather than read directly.
