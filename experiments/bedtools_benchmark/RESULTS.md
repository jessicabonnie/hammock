# bedtools_benchmark — May 2026 results

Two runs against the May-10 optimized `hammock-cpp` binary
(`/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp`),
on Rockfish `shared` partition (sr-class nodes), 16 cpus, 32 GB:

- **2026-05-10 (morning)**: five sweeps with **sequential** pre-sort.
- **2026-05-10 (evening) → 05-11**: same five sweeps with
  **parallel** pre-sort (`ThreadPoolExecutor` sized to `num_threads`).
  Used for the canonical tables below — this is the workflow a real
  HPC user would run.

Headline: the optimized hammock is **~3× faster across the board** vs
the prior build, with **same accuracy properties** and **same scaling
shape** as the orig October 2025 run. The crossover where hammock
beats bedtools sits at **N=128 files** under either sort framing —
parallelizing the sort step shrinks the sort cost ~15× and makes the
"include sort in workflow" comparison effectively the same as the
"pre-sorted inputs" comparison.

## Source files

CSVs/text reports live in `results/` (symlinked to
`/vast/blangme2/jbonnie/hammock_claude_experiments/bedtools_benchmark/results/`).
PNGs in `figures/`.

**Parallel-sort run (canonical):**

| Sweep         | Job ID      | CSV / report stem                            |
| ------------- | ----------- | -------------------------------------------- |
| precision     | 23923204    | `sweep_precision_20260510_213254`            |
| threads       | 23923205    | `sweep_threads_20260510_213254`              |
| intervals     | 23923206    | `sweep_intervals_20260510_213254`            |
| files (t=8)   | 23923207    | `cpp_vs_bedtools_t8_20260510_213254`         |
| files (t=16)  | 23923208    | `cpp_vs_bedtools_t16_20260510_221418`        |

**Sequential-sort run (earlier, kept for reference):**

| Sweep         | Job ID      | CSV / report stem                            |
| ------------- | ----------- | -------------------------------------------- |
| precision     | 23902089    | `sweep_precision_20260510_172628`            |
| threads       | 23902090    | `sweep_threads_20260510_172628`              |
| intervals     | 23902091    | `sweep_intervals_20260510_172628`            |
| files (t=8)   | 23902092    | `cpp_vs_bedtools_t8_20260510_172628`         |
| files (t=16)  | 23902102    | `cpp_vs_bedtools_t16_20260510_181334`        |

The precision sweep also dumps `*_pairs.csv` (per-pair bedtools and
hammock jaccards across all 4096 pairs × 5 precisions × 3 runs =
61,440 rows) for the scatter plot in `figures/`.

## Optimization speedup

Hammock-cpp wall time, before vs after the May-10 build (same inputs):

| Sweep, config                       | Before  | After   | Speedup |
| ----------------------------------- | ------- | ------- | ------- |
| precision p=14, 64 files × 10k      | 82.9 s  | 26.3 s  | 3.15×   |
| threads t=1, 64 files × 10k         | 637.1 s | 198.9 s | 3.20×   |
| threads t=16, 64 files × 10k        | 42.7 s  | 13.7 s  | 3.12×   |
| intervals 1M × 16 files             | 1999 s  | 644 s   | 3.10×   |

Near-constant 3.1–3.2× speedup across configs.

## files sweep at t=16 (the headline comparison)

Parallel-sort run, p=14, 10k intervals/file, N up to 512, 3 runs/config:

| N   | bedtools wall | hammock wall | speedup | sort wall (parallel, 16w) |
| --- | ------------- | ------------ | ------- | ------------------------- |
|   2 |   0.76 s |   0.45 s | 1.70× | 0.10 s |
|   4 |   0.51 s |   0.87 s | 0.58× | 0.09 s |
|   8 |   0.63 s |   1.74 s | 0.36× | 0.11 s |
|  16 |   1.18 s |   3.49 s | 0.34× | 0.22 s |
|  32 |   3.23 s |   6.89 s | 0.47× | 0.41 s |
|  64 |  11.51 s |  14.06 s | 0.82× | 0.74 s |
| 128 |  44.89 s |  27.49 s | **1.63×** | 1.43 s |
| 256 | 184.47 s |  54.99 s | **3.35×** | 2.79 s |
| 512 | 767.51 s | 110.07 s | **6.97×** | 5.52 s |

Crossover at **N=128** (1.63× there, 6.97× at N=512). The orig Oct 2025
run reported 1.45× at N=128 and 5.83× at N=512 — same shape, slightly
larger ratios in our run because the sr-class hardware was slower for
bedtools-via-parallel-fan-out than the orig's c-class hardware.

With parallel sort the bedtools workflow penalty essentially vanishes:

| N   | bt wall | sort | bt + sort | hammock | speedup (with sort) |
| --- | ------- | ---- | --------- | ------- | ------------------- |
|  64 |  11.5 s |  0.7 s |  12.3 s |  14.1 s | 0.87× |
| 128 |  44.9 s |  1.4 s |  46.3 s |  27.5 s | 1.68× |
| 256 | 184.5 s |  2.8 s | 187.3 s |  55.0 s | 3.40× |
| 512 | 767.5 s |  5.5 s | 773.0 s | 110.1 s | 7.02× |

The crossover holds at N=128 either way. See
`docs/bedtools-parallelism-caveat.md` for the fairness framing —
parallelizing sort makes that whole discussion mostly moot.

## threads sweep — scaling preserved

64 files × 10k intervals, p=14:

| threads | bedtools wall | hammock wall | hammock scaling | bedtools scaling |
| ------- | ------------- | ------------ | --------------- | ---------------- |
|  1 | 62.6 s | 198.9 s |  1.00× |  1.00× |
|  2 | 43.2 s | 101.3 s |  1.96× |  1.45× |
|  4 | 21.5 s |  51.3 s |  3.88× |  2.91× |
|  8 | 11.0 s |  26.4 s |  7.54× |  5.71× |
| 16 | 10.8 s |  13.7 s | 14.50× |  5.79× |

Hammock OpenMP scaling is near-linear (efficiency 0.91 at t=16);
bedtools-via-GNU-parallel maxes out around t=8 (the per-pair cost is
small enough that parallel-spawn overhead dominates). At N=64 the
crossover from threads alone doesn't appear — bedtools wins at every
thread count. The crossover lives in the N axis.

## intervals sweep — Mode B is hash-bound at huge files

16 files, p=14, t=8:

| N intervals/file | bedtools wall | hammock wall | hammock/bt |
| ---------------- | ------------- | ------------ | ---------- |
|       1 000 |  1.4 s |   0.72 s | 0.51× |
|      10 000 |  1.2 s |   6.59 s | 5.5× |
|     100 000 |  4.1 s |  64.71 s | 15.8× |
|   1 000 000 | 32.7 s | 643.95 s | 19.7× |

Mode B hashes every base position, so cost is linear in total bp.
With only 16 files, bedtools' O(N²) is trivial (256 pairs) and the
hashing cost dominates. Hammock's win regime is "many files, modest
intervals/file", not "huge files." Worth noting that 1M intervals/file
isn't typical for real BED data — peak callers and DHS catalogs are
usually 10k–500k.

## precision sweep — HLL behaves textbook

64 files × 10k intervals, t=8:

| p  | hammock wall | RSS    | MAE vs bedtools | MAE vs hammock@p=18 |
| -- | ------------ | ------ | --------------- | ------------------- |
| 10 | 26.16 s |  5.6 MB | 0.1667 | 0.0110 |
| 12 | 26.20 s |  6.0 MB | 0.1663 | 0.0056 |
| 14 | 26.33 s |  7.6 MB | 0.1641 | 0.0026 |
| 16 | 27.64 s | 13.9 MB | 0.1639 | 0.0012 |
| 18 | 30.53 s | 39.0 MB | 0.1639 | 0.0000 |

- **MAE vs hammock@p=18** halves per +2 in p, matching the
  `1/√(2^p)` theoretical bound — clean HLL noise behavior.
- **MAE vs bedtools** is flat at ~0.164 across all p — this is the
  definitional gap (register-equality jaccard vs bp set-jaccard),
  *not* HLL noise. Explained in `docs/jaccard-definitional-gap.md`.
- Wall time grows slowly with p until p=18, where the register array
  size (2^18 × 1 byte = 256 KB per sketch × 128 sketches = 32 MB)
  starts to matter for cache.

## Next steps

User decision pending. Options:

1. **Make graphs.** The data is in CSV form across all five sweeps,
   plus per-pair scatter data for the precision sweep. The auto-emitted
   PNGs in `figures/` are working drafts — the next step is a paper-
   quality summary panel: crossover plot (N), parallel-scaling plot
   (threads), accuracy curve (precision, both ground truths), and the
   bp-jaccard vs register-equality scatter from `*_pairs.csv`.
2. **Make hammock even faster first**, then re-run all five sweeps.
   The optimization was a 3× wall-time win — another similar win
   would shift the N-crossover from 128 down to ~32 even on
   pre-sorted inputs, which is a much stronger story.

These aren't mutually exclusive but the order matters: if you
re-optimize first, the graphs you make will be the final ones.

## Caveats and known issues addressed

- **`find_hammock_cpp` was returning the wrong binary.** Old code
  did `sorted(...)[-1]` which lexicographically picked
  `cp38-cp38-…` (older) over `cp310-cp310-…`. Fixed to pick by
  `os.path.getmtime`. Stale `cp38-cp38-…` build directory deleted.
- **Sort time is now captured separately.** All CSVs from this run
  forward include a `sort_time` (sweep.py) or
  `mean_sort_time` (benchmark_cpp_vs_bedtools.py) column. Text
  reports show sort time per config.
- **Sort is now parallelized** (post-2026-05-10-evening commit). Both
  drivers fan the per-file `sort` subprocess calls across a
  `ThreadPoolExecutor` sized to the tool's thread count. The sort_time
  values in the tables above came from the earlier **sequential**
  implementation; future runs with the parallel sort will show
  proportionally smaller sort times (roughly num_threads-fold smaller
  at small N, less at huge N where I/O saturates).

## Related docs

- `docs/jaccard-definitional-gap.md` — why hammock jaccard ≠ bedtools jaccard
- `docs/bedtools-parallelism-caveat.md` — the GNU-parallel framing and the sort-time fairness section
