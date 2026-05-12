# bedtools_benchmark — May 11 (evening) results

Three runs against `hammock-cpp`, on Rockfish `shared` partition (sr-class
nodes), 16 cpus, 32 GB:

- **2026-05-10 (morning)** — five sweeps with **sequential** pre-sort,
  pre-optimization hammock-cpp.
- **2026-05-10 (evening)** — same sweeps with **parallel** pre-sort and the
  first round of inner-loop optimizations (~3.1× speedup).
- **2026-05-11 (evening) — current canonical** — same sweeps with the
  second round of optimizations (inlined xxhash + `hash64_short`, removed
  `thread_local std::string`, incremental ASCII counter for the Mode B/C
  inner loop). Another ~2.1× wall-time win on top of the May-10 evening
  build, so the overall improvement vs the original May-10 morning code
  is ~6.5×.

Headline: with the latest hammock-cpp the **N-axis crossover moves from
N=128 (May-10 evening) to N=32 (where it's a 1.01× tie), with the first
decisive hammock win at N=64 (1.86×)**, and the asymptotic advantage at
N=512 grows from 6.97× to **13.18×**. Bedtools timings are
unchanged across runs (only hammock is moving). Plotting is now done in
R via `make_graphs.R` (ggplot2); the matplotlib `make_graphs.py` is
removed.

## Source files

CSVs/text reports live in `results/` (symlinked to
`/vast/blangme2/jbonnie/hammock_claude_experiments/bedtools_benchmark/results/`).
PNGs in `figures/`.

**Current run (canonical, May 11):**

| Sweep         | Job ID      | CSV / report stem                            |
| ------------- | ----------- | -------------------------------------------- |
| precision     | 23986411    | `sweep_precision_20260511_181918`            |
| threads       | 23986412    | `sweep_threads_20260511_181918`              |
| intervals     | 23986413    | `sweep_intervals_20260511_181919`            |
| files (t=8)   | 23986414    | `cpp_vs_bedtools_t8_20260511_181919`         |
| files (t=16)  | 23986415    | `cpp_vs_bedtools_t16_20260511_181919`        |

**Previous run (May 10 evening, parallel-sort baseline):**

| Sweep         | Job ID      | CSV / report stem                            |
| ------------- | ----------- | -------------------------------------------- |
| precision     | 23923204    | `sweep_precision_20260510_213254`            |
| threads       | 23923205    | `sweep_threads_20260510_213254`              |
| intervals     | 23923206    | `sweep_intervals_20260510_213254`            |
| files (t=8)   | 23923207    | `cpp_vs_bedtools_t8_20260510_213254`         |
| files (t=16)  | 23923208    | `cpp_vs_bedtools_t16_20260510_221418`        |

The precision sweep also dumps `*_pairs.csv` (per-pair bedtools and
hammock jaccards across all 4096 pairs × 5 precisions × 3 runs =
61,440 rows) for the scatter plot in `figures/`.

## Optimization speedup

Hammock-cpp wall time across the three rounds (same inputs):

| Sweep, config                       | May-10 AM | May-10 PM | **May-11 PM** | Total speedup |
| ----------------------------------- | --------- | --------- | ------------- | ------------- |
| precision p=14, 64 files × 10k      | 82.9 s    |  26.3 s   | **12.4 s**    | 6.7×          |
| threads t=1, 64 files × 10k         | 637.1 s   | 198.9 s   | **70.5 s**    | 9.0×          |
| threads t=16, 64 files × 10k        | 42.7 s    |  13.7 s   |  **6.6 s**    | 6.5×          |
| intervals 1M × 16 files             | 1999 s    | 644 s     | **301 s**     | 6.6×          |

This round's win came from three commits on May 11:
`hash64_short` inlining (`6cf430a`), removing `thread_local std::string`
in the stride loop (`9322087`), and incremental ASCII counter to avoid
a redundant `std::count` per interval (`cf796a4`).

## files sweep at t=16 (the headline comparison)

p=14, 10k intervals/file, N up to 512, 3 runs/config:

| N   | bedtools wall | hammock wall | speedup | sort wall (parallel, 16w) |
| --- | ------------- | ------------ | ------- | ------------------------- |
|   2 |   0.95 s |   0.22 s |  4.30× | 0.08 s |
|   4 |   0.53 s |   0.42 s |  1.28× | 0.08 s |
|   8 |   0.99 s |   0.83 s |  1.19× | 0.10 s |
|  16 |   1.45 s |   1.65 s |  0.88× | 0.20 s |
|  32 |   3.34 s |   3.29 s |  1.01× | 0.35 s |
|  64 |  12.23 s |   6.58 s | **1.86×** | 0.69 s |
| 128 |  48.78 s |  13.18 s | **3.70×** | 1.35 s |
| 256 | 200.62 s |  26.40 s | **7.60×** | 2.68 s |
| 512 | 694.26 s |  52.67 s | **13.18×** | 5.34 s |

Persistent crossover is now at **N=32** (1.01× — within run-to-run noise,
treat as a tie); the first **decisive** hammock win is at **N=64** (1.86×),
and the lead grows to 13.18× at N=512. The May-10 evening run had
crossover at N=128 (1.63×) and 6.97× at N=512 — the new optimizations
roughly halved hammock wall time without affecting bedtools, so the
crossover moves two notches leftward and the asymptotic ratio doubles.

The N=2 anomaly (4.30× speedup at N=2 then dropping back) is a bedtools
fixed-cost artifact: at N=2, bedtools spends most of its time in
GNU-parallel spin-up, which doesn't scale down. Hammock starts up
faster but its OpenMP threads don't fully amortize at N=2 either —
both lines bottom out around 0.5–1 s and the ratio is mostly noise.

## threads sweep — scaling preserved, hammock now faster than bedtools at t=16

64 files × 10k intervals, p=14:

| threads | bedtools wall | hammock wall | hammock scaling | bedtools scaling |
| ------- | ------------- | ------------ | --------------- | ---------------- |
|  1 | 57.9 s | 70.5 s |  1.00× | 1.00× |
|  2 | 38.2 s | 41.7 s |  1.69× | 1.51× |
|  4 | 19.2 s | 23.8 s |  2.97× | 3.01× |
|  8 | 10.1 s | 12.3 s |  5.72× | 5.72× |
| 16 |  9.8 s |  6.6 s | 10.74× | 5.91× |

Hammock OpenMP scaling stays near-linear (efficiency 0.67 at t=16, vs
0.91 in the prior round — the sketching kernel got fast enough that
OpenMP overhead and memory bandwidth show up); bedtools-via-GNU-parallel
still maxes out at t=8. **At t=16, hammock now beats bedtools at
N=64**, where in the prior round bedtools won at every thread count.
This brings the t=16 / N=64 cell (1.49×) into the win region for the
first time on this hardware.

## intervals sweep — Mode B is hash-bound at huge files

16 files, p=14, t=8:

| N intervals/file | bedtools wall | hammock wall | hammock/bt |
| ---------------- | ------------- | ------------ | ---------- |
|       1 000 |  1.4 s |   0.34 s | 0.25× |
|      10 000 |  1.6 s |   3.10 s | 1.94× |
|     100 000 |  5.2 s |  30.29 s | 5.83× |
|   1 000 000 | 34.3 s | 300.86 s | 8.76× |

Mode B hashes every base position, so cost is linear in total bp.
With only 16 files, bedtools' O(N²) is trivial (256 pairs) and the
hashing cost dominates. Hammock's win regime is "many files, modest
intervals/file," not "huge files." 1M intervals/file isn't typical for
real BED data (peak callers and DHS catalogs are usually 10k–500k), but
it's a useful upper bound: at 1M, hammock is now 8.76× *slower* than
bedtools (vs 19.7× slower in the May-10 evening run).

## precision sweep — HLL behaves textbook

64 files × 10k intervals, t=8:

| p  | hammock wall | RSS    | MAE vs bedtools | MAE vs hammock@p=18 |
| -- | ------------ | ------ | --------------- | ------------------- |
| 10 | 12.22 s |  5.6 MB | 0.1654 | 0.0110 |
| 12 | 12.25 s |  6.0 MB | 0.1659 | 0.0056 |
| 14 | 12.36 s |  7.6 MB | 0.1640 | 0.0026 |
| 16 | 13.42 s | 13.9 MB | 0.1638 | 0.0012 |
| 18 | 15.82 s | 39.0 MB | 0.1639 | 0.0000 |

- **MAE vs hammock@p=18** halves per +2 in p, matching the
  `1/√(2^p)` theoretical bound — clean HLL noise behavior.
- **MAE vs bedtools** is flat at ~0.164 across all p — this is the
  definitional gap (register-equality jaccard vs bp set-jaccard),
  *not* HLL noise. Explained in `docs/jaccard-definitional-gap.md`.
- Wall time grows slowly with p until p=18, where the register array
  size (2^18 × 1 byte = 256 KB per sketch × 128 sketches = 32 MB)
  starts to matter for cache. RSS is identical to the prior round —
  HLL register layout is unchanged.

## Plotting

R replaces matplotlib for figure generation:

```bash
ml r/4.3.0
Rscript experiments/bedtools_benchmark/make_graphs.R \
  --files-csv results/cpp_vs_bedtools_t16_20260511_181919.csv \
  --pairs-csv results/sweep_precision_20260511_181918_pairs.csv
```

Produces four PNGs in `figures/`:

- `*_sketch_compare_split.png` — bedtools vs hammock(sketch + compare),
  files axis, log-log, with a persistent-crossover annotation.
- `*_cost_per_pair.png` — wall_time / N² vs N.
- `*_jaccard_scatter.png` — per-pair (bedtools, hammock) jaccards
  colored by precision, zoomed to the data cluster.
- `*_jaccard_delta.png` — (hammock − bedtools) vs bedtools across all
  pairs, with per-precision running mean.

Output filenames mirror what the matplotlib version produced. Requires
`ggplot2`, `scales`, `dplyr`, `readr`, `tidyr`, `Cairo` (the system R
build at `r/4.3.0` lacks native PNG/cairo capability, so we render via
`CairoPNG`).

## Caveats and known issues addressed

- **`find_hammock_cpp` was returning the wrong binary.** Old code
  did `sorted(...)[-1]` which lexicographically picked
  `cp38-cp38-…` (older) over `cp310-cp310-…`. Fixed to pick by
  `os.path.getmtime`. Stale `cp38-cp38-…` build directory deleted.
- **Sort time is captured separately.** All CSVs from the May-10
  evening run forward include a `sort_time` (sweep.py) or
  `mean_sort_time` (benchmark_cpp_vs_bedtools.py) column. Text
  reports show sort time per config.
- **Sort is parallelized** (post-2026-05-10-evening commit). Both
  drivers fan the per-file `sort` subprocess calls across a
  `ThreadPoolExecutor` sized to the tool's thread count.

## Related docs

- `docs/jaccard-definitional-gap.md` — why hammock jaccard ≠ bedtools jaccard
- `docs/bedtools-parallelism-caveat.md` — the GNU-parallel framing and the sort-time fairness section
- `docs/mode-b-subsampling-perf.md` — design notes for the inner-loop
  optimizations that drove this round's speedup
