# CLAUDE.md

Guidance for Claude Code working in this repository. Internal notes — keep
divergences and unimplemented-feature info here, not in README.md.

## ⛔ Bedtools-relative performance numbers are retracted (2026-08-09)

Four harness defects, all inflating the bedtools baseline, all in our favour.
The largest: the bedtools module exports `LD_LIBRARY_PATH` with 17 directories
of which 4 are used, and `ld.so` searches the other 13 (on GPFS) on **every one
of N² execs** — worth 2.4–2.8×. N=512 bedtools 1978 s → 716 s; the headline
falls 27.61× → ~10.2×, itself provisional.

Also retracted: the "~123 exec/s process-creation cap" (measured 364 exec/s
clean), and the `md5sum` control that appeared to confirm it — that control ran
in the same polluted environment.

**Do not quote any speedup-over-bedtools number from this file or any other
until it is listed under "Re-measured and gated" in
`docs/bedtools-baseline-retraction.md`.** Hammock-internal ratios (mixed-stride
vs hash-threshold, precision sweeps, the fused-pass A/B, Mode D threading) are
unaffected — they divide two hammock numbers from one run.

The recurring failure was **verifying mechanisms without verifying magnitudes**:
a 0.7 s floor was accepted as the explanation for a 1262 s gap. Check a fix
against the size of the thing it claims to explain.

## What this is

A Python + C++ refactor of [`hammock`](https://github.com/jessicabonnie/hammock).
Same CLI; HyperLogLog backing for all four modes (A/B/C BED-interval, D FASTA);
parity-tested byte-equal CSV output against the orig on tested paths.

## Build / test

```bash
pip install -e . --no-build-isolation   # builds the C++ extension
pytest tests/
```

**Cluster compiler caveat (Mode D depends on it).** Build the extension with the
**conda env's own compiler**, not a spack `gcc/9.3.0` env-module. A loaded gcc
module sets `CC`/`CXX` to a spack gcc that bakes its own (old) `libstdc++` dir
into `_core`'s DT_RPATH *first*; that old libstdc++ loads before the conda one
and breaks `import digest` (`GLIBCXX_… not found`) → Mode D silently degrades
(now a hard error, see divergence #6). CMakeLists adds an `$ORIGIN`-relative
DT_RPATH to the env lib, but it only wins if the spack paths aren't injected
ahead of it — so build like:

```bash
CR=$CONDA_PREFIX   # the refactor env (claude-ref-comparison)
rm -rf build/      # wipe cached CMakeCache when changing compiler
CC=$CR/bin/x86_64-conda-linux-gnu-gcc CXX=$CR/bin/x86_64-conda-linux-gnu-g++ \
  pip install -e . --no-build-isolation
# verify: readelf -d <_core.so> | grep RPATH  → $ORIGIN/../../.. first, no gcc-9.3.0
#         python -c "from hammock.modes.sequence import _DIGEST_AVAILABLE as d; print(d)"  → True
```

See `memory/project_modeD_zero_rpath_digest.md`.

**Re-validate before changing floating-point build flags.** The Ertl
estimator feeds every cardinality-derived column (`containment_*`,
`cosketch_*`, `jaccard_similarity_ie`), and its τ/σ series terminate on exact
float equality, so a reassociating flag could in principle shift them.
Reassuringly, a spot check found bit-identical output across `-O0/-O2/-O3`,
every `-march=` tried (x86-64, v2, haswell, skylake-avx512, native), **and**
`-Ofast`/`-ffast-math` — the τ/σ loops are loop-carried, so there is nothing
to reassociate. Treat that as "no known hazard", not a guarantee: if you touch
the flags, re-run the byte-identity gate (capture a Mode B + Mode D CSV before
and after and diff the whole row). `jaccard_similarity` itself is an exact
ratio of two integer register counts and is immune either way.

The bed2fasta tests (`tests/test_bed2fasta*.py`) and their `--ref` end-to-end
paths need `bedtools` (and `samtools` for indexing) on `PATH` — `ml bedtools
samtools` on the cluster; they self-skip otherwise. Mode D parity needs the
conda-orig env (see Parity environments).

The wheel includes a standalone `hammock-cpp` binary built from the same
`hammock_core` static lib (in `build/`); intended for max-speed Mode B
benchmarking, no Python in the loop. The wheel does install it (to
`<site-packages>/bin/hammock-cpp`), but that directory is **not** on `$PATH`,
so invoke it by full path — `build/<wheel-tag>/hammock-cpp` after a local
build, or the site-packages copy.

**Since v0.7.0 it emits the full 9-column block by default** — `query`,
`reference`, `jaccard_similarity`, `jaccard_similarity_ie`,
`containment_AB/BA`, `cosketch_*` — matching the Python CLI **bit-for-bit** (the
IE derivation is written the same way in both — see
`jaccard_ie_from_containments` in `hammock_cli.cpp` and
`runner._jaccard_ie_from_containments`; keep them in sync or
`tests/test_hammock_cpp_metrics.py` fails on `==`). `--metrics` is still
accepted and is now a no-op.

**Pass `--no-metrics` for timing runs.** It drops back to the 3 columns and
tags the output `_j3`. The block costs one extra cardinality estimate per pair
(the union histogram is accumulated inside the Jaccard pass — see the fused-pass
note under Open seeds), so a timed run with it on is not comparable to the
numbers in `experiments/bedtools_benchmark/RESULTS.md`, which are all
`--no-metrics`. The benchmark harnesses pass the flag explicitly in both
directions, so the shape of a timed run no longer depends on a default.

**Measured cost of the metrics block, settled 2026-08-09** (job 29628907,
`docs/data/fusion_ab_20260808_232422.csv`; sr09, 16 reserved cores, N=64/side,
10k intervals, seeded corpus, 5 replicates, arm order permuted per replicate):

| p | 12 | 14 | 16 | 18 | 20 | 22 | 24 |
|---|---|---|---|---|---|---|---|
| block cost, `pair_time` (metrics ÷ no_metrics) | 1.79× | 1.62× | 1.61× | 1.59× | 1.64× | 1.67× | 1.67× |
| what the fusion bought (pre ÷ post, metrics arm) | 1.60× | 1.63× | 1.66× | 1.67× | 1.72× | 1.77× | 1.78× |

The **control** — the same two binaries on `--no-metrics`, which is
byte-identical code — reads 0.96–1.02, so ±2–4% is this experiment's
measurement error.

**Do not compare absolute timings across benchmark runs.** The earlier version
of this section reported a "flat ≈2.5×" multiplier and a 2.09× fusion speedup by
comparing a 2026-08-04 run to a 2026-08-08 one. The 08-04 run had no SLURM
allocation and was **1.27–1.53× slower on unchanged code**; the same control
above read 0.65–0.79 there. Within-run *ratios* survive that (they are paired);
absolutes do not. Full reasoning, and eleven still-open items:
`docs/seed-benchmark-methodology.md`.

Two things that section used to claim and should not be repeated. The
`comparison_time` multiplier is **not** the estimator's cost — at N=512 it is
now **72% serial `fprintf`** of six extra columns, precision-independent and
Θ(N²). And peak RSS cannot show the removed per-pair union allocation at all:
the per-thread sketching sketch is exactly co-sized with it and lives in a
different phase, so `max(RSS)` is identical either way. Older per-precision
tables: `docs/metrics-by-default.md` (still carries superseded numbers).

## Architecture

- `python/hammock/` — CLI (`cli.py`), orchestrator (`runner.py`), Mode D
  (`modes/sequence.py`). Threading via `concurrent.futures.ThreadPoolExecutor`;
  default `--threads = min(8, cpu_count())` for interval modes, **1 for Mode D**
  (`cli._default_threads`). The C++ extension releases the GIL, so
  interval-mode threading is real parallelism.

  **There are three thread budgets, not two** (`cli.main`):

  | budget | drives | Mode D clamp? |
  |---|---|---|
  | `threads` | the sketching `ThreadPoolExecutor` | yes — 1, GIL convoy |
  | `io_threads` | `bedtools getfasta` extraction | no |
  | `omp_threads` | the C++ OpenMP **pairwise** phase | no |

  `omp_threads` must not inherit Mode D's clamp: that clamp is about the GIL
  while *sketching*, whereas the pairwise loop runs once from the main thread
  with the GIL released, so clamping it to 1 would only make Mode D slower. It
  is `0` ("leave OpenMP's own default") unless `--threads` was given explicitly,
  so the no-flag path is unchanged. Before this existed the pairwise phase
  ignored `--threads` outright — `omp_set_num_threads` was called only by the
  standalone binary — so `hammock --threads 4` in a 4-CPU cgroup still spawned a
  team per core on the node. Applied as a `num_threads()` clause, never
  `omp_set_num_threads`, so it cannot leak into the sketching regions, which are
  already parallel *and* nested inside the pool (see the seed note below).

  **Mode D threading is not — it is a GIL convoy.** Measured 2026-08-03 on a
  48-core node, 4 Maurano DHS FASTAs,
  k=8 w=40 p=20: `--threads 1` = 31.9 s vs `--threads 8` = **68.3 s (2.1×
  slower)**. In-process, per-task time degrades monotonically (T1 1.00×,
  T2 0.61×, T4 0.86×, T8 0.40×) even with private per-thread sequence lists,
  so it is not refcount contention on shared data. `digest.window_minimizer`
  does *not* usefully release the GIL at this call pattern (~11 µs per call on
  a ~466 bp record), and the surrounding per-record Python work holds it
  outright. Eight concurrent *processes* on the same workload scale ~7.1×, so
  the machine is fine — it is the threading model. Corroborated by the archived
  sweep index (`maurano_dhs_validation/results/sweep_d_*.csv`), where every
  8-thread Mode D cell records `cpu_s/wall_s ≈ 0.8`.

  **Fixed in v0.6.1 by defaulting Mode D to `--threads 1`** (interval modes
  keep `min(8, cpu_count())`). An explicit `--threads > 1` is still honored, with
  a one-line stderr note. Re-verified 2026-08-04 on the same 4-FASTA workload:
  31.9 s default vs 63.1 s at `--threads 8`, CSVs byte-identical.

  Two consequences that outlive the default change:

  - **Every archived Mode D wall time is inflated** (they were all run at
    `--threads 8`). Re-baseline; don't compare new timings against them.
  - **Mode D still uses one core.** Making it *actually* parallel — a process
    pool — is worth more than any sketching optimization so far (cf. divergence
    #8: removing `_with_ends` bought 1.5–2.5× single-threaded, and processes
    scale ~7.1× here). That remains open; see the seed.

  bed2fasta extraction is exempt from the clamp: it shells out to `bedtools
  getfasta`, which parallelizes for real, so `runner._run_bed2fasta` reads a
  separate `args.io_threads` budget. Keep that split if you rework threading.

  Full evidence, option space, and reproduction: `docs/seed-mode-d-threading.md`.

## Open seeds (handoff notes for future work)

None is a decision; each is evidence gathered plus what still needs
establishing. Read the seed before re-litigating the question.

- `docs/seed-subsampling-synthetic-supplement.md` — the outline/draft's
  synthetic-corpus subB speedup claim ("4.21× at p=14") didn't trace to any
  CSV in the repo at any precision; the paired Maurano number ("1.83×") was
  real but from a superseded register-equality-era run (current is 1.87×,
  `docs/data/maurano_subB_ie_summary.csv`). Both fixed in the paper text.
  Still open: a p=18, `jaccard_similarity_ie`-era synthetic-corpus
  measurement to restore the corpus-dependence claim with real numbers.
- `docs/seed-hammock-cpp-file-dispatch.md` — two claims at different
  confidence levels. **Part 1 FIXED in v0.7.1 (2026-08-11)**: the Python
  CLI's `--threads` did not bound sketching-phase parallelism at all
  (verified by code read — dispositive, no thread parameter existed to carry
  the flag through — corroborated by a replicated isolation test where
  `--threads 2` used ~11.5 cores and CPU usage didn't track `--threads` at
  all). Fix: `sketch_bed_file_hll` gained a `threads` parameter threaded into
  `process_bed_file_mode_b/c`'s `num_threads()` clause (shared
  `cpp/include/hammock/omp_util.hpp::omp_team_size` helper, also used by the
  pairwise loops); `runner._sketch_many` (`_split_thread_budget`) divides
  `--threads` between its outer file-dispatch pool and each file's inner OMP
  team so the product stays bounded by the user's request.
  `hammock-cpp` was already correct (unaffected). **Part 2 not yet
  confirmed, and now needs re-measurement**: a related observation that the
  CLI was faster than `hammock-cpp` at N≥32 files/side, which would mean the
  headline benchmark figures (all timed on `hammock-cpp`) understate
  hammock's real throughput — direction was code-consistent but magnitude
  was never checked against the observed effect size, so don't treat it as
  settled or re-derive any figure from it. The Part 1 fix removes the
  sketch-phase oversubscription Part 2's candidate mechanism depended on, so
  the crossover must be re-measured post-fix, not assumed to still hold —
  see the seed doc's "Next steps."

- `docs/seed-benchmark-methodology.md` — **read this before quoting or
  generating any performance number.** Comparing one benchmark run to another in
  this repo is confounded by an unmeasured machine factor of up to 1.5×: the
  `--no-metrics` arm is byte-identical code across the fusion and it still
  measured 1.27–1.53× slower on the 2026-08-04 run, which had no SLURM
  allocation. Within-run *ratios* survive; cross-run absolutes do not. The seed
  carries the corrective paired/interleaved design (job 29628907), the reason
  the peak-RSS argument is null by construction, and the fact that the residual
  `+IE` cost is now 72% `fprintf` rather than estimator work.

  **The bedtools baseline is process-creation-bound (2026-08-09).** A pairwise
  `bedtools jaccard` workflow launches one process per pair — N² of them — and
  on these nodes process creation caps near **123 exec/s and does not scale with
  cores**, so "bedtools at t=16" can mean "bedtools at t≈1.5" and inflate a
  quoted speedup by up to ~6×. It is not a bedtools defect and not GNU
  parallel's: `md5sum` on node-local files measures **0.46×** at 16-way,
  `xargs -P16` hits the same ceiling, and moving the binary off GPFS to local
  NVMe changes nothing. Measured achieved efficiency swings **1.17×–2.86×**
  across nodes on identical code, and bedtools *regresses* past t=8 (3.30× at
  t=8 → 1.89× at t=48, wall rising 18.7→32.5 s) while hammock is monotonic to
  18.81× and keeps 43.3 of 48 cores busy.

  Consequences, all already applied: `bedtools.sh` runs **one** process per pair
  (`parallel --tagstring` + one `awk`) rather than three, worth ~2.1×; every
  bedtools CSV row now carries `mean_bedtools_serial_ms` and
  `mean_bedtools_parallel_eff` — **read the latter before quoting any speedup**,
  and only at N ≳ 32, since smaller cells are startup-dominated and read ~0.01.
  Full evidence and the superseded "~5.7× premium" table:
  `docs/bedtools-parallelism-caveat.md`. Figure:
  `paper/figures/threading_supplement.png`.

  **`bedtools` is pinned to 2.30.0**, and this is a correctness fix, not
  hygiene: 2.27.1 computes an **order-dependent union**, so 93 of the 190
  unordered Maurano pairs had J(A,B) ≠ J(B,A) in what is supposed to be the
  *exact* reference. Magnitude is small (max |Δ| 1.3×10⁻⁵, mean 5.6×10⁻⁷) so no
  archived conclusion changes, but it is 6.4% of the IE MAE at p=23 in the worst
  pair — do not use pre-2026-08-09 bedtools columns to argue about differences
  at the 10⁻⁵ level.

  **`docs/data/maurano_bedtools_ref.tsv` specifically is exempt from that
  warning, verified 2026-08-10, not just asserted.** Its actual provenance was
  unknown — it turned out to be a byte-identical carry-forward of a file dated
  2025-07-21 in the pre-refactor repo, with no record of which bedtools build
  produced it, and it had never been checked against the warning above despite
  being named in it. Regenerated it independently against the pinned 2.30.0
  binary (`bedtools jaccard` over all 400 ordered pairs on the same 20 Maurano
  BEDs) and diffed every column: **0/400 rows differ, 0/190 unique pairs
  asymmetric** — bit-identical to the checked-in file, `md5sum` included.
  Whatever produced the original file, it did not exhibit the 2.27.1
  order-dependent-union bug. So this specific file is fine to use at the 10⁻⁵
  level too; the warning above still applies to other pre-2026-08-09 bedtools
  columns whose provenance hasn't been checked the same way.

- `docs/seed-mode-d-threading.md` — the *default* is fixed (Mode D → 1 thread,
  v0.6.1), but Mode D still runs on one core. Still open: making it genuinely
  parallel via a process pool, which needs pickling on `HLLSketch` — a
  binding-only addition, since `hll_sketch.hpp` already exposes `registers()`
  publicly. **Closed 2026-08-08: interval-mode `_parallel_map` threading.** It
  had never been measured, and the guess that `min(8, cpu_count())` was leaving
  a 48-core node idle is **wrong** — `process_bed_file_mode_b/c` are themselves
  `#pragma omp parallel` (`cpp/src/processing_modes.cpp:134,241`), so sketching
  already saturates the machine from inside one C++ call. Measured on 24 BED
  files × 30k intervals, mode B p=16, medians of 3:

  | `--threads` | 1 | 2 | 4 | 8 | 16 | 32 | 48 |
  |---|---|---|---|---|---|---|---|
  | wall (s) | 5.52 | 4.22 | 3.71 | 3.73 | 3.81 | 3.78 | 3.76 |
  | cpu/wall | 34.8 | 39.4 | 43.8 | 44.7 | 44.0 | 44.4 | 44.6 |

  `cpu/wall ≈ 35` at **one** Python thread is the tell — 72% of a 48-core node
  busy before the pool contributes anything. The pool plateaus by ~4 threads and
  total CPU is flat past that (~165 s for 2–48 threads; the 1-thread cell is
  192 s, 16% higher), so the Python layer buys nothing beyond ~4 — raising the
  cap would only deepen the nesting. **Scope: one config** (24 files × 30k
  intervals, mode B, **p=16**, medians of 3, one node). It does not reach the
  p=24 regime where each nested thread holds a 16 MiB sketch, which is where
  oversubscription would bite hardest, so read it as "no headroom at p=16",
  not as a precision-independent result. Note the real
  consequence runs the other way: 8 pool threads × a per-core OpenMP team, each
  allocating a thread-local `HLLSketch` (16 MiB at p=24), is oversubscription,
  not headroom. Sizing *that* is the open question, not raising the cap.
- `docs/seed-mode-d-hash-width.md` — `digest` returns ≤32-bit minimizer hashes
  while the HLL assumes 64, biasing Mode D cardinality −0.5% to −8.3%. The
  `hash_size=32` option became viable once `_with_ends` was removed (divergence
  #8), since there is no longer a second sketch to merge with.
  **The pairwise metrics loop is one fused register pass** (2026-08-08).
  `HLLSketch::jaccard_and_union_cardinality` walks both register arrays once,
  accumulating the matching/active tallies *and* `counts[max(a[i], b[i])]++` —
  the union's register-value histogram — so no union sketch is ever built. It
  replaced `jaccard_similarity()` + `intersection_size()`, which between them
  walked the registers five times per pair and heap-allocated a whole sketch
  (16 MiB at p=24) N·M times, while re-estimating both operands' cardinalities
  that `pairwise_metrics_hll` had already hoisted.

  **Bit-identical by construction, not by measurement.**
  `ertl_improved_estimate` consumes registers *only* as the integer histogram it
  builds at the top, and the fused histogram is the same multiset the
  materialized union would produce — identical integers in, identical doubles
  out of the same τ/σ code. The one register-order-sensitive path is the
  `z < 1e-10` Flajolet fallback, reproduced over `max(a[i], b[i])` in index
  order. Both front ends use the helper, so
  `tests/test_hammock_cpp_metrics.py`'s cross-tool `==` still gates them.

  Measured (medians of 9, `OMP_NUM_THREADS=16`), Python `pairwise_metrics_hll`:

  | | p=20 | p=22 | p=24 |
  |---|---|---|---|
  | 144 pairs, before → after | 102 → 39 ms | 350 → 76 ms | 1320 → 274 ms |
  | 1024 pairs (256 at p=24) | 580 → 268 ms | 2070 → 486 ms | 2342 → 432 ms |

  i.e. **2.2–5.4×**, growing with precision. Part of that is the fusion and part
  is no longer deep-copying the sketch list across the pybind boundary (see
  below); they were measured together.

  **`pairwise_metrics_hll`/`pairwise_jaccard_hll` borrow their sketches.** They
  take `std::vector<const HLLSketch*>`, not `std::vector<HLLSketch>` — the
  by-value form made `pybind11/stl.h` deep-copy every register array on every
  call, serially and under the GIL (peak RSS 275 → 814 MiB for a 64-sketch
  p=22 pool passed as both operands; now 268 → 282). Two consequences:
  both arguments are `.noconvert()` and `require_uniform_precision` scans for
  nulls **before** its `a.empty() || b.empty()` short-circuit, because pybind's
  pointer caster turns `None` into a null that would be dereferenced with the
  GIL released inside an OpenMP region; and mutating a sketch from another
  Python thread *during* one of these calls is now a data race, where the copy
  used to make it a snapshot.

- `cpp/src/` — HLL sketch (parity-rewritten from scratch to match orig
  Python's algorithm exactly: low-bit register index, ctz rho, Ertl 2017
  improved estimator, register-equality Jaccard). Mode procs in
  `processing_modes.cpp`; mixed-stride in `stride.cpp`. xxh32 used for
  subsampling decisions (seed 31337); xxh64 for HLL ingestion (seed=42 by
  default, exposed as `--seed`).
- `bindings/_core.cpp` — pybind11 surface.
- `extern/hll/` — vendored upstream HLL library; we don't actually use its
  algorithm (kept for the optional fast path; current build links it for
  the standalone `hammock-cpp`).

## Intentional divergences from orig

These are deliberate; parity tests that touch them are skipped or projected.

1. **Mode B `--subB` honors the rate.** Orig's `intervals.py:542` silently
   ignored `subsample[1]` in Mode B (`point_subsample = subsample[1] if mode == "C" else 1.0`).
   We honor it. Any prior Mode B + `--subB` run on the orig is not
   parity-comparable.
2. **Containment + co-sketch columns.** Orig's single `containment` column
   was a placeholder. We replace it with a five-column block, computed
   from the **inclusion-exclusion** intersection — `|A| + |B| - |A ∪ B|`,
   Ertl estimator on each, union by register-wise max, clamped to `>= 0`
   (the `>= 0` clamp lives in `pairwise_metrics_hll`; `HLLSketch::intersection_size`
   is the equivalent scalar reference path, kept for tests). Same
   formula as orig's `hyperloglog.py estimate_intersection`. This is **not**
   the register-equality path that `jaccard_similarity` uses — see the
   estimator note below, it matters:
   - `containment_AB = |A ∩ B| / |A|`
   - `containment_BA = |A ∩ B| / |B|`
   - `cosketch_geom  = sqrt(C_AB · C_BA)`
   - `cosketch_arith = (C_AB + C_BA) / 2`
   - `cosketch_max   = max(C_AB, C_BA)`

   No `expA` exponent is applied — the orig's `** expA` semantics on a
   placeholder column were never load-bearing, and `expA` already changes
   what `A` *contains* (interval multiplicity) upstream of the intersection.

   Mode D emits this block once, on the minimizer HLL. (It used to emit a
   second `_with_ends` copy on a merged minimizer ∪ start-end HLL; that whole
   family was removed — see divergence #8.)

   Parity tests project these columns out before comparing.

   **The containment columns are emitted unclamped, so values just over 1.0
   reach the CSV** (`1.0000000000000007` occurs in ordinary runs), and
   `cosketch_max` inherits them. Only `jaccard_similarity_ie` clamps before
   use. Don't assert `<= 1.0` downstream. Magnitudes and mechanism:
   `runner._jaccard_ie_from_containments`'s docstring.

   **Read these CSVs back with `float_precision="round_trip"`.** Python writes
   `repr(float)` and `hammock-cpp` writes `%.17g` — the same double, different
   text, both exact. `pandas.read_csv`'s *default* float parser is not
   round-trip exact and will manufacture spurious ~1e-16 disagreements when you
   recompute `cosketch_geom`/`cosketch_arith` from a row's own containment
   columns, or diff the two front ends. `csv` + `float()` is also fine (what
   `tests/test_parity_against_original.py` uses). `cosketch_max` is immune,
   being a selection rather than arithmetic — if it is the only column that
   matches, you have hit this.

   **The two estimator families in a row are on different scales — do not
   mix them.** `jaccard_similarity` is register-equality, which carries a
   chance-agreement floor `c` (two registers tie when different elements
   have equal ρ). The containment/cosketch block is inclusion-exclusion and
   carries no such floor. Measured on *disjoint* inputs at p=16, n=2×10⁵:
   `jaccard_similarity = 0.168` while `containment_AB = 0.000`. So
   `jaccard_similarity` is *approximately* an affine transform `c + (1−c)·J`
   of set-Jaccard, while the containments estimate the true set quantities.
   Two qualifications, both measured, both load-bearing:

   - **The affine model is only approximate.** Against a `c + (1−c)·J` fit
     with c=0.180 the residual peaks at **+0.025 near J≈0.5**, is −0.010 at
     J=0, and is identically 0 at J=1 (any slope-(1−c) line passes through
     (1,1) by construction), so the error is largest in the *middle* of the
     range, not at the ends. That curvature is why a fitted intercept (0.180)
     exceeds the true disjoint floor (0.1699) — consistent, not a bug.
   - **`c` is set by the load factor λ = n/m *and* by the cardinality ratio
     |A|/|B|**, not by precision as such: 0.1699 as λ→∞ at equal cardinality,
     but 0.152 at ratio 2, 0.097 at ratio 5, 0.058 at ratio 10; and 0.045 at
     p=24 once m > n. It is **not** a constant you can subtract, and
     **`jaccard_similarity` is therefore not rank-faithful across pairs of
     differing cardinality ratio.** On the 20-sample Maurano corpus (ratio up
     to 2.2) it inverts bedtools' ordering for 2.5% of pairs (Kendall
     τ = 0.951). Rank only within comparable pairs.

   Consequence: **`jaccard_similarity_ie` now ships as its own column** (see
   divergence #7), and the same value is recoverable from `C_AB`/`C_BA` in any
   CSV that has them — `J = 1/(1/C_AB + 1/C_BA − 1)`, exact to ~2 ulp, not an
   approximation. **Caveat: CSVs written before 2026-05-14 have the orig
   `containment` placeholder (constant 1.0) instead of the block, so nothing
   is recoverable from them** — that includes all the archived interval-mode
   A/B/C output.

   **Pick by what you need.** `jaccard_similarity_ie` wins on *calibration*
   (MAE vs bedtools 5×10⁻⁴ at p=20 vs 0.15 for `jaccard_similarity`), and it
   is the only one of the two that is comparable across pairs of differing set
   size. Use it for magnitude and for any comparison spanning different set
   sizes; use `jaccard_similarity` only to rank pairs of comparable size.
   Caveat on `_ie`: it is censored at 0 by the `>= 0` clamp (25/90 pairs at
   p=12, all low-J — an exact 0.0 means "clamped or empty", never "measured
   zero") and is uninformative below J ≈ a few/√m.

   **The counter-claim that `jaccard_similarity` wins on *resolution* is
   retracted as stated.** The measurement behind it (error sd 0.0014 vs
   0.0024 at p=16, J<0.05) used a binned statistic that mixes within-bin
   true-J variance into the register-equality column only, and it was taken
   with the size ratio pinned near 1. It was replaced by a rank statistic,
   which is immune to that contamination because Kendall τ is invariant under
   any monotone transform of the estimator. What survives is narrower and
   precision-indexed:

   | J < 0.05, τ vs bedtools | p=12 | p=16 | p=20 | p=24 |
   |---|---|---|---|---|
   | `jaccard_similarity` | 0.335 | **0.658** | **0.905** | 0.907 |
   | `jaccard_similarity_ie` | 0.289 | 0.562 | 0.794 | **0.967** |

   Above J = 0.05 both reach τ = 1 by p=20 and there is nothing to choose.
   (Corrected 2026-08-06: this said "by p=16". In the J ≥ 0.2 stratum
   register-equality is τ = 0.9804 at p=16 — one discordant pair out of 102 —
   and only ties at p=20. It does not change the conclusion, since the cell
   turns on a single comparison, but p=16 was not where they converge.)
   So the register-equality advantage exists, but only for *ranking*, only
   below J ≈ 0.05, only at p ≤ 20, and only among pairs of comparable size —
   and one step of `-p` removes it. **The operative rule is: read
   `jaccard_similarity_ie`; if your corpus is low-J and you need ranking,
   raise `-p` to 24 rather than switching columns.** Generators:
   `experiments/bedtools_benchmark/estimator_rank_by_precision.py`,
   `paper/estimator_crossover/plot_estimator_crossover.R`. Caveats and the
   cross-species checks are in `docs/estimator-analysis-findings.md` §9.
3. **Default `--subB-method=mixed-stride`** — deterministic chr-keyed
   stride sampling. Orig's pipx-installed 0.4.0 didn't accept the flag
   at all (it lived only in WIP changes). We made mixed-stride the
   default in v.X because it is substantially faster than hash-threshold
   at every subB level and statistically well-behaved on genomic data;
   the structured-sampling concern is theoretical for BED inputs and the
   chr-keyed offset breaks cross-chromosome alignment. To get orig parity
   on subB runs, pass `--subB-method=hash-threshold` explicitly.
   Historical mixed-stride results from orig need to be re-run; see
   `memory/project_mixed_stride_rerun.md`.
4. **`--subB-method=single-hash`** is an opt-in parity divergence. Uses
   one `xxh64(point, seed=hll_seed)` to drive both the gate (high 32 bits
   compared to `subB * UINT32_MAX`) and HLL ingestion (full 64 bits). Orig
   uses `xxh32(point, seed=31337)` for the gate, distinct from the `xxh64`
   ingestion hash. Statistically equivalent estimator but a different
   *accepted-position set* per file, so CSV output is not byte-equal to
   orig hash-threshold. Prints a one-line stderr note when used.
   Note: single-hash is **not** automatically faster than hash-threshold —
   it trades one xxh32 (cheap) for one xxh64 (more work) per position.
   At low subB it's slightly slower than hash-threshold; its only real win
   is the **middle** of the range. Measured on Maurano
   (`docs/data/maurano_subB_summary.csv`, wall median, s):

   | subB | hash-threshold | single-hash | |
   |---|---|---|---|
   | 0.01 | 7.26 | 7.97 | slower |
   | 0.10 | 8.95 | 9.21 | slower |
   | 0.25 | 11.18 | 11.06 | 1.1% faster |
   | 0.50 | 14.67 | 13.62 | **7% faster** |
   | 1.00 | 9.555 | 9.546 | identical |

   At subB=1.0 the two are identical because the gate short-circuits and is
   never entered (`stride.cpp:31` — `do_subsample = subsample < 1.0` gates
   both methods), so there is no hash to trade.

   **Do not read the 0.50 row as a reason to use single-hash.** At subB 0.25
   and 0.50 *both* methods are slower than not subsampling at all (11.18 and
   14.67 s against `wall_nosub` 9.54 s — 0.65× and 0.70×). The 7% win is a
   win inside a regime that is already a net loss, so it never makes
   single-hash the right choice; it only shows the trade is not monotone.
   (Corrected 2026-08-06: this entry previously read "only at subB near 1.0
   is it marginally faster", which has the direction backwards.) Kept as a
   research/comparison flag.
5. **`--gate-seed N`** (default 31337) seeds the xxh32 gate hash and
   mixed-stride chr->stride hash. Default 31337 preserves orig contract
   in hash-threshold mode. Lets users sweep gate randomness
   independently of `--seed` (HLL ingestion).
6. **Mode D minimizer hash insertion uses `add_hash64` directly.** Orig's
   `lib/minimizer.py:118-121` has two branches: a FAST_HLL fast path that
   calls `_cpp_sketch.add_hash64(np.uint64(hash_val))` and a Python-fallback
   slow path that does `add_string(str(np.uint64(hash_val)))`. The slow
   path kmer-iterates over the *decimal digits* of the hash and silently
   drops any minimizer where `len(str(hash_val)) < k` — which is most of
   them at typical k≥10 on random-ACGT corpora (selector hashes are
   biased small). In the orig conda env (Py3.12 + bioconda digest),
   `_acceleration_type='Python'` so orig falls through to the buggy slow
   path and returns `jaccard_similarity=0` even self-vs-self on synthetic
   FASTAs. We use `_core.HLLSketch.add_hash64(hash_val)` unconditionally
   — matching the *intent* of orig's fast path. Consequently our Mode D
   `jaccard_similarity` differs from orig **even on `tiny.fa`** (e.g. 0.75
   vs orig's 0.7903 at k=5,w=20): the minimizer *sets* are identical (same
   `digest`), but orig's slow path hashes the *decimal digits* of each
   minimizer hash as k-mers while we ingest the raw 64-bit hash. So Mode D
   parity is **structural, not byte-equal** — `tests/test_mode_d_parity.py`
   asserts matching structural columns + well-formed similarity values, not
   equality. (An earlier byte-equal version only "passed" because
   `shutil.which("hammock")` resolved to the orig binary and self-compared;
   run the refactor env's `bin` first on PATH to actually exercise it.)
   The earlier "always 0" observation in memory turned out to be this bug;
   see `memory/project_modeD_no_ends_zero_bug.md`. For the separate
   silent-zero-from-broken-`digest` failure mode (RPATH shadowing
   libstdc++), see `memory/project_modeD_zero_rpath_digest.md`;
   `sketch_fasta` now raises loudly instead of falling back silently.
7. **Second Jaccard column: `jaccard_similarity_ie`** (v0.4.0). Orig emits one
   Jaccard column, computed by register equality. We emit that column
   unchanged — byte-equal, and parity tests still compare it — plus a second
   inclusion-exclusion column immediately after it. Rationale is divergence #2's
   estimator note: the register-equality column is not set Jaccard and is not
   rank-faithful across pairs of differing size, so a run had no
   bedtools-comparable Jaccard in it even though the underlying quantity was
   already being computed for the containment block.

   Derived Python-side in `runner._jaccard_ie_from_containments` from the
   `c_ab`/`c_ba` arrays `pairwise_metrics_hll` already returns — no C++ or
   binding change, mirroring how `_cosketch_from_containments` works.
   Containments are clamped to 1.0 first (Ertl noise can push them a few ulp
   past it), which forces `denom >= 1` and makes divide-by-zero and
   out-of-range results unreachable rather than merely unlikely. A zero
   containment scores `J_ie = 0.0`.

   `tests/test_parity_against_original.py` adds the new name to
   `_PROJECTED_OUT` so it is dropped before comparison — projecting out a
   *Jaccard* column looks wrong at a glance,
   but orig has no counterpart to be unfaithful to and `jaccard_similarity`
   is still compared byte-for-byte.
8. **The Mode D `_with_ends` column family is removed** (v0.6.0). Orig emits
   `jaccard_similarity_with_ends`; we no longer do, and we also dropped the six
   `_ie`/containment/cosketch twins we had added alongside it. Mode D now emits
   exactly one similarity block, on the minimizer HLL. `MinimizerSketch` is a
   single-HLL class: `startend_hll`, `merged()`, `canonicalize_kmer` and
   `_add_kmers_to` are gone.

   **Why.** The merged sketch was `minimizer_hll ∪ startend_hll`, where the
   ends payload was `canonicalize_kmer(s[:k] + s[-k:])` per record with *every*
   sliding k-mer of that 2k string inserted. Four measured problems:

   - **It was an uncontrollable blend, not a distinct quantity.** Because the
     two hash sets are disjoint, the merged Jaccard is exactly the size-weighted
     mediant of the minimizer and flank Jaccards: `J = (I_M+I_E)/(U_M+U_E)`.
     Verified against the shipped register-equality column to 2×10⁻⁴ at p=24
     (k=10, w=100: predicted 0.4095, actual 0.4097); the identity degrades to
     ±0.02 at p≤14 where the chance-agreement floor bends it. The mixing weight
     φ = |E|/(|M|+|E|) is set by (k, w, record count, length) and **never by the
     user** — on a Maurano file it is 0.55 at w=100 and 0.86 at w=500.
   - **k−1 of the k+1 inserted elements are chimeric.** Only indices 0 and k of
     the 2k payload are real sequence k-mers; the rest span the artificial
     `start||end` splice and exist in no genome. Measured on real peak FASTAs
     they are 77–79% of distinct end elements. Canonicalization was applied to
     the whole 2k concatenation, so a record's start-k-mer orientation depended
     on its *other* end.
   - **The flank component measures exact boundary identity, as a cliff.** At
     k=10, w=50, p=24, jittering record boundaries by **1 bp** takes the flank
     Jaccard from 1.000 to 0.130 (minimizer Jaccard: 1.000 → 0.997). On real
     data that bit is always off — two Maurano DHS BEDs share 15 exactly-matching
     intervals out of 168,322 (0.01%).
   - **The stated rationale was false.** The column was justified as rescuing
     records too short to yield minimizers. But the no-minimizer fallback fed
     the whole sequence into *both* HLLs and returned early, so on a
     constant-length sub-threshold corpus the two columns were **bit-identical**
     — it added nothing exactly where it was supposed to help. Dropout begins at
     `L < k + w − 1`.

   Empirically it also just lost: on the 235 existing Maurano Mode D configs,
   scored by Spearman ρ against bedtools, `no_ends` wins 185/49 (exact
   comparison). Filtering to cells that are unsaturated with <10% record
   dropout — an outcome-independent filter — `no_ends` wins 93.2%, and every
   exception is a cell where *both* columns have ρ < 0.39.

   **Caveats worth keeping.** The one place `_with_ends` won on the
   inclusion-exclusion family was 12 cells at w=500, ρ up to 0.9498 vs 0.8950 —
   but those have 72.5–73.1% record dropout, i.e. the regime where the column
   degenerates to record identity, and they sit well below the global optimum
   (k=20, w=20, p=24 → ρ_no_ie = 0.9997 at 1.9% dropout). Also note the arena
   is partly circular: Mode D FASTAs are `bedtools getfasta` extractions from
   one reference, so minimizer Jaccard nearly restates coordinate overlap
   (ρ = 0.9997). Full analysis and the adversarial review that produced it:
   `docs/mode-d-ends-removal.md`.

   Archived CSVs keep the columns; nothing on disk is invalidated. Analyses that
   consumed `jaccard_similarity_with_ends` were re-pointed to
   `jaccard_similarity` — see that doc for the list.
9. **`hammock-cpp`'s defaults now match the Python CLI's** (v0.7.0), and
   BagMinHash is gone. Three user-visible changes, best read as one sentence
   about the bare invocation: `hammock-cpp a.txt b.txt -o out` used to write
   `out_hll_p18_jaccA.csv` with 3 columns and now writes
   `out_hll_p18_jaccB.csv` with 9.
   - **Mode defaults to B**, not A, matching `_autodetect_mode` for BED input.
     Previously the two front-ends produced a column of the same name from
     different algorithms with no warning. Every in-repo caller passes `--mode`
     explicitly, so nothing archived is affected.
   - **The metrics block is on by default**; `--no-metrics` is the opt-out and
     tags its 3-column output `_j3`. The tag is on the *reduced* shape
     deliberately: that keeps filename ↔ column count one-to-one while leaving
     the default path where every pre-0.7.0 default run already wrote. Note the
     C++ suffix grammar therefore has a `_j3` component with **no Python
     counterpart** — do not "restore parity" by deleting it.
   - **`--peak-height` and `BagMinHashSketch` are deleted**, including
     `hammock.BagMinHashSketch`, `_core.BagMinHashSketch` and
     `_core.sketch_bed_file_bmh`. The flag selected the BagMinHash backend on
     the column index alone regardless of mode, so under the new mode default it
     would have silently discarded its own count weights and run serial — see
     the commit for the trace. Nothing in the repo used it.

   Two silent-overwrite bugs in `outprefix_with_suffix` were fixed at the same
   time: `--subA` never reached the filename (two Mode C runs differing only in
   `--subA` overwrote each other), and `subB` used an unconditional `%.2f` so
   `0.001` and `0.005` collided on `_B0.00`. Both now mirror `outprefix.py`,
   including its strict `< 0.01` boundary for the `.4f` form.

   `--version` exists on both front-ends now, sourced from `pyproject.toml` via
   `SKBUILD_PROJECT_VERSION`. The benchmark harnesses probe it on the resolved
   binary path and refuse anything older than 0.7.0 — they pass `--no-metrics`
   and parse microsecond timings, neither of which an older binary has.

## Mode D BED→FASTA (bed2fasta) — SHIPPED

Two BED lists, each tagged with a reference, are converted to FASTA and run
through Mode D. `python/hammock/refs.py` (reference resolution),
`python/hammock/bed2fasta.py` (extraction), wired into `runner._run_bed2fasta`.

- **Flags:** `--ref` (both lists) XOR `--ref1`/`--ref2` (per list);
  `--ref-cache-dir`, `--fasta-outdir`. Any reference flag forces Mode D and
  reinterprets the two positional lists as BED files. Combining a reference
  flag with `--mode A|B|C` is an error.
- **A reference is a keyword (`hg38`, `mm10`, `hg19`, `mm39`, `hs1`) or a local
  FASTA path.** `resolve_reference` **never downloads during a run** (HPC
  compute nodes are firewalled): a keyword must already be cached, else it
  raises with the exact `hammock fetch-ref <kw> --ref-cache-dir <dir>` command.
  `hammock fetch-ref` is a separate subcommand (run on a login node) that
  downloads + decompresses + indexes into the cache atomically (temp +
  `os.replace` + `.done` sentinel; http(s) only; size-capped).
- **Single extraction backend: `bedtools getfasta`** (default strand, headers
  `chrom:start-end`). `twoBitToFa` was deliberately *not* used — its `-bed`
  mode reverse-complements minus-strand intervals and headers by BED-name,
  which would silently diverge from getfasta and corrupt
  `jaccard_similarity_with_ends`. Needs `bedtools` on PATH (`ml bedtools`);
  `samtools` (`ml samtools`) is used to index references when available.
- **Silent-zero guard:** a chromosome-name mismatch (`chr1` vs `1`) makes
  bedtools warn + skip intervals but exit 0, yielding an empty FASTA →
  Jaccard≈0. `bed_to_fasta` turns an empty output or a "not found"/"WARNING"
  stderr into a hard `ConversionError`. It also warns on high-N extractions
  (assembly gaps → spurious shared minimizers) and record-count shortfalls.
- **CSV additions:** Mode D now emits **always-present trailing `ref1`/`ref2`
  columns** (`"NA"` for plain FASTA runs; matches the `num_hashes="NA"`
  convention). Parity tests project them out (`tests/test_mode_d_parity.py`
  `_PROJECTED_OUT`); the exact-header assertion in `tests/test_mode_d.py` was
  updated. A/B/C output and `_row_prefix` are untouched. The output filename is
  tagged `_<ref1tag>-vs-<ref2tag>` so cross-reference runs to the same `-o`
  don't overwrite each other.
- **CSV labels** are the original BED basenames (captured before the FASTA
  paths are swapped in), honoring `--full-paths`; generated FASTAs use unique
  per-index names so same-basename BEDs in different dirs never collide.
- **Cross-reference caveat:** cross-species Mode D Jaccard (e.g. hg38 vs mm10)
  measures **shared k-mer content — repeat/low-complexity-driven — not
  homology.** The `k=8` default suits within-reference use; prefer larger `k`
  for cross-species (cf. `memory/project_primate_phylogeny_substrate_engineering.md`,
  where deep topology failed without repeat-masking / larger k). Chromosome
  naming must match the reference; N-runs inflate similarity.

Deferred (not v1): in-run downloading / remote streaming (`twoBitToFa` over
HTTP fails on firewalled compute nodes), per-file references within a list
(`--ref1/--ref2` are per-list), and reference chr-prefix auto-fixup (cf.
`experiments/primate-phylogeny` `peak_chr_prefix`).

## Mode naming + default (user-facing)

The CLI presents two top-level modes — **interval** (BED) and **sequence**
(FASTA) — with A/B/C as deprioritized interval *flavors*. Names map to the
canonical single letters used everywhere internally (and in the CSV `mode`
column / output filename, which **keep the letter** for orig parity):

- `interval` / `interval-points` → **B** (base-level points; the default for BED)
- `interval-string` → **A** (exact interval strings)
- `interval-hybrid` → **C** (both; `--subA`/`--subB`/`--expA`)
- `sequence` → **D** (FASTA; the default for FASTA or any `--ref*` flag)

`--mode` accepts names or letters (`_MODE_ALIASES`/`_normalize_mode` in
`cli.py`, normalized to the letter at parse time). **The BED autodetect default
flipped A→B** (`_autodetect_mode`): plain BED → B; `--subA`/`--expA` →
C; `--subB` alone stays B. Parity tests pass explicit `--mode`, so they're
unaffected; `tests/test_autodetect.py` asserts the new default.

## Not implemented (would need work to add)

- Sketch types: `--minhash`, `--exact`, `--minimizer` for A/B/C. Phase 1
  shipped HLL only, and as of v0.7.0 there is no dormant C++ class behind any of
  them — `BagMinHashSketch` and its `--peak-height` flag were deleted (see
  divergence #9). HLL is the only `AbstractSketch` implementation now, so that
  interface is a single-implementation abstraction; re-deriving it when a second
  backend actually lands will be cheaper than carrying the current one, whose
  `jaccard_similarity(const AbstractSketch&)` + `dynamic_cast` shape is exactly
  what made BagMinHash awkward.
- File-level multiprocessing fallback (we use threads only).

## Parity environments

Two different orig installs because bioconda's `digest` is Python ≥ 3.9
while the orig pipx venv is on 3.8:

- `hammock-orig` (pipx, Py 3.8) → `tests/test_parity_against_original.py`,
  Modes A/B/C only.
- `/home/jbonnie1/.conda/envs/hammock/bin/hammock` (Py 3.12 + bioconda
  `digest`) → `tests/test_mode_d_parity.py`.

If a parity test is added for a new mode, prefer the conda env (modern
Python; has `digest`).

**Scope caveat: the Mode B byte-equality claim holds only for the list shape
the test happens to use.** `hammock-orig` is inconsistent *with itself* on
Mode B — for the same file pair at p=12 it emits `0.5478468899521531` when
both positional lists are `[tiny_a, tiny_b]`, but `0.5455607476635514` for
1×1, 1×2, and 2×1 shapes. Ours is shape-invariant and returns the first value
in every shape; an independent pure-Python reimplementation of the register
model (xxh64 seed 42 over the `chr\tpos` point strings) confirms ours is the
correct one, so **orig is the one that's wrong**. It matters because
`test_jaccard_byte_equal` passes the *same* two-file list for both positionals
(`_files_list`), which is exactly the one configuration where orig agrees.
So: do not "fix" a parity failure that appears after widening the test to
genuinely different query/ref lists — that is orig's bug surfacing, not a
drift in this repo. The mechanism was not pinned; orig's alternate value is in
the neighborhood of k-merizing the point strings at the default `kmer_size=8`
(closest probe 0.5453368), i.e. plausibly the same class of dual-path bug as
divergence #6. Established 2026-08-06 by end-to-end verification of every
emitted metric.

## CLI surface

Minimal — every flag does something. Removed during cleanup (see
`memory/`-tracked PRs for context):
`--debug`, `--python-only`, `--exact`, `--minhash`, `--minimizer`,
`--hyperloglog`, `--num_hashes/-n`. The `num_hashes` CSV column is hardcoded
to `"NA"` to preserve orig row layout.

## Conventions

- Python is the program. C++ extension is for hot paths only — no
  business logic in C++.
- Keep `_core` bindings narrow and stateless from the Python side. Sketches
  are passed back as Python objects with the heavy lifting done in C++ under
  GIL release.
- Parity tests over invented goldens, not the other way around. If we find
  a parity drift, fix our code or document the divergence here.
