# Seed: CLI's `--threads` doesn't bound sketching parallelism, and a related hammock-cpp-vs-CLI wall-time crossover (found 2026-08-11, re-measured 2026-08-11)

This seed carries two claims at different confidence levels. Keep them
separate — the first is a confirmed bug, fixed. The second started as a
direction-only observation that didn't meet this repo's own bar for
"confirmed" (see `CLAUDE.md`'s retraction section: "verifying mechanisms
without verifying magnitudes was the recurring failure"); as of the "Part 2
update" section below it is measured, controlled (`--exclusive`), reproduced
on two node/CPU configurations, and turns out to be a **non-monotonic hump**
rather than the simple one-directional crossover first suspected — real in
the N≈8-1024 band, reversed by N=2048. Still job-level n=1 (one allocation
each) and one proposed mechanism (pairwise-phase dilution) is plausible but
not yet checked against its own magnitude — see "Still open" below.

## Part 1 — FIXED in v0.7.1 (2026-08-11): the Python CLI's `--threads` did not bound sketching-phase parallelism

**Fixed.** `sketch_bed_file_hll` (`bindings/_core.cpp`) gained a `threads`
parameter, passed through unresolved into `process_bed_file_mode_b/c`
(`cpp/src/processing_modes.cpp`), which now attach a `num_threads()` clause to
their `#pragma omp parallel` sketching regions via a shared
`omp_team_size()` helper (`cpp/include/hammock/omp_util.hpp` — hoisted out of
`bindings/_core.cpp`, which already had its own copy for the pairwise loops,
so both front ends' clamp formula can't drift apart). `runner._sketch_many`
computes the per-call split (`runner._split_thread_budget`): `outer_workers`
is the concurrency its own `ThreadPoolExecutor`-based file dispatch will
actually use (`min(threads, n)`, or 1 below the threshold where it falls back
to serial), and each dispatched file's inner OMP team gets
`max(1, total_threads // outer_workers)` — floor division, so the product is
always `<= total_threads` (never exceeds the budget; can under-use it by a few
threads when `n` doesn't divide evenly, a documented, accepted trade-off).
`hammock_cli.cpp`'s call sites are unchanged (they still pass no `threads`
arg, defaulting to 0/"ambient", which is correct since that binary already
calls `omp_set_num_threads(args.threads)` globally up front). Verified
empirically post-fix (temporary instrumentation, since reverted): for
`n=1` (single file per side, the isolation-test shape below), the OMP
region's actual `omp_get_num_threads()` matched the requested `--threads`
exactly at 1, 2, and 16; for `n=3, --threads 8` it correctly resolved to a
team size of 2 per file (`8 // 3`); a large synthetic workload (5M intervals)
showed CPU% now tracking `--threads` monotonically (144% → 183% → 272% at
threads 1/2/16) instead of the pre-fix ~11–14 cores regardless of the flag.
Full test suite (`pytest tests/`, including new `tests/test_sketch_threads.py`
covering the split arithmetic and the binding's `threads=` argument) passes,
244 passed / 8 skipped.

Note for anyone re-running the isolation script below post-fix: at the
script's default small workload (200k intervals, sub-second wall time), the
"Percent of CPU" metric is dominated by numpy/BLAS thread-pool spin-up at
Python import time, not by the sketch-phase OMP team — this confounds the
metric enough to mask the fix's effect entirely at that scale. Use a larger
`--num-intervals` (millions) or time the sketch phase in isolation
(`--verbose` on both front-ends) to see the bound take effect; don't conclude
the fix didn't work from a small-workload CPU% reading alone.

**Original bug report, preserved for context — this is bad independent of any
wall-time benefit it happened to produce:** a user who passes `--threads 4` —
e.g. because that's their actual
core allocation on a shared cluster node — should get *at most* 4 threads of
CPU usage, not an unbounded number decided by whatever OpenMP's ambient
default happens to be. Silently using far more than requested can starve
co-tenant jobs, get killed by a cgroup, or just make `--threads` a lie for
anyone trying to control resource usage. It also means **no CLI thread-scaling
benchmark can be trusted until this is fixed** — sweeping `--threads` and
expecting proportional behavior will not show what it appears to show, since
the sketching phase (the dominant cost at most N, see `CLAUDE.md`'s
metrics-block-cost notes) isn't actually responding to the flag.

**Verified two ways, at two different evidentiary strengths — kept separate
below rather than bundled under one "verified" claim, per the same
magnitude-checking standard Part 2 is held to. Pre-fix state, preserved for
context; line numbers below are as of the pre-fix code and no longer point at
the right spot post-fix (`sketch_bed_file_hll` moved and grew a `threads`
param) — see "Fixed" above for current locations:**

1. **The mechanism is proven by code inspection alone — dispositive, not
   statistical.** The pybind11 binding for sketching, `sketch_bed_file_hll`
   (`bindings/_core.cpp:354-367` pre-fix), took no thread argument at all — it
   only wrapped the call in `py::call_guard<py::gil_scoped_release>()`.
   `process_bed_file_mode_b`'s internal `#pragma omp parallel`
   (`cpp/src/processing_modes.cpp:134`) has no `num_threads()` clause, so it
   always uses OpenMP's ambient default team size. Nothing on the Python side
   ever calls `omp_set_num_threads` or sets `OMP_NUM_THREADS` (checked: no
   hits in `cli.py`/`runner.py`), so `--threads` has no code path — not even
   a side channel — by which it could reach this region. Contrast with the
   pairwise-comparison phase, which *does* get a `num_threads()` clause fed
   from `omp_threads` (`_core.cpp:214/265/274/283`) per the existing
   three-thread-budget table in `CLAUDE.md`'s Architecture section. The fix
   applied there (after a documented past bug: "the pairwise phase ignored
   `--threads` outright ... a team per core on the node") was never applied to
   the sketching phase. An absent function parameter cannot be a matter of
   run-to-run noise — this claim doesn't need replication the way an
   empirical measurement would.
2. **Empirical corroboration — illustrative, not the load-bearing evidence.**
   Single file per side, so no file-dispatch (`ThreadPoolExecutor`)
   parallelism is possible for either tool — isolates the sketching-phase OMP
   team specifically. 3 replicates, `devlangmead1`, `/usr/bin/time -v`,
   median `Percent of CPU`; data:
   `experiments/bedtools_benchmark/results/threads_isolation_1786475597.csv`
   (script: `measure_threads_isolation.py`; see "Reproducing" below).

   | `--threads` | hammock-cpp (median of 3) | hammock CLI (median of 3) |
   |---|---|---|
   | 2  | 167%  | **1154%** |
   | 4  | 269%  | **1565%** |
   | 16 | 647%  | **1459%** |

   `hammock-cpp` scales with the requested thread count (it calls
   `omp_set_num_threads(args.threads)` once, globally,
   `cpp/app/hammock_cli.cpp:304`, before its serial per-file loop — see Part 2
   — so every subsequent parallel region inherits that bound with no
   concurrent teams to conflict with it). The CLI does not: even at
   `--threads 2` it used ~11.5 cores, and CPU% doesn't track `--threads` in
   any consistent direction (it's *higher* at 4 than at 16) — consistent with
   "ignores the flag, uses whatever the node happens to give it," not with
   any partial/soft honoring of it. This table is a single-run
   corroboration of a claim already settled by (1); treat the exact
   percentages as illustrative of scale (an order of magnitude beyond what
   was requested), not as precise measurements.

**Fix direction — DONE, see "Fixed" above.** Chose the `sketch_bed_file_hll`
thread-count-parameter route (not the `omp_set_num_threads` route), and did
resolve the outer-pool interaction called out below: rather than naively
clamping the inner team to the full `threads` value while the outer pool is
*also* sized to `threads` (which would have reproduced the old
unbounded-oversubscription bug, just at the file level instead of the OMP
level), `runner._split_thread_budget` divides the budget between the two so
their product is bounded.

## Part 2 — LESS CERTAIN: does this (or something related) explain a hammock-cpp-vs-CLI wall-time crossover?

**Question, and why it isn't answered yet.** Every headline benchmark figure
in this paper (Figure 3A, S6, S7, S8) times the standalone `hammock-cpp`
binary, on the assumption that it is hammock's speed ceiling — same core
(`hammock_core`), minus Python overhead, so strictly faster than the CLI. A
side-benchmark measuring "how much slower is the CLI" instead found the CLI
**faster** at N≥32 files per side. **This section's status, per the review
that caught the gap: direction is plausible and code-consistent; magnitude
against the observed effect size, and generalization to what the headline
figures actually measure, are NOT yet verified. Do not treat this as settled.**

### The measurement

`experiments/bedtools_benchmark/measure_cli_overhead.py` /
`sbatch_cli_overhead.sh`, job 29758101. Median of 5 replicates,
`--mode B -p 18 --threads 16`, synthetic corpus, node sr15
(`--partition=shared`, not exclusive). Data:
`experiments/bedtools_benchmark/results/cli_overhead_1786473951.csv`:

| N (files/side) | hammock-cpp | hammock CLI | cli/cpp |
|---|---|---|---|
| 2   | 0.2797s | 0.4506s† | 1.611 (CLI slower, as expected) |
| 8   | 1.0914s | 1.1370s | 1.042 (CLI slower) |
| 32  | 4.3527s | 3.7613s | **0.864 (CLI faster)** |
| 128 | 17.6513s| 14.2486s| **0.807 (CLI faster)** |

† One of the 5 N=2 CLI replicates was a 5-8× outlier (2.4051s vs 0.43-0.46s
for the other 4) — plausibly a cold-cache first-process-in-job effect. The
median is nominally robust to it, but n=5 is thin enough that this single
point still shapes the reported ratio; don't treat 1.611 as precise.

A "fixed Python startup cost" model predicts the ratio approaches 1.0 from
above as N grows — never crosses under it. It crosses under it at N=32.

**Important scope note, corrected after review**: this was necessarily
measured entirely on hammock-cpp's `--metrics` arm, because the CLI has no
`--metrics`/`--no-metrics` opt-out (it always emits the full 9-column block).
`benchmark_cpp_vs_bedtools.py` — the harness behind every headline
figure — passes `--no-metrics` for its hammock-cpp timing arm. So **the arm
this crossover was measured on is not the arm the headline figures use**, and
the metrics block's own cost is precision- and N-dependent (`CLAUDE.md`:
~1.6-1.8× at N=64, "72% `fprintf`" at N=512), meaning the sketch-phase-vs-
pairwise-phase balance likely differs between the two arms in a way that
hasn't been characterized. This is the reverse of how an earlier draft of
this doc stated the gap (it read as if the metrics arm was the *untested*
extension; it is actually the *only* arm tested, and `--no-metrics` — what
the figures use — is untested).

### Candidate mechanism (code-read, internally consistent, magnitude NOT checked against the observed 13-19% effect)

- `hammock-cpp` sketches its N query files, then its N ref files, in a
  **strict serial `for` loop** (`cpp/app/hammock_cli.cpp:334-343`) — one file
  at a time, each getting a dedicated OMP team (bounded by `--threads`, see
  Part 1). Zero cross-file overlap, including the plain single-threaded
  `read_intervals()` (`processing_modes.cpp:117`) that runs before each
  file's parallel region.
- The CLI's `_sketch_many`/`_parallel_map` (`runner.py:248-291`) dispatches up
  to `--threads` files **concurrently** via `ThreadPoolExecutor`, each worker
  spawning its own unclamped OMP team (Part 1's bug). This is a genuinely
  oversubscribed pattern that nonetheless produced lower wall time at N≥32.
- Pairwise-phase and I/O-format differences were checked and are not
  first-order: both front-ends call the literal same function for this phase
  (`HLLSketch::jaccard_and_union_cardinality`, `hll_sketch.cpp:173`, invoked
  from both `hammock_cli.cpp:417` and `_core.cpp:296`), reaching the same
  effective thread count by different mechanisms — `hammock-cpp`'s pairwise
  loop (`hammock_cli.cpp:405`) has no `num_threads()` clause and relies on the
  ambient team size set once globally by its earlier `omp_set_num_threads()`
  call, while the Python bindings' `pairwise_metrics_hll`/`pairwise_jaccard_hll`
  attach an explicit `num_threads(nt)` clause (`_core.cpp:214/265/274`). Same
  outcome, not "matching clauses." hammock-cpp's own `--verbose` breakdown
  puts the pairwise+write phase at ~1.1% of wall time even at N=32 (1024
  pairs) — independently re-derived and confirmed (~1.157%, and consistent
  with archived scaling data in
  `docs/data/pairwise_cost_by_precision_20260804_164807.csv`).
- **What this mechanism does NOT explain, by its own logic**: overlapping
  serial-I/O-then-parallel-compute phases across files should make *hidden
  idle time productive*, which predicts the CLI's total CPU-seconds should be
  equal to or higher than hammock-cpp's — not lower. An informal, uncommitted
  second-node measurement (see caveat below) reported the CLI's CPU-seconds as
  *lower* at N=32/128, which is not explained by the stated mechanism and was
  not chased down further (candidates not distinguished: OMP team
  setup/teardown cost paid 64 times serially vs. spread across worker
  threads; NUMA effects from 16-wide teams spawned across the node; malloc-
  arena/cache-locality differences from thread-local state). **No magnitude
  estimate was made for the proposed mechanism against the observed
  13-19% wall-time effect** — this is exactly the check `CLAUDE.md`'s
  retraction section demands and it has not been done here.

**Caveat on the "independent reproduction" claim.** A second run on a
different node was described as reproducing the direction and rough magnitude
(0.898/0.868 vs the table's 0.864/0.807). **That run has no artifact in this
repo — no CSV, no log, no script, no node name recorded** — unlike every
other quantitative claim in this repo's seed docs. The two runs' ratios
differ by 4-7.5% relative, which is the same order as the ±2-4% this repo's
own controlled, paired benchmark design (`docs/seed-benchmark-methodology.md`,
job 29628907) attributes to pure machine noise on one allocation, and well
within the 1.27-1.53× machine-factor range that same doc documents for
*uncontrolled* shared-partition comparisons — which is exactly what both of
these runs are (neither used `--exclusive`). Treat the second run as an
informal, unreproducible spot-check that didn't contradict the first, not as
independent confirmation of magnitude.

### What would disconfirm this (falsification criteria, not just extensions)

- **If `--exclusive` brings the N=32/128 ratios to within ±2-4% of 1.0**
  (this repo's own controlled-benchmark noise floor, job 29628907 via
  `docs/seed-benchmark-methodology.md` — not the looser 1.27-1.53×
  uncontrolled-shared-partition figure, since `--exclusive` is the controlled
  case), the effect is shared-partition contention, not a real
  dispatch-strategy difference, and this part of the seed should be closed as
  a false lead.
- **If a `perf stat`/`strace -f` pass shows OMP team setup/teardown cost is
  negligible relative to the CPU-seconds gap**, the teardown-cost candidate
  explanation is ruled out (not just untested).
- **If timing the sketch phase alone (isolated from pairwise/write, e.g. via
  `--verbose` on both front-ends) doesn't reproduce the crossover**, the
  mechanism proposed above is wrong and something else (comparison phase,
  I/O, output writing) needs to be re-examined despite the earlier "~1.1% of
  wall time" ballpark.

### Consequence, if this holds up under the checks above

This would be the **opposite direction** from the bedtools-baseline
retraction (`CLAUDE.md`'s top section: harness bugs that inflated hammock's
apparent advantage). If hammock-cpp's serial-file dispatch really is
suboptimal at the N and on the arm (`--no-metrics`) the headline figures
actually use, those figures would be **understating** achievable hammock
throughput. **Nothing currently published is known to be wrong.** The effect
has only been measured on the `--metrics` arm, only up to N=128, only on
shared (non-exclusive) allocations, and its proposed mechanism has an
unexplained residual (see above). Do not re-derive or restate any existing
figure's numbers from this seed alone — it identifies a question, not a
correction.

## Part 2 update (2026-08-11): re-measured on --exclusive, extended to N=2048, both arms — job 29763124

**Next step 2, below, is done.** `sbatch_cli_overhead.sh` / `measure_cli_overhead.py`
now run two arms back-to-back in one `--exclusive` allocation: Arm B
(`--metrics`, N=2/8/32/128, otherwise reproducing job 29758101's config
exactly) to test the falsification criterion on the same arm it was
calibrated against, and Arm A (`--no-metrics`, N=2..2048, the arm the
headline figures actually use). Node c192 (Xeon Gold 6248R), 6 replicates
per cell, alternated tool order, checkpointed after every N block. Data:
`experiments/bedtools_benchmark/results/cli_overhead_metrics_exclusive_29763124.csv`,
`experiments/bedtools_benchmark/results/cli_overhead_nometrics_exclusive_29763124.csv`.

| N | Arm B (`--metrics`) cli/cpp | Arm A (`--no-metrics`) cli/cpp |
|---|---|---|
| 2 | 1.562 | 1.539 |
| 8 | 0.926 | 0.949 |
| 32 | 0.722 | 0.723 |
| 128 | 0.706 | 0.728 |
| 512 | -- | 0.859 |
| 1024 | -- | 0.959 |
| 2048 | -- | **1.127** |

**1. The "Important scope note" arm-mismatch gap is closed.** Arm A and Arm B
agree closely everywhere they overlap (N=2/8/32/128) -- the crossover is not
an artifact of the `--metrics` arm job 29758101 happened to use; it
reproduces on `--no-metrics`, the arm every headline figure actually runs.

**2. The falsification criterion FAILS to close Part 2 as a false lead.** The
criterion (above): if `--exclusive` brings N=32/128 to within ±2-4% of 1.0,
it's shared-partition contention. Instead both ratios moved *further* from
1.0 than job 29758101's original shared-partition numbers (0.864→0.722 at
N=32, 0.807→0.706 at N=128) -- the opposite of what the contention hypothesis
predicts. **One real, unresolved confound**: this job landed on a different
CPU model (c192, Xeon Gold 6248R) than job 29758101 (sr15, Xeon Gold 6448Y),
so the *magnitude* of the shift cannot be cleanly attributed to exclusivity
alone -- CLAUDE.md already documents `cpu_model` as load-bearing given
`-march=native`. But the *qualitative* result (a real crossover exists, and
now: it reverses at N=2048, see below) does not depend on resolving that
confound, since job 29758101 already established the same-direction effect
on a third node/CPU (sr15/6448Y, shared partition) and this job adds a
second, different node/CPU (c192/6248R, exclusive) -- two different physical
configurations agreeing on the qualitative shape is stronger than either
alone, even with the magnitude still uncertain.

**3. New finding the original seed did not anticipate: the effect is a hump,
not monotonic.** CLI is faster from N≈8 through N≈1024 (peak relative
advantage around N=32-128, ratio ~0.72), converges back to near-parity by
N=1024 (0.959), and **reverses by N=2048 -- hammock-cpp is ~13% faster
(ratio 1.127)**. N=2048 is the point closest to what the paper's real
headline figures (N up to 2048, see `sbatch_fig3_largeN.sh`) actually
measure. Candidate mechanism, not yet verified against the observed
magnitude the same way Part 2's original mechanism claim was flagged as
unverified: the pairwise-comparison phase is literal shared C++ code between
both front-ends (`HLLSketch::jaccard_and_union_cardinality`, already
parallelized identically via OpenMP on both sides -- confirmed by reading
`hammock_cli.cpp:405` and `_core.cpp:206/275`, both
`#pragma omp parallel for collapse(2)`), grows ~quadratically with N, and
costs both tools the same wall-clock regardless of dispatch strategy. As it
becomes a larger share of total wall time, it dilutes -- and past some N,
overwhelms -- a per-file dispatch-strategy difference that only lives in the
sketch phase. Consistent with, but not proof of: the seed doc's own earlier
note that pairwise+write was ~1.1% of wall at N=32 (1024 pairs) and grows
from there.

**4. Consequence, updated from the original "Nothing currently published is
known to be wrong."** Still true, and now on firmer ground specifically
*because* N=2048 reverses: the headline figures time hammock-cpp at
N=512/1024/2048, and at N=2048 -- the point measured here closest to that
regime -- hammock-cpp is faster than the CLI, not slower. So the original
worry (headline figures might be *understating* achievable throughput by
timing the "wrong" tool) does not hold at the N the figures actually use.
It does hold, newly and concretely, in the N≈8-1024 band -- real, reproduced
on two node/CPU configurations, but not the regime any current figure
targets.

**5. Provenance: one benign drift event, verified by hand.** The detector
fired at N=1024 run 3 (`git_head` changed mid-job). Traced by hand: commit
`256931c` ("Clarify legacy register-equality similarity naming") touched
only `README.md` and `paper/outline.md` -- neither `python/hammock/`, `cpp/`,
nor `bindings/` -- so it did not affect either tool's measured behavior. Full
run otherwise clean (no other warnings, no errors).

**Still open**: job-level n=1 on `--exclusive` (one allocation, one node);
the CPU-model confound on the N=32/128 *magnitude* shift specifically; and
the N=1024→2048 reversal mechanism above is a plausible, code-consistent
candidate, not a checked one -- nobody has yet timed the sketch phase and
pairwise phase separately at N=2048 on both tools to confirm the pairwise
share is actually large enough to explain an 11-13-percentage-point swing.

## Next steps

1. ~~**Fix Part 1**~~ — **done, v0.7.1 (2026-08-11).**
2. ~~**Re-run Part 2 on `--exclusive`, extend to N=2048, `--no-metrics`
   arm**~~ — **done, job 29763124 (2026-08-11), see update above.**
3. If the N≈8-1024 crossover is judged worth chasing on its own: prototype
   parallelizing hammock-cpp's per-file sketch loop and re-measure against
   bedtools. Would need to also explain/survive the N=2048 reversal (item 4
   below) to be worth doing -- porting the CLI's dispatch strategy into
   hammock-cpp only helps if the sketch phase is still the bottleneck at the
   N where it would be applied.
4. **Do NOT add the Python CLI as a competing timed arm in
   `benchmark_cpp_vs_bedtools.py`** (the idea behind this update). Considered
   and rejected on the strength of the N=2048 data above: the headline
   figures' regime (N up to 2048) is exactly where hammock-cpp remains
   faster, so a CLI arm would not beat hammock-cpp where it matters, only in
   a band (N≈8-1024) none of the current figures target. Revisit only if a
   future figure specifically needs that N band.
5. **Not yet started.** Confirm or refute the pairwise-phase-dilution
   mechanism (item 3 above) by timing the sketch phase and pairwise phase
   separately at N=1024 and N=2048 on both tools (`--verbose` on
   `hammock-cpp`; the Python CLI would need an equivalent breakdown added, or
   time `_sketch_many` vs the pairwise call directly). This is the check the
   mechanism claim is missing -- don't treat item 4's "hammock-cpp is faster
   at N=2048" as mechanistically understood, only as measured.
6. **Not yet started.** Replicate Arm B (or a smaller N=32/128-only rerun) on
   a second `--exclusive` allocation to separate the CPU-model confound from
   exclusivity itself, if the N=32/128 *magnitude* (as opposed to direction)
   ever becomes load-bearing for a claim.

## Reproducing the measurements above

Part 1's isolation test (single file per side, no file-dispatch parallelism
possible) — the actual script used, defaults matching the table above:
```bash
python3 experiments/bedtools_benchmark/measure_threads_isolation.py \
    --threads-list 2,4,16 --reps 3 --num-intervals 200000 --precision 18
```
Writes `experiments/bedtools_benchmark/results/threads_isolation_<stamp>.csv`
and prints the same median-CPU% table shown above. Binaries default to
`$HAMMOCK_CPP_BIN` (falls back to the local `build/` binary) and
`/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock`; override with
`--hammock-cpp-bin`/`--hammock-cli-bin` if needed.

Part 2's N-sweep — the actual invocation behind job 29758101:
```bash
python3 experiments/bedtools_benchmark/measure_cli_overhead.py \
    --num-files 2,8,32,128 --num-intervals 10000 \
    --precision 18 --threads 16 --runs 5 --corpus-seed 20260811
```
Part 2's re-measurement (done, job 29763124 -- see "Part 2 update" above) --
`sbatch_cli_overhead.sh` as checked in now runs both arms; to reproduce
directly:
```bash
sbatch experiments/bedtools_benchmark/sbatch_cli_overhead.sh
```
or invoke `measure_cli_overhead.py` directly per-arm, e.g. Arm A:
```bash
python3 experiments/bedtools_benchmark/measure_cli_overhead.py \
    --num-files 2,8,32,128,512,1024,2048 --num-intervals 10000 \
    --precision 18 --threads 16 --runs 6 --corpus-seed 20260811 \
    --cpp-metrics-arm no-metrics
```

## Related

- A limitations-section sentence about the CLI-vs-hammock-cpp gap can now be
  written with real numbers instead of being blocked: hammock-cpp is the
  right speed-ceiling choice at the N the headline figures actually use
  (N=2048: hammock-cpp ~13% faster), but the Python CLI is genuinely faster
  in a middle band (N≈8-1024, up to ~28% faster around N=32-128) that no
  current figure targets. A flat "CLI is ~1.5x slower" or "CLI is faster"
  would both be wrong depending on N -- the "Part 2 update" table above is
  the citable source.
- `docs/seed-benchmark-methodology.md` — the shared-partition noise floor;
  Part 2's N=32/128 *magnitude* (as opposed to direction) still hasn't been
  checked against it independently of the CPU-model confound (see "Still
  open" above).
