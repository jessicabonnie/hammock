# Seed: CLI's `--threads` doesn't bound sketching parallelism, and a related, less-certain hammock-cpp-vs-CLI wall-time crossover (found 2026-08-11)

This seed carries two claims at very different confidence levels. Keep them
separate — the first is a confirmed bug with a planned fix; the second is a
direction-only observation that does not yet meet this repo's own bar for
"confirmed" (see `CLAUDE.md`'s retraction section: "verifying mechanisms
without verifying magnitudes was the recurring failure").

## Part 1 — CONFIRMED BUG, TO BE FIXED: the Python CLI's `--threads` does not bound sketching-phase parallelism

**To fix.** This is bad independent of any wall-time benefit it happens to
produce: a user who passes `--threads 4` — e.g. because that's their actual
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
magnitude-checking standard Part 2 is held to:**

1. **The mechanism is proven by code inspection alone — dispositive, not
   statistical.** The pybind11 binding for sketching, `sketch_bed_file_hll`
   (`bindings/_core.cpp:354-367`), takes no thread argument at all — it only
   wraps the call in `py::call_guard<py::gil_scoped_release>()`.
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

**Fix direction (not started).** Thread the CLI's `threads` budget into the
sketching phase the same way `omp_threads` already reaches the pairwise
phase: either give `sketch_bed_file_hll` a thread-count parameter that
becomes a `num_threads()` clause inside `process_bed_file_mode_b`'s pragma, or
call `omp_set_num_threads` once from Python before sketching starts (matching
what `hammock-cpp` already does). Needs a decision on how this interacts with
the outer `ThreadPoolExecutor` — see Part 2's mechanism, since naively
clamping the inner team to `threads` while the outer pool is *also* sized to
`threads` would change the wall-time behavior investigated there, not just
fix the resource-accounting bug.

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

## Next steps (not started)

1. **Fix Part 1** (the confirmed `--threads` bug) — see "Fix direction" above.
   Note this will *change* Part 2's mechanism once shipped (the sketching
   phase's oversubscription goes away), so re-measure the crossover after
   fixing, not just before.
2. Re-run Part 2's N=2/8/32/128 comparison on `--exclusive`, and extend to
   N=512/1024/2048 (where the paper's real claims live) and to the
   `--no-metrics` arm (what the figures actually use) — see the
   falsification criteria above for what a negative result would look like.
3. If the crossover survives all of the above: prototype parallelizing
   hammock-cpp's per-file sketch loop and re-measure against bedtools. Real
   potential improvement to the paper's headline numbers if so — but new
   scope, not something to fold into the current session.

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
For next step 2 (re-run on `--exclusive`, extend N, switch to `--no-metrics`
for the cpp side to match what the headline figures actually use): add
`--exclusive` to `sbatch_cli_overhead.sh`'s header and widen `--num-files` to
`2,8,32,128,512,1024,2048`. Note `measure_cli_overhead.py` currently hardcodes
`--metrics` on the cpp side (see "Important scope note" above) — switching
that to `--no-metrics` needs a small script edit, not just a flag change.

## Related

- A limitations-section sentence about the CLI-vs-hammock-cpp gap was blocked
  on this seed rather than written, since a simple "CLI is ~1.5x slower"
  claim would be actively misleading given the sign flip in Part 2.
- `docs/seed-benchmark-methodology.md` — the shared-partition noise floor
  Part 2's magnitude hasn't been checked against yet, and the model for what
  "independently reproduced" should actually require in this repo (a
  committed artifact, not prose).
