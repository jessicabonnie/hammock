# The bedtools parallelism caveat

## Summary

`bedtools jaccard` is single-threaded by design. Whenever this codebase
shows "bedtools at t=N" we mean **N concurrent independent
`bedtools jaccard` processes orchestrated by GNU parallel**, not one
multi-threaded bedtools call. That's a workflow we built around the tool
(`experiments/bedtools_benchmark/bedtools.sh`), not behavior of bedtools
itself. The same wrapper was used in the original hammock October 2025
benchmark, so existing comparisons are internally consistent — but it
matters how you frame the result, because "hammock parallelizes itself"
and "bedtools parallelizes through a wrapper" aren't the same kind of
parallelism.

## What "threads" means for each tool

| Tool       | Threads = N means…                                        |
| ---------- | --------------------------------------------------------- |
| hammock-cpp | One process, N OpenMP threads sharing memory; sketches built and compared in parallel within the process |
| bedtools   | N independent single-threaded `bedtools jaccard` processes spawned by `parallel --jobs N`, one per (a, b) pair, started/ended N at a time |

Both use ~N cores. They differ in:

- **Process model.** Hammock has one address space; bedtools-via-parallel
  has up to N address spaces alive at once.
- **Per-pair overhead.** bedtools-via-parallel pays process spawn,
  binary load, and file open for every single one of the N² pair calls.
  Hammock loads each input file once, sketches it, then does N² compares
  on cached sketches in memory.
- **Memory accounting.** `/usr/bin/time -v` on the bedtools wrapper
  reports the *largest single child's* peak RSS (~19 MB), not the sum
  across the N concurrent processes. Hammock's reported RSS is the
  whole-process peak. So hammock RSS is higher-fidelity; the bedtools
  number under-counts true memory pressure under parallel fan-out.
- **Cache friendliness.** Hammock's threads share cache lines for
  sketch state. bedtools' processes don't share anything.

## Two ways to read a "hammock vs bedtools" comparison

### Tool-vs-tool (strict): not apples-to-apples

If the question is *"which tool is faster on its own merits?"*, then a
hammock-`--threads 16` vs bedtools-`-`(no flag) comparison is the honest
one. By that yardstick bedtools is single-threaded and hammock has
built-in parallelism, period — the numbers we're plotting are not the
right ones.

### Workflow-vs-workflow: yes, apples-to-apples

If the question is *"what's the fastest way to compute pairwise jaccard
on N files?"* — which is what an HPC user actually cares about — then
GNU parallel is the standard idiom for bedtools, and the comparison
*does* hold. The orig hammock benchmark used the same idiom, so all
hammock-vs-bedtools comparisons in this codebase and its predecessor
are workflow-vs-workflow comparisons.

The one we report is the workflow comparison. We surface the strict
tool-vs-tool baseline alongside it (see below) so the framing is
visible rather than hidden.

## Empirical: the parallelism premium is NOT a constant (updated 2026-08-09)

**Read this before the table below.** The premium has now been measured twice on
the same workload and disagrees by 3×. It is a property of the node's
process-creation throughput on the day, not a property of either tool, and the
larger of the two numbers is the one that flatters hammock's baseline least —
so do not pick one and quote it.

| N=64, 10k intervals/file | t=1 | t=2 | t=4 | t=8 | t=16 | t=32 | t=48 |
|---|---|---|---|---|---|---|---|
| 2026-05 wall (s) | 57.0 | 37.8 | 19.4 | 10.3 | 10.1 | — | — |
| 2026-05 speedup | 1.00 | 1.51 | 2.94 | 5.54 | **5.65** | — | — |
| 2026-08-09 wall (s) | 61.62 | 37.53 | 21.82 | 18.66 | 32.06 | 32.45 | 32.53 |
| 2026-08-09 speedup | 1.00 | 1.64 | 2.82 | **3.30** | 1.92 | 1.90 | 1.89 |
| 2026-08-09 cores busy | 1.4 | 2.0 | 3.9 | 5.5 | 3.6 | 3.6 | 3.6 |

The 2026-08-09 row is job 29652415 (node c707, one exclusive allocation, p=18,
3 replicates, medians) using the **fixed** one-process-per-pair `bedtools.sh`,
which is *faster* per pair than the version behind the 2026-05 row. It still
scales worse. Two findings follow.

**1. bedtools does not merely plateau — it regresses.** Wall time rises from
18.66 s at t=8 to 32.53 s at t=48 and stays there, so 16, 32 and 48 concurrent
jobs are all *slower* than 8. `cpu_time/wall_time` shows why: the workflow keeps
at most ~5.5 cores busy and falls back to 3.6 once oversubscribed. The 2026-05
table stopped at t=16 and so could not see this. For contrast, hammock on the
same node and workload is monotonic to 18.81× and keeps 43.3 of 48 cores busy.

**2. The ceiling is process creation, and it is not bedtools' fault.** A
pairwise workflow launches one process per pair — N² of them — and on these
nodes process creation caps near **123 exec/s** and does not scale with cores.
Measured controls, all on the same node:

- `md5sum` on node-local files: **0.46×** at 16-way — *slower* than serial.
  Nothing to do with bedtools.
- `xargs -P16` hits the same ceiling as GNU `parallel` (8.26 s vs 7.74 s for
  1024 pairs), so it is not the dispatcher.
- Copying the bedtools binary from GPFS to local NVMe changes nothing
  (1.46× vs 1.48×), so it is not filesystem latency on the executable.
- GNU parallel itself is healthy: 32 × `sleep 1` at `-j16` takes 2.19 s
  (14.6× parallelism), and it dispatches shell builtins at 468 jobs/s.

**Consequence for every speedup this repo reports.** "bedtools at t=16" can
silently mean "bedtools at t≈1.5", which inflates a quoted speedup by up to ~6×.
Job 29651772 was cancelled 40 minutes in for exactly this. Two mitigations are
now in place:

- `bedtools.sh` runs **one** process per pair (`parallel --tagstring` + a single
  `awk`) instead of three (`bedtools | tail | cut` inside an exported function).
  Worth ~2.1× consistently: measured 0.81→1.68 (exclusive c599), 1.27→2.86
  (c516), 0.62→1.17 (shared/sr08) in achieved efficiency.
- Every bedtools row in `benchmark_cpp_vs_bedtools.py`'s CSV now carries
  `mean_bedtools_serial_ms` and `mean_bedtools_parallel_eff`, so the achieved
  parallelism is recorded rather than assumed. **Read it before quoting a
  speedup.** It is only meaningful at large N; below about N=32 the bedtools leg
  is a few pairs behind a fixed ~0.5 s of startup and reads ~0.01.

Figure: `paper/figures/threading_supplement.png`
(`paper/pairwise_scaling/plot_threading_supplement.R`, data
`docs/data/sweep_threads_p18.csv`).

### The original (2026-05) measurement, kept for the record

From `experiments/bedtools_benchmark/sweep.py --axis threads` at
N=64 files × 10k intervals/file:

| threads | bedtools wall (s) | bedtools speedup vs t=1 |
| ------- | ----------------- | ----------------------- |
| 1       | 57.0              | 1.00× (raw bedtools)    |
| 2       | 37.8              | 1.51×                   |
| 4       | 19.4              | 2.94×                   |
| 8       | 10.3              | 5.54×                   |
| 16      | 10.1              | 5.65×                   |

Going from raw single-threaded bedtools to GNU-parallel-wrapped at
t=16 buys ~5.7× wall-time speedup on this hardware. That's the
"parallel idiom premium" — it's a real number, but it's a property
of the workflow, not of the bedtools binary. Note also the diminishing
return between t=8 and t=16 (10.3 → 10.1): bedtools wall stops scaling
because per-pair cost is small enough that GNU parallel's spawn
overhead and node I/O bandwidth become the bottleneck, not core count.

> **Superseded as a quotable number (2026-08-09).** "~5.7×" is one observation,
> not the premium. The same measurement on the same workload now reads 1.92× at
> t=16, and the diagnosis in the paragraph above — spawn overhead becoming the
> bottleneck — turned out to be right in kind and badly understated in degree:
> the ceiling is ~123 exec/s and it makes wall time *rise* past t=8. See the
> updated section at the top of this file. The hammock rows below are from the
> same 2026-05 run and are p=14, so they are not comparable to the p=18 numbers
> quoted elsewhere either.

Hammock-cpp at the same N=64, p=14:

| threads | hammock wall (s) | speedup vs t=1 |
| ------- | ---------------- | -------------- |
| 1       | 637.1            | 1.00×          |
| 2       | 321.0            | 1.98×          |
| 4       | 162.2            | 3.93×          |
| 8       | 82.95            | 7.68×          |
| 16      | 42.7             | 14.92×         |

Near-linear scaling all the way out — OpenMP thread parallelism on
shared memory pays off cleanly here.

> **Superseded 2026-08-06 — the hammock table above is from a
> pre-optimization build and has no surviving CSV.** No archived
> `sweep_threads_*.csv` contains a 637 s t=1 point; the earliest one on disk
> (`sweep_threads_20260510_172628.csv`) already reads 198.9 s. The bedtools
> table above matches the May 11/12 sweeps to within run-to-run noise, so the
> two tables in this section are from *different builds* and their ratios should
> not be compared to each other.
>
> Canonical current numbers, `sweep_threads_20260512_150458.csv` (the run
> `experiments/bedtools_benchmark/RESULTS.md` tabulates), N=64, p=14, subB=1.0:
>
> | threads | bedtools | hammock-cpp |
> | --- | --- | --- |
> | 1 | 58.11 s | 72.52 s |
> | 2 | 38.20 s | 42.41 s |
> | 4 | 19.63 s | 24.01 s |
> | 8 | 10.26 s | 12.45 s |
> | 16 | 9.90 s | 6.63 s |
>
> Hammock's OpenMP scaling is **10.9×** at t=16, not 14.9× — still clean, and
> still the point being made, but the optimized binary starts from a much lower
> t=1 so there is less headroom. The qualitative contrast survives unchanged:
> hammock keeps scaling from t=8 to t=16 (12.45 → 6.63 s) while
> GNU-parallel-wrapped bedtools flattens (10.26 → 9.90 s).

## How to read summary plots

Wall-time-vs-N plots in this codebase should show three reference
curves:

1. **Raw bedtools** (single-threaded, t=1) — the tool by itself
2. **Parallel-wrapped bedtools** (GNU parallel, t=16) — the workflow
3. **hammock-cpp Mode B** (OpenMP, t=16) — hammock's built-in parallelism

This makes the parallel idiom an explicit step in the comparison rather
than a hidden assumption. If a plot only shows curves 2 and 3, the
delta between hammock and bedtools-with-parallel-wrapping is what's
being reported — *not* the delta between hammock and "what bedtools
does on its own."

## Where this matters in experiments

- **`experiments/bedtools_benchmark/sweep.py --axis threads`** — every
  data point uses parallel-wrapping at the listed thread count. The
  scaling curve for bedtools is GNU parallel's scaling, not bedtools'.

- **`experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py`
  (num_files axis)** — same caveat. The crossover N where hammock
  starts winning depends heavily on which "bedtools" you compare to.
  The orig October 2025 run showed a crossover at N≈128 against
  parallel-wrapped bedtools at t=16; against raw bedtools the
  crossover would land at much smaller N.

- **t=16 N-up-to-512 follow-up** (job 23901950) — re-runs the orig
  protocol on our hardware. **Completed**; results are in "Empirical
  follow-up: t=16 replication of the orig benchmark" below. That protocol
  has since been re-run again on Aug 4 2026 (job 29552415,
  `docs/data/cpp_vs_bedtools_t16_20260804_172242.csv`); see
  `docs/figure3-panel-a-rebuild.md`. The crossover conclusion is unchanged
  across all three runs.

## Where this *doesn't* matter

- The HLL precision/accuracy story (`docs/jaccard-definitional-gap.md`)
  is tool-internal — bedtools' parallelism doesn't enter.
- Mode D parity tests against the conda `hammock` install
  (`tests/test_mode_d_parity.py`) — both sides are hammock; no
  bedtools.

## Empirical follow-up: t=16 replication of the orig benchmark

Run on 2026-05-10 with the optimized hammock (post-May-10 build) at
t=16, p=14, 10k intervals/file, N up to 512, 3 runs/config. Source:
`experiments/bedtools_benchmark/results/cpp_vs_bedtools_t16_20260510_181334.{txt,csv}`.

| N   | orig bedtools wall | new bedtools wall | orig hammock wall | new hammock wall | orig speedup | new speedup |
| --- | ------------------ | ----------------- | ----------------- | ---------------- | ------------ | ----------- |
| 64  |    33.2 s |   10.4 s |    46.2 s |   13.7 s | 0.72× | 0.76× |
| 128 |   134.9 s |   39.6 s |    92.9 s |   27.3 s | 1.45× | 1.45× |
| 256 |   551.3 s |  157.8 s |   186.2 s |   54.6 s | 2.96× | 2.89× |
| 512 |  2196.4 s |  645.3 s |   376.8 s |  109.4 s | 5.83× | 5.90× |

**The crossover behavior is preserved.** Both hammock and bedtools are
2–3× faster in absolute terms on the May 2026 hardware/build, but the
*ratio* — the thing the comparison plot reports — is essentially
identical. Crossover at N=128, ~6× speedup at N=512. So the orig's
"hammock wins above N≈128" claim replicates with the optimized
hammock and the current bedtools baseline.

## Sort time and the fairness question

Pre-sort time from the 2026-05-10 evening rerun, with the sort step
parallelized across 16 workers
(`ThreadPoolExecutor` of `subprocess` calls to `sort`):

| N   | bedtools wall | sort wall (parallel, 16w) | sort as % of bt workflow | bt + sort | hammock wall | speedup (with sort) |
| --- | ------------- | ------------------------- | ------------------------ | --------- | ------------ | ------------------- |
|  64 |  11.5 s |  0.74 s |  6% |  12.3 s |  14.1 s | 0.87× |
| 128 |  44.9 s |  1.43 s |  3% |  46.3 s |  27.5 s | 1.68× |
| 256 | 184.5 s |  2.79 s |  1% | 187.3 s |  55.0 s | 3.40× |
| 512 | 767.5 s |  5.52 s |  1% | 773.0 s | 110.1 s | 7.02× |

With parallel sort the pre-sort step is **a few percent of the
bedtools workflow at most** — the "fairness gap" between the
pre-sorted and unsorted-input framings has essentially vanished.
Crossover sits at **N=128 in either framing**, with hammock's lead
growing to ~7× at N=512.

For context, here is the same table from the earlier sequential-sort
run (2026-05-10 morning), where sort was a much bigger story:

| N   | bedtools wall | sort wall (sequential) | sort % | bt + sort | hammock wall | speedup |
| --- | ------------- | ---------------------- | ------ | --------- | ------------ | ------- |
|  64 |  10.4 s |  10.4 s | 50% |  20.8 s |  13.7 s | 1.52× |
| 128 |  39.6 s |  20.6 s | 34% |  60.2 s |  27.3 s | 2.20× |
| 256 | 157.8 s |  41.2 s | 21% | 199.0 s |  54.6 s | 3.64× |
| 512 | 645.3 s |  82.5 s | 11% | 727.8 s | 109.4 s | 6.65× |

Sort scales linearly in total intervals (O(N·intervals/file)) and
parallelizes embarrassingly across files, so going from sequential
to 16-way parallel cuts it ~15× (close to the theoretical 16×
ceiling, limited only by per-file `subprocess` setup overhead).

**Takeaway:** the "include sort in the bedtools workflow" comparison
matters only when (a) the user's sort step is sequential, or (b) the
sort dataset is small enough that fan-out overhead dominates. In any
realistic HPC setup, the pre-sorted and unsorted-input framings give
the same answer for which tool wins.

Both views are honest. We capture sort time in every CSV and
text report (`mean_sort_time` / `sort_time`) so downstream analysis
can pick whichever framing fits the use case.

## References

- Wrapper script: `experiments/bedtools_benchmark/bedtools.sh`
- Threads-sweep results: `experiments/bedtools_benchmark/results/sweep_threads_*.csv`
- Orig benchmark for comparison: `/home/jbonnie1/interval_sketch/hammock/benchmarks/`
- Sibling doc: `docs/jaccard-definitional-gap.md`
