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

## Empirical: the parallelism premium for bedtools

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
