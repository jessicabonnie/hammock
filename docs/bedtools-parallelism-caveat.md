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
  protocol at our hardware. Numbers will be filled in here when it
  completes.

## Where this *doesn't* matter

- The HLL precision/accuracy story (`docs/jaccard-definitional-gap.md`)
  is tool-internal — bedtools' parallelism doesn't enter.
- Mode D parity tests against the conda `hammock` install
  (`tests/test_mode_d_parity.py`) — both sides are hammock; no
  bedtools.

## Empirical follow-up: t=16 replication of the orig benchmark

*(To be filled in after job 23901950 completes.)*

The original Oct 2025 benchmark ran t=16, p=16, 10k intervals/file,
N up to 512 (`hammock/benchmarks/results/cpp_vs_bedtools_t16_20251030_101826.txt`):

| N | orig bedtools wall (s) | orig hammock wall (s) | orig speedup |
| - | --------------------- | --------------------- | ------------ |
| 64  |   33.2 |   46.2 | 0.72× |
| 128 |  134.9 |   92.9 | 1.45× |
| 256 |  551.3 |  186.2 | 2.96× |
| 512 | 2196.4 |  376.8 | 5.83× |

We replicate this on our hardware to disentangle "did hammock change?"
from "did the bedtools baseline shift between Oct 2025 and now?". The
files-axis sweep at t=8 already showed our bedtools running ~3× faster
than orig at comparable N, suggesting the baseline did move. The t=16
result will say whether the crossover behavior is preserved.

## References

- Wrapper script: `experiments/bedtools_benchmark/bedtools.sh`
- Threads-sweep results: `experiments/bedtools_benchmark/results/sweep_threads_*.csv`
- Orig benchmark for comparison: `/home/jbonnie1/interval_sketch/hammock/benchmarks/`
- Sibling doc: `docs/jaccard-definitional-gap.md`
