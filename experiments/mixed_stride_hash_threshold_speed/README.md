# Mixed-stride versus hash-threshold speed experiment

This experiment measures the speed gained by hammock's public/default
`mixed-stride` sampler over the original-compatible `hash-threshold` sampler
on the 20-file Maurano fetal-tissue DHS corpus.

It is deliberately isolated from `experiments/subB_mixed_stride/`, whose
checked-in results are historical and should not be silently replaced.

## Predeclared design

- Mode B, HLL precision 18, eight threads.
- HLL seed 42 and sampling/gate seed 31337, passed explicitly.
- Methods: `mixed-stride` and `hash-threshold` only.
- Rates: `0.1, 0.06, 0.03, 0.01, 0.006, 0.003, 0.001`.
- Ten measured rounds after one unmeasured warm-up round.
- Every round contains all 14 method/rate cells in a deterministic shuffled
  order. Comparisons are paired by round to reduce node-level timing drift.
- Run on one exclusive Slurm node with OpenMP threads bound to cores.
- The bare/default IE output shape is used for both methods. It represents the
  public end-to-end path and keeps output work identical between methods.

The reciprocal is integral at 0.1, 0.01, and 0.001. The other four rates
exercise mixed-stride v2's within-chromosome mixed-gap path.

### Endpoints

The primary endpoint is the paired per-round speedup in hammock's reported
sketch-construction phase:

```text
hash-threshold sketch seconds / mixed-stride sketch seconds
```

The secondary endpoint is the analogous end-to-end wall-time speedup. CPU
time, pairwise/write time, maximum RSS, sampled-point count, realized sampling
fraction, output digest, binary digest, commit, node, and CPU model are retained
for diagnosis and provenance.

For each rate, report the median paired speedup and a deterministic 95%
percentile bootstrap confidence interval over rounds. Also report each
method's median raw time. Do not pool the seven rates into one headline ratio.

### Interpretation gates

- A cell is complete only if it has one run per method in every measured round.
- Repeated output digests must be identical within a method/rate cell.
- The parsed total-point count must agree across every run.
- Sampling fractions are recorded to show that a speed advantage is not caused
  by a materially different requested workload.
- The primary claim is computational: mixed-stride avoids visiting and gate-
  hashing rejected positions. Accuracy remains a separate endpoint from the
  existing Maurano validation and must not be inferred from timing alone.

## Run

From the repository root:

```bash
sbatch experiments/mixed_stride_hash_threshold_speed/sbatch.sh
```

For a local plumbing check (not a publishable timing):

```bash
python experiments/mixed_stride_hash_threshold_speed/run_benchmark.py --smoke
```

The Slurm job writes a timestamped run directory under `results/` containing
`metadata.json`, `files.txt`, `raw_runs.csv`, the first deterministic estimate
file for every method/rate cell, and `summary.csv`.

## Archived result

Job 30459651 completed on node c531 (Intel Xeon Gold 6248R) in 14:44. The
repository retains its path and checksum input manifests, machine and binary provenance, all 140
measured timings, and the validated summary under `results/job_30459651/`.
The per-cell estimate tables are omitted because they are regenerable and
their SHA-256 digests are recorded in `raw_runs.csv`.
The corresponding Slurm stdout/stderr are retained under `logs/`.

Across the requested fractions, mixed-stride was 1.736--2.271x faster in the
reported sketch-construction phase and 1.727--2.237x faster end to end. Every
end-to-end comparison at `subB <= 0.01` exceeded 2x. Realized fractions were
closely matched; the largest absolute difference between methods was
4.09e-6 (`subB=0.003`). These measurements establish a computational
advantage at matched sampling density; accuracy is evaluated separately in
the Maurano validation.

The run used a checkout reported as dirty because the benchmark harness and
sampling-validation tests were not yet committed. The executable is pinned by
its SHA-256 digest in `metadata.json`, but the commit alone does not reconstruct
the exact source state. This limitation should remain visible when the result
is cited; a clean-checkout rerun is preferable before final publication.

Generate the supplementary figure from the archived summary with:

```bash
python experiments/mixed_stride_hash_threshold_speed/plot_results.py
```

Suggested supplementary-results text:

> On the 20-file Maurano DHS benchmark (p=18, eight threads), mixed-stride
> subsampling reduced end-to-end runtime by 1.73--2.24-fold relative to
> hash-threshold sampling across subB=0.1--0.001 while producing closely
> matched realized sampling fractions. The advantage exceeded twofold at
> every subB value at or below 0.01.

Suggested caption:

> **Mixed-stride avoids the per-position cost of hash-threshold sampling.**
> (A) Median sketch-construction time for each method. (B) Median paired
> hash-threshold/mixed-stride runtime ratio over
> ten randomized rounds on one exclusive node; values above one favor
> mixed-stride. Error bars are deterministic percentile-bootstrap 95%
> intervals over the ten paired rounds. Sketch construction is the
> predeclared primary endpoint and total wall time is secondary. Both methods
> used the same 20 Maurano DHS files, precision 18, eight threads, HLL seed 42,
> and sampling seed 31337. Intervals quantify within-node timing variability,
> not generalization across nodes or datasets.
