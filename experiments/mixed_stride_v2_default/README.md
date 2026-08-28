# Mixed-stride v2 default verification

This experiment verifies the promotion of the chromosome-anchored
mixed-stride v2 sampler to the public `mixed-stride` default. It is isolated
from the archived `subB_mixed_stride` experiment so its raw outputs cannot be
mistaken for, or overwrite, the historical v1 results.

The runner uses the 20 Maurano fetal-tissue DHS BED files, Mode B, precision
18, eight threads, and `jaccard_similarity_ie`. It runs an unsubsampled
`subB=1.0` accuracy baseline and compares public/default `mixed-stride` with
explicit legacy `mixed-stride-v1` at:

```text
0.3, 0.1, 0.03, 0.01, 0.008, 0.005, 0.003, 0.001
```

Rates `0.1`, `0.01`, `0.008`, `0.005`, and `0.001` have integral reciprocal
strides, so the two methods must emit exactly equal pairwise values. Any
mismatch at those rates is a hard failure. Rates `0.3`, `0.03`, and `0.003`
have non-integral reciprocals and are expected to demonstrate that v2 differs
from v1. The summary reports wall time, mean/max absolute per-pair change from
the unsubsampled baseline, and direct v1-v2 differences for all 190 unique
non-self pairs.

Run from the repository root after building `hammock-cpp`:

```bash
python experiments/mixed_stride_v2_default/run_maurano_parity.py
```

Each invocation creates a timestamped directory under `results/` containing
the file list, raw hammock TSVs, captured stderr, and `summary.csv`.

## Completed run

`results/run_20260828_141330/summary.csv` records the branch verification run.
All 190 pairwise values were exactly equal between v1 and v2 at every
integral-reciprocal rate. All 190 differed at each non-integral rate, as
expected. At the exact-parity rates, the single-run v2/default wall times
differed from v1 by at most 1.1%; these timings are a regression smoke test,
not a replicated performance estimate.
