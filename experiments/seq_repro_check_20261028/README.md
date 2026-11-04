# Sequence-mode reproducibility check (2026-10-28)

Asks whether current `main` reproduces bit-for-bit the
archived Mode D `--sequence-hll-hash rehash-selector64` outputs behind Figs 5-9.

- Source: `main` at 05cb64529178b64b3e16f13ebc09623d356ab43e, cloned to
  `/vast/blangme2/jbonnie/hammock-seq-repro-20261028/hammock`.
- Build: `pip install --no-deps --no-build-isolation --target
  /vast/blangme2/jbonnie/hammock-seq-repro-20261028/runtime/site-packages` with the
  `claude-ref-comparison` conda env python and its conda gcc (now 14.3.0; the archived
  `_core` was built with conda gcc 13.3.0). Log: `.../build.log`. The archived
  `sitecustomize.py` was copied in to drop the env's editable-hammock finder.
- Runtime: `claude-ref-comparison` env python 3.10.20, `LD_PRELOAD` env libstdc++,
  `PYTHONPATH` = isolated site-packages (mirrors `maurano_dhs_rehash_validation/workflow/hammock_runtime.sh`).
- digest: single unversioned extension
  `claude-ref-comparison/lib/python3.10/site-packages/digest.cpython-310-x86_64-linux-gnu.so`
  (mtime 2025-05-28, sha256 db922f2b16dde072220f5829d6c10060998b791285d0745995907af334596426,
  built with Spack GCC 9.3.0). No dist-info/version metadata exists. The archived runs
  (Aug 2026) used the same env and the same path (the archived runtime site-packages contain
  only hammock), and the file predates them, so old and new use the same digest binary.
  xxhash 3.5.0, numpy 2.2.6 in both.
- Archives (read-only): Maurano CSVs under
  `/vast/blangme2/jbonnie/hammock-mode-d-rehash-sweep.sCydsC/hammock/experiments/maurano_dhs_rehash_validation/results/rehash/`
  (commit 008e9f24); Fig 5 CSV
  `/vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/exp_a_rehash_selector64/seed42/broad/k20_w20/`
  (p=24, the precision of `docs/data/exp_a_broad_k20_w20_full.csv`, which is identical to it).
- New outputs: `/vast/blangme2/jbonnie/hammock-seq-repro-20261028/results/task*/`.

Files: `jobs.tsv` (30 tasks; task 29 is a control re-running the *archived* runtime
today), `run_task.sh` (SLURM array, `sbatch --array=0-29 run_task.sh`), `compare.py`
(writes `comparison.tsv`; compares file1/file2 keys, row order, metadata columns and the
seven estimate columns as exact text, plus max abs difference).

SLURM job: 32635054 (array 0-29), all COMPLETED. Result: 30/30 IDENTICAL; see comparison.tsv. Whole files are also byte-identical after CRLF->LF normalization.

## Addendum (2026-10-28): cross-reference at p=18

Figure 5 now uses the p=18 cross-reference run (`exp_a_rehash_selector64_p18`),
so tasks 30 (broad) and 31 (narrow), k=20,w=20,p=18,seed 42, were added
(SLURM 32636312, 32636317). Both are IDENTICAL (81 rows, max abs diff 0).
Runs relevant to the current manuscript: tasks 0-27, 30, 31 (30 runs); task 28
(p=24) matches an input no longer used by any figure; task 29 is the control.
