# CLAUDE.md

Guidance for Claude Code working in this repository. Internal notes — keep
divergences and unimplemented-feature info here, not in README.md.

## What this is

A Python + C++ refactor of [`hammock`](https://github.com/jessicabonnie/hammock).
Same CLI; HyperLogLog backing for all four modes (A/B/C BED-interval, D FASTA);
parity-tested byte-equal CSV output against the orig on tested paths.

## Build / test

```bash
pip install -e . --no-build-isolation   # builds the C++ extension
pytest tests/
```

The wheel includes a standalone `hammock-cpp` binary built from the same
`hammock_core` static lib (in `build/`); intended for max-speed Mode B
benchmarking, no Python in the loop.

## Architecture

- `python/hammock/` — CLI (`cli.py`), orchestrator (`runner.py`), Mode D
  (`modes/sequence.py`). Threading via `concurrent.futures.ThreadPoolExecutor`;
  default `--threads = min(8, cpu_count())`. The C++ extension and `digest`
  both release the GIL, so threads are real parallelism.
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
2. **Containment column is well-defined.** Orig's was a placeholder. Ours:
   `(intersection / cardinality(self)) ** expA`, with `1.0` sentinel when
   `expA == 0`. Parity tests project this column out before comparing.
3. **`--mixed-stride`** is a real flag. Orig's pipx-installed 0.4.0 doesn't
   accept it (it lived only in WIP changes). Historical mixed-stride results
   from the orig need to be re-run; see `memory/project_mixed_stride_rerun.md`.

## Not implemented (would need work to add)

- Sketch types: `--minhash`, `--exact`, `--minimizer` for A/B/C. Phase 1
  shipped HLL only. The `BagMinHashSketch` C++ class exists; the CLI/runner
  glue doesn't.
- Mode D BED→FASTA flags `--bed1/--ref1/--bed2/--ref2`. Phase 5 in the plan.
- Reference-name resolution (`--reference mm10`). Phase 6 (aspirational).
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
