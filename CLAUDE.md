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

The bed2fasta tests (`tests/test_bed2fasta*.py`) and their `--ref` end-to-end
paths need `bedtools` (and `samtools` for indexing) on `PATH` — `ml bedtools
samtools` on the cluster; they self-skip otherwise. Mode D parity needs the
conda-orig env (see Parity environments).

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
2. **Containment + co-sketch columns.** Orig's single `containment` column
   was a placeholder. We replace it with a five-column block, computed
   from the HLL register-equality intersection:
   - `containment_AB = |A ∩ B| / |A|`
   - `containment_BA = |A ∩ B| / |B|`
   - `cosketch_geom  = sqrt(C_AB · C_BA)`
   - `cosketch_arith = (C_AB + C_BA) / 2`
   - `cosketch_max   = max(C_AB, C_BA)`

   No `expA` exponent is applied — the orig's `** expA` semantics on a
   placeholder column were never load-bearing, and `expA` already changes
   what `A` *contains* (interval multiplicity) upstream of the intersection.

   Mode D emits the same block twice: once on the minimizer HLL (paired
   with `jaccard_similarity`) and once on the merged minimizer ∪ start-end
   HLL (paired with `jaccard_similarity_with_ends`, columns suffixed
   `_with_ends`).

   Parity tests project these columns out before comparing.
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
   At low subB it's slightly slower than hash-threshold; only at subB
   near 1.0 is it marginally faster. Kept as a research/comparison flag.
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
  shipped HLL only. The `BagMinHashSketch` C++ class exists; the CLI/runner
  glue doesn't.
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
