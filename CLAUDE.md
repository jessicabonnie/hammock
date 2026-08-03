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

**Re-validate before changing floating-point build flags.** The Ertl
estimator feeds every cardinality-derived column (`containment_*`,
`cosketch_*`, `jaccard_similarity_ie`), and its τ/σ series terminate on exact
float equality, so a reassociating flag could in principle shift them.
Reassuringly, a spot check found bit-identical output across `-O0/-O2/-O3`,
every `-march=` tried (x86-64, v2, haswell, skylake-avx512, native), **and**
`-Ofast`/`-ffast-math` — the τ/σ loops are loop-carried, so there is nothing
to reassociate. Treat that as "no known hazard", not a guarantee: if you touch
the flags, re-run the byte-identity gate (capture a Mode B + Mode D CSV before
and after and diff the whole row). `jaccard_similarity` itself is an exact
ratio of two integer register counts and is immune either way.

The bed2fasta tests (`tests/test_bed2fasta*.py`) and their `--ref` end-to-end
paths need `bedtools` (and `samtools` for indexing) on `PATH` — `ml bedtools
samtools` on the cluster; they self-skip otherwise. Mode D parity needs the
conda-orig env (see Parity environments).

The wheel includes a standalone `hammock-cpp` binary built from the same
`hammock_core` static lib (in `build/`); intended for max-speed Mode B
benchmarking, no Python in the loop. The wheel does install it (to
`<site-packages>/bin/hammock-cpp`), but that directory is **not** on `$PATH`,
so invoke it by full path — `build/<wheel-tag>/hammock-cpp` after a local
build, or the site-packages copy.

By default it emits 3 columns (`query`, `reference`, `jaccard_similarity`);
`--metrics` adds `jaccard_similarity_ie`, `containment_AB/BA` and `cosketch_*`,
matching the Python CLI **bit-for-bit** (the IE derivation is written the same
way in both — see `jaccard_ie_from_containments` in `hammock_cli.cpp` and
`runner._jaccard_ie_from_containments`; keep them in sync or
`tests/test_hammock_cpp_metrics.py` fails on `==`). Off by default because it
costs a union + cardinality per pair and would invalidate the timings in
`experiments/bedtools_benchmark/RESULTS.md`. Needed before any interval-mode
rerun: without it the binary cannot emit output from which set-Jaccard is
recoverable.

## Architecture

- `python/hammock/` — CLI (`cli.py`), orchestrator (`runner.py`), Mode D
  (`modes/sequence.py`). Threading via `concurrent.futures.ThreadPoolExecutor`;
  default `--threads = min(8, cpu_count())`. The C++ extension releases the
  GIL, so interval-mode threading is real parallelism.

  **Mode D threading is not — it is a GIL convoy, and the default makes it
  slower.** Measured 2026-08-03 on a 48-core node, 4 Maurano DHS FASTAs,
  k=8 w=40 p=20: `--threads 1` = 31.9 s vs `--threads 8` = **68.3 s (2.1×
  slower)**. In-process, per-task time degrades monotonically (T1 1.00×,
  T2 0.61×, T4 0.86×, T8 0.40×) even with private per-thread sequence lists,
  so it is not refcount contention on shared data. `digest.window_minimizer`
  does *not* usefully release the GIL at this call pattern (~11 µs per call on
  a ~466 bp record), and the surrounding per-record Python work holds it
  outright. Eight concurrent *processes* on the same workload scale ~7.1×, so
  the machine is fine — it is the threading model. Corroborated by the archived
  sweep index (`maurano_dhs_validation/results/sweep_d_*.csv`), where every
  8-thread Mode D cell records `cpu_s/wall_s ≈ 0.8`.

  Consequence: **pass `--threads 1` for Mode D**, and treat every archived Mode D
  wall time as inflated. Switching Mode D to processes (or dropping the pool) is
  worth more than any sketching optimization so far — see divergence #8's note
  that removing `_with_ends` bought 1.5–2.5× single-threaded.
  Full evidence, option space, and reproduction: `docs/seed-mode-d-threading.md`.

## Open seeds (handoff notes for future work)

Neither is a decision; each is evidence gathered plus what still needs
establishing. Read the seed before re-litigating the question.

- `docs/seed-mode-d-threading.md` — Mode D's thread pool is a GIL convoy and the
  default makes it 2–4.5× slower. Blocked route (process pool) needs pickling
  support on `HLLSketch`, which the bindings do not expose.
- `docs/seed-mode-d-hash-width.md` — `digest` returns ≤32-bit minimizer hashes
  while the HLL assumes 64, biasing Mode D cardinality −0.5% to −8.3%. The
  `hash_size=32` option became viable once `_with_ends` was removed (divergence
  #8), since there is no longer a second sketch to merge with.
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
   from the **inclusion-exclusion** intersection — `|A| + |B| - |A ∪ B|`,
   Ertl estimator on each, union by register-wise max, clamped to `>= 0`
   (`HLLSketch::intersection_size`, `cpp/src/hll_sketch.cpp:169`). Same
   formula as orig's `hyperloglog.py estimate_intersection`. This is **not**
   the register-equality path that `jaccard_similarity` uses — see the
   estimator note below, it matters:
   - `containment_AB = |A ∩ B| / |A|`
   - `containment_BA = |A ∩ B| / |B|`
   - `cosketch_geom  = sqrt(C_AB · C_BA)`
   - `cosketch_arith = (C_AB + C_BA) / 2`
   - `cosketch_max   = max(C_AB, C_BA)`

   No `expA` exponent is applied — the orig's `** expA` semantics on a
   placeholder column were never load-bearing, and `expA` already changes
   what `A` *contains* (interval multiplicity) upstream of the intersection.

   Mode D emits this block once, on the minimizer HLL. (It used to emit a
   second `_with_ends` copy on a merged minimizer ∪ start-end HLL; that whole
   family was removed — see divergence #8.)

   Parity tests project these columns out before comparing.

   **The two estimator families in a row are on different scales — do not
   mix them.** `jaccard_similarity` is register-equality, which carries a
   chance-agreement floor `c` (two registers tie when different elements
   have equal ρ). The containment/cosketch block is inclusion-exclusion and
   carries no such floor. Measured on *disjoint* inputs at p=16, n=2×10⁵:
   `jaccard_similarity = 0.168` while `containment_AB = 0.000`. So
   `jaccard_similarity` is *approximately* an affine transform `c + (1−c)·J`
   of set-Jaccard, while the containments estimate the true set quantities.
   Two qualifications, both measured, both load-bearing:

   - **The affine model is only approximate.** Against a `c + (1−c)·J` fit
     with c=0.180 the residual peaks at **+0.025 near J≈0.5**, is −0.010 at
     J=0, and is identically 0 at J=1 (any slope-(1−c) line passes through
     (1,1) by construction), so the error is largest in the *middle* of the
     range, not at the ends. That curvature is why a fitted intercept (0.180)
     exceeds the true disjoint floor (0.1699) — consistent, not a bug.
   - **`c` is set by the load factor λ = n/m *and* by the cardinality ratio
     |A|/|B|**, not by precision as such: 0.1699 as λ→∞ at equal cardinality,
     but 0.152 at ratio 2, 0.097 at ratio 5, 0.058 at ratio 10; and 0.045 at
     p=24 once m > n. It is **not** a constant you can subtract, and
     **`jaccard_similarity` is therefore not rank-faithful across pairs of
     differing cardinality ratio.** On the 20-sample Maurano corpus (ratio up
     to 2.2) it inverts bedtools' ordering for 2.5% of pairs (Kendall
     τ = 0.951). Rank only within comparable pairs.

   Consequence: **`jaccard_similarity_ie` now ships as its own column** (see
   divergence #7), and the same value is recoverable from `C_AB`/`C_BA` in any
   CSV that has them — `J = 1/(1/C_AB + 1/C_BA − 1)`, exact to ~2 ulp, not an
   approximation. **Caveat: CSVs written before 2026-05-14 have the orig
   `containment` placeholder (constant 1.0) instead of the block, so nothing
   is recoverable from them** — that includes all the archived interval-mode
   A/B/C output.

   **Pick by what you need.** `jaccard_similarity_ie` wins on *calibration*
   (MAE vs bedtools 5×10⁻⁴ at p=20 vs 0.15 for `jaccard_similarity`), and it
   is the only one of the two that is comparable across pairs of differing set
   size. Use it for magnitude and for any comparison spanning different set
   sizes; use `jaccard_similarity` only to rank pairs of comparable size.
   Caveat on `_ie`: it is censored at 0 by the `>= 0` clamp (25/90 pairs at
   p=12, all low-J — an exact 0.0 means "clamped or empty", never "measured
   zero") and is uninformative below J ≈ a few/√m.

   **The counter-claim that `jaccard_similarity` wins on *resolution* is
   under revision — don't rely on it.** The measurement behind it (error sd
   0.0014 vs 0.0024 at p=16, J<0.05) used a binned statistic that mixes
   within-bin true-J variance into the register-equality column only, and it
   was taken with the size ratio pinned near 1. See the boxed note in
   `docs/jaccard-definitional-gap.md`; a dedicated experiment is in progress.
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
7. **Second Jaccard column: `jaccard_similarity_ie`** (v0.4.0). Orig emits one
   Jaccard column, computed by register equality. We emit that column
   unchanged — byte-equal, and parity tests still compare it — plus a second
   inclusion-exclusion column immediately after it. Rationale is divergence #2's
   estimator note: the register-equality column is not set Jaccard and is not
   rank-faithful across pairs of differing size, so a run had no
   bedtools-comparable Jaccard in it even though the underlying quantity was
   already being computed for the containment block.

   Derived Python-side in `runner._jaccard_ie_from_containments` from the
   `c_ab`/`c_ba` arrays `pairwise_metrics_hll` already returns — no C++ or
   binding change, mirroring how `_cosketch_from_containments` works.
   Containments are clamped to 1.0 first (Ertl noise can push them a few ulp
   past it), which forces `denom >= 1` and makes divide-by-zero and
   out-of-range results unreachable rather than merely unlikely. A zero
   containment scores `J_ie = 0.0`.

   `tests/test_parity_against_original.py` adds the new name to
   `_PROJECTED_OUT` so it is dropped before comparison — projecting out a
   *Jaccard* column looks wrong at a glance,
   but orig has no counterpart to be unfaithful to and `jaccard_similarity`
   is still compared byte-for-byte.
8. **The Mode D `_with_ends` column family is removed** (v0.6.0). Orig emits
   `jaccard_similarity_with_ends`; we no longer do, and we also dropped the six
   `_ie`/containment/cosketch twins we had added alongside it. Mode D now emits
   exactly one similarity block, on the minimizer HLL. `MinimizerSketch` is a
   single-HLL class: `startend_hll`, `merged()`, `canonicalize_kmer` and
   `_add_kmers_to` are gone.

   **Why.** The merged sketch was `minimizer_hll ∪ startend_hll`, where the
   ends payload was `canonicalize_kmer(s[:k] + s[-k:])` per record with *every*
   sliding k-mer of that 2k string inserted. Four measured problems:

   - **It was an uncontrollable blend, not a distinct quantity.** Because the
     two hash sets are disjoint, the merged Jaccard is exactly the size-weighted
     mediant of the minimizer and flank Jaccards: `J = (I_M+I_E)/(U_M+U_E)`.
     Verified against the shipped register-equality column to 2×10⁻⁴ at p=24
     (k=10, w=100: predicted 0.4095, actual 0.4097); the identity degrades to
     ±0.02 at p≤14 where the chance-agreement floor bends it. The mixing weight
     φ = |E|/(|M|+|E|) is set by (k, w, record count, length) and **never by the
     user** — on a Maurano file it is 0.55 at w=100 and 0.86 at w=500.
   - **k−1 of the k+1 inserted elements are chimeric.** Only indices 0 and k of
     the 2k payload are real sequence k-mers; the rest span the artificial
     `start||end` splice and exist in no genome. Measured on real peak FASTAs
     they are 77–79% of distinct end elements. Canonicalization was applied to
     the whole 2k concatenation, so a record's start-k-mer orientation depended
     on its *other* end.
   - **The flank component measures exact boundary identity, as a cliff.** At
     k=10, w=50, p=24, jittering record boundaries by **1 bp** takes the flank
     Jaccard from 1.000 to 0.130 (minimizer Jaccard: 1.000 → 0.997). On real
     data that bit is always off — two Maurano DHS BEDs share 15 exactly-matching
     intervals out of 168,322 (0.01%).
   - **The stated rationale was false.** The column was justified as rescuing
     records too short to yield minimizers. But the no-minimizer fallback fed
     the whole sequence into *both* HLLs and returned early, so on a
     constant-length sub-threshold corpus the two columns were **bit-identical**
     — it added nothing exactly where it was supposed to help. Dropout begins at
     `L < k + w − 1`.

   Empirically it also just lost: on the 235 existing Maurano Mode D configs,
   scored by Spearman ρ against bedtools, `no_ends` wins 185/49 (exact
   comparison). Filtering to cells that are unsaturated with <10% record
   dropout — an outcome-independent filter — `no_ends` wins 93.2%, and every
   exception is a cell where *both* columns have ρ < 0.39.

   **Caveats worth keeping.** The one place `_with_ends` won on the
   inclusion-exclusion family was 12 cells at w=500, ρ up to 0.9498 vs 0.8950 —
   but those have 72.5–73.1% record dropout, i.e. the regime where the column
   degenerates to record identity, and they sit well below the global optimum
   (k=20, w=20, p=24 → ρ_no_ie = 0.9997 at 1.9% dropout). Also note the arena
   is partly circular: Mode D FASTAs are `bedtools getfasta` extractions from
   one reference, so minimizer Jaccard nearly restates coordinate overlap
   (ρ = 0.9997). Full analysis and the adversarial review that produced it:
   `docs/mode-d-ends-removal.md`.

   Archived CSVs keep the columns; nothing on disk is invalidated. Analyses that
   consumed `jaccard_similarity_with_ends` were re-pointed to
   `jaccard_similarity` — see that doc for the list.

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
