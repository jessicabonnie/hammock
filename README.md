# hammock

Pairwise sketch similarity and estimated set Jaccard for BED intervals and FASTA sequences, via
HyperLogLog sketches. A clean Python + C++ refactor of the original
[`hammock`](https://github.com/jessicabonnie/hammock); same CLI, faster
sketching, byte-equal `reg_eq_similarity` (orig's `jaccard_similarity`,
renamed — see [Output shapes and columns](#output-shapes-and-columns)) for
modes A/B/C on the parity-tested paths (we add columns orig has no
counterpart for; mode D matches structurally — see [Testing](#testing)).

## What it does

Given two text files, each listing one path per line:

```
$ cat queries.txt
sample1.bed
sample2.bed
$ cat refs.txt
ref.bed
```

`hammock queries.txt refs.txt` produces a **comma-separated** CSV of all
pairwise similarity summaries (BED input → interval mode by default). The bare
invocation writes the `_ie` shape — just the set-Jaccard column (see
[Output shapes and columns](#output-shapes-and-columns) below):

```
file1,file2,sketch_type,mode,precision,num_hashes,kmer_size,window_size,jaccard_similarity_ie
sample1.bed,ref.bed,hyperloglog,B,18,NA,8,40,0.6190
sample2.bed,ref.bed,hyperloglog,B,18,NA,8,40,0.1210
```

(Values rounded for display; the CSV writes full precision.
`kmer_size`/`window_size` show the sequence-mode defaults `8`/`40` even for
interval runs, where they're unused; `num_hashes` is always `NA`.)

### Output shapes and columns

Every invocation writes exactly one of **three mutually exclusive shapes**,
selected by flag, and the output filename is tagged to match so the three
never collide on one path (see [CLI](#cli)):

| Flag | Filename tag | Columns written | Cost |
|------|--------------|------------------|------|
| *(none — default)* | `_ie` | `jaccard_similarity_ie` | needs the same fused union pass as `--metrics` — **not** the cheap arm |
| `--register-equality` / `--re` | `_re` | `reg_eq_similarity` | the cheap arm on `hammock-cpp` (skips the union/containment pass entirely); on the Python CLI it is **not** cheaper than `--metrics` — the Python binding always computes the fused pass regardless of shape |
| `--metrics` | `_full` | `reg_eq_similarity`, `jaccard_similarity_ie`, `containment_AB`, `containment_BA`, `cosketch_geom`, `cosketch_arith`, `cosketch_max` (7 columns) | full cost |

All values are in `[0, 1]`. Column meanings:

| Column | Meaning |
|--------|---------|
| `reg_eq_similarity` | **Register-equality similarity** (renamed from `jaccard_similarity` — see below): the fraction of active HLL registers holding equal values. Despite the *original* `jaccard_similarity` name, this is *not* set Jaccard — it is biased upward, and the bias depends on how loaded the sketches are **and** on `\|A\|/\|B\|`. |
| `jaccard_similarity_ie` | **Set Jaccard** via inclusion-exclusion, `\|A ∩ B\| / \|A ∪ B\|` with `\|A ∩ B\| = \|A\| + \|B\| − \|A ∪ B\|`. Comparable to `bedtools jaccard`. Noisier than the column above at low precision *and* low similarity, and censored at 0 (an exact `0.0` means "estimated ≤ 0", not "measured zero"). |
| `containment_AB` | `\|A ∩ B\| / \|A\|` — fraction of side A (file1/LIST1) covered by B |
| `containment_BA` | `\|A ∩ B\| / \|B\|` — fraction of side B (file2/LIST2) covered by A |
| `cosketch_geom` / `cosketch_arith` / `cosketch_max` | geometric / arithmetic / max mean of the two containments |

Which to use: `jaccard_similarity_ie` (the bare default) when you want a value
comparable to an exact tool or to other corpora; `reg_eq_similarity` only to
rank pairs of similar size against each other. See
[`docs/jaccard-definitional-gap.md`](docs/jaccard-definitional-gap.md).

**Renamed in v0.9.0 (`docs/seed-jaccard-reg-eq-rename.md`): `jaccard_similarity`
is now `reg_eq_similarity`**, and the `register_equality_similarity` column
(a literal duplicate of it, added in v0.8.0 so `--re`/`--metrics` output was
self-describing without relying on column position) is dropped — `reg_eq_similarity`
now does that on its own. Archived CSVs written before this rename still carry
the old `jaccard_similarity` name; code that reads them should prefer
`reg_eq_similarity` if present and fall back to `jaccard_similarity` if not.

**Interval-hybrid (C)** inserts three parameter columns — `subA`,`subB`,`expA` —
between `window_size` and the similarity columns, regardless of shape; modes A
and B do not.

**Sequence mode (D)** emits the same shape, computed on the FASTA's window
minimizers. Every Mode D row also ends with trailing **`ref1`,`ref2`** columns
(the reference each list was extracted against; `NA` for plain-FASTA runs).

Sequence mode used to emit a second `_with_ends` copy of the block, computed on
the minimizers plus each record's start/end k-mers. Those seven columns were
removed in v0.6.0: the merged sketch was a blend of minimizer similarity and
*exact record-boundary identity* at a mixing weight the user could not set, and
the boundary term is destroyed by 1 bp of coordinate jitter. See divergence #8
in `CLAUDE.md` and [`docs/mode-d-ends-removal.md`](docs/mode-d-ends-removal.md).

### Modes

The top-level choice is **interval** mode (compare BED interval sets) vs
**sequence** mode (compare FASTA sequences). Interval mode has three flavors;
`interval` (base-level overlap) is the default.

| `--mode` | Letter | Input | What it sketches |
|----------|--------|-------|------------------|
| **`interval`** / `interval-points` | **B** | BED | Base-level points — every position in every interval. Compare `jaccard_similarity_ie` against `bedtools jaccard`, not `reg_eq_similarity`. **Default for BED.** |
| `interval-string` | A | BED | Intervals as exact `chr\tstart\tend` strings |
| `interval-hybrid` | C | BED | Both, with subsampling (`--subA`, `--subB`, `--expA`) |
| **`sequence`** | D | FASTA | Sliding-window minimizers. **Default for FASTA / `--ref`.** |

The `mode` **column in the CSV keeps the letter** (A/B/C/D) for compatibility;
the names above are just the CLI/`--mode` spelling (letters still accepted).

Modes are auto-detected: `.fa/.fasta/.fna/.ffn/.faa/.frn` (plus `.gz`) or any
`--ref*` flag → **sequence** mode; otherwise → **interval** mode (B), except that
`--subA`/`--expA` (interval-string knobs) select **interval-hybrid** (C).

Interval modes read **plain-text** BED only: the parser does not decompress
`.gz` or decode binary BigBed. hammock detects the common compressed and
binary formats by magic bytes and **refuses to run** (exit 2) with a message
naming the conversion command — decompress or convert first. The check is
deliberately conservative, so an unrecognized binary format can still slip
through and sketch to nothing, scoring `0.0` against every file including
itself. (Sequence mode does handle `.gz` FASTA.)

**BED→FASTA (Mode D from BED input):** pass `--ref`/`--ref1`/`--ref2` and the
two lists are treated as BED files, converted to FASTA with `bedtools getfasta`
against the given reference(s), then compared as sequences. This enables
cross-reference comparisons (e.g. hg38-derived vs mm10-derived peak sequences).
See [Installation](#installation) for the `bedtools`/`samtools` requirements and
[CLI](#cli) for the flags. Note: cross-species Mode D similarity reflects shared
k-mer content (repeat/low-complexity-driven), not homology — prefer a larger
`-k` than the default for cross-species runs.

## Installation

```bash
pip install -e . --no-build-isolation
```

This builds the C++ extension (HLL sketching + BED parser) via
`scikit-build-core` + CMake + pybind11 and installs the `hammock` entry point.
Requires Python ≥ 3.8, CMake ≥ 3.18, and a C++17 compiler with OpenMP. The build
uses a vendored copy of [mindis/hll](https://github.com/mindis/hll); no submodule
init needed.

Build with the compiler that belongs to the environment you will run in. If a
module system or a stale toolchain puts a different `CC`/`CXX` in front, the
extension can bake in an RPATH to an older `libstdc++`, which surfaces later as
`GLIBCXX_… not found` on `import digest` (see [Troubleshooting](#troubleshooting)).
Wipe `build/` when you change compilers — CMake caches the old one.

For mode D you also need `digest` and `xxhash`:

```bash
# Easiest path: bioconda (Python ≥ 3.9)
conda install -c bioconda digest
pip install xxhash
```

Both are in `[project.optional-dependencies]` as `mode_d`, so
`pip install -e '.[mode_d]'` covers them — but bioconda is the recommended
source for `digest`, because the PyPI package of that name is a different
project. `environment.yml` lists both (`digest`, `python-xxhash`).

### BED→FASTA mode (`--ref`) requirements

Converting BED lists to FASTA (the `--ref`/`--ref1`/`--ref2` flags, see
[CLI](#cli)) shells out to external command-line tools — no extra Python
packages are needed:

| Tool | Needed for | Required? |
|------|------------|-----------|
| **bedtools** (`getfasta`) | extracting sequences from a reference | **required** |
| **samtools** (`faidx`)    | indexing references (`.fai`)          | recommended — falls back to `bedtools` building the index if the reference dir is writable; **required** for `fetch-ref`, for a gzipped local reference, and when the reference lives in a read-only directory |

Make sure both are on `PATH` before running. On an HPC cluster with environment
modules:

```bash
ml bedtools samtools          # or: ml UCSC_Genome_Browser/2021 (bundles bedtools)
```

Reference genomes are supplied as a **keyword** (`hg38`, `mm10`, `hg19`, `mm39`,
`hs1`), a **local FASTA path**, or a **URL**. Local paths are used directly;
keywords and URLs are resolved from a local cache and are **never downloaded
during a run** (HPC compute nodes are usually firewalled) — a run errors if they
aren't cached yet. Populate the cache once on a networked (login) node:

```bash
hammock fetch-ref hg38 --ref-cache-dir /shared/refs    # downloads + indexes once
```

`fetch-ref` needs `samtools` on `PATH` and network access (http(s) only).

## CLI

```
hammock <queries.txt> <refs.txt> [--mode MODE] [options]

  --mode MODE               interval (default, BED) | sequence (FASTA/--ref);
                            advanced: interval-string (A), interval-points (B),
                            interval-hybrid (C). Letters A/B/C/D also accepted.
  -p, --precision N         HyperLogLog precision (4..24, default 18)
  --subA F                  interval-hybrid: subsampling rate for interval strings
  --subB F                  interval-points/hybrid: subsampling rate for points
  --expA F                  interval-hybrid: power-of-10 exponent weighting interval
                            strings (e.g. 2 → each string counts ×100; default 0)
  --subB-method M           interval-points/hybrid: point-sampling method —
                            mixed-stride (default; deterministic chromosome-anchored
                            fractional-interval sampling that mixes adjacent integer
                            gap lengths), mixed-stride-v1 (legacy one-stride-per-
                            chromosome behavior), mixed-stride-v2 (alias for the
                            default), hash-threshold (random gate; use this for byte-
                            for-byte parity with the original hammock), or single-hash.
                            NOTE: the mixed-stride default means subB
                            output is NOT byte-equal to orig unless you pass
                            --subB-method hash-threshold.
  -k, --kmer_size N         sequence mode: k-mer length (default 8)
  -w, --window_size N       sequence mode: window size (default 40)
  --sequence-hll-hash M     sequence HLL hash: rehash-selector64 (default) or
                            legacy-selector32 for reproducing pre-v0.12 results
  --seed N                  HLL ingestion seed (xxh64, default 42)
  --gate-seed N             Seed for the subB sampling-decision hash (xxh32) and the
                            mixed-stride grid phase (default 31337 = orig hammock;
                            independent of --seed)
  --threads N               Default: min(8, cpu_count()) for interval modes, 1 for
                            sequence mode (D), whose sketch loop is GIL-bound — an
                            explicit --threads > 1 there is honored but usually slower
  --full-paths              Write normalized input paths in the CSV file1/file2 columns
                            instead of basenames
  --metrics                 Emit the full 7-column block (reg_eq_similarity,
                            jaccard_similarity_ie, containment_AB/BA, cosketch_*);
                            tags the output _full. Mutually exclusive with
                            --register-equality/--re.
  --register-equality, --re Emit only reg_eq_similarity — the cheap
                            register-equality-only arm on hammock-cpp (not on
                            the Python CLI, see Output shapes); tags the
                            output _re. Mutually exclusive with --metrics.
                            Default (neither flag): emit jaccard_similarity_ie alone,
                            tagged _ie.
  -o, --outprefix PREFIX    Output prefix (default "hammock")
  --memory-limit-gb F       Soft memory cap in GiB (default 0 = disabled)
  --verbose                 Per-file sketching progress on stderr
  --version                 Print the version and exit

  BED→FASTA (Mode D) — treat the two lists as BED files, convert via bedtools:
  --ref REF                 Reference for BOTH lists (keyword | local FASTA | URL)
  --ref1 REF                Reference for list 1 (mutually exclusive with --ref)
  --ref2 REF                Reference for list 2 (mutually exclusive with --ref)
  --ref-cache-dir DIR       Cached/indexed references (default: $HAMMOCK_REF_CACHE
                            or ~/.hammock/refs)
  --fasta-outdir DIR        Keep generated FASTAs (default: temp dir, auto-cleaned)

hammock fetch-ref <keyword|url> [--ref-cache-dir DIR] [--force]
                            Download + index a reference into the cache
                            (run once on a networked/login node)
```

Any `--ref*` flag reinterprets the two positional lists as BED files, forces
Mode D, and adds trailing `ref1`/`ref2` columns to the CSV (they are `"NA"` for
plain FASTA runs). Both list references must be given (`--ref`, or both `--ref1`
and `--ref2`); combining a reference flag with an interval `--mode` is an error.
Requires `bedtools` — see [Installation](#installation).

The output filename embeds the parameters so runs with different settings don't
collide — e.g. `-o out` gives `out_hll_p18_jaccB_ie.csv` (interval, bare
default), or `out_hll_p18_jaccC_expA0.50_B0.30_ie.csv` with hybrid subsampling.
Sequence-mode outputs also embed `_k<k>_w<w>_rehash-selector64`
(`out_mnmzr_p18_jaccD_k8_w40_rehash-selector64_ie.csv`); explicit legacy
hashing retains the pre-v0.12 `_k<k>_w<w>` filename. BED→FASTA runs insert
`_<ref1>-vs-<ref2>` right after your prefix, so cross-reference runs to the
same `-o` don't overwrite each other
(`out_hg38-vs-mm10_mnmzr_p18_jaccD_k8_w40_rehash-selector64_ie.csv`). The trailing `_ie`/`_re`/
`_full` always names the output shape (see
[Output shapes and columns](#output-shapes-and-columns)) — `--register-equality`/
`--re` gives `..._re.csv`, `--metrics` gives `..._full.csv`. Run with
`--verbose` to print the exact path.

## Examples

```bash
# Interval mode (default) on two sets of BED files, 8 threads, full paths:
hammock queries.txt refs.txt -p 18 --threads 8 --full-paths -o results

# Interval-string flavor (exact interval boundaries) explicitly:
hammock queries.txt refs.txt --mode interval-string -o results

# Interval-hybrid with 30% point subsampling and an interval-string multiplier:
hammock queries.txt refs.txt --mode interval-hybrid --subB 0.3 --expA 0.5 -o results

# Sequence mode, k=10, window=30:
hammock fastas.txt fastas.txt --mode sequence -k 10 -w 30 -o seq_results

# BED→FASTA: compare two BED-peak lists as sequences, both on hg38
# (requires `ml bedtools`; hg38 must be cached — see fetch-ref below):
hammock peaks1.txt peaks2.txt --ref hg38 --ref-cache-dir /shared/refs -o seq_out

# Cross-reference: list 1 from hg38, list 2 from mm10, using local FASTAs:
hammock human_peaks.txt mouse_peaks.txt \
        --ref1 /refs/hg38.fa --ref2 /refs/mm10.fa -o cross_out

# Populate the reference cache once on a login node:
hammock fetch-ref hg38 --ref-cache-dir /shared/refs

# Quiet sequential run (e.g. inside a Snakemake job):
hammock queries.txt refs.txt --mode interval --threads 1

# Full metrics block (containment + cosketch columns), not just jaccard_similarity_ie:
hammock queries.txt refs.txt --metrics -o results   # -> results_hll_p18_jaccB_full.csv
```

## Troubleshooting

hammock fails loudly rather than emitting misleading zeros. Common errors and
their fixes:

- **"looks like a BED/FASTA file, not a list of paths"** — the two positional
  arguments are *text files listing one input path per line*, not the inputs
  themselves. Build them first: `ls *.bed > queries.txt`.
- **BED→FASTA gives an error about intervals "not found in the fasta"** — your
  BED chromosome names don't match the reference (`chr1` vs `1`). This is caught
  as a hard error, **not** a silent Jaccard≈0. Use a reference whose naming
  matches your BEDs (or re-name the BED chromosomes).
- **"Query input … is gzip-compressed"** (or zstd/bzip2/xz/BigBed/BigWig) —
  interval modes read plain-text BED only, so hammock rejects the file instead
  of sketching it to nothing. The message names a conversion command; run it
  and point the list at the converted file.
- **"Reference keyword 'hg38' (→ hg38) is not cached in <dir>"** — keywords/URLs
  are never downloaded mid-run. The message prints the exact `hammock fetch-ref
  hg38 --ref-cache-dir <dir>` command; run it once on a networked node.
- **"--ref-cache-dir/--fasta-outdir only apply with a reference flag"** or
  **"--ref is mutually exclusive with --ref1/--ref2"** or **"BED→FASTA mode needs
  a reference for both lists"** — reference-flag validation; supply exactly
  `--ref` (both lists) or both of `--ref1`/`--ref2`.
- **Sequence mode raises "requires the `digest` module"** — `digest` isn't
  importable. Verify `python -c "import digest"`; if it fails with `GLIBCXX_…
  not found`, your `_core` extension has a stale RPATH pulling in an old
  libstdc++ — remove `build/` and rebuild with the environment's own compiler
  (see [Installation](#installation)).

## Layout

```
python/hammock/        # the Python program (CLI + orchestration + Mode D)
cpp/                   # C++ extension source
  include/hammock/     #   public headers
  src/                 #   sketches, BED parser, mode procs, mixed-stride
  app/                 #   `hammock-cpp` benchmark binary (no Python)
bindings/_core.cpp     # pybind11 module
extern/hll/            # vendored mindis/hll
tests/                 # pytest: unit + parity tests
```

A standalone `hammock-cpp` binary is built alongside the wheel for max-speed
benchmarking — same algorithms as `hammock` for modes A/B/C, no Python in the
loop (mode D needs `digest` and stays on the Python side). Useful when measuring
Mode B throughput. It writes a **tab**-separated file.

**It is not on `PATH`.** Invoke it by full path: `build/<wheel-tag>/hammock-cpp`
after a local build, or `<site-packages>/bin/hammock-cpp` from an installed
wheel.

Since **0.7.0** its mode default matches the Python CLI's: mode **B** for BED
input. Since **0.8.0** it also matches the Python CLI's output-shape contract
(see [Output shapes and columns](#output-shapes-and-columns)) — three mutually
exclusive shapes, tagged filenames, none bare:

```bash
hammock-cpp queries.txt refs.txt -p 20 -o out      # -> out_hll_p20_jaccB_ie.csv (default: jaccard_similarity_ie only)
```

Those values are **bit-for-bit identical** to the Python CLI's on the same input,
for every shape (`tests/test_hammock_cpp_metrics.py` asserts exact equality).

Pass **`--register-equality`**/**`--re`** for timing runs. It emits only
`query`, `reference`, `reg_eq_similarity` and tags the file `_re`. Unlike
the default/`--metrics` shapes, this one **skips
the union pass entirely on `hammock-cpp`** — it is the cheap arm, worth using
when you are measuring throughput, not similarity. (On the Python CLI the
same flag is **not** cheaper: the binding always computes the fused union
pass regardless of shape, so `--re` there saves output columns, not compute.)

```bash
hammock-cpp queries.txt refs.txt -p 20 --register-equality -o out   # -> out_hll_p20_jaccB_re.csv
```

Pass **`--metrics`** for the full 7-column block — `reg_eq_similarity`,
`jaccard_similarity_ie`, `containment_AB`/`containment_BA`, and the three
`cosketch_*` columns — tagged `_full`:

```bash
hammock-cpp queries.txt refs.txt -p 20 --metrics -o out   # -> out_hll_p20_jaccB_full.csv
```

`--no-metrics` was removed in v0.8.0, not aliased — pass `--register-equality`/
`--re` instead. `--version` prints the version to stdout.

Upgrading from ≤ 0.6.x: a bare invocation used to mean mode A with 3 columns.
Pass `--mode A` if you relied on that, and note that `--peak-height` and its
BagMinHash backend were removed — they were never wired into either CLI.
Upgrading from 0.7.x: a bare invocation used to mean mode B with the full
7-column block (untagged); it now means mode B with `jaccard_similarity_ie`
alone, tagged `_ie`. Pass `--metrics` for the closest equivalent to the old
shape — now 7 columns, tagged `_full` (it briefly gained an 8th,
`register_equality_similarity`, in 0.8.0; that column was dropped again in
0.9.0, see below).
Upgrading from 0.8.x: `jaccard_similarity` is renamed `reg_eq_similarity`,
and `register_equality_similarity` (the 0.8.0-era literal duplicate of it)
is gone — `reg_eq_similarity` is self-describing on its own now, so `re`
drops from 4 columns to 3 and `full` from 8 to 7
(`docs/seed-jaccard-reg-eq-rename.md`).

The current Mode D parallelism limitations and the preferred C++ implementation
boundary are documented in
[`docs/mode-d-parallelism.md`](docs/mode-d-parallelism.md).

## Testing

```bash
pytest tests/
```

The test suite covers:

- **Parity vs orig (modes A/B/C):** byte-equal `reg_eq_similarity` column
  (compared by position against orig's frozen `jaccard_similarity`, since
  v0.9.0 the two are no longer the same literal name) against `hammock-orig`
  (the upstream `hammock`, installed via pipx under that name).
- **Structural parity vs orig (mode D):** matching structural columns and
  well-formed similarity values against the conda-env `hammock` (Python 3.12 +
  bioconda `digest`). Mode D is **not** byte-equal by design — it ingests each
  minimizer's raw 64-bit hash rather than the original's decimal-digit-k-mer
  slow path, so the values differ (e.g. 0.75 vs 0.79 on the fixtures).
- **Threading correctness:** `--threads 4` produces identical CSVs to
  single-threaded for every mode.
- **Functional regressions:** mode B `--subB` actually subsamples, and
  `--subB-method mixed-stride` (the default) is deterministic across runs.
- **Mode D unit tests:** exact CSV header, empty-minimizer fallback,
  determinism, self-similarity = 1.0.

Parity tests skip automatically if their respective `hammock` binary isn't
on disk.

## License

Same as the upstream `hammock`.
