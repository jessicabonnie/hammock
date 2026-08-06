# hammock

Pairwise Jaccard similarity for BED intervals and FASTA sequences, via
HyperLogLog sketches. A clean Python + C++ refactor of the original
[`hammock`](https://github.com/jessicabonnie/hammock); same CLI, faster
sketching, byte-equal `jaccard_similarity` for modes A/B/C on the
parity-tested paths (we add columns orig has no counterpart for; mode D matches
structurally — see [Testing](#testing)).

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
pairwise Jaccard estimates (BED input → interval mode by default):

```
file1,file2,sketch_type,mode,precision,num_hashes,kmer_size,window_size,jaccard_similarity,jaccard_similarity_ie,containment_AB,containment_BA,cosketch_geom,cosketch_arith,cosketch_max
sample1.bed,ref.bed,hyperloglog,B,18,NA,8,40,0.7111,0.6190,0.7650,0.7643,0.7647,0.7647,0.7650
sample2.bed,ref.bed,hyperloglog,B,18,NA,8,40,0.2481,0.1210,0.3827,0.1503,0.2398,0.2665,0.3827
```

(Values rounded for display; the CSV writes full precision.
`kmer_size`/`window_size` show the sequence-mode defaults `8`/`40` even for
interval runs, where they're unused; `num_hashes` is always `NA`.)

### Output columns

Every row ends with **two Jaccard estimates** plus a containment/cosketch
block, all in `[0, 1]`:

| Column | Meaning |
|--------|---------|
| `jaccard_similarity` | **Register-equality** statistic: the fraction of active HLL registers holding equal values. This is *not* set Jaccard — it is biased upward, and the bias depends on how loaded the sketches are **and** on `\|A\|/\|B\|`. Kept for parity with the original `hammock`. |
| `jaccard_similarity_ie` | **Set Jaccard** via inclusion-exclusion, `\|A ∩ B\| / \|A ∪ B\|` with `\|A ∩ B\| = \|A\| + \|B\| − \|A ∪ B\|`. Comparable to `bedtools jaccard`. Noisier than the column above at low precision *and* low similarity, and censored at 0 (an exact `0.0` means "estimated ≤ 0", not "measured zero"). |
| `containment_AB` | `\|A ∩ B\| / \|A\|` — fraction of side A (file1/LIST1) covered by B |
| `containment_BA` | `\|A ∩ B\| / \|B\|` — fraction of side B (file2/LIST2) covered by A |
| `cosketch_geom` / `cosketch_arith` / `cosketch_max` | geometric / arithmetic / max mean of the two containments |

Which to use: `jaccard_similarity_ie` when you want a value comparable to an
exact tool or to other corpora; `jaccard_similarity` only to rank pairs of
similar size against each other. See
[`docs/jaccard-definitional-gap.md`](docs/jaccard-definitional-gap.md).

**Sequence mode (D)** emits the same block, computed on the FASTA's window
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
| **`interval`** / `interval-points` | **B** | BED | Base-level points — every position in every interval. Compare `jaccard_similarity_ie` against `bedtools jaccard`, not `jaccard_similarity`. **Default for BED.** |
| `interval-string` | A | BED | Intervals as exact `chr\tstart\tend` strings |
| `interval-hybrid` | C | BED | Both, with subsampling (`--subA`, `--subB`, `--expA`) |
| **`sequence`** | D | FASTA | Sliding-window minimizers. **Default for FASTA / `--ref`.** |

The `mode` **column in the CSV keeps the letter** (A/B/C/D) for compatibility;
the names above are just the CLI/`--mode` spelling (letters still accepted).

Modes are auto-detected: `.fa/.fasta/.fna/.ffn/.faa/.frn` (plus `.gz`) or any
`--ref*` flag → **sequence** mode; otherwise → **interval** mode (B), except that
`--subA`/`--expA` (interval-string knobs) select **interval-hybrid** (C).

**BED→FASTA (Mode D from BED input):** pass `--ref`/`--ref1`/`--ref2` and the
two lists are treated as BED files, converted to FASTA with `bedtools getfasta`
against the given reference(s), then compared as sequences. This enables
cross-reference comparisons (e.g. hg38-derived vs mm10-derived peak sequences).
See [Installation](#installation) for the `bedtools`/`samtools` requirements and
[CLI](#cli) for the flags. Note: cross-species Mode D Jaccard reflects shared
k-mer content (repeat/low-complexity-driven), not homology — prefer a larger
`-k` than the default for cross-species runs.

## Installation

```bash
pip install -e . --no-build-isolation
```

This builds the C++ extension (HLL sketching + BED parser) via
`scikit-build-core` + CMake + pybind11 and installs the `hammock` entry point.
Requires CMake ≥ 3.18 and a C++17 compiler. The build uses a vendored copy of
[mindis/hll](https://github.com/mindis/hll); no submodule init needed.

For mode D you also need `digest`:

```bash
# Easiest path: bioconda (Python ≥ 3.9)
conda install -c bioconda digest
```

`digest` is in `[project.optional-dependencies]` as `mode_d`, but bioconda is
the recommended source — the PyPI `digest` is a different package.

### BED→FASTA mode (`--ref`) requirements

Converting BED lists to FASTA (the `--ref`/`--ref1`/`--ref2` flags, see
[CLI](#cli)) shells out to external command-line tools — no extra Python
packages are needed:

| Tool | Needed for | Required? |
|------|------------|-----------|
| **bedtools** (`getfasta`) | extracting sequences from a reference | **required** |
| **samtools** (`faidx`)    | indexing references (`.fai`)          | recommended — falls back to `bedtools` building the index if the reference dir is writable |

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
                            mixed-stride (default; deterministic chr-keyed stride,
                            fastest at low subB), hash-threshold (random gate; use
                            this for byte-for-byte parity with the original hammock),
                            or single-hash. NOTE: the mixed-stride default means subB
                            output is NOT byte-equal to orig unless you pass
                            --subB-method hash-threshold.
  -k, --kmer_size N         sequence mode: k-mer length (default 8)
  -w, --window_size N       sequence mode: window size (default 40)
  --seed N                  HLL ingestion seed (xxh64, default 42)
  --gate-seed N             Seed for the subB sampling-decision hash (xxh32) and the
                            mixed-stride stride assignment (default 31337 = orig
                            hammock; independent of --seed)
  --threads N               Default min(8, cpu_count())
  --full-paths              Use full paths in CSV file1/file2 columns instead of basenames
  -o, --outprefix PREFIX    Output prefix (default "hammock")
  --memory-limit-gb F       Soft memory cap (default 0 = disabled)
  --verbose                 Per-file sketching progress on stderr

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
and `--ref2`). Requires `bedtools` — see [Installation](#installation).

The output filename embeds the parameters so runs with different settings don't
collide — e.g. `-o out` gives `out_hll_p18_jaccB.csv` (interval), or
`out_hll_p18_jaccC_expA0.50_B0.30.csv` with hybrid subsampling. Sequence-mode
outputs also embed `_k<k>_w<w>` (`out_mnmzr_p18_jaccD_k8_w40.csv`), and BED→FASTA
runs prepend `_<ref1>-vs-<ref2>` so cross-reference runs to the same `-o` don't
overwrite each other (`out_hg38-vs-mm10_mnmzr_p18_jaccD_k8_w40.csv`). Run with
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
- **"reference keyword 'hg38' is not cached"** — keywords/URLs are never
  downloaded mid-run. The message prints the exact `hammock fetch-ref hg38
  --ref-cache-dir <dir>` command; run it once on a networked node.
- **"--ref-cache-dir/--fasta-outdir only apply with a reference"** or
  **"--ref … mutually exclusive with --ref1/--ref2"** or **"give both --ref1
  and --ref2"** — reference-flag validation; supply exactly `--ref` (both lists)
  or both of `--ref1`/`--ref2`.
- **Sequence mode raises "requires the `digest` module"** — `digest` isn't
  importable. Verify `python -c "import digest"`; if it fails with `GLIBCXX_…
  not found`, your `_core` extension has a stale RPATH pulling in an old
  libstdc++ — rebuild with the environment's own compiler (see the build notes).

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
benchmarking — same algorithms as `hammock`, no Python in the loop. Useful
when measuring Mode B throughput. It writes a **tab**-separated file.

Since **0.7.0** its defaults match the Python CLI's: mode **B** for BED input,
and the full similarity block — `jaccard_similarity`, `jaccard_similarity_ie`,
`containment_AB`/`containment_BA`, and the three `cosketch_*` columns.

```bash
hammock-cpp queries.txt refs.txt -p 20 -o out      # -> out_hll_p20_jaccB.csv
```

Those values are **bit-for-bit identical** to the Python CLI's on the same input
(`tests/test_hammock_cpp_metrics.py` asserts exact equality).

Pass **`--no-metrics`** for timing runs. It emits only `query`, `reference`,
`jaccard_similarity` and tags the file `_j3`, so the two shapes never collide.
The block costs a union plus two cardinality estimates per pair — worth skipping
when you are measuring throughput, not similarity:

```bash
hammock-cpp queries.txt refs.txt -p 20 --no-metrics -o out   # -> out_hll_p20_jaccB_j3.csv
```

`--metrics` is still accepted, and is now a no-op. `--version` prints the
version to stdout.

Upgrading from ≤ 0.6.x: a bare invocation used to mean mode A with 3 columns.
Pass `--mode A` if you relied on that, and note that `--peak-height` and its
BagMinHash backend were removed — they were never wired into either CLI.

## Testing

```bash
pytest tests/
```

The test suite covers:

- **Parity vs orig (modes A/B/C):** byte-equal Jaccard column against
  `hammock-orig` (orig 0.4.0 installed via pipx).
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
  determinism, self-Jaccard = 1.0.

Parity tests skip automatically if their respective `hammock` binary isn't
on disk.

## License

Same as the upstream `hammock`.
