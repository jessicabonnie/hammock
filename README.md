# hammock

Pairwise Jaccard similarity for BED intervals and FASTA sequences, via
HyperLogLog sketches. A clean Python + C++ refactor of the original
[`hammock`](https://github.com/jessicabonnie/hammock); same CLI, faster
sketching, byte-equal output for modes A/B/C/D on the parity-tested paths.

## What it does

Given two text files, each listing one path per line:

```
$ cat queries.txt
sample1.bed
sample2.bed
$ cat refs.txt
ref.bed
```

`hammock queries.txt refs.txt` produces a CSV of all pairwise Jaccard
estimates (BED input → interval mode by default):

```
file1,file2,sketch_type,mode,precision,num_hashes,kmer_size,window_size,jaccard_similarity,containment
sample1.bed,ref.bed,hyperloglog,B,18,NA,NA,NA,0.812...,1.0
sample2.bed,ref.bed,hyperloglog,B,18,NA,NA,NA,0.534...,1.0
```

### Modes

The top-level choice is **interval** mode (compare BED interval sets) vs
**sequence** mode (compare FASTA sequences). Interval mode has three flavors;
`interval` (base-level overlap) is the default.

| `--mode` | Letter | Input | What it sketches |
|----------|--------|-------|------------------|
| **`interval`** / `interval-points` | **B** | BED | Base-level points — every position in every interval (base-pair overlap Jaccard, like `bedtools jaccard`). **Default for BED.** |
| `interval-string` | A | BED | Intervals as exact `chr\tstart\tend` strings |
| `interval-hybrid` | C | BED | Both, with subsampling (`--subA`, `--subB`, `--expA`) |
| **`sequence`** | D | FASTA | Sliding-window minimizers + canonicalized start/end k-mers. **Default for FASTA / `--ref`.** |

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
`hs1`), a **local FASTA path**, or a URL. Keywords are resolved from a local
cache and are **never downloaded during a run** (HPC compute nodes are usually
firewalled). Populate the cache once on a networked (login) node:

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
  --expA F                  interval-hybrid: power-of-10 multiplier for interval strings
  --mixed-stride            interval-points/hybrid: deterministic point subsampling
  -k, --kmer_size N         sequence mode: k-mer length (default 8)
  -w, --window_size N       sequence mode: window size (default 40)
  --seed N                  Hash seed (default 42; xxh64 throughout)
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

The output filename embeds the parameters (e.g.
`hammock_hll_p18_jaccA_B0.50.csv`) so output from different runs doesn't
collide.

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
when measuring Mode B throughput.

## Testing

```bash
pytest tests/
```

The test suite covers:

- **Parity vs orig (modes A/B/C):** byte-equal Jaccard column against
  `hammock-orig` (orig 0.4.0 installed via pipx).
- **Parity vs orig (mode D):** byte-equal CSVs against the conda-env
  `hammock` (Python 3.12 + bioconda `digest`).
- **Threading correctness:** `--threads 4` produces identical CSVs to
  single-threaded for every mode.
- **Functional regressions:** mode B `--subB` actually subsamples, mode B
  `--mixed-stride` is deterministic across runs.
- **Mode D unit tests:** `canonicalize_kmer`, empty-minimizer fallback,
  determinism, self-Jaccard = 1.0.

Parity tests skip automatically if their respective `hammock` binary isn't
on disk.

## License

Same as the upstream `hammock`.
