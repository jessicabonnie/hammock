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

`hammock queries.txt refs.txt --mode A` produces a CSV of all pairwise
Jaccard estimates:

```
file1,file2,sketch_type,mode,precision,num_hashes,kmer_size,window_size,jaccard_similarity,containment
sample1.bed,ref.bed,hyperloglog,A,18,NA,NA,NA,0.812...,1.0
sample2.bed,ref.bed,hyperloglog,A,18,NA,NA,NA,0.534...,1.0
```

### Modes

| Mode | Input    | What it sketches |
|------|----------|------------------|
| **A** | BED     | Intervals (`chr\tstart\tend` formatted strings) |
| **B** | BED     | Points: every position inside every interval |
| **C** | BED     | Combined intervals + points, with subsampling (`--subA`, `--subB`, `--expA`) |
| **D** | FASTA   | Sliding-window minimizers + canonicalized start/end k-mers |

Modes are auto-detected from file extension: `.fa/.fasta/.fna/.ffn/.faa/.frn`
(plus `.gz` variants) → mode D; everything else defaults to mode A unless
`--subA/--subB/--expA` are set, in which case mode C.

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

## CLI

```
hammock <queries.txt> <refs.txt> --mode {A,B,C,D} [options]

  -p, --precision N         HyperLogLog precision (4..24, default 18)
  --subA F                  Mode C: subsampling rate for intervals
  --subB F                  Mode B/C: subsampling rate for points
  --expA F                  Mode C: power-of-10 multiplier for A intervals
  --mixed-stride            Mode B/C: deterministic interval-independent point subsampling
  -k, --kmer_size N         Mode D: k-mer length (default 8)
  -w, --window_size N       Mode D: window size (default 40)
  --seed N                  Hash seed (default 42; xxh64 throughout)
  --threads N               Default min(8, cpu_count())
  --full-paths              Use full paths in CSV file1/file2 columns instead of basenames
  -o, --outprefix PREFIX    Output prefix (default "hammock")
  --memory-limit-gb F       Soft memory cap (default 0 = disabled)
  --verbose                 Per-file sketching progress on stderr
```

The output filename embeds the parameters (e.g.
`hammock_hll_p18_jaccA_B0.50.csv`) so output from different runs doesn't
collide.

## Examples

```bash
# Mode A on two sets of BED files, 8 threads, full paths in output:
hammock queries.txt refs.txt --mode A -p 18 --threads 8 --full-paths -o results

# Mode C with 30% point subsampling and an A-side multiplier:
hammock queries.txt refs.txt --mode C --subB 0.3 --expA 0.5 -o results

# Mode D, k=10, window=30:
hammock fastas.txt fastas.txt --mode D -k 10 -w 30 -o seq_results

# Quiet sequential run (e.g. inside a Snakemake job):
hammock queries.txt refs.txt --mode A --threads 1
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
