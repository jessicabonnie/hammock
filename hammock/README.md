# Hammock

A Python library for sketching and comparing files containing lists of things, with special capabilities for comparing bedfiles.

## Features

- Multiple sketching methods:
  - HyperLogLog (with Rust acceleration)
  - MinHash
  - Minimizer
  - Exact counting
- Support for different file types:
  - BED files
  - Sequence files
  - General tab-delimited files
- Different comparison modes:
  - Mode A: Compare intervals only
  - Mode B: Compare points only
  - Mode C: Compare both intervals and points
  - Mode D: Compare sequences

## Installation

```bash
pip install hammock
```

For development installation:
```bash
git clone https://github.com/jbonnie1/hammock.git
cd hammock
pip install -e ".[dev]"
```

## Usage

### Python API

```python
from hammock import compare_files, compare_bed_files, compare_sequence_files

# Compare two BED files
results = compare_bed_files("file1.bed", "file2.bed")
print(results)  # {'jaccard': 0.75, 'containment': 0.8, ...}

# Compare two sequence files
results = compare_sequence_files("file1.fa", "file2.fa")
print(results)  # {'jaccard': 0.65, 'containment': 0.7, ...}

# Compare any two files with custom parameters
results = compare_files(
    "file1.txt",
    "file2.txt",
    mode='A',  # Compare intervals only
    sketch_type="hyperloglog",
    precision=12,
    num_hashes=64,
    subA=1.0,
    subB=1.0,
    expA=0.5,
    use_rust=True
)
print(results)
```

### Command Line Interface

The command-line interface supports batch processing of multiple files against a set of reference files.

```bash
# Basic usage
hammock files_to_compare.txt reference_files.txt

# With options
hammock files_to_compare.txt reference_files.txt \
    --mode A \
    --precision 12 \
    --num_hashes 128 \
    --threads 4 \
    --outprefix results
```

#### Input Files

1. `files_to_compare.txt`: A text file containing paths to files to be compared, one per line
2. `reference_files.txt`: A text file containing paths to primary/reference files to compare against, one per line

#### Command Line Options

- `--mode`: Comparison mode (default: 'A')
  - A: Compare intervals only (default for BED/BigBed files)
  - B: Compare points only
  - C: Compare both intervals and points
  - D: Compare sequences (auto-detected for sequence files)

- `--outprefix`, `-o`: Output file prefix (default: "hammock")
- `--precision`, `-p`: Precision for HyperLogLog sketching (default: 12)
- `--num_hashes`, `-n`: Number of hashes for MinHash sketching (default: 128)
- `--subA`: Subsampling rate for intervals (0 to 1, default: 1.0)
- `--subB`: Subsampling rate for points (0 to 1, default: 1.0)
- `--expA`: Power of 10 exponent for A-type intervals (default: 0)
- `--threads`: Number of threads to use (default: None)
- `--kmer_size`: Size of k-mers for sequence sketching (default: 8)
- `--window_size`: Size of sliding window for sequence sketching (default: 40)
- `--seed`: Random seed for hashing (default: 42)
- `--verbose`: Print verbose output
- `--debug`: Enable debug mode
- `--rust`: Use Rust implementation for HyperLogLog when available

Sketch type options (mutually exclusive):
- `--hyperloglog`: Use HyperLogLog sketching
- `--exact`: Use exact counting
- `--minhash`: Use MinHash sketching
- `--minimizer`: Use minimizer sketching

- `--hashsize`: Hash size in bits (32 or 64, default: 32)

#### Output

The command generates tab-delimited output files with the following naming pattern:
- `{outprefix}_{sketch_type}_{parameters}.tsv`

Example output files:
- `hammock_hll_p12_jaccA.tsv`: HyperLogLog with precision 12, mode A
- `hammock_mh_n128_jaccB.tsv`: MinHash with 128 hashes, mode B
- `hammock_mnmzr_jaccC.tsv`: Minimizer, mode C

The output file contains columns:
- `query`: Name of the query file
- `reference`: Name of the reference file
- Similarity metrics (jaccard, containment, etc.)

## License

MIT License

## Citation

If you use Hammock in your research, please cite:

```
@software{hammock2024,
  author = {Bonnie, Jessica},
  title = {Hammock: A library for sketching and comparing files containing lists of things},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/jbonnie1/hammock}
}
``` 