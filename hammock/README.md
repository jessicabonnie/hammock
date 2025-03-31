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

### Basic Usage

```python
from hammock import compare_files, compare_bed_files, compare_sequence_files

# Compare two BED files
results = compare_bed_files("file1.bed", "file2.bed")
print(results)

# Compare two sequence files
results = compare_sequence_files("file1.fa", "file2.fa")
print(results)

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

```bash
# Compare BED files
hammock compare-bed file1.bed file2.bed

# Compare sequence files
hammock compare-seq file1.fa file2.fa

# Compare any files with custom parameters
hammock compare file1.txt file2.txt --mode A --sketch-type hyperloglog --precision 12
```

## Parameters

### Common Parameters

- `sketch_type`: Type of sketch to use
  - "hyperloglog" (default)
  - "minhash"
  - "minimizer"
  - "exact"
- `precision`: Precision for HyperLogLog sketching (default: 12)
- `num_hashes`: Number of hashes for MinHash sketching (default: 64)

### BED File Parameters

- `subA`: Subsampling rate for intervals (0 to 1, default: 1.0)
- `subB`: Subsampling rate for points (0 to 1, default: 1.0)
- `expA`: Power of 10 exponent for A-type intervals (default: 0.5)
- `use_rust`: Whether to use Rust implementation for HyperLogLog (default: True)

### Sequence File Parameters

- `kmer_size`: Size of k-mers for sequence sketching (default: 8)
- `window_size`: Size of sliding window for sequence sketching (default: 40)

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