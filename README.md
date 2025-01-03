# hammock

A tool for comparing BED files using sketching algorithms.

## Installation

```bash
git clone https://github.com/jessicabonnie/hammock.git
cd hammock
pip install -e .
```
`hammock`'s mode D makes use of `Digest`, a C++ library that supports various sub-sampling schemes for $k$-mers in DNA sequences.
```bash
cd ..
git clone https://github.com/VeryAmazed/digest.git
cd digest

meson setup --prefix=$(pwd)/build --buildtype=release build
meson install -C build
pip install .
```

## Usage

Calculate pairwise Jaccard similarities between lists of BED or sequence files.
Both input files are text files containing a path per line to the files to be compared. For computational purposes, `<primary_paths_file>` should be the shorter of the two file lists.
```bash
hammock <filepaths_file> <primary_file> --mode {A|B|C|D} [options]
```


Modes:
- `A`: Compare intervals only
- `B`: Compare points only
- `C`: Compare both intervals and points
- `D`: Compare sequences

Options:
- `--hyperloglog`: Use HyperLogLog sketching (default)
- `--minhash`: Use MinHash sketching
- `--exact`: Use exact counting
- `--precision`, `-p`: Precision for HyperLogLog (default: 8)
- `--num_hashes`, `-n`: Number of hashes for MinHash (default: 128)
- `--subsample`: Subsampling rate (default: 1.0)
- `--balance`: When set, type A values will be subsampled by (1-subsample)
- `--debug`: Output debug information including sparsity comparisons

### Output

Results are written to CSV files with the following format:
```bed1, bed2, sketch_type, mode, num_hashes, precision, subsample, balance, jaccardfunc, jaccardcalc, intersect, union```


## Examples

```bash
# Compare intervals using HyperLogLog
hammock files.txt primary.txt --mode A --precision 12
# Compare points using MinHash with subsampling
hammock files.txt primary.txt --mode B --minhash --num_hashes 256 --subsample 0.1
```

## Testing

In order to run the tests requiring BigBed files, you need to have `bedToBigBed` installed and in your PATH.
```bash
#conda install -c bioconda ucsc-bedtools
```

Run all tests:
```bash
python3 -m pytest hammock/tests
```

### Test Files
- `test_hyperloglog.py`: Tests HyperLogLog sketching with various precisions and set sizes
- `test_minhash.py`: Tests MinHash sketching with different numbers of hash functions
- `benchmark_sketches.py`: Performance benchmarks comparing different sketching methods

### Test Cases
Tests include:
- Perfect overlap scenarios
- No overlap scenarios
- Partial overlap with different set sizes
- Sparse vs dense comparisons
- Edge cases with very different sparsity levels

### Running Specific Tests

```bash
# Run HyperLogLog tests
python -m pytest hammock/tests/test_hyperloglog.py
# Run MinHash tests
python -m pytest hammock/tests/test_minhash.py
# Run benchmarks
python -m pytest hammock/tests/benchmark_sketches.py
```