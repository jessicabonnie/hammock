# hammock

A tool for comparing BED files using sketching algorithms.

## Installation

```bash
git clone https://github.com/jessicabonnie/hammock.git
cd hammock
pip install -e .
```

## Usage

Calculate Jaccard similarities between BED files:
```bash
hammock <filepaths_file> <primary_file> --mode {A|B|C} [options]
```


Modes:
- `A`: Compare intervals only
- `B`: Compare points only
- `C`: Compare both intervals and points

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
```bed1,bed2,sketch_type,mode,num_hashes,precision,subsample,balance,jaccardfunc,jaccardcalc,intersect,union```


## Examples

```bash
# Compare intervals using HyperLogLog
hammock files.txt primary.txt --mode A --precision 12
# Compare points using MinHash with subsampling
hammock files.txt primary.txt --mode B --minhash --num_hashes 256 --subsample 0.1
```

## Testing

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