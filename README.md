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