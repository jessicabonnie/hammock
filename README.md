# hammock

A tool for comparing BED files using sketching algorithms.

## Installation

### Basic Installation (Pure Python)

```bash
git clone https://github.com/jessicabonnie/hammock.git
cd hammock
pip install -e .
```

### Optimized Installation (Recommended)

For **2-5x performance improvement** with FastHyperLogLog and optional Cython acceleration:

```bash
# Install with Cython acceleration support
git clone https://github.com/jessicabonnie/hammock.git
cd hammock
pip install -e ".[fast]"
```

This will:
- ✅ Enable FastHyperLogLog with Cython acceleration (2-5x speedup)
- ✅ Automatically fall back to pure Python if Cython build fails
- ✅ Provide identical results with significantly better performance
- ✅ Work transparently - existing code runs faster without changes

**Alternative installation methods:**
```bash
# Option 1: Install from PyPI (when available)
pip install "hammock[fast]"

# Option 2: Development installation with all tools
pip install -e ".[fast,dev]"
python setup.py build_ext --inplace

# Option 3: Standard installation without acceleration
pip install hammock  # Pure Python, works everywhere
```

### Digest for $k$-mers

`hammock`'s mode D makes use of `Digest`, a C++ library that supports various sub-sampling schemes for $k$-mers in DNA sequences. `Digest` is now available on bioconda:    
```bash
conda install -c bioconda digest
```

If you would like to install `Digest` from source, follow the instructions [here](https://github.com/VeryAmazed/digest.git).

### Performance Comparison

Hammock offers multiple implementation options for different performance needs:

| Implementation | Installation | Performance | When to Use |
|----------------|--------------|-------------|-------------|
| **Pure Python** | `pip install hammock` | Baseline (1x) | ✅ Simple installation<br>✅ Works everywhere<br>⚠️ Slower for large datasets |
| **FastHyperLogLog** | `pip install "hammock[fast]"` | **2-5x faster** | 🚀 **Recommended for most users**<br>✅ Significant speedup<br>✅ Automatic fallback<br>✅ Easy installation |
| **Rust HyperLogLog** | Build with maturin | **5-10x faster** | 🔥 Maximum performance<br>⚠️ Requires Rust build tools<br>⚠️ More complex setup |

**Performance results** (based on comprehensive testing):
- **Small datasets** (100-1000 items): FastHyperLogLog provides ~2x speedup
- **Medium datasets** (1000-10000 items): FastHyperLogLog provides ~3x speedup  
- **Large datasets** (10000+ items): FastHyperLogLog provides ~2-5x speedup
- **Batch operations**: Up to 5x speedup with optimized batch processing

**Recommendation**: Install with `pip install "hammock[fast]"` for the best balance of performance and ease of use.


## Usage

Calculate pairwise Jaccard similarities between lists of BED or sequence files.
Both input files are text files containing a path per line to the files to be compared. For computational purposes, `<primary_paths_file>` should be the shorter of the two file lists.
```bash
hammock <filepaths_file> <primary_file> --mode {A|B|C|D} [options]
```

### Modes and Auto-Detection

Hammock supports different comparison modes and includes automatic mode detection:

- `A`: Compare intervals only (default for BED files)
- `B`: Compare points only
- `C`: Compare both intervals and points
- `D`: Compare sequences (auto-detected for sequence files)

### Supported File Formats

Interval Files:
- BED format (`.bed`): Tab-delimited format with chromosome, start, and end positions
- BigBed format (`.bb`): Binary indexed version of BED format, useful for large datasets
- BigWig format (`.bw`): Binary format for continuous data, intervals created from runs of non-zero values
- Any tab-delimited file with at least 3 columns (chr, start, end) in BED-style format

Sequence Files (triggers mode D):
- FASTA format (`.fa`, `.fasta`)
- Nucleotide FASTA (`.fna`, `.ffn`)
- Amino acid FASTA (`.faa`)
- Non-coding RNA FASTA (`.frn`)

Mode Detection:
- If any file path listed in `<filepaths_file>` has a sequence file extension (`.fa`, `.fasta`, `.fna`, `.ffn`, `.faa`, or `.frn`), hammock automatically switches to mode D
- Using C-specific parameters (`--subA`, `--subB`, or `--expA`) automatically switches to mode C
- Otherwise, defaults to mode A for BED files

For example, if your `filepaths_file` contains:
```
data/sample1.fasta
data/sample2.fasta
data/sample3.fa
```
`hammock` will automatically use mode D for sequence comparison, regardless of the specified mode.

### Options

Sketching Options:
- `--hyperloglog`: Use HyperLogLog sketching (default, automatically uses FastHyperLogLog when available)
- `--minhash`: Use MinHash sketching
- `--exact`: Use exact counting
- `--precision`, `-p`: Precision for HyperLogLog (default: 8)
- `--num_hashes`, `-n`: Number of hashes for MinHash (default: 128)
- `--hashsize`: Hash size in bits for HyperLogLog (32 or 64, default: 64)


**Performance Note**: When you install with `pip install "hammock[fast]"`, HyperLogLog sketching automatically uses FastHyperLogLog with Cython acceleration for 2-5x better performance. This works transparently with all existing commands.

Mode C Options:
- `--subA`: Subsampling rate for intervals (default: 1.0)
- `--subB`: Subsampling rate for points (default: 1.0)
- `--expA`: Power of 10 exponent to multiply contribution of A-type intervals (default: 0)

Mode D Options:
- `--window_size`, `-w`: Window size for sequence minimizers (default: 40)
- `--kmer_size`, `-k`: K-mer size for sequence minimizers (default: 8)

General Options:
- `--outprefix`, `-o`: The output file prefix (default: "hammock")
- `--threads`: Number of threads to use for parallel processing (default: number of CPUs)
- `--debug`: Output debug information including sparsity comparisons

### Output

Results are written to CSV files with the following format:
```bed1, bed2, sketch_type, mode, num_hashes, precision, subsample, balance, jaccardfunc, jaccardcalc, intersect, union```

## Examples

```bash
# Compare intervals using HyperLogLog (mode A)
# Automatically uses FastHyperLogLog for 2-5x speedup when installed with [fast]
hammock files.txt primary.txt --precision 12

# Compare points using MinHash (mode B)
hammock files.txt primary.txt --mode B --minhash --num_hashes 256

# Compare both with subsampling (automatically uses mode C)
hammock files.txt primary.txt --subA 0.5 --subB 0.5

# Compare sequences (automatically uses mode D)
hammock fasta_files.txt primary_fastas.txt

# Use 32-bit hashing
hammock files.txt primary.txt --hashsize 32

# Use 64-bit hashing (default)
hammock files.txt primary.txt --hashsize 64
```

### Performance Verification

If you installed with `pip install "hammock[fast]"`, you can verify acceleration is working:

```python
# Check if FastHyperLogLog acceleration is available
from hammock.lib import get_performance_info
info = get_performance_info()
print(f"Cython acceleration: {info['cython_available']}")
print(f"Expected speedup: {info['performance_gain']}")

# Create an optimized sketch
from hammock.lib import create_fast_hyperloglog
sketch = create_fast_hyperloglog(precision=12, debug=True)
# Output shows: "Using Cython acceleration for FastHyperLogLog"
```

## Visualization (Experimental)

**⚠️ Note: This is an experimental feature and is not fully supported. Use at your own discretion.**

Hammock includes an R script for generating heatmaps from comparison results, particularly useful for visualizing ENCODE ChIP-seq data comparisons. After installation, the script is available as a system command:

```bash
# Generate heatmaps with default "encode" prefix
encode_heatmap.R encode.report hammock.results

# Generate heatmaps with custom prefix
encode_heatmap.R encode.report hammock.results my_experiment
```

### Requirements
- R must be installed and available in your system PATH
- R packages: `tidyr`, `dplyr`, `tibble`, `gplots`, `RColorBrewer`, `readr`, `stringr`

### Input Files
- `encode.report`: Tab-delimited ENCODE metadata file with columns for accession, biosample, organism, and target
- `hammock.results`: CSV output from hammock comparison

### Output
Three PDF heatmaps showing Jaccard similarities with border colors indicating:
- `*_biosample.pdf`: Colored by biosample type
- `*_organism.pdf`: Colored by organism
- `*_target.pdf`: Colored by target of assay

### Customization
The R script is located at `scripts/encode_heatmap.R` in the hammock installation directory. Advanced users can modify this script to customize heatmap generation, color schemes, or add additional visualizations.

## Testing

In order to run the tests requiring BigBed files, you need to have `bedToBigBed` installed and in your PATH.
```bash
#conda install -c bioconda ucsc-bedtools
```

Run all tests:
```bash
python3 -m pytest hammock/tests            # Run all tests
python3 -m pytest hammock/tests -m quick   # Run quick tests only
python3 -m pytest hammock/tests -m full    # Run full tests only
```

### Test Files
- `test_hyperloglog.py`: Tests HyperLogLog sketching with various precisions and set sizes
- `test_minhash.py`: Tests MinHash sketching with different numbers of hash functions
- `test_intervals.py`: Tests interval processing and mode detection
- `test_sequences.py`: Tests sequence processing and mode D
- `test_hammock.py`: Tests CLI functionality and mode switching
- `benchmark_sketches.py`: Performance benchmarks comparing different sketching methods

### Test Cases
Tests include:
- Perfect overlap scenarios
- No overlap scenarios
- Partial overlap with different set sizes
- Sparse vs dense comparisons
- Edge cases with very different sparsity levels
- Automatic mode detection and switching
- Parameter validation

### Running Specific Tests

```bash
# Run HyperLogLog tests
python -m pytest hammock/tests/test_hyperloglog.py
# Run MinHash tests
python -m pytest hammock/tests/test_minhash.py
# Run benchmarks
python -m pytest hammock/tests/benchmark_sketches.py
```

## Fast HyperLogLog with Rust

For maximum performance, Hammock includes a native Rust implementation of HyperLogLog called `RustHLLWrapper`. This implementation provides significant speed improvements, especially for large datasets.

### Installation

The Rust implementation requires building the Rust extension:

```bash
cd rust_hll
python -m maturin develop
```

### Usage

```python
from hammock import RustHLLWrapper

# Create a sketch with precision 12 (2^12 = 4096 registers)
sketch = RustHLLWrapper(precision=12)

# Add values
sketch.add("item1")
sketch.add("item2")

# Add in batch (much faster)
sketch.add_batch(["item3", "item4", "item5", "item6"])

# Get cardinality estimate
estimate = sketch.cardinality()
print(f"Estimated cardinality: {estimate}")

# Merge with another sketch
sketch2 = RustHLLWrapper(precision=12)
sketch2.add_batch(["item5", "item6", "item7", "item8"])
sketch.merge(sketch2)

# Calculate Jaccard similarity
jaccard = sketch.jaccard(sketch2)
print(f"Jaccard similarity: {jaccard}")
```

### Performance Features

The Rust implementation includes several performance optimizations:

- Multi-threaded batch processing for large datasets
- SIMD-style vectorized operations
- Optimized hashing and register access
- Fast Maximum Likelihood Estimation

You can control threading with:

```python
# Disable threading
sketch = RustHLLWrapper(precision=12, use_threading=False)

# Enable threading with custom batch size threshold
sketch = RustHLLWrapper(precision=12, use_threading=True, min_thread_batch=100000)
```