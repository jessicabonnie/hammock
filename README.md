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

For **2-15x performance improvement** with FastHyperLogLog and C++ acceleration:

```bash
# Install with Cython acceleration support
git clone https://github.com/jessicabonnie/hammock.git
cd hammock
pip install -e ".[fast]"
```

This will:
- ✅ Enable FastHyperLogLog with C++ acceleration (2-15x speedup)
- ✅ Automatically fall back to Cython (2-5x speedup) or Python if compilation fails
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
| **FastHyperLogLog (C++)** | `pip install "hammock[fast]"` | **2-15x faster** | 🚀 **🚀 Recommended for most users**<br>✅ Maximum speedup<br>✅ Automatic fallback<br>✅ Easy installation |
| **FastHyperLogLog (Cython)** | `pip install "hammock[fast]"` | **2-5x faster** | ✅ Good speedup<br>✅ Automatic fallback<br>✅ Easy installation |

**Performance results** (based on comprehensive testing):
- **Small datasets** (100-1000 items): FastHyperLogLog (C++) provides ~2.5x speedup
- **Medium datasets** (1000-10000 items): FastHyperLogLog (C++) provides ~3.3x speedup  
- **Large datasets** (10000+ items): FastHyperLogLog (C++) provides ~2-15x speedup
- **Batch operations**: Up to 15x speedup with optimized batch processing
- **C++ vs Python**: Consistent 2-15x improvement across all dataset sizes
- **Cython fallback**: 2-5x speedup when C++ compilation fails

**Recommendation**: Install with `pip install "hammock[fast]"` for the best balance of performance and ease of use. FastHyperLogLog automatically uses C++ acceleration when available, providing maximum performance with automatic fallbacks.


## Usage

Calculate pairwise Jaccard similarities between lists of BED or sequence files.
Both input files are text files containing a path per line to the files to be compared. For computational purposes, `<primary_paths_file>` should be the shorter of the two file lists.
```bash
hammock <filepaths_file> <primary_file> --mode {A|B|C|D} [options]
```

## Analysis and Visualization Scripts

See `scripts/README.md` for full details. Notable utilities:

- `scripts/clustering_analysis.py`: Evaluate hierarchical clustering (NMI, silhouette) over multiple cluster counts and linkage methods from a hammock CSV or bedtools TSV. Supports long-table TSV to stdout or file.
- `scripts/draw_dendrograms.py` (new): Create a multi-page PDF with one dendrogram per linkage method, colored by tissue labels, from a hammock CSV or bedtools TSV plus an accession key.

Examples:

```bash
# Long-table clustering analysis
python3 scripts/clustering_analysis.py results/mouseonly_hll_p18_jaccA.csv data/filtered_accession_key.tsv --stdout

# Multi-page dendrogram PDF
python3 scripts/draw_dendrograms.py results/mouseonly_hll_p18_jaccA.csv data/filtered_accession_key.tsv --output results/mouseonly_dendrograms.pdf
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


**Performance Note**: When you install with `pip install "hammock[fast]"`, HyperLogLog sketching automatically uses FastHyperLogLog with C++ acceleration for 2-15x better performance. If C++ compilation fails, it automatically falls back to Cython acceleration (2-5x speedup) or pure Python. This works transparently with all existing commands.

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
# Automatically uses FastHyperLogLog for 2-15x speedup when installed with [fast]
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
- `benchmark_hll.py`: Comprehensive HyperLogLog benchmarks including C++ vs Python performance

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

# Run comprehensive HyperLogLog benchmarks (includes C++ vs Python comparison)
cd benchmarks
python benchmark_hll.py
```



# Fast HyperLogLog with C++ (Recommended)

For **2-15x performance improvement**, Hammock includes a high-performance C++ implementation of HyperLogLog called `FastHyperLogLog`. This implementation provides the best balance of performance and ease of use, automatically falling back to Cython or Python if C++ compilation fails.

## Obtaining the C++ HLL Code

**The C++ HyperLogLog implementation is already included** when you clone the main hammock repository. You don't need to download or clone any additional repositories. Simply clone hammock and the C++ code will be available in the `hll/` directory:

```bash
git clone https://github.com/jessicabonnie/hammock.git
cd hammock
# The C++ code is now available in hll/ directory
ls hll/
```

**Note**: The C++ HyperLogLog implementation is based on the [hll library](https://github.com/mindis/hll) by Daniel Baker, which provides a high-performance C++ implementation with SIMD parallelism. This library is included under the MIT License and is fully open source for both academic and commercial use.

This implementation provides:

- **SIMD optimizations** (SSE2, AVX2, AVX-512) for maximum performance
- **Thread-safe operations** for concurrent access
- **Optimized memory management** for large datasets
- **High-precision cardinality estimation** with multiple correction methods

## Installation and Setup

### Prerequisites

Ensure you have the following build tools installed:

```bash
# Ubuntu/Debian
sudo apt-get install build-essential python3-dev

# CentOS/RHEL/Rocky Linux
sudo yum groupinstall "Development Tools"
sudo yum install python3-devel

# macOS (using Homebrew)
brew install gcc python3

# Windows (using MSVC)
# Install Visual Studio Build Tools or use WSL
```

### Automatic Installation (Recommended)

The C++ implementation is automatically built when you install hammock with the `[fast]` option. The C++ source code is already included in the repository:

```bash
git clone https://github.com/jessicabonnie/hammock.git
cd hammock
pip install -e ".[fast]"
```

This will:
- ✅ Compile the C++ HyperLogLog library
- ✅ Build the Cython wrapper
- ✅ Install FastHyperLogLog with automatic acceleration detection
- ✅ Provide 2-15x performance improvement
- ✅ Automatically fall back to Cython or Python if compilation fails

### Manual Build (Advanced Users)

If you need to customize the build process or troubleshoot compilation issues:

```bash
# 1. Build the standalone C++ library
cd hll
make clean && make

# 2. Build the Python extension
cd ..
python setup.py build_ext --inplace

# 3. Test the installation
python -c "from hammock.lib.hyperloglog_fast import FastHyperLogLog; print('C++ HLL successfully installed!')"
```

### Build Script

Use the provided build script for automated compilation:

```bash
python build_cpp_extension.py
```

This script will:
- Check dependencies
- Build the C++ library
- Compile the Cython extension
- Run basic tests
- Provide detailed error messages if something fails

## Usage

FastHyperLogLog automatically detects and uses the best available acceleration:

```python
from hammock.lib.hyperloglog_fast import FastHyperLogLog

# Create an optimized sketch (automatically uses C++ if available)
sketch = FastHyperLogLog(precision=12, hash_size=64)

# Add items individually
sketch.add_string("item1")
sketch.add_string("item2")

# Add items in batch (much faster)
sketch.add_batch(["item3", "item4", "item5", "item6"])

# Get cardinality estimate
estimate = sketch.estimate_cardinality()
print(f"Estimated cardinality: {estimate}")

# Check which acceleration is being used
print(f"Acceleration type: {sketch._acceleration_type}")
print(f"Performance info: {sketch.get_performance_info()}")
```

## Performance Features

The C++ implementation includes several performance optimizations:

- **SIMD vectorization** for register operations
- **Optimized hashing** using xxHash
- **Memory-efficient storage** with bit-packing
- **Fast set operations** (union, intersection)
- **Automatic precision optimization** based on data size

## Performance Comparison

| Implementation | Speedup | Installation | Best For |
|----------------|---------|--------------|----------|
| **Pure Python** | 1x | `pip install hammock` | Simple setup, small datasets |
| **FastHyperLogLog (C++)** | **2-15x** | `pip install "hammock[fast]"` | **🚀 Best balance** |
| **FastHyperLogLog (Cython)** | 2-5x | `pip install "hammock[fast]"` | Fallback when C++ fails |

## Troubleshooting

### Common Issues

1. **"No module named 'hammock.lib.cpp_hll_wrapper'"**
   - Ensure you have build tools installed
   - Try: `python build_cpp_extension.py`
   - Check: `python setup.py build_ext --inplace`

2. **"C++ compiler not found"**
   - Install build-essential (Ubuntu) or Development Tools (CentOS)
   - On macOS: `xcode-select --install`

3. **"Cython compilation failed"**
   - Update Cython: `pip install --upgrade cython`
   - Check Python version compatibility

4. **Performance not as expected**
   - Verify C++ acceleration: `sketch._acceleration_type == 'C++'`
   - Check precision settings (higher precision = better accuracy, slower speed)

### Verification

Test that C++ acceleration is working:

```python
from hammock.lib.hyperloglog_fast import FastHyperLogLog, CPP_AVAILABLE

print(f"C++ available: {CPP_AVAILABLE}")

sketch = FastHyperLogLog(precision=12)
print(f"Acceleration type: {sketch._acceleration_type}")

# Should show "C++" if working correctly
```

## Benchmarks

Run comprehensive performance benchmarks:

```bash
cd benchmarks
python benchmark_hll.py
```

This will test all implementations and show detailed performance comparisons including:
- Individual add performance
- Batch add performance  
- Cardinality estimation speed
- Memory usage
- Accuracy comparisons

## Citation and Attribution

The C++ HyperLogLog implementation in hammock is based on the following work:

**Baker, D. (2016). hll: C++ HyperLogLog Implementation with SIMD Parallelism.**  
GitHub repository: [https://github.com/mindis/hll](https://github.com/mindis/hll)  
License: MIT License

**If you use the C++ HyperLogLog acceleration in hammock, please cite:**

Baker, D. (2016). hll: C++ HyperLogLog Implementation with SIMD Parallelism. GitHub repository. https://github.com/mindis/hll

**Original HyperLogLog Algorithm:**

Flajolet, P., Fusy, É., Gandouet, O., & Meunier, F. (2007). HyperLogLog: the analysis of a near-optimal cardinality estimation algorithm. In *Discrete Mathematics and Theoretical Computer Science* (pp. 137-156).