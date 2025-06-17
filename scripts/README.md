# Hammock Scripts Directory

This directory contains various utility scripts for running hammock analyses, benchmarking, and processing results. The scripts are organized into several categories based on their functionality.

## Overview

After installation (`pip install .`), all scripts in this directory become globally accessible from any location on your system. Scripts that depend on other scripts in this directory will automatically find their dependencies regardless of where you run them from.

## Script Categories

### 🔬 Analysis and Comparison Scripts

#### `compare_sim_matrices.py`
**Purpose**: Compare similarity matrices from bedtools and hammock outputs, calculating various metrics.

**Usage**:
```bash
compare_sim_matrices.py <hammock_output> <bedtools_output> [--table] [--quiet]
```

**Features**:
- Automatically detects file format (bedtools vs hammock)
- Supports both BED (mode B) and FASTA (mode D) file analyses
- Calculates Frobenius norm, correlation, mean/max absolute error
- Outputs detailed comparison metrics or table format
- Handles different file extensions (.bed.gz, .fa.gz, etc.)

**Example**:
```bash
# Compare matrices and get detailed output
compare_sim_matrices.py hammock_results.csv bedtools_output.txt

# Get table format for parameter sweeps
compare_sim_matrices.py hammock_results.csv bedtools_output.txt --table
```

---

### 🧪 Parameter Sweep Scripts

#### `complete_parameter_sweep.sh`
**Purpose**: Complete parameter sweep for hammock that:
1. Takes a list of BED files and runs pairwise bedtools jaccard on ALL combinations
2. Uses the resulting jaccard matrix as reference for hammock parameter sweep  
3. Compares hammock output against bedtools reference for different parameters

**Usage**:
```bash
complete_parameter_sweep.sh -b <bed_file_list> [-f <fasta_file_list>] [-o <output_prefix>] [-c] [-v] [-q]
```

**Parameters**:
- `-b`: BED file list (required) - text file with one BED file path per line
  - Used for generating bedtools pairwise jaccard reference
- `-f`: FASTA file list (optional) - text file with one FASTA file path per line
  - If not provided, hammock will run on the BED files
- `-o`: Output prefix with path (default: complete_parameter_sweep)
  - Creates: `<prefix>_results/` directory, `<prefix>_bedtools_ref.tsv`, `<prefix>_results.tsv`
- `-c`: Clean up intermediate files after completion (bedtools reference always preserved)
- `-v`: Verbose output (shows progress of bedtools comparisons)
- `-q, --quick`: Quick mode with limited parameter combinations for testing

**Mode Selection**:
- **BED files only**: Runs bedtools jaccard pairwise, then tests precision parameter (mode B)
- **FASTA files + BED**: Uses BED files for reference, tests all parameters on FASTA (mode D)

**Parameter Ranges**:
- **Full mode (default)**:
  - BED files: precision = 16,18,19,20,21,22,23,25
  - FASTA files: klen = 15,20,25; window = 100,200,500; precision = 20,23,25
- **Quick mode (-q)**:
  - BED files: precision = 20,23
  - FASTA files: klen = 20; window = 200; precision = 20,23

**Features**:
- Automatically generates bedtools pairwise jaccard matrix from BED file list
- Reuses existing bedtools reference if found (avoids time-intensive regeneration)
- Always preserves bedtools reference file for future use (even with cleanup enabled)
- Skips bedtools generation if reference already exists

**Examples**:
```bash
# Run complete sweep on BED files only
complete_parameter_sweep.sh -b bed_files_list.txt -c -v

# Run sweep on FASTA files using corresponding BED files for reference
complete_parameter_sweep.sh -b bed_files_list.txt -f fasta_files_list.txt -c -v

# Quick mode for testing
complete_parameter_sweep.sh -b bed_files_list.txt --quick -c -v

# Custom output prefix and path
complete_parameter_sweep.sh -b bed_files_list.txt -o /path/to/my_experiment -c -v
```

#### `complete_parameter_sweep_quick.sh`
**Purpose**: Faster version of complete parameter sweep with limited parameter combinations for testing.

**Usage**: Same as `complete_parameter_sweep.sh`

**Parameter Differences**:
- **BED files**: Tests only precision (20, 23)
- **FASTA files**: Tests only k-mer=20, window=200, precision (20, 23)

#### `parameter_sweep.sh`
**Purpose**: Run hammock parameter sweep against pre-generated bedtools reference matrix.

**Usage**:
```bash
parameter_sweep.sh -b <bedtools_output_file> -f <file_list> [-o <output_dir>] [-r <results_table>] [-q]
```

**Parameters**:
- `-b`: Path to bedtools pairwise jaccard output matrix file (NOT a list of bed files - must be the actual bedtools output)
- `-f`: File containing list of BED or FASTA file paths (one per line)
- `-o`: Output directory for hammock results (default: parameter_sweep_results)
- `-r`: Results table filename (default: parameter_sweep_results.tsv)
- `-q, --quick`: Quick mode with limited parameter combinations for testing

**File Type Detection**:
- BED files (.bed): Only precision parameter is varied (mode B)
- FASTA files (.fa, .fasta, etc.): All parameters varied (mode D)

**Parameter Ranges**:
- **Full mode (default)**:
  - BED files: precision = 16,18,19,20,21,22,23,25
  - FASTA files: klen = 15,20,25; window = 100,200,500; precision = 20,23,25
- **Quick mode (-q)**:
  - BED files: precision = 20,23
  - FASTA files: klen = 20; window = 200; precision = 20,23

**Examples**:
```bash
# Basic parameter sweep
parameter_sweep.sh -b bedtools_reference.txt -f files.txt -o sweep_results/

# Quick mode for testing
parameter_sweep.sh -b bedtools_output.txt -f bed_files_list.txt --quick

# For FASTA files (assuming you have bedtools output)
parameter_sweep.sh -b bedtools_fasta_output.txt -f fasta_files_list.txt
```

#### `parameter_sweep_quick.sh`
**Purpose**: Quick version of parameter sweep with limited combinations.

**Usage**: Same as `parameter_sweep.sh` but with reduced parameter space for faster testing.

---

### 🔧 Utility Scripts

#### `bedtools_pairwise.sh`
**Purpose**: Generate pairwise Jaccard similarity matrix using bedtools for comparison reference.

**Usage**:
```bash
bedtools_pairwise.sh -f <file_list> [-o <output_file>] [-v]
```

**Features**:
- Handles both .bed and .bed.gz files
- Automatically sorts files for bedtools compatibility
- Generates reference similarity matrix for hammock comparison

**Example**:
```bash
bedtools_pairwise.sh -f bed_files.txt -o reference_matrix.txt -v
```

#### `ENCODE_report_to_key.sh`
**Purpose**: Process ENCODE metadata reports to extract file information and create lookup keys.

**Usage**:
```bash
ENCODE_report_to_key.sh <encode_report.tsv>
```

---

### 📊 Visualization Scripts

#### `encode_heatmap.R`
**Purpose**: Generate heatmaps from similarity matrices with ENCODE metadata integration.

**Usage**:
```R
# Run in R environment
source("encode_heatmap.R")
```

**Features**:
- Creates publication-quality heatmaps
- Integrates ENCODE metadata for labeling
- Supports various similarity matrix formats
- Customizable color schemes and annotations

---

### 🖥️ SLURM Job Scripts -- DRAFT, NOT INCLUDED

These scripts are designed for high-performance computing environments using SLURM job scheduler.

#### `run_hammock.slurm`
**Purpose**: Basic SLURM job script for single hammock runs.

**Features**:
- Configurable memory and CPU resources
- Error handling and logging
- Suitable for single dataset analysis

#### `run_hammock_batch.slurm`
**Purpose**: Batch processing of multiple hammock jobs.

**Features**:
- Array job support for processing multiple file lists
- Parallel execution across multiple nodes
- Automated output organization

#### `run_hammock_sweep.slurm`
**Purpose**: SLURM job for running parameter sweeps on HPC systems.

**Features**:
- High-memory allocation for large parameter spaces
- Integrated with parameter sweep scripts
- Automatic result aggregation

#### `run_mode_b_benchmark.slurm`
**Purpose**: Benchmarking script specifically for hammock mode B (BED files).

**Features**:
- Performance profiling
- Memory usage tracking
- Comparative analysis with bedtools

---

## Workflow Examples

### Basic Parameter Sweep Workflow
```bash
# 1. Quick test with limited parameters (automatically generates bedtools reference)
complete_parameter_sweep.sh -b test_files.txt -o test_results -v -q

# 2. Full parameter sweep (automatically generates bedtools reference)
complete_parameter_sweep.sh -b all_files.txt -o full_results -c -v

# 3. Rerun sweep with same data (reuses existing bedtools reference - much faster!)
complete_parameter_sweep.sh -b all_files.txt -o full_results -c -v

# 4. Custom parameter sweep with pre-existing bedtools reference
parameter_sweep.sh -b full_results_bedtools_ref.tsv -f files.txt -o custom_sweep/
```

### Manual Analysis Workflow
```bash
# 1. Generate bedtools pairwise jaccard reference matrix
bedtools_pairwise.sh -f files.txt -o bedtools_ref.txt

# 2. Run hammock with specific parameters
hammock files.txt files.txt -o hammock_results -p 20 --mode B

# 3. Compare results (using bedtools output matrix)
compare_sim_matrices.py hammock_results_hll_p20_jaccB.csv bedtools_ref.txt
```

### HPC Workflow
```bash
# Submit parameter sweep job to SLURM
sbatch run_hammock_sweep.slurm

# Submit batch processing job
sbatch run_hammock_batch.slurm
```

## File Type Auto-Detection

Scripts automatically detect file types based on extensions:

- **BED files**: `.bed`, `.bed.gz` → Uses hammock mode B
- **FASTA files**: `.fa`, `.fasta`, `.fna`, `.ffn`, `.faa`, `.frn` (with optional `.gz`) → Uses hammock mode D

## Dependencies

### Required Software
- **hammock**: Main analysis tool
- **bedtools**: For reference similarity calculations
- **Python 3.6+**: With numpy, matplotlib, scipy, pandas
- **R**: For visualization scripts (ggplot2, pheatmap, etc.)
- **bash**: Shell scripting environment

### Script Dependencies
Scripts that call other scripts will automatically find them after installation:
- `complete_parameter_sweep*.sh` → `parameter_sweep*.sh`
- `parameter_sweep*.sh` → `compare_sim_matrices.py`

## Output Files

### Parameter Sweep Results
Results are saved in TSV format with columns:
- `klen`, `window`, `precision`: Parameter values
- `hammock_file`, `bedtools_file`: Input files
- `format1`, `format2`: File format detection
- `matrix_size`: Number of files compared
- `frobenius_norm`: Matrix difference metric
- `correlation`: Pearson correlation coefficient
- `mean_abs_error`, `max_abs_error`: Error metrics
- `runtime_seconds`: Execution time

### Log Files
- Most scripts generate detailed logs for debugging
- SLURM scripts create separate `.out` and `.err` files
- Use `-v` flag for verbose output where available

## Troubleshooting

### Common Issues
1. **"Script not found"**: Run `pip install .` from the main hammock directory
2. **"Permission denied"**: Ensure scripts have execute permissions
3. **"Bedtools failed"**: Check that input BED files are properly formatted
4. **"Comparison failed"**: Verify that both hammock and bedtools outputs exist

### Getting Help
```bash
# Most scripts support help flags
script_name.sh -h
script_name.py --help
```

---

*This README was generated for hammock v0.2.1. For the latest documentation, visit the [hammock GitHub repository](https://github.com/jessicabonnie/hammock).* 