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

#### `compare_clustering_trees.py`
**Purpose**: Compare hierarchical clustering trees built from Jaccard similarity matrices using Robinson-Foulds distance. This provides a biologically meaningful comparison of how well hammock approximates bedtools clustering patterns.

**Usage**:
```bash
compare_clustering_trees.py <file1> <file2> [--linkage METHOD] [--table] [--verbose] [--save-trees FILE] [--output FILE]
```

**Features**:
- Automatically detects file format (bedtools vs hammock)
- Builds hierarchical clustering trees from similarity matrices
- Calculates accurate Robinson-Foulds distance using ete3 library
- Supports multiple linkage methods (single, complete, average, ward)
- Outputs normalized RF distance (0-1 scale) and tree similarity
- Can save Newick format trees for visualization
- Table format output for integration with parameter sweeps

**Parameters**:
- `--linkage, -l`: Clustering linkage method (default: average)
  - `single`: Single linkage (minimum distance)
  - `complete`: Complete linkage (maximum distance) 
  - `average`: Average linkage (UPGMA)
  - `ward`: Ward's minimum variance method
- `--table, -t`: Output results as tab-delimited line
- `--header`: Print tab-delimited header (use with --table)
- `--verbose, -v`: Show detailed progress and tree information
- `--save-trees FILE`: Save Newick format trees to file
- `--output, -o FILE`: Save detailed comparison results to file

**Output Metrics**:
- **Robinson-Foulds distance**: Number of different splits between trees
- **Maximum RF distance**: Theoretical maximum RF distance for n leaves
- **Normalized RF distance**: RF distance / max RF distance (0-1 scale)
- **Tree similarity**: 1 - normalized RF distance (higher = more similar)

**Dependencies**:
- **Required**: scipy, pandas, numpy
- **Recommended**: ete3 (for accurate RF calculations)
- **Alternative**: dendropy (fallback if ete3 unavailable)

**Example**:
```bash
# Basic tree comparison with verbose output
compare_clustering_trees.py hammock_results.csv bedtools_output.tsv --verbose

# Compare with different linkage methods
compare_clustering_trees.py hammock_results.csv bedtools_output.tsv --linkage complete

# Table format for parameter sweeps
compare_clustering_trees.py hammock_results.csv bedtools_output.tsv --table

# Save trees and detailed results for further analysis
compare_clustering_trees.py hammock_results.csv bedtools_output.tsv \
    --save-trees clustering_trees.newick \
    --output detailed_comparison.txt \
    --verbose

# Print header for parameter sweep integration
compare_clustering_trees.py --header
```

**Interpretation**:
- **RF distance = 0**: Trees are identical (perfect match)
- **RF distance = max**: Trees are maximally different
- **Normalized RF < 0.2**: Trees are very similar (good approximation)
- **Normalized RF > 0.8**: Trees are very different (poor approximation)

**Integration with Parameter Sweeps**:
The table output format makes it easy to integrate with parameter sweep workflows:
```bash
# Header output
file1	file2	format1	format2	matrix_size	rf_distance	max_rf_distance	normalized_rf	linkage_method

# Example result
hammock_p20_k15_w100.csv	bedtools_ref.tsv	hammock	bedtools	68	114	130	0.876923	average
```

#### `draw_clustering_tree.py`
**Purpose**: Generate hierarchical clustering dendrograms from bedtools/hammock similarity matrices OR visualize saved Newick format trees with automatic format detection.

**Usage**:
```bash
# From similarity matrices
draw_clustering_tree.py <matrix_file> [--output OUTPUT] [--linkage METHOD] [--figsize W H] [options]

# From Newick tree files
draw_clustering_tree.py <newick_file> [--output OUTPUT] [--tree-index N] [--list-trees] [options]
```

**Features**:
- **NEW**: Automatically detects file format (bedtools, hammock, or Newick)
- **NEW**: Visualizes saved Newick trees from `compare_clustering_trees.py --save-trees`
- **NEW**: Supports multi-tree Newick files with tree selection by index
- Supports all major linkage methods (single, complete, average, ward) for similarity matrices
- Generates high-quality dendrogram visualizations
- Customizable figure size, labels, and color schemes
- Detailed clustering summary with similarity statistics (for similarity matrices)
- Saves plots in multiple formats (PNG, PDF, SVG, etc.)
- Option to suppress plot display for batch processing

**Parameters**:
- `input_file`: Input file (similarity matrix: bedtools/hammock format, or Newick tree file)
- `--output, -o FILE`: Output file for dendrogram (PNG, PDF, SVG, etc.)
- `--linkage, -l METHOD`: Linkage method (single, complete, average, ward; default: average, ignored for Newick files)
- `--figsize W H`: Figure size as width height (default: 12 8)
- `--title, -t TITLE`: Custom title for the dendrogram
- `--max-label-length N`: Maximum length for sample labels (default: 20, 0 for no truncation)
- `--color-threshold FLOAT`: Distance threshold for coloring clusters (ignored for Newick files)
- `--no-show`: Do not display the plot (useful when saving to file)
- `--summary, -s`: Print detailed clustering summary (ignored for Newick files)
- `--quiet, -q`: Suppress progress messages
- **NEW** `--tree-index N`: Index of tree to visualize from Newick file (default: 0, first tree)
- **NEW** `--list-trees`: List all trees in Newick file and exit

**Examples**:
```bash
# Basic dendrogram generation with display
draw_clustering_tree.py similarity_matrix.csv

# Save dendrogram to file without displaying
draw_clustering_tree.py bedtools_output.tsv --output tree.png --no-show

# Generate with different linkage method and custom size
draw_clustering_tree.py hammock_results.csv --linkage complete --figsize 16 10 --output large_tree.pdf

# Batch processing with summary statistics
draw_clustering_tree.py data.csv --no-show --output tree.png --summary --quiet

# Compare different linkage methods
for method in single complete average ward; do
    draw_clustering_tree.py data.csv --linkage $method --output tree_${method}.png --no-show
done
```

**Output Information**:
- **Dendrogram Plot**: Visual representation of hierarchical clustering
- **Clustering Summary** (with --summary): 
  - Number of samples and matrix dimensions
  - Similarity range, mean, and median
  - Clustering distance statistics
  - Top 10 most similar sample pairs
- **Distance Scale**: Y-axis shows distance (1 - Jaccard similarity)
- **Color Coding**: Clusters colored by distance threshold (customizable)

---

### 🧪 Parameter Sweep Scripts

#### `complete_parameter_sweep.sh`
**Purpose**: Complete parameter sweep for hammock that:
1. Takes a list of BED files and runs pairwise bedtools jaccard on ALL combinations
2. Uses the resulting jaccard matrix as reference for hammock parameter sweep  
3. Compares hammock output against bedtools reference for different parameters
4. Automatically includes clustering tree analysis when `compare_clustering_trees.py` is available

**Usage**:
```bash
complete_parameter_sweep.sh -b <bed_file_list> [-f <fasta_file_list>] [-o <output_prefix>] [-c] [-v] [-q] [-d] [--klen <values>] [--window <values>] [--precision <values>]
```

**Parameters**:
- `-b`: BED file list (required) - text file with one BED file path per line
  - Used for generating bedtools pairwise jaccard reference
- `-f`: FASTA file list (optional) - text file with one FASTA file path per line
  - If not provided, hammock will run on the BED files
- `-o`: Output prefix with path (default: complete_parameter_sweep)
  - Creates: `<prefix>_results/` directory, `<prefix>_bedtools_ref.tsv`, `<prefix>_results.tsv`, `<prefix>_results_clustering.tsv`
- `-c`: Clean up intermediate files after completion (bedtools reference always preserved)
- `-v`: Verbose output (shows progress of bedtools comparisons)
- `-q, --quick`: Quick mode with limited parameter combinations for testing
- `-d, --debug`: Debug mode with detailed output and preserved intermediate files
- `--klen <values>`: **NEW** Custom klen values (comma-separated, e.g., '10,15,20')
- `--window <values>`: **NEW** Custom window values (comma-separated, e.g., '50,100,200')
- `--precision <values>`: **NEW** Custom precision values (comma-separated, e.g., '20,22,24')

**Mode Selection**:
- **BED files only**: Runs bedtools jaccard pairwise, then tests precision parameter (mode B)
- **FASTA files + BED**: Uses BED files for reference, tests all parameters on FASTA (mode D)

**Parameter Ranges**:
- **Custom parameters (--klen, --window, --precision)**: **NEW**
  - Uses your specified values instead of defaults/quick mode
  - Custom parameters override default/quick mode parameters
  - Values must be comma-separated without spaces (e.g., '10,15,20,25')
- **Full mode (default)**:
  - BED files: precision = 16,19,20,21,22,23,24,26,28,30
  - FASTA files: klen = 10,15,20,25,30; window = 25,50,100,200,300,400,500; precision = 16,19,20,21,22,23,24,26,28,30
- **Quick mode (-q)**:
  - BED files: precision = 20,23
  - FASTA files: klen = 20; window = 200; precision = 20,23

**Features**:
- Automatically generates bedtools pairwise jaccard matrix from BED file list
- Reuses existing bedtools reference if found (avoids time-intensive regeneration)
- Always preserves bedtools reference file for future use (even with cleanup enabled)
- Skips bedtools generation if reference already exists
- **NEW**: Automatically runs clustering tree analysis when available
- **NEW**: Generates dual output tables: similarity matrices and clustering comparisons
- **NEW**: Shows best parameter combinations based on lowest normalized Robinson-Foulds distance
- **NEW**: Automatically generates bedtools reference dendrogram visualization
- **NEW**: Debug mode provides detailed output and preserves intermediate files for troubleshooting

**Output Files**:
- `<prefix>_results.tsv`: Traditional similarity matrix comparison results
- `<prefix>_results_clustering.tsv`: Clustering tree topology comparison results (when available)
- `<prefix>_bedtools_ref.tsv`: Bedtools pairwise jaccard reference matrix
- `<prefix>_bedtools_reference_dendrogram.png`: Visual dendrogram of bedtools reference clustering (when available)

**Examples**:
```bash
# Run complete sweep on BED files only (includes clustering analysis)
complete_parameter_sweep.sh -b bed_files_list.txt -c -v

# Run sweep on FASTA files using corresponding BED files for reference
complete_parameter_sweep.sh -b bed_files_list.txt -f fasta_files_list.txt -c -v

# Quick mode for testing (still includes clustering if available)
complete_parameter_sweep.sh -b bed_files_list.txt --quick -c -v

# Debug mode with detailed output and preserved intermediate files
complete_parameter_sweep.sh -b bed_files_list.txt -f fasta_files_list.txt -d -v

# Custom output prefix and path
complete_parameter_sweep.sh -b bed_files_list.txt -o /path/to/my_experiment -c -v

# Custom parameter values for FASTA files (NEW)
complete_parameter_sweep.sh -b bed_files_list.txt -f fasta_files_list.txt --klen 15,20,25 --window 100,200 --precision 20,22,24 -c -v

# Custom precision values for BED files (NEW)
complete_parameter_sweep.sh -b bed_files_list.txt --precision 18,20,22,24 -c -v

# Use default klen and window, but custom precision for FASTA files (NEW)
complete_parameter_sweep.sh -b bed_files_list.txt -f fasta_files_list.txt --precision 19,21,23 -c -v

# View clustering results (best parameters have lowest normalized RF distance)
sort -k8 -n my_experiment_results_clustering.tsv | head -5
```

#### `complete_parameter_sweep_quick.sh`
**Purpose**: Faster version of complete parameter sweep with limited parameter combinations for testing.

**Usage**: Same as `complete_parameter_sweep.sh`

**Parameter Differences**:
- **BED files**: Tests only precision (20, 23)
- **FASTA files**: Tests only k-mer=20, window=200, precision (20, 23)

#### `parameter_sweep.sh`
**Purpose**: Run hammock parameter sweep against pre-generated bedtools reference matrix with optional clustering tree analysis.

**Usage**:
```bash
parameter_sweep.sh -b <bedtools_output_file> -f <file_list> [-o <output_dir>] [-r <results_table>] [-q] [-d] [--klen <values>] [--window <values>] [--precision <values>]
```

**Parameters**:
- `-b`: Path to bedtools pairwise jaccard output matrix file (NOT a list of bed files - must be the actual bedtools output)
- `-f`: File containing list of BED or FASTA file paths (one per line)
- `-o`: Output directory for hammock results (default: parameter_sweep_results)
- `-r`: Results table filename (default: parameter_sweep_results.tsv)
- `-q, --quick`: Quick mode with limited parameter combinations for testing
- `-d, --debug`: Debug mode with detailed output and preserved intermediate files
- `--klen <values>`: **NEW** Custom klen values (comma-separated, e.g., '10,15,20')
- `--window <values>`: **NEW** Custom window values (comma-separated, e.g., '50,100,200')
- `--precision <values>`: **NEW** Custom precision values (comma-separated, e.g., '20,22,24')

**File Type Detection**:
- BED files (.bed): Only precision parameter is varied (mode B)
- FASTA files (.fa, .fasta, etc.): All parameters varied (mode D)

**Features**:
- **NEW**: Automatically detects and uses `compare_clustering_trees.py` when available
- **NEW**: Generates dual output tables: similarity matrices and clustering tree comparisons
- **NEW**: Uses average linkage clustering method for Robinson-Foulds distance calculations
- **NEW**: Automatically generates bedtools reference dendrogram visualization
- **NEW**: Debug mode provides detailed output and preserves intermediate files for troubleshooting
- Enhanced progress reporting with clustering analysis status
- Graceful fallback when clustering comparison tools are unavailable

**Debug Mode Features**:
When enabled with `-d` or `--debug`, the script provides:
- **Detailed clustering output**: Shows clustering exit codes, result lengths, and content
- **File operation tracking**: Shows table paths, existence checks, and write operations
- **Table monitoring**: Reports table sizes before and after each write operation
- **Content visualization**: Displays the actual lines being written to tables
- **File preservation**: Preserves intermediate hammock output files for manual inspection
- **Enhanced error reporting**: Shows complete error messages and debugging information

**When to Use Debug Mode**:
- Investigating empty or incomplete results tables
- Troubleshooting clustering comparison failures
- Verifying parameter combinations are being processed correctly
- Manually inspecting intermediate hammock output files
- Debugging file I/O issues or table generation problems

**Output Files**:
- `<results_table>`: Traditional similarity matrix comparison results
- `<results_table_base>_clustering.tsv`: Clustering tree topology comparison results (when available)
- `bedtools_reference_dendrogram.png`: Visual dendrogram of bedtools reference clustering (when available)

**Parameter Ranges**:
- **Custom parameters (--klen, --window, --precision)**: **NEW**
  - Uses your specified values instead of defaults/quick mode
  - Custom parameters override default/quick mode parameters
  - Values must be comma-separated without spaces (e.g., '10,15,20,25')
- **Full mode (default)**:
  - BED files: precision = 16,19,20,21,22,23,24,26,28,30
  - FASTA files: klen = 10,15,20,25,30; window = 25,50,100,200,300,400,500; precision = 16,19,20,21,22,23,24,26,28,30
- **Quick mode (-q)**:
  - BED files: precision = 20,23
  - FASTA files: klen = 20; window = 200; precision = 20,23

**Examples**:
```bash
# Basic parameter sweep (includes clustering analysis when available)
parameter_sweep.sh -b bedtools_reference.txt -f files.txt -o sweep_results/

# Quick mode for testing
parameter_sweep.sh -b bedtools_output.txt -f bed_files_list.txt --quick

# Debug mode with detailed output and preserved intermediate files
parameter_sweep.sh -b bedtools_reference.txt -f files.txt -o sweep_results/ --debug

# For FASTA files (assuming you have bedtools output)
parameter_sweep.sh -b bedtools_fasta_output.txt -f fasta_files_list.txt

# Custom parameter values for FASTA files (NEW)
parameter_sweep.sh -b bedtools_ref.tsv -f fasta_files_list.txt --klen 12,18,24 --window 75,150,300 --precision 18,20,22

# Custom precision values for BED files (NEW)
parameter_sweep.sh -b bedtools_ref.tsv -f bed_files_list.txt --precision 19,21,25

# Test only specific combinations (NEW)
parameter_sweep.sh -b bedtools_ref.tsv -f fasta_files_list.txt --klen 20 --window 100,200 --precision 22,24

# View clustering results to find best parameters
sort -k8 -n sweep_results/parameter_sweep_results_clustering.tsv | head -5
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
```bash
encode_heatmap.R <encode.report> <hammock.results> [file_prefix]
```

**Features**:
- Creates publication-quality heatmaps
- Integrates ENCODE metadata for labeling
- Supports various similarity matrix formats
- Customizable color schemes and annotations
- Generates three types of heatmaps:
  - Biosample-based clustering
  - Organism-based clustering
  - Target of assay-based clustering

**Output Files**:
- `{file_prefix}_biosample.pdf` - Heatmap colored by biosample
- `{file_prefix}_organism.pdf` - Heatmap colored by organism
- `{file_prefix}_target.pdf` - Heatmap colored by target of assay

#### `encode_PCA.R`
**Purpose**: Generate PCA plots from similarity matrices with ENCODE metadata integration, including both global and organism-specific analyses.

**Usage**:
```bash
encode_PCA.R <encode.report> <hammock.results> [file_prefix]
```

**Features**:
- **Global PCA**: Single PCA analysis of all organisms together
- **Organism-specific PCA**: Separate PCA analysis for each organism
- Large, readable text and point sizes for publication quality
- Color coding by biosample, shape coding by organism
- Variance explained percentages in axis labels
- Summary statistics for each organism

**Output Files**:
- `{file_prefix}_pca.pdf` - Global PCA plot (all organisms together)
- `{file_prefix}_organism_specific_pca.pdf` - Organism-specific PCA plots (one per organism)

**Analysis Types**:
1. **Global PCA**: Uses the full similarity matrix to show how all samples cluster together
2. **Organism-specific PCA**: Calculates separate PCA for each organism using only within-organism similarities

#### `encode_graphs.sh`
**Purpose**: Orchestrate the generation of all ENCODE visualization graphs (heatmaps and PCA plots) in one command.

**Usage**:
```bash
encode_graphs.sh <encode.report> <hammock.results> [file_prefix]
```

**Features**:
- Calls both `encode_heatmap.R` and `encode_PCA.R` automatically
- Validates input files before processing
- Clear progress reporting for each step
- Error handling with informative messages
- Summary of all generated files

**Generated Files**:
- `{file_prefix}_biosample.pdf` - Biosample heatmap
- `{file_prefix}_organism.pdf` - Organism heatmap
- `{file_prefix}_target.pdf` - Target heatmap
- `{file_prefix}_pca.pdf` - Global PCA plot (all organisms)
- `{file_prefix}_organism_specific_pca.pdf` - Organism-specific PCA plots

**Example Workflow**:
```bash
# Generate all ENCODE visualizations in one command
encode_graphs.sh encode_metadata.tsv hammock_results.csv my_experiment

# This produces 5 PDF files with comprehensive visualizations
```

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
# 1. Quick test with limited parameters (automatically generates bedtools reference + dendrogram)
complete_parameter_sweep.sh -b test_files.txt -o test_results -v -q

# 2. Full parameter sweep (automatically generates bedtools reference + dendrogram)
complete_parameter_sweep.sh -b all_files.txt -o full_results -c -v

# 3. Rerun sweep with same data (reuses existing bedtools reference - much faster!)
complete_parameter_sweep.sh -b all_files.txt -o full_results -c -v

# 4. Custom parameter sweep with pre-existing bedtools reference (generates dendrogram)
parameter_sweep.sh -b full_results_bedtools_ref.tsv -f files.txt -o custom_sweep/

# 5. Custom parameter sweep with specific values (NEW)
complete_parameter_sweep.sh -b files.txt -f fasta_files.txt --klen 15,20,25 --precision 20,22,24 -o custom_params/ -c -v

# 6. Narrow parameter sweep based on initial results (NEW)
parameter_sweep.sh -b custom_params_bedtools_ref.tsv -f fasta_files.txt --klen 20 --window 150,200,250 --precision 22 -o narrow_sweep/

# 7. View generated dendrogram to understand reference clustering structure
# full_results_bedtools_reference_dendrogram.png (or custom_sweep/bedtools_reference_dendrogram.png)

# 8. Analyze results - both similarity matrices and clustering trees
# Best numerical accuracy (lowest error)
sort -k7 -n full_results_results.tsv | head -5

# Best biological accuracy (lowest normalized RF distance)
sort -k8 -n full_results_results_clustering.tsv | head -5
```

### Debug Mode Workflow
```bash
# 1. Debug mode for investigating empty or incomplete results
complete_parameter_sweep.sh -b files.txt -o debug_test -d -v

# 2. Debug mode for parameter sweep with existing bedtools reference
parameter_sweep.sh -b bedtools_ref.tsv -f files.txt -o debug_sweep/ --debug

# 3. Manually inspect preserved intermediate files
ls debug_sweep/  # Shows all preserved hammock output files
head -5 debug_sweep/hammock_*_k15_w100_p20.csv  # Check hammock output format

# 4. Compare intermediate files with debug output
# Debug output shows exact table paths, write operations, and content
# Compare with actual table contents to identify issues

# 5. Debug clustering comparison issues
parameter_sweep.sh -b bedtools_ref.tsv -f files.txt -o debug_sweep/ --debug | grep "DEBUG:"
# Shows clustering exit codes, result lengths, and table operations

# 6. Use debug mode to troubleshoot specific parameter combinations
# Debug output shows which combinations succeed and which fail
```

### Custom Parameter Optimization Workflow *(NEW)*
```bash
# 1. Quick exploration with broad custom parameter ranges
complete_parameter_sweep.sh -b files.txt -f fasta_files.txt --klen 10,20,30 --window 50,200,400 --precision 18,22,26 -o broad_sweep/ -c -v

# 2. Analyze initial results to identify promising regions
sort -k11 -n broad_sweep_results_clustering.tsv | head -5  # Best biological accuracy
sort -k10 -n broad_sweep_results.tsv | head -5  # Best numerical accuracy

# 3. Narrow sweep based on promising parameters (e.g., best was klen=20, window=200)
parameter_sweep.sh -b broad_sweep_bedtools_ref.tsv -f fasta_files.txt --klen 18,20,22 --window 150,200,250 --precision 20,21,22,23 -o narrow_sweep/

# 4. Fine-tune around best parameters (e.g., best was klen=20, window=200, precision=22)
parameter_sweep.sh -b broad_sweep_bedtools_ref.tsv -f fasta_files.txt --klen 19,20,21 --window 190,200,210 --precision 21,22,23 -o fine_tune/

# 5. Compare custom results with default parameter ranges
complete_parameter_sweep.sh -b files.txt -f fasta_files.txt -o default_sweep/ -c -v
compare_sim_matrices.py fine_tune/best_result.csv default_sweep_results/best_default.csv

# 6. Test single optimal combination extensively
hammock fasta_files.txt fasta_files.txt -o final_test -k 20 -w 200 -p 22 --mode D
compare_clustering_trees.py final_test_mnmzr_p22_jaccD_k20_w200.csv broad_sweep_bedtools_ref.tsv --verbose
```

### Manual Analysis Workflow
```bash
# 1. Generate bedtools pairwise jaccard reference matrix
bedtools_pairwise.sh -f files.txt -o bedtools_ref.txt

# 2. Run hammock with specific parameters
hammock files.txt files.txt -o hammock_results -p 20 --mode B

# 3. Compare similarity matrices (numerical accuracy)
compare_sim_matrices.py hammock_results_hll_p20_jaccB.csv bedtools_ref.txt

# 4. Compare clustering tree topologies (biological accuracy)
compare_clustering_trees.py hammock_results_hll_p20_jaccB.csv bedtools_ref.txt --verbose

# 5. Compare different parameter combinations
for p in 18 20 22 25; do
    echo "Testing precision $p"
    hammock files.txt files.txt -o hammock_p${p} -p $p --mode B
    echo -n "Numerical: "; compare_sim_matrices.py hammock_p${p}_hll_p${p}_jaccB.csv bedtools_ref.txt --table
    echo -n "Biological: "; compare_clustering_trees.py hammock_p${p}_hll_p${p}_jaccB.csv bedtools_ref.txt --table
done
```

### Clustering Tree Analysis Workflow
```bash
# 1. Generate dendrograms for visual inspection (bedtools reference automatically generated by parameter sweep)
draw_clustering_tree.py bedtools_ref.tsv --output bedtools_tree.png --summary
draw_clustering_tree.py hammock_output.csv --output hammock_tree.png --summary

# 2. Compare dendrograms with different linkage methods
for method in single complete average ward; do
    draw_clustering_tree.py bedtools_ref.tsv --linkage $method --output bedtools_${method}.png --no-show --quiet
    draw_clustering_tree.py hammock_output.csv --linkage $method --output hammock_${method}.png --no-show --quiet
done

# 3. Quantitative tree comparison with detailed output
compare_clustering_trees.py hammock_output.csv bedtools_ref.tsv --verbose \
    --save-trees trees.newick --output detailed_results.txt

# 4. Visualize saved Newick trees from clustering comparison
draw_clustering_tree.py trees.newick --tree-index 0 --output hammock_tree_from_newick.png --title "Hammock Tree"
draw_clustering_tree.py trees.newick --tree-index 1 --output bedtools_tree_from_newick.png --title "Bedtools Reference"

# 5. Test different linkage methods to see which preserves topology best
compare_clustering_trees.py hammock_output.csv bedtools_ref.tsv --linkage single --table
compare_clustering_trees.py hammock_output.csv bedtools_ref.tsv --linkage complete --table
compare_clustering_trees.py hammock_output.csv bedtools_ref.tsv --linkage average --table
compare_clustering_trees.py hammock_output.csv bedtools_ref.tsv --linkage ward --table

# 6. Batch comparison for parameter sweep results
echo "file1	file2	format1	format2	matrix_size	rf_distance	max_rf_distance	normalized_rf	linkage_method" > clustering_results.tsv
for hammock_file in parameter_sweep_results/*.csv; do
    compare_clustering_trees.py "$hammock_file" bedtools_ref.tsv --table >> clustering_results.tsv
done

# 7. Generate dendrograms for best parameter combinations
sort -k8 -n clustering_results.tsv | head -3 | while read line; do
    hammock_file=$(echo "$line" | cut -f1)
    normalized_rf=$(echo "$line" | cut -f8)
    base_name=$(basename "$hammock_file" .csv)
    draw_clustering_tree.py "$hammock_file" --output "best_${base_name}_rf${normalized_rf}.png" --no-show
done

# 8. Find the best parameter combination (lowest normalized RF distance)
sort -k8 -n clustering_results.tsv | head -5
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

### Python Dependencies for Visualization
- **matplotlib**: For dendrogram plotting (required for `draw_clustering_tree.py` and automatic dendrogram generation)
- **scipy**: For hierarchical clustering algorithms
- **pandas**: For data manipulation
- **numpy**: For numerical operations

### Optional Python Libraries
- **ete3**: For accurate Robinson-Foulds distance calculations (recommended)
  ```bash
  conda install ete3
  # or
  pip install ete3
  ```
- **dendropy**: Alternative to ete3 for tree comparisons
  ```bash
  pip install dendropy
  ```

### Script Dependencies
Scripts that call other scripts will automatically find them after installation:
- `complete_parameter_sweep*.sh` → `parameter_sweep*.sh`
- `parameter_sweep*.sh` → `compare_sim_matrices.py`, `compare_clustering_trees.py`, `draw_clustering_tree.py`

## Output Files

### Parameter Sweep Results

#### Similarity Matrix Results (`*_results.tsv`)
Results are saved in TSV format with columns:
- `klen`, `window`, `precision`: Parameter values
- `hammock_file`, `bedtools_file`: Input files
- `format1`, `format2`: File format detection
- `matrix_size`: Number of files compared
- `frobenius_norm`: Matrix difference metric
- `correlation`: Pearson correlation coefficient
- `mean_abs_error`, `max_abs_error`: Error metrics
- `runtime_seconds`: Execution time

#### Clustering Tree Results (`*_results_clustering.tsv`)
**NEW**: Clustering tree topology comparison results (when `compare_clustering_trees.py` is available):
- `file1`, `file2`: Input files (hammock vs bedtools)
- `format1`, `format2`: File format detection
- `matrix_size`: Number of files compared
- `rf_distance`: Robinson-Foulds distance (number of different splits)
- `max_rf_distance`: Maximum possible RF distance for n leaves
- `normalized_rf`: RF distance / max RF distance (0-1 scale)
- `linkage_method`: Clustering method used (typically 'average')

**Interpretation of Clustering Results**:
- **Normalized RF < 0.2**: Excellent biological accuracy (trees very similar)
- **Normalized RF 0.2-0.5**: Good biological accuracy (moderate similarity)
- **Normalized RF 0.5-0.8**: Poor biological accuracy (trees quite different)
- **Normalized RF > 0.8**: Very poor biological accuracy (trees very different)

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
5. **"Clustering comparison unavailable"**: Install ete3 (`conda install ete3`) for accurate Robinson-Foulds calculations
6. **"Different best parameters"**: Numerical accuracy (similarity matrices) and biological accuracy (clustering trees) may suggest different optimal parameters - consider your analysis goals
7. **"No display available"**: Use `--no-show` flag with `draw_clustering_tree.py` when running on headless systems
8. **"Matplotlib backend error"**: Set `MPLBACKEND=Agg` environment variable for non-interactive plotting
9. **"Empty or incomplete results tables"**: Use debug mode (`-d` or `--debug`) to investigate table writing issues and preserve intermediate files
10. **"Clustering results missing"**: Enable debug mode to see clustering comparison exit codes and error messages
11. **"Custom parameter format error"**: Ensure custom parameter values are comma-separated without spaces (e.g., '10,15,20' not '10, 15, 20')
12. **"Too many parameter combinations"**: Use custom parameters to limit the search space instead of full parameter sweeps

### Getting Help
```bash
# Most scripts support help flags
script_name.sh -h
script_name.py --help

# See custom parameter options
complete_parameter_sweep.sh -h
parameter_sweep.sh -h

# Use debug mode for troubleshooting
complete_parameter_sweep.sh -b files.txt -o debug_test -d -v
parameter_sweep.sh -b bedtools_ref.tsv -f files.txt --debug

# Test custom parameters with debug mode
complete_parameter_sweep.sh -b files.txt --precision 20,22 -d -v
parameter_sweep.sh -b bedtools_ref.tsv -f files.txt --klen 15,20 --debug
```

---

*This README was generated for hammock v0.2.1. For the latest documentation, visit the [hammock GitHub repository](https://github.com/jessicabonnie/hammock).* 