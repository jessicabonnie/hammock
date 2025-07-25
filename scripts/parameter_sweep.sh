#!/bin/bash

# Parameter sweep script for hammock
# Runs hammock with different parameter combinations and compares to bedtools output
# Automatically detects BED vs FASTA files and adjusts parameters accordingly

set -e  # Exit on any error

# Initialize error logging variables early
ERROR_LOG=""
SCRIPT_START_TIME=$(date)

# Function to log errors and exit gracefully
error_handler() {
    local exit_code=$?
    local line_number=$1
    local command="$2"
    
    if [[ -n "$ERROR_LOG" ]]; then
        {
            echo ""
            echo "FATAL ERROR - $(date)"
            echo "===========================================" 
            echo "Exit code: $exit_code"
            echo "Line number: $line_number"
            echo "Failed command: $command"
            echo "Current working directory: $(pwd)"
            echo "Script arguments: $0 $*"
            echo "Environment PATH: $PATH"
            echo ""
            echo "Script execution details:"
            echo "  Start time: $SCRIPT_START_TIME"
            echo "  End time: $(date)"
            echo "  Bedtools file: ${BEDTOOLS_FILE:-'not set'}"
            echo "  File list: ${FILE_LIST:-'not set'}"
            echo "  Output directory: ${OUTPUT_DIR:-'not set'}"
            echo "  Quick mode: ${QUICK_MODE:-'not set'}"
            echo ""
            echo "Available commands check:"
            echo "  bedtools: $(command -v bedtools 2>/dev/null || echo 'NOT FOUND')"
            echo "  hammock: $(command -v hammock 2>/dev/null || echo 'NOT FOUND')"
            echo "  compare_sim_matrices.py: $(command -v compare_sim_matrices.py 2>/dev/null || echo 'NOT FOUND')"
            echo ""
            echo "Recent shell history (last 10 commands):"
            history | tail -10 2>/dev/null || echo "History not available"
            echo ""
        } >> "$ERROR_LOG"
    fi
    
    echo "FATAL ERROR: Script failed at line $line_number with exit code $exit_code" >&2
    echo "Check error log for details: ${ERROR_LOG:-'Error log not initialized'}" >&2
    exit $exit_code
}

# Set up error trapping
trap 'error_handler $LINENO "$BASH_COMMAND"' ERR

# Default values - modify these as needed
BEDTOOLS_FILE=""
FILE_LIST=""
OUTPUT_DIR="parameter_sweep_results"
RESULTS_TABLE="parameter_sweep_results.tsv"
QUICK_MODE=false
DEBUG_MODE=false

# Custom parameter ranges (empty means use defaults)
CUSTOM_WINDOW=""
CUSTOM_KLEN=""
CUSTOM_PRECISION=""

# Parameter arrays - will be set based on file type detection
KLEN_VALUES=()
WINDOW_VALUES=()
PRECISION_VALUES=()
HAMMOCK_MODE=""
FILE_TYPE=""

# Function to display usage
usage() {
    echo "Usage: $0 -b <bedtools_output_file> -f <file_list> [-o <output_dir>] [-r <results_table>] [-q] [-d] [--window <values>] [--klen <values>] [--precision <values>]"
    echo ""
    echo "Required arguments:"
    echo "  -b <bedtools_output_file>  Path to bedtools pairwise jaccard output matrix file"
    echo "                             (NOT a list of bed files - must be the actual bedtools output)"
    echo "  -f <file_list>             File containing list of BED or FASTA file paths (one per line)"
    echo ""
    echo "Optional arguments:"
    echo "  -o <output_dir>            Output directory for hammock results (default: parameter_sweep_results)"
    echo "  -r <results_table>         Results table filename (default: parameter_sweep_results.tsv)"
    echo "  -q, --quick                Quick mode with limited parameter combinations for testing"
    echo "  -d, --debug                Debug mode with detailed output and preserved intermediate files"
    echo "  --window <values>          Custom window sizes (comma-separated, e.g., '10,20,50')"
    echo "  --klen <values>            Custom k-mer lengths (comma-separated, e.g., '10,15,20')"
    echo "  --precision <values>       Custom precision values (comma-separated, e.g., '20,22,24')"
    echo "  -h                         Show this help message"
    echo ""
    echo "File type detection:"
    echo "  BED files (.bed): Only precision parameter is varied (mode B)"
    echo "  FASTA files (.fa, .fasta, etc.): All parameters varied (mode D)"
    echo ""
    echo "Parameter ranges:"
    echo "  Full mode (default):"
    echo "    BED files: precision = 16,19,20,21,22,23,24,26,28,30"
    echo "    FASTA files: klen = 10,15,20,25,30; window = 25,50,100,200,300,400,500; precision = 16,19,20,21,22,23,24,26,28,30"
    echo "  Quick mode (-q):"
    echo "    BED files: precision = 20,23"
    echo "    FASTA files: klen = 20; window = 200; precision = 20,23"
    echo "  Custom mode:"
    echo "    Use --window, --klen, and --precision to override defaults"
    echo "    Custom parameters take precedence over quick/full mode settings"
    echo ""
    echo "IMPORTANT: You must first generate bedtools output before running this script:"
    echo "  For BED files:"
    echo "    bedtools jaccard -a file1.bed -b file2.bed > bedtools_output.txt"
    echo "    (repeat for all pairwise combinations)"
    echo "  For FASTA files:"
    echo "    Generate k-mer sketches first, then compute jaccard similarities"
    echo ""
    echo "Examples:"
    echo "  # First generate bedtools output from your bed files"
    echo "  generate_bedtools_matrix.py bed_files_list.txt bedtools_output.txt"
    echo "  # Then run parameter sweep"
    echo "  $0 -b bedtools_output.txt -f bed_files_list.txt"
    echo ""
    echo "  # Quick mode for testing"
    echo "  $0 -b bedtools_output.txt -f bed_files_list.txt --quick"
    echo ""
    echo "  # Custom parameter ranges"
    echo "  $0 -b bedtools_output.txt -f fasta_files_list.txt --window 8,10,20,50 --klen 8,15,20 --precision 20,22,24"
    echo ""
    echo "  # For FASTA files (assuming you have bedtools output)"
    echo "  $0 -b bedtools_fasta_output.txt -f fasta_files_list.txt"
}

# Function to detect file type from file list
detect_file_type() {
    local file_list="$1"
    
    # Read first non-empty, non-comment line
    while IFS= read -r file_path; do
        if [[ -n "$file_path" && ! "$file_path" =~ ^[[:space:]]*# ]]; then
            file_path=$(echo "$file_path" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            
            # Check file extension (including gzipped versions)
            if [[ "$file_path" =~ \.(fa|fasta|fna|ffn|faa|frn)(\.gz)?$ ]]; then
                echo "FASTA"
                return 0
            elif [[ "$file_path" =~ \.bed(\.gz)?$ ]]; then
                echo "BED"
                return 0
            fi
        fi
    done < "$file_list"
    
    echo "UNKNOWN"
    return 1
}

# Function to set parameters based on file type
set_parameters() {
    local file_type="$1"
    local quick_mode="$2"
    
    # First set default parameters based on file type and mode
    if [[ "$file_type" == "FASTA" ]]; then
        # FASTA files: test all parameters (mode D)
        if [[ "$quick_mode" == "true" ]]; then
            # Quick mode: limited parameters for testing
            KLEN_VALUES=(20)
            WINDOW_VALUES=(200)
            PRECISION_VALUES=(20 23)
            echo "Detected FASTA files - using mode D with limited parameter sweep (quick mode)"
        else
            # Full mode: comprehensive parameter sweep
            KLEN_VALUES=(10 15 20 25 30)
            WINDOW_VALUES=(25 50 100 200 300 400 500)
            PRECISION_VALUES=(16 19 20 21 22 23 24 26 28 30)
            echo "Detected FASTA files - using mode D with full parameter sweep"
        fi
        HAMMOCK_MODE="D"
    elif [[ "$file_type" == "BED" ]]; then
        # BED files: only test precision (mode B)
        KLEN_VALUES=(0)  # Not used for BED files
        WINDOW_VALUES=(0)  # Not used for BED files
        if [[ "$quick_mode" == "true" ]]; then
            # Quick mode: limited precision values
            PRECISION_VALUES=(20 23)
            echo "Detected BED files - using mode B with limited precision sweep (quick mode)"
        else
            # Full mode: comprehensive precision sweep
            PRECISION_VALUES=(16 19 20 21 22 23 24 26 28 30)
            echo "Detected BED files - using mode B with precision-only sweep"
        fi
        HAMMOCK_MODE="B"
    else
        echo "Error: Could not detect file type from file list" >&2
        exit 1
    fi
    
    # Override with custom parameters if provided
    if [[ -n "$CUSTOM_KLEN" ]]; then
        IFS=',' read -ra KLEN_VALUES <<< "$CUSTOM_KLEN"
        echo "Using custom k-mer lengths: ${KLEN_VALUES[*]}"
    fi
    if [[ -n "$CUSTOM_WINDOW" ]]; then
        IFS=',' read -ra WINDOW_VALUES <<< "$CUSTOM_WINDOW"
        echo "Using custom window sizes: ${WINDOW_VALUES[*]}"
    fi
    if [[ -n "$CUSTOM_PRECISION" ]]; then
        IFS=',' read -ra PRECISION_VALUES <<< "$CUSTOM_PRECISION"
        echo "Using custom precision values: ${PRECISION_VALUES[*]}"
    fi
}

# Parse command line arguments
while getopts "b:f:o:r:qdh-:" opt; do
    case $opt in
        b)
            BEDTOOLS_FILE="$OPTARG"
            ;;
        f)
            FILE_LIST="$OPTARG"
            ;;
        o)
            OUTPUT_DIR="$OPTARG"
            ;;
        r)
            RESULTS_TABLE="$OPTARG"
            ;;
        q)
            QUICK_MODE=true
            ;;
        d)
            DEBUG_MODE=true
            ;;
        h)
            usage
            exit 0
            ;;
        -)
            case "$OPTARG" in
                quick)
                    QUICK_MODE=true
                    ;;
                debug)
                    DEBUG_MODE=true
                    ;;
                window)
                    CUSTOM_WINDOW="${!OPTIND}"; OPTIND=$(( OPTIND + 1 ))
                    ;;
                klen)
                    CUSTOM_KLEN="${!OPTIND}"; OPTIND=$(( OPTIND + 1 ))
                    ;;
                precision)
                    CUSTOM_PRECISION="${!OPTIND}"; OPTIND=$(( OPTIND + 1 ))
                    ;;
                *)
                    echo "Invalid long option: --$OPTARG" >&2
                    usage
                    exit 1
                    ;;
            esac
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
    esac
done

# Check required arguments
if [[ -z "$BEDTOOLS_FILE" || -z "$FILE_LIST" ]]; then
    echo "Error: Both -b (bedtools file) and -f (file list) are required" >&2
    usage
    exit 1
fi

# Validate inputs
if [[ ! -f "$BEDTOOLS_FILE" ]]; then
    echo "Error: Bedtools file '$BEDTOOLS_FILE' does not exist" >&2
    exit 1
fi

if [[ ! -f "$FILE_LIST" ]]; then
    echo "Error: File list '$FILE_LIST' does not exist" >&2
    exit 1
fi

# Detect file type and set parameters
FILE_TYPE=$(detect_file_type "$FILE_LIST")
set_parameters "$FILE_TYPE" "$QUICK_MODE"

# Validate that all files in the list exist
echo "Validating file list..."
while IFS= read -r file_path; do
    # Skip empty lines and comments
    if [[ -n "$file_path" && ! "$file_path" =~ ^[[:space:]]*# ]]; then
        # Remove leading/trailing whitespace
        file_path=$(echo "$file_path" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        
        if [[ ! -f "$file_path" ]]; then
            echo "Error: File '$file_path' from list does not exist" >&2
            exit 1
        fi
    fi
done < "$FILE_LIST"

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Test if compare_sim_matrices.py is available and working
echo "Testing compare_sim_matrices.py availability..."
if [[ -f "$SCRIPT_DIR/compare_sim_matrices.py" ]]; then
    echo "Found compare_sim_matrices.py in script directory"
    if python "$SCRIPT_DIR/compare_sim_matrices.py" --help >/dev/null 2>&1; then
        echo "compare_sim_matrices.py is working correctly"
    else
        echo "WARNING: compare_sim_matrices.py found but not working properly" >&2
    fi
elif command -v compare_sim_matrices.py >/dev/null 2>&1; then
    echo "Found compare_sim_matrices.py in PATH"
    if compare_sim_matrices.py --help >/dev/null 2>&1; then
        echo "compare_sim_matrices.py is working correctly"
    else
        echo "WARNING: compare_sim_matrices.py found in PATH but not working properly" >&2
    fi
else
    echo "ERROR: compare_sim_matrices.py not found in script directory or PATH" >&2
    echo "Script directory: $SCRIPT_DIR" >&2
    echo "Contents of script directory:" >&2
    ls -la "$SCRIPT_DIR" >&2
    exit 1
fi

# Test if compare_clustering_trees.py is available and working
echo "Testing compare_clustering_trees.py availability..."
CLUSTERING_COMPARISON_AVAILABLE=false
if [[ -f "$SCRIPT_DIR/compare_clustering_trees.py" ]]; then
    echo "Found compare_clustering_trees.py in script directory"
    if python "$SCRIPT_DIR/compare_clustering_trees.py" --help >/dev/null 2>&1; then
        echo "compare_clustering_trees.py is working correctly"
        CLUSTERING_COMPARISON_AVAILABLE=true
    else
        echo "WARNING: compare_clustering_trees.py found but not working properly" >&2
    fi
elif command -v compare_clustering_trees.py >/dev/null 2>&1; then
    echo "Found compare_clustering_trees.py in PATH"
    if compare_clustering_trees.py --help >/dev/null 2>&1; then
        echo "compare_clustering_trees.py is working correctly"
        CLUSTERING_COMPARISON_AVAILABLE=true
    else
        echo "WARNING: compare_clustering_trees.py found in PATH but not working properly" >&2
    fi
else
    echo "WARNING: compare_clustering_trees.py not found - clustering tree comparisons will be skipped" >&2
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Initialize error log
ERROR_LOG="$OUTPUT_DIR/parameter_sweep_error.log"
{
    echo "Parameter Sweep Error Log - $(date)"
    echo "========================================"
    echo "Script: $0"
    echo "Arguments: $*"
    echo "Working directory: $(pwd)"
    echo "User: $(whoami)"
    echo "Host: $(hostname)"
    echo ""
    echo "Configuration:"
    echo "  Bedtools file: ${BEDTOOLS_FILE}"
    echo "  File list: ${FILE_LIST}"
    echo "  Output directory: ${OUTPUT_DIR}"
    echo "  Results table: ${RESULTS_TABLE}"
    echo "  Quick mode: ${QUICK_MODE}"
    echo ""
    echo "Environment check:"
    echo "  PATH: $PATH"
    echo "  bedtools: $(command -v bedtools 2>/dev/null || echo 'NOT FOUND')"
    echo "  hammock: $(command -v hammock 2>/dev/null || echo 'NOT FOUND')"
    echo "  compare_sim_matrices.py: $(command -v compare_sim_matrices.py 2>/dev/null || echo 'NOT FOUND')"
    echo ""
    echo "Log entries:"
    echo "============"
} > "$ERROR_LOG"

# Initialize results table with header
echo "Starting parameter sweep..."
echo "Bedtools file: $BEDTOOLS_FILE"
echo "File list: $FILE_LIST"
echo "File type: $FILE_TYPE"
echo "Hammock mode: $HAMMOCK_MODE"
echo "Quick mode: $QUICK_MODE"
echo "Output directory: $OUTPUT_DIR"
echo "Similarity matrix results table: $RESULTS_TABLE"

# Set clustering results table path if clustering comparison is available
if [[ "$CLUSTERING_COMPARISON_AVAILABLE" == "true" ]]; then
    CLUSTERING_RESULTS_TABLE="${RESULTS_TABLE%.tsv}_clustering.tsv"
    echo "Clustering tree results table: $CLUSTERING_RESULTS_TABLE"
    echo "Clustering comparison: ENABLED (using average linkage)"
else
    echo "Clustering comparison: DISABLED (compare_clustering_trees.py not available)"
fi
echo ""

# Create results table with extended header
echo -e "klen\twindow\tprecision\thammock_file\tbedtools_file\tformat1\tformat2\tmatrix_size\tfrobenius_norm\tcorrelation\tmean_abs_error\tmax_abs_error\truntime_seconds" > "$RESULTS_TABLE"

# Create clustering results table if clustering comparison is available
if [[ "$CLUSTERING_COMPARISON_AVAILABLE" == "true" ]]; then
    echo "Clustering results will be saved to: $CLUSTERING_RESULTS_TABLE"
    echo -e "klen\twindow\tprecision\thammock_file\tbedtools_file\tformat1\tformat2\tmatrix_size\trf_distance\tmax_rf_distance\tnormalized_rf\tlinkage_method\tclustering_method\tn_clusters1\tn_clusters2\tadjusted_rand_index\tnormalized_mutual_info\tcluster_correspondence\tcluster_stability\truntime_seconds" > "$CLUSTERING_RESULTS_TABLE"
fi

# Count total combinations for progress tracking
total_combinations=$((${#KLEN_VALUES[@]} * ${#WINDOW_VALUES[@]} * ${#PRECISION_VALUES[@]}))
current_combination=0

echo "Testing $total_combinations parameter combinations..."
echo ""

# Loop through all parameter combinations
for klen in "${KLEN_VALUES[@]}"; do
    for window in "${WINDOW_VALUES[@]}"; do
        for precision in "${PRECISION_VALUES[@]}"; do
            current_combination=$((current_combination + 1))
            
            if [[ "$FILE_TYPE" == "BED" ]]; then
                echo "[$current_combination/$total_combinations] Testing precision=$precision"
                # Generate base output filename (hammock will add suffix)
                hammock_base="$OUTPUT_DIR/hammock"
                # Run hammock with precision only
                hammock_cmd="hammock '$FILE_LIST' '$FILE_LIST' -o '$hammock_base' -p $precision --mode $HAMMOCK_MODE"
            else
                echo "[$current_combination/$total_combinations] Testing klen=$klen, window=$window, precision=$precision"
                # Generate base output filename (hammock will add suffix)
                hammock_base="$OUTPUT_DIR/hammock"
                # Run hammock with all parameters
                hammock_cmd="hammock '$FILE_LIST' '$FILE_LIST' -o '$hammock_base' -k $klen -w $window -p $precision --mode $HAMMOCK_MODE --minimizer"
            fi
            
            # Record start time
            start_time=$(date +%s)
            
            # Run hammock with current parameters (pairwise comparison using same file list twice)
            echo "  Running hammock..."
            echo "  Command: $hammock_cmd"
            hammock_error_output=$(mktemp)
            hammock_stdout_output=$(mktemp)
            if ! eval "$hammock_cmd" > "$hammock_stdout_output" 2>"$hammock_error_output"; then
                echo "  ERROR: Hammock failed with parameters klen=$klen, window=$window, precision=$precision" >&2
                echo "  ERROR: Hammock command: $hammock_cmd" >&2
                echo "  ERROR: Current working directory: $(pwd)" >&2
                echo "  ERROR: Hammock stdout:" >&2
                cat "$hammock_stdout_output" >&2
                echo "  ERROR: Hammock stderr:" >&2
                cat "$hammock_error_output" >&2
                echo "  ERROR: File list contents (first 5 lines):" >&2
                head -5 "$FILE_LIST" >&2
                echo -e "$klen\t$window\t$precision\tHAMMOCK_FAILED\t-\t-\t-\t-\t-\t-\t-\t-\t$runtime" >> "$RESULTS_TABLE"
                rm -f "$hammock_error_output" "$hammock_stdout_output"
                continue
            fi
            rm -f "$hammock_error_output" "$hammock_stdout_output"
            
            # Record end time and calculate runtime
            end_time=$(date +%s)
            runtime=$((end_time - start_time))
            
            # Determine actual output filename (hammock adds suffix based on parameters)
            if [[ "$HAMMOCK_MODE" == "B" ]]; then
                # For mode B with hyperloglog (default), precision p, the suffix is: _hll_p{precision}_jaccB.csv
                hammock_output="${hammock_base}_hll_p${precision}_jaccB.csv"
            else
                # For mode D with minimizer (default for FASTA), the suffix includes kmer and window
                hammock_output="${hammock_base}_mnmzr_p${precision}_jaccD_k${klen}_w${window}.csv"
            fi
            
            # Check if hammock output was created
            if [[ ! -f "$hammock_output" ]]; then
                echo "  ERROR: Hammock output file was not created: $hammock_output" >&2
                echo -e "$klen\t$window\t$precision\tNO_OUTPUT\t-\t-\t-\t-\t-\t-\t-\t-\t$runtime" >> "$RESULTS_TABLE"
                continue
            fi
            
            # Compare with bedtools output
            echo "  Comparing with bedtools output..."
            echo "  Hammock output file: $hammock_output"
            echo "  Bedtools file: $BEDTOOLS_FILE"
            
            # Check if files exist before comparison
            if [[ ! -f "$hammock_output" ]]; then
                echo "  ERROR: Hammock output file missing: $hammock_output" >&2
                echo -e "$klen\t$window\t$precision\tMISSING_HAMMOCK_OUTPUT\t-\t-\t-\t-\t-\t-\t-\t-\t$runtime" >> "$RESULTS_TABLE"
                continue
            fi
            
            if [[ ! -f "$BEDTOOLS_FILE" ]]; then
                echo "  ERROR: Bedtools file missing: $BEDTOOLS_FILE" >&2
                echo -e "$klen\t$window\t$precision\tMISSING_BEDTOOLS_FILE\t-\t-\t-\t-\t-\t-\t-\t-\t$runtime" >> "$RESULTS_TABLE"
                continue
            fi
            
            # Show file sizes for debugging
            echo "  Hammock output size: $(wc -l < "$hammock_output") lines"
            echo "  Bedtools file size: $(wc -l < "$BEDTOOLS_FILE") lines"
            
            # Try to find compare_sim_matrices.py in the same directory first, then fallback to PATH
            comparison_error_output=$(mktemp)
            comparison_stdout_output=$(mktemp)
            
            if [[ -f "$SCRIPT_DIR/compare_sim_matrices.py" ]]; then
                comparison_cmd="python $SCRIPT_DIR/compare_sim_matrices.py $hammock_output $BEDTOOLS_FILE --table"
                echo "  Comparison command: $comparison_cmd"
                comparison_result=$($comparison_cmd > "$comparison_stdout_output" 2>"$comparison_error_output")
                comparison_exit_code=$?
            else
                comparison_cmd="compare_sim_matrices.py $hammock_output $BEDTOOLS_FILE --table"
                echo "  Comparison command: $comparison_cmd"
                comparison_result=$($comparison_cmd > "$comparison_stdout_output" 2>"$comparison_error_output")
                comparison_exit_code=$?
            fi
            
            # Read the actual stdout content
            if [[ -f "$comparison_stdout_output" ]]; then
                comparison_result=$(cat "$comparison_stdout_output")
            fi
            
            if [[ $comparison_exit_code -eq 0 && -n "$comparison_result" ]]; then
                # Add parameter columns and runtime to the comparison result
                echo -e "$klen\t$window\t$precision\t$comparison_result\t$runtime" >> "$RESULTS_TABLE"
                echo "  ✓ Similarity matrix comparison completed successfully"
                
                # Perform clustering tree comparison if available
                if [[ "$CLUSTERING_COMPARISON_AVAILABLE" == "true" ]]; then
                    echo "  Comparing clustering trees..."
                    clustering_error_output=$(mktemp)
                    clustering_stdout_output=$(mktemp)
                    
                    if [[ -f "$SCRIPT_DIR/compare_clustering_trees.py" ]]; then
                        clustering_cmd="python $SCRIPT_DIR/compare_clustering_trees.py $hammock_output $BEDTOOLS_FILE --table --linkage average"
                    else
                        clustering_cmd="compare_clustering_trees.py $hammock_output $BEDTOOLS_FILE --table --linkage average"
                    fi
                    
                    echo "  Clustering comparison command: $clustering_cmd"
                    clustering_result=$($clustering_cmd > "$clustering_stdout_output" 2>"$clustering_error_output")
                    clustering_exit_code=$?
                    
                    # Read the actual stdout content
                    if [[ -f "$clustering_stdout_output" ]]; then
                        clustering_result=$(cat "$clustering_stdout_output")
                    fi
                    
                    # DEBUG: Show what we captured
                    if [[ "$DEBUG_MODE" == "true" ]]; then
                        echo "  DEBUG: Clustering exit code: $clustering_exit_code"
                        echo "  DEBUG: Clustering result length: ${#clustering_result}"
                        echo "  DEBUG: Clustering result content: '$clustering_result'"
                    fi
                    
                    if [[ $clustering_exit_code -eq 0 && -n "$clustering_result" ]]; then
                        # Add parameter columns and runtime to the clustering comparison result
                        clustering_line="$klen\t$window\t$precision\t$clustering_result\t$runtime"
                        if [[ "$DEBUG_MODE" == "true" ]]; then
                            echo "  DEBUG: Writing clustering line: '$clustering_line'"
                            echo "  DEBUG: Clustering table path: '$CLUSTERING_RESULTS_TABLE'"
                            echo "  DEBUG: Table exists before write: $(test -f "$CLUSTERING_RESULTS_TABLE" && echo "YES" || echo "NO")"
                            echo "  DEBUG: Table size before write: $(wc -l < "$CLUSTERING_RESULTS_TABLE" 2>/dev/null || echo "ERROR") lines"
                        fi
                        echo -e "$clustering_line" >> "$CLUSTERING_RESULTS_TABLE"
                        if [[ "$DEBUG_MODE" == "true" ]]; then
                            echo "  DEBUG: Write operation exit code: $?"
                            echo "  DEBUG: Table size after write: $(wc -l < "$CLUSTERING_RESULTS_TABLE" 2>/dev/null || echo "ERROR") lines"
                            echo "  DEBUG: Last few lines of table:"
                            tail -3 "$CLUSTERING_RESULTS_TABLE" 2>/dev/null | sed 's/^/    /' || echo "    ERROR reading table"
                        fi
                        echo "  ✓ Clustering tree comparison completed successfully"
                    else
                        echo "  WARNING: Clustering tree comparison failed" >&2
                        echo "  WARNING: Clustering command: $clustering_cmd" >&2
                        echo "  WARNING: Clustering command exit code: $clustering_exit_code" >&2
                        echo "  WARNING: Clustering stdout:" >&2
                        cat "$clustering_stdout_output" >&2
                        echo "  WARNING: Clustering stderr:" >&2
                        cat "$clustering_error_output" >&2
                        echo -e "$klen\t$window\t$precision\tCLUSTERING_FAILED\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t$runtime" >> "$CLUSTERING_RESULTS_TABLE"
                    fi
                    rm -f "$clustering_error_output" "$clustering_stdout_output"
                fi
            else
                echo "  ERROR: Similarity matrix comparison failed" >&2
                echo "  ERROR: Comparison command: $comparison_cmd" >&2
                echo "  ERROR: Comparison command exit code: $comparison_exit_code" >&2
                echo "  ERROR: Comparison stdout:" >&2
                cat "$comparison_stdout_output" >&2
                echo "  ERROR: Comparison stderr:" >&2
                cat "$comparison_error_output" >&2
                echo "  ERROR: First 5 lines of hammock output:" >&2
                head -5 "$hammock_output" >&2
                echo "  ERROR: First 5 lines of bedtools file:" >&2
                head -5 "$BEDTOOLS_FILE" >&2
                echo -e "$klen\t$window\t$precision\tCOMPARISON_FAILED\t-\t-\t-\t-\t-\t-\t-\t-\t$runtime" >> "$RESULTS_TABLE"
                
                # Also add failed entry to clustering results if available
                if [[ "$CLUSTERING_COMPARISON_AVAILABLE" == "true" ]]; then
                    echo -e "$klen\t$window\t$precision\tSIM_COMPARISON_FAILED\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t$runtime" >> "$CLUSTERING_RESULTS_TABLE"
                fi
            fi
            rm -f "$comparison_error_output" "$comparison_stdout_output"
            
            # Clean up hammock output file (skip during debug mode)
            if [[ "$DEBUG_MODE" == "true" ]]; then
                echo "  DEBUG MODE: Preserving hammock output file: $hammock_output"
            else
                echo "  Cleaning up hammock output..."
                rm -f "$hammock_output"
            fi
            
            echo "  Runtime: ${runtime}s"
            echo ""
        done
    done
done

echo "Parameter sweep completed!"
echo "Results saved to: $RESULTS_TABLE"
if [[ "$CLUSTERING_COMPARISON_AVAILABLE" == "true" ]]; then
    echo "Clustering tree results saved to: $CLUSTERING_RESULTS_TABLE"
fi
echo ""

# Generate dendrogram visualization of bedtools reference
echo "Generating bedtools reference dendrogram..."
DENDROGRAM_OUTPUT="$OUTPUT_DIR/bedtools_reference_dendrogram.png"

# Test if draw_clustering_tree.py is available
DRAW_TREE_AVAILABLE=false
if [[ -f "$SCRIPT_DIR/draw_clustering_tree.py" ]]; then
    echo "Found draw_clustering_tree.py in script directory"
    if python "$SCRIPT_DIR/draw_clustering_tree.py" --help >/dev/null 2>&1; then
        DRAW_TREE_AVAILABLE=true
        DRAW_TREE_CMD="python $SCRIPT_DIR/draw_clustering_tree.py"
    else
        echo "WARNING: draw_clustering_tree.py found but not working properly" >&2
    fi
elif command -v draw_clustering_tree.py >/dev/null 2>&1; then
    echo "Found draw_clustering_tree.py in PATH"
    if draw_clustering_tree.py --help >/dev/null 2>&1; then
        DRAW_TREE_AVAILABLE=true
        DRAW_TREE_CMD="draw_clustering_tree.py"
    else
        echo "WARNING: draw_clustering_tree.py found in PATH but not working properly" >&2
    fi
else
    echo "WARNING: draw_clustering_tree.py not found - dendrogram generation will be skipped" >&2
fi

if [[ "$DRAW_TREE_AVAILABLE" == "true" ]]; then
    echo "Generating dendrogram for bedtools reference file..."
    dendrogram_error_output=$(mktemp)
    dendrogram_stdout_output=$(mktemp)
    
    dendrogram_cmd="$DRAW_TREE_CMD '$BEDTOOLS_FILE' --output '$DENDROGRAM_OUTPUT' --no-show --title 'Bedtools Reference Dendrogram' --figsize 14 10"
    echo "Dendrogram command: $dendrogram_cmd"
    
    if eval "$dendrogram_cmd" > "$dendrogram_stdout_output" 2>"$dendrogram_error_output"; then
        echo "✓ Bedtools reference dendrogram saved to: $DENDROGRAM_OUTPUT"
    else
        echo "WARNING: Failed to generate bedtools reference dendrogram" >&2
        echo "WARNING: Dendrogram command: $dendrogram_cmd" >&2
        echo "WARNING: Dendrogram stdout:" >&2
        cat "$dendrogram_stdout_output" >&2
        echo "WARNING: Dendrogram stderr:" >&2
        cat "$dendrogram_error_output" >&2
    fi
    rm -f "$dendrogram_error_output" "$dendrogram_stdout_output"
else
    echo "Dendrogram generation skipped (draw_clustering_tree.py not available)"
fi
echo ""

echo "Summary:"
echo "  File type: $FILE_TYPE"
echo "  Hammock mode: $HAMMOCK_MODE"
echo "  Quick mode: $QUICK_MODE"
echo "  Total combinations tested: $total_combinations"
echo "  Similarity matrix results: $RESULTS_TABLE"
if [[ "$CLUSTERING_COMPARISON_AVAILABLE" == "true" ]]; then
    echo "  Clustering tree results: $CLUSTERING_RESULTS_TABLE"
fi
if [[ "$DRAW_TREE_AVAILABLE" == "true" && -f "$DENDROGRAM_OUTPUT" ]]; then
    echo "  Bedtools reference dendrogram: $DENDROGRAM_OUTPUT"
fi

# Show a preview of the results
echo ""
echo "Similarity matrix results preview (first 5 lines):"
head -n 6 "$RESULTS_TABLE" | column -t -s $'\t'

if [[ "$CLUSTERING_COMPARISON_AVAILABLE" == "true" && -f "$CLUSTERING_RESULTS_TABLE" ]]; then
    echo ""
    echo "Clustering tree results preview (first 5 lines):"
    head -n 6 "$CLUSTERING_RESULTS_TABLE" | column -t -s $'\t'
    
    # Show best clustering results (lowest normalized RF distance)
    echo ""
    echo "Best clustering approximations (lowest normalized RF distance):"
    if [[ $(wc -l < "$CLUSTERING_RESULTS_TABLE") -gt 1 ]]; then
        # Skip header and sort by normalized_rf column (column 11), then show top 3
        tail -n +2 "$CLUSTERING_RESULTS_TABLE" | sort -k11 -n | head -3 | column -t -s $'\t'
    else
        echo "No successful clustering comparisons found."
    fi
fi 