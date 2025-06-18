#!/bin/bash

# Complete parameter sweep script for hammock
# First generates bedtools reference, then runs hammock parameter sweep
# Automatically detects BED vs FASTA files and adjusts parameters accordingly

set -e  # Exit on any error

# Default values - modify these as needed
FASTA_FILE_LIST=""
BED_FILE_LIST=""
OUTPUT_PREFIX="complete_parameter_sweep"
CLEANUP=false
VERBOSE=false
QUICK_MODE=false

# Parameter arrays - will be set based on file type detection
KLEN_VALUES=()
WINDOW_VALUES=()
PRECISION_VALUES=()
HAMMOCK_MODE=""
FILE_TYPE=""

# Function to display usage
usage() {
    echo "Usage: $0 -b <bed_file_list> [-f <fasta_file_list>] [-o <output_prefix>] [-c] [-v] [-q]"
    echo ""
    echo "Complete parameter sweep for hammock:"
    echo "  1. Takes a list of BED files and runs pairwise bedtools jaccard on ALL combinations"
    echo "  2. Uses the resulting jaccard matrix as reference for hammock parameter sweep"
    echo "  3. Compares hammock output against bedtools reference for different parameters"
    echo ""
    echo "Required arguments:"
    echo "  -b <bed_file_list>    File containing list of BED file paths (one per line)"
    echo "                        Used for generating bedtools pairwise jaccard reference"
    echo ""
    echo "Optional arguments:"
    echo "  -f <fasta_file_list>  File containing list of FASTA file paths for hammock (one per line)"
    echo "                        If not provided, hammock will run on the BED files"
    echo "  -o <output_prefix>    Output prefix with path (default: complete_parameter_sweep)"
    echo "                        Creates: <prefix>_results/ directory, <prefix>_bedtools_ref.tsv, <prefix>_results.tsv"
    echo "  -c                    Clean up intermediate files after completion (bedtools reference always preserved)"
    echo "  -v                    Verbose output (shows progress of bedtools comparisons)"
    echo "  -q, --quick           Quick mode with limited parameter combinations for testing"
    echo "  -h                    Show this help message"
    echo ""
    echo "Mode selection:"
    echo "  BED files only:    Runs bedtools jaccard pairwise, then tests precision parameter (mode B)"
    echo "  FASTA files + BED: Uses BED files for reference, tests all parameters on FASTA (mode D)"
    echo ""
    echo "Parameter ranges:"
    echo "  Full mode (default):"
    echo "    BED files: precision = 16,18,19,20,21,22,23,25"
    echo "    FASTA files: klen = 15,20,25; window = 100,200,500; precision = 20,23,25"
    echo "  Quick mode (-q):"
    echo "    BED files: precision = 20,23"
    echo "    FASTA files: klen = 20; window = 200; precision = 20,23"
    echo ""
    echo "Examples:"
    echo "  # Run complete sweep on BED files only"
    echo "  $0 -b bed_files_list.txt -c -v"
    echo "  # Run sweep on FASTA files using corresponding BED files for reference"
    echo "  $0 -b bed_files_list.txt -f fasta_files_list.txt -c -v"
    echo "  # Quick mode for testing"
    echo "  $0 -b bed_files_list.txt --quick -c -v"
    echo "  # Custom output prefix and path"
    echo "  $0 -b bed_files_list.txt -o /path/to/my_experiment -c -v"
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
            PRECISION_VALUES=(16 18 19 20 21 22 23 25)
            echo "Detected BED files - using mode B with precision-only sweep"
        fi
        HAMMOCK_MODE="B"
    else
        echo "Error: Could not detect file type from file list" >&2
        exit 1
    fi
}

# Parse command line arguments
while getopts "b:f:o:cvqh-:" opt; do
    case $opt in
        b)
            BED_FILE_LIST="$OPTARG"
            ;;
        f)
            FASTA_FILE_LIST="$OPTARG"
            ;;
        o)
            OUTPUT_PREFIX="$OPTARG"
            ;;
        c)
            CLEANUP=true
            ;;
        v)
            VERBOSE=true
            ;;
        q)
            QUICK_MODE=true
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
if [[ -z "$BED_FILE_LIST" ]]; then
    echo "Error: -b (BED file list) is required" >&2
    usage
    exit 1
fi

# Validate inputs
if [[ ! -f "$BED_FILE_LIST" ]]; then
    echo "Error: BED file list '$BED_FILE_LIST' does not exist" >&2
    exit 1
fi

if [[ -n "$FASTA_FILE_LIST" && ! -f "$FASTA_FILE_LIST" ]]; then
    echo "Error: FASTA file list '$FASTA_FILE_LIST' does not exist" >&2
    exit 1
fi

# Determine file type and mode
if [[ -n "$FASTA_FILE_LIST" ]]; then
    FILE_TYPE="FASTA"
    HAMMOCK_FILE_LIST="$FASTA_FILE_LIST"
    echo "Mode: FASTA files with BED reference"
else
    FILE_TYPE="BED"
    HAMMOCK_FILE_LIST="$BED_FILE_LIST"
    echo "Mode: BED files only"
fi

# Set parameters based on file type
set_parameters "$FILE_TYPE" "$QUICK_MODE"

# Set up output paths
OUTPUT_DIR="${OUTPUT_PREFIX}_results"
BEDTOOLS_OUTPUT="${OUTPUT_PREFIX}_bedtools_ref.tsv"
RESULTS_TABLE="${OUTPUT_PREFIX}_results.tsv"
ERROR_LOG="${OUTPUT_PREFIX}_error.log"

# Validate that all files in the lists exist
echo "Validating BED file list..."
while IFS= read -r file_path; do
    # Skip empty lines and comments
    if [[ -n "$file_path" && ! "$file_path" =~ ^[[:space:]]*# ]]; then
        # Remove leading/trailing whitespace
        file_path=$(echo "$file_path" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        
        if [[ ! -f "$file_path" ]]; then
            echo "Error: File '$file_path' from BED list does not exist" >&2
            exit 1
        fi
    fi
done < "$BED_FILE_LIST"

# Validate FASTA file list if provided
if [[ -n "$FASTA_FILE_LIST" ]]; then
    echo "Validating FASTA file list..."
    while IFS= read -r file_path; do
        # Skip empty lines and comments
        if [[ -n "$file_path" && ! "$file_path" =~ ^[[:space:]]*# ]]; then
            # Remove leading/trailing whitespace
            file_path=$(echo "$file_path" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            
            if [[ ! -f "$file_path" ]]; then
                echo "Error: File '$file_path' from FASTA list does not exist" >&2
                exit 1
            fi
        fi
    done < "$FASTA_FILE_LIST"
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Starting complete parameter sweep..."
echo "BED file list (for bedtools): $BED_FILE_LIST"
echo "FASTA file list (for hammock): ${FASTA_FILE_LIST:-$BED_FILE_LIST}"
echo "File type: $FILE_TYPE"
echo "Hammock mode: $HAMMOCK_MODE"
echo "Quick mode: $QUICK_MODE"
echo "Output directory: $OUTPUT_DIR"
echo "Bedtools reference: $BEDTOOLS_OUTPUT"
echo "Results table: $RESULTS_TABLE"
echo "Error log: $ERROR_LOG"
echo "Cleanup intermediate files: $CLEANUP"
echo "Verbose: $VERBOSE"
echo ""

# Initialize error log (will only be created if errors occur)
# Note: Error log will remain empty if no errors occur

# Step 1: Generate bedtools reference
echo "Step 1: Generating bedtools reference..."
if [[ -f "$BEDTOOLS_OUTPUT" ]]; then
    echo "  ✓ Bedtools reference already exists: $BEDTOOLS_OUTPUT"
    echo "  Skipping bedtools generation (delete file to regenerate)"
else
    # For both BED and FASTA files, use bedtools jaccard on BED files
    echo "  Running bedtools jaccard on BED files..."
    
    # Create temporary script for bedtools pairwise comparison
    temp_script=$(mktemp)
    cat > "$temp_script" << 'EOF'
#!/bin/bash
file_list="$1"
output_file="$2"
verbose="$3"

set -e  # Exit on any error

# Read files into array
files=()
while IFS= read -r file_path; do
    if [[ -n "$file_path" && ! "$file_path" =~ ^[[:space:]]*# ]]; then
        file_path=$(echo "$file_path" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        files+=("$file_path")
    fi
done < "$file_list"

echo "  Found ${#files[@]} files for pairwise comparison"

# Write header matching bedtools jaccard output format
echo -e "file1\tfile2\tintersection\tunion\tjaccard\tn_intersections" > "$output_file"

# Compare all pairs
total_pairs=$(( ${#files[@]} * ${#files[@]} ))
current_pair=0

# Create temporary output file to collect all results
temp_output=$(mktemp)

for file1 in "${files[@]}"; do
    for file2 in "${files[@]}"; do
        current_pair=$((current_pair + 1))
        if [[ "$verbose" == "true" ]]; then
            echo "    [$current_pair/$total_pairs] Comparing $(basename "$file1") vs $(basename "$file2")"
        fi
        
        # Create sorted versions of the files for bedtools
        sorted_file1=$(mktemp --suffix=.bed)
        sorted_file2=$(mktemp --suffix=.bed)
        
        # Sort the files (handle both .bed and .bed.gz)
        if [[ "$file1" == *.gz ]]; then
            if ! zcat "$file1" | sort -k1,1 -k2,2n > "$sorted_file1"; then
                echo "Error: Failed to sort $file1" >&2
                rm -f "$sorted_file1" "$sorted_file2"
                exit 1
            fi
        else
            if ! sort -k1,1 -k2,2n "$file1" > "$sorted_file1"; then
                echo "Error: Failed to sort $file1" >&2
                rm -f "$sorted_file1" "$sorted_file2"
                exit 1
            fi
        fi
        
        if [[ "$file2" == *.gz ]]; then
            if ! zcat "$file2" | sort -k1,1 -k2,2n > "$sorted_file2"; then
                echo "Error: Failed to sort $file2" >&2
                rm -f "$sorted_file1" "$sorted_file2"
                exit 1
            fi
        else
            if ! sort -k1,1 -k2,2n "$file2" > "$sorted_file2"; then
                echo "Error: Failed to sort $file2" >&2
                rm -f "$sorted_file1" "$sorted_file2"
                exit 1
            fi
        fi
        
        # Run bedtools jaccard on sorted files and collect output
        if bedtools_result=$(bedtools jaccard -a "$sorted_file1" -b "$sorted_file2" 2>/dev/null); then
            # Extract the data line (skip header) and add file names
            data_line=$(echo "$bedtools_result" | tail -n +2)
            echo -e "$(basename "$file1")\t$(basename "$file2")\t$data_line" >> "$temp_output"
        else
            echo "Error: bedtools jaccard failed for $file1 vs $file2" >&2
            rm -f "$sorted_file1" "$sorted_file2"
            exit 1
        fi
        
        # Clean up temporary sorted files
        rm -f "$sorted_file1" "$sorted_file2"
    done
done

# Append all results to final output file
cat "$temp_output" >> "$output_file"
rm -f "$temp_output"

echo "  ✓ Generated $(( total_pairs )) pairwise comparisons"
EOF
    
    chmod +x "$temp_script"
    bedtools_error_output=$(mktemp)
    if ! "$temp_script" "$BED_FILE_LIST" "$BEDTOOLS_OUTPUT" "$VERBOSE" 2>"$bedtools_error_output"; then
        echo "ERROR: Failed to generate bedtools reference" >&2
        # Only create/write to error log if there was an actual error
        {
            echo "Complete Parameter Sweep Error Log - $(date)"
            echo "=========================================="
            echo ""
            echo "ERROR: Failed to generate bedtools reference - $(date)"
            echo "Bedtools error output:"
            cat "$bedtools_error_output"
            echo ""
        } > "$ERROR_LOG"
        rm -f "$bedtools_error_output"
        exit 1
    fi
    rm -f "$bedtools_error_output"
    rm -f "$temp_script"
fi

if [[ ! -f "$BEDTOOLS_OUTPUT" ]]; then
    echo "Error: Failed to generate bedtools reference" >&2
    exit 1
fi

echo "  ✓ Bedtools reference generated: $BEDTOOLS_OUTPUT"
echo ""

# Step 2: Run parameter sweep
echo "Step 2: Running hammock parameter sweep..."

# Create parameter sweep subdirectory
SWEEP_DIR="$OUTPUT_DIR/parameter_sweep"
mkdir -p "$SWEEP_DIR"

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Build parameter sweep command
PARAM_SWEEP_CMD="$SCRIPT_DIR/parameter_sweep.sh -b $BEDTOOLS_OUTPUT -f $HAMMOCK_FILE_LIST -o $SWEEP_DIR -r $RESULTS_TABLE"
if [[ "$QUICK_MODE" == "true" ]]; then
    PARAM_SWEEP_CMD="$PARAM_SWEEP_CMD --quick"
fi

# Run parameter sweep with real-time progress output
# We'll filter the output to show progress messages while capturing errors
param_sweep_stderr=$(mktemp)

# Use a named pipe to filter stdout in real-time
progress_pipe=$(mktemp -u)
mkfifo "$progress_pipe"

# Start a background process to filter and display progress messages
{
    while IFS= read -r line; do
        # Show progress lines (lines that contain "Testing" or show completion status)
        if [[ "$line" =~ ^[[:space:]]*\[[0-9]+/[0-9]+\][[:space:]]*Testing ]] || \
           [[ "$line" =~ ^[[:space:]]*✓[[:space:]] ]] || \
           [[ "$line" =~ ^[[:space:]]*Testing[[:space:]] ]] || \
           [[ "$line" =~ parameter[[:space:]]combinations ]] || \
           [[ "$line" =~ Parameter[[:space:]]sweep[[:space:]]completed ]]; then
            echo "  $line"
        elif [[ "$VERBOSE" == "true" ]]; then
            # In verbose mode, show all non-error output
            echo "  $line"
        fi
    done < "$progress_pipe"
} &
filter_pid=$!

# Run the parameter sweep command
if ! $PARAM_SWEEP_CMD > "$progress_pipe" 2>"$param_sweep_stderr"; then
    # Kill the filter process
    kill $filter_pid 2>/dev/null
    wait $filter_pid 2>/dev/null
    rm -f "$progress_pipe"
    
    echo "ERROR: Parameter sweep failed" >&2
    # Only create/write to error log if there was an actual error
    {
        if [[ ! -f "$ERROR_LOG" ]]; then
            echo "Complete Parameter Sweep Error Log - $(date)"
            echo "=========================================="
            echo ""
        fi
        echo "ERROR: Parameter sweep failed - $(date)"
        echo "Running parameter sweep command: $PARAM_SWEEP_CMD"
        echo "Current working directory: $(pwd)"
        echo "SCRIPT_DIR: $SCRIPT_DIR"
        echo "Files exist check:"
        echo "  BEDTOOLS_OUTPUT ($BEDTOOLS_OUTPUT): $(test -f "$BEDTOOLS_OUTPUT" && echo "EXISTS" || echo "MISSING")"
        echo "  HAMMOCK_FILE_LIST ($HAMMOCK_FILE_LIST): $(test -f "$HAMMOCK_FILE_LIST" && echo "EXISTS" || echo "MISSING")"
        echo "  parameter_sweep.sh ($SCRIPT_DIR/parameter_sweep.sh): $(test -f "$SCRIPT_DIR/parameter_sweep.sh" && echo "EXISTS" || echo "MISSING")"
        echo ""
        echo "Parameter sweep stderr:"
        cat "$param_sweep_stderr"
        echo ""
    } >> "$ERROR_LOG"
    rm -f "$param_sweep_stderr"
    exit 1
fi

# Wait for the filter process to finish
wait $filter_pid 2>/dev/null
rm -f "$progress_pipe" "$param_sweep_stderr"

if [[ ! -f "$RESULTS_TABLE" ]]; then
    echo "Error: Parameter sweep failed to generate results table" >&2
    exit 1
fi

echo "  ✓ Parameter sweep completed: $RESULTS_TABLE"

# Check if clustering results were also generated
CLUSTERING_RESULTS_TABLE="${RESULTS_TABLE%.tsv}_clustering.tsv"
if [[ -f "$CLUSTERING_RESULTS_TABLE" ]]; then
    echo "  ✓ Clustering tree comparison results: $CLUSTERING_RESULTS_TABLE"
fi
echo ""

# Step 3: Cleanup if requested
if [[ "$CLEANUP" == "true" ]]; then
    echo "Step 3: Cleaning up intermediate files..."
    echo "  Removing temporary parameter sweep directory: $SWEEP_DIR"
    rm -rf "$SWEEP_DIR"
    echo "  ✓ Intermediate files cleaned up"
    echo "  ✓ Bedtools reference preserved: $BEDTOOLS_OUTPUT"
    echo ""
fi

echo "Complete parameter sweep finished!"
echo ""
echo "Summary:"
echo "  File type: $FILE_TYPE"
echo "  Hammock mode: $HAMMOCK_MODE"
echo "  Quick mode: $QUICK_MODE"
echo "  Bedtools reference: $BEDTOOLS_OUTPUT"
echo "  Similarity matrix results: $RESULTS_TABLE"
if [[ -f "$CLUSTERING_RESULTS_TABLE" ]]; then
    echo "  Clustering tree results: $CLUSTERING_RESULTS_TABLE"
fi
echo "  Error log: $ERROR_LOG"
echo "  Output directory: $OUTPUT_DIR"

# Only add completion summary to error log if results table creation failed
if [[ ! -f "$RESULTS_TABLE" ]]; then
    {
        if [[ ! -f "$ERROR_LOG" ]]; then
            echo "Complete Parameter Sweep Error Log - $(date)"
            echo "=========================================="
            echo ""
        fi
        echo "ERROR: Results table was not created - $(date)"
        echo "Results table path: $RESULTS_TABLE"
        echo ""
    } >> "$ERROR_LOG"
elif [[ -f "$RESULTS_TABLE" ]]; then
    result_count=$(tail -n +2 "$RESULTS_TABLE" | wc -l)
    if [[ $result_count -eq 0 ]]; then
        {
            if [[ ! -f "$ERROR_LOG" ]]; then
                echo "Complete Parameter Sweep Error Log - $(date)"
                echo "=========================================="
                echo ""
            fi
            echo "WARNING: No results generated - $(date)"
            echo "Results table: $RESULTS_TABLE"
            echo "Results table exists but is empty - check for parameter sweep issues"
            echo ""
        } >> "$ERROR_LOG"
    fi
fi

# Show a preview of the results
echo ""
echo "Similarity matrix results preview (first 5 lines):"
head -n 6 "$RESULTS_TABLE" | column -t -s $'\t'

# Show clustering results preview if available
if [[ -f "$CLUSTERING_RESULTS_TABLE" ]]; then
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

# Show error log info
echo ""
if [[ -f "$ERROR_LOG" && -s "$ERROR_LOG" ]]; then
    echo "Error log saved to: $ERROR_LOG"
    echo "Error log size: $(wc -l < "$ERROR_LOG") lines"
    echo "To view errors: cat $ERROR_LOG"
else
    echo "Error log: No errors detected (error log is empty or doesn't exist)"
    # Create empty error log file if it doesn't exist, for consistency
    touch "$ERROR_LOG"
fi 