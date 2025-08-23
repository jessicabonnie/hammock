#!/bin/bash

# Improved bedtools pairwise jaccard script
# Takes a list of BED files and computes pairwise Jaccard similarities

set -e  # Exit on any error

# Default values
OUTPUT_FILE="bedtools_pairwise_jaccard.tsv"
VERBOSE=false

# Function to display usage
usage() {
    echo "Usage: $0 [OPTIONS] <bed_file1> <bed_file2> ... <bed_fileN>"
    echo "   or: $0 [OPTIONS] -f <file_list>"
    echo ""
    echo "Compute pairwise Jaccard similarities between BED files using bedtools."
    echo ""
    echo "Options:"
    echo "  -o <output_file>    Output file (default: bedtools_pairwise_jaccard.txt)"
    echo "  -f <file_list>      File containing list of BED file paths (one per line)"
    echo "  -v                  Verbose output"
    echo "  -h                  Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 file1.bed file2.bed file3.bed"
    echo "  $0 -f bed_files_list.txt -o results.txt"
    echo "  $0 -o custom_output.txt *.bed"
}

# Parse command line arguments
while getopts "o:f:vh" opt; do
    case $opt in
        o)
            OUTPUT_FILE="$OPTARG"
            ;;
        f)
            FILE_LIST="$OPTARG"
            ;;
        v)
            VERBOSE=true
            ;;
        h)
            usage
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
    esac
done

shift $((OPTIND-1))

# Collect BED files
BED_FILES=()

if [[ -n "$FILE_LIST" ]]; then
    # Use files from file list
    if [[ ! -f "$FILE_LIST" ]]; then
        echo "Error: File list '$FILE_LIST' does not exist" >&2
        exit 1
    fi
    
    # Read file paths from the list file
    while IFS= read -r file_path; do
        # Skip empty lines and comments
        if [[ -n "$file_path" && ! "$file_path" =~ ^[[:space:]]*# ]]; then
            # Remove leading/trailing whitespace
            file_path=$(echo "$file_path" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            
            if [[ ! -f "$file_path" ]]; then
                echo "Error: File '$file_path' from list does not exist" >&2
                exit 1
            fi
            BED_FILES+=("$file_path")
        fi
    done < "$FILE_LIST"
    
    if [[ ${#BED_FILES[@]} -eq 0 ]]; then
        echo "Error: No valid BED files found in file list '$FILE_LIST'" >&2
        exit 1
    fi
    
    if [[ "$VERBOSE" == true ]]; then
        echo "Found ${#BED_FILES[@]} BED files in file list '$FILE_LIST'"
    fi
    
else
    # Use files from command line arguments
    if [[ $# -eq 0 ]]; then
        echo "Error: No BED files specified" >&2
        usage
        exit 1
    fi
    
    # Validate that all files exist
    for file in "$@"; do
        if [[ ! -f "$file" ]]; then
            echo "Error: File '$file' does not exist" >&2
            exit 1
        fi
        BED_FILES+=("$file")
    done
fi

# Check if bedtools is available
if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools is not installed or not in PATH" >&2
    exit 1
fi

# Sort files for consistent ordering
IFS=$'\n' BED_FILES=($(sort <<<"${BED_FILES[*]}"))
unset IFS

if [[ "$VERBOSE" == true ]]; then
    echo "Processing ${#BED_FILES[@]} BED files:"
    for file in "${BED_FILES[@]}"; do
        echo "  $(basename "$file")"
    done
    echo "Output file: $OUTPUT_FILE"
    echo ""
fi

# Create output file with header
echo "file1 file2 intersection union jaccard n_intersections" > "$OUTPUT_FILE"

# Calculate total number of comparisons for progress tracking
total_comparisons=$((${#BED_FILES[@]} * ${#BED_FILES[@]}))
current_comparison=0

if [[ "$VERBOSE" == true ]]; then
    echo "Computing $total_comparisons pairwise comparisons..."
    echo ""
fi

# Compute pairwise Jaccard similarities
for file1 in "${BED_FILES[@]}"; do
    for file2 in "${BED_FILES[@]}"; do
        current_comparison=$((current_comparison + 1))
        
        if [[ "$VERBOSE" == true ]]; then
            echo "[$current_comparison/$total_comparisons] Comparing $(basename "$file1") vs $(basename "$file2")"
        fi
        
        # Compute the jaccard statistic for these two files
        bvalues=$(bedtools jaccard -a "$file1" -b "$file2" | awk 'NR==2 {print $0}')
        
        # Extract just the filenames (not full paths) for output
        file1_name=$(basename "$file1")
        file2_name=$(basename "$file2")
        
        # Write results to output file
        echo "$file1_name $file2_name $bvalues" >> "$OUTPUT_FILE"
    done
done

if [[ "$VERBOSE" == true ]]; then
    echo ""
    echo "Pairwise Jaccard computation completed!"
    echo "Results saved to: $OUTPUT_FILE"
    echo "Total comparisons: $total_comparisons"
    
    # Show a preview of results
    echo ""
    echo "Results preview (first 5 lines):"
    head -n 6 "$OUTPUT_FILE"
else
    echo "Results saved to: $OUTPUT_FILE"
fi
