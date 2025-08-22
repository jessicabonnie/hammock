#!/bin/bash

# Script to batch process CSV files with filter_hammock_output.py
# Usage: batch_filter_hammock.sh <input_directory> <accession_list> [tag]

# Function to display usage
usage() {
    echo "Usage: $0 <input_directory> <accession_list> [tag]"
    echo ""
    echo "Arguments:"
    echo "  input_directory  Directory containing CSV files to process"
    echo "  accession_list   Path to accession list file"
    echo "  tag             Optional tag to append to output filenames (before extension)"
    echo ""
    echo "Examples:"
    echo "  $0 /path/to/csvs accessions.txt"
    echo "  $0 /path/to/csvs accessions.txt _filtered"
    echo ""
    echo "Output files are written to the current directory."
    exit 1
}

# Check if at least 2 arguments are provided
if [ $# -lt 2 ]; then
    echo "Error: Insufficient arguments"
    usage
fi

INPUT_DIR="$1"
ACCESSION_LIST="$2"
TAG="${3:-}"  # Optional tag, empty string if not provided

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FILTER_SCRIPT="$SCRIPT_DIR/filter_hammock_output.py"

# Validation
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory '$INPUT_DIR' does not exist"
    exit 1
fi

if [ ! -f "$ACCESSION_LIST" ]; then
    echo "Error: Accession list file '$ACCESSION_LIST' does not exist"
    exit 1
fi

if [ ! -f "$FILTER_SCRIPT" ]; then
    echo "Error: Filter script '$FILTER_SCRIPT' not found"
    echo "Make sure filter_hammock_output.py is in the same directory as this script"
    exit 1
fi

# Check if Python script is executable
if [ ! -x "$FILTER_SCRIPT" ]; then
    echo "Making filter script executable..."
    chmod +x "$FILTER_SCRIPT"
fi

# Count CSV files in input directory
CSV_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -name "*.csv" -type f | wc -l)

if [ "$CSV_COUNT" -eq 0 ]; then
    echo "No CSV files found in directory: $INPUT_DIR"
    exit 1
fi

echo "Found $CSV_COUNT CSV file(s) in $INPUT_DIR"
echo "Using accession list: $ACCESSION_LIST"
if [ -n "$TAG" ]; then
    echo "Output tag: $TAG"
fi
echo "Output directory: $(pwd)"
echo ""

# Process each CSV file
PROCESSED=0
FAILED=0

for csv_file in "$INPUT_DIR"/*.csv; do
    # Skip if no CSV files match (shouldn't happen due to earlier check)
    [ ! -f "$csv_file" ] && continue
    
    # Get the basename without path and extension
    basename=$(basename "$csv_file" .csv)
    
    # Construct output filename
    if [ -n "$TAG" ]; then
        output_file="${basename}${TAG}.csv"
    else
        output_file="${basename}.csv"
    fi
    
    echo "Processing: $(basename "$csv_file") -> $output_file"
    
    # Run the filter script
    if "$FILTER_SCRIPT" "$csv_file" "$ACCESSION_LIST" "$output_file" 2>&1; then
        ((PROCESSED++))
        echo "✓ Successfully processed $(basename "$csv_file")"
    else
        exit_code=$?
        ((FAILED++))
        echo "✗ Failed to process $(basename "$csv_file") (exit code: $exit_code)"
    fi
    echo ""
done

# Summary
echo "=========================="
echo "Processing complete!"
echo "Successfully processed: $PROCESSED files"
echo "Failed: $FAILED files"
echo "Output files written to: $(pwd)"

if [ "$FAILED" -gt 0 ]; then
    exit 1
fi
