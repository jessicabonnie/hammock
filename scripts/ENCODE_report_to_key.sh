#!/bin/bash

# Script to create a long table from ENCODE TSV data
# Converts experiment report with comma-separated files into individual file entries
# Usage: ./ENCODE_report_to_key.sh input_file.tsv output_file.tsv

# Check if correct number of arguments provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 input_file.tsv output_file.tsv"
    echo "Example: $0 experiment_report_2025_6_11_18h_49m.tsv long_table.tsv"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found!"
    exit 1
fi

# Write header to output file
echo -e "Accession\tFile\tTarget_of_Assay\tBiosample_term_name\tOrganism" > "$OUTPUT_FILE"

# Process the TSV file using awk for more reliable field parsing
# Skip the first line (timestamp) and header line, then process each data row
tail -n +3 "$INPUT_FILE" | awk -F'\t' -v output_file="$OUTPUT_FILE" '
{
    # Extract relevant columns (1-based indexing in awk)
    # Accession=2, Target_of_assay=7, Biosample_term_name=10, Files=16, Organism=22
    accession = $2
    target_of_assay = $7
    biosample_term_name = $10
    files = $16
    organism = $22
    
    # Skip empty lines or lines with missing data
    if (accession == "" || files == "") {
        next
    }
    
    # Split files by comma and process each file
    split(files, file_array, ",")
    for (i in file_array) {
        # Remove leading/trailing whitespace
        gsub(/^[ \t]+|[ \t]+$/, "", file_array[i])
        
        # Extract just the filename from the path (remove /files/ prefix and trailing /)
        filename = file_array[i]
        gsub(/^\/files\//, "", filename)
        gsub(/\/$/, "", filename)
        
        # Skip empty filenames
        if (filename != "") {
            # Write to output file
            print accession "\t" filename "\t" target_of_assay "\t" biosample_term_name "\t" organism >> output_file
        }
    }
}'

echo "Long table created successfully: $OUTPUT_FILE"
echo "Total rows (including header): $(wc -l < "$OUTPUT_FILE")"