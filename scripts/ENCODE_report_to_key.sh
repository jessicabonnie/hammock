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
echo -e "Accession\tFile\tTarget_of_Assay\tBiosample_term_name\tOrganism\tLife_stage" > "$OUTPUT_FILE"

# Process the TSV file using awk for more reliable field parsing
# Skip only the first line (timestamp); read the header row to map column names dynamically
tail -n +2 "$INPUT_FILE" | awk -F'\t' -v output_file="$OUTPUT_FILE" -v input_file="$INPUT_FILE" '
NR == 1 {
    # Map header names to indices (1-based)
    for (i = 1; i <= NF; i++) {
        header[$i] = i
    }

    # Resolve required/optional column indices
    acc_i     = header["Accession"]
    target_i  = header["Target of assay"]        # optional in some reports
    bios_i    = header["Biosample term name"]
    files_i   = header["Files"]
    org_i     = header["Organism"]
    stage_i   = header["Life stage"]              # optional in some reports

    # Validate required columns exist
    if (!acc_i || !files_i || !bios_i || !org_i) {
        msg = "Error: One or more required columns are missing in the header. Required: Accession, Biosample term name, Files, Organism. Optional: Target of assay, Life stage"
        print msg > "/dev/stderr"
        exit 1
    }

    # Warn for optional columns missing
    if (!target_i) {
        print "Warning: Optional column Target of assay not found in header for " input_file > "/dev/stderr"
    }
    if (!stage_i) {
        print "Warning: Optional column Life stage not found in header for " input_file > "/dev/stderr"
    }

    next
}

{
    accession = $(acc_i)
    target_of_assay = (target_i ? $(target_i) : "")
    biosample_term_name = $(bios_i)
    files = $(files_i)
    organism = $(org_i)
    life_stage = (stage_i ? $(stage_i) : "")

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
            print accession "\t" filename "\t" target_of_assay "\t" biosample_term_name "\t" organism "\t" life_stage >> output_file
        }
    }
}'

echo "Long table created successfully: $OUTPUT_FILE"
echo "Total rows (including header): $(wc -l < "$OUTPUT_FILE")"