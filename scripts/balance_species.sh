#!/bin/bash

# =============================================================================
# ENCODE Species Balancing Script
# =============================================================================
#
# PURPOSE:
# This script creates a filtered accession key that ensures equal representation of 
# Homo sapiens and Mus musculus samples across all tissue categories. The output
# maintains the original TSV format but only includes accessions that will produce
# a balanced dataset for comparative genomics studies.
#
# FUNCTIONALITY:
# 1. Reads a TSV file containing accession keys with tissue and species information
# 2. Identifies all unique tissue categories in the dataset
# 3. For each tissue category, counts samples from both Homo sapiens and Mus musculus
# 4. Determines the minimum count between the two species for each tissue
# 5. Randomly (with set seed) selects that minimum number of samples from each species
# 6. Outputs a filtered accession key with only the selected accessions
#
# INPUT FORMAT:
# The input TSV file should have the following columns:
# Column 1: ENCSR accession (experiment)
# Column 2: ENCFF accession (file)
# Column 3: (empty)
# Column 4: Tissue category (e.g., "adrenal gland", "heart", "brain")
# Column 5: Species ("Homo sapiens" or "Mus musculus")
#
# OUTPUT FORMAT:
# Same as input format, but filtered to include only accessions needed for balancing
#
# USAGE:
# ./balance_species.sh <input_file> <output_file>
# Example: ./balance_species.sh filtered_accession_key.tsv balanced_accession_key.tsv
#
# FEATURES:
# - Random sampling to avoid selection bias
# - Comprehensive progress reporting
# - Detailed summary statistics
# - Handles tissue names with spaces correctly
# - Skips tissues missing one or both species
# - Maintains original file format and structure
# - Outputs filtered accession key for downstream processing
# - Preserves header row from input file
#
# DEPENDENCIES:
# - bash
# - awk
# - shuf (for random sampling)
# - sort
# - wc (for line counting)
#
# AUTHOR: Generated for DNase-seq comparative analysis
# DATE: 2024
# =============================================================================

# Script to create a filtered accession key with balanced species representation
# Usage: ./balance_species.sh input_file output_file [seed]
# 
# The script uses reproducible random sampling. If no seed is provided, 
# a default seed of 42 is used. To use a different seed:
# ./balance_species.sh input.tsv output.tsv 12345

INPUT_FILE="$1"
OUTPUT_FILE="$2"

# Set default seed if not provided
SEED="${3:-42}"

if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "Usage: $0 input_file output_file [seed]"
    echo "Example: $0 filtered_accession_key.tsv balanced_accession_key.tsv"
    echo "Example with seed: $0 filtered_accession_key.tsv balanced_accession_key.tsv 12345"
    exit 1
fi

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found"
    exit 1
fi

echo "Creating filtered accession key with balanced species representation..."
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Random seed: $SEED"

# Create temporary directory for processing
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

# Extract and save the header row
head -1 "$INPUT_FILE" > "$TEMP_DIR/header.txt"

# Get unique tissue categories (skip header row)
tail -n +2 "$INPUT_FILE" | awk -F'\t' '{print $4}' | sort | uniq > "$TEMP_DIR/tissues.txt"

# Process each tissue category
while read tissue; do
    echo "Processing tissue: '$tissue'"
    
    # Count samples for each species for this tissue (skip header)
    homo_count=$(tail -n +2 "$INPUT_FILE" | awk -F'\t' -v tissue="$tissue" '$4 == tissue && $5 == "Homo sapiens"' | wc -l)
    mus_count=$(tail -n +2 "$INPUT_FILE" | awk -F'\t' -v tissue="$tissue" '$4 == tissue && $5 == "Mus musculus"' | wc -l)
    
    echo "  Homo sapiens: $homo_count, Mus musculus: $mus_count"
    
    # Find minimum count
    if [ $homo_count -eq 0 ] || [ $mus_count -eq 0 ]; then
        echo "  Skipping '$tissue' - missing one or both species"
        continue
    fi
    
    min_count=$((homo_count < mus_count ? homo_count : mus_count))
    echo "  Minimum count for '$tissue': $min_count"
    
    # Extract and randomly select samples for each species (skip header)
    # Use reproducible random sampling with the specified seed
    tail -n +2 "$INPUT_FILE" | awk -F'\t' -v tissue="$tissue" '$4 == tissue && $5 == "Homo sapiens"' | \
        shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:"$SEED" -nosalt </dev/zero 2>/dev/null) | head -n "$min_count" >> "$TEMP_DIR/balanced_samples.txt"
    
    tail -n +2 "$INPUT_FILE" | awk -F'\t' -v tissue="$tissue" '$4 == tissue && $5 == "Mus musculus"' | \
        shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:"$SEED" -nosalt </dev/zero 2>/dev/null) | head -n "$min_count" >> "$TEMP_DIR/balanced_samples.txt"
    
    echo "  Selected $min_count samples from each species for '$tissue'"
    
done < "$TEMP_DIR/tissues.txt"

# Check if we have any balanced samples
if [ ! -f "$TEMP_DIR/balanced_samples.txt" ]; then
    echo "Error: No balanced samples found. Creating output file with header only."
    cat "$TEMP_DIR/header.txt" > "$OUTPUT_FILE"
else
    # Combine header with balanced samples and sort
    cat "$TEMP_DIR/header.txt" "$TEMP_DIR/balanced_samples.txt" | \
        awk 'NR==1 {print; next} {print}' | sort -k4,4 -k5,5 > "$OUTPUT_FILE"
fi

# Print summary statistics (skip header in counts)
echo ""
echo "Filtered accession key created successfully!"
echo "Summary of balanced dataset:"
echo "Tissue Category | Homo sapiens | Mus musculus | Total"
echo "----------------|--------------|--------------|-------"

while read tissue; do
    homo_count=$(tail -n +2 "$OUTPUT_FILE" | awk -F'\t' -v tissue="$tissue" '$4 == tissue && $5 == "Homo sapiens"' | wc -l)
    mus_count=$(tail -n +2 "$OUTPUT_FILE" | awk -F'\t' -v tissue="$tissue" '$4 == tissue && $5 == "Mus musculus"' | wc -l)
    total=$((homo_count + mus_count))
    printf "%-15s | %-12s | %-12s | %s\n" "$tissue" "$homo_count" "$mus_count" "$total"
done < "$TEMP_DIR/tissues.txt"

echo ""
echo "Total accessions in filtered key: $(tail -n +2 "$OUTPUT_FILE" | wc -l)"
echo "Original dataset had: $(tail -n +2 "$INPUT_FILE" | wc -l) accessions"
echo ""
echo "This filtered accession key can now be used to create balanced datasets"
echo "for downstream analysis by filtering your bed files and fasta files accordingly." 