#!/bin/bash

# Script to generate all ENCODE graphs (heatmaps and PCA plots)
# Usage: ./encode_graphs.sh <encode.report> <hammock.results> [file_prefix]

# Check arguments
if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "Usage: $0 <encode.report> <hammock.results> [file_prefix]"
    echo "  encode.report: ENCODE metadata report file"
    echo "  hammock.results: Hammock similarity results file"
    echo "  file_prefix: Optional prefix for output files (default: encode)"
    exit 1
fi

# Get arguments
encode_report="$1"
hammock_results="$2"
file_prefix="${3:-encode}"

# Check if input files exist
if [ ! -f "$encode_report" ]; then
    echo "Error: ENCODE report file '$encode_report' not found"
    exit 1
fi

if [ ! -f "$hammock_results" ]; then
    echo "Error: Hammock results file '$hammock_results' not found"
    exit 1
fi

# Get the directory where this script is located
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Generating ENCODE graphs..."
echo "ENCODE report: $encode_report"
echo "Hammock results: $hammock_results"
echo "File prefix: $file_prefix"
echo "Script directory: $script_dir"
echo ""

# Generate heatmaps and global PCA
echo "Step 1: Generating heatmaps and global PCA plot..."
if Rscript "$script_dir/encode_heatmap.R" "$encode_report" "$hammock_results" "$file_prefix"; then
    echo "✓ Heatmaps and global PCA generated successfully"
else
    echo "✗ Error generating heatmaps and global PCA"
    exit 1
fi

echo ""

# Generate organism-specific PCA plots
echo "Step 2: Generating organism-specific PCA plots..."
if Rscript "$script_dir/encode_PCA.R" "$encode_report" "$hammock_results" "$file_prefix"; then
    echo "✓ Organism-specific PCA plots generated successfully"
else
    echo "✗ Error generating organism-specific PCA plots"
    exit 1
fi

echo ""
echo "All graphs generated successfully!"
echo ""
echo "Generated files:"
echo "  ${file_prefix}_biosample.pdf - Biosample heatmap"
echo "  ${file_prefix}_organism.pdf - Organism heatmap"
echo "  ${file_prefix}_target.pdf - Target heatmap"
echo "  ${file_prefix}_pca.pdf - Global PCA plot (all organisms)"
echo "  ${file_prefix}_organism_specific_pca.pdf - Organism-specific PCA plots" 