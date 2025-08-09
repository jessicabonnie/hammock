#!/usr/bin/env python3
"""
Script to analyze filtered accession key data and create a summary table
of tissues and sample counts by organism.
"""

import pandas as pd
import argparse
import sys
from collections import defaultdict

def analyze_tissue_summary(file_path):
    """
    Analyze the filtered accession key file and create a summary table.
    
    Args:
        file_path (str): Path to the TSV file
        
    Returns:
        dict: Summary data organized by organism and tissue
    """
    # Read the TSV file
    df = pd.read_csv(file_path, sep='\t')
    
    # Create summary dictionary
    summary = defaultdict(lambda: defaultdict(int))
    
    # Count samples by organism and tissue
    for _, row in df.iterrows():
        organism = row['Organism']
        tissue = row['Biosample_term_name']
        summary[organism][tissue] += 1
    
    return summary

def print_summary_table(summary):
    """
    Print a formatted summary table.
    
    Args:
        summary (dict): Summary data from analyze_tissue_summary
    """
    print("=" * 80)
    print("TISSUE AND SAMPLE COUNT SUMMARY BY ORGANISM")
    print("=" * 80)
    
    for organism in sorted(summary.keys()):
        print(f"\n{organism}:")
        print("-" * 40)
        
        # Sort tissues by count (descending)
        tissues = sorted(summary[organism].items(), key=lambda x: x[1], reverse=True)
        
        total_samples = sum(summary[organism].values())
        print(f"Total samples: {total_samples}")
        print(f"Unique tissues: {len(summary[organism])}")
        print()
        
        print(f"{'Tissue':<30} {'Count':<10} {'Percentage':<10}")
        print("-" * 50)
        
        for tissue, count in tissues:
            percentage = (count / total_samples) * 100
            print(f"{tissue:<30} {count:<10} {percentage:>6.1f}%")
    
    print("\n" + "=" * 80)

def create_markdown_table(summary):
    """
    Create a markdown formatted table.
    
    Args:
        summary (dict): Summary data from analyze_tissue_summary
        
    Returns:
        str: Markdown formatted table
    """
    markdown = "# Tissue and Sample Count Summary by Organism\n\n"
    
    for organism in sorted(summary.keys()):
        markdown += f"## {organism}\n\n"
        
        tissues = sorted(summary[organism].items(), key=lambda x: x[1], reverse=True)
        total_samples = sum(summary[organism].values())
        
        markdown += f"**Total samples:** {total_samples}  \n"
        markdown += f"**Unique tissues:** {len(summary[organism])}\n\n"
        
        markdown += "| Tissue | Count | Percentage |\n"
        markdown += "|--------|-------|------------|\n"
        
        for tissue, count in tissues:
            percentage = (count / total_samples) * 100
            markdown += f"| {tissue} | {count} | {percentage:.1f}% |\n"
        
        markdown += "\n"
    
    return markdown

def main():
    """Main function to run the analysis."""
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Analyze accession key file and create tissue summary table by organism"
    )
    parser.add_argument(
        "file_path", 
        help="Path to the TSV accession key file"
    )
    parser.add_argument(
        "--output", "-o",
        default="tissue_summary_table.md",
        help="Output markdown file path (default: tissue_summary_table.md)"
    )
    parser.add_argument(
        "--no-markdown", 
        action="store_true",
        help="Don't save markdown file, only print to console"
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    try:
        # Analyze the data
        summary = analyze_tissue_summary(args.file_path)
        
        # Print summary table
        print_summary_table(summary)
        
        # Create and save markdown table if requested
        if not args.no_markdown:
            markdown_table = create_markdown_table(summary)
            
            with open(args.output, "w") as f:
                f.write(markdown_table)
            
            print(f"\nMarkdown table saved to: {args.output}")
        
    except FileNotFoundError:
        print(f"Error: File '{args.file_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
