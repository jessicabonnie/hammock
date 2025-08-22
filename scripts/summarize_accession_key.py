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
        tuple: (summary_data, has_life_stage) where summary_data is organized by organism, tissue, and optionally life_stage
    """
    # Read the TSV file
    df = pd.read_csv(file_path, sep='\t')
    
    # Check if Life_stage column exists
    has_life_stage = 'Life_stage' in df.columns
    
    # Create summary dictionary
    if has_life_stage:
        # Three-level dictionary: organism -> tissue -> life_stage -> count
        summary = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        
        # Count samples by organism, tissue, and life stage
        for _, row in df.iterrows():
            organism = row['Organism']
            tissue = row['Biosample_term_name']
            life_stage = row['Life_stage'] if pd.notna(row['Life_stage']) and row['Life_stage'].strip() else 'unknown'
            summary[organism][tissue][life_stage] += 1
    else:
        # Two-level dictionary: organism -> tissue -> count (backward compatibility)
        summary = defaultdict(lambda: defaultdict(int))
        
        # Count samples by organism and tissue
        for _, row in df.iterrows():
            organism = row['Organism']
            tissue = row['Biosample_term_name']
            summary[organism][tissue] += 1
    
    return summary, has_life_stage

def print_summary_table(summary, has_life_stage=False):
    """
    Print a formatted summary table.
    
    Args:
        summary (dict): Summary data from analyze_tissue_summary
        has_life_stage (bool): Whether the data includes life stage information
    """
    print("=" * 100)
    if has_life_stage:
        print("TISSUE AND SAMPLE COUNT SUMMARY BY ORGANISM (WITH LIFE STAGE)")
    else:
        print("TISSUE AND SAMPLE COUNT SUMMARY BY ORGANISM")
    print("=" * 100)
    
    for organism in sorted(summary.keys()):
        print(f"\n{organism}:")
        print("-" * 60)
        
        if has_life_stage:
            # Calculate totals for tissues across all life stages
            tissue_totals = defaultdict(int)
            life_stage_totals = defaultdict(int)
            total_samples = 0
            
            for tissue, life_stage_data in summary[organism].items():
                for life_stage, count in life_stage_data.items():
                    tissue_totals[tissue] += count
                    life_stage_totals[life_stage] += count
                    total_samples += count
            
            print(f"Total samples: {total_samples}")
            print(f"Unique tissues: {len(tissue_totals)}")
            print(f"Life stages: {', '.join(sorted(life_stage_totals.keys()))}")
            print()
            
            # Sort tissues by total count (descending)
            tissues = sorted(tissue_totals.items(), key=lambda x: x[1], reverse=True)
            
            print(f"{'Tissue':<35} {'Total':<8} {'%':<6} {'Life Stage Breakdown'}")
            print("-" * 90)
            
            for tissue, total_count in tissues:
                percentage = (total_count / total_samples) * 100
                
                # Get life stage breakdown for this tissue
                life_stages = summary[organism][tissue]
                life_stage_str = ", ".join([f"{ls}:{count}" for ls, count in sorted(life_stages.items())])
                
                print(f"{tissue:<35} {total_count:<8} {percentage:>5.1f}% {life_stage_str}")
                
        else:
            # Original format without life stage
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
    
    print("\n" + "=" * 100)

def create_markdown_table(summary, has_life_stage=False):
    """
    Create a markdown formatted table.
    
    Args:
        summary (dict): Summary data from analyze_tissue_summary
        has_life_stage (bool): Whether the data includes life stage information
        
    Returns:
        str: Markdown formatted table
    """
    if has_life_stage:
        markdown = "# Tissue and Sample Count Summary by Organism (with Life Stage)\n\n"
    else:
        markdown = "# Tissue and Sample Count Summary by Organism\n\n"
    
    for organism in sorted(summary.keys()):
        markdown += f"## {organism}\n\n"
        
        if has_life_stage:
            # Calculate totals for tissues across all life stages
            tissue_totals = defaultdict(int)
            life_stage_totals = defaultdict(int)
            total_samples = 0
            
            for tissue, life_stage_data in summary[organism].items():
                for life_stage, count in life_stage_data.items():
                    tissue_totals[tissue] += count
                    life_stage_totals[life_stage] += count
                    total_samples += count
            
            markdown += f"**Total samples:** {total_samples}  \n"
            markdown += f"**Unique tissues:** {len(tissue_totals)}  \n"
            markdown += f"**Life stages:** {', '.join(sorted(life_stage_totals.keys()))}  \n"
            
            # Add life stage summary
            markdown += f"**Life stage distribution:** "
            life_stage_summary = []
            for life_stage in sorted(life_stage_totals.keys()):
                count = life_stage_totals[life_stage]
                pct = (count / total_samples) * 100
                life_stage_summary.append(f"{life_stage} ({count}, {pct:.1f}%)")
            markdown += ", ".join(life_stage_summary) + "\n\n"
            
            # Sort tissues by total count (descending)
            tissues = sorted(tissue_totals.items(), key=lambda x: x[1], reverse=True)
            
            markdown += "| Tissue | Total Count | Percentage | Life Stage Breakdown |\n"
            markdown += "|--------|-------------|------------|----------------------|\n"
            
            for tissue, total_count in tissues:
                percentage = (total_count / total_samples) * 100
                
                # Get life stage breakdown for this tissue
                life_stages = summary[organism][tissue]
                life_stage_breakdown = []
                for ls in sorted(life_stages.keys()):
                    count = life_stages[ls]
                    life_stage_breakdown.append(f"{ls}: {count}")
                life_stage_str = ", ".join(life_stage_breakdown)
                
                markdown += f"| {tissue} | {total_count} | {percentage:.1f}% | {life_stage_str} |\n"
        
        else:
            # Original format without life stage
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
        summary, has_life_stage = analyze_tissue_summary(args.file_path)
        
        # Print summary table
        print_summary_table(summary, has_life_stage)
        
        # Create and save markdown table if requested
        if not args.no_markdown:
            markdown_table = create_markdown_table(summary, has_life_stage)
            
            with open(args.output, "w") as f:
                f.write(markdown_table)
            
            print(f"\nMarkdown table saved to: {args.output}")
            if has_life_stage:
                print("Note: Life stage information was detected and included in the analysis.")
            else:
                print("Note: No Life_stage column found - analysis performed without life stage information.")
        
    except FileNotFoundError:
        print(f"Error: File '{args.file_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
