#!/usr/bin/env python3
"""
Filter hammock output to only include rows where both file1 and file2 
have basenames that match entries in the accession list.
"""

import sys
import os
import pandas as pd
from pathlib import Path
import argparse


def load_accession_list(accession_file):
    """
    Load accession list from file.
    
    Args:
        accession_file (str): Path to file containing accessions (one per line)
        
    Returns:
        set: Set of accession basenames (without extensions)
    """
    accessions = set()
    
    with open(accession_file, 'r') as f:
        for line in f:
            accession = line.strip()
            if accession:  # Skip empty lines
                # Remove file extension to get basename
                basename = Path(accession).stem
                accessions.add(basename)
    
    return accessions


def filter_hammock_output(hammock_file, accession_list, output_file=None):
    """
    Filter hammock output to only include rows where both file1 and file2 
    have basenames in the accession list.
    
    Args:
        hammock_file (str): Path to hammock CSV output file
        accession_list (set): Set of accession basenames to keep
        output_file (str): Output file path (if None, prints to stdout)
        
    Returns:
        int: Number of rows in filtered output
    """
    # Read hammock output
    df = pd.read_csv(hammock_file)
    
    print(f"Original hammock output: {len(df)} rows")
    print(f"Accession list contains: {len(accession_list)} accessions")
    
    # Get basenames of file1 and file2
    def get_basename(filename):
        """Extract basename from filename, handling various extensions."""
        name = str(filename)
        # Handle common double extensions
        if name.endswith('.bed.gz'):
            name = name[:-7]  # Remove .bed.gz
        elif name.endswith('.fa.gz'):
            name = name[:-6]  # Remove .fa.gz
        else:
            # Use Path.stem for single extensions
            name = Path(name).stem
        return name
    
    # Apply basename extraction
    df['file1_basename'] = df['file1'].apply(get_basename)
    df['file2_basename'] = df['file2'].apply(get_basename)
    
    # Filter rows where both file1 and file2 basenames are in accession list
    filtered_df = df[
        (df['file1_basename'].isin(accession_list)) & 
        (df['file2_basename'].isin(accession_list))
    ]
    
    # Remove the temporary basename columns
    filtered_df = filtered_df.drop(['file1_basename', 'file2_basename'], axis=1)
    
    print(f"Filtered output: {len(filtered_df)} rows")
    print(f"Kept {len(filtered_df) / len(df) * 100:.1f}% of original rows")
    
    # Get unique files in filtered output
    unique_files = set()
    for _, row in filtered_df.iterrows():
        unique_files.add(get_basename(row['file1']))
        unique_files.add(get_basename(row['file2']))
    
    print(f"Unique files in filtered output: {len(unique_files)}")
    
    # Write output
    if output_file:
        filtered_df.to_csv(output_file, index=False)
        print(f"Filtered output saved to: {output_file}")
    else:
        # Print to stdout
        filtered_df.to_csv(sys.stdout, index=False)
    
    return len(filtered_df)


def main():
    """
    Main function to filter hammock output.
    """
    parser = argparse.ArgumentParser(
        description='Filter hammock output to only include rows where both file1 and file2 have basenames in the accession list',
        epilog='Example: python filter_hammock_output.py hammock_output.csv accession_list.txt filtered_output.csv'
    )
    parser.add_argument('hammock_file', help='Path to hammock CSV output file')
    parser.add_argument('accession_list', help='Path to file containing accessions (one per line)')
    parser.add_argument('output_file', nargs='?', help='Output file path (if not provided, prints to stdout)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not os.path.exists(args.hammock_file):
        print(f"Error: Hammock file '{args.hammock_file}' not found")
        sys.exit(1)
    
    if not os.path.exists(args.accession_list):
        print(f"Error: Accession list file '{args.accession_list}' not found")
        sys.exit(1)
    
    if args.verbose:
        print(f"Filtering hammock output: {args.hammock_file}")
        print(f"Using accession list: {args.accession_list}")
        if args.output_file:
            print(f"Output file: {args.output_file}")
        else:
            print("Output: stdout")
        print()
    
    # Load accession list
    accession_set = load_accession_list(args.accession_list)
    
    if args.verbose:
        print(f"Loaded {len(accession_set)} accessions from list")
        print("Sample accessions:", list(accession_set)[:5])
        print()
    
    # Filter hammock output
    try:
        num_rows = filter_hammock_output(args.hammock_file, accession_set, args.output_file)
        
        if args.verbose:
            print(f"\nFiltering completed successfully!")
            print(f"Final output contains {num_rows} rows")
            
    except Exception as e:
        print(f"Error during filtering: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 