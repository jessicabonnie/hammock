#!/usr/bin/env python3
"""
Script to compare similarity matrices from bedtools and hammock outputs.
Detects file format automatically and calculates Frobenius norm between matrices.
"""

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path


def detect_file_format(filepath):
    """
    Detect if the file is in bedtools or hammock format.
    
    Args:
        filepath (str): Path to the input file
        
    Returns:
        str: Either 'bedtools' or 'hammock'
    """
    try:
        # Read first few lines to detect format
        with open(filepath, 'r') as f:
            first_line = f.readline().strip()
            
        # Check for bedtools format (has specific column names, tab or space separated)
        if (first_line.startswith('file1\tfile2\tintersection\tunion\tjaccard') or 
            first_line.startswith('file1 file2 intersection union jaccard')):
            return 'bedtools'
        elif first_line.startswith('file1,file2') and ('jaccard_similarity' in first_line or 'jaccard_similarity_with_ends' in first_line):
            return 'hammock'
        else:
            # Try to determine by content
            df = pd.read_csv(filepath, nrows=5)
            if 'jaccard' in df.columns and 'intersection' in df.columns:
                return 'bedtools'
            elif 'jaccard_similarity' in df.columns or 'jaccard_similarity_with_ends' in df.columns:
                return 'hammock'
            else:
                raise ValueError(f"Cannot determine file format for {filepath}")
                
    except Exception as e:
        raise ValueError(f"Error detecting file format for {filepath}: {e}")


def parse_bedtools_format(filepath):
    """
    Parse bedtools pairwise jaccard output into a similarity matrix.
    
    Args:
        filepath (str): Path to bedtools output file
        
    Returns:
        pandas.DataFrame: Symmetric similarity matrix
    """
    # Read the bedtools output (mixed space/tab separated - use whitespace)
    df = pd.read_csv(filepath, sep=r'\s+')
    
    # Extract similarity values and file pairs
    similarity_data = df[['file1', 'file2', 'jaccard']].copy()
    
    # Convert file names to basenames without extensions for matching
    # Handle double extensions like .bed.gz
    def get_basename(filename):
        # Remove common double extensions first
        name = str(filename)
        if name.endswith('.bed.gz'):
            name = name[:-7]  # Remove .bed.gz
        elif name.endswith('.fa.gz'):
            name = name[:-6]  # Remove .fa.gz
        else:
            # Use Path.stem for single extensions
            name = Path(name).stem
        return name
    
    similarity_data['file1_base'] = similarity_data['file1'].apply(get_basename)
    similarity_data['file2_base'] = similarity_data['file2'].apply(get_basename)
    
    # Get unique file basenames
    all_files = sorted(set(similarity_data['file1_base'].tolist() + similarity_data['file2_base'].tolist()))
    
    # Create empty similarity matrix
    n_files = len(all_files)
    sim_matrix = np.zeros((n_files, n_files))
    
    # Create file name to index mapping
    file_to_idx = {file: idx for idx, file in enumerate(all_files)}
    
    # Fill the similarity matrix
    for _, row in similarity_data.iterrows():
        i = file_to_idx[row['file1_base']]
        j = file_to_idx[row['file2_base']]
        sim_matrix[i, j] = row['jaccard']
        # Make symmetric (though bedtools should already include both directions)
        sim_matrix[j, i] = row['jaccard']
    
    # Convert to DataFrame with proper labels (using basenames)
    sim_df = pd.DataFrame(sim_matrix, index=all_files, columns=all_files)
    
    return sim_df


def parse_hammock_format(filepath):
    """
    Parse hammock CSV output into a similarity matrix.
    
    Args:
        filepath (str): Path to hammock CSV file
        
    Returns:
        pandas.DataFrame: Symmetric similarity matrix
    """
    # Read the hammock CSV output
    df = pd.read_csv(filepath)
    
    # Determine which jaccard column to use
    if 'jaccard_similarity_with_ends' in df.columns:
        # Mode D (FASTA files) - use jaccard_similarity_with_ends
        jaccard_col = 'jaccard_similarity_with_ends'
    elif 'jaccard_similarity' in df.columns:
        # Mode B (BED files) - use jaccard_similarity
        jaccard_col = 'jaccard_similarity'
    else:
        raise ValueError("No jaccard similarity column found in hammock output")
    
    # Extract similarity values and file pairs
    similarity_data = df[['file1', 'file2', jaccard_col]].copy()
    
    # Convert file names to basenames without extensions for matching
    # Handle double extensions like .bed.gz
    def get_basename(filename):
        # Remove common double extensions first
        name = str(filename)
        if name.endswith('.bed.gz'):
            name = name[:-7]  # Remove .bed.gz
        elif name.endswith('.fa.gz'):
            name = name[:-6]  # Remove .fa.gz
        else:
            # Use Path.stem for single extensions
            name = Path(name).stem
        return name
    
    similarity_data['file1_base'] = similarity_data['file1'].apply(get_basename)
    similarity_data['file2_base'] = similarity_data['file2'].apply(get_basename)
    
    # Get unique file basenames
    all_files = sorted(set(similarity_data['file1_base'].tolist() + similarity_data['file2_base'].tolist()))
    
    # Create empty similarity matrix
    n_files = len(all_files)
    sim_matrix = np.zeros((n_files, n_files))
    
    # Create file name to index mapping
    file_to_idx = {file: idx for idx, file in enumerate(all_files)}
    
    # Fill the similarity matrix
    for _, row in similarity_data.iterrows():
        i = file_to_idx[row['file1_base']]
        j = file_to_idx[row['file2_base']]
        sim_matrix[i, j] = row[jaccard_col]
        # Make symmetric
        sim_matrix[j, i] = row[jaccard_col]
    
    # Convert to DataFrame with proper labels (using basenames)
    sim_df = pd.DataFrame(sim_matrix, index=all_files, columns=all_files)
    
    return sim_df


def align_matrices(matrix1, matrix2, quiet=False):
    """
    Align two similarity matrices to have the same dimensions and ordering.
    Only keeps files that are present in both matrices.
    
    Args:
        matrix1 (pandas.DataFrame): First similarity matrix
        matrix2 (pandas.DataFrame): Second similarity matrix
        quiet (bool): If True, suppress status messages
        
    Returns:
        tuple: (aligned_matrix1, aligned_matrix2)
    """
    # Find common files
    common_files = sorted(set(matrix1.index).intersection(set(matrix2.index)))
    
    if len(common_files) == 0:
        raise ValueError("No common files found between the two matrices")
    
    if not quiet:
        print(f"Found {len(common_files)} common files between matrices")
        print(f"Matrix 1 has {len(matrix1)} files, Matrix 2 has {len(matrix2)} files")
    
    # Subset both matrices to common files and same order
    aligned_matrix1 = matrix1.loc[common_files, common_files]
    aligned_matrix2 = matrix2.loc[common_files, common_files]
    
    return aligned_matrix1, aligned_matrix2


def calculate_frobenius_norm(matrix1, matrix2):
    """
    Calculate the Frobenius norm between two matrices.
    
    Args:
        matrix1 (pandas.DataFrame): First matrix
        matrix2 (pandas.DataFrame): Second matrix
        
    Returns:
        float: Frobenius norm of the difference
    """
    # Ensure matrices have same shape
    if matrix1.shape != matrix2.shape:
        raise ValueError(f"Matrix shapes don't match: {matrix1.shape} vs {matrix2.shape}")
    
    # Calculate element-wise difference
    diff_matrix = matrix1.values - matrix2.values
    
    # Calculate Frobenius norm
    frobenius_norm = np.linalg.norm(diff_matrix, 'fro')
    
    return frobenius_norm


def calculate_additional_metrics(matrix1, matrix2):
    """
    Calculate additional comparison metrics between two matrices.
    
    Args:
        matrix1 (pandas.DataFrame): First matrix
        matrix2 (pandas.DataFrame): Second matrix
        
    Returns:
        dict: Dictionary of metrics
    """
    # Convert to numpy arrays
    m1 = matrix1.values
    m2 = matrix2.values
    
    # Flatten for correlation calculation (excluding diagonal)
    mask = ~np.eye(m1.shape[0], dtype=bool)
    m1_flat = m1[mask]
    m2_flat = m2[mask]
    
    # Calculate metrics
    correlation = np.corrcoef(m1_flat, m2_flat)[0, 1]
    mean_abs_error = np.mean(np.abs(m1 - m2))
    max_abs_error = np.max(np.abs(m1 - m2))
    
    metrics = {
        'correlation': correlation,
        'mean_absolute_error': mean_abs_error,
        'max_absolute_error': max_abs_error,
        'matrix_size': m1.shape[0]
    }
    
    return metrics


def print_table_header():
    """
    Print the header for tab-delimited output.
    """
    print("file1\tfile2\tformat1\tformat2\tmatrix_size\tfrobenius_norm\tcorrelation\tmean_abs_error\tmax_abs_error")


def main():
    """
    Main function to compare similarity matrices from two files.
    """
    parser = argparse.ArgumentParser(
        description='Compare similarity matrices from bedtools and hammock outputs',
        epilog='Example: python compare_sim_matrices.py file1.txt file2.csv --table'
    )
    parser.add_argument('file1', nargs='?', help='First input file (bedtools or hammock format)')
    parser.add_argument('file2', nargs='?', help='Second input file (bedtools or hammock format)')
    parser.add_argument('--output', '-o', help='Output file for detailed comparison (optional)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    parser.add_argument('--table', '-t', action='store_true', 
                       help='Output results as tab-delimited line (use --header to print column names)')
    parser.add_argument('--header', action='store_true', 
                       help='Print tab-delimited header (use with --table)')
    
    args = parser.parse_args()
    
    # Print header if requested
    if args.header:
        print_table_header()
        if not args.table:
            sys.exit(0)  # Exit if only header was requested
    
    # Check if files are provided when not just printing header
    if not args.file1 or not args.file2:
        if not args.header:
            parser.error("file1 and file2 are required unless using --header")
        else:
            sys.exit(0)
    
    # Check if files exist
    if not Path(args.file1).exists():
        print(f"Error: File {args.file1} does not exist", file=sys.stderr)
        sys.exit(1)
    if not Path(args.file2).exists():
        print(f"Error: File {args.file2} does not exist", file=sys.stderr)
        sys.exit(1)
    
    try:
        # Detect file formats
        format1 = detect_file_format(args.file1)
        format2 = detect_file_format(args.file2)
        
        if args.verbose:
            print(f"File 1 ({args.file1}) detected as: {format1}")
            print(f"File 2 ({args.file2}) detected as: {format2}")
        
        # Parse files into similarity matrices
        if not args.table:
            print("Parsing similarity matrices...")
        if format1 == 'bedtools':
            matrix1 = parse_bedtools_format(args.file1)
        else:
            matrix1 = parse_hammock_format(args.file1)
            
        if format2 == 'bedtools':
            matrix2 = parse_bedtools_format(args.file2)
        else:
            matrix2 = parse_hammock_format(args.file2)
        
        # Align matrices to common files
        if not args.table:
            print("Aligning matrices...")
        aligned_matrix1, aligned_matrix2 = align_matrices(matrix1, matrix2, quiet=args.table)
        
        # Calculate Frobenius norm
        if not args.table:
            print("Calculating Frobenius norm...")
        frobenius_norm = calculate_frobenius_norm(aligned_matrix1, aligned_matrix2)
        
        # Calculate additional metrics
        additional_metrics = calculate_additional_metrics(aligned_matrix1, aligned_matrix2)
        
        # Output results in table format if requested
        if args.table:
            # Extract just the filename from the full path
            file1_name = Path(args.file1).name
            file2_name = Path(args.file2).name
            
            # Print tab-delimited line
            print(f"{file1_name}\t{file2_name}\t{format1}\t{format2}\t{additional_metrics['matrix_size']}\t{frobenius_norm:.6f}\t{additional_metrics['correlation']:.6f}\t{additional_metrics['mean_absolute_error']:.6f}\t{additional_metrics['max_absolute_error']:.6f}")
        else:
            # Print results in normal format
            print("\n" + "="*50)
            print("SIMILARITY MATRIX COMPARISON RESULTS")
            print("="*50)
            print(f"File 1: {args.file1} ({format1} format)")
            print(f"File 2: {args.file2} ({format2} format)")
            print(f"Matrix size: {additional_metrics['matrix_size']} x {additional_metrics['matrix_size']}")
            print(f"Frobenius norm: {frobenius_norm:.6f}")
            print(f"Correlation: {additional_metrics['correlation']:.6f}")
            print(f"Mean absolute error: {additional_metrics['mean_absolute_error']:.6f}")
            print(f"Max absolute error: {additional_metrics['max_absolute_error']:.6f}")
        
        # Save detailed output if requested
        if args.output:
            with open(args.output, 'w') as f:
                f.write("SIMILARITY MATRIX COMPARISON RESULTS\n")
                f.write("="*50 + "\n")
                f.write(f"File 1: {args.file1} ({format1} format)\n")
                f.write(f"File 2: {args.file2} ({format2} format)\n")
                f.write(f"Matrix size: {additional_metrics['matrix_size']} x {additional_metrics['matrix_size']}\n")
                f.write(f"Frobenius norm: {frobenius_norm:.6f}\n")
                f.write(f"Correlation: {additional_metrics['correlation']:.6f}\n")
                f.write(f"Mean absolute error: {additional_metrics['mean_absolute_error']:.6f}\n")
                f.write(f"Max absolute error: {additional_metrics['max_absolute_error']:.6f}\n")
                f.write("\nDifference Matrix (Matrix1 - Matrix2):\n")
                diff_matrix = aligned_matrix1 - aligned_matrix2
                f.write(diff_matrix.to_string())
            if not args.table:
                print(f"\nDetailed results saved to: {args.output}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
