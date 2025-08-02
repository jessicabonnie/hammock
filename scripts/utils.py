#!/usr/bin/env python3
"""
Utility functions for hammock analysis scripts.
Shared functions used across multiple analysis scripts.
"""

import pandas as pd
import numpy as np
from pathlib import Path


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
        # Mode B (BED files) or other modes - use jaccard_similarity
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


def load_accession_key(filepath):
    """
    Load the accession key file and create a mapping from file to true labels.
    
    Args:
        filepath (str): Path to the accession key TSV file
        
    Returns:
        dict: Mapping from file basename to true label (tissue type)
    """
    # Read the accession key file
    df = pd.read_csv(filepath, sep='\t')
    
    # Create mapping from file to tissue type
    file_to_tissue = {}
    
    for _, row in df.iterrows():
        file_basename = Path(row['File']).stem  # Remove .fa extension
        tissue = row['Biosample_term_name']
        file_to_tissue[file_basename] = tissue
    
    return file_to_tissue


def filter_hammock_by_accessions(hammock_file, accession_list, output_file=None):
    """
    Filter hammock output to only include rows where both file1 and file2 
    have basenames in the accession list.
    
    Args:
        hammock_file (str): Path to hammock CSV output file
        accession_list (set): Set of accession basenames to keep
        output_file (str): Output file path (if None, returns DataFrame)
        
    Returns:
        pandas.DataFrame or int: Filtered DataFrame or number of rows if output_file is specified
    """
    # Read hammock output
    df = pd.read_csv(hammock_file)
    
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
    
    # Write output or return DataFrame
    if output_file:
        filtered_df.to_csv(output_file, index=False)
        return len(filtered_df)
    else:
        return filtered_df 