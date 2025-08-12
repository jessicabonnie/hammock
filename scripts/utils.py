#!/usr/bin/env python3
"""
Utility functions for hammock analysis scripts.
Shared functions used across multiple analysis scripts.
"""

import pandas as pd
import numpy as np
from pathlib import Path


def parse_hammock_format(filepath, preferred_similarity_column=None):
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
    if preferred_similarity_column and preferred_similarity_column in df.columns:
        jaccard_col = preferred_similarity_column
    else:
        if 'jaccard_similarity_with_ends' in df.columns:
            # Mode D (FASTA files) - use jaccard_similarity_with_ends
            jaccard_col = 'jaccard_similarity_with_ends'
        elif 'jaccard_similarity' in df.columns:
            # Modes A/B/C - use jaccard_similarity
            jaccard_col = 'jaccard_similarity'
        else:
            raise ValueError("No jaccard similarity column found in hammock output")
    # Normalize to a canonical column name for downstream use
    df['similarity'] = df[jaccard_col]
    
    # Extract similarity values and file pairs
    similarity_data = df[['file1', 'file2', 'similarity']].copy()
    
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
        sim_matrix[i, j] = row['similarity']
        # Make symmetric
        sim_matrix[j, i] = row['similarity']
    
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


def detect_hammock_similarity_column(filepath):
    """
    Detect which Jaccard similarity column to use from a hammock CSV.
    Modes A/B/C use 'jaccard_similarity', mode D uses 'jaccard_similarity_with_ends'.

    Args:
        filepath (str): Path to hammock CSV file

    Returns:
        str: Column name to use for similarity
    """
    df_head = pd.read_csv(filepath, nrows=1)
    # Prefer column based on filename-indicated mode when possible
    params = extract_hammock_parameters_from_filename(filepath)
    mode = (params.get('mode') or '').upper()
    if mode == 'D' and 'jaccard_similarity_with_ends' in df_head.columns:
        return 'jaccard_similarity_with_ends'
    if mode in {'A', 'B', 'C'} and 'jaccard_similarity' in df_head.columns:
        return 'jaccard_similarity'
    # Fallback to available columns
    if 'jaccard_similarity_with_ends' in df_head.columns:
        return 'jaccard_similarity_with_ends'
    if 'jaccard_similarity' in df_head.columns:
        return 'jaccard_similarity'
    raise ValueError("No recognized Jaccard similarity column found in hammock output")


def detect_hammock_exp(filepath):
    """
    Detect exponent parameter for mode C outputs from the hammock CSV content.
    Looks for columns like 'expA' or 'exp', returns a representative value.

    Args:
        filepath (str): Path to hammock CSV

    Returns:
        int|float|None: Exponent value if detected, else None
    """
    df_head = pd.read_csv(filepath, nrows=10)
    candidate_cols = [c for c in df_head.columns if c.lower() in {'exp', 'expa'}]
    for col in candidate_cols:
        values = df_head[col].dropna().unique()
        if len(values) > 0:
            try:
                return int(values[0])
            except Exception:
                try:
                    return float(values[0])
                except Exception:
                    return None
    return None


def extract_hammock_parameters_from_filename(filename):
    """
    Extract mode (A-D), klen (k), window (w), and precision (p) from a hammock output filename.

    Expected patterns include tokens like: 'jaccD', 'k10', 'w25', 'p16'.

    Args:
        filename (str): Base filename or full path

    Returns:
        dict: {
            'mode': str|None,
            'klen': int|None,
            'window': int|None,
            'precision': int|None,
            'sketch': str|None,              # hll, mnmzr, kmv, minhash, hyperloglog
            'species_scope': str|None,       # e.g., mouseonly, manmouse, humanonly
            'subset_tag': str|None,          # e.g., subtissue, balanced, etc.
            'tags': str|None,                # joined non-parameter tokens
            'subA': float|None,              # parsed from tokens like A0.50
            'subB': float|None,              # parsed from tokens like B0.60
            'expA': float|None,              # parsed from tokens like expA2.00
        }
    """
    import re
    name = Path(filename).name

    mode = None
    mode_match = re.search(r'jacc([A-D])', name, flags=re.IGNORECASE)
    if mode_match:
        mode = mode_match.group(1).upper()

    klen = None
    k_match = re.search(r'k(\d+)', name)
    if k_match:
        klen = int(k_match.group(1))

    window = None
    w_match = re.search(r'w(\d+)', name)
    if w_match:
        window = int(w_match.group(1))

    precision = None
    p_match = re.search(r'p(\d+)', name)
    if p_match:
        precision = int(p_match.group(1))

    # Detect sketch type tokens
    sketch = None
    lowered = name.lower()
    if re.search(r'(^|[_\.-])(hll|hyperloglog)([_\.-]|$)', lowered):
        sketch = 'hyperloglog'
    elif re.search(r'(^|[_\.-])(mnmzr|minimizer)([_\.-]|$)', lowered):
        sketch = 'minimizer'
    elif re.search(r'(^|[_\.-])(kmv|minhash)([_\.-]|$)', lowered):
        sketch = 'kmv'

    # Tokenize to capture leading descriptive tags
    base_no_ext = Path(name).stem
    tokens = re.split(r'[_\.-]+', base_no_ext)

    # Helper to identify parameter tokens
    def is_param_token(tok: str) -> bool:
        return (
            re.fullmatch(r'[pP]\d+', tok) is not None or
            re.fullmatch(r'jacc[ABCDabcd]', tok) is not None or
            re.fullmatch(r'[kK]\d+', tok) is not None or
            re.fullmatch(r'[wW]\d+', tok) is not None or
            tok.lower() in {'hll', 'hyperloglog', 'mnmzr', 'minimizer', 'kmv', 'minhash'}
        )

    # Extract optional subA/subB and expA tokens
    subA = None
    subB = None
    expA = None
    for tok in tokens:
        # expA token like expA2.00
        m_exp = re.fullmatch(r'expA(\d+(?:\.\d+)?)', tok, flags=re.IGNORECASE)
        if m_exp:
            try:
                expA = float(m_exp.group(1))
            except Exception:
                pass
            continue
        # A token like A0.75
        m_a = re.fullmatch(r'A(\d+(?:\.\d+)?)', tok)
        if m_a:
            try:
                subA = float(m_a.group(1))
            except Exception:
                pass
            continue
        # B token like B0.60
        m_b = re.fullmatch(r'B(\d+(?:\.\d+)?)', tok)
        if m_b:
            try:
                subB = float(m_b.group(1))
            except Exception:
                pass
            continue

    non_param_tokens = [t for t in tokens if t and not is_param_token(t)]
    species_scope = non_param_tokens[0] if len(non_param_tokens) >= 1 else None
    subset_tag = non_param_tokens[1] if len(non_param_tokens) >= 2 else None
    tags_joined = '_'.join(non_param_tokens) if non_param_tokens else None

    return {
        'mode': mode,
        'klen': klen,
        'window': window,
        'precision': precision,
        'sketch': sketch,
        'species_scope': species_scope,
        'subset_tag': subset_tag,
        'tags': tags_joined,
        'subA': subA,
        'subB': subB,
        'expA': expA,
    }


def detect_file_format(filepath):
    """
    Detect if the file is in bedtools or hammock format.

    Args:
        filepath (str): Path to the input file

    Returns:
        str: 'bedtools' or 'hammock'
    """
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline().strip()
        if (first_line.startswith('file1\tfile2\tintersection\tunion\tjaccard') or
                first_line.startswith('file1 file2 intersection union jaccard')):
            return 'bedtools'
        elif first_line.startswith('file1,file2') and (
                'jaccard_similarity' in first_line or 'jaccard_similarity_with_ends' in first_line):
            return 'hammock'
        else:
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
    Parse bedtools pairwise jaccard output into a similarity matrix using 'jaccard' column.

    Args:
        filepath (str): Path to bedtools output file

    Returns:
        pandas.DataFrame: Symmetric similarity matrix
    """
    # Read the bedtools output (mixed space/tab separated - use whitespace)
    df = pd.read_csv(filepath, sep=r'\s+')

    # Extract similarity values and file pairs
    if not {'file1', 'file2', 'jaccard'}.issubset(df.columns):
        raise ValueError("Bedtools file missing required columns: file1, file2, jaccard")
    similarity_data = df[['file1', 'file2', 'jaccard']].copy()

    # Convert file names to basenames without extensions for matching
    def get_basename(filename):
        name = str(filename)
        if name.endswith('.bed.gz'):
            name = name[:-7]
        elif name.endswith('.fa.gz'):
            name = name[:-6]
        else:
            name = Path(name).stem
        return name

    similarity_data['file1_base'] = similarity_data['file1'].apply(get_basename)
    similarity_data['file2_base'] = similarity_data['file2'].apply(get_basename)

    # Get unique file basenames
    all_files = sorted(set(similarity_data['file1_base'].tolist() + similarity_data['file2_base'].tolist()))

    # Create empty similarity matrix
    n_files = len(all_files)
    sim_matrix = np.zeros((n_files, n_files))

    # Create mapping
    file_to_idx = {file: idx for idx, file in enumerate(all_files)}

    # Fill matrix symmetrically
    for _, row in similarity_data.iterrows():
        i = file_to_idx[row['file1_base']]
        j = file_to_idx[row['file2_base']]
        sim_matrix[i, j] = row['jaccard']
        sim_matrix[j, i] = row['jaccard']

    return pd.DataFrame(sim_matrix, index=all_files, columns=all_files)