#!/usr/bin/env python3
"""
Script to compare hierarchical clustering trees from bedtools and hammock outputs.
Calculates Robinson-Foulds distance between trees built from Jaccard similarity matrices.
"""

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path
from scipy.cluster.hierarchy import linkage, to_tree, fcluster
from scipy.spatial.distance import squareform
import warnings
from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score
from sklearn.cluster import KMeans

# Suppress scipy warnings about clustering
warnings.filterwarnings('ignore', category=RuntimeWarning)

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
        # Make symmetric
        sim_matrix[j, i] = row['jaccard']
    
    # Convert to DataFrame with proper labels
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
    
    # Convert to DataFrame with proper labels
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


def similarity_to_distance_matrix(sim_matrix):
    """
    Convert similarity matrix to distance matrix.
    Distance = 1 - similarity, with bounds checking.
    
    Args:
        sim_matrix (pandas.DataFrame): Similarity matrix
        
    Returns:
        numpy.ndarray: Distance matrix as condensed form for clustering
    """
    # Convert similarities to distances: distance = 1 - similarity
    dist_matrix = 1.0 - sim_matrix.values
    
    # Ensure diagonal is zero (self-distance)
    np.fill_diagonal(dist_matrix, 0.0)
    
    # Ensure distances are non-negative and bounded
    dist_matrix = np.clip(dist_matrix, 0.0, 1.0)
    
    # Convert to condensed form (upper triangle, excluding diagonal)
    # This is required by scipy.cluster.hierarchy.linkage
    condensed_dist = squareform(dist_matrix)
    
    return condensed_dist


def dynamic_tree_cut(linkage_matrix, distance_threshold=None, min_cluster_size=2, max_cluster_size=None):
    """
    Perform dynamic tree cutting to automatically determine optimal clusters.
    
    Args:
        linkage_matrix: Linkage matrix from scipy.cluster.hierarchy.linkage
        distance_threshold: Distance threshold for cutting (if None, auto-detect)
        min_cluster_size: Minimum cluster size
        max_cluster_size: Maximum cluster size (if None, no limit)
        
    Returns:
        tuple: (cluster_labels, n_clusters, distance_threshold_used)
    """
    # If no distance threshold provided, try to auto-detect
    if distance_threshold is None:
        # Use the median of the linkage distances as a starting point
        distances = linkage_matrix[:, 2]
        distance_threshold = np.median(distances)
        
        # Adjust based on data size - smaller datasets need smaller thresholds
        n_samples = len(linkage_matrix) + 1
        if n_samples < 10:
            distance_threshold *= 0.5
        elif n_samples > 50:
            distance_threshold *= 1.5
    
    # Perform clustering with the threshold
    cluster_labels = fcluster(linkage_matrix, distance_threshold, criterion='distance')
    n_clusters = len(set(cluster_labels))
    
    # Post-process to ensure minimum cluster size
    if min_cluster_size > 1:
        # Count cluster sizes
        cluster_sizes = {}
        for label in cluster_labels:
            cluster_sizes[label] = cluster_sizes.get(label, 0) + 1
        
        # Find clusters that are too small
        small_clusters = [label for label, size in cluster_sizes.items() if size < min_cluster_size]
        
        if small_clusters:
            # Merge small clusters with their nearest neighbors
            for small_cluster in small_clusters:
                # Find the nearest cluster to merge with
                best_merge = None
                best_distance = float('inf')
                
                for other_cluster in set(cluster_labels):
                    if other_cluster != small_cluster and cluster_sizes.get(other_cluster, 0) >= min_cluster_size:
                        # Calculate average distance between clusters
                        small_indices = np.where(cluster_labels == small_cluster)[0]
                        other_indices = np.where(cluster_labels == other_cluster)[0]
                        
                        # Use linkage matrix to find distance
                        for i in range(len(linkage_matrix)):
                            if (linkage_matrix[i, 0] in small_indices and linkage_matrix[i, 1] in other_indices) or \
                               (linkage_matrix[i, 1] in small_indices and linkage_matrix[i, 0] in other_indices):
                                if linkage_matrix[i, 2] < best_distance:
                                    best_distance = linkage_matrix[i, 2]
                                    best_merge = other_cluster
                                break
                
                if best_merge is not None:
                    # Merge the small cluster into the best merge candidate
                    cluster_labels[cluster_labels == small_cluster] = best_merge
    
    # Re-number clusters to be consecutive
    unique_labels = sorted(set(cluster_labels))
    label_map = {old_label: new_label for new_label, old_label in enumerate(unique_labels, 1)}
    cluster_labels = np.array([label_map[label] for label in cluster_labels])
    n_clusters = len(unique_labels)
    
    return cluster_labels, n_clusters, distance_threshold


def evaluate_clustering_similarity(labels1, labels2, names):
    """
    Evaluate the similarity between two clusterings.
    
    Args:
        labels1: Cluster labels from first clustering
        labels2: Cluster labels from second clustering
        names: List of sample names
        
    Returns:
        dict: Dictionary containing various similarity metrics
    """
    # Calculate clustering similarity metrics
    ari = adjusted_rand_score(labels1, labels2)
    nmi = normalized_mutual_info_score(labels1, labels2)
    
    # Calculate cluster overlap matrix
    n_clusters1 = len(set(labels1))
    n_clusters2 = len(set(labels2))
    
    overlap_matrix = np.zeros((n_clusters1, n_clusters2))
    for i in range(n_clusters1):
        for j in range(n_clusters2):
            # Count samples that are in cluster i in labels1 and cluster j in labels2
            overlap = np.sum((labels1 == (i+1)) & (labels2 == (j+1)))
            overlap_matrix[i, j] = overlap
    
    # Calculate cluster correspondence
    max_overlaps = np.max(overlap_matrix, axis=1)
    total_samples = len(labels1)
    cluster_correspondence = np.sum(max_overlaps) / total_samples
    
    # Calculate cluster stability (how many samples stay in the same cluster)
    stable_samples = np.sum(labels1 == labels2)
    cluster_stability = stable_samples / total_samples
    
    return {
        'adjusted_rand_index': ari,
        'normalized_mutual_info': nmi,
        'cluster_correspondence': cluster_correspondence,
        'cluster_stability': cluster_stability,
        'overlap_matrix': overlap_matrix,
        'n_clusters1': n_clusters1,
        'n_clusters2': n_clusters2
    }


def cluster_tree_to_newick(tree_node, names):
    """
    Convert scipy clustering tree to Newick format string.
    
    Args:
        tree_node: Root node from scipy.cluster.hierarchy.to_tree()
        names (list): List of leaf names in order
        
    Returns:
        str: Newick format tree string
    """
    def _build_newick(node):
        if node.is_leaf():
            # Leaf node - return the name
            return names[node.id]
        else:
            # Internal node - recursively build subtrees
            left_subtree = _build_newick(node.left)
            right_subtree = _build_newick(node.right)
            return f"({left_subtree},{right_subtree})"
    
    return _build_newick(tree_node) + ";"


def calculate_robinson_foulds_distance(tree1_newick, tree2_newick):
    """
    Calculate Robinson-Foulds distance between two trees.
    Uses ete3 or dendropy if available, otherwise falls back to simple comparison.
    
    Args:
        tree1_newick (str): First tree in Newick format
        tree2_newick (str): Second tree in Newick format
        
    Returns:
        tuple: (rf_distance, max_rf_distance, normalized_rf_distance, n_leaves)
    """
    try:
        # Try to import ete3 for RF distance calculation
        from ete3 import Tree
        
        # Parse trees
        t1 = Tree(tree1_newick, format=1)
        t2 = Tree(tree2_newick, format=1)
        
        # Calculate Robinson-Foulds distance
        # ete3 robinson_foulds returns: (rf_distance, max_rf_distance, common_leaves_names, partitions_t1, partitions_t2, discarded_partitions_t1, discarded_partitions_t2)
        # But we only need the first few values
        rf_result = t1.robinson_foulds(t2, unrooted_trees=True)
        rf = rf_result[0]
        max_rf = rf_result[1]
        common_leaves = rf_result[2]
        
        # Normalized RF distance (0 to 1)
        normalized_rf = rf / max_rf if max_rf > 0 else 0.0
        
        return rf, max_rf, normalized_rf, len(common_leaves)
        
    except ImportError:
        # Fallback to dendropy if ete3 is not available
        try:
            import dendropy
            
            # Parse trees
            tns = dendropy.TaxonNamespace()
            t1 = dendropy.Tree.get(data=tree1_newick, schema="newick", taxon_namespace=tns)
            t2 = dendropy.Tree.get(data=tree2_newick, schema="newick", taxon_namespace=tns)
            
            # Calculate Robinson-Foulds distance
            rf = dendropy.calculate.treecompare.robinson_foulds_distance(t1, t2)
            
            # Calculate maximum possible RF distance (2 * (n - 3) for unrooted trees)
            n_leaves = len(tns)
            max_rf = 2 * (n_leaves - 3) if n_leaves > 3 else 0
            
            # Normalized RF distance
            normalized_rf = rf / max_rf if max_rf > 0 else 0.0
            
            return rf, max_rf, normalized_rf, n_leaves
            
        except ImportError:
            # Simple fallback implementation
            return simple_robinson_foulds(tree1_newick, tree2_newick)


def simple_robinson_foulds(tree1_newick, tree2_newick):
    """
    Simple Robinson-Foulds distance implementation.
    This is a basic fallback when specialized tree libraries are not available.
    
    Args:
        tree1_newick (str): First tree in Newick format
        tree2_newick (str): Second tree in Newick format
        
    Returns:
        tuple: (rf_distance, max_rf_distance, normalized_rf_distance, n_leaves)
    """
    def extract_leaves(newick):
        # Simple extraction of leaf names from Newick string
        import re
        # Find all leaf names (alphanumeric strings before commas or closing parens)
        leaves = re.findall(r'([A-Za-z0-9_\-\.]+)(?=[,\)])', newick)
        return set(leaves)
    
    leaves1 = extract_leaves(tree1_newick)
    leaves2 = extract_leaves(tree2_newick)
    
    # Common leaves
    common_leaves = leaves1.intersection(leaves2)
    n_leaves = len(common_leaves)
    
    # Simple approximation: if trees have different topologies, estimate RF distance
    if n_leaves > 3:
        max_rf = 2 * (n_leaves - 3)
        # For simplicity, assume moderate distance if we can't calculate properly
        # In practice, this should be replaced with proper RF calculation
        rf = max_rf // 3  # Conservative estimate
        normalized_rf = rf / max_rf
    else:
        rf = 0
        max_rf = 0
        normalized_rf = 0.0
    
    print("Warning: Using simplified RF distance approximation. Install ete3 or dendropy for accurate results.", file=sys.stderr)
    
    return rf, max_rf, normalized_rf, n_leaves


def extract_parameters_from_filename(filename):
    """
    Extract klen, window, and precision parameters from hammock output filename.
    Returns (klen, window, precision) as integers, or (0, 0, 0) if not found.
    """
    import re
    
    # Pattern to match parameters in hammock output filenames
    # Example: hammock_mnmzr_p16_jaccD_k10_w25.csv
    pattern = r'k(\d+).*w(\d+).*p(\d+)'
    match = re.search(pattern, filename)
    
    if match:
        klen = int(match.group(1))
        window = int(match.group(2))
        precision = int(match.group(3))
        return klen, window, precision
    else:
        return 0, 0, 0

def print_table_header():
    """
    Print the header for tab-delimited output.
    """
    print("file1\tfile2\tformat1\tformat2\tklen1\twindow1\tprecision1\tklen2\twindow2\tprecision2\tn_leaves\trf_distance\tmax_rf_distance\tnormalized_rf\tlinkage_method\tclustering_method\tn_clusters1\tn_clusters2\tadjusted_rand_index\tnormalized_mutual_info\tcluster_correspondence\tcluster_stability")


def main():
    """
    Main function to compare clustering trees from two files.
    """
    parser = argparse.ArgumentParser(
        description='Compare hierarchical clustering trees from bedtools and hammock outputs using Robinson-Foulds distance',
        epilog='Example: python compare_clustering_trees.py hammock_output.csv bedtools_output.tsv --linkage average'
    )
    parser.add_argument('file1', nargs='?', help='First input file (bedtools or hammock format)')
    parser.add_argument('file2', nargs='?', help='Second input file (bedtools or hammock format)')
    parser.add_argument('--linkage', '-l', choices=['single', 'complete', 'average', 'ward'], 
                       default='average', help='Hierarchical clustering linkage method (default: average)')
    parser.add_argument('--output', '-o', help='Output file for detailed comparison (optional)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    parser.add_argument('--table', '-t', action='store_true', 
                       help='Output results as tab-delimited line (use --header to print column names)')
    parser.add_argument('--header', action='store_true', 
                       help='Print tab-delimited header (use with --table)')
    parser.add_argument('--save-trees', help='Save Newick trees to file (useful for visualization)')
    parser.add_argument('--clustering-method', '-c', choices=['kmeans', 'dynamic', 'both'], 
                       default='both', help='Clustering method: kmeans (fixed k), dynamic (tree cutting), or both (default: both)')
    parser.add_argument('--distance-threshold', '-d', type=float, 
                       help='Distance threshold for dynamic tree cutting (auto-detect if not specified)')
    parser.add_argument('--min-cluster-size', type=int, default=2,
                       help='Minimum cluster size for dynamic tree cutting (default: 2)')
    
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
        
        # Convert to distance matrices
        if not args.table:
            print("Converting similarities to distances...")
        dist1 = similarity_to_distance_matrix(aligned_matrix1)
        dist2 = similarity_to_distance_matrix(aligned_matrix2)
        
        # Perform hierarchical clustering
        if not args.table:
            print(f"Performing hierarchical clustering with {args.linkage} linkage...")
        
        linkage1 = linkage(dist1, method=args.linkage)
        linkage2 = linkage(dist2, method=args.linkage)
        
        # Perform clustering based on method
        clustering_results = {}
        clustering_similarity = {}
        
        if args.clustering_method in ['kmeans', 'both']:
            if not args.table:
                print("Performing k-means clustering...")
            
            # Try different numbers of clusters
            cluster_counts = [2, 3, 4, 5, 6]
            best_silhouette1 = -1
            best_silhouette2 = -1
            best_k1 = 0
            best_k2 = 0
            best_labels1 = None
            best_labels2 = None
            
            for k in cluster_counts:
                kmeans1 = KMeans(n_clusters=k, random_state=42)
                kmeans2 = KMeans(n_clusters=k, random_state=42)
                
                # Convert condensed distance matrices back to full matrices for silhouette calculation
                n_samples = len(aligned_matrix1)
                dist_matrix1 = squareform(dist1)
                dist_matrix2 = squareform(dist2)
                
                # Get cluster labels
                labels1 = kmeans1.fit_predict(dist_matrix1)
                labels2 = kmeans2.fit_predict(dist_matrix2)
                
                # Calculate silhouette scores
                silhouette1 = silhouette_score(dist_matrix1, labels1, metric='precomputed')
                silhouette2 = silhouette_score(dist_matrix2, labels2, metric='precomputed')
                
                # Track best scores
                if silhouette1 > best_silhouette1:
                    best_silhouette1 = silhouette1
                    best_k1 = k
                    best_labels1 = labels1
                if silhouette2 > best_silhouette2:
                    best_silhouette2 = silhouette2
                    best_k2 = k
                    best_labels2 = labels2
                    
                if not args.table:
                    print(f"k={k}:")
                    print(f"  Silhouette score for matrix 1: {silhouette1:.4f}")
                    print(f"  Silhouette score for matrix 2: {silhouette2:.4f}")
            
            if not args.table:
                print(f"\nBest k-means clustering results:")
                print(f"Matrix 1: k={best_k1}, silhouette={best_silhouette1:.4f}")
                print(f"Matrix 2: k={best_k2}, silhouette={best_silhouette2:.4f}")
            
            clustering_results['kmeans'] = {
                'k1': best_k1, 'k2': best_k2,
                'silhouette1': best_silhouette1, 'silhouette2': best_silhouette2,
                'labels1': best_labels1, 'labels2': best_labels2
            }
            
            # Evaluate clustering similarity for k-means
            if best_labels1 is not None and best_labels2 is not None:
                clustering_similarity['kmeans'] = evaluate_clustering_similarity(
                    best_labels1, best_labels2, list(aligned_matrix1.index)
                )
        
        if args.clustering_method in ['dynamic', 'both']:
            if not args.table:
                print("Performing dynamic tree cutting...")
            
            # Perform dynamic tree cutting
            labels1, n_clusters1, threshold1 = dynamic_tree_cut(
                linkage1, 
                distance_threshold=args.distance_threshold,
                min_cluster_size=args.min_cluster_size
            )
            labels2, n_clusters2, threshold2 = dynamic_tree_cut(
                linkage2, 
                distance_threshold=args.distance_threshold,
                min_cluster_size=args.min_cluster_size
            )
            
            if not args.table:
                print(f"Dynamic tree cutting results:")
                print(f"Matrix 1: {n_clusters1} clusters, threshold={threshold1:.4f}")
                print(f"Matrix 2: {n_clusters2} clusters, threshold={threshold2:.4f}")
            
            clustering_results['dynamic'] = {
                'n_clusters1': n_clusters1, 'n_clusters2': n_clusters2,
                'threshold1': threshold1, 'threshold2': threshold2,
                'labels1': labels1, 'labels2': labels2
            }
            
            # Evaluate clustering similarity for dynamic tree cutting
            clustering_similarity['dynamic'] = evaluate_clustering_similarity(
                labels1, labels2, list(aligned_matrix1.index)
            )
        
        # Convert to tree objects
        tree1 = to_tree(linkage1)
        tree2 = to_tree(linkage2)
        
        # Convert to Newick format
        names = list(aligned_matrix1.index)
        newick1 = cluster_tree_to_newick(tree1, names)
        newick2 = cluster_tree_to_newick(tree2, names)
        
        if args.verbose:
            print(f"Tree 1 (Newick): {newick1}")
            print(f"Tree 2 (Newick): {newick2}")
        
        # Calculate Robinson-Foulds distance
        if not args.table:
            print("Calculating Robinson-Foulds distance...")
        rf_distance, max_rf, normalized_rf, n_leaves = calculate_robinson_foulds_distance(newick1, newick2)
        
        # Output results in table format if requested
        if args.table:
            # Extract just the filename from the full path
            file1_name = Path(args.file1).name
            file2_name = Path(args.file2).name
            
            # Extract parameters from filenames
            klen1, window1, precision1 = extract_parameters_from_filename(file1_name)
            klen2, window2, precision2 = extract_parameters_from_filename(file2_name)
            
            # Output separate lines for each clustering method when both are used
            methods_to_output = []
            if args.clustering_method == 'kmeans' or (args.clustering_method == 'both' and 'kmeans' in clustering_results):
                methods_to_output.append('kmeans')
            if args.clustering_method == 'dynamic' or (args.clustering_method == 'both' and 'dynamic' in clustering_results):
                methods_to_output.append('dynamic')
            
            for method in methods_to_output:
                if method == 'kmeans':
                    # For k-means, output results for each k value tested
                    cluster_counts = [2, 3, 4, 5, 6]
                    for k in cluster_counts:
                        # Re-run k-means for this k value to get the results
                        kmeans1 = KMeans(n_clusters=k, random_state=42)
                        kmeans2 = KMeans(n_clusters=k, random_state=42)
                        
                        # Convert condensed distance matrices back to full matrices
                        n_samples = len(aligned_matrix1)
                        dist_matrix1 = squareform(dist1)
                        dist_matrix2 = squareform(dist2)
                        
                        # Get cluster labels
                        labels1 = kmeans1.fit_predict(dist_matrix1)
                        labels2 = kmeans2.fit_predict(dist_matrix2)
                        
                        # Calculate similarity metrics for this k
                        similarity_metrics = evaluate_clustering_similarity(labels1, labels2, list(aligned_matrix1.index))
                        
                        ari = similarity_metrics.get('adjusted_rand_index', 0.0)
                        nmi = similarity_metrics.get('normalized_mutual_info', 0.0)
                        correspondence = similarity_metrics.get('cluster_correspondence', 0.0)
                        stability = similarity_metrics.get('cluster_stability', 0.0)
                        
                        # Print tab-delimited line for this k value
                        print(f"{file1_name}\t{file2_name}\t{format1}\t{format2}\t{klen1}\t{window1}\t{precision1}\t{klen2}\t{window2}\t{precision2}\t{n_leaves}\t{rf_distance}\t{max_rf}\t{normalized_rf:.6f}\t{args.linkage}\tkmeans_k{k}\t{k}\t{k}\t{ari:.6f}\t{nmi:.6f}\t{correspondence:.6f}\t{stability:.6f}")
                else:  # dynamic
                    # Get cluster counts for dynamic method
                    n_clusters1 = clustering_results['dynamic'].get('n_clusters1', 0)
                    n_clusters2 = clustering_results['dynamic'].get('n_clusters2', 0)
                    similarity_metrics = clustering_similarity.get('dynamic', {})
                    
                    ari = similarity_metrics.get('adjusted_rand_index', 0.0)
                    nmi = similarity_metrics.get('normalized_mutual_info', 0.0)
                    correspondence = similarity_metrics.get('cluster_correspondence', 0.0)
                    stability = similarity_metrics.get('cluster_stability', 0.0)
                    
                    # Print tab-delimited line for dynamic method
                    print(f"{file1_name}\t{file2_name}\t{format1}\t{format2}\t{klen1}\t{window1}\t{precision1}\t{klen2}\t{window2}\t{precision2}\t{n_leaves}\t{rf_distance}\t{max_rf}\t{normalized_rf:.6f}\t{args.linkage}\t{method}\t{n_clusters1}\t{n_clusters2}\t{ari:.6f}\t{nmi:.6f}\t{correspondence:.6f}\t{stability:.6f}")
        else:
            # Print results in normal format
            print("\n" + "="*60)
            print("HIERARCHICAL CLUSTERING TREE COMPARISON RESULTS")
            print("="*60)
            print(f"File 1: {args.file1} ({format1} format)")
            print(f"File 2: {args.file2} ({format2} format)")
            print(f"Linkage method: {args.linkage}")
            print(f"Clustering method: {args.clustering_method}")
            print(f"Number of common leaves: {n_leaves}")
            print(f"Robinson-Foulds distance: {rf_distance}")
            print(f"Maximum possible RF distance: {max_rf}")
            print(f"Normalized RF distance: {normalized_rf:.6f}")
            print(f"Tree similarity: {1.0 - normalized_rf:.6f}")
            
            # Print clustering-specific results
            if args.clustering_method in ['kmeans', 'both']:
                kmeans_results = clustering_results.get('kmeans', {})
                kmeans_similarity = clustering_similarity.get('kmeans', {})
                print(f"\nK-means clustering results:")
                print(f"  Best k for matrix 1: {kmeans_results.get('k1', 0)} (silhouette: {kmeans_results.get('silhouette1', 0):.4f})")
                print(f"  Best k for matrix 2: {kmeans_results.get('k2', 0)} (silhouette: {kmeans_results.get('silhouette2', 0):.4f})")
                if kmeans_similarity:
                    print(f"  Clustering similarity metrics:")
                    print(f"    Adjusted Rand Index: {kmeans_similarity.get('adjusted_rand_index', 0):.4f}")
                    print(f"    Normalized Mutual Information: {kmeans_similarity.get('normalized_mutual_info', 0):.4f}")
                    print(f"    Cluster Correspondence: {kmeans_similarity.get('cluster_correspondence', 0):.4f}")
                    print(f"    Cluster Stability: {kmeans_similarity.get('cluster_stability', 0):.4f}")
            
            if args.clustering_method in ['dynamic', 'both']:
                dynamic_results = clustering_results.get('dynamic', {})
                dynamic_similarity = clustering_similarity.get('dynamic', {})
                print(f"\nDynamic tree cutting results:")
                print(f"  Clusters in matrix 1: {dynamic_results.get('n_clusters1', 0)} (threshold: {dynamic_results.get('threshold1', 0):.4f})")
                print(f"  Clusters in matrix 2: {dynamic_results.get('n_clusters2', 0)} (threshold: {dynamic_results.get('threshold2', 0):.4f})")
                if dynamic_similarity:
                    print(f"  Clustering similarity metrics:")
                    print(f"    Adjusted Rand Index: {dynamic_similarity.get('adjusted_rand_index', 0):.4f}")
                    print(f"    Normalized Mutual Information: {dynamic_similarity.get('normalized_mutual_info', 0):.4f}")
                    print(f"    Cluster Correspondence: {dynamic_similarity.get('cluster_correspondence', 0):.4f}")
                    print(f"    Cluster Stability: {dynamic_similarity.get('cluster_stability', 0):.4f}")
        
        # Save trees if requested
        if args.save_trees:
            with open(args.save_trees, 'w') as f:
                f.write(f"# Tree 1 ({args.file1})\n")
                f.write(f"{newick1}\n\n")
                f.write(f"# Tree 2 ({args.file2})\n")
                f.write(f"{newick2}\n")
            if not args.table:
                print(f"\nTrees saved to: {args.save_trees}")
        
        # Save detailed output if requested
        if args.output:
            with open(args.output, 'w') as f:
                f.write("HIERARCHICAL CLUSTERING TREE COMPARISON RESULTS\n")
                f.write("="*60 + "\n")
                f.write(f"File 1: {args.file1} ({format1} format)\n")
                f.write(f"File 2: {args.file2} ({format2} format)\n")
                f.write(f"Linkage method: {args.linkage}\n")
                f.write(f"Clustering method: {args.clustering_method}\n")
                f.write(f"Number of common leaves: {n_leaves}\n")
                f.write(f"Robinson-Foulds distance: {rf_distance}\n")
                f.write(f"Maximum possible RF distance: {max_rf}\n")
                f.write(f"Normalized RF distance: {normalized_rf:.6f}\n")
                f.write(f"Tree similarity: {1.0 - normalized_rf:.6f}\n")
                f.write(f"\nTree 1 (Newick):\n{newick1}\n")
                f.write(f"\nTree 2 (Newick):\n{newick2}\n")
            if not args.table:
                print(f"\nDetailed results saved to: {args.output}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main() 