#!/usr/bin/env python3
"""
Script to generate hierarchical clustering trees from bedtools or hammock similarity matrices.
Automatically detects file format and creates dendrogram visualizations.
"""

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import matplotlib
import warnings

# Suppress scipy warnings about clustering
warnings.filterwarnings('ignore', category=RuntimeWarning)

def detect_file_format(filepath):
    """
    Detect if the file is in bedtools, hammock, or newick format.
    
    Args:
        filepath (str): Path to the input file
        
    Returns:
        str: Either 'bedtools', 'hammock', or 'newick'
    """
    try:
        # Read first few lines to detect format
        with open(filepath, 'r') as f:
            first_line = f.readline().strip()
            second_line = f.readline().strip()
            
        # Check for Newick format (comment followed by tree string)
        if (first_line.startswith('#') and 
            (second_line.endswith(';') or '(' in second_line)):
            return 'newick'
        
        # Check for bedtools format (has specific column names, tab or space separated)
        if (first_line.startswith('file1\tfile2\tintersection\tunion\tjaccard') or 
            first_line.startswith('file1 file2 intersection union jaccard')):
            return 'bedtools'
        elif first_line.startswith('file1,file2') and ('jaccard_similarity' in first_line or 'jaccard_similarity_with_ends' in first_line):
            return 'hammock'
        else:
            # Try to determine by content
            try:
                df = pd.read_csv(filepath, nrows=5)
                if 'jaccard' in df.columns and 'intersection' in df.columns:
                    return 'bedtools'
                elif 'jaccard_similarity' in df.columns or 'jaccard_similarity_with_ends' in df.columns:
                    return 'hammock'
                else:
                    raise ValueError(f"Cannot determine file format for {filepath}")
            except:
                # If CSV parsing fails, check if it's a simple Newick file
                with open(filepath, 'r') as f:
                    content = f.read().strip()
                if content.endswith(';') and '(' in content:
                    return 'newick'
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


def parse_newick_format(filepath):
    """
    Parse Newick format tree file and extract trees for visualization.
    
    Args:
        filepath (str): Path to Newick file
        
    Returns:
        list: List of dictionaries containing tree info and Newick strings
    """
    trees = []
    current_tree = None
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('#'):
                # Comment line with tree description
                if current_tree is not None:
                    trees.append(current_tree)
                
                # Extract description from comment
                description = line[1:].strip()
                current_tree = {
                    'description': description,
                    'newick': None
                }
            elif line.endswith(';'):
                # Newick tree string
                if current_tree is not None:
                    current_tree['newick'] = line
                else:
                    # Simple Newick file without comments
                    current_tree = {
                        'description': 'Tree',
                        'newick': line
                    }
                trees.append(current_tree)
                current_tree = None
    
    # Handle case where last tree doesn't end with semicolon on same line
    if current_tree is not None and current_tree['newick'] is None:
        # Try to read remaining content as Newick
        with open(filepath, 'r') as f:
            content = f.read()
            # Find the last line that looks like Newick
            for line in content.split('\n'):
                if '(' in line and ')' in line:
                    current_tree['newick'] = line.strip()
                    if not current_tree['newick'].endswith(';'):
                        current_tree['newick'] += ';'
                    break
        if current_tree['newick']:
            trees.append(current_tree)
    
    return trees


def newick_to_dendrogram(newick_string, tree_description="Tree", figsize=(12, 8), 
                        title=None, output_file=None, max_label_length=20, 
                        show_plot=True):
    """
    Create a dendrogram from a Newick format tree string.
    
    Args:
        newick_string (str): Newick format tree string
        tree_description (str): Description of the tree
        figsize (tuple): Figure size for the plot
        title (str): Title for the plot
        output_file (str): Path to save the plot (optional)
        max_label_length (int): Maximum length for leaf labels
        show_plot (bool): Whether to display the plot
        
    Returns:
        dict: Dictionary containing tree information
    """
    try:
        # Try to use ete3 for parsing and convert to scipy format for dendrogram
        from ete3 import Tree
        
        # Parse the Newick tree
        tree = Tree(newick_string, format=1)
        
        # Extract leaf names
        leaf_names = [leaf.name for leaf in tree.get_leaves()]
        
        # Truncate labels if they're too long
        if max_label_length > 0:
            display_labels = [label[:max_label_length] + '...' if len(label) > max_label_length else label 
                             for label in leaf_names]
        else:
            display_labels = leaf_names
        
        # Convert ete3 tree to distance matrix and then to scipy linkage format
        # This is a simplified approach - we'll create a dendrogram-style plot
        try:
            # Try to convert to scipy format for proper dendrogram
            distance_matrix = newick_to_distance_matrix(tree, leaf_names)
            
            # Create linkage matrix from distance matrix
            from scipy.spatial.distance import squareform
            from scipy.cluster.hierarchy import linkage, dendrogram
            
            condensed_dist = squareform(distance_matrix, checks=False)
            linkage_matrix = linkage(condensed_dist, method='average')
            
            # Create the plot
            plt.figure(figsize=figsize)
            
            # Create dendrogram
            dendrogram_info = dendrogram(
                linkage_matrix,
                labels=display_labels,
                orientation='top',
                distance_sort='descending',
                show_leaf_counts=True,
                leaf_rotation=90,
                leaf_font_size=8
            )
            
            # Set title
            if title:
                plt.title(title, fontsize=14, fontweight='bold')
            else:
                plt.title(f'Tree Visualization: {tree_description}', fontsize=14, fontweight='bold')
            
            # Set labels
            plt.xlabel('Samples', fontsize=12)
            plt.ylabel('Distance', fontsize=12)
            
            # Adjust layout to prevent label cutoff
            plt.tight_layout()
            
        except Exception as e:
            # Fallback to ete3's own plotting if scipy conversion fails
            print(f"Warning: Could not convert to scipy dendrogram, using ete3 tree plot: {e}")
            
            # Create matplotlib figure
            fig, ax = plt.subplots(figsize=figsize)
            
            # Use ete3's ASCII representation as a starting point
            # This is a simplified approach - we'll create a tree-like structure
            create_tree_plot_from_ete3(tree, ax, display_labels, max_label_length)
            
            if title:
                plt.title(title, fontsize=14, fontweight='bold')
            else:
                plt.title(f'Tree Visualization: {tree_description}', fontsize=14, fontweight='bold')
        
        # Save plot if requested
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Tree visualization saved to: {output_file}")
        
        # Show plot if requested
        if show_plot:
            plt.show()
        else:
            plt.close()
        
        return {
            'tree': tree,
            'leaf_names': leaf_names,
            'newick': newick_string,
            'description': tree_description
        }
        
    except ImportError:
        # Fallback without ete3 - create a simple tree visualization using manual parsing
        print("Warning: ete3 not available, using simplified tree visualization")
        
        # Extract leaf names using simple parsing
        import re
        leaf_names = re.findall(r'([A-Za-z0-9_\-\.]+)(?=[,\)])', newick_string)
        
        # Truncate labels if they're too long
        if max_label_length > 0:
            display_labels = [label[:max_label_length] + '...' if len(label) > max_label_length else label 
                             for label in leaf_names]
        else:
            display_labels = leaf_names
        
        # Create a simple tree visualization
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create a simple hierarchical layout
        create_simple_tree_plot(newick_string, ax, display_labels)
        
        if title:
            plt.title(title, fontsize=14, fontweight='bold')
        else:
            plt.title(f'Tree Visualization: {tree_description}', fontsize=14, fontweight='bold')
        
        plt.xlabel('Samples', fontsize=12)
        plt.ylabel('Hierarchical Structure', fontsize=12)
        plt.tight_layout()
        
        # Save plot if requested
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Tree visualization saved to: {output_file}")
        
        # Show plot if requested
        if show_plot:
            plt.show()
        else:
            plt.close()
        
        return {
            'tree': None,
            'leaf_names': leaf_names,
            'newick': newick_string,
            'description': tree_description
        }


def newick_to_distance_matrix(tree, leaf_names):
    """
    Convert ete3 tree to distance matrix.
    
    Args:
        tree: ete3 Tree object
        leaf_names: List of leaf names
        
    Returns:
        numpy.ndarray: Distance matrix
    """
    n = len(leaf_names)
    distance_matrix = np.zeros((n, n))
    
    # Calculate pairwise distances between all leaves
    leaves = tree.get_leaves()
    name_to_leaf = {leaf.name: leaf for leaf in leaves}
    
    for i, name1 in enumerate(leaf_names):
        for j, name2 in enumerate(leaf_names):
            if i != j:
                leaf1 = name_to_leaf[name1]
                leaf2 = name_to_leaf[name2]
                distance = leaf1.get_distance(leaf2)
                distance_matrix[i, j] = distance
    
    return distance_matrix


def create_tree_plot_from_ete3(tree, ax, display_labels, max_label_length):
    """
    Create a tree plot using ete3 tree structure.
    
    Args:
        tree: ete3 Tree object
        ax: matplotlib axis
        display_labels: List of display labels
        max_label_length: Maximum label length
    """
    # Get tree layout information
    leaves = tree.get_leaves()
    
    # Simple layout: arrange leaves horizontally and draw connections
    n_leaves = len(leaves)
    
    # Position leaves
    leaf_positions = {}
    for i, leaf in enumerate(leaves):
        x = i
        y = 0
        leaf_positions[leaf] = (x, y)
        
        # Draw leaf label
        label = leaf.name
        if max_label_length > 0 and len(label) > max_label_length:
            label = label[:max_label_length] + '...'
        ax.text(x, y, label, rotation=90, ha='center', va='bottom', fontsize=8)
    
    # Draw tree structure (simplified)
    def draw_node(node, level=1):
        if node.is_leaf():
            return leaf_positions[node]
        
        children = node.get_children()
        child_positions = [draw_node(child, level+1) for child in children]
        
        # Calculate position for internal node
        x_coords = [pos[0] for pos in child_positions]
        y_coord = level
        x_coord = sum(x_coords) / len(x_coords)
        
        # Draw connections to children
        for child_pos in child_positions:
            ax.plot([x_coord, child_pos[0]], [y_coord, child_pos[1]], 'k-', linewidth=1)
        
        # Draw horizontal line connecting children
        if len(child_positions) > 1:
            x_min = min(pos[0] for pos in child_positions)
            x_max = max(pos[0] for pos in child_positions)
            ax.plot([x_min, x_max], [y_coord, y_coord], 'k-', linewidth=1)
        
        return (x_coord, y_coord)
    
    # Draw the tree
    draw_node(tree)
    
    # Set axis properties
    ax.set_xlim(-0.5, n_leaves - 0.5)
    ax.set_ylim(-0.5, tree.get_farthest_leaf()[1] + 1)
    ax.set_aspect('auto')


def create_simple_tree_plot(newick_string, ax, display_labels):
    """
    Create a simple tree plot without ete3.
    
    Args:
        newick_string: Newick format string
        ax: matplotlib axis
        display_labels: List of display labels
    """
    # This is a very simplified tree visualization
    # Just show the leaf names in a hierarchical arrangement
    
    n_leaves = len(display_labels)
    
    # Arrange leaves in a simple layout
    for i, label in enumerate(display_labels):
        x = i
        y = 0
        ax.text(x, y, label, rotation=90, ha='center', va='bottom', fontsize=8)
        
        # Draw a simple vertical line to suggest tree structure
        ax.plot([x, x], [0, 1], 'k-', linewidth=1, alpha=0.5)
    
    # Draw a horizontal line at the top to suggest common ancestor
    if n_leaves > 1:
        ax.plot([0, n_leaves-1], [1, 1], 'k-', linewidth=2)
    
    # Set axis properties
    ax.set_xlim(-0.5, n_leaves - 0.5)
    ax.set_ylim(-0.5, 1.5)
    ax.set_aspect('auto')
    
    # Add note about simplified visualization
    ax.text(n_leaves/2, 1.2, 'Simplified tree structure\n(Install ete3 for detailed visualization)', 
            ha='center', va='bottom', fontsize=10, style='italic',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.7))


def similarity_to_distance_matrix(sim_matrix):
    """
    Convert similarity matrix to distance matrix for hierarchical clustering.
    
    Args:
        sim_matrix (pandas.DataFrame): Symmetric similarity matrix
        
    Returns:
        numpy.ndarray: Condensed distance matrix for scipy.cluster.hierarchy
    """
    # Convert similarities to distances (1 - similarity)
    dist_matrix = 1.0 - sim_matrix.values
    
    # Ensure diagonal is zero (perfect similarity with self)
    np.fill_diagonal(dist_matrix, 0.0)
    
    # Ensure distances are non-negative
    dist_matrix = np.maximum(dist_matrix, 0.0)
    
    # Convert to condensed form for scipy clustering
    # Only upper triangle (excluding diagonal) in row-wise order
    condensed_dist = squareform(dist_matrix, checks=False)
    
    return condensed_dist


def create_dendrogram(sim_matrix, linkage_method='average', figsize=(12, 8), 
                     title=None, output_file=None, max_label_length=20, 
                     color_threshold=None, show_plot=True):
    """
    Create and display/save a dendrogram from a similarity matrix.
    
    Args:
        sim_matrix (pandas.DataFrame): Symmetric similarity matrix
        linkage_method (str): Linkage method for clustering
        figsize (tuple): Figure size for the plot
        title (str): Title for the plot
        output_file (str): Path to save the plot (optional)
        max_label_length (int): Maximum length for leaf labels
        color_threshold (float): Threshold for coloring clusters
        show_plot (bool): Whether to display the plot
        
    Returns:
        dict: Dictionary containing linkage matrix and dendrogram info
    """
    # Convert to distance matrix
    dist_matrix = similarity_to_distance_matrix(sim_matrix)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(dist_matrix, method=linkage_method)
    
    # Truncate labels if they're too long
    labels = list(sim_matrix.index)
    if max_label_length > 0:
        labels = [label[:max_label_length] + '...' if len(label) > max_label_length else label 
                 for label in labels]
    
    # Create the plot
    plt.figure(figsize=figsize)
    
    # Create dendrogram
    dendrogram_info = dendrogram(
        linkage_matrix,
        labels=labels,
        orientation='top',
        distance_sort='descending',
        show_leaf_counts=True,
        color_threshold=color_threshold,
        leaf_rotation=90,
        leaf_font_size=8
    )
    
    # Set title
    if title:
        plt.title(title, fontsize=14, fontweight='bold')
    else:
        plt.title(f'Hierarchical Clustering Dendrogram ({linkage_method} linkage)', 
                 fontsize=14, fontweight='bold')
    
    # Set labels
    plt.xlabel('Samples', fontsize=12)
    plt.ylabel('Distance (1 - Jaccard Similarity)', fontsize=12)
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save plot if requested
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Dendrogram saved to: {output_file}")
    
    # Show plot if requested
    if show_plot:
        plt.show()
    else:
        plt.close()
    
    return {
        'linkage_matrix': linkage_matrix,
        'dendrogram_info': dendrogram_info,
        'labels': labels,
        'distance_matrix': dist_matrix
    }


def print_clustering_summary(sim_matrix, linkage_matrix, linkage_method):
    """
    Print a summary of the clustering results.
    
    Args:
        sim_matrix (pandas.DataFrame): Original similarity matrix
        linkage_matrix (numpy.ndarray): Linkage matrix from clustering
        linkage_method (str): Linkage method used
    """
    n_samples = len(sim_matrix)
    
    print("\nClustering Summary:")
    print("=" * 50)
    print(f"Number of samples: {n_samples}")
    print(f"Linkage method: {linkage_method}")
    print(f"Matrix size: {sim_matrix.shape[0]} x {sim_matrix.shape[1]}")
    
    # Calculate some basic statistics
    similarity_values = sim_matrix.values[np.triu_indices(n_samples, k=1)]
    print(f"Similarity range: {similarity_values.min():.4f} - {similarity_values.max():.4f}")
    print(f"Mean similarity: {similarity_values.mean():.4f}")
    print(f"Median similarity: {np.median(similarity_values):.4f}")
    
    # Distance statistics from linkage matrix
    distances = linkage_matrix[:, 2]
    print(f"Clustering distance range: {distances.min():.4f} - {distances.max():.4f}")
    print(f"Mean clustering distance: {distances.mean():.4f}")
    
    print("\nTop 10 most similar pairs:")
    print("-" * 30)
    
    # Find most similar pairs
    sim_values = sim_matrix.values.copy()
    np.fill_diagonal(sim_values, -1)  # Exclude diagonal
    
    # Get indices of top similarities
    flat_indices = np.argsort(sim_values.flatten())[-10:][::-1]
    row_indices, col_indices = np.unravel_index(flat_indices, sim_values.shape)
    
    for i, (row_idx, col_idx) in enumerate(zip(row_indices, col_indices)):
        if sim_values[row_idx, col_idx] > 0:  # Only show positive similarities
            sample1 = sim_matrix.index[row_idx]
            sample2 = sim_matrix.index[col_idx]
            similarity = sim_values[row_idx, col_idx]
            print(f"{i+1:2d}. {sample1[:25]:25s} <-> {sample2[:25]:25s} ({similarity:.4f})")


def main():
    """
    Main function to generate hierarchical clustering tree from similarity matrix.
    """
    parser = argparse.ArgumentParser(
        description='Generate hierarchical clustering tree from bedtools/hammock similarity matrix or Newick tree file',
        epilog='Examples:\n'
               '  python draw_clustering_tree.py input_matrix.csv --output tree.png --linkage complete\n'
               '  python draw_clustering_tree.py saved_trees.newick --output tree.png --tree-index 0',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('input_file', help='Input file: similarity matrix (bedtools/hammock format) or Newick tree file')
    parser.add_argument('--output', '-o', help='Output file for dendrogram (PNG, PDF, SVG, etc.)')
    parser.add_argument('--linkage', '-l', choices=['single', 'complete', 'average', 'ward'], 
                       default='average', help='Linkage method for clustering (default: average, ignored for Newick files)')
    parser.add_argument('--figsize', nargs=2, type=float, default=[12, 8], 
                       help='Figure size as width height (default: 12 8)')
    parser.add_argument('--title', '-t', help='Custom title for the dendrogram')
    parser.add_argument('--max-label-length', type=int, default=20, 
                       help='Maximum length for sample labels (default: 20, 0 for no truncation)')
    parser.add_argument('--color-threshold', type=float, 
                       help='Distance threshold for coloring clusters (default: auto, ignored for Newick files)')
    parser.add_argument('--no-show', action='store_true', 
                       help='Do not display the plot (useful when saving to file)')
    parser.add_argument('--summary', '-s', action='store_true', 
                       help='Print detailed clustering summary (ignored for Newick files)')
    parser.add_argument('--quiet', '-q', action='store_true', 
                       help='Suppress progress messages')
    parser.add_argument('--tree-index', type=int, default=0,
                       help='Index of tree to visualize from Newick file (default: 0, first tree)')
    parser.add_argument('--list-trees', action='store_true',
                       help='List all trees in Newick file and exit (use with Newick files)')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not Path(args.input_file).exists():
        print(f"Error: Input file {args.input_file} does not exist", file=sys.stderr)
        sys.exit(1)
    
    try:
        # Set matplotlib backend for non-interactive use if needed
        if args.no_show and args.output:
            matplotlib.use('Agg')
        
        # Detect file format
        file_format = detect_file_format(args.input_file)
        if not args.quiet:
            print(f"Detected file format: {file_format}")
        
        if file_format == 'newick':
            # Handle Newick format files
            if not args.quiet:
                print("Parsing Newick tree file...")
            
            trees = parse_newick_format(args.input_file)
            
            if args.list_trees:
                # List all trees and exit
                print(f"Found {len(trees)} tree(s) in file:")
                for i, tree_info in enumerate(trees):
                    print(f"  {i}: {tree_info['description']}")
                return
            
            if args.tree_index >= len(trees):
                print(f"Error: Tree index {args.tree_index} not found. File contains {len(trees)} tree(s).", file=sys.stderr)
                print("Use --list-trees to see available trees.", file=sys.stderr)
                sys.exit(1)
            
            # Get the specified tree
            tree_info = trees[args.tree_index]
            
            if not args.quiet:
                print(f"Visualizing tree {args.tree_index}: {tree_info['description']}")
            
            # Generate title if not provided
            title = args.title
            if not title:
                title = f"Tree {args.tree_index}: {tree_info['description']}"
            
            # Create visualization from Newick string
            results = newick_to_dendrogram(
                newick_string=tree_info['newick'],
                tree_description=tree_info['description'],
                figsize=tuple(args.figsize),
                title=title,
                output_file=args.output,
                max_label_length=args.max_label_length,
                show_plot=not args.no_show
            )
            
            if not args.quiet:
                print("Newick tree visualization completed successfully!")
                print(f"Tree has {len(results['leaf_names'])} leaves")
                if args.output:
                    print(f"Visualization saved to: {args.output}")
        
        else:
            # Handle similarity matrix files (bedtools/hammock)
            if not args.quiet:
                print("Parsing similarity matrix...")
            
            if file_format == 'bedtools':
                sim_matrix = parse_bedtools_format(args.input_file)
            else:
                sim_matrix = parse_hammock_format(args.input_file)
            
            if not args.quiet:
                print(f"Matrix size: {sim_matrix.shape[0]} x {sim_matrix.shape[1]}")
            
            # Generate title if not provided
            title = args.title
            if not title:
                input_name = Path(args.input_file).stem
                title = f'Hierarchical Clustering: {input_name} ({args.linkage} linkage)'
            
            # Create dendrogram
            if not args.quiet:
                print(f"Creating dendrogram with {args.linkage} linkage...")
            
            results = create_dendrogram(
                sim_matrix=sim_matrix,
                linkage_method=args.linkage,
                figsize=tuple(args.figsize),
                title=title,
                output_file=args.output,
                max_label_length=args.max_label_length,
                color_threshold=args.color_threshold,
                show_plot=not args.no_show
            )
            
            # Print summary if requested
            if args.summary:
                print_clustering_summary(sim_matrix, results['linkage_matrix'], args.linkage)
            
            if not args.quiet:
                print("Clustering tree generation completed successfully!")
                if args.output:
                    print(f"Dendrogram saved to: {args.output}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main() 