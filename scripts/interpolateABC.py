#!/usr/bin/env python3
"""
Script to calculate Robinson-Foulds distances between hammock output trees.

Compares mode C output trees (precision 24) against:
1. Mode A output trees (precision 24) 
2. Mode B output trees (precision 24)

For each comparison, tracks the expA value and outputs:
File1, File2, expA, RFA, RFB

Also generates a line graph showing RF distances vs expA values.
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
import warnings
import matplotlib.pyplot as plt
import seaborn as sns

# Add the scripts directory to Python path
scripts_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(scripts_path)

# Import utility functions
from utils import parse_hammock_format, extract_hammock_parameters_from_filename

# Suppress scipy warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)


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
        rf = max_rf // 3  # Conservative estimate
        normalized_rf = rf / max_rf
    else:
        rf = 0
        max_rf = 0
        normalized_rf = 0.0
    
    print("Warning: Using simplified RF distance approximation. Install ete3 or dendropy for accurate results.", file=sys.stderr)
    
    return rf, max_rf, normalized_rf, n_leaves


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


def build_tree_from_similarity_matrix(sim_matrix, linkage_method='average'):
    """
    Build a hierarchical clustering tree from a similarity matrix.
    
    Args:
        sim_matrix (pandas.DataFrame): Similarity matrix
        linkage_method (str): Linkage method for clustering
        
    Returns:
        str: Newick format tree string
    """
    # Convert similarity to distance matrix
    condensed_dist = similarity_to_distance_matrix(sim_matrix)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_dist, method=linkage_method)
    
    # Convert to tree
    tree_node = to_tree(linkage_matrix)
    
    # Convert to Newick format
    names = list(sim_matrix.index)
    newick_tree = cluster_tree_to_newick(tree_node, names)
    
    return newick_tree


def find_matching_files(directory, mode, precision, expA=None):
    """
    Find hammock output files matching specific parameters.
    
    Args:
        directory (str): Directory to search
        mode (str): Mode (A, B, or C)
        precision (int): Precision value
        expA (float, optional): expA value for mode C files
        
    Returns:
        list: List of matching file paths
    """
    matching_files = []
    
    for file_path in Path(directory).glob("*.csv"):
        try:
            params = extract_hammock_parameters_from_filename(str(file_path))
            
            # Check if mode and precision match
            if params.get('mode') == mode and params.get('precision') == precision:
                # For mode C, optionally filter by expA
                if mode == 'C' and expA is not None:
                    if params.get('expA') == expA:
                        matching_files.append(str(file_path))
                else:
                    matching_files.append(str(file_path))
                    
        except Exception as e:
            print(f"Warning: Could not parse parameters from {file_path}: {e}", file=sys.stderr)
            continue
    
    return matching_files


def create_two_page_plot(results_expa, results_subb, output_basename):
    """
    Create a two-page PDF with separate plots for expA and subB parameters.
    
    Args:
        results_expa (list): List of dictionaries with expA results
        results_subb (list): List of dictionaries with subB results
        output_basename (str): Base name for output files (without extension)
    """
    # Set up the plot style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with two pages
    fig = plt.figure(figsize=(12, 8))
    
    # Page 1: expA plot
    if results_expa:
        ax1 = plt.subplot(111)
        df_expa = pd.DataFrame(results_expa)
        df_sorted = df_expa.sort_values('expA')
        
        # Plot baseline (A vs B distance)
        ax1.plot(df_sorted['expA'], df_sorted['RFA'], 
                marker='o', linewidth=2, markersize=6, 
                label='Baseline: Mode A vs Mode B', color='blue')
        
        # Add a horizontal line to show the baseline level
        baseline_value = df_sorted['RFA'].iloc[0]  # All values should be the same
        ax1.axhline(y=baseline_value, color='gray', linestyle='--', alpha=0.7, 
                   label=f'Baseline Level: {baseline_value:.3f}')
        
        # Customize the plot
        ax1.set_xlabel('expA', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Normalized Robinson-Foulds Distance', fontsize=12, fontweight='bold')
        ax1.set_title('Baseline Normalized RF Distance vs expA: Mode A vs Mode B', 
                     fontsize=14, fontweight='bold', pad=20)
        
        # Add grid
        ax1.grid(True, alpha=0.3)
        
        # Add legend
        ax1.legend(fontsize=11, frameon=True, fancybox=True, shadow=True)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save the first page
        plot_filename = f"{output_basename}_expA.png"
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        print(f"expA plot saved as: {plot_filename}")
        
        plt.close()
    
    # Page 2: subB plot
    if results_subb:
        fig2 = plt.figure(figsize=(12, 8))
        ax2 = plt.subplot(111)
        df_subb = pd.DataFrame(results_subb)
        df_sorted = df_subb.sort_values('subB')
        
        # Plot baseline (A vs B distance)
        ax2.plot(df_sorted['subB'], df_sorted['RFA'], 
                marker='s', linewidth=2, markersize=6, 
                label='Baseline: Mode A vs Mode B', color='red')
        
        # Add a horizontal line to show the baseline level
        baseline_value = df_sorted['RFA'].iloc[0]  # All values should be the same
        ax2.axhline(y=baseline_value, color='gray', linestyle='--', alpha=0.7, 
                   label=f'Baseline Level: {baseline_value:.3f}')
        
        # Customize the plot
        ax2.set_xlabel('subB', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Normalized Robinson-Foulds Distance', fontsize=12, fontweight='bold')
        ax2.set_title('Baseline Normalized RF Distance vs subB: Mode A vs Mode B', 
                     fontsize=14, fontweight='bold', pad=20)
        
        # Add grid
        ax2.grid(True, alpha=0.3)
        
        # Add legend
        ax2.legend(fontsize=11, frameon=True, fancybox=True, shadow=True)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save the second page
        plot_filename = f"{output_basename}_subB.png"
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        print(f"subB plot saved as: {plot_filename}")
        
        plt.close()
    
    # Create combined PDF with both pages
    if results_expa or results_subb:
        from matplotlib.backends.backend_pdf import PdfPages
        
        pdf_filename = f"{output_basename}.pdf"
        with PdfPages(pdf_filename) as pdf:
            # Page 1: expA
            if results_expa:
                fig1 = plt.figure(figsize=(12, 8))
                ax1 = plt.subplot(111)
                df_expa = pd.DataFrame(results_expa)
                df_sorted = df_expa.sort_values('expA')
                
                ax1.plot(df_sorted['expA'], df_sorted['RFA'], 
                        marker='o', linewidth=2, markersize=6, 
                        label='Baseline: Mode A vs Mode B', color='blue')
                
                baseline_value = df_sorted['RFA'].iloc[0]
                ax1.axhline(y=baseline_value, color='gray', linestyle='--', alpha=0.7, 
                           label=f'Baseline Level: {baseline_value:.3f}')
                
                ax1.set_xlabel('expA', fontsize=12, fontweight='bold')
                ax1.set_ylabel('Normalized Robinson-Foulds Distance', fontsize=12, fontweight='bold')
                ax1.set_title('Baseline Normalized RF Distance vs expA: Mode A vs Mode B', 
                             fontsize=14, fontweight='bold', pad=20)
                ax1.grid(True, alpha=0.3)
                ax1.legend(fontsize=11, frameon=True, fancybox=True, shadow=True)
                plt.tight_layout()
                pdf.savefig(fig1)
                plt.close()
            
            # Page 2: subB
            if results_subb:
                fig2 = plt.figure(figsize=(12, 8))
                ax2 = plt.subplot(111)
                df_subb = pd.DataFrame(results_subb)
                df_sorted = df_subb.sort_values('subB')
                
                ax2.plot(df_sorted['subB'], df_sorted['RFA'], 
                        marker='s', linewidth=2, markersize=6, 
                        label='Baseline: Mode A vs Mode B', color='red')
                
                baseline_value = df_sorted['RFA'].iloc[0]
                ax2.axhline(y=baseline_value, color='gray', linestyle='--', alpha=0.7, 
                           label=f'Baseline Level: {baseline_value:.3f}')
                
                ax2.set_xlabel('subB', fontsize=12, fontweight='bold')
                ax2.set_ylabel('Normalized Robinson-Foulds Distance', fontsize=12, fontweight='bold')
                ax2.set_title('Baseline Normalized RF Distance vs subB: Mode A vs Mode B', 
                             fontsize=14, fontweight='bold', pad=20)
                ax2.grid(True, alpha=0.3)
                ax2.legend(fontsize=11, frameon=True, fancybox=True, shadow=True)
                plt.tight_layout()
                pdf.savefig(fig2)
                plt.close()
        
        print(f"Combined PDF saved as: {pdf_filename}")


def main():
    parser = argparse.ArgumentParser(
        description='Calculate Robinson-Foulds distances between hammock output trees and create interpolation plots',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python interpolateABC.py /path/to/hammock/outputs --output ABC_interpolation.tsv
  python interpolateABC.py /path/to/hammock/outputs --precision 24 --linkage average
  python interpolateABC.py /path/to/hammock/outputs --output-tag test_run
  python interpolateABC.py /path/to/hammock/outputs --output my_results.tsv --output-tag experiment1
        """
    )
    
    parser.add_argument('directory', help='Directory containing hammock output files')
    parser.add_argument('--output', '-o', default='ABC_interpolation.tsv', 
                       help='Output file for results (default: ABC_interpolation.tsv)')
    parser.add_argument('--output-tag', '-t', 
                       help='Optional tag to add to output basename before extension')
    parser.add_argument('--precision', '-p', type=int, default=24,
                       help='Precision value to filter files (default: 24)')
    parser.add_argument('--linkage', '-l', choices=['single', 'complete', 'average', 'ward'],
                       default='average', help='Linkage method for clustering (default: average)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.directory):
        print(f"Error: Directory {args.directory} does not exist", file=sys.stderr)
        sys.exit(1)
    
    print(f"Searching for hammock output files in {args.directory}")
    print(f"Precision: {args.precision}, Linkage method: {args.linkage}")
    
    # Find all mode C files with the specified precision
    mode_c_files = find_matching_files(args.directory, 'C', args.precision)
    
    if not mode_c_files:
        print(f"No mode C files found with precision {args.precision}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(mode_c_files)} mode C files")
    
    # Results storage for expA and subB
    results_expa = []
    results_subb = []
    
    for mode_c_file in mode_c_files:
        if args.verbose:
            print(f"\nProcessing mode C file: {os.path.basename(mode_c_file)}")
        
        # Extract parameters from mode C file
        try:
            params = extract_hammock_parameters_from_filename(mode_c_file)
            expA = params.get('expA')
            subB = params.get('subB')
            
            # Skip if neither expA nor subB is present
            if expA is None and subB is None:
                print(f"Warning: Could not extract expA or subB from {mode_c_file}, skipping", file=sys.stderr)
                continue
                
        except Exception as e:
            print(f"Warning: Could not parse parameters from {mode_c_file}: {e}", file=sys.stderr)
            continue
        
        # Find corresponding mode A and mode B files
        mode_a_files = find_matching_files(args.directory, 'A', args.precision)
        mode_b_files = find_matching_files(args.directory, 'B', args.precision)
        
        if not mode_a_files:
            print(f"Warning: No mode A files found with precision {args.precision}", file=sys.stderr)
            continue
            
        if not mode_b_files:
            print(f"Warning: No mode B files found with precision {args.precision}", file=sys.stderr)
            continue
        
        # For now, use the first matching file of each mode
        # In a more sophisticated version, you might want to match by other parameters
        mode_a_file = mode_a_files[0]
        mode_b_file = mode_b_files[0]
        
        if args.verbose:
            print(f"  Mode A file: {os.path.basename(mode_a_file)}")
            print(f"  Mode B file: {os.path.basename(mode_b_file)}")
            if expA is not None:
                print(f"  expA: {expA}")
            if subB is not None:
                print(f"  subB: {subB}")
        
        try:
            # Parse similarity matrices
            sim_matrix_c = parse_hammock_format(mode_c_file)
            sim_matrix_a = parse_hammock_format(mode_a_file)
            sim_matrix_b = parse_hammock_format(mode_b_file)
            
            # Build trees
            tree_c = build_tree_from_similarity_matrix(sim_matrix_c, args.linkage)
            tree_a = build_tree_from_similarity_matrix(sim_matrix_a, args.linkage)
            tree_b = build_tree_from_similarity_matrix(sim_matrix_b, args.linkage)
            
            # Calculate RF distances
            rf_c_vs_a, max_rf_a, norm_rf_a, n_leaves_a = calculate_robinson_foulds_distance(tree_c, tree_a)
            rf_c_vs_b, max_rf_b, norm_rf_b, n_leaves_b = calculate_robinson_foulds_distance(tree_c, tree_b)
            rf_a_vs_b, max_rf_ab, norm_rf_ab, n_leaves_ab = calculate_robinson_foulds_distance(tree_a, tree_b)
            
            # Store results based on parameter type
            result_entry = {
                'File': os.path.basename(mode_c_file),
                'RFA': norm_rf_ab,  # Normalized RF distance between mode A and mode B
                'RFB': norm_rf_ab   # Normalized RF distance between mode A and mode B
            }
            
            if expA is not None:
                result_entry['expA'] = expA
                results_expa.append(result_entry.copy())
                
            if subB is not None:
                result_entry['subB'] = subB
                results_subb.append(result_entry.copy())
            
            if args.verbose:
                print(f"  RF distance C vs A: {rf_c_vs_a} (normalized: {norm_rf_a:.4f})")
                print(f"  RF distance C vs B: {rf_c_vs_b} (normalized: {norm_rf_b:.4f})")
                print(f"  RF distance A vs B: {rf_a_vs_b} (normalized: {norm_rf_ab:.4f})")
                print(f"  Using normalized RF distance: {norm_rf_ab:.4f}")
                
        except Exception as e:
            print(f"Error processing {mode_c_file}: {e}", file=sys.stderr)
            continue
    
    # Create results DataFrames and save
    if results_expa or results_subb:
        
        # Handle output tag if provided
        if args.output_tag:
            output_path = Path(args.output)
            # Insert tag before extension
            new_basename = f"{output_path.stem}_{args.output_tag}{output_path.suffix}"
            output_file = output_path.parent / new_basename
        else:
            output_file = args.output
            
        # Save expA results if any
        if results_expa:
            df_expa = pd.DataFrame(results_expa)
            expa_output = output_file.parent / f"{Path(output_file).stem}_expA{Path(output_file).suffix}"
            df_expa.to_csv(expa_output, sep='\t', index=False)
            print(f"\nexpA results saved to {expa_output}")
            print(f"Processed {len(results_expa)} expA comparisons")
            
        # Save subB results if any
        if results_subb:
            df_subb = pd.DataFrame(results_subb)
            subb_output = output_file.parent / f"{Path(output_file).stem}_subB{Path(output_file).suffix}"
            df_subb.to_csv(subb_output, sep='\t', index=False)
            print(f"subB results saved to {subb_output}")
            print(f"Processed {len(results_subb)} subB comparisons")
        
        # Create the two-page plot
        output_basename = Path(output_file).stem
        create_two_page_plot(results_expa, results_subb, output_basename)
        
        # Print summary
        print("\nSummary:")
        if results_expa:
            df_expa = pd.DataFrame(results_expa)
            print(f"expA results: {len(results_expa)} comparisons")
            print(f"  Average normalized RF distance A vs B: {df_expa['RFA'].mean():.4f}")
            print(f"  expA range: {df_expa['expA'].min()} to {df_expa['expA'].max()}")
        if results_subb:
            df_subb = pd.DataFrame(results_subb)
            print(f"subB results: {len(results_subb)} comparisons")
            print(f"  Average normalized RF distance A vs B: {df_subb['RFA'].mean():.4f}")
            print(f"  subB range: {df_subb['subB'].min()} to {df_subb['subB'].max()}")
        print(f"Note: RFA and RFB columns both contain the normalized RF distance between mode A and mode B files")
        
    else:
        print("No valid comparisons found", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
