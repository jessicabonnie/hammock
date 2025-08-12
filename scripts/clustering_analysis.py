#!/usr/bin/env python3
"""
Script for performing agglomerative clustering analysis on hammock similarity matrices
and evaluating clustering quality using Normalized Mutual Information (NMI).
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from pathlib import Path
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import normalized_mutual_info_score, silhouette_score

# Add the scripts directory to Python path
scripts_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(scripts_path)

# Import utility functions
from utils import (
    parse_hammock_format,
    load_accession_key,
    detect_hammock_similarity_column,
    extract_hammock_parameters_from_filename,
    detect_file_format,
    parse_bedtools_format,
    detect_hammock_exp,
)


def perform_agglomerative_clustering(similarity_matrix, n_clusters, linkage='ward'):
    """
    Perform agglomerative clustering on a similarity matrix.
    
    Args:
        similarity_matrix (pandas.DataFrame): Similarity matrix
        n_clusters (int): Number of clusters to create
        linkage (str): Linkage method ('ward', 'complete', 'average', 'single')
        
    Returns:
        tuple: (cluster_labels, clustering_model)
    """
    # Convert similarity to distance (1 - similarity)
    distance_matrix = 1 - similarity_matrix.values
    
    # For ward linkage, we need to use euclidean metric with feature matrix
    if linkage == 'ward':
        # Use complete linkage instead of ward for precomputed distances
        # Ward linkage requires euclidean metric and doesn't work well with precomputed distances
        clustering = AgglomerativeClustering(
            n_clusters=n_clusters,
            linkage='complete',
            metric='precomputed'
        )
        cluster_labels = clustering.fit_predict(distance_matrix)
    else:
        # For other linkage methods, we can use precomputed distances
        clustering = AgglomerativeClustering(
            n_clusters=n_clusters,
            linkage=linkage,
            metric='precomputed'
        )
        cluster_labels = clustering.fit_predict(distance_matrix)
    
    return cluster_labels, clustering


def calculate_nmi(true_labels, predicted_labels):
    """
    Calculate Normalized Mutual Information between true and predicted labels.
    
    Args:
        true_labels (list): True labels
        predicted_labels (list): Predicted cluster labels
        
    Returns:
        float: Normalized Mutual Information score
    """
    return normalized_mutual_info_score(true_labels, predicted_labels)


def evaluate_clustering_with_nmi(similarity_matrix, true_labels_dict, n_clusters_range=None):
    """
    Evaluate agglomerative clustering using NMI across different numbers of clusters.
    
    Args:
        similarity_matrix (pandas.DataFrame): Similarity matrix
        true_labels_dict (dict): Mapping from file basename to true label
        n_clusters_range (list): Range of cluster numbers to test
        
    Returns:
        dict: Results with NMI scores for each number of clusters
    """
    if n_clusters_range is None:
        # Default to the range 2..30, capped by number of available samples
        upper = min(30, max(2, len(similarity_matrix)))
        n_clusters_range = range(2, upper + 1)
    
    results = {}
    
    # Get true labels for files in the similarity matrix
    true_labels = []
    for file in similarity_matrix.index:
        if file in true_labels_dict:
            true_labels.append(true_labels_dict[file])
        else:
            true_labels.append('unknown')
    
    # Test different numbers of clusters
    for n_clusters in n_clusters_range:
        try:
            cluster_labels, _ = perform_agglomerative_clustering(similarity_matrix, n_clusters)
            
            # Calculate NMI
            nmi_score = calculate_nmi(true_labels, cluster_labels)
            
            results[n_clusters] = {
                'nmi_score': nmi_score,
                'cluster_labels': cluster_labels
            }
            
        except Exception as e:
            print(f"Error with {n_clusters} clusters: {e}")
            results[n_clusters] = {
                'nmi_score': None,
                'cluster_labels': None
            }
    
    return results


def evaluate_nmi_long_table(similarity_matrix, true_labels_dict, n_clusters_range, linkage_methods):
    """
    Compute NMI for each combination of number of clusters and linkage method.

    Args:
        similarity_matrix (pandas.DataFrame): Similarity matrix filtered to files with labels
        true_labels_dict (dict): Mapping from file basename to true label
        n_clusters_range (Iterable[int]): Cluster counts to evaluate
        linkage_methods (Iterable[str]): Linkage methods to evaluate

    Returns:
        pandas.DataFrame: Long-form table with columns [n_clusters, linkage, nmi, silhouette]
    """
    # Build true labels aligned to the matrix order
    aligned_true_labels = [true_labels_dict.get(file, 'unknown') for file in similarity_matrix.index]

    rows = []
    # Precompute full distance matrix for silhouette (0 diagonal, clipped)
    distance_square = 1.0 - similarity_matrix.values
    np.fill_diagonal(distance_square, 0.0)
    distance_square = np.clip(distance_square, 0.0, 1.0)
    for linkage in linkage_methods:
        for n_clusters in n_clusters_range:
            try:
                predicted_labels, _ = perform_agglomerative_clustering(
                    similarity_matrix, n_clusters, linkage
                )
                nmi_score = calculate_nmi(aligned_true_labels, predicted_labels)
                # Compute silhouette score using precomputed distances
                silhouette = np.nan
                try:
                    silhouette = silhouette_score(distance_square, predicted_labels, metric='precomputed')
                except Exception:
                    silhouette = np.nan
                rows.append({
                    'n_clusters': int(n_clusters),
                    'linkage': str(linkage),
                    'nmi': float(nmi_score),
                    'silhouette': float(silhouette) if not pd.isna(silhouette) else np.nan,
                })
            except Exception as e:
                rows.append({
                    'n_clusters': int(n_clusters),
                    'linkage': str(linkage),
                    'nmi': np.nan,
                    'silhouette': np.nan,
                    'error': str(e)
                })

    return pd.DataFrame(rows)


def plot_nmi_results(nmi_results, title="NMI Scores vs Number of Clusters"):
    """
    Plot NMI scores against number of clusters.
    
    Args:
        nmi_results (dict): Results from evaluate_clustering_with_nmi
        title (str): Plot title
    """
    n_clusters = list(nmi_results.keys())
    nmi_scores = [results['nmi_score'] for results in nmi_results.values()]
    
    plt.figure(figsize=(10, 6))
    plt.plot(n_clusters, nmi_scores, 'bo-', linewidth=2, markersize=8)
    plt.xlabel('Number of Clusters')
    plt.ylabel('Normalized Mutual Information Score')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    return plt.gcf()


def create_confusion_matrix(similarity_matrix, true_labels_dict, n_clusters, linkage='ward'):
    """
    Create a confusion matrix comparing true labels with predicted clusters.
    
    Args:
        similarity_matrix (pandas.DataFrame): Similarity matrix
        true_labels_dict (dict): Mapping from file basename to true label
        n_clusters (int): Number of clusters
        linkage (str): Linkage method
        
    Returns:
        pandas.DataFrame: Confusion matrix
    """
    # Get true labels
    true_labels = []
    for file in similarity_matrix.index:
        if file in true_labels_dict:
            true_labels.append(true_labels_dict[file])
        else:
            true_labels.append('unknown')
    
    # Perform clustering
    cluster_labels, _ = perform_agglomerative_clustering(similarity_matrix, n_clusters, linkage)
    
    # Create confusion matrix
    confusion_df = pd.DataFrame({
        'True_Label': true_labels,
        'Predicted_Cluster': cluster_labels
    })
    
    # Create pivot table
    confusion_matrix = pd.crosstab(
        confusion_df['True_Label'], 
        confusion_df['Predicted_Cluster'], 
        margins=True
    )
    
    return confusion_matrix


def plot_confusion_matrix(confusion_matrix, title="Confusion Matrix"):
    """
    Plot confusion matrix as a heatmap.
    
    Args:
        confusion_matrix (pandas.DataFrame): Confusion matrix
        title (str): Plot title
    """
    plt.figure(figsize=(10, 8))
    sns.heatmap(confusion_matrix, annot=True, fmt='d', cmap='Blues', cbar=True)
    plt.title(title)
    plt.xlabel('Predicted Cluster')
    plt.ylabel('True Label')
    plt.tight_layout()
    plt.show()
    
    return plt.gcf()


def plot_nmi_vs_clusters_colored_by_silhouette(long_table_df: pd.DataFrame, title_prefix: str = "NMI vs #Clusters") -> list:
    """
    Scatter plot of NMI vs number of clusters, colored by silhouette score, per linkage method.

    Args:
        long_table_df (pd.DataFrame): Must contain columns 'n_clusters', 'nmi', 'silhouette', 'linkage'
        title_prefix (str): Title prefix for figures
    """
    figs = []
    methods = sorted([m for m in long_table_df['linkage'].dropna().unique()])
    for method in methods:
        sub = long_table_df[long_table_df['linkage'] == method].copy()
        # Drop rows with missing NMI
        sub = sub.dropna(subset=['nmi'])
        plt.figure(figsize=(8, 6))
        sc = plt.scatter(
            sub['n_clusters'],
            sub['nmi'],
            c=sub['silhouette'],
            cmap='viridis',
            s=60,
            edgecolors='k',
            linewidths=0.3,
        )
        cbar = plt.colorbar(sc)
        cbar.set_label('Silhouette')
        plt.xlabel('Number of Clusters')
        plt.ylabel('Normalized Mutual Information (NMI)')
        plt.title(f"{title_prefix} (linkage: {method})")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        figs.append(plt.gcf())
    return figs

def main():
    """
    Main function to perform clustering analysis.
    """
    import argparse
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Perform agglomerative clustering analysis on similarity matrices (hammock CSV or bedtools TSV)',
        epilog='Example: python clustering_analysis.py input.csv accession_key.tsv --stdout'
    )
    parser.add_argument('hammock_output', help='Path to input file (hammock CSV or bedtools pairwise TSV)')
    parser.add_argument('accession_key', help='Path to accession key TSV file')
    parser.add_argument('--clusters', type=int, nargs='+', default=range(2, 31),
                       help='Number of clusters to test (default: 2-30)')
    parser.add_argument('--linkage', choices=['ward', 'complete', 'average', 'single'], 
                       default='ward', help='Linkage method (default: ward)')
    parser.add_argument('--linkage-methods', nargs='+', choices=['ward', 'complete', 'average', 'single'],
                       default=['average', 'complete'],
                       help='Linkage methods to include in long-table output (default: average complete)')
    # Preferred: -o/--out. Keep --long-table-out for backward compatibility
    parser.add_argument('-o', '--out', dest='out', metavar='TSV',
                       help='Write long-table NMI results to this TSV file (default: <input_basename>_cluster.tsv)')
    parser.add_argument('--long-table-out', dest='long_table_out', metavar='TSV',
                       help=argparse.SUPPRESS)
    parser.add_argument('--stdout', action='store_true',
                       help='Write the long-table NMI results to standard out (TSV)')
    parser.add_argument('--plots-pdf', metavar='PDF',
                       help='Write plots (NMI vs clusters and scatter colored by silhouette) to a multi-page PDF')
    
    args = parser.parse_args()
    quiet = bool(args.stdout)
    
    # File paths from command line arguments
    hammock_output = args.hammock_output
    accession_key = args.accession_key
    
    if not quiet:
        print("Loading data...")
    
    # Detect file format and parse similarity matrix accordingly
    file_format = detect_file_format(hammock_output)
    if file_format == 'hammock':
        similarity_column = detect_hammock_similarity_column(hammock_output)
        if not quiet:
            print(f"Detected similarity column: {similarity_column}")
        similarity_matrix = parse_hammock_format(hammock_output, preferred_similarity_column=similarity_column)
    else:  # bedtools
        similarity_column = 'jaccard'
        if not quiet:
            print("Detected bedtools format; using 'jaccard' column")
        similarity_matrix = parse_bedtools_format(hammock_output)
    if not quiet:
        print(f"Loaded similarity matrix with {len(similarity_matrix)} files")
    
    # Load true labels (tissue types)
    true_labels_dict = load_accession_key(accession_key)
    if not quiet:
        print(f"Loaded {len(true_labels_dict)} true labels (tissue types)")
    
    # Check overlap between files in similarity matrix and true labels
    files_in_matrix = set(similarity_matrix.index)
    files_with_labels = set(true_labels_dict.keys())
    common_files = files_in_matrix.intersection(files_with_labels)

    # If no overlap, fall back to inferring labels from file names
    def infer_label_from_basename(basename: str) -> str:
        token = basename.split('-', 1)[0]
        if token.startswith('f') and len(token) > 1:
            token = token[1:]
        return token

    if len(common_files) == 0:
        if not quiet:
            print("Warning: No overlap between similarity matrix files and accession key. Inferring labels from filenames.")
        filtered_matrix = similarity_matrix
        effective_labels = {name: infer_label_from_basename(name) for name in filtered_matrix.index}
    else:
        if not quiet:
            print(f"Files in similarity matrix: {len(files_in_matrix)}")
            print(f"Files with tissue labels: {len(files_with_labels)}")
            print(f"Common files: {len(common_files)}")
        # Filter similarity matrix to only include files with true labels
        filtered_matrix = similarity_matrix.loc[list(common_files), list(common_files)]
        effective_labels = true_labels_dict
        if not quiet:
            print(f"Filtered matrix size: {len(filtered_matrix)} x {len(filtered_matrix)}")
    
    # Evaluate clustering with different numbers of clusters
    if not quiet:
        print("\nEvaluating clustering with NMI...")
    nmi_results = evaluate_clustering_with_nmi(
        filtered_matrix,
        effective_labels,
        n_clusters_range=args.clusters
    )

    # Compute long-table across all linkage methods and cluster counts
    long_table_df = evaluate_nmi_long_table(
        filtered_matrix,
        effective_labels,
        n_clusters_range=args.clusters,
        linkage_methods=args.linkage_methods,
    )

    # Add filename-derived parameters to long table for traceability
    filename_params = extract_hammock_parameters_from_filename(hammock_output) if file_format == 'hammock' else {
        'mode': None, 'klen': None, 'window': None, 'precision': None, 'subA': None, 'subB': None, 'expA': None
    }
    # Only include requested parameter fields in the output
    mode_token = filename_params.get('mode')
    expA_value = filename_params.get('expA')
    if expA_value is None and file_format == 'hammock' and mode_token == 'C':
        expA_value = detect_hammock_exp(hammock_output)

    long_table_df.insert(0, 'sketch', filename_params.get('sketch'))
    long_table_df.insert(1, 'mode', mode_token)
    # For mode C, include 'exp' alongside precision; otherwise, leave exp blank
    long_table_df.insert(2, 'precision', filename_params.get('precision'))
    long_table_df.insert(3, 'expA', expA_value)
    long_table_df.insert(4, 'subA', filename_params.get('subA'))
    long_table_df.insert(5, 'subB', filename_params.get('subB'))
    long_table_df.insert(6, 'window', filename_params.get('window'))
    long_table_df.insert(7, 'klen', filename_params.get('klen'))

    # Output long table (default to <input_basename>_cluster.tsv in the same directory)
    output_path = args.out or args.long_table_out
    if not output_path:
        in_path = Path(hammock_output)
        default_name = f"{in_path.stem}_cluster.tsv"
        output_path = str(in_path.with_name(default_name))
    long_table_df.sort_values(['n_clusters', 'linkage']).to_csv(
        output_path, sep='\t', index=False
    )
    if not quiet:
        print(f"Long-table NMI results written to: {output_path}")
    if args.stdout:
        # Print header and rows as TSV
        # Print header and rows including parameter columns
        ordered_cols = ['sketch', 'mode', 'precision', 'expA', 'subA', 'subB', 'window', 'klen', 'n_clusters', 'linkage', 'nmi', 'silhouette']
        print("\t".join(ordered_cols))
        for _, row in long_table_df.sort_values(['n_clusters', 'linkage']).iterrows():
            nmi_display = '' if pd.isna(row.get('nmi')) else f"{row['nmi']:.6f}"
            sil_display = '' if pd.isna(row.get('silhouette')) else f"{row['silhouette']:.6f}"
            values = [
                str(row.get('sketch') or ''),
                str(row.get('mode') or ''),
                str(row.get('precision') or ''),
                str(row.get('expA') or ''),
                str(row.get('subA') or ''),
                str(row.get('subB') or ''),
                str(row.get('window') or ''),
                str(row.get('klen') or ''),
                str(int(row['n_clusters'])),
                str(row['linkage']),
                nmi_display,
                sil_display,
            ]
            print("\t".join(values))
        # If stdout table requested, stop after printing
        return
    
    # Print results
    if not quiet:
        print("\nNMI Results:")
        print("Clusters\tNMI Score")
        print("-" * 20)
        for n_clusters, result in nmi_results.items():
            if result['nmi_score'] is not None:
                print(f"{n_clusters}\t\t{result['nmi_score']:.4f}")
            else:
                print(f"{n_clusters}\t\tError")
    
    # Find best number of clusters
    valid_results = {k: v for k, v in nmi_results.items() if v['nmi_score'] is not None}
    if valid_results:
        best_n_clusters = max(valid_results.keys(), key=lambda x: valid_results[x]['nmi_score'])
        best_nmi = valid_results[best_n_clusters]['nmi_score']
        if not quiet:
            print(f"\nBest clustering: {best_n_clusters} clusters with NMI = {best_nmi:.4f}")
        
        # Plot NMI results and NMI vs clusters colored by silhouette
        if not quiet or args.plots_pdf:
            if not quiet:
                print("\nPlotting NMI results...")
            figs = []
            figs.append(plot_nmi_results(nmi_results, "NMI Scores vs Number of Clusters"))
            try:
                if not quiet:
                    print("Plotting NMI vs #clusters colored by silhouette...")
                figs.extend(plot_nmi_vs_clusters_colored_by_silhouette(long_table_df, "NMI vs Number of Clusters"))
            except Exception as e:
                print(f"Warning: failed to plot NMI vs clusters colored by silhouette: {e}")
            if args.plots_pdf:
                with PdfPages(args.plots_pdf) as pdf:
                    for fig in figs:
                        pdf.savefig(fig)
                if not quiet:
                    print(f"Plots saved to: {args.plots_pdf}")
        
        # Confusion matrix generation is currently disabled; keep code for future use
        generate_confusion_outputs = False
        if generate_confusion_outputs:
            if not quiet:
                print(f"\nCreating confusion matrix for {best_n_clusters} clusters...")
            confusion_matrix = create_confusion_matrix(
                filtered_matrix,
                true_labels_dict,
                best_n_clusters
            )
            if not quiet:
                plot_confusion_matrix(confusion_matrix, f"Confusion Matrix ({best_n_clusters} clusters)")
        
        # Additional analysis: compare different linkage methods
        if not quiet:
            print("\nComparing different linkage methods...")
        linkage_methods = ['ward', 'complete', 'average', 'single']
        linkage_results = {}
        
        for linkage in linkage_methods:
            try:
                cluster_labels, _ = perform_agglomerative_clustering(
                    filtered_matrix, best_n_clusters, linkage
                )
                
                # Get true labels (tissue types) for files in filtered matrix
                true_labels = [true_labels_dict[file] for file in filtered_matrix.index]
                
                nmi_score = calculate_nmi(true_labels, cluster_labels)
                linkage_results[linkage] = nmi_score
                if not quiet:
                    print(f"{linkage}: NMI = {nmi_score:.4f}")
                
            except Exception as e:
                if not quiet:
                    print(f"{linkage}: Error - {e}")
                linkage_results[linkage] = None
        
        # Plot linkage comparison
        if not quiet and any(v is not None for v in linkage_results.values()):
            plt.figure(figsize=(10, 6))
            methods = list(linkage_results.keys())
            scores = [linkage_results[m] if linkage_results[m] is not None else 0 for m in methods]
            
            plt.bar(methods, scores)
            plt.ylabel('NMI Score')
            plt.title(f'Linkage Method Comparison ({best_n_clusters} clusters)')
            plt.ylim(0, 1)
            plt.tight_layout()
            plt.show()
    
    else:
        if not quiet:
            print("No valid clustering results found.")
    
    if not quiet:
        print("\nAnalysis complete!")


if __name__ == "__main__":
    main() 