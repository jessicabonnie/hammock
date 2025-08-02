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
import seaborn as sns
from pathlib import Path
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import normalized_mutual_info_score

# Add the scripts directory to Python path
scripts_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(scripts_path)

# Import utility functions
from utils import parse_hammock_format, load_accession_key


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
        n_clusters_range = range(2, min(11, len(similarity_matrix) // 2))
    
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
    plt.show()
    
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


def main():
    """
    Main function to perform clustering analysis.
    """
    import argparse
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Perform agglomerative clustering analysis on hammock similarity matrices',
        epilog='Example: python clustering_analysis.py hammock_output.csv accession_key.tsv'
    )
    parser.add_argument('hammock_output', help='Path to hammock CSV output file')
    parser.add_argument('accession_key', help='Path to accession key TSV file')
    parser.add_argument('--clusters', type=int, nargs='+', default=range(2, 11),
                       help='Number of clusters to test (default: 2-10)')
    parser.add_argument('--linkage', choices=['ward', 'complete', 'average', 'single'], 
                       default='ward', help='Linkage method (default: ward)')
    
    args = parser.parse_args()
    
    # File paths from command line arguments
    hammock_output = args.hammock_output
    accession_key = args.accession_key
    
    print("Loading data...")
    
    # Load hammock similarity matrix
    similarity_matrix = parse_hammock_format(hammock_output)
    print(f"Loaded similarity matrix with {len(similarity_matrix)} files")
    
    # Load true labels (tissue types)
    true_labels_dict = load_accession_key(accession_key)
    print(f"Loaded {len(true_labels_dict)} true labels (tissue types)")
    
    # Check overlap between files in similarity matrix and true labels
    files_in_matrix = set(similarity_matrix.index)
    files_with_labels = set(true_labels_dict.keys())
    common_files = files_in_matrix.intersection(files_with_labels)
    
    print(f"Files in similarity matrix: {len(files_in_matrix)}")
    print(f"Files with tissue labels: {len(files_with_labels)}")
    print(f"Common files: {len(common_files)}")
    
    # Filter similarity matrix to only include files with true labels
    filtered_matrix = similarity_matrix.loc[list(common_files), list(common_files)]
    print(f"Filtered matrix size: {len(filtered_matrix)} x {len(filtered_matrix)}")
    
    # Evaluate clustering with different numbers of clusters
    print("\nEvaluating clustering with NMI...")
    nmi_results = evaluate_clustering_with_nmi(
        filtered_matrix, 
        true_labels_dict,
        n_clusters_range=args.clusters
    )
    
    # Print results
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
        print(f"\nBest clustering: {best_n_clusters} clusters with NMI = {best_nmi:.4f}")
        
        # Plot NMI results
        print("\nPlotting NMI results...")
        plot_nmi_results(nmi_results, "NMI Scores vs Number of Clusters")
        
        # Create confusion matrix for best clustering
        print(f"\nCreating confusion matrix for {best_n_clusters} clusters...")
        confusion_matrix = create_confusion_matrix(
            filtered_matrix, 
            true_labels_dict, 
            best_n_clusters
        )
        
        print("\nConfusion Matrix:")
        print(confusion_matrix)
        
        # Plot confusion matrix
        plot_confusion_matrix(confusion_matrix, f"Confusion Matrix ({best_n_clusters} clusters)")
        
        # Additional analysis: compare different linkage methods
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
                print(f"{linkage}: NMI = {nmi_score:.4f}")
                
            except Exception as e:
                print(f"{linkage}: Error - {e}")
                linkage_results[linkage] = None
        
        # Plot linkage comparison
        if any(v is not None for v in linkage_results.values()):
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
        print("No valid clustering results found.")
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main() 