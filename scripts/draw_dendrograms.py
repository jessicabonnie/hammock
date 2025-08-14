#!/usr/bin/env python3
"""
Draw dendrograms for a similarity matrix (hammock CSV or bedtools TSV), one page per
hierarchical clustering linkage method, with leaf labels colored by true tissue labels.

Usage examples:
  - python3 scripts/draw_dendrograms.py <input_matrix.csv|tsv> <accession_key.tsv>
  - python3 scripts/draw_dendrograms.py input.csv key.tsv --output results/dendrograms.pdf
  - python3 scripts/draw_dendrograms.py input.csv key.tsv --linkage-methods average complete single
"""

import os
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# Ensure we can import shared utilities
scripts_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(scripts_path)

from utils import (
    parse_hammock_format,
    parse_bedtools_format,
    detect_file_format,
    load_accession_key,
)


def convert_similarity_to_condensed_distance(similarity_matrix: pd.DataFrame) -> np.ndarray:
    """
    Convert a similarity matrix to a condensed distance matrix suitable for scipy.linkage.

    Distance = 1 - similarity, diagonal forced to 0, values clipped to [0, 1].
    """
    distance_square = 1.0 - similarity_matrix.values.astype(float)
    np.fill_diagonal(distance_square, 0.0)
    distance_square = np.clip(distance_square, 0.0, 1.0)
    return squareform(distance_square)


def infer_label_from_basename(basename: str) -> str:
    """
    Infer a label from a file basename if not present in the accession key.
    Heuristic: take token before the first '-', drop leading 'f' if present.
    """
    token = str(basename).split('-', 1)[0]
    if token.startswith('f') and len(token) > 1:
        token = token[1:]
    return token


def align_matrix_and_labels(
    similarity_matrix: pd.DataFrame,
    accession_labels: Dict[str, str],
    quiet: bool = False,
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """
    Subset the similarity matrix to files present in the accession key. If there is no
    overlap, return the original matrix and infer labels from basenames.
    """
    files_in_matrix = set(similarity_matrix.index)
    files_with_labels = set(accession_labels.keys())
    common_files = sorted(files_in_matrix.intersection(files_with_labels))

    if len(common_files) == 0:
        if not quiet:
            print("Warning: No overlap with accession key; inferring labels from basenames.")
        inferred = {name: infer_label_from_basename(name) for name in similarity_matrix.index}
        return similarity_matrix, inferred

    filtered = similarity_matrix.loc[common_files, common_files]
    return filtered, accession_labels


def build_label_color_map(labels: List[str]) -> Dict[str, str]:
    """
    Assign a distinct color per unique label using a categorical palette.
    """
    unique_labels = sorted(set(labels))
    # Use matplotlib tab20 palette (repeats if more than 20)
    base_colors = plt.get_cmap('tab20').colors
    color_map: Dict[str, str] = {}
    for idx, lab in enumerate(unique_labels):
        color = base_colors[idx % len(base_colors)]
        # Convert RGBA tuple to hex for consistent use
        color_map[lab] = '#%02x%02x%02x' % (
            int(color[0] * 255), int(color[1] * 255), int(color[2] * 255)
        )
    return color_map


def plot_dendrogram_colored(
    ax: plt.Axes,
    condensed_distance: np.ndarray,
    sample_names: List[str],
    name_to_label: Dict[str, str],
    method: str,
    title: str,
) -> None:
    """
    Plot a dendrogram on the given axes with leaf labels colored according to name_to_label.
    """
    # Ward is not compatible with condensed distance input; map to complete with a note
    effective_method = method
    if method.lower() == 'ward':
        effective_method = 'complete'

    Z = linkage(condensed_distance, method=effective_method)

    # Build mapping from cluster ids to tissue-uniform colors for branch coloring
    n = len(sample_names)
    sample_labels_by_index = [name_to_label.get(name, 'unknown') for name in sample_names]
    color_map = build_label_color_map(sample_labels_by_index)

    # Prepare label sets per node id (leaves: 0..n-1; internal: n..n+(n-2))
    label_sets: Dict[int, set] = {i: {sample_labels_by_index[i]} for i in range(n)}
    link_color_map: Dict[int, str] = {}
    for i in range(Z.shape[0]):
        left = int(Z[i, 0])
        right = int(Z[i, 1])
        curr_id = n + i
        combined = label_sets[left] | label_sets[right]
        label_sets[curr_id] = combined
        if len(combined) == 1:
            only_label = next(iter(combined))
            link_color_map[curr_id] = color_map.get(only_label, '#444444')
        else:
            link_color_map[curr_id] = '#888888'

    def link_color_func(k: int) -> str:
        return link_color_map.get(k, '#888888')

    dendro = dendrogram(
        Z,
        labels=sample_names,
        leaf_rotation=90,
        leaf_font_size=7,
        color_threshold=0,  # disable distance threshold coloring
        above_threshold_color='#888888',
        link_color_func=link_color_func,
        ax=ax,
        no_labels=False,
    )

    ax.set_title(title, fontsize=12)
    ax.set_ylabel('Distance (1 - similarity)')

    # Also color tick labels by tissue
    for tick_label in ax.get_xmajorticklabels():
        sample_name = tick_label.get_text()
        tissue = name_to_label.get(sample_name, 'unknown')
        tick_label.set_color(color_map.get(tissue, '#444444'))

    # Build legend for tissues present
    labels_in_plot = dendro.get('ivl', [])
    label_values = [name_to_label.get(name, 'unknown') for name in labels_in_plot]
    unique_tissues: List[str] = []
    for tissue in label_values:
        if tissue not in unique_tissues:
            unique_tissues.append(tissue)
    handles = [
        plt.Line2D([0], [0], marker='o', color='w', label=t,
                   markerfacecolor=color_map.get(t, '#444444'), markersize=6)
        for t in unique_tissues
    ]
    ax.legend(handles=handles, title='Tissue', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Draw dendrograms (one per linkage) from a similarity matrix (hammock CSV or bedtools TSV), colored by tissue labels',
        epilog='Example: python3 scripts/draw_dendrograms.py input.csv accession_key.tsv --output dendrograms.pdf'
    )
    parser.add_argument('input_matrix', help='Path to hammock CSV output file or bedtools pairwise TSV')
    parser.add_argument('accession_key', help='Path to accession key TSV file')
    parser.add_argument('--linkage-methods', nargs='+', choices=['single', 'complete', 'average', 'ward'],
                        default=['average', 'complete', 'single'],
                        help='Linkage methods to render (default: average complete single)')
    parser.add_argument('--output', '-o', help='Output PDF path (default: <input_basename>_dendrograms.pdf)')
    parser.add_argument('--figsize', nargs=2, type=float, metavar=('W', 'H'), default=[12.0, 8.0],
                        help='Figure size in inches (default: 12 8)')
    parser.add_argument('--quiet', action='store_true', help='Suppress status prints')

    args = parser.parse_args()

    input_path = args.input_matrix
    key_path = args.accession_key
    methods = args.linkage_methods
    figsize = tuple(args.figsize)
    quiet = bool(args.quiet)

    if not Path(input_path).exists():
        print(f"Error: input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    if not Path(key_path).exists():
        print(f"Error: accession key not found: {key_path}", file=sys.stderr)
        sys.exit(1)

    if args.output:
        output_pdf = args.output
    else:
        # Default: same basename (including extension) with "_dend.pdf" appended
        base_with_ext = Path(input_path).name
        output_pdf = str(Path(input_path).parent / f"{base_with_ext}_dend.pdf")

    # Parse input matrix
    file_format = detect_file_format(input_path)
    if not quiet:
        print(f"Detected input format: {file_format}")

    if file_format == 'hammock':
        sim_matrix = parse_hammock_format(input_path)
    else:
        sim_matrix = parse_bedtools_format(input_path)

    if not quiet:
        print(f"Loaded similarity matrix with {len(sim_matrix)} samples")

    # Load labels and align
    accession_labels = load_accession_key(key_path)
    sim_matrix_aligned, label_map = align_matrix_and_labels(sim_matrix, accession_labels, quiet=quiet)

    if sim_matrix_aligned.shape[0] < 2:
        print("Error: need at least 2 samples to draw a dendrogram", file=sys.stderr)
        sys.exit(1)

    # Convert to condensed distance
    condensed = convert_similarity_to_condensed_distance(sim_matrix_aligned)
    names = list(sim_matrix_aligned.index)

    # Render one dendrogram per method into a single PDF
    if not quiet:
        print(f"Writing dendrograms to: {output_pdf}")

    with PdfPages(output_pdf) as pdf:
        for method in methods:
            fig, ax = plt.subplots(figsize=figsize)
            title = f"Hierarchical Clustering Dendrogram ({method})"
            if method.lower() == 'ward':
                title += " [ward uses 'complete' on precomputed distances]"
            plot_dendrogram_colored(
                ax=ax,
                condensed_distance=condensed,
                sample_names=names,
                name_to_label=label_map,
                method=method,
                title=title,
            )
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    if not quiet:
        print("Done.")


if __name__ == '__main__':
    main()


