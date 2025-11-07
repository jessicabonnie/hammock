#!/usr/bin/env python3
"""
Benchmark C++ hammock subB strategies (hash-threshold vs mixed-stride)

This script compares:
1. Hash-threshold subsampling (default)
2. Mixed-stride subsampling (--mixed-stride)

Tests accuracy against bedtools jaccard (ground truth) and measures performance.
"""

import argparse
import subprocess
import os
import sys
import tempfile
import time
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import List, Tuple, Dict
import hashlib

# Set style for better-looking plots
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Benchmark C++ hammock subB strategies vs bedtools jaccard"
    )
    parser.add_argument(
        'bed_files',
        nargs='+',
        help='BED files to compare OR a text file containing BED file paths (one per line). '
             'All pairwise comparisons will be made.'
    )
    parser.add_argument(
        '--hammock-bin',
        default='hammock_cpp/bin/hammock',
        help='Path to hammock C++ binary (default: hammock_cpp/bin/hammock)'
    )
    parser.add_argument(
        '--precision',
        type=int,
        default=18,
        help='HyperLogLog precision (default: 18)'
    )
    parser.add_argument(
        '--runs',
        type=int,
        default=3,
        help='Number of runs for timing (default: 3)'
    )
    parser.add_argument(
        '--output-prefix',
        default='cpp_subB_benchmark',
        help='Output file prefix (default: cpp_subB_benchmark)'
    )
    parser.add_argument(
        '--cache-dir',
        default='benchmarks/cache',
        help='Directory for caching bedtools results (default: benchmarks/cache)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed for hammock (default: 42)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=None,
        help='Number of threads to use (default: auto-select by hammock)'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Print verbose output'
    )
    
    return parser.parse_args()


def ensure_cache_dir(cache_dir: str) -> Path:
    """Ensure cache directory exists."""
    cache_path = Path(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)
    return cache_path


def read_bed_files_list(file_or_files: List[str]) -> List[str]:
    """
    Read BED file paths from arguments.
    
    If a single file is provided and it's a text file (not .bed/.bed.gz),
    read paths from it (one per line). Otherwise, treat as list of BED files.
    
    Args:
        file_or_files: List of file paths from command line
        
    Returns:
        List of BED file paths
    """
    # If single argument and it looks like a text file list
    if len(file_or_files) == 1:
        filepath = file_or_files[0]
        # Check if it's a text file (not a BED file)
        if not filepath.endswith(('.bed', '.bed.gz', '.bedgraph', '.bg')):
            # Try to read as file list
            try:
                with open(filepath, 'r') as f:
                    paths = []
                    for line in f:
                        line = line.strip()
                        # Skip empty lines and comments
                        if line and not line.startswith('#'):
                            paths.append(line)
                    if paths:
                        print(f"Read {len(paths)} BED file paths from {filepath}")
                        return paths
            except (IOError, OSError):
                # If can't read as file, treat as single BED path
                pass
    
    # Otherwise, treat as direct list of BED files
    return file_or_files


def get_file_hash(filepath: str) -> str:
    """Get hash of file for caching."""
    with open(filepath, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()[:8]


def get_bedtools_cache_key(file1: str, file2: str) -> str:
    """Get cache key for bedtools result."""
    hash1 = get_file_hash(file1)
    hash2 = get_file_hash(file2)
    # Sort to ensure same key regardless of order
    sorted_hashes = tuple(sorted([hash1, hash2]))
    return f"bedtools_{sorted_hashes[0]}_{sorted_hashes[1]}.json"


def run_bedtools_jaccard(file1: str, file2: str, cache_dir: Path, verbose: bool = False) -> float:
    """Run bedtools jaccard with caching."""
    cache_file = cache_dir / get_bedtools_cache_key(file1, file2)
    
    # Check cache first
    if cache_file.exists():
        if verbose:
            print(f"  Using cached bedtools result: {cache_file}")
        with open(cache_file, 'r') as f:
            data = json.load(f)
            return data['jaccard']
    
    # Run bedtools
    if verbose:
        print(f"  Running bedtools jaccard: {file1} vs {file2}")
    
    try:
        result = subprocess.run(
            ['bedtools', 'jaccard', '-a', file1, '-b', file2],
            capture_output=True,
            text=True,
            check=True
        )
        
        # Parse output (header + data line)
        lines = result.stdout.strip().split('\n')
        if len(lines) >= 2:
            fields = lines[1].split('\t')
            jaccard = float(fields[2])  # Third column is jaccard
            
            # Cache result
            with open(cache_file, 'w') as f:
                json.dump({
                    'file1': file1,
                    'file2': file2,
                    'jaccard': jaccard,
                    'timestamp': time.time()
                }, f)
            
            return jaccard
        else:
            raise ValueError(f"Unexpected bedtools output: {result.stdout}")
            
    except subprocess.CalledProcessError as e:
        print(f"Error running bedtools: {e}")
        print(f"stderr: {e.stderr}")
        return None
    except FileNotFoundError:
        print("Error: bedtools not found. Please install bedtools.")
        return None


def run_hammock_cpp(file1: str, file2: str, hammock_bin: str, subB: float, 
                    mixed_stride: bool, precision: int, seed: int, 
                    threads: int = None, verbose: bool = False) -> Tuple[float, float, float]:
    """Run C++ hammock and return (jaccard, wall_time, cpu_time)."""
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create file lists
        files_list = os.path.join(tmpdir, 'files.txt')
        primary_list = os.path.join(tmpdir, 'primary.txt')
        
        with open(files_list, 'w') as f:
            f.write(f"{file1}\n")
        
        with open(primary_list, 'w') as f:
            f.write(f"{file2}\n")
        
        # Build command
        output_prefix = os.path.join(tmpdir, 'result')
        cmd = [
            hammock_bin,
            files_list,
            primary_list,
            '--mode', 'B',
            '--subB', str(subB),
            '--precision', str(precision),
            '--seed', str(seed),
            '-o', output_prefix
        ]
        
        if mixed_stride:
            cmd.append('--mixed-stride')
        
        if threads is not None:
            cmd.extend(['--threads', str(threads)])
        
        # Run hammock
        cpu_start = time.process_time()
        wall_start = time.time()
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        wall_time = time.time() - wall_start
        cpu_time = time.process_time() - cpu_start
        
        # Parse timing from stderr
        timing_line = [line for line in result.stderr.split('\n') if line.startswith('TIMING:')]
        if timing_line:
            # Extract total time from TIMING line
            parts = timing_line[0].split()
            for part in parts:
                if part.startswith('total='):
                    wall_time = float(part.split('=')[1])
                    break
        
        # Read result CSV
        csv_files = list(Path(tmpdir).glob('*.csv'))
        if csv_files:
            df = pd.read_csv(csv_files[0])
            jaccard = df['jaccard'].iloc[0]
            return jaccard, wall_time, cpu_time
        else:
            raise ValueError("No output CSV found")


def run_benchmark(args):
    """Run the benchmark."""
    
    # Read BED files (from list or file)
    bed_files = read_bed_files_list(args.bed_files)
    
    # Validate inputs
    if len(bed_files) < 2:
        print("Error: Need at least 2 BED files for comparison")
        sys.exit(1)
    
    # Check hammock binary exists
    if not os.path.exists(args.hammock_bin):
        print(f"Error: Hammock binary not found at {args.hammock_bin}")
        sys.exit(1)
    
    # Ensure cache directory
    cache_dir = ensure_cache_dir(args.cache_dir)
    
    # Generate all pairwise comparisons
    pairs = []
    for i, file1 in enumerate(bed_files):
        for file2 in bed_files[i+1:]:
            pairs.append((file1, file2))
    
    if args.verbose:
        print(f"Benchmarking {len(pairs)} pairwise comparisons")
        print(f"Files: {bed_files}")
        print(f"Precision: {args.precision}")
        print(f"Runs per configuration: {args.runs}")
        if args.threads is not None:
            print(f"Threads: {args.threads}")
        else:
            print(f"Threads: auto-select")
        print(f"Cache directory: {cache_dir}")
    
    # Get bedtools ground truth for all pairs
    bedtools_results = {}
    print("\n=== Getting bedtools ground truth ===")
    for file1, file2 in pairs:
        jaccard = run_bedtools_jaccard(file1, file2, cache_dir, args.verbose)
        if jaccard is None:
            print(f"Warning: Could not get bedtools result for {file1} vs {file2}")
            continue
        bedtools_results[(file1, file2)] = jaccard
        print(f"  {os.path.basename(file1)} vs {os.path.basename(file2)}: {jaccard:.6f}")
    
    # Test subB values from 0.1 to 1.0 in steps of 0.1 (excluding 0.0)
    subB_values = [round(x * 0.1, 1) for x in range(1, 11)]  # 0.1, 0.2, ..., 1.0
    
    # Results storage
    results = []
    
    print("\n=== Running hammock benchmarks ===")
    
    # Test both strategies
    strategies = [
        ('hash-threshold', False),
        ('mixed-stride', True)
    ]
    
    total_tests = len(pairs) * len(subB_values) * len(strategies) * args.runs
    test_count = 0
    
    for file1, file2 in pairs:
        pair_name = f"{os.path.basename(file1)}_vs_{os.path.basename(file2)}"
        
        if (file1, file2) not in bedtools_results:
            continue
        
        bedtools_jaccard = bedtools_results[(file1, file2)]
        
        for subB in subB_values:
            for strategy_name, mixed_stride in strategies:
                # Run multiple times for statistics
                jaccards = []
                wall_times = []
                cpu_times = []
                
                for run in range(args.runs):
                    test_count += 1
                    if args.verbose:
                        print(f"  [{test_count}/{total_tests}] {pair_name}, subB={subB}, "
                              f"strategy={strategy_name}, run={run+1}")
                    
                    try:
                        jaccard, wall_time, cpu_time = run_hammock_cpp(
                            file1, file2, args.hammock_bin, subB, mixed_stride,
                            args.precision, args.seed, args.threads, args.verbose
                        )
                        jaccards.append(jaccard)
                        wall_times.append(wall_time)
                        cpu_times.append(cpu_time)
                    except Exception as e:
                        print(f"    Error: {e}")
                        continue
                
                if jaccards and wall_times:
                    # Calculate statistics
                    mean_jaccard = np.mean(jaccards)
                    std_jaccard = np.std(jaccards)
                    mean_wall_time = np.mean(wall_times)
                    std_wall_time = np.std(wall_times)
                    mean_cpu_time = np.mean(cpu_times)
                    std_cpu_time = np.std(cpu_times)
                    
                    # Calculate error vs bedtools
                    abs_error = abs(mean_jaccard - bedtools_jaccard)
                    rel_error = abs_error / bedtools_jaccard if bedtools_jaccard > 0 else 0
                    
                    results.append({
                        'file1': os.path.basename(file1),
                        'file2': os.path.basename(file2),
                        'pair': pair_name,
                        'subB': subB,
                        'strategy': strategy_name,
                        'threads': args.threads if args.threads is not None else 'auto',
                        'hammock_jaccard_mean': mean_jaccard,
                        'hammock_jaccard_std': std_jaccard,
                        'bedtools_jaccard': bedtools_jaccard,
                        'abs_error': abs_error,
                        'rel_error_pct': rel_error * 100,
                        'wall_time_mean': mean_wall_time,
                        'wall_time_std': std_wall_time,
                        'cpu_time_mean': mean_cpu_time,
                        'cpu_time_std': std_cpu_time,
                        'runs': args.runs
                    })
                    
                    if not args.verbose:
                        # Print progress for non-verbose mode
                        print(f"  [{test_count}/{total_tests}] {strategy_name:15s} "
                              f"subB={subB:3.1f} jaccard={mean_jaccard:.4f} "
                              f"rel_err={rel_error * 100:.2f}% wall={mean_wall_time:.4f}s cpu={mean_cpu_time:.4f}s", 
                              end='\r')
    
    if not args.verbose:
        print()  # New line after progress
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Save CSV
    csv_file = f"{args.output_prefix}.csv"
    df.to_csv(csv_file, index=False)
    print(f"\n=== Results saved to {csv_file} ===")
    
    # Print summary statistics
    print("\n=== Summary Statistics ===")
    summary = df.groupby(['strategy', 'subB']).agg({
        'rel_error_pct': ['mean', 'std'],
        'wall_time_mean': ['mean', 'std'],
        'cpu_time_mean': ['mean', 'std']
    }).round(6)
    print(summary)
    
    return df


def create_plots(df: pd.DataFrame, output_prefix: str):
    """Create visualization plots."""
    
    print("\n=== Creating plots ===")
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Color palette
    colors = {'hash-threshold': '#1f77b4', 'mixed-stride': '#ff7f0e'}
    
    # Plot 1: CPU Time vs subB
    ax1 = axes[0, 0]
    for strategy in df['strategy'].unique():
        strategy_data = df[df['strategy'] == strategy]
        grouped = strategy_data.groupby('subB').agg({
            'cpu_time_mean': 'mean',
            'cpu_time_std': 'mean'
        }).reset_index()
        
        ax1.errorbar(
            grouped['subB'],
            grouped['cpu_time_mean'],
            yerr=grouped['cpu_time_std'],
            label=strategy,
            marker='o',
            capsize=5,
            color=colors[strategy]
        )
    
    ax1.set_xlabel('Subsampling Rate (subB)', fontsize=12)
    ax1.set_ylabel('CPU Time (seconds)', fontsize=12)
    ax1.set_title('Performance: CPU Time vs Subsampling Rate', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Relative Error vs subB
    ax2 = axes[0, 1]
    for strategy in df['strategy'].unique():
        strategy_data = df[df['strategy'] == strategy]
        grouped = strategy_data.groupby('subB').agg({
            'rel_error_pct': 'mean'
        }).reset_index()
        
        ax2.plot(
            grouped['subB'],
            grouped['rel_error_pct'],
            label=strategy,
            marker='o',
            color=colors[strategy]
        )
    
    ax2.set_xlabel('Subsampling Rate (subB)', fontsize=12)
    ax2.set_ylabel('Relative Error (%) vs bedtools', fontsize=12)
    ax2.set_title('Accuracy: Relative Error vs Subsampling Rate', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_yscale('log')
    
    # Plot 3: Wall Time vs subB (with error bars)
    ax3 = axes[1, 0]
    for strategy in df['strategy'].unique():
        strategy_data = df[df['strategy'] == strategy]
        grouped = strategy_data.groupby('subB').agg({
            'wall_time_mean': 'mean',
            'wall_time_std': 'mean'
        }).reset_index()
        
        ax3.errorbar(
            grouped['subB'],
            grouped['wall_time_mean'],
            yerr=grouped['wall_time_std'],
            label=strategy,
            marker='o',
            capsize=5,
            color=colors[strategy]
        )
    
    ax3.set_xlabel('Subsampling Rate (subB)', fontsize=12)
    ax3.set_ylabel('Wall Time (seconds)', fontsize=12)
    ax3.set_title('Performance: Wall Time vs Subsampling Rate', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Speedup of mixed-stride vs hash-threshold (wall time)
    ax4 = axes[1, 1]
    
    # Calculate speedup
    hash_data = df[df['strategy'] == 'hash-threshold'].groupby('subB')['wall_time_mean'].mean()
    mixed_data = df[df['strategy'] == 'mixed-stride'].groupby('subB')['wall_time_mean'].mean()
    speedup = hash_data / mixed_data
    
    ax4.plot(speedup.index, speedup.values, marker='o', color='green', linewidth=2)
    ax4.axhline(y=1.0, color='k', linestyle='--', alpha=0.5, label='No speedup')
    ax4.fill_between(speedup.index, 1.0, speedup.values, 
                     where=(speedup.values > 1.0), alpha=0.3, color='green',
                     label='Mixed-stride faster')
    ax4.fill_between(speedup.index, speedup.values, 1.0,
                     where=(speedup.values < 1.0), alpha=0.3, color='red',
                     label='Hash-threshold faster')
    
    ax4.set_xlabel('Subsampling Rate (subB)', fontsize=12)
    ax4.set_ylabel('Speedup (wall time, hash-threshold / mixed-stride)', fontsize=12)
    ax4.set_title('Performance: Mixed-stride Speedup', fontsize=14, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    plot_file = f"{output_prefix}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {plot_file}")
    
    # Also create a focused comparison plot
    fig2, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # Scatter plot: error vs time, colored by strategy, sized by subB
    for strategy in df['strategy'].unique():
        strategy_data = df[df['strategy'] == strategy]
        scatter = ax.scatter(
            strategy_data['wall_time_mean'],
            strategy_data['rel_error_pct'],
            c=strategy_data['subB'],
            s=100,
            alpha=0.6,
            cmap='viridis',
            marker='o' if strategy == 'hash-threshold' else '^',
            label=strategy,
            edgecolors='black',
            linewidth=0.5
        )
    
    ax.set_xlabel('Wall Time (seconds)', fontsize=12)
    ax.set_ylabel('Relative Error (%) vs bedtools', fontsize=12)
    ax.set_title('Accuracy vs Performance Trade-off', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Subsampling Rate (subB)', fontsize=10)
    
    plt.tight_layout()
    tradeoff_file = f"{output_prefix}_tradeoff.png"
    plt.savefig(tradeoff_file, dpi=300, bbox_inches='tight')
    print(f"Trade-off plot saved to {tradeoff_file}")
    
    plt.close('all')


def main():
    """Main function."""
    args = parse_args()
    
    # Read BED files to determine pair count for descriptive filename
    bed_files = read_bed_files_list(args.bed_files)
    if len(bed_files) < 2:
        print("Error: Need at least 2 BED files for comparison")
        sys.exit(1)
    
    # Calculate number of pairs
    num_pairs = len(bed_files) * (len(bed_files) - 1) // 2
    
    # Build descriptive output prefix
    # Format: <base>_p<precision>_npairs<num>[_t<threads>]
    parts = [args.output_prefix]
    parts.append(f"p{args.precision}")
    parts.append(f"npairs{num_pairs}")
    if args.threads is not None:
        parts.append(f"t{args.threads}")
    
    output_prefix = "_".join(parts)
    
    # Update args.output_prefix for use in run_benchmark
    args.output_prefix = output_prefix
    
    print("=" * 70)
    print("C++ Hammock subB Strategy Benchmark")
    print("=" * 70)
    
    # Run benchmark
    df = run_benchmark(args)
    
    if df.empty:
        print("Error: No results collected")
        sys.exit(1)
    
    # Create plots
    create_plots(df, args.output_prefix)
    
    print("\n" + "=" * 70)
    print("Benchmark complete!")
    print("=" * 70)
    
    # Print key findings
    print("\n=== Key Findings ===")
    
    # Find optimal subB for each strategy
    for strategy in df['strategy'].unique():
        strategy_data = df[df['strategy'] == strategy]
        
        # Best accuracy (lowest error)
        best_acc = strategy_data.loc[strategy_data['rel_error_pct'].idxmin()]
        print(f"\n{strategy}:")
        print(f"  Best accuracy: subB={best_acc['subB']}, rel_error={best_acc['rel_error_pct']:.2f}%")
        
        # Best performance (lowest wall time)
        best_perf = strategy_data.loc[strategy_data['wall_time_mean'].idxmin()]
        print(f"  Best performance: subB={best_perf['subB']}, wall_time={best_perf['wall_time_mean']:.4f}s, cpu_time={best_perf['cpu_time_mean']:.4f}s")
        
        # Best trade-off (normalized)
        # Normalize error and time to [0, 1], then find minimum sum
        norm_error = (strategy_data['rel_error_pct'] - strategy_data['rel_error_pct'].min()) / \
                     (strategy_data['rel_error_pct'].max() - strategy_data['rel_error_pct'].min() + 1e-10)
        norm_time = (strategy_data['wall_time_mean'] - strategy_data['wall_time_mean'].min()) / \
                    (strategy_data['wall_time_mean'].max() - strategy_data['wall_time_mean'].min() + 1e-10)
        tradeoff_score = norm_error + norm_time
        best_tradeoff_idx = tradeoff_score.idxmin()
        best_tradeoff = strategy_data.loc[best_tradeoff_idx]
        print(f"  Best trade-off: subB={best_tradeoff['subB']}, "
              f"rel_error={best_tradeoff['rel_error_pct']:.2f}%, "
              f"wall_time={best_tradeoff['wall_time_mean']:.4f}s, "
              f"cpu_time={best_tradeoff['cpu_time_mean']:.4f}s")


if __name__ == '__main__':
    main()

