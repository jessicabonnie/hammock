#!/usr/bin/env python3

import os
import time
import random
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from datetime import datetime
import tempfile
from typing import List
import pybedtools # type: ignore
from hammock.lib.intervals import IntervalSketch
from hammock.lib.sequences import SequenceSketch

# Constants for benchmarking
NUM_INTERVALS = 1000  # Fixed number of intervals per file
NUM_FILES_LIST = [2, 4, 8, 16, 32]  # Number of files to compare
NUM_RUNS = 5

def generate_bed_file(num_intervals: int, output_file: str):
    """Generate a test BED file with random intervals."""
    with open(output_file, 'w') as f:
        for i in range(num_intervals):
            # Generate random chromosome (1-22, X, Y)
            chrom = random.choice([str(i) for i in range(1, 23)] + ['X', 'Y'])
            # Generate random start position (0-1000000)
            start = random.randint(0, 1000000)
            # Generate random end position (start+1 to start+1000)
            end = random.randint(start + 1, start + 1000)
            f.write(f"{chrom}\t{start}\t{end}\n")

def run_pybedtools_benchmark(file1_list: List[str], file2_list: List[str]) -> float:
    """Run pybedtools benchmark with the given file lists."""
    start_time = time.time()
    
    # Calculate pairwise Jaccard similarities
    for file1 in file1_list:
        # Sort file for pybedtools
        sorted_file1 = file1 + '.sorted.bed'
        os.system(f'sort -k1,1 -k2,2n {file1} > {sorted_file1}')
        bed1 = pybedtools.BedTool(sorted_file1)
        for file2 in file2_list:
            # Sort file for pybedtools
            sorted_file2 = file2 + '.sorted.bed'
            os.system(f'sort -k1,1 -k2,2n {file2} > {sorted_file2}')
            bed2 = pybedtools.BedTool(sorted_file2)
            # Calculate Jaccard similarity
            jaccard = bed1.jaccard(bed2)
            # Clean up sorted files
            os.remove(sorted_file2)
        # Clean up sorted files
        os.remove(sorted_file1)
    
    end_time = time.time()
    return end_time - start_time

def run_hammock_benchmark(file1_list: List[str], file2_list: List[str], 
                        mode='A', sketch_type='hyperloglog', precision=12, 
                        use_rust=False) -> float:
    """Run hammock benchmark with specified parameters."""
    start_time = time.time()
    
    # Create sketches for first set of files
    sketches1 = {}
    for file1 in file1_list:
        sketch1 = IntervalSketch.from_file(
            filename=file1,
            mode=mode,
            precision=precision,
            sketch_type=sketch_type,
            use_rust=use_rust
        )
        sketches1[file1] = sketch1
    
    # Create sketches for second set of files
    sketches2 = {}
    for file2 in file2_list:
        sketch2 = IntervalSketch.from_file(
            filename=file2,
            mode=mode,
            precision=precision,
            sketch_type=sketch_type,
            use_rust=use_rust
        )
        sketches2[file2] = sketch2
    
    # Calculate pairwise similarities
    for sketch1 in sketches1.values():
        for sketch2 in sketches2.values():
            similarity = sketch1.similarity_values(sketch2)
    
    end_time = time.time()
    return end_time - start_time

def run_benchmark(num_files_list: List[int] = NUM_FILES_LIST, num_runs: int = NUM_RUNS):
    """Run benchmarks with different numbers of files."""
    results = []
    
    for num_files in num_files_list:
        print(f"\nBenchmarking with {num_files} files...")
        
        # Create temporary directory for test files
        with tempfile.TemporaryDirectory() as temp_dir:
            # Generate two sets of BED files
            file1_list = []
            file2_list = []
            
            for i in range(num_files):
                file1 = os.path.join(temp_dir, f"set1_file{i}.bed")
                file2 = os.path.join(temp_dir, f"set2_file{i}.bed")
                generate_bed_file(NUM_INTERVALS, file1)
                generate_bed_file(NUM_INTERVALS, file2)
                file1_list.append(file1)
                file2_list.append(file2)
            
            # Run pybedtools benchmark
            pybedtools_times = []
            for run in range(num_runs):
                time_taken = run_pybedtools_benchmark(file1_list, file2_list)
                pybedtools_times.append(time_taken)
            
            # Run hammock benchmarks with different modes and sketching options
            hammock_results = {}
            for mode in ['A', 'B', 'C']:
                for sketch_type in ['hyperloglog', 'minhash']:
                    for precision in [8, 12, 16]:
                        # Run with Python implementation
                        key = f"hammock_{mode}_{sketch_type}_p{precision}"
                        times = []
                        for run in range(num_runs):
                            time_taken = run_hammock_benchmark(
                                file1_list, file2_list,
                                mode=mode, sketch_type=sketch_type,
                                precision=precision, use_rust=False
                            )
                            times.append(time_taken)
                        hammock_results[key] = {
                            'mean_time': np.mean(times),
                            'std_time': np.std(times),
                            'min_time': np.min(times),
                            'max_time': np.max(times)
                        }
                        
                        # Run with Rust implementation for HyperLogLog
                        if sketch_type == 'hyperloglog':
                            key = f"hammock_{mode}_{sketch_type}_rust_p{precision}"
                            times = []
                            for run in range(num_runs):
                                time_taken = run_hammock_benchmark(
                                    file1_list, file2_list,
                                    mode=mode, sketch_type=sketch_type,
                                    precision=precision, use_rust=True
                                )
                                times.append(time_taken)
                            hammock_results[key] = {
                                'mean_time': np.mean(times),
                                'std_time': np.std(times),
                                'min_time': np.min(times),
                                'max_time': np.max(times)
                            }
            
            # Store results
            result = {
                'num_files': num_files,
                'pybedtools': {
                    'mean_time': np.mean(pybedtools_times),
                    'std_time': np.std(pybedtools_times),
                    'min_time': np.min(pybedtools_times),
                    'max_time': np.max(pybedtools_times)
                }
            }
            result.update(hammock_results)
            results.append(result)
    
    return results

def plot_results(results):
    """Plot the benchmark results."""
    num_files = [r['num_files'] for r in results]
    
    # Create a figure with multiple subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Plot pybedtools results
    pybedtools_mean = [r['pybedtools']['mean_time'] for r in results]
    pybedtools_std = [r['pybedtools']['std_time'] for r in results]
    ax1.errorbar(num_files, pybedtools_mean, yerr=pybedtools_std, 
                fmt='o-', capsize=5, label='pybedtools')
    
    # Plot hammock results for different modes and sketching options
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    markers = ['o', 's', '^', 'D', 'v', '<']
    
    # Plot Python implementation results
    for idx, (mode, sketch_type) in enumerate([('A', 'hyperloglog'), ('B', 'hyperloglog'),
                                              ('C', 'hyperloglog'), ('A', 'minhash'),
                                              ('B', 'minhash'), ('C', 'minhash')]):
        for precision in [8, 12, 16]:
            key = f"hammock_{mode}_{sketch_type}_p{precision}"
            mean_times = [r[key]['mean_time'] for r in results]
            std_times = [r[key]['std_time'] for r in results]
            label = f"{key} (p={precision})"
            ax2.errorbar(num_files, mean_times, yerr=std_times,
                        fmt=f'{markers[idx]}-', capsize=5, color=colors[idx],
                        label=label)
    
    # Plot Rust implementation results
    for idx, mode in enumerate(['A', 'B', 'C']):
        for precision in [8, 12, 16]:
            key = f"hammock_{mode}_hyperloglog_rust_p{precision}"
            mean_times = [r[key]['mean_time'] for r in results]
            std_times = [r[key]['std_time'] for r in results]
            label = f"{key} (p={precision})"
            ax2.errorbar(num_files, mean_times, yerr=std_times,
                        fmt=f'{markers[idx]}--', capsize=5, color=colors[idx],
                        label=label)
    
    # Customize plots
    for ax in [ax1, ax2]:
        ax.set_xlabel('Number of Files')
        ax.set_ylabel('Execution Time (seconds)')
        ax.grid(True)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    ax1.set_title('pybedtools Performance')
    ax2.set_title('hammock Performance (solid=Python, dashed=Rust)')
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save plot
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    plot_path = os.path.join('benchmarks/results', f'pybedtools_hammock_benchmark_{timestamp}.png')
    plt.savefig(plot_path, bbox_inches='tight')
    plt.close()
    
    # Save numerical results
    results_path = os.path.join('benchmarks/results', f'pybedtools_hammock_benchmark_{timestamp}.txt')
    with open(results_path, 'w') as f:
        f.write("Benchmark Results:\n")
        f.write("================\n\n")
        for r in results:
            f.write(f"Number of files: {r['num_files']}\n")
            f.write("\npybedtools:\n")
            f.write(f"  Mean time: {r['pybedtools']['mean_time']:.3f} seconds\n")
            f.write(f"  Standard deviation: {r['pybedtools']['std_time']:.3f} seconds\n")
            f.write(f"  Min time: {r['pybedtools']['min_time']:.3f} seconds\n")
            f.write(f"  Max time: {r['pybedtools']['max_time']:.3f} seconds\n")
            f.write("\nhammock:\n")
            for key in r:
                if key.startswith('hammock_'):
                    f.write(f"\n  {key}:\n")
                    f.write(f"    Mean time: {r[key]['mean_time']:.3f} seconds\n")
                    f.write(f"    Standard deviation: {r[key]['std_time']:.3f} seconds\n")
                    f.write(f"    Min time: {r[key]['min_time']:.3f} seconds\n")
                    f.write(f"    Max time: {r[key]['max_time']:.3f} seconds\n")
            f.write("\n" + "="*50 + "\n\n")

def main():
    # Create results directory if it doesn't exist
    os.makedirs('benchmarks/results', exist_ok=True)
    
    print("Starting pybedtools and hammock benchmark...")
    results = run_benchmark()
    plot_results(results)
    print("\nBenchmark completed. Results saved in benchmarks/results/")

if __name__ == "__main__":
    main() 