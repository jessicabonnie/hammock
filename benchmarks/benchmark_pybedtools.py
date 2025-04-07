#!/usr/bin/env python3

import os
import time
import random
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from datetime import datetime
import tempfile
from typing import List, Tuple, Dict
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
            end = random.randint(start + 200, start + 2000)
            f.write(f"{chrom}\t{start}\t{end}\n")

def run_pybedtools_benchmark(file1_list: List[str], file2_list: List[str]) -> Tuple[float, Dict[str, float]]:
    """Run pybedtools benchmark with the given file lists."""
    start_time = time.process_time()
    
    # Sort files for pybedtools
    sorted_files1 = []
    sorted_files2 = []
    for f in file1_list:
        sorted_f = f + '.sorted.bed'
        os.system(f'sort -k1,1 -k2,2n {f} > {sorted_f}')
        sorted_files1.append(sorted_f)
    for f in file2_list:
        sorted_f = f + '.sorted.bed'
        os.system(f'sort -k1,1 -k2,2n {f} > {sorted_f}')
        sorted_files2.append(sorted_f)
    
    # Create BedTool objects for each file
    bed_files1 = [pybedtools.BedTool(f) for f in sorted_files1]
    bed_files2 = [pybedtools.BedTool(f) for f in sorted_files2]
    
    # Calculate Jaccard similarities
    similarities = {}
    for i, bed1 in enumerate(bed_files1):
        for j, bed2 in enumerate(bed_files2):
            key = f"{os.path.basename(file1_list[i])}_{os.path.basename(file2_list[j])}"
            similarities[key] = bed1.jaccard(bed2)
    
    # Clean up sorted files
    for f in sorted_files1 + sorted_files2:
        os.remove(f)
    
    end_time = time.process_time()
    return end_time - start_time, similarities

def run_hammock_benchmark(file1_list: List[str], file2_list: List[str], mode='A', sketch_type='hyperloglog', precision=12, use_rust=False) -> Tuple[float, Dict[str, Dict[str, float]]]:
    """Run hammock benchmark with specified parameters."""
    start_time = time.process_time()
    
    # Create sketches for first set of files
    sketches1 = []
    for file1 in file1_list:
        sketch = IntervalSketch.from_file(
            filename=file1,
            mode=mode,
            precision=precision,
            sketch_type=sketch_type,
            use_rust=use_rust
        )
        sketches1.append(sketch)
    
    # Create sketches for second set of files
    sketches2 = []
    for file2 in file2_list:
        sketch = IntervalSketch.from_file(
            filename=file2,
            mode=mode,
            precision=precision,
            sketch_type=sketch_type,
            use_rust=use_rust
        )
        sketches2.append(sketch)
    
    # Calculate similarities
    similarities = {}
    for i, sketch1 in enumerate(sketches1):
        for j, sketch2 in enumerate(sketches2):
            key = f"{os.path.basename(file1_list[i])}_{os.path.basename(file2_list[j])}"
            similarities[key] = sketch1.similarity_values(sketch2)
    
    end_time = time.process_time()
    return end_time - start_time, similarities

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
                time_taken, similarities = run_pybedtools_benchmark(file1_list, file2_list)
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
                            time_taken, similarities = run_hammock_benchmark(
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
                                time_taken, similarities = run_hammock_benchmark(
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
    
    # Create a figure with three subplots, one for each mode
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 15))
    
    # Colors for different implementations
    colors = {
        'pybedtools': 'black',
        'hammock_python': '#1f77b4',
        'hammock_rust': '#ff7f0e'
    }
    
    # Debug print to verify data
    print("\nVerifying data for plotting:")
    for r in results:
        print(f"\nFor {r['num_files']} files:")
        print("pybedtools:", r['pybedtools']['mean_time'])
        for mode in ['A', 'B', 'C']:
            key = f"hammock_{mode}_hyperloglog_p16"
            print(f"{key}:", r[key]['mean_time'])
            key = f"hammock_{mode}_hyperloglog_rust_p16"
            print(f"{key}:", r[key]['mean_time'])
    
    # Plot for each mode
    for mode, ax in [('A', ax1), ('B', ax2), ('C', ax3)]:
        # Plot pybedtools results
        pybedtools_mean = [r['pybedtools']['mean_time'] for r in results]
        pybedtools_std = [r['pybedtools']['std_time'] for r in results]
        ax.errorbar(num_files, pybedtools_mean, yerr=pybedtools_std,
                   fmt='o-', color=colors['pybedtools'], capsize=5,
                   label='pybedtools', linewidth=2)
        
        # Plot hammock Python implementation (precision 16)
        key = f"hammock_{mode}_hyperloglog_p16"
        mean_times = [r[key]['mean_time'] for r in results]
        std_times = [r[key]['std_time'] for r in results]
        ax.errorbar(num_files, mean_times, yerr=std_times,
                   fmt='--', color=colors['hammock_python'], capsize=5,
                   label='hammock (Python)', linewidth=2)
        
        # Plot hammock Rust implementation (precision 16)
        key = f"hammock_{mode}_hyperloglog_rust_p16"
        mean_times = [r[key]['mean_time'] for r in results]
        std_times = [r[key]['std_time'] for r in results]
        ax.errorbar(num_files, mean_times, yerr=std_times,
                   fmt='-', color=colors['hammock_rust'], capsize=5,
                   label='hammock (Rust)', linewidth=2)
        
        # Customize each subplot
        ax.set_xlabel('Number of Files')
        ax.set_ylabel('CPU Time (seconds)')
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend()
        ax.set_title(f'Mode {mode} Performance Comparison')
        
        # Use log scale for y-axis to better show differences
        ax.set_yscale('log')
        
        # Add y-axis limits to ensure consistent scale across plots
        ax.set_ylim(bottom=0.001, top=1000)  # Adjust these values based on your data
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save plot
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    plot_path = os.path.join('benchmarks/results', f'pybedtools_hammock_benchmark_{timestamp}.png')
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    
    # Save numerical results
    results_path = os.path.join('benchmarks/results', f'pybedtools_hammock_benchmark_{timestamp}.txt')
    with open(results_path, 'w') as f:
        f.write("Benchmark Results:\n")
        f.write("================\n\n")
        for r in results:
            f.write(f"Number of files: {r['num_files']}\n")
            f.write("\npybedtools:\n")
            f.write(f"  Mean CPU time: {r['pybedtools']['mean_time']:.3f} seconds\n")
            f.write(f"  Standard deviation: {r['pybedtools']['std_time']:.3f} seconds\n")
            f.write(f"  Min CPU time: {r['pybedtools']['min_time']:.3f} seconds\n")
            f.write(f"  Max CPU time: {r['pybedtools']['max_time']:.3f} seconds\n")
            f.write("\nhammock (precision 16):\n")
            for mode in ['A', 'B', 'C']:
                f.write(f"\n  Mode {mode}:\n")
                # Python implementation
                key = f"hammock_{mode}_hyperloglog_p16"
                f.write(f"    Python:\n")
                f.write(f"      Mean time: {r[key]['mean_time']:.3f} seconds\n")
                f.write(f"      Standard deviation: {r[key]['std_time']:.3f} seconds\n")
                f.write(f"      Min time: {r[key]['min_time']:.3f} seconds\n")
                f.write(f"      Max time: {r[key]['max_time']:.3f} seconds\n")
                # Rust implementation
                key = f"hammock_{mode}_hyperloglog_rust_p16"
                f.write(f"    Rust:\n")
                f.write(f"      Mean time: {r[key]['mean_time']:.3f} seconds\n")
                f.write(f"      Standard deviation: {r[key]['std_time']:.3f} seconds\n")
                f.write(f"      Min time: {r[key]['min_time']:.3f} seconds\n")
                f.write(f"      Max time: {r[key]['max_time']:.3f} seconds\n")
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