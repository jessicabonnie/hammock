#!/usr/bin/env python3

import os
import time
import random
import subprocess
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from datetime import datetime
import tempfile
from typing import List

# Get the absolute path to the bedtools.sh script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BEDTOOLS_SCRIPT = os.path.join(SCRIPT_DIR, 'bedtools.sh')

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

def run_bedtools_benchmark(file1_list_path: str, file2_list_path: str) -> float:
    """Run bedtools benchmark with the given file lists."""
    start_time = time.process_time()
    subprocess.run(['bash', BEDTOOLS_SCRIPT, file1_list_path, file2_list_path], 
                  check=True, capture_output=True)
    end_time = time.process_time()
    return end_time - start_time

def run_hammock_benchmark(file1_list_path: str, file2_list_path: str, mode='A', sketch_type='hyperloglog', precision=12, use_rust=False) -> float:
    """Run hammock benchmark with specified parameters."""
    cmd = ['hammock', file1_list_path, file2_list_path, 
           '--mode', mode,
           f'--{sketch_type}',
           '--precision', str(precision)]
    
    if use_rust and sketch_type == 'hyperloglog':
        cmd.append('--rust')
    
    start_time = time.process_time()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        end_time = time.process_time()
        return end_time - start_time
    except subprocess.CalledProcessError as e:
        print(f"\nError running hammock command: {' '.join(cmd)}")
        print(f"Return code: {e.returncode}")
        print(f"Error output: {e.stderr}")
        print(f"Standard output: {e.stdout}")
        raise

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
            
            # Create input files for both tools
            file1_list_path = os.path.join(temp_dir, "file1_list.txt")
            file2_list_path = os.path.join(temp_dir, "file2_list.txt")
            
            with open(file1_list_path, 'w') as f:
                f.write('\n'.join(file1_list))
            with open(file2_list_path, 'w') as f:
                f.write('\n'.join(file2_list))
            
            # Run bedtools benchmark
            bedtools_times = []
            for run in range(num_runs):
                time_taken = run_bedtools_benchmark(file1_list_path, file2_list_path)
                bedtools_times.append(time_taken)
            
            # Run hammock benchmarks with different modes and sketching options
            hammock_results = {}
            for mode in ['A', 'B', 'C']:
                for sketch_type in ['hyperloglog', 'minhash']:
                    for precision in [8, 12, 16]:
                        # Run with Python implementation
                        key = f"hammock_{mode}_{sketch_type}_p{precision}"
                        times = []
                        for run in range(num_runs):
                            try:
                                time_taken = run_hammock_benchmark(
                                    file1_list_path, file2_list_path,
                                    mode=mode, sketch_type=sketch_type,
                                    precision=precision, use_rust=False
                                )
                                times.append(time_taken)
                            except Exception as e:
                                print(f"Skipping {key} run {run+1} due to error: {e}")
                        
                        if times:
                            hammock_results[key] = {
                                'mean_time': np.mean(times),
                                'std_time': np.std(times) if len(times) > 1 else 0,
                                'min_time': np.min(times),
                                'max_time': np.max(times)
                            }
                        else:
                            print(f"All runs failed for {key}, skipping this configuration")
                            hammock_results[key] = {
                                'mean_time': float('nan'),
                                'std_time': float('nan'),
                                'min_time': float('nan'),
                                'max_time': float('nan')
                            }
                        
                        # Run with Rust implementation for HyperLogLog
                        if sketch_type == 'hyperloglog':
                            key = f"hammock_{mode}_{sketch_type}_rust_p{precision}"
                            times = []
                            for run in range(num_runs):
                                try:
                                    time_taken = run_hammock_benchmark(
                                        file1_list_path, file2_list_path,
                                        mode=mode, sketch_type=sketch_type,
                                        precision=precision, use_rust=True
                                    )
                                    times.append(time_taken)
                                except Exception as e:
                                    print(f"Skipping {key} run {run+1} due to error: {e}")
                            
                            if times:
                                hammock_results[key] = {
                                    'mean_time': np.mean(times),
                                    'std_time': np.std(times) if len(times) > 1 else 0,
                                    'min_time': np.min(times),
                                    'max_time': np.max(times)
                                }
                            else:
                                print(f"All runs failed for {key}, skipping this configuration")
                                hammock_results[key] = {
                                    'mean_time': float('nan'),
                                    'std_time': float('nan'),
                                    'min_time': float('nan'),
                                    'max_time': float('nan')
                                }
            
            # Store results
            result = {
                'num_files': num_files,
                'bedtools': {
                    'mean_time': np.mean(bedtools_times),
                    'std_time': np.std(bedtools_times),
                    'min_time': np.min(bedtools_times),
                    'max_time': np.max(bedtools_times)
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
    
    # Plot bedtools results
    bedtools_mean = [r['bedtools']['mean_time'] for r in results]
    bedtools_std = [r['bedtools']['std_time'] for r in results]
    ax1.errorbar(num_files, bedtools_mean, yerr=bedtools_std, 
                fmt='o-', capsize=5, label='bedtools')
    
    # Plot hammock results for different modes and sketching options
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    markers = ['o', 's', '^', 'D', 'v', '<']
    
    # Plot Python implementation results
    for idx, (mode, sketch_type) in enumerate([('A', 'hyperloglog'), ('B', 'hyperloglog'),
                                              ('C', 'hyperloglog'), ('A', 'minhash'),
                                              ('B', 'minhash'), ('C', 'minhash')]):
        for precision in [8, 12, 16]:
            key = f"hammock_{mode}_{sketch_type}_p{precision}"
            if all(key in r and not np.isnan(r[key]['mean_time']) for r in results):
                mean_times = [r[key]['mean_time'] for r in results]
                std_times = [r[key]['std_time'] for r in results]
                label = f"{key} (p={precision})"
                ax2.errorbar(num_files, mean_times, yerr=std_times,
                            fmt=f'{markers[idx]}-', capsize=5, color=colors[idx],
                            label=label)
            else:
                print(f"Skipping plot for {key} due to missing or NaN results")
    
    # Plot Rust implementation results
    for idx, mode in enumerate(['A', 'B', 'C']):
        for precision in [8, 12, 16]:
            key = f"hammock_{mode}_hyperloglog_rust_p{precision}"
            if all(key in r and not np.isnan(r[key]['mean_time']) for r in results):
                mean_times = [r[key]['mean_time'] for r in results]
                std_times = [r[key]['std_time'] for r in results]
                label = f"{key} (p={precision})"
                ax2.errorbar(num_files, mean_times, yerr=std_times,
                            fmt=f'{markers[idx]}--', capsize=5, color=colors[idx],
                            label=label)
            else:
                print(f"Skipping plot for {key} due to missing or NaN results")
    
    # Customize plots
    for ax in [ax1, ax2]:
        ax.set_xlabel('Number of Files')
        ax.set_ylabel('CPU Time (seconds)')
        ax.grid(True)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    ax1.set_title('bedtools.sh Performance')
    ax2.set_title('hammock Performance (solid=Python, dashed=Rust)')
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save plot
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    plot_path = os.path.join('benchmarks/results', f'bedtools_hammock_benchmark_{timestamp}.png')
    plt.savefig(plot_path, bbox_inches='tight')
    plt.close()
    
    # Save numerical results
    results_path = os.path.join('benchmarks/results', f'bedtools_hammock_benchmark_{timestamp}.txt')
    with open(results_path, 'w') as f:
        f.write("Benchmark Results:\n")
        f.write("================\n\n")
        for r in results:
            f.write(f"Number of files: {r['num_files']}\n")
            f.write("\nbedtools:\n")
            f.write(f"  Mean CPU time: {r['bedtools']['mean_time']:.3f} seconds\n")
            f.write(f"  Standard deviation: {r['bedtools']['std_time']:.3f} seconds\n")
            f.write(f"  Min CPU time: {r['bedtools']['min_time']:.3f} seconds\n")
            f.write(f"  Max CPU time: {r['bedtools']['max_time']:.3f} seconds\n")
            f.write("\nhammock:\n")
            for key in r:
                if key.startswith('hammock_'):
                    f.write(f"\n  {key}:\n")
                    f.write(f"    Mean CPU time: {r[key]['mean_time']:.3f} seconds\n")
                    f.write(f"    Standard deviation: {r[key]['std_time']:.3f} seconds\n")
                    f.write(f"    Min CPU time: {r[key]['min_time']:.3f} seconds\n")
                    f.write(f"    Max CPU time: {r[key]['max_time']:.3f} seconds\n")
            f.write("\n" + "="*50 + "\n\n")

def main():
    # Create results directory if it doesn't exist
    os.makedirs('benchmarks/results', exist_ok=True)
    
    print("Starting bedtools.sh and hammock benchmark...")
    results = run_benchmark()
    plot_results(results)
    print("\nBenchmark completed. Results saved in benchmarks/results/")

if __name__ == "__main__":
    main() 