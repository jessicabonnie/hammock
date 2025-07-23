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

def run_hammock_benchmark(file1_list_path: str, file2_list_path: str, mode='A', sketch_type='hyperloglog', precision=12, use_rust=False, python_only=False) -> float:
    """Run hammock benchmark with specified parameters."""
    cmd = ['hammock', file1_list_path, file2_list_path, 
           '--mode', mode,
           f'--{sketch_type}',
           '--precision', str(precision)]
    
    # Add --python-only flag if requested
    if python_only:
        cmd.append('--python-only')
    
    # Note: The --rust flag has been removed. HyperLogLog now automatically uses
    # Cython acceleration when available (FastHyperLogLog), falling back to
    # pure Python implementation when not available.
    
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
                        # Run with automatic acceleration (Cython when available, Python fallback)
                        key = f"hammock_{mode}_{sketch_type}_p{precision}"
                        times = []
                        for run in range(num_runs):
                            try:
                                time_taken = run_hammock_benchmark(
                                    file1_list_path, file2_list_path,
                                    mode=mode, sketch_type=sketch_type,
                                    precision=precision, use_rust=False, python_only=False
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
                        
                        # Run with Python-only implementation for HyperLogLog
                        if sketch_type == 'hyperloglog':
                            key = f"hammock_{mode}_{sketch_type}_python_p{precision}"
                            times = []
                            for run in range(num_runs):
                                try:
                                    time_taken = run_hammock_benchmark(
                                        file1_list_path, file2_list_path,
                                        mode=mode, sketch_type=sketch_type,
                                        precision=precision, use_rust=False, python_only=True
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
    """Plot the benchmark results with separate graphs for each mode and a comparison graph."""
    num_files = [r['num_files'] for r in results]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Colors and markers for different configurations
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']
    
    # Create separate graphs for each mode (A, B, C)
    for mode in ['A', 'B', 'C']:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
        # Plot bedtools for reference
        bedtools_mean = [r['bedtools']['mean_time'] for r in results]
        bedtools_std = [r['bedtools']['std_time'] for r in results]
        ax.errorbar(num_files, bedtools_mean, yerr=bedtools_std, 
                   fmt='ko-', capsize=5, label='bedtools', linewidth=2, markersize=8)
        
        # Plot hammock results for this mode
        color_idx = 0
        for sketch_type in ['hyperloglog', 'minhash']:
            for precision in [8, 12, 16]:
                # Cython version (automatic acceleration)
                key = f"hammock_{mode}_{sketch_type}_p{precision}"
                if all(key in r and not np.isnan(r[key]['mean_time']) for r in results):
                    mean_times = [r[key]['mean_time'] for r in results]
                    std_times = [r[key]['std_time'] for r in results]
                    label = f"hammock {mode} {sketch_type} p{precision} (Cython)"
                    ax.errorbar(num_files, mean_times, yerr=std_times,
                               fmt=f'{markers[color_idx]}-', capsize=5, 
                               color=colors[color_idx], label=label, linewidth=2)
                    color_idx += 1
                
                # Python-only version (only for HyperLogLog)
                if sketch_type == 'hyperloglog':
                    key = f"hammock_{mode}_{sketch_type}_python_p{precision}"
                    if all(key in r and not np.isnan(r[key]['mean_time']) for r in results):
                        mean_times = [r[key]['mean_time'] for r in results]
                        std_times = [r[key]['std_time'] for r in results]
                        label = f"hammock {mode} {sketch_type} p{precision} (Python)"
                        ax.errorbar(num_files, mean_times, yerr=std_times,
                                   fmt=f'{markers[color_idx-1]}--', capsize=5, 
                                   color=colors[color_idx-1], label=label, linewidth=2)
        
        ax.set_xlabel('Number of Files')
        ax.set_ylabel('CPU Time (seconds)')
        ax.set_title(f'hammock Mode {mode} Performance Comparison')
        ax.grid(True, alpha=0.3)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        # Save plot
        plot_path = os.path.join('benchmarks/results', f'hammock_mode_{mode}_benchmark_{timestamp}.png')
        plt.savefig(plot_path, bbox_inches='tight', dpi=300)
        plt.close()
    
    # Create comparison graph: bedtools vs hammock mode B HyperLogLog (Python vs Cython)
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    # Plot bedtools
    bedtools_mean = [r['bedtools']['mean_time'] for r in results]
    bedtools_std = [r['bedtools']['std_time'] for r in results]
    ax.errorbar(num_files, bedtools_mean, yerr=bedtools_std, 
               fmt='ko-', capsize=5, label='bedtools', linewidth=3, markersize=10)
    
    # Calculate average across precision levels for hammock mode B HyperLogLog
    cython_times = []
    python_times = []
    
    for file_idx in range(len(num_files)):
        cython_file_times = []
        python_file_times = []
        
        for precision in [8, 12, 16]:
            # Cython version
            key = f"hammock_B_hyperloglog_p{precision}"
            if key in results[file_idx] and not np.isnan(results[file_idx][key]['mean_time']):
                cython_file_times.append(results[file_idx][key]['mean_time'])
            
            # Python-only version
            key = f"hammock_B_hyperloglog_python_p{precision}"
            if key in results[file_idx] and not np.isnan(results[file_idx][key]['mean_time']):
                python_file_times.append(results[file_idx][key]['mean_time'])
        
        if cython_file_times:
            cython_times.append(np.mean(cython_file_times))
        else:
            cython_times.append(np.nan)
            
        if python_file_times:
            python_times.append(np.mean(python_file_times))
        else:
            python_times.append(np.nan)
    
    # Plot averaged hammock results
    if not all(np.isnan(cython_times)):
        ax.errorbar(num_files, cython_times, fmt='o-', capsize=5, 
                   label='hammock B HyperLogLog (Cython, avg)', linewidth=2, color='blue')
    
    if not all(np.isnan(python_times)):
        ax.errorbar(num_files, python_times, fmt='s--', capsize=5, 
                   label='hammock B HyperLogLog (Python, avg)', linewidth=2, color='red')
    
    ax.set_xlabel('Number of Files')
    ax.set_ylabel('CPU Time (seconds)')
    ax.set_title('bedtools vs hammock Mode B HyperLogLog Comparison')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    # Save comparison plot
    plot_path = os.path.join('benchmarks/results', f'bedtools_vs_hammock_B_comparison_{timestamp}.png')
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
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
    
    # Save results table
    table_path = os.path.join('benchmarks/results', f'benchmark_table_{timestamp}.csv')
    save_results_table(results, table_path)

def save_results_table(results, output_path):
    """Save benchmark results as a CSV table."""
    import csv
    
    # Define all possible columns
    columns = ['num_files', 'bedtools_mean', 'bedtools_std', 'bedtools_min', 'bedtools_max']
    
    # Add hammock columns for all configurations
    for mode in ['A', 'B', 'C']:
        for sketch_type in ['hyperloglog', 'minhash']:
            for precision in [8, 12, 16]:
                base_key = f"hammock_{mode}_{sketch_type}_p{precision}"
                columns.extend([
                    f"{base_key}_mean",
                    f"{base_key}_std", 
                    f"{base_key}_min",
                    f"{base_key}_max"
                ])
                
                # Add Python-only version for HyperLogLog
                if sketch_type == 'hyperloglog':
                    python_key = f"hammock_{mode}_{sketch_type}_python_p{precision}"
                    columns.extend([
                        f"{python_key}_mean",
                        f"{python_key}_std",
                        f"{python_key}_min", 
                        f"{python_key}_max"
                    ])
    
    # Write CSV file
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=columns)
        writer.writeheader()
        
        for r in results:
            row = {
                'num_files': r['num_files'],
                'bedtools_mean': f"{r['bedtools']['mean_time']:.6f}",
                'bedtools_std': f"{r['bedtools']['std_time']:.6f}",
                'bedtools_min': f"{r['bedtools']['min_time']:.6f}",
                'bedtools_max': f"{r['bedtools']['max_time']:.6f}"
            }
            
            # Add hammock results
            for mode in ['A', 'B', 'C']:
                for sketch_type in ['hyperloglog', 'minhash']:
                    for precision in [8, 12, 16]:
                        base_key = f"hammock_{mode}_{sketch_type}_p{precision}"
                        if base_key in r and not np.isnan(r[base_key]['mean_time']):
                            row[f"{base_key}_mean"] = f"{r[base_key]['mean_time']:.6f}"
                            row[f"{base_key}_std"] = f"{r[base_key]['std_time']:.6f}"
                            row[f"{base_key}_min"] = f"{r[base_key]['min_time']:.6f}"
                            row[f"{base_key}_max"] = f"{r[base_key]['max_time']:.6f}"
                        else:
                            row[f"{base_key}_mean"] = "NaN"
                            row[f"{base_key}_std"] = "NaN"
                            row[f"{base_key}_min"] = "NaN"
                            row[f"{base_key}_max"] = "NaN"
                        
                        # Add Python-only version for HyperLogLog
                        if sketch_type == 'hyperloglog':
                            python_key = f"hammock_{mode}_{sketch_type}_python_p{precision}"
                            if python_key in r and not np.isnan(r[python_key]['mean_time']):
                                row[f"{python_key}_mean"] = f"{r[python_key]['mean_time']:.6f}"
                                row[f"{python_key}_std"] = f"{r[python_key]['std_time']:.6f}"
                                row[f"{python_key}_min"] = f"{r[python_key]['min_time']:.6f}"
                                row[f"{python_key}_max"] = f"{r[python_key]['max_time']:.6f}"
                            else:
                                row[f"{python_key}_mean"] = "NaN"
                                row[f"{python_key}_std"] = "NaN"
                                row[f"{python_key}_min"] = "NaN"
                                row[f"{python_key}_max"] = "NaN"
            
            writer.writerow(row)
    
    print(f"Results table saved to: {output_path}")

def quick_test():
    """Run a quick test to verify graphing functionality."""
    print("Running quick test to verify graphing functionality...")
    
    # Create results directory if it doesn't exist
    os.makedirs('benchmarks/results', exist_ok=True)
    
    # Run with minimal parameters for quick test
    results = run_benchmark(num_files_list=[2, 4], num_runs=2)
    plot_results(results)
    print("\nQuick test completed. Check benchmarks/results/ for generated graphs.")

def main():
    # Create results directory if it doesn't exist
    os.makedirs('benchmarks/results', exist_ok=True)
    
    print("Starting bedtools.sh and hammock benchmark...")
    results = run_benchmark()
    plot_results(results)
    print("\nBenchmark completed. Results saved in benchmarks/results/")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--test":
        quick_test()
    else:
        main() 