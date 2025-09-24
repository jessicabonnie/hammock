#!/usr/bin/env python3
"""
Benchmark script comparing hammock performance against bedtools.

This version includes a dedicated graph comparing comparison times between
C++ implementation, bedtools, and Python implementation.

Note: MinHash benchmarking is currently disabled but the code is preserved
for future use. To re-enable MinHash benchmarking, remove the conditional
check in the main benchmarking loop.
"""

import os
import time
import random
import subprocess
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from datetime import datetime
import tempfile
from typing import List, Dict, Any
import psutil
import platform
import threading
import queue

# Get the absolute path to the bedtools.sh script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BEDTOOLS_SCRIPT = os.path.join(SCRIPT_DIR, 'bedtools.sh')

def get_memory_usage():
    """Get current memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024

def monitor_memory_usage(stop_event, memory_queue):
    """Monitor memory usage in a separate thread."""
    process = psutil.Process(os.getpid())
    peak_memory = 0
    
    while not stop_event.is_set():
        try:
            current_memory = process.memory_info().rss / 1024 / 1024
            peak_memory = max(peak_memory, current_memory)
            time.sleep(0.01)  # Sample every 10ms
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            break
    
    memory_queue.put(peak_memory)

def get_file_sizes(file_list_path: str) -> Dict[str, int]:
    """Get file sizes for all files in the list."""
    sizes = {}
    with open(file_list_path, 'r') as f:
        for file_path in f:
            file_path = file_path.strip()
            if os.path.exists(file_path):
                sizes[file_path] = os.path.getsize(file_path)
    return sizes

def get_system_info() -> Dict[str, Any]:
    """Get system information for reproducibility."""
    return {
        'cpu_count': psutil.cpu_count(),
        'memory_total_gb': psutil.virtual_memory().total / (1024**3),
        'platform': platform.platform(),
        'python_version': platform.python_version()
    }

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

def run_bedtools_benchmark(file1_list_path: str, file2_list_path: str) -> Dict[str, Any]:
    """Run bedtools benchmark with the given file lists."""
    # Set up memory monitoring
    stop_event = threading.Event()
    memory_queue = queue.Queue()
    memory_thread = threading.Thread(target=monitor_memory_usage, args=(stop_event, memory_queue))
    memory_thread.start()
    
    # Capture initial memory
    memory_start = get_memory_usage()
    
    # Capture timing
    wall_start = time.time()
    cpu_start = time.process_time()
    
    # Run the command
    result = subprocess.run(['bash', BEDTOOLS_SCRIPT, file1_list_path, file2_list_path], 
                          check=True, capture_output=True, text=True)
    
    # Capture final timing
    cpu_end = time.process_time()
    wall_end = time.time()
    memory_end = get_memory_usage()
    
    # Stop memory monitoring and get peak memory
    stop_event.set()
    memory_thread.join()
    memory_peak = memory_queue.get()
    
    # Get file sizes
    file1_sizes = get_file_sizes(file1_list_path)
    file2_sizes = get_file_sizes(file2_list_path)
    
    return {
        'cpu_time': cpu_end - cpu_start,
        'wall_time': wall_end - wall_start,
        'memory_start_mb': memory_start,
        'memory_end_mb': memory_end,
        'memory_peak_mb': memory_peak,
        'memory_delta_mb': memory_end - memory_start,
        'file1_sizes': file1_sizes,
        'file2_sizes': file2_sizes,
        'stdout': result.stdout,
        'stderr': result.stderr,
        'return_code': result.returncode
    }

def run_hammock_benchmark(file1_list_path: str, file2_list_path: str, mode='A', sketch_type='hyperloglog', precision=12, python_only=False) -> Dict[str, Any]:
    """Run hammock benchmark with specified parameters."""
    cmd = ['hammock', file1_list_path, file2_list_path, 
           '--mode', mode,
           f'--{sketch_type}',
           '--precision', str(precision)]
    
    # Add --python-only flag if requested
    if python_only:
        cmd.append('--python-only')
    
    # Note: HyperLogLog now automatically uses C++ acceleration when available 
    # (FastHyperLogLog), falling back to Cython or pure Python implementation 
    # when not available.
    
    # Set up memory monitoring
    stop_event = threading.Event()
    memory_queue = queue.Queue()
    memory_thread = threading.Thread(target=monitor_memory_usage, args=(stop_event, memory_queue))
    memory_thread.start()
    
    # Capture initial memory
    memory_start = get_memory_usage()
    
    # Capture timing
    wall_start = time.time()
    cpu_start = time.process_time()
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Capture final timing
        cpu_end = time.process_time()
        wall_end = time.time()
        memory_end = get_memory_usage()
        
        # Stop memory monitoring and get peak memory
        stop_event.set()
        memory_thread.join()
        memory_peak = memory_queue.get()
        
        # Get file sizes
        file1_sizes = get_file_sizes(file1_list_path)
        file2_sizes = get_file_sizes(file2_list_path)
        
        return {
            'cpu_time': cpu_end - cpu_start,
            'wall_time': wall_end - wall_start,
            'memory_start_mb': memory_start,
            'memory_end_mb': memory_end,
            'memory_peak_mb': memory_peak,
            'memory_delta_mb': memory_end - memory_start,
            'file1_sizes': file1_sizes,
            'file2_sizes': file2_sizes,
            'stdout': result.stdout,
            'stderr': result.stderr,
            'return_code': result.returncode
        }
    except subprocess.CalledProcessError as e:
        # Stop memory monitoring
        stop_event.set()
        memory_thread.join()
        
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
            bedtools_runs = []
            for run in range(num_runs):
                run_data = run_bedtools_benchmark(file1_list_path, file2_list_path)
                bedtools_runs.append(run_data)
            
            # Extract timing and memory data for bedtools
            bedtools_cpu_times = [run['cpu_time'] for run in bedtools_runs]
            bedtools_wall_times = [run['wall_time'] for run in bedtools_runs]
            bedtools_memory_peaks = [run['memory_peak_mb'] for run in bedtools_runs]
            bedtools_memory_deltas = [run['memory_delta_mb'] for run in bedtools_runs]
            
            # Run hammock benchmarks with different modes and sketching options
            hammock_results = {}
            for mode in ['A', 'B', 'C']:
                for sketch_type in ['hyperloglog', 'minhash']:
                    # Skip MinHash benchmarking for now (keep code for future use)
                    if sketch_type == 'minhash':
                        continue
                        
                    for precision in [8, 12, 16]:
                        # Run with automatic acceleration (Cython when available, Python fallback)
                        key = f"hammock_{mode}_{sketch_type}_p{precision}"
                        runs = []
                        for run in range(num_runs):
                            try:
                                run_data = run_hammock_benchmark(
                                    file1_list_path, file2_list_path,
                                    mode=mode, sketch_type=sketch_type,
                                    precision=precision, python_only=False
                                )
                                runs.append(run_data)
                            except Exception as e:
                                print(f"Skipping {key} run {run+1} due to error: {e}")
                        
                        if runs:
                            # Extract timing and memory data
                            cpu_times = [run['cpu_time'] for run in runs]
                            wall_times = [run['wall_time'] for run in runs]
                            memory_peaks = [run['memory_peak_mb'] for run in runs]
                            memory_deltas = [run['memory_delta_mb'] for run in runs]
                            
                            hammock_results[key] = {
                                'mean_cpu_time': np.mean(cpu_times),
                                'std_cpu_time': np.std(cpu_times) if len(cpu_times) > 1 else 0,
                                'min_cpu_time': np.min(cpu_times),
                                'max_cpu_time': np.max(cpu_times),
                                'mean_wall_time': np.mean(wall_times),
                                'std_wall_time': np.std(wall_times) if len(wall_times) > 1 else 0,
                                'min_wall_time': np.min(wall_times),
                                'max_wall_time': np.max(wall_times),
                                'mean_memory_peak': np.mean(memory_peaks),
                                'std_memory_peak': np.std(memory_peaks) if len(memory_peaks) > 1 else 0,
                                'min_memory_peak': np.min(memory_peaks),
                                'max_memory_peak': np.max(memory_peaks),
                                'mean_memory_delta': np.mean(memory_deltas),
                                'std_memory_delta': np.std(memory_deltas) if len(memory_deltas) > 1 else 0,
                                'min_memory_delta': np.min(memory_deltas),
                                'max_memory_delta': np.max(memory_deltas),
                                'detailed_runs': runs
                            }
                        else:
                            print(f"All runs failed for {key}, skipping this configuration")
                            hammock_results[key] = {
                                'mean_cpu_time': float('nan'),
                                'std_cpu_time': float('nan'),
                                'min_cpu_time': float('nan'),
                                'max_cpu_time': float('nan'),
                                'mean_wall_time': float('nan'),
                                'std_wall_time': float('nan'),
                                'min_wall_time': float('nan'),
                                'max_wall_time': float('nan'),
                                'mean_memory_peak': float('nan'),
                                'std_memory_peak': float('nan'),
                                'min_memory_peak': float('nan'),
                                'max_memory_peak': float('nan'),
                                'mean_memory_delta': float('nan'),
                                'std_memory_delta': float('nan'),
                                'min_memory_delta': float('nan'),
                                'max_memory_delta': float('nan'),
                                'detailed_runs': []
                            }
                        
                        # Run with Python-only implementation for HyperLogLog
                        if sketch_type == 'hyperloglog':
                            key = f"hammock_{mode}_{sketch_type}_python_p{precision}"
                            runs = []
                            for run in range(num_runs):
                                try:
                                    run_data = run_hammock_benchmark(
                                        file1_list_path, file2_list_path,
                                        mode=mode, sketch_type=sketch_type,
                                        precision=precision, python_only=True
                                    )
                                    runs.append(run_data)
                                except Exception as e:
                                    print(f"Skipping {key} run {run+1} due to error: {e}")
                            
                            if runs:
                                # Extract timing and memory data
                                cpu_times = [run['cpu_time'] for run in runs]
                                wall_times = [run['wall_time'] for run in runs]
                                memory_peaks = [run['memory_peak_mb'] for run in runs]
                                memory_deltas = [run['memory_delta_mb'] for run in runs]
                                
                                hammock_results[key] = {
                                    'mean_cpu_time': np.mean(cpu_times),
                                    'std_cpu_time': np.std(cpu_times) if len(cpu_times) > 1 else 0,
                                    'min_cpu_time': np.min(cpu_times),
                                    'max_cpu_time': np.max(cpu_times),
                                    'mean_wall_time': np.mean(wall_times),
                                    'std_wall_time': np.std(wall_times) if len(wall_times) > 1 else 0,
                                    'min_wall_time': np.min(wall_times),
                                    'max_wall_time': np.max(wall_times),
                                    'mean_memory_peak': np.mean(memory_peaks),
                                    'std_memory_peak': np.std(memory_peaks) if len(memory_peaks) > 1 else 0,
                                    'min_memory_peak': np.min(memory_peaks),
                                    'max_memory_peak': np.max(memory_peaks),
                                    'mean_memory_delta': np.mean(memory_deltas),
                                    'std_memory_delta': np.std(memory_deltas) if len(memory_deltas) > 1 else 0,
                                    'min_memory_delta': np.min(memory_deltas),
                                    'max_memory_delta': np.max(memory_deltas),
                                    'detailed_runs': runs
                                }
                            else:
                                print(f"All runs failed for {key}, skipping this configuration")
                                hammock_results[key] = {
                                    'mean_cpu_time': float('nan'),
                                    'std_cpu_time': float('nan'),
                                    'min_cpu_time': float('nan'),
                                    'max_cpu_time': float('nan'),
                                    'mean_wall_time': float('nan'),
                                    'std_wall_time': float('nan'),
                                    'min_wall_time': float('nan'),
                                    'max_wall_time': float('nan'),
                                    'mean_memory_peak': float('nan'),
                                    'std_memory_peak': float('nan'),
                                    'min_memory_peak': float('nan'),
                                    'max_memory_peak': float('nan'),
                                    'mean_memory_delta': float('nan'),
                                    'std_memory_delta': float('nan'),
                                    'min_memory_delta': float('nan'),
                                    'max_memory_delta': float('nan'),
                                    'detailed_runs': []
                                }
            
            # Store results with comprehensive bedtools data
            result = {
                'num_files': num_files,
                'bedtools': {
                    'mean_cpu_time': np.mean(bedtools_cpu_times),
                    'std_cpu_time': np.std(bedtools_cpu_times),
                    'min_cpu_time': np.min(bedtools_cpu_times),
                    'max_cpu_time': np.max(bedtools_cpu_times),
                    'mean_wall_time': np.mean(bedtools_wall_times),
                    'std_wall_time': np.std(bedtools_wall_times),
                    'min_wall_time': np.min(bedtools_wall_times),
                    'max_wall_time': np.max(bedtools_wall_times),
                    'mean_memory_peak': np.mean(bedtools_memory_peaks),
                    'std_memory_peak': np.std(bedtools_memory_peaks),
                    'min_memory_peak': np.min(bedtools_memory_peaks),
                    'max_memory_peak': np.max(bedtools_memory_peaks),
                    'mean_memory_delta': np.mean(bedtools_memory_deltas),
                    'std_memory_delta': np.std(bedtools_memory_deltas),
                    'min_memory_delta': np.min(bedtools_memory_deltas),
                    'max_memory_delta': np.max(bedtools_memory_deltas),
                    'detailed_runs': bedtools_runs
                }
            }
            result.update(hammock_results)
            results.append(result)
    
    return results

def plot_results(results, timestamp=None):
    """Plot the benchmark results with separate graphs for each mode and a comparison graph."""
    num_files = [r['num_files'] for r in results]
    if timestamp is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Colors and markers for different configurations
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']
    
    # Create separate graphs for each mode (A, B, C)
    for mode in ['A', 'B', 'C']:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
        # Plot bedtools for reference
        bedtools_mean = [r['bedtools']['mean_cpu_time'] for r in results]
        bedtools_std = [r['bedtools']['std_cpu_time'] for r in results]
        ax.errorbar(num_files, bedtools_mean, yerr=bedtools_std, 
                   fmt='ko-', capsize=5, label='bedtools', linewidth=2, markersize=8)
        
        # Plot hammock results for this mode
        color_idx = 0
        for sketch_type in ['hyperloglog', 'minhash']:
            # Note: MinHash benchmarking is currently disabled, so this loop will only
            # process HyperLogLog data. The MinHash code is preserved for future use.
            for precision in [8, 12, 16]:
                # Cython version (automatic acceleration)
                key = f"hammock_{mode}_{sketch_type}_p{precision}"
                if all(key in r and not np.isnan(r[key]['mean_cpu_time']) for r in results):
                    mean_times = [r[key]['mean_cpu_time'] for r in results]
                    std_times = [r[key]['std_cpu_time'] for r in results]
                    label = f"hammock {mode} {sketch_type} p{precision} (Cython)"
                    ax.errorbar(num_files, mean_times, yerr=std_times,
                               fmt=f'{markers[color_idx]}:', capsize=5, 
                               color=colors[color_idx], label=label, linewidth=2)
                    color_idx += 1
                
                # Python-only version (only for HyperLogLog)
                if sketch_type == 'hyperloglog':
                    key = f"hammock_{mode}_{sketch_type}_python_p{precision}"
                    if all(key in r and not np.isnan(r[key]['mean_cpu_time']) for r in results):
                        mean_times = [r[key]['mean_cpu_time'] for r in results]
                        std_times = [r[key]['std_cpu_time'] for r in results]
                        label = f"hammock {mode} {sketch_type} p{precision} (Python)"
                        ax.errorbar(num_files, mean_times, yerr=std_times,
                                   fmt=f'{markers[color_idx-1]}-', capsize=5, 
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
    bedtools_mean = [r['bedtools']['mean_cpu_time'] for r in results]
    bedtools_std = [r['bedtools']['std_cpu_time'] for r in results]
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
            if key in results[file_idx] and not np.isnan(results[file_idx][key]['mean_cpu_time']):
                cython_file_times.append(results[file_idx][key]['mean_cpu_time'])
            
            # Python-only version
            key = f"hammock_B_hyperloglog_python_p{precision}"
            if key in results[file_idx] and not np.isnan(results[file_idx][key]['mean_cpu_time']):
                python_file_times.append(results[file_idx][key]['mean_cpu_time'])
        
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
        ax.errorbar(num_files, cython_times, fmt='o:', capsize=5, 
                   label='hammock B HyperLogLog (Cython, avg)', linewidth=2, color='blue')
    
    if not all(np.isnan(python_times)):
        ax.errorbar(num_files, python_times, fmt='s-', capsize=5, 
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
    
    # Create separate graphs for mode B HyperLogLog and MinHash
    create_mode_B_separate_graphs(results, timestamp)
    
    # Create comparison-focused graph: bedtools vs hammock comparison times
    create_comparison_focused_graph(results, timestamp)
    
    # Save numerical results
    results_path = os.path.join('benchmarks/results', f'bedtools_hammock_benchmark_{timestamp}.txt')
    with open(results_path, 'w') as f:
        f.write("Benchmark Results:\n")
        f.write("================\n\n")
        for r in results:
            f.write(f"Number of files: {r['num_files']}\n")
            f.write("\nbedtools:\n")
            f.write(f"  CPU Time - Mean: {r['bedtools']['mean_cpu_time']:.3f} ± {r['bedtools']['std_cpu_time']:.3f} seconds (min: {r['bedtools']['min_cpu_time']:.3f}, max: {r['bedtools']['max_cpu_time']:.3f})\n")
            f.write(f"  Wall Time - Mean: {r['bedtools']['mean_wall_time']:.3f} ± {r['bedtools']['std_wall_time']:.3f} seconds (min: {r['bedtools']['min_wall_time']:.3f}, max: {r['bedtools']['max_wall_time']:.3f})\n")
            f.write(f"  Memory Peak - Mean: {r['bedtools']['mean_memory_peak']:.1f} ± {r['bedtools']['std_memory_peak']:.1f} MB (min: {r['bedtools']['min_memory_peak']:.1f}, max: {r['bedtools']['max_memory_peak']:.1f})\n")
            f.write(f"  Memory Delta - Mean: {r['bedtools']['mean_memory_delta']:.1f} ± {r['bedtools']['std_memory_delta']:.1f} MB (min: {r['bedtools']['min_memory_delta']:.1f}, max: {r['bedtools']['max_memory_delta']:.1f})\n")
            f.write("\nhammock:\n")
            for key in r:
                if key.startswith('hammock_'):
                    f.write(f"\n  {key}:\n")
                    f.write(f"    CPU Time - Mean: {r[key]['mean_cpu_time']:.3f} ± {r[key]['std_cpu_time']:.3f} seconds (min: {r[key]['min_cpu_time']:.3f}, max: {r[key]['max_cpu_time']:.3f})\n")
                    f.write(f"    Wall Time - Mean: {r[key]['mean_wall_time']:.3f} ± {r[key]['std_wall_time']:.3f} seconds (min: {r[key]['min_wall_time']:.3f}, max: {r[key]['max_wall_time']:.3f})\n")
                    f.write(f"    Memory Peak - Mean: {r[key]['mean_memory_peak']:.1f} ± {r[key]['std_memory_peak']:.1f} MB (min: {r[key]['min_memory_peak']:.1f}, max: {r[key]['max_memory_peak']:.1f})\n")
                    f.write(f"    Memory Delta - Mean: {r[key]['mean_memory_delta']:.1f} ± {r[key]['std_memory_delta']:.1f} MB (min: {r[key]['min_memory_delta']:.1f}, max: {r[key]['max_memory_delta']:.1f})\n")
            f.write("\n" + "="*50 + "\n\n")
    
    # Save results table
    table_path = os.path.join('benchmarks/results', f'benchmark_table_{timestamp}.csv')
    save_results_table(results, table_path)

def save_results_table(results, output_path):
    """Save benchmark results as a CSV table."""
    import csv
    
    # Create a comprehensive table with all metrics
    columns = ['num_files', 'tool', 'mode', 'sketch_type', 'precision', 'implementation',
               'cpu_mean', 'cpu_std', 'cpu_min', 'cpu_max',
               'wall_mean', 'wall_std', 'wall_min', 'wall_max',
               'memory_peak_mean', 'memory_peak_std', 'memory_peak_min', 'memory_peak_max',
               'memory_delta_mean', 'memory_delta_std', 'memory_delta_min', 'memory_delta_max']
    
    # Write CSV file
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=columns)
        writer.writeheader()
        
        for r in results:
            # Add bedtools results
            bedtools_runs = r['bedtools']['detailed_runs']
            bedtools_cpu_times = [run['cpu_time'] for run in bedtools_runs]
            bedtools_wall_times = [run['wall_time'] for run in bedtools_runs]
            bedtools_memory_peaks = [run['memory_peak_mb'] for run in bedtools_runs]
            bedtools_memory_deltas = [run['memory_delta_mb'] for run in bedtools_runs]
            
            row = {
                'num_files': r['num_files'],
                'tool': 'bedtools',
                'mode': 'N/A',
                'sketch_type': 'N/A',
                'precision': 'N/A',
                'implementation': 'N/A',
                'cpu_mean': f"{np.mean(bedtools_cpu_times):.6f}",
                'cpu_std': f"{np.std(bedtools_cpu_times):.6f}",
                'cpu_min': f"{np.min(bedtools_cpu_times):.6f}",
                'cpu_max': f"{np.max(bedtools_cpu_times):.6f}",
                'wall_mean': f"{np.mean(bedtools_wall_times):.6f}",
                'wall_std': f"{np.std(bedtools_wall_times):.6f}",
                'wall_min': f"{np.min(bedtools_wall_times):.6f}",
                'wall_max': f"{np.max(bedtools_wall_times):.6f}",
                'memory_peak_mean': f"{np.mean(bedtools_memory_peaks):.6f}",
                'memory_peak_std': f"{np.std(bedtools_memory_peaks):.6f}",
                'memory_peak_min': f"{np.min(bedtools_memory_peaks):.6f}",
                'memory_peak_max': f"{np.max(bedtools_memory_peaks):.6f}",
                'memory_delta_mean': f"{np.mean(bedtools_memory_deltas):.6f}",
                'memory_delta_std': f"{np.std(bedtools_memory_deltas):.6f}",
                'memory_delta_min': f"{np.min(bedtools_memory_deltas):.6f}",
                'memory_delta_max': f"{np.max(bedtools_memory_deltas):.6f}"
            }
            writer.writerow(row)
            
            # Add hammock results
            for mode in ['A', 'B', 'C']:
                for sketch_type in ['hyperloglog', 'minhash']:
                    # Note: MinHash benchmarking is currently disabled, so this loop will only
                    # process HyperLogLog data. The MinHash code is preserved for future use.
                    for precision in [8, 12, 16]:
                        base_key = f"hammock_{mode}_{sketch_type}_p{precision}"
                        if base_key in r and not np.isnan(r[base_key]['mean_cpu_time']):
                            detailed_runs = r[base_key]['detailed_runs']
                            cpu_times = [run['cpu_time'] for run in detailed_runs]
                            wall_times = [run['wall_time'] for run in detailed_runs]
                            memory_peaks = [run['memory_peak_mb'] for run in detailed_runs]
                            memory_deltas = [run['memory_delta_mb'] for run in detailed_runs]
                            
                            row = {
                                'num_files': r['num_files'],
                                'tool': 'hammock',
                                'mode': mode,
                                'sketch_type': sketch_type,
                                'precision': precision,
                                'implementation': 'Cython',
                                'cpu_mean': f"{np.mean(cpu_times):.6f}",
                                'cpu_std': f"{np.std(cpu_times):.6f}",
                                'cpu_min': f"{np.min(cpu_times):.6f}",
                                'cpu_max': f"{np.max(cpu_times):.6f}",
                                'wall_mean': f"{np.mean(wall_times):.6f}",
                                'wall_std': f"{np.std(wall_times):.6f}",
                                'wall_min': f"{np.min(wall_times):.6f}",
                                'wall_max': f"{np.max(wall_times):.6f}",
                                'memory_peak_mean': f"{np.mean(memory_peaks):.6f}",
                                'memory_peak_std': f"{np.std(memory_peaks):.6f}",
                                'memory_peak_min': f"{np.min(memory_peaks):.6f}",
                                'memory_peak_max': f"{np.max(memory_peaks):.6f}",
                                'memory_delta_mean': f"{np.mean(memory_deltas):.6f}",
                                'memory_delta_std': f"{np.std(memory_deltas):.6f}",
                                'memory_delta_min': f"{np.min(memory_deltas):.6f}",
                                'memory_delta_max': f"{np.max(memory_deltas):.6f}"
                            }
                            writer.writerow(row)
                        
                        # Add Python-only version for HyperLogLog
                        if sketch_type == 'hyperloglog':
                            python_key = f"hammock_{mode}_{sketch_type}_python_p{precision}"
                            if python_key in r and not np.isnan(r[python_key]['mean_cpu_time']):
                                detailed_runs = r[python_key]['detailed_runs']
                                cpu_times = [run['cpu_time'] for run in detailed_runs]
                                wall_times = [run['wall_time'] for run in detailed_runs]
                                memory_peaks = [run['memory_peak_mb'] for run in detailed_runs]
                                memory_deltas = [run['memory_delta_mb'] for run in detailed_runs]
                                
                                row = {
                                    'num_files': r['num_files'],
                                    'tool': 'hammock',
                                    'mode': mode,
                                    'sketch_type': sketch_type,
                                    'precision': precision,
                                    'implementation': 'Python',
                                    'cpu_mean': f"{np.mean(cpu_times):.6f}",
                                    'cpu_std': f"{np.std(cpu_times):.6f}",
                                    'cpu_min': f"{np.min(cpu_times):.6f}",
                                    'cpu_max': f"{np.max(cpu_times):.6f}",
                                    'wall_mean': f"{np.mean(wall_times):.6f}",
                                    'wall_std': f"{np.std(wall_times):.6f}",
                                    'wall_min': f"{np.min(wall_times):.6f}",
                                    'wall_max': f"{np.max(wall_times):.6f}",
                                    'memory_peak_mean': f"{np.mean(memory_peaks):.6f}",
                                    'memory_peak_std': f"{np.std(memory_peaks):.6f}",
                                    'memory_peak_min': f"{np.min(memory_peaks):.6f}",
                                    'memory_peak_max': f"{np.max(memory_peaks):.6f}",
                                    'memory_delta_mean': f"{np.mean(memory_deltas):.6f}",
                                    'memory_delta_std': f"{np.std(memory_deltas):.6f}",
                                    'memory_delta_min': f"{np.min(memory_deltas):.6f}",
                                    'memory_delta_max': f"{np.max(memory_deltas):.6f}"
                                }
                                writer.writerow(row)
    
    print(f"Results table saved to: {output_path}")
    
    # Save detailed run data as JSON
    detailed_path = output_path.replace('.csv', '_detailed.json')
    import json
    with open(detailed_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Detailed run data saved to: {detailed_path}")

def create_mode_B_separate_graphs(results, timestamp):
    """Create separate graphs for mode B HyperLogLog and MinHash implementations."""
    num_files = [r['num_files'] for r in results]
    
    # Create separate graph for mode B HyperLogLog
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    # Plot bedtools for reference
    bedtools_mean = [r['bedtools']['mean_cpu_time'] for r in results]
    bedtools_std = [r['bedtools']['std_cpu_time'] for r in results]
    ax.errorbar(num_files, bedtools_mean, yerr=bedtools_std, 
               fmt='ko-', capsize=5, label='bedtools', linewidth=3, markersize=10)
    
    # Plot mode B HyperLogLog with different precisions
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    markers = ['o', 's', '^']
    
    for idx, precision in enumerate([8, 12, 16]):
        # Cython version
        key = f"hammock_B_hyperloglog_p{precision}"
        if all(key in r and not np.isnan(r[key]['mean_cpu_time']) for r in results):
            mean_times = [r[key]['mean_cpu_time'] for r in results]
            std_times = [r[key]['std_cpu_time'] for r in results]
            label = f"hammock B HyperLogLog p{precision} (Cython)"
            ax.errorbar(num_files, mean_times, yerr=std_times,
                       fmt=f'{markers[idx]}:', capsize=5, color=colors[idx],
                       label=label, linewidth=2)
        
        # Python-only version
        key = f"hammock_B_hyperloglog_python_p{precision}"
        if all(key in r and not np.isnan(r[key]['mean_cpu_time']) for r in results):
            mean_times = [r[key]['mean_cpu_time'] for r in results]
            std_times = [r[key]['std_cpu_time'] for r in results]
            label = f"hammock B HyperLogLog p{precision} (Python)"
            ax.errorbar(num_files, mean_times, yerr=std_times,
                       fmt=f'{markers[idx]}-', capsize=5, color=colors[idx],
                       label=label, linewidth=2)
    
    ax.set_xlabel('Number of Files')
    ax.set_ylabel('CPU Time (seconds)')
    ax.set_title('Mode B HyperLogLog Performance Comparison')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    # Save HyperLogLog plot
    plot_path = os.path.join('benchmarks/results', f'hammock_B_hyperloglog_comparison_{timestamp}.png')
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    
    # Create separate graph for mode B MinHash
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    # Plot bedtools for reference
    ax.errorbar(num_files, bedtools_mean, yerr=bedtools_std, 
               fmt='ko-', capsize=5, label='bedtools', linewidth=3, markersize=10)
    
    # Plot mode B MinHash with different precisions
    for idx, precision in enumerate([8, 12, 16]):
        key = f"hammock_B_minhash_p{precision}"
        if all(key in r and not np.isnan(r[key]['mean_cpu_time']) for r in results):
            mean_times = [r[key]['mean_cpu_time'] for r in results]
            std_times = [r[key]['std_cpu_time'] for r in results]
            label = f"hammock B MinHash p{precision} (Cython)"
            ax.errorbar(num_files, mean_times, yerr=std_times,
                       fmt=f'{markers[idx]}:', capsize=5, color=colors[idx],
                       label=label, linewidth=2)
    
    ax.set_xlabel('Number of Files')
    ax.set_ylabel('CPU Time (seconds)')
    ax.set_title('Mode B MinHash Performance Comparison')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    # Save MinHash plot
    plot_path = os.path.join('benchmarks/results', f'hammock_B_minhash_comparison_{timestamp}.png')
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    plt.close()

def create_comparison_focused_graph(results, timestamp):
    """Create a graph focused on comparison times between different implementations."""
    num_files = [r['num_files'] for r in results]
    
    # Create comparison-focused graph
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # Plot bedtools for reference (this represents the "comparison" time for bedtools)
    bedtools_mean = [r['bedtools']['mean_cpu_time'] for r in results]
    bedtools_std = [r['bedtools']['std_cpu_time'] for r in results]
    ax.errorbar(num_files, bedtools_mean, yerr=bedtools_std, 
               fmt='ko-', capsize=5, label='bedtools (total time)', linewidth=3, markersize=10)
    
    # Colors and markers for different configurations
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']
    
    # Plot hammock total times for different configurations (since we don't have separated timing in the original)
    color_idx = 0
    
    # Focus on Mode B HyperLogLog for comparison (precision=16 only)
    precision = 16
    # C++/Cython version total times
    key = f"hammock_B_hyperloglog_p{precision}"
    if all(key in r and not np.isnan(r[key]['mean_cpu_time']) for r in results):
        total_times = [r[key]['mean_cpu_time'] for r in results]
        total_stds = [r[key]['std_cpu_time'] for r in results]
        label = f"hammock B HyperLogLog p{precision} (C++/Cython, total time)"
        ax.errorbar(num_files, total_times, yerr=total_stds,
                   fmt=f'{markers[color_idx]}:', capsize=5, 
                   color=colors[color_idx], label=label, linewidth=2, markersize=8)
        color_idx += 1
    
    # Python-only version total times
    key_python = f"hammock_B_hyperloglog_python_p{precision}"
    if all(key_python in r and not np.isnan(r[key_python]['mean_cpu_time']) for r in results):
        total_times = [r[key_python]['mean_cpu_time'] for r in results]
        total_stds = [r[key_python]['std_cpu_time'] for r in results]
        label = f"hammock B HyperLogLog p{precision} (Python, total time)"
        ax.errorbar(num_files, total_times, yerr=total_stds,
                   fmt=f'{markers[color_idx-1]}-', capsize=5, 
                   color=colors[color_idx-1], label=label, linewidth=2, markersize=8)
    
    ax.set_xlabel('Number of Files', fontsize=12)
    ax.set_ylabel('CPU Time (seconds)', fontsize=12)
    ax.set_title('Performance Comparison: bedtools vs hammock (C++ vs Python Implementation)', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    ax.set_yscale('log')  # Use log scale to better show the differences
    
    plt.tight_layout()
    
    # Save comparison plot
    plot_path = os.path.join('benchmarks/results', f'implementation_comparison_{timestamp}.png')
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"Implementation comparison graph saved: implementation_comparison_{timestamp}.png")

def quick_test():
    """Run a quick test to verify graphing functionality."""
    print("Running quick test to verify graphing functionality...")
    
    # Create results directory if it doesn't exist
    results_dir = 'benchmarks/results'
    os.makedirs(results_dir, exist_ok=True)
    print(f"Output will be saved to: {os.path.abspath(results_dir)}")
    
    # Generate timestamp for this run
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Run with minimal parameters for quick test
    results = run_benchmark(num_files_list=[2, 4], num_runs=2)
    plot_results(results, timestamp)
    print("\nQuick test completed. Check benchmarks/results/ for generated graphs.")

def main():
    # Create results directory if it doesn't exist
    results_dir = 'benchmarks/results'
    os.makedirs(results_dir, exist_ok=True)
    print(f"Output will be saved to: {os.path.abspath(results_dir)}")
    
    # Generate timestamp for this run
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    print("Starting bedtools.sh and hammock benchmark...")
    results = run_benchmark()
    plot_results(results, timestamp)
    print("\nBenchmark completed. Results saved in benchmarks/results/")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--test":
        quick_test()
    else:
        main() 