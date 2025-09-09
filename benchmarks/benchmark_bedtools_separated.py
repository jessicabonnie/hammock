#!/usr/bin/env python3
"""
Benchmark script comparing hammock performance against bedtools with separated timing.

This version separates sketch creation timing from comparison timing to provide
more accurate performance analysis. It also includes a dedicated graph comparing
comparison times between C++ implementation, bedtools, and Python implementation.

"""

import os
import sys
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

# Add the parent directory to the path so we can import hammock modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from benchmarks.sketch_comparison_timer import SketchComparisonTimer
from hammock.lib.intervals import IntervalSketch

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
        'return_code': result.returncode,
        'system_info': get_system_info()
    }

def run_hammock_benchmark_separated(file1_list_path: str, file2_list_path: str, mode='A', sketch_type='hyperloglog', precision=12, python_only=False) -> Dict[str, Any]:
    """Run hammock benchmark with separated creation and comparison timing."""
    
    # === PHASE 1: MEASURE SKETCH CREATION TIME ===

    
    creation_start = time.time()
    creation_cpu_start = time.process_time()
    
    try:
        # Create sketches from files
        sketches1 = []
        sketches2 = []
        
        with open(file1_list_path, 'r') as f:
            file1_list = [line.strip() for line in f if line.strip()]
        with open(file2_list_path, 'r') as f:
            file2_list = [line.strip() for line in f if line.strip()]
        
        # Create sketches for each file pair
        for file1, file2 in zip(file1_list, file2_list):
            sketch1 = IntervalSketch.from_file(
                file1, mode=mode, precision=precision, 
                sketch_type=sketch_type, use_cpp=not python_only
            )
            sketch2 = IntervalSketch.from_file(
                file2, mode=mode, precision=precision, 
                sketch_type=sketch_type, use_cpp=not python_only
            )
            sketches1.append(sketch1)
            sketches2.append(sketch2)
        
        creation_time = time.time() - creation_start
        creation_cpu_time = time.process_time() - creation_cpu_start
        
        
    except Exception as e:
        print(f"        Error creating sketches: {e}")
        return {
            'creation_cpu_time': float('nan'),
            'creation_wall_time': float('nan'),
            'comparison_cpu_time': float('nan'),
            'comparison_wall_time': float('nan'),
            'total_cpu_time': float('nan'),
            'total_wall_time': float('nan'),
            'memory_start_mb': 0,
            'memory_end_mb': 0,
            'memory_peak_mb': 0,
            'memory_delta_mb': 0,
            'error': str(e)
        }
    
    # === PHASE 2: MEASURE COMPARISON TIME ONLY ===
    
    # Initialize comparison timer
    timer = SketchComparisonTimer()
    
    # Time all pairwise comparisons
    total_comparison_time = 0.0
    total_comparison_cpu_time = 0.0
    num_comparisons = 0
    
    for i, sketch1 in enumerate(sketches1):
        for j, sketch2 in enumerate(sketches2):
            if i != j:  # Don't compare sketch with itself
                comparison_results = timer.time_jaccard_comparison(sketch1, sketch2, num_runs=1)
                if comparison_results:
                    total_comparison_time += comparison_results[0].wall_time
                    total_comparison_cpu_time += comparison_results[0].cpu_time
                    num_comparisons += 1
    
    if num_comparisons > 0:
        avg_comparison_time = total_comparison_time / num_comparisons
        avg_comparison_cpu_time = total_comparison_cpu_time / num_comparisons
    else:
        avg_comparison_time = 0.0
        avg_comparison_cpu_time = 0.0
    

    
    # Calculate total times
    total_wall_time = creation_time + total_comparison_time
    total_cpu_time = creation_cpu_time + total_comparison_cpu_time
    
    return {
        'creation_cpu_time': creation_cpu_time,
        'creation_wall_time': creation_time,
        'comparison_cpu_time': avg_comparison_cpu_time,
        'comparison_wall_time': avg_comparison_time,
        'total_cpu_time': total_cpu_time,
        'total_wall_time': total_wall_time,
        'num_comparisons': num_comparisons,
        'memory_start_mb': 0,
        'memory_end_mb': 0,
        'memory_peak_mb': 0,
        'memory_delta_mb': 0,
        'system_info': get_system_info()
    }

def run_benchmark_separated(num_files_list: List[int] = NUM_FILES_LIST, num_runs: int = NUM_RUNS):
    """Run benchmarks with separated creation and comparison timing."""
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
            print("  Running bedtools benchmark...")
            bedtools_runs = []
            for run in range(num_runs):
                run_data = run_bedtools_benchmark(file1_list_path, file2_list_path)
                bedtools_runs.append(run_data)
            
            # Extract timing and memory data for bedtools
            bedtools_cpu_times = [run['cpu_time'] for run in bedtools_runs]
            bedtools_wall_times = [run['wall_time'] for run in bedtools_runs]
            bedtools_memory_peaks = [run['memory_peak_mb'] for run in bedtools_runs]
            bedtools_memory_deltas = [run['memory_delta_mb'] for run in bedtools_runs]
            
            # Run hammock benchmarks with separated timing
            print("  Running hammock benchmarks with separated timing...")
            hammock_results = {}
            for mode in ['A', 'B', 'C']:
                for sketch_type in ['hyperloglog']:  # Focus on HyperLogLog for now
                    for precision in [8, 12, 16]:
                        print(f"    Testing {mode} mode, {sketch_type} p{precision}...")
                        # Run with C++/Cython acceleration
                        key = f"hammock_{mode}_{sketch_type}_p{precision}"
                        runs = []
                        for run in range(num_runs):
                            try:
                                run_data = run_hammock_benchmark_separated(
                                    file1_list_path, file2_list_path,
                                    mode=mode, sketch_type=sketch_type,
                                    precision=precision, python_only=False
                                )
                                runs.append(run_data)
                            except Exception as e:
                                print(f"Skipping {key} run {run+1} due to error: {e}")
                        
                        if runs:
                            # Extract separated timing data
                            creation_cpu_times = [run['creation_cpu_time'] for run in runs if not np.isnan(run['creation_cpu_time'])]
                            creation_wall_times = [run['creation_wall_time'] for run in runs if not np.isnan(run['creation_wall_time'])]
                            comparison_cpu_times = [run['comparison_cpu_time'] for run in runs if not np.isnan(run['comparison_cpu_time'])]
                            comparison_wall_times = [run['comparison_wall_time'] for run in runs if not np.isnan(run['comparison_wall_time'])]
                            total_cpu_times = [run['total_cpu_time'] for run in runs if not np.isnan(run['total_cpu_time'])]
                            total_wall_times = [run['total_wall_time'] for run in runs if not np.isnan(run['total_wall_time'])]
                            
                            hammock_results[key] = {
                                'creation_cpu_mean': np.mean(creation_cpu_times) if creation_cpu_times else float('nan'),
                                'creation_cpu_std': np.std(creation_cpu_times) if len(creation_cpu_times) > 1 else 0,
                                'creation_wall_mean': np.mean(creation_wall_times) if creation_wall_times else float('nan'),
                                'creation_wall_std': np.std(creation_wall_times) if len(creation_wall_times) > 1 else 0,
                                'comparison_cpu_mean': np.mean(comparison_cpu_times) if comparison_cpu_times else float('nan'),
                                'comparison_cpu_std': np.std(comparison_cpu_times) if len(comparison_cpu_times) > 1 else 0,
                                'comparison_wall_mean': np.mean(comparison_wall_times) if comparison_wall_times else float('nan'),
                                'comparison_wall_std': np.std(comparison_wall_times) if len(comparison_wall_times) > 1 else 0,
                                'total_cpu_mean': np.mean(total_cpu_times) if total_cpu_times else float('nan'),
                                'total_cpu_std': np.std(total_cpu_times) if len(total_cpu_times) > 1 else 0,
                                'total_wall_mean': np.mean(total_wall_times) if total_wall_times else float('nan'),
                                'total_wall_std': np.std(total_wall_times) if len(total_wall_times) > 1 else 0,
                                'detailed_runs': runs
                            }
                        else:
                            print(f"All runs failed for {key}, skipping this configuration")
                            hammock_results[key] = {
                                'creation_cpu_mean': float('nan'),
                                'creation_cpu_std': float('nan'),
                                'creation_wall_mean': float('nan'),
                                'creation_wall_std': float('nan'),
                                'comparison_cpu_mean': float('nan'),
                                'comparison_cpu_std': float('nan'),
                                'comparison_wall_mean': float('nan'),
                                'comparison_wall_std': float('nan'),
                                'total_cpu_mean': float('nan'),
                                'total_cpu_std': float('nan'),
                                'total_wall_mean': float('nan'),
                                'total_wall_std': float('nan'),
                                'detailed_runs': []
                            }
                        
                        # Run with Python-only implementation
                        print(f"      Testing Python-only implementation...")
                        key_python = f"hammock_{mode}_{sketch_type}_python_p{precision}"
                        runs_python = []
                        for run in range(num_runs):
                            try:
                                run_data = run_hammock_benchmark_separated(
                                    file1_list_path, file2_list_path,
                                    mode=mode, sketch_type=sketch_type,
                                    precision=precision, python_only=True
                                )
                                runs_python.append(run_data)
                            except Exception as e:
                                print(f"Skipping {key_python} run {run+1} due to error: {e}")
                        
                        if runs_python:
                            # Extract separated timing data for Python version
                            creation_cpu_times = [run['creation_cpu_time'] for run in runs_python if not np.isnan(run['creation_cpu_time'])]
                            creation_wall_times = [run['creation_wall_time'] for run in runs_python if not np.isnan(run['creation_wall_time'])]
                            comparison_cpu_times = [run['comparison_cpu_time'] for run in runs_python if not np.isnan(run['comparison_cpu_time'])]
                            comparison_wall_times = [run['comparison_wall_time'] for run in runs_python if not np.isnan(run['comparison_wall_time'])]
                            total_cpu_times = [run['total_cpu_time'] for run in runs_python if not np.isnan(run['total_cpu_time'])]
                            total_wall_times = [run['total_wall_time'] for run in runs_python if not np.isnan(run['total_wall_time'])]
                            
                            hammock_results[key_python] = {
                                'creation_cpu_mean': np.mean(creation_cpu_times) if creation_cpu_times else float('nan'),
                                'creation_cpu_std': np.std(creation_cpu_times) if len(creation_cpu_times) > 1 else 0,
                                'creation_wall_mean': np.mean(creation_wall_times) if creation_wall_times else float('nan'),
                                'creation_wall_std': np.std(creation_wall_times) if len(creation_wall_times) > 1 else 0,
                                'comparison_cpu_mean': np.mean(comparison_cpu_times) if comparison_cpu_times else float('nan'),
                                'comparison_cpu_std': np.std(comparison_cpu_times) if len(comparison_cpu_times) > 1 else 0,
                                'comparison_wall_mean': np.mean(comparison_wall_times) if comparison_wall_times else float('nan'),
                                'comparison_wall_std': np.std(comparison_wall_times) if len(comparison_wall_times) > 1 else 0,
                                'total_cpu_mean': np.mean(total_cpu_times) if total_cpu_times else float('nan'),
                                'total_cpu_std': np.std(total_cpu_times) if len(total_cpu_times) > 1 else 0,
                                'total_wall_mean': np.mean(total_wall_times) if total_wall_times else float('nan'),
                                'total_wall_std': np.std(total_wall_times) if len(total_wall_times) > 1 else 0,
                                'detailed_runs': runs_python
                            }
                        else:
                            print(f"All runs failed for {key_python}, skipping this configuration")
                            hammock_results[key_python] = {
                                'creation_cpu_mean': float('nan'),
                                'creation_cpu_std': float('nan'),
                                'creation_wall_mean': float('nan'),
                                'creation_wall_std': float('nan'),
                                'comparison_cpu_mean': float('nan'),
                                'comparison_cpu_std': float('nan'),
                                'comparison_wall_mean': float('nan'),
                                'comparison_wall_std': float('nan'),
                                'total_cpu_mean': float('nan'),
                                'total_cpu_std': float('nan'),
                                'total_wall_mean': float('nan'),
                                'total_wall_std': float('nan'),
                                'detailed_runs': []
                            }
            
            # Store results for this number of files
            result = {
                'num_files': num_files,
                'bedtools': {
                    'mean_cpu_time': np.mean(bedtools_cpu_times),
                    'std_cpu_time': np.std(bedtools_cpu_times) if len(bedtools_cpu_times) > 1 else 0,
                    'min_cpu_time': np.min(bedtools_cpu_times),
                    'max_cpu_time': np.max(bedtools_cpu_times),
                    'mean_wall_time': np.mean(bedtools_wall_times),
                    'std_wall_time': np.std(bedtools_wall_times) if len(bedtools_wall_times) > 1 else 0,
                    'min_wall_time': np.min(bedtools_wall_times),
                    'max_wall_time': np.max(bedtools_wall_times),
                    'mean_memory_peak': np.mean(bedtools_memory_peaks),
                    'std_memory_peak': np.std(bedtools_memory_peaks) if len(bedtools_memory_peaks) > 1 else 0,
                    'min_memory_peak': np.min(bedtools_memory_peaks),
                    'max_memory_peak': np.max(bedtools_memory_peaks),
                    'mean_memory_delta': np.mean(bedtools_memory_deltas),
                    'std_memory_delta': np.std(bedtools_memory_deltas) if len(bedtools_memory_deltas) > 1 else 0,
                    'min_memory_delta': np.min(bedtools_memory_deltas),
                    'max_memory_delta': np.max(bedtools_memory_deltas),
                    'detailed_runs': bedtools_runs
                },
                **hammock_results
            }
            
            results.append(result)
            print(f"  Completed benchmarking with {num_files} files")
    
    return results

def plot_results_separated(results, timestamp=None):
    """Plot the benchmark results with separated timing and comparison focus."""
    num_files = [r['num_files'] for r in results]
    if timestamp is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Colors and markers for different configurations
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']
    
    # Create comparison-focused graph: bedtools vs hammock comparison times
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # Plot bedtools for reference (this represents the "comparison" time for bedtools)
    bedtools_mean = [r['bedtools']['mean_cpu_time'] for r in results]
    bedtools_std = [r['bedtools']['std_cpu_time'] for r in results]
    ax.errorbar(num_files, bedtools_mean, yerr=bedtools_std, 
               fmt='ko-', capsize=5, label='bedtools (total time)', linewidth=3, markersize=10)
    
    # Plot hammock comparison times for different configurations
    color_idx = 0
    
    # Focus on Mode B HyperLogLog for comparison
    for precision in [8, 12, 16]:
        # C++/Cython version comparison times
        key = f"hammock_B_hyperloglog_p{precision}"
        if all(key in r and not np.isnan(r[key]['comparison_cpu_mean']) for r in results):
            comparison_times = [r[key]['comparison_cpu_mean'] for r in results]
            comparison_stds = [r[key]['comparison_cpu_std'] for r in results]
            label = f"hammock B HyperLogLog p{precision} (C++/Cython, comparison only)"
            ax.errorbar(num_files, comparison_times, yerr=comparison_stds,
                       fmt=f'{markers[color_idx]}:', capsize=5, 
                       color=colors[color_idx], label=label, linewidth=2, markersize=8)
            color_idx += 1
        
        # Python-only version comparison times
        key_python = f"hammock_B_hyperloglog_python_p{precision}"
        if all(key_python in r and not np.isnan(r[key_python]['comparison_cpu_mean']) for r in results):
            comparison_times = [r[key_python]['comparison_cpu_mean'] for r in results]
            comparison_stds = [r[key_python]['comparison_cpu_std'] for r in results]
            label = f"hammock B HyperLogLog p{precision} (Python, comparison only)"
            ax.errorbar(num_files, comparison_times, yerr=comparison_stds,
                       fmt=f'{markers[color_idx-1]}-', capsize=5, 
                       color=colors[color_idx-1], label=label, linewidth=2, markersize=8)
    
    ax.set_xlabel('Number of Files', fontsize=12)
    ax.set_ylabel('CPU Time (seconds)', fontsize=12)
    ax.set_title('Comparison Time Analysis: bedtools vs hammock Sketch Comparison Operations', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    ax.set_yscale('log')  # Use log scale to better show the differences
    
    plt.tight_layout()
    
    # Save comparison plot
    plot_path = os.path.join('benchmarks/results', f'comparison_times_analysis_{timestamp}.png')
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    
    # Create creation vs comparison breakdown graph
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Left plot: Creation times
    ax1.set_title('Sketch Creation Times', fontsize=14)
    color_idx = 0
    
    for precision in [8, 12, 16]:
        # C++/Cython version creation times
        key = f"hammock_B_hyperloglog_p{precision}"
        if all(key in r and not np.isnan(r[key]['creation_cpu_mean']) for r in results):
            creation_times = [r[key]['creation_cpu_mean'] for r in results]
            creation_stds = [r[key]['creation_cpu_std'] for r in results]
            label = f"hammock B HyperLogLog p{precision} (C++/Cython)"
            ax1.errorbar(num_files, creation_times, yerr=creation_stds,
                        fmt=f'{markers[color_idx]}:', capsize=5, 
                        color=colors[color_idx], label=label, linewidth=2)
            color_idx += 1
        
        # Python-only version creation times
        key_python = f"hammock_B_hyperloglog_python_p{precision}"
        if all(key_python in r and not np.isnan(r[key_python]['creation_cpu_mean']) for r in results):
            creation_times = [r[key_python]['creation_cpu_mean'] for r in results]
            creation_stds = [r[key_python]['creation_cpu_std'] for r in results]
            label = f"hammock B HyperLogLog p{precision} (Python)"
            ax1.errorbar(num_files, creation_times, yerr=creation_stds,
                        fmt=f'{markers[color_idx-1]}-', capsize=5, 
                        color=colors[color_idx-1], label=label, linewidth=2)
    
    ax1.set_xlabel('Number of Files', fontsize=12)
    ax1.set_ylabel('CPU Time (seconds)', fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    ax1.set_yscale('log')
    
    # Right plot: Comparison times
    ax2.set_title('Sketch Comparison Times', fontsize=14)
    color_idx = 0
    
    for precision in [8, 12, 16]:
        # C++/Cython version comparison times
        key = f"hammock_B_hyperloglog_p{precision}"
        if all(key in r and not np.isnan(r[key]['comparison_cpu_mean']) for r in results):
            comparison_times = [r[key]['comparison_cpu_mean'] for r in results]
            comparison_stds = [r[key]['comparison_cpu_std'] for r in results]
            label = f"hammock B HyperLogLog p{precision} (C++/Cython)"
            ax2.errorbar(num_files, comparison_times, yerr=comparison_stds,
                        fmt=f'{markers[color_idx]}:', capsize=5, 
                        color=colors[color_idx], label=label, linewidth=2)
            color_idx += 1
        
        # Python-only version comparison times
        key_python = f"hammock_B_hyperloglog_python_p{precision}"
        if all(key_python in r and not np.isnan(r[key_python]['comparison_cpu_mean']) for r in results):
            comparison_times = [r[key_python]['comparison_cpu_mean'] for r in results]
            comparison_stds = [r[key_python]['comparison_cpu_std'] for r in results]
            label = f"hammock B HyperLogLog p{precision} (Python)"
            ax2.errorbar(num_files, comparison_times, yerr=comparison_stds,
                        fmt=f'{markers[color_idx-1]}-', capsize=5, 
                        color=colors[color_idx-1], label=label, linewidth=2)
    
    ax2.set_xlabel('Number of Files', fontsize=12)
    ax2.set_ylabel('CPU Time (seconds)', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    ax2.set_yscale('log')
    
    plt.tight_layout()
    
    # Save breakdown plot
    plot_path = os.path.join('benchmarks/results', f'creation_vs_comparison_breakdown_{timestamp}.png')
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    
    # Save numerical results
    results_path = os.path.join('benchmarks/results', f'separated_timing_results_{timestamp}.txt')
    with open(results_path, 'w') as f:
        f.write("Separated Timing Benchmark Results:\n")
        f.write("==================================\n\n")
        for r in results:
            f.write(f"Number of files: {r['num_files']}\n")
            f.write("\nbedtools:\n")
            f.write(f"  Total CPU Time - Mean: {r['bedtools']['mean_cpu_time']:.3f} ± {r['bedtools']['std_cpu_time']:.3f} seconds\n")
            f.write(f"  Total Wall Time - Mean: {r['bedtools']['mean_wall_time']:.3f} ± {r['bedtools']['std_wall_time']:.3f} seconds\n")
            f.write("\nhammock (separated timing):\n")
            for key in r:
                if key.startswith('hammock_'):
                    f.write(f"\n  {key}:\n")
                    f.write(f"    Creation CPU Time - Mean: {r[key]['creation_cpu_mean']:.6f} ± {r[key]['creation_cpu_std']:.6f} seconds\n")
                    f.write(f"    Creation Wall Time - Mean: {r[key]['creation_wall_mean']:.6f} ± {r[key]['creation_wall_std']:.6f} seconds\n")
                    f.write(f"    Comparison CPU Time - Mean: {r[key]['comparison_cpu_mean']:.6f} ± {r[key]['comparison_cpu_std']:.6f} seconds\n")
                    f.write(f"    Comparison Wall Time - Mean: {r[key]['comparison_wall_mean']:.6f} ± {r[key]['comparison_wall_std']:.6f} seconds\n")
                    f.write(f"    Total CPU Time - Mean: {r[key]['total_cpu_mean']:.6f} ± {r[key]['total_cpu_std']:.6f} seconds\n")
                    f.write(f"    Total Wall Time - Mean: {r[key]['total_wall_mean']:.6f} ± {r[key]['total_wall_std']:.6f} seconds\n")
            f.write("\n" + "="*50 + "\n\n")
    
    # Save results as CSV table
    csv_path = os.path.join('benchmarks/results', f'separated_timing_results_{timestamp}.csv')
    save_separated_results_table(results, csv_path)
    
    print(f"\nSeparated timing analysis completed!")
    print(f"Results saved in benchmarks/results/")
    print(f"Key findings:")
    print(f"- Comparison-focused graph: comparison_times_analysis_{timestamp}.png")
    print(f"- Creation vs comparison breakdown: creation_vs_comparison_breakdown_{timestamp}.png")
    print(f"- Detailed results: separated_timing_results_{timestamp}.txt")
    print(f"- CSV table: separated_timing_results_{timestamp}.csv")

def save_separated_results_table(results, output_path):
    """Save separated timing benchmark results as a CSV table."""
    import csv
    
    # Create a comprehensive table with all metrics including separated timing
    columns = ['num_files', 'tool', 'mode', 'sketch_type', 'precision', 'implementation',
               'cpu_mean', 'cpu_std', 'cpu_min', 'cpu_max',
               'wall_mean', 'wall_std', 'wall_min', 'wall_max',
               'creation_cpu_mean', 'creation_cpu_std', 'creation_cpu_min', 'creation_cpu_max',
               'creation_wall_mean', 'creation_wall_std', 'creation_wall_min', 'creation_wall_max',
               'comparison_cpu_mean', 'comparison_cpu_std', 'comparison_cpu_min', 'comparison_cpu_max',
               'comparison_wall_mean', 'comparison_wall_std', 'comparison_wall_min', 'comparison_wall_max',
               'memory_peak_mean', 'memory_peak_std', 'memory_peak_min', 'memory_peak_max',
               'memory_delta_mean', 'memory_delta_std', 'memory_delta_min', 'memory_delta_max']
    
    # Write CSV file
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=columns)
        writer.writeheader()
        
        for r in results:
            # Add bedtools results (no separated timing available)
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
                'creation_cpu_mean': 'N/A',
                'creation_cpu_std': 'N/A',
                'creation_cpu_min': 'N/A',
                'creation_cpu_max': 'N/A',
                'creation_wall_mean': 'N/A',
                'creation_wall_std': 'N/A',
                'creation_wall_min': 'N/A',
                'creation_wall_max': 'N/A',
                'comparison_cpu_mean': 'N/A',
                'comparison_cpu_std': 'N/A',
                'comparison_cpu_min': 'N/A',
                'comparison_cpu_max': 'N/A',
                'comparison_wall_mean': 'N/A',
                'comparison_wall_std': 'N/A',
                'comparison_wall_min': 'N/A',
                'comparison_wall_max': 'N/A',
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
            
            # Add hammock results with separated timing
            for key in r:
                if key.startswith('hammock_'):
                    # Parse the key to extract mode, sketch_type, precision, and implementation
                    parts = key.split('_')
                    if len(parts) >= 4:
                        mode = parts[1]
                        sketch_type = parts[2]
                        
                        # Handle both formats: hammock_A_hyperloglog_p8 and hammock_A_hyperloglog_python_p8
                        if 'python' in key:
                            # For Python implementations: hammock_A_hyperloglog_python_p8
                            implementation = 'Python'
                            # Find the part that starts with 'p' (precision)
                            precision = 'N/A'
                            for part in parts:
                                if part.startswith('p') and part[1:].isdigit():
                                    precision = part[1:]
                                    break
                        else:
                            # For C++/Cython implementations: hammock_A_hyperloglog_p8
                            implementation = 'C++/Cython'
                            precision_part = parts[3]
                            if precision_part.startswith('p'):
                                precision = precision_part[1:]
                            else:
                                precision = 'N/A'
                        
                        # Get the separated timing data
                        data = r[key]
                        detailed_runs = data.get('detailed_runs', [])
                        
                        if detailed_runs:
                            # Extract timing data from detailed runs
                            creation_cpu_times = [run.get('creation_cpu_time', 0) for run in detailed_runs if not np.isnan(run.get('creation_cpu_time', 0))]
                            creation_wall_times = [run.get('creation_wall_time', 0) for run in detailed_runs if not np.isnan(run.get('creation_wall_time', 0))]
                            comparison_cpu_times = [run.get('comparison_cpu_time', 0) for run in detailed_runs if not np.isnan(run.get('comparison_cpu_time', 0))]
                            comparison_wall_times = [run.get('comparison_wall_time', 0) for run in detailed_runs if not np.isnan(run.get('comparison_wall_time', 0))]
                            
                            row = {
                                'num_files': r['num_files'],
                                'tool': 'hammock',
                                'mode': mode,
                                'sketch_type': sketch_type,
                                'precision': precision,
                                'implementation': implementation,
                                'cpu_mean': f"{data.get('total_cpu_mean', 0):.6f}",
                                'cpu_std': f"{data.get('total_cpu_std', 0):.6f}",
                                'cpu_min': f"{data.get('total_cpu_mean', 0):.6f}",  # Approximate
                                'cpu_max': f"{data.get('total_cpu_mean', 0):.6f}",  # Approximate
                                'wall_mean': f"{data.get('total_wall_mean', 0):.6f}",
                                'wall_std': f"{data.get('total_wall_std', 0):.6f}",
                                'wall_min': f"{data.get('total_wall_mean', 0):.6f}",  # Approximate
                                'wall_max': f"{data.get('total_wall_mean', 0):.6f}",  # Approximate
                                'creation_cpu_mean': f"{data.get('creation_cpu_mean', 0):.6f}",
                                'creation_cpu_std': f"{data.get('creation_cpu_std', 0):.6f}",
                                'creation_cpu_min': f"{data.get('creation_cpu_mean', 0):.6f}",  # Approximate
                                'creation_cpu_max': f"{data.get('creation_cpu_mean', 0):.6f}",  # Approximate
                                'creation_wall_mean': f"{data.get('creation_wall_mean', 0):.6f}",
                                'creation_wall_std': f"{data.get('creation_wall_std', 0):.6f}",
                                'creation_wall_min': f"{data.get('creation_wall_mean', 0):.6f}",  # Approximate
                                'creation_wall_max': f"{data.get('creation_wall_mean', 0):.6f}",  # Approximate
                                'comparison_cpu_mean': f"{data.get('comparison_cpu_mean', 0):.6f}",
                                'comparison_cpu_std': f"{data.get('comparison_cpu_std', 0):.6f}",
                                'comparison_cpu_min': f"{data.get('comparison_cpu_mean', 0):.6f}",  # Approximate
                                'comparison_cpu_max': f"{data.get('comparison_cpu_mean', 0):.6f}",  # Approximate
                                'comparison_wall_mean': f"{data.get('comparison_wall_mean', 0):.6f}",
                                'comparison_wall_std': f"{data.get('comparison_wall_std', 0):.6f}",
                                'comparison_wall_min': f"{data.get('comparison_wall_mean', 0):.6f}",  # Approximate
                                'comparison_wall_max': f"{data.get('comparison_wall_mean', 0):.6f}",  # Approximate
                                'memory_peak_mean': '0.000000',
                                'memory_peak_std': '0.000000',
                                'memory_peak_min': '0.000000',
                                'memory_peak_max': '0.000000',
                                'memory_delta_mean': '0.000000',
                                'memory_delta_std': '0.000000',
                                'memory_delta_min': '0.000000',
                                'memory_delta_max': '0.000000'
                            }
                            writer.writerow(row)
    
    print(f"CSV table saved to: {output_path}")

def quick_test():
    """Run a quick test to verify separated timing functionality."""
    print("Running quick test to verify separated timing functionality...")
    
    # Create results directory if it doesn't exist
    os.makedirs('benchmarks/results', exist_ok=True)
    
    # Generate timestamp for this run
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Run with minimal parameters for quick test
    results = run_benchmark_separated(num_files_list=[2, 4], num_runs=2)
    plot_results_separated(results, timestamp)
    print("\nQuick test completed. Check benchmarks/results/ for generated graphs.")

def main():
    # Create results directory if it doesn't exist
    os.makedirs('benchmarks/results', exist_ok=True)
    
    # Generate timestamp for this run
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    print("Starting separated timing benchmark...")
    print("This benchmark separates sketch creation time from comparison time")
    print("to provide more accurate performance analysis.\n")
    
    results = run_benchmark_separated()
    plot_results_separated(results, timestamp)
    print("\nSeparated timing benchmark completed. Results saved in benchmarks/results/")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--test":
        quick_test()
    else:
        main()
