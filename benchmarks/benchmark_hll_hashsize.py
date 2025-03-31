#!/usr/bin/env python
import time
import random
import numpy as np
from typing import List, Tuple, Dict
import matplotlib.pyplot as plt
from datetime import datetime
import os
import sys
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.rusthll import RustHyperLogLog

# Create results directory if it doesn't exist
RESULTS_DIR = os.path.join(os.path.dirname(__file__), 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)

# Get script name without extension and create timestamp
SCRIPT_NAME = os.path.splitext(os.path.basename(__file__))[0]
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

# Create output file paths with script name and timestamp
DEBUG_LOG = os.path.join(RESULTS_DIR, f"{SCRIPT_NAME}_{TIMESTAMP}.log")
RESULTS_PNG = os.path.join(RESULTS_DIR, f"{SCRIPT_NAME}_{TIMESTAMP}.png")
RESULTS_PNG_PRECISION = os.path.join(RESULTS_DIR, f"{SCRIPT_NAME}_{TIMESTAMP}_precision.png")
RESULTS_TXT = os.path.join(RESULTS_DIR, f"{SCRIPT_NAME}_{TIMESTAMP}.txt")

# Redirect stdout to the debug log
sys.stdout = open(DEBUG_LOG, 'w')

# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)

# Constants for benchmarking
PRECISION_VALUES = [6, 8, 10, 12, 14, 16, 20, 24, 28, 30]
DATA_SIZES = [1000, 10000, 100000, 1000000]
BATCH_SIZES = [100, 1000, 10000]
NUM_RUNS = 5

def generate_test_data(size: int) -> List[str]:
    """Generate test data of specified size."""
    return [f"item_{i}" for i in range(size)]

def benchmark_add_string(sketch, data: List[str]) -> float:
    """Benchmark adding values one at a time."""
    start_time = time.time()
    for item in data:
        sketch.add_string(item)
    end_time = time.time()
    return end_time - start_time

def benchmark_add_batch(sketch, data: List[str], batch_size: int = 1000) -> float:
    """Benchmark adding values in batches."""
    start_time = time.time()
    for i in range(0, len(data), batch_size):
        batch = data[i:i + batch_size]
        sketch.add_batch(batch)
    end_time = time.time()
    return end_time - start_time

def benchmark_merge(sketch1, sketch2) -> float:
    """Benchmark merging two sketches."""
    start_time = time.time()
    sketch1.merge(sketch2)
    end_time = time.time()
    return end_time - start_time

def benchmark_estimate(sketch) -> float:
    """Benchmark cardinality estimation."""
    start_time = time.time()
    sketch.estimate_cardinality()
    end_time = time.time()
    return end_time - start_time

def run_benchmarks(data_sizes: List[int] = DATA_SIZES, 
                  batch_sizes: List[int] = BATCH_SIZES,
                  precision_values: List[int] = PRECISION_VALUES,
                  num_runs: int = NUM_RUNS) -> Dict:
    """Run comprehensive benchmarks."""
    results = {
        'add_string': [],
        'add_batch': [],
        'merge': [],
        'estimate': [],
        'accuracy': []
    }
    
    for precision in precision_values:
        print(f"\nBenchmarking with precision {precision}...")
        
        for size in data_sizes:
            print(f"  Testing with {size} items...")
            data = generate_test_data(size)
            
            # Create sketches
            sketches = {
                'py_32': HyperLogLog(precision=precision, hash_size=32),
                'py_64': HyperLogLog(precision=precision, hash_size=64),
                'rust_32': RustHyperLogLog(precision=precision, hash_size=32),
                'rust_64': RustHyperLogLog(precision=precision, hash_size=64)
            }
            
            # Benchmark add_string
            for name, sketch in sketches.items():
                times = []
                for _ in range(num_runs):
                    # Create a fresh sketch for each run
                    if name.startswith('py'):
                        sketch = HyperLogLog(precision=precision, hash_size=32 if name.endswith('32') else 64)
                    else:
                        sketch = RustHyperLogLog(precision=precision, hash_size=32 if name.endswith('32') else 64)
                    time_taken = benchmark_add_string(sketch, data)
                    times.append(time_taken)
                results['add_string'].append({
                    'size': size,
                    'precision': precision,
                    'implementation': name,
                    'mean_time': np.mean(times),
                    'std_time': np.std(times)
                })
            
            # Benchmark add_batch
            for batch_size in batch_sizes:
                for name, sketch in sketches.items():
                    times = []
                    for _ in range(num_runs):
                        if name.startswith('py'):
                            sketch = HyperLogLog(precision=precision, hash_size=32 if name.endswith('32') else 64)
                        else:
                            sketch = RustHyperLogLog(precision=precision, hash_size=32 if name.endswith('32') else 64)
                        time_taken = benchmark_add_batch(sketch, data, batch_size)
                        times.append(time_taken)
                    results['add_batch'].append({
                        'size': size,
                        'precision': precision,
                        'batch_size': batch_size,
                        'implementation': name,
                        'mean_time': np.mean(times),
                        'std_time': np.std(times)
                    })
            
            # Benchmark merge
            for name, sketch in sketches.items():
                times = []
                for _ in range(num_runs):
                    if name.startswith('py'):
                        sketch1 = HyperLogLog(precision=precision, hash_size=32 if name.endswith('32') else 64)
                        sketch2 = HyperLogLog(precision=precision, hash_size=32 if name.endswith('32') else 64)
                    else:
                        sketch1 = RustHyperLogLog(precision=precision, hash_size=32 if name.endswith('32') else 64)
                        sketch2 = RustHyperLogLog(precision=precision, hash_size=32 if name.endswith('32') else 64)
                    
                    # Add data to both sketches
                    for item in data[:size//2]:
                        sketch1.add_string(item)
                    for item in data[size//2:]:
                        sketch2.add_string(item)
                    
                    time_taken = benchmark_merge(sketch1, sketch2)
                    times.append(time_taken)
                results['merge'].append({
                    'size': size,
                    'precision': precision,
                    'implementation': name,
                    'mean_time': np.mean(times),
                    'std_time': np.std(times)
                })
            
            # Benchmark estimate and measure accuracy
            for name, sketch in sketches.items():
                times = []
                estimates = []
                for _ in range(num_runs):
                    if name.startswith('py'):
                        sketch = HyperLogLog(precision=precision, hash_size=32 if name.endswith('32') else 64)
                    else:
                        sketch = RustHyperLogLog(precision=precision, hash_size=32 if name.endswith('32') else 64)
                    
                    # Add data
                    for item in data:
                        sketch.add_string(item)
                    
                    time_taken = benchmark_estimate(sketch)
                    estimate = sketch.estimate_cardinality()
                    
                    times.append(time_taken)
                    estimates.append(estimate)
                
                results['estimate'].append({
                    'size': size,
                    'precision': precision,
                    'implementation': name,
                    'mean_time': np.mean(times),
                    'std_time': np.std(times)
                })
                
                results['accuracy'].append({
                    'size': size,
                    'precision': precision,
                    'implementation': name,
                    'mean_estimate': np.mean(estimates),
                    'std_estimate': np.std(estimates),
                    'relative_error': abs(np.mean(estimates) - size) / size
                })
    
    return results

def plot_results(results: Dict):
    """Plot benchmark results."""
    # Plot 1: Performance by data size
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('HyperLogLog Implementation Comparison (Performance by Data Size)')
    
    # Plot add_string times
    ax = axes[0, 0]
    for impl in ['py_32', 'py_64', 'rust_32', 'rust_64']:
        data = [r for r in results['add_string'] if r['implementation'] == impl and r['precision'] == 14]
        sizes = [r['size'] for r in data]
        times = [r['mean_time'] for r in data]
        stds = [r['std_time'] for r in data]
        ax.errorbar(sizes, times, yerr=stds, label=impl)
    ax.set_xlabel('Number of Items')
    ax.set_ylabel('Time (seconds)')
    ax.set_title('Add Value Performance (precision=14)')
    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Plot add_batch times for largest batch size
    ax = axes[0, 1]
    batch_size = max(r['batch_size'] for r in results['add_batch'])
    for impl in ['py_32', 'py_64', 'rust_32', 'rust_64']:
        data = [r for r in results['add_batch'] 
                if r['implementation'] == impl and r['batch_size'] == batch_size and r['precision'] == 14]
        sizes = [r['size'] for r in data]
        times = [r['mean_time'] for r in data]
        stds = [r['std_time'] for r in data]
        ax.errorbar(sizes, times, yerr=stds, label=impl)
    ax.set_xlabel('Number of Items')
    ax.set_ylabel('Time (seconds)')
    ax.set_title(f'Add Batch Performance (precision=14, batch_size={batch_size})')
    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Plot merge times
    ax = axes[1, 0]
    for impl in ['py_32', 'py_64', 'rust_32', 'rust_64']:
        data = [r for r in results['merge'] if r['implementation'] == impl and r['precision'] == 14]
        sizes = [r['size'] for r in data]
        times = [r['mean_time'] for r in data]
        stds = [r['std_time'] for r in data]
        ax.errorbar(sizes, times, yerr=stds, label=impl)
    ax.set_xlabel('Number of Items')
    ax.set_ylabel('Time (seconds)')
    ax.set_title('Merge Performance (precision=14)')
    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Plot relative errors
    ax = axes[1, 1]
    for impl in ['py_32', 'py_64', 'rust_32', 'rust_64']:
        data = [r for r in results['accuracy'] if r['implementation'] == impl and r['precision'] == 14]
        sizes = [r['size'] for r in data]
        errors = [r['relative_error'] for r in data]
        ax.plot(sizes, errors, label=impl)
    ax.set_xlabel('Number of Items')
    ax.set_ylabel('Relative Error')
    ax.set_title('Estimation Accuracy (precision=14)')
    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(RESULTS_PNG)
    plt.close()
    
    # Plot 2: Performance by precision
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('HyperLogLog Implementation Comparison (Performance by Precision)')
    
    # Plot add_string times by precision
    ax = axes[0, 0]
    size = max(r['size'] for r in results['add_string'])
    for impl in ['py_32', 'py_64', 'rust_32', 'rust_64']:
        data = [r for r in results['add_string'] if r['implementation'] == impl and r['size'] == size]
        precisions = [r['precision'] for r in data]
        times = [r['mean_time'] for r in data]
        stds = [r['std_time'] for r in data]
        ax.errorbar(precisions, times, yerr=stds, label=impl)
    ax.set_xlabel('Precision')
    ax.set_ylabel('Time (seconds)')
    ax.set_title(f'Add Value Performance (size={size})')
    ax.legend()
    
    # Plot add_batch times by precision
    ax = axes[0, 1]
    batch_size = max(r['batch_size'] for r in results['add_batch'])
    for impl in ['py_32', 'py_64', 'rust_32', 'rust_64']:
        data = [r for r in results['add_batch'] 
                if r['implementation'] == impl and r['size'] == size and r['batch_size'] == batch_size]
        precisions = [r['precision'] for r in data]
        times = [r['mean_time'] for r in data]
        stds = [r['std_time'] for r in data]
        ax.errorbar(precisions, times, yerr=stds, label=impl)
    ax.set_xlabel('Precision')
    ax.set_ylabel('Time (seconds)')
    ax.set_title(f'Add Batch Performance (size={size}, batch_size={batch_size})')
    ax.legend()
    
    # Plot merge times by precision
    ax = axes[1, 0]
    for impl in ['py_32', 'py_64', 'rust_32', 'rust_64']:
        data = [r for r in results['merge'] if r['implementation'] == impl and r['size'] == size]
        precisions = [r['precision'] for r in data]
        times = [r['mean_time'] for r in data]
        stds = [r['std_time'] for r in data]
        ax.errorbar(precisions, times, yerr=stds, label=impl)
    ax.set_xlabel('Precision')
    ax.set_ylabel('Time (seconds)')
    ax.set_title(f'Merge Performance (size={size})')
    ax.legend()
    
    # Plot relative errors by precision
    ax = axes[1, 1]
    for impl in ['py_32', 'py_64', 'rust_32', 'rust_64']:
        data = [r for r in results['accuracy'] if r['implementation'] == impl and r['size'] == size]
        precisions = [r['precision'] for r in data]
        errors = [r['relative_error'] for r in data]
        ax.plot(precisions, errors, label=impl)
    ax.set_xlabel('Precision')
    ax.set_ylabel('Relative Error')
    ax.set_title(f'Estimation Accuracy (size={size})')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(RESULTS_PNG_PRECISION)
    plt.close()

def print_results(results: Dict):
    """Print benchmark results in a readable format."""
    # Create a text file for the results
    with open(RESULTS_TXT, 'w') as f:
        f.write("\nBenchmark Results:\n")
        f.write("=" * 80 + "\n")
        
        # Print add_string results
        f.write("\nAdd Value Performance:\n")
        f.write("-" * 80 + "\n")
        for precision in sorted(set(r['precision'] for r in results['add_string'])):
            f.write(f"\nPrecision: {precision}\n")
            for size in sorted(set(r['size'] for r in results['add_string'])):
                f.write(f"\n  Data size: {size}\n")
                for impl in ['py_32', 'py_64', 'rust_32', 'rust_64']:
                    data = [r for r in results['add_string'] 
                           if r['implementation'] == impl and r['size'] == size and r['precision'] == precision][0]
                    f.write(f"    {impl:8s}: {data['mean_time']:.6f} ± {data['std_time']:.6f} seconds\n")
        
        # Print add_batch results
        f.write("\nAdd Batch Performance:\n")
        f.write("-" * 80 + "\n")
        for precision in sorted(set(r['precision'] for r in results['add_batch'])):
            f.write(f"\nPrecision: {precision}\n")
            for size in sorted(set(r['size'] for r in results['add_batch'])):
                f.write(f"\n  Data size: {size}\n")
                for batch_size in sorted(set(r['batch_size'] for r in results['add_batch'])):
                    f.write(f"    Batch size: {batch_size}\n")
                    for impl in ['py_32', 'py_64', 'rust_32', 'rust_64']:
                        data = [r for r in results['add_batch'] 
                               if r['implementation'] == impl and r['size'] == size 
                               and r['batch_size'] == batch_size and r['precision'] == precision][0]
                        f.write(f"      {impl:8s}: {data['mean_time']:.6f} ± {data['std_time']:.6f} seconds\n")
        
        # Print merge results
        f.write("\nMerge Performance:\n")
        f.write("-" * 80 + "\n")
        for precision in sorted(set(r['precision'] for r in results['merge'])):
            f.write(f"\nPrecision: {precision}\n")
            for size in sorted(set(r['size'] for r in results['merge'])):
                f.write(f"\n  Data size: {size}\n")
                for impl in ['py_32', 'py_64', 'rust_32', 'rust_64']:
                    data = [r for r in results['merge'] 
                           if r['implementation'] == impl and r['size'] == size and r['precision'] == precision][0]
                    f.write(f"    {impl:8s}: {data['mean_time']:.6f} ± {data['std_time']:.6f} seconds\n")
        
        # Print accuracy results
        f.write("\nEstimation Accuracy:\n")
        f.write("-" * 80 + "\n")
        for precision in sorted(set(r['precision'] for r in results['accuracy'])):
            f.write(f"\nPrecision: {precision}\n")
            for size in sorted(set(r['size'] for r in results['accuracy'])):
                f.write(f"\n  Data size: {size}\n")
                for impl in ['py_32', 'py_64', 'rust_32', 'rust_64']:
                    data = [r for r in results['accuracy'] 
                           if r['implementation'] == impl and r['size'] == size and r['precision'] == precision][0]
                    f.write(f"    {impl:8s}: {data['mean_estimate']:.2f} ± {data['std_estimate']:.2f} "
                           f"(relative error: {data['relative_error']:.4f})\n")
    
    # Also print to console
    with open(RESULTS_TXT, 'r') as f:
        print(f.read())

def main():
    """Main entry point."""
    print("\nNOTE: Before running benchmarks, it is recommended to run the tests first:")
    print("      python -m pytest tests/test_hll_hashsize.py -v\n")
    
    print(f"Output files will be saved with prefix: {SCRIPT_NAME}_{TIMESTAMP}")
    print(f"Results directory: {RESULTS_DIR}\n")
    
    # Run benchmarks
    results = run_benchmarks()
    
    # Plot results
    plot_results(results)
    
    # Print results
    print_results(results)
    
    print(f"\nBenchmark results saved in directory: {RESULTS_DIR}")
    print(f"Debug log: {DEBUG_LOG}")
    print(f"Results plot (by size): {RESULTS_PNG}")
    print(f"Results plot (by precision): {RESULTS_PNG_PRECISION}")
    print(f"Results text: {RESULTS_TXT}")

if __name__ == "__main__":
    main() 