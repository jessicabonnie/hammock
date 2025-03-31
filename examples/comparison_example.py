#!/usr/bin/env python3
"""
Comparison of HyperLogLog implementations in Hammock.

This example compares the performance and accuracy of three implementations:
1. HyperLogLog (Pure Python)
2. FastHyperLogLog (Python with Rust via Py03)
3. RustHLLWrapper (Native Rust implementation)
"""

import time
import random
import numpy as np
import matplotlib.pyplot as plt
from hammock import HyperLogLog, FastHyperLogLog

# Check if RustHLLWrapper is available
try:
    from hammock import RustHLLWrapper
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False
    print("Native Rust implementation not available.")
    print("Build it with: cd rust_hll && python -m maturin develop")

def generate_data(size, unique_ratio=1.0):
    """Generate data with controlled uniqueness."""
    unique_count = int(size * unique_ratio)
    data = []
    for _ in range(size):
        data.append(f"item_{random.randint(0, unique_count)}")
    return data

def benchmark_add(implementations, data):
    """Benchmark adding items one by one."""
    results = {}
    
    for name, sketch in implementations.items():
        start_time = time.time()
        
        if name == "RustHLLWrapper (Native)":
            for item in data:
                sketch.add(item)
        else:
            for item in data:
                sketch.add_string(item)
                
        elapsed = time.time() - start_time
        results[name] = elapsed
        
        print(f"{name}: Added {len(data)} items in {elapsed:.4f}s")
    
    return results

def benchmark_batch(implementations, data):
    """Benchmark adding items in batch."""
    results = {}
    
    for name, sketch in implementations.items():
        if not hasattr(sketch, "add_batch"):
            results[name] = float('inf')
            print(f"{name}: Batch add not supported")
            continue
            
        start_time = time.time()
        sketch.add_batch(data)
        elapsed = time.time() - start_time
        results[name] = elapsed
        
        print(f"{name}: Batch added {len(data)} items in {elapsed:.4f}s")
    
    return results

def benchmark_estimation(implementations, actual_unique):
    """Benchmark cardinality estimation."""
    results = {}
    
    for name, sketch in implementations.items():
        start_time = time.time()
        
        if name == "RustHLLWrapper (Native)":
            estimate = sketch.cardinality()
        else:
            estimate = sketch.estimate_cardinality(method="ertl_mle")
            
        elapsed = time.time() - start_time
        error = abs(estimate - actual_unique) / actual_unique * 100
        
        results[name] = {
            "time": elapsed,
            "estimate": estimate,
            "error": error
        }
        
        print(f"{name}: Estimated {estimate:.1f} (actual: {actual_unique}) in {elapsed:.6f}s, error: {error:.2f}%")
    
    return results

def plot_results(title, results, ylabel, log_scale=False):
    """Plot benchmark results."""
    plt.figure(figsize=(10, 6))
    
    names = list(results.keys())
    values = list(results.values())
    
    # Sort by value
    sorted_indices = np.argsort(values)
    sorted_names = [names[i] for i in sorted_indices]
    sorted_values = [values[i] for i in sorted_indices]
    
    # Plot
    bars = plt.bar(sorted_names, sorted_values)
    
    # Add values on top of bars
    for i, bar in enumerate(bars):
        value = sorted_values[i]
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + (max(sorted_values) * 0.01),
            f"{value:.4f}",
            ha='center'
        )
    
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xticks(rotation=15, ha='right')
    plt.tight_layout()
    
    if log_scale:
        plt.yscale('log')
    
    plt.savefig(f"{title.lower().replace(' ', '_')}.png")
    plt.show()

def main():
    """Run the benchmark comparison."""
    print("Hammock HyperLogLog Implementations Comparison")
    print("=" * 45)
    
    # Create implementations
    implementations = {
        "HyperLogLog (Python)": HyperLogLog(precision=12),
        "FastHyperLogLog (Py03)": FastHyperLogLog(precision=12)
    }
    
    if RUST_AVAILABLE:
        implementations["RustHLLWrapper (Native)"] = RustHLLWrapper(precision=12)
    
    # Generate data
    print("\nGenerating test data...")
    data_size = 100000
    data = generate_data(data_size, unique_ratio=0.8)
    unique_count = len(set(data))
    print(f"Generated {data_size} items with {unique_count} unique values")
    
    # Benchmark adding items
    print("\nBenchmarking adding items one by one...")
    add_results = benchmark_add(implementations, data[:10000])  # Use smaller set for single adds
    
    # Benchmark batch adding
    print("\nBenchmarking batch adding...")
    batch_results = benchmark_batch(implementations, data)
    
    # Benchmark estimation
    print("\nBenchmarking cardinality estimation...")
    estimate_results = benchmark_estimation(implementations, unique_count)
    
    # Extract estimation times
    estimation_times = {name: results["time"] for name, results in estimate_results.items()}
    
    # Plot results
    if not plt.isinteractive():
        # Only generate plots if we can display them
        plot_results("Item Addition Time (seconds)", add_results, "Time (s)")
        plot_results("Batch Addition Time (seconds)", batch_results, "Time (s)")
        plot_results("Estimation Time (seconds)", estimation_times, "Time (s)")
    
    # Print summary
    print("\nPerformance Summary:")
    print("=" * 45)
    
    fastest_add = min(add_results, key=add_results.get)
    fastest_batch = min(batch_results, key=batch_results.get) if any(v != float('inf') for v in batch_results.values()) else "N/A"
    fastest_estimate = min(estimation_times, key=estimation_times.get)
    
    most_accurate = min(estimate_results, key=lambda k: estimate_results[k]["error"])
    
    print(f"Fastest for single adds: {fastest_add}")
    print(f"Fastest for batch adds: {fastest_batch}")
    print(f"Fastest for estimation: {fastest_estimate}")
    print(f"Most accurate: {most_accurate} (error: {estimate_results[most_accurate]['error']:.2f}%)")
    
    # Speed comparison
    if RUST_AVAILABLE:
        rust_add = add_results["RustHLLWrapper (Native)"]
        python_add = add_results["HyperLogLog (Python)"]
        speedup = python_add / rust_add
        print(f"\nRust is {speedup:.1f}x faster than Python for single adds")
        
        if "RustHLLWrapper (Native)" in batch_results and batch_results["RustHLLWrapper (Native)"] != float('inf'):
            rust_batch = batch_results["RustHLLWrapper (Native)"]
            if "FastHyperLogLog (Py03)" in batch_results and batch_results["FastHyperLogLog (Py03)"] != float('inf'):
                py03_batch = batch_results["FastHyperLogLog (Py03)"]
                batch_speedup = py03_batch / rust_batch
                print(f"Rust is {batch_speedup:.1f}x faster than Py03 for batch adds")

if __name__ == "__main__":
    main() 