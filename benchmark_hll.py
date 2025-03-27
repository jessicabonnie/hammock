#!/usr/bin/env python3
import time
import random
import numpy as np
import matplotlib.pyplot as plt
from hammock.lib import FastHyperLogLog
from hammock.lib import HyperLogLog

# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)

def generate_data(size, unique_ratio=1.0):
    """Generate test data with controlled uniqueness."""
    unique_count = int(size * unique_ratio)
    unique_items = [f"unique_item_{i}" for i in range(unique_count)]
    
    if unique_ratio < 1.0:
        # Generate data with duplicates
        result = []
        for _ in range(size):
            if random.random() < unique_ratio:
                # Add a unique item
                result.append(random.choice(unique_items))
            else:
                # Add a duplicate (from the first half of unique items)
                result.append(random.choice(unique_items[:unique_count//2]))
        return result
    else:
        # All unique items
        return [f"item_{i}" for i in range(size)]

def benchmark_add(data_sizes, sketch_classes):
    """Benchmark adding items to sketches."""
    results = {name: [] for name in sketch_classes.keys()}
    
    for size in data_sizes:
        data = generate_data(size)
        
        for name, cls in sketch_classes.items():
            # Create sketch
            sketch = cls(precision=12)
            
            # Measure add time
            start_time = time.time()
            
            if name == "FastHyperLogLog (Rust)" and hasattr(sketch, "add_batch"):
                sketch.add_batch(data)
            else:
                for item in data:
                    if name == "HyperLogLog (Python)":
                        sketch.add_string(item)
                    else:
                        sketch.add(item)
            
            elapsed = time.time() - start_time
            results[name].append(elapsed)
            
            print(f"{name}: Added {size} items in {elapsed:.4f}s")
    
    return results

def benchmark_cardinality(data_sizes, methods, sketch_classes):
    """Benchmark cardinality estimation with different methods."""
    results = {f"{name} ({method})": [] 
               for name in sketch_classes.keys() 
               for method in methods}
    
    for size in data_sizes:
        data = generate_data(size)
        
        for name, cls in sketch_classes.items():
            # Create and fill sketch
            sketch = cls(precision=12)
            
            if name == "FastHyperLogLog (Rust)" and hasattr(sketch, "add_batch"):
                sketch.add_batch(data)
            else:
                for item in data:
                    if name == "HyperLogLog (Python)":
                        sketch.add_string(item)
                    else:
                        sketch.add(item)
            
            # Benchmark each estimation method
            for method in methods:
                start_time = time.time()
                
                for _ in range(10):  # Run multiple times for more stable measurements
                    if name == "HyperLogLog (Python)":
                        sketch.estimate_cardinality(method=method)
                    else:
                        sketch.estimate_cardinality(method=method)
                
                elapsed = (time.time() - start_time) / 10
                results[f"{name} ({method})"].append(elapsed)
                
                print(f"{name} ({method}): Estimated {size} items in {elapsed:.6f}s")
    
    return results

def benchmark_merge(data_sizes, sketch_classes):
    """Benchmark merging sketches."""
    results = {name: [] for name in sketch_classes.keys()}
    
    for size in data_sizes:
        for name, cls in sketch_classes.items():
            # Create two sketches
            sketch1 = cls(precision=12)
            sketch2 = cls(precision=12)
            
            # Fill sketches with different data
            data1 = generate_data(size, unique_ratio=0.8)
            data2 = generate_data(size, unique_ratio=0.8)
            
            if name == "FastHyperLogLog (Rust)" and hasattr(sketch1, "add_batch"):
                sketch1.add_batch(data1)
                sketch2.add_batch(data2)
            else:
                for item in data1:
                    if name == "HyperLogLog (Python)":
                        sketch1.add_string(item)
                    else:
                        sketch1.add(item)
                
                for item in data2:
                    if name == "HyperLogLog (Python)":
                        sketch2.add_string(item)
                    else:
                        sketch2.add(item)
            
            # Measure merge time
            start_time = time.time()
            
            for _ in range(10):  # Run multiple times for stability
                if name == "HyperLogLog (Python)":
                    merged = HyperLogLog(precision=12)
                    merged.merge(sketch1)
                    merged.merge(sketch2)
                else:
                    merged = FastHyperLogLog(precision=12)
                    merged.merge(sketch1)
                    merged.merge(sketch2)
            
            elapsed = (time.time() - start_time) / 10
            results[name].append(elapsed)
            
            print(f"{name}: Merged {size} items in {elapsed:.6f}s")
    
    return results

def benchmark_accuracy(data_sizes, sketch_classes):
    """Benchmark estimation accuracy."""
    results = {name: {"error": [], "estimate": []} for name in sketch_classes.keys()}
    
    for size in data_sizes:
        data = generate_data(size)
        
        for name, cls in sketch_classes.items():
            # Create and fill sketch
            sketch = cls(precision=12)
            
            if name == "FastHyperLogLog (Rust)" and hasattr(sketch, "add_batch"):
                sketch.add_batch(data)
            else:
                for item in data:
                    if name == "HyperLogLog (Python)":
                        sketch.add_string(item)
                    else:
                        sketch.add(item)
            
            # Get estimate
            if name == "HyperLogLog (Python)":
                estimate = sketch.estimate_cardinality(method="ertl_mle")
            else:
                estimate = sketch.estimate_cardinality(method="ertl_mle")
            
            # Calculate error
            error = abs(estimate - size) / size * 100
            
            results[name]["error"].append(error)
            results[name]["estimate"].append(estimate)
            
            print(f"{name}: Size {size}, Estimate {estimate:.2f}, Error {error:.2f}%")
    
    return results

def plot_results(title, x_data, y_data, x_label, y_label, legend_loc='upper left'):
    """Plot benchmark results."""
    plt.figure(figsize=(10, 6))
    
    for name, data in y_data.items():
        plt.plot(x_data, data, marker='o', linewidth=2, label=name)
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc=legend_loc)
    plt.tight_layout()
    
    # Save to file
    filename = title.lower().replace(' ', '_') + '.png'
    plt.savefig(filename)
    print(f"Saved plot to {filename}")
    
    plt.close()

def run_benchmarks():
    """Run all benchmarks."""
    # Check if Rust implementation is available
    rust_available = FastHyperLogLog(precision=12).is_using_rust()
    
    if rust_available:
        sketch_classes = {
            "FastHyperLogLog (Rust)": FastHyperLogLog,
            "HyperLogLog (Python)": HyperLogLog
        }
    else:
        print("Rust implementation not available, using Python implementation only")
        sketch_classes = {
            "FastHyperLogLog (Python fallback)": FastHyperLogLog,
            "HyperLogLog (Python)": HyperLogLog
        }
    
    # Define data sizes for benchmarks
    small_sizes = [1000, 5000, 10000, 50000, 100000]
    large_sizes = [100000, 500000, 1000000, 5000000]
    
    # Run add benchmark
    print("\n=== Benchmarking Add Operation ===")
    add_results = benchmark_add(small_sizes, sketch_classes)
    plot_results(
        "HyperLogLog Add Performance", 
        small_sizes, 
        add_results, 
        "Number of Items", 
        "Time (seconds)"
    )
    
    # Run cardinality estimation benchmark
    print("\n=== Benchmarking Cardinality Estimation ===")
    methods = ["original", "ertl_improved", "ertl_mle"]
    cardinality_results = benchmark_cardinality(small_sizes, methods, sketch_classes)
    plot_results(
        "HyperLogLog Cardinality Estimation Performance", 
        small_sizes, 
        cardinality_results, 
        "Number of Items", 
        "Time (seconds)"
    )
    
    # Run merge benchmark
    print("\n=== Benchmarking Merge Operation ===")
    merge_results = benchmark_merge(small_sizes, sketch_classes)
    plot_results(
        "HyperLogLog Merge Performance", 
        small_sizes, 
        merge_results, 
        "Number of Items per Sketch", 
        "Time (seconds)"
    )
    
    # Run accuracy benchmark
    print("\n=== Benchmarking Estimation Accuracy ===")
    accuracy_results = benchmark_accuracy(small_sizes, sketch_classes)
    
    # Plot accuracy results
    plt.figure(figsize=(10, 6))
    
    for name, data in accuracy_results.items():
        plt.plot(small_sizes, data["error"], marker='o', linewidth=2, label=name)
    
    plt.title("HyperLogLog Estimation Error")
    plt.xlabel("Number of Items")
    plt.ylabel("Error (%)")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='upper right')
    plt.tight_layout()
    
    plt.savefig("hyperloglog_estimation_error.png")
    print("Saved plot to hyperloglog_estimation_error.png")
    
    # If Rust is available, run large-scale benchmark
    if rust_available:
        print("\n=== Running Large-Scale Benchmark (Rust only) ===")
        large_sketch_classes = {"FastHyperLogLog (Rust)": FastHyperLogLog}
        large_add_results = benchmark_add(large_sizes, large_sketch_classes)
        plot_results(
            "HyperLogLog Large-Scale Add Performance", 
            large_sizes, 
            large_add_results, 
            "Number of Items", 
            "Time (seconds)"
        )

if __name__ == "__main__":
    run_benchmarks() 