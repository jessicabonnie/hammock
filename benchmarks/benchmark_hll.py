#!/usr/bin/env python3
import time
import random
import numpy as np
import matplotlib.pyplot as plt
import os
from hammock.lib.rusthll import RustHyperLogLog
from hammock.lib.hyperloglog import HyperLogLog
from datetime import datetime

# Create results directory if it doesn't exist
RESULTS_DIR = os.path.join(os.path.dirname(__file__), 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)

# Try importing the rust_hll module
try:
    from hammock.lib.rusthll_compat import RustHLLWrapper, RUST_AVAILABLE
    RUST_HLL_AVAILABLE = True
except ImportError:
    RUST_HLL_AVAILABLE = False

# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)

# Constants
PRECISION_VALUES = [8, 10, 12, 14, 16]
NUM_ITEMS = [1000, 10000, 100000, 1000000]
METHODS = ["original", "ertl_improved", "ertl_mle"]

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

def benchmark_add(sketch_classes, precision=12):
    """Benchmark adding items to different HLL implementations."""
    results = {}
    
    for sketch_class in sketch_classes:
        sketch_name = sketch_class.__name__
        sketch = sketch_class(precision)
        results[sketch_name] = {}
        
        for num in NUM_ITEMS:
            items = list(range(num))
            start_time = time.time()
            for item in items:
                if sketch_name == "HyperLogLog":
                    sketch.add_string(str(item))
                elif sketch_name == "RustHLLWrapper":
                    sketch.add(item)
                else:  # FastHyperLogLog
                    sketch.add_string(str(item))
            end_time = time.time()
            
            time_taken = end_time - start_time
            items_per_second = num / time_taken if time_taken > 0 else float('inf')
            
            results[sketch_name][num] = {
                'time': time_taken,
                'items_per_second': items_per_second
            }
    
    return results

def benchmark_add_batch(sketch_classes, precision=12):
    """Benchmark adding items in batch to different HLL implementations."""
    results = {}
    
    for sketch_class in sketch_classes:
        sketch_name = sketch_class.__name__
        sketch = sketch_class(precision)
        results[sketch_name] = {}
        
        for num in NUM_ITEMS:
            items = [str(item) for item in range(num)]  # Convert integers to strings
            
            if hasattr(sketch, 'add_batch'):
                start_time = time.time()
                sketch.add_batch(items)
                end_time = time.time()
                
                time_taken = end_time - start_time
                items_per_second = num / time_taken if time_taken > 0 else float('inf')
            else:
                # Fallback for implementations without batch support
                time_taken = float('inf')
                items_per_second = 0
            
            results[sketch_name][num] = {
                'time': time_taken,
                'items_per_second': items_per_second
            }
    
    return results

def benchmark_cardinality(sketch_classes, precision=12):
    """Benchmark cardinality estimation using different methods."""
    results = {}
    
    for sketch_class in sketch_classes:
        sketch_name = sketch_class.__name__
        results[sketch_name] = {}
        
        for num in NUM_ITEMS:
            # Only test one of the larger datasets
            if num < 100000:
                continue
                
            items = list(range(num))
            sketch = sketch_class(precision)
            
            # Add items
            if hasattr(sketch, 'add_batch'):
                sketch.add_batch([str(item) for item in items])
            else:
                for item in items:
                    if sketch_name == "HyperLogLog":
                        sketch.add_string(str(item))
                    elif sketch_name == "RustHLLWrapper":
                        sketch.add_string(str(item))
                    else:  # FastHyperLogLog
                        sketch.add_string(str(item))
            
            times = {}
            estimates = {}
            
            # For Rust implementation, we only have one estimation method
            if sketch_name == "RustHLLWrapper":
                start_time = time.time()
                estimates["fast_mle"] = sketch.estimate_cardinality()
                end_time = time.time()
                times["fast_mle"] = end_time - start_time
            else:
                # For other implementations, test all methods
                for method in METHODS:
                    try:
                        start_time = time.time()
                        estimates[method] = sketch.estimate_cardinality(method=method)
                        end_time = time.time()
                        times[method] = end_time - start_time
                    except ValueError as e:
                        # Skip methods not supported by this implementation
                        print(f"Skipping method {method} for {sketch_name}: {e}")
                        continue
            
            results[sketch_name][num] = {
                'times': times,
                'estimates': estimates,
                'actual': num
            }
    
    return results

def benchmark_merge(sketch_classes, precision=12):
    """Benchmark merging two sketches."""
    results = {}
    
    for sketch_class in sketch_classes:
        sketch_name = sketch_class.__name__
        results[sketch_name] = {}
        
        for num in NUM_ITEMS:
            # Create two sketches with some overlap
            sketch1 = sketch_class(precision)
            sketch2 = sketch_class(precision)
            
            overlap = num // 2
            items1 = list(range(num))
            items2 = list(range(overlap, overlap + num))
            
            # Add items to sketches
            if hasattr(sketch1, 'add_batch'):
                sketch1.add_batch([str(item) for item in items1])
                sketch2.add_batch([str(item) for item in items2])
            else:
                for item in items1:
                    if sketch_name == "HyperLogLog":
                        sketch1.add_string(str(item))
                    elif sketch_name == "RustHLLWrapper":
                        sketch1.add(item)
                    else:  # FastHyperLogLog
                        sketch1.add_string(str(item))
                for item in items2:
                    if sketch_name == "HyperLogLog":
                        sketch2.add_string(str(item))
                    elif sketch_name == "RustHLLWrapper":
                        sketch2.add(item)
                    else:  # FastHyperLogLog
                        sketch2.add_string(str(item))
            
            # Measure merge time
            sketch_copy = sketch_class(precision)
            if hasattr(sketch1, 'add_batch'):
                sketch_copy.add_batch([str(item) for item in items1])
            else:
                for item in items1:
                    if sketch_name == "HyperLogLog":
                        sketch_copy.add_string(str(item))
                    elif sketch_name == "RustHLLWrapper":
                        sketch_copy.add(item)
                    else:  # FastHyperLogLog
                        sketch_copy.add_string(str(item))
                    
            start_time = time.time()
            sketch_copy.merge(sketch2)
            end_time = time.time()
            
            time_taken = end_time - start_time
            est_cardinality = sketch_copy.estimate_cardinality() if hasattr(sketch_copy, 'estimate_cardinality') else sketch_copy.cardinality()
            actual_cardinality = len(set(items1 + items2))
            
            results[sketch_name][num] = {
                'time': time_taken,
                'estimated_cardinality': est_cardinality,
                'actual_cardinality': actual_cardinality,
                'error': abs(est_cardinality - actual_cardinality) / actual_cardinality * 100
            }
    
    return results

def benchmark_accuracy(sketch_classes, precision_values=PRECISION_VALUES):
    """Benchmark the accuracy of different HLL implementations across precisions."""
    results = {}
    
    for precision in precision_values:
        results[precision] = {}
        
        for sketch_class in sketch_classes:
            sketch_name = sketch_class.__name__
            results[precision][sketch_name] = {}
            
            for num in NUM_ITEMS:
                items = list(range(num))
                sketch = sketch_class(precision)
                
                # Add items
                if hasattr(sketch, 'add_batch'):
                    sketch.add_batch([str(item) for item in items])
                else:
                    for item in items:
                        if sketch_name == "HyperLogLog":
                            sketch.add_string(str(item))
                        elif sketch_name == "RustHLLWrapper":
                            sketch.add(item)
                        else:  # FastHyperLogLog
                            sketch.add_string(str(item))
                
                # Get estimate
                if sketch_name == "RustHLLWrapper":
                    estimate = sketch.estimate_cardinality()
                else:
                    # Use the MLE estimator for fairest comparison
                    estimate = sketch.estimate_cardinality(method="ertl_mle")
                
                # Calculate error
                error = abs(estimate - num) / num * 100
                
                results[precision][sketch_name][num] = {
                    'estimate': estimate,
                    'actual': num,
                    'error': error
                }
    
    return results

def save_results_to_file(results, filename):
    """Save benchmark results to a file in the results directory."""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    base, ext = os.path.splitext(filename)
    filepath = os.path.join(RESULTS_DIR, f"{base}_{timestamp}{ext}")
    with open(filepath, 'w') as f:
        for category, data in results.items():
            f.write(f"\n{category.upper()}:\n")
            f.write("-" * 40 + "\n")
            for key, value in data.items():
                f.write(f"{key}: {value}\n")
    print(f"Results saved to {filepath}")

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
    
    # Save to file in results directory with timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    filename = title.lower().replace(' ', '_') + f'_{timestamp}.png'
    filepath = os.path.join(RESULTS_DIR, filename)
    plt.savefig(filepath)
    print(f"Saved plot to {filepath}")
    
    plt.close()

def run_benchmarks():
    """Run all benchmarks and return the results."""
    sketch_classes = [HyperLogLog, RustHyperLogLog]
    
    if RUST_HLL_AVAILABLE:
        sketch_classes.append(RustHLLWrapper)
    
    results = {
        'add': benchmark_add(sketch_classes),
        'add_batch': benchmark_add_batch(sketch_classes),
        'cardinality': benchmark_cardinality(sketch_classes),
        'merge': benchmark_merge(sketch_classes),
        'accuracy': benchmark_accuracy(sketch_classes)
    }
    
    # Save results to file
    save_results_to_file(results, 'benchmark_results.txt')
    
    # Plot results
    # Plot add performance
    add_data = {name: [data[num]['items_per_second'] for num in NUM_ITEMS] 
                for name, data in results['add'].items()}
    plot_results('Add Performance', NUM_ITEMS, add_data, 'Number of Items', 'Items per Second')
    
    # Plot batch add performance
    batch_data = {name: [data[num]['items_per_second'] for num in NUM_ITEMS] 
                  for name, data in results['add_batch'].items()}
    plot_results('Batch Add Performance', NUM_ITEMS, batch_data, 'Number of Items', 'Items per Second')
    
    # Plot cardinality estimation times
    cardinality_data = {}
    for name, data in results['cardinality'].items():
        times = []
        for num in NUM_ITEMS:
            if num >= 100000:
                # For Rust implementation, we only have one estimation method
                if name == "RustHLLWrapper":
                    times.append(data[num]['times']['fast_mle'])
                else:
                    # For other implementations, use the fast_mle method if available
                    if 'fast_mle' in data[num]['times']:
                        times.append(data[num]['times']['fast_mle'])
                    else:
                        # Fallback to ertl_mle if fast_mle not available
                        times.append(data[num]['times']['ertl_mle'])
        cardinality_data[name] = times
    
    plot_results('Cardinality Estimation Time', [num for num in NUM_ITEMS if num >= 100000], 
                cardinality_data, 'Number of Items', 'Time (seconds)')
    
    # Plot merge performance
    merge_data = {name: [data[num]['time'] for num in NUM_ITEMS] 
                  for name, data in results['merge'].items()}
    plot_results('Merge Performance', NUM_ITEMS, merge_data, 'Number of Items', 'Time (seconds)')
    
    # Plot accuracy across precisions
    for num in NUM_ITEMS:
        accuracy_data = {cls.__name__: [results['accuracy'][p][cls.__name__][num]['error'] for p in PRECISION_VALUES] 
                        for cls in sketch_classes}
        plot_results(f'Accuracy at {num} Items', PRECISION_VALUES, accuracy_data, 
                    'Precision', 'Error (%)', legend_loc='upper right')
    
    return results

if __name__ == "__main__":
    run_benchmarks() 