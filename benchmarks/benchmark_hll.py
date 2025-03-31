#!/usr/bin/env python3
import time
import random
import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import os
import sys
from hammock.lib.rusthll import RustHLL
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
PRECISION_VALUES = [6, 8, 10, 12, 14, 16]
NUM_ITEMS = [1000, 10000, 100000, 1000000, 10000000]
METHODS = ["original", "ertl_improved", "ertl_mle"]

# Create a debug log file
DEBUG_LOG = os.path.join(os.path.dirname(__file__), 'benchmark_debug.log')
# Redirect stdout to the debug log
sys.stdout = open(DEBUG_LOG, 'w')

print("\nNOTE: Before running benchmarks, it is recommended to run the tests first:")
print("      python -m pytest tests/test_rusthll.py::test_high_precision_support -v\n")

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

def benchmark_add(sketch_classes, precision=14):
    """Benchmark adding items to different HLL implementations."""
    results = {}
    
    for sketch_class in sketch_classes:
        sketch_name = sketch_class.__name__
        sketch = sketch_class(precision, hash_size=32)
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

def benchmark_add_batch(sketch_classes, precision=14):
    """Benchmark adding items in batch to different HLL implementations."""
    results = {}
    
    for sketch_class in sketch_classes:
        sketch_name = sketch_class.__name__
        sketch = sketch_class(precision, hash_size=32)
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

def benchmark_cardinality(sketch_classes, precision=14):
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
            sketch = sketch_class(precision, hash_size=32)
            
            # Check if using native or fallback implementation
            using_native = False
            if hasattr(sketch, 'is_using_rust'):
                using_native = sketch.is_using_rust()
            elif hasattr(sketch, 'using_rust'):
                using_native = sketch.using_rust
            implementation = "native" if using_native else "fallback"
            
            # Add items
            if hasattr(sketch, 'add_batch'):
                sketch.add_batch([str(item) for item in items])
            else:
                for item in items:
                    if sketch_name == "HyperLogLog":
                        sketch.add_string(str(item))
                    elif sketch_name == "RustHLLWrapper":
                        sketch.add_string(str(item))
                    else:  # FastHyperLogLog or RustHLL
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
                    except TypeError:
                        # Some implementations might not accept the method parameter
                        start_time = time.time()
                        estimates[method] = sketch.estimate_cardinality()
                        end_time = time.time()
                        times[method] = end_time - start_time
            
            results[sketch_name][num] = {
                'times': times,
                'estimates': estimates,
                'actual': num,
                'implementation': implementation
            }
    
    return results

def benchmark_merge(sketch_classes, precision=14):
    """Benchmark merging two sketches."""
    results = {}
    
    for sketch_class in sketch_classes:
        sketch_name = sketch_class.__name__
        results[sketch_name] = {}
        
        for num in NUM_ITEMS:
            # Create two sketches with some overlap
            sketch1 = sketch_class(precision, hash_size=32)
            sketch2 = sketch_class(precision, hash_size=32)
            
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
            sketch_copy = sketch_class(precision, hash_size=32)
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

def benchmark_accuracy(sketch_classes, precision_values=None):
    """Benchmark the accuracy of different HLL implementations across precisions."""
    if precision_values is None:
        precision_values = PRECISION_VALUES
        
    results = {}
    
    # Make sure results dictionary is properly initialized with all precisions
    for precision in precision_values:
        results[precision] = {}
        for sketch_class in sketch_classes:
            sketch_name = sketch_class.__name__
            results[precision][sketch_name] = {}
    
    # Run benchmarks
    for precision in precision_values:
        print(f"\nTesting precision {precision}...")
        
        for sketch_class in sketch_classes:
            sketch_name = sketch_class.__name__
            
            for num in NUM_ITEMS:
                print(f"  {sketch_name} with {num} items...", end="", flush=True)
                try:
                    items = list(range(num))
                    sketch = sketch_class(precision, hash_size=32)
                    
                    # Check if using native or fallback implementation
                    using_native = False
                    if hasattr(sketch, 'is_using_rust'):
                        using_native = sketch.is_using_rust()
                    elif hasattr(sketch, 'using_rust'):
                        using_native = sketch.using_rust
                    implementation = "native" if using_native else "fallback"
                    
                    # Add items
                    if hasattr(sketch, 'add_batch'):
                        sketch.add_batch([str(item) for item in items])
                    else:
                        for item in items:
                            if sketch_name == "HyperLogLog":
                                sketch.add_string(str(item))
                            elif sketch_name == "RustHLLWrapper":
                                sketch.add(item)
                            else:  # FastHyperLogLog or RustHLL
                                sketch.add_string(str(item))
                    
                    # Get estimate
                    if sketch_name == "RustHLLWrapper":
                        estimate = sketch.estimate_cardinality()
                    else:
                        # Try the MLE estimator for fairest comparison
                        try:
                            estimate = sketch.estimate_cardinality(method="ertl_mle")
                        except (ValueError, TypeError):
                            # Fallback if method parameter not supported
                            estimate = sketch.estimate_cardinality()
                    
                    # Calculate error
                    error = abs(estimate - num) / num * 100
                    
                    results[precision][sketch_name][num] = {
                        'estimate': float(estimate),
                        'actual': num,
                        'error': float(error),
                        'implementation': implementation
                    }
                    print(f" error: {error:.2f}%")
                except Exception as e:
                    # Record the error but keep going
                    print(f" ERROR: {e}")
                    results[precision][sketch_name][num] = {
                        'error': "Failed",
                        'exception': str(e)
                    }
    
    return results

def save_results_to_file(results, filename):
    """Save benchmark results to a file in the results directory."""
    # Use the filename as provided without adding another timestamp
    with open(filename, 'w') as f:
        for category, data in results.items():
            f.write(f"\n{category.upper()}:\n")
            f.write("-" * 40 + "\n")
            for key, value in data.items():
                f.write(f"{key}: {value}\n")
    print(f"Results saved to {filename}")

def plot_results(title, x_data, y_data, x_label, y_label, legend_loc='upper left'):
    """Plot benchmark results."""
    plt.figure(figsize=(10, 6))
    has_data = False
    for name, data in y_data.items():
        x = []
        y = []
        for num in x_data:
            if num in data:
                x.append(num)
                y.append(data[num])
        if x and y:  # Only plot if we have data
            plt.plot(x, y, marker='o', label=name)
            has_data = True
    
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.grid(True)
    if has_data:
        plt.legend(loc=legend_loc)
    plt.xscale('log')
    plt.yscale('log')
    return plt.gcf()

def plot_all_results(add_results, batch_results, cardinality_results, merge_results, accuracy_results, base_filename):
    """Create a single figure with all benchmark plots."""
    # Create a figure with subplots
    fig = plt.figure(figsize=(20, 15))
    
    # Add performance subplot
    ax1 = fig.add_subplot(231)
    has_data1 = False
    for name, data in add_results.items():
        x = []
        y = []
        for num in NUM_ITEMS:
            if num in data:
                x.append(num)
                y.append(data[num]['items_per_second'])
        if x and y:  # Only plot if we have data
            ax1.plot(x, y, marker='o', label=name)
            has_data1 = True
    ax1.set_title('Add Performance')
    ax1.set_xlabel('Number of Items')
    ax1.set_ylabel('Items per Second')
    ax1.grid(True)
    if has_data1:
        ax1.legend(loc='upper left')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    
    # Batch add performance subplot
    ax2 = fig.add_subplot(232)
    has_data2 = False
    for name, data in batch_results.items():
        x = []
        y = []
        for num in NUM_ITEMS:
            if num in data:
                x.append(num)
                y.append(data[num]['items_per_second'])
        if x and y:  # Only plot if we have data
            ax2.plot(x, y, marker='o', label=name)
            has_data2 = True
    ax2.set_title('Batch Add Performance')
    ax2.set_xlabel('Number of Items')
    ax2.set_ylabel('Items per Second')
    ax2.grid(True)
    if has_data2:
        ax2.legend(loc='upper left')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    
    # Cardinality estimation time subplot
    ax3 = fig.add_subplot(233)
    has_data3 = False
    for name, data in cardinality_results.items():
        x = []
        y = []
        for num in NUM_ITEMS:
            if num >= 100000 and num in data:
                if 'fast_mle' in data[num]['times']:
                    x.append(num)
                    y.append(data[num]['times']['fast_mle'])
                elif 'ertl_mle' in data[num]['times']:
                    x.append(num)
                    y.append(data[num]['times']['ertl_mle'])
        if x and y:  # Only plot if we have data
            ax3.plot(x, y, marker='o', label=name)
            has_data3 = True
    ax3.set_title('Cardinality Estimation Time')
    ax3.set_xlabel('Number of Items')
    ax3.set_ylabel('Time (seconds)')
    ax3.grid(True)
    if has_data3:
        ax3.legend(loc='upper left')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    
    # Merge performance subplot
    ax4 = fig.add_subplot(234)
    has_data4 = False
    for name, data in merge_results.items():
        x = []
        y = []
        for num in NUM_ITEMS:
            if num in data:
                x.append(num)
                y.append(data[num]['time'])
        if x and y:  # Only plot if we have data
            ax4.plot(x, y, marker='o', label=name)
            has_data4 = True
    ax4.set_title('Merge Performance')
    ax4.set_xlabel('Number of Items')
    ax4.set_ylabel('Time (seconds)')
    ax4.grid(True)
    if has_data4:
        ax4.legend(loc='upper left')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    
    # Accuracy subplot - use precision 14 as before
    ax5 = fig.add_subplot(235)
    has_data5 = False
    
    # Check if precision 14 is available
    if 14 in accuracy_results:
        print(f"Found precision 14 data, plotting accuracy...")
        for sketch_name in accuracy_results[14].keys():
            x = []
            y = []
            for num in NUM_ITEMS:
                if num in accuracy_results[14][sketch_name]:
                    error_value = accuracy_results[14][sketch_name][num].get('error')
                    if isinstance(error_value, (int, float)):
                        x.append(num)
                        y.append(error_value)
                        print(f"  {sketch_name} at {num} items: error = {error_value:.2f}%")
            
            if x and y:  # Only plot if we have data
                ax5.plot(x, y, marker='o', label=sketch_name)
                has_data5 = True
                print(f"  Plotted {len(x)} points for {sketch_name}")
            else:
                print(f"  No valid data to plot for {sketch_name}")
        
        ax5.set_title('Accuracy (precision=14)')
        ax5.set_xlabel('Number of Items')
        ax5.set_ylabel('Error (%)')
        ax5.grid(True)
        if has_data5:
            ax5.legend(loc='upper left')
        ax5.set_xscale('log')
    else:
        print(f"ERROR: Precision 14 not found in accuracy_results! Available: {list(accuracy_results.keys())}")
        ax5.set_title('Accuracy data for precision 14 not available')
        ax5.text(0.5, 0.5, "Precision 14 data not available", 
                horizontalalignment='center',
                verticalalignment='center',
                transform=ax5.transAxes)
        ax5.grid(True)
    
    # Adjust layout to prevent overlap
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(os.path.join(RESULTS_DIR, f'{base_filename}_results.png'))
    plt.close()

def plot_accuracy_results(results, base_filename):
    """Create a single figure with accuracy subplots for different item counts."""
    # First check if we have valid results
    print("\nPreparing accuracy plots...")
    
    if not results:
        print("ERROR: No accuracy results available!")
        return
        
    # Print available precision keys
    print(f"Available precision keys: {list(results.keys())}")
    
    # Create a figure with subplots
    fig = plt.figure(figsize=(20, 15))
    
    # Get the number of item counts and create a grid of subplots
    num_items = len(NUM_ITEMS)
    num_cols = 2  # Number of columns in the grid
    num_rows = (num_items + num_cols - 1) // num_cols  # Calculate number of rows needed
    
    for idx, num in enumerate(NUM_ITEMS):
        ax = fig.add_subplot(num_rows, num_cols, idx + 1)
        
        # Track if we plotted any data for this subplot
        has_data = False
        
        # Plot accuracy for each sketch type
        for sketch_name in results[list(results.keys())[0]].keys():
            print(f"Processing {sketch_name} for item count {num}...")
            
            precisions = []
            errors = []
            
            for p in results.keys():
                # Check if we have data for this precision, sketch, and item count
                if sketch_name in results[p] and num in results[p][sketch_name]:
                    error_value = results[p][sketch_name][num].get('error')
                    
                    # Only use if it's a numeric error
                    if isinstance(error_value, (int, float)):
                        precisions.append(p)
                        errors.append(error_value)
                        print(f"  Precision {p}: error = {error_value:.2f}%")
                    else:
                        print(f"  Precision {p}: INVALID error value: {error_value}")
            
            if precisions and errors:  # Only plot if we have data
                ax.plot(precisions, errors, marker='o', label=sketch_name)
                has_data = True
                print(f"  Plotted {len(precisions)} points for {sketch_name}")
            else:
                print(f"  No valid data to plot for {sketch_name}")
        
        ax.set_title(f'Accuracy at {num:,} Items')
        ax.set_xlabel('Precision')
        ax.set_ylabel('Error (%)')
        ax.grid(True)
        
        # Only add a legend if we actually plotted data
        if has_data:
            ax.legend(loc='upper right')
        else:
            ax.text(0.5, 0.5, "No data available", 
                   horizontalalignment='center',
                   verticalalignment='center',
                   transform=ax.transAxes)
    
    # Adjust layout to prevent overlap
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(os.path.join(RESULTS_DIR, f'{base_filename}_accuracy.png'))
    plt.close()

def run_benchmarks():
    """Run all benchmarks and save results."""
    print("\n========================================================")
    print("    HyperLogLog Benchmark")
    print("========================================================")
    print(f"Precision values: {PRECISION_VALUES}")
    print(f"Data sizes: {NUM_ITEMS}")
    
    # Define sketch classes to benchmark
    sketch_classes = [HyperLogLog]
    if RUST_HLL_AVAILABLE:
        sketch_classes.extend([RustHLL, RustHLLWrapper])
        print(f"Benchmarking {len(sketch_classes)} sketch implementations")
    else:
        print("Benchmarking Python HyperLogLog implementation only")
    
    print("Running benchmarks...")
    
    # Run benchmarks
    add_results = benchmark_add(sketch_classes)
    batch_results = benchmark_add_batch(sketch_classes)
    cardinality_results = benchmark_cardinality(sketch_classes)
    merge_results = benchmark_merge(sketch_classes)
    
    # Debug print for benchmark functions
    print("\nResults structure check:")
    print(f"add_results keys: {list(add_results.keys())}")
    if add_results:
        first_key = list(add_results.keys())[0]
        print(f"Sample add_results[{first_key}] keys: {list(add_results[first_key].keys())}")
    
    print(f"batch_results keys: {list(batch_results.keys())}")
    print(f"cardinality_results keys: {list(cardinality_results.keys())}")
    print(f"merge_results keys: {list(merge_results.keys())}")
    
    # Run accuracy benchmark
    print("\nRunning accuracy benchmark...")
    accuracy_results = benchmark_accuracy(sketch_classes)
    
    # Debug print for accuracy results
    print("\nAccuracy results structure:")
    print(f"accuracy_results keys: {list(accuracy_results.keys())}")
    if 14 in accuracy_results:
        print(f"Precision 14 keys: {list(accuracy_results[14].keys())}")
        for name in accuracy_results[14]:
            print(f"  {name} keys: {list(accuracy_results[14][name].keys())}")
    else:
        print(f"ERROR: Precision 14 not in accuracy_results!")
        print(f"Available precisions: {list(accuracy_results.keys())}")
    
    print("Benchmarks completed")
    
    # Generate timestamp once for both outputs
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    script_name = os.path.splitext(os.path.basename(__file__))[0]
    base_filename = f'{script_name}_{timestamp}'
    
    # Save results to file
    save_results_to_file({
        'ADD': add_results,
        'ADD_BATCH': batch_results,
        'CARDINALITY': cardinality_results,
        'MERGE': merge_results,
        'ACCURACY': accuracy_results
    }, os.path.join(RESULTS_DIR, f'{base_filename}.txt'))
    
    # Plot all results in a single figure
    plot_all_results(add_results, batch_results, cardinality_results, merge_results, accuracy_results, base_filename)
    
    # Plot accuracy results in a separate figure with subplots
    plot_accuracy_results(accuracy_results, base_filename)

if __name__ == "__main__":
    run_benchmarks() 