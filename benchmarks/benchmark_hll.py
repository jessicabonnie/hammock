#!/usr/bin/env python3
import time
import random
import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import os
import sys
from hammock.lib.rusthll_compat import RustHLLWrapper
from hammock.lib.hyperloglog import HyperLogLog
from datetime import datetime

# Create results directory if it doesn't exist
RESULTS_DIR = os.path.join(os.path.dirname(__file__), 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)

# Try importing the rust_hll module
try:
    import rust_hll
    from hammock.lib.rusthll_compat import RUST_AVAILABLE, RustHLLWrapper
    RUST_HLL_AVAILABLE = True
except ImportError:
    RUST_HLL_AVAILABLE = False

# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)

# Constants
PRECISION_VALUES = [8, 10, 12, 14, 16, 18, 20]
NUM_ITEMS = [1000, 10000, 100000, 1000000, 10000000]
METHODS = ["original", "ertl_mle", "fast_mle"]

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
        sketch = sketch_class(precision=precision, hash_size=32, debug=False)
        results[sketch_name] = {}
        
        for num in NUM_ITEMS:
            items = list(range(num))
            start_time = time.time()
            for item in items:
                if sketch_name == "HyperLogLog":
                    sketch.add_string(str(item))
                elif sketch_name == "RustHLLWrapper":
                    sketch.add(str(item))
                elif sketch_name == "RustHLL":
                    sketch.add_value(str(item))
                else:  # Fallback
                    sketch.add(str(item))
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
        sketch = sketch_class(precision, hash_size=32, debug=False)
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

def benchmark_cardinality(sketch_classes, precision=14, num_runs=5):
    """Benchmark cardinality estimation using different methods."""
    results = {}
    
    for sketch_class in sketch_classes:
        sketch_name = sketch_class.__name__
        results[sketch_name] = {}
        
        for num in NUM_ITEMS:
            # Only test one of the larger datasets
            if num < 100000:
                continue
                
            # Initialize results for this size
            results[sketch_name][num] = {
                'times': {},
                'estimates': {},
                'actual': num,
                'implementation': None,
                'errors': {}
            }
            
            # Run multiple times with different seeds
            for run in range(num_runs):
                seed = run + 1  # Use different seeds for each run
                items = list(range(num))
                sketch = sketch_class(precision=precision, hash_size=32, seed=seed, debug=False)
                
                # Check if using native or fallback implementation
                using_native = False
                if hasattr(sketch, 'is_using_rust'):
                    using_native = sketch.is_using_rust()
                elif hasattr(sketch, 'using_rust'):
                    using_native = sketch.using_rust
                elif sketch_name == "RustHLL":
                    using_native = True
                implementation = "native" if using_native else "fallback"
                
                # Store implementation type (should be same for all runs)
                if results[sketch_name][num]['implementation'] is None:
                    results[sketch_name][num]['implementation'] = implementation
                
                # Add items
                if hasattr(sketch, 'add_batch'):
                    sketch.add_batch([str(item) for item in items])
                else:
                    for item in items:
                        if sketch_name == "HyperLogLog":
                            sketch.add_string(str(item))
                        elif sketch_name == "RustHLLWrapper":
                            sketch.add_string(str(item))
                        elif sketch_name == "RustHLL":
                            sketch.add_value(str(item))
                        else:  # Fallback
                            sketch.add_string(str(item))
                
                # For Rust implementations, we only have one estimation method
                if sketch_name == "RustHLL":
                    start_time = time.time()
                    estimate = sketch.estimate_cardinality()
                    end_time = time.time()
                    time_taken = end_time - start_time
                    
                    # Store results under 'fast_mle' method for consistency
                    if 'fast_mle' not in results[sketch_name][num]['times']:
                        results[sketch_name][num]['times']['fast_mle'] = []
                        results[sketch_name][num]['estimates']['fast_mle'] = []
                        results[sketch_name][num]['errors']['fast_mle'] = []
                    
                    results[sketch_name][num]['times']['fast_mle'].append(time_taken)
                    results[sketch_name][num]['estimates']['fast_mle'].append(estimate)
                    results[sketch_name][num]['errors']['fast_mle'].append(abs(estimate - num) / num)
                
                # For RustHLLWrapper, store in fast_mle slot
                elif sketch_name == "RustHLLWrapper":
                    start_time = time.time()
                    estimate = sketch.estimate_cardinality()
                    end_time = time.time()
                    time_taken = end_time - start_time
                    
                    if 'fast_mle' not in results[sketch_name][num]['times']:
                        results[sketch_name][num]['times']['fast_mle'] = []
                        results[sketch_name][num]['estimates']['fast_mle'] = []
                        results[sketch_name][num]['errors']['fast_mle'] = []
                    
                    results[sketch_name][num]['times']['fast_mle'].append(time_taken)
                    results[sketch_name][num]['estimates']['fast_mle'].append(estimate)
                    results[sketch_name][num]['errors']['fast_mle'].append(abs(estimate - num) / num)
                
                # For Python implementation, test all methods
                else:
                    for method in METHODS:
                        try:
                            start_time = time.time()
                            estimate = sketch.estimate_cardinality(method=method)
                            end_time = time.time()
                            time_taken = end_time - start_time
                            
                            if method not in results[sketch_name][num]['times']:
                                results[sketch_name][num]['times'][method] = []
                                results[sketch_name][num]['estimates'][method] = []
                                results[sketch_name][num]['errors'][method] = []
                            
                            results[sketch_name][num]['times'][method].append(time_taken)
                            results[sketch_name][num]['estimates'][method].append(estimate)
                            results[sketch_name][num]['errors'][method].append(abs(estimate - num) / num)
                        except (ValueError, TypeError) as e:
                            # Skip methods not supported by this implementation
                            print(f"Skipping method {method} for {sketch_name}: {e}")
                            continue
    
    return results

def benchmark_merge(sketch_classes, precision=14):
    """Benchmark merging two sketches."""
    results = {}
    
    for sketch_class in sketch_classes:
        sketch_name = sketch_class.__name__
        results[sketch_name] = {}
        
        for num in NUM_ITEMS:
            # Create two sketches with some overlap
            sketch1 = sketch_class(precision, hash_size=32, debug=False)
            sketch2 = sketch_class(precision, hash_size=32, debug=False)
            
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
                        sketch1.add(str(item))
                    else:  # RustHyperLogLog
                        sketch1.add(str(item))
                for item in items2:
                    if sketch_name == "HyperLogLog":
                        sketch2.add_string(str(item))
                    elif sketch_name == "RustHLLWrapper":
                        sketch2.add(str(item))
                    else:  # RustHyperLogLog
                        sketch2.add(str(item))
            
            # Measure merge time
            start_time = time.time()
            sketch1.merge(sketch2)  # Merge directly into sketch1
            end_time = time.time()
            
            time_taken = end_time - start_time
            est_cardinality = sketch1.estimate_cardinality() if hasattr(sketch1, 'estimate_cardinality') else sketch1.cardinality()
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
                    sketch = sketch_class(precision, hash_size=32, debug=False)
                    
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

def plot_all_results(add_results, batch_results, cardinality_results, merge_results, accuracy_results, run_dir):
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
        for num in NUM_ITEMS:
            if num >= 100000 and num in data:
                # Plot each estimation method separately
                for method in data[num]['times'].keys():
                    x = []
                    y = []
                    x.append(num)
                    if isinstance(data[num]['times'][method], (list, np.ndarray)):
                        y.append(np.mean(data[num]['times'][method]))
                    else:
                        y.append(data[num]['times'][method])
                    if x and y:  # Only plot if we have data
                        ax3.plot(x, y, marker='o', label=f"{name} ({method})")
                        has_data3 = True
    ax3.set_title('Cardinality Estimation Time')
    ax3.set_xlabel('Number of Items')
    ax3.set_ylabel('Time (seconds)')
    ax3.grid(True)
    if has_data3:
        ax3.legend(loc='upper left', fontsize='small')
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
    plt.savefig(os.path.join(run_dir, 'results.png'))
    plt.close()

def plot_accuracy_results(results, run_dir):
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
    plt.savefig(os.path.join(run_dir, 'accuracy.png'))
    plt.close()

def plot_estimation_methods(cardinality_results, run_dir):
    """Plot the performance of different estimation methods for cardinality estimation."""
    print("\nPreparing estimation methods plots...")
    
    if not cardinality_results:
        print("ERROR: No cardinality results available!")
        return
    
    # Create a new figure with two subplots (one for accuracy, one for speed)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # For each sketch class and dataset size, plot the error percentage and time taken for each method
    sketch_names = list(cardinality_results.keys())
    dataset_sizes = sorted([size for size in NUM_ITEMS if size >= 100000])
    
    # Dictionary to hold data for plotting
    method_data = {
        'errors': {},
        'times': {}
    }
    
    # Parse the results
    for sketch_name in sketch_names:
        for size in dataset_sizes:
            if size not in cardinality_results[sketch_name]:
                continue
                
            result = cardinality_results[sketch_name][size]
            actual = result['actual']
            
            # Calculate mean and std of errors across runs
            for method in result['estimates'].keys():
                # Initialize method in dictionaries if not present
                if method not in method_data['errors']:
                    method_data['errors'][method] = []
                    method_data['times'][method] = []
                
                # Calculate mean error across runs
                errors = result['errors'][method]
                if isinstance(errors, (list, np.ndarray)):
                    mean_error = np.mean(errors)
                    std_error = np.std(errors)
                else:
                    # If we have a scalar error value
                    mean_error = errors
                    std_error = 0
                
                method_data['errors'][method].append((sketch_name, size, mean_error, std_error))
                
                # Calculate mean time across runs
                times = result['times'][method]
                if isinstance(times, (list, np.ndarray)):
                    mean_time = np.mean(times) * 1000  # Convert to milliseconds
                else:
                    mean_time = times * 1000  # Convert scalar time to milliseconds
                method_data['times'][method].append((sketch_name, size, mean_time))
    
    # Plot error rates with error bars
    for method in method_data['errors']:
        data = method_data['errors'][method]
        sizes = [d[1] for d in data]
        errors = [d[2] for d in data]
        stds = [d[3] for d in data]
        
        ax1.errorbar(sizes, errors, yerr=stds, label=f"{method}", marker='o')
    
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Dataset Size')
    ax1.set_ylabel('Mean Error (%)')
    ax1.set_title('Cardinality Estimation Error by Method')
    ax1.grid(True)
    ax1.legend()
    
    # Plot execution times
    for method in method_data['times']:
        data = method_data['times'][method]
        sizes = [d[1] for d in data]
        times = [d[2] for d in data]
        
        ax2.plot(sizes, times, label=f"{method}", marker='o')
    
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Dataset Size')
    ax2.set_ylabel('Mean Time (ms)')
    ax2.set_title('Cardinality Estimation Time by Method')
    ax2.grid(True)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(run_dir, 'cardinality_estimates.png'))
    plt.close()

def plot_estimation_bias(cardinality_results, run_dir):
    """Plot bias of different estimation methods by analyzing their estimated-to-actual ratio."""
    print("\nPreparing estimation bias plots...")
    
    if not cardinality_results:
        print("ERROR: No cardinality results available!")
        return
    
    # Create a figure
    fig, ax = plt.subplots(figsize=(15, 10))
    
    # For each sketch class and dataset size, calculate the bias ratio (estimate/actual)
    sketch_names = list(cardinality_results.keys())
    dataset_sizes = sorted([size for size in NUM_ITEMS if size >= 100000])
    
    # Dictionary to hold data for plotting
    method_data = {}
    
    # Parse the results
    for sketch_name in sketch_names:
        for size in dataset_sizes:
            if size not in cardinality_results[sketch_name]:
                continue
                
            result = cardinality_results[sketch_name][size]
            actual = result['actual']
            
            # For each method that has data
            for method in result['estimates'].keys():
                # Initialize method in dictionary if not present
                if method not in method_data:
                    method_data[method] = []
                
                # Get the estimate (might be a list for multiple runs)
                estimates = result['estimates'][method]
                if isinstance(estimates, (list, np.ndarray)) and len(estimates) > 0:
                    estimate = np.mean(estimates)
                    ratio = float(estimate) / actual
                    method_data[method].append((sketch_name, size, ratio))
    
    # Check if we have any data to plot
    if not method_data:
        print("No valid data available for estimation bias plot")
        plt.close()
        return
    
    # Plot bias ratio (logarithmic, so 0 means unbiased, > 0 means overestimation, < 0 means underestimation)
    markers = ['o', 's', '^', 'v', 'D', '*']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    # Plot a horizontal line at y=1 (unbiased)
    unbiased_line = ax.axhline(y=1.0, color='k', linestyle='--', alpha=0.3, label="Unbiased (ratio=1)")
    
    has_data = False  # Track if we've plotted any data
    all_ratios = []  # Collect all ratios for y-axis limits
    plotted_lines = [unbiased_line]  # Keep track of plotted lines for legend
    
    # Group by method
    for i, (method, data_points) in enumerate(method_data.items()):
        # Group by sketch name
        sketch_dict = {}
        for sketch_name, size, ratio in data_points:
            if sketch_name not in sketch_dict:
                sketch_dict[sketch_name] = {'sizes': [], 'ratios': []}
            sketch_dict[sketch_name]['sizes'].append(size)
            sketch_dict[sketch_name]['ratios'].append(ratio)
            all_ratios.append(ratio)
        
        # Plot each sketch as a separate line
        for j, (sketch_name, data) in enumerate(sketch_dict.items()):
            if data['sizes'] and data['ratios']:
                marker_idx = j % len(markers)
                color_idx = i % len(colors)
                label = f"{sketch_name} ({method})"
                
                # Sort by size to ensure correct line order
                points = sorted(zip(data['sizes'], data['ratios']))
                sizes, ratios = zip(*points)
                
                line = ax.plot(sizes, ratios, marker=markers[marker_idx], color=colors[color_idx], 
                             label=label, linestyle='-')[0]
                plotted_lines.append(line)
                has_data = True
    
    if not has_data:
        print("No valid data to plot in estimation bias plot")
        plt.close()
        return
    
    ax.set_title('Estimation Method Bias (Estimate/Actual Ratio)')
    ax.set_xlabel('Dataset Size')
    ax.set_ylabel('Estimate/Actual Ratio')
    ax.set_xscale('log')
    ax.grid(True)
    
    # Add region highlighting the 0.95-1.05 range (5% error margin)
    error_band = ax.axhspan(0.95, 1.05, alpha=0.2, color='green', label="±5% error band")
    plotted_lines.append(error_band)
    
    # Calculate y-axis limits from actual data
    if all_ratios:
        y_min = min(all_ratios)
        y_max = max(all_ratios)
        padding = (y_max - y_min) * 0.1
        ax.set_ylim([max(0, y_min - padding), y_max + padding])
    
    # Only add legend if we have plotted lines
    if plotted_lines:
        ax.legend(handles=plotted_lines, loc='upper right', fontsize='small')
    
    plt.tight_layout()
    plt.savefig(os.path.join(run_dir, 'estimation_bias.png'))
    plt.close()
    
    # Create a boxplot to show the distribution of bias across all dataset sizes
    if not all_ratios:
        print("No data available for boxplot")
        return
        
    fig, ax = plt.subplots(figsize=(15, 8))
    
    # Prepare data for boxplot
    boxplot_data = []
    boxplot_labels = []
    
    # Combine data across all dataset sizes for each method/sketch combination
    for method, data_points in method_data.items():
        sketch_dict = {}
        for sketch_name, _, ratio in data_points:
            if f"{sketch_name}_{method}" not in sketch_dict:
                sketch_dict[f"{sketch_name}_{method}"] = []
            sketch_dict[f"{sketch_name}_{method}"].append(ratio)
        
        # Add each combination as a separate box
        for label, ratios in sketch_dict.items():
            if ratios:
                boxplot_data.append(ratios)
                boxplot_labels.append(label)
    
    if not boxplot_data:
        print("No data available for boxplot")
        plt.close()
        return
    
    # Create the boxplot
    bp = ax.boxplot(boxplot_data, patch_artist=True, showfliers=True)
    
    # Color the boxes according to the method
    method_names = list(method_data.keys())
    for i, box in enumerate(bp['boxes']):
        # Extract the method name from the label
        for method in method_names:
            if method in boxplot_labels[i]:
                method_idx = method_names.index(method)
                box.set_facecolor(colors[method_idx % len(colors)])
                break
    
    # Add a horizontal line at y=1 (unbiased)
    ax.axhline(y=1.0, color='k', linestyle='--', alpha=0.7)
    
    # Add region highlighting the 0.95-1.05 range (5% error margin)
    ax.axhspan(0.95, 1.05, alpha=0.2, color='green', label="±5% error band")
    
    ax.set_title('Distribution of Estimation Bias by Method and Implementation')
    ax.set_ylabel('Estimate/Actual Ratio')
    ax.set_xticklabels(boxplot_labels, rotation=45, ha='right')
    ax.grid(True, axis='y')
    
    plt.tight_layout()
    plt.savefig(os.path.join(run_dir, 'estimation_bias_boxplot.png'))
    plt.close()

def plot_performance_vs_accuracy(cardinality_results, run_dir):
    """Create a scatter plot showing the trade-off between performance and accuracy for different methods."""
    print("\nPreparing performance vs accuracy tradeoff plots...")
    
    if not cardinality_results:
        print("ERROR: No cardinality results available!")
        return
    
    # Create the figure
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # For each sketch class, dataset size, and method, calculate performance and accuracy metrics
    sketch_names = list(cardinality_results.keys())
    dataset_sizes = sorted([size for size in NUM_ITEMS if size >= 100000])
    
    # Collect data for plotting
    plot_data = []
    
    # Parse the results
    for sketch_name in sketch_names:
        for size in dataset_sizes:
            if size not in cardinality_results[sketch_name]:
                continue
                
            result = cardinality_results[sketch_name][size]
            actual = result['actual']
            
            # For each method that has data
            for method in result['estimates'].keys():
                # Get the estimates and times (might be lists for multiple runs)
                estimates = result['estimates'][method]
                times = result['times'][method]
                
                # Skip if not valid data
                if not isinstance(estimates, (list, np.ndarray)) or not isinstance(times, (list, np.ndarray)):
                    continue
                
                # Calculate mean values
                mean_estimate = np.mean(estimates)
                mean_time = np.mean(times)
                
                # Calculate error percentage
                error_pct = abs(float(mean_estimate) - actual) / actual * 100
                time_ms = mean_time * 1000  # Convert to milliseconds
                ops_per_second = 1 / mean_time if mean_time > 0 else float('inf')
                
                # Add to plot data
                plot_data.append({
                    'sketch_name': sketch_name,
                    'method': method,
                    'size': size,
                    'error_pct': error_pct,
                    'time_ms': time_ms,
                    'ops_per_second': ops_per_second
                })
    
    # Group data by method and sketch name
    method_sketch_pairs = set([(d['method'], d['sketch_name']) for d in plot_data])
    
    # Plot the data as a scatter plot
    markers = ['o', 's', '^', 'v', 'D', '*']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    sizes = [50, 100, 200, 300, 400]  # Size points by dataset size
    
    # Add text labels for each point
    for i, (method, sketch_name) in enumerate(method_sketch_pairs):
        # Filter data for this method and sketch
        method_data = [d for d in plot_data if d['method'] == method and d['sketch_name'] == sketch_name]
        
        # Continue if no data
        if not method_data:
            continue
            
        # Calculate averages across all sizes
        avg_error = sum(d['error_pct'] for d in method_data) / len(method_data)
        avg_time = sum(d['time_ms'] for d in method_data) / len(method_data)
        
        # Plot with different marker for each sketch and color for each method
        method_idx = [d['method'] for d in plot_data].index(method) if method in [d['method'] for d in plot_data] else 0
        sketch_idx = [d['sketch_name'] for d in plot_data].index(sketch_name) if sketch_name in [d['sketch_name'] for d in plot_data] else 0
        
        marker_idx = sketch_idx % len(markers)
        color_idx = method_idx % len(colors)
        
        label = f"{sketch_name} ({method})"
        
        # Plot individual points
        for d in method_data:
            size_idx = dataset_sizes.index(d['size']) if d['size'] in dataset_sizes else 0
            point_size = sizes[size_idx % len(sizes)]
            
            ax.scatter(d['time_ms'], d['error_pct'], 
                     s=point_size, 
                     marker=markers[marker_idx], 
                     color=colors[color_idx],
                     alpha=0.7)
        
        # Plot average with a larger marker and label
        ax.scatter(avg_time, avg_error, 
                 s=250, 
                 marker=markers[marker_idx], 
                 color=colors[color_idx],
                 edgecolors='black',
                 linewidth=1.5,
                 label=label)
        
        # Add method name as annotation
        ax.annotate(label, 
                  (avg_time, avg_error),
                  textcoords="offset points",
                  xytext=(5, 5),
                  fontsize=8)
    
    # Set axis labels and title
    ax.set_xlabel('Average Execution Time (ms)')
    ax.set_ylabel('Average Error (%)')
    ax.set_title('Performance vs Accuracy Tradeoff for Different Estimation Methods')
    
    # Set axis scales to logarithmic
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Add grid
    ax.grid(True, which='both', linestyle='--', alpha=0.5)
    
    # Add legend
    ax.legend(loc='upper right', fontsize='small')
    
    # Add annotations
    ax.text(0.02, 0.98, "Better", transform=ax.transAxes, fontsize=14, 
            verticalalignment='top', horizontalalignment='left',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='green', alpha=0.3))
    
    ax.text(0.98, 0.02, "Worse", transform=ax.transAxes, fontsize=14, 
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='red', alpha=0.3))
    
    # Add curved arrow from "Better" to "Worse"
    ax.annotate("", xy=(0.98, 0.02), xytext=(0.02, 0.98), 
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=dict(arrowstyle="fancy", connectionstyle="arc3,rad=0.3", 
                                facecolor='gray', edgecolor='gray', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(os.path.join(run_dir, 'performance_vs_accuracy.png'))
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
        sketch_classes.extend([RustHLLWrapper, rust_hll.RustHLL])
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
    
    # Generate timestamp and create results subdirectory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    script_name = os.path.splitext(os.path.basename(__file__))[0]
    run_dir = os.path.join(RESULTS_DIR, f'{script_name}_{timestamp}')
    os.makedirs(run_dir, exist_ok=True)
    print(f"\nSaving results to: {run_dir}")
    
    # Save results to file
    save_results_to_file({
        'ADD': add_results,
        'ADD_BATCH': batch_results,
        'CARDINALITY': cardinality_results,
        'MERGE': merge_results,
        'ACCURACY': accuracy_results
    }, os.path.join(run_dir, 'results.txt'))
    
    # Plot all results in a single figure
    plot_all_results(add_results, batch_results, cardinality_results, merge_results, accuracy_results, run_dir)
    
    # Plot accuracy results in a separate figure with subplots
    plot_accuracy_results(accuracy_results, run_dir)
    
    # Plot estimation methods comparison
    plot_estimation_methods(cardinality_results, run_dir)
    
    # Plot estimation bias
    plot_estimation_bias(cardinality_results, run_dir)
    
    # Plot performance vs accuracy tradeoff
    plot_performance_vs_accuracy(cardinality_results, run_dir)

if __name__ == "__main__":
    run_benchmarks() 