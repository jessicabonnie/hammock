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
PRECISION_VALUES = [8, 10, 12, 14, 16, 18, 20]
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
                    sketch.add(str(item))
                else:  # RustHyperLogLog
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

def plot_estimation_methods(cardinality_results, base_filename):
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
            
            for method, estimate in result['estimates'].items():
                # Initialize method in dictionaries if not present
                if method not in method_data['errors']:
                    method_data['errors'][method] = []
                    method_data['times'][method] = []
                
                # Calculate error
                if isinstance(estimate, (int, float, np.number)):
                    error_pct = abs(float(estimate) - actual) / actual * 100
                    method_data['errors'][method].append((sketch_name, size, error_pct))
                
                # Get time
                if method in result['times']:
                    time_taken = result['times'][method] * 1000  # Convert to milliseconds
                    method_data['times'][method].append((sketch_name, size, time_taken))
    
    # Plot accuracy (error percentage)
    markers = ['o', 's', '^', 'v', 'D', '*']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    # Group by method
    for i, (method, data_points) in enumerate(method_data['errors'].items()):
        # Group by sketch name
        sketch_dict = {}
        for sketch_name, size, error in data_points:
            if sketch_name not in sketch_dict:
                sketch_dict[sketch_name] = {'sizes': [], 'errors': []}
            sketch_dict[sketch_name]['sizes'].append(size)
            sketch_dict[sketch_name]['errors'].append(error)
        
        # Plot each sketch as a separate line
        for j, (sketch_name, data) in enumerate(sketch_dict.items()):
            marker_idx = j % len(markers)
            color_idx = i % len(colors)
            label = f"{sketch_name} ({method})"
            
            if data['sizes'] and data['errors']:
                # Sort by size to ensure correct line order
                points = sorted(zip(data['sizes'], data['errors']))
                sizes, errors = zip(*points)
                
                ax1.plot(sizes, errors, marker=markers[marker_idx], color=colors[color_idx], 
                          label=label, linestyle='-')
    
    ax1.set_title('Estimation Method Accuracy')
    ax1.set_xlabel('Dataset Size')
    ax1.set_ylabel('Error (%)')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.grid(True)
    ax1.legend(loc='upper right', fontsize='small')
    
    # Plot speed (time taken)
    for i, (method, data_points) in enumerate(method_data['times'].items()):
        # Group by sketch name
        sketch_dict = {}
        for sketch_name, size, time_taken in data_points:
            if sketch_name not in sketch_dict:
                sketch_dict[sketch_name] = {'sizes': [], 'times': []}
            sketch_dict[sketch_name]['sizes'].append(size)
            sketch_dict[sketch_name]['times'].append(time_taken)
        
        # Plot each sketch as a separate line
        for j, (sketch_name, data) in enumerate(sketch_dict.items()):
            marker_idx = j % len(markers)
            color_idx = i % len(colors)
            label = f"{sketch_name} ({method})"
            
            if data['sizes'] and data['times']:
                # Sort by size to ensure correct line order
                points = sorted(zip(data['sizes'], data['times']))
                sizes, times = zip(*points)
                
                ax2.plot(sizes, times, marker=markers[marker_idx], color=colors[color_idx], 
                          label=label, linestyle='-')
    
    ax2.set_title('Estimation Method Performance')
    ax2.set_xlabel('Dataset Size')
    ax2.set_ylabel('Execution Time (ms)')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.grid(True)
    ax2.legend(loc='upper left', fontsize='small')
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, f'{base_filename}_estimation_methods.png'))
    plt.close()
    
    # Create a more detailed plot showing actual vs estimated values
    fig, ax = plt.subplots(figsize=(15, 10))
    
    bar_width = 0.15
    opacity = 0.8
    index = np.arange(len(dataset_sizes))
    
    # Plot actual values as a line
    ax.plot(index, dataset_sizes, 'k--', linewidth=2, label='Actual')
    
    # Plot estimated values as bars grouped by sketch
    offset = -bar_width * (len(method_data['errors']) / 2)
    
    for i, (method, data_points) in enumerate(method_data['errors'].items()):
        # Group by sketch name
        sketch_dict = {}
        for sketch_name, size, _ in data_points:
            if sketch_name not in sketch_dict:
                sketch_dict[sketch_name] = {'sizes': [], 'estimates': []}
            
            # Find the estimate value from the original results
            if size in cardinality_results[sketch_name]:
                estimate = cardinality_results[sketch_name][size]['estimates'].get(method)
                if estimate is not None:
                    sketch_dict[sketch_name]['sizes'].append(size)
                    sketch_dict[sketch_name]['estimates'].append(float(estimate))
        
        # Plot each sketch as a separate set of bars
        for j, (sketch_name, data) in enumerate(sketch_dict.items()):
            if not data['sizes'] or not data['estimates']:
                continue
                
            # Create mapping from size to index position
            size_to_idx = {size: idx for idx, size in enumerate(dataset_sizes)}
            
            # Create arrays for positions and values
            positions = [size_to_idx[size] for size in data['sizes']]
            values = data['estimates']
            
            color_idx = i % len(colors)
            label = f"{sketch_name} ({method})"
            
            current_offset = offset + i * bar_width
            ax.bar(index[positions] + current_offset, values, bar_width,
                   alpha=opacity, color=colors[color_idx], label=label)
    
    ax.set_xlabel('Dataset Size')
    ax.set_ylabel('Cardinality Estimate')
    ax.set_title('Actual vs Estimated Cardinality by Method')
    ax.set_xticks(index)
    ax.set_xticklabels([f"{size:,}" for size in dataset_sizes], rotation=45)
    ax.set_yscale('log')
    ax.legend(loc='upper left', fontsize='small')
    
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, f'{base_filename}_cardinality_estimates.png'))
    plt.close()

def plot_estimation_bias(cardinality_results, base_filename):
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
            
            for method, estimate in result['estimates'].items():
                # Initialize method in dictionary if not present
                if method not in method_data:
                    method_data[method] = []
                
                # Calculate bias ratio (> 1 means overestimation, < 1 means underestimation)
                if isinstance(estimate, (int, float, np.number)):
                    ratio = float(estimate) / actual
                    method_data[method].append((sketch_name, size, ratio))
    
    # Plot bias ratio (logarithmic, so 0 means unbiased, > 0 means overestimation, < 0 means underestimation)
    markers = ['o', 's', '^', 'v', 'D', '*']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    # Plot a horizontal line at y=1 (unbiased)
    ax.axhline(y=1.0, color='k', linestyle='--', alpha=0.3, label="Unbiased (ratio=1)")
    
    # Group by method
    for i, (method, data_points) in enumerate(method_data.items()):
        # Group by sketch name
        sketch_dict = {}
        for sketch_name, size, ratio in data_points:
            if sketch_name not in sketch_dict:
                sketch_dict[sketch_name] = {'sizes': [], 'ratios': []}
            sketch_dict[sketch_name]['sizes'].append(size)
            sketch_dict[sketch_name]['ratios'].append(ratio)
        
        # Plot each sketch as a separate line
        for j, (sketch_name, data) in enumerate(sketch_dict.items()):
            marker_idx = j % len(markers)
            color_idx = i % len(colors)
            label = f"{sketch_name} ({method})"
            
            if data['sizes'] and data['ratios']:
                # Sort by size to ensure correct line order
                points = sorted(zip(data['sizes'], data['ratios']))
                sizes, ratios = zip(*points)
                
                ax.plot(sizes, ratios, marker=markers[marker_idx], color=colors[color_idx], 
                          label=label, linestyle='-')
    
    ax.set_title('Estimation Method Bias (Estimate/Actual Ratio)')
    ax.set_xlabel('Dataset Size')
    ax.set_ylabel('Estimate/Actual Ratio')
    ax.set_xscale('log')
    ax.grid(True)
    
    # Add region highlighting the 0.95-1.05 range (5% error margin)
    ax.axhspan(0.95, 1.05, alpha=0.2, color='green', label="±5% error band")
    
    # Refine y-axis limits to focus on the relevant range
    y_min = min([min(d['ratios']) for d in [data for data in sketch_dict.values() 
                                           for sketch_dict in method_data.values()]])
    y_max = max([max(d['ratios']) for d in [data for data in sketch_dict.values() 
                                           for sketch_dict in method_data.values()]])
    
    # Add some padding to the limits
    padding = (y_max - y_min) * 0.1
    ax.set_ylim([max(0, y_min - padding), y_max + padding])
    
    ax.legend(loc='upper right', fontsize='small')
    
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, f'{base_filename}_estimation_bias.png'))
    plt.close()
    
    # Create a boxplot to show the distribution of bias across all dataset sizes
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
    ax.axhspan(0.95, 1.05, alpha=0.2, color='green')
    
    ax.set_title('Distribution of Estimation Bias by Method and Implementation')
    ax.set_ylabel('Estimate/Actual Ratio')
    ax.set_xticklabels(boxplot_labels, rotation=45, ha='right')
    ax.grid(True, axis='y')
    
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, f'{base_filename}_estimation_bias_boxplot.png'))
    plt.close()

def plot_performance_vs_accuracy(cardinality_results, base_filename):
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
            
            for method, estimate in result['estimates'].items():
                # Skip if not a valid estimate or time
                if not isinstance(estimate, (int, float, np.number)) or method not in result['times']:
                    continue
                
                # Calculate error percentage and performance (ops per second)
                error_pct = abs(float(estimate) - actual) / actual * 100
                time_ms = result['times'][method] * 1000  # Convert to milliseconds
                ops_per_second = 1 / result['times'][method] if result['times'][method] > 0 else float('inf')
                
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
    plt.savefig(os.path.join(RESULTS_DIR, f'{base_filename}_performance_vs_accuracy.png'))
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
    
    # Plot estimation methods comparison
    plot_estimation_methods(cardinality_results, base_filename)
    
    # Plot estimation bias
    plot_estimation_bias(cardinality_results, base_filename)
    
    # Plot performance vs accuracy tradeoff
    plot_performance_vs_accuracy(cardinality_results, base_filename)

if __name__ == "__main__":
    run_benchmarks() 