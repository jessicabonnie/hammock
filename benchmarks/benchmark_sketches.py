import time
import csv
import os
from datetime import datetime
from memory_profiler import memory_usage
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash
from hammock.lib.minimizer import MinimizerSketch
from hammock.lib.exact import ExactCounter
from hammock.lib.setsketch import SetSketch
import random
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore

# Create results directory if it doesn't exist
RESULTS_DIR = os.path.join(os.path.dirname(__file__), 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)

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

def run_test_case(sketch_class, param_name: str, param_value: int, 
                  set1_size: int, set2_size: int, set2_offset: int = 0):
    """Run a single test case and return results."""
    # Initialize sketches with the given parameter
    sketch1 = sketch_class(**{param_name: param_value})
    sketch2 = sketch_class(**{param_name: param_value})
    
    # Time the data addition
    start_add = time.time()
    for i in range(set1_size):
        sketch1.add_string(str(i))
    for i in range(set2_size):
        sketch2.add_string(str(i + set2_offset))
    end_add = time.time()
    
    # Calculate expected Jaccard
    intersection = min(set1_size, max(0, set2_size - set2_offset))
    union = set1_size + set2_size - intersection
    expected_jaccard = intersection / union if union > 0 else 0.0
    
    # Time the Jaccard calculation
    start_jaccard = time.time()
    jaccard = sketch1.estimate_jaccard(sketch2)
    end_jaccard = time.time()
    
    return {
        'add_time': end_add - start_add,
        'jaccard_time': end_jaccard - start_jaccard,
        'expected_jaccard': expected_jaccard,
        'calculated_jaccard': jaccard,
        'absolute_error': abs(jaccard - expected_jaccard)
    }

def benchmark_sketch(sketch_class, param_name: str, param_values: list, test_cases: list):
    """Benchmark a sketch type across different parameter values."""
    results = []
    
    for param_value in param_values:
        print(f"\n{'='*60}")
        print(f"Testing {sketch_class.__name__} with {param_name}={param_value}")
        print('='*60)
        
        for name, desc, set1_size, set2_size, set2_offset in test_cases:
            print(f"\nRunning: {name}")
            print(desc)
            
            # Measure memory and run test
            start_mem = memory_usage()[0]
            test_results = run_test_case(
                sketch_class, param_name, param_value,
                set1_size, set2_size, set2_offset
            )
            end_mem = memory_usage()[0]
            
            # Combine results
            result = {
                'sketch_type': sketch_class.__name__,
                param_name: param_value,
                'test_name': name,
                'set1_size': set1_size,
                'set2_size': set2_size,
                'memory_delta': end_mem - start_mem,
                **test_results
            }
            
            results.append(result)
            print(f"Memory usage: {result['memory_delta']:.2f} MiB")
            print(f"Add time: {result['add_time']:.3f}s")
            print(f"Jaccard time: {result['jaccard_time']:.3f}s")
            print(f"Expected Jaccard: {result['expected_jaccard']:.3f}")
            print(f"Calculated Jaccard: {result['calculated_jaccard']:.3f}")
            print(f"Absolute error: {result['absolute_error']:.3f}")
    
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
    
    # Save to file in results directory with timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    filename = title.lower().replace(' ', '_') + f'_{timestamp}.png'
    filepath = os.path.join(RESULTS_DIR, filename)
    plt.savefig(filepath)
    print(f"Saved plot to {filepath}")
    
    plt.close()

def plot_accuracy_results(results, base_filename):
    """Create a single figure with accuracy subplots for different item counts."""
    # Create a figure with subplots
    fig = plt.figure(figsize=(20, 15))
    
    # Get the number of item counts and create a grid of subplots
    num_items = len(NUM_ITEMS)
    num_cols = 3  # Number of columns in the grid
    num_rows = (num_items + num_cols - 1) // num_cols  # Calculate number of rows needed
    
    for idx, num in enumerate(NUM_ITEMS):
        ax = fig.add_subplot(num_rows, num_cols, idx + 1)
        
        # Plot accuracy for each sketch type
        for sketch_name, data in results.items():
            precisions = []
            errors = []
            for p in PRECISION_VALUES:
                if p in data and num in data[p]:
                    precisions.append(p)
                    errors.append(data[p][num]['error'])
            
            if precisions and errors:  # Only plot if we have data
                ax.plot(precisions, errors, marker='o', label=sketch_name)
        
        ax.set_title(f'Accuracy at {num:,} Items')
        ax.set_xlabel('Precision')
        ax.set_ylabel('Error (%)')
        ax.grid(True)
        ax.legend(loc='upper right')
    
    # Adjust layout to prevent overlap
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(os.path.join(RESULTS_DIR, f'{base_filename}_accuracy.png'))
    plt.close()

def run_benchmarks():
    """Run all benchmarks and save results."""
    # Define sketch classes to benchmark
    sketch_classes = [HyperLogLog, MinHash, SetSketch]
    
    # Run benchmarks
    add_results = benchmark_add(sketch_classes)
    batch_results = benchmark_add_batch(sketch_classes)
    cardinality_results = benchmark_cardinality(sketch_classes)
    merge_results = benchmark_merge(sketch_classes)
    accuracy_results = benchmark_accuracy(sketch_classes)
    
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