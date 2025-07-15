#!/usr/bin/env python3
"""
Enhanced HyperLogLog benchmark with integrated cProfile profiling.
"""
import time
import random
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from hammock.lib.rusthll_compat import RustHLLWrapper
from hammock.lib.hyperloglog import HyperLogLog
from datetime import datetime

# Import profiling utilities
from profiling_utils import ProfiledBenchmark, profile_function, BenchmarkProfiler, analyze_profile_file

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
NUM_ITEMS = [1000, 10000, 100000, 1000000]
METHODS = ["original", "ertl_mle", "fast_mle"]

def generate_data(size, unique_ratio=1.0):
    """Generate test data with controlled uniqueness."""
    unique_count = int(size * unique_ratio)
    unique_items = [f"unique_item_{i}" for i in range(unique_count)]
    
    if unique_ratio < 1.0:
        # Generate data with duplicates
        result = []
        for _ in range(size):
            if random.random() < unique_ratio:
                result.append(random.choice(unique_items))
            else:
                result.append(random.choice(unique_items[:unique_count//2]))
        return result
    else:
        return [f"item_{i}" for i in range(size)]

@profile_function(save_profile=True, output_dir=RESULTS_DIR)
def benchmark_add_profiled(sketch_classes, precision=14, num_items=100000):
    """Benchmark adding items to different HLL implementations with profiling."""
    results = {}
    
    for sketch_class in sketch_classes:
        sketch_name = sketch_class.__name__
        
        # Create sketch
        if sketch_name == "RustHLLWrapper":
            sketch = sketch_class(precision=precision)
        else:
            sketch = sketch_class(precision=precision, hash_size=32, debug=False)
        
        results[sketch_name] = {}
        
        # Generate data
        items = generate_data(num_items)
        
        # Benchmark adding items
        start_time = time.time()
        for item in items:
            if sketch_name == "HyperLogLog":
                sketch.add_string(item)
            elif sketch_name == "RustHLLWrapper":
                sketch.add(item)
            else:
                sketch.add(item)
        end_time = time.time()
        
        time_taken = end_time - start_time
        items_per_second = num_items / time_taken if time_taken > 0 else float('inf')
        
        results[sketch_name] = {
            'time': time_taken,
            'items_per_second': items_per_second,
            'cardinality_estimate': sketch.estimate_cardinality() if hasattr(sketch, 'estimate_cardinality') else sketch.cardinality()
        }
    
    return results

def benchmark_add_with_detailed_profiling(sketch_classes, precision=14, num_items=100000):
    """Benchmark with detailed profiling using context manager."""
    results = {}
    
    for sketch_class in sketch_classes:
        sketch_name = sketch_class.__name__
        
        with ProfiledBenchmark(f"add_benchmark_{sketch_name}", output_dir=RESULTS_DIR) as profiler:
            # Create sketch
            if sketch_name == "RustHLLWrapper":
                sketch = sketch_class(precision=precision)
            else:
                sketch = sketch_class(precision=precision, hash_size=32, debug=False)
            
            # Generate data
            items = generate_data(num_items)
            
            # Benchmark adding items
            start_time = time.time()
            for item in items:
                if sketch_name == "HyperLogLog":
                    sketch.add_string(item)
                elif sketch_name == "RustHLLWrapper":
                    sketch.add(item)
                else:
                    sketch.add(item)
            end_time = time.time()
            
            # Get cardinality estimate
            cardinality = sketch.estimate_cardinality() if hasattr(sketch, 'estimate_cardinality') else sketch.cardinality()
        
        # Get function-level statistics
        add_stats = profiler.get_function_stats('add')
        hash_stats = profiler.get_function_stats('hash')
        
        results[sketch_name] = {
            'time': end_time - start_time,
            'items_per_second': num_items / (end_time - start_time) if (end_time - start_time) > 0 else float('inf'),
            'cardinality_estimate': cardinality,
            'wall_time': profiler.wall_time,
            'memory_delta': profiler.memory_delta,
            'add_function_stats': add_stats,
            'hash_function_stats': hash_stats,
            'profile_stats': profiler.profile_stats
        }
    
    return results

def benchmark_cardinality_profiled(sketch_classes, precision=14, num_runs=5):
    """Benchmark cardinality estimation with profiling."""
    results = {}
    
    for sketch_class in sketch_classes:
        sketch_name = sketch_class.__name__
        results[sketch_name] = {}
        
        for method in METHODS:
            if sketch_name == "RustHLLWrapper" and method != "original":
                continue  # Rust wrapper might not support all methods
            
            method_results = []
            
            for run in range(num_runs):
                test_name = f"cardinality_{sketch_name}_{method}_run_{run}"
                
                with ProfiledBenchmark(test_name, output_dir=RESULTS_DIR) as profiler:
                    # Create sketch
                    if sketch_name == "RustHLLWrapper":
                        sketch = sketch_class(precision=precision)
                    else:
                        sketch = sketch_class(precision=precision, hash_size=32, debug=False)
                    
                    # Add data
                    items = generate_data(50000, unique_ratio=0.8)
                    for item in items:
                        if sketch_name == "HyperLogLog":
                            sketch.add_string(item)
                        elif sketch_name == "RustHLLWrapper":
                            sketch.add(item)
                        else:
                            sketch.add(item)
                    
                    # Benchmark estimation
                    start_time = time.time()
                    if method == "original":
                        estimate = sketch.estimate_cardinality() if hasattr(sketch, 'estimate_cardinality') else sketch.cardinality()
                    elif method == "ertl_mle" and hasattr(sketch, 'estimate_cardinality_ertl_mle'):
                        estimate = sketch.estimate_cardinality_ertl_mle()
                    elif method == "fast_mle" and hasattr(sketch, 'estimate_cardinality_fast_mle'):
                        estimate = sketch.estimate_cardinality_fast_mle()
                    else:
                        estimate = sketch.estimate_cardinality() if hasattr(sketch, 'estimate_cardinality') else sketch.cardinality()
                    end_time = time.time()
                
                true_cardinality = len(set(items))
                error = abs(estimate - true_cardinality) / true_cardinality * 100
                
                method_results.append({
                    'estimate': estimate,
                    'true_cardinality': true_cardinality,
                    'error_percent': error,
                    'estimation_time': end_time - start_time,
                    'wall_time': profiler.wall_time,
                    'memory_delta': profiler.memory_delta,
                    'profile_stats': profiler.profile_stats
                })
            
            results[sketch_name][method] = method_results
    
    return results

def run_comprehensive_hll_profiling():
    """Run comprehensive HLL profiling across different configurations."""
    
    # Available sketch classes
    sketch_classes = [HyperLogLog]
    if RUST_HLL_AVAILABLE:
        sketch_classes.append(RustHLLWrapper)
    
    # Create profiler
    profiler = BenchmarkProfiler(output_dir=RESULTS_DIR)
    
    # Define test configurations
    configurations = []
    for precision in [12, 16, 20]:
        for num_items in [10000, 100000]:
            configurations.append({
                'name': f'precision_{precision}_items_{num_items}',
                'precision': precision,
                'num_items': num_items,
                'sketch_classes': sketch_classes
            })
    
    def run_add_benchmark(precision, num_items, sketch_classes):
        """Run add benchmark for a specific configuration."""
        return benchmark_add_with_detailed_profiling(sketch_classes, precision, num_items)
    
    # Run profiled benchmarks
    all_results = []
    for config in configurations:
        results = profiler.run_profiled_benchmark(
            benchmark_func=run_add_benchmark,
            name=f"hll_add_{config['name']}",
            configurations=[config]
        )
        all_results.extend(results)
    
    # Generate summary report
    profiler.save_summary_report("hll_benchmark_summary.csv")
    
    # Generate plots
    profiler.plot_performance_comparison('wall_time')
    profiler.plot_performance_comparison('memory_delta')
    
    return all_results

def analyze_performance_bottlenecks(sketch_class, precision=14, num_items=100000):
    """Analyze performance bottlenecks in a specific sketch implementation."""
    
    sketch_name = sketch_class.__name__
    print(f"\nAnalyzing performance bottlenecks for {sketch_name}")
    print("="*60)
    
    with ProfiledBenchmark(f"bottleneck_analysis_{sketch_name}", output_dir=RESULTS_DIR) as profiler:
        # Create sketch
        if sketch_name == "RustHLLWrapper":
            sketch = sketch_class(precision=precision)
        else:
            sketch = sketch_class(precision=precision, hash_size=32, debug=False)
        
        # Generate data
        items = generate_data(num_items)
        
        # Add items with periodic profiling checkpoints
        chunk_size = num_items // 10
        for i in range(0, num_items, chunk_size):
            chunk = items[i:i + chunk_size]
            
            chunk_start = time.time()
            for item in chunk:
                if sketch_name == "HyperLogLog":
                    sketch.add_string(item)
                elif sketch_name == "RustHLLWrapper":
                    sketch.add(item)
                else:
                    sketch.add(item)
            chunk_end = time.time()
            
            print(f"  Chunk {i//chunk_size + 1}: {chunk_end - chunk_start:.4f}s")
        
        # Test cardinality estimation
        estimation_start = time.time()
        cardinality = sketch.estimate_cardinality() if hasattr(sketch, 'estimate_cardinality') else sketch.cardinality()
        estimation_end = time.time()
        
        print(f"  Cardinality estimation: {estimation_end - estimation_start:.4f}s")
        print(f"  Estimated cardinality: {cardinality}")
    
    print(f"\nTotal wall time: {profiler.wall_time:.4f}s")
    print(f"Memory delta: {profiler.memory_delta:.2f} MiB")
    
    # Analyze top functions
    if profiler.profile_stats:
        print("\nTop 10 functions by cumulative time:")
        sorted_stats = sorted(profiler.profile_stats.stats.items(), key=lambda x: x[1][3], reverse=True)
        for i, (func_name, (cc, nc, tt, ct, callers)) in enumerate(sorted_stats[:10]):
            filename, line_num, func = func_name
            print(f"  {i+1}. {func} ({filename}:{line_num})")
            print(f"     Calls: {cc}, Total time: {tt:.4f}s, Cumulative: {ct:.4f}s")
    
    return profiler.profile_stats

def compare_implementations():
    """Compare different HLL implementations with detailed profiling."""
    
    if not RUST_HLL_AVAILABLE:
        print("Rust HLL not available, skipping comparison")
        return
    
    print("\nComparing HLL implementations...")
    print("="*60)
    
    sketch_classes = [HyperLogLog, RustHLLWrapper]
    precision = 14
    num_items = 50000
    
    results = {}
    profile_files = []
    
    for sketch_class in sketch_classes:
        sketch_name = sketch_class.__name__
        
        # Run detailed profiling
        profile_stats = analyze_performance_bottlenecks(sketch_class, precision, num_items)
        
        # Save profile for comparison
        profile_file = os.path.join(RESULTS_DIR, f"comparison_{sketch_name}.prof")
        if profile_stats:
            profile_stats.dump_stats(profile_file)
            profile_files.append(profile_file)
    
    # Compare profiles
    if len(profile_files) >= 2:
        print("\nComparing profile files...")
        from profiling_utils import compare_profiles
        compare_profiles(profile_files)

if __name__ == "__main__":
    print("Running comprehensive HLL profiling...")
    
    # Run comprehensive profiling
    results = run_comprehensive_hll_profiling()
    
    # Analyze bottlenecks
    print("\n" + "="*60)
    print("Analyzing performance bottlenecks...")
    analyze_performance_bottlenecks(HyperLogLog, precision=14, num_items=100000)
    
    # Compare implementations if available
    if RUST_HLL_AVAILABLE:
        compare_implementations()
    
    print("\nHLL profiling completed! Check the results directory for detailed analysis.")
    print("Use analyze_profile_file() to examine specific profile files.") 