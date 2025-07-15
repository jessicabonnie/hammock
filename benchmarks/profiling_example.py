#!/usr/bin/env python3
"""
Simple example demonstrating cProfile integration with existing benchmarks.
"""

import os
import sys
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash
from profiling_utils import ProfiledBenchmark, profile_function, BenchmarkProfiler

def simple_sketch_test():
    """Simple sketch test for demonstration."""
    # Create two sketches
    sketch1 = HyperLogLog(precision=12)
    sketch2 = HyperLogLog(precision=12)
    
    # Add some data
    for i in range(10000):
        sketch1.add_string(f"item_{i}")
        sketch2.add_string(f"item_{i + 5000}")  # 50% overlap
    
    # Calculate Jaccard similarity
    jaccard = sketch1.estimate_jaccard(sketch2)
    
    return {
        'jaccard': jaccard,
        'cardinality1': sketch1.estimate_cardinality(),
        'cardinality2': sketch2.estimate_cardinality()
    }

@profile_function(save_profile=True, output_dir="results")
def profiled_sketch_test():
    """The same test but with profiling decorator."""
    return simple_sketch_test()

def main():
    print("cProfile Integration Example")
    print("="*40)
    
    # Example 1: Using the context manager
    print("\n1. Using ProfiledBenchmark context manager:")
    
    with ProfiledBenchmark("hyperloglog_test", output_dir="results") as profiler:
        result = simple_sketch_test()
    
    print(f"   Jaccard similarity: {result['jaccard']:.3f}")
    print(f"   Wall time: {profiler.wall_time:.4f}s")
    print(f"   Memory delta: {profiler.memory_delta:.2f} MiB")
    
    # Example 2: Using the decorator
    print("\n2. Using @profile_function decorator:")
    
    result = profiled_sketch_test()
    if 'profile_info' in result:
        print(f"   Jaccard similarity: {result['jaccard']:.3f}")
        print(f"   Wall time: {result['profile_info']['wall_time']:.4f}s")
        print(f"   Memory delta: {result['profile_info']['memory_delta']:.2f} MiB")
    
    # Example 3: Using BenchmarkProfiler for multiple configurations
    print("\n3. Using BenchmarkProfiler for multiple configurations:")
    
    profiler = BenchmarkProfiler(output_dir="results")
    
    # Define test configurations
    configurations = [
        {
            'name': 'precision_12',
            'sketch_class': HyperLogLog,
            'precision': 12,
            'num_items': 5000
        },
        {
            'name': 'precision_16',
            'sketch_class': HyperLogLog,
            'precision': 16,
            'num_items': 5000
        },
        {
            'name': 'minhash_k100',
            'sketch_class': MinHash,
            'num_hashes': 100,
            'num_items': 5000
        }
    ]
    
    def run_config_test(sketch_class, num_items, **kwargs):
        """Run test for a specific configuration."""
        # Create sketch with parameters
        if sketch_class == HyperLogLog:
            sketch = sketch_class(precision=kwargs['precision'])
        else:
            sketch = sketch_class(num_hashes=kwargs['num_hashes'])
        
        # Add data
        for i in range(num_items):
            sketch.add_string(f"item_{i}")
        
        # Get cardinality estimate
        cardinality = sketch.estimate_cardinality()
        
        return {
            'cardinality': cardinality,
            'error': abs(cardinality - num_items) / num_items * 100
        }
    
    # Run profiled benchmarks
    results = profiler.run_profiled_benchmark(
        benchmark_func=run_config_test,
        name="sketch_comparison",
        configurations=configurations
    )
    
    print(f"   Ran {len(results)} configurations")
    for result in results:
        config = result['config']
        print(f"   {config['name']}: {result['wall_time']:.4f}s, "
              f"error: {result['benchmark_result']['error']:.2f}%")
    
    # Generate summary report
    profiler.save_summary_report("example_summary.csv")
    
    # Generate performance plots
    profiler.plot_performance_comparison('wall_time', save_plot=True)
    
    print("\n4. Profile files generated:")
    results_dir = "results"
    if os.path.exists(results_dir):
        profile_files = [f for f in os.listdir(results_dir) if f.endswith('.prof')]
        for pf in profile_files:
            print(f"   {pf}")
    
    print("\nExample completed! Check the 'results' directory for detailed profiling data.")
    print("You can analyze profile files using:")
    print("   from profiling_utils import analyze_profile_file")
    print("   analyze_profile_file('results/your_profile.prof')")

if __name__ == "__main__":
    main() 