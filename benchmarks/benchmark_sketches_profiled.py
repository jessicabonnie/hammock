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
import numpy as np
import matplotlib.pyplot as plt

# Import our new profiling utilities
from profiling_utils import ProfiledBenchmark, profile_function, BenchmarkProfiler

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

@profile_function(save_profile=True, output_dir=RESULTS_DIR)
def run_test_case_profiled(sketch_class, param_name: str, param_value: int, 
                          set1_size: int, set2_size: int, set2_offset: int = 0):
    """Run a single test case with profiling."""
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

def run_test_case_with_profiling(sketch_class, param_name: str, param_value: int, 
                                set1_size: int, set2_size: int, set2_offset: int = 0):
    """Run a single test case with detailed profiling using context manager."""
    
    test_name = f"{sketch_class.__name__}_{param_name}_{param_value}"
    
    with ProfiledBenchmark(test_name, output_dir=RESULTS_DIR) as profiler:
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
    
    # Get function-level statistics
    add_stats = profiler.get_function_stats('add_string')
    jaccard_stats = profiler.get_function_stats('estimate_jaccard')
    
    return {
        'add_time': end_add - start_add,
        'jaccard_time': end_jaccard - start_jaccard,
        'expected_jaccard': expected_jaccard,
        'calculated_jaccard': jaccard,
        'absolute_error': abs(jaccard - expected_jaccard),
        'wall_time': profiler.wall_time,
        'memory_delta': profiler.memory_delta,
        'add_function_stats': add_stats,
        'jaccard_function_stats': jaccard_stats,
        'profile_stats': profiler.profile_stats
    }

def benchmark_sketch_profiled(sketch_class, param_name: str, param_values: list, test_cases: list):
    """Enhanced benchmark with profiling capabilities."""
    results = []
    
    for param_value in param_values:
        print(f"\n{'='*60}")
        print(f"Testing {sketch_class.__name__} with {param_name}={param_value}")
        print('='*60)
        
        for name, desc, set1_size, set2_size, set2_offset in test_cases:
            print(f"\nRunning: {name}")
            print(desc)
            
            # Run with detailed profiling
            result = run_test_case_with_profiling(
                sketch_class, param_name, param_value,
                set1_size, set2_size, set2_offset
            )
            
            # Combine results
            combined_result = {
                'sketch_type': sketch_class.__name__,
                param_name: param_value,
                'test_name': name,
                'set1_size': set1_size,
                'set2_size': set2_size,
                **result
            }
            
            results.append(combined_result)
            
            # Print results
            print(f"Wall time: {result['wall_time']:.4f}s")
            print(f"Memory usage: {result['memory_delta']:.2f} MiB")
            print(f"Add time: {result['add_time']:.3f}s")
            print(f"Jaccard time: {result['jaccard_time']:.3f}s")
            print(f"Expected Jaccard: {result['expected_jaccard']:.3f}")
            print(f"Calculated Jaccard: {result['calculated_jaccard']:.3f}")
            print(f"Absolute error: {result['absolute_error']:.3f}")
            
            # Print top functions from profiling
            if result['profile_stats']:
                top_funcs = result['profile_stats'].stats
                if top_funcs:
                    sorted_funcs = sorted(top_funcs.items(), key=lambda x: x[1][3], reverse=True)
                    print(f"Top function by cumulative time: {sorted_funcs[0][0]}")
    
    return results

def run_comprehensive_profiled_benchmark():
    """Run comprehensive benchmark with profiling for all sketch types."""
    
    # Create benchmark profiler
    profiler = BenchmarkProfiler(output_dir=RESULTS_DIR)
    
    # Define sketch types and their parameters
    sketch_configs = [
        {
            'name': 'HyperLogLog_precision_12',
            'sketch_class': HyperLogLog,
            'param_name': 'precision',
            'param_value': 12
        },
        {
            'name': 'HyperLogLog_precision_16',
            'sketch_class': HyperLogLog,
            'param_name': 'precision',
            'param_value': 16
        },
        {
            'name': 'MinHash_k_100',
            'sketch_class': MinHash,
            'param_name': 'k',
            'param_value': 100
        },
        {
            'name': 'MinHash_k_200',
            'sketch_class': MinHash,
            'param_name': 'k',
            'param_value': 200
        }
    ]
    
    # Test cases: (name, description, set1_size, set2_size, set2_offset)
    test_cases = [
        ("identical_sets", "Two identical sets of 1000 elements", 1000, 1000, 0),
        ("disjoint_sets", "Two disjoint sets of 1000 elements", 1000, 1000, 1000),
        ("half_overlap", "Two sets with 50% overlap", 1000, 1000, 500),
        ("small_large", "Small set vs large set", 100, 2000, 0),
    ]
    
    def run_single_config(sketch_class, param_name, param_value, test_case):
        """Run a single configuration."""
        name, desc, set1_size, set2_size, set2_offset = test_case
        return run_test_case_with_profiling(
            sketch_class, param_name, param_value,
            set1_size, set2_size, set2_offset
        )
    
    # Run all configurations
    all_results = []
    for config in sketch_configs:
        for test_case in test_cases:
            config_with_test = {
                **config,
                'test_case': test_case
            }
            
            results = profiler.run_profiled_benchmark(
                benchmark_func=run_single_config,
                name=f"{config['name']}_{test_case[0]}",
                configurations=[config_with_test]
            )
            all_results.extend(results)
    
    # Generate summary report
    profiler.save_summary_report()
    
    # Generate performance comparison plots
    profiler.plot_performance_comparison('wall_time')
    profiler.plot_performance_comparison('memory_delta')
    
    return all_results

def analyze_sketch_performance(sketch_class, precision_values=[12, 16, 20]):
    """Analyze performance characteristics of a specific sketch type."""
    
    # Create configurations for different precision values
    configurations = []
    for precision in precision_values:
        configurations.append({
            'name': f'precision_{precision}',
            'sketch_class': sketch_class,
            'param_name': 'precision',
            'param_value': precision,
            'test_case': ("performance_test", "Performance analysis", 10000, 10000, 5000)
        })
    
    # Run profiled benchmark
    profiler = BenchmarkProfiler(output_dir=RESULTS_DIR)
    
    def run_performance_test(sketch_class, param_name, param_value, test_case):
        """Run performance test."""
        name, desc, set1_size, set2_size, set2_offset = test_case
        return run_test_case_with_profiling(
            sketch_class, param_name, param_value,
            set1_size, set2_size, set2_offset
        )
    
    results = profiler.run_profiled_benchmark(
        benchmark_func=run_performance_test,
        name=f"{sketch_class.__name__}_performance_analysis",
        configurations=configurations
    )
    
    # Print detailed analysis
    print(f"\nPerformance Analysis for {sketch_class.__name__}")
    print("="*60)
    
    for result in results:
        config = result['config']
        print(f"\nPrecision {config['param_value']}:")
        print(f"  Wall time: {result['wall_time']:.4f}s")
        print(f"  Memory delta: {result['memory_delta']:.2f} MiB")
        
        # Get function-level statistics
        benchmark_result = result['benchmark_result']
        if benchmark_result.get('add_function_stats'):
            add_stats = benchmark_result['add_function_stats']
            print(f"  Add function calls: {sum(stats['call_count'] for stats in add_stats.values())}")
            
        if benchmark_result.get('jaccard_function_stats'):
            jaccard_stats = benchmark_result['jaccard_function_stats']
            print(f"  Jaccard function calls: {sum(stats['call_count'] for stats in jaccard_stats.values())}")
    
    return results

if __name__ == "__main__":
    print("Running comprehensive profiled benchmark...")
    
    # Run comprehensive benchmark
    results = run_comprehensive_profiled_benchmark()
    
    # Run detailed analysis for HyperLogLog
    print("\n" + "="*60)
    print("Running detailed HyperLogLog analysis...")
    hll_results = analyze_sketch_performance(HyperLogLog, [12, 16, 20])
    
    print("\nBenchmark completed! Check the results directory for detailed profiling data.") 