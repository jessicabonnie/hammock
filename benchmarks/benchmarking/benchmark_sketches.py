import time
import csv
import os
from datetime import datetime
from memory_profiler import memory_usage
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash

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

def run_benchmarks():
    """Run all benchmarks and save results."""
    # Define test cases
    test_cases = [
        (
            "Very sparse perfect overlap",
            "10 integers with perfect overlap",
            10, 10, 0
        ),
        (
            "Sparse vs Dense",
            "100 vs 1000 integers with partial overlap",
            100, 1000, 0
        ),
        (
            "Medium density overlap",
            "1000 integers with 90% overlap",
            1000, 1000, 100
        ),
        (
            "Dense high overlap",
            "10000 integers with 90% overlap",
            10000, 10000, 1000
        ),
        (
            "Very dense partial overlap",
            "100000 integers with 50% overlap",
            100000, 100000, 50000
        ),
        (
            "Extremely dense minimal overlap",
            "1000000 integers with 5% overlap",
            1000000, 1000000, 900000
        )
    ]
    
    # Run benchmarks for each sketch type
    results = []
    results.extend(benchmark_sketch(
        MinHash, 'num_hashes', 
        [16, 32, 64, 128, 256, 512, 1024],
        test_cases
    ))
    results.extend(benchmark_sketch(
        HyperLogLog, 'precision',
        [4, 6, 8, 10, 12, 14],
        test_cases
    ))
    
    # Create results directory if it doesn't exist
    os.makedirs('test_results', exist_ok=True)
    
    # Write results to CSV
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    filename = f'test_results/sketch_benchmark_{timestamp}.csv'
    
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    
    print(f"\nResults written to {filename}")

if __name__ == "__main__":
    run_benchmarks() 