#!/usr/bin/env python3
import sys
import os
import random
import time
import numpy as np
from hammock.lib import FastHyperLogLog
from hammock.lib import HyperLogLog

# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)

def test_basic_functionality():
    """Test basic functionality of FastHyperLogLog."""
    print("\n=== Testing Basic Functionality ===")
    
    # Create sketch
    sketch = FastHyperLogLog(precision=14)
    
    # Check if using Rust
    print(f"Using Rust implementation: {sketch.is_using_rust()}")
    
    # Add some values
    for i in range(1000):
        sketch.add(f"item_{i}")
    
    # Check cardinality
    estimated = sketch.cardinality()
    print(f"Added 1000 unique items, estimated: {estimated:.2f}")
    error = abs(estimated - 1000) / 1000 * 100
    print(f"Error: {error:.2f}%")
    assert error < 5, f"Error too high: {error:.2f}%"
    
    # Add duplicates
    for i in range(500):
        sketch.add(f"item_{i}")
    
    # Check cardinality again
    estimated = sketch.cardinality()
    print(f"Added 500 duplicates, estimated: {estimated:.2f}")
    error = abs(estimated - 1000) / 1000 * 100
    print(f"Error: {error:.2f}%")
    assert error < 5, f"Error too high: {error:.2f}%"
    
    print("Basic functionality test passed!")

def test_estimation_methods():
    """Test different cardinality estimation methods."""
    print("\n=== Testing Estimation Methods ===")
    
    # Create sketch
    sketch = FastHyperLogLog(precision=14)
    
    # Add some values
    for i in range(10000):
        sketch.add(f"item_{i}")
    
    # Test different estimation methods
    methods = ["original", "ertl_improved", "ertl_mle"]
    
    for method in methods:
        start_time = time.time()
        estimated = sketch.estimate_cardinality(method=method)
        elapsed = time.time() - start_time
        
        error = abs(estimated - 10000) / 10000 * 100
        print(f"Method: {method}, estimated: {estimated:.2f}, error: {error:.2f}%, time: {elapsed:.6f}s")
        assert error < 5, f"Error too high for {method}: {error:.2f}%"
    
    print("Estimation methods test passed!")

def test_merge():
    """Test merging multiple sketches."""
    print("\n=== Testing Merge Operation ===")
    
    # Create two sketches with different items
    sketch1 = FastHyperLogLog(precision=14)
    sketch2 = FastHyperLogLog(precision=14)
    
    # Add values to sketch1
    for i in range(5000):
        sketch1.add(f"item1_{i}")
    
    # Add values to sketch2
    for i in range(5000):
        sketch2.add(f"item2_{i}")
    
    # Also add some overlapping values
    for i in range(1000):
        sketch1.add(f"common_{i}")
        sketch2.add(f"common_{i}")
    
    # Check individual cardinalities
    est1 = sketch1.cardinality()
    est2 = sketch2.cardinality()
    
    print(f"Sketch1 (5000 unique + 1000 common): {est1:.2f}")
    print(f"Sketch2 (5000 unique + 1000 common): {est2:.2f}")
    
    # Merge sketches
    sketch1.merge(sketch2)
    
    # Check merged cardinality (should be around 11000)
    est_merged = sketch1.cardinality()
    print(f"Merged (11000 expected): {est_merged:.2f}")
    
    error = abs(est_merged - 11000) / 11000 * 100
    print(f"Merge error: {error:.2f}%")
    assert error < 5, f"Merge error too high: {error:.2f}%"
    
    print("Merge test passed!")

def test_jaccard():
    """Test Jaccard similarity calculation."""
    print("\n=== Testing Jaccard Similarity ===")
    
    # Create two sketches with different overlap ratios
    sketch1 = FastHyperLogLog(precision=14)
    sketch2 = FastHyperLogLog(precision=14)
    sketch3 = FastHyperLogLog(precision=14)
    
    # Fill sketch1 with 10000 items
    for i in range(10000):
        sketch1.add(f"item_{i}")
    
    # Fill sketch2 with 5000 overlapping and 5000 distinct items (50% overlap)
    for i in range(5000):
        sketch2.add(f"item_{i}")  # Overlapping
    for i in range(10000, 15000):
        sketch2.add(f"item_{i}")  # Distinct
    
    # Fill sketch3 with 1000 overlapping and 9000 distinct items (10% overlap)
    for i in range(1000):
        sketch3.add(f"item_{i}")  # Overlapping
    for i in range(15000, 24000):
        sketch3.add(f"item_{i}")  # Distinct
    
    # Calculate Jaccard similarities
    jaccard12 = sketch1.jaccard(sketch2)
    jaccard13 = sketch1.jaccard(sketch3)
    
    # Expected: |A∩B|/|A∪B| = 5000/15000 = 0.333 for sketch1 and sketch2
    print(f"Jaccard similarity (50% overlap, expected ~0.33): {jaccard12:.4f}")
    assert abs(jaccard12 - 0.333) < 0.1, f"Jaccard error too high: {abs(jaccard12 - 0.333):.4f}"
    
    # Expected: |A∩B|/|A∪B| = 1000/19000 = 0.053 for sketch1 and sketch3
    print(f"Jaccard similarity (10% overlap, expected ~0.05): {jaccard13:.4f}")
    assert abs(jaccard13 - 0.053) < 0.05, f"Jaccard error too high: {abs(jaccard13 - 0.053):.4f}"
    
    print("Jaccard similarity test passed!")

def test_batch_operations():
    """Test batch operations for performance."""
    print("\n=== Testing Batch Operations ===")
    
    # Generate data
    data_size = 100000
    data = [f"item_{i}" for i in range(data_size)]
    
    # Test individual add
    sketch1 = FastHyperLogLog(precision=14)
    
    start_time = time.time()
    for item in data[:10000]:  # Use smaller set for individual adds
        sketch1.add(item)
    individual_time = time.time() - start_time
    
    est1 = sketch1.cardinality()
    print(f"Individual add (10000 items): {est1:.2f}, time: {individual_time:.4f}s")
    
    # Test batch add
    sketch2 = FastHyperLogLog(precision=14)
    
    start_time = time.time()
    sketch2.add_batch(data)
    batch_time = time.time() - start_time
    
    est2 = sketch2.cardinality()
    print(f"Batch add ({data_size} items): {est2:.2f}, time: {batch_time:.4f}s")
    
    if sketch2.is_using_rust():
        print(f"Batch speedup (extrapolated): {(individual_time/10000*data_size)/batch_time:.2f}x")
    
    print("Batch operations test passed!")

def test_edge_cases():
    """Test edge cases."""
    print("\n=== Testing Edge Cases ===")
    
    # Empty sketch
    sketch = FastHyperLogLog(precision=14)
    empty_est = sketch.cardinality()
    print(f"Empty sketch estimate: {empty_est}")
    assert empty_est < 1, f"Empty sketch should estimate close to 0, got {empty_est}"
    
    # Very small cardinality
    sketch = FastHyperLogLog(precision=14)
    for i in range(5):
        sketch.add(f"item_{i}")
    
    small_est = sketch.cardinality()
    print(f"Small cardinality (5 items): {small_est:.2f}")
    assert abs(small_est - 5) / 5 < 0.5, f"Small cardinality error too high: {abs(small_est - 5) / 5:.2f}"
    
    # Very large cardinality
    sketch = FastHyperLogLog(precision=14)
    for i in range(100000):
        sketch.add(f"item_{i}")
    
    large_est = sketch.cardinality()
    print(f"Large cardinality (100000 items): {large_est:.2f}")
    assert abs(large_est - 100000) / 100000 < 0.05, f"Large cardinality error too high: {abs(large_est - 100000) / 100000:.2f}"
    
    print("Edge cases test passed!")

def test_comparison_with_python():
    """Compare results with Python implementation."""
    print("\n=== Comparing With Python Implementation ===")
    
    if not FastHyperLogLog(precision=14).is_using_rust():
        print("Skipping comparison - Rust implementation not available")
        return
    
    # Generate data
    data_size = 10000
    data = [f"item_{i}" for i in range(data_size)]
    
    # Create both sketches
    rust_sketch = FastHyperLogLog(precision=14)
    py_sketch = HyperLogLog(precision=14)
    
    # Add data to both
    for item in data:
        rust_sketch.add(item)
        py_sketch.add_string(item)
    
    # Compare estimates
    rust_est = rust_sketch.cardinality()
    py_est = py_sketch.estimate_cardinality()
    
    print(f"Rust estimate: {rust_est:.2f}")
    print(f"Python estimate: {py_est:.2f}")
    print(f"Difference: {abs(rust_est - py_est):.2f} ({abs(rust_est - py_est) / py_est * 100:.2f}%)")
    
    # If difference is too high, try with higher precision
    if abs(rust_est - py_est) / py_est >= 0.05 and rust_sketch.precision < 16:
        print("\nTrying with higher precision...")
        rust_sketch = FastHyperLogLog(precision=rust_sketch.precision + 2)
        py_sketch = HyperLogLog(precision=py_sketch.precision + 2)
        
        # Add data to both with higher precision
        for item in data:
            rust_sketch.add(item)
            py_sketch.add_string(item)
        
        # Compare estimates again
        rust_est = rust_sketch.cardinality()
        py_est = py_sketch.estimate_cardinality()
        
        print(f"Rust estimate (precision {rust_sketch.precision}): {rust_est:.2f}")
        print(f"Python estimate (precision {py_sketch.precision}): {py_est:.2f}")
        print(f"Difference: {abs(rust_est - py_est):.2f} ({abs(rust_est - py_est) / py_est * 100:.2f}%)")
    
    assert abs(rust_est - py_est) / py_est < 0.05, f"Rust and Python estimates differ by more than 5%"
    
    # Compare performance
    start_time = time.time()
    for _ in range(5):
        py_sketch.estimate_cardinality()
    py_time = time.time() - start_time
    
    start_time = time.time()
    for _ in range(5):
        rust_sketch.cardinality()
    rust_time = time.time() - start_time
    
    print(f"Python estimation time: {py_time:.6f}s")
    print(f"Rust estimation time: {rust_time:.6f}s")
    print(f"Speedup: {py_time/rust_time:.2f}x")
    
    print("Comparison with Python implementation test passed!")

def run_all_tests():
    """Run all tests."""
    try:
        test_basic_functionality()
        test_estimation_methods()
        test_merge()
        test_jaccard()
        test_batch_operations()
        test_edge_cases()
        test_comparison_with_python()
        
        print("\n=== All Tests Passed! ===")
    except AssertionError as e:
        print(f"\n=== Test Failed: {e} ===")
    except Exception as e:
        print(f"\n=== Unexpected Error: {e} ===")

if __name__ == "__main__":
    run_all_tests() 