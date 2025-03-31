#!/usr/bin/env python3
"""
Test the optimized Rust HyperLogLog implementation.
"""
import time
import signal
import sys

# Try to import the rust_hll module
try:
    import rust_hll
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False
    print("Rust HLL module not available. Please build it with 'cd rust_hll && python -m maturin develop'")
    sys.exit(1)

def main():
    # Print module information
    print(f"Available classes in rust_hll: {dir(rust_hll)}")
    
    # Create a RustHLL sketch
    try:
        RustHLL = getattr(rust_hll, 'RustHLL')
        print(f"Found Rust HyperLogLog class: {RustHLL.__name__}")
    except AttributeError:
        print("Could not find RustHLL class in the rust_hll module")
        sys.exit(1)
    
    # Create a new sketch
    print("Creating FastHyperLogLog sketch...")
    print("Using Rust HyperLogLog implementation")
    sketch = RustHLL(precision=12, use_threading=True, min_thread_batch=10000)
    
    # Print available methods
    methods = dir(sketch)
    print(f"Available methods: {methods}")
    print(f"Using Rust implementation: {True}")
    
    # Add some values
    print("\nAdding values...")
    start_time = time.time()
    values = [f"item_{i}" for i in range(10000)]
    sketch.add_batch(values)
    end_time = time.time()
    print(f"Added {len(values)} values in {end_time - start_time:.4f}s")
    
    # Get cardinality estimate
    print("\nTesting cardinality estimation:")
    est = sketch.estimate()
    print(f"Cardinality estimate: {est:.2f}")
    print(f"Actual cardinality: {len(set(values))}")
    print(f"Error: {abs(est - len(set(values))) / len(set(values)) * 100:.2f}%")
    
    # Test Jaccard similarity
    print("\nTesting Jaccard similarity...")
    # Create two sketches with overlapping values
    sketch1 = RustHLL(precision=14)
    sketch2 = RustHLL(precision=14)
    
    # Add values with 50% overlap
    values1 = [f"item_{i}" for i in range(1000)]
    values2 = [f"item_{i}" for i in range(500, 1500)]
    
    sketch1.add_batch(values1)
    sketch2.add_batch(values2)
    
    # Calculate Jaccard similarity
    jaccard = sketch1.jaccard(sketch2)
    print(f"Jaccard similarity: {jaccard:.4f}")
    print(f"Expected similarity: ~0.3333 (500 / 1500)")
    
    # Calculate actual Jaccard for comparison
    set1 = set(values1)
    set2 = set(values2)
    actual_jaccard = len(set1.intersection(set2)) / len(set1.union(set2))
    print(f"Actual Jaccard similarity: {actual_jaccard:.4f}")
    print(f"Error: {abs(jaccard - actual_jaccard):.4f}")
    
    # Benchmark
    print("\nBenchmarking batch add...")
    large_values = [f"large_item_{i}" for i in range(100000)]
    
    start_time = time.time()
    sketch = RustHLL(precision=14, use_threading=True)
    sketch.add_batch(large_values)
    threaded_time = time.time() - start_time
    print(f"Threaded batch add time: {threaded_time:.4f} seconds")
    
    start_time = time.time()
    sketch = RustHLL(precision=14, use_threading=False)
    sketch.add_batch(large_values)
    non_threaded_time = time.time() - start_time
    print(f"Non-threaded batch add time: {non_threaded_time:.4f} seconds")
    print(f"Speedup: {non_threaded_time / threaded_time:.2f}x")

if __name__ == "__main__":
    main() 