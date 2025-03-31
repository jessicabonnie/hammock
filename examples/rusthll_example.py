#!/usr/bin/env python3
"""
Example usage of the optimized RustHLL implementation.
This demonstrates the key features and performance characteristics.
"""

import time
import random
import sys
from typing import Set, List

# Try to import the rust_hll module
try:
    import rust_hll
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False
    print("Error: rust_hll module not found.")
    print("Please build the Rust extension first:")
    print("  cd rust_hll")
    print("  python -m maturin develop")
    sys.exit(1)

def generate_data(n: int, seed: int = 42) -> List[str]:
    """Generate n random data items with controlled uniqueness."""
    random.seed(seed)
    return [f"item_{random.randint(0, int(n*1.2))}" for _ in range(n)]

def main():
    print("RustHLL Example - Optimized HyperLogLog implementation")
    print("-" * 60)
    
    # Create a new RustHLL sketch with precision 12 (4096 registers)
    sketch = rust_hll.RustHLL(12)
    print(f"Created sketch: {sketch.debug_info()}")
    
    # Generate some test data
    n_items = 100_000
    print(f"\nGenerating {n_items:,} random items...")
    data = generate_data(n_items)
    
    # Calculate actual unique items
    unique_count = len(set(data))
    print(f"Actual unique items: {unique_count:,}")
    
    # Add data to the sketch
    print("\nAdding items to sketch...")
    start_time = time.time()
    sketch.add_batch(data)
    add_time = time.time() - start_time
    print(f"Added {n_items:,} items in {add_time:.4f} seconds")
    print(f"Processing rate: {n_items/add_time:,.0f} items/second")
    
    # Estimate cardinality
    print("\nEstimating cardinality...")
    start_time = time.time()
    estimate = sketch.estimate()
    est_time = time.time() - start_time
    error = abs(estimate - unique_count) / unique_count * 100
    print(f"Estimated unique items: {estimate:,.0f}")
    print(f"Actual unique items: {unique_count:,}")
    print(f"Error: {error:.2f}%")
    print(f"Estimation time: {est_time:.6f} seconds")
    
    # Create a second sketch with some overlap
    print("\nCreating second sketch with ~50% overlap...")
    sketch2 = rust_hll.RustHLL(12)
    data2 = generate_data(n_items, seed=84)  # Different seed gives partial overlap
    sketch2.add_batch(data2)
    
    # Calculate Jaccard similarity
    print("\nCalculating Jaccard similarity...")
    start_time = time.time()
    jaccard = sketch.jaccard(sketch2)
    jaccard_time = time.time() - start_time
    
    # Calculate actual Jaccard for comparison
    set1 = set(data)
    set2 = set(data2)
    actual_jaccard = len(set1.intersection(set2)) / len(set1.union(set2))
    
    print(f"Estimated Jaccard: {jaccard:.4f}")
    print(f"Actual Jaccard: {actual_jaccard:.4f}")
    print(f"Jaccard calculation time: {jaccard_time:.6f} seconds")
    
    # Test threading options with larger dataset
    print("\nTesting threading options with larger dataset...")
    big_data = generate_data(1_000_000, seed=100)
    
    # Threading enabled (default)
    sketch_threaded = rust_hll.RustHLL(12, use_threading=True, min_thread_batch=50000)
    start_time = time.time()
    sketch_threaded.add_batch(big_data)
    threaded_time = time.time() - start_time
    
    # Threading disabled
    sketch_single = rust_hll.RustHLL(12, use_threading=False)
    start_time = time.time()
    sketch_single.add_batch(big_data)
    single_time = time.time() - start_time
    
    print(f"With threading:    {threaded_time:.4f} seconds")
    print(f"Without threading: {single_time:.4f} seconds")
    print(f"Speedup factor: {single_time/threaded_time:.2f}x")
    
    # Debug info
    print("\nSketch debug info:")
    print(sketch.debug_info())
    
    print("\nExample completed successfully!")

if __name__ == "__main__":
    main() 