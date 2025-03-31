#!/usr/bin/env python3
"""
Example usage of FastHyperLogLog for cardinality estimation.
"""
from hammock.lib import FastHyperLogLog

def main():
    """Demonstrate basic usage of FastHyperLogLog."""
    # Create a new sketch with precision 12 (2^12 = 4096 registers)
    sketch = FastHyperLogLog(precision=12)
    
    # Check if using Rust implementation
    using_rust = sketch.is_using_rust()
    print(f"Using Rust implementation: {using_rust}")
    
    # Add some items to the sketch
    print("\nAdding items to sketch...")
    for i in range(10000):
        sketch.add(f"item_{i}")
    
    # Add some duplicates
    for i in range(2000):
        sketch.add(f"item_{i}")
    
    # Get cardinality estimate with different methods
    print("\nCardinality estimates:")
    original_est = sketch.estimate_cardinality(method="original")
    improved_est = sketch.estimate_cardinality(method="ertl_improved")
    mle_est = sketch.estimate_cardinality(method="ertl_mle")
    
    print(f"Original method:   {original_est:.2f}")
    print(f"Improved method:   {improved_est:.2f}")
    print(f"MLE method:        {mle_est:.2f}")
    print(f"Truth:             10000")
    
    # Create another sketch for similarity calculation
    print("\nCreating second sketch with 50% overlap...")
    sketch2 = FastHyperLogLog(precision=12)
    
    # Add 5000 overlapping items
    for i in range(5000):
        sketch2.add(f"item_{i}")
    
    # Add 5000 distinct items
    for i in range(10000, 15000):
        sketch2.add(f"item_{i}")
    
    # Calculate Jaccard similarity
    jaccard = sketch.jaccard(sketch2)
    print(f"Jaccard similarity: {jaccard:.4f} (expected ~0.33)")
    
    # Demonstrate batch operations
    print("\nDemonstrating batch operations...")
    batch_sketch = FastHyperLogLog(precision=12)
    
    # Create batch data
    batch_data = [f"batch_item_{i}" for i in range(20000)]
    
    # Use batch add if available
    if using_rust:
        batch_sketch.add_batch(batch_data)
    else:
        for item in batch_data:
            batch_sketch.add(item)
    
    # Get cardinality
    batch_est = batch_sketch.cardinality()
    print(f"Batch cardinality: {batch_est:.2f} (expected ~20000)")
    
    # Demonstrate merging
    print("\nDemonstrating merge operation...")
    
    # Merge the two sketches
    sketch.merge(sketch2)
    
    # Get cardinality of merged sketch (should be around 15000)
    merged_est = sketch.cardinality()
    print(f"Merged cardinality: {merged_est:.2f} (expected ~15000)")

if __name__ == "__main__":
    main() 