"""
Example usage of RustHyperLogLog for cardinality estimation.
"""
import numpy as np
from hammock.lib import RustHyperLogLog

def main():
    """Demonstrate basic usage of RustHyperLogLog."""
    
    # Create a sketch with precision 12 (default)
    sketch = RustHyperLogLog(precision=12)
    
    # Add some strings
    sketch.add("hello")
    sketch.add("world")
    sketch.add("hello")  # Duplicate, should be counted only once
    
    # Get cardinality estimate
    est = sketch.estimate_cardinality()
    print(f"Estimated cardinality: {est:.2f}")
    
    # Add more strings
    for i in range(1000):
        sketch.add(f"item{i}")
    
    # Get updated estimate
    est = sketch.estimate_cardinality()
    print(f"Estimated cardinality after adding 1000 items: {est:.2f}")
    
    # Create another sketch and merge
    sketch2 = RustHyperLogLog(precision=12)
    sketch2.add("world")
    sketch2.add("python")
    
    # Merge sketches
    sketch.merge(sketch2)
    
    # Get merged estimate
    est = sketch.estimate_cardinality()
    print(f"Estimated cardinality after merge: {est:.2f}")
    
    # Demonstrate batch operations
    batch_sketch = RustHyperLogLog(precision=12)
    strings = [f"batch{i}" for i in range(100)]
    batch_sketch.add_batch(strings)
    
    # Get batch estimate
    est = batch_sketch.estimate_cardinality()
    print(f"Estimated cardinality of batch: {est:.2f}")
    
    # Save and load sketch
    sketch.write("test_sketch.npy")
    loaded_sketch = RustHyperLogLog.load("test_sketch.npy")
    est = loaded_sketch.estimate_cardinality()
    print(f"Estimated cardinality of loaded sketch: {est:.2f}")

if __name__ == "__main__":
    main() 