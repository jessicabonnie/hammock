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

class TimeoutError(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutError("Function took too long to complete")

def run_with_timeout(func, timeout=5):
    """Run a function with a timeout."""
    # Set the timeout handler
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(timeout)
    
    try:
        result = func()
        signal.alarm(0)  # Disable the alarm
        return result
    except TimeoutError as e:
        print(f"TIMEOUT: {e}")
        return None

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
    print("Using Rust HyperLogLog implementation")
    sketch = RustHLL(precision=12, use_threading=True, min_thread_batch=10000)
    
    # Print available methods
    methods = dir(sketch)
    print(f"Available methods: {methods}")
    
    # Report if Rust is available
    print(f"Rust HyperLogLog available: {RUST_AVAILABLE}")
    print(f"Using Rust implementation: {True}")
    
    # Add some values
    print("\nAdding values...")
    values = [f"item_{i}" for i in range(1000)]
    
    try:
        sketch.add_batch(values)
        print("Successfully added batch of values")
    except Exception as e:
        print(f"Error adding batch: {e}")
        print("Falling back to individual adds")
        for value in values[:10]:  # Just add a few values for testing
            try:
                sketch.add_value(value)
                print(f"Successfully added value: {value}")
            except Exception as e:
                print(f"Error adding value {value}: {e}")
    
    # Get cardinality estimate with timeout
    print("\nGetting cardinality estimate...")
    def get_cardinality():
        return sketch.estimate()
    
    est = run_with_timeout(get_cardinality)
    if est is not None:
        print(f"Cardinality estimate: {est:.2f}")
    
    # Test optimization
    print("\nTesting optimization...")
    try:
        sketch.optimize()
        print("Successfully optimized sketch")
    except Exception as e:
        print(f"Error optimizing sketch: {e}")
    
    # Test debug info
    print("\nGetting debug info...")
    try:
        debug_info = sketch.debug_info()
        print(f"Debug info: {debug_info}")
    except Exception as e:
        print(f"Error getting debug info: {e}")
    
    # Test merging
    print("\nTesting merging...")
    try:
        # Create another sketch with different values
        other_sketch = RustHLL(precision=12)
        other_values = [f"other_item_{i}" for i in range(500, 1500)]
        other_sketch.add_batch(other_values)
        
        # Merge the sketches
        sketch.merge(other_sketch)
        print("Successfully merged sketches")
        
        # Get cardinality estimate of merged sketch
        merged_est = sketch.estimate()
        print(f"Merged cardinality estimate: {merged_est:.2f}")
    except Exception as e:
        print(f"Error during merge test: {e}")
    
    # Test Jaccard similarity
    print("\nTesting Jaccard similarity...")
    try:
        # Create two sketches with overlapping values
        sketch1 = RustHLL(precision=12)
        sketch2 = RustHLL(precision=12)
        
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
    except Exception as e:
        print(f"Error during Jaccard test: {e}")
    
    # Benchmark
    print("\nBenchmarking batch add...")
    large_values = [f"large_item_{i}" for i in range(100000)]
    
    start_time = time.time()
    sketch = RustHLL(precision=12, use_threading=True)
    sketch.add_batch(large_values)
    threaded_time = time.time() - start_time
    print(f"Threaded batch add time: {threaded_time:.4f} seconds")
    
    start_time = time.time()
    sketch = RustHLL(precision=12, use_threading=False)
    sketch.add_batch(large_values)
    non_threaded_time = time.time() - start_time
    print(f"Non-threaded batch add time: {non_threaded_time:.4f} seconds")
    print(f"Speedup: {non_threaded_time / threaded_time:.2f}x")

if __name__ == "__main__":
    main() 