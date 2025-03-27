#!/usr/bin/env python3
"""
Test the Rust HyperLogLog implementation.
"""
from hammock.lib.rusthll import FastHyperLogLog, RUST_AVAILABLE
import time
import signal

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
    # Create a sketch with debug enabled
    sketch = FastHyperLogLog(precision=12, debug=True)
    
    # Report if Rust is available
    print(f"Rust HyperLogLog available: {RUST_AVAILABLE}")
    print(f"Using Rust implementation: {sketch.is_using_rust()}")
    
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
                sketch.add(value)
                print(f"Successfully added value: {value}")
            except Exception as e:
                print(f"Error adding value {value}: {e}")
    
    # Get cardinality estimate with timeout
    print("\nGetting cardinality estimate...")
    def get_cardinality():
        return sketch.cardinality()
    
    est = run_with_timeout(get_cardinality)
    if est is not None:
        print(f"Cardinality estimate: {est:.2f}")
    
    # Test different estimation methods with timeout
    for method in ["original", "ertl_improved", "ertl_mle"]:
        print(f"\nTesting {method} method...")
        def get_estimate():
            return sketch.estimate_cardinality(method=method)
        
        est = run_with_timeout(get_estimate)
        if est is not None:
            print(f"{method} estimate: {est:.2f}")

if __name__ == "__main__":
    main() 