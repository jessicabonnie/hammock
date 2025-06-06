#!/usr/bin/env python3
"""
Compatibility wrapper for native RustHLL to match the API expected by benchmark scripts.
"""

try:
    import rust_hll
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False

# Import Python HLL implementation for fallback
try:
    from .hyperloglog import HyperLogLog
    PYTHON_HLL_AVAILABLE = True
except ImportError:
    PYTHON_HLL_AVAILABLE = False

class RustHLLWrapper:
    """Wrapper for native RustHLL to make it API-compatible with the benchmark scripts."""
    
    def __init__(self, precision=16, hash_size=32, debug=False, seed=None):
        """Initialize RustHLL with specified precision and hash size.
        
        Args:
            precision: Precision parameter for the HLL sketch
            hash_size: Size of the hash in bits (32 or 64)
            debug: Debug flag (not used, provided for backward compatibility)
            seed: Random seed for hashing
        """
        self._precision = precision
        self._hash_size = hash_size
        self._seed = seed if seed is not None else 42
        
        try:
            self._using_rust = True
            # Import rust_hll inside the method to ensure it's available in this scope
            import rust_hll
            # Create RustHLL with the correct parameters
            self._sketch = rust_hll.RustHLL(
                precision=precision,
                hash_size=hash_size,
                debug=debug,
                seed=self._seed
            )
        except (ImportError, ValueError) as e:
            self._using_rust = False
            if isinstance(e, ValueError) and "Precision must be between" in str(e):
                print(f"Rust implementation error: {e}")
                print(f"Falling back to Python implementation")
            else:
                print(f"Failed to load Rust implementation: {e}")
                print(f"Falling back to Python implementation")
                
            # Fall back to Python implementation
            from hammock.lib.hyperloglog import HyperLogLog
            self._sketch = HyperLogLog(hash_size=hash_size, precision=precision)
    
    def add(self, value):
        """Add a value to the sketch."""
        self._sketch.add_value(str(value))
    
    def add_string(self, value):
        """Add a string value to the sketch."""
        self._sketch.add_value(value)
    
    def add_batch(self, values):
        """Add a batch of values to the sketch."""
        # Convert all values to strings for consistent handling
        string_values = [str(value) for value in values]
        self._sketch.add_batch(string_values)
    
    def cardinality(self):
        """Estimate the cardinality."""
        return self._sketch.estimate_cardinality()
    
    def estimate_cardinality(self):
        """Estimate the cardinality of the set."""
        if self._using_rust:
            return self._sketch.estimate_cardinality()
        else:
            return self._sketch.estimate_cardinality()
    
    def merge(self, other):
        """Merge with another sketch."""
        if isinstance(other, RustHLLWrapper):
            self._sketch.merge(other._sketch)
        else:
            # Try to access the underlying native sketch if possible
            try:
                self._sketch.merge(other._sketch)
            except (AttributeError, TypeError):
                raise TypeError("Cannot merge with incompatible sketch type")
    
    def jaccard(self, other):
        """Calculate Jaccard similarity with another sketch."""
        if not hasattr(self._sketch, 'jaccard'):
            # If using Python implementation that doesn't support jaccard,
            # compute it from the cardinality estimates
            my_card = self.estimate_cardinality()
            other_card = other.estimate_cardinality()
            
            # Clone self using our clone method
            merged = self.clone()
            
            # Merge with other
            merged.merge(other)
            union_card = merged.estimate_cardinality()
            
            # Calculate Jaccard index
            if union_card == 0:
                return 1.0  # Both are empty
            
            # Jaccard index = size of intersection / size of union
            # Intersection size = sum of individual sizes - union size
            intersection_card = my_card + other_card - union_card
            # Ensure non-negative intersection (can happen due to estimation errors)
            intersection_card = max(0, intersection_card)
            
            return intersection_card / union_card
        
        # If jaccard is supported directly
        if isinstance(other, RustHLLWrapper):
            return self._sketch.jaccard(other._sketch)
        else:
            # Try to access the underlying native sketch if possible
            try:
                return self._sketch.jaccard(other._sketch)
            except (AttributeError, TypeError):
                raise TypeError("Cannot calculate Jaccard with incompatible sketch type")
    
    def __str__(self):
        """Return a string representation."""
        return f"RustHLLWrapper(precision={self._precision})"
    
    @property
    def precision(self):
        """Get the sketch precision."""
        return self._precision
    
    @property
    def hash_size(self):
        """Get the hash size."""
        return self._hash_size
    
    @property
    def using_rust(self):
        """Get whether the Rust implementation is being used."""
        return self._using_rust
    
    def clone(self):
        """Create a clone of this sketch."""
        result = RustHLLWrapper(precision=self._precision, hash_size=self._hash_size, debug=False)
        
        if self._using_rust:
            # If using Rust and it has clone method
            try:
                import rust_hll
                result._using_rust = True
                result._sketch = self._sketch.clone()
            except (AttributeError, ImportError):
                # If no clone method, create new and merge
                result._using_rust = True
                import rust_hll
                result._sketch = rust_hll.RustHLL(self._precision)
                result._sketch.merge(self._sketch)
        else:
            # For Python implementation
            from hammock.lib.hyperloglog import HyperLogLog
            result._using_rust = False
            result._sketch = HyperLogLog(hash_size=self._hash_size, precision=self._precision)
            result._sketch.registers = self._sketch.registers.copy()
            
        return result 