"""
Fast HyperLogLog implementation using Rust.
"""
from typing import Optional, Dict, Any, List, Set, Tuple, Union
import sys
import os

# Try to import the Rust extension
try:
    # Add the rust_hll directory to the Python path
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "rust_hll"))
    from rust_hll import RustHLL
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False

# If the Rust extension is not available, use the Python implementation
if not RUST_AVAILABLE:
    from .hyperloglog import HyperLogLog as PyHyperLogLog
    print("Rust HyperLogLog extension not found. Using Python implementation.")


class FastHyperLogLog:
    """
    A fast HyperLogLog implementation that uses Rust if available,
    otherwise falls back to the Python implementation.
    """
    
    def __init__(self, precision: int = 12):
        """
        Initialize a new HyperLogLog sketch.
        
        Args:
            precision: The precision parameter (number of registers = 2^precision)
        """
        self.precision = precision
        
        if RUST_AVAILABLE:
            self._sketch = RustHLL(precision)
            self._using_rust = True
        else:
            self._sketch = PyHyperLogLog(p=precision)
            self._using_rust = False
    
    def add(self, value: str) -> None:
        """
        Add a value to the sketch.
        
        Args:
            value: The string value to add
        """
        self._sketch.add(value)
    
    def add_batch(self, values: List[str]) -> None:
        """
        Add multiple values to the sketch.
        
        Args:
            values: List of string values to add
        """
        if self._using_rust:
            self._sketch.add_batch(values)
        else:
            for value in values:
                self._sketch.add(value)
    
    def cardinality(self) -> float:
        """
        Get the estimated cardinality (number of unique values).
        
        Returns:
            Estimated number of unique values
        """
        if self._using_rust:
            return self._sketch.estimate()
        else:
            return self._sketch.cardinality()
    
    def merge(self, other: 'FastHyperLogLog') -> None:
        """
        Merge another HyperLogLog sketch into this one.
        
        Args:
            other: Another HyperLogLog sketch to merge
        """
        if self.precision != other.precision:
            raise ValueError("Cannot merge sketches with different precisions")
        
        if self._using_rust and other._using_rust:
            self._sketch.merge(other._sketch)
        elif not self._using_rust and not other._using_rust:
            self._sketch.merge(other._sketch)
        else:
            # Different implementations, need to convert
            raise NotImplementedError("Merging different HyperLogLog implementations not supported yet")
    
    def jaccard(self, other: 'FastHyperLogLog') -> float:
        """
        Calculate the Jaccard similarity with another sketch.
        
        Args:
            other: Another HyperLogLog sketch to compare with
            
        Returns:
            Jaccard similarity [0, 1]
        """
        if self.precision != other.precision:
            raise ValueError("Cannot compare sketches with different precisions")
        
        if self._using_rust and other._using_rust:
            return self._sketch.jaccard(other._sketch)
        else:
            # For Python implementation, calculate Jaccard similarity manually
            clone = FastHyperLogLog(self.precision)
            if self._using_rust:
                # Convert to Python implementation
                raise NotImplementedError("Jaccard similarity between different implementations not supported yet")
            else:
                # Both are Python
                union_card = self._sketch.cardinality() + other._sketch.cardinality()
                
                # Create a copy for merging
                merged = PyHyperLogLog(p=self.precision)
                merged.merge(self._sketch)
                merged.merge(other._sketch)
                
                union_est = merged.cardinality()
                if union_est == 0:
                    return 1.0  # Both are empty, consider them identical
                
                # Estimate intersection using inclusion-exclusion principle
                intersection_est = union_card - union_est
                return max(0.0, min(1.0, intersection_est / union_est))
    
    def is_using_rust(self) -> bool:
        """
        Check if this instance is using the Rust implementation.
        
        Returns:
            True if using Rust, False otherwise
        """
        return self._using_rust
    
    def __str__(self) -> str:
        """Get string representation."""
        impl = "Rust" if self._using_rust else "Python"
        return f"FastHyperLogLog(precision={self.precision}, implementation={impl})"
    
    def __repr__(self) -> str:
        """Get representation string."""
        return self.__str__() 