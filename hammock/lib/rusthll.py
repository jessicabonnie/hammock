"""
Fast HyperLogLog implementation using Rust.
"""
from __future__ import annotations
from typing import Optional, Dict, Any, List, Set, Tuple, Union
import sys
import os

from hammock.lib.abstractsketch import AbstractSketch
import numpy as np  # type: ignore
import xxhash  # type: ignore

# Try to import the Rust extension
try:
    # Add the rust_hll directory to the Python path
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "rust_hll"))
    # Also add the target/release directory
    rust_target_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "rust_hll", "target", "release")
    sys.path.append(rust_target_dir)
    
    # Import the Rust module
    import rust_hll
    RUST_AVAILABLE = True
    print(f"Found Rust HyperLogLog module")
except ImportError as e:
    print(f"Error importing Rust HyperLogLog module: {e}")
    RUST_AVAILABLE = False
    rust_hll = None

# If the Rust extension is not available, use the Python implementation
if not RUST_AVAILABLE:
    from .hyperloglog import HyperLogLog as PyHyperLogLog
    print("Rust HyperLogLog extension not found. Using Python implementation.")


class RustHyperLogLog(AbstractSketch):
    """Rust implementation of HyperLogLog."""
    
    def __init__(self, precision: int = 12, debug: bool = False, expected_cardinality: Optional[int] = None):
        """Initialize the sketch.
        
        Args:
            precision: Base precision value (4-24)
            debug: Enable debug output
            expected_cardinality: Expected number of unique items. If provided, precision will be adjusted.
        """
        if not RUST_AVAILABLE:
            raise ImportError("Rust HyperLogLog module is not available")
            
        if expected_cardinality is not None:
            # Adjust precision based on expected cardinality
            if expected_cardinality < 50:
                precision = 4 
            elif expected_cardinality < 100:
                precision = 6  # For smallest sets, use lower precision
            elif expected_cardinality < 1000:
                precision = 8  # For small sets, use lower precision
            elif expected_cardinality < 10000:
                precision = 12  # For medium sets
            elif expected_cardinality < 100000:
                precision = 18  # For larger sets
            else:
                precision = 22  # For very large sets
        
        if not 4 <= precision <= 24:
            raise ValueError("Precision must be between 4 and 24")
        
        self.precision = precision
        self.num_registers = 1 << precision
        self.debug = debug
        self.sketch = rust_hll.RustHLL(precision)
        
    def is_using_rust(self) -> bool:
        """Check if using Rust implementation."""
        return True  # This class always uses Rust
        
    def add(self, value: Union[str, int, float]) -> None:
        """Add a value to the sketch."""
        # Convert to string and use add for Rust implementation
        if isinstance(value, (int, float)):
            value = str(value)
        self.sketch.add_value(value)
        
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        self.sketch.add_value(s)
        
    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch."""
        # Convert any non-string values to strings
        string_values = [str(s) if not isinstance(s, str) else s for s in strings]
        self.sketch.add_batch(string_values)
        
    def estimate_cardinality(self, method: str = "fast_mle") -> float:
        """Estimate the cardinality of the set.
        
        Args:
            method: The estimation method to use. One of:
                - "original": Original HyperLogLog estimation
                - "ertl_mle": ERTL Maximum Likelihood Estimation
                - "fast_mle": Fast MLE implementation (default)
        """
        if method not in ["original", "ertl_mle", "fast_mle"]:
            raise ValueError("Method must be one of: original, ertl_mle, fast_mle")
        return float(self.sketch.estimate_cardinality())
        
    def estimate_intersection(self, other: 'RustHyperLogLog') -> float:
        """Estimate intersection cardinality with another sketch."""
        if not isinstance(other, RustHyperLogLog):
            raise ValueError("Can only compare with another RustHyperLogLog")
        if self.precision != other.precision:
            raise ValueError("Cannot compute intersection of HLLs with different precision")
        
        # Basic inclusion-exclusion
        a = self.estimate_cardinality()
        b = other.estimate_cardinality()
        union = self.estimate_union(other)
        intersection = max(0.0, a + b - union)
        
        if self.debug:
            print(f"DEBUG: a={a:.1f}, b={b:.1f}, union={union:.1f}, intersection={intersection:.1f}")
        return float(intersection)
        
    def estimate_union(self, other: 'RustHyperLogLog') -> float:
        """Estimate union cardinality with another sketch."""
        if not isinstance(other, RustHyperLogLog):
            raise ValueError("Can only compare with another RustHyperLogLog")
        if self.precision != other.precision:
            raise ValueError("Cannot compute union of HLLs with different precision")
        
        # Create a temporary sketch for the union
        temp_sketch = rust_hll.RustHLL(self.precision)
        temp_sketch.merge(self.sketch)
        temp_sketch.merge(other.sketch)
        return float(temp_sketch.estimate_cardinality())
        
    def estimate_jaccard(self, other: 'RustHyperLogLog') -> float:
        """Estimate Jaccard similarity with another sketch."""
        if not isinstance(other, RustHyperLogLog):
            raise ValueError("Can only compare with another RustHyperLogLog")
        if self.precision != other.precision:
            raise ValueError("Cannot compare HLLs with different precision")
        
        if self.is_empty() or other.is_empty():
            return 0.0
        
        intersection = self.estimate_intersection(other)
        union = self.estimate_union(other)
        
        if union == 0:
            return 0.0
            
        return float(intersection / union)
        
    def is_empty(self) -> bool:
        """Check if the sketch is empty."""
        return self.sketch.is_empty()
        
    def merge(self, other: 'RustHyperLogLog') -> None:
        """Merge another sketch into this one."""
        if not isinstance(other, RustHyperLogLog):
            raise ValueError("Can only merge with another RustHyperLogLog")
        if self.precision != other.precision:
            raise ValueError("Cannot merge HLLs with different precision")
        
        self.sketch.merge(other.sketch)
        
    def save(self, filepath: str) -> None:
        """Write sketch to file."""
        # Convert to absolute path and ensure directory exists
        filepath = os.path.abspath(filepath)
        dir_path = os.path.dirname(filepath)
        if dir_path:
            os.makedirs(dir_path, exist_ok=True)
        try:
            self.sketch.write(filepath)
        except Exception as e:
            print(f"Error saving sketch to {filepath}: {e}", file=sys.stderr)
            raise
        
    @classmethod
    def load(cls, filepath: str) -> 'RustHyperLogLog':
        """Load sketch from file."""
        # Convert to absolute path
        filepath = os.path.abspath(filepath)
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Sketch file not found: {filepath}")
        try:
            sketch = cls()
            sketch.sketch = rust_hll.RustHLL.load(filepath)
            return sketch
        except Exception as e:
            print(f"Error loading sketch from {filepath}: {e}", file=sys.stderr)
            raise
    
    def cardinality(self) -> float:
        """Get the estimated cardinality (number of unique values)."""
        return float(self.sketch.estimate_cardinality())
    
    def jaccard(self, other: 'RustHyperLogLog') -> float:
        """Calculate the Jaccard similarity with another sketch."""
        return float(self.estimate_jaccard(other))
    
    def similarity_values(self, other: 'RustHyperLogLog') -> Dict[str, float]:
        """Get similarity values between this sketch and another."""
        return {"jaccard": float(self.jaccard(other))}
    
    def __str__(self) -> str:
        """Get string representation."""
        return f"RustHyperLogLog(precision={self.precision})"
    
    def __repr__(self) -> str:
        """Get representation string."""
        return self.__str__() 