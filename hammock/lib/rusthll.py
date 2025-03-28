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
    # Also add the target/release directory
    rust_target_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "rust_hll", "target", "release")
    sys.path.append(rust_target_dir)
    
    # First import the module to check available classes
    import rust_hll
    
    # Check what classes are available
    rust_classes = [name for name in dir(rust_hll) if not name.startswith('_')]
    print(f"Available classes in rust_hll: {rust_classes}")
    
    # Try to find a suitable HyperLogLog class
    RustHLLClass = None
    for class_name in rust_classes:
        try:
            class_obj = getattr(rust_hll, class_name)
            # Check if this class has the methods we need
            if (hasattr(class_obj, 'new') or hasattr(class_obj, '__new__') or 
                hasattr(class_obj, '__init__')):
                # This looks like a class
                if hasattr(class_obj, 'add') or hasattr(class_obj, 'estimate'):
                    # This looks like a HyperLogLog class
                    RustHLLClass = class_obj
                    break
        except Exception:
            pass
    
    if RustHLLClass:
        RUST_AVAILABLE = True
        print(f"Found Rust HyperLogLog class: {RustHLLClass.__name__}")
    else:
        print("No suitable HyperLogLog class found in rust_hll module")
        RUST_AVAILABLE = False
    
except ImportError as e:
    print(f"Error importing Rust HyperLogLog module: {e}")
    RUST_AVAILABLE = False
    RustHLLClass = None

# If the Rust extension is not available, use the Python implementation
if not RUST_AVAILABLE:
    from .hyperloglog import HyperLogLog as PyHyperLogLog
    print("Rust HyperLogLog extension not found. Using Python implementation.")


class FastHyperLogLog:
    """
    A fast HyperLogLog implementation that uses Rust if available,
    otherwise falls back to the Python implementation.
    """
    
    def __init__(self, precision: int = 12, kmer_size: int = 0, window_size: int = 0, 
                 seed: Optional[int] = None, debug: bool = False):
        """
        Initialize a new HyperLogLog sketch.
        
        Args:
            precision: The precision parameter (number of registers = 2^precision)
            kmer_size: Size of k-mers (0 for whole string mode)
            window_size: Size of sliding window (0 or == kmer_size for no windowing)
            seed: Random seed for hashing
            debug: Whether to print debug information
        """
        self.precision = precision
        self.kmer_size = kmer_size
        self.window_size = window_size
        self.seed = seed
        self.debug = debug
        
        if RUST_AVAILABLE:
            try:
                self._sketch = RustHLLClass(precision)
                self._using_rust = True
                if debug:
                    print(f"Using Rust HyperLogLog implementation")
                    print(f"Available methods: {dir(self._sketch)}")
            except Exception as e:
                print(f"Error initializing Rust HyperLogLog: {e}, falling back to Python implementation")
                self._using_rust = False
                self._sketch = PyHyperLogLog(
                    precision=precision,
                    kmer_size=kmer_size,
                    window_size=window_size,
                    seed=seed,
                    debug=debug
                )
        else:
            self._sketch = PyHyperLogLog(
                precision=precision,
                kmer_size=kmer_size,
                window_size=window_size,
                seed=seed,
                debug=debug
            )
            self._using_rust = False
            if debug:
                print(f"Using Python HyperLogLog implementation")
    
    def add(self, value):
        """Add a value to the sketch."""
        if self._using_rust:
            # Ensure value is a string for rust_hll
            value_str = str(value)
            self._sketch.add_value(value_str)
        else:
            self._sketch.add(value)
    
    def add_batch(self, values: List[str]) -> None:
        """
        Add a batch of values to the sketch.
        
        Args:
            values: The values to add
        """
        if self._using_rust:
            # Ensure all values are strings for rust_hll
            string_values = [str(v) for v in values]
            self._sketch.add_batch(string_values)
        else:
            for value in values:
                self._sketch.add_string(value)
    
    def add_int(self, value: int) -> None:
        """
        Add an integer value to the sketch.
        
        Args:
            value: Integer value to add
        """
        if self._using_rust:
            self._sketch.add(str(value))
        else:
            self._sketch.add_int(value)
    
    def add_string(self, s: str) -> None:
        """
        Add a string value to the sketch.
        
        Args:
            s: The string value to add
        """
        if self._using_rust:
            self._sketch.add_value(s)
        else:
            self._sketch.add_string(s)
    
    def cardinality(self) -> float:
        """
        Get the estimated cardinality (number of unique values).
        
        Returns:
            Estimated number of unique values
        """
        if self._using_rust:
            return self._sketch.estimate()
        else:
            return self._sketch.estimate_cardinality()
    
    def estimate_cardinality(self, method="ertl_mle"):
        """
        Estimate the cardinality of the multiset.
        
        Args:
            method: Estimation method ('original', 'ertl_improved', or 'ertl_mle')
                    Rust implementation supports all three methods
        
        Returns:
            Estimated cardinality
        """
        if self._using_rust:
            try:
                # Call the estimate method directly for the Rust implementation
                return self._sketch.estimate()
            except ValueError as e:
                # Handle unknown method errors
                print(f"Warning: {str(e)}, falling back to default method")
                return self._sketch.estimate()
        else:
            # For the Python fallback, call estimate_cardinality directly
            return self._sketch.estimate_cardinality(method)
    
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
            if self._using_rust:
                # Convert to Python implementation
                raise NotImplementedError("Jaccard similarity between different implementations not supported yet")
            else:
                # Use the Python implementation's method
                return self._sketch.estimate_jaccard(other._sketch)
    
    def estimate_jaccard(self, other: 'FastHyperLogLog') -> float:
        """
        Alias for jaccard() to maintain compatibility with existing code.
        
        Args:
            other: Another HyperLogLog sketch to compare with
            
        Returns:
            Jaccard similarity [0, 1]
        """
        return self.jaccard(other)
    
    def similarity_values(self, other: 'FastHyperLogLog') -> Dict[str, float]:
        """
        Get similarity values between this sketch and another.
        
        Args:
            other: Another HyperLogLog sketch to compare with
            
        Returns:
            Dictionary of similarity measures
        """
        if self._using_rust:
            return {"jaccard": self.jaccard(other)}
        else:
            # Use the Python implementation's method if available
            if hasattr(self._sketch, 'similarity_values'):
                return self._sketch.similarity_values(other._sketch)
            else:
                return {"jaccard": self.jaccard(other)}
    
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