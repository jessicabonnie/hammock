"""
This module provides the main interface for the Rust HyperLogLog implementation.

The module structure is:
1. rust_hll_ext.py (this file) - Main interface that provides a unified API
2. rust_hll_ext_rs.py - Contains the actual Rust implementation
3. rust_hll.py - Python wrapper that uses this interface

This file acts as a bridge between the Python and Rust implementations:
- If the Rust implementation is available (from rust_hll_ext_rs.py), it uses that
- If not, it falls back to a pure Python implementation
- The user doesn't need to know which implementation is being used

The fallback implementation is provided here rather than in rust_hll_ext_rs.py
to avoid circular imports and to keep the Rust-specific code separate.
"""

from __future__ import annotations
from typing import Optional, List, Tuple, TextIO, Union, Dict
import numpy as np  # type: ignore
import xxhash  # type: ignore
import warnings

# Try to import the Rust extension
try:
    from hammock.lib.rust_hll_ext_rs import RustHLLExtRS as RustHLLClass
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False

class RustHLLWrapper:
    """Rust implementation of HyperLogLog with fallback to pure Python."""
    
    def __init__(self, precision: int = 12, hash_size: int = 32):
        """Initialize the sketch.
        
        Args:
            precision: Number of bits to use for register addressing.
            hash_size: Size of hash in bits (32 or 64)
        """
        self.precision = precision
        self.num_registers = 1 << precision
        self.hash_size = hash_size
        
        if RUST_AVAILABLE:
            self._sketch = RustHLLClass(precision, hash_size=hash_size)
            self._using_rust = True
        else:
            self.registers = np.zeros(self.num_registers, dtype=np.float64)
            self._using_rust = False
        
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        self.add(s)
        
    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch."""
        if self._using_rust:
            self._sketch.add_batch(strings)
        else:
            for s in strings:
                self.add(s)
        
    def add(self, value: str) -> None:
        """Add a value to the sketch."""
        if self._using_rust:
            self._sketch.add_value(value)
        else:
            if self.hash_size == 32:
                hash_val = xxhash.xxh32(value.encode()).intdigest()
            else:
                hash_val = xxhash.xxh64(value.encode()).intdigest()
            idx = hash_val & (self.num_registers - 1)
            rank = self._rho(hash_val)
            self.registers[idx] = max(self.registers[idx], rank)
        
    def _rho(self, value: int) -> int:
        """Calculate the position of the leftmost 1-bit."""
        if value == 0:
            return self.hash_size
        return 1 + (value & -value).bit_length() - 1
        
    def estimate_cardinality(self) -> float:
        """Estimate the cardinality of the set."""
        if self._using_rust:
            return self._sketch.estimate_cardinality()
            
        # Calculate raw estimate
        raw_estimate = self._raw_estimate()
        
        # Apply bias correction for small cardinalities
        if raw_estimate <= 5 * self.num_registers:
            raw_estimate = self._bias_correct(raw_estimate)
            
        return raw_estimate
        
    def _raw_estimate(self) -> float:
        """Calculate raw estimate of cardinality."""
        # Calculate harmonic mean of 2^(-M[j])
        sum_inv = np.sum(2.0 ** -self.registers)
        alpha = self._get_alpha()
        return alpha * self.num_registers * self.num_registers / sum_inv
        
    def _bias_correct(self, raw_estimate: float) -> float:
        """Apply bias correction for small cardinalities."""
        # Use linear interpolation between known bias values
        if raw_estimate <= 0:
            return 0.0
            
        # Get bias correction table
        bias_table = self._get_bias_table()
        
        # Find the closest cardinality in the table
        closest_card = min(bias_table.keys(), key=lambda x: abs(x - raw_estimate))
        bias = bias_table[closest_card]
        
        return raw_estimate * (1 - bias)
        
    def _get_alpha(self) -> float:
        """Get alpha constant based on precision."""
        if self.precision == 4:
            return 0.673
        elif self.precision == 5:
            return 0.697
        elif self.precision == 6:
            return 0.709
        else:
            return 0.7213 / (1 + 1.079 / self.num_registers)
            
    def _get_bias_table(self) -> Dict[float, float]:
        """Get bias correction table based on precision."""
        if self.precision == 4:
            return {
                0: 0.673,
                1: 0.697,
                2: 0.709,
                3: 0.715,
                4: 0.718,
                5: 0.719
            }
        elif self.precision == 5:
            return {
                0: 0.697,
                1: 0.709,
                2: 0.715,
                3: 0.718,
                4: 0.719,
                5: 0.720
            }
        elif self.precision == 6:
            return {
                0: 0.709,
                1: 0.715,
                2: 0.718,
                3: 0.719,
                4: 0.720,
                5: 0.721
            }
        else:
            return {
                0: 0.7213,
                1: 0.7213,
                2: 0.7213,
                3: 0.7213,
                4: 0.7213,
                5: 0.7213
            }
            
    def is_empty(self) -> bool:
        """Check if the sketch is empty."""
        if self._using_rust:
            return self._sketch.is_empty()
        return np.all(self.registers == 0)
        
    def merge(self, other: 'RustHLLWrapper') -> None:
        """Merge another sketch into this one."""
        if not isinstance(other, RustHLLWrapper):
            raise ValueError("Can only merge with another RustHLLWrapper")
        if self.precision != other.precision:
            raise ValueError("Cannot merge HLLs with different precision")
        
        if self._using_rust and other._using_rust:
            self._sketch.merge(other._sketch)
        elif not self._using_rust and not other._using_rust:
            self.registers = np.maximum(self.registers, other.registers)
        else:
            raise ValueError("Cannot merge sketches with different implementations")
        
    def write(self, filepath: str) -> None:
        """Write sketch to file."""
        if self._using_rust:
            self._sketch.write(filepath)
        else:
            np.save(filepath, self.registers)
        
    @classmethod
    def load(cls, filepath: str) -> 'RustHLLWrapper':
        """Load sketch from file."""
        if RUST_AVAILABLE:
            sketch = cls()
            sketch._sketch = RustHLLClass.load(filepath)
            sketch._using_rust = True
            return sketch
        else:
            registers = np.load(filepath)
            sketch = cls(precision=int(np.log2(len(registers))))
            sketch.registers = registers
            sketch._using_rust = False
            return sketch 