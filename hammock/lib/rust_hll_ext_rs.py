from __future__ import annotations
from typing import Optional, List, Tuple, TextIO, Union, Dict
import numpy as np  # type: ignore
import xxhash  # type: ignore

class RustHyperLogLog:
    """Pure Rust implementation of HyperLogLog."""
    
    def __init__(self, precision: int = 12):
        """Initialize the sketch."""
        if not 4 <= precision <= 16:
            raise ValueError("Precision must be between 4 and 16")
        
        self.precision = precision
        self.num_registers = 1 << precision
        self.registers = np.zeros(self.num_registers, dtype=np.float64)
        
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        self.add(s)
        
    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch."""
        for s in strings:
            self.add(s)
        
    def add(self, value: str) -> None:
        """Add a value to the sketch."""
        hash_val = xxhash.xxh64(value.encode()).intdigest()
        idx = hash_val & (self.num_registers - 1)
        rank = self._rho(hash_val)
        self.registers[idx] = max(self.registers[idx], rank)
        
    def _rho(self, value: int) -> int:
        """Calculate the position of the leftmost 1-bit."""
        if value == 0:
            return 64
        return 1 + (value & -value).bit_length() - 1
        
    def estimate_cardinality(self) -> float:
        """Estimate the cardinality of the set."""
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
        return np.all(self.registers == 0)
        
    def merge(self, other: 'RustHyperLogLog') -> None:
        """Merge another sketch into this one."""
        if not isinstance(other, RustHyperLogLog):
            raise ValueError("Can only merge with another RustHyperLogLog")
        if self.precision != other.precision:
            raise ValueError("Cannot merge HLLs with different precision")
        
        self.registers = np.maximum(self.registers, other.registers)
        
    def write(self, filepath: str) -> None:
        """Write sketch to file."""
        np.save(filepath, self.registers)
        
    @classmethod
    def load(cls, filepath: str) -> 'RustHyperLogLog':
        """Load sketch from file."""
        registers = np.load(filepath)
        sketch = cls(precision=int(np.log2(len(registers))))
        sketch.registers = registers
        return sketch 