from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Optional, Dict, List
import xxhash # type: ignore

class AbstractSketch(ABC):
    """Base class for all sketch types."""
    
    @abstractmethod
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        pass
    
    @abstractmethod
    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch.
        
        Args:
            strings: List of strings to add to the sketch
        """
        pass
    
    @abstractmethod
    def similarity_values(self, other: 'AbstractSketch') -> Dict[str, float]:
        """Estimate similarity with another sketch.
        
        Returns:
            Dictionary of similarity measures
        """
        pass
    
    # Hash functions - static methods for use by all sketch implementations
    @staticmethod
    def _hash64_int(x: int, seed: int = 0) -> int:
        """64-bit hash function for integers.
        
        Args:
            x: Integer value to hash
            seed: Random seed for hashing
            
        Returns:
            64-bit hash value as integer
        """
        hasher = xxhash.xxh64(seed=seed)
        hasher.update(x.to_bytes(8, byteorder='little'))
        return hasher.intdigest()
    
    @staticmethod
    def _hash32_int(x: int, seed: int = 0) -> int:
        """32-bit hash function for integers.
        
        Args:
            x: Integer value to hash
            seed: Random seed for hashing
            
        Returns:
            32-bit hash value as integer
        """
        hasher = xxhash.xxh32(seed=seed)
        hasher.update(x.to_bytes(8, byteorder='little'))
        return hasher.intdigest()

    @staticmethod
    def _hash_str(s: bytes, seed: int = 0, hash_size: int = 32) -> int:
        """Hash a string using xxhash.
        
        Args:
            s: String to hash
            seed: Random seed for hashing
            hash_size: Size of hash in bits (32 or 64)
            
        Returns:
            Hash value as integer
        """
        if hash_size == 32:
            hasher = xxhash.xxh32(seed=seed)
        else:
            hasher = xxhash.xxh64(seed=seed)
        hasher.update(s)
        return hasher.intdigest()

    # Instance methods that use sketch's seed and hash_size (if available)
    def hash_str(self, s: bytes) -> int:
        """Instance method to hash a string using the instance's seed and hash_size.
        
        Args:
            s: String to hash
            
        Returns:
            Hash value as integer
        """
        seed = getattr(self, 'seed', 0)
        hash_size = getattr(self, 'hash_size', 32)
        return self._hash_str(s, seed=seed, hash_size=hash_size)

    def hash64_int(self, x: int) -> int:
        """Instance method to hash an integer using the instance's seed and 64-bit hash.
        
        Args:
            x: Integer value to hash
            
        Returns:
            64-bit hash value as integer
        """
        seed = getattr(self, 'seed', 0)
        return self._hash64_int(x, seed=seed)
    
    def hash32_int(self, x: int) -> int:
        """Instance method to hash an integer using the instance's seed and 32-bit hash.
        
        Args:
            x: Integer value to hash
            
        Returns:
            32-bit hash value as integer
        """
        seed = getattr(self, 'seed', 0)
        return self._hash32_int(x, seed=seed)

    # @abstractmethod
    # def merge(self, other: 'AbstractSketch') -> None:
    #     """Merge another sketch into this one."""
    #     pass
        
    # @abstractmethod
    # def write(self, filepath: str) -> None:
    #     """Write sketch to file."""
    #     pass
        
    # @classmethod
    # @abstractmethod
    # def load(cls, filepath: str) -> 'AbstractSketch':
    #     """Load sketch from file."""
    #     pass

# class AbstractDataSketch(ABC):
#     """Abstract base class for data-specific sketch wrappers."""
    
#     @abstractmethod
#     def __init__(self, sketch_type: str, **kwargs):
#         """Initialize with specified sketch type."""
#         pass
    
#     @classmethod
#     @abstractmethod
#     def from_file(cls, filename: str, sketch_type: str, **kwargs) -> Optional['AbstractDataSketch']:
#         """Create sketch from file."""
#         pass 