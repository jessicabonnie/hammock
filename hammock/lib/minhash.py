from __future__ import annotations
import numpy as np # type: ignore
import xxhash # type: ignore
from typing import Optional, List, Union, Dict
from hammock.lib.abstractsketch import AbstractSketch

class MinHash(AbstractSketch):
    def __init__(self, 
                 num_hashes: int = 128,
                 kmer_size: int = 0,
                 window_size: int = 0,
                 seed: int = 42,
                 debug: bool = False):
        """Initialize MinHash sketch.
        
        Args:
            num_hashes: Number of hash functions
            kmer_size: Size of k-mers (0 for whole string mode)
            window_size: Size of sliding window (0 or == kmer_size for no windowing)
            seed: Random seed for hashing
            debug: Whether to print debug information
        """
        super().__init__()
        self.num_hashes = num_hashes
        self.kmer_size = kmer_size
        self.window_size = window_size
        self.seed = seed
        self.debug = debug
        # Initialize with maximum possible value
        max_val = np.iinfo(np.uint64).max
        self.min_hashes = np.full(num_hashes, max_val, dtype=np.uint64)
        # self.hashers = [xxhash.xxh64(seed=seed + i * 31337) for i in range(0,num_hashes)]
        self.hashers=[None]*num_hashes
        if self.debug:
            print(f"\nDebug MinHash init:")
            print(f"Max uint64 value: {max_val}")
            print(f"Initial min_hashes type: {self.min_hashes.dtype}")
            print(f"Initial min_hashes first value: {self.min_hashes[0]}")
        
        # Create independent hash functions with different seeds
        # self.hashers = [xxhash.xxh64(seed=seed + i * 31337) for i in range(0,num_hashes)]
    
    def _hash_str(self, s: str) -> np.ndarray:
        """Generate all hash values for a string."""
        if not s:  # Guard against empty strings
            raise ValueError("Cannot hash empty string")
        
        hashes = np.zeros(self.num_hashes, dtype=np.uint64)
        s_bytes = s.encode()
        
        # Use a large prime number for better distribution
        prime = 2**61 - 1  # Mersenne prime
        
        for i in range(self.num_hashes):
            # Create and update hashers separately
            hasher1 = xxhash.xxh64(seed=self.seed + i)
            hasher2 = xxhash.xxh64(seed=self.seed + i + self.num_hashes)
            
            hasher1.update(s_bytes)
            hasher2.update(s_bytes)
            
            h1 = hasher1.intdigest()
            h2 = hasher2.intdigest()
            
            # Combine hashes using linear combination
            combined = (h1 + i * h2) % prime
            
            # Map to uint64 range
            hashes[i] = combined % (1 << 64)
        
        return hashes
    
    def _process_kmer(self, kmer: str) -> None:
        """Process a k-mer and update minimum hash values."""
        if not kmer:
            return
        
        # Get hash values for this k-mer
        try:
            hashes = self._hash_str(kmer)
            # Update min_hashes with element-wise minimum
            self.min_hashes = np.minimum(self.min_hashes, hashes)
            
            if self.debug and int(kmer) % 2000 == 0:
                print(f"\nProcessing k-mer: {kmer}")
                print(f"Hash values: {hashes[:5]}")
                print(f"Updated min_hashes: {self.min_hashes[:5]}")
        except Exception as e:
            if self.debug:
                print(f"Error processing k-mer {kmer}: {str(e)}")
            raise

    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        if not s:
            return
        
        if self.kmer_size == 0:
            # Whole string mode
            self._process_kmer(s)
        else:
            # K-mer mode - process all k-mers
            for i in range(len(s) - self.kmer_size + 1):
                kmer = s[i:i + self.kmer_size]
                self._process_kmer(kmer)

    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch.
        
        Args:
            strings: List of strings to add to the sketch
        """
        for s in strings:
            self.add_string(s)

    def merge(self, other: 'MinHash') -> 'MinHash':
        """Merge another MinHash sketch into this one.
        
        Returns:
            A new MinHash sketch containing the merged values
        """
        if self.num_hashes != other.num_hashes:
            raise ValueError("Cannot merge MinHash sketches with different sizes")
        
        # Create a new sketch with same parameters
        merged = MinHash(
            num_hashes=self.num_hashes,
            kmer_size=self.kmer_size,
            window_size=self.window_size,
            seed=self.seed,
            debug=self.debug
        )
        
        # Set its min_hashes to the minimum of both sketches
        merged.min_hashes = np.minimum(self.min_hashes, other.min_hashes)
        
        return merged

    def estimate_cardinality(self) -> float:
        """Estimate cardinality of the set using MinHash signatures."""
        mask = ~np.isinf(self.min_hashes)
        if not np.any(mask):
            return 0.0
            
        k = max(1, int(0.7 * self.num_hashes))
        kth_min = np.partition(self.min_hashes[mask], k)[k]
        
        return (k / kth_min) * (2**64) if kth_min > 0 else 0.0

    def estimate_intersection(self, other: 'MinHash') -> float:
        """Estimate intersection cardinality with another MinHash sketch."""
        if self.num_hashes != other.num_hashes:
            raise ValueError("Cannot compute intersection of MinHash sketches with different sizes")
        
        jaccard = self.estimate_jaccard(other)
        union = self.estimate_union(other)
        
        return max(0.0, jaccard * union)

    def estimate_jaccard(self, other: 'MinHash') -> float:
        """Estimate Jaccard similarity with another MinHash sketch."""
        if not isinstance(other, MinHash):
            raise TypeError("Can only compare with another MinHash sketch")
        if self.num_hashes != other.num_hashes:
            raise ValueError("Cannot compare MinHash sketches with different numbers of hashes")
        if self.seed != other.seed:
            raise ValueError("Cannot compare MinHash sketches with different seeds")
        
        max_val = np.iinfo(np.uint64).max
        # Return 0 if either sketch is empty
        if np.all(self.min_hashes == max_val) or np.all(other.min_hashes == max_val):
            return 0.0
        
        # Count matches using element-wise comparison
        matches = np.sum(self.min_hashes == other.min_hashes)
        
        if self.debug:
            print(f"\nJaccard estimation:")
            print(f"Number of matching hashes: {matches}")
            print(f"Total number of hashes: {self.num_hashes}")
            print(f"First few hashes sketch1: {self.min_hashes[:5]}")
            print(f"First few hashes sketch2: {other.min_hashes[:5]}")
        
        # Jaccard similarity is the fraction of matching minimums
        return float(matches) / self.num_hashes

    def estimate_union(self, other: 'MinHash') -> float:
        """Estimate union cardinality with another MinHash sketch."""
        if self.num_hashes != other.num_hashes:
            raise ValueError("Cannot compute union of MinHash sketches with different sizes")
        
        combined = np.minimum(self.min_hashes, other.min_hashes)
        mask = ~np.isinf(combined)
        if not np.any(mask):
            return 0.0
            
        k = max(1, int(0.7 * self.num_hashes))
        kth_min = np.partition(combined[mask], k)[k]
        
        return (k / kth_min) * (2**64) if kth_min > 0 else 0.0

    def write(self, filepath: str) -> None:
        """Write the MinHash sketch to a file.
        
        Args:
            filepath: Path to write the sketch to
        """
        import pickle
        with open(filepath, 'wb') as f:
            pickle.dump({
                'num_hashes': self.num_hashes,
                'kmer_size': self.kmer_size,
                'window_size': self.window_size,
                'seed': self.seed,
                'debug': self.debug,
                'min_hashes': self.min_hashes
            }, f)

    @classmethod
    def load(cls, filepath: str) -> 'MinHash':
        """Load a MinHash sketch from a file.
        
        Args:
            filepath: Path to load the sketch from
            
        Returns:
            The loaded MinHash sketch
        """
        import pickle
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
            
        sketch = cls(
            num_hashes=data['num_hashes'],
            kmer_size=data['kmer_size'],
            window_size=data['window_size'],
            seed=data['seed'],
            debug=data['debug']
        )
        
        sketch.min_hashes = data['min_hashes']
        return sketch

    def add_int(self, value: int) -> None:
        """Add an integer to the sketch."""
        # Convert integer to string to ensure consistent hashing
        self._process_kmer(str(value))

    def add_int_old(self, value: int) -> None:
        hashes = np.array([h.intdigest() for h in self.hashers])
        self.min_hashes = np.minimum(self.min_hashes, hashes)

    def similarity_values(self, other: 'AbstractSketch') -> Dict[str, float]:
        """Estimate similarity with another sketch.
        
        Returns:
            Dictionary containing similarity metrics
        """
        if not isinstance(other, MinHash):
            raise TypeError("Can only compare MinHash sketches")
        if self.num_hashes != other.num_hashes:
            raise ValueError("Cannot compare MinHash sketches with different numbers of hashes")
        if self.seed != other.seed:
            raise ValueError("Cannot compare MinHash sketches with different seeds")
        
        max_val = np.iinfo(np.uint64).max
        # Return 0 if either sketch is empty
        if np.all(self.min_hashes == max_val) or np.all(other.min_hashes == max_val):
            return {
                'jaccard': 0.0,
                'containment': 0.0,
                'intersection': 0.0,
                'union': 0.0
            }
        
        # Count matches using element-wise comparison
        matches = np.sum(self.min_hashes == other.min_hashes)
        jaccard = float(matches) / self.num_hashes
        
        return {
            'jaccard': jaccard,
            'containment': jaccard * self.estimate_union(other) / self.estimate_cardinality(),
            'intersection': jaccard * self.estimate_union(other),
            'union': self.estimate_union(other)
        }
