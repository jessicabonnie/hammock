from typing import Literal, Optional, Union
from hammock.lib.hyperloglog import HyperLogLog
from minhash import MinHash
from hammock.lib.exact import ExactCounter

class Sketch:
    def __init__(self, 
                 sketch_type: Literal["hyperloglog", "minhash", "exact"],
                 precision: int = 14,
                 num_hashes: int = 128,
                 kmer_size: int = 0,
                 window_size: int = 0,
                 seed: int = 0):
        """Initialize sketch with specified type.
        
        Args:
            sketch_type: "hyperloglog", "minhash", or "exact"
            precision: Precision for HyperLogLog
            num_hashes: Number of hashes for MinHash
            kmer_size: Size of k-mers
            window_size: Size of sliding window
            seed: Random seed
        """
        if sketch_type == "hyperloglog":
            self.sketch = HyperLogLog(precision, kmer_size, window_size, seed)
        elif sketch_type == "minhash":
            self.sketch = MinHash(num_hashes, kmer_size, window_size, seed)
            precision = None  # Set precision to None for MinHash
        elif sketch_type == "exact":
            self.sketch = ExactCounter(None, kmer_size, window_size, seed)
            precision = None
            num_hashes = None
        else:
            raise ValueError("Invalid sketch type. Use 'hyperloglog', 'minhash', or 'exact'")
        
        self.type = sketch_type

    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        self.sketch.add_string(s)

    def merge(self, other: 'Sketch') -> None:
        """Merge another sketch into this one."""
        if self.type != other.type:
            raise ValueError("Cannot merge different sketch types")
        self.sketch.merge(other.sketch)

    def estimate_similarity(self, other: 'Sketch') -> float:
        """Estimate similarity with another sketch.
        
        Returns Jaccard similarity for MinHash or
        intersection cardinality for HyperLogLog.
        """
        if self.type != other.type:
            raise ValueError("Cannot compare different sketch types")
            
        # if self.type == "minhash":
        return self.sketch.estimate_jaccard(other.sketch)
        # else:
        #     return self.sketch.estimate_intersection(other.sketch) 

    def estimate_cardinality(self) -> float:
        """Estimate cardinality of the sketch."""
        return self.sketch.estimate_cardinality()

    def estimate_intersection(self, other: 'Sketch') -> float:
        """Estimate intersection cardinality with another sketch."""
        if self.type != other.type:
            raise ValueError("Cannot compare different sketch types")
        return self.sketch.estimate_intersection(other.sketch)

    def estimate_union(self, other: 'Sketch') -> float:
        """Estimate union cardinality with another sketch."""
        if self.type != other.type:
            raise ValueError("Cannot compute union of different sketch types")
            
        return self.sketch.estimate_union(other.sketch)

    def estimate_jaccard(self, other: 'Sketch') -> float:
        """Estimate Jaccard similarity with another sketch."""
        if self.type != other.type:
            raise ValueError("Cannot compare different sketch types")
            
        return self.sketch.estimate_jaccard(other.sketch)