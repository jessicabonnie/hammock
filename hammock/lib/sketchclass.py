from typing import Literal, Optional, Union
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash
from hammock.lib.exact import ExactCounter

class Sketch:
    def __init__(self, 
                 sketch_type: Literal["hyperloglog", "minhash", "exact"],
                 precision: int = 14,
                 num_hashes: int = 128,
                 kmer_size: int = 0,
                 window_size: int = 0,
                 seed: int = 0,
                 debug: bool = False):
        """Initialize sketch with specified type.
        
        Args:
            sketch_type: "hyperloglog", "minhash", or "exact"
            precision: Precision for HyperLogLog
            num_hashes: Number of hashes for MinHash
            kmer_size: Size of k-mers
            window_size: Size of sliding window
            seed: Random seed
            debug: Whether to print debug information
        """
        if sketch_type == "hyperloglog":
            self.sketch = HyperLogLog(precision, kmer_size, window_size, seed, debug)
        elif sketch_type == "minhash":
            self.sketch = MinHash(num_hashes, kmer_size, window_size, seed, debug)
            precision = None  # Set precision to None for MinHash
        elif sketch_type == "exact":
            self.sketch = ExactCounter(None, kmer_size, window_size, seed, debug)
            precision = None
            num_hashes = None
        else:
            raise ValueError("Invalid sketch type. Use 'hyperloglog', 'minhash', or 'exact'")
        
        self.type = sketch_type
        self.total_interval_size = 0
        self.num_intervals = 0
        self.avg_interval_size = 0


    def add_interval_size(self, size: int) -> None:
        """Add an interval size to the running statistics.
        
        Args:
            size: Size of the interval being added
        """
        self.total_interval_size += size
        self.num_intervals += 1
        self.avg_interval_size = self.total_interval_size / self.num_intervals if self.num_intervals > 0 else 0

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
    
    def similarity_values(self, other: 'Sketch') -> tuple[float, float, float, float]:
        """Calculate similarity values between two sketches.
        """
        union_size = self.estimate_union(other)
        intersect = self.estimate_intersection(other)
        jaccard = self.estimate_jaccard(other)
        
        # Handle potential division by zero
        jacc_calc = intersect / union_size if union_size != 0 else 0
        
        return union_size, intersect, jaccard, jacc_calc

    def write_sketch(self, filepath: str) -> None:
        """Write sketch to file.
        
        The file format depends on the sketch type:
        - hyperloglog/minhash: Binary .npz format
        - exact: Text format
        
        Args:
            filepath: Path to output file
            
        Raises:
            NotImplementedError: This method is not yet implemented
        """
        # self.sketch.write_sketch(filepath)
        raise NotImplementedError("write_sketch() is not yet implemented")

    @classmethod
    def read_sketch(cls, filepath: str, sketch_type: Literal["hyperloglog", "minhash", "exact"]) -> 'Sketch':
        """Read sketch from file.
        
        Args:
            filepath: Path to input file
            sketch_type: Type of sketch to read ("hyperloglog", "minhash", or "exact")
            
        Returns:
            Sketch object loaded from file
            
        Raises:
            NotImplementedError: This method is not yet implemented
            ValueError: If sketch_type is invalid
        """
        # if sketch_type == "hyperloglog":
        #     inner_sketch = HyperLogLog.read_sketch(filepath)
        # elif sketch_type == "minhash":
        #     inner_sketch = MinHash.read_sketch(filepath)
        # elif sketch_type == "exact":
        #     inner_sketch = ExactCounter.read_sketch(filepath)
        # else:
        if sketch_type not in ["hyperloglog", "minhash", "exact"]:
            raise ValueError("Invalid sketch type. Use 'hyperloglog', 'minhash', or 'exact'")
        # sketch = cls(sketch_type=sketch_type)
        # sketch.sketch = inner_sketch
        # return sketch
        raise NotImplementedError("read_sketch() is not yet implemented")