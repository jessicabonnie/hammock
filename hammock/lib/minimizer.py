#!/usr/bin/env python
from Bio import SeqIO # type: ignore
from Digest import window_minimizer # type: ignore
from typing import Optional, Union, Literal, Dict
from hammock.lib.abstractsketch import AbstractSketch
from hammock.lib.hyperloglog import HyperLogLog

class MinimizerSketch(AbstractSketch):
    def __init__(self, k: int, w: int, gapk: int):
        if k <= 3:
            raise ValueError("k must be greater than 3")
        self.k = k
        self.w = w
        self.gapk = gapk
        # Use HyperLogLog directly for both sketches
        self.minimizer_sketch = HyperLogLog(kmer_size=0)
        self.gap_sketch = HyperLogLog(kmer_size=gapk)
    
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        if not s or len(s) < self.k:
            return
            
        # Get minimizers and their positions
        try:
            minimizers = window_minimizer(s, self.k, self.w, include_hash=True)
            
            # 1. Add hashes and end kmers to first sketch
            hashes = [h for _, h in minimizers]
            if len(hashes) >= 2:
                first_last = s[:self.k] + s[-self.k:]
                self.minimizer_sketch.add_string(first_last)
            for h in hashes:
                self.minimizer_sketch.add_string(str(h))
                
            # 2. Create and add gap pattern to second sketch
            if len(minimizers) >= 2:
                positions = [pos for pos, _ in minimizers]
                gaps = [str(positions[i+1] - positions[i]) for i in range(len(positions)-1)]
                gap_pattern = ",".join(gaps)
                self.gap_sketch.add_string(gap_pattern)
                
        except RuntimeError:
            return  # Handle case where minimizer generation fails
            
    def compare_overlaps(self, other: 'MinimizerSketch') -> tuple[float, float]:
        """Returns (hash_similarity, gap_similarity)"""
        if not isinstance(other, MinimizerSketch):
            raise ValueError("Can only compare with another MinimizerSketch")
        if self.k != other.k or self.w != other.w or self.gapk != other.gapk:
            raise ValueError("Cannot compare sketches with different parameters")
            
        # Check if either sketch is empty (no items added)
        if self.minimizer_sketch.item_count == 0 or other.minimizer_sketch.item_count == 0:
            return 0.0, 0.0
        
        hash_sim = self.minimizer_sketch.estimate_jaccard(other.minimizer_sketch)
        gap_sim = self.gap_sketch.estimate_jaccard(other.gap_sketch)
        return hash_sim, gap_sim
        
    def estimate_similarity(self, other: AbstractSketch) -> Dict[str, float]:
        """Estimate similarity with another sketch.
        
        Returns:
            Dictionary containing 'hash_similarity', 'gap_similarity', and 'combined_similarity'
        """
        if not isinstance(other, MinimizerSketch):
            raise ValueError("Can only compare with another MinimizerSketch")
        
        hash_sim, gap_sim = self.compare_overlaps(other)
        combined_sim = (hash_sim + gap_sim) / 2
        
        return {
            'hash_similarity': hash_sim,
            'gap_similarity': gap_sim,
            'combined_similarity': combined_sim
        }

# class MinimizerSketch(AbstractSketch):
#     """Sketch class for sequence data using minimizers with an underlying sketch type."""
    
#     def __init__(self, 
#                  window_size: int = 40, 
#                  kmer_size: int = 8,
#                  sketch_type: Literal["hyperloglog", "minhash"] = "hyperloglog",
#                  precision: int = 8,
#                  num_hashes: int = 128,
#                  seed: int = 0):
#         """Initialize MinimizerSketch.
        
#         Args:
#             window_size: Size of sliding window for minimizer selection
#             kmer_size: Size of kmers to use
#             sketch_type: Type of underlying sketch ("hyperloglog" or "minhash")
#             precision: Precision for HyperLogLog sketching
#             num_hashes: Number of hash functions for MinHash sketching
#             seed: Random seed for hash functions
#         """
#         self.window_size = window_size
#         self.kmer_size = kmer_size
#         self.seed = seed
        
#         # Initialize underlying sketch
#         if sketch_type == "hyperloglog":
#             self.sketch = HyperLogLog(precision=precision, seed=seed)
#         elif sketch_type == "minhash":
#             self.sketch = MinHash(num_hashes=num_hashes, seed=seed)
#         else:
#             raise ValueError(f"Invalid sketch type: {sketch_type}")
    
#     def add_string(self, s: str) -> None:
#         """Add a string to the sketch."""
#         minimizers = window_minimizer(s, 
#                                     w=self.window_size, 
#                                     k=self.kmer_size, 
#                                     include_hash=True)
#         for _, hash_val in minimizers:
#             self.sketch.add_int(hash_val)
            
#     def add_sequence(self, sequence: str) -> None:
#         """Add a sequence using the minimizer approach."""
#         self.add_string(sequence)
#         # Add concatenated flanking kmers from start and end of sequence
#         self.sketch.add_string(sequence[:self.kmer_size] + sequence[-self.kmer_size:])
#         # For each window:
#         #   1. Find minimizer
#         #   2. Store minimizer with its position
#         # ... implementation ...
    
#     def find_overlaps(self, other_sketch: 'MinimizerSketch') -> list:
#         """Find overlapping regions between two sketches"""
#         # Compare minimizers and their positions to identify overlaps
#         # Return list of overlap regions
#         # ... implementation ...
        
#     def estimate_cardinality(self) -> float:
#         """Estimate cardinality using underlying sketch."""
#         return self.sketch.estimate_cardinality()
        
#     def estimate_jaccard(self, other: 'AbstractSketch') -> float:
#         """Estimate Jaccard similarity with another sketch."""
#         if not isinstance(other, MinimizerSketch):
#             raise ValueError("Can only compare with another MinimizerSketch")
#         return self.sketch.estimate_jaccard(other.sketch)
        
#     def merge(self, other: 'AbstractSketch') -> None:
#         """Merge another sketch into this one."""
#         if not isinstance(other, MinimizerSketch):
#             raise ValueError("Can only merge with another MinimizerSketch")
#         self.sketch.merge(other.sketch)

#     def write(self, filepath: str) -> None:
#         """Write sketch to file."""
#         self.sketch.write(filepath)

#     @classmethod
#     def load(cls, filepath: str) -> 'MinimizerSketch':
#         """Load sketch from file."""
#         # Load parameters from the underlying sketch file
#         data = np.load(filepath)
#         sketch = cls(
#             kmer_size=int(data['kmer_size'][0]),
#             window_size=int(data['window_size'][0]),
#             precision=int(data['precision'][0]) if 'precision' in data else 8,
#             num_hashes=int(data['num_hashes'][0]) if 'num_hashes' in data else 128,
#             seed=int(data['seed'][0])
#         )
#         sketch.sketch = HyperLogLog.load(filepath)
#         return sketch