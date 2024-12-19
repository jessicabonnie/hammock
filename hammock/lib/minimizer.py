#!/usr/bin/env python
from Bio import SeqIO
from Digest import window_minimizer
from typing import Optional, Union, Literal
from hammock.lib.abstractsketch import AbstractSketch
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash

class MinimizerSketch(AbstractSketch):
    """Sketch class for sequence data using minimizers with an underlying sketch type."""
    
    def __init__(self, 
                 window_size: int = 40, 
                 kmer_size: int = 8,
                 sketch_type: Literal["hyperloglog", "minhash"] = "hyperloglog",
                 precision: int = 8,
                 num_hashes: int = 128,
                 seed: int = 0):
        """Initialize MinimizerSketch.
        
        Args:
            window_size: Size of sliding window for minimizer selection
            kmer_size: Size of kmers to use
            sketch_type: Type of underlying sketch ("hyperloglog" or "minhash")
            precision: Precision for HyperLogLog sketching
            num_hashes: Number of hash functions for MinHash sketching
            seed: Random seed for hash functions
        """
        self.window_size = window_size
        self.kmer_size = kmer_size
        self.seed = seed
        
        # Initialize underlying sketch
        if sketch_type == "hyperloglog":
            self.sketch = HyperLogLog(precision=precision, seed=seed)
        elif sketch_type == "minhash":
            self.sketch = MinHash(num_hashes=num_hashes, seed=seed)
        else:
            raise ValueError(f"Invalid sketch type: {sketch_type}")
    
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        minimizers = window_minimizer(s, 
                                    w=self.window_size, 
                                    k=self.kmer_size, 
                                    include_hash=True)
        for _, hash_val in minimizers:
            self.sketch.add_int(hash_val)
            
    def add_sequence(self, sequence: str) -> None:
        """Add a sequence using the minimizer approach."""
        self.add_string(sequence)
        # Add concatenated flanking kmers from start and end of sequence
        self.sketch.add_string(sequence[:self.kmer_size] + sequence[-self.kmer_size:])
        
    def estimate_cardinality(self) -> float:
        """Estimate cardinality using underlying sketch."""
        return self.sketch.estimate_cardinality()
        
    def estimate_jaccard(self, other: 'AbstractSketch') -> float:
        """Estimate Jaccard similarity with another sketch."""
        if not isinstance(other, MinimizerSketch):
            raise ValueError("Can only compare with another MinimizerSketch")
        return self.sketch.estimate_jaccard(other.sketch)
        
    def merge(self, other: 'AbstractSketch') -> None:
        """Merge another sketch into this one."""
        if not isinstance(other, MinimizerSketch):
            raise ValueError("Can only merge with another MinimizerSketch")
        self.sketch.merge(other.sketch)