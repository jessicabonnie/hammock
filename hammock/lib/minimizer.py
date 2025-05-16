#!/usr/bin/env python
from Bio import SeqIO # type: ignore
from digest import window_minimizer # type: ignore
from typing import Optional, Union, Literal, Dict, List
from hammock.lib.abstractsketch import AbstractSketch
from hammock.lib.hyperloglog import HyperLogLog
import numpy as np # type: ignore

class MinimizerSketch(AbstractSketch):
    def __init__(self, 
                 kmer_size: int = 4, 
                 window_size: int = 8, 
                 seed: int = 42,
                 debug: bool = False):
        """Initialize minimizer sketch.
        
        Args:
            window_size: Size of sliding window
            kmer_size: Size of k-mers
            seed: Random seed for hashing
            debug: Whether to print debug information
        """
        super().__init__()
        self.window_size = window_size
        self.kmer_size = kmer_size
        self.seed = seed
        self.debug = debug
        # Initialize HyperLogLog sketches with proper parameters
        self.minimizer_sketch = HyperLogLog(kmer_size=kmer_size, window_size=window_size, seed=seed)
        self.startend_sketch = HyperLogLog(kmer_size=kmer_size, window_size=window_size, seed=seed)
        self.startend_kmers = set()
        self.minimizers = set()
    
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        # Get minimizers from the string
        minimizers = window_minimizer(s, self.window_size, self.kmer_size, self.seed)
        
        # Add minimizers to the sketch
        for _, hash_val in minimizers:
            self.minimizers.add(hash_val)
        
        # Add start and end k-mers
        if len(s) >= self.kmer_size:
            start_kmer = s[:self.kmer_size]
            end_kmer = s[-self.kmer_size:]
            self.startend_kmers.add(start_kmer)
            self.startend_kmers.add(end_kmer)
    
    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch.
        
        Args:
            strings: List of strings to add to the sketch
        """
        for s in strings:
            self.add_string(s)

    def estimate_jaccard(self, other: 'MinimizerSketch') -> float:
        """Estimate Jaccard similarity with another minimizer sketch."""
        if not isinstance(other, MinimizerSketch):
            raise TypeError("Can only compare with another minimizer sketch")
        
        # Return 0 if either sketch is empty
        if not self.minimizer_sketch or not other.minimizer_sketch:
            return 0.0
        
        return self.minimizer_sketch.estimate_jaccard(other.minimizer_sketch)

    def similarity_values(self, other: 'AbstractSketch') -> Dict[str, float]:
        """Calculate similarity values between minimizer sketches.
        
        Returns:
            Dictionary containing similarity metrics
        """
        if not isinstance(other, MinimizerSketch):
            raise ValueError("Can only compare with another MinimizerSketch")
        if self.kmer_size != other.kmer_size or self.window_size != other.window_size:
            raise ValueError("Cannot compare sketches with different parameters")
            
        # Calculate hash similarity using minimizer sets only
        intersection = len(self.minimizers & other.minimizers)
        union = len(self.minimizers | other.minimizers)
        hash_sim = float(intersection) / union if union > 0 else 0.0
        
        if self.debug:
            print(f"Hash similarity - intersection: {intersection}, union: {union}, similarity: {hash_sim:.4f}")
        
        # Calculate hash similarity including end k-mers
        combined_self = self.minimizers | self.startend_kmers
        combined_other = other.minimizers | other.startend_kmers
        combined_self_sketch = self.minimizer_sketch.merge(self.startend_sketch)
        combined_other_sketch = other.minimizer_sketch.merge(other.startend_sketch)
        intersection = len(combined_self & combined_other)
        union = len(combined_self | combined_other)
        hash_with_ends_sim = float(intersection) / union if union > 0 else 0.0
        
        if self.debug:
            print(f"Hash+ends similarity - intersection: {intersection}, union: {union}, similarity: {hash_with_ends_sim:.4f}")
        
        # Calculate overall Jaccard similarity
        jaccard_sim = self.estimate_jaccard(other)
        
        if self.debug:
            print(f"\nSimilarity metrics:")
            print(f"Hash similarity: {hash_sim:.4f}")
            print(f"Hash+ends similarity: {hash_with_ends_sim:.4f}")
            print(f"Jaccard similarity: {jaccard_sim:.4f}")
        
        return {
            'hash_similarity': hash_sim,
            'hash_with_ends_similarity': hash_with_ends_sim,
            'jaccard_similarity': jaccard_sim
        }
        
    def estimate_cardinality(self) -> float:
        """Estimate the number of unique minimizers in the sequences."""
        return float(len(self.minimizers))
    
    def merge(self, other: 'MinimizerSketch') -> 'MinimizerSketch':
        """Merge another minimizer sketch into this one."""
        if not isinstance(other, MinimizerSketch):
            raise TypeError("Can only merge with another MinimizerSketch")
        
        result = MinimizerSketch(
            kmer_size=self.kmer_size,
            window_size=self.window_size,
            seed=self.seed,
            debug=self.debug
        )
        
        # Merge minimizer sets
        result.minimizers = self.minimizers.union(other.minimizers)
        
        # Merge startend kmers
        result.startend_kmers = self.startend_kmers.union(other.startend_kmers)
        
        # Merge sketches
        result.minimizer_sketch = self.minimizer_sketch.merge(other.minimizer_sketch)
        result.startend_sketch = self.startend_sketch.merge(other.startend_sketch)
        
        return result

    def write(self, filepath: str) -> None:
        """Write sketch to file."""
        # Save the sketches to a single npz file
        np.savez(filepath,
                 # Save HyperLogLog sketches
                 minimizer_sketch_registers=self.minimizer_sketch.registers,
                 startend_sketch_registers=self.startend_sketch.registers,
                 
                 # Save sets
                 minimizers=list(self.minimizers),
                 startend_kmers=list(self.startend_kmers),
                 
                 # Save parameters
                 kmer_size=self.kmer_size,
                 window_size=self.window_size,
                 seed=self.seed,
                 debug=self.debug)

    @classmethod
    def load(cls, filepath: str) -> 'MinimizerSketch':
        """Load sketch from file."""
        try:
            # Load the npz file
            data = np.load(filepath + '.npz', allow_pickle=True)
            
            # Create new sketch with loaded parameters
            sketch = cls(
                kmer_size=int(data['kmer_size']),
                window_size=int(data['window_size']),
                seed=int(data['seed']),
                debug=bool(data['debug'])
            )
            
            # Restore the sets
            sketch.minimizers = set(data['minimizers'])
            sketch.startend_kmers = set(data['startend_kmers'])
            
            # Restore the registers for each sketch
            sketch.minimizer_sketch.registers = data['minimizer_sketch_registers']
            sketch.startend_sketch.registers = data['startend_sketch_registers']
            
            return sketch
            
        except Exception as e:
            raise ValueError(f"Could not load sketch from file: {str(e)}")