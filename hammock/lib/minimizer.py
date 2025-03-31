#!/usr/bin/env python
from Bio import SeqIO # type: ignore
from digest import window_minimizer # type: ignore
from typing import Optional, Union, Literal, Dict, List
from hammock.lib.abstractsketch import AbstractSketch
from hammock.lib.hyperloglog import HyperLogLog
import numpy as np

class MinimizerSketch(AbstractSketch):
    def __init__(self, 
                 kmer_size: int = 4, 
                 window_size: int = 8, 
                 gapn: int = 0,
                 seed: int = 42,
                 debug: bool = False):
        """Initialize minimizer sketch.
        
        Args:
            window_size: Size of sliding window
            kmer_size: Size of k-mers
            gapn: Number of consecutive gaps to group at a time
            seed: Random seed for hashing
            debug: Whether to print debug information
        """
        super().__init__()
        self.window_size = window_size
        self.kmer_size = kmer_size
        self.gapn = gapn
        self.seed = seed
        self.debug = debug
        # Initialize HyperLogLog sketches with proper parameters
        self.minimizer_sketch = HyperLogLog(kmer_size=kmer_size, window_size=window_size, seed=seed)
        self.gap_sketch = HyperLogLog(kmer_size=kmer_size, window_size=window_size, seed=seed)
        self.startend_sketch = HyperLogLog(kmer_size=kmer_size, window_size=window_size, seed=seed)
        self.startend_kmers = set()
        self.minimizers = set()
    
    def _process_gap_patterns(self, minimizers: list, debug: bool = False) -> None:
        """Process gap patterns from minimizer positions and add them to the gap sketch.
        
        Takes consecutive gaps between minimizer positions and groups them into patterns
        of length gapn. Each pattern is added to the gap sketch for similarity comparison.
        
        Args:
            minimizers: List of (position, hash) tuples from window_minimizer
            debug: Whether to print debug information
        """
        if len(minimizers) < self.gapn:
            if self.debug or debug:
                print("Not enough minimizers to calculate gap patterns")
            return 
        
        # Extract positions and calculate gaps
        positions = [pos for pos, _ in minimizers]
        gaps = [positions[i+1] - positions[i] for i in range(len(positions)-1)]
        
        if debug:
            print(f"Minimizer positions: {positions[:5]}...")
            print(f"Gaps between minimizers: {gaps[:5]}...")
        
        # Generate gap patterns of length gapn
        for i in range(len(gaps) - self.gapn + 1):
            pattern = ','.join(str(gaps[i:i + self.gapn]))
            self.gap_sketch.add_string(pattern)
        return
        
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        # Get minimizers from the string
        minimizers = window_minimizer(s, self.window_size, self.kmer_size, self.seed)
        
        # Add minimizers to the sketch
        for _, hash_val in minimizers:
            self.minimizer_sketch.add_string(str(hash_val))
            self.minimizers.add(hash_val)
        
        # Process gap patterns if gapn > 0
        if self.gapn > 0:
            self._process_gap_patterns(minimizers, self.debug)
        
        # Add start and end k-mers
        if len(s) >= self.kmer_size:
            start_kmer = s[:self.kmer_size]
            end_kmer = s[-self.kmer_size:]
            self.startend_sketch.add_string(start_kmer)
            self.startend_sketch.add_string(end_kmer)
            self.startend_kmers.add(start_kmer)
            self.startend_kmers.add(end_kmer)
    
    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch.
        
        Args:
            strings: List of strings to add to the sketch
        """
        for s in strings:
            self.add_string(s)
    
    def compare_overlaps(self, other: 'MinimizerSketch') -> Dict[str, float]:
        """Returns dictionary of similarity metrics between two minimizer sketches.
        
        Returns:
            Dictionary containing:
                - 'hash_similarity': Similarity based on minimizer hashes only
                - 'hash_with_ends_similarity': Similarity including end k-mers
                - 'gap_similarity': Similarity of gap patterns
                - 'jaccard_similarity': Overall Jaccard similarity
        """
        if not isinstance(other, MinimizerSketch):
            raise ValueError("Can only compare with another MinimizerSketch")
        if self.kmer_size != other.kmer_size or self.window_size != other.window_size or self.gapn != other.gapn:
            raise ValueError("Cannot compare sketches with different parameters")
            
        # Calculate hash similarity using minimizer sets only
        # NOTE: This should probably not be part of this function
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
        
        # Calculate gap similarity using HyperLogLog
        gap_sim = self.gap_sketch.estimate_jaccard(other.gap_sketch)
        
        # Calculate overall Jaccard similarity
        jaccard_sim = self.estimate_jaccard(other)
        
        if self.debug:
            print(f"\nSimilarity metrics:")
            print(f"Hash similarity: {hash_sim:.4f}")
            print(f"Hash+ends similarity: {hash_with_ends_sim:.4f}")
            print(f"Gap similarity: {gap_sim:.4f}")
            print(f"Jaccard similarity: {jaccard_sim:.4f}")
        
        return {
            'hash_similarity': hash_sim,
            'hash_with_ends_similarity': hash_with_ends_sim,
            'gap_similarity': gap_sim,
            'jaccard_similarity': jaccard_sim
        }

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
        if self.kmer_size != other.kmer_size or self.window_size != other.window_size or self.gapn != other.gapn:
            raise ValueError("Cannot compare sketches with different parameters")
        
        return self.compare_overlaps(other)

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
            gapn=self.gapn,
            seed=self.seed,
            debug=self.debug
        )
        
        # Merge minimizer sets
        result.minimizers = self.minimizers.union(other.minimizers)
        
        # Merge startend kmers
        result.startend_kmers = self.startend_kmers.union(other.startend_kmers)
        
        # Merge sketches
        result.minimizer_sketch = self.minimizer_sketch.merge(other.minimizer_sketch)
        result.gap_sketch = self.gap_sketch.merge(other.gap_sketch)
        result.startend_sketch = self.startend_sketch.merge(other.startend_sketch)
        
        return result

    def write(self, filepath: str) -> None:
        """Write the MinimizerSketch to a file.
        
        Args:
            filepath: Path to write the sketch to
        """
        import pickle
        with open(filepath, 'wb') as f:
            pickle.dump({
                'kmer_size': self.kmer_size,
                'window_size': self.window_size,
                'gapn': self.gapn,
                'seed': self.seed,
                'debug': self.debug,
                'minimizers': self.minimizers,
                'startend_kmers': self.startend_kmers,
                'minimizer_sketch': self.minimizer_sketch,
                'gap_sketch': self.gap_sketch,
                'startend_sketch': self.startend_sketch
            }, f)

    @classmethod
    def load(cls, filepath: str) -> 'MinimizerSketch':
        """Load a MinimizerSketch from a file.
        
        Args:
            filepath: Path to load the sketch from
            
        Returns:
            The loaded MinimizerSketch
        """
        import pickle
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
            
        sketch = cls(
            kmer_size=data['kmer_size'],
            window_size=data['window_size'],
            gapn=data['gapn'],
            seed=data['seed'],
            debug=data['debug']
        )
        
        # Restore the sets and sketches
        sketch.minimizers = data['minimizers']
        sketch.startend_kmers = data['startend_kmers']
        sketch.minimizer_sketch = data['minimizer_sketch']
        sketch.gap_sketch = data['gap_sketch']
        sketch.startend_sketch = data['startend_sketch']
        
        return sketch