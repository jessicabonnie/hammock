#!/usr/bin/env python
from Bio import SeqIO # type: ignore
from Digest import window_minimizer # type: ignore
from typing import Optional, Union, Literal, Dict
from hammock.lib.abstractsketch import AbstractSketch
from hammock.lib.hyperloglog import HyperLogLog

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
        # Use HyperLogLog directly for both sketches
        self.minimizer_sketch = HyperLogLog(kmer_size=0)
        self.gap_sketch = HyperLogLog(kmer_size=0)
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
        if not s or len(s) < self.kmer_size or len(s)==0:
            return
            
        # Get minimizers and their positions
        try:
            minimizers = window_minimizer(s, self.kmer_size, self.window_size, include_hash=True)
            
            if self.debug:
                print(f"\nDebug: add_string")
                print(f"Number of minimizers found: {len(minimizers)}")
                print(f"First few minimizers (pos, hash): {minimizers[:5]}")
            
            # 1. Add minimizer hashes to sketch
            for _, h in minimizers:
                self.minimizer_sketch.add_string(str(h))
                self.minimizers.add(h)
            # 2. Store start+end k-mers separately
            if len(s) >= self.kmer_size:
                self.startend_kmers.add(s[:self.kmer_size]+
                                    s[-self.kmer_size:])
                
                if self.debug:
                    print(f"End k-mers added: {s[:self.kmer_size]}, {s[-self.kmer_size:]}")
            
            # Calculate and store gap patterns
            # NOTE: currently cycling through minimizers more than once
            self._process_gap_patterns(minimizers)
            
        except RuntimeError as e:
            if self.debug:
                print(f"Error in add_string: {str(e)}")
                print(f"Parameters: k={self.kmer_size}, w={self.window_size}, gapn={self.gapn}")
            return
            
    def compare_overlaps(self, other: 'MinimizerSketch') -> tuple[float, float, float, float]:
        """Returns (hash_similarity, hash_with_ends_similarity, gap_similarity, jaccard_similarity)"""
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
        intersection = len(combined_self & combined_other)
        union = len(combined_self | combined_other)
        hash_with_ends_sim = float(intersection) / union if union > 0 else 0.0
        
        if self.debug:
            print(f"Hash+ends similarity - intersection: {intersection}, union: {union}, similarity: {hash_with_ends_sim:.4f}")
        
        # 3. Calculate gap similarity using HyperLogLog
        gap_sim = self.gap_sketch.estimate_jaccard(other.gap_sketch)
        
        # 4. Calculate overall Jaccard similarity
        jaccard_sim = self.estimate_jaccard(other)
        
        if self.debug:
            print(f"\nSimilarity metrics:")
            print(f"Hash similarity: {hash_sim:.4f}")
            print(f"Hash+ends similarity: {hash_with_ends_sim:.4f}")
            print(f"Gap similarity: {gap_sim:.4f}")
            print(f"Jaccard similarity: {jaccard_sim:.4f}")
        
        return hash_sim, hash_with_ends_sim, gap_sim, jaccard_sim
        
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
            Dictionary containing 'hash_similarity', 'hash_with_ends_similarity', 
            'gap_similarity', and 'jaccard_similarity'
        """
        if not isinstance(other, MinimizerSketch):
            raise ValueError("Can only compare with another MinimizerSketch")
        if self.kmer_size != other.kmer_size or self.window_size != other.window_size or self.gapn != other.gapn:
            raise ValueError("Cannot compare sketches with different parameters")
        # NOTE: Perhaps the compare_overlaps function should just return a dictionary
        hash_sim, hash_with_ends_sim, gap_sim, jaccard_sim = self.compare_overlaps(other)
        
        return {
            'hash_similarity': hash_sim,
            'hash_with_ends_similarity': hash_with_ends_sim,
            'gap_similarity': gap_sim,
            'jaccard_similarity': jaccard_sim
        }