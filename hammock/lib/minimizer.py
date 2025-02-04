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
            
        hash_sim = self.minimizer_sketch.estimate_jaccard(other.minimizer_sketch)
        gap_sim = self.gap_sketch.estimate_jaccard(other.gap_sketch)
        return hash_sim, gap_sim
        
    def estimate_jaccard(self, other: 'MinimizerSketch') -> float:
        """Estimate Jaccard similarity with another minimizer sketch."""
        if not isinstance(other, MinimizerSketch):
            raise TypeError("Can only compare with another minimizer sketch")
        
        # Return 0 if either sketch is empty
        if not self.minimizers or not other.minimizers:
            return 0.0
        
        # Calculate intersection and union sizes
        intersection = len(set(self.minimizers) & set(other.minimizers))
        union = len(set(self.minimizers) | set(other.minimizers))
        
        return float(intersection) / union if union > 0 else 0.0

    def similarity_values(self, other: 'AbstractSketch') -> Dict[str, float]:
        """Calculate similarity values between minimizer sketches.
        
        Returns:
            Dictionary containing 'jaccard_similarity' and 'overlap_similarity'
        """
        if not isinstance(other, MinimizerSketch):
            raise ValueError("Can only compare with another MinimizerSketch")
        if self.k != other.k or self.w != other.w or self.gapk != other.gapk:
            raise ValueError("Cannot compare sketches with different parameters")
        
        gap_sim, jaccard = self.compare_overlaps(other)
        
        return {
            'jaccard_similarity': jaccard,
            'gap_similarity': gap_sim
        }