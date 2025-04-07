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
        # Check if the string is empty
        if not s:
            return
        
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
        
        # Check if sketches are None before merging
        if self.minimizer_sketch is None or self.startend_sketch is None:
            hash_with_ends_sim = 0.0
        else:
            try:
                combined_self_sketch = self.minimizer_sketch.merge(self.startend_sketch)
                combined_other_sketch = other.minimizer_sketch.merge(other.startend_sketch)
                intersection = len(combined_self & combined_other)
                union = len(combined_self | combined_other)
                hash_with_ends_sim = float(intersection) / union if union > 0 else 0.0
            except (AttributeError, TypeError):
                # If merge fails, fall back to set-based calculation
                intersection = len(combined_self & combined_other)
                union = len(combined_self | combined_other)
                hash_with_ends_sim = float(intersection) / union if union > 0 else 0.0
        
        if self.debug:
            print(f"Hash+ends similarity - intersection: {intersection}, union: {union}, similarity: {hash_with_ends_sim:.4f}")
        
        # Calculate gap similarity using HyperLogLog
        if self.gap_sketch is None or other.gap_sketch is None:
            gap_sim = 0.0
        else:
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
        
        # Create new sketches for the result
        result.minimizer_sketch = HyperLogLog(
            precision=self.minimizer_sketch.precision,
            kmer_size=self.minimizer_sketch.kmer_size,
            window_size=self.minimizer_sketch.window_size,
            seed=self.minimizer_sketch.seed,
            hash_size=self.minimizer_sketch.hash_size
        )
        result.gap_sketch = HyperLogLog(
            precision=self.gap_sketch.precision,
            kmer_size=self.gap_sketch.kmer_size,
            window_size=self.gap_sketch.window_size,
            seed=self.gap_sketch.seed,
            hash_size=self.gap_sketch.hash_size
        )
        result.startend_sketch = HyperLogLog(
            precision=self.startend_sketch.precision,
            kmer_size=self.startend_sketch.kmer_size,
            window_size=self.startend_sketch.window_size,
            seed=self.startend_sketch.seed,
            hash_size=self.startend_sketch.hash_size
        )
        
        # Copy registers from self
        result.minimizer_sketch.registers = self.minimizer_sketch.registers.copy()
        result.gap_sketch.registers = self.gap_sketch.registers.copy()
        result.startend_sketch.registers = self.startend_sketch.registers.copy()
        
        # Merge with other
        result.minimizer_sketch.merge(other.minimizer_sketch)
        result.gap_sketch.merge(other.gap_sketch)
        result.startend_sketch.merge(other.startend_sketch)
        
        return result

    def write(self, filepath: str) -> None:
        """Write sketch to file."""
        # Save the sketches to a single npz file
        np.savez(filepath,
                 # Save HyperLogLog sketches
                 minimizer_sketch_registers=self.minimizer_sketch.registers,
                 gap_sketch_registers=self.gap_sketch.registers,
                 startend_sketch_registers=self.startend_sketch.registers,
                 
                 # Save sets
                 minimizers=list(self.minimizers),
                 startend_kmers=list(self.startend_kmers),
                 
                 # Save parameters
                 kmer_size=self.kmer_size,
                 window_size=self.window_size,
                 gapn=self.gapn,
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
                gapn=int(data['gapn']),
                seed=int(data['seed']),
                debug=bool(data['debug'])
            )
            
            # Restore the sets
            sketch.minimizers = set(data['minimizers'])
            sketch.startend_kmers = set(data['startend_kmers'])
            
            # Restore the registers for each sketch
            sketch.minimizer_sketch.registers = data['minimizer_sketch_registers']
            sketch.gap_sketch.registers = data['gap_sketch_registers']
            sketch.startend_sketch.registers = data['startend_sketch_registers']
            
            return sketch
            
        except Exception as e:
            raise ValueError(f"Could not load sketch from file: {str(e)}")