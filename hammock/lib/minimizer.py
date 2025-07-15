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
                 precision: int = 16,
                 debug: bool = False,
                 use_sets: bool = False):
        """Initialize minimizer sketch.
        
        Args:
            window_size: Size of sliding window
            kmer_size: Size of k-mers
            seed: Random seed for hashing
            precision: Precision of the HyperLogLog sketchs
            debug: Whether to print debug information
            use_sets: Whether to store minimizers in sets (exact) in addition to sketches
        """
        super().__init__()
        self.window_size = window_size
        self.kmer_size = kmer_size
        self.seed = seed
        self.precision = precision
        self.debug = debug
        self.use_sets = use_sets
        # Initialize HyperLogLog sketches with proper parameters
        self.minimizer_sketch = HyperLogLog(precision=precision, kmer_size=kmer_size, window_size=window_size, seed=seed)
        self.startend_sketch = HyperLogLog(precision=precision, kmer_size=kmer_size, window_size=window_size, seed=seed)
        # Initialize sets only if use_sets is True
        if self.use_sets:
            self.startend_kmers = set()
            self.minimizers = set()
    
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        # Get minimizers from the string
        minimizers = window_minimizer(s, self.window_size, self.kmer_size, self.seed)
        # TODO: what happens if the string is shorter than the window size?
        # Add minimizers to the sketch
        for _, hash_val in minimizers:
            if self.use_sets:
                self.minimizers.add(hash_val)
            self.minimizer_sketch.add_string(str(hash_val))
        
        # Add start and end k-mers
        if len(s) >= self.kmer_size:
            start_kmer = s[:self.kmer_size]
            end_kmer = s[-self.kmer_size:]
            self.startend_sketch.add_string(start_kmer + end_kmer)
            if self.use_sets:
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
            
        # Calculate set similarity using minimizer sets only (if using sets)
        if self.use_sets and other.use_sets:
            intersection = len(self.minimizers & other.minimizers)
            union = len(self.minimizers | other.minimizers)
            set_sim = float(intersection) / union if union > 0 else 0.0
            
            if self.debug:
                print(f"Set similarity - intersection: {intersection}, union: {union}, similarity: {set_sim:.4f}")
        else:
            set_sim = None
        
        # Calculate set similarity including end k-mers (if using sets)
        if self.use_sets and other.use_sets:
            combined_self = self.minimizers | self.startend_kmers
            combined_other = other.minimizers | other.startend_kmers
            intersection = len(combined_self & combined_other)
            union = len(combined_self | combined_other)
            set_with_ends_sim = float(intersection) / union if union > 0 else 0.0
            
            if self.debug:
                print(f"Set+ends similarity - intersection: {intersection}, union: {union}, similarity: {set_with_ends_sim:.4f}")
        else:
            set_with_ends_sim = None
        # Calculate sketch-based similarity (always available)
        combined_self_sketch = self.minimizer_sketch.merge_new(self.startend_sketch)
        combined_other_sketch = other.minimizer_sketch.merge_new(other.startend_sketch)
        
        # Calculate overall Jaccard similarity
        jaccard_sim = self.estimate_jaccard(other)
        jaccard_sim_with_ends = combined_self_sketch.estimate_jaccard(combined_other_sketch)
        
        if self.debug:
            print(f"\nSimilarity metrics:")
            if set_sim is not None:
                print(f"Set similarity: {set_sim:.4f}")
            if set_with_ends_sim is not None:
                print(f"Set+ends similarity: {set_with_ends_sim:.4f}")
            print(f"Jaccard similarity: {jaccard_sim:.4f}")
            print(f"Jaccard similarity with ends: {jaccard_sim_with_ends:.4f}")
        
        # Build return dictionary - include set-based metrics only if available
        result = {
            'jaccard_similarity': jaccard_sim,
            'jaccard_similarity_with_ends': jaccard_sim_with_ends
        }
        
        if set_sim is not None:
            result['set_similarity'] = set_sim
        if set_with_ends_sim is not None:
            result['set_with_ends_similarity'] = set_with_ends_sim
            
        return result
        
    def estimate_cardinality(self) -> float:
        """Estimate the number of unique minimizers in the sequences."""
        if self.use_sets:
            return float(len(self.minimizers))
        else:
            return self.minimizer_sketch.estimate_cardinality()
    
    def merge(self, other: 'MinimizerSketch') -> 'MinimizerSketch':
        """Merge another minimizer sketch into this one."""
        if not isinstance(other, MinimizerSketch):
            raise TypeError("Can only merge with another MinimizerSketch")
        
        result = MinimizerSketch(
            kmer_size=self.kmer_size,
            window_size=self.window_size,
            seed=self.seed,
            precision=self.precision,
            debug=self.debug,
            use_sets=self.use_sets
        )
        
        # Merge minimizer sets (if using sets)
        if self.use_sets and other.use_sets:
            result.minimizers = self.minimizers.union(other.minimizers)
            result.startend_kmers = self.startend_kmers.union(other.startend_kmers)
        
        # Merge sketches
        print(f"Merging {self.minimizer_sketch.estimate_cardinality()} and {other.minimizer_sketch.estimate_cardinality()} cardinality sketches")
        result.minimizer_sketch = self.minimizer_sketch.merge_new(other.minimizer_sketch)
        result.startend_sketch = self.startend_sketch.merge_new(other.startend_sketch)
        
        return result

    def write(self, filepath: str) -> None:
        """Write sketch to file."""
        # Save the sketches to a single npz file
        save_data = {
            # Save HyperLogLog sketches
            'minimizer_sketch_registers': self.minimizer_sketch.registers,
            'startend_sketch_registers': self.startend_sketch.registers,
            
            # Save parameters
            'kmer_size': self.kmer_size,
            'window_size': self.window_size,
            'seed': self.seed,
            'precision': self.precision,
            'debug': self.debug,
            'use_sets': self.use_sets
        }
        
        # Save sets only if using sets
        if self.use_sets:
            save_data['minimizers'] = list(self.minimizers)
            save_data['startend_kmers'] = list(self.startend_kmers)
        
        np.savez(filepath, **save_data)

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
                precision=int(data['precision']) if 'precision' in data else 16,
                debug=bool(data['debug']),
                use_sets=bool(data['use_sets']) if 'use_sets' in data else False
            )
            
            # Restore the sets (if they were saved)
            if 'minimizers' in data and 'startend_kmers' in data:
                sketch.minimizers = set(data['minimizers'])
                sketch.startend_kmers = set(data['startend_kmers'])
            
            # Restore the registers for each sketch
            sketch.minimizer_sketch.registers = data['minimizer_sketch_registers']
            sketch.startend_sketch.registers = data['startend_sketch_registers']
            
            return sketch
            
        except Exception as e:
            raise ValueError(f"Could not load sketch from file: {str(e)}")