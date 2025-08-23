#!/usr/bin/env python
from __future__ import annotations
import gc
import gzip
from typing import Dict, Iterator, List, Optional, Union, TYPE_CHECKING

from hammock.lib.abstractsketch import AbstractSketch
from hammock.lib.utils import read_sequences
from hammock.lib.minimizer import MinimizerSketch
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash

# Try to import FastHyperLogLog for better performance
try:
    from hammock.lib.hyperloglog_fast import FastHyperLogLog
    FAST_HLL_AVAILABLE = True
except ImportError:
    FAST_HLL_AVAILABLE = False
    FastHyperLogLog = None

# Use TYPE_CHECKING to avoid circular imports
if TYPE_CHECKING:
    from hammock.lib.minimizer import MinimizerSketch
    from hammock.lib.hyperloglog import HyperLogLog

class SequenceSketch(AbstractSketch):
    """Base class for sequence sketches."""
    
    def __init__(self, 
                 sketch_type: str = "minimizer",
                 kmer_size: int = 8,
                 window_size: int = 40,
                #  gapn: int = 0,
                 precision: int = 12,
                 num_hashes: int = 128,
                 seed: int = 42,
                 debug: bool = False,
                 use_sets: bool = False):
        """Initialize a SequenceSketch.
        
        Args:
            sketch_type: Type of sketch to use ("minimizer", "hyperloglog", or "minhash")
            kmer_size: Size of k-mers
            window_size: Size of sliding window
            gapn: Number of consecutive gaps to group (for minimizer sketch)
            precision: Precision parameter for HyperLogLog
            num_hashes: Number of hash functions for MinHash
            seed: Random seed for hashing
            debug: Whether to print debug information
            use_sets: Whether to store exact sets (only applies to minimizer sketch)
        """
        super().__init__()
        self.sketch_type = sketch_type
        self.kmer_size = kmer_size
        self.window_size = window_size
        self.seed = seed
        self.debug = debug
        self.use_sets = use_sets
        
        # Create the appropriate sketch based on sketch_type
        if sketch_type == "minimizer":
            self.sketch = MinimizerSketch(
                kmer_size=kmer_size,
                window_size=window_size,
                #  gapn=gapn,
                seed=seed,
                precision=precision,
                debug=debug,
                use_sets=use_sets
            )
        elif sketch_type == "hyperloglog":
            # Use FastHyperLogLog if available for better performance
            if FAST_HLL_AVAILABLE:
                self.sketch = FastHyperLogLog(
                    precision=precision,
                    kmer_size=kmer_size,
                    window_size=window_size,
                    seed=seed,
                    debug=debug
                )
            else:
                self.sketch = HyperLogLog(
                    precision=precision,
                    kmer_size=kmer_size,
                    window_size=window_size,
                    seed=seed,
                    debug=debug
                )
        elif sketch_type == "minhash":
            self.sketch = MinHash(
                num_hashes=num_hashes,
                kmer_size=kmer_size,
                window_size=window_size,
                seed=seed,
                debug=debug
            )
        else:
            raise ValueError(f"Unknown sketch type: {sketch_type}")
    
    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        self.sketch.add_string(s)
    
    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch.
        
        Args:
            strings: List of strings to add to the sketch
        """
        for s in strings:
            self.sketch.add_string(s)
    
    def estimate_jaccard(self, other: 'AbstractSketch') -> float:
        """Estimate Jaccard similarity with another sketch."""
        if not isinstance(other, SequenceSketch):
            raise TypeError("Can only compare with another SequenceSketch")
        
        return self.sketch.estimate_jaccard(other.sketch)
    
    def estimate_cardinality(self) -> float:
        """Estimate the number of unique k-mers in the sequences."""
        return self.sketch.estimate_cardinality()
    
    def similarity_values(self, other: 'AbstractSketch') -> Dict[str, float]:
        """Calculate similarity values between sketches.
        
        Returns:
            Dictionary containing similarity metrics
        """
        if not isinstance(other, SequenceSketch):
            raise ValueError("Can only compare with another SequenceSketch")
        
        # For MinimizerSketch, use its specialized similarity metrics
        if self.sketch_type == "minimizer" and other.sketch_type == "minimizer":
            return self.sketch.similarity_values(other.sketch)
        
        # For other sketch types, just return Jaccard similarity
        return {"jaccard_similarity": self.estimate_jaccard(other)}

    @classmethod
    def from_file(cls, filename: str, **kwargs) -> Optional['SequenceSketch']:
        """Create a sketch from a FASTA/FASTQ file."""
        try:
            # Create a SequenceSketch with the appropriate sketch type
            sketch = cls(
                sketch_type=kwargs.get('sketch_type', 'minimizer'),
                kmer_size=kwargs.get('kmer_size', 8),
                window_size=kwargs.get('window_size', 40),
                #  gapn=kwargs.get('gapn', 0),
                precision=kwargs.get('precision', 12),
                num_hashes=kwargs.get('num_hashes', 128),
                seed=kwargs.get('seed', 42),
                debug=kwargs.get('debug', False),
                use_sets=kwargs.get('use_sets', False)
            )
            
            # Process the file
            for records in read_sequences(filename, kwargs.get('chunk_size', 1000)):
                for record in records:
                    sketch.add_string(str(record.seq))
                if kwargs.get('verbose', False):
                    print(f"Processed {len(records)} sequences from {filename}")
                gc.collect()
            
            return sketch
            
        except Exception as e:
            print(f"Error processing file {filename}: {str(e)}")
            return None

    def write(self, filepath: str) -> None:
        """Write sketch to file."""
        self.sketch.write(filepath)

    @classmethod
    def load(cls, filepath: str) -> 'SequenceSketch':
        """Load sketch from file."""
        # Try each sketch type and catch specific exceptions
        try:
            sketch = HyperLogLog.load(filepath)
            seq_sketch = cls(sketch_type="hyperloglog")
            seq_sketch.sketch = sketch
            return seq_sketch
        except (ValueError, OSError):
            pass

        try:
            sketch = MinHash.load(filepath)
            seq_sketch = cls(sketch_type="minhash")
            seq_sketch.sketch = sketch
            return seq_sketch
        except (ValueError, OSError):
            pass

        try:
            sketch = MinimizerSketch.load(filepath)
            seq_sketch = cls(sketch_type="minimizer")
            seq_sketch.sketch = sketch
            return seq_sketch
        except (ValueError, OSError):
            pass

        raise ValueError(f"Could not load sketch from file: {filepath}")