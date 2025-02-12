#!/usr/bin/env python
from typing import Optional, List, Iterator, Literal, Dict
from Bio import SeqIO # type: ignore
from hammock.lib.abstractsketch import AbstractSketch
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash
from hammock.lib.minimizer import MinimizerSketch
import gc

class SequenceSketch(AbstractSketch):
    """Sketch class for sequence data using various sketching methods."""
    
    def __init__(self, 
                 sketch_type: str = "minimizer",
                 kmer_size: int = 8,
                 window_size: int = 40,
                 precision: int = 8,
                 num_hashes: int = 128,
                 seed: int = 42,
                 debug: bool = False):
        """Initialize a SequenceSketch.
        
        Args:
            sketch_type: Type of sketch to use
            kmer_size: Size of kmers
            window_size: Size of sliding window
            precision: Precision for HyperLogLog
            num_hashes: Number of hash functions for MinHash
            seed: Seed for random number generation
            debug: Whether to print debug information
        """
        super().__init__()
        self.sketch_type = sketch_type
        self.kmer_size = kmer_size
        self.debug = debug
        
        if sketch_type == "hyperloglog":
            self.sketch = HyperLogLog(kmer_size=kmer_size, precision=precision, seed=seed, debug=debug)
        elif sketch_type == "minhash":
            self.sketch = MinHash(kmer_size=kmer_size, num_hashes=num_hashes, seed=seed, debug=debug)
        elif sketch_type == "minimizer":
            self.sketch = MinimizerSketch(kmer_size=kmer_size, window_size=window_size, seed=seed, debug=debug)
        else:
            raise ValueError(f"Invalid sketch type: {sketch_type}")
        
        self.window_size = window_size
        self.total_sequence_length = 0
        self.num_sequences = 0
        self.seed = seed
        
    def add_sequence(self, sequence: str) -> None:
        """Add a sequence to the sketch.
        
        Args:
            sequence: DNA/RNA sequence string
        """
        self.total_sequence_length += len(sequence)
        self.num_sequences += 1
        
        # Process sequence based on sketch type
        if isinstance(self.sketch, MinimizerSketch):
            self.sketch.add_sequence(sequence)
        else:
            # For HyperLogLog and MinHash, use k-mer approach
            for i in range(len(sequence) - self.kmer_size + 1):
                kmer = sequence[i:i + self.kmer_size]
                self.sketch.add_string(kmer)
    
    @classmethod
    def from_file(cls,
                  filename: str,
                  sketch_type: str = "minimizer",
                  kmer_size: int = 8,
                  window_size: int = 40,
                  precision: int = 8,
                  num_hashes: int = 128,
                  chunk_size: int = 1000,
                  verbose: bool = False,
                  **kwargs) -> Optional['SequenceSketch']:
        """Create a SequenceSketch from a FASTA/FASTQ file.
        
        Args:
            filename: Path to FASTA/FASTQ file
            sketch_type: Type of sketch to use
            kmer_size: Size of kmers
            window_size: Size of sliding window
            precision: Precision for HyperLogLog
            num_hashes: Number of hash functions for MinHash
            chunk_size: Number of sequences to process at once
            verbose: Whether to print progress
            
        Returns:
            SequenceSketch object or None if file processing fails
        """
        try:
            sketch = cls(
                sketch_type=sketch_type,
                kmer_size=kmer_size,
                window_size=window_size,
                precision=precision,
                num_hashes=num_hashes,
            )
            
            for records in read_sequences(filename, chunk_size):
                for record in records:
                    sketch.add_sequence(str(record.seq))
                if verbose:
                    print(f"Processed {chunk_size} sequences from {filename}")
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
        # First try loading as HyperLogLog
        try:
            sketch = HyperLogLog.load(filepath)
            seq_sketch = cls(sketch_type="hyperloglog")
            seq_sketch.sketch = sketch
            return seq_sketch
        except:
            pass

        # Then try as MinHash
        try:
            sketch = MinHash.load(filepath)
            seq_sketch = cls(sketch_type="minhash")
            seq_sketch.sketch = sketch
            return seq_sketch
        except:
            pass

        # Finally try as Minimizer
        try:
            sketch = MinimizerSketch.load(filepath)
            seq_sketch = cls(sketch_type="minimizer")
            seq_sketch.sketch = sketch
            return seq_sketch
        except:
            raise ValueError("Could not load sketch from file")

    def estimate_cardinality(self) -> float:
        """Estimate the number of unique k-mers in the sequences.
        
        Returns:
            Estimated number of unique k-mers
        """
        return self.sketch.estimate_cardinality()

    def estimate_jaccard(self, other: 'SequenceSketch') -> float:
        """Estimate Jaccard similarity with another sketch."""
        if not isinstance(other, SequenceSketch):
            raise TypeError("Can only compare with another SequenceSketch")
        if self.sketch_type != other.sketch_type:
            raise ValueError("Cannot compare sketches of different types")
        return self.sketch.estimate_jaccard(other.sketch)

    def add_string(self, s: str) -> None:
        """Add a string to the sketch."""
        self.sketch.add_string(s)
    
    def similarity_values(self, other: 'AbstractSketch') -> Dict[str, float]:
        """Calculate similarity values between sequence sketches.
        
        Returns:
            Dictionary containing all similarity measures from the underlying sketch
        """
        if not isinstance(other, SequenceSketch):
            raise ValueError("Can only compare with another SequenceSketch")
        
        # Pass through all similarity values from the underlying sketch
        return self.sketch.similarity_values(other.sketch)

def read_sequences(filename: str, chunk_size: int = 1000) -> Iterator[List[SeqIO.SeqRecord]]:
    """Read sequences from a FASTA/FASTQ file in chunks."""
    formatx = "fasta" if filename.endswith((".fa", ".fasta")) else "fastq"
    records = []
    
    with open(filename, "r") as file:
        for record in SeqIO.parse(file, formatx):
            records.append(record)
            if len(records) >= chunk_size:
                yield records
                records = []
        if records:  # Yield any remaining records
            yield records