#!/usr/bin/env python
from typing import Optional, List, Iterator
from Bio import SeqIO # type: ignore
from hammock.lib.abstractsketch import AbstractDataSketch
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash
from hammock.lib.minimizer import MinimizerSketch
import gc

class SequenceSketch(AbstractDataSketch):
    """Sketch class for sequence data using various sketching methods."""
    
    def __init__(self, 
                 sketch_type: str = "minimizer",
                 kmer_size: int = 8,
                 window_size: int = 40,
                 precision: int = 8,
                 num_hashes: int = 128,
                 seed: int = 0):
        """Initialize sequence sketch.
        
        Args:
            sketch_type: Type of sketch to use (minimizer/hyperloglog/minhash)
            kmer_size: Size of kmers
            window_size: Size of sliding window
            precision: Precision for HyperLogLog sketching
            num_hashes: Number of hash functions for MinHash sketching
            seed: Random seed for hash functions
        """
        if sketch_type not in ["minimizer", "hyperloglog", "minhash"]:
            raise ValueError(f"Invalid sketch type for sequences: {sketch_type}")
            
        if sketch_type == "minimizer":
            self.sketch = MinimizerSketch(kmer_size=kmer_size, 
                                        window_size=window_size,
                                        seed=seed)
        elif sketch_type == "hyperloglog":
            self.sketch = HyperLogLog(precision=precision,
                                    kmer_size=kmer_size,
                                    seed=seed)
        elif sketch_type == "minhash":
            self.sketch = MinHash(num_hashes=num_hashes,
                                kmer_size=kmer_size,
                                seed=seed)
        
        self.kmer_size = kmer_size
        self.window_size = window_size
        self.total_sequence_length = 0
        self.num_sequences = 0
        
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
                **kwargs
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