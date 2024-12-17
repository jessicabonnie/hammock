#!/usr/bin/env python
from typing import Optional, List, Iterator
from Bio import SeqIO # type: ignore
from hammock.lib.sketchclass import Sketch
from Digest import window_minimizer # type: ignore  
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash
import gc

class SequenceSketch(Sketch):
    """Sketch class for sequence data using various sketching methods."""
    
    def __init__(self, 
                 kmer_size: int = 8,
                 window_size: int = 40,
                 precision: int = 8,
                 num_hashes: int = 128,
                 sketch_type: str = "hyperloglog",
                 seed: int = 0):
        """Initialize SequenceSketch.
        
        Args:
            kmer_size: Size of kmers to use
            window_size: Size of sliding window
            precision: Precision for HyperLogLog sketching
            num_hashes: Number of hash functions for MinHash sketching
            sketch_type: Type of sketch to use (hyperloglog/minhash/minimizer)
            seed: Random seed for hash functions
        """
        if sketch_type not in ["hyperloglog", "minhash", "minimizer"]:
            raise ValueError(f"Invalid sketch type: {sketch_type}")
            
        super().__init__(
            sketch_type=sketch_type,
            precision=precision,
            num_hashes=num_hashes,
            kmer_size=kmer_size,
            seed=seed
        )
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
        
        if self.sketch_type == "minimizer":
            # Use minimizer approach
            minimizers = window_minimizer(sequence, 
                                        w=self.window_size, 
                                        k=self.kmer_size, 
                                        include_hash=True)
            for _, hash_val in minimizers:
                self.add_hash(hash_val)
        else:
            # Use regular k-mer approach for hyperloglog/minhash
            for i in range(len(sequence) - self.kmer_size + 1):
                kmer = sequence[i:i + self.kmer_size]
                self.add_string(kmer)
    
    @classmethod
    def from_file(cls,
                  filename: str,
                  kmer_size: int = 8,
                  window_size: int = 40,
                  precision: int = 8,
                  num_hashes: int = 128,
                  sketch_type: str = "minimizer",
                  chunk_size: int = 1000,
                  verbose: bool = False) -> Optional['SequenceSketch']:
        """Create a SequenceSketch from a FASTA/FASTQ file.
        
        Args:
            filename: Path to FASTA/FASTQ file
            kmer_size: Size of kmers
            window_size: Size of sliding window
            precision: Precision for HyperLogLog
            num_hashes: Number of hash functions for MinHash
            sketch_type: Type of sketch to use
            chunk_size: Number of sequences to process at once
            verbose: Whether to print progress
            
        Returns:
            SequenceSketch object or None if file processing fails
        """
        sketch = cls(
            kmer_size=kmer_size,
            window_size=window_size,
            precision=precision,
            num_hashes=num_hashes,
            sketch_type=sketch_type
        )
        
        try:
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
    """Read sequences from a FASTA/FASTQ file in chunks.
    
    Args:
        filename: Path to FASTA/FASTQ file
        chunk_size: Number of sequences to read at a time
        
    Returns:
        Iterator yielding lists of SeqRecord objects, with each list containing
        up to chunk_size sequences
    """
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