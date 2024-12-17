#!/usr/bin/env python
from Bio import SeqIO # type: ignore
from Digest import window_minimizer # type: ignore
from hammock.lib.sketchclass import Sketch
from typing import Optional

class MinimizerSketch(Sketch):
    """Sketch class for sequence data using minimizers."""
    
    def __init__(self, 
                 window_size: int = 40, 
                 kmer_size: int = 8,
                 precision: int = 8,
                 num_hashes: int = 128,
                 seed: int = 0):
        """Initialize MinimizerSketch.
        
        Args:
            window_size: Size of sliding window for minimizer selection
            kmer_size: Size of kmers to use
            precision: Precision for HyperLogLog sketching
            num_hashes: Number of hash functions for MinHash sketching
            seed: Random seed for hash functions
        """
        super().__init__(
            sketch_type="hyperloglog",
            precision=precision,
            num_hashes=num_hashes,
            kmer_size=kmer_size,
            seed=seed
        )
        self.window_size = window_size
        self.kmer_size = kmer_size
    
    # def _hash64(self, x: int) -> int:
    #     """64-bit hash function.
        
    #     Args:
    #         x: Integer value to hash
            
    #     Returns:
    #         64-bit hash value as integer
    #     """
    #     hasher = xxhash.xxh64(seed=self.seed)
    #     hasher.update(x.to_bytes(8, byteorder='little'))
    #     return hasher.intdigest()
        
    def add_sequence(self, sequence: str) -> None:
        """Add a sequence to the sketch using minimizers.
        
        Args:
            sequence: DNA/RNA sequence string
        """
        minimizers = window_minimizer(sequence, 
                                    w=self.window_size, 
                                    k=self.kmer_size, 
                                    include_hash=True)
        for _, hash_val in minimizers:
            # Add the minimizer hash value to the underlying HyperLogLog sketch
            self.add_int(hash_val)
        # add concatenated flanking kmers from start and end of sequence
        self.add_string(sequence[:self.kmer_size]+sequence[-self.kmer_size:])
            
    @classmethod
    def from_file(cls, 
                  filename: str, 
                  window_size: int = 40, 
                  kmer_size: int = 8,
                  precision: int = 8,
                  num_hashes: int = 128,
                  seed: int = 0,
                  chunk_size: int = 1000,
                  verbose: bool = False) -> Optional['MinimizerSketch']:
        """Create a MinimizerSketch from a FASTA/FASTQ file. 
        
        This function will read in a FASTA/FASTQ file and add the sequences to the sketch.
        The sequences are processed in chunks of `chunk_size` to avoid memory issues.
        
        Args:
            filename: Path to FASTA/FASTQ file
            window_size: Size of sliding window
            kmer_size: Size of kmers
            precision: Precision for HyperLogLog
            num_hashes: Number of hash functions for MinHash
            seed: Random seed
            chunk_size: Number of sequences to process at once
            verbose: Whether to print progress
            
        Returns:
            MinimizerSketch object or None if file processing fails
        """
        sketch = cls(
            window_size=window_size,
            kmer_size=kmer_size,
            precision=precision,
            num_hashes=num_hashes,
            seed=seed
        )
        
        try:
            formatx = "fasta" if filename.endswith((".fa", ".fasta")) else "fastq"
            with open(filename, "r") as file:
                records = []
                for record in SeqIO.parse(file, formatx):
                    records.append(record)
                    if len(records) >= chunk_size:
                        for rec in records:
                            sketch.add_sequence(str(rec.seq))
                        records = []
                        if verbose:
                            print(f"Processed {chunk_size} sequences from {filename}")
                # Process remaining records
                for rec in records:
                    sketch.add_sequence(str(rec.seq))
            return sketch
            
        except Exception as e:
            print(f"Error processing file {filename}: {str(e)}")
            return None