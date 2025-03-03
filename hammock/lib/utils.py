import gzip
from typing import Iterator, List, Optional
from Bio import SeqIO

def read_sequences(filename: str, chunk_size: int = 1000) -> Iterator[List[SeqIO.SeqRecord]]:
    """Read sequences from a FASTA/FASTQ file in chunks.
    
    FASTA format:
    >sequence_identifier [optional description]
    ACGTACGTACGT...
    ACGTACGTACGT... (sequence can span multiple lines)
    
    FASTQ format:
    @sequence_identifier [optional description]
    ACGTACGTACGT...
    +[optional repeat of identifier]
    !@#$%^&*... (quality scores)
    """
    # Expanded list of FASTA file extensions
    fasta_extensions = (".fa", ".fasta", ".fna", ".ffn", ".frn", ".faa")
    fastq_extensions = (".fq", ".fastq")
    
    # Check if file is gzipped
    is_gzipped = filename.endswith(".gz")
    base_filename = filename[:-3] if is_gzipped else filename
    
    # Determine format
    formatx = "fasta" if any(base_filename.endswith(ext) for ext in fasta_extensions) else "fastq"
    
    records = []
    
    # Open with appropriate method
    opener = gzip.open if is_gzipped else open
    mode = "rt" if is_gzipped else "r"  # text mode for gzip
    
    try:
        with opener(filename, mode) as file:
            # Use BioPython's SeqIO.parse which correctly handles multi-line FASTA/FASTQ
            for record in SeqIO.parse(file, formatx):
                records.append(record)
                if len(records) >= chunk_size:
                    yield records
                    records = []
            if records:  # Yield any remaining records
                yield records
    except Exception as e:
        print(f"Error reading file {filename}: {str(e)}")
        yield []  # Return empty list on error 