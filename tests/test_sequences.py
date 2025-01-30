from __future__ import annotations
import pytest # type: ignore
from hammock.lib.sequences import SequenceSketch
from Bio import SeqIO
from io import StringIO

class TestSequenceSketchQuick:
    """Quick tests for sequence sketching functionality."""
    
    def test_initialization(self):
        """Test that SequenceSketch initializes correctly."""
        # Test default initialization
        sketch = SequenceSketch()
        assert sketch.kmer_size == 8
        assert sketch.window_size == 40
        assert sketch.num_sequences == 0
        assert sketch.total_sequence_length == 0
        
        # Test invalid sketch type
        with pytest.raises(ValueError):
            SequenceSketch(sketch_type="invalid")

class TestSequenceSketchFull:
    """Full test suite for sequence sketching functionality."""
    
    def test_custom_initialization(self):
        """Test custom initialization parameters."""
        sketch = SequenceSketch(
            sketch_type="hyperloglog",
            kmer_size=10,
            window_size=50,
            precision=12,
            num_hashes=256,
            seed=42
        )
        assert sketch.kmer_size == 10
        assert sketch.window_size == 50
        assert sketch.num_sequences == 0
        assert sketch.total_sequence_length == 0
    
    def test_sequence_cardinality(self):
        """Test cardinality estimation for sequence sketches."""
        # Create test sequences with known cardinality
        sequences = [
            "ATCGATCGATCG",  # 12 kmers with k=4
            "GGGGGGGGGGGG",  # 9 kmers with k=4 (all G's)
            "ATATATATATAT"   # 9 kmers with k=4 (alternating A/T)
        ]
        
        # Test with different k-mer sizes
        k_sizes = [4, 6, 8]
        
        for k in k_sizes:
            # Create sketch
            sketch = SequenceSketch(kmer_size=k, precision=12)  # Higher precision for better accuracy
            
            # Calculate expected number of k-mers
            expected_kmers = 0
            for seq in sequences:
                if len(seq) >= k:
                    # Count unique k-mers in this sequence
                    kmers = set()
                    for i in range(len(seq) - k + 1):
                        kmers.add(seq[i:i+k])
                    expected_kmers += len(kmers)
            
            # Add sequences to sketch
            for seq in sequences:
                sketch.add_sequence(seq)
            
            # Get cardinality estimate
            estimated = sketch.estimate_cardinality()
            
            # Allow for some error in the estimate
            error = abs(estimated - expected_kmers) / expected_kmers
            assert error < 0.15, f"Cardinality estimate for k={k} too far off: " \
                               f"expected {expected_kmers}, got {estimated:.1f} (error: {error:.3f})" 