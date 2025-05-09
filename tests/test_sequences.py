from __future__ import annotations
import pytest # type: ignore
from hammock.lib.sequences import SequenceSketch
from Bio import SeqIO # type: ignore
from io import StringIO

class TestSequenceSketchQuick:
    """Quick tests for sequence sketching functionality."""
    
    def test_initialization(self):
        """Test that SequenceSketch initializes correctly."""
        sketch = SequenceSketch(sketch_type="hyperloglog", kmer_size=3)
        assert sketch.kmer_size == 3
        assert sketch.sketch_type == "hyperloglog"
        
        # Test invalid sketch type
        with pytest.raises(ValueError):
            SequenceSketch(sketch_type="invalid", kmer_size=3)

class TestSequenceSketchFull:
    """Full test suite for sequence sketching functionality."""
    
    def test_custom_initialization(self):
        """Test custom initialization parameters."""
        sketch = SequenceSketch(
            sketch_type="hyperloglog",
            kmer_size=10,
            precision=12
        )
        assert sketch.kmer_size == 10
    
    def test_sequence_similarity(self):
        """Test similarity estimation for sequence sketches."""
        sketch1 = SequenceSketch(sketch_type="hyperloglog", kmer_size=3)
        sketch2 = SequenceSketch(sketch_type="hyperloglog", kmer_size=3)
        
        sketch1.add_string("ACGTACGT")
        sketch2.add_string("ACGTACGT")
        
        sim = sketch1.similarity_values(sketch2)
        assert sim['jaccard_similarity'] == 1.0

def test_sequence_sketch_initialization():
    sketch = SequenceSketch(sketch_type="hyperloglog", kmer_size=3)
    assert sketch.kmer_size == 3
    assert sketch.sketch_type == "hyperloglog"

def test_invalid_sketch_type():
    with pytest.raises(ValueError):
        SequenceSketch(sketch_type="invalid", kmer_size=3)

def test_add_string():
    sketch = SequenceSketch(sketch_type="hyperloglog", kmer_size=3)
    sketch.add_string("ACGTACGT")
    
    # Add same sequence to another sketch
    sketch2 = SequenceSketch(sketch_type="hyperloglog", kmer_size=3)
    sketch2.add_string("ACGTACGT")
    
    sim = sketch.similarity_values(sketch2)
    assert sim['jaccard_similarity'] == 1.0

def test_sequence_similarity():
    """Test sequence similarity calculation."""
    sketch1 = SequenceSketch(sketch_type="minimizer", kmer_size=4)
    sketch2 = SequenceSketch(sketch_type="minimizer", kmer_size=4)
    
    sketch1.add_string("ACGTACGT")
    sketch2.add_string("ACGTACGT")
    
    result = sketch1.similarity_values(sketch2)
    assert 'jaccard_similarity' in result
    assert result['jaccard_similarity'] == 1.0
    assert 'gap_similarity' in result
    assert result['gap_similarity'] == 1.0

def test_different_sequences():
    sketch1 = SequenceSketch(sketch_type="hyperloglog", kmer_size=3)
    sketch2 = SequenceSketch(sketch_type="hyperloglog", kmer_size=3)
    
    sketch1.add_string("AAAAAAA")
    sketch2.add_string("TTTTTTT")
    
    sim = sketch1.similarity_values(sketch2)
    assert 0.0 <= sim['jaccard_similarity'] <= 1.0

@pytest.mark.quick
def test_empty_sequence():
    """Test handling of empty sequences."""
    sketch1 = SequenceSketch(sketch_type="hyperloglog", kmer_size=4)
    sketch2 = SequenceSketch(sketch_type="hyperloglog", kmer_size=4)
    
    sketch1.add_string("")
    sketch2.add_string("ACGT")
    
    # Use similarity_values
    result = sketch1.similarity_values(sketch2)
    assert result['jaccard_similarity'] == 0.0

@pytest.mark.quick
def test_sequence_shorter_than_kmer():
    """Test handling of sequences shorter than k-mer size."""
    sketch1 = SequenceSketch(sketch_type="hyperloglog", kmer_size=4)
    sketch2 = SequenceSketch(sketch_type="hyperloglog", kmer_size=4)
    
    sketch1.add_string("ACG")  # shorter than k=4
    sketch2.add_string("ACGT")
    
    # Use similarity_values
    result = sketch1.similarity_values(sketch2)
    assert result['jaccard_similarity'] == 0.0

def test_different_kmer_sizes():
    sketch1 = SequenceSketch(sketch_type="hyperloglog", kmer_size=3)
    sketch2 = SequenceSketch(sketch_type="hyperloglog", kmer_size=4)
    
    with pytest.raises(ValueError):
        sketch1.similarity_values(sketch2) 