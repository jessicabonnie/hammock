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

def test_sequence_sketch_use_sets_parameter():
    """Test that use_sets parameter is properly passed to MinimizerSketch."""
    # Test default behavior (use_sets=False)
    sketch1 = SequenceSketch(sketch_type="minimizer", kmer_size=4, use_sets=False)
    assert sketch1.use_sets == False
    assert sketch1.sketch.use_sets == False
    assert not hasattr(sketch1.sketch, 'minimizers')
    
    # Test with use_sets=True
    sketch2 = SequenceSketch(sketch_type="minimizer", kmer_size=4, use_sets=True)
    assert sketch2.use_sets == True
    assert sketch2.sketch.use_sets == True
    assert hasattr(sketch2.sketch, 'minimizers')
    assert hasattr(sketch2.sketch, 'startend_kmers')

def test_sequence_sketch_use_sets_only_affects_minimizer():
    """Test that use_sets parameter only affects minimizer sketch type."""
    # use_sets should not affect non-minimizer sketch types
    sketch_hll = SequenceSketch(sketch_type="hyperloglog", kmer_size=4, use_sets=True)
    assert sketch_hll.use_sets == True
    assert not hasattr(sketch_hll.sketch, 'use_sets')  # HyperLogLog doesn't have use_sets
    
    sketch_mh = SequenceSketch(sketch_type="minhash", kmer_size=4, use_sets=True)
    assert sketch_mh.use_sets == True
    assert not hasattr(sketch_mh.sketch, 'use_sets')  # MinHash doesn't have use_sets

def test_sequence_sketch_use_sets_similarity_metrics():
    """Test that SequenceSketch returns correct similarity metrics with use_sets=True."""
    sketch1 = SequenceSketch(sketch_type="minimizer", kmer_size=4, use_sets=True)
    sketch2 = SequenceSketch(sketch_type="minimizer", kmer_size=4, use_sets=True)
    
    sequence = "ACGTACGTACGTACGT"
    sketch1.add_string(sequence)
    sketch2.add_string(sequence)
    
    sim = sketch1.similarity_values(sketch2)
    
    # Should have both sketch-based and set-based metrics
    assert 'jaccard_similarity' in sim
    assert 'jaccard_similarity_with_ends' in sim
    assert 'set_similarity' in sim
    assert 'set_with_ends_similarity' in sim
    
    # For identical sequences, set-based should be perfect
    assert sim['set_similarity'] == 1.0
    assert sim['set_with_ends_similarity'] == 1.0

def test_sequence_sketch_from_file_use_sets(tmp_path):
    """Test that from_file method supports use_sets parameter."""
    # Create a simple FASTA file
    fasta_content = ">seq1\nACGTACGTACGTACGT\n>seq2\nTGCATGCATGCATGCA\n"
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)
    
    # Test with use_sets=True
    sketch = SequenceSketch.from_file(str(fasta_file), 
                                      sketch_type="minimizer", 
                                      kmer_size=4, 
                                      use_sets=True)
    
    assert sketch is not None
    assert sketch.sketch_type == "minimizer"
    assert sketch.use_sets == True
    assert sketch.sketch.use_sets == True
    assert hasattr(sketch.sketch, 'minimizers')
    assert hasattr(sketch.sketch, 'startend_kmers')
    
    # Should have non-empty sets after processing sequences
    assert len(sketch.sketch.minimizers) > 0
    assert len(sketch.sketch.startend_kmers) > 0

def test_sequence_sketch_cardinality_with_sets():
    """Test that cardinality estimation uses sets when use_sets=True."""
    sequence = "ACGTACGTACGTACGT"
    
    # Test with use_sets=False (sketch-based)
    sketch_no_sets = SequenceSketch(sketch_type="minimizer", kmer_size=4, use_sets=False)
    sketch_no_sets.add_string(sequence)
    cardinality_sketch = sketch_no_sets.estimate_cardinality()
    
    # Test with use_sets=True (exact count)
    sketch_with_sets = SequenceSketch(sketch_type="minimizer", kmer_size=4, use_sets=True)
    sketch_with_sets.add_string(sequence)
    cardinality_exact = sketch_with_sets.estimate_cardinality()
    
    # Both should be positive
    assert cardinality_sketch > 0
    assert cardinality_exact > 0
    
    # Exact count should equal the number of unique minimizers
    assert cardinality_exact == float(len(sketch_with_sets.sketch.minimizers)) 