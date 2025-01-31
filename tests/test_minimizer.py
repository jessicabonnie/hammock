import pytest # type: ignore
from hammock.lib.minimizer import MinimizerSketch

def test_minimizer_sketch_initialization():
    sketch = MinimizerSketch(k=5, w=10, gapk=3)
    assert sketch.k == 5
    assert sketch.w == 10
    assert sketch.gapk == 3

def test_add_sequence():
    sketch = MinimizerSketch(k=4, w=6, gapk=2)
    # Simple sequence with known minimizers
    sequence = "ACGTACGTACGT"
    sketch.add_string(sequence)
    
    # Add same sequence again - should not affect similarity
    sketch2 = MinimizerSketch(k=4, w=6, gapk=2)
    sketch2.add_string(sequence)
    
    sim = sketch.estimate_similarity(sketch2)
    assert sim['hash_similarity'] == 1.0
    assert sim['gap_similarity'] == 1.0
    assert sim['combined_similarity'] == 1.0

def test_compare_different_sequences():
    sketch1 = MinimizerSketch(k=4, w=6, gapk=2)
    sketch2 = MinimizerSketch(k=4, w=6, gapk=2)
    
    sketch1.add_string("ACGTACGTACGT")
    sketch2.add_string("TGCATGCATGCA")
    
    sim = sketch1.estimate_similarity(sketch2)
    assert 0 <= sim['hash_similarity'] <= 1
    assert 0 <= sim['gap_similarity'] <= 1
    assert 0 <= sim['combined_similarity'] <= 1

def test_empty_sequence():
    sketch = MinimizerSketch(k=4, w=6, gapk=2)
    sketch.add_string("")  # empty sequence
    
    sketch2 = MinimizerSketch(k=4, w=6, gapk=2)
    sim = sketch.estimate_similarity(sketch2)
    # Both sketches are empty, so similarity should be 0
    assert sim['hash_similarity'] == 0
    assert sim['gap_similarity'] == 0
    assert sim['combined_similarity'] == 0

def test_sequence_shorter_than_k():
    sketch = MinimizerSketch(k=5, w=10, gapk=2)
    sketch.add_string("ACG")  # shorter than k
    
    sketch2 = MinimizerSketch(k=5, w=10, gapk=2)
    sim = sketch.estimate_similarity(sketch2)
    # Sequence too short to generate minimizers
    assert sim['hash_similarity'] == 0
    assert sim['gap_similarity'] == 0
    assert sim['combined_similarity'] == 0

def test_different_parameters():
    sketch1 = MinimizerSketch(k=4, w=6, gapk=2)
    sketch2 = MinimizerSketch(k=5, w=6, gapk=2)
    
    sketch1.add_string("ACGTACGTACGT")
    sketch2.add_string("ACGTACGTACGT")
    
    with pytest.raises(ValueError):
        sketch1.estimate_similarity(sketch2) 

class TestMinimizerSketchQuick:
    """Quick tests for minimizer sketching functionality."""
    
    def test_initialization(self):
        """Test that MinimizerSketch initializes correctly."""
        sketch = MinimizerSketch(k=5, w=10, gapk=3)
        assert sketch.k == 5
        assert sketch.w == 10
        assert sketch.gapk == 3
    
    def test_basic_similarity(self):
        """Test basic similarity calculation."""
        sketch1 = MinimizerSketch(k=4, w=6, gapk=2)
        sketch2 = MinimizerSketch(k=4, w=6, gapk=2)
        
        sketch1.add_string("ACGTACGT")
        sketch2.add_string("ACGTACGT")
        
        sim = sketch1.estimate_similarity(sketch2)
        assert sim['combined_similarity'] == 1.0 