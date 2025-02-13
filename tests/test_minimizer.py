import pytest # type: ignore
from hammock.lib.minimizer import MinimizerSketch

def test_minimizer_sketch_initialization():
    sketch = MinimizerSketch(kmer_size=5, window_size=10, gapk=3)
    assert sketch.kmer_size == 5
    assert sketch.window_size == 10
    assert sketch.gapk == 3

def test_add_string():
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6, gapk=2)
    sequence = "ACGTACGTATTAGATCCG"
    sketch1.add_string(sequence)
    
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6, gapk=2)
    sketch2.add_string(sequence)
    
    hash_sim, hash_ends_sim, gap_sim, jaccard_sim = sketch1.compare_overlaps(sketch2)
    assert abs(1.0 - hash_sim) < 0.1
    assert abs(1.0 - hash_ends_sim) < 0.1
    assert abs(1.0 - gap_sim) < 0.1
    assert abs(1.0 - jaccard_sim) < 0.1

def test_compare_different_sequences():
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6, gapk=2)
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6, gapk=2)
    
    sketch1.add_string("ACGTACGTACGT")
    sketch2.add_string("TGCATGCATGCA")
    
    hash_sim, hash_ends_sim, gap_sim, jaccard_sim = sketch1.compare_overlaps(sketch2)
    assert 0 <= hash_sim <= 1
    assert 0 <= hash_ends_sim <= 1
    assert 0 <= gap_sim <= 1
    assert 0 <= jaccard_sim <= 1

def test_empty_sequence():
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6, gapk=2)
    sketch1.add_string("")  # empty sequence
    
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6, gapk=2)
    hash_sim, hash_ends_sim, gap_sim, jaccard_sim = sketch1.compare_overlaps(sketch2)
    assert hash_sim == 0
    assert hash_ends_sim == 0
    assert gap_sim == 0
    assert jaccard_sim == 0

def test_sequence_shorter_than_k():
    sketch1 = MinimizerSketch(kmer_size=5, window_size=10, gapk=2)
    sketch1.add_string("ACG")  # shorter than k
    
    sketch2 = MinimizerSketch(kmer_size=5, window_size=10, gapk=2)
    hash_sim, hash_ends_sim, gap_sim, jaccard_sim = sketch1.compare_overlaps(sketch2)
    assert hash_sim == 0
    assert hash_ends_sim == 0
    assert gap_sim == 0
    assert jaccard_sim == 0

def test_different_parameters():
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6, gapk=2)
    sketch2 = MinimizerSketch(kmer_size=5, window_size=6, gapk=2)
    
    sketch1.add_string("ACGTACGTACGT")
    sketch2.add_string("ACGTACGTACGT")
    
    with pytest.raises(ValueError):
        sketch1.compare_overlaps(sketch2)

def test_gap_patterns():
    sketch1 = MinimizerSketch(kmer_size=4, window_size=5, gapk=2)
    sequence = "ACGTAAAGTACGTAAGG"
    sketch1.add_string(sequence)
    
    sketch2 = MinimizerSketch(kmer_size=4, window_size=5, gapk=2)
    sketch2.add_string(sequence)
    
    hash_sim, hash_ends_sim, gap_sim, jaccard_sim = sketch1.compare_overlaps(sketch2)
    assert abs(1.0 - gap_sim) < 0.1

def test_end_kmers():
    sketch = MinimizerSketch(kmer_size=4, window_size=6, gapk=2)
    sequence = "ACGTACGTACGT"
    sketch.add_string(sequence)
    
    # Create sketch with same start/end but different middle
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6, gapk=2)
    modified_sequence = "ACGT" + "TTTT" + sequence[-4:]  # Same ends, different middle
    sketch2.add_string(modified_sequence)
    
    hash_sim, hash_ends_sim, gap_sim, jaccard_sim = sketch.compare_overlaps(sketch2)
    assert hash_ends_sim > hash_sim  # End k-mer similarity should increase the score

class TestMinimizerSketchQuick:
    """Quick tests for minimizer sketching functionality."""
    
    def test_initialization(self):
        """Test that MinimizerSketch initializes correctly."""
        sketch = MinimizerSketch(kmer_size=5, window_size=10, gapk=3)
        assert sketch.kmer_size == 5
        assert sketch.window_size == 10
        assert sketch.gapk == 3
    
    def test_basic_similarity(self):
        """Test basic similarity calculation."""
        sketch1 = MinimizerSketch(kmer_size=4, window_size=5, gapk=2)
        sketch2 = MinimizerSketch(kmer_size=4, window_size=5, gapk=2)
        
        sketch1.add_string("ACGTACGTAAGG")
        sketch2.add_string("ACGTACGTAAGG")
        
        sim = sketch1.similarity_values(sketch2)
        assert sim['jaccard_similarity'] == 1.0

def test_basic_similarity():
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6, gapk=2)
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6, gapk=2)
    
    # Identical sequences should have similarity close to 1
    sequence = "ACCGTACTCGTACGTAA"
    sketch1.add_string(sequence)
    sketch2.add_string(sequence)
    
    hash_sim, hash_ends_sim, gap_sim, jaccard_sim = sketch1.compare_overlaps(sketch2)
    assert abs(1.0 - hash_sim) < 0.1 