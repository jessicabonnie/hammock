import pytest # type: ignore
from hammock.lib.minimizer import MinimizerSketch

def test_minimizer_sketch_initialization():
    sketch = MinimizerSketch(kmer_size=5, window_size=10)
    assert sketch.kmer_size == 5
    assert sketch.window_size == 10

def test_add_string():
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6)
    sequence = "ACGTACGTATTAGATCCG"
    sketch1.add_string(sequence)
    
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6)
    sketch2.add_string(sequence)
    
    sim = sketch1.similarity_values(sketch2)
    # assert abs(1.0 - sim['hash_similarity']) < 0.1
    # assert abs(1.0 - sim['hash_with_ends_similarity']) < 0.1
    assert abs(1.0 - sim['jaccard_similarity']) < 0.1
    assert abs(1.0 - sim['jaccard_similarity_with_ends']) < 0.1

def test_compare_different_sequences():
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6)
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6)
    
    sketch1.add_string("ACGTACGTACGT")
    sketch2.add_string("TGCATGCATGCA")
    
    sim = sketch1.similarity_values(sketch2)
    # assert 0 <= sim['hash_similarity'] <= 1
    # assert 0 <= sim['hash_with_ends_similarity'] <= 1
    assert 0 <= sim['jaccard_similarity'] <= 1
    assert 0 <= sim['jaccard_similarity_with_ends'] <= 1

def test_empty_sequence():
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6)
    sketch1.add_string("")  # empty sequence
    
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6)
    sim = sketch1.similarity_values(sketch2)
    # assert sim['hash_similarity'] == 0
    # assert sim['hash_with_ends_similarity'] == 0
    assert sim['jaccard_similarity'] == 0
    assert sim['jaccard_similarity_with_ends'] == 0

def test_sequence_shorter_than_k():
    sketch1 = MinimizerSketch(kmer_size=5, window_size=10)
    sketch1.add_string("ACG")  # shorter than k
    
    sketch2 = MinimizerSketch(kmer_size=5, window_size=10)
    sim = sketch1.similarity_values(sketch2)
    # assert sim['hash_similarity'] == 0
    # assert sim['hash_with_ends_similarity'] == 0
    assert sim['jaccard_similarity'] == 0
    assert sim['jaccard_similarity_with_ends'] == 0

def test_different_parameters():
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6)
    sketch2 = MinimizerSketch(kmer_size=5, window_size=6)
    
    sketch1.add_string("ACGTACGTACGT")
    sketch2.add_string("ACGTACGTACGT")
    
    with pytest.raises(ValueError):
        sketch1.similarity_values(sketch2)

def test_end_kmers():
    sketch = MinimizerSketch(kmer_size=4, window_size=6)
    sequence = "ACGTACGTACGT"
    sketch.add_string(sequence)
    
    # Create sketch with same start/end but different middle
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6)
    modified_sequence = "ACGT" + "TTTT" + sequence[-4:]  # Same ends, different middle
    sketch2.add_string(modified_sequence)
    
    sim = sketch.similarity_values(sketch2)
    assert sim['jaccard_similarity_with_ends'] > sim['jaccard_similarity']  # End k-mer similarity should increase the score

class TestMinimizerSketchQuick:
    """Quick tests for minimizer sketching functionality."""
    
    def test_initialization(self):
        """Test that MinimizerSketch initializes correctly."""
        sketch = MinimizerSketch(kmer_size=5, window_size=10)
        assert sketch.kmer_size == 5
        assert sketch.window_size == 10
    
    def test_basic_similarity(self):
        """Test basic similarity calculation."""
        sketch1 = MinimizerSketch(kmer_size=4, window_size=5)
        sketch2 = MinimizerSketch(kmer_size=4, window_size=5)
        
        sketch1.add_string("ACGTACGTAAGG")
        sketch2.add_string("ACGTACGTAAGG")
        
        sim = sketch1.similarity_values(sketch2)
        assert sim['jaccard_similarity'] == 1.0

def test_basic_similarity():
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6)
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6)
    
    # Identical sequences should have similarity close to 1
    sequence = "ACCGTACTCGTACGTAA"
    sketch1.add_string(sequence)
    sketch2.add_string(sequence)
    
    sim = sketch1.similarity_values(sketch2)
    assert abs(1.0 - sim['jaccard_similarity']) < 0.1 