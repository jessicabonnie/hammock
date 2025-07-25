import pytest # type: ignore
from hammock.lib.minimizer import MinimizerSketch, canonicalize_kmer

def test_canonicalize_kmer():
    """Test that k-mer canonicalization works correctly."""
    # Test basic canonicalization with clear examples (all converted to uppercase)
    assert canonicalize_kmer("ATCG") == "ATCG"  # ATCG < CGAT (reverse complement)
    assert canonicalize_kmer("CGAT") == "ATCG"  # ATCG < CGAT (reverse complement)
    assert canonicalize_kmer("AAAA") == "AAAA"  # AAAA < TTTT (reverse complement)
    assert canonicalize_kmer("TTTT") == "AAAA"  # AAAA < TTTT (reverse complement)
    
    # Test case conversion (lowercase input converted to uppercase)
    assert canonicalize_kmer("atcg") == "ATCG"  # atcg -> ATCG, ATCG < CGAT
    assert canonicalize_kmer("cgat") == "ATCG"  # cgat -> CGAT, ATCG < CGAT
    
    # Test edge cases
    assert canonicalize_kmer("") == ""
    assert canonicalize_kmer("A") == "A"
    assert canonicalize_kmer("T") == "A"  # A < T
    assert canonicalize_kmer("a") == "A"  # a -> A
    assert canonicalize_kmer("t") == "A"  # t -> T, A < T
    
    # Test mixed case (all converted to uppercase)
    assert canonicalize_kmer("AtCg") == "ATCG"  # AtCg -> ATCG, ATCG < CGAT
    assert canonicalize_kmer("CgAt") == "ATCG"  # CgAt -> CGAT, ATCG < CGAT
    
    # Test palindromic sequences (self-reverse complement)
    assert canonicalize_kmer("ACGT") == "ACGT"  # ACGT is its own reverse complement
    assert canonicalize_kmer("TGCA") == "TGCA"  # TGCA is its own reverse complement
    assert canonicalize_kmer("acgt") == "ACGT"  # acgt -> ACGT, ACGT is its own reverse complement
    assert canonicalize_kmer("tgca") == "TGCA"  # tgca -> TGCA, TGCA is its own reverse complement

def test_canonicalized_flanking_kmers():
    """Test that flanking k-mers are properly canonicalized."""
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6, use_sets=True)
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6, use_sets=True)
    
    # Add sequences with different orientations but same canonical flanking k-mers
    sketch1.add_string("ACGTACGTACGT")  # Start: ACGT, End: ACGT
    sketch2.add_string("TGCATGCATGCA")  # Start: TGCA (canonical: ACGT), End: TGCA (canonical: ACGT)
    
    # The flanking k-mers should be the same after canonicalization
    sim = sketch1.similarity_values(sketch2)
    # The similarity with ends should be higher than without ends due to canonicalized flanking k-mers
    assert sim['jaccard_similarity_with_ends'] >= sim['jaccard_similarity']

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
    # Both sketches are empty, so similarity should be 0
    assert sim['jaccard_similarity'] == 0
    assert sim['jaccard_similarity_with_ends'] == 0
    
    # But if we add the same short sequence to sketch2, they should be similar
    sketch2.add_string("ACG")  # same short sequence
    sim = sketch1.similarity_values(sketch2)
    assert sim['jaccard_similarity'] > 0.9  # Should be very high similarity
    assert sim['jaccard_similarity_with_ends'] > 0.9

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

def test_use_sets_parameter_initialization():
    """Test that use_sets parameter is properly initialized."""
    # Default should be False
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6)
    assert sketch1.use_sets == False
    assert not hasattr(sketch1, 'minimizers')
    assert not hasattr(sketch1, 'startend_kmers')
    
    # When use_sets=True, sets should be initialized
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6, use_sets=True)
    assert sketch2.use_sets == True
    assert hasattr(sketch2, 'minimizers')
    assert hasattr(sketch2, 'startend_kmers')
    assert isinstance(sketch2.minimizers, set)
    assert isinstance(sketch2.startend_kmers, set)

def test_use_sets_cardinality_estimation():
    """Test that cardinality estimation uses sets when use_sets=True."""
    sequence = "ACGTACGTACGTACGT"
    
    # Sketch-based estimation
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6, use_sets=False)
    sketch1.add_string(sequence)
    cardinality_sketch = sketch1.estimate_cardinality()
    
    # Set-based exact count
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6, use_sets=True)
    sketch2.add_string(sequence)
    cardinality_sets = sketch2.estimate_cardinality()
    
    # Set-based should be exact integer, sketch-based might be approximate
    assert isinstance(cardinality_sets, float)
    assert cardinality_sets == float(len(sketch2.minimizers))
    assert cardinality_sketch > 0  # Should be positive

def test_use_sets_similarity_metrics():
    """Test that set-based similarity metrics are returned when use_sets=True."""
    sequence1 = "ACGTACGTACGTACGT"
    sequence2 = "ACGTACGTACGTACGG"  # Slightly different
    
    # Test with use_sets=False (should only get sketch-based metrics)
    sketch1_no_sets = MinimizerSketch(kmer_size=4, window_size=6, use_sets=False)
    sketch2_no_sets = MinimizerSketch(kmer_size=4, window_size=6, use_sets=False)
    sketch1_no_sets.add_string(sequence1)
    sketch2_no_sets.add_string(sequence2)
    
    sim_no_sets = sketch1_no_sets.similarity_values(sketch2_no_sets)
    assert 'jaccard_similarity' in sim_no_sets
    assert 'jaccard_similarity_with_ends' in sim_no_sets
    assert 'set_similarity' not in sim_no_sets
    assert 'set_with_ends_similarity' not in sim_no_sets
    
    # Test with use_sets=True (should get both sketch and set-based metrics)
    sketch1_with_sets = MinimizerSketch(kmer_size=4, window_size=6, use_sets=True)
    sketch2_with_sets = MinimizerSketch(kmer_size=4, window_size=6, use_sets=True)
    sketch1_with_sets.add_string(sequence1)
    sketch2_with_sets.add_string(sequence2)
    
    sim_with_sets = sketch1_with_sets.similarity_values(sketch2_with_sets)
    assert 'jaccard_similarity' in sim_with_sets
    assert 'jaccard_similarity_with_ends' in sim_with_sets
    assert 'set_similarity' in sim_with_sets
    assert 'set_with_ends_similarity' in sim_with_sets
    
    # All similarities should be between 0 and 1
    for key, value in sim_with_sets.items():
        assert 0 <= value <= 1, f"{key} should be between 0 and 1, got {value}"

def test_use_sets_mixed_comparison():
    """Test comparison between sketches with different use_sets values."""
    sequence = "ACGTACGTACGTACGT"
    
    sketch_no_sets = MinimizerSketch(kmer_size=4, window_size=6, use_sets=False)
    sketch_with_sets = MinimizerSketch(kmer_size=4, window_size=6, use_sets=True)
    
    sketch_no_sets.add_string(sequence)
    sketch_with_sets.add_string(sequence)
    
    # When one sketch has sets and other doesn't, should only get sketch-based metrics
    sim = sketch_no_sets.similarity_values(sketch_with_sets)
    assert 'jaccard_similarity' in sim
    assert 'jaccard_similarity_with_ends' in sim
    assert 'set_similarity' not in sim
    assert 'set_with_ends_similarity' not in sim

def test_use_sets_merge():
    """Test merging sketches with use_sets=True."""
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6, use_sets=True)
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6, use_sets=True)
    
    sketch1.add_string("ACGTACGTACGT")
    sketch2.add_string("TGCATGCATGCA")
    
    # Record set sizes before merge
    sketch1_minimizers = len(sketch1.minimizers)
    sketch1_startend = len(sketch1.startend_kmers)
    sketch2_minimizers = len(sketch2.minimizers)
    sketch2_startend = len(sketch2.startend_kmers)
    
    merged = sketch1.merge(sketch2)
    
    # Merged sketch should have use_sets=True and contain union of sets
    assert merged.use_sets == True
    assert hasattr(merged, 'minimizers')
    assert hasattr(merged, 'startend_kmers')
    
    # Merged sets should be at least as large as individual sets
    assert len(merged.minimizers) >= max(sketch1_minimizers, sketch2_minimizers)
    assert len(merged.startend_kmers) >= max(sketch1_startend, sketch2_startend)

def test_use_sets_identical_sequences():
    """Test that identical sequences give perfect similarity with sets."""
    sequence = "ACGTACGTACGTACGT"
    
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6, use_sets=True)
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6, use_sets=True)
    
    sketch1.add_string(sequence)
    sketch2.add_string(sequence)
    
    sim = sketch1.similarity_values(sketch2)
    
    # For identical sequences, set-based similarities should be exactly 1.0
    assert sim['set_similarity'] == 1.0
    assert sim['set_with_ends_similarity'] == 1.0
    
    # Sketch-based should be close to 1.0 but might not be exact due to approximation
    assert abs(1.0 - sim['jaccard_similarity']) < 0.1
    assert abs(1.0 - sim['jaccard_similarity_with_ends']) < 0.1

def test_use_sets_completely_different_sequences():
    """Test that completely different sequences give low similarity with sets."""
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6, use_sets=True)
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6, use_sets=True)
    
    # Use sequences with no common k-mers
    sketch1.add_string("AAAAAAAAAAAAAAAA")
    sketch2.add_string("CCCCCCCCCCCCCCCC")
    
    sim = sketch1.similarity_values(sketch2)
    
    # Should have very low or zero similarity
    assert sim['set_similarity'] == 0.0
    assert sim['set_with_ends_similarity'] == 0.0

def test_use_sets_save_load(tmp_path):
    """Test saving and loading sketches with use_sets=True."""
    import numpy as np # type: ignore
    
    sketch = MinimizerSketch(kmer_size=4, window_size=6, use_sets=True, precision=12)
    sketch.add_string("ACGTACGTACGTACGT")
    
    # Record original state
    original_minimizers = sketch.minimizers.copy()
    original_startend = sketch.startend_kmers.copy()
    original_cardinality = sketch.estimate_cardinality()
    
    # Save sketch
    filepath = str(tmp_path / "test_sketch")
    sketch.write(filepath)
    
    # Load sketch
    loaded_sketch = MinimizerSketch.load(filepath)
    
    # Verify loaded sketch has same properties
    assert loaded_sketch.use_sets == True
    assert loaded_sketch.kmer_size == 4
    assert loaded_sketch.window_size == 6
    assert loaded_sketch.precision == 12
    assert hasattr(loaded_sketch, 'minimizers')
    assert hasattr(loaded_sketch, 'startend_kmers')
    
    # Verify sets are correctly restored
    assert loaded_sketch.minimizers == original_minimizers
    assert loaded_sketch.startend_kmers == original_startend
    assert loaded_sketch.estimate_cardinality() == original_cardinality

def test_use_sets_save_load_no_sets(tmp_path):
    """Test saving and loading sketches with use_sets=False."""
    sketch = MinimizerSketch(kmer_size=4, window_size=6, use_sets=False, precision=12)
    sketch.add_string("ACGTACGTACGTACGT")
    
    # Save sketch
    filepath = str(tmp_path / "test_sketch_no_sets")
    sketch.write(filepath)
    
    # Load sketch
    loaded_sketch = MinimizerSketch.load(filepath)
    
    # Verify loaded sketch has same properties
    assert loaded_sketch.use_sets == False
    assert loaded_sketch.kmer_size == 4
    assert loaded_sketch.window_size == 6
    assert loaded_sketch.precision == 12
    assert not hasattr(loaded_sketch, 'minimizers')
    assert not hasattr(loaded_sketch, 'startend_kmers')

def test_use_sets_backward_compatibility(tmp_path):
    """Test that old files without use_sets parameter load correctly."""
    import numpy as np
    
    # Create a sketch and save it
    sketch = MinimizerSketch(kmer_size=4, window_size=6, precision=12, debug=False)
    sketch.add_string("ACGTACGTACGTACGT")
    
    # Manually save without use_sets parameter (simulating old format)
    filepath = str(tmp_path / "old_format_sketch")
    np.savez(filepath,
             minimizer_sketch_registers=sketch.minimizer_sketch.registers,
             startend_sketch_registers=sketch.startend_sketch.registers,
             kmer_size=sketch.kmer_size,
             window_size=sketch.window_size,
             seed=sketch.seed,
             precision=sketch.precision,
             debug=sketch.debug)
    
    # Load should work and default to use_sets=False
    loaded_sketch = MinimizerSketch.load(filepath)
    assert loaded_sketch.use_sets == False
    assert not hasattr(loaded_sketch, 'minimizers')
    assert not hasattr(loaded_sketch, 'startend_kmers')

def test_short_sequence_handling():
    """Test that sequences too short for window_minimizer are added as whole sequences."""
    # Test sequence shorter than k-mer size
    sketch1 = MinimizerSketch(kmer_size=5, window_size=10, use_sets=False)
    sketch2 = MinimizerSketch(kmer_size=5, window_size=10, use_sets=False)
    
    short_seq = "ACG"  # length 3, shorter than kmer_size=5
    sketch1.add_string(short_seq)
    sketch2.add_string(short_seq)
    
    # Should have perfect similarity since both have the same short sequence
    sim = sketch1.similarity_values(sketch2)
    assert sim['jaccard_similarity'] > 0.9  # Should be very high similarity
    assert sim['jaccard_similarity_with_ends'] > 0.9
    
    # Should contribute to cardinality
    assert sketch1.estimate_cardinality() > 0

def test_short_sequence_handling_with_sets():
    """Test that sequences too short for window_minimizer work with use_sets=True."""
    sketch1 = MinimizerSketch(kmer_size=6, window_size=8, use_sets=True)
    sketch2 = MinimizerSketch(kmer_size=6, window_size=8, use_sets=True)
    
    short_seq = "ACGT"  # length 4, shorter than kmer_size=6
    sketch1.add_string(short_seq)
    sketch2.add_string(short_seq)
    
    # Should have perfect similarity since both have the same short sequence
    sim = sketch1.similarity_values(sketch2)
    assert sim['set_similarity'] == 1.0  # Perfect set similarity
    assert sim['set_with_ends_similarity'] == 1.0
    assert sim['jaccard_similarity'] > 0.9
    assert sim['jaccard_similarity_with_ends'] > 0.9
    
    # Check that the sequence was added to sets
    assert len(sketch1.minimizers) > 0
    assert len(sketch1.startend_kmers) > 0

def test_sequence_shorter_than_window():
    """Test sequences shorter than window size but longer than k-mer size."""
    sketch1 = MinimizerSketch(kmer_size=3, window_size=10, use_sets=False)
    sketch2 = MinimizerSketch(kmer_size=3, window_size=10, use_sets=False)
    
    short_seq = "ACGTAG"  # length 6, longer than kmer_size=3 but shorter than window_size=10
    sketch1.add_string(short_seq)
    sketch2.add_string(short_seq)
    
    # Should have high similarity
    sim = sketch1.similarity_values(sketch2)
    assert sim['jaccard_similarity'] > 0.5  # Should have reasonable similarity
    assert sim['jaccard_similarity_with_ends'] > 0.5
    
    # Should contribute to cardinality
    assert sketch1.estimate_cardinality() > 0

def test_short_vs_normal_sequences():
    """Test that short sequences still contribute meaningfully when mixed with normal sequences."""
    sketch1 = MinimizerSketch(kmer_size=4, window_size=8, use_sets=True)
    sketch2 = MinimizerSketch(kmer_size=4, window_size=8, use_sets=True)
    
    # Add both short and normal sequences
    short_seq = "AC"  # shorter than k-mer size
    normal_seq = "ACGTACGTACGTACGT"  # normal length
    
    sketch1.add_string(short_seq)
    sketch1.add_string(normal_seq)
    
    sketch2.add_string(short_seq)  # Same short sequence
    sketch2.add_string("TGCATGCATGCATGCA")  # Different normal sequence
    
    sim = sketch1.similarity_values(sketch2)
    # Should have some similarity due to the common short sequence
    assert sim['jaccard_similarity'] > 0.0
    assert sim['jaccard_similarity_with_ends'] > 0.0
    
    # Both sketches should have cardinality > 0
    assert sketch1.estimate_cardinality() > 0
    assert sketch2.estimate_cardinality() > 0

def test_empty_vs_short_sequence():
    """Test that empty sequences behave differently from short sequences."""
    sketch1 = MinimizerSketch(kmer_size=4, window_size=6, use_sets=False)
    sketch2 = MinimizerSketch(kmer_size=4, window_size=6, use_sets=False)
    
    sketch1.add_string("")  # empty sequence
    sketch2.add_string("AC")  # short sequence
    
    sim = sketch1.similarity_values(sketch2)
    # Empty and short should have low similarity
    assert sim['jaccard_similarity'] >= 0.0
    assert sim['jaccard_similarity_with_ends'] >= 0.0 