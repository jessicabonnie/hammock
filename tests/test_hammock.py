#!/usr/bin/env python
from __future__ import annotations
from hammock.lib.sequences import SequenceSketch
from hammock.lib.intervals import IntervalSketch
from unittest.mock import patch, mock_open
import pytest # type: ignore
from hammock.hammock import main
import os

@pytest.mark.full
def test_sequences():
    """Test sequence sketching with different methods."""
    # Create two similar sequences
    seq1 = "ATCGATCGATCG" * 10
    seq2 = "ATCGATCGATCG" * 8 + "GTTCATGCATAT" * 2
    
    # Test with different sketch types
    for sketch_type in ["minimizer", "hyperloglog", "minhash"]:
        print(f"\nTesting {sketch_type} sketch for sequences:")
        
        # Create sketches
        sketch1 = SequenceSketch(sketch_type=sketch_type, kmer_size=8)
        sketch2 = SequenceSketch(sketch_type=sketch_type, kmer_size=8)
        
        # Add sequences
        sketch1.add_string(seq1)
        sketch2.add_string(seq2)
        
        # Get similarity values
        result = sketch1.similarity_values(sketch2)
        jaccard = result['jaccard_similarity']
        print(f"Jaccard similarity: {jaccard:.3f}")
        if sketch_type != "minimizer":
            print(f"Cardinality 1: {sketch1.sketch.estimate_cardinality():.0f}")
            print(f"Cardinality 2: {sketch2.sketch.estimate_cardinality():.0f}")
        else:
            print(f"Gap Similarity: {result['gap_similarity']:.3f}")

@pytest.mark.full
def test_intervals():
    """Test interval sketching with different methods."""
    # Create more test BED lines with distinct intervals
    bed1 = []
    bed2 = []
    
    # Create smaller test sets first for debugging
    n_intervals = 20  # Start with smaller sets
    n_shared = 4     # 20% overlap
    
    # Create test intervals
    for i in range(n_intervals):
        # First set: even positions
        bed1.append(f"chr1\t{2*i*100}\t{(2*i+1)*100}")
        # Second set: with 20% overlap with bed1
        if i < n_shared:  # shared intervals
            bed2.append(f"chr1\t{2*i*100}\t{(2*i+1)*100}")  # Same as bed1
        else:
            bed2.append(f"chr1\t{(2*i+1)*100}\t{(2*i+2)*100}")  # Different positions
    
    # Test each sketch type separately
    for sketch_type in ["minhash", "hyperloglog"]:  # Test MinHash first
        print(f"\n{'='*50}")
        print(f"Testing {sketch_type} sketch for intervals:")
        print(f"{'='*50}")
        
        # Create sketches
        sketch_args = {
            'kmer_size': 0,
            'mode': "A",
            'sketch_type': sketch_type
        }
        if sketch_type == "hyperloglog":
            sketch_args['precision'] = 14
            print(f"HyperLogLog precision: {sketch_args['precision']}")
        elif sketch_type == "minhash":
            sketch_args['num_hashes'] = 256
            print(f"MinHash num_hashes: {sketch_args['num_hashes']}")
        
        sketch1 = IntervalSketch(**sketch_args)
        sketch2 = IntervalSketch(**sketch_args)
        
        # Add intervals with debugging output
        print("\nFirst 3 intervals from each set:")
        print("Set 1:")
        for i, line in enumerate(bed1):
            interval, points, size = sketch1.bedline(line, mode="A", sep="-")
            if interval and i < 3:
                print(f"  {interval}")
            if interval:
                sketch1.sketch.add_string(interval)
        
        print("\nSet 2:")
        for i, line in enumerate(bed2):
            interval, points, size = sketch2.bedline(line, mode="A", sep="-")
            if interval and i < 3:
                print(f"  {interval}")
            if interval:
                sketch2.sketch.add_string(interval)
        
        # Get cardinality estimates
        # if sketch_type == "hyperloglog":
        print("\nCardinality estimates:")
        print(f"Set 1: {sketch1.sketch.estimate_cardinality()}")
        print(f"Set 2: {sketch2.sketch.estimate_cardinality()}")
    
        # Get similarity values
        result = sketch1.similarity_values(sketch2)
        jaccard = result['jaccard_similarity']
        
        # Calculate expected Jaccard
        total_unique = 2 * n_intervals - n_shared  # Total unique intervals
        expected_jaccard = n_shared / total_unique
        
        print(f"\nJaccard Similarity:")
        print(f"Expected: {expected_jaccard:.3f}")
        print(f"Got:      {jaccard:.3f}")
        
        # Verify with appropriate tolerance
        tolerance = 0.1 if sketch_type == "hyperloglog" else 0.05
        assert abs(jaccard - expected_jaccard) < tolerance, \
            f"Expected Jaccard ~{expected_jaccard:.3f}, got {jaccard:.3f} (type: {sketch_type})"

@pytest.mark.quick
def test_mode_switching(capsys):
    """Test automatic mode switching based on parameters."""
    # Setup basic arguments
    base_args = [
        "hammock",
        "test_paths.txt",
        "test_primary.txt",
        "--mode", "A"
    ]
    
    # Test default mode A
    with patch('sys.argv', base_args):
        with patch('builtins.open', mock_open(read_data='test.bed\n')):
            main()
    captured = capsys.readouterr()
    assert "Using default mode A for interval comparison" in captured.out
    
    # Test switching to mode C with subA
    c_mode_args = base_args + ["--subA", "0.5"]
    with patch('sys.argv', c_mode_args):
        with patch('builtins.open', mock_open(read_data='test.bed\n')):
            main()
    captured = capsys.readouterr()
    assert "C-mode parameters detected" in captured.out
    assert "switching from mode A to mode C" in captured.out
    
    # Test switching to mode C with subB
    c_mode_args = base_args + ["--subB", "0.5"]
    with patch('sys.argv', c_mode_args):
        with patch('builtins.open', mock_open(read_data='test.bed\n')):
            main()
    captured = capsys.readouterr()
    assert "C-mode parameters detected" in captured.out
    assert "switching from mode A to mode C" in captured.out
    
    # Test switching to mode C with expA
    c_mode_args = base_args + ["--expA", "1"]
    with patch('sys.argv', c_mode_args):
        with patch('builtins.open', mock_open(read_data='test.bed\n')):
            main()
    captured = capsys.readouterr()
    assert "C-mode parameters detected" in captured.out
    assert "switching from mode A to mode C" in captured.out

@pytest.mark.quick
def test_mode_validation():
    """Test validation of mode-specific parameters."""
    # Test that mode A with expA switches to mode C
    args = [
        "hammock",
        "test_paths.txt",
        "test_primary.txt",
        "--mode", "A",
        "--expA", "1"
    ]
    with patch('sys.argv', args):
        with patch('builtins.open', mock_open(read_data='test.bed\n')):
            main()  # Should switch to mode C instead of raising error
    
    # Test invalid subsampling in mode A
    args = [
        "hammock",
        "test_paths.txt",
        "test_primary.txt",
        "--mode", "A",
        "--subA", "0.5"
    ]
    with patch('sys.argv', args):
        with patch('builtins.open', mock_open(read_data='test.bed\n')):
            main()  # Should switch to mode C instead of raising error
    
    # Test negative expA in mode C
    args = [
        "hammock",
        "test_paths.txt",
        "test_primary.txt",
        "--mode", "C",
        "--expA", "-1"
    ]
    with patch('sys.argv', args):
        with patch('builtins.open', mock_open(read_data='test.bed\n')):
            with pytest.raises(ValueError, match="--expA parameter must be non-negative for mode C"):
                main()
    
    # Test invalid subsample rates
    args = [
        "hammock",
        "test_paths.txt",
        "test_primary.txt",
        "--mode", "C",
        "--subA", "1.5"
    ]
    with patch('sys.argv', args):
        with patch('builtins.open', mock_open(read_data='test.bed\n')):
            with pytest.raises(ValueError, match="Subsample rates must be between 0 and 1"):
                main()

@pytest.mark.quick
def test_sequence_mode_switching(capsys):
    """Test automatic switching to mode D for sequence files."""
    # Test with fasta file
    args = [
        "hammock",
        "test_paths.txt",
        "test_primary.txt",
        "--mode", "A"
    ]
    
    # Test switching to mode D with .fasta file
    with patch('sys.argv', args):
        with patch('builtins.open', mock_open(read_data='test.fasta\n')):
            main()
    captured = capsys.readouterr()
    assert "Detected sequence file format, switching to mode D" in captured.out
    
    # Test switching to mode D with other sequence formats
    for ext in ['.fa', '.fna', '.ffn', '.faa', '.frn']:
        with patch('sys.argv', args):
            with patch('builtins.open', mock_open(read_data=f'test{ext}\n')):
                main()
        captured = capsys.readouterr()
        assert "Detected sequence file format, switching to mode D" in captured.out
    
    # Test no switch with .bed file
    with patch('sys.argv', args):
        with patch('builtins.open', mock_open(read_data='test.bed\n')):
            main()
    captured = capsys.readouterr()
    assert "Detected sequence file format, switching to mode D" not in captured.out

if __name__ == "__main__":
    print("Testing sequence sketching...")
    test_sequences()
    
    print("\nTesting interval sketching...")
    test_intervals() 