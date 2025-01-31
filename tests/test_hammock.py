#!/usr/bin/env python
from __future__ import annotations
from hammock.lib.sequences import SequenceSketch
from hammock.lib.intervals import IntervalSketch
from unittest.mock import patch, mock_open
import pytest # type: ignore
from hammock.hammock import main

@pytest.mark.full
def test_sequences():
    """Test sequence sketching with different methods."""
    # Create two similar sequences
    seq1 = "ATCGATCGATCG" * 10
    seq2 = "ATCGATCGATCG" * 8 + "GCATGCATGCAT" * 2
    
    # Test with different sketch types
    for sketch_type in ["minimizer", "hyperloglog", "minhash"]:
        print(f"\nTesting {sketch_type} sketch for sequences:")
        
        # Create sketches
        sketch1 = SequenceSketch(sketch_type=sketch_type, kmer_size=8)
        sketch2 = SequenceSketch(sketch_type=sketch_type, kmer_size=8)
        
        # Add sequences
        sketch1.add_sequence(seq1)
        sketch2.add_sequence(seq2)
        
        # Compare
        jaccard = sketch1.sketch.estimate_jaccard(sketch2.sketch)
        print(f"Jaccard similarity: {jaccard:.3f}")
        print(f"Cardinality 1: {sketch1.sketch.estimate_cardinality():.0f}")
        print(f"Cardinality 2: {sketch2.sketch.estimate_cardinality():.0f}")

@pytest.mark.full
def test_intervals():
    """Test interval sketching with different methods."""
    # Create some test BED lines with more distinct intervals
    bed1 = [
        "chr1\t100\t200",
        "chr1\t300\t400",
        "chr2\t100\t200",
        "chr3\t100\t200",
        "chr4\t100\t200"
    ]
    bed2 = [
        "chr1\t100\t200",
        "chr1\t300\t400",
        "chr5\t100\t200",
        "chr6\t100\t200",
        "chr7\t100\t200"
    ]
    
    # Test with different sketch types
    for sketch_type in ["hyperloglog", "minhash"]:
        print(f"\nTesting {sketch_type} sketch for intervals:")
        
        # Create sketches with higher precision/num_hashes for better accuracy
        sketch_args = {
            'precision': 14 if sketch_type == "hyperloglog" else 256,  # Increase precision
            'mode': "A",
            'sketch_type': sketch_type
        }
        
        sketch1 = IntervalSketch(**sketch_args)
        sketch2 = IntervalSketch(**sketch_args)
        
        # Add intervals
        for line in bed1:
            interval, points, size = sketch1.bedline(line, mode="A", sep="-")
            if interval:
                sketch1.sketch.add_string(interval)
                
        for line in bed2:
            interval, points, size = sketch2.bedline(line, mode="A", sep="-")
            if interval:
                sketch2.sketch.add_string(interval)
        
        # Compare
        jaccard = sketch1.sketch.estimate_jaccard(sketch2.sketch)
        
        # Verify expected Jaccard similarity (2 shared intervals out of 8 total)
        expected_jaccard = 2/8  # 2 shared intervals out of 8 unique intervals
        assert abs(jaccard - expected_jaccard) < 0.1, f"Expected Jaccard ~{expected_jaccard}, got {jaccard:.3f}"
        
        # Verify cardinalities
        assert abs(sketch1.sketch.estimate_cardinality() - 5) < 0.5
        assert abs(sketch2.sketch.estimate_cardinality() - 5) < 0.5

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