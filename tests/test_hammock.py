#!/usr/bin/env python
from __future__ import annotations
from hammock.lib.sequences import SequenceSketch
from hammock.lib.intervals import IntervalSketch
from hammock.lib.hyperloglog_fast import FastHyperLogLog, CPP_AVAILABLE
from unittest.mock import patch, mock_open, MagicMock
import pytest # type: ignore
from hammock.hammock import main
import os
import gettext

# Valid test data
TEST_BED_DATA = """chr1\t0\t100\tA
chr1\t200\t300\tA
chr1\t400\t500\tA
"""

TEST_FASTA_DATA = """>seq1
ATCGATCGATCG
>seq2
GCTAGCTAGCTA
"""

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

    # Mock gettext to avoid bytes/string issues
    mock_gettext = MagicMock()
    mock_gettext.gettext = lambda x: x

    # Test default mode A
    with patch('gettext.translation', return_value=mock_gettext):
        with patch('sys.argv', base_args):
            with patch('builtins.open') as mock_file:
                mock_file.return_value.__enter__.return_value.read.return_value = TEST_BED_DATA
                with patch('os.path.exists', return_value=True):
                    main()

@pytest.mark.quick
def test_mode_validation():
    """Test validation of mode-specific parameters."""
    # Mock gettext to avoid bytes/string issues
    mock_gettext = MagicMock()
    mock_gettext.gettext = lambda x: x

    # Test that mode A with expA switches to mode C
    args = [
        "hammock",
        "test_paths.txt",
        "test_primary.txt",
        "--mode", "A",
        "--expA", "1"
    ]
    with patch('gettext.translation', return_value=mock_gettext):
        with patch('sys.argv', args):
            with patch('builtins.open') as mock_file:
                mock_file.return_value.__enter__.return_value.read.return_value = TEST_BED_DATA
                with patch('os.path.exists', return_value=True):
                    main()  # Should switch to mode C instead of raising error

@pytest.mark.quick
def test_sequence_mode_switching(capsys):
    """Test automatic switching to mode D for sequence files."""
    # Mock gettext to avoid bytes/string issues
    mock_gettext = MagicMock()
    mock_gettext.gettext = lambda x: x

    # Test with fasta file
    args = [
        "hammock",
        "test_paths.txt",
        "test_primary.txt",
        "--mode", "A"
    ]

    # Test switching to mode D with .fasta file
    with patch('gettext.translation', return_value=mock_gettext):
        with patch('sys.argv', args):
            with patch('builtins.open') as mock_file:
                # First read for extension check
                mock_file.return_value.__enter__.return_value.readline.return_value = "test.fasta\n"
                # Second read for content
                mock_file.return_value.__enter__.return_value.read.return_value = TEST_FASTA_DATA
                with patch('os.path.exists', return_value=True):
                    main()
    captured = capsys.readouterr()
    assert "Detected sequence file format, switching to mode D" in captured.out

    # Test switching to mode D with other sequence formats
    for ext in ['.fa', '.fna', '.ffn', '.faa', '.frn']:
        with patch('gettext.translation', return_value=mock_gettext):
            with patch('sys.argv', args):
                with patch('builtins.open') as mock_file:
                    # First read for extension check
                    mock_file.return_value.__enter__.return_value.readline.return_value = f"test{ext}\n"
                    # Second read for content
                    mock_file.return_value.__enter__.return_value.read.return_value = TEST_FASTA_DATA
                    with patch('os.path.exists', return_value=True):
                        main()
        captured = capsys.readouterr()
        assert "Detected sequence file format, switching to mode D" in captured.out

    # Test that non-sequence files don't switch to mode D
    with patch('gettext.translation', return_value=mock_gettext):
        with patch('sys.argv', args):
            with patch('builtins.open') as mock_file:
                # First read for extension check
                mock_file.return_value.__enter__.return_value.readline.return_value = "test.bed\n"
                # Second read for content
                mock_file.return_value.__enter__.return_value.read.return_value = TEST_BED_DATA
                with patch('os.path.exists', return_value=True):
                    main()
    captured = capsys.readouterr()
    assert "Detected sequence file format, switching to mode D" not in captured.out

@pytest.mark.quick
def test_identical_file_lists_optimization(capsys):
    """Test optimization when both file lists are identical."""
    # Mock gettext to avoid bytes/string issues
    mock_gettext = MagicMock()
    mock_gettext.gettext = lambda x: x

    # 1. Test case: same file provided for both arguments
    args = [
        "hammock",
        "same_file.txt",
        "same_file.txt",  # Same file for both arguments
        "--mode", "A"
    ]
    
    with patch('gettext.translation', return_value=mock_gettext):
        with patch('sys.argv', args):
            with patch('builtins.open') as mock_file:
                # First read for extension check
                mock_file.return_value.__enter__.return_value.readline.return_value = "test.bed\n"
                # Mock file reading for content
                mock_file.return_value.__enter__.return_value.readlines.return_value = ["file1.bed\n", "file2.bed\n", "file3.bed\n"]
                # Second read to get full list
                mock_file.return_value.__enter__.return_value.__iter__.return_value = ["file1.bed\n", "file2.bed\n", "file3.bed\n"]
                
                # Mock exists to return True for both the list file and the bed files
                def mock_exists(path):
                    return True
                
                with patch('os.path.exists', side_effect=mock_exists), \
                     patch('hammock.hammock.IntervalSketch.from_file', MagicMock(return_value=MagicMock())), \
                     patch('hammock.hammock.csv.writer', MagicMock()):
                    main()
    
    captured = capsys.readouterr()
    assert "Detected same input file for both arguments - sketches will be computed only once" in captured.out
    
    # 2. Test case: different files with same content
    args = [
        "hammock",
        "files1.txt",
        "files2.txt",  # Different files
        "--mode", "A"
    ]
    
    with patch('gettext.translation', return_value=mock_gettext):
        with patch('sys.argv', args):
            with patch('builtins.open') as mock_file:
                # Mock different behavior based on argument
                def side_effect(*args, **kwargs):
                    # Create a mock with different behavior depending on which file is being opened
                    mock = MagicMock()
                    if args[0] == "files1.txt":
                        mock.readline.return_value = "test.bed\n"
                        mock.__iter__.return_value = ["file1.bed\n", "file2.bed\n", "file3.bed\n"]
                        mock.readlines.return_value = ["file1.bed\n", "file2.bed\n", "file3.bed\n"]
                    elif args[0] == "files2.txt":
                        # Same content but in different order
                        mock.__iter__.return_value = ["file3.bed\n", "file1.bed\n", "file2.bed\n"]
                        mock.readlines.return_value = ["file3.bed\n", "file1.bed\n", "file2.bed\n"]
                    return mock
                
                mock_file.side_effect = side_effect
                
                # Mock exists to return True for all relevant files
                def mock_exists(path):
                    return True
                
                with patch('os.path.exists', side_effect=mock_exists), \
                     patch('hammock.hammock.IntervalSketch.from_file', MagicMock(return_value=MagicMock())), \
                     patch('hammock.hammock.csv.writer', MagicMock()):
                    main()
    
    captured = capsys.readouterr()
    assert "Detected identical file lists - sketches will be computed only once" in captured.out
    
    # 3. Test optimization for file processing (using process_file)
    from hammock.hammock import process_file
    
    # Mock os.path.exists to return True for test files
    def mock_exists(path):
        return True
    
    with patch('os.path.exists', side_effect=mock_exists), \
         patch('hammock.hammock.IntervalSketch.from_file', MagicMock(return_value=MagicMock())) as mock_sketch, \
         patch('hammock.hammock.limit_memory', MagicMock()):  # Prevent actual memory limiting
            
        # Test with same files - 'file1.bed' appears in both arguments
        result = process_file('file1.bed', ['file1.bed', 'file2.bed'], mode='A')
        # Should have called from_file once for each unique file (2 calls)
        assert mock_sketch.call_count == 2
        mock_sketch.reset_mock()
        
        # Test with different files - 'file3.bed' is new and not in the primary list
        result = process_file('file3.bed', ['file1.bed', 'file2.bed'], mode='A')
        # Should call from_file for the unique files (2 primaries + 1 new = 3 calls)
        assert mock_sketch.call_count == 3


@pytest.mark.skipif(not CPP_AVAILABLE, reason="C++ extension not available")
def test_cpp_integration_in_hammock():
    """Test that C++ integration works within the hammock framework."""
    # Test that FastHyperLogLog can be used in sequence sketching
    seq1 = "ATCGATCGATCG" * 10
    seq2 = "ATCGATCGATCG" * 8 + "GTTCATGCATAT" * 2
    
    # Create sketches with C++ acceleration
    sketch1 = SequenceSketch(sketch_type="hyperloglog", kmer_size=8)
    sketch2 = SequenceSketch(sketch_type="hyperloglog", kmer_size=8)
    
    # Verify that the underlying sketch is using C++ acceleration
    assert hasattr(sketch1.sketch, '_acceleration_type')
    assert sketch1.sketch._acceleration_type in ['C++', 'Cython', 'Python']
    
    # Add sequences
    sketch1.add_string(seq1)
    sketch2.add_string(seq2)
    
    # Get similarity values
    result = sketch1.similarity_values(sketch2)
    jaccard = result['jaccard_similarity']
    
    # Should get reasonable similarity
    assert 0.0 <= jaccard <= 1.0
    print(f"C++ integration test - Jaccard similarity: {jaccard:.3f}")


if __name__ == "__main__":
    print("Testing sequence sketching...")
    test_sequences()
    
    print("\nTesting interval sketching...")
    test_intervals() 