#!/usr/bin/env python
from __future__ import annotations
import tempfile
import pytest # type: ignore
import pyBigWig # type: ignore
from hammock.lib.intervals import IntervalSketch
import os
import pysam # type: ignore

@pytest.mark.quick
def test_interval_sketch_init():
    """Test IntervalSketch initialization with different modes."""
    # Test mode A
    sketch = IntervalSketch(mode="A")
    assert sketch.mode == "A"
    assert sketch.expA == 0
    assert sketch.subsample == (1.0, 1.0)
    assert sketch.num_intervals == 0  # Start with 0 intervals

    # Test mode B
    sketch = IntervalSketch(mode="B")
    assert sketch.mode == "B"
    assert sketch.num_intervals == 0  # Start with 0 intervals

    # Test mode C
    sketch = IntervalSketch(mode="C")
    assert sketch.mode == "C"
    assert sketch.num_intervals == 0  # Start with 0 intervals

@pytest.mark.quick
def test_invalid_modes():
    """Test that invalid modes raise appropriate errors."""
    with pytest.raises(ValueError):
        IntervalSketch(mode="X")

@pytest.mark.quick
def test_expA_validation():
    """Test expA parameter validation."""
    # Test valid expA in mode C
    sketch = IntervalSketch(mode="C", expA=1.0)
    assert sketch.expA == 1.0

    # Test expA > 0 in mode A should raise error
    with pytest.raises(ValueError, match="Multiplicity .expA. can only be used with mode C"):
        IntervalSketch(mode="A", expA=1.0)

    # Test negative expA should raise error
    with pytest.raises(ValueError):
        IntervalSketch(mode="C", expA=-1.0)

@pytest.mark.quick
def test_subsample_validation():
    """Test subsample parameter validation."""
    # Test valid subsample rates
    sketch = IntervalSketch(mode="C", subsample=(0.5, 0.5))
    assert sketch.subsample == (0.5, 0.5)

    # Test invalid subsample rates
    with pytest.raises(ValueError):
        IntervalSketch(mode="C", subsample=(1.5, 0.5))
    with pytest.raises(ValueError):
        IntervalSketch(mode="C", subsample=(0.5, -0.1))

@pytest.mark.full
def test_expA_from_file():
    """Test that expA parameter is properly passed through from_file."""
    # Create a temporary bed file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed') as f:
        f.write("chr1\t100\t200\n")
        f.flush()
        
        # Test with expA > 0
        sketch = IntervalSketch.from_file(
            filename=f.name,
            mode="C",
            expA=1.0
        )
        assert sketch.expA == 1.0
        
        # Test with default expA
        sketch = IntervalSketch.from_file(
            filename=f.name,
            mode="C"
        )
        assert sketch.expA == 0

@pytest.mark.full
def test_bed_file_processing():
    """Test processing of BED files."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed') as f:
        # Write test intervals
        f.write("chr1\t100\t200\n")
        f.write("chr1\t150\t250\n")
        f.write("chr2\t300\t400\n")
        f.flush()

        # Test mode A (intervals)
        sketch_a = IntervalSketch.from_file(filename=f.name, mode="A")
        assert sketch_a is not None
        assert sketch_a.num_intervals == 3  # Should count all intervals
        assert sketch_a.sketch.estimate_cardinality() > 0  # Should have interval hashes

        # Test mode B (points)
        sketch_b = IntervalSketch.from_file(filename=f.name, mode="B")
        assert sketch_b is not None
        assert sketch_b.num_intervals == 3  # Should count number of intervals in file
        assert sketch_b.sketch.estimate_cardinality() > 0  # Should have points

        # Test mode C (both)
        sketch_c = IntervalSketch.from_file(filename=f.name, mode="C")
        assert sketch_c is not None
        assert sketch_c.num_intervals == 3  # Should count intervals
        assert sketch_c.sketch.estimate_cardinality() > 0  # Should have both intervals and points

@pytest.mark.full
def test_invalid_file():
    """Test handling of invalid files."""
    # Test with non-existent file
    sketch = IntervalSketch.from_file(filename="nonexistent.bed", mode="A")
    assert sketch is None

    # Test with invalid format
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed') as f:
        f.write("invalid format\n")
        f.flush()
        sketch = IntervalSketch.from_file(filename=f.name, mode="A")
        assert sketch is None

@pytest.mark.full
def test_bigwig_file_processing():
    """Test processing of BigWig files."""
    with tempfile.NamedTemporaryFile(suffix='.bw', delete=False) as f:
        # Create a test BigWig file
        bw = pyBigWig.open(f.name, 'w')
        # Add a header with chromosome sizes
        chroms = [("chr1", 1000), ("chr2", 1000)]
        bw.addHeader(list(chroms))
        
        # Add some test data
        # chr1: two intervals with non-zero values
        chroms_values = [
            ("chr1", 100, 200, 1.0),
            ("chr1", 300, 400, 1.0),
            ("chr2", 500, 600, 1.0)
        ]
        for chrom, start, end, value in chroms_values:
            bw.addEntries([chrom], [start], ends=[end], values=[value])
        bw.close()

        # Test mode A (intervals)
        sketch_a = IntervalSketch.from_file(filename=f.name, mode="A")
        assert sketch_a is not None
        assert sketch_a.num_intervals == 3  # Should have 3 intervals
        assert sketch_a.total_interval_size == 300  # Total size: (200-100) + (400-300) + (600-500)

        # Test mode B (points)
        sketch_b = IntervalSketch.from_file(filename=f.name, mode="B")
        assert sketch_b is not None
        assert sketch_b.num_intervals == 3  # Should count number of intervals in file
        assert sketch_b.sketch.estimate_cardinality() > 0  # Should have points in the sketch

        # Test mode C (both)
        sketch_c = IntervalSketch.from_file(filename=f.name, mode="C")
        assert sketch_c is not None
        assert sketch_c.num_intervals == 3  # Should count intervals in mode C
@pytest.mark.full
def test_bigwig_invalid_file():
    """Test handling of invalid BigWig files."""
    # Test with non-existent file
    sketch = IntervalSketch.from_file(filename="nonexistent.bw", mode="A")
    assert sketch is None

    # Test with invalid format
    with tempfile.NamedTemporaryFile(suffix='.bw') as f:
        f.write(b"invalid bigwig format")
        f.flush()
        sketch = IntervalSketch.from_file(filename=f.name, mode="A")
        assert sketch is None

@pytest.mark.full
def test_bigwig_with_expA():
    """Test bigWig comparison with expA."""
    # Check if test files exist
    test1_path = "tests/data/test1.bigwig"
    test2_path = "tests/data/test2.bigwig"
    
    assert os.path.exists(test1_path), f"Test file not found: {test1_path}"
    assert os.path.exists(test2_path), f"Test file not found: {test2_path}"
    
    sketch1 = IntervalSketch.from_file(
        filename=test1_path,
        mode="C",
        expA=2.0,
        sketch_type="hyperloglog"
    )
    assert sketch1 is not None, "Failed to create first sketch"
    
    sketch2 = IntervalSketch.from_file(
        filename=test2_path,
        mode="C",
        expA=2.0,
        sketch_type="hyperloglog"
    )
    assert sketch2 is not None, "Failed to create second sketch"
    
    # Use similarity_values
    result = sketch1.similarity_values(sketch2)
    assert 'jaccard_similarity' in result
    assert 'containment' in result
    assert result['containment'] is not None

@pytest.mark.full
def test_bigwig_with_subsampling():
    """Test BigWig processing with subsampling."""
    with tempfile.NamedTemporaryFile(suffix='.bw', delete=False) as f:
        # Create a test BigWig file
        bw = pyBigWig.open(f.name, 'w')
        # Add a header with chromosome sizes
        chroms = [("chr1", 1000)]
        bw.addHeader(list(chroms))
        
        # Add multiple intervals
        for i in range(5):
            bw.addEntries(["chr1"], [i*200], ends=[(i+1)*200], values=[1.0])
        bw.close()

        # Test with subsampling in mode C
        sketch = IntervalSketch.from_file(
            filename=f.name,
            mode="C",
            subsample=(0.5, 0.5)  # Subsample both intervals and points
        )
        assert sketch is not None
        assert sketch.num_intervals <= 5  # Should have fewer intervals due to subsampling
        assert sketch.sketch.estimate_cardinality() > 0  # Should still have some points

@pytest.mark.full
def test_subsampling():
    """Test interval subsampling."""
    total_intervals = 10
    subsample_ratio = 0.5  # 50% sampling
    
    # Create test bed file with known number of intervals
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        for i in range(total_intervals):
            f.write(f"chr1\t{i*100}\t{(i+1)*100}\n")
        f.flush()
        
        try:
            # Test mode A with subsampling
            sketch_a = IntervalSketch.from_file(
                filename=f.name,
                mode="A",
                subsample=(subsample_ratio, 0.0)  # Sample 50% of intervals, no points
            )
            assert sketch_a is not None
            expected = total_intervals * subsample_ratio
            tolerance = 2  # Allow some variation due to random sampling
            assert abs(sketch_a.num_intervals - expected) <= tolerance, \
                f"Expected ~{expected} intervals (±{tolerance}), got {sketch_a.num_intervals}"
            
            # Test mode B with subsampling
            sketch_b = IntervalSketch.from_file(
                filename=f.name,
                mode="B",
                subsample=(0.0, subsample_ratio)  # No intervals, sample 50% of points
            )
            assert sketch_b is not None
            assert sketch_b.num_intervals == total_intervals  # Should count all intervals in file
            assert sketch_b.sketch.estimate_cardinality() > 0  # Should have some points
            
        finally:
            os.unlink(f.name)

@pytest.mark.full
def test_bam_file_processing():
    """Test processing of BAM files."""
    # Create a temporary BAM file
    with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as f:
        # Create a test BAM file
        header = {'HD': {'VN': '1.0'},
                 'SQ': [{'LN': 1000, 'SN': 'chr1'},
                       {'LN': 1000, 'SN': 'chr2'}]}
        
        with pysam.AlignmentFile(f.name, 'wb', header=header) as bam:
            # Create test alignments
            a1 = pysam.AlignedSegment()
            a1.query_name = "read1"
            a1.reference_id = 0  # chr1
            a1.reference_start = 100
            a1.reference_end = 200
            a1.query_sequence = "A" * 100
            a1.flag = 0
            a1.mapping_quality = 20
            a1.cigar = [(0, 100)]  # 100M (match)
            
            a2 = pysam.AlignedSegment()
            a2.query_name = "read2"
            a2.reference_id = 0  # chr1
            a2.reference_start = 150
            a2.reference_end = 250
            a2.query_sequence = "A" * 100
            a2.flag = 0
            a2.mapping_quality = 20
            a2.cigar = [(0, 100)]  # 100M (match)
            
            a3 = pysam.AlignedSegment()
            a3.query_name = "read3"
            a3.reference_id = 1  # chr2
            a3.reference_start = 300
            a3.reference_end = 400
            a3.query_sequence = "A" * 100
            a3.flag = 0
            a3.mapping_quality = 20
            a3.cigar = [(0, 100)]  # 100M (match)
            
            # Write alignments
            bam.write(a1)
            bam.write(a2)
            bam.write(a3)
        
        # Create index for the BAM file
        pysam.index(f.name)
        
        try:
            # Test mode A (intervals)
            sketch_a = IntervalSketch.from_file(filename=f.name, mode="A")
            assert sketch_a is not None
            assert sketch_a.num_intervals == 3  # Should count all intervals
            assert sketch_a.sketch.estimate_cardinality() > 0  # Should have interval hashes
            
            # Test mode B (points)
            sketch_b = IntervalSketch.from_file(filename=f.name, mode="B")
            assert sketch_b is not None
            assert sketch_b.num_intervals == 3  # Should count number of intervals
            assert sketch_b.sketch.estimate_cardinality() > 0  # Should have points
            
            # Test mode C (both)
            sketch_c = IntervalSketch.from_file(filename=f.name, mode="C")
            assert sketch_c is not None
            assert sketch_c.num_intervals == 3  # Should count intervals
            assert sketch_c.sketch.estimate_cardinality() > 0  # Should have both intervals and points
            
        finally:
            # Clean up
            os.unlink(f.name)
            if os.path.exists(f.name + '.bai'):
                os.unlink(f.name + '.bai')

@pytest.mark.full
def test_bam_invalid_file():
    """Test handling of invalid BAM files."""
    # Test with non-existent file
    sketch = IntervalSketch.from_file(filename="nonexistent.bam", mode="A")
    assert sketch is None
    
    # Test with invalid format
    with tempfile.NamedTemporaryFile(suffix='.bam') as f:
        f.write(b"invalid bam format")
        f.flush()
        sketch = IntervalSketch.from_file(filename=f.name, mode="A")
        assert sketch is None

@pytest.mark.full
def test_bam_with_subsampling():
    """Test BAM processing with subsampling."""
    with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as f:
        # Create header
        header = {'HD': {'VN': '1.0'},
                 'SQ': [{'LN': 1000, 'SN': 'chr1'}]}
        
        with pysam.AlignmentFile(f.name, 'wb', header=header) as bam:
            # Create 10 test alignments
            for i in range(10):
                a = pysam.AlignedSegment()
                a.query_name = f"read{i}"
                a.reference_id = 0
                a.reference_start = i * 100
                a.reference_end = (i + 1) * 100
                a.query_sequence = "A" * 100
                a.flag = 0
                a.mapping_quality = 20
                a.cigar = [(0, 100)]
                bam.write(a)
        
        # Create index
        pysam.index(f.name)
        
        try:
            # Test with subsampling in mode C
            sketch = IntervalSketch.from_file(
                filename=f.name,
                mode="C",
                subsample=(0.5, 0.5)  # Subsample both intervals and points
            )
            assert sketch is not None
            assert sketch.num_intervals <= 10  # Should have fewer intervals due to subsampling
            assert sketch.sketch.estimate_cardinality() > 0  # Should still have some points
            
        finally:
            # Clean up
            os.unlink(f.name)
            if os.path.exists(f.name + '.bai'):
                os.unlink(f.name + '.bai')
