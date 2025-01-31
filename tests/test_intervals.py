#!/usr/bin/env python
from __future__ import annotations
import tempfile
import pytest # type: ignore
import pyBigWig # type: ignore
from hammock.lib.intervals import IntervalSketch

@pytest.mark.quick
def test_interval_sketch_init():
    """Test IntervalSketch initialization with different modes."""
    # Test mode A
    sketch = IntervalSketch(mode="A")
    assert sketch.mode == "A"
    assert sketch.expA == 0
    assert sketch.subsample == (1.0, 1.0)

    # Test mode B
    sketch = IntervalSketch(mode="B")
    assert sketch.mode == "B"

    # Test mode C
    sketch = IntervalSketch(mode="C")
    assert sketch.mode == "C"

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

        # Test mode B (points)
        sketch_b = IntervalSketch.from_file(filename=f.name, mode="B")
        assert sketch_b is not None

        # Test mode C (both)
        sketch_c = IntervalSketch.from_file(filename=f.name, mode="C")
        assert sketch_c is not None

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
        assert sketch_b.num_intervals == 300  # Should have one point per position

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
    """Test BigWig processing with expA parameter."""
    with tempfile.NamedTemporaryFile(suffix='.bw', delete=False) as f:
        # Create a test BigWig file
        bw = pyBigWig.open(f.name, 'w')
        # Add a header with chromosome sizes
        chroms = [("chr1", 1000)]
        bw.addHeader(list(chroms))
        
        # Add one test interval
        bw.addEntries(["chr1"], [100], ends=[200], values=[1.0])
        bw.close()

        # Test with expA > 0 in mode C
        sketch = IntervalSketch.from_file(
            filename=f.name,
            mode="C",
            expA=1.0  # This should create multiple copies of each interval
        )
        assert sketch is not None
        assert sketch.expA == 1.0

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
        # The actual number of intervals/points will be probabilistic due to subsampling
        assert sketch.num_intervals <= 5  # Should have fewer than total intervals due to subsampling
