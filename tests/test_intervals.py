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
def test_subsampling():
    """Test interval subsampling."""
    total_intervals = 1000  # Increase number of intervals for more reliable sampling
    subsample_ratio = 0.5  # 50% sampling
    
    # Create test bed file with known number of intervals
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        for i in range(total_intervals):
            f.write(f"chr1\t{i*100}\t{(i+1)*100}\n")
        f.flush()
        
        try:
            # Test mode A (no subsampling allowed)
            sketch_a = IntervalSketch.from_file(
                filename=f.name,
                mode="A",
                subsample=(1.0, 1.0)  # No subsampling in mode A
            )
            assert sketch_a is not None
            assert sketch_a.num_intervals == total_intervals  # Should have all intervals
            
            # Test mode B (only point subsampling)
            sketch_b = IntervalSketch.from_file(
                filename=f.name,
                mode="B",
                subsample=(1.0, subsample_ratio)  # Only subsample points
            )
            assert sketch_b is not None
            assert sketch_b.sketch.estimate_cardinality() > 0  # Should have some points
            
            # Test mode C (both interval and point subsampling)
            sketch_c = IntervalSketch.from_file(
                filename=f.name,
                mode="C",
                subsample=(subsample_ratio, subsample_ratio)  # Sample both intervals and points
            )
            assert sketch_c is not None
            expected = total_intervals * subsample_ratio
            tolerance = expected * 0.2  # Allow 20% variation due to random sampling
            assert abs(sketch_c.num_intervals - expected) <= tolerance, \
                f"Expected ~{expected} intervals (±{tolerance}), got {sketch_c.num_intervals}"
            assert sketch_c.sketch.estimate_cardinality() > 0
            
        finally:
            os.unlink(f.name)


@pytest.mark.full
def test_unsupported_file():
    """Test handling of unsupported file formats."""
    with tempfile.NamedTemporaryFile(suffix='.cram', delete=False) as f:
        sketch = IntervalSketch.from_file(filename=f.name, mode="A")
        assert sketch is None

@pytest.mark.full
def test_error_handling():
    """Test error handling for malformed files."""
    # Test non-existent file
    sketch = IntervalSketch.from_file(filename="nonexistent.bed", mode="A")
    assert sketch is None

    # Test malformed Bed file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write("ch1\t100\t200\t300\n")
        f.flush()
        sketch = IntervalSketch.from_file(filename=f.name, mode="A")
        assert sketch is None
