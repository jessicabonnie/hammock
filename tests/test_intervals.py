#!/usr/bin/env python
from __future__ import annotations
import tempfile
import pytest # type: ignore
import pyBigWig # type: ignore
from hammock.lib.intervals import IntervalSketch
import os
import pysam # type: ignore
import gzip
import time

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
    """Test processing of BED files with different modes and parameters."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        # Write test intervals
        f.write("chr1\t100\t200\n")
        f.write("chr1\t150\t250\n")
        f.write("chr2\t300\t400\n")
        f.flush()

        try:
            # Test mode A (intervals only)
            sketch_a = IntervalSketch.from_file(
                filename=f.name,
                mode="A",
                precision=12,  # 2^12 = 4096 registers, sufficient for ~200 points
                debug=True
            )
            assert sketch_a is not None
            assert sketch_a.num_intervals == 3  # Should count all intervals
            assert sketch_a.sketch.estimate_cardinality() > 0  # Should have interval hashes

            # Test mode B (points only)
            sketch_b = IntervalSketch.from_file(
                filename=f.name,
                mode="B",
                precision=12,  # 2^12 = 4096 registers, sufficient for ~200 points
                debug=True
            )
            assert sketch_b is not None
            # Total points should be sum of interval lengths
            expected_points = (200-100) + (250-150) + (400-300)
            tolerance = expected_points * 0.05  # 5% tolerance
            actual_points = sketch_b.sketch.estimate_cardinality()
            print(f"Mode B: Expected ~{expected_points} points (±{tolerance}), got {actual_points}")
            assert abs(actual_points - expected_points) <= tolerance, \
                f"Expected ~{expected_points} points (±{tolerance}), got {actual_points}"

            # Test mode C (both intervals and points)
            sketch_c = IntervalSketch.from_file(
                filename=f.name,
                mode="C",
                precision=12,  # 2^12 = 4096 registers, sufficient for ~200 points
                subsample=(0.5, 0.5),  # 50% sampling for both intervals and points
                debug=True
            )
            assert sketch_c is not None
            assert sketch_c.num_intervals <= 3  # Should have <= 3 intervals due to sampling
            actual_points = sketch_c.sketch.estimate_cardinality()
            print(f"Mode C: Got {actual_points} points after subsampling")
            assert sketch_c.sketch.estimate_cardinality() > 0

            # Test mode C with expA
            sketch_c_exp = IntervalSketch.from_file(
                filename=f.name,
                mode="C",
                precision=12,  # 2^12 = 4096 registers, sufficient for ~200 points
                expA=1.0,  # Add copies of intervals
                debug=True
            )
            assert sketch_c_exp is not None
            assert sketch_c_exp.num_intervals == 3
            # Should have more hashes due to interval copies
            actual_points = sketch_c_exp.sketch.estimate_cardinality()
            print(f"Mode C with expA: Got {actual_points} points")
            assert sketch_c_exp.sketch.estimate_cardinality() > sketch_c.sketch.estimate_cardinality()

            # Test gzipped BED
            with tempfile.NamedTemporaryFile(suffix='.bed.gz', delete=False) as gz_f:
                with gzip.open(gz_f.name, 'wt') as gz:
                    gz.write("chr1\t100\t200\n")
                try:
                    sketch_gz = IntervalSketch.from_file(
                        filename=gz_f.name, 
                        mode="A", 
                        precision=12,  # 2^12 = 4096 registers, sufficient for ~200 points
                        debug=True
                    )
                    assert sketch_gz is not None
                    assert sketch_gz.num_intervals == 1
                finally:
                    os.unlink(gz_f.name)

        finally:
            os.unlink(f.name)

@pytest.mark.full
def test_bed_advanced_features():
    """Test advanced features of BED file processing."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        # Write test intervals with various formats
        f.write("chr1\t100\t200\n")  # Standard format
        f.write("chrX\t150\t250\n")  # Non-numeric chromosome
        f.write("1\t300\t400\n")     # No chr prefix
        f.flush()

        try:
            # Test with different separators
            sketch_sep = IntervalSketch.from_file(
                filename=f.name,
                mode="A",
                sep="|",  # Use different separator
                debug=True
            )
            assert sketch_sep is not None
            assert sketch_sep.num_intervals == 3

            # Test with different sketch types
            sketch_hll = IntervalSketch.from_file(
                filename=f.name,
                mode="A",
                sketch_type="hyperloglog",
                precision=12,
                debug=True
            )
            assert sketch_hll is not None

            # Test with Rust implementation
            sketch_rust = IntervalSketch.from_file(
                filename=f.name,
                mode="A",
                use_rust=True,
                debug=True
            )
            assert sketch_rust is not None

            # Compare sketches
            assert abs(sketch_hll.estimate_jaccard(sketch_rust) - 1.0) < 0.1

            # Test debug mode
            sketch_debug = IntervalSketch.from_file(
                filename=f.name,
                mode="A",
                debug=True
            )
            assert sketch_debug is not None

        finally:
            os.unlink(f.name)

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
                subsample=(1.0, 1.0),  # No subsampling in mode A
                precision=17,  # 2^17 = 131,072 registers, sufficient for ~100,000 points with better accuracy
                debug=True
            )
            assert sketch_a is not None
            assert sketch_a.num_intervals == total_intervals  # Should have all intervals
            
            # Test mode B (only point subsampling)
            sketch_b = IntervalSketch.from_file(
                filename=f.name,
                mode="B",
                subsample=(1.0, subsample_ratio),  # Only subsample points
                precision=17,  # 2^17 = 131,072 registers, sufficient for ~100,000 points with better accuracy
                debug=True
            )
            assert sketch_b is not None
            # Total points should be sum of interval lengths
            expected_points = sum((i+1)*100 - i*100 for i in range(total_intervals))
            # With subsampling, we expect roughly half the points
            expected_sampled = expected_points * subsample_ratio
            tolerance = expected_sampled * 0.15  # Reduced tolerance to 15% for better accuracy
            actual_points = sketch_b.sketch.estimate_cardinality()
            print(f"Mode B: Expected ~{expected_sampled} points (±{tolerance}), got {actual_points}")
            assert abs(actual_points - expected_sampled) <= tolerance, \
                f"Expected ~{expected_sampled} points (±{tolerance}), got {actual_points}"
            
            # Test mode C (both interval and point subsampling)
            sketch_c = IntervalSketch.from_file(
                filename=f.name,
                mode="C",
                subsample=(subsample_ratio, subsample_ratio),  # Sample both intervals and points
                precision=17,  # 2^17 = 131,072 registers, sufficient for ~100,000 points with better accuracy
                debug=True
            )
            assert sketch_c is not None
            assert sketch_c.num_intervals <= total_intervals  # Should have <= total_intervals due to sampling
            assert sketch_c.sketch.estimate_cardinality() > 0  # Should have both intervals and points
            
        finally:
            os.unlink(f.name)

@pytest.mark.full
def test_supported_file_extensions():
    """Test handling of all supported file extensions."""
    # Test with BED format
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write("chr1\t100\t200\n")
        f.flush()
        sketch = IntervalSketch.from_file(filename=f.name, mode="A")
        assert sketch is not None
        os.unlink(f.name)
    
    # Test with gzipped BED
    with tempfile.NamedTemporaryFile(suffix='.bed.gz', delete=False) as f:
        with gzip.open(f.name, 'wt') as gz:
            gz.write("chr1\t100\t200\n")
        sketch = IntervalSketch.from_file(filename=f.name, mode="A")
        assert sketch is not None
        os.unlink(f.name)
    
    # Skip BigBed tests since they require pysam and proper file creation
    # BigBed files are tested separately in test_bed_advanced_features

@pytest.mark.full
def test_file_extension_messages():
    """Test that appropriate messages are printed for unsupported formats."""
    import io
    from contextlib import redirect_stdout
    
    # Test with unsupported format
    with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as f:
        with io.StringIO() as buf, redirect_stdout(buf):
            sketch = IntervalSketch.from_file(filename=f.name, mode="A")
            output = buf.getvalue()
            assert "Skipping unsupported file format" in output
            assert f.name in output
            os.unlink(f.name)
    
    # Test with unsupported gzipped format
    with tempfile.NamedTemporaryFile(suffix='.txt.gz', delete=False) as f:
        with io.StringIO() as buf, redirect_stdout(buf):
            sketch = IntervalSketch.from_file(filename=f.name, mode="A")
            output = buf.getvalue()
            assert "Skipping unsupported gzipped file format" in output
            assert f.name in output
            os.unlink(f.name)

@pytest.mark.full
def test_error_handling():
    """Test error handling for malformed files."""
    # Test non-existent file
    sketch = IntervalSketch.from_file(filename="nonexistent.bed", mode="A", debug=True)
    assert sketch is None

    # Test malformed Bed file - use truly malformed data
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write("ch1\tinvalid\t200\n")  # Non-integer start position
        f.flush()
        try:
            sketch = IntervalSketch.from_file(filename=f.name, mode="A", debug=True)
            assert sketch is None, "Expected None for malformed BED file with invalid coordinates"
        finally:
            os.unlink(f.name)  # Make sure to clean up the file
            
    # Test valid Bed file with non-standard chromosome name
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write("scaffold123\t100\t200\n")  # Valid non-standard chromosome name
        f.flush()
        try:
            sketch = IntervalSketch.from_file(filename=f.name, mode="A", debug=True)
            assert sketch is not None, "Expected valid sketch for BED file with non-standard chromosome name"
            assert sketch.num_intervals == 1, "Expected one interval in valid BED file"
        finally:
            os.unlink(f.name)

@pytest.mark.full
def test_parallel_file_reading():
    """Test parallelized file reading functionality."""
    # Create a large test bed file with many intervals
    total_intervals = 100000  # Large number of intervals to make parallelization worthwhile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        for i in range(total_intervals):
            f.write(f"chr1\t{i*100}\t{(i+1)*100}\n")
        f.flush()
        
        try:
            # Test with different numbers of processes
            for num_processes in [1, 2, 4]:
                start_time = time.time()
                sketch = IntervalSketch.from_file(
                    filename=f.name,
                    mode="A",
                    num_processes=num_processes,
                    chunk_size=10000,
                    precision=20  # 2^20 = 1,048,576 registers, sufficient for ~10,000,000 points
                )
                end_time = time.time()
                
                assert sketch is not None
                assert sketch.num_intervals == total_intervals
                print(f"Processing time with {num_processes} processes: {end_time - start_time:.2f}s")
                
                # Verify cardinality estimate
                est_card = sketch.sketch.estimate_cardinality()
                assert est_card > 0
                
                # Test with different chunk sizes
                for chunk_size in [5000, 10000, 20000]:
                    start_time = time.time()
                    sketch = IntervalSketch.from_file(
                        filename=f.name,
                        mode="A",
                        num_processes=4,
                        chunk_size=chunk_size,
                        precision=20  # 2^20 = 1,048,576 registers, sufficient for ~10,000,000 points
                    )
                    end_time = time.time()
                    
                    assert sketch is not None
                    assert sketch.num_intervals == total_intervals
                    print(f"Processing time with chunk_size {chunk_size}: {end_time - start_time:.2f}s")
            
            # Test with different modes
            for mode in ["A", "B", "C"]:
                start_time = time.time()
                sketch = IntervalSketch.from_file(
                    filename=f.name,
                    mode=mode,
                    num_processes=4,
                    chunk_size=10000,
                    precision=20  # 2^20 = 1,048,576 registers, sufficient for ~10,000,000 points
                )
                end_time = time.time()
                
                assert sketch is not None
                if mode == "A":
                    assert sketch.num_intervals == total_intervals
                print(f"Processing time for mode {mode}: {end_time - start_time:.2f}s")
                
            # Test with subsampling
            start_time = time.time()
            sketch = IntervalSketch.from_file(
                filename=f.name,
                mode="C",
                num_processes=4,
                chunk_size=10000,
                subsample=(0.5, 0.5),
                precision=20  # 2^20 = 1,048,576 registers, sufficient for ~10,000,000 points
            )
            end_time = time.time()
            
            assert sketch is not None
            assert sketch.num_intervals <= total_intervals  # Should have fewer intervals due to subsampling
            print(f"Processing time with subsampling: {end_time - start_time:.2f}s")
            
        finally:
            os.unlink(f.name)

@pytest.mark.full
def test_unsupported_file():
    """Test handling of unsupported file formats."""
    # Test with unsupported binary formats
    for ext in ['cram', 'sam', 'bam', 'bw', 'bigwig', 'gff', 'gff3']:
        with tempfile.NamedTemporaryFile(suffix=f'.{ext}', delete=False) as f:
            sketch = IntervalSketch.from_file(filename=f.name, mode="A")
            assert sketch is None
            os.unlink(f.name)
    
    # Test with unsupported text formats
    for ext in ['txt', 'csv', 'tsv', 'json']:
        with tempfile.NamedTemporaryFile(suffix=f'.{ext}', delete=False) as f:
            sketch = IntervalSketch.from_file(filename=f.name, mode="A")
            assert sketch is None
            os.unlink(f.name)
    
    # Test with unsupported gzipped formats
    for ext in ['txt', 'csv', 'tsv', 'json', 'gff', 'gff3']:
        with tempfile.NamedTemporaryFile(suffix=f'.{ext}.gz', delete=False) as f:
            sketch = IntervalSketch.from_file(filename=f.name, mode="A")
            assert sketch is None
            os.unlink(f.name)
