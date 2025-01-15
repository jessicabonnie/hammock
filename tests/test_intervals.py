import pytest # type: ignore
from hammock.lib.intervals import IntervalSketch
from hammock.lib.abstractsketch import AbstractSketch
import tempfile
import os
import subprocess

@pytest.fixture
def small_bed_file():
    """Fixture for small BED file used in quick tests"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        f.write("chr1\t100\t200\n")
        f.write("chr2\t150\t300\n")
        f.write("chr3\t1000\t2000\n")
    yield f.name
    os.unlink(f.name)

@pytest.fixture
def large_bed_file():
    """Fixture for large BED file used in full tests"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        # Write 1000 intervals
        for i in range(1000):
            f.write(f"chr1\t{i*1000}\t{(i+1)*1000}\n")
    yield f.name
    os.unlink(f.name)

@pytest.fixture
def small_bigbed_file():
    """Fixture for small BigBed file used in quick tests"""
    # Create a temporary BED file
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        f.write("chr1\t100\t200\n")
        f.write("chr2\t150\t300\n")
        f.write("chr3\t1000\t2000\n")
        bed_file = f.name

    # Convert to BigBed using bedToBigBed
    bb_file = bed_file + '.bb'
    try:
        # Check if bedToBigBed is available
        result = subprocess.run(['which', 'bedToBigBed'], 
                              capture_output=True, 
                              text=True)
        if result.returncode != 0:
            pytest.skip("bedToBigBed not found in PATH")
            
        # Create chrom.sizes file
        chrom_sizes = bed_file + '.sizes'
        with open(chrom_sizes, 'w') as f:
            f.write("chr1\t1000000\n")
            f.write("chr2\t1000000\n")
            f.write("chr3\t1000000\n")
            
        # Convert to BigBed
        subprocess.run(['bedToBigBed', bed_file, chrom_sizes, bb_file], 
                      check=True,
                      capture_output=True)
        os.unlink(chrom_sizes)
        
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        pytest.skip(f"bedToBigBed conversion failed: {str(e)}")
        
    yield bb_file
    
    # Cleanup
    try:
        os.unlink(bed_file)
        if os.path.exists(bb_file):
            os.unlink(bb_file)
    except OSError:
        pass  # Ignore cleanup errors

@pytest.mark.quick
class TestIntervalSketchQuick:
    """Quick tests for IntervalSketch class."""
    
    def test_init(self):
        """Test basic initialization."""
        sketch = IntervalSketch(mode="A")
        assert sketch.mode == "A"
        assert sketch.sketch is not None
        assert sketch.expA == 0
        assert sketch.subsample == (1.0, 1.0)
    
    def test_bedline_mode_A(self):
        """Test mode A (intervals only) with small input"""
        sketch = IntervalSketch(mode="A")
        interval, points, size = sketch.bedline("chr1\t100\t200", mode="A", sep="-")
        assert interval == "1-100-200-A"
        assert points == []
        assert size == 100

    def test_bedline_mode_B(self):
        """Test mode B (points only) with small input"""
        sketch = IntervalSketch(mode="B")
        interval, points, size = sketch.bedline("chr1\t100\t103", mode="B", sep="-")
        print(len(points))
        print([p for p in points])
        assert interval is None
        assert len(points) == 3
        assert all(isinstance(p, str) or p is None for p in points)
        assert size == 3

    def test_bedline_mode_C(self):
        """Test mode C (both) with small input"""
        sketch = IntervalSketch(mode="C")
        interval, points, size = sketch.bedline("chr1\t100\t103", mode="C", sep="-")
        if interval is not None:
            assert isinstance(interval, str)
        assert len(points) == 3
        assert all(isinstance(p, str) or p is None for p in points)
        assert size == 3

    def test_from_file_small(self, small_bed_file):
        """Test file processing with small file"""
        sketch = IntervalSketch.from_file(
            filename=small_bed_file,
            mode="A",
            precision=8,
            sketch_type="hyperloglog",
            expA=0
        )
        assert isinstance(sketch, IntervalSketch)
        assert sketch.num_intervals == 3
        assert sketch.total_interval_size == 1250

    def test_add_string(self):
        """Test adding strings to sketch"""
        sketch = IntervalSketch(mode="A")
        sketch.add_string("test1")
        sketch.add_string("test2")
        # Just verify it doesn't raise an exception
        assert True

    def test_bedline_mode_C_with_subsampling(self):
        """Test mode C with different subsampling rates"""
        sketch = IntervalSketch(mode="C")
        
        # Test with no subsampling
        interval, points, size = sketch.bedline(
            "chr1\t100\t103", 
            mode="C", 
            sep="-",
            subsample=(1.0, 1.0)
        )
        if interval is not None:
            assert isinstance(interval, str)
        assert len(points) == 3
        assert all(isinstance(p, str) or p is None for p in points)
        
        # Test with interval subsampling only
        interval, points, size = sketch.bedline(
            "chr1\t100\t103", 
            mode="C", 
            sep="-",
            subsample=(0.0, 1.0)
        )
        assert interval is None
        assert len(points) == 3
        assert all(isinstance(p, str) or p is None for p in points)
        
        # Test with point subsampling only
        interval, points, size = sketch.bedline(
            "chr1\t100\t103", 
            mode="C", 
            sep="-",
            subsample=(1.0, 0.0)
        )
        if interval is not None:
            assert isinstance(interval, str)
        assert all(p is None for p in points)

    def test_bigbed_processing(self, small_bigbed_file):
        """Test processing of BigBed files"""
        if small_bigbed_file is None:
            pytest.skip("BigBed file creation failed")
        
        sketch = IntervalSketch.from_file(
            filename=small_bigbed_file,
            mode="A",
            precision=8
        )
        assert sketch is not None
        assert sketch.num_intervals == 3

    def test_multiplicity_constraints(self):
        """Test that multiplicity constraints are enforced."""
        # Create a temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write("chr1\t100\t200\n")  # This creates an interval of size 100
            temp_file = f.name

        try:
            # Test that multiplicity only works in mode C
            with pytest.raises(ValueError, match="Multiplicity .* mode C"):
                IntervalSketch.from_file(
                    filename=temp_file,
                    mode="A",
                    expA=2
                )
            
            with pytest.raises(ValueError, match="Multiplicity .* mode C"):
                IntervalSketch.from_file(
                    filename=temp_file,
                    mode="B",
                    expA=2
                )

            # Test that multiplicity cannot be used with subsampling
            with pytest.raises(ValueError, match="Multiplicity .* subsampling"):
                IntervalSketch.from_file(
                    filename=temp_file,
                    mode="C",
                    expA=2,
                    subsample=(0.5, 1.0)
                )

            # Test that valid usage works
            sketch = IntervalSketch.from_file(
                filename=temp_file,
                mode="C",
                expA=2,
                subsample=(1.0, 1.0),
                sketch_type="exacttest"  # Use exacttest instead of exact
            )
            assert sketch is not None
            assert sketch.num_intervals == 1
            cardinality = sketch.sketch.estimate_cardinality()
            # Expected: 1 original interval + 99 copies + 100 points = 200
            expected = 200
            assert cardinality == expected, f"Expected cardinality {expected}, got {cardinality}"

        finally:
            os.unlink(temp_file)

@pytest.mark.full
class TestIntervalSketchFull:
    """Full tests for IntervalSketch class."""
    
    def test_large_file_processing(self, large_bed_file):
        """Test processing of large BED file"""
        sketch = IntervalSketch.from_file(
            filename=large_bed_file,
            mode="A",
            precision=12,
            sketch_type="hyperloglog"
        )
        assert sketch is not None
        assert sketch.num_intervals == 1000
        assert sketch.total_interval_size == 1000000

    def test_subsampling_modes(self, large_bed_file):
        """Test different subsampling rate combinations"""
        test_cases = [
            (1.0, 1.0),  # No subsampling
            (0.5, 1.0),  # Interval subsampling only
            (1.0, 0.5),  # Point subsampling only
            (0.5, 0.5),  # Both subsampled
            (0.0, 1.0),  # No intervals
            (1.0, 0.0),  # No points
        ]
        
        for subA, subB in test_cases:
            sketch = IntervalSketch.from_file(
                filename=large_bed_file,
                mode="C",
                subsample=(subA, subB),
                expA=0
            )
            assert sketch is not None
            
            # For complete subsampling, verify no elements
            if subA == 0.0:
                assert sketch.num_intervals == 0
            if subB == 0.0 and subA == 0.0:
                assert sketch.sketch.estimate_cardinality() == 0

    def test_mode_specific_subsampling(self, large_bed_file):
        """Test that subsampling only applies in mode C"""
        # Mode A should ignore subsampling
        sketch_A = IntervalSketch.from_file(
            filename=large_bed_file,
            mode="A",
            subsample=(0.0, 0.0)  # Should be ignored
        )
        assert sketch_A.num_intervals > 0
        
        # Mode B should ignore subsampling
        sketch_B = IntervalSketch.from_file(
            filename=large_bed_file,
            mode="B",
            subsample=(0.0, 0.0)  # Should be ignored
        )
        assert sketch_B.num_intervals > 0
        
        # Mode C should respect subsampling
        sketch_C = IntervalSketch.from_file(
            filename=large_bed_file,
            mode="C",
            subsample=(0.0, 0.0)
        )
        assert sketch_C.num_intervals == 0

    def test_generate_points_large(self):
        """Test point generation with large intervals"""
        sketch = IntervalSketch(mode="B")
        points = sketch.generate_points("chr1", 1000, 2000, subsample=0.5)
        assert len([p for p in points if p is not None]) < 1001
        assert all(isinstance(p, str) or p is None for p in points)

    def test_chunk_processing(self, large_bed_file):
        """Test chunk processing with large BED file"""
        sketch = IntervalSketch.from_file(
            filename=large_bed_file,
            mode="A",
            precision=12,
            sketch_type="hyperloglog"
        )
        assert sketch is not None
        assert sketch.num_intervals == 1000
        assert sketch.total_interval_size == 1000000

if __name__ == "__main__":
    pytest.main([__file__])