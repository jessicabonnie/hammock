import pytest # type: ignore
from hammock.lib.intervals import IntervalSketch
from hammock.lib.abstractsketch import AbstractSketch
import tempfile
import os

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

@pytest.mark.quick
class TestIntervalSketchQuick:
    """Quick tests for IntervalSketch class."""
    
    def test_init(self):
        """Test basic initialization."""
        sketch = IntervalSketch(mode="A", sketch_type="hyperloglog")
        assert sketch.mode == "A"
        assert sketch.sketch is not None
    
    def test_bedline_mode_A(self):
        """Test mode A (intervals only) with small input"""
        sketch = IntervalSketch(mode="A")
        interval, points, size = sketch.bedline("chr1\t100\t200", mode="A", sep="-")
        assert interval == b"1-100-200-A"
        assert points == []
        assert size == 100

    def test_bedline_mode_B(self):
        """Test mode B (points only) with small input"""
        sketch = IntervalSketch(mode="B")
        interval, points, size = sketch.bedline("chr1\t100\t103", mode="B", sep="-")
        assert interval is None
        assert len(points) == 3
        assert size == 3

    def test_bedline_mode_C(self):
        """Test mode C (both) with small input"""
        sketch = IntervalSketch(mode="C")
        interval, points, size = sketch.bedline("chr1\t100\t103", mode="C", sep="-")
        assert interval == b"1-100-103-A"
        assert len(points) == 3
        assert size == 3

    def test_from_file_small(self, small_bed_file):
        """Test file processing with small file"""
        sketch = IntervalSketch.from_file(
            filename=small_bed_file,
            mode="A",
            precision=8,
            sketch_type="hyperloglog"
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
        assert interval == b"1-100-103-A"
        assert len(points) == 3
        
        # Test with interval subsampling only
        interval, points, size = sketch.bedline(
            "chr1\t100\t103", 
            mode="C", 
            sep="-",
            subsample=(0.0, 1.0)
        )
        assert interval is None  # Should be subsampled out
        assert len(points) == 3
        
        # Test with point subsampling only
        interval, points, size = sketch.bedline(
            "chr1\t100\t103", 
            mode="C", 
            sep="-",
            subsample=(1.0, 0.0)
        )
        assert interval == b"1-100-103-A"
        assert all(p is None for p in points)  # All points should be subsampled out

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
                subsample=(subA, subB)
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
        assert len(sketch.generate_points("chr1", 1000, 2000, subsample=0.5)) < 1001

if __name__ == "__main__":
    pytest.main([__file__])