import pytest
from hammock.lib.intervals import IntervalSketch
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
    """Quick tests for IntervalSketch class"""
    
    def test_init(self):
        """Test basic initialization"""
        sketch = IntervalSketch(mode="A")
        assert sketch.mode == "A"
        assert sketch.total_interval_size == 0
        assert sketch.num_intervals == 0

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
        interval, points, size = sketch.bedline("chr1\t100\t102", mode="B", sep="-")
        assert interval is None
        assert len(points) == 3
        assert size == 2

    def test_bedline_mode_C(self):
        """Test mode C (both) with small input"""
        sketch = IntervalSketch(mode="C")
        interval, points, size = sketch.bedline("chr1\t100\t102", mode="C", sep="-")
        assert interval == b"1-100-102-A"
        assert len(points) == 3
        assert size == 2

    def test_from_file_small(self, small_bed_file):
        """Test file processing with small file"""
        sketch = IntervalSketch.from_file(
            filename=small_bed_file,
            mode="A",
            precision=8,
            sketch_type="hyperloglog"
        )
        assert sketch is not None
        assert sketch.num_intervals == 3
        assert sketch.total_interval_size == 1250

    def test_add_string(self):
        """Test adding strings to sketch"""
        sketch = IntervalSketch(mode="A")
        sketch.add_string("test1")
        sketch.add_string("test2")
        # Just verify it doesn't raise an exception
        assert True

@pytest.mark.full
class TestIntervalSketchFull:
    """Full test suite for IntervalSketch class"""
    
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
        """Test different subsampling rates"""
        rates = [-0.5, 0.5, 0.1]
        for rate in rates:
            sketch = IntervalSketch.from_file(
                filename=large_bed_file,
                mode="C",
                subsample=rate
            )
            assert sketch is not None

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