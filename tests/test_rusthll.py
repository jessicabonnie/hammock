import pytest
import numpy as np
import time
import signal
import os
from hammock.lib.rusthll import RustHyperLogLog
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.rusthll_compat import RustHLLWrapper, RUST_AVAILABLE

class TimeoutError(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutError("Function took too long to complete")

def run_with_timeout(func, timeout=5):
    """Run a function with a timeout."""
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(timeout)
    try:
        result = func()
        signal.alarm(0)
        return result
    except TimeoutError as e:
        print(f"TIMEOUT: {e}")
        return None

def test_rust_import():
    """Test that Rust implementation is available."""
    assert RUST_AVAILABLE, "Rust implementation should be available"

def test_rusthll_initialization():
    """Test initialization of RustHyperLogLog with different precisions."""
    for precision in range(4, 25):
        sketch = RustHyperLogLog(precision)
        assert sketch.precision == precision
        assert sketch.num_registers == 2 ** precision
        assert sketch.is_using_rust()  # Check that Rust implementation is being used

def test_rusthll_add():
    """Test adding items to RustHyperLogLog."""
    sketch = RustHyperLogLog(14)
    items = ["item1", "item2", "item3"]
    
    for item in items:
        sketch.add_string(item)
    
    # Test adding a batch
    sketch.add_batch(items)
    
    # Test adding integers
    sketch.add(42)
    sketch.add_batch([1, 2, 3])

def test_rusthll_cardinality():
    """Test cardinality estimation."""
    # For small sets, use lower precision
    sketch = RustHyperLogLog(expected_cardinality=1000)
    items = [f"item{i}" for i in range(1000)]
    
    # Add items
    sketch.add_batch(items)
    
    # Get estimate
    estimate = sketch.estimate_cardinality()
    
    # Should be close to actual cardinality
    assert abs(estimate - 1000) / 1000 < 0.1  # Within 10% error

def test_rusthll_merge():
    """Test merging two RustHyperLogLog sketches."""
    # For medium sets, use medium precision
    sketch1 = RustHyperLogLog(expected_cardinality=5000)
    sketch2 = RustHyperLogLog(expected_cardinality=5000)
    
    # Add different items to each sketch
    items1 = [f"item1_{i}" for i in range(500)]
    items2 = [f"item2_{i}" for i in range(500)]
    
    sketch1.add_batch(items1)
    sketch2.add_batch(items2)
    
    # Merge sketches
    sketch1.merge(sketch2)
    
    # Estimate should be close to union cardinality
    estimate = sketch1.estimate_cardinality()
    assert abs(estimate - 1000) / 1000 < 0.1  # Within 10% error

def test_rusthll_intersection():
    """Test intersection estimation."""
    sketch1 = RustHyperLogLog(12)
    sketch2 = RustHyperLogLog(12)
    
    # Add overlapping items
    items1 = [f"item_{i}" for i in range(10000)]
    items2 = [f"item_{i}" for i in range(5000, 15000)]  # 5000 items overlap
    
    sketch1.add_batch(items1)
    sketch2.add_batch(items2)
    
    # Estimate intersection
    intersection = sketch1.estimate_intersection(sketch2)
    assert abs(intersection - 5000) / 5000 < 0.2  # Within 20% error

def test_rusthll_union():
    """Test union estimation."""
    sketch1 = RustHyperLogLog(12)
    sketch2 = RustHyperLogLog(12)
    
    # Add overlapping items
    items1 = [f"item_{i}" for i in range(10000)]
    items2 = [f"item_{i}" for i in range(5000, 15000)]  # 5000 items overlap
    
    sketch1.add_batch(items1)
    sketch2.add_batch(items2)
    
    # Estimate union
    union = sketch1.estimate_union(sketch2)
    assert abs(union - 15000) / 15000 < 0.2  # Within 20% error

def test_rusthll_jaccard():
    """Test Jaccard similarity estimation."""
    sketch1 = RustHyperLogLog(12)
    sketch2 = RustHyperLogLog(12)
    
    # Add overlapping items
    items1 = [f"item_{i}" for i in range(10000)]
    items2 = [f"item_{i}" for i in range(5000, 15000)]  # 5000 items overlap
    
    sketch1.add_batch(items1)
    sketch2.add_batch(items2)
    
    # Calculate Jaccard similarity
    jaccard = sketch1.estimate_jaccard(sketch2)
    expected_jaccard = 5000 / 15000  # intersection / union
    assert abs(jaccard - expected_jaccard) < 0.1  # Within 0.1 error

def test_rusthll_serialization():
    """Test saving and loading sketches."""
    import tempfile
    
    sketch = RustHyperLogLog(14)
    items = [f"item_{i}" for i in range(1000)]
    sketch.add_batch(items)
    
    # Create a temporary directory for the test
    with tempfile.TemporaryDirectory() as temp_dir:
        # Use a temporary file with .npy extension
        test_file = os.path.join(temp_dir, "test_sketch.npy")
        os.makedirs(os.path.dirname(test_file), exist_ok=True)
        sketch.save(test_file)
        
        # Load sketch
        loaded_sketch = RustHyperLogLog.load(test_file)
        
        # Compare estimates
        original_estimate = sketch.estimate_cardinality()
        loaded_estimate = loaded_sketch.estimate_cardinality()
        assert abs(original_estimate - loaded_estimate) / original_estimate < 0.01

def test_rusthll_precision_compatibility():
    """Test that sketches with different precisions cannot be merged."""
    sketch1 = RustHyperLogLog(14)
    sketch2 = RustHyperLogLog(16)
    
    with pytest.raises(ValueError):
        sketch1.merge(sketch2)

def test_rusthll_wrapper():
    """Test the RustHLLWrapper class."""
    if not RUST_AVAILABLE:
        pytest.skip("Rust implementation not available")
    
    wrapper = RustHLLWrapper(14)
    items = [f"item_{i}" for i in range(1000)]
    
    # Test adding items
    for item in items:
        wrapper.add(item)
    
    # Test batch adding
    wrapper.add_batch(items)
    
    # Test cardinality estimation
    estimate = wrapper.estimate_cardinality()
    assert abs(estimate - 1000) / 1000 < 0.1

def test_rust_flag():
    """Test that the --rust flag correctly selects the Rust implementation."""
    from hammock.hammock import main
    import sys
    import io
    import tempfile
    
    # Create temporary test files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write("chr1\t0\t100\n")
        test_file = f.name
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write("chr1\t0\t100\n")
        primary_file = f.name
    
    # Create temporary files for filepaths and primary files
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        f.write(test_file + "\n")
        filepaths_file = f.name
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        f.write(primary_file + "\n")
        primary_file_list = f.name
    
    # Capture stdout
    stdout = io.StringIO()
    sys.stdout = stdout
    
    # Test with --rust flag
    sys.argv = ["hammock", "--rust", "--hyperloglog", "--mode", "A", "--precision", "12", filepaths_file, primary_file_list]
    try:
        main()
    except SystemExit as e:
        assert e.code == 0  # Should succeed
    
    # Check that Rust implementation was used
    output = stdout.getvalue()
    assert "Using Rust implementation" in output
    
    # Clean up
    os.unlink(test_file)
    os.unlink(primary_file)
    os.unlink(filepaths_file)
    os.unlink(primary_file_list)

def test_rusthll_large_batch():
    """Test handling of very large batches."""
    if not RUST_AVAILABLE:
        pytest.skip("Rust implementation not available")
    
    # For very large sets, use higher precision
    sketch = RustHyperLogLog(expected_cardinality=1000000)
    large_values = [f"large_item_{i}" for i in range(1000000)]
    
    # Test adding large batch
    sketch.add_batch(large_values)
    
    # Get estimate
    estimate = sketch.estimate_cardinality()
    assert abs(estimate - 1000000) / 1000000 < 0.1

def test_rusthll_basic_operations():
    """Test basic operations of RustHyperLogLog."""
    # For very small sets, use lowest precision
    sketch = RustHyperLogLog(expected_cardinality=10)
    assert sketch.precision == 4
    
    # Test adding strings
    sketch.add("test1")
    sketch.add("test2")
    sketch.add("test3")
    
    
    # Test cardinality estimation
    est = sketch.estimate_cardinality()
    assert est > 0
    assert est < 5  # Should be close to 3, allow for some error
    
    # Test batch operations
    strings = ["batch1", "batch2", "batch3"]
    sketch.add_batch(strings)
    est = sketch.estimate_cardinality()
    assert est > 0
    # assert est < 8  # Should be close to 6, allow for some error
    
    # Test merging
    sketch1 = RustHyperLogLog(precision=4)
    sketch2 = RustHyperLogLog(precision=4)
    
    sketch1.add("merge1")
    sketch1.add("merge2")
    
    sketch2.add("merge2")
    sketch2.add("merge3")
    
    sketch1.merge(sketch2)
    est = sketch1.estimate_cardinality()
    assert est > 0
    assert est < 5  # Should be close to 3, allow for some error
    
    # Test file operations
    sketch1 = RustHyperLogLog(precision=8)
    sketch2 = RustHyperLogLog(precision=8)
    sketch3 = RustHyperLogLog(precision=8)
    
    sketch1.add("file1")
    sketch2.add("file2")
    sketch3.add("file3")
    
    # Test writing and loading
    test_file = os.path.join("test_results", "test_sketch1.npy")
    os.makedirs("test_results", exist_ok=True)
    try:
        sketch1.save(test_file)
        loaded_sketch = RustHyperLogLog.load(test_file)
        assert abs(sketch1.estimate_cardinality() - loaded_sketch.estimate_cardinality()) < 0.2
    finally:
        if os.path.exists(test_file):
            os.remove(test_file)
    
    # Test precision validation
    with pytest.raises(ValueError):
        RustHyperLogLog(precision=3)  # Too low
    with pytest.raises(ValueError):
        RustHyperLogLog(precision=25)  # Too high
        
    # Test merge validation
    sketch1 = RustHyperLogLog(precision=14)
    sketch2 = RustHyperLogLog(precision=15)  # Different precision
    with pytest.raises(ValueError):
        sketch1.merge(sketch2)
        
    # Test empty sketch
    sketch = RustHyperLogLog(precision=14)
    assert sketch.is_empty()
    sketch.add("test")
    assert not sketch.is_empty()
    
    # Test Jaccard similarity
    sketch1 = RustHyperLogLog(precision=14)
    sketch2 = RustHyperLogLog(precision=14)
    
    sketch1.add("jaccard1")
    sketch1.add("jaccard2")
    
    sketch2.add("jaccard2")
    sketch2.add("jaccard3")
    
    jaccard = sketch1.estimate_jaccard(sketch2)
    assert jaccard > 0
    assert jaccard < 1
    
    # Test retry with higher precision
    if not RustHyperLogLog(precision=14).is_using_rust():
        # Create a sketch with low precision
        rust_sketch = RustHyperLogLog(precision=14)
        rust_sketch.add("test1")
        rust_sketch.add("test2")
        
        # Create a sketch with higher precision
        high_precision_sketch = RustHyperLogLog(precision=rust_sketch.precision + 2)
        high_precision_sketch.add("test1")
        high_precision_sketch.add("test2")
        
        # Compare estimates
        low_est = rust_sketch.estimate_cardinality()
        high_est = high_precision_sketch.estimate_cardinality()
        
        # Higher precision should give more accurate estimate
        assert high_est > 0
        assert high_est < 3  # Should be close to 2 