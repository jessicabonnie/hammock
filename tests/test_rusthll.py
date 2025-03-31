import pytest
import numpy as np
import time
import signal
import os
from hammock.lib.rusthll import RustHLL
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

def test_initialization():
    """Test initialization of RustHLL with different precisions."""
    for precision in [4, 8, 12, 16]:
        sketch = RustHLL(precision)
        assert sketch.precision == precision

def test_add_items():
    """Test adding items to RustHLL."""
    sketch = RustHLL(14)
    items = ["item1", "item2", "item3"]
    for item in items:
        sketch.add(item)
    assert not sketch.is_empty()

def test_expected_cardinality():
    """Test initialization with expected cardinality."""
    sketch = RustHLL(expected_cardinality=1000)
    assert sketch.precision >= 8  # Should adjust precision based on expected cardinality

def test_merge():
    """Test merging two RustHLL sketches."""
    sketch1 = RustHLL(expected_cardinality=5000)
    sketch2 = RustHLL(expected_cardinality=5000)
    
    # Add different items to each sketch
    for i in range(1000):
        sketch1.add(f"item1_{i}")
        sketch2.add(f"item2_{i}")
    
    # Merge sketches
    sketch1.merge(sketch2)
    
    # Check that merged sketch has higher cardinality
    assert sketch1.estimate_cardinality() > 1000

def test_precision_validation():
    """Test precision validation."""
    # Test valid precisions
    for precision in [4, 8, 12, 16]:
        sketch = RustHLL(precision)
        assert sketch.precision == precision
    
    # Test invalid precisions
    with pytest.raises(ValueError):
        RustHLL(precision=3)  # Too low

def test_merge_different_precision():
    """Test merging sketches with different precision."""
    sketch1 = RustHLL(precision=14)
    sketch2 = RustHLL(precision=15)  # Different precision
    
    with pytest.raises(ValueError):
        sketch1.merge(sketch2)

def test_empty_sketch():
    """Test empty sketch operations."""
    sketch = RustHLL(precision=14)
    assert sketch.is_empty()
    assert sketch.estimate_cardinality() == 0

def test_rust_implementation():
    """Test Rust implementation availability."""
    if not RustHLL(precision=14).is_using_rust():
        pytest.skip("Rust implementation not available")
    
    rust_sketch = RustHLL(precision=14)
    assert rust_sketch.is_using_rust()
    
    # Test with high precision
    high_precision_sketch = RustHLL(precision=rust_sketch.precision + 2)
    assert high_precision_sketch.is_using_rust()

def test_rusthll_cardinality():
    """Test cardinality estimation."""
    # For small sets, use lower precision
    sketch = RustHLL(expected_cardinality=1000)
    items = [f"item{i}" for i in range(1000)]
    
    # Add items
    sketch.add_batch(items)
    
    # Get estimate
    estimate = sketch.estimate_cardinality()
    
    # Should be close to actual cardinality
    assert abs(estimate - 1000) / 1000 < 0.1  # Within 10% error

def test_rusthll_intersection():
    """Test intersection estimation."""
    sketch1 = RustHLL(12)
    sketch2 = RustHLL(12)
    
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
    sketch1 = RustHLL(12)
    sketch2 = RustHLL(12)
    
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
    sketch1 = RustHLL(12)
    sketch2 = RustHLL(12)
    
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
    
    sketch = RustHLL(14)
    items = [f"item_{i}" for i in range(1000)]
    sketch.add_batch(items)
    
    # Create a temporary directory for the test
    with tempfile.TemporaryDirectory() as temp_dir:
        # Use a temporary file with .npy extension
        test_file = os.path.join(temp_dir, "test_sketch.npy")
        os.makedirs(os.path.dirname(test_file), exist_ok=True)
        sketch.save(test_file)
        
        # Load sketch
        loaded_sketch = RustHLL.load(test_file)
        
        # Compare estimates
        original_estimate = sketch.estimate_cardinality()
        loaded_estimate = loaded_sketch.estimate_cardinality()
        assert abs(original_estimate - loaded_estimate) / original_estimate < 0.01

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
    sketch = RustHLL(expected_cardinality=1000000)
    large_values = [f"large_item_{i}" for i in range(1000000)]
    
    # Test adding large batch
    sketch.add_batch(large_values)
    
    # Get estimate
    estimate = sketch.estimate_cardinality()
    assert abs(estimate - 1000000) / 1000000 < 0.1

def test_rusthll_basic_operations():
    """Test basic operations of RustHLL."""
    # For very small sets, use lowest precision
    sketch = RustHLL(expected_cardinality=10)
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
    sketch1 = RustHLL(precision=4)
    sketch2 = RustHLL(precision=4)
    
    sketch1.add("merge1")
    sketch1.add("merge2")
    
    sketch2.add("merge2")
    sketch2.add("merge3")
    
    sketch1.merge(sketch2)
    est = sketch1.estimate_cardinality()
    assert est > 0
    assert est < 5  # Should be close to 3, allow for some error
    
    # Test file operations
    sketch1 = RustHLL(precision=8)
    sketch2 = RustHLL(precision=8)
    sketch3 = RustHLL(precision=8)
    
    sketch1.add("file1")
    sketch2.add("file2")
    sketch3.add("file3")
    
    # Test writing and loading
    test_file = os.path.join("test_results", "test_sketch1.npy")
    os.makedirs("test_results", exist_ok=True)
    try:
        sketch1.save(test_file)
        loaded_sketch = RustHLL.load(test_file)
        assert abs(sketch1.estimate_cardinality() - loaded_sketch.estimate_cardinality()) < 0.2
    finally:
        if os.path.exists(test_file):
            os.remove(test_file)
    
    # Test precision validation
    with pytest.raises(ValueError):
        RustHLL(precision=3)  # Too low
        
    # Test merge validation
    sketch1 = RustHLL(precision=14)
    sketch2 = RustHLL(precision=15)  # Different precision
    with pytest.raises(ValueError):
        sketch1.merge(sketch2)
        
    # Test empty sketch
    sketch = RustHLL(precision=14)
    assert sketch.is_empty()
    sketch.add("test")
    assert not sketch.is_empty()
    
    # Test Jaccard similarity
    sketch1 = RustHLL(precision=14)
    sketch2 = RustHLL(precision=14)
    
    sketch1.add("jaccard1")
    sketch1.add("jaccard2")
    
    sketch2.add("jaccard2")
    sketch2.add("jaccard3")
    
    jaccard = sketch1.estimate_jaccard(sketch2)
    assert jaccard > 0
    assert jaccard < 1
    
    # Test retry with higher precision
    if not RustHLL(precision=14).is_using_rust():
        # Create a sketch with low precision
        rust_sketch = RustHLL(precision=14)
        rust_sketch.add("test1")
        rust_sketch.add("test2")
        
        # Create a sketch with higher precision
        high_precision_sketch = RustHLL(precision=rust_sketch.precision + 2)
        high_precision_sketch.add("test1")
        high_precision_sketch.add("test2")
        
        # Compare estimates
        low_est = rust_sketch.estimate_cardinality()
        high_est = high_precision_sketch.estimate_cardinality()
        
        # Higher precision should give more accurate estimate
        assert high_est > 0
        assert high_est < 3  # Should be close to 2

def test_high_precision_support():
    """Test support for high precision values (over 16)."""
    # Test with high precision values
    high_precision_values = [18, 20, 22, 24]
    
    for precision in high_precision_values:
        # Test Python implementation
        py_sketch = HyperLogLog(precision=precision, hash_size=32)
        assert py_sketch.precision == precision
        
        # Test Rust implementation
        rust_sketch = RustHLL(precision=precision, hash_size=32)
        assert rust_sketch.precision == precision
        
        # Verify Rust implementation is actually being used
        if RUST_AVAILABLE:
            assert rust_sketch.is_using_rust(), f"RustHLL should use Rust implementation for precision={precision}"
        
        # Test adding items with high precision
        for i in range(100):
            rust_sketch.add(f"item_{i}")
        
        # Test estimation
        cardinality = rust_sketch.estimate_cardinality()
        assert 90 <= cardinality <= 110, f"Estimated cardinality {cardinality} should be close to 100"
        
        # Test with RustHLLWrapper if available
        if RUST_AVAILABLE:
            wrapper_sketch = RustHLLWrapper(precision=precision, hash_size=32)
            for i in range(100):
                wrapper_sketch.add(f"item_{i}")
            wrapper_cardinality = wrapper_sketch.estimate_cardinality()
            assert 90 <= wrapper_cardinality <= 110, f"Wrapper cardinality {wrapper_cardinality} should be close to 100" 