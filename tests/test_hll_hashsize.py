import pytest # type: ignore
import numpy as np # type: ignore
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.rusthll import RustHyperLogLog

def test_hash_size_initialization():
    """Test that both implementations can be initialized with both hash sizes."""
    # Test Python implementation
    hll_py_32 = HyperLogLog(precision=12, hash_size=32)
    hll_py_64 = HyperLogLog(precision=12, hash_size=64)
    
    # Test Rust implementation
    hll_rust_32 = RustHyperLogLog(precision=12, hash_size=32)
    hll_rust_64 = RustHyperLogLog(precision=12, hash_size=64)
    
    # Verify all instances were created successfully
    assert hll_py_32 is not None
    assert hll_py_64 is not None
    assert hll_rust_32 is not None
    assert hll_rust_64 is not None

def test_invalid_hash_size():
    """Test that invalid hash sizes are rejected."""
    with pytest.raises(ValueError, match="hash_size must be 32 or 64"):
        HyperLogLog(precision=12, hash_size=16)
    
    with pytest.raises(ValueError, match="hash_size must be 32 or 64"):
        RustHyperLogLog(precision=12, hash_size=16)

def test_basic_functionality():
    """Test basic functionality with both hash sizes."""
    # Create test data
    test_data = [f"item_{i}" for i in range(1000)]
    
    # Test Python implementation
    hll_py_32 = HyperLogLog(precision=12, hash_size=32)
    hll_py_64 = HyperLogLog(precision=12, hash_size=64)
    
    # Test Rust implementation
    hll_rust_32 = RustHyperLogLog(precision=12, hash_size=32)
    hll_rust_64 = RustHyperLogLog(precision=12, hash_size=64)
    
    # Add data to all sketches
    for item in test_data:
        hll_py_32.add_string(item)
        hll_py_64.add_string(item)
        hll_rust_32.add_string(item)
        hll_rust_64.add_string(item)
    
    # Get estimates
    est_py_32 = hll_py_32.estimate_cardinality()
    est_py_64 = hll_py_64.estimate_cardinality()
    est_rust_32 = hll_rust_32.estimate_cardinality()
    est_rust_64 = hll_rust_64.estimate_cardinality()
    
    # Verify estimates are reasonable
    assert 900 < est_py_32 < 1100
    assert 900 < est_py_64 < 1100
    assert 900 < est_rust_32 < 1100
    assert 900 < est_rust_64 < 1100

def test_batch_operations():
    """Test batch operations with both hash sizes."""
    # Create test data
    test_data = [f"item_{i}" for i in range(1000)]
    
    # Test Python implementation
    hll_py_32 = HyperLogLog(precision=12, hash_size=32)
    hll_py_64 = HyperLogLog(precision=12, hash_size=64)
    
    # Test Rust implementation
    hll_rust_32 = RustHyperLogLog(precision=12, hash_size=32)
    hll_rust_64 = RustHyperLogLog(precision=12, hash_size=64)
    
    # Add data in batches
    batch_size = 100
    for i in range(0, len(test_data), batch_size):
        batch = test_data[i:i + batch_size]
        hll_py_32.add_batch(batch)
        hll_py_64.add_batch(batch)
        hll_rust_32.add_batch(batch)
        hll_rust_64.add_batch(batch)
    
    # Get estimates
    est_py_32 = hll_py_32.estimate_cardinality()
    est_py_64 = hll_py_64.estimate_cardinality()
    est_rust_32 = hll_rust_32.estimate_cardinality()
    est_rust_64 = hll_rust_64.estimate_cardinality()
    
    # Verify estimates are reasonable
    assert 900 < est_py_32 < 1100
    assert 900 < est_py_64 < 1100
    assert 900 < est_rust_32 < 1100
    assert 900 < est_rust_64 < 1100

def test_merge_operations():
    """Test merge operations with both hash sizes."""
    # Create two sets of test data
    data1 = [f"item_{i}" for i in range(500)]
    data2 = [f"item_{i}" for i in range(500, 1000)]
    
    # Create sketches for first set
    hll_py_32_1 = HyperLogLog(precision=12, hash_size=32)
    hll_py_64_1 = HyperLogLog(precision=12, hash_size=64)
    hll_rust_32_1 = RustHyperLogLog(precision=12, hash_size=32)
    hll_rust_64_1 = RustHyperLogLog(precision=12, hash_size=64)
    
    # Create sketches for second set
    hll_py_32_2 = HyperLogLog(precision=12, hash_size=32)
    hll_py_64_2 = HyperLogLog(precision=12, hash_size=64)
    hll_rust_32_2 = RustHyperLogLog(precision=12, hash_size=32)
    hll_rust_64_2 = RustHyperLogLog(precision=12, hash_size=64)
    
    # Add data to respective sketches
    for item in data1:
        hll_py_32_1.add_string(item)
        hll_py_64_1.add_string(item)
        hll_rust_32_1.add_string(item)
        hll_rust_64_1.add_string(item)
    
    for item in data2:
        hll_py_32_2.add_string(item)
        hll_py_64_2.add_string(item)
        hll_rust_32_2.add_string(item)
        hll_rust_64_2.add_string(item)
    
    # Merge sketches
    hll_py_32_1.merge(hll_py_32_2)
    hll_py_64_1.merge(hll_py_64_2)
    hll_rust_32_1.merge(hll_rust_32_2)
    hll_rust_64_1.merge(hll_rust_64_2)
    
    # Get estimates
    est_py_32 = hll_py_32_1.estimate_cardinality()
    est_py_64 = hll_py_64_1.estimate_cardinality()
    est_rust_32 = hll_rust_32_1.estimate_cardinality()
    est_rust_64 = hll_rust_64_1.estimate_cardinality()
    
    # Verify estimates are reasonable
    assert 900 < est_py_32 < 1100
    assert 900 < est_py_64 < 1100
    assert 900 < est_rust_32 < 1100
    assert 900 < est_rust_64 < 1100

def test_edge_cases():
    """Test edge cases with both hash sizes."""
    # Test empty sketches
    hll_py_32 = HyperLogLog(precision=12, hash_size=32)
    hll_py_64 = HyperLogLog(precision=12, hash_size=64)
    hll_rust_32 = RustHyperLogLog(precision=12, hash_size=32)
    hll_rust_64 = RustHyperLogLog(precision=12, hash_size=64)
    
    assert hll_py_32.is_empty()
    assert hll_py_64.is_empty()
    assert hll_rust_32.is_empty()
    assert hll_rust_64.is_empty()
    
    # Test with single item
    test_item = "test_item"
    hll_py_32.add_string(test_item)
    hll_py_64.add_string(test_item)
    hll_rust_32.add_string(test_item)
    hll_rust_64.add_string(test_item)
    
    assert not hll_py_32.is_empty()
    assert not hll_py_64.is_empty()
    assert not hll_rust_32.is_empty()
    assert not hll_rust_64.is_empty()
    
    # Verify estimates are reasonable for single item
    assert 0.9 < hll_py_32.estimate_cardinality() < 1.1
    assert 0.9 < hll_py_64.estimate_cardinality() < 1.1
    assert 0.9 < hll_rust_32.estimate_cardinality() < 1.1
    assert 0.9 < hll_rust_64.estimate_cardinality() < 1.1

def test_precision_limits():
    """Test precision limits with both hash sizes."""
    # Test maximum precision for 32-bit hash
    with pytest.raises(ValueError, match="Precision must be less than hash_size"):
        HyperLogLog(precision=32, hash_size=32)
    
    with pytest.raises(ValueError, match="Precision must be less than hash_size"):
        RustHyperLogLog(precision=32, hash_size=32)
    
    # Test maximum precision for 64-bit hash
    with pytest.raises(ValueError, match="Precision must be less than hash_size"):
        HyperLogLog(precision=64, hash_size=64)
    
    with pytest.raises(ValueError, match="Precision must be less than hash_size"):
        RustHyperLogLog(precision=64, hash_size=64)
    
    # Test minimum precision
    with pytest.raises(ValueError, match="Precision must be at least 4"):
        HyperLogLog(precision=3, hash_size=32)
    
    with pytest.raises(ValueError, match="Precision must be at least 4"):
        RustHyperLogLog(precision=3, hash_size=32) 