import os
import tempfile
import pytest
import numpy as np
from hammock.lib.sketchclass import Sketch
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash
from hammock.lib.exact import ExactCounter

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        yield tmpdirname

def test_hyperloglog_io(temp_dir):
    """Test HyperLogLog read/write functionality."""
    # Create and populate a HLL sketch
    hll = HyperLogLog(precision=8, kmer_size=3, window_size=5, seed=42)
    hll.add_string("ACGTACGT")
    hll.add_string("TGCATGCA")
    
    # Write to file
    filepath = os.path.join(temp_dir, "test_hll.npz")
    hll.write_sketch(filepath)
    
    # Read back and compare
    hll2 = HyperLogLog.read_sketch(filepath)
    
    assert hll.precision == hll2.precision
    assert hll.kmer_size == hll2.kmer_size
    assert hll.window_size == hll2.window_size
    assert hll.seed == hll2.seed
    np.testing.assert_array_equal(hll.registers, hll2.registers)
    assert hll.estimate_cardinality() == hll2.estimate_cardinality()

def test_minhash_io(temp_dir):
    """Test MinHash read/write functionality."""
    # Create and populate a MinHash sketch
    mh = MinHash(num_hashes=128, kmer_size=3, window_size=5, seed=42)
    mh.add_string("ACGTACGT")
    mh.add_string("TGCATGCA")
    
    # Write to file
    filepath = os.path.join(temp_dir, "test_mh.npz")
    mh.write_sketch(filepath)
    
    # Read back and compare
    mh2 = MinHash.read_sketch(filepath)
    
    assert mh.num_hashes == mh2.num_hashes
    assert mh.kmer_size == mh2.kmer_size
    assert mh.window_size == mh2.window_size
    assert mh.seed == mh2.seed
    np.testing.assert_array_equal(mh.signatures, mh2.signatures)
    assert mh.estimate_cardinality() == mh2.estimate_cardinality()

def test_exact_io(temp_dir):
    """Test ExactCounter read/write functionality."""
    # Create and populate an ExactCounter sketch with minimal parameters
    ec = ExactCounter(kmer_size=0, window_size=0, seed=42)  # No k-mer processing
    ec.add_string("ACGT")  # Simple test string
    
    # Write to file
    filepath = os.path.join(temp_dir, "test_exact.txt")
    ec.write_sketch(filepath)
    
    # Read back and compare
    ec2 = ExactCounter.read_sketch(filepath)
    
    assert ec.kmer_size == ec2.kmer_size
    assert ec.window_size == ec2.window_size
    assert ec.seed == ec2.seed
    assert ec.elements == ec2.elements
    assert ec.estimate_cardinality() == ec2.estimate_cardinality()

def test_sketch_wrapper_io(temp_dir):
    """Test Sketch wrapper class read/write functionality for all sketch types."""
    test_data = [
        ("hyperloglog", "test_hll.npz"),
        ("minhash", "test_mh.npz"),
        ("exact", "test_exact.txt")
    ]
    
    for sketch_type, filename in test_data:
        # Create and populate a sketch with minimal parameters
        sketch = Sketch(
            sketch_type=sketch_type,
            precision=8,
            num_hashes=128,
            kmer_size=0,  # No k-mer processing
            window_size=0
        )
        sketch.add_string("ACGT")
        
        # Write to file
        filepath = os.path.join(temp_dir, filename)
        sketch.write_sketch(filepath)
        
        # Read back and compare
        sketch2 = Sketch.read_sketch(filepath, sketch_type=sketch_type)
        
        # Compare cardinality estimates with tolerance
        assert abs(sketch.estimate_cardinality() - sketch2.estimate_cardinality()) < 1e-10
        
        # Compare Jaccard similarity with self (should be 1.0)
        assert abs(sketch.estimate_jaccard(sketch2) - 1.0) < 1e-10

def test_invalid_sketch_type(temp_dir):
    """Test that invalid sketch types raise appropriate errors."""
    with pytest.raises(ValueError):
        Sketch.read_sketch("test.npz", sketch_type="invalid_type")

def test_file_not_found():
    """Test that attempting to read non-existent files raises appropriate errors."""
    with pytest.raises(FileNotFoundError):
        HyperLogLog.read_sketch("nonexistent.npz")
    with pytest.raises(FileNotFoundError):
        MinHash.read_sketch("nonexistent.npz")
    with pytest.raises(FileNotFoundError):
        ExactCounter.read_sketch("nonexistent.txt")

def test_corrupted_files(temp_dir):
    """Test handling of corrupted files."""
    # Create corrupted npz file
    filepath = os.path.join(temp_dir, "corrupted.npz")
    with open(filepath, 'wb') as f:
        f.write(b'corrupted data')
    
    with pytest.raises((ValueError, OSError)):  # NPZ files raise ValueError or OSError
        HyperLogLog.read_sketch(filepath)
    
    # Create corrupted text file
    filepath = os.path.join(temp_dir, "corrupted.txt")
    with open(filepath, 'w') as f:
        f.write('corrupted data\n')  # No metadata at all
    
    with pytest.raises((ValueError, KeyError)):  # Should raise when metadata is missing
        ExactCounter.read_sketch(filepath)