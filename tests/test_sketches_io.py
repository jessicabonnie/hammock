import os
import tempfile
import pytest # type: ignore
import numpy as np # type: ignore
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash
from hammock.lib.exact import ExactCounter
from hammock.lib.sequences import SequenceSketch

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        yield tmpdirname

@pytest.mark.quick
class TestSketchesIOQuick:
    """Quick tests for sketch I/O functionality."""
    
    def test_hyperloglog_io(self, temp_dir):
        """Test HyperLogLog read/write functionality."""
        # Create and populate a HLL sketch
        hll = HyperLogLog(precision=8)
        hll.add_string("ACGTACGT")
        
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
    
    def test_minhash_io(self, temp_dir):
        """Test MinHash read/write functionality."""
        # Create and populate a MinHash sketch
        mh = MinHash(num_hashes=128)
        mh.add_string("ACGTACGT")
        
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
    
    def test_exact_io(self, temp_dir):
        """Test ExactCounter read/write functionality."""
        # Create and populate an ExactCounter sketch
        ec = ExactCounter()
        ec.add_string("ACGT")
        
        # Write to file
        filepath = os.path.join(temp_dir, "test_exact.txt")
        ec.write(filepath)
        
        # Read back and compare
        ec2 = ExactCounter.read(filepath)
        
        assert ec.kmer_size == ec2.kmer_size
        assert ec.seed == ec2.seed
        assert ec.elements == ec2.elements

@pytest.mark.full
class TestSketchesIOFull:
    """Full tests for sketch I/O functionality."""
    
    def test_sequence_sketch_io(self, temp_dir):
        """Test SequenceSketch read/write functionality."""
        # Create and populate a sequence sketch
        seq = SequenceSketch(kmer_size=8, window_size=12)
        seq.add_string("ACGTACGTACGTACGT")
        
        # Write to file
        filepath = os.path.join(temp_dir, "test_seq.npz")
        seq.write_sketch(filepath)
        
        # Read back and compare
        seq2 = SequenceSketch.read_sketch(filepath)
        
        assert seq.kmer_size == seq2.kmer_size
        assert seq.window_size == seq2.window_size
        assert seq.seed == seq2.seed
        assert seq.estimate_cardinality() == seq2.estimate_cardinality()
        assert seq.estimate_jaccard(seq2) > 0.99
    
    def test_invalid_sketch_type(self, temp_dir):
        """Test that invalid sketch types raise appropriate errors."""
        filepath = os.path.join(temp_dir, "nonexistent.npz")
        with pytest.raises(ValueError):
            SequenceSketch.read_sketch(filepath)
    
    def test_file_not_found(self):
        """Test that attempting to read non-existent files raises appropriate errors."""
        with pytest.raises(FileNotFoundError):
            HyperLogLog.read_sketch("nonexistent.npz")
        with pytest.raises(FileNotFoundError):
            MinHash.read_sketch("nonexistent.npz")
        with pytest.raises(FileNotFoundError):
            ExactCounter.read_sketch("nonexistent.txt")
    
    def test_corrupted_files(self, temp_dir):
        """Test handling of corrupted files."""
        # Create corrupted npz file
        filepath = os.path.join(temp_dir, "corrupted.npz")
        with open(filepath, 'wb') as f:
            f.write(b'corrupted data')
        
        with pytest.raises((ValueError, OSError)):
            HyperLogLog.read_sketch(filepath)
        
        # Create corrupted text file
        filepath = os.path.join(temp_dir, "corrupted.txt")
        with open(filepath, 'w') as f:
            f.write('corrupted data\n')
        
        with pytest.raises((ValueError, KeyError)):
            ExactCounter.read_sketch(filepath)

if __name__ == "__main__":
    pytest.main([__file__])