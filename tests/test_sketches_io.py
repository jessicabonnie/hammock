from __future__ import annotations
import os
import tempfile
import pytest # type: ignore
import numpy as np # type: ignore
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.minhash import MinHash
from hammock.lib.exact import ExactCounter
from hammock.lib.sequences import SequenceSketch
import gzip
from hammock.lib.intervals import IntervalSketch

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
        hll.write(filepath)
        
        # Read back and compare
        hll2 = HyperLogLog.load(filepath)
        
        assert hll.precision == hll2.precision
        assert hll.kmer_size == hll2.kmer_size
        assert hll.window_size == hll2.window_size
        assert hll.seed == hll2.seed
        np.testing.assert_array_equal(hll.registers, hll2.registers)

    def test_hyperloglog_hash_size_io(self, temp_dir):
        """Test HyperLogLog read/write functionality with different hash sizes."""
        # Test with 32-bit hash
        hll_32 = HyperLogLog(precision=8, hash_size=32)
        hll_32.add_string("ACGTACGT")
        
        # Test with 64-bit hash
        hll_64 = HyperLogLog(precision=8, hash_size=64)
        hll_64.add_string("ACGTACGT")
        
        # Write to files
        filepath_32 = os.path.join(temp_dir, "test_hll_32.npz")
        filepath_64 = os.path.join(temp_dir, "test_hll_64.npz")
        hll_32.write(filepath_32)
        hll_64.write(filepath_64)
        
        # Read back and compare
        hll_32_loaded = HyperLogLog.load(filepath_32)
        hll_64_loaded = HyperLogLog.load(filepath_64)
        
        # Check hash sizes are preserved
        assert hll_32.hash_size == hll_32_loaded.hash_size == 32
        assert hll_64.hash_size == hll_64_loaded.hash_size == 64
        
        # Check other parameters are preserved
        assert hll_32.precision == hll_32_loaded.precision
        assert hll_64.precision == hll_64_loaded.precision
        assert hll_32.kmer_size == hll_32_loaded.kmer_size
        assert hll_64.kmer_size == hll_64_loaded.kmer_size
        assert hll_32.window_size == hll_32_loaded.window_size
        assert hll_64.window_size == hll_64_loaded.window_size
        assert hll_32.seed == hll_32_loaded.seed
        assert hll_64.seed == hll_64_loaded.seed
        
        # Check registers are preserved
        np.testing.assert_array_equal(hll_32.registers, hll_32_loaded.registers)
        np.testing.assert_array_equal(hll_64.registers, hll_64_loaded.registers)
        
        # Check cardinality estimates are reasonable
        assert abs(hll_32.estimate_cardinality() - hll_32_loaded.estimate_cardinality()) < 1e-10
        assert abs(hll_64.estimate_cardinality() - hll_64_loaded.estimate_cardinality()) < 1e-10
    
    def test_minhash_io(self, temp_dir):
        """Test MinHash read/write functionality."""
        # Create and populate a MinHash sketch
        mh = MinHash(num_hashes=128)
        mh.add_string("ACGTACGT")
        
        # Write to file
        filepath = os.path.join(temp_dir, "test_mh.npz")
        mh.write(filepath)
        
        # Read back and compare
        mh2 = MinHash.load(filepath)
        
        assert mh.num_hashes == mh2.num_hashes
        assert mh.kmer_size == mh2.kmer_size
        assert mh.window_size == mh2.window_size
        assert mh.seed == mh2.seed
        np.testing.assert_array_equal(mh.hashers, mh2.hashers)
    
    # def test_exact_io(self, temp_dir):
    #     """Test ExactCounter read/write functionality."""
    #     # Create and populate an ExactCounter sketch
    #     ec = ExactCounter()
    #     ec.add_string("ACGT")
        
    #     # Write to file
    #     filepath = os.path.join(temp_dir, "test_exact.txt")
    #     ec.write(filepath)
        
    #     # Read back and compare
    #     ec2 = ExactCounter.read(filepath)
        
    #     assert ec.kmer_size == ec2.kmer_size
    #     assert ec.seed == ec2.seed
    #     assert ec.elements == ec2.elements

@pytest.mark.full
class TestSketchesIOFull:
    """Full tests for sketch I/O functionality."""
    
    def test_sequence_sketch_io(self, temp_dir):
        """Test SequenceSketch read/write functionality."""
        # Create and populate a sequence sketch with explicit parameters
        window_size = 40  # Match the default window size
        seq = SequenceSketch(
            kmer_size=8,
            window_size=window_size,
            sketch_type="hyperloglog"
        )
        # Use a longer sequence that's at least a few windows long
        test_seq = "ACGTACGTACGTACGT" * 8  # 128 bp sequence, longer than window size
        seq.add_string(test_seq)
        
        # Write to file
        filepath = os.path.join(temp_dir, "test_seq.npz")
        seq.write(filepath)
        
        # Read back and compare
        seq2 = SequenceSketch.load(filepath)
        
        assert seq.kmer_size == seq2.kmer_size
        assert seq.window_size == seq2.window_size  # Now should be 40 for both
        assert seq.seed == seq2.seed
        assert seq.estimate_cardinality() == seq2.estimate_cardinality()
        assert seq.estimate_jaccard(seq2) > 0.99
    
    def test_invalid_sketch_type(self, temp_dir):
        """Test that invalid sketch types raise appropriate errors."""
        filepath = os.path.join(temp_dir, "nonexistent.npz")
        with pytest.raises(ValueError):
            SequenceSketch.load(filepath)
    
    def test_file_not_found(self):
        """Test that attempting to read non-existent files raises appropriate errors."""
        with pytest.raises(FileNotFoundError):
            HyperLogLog.load("nonexistent.npz")
        with pytest.raises(FileNotFoundError):
            MinHash.load("nonexistent.npz")
        # with pytest.raises(FileNotFoundError):
        #     ExactCounter.load("nonexistent.txt")
    
    def test_corrupted_files(self, temp_dir):
        """Test handling of corrupted files."""
        # Create corrupted npz file
        filepath = os.path.join(temp_dir, "corrupted.npz")
        with open(filepath, 'wb') as f:
            f.write(b'corrupted data')
        
        with pytest.raises((ValueError, OSError)):
            HyperLogLog.load(filepath)
    
    def test_exact_counter_io_disabled(self):
        """Test that ExactCounter I/O operations raise NotImplementedError."""
        counter = ExactCounter()
        with pytest.raises(NotImplementedError):
            counter.write("test.txt")
        with pytest.raises(NotImplementedError):
            ExactCounter.load("test.txt")

    def test_gff_file_processing(self, temp_dir):
        """Test processing of GFF format files."""
        filepath = os.path.join(temp_dir, "test.gff")
        with open(filepath, 'w') as f:
            f.write("""##gff-version 3
#!genome-build GRCh38.p13
chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tID=ENSG00000223972.5
chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tParent=ENSG00000223972.5
chr1\tHAVANA\tCDS\t12010\t12057\t.\t+\t0\tParent=ENSG00000223972.5
chr1\tHAVANA\tgene\t14404\t29570\t.\t-\t.\tID=ENSG00000227232.5
""")
        
        # Test basic GFF processing
        sketch = IntervalSketch.from_file(filepath, mode="A")
        assert sketch is not None
        assert sketch.num_intervals == 4
        
        # Test feature filtering
        sketch_genes = IntervalSketch.from_file(filepath, mode="A", feature_types=['gene'])
        assert sketch_genes is not None
        assert sketch_genes.num_intervals == 2

    def test_gff_gz_processing(self, temp_dir):
        """Test processing of gzipped GFF files."""
        filepath = os.path.join(temp_dir, "test.gff.gz")
        with gzip.open(filepath, 'wt') as f:
            f.write("""##gff-version 3
chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tID=ENSG00000223972.5
""")
        
        sketch = IntervalSketch.from_file(filepath, mode="A")
        assert sketch is not None
        assert sketch.num_intervals == 1

    def test_gff_coordinate_conversion(self, temp_dir):
        """Test conversion from 1-based GFF coordinates to 0-based internal coordinates."""
        gff_path = os.path.join(temp_dir, "test.gff")
        bed_path = os.path.join(temp_dir, "test.bed")
        
        # Write equivalent GFF and BED files
        with open(gff_path, 'w') as f:
            f.write("""##gff-version 3
chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=test1
""")
        
        with open(bed_path, 'w') as f:
            f.write("chr1\t0\t100\n")  # BED coordinates are 0-based
        
        sketch_gff = IntervalSketch.from_file(gff_path, mode="A")
        sketch_bed = IntervalSketch.from_file(bed_path, mode="A")
        
        assert sketch_gff is not None
        assert sketch_bed is not None
        assert sketch_gff.estimate_jaccard(sketch_bed) == 1.0

    def test_invalid_gff(self, temp_dir):
        """Test handling of invalid GFF files."""
        filepath = os.path.join(temp_dir, "invalid.gff")
        with open(filepath, 'w') as f:
            f.write("invalid\tgff\tdata\n")
        
        sketch = IntervalSketch.from_file(filepath, mode="A")
        assert sketch is not None
        assert sketch.num_intervals == 0

if __name__ == "__main__":
    pytest.main([__file__])