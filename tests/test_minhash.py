import pytest # type: ignore
from hammock.lib.minhash import MinHash
import csv
from datetime import datetime
import os

timelimit = 480

def run_test_case(num_hashes: int, name: str, desc: str, expected: float, 
                  set1_size: int, set2_size: int, set2_offset: int = 0):
    """Run a single test case with given parameters and return results."""
    print(f"\nTest case with num_hashes={num_hashes}: {name}")
    print(desc)
    print(f"Expected Jaccard: {expected:.3f}")
    
    sketch1 = MinHash(num_hashes=num_hashes, kmer_size=3, debug=True)
    sketch2 = MinHash(num_hashes=num_hashes, kmer_size=3, debug=True)
    
    # Add more elements for better statistical accuracy
    for i in range(set1_size * 10):  # Multiply by 10 for more elements
        sketch1.add_string(str(i).zfill(3))
    for i in range(set2_size * 10):  # Multiply by 10 for more elements
        sketch2.add_string(str(i + set2_offset).zfill(3))
        
    jaccard = sketch1.similarity_values(sketch2)['jaccard_similarity']
    error = abs(jaccard - expected)
    
    print(f"Calculated Jaccard similarity: {jaccard:.3f}")
    print(f"Difference from expected: {error:.3f}")
    
    return {
        'num_hashes': num_hashes,
        'test_name': name,
        'set1_size': set1_size * 10,  # Adjusted size
        'set2_size': set2_size * 10,  # Adjusted size
        'expected_jaccard': expected,
        'calculated_jaccard': jaccard,
        'absolute_error': error
    }

@pytest.mark.quick
class TestMinHashQuick:
    """Quick tests for MinHash class."""
    
    def test_init(self):
        """Test MinHash initialization."""
        sketch = MinHash(num_hashes=16, debug=False)
        assert sketch.num_hashes == 16
    
    def test_add_string(self):
        """Test adding strings to MinHash."""
        sketch = MinHash(num_hashes=16, debug=False)
        sketch.add_string("test")
        assert len(sketch.min_hashes) == 16
    
    @pytest.mark.quick
    def test_jaccard(self):
        """Test Jaccard similarity calculation."""
        # Test identical sketches
        sketch1 = MinHash(num_hashes=128, debug=False)
        sketch2 = MinHash(num_hashes=128, debug=False)
        for i in range(100):
            sketch1.add_int(i)
            sketch2.add_int(i)
        
        result = sketch1.similarity_values(sketch2)
        assert result['jaccard_similarity'] == 1.0
        
        # Test partial overlap with fresh sketches
        sketch1 = MinHash(num_hashes=128, debug=False)
        sketch2 = MinHash(num_hashes=128, debug=False)
        
        # Add elements to first sketch: {0, 1, ..., 99}
        for i in range(100):
            sketch1.add_int(i)
            
        # Add elements to second sketch: {50, 51, ..., 149}
        for i in range(50, 150):
            sketch2.add_int(i)
        
        # Expected Jaccard: |intersection| / |union| = 50 / 150 = 0.333...
        result = sketch1.similarity_values(sketch2)
        assert 0.25 <= result['jaccard_similarity'] <= 0.45  # Allow for estimation error

@pytest.mark.full
@pytest.mark.timeout(timelimit)
class TestMinHashFull:
    """Full tests for MinHash class."""
    
    @pytest.mark.full
    @pytest.mark.timeout(timelimit)
    def test_minhash_accuracy(self):
        """Test MinHash accuracy with various set sizes and overlaps."""
        results = []
        for num_hashes in [128, 256, 512]:
            results.extend([
                run_test_case(
                    num_hashes=num_hashes,
                    name="No overlap - small",
                    desc="Small sets (1k elements) with no overlap",
                    expected=0.0,
                    set1_size=1000,
                    set2_size=1000,
                    set2_offset=1000
                ),
                run_test_case(
                    num_hashes=num_hashes,
                    name="Complete overlap - small",
                    desc="Small sets (1k elements) with complete overlap",
                    expected=1.0,
                    set1_size=1000,
                    set2_size=1000,
                    set2_offset=0
                ),
                run_test_case(
                    num_hashes=num_hashes,
                    name="Partial overlap - small",
                    desc="Small sets (1k elements) with 50% overlap",
                    expected=0.5,
                    set1_size=1000,
                    set2_size=1000,
                    set2_offset=500
                )
            ])
            
        # Verify results
        for result in results:
            assert result['absolute_error'] < 0.1, \
                f"Error too large for {result['test_name']}"
        
        save_results(results, "minhash_accuracy_test")
    
    @pytest.mark.full
    @pytest.mark.timeout(timelimit)
    def test_minhash_full(self):
        """Full test suite with larger sets and more hash functions"""
        results = []
        
        for num_hashes in [16, 32, 64, 128, 256, 512, 1024]:
            results.extend([
                # Large set tests
                run_test_case(
                    num_hashes=num_hashes,
                    name="Large sets - high overlap",
                    desc="10K elements with 90% overlap",
                    expected=0.90,
                    set1_size=10000,
                    set2_size=10000,
                    set2_offset=1000
                ),
                run_test_case(
                    num_hashes=num_hashes,
                    name="Very large sets - small overlap",
                    desc="100K elements with 0.1% overlap",
                    expected=0.001,
                    set1_size=100000,
                    set2_size=100000,
                    set2_offset=99900
                ),
                # Different size sets
                run_test_case(
                    num_hashes=num_hashes,
                    name="Different sizes - medium",
                    desc="1K vs 10K elements",
                    expected=0.1,
                    set1_size=1000,
                    set2_size=10000
                ),
                # Edge cases with large sets
                run_test_case(
                    num_hashes=num_hashes,
                    name="Large identical sets",
                    desc="50K identical elements",
                    expected=1.0,
                    set1_size=50000,
                    set2_size=50000
                )
            ])
        
        save_results(results, "minhash_full_test")

@pytest.mark.full
@pytest.mark.timeout(timelimit)
def test_large_set_accuracy():
    """Test MinHash accuracy with very large sets."""
    results = []
    for num_hashes in [128, 256, 512, 1024]:
        results.extend([
            run_test_case(
                num_hashes=num_hashes,
                name="Large set comparison",
                desc="Large sets (100k elements) with 50% overlap",
                expected=0.5,
                set1_size=100000,
                set2_size=100000,
                set2_offset=50000
            )
        ])
        
    # Verify results
    for result in results:
        assert result['absolute_error'] < 0.1, f"Error too large for {result['test_name']}"
    
    save_results(results, "minhash_large_test")

def save_results(results: list, test_name: str):
    """Save test results to CSV file"""
    os.makedirs('test_results', exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    filename = f'test_results/{test_name}_{timestamp}.csv'
    
    with open(filename, 'w', newline='') as f:
        if results:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
    
    print(f"\nResults written to {filename}")

@pytest.mark.quick
def test_minhash_empty():
    """Test MinHash with empty sketches."""
    sketch1 = MinHash(num_hashes=128, debug=True)
    sketch2 = MinHash(num_hashes=128, debug=True)
    
    # Add data to only one sketch
    sketch1.update("test")
    
    # Jaccard similarity between non-empty and empty should be 0
    assert sketch1.jaccard(sketch2) == 0.0
    
    # Empty sketches should have 0 similarity
    sketch3 = MinHash(num_hashes=128, debug=True)
    sketch4 = MinHash(num_hashes=128, debug=True)
    assert sketch3.jaccard(sketch4) == 0.0

def test_minhash_identical():
    mh1 = MinHash(kmer_size=3, debug=True)
    mh2 = MinHash(kmer_size=3, debug=True)
    
    mh1.add_string("test string")
    mh2.add_string("test string")
    
    sim = mh1.similarity_values(mh2)
    assert sim['jaccard_similarity'] == 1.0

def test_minhash_different():
    mh1 = MinHash(kmer_size=3, debug=True)
    mh2 = MinHash(kmer_size=3, debug=True)
    
    mh1.add_string("test string 1")
    mh2.add_string("test string 2")
    
    sim = mh1.similarity_values(mh2)
    assert 0.0 <= sim['jaccard_similarity'] <= 1.0

def test_minhash_different_sizes():
    mh1 = MinHash(num_hashes=128, debug=True)
    mh2 = MinHash(num_hashes=64, debug=True)
    
    mh1.add_string("test")
    mh2.add_string("test")
    
    with pytest.raises(ValueError):
        mh1.similarity_values(mh2)

def test_minhash():
    """Test MinHash sketching."""
    # Create sketches with debug enabled and non-zero seed
    seed = 42  # Any non-zero value would work
    sketch1 = MinHash(num_hashes=128, seed=seed, debug=True)
    sketch2 = MinHash(num_hashes=128, seed=seed, debug=True)
    
    # Add elements to sketch1 (0-999)
    for i in range(1000):
        sketch1.add_string(str(i).zfill(3))
    
    # Add elements to sketch2 (500-1499)
    # This creates 50% overlap with sketch1
    for i in range(500, 1500):
        sketch2.add_string(str(i).zfill(3))
    
    # Calculate Jaccard similarity
    jaccard = sketch1.estimate_jaccard(sketch2)
    
    # Expected Jaccard is 500/1500 = 1/3
    expected = 1/3
    assert abs(jaccard - expected) < 0.1, f"Expected Jaccard ~{expected:.3f}, got {jaccard:.3f}"

if __name__ == "__main__":
    # Run quick tests
    test = TestMinHashQuick()
    test.test_init()
    test.test_add_string()
    test.test_jaccard()
    
    # Run accuracy tests
    test = TestMinHashFull()
    test.test_minhash_accuracy()