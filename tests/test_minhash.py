import pytest # type: ignore
from hammock.lib.minhash import MinHash
import csv
from datetime import datetime
import os

# Add the slow marker at the top of the file
pytest.mark.slow = pytest.mark.skipif(
    not pytest.config.getoption("--runslow"),
    reason="need --runslow option to run"
)

def run_test_case(num_hashes: int, name: str, desc: str, expected: float, 
                  set1_size: int, set2_size: int, set2_offset: int = 0):
    """Run a single test case with given parameters and return results."""
    print(f"\nTest case with num_hashes={num_hashes}: {name}")
    print(desc)
    print(f"Expected Jaccard: {expected:.3f}")
    
    sketch1 = MinHash(num_hashes=num_hashes)
    sketch2 = MinHash(num_hashes=num_hashes)
    
    # Add more elements for better statistical accuracy
    for i in range(set1_size * 10):  # Multiply by 10 for more elements
        sketch1.add_string(str(i))
    for i in range(set2_size * 10):  # Multiply by 10 for more elements
        sketch2.add_string(str(i + set2_offset))
        
    jaccard = sketch1.estimate_jaccard(sketch2)
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
        sketch = MinHash(num_hashes=16)
        assert sketch.num_hashes == 16
    
    def test_add_string(self):
        """Test adding strings to MinHash."""
        sketch = MinHash(num_hashes=16)
        sketch.add_string("test")
        assert len(sketch.signatures) == 16
    
    def test_jaccard(self):
        """Test MinHash estimate_jaccard method."""
        sketch1 = MinHash(num_hashes=16)
        sketch2 = MinHash(num_hashes=16)
        sketch1.add_string("test")
        sketch2.add_string("test")
        assert sketch1.estimate_jaccard(sketch2) == 1.0

@pytest.mark.full
class TestMinHashFull:
    """Full tests for MinHash class."""
    
    def test_minhash_accuracy(self):
        """Test MinHash accuracy with various set sizes and overlaps."""
        results = []
        
        for num_hashes in [128, 256, 512, 1024]:  # Increased number of hashes
            results.extend([
                run_test_case(
                    num_hashes=num_hashes,
                    name="Perfect overlap - small",
                    desc="Small sets with perfect overlap",
                    expected=1.0,
                    set1_size=100,
                    set2_size=100
                ),
                run_test_case(
                    num_hashes=num_hashes,
                    name="No overlap - small",
                    desc="Small sets with no overlap",
                    expected=0.0,
                    set1_size=100,
                    set2_size=100,
                    set2_offset=1000  # Increased offset
                ),
                run_test_case(
                    num_hashes=num_hashes,
                    name="Partial overlap - small",
                    desc="Small sets with 50% overlap",
                    expected=0.5,
                    set1_size=1000,  # Increased size
                    set2_size=1000,  # Increased size
                    set2_offset=500
                )
            ])
        
        # Verify results
        for result in results:
            assert result['absolute_error'] < 0.1, f"Error too large for {result['test_name']}"
    
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

@pytest.mark.slow
def test_large_set_accuracy():
    """Test MinHash accuracy with very large sets."""
    # This test takes a long time to run
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
    # ... rest of test ...

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

if __name__ == "__main__":
    # Run quick tests
    test = TestMinHashQuick()
    test.test_init()
    test.test_add_string()
    test.test_jaccard()
    
    # Run accuracy tests
    test = TestMinHashFull()
    test.test_minhash_accuracy()