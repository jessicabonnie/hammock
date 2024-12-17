from hammock.lib.sketchclass import Sketch
import csv
from datetime import datetime
import os
import pytest

def run_test_case(precision: int, name: str, desc: str, expected: float, 
                  set1_size: int, set2_size: int, set2_offset: int = 0):
    """Run a single test case with given parameters and return results."""
    print(f"\nTest case with precision={precision}: {name}")
    print(desc)
    print(f"Expected Jaccard: {expected:.3f}")
    
    sketch1 = Sketch(precision=precision, sketch_type="hyperloglog")
    sketch2 = Sketch(precision=precision, sketch_type="hyperloglog")
    
    for i in range(set1_size):
        sketch1.add_int(i)
    for i in range(set2_size):
        sketch2.add_int(i + set2_offset)
        
    jaccard = sketch1.estimate_jaccard(sketch2)
    error = abs(jaccard - expected)
    
    print(f"Calculated Jaccard similarity: {jaccard:.3f}")
    print(f"Difference from expected: {error:.3f}")
    
    return {
        'precision': precision,
        'num_registers': 2**precision,
        'test_name': name,
        'set1_size': set1_size,
        'set2_size': set2_size,
        'expected_jaccard': expected,
        'calculated_jaccard': jaccard,
        'absolute_error': error
    }

@pytest.mark.quick
def test_hll_quick():
    """Quick tests with small sets and lower precision"""
    results = []
    
    for precision in [4, 6, 8]:
        results.extend([
            # Basic functionality tests
            run_test_case(
                precision=precision,
                name="Perfect overlap - small",
                desc="Small sets with perfect overlap",
                expected=1.0,
                set1_size=10,
                set2_size=10
            ),
            run_test_case(
                precision=precision,
                name="No overlap - small",
                desc="Small sets with no overlap",
                expected=0.0,
                set1_size=10,
                set2_size=10,
                set2_offset=10
            ),
            run_test_case(
                precision=precision,
                name="Partial overlap - small",
                desc="Small sets with 50% overlap",
                expected=0.5,
                set1_size=100,
                set2_size=100,
                set2_offset=50
            ),
            # Edge cases
            run_test_case(
                precision=precision,
                name="Empty sets",
                desc="Testing empty set handling",
                expected=1.0,
                set1_size=0,
                set2_size=0
            ),
            run_test_case(
                precision=precision,
                name="Single element - same",
                desc="Single element equality",
                expected=1.0,
                set1_size=1,
                set2_size=1
            )
        ])
    
    save_results(results, "hll_quick_test")

@pytest.mark.full
def test_hll_full():
    """Full test suite with larger sets and higher precision"""
    results = []
    
    for precision in [4, 6, 8, 10, 12, 14]:
        results.extend([
            # Large set tests
            run_test_case(
                precision=precision,
                name="Large sets - high overlap",
                desc="10K elements with 90% overlap",
                expected=0.90,
                set1_size=10000,
                set2_size=10000,
                set2_offset=1000
            ),
            run_test_case(
                precision=precision,
                name="Very large sets - small overlap",
                desc="100K elements with 0.1% overlap",
                expected=0.001,
                set1_size=100000,
                set2_size=100000,
                set2_offset=99900
            ),
            # Different size sets
            run_test_case(
                precision=precision,
                name="Different sizes - medium",
                desc="1K vs 10K elements",
                expected=0.1,
                set1_size=1000,
                set2_size=10000
            ),
            # Precision-sensitive cases
            run_test_case(
                precision=precision,
                name="Dense set test",
                desc="Testing precision impact on dense sets",
                expected=0.82,
                set1_size=10000,
                set2_size=10000,
                set2_offset=1000
            ),
            run_test_case(
                precision=precision,
                name="Sparse set test",
                desc="Testing precision impact on sparse sets",
                expected=0.001,
                set1_size=100,
                set2_size=10000
            )
        ])
    
    save_results(results, "hll_full_test")

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
    test_hll_quick() 