from hammock.lib.hyperloglog import HyperLogLog
import csv
from datetime import datetime
import os

def run_test_case(precision: int, name: str, desc: str, expected: float, 
                  set1_size: int, set2_size: int, set2_offset: int = 0):
    """Run a single test case with given parameters and return results."""
    print(f"\nTest case with precision={precision}: {name}")
    print(desc)
    print(f"Expected Jaccard: {expected:.3f}")
    
    hll1 = HyperLogLog(precision=precision)
    hll2 = HyperLogLog(precision=precision)
    
    for i in range(set1_size):
        hll1.add_int(i)
    for i in range(set2_size):
        hll2.add_int(i + set2_offset)
        
    jaccard = hll1.estimate_jaccard(hll2)
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

def test_hll_estimates():
    results = []
    
    # Test different precision values
    for precision in [4, 6, 8, 10, 12, 14]:
        print(f"\n{'='*60}")
        print(f"Testing with precision {precision}")
        print(f"Number of registers: {2**precision}")
        print('='*60)
        
        # Test 1: Similar sized sets with high overlap
        results.append(run_test_case(
            precision=precision,
            name="Similar counts with overlap",
            desc="Sketch 1: 1000 integers (0-999)\n"
                 "Sketch 2: 1000 integers (100-1099)\n"
                 "900 integers overlap",
            expected=0.82,  # 900/1100
            set1_size=1000,
            set2_size=1000,
            set2_offset=100
        ))
        
        # Test 2: Sparse vs Dense
        results.append(run_test_case(
            precision=precision,
            name="Sparse vs Dense",
            desc="Sketch 1: 100 integers (0-99)\n"
                 "Sketch 2: 1000 integers (0-999)\n"
                 "100 integers overlap",
            expected=0.1,  # 100/1000
            set1_size=100,
            set2_size=1000
        ))
        
        # Test 3: Very sparse with perfect overlap
        results.append(run_test_case(
            precision=precision,
            name="Very sparse with perfect overlap",
            desc="Sketch 1: 10 integers (0-9)\n"
                 "Sketch 2: 10 integers (0-9)\n"
                 "Perfect overlap",
            expected=1.0,
            set1_size=10,
            set2_size=10
        ))
    
    # Create results directory if it doesn't exist
    os.makedirs('test_results', exist_ok=True)
    
    # Write results to CSV
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    filename = f'test_results/hll_precision_test_{timestamp}.csv'
    
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    
    print(f"\nResults written to {filename}")

if __name__ == "__main__":
    test_hll_estimates() 