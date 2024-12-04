from hammock.lib.minhash import MinHash
import csv
from datetime import datetime
import os

def run_test_case(num_hashes: int, name: str, desc: str, expected: float, 
                  set1_size: int, set2_size: int, set2_offset: int = 0):
    """Run a single test case with given parameters and return results."""
    print(f"\nTest case with num_hashes={num_hashes}: {name}")
    print(desc)
    print(f"Expected Jaccard: {expected:.3f}")
    
    mh1 = MinHash(num_hashes=num_hashes)
    mh2 = MinHash(num_hashes=num_hashes)
    
    for i in range(set1_size):
        mh1.add_string(str(i))
    for i in range(set2_size):
        mh2.add_string(str(i + set2_offset))
        
    jaccard = mh1.estimate_jaccard(mh2)
    error = abs(jaccard - expected)
    
    print(f"Calculated Jaccard similarity: {jaccard:.3f}")
    print(f"Difference from expected: {error:.3f}")
    
    return {
        'num_hashes': num_hashes,
        'test_name': name,
        'set1_size': set1_size,
        'set2_size': set2_size,
        'expected_jaccard': expected,
        'calculated_jaccard': jaccard,
        'absolute_error': error
    }

def test_minhash_estimates():
    results = []
    
    # Test different numbers of hash functions
    for num_hashes in [16, 32, 64, 128, 256, 512, 1024]:
        print(f"\n{'='*60}")
        print(f"Testing with num_hashes={num_hashes}")
        print('='*60)
        
        # Test 1: Very sparse with perfect overlap
        results.append(run_test_case(
            num_hashes=num_hashes,
            name="Very sparse with perfect overlap",
            desc="Sketch 1: 10 integers (0-9)\n"
                 "Sketch 2: 10 integers (0-9)\n"
                 "Perfect overlap",
            expected=1.0,
            set1_size=10,
            set2_size=10
        ))
        
        # Test 2: Sparse vs Dense
        results.append(run_test_case(
            num_hashes=num_hashes,
            name="Sparse vs Dense",
            desc="Sketch 1: 100 integers (0-99)\n"
                 "Sketch 2: 1000 integers (0-999)\n"
                 "100 integers overlap",
            expected=0.1,  # 100/1000
            set1_size=100,
            set2_size=1000
        ))
        
        # Test 3: Medium density with overlap
        results.append(run_test_case(
            num_hashes=num_hashes,
            name="Medium density with overlap",
            desc="Sketch 1: 1000 integers (0-999)\n"
                 "Sketch 2: 1000 integers (100-1099)\n"
                 "900 integers overlap",
            expected=0.82,  # 900/1100
            set1_size=1000,
            set2_size=1000,
            set2_offset=100
        ))
        
        # Test 4: Dense with high overlap
        results.append(run_test_case(
            num_hashes=num_hashes,
            name="Dense with high overlap",
            desc="Sketch 1: 10000 integers (0-9999)\n"
                 "Sketch 2: 10000 integers (1000-10999)\n"
                 "9000 integers overlap",
            expected=0.82,  # 9000/11000
            set1_size=10000,
            set2_size=10000,
            set2_offset=1000
        ))
        
        # Test 5: Very dense with partial overlap
        results.append(run_test_case(
            num_hashes=num_hashes,
            name="Very dense with partial overlap",
            desc="Sketch 1: 100000 integers (0-99999)\n"
                 "Sketch 2: 100000 integers (50000-149999)\n"
                 "50000 integers overlap",
            expected=0.33,  # 50000/150000
            set1_size=100000,
            set2_size=100000,
            set2_offset=50000
        ))
        
        # Test 6: Extremely dense with minimal overlap
        results.append(run_test_case(
            num_hashes=num_hashes,
            name="Extremely dense with minimal overlap",
            desc="Sketch 1: 1000000 integers (0-999999)\n"
                 "Sketch 2: 1000000 integers (900000-1899999)\n"
                 "100000 integers overlap",
            expected=0.053,  # 100000/1900000
            set1_size=1000000,
            set2_size=1000000,
            set2_offset=900000
        ))
        
        # Test 7: Drastically different sparsity
        results.append(run_test_case(
            num_hashes=num_hashes,
            name="Drastically different sparsity",
            desc="Sketch 1: 10 integers (0-9)\n"
                 "Sketch 2: 10000 integers (0-9999)\n"
                 "10 integers overlap",
            expected=0.001,  # 10/10000
            set1_size=10,
            set2_size=10000
        ))
        
        # Test 8: Sparse vs Very Dense
        results.append(run_test_case(
            num_hashes=num_hashes,
            name="Sparse vs Very Dense",
            desc="Sketch 1: 100 integers (0-99)\n"
                 "Sketch 2: 100000 integers (0-99999)\n"
                 "100 integers overlap",
            expected=0.001,  # 100/100000
            set1_size=100,
            set2_size=100000
        ))
    
    # Create results directory if it doesn't exist
    os.makedirs('test_results', exist_ok=True)
    
    # Write results to CSV
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    filename = f'test_results/minhash_numhashes_test_{timestamp}.csv'
    
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    
    print(f"\nResults written to {filename}")

if __name__ == "__main__":
    test_minhash_estimates() 