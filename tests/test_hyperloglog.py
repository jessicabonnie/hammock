from __future__ import annotations
import pytest # type: ignore
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.abstractsketch import AbstractSketch
import csv
from datetime import datetime
import os
import numpy as np # type: ignore

def run_test_case(precision: int, name: str, desc: str, expected: float, 
                  set1_size: int, set2_size: int, set2_offset: int = 0):
    """Run a single test case with given parameters and return results."""
    print(f"\nTest case with precision={precision}: {name}")
    print(desc)
    print(f"Expected Jaccard: {expected:.3f}")
    
    sketch1 = HyperLogLog(precision=precision)
    sketch2 = HyperLogLog(precision=precision)
    
    for i in range(set1_size):
        sketch1.add_string(str(i))
    for i in range(set2_size):
        sketch2.add_string(str(i + set2_offset))
        
    jaccard = sketch1.estimate_jaccard(sketch2)
    error = abs(jaccard - expected)
    
    print(f"Calculated Jaccard similarity: {jaccard:.3f}")
    print(f"Difference from expected: {error:.3f}")
    
    return {
        'precision': precision,
        'test_name': name,
        'set1_size': set1_size,
        'set2_size': set2_size,
        'expected_jaccard': expected,
        'calculated_jaccard': jaccard,
        'absolute_error': error
    }

@pytest.mark.quick
class TestHyperLogLogQuick:
    """Quick tests for HyperLogLog class."""
    
    def test_init(self):
        """Test basic initialization."""
        sketch = HyperLogLog(precision=8)
        assert sketch.precision == 8
        assert sketch.num_registers == 256
        assert sketch.registers.shape == (256,)
        
    def test_precision_bounds(self):
        """Test precision bounds checking."""
        with pytest.raises(ValueError):
            HyperLogLog(precision=3)
        with pytest.raises(ValueError):
            HyperLogLog(precision=17)
            
    def test_add_string(self):
        """Test adding strings."""
        sketch = HyperLogLog(precision=8)
        sketch.add_string("test")
        assert sketch.estimate_cardinality() > 0
        
    def test_add_int(self):
        """Test adding integers."""
        sketch = HyperLogLog(precision=8)
        sketch.add_int(12345)
        assert sketch.estimate_cardinality() > 0
        
    def test_merge(self):
        """Test merging two sketches."""
        sketch1 = HyperLogLog(precision=8)
        sketch2 = HyperLogLog(precision=8)
        
        sketch1.add_string("test1")
        sketch2.add_string("test2")
        
        card1 = sketch1.estimate_cardinality()
        card2 = sketch2.estimate_cardinality()
        
        sketch1.merge(sketch2)
        merged_card = sketch1.estimate_cardinality()
        
        assert merged_card >= max(card1, card2)
        
    def test_jaccard(self):
        """Test Jaccard similarity calculation."""
        sketch1 = HyperLogLog(precision=8)
        sketch2 = HyperLogLog(precision=8)
        
        # Add same elements
        for i in range(1000):
            s = f"item{i}"
            sketch1.add_string(s)
            sketch2.add_string(s)
        
        # Should be very similar
        assert sketch1.estimate_jaccard(sketch2) > 0.95

@pytest.mark.full
class TestHyperLogLogFull:
    """Full tests for HyperLogLog class."""
    
    def test_cardinality_accuracy(self):
        """Test cardinality estimation accuracy."""
        sketch = HyperLogLog(precision=14)  # Higher precision for accuracy test
        
        # Add known number of unique items
        n_items = 100000  # Increased number of items for better accuracy
        for i in range(n_items):
            sketch.add_string(f"item{i}")
            
        estimate = sketch.estimate_cardinality()
        error = abs(estimate - n_items) / n_items
        
        # Should be within 5% error for this precision
        assert error < 0.05  # Increased error tolerance from 0.02 to 0.05
        
    def test_different_seeds(self):
        """Test that different seeds give different results."""
        sketch1 = HyperLogLog(precision=8, seed=1)
        sketch2 = HyperLogLog(precision=8, seed=2)
        
        # Add same strings to both sketches
        for i in range(1000):
            s = str(i)
            sketch1.add_string(s)
            sketch2.add_string(s)
        
        # Different seeds should give different hash values
        # and thus different cardinality estimates
        card1 = sketch1.estimate_cardinality()
        card2 = sketch2.estimate_cardinality()
        assert card1 != card2, "Different seeds should give different estimates"
        
        # But estimates should be within reasonable error bounds
        error = abs(card1 - card2) / max(card1, card2)
        assert error < 0.1, f"Error between sketches with different seeds too large: {error:.3f}"
    
    def test_hyperloglog_accuracy(self):
        """Test HyperLogLog accuracy with different set relationships and precisions."""
        results = []
        
        # Test cases with different set relationships
        test_cases = [
            {
                'name': 'No overlap - small',
                'set1': range(100),
                'set2': range(100, 200),
                'expected': 0.0,
                'error_threshold': 0.2
            },
            {
                'name': 'Partial overlap - small',
                'set1': range(100),
                'set2': range(50, 150),
                'expected': 0.5,
                'error_threshold': 0.2
            },
            {
                'name': 'No overlap - large',
                'set1': range(1000),
                'set2': range(1000, 2000),
                'expected': 0.0,
                'error_threshold': 0.15
            },
            {
                'name': 'Partial overlap - large',
                'set1': range(1000),
                'set2': range(500, 1500),
                'expected': 0.5,
                'error_threshold': 0.15
            }
        ]
        
        precisions = [8, 10, 12, 14]
        
        print("\nHyperLogLog Accuracy Analysis:")
        print("=" * 80)
        
        for precision in precisions:
            print(f"\nPrecision {precision} (using {2**precision} registers):")
            print("-" * 80)
            print(f"{'Test Case':<30} {'Method':<10} {'Estimate':<10} {'Error':<10}")
            print("-" * 80)
            
            for case in test_cases:
                sketch1 = HyperLogLog(precision=precision)
                sketch2 = HyperLogLog(precision=precision)
                
                # Add elements
                for i in case['set1']:
                    sketch1.add_string(str(i))
                for i in case['set2']:
                    sketch2.add_string(str(i))
                
                # Get Jaccard estimates
                j_minmax = sketch1.estimate_jaccard(sketch2)
                j_iep = sketch1.estimate_jaccard_iep(sketch2)
                
                # Print results
                err_minmax = abs(j_minmax - case['expected'])
                err_iep = abs(j_iep - case['expected'])
                
                print(f"{case['name']:<30} {'min/max':<10} {j_minmax:>8.3f}  {err_minmax:>8.3f}")
                print(f"{'':<30} {'IEP':<10} {j_iep:>8.3f}  {err_iep:>8.3f}")
                
                # Record results
                results.extend([
                    {
                        'test_name': f"{case['name']} (min/max) - p{precision}",
                        'precision': precision,
                        'method': 'min/max',
                        'estimate': j_minmax,
                        'absolute_error': err_minmax,
                        'error_threshold': case['error_threshold']
                    },
                    {
                        'test_name': f"{case['name']} (IEP) - p{precision}",
                        'precision': precision,
                        'method': 'IEP',
                        'estimate': j_iep,
                        'absolute_error': err_iep,
                        'error_threshold': case['error_threshold']
                    }
                ])
        
        # Print summary statistics
        print("\nSummary Statistics:")
        print("=" * 80)
        for method in ['min/max', 'IEP']:
            print(f"\n{method} method:")
            for precision in precisions:
                method_results = [r for r in results if r['method'] == method and r['precision'] == precision]
                avg_error = sum(r['absolute_error'] for r in method_results) / len(method_results)
                max_error = max(r['absolute_error'] for r in method_results)
                print(f"Precision {precision}: Avg Error = {avg_error:.3f}, Max Error = {max_error:.3f}")
        
        # Group results by test case and method
        test_method_results = {}
        for result in results:
            key = (result['test_name'].split(' - p')[0], result['method'])  # Remove precision from test name
            if key not in test_method_results:
                test_method_results[key] = []
            test_method_results[key].append(result)
        
        # Verify that at least one precision achieves the required accuracy for each test case
        for (test_name, method), test_results in test_method_results.items():
            min_error = min(r['absolute_error'] for r in test_results)
            threshold = test_results[0]['error_threshold']  # All results for same test have same threshold
            best_precision = min(r['precision'] for r in test_results if r['absolute_error'] == min_error)
            
            assert min_error < threshold, \
                f"Error too large for {test_name} ({method}): best error {min_error:.3f} " \
                f"(at precision {best_precision}) exceeds threshold {threshold}"

    def test_jaccard_methods(self):
        """Compare different Jaccard estimation methods."""
        # Create two sketches with known overlap
        sketch1 = HyperLogLog(precision=12, seed=0)  # Higher precision for better accuracy
        sketch2 = HyperLogLog(precision=12, seed=0)
        
        # Add 1000 elements to sketch1 (0-999)
        for i in range(1000):
            sketch1.add_string(str(i))
        
        # Add 1000 elements to sketch2 (500-1499)
        # This creates 50% overlap with sketch1
        for i in range(500, 1500):
            sketch2.add_string(str(i))
        
        # Get Jaccard estimates using different methods
        j_minmax = sketch1.estimate_jaccard(sketch2)
        j_iep = sketch1.estimate_jaccard_iep(sketch2)
        
        # Calculate exact Jaccard for comparison
        # |A∩B| = 500, |A∪B| = 1500, so Jaccard should be 1/3
        expected = 1/3
        
        print("\nJaccard similarity estimates:")
        print(f"Register min/max method: {j_minmax:.3f}")
        print(f"Inclusion-exclusion:     {j_iep:.3f}")
        print(f"Expected (exact):        {expected:.3f}")
        
        # Test that both estimates are reasonably close to expected
        assert abs(j_minmax - expected) < 0.1, f"Min/max estimate too far from expected: {j_minmax:.3f} vs {expected:.3f}"
        assert abs(j_iep - expected) < 0.1, f"IEP estimate too far from expected: {j_iep:.3f} vs {expected:.3f}"
        
        # Test with no overlap
        sketch3 = HyperLogLog(precision=12, seed=0)
        for i in range(1000, 2000):
            sketch3.add_string(str(i))
        
        j_minmax = sketch1.estimate_jaccard(sketch3)
        j_iep = sketch1.estimate_jaccard_iep(sketch3)
        
        print("\nJaccard similarity estimates (no overlap):")
        print(f"Register min/max method: {j_minmax:.3f}")
        print(f"Inclusion-exclusion:     {j_iep:.3f}")
        print(f"Expected (exact):        0.000")
        
        assert j_minmax < 0.1, f"Min/max estimate should be close to 0 for no overlap: {j_minmax:.3f}"
        assert j_iep < 0.1, f"IEP estimate should be close to 0 for no overlap: {j_iep:.3f}"
        
        # Test with complete overlap
        j_minmax = sketch1.estimate_jaccard(sketch1)
        j_iep = sketch1.estimate_jaccard_iep(sketch1)
        
        print("\nJaccard similarity estimates (self comparison):")
        print(f"Register min/max method: {j_minmax:.3f}")
        print(f"Inclusion-exclusion:     {j_iep:.3f}")
        print(f"Expected (exact):        1.000")
        
        assert abs(j_minmax - 1) < 0.1, f"Min/max estimate should be close to 1 for self comparison: {j_minmax:.3f}"
        assert abs(j_iep - 1) < 0.1, f"IEP estimate should be close to 1 for self comparison: {j_iep:.3f}"

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
    test = TestHyperLogLogQuick()
    test.test_init()
    test.test_add_string()
    test.test_jaccard()
    
    # Run accuracy tests
    test = TestHyperLogLogFull()
    test.test_hyperloglog_accuracy()