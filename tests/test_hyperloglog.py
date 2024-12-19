import pytest
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.abstractsketch import AbstractSketch
import csv
from datetime import datetime
import os

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
        """Test Jaccard similarity estimation."""
        sketch1 = HyperLogLog(precision=8)
        sketch2 = HyperLogLog(precision=8)
        
        # Add same strings
        sketch1.add_string("test")
        sketch2.add_string("test")
        
        # Should be similar
        assert sketch1.estimate_jaccard(sketch2) > 0.9
        
        # Add different strings
        sketch2.add_string("different")
        
        # Should be less similar
        assert sketch1.estimate_jaccard(sketch2) < 0.9

@pytest.mark.full
class TestHyperLogLogFull:
    """Full tests for HyperLogLog class."""
    
    def test_cardinality_accuracy(self):
        """Test cardinality estimation accuracy."""
        sketch = HyperLogLog(precision=14)  # Higher precision for accuracy test
        
        # Add known number of unique items
        n_items = 10000
        for i in range(n_items):
            sketch.add_string(f"item{i}")
            
        estimate = sketch.estimate_cardinality()
        error = abs(estimate - n_items) / n_items
        
        # Should be within 2% error for this precision
        assert error < 0.02
        
    def test_different_seeds(self):
        """Test that different seeds give different results."""
        sketch1 = HyperLogLog(precision=8, seed=1)
        sketch2 = HyperLogLog(precision=8, seed=2)
        
        # Add same strings
        for i in range(1000):
            s = f"item{i}"
            sketch1.add_string(s)
            sketch2.add_string(s)
            
        # Should give different register values
        assert not (sketch1.registers == sketch2.registers).all() 
    
    def test_hyperloglog_accuracy(self):
        """Test HyperLogLog accuracy with various set sizes and overlaps."""
        results = []
        
        for precision in [8, 10, 12, 14]:
            results.extend([
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
                )
            ])
        
        # Verify results
        for result in results:
            assert result['absolute_error'] < 0.1, f"Error too large for {result['test_name']}"

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