#!/usr/bin/env python3
"""
Simple test script to verify the Rust Jaccard methods work correctly.
"""
import time
import random
import sys

# Try to import the rust_hll module
try:
    import rust_hll
    print("rust_hll module imported successfully")
except ImportError as e:
    print(f"Error importing rust_hll: {e}")
    print("Please build the Rust extension first:")
    print("  cd rust_hll")
    print("  python -m maturin develop")
    sys.exit(1)

def generate_test_sets(size1, size2, overlap_ratio):
    """Generate two sets with controlled overlap for testing."""
    overlap_size = int(min(size1, size2) * overlap_ratio)
    universe = list(range(max(size1, size2) * 3))  # Larger universe to avoid accidental overlaps
    
    # Shuffle the universe to ensure randomness
    random.shuffle(universe)
    
    # Create the first set
    set1 = set(universe[:size1])
    
    # Create the second set with controlled overlap
    overlap = list(set1)[:overlap_size]
    remaining = universe[size1:size1+(size2-overlap_size)]
    set2 = set(overlap + remaining)
    
    # Calculate true Jaccard similarity
    true_jaccard = len(set1.intersection(set2)) / len(set1.union(set2))
    
    return set1, set2, true_jaccard

def test_jaccard_methods():
    """Test all Jaccard similarity methods."""
    print("\n=== Testing Jaccard Similarity Methods ===")
    
    # Test parameters
    precisions = [12, 16]
    set_sizes = [(1000, 1000), (10000, 5000)]
    overlap_ratios = [0.1, 0.5, 0.8]
    
    for precision in precisions:
        print(f"\nPrecision: {precision}")
        
        for size1, size2 in set_sizes:
            print(f"\nSet sizes: {size1}, {size2}")
            
            for overlap_ratio in overlap_ratios:
                # Generate test sets
                set1, set2, true_jaccard = generate_test_sets(size1, size2, overlap_ratio)
                print(f"\nOverlap ratio: {overlap_ratio:.1f}")
                print(f"True Jaccard: {true_jaccard:.4f}")
                
                # Create sketches
                sketch1 = rust_hll.RustHLL(precision)
                sketch2 = rust_hll.RustHLL(precision)
                
                # Add elements to sketches
                for item in set1:
                    sketch1.add_value(str(item))
                for item in set2:
                    sketch2.add_value(str(item))
                
                # Get cardinality estimates
                card1 = sketch1.estimate_cardinality()
                card2 = sketch2.estimate_cardinality()
                print(f"Cardinality estimates: {card1:.1f}, {card2:.1f}")
                
                # Test original jaccard method
                try:
                    original_jaccard = sketch1.jaccard(sketch2)
                    error = abs(original_jaccard - true_jaccard)
                    rel_error = error / true_jaccard if true_jaccard > 0 else float('inf')
                    print(f"Original jaccard:       {original_jaccard:.4f} (error: {error:.4f}, rel_error: {rel_error:.2%})")
                except Exception as e:
                    print(f"Original jaccard failed: {e}")
                
                # Test register match method
                try:
                    register_jaccard = sketch1.jaccard_register_match(sketch2)
                    error = abs(register_jaccard - true_jaccard)
                    rel_error = error / true_jaccard if true_jaccard > 0 else float('inf')
                    print(f"Register match jaccard: {register_jaccard:.4f} (error: {error:.4f}, rel_error: {rel_error:.2%})")
                except Exception as e:
                    print(f"Register match jaccard failed: {e}")
                
                # Test minmax method
                try:
                    minmax_jaccard = sketch1.jaccard_minmax(sketch2)
                    error = abs(minmax_jaccard - true_jaccard)
                    rel_error = error / true_jaccard if true_jaccard > 0 else float('inf')
                    print(f"Minmax jaccard:        {minmax_jaccard:.4f} (error: {error:.4f}, rel_error: {rel_error:.2%})")
                except Exception as e:
                    print(f"Minmax jaccard failed: {e}")
                
                # Test minhash method
                try:
                    minhash_jaccard = sketch1.jaccard_with_method(sketch2, "minhash")
                    error = abs(minhash_jaccard - true_jaccard)
                    rel_error = error / true_jaccard if true_jaccard > 0 else float('inf')
                    print(f"Minhash jaccard:       {minhash_jaccard:.4f} (error: {error:.4f}, rel_error: {rel_error:.2%})")
                except Exception as e:
                    print(f"Minhash jaccard failed: {e}")

if __name__ == "__main__":
    random.seed(42)  # For reproducibility
    test_jaccard_methods() 