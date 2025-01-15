#!/usr/bin/env python

from hammock.lib.sequences import SequenceSketch
from hammock.lib.intervals import IntervalSketch

def test_sequences():
    """Test sequence sketching with different methods."""
    # Create two similar sequences
    seq1 = "ATCGATCGATCG" * 10
    seq2 = "ATCGATCGATCG" * 8 + "GCATGCATGCAT" * 2
    
    # Test with different sketch types
    for sketch_type in ["minimizer", "hyperloglog", "minhash"]:
        print(f"\nTesting {sketch_type} sketch for sequences:")
        
        # Create sketches
        sketch1 = SequenceSketch(sketch_type=sketch_type, kmer_size=8)
        sketch2 = SequenceSketch(sketch_type=sketch_type, kmer_size=8)
        
        # Add sequences
        sketch1.add_sequence(seq1)
        sketch2.add_sequence(seq2)
        
        # Compare
        jaccard = sketch1.sketch.estimate_jaccard(sketch2.sketch)
        print(f"Jaccard similarity: {jaccard:.3f}")
        print(f"Cardinality 1: {sketch1.sketch.estimate_cardinality():.0f}")
        print(f"Cardinality 2: {sketch2.sketch.estimate_cardinality():.0f}")

def test_intervals():
    """Test interval sketching with different methods."""
    # Create some test BED lines
    bed1 = ["chr1\t100\t200", "chr1\t300\t400", "chr2\t100\t200"]
    bed2 = ["chr1\t100\t200", "chr1\t300\t400", "chr3\t100\t200"]
    
    # Test with different sketch types
    for sketch_type in ["hyperloglog", "minhash", "exact"]:
        print(f"\nTesting {sketch_type} sketch for intervals:")
        
        # Create sketches
        sketch1 = IntervalSketch(mode="A", sketch_type=sketch_type)
        sketch2 = IntervalSketch(mode="A", sketch_type=sketch_type)
        
        # Add intervals
        for line in bed1:
            interval, points, size = sketch1.bedline(line, mode="A", sep="-")
            if interval:
                sketch1.sketch.add_string(interval)
                
        for line in bed2:
            interval, points, size = sketch2.bedline(line, mode="A", sep="-")
            if interval:
                sketch2.sketch.add_string(interval)
        
        # Compare
        jaccard = sketch1.sketch.estimate_jaccard(sketch2.sketch)
        print(f"Jaccard similarity: {jaccard:.3f}")
        print(f"Cardinality 1: {sketch1.sketch.estimate_cardinality():.0f}")
        print(f"Cardinality 2: {sketch2.sketch.estimate_cardinality():.0f}")
        
        # Verify expected Jaccard similarity (2 shared intervals out of 4 total)
        assert abs(jaccard - 0.5) < 0.1, f"Expected Jaccard ~0.5, got {jaccard:.3f}"
        
        # Verify cardinalities
        assert abs(sketch1.sketch.estimate_cardinality() - 3) < 0.5
        assert abs(sketch2.sketch.estimate_cardinality() - 3) < 0.5

if __name__ == "__main__":
    print("Testing sequence sketching...")
    test_sequences()
    
    print("\nTesting interval sketching...")
    test_intervals() 