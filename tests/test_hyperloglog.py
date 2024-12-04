from hammock.lib.hyperloglog import HyperLogLog

def test_hll_estimates():
    # Create HLLs with precision 8
    hll1 = HyperLogLog(precision=8)
    hll2 = HyperLogLog(precision=8)
    
    # Test case 1: Similar number of items
    print("\nTest 1: Similar counts (1000 vs 1000)")
    for i in range(1000):
        hll1.add_int(i)
        hll2.add_int(i + 100)  # 900 items overlap
    print(f"Jaccard similarity: {hll1.estimate_jaccard(hll2)}")
    
    # Test case 2: Sparse vs Dense
    hll3 = HyperLogLog(precision=8)
    hll4 = HyperLogLog(precision=8)
    print("\nTest 2: Sparse vs Dense (100 vs 1000)")
    for i in range(100):
        hll3.add_int(i)
    for i in range(1000):
        hll4.add_int(i)  # Complete overlap with hll3
    print(f"Jaccard similarity: {hll3.estimate_jaccard(hll4)}")
    
    # Test case 3: Very sparse
    hll5 = HyperLogLog(precision=8)
    hll6 = HyperLogLog(precision=8)
    print("\nTest 3: Very sparse (10 vs 10)")
    for i in range(10):
        hll5.add_int(i)
        hll6.add_int(i)  # Complete overlap
    print(f"Jaccard similarity: {hll5.estimate_jaccard(hll6)}")

if __name__ == "__main__":
    test_hll_estimates() 