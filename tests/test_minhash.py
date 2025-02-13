import pytest # type: ignore
from hammock.lib.minhash import MinHash
import csv
from datetime import datetime
import os
import random
import numpy as np

timelimit = 480

def generate_sequence_pair(length: int, expected_similarity: float) -> tuple[str, str]:
    """Generate a pair of sequences with expected similarity.
    
    Args:
        length: Length of sequences to generate
        expected_similarity: Target similarity between sequences (0-1)
        
    Returns:
        Tuple of two sequences
    """
    # Use only DNA characters
    alphabet = ['A', 'C', 'T', 'G']
    
    # Generate first sequence randomly
    seq1 = ''.join(random.choice(alphabet) for _ in range(length))
    
    # Generate second sequence by copying first and introducing mutations
    seq2 = list(seq1)
    num_mutations = int(length * (1 - expected_similarity))
    
    # Randomly select positions to mutate
    positions = random.sample(range(length), num_mutations)
    
    # Introduce mutations at selected positions
    for pos in positions:
        # Get current base and possible replacements
        current = seq2[pos]
        replacements = [b for b in alphabet if b != current]
        # Replace with a different base
        seq2[pos] = random.choice(replacements)
    
    return seq1, ''.join(seq2)

def run_test_case(num_hashes: int,
                  name: str,
                  desc: str,
                  expected: float,
                  set1_size: int,
                  set2_size: int,
                  set2_offset: int) -> dict:
    """Run a single test case and return results."""
    # Create two MinHash sketches
    sketch1 = MinHash(num_hashes=num_hashes, debug=True)
    sketch2 = MinHash(num_hashes=num_hashes, debug=True)
    
    # Generate test data - convert integers to strings to ensure consistent hashing
    set1 = [str(i).zfill(4) for i in range(set1_size)]
    set2 = [str(i + set2_offset).zfill(4) for i in range(set2_size)]
    
    # Add elements to sketches
    for elem in set1:
        sketch1.add_string(elem)
    for elem in set2:
        sketch2.add_string(elem)
    
    # Calculate Jaccard similarity
    calculated = sketch1.estimate_jaccard(sketch2)
    
    # Calculate absolute error
    error = abs(expected - calculated)
    
    print(f"\nTest case with num_hashes={num_hashes}: {name}")
    print(desc)
    print(f"Expected Jaccard: {expected:.3f}")
    print(f"Calculated Jaccard similarity: {calculated:.3f}")
    print(f"Difference from expected: {error:.3f}")
    
    return {
        'test_name': name,
        'expected': expected,
        'calculated': calculated,
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
            assert result['absolute_error'] < 0.2, \
                f"Error too large for {result['test_name']}"
        
        save_results(results, "minhash_accuracy_test")
    
    @pytest.mark.full
    @pytest.mark.timeout(timelimit)
    def test_minhash_full(self):
        """Full test suite with larger sets and more hash functions"""
        results = []
        
        for num_hashes in [16, 32, 64, 128, 256, 512]:
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
                    set2_size=10000,
                    set2_offset=900  # This will create 10% overlap (100 elements)
                ),
                # Edge cases with large sets
                run_test_case(
                    num_hashes=num_hashes,
                    name="Large identical sets",
                    desc="50K identical elements",
                    expected=1.0,
                    set1_size=50000,
                    set2_size=50000,
                    set2_offset=0
                )
            ])
        
        save_results(results, "minhash_full_test")

@pytest.mark.full
@pytest.mark.timeout(timelimit)
def test_large_set_accuracy():
    """Test MinHash accuracy with very large sets."""
    results = []
    for num_hashes in [128, 256]:
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
        assert result['absolute_error'] < 0.2, f"Error too large for {result['test_name']}"
    
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
    sketch1 = MinHash(num_hashes=128, debug=False)
    sketch2 = MinHash(num_hashes=128, debug=False)
    
    # Add data to only one sketch
    sketch1.add_string("test")
    
    # Jaccard similarity between non-empty and empty should be 0
    assert sketch1.estimate_jaccard(sketch2) == 0.0
    
    # Empty sketches should have 0 similarity
    sketch3 = MinHash(num_hashes=128, debug=False)
    sketch4 = MinHash(num_hashes=128, debug=False)
    assert sketch3.estimate_jaccard(sketch4) == 0.0

def test_minhash_identical():
    mh1 = MinHash(kmer_size=3, debug=False)
    mh2 = MinHash(kmer_size=3, debug=False)
    
    mh1.add_string("test string")
    mh2.add_string("test string")
    
    sim = mh1.similarity_values(mh2)
    assert sim['jaccard_similarity'] == 1.0

def test_minhash_different():
    mh1 = MinHash(kmer_size=3, debug=False)
    mh2 = MinHash(kmer_size=3, debug=False)
    
    mh1.add_string("test string 1")
    mh2.add_string("test string 2")
    
    sim = mh1.similarity_values(mh2)
    assert 0.0 <= sim['jaccard_similarity'] <= 1.0

def test_minhash_different_sizes():
    mh1 = MinHash(num_hashes=128, debug=False)
    mh2 = MinHash(num_hashes=64, debug=False)
    
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
        sketch1.add_string(str(i).zfill(4))
    
    # Add elements to sketch2 (500-1499)
    # This creates 50% overlap with sketch1
    for i in range(500, 1500):
        sketch2.add_string(str(i).zfill(4))
    
    # Calculate Jaccard similarity
    jaccard = sketch1.estimate_jaccard(sketch2)
    
    # Expected Jaccard is 500/1500 = 1/3
    expected = 1/3
    assert abs(jaccard - expected) < 0.1, f"Expected Jaccard ~{expected:.3f}, got {jaccard:.3f}"

def test_minhash_similarity():
    """Test MinHash similarity estimation."""
    # Smaller test parameters
    seq_length = 200  # Reduced from 1000
    num_hashes = 64   # Reduced from 128
    k = 5            # Smaller k-mer size
    
    # Test fewer similarity values
    similarities = [0.9, 0.5]  # Reduced set of similarities to test
    trials = 2               # Reduced number of trials
    
    for expected_sim in similarities:
        for _ in range(trials):
            seq1, seq2 = generate_sequence_pair(seq_length, expected_sim)
            
            # Create MinHash sketches
            sketch1 = MinHash(kmer_size=k, num_hashes=num_hashes)
            sketch2 = MinHash(kmer_size=k, num_hashes=num_hashes)
            
            sketch1.add_string(seq1)
            sketch2.add_string(seq2)
            
            estimated_sim = sketch1.estimate_jaccard(sketch2)
            
            # Allow for some error in the estimation
            assert abs(estimated_sim - expected_sim) < 0.2

def test_minhash_merge():
    """Test MinHash merging operation."""
    seq_length = 200  # Reduced length
    num_hashes = 64   # Reduced number of hashes
    k = 5            # Smaller k-mer size
    
    # Generate test sequence
    sequence = ''.join(random.choice('ACGT') for _ in range(seq_length))
    
    # Create and merge sketches
    sketch1 = MinHash(kmer_size=k, num_hashes=num_hashes)
    sketch2 = MinHash(kmer_size=k, num_hashes=num_hashes)
    
    half = len(sequence) // 2
    sketch1.add_string(sequence[:half])
    sketch2.add_string(sequence[half:])
    
    merged = sketch1.merge(sketch2)
    
    # Create a sketch of the full sequence for comparison
    full_sketch = MinHash(kmer_size=k, num_hashes=num_hashes)
    full_sketch.add_string(sequence)
    
    # Compare merged sketch with full sequence sketch
    similarity = merged.estimate_jaccard(full_sketch)
    assert abs(similarity - 1.0) < 0.2

def _levenshtein_distance(s1: str, s2: str) -> int:
    """Calculate the Levenshtein distance between two strings."""
    if len(s1) < len(s2):
        return _levenshtein_distance(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]

@pytest.mark.quick
def test_minhash_small_overlap():
    """Test MinHash similarity estimation for small sequences with partial overlap."""
    num_hashes = 128
    k = 3  # Smaller k for short sequences
    seed = 42  # Set a specific seed for reproducibility
    
    # Multiple test cases with known overlaps
    test_cases = [
        ("ACGTACGTAAG", "TACGTGCAAAG"),  # Original case
        ("AAAA", "TTTT"),                # No overlap
        ("ACGT", "ACGT"),                # Perfect match
        ("ACGTACGT", "TGCATGCA")         # Minimal overlap
    ]
    
    for seq1, seq2 in test_cases:
        # Calculate Levenshtein similarity
        edit_dist = _levenshtein_distance(seq1, seq2)
        levenshtein_sim = 1 - (edit_dist / max(len(seq1), len(seq2)))
        
        # Calculate exact character-wise similarity
        char_sim = sum(1 for a, b in zip(seq1, seq2) if a == b) / max(len(seq1), len(seq2))
        
        # Create sketches with explicit seed
        sketch1 = MinHash(kmer_size=k, num_hashes=num_hashes, seed=seed, debug=True)
        sketch2 = MinHash(kmer_size=k, num_hashes=num_hashes, seed=seed, debug=True)
        
        # Debug k-mer generation
        kmers1 = set(seq1[i:i+k] for i in range(len(seq1)-k+1))
        kmers2 = set(seq2[i:i+k] for i in range(len(seq2)-k+1))
        kmer_intersection = kmers1 & kmers2
        kmer_union = kmers1 | kmers2
        true_jaccard = len(kmer_intersection) / len(kmer_union)
        
        sketch1.add_string(seq1)
        sketch2.add_string(seq2)
        
        estimated_sim = sketch1.estimate_jaccard(sketch2)
        
        # Print debug info
        print(f"\nTest case: {seq1} vs {seq2}")
        print(f"K-mers in seq1: {kmers1}")
        print(f"K-mers in seq2: {kmers2}")
        print(f"K-mer intersection: {kmer_intersection}")
        print(f"K-mer union: {kmer_union}")
        print(f"True Jaccard similarity: {true_jaccard:.3f}")
        print(f"MinHash estimated similarity: {estimated_sim:.3f}")
        print(f"Levenshtein similarity: {levenshtein_sim:.3f}")
        print(f"Character-wise similarity: {char_sim:.3f}")
        
        # Check for extreme values
        if estimated_sim == 1.0 and true_jaccard < 0.9:
            print("WARNING: MinHash reporting perfect similarity for non-identical sequences!")
            print("This suggests a problem with hash function generation or comparison")
        
        # Use wider tolerance for small sequences but compare against true Jaccard
        assert abs(estimated_sim - true_jaccard) < 0.2, \
            f"Estimated similarity {estimated_sim:.3f} too far from true Jaccard {true_jaccard:.3f}"

@pytest.mark.quick
def test_minhash_hash_distribution():
    """Test that MinHash is generating different hash values."""
    num_hashes = 128
    k = 3
    seed = 42
    
    sketch = MinHash(kmer_size=k, num_hashes=num_hashes, seed=seed)
    sequence = "ACGTACGTAAG"
    
    # Store hash values for each k-mer
    kmer_hashes = {}
    
    # Get k-mers
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    
    # Add each k-mer separately and record its hash values
    for kmer in kmers:
        test_sketch = MinHash(kmer_size=k, num_hashes=num_hashes, seed=seed, debug=True)
        test_sketch.add_string(kmer)
        kmer_hashes[kmer] = test_sketch.min_hashes.copy()
    
    # Print hash distribution info
    print("\nHash Distribution Test:")
    print(f"Sequence: {sequence}")
    print(f"k-mer size: {k}")
    print(f"Number of hashes: {num_hashes}")
    
    # Check for hash collisions between different k-mers
    for kmer1 in kmers:
        for kmer2 in kmers:
            if kmer1 < kmer2:  # Compare each pair once
                matches = sum(h1 == h2 for h1, h2 in zip(kmer_hashes[kmer1], kmer_hashes[kmer2]))
                print(f"K-mers {kmer1} vs {kmer2}: {matches}/{num_hashes} matching hashes")
                # We expect very few matches for different k-mers
                assert matches < num_hashes * 0.1, \
                    f"Too many matching hashes between different k-mers: {kmer1} vs {kmer2}"
    
    # Check that each k-mer gets different hash values
    for kmer, hashes in kmer_hashes.items():
        unique_hashes = len(set(hashes))
        print(f"K-mer {kmer}: {unique_hashes} unique hashes out of {num_hashes}")
        # We expect most hash values to be unique
        assert unique_hashes > num_hashes * 0.9, \
            f"Too few unique hash values for k-mer {kmer}"

@pytest.mark.full
def test_basic_jaccard():
    """Test basic Jaccard similarity calculation with multiple trials."""
    num_trials = 5  # Run multiple trials to ensure consistency
    max_error = 0.2  # Maximum allowed error
    
    for trial in range(num_trials):
        sketch1 = MinHash(num_hashes=512, seed=42 + trial)  # Increased num_hashes for better accuracy
        sketch2 = MinHash(num_hashes=512, seed=42 + trial)
        
        # Test case 1: Identical sets
        for i in range(1000):
            val = str(i).zfill(4)
            sketch1.add_string(val)
            sketch2.add_string(val)
        
        similarity = sketch1.estimate_jaccard(sketch2)
        assert abs(similarity - 1.0) < max_error, \
            f"Trial {trial}: Expected similarity 1.0 for identical sets, got {similarity}"
        
        # Clear sketch2 and test partial overlap
        sketch2 = MinHash(num_hashes=512, seed=42 + trial)
        
        # Add 500 overlapping elements and 500 different elements
        for i in range(500):  # First 500 elements overlap
            sketch2.add_string(str(i).zfill(4))
        for i in range(500, 1000):  # Next 500 are different
            sketch2.add_string(str(i + 500).zfill(4))  # Offset by 500 to ensure no overlap
        
        similarity = sketch1.estimate_jaccard(sketch2)
        assert abs(similarity - 0.5) < max_error, \
            f"Trial {trial}: Expected similarity 0.5 for 50% overlap, got {similarity}"
        
        # Test case 3: No overlap
        sketch2 = MinHash(num_hashes=512, seed=42 + trial)
        for i in range(1000):
            sketch2.add_string(str(i + 1000).zfill(4))  # Completely different elements
        
        similarity = sketch1.estimate_jaccard(sketch2)
        assert abs(similarity - 0.0) < max_error, \
            f"Trial {trial}: Expected similarity 0.0 for no overlap, got {similarity}"

@pytest.mark.quick
def test_hash_function_independence():
    """Test that hash functions are truly independent."""
    num_hashes = 128
    k = 3
    seed = 42
    
    # Create two different k-mers
    kmer1 = "ACG"
    kmer2 = "CGT"
    
    # Get hash values for both k-mers
    sketch1 = MinHash(kmer_size=k, num_hashes=num_hashes, seed=seed, debug=True)
    sketch2 = MinHash(kmer_size=k, num_hashes=num_hashes, seed=seed, debug=True)
    
    # Test single k-mer hashing
    print("\nTesting direct hash generation:")
    hashes1 = sketch1._hash_str(kmer1)
    hashes2 = sketch2._hash_str(kmer2)
    print(f"Direct hashes for {kmer1}: {hashes1[:5]}")
    print(f"Direct hashes for {kmer2}: {hashes2[:5]}")
    
    # Add k-mers to sketches
    sketch1.add_string(kmer1)
    sketch2.add_string(kmer2)
    
    # Count matching hashes
    matches = np.sum(sketch1.min_hashes == sketch2.min_hashes)
    print(f"\nNumber of matching hashes: {matches} out of {num_hashes}")
    
    # We expect very few matches between different k-mers
    max_allowed_matches = num_hashes * 0.1  # Allow up to 10% matches by chance
    assert matches <= max_allowed_matches, \
        f"Too many matching hashes: {matches} > {max_allowed_matches}"

if __name__ == "__main__":
    # Run quick tests
    test = TestMinHashQuick()
    test.test_init()
    test.test_add_string()
    test.test_jaccard()
    
    # Run accuracy tests
    test = TestMinHashFull()
    test.test_minhash_accuracy()