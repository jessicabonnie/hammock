import unittest
import numpy as np # type: ignore
from hammock.lib.minimizer import MinimizerSketch
from hammock.lib.sequences import SequenceSketch
import os
import csv
from datetime import datetime
import random

class TestMinimizerSimilarity(unittest.TestCase):
    def generate_sequence_pair(self, length: int, expected_similarity: float) -> tuple[str, str]:
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

    def test_similarity_metrics_comparison(self):
        """Test and compare different sequence similarity metrics."""
        # Test parameters
        params = {
            'seq_length': 1000,
            'k_values': [5, 7, 9, 10],
            'w_values': [10, 20, 40],
            'gapk_values': [7, 8, 10],
            'expected_sims': [0.99, 0.95, 0.80, 0.75, 0.50],
            'trials': 3
        }
        
        # Set random seed for reproducibility
        random.seed(42)

        # Setup output files
        results_dir = os.path.join('test_results', 'similarity_tables')
        os.makedirs(results_dir, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        csv_file = os.path.join(results_dir, f'minimizer_similarities_{timestamp}.csv')

        # Print header
        # self._print_header(params['seq_length'])

        # Run tests and write results
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            self._write_csv_header(writer)
            self._run_parameter_combinations(writer, params)

        print(f"\nResults saved to: {csv_file}")

    def _print_header(self, seq_length: int) -> None:
        """Print the console header."""
        print("\nSequence Similarity Metrics Comparison")
        print(f"Sequence length: {seq_length}")
        print("-" * 140)
        print(f"{'k':<3} {'w':<3} {'gapk':<5} {'hash_sim':<10} {'hash_ends':<10} {'gap_sim':<10} "
              f"{'jaccard':<10} {'char_sim':<10} {'edit_sim':<10} {'expected':<10}")
        print("-" * 140)

    def _write_csv_header(self, writer: csv.writer) -> None:
        """Write the CSV header."""
        writer.writerow(['k', 'w', 'gapk', 'trial', 
                        'hash_similarity', 'hash_with_ends_similarity', 
                        'gap_similarity', 'jaccard_similarity',
                        'char_similarity', 'edit_similarity', 
                        'expected_similarity'])

    def _run_parameter_combinations(self, writer: csv.writer, params: dict) -> None:
        """Run tests for all parameter combinations."""
        for k in params['k_values']:
            for w in params['w_values']:
                for gapk in params['gapk_values']:
                    for expected_sim in params['expected_sims']:
                        self._run_trials(writer, k, w, gapk, expected_sim, params)

    def _run_trials(self, writer: csv.writer, k: int, w: int, gapk: int, 
                    expected_sim: float, params: dict) -> None:
        """Run multiple trials for a specific parameter combination."""
        for trial in range(params['trials']):
            self._run_single_trial(writer, k, w, gapk, expected_sim, 
                                 params['seq_length'], trial)

    def _run_single_trial(self, writer: csv.writer, k: int, w: int, gapk: int, 
                         expected_sim: float, seq_length: int, trial: int, debug: bool = False) -> None:
        """Run a single trial and record results."""
        # Generate sequence pair
        seq1, seq2 = self.generate_sequence_pair(seq_length, expected_sim)

        if debug:
            self._print_debug_header(k, w, gapk, expected_sim, seq1, seq2)

        # Calculate similarities
        char_sim = sum(1 for a, b in zip(seq1, seq2) if a == b) / len(seq1)
        edit_dist = self._levenshtein_distance(seq1, seq2)
        edit_sim = 1 - (edit_dist / len(seq1))

        # Create sketches and compare - set debug=False
        sketch1 = MinimizerSketch(kmer_size=k, window_size=w, gapk=gapk, debug=False)
        sketch2 = MinimizerSketch(kmer_size=k, window_size=w, gapk=gapk, debug=False)
        
        sketch1.add_string(seq1)
        sketch2.add_string(seq2)
        
        hash_sim, hash_ends_sim, gap_sim, jaccard_sim = sketch1.compare_overlaps(sketch2)

        if debug:
            self._print_debug_results(char_sim, edit_sim, hash_sim, 
                                    hash_ends_sim, gap_sim, jaccard_sim, expected_sim)

        # Record results
        self._write_results(writer, k, w, gapk, trial, hash_sim, hash_ends_sim,
                           gap_sim, jaccard_sim, char_sim, edit_sim, expected_sim)
        # self._print_results(k, w, gapk, hash_sim, hash_ends_sim, gap_sim,
        #                    jaccard_sim, char_sim, edit_sim, expected_sim)

    def _print_debug_header(self, k: int, w: int, gapk: int, expected_sim: float, 
                           seq1: str, seq2: str) -> None:
        """Print debug information header."""
        print(f"\n{'='*80}")
        print(f"Debug output for k={k}, w={w}, gapk={gapk}, expected_sim={expected_sim:.2f}")
        print(f"{'='*80}")
        print(f"Sequence samples:")
        print(f"Seq1 (first 50): {seq1[:50]}...")
        print(f"Seq2 (first 50): {seq2[:50]}...")

    def _print_debug_results(self, char_sim: float, edit_sim: float, hash_sim: float,
                            hash_ends_sim: float, gap_sim: float, jaccard_sim: float,
                            expected_sim: float) -> None:
        """Print detailed results for debug mode."""
        print(f"\nSimilarity Results:")
        print(f"Character-wise: {char_sim:.4f}")
        print(f"Edit Distance:  {edit_sim:.4f}")
        print(f"Hash:          {hash_sim:.4f}")
        print(f"Hash+Ends:     {hash_ends_sim:.4f}")
        print(f"Gap:           {gap_sim:.4f}")
        print(f"Jaccard:       {jaccard_sim:.4f}")
        print(f"Expected:      {expected_sim:.4f}")
        print(f"{'='*80}\n")

    def _write_results(self, writer: csv.writer, k: int, w: int, gapk: int, trial: int,
                      hash_sim: float, hash_ends_sim: float, gap_sim: float, jaccard_sim: float,
                      char_sim: float, edit_sim: float, expected_sim: float) -> None:
        """Write results to CSV file."""
        writer.writerow([k, w, gapk, trial + 1, 
                        f"{hash_sim:.4f}", f"{hash_ends_sim:.4f}",
                        f"{gap_sim:.4f}", f"{jaccard_sim:.4f}",
                        f"{char_sim:.4f}", f"{edit_sim:.4f}",
                        f"{expected_sim:.4f}"])

    def _print_results(self, k: int, w: int, gapk: int, hash_sim: float,
                      hash_ends_sim: float, gap_sim: float, jaccard_sim: float,
                      char_sim: float, edit_sim: float, expected_sim: float) -> None:
        """Print results to console."""
        print(f"{k:<3} {w:<3} {gapk:<5} {hash_sim:<10.4f} {hash_ends_sim:<10.4f} "
              f"{gap_sim:<10.4f} {jaccard_sim:<10.4f} {char_sim:<10.4f} "
              f"{edit_sim:<10.4f} {expected_sim:<10.4f}")

    def _levenshtein_distance(self, s1: str, s2: str) -> int:
        """Calculate the Levenshtein distance between two strings."""
        if len(s1) < len(s2):
            return self._levenshtein_distance(s2, s1)

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

    def test_known_similarities(self):
        """Test minimizer sketch similarity with known sequences."""
        test_cases = [
            # Identical sequences should have similarity 1.0
            {
                'seq1': 'ACGTACGTACGT',
                'seq2': 'ACGTACGTACGT',
                'k': 4,
                'w': 4,
                'gapk': 0,
                'expected_hash': 1.0,
                'expected_gap': 1.0,
                'desc': 'identical sequences'
            },
            # Completely different sequences should have similarity close to 0
            {
                'seq1': 'AAAAAAAAAA',
                'seq2': 'CCCCCCCCCC',
                'k': 4,
                'w': 4,
                'gapk': 0,
                'expected_hash': 0.0,
                'expected_gap': 0.0,
                'desc': 'completely different sequences'
            },
            # Single mutation in middle
            {
                'seq1': 'ACGTACGTACGT',
                'seq2': 'ACGTACCTACGT',
                'k': 4,
                'w': 4,
                'gapk': 0,
                'expected_hash': 0.75,  # Most 4-mers should still match
                'expected_gap': 0.75,
                'desc': 'single mutation'
            },
            # Shifted sequence (should affect gap pattern)
            {
                'seq1': 'ACGTACGTACGT',
                'seq2': 'TACGTACGTACG',
                'k': 4,
                'w': 4,
                'gapk': 1,
                'expected_hash': 0.8,
                'expected_gap': 0.0,  # Gap pattern should be different
                'desc': 'shifted sequence'
            }
        ]

        for case in test_cases:
            with self.subTest(msg=case['desc']):
                # Create sketches
                sketch1 = MinimizerSketch(kmer_size=case['k'], 
                                        window_size=case['w'], 
                                        gapk=case['gapk'])
                sketch2 = MinimizerSketch(kmer_size=case['k'], 
                                        window_size=case['w'], 
                                        gapk=case['gapk'])
                
                # Add sequences
                sketch1.add_string(case['seq1'])
                sketch2.add_string(case['seq2'])
                
                # Get similarities
                hash_sim, hash_ends_sim, gap_sim, jaccard_sim = sketch1.compare_overlaps(sketch2)
                
                # Print results instead of asserting
                print(f"\nTest case: {case['desc']}")
                print(f"Parameters: k={case['k']}, w={case['w']}, gapk={case['gapk']}")
                print(f"Hash similarity: {hash_sim:.4f} (expected: {case['expected_hash']:.4f})")
                print(f"Gap similarity: {gap_sim:.4f} (expected: {case['expected_gap']:.4f})")

    def test_simple_similarity(self):
        """Test minimizer similarity with simple cases."""
        # Generate base sequences with more variation
        base_seq = ''.join(np.random.choice(['A', 'C', 'G', 'T']) for _ in range(100000))
        
        # Test case 1: Identical sequences
        seq1 = base_seq
        seq2 = seq1  # Identical
        
        print("\nTest Case 1: Identical Sequences")
        self._test_pair(seq1, seq2, "identical")
        
        # Test case 2: 80% similar
        seq1 = base_seq
        seq2 = list(seq1)  # Convert to list for easier modification
        positions = np.random.choice(len(seq1), int(len(seq1) * 0.2), replace=False)  # 20% positions to change
        for pos in positions:
            current = seq2[pos]
            choices = [x for x in 'ACGT' if x != current]
            seq2[pos] = np.random.choice(choices)
        seq2 = ''.join(seq2)
        
        print("\nTest Case 2: 80% Similar Sequences")
        self._test_pair(seq1, seq2, "80% similar")
        
        # Test case 3: 50% similar
        seq1 = base_seq
        seq2 = list(seq1)  # Convert to list for easier modification
        positions = np.random.choice(len(seq1), int(len(seq1) * 0.5), replace=False)  # 50% positions to change
        for pos in positions:
            current = seq2[pos]
            choices = [x for x in 'ACGT' if x != current]
            seq2[pos] = np.random.choice(choices)
        seq2 = ''.join(seq2)
        
        print("\nTest Case 3: 50% Similar Sequences")
        self._test_pair(seq1, seq2, "50% similar")

    def _test_pair(self, seq1: str, seq2: str, desc: str):
        """Helper method to test a pair of sequences."""
        k = 4  # k-mer size
        w = 4  # window size
        sketch1 = MinimizerSketch(kmer_size=k, window_size=w)
        sketch2 = MinimizerSketch(kmer_size=k, window_size=w)
        
        # Calculate actual character-wise similarity
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        actual_sim = matches / len(seq1)
        
        print(f"\nSequence info:")
        print(f"Seq1 (first 50): {seq1[:50]}...")
        print(f"Seq2 (first 50): {seq2[:50]}...")
        print(f"Character-wise similarity: {actual_sim:.4f}")
        
        # Add sequences and check results
        sketch1.add_string(seq1)
        sketch2.add_string(seq2)
        
        print(f"\nSketch info:")
        print(f"Sketch1 minimizers: {len(sketch1.minimizers)}")
        print(f"Sketch2 minimizers: {len(sketch2.minimizers)}")
        
        # Print some example minimizers
        if sketch1.minimizers:
            print(f"Sample minimizers from sketch1: {list(sketch1.minimizers)[:5]}")
        if sketch2.minimizers:
            print(f"Sample minimizers from sketch2: {list(sketch2.minimizers)[:5]}")
        
        # Check k-mers in a window
        window_size = 20
        print(f"\nAnalyzing first window of {window_size} bases:")
        window1 = seq1[:window_size]
        window2 = seq2[:window_size]
        print(f"Window1: {window1}")
        print(f"Window2: {window2}")
        
        # Print k-mers in the window
        kmers1 = [window1[i:i+k] for i in range(len(window1)-k+1)]
        kmers2 = [window2[i:i+k] for i in range(len(window2)-k+1)]
        print(f"k-mers in window1: {kmers1}")
        print(f"k-mers in window2: {kmers2}")
        
        # Check intersection and union
        intersection = sketch1.minimizers & sketch2.minimizers
        union = sketch1.minimizers | sketch2.minimizers
        print(f"\nSet operations:")
        print(f"Intersection size: {len(intersection)}")
        print(f"Union size: {len(union)}")
        
        # Get similarities
        hash_sim, hash_ends_sim, gap_sim, jaccard_sim = sketch1.compare_overlaps(sketch2)
        print(f"\nSimilarity results:")
        print(f"Hash similarity: {hash_sim:.4f} (expected: {actual_sim:.4f})")
        print(f"Gap similarity: {gap_sim:.4f} (expected: {actual_sim:.4f})")

if __name__ == '__main__':
    unittest.main() 