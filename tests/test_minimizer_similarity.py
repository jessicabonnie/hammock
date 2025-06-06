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

        # Run tests and write results
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            self._write_csv_header(writer)
            self._run_parameter_combinations(writer, params)

        print(f"\nResults saved to: {csv_file}")

    def _write_csv_header(self, writer: csv.writer) -> None:
        """Write the CSV header."""
        writer.writerow(['k', 'w', 'trial', 
                        # 'hash_similarity', 'hash_with_ends_similarity', 
                        'jaccard_similarity',
                        'jaccard_similarity_with_ends',
                        'char_similarity', 'edit_similarity', 
                        'expected_similarity'])

    def _run_parameter_combinations(self, writer: csv.writer, params: dict) -> None:
        """Run tests for all parameter combinations."""
        for k in params['k_values']:
            for w in params['w_values']:
                for expected_sim in params['expected_sims']:
                    self._run_trials(writer, k, w, expected_sim, params)

    def _run_trials(self, writer: csv.writer, k: int, w: int, 
                    expected_sim: float, params: dict) -> None:
        """Run multiple trials for a specific parameter combination."""
        for trial in range(params['trials']):
            self._run_single_trial(writer, k, w, expected_sim, 
                                 params['seq_length'], trial)

    def _run_single_trial(self, writer: csv.writer, k: int, w: int, 
                         expected_sim: float, seq_length: int, trial: int, debug: bool = False) -> None:
        """Run a single trial and record results."""
        # Generate sequence pair
        seq1, seq2 = self.generate_sequence_pair(seq_length, expected_sim)

        if debug:
            self._print_debug_header(k, w, expected_sim, seq1, seq2)

        # Calculate similarities
        char_sim = sum(1 for a, b in zip(seq1, seq2) if a == b) / len(seq1)
        edit_dist = self._levenshtein_distance(seq1, seq2)
        edit_sim = 1 - (edit_dist / len(seq1))

        # Create sketches and compare
        sketch1 = MinimizerSketch(kmer_size=k, window_size=w, debug=False)
        sketch2 = MinimizerSketch(kmer_size=k, window_size=w, debug=False)
        
        sketch1.add_string(seq1)
        sketch2.add_string(seq2)
        
        sim = sketch1.similarity_values(sketch2)

        if debug:
            self._print_debug_results(char_sim, edit_sim, sim, expected_sim)

        # Add assertions to verify similarity values
        # assert 0 <= sim['hash_similarity'] <= 1, f"Hash similarity {sim['hash_similarity']} out of range [0,1]"
        # assert 0 <= sim['hash_with_ends_similarity'] <= 1, f"Hash+ends similarity {sim['hash_with_ends_similarity']} out of range [0,1]"
        assert 0 <= sim['jaccard_similarity'] <= 1, f"Jaccard similarity {sim['jaccard_similarity']} out of range [0,1]"
        assert 0 <= sim['jaccard_similarity_with_ends'] <= 1, f"Jaccard+ends similarity {sim['jaccard_similarity_with_ends']} out of range [0,1]"
        # Verify character similarity matches expected similarity within tolerance
        tolerance = 0.1
        assert abs(char_sim - expected_sim) < tolerance, \
            f"Character similarity {char_sim} too far from expected {expected_sim}"
        
        # Verify similarity metrics are correlated
        # Jaccard similarity should be roughly similar to character similarity
        assert abs(sim['jaccard_similarity'] - char_sim) < 0.3, \
            f"Jaccard similarity {sim['jaccard_similarity']} too different from char similarity {char_sim}"
        
        # Hash with ends similarity should be at least as high as hash similarity
        tolerance = 0.1
        assert abs(sim['jaccard_similarity_with_ends'] - sim['jaccard_similarity']) < tolerance, \
            f"Jaccard+ends similarity {sim['jaccard_similarity_with_ends']} too different from jaccard similarity {sim['jaccard_similarity']}"

        # Record results
        self._write_results(writer, k, w, trial, sim, char_sim, edit_sim, expected_sim)

    def _print_debug_header(self, k: int, w: int, expected_sim: float, 
                           seq1: str, seq2: str) -> None:
        """Print debug information header."""
        print(f"\n{'='*80}")
        print(f"Debug output for k={k}, w={w}, expected_sim={expected_sim:.2f}")
        print(f"{'='*80}")
        print(f"Sequence samples:")
        print(f"Seq1 (first 50): {seq1[:50]}...")
        print(f"Seq2 (first 50): {seq2[:50]}...")

    def _print_debug_results(self, char_sim: float, edit_sim: float, sim: dict,
                            expected_sim: float) -> None:
        """Print detailed results for debug mode."""
        print(f"\nSimilarity Results:")
        print(f"Character-wise: {char_sim:.4f}")
        print(f"Edit Distance:  {edit_sim:.4f}")
        # print(f"Hash:          {sim['hash_similarity']:.4f}")
        # print(f"Hash+Ends:     {sim['hash_with_ends_similarity']:.4f}")
        print(f"Jaccard:       {sim['jaccard_similarity']:.4f}")
        print(f"Jaccard+Ends:  {sim['jaccard_similarity_with_ends']:.4f}")
        print(f"Expected:      {expected_sim:.4f}")
        print(f"{'='*80}\n")

    def _write_results(self, writer: csv.writer, k: int, w: int, trial: int,
                      sim: dict, char_sim: float, edit_sim: float, expected_sim: float) -> None:
        """Write results to CSV file."""
        writer.writerow([k, w, trial + 1, 
                        # f"{sim['hash_similarity']:.4f}", f"{sim['hash_with_ends_similarity']:.4f}",
                        f"{sim['jaccard_similarity']:.4f}",
                        f"{sim['jaccard_similarity_with_ends']:.4f}",
                        f"{char_sim:.4f}", f"{edit_sim:.4f}",
                        f"{expected_sim:.4f}"])

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
                'expected_hash': 1.0,
                'expected_hash_ends': 1.0,
                'desc': 'identical sequences'
            },
            # Completely different sequences should have similarity close to 0
            {
                'seq1': 'AAAAAAAAAA',
                'seq2': 'CCCCCCCCCC',
                'k': 4,
                'w': 4,
                'expected_hash': 0.0,
                'expected_hash_ends': 0.0,
                'desc': 'completely different sequences'
            },
            # Single mutation in middle
            {
                'seq1': 'ACGTACGTACGT',
                'seq2': 'ACGTACCTACGT',
                'k': 4,
                'w': 4,
                'expected_hash': 0.8,
                'expected_hash_ends': 0.8,
                'desc': 'single mutation'
            },
            # Same ends, different middle
            {
                'seq1': 'ACGTACGTACGT',
                'seq2': 'ACGT' + 'TTTT' + 'ACGT',
                'k': 4,
                'w': 4,
                'expected_hash': 0.5,
                'expected_hash_ends': 0.8,  # Higher due to matching ends
                'desc': 'same ends, different middle'
            }
        ]

        for case in test_cases:
            with self.subTest(desc=case['desc']):
                sketch1 = MinimizerSketch(kmer_size=case['k'], window_size=case['w'])
                sketch2 = MinimizerSketch(kmer_size=case['k'], window_size=case['w'])
                
                sketch1.add_string(case['seq1'])
                sketch2.add_string(case['seq2'])
                
                sim = sketch1.similarity_values(sketch2)
                # self.assertAlmostEqual(sim['hash_similarity'], case['expected_hash'], delta=0.2)
                # self.assertAlmostEqual(sim['hash_with_ends_similarity'], case['expected_hash_ends'], delta=0.2)
                self.assertAlmostEqual(sim['jaccard_similarity'], case['expected_hash'], delta=0.2)
                self.assertAlmostEqual(sim['jaccard_similarity_with_ends'], case['expected_hash_ends'], delta=0.2)

    def test_simple_similarity(self):
        """Test basic similarity calculation with simple sequences."""
        sketch1 = MinimizerSketch(kmer_size=4, window_size=5)
        sketch2 = MinimizerSketch(kmer_size=4, window_size=5)
        
        # Test identical sequences
        seq = "ACGTACGTACGT"
        sketch1.add_string(seq)
        sketch2.add_string(seq)
        
        sim = sketch1.similarity_values(sketch2)
        # self.assertAlmostEqual(sim['hash_similarity'], 1.0)
        # self.assertAlmostEqual(sim['hash_with_ends_similarity'], 1.0)
        self.assertAlmostEqual(sim['jaccard_similarity'], 1.0)
        self.assertAlmostEqual(sim['jaccard_similarity_with_ends'], 1.0)

        # Test completely different sequences
        sketch1 = MinimizerSketch(kmer_size=4, window_size=5)
        sketch2 = MinimizerSketch(kmer_size=4, window_size=5)
        
        sketch1.add_string("AAAAAAAAAAAA")
        sketch2.add_string("CCCCCCCCCCCC")
        
        sim = sketch1.similarity_values(sketch2)
        # self.assertAlmostEqual(sim['hash_similarity'], 0.0, delta=0.1)
        # self.assertAlmostEqual(sim['hash_with_ends_similarity'], 0.0, delta=0.1)
        self.assertAlmostEqual(sim['jaccard_similarity'], 0.0, delta=0.1)
        self.assertAlmostEqual(sim['jaccard_similarity_with_ends'], 0.0, delta=0.1)

    def _test_pair(self, seq1: str, seq2: str, desc: str):
        """Test a pair of sequences with different similarity metrics."""
        sketch1 = MinimizerSketch(kmer_size=4, window_size=5)
        sketch2 = MinimizerSketch(kmer_size=4, window_size=5)
        
        sketch1.add_string(seq1)
        sketch2.add_string(seq2)
        
        sim = sketch1.similarity_values(sketch2)
        
        # Print results for analysis
        print(f"\nTest case: {desc}")
        print(f"Sequence 1: {seq1}")
        print(f"Sequence 2: {seq2}")
        # print(f"Hash similarity: {sim['hash_similarity']:.4f}")
        # print(f"Hash+ends similarity: {sim['hash_with_ends_similarity']:.4f}")
        print(f"Jaccard similarity: {sim['jaccard_similarity']:.4f}")
        print(f"Jaccard+ends similarity: {sim['jaccard_similarity_with_ends']:.4f}")

if __name__ == '__main__':
    unittest.main() 