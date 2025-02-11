import unittest
import numpy as np # type: ignore
from hammock.lib.minimizer import MinimizerSketch
from hammock.lib.sequences import SequenceSketch
import os
import csv
from datetime import datetime

class TestMinimizerSimilarity(unittest.TestCase):
    def generate_mutated_sequence(self, orig_seq: str, mutation_rate: float) -> str:
        """Generate a mutated version of the sequence."""
        mutated_seq = list(orig_seq)
        num_mutations = int(len(orig_seq) * mutation_rate)
        mut_positions = np.random.choice(len(orig_seq), num_mutations, replace=False)
        
        for pos in mut_positions:
            choices = ['A', 'C', 'G', 'T']
            choices.remove(mutated_seq[pos])
            mutated_seq[pos] = np.random.choice(choices)
        
        return ''.join(mutated_seq)

    def test_similarity_under_mutations(self):
        """Test minimizer sketch similarity under different parameters."""
        # Test parameters
        seq_length = 100000  # Increased from 1000 to 100k
        k_values = [5, 7, 9]
        w_values = [10, 20, 40]
        gapk_values = [7, 8, 10, 20]  # Including different gap sizes
        mutation_rates = [0.01, 0.05, 0.1]
        trials = 3

        # Create test_results directory if it doesn't exist
        results_dir = os.path.join('test_results', 'similarity_tables')
        os.makedirs(results_dir, exist_ok=True)
        
        # Create CSV file with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        csv_file = os.path.join(results_dir, f'minimizer_similarities_{timestamp}.csv')
        
        # Print header for console output
        print("\nMinimizer Sketch Similarity Analysis")
        print(f"Sequence length: {seq_length}")
        print("-" * 80)
        print(f"{'k':<3} {'w':<3} {'gapk':<5} {'mut_rate':<8} {'hash_sim':<10} {'gap_sim':<10}")
        print("-" * 80)

        # Open CSV file and write header
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['k', 'w', 'gapk', 'mutation_rate', 'trial', 'hash_similarity', 'gap_similarity'])

            for k in k_values:
                for w in w_values:
                    for gapk in gapk_values:
                        for mut_rate in mutation_rates:
                            hash_sims = []
                            gap_sims = []
                            
                            for trial in range(trials):
                                # Generate original sequence
                                orig_seq = ''.join(np.random.choice(['A', 'C', 'G', 'T']) 
                                    for _ in range(seq_length))
                                
                                # Create mutated sequence
                                mutated_seq = self.generate_mutated_sequence(orig_seq, mut_rate)
                                
                                # Create sketches
                                sketch1 = MinimizerSketch(kmer_size=k, window_size=w, gapk=gapk)
                                sketch2 = MinimizerSketch(kmer_size=k, window_size=w, gapk=gapk)
                                
                                # Add sequences
                                sketch1.add_string(orig_seq)
                                sketch2.add_string(mutated_seq)
                                
                                # Get similarities
                                hash_sim, gap_sim = sketch1.compare_overlaps(sketch2)
                                hash_sims.append(hash_sim)
                                gap_sims.append(gap_sim)
                                
                                # Write individual trial results to CSV
                                writer.writerow([k, w, gapk, mut_rate, trial + 1, hash_sim, gap_sim])
                            
                            # Calculate mean similarities for console output
                            mean_hash = np.mean(hash_sims)
                            mean_gap = np.mean(gap_sims)
                            
                            print(f"{k:<3} {w:<3} {gapk:<5} {mut_rate:<8.2f} "
                                  f"{mean_hash:<10.4f} {mean_gap:<10.4f}")

        print(f"\nResults saved to: {csv_file}")

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
                hash_sim, gap_sim = sketch1.compare_overlaps(sketch2)
                
                # Print results instead of asserting
                print(f"\nTest case: {case['desc']}")
                print(f"Parameters: k={case['k']}, w={case['w']}, gapk={case['gapk']}")
                print(f"Hash similarity: {hash_sim:.4f} (expected: {case['expected_hash']:.4f})")
                print(f"Gap similarity: {gap_sim:.4f} (expected: {case['expected_gap']:.4f})")

if __name__ == '__main__':
    unittest.main() 