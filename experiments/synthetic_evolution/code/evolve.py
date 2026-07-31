#!/usr/bin/env python3
"""
evolve.py - Simulate evolution of BED file intervals

This program simulates evolutionary changes to a BED file through multiple generations.
Four types of evolution are supported:
- Type A: Intervals appear or disappear (insertion/deletion)
- Type A2: Intervals only disappear (deletion only)
- Type B: Intervals jitter (small position changes)
- Type C: Both A and B occur with equal probability
"""

import sys
import random
import argparse
from pathlib import Path


def read_bed_file(input_file):
    """Read a BED file and return a list of intervals."""
    intervals = []
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                # Keep any additional columns
                extra = parts[3:] if len(parts) > 3 else []
                intervals.append((chrom, start, end, extra))
    return intervals


def write_bed_file(output_file, intervals):
    """Write intervals to a BED file."""
    with open(output_file, 'w') as f:
        for interval in intervals:
            chrom, start, end, extra = interval
            if extra:
                f.write(f"{chrom}\t{start}\t{end}\t{chr(9).join(extra)}\n")
            else:
                f.write(f"{chrom}\t{start}\t{end}\n")


def evolve_type_a(interval, change_prob, chromosome_options,
                  genome_range=1000000):
    """
    Type A evolution: Intervals appear or disappear.
    Returns None if interval should disappear, or a new random interval if it should change.
    Replacement intervals are drawn uniformly from [0, genome_range].
    """
    if random.random() < change_prob:
        # 50% chance to disappear, 50% chance to mutate to a new interval
        if random.random() < 0.5:
            return None  # Disappear
        else:
            # Mutate to a completely new interval
            chrom = random.choice(chromosome_options)
            interval_length = random.choice([100, 500, 1000, 5000, 10000])
            start = random.randint(0, genome_range)
            end = start + interval_length
            return (chrom, start, end, interval[3])
    return interval


def evolve_type_a2(interval, change_rate):
    """
    Type A2 evolution: Intervals only disappear (deletion only).
    Returns None if interval should disappear, otherwise returns unchanged interval.
    """
    if random.random() < change_rate:
        return None  # Disappear
    return interval


def evolve_type_b(interval, change_rate, jitter_amount=100):
    """
    Type B evolution: Intervals jitter by a small amount.
    Both start and end positions are adjusted independently.
    """
    if random.random() < change_rate:
        chrom, start, end, extra = interval
        
        # Jitter the start and end positions independently
        start_jitter = random.randint(-jitter_amount, jitter_amount)
        end_jitter = random.randint(-jitter_amount, jitter_amount)
        
        new_start = max(0, start + start_jitter)
        new_end = max(new_start + 1, end + end_jitter)  # Ensure end > start
        
        return (chrom, new_start, new_end, extra)
    return interval


def evolve_type_c(interval, change_rate, chromosome_options,
                  jitter_amount=100, genome_range=1000000):
    """
    Type C evolution: Both A and B occur with equal probability.
    """
    if random.random() < change_rate:
        # Equal probability of type A or type B
        if random.random() < 0.5:
            return evolve_type_a(interval, 1.0, chromosome_options,
                                 genome_range=genome_range)
        else:
            return evolve_type_b(interval, 1.0, jitter_amount)
    return interval


def evolve_generation(intervals, evolution_type, change_rate,
                      chromosome_options, jitter_amount, genome_range=1000000):
    """
    Evolve all intervals for one generation.
    """
    new_intervals = []

    for interval in intervals:
        if evolution_type == 'A':
            evolved = evolve_type_a(interval, change_rate, chromosome_options,
                                    genome_range=genome_range)
        elif evolution_type == 'A2':
            evolved = evolve_type_a2(interval, change_rate)
        elif evolution_type == 'B':
            evolved = evolve_type_b(interval, change_rate, jitter_amount)
        elif evolution_type == 'C':
            evolved = evolve_type_c(interval, change_rate, chromosome_options,
                                    jitter_amount, genome_range=genome_range)
        else:
            raise ValueError(f"Unknown evolution type: {evolution_type}")
        
        # Only add if interval didn't disappear (type A/A2 can return None)
        if evolved is not None:
            new_intervals.append(evolved)
    
    return new_intervals


def evolve(input_file, output_file, num_generations, change_rate, evolution_type,
           jitter_amount=100, chromosome_options=None, step=None, verbose=False,
           genome_range=1000000):
    """
    Main evolution function.
    
    Args:
        input_file: Path to input BED file
        output_file: Path to output BED file
        num_generations: Number of generations to evolve
        change_rate: Probability of a line being changed per generation (0.0 to 1.0)
        evolution_type: Type of evolution ('A', 'B', or 'C')
        jitter_amount: Amount of jitter for type B evolution (default 100bp)
        chromosome_options: List of valid chromosomes for type A evolution
        step: If provided, output intermediate files every 'step' generations
        verbose: Print progress information
    """
    if chromosome_options is None:
        chromosome_options = ["chr1", "chr2", "chr3", "chr4", "chr5"]
    
    # Validate parameters
    if not 0.0 <= change_rate <= 1.0:
        raise ValueError("change_rate must be between 0.0 and 1.0")
    
    if evolution_type not in ['A', 'A2', 'B', 'C']:
        raise ValueError("evolution_type must be 'A', 'A2', 'B', or 'C'")
    
    if num_generations < 0:
        raise ValueError("num_generations must be non-negative")
    
    if step is not None and step <= 0:
        raise ValueError("step must be positive")
    
    # Read initial intervals
    intervals = read_bed_file(input_file)
    initial_count = len(intervals)
    
    if verbose:
        print(f"Starting evolution with {initial_count} intervals")
        print(f"Evolution type: {evolution_type}")
        print(f"Change rateability: {change_rate}")
        print(f"Generations: {num_generations}")
        if evolution_type in ['B', 'C']:
            print(f"Jitter amount: ±{jitter_amount}bp")
        if step is not None:
            print(f"Saving intermediate files every {step} generation(s)")
    
    # Determine output files
    output_files = []
    if step is not None:
        # Generate intermediate output filenames
        # Extract base and extension from output_file
        if output_file.endswith('.bed'):
            base = output_file[:-4]
        else:
            base = output_file
        
        for gen in range(step, num_generations + 1, step):
            intermediate_file = f"{base}_{gen:03d}.bed"
            output_files.append((gen, intermediate_file))
    else:
        # Only output final file
        output_files.append((num_generations, output_file))
    
    # Evolve for specified number of generations
    for gen in range(1, num_generations + 1):
        intervals = evolve_generation(intervals, evolution_type, change_rate,
                                     chromosome_options, jitter_amount,
                                     genome_range=genome_range)
        if verbose:
            print(f"Generation {gen}: {len(intervals)} intervals")
        
        # Check if we need to save this generation
        for target_gen, target_file in output_files:
            if gen == target_gen:
                write_bed_file(target_file, intervals)
                if verbose:
                    print(f"  → Saved to: {target_file}")
    
    if verbose:
        print(f"\nEvolution complete!")
        print(f"Initial intervals: {initial_count}")
        print(f"Final intervals: {len(intervals)}")
        print(f"Change: {len(intervals) - initial_count:+d} ({100 * (len(intervals) - initial_count) / initial_count:+.1f}%)")
        if step is None:
            print(f"Output written to: {output_file}")
        else:
            print(f"Output files written: {len(output_files)} files")


def main():
    parser = argparse.ArgumentParser(
        description="Simulate evolution of BED file intervals",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Evolution types:
  A  - Intervals appear or disappear (insertion/deletion)
  A2 - Intervals only disappear (deletion only)
  B  - Intervals jitter by small amounts (position changes)
  C  - Both A and B occur with equal rate

Output filename format:
  The output file will be automatically named as:
    {output_prefix}_g{generations}_r{rate}_{type}.bed
  
  With --step, intermediate files are named:
    {output_prefix}_g{generations}_r{rate}_{type}_{gen:03d}.bed
  
  Example: "evolved_g10_r0.1_A.bed" for prefix "evolved", 10 generations, 0.1 rate, type A

Examples:
  # Type A evolution: appearance/disappearance (creates evolved_g5_r0.1_A.bed)
  python evolve.py input.bed evolved --generations 5 --prob 0.1 --type A
  
  # Type A2 evolution: deletion only (creates evolved_g5_r0.1_A2.bed)
  python evolve.py input.bed evolved --generations 5 --prob 0.1 --type A2
  
  # Type B evolution with custom jitter (creates test_g10_r0.05_B.bed)
  python evolve.py input.bed test --generations 10 --prob 0.05 --type B --jitter 50
  
  # Type C evolution with intermediate snapshots
  # Creates: evolved_g20_r0.2_C_005.bed, evolved_g20_r0.2_C_010.bed, ..., evolved_g20_r0.2_C_020.bed
  python evolve.py input.bed evolved --generations 20 --prob 0.2 --type C --step 5 --verbose
        """
    )
    
    parser.add_argument("input_file", 
                       help="Input BED file to evolve")
    parser.add_argument("output_prefix", 
                       help="Output filename prefix (full path will be {prefix}_g{gen}_r{prob}_{type}.bed)")
    parser.add_argument("-g", "--generations", 
                       type=int, required=True,
                       help="Number of generations to evolve")
    parser.add_argument("-p", "--prob", "--rate", "-r",
                       type=float, required=True,
                       help="Probability of a line being changed per generation (0.0-1.0)")
    parser.add_argument("-t", "--type", 
                       choices=['A', 'A2', 'B', 'C'], required=True,
                       help="Type of evolution: A (appear/disappear), A2 (deletion only), B (jitter), C (both A+B)")
    parser.add_argument("-j", "--jitter", 
                       type=int, default=100,
                       help="Jitter amount in base pairs for type B/C evolution (default: 100)")
    parser.add_argument("-c", "--chromosomes",
                       nargs='+',
                       default=["chr1", "chr2", "chr3", "chr4", "chr5"],
                       help="Valid chromosome names for type A evolution")
    parser.add_argument("--step",
                       type=int,
                       help="Output intermediate files every STEP generations (generation numbers padded to 3 digits)")
    parser.add_argument("--genome_range",
                       type=int, default=1000000,
                       help="Max start position per chromosome for replacement "
                            "intervals (must match createBED.py --genome_range). "
                            "Default 1e6.")
    parser.add_argument("-s", "--seed",
                       type=int,
                       help="Random seed for reproducibility")
    parser.add_argument("-v", "--verbose",
                       action="store_true",
                       help="Print progress information")

    args = parser.parse_args()
    
    # Set random seed if provided
    if args.seed is not None:
        random.seed(args.seed)
    
    # Check if input file exists
    if not Path(args.input_file).exists():
        print(f"Error: Input file '{args.input_file}' does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Construct output filename from prefix and parameters
    output_file = f"{args.output_prefix}_g{args.generations}_r{args.prob}_{args.type}.bed"
    
    # Run evolution
    try:
        evolve(args.input_file, output_file, args.generations, args.prob,
               args.type, args.jitter, args.chromosomes, args.step, args.verbose,
               genome_range=args.genome_range)
    except Exception as e:
        print(f"Error during evolution: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

