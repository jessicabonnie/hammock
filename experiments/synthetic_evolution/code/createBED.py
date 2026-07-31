#!/usr/bin/env python3
"""
createBED.py - Create a synthetic BED file with various interval sizes

This program creates a synthetic BED file that is large and contains both small 
and large intervals by randomly sampling from different interval lengths.
"""

import sys
import random
import argparse


# Interval length options: includes both small and large intervals
interval_length_options = [10, 50, 100, 500, 1000, 5000, 10000]

# Chromosome options for synthetic data
chromosome_options = ["chr1", "chr2", "chr3", "chr4", "chr5"]


def generate_random_bed_file(num_lines, output_file, genome_range=1000000):
    """
    Generate a synthetic BED file with random intervals.

    Args:
        num_lines:    Number of intervals to generate
        output_file:  Path to output BED file
        genome_range: Max start position per chromosome (default 1e6).
                      Coverage saturation: with N intervals of average
                      length L over C chromosomes of range R, the average
                      coverage fraction is N*L / (C*R). The original
                      1e6 default gives ~47× over-coverage at N=100k —
                      use 1e8 (~0.5%) for a biologically realistic
                      sparse setup.
    """
    with open(output_file, 'w') as f:
        for _ in range(num_lines):
            chromosome = random.choice(chromosome_options)
            interval_length = random.choice(interval_length_options)
            start = random.randint(0, genome_range)
            end = start + interval_length
            f.write(f"{chromosome}\t{start}\t{end}\n")

    print(f"Generated {num_lines} intervals in {output_file} "
          f"(genome_range={genome_range})")


def main():
    parser = argparse.ArgumentParser(
        description="Create a synthetic BED file with various interval sizes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create default synthetic BED file (10000 lines)
  python createBED.py
  
  # Create custom-sized BED file
  python createBED.py --num_lines 100000 --output_file large_synthetic.bed
  
  # Create small test file
  python createBED.py --num_lines 100 --output_file test.bed
        """
    )
    
    parser.add_argument("--num_lines",
                       type=int,
                       default=10000,
                       help="Number of lines/intervals to generate (default: 10000)")
    parser.add_argument("--output_file",
                       type=str,
                       default="synthetic-c.bed",
                       help="Output BED file name (default: synthetic-c.bed)")
    parser.add_argument("--genome_range",
                       type=int,
                       default=1000000,
                       help="Max start position per chromosome (default: 1e6). "
                            "Lower => more saturated; higher => sparser.")
    parser.add_argument("-s", "--seed",
                       type=int,
                       help="Random seed for reproducibility")
    
    args = parser.parse_args()
    
    # Set random seed if provided
    if args.seed is not None:
        random.seed(args.seed)
        print(f"Using random seed: {args.seed}")
    
    # Validate num_lines
    if args.num_lines <= 0:
        print(f"Error: num_lines must be positive, got {args.num_lines}", file=sys.stderr)
        sys.exit(1)
    
    # Generate the BED file
    generate_random_bed_file(args.num_lines, args.output_file,
                             genome_range=args.genome_range)


if __name__ == "__main__":
    main()

