#!/usr/bin/env python
import random
import sys
import math

# given the starting position of the first interval in the bedfile (start_pos), the number of intervals/lines desired (count), the distance between the start and end of each interval (size), the chromosome value (chrom), and the distance between the end position of one interval and the start position of the next interval, generate a bedfile


def generate_bed_file(start_pos, interval_count, size, distance, output_file, chrom=10, seed=0):
    """
    Generate a BED file with specified intervals.

    Args:
        start_pos (int): The starting position of the first interval.
        interval_count (int): The number of intervals to generate.
        size (int): The size of each interval.
        distance (int): The distance between the start position of each interval.
        output_file (str): The path to the output BED file.
        chrom (int, optional): The chromosome number. Defaults to 10.
        seed (int, optional): The random seed. Defaults to 0.

    Returns:
        None
    """
    with open(output_file, 'w') as file:
        file.write("#chrom\tstart\tend\n")
        for _ in range(interval_count):
            random.seed( _ + seed)
            start = start_pos + (_ * random.randrange(math.floor(.75*distance), math.floor(1.25*distance)))
            end = start + size - random.randrange(0,math.floor(.5*size))
            file.write(f"{chrom}\t{start}\t{end}\n")


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage: python determined_bed.py <input csv>")
        sys.exit(1)
    with open(sys.argv[1], "r") as filein:
        # Assume header
        lines = filein.readlines()[1:]
        for line in lines:
            interval_vals = line.strip().split(',')
            generate_bed_file(
                start_pos=int(interval_vals[0]),
                interval_count=int(interval_vals[1]),
                size=int(interval_vals[2]),
                distance=int(interval_vals[3]),
                output_file=interval_vals[4],
                chrom=int(interval_vals[5]),
                seed=int(interval_vals[6])
            )


# <starting position of 1st interval> <number of intervals> <uniform size of each interval> <uniform distance between each interval start position> <output file name>  [chromosome_number]")
# if len(sys.argv) > 6:
    #     chrom = sys.argv[6]

 # start_pos = int(sys.argv[1])
    # interval_count = sys.argv[2]
    # size = sys.argv[3]
    # distance = sys.argv[4]
    # output_file = sys.argv[5]
    # chrom = 10

    