#!/usr/bin/env python
import sys
import os

def basic_bedline(line):
    """Parse a single line from a BED file into chromosome, start, and end coordinates.
    
    Args:
        line: A string containing a single line from a BED file
        
    Returns:
        Tuple of (chrom, start, end) where:
            chrom: Chromosome name with 'chr' prefix removed if present 
            start: Integer start coordinate
            end: Integer end coordinate
            
    Raises:
        ValueError: If line has fewer than 3 tab or space-separated columns
    """
    columns = line.strip().split('\t')
    if len(columns) < 2:
        columns = line.strip().split(" ")
        if len(columns) < 2:
            raise ValueError("bedline: one of the lines in malformed")
    if columns[0].startswith('chr'):
        columns[0] = columns[0][3:]
    return columns[0], int(columns[1]), int(columns[2])

def bed_to_intervals(bedfile):
    """Read a BED file and convert to a list of interval tuples.
    
    Args:
        bedfile: Path to BED format file
        
    Returns:
        List of tuples, each containing (chrom, start, end) for an interval.
        Only includes intervals where chromosome is a number.
    """
    intervals = []
    with open(bedfile, 'r') as file:
        for line in file:
            interval = basic_bedline(line)
            if not is_int(interval[0]):
                continue
            intervals.append(basic_bedline(line))
    return intervals

def no_overlap(intervals):
    # print(len(intervals))
    new_intervals = []
    for interval in intervals:
        new_intervals.append(["88", interval[1], interval[2]])
    return new_intervals

def shrink_intervals(intervals, denominator=2, padding=True):
    """Shrink intervals by dividing their length by a denominator.
    
    For each interval, creates a new interval with the same start position but
    an end position that is (end-start)/denominator distance from the start.
    Optionally adds padding intervals to maintain the original number of points.
    
    Args:
        intervals: List of [chrom, start, end] interval lists
        denominator: Integer to divide interval lengths by. Defaults to 2
        padding: If True, adds padding intervals with chrom="55" to maintain 
                original number of points. Defaults to True
                
    Returns:
        List of [chrom, start, end] lists representing the shrunk intervals,
        with optional padding intervals appended
    """
    shrunk_intervals = []
    for i in range(len(intervals)):
        interval = intervals[i]
        chrom = interval[0]
        start = int(interval[1])
        end = int(interval[2])
        increment = (end-start) // denominator
        partial = start + increment
        shrunk_intervals.append([chrom, str(start), str(partial)])
        if padding:
            shrunk_intervals.append(["55", str(partial), str(end)])
    return shrunk_intervals
def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def subselect_intervals(intervals, n, padding=True):
    subselected_intervals = []
    for i in range(len(intervals)):
        if i % n == 0:
            # Add the interval to the subselected list if its index is a multiple of n
            subselected_intervals.append(intervals[i])
        else:
            if padding:
                # Add a modified interval with a 20 base pair shift if padding is enabled
                new_interval = (str(int(intervals[i][0])+30), intervals[i][1], intervals[i][2])
                subselected_intervals.append(new_interval)
    return subselected_intervals
# def every_nth(lst):
#     return lst[(n-1)::n]

def write_bedfile(intervals, output_file, sep='\t'):
    with open(output_file, 'w') as file:
        for interval in intervals:
            file.write(f"{interval[0]}{sep}{interval[1]}{sep}{interval[2]}\n")
    return output_file

def print_filenames(filenames):
    for filename in filenames:
        sys.stdout.write(os.path.abspath(filename) + "\n")
        sys.stdout.flush()
def generate_variant_bed_files(bedfile, output_prefix, ratio_denominators, output_dir='.'):
    """Generate multiple variant BED files from an input BED file.
    
    This function creates several modified versions of the input BED file:
    - Mode B files: Intervals are shrunk by dividing their length by each ratio denominator
      - With and without padding intervals to maintain the original number of points
    - Mode A files: Every nth interval is selected, where n is each ratio denominator
      - With and without padding with nonesense intervals to maintain the original number of intervals
    - No overlap file: All intervals are modified to ensure no overlaps
    
    Args:
        bedfile: Path to input BED file
        output_prefix: Prefix for output filenames
        ratio_denominators: List of integers to use as denominators for shrinking/selecting intervals
        output_dir: Directory to write output files to. Defaults to current directory.
        
    Returns:
        None. Prints absolute paths of all generated files.
    """
    intervals = bed_to_intervals(bedfile)
    output_files = []
    for i in ratio_denominators:
        output_files.append(write_bedfile(shrink_intervals(intervals, i), 
                                          os.path.join(output_dir, f"{output_prefix}_modeB_{i}.bed")))
        output_files.append(write_bedfile(shrink_intervals(intervals, i, padding=False), os.path.join(output_dir, f"{output_prefix}_nopad_modeB_{i}.bed")))
        output_files.append(write_bedfile(subselect_intervals(intervals, i), 
                                          os.path.join(output_dir, f"{output_prefix}_modeA_{i}.bed")))
        output_files.append(
            write_bedfile(subselect_intervals(intervals, i, padding=False),
                          os.path.join(output_dir, 
                                       f"{output_prefix}_nopad_modeA_{i}.bed")))

    output_files.append(write_bedfile(no_overlap(intervals), 
                                      os.path.join(output_dir, f"{output_prefix}_no_overlap.bed")))
    print_filenames(output_files)
if __name__ == "__main__":
    if len(sys.argv) not in [3, 4]:
        print("Usage: python script.py <input_bedfile> <ratio_denominators> [output_directory]")
        print("Example: python script.py input.bed 2,3,4,5 /path/to/output")
        sys.exit(1)
    
    input_bed = sys.argv[1]
    output_prefix = os.path.basename(os.path.splitext(input_bed)[0])
    ratio_denominators = [int(x) for x in sys.argv[2].split(',')]
    
    output_dir = '.'  # Default to current working directory
    if len(sys.argv) == 4:
        output_dir = sys.argv[3]
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    generate_variant_bed_files(input_bed, output_prefix, ratio_denominators, output_dir)
