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

def shrink_intervals(intervals, denominator=2, numerator=0, padding=True):
    """Shrink intervals by dividing their length by a denominator.
    
    For each interval, creates a new interval with the same start position but
    an end position that is numerator*(end-start)/denominator distance from the start.
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
        # NOTE: theoretically the start of the new interval could be varied from the original start
        # but this is not implemented
        new_length = numerator*increment
        partial = start + new_length
        shrunk_intervals.append([chrom, str(start), str(partial)])
        if padding:
            shrunk_intervals.append([str(int(chrom)+37), str(partial), str(end)])
    return shrunk_intervals
def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def subselect_intervals(intervals, numerator=1, denominator=2, padding=True):
    subselected_intervals = []
    for i in range(len(intervals)):
        if i % denominator <= numerator-1:
            # Add the interval to the subselected list if its index is a multiple of n
            subselected_intervals.append(intervals[i])
        else:
            if padding:
                # Add a modified interval with a 20 base pair shift if padding is enabled
                new_interval = (str(int(intervals[i][0])+23), intervals[i][1], intervals[i][2])
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
    - Mode A files: Every fraction of n intervals is selected, where n is each ratio denominator
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
    ratios_computed = []
    pad_dict = {True: "pad", False: "nopad"}
    intervals = bed_to_intervals(bedfile)
    output_files = []
    for i in ratio_denominators:
        for j in range(i):
            j_over_i = f"{j/i:.2}".rstrip('0').rstrip('.')
            if j_over_i in ratios_computed:
                continue
            ratios_computed.append(j_over_i)
            for pad in [True, False]:
                output_files.append(
                    write_bedfile(
                        shrink_intervals(intervals, numerator=j, denominator=i, padding=pad), 
                        os.path.join(output_dir, f"{output_prefix}_{pad_dict[pad]}_modeB_{j_over_i}.bed")))
                # output_files.append(
                #     write_bedfile(
                #         shrink_intervals(intervals, numerator=j, denominator=i, padding=False), 
                #         os.path.join(output_dir, f"{output_prefix}_{pad_dict[pad]}_modeB_{j_over_i}.bed")))
                output_files.append(
                    write_bedfile(
                        subselect_intervals(intervals, numerator=j, denominator=i, padding=pad), 
                        os.path.join(output_dir, f"{output_prefix}_{pad_dict[pad]}_modeA_{j_over_i}.bed")))
            # output_files.append(
            #     write_bedfile(
            #         subselect_intervals(intervals, numerator=j, denominator=i, padding=False),
            #               os.path.join(output_dir, 
            #                            f"{output_prefix}_nopad_modeA_{j_over_i}.bed")))

    output_files.append(write_bedfile(no_overlap(intervals), 
                                      os.path.join(output_dir, f"{output_prefix}_no_overlap.bed")))
    print_filenames(output_files)
if __name__ == "__main__":
    if len(sys.argv) not in [3, 4]:
        print("Usage: python script.py <input_bedfile> <ratio_denominators> [output_directory]")
        print("Example: python script.py input.bed 2,3,4,5 /path/to/output")
        sys.exit(1)

    # Convert string representation of list back to actual arguments if needed
    if isinstance(sys.argv, str):
        import ast
        try:
            args = ast.literal_eval(sys.argv)
            input_bed = args[1]
            ratio_denominators = [int(x) for x in args[2].split(',')]
            output_dir = args[3] if len(args) == 4 else '.'
        except:
            print(f"Error parsing arguments: {sys.argv}")
            sys.exit(1)
    else:
        input_bed = sys.argv[1]
        ratio_denominators = [int(x) for x in sys.argv[2].split(',')]
        output_dir = sys.argv[3] if len(sys.argv) == 4 else '.'

    # Resolve relative paths
    input_bed = os.path.abspath(os.path.expanduser(input_bed))
    output_dir = os.path.abspath(os.path.expanduser(output_dir))

    if not os.path.exists(input_bed):
        print(f"Error: Input bed file not found: {input_bed}")
        sys.exit(1)

    output_prefix = os.path.basename(os.path.splitext(input_bed)[0])
    
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except OSError as e:
            print(f"Error creating output directory: {e}")
            sys.exit(1)
    
    try:
        generate_variant_bed_files(input_bed, output_prefix, ratio_denominators, output_dir)
    except Exception as e:
        print(f"Error generating variant bed files: {e}")
        sys.exit(1)
