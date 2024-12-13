#!/usr/bin/env python
# import random
import sys
import numpy as np # type: ignore
import os
# import xxhash
from multiprocessing import Pool
from hammock.lib.bedprocessing import basic_bedline
import argparse
import resource
from itertools import islice
import gc
from hammock.lib.sketchclass import Sketch
from hammock.lib.hyperloglog import HyperLogLog
# from hammock.lib.exact import ExactCounter
from typing import Optional

# Set memory limit to 28GB (adjust as needed)
def limit_memory():
    """Set memory usage limit for the process.
    
    Sets a soft memory limit of 28GB while preserving the existing hard limit.
    This helps prevent out-of-memory errors by limiting memory allocation.
    
    Uses resource.RLIMIT_AS to set the maximum area (in bytes) of address space 
    that may be taken by the process.
    """
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (28 * 1024 * 1024 * 1024, hard))

def bed_to_sets(filename: str, 
                mode: str, 
                num_hashes: int, 
                precision: int, 
                sep: str = "-", 
                subsample: float = 1, 
                verbose: bool = False, 
                chunk_size: int = 1000, 
                sketch_type: str = "hyperloglog") -> Sketch:
    """Process bed file into a sketch of points or intervals or both. Default sketch type is HyperLogLog.
    
    Args:
        filename: Path to the BED file to process
        mode: A/B/C for interval/point/both comparison types
            A: Compare intervals only
            B: Compare points only 
            C: Compare both intervals and points
        num_hashes: Number of hash functions for MinHash sketching
        precision: Precision parameter for HyperLogLog sketching
        sep: Separator for combining fields. Defaults to "-"
        subsample: Fraction of points (type B) to sample. Defaults to 1. 
                  If negative, interval (type A) values will be subsampled by (1-subsample)
        verbose: Whether to print progress. Defaults to False
        chunk_size: Number of lines to process at once. Defaults to 1000
        sketch_type: Type of sketch to use:
            hyperloglog: HyperLogLog sketch for cardinality estimation
            minhash: MinHash sketch for Jaccard similarity
            exact: Exact set representation
            
    Returns:
        Sketch object containing the processed data
    """
    outset = Sketch(
        sketch_type=sketch_type,
        precision=precision if sketch_type == "hyperloglog" else 14,
        num_hashes=num_hashes if sketch_type == "minhash" else 128,
        kmer_size=0,
        seed=0
    )
    
    def process_chunk(chunk: list[str]) -> None:
        bedline_results = [bedline(line, mode=mode, sep=sep, subsample=subsample) for line in chunk if not line.startswith('#')]
        #NOTE: this would be more efficient if it was in the loop
        outset.total_interval_size += sum(isize for _, _, isize in bedline_results)
        outset.num_intervals += len(bedline_results)
        if mode in ["B", "C"]:
            for _, b, _ in bedline_results:
                if b is not None:
                    for point in b:
                        if point is not None:
                            outset.add_string(point.decode())
        if mode in ["A", "C"]:
            for a, _, _ in bedline_results:
                if a is not None:
                    outset.add_string(a.decode())
        gc.collect()

    with open(filename, "r") as file:
        while True:
            chunk = list(islice(file, chunk_size))
            if not chunk:
                break
            process_chunk(chunk)
            if verbose:
                print(f"Processed {chunk_size} lines from {filename}")
            gc.collect()

    return outset

# def bed_to_sets_parallel(args):
#     """Helper function for parallel processing of bed files into sketches.
    
#     Args:
#         args: Tuple containing:
#             basename: Base name for the output file
#             filepath: Path to the input BED file
#             mode: A/B/C for interval/point/both comparison types 
#             num_hashes: Number of hash functions for MinHash sketching
#             precision: Precision parameter for HyperLogLog sketching
#             subsample: Fraction of points to sample
#             sketch_type: Type of sketch to use (hyperloglog/minhash/exact)
            
#     Returns:
#         Tuple of (basename, sketch) where sketch is the processed Sketch object
#     """
#     basename, filepath, mode, num_hashes, precision, subsample, sketch_type = args
#     return (basename, bed_to_sets(filename=filepath, mode=mode, num_hashes=num_hashes, 
#                                 precision=precision, subsample=subsample, sketch_type=sketch_type))

def bedline(line: str, 
            mode: str, 
            sep: str, 
            subsample: float = 1) -> tuple[Optional[bytes], list[Optional[bytes]]]:
    """Process a single line describing an interval from a BED file. If subsample is provided points will be sampled from the interval.
    
    Args:
        line: Single line from BED file
        mode: A/B/C for interval/point/both comparison types
        sep: Separator for combining fields
        subsample: Subsampling rate. If negative, interval values are subsampled
        
    Returns:
        Tuple of (interval, points) where:
            interval: Encoded interval string or None
            points: List of encoded point strings or empty list
    """
    # NOTE: CHECK THIS: right now we are counting i:i+1 to be the position of the ith bp, meaning chrx 2:4 contains bp 2,3 only. see here https://bedtools.readthedocs.io/en/latest/content/general-usage.html
    interval = ''
    points = []
    chrval, startx, endx = basic_bedline(line)
    start = min(startx, endx)
    end = max(startx, endx)
    interval_size = end - start
    
    if mode in ["A","C"]:
        interval = sep.join([chrval, str(start), str(end), "A"]).encode('utf-8')
        if subsample < 0:
            hashv = HyperLogLog._hash_str(interval, seed=777)
            if hashv % (2**32) > int((1+subsample) * (2**32)): # Modulo to ensure we're in 32-bit space
                interval = None
    if mode in ["B","C"]:
        points = generate_points(chrval, start, end, sep=sep, subsample=abs(subsample))
        
        # for val in range(start, end):
        #     points.append(sep.join([str(chrval), str(val), str(val+1), "B"]).encode('utf-8'))
    return interval, points, interval_size #[g for g in points]

# def parallel_criterion(args):
#     chrval, start, sep, maximum, seed = args
#     outstr = sep.join([str(chrval), str(start), str(start+1)])
#     hashv = HyperLogLog._hash_str(outstr.encode('utf-8'), seed)
#     if hashv <= maximum:
#         return outstr.encode('utf-8')
    
# def generate_points_parallel(chrval, start, end, sep="-", subsample=1, seed=23):
#     args = ((chrval, val, sep, abs(subsample)*sys.maxsize, seed) for val in range(start, end+1))
#     with Pool() as pool:
#         points = pool.map_async(parallel_criterion, args, chunksize=1000)
#     return points


def generate_points(chrval: str, 
                   start: int, 
                   end: int, 
                   sep: str = "-", 
                   subsample: float = 1, 
                   seed: int = 23) -> list[Optional[bytes]]:
    """Generate points from a BED interval. If subsample is provided, points will be sampled from the interval.
    
    Args:
        chrval: Chromosome value
        start: Start position
        end: End position
        sep: Separator for combining fields
        subsample: Fraction of points to sample (0-1)
        seed: Random seed for sampling
        
    Returns:
        List of encoded point strings that passed sampling
    """
    # Convert subsample ratio to the correct threshold
    maximum = int(subsample * (2**32))  # Using 32-bit space for better distribution
    def gp(x):
        outstr = sep.join([str(chrval), str(x), str(x+1)])
        hashv = HyperLogLog._hash_str(outstr.encode('utf-8'), seed)
        if hashv % (2**32) <= maximum:  # Modulo to ensure we're in 32-bit space
            return outstr.encode('utf-8')
    return [gp(x) for x in range(start, end+1)]

def similarity_values(sketch1: Sketch, 
                     sketch2: Sketch) -> tuple[float, float, float, float]:
    """Calculate similarity values between two sketches.
    
    Args:
        sketch1: First sketch
        sketch2: Second sketch
        
    Returns:
        Tuple of (union_size, intersection_size, jaccard, jaccard_calc) where:
            union_size: Estimated size of union
            intersection_size: Estimated size of intersection  
            jaccard: Jaccard similarity from sketch's estimate
            jaccard_calc: Jaccard calculated from union/intersection
    """
    union_size = sketch1.estimate_union(sketch2)
    intersect = sketch1.estimate_intersection(sketch2)
    jaccard = sketch1.estimate_jaccard(sketch2)
    
    # Handle potential division by zero
    jacc_calc = intersect / union_size if union_size != 0 else 0
    
    return union_size, intersect, jaccard, jacc_calc


def process_file(args: tuple[str, dict, list, str, int, int, float, str]) -> tuple[Optional[str], Optional[dict]]:
    """Process a single BED file and calculate similarity values against primary sets.
    
    Args:
        args: Tuple containing:
            filepath: Path to BED file to process
            primary_sets: Dictionary of primary Sketch objects to compare against
            primary_keys: List of keys for primary_sets
            mode: Mode for bed_to_sets ('A', 'B', or 'C')
            num_hashes: Number of hashes for MinHash sketching
            precision: Precision for HyperLogLog sketching
            subsample: Subsampling rate for points
            sketch_type: Type of sketching to use
            
    Returns:
        Tuple of (basename, output_dict) where:
            basename: Base filename of processed file
            output_dict: Dictionary mapping primary keys to similarity values
            Returns (None, None) if processing fails
    """
    filepath, primary_sets, primary_keys, mode, num_hashes, precision, subsample, sketch_type = args
    basename = os.path.basename(filepath)
    if not filepath or not os.path.exists(filepath):
        print(f"Error: Invalid or non-existent file path: '{filepath}'", file=sys.stderr)
        return None, None

    basename = os.path.basename(filepath)
    comparator = bed_to_sets(
        filename=filepath, 
        mode=mode, 
        num_hashes=num_hashes, 
        precision=precision, 
        subsample=subsample,
        sketch_type=sketch_type
    )
    if comparator is None:
        return None, None

    output = {}

    for i in range(len(primary_keys)):
        primary = primary_sets[primary_keys[i]]
        if not isinstance(primary, Sketch):
            print(f"Error: Invalid type for primary set {primary_keys[i]}: {type(primary)}")
            continue
        sim_values = similarity_values(primary, comparator)
        if sim_values is not None:
            output[primary_keys[i]] = sim_values
    
    return basename, output

def get_parser():
    """Create and configure the argument parser for hammock.
    
    Returns:
        argparse.ArgumentParser: Configured parser with all command line arguments
    """ 
    parser = argparse.ArgumentParser(
        description="Calculate similarity between BED files using sketching",
        epilog="\n\n"
    )
    parser.add_argument("filepaths_file", type=str, help="File containing paths of BED files to compare.")
    parser.add_argument("primary_file", type=str, help="File containing paths of primary comparator BED files.")
    parser.add_argument('--outprefix', '-o', type=str, default="hammock", help='The output file prefix')
    parser.add_argument(
        "--mode",
        type=str,
        required=True,
        choices=["A", "B", "C"],
        help="Mode to indicate comparison type: A=interval, B=point, C=both"
    )
    parser.add_argument("--precision", "-p", type=int, help="Precision for HyperLogLog sketching", default=8)
    parser.add_argument("--num_hashes", "-n", type=int, help="Number of hashes for MinHash sketching", default=128)
    parser.add_argument("--subsample", type=float, default=1, help="Subsampling rate for points: provide decimal ratio.")
    parser.add_argument("--balance", action="store_true", help="When set, type A values will be subsampled by (1-subsample)")
    
    # Add mutually exclusive sketch type flags
    sketch_group = parser.add_mutually_exclusive_group()
    sketch_group.add_argument(
        "--hyperloglog",
        action="store_const",
        const="hyperloglog",
        dest="sketch_type",
        help="Use HyperLogLog sketching (default)"
    )
    sketch_group.add_argument(
        "--exact",
        action="store_const",
        const="exact",
        dest="sketch_type",
        help="Use exact counting"
    )
    sketch_group.add_argument(
        "--minhash",
        action="store_const",
        const="minhash",
        dest="sketch_type",
        help="Use MinHash sketching"
    )
    
    # Set default sketch type
    parser.set_defaults(sketch_type="hyperloglog")
    
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Output debug information including sparsity comparisons"
    )
    
    return parser

def get_new_prefix(outprefix: str, 
                   sketch_type: str, 
                   num_hashes: Optional[int] = None, 
                   precision: Optional[int] = None) -> str:
    """Get output prefix with appropriate suffix based on sketch type.
    
    Args:
        outprefix: Base output prefix
        sketch_type: Type of sketch (hyperloglog/minhash/exact)
        num_hashes: Number of hashes for MinHash
        precision: Precision for HyperLogLog
        
    Returns:
        String with appropriate suffix for sketch type
    """
    if sketch_type == "minhash":
        suffix = f"_n{num_hashes}_mh"
    elif sketch_type == "hyperloglog":
        suffix = f"_p{precision}_hll"
    elif sketch_type == "exact":
        suffix = "_exact"
    return f"{outprefix}{suffix}"

def main():
    args = get_parser().parse_args()
    
    if args.balance and args.mode != 'C':
        print(f"\nError: --balance flag is only valid with mode C.")
        print("This flag controls subsampling between interval (A) and point (B) comparisons.")
        print(f"It has no effect in mode {args.mode} which only uses {'intervals' if args.mode == 'A' else 'points'}.")
        sys.exit(1)
    
    limit_memory()
    subsample = args.subsample
    mode = args.mode
    outprefix = args.outprefix
    if os.path.isdir(outprefix):
        outprefix = os.path.join(os.path.dirname(outprefix), os.path.basename(outprefix))

    # If balance is set, make subsample negative to indicate that type A values should be subsampled by (1-subsample)
    if args.balance:
        subsample = - subsample

    bed_sets = {}
    filepaths_file = args.filepaths_file
    pfile = args.primary_file
    # pfile = "/home/jbonnie1/interval_sketch/hammock/cobind_repro/data/TFbeds_primary.txt"
    num_hashes = args.num_hashes
    precision = args.precision
    with open(filepaths_file, "r") as filein:
        filepaths = [f.strip() for f in filein.readlines()]
    
    prime_sets = {}
    with open(pfile, "r") as filein:
        for line in filein:
            filepath = line.strip()
            basename = os.path.basename(filepath)
            prime_sets[basename] = bed_to_sets(
                filename=filepath, 
                mode=mode, 
                num_hashes=num_hashes, 
                precision=precision, 
                subsample=subsample,
                sketch_type=args.sketch_type
            )
    
    prime_keys = list(prime_sets.keys())
    prime_keys.sort()

    jacc_long = [",".join(["bed1", "bed2", "sketch_type", "mode", "num_hashes", "precision", "subsample","balance", "jaccardfunc", "jaccardcalc","intersect", "union" ])]
    # Create a list of arguments for each worker
    args_list = [(f, prime_sets, prime_keys, args.mode, args.num_hashes, 
                  args.precision, subsample, args.sketch_type) for f in filepaths]
    results = []
    # with Pool() as pool:
        # Map the worker function to the arguments list and get the results
        # results = pool.map(process_file, args_list)
    for argy in args_list:
        results.append(process_file(argy))

    result_dict = {b:out for (b, out) in results if b is not None and out is not None}
    result_keys = list(result_dict.keys())
    result_keys.sort()
    for j in range(len(result_keys)):
        secset = result_keys[j]
        for i in range(len(prime_keys)):
            primeset = prime_keys[i]
            union, intersect, jaccard, jacc_calc = result_dict[secset][primeset]
            if args.debug:
                print(f"\nComparing {primeset} with {secset}:")
                print(f"Set 1 size: {prime_sets[primeset].get_size()}")
                print(f"Set 2 size: {result_dict[secset][primeset][0]}")  # union size
                print(f"Sparsity ratio: {min(prime_sets[primeset].get_size(), result_dict[secset][primeset][0]) / max(prime_sets[primeset].get_size(), result_dict[secset][primeset][0]):.3f}")
                print(f"Jaccard similarity: {jaccard:.3f}")
            
            if primeset == secset:
                secset = "-"
            prec_str = "-" if args.sketch_type in ["exact", "minhash"] else str(precision)
            hash_str = "-" if args.sketch_type == "exact" else str(num_hashes)
            jacc_long.append(",".join(
                [primeset,
                 secset,
                 args.sketch_type,
                 mode,
                 hash_str,
                 prec_str,
                 str(args.subsample),
                 str(args.balance),
                 f"{jaccard:>5f}",
                 f"{jacc_calc:>5f}",
                 f"{intersect:>5f}",
                 f"{union:>5f}"]))
    new_prefix = get_new_prefix(outprefix, args.sketch_type, num_hashes, precision)
    if mode == "C":
        # Recall that subsample is negative if balance is set
        new_prefix = f"{new_prefix}_sub{subsample}"
    # with open(f"{new_prefix}_card{mode}.csv", mode="w") as outfile:
    #     for line in outlines:
    #         outfile.write(line + "\n")
    with open(f"{new_prefix}_jacc{mode}.csv", mode="w") as outfile:
        for line in jacc_long:
            outfile.write(line + "\n")

    limit_memory()
    return 0

if __name__ == "__main__":
    main()