#!/usr/bin/env python
import sys
import os
from multiprocessing import Pool
import argparse
import resource
from itertools import islice
import gc
from typing import Optional
from Bio import SeqIO # type: ignore
from hammock.lib.intervals import IntervalSketch
from hammock.lib.minimizer import MinimizerSketch

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

def process_file(args: tuple[str, dict, list, str, int, int, tuple[float, float], str]) -> tuple[Optional[str], Optional[dict]]:
    """Process a single file and calculate similarity values against primary sets.
    
    Args:
        args: Tuple containing:
            filepath: Path to file to process
            primary_sets: Dictionary of primary Sketch objects to compare against
            primary_keys: List of keys for primary_sets
            mode: Mode for comparison (A/B/C/D)
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

    # Choose sketch class based on mode
    if mode == "D":
        sketch_args = {
            'window_size': args.window_size,
            'kmer_size': args.kmer_size,
            'precision': precision,
            'num_hashes': num_hashes
        }
        comparator = MinimizerSketch.from_file(
            filename=filepath,
            **sketch_args
        )
    else:
        sketch_args = {
            'mode': mode,
            'precision': precision,
            'num_hashes': num_hashes,
            'sketch_type': sketch_type,
            'subsample': subsample
        }
        comparator = IntervalSketch.from_file(
            filename=filepath,
            **sketch_args
        )
    
    if comparator is None:
        return None, None

    output = {}
    for i in range(len(primary_keys)):
        primary = primary_sets[primary_keys[i]]
        if not isinstance(primary, (IntervalSketch, MinimizerSketch)):
            print(f"Error: Invalid type for primary set {primary_keys[i]}: {type(primary)}")
            continue
        sim_values = primary.similarity_values(comparator)
        if sim_values is not None:
            output[primary_keys[i]] = sim_values
    
    return basename, output

def get_parser():
    """Create and configure the argument parser for hammock.
    
    Returns:
        argparse.ArgumentParser: Configured parser with all command line arguments
    """ 
    parser = argparse.ArgumentParser(
        description="Calculate similarity between BED or sequence files using sketching",
        epilog="\n\n"
    )
    parser.add_argument("filepaths_file", type=str, help="File containing paths of BED files to compare.")
    parser.add_argument("primary_file", type=str, help="File containing paths of primary comparator BED files. If one of the file lists is shorter, it should be this one.")
    parser.add_argument('--outprefix', '-o', type=str, default="hammock", help='The output file prefix')
    parser.add_argument(
        "--mode",
        type=str,
        required=True,
        choices=["A", "B", "C", "D"],
        help="Mode to indicate comparison type: A=interval, B=point, C=both, D=sequence"
    )
    parser.add_argument("--precision", "-p", type=int, help="Precision for HyperLogLog sketching", default=8)
    parser.add_argument("--num_hashes", "-n", type=int, help="Number of hashes for MinHash sketching", default=128)
    parser.add_argument("--subA", type=float, default=1.0, help="Subsampling rate for intervals (0 to 1)")
    parser.add_argument("--subB", type=float, default=1.0, help="Subsampling rate for points (0 to 1)")
    
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

    # Add sequence-specific arguments
    parser.add_argument(
        "--window_size", "-w",
        type=int,
        default=40,
        help="Window size for sequence minimizers (mode D only)"
    )
    parser.add_argument(
        "--kmer_size", "-k",
        type=int,
        default=8,
        help="K-mer size for sequence minimizers (mode D only)"
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
        return f"{outprefix}_minhash_{num_hashes}"
    elif sketch_type == "hyperloglog":
        return f"{outprefix}_hll_{precision}"
    else:
        return f"{outprefix}_exact"

def main():
    """Main entry point for hammock."""
    parser = get_parser()
    args = parser.parse_args()
    # Package subA and subB into a tuple for processing
    subsample = (args.subA, args.subB)
    # Validate subsample rates
    if not (0 <= subsample[0] <= 1 and 0 <= subsample[1] <= 1):
        raise ValueError("Subsample rates must be between 0 and 1")
        
    
    # Set memory limit
    limit_memory()
    
    # Read file paths
    with open(args.filepaths_file) as f:
        filepaths = [line.strip() for line in f if line.strip()]
    with open(args.primary_file) as f:
        primary_paths = [line.strip() for line in f if line.strip()]
        
    # Process primary files first
    primary_sets = {}
    primary_keys = []
    
    for filepath in primary_paths:
        if not os.path.exists(filepath):
            print(f"Error: Primary file not found: {filepath}", file=sys.stderr)
            continue
            
        basename = os.path.basename(filepath)
        primary_keys.append(basename)
        
        # Choose sketch class based on mode
        if args.mode == "D":
            sketch_args = {
                'window_size': args.window_size,
                'kmer_size': args.kmer_size,
                'precision': args.precision,
                'num_hashes': args.num_hashes
            }
            primary_sets[basename] = MinimizerSketch.from_file(
                filename=filepath,
                **sketch_args
            )
        else:
            sketch_args = {
                'mode': args.mode,
                'precision': args.precision,
                'num_hashes': args.num_hashes,
                'sketch_type': args.sketch_type,
                'subsample': subsample # Pass the tuple instead of separate values
            }
            primary_sets[basename] = IntervalSketch.from_file(
                filename=filepath,
                **sketch_args
            )
    
    # Process remaining files in parallel
    pool_args = [
        (filepath, primary_sets, primary_keys, args.mode, args.num_hashes, 
         args.precision, subsample, args.sketch_type)  # Pass the tuple here
        for filepath in filepaths
    ]
    
    with Pool() as pool:
        results = pool.map(process_file, pool_args)
    
    # Write results
    outprefix = get_new_prefix(args.outprefix, args.sketch_type, args.num_hashes, args.precision)
    with open(f"{outprefix}.tsv", "w") as f:
        # Write header
        f.write("file\t" + "\t".join(primary_keys) + "\n")
        
        # Write results
        for basename, output in results:
            if basename and output:
                f.write(basename)
                for key in primary_keys:
                    if key in output:
                        f.write(f"\t{output[key]}")
                    else:
                        f.write("\tNA")
                f.write("\n")

if __name__ == "__main__":
    main()