#!/usr/bin/env python
from __future__ import annotations
import sys
import os
from multiprocessing import Pool
import argparse
import resource
import csv
import gc
from typing import Optional, Dict, List, Tuple, Any
from hammock.lib.abstractsketch import AbstractSketch
from hammock.lib.intervals import IntervalSketch
from hammock.lib.sequences import SequenceSketch
from hammock.lib.minimizer import MinimizerSketch
# from hammock.lib.exact import ExactCounter
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

def process_file(args: tuple[str, dict, list, str, int, int, int, int,tuple[float, float], float, str]) -> tuple[Optional[str], Optional[dict]]:
    """Process a single file and calculate similarity values against primary sets.
    
    Args:
        args: Tuple containing:
            filepath: Path to file to process
            primary_sets: Dictionary of primary Sketch objects to compare against
            primary_keys: List of keys for primary_sets
            mode: Mode for comparison (A/B/C/D)
            num_hashes: Number of hashes for MinHash sketching
            precision: Precision for HyperLogLog sketching
            kmer_size: Size of k-mers for sequence sketching
            window_size: Size of sliding window for sequence sketching
            subsample: Subsampling rate for points
            expA: Power of 10 exponent to use to multiply contribution of A-type intervals
            sketch_type: Type of sketching to use
            
    Returns:
        Tuple of (basename, output_dict) where:
            basename: Base filename of processed file
            output_dict: Dictionary mapping primary keys to similarity values
            Returns (None, None) if processing fails
    """
    filepath, primary_sets, primary_keys, mode, num_hashes, precision, kmer_size, window_size, subsample, expA, sketch_type = args
    basename = os.path.basename(filepath)
    if not filepath or not os.path.exists(filepath):
        print(f"Error: Invalid or non-existent file path: '{filepath}'", file=sys.stderr)
        return None, None

    # Choose sketch class based on mode
    if mode == "D":
        sketch_args = {
            'sketch_type': sketch_type,
            'window_size': window_size,
            'kmer_size': kmer_size,
            'precision': precision,
            'num_hashes': num_hashes
        }
        comparator = SequenceSketch.from_file(
            filename=filepath,
            **sketch_args
        )
    else:
        sketch_args = {
            'mode': mode,
            'precision': precision,
            'num_hashes': num_hashes,
            'sketch_type': sketch_type,
            'subsample': subsample,
            'expA': expA
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
        if not isinstance(primary, (IntervalSketch, MinimizerSketch, SequenceSketch)):
            print(f"Error: Invalid type for primary set {primary_keys[i]}: {type(primary)}")
            continue
        sim_values = primary.similarity_values(comparator)
        if sim_values is not None:
            output[primary_keys[i]] = sim_values
    
    return basename, output

def parse_args():
    """Parse command line arguments."""
    arg_parser = argparse.ArgumentParser(
        description="""Calculate pairwise Jaccard similarities between lists of BED or sequence files.
        
        Supported file formats:
        - BED format (.bed): Tab-delimited format with chromosome, start, and end positions
        - BigBed format (.bb): Binary indexed version of BED format
        - Any tab-delimited file with at least 3 columns (chr, start, end) in BED-style format
        - Sequence files (.fa, .fasta, .fna, .ffn, .faa, .frn) - automatically uses mode D
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    arg_parser.add_argument('filepaths_file',
                       help='Text file containing paths to files to be compared')
    arg_parser.add_argument('primary_file',
                       help='Text file containing paths to primary files to compare against')
    
    arg_parser.add_argument('--mode', choices=['A', 'B', 'C', 'D'], default='A',
                       help='''Mode for comparison:
                       A: Compare intervals only (default for BED/BigBed files)
                       B: Compare points only
                       C: Compare both intervals and points
                       D: Compare sequences (auto-detected for sequence files)''')
    
    arg_parser.add_argument('--outprefix', '-o', type=str, default="hammock", help='The output file prefix')
    arg_parser.add_argument("--precision", "-p", type=int, help="Precision for HyperLogLog sketching", default=12)
    arg_parser.add_argument("--num_hashes", "-n", type=int, help="Number of hashes for MinHash sketching", default=128)
    arg_parser.add_argument("--subA", type=float, default=1.0, help="Subsampling rate for intervals (0 to 1)")
    arg_parser.add_argument("--subB", type=float, default=1.0, help="Subsampling rate for points (0 to 1)")
    arg_parser.add_argument( "--expA", type=float, default=0, help="Power of 10 exponent to use to multiply contribution of A-type intervals (only valid for mode C)",
    )
    
    # Add mutually exclusive sketch type flags
    sketch_group = arg_parser.add_mutually_exclusive_group()
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
    sketch_group.add_argument(
        "--minimizer",
        action="store_const",
        const="minimizer",
        dest="sketch_type",
        help="Use Minimizer sketching (default for mode D)"
    )
    
    # Set default sketch type
    arg_parser.set_defaults(sketch_type="hyperloglog")
    
    arg_parser.add_argument(
        "--debug",
        action="store_true",
        help="Output debug information including sparsity comparisons"
    )

    # Add sequence-specific arguments
    arg_parser.add_argument('--kmer_size', "-k",type=int,
                       help='Size of k-mers (0 for whole string mode, defaults to 0 for modes A/B/C, 14 for mode D)')
    arg_parser.add_argument('--window_size', '-w', type=int,
                       help='Size of sliding window (defaults to 0 if kmer_size is 0, otherwise 40)')
    
    # Add threads argument
    arg_parser.add_argument("--threads", type=int, default=None,
                       help="Number of threads to use for parallel processing")
    
    # Add seed argument
    arg_parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for sketching (default: 42)")
    
    # Add verbose argument
    arg_parser.add_argument("--verbose", action="store_true",
                        help="Print verbose output")
    
    args = arg_parser.parse_args()
    
    # Set kmer_size default based on mode
    if args.kmer_size is None:
        args.kmer_size = 14 if args.mode == "D" else 0
    
    # Set window_size default based on kmer_size
    if args.window_size is None:
        args.window_size = 0 if args.kmer_size == 0 else 40
    
    return args

def get_new_prefix(outprefix: str, 
                   sketch_type: str,
                   mode: str,
                   num_hashes: Optional[int] = None, 
                   precision: Optional[int] = None,
                   subA: Optional[float] = None,
                   subB: Optional[float] = None,
                   expA: Optional[float] = None) -> str:
    """Get output prefix with appropriate suffix based on sketch type.
    
    Args:
        outprefix: Base output prefix
        sketch_type: Type of sketch (hyperloglog/minhash/exact)
        mode: Mode (A/B/C/D)
        num_hashes: Number of hashes for MinHash
        precision: Precision for HyperLogLog
        
    Returns:
        String with appropriate suffix for sketch type and mode
    """

    if sketch_type == "minhash":
        outprefix = f"{outprefix}_mh_n{num_hashes}_jacc{mode}"
    elif sketch_type == "hyperloglog":
        outprefix = f"{outprefix}_hll_p{precision}_jacc{mode}"
    elif sketch_type == "minimizer":
        outprefix = f"{outprefix}_mnmzr_jacc{mode}"
    else:
        outprefix = f"{outprefix}_exact_jacc{mode}"
        
    # Add either subA or expA (mutually exclusive)
    if expA is not None and expA > 0:
        outprefix = f"{outprefix}_expA{expA:.2f}"
    elif subA is not None and subA != 1.0:
        outprefix = f"{outprefix}_A{subA:.2f}"
        
    # Add subB if present
    if subB is not None and subB != 1.0:
        outprefix = f"{outprefix}_B{subB:.2f}"
        
    return outprefix

def write_results(results: List[Dict[str, Any]], output: str) -> None:
    """Write comparison results to file.
    
    Args:
        results: List of dictionaries containing comparison results
        output: Output file path
    """
    # Get all possible column names from the similarity dictionaries
    similarity_columns = set()
    for result in results:
        if 'similarity' in result:
            similarity_columns.update(result['similarity'].keys())
    
    # Define column order
    columns = ['query', 'reference']
    columns.extend(sorted(similarity_columns))  # Add all similarity measures
    if 'metadata' in results[0]:
        columns.append('metadata')
    
    # Write results
    with open(output, 'w') as f:
        # Write header
        f.write('\t'.join(columns) + '\n')
        
        # Write data
        for result in results:
            row = [result['query'], result['reference']]
            # Add similarity values, defaulting to 0.0 if not present
            for col in sorted(similarity_columns):
                row.append(str(result['similarity'].get(col, 0.0)))
            if 'metadata' in result:
                row.append(result['metadata'])
            f.write('\t'.join(row) + '\n')

def process_sketches(sketches: List[Tuple[str, AbstractSketch]], 
                    output: str,
                    metadata: Optional[Dict[str, str]] = None) -> None:
    """Process all sketch comparisons and write results.
    
    Args:
        sketches: List of (name, sketch) tuples
        output: Output file path
        metadata: Optional metadata dictionary
    """
    results = []
    
    for i, (query_name, query_sketch) in enumerate(sketches):
        for ref_name, ref_sketch in sketches[i:]:
            if query_name == ref_name:
                continue
                
            # Get similarity measures as dictionary
            similarity = query_sketch.similarity_values(ref_sketch)
            
            # Create result entry
            result = {
                'query': query_name,
                'reference': ref_name,
                'similarity': similarity
            }
            
            # Add metadata if provided
            if metadata and query_name in metadata:
                result['metadata'] = metadata[query_name]
                
            results.append(result)
            
            # Add reverse comparison if needed
            if query_name != ref_name:
                reverse_result = {
                    'query': ref_name,
                    'reference': query_name,
                    'similarity': similarity
                }
                if metadata and ref_name in metadata:
                    reverse_result['metadata'] = metadata[ref_name]
                results.append(reverse_result)
    
    write_results(results, output)

def main():
    """Main entry point for hammock."""
    args = parse_args()
    
    # Read first filepath to check extension
    with open(args.filepaths_file) as f:
        first_file = f.readline().strip()
        if first_file.endswith(('.fa', '.fasta', '.fna', '.ffn', '.faa', '.frn')):
            if args.mode != "D":
                args.mode = "D"
                print(f"Detected sequence file format, switching to mode D")
                
                # Check if any sketch type flag was explicitly provided in command line
                sketch_flags = ["--hyperloglog", "--exact", "--minhash", "--minimizer"]
                sketch_type_explicitly_set = any(flag in sys.argv for flag in sketch_flags)
                
                # Only change default sketch type if not explicitly set by user
                if not sketch_type_explicitly_set:
                    args.sketch_type = "minimizer"
                    print("Changing default sketch type to 'minimizer' for sequence data")
    
    # Check if C-specific parameters are used and switch mode if needed
    c_mode_params_used = args.subA != 1.0 or args.subB != 1.0 or args.expA > 0
    if c_mode_params_used and args.mode == "A":
        print(f"C-mode parameters detected (subA={args.subA}, subB={args.subB}, expA={args.expA}), switching from mode A to mode C")
        args.mode = "C"
    elif args.mode == "A":
        print("Using default mode A for interval comparison")
            
    # Validate expA is only used with mode C
    if args.mode != "C":
        if args.expA > 0:
            raise ValueError("--expA parameter is invalid outside of mode C")
        if args.subA != 1.0 or args.subB != 1.0:
            raise ValueError("Mode C is the only mode that allows subsampling. Please change mode to C or remove --subA and --subB from your command.")
    else:
        if args.expA < 0:
            raise ValueError("--expA parameter must be non-negative for mode C")
        # Validate subsample rates
        if not (0 <= args.subA <= 1 and 0 <= args.subB <= 1):
            raise ValueError("Subsample rates must be between 0 and 1")

    # Package subA and subB into a tuple for processing
    subsample = (args.subA, args.subB)
    
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
            primary_sets[basename] = SequenceSketch.from_file(
                filename=filepath,
                sketch_type=args.sketch_type,
                kmer_size=args.kmer_size,
                window_size=args.window_size,
                precision=args.precision,
                num_hashes=args.num_hashes,
                seed=args.seed,
                verbose=args.verbose
            )
        else:
            sketch_args = {
                'mode': args.mode,
                'precision': args.precision,
                'num_hashes': args.num_hashes,
                'kmer_size': args.kmer_size,
                'window_size': args.window_size,
                'sketch_type': args.sketch_type,
                'subsample': subsample,
                'expA': args.expA,
                'debug': args.debug
            }
            primary_sets[basename] = IntervalSketch.from_file(
                filename=filepath,
                **sketch_args
            )
    
    # Use specified threads or default to os.cpu_count()
    num_threads = args.threads if args.threads is not None else os.cpu_count()
    
    # Process remaining files in parallel
    pool_args = [
        (filepath, primary_sets, primary_keys, args.mode, args.num_hashes, 
         args.precision, args.kmer_size, args.window_size, subsample, args.expA, args.sketch_type)
        for filepath in filepaths
    ]
    
    with Pool(processes=num_threads) as pool:
        results = pool.map(process_file, pool_args)
    
    # Write results
    outprefix = get_new_prefix(args.outprefix, args.sketch_type, args.mode, args.num_hashes, args.precision, args.subA, args.subB, args.expA)
    with open(f"{outprefix}.csv", "w", newline='') as f:
        writer = csv.writer(f)
        
        # Get all similarity measure names from first result
        similarity_measures = []
        for basename, output in results:
            if basename and output and output.get(primary_keys[0]):
                if isinstance(output[primary_keys[0]], dict):
                    similarity_measures = list(output[primary_keys[0]].keys())
                else:
                    similarity_measures = ["jaccard"]
                break
        
        # Write header: file1, file2, sketch_type, mode, [sketch params], [C-mode params], similarity measures
        header = ["file1", "file2", "sketch_type", "mode"]
        if args.sketch_type in ["hyperloglog", "minhash", "minimizer"]:
            header.extend(["precision", "num_hashes", "kmer_size", "window_size"])
        if args.mode == "C":
            header.extend(["subA", "subB", "expA"])
        header.extend(similarity_measures)
        writer.writerow(header)
        
        # Write results
        for basename, output in results:
            if basename and output:
                for key in primary_keys:
                    if key in output:
                        row = [basename, key, args.sketch_type, args.mode]  # files and sketch info
                        if args.sketch_type in ["hyperloglog", "minhash"]:
                            precision = args.precision if args.sketch_type == "hyperloglog" else "NA"
                            num_hashes = args.num_hashes if args.sketch_type == "minhash" else "NA"
                            window_size = "NA" if args.kmer_size == 0 else args.window_size
                            row.extend([precision, num_hashes, args.kmer_size, window_size])
                        if args.mode == "C":
                            row.extend([args.subA, args.subB, args.expA])  # C-mode parameters
                        # Add similarity values
                        if isinstance(output[key], dict):
                            row.extend(output[key].values())
                        else:
                            row.append(output[key])
                        writer.writerow(row)

if __name__ == "__main__":
    main()