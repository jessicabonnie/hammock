#!/usr/bin/env python
from __future__ import annotations
import sys
import os
from multiprocessing import Pool, cpu_count
import argparse
import resource
import csv
import gc
from typing import Optional, Dict, List, Tuple, Any
from hammock.lib.abstractsketch import AbstractSketch
from hammock.lib.intervals import IntervalSketch
from hammock.lib.sequences import SequenceSketch
# from hammock.lib.minimizer import MinimizerSketch
# from hammock.lib.exact import ExactCounter

# Maximum precision for HyperLogLog sketches
MAX_PRECISION = 24

# Set memory limit to 28GB (adjust as needed)
def limit_memory():
    """Set memory usage limit for the process.
    
    Sets a soft memory limit of 28GB while preserving the existing hard limit.
    This helps prevent out-of-memory errors by limiting memory allocation.
    
    Uses resource.RLIMIT_AS to set the maximum area (in bytes) of address space 
    that may be taken by the process.
    """
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    # Set soft limit to 28GB, preserve hard limit
    resource.setrlimit(resource.RLIMIT_AS, (28 * 1024 * 1024 * 1024, hard))
    print(f"Set memory limit to 28GB (soft limit)")

def process_file(filepath: str, primary_files: List[str], mode: str = 'A',
                num_hashes: int = 64, precision: int = 12, kmer_size: int = 0,
                window_size: int = 0, subA: float = 1.0, subB: float = 1.0,
                expA: float = 0.5, use_rust: bool = False, sketch_type: str = "hyperloglog",
                hash_size: int = 32) -> Dict[str, float]:
    """Process a file and calculate similarity against primary sets."""
    # Set memory limit before processing
    limit_memory()
    
    if not os.path.exists(filepath):
        print(f"Error: File {filepath} does not exist", file=sys.stderr)
        sys.exit(2)
        
    # Print implementation info
    if sketch_type == "hyperloglog":
        if use_rust:
            print(f"Using Rust implementation with {hash_size}-bit hashing")
        else:
            print(f"Using Python implementation with {hash_size}-bit hashing")
    elif sketch_type == "minhash":
        print("Using MinHash implementation")
    elif sketch_type == "minimizer":
        print("Using Minimizer implementation")
    elif sketch_type == "exact":
        print("Using Exact implementation")
        
    # Force kmer_size to 0 for mode B
    if mode == "B":
        kmer_size = 0
        window_size = 0
        
    # Create primary sketches
    primary_sets = {}
    for primary_file in primary_files:
        if not os.path.exists(primary_file):
            print(f"Error: Primary file {primary_file} does not exist", file=sys.stderr)
            sys.exit(2)
            
        # Choose sketch class based on mode
        if mode == "D":
            sketch = SequenceSketch.from_file(
                filename=primary_file,
                sketch_type=sketch_type,
                kmer_size=kmer_size,
                window_size=window_size,
                precision=precision,
                num_hashes=num_hashes
            )
        else:
            sketch = IntervalSketch.from_file(
                filename=primary_file,
                mode=mode,
                precision=precision,
                num_hashes=num_hashes,
                kmer_size=kmer_size,
                window_size=window_size,
                sketch_type=sketch_type,
                subsample=(subA, subB),
                expA=expA,
                use_rust=use_rust
            )
                        
        primary_sets[os.path.basename(primary_file)] = sketch
    
    # Check if the current file is already in primary_sets to avoid recreating sketch
    current_basename = os.path.basename(filepath)
    if current_basename in primary_sets:
        sketch = primary_sets[current_basename]
    else:
        # Create comparator sketch from file
        if mode == "D":
            sketch = SequenceSketch.from_file(
                filename=filepath,
                sketch_type=sketch_type,
                kmer_size=kmer_size,
                window_size=window_size,
                precision=precision,
                num_hashes=num_hashes
            )
        else:
            sketch = IntervalSketch.from_file(
                filename=filepath,
                mode=mode,
                precision=precision,
                num_hashes=num_hashes,
                kmer_size=kmer_size,
                window_size=window_size,
                sketch_type=sketch_type,
                subsample=(subA, subB),
                expA=expA,
                use_rust=use_rust
            )
                    
    # Calculate similarity values
    similarity_values = {}
    for primary_name, primary_sketch in primary_sets.items():
        similarity_values[primary_name] = sketch.similarity_values(primary_sketch)
            
    return similarity_values

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
    
    arg_parser.add_argument('--outprefix', '-o','--out', type=str, default="hammock", help='The output file prefix')
    arg_parser.add_argument("--precision", "-p", type=int, help="Precision for HyperLogLog sketching", default=18)
    arg_parser.add_argument("--num_hashes", "-n", type=int, help="Number of hashes for MinHash sketching", default=128)
    arg_parser.add_argument("--subA", type=float, default=1.0, help="Subsampling rate for intervals (0 to 1)")
    arg_parser.add_argument("--subB", type=float, default=1.0, help="Subsampling rate for points (0 to 1)")
    arg_parser.add_argument("--expA", type=float, default=0, help="Power of 10 exponent to multiply contribution of A-type intervals")
    arg_parser.add_argument("--threads", type=int, help="Number of threads to use", default=None)
    arg_parser.add_argument("--kmer_size", '-k', type=int, default=8, help="Size of k-mers for sequence sketching")
    arg_parser.add_argument("--window_size", '-w', '--window',type=int, default=40, help="Size of sliding window for sequence sketching")
    arg_parser.add_argument("--seed", type=int, default=42, help="Random seed for hashing")
    arg_parser.add_argument("--verbose", action="store_true", help="Print verbose output")
    arg_parser.add_argument("--debug", action="store_true", help="Enable debug mode")
    arg_parser.add_argument("--rust", action="store_true", help="Use Rust implementation for HyperLogLog when available")
    
    # Sketch type arguments
    sketch_group = arg_parser.add_mutually_exclusive_group()
    sketch_group.add_argument("--hyperloglog", action="store_true", help="Use HyperLogLog sketching")
    sketch_group.add_argument("--exact", action="store_true", help="Use exact counting")
    sketch_group.add_argument("--minhash", action="store_true", help="Use MinHash sketching")
    sketch_group.add_argument("--minimizer", action="store_true", help="Use minimizer sketching")
    
    arg_parser.add_argument('--hashsize', type=int, default=32, choices=[32, 64], help='Hash size in bits (32 or 64, default: 32)', dest='hash_size')
    
    return arg_parser.parse_args()

def get_new_prefix(outprefix: str, 
                   sketch_type: str,
                   mode: str,
                   num_hashes: Optional[int] = None, 
                   precision: Optional[int] = None,
                   subA: Optional[float] = None,
                   subB: Optional[float] = None,
                   expA: Optional[float] = None,
                   kmer_size: Optional[int] = None,
                   window_size: Optional[int] = None) -> str:
    """Get output prefix with appropriate suffix based on sketch type.
    
    Args:
        outprefix: Base output prefix
        sketch_type: Type of sketch (hyperloglog/minhash/exact)
        mode: Mode (A/B/C/D)
        num_hashes: Number of hashes for MinHash
        precision: Precision for HyperLogLog
        subA: Subsampling rate for intervals
        subB: Subsampling rate for points
        expA: Power of 10 exponent for A-type intervals
        kmer_size: Size of k-mers for sequence sketching
        window_size: Size of sliding window for sequence sketching
        
    Returns:
        String with appropriate suffix for sketch type and mode
    """

    if sketch_type == "minhash":
        outprefix = f"{outprefix}_mh_n{num_hashes}_jacc{mode}"
    elif sketch_type == "hyperloglog":
        outprefix = f"{outprefix}_hll_p{precision}_jacc{mode}"
    elif sketch_type == "minimizer":
        outprefix = f"{outprefix}_mnmzr_p{precision}_jacc{mode}"
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
        
    # Add kmer_size and window_size for mode D
    if mode == "D" and kmer_size is not None and window_size is not None:
        outprefix = f"{outprefix}_k{kmer_size}_w{window_size}"
        
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

def process_files(files: List[str], **kwargs) -> Dict[str, IntervalSketch]:
    """Process multiple files in parallel."""
    mode = kwargs.get('mode', 'A')
    sep = kwargs.get('sep', '-')
    subsample = kwargs.get('subsample', (1.0, 1.0))
    expA = kwargs.get('expA', 0)
    debug = kwargs.get('debug', False)
    use_rust = kwargs.get('use_rust', False)
    
    # If using Rust, be more conservative with Python-level threading
    if use_rust:
        num_threads = 1  # Let Rust handle the threading
    else:
        num_threads = min(cpu_count(), 4)  # Limit Python threads
    
    if debug:
        print(f"Processing {len(files)} files with mode={mode}, subsample={subsample}, debug={debug}")
        print(f"Using {num_threads} Python threads (Rust threading: {use_rust})")
    
    # Create a sketch for each file
    sketches = {}
    for filename in files:
        try:
            sketch = IntervalSketch.from_bed(
                filename,
                mode=mode,
                sep=sep,
                subsample=subsample,
                expA=expA,
                debug=debug,
                use_rust=use_rust,
                num_threads=num_threads
            )
            if sketch is not None:
                sketches[filename] = sketch
        except Exception as e:
            if debug:
                print(f"Error processing {filename}: {str(e)}")
            continue
    
    return sketches

def main():
    """Main entry point for hammock."""
    import warnings
    
    # Parse arguments
    args = parse_args()
    
    # Set memory limit
    limit_memory()
    
    # Validate precision value
    if args.precision < 4 or args.precision >= args.hash_size:
        warnings.warn(f"Precision {args.precision} is outside recommended range (4-{args.hash_size-10}). "
                     "This may reduce accuracy.", RuntimeWarning)
    
    # Determine sketch type
    if args.hyperloglog:
        args.sketch_type = "hyperloglog"
    elif args.minhash:
        args.sketch_type = "minhash"
    elif args.minimizer:
        args.sketch_type = "minimizer"
    elif args.exact:
        args.sketch_type = "exact"
    else:
        # Default to hyperloglog if no sketch type specified
        args.sketch_type = "hyperloglog"
    
    # Print implementation info first
    if args.sketch_type == "hyperloglog":
        if args.rust:
            print("Using Rust implementation")
        else:
            print("Using Python implementation")
    elif args.sketch_type == "minhash":
        print("Using MinHash implementation")
    elif args.sketch_type == "minimizer":
        print("Using Minimizer implementation")
    elif args.sketch_type == "exact":
        print("Using Exact implementation")
    
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
    if args.mode != "C" and args.expA > 0:
        print(f"C-mode parameters detected (subA={args.subA}, subB={args.subB}, expA={args.expA}), switching from mode A to mode C")
        args.mode = "C"

    # Validate subsample rates
    if args.subA < 0 or args.subA > 1 or args.subB < 0 or args.subB > 1:
        raise ValueError("Subsample rates must be between 0 and 1")

    # Validate expA is non-negative for mode C
    if args.mode == "C" and args.expA < 0:
        raise ValueError("--expA parameter must be non-negative for mode C")

    # Package subA and subB into a tuple for processing
    subsample = (args.subA, args.subB)
    
    # Read file paths
    with open(args.filepaths_file) as f:
        filepaths = [line.strip() for line in f if line.strip()]
    with open(args.primary_file) as f:
        primary_paths = [line.strip() for line in f if line.strip()]
    
    # Check if both input files are the same
    same_input_lists = False
    if args.filepaths_file == args.primary_file:
        same_input_lists = True
        print("Detected same input file for both arguments - sketches will be computed only once")
    else:
        # Check if the contents match
        if set(filepaths) == set(primary_paths):
            same_input_lists = True
            print("Detected identical file lists - sketches will be computed only once")
        
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
                'debug': args.debug,
                'use_rust': args.rust  # Add use_rust flag to sketch arguments
            }
            primary_sets[basename] = IntervalSketch.from_file(
                filename=filepath,
                **sketch_args
            )
    
    # Use specified threads or default to os.cpu_count()
    num_threads = args.threads if args.threads is not None else os.cpu_count()
    
    # Process remaining files in parallel
    pool_args = []
    
    # If using the same input lists, compute pairwise comparisons more efficiently
    if same_input_lists:
        print("Computing all pairwise comparisons directly")
        
        # Create a cache to store already computed similarities
        similarity_cache = {}
        
        # We already have all the sketches in primary_sets
        # For each file, we need to compute the similarity with every other file
        results = []
        for filepath in filepaths:
            basename = os.path.basename(filepath)
            if basename in primary_sets:
                # Get the sketch for this file
                sketch = primary_sets[basename]
                
                # Calculate similarity with all primary files (including self)
                similarity_values = {}
                for primary_name, primary_sketch in primary_sets.items():
                    # Create a cache key that's consistent regardless of order
                    cache_key = tuple(sorted([basename, primary_name]))
                    
                    # Check if we've already computed this similarity
                    if cache_key in similarity_cache:
                        similarity_values[primary_name] = similarity_cache[cache_key]
                    else:
                        # Compute the similarity and cache it, including for self-comparison
                        similarity_result = sketch.similarity_values(primary_sketch)
                        similarity_values[primary_name] = similarity_result
                        similarity_cache[cache_key] = similarity_result
                
                # Add to results
                results.append((basename, similarity_values))
    else:
        # Standard processing for different file lists
        pool_args = [
            (filepath, primary_paths, args.mode, args.num_hashes, args.precision, args.kmer_size, args.window_size, args.subA, args.subB, args.expA, args.rust, args.sketch_type, args.hash_size)
            for filepath in filepaths
        ]
        
        # Process files sequentially instead of in parallel to avoid pickling issues
        results = []
        for pool_arg in pool_args:
            filepath, primary_paths, mode, num_hashes, precision, kmer_size, window_size, subA, subB, expA, use_rust, sketch_type, hash_size = pool_arg
            result = process_file(filepath, primary_paths, mode, num_hashes, precision, kmer_size, window_size, subA, subB, expA, use_rust, sketch_type, hash_size)
            results.append((os.path.basename(filepath), result))
    
    # Write results
    outprefix = get_new_prefix(args.outprefix, args.sketch_type, args.mode, args.num_hashes, args.precision, args.subA, args.subB, args.expA, args.kmer_size, args.window_size)
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
                        if args.sketch_type in ["hyperloglog", "minhash", "minimizer"]:
                            precision = args.precision if args.sketch_type in ["hyperloglog", "minimizer"] else "NA"
                            num_hashes = args.num_hashes if args.sketch_type == "minhash" else "NA"
                            # Only include kmer size and window size if they are relevant
                            if args.mode == "D" or args.kmer_size > 0:
                                row.extend([precision, num_hashes, args.kmer_size, args.window_size])
                            else:
                                row.extend([precision, num_hashes, "NA", "NA"])
                        if args.mode == "C":
                            row.extend([args.subA, args.subB, args.expA])  # C-mode parameters
                        # Add similarity values in the same order as similarity_measures
                        if isinstance(output[key], dict):
                            # Ensure order matches similarity_measures order
                            for measure in similarity_measures:
                                row.append(output[key].get(measure, 0.0))
                        else:
                            row.append(output[key])
                        writer.writerow(row)

if __name__ == "__main__":
    main()