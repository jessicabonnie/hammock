#!/usr/bin/env python
# import random
import sys
import numpy as np
import os
import xxhash
from multiprocessing import Pool
from datasketch import MinHash, LeanMinHash
from bedprocessing import basic_bedline
import argparse
import resource
from itertools import islice
import gc

# Set memory limit to 28GB (adjust as needed)
def limit_memory():
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (28 * 1024 * 1024 * 1024, hard))

# class HashXX32(object):
#     def __init__(self, seed):
#         self.h = xxhash.xxh32(seed=seed)

#     def hash(self, o):
#         self.h.reset()
#         self.h.update(o)
#         return self.h.intdigest() % sys.maxsize
def _hash_64(b, seed=0):
    return xxhash.xxh64(b, seed=seed).intdigest() % sys.maxsize

def _hash_32(b, seed=0):
    return xxhash.xxh32_intdigest(b, seed=0) % sys.maxsize

def hashfunc(x):
    if isinstance(x, list):
        # If x is a list, hash each element separately and combine the results
        return xxhash.xxh32(b''.join(str(item).encode('utf-8') for item in x)).intdigest()
    else:
        # If x is not a list, proceed as before
        return xxhash.xxh32(str(x).encode('utf-8')).intdigest()

def bed_to_sets(filename, mode, parts, sep="-", subsample=1, verbose=False, chunk_size=1000):
    outset = MinHash(num_perm=int(parts), seed=0, hashfunc=hashfunc)
    
    def process_chunk(chunk):
        bedline_results = [bedline(line, mode=mode, sep=sep, subsample=subsample) for line in chunk if not line.startswith('#')]
        if mode in ["B", "C"]:
            outset.update_batch( [b for (_, b) in bedline_results if b is not None] )
        if mode in ["A", "C"]:
            outset.update_batch([a for (a, _) in bedline_results if a != ''])
        # del bedline_results
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

    return LeanMinHash(outset)

def bed_to_sets_parallel(args):
    basename, filepath, mode, parts, subsample = args
    return (basename, bed_to_sets(filename=filepath, mode=mode, parts=parts, subsample=subsample))

def bedline(line, mode, sep, subsample=1):
    # NOTE: right now we are counting i:i+1 to be the position of the ith bp, meaning chrx 2:4 contains bp 2,3 only. see here https://bedtools.readthedocs.io/en/latest/content/general-usage.html
    interval = ''
    points = []
    chrval, startx, endx = basic_bedline(line)
    start = min(startx, endx)
    end = max(startx, endx)
    if mode in ["A","C"]:
        interval = sep.join([chrval, str(start), str(end), "A"]).encode('utf-8')
    if mode in ["B","C"]:
        points = generate_points(chrval, start, end, sep=sep, subsample=subsample)
        
        # for val in range(start, end):
        #     points.append(sep.join([str(chrval), str(val), str(val+1), "B"]).encode('utf-8'))
    return interval, points #[g for g in points]

def parallel_criterion(args):
    chrval, start, sep, maximum, seed = args
    outstr = sep.join([str(chrval), str(start), str(start+1)])
    hashv = _hash_32(outstr.encode('utf-8'), seed=seed)
    if hashv <= maximum:
        return outstr.encode('utf-8')
    
def generate_points_parallel(chrval, start, end, sep="-", subsample=1, seed=23):
    args = ((chrval, val, sep, subsample*sys.maxsize, seed) for val in range(start, end+1))
    with Pool() as pool:
        points = pool.map_async(parallel_criterion, args, chunksize=1000)
    return points


def generate_points(chrval, start, end, sep="-", subsample=1,seed=23):
    '''
    Generate points from a BED interval.
    '''
    maximum = subsample * sys.maxsize
    def gp(x):
        outstr = sep.join([str(chrval), str(x), str(x+1)])
        hashv = _hash_32(outstr.encode('utf-8'), seed=seed)
        if hashv <= maximum:
            return outstr.encode('utf-8')
    # if not generator:
    #     points = []
    # points = []
    # for val in range(start, end):
    #     points.append(gp(val))
        # points.append(outstr.encode('utf-8'))
    # return points
    return [gp(x) for x in range(start, end+1)]
    


def similarity_values(set1, set2):
    union = LeanMinHash.union(set1, set2)
    unionc = union.count()
    intersect = abs(set1.count() + set2.count() - unionc)
    jaccard = set1.jaccard(set2)
    
    # Handle potential division by zero
    jacc_calc = intersect / unionc if unionc != 0 else 0
    
    return unionc, intersect, jaccard, jacc_calc


def process_file(args):
    filepath, primary_sets, primary_keys, mode, parts, subsample = args
    basename = os.path.basename(filepath)
    if not filepath or not os.path.exists(filepath):
        print(f"Error: Invalid or non-existent file path: '{filepath}'", file=sys.stderr)
        return None, None

    basename = os.path.basename(filepath)
    comparator = bed_to_sets(filename=filepath, mode=mode, parts=parts, subsample=subsample)
    if comparator is None:
        return None, None

    output = {}

    for i in range(len(primary_keys)):
        primary = primary_sets[primary_keys[i]]
        if not isinstance(primary, LeanMinHash):
            print(f"Error: Invalid type for primary set {primary_keys[i]}: {type(primary)}")
            continue
        sim_values = similarity_values(primary, comparator)
        if sim_values is not None:
            output[primary_keys[i]] = sim_values
    
    return basename, output


def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate Jaccard similarity between BED files.")
    parser.add_argument("filepaths_file", type=str, help="File containing paths of BED files to compare.")
    parser.add_argument("primary_file", type=str, help="File containing paths of primary comparator BED files.")
    parser.add_argument('--output', '-o', type=str, default="hammock", help='The output file prefix')
    parser.add_argument("--mode", type=str, required=True, help="Mode to indicate what kind of similarity comparison is desired.")
    parser.add_argument("--perm", "-p", type=int, help="Number of permutations for MinHash.", required=True)
    parser.add_argument("--subsample", type=float, default=1, help="Subsampling rate for points: provide decimal ratio.")
    return parser.parse_args()


if __name__ == "__main__":
    limit_memory()
    args = parse_arguments()
    subsample = args.subsample
    mode = args.mode
    outprefix = args.output

    bed_sets = {}
    filepaths_file = args.filepaths_file
    pfile = args.primary_file
    # pfile = "/home/jbonnie1/interval_sketch/hammock/cobind_repro/data/TFbeds_primary.txt"
    parts = args.perm
    with open(filepaths_file, "r") as filein:
        filepaths = [f.strip() for f in filein.readlines()]
    
    prime_sets = {}
    with open(pfile, "r") as filein:
        for line in filein:
            filepath = line.strip()
            basename = os.path.basename(filepath)
            prime_sets[basename] = bed_to_sets(filename=filepath, mode=mode, parts=parts, subsample=subsample)
    
    prime_keys = list(prime_sets.keys())
    prime_keys.sort()

    # outlines = [",".join(["bed1", "bed2", "mode", "permutations", "subsample", "union", "intersect_inexact"])]
    jacc_long = [",".join(["bed1", "bed2", "mode", "permutations", "subsample", "jaccardfunc", "jaccardcalc","intersect_inexact", "union" ])]
    # Create a list of arguments for each worker
    args_list = [(f, prime_sets, prime_keys, mode, parts, subsample) for f in filepaths]
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
            # primename = prime_keys[i]
            primeset = prime_keys[i]
            
            union, intersect, jaccard, jacc_calc = result_dict[secset][primeset]
            # calc_jacc = intersect/union
            if primeset == secset:
                secset = "-"

            # outlines.append(",".join(
            #     [primeset,
            #      secset,
            #      mode,
            #      str(parts),
            #      str(subsample),
            #      f"{intersect:>5f}",
            #      f"{union:>5f}"]))
            jacc_long.append(",".join(
                [primeset,
                 secset,
                 mode,
                 str(parts),
                 str(subsample),
                 f"{jaccard:>5f}",
                 f"{jacc_calc:>5f}",
                 f"{intersect:>5f}",
                 f"{union:>5f}"]))
    new_prefix = f"{outprefix}_p{parts}"
    if mode == "C":
        new_prefix = f"{new_prefix}_sub{subsample}"
    # with open(f"{new_prefix}_card{mode}.csv", mode="w") as outfile:
    #     for line in outlines:
    #         outfile.write(line + "\n")
    with open(f"{new_prefix}_jacc{mode}.csv", mode="w") as outfile:
        for line in jacc_long:
            outfile.write(line + "\n")

    limit_memory()
