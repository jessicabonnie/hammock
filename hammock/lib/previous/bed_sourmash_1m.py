#!/usr/bin/env python
# import random
import sys
import numpy as np
import os
import xxhash
from multiprocessing import Pool
from datasketch import MinHash
import sourmash
from bedprocessing import basic_bedline
import argparse

def _hash_32(b, seed=0, **kwargs):
    return xxhash.xxh32_intdigest(b, seed=0, **kwargs) % sys.maxsize

# def subsample_points(points, subsample=1, seed=0):
#     '''
#     Subsample a list of points.
#     '''
#     if subsample == 1:
#         return points
#     else:
#         return [ b for b in points if _hash_64(b,seed=seed) % 10 <= subsample*10]
    
# class ExtendMinHash(MinHash):
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
    
#     def cardinality(self):
#         return self.count()
    
#     def union(self, b):
#         """
#         Merges this sketch with one or more other sketches.
#         """
        
#         return( MinHash.union(self,b))
#         # for other in others:
#         #     if not isinstance(other, HyperLogLog):
#         #         raise TypeError(f"Cannot merge ExtendedHyperLogLog with {type(other)}")
#             # self.merge(other)
#         # return self
    
#     def size_union(self, others):
#         """
#         Returns the cardinality of the union of this MinHash sketch
#         with one or more other MinHash sketches.
#         """
        
#         return self.union(others).count()
    
#     def size_intersection(self, b):
#         '''Return the estimated size of the intersection
#             intersection(a,b)=cardinality(a)+cardinality(b)-union(a,b) '''
#         # return self.intersect(b).count()
#         return abs(self.count() + b.count() - self.size_union(b))


def bed_to_sets(filename, mode, parts, sep="-", subsample=1, verbose=False):
    """
    Convert a BED file to sets of intervals and points.

    Args:
        filename (str): The path to the BED file.
        mode (str, optional): The mode to determine how to parse the BED file. Defaults to "A".
        sep (str, optional): The separator used to split the BED file line. Defaults to "-".
        verbose (bool, optional): Whether to print additional information during conversion. Defaults to False.

    Returns:
        list: A list containing two sets - Aset (intervals) and Bset (points).
    """
    outset = sourmash.MinHash(n=parts, ksize=21, is_protein=False, dayhoff=False, seed=0 )
    
    with open(filename, "r") as file:
        for line in file:
            if not line.startswith('#'):
                interval, points = bedline(line, mode=mode, sep=sep, subsample=subsample, returnhash=True)
                if interval == '':
                    break
                if mode in ["A", "C"]:
                    outset.add_hash(interval)
                    # outset.update(interval)
                if mode in ["B", "C"]:
                    outset.add_many(points)
                    # outset.update_batch(points)
                # if mode in ["C"]:
                #     outset.update_batch(points, subsample=subsample)
                    # outset.update_batch(subsample_points(points, subsample))
    print(f"finished bed_to_sets for {filename}")
    return outset

def bedline(line, mode, sep, subsample=1, returnhash=False):
    # NOTE: right now we are counting i:i+1 to be the position of the ith bp, meaning chrx 2:4 contains bp 2,3 only. see here https://bedtools.readthedocs.io/en/latest/content/general-usage.html
    interval = ''
    points = []
    chrval, start, end = basic_bedline(line)
    if mode in ["A","C"]:
        interval = sep.join([chrval, str(start), str(end), "A"]).encode('utf-8')
    if mode in ["B","C"]:
        points = generate_points(chrval, start, end, sep=sep, subsample=subsample, returnhash=returnhash, generator=False)    
        # for val in range(start, end):
        #     points.append(sep.join([str(chrval), str(val), str(val+1), "B"]).encode('utf-8'))
    if returnhash:
        interval = _hash_32(interval)
        # interval = sourmash.minhash.hash_murmur(interval, 0)
    return interval, [g for g in points]

def parallel_criterion(args):
    chrval, start, sep, maximum, seed = args
    outstr = sep.join([str(chrval), str(start), str(start+1)])
    hashv = _hash_32(outstr.encode('utf-8'), seed=seed)
    if hashv <= maximum:
        return outstr.encode('utf-8')

def generate_points(chrval, start, end, sep="-", subsample=1, returnhash=False, seed=0, generator=True):
    '''
    Generate points from a BED interval.
    '''
    maximum = sys.maxsize * subsample
    if not generator:
        points = []
    for val in range(start, end):
        outstr = sep.join([str(chrval), str(val), str(val+1)])
        # hashv = sourmash.minhash.hash_murmur(sourmash.minhash.to_bytes(outstr), seed=seed)
        # hashv = sourmash.minhash.hash_murmur(outstr.encode('utf-8'), seed=seed)
        hashv = _hash_32(outstr.encode('utf-8'), seed=seed)
        # use max value sys.maxsize to get without mod
        if hashv <= maximum:
            if returnhash:
                if generator:
                    yield hashv
                else:
                    points.append(hashv)
            if generator:
                yield outstr.encode('utf-8')
            else:
                points.append(outstr.encode('utf-8'))
    if not generator:
        return points

def process_file(args):
    '''
    Process a file to calculate the Jaccard similarity between the primary sets and the comparator sets.
    '''
    filepath, primary_sets, primary_keys, mode, parts, subsample = args
    basename = os.path.basename(filepath)
    comparator = bed_to_sets(filename=filepath, mode=mode, parts=parts, subsample=subsample)
    output = {}
    
    for i in range(len(primary_keys)):
        primary = primary_sets[primary_keys[i]]
        union = primary.copy_and_clear()
        union.merge(primary)
        union.merge(comparator)
        
        output[primary_keys[i]] = (primary.similarity(comparator), primary.jaccard(comparator), primary.count_common(comparator), len(union))

    return basename, output

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate Jaccard similarity between BED files.")
    parser.add_argument("filepaths_file", type=str, help="File containing paths of BED files to compare.")
    parser.add_argument("primary_file", type=str, help="File containing paths of primary comparator BED files.")
    parser.add_argument('--output', '-o', type=str, default=sys.stdout, help='The output file path')
    parser.add_argument("--mode", type=str, required=True, help="Mode to indicate what kind of similarity comparison is desired.")
    parser.add_argument("--perm", "-p", type=int, help="Number of permutations for MinHash.")
    parser.add_argument("--subsample", type=float, default=1, help="Subsampling rate for points given as decimal ratio.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    subsample = args.subsample
    mode = args.mode
    outprefix = args.output

    bed_sets = {}
    filepaths_file = args.filepaths_file
    pfile = args.primary_file
    parts = args.perm
    with open(filepaths_file, "r") as filein:
        filepaths = [f.strip() for f in filein.readlines()]
    
    prime_sets = {}
    prime_args = {}
    with open(pfile, "r") as filein:
        # primepaths = [f.strip() for f in filein.readlines()]
        for line in filein.readlines():
            filepath = line.strip()
            basename = os.path.basename(filepath)
            prime_sets[basename] = bed_to_sets(filename=filepath, mode=mode, parts=parts, subsample=subsample)
            
            # prime_args[basename] = (filepath, mode, parts, subsample)
    prime_keys = list(prime_sets.keys())
    prime_keys.sort()

    outlines = [",".join(["bed1", "bed2", "intersect","union"])]
    jacc_long = [",".join(["bed1", "bed2", "similarity", "jaccard"])]

    with Pool() as pool:
        # Create a list of arguments for each worker
        args_list = [(f, prime_sets, prime_keys, mode, parts, subsample) for f in filepaths]
        
        # Map the worker function to the arguments list and get the results
        results = pool.map(process_file, args_list)
        # print(results[:3])

    result_dict = {b:out for b, out in results}
    result_keys = list(result_dict.keys())
    result_keys.sort()
    # Update the matrix with the results
    print(result_keys)
    for j in range(len(result_keys)):
        for i in range(len(prime_keys)):
            primeset = prime_keys[i]
            secset = result_keys[j]
            similarity, jaccard, intersect, union = result_dict[secset][primeset]
            # calc_jacc = intersect/union
            if primeset == secset:
                secset = "-"

            outlines.append(",".join([primeset,
                                    secset,
                                    f"{intersect:>5f}",
                                    f"{union:>5f}"]))
            jacc_long.append(",".join([primeset,
                                    secset,
                                    f"{similarity:>5f}",
                                    f"{jaccard:>5f}"]))
    new_prefix = f"{outprefix}_p{parts}"
    if mode == "C" and subsample != 1:
        new_prefix = f"{new_prefix}_sub{subsample}"
    with open(f"{new_prefix}_card{mode}.csv", mode="w") as outfile:
        for line in outlines:
            outfile.write(line + "\n")
    with open(f"{new_prefix}_jacc{mode}.csv", mode="w") as outfile:
        for line in jacc_long:
            outfile.write(line + "\n")
   