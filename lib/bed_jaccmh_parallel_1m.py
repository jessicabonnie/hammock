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
    return xxhash.xxh32(x.encode('utf-8')).intdigest()

# def subsample_points(points, subsample=1, seed=0):
#     '''
#     Subsample a list of points.
#     '''
#     if subsample == 1:
#         return points
#     else:
# #         return [ b for b in points if _hash_64(b,seed=seed) % 10 <= subsample*10]
# class ExtendMinHash(LeanMinHash):
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         # self.hashfunc = _hash_32
#     # def cardinality(self):
#     #     return self.count()
#     # def union(self, b):
#     #     """
#     #     Merges this sketch with one or more other sketches.
#     #     """
#     #     return( MinHash.union(self,b))
#         # for other in others:
#         #     if not isinstance(other, HyperLogLog):
#         #         raise TypeError(f"Cannot merge ExtendedHyperLogLog with {type(other)}")
#             # self.merge(other)
#         # return self
#     # def update_batch(self, new_vals):
#     #     newMH = 
    
#     def size_union(self, others):
#         """
#         Returns the cardinality of the union of this MinHash sketch
#         with one or more other MinHash sketches.
#         """
#         return MinHash.union(self,others).count()
    
#     def size_intersection(self, b):
#         '''Return the estimated size of the intersection
#             intersection(a,b)=cardinality(a)+cardinality(b)-union(a,b) '''
#         # return self.intersect(b).count()
#         return abs(self.count() + b.count() - self.size_union(b))

#     # # def add(self,values):
#     #     '''
#     #     Add a value to the sketch.
#     #     '''
#     #     self.update(values)
#     # def jaccard(self, b):
#     #     '''
#     #     Estimate the Jaccard similarity between A and B.
#     #     Calculate: |sketch(A).insertection(sketch(B))| / |sketch(A).union(sketch(B))|
#     #     '''
#     #     return( MinHash.jaccard(self,b))
#         # return self.size_intersection(b)/self.size_union(b)

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
    
    outset = MinHash(num_perm=int(parts), seed=0, hashfunc=hashfunc)
    #, hashfunc=xxhash.xxh32) #, hashfunc=_hash_32)
    listA = []
    listB = []
    bedline_args = []
    with open(filename, "r") as file:
        for line in file:
            # print(line)
            if not line.startswith('#'):
                bedline_args.append((line, mode, sep, subsample))
                # interval, points = bedline(line, mode=mode, sep=sep, subsample=subsample)
                # print(len(points))
                # if interval == '':
                #     break
                # listA.append(interval)
                # listB.append(points)
    with Pool() as pool:
        bedline_results = pool.starmap(bedline, bedline_args)
    if mode in ["B", "C"]:
        # outset.update_batch(listB)
        outset.update_batch([b for (_,b) in bedline_results if b != []])
    if mode in ["A", "C"]:
        # outset.update_batch(listA)
        outset.update_batch([a for (a,_) in bedline_results if a != ''])

    print(f"finished bed_to_sets for {filename}")
    # print(outset.count())
    return LeanMinHash(outset)
    # return ExtendMinHash(outset)

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

# def subsample_max(subsample):
#     return subsample*sys.maxsize
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
    union = LeanMinHash.union(set1,set2)
    unionc = union.count()
    intersect = set1.count() + set2.count() - unionc
    print(unionc, intersect)
    # intersect = float(len(set1.intersection(set2)))
    # return union, intersect, intersect/union
    return (unionc, intersect, set1.jaccard(set2), intersect/unionc)

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
        # union = MinHash.union(primary, comparator)
        output[primary_keys[i]] = similarity_values(primary, comparator)
    
    return basename, output

# def calculate_jaccard_similarity(set1, set2):
#     '''
#     Calculate the Jaccard similarity between two sets.
#     '''
#     # unionA, unionB = 0, 0
#     outunion = 0
#     outintersect = 0
#     outjacc = 0
#     # intersectA, intersectB = 0, 0
#     # jaccA, jaccB = 0, 0
#     outunion, outintersect, outjacc1, outjacc2 = jaccard_similarity(set1, set2)
#     return outunion, outjacc1, outjacc2

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate Jaccard similarity between BED files.")
    parser.add_argument("filepaths_file", type=str, help="File containing paths of BED files to compare.")
    parser.add_argument("primary_file", type=str, help="File containing paths of primary comparator BED files.")
    parser.add_argument('--output', '-o', type=str, default=sys.stdout, help='The output file path')
    parser.add_argument("--mode", type=str, required=True, help="Mode to indicate what kind of similarity comparison is desired.")
    parser.add_argument("--perm", "-p", type=int, help="Number of permutations for MinHash.")
    parser.add_argument("--subsample", type=float, default=1, help="Subsampling rate for points: provide decimal ratio.")
    return parser.parse_args()


if __name__ == "__main__":
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
    # readargs = []
    with open(pfile, "r") as filein:
        # primepaths = [f.strip() for f in filein.readlines()]
        for line in filein.readlines():
            filepath = line.strip()
            basename = os.path.basename(filepath)
            # readargs.append( (basename, filepath, mode, parts, subsample))
            prime_sets[basename] = bed_to_sets(filename=filepath, mode=mode, parts=parts, subsample=subsample)
    # with Pool() as xpool:
    #     linesets = xpool.map(bed_to_sets_parallel, readargs)
    # prime_sets = {b:out for b, out in linesets}
    prime_keys = list(prime_sets.keys())
    prime_keys.sort()

    outlines = [",".join(["bed1", "bed2", "union", "intersect_inex"])]
    jacc_long = [",".join(["bed1", "bed2", "jaccardfunc", "jaccardcalc"])]
    # Create a list of arguments for each worker
    args_list = [(f, prime_sets, prime_keys, mode, parts, subsample) for f in filepaths]
    results = []
    # with Pool() as pool:
        # Map the worker function to the arguments list and get the results
        # results = pool.map(process_file, args_list)
    for argy in args_list:
        results.append(process_file(argy))

    result_dict = {b:out for (b, out) in results}
    result_keys = list(result_dict.keys())
    result_keys.sort()
    for j in range(len(result_keys)):
            for i in range(len(prime_keys)):
                primeset = prime_keys[i]
                secset = result_keys[j]
                union, intersect, jaccard, jacc_calc = result_dict[secset][primeset]
                # calc_jacc = intersect/union
                if primeset == secset:
                    secset = "-"

                outlines.append(",".join([primeset,
                                        secset,
                                        f"{intersect:>5f}",
                                        f"{union:>5f}"]))
                jacc_long.append(",".join([primeset,
                                        secset,
                                        f"{jaccard:>5f}",
                                        f"{jacc_calc:>5f}"]))
    new_prefix = f"{outprefix}_p{parts}"
    if mode == "C" and subsample != 1:
        new_prefix = f"{new_prefix}_sub{subsample}"
    with open(f"{new_prefix}_card{mode}.csv", mode="w") as outfile:
        for line in outlines:
            outfile.write(line + "\n")
    with open(f"{new_prefix}_jacc{mode}.csv", mode="w") as outfile:
        for line in jacc_long:
            outfile.write(line + "\n")
