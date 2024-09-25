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


# class HashXX32(object):
#     def __init__(self, seed):
#         self.h = xxhash.xxh32(seed=seed)

#     def hash(self, o):
#         self.h.reset()
#         self.h.update(o)
#         return self.h.intdigest() % sys.maxsize
def _hash_64(b, seed=0):
    return xxhash.xxh64(b, seed=seed).intdigest() % sys.maxsize

def _hash_32(b):
    return xxhash.xxh32(b).intdigest() % sys.maxsize

def subsample_points(points, subsample=1, seed=0):
    '''
    Subsample a list of points.
    '''
    if subsample == 1:
        return points
    else:
        return [ b for b in points if _hash_64(b,seed=seed) % 10 <= subsample*10]
class ExtendMinHash(MinHash):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    
    def cardinality(self):
        return self.count()
    
    def union(self, b):
        """
        Merges this sketch with one or more other sketches.
        """
        
        return( MinHash.union(self,b))
        # for other in others:
        #     if not isinstance(other, HyperLogLog):
        #         raise TypeError(f"Cannot merge ExtendedHyperLogLog with {type(other)}")
            # self.merge(other)
        # return self
    
    def size_union(self, others):
        """
        Returns the cardinality of the union of this MinHash sketch
        with one or more other MinHash sketches.
        """
        
        return self.union(others).count()
    
    def size_intersection(self, b):
        '''Return the estimated size of the intersection
            intersection(a,b)=cardinality(a)+cardinality(b)-union(a,b) '''
        # return self.intersect(b).count()
        return abs(self.count() + b.count() - self.size_union(b))

    # def add(self,values):
        '''
        Add a value to the sketch.
        '''
        self.update(values)
    # def jaccard(self, b):
    #     '''
    #     Estimate the Jaccard similarity between A and B.
    #     Calculate: |sketch(A).insertection(sketch(B))| / |sketch(A).union(sketch(B))|
    #     '''
    #     return( MinHash.jaccard(self,b))
        # return self.size_intersection(b)/self.size_union(b)

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
    # outset = ExtendMinHash(num_perm=int(parts), hashfunc=_hash_64)
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

# def basic_bedline(line):
#     columns = line.strip().split('\t')
#     if 0 < len(columns) <= 2:
#         raise ValueError("bedline: one of the lines in malformed")
#     if columns[0].startswith('chr'):
#         columns[0] = columns[0][3:]
#     return str(columns[0]), int(columns[1]), int(columns[2])

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
        # interval = _hash_64(interval)
        interval = sourmash.minhash.hash_murmur(interval, 0)
    return interval, [g for g in points]


def generate_points(chrval, start, end, sep="-", subsample=1, returnhash=False, seed=0, generator=True):
    '''
    Generate points from a BED interval.
    '''
    if not generator:
        points = []
    for val in range(start, end):
        outstr = sep.join([str(chrval), str(val), str(val+1)])
        # hashv = sourmash.minhash.hash_murmur(sourmash.minhash.to_bytes(outstr), seed=seed)
        hashv = sourmash.minhash.hash_murmur(outstr.encode('utf-8'), seed=seed)
        # hashv = _hash_64(outstr, seed=seed)
        if hashv % 10 <= subsample*10:
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
    


# def jaccard_similarity(set1, set2):
#     union = set1.size_union(set2)
#     intersect = set1.size_intersection(set2)

#     # intersect = float(len(set1.intersection(set2)))
#     # return union, intersect, intersect/union
#     return union, intersect, set1.jaccard(set2), intersect/union

# def jaccard_similarity(set1, set2):
#     # union = set1.update(set2)
#     intersect, _ = set1.intersection(set2)
#     return set1.compare(set2), set1.jaccard(set2), set1.angular_similarity(set2), intersect


# def write_pretty_matrix(matrix, row_names, outhandle=sys.stdout):
#     """
#     Prints a matrix in a pretty format with row names.

#     Args:
#         matrix (np.ndarray): The matrix to be printed.
#         row_names (list): A list of row names, one for each row in the matrix.
#         outhandle (file-like object, optional): The output handle to write the matrix to. Defaults to sys.stdout.
#     """
#     # Get the number of rows and columns in the matrix
#     num_rows, num_cols = matrix.shape
#     with outhandle as matout:
#         # Print the header row with column names
#         matout.write(",".join([""] + row_names)+"\n")
#         for i in range(num_rows):
#             outrow = [row_names[i]] + [f"{matrix[i, j]:>5.4f}" for j in range(num_cols)]
#             matout.write(",".join(outrow)+"\n")

def process_file(args):
    '''
    Process a file to calculate the Jaccard similarity between the primary sets and the comparator sets.
    '''
    # def jaccard_similarity(set1, set2):
    #     # union = set1.update(set2)
    #     intersect, _ = set1.intersection(set2)
    #     return set1.compare(set2), set1.jaccard(set2), set1.angular_similarity(set2), intersect

    filepath, primary_args, primary_keys, mode, parts, subsample = args
    basename = os.path.basename(filepath)
    primary_sets = {}
    # for basename in primary_keys:
    #     basename = os.path.basename(pfilepath)
    #     pfilepath, pmode, pparts, psubsample = primary_args[basename]
    #     primary_sets[basename] = bed_to_sets(filename=pfilepath, mode=pmode, parts=pparts, subsample=psubsample)
    
    comparator = bed_to_sets(filename=filepath, mode=mode, parts=parts, subsample=subsample)
    
    output = {}
    for i in range(len(primary_keys)):
        pfilepath, pmode, pparts, psubsample = primary_args[primary_keys[i]]
        primary = bed_to_sets(filename=pfilepath, mode=pmode, parts=pparts, subsample=psubsample)
        print(primary.jaccard(comparator))
        # intersect, _ = primary.intersection(comparator)
        
        # args = (primary_sets[i], comparator, mode)
        # union, jaccard1, jaccard2, _ = jaccard_similarity(primary, comparator)
        # output[primary_keys[i]] = (union, jaccard1, jaccard2)
       
        union = primary.copy_and_clear()
        union.merge(primary)
        union.merge(comparator)
        primary.merge(comparator)
        output[primary_keys[i]] = (primary.similarity(comparator), primary.jaccard(comparator), primary.count_common(comparator), len(union))

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
    parser.add_argument("--subsample", type=int, default=1, help="Subsampling rate for points.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    subsample = args.subsample
    mode = args.mode
    outprefix = args.output
    # if len(sys.argv) < 5:
    #     print("Usage: python bed_jaccards.py <filepaths file> <file of paths of primary comparitor files> <mode> <parts> <prefix>")
    #     sys.exit(1)

    bed_sets = {}
    filepaths_file = args.filepaths_file
    pfile = args.primary_file
    # pfile = "/home/jbonnie1/interval_sketch/hammock/cobind_repro/data/TFbeds_primary.txt"
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
            # prime_sets[basename] = bed_to_sets(filename=filepath, mode=mode, parts=parts, subsample=subsample)
            prime_args[basename] = (filepath, mode, parts, subsample)
    prime_keys = list(prime_args.keys())
    prime_keys.sort()

    # matrix = np.zeros((len(bed_keys), len(bed_keys)))
    outlines = [",".join(["bed1", "bed2", "intersect","union"])]
    jacc_long = [",".join(["bed1", "bed2", "similarity", "jaccard"])]

    with Pool() as pool:
        # Create a list of arguments for each worker
        args_list = [(f, prime_args, prime_keys, mode, parts, subsample) for f in filepaths]
        # Map the worker function to the arguments list and get the results
        results = pool.map(process_file, args_list)

    result_dict = {b:out for b, out in results}
    result_keys = list(result_dict.keys())
    result_keys.sort()
    # Update the matrix with the results
    for j in range(len(result_keys)):
        for i in range(len(prime_keys)):
            primeset = prime_keys[i]
            secset = result_keys[j]
            similarity, jaccard, intersect, union = result_dict[secset][primeset]
            # matrix[i, j] = jaccA
            # matrix[j, i] = jaccB

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
   
    with open(f"{outprefix}_p{parts}_card{mode}.csv", mode="w") as outfile:
        for line in outlines:
            outfile.write(line + "\n")
    with open(f"{outprefix}_p{parts}_jacc{mode}.csv", mode="w") as outfile:
        for line in jacc_long:
            outfile.write(line + "\n")
    # matrix = np.zeros((len(bed_keys), len(bed_keys)))
    # outlines = [",".join(["bed1", "bed2", "unionA", "unionB"])]
    
    #     # freeze the order of keys for bed_sets in a list
    # bed_keys = list(bed_sets.keys())
    # bed_keys.sort()
    # matrix = np.zeros((len(bed_keys), len(bed_keys)))
    # outlines = [",".join(["bed1", "bed2", "unionA", "unionB"])]
    # jacc_long = [",".join(["bed1", "bed2", "jaccA", "jaccB"])]

    # # Update the matrix with the results
    # for result in results:
    #     i, j, unionA, unionB, jaccA, jaccB = result
    #     matrix[i, j] = jaccA
    #     matrix[j, i] = jaccB
    # pretty_print_matrix(matrix,bed_keys)
    # with open(f"{outprefix}_matrix{mode}.csv", mode="w") as matout:
    #     write_pretty_matrix(matrix=matrix, row_names=bed_keys, outhandle=matout)
