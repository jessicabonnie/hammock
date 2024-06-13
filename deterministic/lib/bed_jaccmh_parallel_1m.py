#!/usr/bin/env python
# import random
import sys
import numpy as np
import os
from multiprocessing import Pool
from datasketch import MinHash


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
        Returns the cardinality of the union of this HyperLogLog sketch
        with one or more other HyperLogLog sketches.
        """
        
        return self.union(others).count()
    
    def size_intersection(self, b):
        '''Return the estimated size of the intersection
            intersection(a,b)=cardinality(a)+cardinality(b)-union(a,b) '''
        # return self.intersection(b).count()
        return abs(self.count() + b.count() - self.size_union(b))

    def add(self,values):
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

def bed_to_sets(filename, mode, parts, sep="-", verbose=False):
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
    Aset = ExtendMinHash(num_perm=int(parts)) #set()
    Bset = ExtendMinHash(num_perm=int(parts))#set()
    with open(filename, "r") as file:
        for line in file:
            if not line.startswith('#'):
                interval, points = bedline(line, mode=mode, sep=sep)
                if verbose:
                    print(line)
                    if len(points) > 100:
                        print(points[:10])
                if interval != '':
                    Aset.update(interval)
                if mode in ["C", "B"]:
                    Bset.batch_update(points)
    return [Aset, Bset]


def basic_bedline(line):
    columns = line.strip().split('\t')
    if 0 < len(columns) <= 2:
        raise ValueError("bedline: one of the lines in malformed")
    if columns[0].startswith('chr'):
        columns[0] = columns[0][3:]
    return str(columns[0]), int(columns[1]), int(columns[2])


def bedline(line, mode, sep):
    # NOTE: right now we are counting i:i+1 to be the position of the ith bp, meaning chrx 2:4 contains bp 2,3 only. see here https://bedtools.readthedocs.io/en/latest/content/general-usage.html
    interval = ''
    points = []
    chrval, start, end = basic_bedline(line)
    if mode in ["A","C"]:
        interval = sep.join([chrval, str(start), str(end), "A"]).encode('utf-8')
    if mode in ["B","C"]:
        for val in range(start, end):
            points.append(sep.join([str(chrval), str(val), str(val+1), "B"]).encode('utf-8'))
    return interval, points


def jaccard_similarity(set1, set2):
    # union = set1.size_union(set2)
    # intersect = float(len(set1.intersection(set2)))
    # return union, intersect, intersect/union
    return set1.size_union(set2), set1.size_intersection(set2), set1.jaccard(set2)


def write_pretty_matrix(matrix, row_names, outhandle=sys.stdout):
    """
    Prints a matrix in a pretty format with row names.

    Args:
        matrix (np.ndarray): The matrix to be printed.
        row_names (list): A list of row names, one for each row in the matrix.
        outhandle (file-like object, optional): The output handle to write the matrix to. Defaults to sys.stdout.
    """
    # Get the number of rows and columns in the matrix
    num_rows, num_cols = matrix.shape
    with outhandle as matout:
        # Print the header row with column names
        matout.write(",".join([""] + row_names)+"\n")
        for i in range(num_rows):
            outrow = [row_names[i]] + [f"{matrix[i, j]:>5.4f}" for j in range(num_cols)]
            matout.write(",".join(outrow)+"\n")

def process_file(args):
    '''
    Process a file to calculate the Jaccard similarity between the primary sets and the comparator sets.
    '''
    filepath, primary_sets, primary_keys, mode, parts = args
    basename = os.path.basename(filepath)
    comparator = bed_to_sets(filename=filepath, mode=mode, parts=parts)
    # print(basename)
    output = {}
    for i in range(len(primary_keys)):
        # args = (primary_sets[i], comparator, mode)
        unionA, unionB, jaccA, jaccB = calculate_jaccard_similarity(primary_sets[primary_keys[i]], comparator, mode)
        output[primary_keys[i]] = (unionA, unionB, jaccA, jaccB)
    return basename, output

def calculate_jaccard_similarity(set1, set2, mode):
    '''
    Calculate the Jaccard similarity between two sets.
    '''
    unionA, unionB = 0, 0
    intersectA, intersectB = 0, 0
    jaccA, jaccB = 0, 0
    if mode in ["A", "C"]:
        unionA, intersectA, jaccA = jaccard_similarity(set1[0], set2[0])
    if mode in ["B", "C"]:
        unionB, intersectB, jaccB = jaccard_similarity(set1[1], set2[1])
    return unionA, unionB, jaccA, jaccB

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python bed_jaccards.py <filepaths file> <mode> <parts> <prefix>")
        sys.exit(1)
    outprefix = ""#os.getcwd()
    if len(sys.argv) == 5:
        if os.path.exists(os.path.dirname(sys.argv[4])):
            outprefix = sys.argv[4]
        else: 
            print("The directory given in the output prefix does not exist. Output files will be written to the current working directory using the basename.")
            outprefix=os.path.basename(sys.argv[3])
    bed_sets = {}
    mode = sys.argv[2]
    filepaths_file = sys.argv[1]
    parts = sys.argv[3]
    with open(filepaths_file, "r") as filein:
        filepaths = [f.strip() for f in filein.readlines()]
    
    pfile = "/home/jbonnie1/interval_sketch/hammock/cobind_repro/data/TFbeds_primary.txt"
    prime_sets = {}
    with open(pfile, "r") as filein:
        # primepaths = [f.strip() for f in filein.readlines()]
        for line in filein.readlines():
            filepath = line.strip()
            basename = os.path.basename(filepath)
            prime_sets[basename] = bed_to_sets(filename=filepath, mode=mode, parts=parts)
    prime_keys = list(prime_sets.keys())
    prime_keys.sort()

    # matrix = np.zeros((len(bed_keys), len(bed_keys)))
    outlines = [",".join(["bed1", "bed2", "unionA", "unionB"])]
    jacc_long = [",".join(["bed1", "bed2", "jaccA", "jaccB"])]

    with Pool() as pool:
        # Create a list of arguments for each worker
        args_list = [(f, prime_sets, prime_keys, mode, parts) for f in filepaths]
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
            unionA, unionB, jaccA, jaccB = result_dict[secset][primeset]
            # matrix[i, j] = jaccA
            # matrix[j, i] = jaccB

            if primeset == secset:
                secset = "-"

            outlines.append(",".join([primeset,
                                    secset,
                                    f"{unionA:>5f}",
                                    f"{unionB:>5f}"]))
            jacc_long.append(",".join([primeset,
                                    secset,
                                    f"{jaccA:>5f}",
                                    f"{jaccB:>5f}"]))
   
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
