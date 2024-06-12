#!/usr/bin/env python
import xxhash
import sys
import heapq
import os
from multiprocessing import Pool
import numpy as np
from bedprocessing import write_pretty_matrix, basic_bedline
from datasketch import HyperLogLog

VALID_MODES = {"A", "B", "C", "D"}



def parse_bed(filename, parts, mode="A", sep="-"):
    hll_intervals = HyperLogLog(p=parts)
    hll_points = HyperLogLog(p=parts)
    with open(filename, "r") as file:
        for line in file:
            interval = ''
            # points = []
            if not line.startswith('#'):
                chrval, start, end = basic_bedline(line)
                if mode in ["C", "A"]:
                    interval = sep.join([str(chrval), str(start), str(end), "A"])
                    hll_intervals.update(interval.encode('utf-8'))
                if mode in ["C", "B"]:
                    for bp in range(start, end):
                        point=sep.join([str(chrval), str(bp), str(bp + 1), "B"])
                        hll_points.update(point.encode('utf-8'))
    return [hll_intervals, hll_points]

def process_file(args):
    filepath, parts, mode = args
    basename = os.path.basename(filepath)
    return basename, parse_bed(filename=filepath, parts=parts, mode=mode)

def bedline(line, mode, sep):
    # NOTE: right now we are counting i:i+1 to be the position of the ith bp, meaning chrx 2:4 contains bp 2,3 only. see here https://bedtools.readthedocs.io/en/latest/content/general-usage.html
    interval = ''
    points = []
    chrval, start, end = basic_bedline(line)
    if mode in ["A","C"]:
        interval = sep.join([chrval, start, end, "A"])
    if mode in ["B","C"]:
        for val in range(start, end):
            points.append(sep.join([str(chrval), str(val), str(val+1), "B"]))
    return interval, points

#           {"CardA"} {"CardB"} {"Jaccard"} {"AcontainB"} {"BcontainA"}\n\


def print_output(sketch_a, sketch_b, idvals='A: fileA B: fileB IorP'):
    print(f'{idvals} \
          CardA: {sketch_a.cardinality():.2f} \
          CardB: {sketch_b.cardinality():.2f} \
          Jaccard: {sketch_a.jaccard(sketch_b):.2f} \
            AcontainB: {sketch_a.contain(sketch_b):.2f} \
            BcontainA: {sketch_b.contain(sketch_a):.2f}')
    # print(f'{mh_a.get_heapq()}')
    # print(f'{mh_b.get_heapq()}')

def calculate_jaccard_similarity(args):
    i, j, bed_keys, bed_hashes, mode = args
    unionA, unionB = 0, 0
    intersectA, intersectB = 0, 0
    jaccA, jaccB = 0, 0
    if mode in ["A", "C"]:
        # print(bed_keys[j])
        unionA = bed_hashes[bed_keys[i]][0].size_union(bed_hashes[bed_keys[j]][0])  # union of A and B
        jaccA = bed_hashes[bed_keys[i]][0].jaccard(bed_hashes[bed_keys[j]][0])
    if mode in ["B", "C"]:
        # print(bed_keys[j])
        unionB = bed_hashes[bed_keys[i]][1].size_union(bed_hashes[bed_keys[j]][1])  # union of A and B
        jaccB = bed_hashes[bed_keys[i]][1].jaccard(bed_hashes[bed_keys[j]][1])
    return i, j, unionA, unionB, jaccA, jaccB


class HyperLogLog(HyperLogLog):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def cardinality(self):
        return self.count()
    
    def union(self,*others):
        """
        Merges this HyperLogLog sketch with one or more other HyperLogLog sketches.
        """
        for other in others:
            if not isinstance(other, HyperLogLog):
                raise TypeError(f"Cannot merge ExtendedHyperLogLog with {type(other)}")
            self.merge(other)
        return self
    
    def size_union(self, others):
        """
        Returns the cardinality of the union of this HyperLogLog sketch
        with one or more other HyperLogLog sketches.
        """
        
        return self.union(others).count()
    
    def size_intersection(self, b):
        '''Return the estimated size of the intersection
            intersection(a,b)=cardinality(a)+cardinality(b)-union(a,b) '''
        return abs(self.count() + b.count() - self.size_union(b))
    
    def jaccard(self, b):
        '''
        Estimate the Jaccard similarity between A and B.
        Calculate: |sketch(A).insertection(sketch(B))| / |sketch(A).union(sketch(B))|
        '''
        
        return self.size_intersection(b)/self.size_union(b)
    

if __name__ == '__main__':
    filepaths_file = sys.argv[1]
    with open(filepaths_file, "r") as filein:
        filepaths = [f.strip() for f in filein.readlines()]
    parts = int(sys.argv[2])
    mode = sys.argv[3]
    outprefix=""
    if len(sys.argv) == 5:
        if os.path.exists(os.path.dirname(sys.argv[4])):
            outprefix = sys.argv[4]
        else:
            print("The directory given in the output prefix does not exist. Output files will be written to the current working directory using the basename.")
            outprefix=os.path.basename(sys.argv[4])
    # validate_mode(mode=mode)
    verbose = False
    bed_hashes = {}
    onetomany = False
    pfile = "/home/jbonnie1/interval_sketch/hammock/cobind_repro/data/TFbeds_primary.txt"
    with open(pfile, "r") as filein:
        primepaths = [f.strip() for f in filein.readlines()]
    prime_keys = [ os.path.basename(f) for f in primepaths]
    prime_keys.sort()

    with Pool() as pool:
        args_list = [(f, parts, mode) for f in filepaths]
        results = pool.map(process_file, args_list)
    # args_list = [(f, parts, mode) for f in filepaths]
    # results = []
    # for args in args_list:
    #     results.append(process_file(args))

    for result in results:
        basename, bed_set = result
        bed_hashes[basename] = bed_set
    
    secondary_keys = list(bed_hashes.keys())
    for p in prime_keys:
        secondary_keys.remove(p)
    secondary_keys.sort()
    bed_keys = prime_keys + secondary_keys

    matrix = np.zeros((len(bed_keys), len(bed_keys)))
    # outlines = [",".join(["bed1", "bed2", "unionA", "unionB"])]
    # jacc_long = [",".join(["bed1", "bed2", "jaccA", "jaccB"])]


    # Create a pool of workers
    # args_list = [(i, j, bed_keys, bed_hashes, mode) for i in range(len(bed_keys)) for j in range(i, len(bed_keys))]
    results = []
    # for args in args_list:
    #     # if bed_keys[args[0]] != bed_keys[args[1]]:
    #     #     print(bed_keys[args[0]], bed_keys[args[1]])
    #     results.append(calculate_jaccard_similarity(args))

    with Pool() as pool:
        # Create a list of arguments for each worker
        args_list = [(i, j, bed_keys, bed_hashes, mode) for i in range(len(prime_keys)) for j in range(i, len(bed_keys))]
        # Map the worker function to the arguments list and get the results
        results = pool.map(calculate_jaccard_similarity, args_list)
    

    # Get filehandles for long output files
    cardout=open(f"{outprefix}_hll_p{parts}_card{mode}.csv", mode="w")
    cardout.write(",".join(["bed1", "bed2", "unionA", "unionB"]) + "\n")
    jaccout=open(f"{outprefix}_hll_p{parts}_jacc{mode}.csv", mode="w")
    jaccout.write(",".join(["bed1", "bed2", "jaccA", "jaccB"]) + "\n")
    # Update the matrix with the results
    for result in results:
        i, j, unionA, unionB, jaccA, jaccB = result
        matrix[i, j] = jaccA
        matrix[j, i] = jaccB

        set1 = bed_keys[i]
        set2 = bed_keys[j]
        if set1 == set2:
            set2 = "-"

        # outlines.append(",".join([set1,
        #                         set2,
        #                         f"{unionA:>5f}",
        #                         f"{unionB:>5f}"]))
        cardout.write(",".join([set1,
                                set2,
                                f"{unionA:>5f}",
                                f"{unionB:>5f}"])+ "\n")
        jaccout.write(",".join([set1,
                                set2,
                                f"{jaccA:>5f}",
                                f"{jaccB:>5f}"])+ "\n")
        # jacc_long.append(",".join([set1,
        #                         set2,
        #                         f"{jaccA:>5f}",
        #                         f"{jaccB:>5f}"]))

    with open(f"{outprefix}_hll_p{parts}_matrix{mode}.csv", mode="w") as matout:
        write_pretty_matrix(matrix=matrix, row_names=bed_keys, outhandle=matout)
    
    cardout.close()
    jaccout.close()
    # with open(f"{outprefix}minhash_h{num_hash}_card{mode}.csv", mode="w") as outfile:
    #     for line in outlines:
    #         outfile.write(line + "\n")
    # with open(f"{outprefix}minhash_h{num_hash}_jacc{mode}.csv", mode="w") as outfile:
    #     for line in jacc_long:
    #         outfile.write(line + "\n")

    # if mode in {"A", "C"}:
    #     # print(list(sorted(mh_intervals_a._heapq)))
    #     # print(list(sorted(mh_intervals_b._heapq)))
    #     union_sketch_intervals = mh_intervals_a.size_union(mh_intervals_b)
    #     idvals = f'A: {os.path.basename(fn_a)} B: {os.path.basename(fn_b)} {"Intervals"}'
    #     print_output(mh_intervals_a, mh_intervals_b, idvals)
    # if mode in {"B", "C"}:
    #     print(list(sorted(mh_points_a._heapq)))
    #     print(list(sorted(mh_points_b._heapq)))
    #     idvals = f'A: {os.path.basename(fn_a)}  B: {os.path.basename(fn_b)} {"Points"}'
    #     union_sketch_points = mh_points_a.size_union(mh_points_b)
    #     print_output(mh_points_a, mh_points_b, idvals)
    # print(f'{mh_a.cardinality():.2f} {mh_b.cardinality():.2f} \
    #       {mh_a.jaccard(mh_b):.2f} {mh_a.contain(mh_b):.2f} \
    #         {mh_b.contain(mh_a):.2f}')
    # print(f'{mh_a.get_heapq()}')
    # print(f'{mh_b.get_heapq()}')
