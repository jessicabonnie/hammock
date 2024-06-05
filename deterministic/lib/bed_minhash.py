#!/usr/bin/env python
import xxhash
import sys
import heapq
import os
from multiprocessing import Pool
import numpy as np
from bedprocessing import write_pretty_matrix, basic_bedline

VALID_MODES = {"A", "B", "C", "D"}


def get_words(filename, mode):
    return open(filename, "rb").read().decode("utf8", "ignore").strip().split()


def validate_mode(mode):
    if mode not in VALID_MODES:
        raise ValueError("validate_mode: mode must be one of %r." % VALID_MODES)

# def basic_bedline(line):
#     columns = line.strip().split('\t')
#     if 0 < len(columns) <= 2:
#         raise ValueError("bedline: one of the lines in malformed")
#     if columns[0].startswith('chr'):
#         columns[0] = columns[0][3:]
#     return int(columns[0]), int(columns[1]), int(columns[2])

def parse_bed(filename, num_hash, mode="A", sep="-"):
    mh_intervals = MinHash(num_hash, 0)
    mh_points = MinHash(num_hash, 0)
    with open(filename, "r") as file:
        for line in file:
            interval = ''
            points = []
            if not line.startswith('#'):
                chrval, start, end = basic_bedline(line)
                if mode in ["C", "A"]:
                    interval = sep.join([str(chrval), str(start), str(end), "A"])
                    mh_intervals.update(interval)
                if mode in ["C", "B"]:
                    for bp in range(start, end):
                        point=sep.join([str(chrval), str(bp), str(bp + 1), "B"])
                        mh_points.update(point)
        return [mh_intervals, mh_points]

def process_file(args):
    filepath, num_hash, mode = args
    basename = os.path.basename(filepath)
    return basename, parse_bed(filename=filepath, num_hash=num_hash, mode=mode)

# def bedline(line, mode, sep):
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


def print_output(mh_a, mh_b, idvals='A: fileA B: fileB IorP'):
    print(f'{idvals} \
          CardA: {mh_a.cardinality():.2f} \
          CardB: {mh_b.cardinality():.2f} \
          Jaccard: {mh_a.jaccard(mh_b):.2f} \
            AcontainB: {mh_a.contain(mh_b):.2f} \
            BcontainA: {mh_b.contain(mh_a):.2f}')
    # print(f'{mh_a.get_heapq()}')
    # print(f'{mh_b.get_heapq()}')

def calculate_jaccard_similarity(args):
    i, j, bed_keys, bed_hashes, mode = args
    unionA, unionB = 0, 0
    intersectA, intersectB = 0, 0
    jaccA, jaccB = 0, 0
    if mode in ["A", "C"]:
        unionA = bed_hashes[bed_keys[i]][0].size_union(bed_hashes[bed_keys[j]][0])  # union of A and B
        jaccA = bed_hashes[bed_keys[i]][0].jaccard(bed_hashes[bed_keys[j]][0])
    if mode in ["B", "C"]:
        unionB = bed_hashes[bed_keys[i]][1].size_union(bed_hashes[bed_keys[j]][1])  # union of A and B
        jaccB = bed_hashes[bed_keys[i]][1].jaccard(bed_hashes[bed_keys[j]][1])
    return i, j, unionA, unionB, jaccA, jaccB

class HashXX32(object):
    def __init__(self, seed):
        self.h = xxhash.xxh32(seed=seed)

    def hash(self, o):
        self.h.reset()
        self.h.update(o)
        return self.h.intdigest() % sys.maxsize


class MinHash(object):
    '''
    MinHash data structure.
    The hash values should be stored in a heapq data structure for fast update.
    See https://docs.python.org/3/library/heapq.html
    '''

    def __init__(self, num_item, seed):
        self._hasher = HashXX32(seed)
        self._maxhash = 2**32
        # Multiply every item in the heapq data structure by -1 to form a max heap.
        # Thereby, self._heapq[0] * -1 = max([i * -1 for i in self._heapq])
        self._heapq = [-self._maxhash for i in range(num_item)]
        heapq.heapify(self._heapq)
        self._k = num_item

    def update(self, obj):
        '''
        Update obj to MinHash.

        We store -1 * (hashed result) in `self._heapq`, by which
        (-1 * self._heapq[0]) is the maximal value of the sketch.
        '''
        new_hash = -1*self._hasher.hash(obj)
        if new_hash > self._heapq[0]:
            if new_hash not in self._heapq:
                heapq.heappushpop(self._heapq, new_hash)

    def copy(self):
        new = MinHash(self._k, 0)
        new._heapq = self._heapq[:]
        return new

    def get_heapq(self):
        ''' Return the MinHash sketch. '''
        return self._heapq

    def cardinality(self):
        ''' Return the estimated cardinality of the sketch. '''
        kmv = -1*self._heapq[0]
        numerator = self._k * self._maxhash
        return (numerator/kmv)-1

    def union(self, b):
        ''' Return the union sketch (A and B). '''
        out = self.copy()

        reject = min(out._heapq)
        for i in range(1, b._k):
            # if reject < b._heapq[b._k - i]:
            if b._heapq[b._k - i] not in out._heapq:
                if b._heapq[b._k - i] > reject:
                    heapq.heappushpop(out._heapq, b._heapq[b._k - i])
                    reject = min(out._heapq)
            # else:
                # break
        return out

    def size_intersection(self, b):
        '''Return the estimated size of the intersection
            intersection(a,b)=cardinality(a)+cardinality(b)-union(a,b) '''
        return abs(self.cardinality() + b.cardinality() - self.union(b).cardinality())

    def size_union(self, b):
        ''' Return the size of union sketch (|A and B|). '''
        return self.union(b).cardinality()

    def contain(self, b):
        '''
        Estimate the containment of A in B.
        Calculate: |sketch(A).insertection(sketch(B))| / |sketch(A)|
        '''

        return self.size_intersection(b) / self.cardinality()

    # def jaccard(self, b):
    #     '''
    #     Estimate the Jaccard similarity between A and B.
    #     Calculate: |sketch(A).insertection(sketch(B))| / |sketch(A).union(sketch(B))|
    #     '''
    #     return self.size_intersection(b)/self.size_union(b)
    def jaccard(self, b):
        '''
        Estimate the Jaccard similarity between A and B.
        Calculate: |sketch(A).insertection(sketch(B))| / |sketch(A).union(sketch(B))|
        '''
        # print(list(sorted(list(self._heapq)))[:10])
        # print(list(sorted(list(b._heapq)))[:10])
        # print(list(sorted(self._heapq))[:10])
        # print(list(sorted(b._heapq))[:10])
    #     union_sketch_intervals = self.size_union(b)
    #     print_output(self, b, "blah")
        set_self = set(self._heapq)
        set_other = set(b._heapq)
        # print(list(sorted(set(self._heapq)))[:10])
        # print(list(sorted(set(b._heapq)))[:10])
        # return len(set_self.intersection(set_other))/len(set_self.union(set_other))
        return len(set_self.intersection(set_other))/len(self._heapq) 
        # return self.size_intersection(b)/self.size_union(b)


def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size


if __name__ == '__main__':
    filepaths_file = sys.argv[1]
    with open(filepaths_file, "r") as filein:
        filepaths = [f.strip() for f in filein.readlines()]
    num_hash = int(sys.argv[2])
    mode = sys.argv[3]

    outprefix=""
    validate_mode(mode=mode)
    verbose = False
    bed_hashes = {}

    # with Pool() as pool:
    #     args_list = [(f, num_hash, mode) for f in filepaths]
    #     results = pool.map(process_file, args_list)
    args_list = [(f, num_hash, mode) for f in filepaths]
    results = []
    for args in args_list:
        results.append(process_file(args))

    for result in results:
        basename, bed_set = result
        bed_hashes[basename] = bed_set
    
    bed_keys = list(bed_hashes.keys())
    bed_keys.sort()
    matrix = np.zeros((len(bed_keys), len(bed_keys)))
    outlines = [",".join(["bed1", "bed2", "unionA", "unionB"])]


    # Create a pool of workers
    args_list = [(i, j, bed_keys, bed_hashes, mode) for i in range(len(bed_keys)) for j in range(i, len(bed_keys))]
    results = []
    for args in args_list:
        # if bed_keys[args[0]] != bed_keys[args[1]]:
        #     print(bed_keys[args[0]], bed_keys[args[1]])
        results.append(calculate_jaccard_similarity(args))

    # with Pool() as pool:
    #     # Create a list of arguments for each worker
    #     args_list = [(i, j, bed_keys, bed_hashes, mode) for i in range(len(bed_keys)) for j in range(i, len(bed_keys))]
    #     # Map the worker function to the arguments list and get the results
    #     results = pool.map(calculate_jaccard_similarity, args_list)

    # Update the matrix with the results
    for result in results:
        i, j, unionA, unionB, jaccA, jaccB = result
        matrix[i, j] = jaccA
        matrix[j, i] = jaccB

        set1 = bed_keys[i]
        set2 = bed_keys[j]
        if set1 == set2:
            set2 = "-"

        outlines.append(",".join([set1,
                                set2,
                                f"{unionA:>5f}",
                                f"{unionB:>5f}"]))

    with open(f"{outprefix}minhash_h{num_hash}_matrix{mode}.csv", mode="w") as matout:
        write_pretty_matrix(matrix=matrix, row_names=bed_keys, outhandle=matout)
    with open(f"{outprefix}minhash_h{num_hash}_card{mode}.csv", mode="w") as outfile:
        for line in outlines:
            outfile.write(line + "\n")

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
