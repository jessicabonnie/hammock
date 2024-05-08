#!/usr/bin/env python
# import random
import sys
import numpy as np
import os

# given the starting position of the first interval in the bedfile (start_pos), the number of intervals/lines desired (count), the distance between the start and end of each interval (size), the chromosome value (chrom), and the distance between the end position of one interval and the start position of the next interval, generate a bedfile
# def true_cardinality(obj_list,seed=10):
#     hasher = xxh32(seed)
#     maxhash = 2**32
#     all_val = set()
#     for obj in obj_list:
#         all_val.add(hasher.hash(obj))
#     return(len(all_val))


def bed_to_sets(filename, mode="A", sep="-", verbose=False):
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
    Aset = set()
    Bset = set()
    with open(filename, "r") as file:
        for line in file:
            if not line.startswith('#'):
                interval, points = bedline(line, mode=mode, sep=sep)
                if verbose:
                    print(line)
                    if len(points) > 100:
                        print(points[:10])
                if interval != '':
                    Aset.add(interval)
                if mode in ["C", "B"]:
                    Bset.update(points)
    return [Aset, Bset]


def basic_bedline(line):
    columns = line.strip().split('\t')
    if 0 < len(columns) <= 2:
        raise ValueError("bedline: one of the lines in malformed")
    if columns[0].startswith('chr'):
        columns[0] = columns[0][3:]
    return int(columns[0]), int(columns[1]), int(columns[2])


def bedline(line, mode, sep) -> tuple[str, list]:
    # NOTE: right now we are counting i:i+1 to be the position of the ith bp, meaning chrx 2:4 contains bp 2,3 only. see here https://bedtools.readthedocs.io/en/latest/content/general-usage.html
    interval = ''
    points = []
    columns = line.strip().split('\t')
    if 0 < len(columns) <= 2:
        print(line)
        raise ValueError("bedline: one of the lines in malformed")
    if columns[0].startswith('chr'):
        columns[0] = columns[0][3:]
    if mode in ["A","C"]:
        interval = sep.join([columns[0], columns[1], columns[2], "A"])
    if mode in ["B","C"]:
        for val in range(int(columns[1]), int(columns[2])):
            points.append(sep.join([columns[0], str(val), str(val+1), "B"]))
    # elif mode == "C":
    #     interval = sep.join([columns[0], columns[1], columns[2], "A"])
    #     for val in range(int(columns[1]), int(columns[2])):
    #         points.append(sep.join([columns[0], str(val), str(val+1), "B"]))
    else:
        raise ValueError("bedline: invalid mode")
    return interval, points


def jaccard_similarity(set1, set2):
    union = float(len(set1.union(set2)))
    intersect = float(len(set1.intersection(set2)))
    return union, intersect, intersect/union


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


if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Usage: python bed_jaccards.py <filepaths file> <mode>")
        sys.exit(1)
    bed_sets = {}
    mode = sys.argv[2]
    with open(sys.argv[1], "r") as filein:
        for line in filein.readlines():
            filepath = line.strip()
            basename = os.path.basename(filepath)
            bed_sets[basename] = bed_to_sets(filename=filepath, mode=mode)
    # freeze the order of keys for bed_sets in a list
    bed_keys = list(bed_sets.keys())
    matrix = np.zeros((len(bed_keys), len(bed_keys)))
    outlines = [",".join(["bed1", "bed2", "unionA", "unionB"])]
    for i in range(len(bed_keys)):
        for j in range(i, len(bed_keys)):
            unionA, unionB = 0, 0
            intersectA, intersectB = 0, 0
            jaccA, jaccB = 0, 0
            if mode in ["A", "C"]:
                unionA, intersectA, jaccA = jaccard_similarity(bed_sets[bed_keys[i]][0], bed_sets[bed_keys[j]][0])
            if mode in ["B", "C"]:
                unionB, intersectB, jaccB = jaccard_similarity(bed_sets[bed_keys[i]][1], bed_sets[bed_keys[j]][1])
            # upper triangular matrix holds jaccard similarities for mode A
            matrix[i, j] = jaccA
            # lower triangular matrix holds jaccard similarities for mode B
            matrix[j, i] = jaccB
            set1 = bed_keys[i]
            set2 = bed_keys[j]
            if set1 == set2:
                set2 = "-"

            outlines.append(",".join([set1,
                                      set2,
                                      f"{unionA:>5f}",
                                      f"{unionB:>5f}"]))
            # print(bed_keys[i], bed_keys[j], jaccA, jaccB)
    # pretty_print_matrix(matrix,bed_keys)
    write_pretty_matrix(matrix=matrix, row_names=bed_keys)
    with open("card_out.csv", mode="w") as outfile:
        for line in outlines:
            outfile.write(line + "\n")
