#!/usr/bin/env python
import sys
import numpy as np

def write_pretty_matrix(matrix:np.ndarray, row_names, outhandle=sys.stdout):
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

def basic_bedline(line):
    columns = line.strip().split('\t')
    if 0 < len(columns) <= 2:
        raise ValueError("bedline: one of the lines in malformed")
    if columns[0].startswith('chr'):
        columns[0] = columns[0][3:]
    return int(columns[0]), int(columns[1]), int(columns[2])