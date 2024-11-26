#!/usr/bin/env python
import sys
import numpy as np

def basic_bedline(line):
    """Parse a single line from a BED file into chromosome, start, and end coordinates.
    
    Args:
        line: A string containing a single line from a BED file
        
    Returns:
        Tuple of (chrom, start, end) where:
            chrom: Chromosome name with 'chr' prefix removed if present 
            start: Integer start coordinate
            end: Integer end coordinate
            
    Raises:
        ValueError: If line has fewer than 3 tab or space-separated columns
    """
    columns = line.strip().split('\t')
    if len(columns) < 2:
        columns = line.strip().split(" ")
        if len(columns) < 2:
            raise ValueError("bedline: one of the lines in malformed")
    if columns[0].startswith('chr'):
        columns[0] = columns[0][3:]
    return columns[0], int(columns[1]), int(columns[2])

def read_csv_matrix(file):
    """
    Read a matrix from a CSV file.

    Args:
        file (str): The path to the CSV file containing the matrix.

    Returns:
        np.ndarray: The matrix read from the CSV file.
    """
    with open (file, mode="r") as fh:
        data = np.genfromtxt(fh, delimiter=',', skip_header=1)
        row_names = np.genfromtxt(fh, delimiter=',', dtype=str, usecols=0, skip_header=1)
        col_names = np.genfromtxt(fh, delimiter=',', dtype=str, usecols=range(1, data.shape[1]+1), skip_header=1)
   
    return data, row_names, col_names

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
