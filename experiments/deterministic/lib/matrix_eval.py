#!/usr/bin/env python
import sys
import numpy as np
from bedprocessing import write_pretty_matrix, basic_bedline

# def write_pretty_matrix(matrix:np.ndarray, row_names, outhandle=sys.stdout):
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

# def basic_bedline(line):
#     columns = line.strip().split('\t')
#     if 0 < len(columns) <= 2:
#         raise ValueError("bedline: one of the lines in malformed")
#     if columns[0].startswith('chr'):
#         columns[0] = columns[0][3:]
#     return int(columns[0]), int(columns[1]), int(columns[2])

def matrix_compare(mat1,mat2):
    """
    Compare two matrices, print their difference matrix and compute their euclidean distance.

    Args:
        mat1 (np.ndarray): The first matrix to compare.
        mat2 (np.ndarray): The second matrix to compare.

    Returns:
        the distance Matrix
        euclidean distances.
    """
    mat1 = mat1[:, 1:]
    mat2 = mat2[:, 1:]
    diff = mat1 - mat2
    distances={}
    
    distances["cosine"] = np.dot(diff.flatten(), diff.flatten()) / (np.linalg.norm(diff) * np.linalg.norm(diff))
    # dist = np.linalg.norm(diff)
    # Calculate the Frobenius distance between the matrices
    distances["frobenius"] = np.linalg.norm(diff)
    # Calculate the correlation matrix distance between the matrices
    distances["correlation"] = np.linalg.norm(np.corrcoef(mat1.T) - np.corrcoef(mat2.T))
    # Calculate the Kullback-Leibler divergence between the matrices
    distances["kl_divergence"] = np.sum(mat1 * np.log(mat1 / mat2, where=(mat1 != 0) & (mat2 != 0)))
    # Calculate the eigenvalue-based distance between the matrices
    eigenvalues_diff = np.linalg.eigvals(mat1) - np.linalg.eigvals(mat2)
    distances["eigenvalue"] = np.linalg.norm(eigenvalues_diff)

    return diff, distances

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


if __name__ == "__main__":
    # Example usage of the functions
    # matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    # row_names = ["Row 1", "Row 2", "Row 3"]
    
    # Write the matrix to a file
    mat1, row1, col1 = read_csv_matrix( sys.argv[1])
    mat2, row2, col2 = read_csv_matrix( sys.argv[2])
    
    print("Matrix 1:\n",mat1)
    print("Matrix 2:\n",mat2)
    # Read the matrix from the file
    # data, row_names, col_names = read_csv_matrix("matrix.csv")
    
    # Compare two matrices
    diff, distances = matrix_compare(mat1, mat2)
    print("Difference Matrix:")
    print(diff)
    print("Cosine Similarity:", distances["cosine"])
    print("Frobenius Distance:", distances["frobenius"])
    print("Correlation Distance:", distances["correlation"])
    print("KL Divergence:", distances["kl_divergence"])
    print("Eigenvalue Distance:", distances["eigenvalue"])