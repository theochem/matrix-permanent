"""

===== Module Desciptor =====

This module contains the basic structure to allow you to generate the permanent of a
given matrix using either the combinatoric, Ryser, or Glynn algorithms. When run, if
invalid inputs are given, instructions on how to properly use it will be printed out
on the terminal. 

"""

import numpy as np

import permanent

def solve_permanent(matrix: np.ndarray, algorithm="Combinatoric") -> float:
    m_rows, n_cols = matrix.shape
    if (not m_rows) or (not n_cols):
        print("You have entered a 1x1 matrix! \n")
        print("The permanent is: ", 1)
        return matrix.dtype.type(1)
    if n_cols < m_rows:
        print("n < m.... transposing matrix! \n")
        matrix = matrix.transpose()
    if algorithm == "Combinatoric":
        print("The permanent of your matrix is: ", permanent.combinatoric(matrix))
        return permanent.combinatoric(matrix)
    elif algorithm == "Glynn":
        print("The permanent of your matrix is: ", permanent.glynn(matrix))
        return permanent.glynn(matrix)
    elif algorithm == "Ryser":
        print("The permanent of your matrix is: ", permanent.ryser(matrix))
        return permanent.ryser(matrix)
    else:
        print("Valid entries are of the form solve_permannt(matrix: np.array, algorithm: str) \n where algorithm = Combinatoric, Ryser, or Glynn. \n Check spelling!\n")
        print("Note: algorithm is an optional argument and is set to Combinatoric as default.")


if __name__ == "__main__":
    solve_permanent(np.arange(1, 5, dtype=np.double).reshape(2,2), "Rysr")