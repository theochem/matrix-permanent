r"""Permanent algorithms in Python."""


__all__ = [
    "comb",
    "glynn",
]


import math

import itertools as it

import numpy as np


def comb(matrix: np.ndarray) -> float:
    r"""
    Compute the permanent of a matrix using the definition of the permanent.

    Parameters
    ----------
    matrix : np.ndarray
        Square matrix.

    Returns
    -------
    result : matrix.dtype
        Permanent of the matrix.

    """
    m, n = matrix.shape
    if (not m) or (not n):
        return 1
    if n < m:
        return comb(matrix.T)
    rows = np.arange(n)
    return sum(np.prod(matrix[rows, cols]) for cols in it.permutations(rows, m))


def glynn(matrix: np.ndarray) -> float:
    r"""
    Compute the permanent of a matrix using Glynn's algorithm.

    Gray code generation from Knuth, D. E. (2005). The Art of Computer Programming,
    Volume 4, Fascicle 2: Generating All Tuples and Permutations.

    Glynn's algorithm from Glynn, D. G. (2010). The permanent of a square matrix.
    European Journal of Combinatorics, 31(7), 1887-1891.

    Parameters
    ----------
    matrix : np.ndarray
        Square matrix.

    Returns
    -------
    result : matrix.dtype
        Permanent of the matrix.

    """
    # Permanent of zero-by-zero matrix is 1
    m, n = matrix.shape
    if (not m) or (not n):
        return 1
    if n < m:
        return glynn(matrix.T)

    # Initialize gray code
    pos = 0
    sign = 1
    bound = n - 1
    delta = np.ones(n, dtype=int)
    graycode = np.arange(n, dtype=int)

    # Iterate over every delta
    result = np.prod(np.sum(matrix, axis=0))
    while pos < bound:
        # Update delta and add term to permanent
        sign *= -1
        delta[bound - pos] *= -1
        result += sign * np.prod(delta[:m].dot(matrix) + np.sum(delta[m:]))
        # Update gray code and position
        graycode[0] = 0
        graycode[pos] = graycode[pos + 1]
        graycode[pos + 1] = pos + 1
        pos = graycode[0]

    # Divide by constant factor
    return result / (math.factorial(n - m + 1) * 2 ** bound)
