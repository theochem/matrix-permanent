import numpy as np

import permanent


def test_placeholder():

    matrix = np.zeros((3, 3), dtype=np.double)
    assert permanent.combinatoric(matrix) == 0.

    matrix = np.ones((3, 3), dtype=np.double)
    assert permanent.combinatoric(matrix) == 1.
