import numpy as np

import permanent

def test_3by3_permanent():

    matrix = np.arange(1, 10, dtype = np.double).reshape(3,3)
    assert permanent.combinatoric(matrix) == 450

def test_4by4_permanent():

    matrix = np.arange(1, 17, dtype = np.double).reshape(4,4)
    assert permanent.combinatoric(matrix) == 55456
