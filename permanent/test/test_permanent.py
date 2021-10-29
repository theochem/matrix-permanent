import numpy as np

import permanent

def test_2by2_permanent():
    matrix = np.arange(1, 5, dtype=np.double).reshape(2,2)
    assert permanent.combinatoric(matrix)==10

def test_3by3_permanent():
    matrix = np.arange(1, 10, dtype=np.double).reshape(3,3)
    assert permanent.combinatoric(matrix)==450

def test_4by4_permanent():
    matrix = np.arange(1, 17, dtype=np.double).reshape(4,4)
    assert permanent.combinatoric(matrix)==55456

def test_2by3_permanent():
    matrix = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.double)
    assert permanent.combinatoric(matrix)==58.0