import numpy as np

import permanent

def test_2by2_comb():
    matrix = np.arange(1, 5, dtype=np.double).reshape(2,2)
    assert permanent.combinatoric(matrix)==10

def test_3by3_comb():
    matrix = np.arange(1, 10, dtype=np.double).reshape(3,3)
    assert permanent.combinatoric(matrix)==450

def test_4by4_comb():
    matrix = np.arange(1, 17, dtype=np.double).reshape(4,4)
    assert permanent.combinatoric(matrix)==55456

def test_2by3_comb():
    matrix = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.double)
    assert permanent.combinatoric(matrix)==58.0

def test_2by2_glynn():
    matrix = np.arange(1, 5, dtype=np.double).reshape(2,2)
    assert permanent.combinatoric(matrix)==10

def test_3by3_glynn():
    matrix = np.arange(1, 10, dtype=np.double).reshape(3,3)
    assert permanent.ryser(matrix)==450

def test_4by4_glynn():
    matrix = np.arange(1, 17, dtype=np.double).reshape(4,4)
    assert permanent.ryser(matrix)==55456

def test_2by2_ryser():
    matrix = np.arange(1, 5, dtype=np.double).reshape(2,2)
    assert permanent.combinatoric(matrix)==10

def test_3by3_ryser():
    matrix = np.arange(1, 10, dtype=np.double).reshape(3,3)
    assert permanent.ryser(matrix)==450

def test_4by4_ryser():
    matrix = np.arange(1, 17, dtype=np.double).reshape(4,4)
    assert permanent.ryser(matrix)==55456

def test_2by3_ryser():
    matrix = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.double)
    assert permanent.ryser(matrix)==58.0