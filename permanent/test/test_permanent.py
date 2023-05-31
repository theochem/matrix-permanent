import numpy as np
import math

import permanent

## Test square matrices
def test_2by2_comb():
    matrix = np.arange(1, 5, dtype=np.double).reshape(2,2)
    assert permanent.combinatoric(matrix)==10

def test_2by2_glynn():
    matrix = np.arange(1, 5, dtype=np.double).reshape(2,2)
    assert permanent.combinatoric(matrix)==10

def test_2by2_ryser():
    matrix = np.arange(1, 5, dtype=np.double).reshape(2,2)
    assert permanent.ryser(matrix)==10

def test_3by3_comb():
    matrix = np.arange(1, 10, dtype=np.double).reshape(3,3)
    assert permanent.combinatoric(matrix)==450

def test_3by3_glynn():
    matrix = np.arange(1, 10, dtype=np.double).reshape(3,3)
    assert permanent.ryser(matrix)==450

def test_3by3_ryser():
    matrix = np.arange(1, 10, dtype=np.double).reshape(3,3)
    assert permanent.ryser(matrix)==450

def test_4by4_comb():
    matrix = np.arange(1, 17, dtype=np.double).reshape(4,4)
    assert permanent.combinatoric(matrix)==55456

def test_4by4_glynn():
    matrix = np.arange(1, 17, dtype=np.double).reshape(4,4)
    assert permanent.ryser(matrix)==55456

def test_4by4_ryser():
    matrix = np.arange(1, 17, dtype=np.double).reshape(4,4)
    assert permanent.ryser(matrix)==55456

def test_7by7_comb():
    matrix = np.arange(1, 50, dtype=np.double).reshape(7,7)
    assert permanent.combinatoric(matrix)==5373548250000

def test_7by7_glynn():
    matrix = np.arange(1, 50, dtype=np.double).reshape(7,7)
    assert permanent.ryser(matrix)==5373548250000

def test_7by7_ryser():
    matrix = np.arange(1, 50, dtype=np.double).reshape(7,7)
    assert permanent.ryser(matrix)==5373548250000

## Test rectangular matrices
def test_2by3_comb():
    matrix = np.arange(1, 7, dtype=np.double).reshape(2,3)
    assert permanent.combinatoric(matrix)==58.0

def test_2by3_glynn():
    matrix = np.arange(1, 7, dtype=np.double).reshape(2,3)
    assert permanent.glynn(matrix)==58.0

def test_2by3_ryser():
    matrix = np.arange(1, 7, dtype=np.double).reshape(2,3)
    assert permanent.ryser(matrix)==58.0

def test_2by4_comb():
    matrix = np.arange(1, 9, dtype=np.double).reshape(2,4)
    assert permanent.combinatoric(matrix)==190

def test_2by4_glynn():
    matrix = np.arange(1, 9, dtype=np.double).reshape(2,4)
    assert permanent.glynn(matrix)==190

def test_2by4_ryser():
    matrix = np.arange(1, 9, dtype=np.double).reshape(2,4)
    assert permanent.ryser(matrix)==190

def test_2by7_comb():
    matrix = np.arange(1, 15, dtype=np.double).reshape(2,7)
    assert permanent.combinatoric(matrix)==1820

def test_2by7_glynn():
    matrix = np.arange(1, 15, dtype=np.double).reshape(2,7)
    assert permanent.glynn(matrix)==1820

def test_2by7_ryser():
    matrix = np.arange(1, 15, dtype=np.double).reshape(2,7)
    assert permanent.ryser(matrix)==1820

def test_5by7_comb():
    matrix = np.arange(1, 36, dtype=np.double).reshape(5,7)
    assert permanent.combinatoric(matrix)==1521238320

def test_5by7_glynn():
    matrix = np.arange(1, 36, dtype=np.double).reshape(5,7)
    assert permanent.glynn(matrix)==1521238320

def test_5by7_ryser():
    matrix = np.arange(1, 36, dtype=np.double).reshape(5,7)
    assert permanent.ryser(matrix)==1521238320

def test_6by7_comb():
    matrix = np.arange(1, 43, dtype=np.double).reshape(6,7)
    assert permanent.combinatoric(matrix)==117681979920

def test_6by7_glynn():
    matrix = np.arange(1, 43, dtype=np.double).reshape(6,7)
    assert permanent.glynn(matrix)==117681979920

def test_6by7_ryser():
    matrix = np.arange(1, 43, dtype=np.double).reshape(6,7)
    assert permanent.ryser(matrix)==117681979920

def test_ones_comb():
    matrix = np.ones((10,10), dtype=np.double)
    assert permanent.combinatoric(matrix)==math.factorial(10)

def test_ones_ryser():
    matrix = np.ones((10,10), dtype=np.double)
    assert permanent.ryser(matrix)==math.factorial(10)

def test_ones_glynn():
    matrix = np.ones((10,10), dtype=np.double)
    assert permanent.glynn(matrix)==math.factorial(10)

def test_ones_comb_big():
    matrix = np.ones((12,12), dtype=np.double)
    assert permanent.combinatoric(matrix)==math.factorial(12)

def test_ones_ryser_big():
    matrix = np.ones((12,12), dtype=np.double)
    assert permanent.ryser(matrix)==math.factorial(12)

def test_ones_glynn_big():
    matrix = np.ones((12,12), dtype=np.double)
    assert permanent.glynn(matrix)==math.factorial(12)

def test_identity_comb():
    matrix = np.identity(10, dtype=np.double)
    assert permanent.combinatoric==1.0

def test_identity_ryser():
    matrix = np.identity(10, dtype=np.double)
    assert permanent.ryser==1.0

def test_identity_glynn():
    matrix = np.identity(10, dtype=np.double)
    assert permanent.glynn==1.0

def test_identity_ryser_odd():
    matrix = np.identity(5, dtype=np.double)
    assert permanent.ryser==1.0

def test_identity_glynn_odd():
    matrix = np.identity(5, dtype=np.double)
    assert permanent.glynn==1.0
