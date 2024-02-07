import numpy as np
import math
import permanent

REL_ERROR = 0.0001;

## Test square matrices
def test_2by2_comb():
    matrix = np.arange(1, 5, dtype=np.double).reshape(2,2)
    assert abs((permanent.combinatoric(matrix) - 10) / 10) <= REL_ERROR

def test_2by2_glynn():
    matrix = np.arange(1, 5, dtype=np.double).reshape(2,2)
    assert abs((permanent.glynn(matrix) - 10) / 10) <= REL_ERROR

def test_2by2_ryser():
    matrix = np.arange(1, 5, dtype=np.double).reshape(2,2)
    assert abs((permanent.ryser(matrix) - 10) / 10) <= REL_ERROR

def test_3by3_comb():
    matrix = np.arange(1, 10, dtype=np.double).reshape(3,3)
    assert abs((permanent.combinatoric(matrix) - 450) / 450) <= REL_ERROR

def test_3by3_glynn():
    matrix = np.arange(1, 10, dtype=np.double).reshape(3,3)
    assert abs((permanent.glynn(matrix) - 450) / 450) <= REL_ERROR

def test_3by3_ryser():
    matrix = np.arange(1, 10, dtype=np.double).reshape(3,3)
    assert abs((permanent.ryser(matrix) - 450) / 450) <= REL_ERROR

def test_4by4_comb():
    matrix = np.arange(1, 17, dtype=np.double).reshape(4,4)
    assert abs((permanent.combinatoric(matrix) - 55456) / 55456) <= REL_ERROR

def test_4by4_glynn():
    matrix = np.arange(1, 17, dtype=np.double).reshape(4,4)
    assert abs((permanent.glynn(matrix) - 55456) / 55456) <= REL_ERROR

def test_4by4_ryser():
    matrix = np.arange(1, 17, dtype=np.double).reshape(4,4)
    assert abs((permanent.ryser(matrix) - 55456) / 55456) <= REL_ERROR

def test_7by7_comb():
    matrix = np.arange(1, 50, dtype=np.double).reshape(7,7)
    assert abs((permanent.combinatoric(matrix) - 5373548250000) / 5373548250000) <= REL_ERROR

def test_7by7_glynn():
    matrix = np.arange(1, 50, dtype=np.double).reshape(7,7)
    assert abs((permanent.glynn(matrix) - 5373548250000) / 5373548250000) <= REL_ERROR

def test_7by7_ryser():
    matrix = np.arange(1, 50, dtype=np.double).reshape(7,7)
    assert abs((permanent.ryser(matrix) - 5373548250000) / 5373548250000) <= REL_ERROR

## Test rectangular matrices
def test_2by3_comb():
    matrix = np.arange(1, 7, dtype=np.double).reshape(2,3)
    assert abs((permanent.combinatoric(matrix) - 58) / 58) <= REL_ERROR

def test_2by3_glynn():
    matrix = np.arange(1, 7, dtype=np.double).reshape(2,3)
    assert abs((permanent.glynn(matrix) - 58) / 58) <= REL_ERROR

def test_2by3_ryser():
    matrix = np.arange(1, 7, dtype=np.double).reshape(2,3)
    assert abs((permanent.ryser(matrix) - 58) / 58) <= REL_ERROR

def test_2by4_comb():
    matrix = np.arange(1, 9, dtype=np.double).reshape(2,4)
    assert abs((permanent.combinatoric(matrix) - 190) / 190) <= REL_ERROR

def test_2by4_glynn():
    matrix = np.arange(1, 9, dtype=np.double).reshape(2,4)
    assert abs((permanent.glynn(matrix) - 190) / 190) <= REL_ERROR

def test_2by4_ryser():
    matrix = np.arange(1, 9, dtype=np.double).reshape(2,4)
    assert abs((permanent.ryser(matrix) - 190) / 190) <= REL_ERROR

def test_2by7_comb():
    matrix = np.arange(1, 15, dtype=np.double).reshape(2,7)
    assert abs((permanent.combinatoric(matrix) - 1820) / 1820) <= REL_ERROR

def test_2by7_glynn():
    matrix = np.arange(1, 15, dtype=np.double).reshape(2,7)
    assert abs((permanent.glynn(matrix) - 1820) / 1820) <= REL_ERROR

def test_2by7_ryser():
    matrix = np.arange(1, 15, dtype=np.double).reshape(2,7)
    assert abs((permanent.ryser(matrix) - 1820) / 1820) <= REL_ERROR

def test_5by7_comb():
    matrix = np.arange(1, 36, dtype=np.double).reshape(5,7)
    assert abs((permanent.combinatoric(matrix) - 1521238320) / 1521238320) <= REL_ERROR

def test_5by7_glynn():
    matrix = np.arange(1, 36, dtype=np.double).reshape(5,7)
    assert abs((permanent.glynn(matrix) - 1521238320) / 1521238320) <= REL_ERROR

def test_5by7_ryser():
    matrix = np.arange(1, 36, dtype=np.double).reshape(5,7)
    assert abs((permanent.ryser(matrix) - 1521238320) / 1521238320) <= REL_ERROR

def test_6by7_comb():
    matrix = np.arange(1, 43, dtype=np.double).reshape(6,7)
    assert abs((permanent.combinatoric(matrix) - 117681979920) / 117681979920) <= REL_ERROR

def test_6by7_glynn():
    matrix = np.arange(1, 43, dtype=np.double).reshape(6,7)
    assert abs((permanent.glynn(matrix) - 117681979920) / 117681979920) <= REL_ERROR

def test_6by7_ryser():
    matrix = np.arange(1, 43, dtype=np.double).reshape(6,7)
    assert abs((permanent.ryser(matrix) - 117681979920) / 117681979920) <= REL_ERROR

def test_ones_comb():
    matrix = np.ones((10,10), dtype=np.double)
    assert abs((permanent.combinatoric(matrix) - math.factorial(10)) / math.factorial(10)) <= REL_ERROR

def test_ones_ryser():
    matrix = np.ones((10,10), dtype=np.double)
    assert abs((permanent.ryser(matrix) - math.factorial(10)) / math.factorial(10)) <= REL_ERROR

def test_ones_glynn():
    matrix = np.ones((10,10), dtype=np.double)
    assert abs((permanent.glynn(matrix) - math.factorial(10)) / math.factorial(10)) <= REL_ERROR

def test_ones_comb_big():
    matrix = np.ones((12,12), dtype=np.double)
    assert abs((permanent.combinatoric(matrix) - math.factorial(12)) / math.factorial(12)) <= REL_ERROR

def test_ones_ryser_big():
    matrix = np.ones((12,12), dtype=np.double)
    assert abs((permanent.ryser(matrix) - math.factorial(12)) / math.factorial(12)) <= REL_ERROR

def test_ones_glynn_big():
    matrix = np.ones((12,12), dtype=np.double)
    assert abs((permanent.glynn(matrix) - math.factorial(12)) / math.factorial(12)) <= REL_ERROR

def test_identity_comb():
    matrix = np.identity(10, dtype=np.double)
    assert abs((permanent.combinatoric(matrix) - 1) / 1) <= REL_ERROR

def test_identity_ryser():
    matrix = np.identity(10, dtype=np.double)
    assert abs((permanent.ryser(matrix) - 1) / 1) <= REL_ERROR

def test_identity_glynn():
    matrix = np.identity(10, dtype=np.double)
    assert abs((permanent.glynn(matrix) - 1) / 1) <= REL_ERROR

def test_identity_ryser_odd():
    matrix = np.identity(5, dtype=np.double)
    assert abs((permanent.ryser(matrix) - 1) / 1) <= REL_ERROR

def test_identity_glynn_odd():
    matrix = np.identity(5, dtype=np.double)
    assert abs((permanent.glynn(matrix) - 1) / 1) <= REL_ERROR

def test_identity_comb_diag():
    matrix = np.ones((3,7), dtype=np.double)
    diag_matrix = np.diag(np.diag(matrix))
    assert abs((permanent.combinatoric(diag_matrix) - 1) / 1) <= REL_ERROR

def test_identity_ryser_diag():
    matrix = np.ones((3,7), dtype=np.double)
    diag_matrix = np.diag(np.diag(matrix))
    assert abs((permanent.ryser(diag_matrix) - 1) / 1) <= REL_ERROR

def test_identity_glynn_diag():
    matrix = np.ones((3,7), dtype=np.double)
    diag_matrix = np.diag(np.diag(matrix))
    assert abs((permanent.glynn(diag_matrix) - 1) / 1) <= REL_ERROR
