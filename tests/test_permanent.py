from math import factorial

import numpy as np
import numpy.testing as npt
import permanent
import pytest

ATOL = 0e0


RTOL = 1e-4


FNS = [
    permanent.opt,
    permanent.combinatoric,
    permanent.glynn,
    permanent.ryser,
]


@pytest.mark.parametrize(
    "args",
    [
        (np.arange(1, 5, dtype=float).reshape(2, 2), 10),
        (np.arange(1, 5, dtype=int).reshape(2, 2), 10),
        (np.arange(1, 5, dtype=complex).reshape(2, 2), 10),
        (np.arange(1, 10, dtype=float).reshape(3, 3), 450),
        (np.arange(1, 17, dtype=float).reshape(4, 4), 55456),
        (np.arange(1, 50, dtype=float).reshape(7, 7), 5373548250000),
        # Rectangular matrices
        (np.arange(1, 7, dtype=float).reshape(2, 3), 58),
        (np.arange(1, 9, dtype=float).reshape(2, 4), 190),
        (np.arange(1, 15, dtype=float).reshape(2, 7), 1820),
        (np.arange(1, 36, dtype=float).reshape(5, 7), 1521238320),
        (np.arange(1, 43, dtype=float).reshape(6, 7), 117681979920),
        # Special matrices
        (np.ones((10, 10), dtype=float), factorial(10)),
        (np.ones((12, 12), dtype=float), factorial(12)),
        (np.identity(10, dtype=float), 1),
        (np.identity(5, dtype=float), 1),
        (np.diag(np.diag(np.ones((3, 7), dtype=float))), 1),
    ],
)
@pytest.mark.parametrize("fn", FNS)
def test_compute_permanent(fn, args):
    r"""Ensure that the correct permanent is computed."""
    npt.assert_allclose(fn(args[0]), args[1], atol=ATOL, rtol=RTOL)


@pytest.mark.parametrize(
    "arg",
    [
        np.ones((2, 65)),
        np.ones((65, 2)),
        np.ones((65, 65)),
    ],
)
@pytest.mark.parametrize("fn", FNS)
def test_dim_gt_64_raises(fn, arg):
    r"""Ensure that matrices with any dimension > 64 raise a ValueError."""
    with npt.assert_raises(ValueError):
        fn(arg)


@pytest.mark.parametrize(
    "arg",
    [
        np.ones((2, 23), dtype=int),
        np.ones((23, 2), dtype=int),
        np.ones((43, 64), dtype=int),
        np.ones((64, 43), dtype=int),
    ],
)
@pytest.mark.parametrize("fn", FNS)
def test_int_dim_diff_gt_20_raises(fn, arg):
    r"""Ensure that integer matrices with difference in dimensions > 20 raise a ValueError."""
    with npt.assert_raises(ValueError):
        fn(arg)


@pytest.mark.parametrize(
    "arg",
    [
        np.ones((4, 4), dtype=bool),
        np.array([[object() for _ in range(4)] for _ in range(4)]),
        np.array([[{"hi": 1} for _ in range(4)] for _ in range(4)]),
        np.array([[chr(x) for x in range(65, 69)] for _ in range(4)]),
    ],
)
@pytest.mark.parametrize("fn", FNS)
def test_invalid_type_raises(fn, arg):
    r"""Ensure that matrices with an invalid dtype raise a ValueError."""
    with npt.assert_raises(TypeError):
        fn(arg)
