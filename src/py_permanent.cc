#include <cstdlib>
#include <cstdint>

#include <complex>

#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/ndarraytypes.h>

#include "permanent.h"


static const char DOCSTRING_MODULE[] = R"""(
Permanent module C extension.

)""";


static const char DOCSTRING_PERMANENT[] = R"""(
Compute the permanent of a matrix using the best algorithm.

Parameters
----------
matrix : np.ndarray(M, N, dtype=(np.double|np.complex))

Returns
-------
permanent : (np.double|np.complex)
    Permanent of matrix.

)""";


static const char DOCSTRING_COMBINATORIC[] = R"""(
Compute the permanent of a matrix combinatorically.

.. math::

    \text{per}(A) = \sum_{\sigma \in P(N,M)}{
        \prod_{i=1}^M{a_{i,{\sigma(i)}}}
    }

Parameters
----------
matrix : np.ndarray(M, N, dtype=(np.double|np.complex))

Returns
-------
permanent : (np.double|np.complex)
    Permanent of matrix.

)""";


static const char DOCSTRING_GLYNN[] = R"""(
Compute the permanent of a matrix via Glynn's algorithm [Glynn]_.

.. math::

    \text{per}(A) = \frac{1}{2^{N-1}} \cdot \sum_{\delta}{
        \left(\sum_{k=1}^N{\delta_k}\right)
        \prod_{j=1}^N{\sum_{i=1}^N{\delta_i a_{i,j}}}
    }

The original formula has been generalized here to work with
:math:`M`-by-:math:`N` rectangular permanents with :math:`M \leq N`
by use of the following identity (shown here for :math:`M \geq N`):

.. math::

    {\text{per}}\left(
        \begin{matrix}
            a_{1,1} & \cdots & a_{1,N} \\
            \vdots & \ddots & \vdots \\
            a_{M,1} & \cdots & a_{M,N} \\
        \end{matrix}
    \right)
    = \frac{1}{\left(M - N + 1\right)!} \cdot {\text{per}}\left(
        \begin{matrix}
            a_{1,1} & \cdots & a_{1,N} & 1_{1,N+1} & \cdots & 1_{1,M} \\
            \vdots & \ddots & \vdots & \vdots & \ddots & \vdots \\
            a_{M,1} & \cdots & a_{M,N} & 1_{M,N+1} & \cdots & 1_{M,M} \\
        \end{matrix}
    \right)

This can be neatly fit into the original formula by extending the inner sums
over :math:`\delta` from :math:`\left[1,M\right]` to :math:`\left[1,N\right]`:

.. math::

    \text{per}(A) = \frac{1}{2^{N-1}} \cdot \frac{1}{\left(N - M + 1\right)!}
    \cdot \sum_{\delta}{
        \left(\sum_{k=1}^N{\delta_k}\right)
        \prod_{j=1}^N{\left(
            \sum_{i=1}^M{\delta_i a_{i,j}} + \sum_{i=M+1}^N{\delta_i}
        \right)}
    }

.. [Glynn] Glynn, D. G. (2010). The permanent of a square matrix.
           *European Journal of Combinatorics*, 31(7), 1887-1891.

Parameters
----------
matrix : np.ndarray(M, N, dtype=(np.double|np.complex))

Returns
-------
permanent : (np.double|np.complex)
    Permanent of matrix.

)""";


static const char DOCSTRING_RYSER[] = R"""(
Compute the permanent of a matrix via Ryser's algorithm [Ryser]_.

.. math::

    \text{per}(A) = \sum_{k=0}^{M-1}{
        {\left(-1\right)}^k
        \left(\begin{matrix}N - M + k\\ k\end{matrix}\right)
        \sum_{\sigma \in P(N,M-k)}{
            \prod_{i=1}^M{
                \sum_{j=1}^{M-k}{a_{i,{\sigma(j)}}}
            }
        }
    }

.. [Ryser] Ryser, H. J. (1963). *Combinatorial Mathematics* (Vol. 14).
           American Mathematical Soc..

Parameters
----------
matrix : np.ndarray(M, N, dtype=(np.double|np.complex))

Returns
-------
permanent : (np.double|np.complex)
    Permanent of matrix.

)""";


static PyObject *py_opt(PyObject *module, PyObject *object)
{
    (void)module;

    PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, nullptr, 2, 2, NPY_ARRAY_CARRAY_RO, nullptr);
    if (!matrix) {
        PyErr_SetString(PyExc_TypeError, "Argument must be a 2-dimensional array");
        return nullptr;
    }

    std::size_t m = PyArray_DIMS(matrix)[0];
    std::size_t n = PyArray_DIMS(matrix)[1];
    if (m > n) {
        return py_opt(module, PyArray_Transpose(matrix, nullptr));
    } else if (m > 64 || n > 64) {
        PyErr_SetString(PyExc_ValueError, "Argument array must have <=64 rows and columns");
        return nullptr;
    }

    int type = PyArray_TYPE(matrix);
    if (type == NPY_DOUBLE) {
        double *ptr = (double *)PyArray_GETPTR2(matrix, 0, 0);
        return PyFloat_FromDouble(opt<double>(m, n, ptr));
    }
    else if (type == NPY_COMPLEX128) {
        std::complex<double> *ptr = (std::complex<double> *)PyArray_GETPTR2(matrix, 0, 0);
        std::complex<double> out = opt<std::complex<double>>(m, n, ptr);
        return PyComplex_FromDoubles(out.real(), out.imag());
    } else {
        PyErr_SetString(PyExc_TypeError, "Array must have dtype (double|complex)");
        return nullptr;
    }
}


static PyObject *py_combinatoric(PyObject *module, PyObject *object)
{
    (void)module;

    PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, nullptr, 2, 2, NPY_ARRAY_CARRAY_RO, nullptr);
    if (!matrix) {
        PyErr_SetString(PyExc_TypeError, "Argument must be a 2-dimensional array");
        return nullptr;
    }

    std::size_t m = PyArray_DIMS(matrix)[0];
    std::size_t n = PyArray_DIMS(matrix)[1];
    if (m > n) {
        return py_combinatoric(module, PyArray_Transpose(matrix, nullptr));
    }

    int type = PyArray_TYPE(matrix);
    if (type == NPY_DOUBLE) {
        double *ptr = (double *)PyArray_GETPTR2(matrix, 0, 0);
        return PyFloat_FromDouble(
            (m == n) ? combinatoric<double>(m, n, ptr)
                     : combinatoric_rectangular<double>(m, n, ptr)
        );
    }
    else if (type == NPY_COMPLEX128) {
        std::complex<double> *ptr = (std::complex<double> *)PyArray_GETPTR2(matrix, 0, 0);
        std::complex<double> out =
            (m == n) ? combinatoric<std::complex<double>>(m, n, ptr)
                     : combinatoric_rectangular<std::complex<double>>(m, n, ptr);
        return PyComplex_FromDoubles(out.real(), out.imag());
    } else {
        PyErr_SetString(PyExc_TypeError, "Array must have dtype (double|complex)");
        return nullptr;
    }
}


static PyObject *py_glynn(PyObject *module, PyObject *object)
{
    (void)module;

    PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, nullptr, 2, 2, NPY_ARRAY_CARRAY_RO, nullptr);
    if (!matrix) {
        PyErr_SetString(PyExc_TypeError, "Argument must be a 2-dimensional array");
        return nullptr;
    }

    std::size_t m = PyArray_DIMS(matrix)[0];
    std::size_t n = PyArray_DIMS(matrix)[1];
    if (m > n) {
        return py_glynn(module, PyArray_Transpose(matrix, nullptr));
    }

    int type = PyArray_TYPE(matrix);
    if (type == NPY_DOUBLE) {
        double *ptr = (double *)PyArray_GETPTR2(matrix, 0, 0);
        return PyFloat_FromDouble(
            (m == n) ? glynn<double>(m, n, ptr)
                     : glynn_rectangular<double>(m, n, ptr)
        );
    }
    else if (type == NPY_COMPLEX128) {
        std::complex<double> *ptr = (std::complex<double> *)PyArray_GETPTR2(matrix, 0, 0);
        std::complex<double> out =
            (m == n) ? glynn<std::complex<double>>(m, n, ptr)
                     : glynn_rectangular<std::complex<double>>(m, n, ptr);
        return PyComplex_FromDoubles(out.real(), out.imag());
    } else {
        PyErr_SetString(PyExc_TypeError, "Array must have dtype (double|complex)");
        return nullptr;
    }
}


static PyObject *py_ryser(PyObject *module, PyObject *object)
{
    (void)module;

    PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, nullptr, 2, 2, NPY_ARRAY_CARRAY_RO, nullptr);
    if (!matrix) {
        PyErr_SetString(PyExc_TypeError, "Argument must be a 2-dimensional array");
        return nullptr;
    }

    std::size_t m = PyArray_DIMS(matrix)[0];
    std::size_t n = PyArray_DIMS(matrix)[1];
    if (m > n) {
        return py_ryser(module, PyArray_Transpose(matrix, nullptr));
    }

    int type = PyArray_TYPE(matrix);
    if (type == NPY_DOUBLE) {
        double *ptr = (double *)PyArray_GETPTR2(matrix, 0, 0);
        return PyFloat_FromDouble(
            (m == n) ? ryser<double>(m, n, ptr)
                     : ryser_rectangular<double>(m, n, ptr)
        );
    }
    else if (type == NPY_COMPLEX128) {
        std::complex<double> *ptr = (std::complex<double> *)PyArray_GETPTR2(matrix, 0, 0);
        std::complex<double> out =
            (m == n) ? ryser<std::complex<double>>(m, n, ptr)
                     : ryser_rectangular<std::complex<double>>(m, n, ptr);
        return PyComplex_FromDoubles(out.real(), out.imag());
    } else {
        PyErr_SetString(PyExc_TypeError, "Array must have dtype (double|complex)");
        return nullptr;
    }
}


/* Define the Python methods that will go into the C extension module.       *
 *                                                                           *
 * Note:  METH_O indicates that the Python function takes a single argument. *
 *        On the C side, the function takes two PyObject* arguments;         *
 *        the first one is the C extension module itself,                    *
 *        and the second one is the argument to the Python function.         */

static PyMethodDef methods[] = {
    /* Python function name     C function          Args flag   Docstring */
    { "opt",                    py_opt,             METH_O,     DOCSTRING_PERMANENT },
    { "combinatoric",           py_combinatoric,    METH_O,     DOCSTRING_COMBINATORIC },
    { "glynn",                  py_glynn,           METH_O,     DOCSTRING_GLYNN },
    { "ryser",                  py_ryser,           METH_O,     DOCSTRING_RYSER },
    { nullptr,                  nullptr,            0,          nullptr } /* sentinel value */
};


/* Define the C extension module. */

static struct PyModuleDef definition = {
    PyModuleDef_HEAD_INIT, "permanent", DOCSTRING_MODULE, -1, methods, nullptr, nullptr, nullptr, nullptr
};


/* Initialize the C extension module. */

PyMODINIT_FUNC PyInit_permanent(void) {
    Py_Initialize();                      /* Initialize Python API */
    import_array();                       /* Initialize NumPy NDArray API */
    return PyModule_Create(&definition);  /* Create module. */
}
