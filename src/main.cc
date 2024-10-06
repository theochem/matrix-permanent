/* Copyright 2024 QC-Devs (GPLv3) */

#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/ndarraytypes.h>
#include <permanent.h>

#include <complex>
#include <cstdlib>

#define LITERAL(S) #S
#define STRINGIZE(S) LITERAL(S)

#define PERMANENT_VERSION STRINGIZE(_PERMANENT_VERSION)
#define PERMANENT_COMPILER_VERSION STRINGIZE(_PERMANENT_COMPILER_VERSION)
#define PERMANENT_GIT_BRANCH STRINGIZE(_PERMANENT_GIT_BRANCH)
#define PERMANENT_GIT_COMMIT_HASH STRINGIZE(_PERMANENT_GIT_COMMIT_HASH)

static const char DOCSTRING_MAIN_MODULE[] = R"""(
Permanent module.

)""";

static const char DOCSTRING_VERSION_INFO_MODULE[] = R"""(
version information.

)""";

static const char DOCSTRING_PERMANENT[] = R"""(
Compute the permanent of a matrix using the best algorithm for the shape of the given matrix.

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
for
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

#define PERMANENT_DISPATCH_ARRAY_TYPE(FN, MATRIX)                                          \
  {                                                                                        \
    int type = PyArray_TYPE(MATRIX);                                                       \
    if (type == NPY_INT8) {                                                                \
      if (std::abs<ptrdiff_t>(static_cast<ptrdiff_t>(m) - n) > 20) {                       \
        PyErr_SetString(PyExc_ValueError,                                                  \
                        "Difference between # cols and # rows cannot exceed 20");          \
        return nullptr;                                                                    \
      }                                                                                    \
      const std::int8_t *ptr =                                                             \
          reinterpret_cast<const std::int8_t *>(PyArray_GETPTR2(MATRIX, 0, 0));            \
      return PyLong_FromSsize_t(FN<std::int8_t, Py_ssize_t>(m, n, ptr));                   \
    } else if (type == NPY_INT16) {                                                        \
      if (std::abs<ptrdiff_t>(static_cast<ptrdiff_t>(m) - n) > 20) {                       \
        PyErr_SetString(PyExc_ValueError,                                                  \
                        "Difference between # cols and # rows cannot exceed 20");          \
        return nullptr;                                                                    \
      }                                                                                    \
      const std::int16_t *ptr =                                                            \
          reinterpret_cast<const std::int16_t *>(PyArray_GETPTR2(MATRIX, 0, 0));           \
      return PyLong_FromSsize_t(FN<std::int16_t, Py_ssize_t>(m, n, ptr));                  \
    } else if (type == NPY_INT32) {                                                        \
      if (std::abs<ptrdiff_t>(static_cast<ptrdiff_t>(m) - n) > 20) {                       \
        PyErr_SetString(PyExc_ValueError,                                                  \
                        "Difference between # cols and # rows cannot exceed 20");          \
        return nullptr;                                                                    \
      }                                                                                    \
      const std::int32_t *ptr =                                                            \
          reinterpret_cast<const std::int32_t *>(PyArray_GETPTR2(MATRIX, 0, 0));           \
      return PyLong_FromSsize_t(FN<std::int32_t, Py_ssize_t>(m, n, ptr));                  \
    } else if (type == NPY_INT64) {                                                        \
      if (std::abs<ptrdiff_t>(static_cast<ptrdiff_t>(m) - n) > 20) {                       \
        PyErr_SetString(PyExc_ValueError,                                                  \
                        "Difference between # cols and # rows cannot exceed 20");          \
        return nullptr;                                                                    \
      }                                                                                    \
      const std::int64_t *ptr =                                                            \
          reinterpret_cast<const std::int64_t *>(PyArray_GETPTR2(MATRIX, 0, 0));           \
      return PyLong_FromSsize_t(FN<std::int64_t, Py_ssize_t>(m, n, ptr));                  \
    } else if (type == NPY_UINT8) {                                                        \
      if (std::abs<ptrdiff_t>(static_cast<ptrdiff_t>(m) - n) > 20) {                       \
        PyErr_SetString(PyExc_ValueError,                                                  \
                        "Difference between # cols and # rows cannot exceed 20");          \
        return nullptr;                                                                    \
      }                                                                                    \
      const std::uint8_t *ptr =                                                            \
          reinterpret_cast<const std::uint8_t *>(PyArray_GETPTR2(MATRIX, 0, 0));           \
      return PyLong_FromSsize_t(FN<std::uint8_t, Py_ssize_t>(m, n, ptr));                  \
    } else if (type == NPY_UINT16) {                                                       \
      if (std::abs<ptrdiff_t>(static_cast<ptrdiff_t>(m) - n) > 20) {                       \
        PyErr_SetString(PyExc_ValueError,                                                  \
                        "Difference between # cols and # rows cannot exceed 20");          \
        return nullptr;                                                                    \
      }                                                                                    \
      const std::uint16_t *ptr =                                                           \
          reinterpret_cast<const std::uint16_t *>(PyArray_GETPTR2(MATRIX, 0, 0));          \
      return PyLong_FromSsize_t(FN<std::uint16_t, Py_ssize_t>(m, n, ptr));                 \
    } else if (type == NPY_UINT32) {                                                       \
      if (std::abs<ptrdiff_t>(static_cast<ptrdiff_t>(m) - n) > 20) {                       \
        PyErr_SetString(PyExc_ValueError,                                                  \
                        "Difference between # cols and # rows cannot exceed 20");          \
        return nullptr;                                                                    \
      }                                                                                    \
      const std::uint32_t *ptr =                                                           \
          reinterpret_cast<const std::uint32_t *>(PyArray_GETPTR2(MATRIX, 0, 0));          \
      return PyLong_FromSsize_t(FN<std::uint32_t, Py_ssize_t>(m, n, ptr));                 \
    } else if (type == NPY_UINT64) {                                                       \
      if (std::abs<ptrdiff_t>(static_cast<ptrdiff_t>(m) - n) > 20) {                       \
        PyErr_SetString(PyExc_ValueError,                                                  \
                        "Difference between # cols and # rows cannot exceed 20");          \
        return nullptr;                                                                    \
      }                                                                                    \
      const std::uint64_t *ptr =                                                           \
          reinterpret_cast<const std::uint64_t *>(PyArray_GETPTR2(MATRIX, 0, 0));          \
      return PyLong_FromSsize_t(FN<std::uint64_t, Py_ssize_t>(m, n, ptr));                 \
    } else if (type == NPY_FLOAT) {                                                        \
      const float *ptr = reinterpret_cast<const float *>(PyArray_GETPTR2(MATRIX, 0, 0));   \
      return PyFloat_FromDouble(FN<float>(m, n, ptr));                                     \
    } else if (type == NPY_DOUBLE) {                                                       \
      const double *ptr = reinterpret_cast<const double *>(PyArray_GETPTR2(MATRIX, 0, 0)); \
      return PyFloat_FromDouble(FN<double>(m, n, ptr));                                    \
    } else if (type == NPY_COMPLEX64) {                                                    \
      std::complex<float> *ptr =                                                           \
          reinterpret_cast<std::complex<float> *>(PyArray_GETPTR2(MATRIX, 0, 0));          \
      std::complex<double> out = FN<std::complex<float>>(m, n, ptr);                       \
      return PyComplex_FromDoubles(out.real(), out.imag());                                \
    } else if (type == NPY_COMPLEX128) {                                                   \
      std::complex<double> *ptr =                                                          \
          reinterpret_cast<std::complex<double> *>(PyArray_GETPTR2(MATRIX, 0, 0));         \
      std::complex<double> out = FN<std::complex<double>>(m, n, ptr);                      \
      return PyComplex_FromDoubles(out.real(), out.imag());                                \
    } else {                                                                               \
      PyErr_SetString(PyExc_TypeError, "Array has unsupported dtype");                     \
      return nullptr;                                                                      \
    }                                                                                      \
  }

static PyObject *py_opt(PyObject *module, PyObject *object)
{
  (void)module;

  PyArrayObject *matrix = reinterpret_cast<PyArrayObject *>(
      PyArray_FromAny(object, nullptr, 2, 2, NPY_ARRAY_CARRAY_RO, nullptr));
  if (!matrix) {
    PyErr_SetString(PyExc_TypeError, "Argument must be a 2-dimensional array");
    return nullptr;
  }

  size_t m = PyArray_DIMS(matrix)[0];
  size_t n = PyArray_DIMS(matrix)[1];
  if (m > n) {
    return py_opt(module, PyArray_Transpose(matrix, nullptr));
  } else if (m > 64 || n > 64) {
    PyErr_SetString(PyExc_ValueError, "Argument array must have <=64 rows and columns");
    return nullptr;
  }

  PERMANENT_DISPATCH_ARRAY_TYPE(permanent::opt, matrix)
}

static PyObject *py_combinatoric(PyObject *module, PyObject *object)
{
  (void)module;

  PyArrayObject *matrix = reinterpret_cast<PyArrayObject *>(
      PyArray_FromAny(object, nullptr, 2, 2, NPY_ARRAY_CARRAY_RO, nullptr));
  if (!matrix) {
    PyErr_SetString(PyExc_TypeError, "Argument must be a 2-dimensional array");
    return nullptr;
  }

  size_t m = PyArray_DIMS(matrix)[0];
  size_t n = PyArray_DIMS(matrix)[1];
  if (m > n) {
    return py_combinatoric(module, PyArray_Transpose(matrix, nullptr));
  }

  PERMANENT_DISPATCH_ARRAY_TYPE(permanent::combinatoric, matrix)
}

static PyObject *py_glynn(PyObject *module, PyObject *object)
{
  (void)module;

  PyArrayObject *matrix = reinterpret_cast<PyArrayObject *>(
      PyArray_FromAny(object, nullptr, 2, 2, NPY_ARRAY_CARRAY_RO, nullptr));
  if (!matrix) {
    PyErr_SetString(PyExc_TypeError, "Argument must be a 2-dimensional array");
    return nullptr;
  }

  size_t m = PyArray_DIMS(matrix)[0];
  size_t n = PyArray_DIMS(matrix)[1];
  if (m > n) {
    return py_glynn(module, PyArray_Transpose(matrix, nullptr));
  }

  PERMANENT_DISPATCH_ARRAY_TYPE(permanent::glynn, matrix)
}

static PyObject *py_ryser(PyObject *module, PyObject *object)
{
  PyArrayObject *matrix = reinterpret_cast<PyArrayObject *>(
      PyArray_FromAny(object, nullptr, 2, 2, NPY_ARRAY_CARRAY_RO, nullptr));
  if (!matrix) {
    PyErr_SetString(PyExc_TypeError, "Argument must be a 2-dimensional array");
    return nullptr;
  }

  size_t m = PyArray_DIMS(matrix)[0];
  size_t n = PyArray_DIMS(matrix)[1];
  if (m > n) {
    return py_ryser(module, PyArray_Transpose(matrix, nullptr));
  }

  PERMANENT_DISPATCH_ARRAY_TYPE(permanent::ryser, matrix)
}

/* Define the Python methods that will go into the C extension module.       *
 *                                                                           *
 * Note:  METH_O indicates that the Python function takes a single argument. *
 *        On the C side, the function takes two PyObject* arguments;         *
 *        the first one is the C extension module itself,                    *
 *        and the second one is the argument to the Python function.         */

static PyMethodDef main_methods[] = {
    /* Python function name     C function          Args flag   Docstring */
    {"opt", py_opt, METH_O, DOCSTRING_PERMANENT},
    {"combinatoric", py_combinatoric, METH_O, DOCSTRING_COMBINATORIC},
    {"glynn", py_glynn, METH_O, DOCSTRING_GLYNN},
    {"ryser", py_ryser, METH_O, DOCSTRING_RYSER},
    {nullptr, nullptr, 0, nullptr} /* sentinel value */
};

static PyMethodDef version_info_methods[] = {
    {nullptr, nullptr, 0, nullptr} /* sentinel value */
};

/* Define the C extension module. */

static struct PyModuleDef main_def = {
    PyModuleDef_HEAD_INIT,
    "permanent",
    DOCSTRING_MAIN_MODULE,
    -1,
    main_methods,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
};

static struct PyModuleDef version_info_def = {
    PyModuleDef_HEAD_INIT,
    "version_info",
    DOCSTRING_VERSION_INFO_MODULE,
    -1,
    version_info_methods,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
};

/* Initialize the C extension module. */

// cppcheck-suppress unusedFunction
PyMODINIT_FUNC PyInit_permanent(void)
{
  /* Initialize Python API */
  Py_Initialize();

  /* Initialize NumPy NDArray API */
  import_array();

  /* Create version info module. */
  auto *version_info_module = PyModule_Create(&version_info_def);
  PyModule_AddStringConstant(version_info_module, "version", PERMANENT_VERSION);
  PyModule_AddStringConstant(version_info_module, "compiler_version", PERMANENT_COMPILER_VERSION);
  PyModule_AddStringConstant(version_info_module, "git_branch", PERMANENT_GIT_BRANCH);
  PyModule_AddStringConstant(version_info_module, "git_commit_hash", PERMANENT_GIT_COMMIT_HASH);

  /* Create main module. */
  auto *main_module = PyModule_Create(&main_def);
  PyModule_AddStringConstant(main_module, "__version__", PERMANENT_VERSION);
  PyModule_AddObject(main_module, "version_info", version_info_module);
  return main_module;
}

#undef PERMANENT_DISPATCH_ARRAY_TYPE
#undef PERMANENT_VERSION
#undef PERMANENT_COMPILER_VERSION
#undef PERMANENT_GIT_BRANCH
#undef PERMANENT_GIT_COMMIT_BRANCH
#undef PERMANENT_BUILD_TIME
#undef PERMANENT_COMPILER_VERSION
