/* Preamble */
/* ******** */


/* Python C header. */
#include <Python.h>


/* NumPy C headers. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/ndarraytypes.h>


/* Python documentation strings. */
#define DOCSTRING_MODULE        "Permanent C extension module."
#define DOCSTRING_PERMANENT     "Compute the permanent of a matrix using the best algorithm."
#define DOCSTRING_COMBINATORIC  "Compute the permanent of a matrix combinatorically."
#define DOCSTRING_GLYNN         "Compute the permanent of a matrix via Glynn's algorithm."
#define DOCSTRING_RYSER         "Compute the permanent of a matrix via Ryser's algorithm."


/* C functions for the Python module */
/* ********************************* */


void swap2(uint64_t *const curr_perm, uint64_t i, uint64_t j)
// A function for swapping the values of two entries in an array while keeping the pointer the same
{
    const uint64_t temp = curr_perm[i];
    curr_perm[i] = curr_perm[j];
    curr_perm[j] = temp;
}

void re_sort(uint64_t *const curr_perm, uint64_t start)
// Re-sort the set to be in ascending order - this occurs after you have swapped adjacent values
{
    uint64_t i = start;
    uint64_t j = 0;
    while (i > j)
    {
        swap2(curr_perm, i, j);
        i--;
        j++;
    }
}


static PyObject *combinatoric(PyObject *module, PyObject *object)
{
    /* Cast PyObject* to PyArrayObject*. */
    PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, NULL, 2, 2, NPY_ARRAY_ALIGNED, NULL);

    /* Return element (0, 0) of the matrix. */
    double *ptr = (double *)PyArray_GETPTR2(matrix, 0, 0);
    return PyFloat_FromDouble(*ptr);

    /* Return the permanent of the matrix. */ 
    uint64_t m_rows = PyArray_DIMS(matrix)[0];
    uint64_t n_cols = PyArray_DIMS(matrix)[1];
    double sum_permanent = 0.0;
    double prod_permanent = 1.0;


    uint64_t init_perm(const uint64_t N)
    {
        // generate the first permutation of the set {1, 2, .., N}
        // return a pointer to the array
        // first set is in ascending order

        uint64_t the_set[N];
        for (uint64_t i = 0; i < N; i++)
        {
            the_set[i] = i + 1
        }
        return the_set
    }

    void next_perm(uint64_t *const curr_perm, const uint64_t length)
    {
        // use the current permutation to generate the next permutation lexicographically
        // save memory by changing contents in init_perm while keeping pointer the same
        // has time complexity O(n) in the worst case
        uint64_t i = 1;
        while (i < length - 2 && curr_perm[i] < curr_perm[i + 1])
        {
            i++;
        }
        if (i = length - 1)
        {
            uint64_t j = 0;
            while (curr_perm[j] > curr_perm[i])
            {
                j--;
            }
            swap2(curr_perm, i, j);
        }
        re_sort(curr_perm, i - 1);
    }

    // generate a pointer to the set of elements {1, 2, .., N} in ascending order
    // store the pointer to the array in variable curr_perm
    curr_perm = init_perm(n_cols);
    while (curr_perm[0] < (n_cols - m_rows + 1))
    {
        for (uint64_t i = 0; i < m_rows; ++i)
        {
            prod_permanent *= (matrix[i, curr_perm[i]]);
        }

        sum_permanent += prod_permanent;
        // generate next permutation lexicographically and update the array pointed to by curr_perm
        next_perm(curr_perm, n_cols);
    }
    return sum_permanent;
}


static PyObject *glynn(PyObject *module, PyObject *object)
{
    /* Cast PyObject* to PyArrayObject*. */
    PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, NULL, 2, 2, NPY_ARRAY_ALIGNED, NULL);

    /* Return element (0, 0) of the matrix. */
    double *ptr = (double *)PyArray_GETPTR2(matrix, 0, 0);
    return PyFloat_FromDouble(*ptr);
}


static PyObject *ryser(PyObject *module, PyObject *object)
{
    /* Cast PyObject* to PyArrayObject*. */
    PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, NULL, 2, 2, NPY_ARRAY_ALIGNED, NULL);

    /* Return element (0, 0) of the matrix. */
    double *ptr = (double *)PyArray_GETPTR2(matrix, 0, 0);
    return PyFloat_FromDouble(*ptr);
}


static PyObject *permanent(PyObject *module, PyObject *object)
{
    /* Cast PyObject* to PyArrayObject*. */
    PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, NULL, 2, 2, NPY_ARRAY_ALIGNED, NULL);

    /* Return element (0, 0) of the matrix. */
    double *ptr = (double *)PyArray_GETPTR2(matrix, 0, 0);
    return PyFloat_FromDouble(*ptr);
}


/* Defining the Python module */
/* ************************** */


/* Define the Python methods that will go into the C extension module. */
static PyMethodDef methods[] = {
    /* Python function name     C function      Args flag   Docstring */
    { "permanent",              permanent,      METH_O,     DOCSTRING_PERMANENT },
    { "combinatoric",           combinatoric,   METH_O,     DOCSTRING_COMBINATORIC },
    { "glynn",                  glynn,          METH_O,     DOCSTRING_GLYNN },
    { "ryser",                  ryser,          METH_O,     DOCSTRING_RYSER },
    { NULL,                     NULL,           0,          NULL } /* sentinel value */
};
/* Note:  METH_O indicates that the Python function takes a single argument.
 *        On the C side, the function takes two PyObject* arguments;
 *        the first one is the C extension module itself,
 *        and the second one is the argument to the Python function. */


/* Define the C extension module. */
static struct PyModuleDef definition = {
    /* Shouldn't need to change this. */
    PyModuleDef_HEAD_INIT, "permanent", DOCSTRING_MODULE, -1, methods
};


/* Initialize the C extension module. */
PyMODINIT_FUNC PyInit_permanent(void) {
    Py_Initialize();    /* Initialize Python API */
    import_array();     /* Initialize NumPy NDArray API */
    return PyModule_Create(&definition); /* Create module. */
}
