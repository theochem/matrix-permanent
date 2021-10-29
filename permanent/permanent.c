/* Preamble */
/* ******** */


/* Python C header. */
#include <Python.h>
#include <stdbool.h>


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


void swap2(uint64_t *perm, uint64_t i, uint64_t j)
// A function for swapping the values of two entries in an array while keeping the pointer the same
{
    const uint64_t temp = perm[i];
    perm[i] = perm[j];
    perm[j] = temp;
}

void init_perm(const uint64_t N, const uint64_t M, uint64_t *the_fact_set, uint64_t *the_perm_set, uint64_t *the_inv_set)
{
    // generate the set {1, 2, .., N} to permute
    // update all relevent sets we must keep track up in order to siplify generation of next permutation
    // update values pointed to by respective array pointers
    // first set is in ascending order to ensure no smaller permutation is possible

    for (uint64_t i = 0; i < N; i++)
    {
        the_fact_set[i + 1] = 0;
        the_perm_set[i] = i;
        the_inv_set[i] = i;
    }
}


bool gen_next_perm(uint64_t *const falling_fact, uint64_t *const perm_, uint64_t *const inv_perm, const uint64_t rows, const uint64_t cols, const uint64_t u_)
{
    /* use the current permutation to check for next permutation lexicographically, if it exists
    update the curr_array by swapping the leftmost changed element (at position i < k).
    Replace the elements up to position k by the smallest element lying to the right of i.
    Do not put remaining k, .., n - 1 positions in ascending order to improve efficiency
    has time complexity O(n) in the worst case */
    int64_t i = rows - 1;
    int64_t m1 = cols - i - 1;
    while (falling_fact[i + 1] == m1)
    // begin update of falling_fact and generate next permutation
    {
        falling_fact[i + 1] = 0;
        ++m1;
        --i;
    }

    if (i == -1)
    {
        return false;
        // we have already generated all permutations
    }
    ++falling_fact[i + 1];
    // find smallest element p[i] < p[j] that lies to the right of pos i
    uint64_t z = perm_[i];
    // update the permutation
    do 
    {
        ++z;
    } while (inv_perm[z] <= i);
    const uint64_t j = inv_perm[z];
    swap2(perm_, i, j);
    swap2(inv_perm, perm_[i], perm_[j]);
    ++i;
    uint64_t y = 0;
    while (i < u_)
    {
        // find smallest elements to the right of position i
        while (inv_perm[y] < i)
        {
            ++y;
        }
        const uint64_t j = inv_perm[y];
        swap2(perm_, i, j);
        swap2(inv_perm, perm_[i], perm_[j]);
        ++i;
    }
    return true;
}


static PyObject *combinatoric(PyObject *module, PyObject *object)
{
    /* Cast PyObject* to PyArrayObject*. */
    PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, NULL, 2, 2, NPY_ARRAY_ALIGNED, NULL);

    /* Return element (0, 0) of the matrix. */
    double *ptr = (double *)PyArray_GETPTR2(matrix, 0, 0);
    /* return PyFloat_FromDouble(*ptr); */

    /* Return the permanent of the matrix. */
    uint64_t m_rows = PyArray_DIMS(matrix)[0];
    uint64_t n_cols = PyArray_DIMS(matrix)[1];
    uint64_t sort_up_to;
    // sort up to position u + 1 where u = min(k, n_cols - 1)

    if (m_rows < n_cols)
    {
        sort_up_to = m_rows;
    }
    else
    {
        sort_up_to = n_cols - 1;
    }

    double sum_permanent = 0.0;
    double prod_permanent = 1.0;


    // generate a pointer to the set of elements {1, 2, .., N} and determing the position of the leftmost change
    // store the pointer to the array in variable falling_fact
    uint64_t falling_fact[128];
    falling_fact[0] = 0;
    // this value will be used to check if we have generated all permutations - allows i to become -1 in th generator
    uint64_t perm_[128];
    uint64_t inv_perm[128];
    // inverse permutation is used to simplify the update of the permanent

    init_perm(n_cols, m_rows, falling_fact, perm_, inv_perm);
    // while (falling_fact[1] < (n_cols - m_rows + 1))
    for (uint64_t i = 0; i < 4; i++)
            printf("perm_[%" PRIu64 "] = %" PRIu64 "\n", i, perm_[i]);

    bool gen_next_perm();
    // get value from initialization before trying to generate next perm.
    for (uint64_t i = 0; i < m_rows; ++i)
        {
            prod_permanent *= (ptr[i * n_cols + perm_[i]]);
        }
    sum_permanent = prod_permanent;
    printf("%f", sum_permanent);
    while (gen_next_perm(falling_fact, perm_, inv_perm, m_rows, n_cols, sort_up_to))
    {
        printf("We made it into the loop!");
        prod_permanent = 1.0;
        for (uint64_t i = 0; i < m_rows; ++i)
        {
            prod_permanent *= (ptr[i * n_cols + perm_[i]]);
        }

        sum_permanent += prod_permanent;
        printf("%f", sum_permanent);        

        for (uint64_t i = 0; i < 3; i++)
            printf("perm_[%" PRIu64 "] = %" PRIu64 "\n", i, perm_[i]);

        for (uint64_t i = 0; i < 3; i++)
            printf("falling_fact[%" PRIu64 "] = %" PRIu64 "\n", i, falling_fact[i + 1]);
    }
    return PyFloat_FromDouble(sum_permanent);
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
