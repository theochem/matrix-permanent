/* The permanent commonly appears in problems related to quantum mechanics, and the most common
brute-force combinatorical method has time complexity O(N!N), thus it is useful to look
for more efficient algorithms. The two algorithms considered to be the fastest are one by
Ryser (based on the inclusion-exclusion principle), and one by Glynn (based on invariant theroy). 
All algorithms work for square NxN matrices, and are generalizable for MxN matrices.

The goal is to optimize the code and find the best algorithm for each value of M and N, and 
have a C++ function that will automatically find the best algorithm based on the size of the 
input matrix.

This code works with an input matrix in the form of Python NumPy array.

The function for generating k-permtations lexicographically was adapted from
"Matters Computational; Ideas, Algorithms, Source Code" by Jorg Arndt" Ch 12.1*/
/* ******** */


/* Python C header. */
#include <Python.h>


/* C headers. */
#include <stdbool.h>
#include <stdio.h>
#include <tgmath.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>


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

/* A function to determine the greatest common denominator.*/
int64_t gcd(int64_t x, int64_t y) 
{
    int64_t t;
    if (y < x) 
    {
        t = x;
        x = y;
        y = t;
    }
    while (y > 0) 
    {
        t = y;
        y = x % y;
        x = t;
    }
    return x;
}

/* A function to determine the binomial coefficient. Adapted from Micheal Richer's code as reference. */
int64_t bin_coeff(int64_t n, int64_t k)
{
    if (k == 0)
    {
        return 1;
    }
    else if (k == 1)
    {
        return n;
    }
    else if (k >= n)
    {
        return (int64_t) k == n;
    }
    else if (k > n / 2)
    {
        k = n - k;
    }
    int64_t binom = 1, g, nr, dr;
    int64_t gcd();
    for (int64_t d = 1; d <= k; ++d) 
    {
        if (binom >= LONG_MAX / n) {
            g = gcd(n, d);
            nr = n / g;
            dr = d / g;
            g = gcd(binom, dr);
            binom = binom / g;
            dr = dr / g;
            if (binom >= LONG_MAX / nr)
            {
                return LONG_MAX;
            }
            binom = binom * nr / dr;
        } 
        else
        {
            binom = binom * n-- / d;
        }
    }
    return binom;
}

/* A function for swapping the values of two entries in an array while maintaining value of the pointer. */
void swap2(int64_t *perm, int64_t i, int64_t j)
{
    const int64_t temp = perm[i];
    perm[i] = perm[j];
    perm[j] = temp;
}

/* A function that will initialize the set to permute. This first set is used as the first permutation. The
set is generated in ascedning order to ensure there are no smaller permutations possible. We keep track of the 
inverse permutation in order to simplify the swap update when generating the next permutation. The ffactorial
set is initialized to all zeroes and is used to know when we have generated all possible permutations.*/
void init_perm(const int64_t N, const int64_t M, int64_t *the_fact_set, int64_t *the_perm_set, int64_t *the_inv_set)
{
    for (int64_t i = 0; i < N; i++)
    {
        the_fact_set[i + 1] = 0;
        the_perm_set[i] = i;
        the_inv_set[i] = i;
    }
}

/* A function that will use the current state of the permutation perm_ and update it to reflect that of the next permutation.
If the largest permutation has been generated, the function will return false, else it will determine the next permutation and 
update all three arrays (falling_fact, perm_, inv_perm). This was adapted from Jorgs "Matters Computational" C++ reference.*/
bool gen_next_perm(int64_t *const falling_fact, int64_t *const perm_, int64_t *const inv_perm, const int64_t rows, const int64_t cols, const int64_t u_)
{
    /* Use the current permutation to check for next permutation lexicographically, if it exists
    update the curr_array by swapping the leftmost changed element (at position i < k).
    Replace the elements up to position k by the smallest element lying to the right of i.
    Do not put remaining k, .., n - 1 positions in ascending order to improve efficiency
    has time complexity O(n) in the worst case. */
    int64_t i = rows - 1;
    int64_t m1 = cols - i - 1;
    /* begin update of falling_fact - recall falling_fact[0] = 0, so increment index accordingly.
    If i becomes negative during the check, you are pointing to the sentinal value, so you are
    done generating permutations. */
    while (falling_fact[i + 1] == m1)
    {
        falling_fact[i + 1] = 0;
        ++m1;
        --i;
    }
    if (i == -1)
    {
        return false;
    }
    ++falling_fact[i + 1];
    /* Find smallest element perm_[i] < perm_[j] that lies to the right of pos i,
    and then update the state of the permutation using its inverse to generate next. */
    int64_t z = perm_[i];
    do 
    {
        ++z;
    } while (inv_perm[z] <= i);
    const int64_t j = inv_perm[z];
    swap2(perm_, i, j);
    swap2(inv_perm, perm_[i], perm_[j]);
    ++i;
    int64_t y = 0;
    /* Find the smallest elements to the right of position i. */
    while (i < u_)
    {
        while (inv_perm[y] < i)
        {
            ++y;
        }
        const int64_t j = inv_perm[y];
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
    double *ptr = (double *)PyArray_GETPTR2(matrix, 0, 0);

    /* Return the permanent of the matrix. */
    int64_t m_rows = PyArray_DIMS(matrix)[0];
    int64_t n_cols = PyArray_DIMS(matrix)[1];
    int64_t sort_up_to;
    /* sort up to position u + 1 where u = min(k, n_cols - 1). */
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
    /* Allocate falling_fact, perm_, and inv_perm arrays. */
    int64_t falling_fact[128];
    /* Set sentinal value. */
    falling_fact[0] = 0;
    int64_t perm_[128];
    int64_t inv_perm[128];
    init_perm(n_cols, m_rows, falling_fact, perm_, inv_perm);
    bool gen_next_perm();
    /* Handle first permutation. */
    for (int64_t i = 0; i < m_rows; ++i)
        {
            prod_permanent *= (ptr[i * n_cols + perm_[i]]);
        }
    sum_permanent = prod_permanent;
    /* Iterate over second to last permutations. */
    while (gen_next_perm(falling_fact, perm_, inv_perm, m_rows, n_cols, sort_up_to))
    {
        prod_permanent = 1.0;
        for (int64_t i = 0; i < m_rows; ++i)
        {
            prod_permanent *= (ptr[i * n_cols + perm_[i]]);
        }
        sum_permanent += prod_permanent;
    }
    return PyFloat_FromDouble(sum_permanent);
}

/* This code for generating the permanent of a square NxN matrix was adapted from Michael Richer's draft code. */
static PyObject *glynn(PyObject *module, PyObject *object)
{
    /* Cast PyObject* to PyArrayObject*. */
    PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, NULL, 2, 2, NPY_ARRAY_ALIGNED, NULL);

    /* Return the permanent of a square NxN matrix. */
    double *ptr = (double *)PyArray_GETPTR2(matrix, 0, 0);
    int64_t m_rows = PyArray_DIMS(matrix)[0];
    int64_t n_cols = PyArray_DIMS(matrix)[1];

    /* Initialize gray code. */
    int64_t pos = 0;
    int64_t sign = 1;
    int64_t bound = n_cols - 1;
    int64_t delta[n_cols];
    /* Allocate and fill delta array (all +1 to start). */
    for (int64_t i = 0; i < n_cols; ++i)
    {
        delta[i] = 1;
    }
    int64_t gray[n_cols];
    /* Allocate and fill gray code from 0 to n_cols. */
    for (int64_t i = 0; i < n_cols; ++i)
    {
        gray[i] = i;
    }
    /* Iterate over every delta. */
    /* Allocate matrix for result of manual multiplies. */
    double vec[n_cols];
    double *vec_end = vec + n_cols;
    /* Handle first Gray code. */
    double result = 1.0, sum_ = 0.0, *vec_it;
    for (int64_t i = 0; i < n_cols; ++i)
    {
        for (int64_t j = 0; j < m_rows; ++j)
        {
            sum_ += (ptr[i * n_cols + j] * delta[i]);
        }
        vec[i] = sum_;
    }
    for (vec_it = vec; vec_it != vec_end; ++vec_it)
    {
        result *= *vec_it;
    }
    /* Iterate over second to last Gray codes. */
    double prod;
    while (pos != bound) 
    {
        /* Update sign and delta. */
        sign *= -1;
        *(delta + bound - pos) *= -1;
        /* Compute each Gray code term. */
        sum_ = 0;
        for (int64_t i = 0; i < n_cols; ++i)
        {
            for (int64_t j = 0; j < m_rows; ++j)
            {
                sum_ += (ptr[i * n_cols + j] * delta[i]);
            }
            vec[i] = sum_;
        }
        prod = 1.0;
        for (vec_it = vec; vec_it != vec_end; ++vec_it)
        {
            prod *= *vec_it;
        }
        result += sign * prod;
        /* Go to next Gray code. */
        *gray = 0;
        *(gray + pos) = *(gray + pos + 1);
        ++pos;
        *(gray + pos) = pos;
        pos = *gray;
    }
    /* Divide by external factor and return permanent. */
    return PyFloat_FromDouble(result / (double) pow(2, bound));
}


static PyObject *ryser(PyObject *module, PyObject *object)
{
    /* Cast PyObject* to PyArrayObject*. */
    PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, NULL, 2, 2, NPY_ARRAY_ALIGNED, NULL);

    /* Return the permanent of a rectangular matrix, where M < N. */
    double *ptr = (double *)PyArray_GETPTR2(matrix, 0, 0);
    int64_t m_rows = PyArray_DIMS(matrix)[0];
    int64_t n_cols = PyArray_DIMS(matrix)[1];

    /* Notice we generate permutations of sets the same way as in combinatoric algorithm.
    See above algorithm for comments on declared arrays and permutation generator. */
    int64_t sort_up_to;
    double sum_of_matrix_vals = 0.0;
    double sum_over_k_vals = 0.0;
    double sum_over_perms = 0.0;
    double prod_of_cols = 1.0;
    int64_t falling_fact[128];
    falling_fact[0] = 0;
    int64_t perm_[128];
    int64_t inv_perm[128];

    for (int64_t k = 0; k < m_rows - 1; ++k)
    {
        init_perm(n_cols, m_rows - k, falling_fact, perm_, inv_perm);
        bool gen_next_perm();
        sort_up_to = m_rows - k;

        for (int64_t i = 0; i < m_rows; ++i)
        {
            for (int64_t j = 0; j < m_rows - k; ++j)
            {
                sum_of_matrix_vals += (ptr[i * n_cols + perm_[j]]);
            }
            prod_of_cols *= sum_of_matrix_vals;
        }
        int64_t value_sign = pow(-1, k);
        int64_t bin_coeff();
        int64_t bin_c = bin_coeff(n_cols - m_rows + k, k); 
        sum_over_perms += value_sign * bin_c * prod_of_cols;
        while (gen_next_perm(falling_fact, perm_, inv_perm, m_rows - k, n_cols, sort_up_to))
        {
            sum_of_matrix_vals = 0.0;
            for (int64_t i = 0; i < m_rows; ++i)
            {
                for (int64_t j = 0; j < m_rows - k; ++j)
                {
                    sum_of_matrix_vals += (ptr[i * n_cols + perm_[j]]);
                }
                prod_of_cols *= sum_of_matrix_vals;
            }
            int64_t value_sign = pow(-1, k);
            int64_t bin_coeff();
            int64_t bin_c = bin_coeff(n_cols - m_rows + k, k); 
            sum_over_perms += value_sign * bin_c * prod_of_cols;
        }
        sum_over_k_vals += sum_over_perms;
    }
    return PyFloat_FromDouble(sum_over_k_vals);
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
