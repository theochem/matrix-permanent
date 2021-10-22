/* Preamble */
/* ******** */


/* Python C header. */
#include <Python.h>


/* Include whatever NumPy C headers you'll need here. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/ndarraytypes.h"
#include "numpy/npy_3kcompat.h"
#include "numpy/ufuncobject.h"


/* Python documentation strings. */
#define DOCSTRING_MODULE        "Permanent C extension module."
#define DOCSTRING_PERMANENT     "Compute the permanent of a matrix using the best algorithm."
#define DOCSTRING_COMBINATORIC  "Compute the permanent of a matrix combinatorically."
#define DOCSTRING_GLYNN         "Compute the permanent of a matrix via Glynn's algorithm."
#define DOCSTRING_RYSER         "Compute the permanent of a matrix via Ryser's algorithm."


/* C functions for the Python module */
/* ********************************* */


static PyObject* combinatoric(PyObject *module, PyObject *matrix)
{
    return PyInt_FromLong(0L);
}


static PyObject* glynn(PyObject *module, PyObject *matrix)
{
    return PyInt_FromLong(0L);
}


static PyObject* ryser(PyObject *module, PyObject *matrix)
{
    return PyInt_FromLong(0L);
}


static PyObject* permanent(PyObject *module, PyObject *matrix)
{
    return PyInt_FromLong(0L);
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
