#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/ndarraytypes.h>

#include "permanent.h"


#define DOCSTRING_MODULE        "Permanent C extension module."
#define DOCSTRING_PERMANENT     "Compute the permanent of a matrix using the best algorithm."
#define DOCSTRING_COMBINATORIC  "Compute the permanent of a matrix combinatorically."
#define DOCSTRING_GLYNN         "Compute the permanent of a matrix via Glynn's algorithm."
#define DOCSTRING_RYSER         "Compute the permanent of a matrix via Ryser's algorithm."


static PyObject *py_combinatoric(PyObject *module, PyObject *object)
{
    PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, NULL, 2, 2, NPY_ARRAY_ALIGNED, NULL);
    return PyFloat_FromDouble(combinatoric(PyArray_DIMS(matrix)[0], PyArray_DIMS(matrix)[1], (double *)PyArray_GETPTR2(matrix, 0, 0)));
}


static PyObject *py_glynn(PyObject *module, PyObject *object)
{
    PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, NULL, 2, 2, NPY_ARRAY_ALIGNED, NULL);
    return PyFloat_FromDouble(glynn(PyArray_DIMS(matrix)[0], PyArray_DIMS(matrix)[1], (double *)PyArray_GETPTR2(matrix, 0, 0)));
}


static PyObject *py_ryser(PyObject *module, PyObject *object)
{
    PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, NULL, 2, 2, NPY_ARRAY_ALIGNED, NULL);
    return PyFloat_FromDouble(ryser(PyArray_DIMS(matrix)[0], PyArray_DIMS(matrix)[1], (double *)PyArray_GETPTR2(matrix, 0, 0)));
}


/* Define the Python methods that will go into the C extension module. */
static PyMethodDef methods[] = {
    /* Python function name     C function          Args flag   Docstring */
    { "combinatoric",           py_combinatoric,    METH_O,     DOCSTRING_COMBINATORIC },
    { "glynn",                  py_glynn,           METH_O,     DOCSTRING_GLYNN },
    { "ryser",                  py_ryser,           METH_O,     DOCSTRING_RYSER },
    { NULL,                     NULL,               0,          NULL } /* sentinel value */
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
