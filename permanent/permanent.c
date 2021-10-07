#include <Python.h>

/* Include whatever NumPy headers you'll need here. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/ndarraytypes.h"
#include "numpy/npy_3kcompat.h"
#include "numpy/ufuncobject.h"


static PyObject* combinatoric(PyObject *dummy, PyObject *args)
{
    return PyInt_FromLong(0L);
}


static PyObject* glynn(PyObject *dummy, PyObject *args)
{
    return PyInt_FromLong(0L);
}


static PyObject* ryser(PyObject *dummy, PyObject *args)
{
    return PyInt_FromLong(0L);
}


static PyObject* permanent(PyObject *dummy, PyObject *args)
{
    return PyInt_FromLong(0L);
}


static PyMethodDef methods[] = {
    { "combinatoric", combinatoric, METH_VARARGS, "numpy function tester" },
    { "glynn", glynn, METH_VARARGS, "numpy function tester" },
    { "ryser", ryser, METH_VARARGS, "numpy function tester" },
    { "permanent", permanent, METH_VARARGS, "numpy function tester" },
    { NULL, NULL, 0, NULL }
};


static struct PyModuleDef definition = {
    PyModuleDef_HEAD_INIT, "permanent", "Permanent C extension module.", -1, methods
};


PyMODINIT_FUNC PyInit_permanent(void) {
    Py_Initialize();
    import_array();
    return PyModule_Create(&definition);
}
