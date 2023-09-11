#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyObject* hamming(PyObject* self, PyObject* args) {
    const char *s1;
    const char *s2;
    Py_ssize_t len;
    int distance = 0;

    if (!PyArg_ParseTuple(args, "s#s#", &s1, &len, &s2, &len)) {
        return NULL;
    }

    const char *end = s1 + len;
    while (s1 < end) {
        if (*s1 != *s2) {
            distance++;
        }
        s1++;
        s2++;
    }

    return PyLong_FromLong(distance);
}

static PyMethodDef methods[] = {
    {"hamming", hamming, METH_VARARGS, "Calculate the Hamming distance between two strings"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "hamming",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC PyInit_hamming(void) {
    return PyModule_Create(&module);
}
