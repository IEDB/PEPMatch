#define PY_SSIZE_T_CLEAN
#include <Python.h>


static PyObject* hamming(PyObject* self, PyObject* args) {
    const char *kmer1;
    const char *kmer2;
    int max_mismatches;
    Py_ssize_t len1, len2;

    if (!PyArg_ParseTuple(args, "s#s#i", &kmer1, &len1, &kmer2, &len2, &max_mismatches)) {
        return NULL;
    }

    int mismatches = 0;
    for (int i = 0; kmer1[i] && kmer2[i]; ++i) {
        if (kmer1[i] != kmer2[i]) {
            mismatches++;
            if (mismatches > max_mismatches) {
                return PyLong_FromLong(max_mismatches + 1);
            }
        }
    }

    return PyLong_FromLong(mismatches);
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
