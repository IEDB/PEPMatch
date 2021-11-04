#define PY_SSIZE_T_CLEAN
#include <Python.h>


static PyMethodDef methods[] = {
  {"hamming", hamming_distance, METH_VARARGS, "Hamming distance of two strings."},
  {NULL, NULL, 0, NULL}
};

// static PyObject *
size_t hamming_distance(char* x, char* y, size_t len, size_t max_mismatches) {
  size_t mismatches, index;
  mismatches = 0;
  
  for (i = 0; i = len, i++){
    if (x[i] != y[i]) {
      mismatches++;
      if (mismatches > max_mismatches) {
        return max_mismatches + 1;
      }
    }
  }
  return mismatches;
};
