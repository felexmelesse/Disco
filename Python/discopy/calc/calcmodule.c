#include <stdio.h>
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_11_API_VERSION
#include <numpy/arrayobject.h>
#include "../../../Calc/calc.h"


static char module_docstring[] = 
    "This module provides an interface for calculating various hydro solutions";
static char bondi_newt_docstring[] = 
    "Calculate the transonic Newtonian bondi solution";

static char bondi_rel_docstring[] = 
    "Calculate the Relativistic bondi solution";

static char magnetosonic_cf_int_newt_docstring[] = 
    "Calculate the integral for the magnetosonic Riemann invariant";

static PyObject *calc_bondi_newt(PyObject *self, PyObject *args);
static PyObject *calc_bondi_rel(PyObject *self, PyObject *args);
static PyObject *calc_magnetosonic_cf_int_newt(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"bondi_newt", calc_bondi_newt, METH_VARARGS, bondi_newt_docstring},
    {"bondi_rel", calc_bondi_rel, METH_VARARGS, bondi_rel_docstring},
    {"magnetosonic_cf_int_newt", calc_magnetosonic_cf_int_newt, METH_VARARGS,
        magnetosonic_cf_int_newt_docstring},
    {NULL, NULL, 0, NULL}};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef calcModule = {
    PyModuleDef_HEAD_INIT,
    "calc", //Module Name
    module_docstring,
    0,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

#define INITERROR return NULL
    
PyMODINIT_FUNC PyInit_calc(void)

#else

#define INITERROR return

PyMODINIT_FUNC initcalc(void)
#endif
{

#if PY_MAJOR_VERSION >= 3
    PyObject *m = PyModule_Create(&calcModule);
#else
    PyObject *m = Py_InitModule3("calc", module_methods, module_docstring);
#endif

    // return if there was a problem
    if(m == NULL)
        INITERROR;

    // Load numpy stuff!
    import_array();

#if PY_MAJOR_VERSION >= 3
    return m;
#endif
}

static PyObject *calc_bondi_newt(PyObject *self, PyObject *args)
{
    double Mdot, GM, gam, rho0;
    PyObject *r_obj = NULL;

    //Parse arguments
    if(!PyArg_ParseTuple(args, "ddddO", &Mdot, &GM, &gam, &rho0, &r_obj))
    {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't parse arguments.");
        return NULL;
    }

    //Grab numpy array
    PyArrayObject *r_arr;
    r_arr = (PyArrayObject *) PyArray_FROM_OTF(r_obj, 
                                        NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    // Throw exception
    if(r_arr == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't read r array.");
        Py_XDECREF(r_arr);
        return NULL;
    }

    //Check r is 1D
    int r_dim = (int)PyArray_NDIM(r_arr);
    if(r_dim != 1)
    {
        PyErr_SetString(PyExc_TypeError, "r must be 1-D");
        Py_DECREF(r_arr);
        return NULL;
    }

    int N = (int)PyArray_DIM(r_arr, 0);

    PyObject *rho_arr;
    PyObject *u_arr;
    PyObject *P_arr;
    npy_intp dims[1];
    dims[0] = N;
    rho_arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    u_arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    P_arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    
    // Throw exception
    if(rho_arr == NULL || u_arr == NULL || P_arr == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't make hydro arrays.");
        Py_DECREF(r_arr);
        Py_XDECREF(rho_arr);
        Py_XDECREF(u_arr);
        Py_XDECREF(P_arr);
        return NULL;
    }
    
    // Here's the actual array!
    double *r = (double *) PyArray_DATA(r_arr);
    double *rho = (double *) PyArray_DATA((PyArrayObject *) rho_arr);
    double *u = (double *) PyArray_DATA((PyArrayObject *) u_arr);
    double *P = (double *) PyArray_DATA((PyArrayObject *) P_arr);

    //Here's the function!
    int err = bondi_newt(Mdot, GM, gam, rho0, r, rho, u, P, N);

    //Clean!
    Py_DECREF(r_arr);

    //Build output
    PyObject *ret = Py_BuildValue("NNN", rho_arr, u_arr, P_arr);
    return ret;
}

static PyObject *calc_bondi_rel(PyObject *self, PyObject *args)
{
    double Mdot, GM, gam, a0;
    PyObject *r_obj = NULL;

    //Parse arguments
    if(!PyArg_ParseTuple(args, "ddddO", &Mdot, &GM, &gam, &a0, &r_obj))
    {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't parse arguments.");
        return NULL;
    }

    //Grab numpy array
    PyArrayObject *r_arr;
    r_arr = (PyArrayObject *) PyArray_FROM_OTF(r_obj, 
                                        NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    // Throw exception
    if(r_arr == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't read r array.");
        Py_XDECREF(r_arr);
        return NULL;
    }

    //Check r is 1D
    int r_dim = (int)PyArray_NDIM(r_arr);
    if(r_dim != 1)
    {
        PyErr_SetString(PyExc_TypeError, "r must be 1-D");
        Py_DECREF(r_arr);
        return NULL;
    }

    int N = (int)PyArray_DIM(r_arr, 0);

    PyObject *rho_arr;
    PyObject *u_arr;
    PyObject *P_arr;
    npy_intp dims[1];
    dims[0] = N;
    rho_arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    u_arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    P_arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    
    // Throw exception
    if(rho_arr == NULL || u_arr == NULL || P_arr == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't make hydro arrays.");
        Py_DECREF(r_arr);
        Py_XDECREF(rho_arr);
        Py_XDECREF(u_arr);
        Py_XDECREF(P_arr);
        return NULL;
    }
    
    // Here's the actual array!
    double *r = (double *) PyArray_DATA(r_arr);
    double *rho = (double *) PyArray_DATA((PyArrayObject *) rho_arr);
    double *u = (double *) PyArray_DATA((PyArrayObject *) u_arr);
    double *P = (double *) PyArray_DATA((PyArrayObject *) P_arr);

    //Here's the function!
    int err = bondi_rel(Mdot, GM, gam, a0, r, rho, u, P, N);

    //Clean!
    Py_DECREF(r_arr);

    //Build output
    PyObject *ret = Py_BuildValue("NNN", rho_arr, u_arr, P_arr);
    return ret;
}

static PyObject *calc_magnetosonic_cf_int_newt(PyObject *self, PyObject *args)
{
    double rho0, cs0, cA0, gam;
    PyObject *rho_obj = NULL;

    //Parse arguments
    if(!PyArg_ParseTuple(args, "Odddd", &rho_obj, &rho0, &cs0, &cA0, &gam))
    {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't parse arguments.");
        return NULL;
    }

    //Grab numpy array
    PyArrayObject *rho_arr;
    rho_arr = (PyArrayObject *) PyArray_FROM_OTF(rho_obj, 
                                        NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    // Throw exception
    if(rho_arr == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't read r array.");
        Py_XDECREF(rho_arr);
        return NULL;
    }

    //Check r is 1D
    int rho_dim = (int)PyArray_NDIM(rho_arr);
    if(rho_dim != 1)
    {
        PyErr_SetString(PyExc_TypeError, "rho must be 1-D");
        Py_DECREF(rho_arr);
        return NULL;
    }

    int N = (int)PyArray_DIM(rho_arr, 0);

    PyObject *res_arr;
    npy_intp dims[1];
    dims[0] = N;
    res_arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    
    // Throw exception
    if(res_arr == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't make array.");
        Py_DECREF(rho_arr);
        Py_XDECREF(res_arr);
        return NULL;
    }
    
    // Here's the actual array!
    double *rho = (double *) PyArray_DATA(rho_arr);
    double *res = (double *) PyArray_DATA((PyArrayObject *) res_arr);

    //Here's the function!
    int i;
    for(i=0; i<N; i++)
        res[i] = magnetosonic_cf_int_newt(rho[i], rho0, cs0, cA0, gam);

    //Clean!
    Py_DECREF(rho_arr);

    //Build output
    PyObject *ret = Py_BuildValue("N", res_arr);
    return ret;
}

