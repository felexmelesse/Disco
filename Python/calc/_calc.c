#include <stdio.h>
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_11_API_VERSION
#include <numpy/arrayobject.h>
#include "../../Calc/calc.h"


static char module_docstring[] = 
    "This module provides an interface for calculating various hydro solutions";
static char bondi_newt_docstring[] = 
    "Calculate the transonic Newtonian bondi solution";

static char bondi_rel_docstring[] = 
    "Calculate the Relativistic bondi solution";

static PyObject *calc_bondi_newt(PyObject *self, PyObject *args);
static PyObject *calc_bondi_rel(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"bondi_newt", calc_bondi_newt, METH_VARARGS, bondi_newt_docstring},
    {"bondi_rel", calc_bondi_rel, METH_VARARGS, bondi_rel_docstring},
    {NULL, NULL, 0, NULL}};

PyMODINIT_FUNC init_calc(void)
{
    PyObject *m = Py_InitModule3("_calc", module_methods, module_docstring);

    // return if there was a problem
    if(m == NULL)
        return;

    // Load numpy stuff!
    import_array();
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

