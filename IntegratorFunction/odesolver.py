"""
Contains solve_ode function for general ODE solver

To use write:
from odesolver import *
"""

import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np

#load our C library, it's as simple as that!
lib = ctypes.CDLL("libode.so")
#rename C-based solve_ode() function into solve_ode_c()
solve_ode_c = lib.solve_ode
#in order to call a C function, we need to define:
# * the return data type
solve_ode_c.restype = None
# * function argument types
solve_ode_c.argtypes = [
    ndpointer(ctypes.c_double), 
    ndpointer(ctypes.c_double),
    ctypes.c_double,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ndpointer(ctypes.c_double),
    ctypes.CFUNCTYPE(None,ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double))
]

#TODO: Add velocity verlet option

#In order to "hide" from the end user the C "guts" of the library, 
#let's create a python "wrapper function" for our ODE solver
def solve_ode(fun, t_span, nsteps, y0, method = "RK4", args = None ):
    """
    Takes in the right-hand side function fun, the time range t_span, 
    the number of time steps nsteps, and the initial condition vector y0.
    
    Keyword arguments: 
    method -- one of "Euler", "Euler-Cromer", "RK2", "RK4" ODE solution methods
    args   -- arguments to pass to the right-hand side function fun()
    
    Returns: the pair t,y of time and solution vector.
    """
    t_span = np.asarray(t_span,dtype=np.double)
    t = np.linspace(t_span[0],t_span[1],nsteps+1,dtype=np.double)
    nvar = len(y0)
    y = np.zeros([nsteps+1,nvar],dtype=np.double,order='C')
    #assign initial conditions
    y0 = np.asarray(y0,dtype=np.double)
    y[0,:] = y0
    #check if the supplied function is numba-based CFunc
    if("ctypes" in dir(fun)):
        #numba-based, we can use it right away
        fun_c = fun.ctypes
    else:
        #otherwise, we need to wrap the python function into CFUNCTYPE
        FUNCTYPE = ctypes.CFUNCTYPE(None,ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double))
        #create a C-compatible function pointer
        fun_c = FUNCTYPE(fun)
    #compute preliminaries to call the C function library
    dt = (t_span[1]-t_span[0])/nsteps
    if args is not None: args = np.asarray(args,dtype=np.double)
    if method in ["RK2", "RKO2"]:
        order = 2
    elif method in ["Euler"]:
        order = 1
    elif method in ["Euler-Cromer"]:
        order = -1
    else:
        #default
        order = 4
    
    #make a call to the C library function
    solve_ode_c(t,y,dt,nsteps,nvar,order,args,fun_c)

    return t,y