"""
Contains integrator functions for bodies orbiting around central object such as a star.

To use write:
from interfunc import *
"""

from numba import cfunc, types
import numpy as np

n_var_per_object = 4 #x, v_x, y, v_y

c_sig_index_funcs = types.int32(types.int32)

@cfunc(c_sig_index_funcs)
def ind_x(i):
    """ Gives index for x in sol/dydt vectors corresponding to object i """
    return n_var_per_object * i + 0

@cfunc(c_sig_index_funcs)
def ind_v_x(i):
    """ Gives index for v_x in sol/dydt vectors corresponding to object i """
    return n_var_per_object * i + 1

@cfunc(c_sig_index_funcs)
def ind_y(i):
    """ Gives index for y in sol/dydt vectors corresponding to object i """
    return n_var_per_object * i + 2

@cfunc(c_sig_index_funcs)
def ind_v_y(i):
    """ Gives index for v_y in sol/dydt vectors corresponding to object i """
    return n_var_per_object * i + 3

c_sig = types.void(types.double,
                   types.CPointer(types.double),
                   types.CPointer(types.double),
                   types.CPointer(types.double))

@cfunc(c_sig)
def func_n_body(t,sol,params,dydt):
    """
    For n objects orbiting around a central object (e.g. a star)
    sol = [x, v_x, y, v_y]*(number of objects)
    params[0] = GM_S, the standard gravitational parameter of the star/object being orbited
    params[1] = n, the number of objects orbiting the central object
    params[2,...] = the standard gravitation parameters for each other object
    
    Takes the time and solution vector sol, as well as parameters in the params array,
    and populates the dydt array with the values of the right-hand side function 
    f(t,y).
    Returns: nothing
    """

    GM_S = params[0]
    n = params[1]
    
    for i in range(n):
        #Object input values
        x = sol[ind_x(i)]
        v_x = sol[ind_v_x(i)]
        y = sol[ind_y(i)]
        v_y = sol[ind_v_y(i)]
        r = np.sqrt(x**2 + y**2)

        #Initial object output values
        dydt[ind_x(i)] = v_x #dx/dt
        dydt[ind_v_x(i)] = -GM_S/(r**3)*x #dv_x/dt
        dydt[ind_y(i)] = v_y #dy/dt
        dydt[ind_v_y(i)] = -GM_S/(r**3)*y #dv_y/dt

        #Adding on forces from other orbiting objects
        for j in range(n):
            if (j == i): continue
            GM_j = params[j+2] # Skipping GM_S and n in params
            x_j = sol[ind_x(j)]
            y_j = sol[ind_y(j)]
            x_diff = x - x_j
            y_diff = y - y_j
            r_diff = np.sqrt(x_diff**2 + y_diff**2)
            dydt[ind_v_x(i)] = dydt[ind_v_x(i)] - GM_j/(r_diff**3)*x_diff #dv_x/dt
            dydt[ind_v_y(i)] = dydt[ind_v_y(i)] - GM_j/(r_diff**3)*y_diff #dv_y/dt

@cfunc(c_sig)
def func_2_body(t,sol,params,dydt):
    """
    For 1 object orbiting around another (e.g. a star)
    sol = x, v_x, y, v_y
    params[0] = GM_S, the standard gravitational parameter of the star/object being orbited
    
    Takes the time and solution vector sol, as well as parameters in the params array,
    and populates the dydt array with the values of the right-hand side function 
    f(t,y).
    Returns: nothing
    """

    #Particle Input
    x = sol[0]
    v_x = sol[1]
    y = sol[2]
    v_y = sol[3]
    GM_S = params[0]
    r = np.sqrt(x**2 + y**2)
    
    #Particle Output
    dydt[0] = v_x #dx/dt
    dydt[1] = -GM_S/(r**3)*x #dv_x/dt
    dydt[2] = v_y #dy/dt
    dydt[3] = -GM_S/(r**3)*y #dv_y/dt