"""
Contains integrator functions for bodies orbiting around central object such as a star.

To use write:
from interfunc import *
"""

from numba import cfunc, types
import numpy as np

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
    params[1,...] = the standard gravitation parameters for each other object
    
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

    #Jupiter
    GM_J = params[1]
    r_J = params[2]
    period_J = params[3]
    theta_J = 2*np.pi*t/period_J
    x_J = r_J*np.cos(theta_J)
    y_J = r_J*np.sin(theta_J)
    x_diff = x - x_J
    y_diff = y - y_J
    r_diff = np.sqrt(x_diff**2 + y_diff**2)
    
    #Particle Output
    dydt[0] = v_x #dx/dt
    dydt[1] = -GM_S/(r**3)*x - GM_J/(r_diff**3)*x_diff #dv_x/dt
    dydt[2] = v_y #dy/dt
    dydt[3] = -GM_S/(r**3)*y - GM_J/(r_diff**3)*y_diff #dv_y/dt

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