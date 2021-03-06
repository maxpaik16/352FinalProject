import numpy as np
import math

# Gives orbital period for orbit of semi-major axis "a" around a star (or any body being orbited)
def orbital_period(a, GM_star):
    return 2*np.pi*np.sqrt(a**3/GM_star)

# Removes horizontal jumps in plots of x and y
# For example for chaotic phase space plots
def remove_horizontal_discontinuities(x,y):
    pos = np.where(np.abs(np.diff(x)) >= 6)[0]+1
    x = np.insert(x, pos, np.nan)
    y = np.insert(y, pos, np.nan)
    return x,y

# Converts from ellipital orbit parameters to 2D orbit parameters
def ellipse_to_xy(a, e, theta, theta_E, GM_star):
    r = a*(1-e**2) / (1 + e*np.cos(theta - theta_E))
    u = -GM_star/(2*a)
    h = np.sqrt(GM_star*a*(1-e**2))
    v = np.sqrt(2*(u + GM_star/r))
    alpha = theta + np.arcsin(h/(r*v))
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    v_x = v*np.cos(alpha)
    v_y = v*np.sin(alpha)
    return x, v_x, y, v_y

# Converts from 2D orbit parameters to ellipital orbit parameters
def xy_to_ellipse(x, v_x, y, v_y, GM_star):
    r = np.sqrt(x**2 + y**2)
    v_squared = v_x**2 + v_y**2
    u = v_squared/2 - GM_star/r
    a = -GM_star/(2*u)
    h = x*v_y - y*v_x
    #GM_star*a*(1-e**2) = h^2
    #1-e**2 = h**2/(GM_star*a)
    e = np.sqrt(1. - h**2/(GM_star*a))
    theta = np.arctan2(y,x)
    #1 + e*cos(theta-theta_E) = a*(1-e**2)/r
    #e*cos(theta-theta_E) = a*(1-e**2)/r - 1
    #cos(theta-theta_E) = (a*(1-e**2)/r - 1)/e
    #theta-theta_E = arccos(a(1-e^2)/r - 1)
    acos_term = np.arccos((a*(1-e**2)/r - 1.)/e)
    if(x*v_x + y*v_y < 0):
        theta_E = theta + acos_term
    else:
        theta_E = theta - acos_term
    if(theta_E > np.pi): theta_E = theta_E - 2*np.pi
    if(theta_E < -np.pi): theta_E = theta_E + 2*np.pi
    if(e == 0): theta_E = 0
    return a, e, theta, theta_E