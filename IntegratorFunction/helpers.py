"""
Contains helper functions for misc orbital calculations

To use write:
from helpers import *
"""

import numpy as np
import math

# Converts angle in radians into range (-pi,pi]
def angle_wrap(a):
    b = a/(2.*np.pi)
    b -= np.floor(b)
    b -= (b>0.5)
    b *= 2.*np.pi
    return(b)

def orbital_period(a, GM_star):
    """
    Gives orbital period for orbit of semi-major axis "a" around a star (or any body being orbited)
    """
    return 2*np.pi*np.sqrt(a**3/GM_star)

# Removes horizontal jumps in plots of x and y
# For example for chaotic phase space plots
def remove_horizontal_discontinuities(x,y):
    pos = np.where(np.abs(np.diff(x)) >= 6)[0]+1
    x = np.insert(x, pos, np.nan)
    y = np.insert(y, pos, np.nan)
    return x,y

def ellipse_to_xy(a,e,theta,thetaE):
    """
    Takes the particle's position relative to an ellipse and parameters of the ellipse a,e,theta,theta_E.
    This function returns the Cartesian variables x,V_x,y,V_y.
    
    Returns x,Vx,y,Vy
    """
    # radius using angle theta
    r = a * (1 - e**2) / (1 + e * np.cos(theta - thetaE))
    
    # angular momentum per mass
    h = 2. * np.pi * np.sqrt(np.abs(a * (1. - e **2)))
    
    # energy per mass
    u = - 2. * (np.pi ** 2) / a 
    
    # speed of the particle
    V = np.sqrt(np.abs(2. * u + 8. * (np.pi ** 2) / r)) 
    
    # let Vx = V cos alpha, Vy = V sin alpha
    # buff = alpha - theta
    # when the radial velocity is positive (the planet goes from its periapse to apoapse = sin(theta-theta_E) > 0)
    # alpha - theta should be less then pi/2
    buff_sin = h / (r * V)
    buff_sin[buff_sin < -1.] = -1.
    buff_sin[buff_sin > 1.] = 1.
    #assert -1 <= buff_sin <= 1
    buff = np.pi*(np.sin(theta - thetaE) < 0.) + np.power(-1., np.sin(theta - thetaE) < 0.) * np.arcsin(buff_sin)
    alpha = theta + buff
        
    # x and y
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    
    # Vx and Vy
    Vx = V * np.cos(alpha)
    Vy = V * np.sin(alpha)
    
    return x,Vx,y,Vy

def xy_to_ellipse(x,Vx,y,Vy):
    """
    Takes the Cartesian variables.
    This function returns the particle's position relative to an ellipse and parameters of the ellipse.
    
    Returns a,e,theta,theta_E
    """
    # radius using x and y
    r = np.sqrt(x ** 2 + y ** 2)
    
    # speed of the particle
    V = np.sqrt(Vx ** 2 + Vy ** 2)
    
    # angular momentum per mass
    h = x * Vy - y * Vx
    
    # energy per mass
    u = (V ** 2) / 2. - 4. * (np.pi ** 2) / r
    
    # semi-major axis
    a = -2. * ((np.pi) ** 2) / u
    
    # eccentricity of the elliptical orbit, added absolute value
    e = np.sqrt(np.abs(1 - ((h / (2. * np.pi)) ** 2 )/ a))
    
    # theta
    theta = np.arctan2(y,x)
    
    # theta_E, compute e*cos(theta - thetaE) first
    buff = a * (1. - e ** 2) / r - 1.
    
    # divide buff/e and output 0 if it is a circular orbit
    buff_cos = np.divide(buff, e, out=np.zeros_like(buff), where=(e > np.power(10.,-5.)))
    
    #to make sure that arccos takes values less than 1 and greater than -1
    buff_cos[buff_cos < -1.] = -1.
    buff_cos[buff_cos > 1.] = 1.
    
    delta = np.arccos(buff_cos)
    
    # change the sign if the radial velocity is negative
    delta *= np.power(-1.,(x * Vx + y * Vy) < 0.)
    thetaE = theta - delta
    
    # set thetaE to 0 if it is a circular orbit
    thetaE *= (e > np.power(10.,-5.))
    
    # fix to add 2pi or subtract 2pi if thetaE isn't between -pi and pi
    thetaE -= (thetaE > np.pi) * 2 * np.pi
    thetaE += (thetaE < -np.pi) * 2 * np.pi
    
    return a,e,theta,thetaE

"""

# Converts from ellipital orbit parameters to 2D orbit parameters
def ellipse_to_xy(a, e, theta, theta_E):
    GM_star = 4*np.pi**2
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

"""