
# python script for running simulations

import csv
import ctypes
from ctypes import *
from numpy.ctypeslib import ndpointer
import numpy as np
import matplotlib.pyplot as plt
import numba
import math
import random

import importlib
import odesolver
importlib.reload(odesolver)
from odesolver import *
import helpers
importlib.reload(helpers)
from helpers import *
import interfunc
importlib.reload(interfunc)
from interfunc import *

GM_S = 39.4229 #AU^3.yr^{-2} solar mass parameter

a_0 = 0.39
e_0 = 0.206
theta_E_0 = -3*np.pi/4
theta_0 = theta_E_0
x_0,v_x_0,y_0,v_y_0 = ellipse_to_xy(a_0, e_0, theta_0, theta_E_0)
initial_mercury = [x_0,v_x_0,y_0,v_y_0]

GM_Sun = 39.4229 #AU^3.yr^{-2} solar mass parameter

# Earth's standard gravitational parameter
GM_Ear = 0.00011841685 #AU^3/yr^2

# https://nssdc.gsfc.nasa.gov/planetary/factsheet/planet_table_ratio.html
GM_Mer = GM_Ear*0.0553
GM_Ven = GM_Ear*0.815
GM_Mar = GM_Ear*0.107
GM_Jup = GM_Ear*317.8
GM_Sat = GM_Ear*95.2
GM_Ura = GM_Ear*14.5
GM_Nep = GM_Ear*17.1

# run simulations
for i in range(100):
    # add a bit of variation to all planetary initial conditions
    init_Mer = np.array(ellipse_to_xy(random.gauss(0.3870993, .0001), random.gauss(0.20564, .0001), 0., 0.))
    init_Ven = np.array(ellipse_to_xy(random.gauss(0.723336, .0001), random.gauss(0.00678, .0001), 0., 0.))
    init_Ven[ind_v_y(0)] *= -1
    init_Ear = np.array(ellipse_to_xy(random.gauss(1.000003, .0001), random.gauss(0.01671, .0001), 0., 0.))
    init_Mar = np.array(ellipse_to_xy(random.gauss(1.52371, .0001), random.gauss(0.09339, .0001), 0., 0.))
    init_Jup = np.array(ellipse_to_xy(random.gauss(5.2029, .0001), random.gauss(0.0484, .0001), 0., 0.))
    init_Sat = np.array(ellipse_to_xy(random.gauss(9.537, .0001), random.gauss(0.0539, .0001), 0., 0.))
    init_Ura = np.array(ellipse_to_xy(random.gauss(19.189, .0001), random.gauss(0.04726, .0001), 0., 0.))
    init_Nep = np.array(ellipse_to_xy(random.gauss(30.0699, .0001), random.gauss(0.00859, .0001), 0., 0.))
    n_planets = 8

    # add more significant noise to Mercury's initial conditions
    init_Mer = list(map(lambda x : x + random.gauss(0, .05), init_Mer))
    params = [GM_Sun, n_planets, GM_Mer, GM_Ven, GM_Ear, GM_Mar, GM_Jup, GM_Sat, GM_Ura, GM_Nep]
    noisy_planets = np.concatenate((init_Mer,init_Ven,init_Ear,init_Mar,init_Jup,init_Sat,init_Ura,init_Nep))

    a_Nep = 30.0699
    total_time = 1000 * orbital_period(a_Nep,GM_Sun) # 1 Neptune period
    step_size = orbital_period(a_0,GM_S)/100 # 1/100 of Mercury period
    n_steps = int(total_time/step_size)

    # run simulations
    t,sol_untransposed = solve_ode(func_n_body,[0.,total_time], n_steps, noisy_planets, args=params, method="Yoshida4")
    sol = sol_untransposed.T

    # check to see whether or not Mercury has been ejected from its orbit
    ejected = 0
    if np.amax(sol[ind_x(0)]) > 10:
        ejected = 1

    # write results to CSV file
    with open('training.csv', mode='a') as file:
        writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        data = np.ndarray.tolist(noisy_planets)
        data.append(ejected)
        writer.writerow(data)
