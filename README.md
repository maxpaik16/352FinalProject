# 352FinalProject
Jason Phelan, Orion Forowycz, Valeriia Rohoza, and Max Paik's 352 final project on the stability of planetary orbits

## Description of Contents

Many files are copied over in multiple places to make import statements work. Apologies for the clutter.

ode.c is the source file containing our integrators

ode.h is the corresponding header file

odesolver.py contains Python wrappers over our integrators

helpers.py contains helper functions used in other scripts (indexing or coordinate conversion, for example)

interfunc.py contains right hand side equations for our integrators to solve

in Solver

  both Jupyter Notebooks contain code that we used to test our integrators.
  
in IntegratorFunction

  func_testing.ipynb contains code used to test our Python scripts
  
in Driver

  \tsimulate.py is the script use to simulate many different systems
  
  results.py plots the results of simulate.py
  
  model.py is an attempt at a neural network to predict stability
  
  
