# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 09:04:34 2013

@author: c0918140
"""

import matplotlib.pyplot as plt
import numpy

#Euler method for solving a first order time dependent ODE
function = raw_input("Please enter an ODE to solve, dx/dt. Do not include the differential operator\n")

method = raw_input("Please enter a method to use (euler, runge-kutta)\n")

time0 = input("Please enter an initial time as an integer\n")

time1 = input("Please enter an ending time as an integer\n")

timesteps = input("Please input the number of steps to use as an integer\n")

initialcondition = input("Please input the value of the ODE at t = 0\n")

fsteps = float(timesteps)
ftime0 = float(time0)
ftime1 = float(time1)
finit = float(initialcondition)

def f(x):
    y = eval(function)
    return y

#Overcomplicated time function. numpy linspace achieves the same result in one line
def time(t0, t1, steps):
    t = []
    dt = (t1-t0)/steps
    t.append(t0)
    i = 0
    while i < steps:
        t.append(t[i] + dt)
        i += 1
    return t
    
#Finds the error by taking the difference between the analytical and numerical solutions.
#Could be one line as return anax-numx if the two are numpy arrays
def error(numx, anax):
    imax = len(numx)
    i = 0
    y = []
    while i < imax:
        y.append(anax[i] - numx[i])
        i += 1
    return y
    
#Euler solver. Extremely simple. Could also be reduced to a couple of array multiplications.
def euler(t0, t1, steps, init):
    dt = (t1-t0)/steps
    y = []
    y.append(init)
    i = 0
    t = t0
    while i < steps:
        y.append(y[i] + dt*f(y[i]))
        t += dt
        i += 1    
    return y

#Analytical solution, exponential, again could be a single line
def analytic(t0, t1, steps, init):
    dt = (t1-t0)/steps
    y = []
    i = 0
    t = t0
    while i <= steps:
        y.append(numpy.exp(t))
        t += dt
        i += 1
    return y

#Run everything and plot it. Has basic error handling.
def solver(t0, t1, method, steps, init):
    errortext = "Unknown method!"
    if method == "euler":     
        print "The results are:\n"
        numx = euler(t0, t1, steps, init)
        anax = analytic(t0, t1, steps, init)
        t = time(t0, t1, steps)
        plt.figure("x' = " + function)
        plt.plot(t, numx)
        plt.plot(t, anax, color="r")
        plt.figure("Error")
        plt.plot(t, error(numx, anax), color="y")
        plt.show()
    elif method == "runge-kutta":
        print errortext
    else:
        print errortext
    
solver(ftime0, ftime1, method, fsteps, finit)