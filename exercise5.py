# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 09:31:42 2013

@author: c0918140
"""
from numpy import *
from math import *
from matplotlib.pyplot import *

def gaussian(x):
    y = exp(-(x - 0.5)**2./(2.*0.1**2.))
    return y
    
def gaussAnalytic(x):
    z = exp(-(mod(x - 0.5, 1)**2./(0.1**2.)))
    return z

def solver(v, startTime, endTime, startSpace, endSpace, spaceSteps):

    dx = (endSpace - startSpace)/spaceSteps
    
    dt = dx/v
    
    timeSteps = (endTime - startTime)/dt    
    
    solution = zeros((timeSteps + 1, spaceSteps + 2))

    i = 1
    
    while i < spaceSteps + 1:
        solution[0, i] = gaussian(dx*i)
        i += 1

    solution[0, 0] = solution[0, spaceSteps]
    solution[0, spaceSteps + 1] = solution[0, 1]     
     
    i = 1        
    n = 0    
    
    while n < timeSteps:
        i = 1
        while i < spaceSteps + 1:
            solution[n + 1, i] = solution[n, i - 1]
            i += 1
        n += 1
        solution[n, 0] = solution[n, spaceSteps]
        solution[n, spaceSteps + 1] = solution[n, 1]
        
    k = 0
    x = linspace(startSpace, endSpace, spaceSteps)
    figure("Velocity: " + str(v) + " | Spatial steps: " + str(spaceSteps))    
    while k < timeSteps:
        y = solution[k]

        plot(x, y[1:spaceSteps + 1])
        k += 1    
    return 0

solver(0.5, 0., 2., 0., 1., 20.)
solver(1.0, 0., 2., 0., 1., 20.)
solver(2.0, 0., 2., 0., 1., 20.)

def L2Norm(v, startTime, endTime, startSpace, endSpace, spaceSteps):
    
    dx = (endSpace - startSpace)/spaceSteps
    
    dt = dx/v
    
    timeSteps = (endTime - startTime)/dt    
    
    solution = zeros((timeSteps + 1, spaceSteps + 2))
    solution2 = zeros((timeSteps + 1, spaceSteps + 2))

    i = 1
    
    while i < spaceSteps + 1:
        solution[0, i] = gaussian(dx*i)
        i += 1
        
    i = 1
    
    while i < spaceSteps + 1:
        solution2[0, i] = gaussAnalytic(dx*i)
        i += 1
     

    solution[0, 0] = solution[0, spaceSteps]
    solution[0, spaceSteps + 1] = solution[0, 1]

    solution2[0, 0] = solution2[0, spaceSteps]
    solution2[0, spaceSteps + 1] = solution2[0, 1]         
     
    i = 1        
    n = 0
    difference = 0
    norm = zeros(timeSteps)
    
    while n < timeSteps:
        i = 1
        while i < spaceSteps + 1:
            solution[n + 1, i] = solution[n, i - 1]
            solution2[n + 1, i] = solution2[n, i - 1]
            difference += (solution[n, i] - solution2[n, i])
            
            i += 1
        norm[n] = pow((difference*difference)/spaceSteps, 0.5)
        n += 1
        solution[n, 0] = solution[n, spaceSteps]
        solution[n, spaceSteps + 1] = solution[n, 1]
        solution2[n, 0] = solution2[n, spaceSteps]
        solution2[n, spaceSteps + 1] = solution2[n, 1]
    
    figure("L2 Norm, "+str(spaceSteps)+" spatial steps")
    t = linspace(startTime, endTime, timeSteps)
    plot(t, norm)
    xlabel("Time")
    ylabel("norm")
    
    return 0
    
L2Norm(2., 0., 5., 0., 1., 10.)
L2Norm(2., 0., 5., 0., 1., 100.)
L2Norm(2., 0., 5., 0., 1., 1000.)