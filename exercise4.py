# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 09:04:34 2013

@author: c0918140
"""

import matplotlib.pyplot as plt
import numpy

method = raw_input("Please enter a method to use (euler, runge-kutta-order-2, runge-kutta-order-4)\n")

ftime0 = 0.
ftime1 = 1.
finit = 1.

def f(x):
    return x
    
def f2(x):
    return -x


def error(numx, anax):
    imax = len(numx)
    i = 0
    y = []
    while i < imax:
        y.append(anax[i] - numx[i])
        i += 1
    return y
    
def euler(t0, t1, steps, init):
    print t1, t0, steps
    dt = (t1 - t0)/steps
    y = []
    t = []
    y.append(init)
    i = 0
    t.append(t0)
    while i < steps:
        y.append((y[i] + dt*f(y[i])))
        t.append(t[i]+ dt)
        i += 1  
    return t, y
    
def rungeKuttaO2(t0, t1, steps, init):
    dt = (t1-t0)/steps
    y = []
    t = []
    y.append(init)
    i = 0
    t.append(t0)
    while i < steps:
        k1 = f(y[i])
        k2 = f(y[i] + k1*dt)
        y.append(y[i] + (dt/2.)*(k1 + k2))
        t.append(t[i] + dt)
        i += 1    
    return t, y

def rungeKuttaO4(t0, t1, steps, init):
    dt = (t1-t0)/steps
    y = []
    t = []
    y.append(init)
    i = 0
    t.append(t0)
    while i < steps:
        k1 = f(y[i])
        k2 = f(y[i] + k1*(dt/2.))
        k3 = f(y[i] + k2*(dt/2.))
        k4 = f(y[i] + k3*(dt))
        y.append(y[i] + (dt/6.)*(k1 + (2.*k2) + (2.*k3) + k4))
        t.append(t[i] + dt)
        i += 1    
    return t, y

def eulerConvergence(t0, t1, steps, init):
    time1, result1 = euler(t0, t1, steps, init)
    time2, result2 = euler(t0, t1, 2*steps, init)
    time3, result3 = euler(t0, t1, 4*steps, init)
    
    r1 = numpy.array(result1)    
    r2 = numpy.array(result2)    
    r3 = numpy.array(result3)
    
    diff1 = r1 - r2[::2]
    diff2 = r2[::2] - r3[::4]
    return time1, diff1, diff2

def solver(t0, t1, method, init):
    errortext = "Unknown method!"
    if method == "euler":
        #ordinary solve
        t, num1 = euler(t0, t1, 10., init)
        t2, num2 = euler(t0, t1, 50., init)
        ana1 = numpy.exp(t)
        plt.figure("x' = x")
        plt.plot(t, num1, "r+", label = "Numerical, 10 steps")
        plt.plot(t2, num2, "b+", label = "Numerical, 50 steps")
        plt.plot(t, ana1, color="r", label = "Analytical, exp(t)")
        plt.xlabel("time")
        plt.ylabel("y(t)")
        plt.legend()
        
        #convergence
        convt1, conv1 = euler(t0, t1, 100., init)
        convt2, conv2 = euler(t0, t1, 200., init)
        anaconv1 = numpy.exp(convt1)
        anaconv2 = numpy.exp(convt2)
        
        convdiff1 = conv1 - anaconv1
        convdiff2 = 2.0*(conv2 - anaconv2)
        
        plt.figure("Convergence error")
        plt.plot(convt1, convdiff1, "r", label = "Error, 100 steps")
        plt.plot(convt2, convdiff2, "b", label = "Error, 200 steps")
        plt.legend()
        
        time1, d1, d2 = eulerConvergence(t0, t1, 40., init)
        time2, d3, d4 = eulerConvergence(t0, t1, 80., init)
        time3, d5, d6 = eulerConvergence(t0, t1, 160., init)        
        
        plt.figure("Self-Convergence")
        plt.plot(time1,d1,'r-', label = '(dt=1/40) - (dt=1/80)')
        plt.plot(time1,2.*d2,'b-', label = '(dt=1/80) - (dt=1/160) \n scaled for 1st-order convergence')
        plt.xlabel('t')
        plt.ylabel('Differences in y(t)')
        #plt.plot(t, error(numx, anax), color="y")
    elif method == "runge-kutta-order-2":
        #ordinary solve
        t, num1 = rungeKuttaO2(t0, t1, 10., init)
        t2, num2 = rungeKuttaO2(t0, t1, 50., init)
        ana1 = numpy.exp(t)
        plt.figure("x' = x")
        plt.plot(t, num1, "r+", label = "Numerical, 10 steps")
        plt.plot(t2, num2, "b+", label = "Numerical, 50 steps")
        plt.plot(t, ana1, color="r", label = "Analytical, exp(t)")
        plt.xlabel("time")
        plt.ylabel("y(t)")
        plt.legend()
        
        #convergence
        convt1, conv1 = rungeKuttaO2(t0, t1, 100., init)
        convt2, conv2 = rungeKuttaO2(t0, t1, 200., init)
        anaconv1 = numpy.exp(convt1)
        anaconv2 = numpy.exp(convt2)
        
        convdiff1 = conv1 - anaconv1
        convdiff2 = 2.0*(conv2 - anaconv2)
        
        plt.figure("Convergence error")
        plt.plot(convt1, convdiff1, "r", label = "Error, 100 steps")
        plt.plot(convt2, convdiff2, "b", label = "Error, 200 steps")
        plt.legend()
    elif method == "runge-kutta-order-4":
        #ordinary solve
        t, num1 = rungeKuttaO4(t0, t1, 10., init)
        t2, num2 = rungeKuttaO4(t0, t1, 50., init)
        ana1 = numpy.exp(t)
        plt.figure("x' = x")
        plt.plot(t, num1, "r+", label = "Numerical, 10 steps")
        plt.plot(t2, num2, "b+", label = "Numerical, 50 steps")
        plt.plot(t, ana1, color="r", label = "Analytical, exp(t)")
        plt.xlabel("time")
        plt.ylabel("y(t)")
        plt.legend()
        
        convt1, conv1 = rungeKuttaO4(t0, t1, 100., init)
        convt2, conv2 = rungeKuttaO2(t0, t1, 100., init)
        convt3, conv3 = euler(t0, t1, 100., init)         
        anaconv1 = numpy.exp(convt1)
        
        convdiff1 = abs(conv1 - anaconv1)
        convdiff2 = abs(conv2 - anaconv1)
        convdiff3 = abs(conv3 - anaconv1)
       
        plt.figure("Error comparison")
        plt.plot(convt1, convdiff1, "r", label = "RK Order 4, 100 steps")
        plt.plot(convt2, convdiff2, "b", label = "RK Order 2, 100 steps")
        plt.plot(convt3, convdiff3, "g", label = "Euler, 100 steps")
        plt.legend(loc = "lower left")        
    
    else:
        print errortext
    
solver(ftime0, ftime1, method, finit)