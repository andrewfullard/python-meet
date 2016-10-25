# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 09:04:34 2013

@author: c0918140 Andrew Fullard
"""

#==============================================================================
# 
#  1) Splitting the equation y" = -w*w*y:
#      
#      Let u = y, v = y'
#      
#      => u' = v
#      => v' = -w*w*y
# 
#     Taking an analytic solution of y(t) = cos(wt)
#     
#     => y(t) = cos(2*pi*t) for w = 2pi
#     
#     y' = -w*sin(wt)
#     
#     y" = -w*w*cos(wt) = -w*w*y 
# 
#     y(0) = 1
#     
#     y(0) = cos(0) = 1
#     
#     y'(0) = 0
#     
#     y'(0) = -w*sin(0) = 0
#    
#     => cos(wt) is the solution
#==============================================================================
    
#import necessary modules
import matplotlib.pyplot as plt
import numpy

#define differentials
function1 = "y"

function2 = "-4*3.1415*3.1415*y"

#user input for the method to use
method = raw_input("Please enter a method to use (euler, runge-kutta)\n")

#time info (initial, final, number of steps)
ftime0 = 0.0

ftime1 = 1.0

fsteps = 1.0

#maximum power of 2 to divide the number of steps by
stepmultiple = 10

#initial conditions for y(0), y'(0)
finit1 = 1.0
finit2 = 0.0

#evaluates function 1 to be machine-readable
def y1(y):
    y = eval(function1)
    return y

#evaluates function 2 to be machine-readable
def y2(y):
    y = eval(function2)
    return y

#==============================================================================
# 2) The Euler method is applied to the two differential equations to find the
# solution. It requires the use of an intermediate step that depends on w.
#==============================================================================

#euler solver
def euler(t0, t1, steps, init1, init2):
    #determine time step
    dt = (t1 - t0)/steps
    #initialise lists; w is for w' = z, z is for z' = -w*w*y
    u = []
    v = []
    t = []
    #initial conditions in the lists
    v.append(init1)
    u.append(init2) 
    #loop index
    i = 0
    #initial time
    t.append(t0)
    while i < steps:
        #euler for two ODEs
        v.append(v[i] + dt*u[i])
        u.append(u[i] + dt*y2(v[i]))
        t.append(t[i] + dt)
        i += 1 
    return t, v    
       
#==============================================================================
# 2) Self-convergence testing further confirms the first-order convergence       
#==============================================================================
       
#euler convergence testing
def eulerConvergence(t0, t1, steps, init1, init2):
    #runs euler method for dt, dt/2, dt/4
    time1, result1 = euler(t0, t1, steps, init1, init2)
    time2, result2 = euler(t0, t1, 2*steps, init1, init2)
    time3, result3 = euler(t0, t1, 4*steps, init1, init2)
    
    #array results so they can be used for calculations
    r1 = numpy.array(result1)    
    r2 = numpy.array(result2)    
    r3 = numpy.array(result3)
    
    #calculate differences
    diff1 = r1 - r2[::2]
    diff2 = r2[::2] - r3[::4]
    return time1, diff1, diff2   
    
#error calculation(numerical solution, analytical solution, number of solutions to test)
def errorEuler(stepmultiple):
    #initialise list
    y = []
    #loop index
    i = 0
    timestep = []
    print "Timestep, powers of 2 (Euler)"
    while i <= stepmultiple:
        #print progress through powers of 2
        print i
        #decrease size of timestep by 1/2^i as i runs from 0 to the power defined by stepmultiple
        timestep.append(fsteps*(pow(2, i)))
        #calculate the euler and analytic results for the timesteps
        t, numx = euler(ftime0, ftime1, timestep[i], finit1, finit2)
        #get error value at the end of the time limit (only part of interest)
        nmax = len(numx) - 1
        #get difference between numerical and analytical results. Analytic result is known to be
        #cos(2*pi) = 1 so we do not need to calculate it every time
        diff = 1.0 - numx[nmax]
        #add the modulus of the difference to the error list
        y.append(abs(diff))
        #increment ith power of 2
        i += 1
    #return the timesteps and error list
    return timestep, y

#==============================================================================
# 5) 4th-order Runge-Kutta method, adapted for two ODEs in a similar way to the Euler
# method.
#==============================================================================
   
#runge-kutta solver
def rungeKutta(t0, t1, steps, init1, init2):
    #determine time step
    dt = (t1-t0)/steps
    #lists; w is w' = y, z is z' = -4*pi*pi*y
    u = []
    v = []
    t = []
    v.append(init1)
    u.append(init2)
    i = 0
    t.append(t0)
    while i < steps:
        #solver for 2nd order ODE; calculates constants based on the split 1st order ODEs
        #kwX are the constants for w based on z and kzX are the constants for z based on w
        #z is the final result based of w and all the constants
        ku1 = y2(u[i])
        kv1 = y1(v[i])
        ku2 = y2(u[i] + kv1*(dt/2.))
        kv2 = y1(v[i] + ku1*(dt/2.))
        ku3 = y2(u[i] + kv2*(dt/2.))
        kv3 = y1(v[i] + ku2*(dt/2.))
        ku4 = y2(u[i] + kv3*(dt))
        kv4 = y1(v[i] + ku3*(dt))
        u.append(u[i] + (dt/6.)*(kv1 + (2.*kv2) + (2.*kv3) + kv4))
        v.append(v[i] + (dt/6.)*(ku1 + (2.*ku2) + (2.*ku3) + ku4))
        t.append(t[i] + dt)
        i += 1    
    return t, v

#error calculation(numerical solution, analytical solution, number of solutions to test)
def errorRunge(stepmultiple):
    #initialise list
    y = []
    #loop index
    i = 0
    timestep = []
    print "Timestep, powers of 2 (Runge-Kutta)"
    while i <= stepmultiple:
        #print progress through powers of 2
        print i
        #decrease size of timestep by 1/2^i as i runs from 0 to the power defined by stepmultiple
        timestep.append(fsteps*(pow(2, i)))
        #calculate the euler and analytic results for the timesteps
        t, numx = rungeKutta(ftime0, ftime1, timestep[i], finit1, finit2)
        #get error value at the end of the time limit (only part of interest)
        nmax = len(numx) - 1
        #get difference between numerical and analytical results. Analytical result is known to be
        #cos(2*pi) = 1 so we do not need to calculate it every time
        diff = 1.0 - numx[nmax]
        #add the modulus of the difference to the error list
        y.append(abs(diff))
        #increment ith power of 2
        i += 1
    #return the timesteps and error list
    return timestep, y

#input reader
def solver(t0, t1, method, steps, init1, init2, stepmult):
    #default error message
    errortext = "Unknown method! Enter \"euler\" or \"runge-kutta\""
    #euler input
    if method == "euler":
        #calculate solution with 50 steps, compare to analytical solution
        t, num1 = euler(t0, t1, 50., init1, init2)
        t = numpy.array(t)
        ana1 = numpy.cos(2*3.1415*t)
        plt.figure("y\" = -w*w*y")
        plt.plot(t, num1, "b+", label = "Numerical, 50 steps")
        plt.plot(t, ana1, color="r", label = "Analytical, cos(wt)")
        plt.xlabel("time")
        plt.ylabel("y(t)")
        plt.legend()
        
#==============================================================================
#         2) Here the Euler method is tested against the analytical result and 
#         displays first order convergence because the curves lie on top of one
#         another.
#==============================================================================
        
        #first order convergence test for 100 and 200 steps
        convt1, conv1 = euler(t0, t1, 100., init1, init2)
        convt2, conv2 = euler(t0, t1, 200., init1, init2)
        convt1 = numpy.array(convt1)
        convt2 = numpy.array(convt2)
        anaconv1 = numpy.cos(2*3.1415*convt1)
        anaconv2 = numpy.cos(2*3.1415*convt2)
        
        convdiff1 = conv1 - anaconv1
        convdiff2 = 2.0*(conv2 - anaconv2)
        
        #plotting
        plt.figure("Convergence error")
        plt.plot(convt1, convdiff1, "r", label = "Error, 100 steps")
        plt.plot(convt2, convdiff2, "b", label = "Error, 200 steps")
        plt.xlabel("t")
        plt.legend()

#==============================================================================
#       2) Self-convergence testing further confirms the first-order convergence
#           as the curves lie almost on top of one another again.      
#==============================================================================
       
        #Self-convergence test
        time1, d1, d2 = eulerConvergence(t0, t1, 40., init1, init2)
        time2, d3, d4 = eulerConvergence(t0, t1, 80., init1, init2)
        time3, d5, d6 = eulerConvergence(t0, t1, 160., init1, init2)        
        
        #plotting
        plt.figure("Self-Convergence")
        plt.plot(time1,d1,'r-', label = '(dt=1/40) - (dt=1/80)')
        plt.plot(time1,2.*d2,'b-', label = '(dt=1/80) - (dt=1/160) \n scaled for 1st-order convergence')
        plt.xlabel('time')
        plt.ylabel('Differences in y(t)')
        plt.legend(loc = "best")
        
#==============================================================================
#         2) The ratio of the convergence behaves badly in multiple places, 
#         because the ratios are equivalent to a tangent curve (sin/cos).
#         The asymptotes are occuring at pi and 3pi/4 (i.e. w*0.5 and w*0.75) which
#         is where tan(wt) tends to infinity
#==============================================================================
        
        #convergence ratios
        plt.figure("Convergence ratios")
        plt.plot(time1[1:],d1[1:]/d2[1:], 'r-', label='dt = 1/40, 1/80, 1/160')
        plt.plot(time2[1:],d3[1:]/d4[1:], 'b-', label='dt = 1/80, 1/160, 1/320')
        plt.plot(time3[1:],d5[1:]/d6[1:], 'g-', label='dt = 1/160, 1/320, 1/640')
        plt.title('Self-convergence test (as ratios) for y\" = w*w*y')
        plt.xlabel('t')
        plt.ylabel('Ratios of differences in y(t)')
        plt.legend(loc='lower left')
        
#==============================================================================
#         3) Error calculations for various timesteps dt = 1/2...1/2^18.
#         Plotted as log-log. Optimised by recognising that the value is always 1
#         at t = 1 for the analytic solution, so it does not need to be re-calculated
#         each time.
#==============================================================================
        
        #plot log-log error vs number of timesteps
        plt.figure("Error")
        timestep, errval = errorEuler(stepmult)
        plt.loglog(timestep, errval, color="y")           
        plt.xlabel("Number of Time Steps")
        plt.ylabel("Error in Method")
        
#==============================================================================
#         4) The plot shows that for dt > 1/16, the convergence of the method is
#         not linear. Otherwise first-order convergence is displayed. Where the 
#         convergence is not clean, the method is unstable.        
#==============================================================================
        
    #runge-kutta input
    elif method == "runge-kutta":

#==============================================================================
#         5) Runge-Kutta method. Much closer results to the analytic solution.
#==============================================================================
        
        #ordinary solve
        t, num1 = rungeKutta(t0, t1, 50., init1, init2)
        t = numpy.array(t)
        ana1 = numpy.cos(2*3.1415*t)
        plt.figure("y\" = -w*w*y")
        plt.plot(t, num1, "r+", label = "Numerical, 50 steps")
        plt.plot(t, ana1, color="r", label = "Analytical, cos(wt)")
        plt.xlabel("time")
        plt.ylabel("y(t)")
        plt.legend()
        
#==============================================================================
#         5) The method shows fourth-order convergence as doubling the number of
#         steps decreases error by a factor of ~8.
#==============================================================================
        
        #convergence
        convt1, conv1 = rungeKutta(t0, t1, 100., init1, init2)
        convt2, conv2 = rungeKutta(t0, t1, 200., init1, init2)
        convt1 = numpy.array(convt1)
        convt2 = numpy.array(convt2)
        anaconv1 = numpy.cos(2*3.1415*convt1)
        anaconv2 = numpy.cos(2*3.1415*convt2)
        
        convdiff1 = conv1 - anaconv1
        convdiff2 = 2.0*(conv2 - anaconv2)
       
        plt.figure("Convergence")
        plt.plot(convt1, convdiff1, "r", label = "Error, 100 steps")
        plt.plot(convt2, convdiff2, "b", label = "Error, 200 steps")
        plt.xlabel("time")
        plt.legend(loc = "best")
        
#==============================================================================eul
#         5) The error change with timestep is much steeper for the Runge-Kutta
#         method, displaying 4th-order convergence. Furthermore, for dt > 1/256
#         the error no longer decreases. This is because for such small timesteps
#         float precision is no longer enough to resolve the differences between
#         the numerical and analytical solutions.
#==============================================================================
        
        #plot log-log error vs number of timesteps
        plt.figure("Error")
        timestep, erreul = errorEuler(stepmult)
        timestep, errrun = errorRunge(stepmult)
        plt.loglog(timestep, erreul, color="r", label = "Euler method")
        plt.loglog(timestep, errrun, color="b", label = "Runge-Kutta method")
        plt.xlabel("Number of Time Steps")
        plt.ylabel("Error in Method")          
        plt.legend()
        
    #other input
    else:
        print errortext
        
#run main function
solver(ftime0, ftime1, method, fsteps, finit1, finit2, stepmultiple) 
