# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 13:08:45 2013

@author: Andrew
"""

from pylab import *

#Euler Method
def Euler(u, udot, udash, dx, dt):
    v = zeros(len(u))
    w = zeros(len(udot))
    z = zeros(len(udash))
    for i in range(1, len(w) - 1):
        v[i] = u[i] + dt*udot[i]
        w[i] = udot[i] + dt*f(udash, i, dx)
        z[i] = udash[i] + dt*f(udot, i, dx)
    return v, w, z

#Runge-Kutta 4th order method
def RK4Step(u, udot, udash, dx, dt):
    v = zeros(len(u))
    w = zeros(len(udot))
    z = zeros(len(udash))
    kv1 = zeros(len(u))
    kv2 = zeros(len(u))
    kv3 = zeros(len(u))
    kv4 = zeros(len(u))
    kw1 = zeros(len(u))
    kw2 = zeros(len(u))
    kw3 = zeros(len(u))
    kw4 = zeros(len(u))
    kz1 = zeros(len(u))
    kz2 = zeros(len(u))
    kz3 = zeros(len(u))
    kz4 = zeros(len(u))

    for i in range(1, len(w) - 1):
        kv1[i] = udot[i]
        kw1[i] = f(udash, i, dx)
        kz1[i] = f(udot, i, dx)
    kv1 = ApplyBoundaryConditions(kv1)
    kw1 = ApplyBoundaryConditions(kw1)
    kz1 = ApplyBoundaryConditions(kz1)
   
    for i in range(1, len(w) - 1):
        kv2[i] = udot[i] + kw1[i]*(dt/2.)
        kw2[i] = f(udash + kz1*(dt/2.), i, dx)
        kz2[i] = f(udot + kw1*(dt/2.), i, dx)
    kv2 = ApplyBoundaryConditions(kv2)
    kw2 = ApplyBoundaryConditions(kw2)
    kz2 = ApplyBoundaryConditions(kz2)        
    
    for i in range(1, len(w) - 1):  
        kv3[i] = udot[i] + kw2[i]*(dt/2.)
        kw3[i] = f(udash + kz2*(dt/2.), i, dx)
        kz3[i] = f(udot + kw2*(dt/2.), i, dx)
    kv3 = ApplyBoundaryConditions(kv3)
    kw3 = ApplyBoundaryConditions(kw3)
    kz3 = ApplyBoundaryConditions(kz3)        
    
    for i in range(1, len(w) - 1):    
        kv4[i] = udot[i] + kw3[i]*dt
        kw4[i] = f(udash + kz3*dt, i, dx)
        kz4[i] = f(udot + kw3*dt, i, dx)
    
    kv4 = ApplyBoundaryConditions(kv4)
    kw4 = ApplyBoundaryConditions(kw4)
    kz4 = ApplyBoundaryConditions(kz4)
    
 
    v = u + (dt/6.)*(kv1 + (2.*kv2) + (2.*kv3) + kv4)
    w = udot + (dt/6.)*(kw1 + (2.*kw2) + (2.*kw3) + kw4)
    z = udash + (dt/6.)*(kz1 + (2.*kz2) + (2.*kz3) + kz4)
    return v, w, z

#finite difference function    
def f(u, i, dx):
    return (u[i + 1] - u[i - 1])/(2.*dx)
 
#Applies boundary conditions (from example code by M Hannam)       
def ApplyBoundaryConditions(y):
    Np = len(y)-2
    y[0] = y[Np]
    y[Np+1] = y[1]
    return y

#Calculates the L2 Norm from two data sets, dat2 having 2 times the points
#of dat1. From example code by M. Hannam.
def CvgNorm(times1, dat1, dat2):
    norm = zeros(len(times1))
    for i in range(0, len(times1)):
       norm[i] = sqrt(sum((dat1[i,1:-1] - dat2[2*i,2:-1:2])**2)/len(norm))
    return norm
    
#Main solution loop, parts from example code by M. Hannam
def SolveAdvectionEquation(a, b, t0, t1, Nx, c, method):
    #Calculate dx and dt
    dx = (float) (b - a)/Nx
    dt = c*dx
    
    #get number of time steps    
    Nt = (t1 - t0)/dt

    # We have Nx+1 points for Nx intervals, plus an extra ghost point.
    x = linspace(a - dx, b, Nx + 2)
    times = linspace(t0, t1, Nt + 1)
    
    u = zeros((Nt + 1, Nx + 2))
    #udash space derivative, udot time derivative
    udot = zeros((Nt + 1, Nx + 2))
    udash = zeros((Nt + 1, Nx + 2))
    sigma2 = 2.*0.1**2   
    
    #Initial conditions
    u[0] = exp(-(x**2)/(sigma2))
    u[0] = ApplyBoundaryConditions(u[0])
    udash[0] = -100.*x*u[0]
    udash[0] = ApplyBoundaryConditions(udash[0])
    
    #Time evolution of solution
    for i in range(0, int(Nt)):
        if (method == 'euler'):
            u[i+1], udot[i+1], udash[i+1] = Euler(u[i], udot[i], udash[i], dx, dt)
        if (method == 'rk4'):
            u[i+1], udot[i+1], udash[i+1] = RK4Step(u[i], udot[i], udash[i], dx, dt)
        u[i+1] = ApplyBoundaryConditions(u[i+1])
        udot[i+1] = ApplyBoundaryConditions(udot[i+1])
        udash[i+1] = ApplyBoundaryConditions(udash[i+1])
    return times, x, u
    
#Various plotting for Euler method
times1, grid1, soln1 = SolveAdvectionEquation(-1., 1., 0., 1.2, 100., 0.5, 'euler')
times2, grid2, soln2 = SolveAdvectionEquation(-1., 1., 0., 1.2, 200., 0.5, 'euler')
times3, grid3, soln3 = SolveAdvectionEquation(-1., 1., 0., 1.2, 400., 0.5, 'euler')
times4, grid4, soln4 = SolveAdvectionEquation(-1., 1., 0., 1.2, 800., 0.5, 'euler')
times5, grid5, soln5 = SolveAdvectionEquation(-1., 1., 0., 1.2, 1600., 0.5, 'euler')
times6, grid6, soln6 = SolveAdvectionEquation(-1., 1., 0., 1.2, 400., 0.25, 'euler')
times7, grid7, soln7 = SolveAdvectionEquation(-1., 1., 0., 1.2, 800., 0.25, 'euler')
times8, grid8, soln8 = SolveAdvectionEquation(-1., 1., 0., 1.2, 1600., 0.25, 'euler')
times9, grid9, soln9 = SolveAdvectionEquation(-1., 1., 0., 1.2, 400., 0.75, 'euler')
times10, grid10, soln10 = SolveAdvectionEquation(-1., 1., 0., 1.2, 800., 0.75, 'euler')
times11, grid11, soln11 = SolveAdvectionEquation(-1., 1., 0., 1.2, 1600., 0.75, 'euler')

#RK method plotting
timesrk1, gridrk1, solnrk1 = SolveAdvectionEquation(-1., 1., 0., 1.2, 100., 0.5, 'rk4')
timesrk2, gridrk2, solnrk2 = SolveAdvectionEquation(-1., 1., 0., 1.2, 200., 0.5, 'rk4')
timesrk3, gridrk3, solnrk3 = SolveAdvectionEquation(-1., 1., 0., 1.2, 400., 0.5, 'rk4')
timesrk4, gridrk4, solnrk4 = SolveAdvectionEquation(-1., 1., 0., 1.2, 800., 0.5, 'rk4')

#L2 norm calculations
norm = CvgNorm(times1, soln1, soln2)
norm2 = CvgNorm(times3, soln3, soln4)
norm3 = CvgNorm(times2, soln2, soln3)
norm4 = CvgNorm(times4, soln4, soln5)
#c = 0.25
norm5 = CvgNorm(times6, soln6, soln7)
norm6 = CvgNorm(times7, soln7, soln8)
#c = 0.75
norm8 = CvgNorm(times9, soln9, soln10)
norm9 = CvgNorm(times10, soln10, soln11)

normrk4 = CvgNorm(timesrk1, solnrk1, solnrk2)
normrk42 = CvgNorm(timesrk3, solnrk3, solnrk4)

#Plotting, ignoring ghost points where relevant
figure()
for i in range(0, len(times1)):
    if times1[i] == 0 or times1[i] == 0.25 or times1[i] == 0.5:
        plot(grid1[1:], soln1[i,1:], label="t = " + str(times1[i]))
xlabel("x")
ylabel("u(x, t)")
title("Wave Equation Solution at t = 0, t = 0.25, t = 0.5, Euler Method")
legend(loc=0)

figure()
for i in range(0, len(timesrk1), 10):
    plot(gridrk1[1:], solnrk1[i,1:])
xlabel("x")
ylabel("u(x, t)")
title("Wave Equation Solution, RK4 Method")

figure()
for i in range(0, len(times1)):
    if times1[i] == 0.9 or times1[i] == 1.0 or times1[i] == 1.1:
        plot(grid1[1:], soln1[i,1:], label="t = " + str(times1[i]))
xlabel("x")
ylabel("u(x, t)")
title("Wave Equation Solution at t = 0.9, t = 1.0, t = 1.1, Euler Method")
legend(loc=0)

figure()
plot(times1, log(norm), label = "100-200 steps")
plot(times3, log(norm2), label = "400-800 steps")
plot(times4, log(norm4), label = "800-1600 steps")
xlabel("Time (s)")
ylabel("L2 Norm")
title("L2 Norm for the Euler Method")

figure()
plot(times1, log(norm), label="Euler Method 100-200")
plot(timesrk1, log(normrk4), label="Runge-Kutta 100-200")
plot(timesrk3, log(normrk42), label="Runge-Kutta 400-800")
xlabel("Time (s)")
ylabel("L2 Norm")
title("L2 Norm Comparing 100 and 200 steps")
legend(loc = 0)

figure()
plot(times2, log(norm3), label = "c = 0.5")
xlabel("Time (s)")
ylabel("L2 Norm")
title("L2 Norm Comparing 200 and 400 steps")
legend(loc = 0)

figure()
plot(times3, log(norm2), label = "c = 0.5")
plot(times6, log(norm5), label = "c = 0.25")
plot(times9, log(norm8), label = "c = 0.75")
xlabel("Time (s)")
ylabel("L2 Norm")
title("L2 Norm Comparing 400 and 800 steps")
legend(loc = 0)

figure()
plot(times4, log(norm4), label = "c = 0.5")
plot(times7, log(norm6), label = "c = 0.25")
plot(times10, log(norm9), label = "c = 0.75")
xlabel("Time (s)")
ylabel("L2 Norm")
title("L2 Norm Comparing 800 and 1600 steps")
legend(loc = 0)#