# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 09:47:03 2012

@author: marko
"""

from pylab import *
        
def Advect(y0):
    y1 = zeros(len(y0))
    for i in range(1,len(y0)):
        y1[i] = y0[i-1]
    return y1
        
def ApplyBoundaryConditions(y):
    Np = len(y)-2
    y[0] = y[Np]
    y[Np+1] = y[1]
    return y
        
def ErrorNorm(xpts,v,t,dat):
    norm = zeros(len(dat))
    for i in range(0,len(dat)):
        soln = exp(-((mod(xpts-v*t[i],1)-0.5)**2)/(2.*0.1**2))
        # Don't include the boundary points in the norm.
        norm[i] = sqrt(sum((dat[i,1:-1] - soln[1:-1])**2)/len(soln[1:-1]))
    return norm
    
def CvgNorm(times1,dat1,dat2):
    norm = zeros(len(times1))
    for i in range(0,len(times1)):
       norm[i] = sqrt(sum((dat1[i,1:-1] - dat2[2*i,2:-1:2])**2)/len(norm))
    return norm
    
def SolveAdvectionEquation(a, b, t0, t1, Nx, Nt, v, method):
    dx = (float) (b-a)/Nx
    dt = dx/v
    if (method == 'advect' and abs(Nt*dt - (t1-t0))> dt/2.0):
        print('Final time is not multiple of dt. Try Tmax = ',Nt*dt)
        return 0
    # We have Nx+1 points for Nx intervals, plus an extra ghost point.
    pts = linspace(a-dx,b,Nx+2)
    times = linspace(t0,t1,Nt+1)
    solndata = zeros((Nt+1,Nx+2))
    sigma2 = 2.*0.1**2
    solndata[0] = exp(-((pts-0.5)**2)/(sigma2))
    solndata[0] = ApplyBoundaryConditions(solndata[0])
    for tn in range(0,Nt):
        if (method == 'advect'):
            solndata[tn+1] = Advect(solndata[tn])
            solndata[tn+1] = ApplyBoundaryConditions(solndata[tn+1])
    return times, pts[1:], solndata[:,1:]

times1, grid1, soln1 = SolveAdvectionEquation(0, 1, 0, 5.0, 10, 100, 2.0, 'advect')
times2, grid2, soln2 = SolveAdvectionEquation(0, 1, 0, 5.0, 100, 1000, 2.0, 'advect')
times3, grid3, soln3 = SolveAdvectionEquation(0, 1, 0, 5.0, 1000, 10000, 2.0, 'advect')

norm1 = ErrorNorm(grid1,2.0,times1,soln1)
norm2 = ErrorNorm(grid2,2.0,times2,soln2)
norm3 = ErrorNorm(grid3,2.0,times3,soln3)    
    
figure(1)
clf()
plot(times1,norm1,'r-')
plot(times2,norm2,'g-')
plot(times3,norm3,'b-')