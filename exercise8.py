# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 09:10:22 2013

@author: c0918140
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

def f(x):
    y = -np.exp(x)*(-2. + 2.*x + 5.*x**2. + x**3.)
    return y

def R(r):
    y = -r**2./(1 + r**6)
    return y
    
def solver(N):
    dx = 1./(N + 2)
    xspace = np.linspace(0. - dx, 1. + dx, N + 2)
    x = np.linspace(0., 1., N)
    A = 1./(np.ones(N + 2)*dx*dx)
    B = -2./(np.ones(N + 2)*dx*dx)
    B[0] = -1./(dx*dx)
    B[N+1] = -3./(dx*dx)
    C = 1./(np.ones(N + 2)*dx*dx)
    Cprime = np.zeros(N + 2)
    D = f(xspace)
    Dprime = np.zeros(N + 2)
    
    u = np.zeros(N + 2)

    Cprime[1] = C[1]/B[1]
    Dprime[1] = D[1]/B[1]    
    
    i = 2
    
    while i < N:
        Cprime[i] = C[i]/(B[i] - (Cprime[i - 1]*A[i]))
        Dprime[i] = (D[i] - Dprime[i - 1]*A[i])/(B[i] - Cprime[i - 1]*A[i])
        i += 1

    i = N
    
    while i > 0:
        u[i] = Dprime[i] - Cprime[i]*u[i+1]
        i -= 1
    
    return x, u[1:N+1]
    
def solver_elliptic(N):
    dx = 1./(N + 2)
    xspace = np.linspace(0. - dx, 1. + dx, N + 2)
    x = np.linspace(0. + dx, 1. - dx, N)
    A = 1./(np.ones(N + 2)*dx*dx)
    B = -2./(np.ones(N + 2)*dx*dx)
    B[0] = -1./(dx*dx)
    B[N+1] = -3./(dx*dx)
    C = 1./(np.ones(N + 2)*dx*dx)
    Cprime = np.zeros(N + 2)
    D = R(xspace)
    Dprime = np.zeros(N + 2)
    
    u = np.zeros(N + 2)

    Cprime[1] = C[1]/B[1]
    Dprime[1] = D[1]/B[1]    
    
    i = 2
    
    while i < N:
        Cprime[i] = C[i]/(B[i] - (Cprime[i - 1]*A[i]))
        Dprime[i] = (D[i] - Dprime[i - 1]*A[i])/(B[i] - Cprime[i - 1]*A[i])
        i += 1

    i = N
    
    while i > 0:
        u[i] = Dprime[i] - Cprime[i]*u[i+1]
        i -= 1
    
    return x, u[1:N+1]

x1, u1 = solver(100.)
x2, u2 = solver(200.)

u3 = interpolate.interp1d(x2,u2,'cubic') 
u3c = u3(x1)

conv = u3c - u1
       
plt.figure("Convergence")
plt.plot(x1, conv)

plt.figure()
plt.plot(x1, u1, label="100 points")
plt.plot(x2, u2, label="200 points")
plt.title("Convergence test")
plt.legend(loc = 0)
