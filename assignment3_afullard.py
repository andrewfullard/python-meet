# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 13:58:21 2013

@author: Andrew
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

alpha = 500
sigma = 100

# Step 1
# Finite difference vectors

def finiteDifference(N, dx):
    A = 1./(np.ones(N + 2)*dx*dx)
    B = -2./(np.ones(N + 2)*dx*dx)
    B[0] = -1./(dx*dx)
    B[N+1] = -3./(dx*dx)
    C = 1./(np.ones(N + 2)*dx*dx)
    return A, B, C
    
# Step 2
# Matrix vector product
    
def matrixVectorProduct(a, b, c, u):
    
    v = np.zeros(len(u))
    v[0] = b[0]*u[0] + c[0]*u[1]
    
    for i in range(1, len(u) - 2):
        v[i] = a[i]*u[i-1] + b[i]*u[i] + c[i]*u[i+1]
        
    v[len(u) - 1] = a[len(u)-1]*u[len(u)-2] + b[len(u)-1]*u[len(u)-1]
    
    return v

# Step 3
# LHS of the elliptic equation  

def ellipticLHS(a, b, c, x, u, N):
    
    D2u = matrixVectorProduct(a, b, c, u)
    
    F = D2u + alpha*np.exp(-sigma*x**2) - u**3
    
    return F    

# Step 4
# Derivative vectors of the LHS side of the elliptic equation

def derivativeLHS(dx, u):
    
    a = 1./(np.ones(len(u))*dx*dx)
    c = a
    
    b = (-2./(np.ones(len(u))*dx*dx)) - (3*u*u)

    b[0] = (-1./(dx*dx)) - (3*u[0]*u[0])
    b[len(u) - 1] = (-1./(dx*dx)) - (3*u[len(u)-1]*u[len(u)-1])
    
    return a, b, c
    
# Step 5
# Standard tridiagonal method
    
def tridiagonalMethod(A, B, C, D, x):
    
    Cprime = np.zeros(len(A))
    Dprime = np.zeros(len(A))
    
    u = np.zeros(len(A))

    Cprime[1] = C[1]/B[1]
    Dprime[1] = D[1]/B[1]    
    
    for i in range(2, len(A) - 2):
        Cprime[i] = C[i]/(B[i] - (Cprime[i - 1]*A[i]))
        Dprime[i] = (D[i] - Dprime[i - 1]*A[i])/(B[i] - Cprime[i - 1]*A[i])
    
    for i in range(len(A) - 2, 0, -1):
        u[i] = Dprime[i] - Cprime[i]*u[i+1]
    
    return u
   
# Step 6
# Newton-Raphson method incorporating all previous steps
   
def newtonRaphson(lbound, rbound, N, E):
    dx = (rbound - lbound)/N
    x = np.linspace(lbound - dx, rbound + dx, N + 2)
    uguess = np.ones(N + 2)
    delta = np.zeros(20)
    for i in range(0, 19):
        a, b, c = finiteDifference(N, dx)
        
        F = ellipticLHS(a, b, c, x, uguess, N)
        
        A, B, C = derivativeLHS(dx, uguess)
        
        d = tridiagonalMethod(A, B, C, F, x)   
        
        uguess = uguess - d
        squaresum = 0
        for j in range(0, len(d)):
            squaresum += d[j]*d[j]
        delta[i] = np.sqrt(squaresum)
        if delta[i] < E:
            break
        
    return x[1:N+1], uguess[1:N+1], delta

# Step 7
# Running the code over the required range and accuracy            
            
x1, u1, delta1 = newtonRaphson(-1., 1., 40, 1e-10)
x2, u2, delta2 = newtonRaphson(-1., 1., 80, 1e-10)
x3, u3, delta3 = newtonRaphson(-1., 1., 160, 1e-10)

# Step 8
# Higher accuracy

x4, u4, delta4 = newtonRaphson(-1., 1., 160, 1e-20)

# Interpolation for convergence testing

u3i = interpolate.interp1d(x2,u2,'cubic') 
u3c = u3i(x1)

u4i = interpolate.interp1d(x3,u3,'cubic') 
u4c = u4i(x1)

# Convergence calculation

conv1 = abs(u3c - u1)
conv2 = abs(u4c - u3c)

# Plotting

plt.figure("Convergence")
plt.plot(x1, conv1, label = "40/80 points")
plt.plot(x1, conv2, label = "80/160 points")
plt.legend(loc="upper right")
plt.xlabel("x")
plt.ylabel("Convergence")

plt.figure("Solutions")
plt.plot(x1, u1, label = "40 points")
plt.plot(x2, u2, label = "80 points")
plt.plot(x3, u3, label ="160 points")
plt.legend(loc="upper right")
plt.xlabel("x")
plt.ylabel("u(x)")

plt.figure()
plt.semilogy(delta1, label = "40 points")
plt.semilogy(delta2, label = "80 points")
plt.semilogy(delta3, label = "160 points")
plt.legend(loc="upper right")
plt.xlabel("Iteration")
plt.ylabel("delta")

plt.figure("10^-20 tolerance")
plt.semilogy(delta4)
plt.xlabel("Iteration")
plt.ylabel("delta")
