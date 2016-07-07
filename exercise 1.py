# -*- coding: utf-8 -*-
"""
Created on Wed Oct 02 09:04:01 2013

@author: c0918140
"""

inFunc = raw_input("Please enter a function to integrate where x is the variable\n")

lowFunc = inFunc.lower()

upperBoundInt = raw_input("Please enter an upper bound for the integration\n")

lowerBoundInt = raw_input("Please enter a lower bound for the integration\n")

integrationPoints = input("Please enter the number of integration points as an integer greater than 1\n")

lowerBound = float(lowerBoundInt)

upperBound = float(upperBoundInt)

def f(x):
    y = eval(lowFunc)
    return y

#def fderiv(x):
#    y = 6*x
#    return y
#
#def fderivlist(u, l, w):
#    result = []
#    while l < u:
#        result.append(fderiv(l))
#        l += w
#    return result

def fderivnumerical(u, l, w):
    result = []
    while l < u:
        result.append(((f(l+w)+f(l-w)-2*f(l))/(w*w)))
        l += w
    return result
    
width = float(((upperBound - lowerBound)/integrationPoints))

def trapezoid(u, l, i, h):
    N = 1
    summed = 0
    z = l + h
    while N < i:
        summed += (f(z)*h)
        z += h
        N += 1     
        
    result = (h/2)*f(l) + (h/2)*f(u) + summed
    return result

def selfConvergenceTest(u, l, i, w):
    result = (trapezoid(u, l, i, w) - trapezoid(u, l, i, (w/2)))/(trapezoid(u, l, i, (w/2))-trapezoid(u, l, i, (w/4)))
    return result    
    
#def error(u, l, w):
#    result = ((w*w)/12)*(u-l)*max(fderivlist(u, l, w))
#    return result
    
def errorNumerical(u, l, w):
    result = ((w*w)/12)*(u-l)*max(fderivnumerical(u, l, w))
    return result
    
#def residual(u, l, i, w):
#    r = 40000 - trapezoid(u, l, i, w)
#    R = log10(r)
#    return R
    
print "Result = "
print trapezoid(upperBound, lowerBound, integrationPoints, width)
print "Convergence test = "
print selfConvergenceTest(upperBound, lowerBound, integrationPoints, width)
#print "Error = "
#print error(upperBound, lowerBound, width)
print "Error (numerical derivative) = "
print errorNumerical(upperBound, lowerBound, width)