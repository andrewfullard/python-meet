# -*- coding: utf-8 -*-
"""
Created on Wed Oct 02 09:04:01 2013

@author: c0918140
"""

#Numerical integration methods

#User function input using python mathematical format e.g. powers are ** not ^
inFunc = raw_input("Please enter a function to integrate where x is the variable\n")

#Make all variables lower case
lowFunc = inFunc.lower()

#User input for upper and lower bounds
upperBoundInt = raw_input("Please enter an upper bound for the integration\n")

#User input for lower bound
lowerBoundInt = raw_input("Please enter a lower bound for the integration\n")

#User input for number of integration points
integrationPoints = input("Please enter the number of integration points as an integer greater than 1\n")

#Float conversion for bounds
lowerBound = float(lowerBoundInt)

upperBound = float(upperBoundInt)

#Define function f(x) using eval on the user input function
def f(x):
    y = eval(lowFunc)
    return y

#Old derivation tests
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

#Numerical derivative of the function f(x)
def fderivnumerical(u, l, w):
	#Python list
    result = []
    while l < u:
        result.append(((f(l+w)+f(l-w)-2*f(l))/(w*w)))
        l += w
    return result
 
#Width of each trapezoidal bin 
width = float(((upperBound - lowerBound)/integrationPoints))

#Trapezoidal integration function
def trapezoid(u, l, i, h):
    N = 1
    summed = 0
    z = l + h
	#iterate over trapezoids, summing areas, ignoring first and last
    while N < i:
        summed += (f(z)*h)
        z += h
        N += 1     
       
	#Final result with first and last trapezoids included
    result = (h/2)*f(l) + (h/2)*f(u) + summed
    return result

#Convergence testing: does the result converge, using up to 4x as many points
def selfConvergenceTest(u, l, i, w):
    result = (trapezoid(u, l, i, w) - trapezoid(u, l, i, (w/2)))/(trapezoid(u, l, i, (w/2))-trapezoid(u, l, i, (w/4)))
    return result    
    
#Old error function
#def error(u, l, w):
#    result = ((w*w)/12)*(u-l)*max(fderivlist(u, l, w))
#    return result
    
#numerical error on result
def errorNumerical(u, l, w):
    result = ((w*w)/12)*(u-l)*max(fderivnumerical(u, l, w))
    return result
 
#Old residual function 
#def residual(u, l, i, w):
#    r = 40000 - trapezoid(u, l, i, w)
#    R = log10(r)
#    return R
    
#print final results to screen
print "Result = "
print trapezoid(upperBound, lowerBound, integrationPoints, width)
print "Convergence test = "
print selfConvergenceTest(upperBound, lowerBound, integrationPoints, width)
#print "Error = "
#print error(upperBound, lowerBound, width)
print "Error (numerical derivative) = "
print errorNumerical(upperBound, lowerBound, width)
