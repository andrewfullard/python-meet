# -*- coding: utf-8 -*-
"""
Created on Wed Oct 02 09:04:01 2013

@author: c0918140
"""
import numpy

inFunc = raw_input("Please enter a function to integrate where x is the variable\n")

upperBoundInt = raw_input("Please enter an upper bound for the integration\n")

lowerBoundInt = raw_input("Please enter a lower bound for the integration\n")

lowerBound = float(lowerBoundInt)

upperBound = float(upperBoundInt)

lowFunc = inFunc.lower()

def f(x):
    y = eval(lowFunc)
    return y
  
def gaussLegendre3(a, b):
    mod = (b-a)/2
    mod2 = (a+b)/2
    R1 = mod*-1*(numpy.sqrt(15.0)/5.0) + mod2
    R2 = mod2
    R3 = mod*(numpy.sqrt(15.0)/5.0) + mod2
    result = mod*((5.0/9.0)*f(R1) + (8.0/9.0)*f(R2) + (5.0/9.0)*f(R3))
    return result
    
def gaussLegendre4(a, b):
    mod = (b-a)/2
    mod2 = (a+b)/2
    result = mod*((0.65214515)*f(mod*-0.33998104+mod2) + (0.65214515)*f(mod*0.33998104+mod2) + (0.34785485)*f(mod*-0.86113631+mod2) + (0.34785485)*f(mod*0.86113631+mod2))
    return result
    
    
print "Result for n = 3 = "
print gaussLegendre3(lowerBound, upperBound)
print "Result for n = 4 = "
print gaussLegendre4(lowerBound, upperBound)
