# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 13:22:09 2018

@author: Andrew
"""

import numpy as np
from scipy.interpolate import interpolate

#Get data info from user
datapath = raw_input("Enter the name of your data file: ")
skipnum = raw_input("Enter integer number of lines to skip in the header: ")

skipnum = int(skipnum)
#load data (skip initial header)
data = np.genfromtxt(datapath, skip_header=skipnum)

#choose filter
while True:
    filter_type = raw_input("Enter a filter (B, V, R): ")
    
    if filter_type == 'B':
        filter_array = np.array([[3600.,  0.000],
                                [3700.,  0.030],
                                [3800.,  0.134],
                                [3900.,  0.567],
                                [4000.,  0.920],
                                [4100.,  0.978],
                                [4200.,  1.000],
                                [4300.,  0.978],
                                [4400.,  0.935],
                                [4500.,  0.853],
                                [4600.,  0.740],
                                [4700.,  0.640],
                                [4800.,  0.536],
                                [4900.,  0.424],
                                [5000.,  0.325],
                                [5100.,  0.235],
                                [5200.,  0.150],
                                [5300.,  0.095],
                                [5400.,  0.043],
                                [5500.,  0.009],
                                [5600.,  0.000]])
        break
    
    elif filter_type == 'V':
        filter_array = np.array([[4700.,  0.000],
                                 [4800.,  0.030],
                                 [4900.,  0.163],
                                 [5000.,  0.458],
                                 [5100.,  0.780],
                                 [5200.,  0.967],
                                 [5300.,  1.000],
                                 [5400.,  0.973],
                                 [5500.,  0.898],
                                 [5600.,  0.792],
                                 [5700.,  0.684],
                                 [5800.,  0.574],
                                 [5900.,  0.461],
                                 [6000.,  0.359],
                                 [6100.,  0.270],
                                 [6200.,  0.197],
                                 [6300.,  0.135],
                                 [6400.,  0.081],
                                 [6500.,  0.045],
                                 [6600.,  0.025],
                                 [6700.,  0.017],
                                 [6800.,  0.013],
                                 [6900.,  0.009],
                                 [7000.,  0.000]])
        break
    
    elif filter_type == 'R':
        filter_array = np.array([[5500,  0.00],
                                 [5600,  0.23],
                                 [5700,  0.74],
                                 [5800,  0.91],
                                 [5900,  0.98],
                                 [6000,  1.00],
                                 [6100,  0.98],
                                 [6200,  0.96],
                                 [6300,  0.93],
                                 [6400,  0.90],
                                 [6500,  0.86],
                                 [6600,  0.81],
                                 [6700,  0.78],
                                 [6800,  0.72],
                                 [6900,  0.67],
                                 [7000,  0.61],
                                 [7100,  0.56],
                                 [7200,  0.51],
                                 [7300,  0.46],
                                 [7400,  0.40],
                                 [7500,  0.35],
                                 [8000,  0.14],
                                 [8500,  0.03],
                                 [9000,  0.00]])
        break
    
    else:
         print "Error: choose the filter by typing B V or R"

#extract data columns
lam= data[:,0]
flux = data[:,1]

#extract filter info
wave = filter_array[:,0]
weight = filter_array[:,1]

#Interpolate filter
interp = interpolate.interp1d(wave, weight, bounds_error=False, fill_value=0.0)
weightpol = interp(lam)
filterregion = np.where(weightpol > 0)

#standard columns from polsalt reduction output
cols = [('%Q', 2), ('%U', 3), ('%Qerr', 4), ('%Uerr', 5)]

#cols = [('%Q', 1), ('%U', 2), ('%err', 3)]

#repeat for each stokes value
for stokes, i in cols:
    
    #load appropriate column
    pol = data[:,i]
    
    #integrate and convolve
    top = np.trapz((flux[filterregion]*weightpol[filterregion]*pol[filterregion]), x=lam[filterregion])
    bottom = np.trapz((flux[filterregion]*weightpol[filterregion]), x=lam[filterregion])
    measure = top/bottom
    
    #error calculation for error columns
    if i > 3:
        measure = measure/np.sqrt(len(lam[filterregion]))
    
    #print results
    print stokes, measure