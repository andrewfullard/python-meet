# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 14:55:09 2017

@author: Andrew
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 15:17:48 2017

@author: Andrew
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Oct 09 10:30:00 2017

@author: Andrew
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Column
import seaborn as sns

def calculate_P(Q,U,name):
    return Column(np.sqrt(Q**2 + U**2), name=name)

def phase_wrap(a):
    return np.concatenate((a-1,a,a+1))

def wrap(a):
    return np.concatenate((a,a,a))

def pa_angle(x, a, b):
    return a*x + b
    
def calculate_PA(Q, U, name="PA"):
    return Column(np.rad2deg(0.5*np.arctan2(U, Q)), name=name)   
    
def rotate_pa(VPA, VP, VE, LineP, LinePA, adjust=0):
    mean_pa = np.average(VPA + adjust, weights = 1./np.rad2deg(VE/2*VP))
    
    print "Mean PA: ", mean_pa
    q_rot_v = VP*np.cos(2*np.deg2rad(VPA - mean_pa))
    u_rot_v = VP*np.sin(2*np.deg2rad(VPA - mean_pa))
    
    q_rot_line = LineP*np.cos(2*np.deg2rad(LinePA - mean_pa))
    u_rot_line = LineP*np.sin(2*np.deg2rad(LinePA - mean_pa))
    
    return q_rot_v, u_rot_v, q_rot_line
    
def plotting_BVR(phase, B, BE, V, VE, R, RE, new):
    fig, ax1 = plt.subplots(figsize=(12,6))

    phase_wrapped_v = phase_wrap(phase)
    pol_wrapped_b = wrap(B)
    yerr_wrapped_b = wrap(BE)
    pol_wrapped_v = wrap(V)
    yerr_wrapped_v = wrap(VE)
    pol_wrapped_r = wrap(R)
    yerr_wrapped_r = wrap(RE)
    
    
    ax1.errorbar(phase_wrapped_v, pol_wrapped_b, yerr = yerr_wrapped_b, fmt = "b.-", label="B")
    if new:
        ax1.errorbar(phase_wrapped_v[new], pol_wrapped_b[new], yerr = yerr_wrapped_b[new], fmt = "bo")
    
    ax1.errorbar(phase_wrapped_v, pol_wrapped_v, yerr = yerr_wrapped_v, fmt = "k.-", label="V")
    if new:
        ax1.errorbar(phase_wrapped_v[new], pol_wrapped_v[new], yerr = yerr_wrapped_v[new], fmt = "ko")
    
    ax1.errorbar(phase_wrapped_v, pol_wrapped_r, yerr = yerr_wrapped_r, fmt = "r.-", label="R")
    if new:
        ax1.errorbar(phase_wrapped_v[new], pol_wrapped_r[new], yerr = yerr_wrapped_r[new], fmt = "ro")
    
    ax1.set_xlabel("Phase")
    ax1.set_ylabel("%q")
    plt.ylim((int(np.min(V)), int(np.max(V)+1)))
    ax1.set_xlim((-0.2, 1.2))
    ax1.vlines([0,1], np.min(V)-2, np.max(V)+2, colors="gray", linestyles = "dashed")
    ax1.legend(loc = "right center")
    ax1.tick_params(direction="in", which='both')
    ax1.legend()
    ax1.minorticks_on()
    #fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    
def plotting_lines(phase, V, VE, line1, line1err, line1label="Line 1", line2=[0], line2err=[0], line2label="Line 2", highlight_phases=[]):
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(3.2, 3), dpi=300)
    
    ax1.plot(wr47_light["phase"], wr47_light["c2"],"k.")
    ax1.plot(wr47_light["phase"]-1, wr47_light["c2"],"k.")
    ax1.plot(wr47_light["phase"]+1, wr47_light["c2"],"k.")
    ax1.set_xlim((-0.2, 1.2))
    ax1.invert_yaxis()
    ax1.set_ylabel("Differential magnitude")
    ax1.set_ylim((2.41, 2.2))
    ax1.vlines((0, 1), 2.2, 2.41, colors="gray", linewidth= 0.8, linestyles = 'dashed')
    ax1.tick_params(direction="in", which='both')
    ax1.minorticks_on()
    
    phase_wrapped = phase_wrap(phase)
    pol_wrapped_v = wrap(V)
    pol_wrapped_line1 = wrap(line1)
    pol_wrapped_line2 = wrap(line2)
    
    #plt.figure(figsize=(8,  6))
    #plt.title("WR 47 Line Polarization",y=1.08)
    ax2.errorbar(phase_wrapped, pol_wrapped_line1, yerr = 0, fmt = "bo--", markersize=4, label = line1label, linewidth= 0.8)
    ax2.errorbar(phase_wrapped, pol_wrapped_line1, yerr = 0, fmt = "bo",  markersize=6, markevery = np.in1d(phase_wrapped, highlight_phases))
    ax2.errorbar(phase_wrapped, pol_wrapped_line2, yerr = 0, fmt = "ro--", markersize=4, label = line2label, linewidth= 0.8)
    ax2.errorbar(phase_wrapped, pol_wrapped_line2, yerr = 0, fmt = "ro", markersize=6, markevery = np.in1d(phase_wrapped, highlight_phases))
    #ax2.errorbar(phase_wrapped, pol_wrapped_line3, yerr = 0, fmt = "o-", label = line3label)
    ax2.errorbar(phase_wrapped, pol_wrapped_v, yerr = 0, fmt = "ko-", markersize=4, label = "V-Band")
    ax2.errorbar(phase_wrapped, pol_wrapped_v, yerr = 0, fmt = "ko", markersize=6, markevery = np.in1d(phase_wrapped, highlight_phases))
    
    ax2.errorbar(0.02, 3.5, yerr = np.mean(line1err), fmt = "bo--", markersize=4)
    ax2.errorbar(0.04, 3.5, yerr = np.mean(line2err), fmt = "ro--", markersize=4)
    #ax2.errorbar(0.17, 3.4, yerr = np.mean(line3err), fmt = "o-")
    ax2.errorbar(0.06, 3.5, yerr = np.mean(VE), fmt = "ko-", markersize=4)
    ax2.text(0.08, 3.32, "Representative\nerror")#, fontsize=16)
    
    ax2.set_xlabel("Phase")
    ax2.set_ylabel("Rotated q (%)")
    plt.ylim((int(np.min(V)), int(np.max(V)+1)+0.1))
    ax2.set_xlim((-0.2, 1.2))
    ax2.vlines((0,1), int(np.min(V)), int(np.max(V)+1), linewidth= 0.8, colors="gray", linestyles = "dashed")
    ax2.legend(loc = "lower right", frameon=True)
    ax2.tick_params(direction="in", which='both')
    ax2.minorticks_on()
    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    #plt.tight_layout(pad=0.1,h_pad=0)
    
    plt.savefig('wr47_final.png', dpi='figure')


path="WR_Star_Plots/"

#wr30_light = ascii.read(path+"wr30_light.csv")
wr47_light = ascii.read(path+"wr47_light.csv", data_start=1)

#wr30_light.sort('phase')
wr47_light.sort('phase')

wr47_4686 = ascii.read(path+"wr47-4686-line.txt", 
                  names = ["Date", "Phase", "Q", "U", "Err", "P", "PA"])

wr47_6563 = ascii.read(path+"wr47-6563-line.txt", 
                  names = ["Date", "Phase", "Q", "U", "Err", "P", "PA"])

wr47 = ascii.read(path+"filterpol-phase_wr47.txt",data_start=3, 
                  names = ["Date", "Phase", "BQ", "BU", "BE", "VQ", "VU", "VE",	"RQ", "RU", "RE"])

BP = Column(np.sqrt(wr47["BQ"]**2 + wr47["BU"]**2), name="BP")

VP = Column(np.sqrt(wr47["VQ"]**2 + wr47["VU"]**2), name="VP")

RP = Column(np.sqrt(wr47["RQ"]**2 + wr47["RU"]**2), name="RP")

wr47.add_column(BP, index = 8)
wr47.add_column(VP, index = 9)
wr47.add_column(RP, index = 10)

wr12_4686 = ascii.read(path+"wr12_4686.txt", 
                  names = ["Date", "Phase", "Q", "U", "Err", "P", "PA"])

wr12_6563 = ascii.read(path+"wr12_6563.txt", 
                  names = ["Date", "Phase", "Q", "U", "Err", "P", "PA"])

wr12 = ascii.read(path+"filterpol-phase_wr12.txt",data_start=3, 
                  names = ["Date", "Phase", "BQ", "BU", "BE", "VQ", "VU", "VE",	"RQ", "RU", "RE"])

BP = calculate_P(wr12["BQ"], wr12["BU"], "BP")

VP = calculate_P(wr12["VQ"], wr12["VU"], "VP")

RP = calculate_P(wr12["RQ"], wr12["RU"], "RP")

wr12.add_column(BP, index = 8)
wr12.add_column(VP, index = 9)
wr12.add_column(RP, index = 10)

sns.set()
sns.set_style("white")
sns.set_context("paper", font_scale=0.8)
sns.set_style("ticks")

wr47_4686.sort('Phase')
wr47_6563.sort('Phase')
wr47.sort('Phase')

wr12.sort('Phase')
wr12_4686.sort('Phase')
wr12_6563.sort('Phase')
#Line polarization

wr47_VPA = calculate_PA(wr47["VQ"], wr47["VU"], name="VPA")
wr12_VPA = calculate_PA(wr12["VQ"], wr12["VU"], name="VPA")


wr47_q_rot_v, wr47_u_rot_v, wr47_q_rot_4686 = rotate_pa(wr47_VPA, wr47["VP"],  wr47["VE"], wr47_4686["P"], wr47_4686["PA"])
wr47_q_rot_v, wr47_u_rot_v, wr47_q_rot_6563 = rotate_pa(wr47_VPA, wr47["VP"],  wr47["VE"], wr47_6563["P"], wr47_6563["PA"])

wr12_q_rot_v, wr12_u_rot_v, wr12_q_rot_4686 = rotate_pa(wr12_VPA, wr12["VP"],  wr12["VE"], wr12_4686["P"], wr12_4686["PA"])
wr12_q_rot_v, wr12_u_rot_v, wr12_q_rot_6563 = rotate_pa(wr12_VPA, wr12["VP"],  wr12["VE"], wr12_6563["P"], wr12_6563["PA"])

print calculate_PA(wr12_q_rot_v, wr12_u_rot_v)

#custom_palette = sns.diverging_palette(220, 20, n=2, l=60, s=95, center = "light")

#sns.set_palette(custom_palette)
#palette = sns.husl_palette(l=1,s=1)
#sns.set_palette(palette)

plotting_BVR(wr47["Phase"], wr47["BP"], wr47["BE"], wr47["VP"], wr47["VE"], wr47["RP"], wr47["RE"], False)

plotting_BVR(wr12["Phase"], wr12["BP"], wr12["BE"], wr12["VP"], wr12["VE"], wr12["RP"], wr12["RE"], False)

#plotting_lines(wr47["Phase"], wr47["VP"], wr47["VE"], wr47_4686["P"], wr47_4686["Err"], line2=wr47_6563["P"], line2err=wr47_6563["Err"])

plotting_lines(wr47["Phase"], wr47_q_rot_v, wr47["VE"], wr47_q_rot_4686, wr47_4686["Err"], line1label = "He II 4686", line2=wr47_q_rot_6563, line2err=wr47_6563["Err"], line2label = "HE II/H$\\alpha$ 6563", highlight_phases = [0.996-1, 0.314, 0.674, 0.996])

#plotting_lines(wr12["Phase"], wr12_q_rot_v, wr12["VE"], wr12_q_rot_4686, wr12_4686["Err"], line1label = "He II 4686", line2=wr12_q_rot_6563, line2err=wr12_6563["Err"], line2label = "HE II/H$\\alpha$ 6563")