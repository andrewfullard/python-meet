# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 12:55:56 2016

@author: WolfeTM
"""

from matplotlib import pyplot as plt
import numpy as np

### --------------------------------------------------------------------------
# plt.annotate() - useful for annotating plots!

x = np.linspace(1.,10.,10) # Sample space of x values
y = 3.*x-1 # Make y a linear function because why not?
plt.figure(1)
plt.title('Example use of plt.annotate()',fontsize=16) # Plot title command
plt.xlabel('x axis label',fontsize=16) # Plot x label command
plt.ylabel('y axis label',fontsize=16) # Plot y label command
plt.plot(x,y,'ko') # Plot data as black circles
plt.annotate('check out that \n data point',(x[4],y[4]),xycoords='data',\
             color='red',xytext=(x[4]+2,y[4]+2),\
             arrowprops=dict(color='blue',arrowstyle='->'),zorder=1,\
             horizontalalignment='left',verticalalignment='bottom',fontsize=16)
             
# You can also iterate over a data set to program multiple annotations: let's
# define a list of strings to annotate each data point with.
names = ['check out this \n data point', 'what about this one?!',\
         'this one is \n cool too!', 'this one is lame',\
         'check out that \n data point', 'don\'t forget about \n this one',
         'this is the \n best data point', 'this is the worst',\
         'this one broke physics', 'this is the last one \n I promise']

# We can even iterate over the colors!
colors = ['red', 'green', 'blue', 'black', 'cyan', 'black', 'blue', 'yellow',\
          'orange', 'red']

# Make the plot:
plt.figure(2)
plt.title('Example use of iterating plt.annotate()',fontsize=16) # Plot title command
plt.xlabel('x axis label',fontsize=16) # Plot x label command
plt.ylabel('y axis label',fontsize=16) # Plot y label command
plt.plot(x,y,'ko') # Plot data as black circles

# Now, iterate:
for i,txt in enumerate(names):
    plt.annotate(names[i],(x[i],y[i]),xycoords='data',\
             color=colors[i],xytext=(x[i]+1,y[i]-1),\
             arrowprops=dict(color=colors[i],arrowstyle='->'),zorder=1,\
             horizontalalignment='left',verticalalignment='bottom',fontsize=16)
# leave this line blank for some reason
         
         
         
###---------------------------------------------------------------------------
# Example use of making legends in your plots
x = np.linspace(1.,10.,10) # Sample space of x values
y = 100.*x-1 
y2 = x**2
y3 = x**3
plt.figure(3)
plt.title('Example use of plt.legend()',fontsize=16) # Plot title command
plt.xlabel('x axis label',fontsize=16) # Plot x label command
plt.ylabel('y axis label',fontsize=16) # Plot y label command
# The most basic way to create a legend is to use the label option in plt.plot()
# The labels will show up once the plt.legend() command is called.
plt.plot(x,y,'ko',label='Linear!')
plt.plot(x,y2,'b-',label='Squared!')
plt.plot(x,y3,'r--',label='Cubed!')
plt.legend()
# plt.legend can be used much more flexibly; consult matplotlib's documentation
# for more details



###---------------------------------------------------------------------------
# Example use of plt.quiver for vector field plots, taken from 
# http://matplotlib.org/examples/pylab_examples/quiver_demo.html

# Create the field data from X & Y values.
X, Y = np.meshgrid(np.arange(0, 2 * np.pi, .2), np.arange(0, 2 * np.pi, .2))
U = np.cos(X)
V = np.sin(Y)

plt.figure(4)

Q = plt.quiver(U, V)
# plt.quiverkey() can be used to help show arrow properties.
qk = plt.quiverkey(Q, 0.5, 0.92, 2, r'$2 \frac{m}{s}$', labelpos='W',
                   fontproperties={'weight': 'bold'})
                   
# Configure axis tick settings:
l, r, b, t = plt.axis()
dx, dy = r - l, t - b
plt.axis([l - 0.05*dx, r + 0.05*dx, b - 0.05*dy, t + 0.05*dy])

plt.title('Minimal arguments, no kwargs')

###----------------------------------------------------------------------------
# Putting it all together: what if you need multiple layers of information,
# using regular image and data  plotting techniques along with annotations and
# arrows from plt.quiver?

# Below is an example plot I made for a paper that'll never be published. This 
# involved over-plotting the positions of stars onto an image of part of the sky
# in infrared. Along with annotating the positions, I needed to represent the 
# apparent motions of the stars with plt.quiver, to demonstrate if the stars 
# were moving in similar directions.

# The image itself is in black and white, but enhanced by color-contours in
# order to demonstrate how complicated the dust structure in that region of
# space is.

# Because of all the info needed, the zorder option in the plot commands is
# used to specify how to order each type of plot. That way, the contours aren't
# overlayed ontop of the annotations, for example.

# Number designations of the stars whos positions will be identified.
numbers = np.array(['1', '2', '3', '4', '6', '5', '7', '8', '9', '10', '11'])

# Pixel coordinates of each star in the image.
xpix = np.array([348, 354, 351, 368, 424, 321, 277, 444, 274, 285, 289])
ypix = np.array([401, 392, 386, 382, 395, 450, 333, 385, 295, 287, 273])
# A coordinate transform will be used for clipping the image to a specific
# region.
xpix_sample = xpix-250
ypix_sample = ypix-250

# Proper motion data for each star, to be represented by the use of plt.quiver()
pmD = np.array([-2.66, -4.80, -4.17, -3.88, -6.00, -2.81, -1.31, -1.33, -8.22, -3.13, -4.70])
pmR = np.array([-0.86, -2.00, -3.00, 0.26, 0.50, 1.86, -1.34, 0.65, -0.23, 0.43, 1.60])

# Use astropy.io.fits to load a fits-format image file
from astropy.io import fits

# Open the image file
hdulist = fits.open('i343b3h0.fit')
# Define an object containing the image data
iras60 = hdulist[0].data[0]
# Get rid of negative values because they're useless and can screw up the
# color scaling of the image and the contours we'll use.
iras60[np.where(iras60 < 0)] = 0.
# Flatten the scale of the image data to enhance smaller details.
iras60_scaled = np.sqrt(iras60)
# Crop the image to the relevant region
iras60_scaled_sample = iras60_scaled[250:,250:]

# Create a meshgrid to be used with the plots
x = np.arange(250)
y = np.arange(250)
X,Y = np.meshgrid(x,y)

# Create figure, use imshow to make the first layer, representing the image
# itself.
plt.figure(5)
plt.title('Complicated star field with \n multiple layers of different info',\
           fontsize=16)
fig = plt.imshow(iras60_scaled_sample, cmap=plt.cm.get_cmap('bone'),\
                 origin='lower',zorder=0)

# Define levels for the contours to be plotted (I defined this manually for
# better control over the contours)
levels = np.linspace(0.,2.,20)
CS = plt.contour(X,Y,iras60_scaled_sample,levels,linewidths=(0.75,),\
                 cmap=plt.cm.get_cmap('Set1'),zorder=1)

# Hide x and y axis ticks, values, and labels (not needed in this case)
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)

# Use plt.colorbar() to specifically demonstrate the conour level values.
plt.colorbar(CS, shrink=0.8, extend='both',label='relative intensity units')

# Truncate the plot to the relevant regions
plt.xlim((15,215))
plt.ylim((10,220))

# Iterate over each star number for annotation. Star # 2 is skipped because
# if plotted with the same iteration algorithm, the text will overlap with 
# Star # 1; so we annotate it separately.
for i,txt in enumerate(numbers):
	if i != 2:
		plt.annotate(txt,(xpix_sample[i],ypix_sample[i]),xycoords='data',\
                    color='white',xytext=(xpix_sample[i]+10,ypix_sample[i]+10),\
                    arrowprops=dict(color='white',arrowstyle='->'),zorder=2,\
                    horizontalalignment='left',verticalalignment='bottom',\
                    fontsize=16)
# Annotate Star 2 separately
plt.annotate(numbers[2],(xpix_sample[2],ypix_sample[2]),xycoords='data',\
             color='white',xytext=(xpix_sample[2]+15,ypix_sample[2]+7),\
             arrowprops=dict(color='white',arrowstyle='->'),zorder=3,\
             horizontalalignment='left',verticalalignment='bottom',fontsize=16)
# Create quivers to represent proper motion data.
QVM = plt.quiver(xpix_sample,ypix_sample,pmR,pmD,angles='uv',scale=50,\
                 color='w',linewidth=0.5,headwidth=1,zorder=4)
