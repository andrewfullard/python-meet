# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 17:11:53 2016

@author: rachelbennet
"""

# The purpose of this code is to demonstrate the import and plotting of basic
# ASCII tables using Astropy, Matplotlib, and LaTeX formatting.

#---------------------------------------------------------------------------

# The code below contains important packages to import:

from matplotlib import pyplot as plt # Matplotlib is the most commonly used
    # plotting package for my purposes and functions similarly to MATLAB
    # plotting (hence "MATplotlib").  Since Matplotlib is a very large 
    # package, we will import only the subpackage we need: Pyplot.  We
    # import it as plt to reduce our typing load.

import os # OS allows us to do fun things like changing our working directory.

from astropy.io import ascii # Astropy is a user-created library that contains
    # (most interesting for our purposes) the powerful ASCII package.  This
    # package allows us to import ASCII tables, such as tab-delineated .txt
    # files, as Python arrays.

#---------------------------------------------------------------------------

# WARNING: LaTeX APPROACHING!  YOU MUST HAVE A LaTeX COMPILER INSTALLED ON
# YOUR COMPUTER TO USE THIS CODE!

# On my computer (MacBook Pro (Retina, 13-inch, Late 2013), 2.6 GHz Intel
    # Core i5), 8 GB 1600 MHz DDR3, OSX Yosemite, Python 3.5.1), attempt-
    # ing to use this code inside Spyder (including Run->Run (F5), the
    # iPython console, and the Console) causes the code to terminate with
    # irreconcilable errors.  If you encounter this, run the routine from
    # the Terminal (Applications->Utilities->Terminal).

from matplotlib import rc # Again, Matplotlib is massive, so we only want
    # the RC subpackage, which contains our nice fonts and LaTeX syntax
    # for special characters.

# Below are two options for the type of font family you want to use.
    # Typically, for publication, you will use a serif font like Palatino.

# Sans-serif (Helvetica):
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

# Serif (Palatino):
rc('font',**{'family':'serif','serif':['Palatino'],'size':'14'}) # This
    # summons the font Palatino from the serif font family and designates
    # the font size as 14 pt.

rc('text', usetex=True) # This tells Python that we would like to be able
    # to type LaTeX commands to create math and other symbols in our text
    # strings.

#-------------------------------------------------------------------------

# Now we would like to import our data file into this routine.  First,
    # download the expdata.txt, lndata.txt, and labeldata.txt files for
    # this demo and place them in an easy-to-find folder.  Then copy that
    # file path (do not include the file name!) and paste it into the
    # indicated place below:

path = "INSERT YOUR FILE PATH HERE" # We are giving the text string that
    # is your file path the object designation "path."

os.chdir( path ) # The OS "change directory" command will use whatever
    # location you have specified as "path" and change your Python working
    # directory to that location.  From here on out, it will assume that
    # any files you wish to import come from this folder, and any files
    # you wish to save will save to this folder as well.

# Now we can import a data file.  Let's import the file expdata.txt from
    # the folder you specified above:

ExpDataTable = ascii.read("expdata.txt") # You must have the file extension
    # at the end and the quotation marks (to tell Python to read this as a
    # string of letters) in order for ascii.read to work properly.  This
    # command can read a wide variety of ASCII file types; see the Astropy
    # ascii.read documentation for more info.  This command imports our
    # .txt file as a Python array, which we are naming ExpDataTable for
    # clarity.

#-------------------------------------------------------------------------

# EXERCISE: Using the same procedure, import lndata.txt and labeldata.txt
    # as the objects LnDataTable and LabeledDataTable. Write your code
    # below:

#-------------------------------------------------------------------------

# In order to plot this, we'll need to separate the x- and y- values.
    # Copy the commands on lines 15, 21, 23, 61, 70, and 79 into your
    # IPython console independently (striking "Enter" after each paste),
    # then view the imported Python array by copying in the command:

ExpDataTable

# Two columns should appear in your browser.  On the left is a string of
    # inputs - x=1-10 in steps of 1.  On the right is the output of
    # exp(x) performed on these inputs.

# Python numbers its columns from left to right, 1-N, where N is the
    # total number of columns.  So if we wish to create a new array that
    # is just the x-values, we can "slice" the original array using the
    # command

ExpXValues = ExpDataTable['col1'] # This creates a new array that is only
    # the values in Column 1 (the leftmost column) of ExpDataTable and
    # gives it the object designation ExpXValues.

# Copy this command into your iPython console, then type ExpXValues to
    # view your new array.

#-------------------------------------------------------------------------

# EXERCISE: Using the same procedure, create an object named ExpYValues
    # that is the array of y-values from ExpDataTable.

#-------------------------------------------------------------------------

# If your data columns have headers (as the data in labeldata.txt does),
    # you can import them by column name like so:

LabeledXValues = LabeledDataTable['X Values']

# This is especially useful for data sets that contain more than one set of
    # x- and/or y-values, as labeldata.txt does.

#-------------------------------------------------------------------------

# EXERCISE: Import the two columns of y-values from LabeledDataTable and
    # give them the object designations LabeledExpValues and
    # LabeledLnValues.

#-------------------------------------------------------------------------

# Now we can plot this data:

plt.figure(1, figsize=(3,3)) # Designates the start of Figure 1, which
    # will be 3"x3".

plt.plot(ExpXValues, ExpYValues) # The x- and y-arrays are given to the
    # plotting routine with this command.

plt.xlabel("Time (years)") # Axis labels default to horizontal for the x-
    # axis and vertical for the y-axis.  The words must be given to
    # Python as a string (in quotation marks).

plt.ylabel("$\beta-\tau$ Emission (mrad)") # Since Python does not
    # contain any native support for special characters like Greek
    # letters, we use LaTeX to create them.  In LaTeX, Greek letters are
    # considered mathematical symbols, and so to use them in-line with
    # other text, the commands must be contained between two $ signs as
    # shown.  The command to create the character β is \beta in LaTeX,
    # and similarly the command to create τ is \tau.  To create an upper-
    # case Greek letter, capitalize the first letter of the
    # command - e.g., \omega gives ω, while \Omega gives Ω.

plt.savefig("ExpDataSamplePlot.png", dpi=300) # If you wish to save your
    # plot right away (particularly useful when you have to run the
    # routine through Mac's Terminal), use the SaveFig command.  Your
    # file name and extension must be specified as a string, and here
    # the resolution is defined as 300 dpi.

plt.show() # This will make the routine show your plot at the end, so
    # you can view it right away.

#-------------------------------------------------------------------------

# EXERCISE: How do you place more than one curve on the same plot?  Use
    # the x- and y-values from LabeledDataTable in addition to ExpData
    # to do this.