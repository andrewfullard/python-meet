# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 19:48:38 2016

@author: rachelbennet
"""

# The purpose of this code is to demonstrate the import of basic ASCII
    # tables using Astropy, and the export of data tables for LaTeX using
    # Astropy.

#----------------------------------------------------------------------------

# The code below contains important packages to import:

import numpy as np # Numpy includes most math functions and constants; we
    # import this as "np" to reduce our typing load.

from matplotlib import pyplot as plt # Matplotlib is the most commonly used
    # plotting package for my purposes and functions similarly to MATLAB
    # plotting (hence "MATplotlib").  Since Matplotlib is a very large 
    # package, we will import only the subpackage we need: Pyplot.  We
    # import it as plt to reduce our typing load again.

import os # OS allows us to do fun things like changing our working directory.

from astropy.io import ascii # Astropy is a user-created library that contains
    # (most interesting for our purposes) the powerful ASCII package.  This
    # package allows us to import ASCII tables, such as tab-delineated .txt
    # files, as Python arrays.

#---------------------------------------------------------------------------

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

# In order to manipulate this, we'll need to separate the x- and y- values.
    # Copy the commands on lines 24, 26, 70, 73, and 82 into your IPython
    # console independently (striking "Enter" after each paste), then view
    # the imported Python array by copying in the command:

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

# Say we want to do some manipulations to the values from our data, such
    # as a unit conversion.  We can perform basic arithmetic on all values
    # in an array using Numpy:

# To convert the y-values from radians to degrees, for instance:

ExpYValuesConverted = (360/(2*np.pi))*ExpYValues # Note that Numpy
    # contains the value of π to many more digits of precision than we
    # care to type!  np.pi calls on this float value for you, and it can
    # be treated just like any other float in arithmetic operations.

#-------------------------------------------------------------------------

# Now, we want to make a table of this data in an existing LaTeX
    # document - but who wants to read this array and type in all of the
    # values ourselves?  Let's make Python format the LaTeX code of the
    # table for us.

# Start by creating a new array out of the existing arrays.  If we only
    # want our x-values and the converted y-values in our table, then we
    # do so by creating the object ToExport:

ToExport = [ExpXValues, ExpYValuesConverted] # This syntax appends the
    # array ExpYValuesConverted onto the end of the array ExpXValues - 
    # essentially squishing them together to make a 2D array out of two
    # 1D arrays.

# Just as there is the ascii.read function from Astropy, there is also
    # an ascii.write.  This can be used to export your data file in
    # various formats - including new ASCII tables - but here, we want
    # to use the LaTeX writer:

ascii.write(ToExport, 'Example LaTeX Table 1.txt', names=['Time (s)','Angle (deg)'], Writer=ascii.Latex)

# This will save a basic LaTeX table code under the file name "Example
    # LaTeX Table 1.txt" with the headers "Time (s)" over the x-values
    # and "Angle (deg)" over the y-values.  You can then copy-paste the
    # code out of the .txt file and into your LaTeX document.

# There are further modifications you can make to this - you can use
    # other preferences within ascii.write to do special formatting on
    # your LaTeX table to make it look a certain way, or you can save the
    # file in any of the supported formats shown in the ascii.write
    # documentation from Astropy.