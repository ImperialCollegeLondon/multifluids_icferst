#!/usr/bin/env python

# arguments:: project vtu
# extracts flow parameters for a number of points
# from a vtu file

import vtk
import sys
from math import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.interpolate import interp1d
import os
import csv

print 'Running the model'
path = os.getcwd()
binpath = path[:path.index('legacy_reservoir_prototype')] + 'bin/icferst'
os.system('rm -f ' + path+ '/*.vtu')
os.system(binpath + ' ' + path + '/*mpml')


#TOLERANCE OF THE CHECKING
#The present values are just above the values I got when writing the script


################################AUTOMATIC STUFF###############################
Passed = False

filename = 'z_optimodel_outfluxes.csv'
phase2_out = []
with open(filename, 'rb') as csvfile:
    datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in datareader:
        try:
            phase2_out.append(float(row[18]))#Cumulative production of oil
        except:
            continue

#Check last cumulative production
diff = (abs(phase2_out[-1]) - 40051615.959)/40051615.959


print 'Compare production with stored value: ' + str(diff)
Passed = False

if (abs(diff) < 1e-2): Passed = True

#print time, temp

if (Passed): 
    print 'Well production works OK'
else:
    print 'Well production does NOT work'

