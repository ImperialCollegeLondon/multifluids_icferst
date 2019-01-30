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

filename = 'two_well_test_outfluxes.csv'
phase1_in = []
phase2_out = []
with open(filename, 'rb') as csvfile:
    datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in datareader:
        try:
            phase1_in.append(float(row[10]))#Cumulative injection of water
            phase2_out.append(float(row[17]))#Cumulative production of oil
        except:
            continue

#Check last cumulative production
diff = abs(phase1_in[2] + phase2_out[2])/abs(phase2_out[2])


print 'In-out difference after 15 years: ' + str(diff)
Passed = False
#Check time to produce water with lower temperature than the reservoir
if (abs(diff) < 1e-5): Passed = True

#print time, temp

if (Passed): 
    print 'Well production works OK'
else:
    print 'Well production does NOT work'

