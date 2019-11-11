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
phase1_out = []
phase1_in = []
with open(filename, 'rb') as csvfile:
    datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in datareader:
        try:
            phase2_out.append(float(row[18]))#Cumulative production of oil
            phase1_out.append(float(row[17]))#Cumulative production of water
            phase1_in.append(float(row[9]))#Cumulative injection of water
        except:
            continue

#Check last cumulative production
#diff = abs(phase2_out[-1] - 41048807.176)/41048807.176
diff = abs(phase2_out[-1] + phase1_out[-1] + phase1_in[-1])/abs(phase1_in[-1]) * 100

print 'Compare production with injection: ' + str(diff)
Passed = False
#Below 0.5% we are happy
if (abs(diff) < 0.5 and phase2_out > 0 ): Passed = True

#print time, temp

if (Passed): 
    print 'Well production works OK'
else:
    print 'Well production does NOT work'

