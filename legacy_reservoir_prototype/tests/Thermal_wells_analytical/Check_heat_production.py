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
Lifetime = 22

showPlot = False
################################AUTOMATIC STUFF###############################
Passed = False

filename = 'two_well_test_outfluxes.csv'
time = []
temp = []
with open(filename, 'rb') as csvfile:
    datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in datareader:
        try:
            time.append(float(row[1]))#time in years
            temp.append(float(row[8]))#this is the closest BC to the production
        except:
            continue
pos = 0
for i in range(len(temp)):
    if (temp[i] < (max(temp) -1.0) ):
        pos = i
        break

print 'Lifetime in years: ' + str(time[pos])
#Check time to produce water with lower temperature than the reservoir
if (time[pos] >= Lifetime): Passed = True
#Check the experiment has finished
if (len(temp) < 28): Passed = False

#print time, temp

if (Passed): 
    print 'Geothermal well production works OK'
else:
    print 'Geothermal well production does NOT work'

if (showPlot):
    fig, ax = plt.subplots()

    line = plt.Line2D(time, temp, color='red', linewidth=2)
    ax.add_line(line)
    ax.autoscale(enable=True, axis='both', tight=None)
    plt.show()
