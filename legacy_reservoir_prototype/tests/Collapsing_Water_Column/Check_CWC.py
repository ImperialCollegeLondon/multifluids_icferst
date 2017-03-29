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


def getAnalytical_interpolated( Analytical_X, Analytical_Y, position):
    #Returns the physical line and line to which certain edge belong
    getAnalytical_interpolated = -1
    k = 0
    #Values are ordered in an increase fashion
    for i in range(len(Analytical_X)):
        if (Analytical_X[i]>=position):
            k = i
            break

    a = (Analytical_Y[k-1] - Analytical_Y[k])/(Analytical_X[k-1] - Analytical_X[k])    
    getAnalytical_interpolated = a * (position-Analytical_X[k-1]) + Analytical_Y[k-1]
    
    return getAnalytical_interpolated 

print 'Running the model'

#Get path

path = os.getcwd()
binpath = path[:path.index('legacy_reservoir_prototype')] + 'bin/icferst'
os.system('rm -f ' + path+ '/*.vtu')
os.system(binpath + ' ' + path + '/*mpml')
#THIS SCRIPT CHECKS THE SOLUTION OBTAINED USING IC-FERST USING P2DGP1DG AND 
#A STRUCTURED MESH OF 30 ELEMENTS IN THE X-DIRECTION
#IT COMPARES THE SOLUTION AGAINST AN ACTUAL ANALYTICAL SOLUTION

#TOLERANCE OF THE CHECKING
#The present values are just above the values I got when writing the script
#The errors seem big but that is 
#because the MAXIMUM pressure is about 10^6
Tolerance_L1_NORM = 0.00175
Tolerance_L2_NORM = 0.000175


#The name of the file and number can be introduced here
#To use this, don't introduce a command argument
AutomaticFile = 'cwc'
AutomaticVTU_Number = 40

#Plot the results in 2d?
showPlot = False

#NAME OF THE VARIABLE YOU WANT TO EXTRACT DATA FROM
data_name = 'phase1::Pressure'

#Initial and last coordinate of the probe
x0 = 0.0
x1 = 1.0

y0 = 0.1
y1 = 0.1

z0 = 0.0
z1 = 0.0
#Resolution of the probe
resolution = 500


#TO EXTRACT VECTORIAL VARIABLES,
# PART OF THE CODE HAS TO BE ACTIVATED AND MODIFIED

   

################################AUTOMATIC STUFF###############################

if (len(sys.argv)>1):
    filename   = sys.argv[1]
    vtu_number = int(sys.argv[2])
else:
    filename = AutomaticFile
    vtu_number = int(AutomaticVTU_Number)
    
#print 'reading data...'

U=[]
T=[]
S=[]
FS=[]

# parallel
#reader = vtk.vtkXMLPUnstructuredGridReader()
#reader.SetFileName(filename+'_'+str(vtu_number)+'.pvtu')

# serial
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(filename+'_'+str(vtu_number)+'.vtu')

#reader.Update()

ugrid = reader.GetOutputPort()
#ugrid.Update()

###########Create the probe line#############
        
detector = []

hx = (x1 - x0) / resolution
hy = (y1 - y0) / resolution
hz = (z1 - z0) / resolution

Experimental_X = []
for i in range(resolution+1):
    detector.append([hx * i + x0, hy * i + y0, hz * i + z0])
    Experimental_X.append(hx * i + x0)

#print 'using',len(detector),'detectors'

points = vtk.vtkPoints()
points.SetDataTypeToDouble()

for i in range(len(detector)):
    points.InsertNextPoint(detector[i][0], detector[i][1], detector[i][2])

detectors = vtk.vtkPolyData()
detectors.SetPoints(points)    

###########Create the probe line end#############

for i in range(len(detector)):
    points.InsertNextPoint(detector[i][0], detector[i][1], detector[i][2])


detectors = vtk.vtkPolyData()
detectors.SetPoints(points)


probe = vtk.vtkProbeFilter()
probe.SetInputConnection(ugrid)


probe.SetSourceConnection(ugrid)
probe.SetInputData(detectors)
probe.Update()

data = probe.GetOutput()

for j in range(points.GetNumberOfPoints()):
    FS.append(  data.GetPointData().GetScalars(data_name).GetTuple(j))

#So far we have the information from the analytical result


Analytical_X = []
Analytical_Y = []
Analytical=file('Semi-Analytical','r')


while True:
    cadena=Analytical.readline()

    if len(cadena) ==0:
        break # EOF
    if len(cadena) <2:
        continue # If a line is empty       
    lista = cadena.split()
    Analytical_X.append(float(lista[0]))
    Analytical_Y.append(float(lista[1]))

Analytical.close

#Create spline curve
#tck = interpolate.splrep(Analytical_X, Analytical_Y, s=0.08)
f = interp1d(Analytical_X, Analytical_Y,kind ='linear')

#COnvert tuple to array
Experimental_Y = []
for item in FS:
    Experimental_Y.extend(item)



L1_sum = 0.0
L2_sum = 0.0
L1_sum_shock_front = 0.0
L2_sum_shock_front = 0.0
N_shock = 0
Infinite_Norm = 0.0
for i in range(len(Experimental_X)):
    if (i==0):#The first position is exact, so no need to interpolate
        L1_sum = L1_sum + abs(Analytical_Y[i] - Experimental_Y[i])
        L2_sum = L2_sum + (Analytical_Y[i] - Experimental_Y[i])**2
        continue
    
    position = Experimental_X[i]
#    x = getAnalytical_interpolated( Analytical_X, Analytical_Y, position)
    x = f(position)
    if (x==-1):
        print 'The size of the Experimental and Analytical experiments is different'
        quit

    if (abs(x - Experimental_Y[i])> Infinite_Norm):
        Infinite_Norm = abs(x - Experimental_Y[i])
    L1_sum = L1_sum + abs(x - Experimental_Y[i])
    L2_sum = L2_sum + (x - Experimental_Y[i])**2
    if (abs(x - Experimental_Y[i])>1/100000000):
        N_shock = N_shock + 1
        L1_sum_shock_front = L1_sum_shock_front + abs(x - Experimental_Y[i])
        L2_sum_shock_front = L2_sum_shock_front + (x - Experimental_Y[i])**2      
        
        
L1_norm= L1_sum / (len(Experimental_X)*max(Analytical_Y))
L2_norm = L2_sum**0.5 / (len(Experimental_X)*max(Analytical_Y))

Passed = True

if (L1_norm > Tolerance_L1_NORM): Passed = False
if (L2_norm > Tolerance_L2_NORM): Passed = False
#print L1_norm, L2_norm
if (Passed): 
    print 'CWC works OK'
else:
    print 'CWC does NOT work'

if (showPlot):
    fig, ax = plt.subplots()
    x = []
    y = []
    for i in range(len(detector)):
        x.append(float(detector[i][0]))
        y.append(float(FS[i][0]))
    line = plt.Line2D(x, y, color='red', linewidth=4)
    line2 = plt.Line2D(Analytical_X, Analytical_Y, color='blue', linewidth=2)
    #line.text.set_color('red')
    #line.text.set_fontsize(16)
    ax.add_line(line)
    ax.add_line(line2)
    plt.show()
