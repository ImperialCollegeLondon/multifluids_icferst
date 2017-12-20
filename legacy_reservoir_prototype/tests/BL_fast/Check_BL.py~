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


print 'Running the model'

#Get path

path = os.getcwd()
binpath = path[:path.index('legacy_reservoir_prototype')] + 'bin/multiphase_prototype'
os.system('rm -f ' + path+ '/*.vtu')
os.system(binpath + ' ' + path + '/BL_fast.mpml')
#THIS SCRIPT CHECKS THE SOLUTION OBTAINED USING IC-FERST USING P2DGP1DG AND 
#A STRUCTURED MESH OF 30 ELEMENTS IN THE X-DIRECTION
#IT COMPARES THE SOLUTION AGAINST AN ACTUAL ANALYTICAL SOLUTION

#TOLERANCE OF THE CHECKING
#The present values are just above the values I got when writing the script
Tolerance_L1_NORM = 0.0183
Tolerance_L2_NORM = 0.001125

AutomaticLine = 0

#RETRIEVE AUTOMATICALLY THE LAST VTU FILE
AutoNumber = 0
for files in os.listdir(path):
    if files.endswith(".vtu"):
        pos = files.rfind('_')
        pos2 = files.rfind('.')
        AutoFile = files[:pos]
        AutoNumber = max(AutoNumber, int(files[pos+1:pos2]))


AutomaticFile = AutoFile
AutomaticVTU_Number = AutoNumber

#Plot the results in 2d?
showPlot = False

#NAME OF THE VARIABLE YOU WANT TO EXTRACT DATA FROM
data_name = 'phase1::PhaseVolumeFraction'

#Initial and last coordinate of the probe
x0 = -0.5
x1 = 0.5

y0 = 0.0 # 1.0/float(NUMBER)
y1 = y0 #<==Temporary, it can handle different values

z0 = 0.0
z1 = 0.0
#Resolution of the probe
resolution = 1000


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
ugrid = reader.GetOutput()
ugrid.Update()
###########Create the probe line#############
#Get bounds of the domain
Aux = ugrid.GetBounds()

if (AutomaticLine > 0):
    x0 = float(Aux[0])
    x1 = float(Aux[1])
    y0 = float(Aux[2])
    y1 = float(Aux[3])
    z0 = float(Aux[4])
    z1 = float(Aux[5])  

if (AutomaticLine == 1):#Straight line across the middle
    x0 = float(Aux[0])
    x1 = float(Aux[1])
    y0 = max(float(Aux[2])/2., float(Aux[3])/2.)
    y1 = y0
    z0 = max(float(Aux[4])/2., float(Aux[5])/2.)
    z1 = z0

#Nothing to do for diagonal, it is already prepared
        
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
probe.SetSource(ugrid)
probe.SetInput(detectors)
probe.Update()

data = probe.GetOutput()

for j in range(points.GetNumberOfPoints()):
    FS.append(  data.GetPointData().GetScalars(data_name).GetTuple(j))

#So far we have the information from the analytical result


Analytical_X = []
Analytical_Y = []
Analytical=file('Analytical','r')


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
    Experimental_X[i] = Experimental_X[i] + 0.5#In this test case the origin is in -0.5

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
        
        
L1_norm= L1_sum / len(Experimental_X) 
L2_norm = L2_sum**0.5 / len(Experimental_X)    

Passed = True

if (L1_norm > Tolerance_L1_NORM): Passed = False
if (L2_norm > Tolerance_L2_NORM): Passed = False
#print L1_norm, L2_norm
if (Passed): 
    print 'BL works OK'
else:
    print 'BL does NOT work'


if (showPlot):
    fig, ax = plt.subplots()
    x = []
    y = []
    for i in range(len(detector)):
        x.append(float(detector[i][0])+0.5)#In this test case the origin is in -0.5
        y.append(float(FS[i][0]))
    line = plt.Line2D(x, y, color='red', linewidth=2)
    line2 = plt.Line2D(Analytical_X, Analytical_Y, color='blue', linewidth=2)
    #line.text.set_color('red')
    #line.text.set_fontsize(16)
    ax.add_line(line)
    ax.add_line(line2)
    plt.show()
