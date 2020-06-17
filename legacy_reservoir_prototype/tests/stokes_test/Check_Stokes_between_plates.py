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
binpath = path[:path.index('legacy_reservoir_prototype')] + 'bin/icferst'
os.system('rm -f ' + path+ '/*.vtu')
os.system(binpath + ' ' + path + '/stokes_test.mpml')


#TOLERANCE OF THE CHECKING
#The present values are just above the values I got when writing the script
Tolerance_L1_NORM = 0.00012

#RETRIEVE AUTOMATICALLY THE LAST VTU FILE
#AutoNumber = 0
#for files in os.listdir(path):
#    if files.endswith(".vtu"):
#        pos = files.rfind('_')
#        pos2 = files.rfind('.')
#        AutoFile = files[:pos]
#        AutoNumber = max(AutoNumber, int(files[pos+1:pos2]))


AutomaticFile = "stokes_test"
AutomaticVTU_Number = 1

#Plot the results in 2d?
showPlot = False

#NAME OF THE VARIABLE YOU WANT TO EXTRACT DATA FROM
data_name = 'Velocity'

#Initial and last coordinate of the probe
x0 = 0.5
x1 = 0.5

y0 = 0.0 #
y1 = 0.1 #

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

ugrid = reader.GetOutputPort()
#ugrid.Update()

###########Create the probe line#############
        
detector = []

hx = (x1 - x0) / resolution
hy = (y1 - y0) / resolution
hz = (z1 - z0) / resolution

Experimental_X = []
Experimental_X = []
for i in range(resolution):
    detector.append([hx * i + x0, hy * i + y0 , hz * i + z0])
    Experimental_X.append(hy * i + y0)

#print 'using',len(detector),'detectors'


###########Create the probe line end#############


points = vtk.vtkPoints()
points.SetDataTypeToDouble()

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
#    FS.append(  data.GetPointData().GetScalars(data_name).GetTuple(j))

#So far we have the information from the analytical result


Analytical_X = []
Analytical_Y = []

Analytical_X1 = np.linspace(0.0, 0.1, num=resolution, retstep=True)
Analytical_X = Analytical_X1[0][:]
#Flow between plates formula for pressure drop of 1, Lenght = 1, viscosity = 1, H = 0.2
for i in range(resolution):
    Analytical_Y.append( 0.5*float(Analytical_X[i])*(0.1 - float(Analytical_X[i])))

#Create spline curve
#tck = interpolate.splrep(Analytical_X, Analytical_Y, s=0.08)
f = interp1d(Analytical_X, Analytical_Y,kind ='linear')

#COnvert tuple to array
Experimental_Y = []
for item in FS:
    auxR =  item[0]
    Experimental_Y.append(auxR)
    #Experimental_Y.extend(auxR)

L1_sum = 0.0
L2_sum = 0.0
L1_sum_shock_front = 0.0
L2_sum_shock_front = 0.0
N_shock = 0
Infinite_Norm = 0.0
for i in range(len(Experimental_X)):
    if (i==0):#The first position is exact, so no need to interpolate
        L1_sum = L1_sum + abs(Analytical_Y[i] - Experimental_Y[i])
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
    if (abs(x - Experimental_Y[i])>1/100000000):
        N_shock = N_shock + 1
        L1_sum_shock_front = L1_sum_shock_front + abs(x - Experimental_Y[i])
        
L1_norm= L1_sum / len(Experimental_X) 

Passed = True

if (L1_norm > Tolerance_L1_NORM): Passed = False
#Check the experiment has finished
#if (AutoNumber < 20): Passed = False

#print L1_norm
if (Passed): 
    print 'BL works OK'
else:
    print 'BL does NOT work'


if (showPlot):
    fig, ax = plt.subplots()
    x = []
    y = []
    for i in range(len(detector)):
        x.append(float(detector[i][1]))
        y.append(float(FS[i][0]))
    line = plt.Line2D(x, y, color='red', linewidth=2)
    line2 = plt.Line2D(Analytical_X, Analytical_Y, color='blue', linewidth=2)
    #line.text.set_color('red')
    #line.text.set_fontsize(16)
    ax.add_line(line)
    ax.add_line(line2)
    ax.set(xlim=(0, 0.1), ylim=(0, 2e-3))
    plt.show()
