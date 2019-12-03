#!/usr/bin/env python

# arguments:: project vtu
# extracts flow parameters for a number of points
# from a vtu file

import vtk
import sys
from math import *
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.special import erfc
from scipy import interpolate
from scipy.interpolate import interp1d
import os



#SETTINGS ANALYTICAL SOLUTION
t = 0.1
nodes = 202
gamma = 1.
pi = 3.141596
Tolerance_L1_NORM = 0.01

print 'Running the model'
path = os.getcwd()
binpath = path[:path.index('legacy_reservoir_prototype')] + 'bin/icferst'
os.system('rm -f ' + path+ '/*.vtu')
os.system(binpath + ' ' + path + '/*mpml')


def analytical_solution(X_coords):
    coords = np.asarray(X_coords)
    term1=erfc((coords-(gamma*t))/(2.*sqrt(t)))
    term2=(np.exp(gamma*coords))*erfc((coords+(gamma*t))/(2.*sqrt(t)))
    term3=1.+(0.5*gamma*(2.-coords+(gamma*t)))
    term4=erfc((2.-coords+(gamma*t))/(2.*sqrt(t)))
    term5=gamma*(sqrt(t/pi))*np.exp(-((2.-coords+(gamma*t))**2.)/(4.*t))
    analytical=(0.5*(term1+term2))+(np.exp(gamma)*((term3*term4)-term5))

    return analytical

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
data_name = 'phase1::Temperature'
#Initial and last coordinate of the probe
x0 = 0.0
x1 = 1.0

y0 = 0.0001 # 1.0/float(NUMBER)
y1 = y0 #<==Temporary, it can handle different values

z0 = 0.0
z1 = 0.0
#Resolution of the probe
resolution = nodes


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

for i in range(len(detector)-1):
    points.InsertNextPoint(detector[i][0], detector[i][1], detector[i][2])

detectors = vtk.vtkPolyData()
detectors.SetPoints(points)    

###########Create the probe line end#############

#for i in range(len(detector)):
#    points.InsertNextPoint(detector[i][0], detector[i][1], detector[i][2])


#detectors = vtk.vtkPolyData()
#detectors.SetPoints(points)


probe = vtk.vtkProbeFilter()
probe.SetInputConnection(ugrid)


probe.SetSourceConnection(ugrid)
probe.SetInputData(detectors)
probe.Update()

data = probe.GetOutput()

for j in range(points.GetNumberOfPoints()):
    #print data.GetPointData().GetScalars(data_name)
    FS.append(  data.GetPointData().GetScalars(data_name).GetTuple(j))


#Calculate analytical solution
Experimental_X = np.linspace(0, 1.0, num=nodes, endpoint=True, retstep=False, dtype=None)


#Analytical_X = []
#Analytical_Y = []
#Analytical=file('Analytical','r')
#while True:
#    cadena=Analytical.readline()
#    if len(cadena) ==0:
#        break # EOF
#    if len(cadena) <2:
#        continue # If a line is empty       
#    lista = cadena.split()
#    Analytical_X.append(float(lista[0]))
#    Analytical_Y.append(float(lista[1]))
#Analytical.close

Analytical_Y = analytical_solution(Experimental_X)

#Compare results
#Convert tuple to array
Experimental_Y = []
for item in FS:
    Experimental_Y.extend(item)
    
#Create spline curve
f = interp1d(Experimental_X, Analytical_Y,kind ='cubic')
L1_sum = 0.0
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
     
        
        
L1_norm= L1_sum / len(Experimental_X) 
  
Passed = True

if (L1_norm > Tolerance_L1_NORM): Passed = False

if (Passed): 
    print 'Thermal_analytical_validation_two_phases OK'
else:
    print 'Thermal_analytical_validation_two_phases NOT work'

#print "L1_norm:", L1_norm

if (showPlot):
    fig, ax = plt.subplots()
#    x = []
#    Analytical = []
#    for i in range(len(detector)):
#        x.append(float(detector[i][0]))
#        Analytical.append(float(FS[i][0]))
    line = plt.Line2D(Experimental_X, Analytical_Y, color='black', linewidth=2)
    line2 = plt.Line2D(Experimental_X, FS, color='blue', linewidth=2)
    #line.text.set_color('red')
    #line.text.set_fontsize(16)
    ax.add_line(line)
    ax.add_line(line2)
    plt.autoscale(enable=True, axis='both', tight=None)
    plt.show()
