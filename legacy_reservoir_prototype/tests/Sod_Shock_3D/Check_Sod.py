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
import sod

print 'First solving the analytical'
#Plot analytical results#############
gamma = 1.4
dustFrac = 0.0
npts = 1000
t = 0.075
left_state = (1,1,0)
right_state = (0.1, 0.125, 0.)

# left_state and right_state set pressure, density and u (velocity)
# geometry sets left boundary on 0., right boundary on 1 and initial
# position of the shock xi on 0.5
# t is the time evolution for which positions and states in tube should be
# calculated
# gamma denotes specific heat
positions, regions, values = sod.solve(left_state=left_state, \
	right_state=right_state, geometry=(0.0, 1.0, 0.5), t=t,
	gamma=gamma, npts=npts, dustFrac=dustFrac)

# Printing positions
print('Positions:')
for desc, vals in positions.items():
	print('{0:10} : {1}'.format(desc, vals))

# Printing p, rho and u for regions
print('Regions:')
for region, vals in sorted(regions.items()):
	print('{0:10} : {1}'.format(region, vals))

Analytical_X = []
Analytical_Yrho = []
Analytical_Yp = []
Analytical_Yv = []


Analytical_X=values['x']
Analytical_Yrho=values['rho']
Analytical_Yp=values['p']
Analytical_Yv=values['u']


print 'Now running the numerical model'

#Get path
path = os.getcwd()
binpath = path[:path.index('legacy_reservoir_prototype')] + 'bin/icferst'

os.system('make clean')
os.system('make')
os.system(binpath + ' ' + path + '/sod_shock3d_test.mpml')
#THIS SCRIPT CHECKS THE SOLUTION OBTAINED USING IC-FERST USING P2DGP1DG AND
#IT COMPARES THE SOLUTION AGAINST AN ACTUAL ANALYTICAL SOLUTION

#TOLERANCE OF THE CHECKING
#The present values are just above the values I got when writing the script
Tolerance_L1_NORM = 0.03
Tolerance_L2_NORM = 0.002

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
data_name_rho = 'Density'
data_name_p = 'Pressure'
data_name_v = 'Velocity'
#Initial and last coordinate of the probe
x0 = 0.0
x1 = 1.0

y0 = 0.0125 # 1.0/float(NUMBER)
y1 = y0

z0 = y0
z1 = y0
#Resolution of the probe
resolution = npts

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
FSrho=[]
FSP=[]
FSV=[]

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
    FSrho.append(  data.GetPointData().GetScalars(data_name_rho).GetTuple(j))

for j in range(points.GetNumberOfPoints()):
    FSP.append(  data.GetPointData().GetScalars(data_name_p).GetTuple(j))

for j in range(points.GetNumberOfPoints()):
    FSV.append(  data.GetPointData().GetScalars(data_name_v).GetTuple(j))

#Create spline curve
#tck = interpolate.splrep(Analytical_X, Analytical_Y, s=0.08)
f = interp1d(Analytical_X, Analytical_Yrho,kind ='linear')

#COnvert tuple to array
Experimental_Yrho = []
Experimental_Yp = []
Experimental_Yv = []
for item in FSrho:
    Experimental_Yrho.extend(item)
for item in FSP:
    Experimental_Yp.extend(item)
for item in FSV:
    Experimental_Yv.extend(item)


####################################### RHO #######################################

L1_sum = 0.0
L2_sum = 0.0
Infinite_Norm = 0.0


for i in range(len(Experimental_X)):
    if (i==0):#The first position is exact, so no need to interpolate
        L1_sum = L1_sum + abs(Analytical_Yrho[i] - Experimental_Yrho[i])
        L2_sum = L2_sum + (Analytical_Yrho[i] - Experimental_Yrho[i])**2
        continue
    Experimental_X[i] = Experimental_X[i] #+ 0.5#In this test case the origin is in -0.5

    position = Experimental_X[i]
    x = f(position)
    if (x==-1):
        print 'The size of the Experimental and Analytical experiments is different'
        quit

    if (abs(x - Experimental_Yrho[i])> Infinite_Norm):
        Infinite_Norm = abs(x - Experimental_Yrho[i])
    L1_sum = L1_sum + abs(x - Experimental_Yrho[i])
    L2_sum = L2_sum + (x - Experimental_Yrho[i])**2


L1_norm= L1_sum / len(Experimental_X)
L2_norm = L2_sum**0.5 / len(Experimental_X)

Passed = True

if (L1_norm > Tolerance_L1_NORM): Passed = False
if (L2_norm > Tolerance_L2_NORM): Passed = False
#Check the experiment has finished
if (AutoNumber < 16): Passed = False


if (Passed):
    print 'Sod shock works OK'
else:
    print 'Sod shock NOT work'

#######################################################
if (showPlot):
	fig, ax = plt.subplots(3, sharex=True)
	x = []
	y = []
	yp = []
	yv = []
	for i in range(len(detector)):
		x.append(float(detector[i][0]))#+0.5)#In this test case the origin is in -0.5
		y.append(float(FSrho[i][0]))
		yp.append(float(FSP[i][0]))
		yv.append(float(FSV[i][0]))
	line = plt.Line2D(x, y, color='red', linewidth=2)
	line2 = plt.Line2D(Analytical_X, Analytical_Yrho, color='blue', linewidth=2)
	#line.text.set_color('red')
	#line.text.set_fontsize(16)
	ax[0].add_line(line)
	ax[0].add_line(line2)

	line3 = plt.Line2D(x, yp, color='red', linewidth=2)
	line4 = plt.Line2D(Analytical_X, Analytical_Yp, color='blue', linewidth=2)
	ax[1].add_line(line3)
	ax[1].add_line(line4)

	line5 = plt.Line2D(x, yv, color='red', linewidth=2)
	line6 = plt.Line2D(Analytical_X, Analytical_Yv, color='blue', linewidth=2)
	ax[2].add_line(line5)
	ax[2].add_line(line6)

	plt.show()
