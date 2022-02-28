#!/usr/bin/env python3

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


print('Running the model')

#Get path

path = os.getcwd()
binpath = path[:path.index('ICFERST')] + 'bin/icferst'
os.system('rm -f ' + path+ '/*.vtu')
os.system(binpath + ' ' + path + '/*.mpml')


#TOLERANCE OF THE CHECKING
Tolerance = 1e-5

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
Passed = True
for test in range(2):
  #NAME OF THE VARIABLE YOU WANT TO EXTRACT DATA FROM
  if test == 0:
    data_name = 'phase1::PhaseVolumeFraction'
  else:
    data_name = 'phase1::Concentration'

  #Initial and last coordinate of the probe
  x0 = -0.5
  x1 = 0.5

  y0 = 0.0 # 1.0/float(NUMBER)
  y1 = y0 #<==Temporary, it can handle different values

  z0 = 0.0
  z1 = 0.0
  #Resolution of the probe
  resolution = 1000


  ################################AUTOMATIC STUFF###############################

  try:
      filename   = sys.argv[1]
      vtu_number = int(sys.argv[2])
  except:
      filename = AutomaticFile
      vtu_number = int(AutomaticVTU_Number)
      
  #print 'reading data...'

  U=[]
  T=[]
  S=[]
  FS=[]

  # serial
  reader = vtk.vtkXMLUnstructuredGridReader()
  reader.SetFileName(filename+'_'+str(vtu_number)+'.vtu')


  ugrid = reader.GetOutputPort()

  ###########Create the probe line#############
  detector = []
  hx = (x1 - x0) / resolution
  hy = (y1 - y0) / resolution
  hz = (z1 - z0) / resolution
  Experimental_X = []
  for i in range(resolution+1):
      detector.append([hx * i + x0, hy * i + y0, hz * i + z0])
      Experimental_X.append(hx * i + x0)

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


  #Convert tuple to array
  Experimental_Y = []
  for item in FS:
      Experimental_Y.extend(item)
  #We test here that every value  with a lower X coordinate has a higher or equal value than the next
  for i in range(resolution-1):
    if (Experimental_Y[i]< 1e-10): continue 
    if (Experimental_Y[i] - Experimental_Y[i+1] < -Tolerance ): 
      #print(Experimental_Y[i], Experimental_Y[i+1], Experimental_Y[i] - Experimental_Y[i+1], i)
      Passed = False

#Check the experiment has finished
if (AutoNumber < 20): Passed = False


if (Passed): 
    print('Monotonicity for Tracer and Saturation works OK')
else:
    print('Monotonicity for Tracer and Saturation does NOT work')


