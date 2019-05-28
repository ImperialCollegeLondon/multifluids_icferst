#!/usr/bin/env python

import vtk
import glob
import sys
import os
import vtktools
import operator
import numpy
import pylab
import fluidity_tools

TOLERANCE_H=0.05
TOLERANCE_P=150
plotting=False

# first extract the water gauge data from the vtus
def get_water_depths(filelist, xarray, delta):
  results = []
  for f in filelist:
    try:
      os.stat(f)
    except:
      print "No such file: %s" % f
      sys.exit(1)
    
    y = numpy.arange(delta/2.0,2.0+delta/2.0,delta)[:,numpy.newaxis]
    
    num = int(f.split(".vtu")[0].split('_')[-1])
    vtu = vtktools.vtu(f)
    for name in vtu.GetFieldNames(): 
      if name.endswith("Time"): 
        time = max(vtu.GetScalarRange(name))
        break
    waterdepths = []
    waterdepths.append(num)
    waterdepths.append(time)
    for x in range(len(xarray)):
      coordinates = numpy.concatenate((numpy.ones((len(y),1))*xarray[x], numpy.zeros((len(y),1))+0.5, y),1)
      waterdepths.append(sum(vtu.ProbeData(coordinates, "Component1::ComponentMassFractionPhase1"))[0]*delta)
    
    results.append(waterdepths)
  
  results.sort(key=operator.itemgetter(1))
  results = numpy.array(results)
  return results

print 'Running the model'

#Get path

path = os.getcwd()
binpath = path[:path.index('legacy_reservoir_prototype')] + 'bin/icferst'
os.system('rm -f ' + path+ '/*.vtu')
os.system(binpath + ' ' + path + '/*mpml')

xarray = [0.5, 1.5, 2.15]

filelist = glob.glob("cwc_[0-9]*.vtu")
results = get_water_depths(filelist, xarray, 0.01)

# now let's plot the data
warray = ["H1", "H2", "H3"]
parray = ["P2"]

time = results[:,1]

H_check=True
P_check=True

# first the water gauges
for x in range(len(xarray)):
  experiment = numpy.load(warray[x]+".npy")

  if plotting==True:
	  
	  pylab.figure(x)
	  pylab.title(warray[x]+" water gauge at "+str(xarray[x])+"m.")
	  pylab.xlabel('Time (s)')
	  pylab.ylabel('Water Depth (m)')
	  pylab.plot(experiment[:,1], experiment[:,2], marker = 'o', markerfacecolor='white', markersize=6, markeredgecolor='black', linestyle="None")
	  pylab.plot(time, results[:,2+x], color='black', linestyle="dashed")
	  pylab.axis([0.0, 1.0, 0.0, 0.75])
	  pylab.legend(("Experiment", "Model"), loc="upper left")
	  pylab.savefig("water_gauge_"+warray[x]+".png")

  #print str(warray[x]),"TOL=",TOLERANCE_H,", ERROR=",abs(numpy.std(numpy.array(experiment[:,2])-numpy.array(results[:,2+x])))
  H_check=abs(numpy.std(numpy.array(experiment[:,2])-numpy.array(results[:,2+x])))<TOLERANCE_H
  if H_check==False:
      print "H_check=",H_check
      break
  
ts=time[1:] #ignore 0

# then the pressure gauges - this takes it data from the detectors so no
# need for extraction from the vtus
results = fluidity_tools.stat_parser("cwc.detectors")
time = results["ElapsedTime"]["value"]

for p in range(len(parray)):
	
  data_o = results["phase1"]["Pressure"][parray[p]]
  if "HydrostaticPressure" in results["phase1"]:
      data_o+=results["phase1"]["HydrostaticPressure"][parray[p]]
  experiment = numpy.load(parray[p]+".npy")
  data=experiment.item(0)["phase1"]["Pressure"][parray[p]]
  if "HydrostaticPressure" in experiment.item(0)["phase1"]:
      data+=experiment.item(0)["phase1"]["HydrostaticPressure"][parray[p]]

# We need to calculate pressures at dump file steps
  data_o2=[]
  data_2=[]
  tol=0.01
  for i in ts:
      io2=next(j for j, _ in enumerate(time) if numpy.isclose(_, i, tol))
      i2=next(j for j, _ in enumerate(experiment.item(0)["ElapsedTime"]["value"]) if numpy.isclose(_, i, tol))
      data_o2.append(data_o[io2])
      data_2.append(data[i2])
      
  if plotting==True:
	  pylab.figure(p+len(xarray))
	  pylab.title(parray[p]+' pressure gauge.')
	  pylab.xlabel('Time (s)')
	  pylab.ylabel('Pressure (Pa)')
	  pylab.plot(ts, data_2, marker = 'o', markerfacecolor='white', markersize=6, markeredgecolor='black', linestyle="None")
	  pylab.plot(ts, data_o2, color='black', linestyle="dashed")
	  pylab.axis([0.0, 1.0, -1000., 10000.])
	  pylab.legend(("Experiment", "Model"), loc="upper left")
	  pylab.savefig("pressure_gauge_"+parray[p]+".png")
  
  print str(parray[p]),"TOL=",TOLERANCE_P,", ERROR=",abs(numpy.std(numpy.array(data_o2)-numpy.array(data_2)))
  P_check=abs(numpy.std(numpy.array(data_o2)-numpy.array(data_2)))<TOLERANCE_P
  if P_check==False:
      print "P_check=",P_check
      break

#pylab.show()
Passed=True
if H_check==False: Passed=False
if P_check==False: Passed=False

if (Passed): 
    print 'CWC works OK'
else:
    print 'CWC does NOT work'
