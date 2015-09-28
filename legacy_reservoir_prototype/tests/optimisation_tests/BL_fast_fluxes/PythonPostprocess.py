import csv
import matplotlib.pyplot as plt
import pylab 
import numpy
import os
from scipy.interpolate import interp1d

tC,f1C,f2C,t1C,t2C,ssum = [], [], [], [], [], []
tCA,f1CA,f2CA,t1CA,t2CA,ssumA = [], [], [], [], [], []

tolerance1 = 0.15
tolerance2 = 0.015


path = os.getcwd()
binpath = path[:path.index('legacy_reservoir_prototype')] + 'bin/multiphase_prototype'
os.system('rm -f ' + path+ '/*.vtu')
os.system('rm -f ' + path+ '/outfluxes.csv')
os.system(binpath + ' ' + path + '/*mpml')

with open('outfluxes.csv','r') as f1:
# Need a slightly modified script to correctly read this file in       
    f1.readline()
    for row in f1:
        columns = row.split(",")
        try:
            tC.append(float(columns[1]))
            f1C.append(float(columns[2]))
            f2C.append(float(columns[3]))
            t1C.append(float(columns[4]))
            t2C.append(float(columns[5]))        
        except:
            print "some error"

f1.close()

with open('outfluxesHighRes.txt','r') as f2:
# This file is still in the old format, so we read it in using the old method       
    reader = csv.reader(f2)
    for row in reader:    
        columns = row[0].split()
        try:
            tCA.append(float(columns[0]))
            f1CA.append(float(columns[1]))
            f2CA.append(float(columns[2]))
            t1CA.append(float(columns[3]))
            t2CA.append(float(columns[4]))        
        except:
            print "some error"

f2.close()

#a=numpy.array(t1C)
#b=numpy.array(t2C)

fn1 = interp1d(tCA,f1CA)
fn2 = interp1d(tCA,f2CA)
fn3 = interp1d(tCA,t1CA)
fn4 = interp1d(tCA,t2CA)

L1_sum1 = 0
L1_sum2 = 0
L1_sum3 = 0
L1_sum4 = 0

for i in range(len(tC)):
# This line was to omit the first row from the calculation (since it's text) but we're no longer reading it in.
    #if (i==0): 
        #continue    
    L1_sum1 = L1_sum1 + abs(fn1(tC[i]) - f1C[i])/len(tC)
    L1_sum2 = L1_sum2 + abs(fn2(tC[i]) - f2C[i])/len(tC)
    L1_sum3 = L1_sum3 + abs(fn3(tC[i]) - t1C[i])/len(tC)
    L1_sum4 = L1_sum4 + abs(fn4(tC[i]) - t2C[i])/len(tC)


#print  L1_sum1, L1_sum2, L1_sum3, L1_sum4

if(L1_sum1 < tolerance1 and L1_sum2 < tolerance1 and  L1_sum3 < tolerance2 and L1_sum4 < tolerance2) :
    print "BL with fluxes works OK"
else:
    print "BL with fluxes failed" 

# REMEMBER TO CHECK IN THE OUTFLUXES FILES THAT YOU DON'T HAVE REPEATED SIMULATIONS PRINTING TO THE SAME FILE!!!!!!

#plt.xlabel('Time (s)',fontsize=25)
#plt.ylabel('Phase 1 Flux Across Outlet',fontsize=25)
#plt.xticks(fontsize = 20)
#plt.yticks(fontsize = 20)
#plt.plot(tC, f1C, 'ro')
#plt.show()


#plt.xlabel('Time (s)',fontsize=25)
#plt.ylabel('Phase 2 Flux Across Outlet',fontsize=25)
#plt.xticks(fontsize = 20)
#plt.yticks(fontsize = 20)
#plt.plot(tC, f2C, 'ro')
#plt.show()


#plt.xlabel('Time (s)',fontsize=25)
#plt.ylabel('Time Integrated Phase 1 Outflow',fontsize=25)
#plt.xticks(fontsize = 20)
#plt.yticks(fontsize = 20)
#plt.plot(tC, t1C, 'ro')
#plt.show()


#plt.xlabel('Time (s)',fontsize=25)
#plt.ylabel('Time Integrated Phase 2 Outflow',fontsize=25)
#plt.xticks(fontsize = 20)
#plt.yticks(fontsize = 20)
#plt.plot(tC, t2C, 'ro')
#plt.show()

#plt.xlabel('Time (s)',fontsize=25)
#plt.ylabel('Total flux across Outflow',fontsize=25)
#plt.xticks(fontsize = 20)
#plt.yticks(fontsize = 20)
#plt.plot(tC, a+b, 'go')
#plt.show()

#plt.xlabel('Time (s)',fontsize=25)
#plt.xticks(fontsize = 20)
#plt.yticks(fontsize = 20)
#plt.plot(tC, a+b, 'go', label = "sum")
#plt.plot(tC, t1C, 'ro', label = "phase1 outflux")
#plt.plot(tC, t2C, 'bo', label = "phase2 outflux")
#pylab.legend(loc='upper left')
#plt.show()

