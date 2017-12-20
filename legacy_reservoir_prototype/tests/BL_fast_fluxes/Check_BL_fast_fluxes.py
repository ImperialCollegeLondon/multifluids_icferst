import csv
import matplotlib.pyplot as plt
import pylab 
import numpy
import os
from scipy.interpolate import interp1d

tC,f1C,f2C,t1C,t2C,ssum = [], [], [], [], [], []
tCA,f1CA,f2CA,t1CA,t2CA,ssumA = [], [], [], [], [], []

tolerance1 = 1e-4
tolerance2 = 1e-4
toleranceMass = 2e-5

path = os.getcwd()
binpath = path[:path.index('legacy_reservoir_prototype')] + 'bin/icferst'
os.system('rm -f ' + path+ '/*.vtu')
os.system('rm -f ' + path+ '/QuickTest_DG_outfluxes.csv')
os.system(binpath + ' ' + path + '/*mpml')

with open('QuickTest_DG_outfluxes.csv','r') as f1:
# Need a slightly modified script to correctly read this file in       
    f1.readline()
    for row in f1:
        columns = row.split(",")
        try:
           tC.append(float(columns[0]))
           f1C.append(float(columns[3]))
           f2C.append(float(columns[4]))
           t1C.append(float(columns[5]))
           t2C.append(float(columns[6]))        
        except:
           print "some error1"

f1.close()

with open('outfluxes_Control.csv','r') as f2:    
    f2.readline()
    for row in f2:
        columns = row.split(",")
        try:
            tCA.append(float(columns[0]))
            f1CA.append(float(columns[3]))
            f2CA.append(float(columns[4]))
            t1CA.append(float(columns[5]))
            t2CA.append(float(columns[6]))        
        except:
            print "some error2"

f2.close()

fn1 = interp1d(tCA,f1CA)
fn2 = interp1d(tCA,f2CA)
fn3 = interp1d(tCA,t1CA)
fn4 = interp1d(tCA,t2CA)


#### L1 ERROR NORMS

L1_sum1 = 0.
L1_sum2 = 0.
L1_sum3 = 0.
L1_sum4 = 0.

for i in range(len(tC)):
# IMPORTANT: This line was to omit the first row from the calculation (since it's text) but we're no longer reading it in.
    #if (i==0): 
        #continue    
    L1_sum1 = L1_sum1 + abs(fn1(tC[i]) - f1C[i])/len(tC)
    L1_sum2 = L1_sum2 + abs(fn2(tC[i]) - f2C[i])/len(tC)
    L1_sum3 = L1_sum3 + abs(fn3(tC[i]) - t1C[i])/len(tC)
    L1_sum4 = L1_sum4 + abs(fn4(tC[i]) - t2C[i])/len(tC)

#### L2 ERROR NORMS

L1_sum5 = 0.
L1_sum6 = 0.
L1_sum7 = 0.
L1_sum8 = 0.

for i in range(len(tC)):
# IMPORTANT: This line was to omit the first row from the calculation (since it's text) but we're no longer reading it in.
    #if (i==0): 
        #continue    
    L1_sum5 = L1_sum1 + (fn1(tC[i]) - f1C[i])**2
    L1_sum6 = L1_sum2 + (fn2(tC[i]) - f2C[i])**2
    L1_sum7 = L1_sum3 + (fn3(tC[i]) - t1C[i])**2
    L1_sum8 = L1_sum4 + (fn4(tC[i]) - t2C[i])**2

L1_sum5 = (L1_sum5**0.5)/len(tC)
L1_sum6 = (L1_sum5**0.5)/len(tC)
L1_sum7 = (L1_sum5**0.5)/len(tC)
L1_sum8 = (L1_sum5**0.5)/len(tC)


pos = len(tC)-1
totalMass = t2C[pos]+t1C[pos]
totalMass_reference = fn3(tC[pos])+fn4(tC[pos])
InjectedMass = 0.008
#print totalMass - InjectedMass
#print  L1_sum1, L1_sum2, L1_sum3, L1_sum4, L1_sum5, L1_sum6, L1_sum7, L1_sum8

#First and most important thing of this test is to check if we conserve Mass
if (abs(totalMass-InjectedMass) < toleranceMass):
    print "BL with fluxes conserves mass OK"
else:
    if(L1_sum1 < tolerance1 and L1_sum2 < tolerance1 and  L1_sum3 < tolerance1 and L1_sum4 < tolerance1 and L1_sum5 < tolerance2 and L1_sum6 < tolerance2 and L1_sum7 < tolerance2 and L1_sum8 < tolerance2) :
        print "BL with fluxes works OK"
    else:
        print "BL with fluxes does not work" 



