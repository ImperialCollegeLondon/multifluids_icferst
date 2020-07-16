#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 12:00:00 2019

@author: lv216

Script to read stat files from IC-FERST (tested for IC-FERST only)

Written to be run from terminal together with the stat files to be read.
i.e. > python read_stat.py file1.stat file2.stat etc
Temporary limit set to 5 files (set arbitrarily)

Other inputs are set as pop-up questions in terminal, to be answered by users.

Any bugs/suggestions/questions please address them to Lluis Via-Estrem (lv216@ic.ac.uk)

!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

"""

#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys
import numpy as np

# External argument call and error messages
if len(sys.argv[1])<1:
    exit('MISSING ARGUMENT: 1: *.stat input filename')

# Checking number of input stat files
nfiles=len(sys.argv)
if nfiles>5:
    exit('too many files given (Max is 5!) \n 5 is a completely arbitrary number,if you want more just get into the code...)')
alldata={}
for i in range(1,nfiles):
    
    # Checking is a stat file
    if 'stat' not in sys.argv[i]:
        exit('input file has to be *.stat')

    infilename=str(sys.argv[i])
    
    # Reading stat file
    infile=open(infilename,'r')
    inlines= infile.readlines()
    infile.close()

    k=0
    data={}
    data_vars={}

    print "List of variables in stat file:"

    for line in inlines:
        if '<' in line:
            if 'field column' in line:
                var_id=int(line.split()[1].split('"')[1])
                var_name=line.split()[2].split('"')[1]
                var_stat=line.split()[3].split('"')[1]
                name=str(var_name+'_'+var_stat)
                if 'material_phase' in line:
                    phase_name=line.split()[4].split('"')[1]
                    name=name+'_'+phase_name
                data_vars[var_id]=name
                data[name]=[]
                if i==1:
                    print var_id, name,data_vars[var_id]
        else:
            break
    
    for line in inlines:
        if '<' not in line:
            vals=line.split()
            if len(vals)==len(data):
                k+=1
                for j in range(len(data_vars)):
                    if 'ElapsedWallTime' in data_vars[j+1]:
                        data[data_vars[j+1]].append(float(vals[j])/3600.)
                    else:
                        data[data_vars[j+1]].append(float(vals[j]))
                                
    counts=range(1,k+1)
    data['counter']=counts
    print i
    alldata[i]=data
    print i,'alldata keys',alldata.keys()

flag=False

# Pop-up-in-terminal question on what variables to use for plotting
while flag==False:
    
    input_vars = raw_input("Enter the variables id you want to plot (separated by space). Max 8: \t").split()
    if len(input_vars)>=9:
        input_vars = raw_input("Ehem... Did I forget to say that max 8 variables!?!? \n Enter the variables id you want to plot (separated by space): \t").split()
        
    list_vars=[]
    for i in input_vars:
        list_vars.append(data_vars[int(i)])
    
    print "List of variables you asked for is:\n"
    for i in list_vars:
        print i
        
    print '\n'
    
    inputs = raw_input("Are you happy with that? (YES/no): \t")
    if 'no' in inputs or 'No' in inputs:
        print 'Going to ask again...'
    else:
        flag=True


inputs = raw_input("What do you want to use for the horizontal axis? ElapsedTime (1) or numer of time steps (2): (Default is 2) \t")

base_t='counter'

if '1' in inputs:
    base_t='ElapsedTime_value'

    
# Default limits for x axis
# NOTE: we are using number of time steps instead of simulation time!
# This facilitates finding mesh adapts (given that we know how often they are set to happen)
xlims=[]
for i in range(1,nfiles):
    print [min(alldata[i][base_t]),max(alldata[i][base_t])]
    xlims.append([min(alldata[i][base_t]),max(alldata[i][base_t])])

xlims=np.max(xlims,axis=0)
lims=[]
for i in range(len(list_vars)):
    lims.append([])
    
# Defining number of figures
n=len(list_vars) # n has to be <9
plot_ids=range(int(str(n)+'11'),int(str(n)+'1'+str(n+1)))

plot_flag=True

# Setting up an interactive loop to edit axis
while plot_flag==True:
    
    fig = plt.figure()
    
    for i in range(len(list_vars)):
        print list_vars
        print list_vars[i]
        # Setting up figures
        ax1=fig.add_subplot(plot_ids[i])
        ax1.set_xlabel(base_t)
        ax1.set_ylabel(list_vars[i].split('_')[0])
        ax1.set_xlim(xlims)
        if len(lims[i])>1:
            ax1.set_ylim(lims[i])
            
        ylims=[]
        for ii in range(1,nfiles):
            ax1.plot(alldata[ii][base_t], alldata[ii][list_vars[i]])
        
        # Adjusting lims for y axis according to x lims
            if plot_flag==True and xlims[0]!=min(alldata[ii][base_t]) and xlims[1]!=max(alldata[ii][base_t]):
                alldata[ii][base_t]=np.asarray(alldata[ii][base_t])
                var_0=len(alldata[ii][base_t][alldata[ii][base_t]<xlims[0]])
                var_1=len(alldata[ii][base_t][alldata[ii][base_t]<xlims[1]])
                maxY=max(alldata[ii][list_vars[i]][var_0:var_1])
                minY=min(alldata[ii][list_vars[i]][var_0:var_1])
                ydif=float(maxY)-float(minY)
                ylims.append([float(minY)-ydif*0.1,float(maxY)+ydif*0.1])
        
        if plot_flag==True and xlims[0]!=min(alldata[ii][base_t]) and xlims[1]!=max(alldata[ii][base_t]):
            print ylims
            ylims=np.max(ylims,axis=0)
            print ylims
            ax1.set_ylim(ylims)
    plt.show()
    
    # If YES, loop/script will finish, otherwise will keep iterating on whatever input is given for X axis range
    inputs = raw_input("Are you happy with the axis? (YES/no): \t")
    if 'no' in inputs or 'No' in inputs:
        X_lims = raw_input('What X axis boundaries would you like... ? \n It MUST be within the range ['+str(min(data[base_t]))+' '+str(max(data[base_t]))+'] inputs sep with space. i.e. 150 300: \t').split()
        xlims=[float(X_lims[0]),float(X_lims[1])]
    else:
        plot_flag=False   

'''
Idea from SE Zve, find derivative from data sets and spot possible peaks from it (l2norm?)
'''

