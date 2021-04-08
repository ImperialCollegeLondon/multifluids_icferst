#!/usr/bin/env python3

import sys
import os


# Author: Pablo Salinas
# This script has been created to merge all the different parts of a well that might be created by trellis
# To make it work place all the parts of the well with the nastran extension into the same folder and run the script
# this will create a file with all the wells merged.

class Nodes():
    def __init__(self):
        self.X = 0.
        self.Y = 0.
        self.Z = 0.
class Edges():
    def __init__(self):
        self.node1 = 0
        self.node2 = 0


print('Converting file...')

#Get path

path = os.getcwd()
########INPUT##############
input_files = []
#By default merge all the files inside one folder
#First we create the list of files
for files in os.listdir(path):
    if files.endswith(".bdf"):
        cadena = os.path.join(path+'/', files)
        if (not (('Merged_well.bdf') in cadena)):
            input_files.append(cadena)
###########AUTOMATIC############


##This script is to introduce variables into a structure type variable


#Final list with everything
Final_Lnodes = [Nodes()]
del Final_Lnodes[0]
Final_Ledges = [Edges()]
del Final_Ledges[0]
for files in input_files:
    File=open(files,'r')
    print('Reading file: ', files)
    #Create nodes list and edges list from all the files
    Lnodes = [Nodes()]
    del Lnodes[0]
    Ledges = [Edges()]
    del Ledges[0]
    conversor = []

    while True:
        auxNodes = Nodes()
        auxEdges = Edges()
        cadena=File.readline()
        #Skip the header
        if (cadena[0:4] == "GRID"):
            auxNodes.X = float(cadena[24:32])
            auxNodes.Y = float(cadena[32:40])
            auxNodes.Z = float(cadena[40:49])
            Lnodes.append(auxNodes)
            conversor.append(int(cadena[8:12]))
 
        #Now read the edges
        if (cadena[0:4] == "CROD"): 
            #Only interested in the last two
            aux = cadena[6:len(cadena)].split()
            auxEdges.node1 = int(aux[-1])
            auxEdges.node2 = int(aux[-2])
            Ledges.append(auxEdges)
        #When reaching end of file exit
        if len(cadena) ==0:
            break # EOF
    #Now normalise the edge list so it is consequent with the position in the Nodes list
    start = len(Final_Ledges)
    for edge in Ledges:
        for i in range(len(conversor)):
            if (conversor[i] == edge.node1):
                edge.node1 = i + start
                continue
            if (conversor[i] == edge.node2):
                edge.node2 = i + start
                continue



    #Copy values to the final list
    for n in Lnodes:
        Final_Lnodes.append(n)
        #print n.X, n.Y, n.Z
    for edge in Ledges:
        Final_Ledges.append(edge)
        #print edge.node1, edge.node2

        #Change general integers
        
    #    for var in  names:
            #This search and replace by ignoring case and also looks for the whole word
            #First in global integer
    #        auxC = re.compile("\\b" + re.escape(var)+"\\b", re.IGNORECASE)
    #        cadena = auxC.sub(Struc_name+var, cadena)
    File.close()


#for n in Final_Lnodes:
#    print n.X, n.Y, n.Z
#for edge in Final_Ledges:
#    print edge.node1, edge.node2
#Output.write(cadena)

#Use one file as template to mimic it and create the merged file
Output = open('Merged_well.bdf', "w")
File=open(input_files[0],'r')
first_GRID = True
first_CROD = True
while True:
    auxNodes = Nodes()
    auxEdges = Edges()
    cadena=File.readline()
    #Skip the header
    if (cadena[0:4] == "GRID" and first_GRID):
        first_GRID = False
        for n in Final_Lnodes:
            ##This is what happens when you use a retard file format...
            cadena2 = cadena[0:24]
            cad = str(n.X)
            #Need to fill up up to 8 digits
            while len(cad)<8:
                cad +=' '
            cadena2+= cad[0:8]
            cad = str(n.Y)
            while len(cad)<8:
                cad +=' '
            cadena2+= cad[0:8]
            cad = str(n.Z)
            while len(cad)<8:
                cad +=' '
            cadena2+= cad[0:8]
            Output.write(cadena2+'\n')

    #Now read the edges
    if (cadena[0:4] == "CROD" and first_CROD): 
        first_CROD = False
        for edge in Final_Ledges:
            ##This is what happens when you use a retard file format...
            cadena2 = cadena[0:24]
            cad = str(edge.node1)
            #Need to fill up up to 8 digits
            while len(cad)<8:
                cad +=' '
            cadena2+= cad[0:8]
            cad = str(edge.node2)
            while len(cad)<8:
                cad +=' '
            cadena2+= cad[0:8]
            Output.write(cadena2+'\n')


        #Only interested in the last two
        aux = cadena[6:len(cadena)].split()
        auxEdges.node1 = int(aux[-1])
        auxEdges.node2 = int(aux[-2])
        Ledges.append(auxEdges)
    if (cadena[0:4] != "GRID" and cadena[0:4] != "CROD"):
        Output.write(cadena)
    #When reaching end of file exit
    if len(cadena) ==0:
        break # EOF
File.close()
Output.close()
