#!/usr/bin/env python

#Created by James Percival
#modified by:
#Carl Jacquemyn
#Pablo Salinas
#This code has LGPL license
import vtk
import sys

fname=sys.argv[1]

r=vtk.vtkExodusIIReader()
if (fname[-2:] == '-h'):
    print 'This script converts a general exodus ii file into ASCII .msh file format. It requires the vtk python library to be installed on the system.'
    print 'The input should be the name of the file, and the output is that same name .msh. Example: python exodus2msh.py test.e'
    print 'To convert to .msh binary format use the option /geometry/create_binary_msh in Diamond.'
    exit()
elif (fname[-2:] != '.e'):
    fname+='.e'
print 'Converting the input exodus ii file into .msh format...'
r.SetFileName(fname)
r.UpdateInformation()
r.GenerateGlobalNodeIdArrayOn()
r.GenerateGlobalElementIdArrayOn()
#r.ExodusModelMetadataOn()
#r.PackExodusModelOntoOutputOn()
for i in range(r.GetNumberOfSideSetArrays()):
    name=r.GetSideSetArrayName(i)
    r.SetSideSetArrayStatus(name,1)
r.Update()


data=r.GetOutput()

node_dict={}
ele_face_dict=[]
ele_vol_dict=[]

n=1
def f():
    global n
    a=n
    n=n+1
    return a

for j in range(data.GetBlock(4).GetNumberOfBlocks()):
    ug=data.GetBlock(4).GetBlock(j)
    for k in range(ug.GetNumberOfPoints()):
        lnodes={}
        p=ug.GetPoint(k)
        N=node_dict.get(p)
        if N==None:
            node_dict[p]=n
            n+=1
        lnodes[k]=node_dict.get(p)
    for k in range(ug.GetNumberOfCells()):
        cp=ug.GetCell(k).GetPoints()
        ele_face_dict.append((j+1,
                              node_dict[cp.GetPoint(0)],
                              node_dict[cp.GetPoint(1)],
                              node_dict[cp.GetPoint(2)]))
                                        

for j in range(data.GetBlock(0).GetNumberOfBlocks()):
    ug=data.GetBlock(0).GetBlock(j)
    for k in range(ug.GetNumberOfPoints()):
        lnodes={}
        p=ug.GetPoint(k)
        N=node_dict.get(p)
        if N==None:
            node_dict[p]=n
            n+=1
        lnodes[k]=node_dict.get(p)
    for k in range(ug.GetNumberOfCells()):
        cp=ug.GetCell(k).GetPoints()
        ele_vol_dict.append((j+1,
                              node_dict[cp.GetPoint(0)],
                              node_dict[cp.GetPoint(1)],
                              node_dict[cp.GetPoint(2)],
                              node_dict[cp.GetPoint(3)]))
        
rnode_dict={}
for k,v in node_dict.items():
    rnode_dict[v]=k
fname = fname[:-2]+'.msh'
file=open(fname,'w')

file.writelines(("$MeshFormat\n",
                 "2.2 0 8\n",
                 "$EndMeshFormat\n",
                 "$Nodes\n",
                "%d\n"%len(node_dict)))
for k in range(len(node_dict)):
    p=rnode_dict[k+1]
    file.write("%d %f %f %f\n"%(k+1,p[0],p[1],p[2]))
file.write("$EndNodes\n")
file.write("$Elements\n")
file.write("%d\n"%(len(ele_face_dict)+len(ele_vol_dict)))
for k,ele in enumerate(ele_face_dict):
    file.write("%d 2 2 %d %d %d %d %d\n"%(k+1,ele[0],ele[0],ele[1],ele[2],ele[3]))
for k,ele in enumerate(ele_vol_dict):
    file.write("%d 4 2 %d %d %d %d %d %d\n"%(k+1,ele[0],ele[0],ele[1],ele[2],ele[3],ele[4]))
file.write("$EndElements\n")
file.close()
print '...file created => '+ fname
