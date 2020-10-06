import vtk
import numpy

def AddOrExtract(P,pdict,*args):
    pt=[]
    for arg in args:
        pt.append(P.GetPoint(arg))
    p=tuple(numpy.mean(numpy.array(pt),axis=0))

    if p in pdict:
        return pdict[p]
    else:
        pdict[p]=P.InsertNextPoint(p)
        return pdict[p]

    return p

def MakeHex(*IdList):
    T=vtk.vtkHexahedron()
    for k,id in enumerate(IdList):
        T.GetPointIds().SetId(k,id)

    return T

def MakeQuad(*IdList):
    T=vtk.vtkQuad()
    for k,id in enumerate(IdList):
        T.GetPointIds().SetId(k,id)

    return T

def vtu2cv(BaseName,FieldNames):

    infile=vtk.vtkXMLUnstructuredGridReader()
    infile.SetFileName(BaseName+'.vtu')
    infile.Update()

    Vout=vtk.vtkUnstructuredGrid()
    Vout.Allocate(1,1)
    P=vtk.vtkPoints()
    Vout.SetPoints(P)

    Nele=infile.GetNumberOfCells()
    Npts=infile.GetNumberOfPoints()
    VIn=infile.GetOutput()
    Ids=vtk.vtkIdList()
    pts=VIn.GetPoints()
    
    CellData=vtk.vtkCellData()
    
    for name in FieldNames:
        data=vtk.vtkDoubleArray()
        data.SetNumberOfValues(4*Nele)
        data.SetName(name)
        Vout.GetCellData().AddArray(data)

    point_dict={}
    cell_list=[]

    for k in range(Npts):
        p=pts.GetPoint(k)
        point_dict[p]=P.InsertNextPoint(p)

    for k in range(Nele):
        type=VIn.GetCellType(k)
        VIn.GetCellPoints(k,Ids)

        TT=[]

        if type==vtk.VTK_TRIANGLE:
            p1=Ids.GetId(0)
            p2=Ids.GetId(1)
            p3=Ids.GetId(2)

	    p4=AddOrExtract(P,point_dict,p1,p2)
            p5=AddOrExtract(P,point_dict,p1,p3)
            p6=AddOrExtract(P,point_dict,p2,p3)
	    
	    p7=AddOrExtract(P,point_dict,p1,p2,p3)

	    TT.append(MakeQuad(p1,p4,p7,p5))	
	    TT.append(MakeQuad(p2,p6,p7,p4))	
	    TT.append(MakeQuad(p3,p5,p7,p6))	

	    PD=VIn.GetPointData()		
            for kk,name in enumerate(FieldNames):
                PD.SetActiveScalars(name)
                CD=Vout.GetCellData().GetArray(kk)
                for kkk,p in enumerate((p1,p2,p3)):
                    CD.SetValue(3*k+kkk,
                                PD.GetScalars().GetValue(p))

        elif type==vtk.VTK_TETRA:
            p1=Ids.GetId(0)
            p2=Ids.GetId(1)
            p3=Ids.GetId(2)
            p4=Ids.GetId(3)

            p5=AddOrExtract(P,point_dict,p1,p2)
            p6=AddOrExtract(P,point_dict,p1,p3)
            p7=AddOrExtract(P,point_dict,p1,p4)
            p8=AddOrExtract(P,point_dict,p2,p3)
            p9=AddOrExtract(P,point_dict,p2,p4)
            p10=AddOrExtract(P,point_dict,p3,p4)

            p11=AddOrExtract(P,point_dict,p1,p2,p3)
            p12=AddOrExtract(P,point_dict,p1,p2,p4)
            p13=AddOrExtract(P,point_dict,p1,p3,p4)
            p14=AddOrExtract(P,point_dict,p2,p3,p4)

            p15=AddOrExtract(P,point_dict,p1,p2,p3,p4)
            

            TT.append(MakeHex(p1,p5,p11,p6,p7,p12,p15,p13))
            TT.append(MakeHex(p2,p8,p11,p5,p9,p14,p15,p12))
            TT.append(MakeHex(p3,p6,p11,p8,p10,p13,p15,p14))
            TT.append(MakeHex(p4,p7,p12,p9,p10,p13,p15,p14))


            PD=VIn.GetPointData()
            for kk,name in enumerate(FieldNames):
                PD.SetActiveScalars(name)
                CD=Vout.GetCellData().GetArray(kk)
                for kkk,p in enumerate((p1,p2,p3,p4)):
                    CD.SetValue(4*k+kkk,
                                PD.GetScalars().GetValue(p))

        for T in TT:
            Vout.InsertNextCell(T.GetCellType(),T.GetPointIds())
            cell_list.append(T.GetPointId(0))

    outfile=vtk.vtkXMLUnstructuredGridWriter()
    outfile.SetFileName('cv_'+BaseName+'.vtu')
    outfile.SetInputData(Vout)
    outfile.Write()

    return VIn, Vout


if __name__=="__main__":
    import sys
    import glob

    for fname in glob.glob(sys.argv[1]):
        fname,type=fname.rsplit('.',1)
        if type in ('vtu',) and not fname[0:2]=='cv_':
            vtu2cv(fname,sys.argv[2:])
