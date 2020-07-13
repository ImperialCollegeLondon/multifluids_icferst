/**********************************************************************/
/** Copyright (C) 2008,                                              **/
/** Queen Mary University of London (QMUL) & Imperial College        **/
/** of Science, Technology and Medicine (ICSTM). All rights reserved.**/
/** Implemented for you by Prof Antonio Munjiza & Dr Jiansheng Xiang **
 
 
 
* This code is part of the Virtual Geoscience Workbench (VGW) developed
* jointly by ICSTM and QMUL through two related parallel projects at 
* ICSTM and QMUL respectively funded by EPSRC. 
*
* This code is provided by copyright holders under the GNU Lesser 
* General Public License (LGPL). It is open source code; you can 
* redistribute it and/or modify it under the terms of the GNU Lesser 
* General Public License version 3.  
*  
* This code is distributed in the hope that it will be useful, 
* but WITHOUT ANY WARRANTY; without even the implied warranty 
* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
* the GNU Lesser General Public License for more details,
* http://www.gnu.org/licenses/lgpl.txt. 
*  
* You should have received a copy of the GNU Lesser General Public 
* License along with this code; if not, write to 
* Dr Jiansheng Xiang Prof Antonio Munjiza or Dr John-Paul Latham 
* j.xiang@imperial.ac.uk a.munjiza@qmul.ac.uk 
* or j.p.latham@imperial.ac.uk 
* ******************************************************************* */   

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkZLibDataCompressor.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkSystemIncludes.h>
#include "cubicsolver.h"

extern "C" {
#include "frame.h"
#include "Ytypes.h"
}


static void readppm(FILE* fdata,vtkPoints* pts, vtkCellArray* polys, 
		    vtkDoubleArray* vects, vtkDoubleArray* fcontact, vtkDoubleArray* ftotal,
                    vtkDoubleArray* fcsum, vtkDoubleArray* ftsum, vtkDoubleArray* vload,
		    vtkDoubleArray* sigma, vtkDoubleArray* tensor,
		    vtkDoubleArray* scalar1, vtkDoubleArray* scalar2, vtkDoubleArray* scalar3, vtkDoubleArray* scalar4,
		    INT npoint,
			 DBL dcsizc,DBL dcsizs,DBL dcsizf,DBL dcsizv)
{ 
  INT i1num[100],i;
  DBL d1num[100];
  DBL x[4],y[4],z[4],v[3],xc,yc,zc;
  DBL fc[3],ft[3];
  DBL fcs[3],fts[3];
  DBL vl[3],cals[3];
  DBL *v1,*v2,*v3,a1,a2,a3,tmp1,min,max,med;
  DBL **T;
  vtkIdType npts[4];
  CHR c1code[500000];
 
  T=TalDBL2(3,3);  

  v1=TalDBL1(3);
  v2=TalDBL1(3);
  v3=TalDBL1(3);

  a1=R0;
  a2=R0;
  a3=R0;

  for(i=0;i<3;i++)
    {
      v1[i]=R0;
      v2[i]=R0;
      v3[i]=R0;
    }

  tmp1=R0;
  min=R0;
  max=R0;
  med=R0;

  if(fdata!=FILENULL)
  { CHRr(fdata,c1code);
       CHRr(fdata,c1code);
    while((FILEND(fdata)==0))
    { codeCHRtoINT(c1code,i1num);
      codeINTtoDBL(d1num,i1num);
      if(i1num[2]==YTE3TET4ELS)
      { polys->InsertNextCell(4);
        for(i=0;i<4;i++)
	{ x[i]=d1num[3+i]*dcsizc;
	  y[i]=d1num[7+i]*dcsizc;
	  z[i]=d1num[11+i]*dcsizc;;
   
	  pts->InsertNextPoint(x[i],y[i],z[i]);
          npts[i]=npoint;
          polys->InsertCellPoint(npts[i]);
          npoint++;
	}
    
        T[0][0]=d1num[18]*dcsizs;
        T[0][1]=d1num[21]*dcsizs;
        T[1][0]=d1num[21]*dcsizs;
        T[1][1]=d1num[19]*dcsizs;
	T[0][2]=d1num[22]*dcsizs;
	T[2][0]=d1num[22]*dcsizs;
	T[1][2]=d1num[23]*dcsizs;
	T[2][1]=d1num[23]*dcsizs;
	T[2][2]=d1num[20]*dcsizs;
	       
        solveSymetricEigenProblem(T,v1,v2,v3,&a1,&a2,&a3);

	max=MAXIM(a1,a2);
	max=MAXIM(max,a3);
	min=MINIM(a1,a2);
	min=MINIM(min,a3);

	if((DABS(a1-max)>EPSILON)&&(DABS(a1-min)>EPSILON))
	  {
	    med=a1;
	  }
	if((DABS(a2-max)>EPSILON)&&(DABS(a2-min)>EPSILON))
	  {
	    med=a2;
	  }
	if((DABS(a3-max)>EPSILON)&&(DABS(a3-min)>EPSILON))
	  {
	    med=a3;
	  }

	tmp1=max-min;

	scalar1->InsertNextTuple1(max);
	scalar1->InsertNextTuple1(max);
	scalar1->InsertNextTuple1(max);
	scalar1->InsertNextTuple1(max);

	scalar2->InsertNextTuple1(med);
	scalar2->InsertNextTuple1(med);
	scalar2->InsertNextTuple1(med);
	scalar2->InsertNextTuple1(med);

	scalar3->InsertNextTuple1(min);
	scalar3->InsertNextTuple1(min);
	scalar3->InsertNextTuple1(min);
	scalar3->InsertNextTuple1(min);

	scalar4->InsertNextTuple1(tmp1);
	scalar4->InsertNextTuple1(tmp1);
	scalar4->InsertNextTuple1(tmp1);
	scalar4->InsertNextTuple1(tmp1);

        for(i=0;i<3;i++)
        { 
          v[i]=d1num[15+i]*dcsizv;
          fc[i]=d1num[24+i]*dcsizf;
          ft[i]=d1num[27+i]*dcsizf;
          fcs[i]=d1num[30+i]*dcsizc;
	  fts[i]=d1num[34+i];
        }
	
	vl[0]=R0;
	vl[1]=R0;
	vl[2]=R0;

	cals[0]=R0;
	cals[1]=d1num[33]*dcsizs;
	cals[2]=R0;
       
        vects->InsertNextTuple3(v[0],v[1],v[2]);
        vects->InsertNextTuple3(v[0],v[1],v[2]);
        vects->InsertNextTuple3(v[0],v[1],v[2]);
	vects->InsertNextTuple3(v[0],v[1],v[2]);

        fcontact->InsertNextTuple3(fc[0],fc[1],fc[2]);
        fcontact->InsertNextTuple3(fc[0],fc[1],fc[2]);
        fcontact->InsertNextTuple3(fc[0],fc[1],fc[2]);
	fcontact->InsertNextTuple3(fc[0],fc[1],fc[2]);

        ftotal->InsertNextTuple3(ft[0],ft[1],ft[2]);
        ftotal->InsertNextTuple3(ft[0],ft[1],ft[2]);
        ftotal->InsertNextTuple3(ft[0],ft[1],ft[2]);
	ftotal->InsertNextTuple3(ft[0],ft[1],ft[2]);

        fcsum->InsertNextTuple3(fcs[0],fcs[1],fcs[2]);
	fcsum->InsertNextTuple3(fcs[0],fcs[1],fcs[2]);
        fcsum->InsertNextTuple3(fcs[0],fcs[1],fcs[2]);
	fcsum->InsertNextTuple3(fcs[0],fcs[1],fcs[2]);

        ftsum->InsertNextTuple3(fts[0],fts[1],fts[2]);
        ftsum->InsertNextTuple3(fts[0],fts[1],fts[2]);
        ftsum->InsertNextTuple3(fts[0],fts[1],fts[2]);
	ftsum->InsertNextTuple3(fts[0],fts[1],fts[2]);

        vload->InsertNextTuple3(vl[0],vl[1],vl[2]);
        vload->InsertNextTuple3(vl[0],vl[1],vl[2]);
        vload->InsertNextTuple3(vl[0],vl[1],vl[2]);
	vload->InsertNextTuple3(vl[0],vl[1],vl[2]);

        sigma->InsertNextTuple3(cals[0],cals[1],cals[2]);
        sigma->InsertNextTuple3(cals[0],cals[1],cals[2]);
        sigma->InsertNextTuple3(cals[0],cals[1],cals[2]);
	sigma->InsertNextTuple3(cals[0],cals[1],cals[2]);

        tensor->InsertNextTuple9(T[0][0],T[0][1],T[0][2],
                                 T[1][0],T[1][1],T[1][2],
                                 T[2][0],T[2][1],T[2][2]); 
        tensor->InsertNextTuple9(T[0][0],T[0][1],T[0][2],
                                 T[1][0],T[1][1],T[1][2],
                                 T[2][0],T[2][1],T[2][2]); 
        tensor->InsertNextTuple9(T[0][0],T[0][1],T[0][2],
                                 T[1][0],T[1][1],T[1][2],
                                 T[2][0],T[2][1],T[2][2]); 
        tensor->InsertNextTuple9(T[0][0],T[0][1],T[0][2],
                                 T[1][0],T[1][1],T[1][2],
                                 T[2][0],T[2][1],T[2][2]); 
      }
     
      CHRr(fdata,c1code);

} } }

int main(int argc, char **argv)
{
  FILE *fdata=FILENULL;
  FILE *fcoeff=FILENULL;

  CHR name[300];
  CHR c1tmp[100],c1tmp1[100];
  CHR c1name[100];
  CHR c1name1[100];
  CHR c1name_vtk[100];
  DBL dcsizc,dcsizs,dcsizv,dcsizf;

  INT i,npoint,ini;
  INT start,finish;

  CHRw(stderr,"Start to convert files\n");
  CHRcpy(c1name1,argv[2]);
  fcoeff=fopen(c1name1,"r");
  if(fcoeff!=FILENULL)
  {
  CHRr(fcoeff,name);

  while(FILEND(fcoeff)==0) 
  { 
      if (CHRcmp(name,"/YD/YDC/DCSIZC",14)==0)
	  {DBLr(fcoeff,&dcsizc);
	  }
	  if (CHRcmp(name,"/YD/YDC/DCSIZF",14)==0)
	  {DBLr(fcoeff,&dcsizf);
	  }
	  if (CHRcmp(name,"/YD/YDC/DCSIZS",14)==0)
	  {DBLr(fcoeff,&dcsizs);
	  }
	  if (CHRcmp(name,"/YD/YDC/DCSIZV",14)==0)
	  {DBLr(fcoeff,&dcsizv);
	  }
      
    CHRr(fcoeff,name);
  }
  }
  else
  {
	  CHRw(stderr,"Could not open input file - usage -i inputfile"); 
	  CHRwcr(stderr);   
	  return 0;
  }

  start=atoi(argv[3]);
  finish=atoi(argv[4]);

  for(i=start;i<finish;i++)
  { vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
    vtkPoints *pts = vtkPoints::New();
    vtkCellArray *polys = vtkCellArray::New();

    vtkDoubleArray *scalar1 = vtkDoubleArray::New();
    scalar1->SetName("Sigma1");
    scalar1->SetNumberOfComponents(1);

    vtkDoubleArray *scalar2 = vtkDoubleArray::New();
    scalar2->SetName("Sigma2");
    scalar2->SetNumberOfComponents(1);

    vtkDoubleArray *scalar3 = vtkDoubleArray::New();
    scalar3->SetName("Sigma3");
    scalar3->SetNumberOfComponents(1);

    vtkDoubleArray *scalar4 = vtkDoubleArray::New();
    scalar4->SetName("DiffStress");
    scalar4->SetNumberOfComponents(1);

    vtkDoubleArray *vects = vtkDoubleArray::New();
    vects->SetName("VelocityVectors");
    vects->SetNumberOfComponents(3);

    vtkDoubleArray *fcontact = vtkDoubleArray::New();
    fcontact->SetName("ContactForces");
    fcontact->SetNumberOfComponents(3);

    vtkDoubleArray *ftotal = vtkDoubleArray::New();
    ftotal->SetName("TotalForces");
    ftotal->SetNumberOfComponents(3);

    vtkDoubleArray *fcsum = vtkDoubleArray::New();
    fcsum->SetName("Displacement");
    fcsum->SetNumberOfComponents(3);

    vtkDoubleArray *ftsum = vtkDoubleArray::New();
    ftsum->SetName("Strain");
    ftsum->SetNumberOfComponents(3);

    vtkDoubleArray *vload = vtkDoubleArray::New();
    vload->SetName("LoadingVelocity");
    vload->SetNumberOfComponents(3);

    vtkDoubleArray *sigma = vtkDoubleArray::New();
    sigma->SetName("CalStress");
    sigma->SetNumberOfComponents(3);

    vtkDoubleArray *tensor = vtkDoubleArray::New();
    tensor->SetName("STRESS");
    tensor->SetNumberOfComponents(9);

    dataSet->Allocate();

    vtkZLibDataCompressor* myZlibCompressor  = vtkZLibDataCompressor::New();
    myZlibCompressor->SetCompressionLevel(9);

    vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
    npoint=0; 

    CHRcpy(c1name,argv[1]);
    SINTw(c1tmp,i,0);
    CHRcat(c1name,c1tmp);
    CHRcpy(c1tmp,c1name);
    CHRcpy(c1tmp1,c1name);
    CHRcat(c1name,".ym");
    fdata=fopen(c1name,"r");
    if(fdata!=FILENULL)
      { readppm(fdata,pts,polys,vects,fcontact,ftotal,fcsum,ftsum,vload,sigma,tensor,scalar1,scalar2,scalar3,scalar4,npoint,dcsizc,dcsizs,dcsizf,dcsizv);
      fclose(fdata);
    }
    else break;
    dataSet->SetPoints(pts);
    pts->Delete();

    dataSet->GetPointData()->AddArray(vects);
    dataSet->GetPointData()->SetActiveAttribute("VELOCITY",vtkDataSetAttributes::VECTORS);
    dataSet->GetPointData()->SetVectors(vects);
    vects->Delete();

    dataSet->GetPointData()->AddArray(fcontact);
    dataSet->GetPointData()->SetActiveAttribute("CONTACT",vtkDataSetAttributes::VECTORS);
    dataSet->GetPointData()->SetVectors(fcontact);
    fcontact->Delete();

    dataSet->GetPointData()->AddArray(ftotal);
    dataSet->GetPointData()->SetActiveAttribute("TOTAL",vtkDataSetAttributes::VECTORS);
    dataSet->GetPointData()->SetVectors(ftotal);
    ftotal->Delete();

    dataSet->GetPointData()->AddArray(fcsum);
    dataSet->GetPointData()->SetActiveAttribute("CFSUM",vtkDataSetAttributes::VECTORS);
    dataSet->GetPointData()->SetVectors(fcsum);
    fcsum->Delete();

    dataSet->GetPointData()->AddArray(ftsum);
    dataSet->GetPointData()->SetActiveAttribute("TFSUM",vtkDataSetAttributes::VECTORS);
    dataSet->GetPointData()->SetVectors(ftsum);
    ftsum->Delete();

    dataSet->GetPointData()->AddArray(vload);
    dataSet->GetPointData()->SetActiveAttribute("VELLOAD",vtkDataSetAttributes::VECTORS);
    dataSet->GetPointData()->SetVectors(vload);
    vload->Delete();

    dataSet->GetPointData()->AddArray(sigma);
    dataSet->GetPointData()->SetActiveAttribute("CSTRESS",vtkDataSetAttributes::VECTORS);
    dataSet->GetPointData()->SetVectors(sigma);
    sigma->Delete();

    dataSet->GetPointData()->AddArray(tensor);
    dataSet->GetPointData()->SetActiveAttribute("STRESS",vtkDataSetAttributes::TENSORS);
    dataSet->GetPointData()->SetTensors(tensor);
    tensor->Delete();

    dataSet->GetPointData()->AddArray(scalar1);
    dataSet->GetPointData()->SetActiveAttribute("Sigma1",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(scalar1);
    scalar1->Delete();

    dataSet->GetPointData()->AddArray(scalar2);
    dataSet->GetPointData()->SetActiveAttribute("Sigma2",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(scalar2);
    scalar2->Delete();

    dataSet->GetPointData()->AddArray(scalar3);
    dataSet->GetPointData()->SetActiveAttribute("Sigma3",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(scalar3);
    scalar3->Delete();

    dataSet->GetPointData()->AddArray(scalar4);
    dataSet->GetPointData()->SetActiveAttribute("DiffStress",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(scalar4);
    scalar4->Delete();

    dataSet->SetCells(VTK_TETRA,polys);
    polys->Delete();
    CHRcpy(c1name_vtk,c1tmp1);
    CHRcat(c1name_vtk,".vtu");
    CHRw(stderr, "Writing vtu xml file ");
    CHRw(stderr, c1name_vtk);
    CHRwcr(stderr);
    writer->SetInput(dataSet);
    writer->SetFileName(c1name_vtk);
    writer->SetCompressor(myZlibCompressor);
    writer->Write();
    dataSet->Delete();
    myZlibCompressor->Delete();
    writer->Delete();
  }
}
