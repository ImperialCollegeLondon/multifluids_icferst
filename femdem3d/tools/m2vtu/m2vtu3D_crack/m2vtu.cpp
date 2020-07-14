
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkZLibDataCompressor.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkSystemIncludes.h>

extern "C" {
#include "frame.h"
#include "Ytypes.h"
}


static void readppm(FILE* fdata,vtkPoints* pts, vtkCellArray* polys, 
		    vtkDoubleArray* vects, vtkDoubleArray* tensor, vtkDoubleArray* scalar,
		    INT npoint, DBL dcsizc,DBL dcsizs,DBL dcsizf,DBL dcsizv)
{ CHR c1code[5000000];
  INT i1num[100],i,j;
  DBL d1num[100];
  DBL x[4],y[4],z[4];
  DBL xg,yg;
  DBL Ncoef;   /* Ncoef is the Normalisation factor */

  DBL fr, fr_i;
    INT icount_t,icount_s;
    
    FILE *fda=FILENULL;
    fda=fopen("fract_ana.txt","a");
  
  vtkIdType npts[3];

  j=0;
    icount_t=0;
    icount_s=0;

  if(fdata!=FILENULL)
  { CHRr(fdata,c1code); CHRr(fdata,c1code);
    while((FILEND(fdata)==0))
    { codeCHRtoINT(c1code,i1num);
      codeINTtoDBL(d1num,i1num);
  //    if((i1num[2]==YTE3TET4JOINT)&&(d1num[21]<R0))
    if((i1num[2]==YTE3TET4JOINT))
      { for(i=0;i<3;i++)
	{ x[i]=d1num[3+i]*dcsizc;
	  y[i]=d1num[9+i]*dcsizc;
	  z[i]=d1num[15+i]*dcsizc;
	  pts->InsertNextPoint(x[i],y[i],z[i]);
          npts[i]=npoint;
          npoint++;
	}
	fr=d1num[21];
        fr_i=R0;
        if(fr<R0+0.01)
          {
           fr_i=-1.0;
           }
        else if(fr<0.1+0.01 && fr>0.1-0.01)
           {
             fr_i=1.0;
               icount_t=icount_t+1;
           }
        else if(fr<0.2+0.01 && fr>0.2-0.01)
           {
             fr_i=2.0;
                              icount_s=icount_s+1;
           }
        else if(fr<0.3+0.01 && fr>0.3-0.01)
           {
             fr_i=1.0;
                              icount_t=icount_t+1;
           }


//	DBLw(stdout,fr,10); CHRwsp(stdout);
        scalar->InsertNextTuple1(fr_i);
        scalar->InsertNextTuple1(fr_i);
        scalar->InsertNextTuple1(fr_i);
        polys->InsertNextCell(3,npts);

        for(i=0;i<3;i++)
        { x[i]=d1num[6+i]*dcsizc;
          y[i]=d1num[12+i]*dcsizc;
          z[i]=d1num[18+i]*dcsizc;
          pts->InsertNextPoint(x[i],y[i],z[i]);
          npts[i]=npoint;
          npoint++;
        }

        scalar->InsertNextTuple1(fr_i);
        scalar->InsertNextTuple1(fr_i);
        scalar->InsertNextTuple1(fr_i);
        polys->InsertNextCell(3,npts);

      }
      else if(i1num[2]==YTE2JOINTS)
      {// CHRw(stdout,"drawjoiint");
      }
      else if(i1num[2]==YTEG2RAD5)
      {// CHRw(stdout,"drawG2RP5");
      }
      CHRr(fdata,c1code);
      j=j+1;
} }


    INTw(fda,icount_t,10);
    CHRwsp(fda);
    INTw(fda,icount_s,10);
    CHRwcr(fda);
    fclose(fda);
}

main(INT argc, char **argv)
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

  //  ini=0;
  for(i=start;i<finish;i++)
  { vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
    vtkPoints *pts = vtkPoints::New();
    vtkCellArray *polys = vtkCellArray::New();

    vtkDoubleArray *vects = vtkDoubleArray::New();
    vects->SetName("VelocityVectors");
    vects->SetNumberOfComponents(3);

    vtkDoubleArray *fcontact = vtkDoubleArray::New();
    fcontact->SetName("ContactForces");
    fcontact->SetNumberOfComponents(3);

    vtkDoubleArray *tensor = vtkDoubleArray::New();
    tensor->SetName("STRESS");
    tensor->SetNumberOfComponents(9);

    vtkDoubleArray *scalar = vtkDoubleArray::New();
    scalar->SetName("Fracture");
    scalar->SetNumberOfComponents(1);

    dataSet->Allocate();

    vtkZLibDataCompressor* myZlibCompressor  = vtkZLibDataCompressor::New();
    myZlibCompressor->SetCompressionLevel(9);

    vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
//    vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();
  
    npoint=0; 

    CHRcpy(c1name,argv[1]);
    SINTw(c1tmp,i,0);
    CHRcat(c1name,c1tmp);
    CHRcpy(c1tmp,c1name);
    CHRcpy(c1tmp1,c1name);
    CHRcat(c1name,".ym");
    fdata=fopen(c1name,"r");
    if(fdata!=FILENULL)
      { readppm(fdata,pts,polys,vects,tensor,scalar,npoint,dcsizc,dcsizs,dcsizf,dcsizv);
      fclose(fdata);
    }
    else break;
    dataSet->SetPoints(pts);
    pts->Delete();

    dataSet->GetPointData()->AddArray(scalar);
    dataSet->GetPointData()->SetActiveAttribute("Fracture",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(scalar);
    scalar->Delete();

    dataSet->SetCells(VTK_TRIANGLE,polys);
    polys->Delete();

    CHRcpy(c1name_vtk,argv[1]);
    CHRcat(c1name_vtk,"_crack");
    SINTw(c1tmp,i,0);
    CHRcat(c1name_vtk,c1tmp);
    CHRcat(c1name_vtk,".vtu");
    CHRw(stderr, "Writing vtk xml file ");
    CHRw(stderr, c1name_vtk);
    CHRwcr(stderr);
    writer->SetInput(dataSet);
    writer->SetFileName(c1name_vtk);
//    writer->SetCompressor(myZlibCompressor);
    writer->Write();
    dataSet->Delete();
 //   myZlibCompressor->Delete();
    writer->Delete();
  }
}
