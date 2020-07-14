/*! \file Y3Dod.c
 *  \brief output database for Y3D model
 *
 *  \todo
 *  finish documentation
 *
 *  Copyright (C) 2008, Queen Mary University of London (QMUL) & 
 *  Imperial College of Science, Technology and Medicine (ICSTM).
 *  All rights reserved. Implemented by Prof Antonio Munjiza & 
 *  Dr Jiansheng Xiang.
 *
 *  This code is part of the Virtual Geoscience Workbench (VGW) 
 *  developed jointly by ICSTM and QMUL through two related parallel 
 *  projects at ICSTM and QMUL respectively funded by EPSRC. 
 *
 *  This code is provided by copyright holders under the GNU Lesser 
 *  General Public License (LGPL). It is open source code; you can 
 *  redistribute it and/or modify it under the terms of the GNU Lesser 
 *  General Public License version 3.  
 *  
 *  This code is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty 
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
 *  the GNU Lesser General Public License for more details,
 *  http://www.gnu.org/licenses/lgpl-3.0.txt. 
 *
 *  You should have received a copy of the GNU Lesser General Public 
 *  License along with this code; if not, write to:
 *
 *  Dr Jiansheng Xiang   <j.xiang@imperial.ac.uk>    \n
 *  Prof Antonio Munjiza <a.munjiza@qmul.ac.uk>      \n
 *  Dr John-Paul Latham  <j.p.latham@imperial.ac.uk> \n
 *
 */


/**********************************************************************/
/**********************************************************************/


#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkZLibDataCompressor.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkSystemIncludes.h>
#include "petscksp.h"
//
extern "C" {
#include "Yproto.h"
#include "frame.h"
}
static INT  i1num[100];    /* numbers for space saving format     */
static DBL  d1num[100];    /* numbers for space saving format     */
static CHR c1code[500];    /* coded i1para in space saving format */

/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/


/* small strain softening quadratic tetrahedron output
*/
static void Yod3TET10_T(  
            INT nelem,   INT nnopo, 
            CHR *fout,
            DBL dcsizc,  DBL dcsizs,   DBL dcsizv, DBL dcsizf,
            INT icoutp,  INT iprop, 
            DBL **d2ncc, DBL **d2nvc,  DBL **d2nft,
            INT *i1elpr, INT *i1nobf, INT **i2elto, DBL ***d3tcs,INT *i1ptyp, DBL *d1ntc, DBL dctime
            ) 
{ 
  INT  i1num[100];    /*Z moved - numbers for space saving format     */
  DBL  d1num[100];    /*Z moved - numbers for space saving format     */
  CHR c1code[500];    /*Z moved - coded i1para in space saving format */
  INT ielem, i, j, k, ip;
  DMX *dm1tn;
  DBL **d2tcs, dtnn;
  DBL tmp,tmp1,tmp2;  
  INT *i1elto, *i1tnn;

  dm1tn = TalzDMX1(nnopo); 
  i1tnn = TalzINT1(nnopo);

  for(ielem=0; ielem<nelem; ielem++)
  {
    i1elto = i2elto[ielem];
    d2tcs = d3tcs[ielem];

	if(i1ptyp[i1elpr[ielem]]==YTE3TET4ELS)
	  {
    for(i=0; i<NNODEX; i++)
    {
      ip = i1elto[i];

      for(j=0; j<NDIME; j++)
	    {
	      for(k=0; k<NDIME; k++)
	      {
          dm1tn[ip][j][k] += d2tcs[j][k];
	      }
      }

      i1tnn[ip]++;
    }

	}
	else
	{
    for(i=0; i<NNODE; i++)
    {
      ip = i1elto[i];

      for(j=0; j<NDIME; j++)
	    {
	      for(k=0; k<NDIME; k++)
	      {
          dm1tn[ip][j][k] += d2tcs[j][k];
	      }
      }

      i1tnn[ip]++;
    }

	}

  }

    vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
    vtkPoints *pts = vtkPoints::New();
    vtkCellArray *polys = vtkCellArray::New();
    
    vtkDoubleArray *vects = vtkDoubleArray::New();
    vects->SetName("VelocityVectors");
    vects->SetNumberOfComponents(3);
    
    vtkDoubleArray *fcontact = vtkDoubleArray::New();
    fcontact->SetName("NodalForces");
    fcontact->SetNumberOfComponents(3);
    	
    vtkDoubleArray *tensor = vtkDoubleArray::New();
    tensor->SetName("STRESS");
    tensor->SetNumberOfComponents(9);
 
    vtkDoubleArray *temp = vtkDoubleArray::New();
    temp->SetName("Temperature");
    temp->SetNumberOfComponents(1);

    vtkDoubleArray *index = vtkDoubleArray::New();
    index->SetName("Index");
    index->SetNumberOfComponents(1);

    vtkDoubleArray *time = vtkDoubleArray::New();
    time->SetName("Time");
    time->SetNumberOfComponents(1);
    
    vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
    //    vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();
    
    vtkZLibDataCompressor* myZlibCompressor  = vtkZLibDataCompressor::New();
    myZlibCompressor->SetCompressionLevel(1);
    
    dataSet->Allocate();

	for(i=0; i<nnopo; i++)
	{
	dtnn = (DBL)i1tnn[i];

	for(j=0; j<NDIME; j++)
	{
		for(k=0; k<NDIME; k++)
		{
			dm1tn[i][j][k] /= dtnn;
		}
	}

	tmp=d2ncc[0][i];
	tmp1=d2ncc[1][i];
	tmp2=d2ncc[2][i];
	pts->InsertNextPoint(tmp,tmp1,tmp2);

	tmp=d2nvc[0][i];
	tmp1=d2nvc[1][i];
	tmp2=d2nvc[2][i];
	vects->InsertNextTuple3(tmp,tmp1,tmp2);

	tmp=d2nft[0][i];
	tmp1=d2nft[1][i];
	tmp2=d2nft[2][i];
	fcontact->InsertNextTuple3(tmp,tmp1,tmp2);

	tensor->InsertNextTuple9(dm1tn[i][0][0],dm1tn[i][0][1],dm1tn[i][0][2],
	                       dm1tn[i][0][1],dm1tn[i][1][1],dm1tn[i][1][2],
	                       dm1tn[i][0][2],dm1tn[i][1][2],dm1tn[i][2][2]);

	tmp=d1ntc[i];
	temp->InsertNextTuple1(tmp);
	tmp=i1nobf[i];
	index->InsertNextTuple1(tmp);

	tmp=dctime;
	time->InsertNextTuple1(tmp);
}

      for(ielem=0; ielem<nelem; ielem++)
      {
          i1elto = i2elto[ielem];
          polys->InsertNextCell(4);
          for(i=0;i<4;i++)
          {
              polys->InsertCellPoint(i1elto[i]);
          }
      }
      dataSet->SetPoints(pts);
      pts->Delete();
      
      dataSet->GetPointData()->AddArray(vects);
      dataSet->GetPointData()->SetActiveAttribute("VELOCITY",vtkDataSetAttributes::VECTORS);
      dataSet->GetPointData()->SetVectors(vects);
      vects->Delete();
      
      dataSet->GetPointData()->AddArray(fcontact);
      dataSet->GetPointData()->SetActiveAttribute("NODALFORCES",vtkDataSetAttributes::VECTORS);
      dataSet->GetPointData()->SetVectors(fcontact);
      fcontact->Delete();
 
      dataSet->GetPointData()->AddArray(tensor);
      dataSet->GetPointData()->SetActiveAttribute("STRESS",vtkDataSetAttributes::TENSORS);
      dataSet->GetPointData()->SetTensors(tensor);
      tensor->Delete();

      dataSet->GetPointData()->AddArray(temp);
      dataSet->GetPointData()->SetActiveAttribute("TEMP",vtkDataSetAttributes::SCALARS);
      dataSet->GetPointData()->SetScalars(temp);
      temp->Delete();

      dataSet->GetPointData()->AddArray(index);
      dataSet->GetPointData()->SetActiveAttribute("INDEX",vtkDataSetAttributes::SCALARS);
      dataSet->GetPointData()->SetScalars(index);
      index->Delete();

      dataSet->GetPointData()->AddArray(time);
      dataSet->GetPointData()->SetActiveAttribute("TIME",vtkDataSetAttributes::SCALARS);
      dataSet->GetPointData()->SetScalars(index);
      time->Delete();

      dataSet->SetCells(VTK_TETRA,polys);
      
      polys->Delete();
      writer->GetInput();
      writer->SetFileName(fout);
      writer->SetCompressor(myZlibCompressor);
      writer->Write();
      dataSet->Delete();
      myZlibCompressor->Delete();
      writer->Delete();


  FREE(i1tnn);
  FREE(dm1tn);
}

/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/


/* small strain softening quadratic tetrahedron output */
/*
static void Yod3TET10(  
     INT nelem,   INT nnopo, 
     FILE *fout,  FILE *fvtk,
     DBL dcsizc,  DBL dcsizs,   DBL dcsizv, DBL dcsizf,
     INT icoutp,  INT iprop, 
     DBL **d2ncc, DBL **d2nvc,  DBL **d2nfc,
     INT *i1elpr, INT **i2elto, DBL ***d3tcs,INT *i1ptyp
     ) 
{ 
}
*/

/*   4-node tetrahedra element output   */

static void Yod3TET4ELS(  
     INT nelem, 
     FILE *fout,
     DBL dcsizc,  DBL dcsizs,   DBL dcsizv, DBL dcsizf,
     INT icoutp,  INT iprop, 
     DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
//     DBL *d1ncix, DBL *d1nciy, DBL *d1nciz,
     DBL *d1nvcx, DBL *d1nvcy, DBL *d1nvcz,
     DBL *d1nfcx, DBL *d1nfcy, DBL *d1nfcz,
     INT *i1elpr, INT **i2elto, DBL ***d3tcs
     ) 
{ 
     INT ielem, i;
     DBL **d2tcs;
     INT *i1elto;
//     FILE *fp;

//     fp=fopen("Y3Dod_output.txt","w+");

     for(ielem=0; ielem<nelem; ielem++)
     {
	  if(i1elpr[ielem]==iprop)
	  {
	       i1elto=i2elto[ielem];
	       d2tcs=d3tcs[ielem];
	       
	       /* prepare output */
	       i1num[0]=icoutp;
	       i1num[1]=27;
	       //d1num[2]=R0; /*Z filler, overwrite */
	       d1num[3]=d1nccx[i1elto[0]]/dcsizc;
	       d1num[4]=d1nccx[i1elto[1]]/dcsizc;
	       d1num[5]=d1nccx[i1elto[2]]/dcsizc;
	       d1num[6]=d1nccx[i1elto[3]]/dcsizc;
	       d1num[7]=d1nccy[i1elto[0]]/dcsizc;
	       d1num[8]=d1nccy[i1elto[1]]/dcsizc;
	       d1num[9]=d1nccy[i1elto[2]]/dcsizc;
	       d1num[10]=d1nccy[i1elto[3]]/dcsizc;
	       d1num[11]=d1nccz[i1elto[0]]/dcsizc;
	       d1num[12]=d1nccz[i1elto[1]]/dcsizc;
	       d1num[13]=d1nccz[i1elto[2]]/dcsizc;
	       d1num[14]=d1nccz[i1elto[3]]/dcsizc;

//             fprintf(fp,"ielem=%ld\nNodeNo=%ld\t%ld\t%ld\t%ld\n",ielem,i1elto[0],i1elto[1],i1elto[2],i1elto[3]);

//	       fprintf(fp,"ielem=%ld\nX=%f\t%f\t%f\t%f\n",ielem,d1ncix[i1elto[0]],d1ncix[i1elto[1]],
//		       d1ncix[i1elto[2]],d1ncix[i1elto[3]]);
//	       fprintf(fp,"ielem=%ld\nY=%f\t%f\t%f\t%f\n",ielem,d1nciy[i1elto[0]],d1nciy[i1elto[1]],
//		       d1nciy[i1elto[2]],d1nciy[i1elto[3]]);
//	       fprintf(fp,"ielem=%ld\nZ=%f\t%f\t%f\t%f\n",ielem,d1nciz[i1elto[0]],d1nciz[i1elto[1]],
//		       d1nciz[i1elto[2]],d1nciz[i1elto[3]]);
	       
	       d1num[15]=(d1nvcx[i1elto[0]]+d1nvcx[i1elto[1]]+d1nvcx[i1elto[2]]+d1nvcx[i1elto[3]])
		    /(R4*dcsizv);
	       d1num[16]=(d1nvcy[i1elto[0]]+d1nvcy[i1elto[1]]+d1nvcy[i1elto[2]]+d1nvcy[i1elto[3]])
		    /(R4*dcsizv);
	       d1num[17]=(d1nvcz[i1elto[0]]+d1nvcz[i1elto[1]]+d1nvcz[i1elto[2]]+d1nvcz[i1elto[3]])
		    /(R4*dcsizv);
	       d1num[18]=d2tcs[0][0]/dcsizs;
	       d1num[19]=d2tcs[1][1]/dcsizs;
	       d1num[20]=d2tcs[2][2]/dcsizs;
	       d1num[21]=d2tcs[0][1]/dcsizs;
	       d1num[22]=d2tcs[0][2]/dcsizs;
	       d1num[23]=d2tcs[1][2]/dcsizs;

	       d1num[24]=(d1nfcx[i1elto[0]]+d1nfcx[i1elto[1]]+d1nfcx[i1elto[2]]+d1nfcx[i1elto[3]])/4.0/dcsizf;
	       d1num[25]=(d1nfcy[i1elto[0]]+d1nfcy[i1elto[1]]+d1nfcy[i1elto[2]]+d1nfcy[i1elto[3]])/4.0/dcsizf;
	       d1num[26]=(d1nfcz[i1elto[0]]+d1nfcz[i1elto[1]]+d1nfcz[i1elto[2]]+d1nfcz[i1elto[3]])/4.0/dcsizf;
	       d1num[27]=R0;
/*
	       d1num[25]=d1nfcx[i1elto[0]]/dcsizf;
	       d1num[26]=d1nfcx[i1elto[1]]/dcsizf;
	       d1num[27]=d1nfcx[i1elto[2]]/dcsizf;
	       d1num[28]=d1nfcx[i1elto[3]]/dcsizf;
	       d1num[29]=d1nfcy[i1elto[0]]/dcsizf;
	       d1num[30]=d1nfcy[i1elto[1]]/dcsizf;
	       d1num[31]=d1nfcy[i1elto[2]]/dcsizf;
	       d1num[32]=d1nfcy[i1elto[3]]/dcsizf;
	       d1num[33]=d1nfcz[i1elto[0]]/dcsizf;
	       d1num[34]=d1nfcz[i1elto[1]]/dcsizf;
	       d1num[35]=d1nfcz[i1elto[2]]/dcsizf;
	       d1num[36]=d1nfcz[i1elto[3]]/dcsizf;*/
	       //d1num[37]=R0; /* elastic damage */

	       for(i=3; i<27; i++)
	       { 
		    d1num[i]=MAXIM((-R1),MINIM(d1num[i],R1));
	       }
	       
	       /* translate into INT */
	       codeDBLtoINT(d1num, i1num); 
	       i1num[2] = YTE3TET4ELS; /* overwrite filler */
	       codeINTtoCHR(c1code, i1num);
	       CHRw(fout, c1code);
	       CHRwcr(fout);
	  }
     }

//     fclose(fp);
}



static void Yod3TET4JOINTSINTACT(  /* 3D joint output for 4-node tetrahedra */
     INT nelem,
     FILE *fout,
     DBL dcsizc,DBL dcsizv,
     DBL dpefs, DBL dpeft, DBL dpegf, DBL dpeks, DBL dpepe,
     INT icoutp,INT iprop ,
     DBL *d1nccx,DBL *d1nccy,DBL *d1nccz,
     DBL *d1nvcx,DBL *d1nvcy,DBL *d1nvcz,
     DBL *d1sdel,
     INT *i1elpr,INT **i2elto
     )
{
}

static void Yod3TET4JOINTSBROKEN(  /* 3D joint output for 4-node tetrahedra */
     INT nelem,
     FILE *fout,
     DBL dcsizc, DBL dcsizv,
     INT icoutp, INT iprop ,
     DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
     DBL *d1nvcx, DBL *d1nvcy, DBL *d1nvcz,
     DBL *d1sdel,
     INT *i1elpr, INT **i2elto, INT *i1elty
     )
{    INT ielem;
     INT i;
     INT ipropc;
     INT *i1elto;
     
     for(ielem=0;ielem<nelem;ielem++)
     {
          i1elto=i2elto[ielem]; 
	  ipropc=i1elpr[ielem];
	  if(ipropc==iprop)
	  {
	       if((i1elto[0]==i1elto[5])&&(i1elto[1]==i1elto[4])&&(i1elto[2]==i1elto[3]))
	       {
		    ipropc=ipropc-YIPROPMAX*2;
	       }
	  }
	  
	  if((ipropc<0)&&(((ipropc+YIPROPMAX)==iprop)||((ipropc+YIPROPMAX*2)==iprop)||((ipropc+YIPROPMAX*2)==0)||(i1elty[ielem]>0)))
	  {
               /* prepare output */
	       i1num[0]=icoutp;
	       i1num[1]=22;
	       d1num[3]=d1nccx[i1elto[0]]/dcsizc;
	       d1num[4]=d1nccx[i1elto[1]]/dcsizc;
	       d1num[5]=d1nccx[i1elto[2]]/dcsizc;
	       d1num[6]=d1nccx[i1elto[3]]/dcsizc;
	       d1num[7]=d1nccx[i1elto[4]]/dcsizc;
	       d1num[8]=d1nccx[i1elto[5]]/dcsizc;
	       d1num[9]=d1nccy[i1elto[0]]/dcsizc;
	       d1num[10]=d1nccy[i1elto[1]]/dcsizc;
	       d1num[11]=d1nccy[i1elto[2]]/dcsizc;
	       d1num[12]=d1nccy[i1elto[3]]/dcsizc;
	       d1num[13]=d1nccy[i1elto[4]]/dcsizc;
	       d1num[14]=d1nccy[i1elto[5]]/dcsizc;
	       d1num[15]=d1nccz[i1elto[0]]/dcsizc;
	       d1num[16]=d1nccz[i1elto[1]]/dcsizc;
	       d1num[17]=d1nccz[i1elto[2]]/dcsizc;
	       d1num[18]=d1nccz[i1elto[3]]/dcsizc;
	       d1num[19]=d1nccz[i1elto[4]]/dcsizc;
	       d1num[20]=d1nccz[i1elto[5]]/dcsizc;
	       
	       if(((ipropc+YIPROPMAX)==iprop)||((ipropc+YIPROPMAX*2)==0)) /* fracture surface */
		 {
		   d1num[21]=-0.2;
		 }
	       else		/* boundary surface */
		 {
		   d1num[21]=-0.1;
		 }

               if(i1elty[ielem]==2)              // pre-ex
                 {
                   d1num[21]=0.2;
                 }
            	else if(i1elty[ielem]==3)              // new-normal
                 {
                   d1num[21]=0.3;
                 }
               else if(i1elty[ielem]==4)              // new-shear
                 {
                   d1num[21]=0.4;
                 }
               else if(i1elty[ielem]==5)              // new-mixed rotational
                 {
                   d1num[21]=0.5;
                 }
	//       d1num[21]=i1elty[ielem]/10.0;


//	       if(d1sdel!=DBL1NULL)
//		    d1num[21]=d1sdel[ielem];
	       for(i=3;i<22;i++)
	       {
		    d1num[i]=MAXIM((-R1),MINIM(d1num[i],R1));
	       }
	       /* translate into INT */
	       codeDBLtoINT(d1num,i1num);
	       i1num[2]=YTE3TET4JOINT;
	       codeINTtoCHR(c1code,i1num);
	       CHRw(fout,c1code);
	       CHRwcr(fout);
	  }
     }
}

//-PY changed it for 3D_fracture_coupling_with_multiphase
static void Yod3gid(  
            INT nelem,   INT nnopo,INT nspd, 
            FILE *fout,  FILE *fvtk,
            DBL dcsizc,  DBL dcsizs,   DBL dcsizv, DBL dcsizf,
            INT icoutp,  INT iprop, 
            DBL **d2ncc, DBL **d2nvc,  DBL **d2nfc,
            INT *i1elpr, INT **i2elto, INT *i1ptyp,
			DBL *d1svrx,DBL *d1svry,DBL *d1svrz,
            DBL *d1p0x,DBL *d1p0y,DBL *d1p0z,DBL *d1p1x,DBL *d1p1y,DBL *d1p1z,
            DBL *d1p2x, DBL *d1p2y,DBL *d1p2z,DBL *cenx,DBL *ceny,DBL *cenz,
			DBL *d1svx,DBL *d1svy,DBL *d1svz,INT *i1nind,
		    DBL *d1pax, DBL *d1pay, DBL *d1paz
            ) 
{ 
  INT  i1num[100];    /*Z moved - numbers for space saving format     */
  DBL  d1num[100];    /*Z moved - numbers for space saving format     */
  CHR c1code[500];    /*Z moved - coded i1para in space saving format */
  INT ielem, i, j, k, ip;
  INT *i1elto, ia;
  DBL tmpx,tmpy,tmpz,dxc,dyc,dzc,vrxc,vryc,vrzc;

  for(i=0; i<nnopo; i++)
  {
      ia = i1nind[i] - 1;

      dxc = d2ncc[0][i] - cenx[ia];
      dyc = d2ncc[1][i] - ceny[ia];
      dzc = d2ncc[2][i] - cenz[ia];

      tmpx = d1svrx[ia] * (d1p0y[ia] * dzc - d1p0z[ia] * dyc)
	  + d1svry[ia] * (d1p1y[ia] * dzc - d1p1z[ia] * dyc)
	  + d1svrz[ia] * (d1p2y[ia] * dzc - d1p2z[ia] * dyc);
      tmpy = d1svrx[ia] * (d1p0z[ia] * dxc - d1p0x[ia] * dzc)
	  + d1svry[ia] * (d1p1z[ia] * dxc - d1p1x[ia] * dzc)
	  + d1svrz[ia] * (d1p2z[ia] * dxc - d1p2x[ia] * dzc);
      tmpz = d1svrx[ia] * (d1p0x[ia] * dyc - d1p0y[ia] * dxc)
	  + d1svry[ia] * (d1p1x[ia] * dyc - d1p1y[ia] * dxc)
	  + d1svrz[ia] * (d1p2x[ia] * dyc - d1p2y[ia] * dxc);

      vrxc = tmpx + d1svx[ia];
      vryc = tmpy + d1svy[ia];
      vrzc = tmpz + d1svz[ia];

      /* prepare output */
      i1num[0] = icoutp;
      i1num[1] = 19;
      //INTw(stdout,ia,5);
      //CHRwcr(stdout);
      d1num[2] = R0; /*Z filler, overwrite */
      d1num[3] = d2ncc[0][i] / dcsizc;
      d1num[4] = d2ncc[1][i] / dcsizc;
      d1num[5] = d2ncc[2][i] / dcsizc;
      d1num[6] = vrxc / dcsizv;
      d1num[7] = vryc / dcsizv;
      d1num[8] = vrzc / dcsizv;
      d1num[9] = R0;
      d1num[10] = R0;
      d1num[11] = R0;
      d1num[12] = R0;
      d1num[13] = R0;
      d1num[14] = R0;
      d1num[15] = d2nfc[0][i] / dcsizf;
      d1num[16] = d2nfc[1][i] / dcsizf;
      d1num[17] = d2nfc[2][i] / dcsizf;
      d1num[18] = R0; /* elastic damage */
      d1num[19] = R0;
    for(j=3; j<17; j++)
    { 
      /*Z faster than either ternary or if/elseif */
      if(d1num[j] < -R1) d1num[j] = -R1;
      if(d1num[j] >  R1) d1num[j] =  R1;
    }
    /* translate into INT */
    codeDBLtoINT(d1num, i1num); 
    i1num[2] = YTE3TET10ELS; /* overwrite filler */
//    i1num[1]=ia;
    codeINTtoCHR(c1code, i1num);
    CHRw(fout, c1code); CHRwcr(fout);
  }
	  
	for(ielem=0; ielem<nelem; ielem++)
	{
    i1elto = i2elto[ielem];

    i1num[ 0] = icoutp;
	  i1num[ 1] = 15;
  	  i1num[2]=i1ptyp[i1elpr[ielem]]+1;
	  //i1num[ 2] = YTE3TETELS1;
          i1num[ 3] = i1elto[0];
          i1num[ 4] = i1elto[1];
          i1num[ 5] = i1elto[2];
          i1num[ 6] = i1elto[3];
          i1num[ 7] = i1elto[4];

          if(i1ptyp[i1elpr[ielem]]==YTE3TET4ELS)
          {
          i1num[ 8] = 0;
          i1num[ 9] = 0;
          i1num[10] = 0;
          i1num[11] = 0;
          i1num[12] = 0;
          i1num[13] = 0;
          }
          else
          {
          i1num[ 7] = i1elto[4];
          i1num[ 8] = i1elto[5];
          i1num[ 9] = i1elto[6];
          i1num[10] = i1elto[7];
          i1num[11] = i1elto[8];
          i1num[12] = i1elto[9];
          i1num[13] = i1elto[10];
          }
          i1num[14] = 0;

    codeINTtoCHR(c1code, i1num);
    CHRw(fout, c1code); CHRwcr(fout);
	}
//INTw(stdout,nspd,5);
//CHRwcr(stdout);
  /*
   for(k=0; k<nspd; k++)
   {


   i1num[ 0] = icoutp;
   i1num[ 1] = 25;
   d1num[ 2] = R0;
   d1num[ 3] = cenx[k]/dcsizc;
   d1num[ 4] = ceny[k]/dcsizc;
   d1num[ 5] = cenz[k]/dcsizc;
   d1num[ 6] = fmod(d1pax[k],R2*MYPI)/(R2*MYPI);
   d1num[ 7] = fmod(d1pay[k],R2*MYPI)/(R2*MYPI);
   d1num[ 8] = fmod(d1paz[k],R2*MYPI)/(R2*MYPI);

   d1num[9] = R0;
   codeDBLtoINT(d1num, i1num);
   i1num[2]=25;
   codeINTtoCHR(c1code, i1num);
   CHRw(fout, c1code); CHRwcr(fout);
   }
   */
  
}


/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/


/*! \brief output the Y database
 *  \param[in] namep base file name
 *  \param[in] yd pointer to Y database structure
 *  \par Details:
 *  Yod() writes to the history file every 100 time steps, and outputs frames 
 *  according to the model output frequency.
 *
 *  Z - needs polished
 *
 *  \see #YDC_struct.icoutf
 */
 void Yod(CHR *namep, YD yd)
 { 
      YDC ydc=&(yd->ydc);
      YDE yde=&(yd->yde);
      YDI ydi=&(yd->ydi);
      YDN ydn=&(yd->ydn);
      YDO ydo=&(yd->ydo);
      YDP ydp=&(yd->ydp);
YDX ydx=&(yd->ydx);
  YDXN ydxn=&(yd->ydxn);    

//-PY changed it for 3D_fracture_coupling_with_multiphase
//      YSP ysp=&(yd->ysp);

      INT iprop,ihys,i;
      CHR namef[300];
      CHR cindex[50];
      static INT ncall=0;
      static INT ndefinitions=0;
      FILE *fout=FILENULL;
      FILE *fragm=FILENULL;
//-PY changed it for 3D_fracture_coupling_with_multiphase
//      FILE *fvtk=FILENULL;

      DBL tmp,tmass;

      ncall=ncall+1;

      /* output history */
      if((ydc->ncstep%100)==0)
      {
	   for(ihys=0;ihys<(ydo->nohys);ihys++)
	   { 
		ydo->d1ohyt[ihys]=ydc->dctime;
		if((ydo->d1ohyt[ihys])==(ydc->dctime))
		{ 
		     if((ydo->f2ohyf[ihys])==FILENULL)
		     { 
			  //CHRcpy(namef,namep);
			  //SINTw(cindex,ihys,0);
			  //CHRcat(namef,"h");
			  //CHRcat(namef,cindex); 
			  //ydo->f2ohyf[ihys]=fopen(namef,"a");
		     }
		     if((ydo->f2ohyf[ihys])!=FILENULL)
		     { 
			  tmp=(ydo->d1ohyt[ihys])*(ydo->d1ohyc[ihys]);
			  DBLw((ydo->f2ohyf[ihys]),tmp,17); 
			  CHRwsp(ydo->f2ohyf[ihys]);
			  tmp=(ydo->d1ohys[ihys])*(ydo->d1ohyf[ihys]);
			  DBLw((ydo->f2ohyf[ihys]),tmp,17);
			  CHRwcr(ydo->f2ohyf[ihys]);       
		     } 
		} 
	   }
      }

      //CHRw(stdout,"OK");
      //CHRwcr(stdout);

      /*Z note - odd to hard-code this */
      if((ncall>100)||((ydc->ncstep)>=(ydc->mcstep-2)))
      { 
	   ncall=0;
	   for(ihys=0;ihys<(ydo->nohys);ihys++)
	   { 
		if((ydo->f2ohyf[ihys])!=FILENULL)
		{
		     fclose(ydo->f2ohyf[ihys]);
		}
		ydo->f2ohyf[ihys]=FILENULL;
	   } 
      }

      //CHRw(stdout,"OK");
      //CHRwcr(stdout);

      /* output animation */
      fout=FILENULL;
      if((ydc->ncstep%ydc->icoutf)==0)   
      { 
		if(ydc->iflag[2]>0)
		{
			CHRcpynoext(namef,namep);
			SINTw(cindex,ydc->icouti,0);
			CHRcat(namef,"_continnum_");
			CHRcat(namef,cindex);
			CHRcat(namef,".vtu");

			iprop=0;
			Yod3TET10_T(  
			      ydx->ncmel, ydxn->ncmno,
			      namef,
			      ydc->dcsizc, ydc->dcsizs, ydc->dcsizv, ydc->dcsizf, 
			      ydc->icoutp, iprop,
			      ydn->d2ncc,  ydn->d2nvc,  ydn->d2nft,
			      yde->i1elpr, ydn->i1nobf, ydx->i2cvto, yde->d3tcs,ydp->i1ptyp, ydxn->d1xntc,ydc->dctime);
		}
	   
	   CHRcpynoext(namef,namep);
	   SINTw(cindex,ydc->icouti,0);
	   CHRcat(namef,cindex);
	   CHRcat(namef,".ym");
	   fout=fopen(namef,"w");
          
       fragm=fopen("fragment.txt","a");
          tmass=R0;
          for(i=0;i<yde->nelem;i++)
          {
              if(ydp->i1ptyp[yde->i1elpr[i]]==(YTE3TET4ELS))
              {
                  if((ydn->d2ncc[1][yde->i2elto[i][0]]>0.0205)||
                     (ydn->d2ncc[1][yde->i2elto[i][1]]>0.0205)||
                     (ydn->d2ncc[1][yde->i2elto[i][2]]>0.0205)||
                     (ydn->d2ncc[1][yde->i2elto[i][3]]>0.0205)) tmass=tmass+yde->d1emct[i];
                     
              }
              
          }
          DBLw(fragm,ydc->dctime,15);
          CHRwsp(fragm);
          DBLw(fragm,tmass,15);
          CHRwcr(fragm);
          fclose(fragm);
          
          
	   //CHRwcr(stdout);

	   if(fout!=FILENULL)
	   { 
		CHRw(fout,"CODED");
		CHRwcr(fout);

		for(iprop=0;iprop<ydp->nprop;iprop++)
		{
		     if((ydp->i1ptyp[iprop])==(YTE3TET4ELS))
		     {
			  Yod3TET4ELS(   /*   4-node tetrahedra element output   */  
			       yde->nelem,
			       fout,
			       ydc->dcsizc, ydc->dcsizs, ydc->dcsizv, ydc->dcsizf, 
			       ydc->icoutp, iprop,
                               ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2ncc[2],
			       //ydn->d2nci[0],ydn->d2nci[1],ydn->d2nci[2],
			       ydn->d2nvc[0],ydn->d2nvc[1],ydn->d2nvc[2],
			       ydn->d2nfp[0],ydn->d2nfp[1],ydn->d2nfp[2],
			       yde->i1elpr, yde->i2elto, yde->d3tcs
			       );
		     }

		     else if((ydp->i1ptyp[iprop])==(YTE3TET4JOINT))
		     {
		       /* 3D joint output for 4-node tetrahedra */
		       /*
    			  Yod3TET4JOINTSINTACT(  
		     	       yde->nelem,
			       fout,
			       ydc->dcsizc,ydc->dcsizv,
			       ydp->d1pefs[iprop],ydp->d1peft[iprop],ydp->d1pegf[iprop],
			       ydp->d1peks[iprop],ydp->d1pepe[iprop],
			       ydc->icoutp, iprop,
			       ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2ncc[2],
			       ydn->d2nvc[0],ydn->d2nvc[1],ydn->d2nvc[2],
			       yde->d2elst[ydp->i1psde[iprop]],
			       yde->i1elpr,yde->i2elto
			       );
		       */		    
			  Yod3TET4JOINTSBROKEN(
			       yde->nelem,
			       fout,
			       ydc->dcsizc,ydc->dcsizv,
			       ydc->icoutp, iprop,
			       ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2ncc[2],
			       ydn->d2nvc[0],ydn->d2nvc[1],ydn->d2nvc[2],
			       yde->d2elst[ydp->i1psde[iprop]],
			       yde->i1elpr,yde->i2elto, yde->i1elty
			       );
		     }
		}


//-PY changed it for 3D_fracture_coupling_with_multiphase //PROBLEM HERE
/*
		  Yod3gid(  
			yde->nelem, ydn->nnopo,ysp->nspd,
			fout, fvtk,
			ydc->dcsizc, ydc->dcsizs, ydc->dcsizv, ydc->dcsizf, 
			ydc->icoutp, iprop,
			ydn->d2ncc,  ydn->d2nvc,  ydn->d2nfc,
			yde->i1elpr, yde->i2elto, ydp->i1ptyp,
			ysp->d2pvr[0],ysp->d2pvr[1],ysp->d2pvr[2],
			ysp->d2prinn1[0],ysp->d2prinn1[1],ysp->d2prinn1[2],
			ysp->d2prinn2[0],ysp->d2prinn2[1],ysp->d2prinn2[2],
			ysp->d2prinn3[0],ysp->d2prinn3[1],ysp->d2prinn3[2],
			ysp->d2pp[0],ysp->d2pp[1],ysp->d2pp[2],
			ysp->d2pvt[0],ysp->d2pvt[1],ysp->d2pvt[2],ydn->i1nind,
			ysp->d2pa[0],ysp->d2pa[1],ysp->d2pa[2]);
*/

      ydc->icouti=ydc->icouti+1;
      fclose(fout);

      /*Z progress */
      CHRw(stderr,"frame ");
      INTw(stderr,ydc->icouti,0);
      CHRw(stderr,"/");
      INTw(stderr,ydc->mcstep/ydc->icoutf,0);
      CHRwcr(stderr);
     
	   }
   }
}


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/



