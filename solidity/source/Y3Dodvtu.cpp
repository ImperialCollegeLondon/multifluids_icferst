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

/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/


/* small strain softening quadratic tetrahedron output
*/
static void Yod3TET10(  
            INT nelem,   INT nnopo, 
            CHR *fout,
            DBL dcsizc,  DBL dcsizs,   DBL dcsizv, DBL dcsizf,
            INT icoutp,  INT iprop, 
            DBL **d2ncc, DBL **d2nvc,  DBL **d2nft,
            INT *i1elpr, INT *i1nobf, INT **i2elto, DBL ***d3tcs,INT *i1ptyp, DBL *d1ntc
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
void Yodvtu(CHR *namep, YD yd)
{ 
  YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
  YDO ydo=&(yd->ydo);
  YDP ydp=&(yd->ydp);
  YDX ydx=&(yd->ydx);
  YDXN ydxn=&(yd->ydxn);  
  INT iprop,ihys,i;
  CHR namef[300];
  CHR name_avg_temp[13];
  CHR cindex[50];
  DBL tmp;
  FILE *fout=FILENULL;

  /* output animation */
  if((ydc->ncstep%ydc->icoutf)==0)   
  { 
//  	INTw(stdout,ydc->ncstep/ydc->icoutf,10);CHRw(stdout,"/");INTw(stdout,ydc->mcstep/ydc->icoutf,5);CHRwcr(stdout);
    CHRcpynoext(namef,namep);
    SINTw(cindex,ydc->icouti,0);
    CHRcat(namef,"_continnum_");
    CHRcat(namef,cindex);
    CHRcat(namef,".vtu");

    printf("IN Yodvtu \n");

      iprop=0;
	Yod3TET10(  
          ydx->ncmel, ydxn->ncmno,
          namef,
          ydc->dcsizc, ydc->dcsizs, ydc->dcsizv, ydc->dcsizf, 
          ydc->icoutp, iprop,
          ydn->d2ncc,  ydn->d2nvc,  ydn->d2nft,
          yde->i1elpr, ydn->i1nobf, ydx->i2cvto, yde->d3tcs,ydp->i1ptyp, ydn->d1ntc
       );
  }
}


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/



