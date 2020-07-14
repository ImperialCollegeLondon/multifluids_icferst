/*! \file Y3Dfd.c^M
 *  \brief Y nodal forces^M
 *^M
 *  \todo^M
 *  finish documentation^M
 *^M
 *  Copyright (C) 2008, Queen Mary University of London (QMUL) & ^M
 *  Imperial College of Science, Technology and Medicine (ICSTM).^M
 *  All rights reserved. Implemented by Prof Antonio Munjiza & ^M
 *  Dr Jiansheng Xiang.^M
 *^M
 *  This code is part of the Virtual Geoscience Workbench (VGW) 
 *  developed jointly by ICSTM and QMUL through two related parallel 
 *  projects at ICSTM and QMUL respectively funded by EPSRC. 
 *
 *  This code is provided by copyright holders under the GNU Lesser ^M
 *  General Public License (LGPL). It is open source code; you can 
 *  redistribute it and/or modify it under the terms of the GNU Lesser 
 *  General Public License version 3.  
 *  
 *  This code is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty 
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
 *  the GNU Lesser General Public License for more details,
 *  http://www.gnu.org/licenses/lgpl-3.0.txt. 
 *^M
 *  You should have received a copy of the GNU Lesser General Public ^M
 *  License along with this code; if not, write to:^M
 *^M
 *  Dr Jiansheng Xiang   <j.xiang@imperial.ac.uk>    \n^M
 *  Prof Antonio Munjiza <a.munjiza@qmul.ac.uk>      \n^M
 *  Dr John-Paul Latham  <j.p.latham@imperial.ac.uk> \n^M
 *^M
 */


/**********************************************************************/
/**********************************************************************/


#include "Yproto.h"

static void Y3DPLACEPLATEN(
			   INT nelem, INT nnopo,
			   INT *i1elpr, INT **i2elto,
			   DBL *d1nciy, DBL *d1nccy,
			   INT *i1nopr, INT *i1bnvx, INT *i1bnvy, INT *i1bnvz, DBL *d1bnap
			   )
{
  INT upp,low,specimen;
  DBL ymax,ymin,uppy,lowy;
  INT ielem,inopo;
  INT i;
  INT *i1elto;

  upp=4;
  low=5;
  specimen=0;
  ymax=0.0;
  ymin=0.0;

  for(ielem=0;ielem<nelem;ielem++)
    {
      if(i1elpr[ielem]==specimen)
	{
	  i1elto=i2elto[ielem];
	  for(i=0;i<4;i++)
	    {
	      if(d1nccy[i1elto[i]]>ymax)
		{
		  ymax=d1nccy[i1elto[i]];
		}
	      if(d1nccy[i1elto[i]]<ymin)
		{
		  ymin=d1nccy[i1elto[i]];
		}
	    }
	}
    }

  uppy=0.0195-ymax;
  lowy=-0.0195-ymin;

  for(inopo=0;inopo<nnopo;inopo++)
    {
      if(i1nopr[inopo]==upp)   //upper platen
	{
	  d1nciy[inopo]=d1nciy[inopo]-uppy;
	  d1nccy[inopo]=d1nccy[inopo]-uppy;
	}
      if(i1nopr[inopo]==low)   //lower platen
	{
	  d1nciy[inopo]=d1nciy[inopo]-lowy;
	  d1nccy[inopo]=d1nccy[inopo]-lowy;
	}
    }

  for(i=1;i<4;i++)
    {
      i1bnvx[i]=0;
      i1bnvy[i]=0;
      i1bnvz[i]=0;
    }
  for(i=4;i<6;i++)
    {
      i1bnvx[i]=1;
      i1bnvy[i]=1;
      i1bnvz[i]=1;
    }

  //d1bnap[2]=0.0;
}

void Yreset(YDE yde, YDN ydn, YDB ydb)
{
  Y3DPLACEPLATEN(
		 yde->nelem,ydn->nnopo,
		 yde->i1elpr,yde->i2elto,
		 ydn->d2nci[1],ydn->d2ncc[1],
		 ydn->i1nopr,ydb->i1bnvx,ydb->i1bnvy,ydb->i1bnvz,ydb->d1bnap
		 );
}
