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

static void Y3DBNPRESSURE(
			  INT nelem, INT *i1elpr, INT **i2elto,
			  INT **i2elbnp,
			  DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
			  INT *i1nopr, INT nelemi
			  )
{
  INT ielem,inopo,jnopo,knopo;
  INT i,j,k;
  INT noprx,nopry,noprz;

  noprx=9;
  nopry=1;
  noprz=2;

  //noprx=0;
  //nopry=0;
  //noprz=0;

  for(ielem=0;ielem<nelem;ielem++)
    {
      if(i1elpr[ielem]==0)
	{
	  for(i=0;i<4;i++)
	    {
	      j=i+1;
	      if(j>3)
		{
		  j=0;
		}
	      k=j+1;
	      if(k>3)
		{
		  k=0;
		}

	      inopo=i2elto[ielem][i];
	      jnopo=i2elto[ielem][j];
	      knopo=i2elto[ielem][k];

	      if((i1nopr[inopo]==16)&&(i1nopr[jnopo]==16)&&(i1nopr[knopo]==16))
		{
		  i2elbnp[ielem][i+j+k-3]=16;
		}
	      else if((DABS(d1nccx[inopo]-0.001)<1.0e-4)&&(DABS(d1nccx[jnopo]-0.001)<1.0e-4)&&(DABS(d1nccx[knopo]-0.001)<1.0e-4)) //x=0 surface
		{
		  i2elbnp[ielem][i+j+k-3]=noprx;
		}
	      //else if((DABS(d1nccx[inopo]-0.375)<1.0e-4)&&(DABS(d1nccx[jnopo]-0.375)<1.0e-4)&&(DABS(d1nccx[knopo]-0.375)<1.0e-4)) //x=1.5 surface
	      //{
	      //i2elbnp[ielem][i+j+k-3]=noprx;
	      //}
	      else if((DABS(d1nccy[inopo])<1.0e-4)&&(DABS(d1nccy[jnopo])<1.0e-4)&&(DABS(d1nccy[knopo])<1.0e-4)) //y=0 surface
		{
		  i2elbnp[ielem][i+j+k-3]=nopry;
		}
	      else if((DABS(d1nccy[inopo]-0.1)<1.0e-4)&&(DABS(d1nccy[jnopo]-0.1)<1.0e-4)&&(DABS(d1nccy[knopo]-0.1)<1.0e-4)) //y=0 surface
		{
		  i2elbnp[ielem][i+j+k-3]=nopry;
		}
	      else if((DABS(d1nccz[inopo]-0.05)<1.0e-4)&&(DABS(d1nccz[jnopo]-0.05)<1.0e-4)&&(DABS(d1nccz[knopo]-0.05)<1.0e-4)) //z=0 surface
		{
		  i2elbnp[ielem][i+j+k-3]=noprz;
		}
	      //else if((DABS(d1nccz[inopo]+0.25)<1.0e-4)&&(DABS(d1nccz[jnopo]+0.25)<1.0e-4)&&(DABS(d1nccz[knopo]+0.25)<1.0e-4)) //z=0 surface
	      //{
	      //i2elbnp[ielem][i+j+k-3]=noprz;
	      //}
	    }
	}
    }
/*
    for(ielem=0;ielem<nelem;ielem++)
    {
      if(i1elpr[ielem]==0)
	{
	  for(i=0;i<4;i++)
	    {
	      inopo=i2elto[ielem][i];

	      if((DABS(d1nccz[inopo])<1.0e-4)&&(i1nopr[inopo]==16))   //z=0 surface
		{
		  i1nopr[inopo]=13;
		}
	    }
	}
    }
*/
}

void Ybnp(YDE yde, YDN ydn)
{
  Y3DBNPRESSURE(
		yde->nelem,yde->i1elpr,yde->i2elto,
		yde->i2elbnp,
		ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2ncc[2],
		ydn->i1nopr, yde->nelemi
		);
}
