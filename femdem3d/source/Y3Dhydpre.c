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
* Prof A.Munjiza or Dr J.P.Latham 
* a.munjiza@qmul.ac.uk or j.p.latham@imperial.ac.uk 
* ******************************************************************* */ 

/*1 File   Yhpd.c */

#include "Yproto.h" 
/****************** PUBLIC ************************************/
void Yhpd(ydx,ydn,ydxn,d1hystp)

/*** integrate hydrostatic pressure over surface mesh  ***/
YDX ydx; YDN ydn; YDXN ydxn; DBL *d1hystp;
{ 
  DBL A[3],x[3],y[3],z[3],area_total;
  DBL xnc,ync,znc,norm,xc,yc,zc,nx,ny,nz;
  DBL xe,ye,ze,fx,fy,fz;
  INT i,j,k,n0,n1,n2,l,m,kc;
  area_total=R0; fx=R0; fy=R0; fz=R0;

  for (i=0; i<ydn->nnopo; i++) 
  {
      for(j=0;j<3;j++)
      {
          ydn->d2nfp[j][i]=R0;
      }
  }
  
  for (i=0; i<ydx->nsfel; i++) 
  {
      n0=ydx->i2xbto[i][0];
      n1=ydx->i2xbto[i][1];
      n2=ydx->i2xbto[i][2];

      x[0]=ydxn->d2xncc[0][n1]-ydxn->d2xncc[0][n0];
      y[0]=ydxn->d2xncc[1][n1]-ydxn->d2xncc[1][n0];
      z[0]=ydxn->d2xncc[2][n1]-ydxn->d2xncc[2][n0];
      x[2]=ydxn->d2xncc[0][n2]-ydxn->d2xncc[0][n0];
      y[2]=ydxn->d2xncc[1][n2]-ydxn->d2xncc[1][n0];
      z[2]=ydxn->d2xncc[2][n2]-ydxn->d2xncc[2][n0];
      x[1]=ydxn->d2xncc[0][n2]-ydxn->d2xncc[0][n1];
      y[1]=ydxn->d2xncc[1][n2]-ydxn->d2xncc[1][n1];
      z[1]=ydxn->d2xncc[2][n2]-ydxn->d2xncc[2][n1];
      V3DCro(xnc,ync,znc,x[0],y[0],z[0],x[2],y[2],z[2]);
      V3DNor(norm,xnc,ync,znc);
      ydx->d2normal[0][i]=xnc;
      ydx->d2normal[1][i]=ync;
      ydx->d2normal[2][i]=znc;
	
      xc=R0;yc=R0;zc=R0;
      for (j=0; j<3; j++) 
      {
          n0=ydx->i2xbto[i][j];
          xc=xc+ydxn->d2xncc[0][n0]/R3;
          yc=yc+ydxn->d2xncc[1][n0]/R3;
          zc=zc+ydxn->d2xncc[2][n0]/R3; 

 
     }
      for(j=0;j<3;j++)
      {
          n0=ydx->i2xbto[i][j];
          xe=xc-ydxn->d2xncc[0][n0];
          ye=yc-ydxn->d2xncc[1][n0];
          ze=zc-ydxn->d2xncc[2][n0];


          V3DCro(xnc,ync,znc,x[j],y[j],z[j],xe,ye,ze);
          V3DNor(norm,xnc,ync,znc);
          A[j]=norm/R2;			
      }
		
      if (ydx->nelno==6)
      {
          for (j=0; j<6; j++)
          {
              if(j<=2)
              {
                  k=j-1;
                  if (j==0)k=2;
                  ydx->d2area[j][i]=(A[j]+A[k])/R4;
              }
              if(j>2)
              {
                  k=j-3;
                  ydx->d2area[j][i]=A[k]/R2;	
              }
          }		
      }	
      else if (ydx->nelno==3)
      {
          for (j=0; j<3; j++)
          {
              k=j-1;
              if (j==0)k=2;
              ydx->d2area[j][i]=(A[j]+A[k])/R2;
          }			
      }		
      k=ydx->i2xbtojn[i][0];
      l=ydx->i2xbtojn[i][1];
      m=ydx->i2xbtojn[i][2];

      nx=ydx->d2normal[0][i];
      ny=ydx->d2normal[1][i];
      nz=ydx->d2normal[2][i];
		
      for (j=0; j<ydx->nelno; j++)
      {  



          // this is for fracture body
          k=ydx->i2xbtojn[i][j];
				kc=ydx->i2xbto[i][j];

				ydn->d2nfp[0][k]=ydn->d2nfp[0][k]+d1hystp[kc]*nx*ydx->d2area[j][i];
				ydn->d2nfp[1][k]=ydn->d2nfp[1][k]+d1hystp[kc]*ny*ydx->d2area[j][i];
				ydn->d2nfp[2][k]=ydn->d2nfp[2][k]+d1hystp[kc]*nz*ydx->d2area[j][i];


	  area_total=area_total+ydx->d2area[j][i];
          fx=fx+d1hystp[kc]*nx*ydx->d2area[j][i];
          fy=fy+d1hystp[kc]*ny*ydx->d2area[j][i];
          fz=fz+d1hystp[kc]*nz*ydx->d2area[j][i];

      }      		
  }

}

