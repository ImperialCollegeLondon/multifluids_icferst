/*! \file Y3Drd.c
 *  \brief read model parameters
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


#include "Yproto.h"


static  CHR *cdig="0123456789";


/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/

/* random float between 0 and 1 */
double mrnds(r)   
  double *r;   
//  int n;   
  { int m;   
    double s,u,v,p;   
    s=65536.0; u=2053.0; v=13849.0;   
       *r=u*(*r)+v; m=(int)(*r/s);   
        *r=*r-m*s; p=*r/s;   
         
    return p;   
  }   

/* random integer between a and b  */
int mrabs(a,b,r)
  int a,b,*r;
  { int k,l,m,i;
    int p;
    k=b-a+1; l=2;
    while (l<k) l=l+l;
    m=4*l; k=*r; i=1;
	while(i<=1)
    { k=k+k+k+k+k;
	k=k%m; l=k/4+a;
        if (l<=b) { p=l;i=i+1;}
      }
    *r=k;
    return p;
  }



void Ypacking(YD yd)
{ 
//  YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
//  YDO ydo=&(yd->ydo);
  YDP ydp=&(yd->ydp);
  YDB ydb=&(yd->ydb);
  YPAR ypar=&(yd->ypar);
  YSP ysp=&(yd->ysp);


  DBL osp[3],max,disp[3];
  DBL d1u,d1v,d1w,d1cosine,d1sine;
  DBL zelta,d1ang[3],tmp;//,xmax,ymax,zmax;
  INT ipos,nnet,i,j,index,k,m,n,ind;
  const double mult[10]={1./6.,1./24.,1./24.,1./24.,1./60.,1./60.,1./60.,1./120.,1./120.,1./120.};
  DBL **inertia,cm[3],**points;
  DBL *intg;
    DBL *v1,*v2,*v3,ro;



  nnet=(yde->melem-yde->nelem-1)/ypar->nelem;
  if(nnet>=1)
  {
        v1=TalDBL1(3);
        v2=TalDBL1(3);
        v3=TalDBL1(3);
        inertia=TalDBL2(3,3);
        points=TalDBL2(4,3);
        intg=TalDBL1(10);

	  if(nnet>ypar->nlpar)nnet=ypar->nlpar;

  max=-BEPSILON;
  tmp=-BEPSILON;
/* maximum projection of existing packs    */
for(i=0;i<yde->nelem;i++)
{
	  if(yde->i1elbe[i]>1)
	  { 
	  for(j=0;j<3;j++)
	  {   tmp=ydn->d2ncc[0][yde->i2elto[i][j]]*ypar->d1lv[0]+
	          ydn->d2ncc[1][yde->i2elto[i][j]]*ypar->d1lv[1]+
			  ydn->d2ncc[2][yde->i2elto[i][j]]*ypar->d1lv[2];
		  if(tmp>max)max=tmp;
		  
	    }
    }
}
/* change coordinates according to position of exisiting packs  */
	  disp[0]=R0;
	  disp[1]=R0;
	  disp[2]=R0;
	  index=0;
if(max>(ypar->max-ypar->radius))
{
	disp[0]=(max-(ypar->max-ypar->radius))*ypar->d1lv[0];
	disp[1]=(max-(ypar->max-ypar->radius))*ypar->d1lv[1];
	disp[2]=(max-(ypar->max-ypar->radius))*ypar->d1lv[2];
	index=0;

	  for(i=0;i<ypar->ngrid;i++)
		  {
			  tmp=(ypar->d2gxyz[0][i]+ypar->radius+disp[0])*ypar->d1lv[0]+
               (ypar->d2gxyz[1][i]+ypar->radius+disp[1])*ypar->d1lv[1]+
               (ypar->d2gxyz[2][i]+ypar->radius+disp[2])*ypar->d1lv[2];
		  if(tmp>(ypar->d1dst[0]*ypar->d1lv[0]+
			  ypar->d1dst[1]*ypar->d1lv[1]+ypar->d1dst[2]*ypar->d1lv[2]))
		  {
			  index=1;
		      break;
		  }

	  
  }
}
/* place particles in random grid and rotate it in random angle */
if(index==0)
{	ydi->diedi=R2*ydi->diezon;
	for(i=0;i<(ypar->ngrid);i++)
	{
	ipos=mrabs(1,ypar->ngrid,&ypar->ir);
	if(ipos<=nnet)
	{
	osp[0]=(ypar->d2gxyz[0][i]+disp[0]);
	osp[1]=(ypar->d2gxyz[1][i]+disp[1]);
	osp[2]=(ypar->d2gxyz[2][i]+disp[2]);
	for(j=0;j<3;j++)
	{
	d1ang[j]=MYPI*mrnds(&ypar->r);
	}

	zelta=SQRT(d1ang[0]*d1ang[0]+d1ang[1]*d1ang[1]+d1ang[2]*d1ang[2]);
	if(zelta>EPSILON)
	{
	d1u=d1ang[0]/zelta;
	d1v=d1ang[1]/zelta;
	d1w=d1ang[2]/zelta;
	d1cosine=COS(zelta);
	d1sine=SIN(zelta);
	}
	

	/*node    */

	 for(k=0;k<ypar->nnode;k++)
	 {

	 /*
	V3DRot(&ydn->d2ncc[0][ydn->nnopo],&ydn->d2ncc[1][ydn->nnopo],&ydn->d2ncc[2][ydn->nnopo],
	ypar->d2xyz[0][k]-ypar->d1cg[0],ypar->d2xyz[1][k]-ypar->d1cg[1],ypar->d2xyz[2][k]-ypar->d1cg[2],
	d1u,d1v,d1w,d1cosine,d1sine);
	*/

		 
       d1u=ypar->d2gxyz[3][i];
       d1v=ypar->d2gxyz[4][i];
       d1w=ypar->d2gxyz[5][i];

    
	V3DRotAxes(&ydn->d2ncc[0][ydn->nnopo],&ydn->d2ncc[1][ydn->nnopo],&ydn->d2ncc[2][ydn->nnopo],
    ypar->d2xyz[0][k]-ypar->d1cg[0],ypar->d2xyz[1][k]-ypar->d1cg[1],ypar->d2xyz[2][k]-ypar->d1cg[2],d1u,d1v,d1w);

	
	ydn->d2ncc[0][ydn->nnopo]=ydn->d2ncc[0][ydn->nnopo]+osp[0];
	ydn->d2ncc[1][ydn->nnopo]=ydn->d2ncc[1][ydn->nnopo]+osp[1];
	ydn->d2ncc[2][ydn->nnopo]=ydn->d2ncc[2][ydn->nnopo]+osp[2];
	
    
	ydn->d2nci[0][ydn->nnopo]=ydn->d2ncc[0][ydn->nnopo];
	ydn->d2nci[1][ydn->nnopo]=ydn->d2ncc[1][ydn->nnopo];
	ydn->d2nci[2][ydn->nnopo]=ydn->d2ncc[2][ydn->nnopo];
	
	ydn->nnopo=ydn->nnopo+1;

	 }

	 /*nelem  */
	 ypar->nindex=ypar->nindex+1; //particle index
	 ysp->nspd=ysp->nspd+1;
 	 ydn->nnopo=ydn->nnopo-ypar->nnode;

	 for(k=0;k<ypar->nelem;k++)
	 {
		 for(m=0;m<(yde->nelno-1);m++)
		 {
			 yde->i2elto[yde->nelem][m]=ypar->i2elto[k][m]+ydn->nnopo;
			 ydn->i1nind[yde->i2elto[yde->nelem][m]]=ypar->nindex;
		 }
		 yde->i2elto[yde->nelem][yde->nelno-1]=ypar->nindex;
		 if(ypar->i1elbe[k]!=0)yde->i1elbe[yde->nelem]=ypar->nindex;
		 yde->nelem=yde->nelem+1;
	 }


	 ydn->nnopo=ydn->nnopo+ypar->nnode;

	 for(j=0;j<10;j++)intg[j]=R0;
	 
	 for(k=0;k<ypar->nelem;k++)
	 {	 
		 for(m=0;m<4;m++)
		 {
	     points[m][0]=ydn->d2ncc[0][yde->i2elto[yde->nelem-ypar->nelem+k][m]];
		 points[m][1]=ydn->d2ncc[1][yde->i2elto[yde->nelem-ypar->nelem+k][m]];
		 points[m][2]=ydn->d2ncc[2][yde->i2elto[yde->nelem-ypar->nelem+k][m]];
		 }
		 ro=ydp->d1pero[yde->i1elpr[yde->nelem-ypar->nelem+k]];
		 mirtichRoutine(points,intg,ro);
	 }

     for(n=0;n<10;n++)
      intg[n]=intg[n]*mult[n];

     //this is an error in the original publication!!!!!!!!!!!!!!!!!
     //they dont have the three!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 ind=ysp->nspd-1;
    ysp->d1mass[ind] = intg[0];  //density = 1.0 assumed 

     //center of mass 
    ysp->d2pp[0][ind]=intg[1]/ysp->d1mass[ind];
    ysp->d2pp[1][ind]=intg[2]/ysp->d1mass[ind];
    ysp->d2pp[2][ind]=intg[3]/ysp->d1mass[ind];
	cm[0]=ysp->d2pp[0][ind];
	cm[1]=ysp->d2pp[1][ind];
	cm[2]=ysp->d2pp[2][ind];

    
    inertia[0][0] = intg[5] + intg[6] -ysp->d1mass[ind]*(cm[1]*cm[1]+cm[2]*cm[2]);
    inertia[1][1] = intg[4] + intg[6] -ysp->d1mass[ind]*(cm[0]*cm[0]+cm[2]*cm[2]);
    inertia[2][2] = intg[4] + intg[5] -ysp->d1mass[ind]*(cm[1]*cm[1]+cm[0]*cm[0]);

    inertia[0][1] = -(intg[7] -ysp->d1mass[ind]*cm[1]*cm[0]);
    inertia[1][0]=inertia[0][1];

    inertia[1][2] = -(intg[8] -ysp->d1mass[ind]*cm[1]*cm[2]);
    inertia[2][1]=inertia[1][2];

    inertia[0][2] = -(intg[9] -ysp->d1mass[ind]*cm[2]*cm[0]);
    inertia[2][0]=inertia[0][2];

//    inertia*=((REAL)(1.0/3.0));

	for(m=0;m<3;m++)
	  {

		  for(n=0;n<3;n++)
	  {
       inertia[m][n]=inertia[m][n]*(R1/R3);
	  }
	  }

    //this is because the vertices could have been DoubleLinkedListed in the wrong direction
    if((inertia[0][0]<0.0)||(inertia[1][1]<0.0)||(inertia[2][2]<0.0))
      {
      //inertia*=(-1.0);
	  for(m=0;m<3;m++)
	  {
		  for(n=0;n<3;n++)
	  {
       inertia[m][n]=-inertia[m][n];
	  }
	  }
       ysp->d1mass[ind]*=(-1.0);
      }


	solveSymetricEigenProblem(inertia, v1,v2,v3,
						&ysp->d2prine[0][ind],&ysp->d2prine[1][ind],&ysp->d2prine[2][ind]);
		for(m=0;m<3;m++)
	  {
		ysp->d2prinn1[m][ind]=v1[m];
		ysp->d2prinn2[m][ind]=v2[m];
		ysp->d2prinn3[m][ind]=v3[m];

		ysp->d2pvt[m][ind]=ypar->d1vc[m];
		}

	}
	}
}

    FREE(v1);
	FREE(v2);
	FREE(v3);
	FREE(inertia);
	FREE(points);
	FREE(intg);
}
	return;
}



/* run packing test   */


