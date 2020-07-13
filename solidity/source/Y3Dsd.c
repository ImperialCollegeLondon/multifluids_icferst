/*! \file Y3Dsd.c
 *  \brief solver database for Y3D model
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
 *  Dr Jiansheng Xiang   <j.xiang@imperial.ac.uk>    \nd1bnfx
 *  Prof Antonio Munjiza <a.munjiza@qmul.ac.uk>      \n
 *  Dr John-Paul Latham  <j.p.latham@imperial.ac.uk> \n
 *
 */


/**********************************************************************/
/**********************************************************************/

#include "Yproto.h" 
//-PY changed it for 3D_fracture_coupling_with_multiphase
static INLINE void nodevel(
						   INT nelem,   INT nnopo,INT nspd,
						   DBL **d2ncc, DBL **d2nvc,  DBL **d2nfc,
						   INT *i1elpr, INT **i2elto, INT *i1ptyp,
						   DBL *d1svrx,DBL *d1svry,DBL *d1svrz,
						   DBL *d1p0x,DBL *d1p0y,DBL *d1p0z,DBL *d1p1x,DBL *d1p1y,DBL *d1p1z,
						   DBL *d1p2x, DBL *d1p2y,DBL *d1p2z,DBL *cenx,DBL *ceny,DBL *cenz,
						   DBL *d1svx,DBL *d1svy,DBL *d1svz,INT *i1nind,
						   DBL *d1pax, DBL *d1pay, DBL *d1paz
						   )
{
	INT ielem, i, j, k, ip;
	INT *i1elto, ia;
	DBL tmpx,tmpy,tmpz,dxc,dyc,dzc,vrxc,vryc,vrzc;
	/*
	 i1num=TalINT1(100);
	 d1num=TalDBL1(100);
	 c1code=TalCHR1(500);  
	 */

//     FILE *fn;
//     fn=fopen("nodevel_test.txt","w+");

	for(i=0; i<nnopo; i++)
	{
		ia=i1nind[i]-1;
		
		dxc=d2ncc[0][i]-cenx[ia];
		dyc=d2ncc[1][i]-ceny[ia];
		dzc=d2ncc[2][i]-cenz[ia];
		
        tmpx=d1svrx[ia]*(d1p0y[ia]*dzc-d1p0z[ia]*dyc)
		+d1svry[ia]*(d1p1y[ia]*dzc-d1p1z[ia]*dyc)
		+d1svrz[ia]*(d1p2y[ia]*dzc-d1p2z[ia]*dyc);
        tmpy=d1svrx[ia]*(d1p0z[ia]*dxc-d1p0x[ia]*dzc)
		+d1svry[ia]*(d1p1z[ia]*dxc-d1p1x[ia]*dzc)
		+d1svrz[ia]*(d1p2z[ia]*dxc-d1p2x[ia]*dzc);
        tmpz=d1svrx[ia]*(d1p0x[ia]*dyc-d1p0y[ia]*dxc)
		+d1svry[ia]*(d1p1x[ia]*dyc-d1p1y[ia]*dxc)
		+d1svrz[ia]*(d1p2x[ia]*dyc-d1p2y[ia]*dxc);
		
        d2nvc[0][i]=tmpx+d1svx[ia];
        d2nvc[1][i]=tmpy+d1svy[ia];
        d2nvc[2][i]=tmpz+d1svz[ia];

//fprintf(fn, "d2nvc[0][i], d2nvc[1][i], d2nvc[2][i]   %lf %lf %lf \n",  d2nvc[0][i], d2nvc[1][i], d2nvc[2][i]);


	}
//fclose(fn);
}

/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/

//-PY changed it for 3D_fracture_coupling_with_multiphase
int sign(DBL v)
{
return v > R0 ? 1 : (v < R0 ? -1 : R0);
}



/*Z calculate current nodal force and change in velocity
*/
static INLINE void nfvel(YDN ydn,
     DBL *d1bnfx, DBL *d1bnfy, DBL *d1bnfz, 
     DBL *d1nfcx, DBL *d1nfcy, DBL *d1nfcz, 
     DBL *d1nftx, DBL *d1nfty, DBL *d1nftz,  
     //DBL *d1kax,  DBL *d1kay,  DBL *d1kaz,
     INT *i1bnvx, INT *i1bnvy, INT *i1bnvz,
     DPT *dp1dv,  
     INT nnopo, 
     DBL *d1nmct, INT *i1nopr,  
	 DBL *d1nfvx,  DBL *d1nfvy,  DBL *d1nfvz,
	 DBL *d1nfpx,  DBL *d1nfpy,  DBL *d1nfpz,
	 DBL *d1nvcx,  DBL *d1nvcy,  DBL *d1nvcz,
	 DBL *d1nvct,  DBL *absorp,
	 DBL *denf,     INT index,
	 DBL *f1dt, DBL *f1dt1, DBL *f1dt2
	)
{
     INT inopo, iprop;
	DBL nmct, f1bf[3],factor;


//     FILE *fp;
//     fp=fopen("nfvel_output.txt","w+");
	if(index==2)
	{
		f1dt[0]=R0; f1dt1[0]=R0; f1dt2[0]=R0;
		f1dt[1]=R0; f1dt1[1]=R0; f1dt2[1]=R0;
		f1dt[2]=R0; f1dt1[2]=R0; f1dt2[2]=R0;
		f1dt[3]=R0; f1bf[0]=R0; f1bf[1]=R0; f1bf[2]=R0; 
	}

     for(inopo=0; inopo<nnopo; inopo++)
     {  
	  iprop = i1nopr[inopo];
	  nmct  = d1nmct[inopo];
//factor=0.25;  
factor=0.35;
/*
		 d1nftx[inopo]=d1nftx[inopo]+factor*(d1nfvx[inopo]-d1nvcx[inopo]*absorp[inopo])*d1nvct[inopo]
												+d1nfpx[inopo];
		 d1nfty[inopo]=d1nfty[inopo]+factor*(d1nfvy[inopo]-d1nvcy[inopo]*absorp[inopo])*d1nvct[inopo]
												+d1nfpy[inopo];
		 d1nftz[inopo]=d1nftz[inopo]+factor*(d1nfvz[inopo]-d1nvcz[inopo]*absorp[inopo])*d1nvct[inopo]
												+d1nfpz[inopo];
*/
//		 fprintf(fp, "********************************* \n" );
//		 fprintf(fp,"inopo %ld %f %f %f %f %f\n",
//			 inopo,d1nvcx[inopo],d1nvcy[inopo],d1nvcz[inopo],absorp[inopo],d1nvct[inopo]);
//		 fprintf(fp,"Fx[%ld]=%f\tFy[%ld]=%f\tFz[%ld]=%f\n",
//			 inopo,d1nftx[inopo],inopo,d1nfty[inopo],inopo,d1nftz[inopo]);
//-ao getting rid of absorb


		 d1nftx[inopo]=d1nftx[inopo]+(d1nfvx[inopo]-d1nvcx[inopo])*factor*absorp[inopo]*d1nvct[inopo]+d1nfpx[inopo]+ d1nfcx[inopo];;
		 d1nfty[inopo]=d1nfty[inopo]+(d1nfvy[inopo]-d1nvcy[inopo])*factor*absorp[inopo]*d1nvct[inopo]+d1nfpy[inopo]+ d1nfcy[inopo];;
		 d1nftz[inopo]=d1nftz[inopo]+(d1nfvz[inopo]-d1nvcz[inopo])*factor*absorp[inopo]*d1nvct[inopo]+d1nfpz[inopo]+ d1nfcz[inopo];;

//-PY changed it for 3D solid-fluid coupling
/*
                 d1nftx[inopo]=d1nftx[inopo]+d1nfpx[inopo] + d1nfcx[inopo];
		 d1nfty[inopo]=d1nfty[inopo]+d1nfpy[inopo] + d1nfcy[inopo];
		 d1nftz[inopo]=d1nftz[inopo]+d1nfpz[inopo] + d1nfcz[inopo];
*/
//-PY d1nfpx: pressure force   
// d1nfcx: drag force


//-PY pressure force
/*
		d1nfpx[inopo]=(ydn->d2nfv[0][inopo]-ydn->d2nvc[0][inopo])*ydn->d1nvct[inopo]+ydn->d2nfp[0][inopo];
		d1nfpy[inopo]=(ydn->d2nfv[1][inopo]-ydn->d2nvc[1][inopo])*ydn->d1nvct[inopo]+ydn->d2nfp[1][inopo];
		d1nfpz[inopo]=(ydn->d2nfv[2][inopo]-ydn->d2nvc[2][inopo])*ydn->d1nvct[inopo]+ydn->d2nfp[2][inopo];
*/
/*
		d1nfpx[inopo]=ydn->d2nfp[0][inopo];
		d1nfpy[inopo]=ydn->d2nfp[1][inopo];
		d1nfpz[inopo]=ydn->d2nfp[2][inopo];
  
*/		//PROBLEM is that d1nvcx > [18] is NaN

/*		 fprintf(fp,"Fx[%ld]=%f\tFy[%ld]=%f\tFz[%ld]=%f\n",
			 inopo,d1nftx[inopo],inopo,d1nfty[inopo],inopo,d1nftz[inopo]);
		 fprintf(fp,"nmct=%f\n",nmct);
		 fprintf(fp,"Vx[%ld]=%f\tVy[%ld]=%f\tVz[%ld]=%f\n",
			 inopo,dp1dv[inopo][0],inopo,dp1dv[inopo][1],inopo,dp1dv[inopo][2]);
		 fprintf(fp, "*********************************\n" );*/

		 if(index==2)
		   {
			 f1dt[0]=f1dt[0]+(d1nfvx[inopo]-d1nvcx[inopo]*absorp[inopo])*d1nvct[inopo];
			 f1dt[1]=f1dt[1]+(d1nfvy[inopo]-d1nvcy[inopo]*absorp[inopo])*d1nvct[inopo];
			 f1dt[2]=f1dt[2]+(d1nfvz[inopo]-d1nvcz[inopo]*absorp[inopo])*d1nvct[inopo];
			 f1dt1[0]=f1dt1[0]+(d1nfvx[inopo])*d1nvct[inopo];
			 f1dt1[1]=f1dt1[1]+(d1nfvy[inopo])*d1nvct[inopo];
			 f1dt1[2]=f1dt1[2]+(d1nfvz[inopo])*d1nvct[inopo];
			 f1dt2[0]=f1dt2[0]+(-d1nvcx[inopo])*absorp[inopo]*d1nvct[inopo];
			 f1dt2[1]=f1dt2[1]+(-d1nvcy[inopo])*absorp[inopo]*d1nvct[inopo];
			 f1dt2[2]=f1dt2[2]+(-d1nvcz[inopo])*absorp[inopo]*d1nvct[inopo];

			 f1dt[3]=f1dt[3]+denf[inopo]*d1nvct[inopo]*9.81;

			 f1bf[0]=f1bf[0]+d1nfpx[inopo];
			 f1bf[1]=f1bf[1]+d1nfpy[inopo];
			 f1bf[2]=f1bf[2]+d1nfpz[inopo];
			 /*
			 fprintf(fp,"inopo %ld %f %f %f %f %f\n",
				 inopo,d1nvcx[inopo],d1nvcy[inopo],d1nvcz[inopo],absorp[inopo],d1nvct[inopo]);
			 fprintf(fp,"Fx[%ld]=%f\tFy[%ld]=%f\tFz[%ld]=%f\n",
				 inopo,d1nftx[inopo],inopo,d1nfty[inopo],inopo,d1nftz[inopo]);
			 fprintf(fp,"nmct=%f\n",nmct);
			 fprintf(fp,"Vx[%ld]=%f\tVy[%ld]=%f\tVz[%ld]=%f\n",
				 inopo,dp1dv[inopo][0],inopo,dp1dv[inopo][1],inopo,dp1dv[inopo][2]);
		 */
		 }


//		 fprintf(fp, "1*********************************1 \n" ); //NaN happening after here
//		 fprintf(fp,"inopo %ld %f %f %f %f %f %f %f\n",
//			 inopo,d1nvcx[inopo],d1nvcy[inopo],d1nvcz[inopo],absorp[inopo],d1nvct[inopo],dp1dv[inopo][2], nmct);
//		 fprintf(fp,"Fx[%ld]=%f\tFy[%ld]=%f\tFz[%ld]=%f\n",
//			 inopo,d1nftx[inopo],inopo,d1nfty[inopo],inopo,d1nftz[inopo]);

	  if((i1bnvx[iprop])||(ABS(nmct)<EPSILON)) 
	  {
	       d1nftx[inopo]=R0;
	       dp1dv[inopo][0] =R0;
	  }
	  else 
	  {
//-PY changed it for 3D_fracture_coupling_with_multiphase	   
              d1nftx[inopo] += (d1bnfx[iprop] );
          //    d1nftx[inopo] += (d1bnfx[iprop] + d1nfcx[inopo]+ d1nfpx[inopo]);// + d1kax[iprop]*nmct);
	       dp1dv[inopo][0] = d1nftx[inopo]/nmct;


//fprintf(fp,"d1bnfx[iprop] d1nfcx[inopo] d1kax[iprop]    %f %f %f\n", d1bnfx[iprop], d1nfcx[inopo], d1kax[iprop]);

	  }

	  if((i1bnvy[iprop])||(ABS(nmct)<EPSILON)) 
	  {
	       d1nfty[inopo]=R0;
	       dp1dv[inopo][1] =R0;
	  }
	  else 
	  {
//-PY changed it for 3D_fracture_coupling_with_multiphase
	       d1nfty[inopo] += (d1bnfy[iprop] );
   //            d1nfty[inopo] += (d1bnfy[iprop] + d1nfcy[inopo]+ d1nfpy[inopo]);// + d1kay[iprop]*nmct);
	       dp1dv[inopo][1] = d1nfty[inopo]/nmct;

//fprintf(fp,"d1bnfy[iprop] d1nfcy[inopo]  d1kay[iprop] %f %f %f\n", d1bnfy[iprop], d1nfcy[inopo], d1kay[iprop]);
	  }

	  if((i1bnvz[iprop])||(ABS(nmct)<EPSILON))
	  {
	       d1nftz[inopo]=R0;
	       dp1dv[inopo][2] =R0;
	  }
	  else
	  {
//-PY changed it for 3D_fracture_coupling_with_multiphase
	       d1nftz[inopo] += (d1bnfz[iprop] );   //-PY divided by 3 to make gravity correct.
//              d1nftz[inopo] += (d1bnfz[iprop] + d1nfcz[inopo]+ d1nfpz[inopo]);// + d1kaz[iprop]*nmct);
	       dp1dv[inopo][2] = d1nftz[inopo]/nmct;

//fprintf(fp,"d1bnfz[iprop] d1nfcz[inopo]  d1kaz[iprop] %f %f %f\n", d1bnfz[iprop], d1nfcz[inopo], d1kaz[iprop]);
	  }

//		 fprintf(fp, "2*********************************2 \n" );
//		 fprintf(fp,"inopo %ld %f %f %f %f %f %f %f\n",
//			 inopo,d1nvcx[inopo],d1nvcy[inopo],d1nvcz[inopo],absorp[inopo],d1nvct[inopo], dp1dv[inopo][2], nmct);
//		 fprintf(fp,"Fx[%ld]=%f\tFy[%ld]=%f\tFz[%ld]=%f\n",
//			 inopo,d1nftx[inopo],inopo,d1nfty[inopo],inopo,d1nftz[inopo]);
     }


//    fclose(fp);
}


/**********************************************************************/
/**********************************************************************/


/*Z one-pass solver
*/
static INLINE void solve(
     DBL *d1nftx0, DBL *d1nfty0, DBL *d1nftz0,
     DBL *d1nftx1, DBL *d1nfty1, DBL *d1nftz1,
     DBL *d1nftx2, DBL *d1nfty2, DBL *d1nftz2,
     INT *i1bnvx, INT *i1bnvy, INT *i1bnvz,
     DPT *dp1dv,  
     DPT *dp1am,  
     INT nnopo,
     DBL *d1nmct,    
     INT nelem, 
     DBL *d1emct,   
     INT **i2elto, 
     DBL **d2csmm, 
     DBL dcurelx,
     INT *i1nopr,INT iqua)
{
     INT ielem, inopo, i, j, ip, jp, *i1elto,iprop;
     DBL emct, mass, *d1csmm;


     /*Z zero acceleration work array (old method was less cohesive)
      */
     TzDPT1(dp1am, nnopo);

     /*Z calculate force
      */
     for(ielem=0; ielem<nelem; ielem++)
     {
	  i1elto = i2elto[ielem];
	  emct   = d1emct[ielem];
	
	  if(iqua==1)
	  {
	       for(i=0; i<NELNO; i++)
	       {  
		    ip     = i1elto[i];
		    d1csmm = d2csmm[i];

		    for(j=0; j<NELNO; j++)
		    {
			 jp   = i1elto[j];
			 mass = d1csmm[j]*emct;

			 dp1am[ip][0] += mass*dp1dv[jp][0];
			 dp1am[ip][1] += mass*dp1dv[jp][1];
			 dp1am[ip][2] += mass*dp1dv[jp][2];
		    }
	       }
	  }
	  else
	  {
	       for(i=0; i<NNODEX; i++)
	       {  
		    ip     = i1elto[i];
      
		    for(j=0; j<NNODEX; j++)
		    {
			 jp   = i1elto[j];
			 if(i==j){mass = 0.01*emct;}
			 else
			 {mass = 0.005*emct;}



			 dp1am[ip][0] += mass*dp1dv[jp][0];
			 dp1am[ip][1] += mass*dp1dv[jp][1];
			 dp1am[ip][2] += mass*dp1dv[jp][2];
		    }
	       }
	  }
     }

     /*Z calculate velocity*/
     for(inopo=0; inopo<nnopo; inopo++)
     {
	  iprop=i1nopr[inopo];
	  if(i1bnvx[iprop]) /* supplied velocity     */
	       dp1am[inopo][0]=R0;
 
	  if(i1bnvy[iprop]) /* supplied velocity     */
	       dp1am[inopo][1]=R0;
      
	  if(i1bnvz[iprop]) /* supplied velocity     */
	       dp1am[inopo][2]=R0;
 
	  d1nftx2[inopo] = d1nftx0[inopo] + (d1nftx1[inopo] - dp1am[inopo][0])*dcurelx; 
	  d1nfty2[inopo] = d1nfty0[inopo] + (d1nfty1[inopo] - dp1am[inopo][1])*dcurelx; 
	  d1nftz2[inopo] = d1nftz0[inopo] + (d1nftz1[inopo] - dp1am[inopo][2])*dcurelx;





		if(ABS(d1nmct[inopo])<EPSILON)
		{     
		dp1dv[inopo][0]=R0;
		dp1dv[inopo][1]=R0;
		dp1dv[inopo][2]=R0;
		}
		else
		{  	       
	  dp1dv[inopo][0] = d1nftx2[inopo]/d1nmct[inopo];
	  dp1dv[inopo][1] = d1nfty2[inopo]/d1nmct[inopo];
	  dp1dv[inopo][2] = d1nftz2[inopo]/d1nmct[inopo];
	  }
     }


}


/**********************************************************************/
/**********************************************************************/


/*Z calculate updated nodal velocity and coordinates
*/
static INLINE void nupdate(
     DBL *d1nccx, DBL *d1nccy, DBL *d1nccz, 
     DBL *d1nvcx, DBL *d1nvcy, DBL *d1nvcz, 
     DBL *d1bnvx, DBL *d1bnvy, DBL *d1bnvz, 
     INT *i1bnvx, INT *i1bnvy, INT *i1bnvz, 
     DBL *d1bnfx, DBL *d1bnfy, DBL *d1bnfz, 
     DBL *d1kax,  DBL *d1kay,  DBL *d1kaz, 
     DPT *dp1dv, 
     INT nnopo,
     INT *i1nopr, 
     DBL dcstec,DBL dctime)
{
     INT inopo, iprop;



     for(inopo=0; inopo<nnopo; inopo++)
     {    
	  iprop = i1nopr[inopo];

	  /*Z current velocity */
	  d1nvcx[inopo] += dp1dv[inopo][0]*dcstec;
	  d1nvcy[inopo] += dp1dv[inopo][1]*dcstec;
	  d1nvcz[inopo] += dp1dv[inopo][2]*dcstec;


      d1nvcx[inopo] += d1kax[iprop]*dcstec;
	  d1nvcy[inopo] += d1kay[iprop]*dcstec;
	  d1nvcz[inopo] += d1kaz[iprop]*dcstec;




	  if(i1bnvx[iprop]==1)
	       d1nvcx[inopo] = d1bnvx[iprop];
	  else if(i1bnvx[iprop]==2)
	       d1nvcx[inopo] = d1bnvx[iprop]*COS(d1bnfx[iprop]*dctime*R2*MYPI);

	  if(i1bnvy[iprop]==1)
	       d1nvcy[inopo] = d1bnvy[iprop];
	  else if(i1bnvy[iprop]==2)
	       d1nvcy[inopo] = d1bnvy[iprop]*COS(d1bnfy[iprop]*dctime*R2*MYPI);
		    
	  if(i1bnvz[iprop]==1)
	       d1nvcz[inopo] = d1bnvz[iprop];
	  else if(i1bnvz[iprop]==2)
	       d1nvcz[inopo] = d1bnvz[iprop]*COS(d1bnfz[iprop]*dctime*R2*MYPI);



	  /*Z new coordinates */
	  d1nccx[inopo] += d1nvcx[inopo]*dcstec; 
	  d1nccy[inopo] += d1nvcy[inopo]*dcstec;
	  d1nccz[inopo] += d1nvcz[inopo]*dcstec;
     }


}


/**********************************************************************/
/**********************************************************************/


/* mechanical solver for 3D nodes with x,y,z d.o.f.
*/
static void Ysd3MEC( YDN ydn,
     DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
     DBL *d1nfcx, DBL *d1nfcy, DBL *d1nfcz,
     DBL *d1nftx, DBL *d1nfty, DBL *d1nftz,
     DBL *d1nvcx, DBL *d1nvcy, DBL *d1nvcz,
     INT *i1bnvx, INT *i1bnvy, INT *i1bnvz,
     DBL *d1bnvx, DBL *d1bnvy, DBL *d1bnvz,
     DBL *d1bnfx, DBL *d1bnfy, DBL *d1bnfz,
     DBL *d1kax,  DBL *d1kay,  DBL *d1kaz,   //-PY gravity
     DBL **d2csmm, 
     INT nnopo, 
     DBL *d1nmct,
     INT *i1nopr,
     INT *i1nobf,
     DBL *d1emct,
     INT nelem,
     INT **i2elto,
     INT initer,
     DBL dcstec,DBL dctime, 
     DBL dcurelx,INT iqua,
	 DBL *d1nfpx,  DBL *d1nfpy,  DBL *d1nfpz,
	 DBL *d1nfvx,  DBL *d1nfvy,  DBL *d1nfvz,
	 DBL *d1nvct,  DBL *absorp,
	 DBL *denf,     INT index,
	DBL *f1dt, DBL *f1dt1, DBL *f1dt2
     ) 
{ 


     DBL *d1nftx1, *d1nfty1, *d1nftz1;
     DPT *dp1dv, *dp1am;

     /*Z allocate temporary working arrays (could be static in YDK) */
     d1nftx1 = TalDBL1(nnopo);
     d1nfty1 = TalDBL1(nnopo);
     d1nftz1 = TalDBL1(nnopo);
     dp1dv   = TalDPT1(nnopo);
     dp1am   = TalDPT1(nnopo);
 
     /*Z calculate current nodal force nft, and change in velocity dv */
     nfvel(ydn, d1bnfx, d1bnfy, d1bnfz, 
	   d1nfcx, d1nfcy, d1nfcz, 
	   d1nftx, d1nfty, d1nftz, 
	   //d1kax,  d1kay,  d1kaz, 
	   i1bnvx, i1bnvy, i1bnvz,
	   dp1dv, 
	   nnopo,
	   d1nmct, 
	   i1nopr,
	   d1nfvx,d1nfvy,d1nfvz,
		   d1nfpx,d1nfpy,d1nfpz,
		   d1nvcx, d1nvcy, d1nvcz,
		   d1nvct,absorp,denf,index,
		f1dt,f1dt1,f1dt2);

     /*Z solver pass 1 (or could have simpler solve1() and solve2()) */
     if(initer>0)
     {
	  solve(d1nftx,  d1nfty,  d1nftz, 
		d1nftx,  d1nfty,  d1nftz,  
		d1nftx1, d1nfty1, d1nftz1,
		i1bnvx, i1bnvy, i1bnvz,
		dp1dv, 
		dp1am, 
		nnopo,
		d1nmct,
		nelem,
		d1emct, 
		i2elto, 
		d2csmm, 
		dcurelx,
		i1nopr,iqua);

	  /*Z solver pass 2 (this could be compile-time constant: #if(INITER==2) */
      if(initer == 2)
      {
	   solve(d1nftx,  d1nfty,  d1nftz,
		 d1nftx1, d1nfty1, d1nftz1,
		 d1nftx1, d1nfty1, d1nftz1,
		 i1bnvx, i1bnvy, i1bnvz,
		 dp1dv,
		 dp1am,
		 nnopo,
		 d1nmct,
		 nelem,
		 d1emct,
		 i2elto,
		 d2csmm,
		 dcurelx,
		 i1nopr,iqua);
	  }
     }
     /*Z update nodal velocity nvc and coordinates ncc */
     nupdate(d1nccx, d1nccy, d1nccz,
	     d1nvcx, d1nvcy, d1nvcz, 
	     d1bnvx, d1bnvy, d1bnvz, 
	     i1bnvx, i1bnvy, i1bnvz, 
	     d1bnfx, d1bnfy, d1bnfz, 
	     d1kax,  d1kay,  d1kaz, 
	     dp1dv, 
	     nnopo,
	     i1nopr,
	     dcstec,dctime);

     /*Z release temporary working arrays (could be static in YDK) */
     FREE(dp1am);
     FREE(dp1dv);
     FREE(d1nftz1);
     FREE(d1nfty1);
     FREE(d1nftx1);
}


static void Ysd3Rig( ypar,ydb,ydc,ydn,ysp)
YPAR ypar; YDB ydb; YDC ydc; YDN ydn; YSP ysp;
{

//     FILE *fr;
//     fr=fopen("rigid_test.txt","w+");

/* define variable     */
  INT i,j,k,id;
  DBL smass;//,smoment,spera;
  DBL osp[3];
  DBL dis[3];
  DBL *d1u,*d1v,*d1w,*d1cosine,*d1sine;
  INT *i1tran,icon[3],indb;
  DBL ninertia,nvr,zelta,d1ang[3],tmp;
  DBL t1,t2,t3,d1sgray[3];
  DBL f[3],a[3],damp,d1nvcx,d1nvcy,d1nvcz;
  DBL tmpx,tmpy,tmpz,fx,fy,fz;
  INT ia;
  DBL dxc,dyc,dzc;
  DBL mxc,myc,mzc;
  DBL mxc_t,myc_t,mzc_t;
//  FILE *FOUT=FILENULL;
  
//  FOUT=fopen("sdw.txt","w");
  d1u=TalDBL1(ysp->nspd);
  d1v=TalDBL1(ysp->nspd);
  d1w=TalDBL1(ysp->nspd);
  d1cosine=TalDBL1(ysp->nspd);
  d1sine=TalDBL1(ysp->nspd);
  i1tran=TalINT1(ysp->nspd);

  //ysd->dsengy=R0;



 //-PY initialize the ysp->d2pfd
int ia1;
    for(i=0; i<ydn->nnopo; i++)
    ia1=ydn->i1nind[i]-1;
	{
         ysp->d2pfd[0][ia1]=0.0;
         ysp->d2pfd[1][ia1]=0.0;
         ysp->d2pfd[2][ia1]=0.0;

         }
  

  
  d1sgray[0]=ydc->dcgrax;
  d1sgray[1]=ydc->dcgray;
  d1sgray[2]=ydc->dcgraz;
	

	for(i=0;i<ydn->nnopi;i++)
	{
		
	ia=ydn->i1nind[i]-1;
	//distance between node and centroid	
	dxc=ydn->d2ncc[0][i]-ysp->d2pp[0][ia];
	dyc=ydn->d2ncc[1][i]-ysp->d2pp[1][ia];
	dzc=ydn->d2ncc[2][i]-ysp->d2pp[2][ia];

//t=r*f
	
	tmpx=ysp->d2pvr[0][ia]*(ysp->d2prinn1[1][ia]*dzc-ysp->d2prinn1[2][ia]*dyc)
		+ysp->d2pvr[1][ia]*(ysp->d2prinn2[1][ia]*dzc-ysp->d2prinn2[2][ia]*dyc)
		+ysp->d2pvr[2][ia]*(ysp->d2prinn3[1][ia]*dzc-ysp->d2prinn3[2][ia]*dyc);
	tmpy=ysp->d2pvr[0][ia]*(ysp->d2prinn1[2][ia]*dxc-ysp->d2prinn1[0][ia]*dzc)
		+ysp->d2pvr[1][ia]*(ysp->d2prinn2[2][ia]*dxc-ysp->d2prinn2[0][ia]*dzc)
		+ysp->d2pvr[2][ia]*(ysp->d2prinn3[2][ia]*dxc-ysp->d2prinn3[0][ia]*dzc);
	tmpz=ysp->d2pvr[0][ia]*(ysp->d2prinn1[0][ia]*dyc-ysp->d2prinn1[1][ia]*dxc)
		+ysp->d2pvr[1][ia]*(ysp->d2prinn2[0][ia]*dyc-ysp->d2prinn2[1][ia]*dxc)
		+ysp->d2pvr[2][ia]*(ysp->d2prinn3[0][ia]*dyc-ysp->d2prinn3[1][ia]*dxc);
	
	d1nvcx=tmpx+ysp->d2pvt[0][ia];
	d1nvcy=tmpy+ysp->d2pvt[1][ia];
	d1nvcz=tmpz+ysp->d2pvt[2][ia];

//fprintf(fr, " tmpx,  tmpy,  tmpz   %lf  %lf  %lf \n", tmpx,  tmpy,  tmpz  );

//fprintf(fr, " 1ysp->d2pvt[0][ia],  ysp->d2pvt[1][ia],  ysp->d2pvt[2][ia]   %lf  %lf  %lf \n", ysp->d2pvt[0][ia],  ysp->d2pvt[1][ia],  ysp->d2pvt[2][ia] );

        ydn->d2nvc[0][i]= d1nvcx;
        ydn->d2nvc[1][i]= d1nvcy;
        ydn->d2nvc[2][i]= d1nvcz;

//-PY changed it for 3D_fracture_coupling_with_multiphase
		
//		fx=(ydn->d2nfv[0][i]-ydn->d2nvc[0][i])*ydn->d1nvct[i]+ydn->d2nfp[0][i];
//		fy=(ydn->d2nfv[1][i]-ydn->d2nvc[1][i])*ydn->d1nvct[i]+ydn->d2nfp[1][i];
//		fz=(ydn->d2nfv[2][i]-ydn->d2nvc[2][i])*ydn->d1nvct[i]+ydn->d2nfp[2][i];
		


//fprintf(fr, " 1ydn->d2nvc[0][i],  ydn->d2nvc[1][i],  ydn->d2nvc[2][i]   %lf  %lf  %lf \n", ydn->d2nvc[0][i],  ydn->d2nvc[1][i],  ydn->d2nvc[2][i] );


//-PY changed it.

       		
//		fx=(ydn->d2nfv[0][i]-ydn->d2nvc[0][i])*ydn->d1nvct[i]+ydn->d2nfp[0][i]+ydn->d2nfc[0][i];
//		fy=(ydn->d2nfv[1][i]-ydn->d2nvc[1][i])*ydn->d1nvct[i]+ydn->d2nfp[1][i]+ydn->d2nfc[1][i];
//		fz=(ydn->d2nfv[2][i]-ydn->d2nvc[2][i])*ydn->d1nvct[i]+ydn->d2nfp[2][i]+ydn->d2nfc[2][i];




		fx=ydn->d2nfp[0][i]+ydn->d2nfc[0][i]; //  +ydn->d1nmct[i]*d1sgray[0];  //-PY d2nfp is the pressure force 
		fy=ydn->d2nfp[1][i]+ydn->d2nfc[1][i]; //  +ydn->d1nmct[i]*d1sgray[1];
		fz=ydn->d2nfp[2][i]+ydn->d2nfc[2][i]; //  +ydn->d1nmct[i]*d1sgray[2];	


//fprintf(fr, " d1sgray[0],  d1sgray[1],  d1sgray[2]   %lf  %lf  %lf \n", d1sgray[0],  d1sgray[1],  d1sgray[2] );
//fprintf(fr, " ydn->d2nfp[0][i],  ydn->d2nfp[1][i],  ydn->d2nfp[2][i]   %lf  %lf  %lf \n", ydn->d2nfp[0][i],  ydn->d2nfp[1][i],  ydn->d2nfp[2][i] );
	mxc_t=fz*dyc-fy*dzc;
	myc_t=fx*dzc-fz*dxc;
	mzc_t=fy*dxc-fx*dyc;
						
	mxc=mxc_t*ysp->d2prinn1[0][ia]+myc_t*ysp->d2prinn1[1][ia]+mzc_t*ysp->d2prinn1[2][ia];
	myc=mxc_t*ysp->d2prinn2[0][ia]+myc_t*ysp->d2prinn2[1][ia]+mzc_t*ysp->d2prinn2[2][ia];
	mzc=mxc_t*ysp->d2prinn3[0][ia]+myc_t*ysp->d2prinn3[1][ia]+mzc_t*ysp->d2prinn3[2][ia];						   
	
/*temporarily block torque*/

	ysp->d2pmd[0][ia]=ysp->d2pmd[0][ia]+mxc;
	ysp->d2pmd[1][ia]=ysp->d2pmd[1][ia]+myc;
	ysp->d2pmd[2][ia]=ysp->d2pmd[2][ia]+mzc;

						   

        ysp->d2pfd[0][ia]=ysp->d2pfd[0][ia]+fx;
        ysp->d2pfd[1][ia]=ysp->d2pfd[1][ia]+fy;
        ysp->d2pfd[2][ia]=ysp->d2pfd[2][ia]+fz;



/*
        ysp->d2pfd[0][ia]+=fx;
        ysp->d2pfd[1][ia]+=fy;
        ysp->d2pfd[2][ia]+=fz;
*/					   
    }

	
 for(j=0;j<ysp->nspd;j++)
 {
  /*calculate the volume of sphere, mass and moment of sphere */
 /*V=4/3*pi*r*r*r                                            */
 /*M=V*RO                                                    */
 /*Moment=2/5*M*r*r                                          */
	 i1tran[j]=0;
     indb=ysp->i1sppr[j];
     damp=0.5;  //hardwired value, to be changed later
	 f[0]=ydb->d1bnfx[indb];
	 f[1]=ydb->d1bnfy[indb];
	 f[2]=ydb->d1bnfz[indb];
	 a[0]=ydb->d1bnax[indb];
	 a[1]=ydb->d1bnay[indb];
	 a[2]=ydb->d1bnaz[indb];

   smass=ysp->d1mass[j];
   //   spera=
 for(i=0;i<ydn->nnodim;i++)
   {
   ysp->d2poldvr[i][j]=ysp->d2pvr[i][j];
   ysp->d2poldp[i][j]=ysp->d2pp[i][j];
   icon[i]=0;
 }

 if(ydb->i1bnvx[indb]>0)icon[0]=ydb->i1bnvx[indb];
 if(ydb->i1bnvy[indb]>0)icon[1]=ydb->i1bnvy[indb];
 if(ydb->i1bnvz[indb]>0)icon[2]=ydb->i1bnvz[indb];


 for(i=0;i<ydn->nnodim;i++)
   {
   osp[i]=ysp->d2pp[i][j];
          
/* calculate updated translational velocity of sphere    */
/* u=u0+dt*force/mass;                                    */

   	 if(ysp->i1con[j]==0)
	 {
		ysp->d2pvt[i][j]=ypar->d1vc[i];
	    
	 }
	 else
	 {

   if(icon[i]>0)
   {
	   if(icon[i]==1)
	   {
	   if(i==0)ysp->d2pvt[i][j]=ydb->d1bnvx[indb];
	   if(i==1)ysp->d2pvt[i][j]=ydb->d1bnvy[indb];
	   if(i==2)ysp->d2pvt[i][j]=ydb->d1bnvz[indb];
	   }
	   else if(icon[i]==2)
	   {

	    if(i==0) ysp->d2pvt[i][j] = ydb->d1bnvx[indb]*SIN(ydb->d1bnfx[indb]*ydc->dctime*R2*MYPI);
	    if(i==1) ysp->d2pvt[i][j] = ydb->d1bnvy[indb]*SIN(ydb->d1bnfy[indb]*ydc->dctime*R2*MYPI);
	    if(i==2) ysp->d2pvt[i][j] = ydb->d1bnvz[indb]*SIN(ydb->d1bnfz[indb]*ydc->dctime*R2*MYPI);
	   }

   }
   else
   {
	   ysp->d2ptf[i][j]=ysp->d2ptfc[i][j]+f[i]+ysp->d2pfd[i][j]+smass*(a[i])-damp*ABS(ysp->d2ptfc[i][j])*sign(ysp->d2pvt[i][j]);
   
   ysp->d2pvt[i][j]=ysp->d2pvt[i][j]+ydc->dcstec*(ysp->d2ptf[i][j]/smass+d1sgray[i]);

//fprintf(fr, " ysp->d2ptf[i][j],  %lf \n", ysp->d2ptf[i][j] );
   }






      if(ABS(ysp->d2pvt[i][j])>=10.0)
   {
   CHRw(stderr,"MOVE TOO FAST IN Y DIRECTION");
   }

   
//   CHRwcr(stderr);
//DBLw(stderr,ysp->d2ptf[i][j],10);
//CHRwcr(stderr);
//DBLw(stderr,ysp->d2pvt[i][j],10);

//   DBLw(stdout, ysp->d2pvt[i][j],10);
/* calculate updated rotational velocity of sphere    */
/* Angel velocity=old velocity+torque/moment          */
   k=i+1;
   if(k>=ydn->nnodim)k=0;
   if(k==(ydn->nnodim-1))
   {
      ninertia=ysp->d2prine[0][j];
   nvr=ysp->d2poldvr[0][j];
   }
   else
   {
   ninertia=ysp->d2prine[k+1][j];
   nvr=ysp->d2poldvr[k+1][j];
   } 
   if(icon[i]>0)
   {ysp->d2pvr[i][j]=R0;
   }
   else
   {
   ysp->d2pvr[i][j]=ysp->d2pvr[i][j]+ydc->dcstec*(ysp->d2pm[i][j]+ysp->d2pmd[i][j]-damp*ABS(ysp->d2pm[i][j])*sign(ysp->d2pvr[i][j])+
   (ysp->d2prine[k][j]-ninertia)*ysp->d2poldvr[k][j]*nvr)/ysp->d2prine[i][j];
   }
	 }


//fprintf(fr, " 2ysp->d2pvt[0][ia],  ysp->d2pvt[1][ia],  ysp->d2pvt[2][ia]   %lf  %lf  %lf \n", ysp->d2pvt[0][ia],  ysp->d2pvt[1][ia],  ysp->d2pvt[2][ia] );

   
/* calculate updated position of sphere     */
/* x=x0+dt*u                                */
   ysp->d2pp[i][j]=ysp->d2pp[i][j]+ydc->dcstec*ysp->d2pvt[i][j];

//CHRwcr(stderr);
//DBLw(stderr,ysp->d2pp[i][j],10);
//CHRwcr(stderr);
//DBLw(stderr,ysp->d2pvt[i][j],10);

 
/* calculate updated angel of sphere     */
/* Angel=old angel+angle velocity by time step   */
   ysp->d2pa[i][j]=ysp->d2pa[i][j]+ydc->dcstec*ysp->d2pvr[i][j];

/* calculate the moving distance of sphere during last time step   */
   d1ang[i]=ydc->dcstec*ysp->d2pvr[i][j];
   dis[i]=ABS(ysp->d2pp[i][j]-osp[i]);

   //   ysp->dsengy=ysp->dsengy+0.5*smass*ysp->d2pvt[i][j]*ysp->d2pvt[i][j]+
   //   0.5*smoment*ysd->d2svr[i][j]*ysd->d2svr[i][j];

 } 
 //   if(dis[0]>=spera)CHRw(stderr,"MOVE TOO FAST IN X DIRECTION");
 //   if(dis[1]>=spera)
 //   {
 //   CHRw(stderr,"MOVE TOO FAST IN Y DIRECTION");
 //   }
 //   if(dis[2]>=spera)CHRw(stderr,"MOVE TOO FAST IN Z DIRECTION");

zelta=SQRT(d1ang[0]*d1ang[0]+d1ang[1]*d1ang[1]+d1ang[2]*d1ang[2]);
if(zelta>EPSILON)
{
i1tran[j]=1;
t1=ysp->d2prinn1[0][j];
t2=ysp->d2prinn2[0][j];
t3=ysp->d2prinn3[0][j];

d1u[j]=(d1ang[0]*ysp->d2prinn1[0][j]+d1ang[1]*ysp->d2prinn2[0][j]+
        d1ang[2]*ysp->d2prinn3[0][j])/zelta;
        
t1=ysp->d2prinn1[1][j];
t2=ysp->d2prinn2[1][j];
t3=ysp->d2prinn3[1][j];
        
d1v[j]=(d1ang[0]*ysp->d2prinn1[1][j]+d1ang[1]*ysp->d2prinn2[1][j]+
        d1ang[2]*ysp->d2prinn3[1][j])/zelta;
        
t1=ysp->d2prinn1[2][j];
t2=ysp->d2prinn2[2][j];
t3=ysp->d2prinn3[2][j];
        
d1w[j]=(d1ang[0]*ysp->d2prinn1[2][j]+d1ang[1]*ysp->d2prinn2[2][j]+
        d1ang[2]*ysp->d2prinn3[2][j])/zelta;
d1cosine[j]=COS(zelta);
d1sine[j]=SIN(zelta);

V3DNor(tmp,d1u[j],d1v[j],d1w[j]);
V3DRot(&ysp->d2prinn1[0][j],&ysp->d2prinn1[1][j],&ysp->d2prinn1[2][j],
       ysp->d2prinn1[0][j],ysp->d2prinn1[1][j],ysp->d2prinn1[2][j],
       d1u[j],d1v[j],d1w[j],d1cosine[j],d1sine[j]);
V3DRot(&ysp->d2prinn2[0][j],&ysp->d2prinn2[1][j],&ysp->d2prinn2[2][j],
       ysp->d2prinn2[0][j],ysp->d2prinn2[1][j],ysp->d2prinn2[2][j],
       d1u[j],d1v[j],d1w[j],d1cosine[j],d1sine[j]);
V3DRot(&ysp->d2prinn3[0][j],&ysp->d2prinn3[1][j],&ysp->d2prinn3[2][j],
       ysp->d2prinn3[0][j],ysp->d2prinn3[1][j],ysp->d2prinn3[2][j],
       d1u[j],d1v[j],d1w[j],d1cosine[j],d1sine[j]);
       
V3DNor(tmp,ysp->d2prinn1[0][j],ysp->d2prinn1[1][j],ysp->d2prinn1[2][j]);
V3DNor(tmp,ysp->d2prinn2[0][j],ysp->d2prinn2[1][j],ysp->d2prinn2[2][j]);
V3DNor(tmp,ysp->d2prinn3[0][j],ysp->d2prinn3[1][j],ysp->d2prinn3[2][j]);
}
}



 for(j=0;j<ydn->nnopo;j++)
 {
  /*calculate the volume of sphere, mass and moment of sphere */
 /*V=4/3*pi*r*r*r                                            */
 /*M=V*RO                                                    */
 /*Moment=2/5*M*r*r                                          */

	 /*
 spera=ysdp->d1spera[ysd->i1nopr[j]];
 smass=(R4/R3)*MYPI*spera*spera*spera*ysdp->d1spero[ysd->i1nopr[j]];
 smoment=(R2/R5)*smass*spera*spera;*/
 id=ydn->i1nind[j]-1;
 for(i=0;i<ydn->nnodim;i++)
 {
 
/* store the old position of sphere   */
   osp[i]=ydn->d2ncc[i][j];
   ydn->d2ncc[i][j]=ydn->d2ncc[i][j]-ysp->d2poldp[i][id];

          
/* calculate updated translational velocity of sphere    */
/* u=u0+dt*force/mass;                                    */

//   ysd->d2stf[i][j]=ysd->d2stfc[i][j]+smass*ysdc->d1sgray[i];
   
   ydn->d2nvc[i][j]=ysp->d2pvt[i][id];

/* calculate updated rotational velocity of sphere    */
/* Angel velocity=old velocity+torque/moment          */

//   ysd->d2svr[i][j]=ysd->d2svr[i][j]+ysdc->dscstec*ysd->d2sm[i][j]/smoment;
   
   
/* calculate updated position of sphere     */
/* x=x0+dt*u                                */
//   ysd->d2sp[i][j]=ysd->d2sp[i][j]+ysdc->dscstec*ysd->d2svt[i][j];
 
/* calculate updated angel of sphere     */
/* Angel=old angel+angle velocity by time step   */
//   ysd->d2sa[i][j]=ysd->d2sa[i][j]+ysdc->dscstec*ysd->d2svr[i][j];

/* calculate the moving distance of sphere during last time step   */

   //   ysd->dsengy=ysd->dsengy+0.5*smass*ysd->d2svt[i][j]*ysd->d2svt[i][j]+
   //   0.5*smoment*ysd->d2svr[i][j]*ysd->d2svr[i][j];

//CHRwcr(stderr);
//DBLw(stderr,ysd->d2sp[i][j],10);

 }
if(i1tran[id]==1)
{
V3DRot(&ydn->d2ncc[0][j],&ydn->d2ncc[1][j],&ydn->d2ncc[2][j],
ydn->d2ncc[0][j],ydn->d2ncc[1][j],ydn->d2ncc[2][j],
d1u[id],d1v[id],d1w[id],d1cosine[id],d1sine[id]);
}
ydn->d2ncc[0][j]=ydn->d2ncc[0][j]+ysp->d2pp[0][id];
ydn->d2ncc[1][j]=ydn->d2ncc[1][j]+ysp->d2pp[1][id];
ydn->d2ncc[2][j]=ydn->d2ncc[2][j]+ysp->d2pp[2][id];




/*
   dis[0]=ABS(ydn->d2ncc[0][j]-osp[0]);
   dis[1]=ABS(ydn->d2ncc[1][j]-osp[1]);
   dis[2]=ABS(ydn->d2ncc[2][j]-osp[2]);


   if(dis[0]>=spera)CHRw(stderr,"MOVE TOO FAST IN X DIRECTION");
   if(dis[1]>=spera)
   {
   CHRw(stderr,"MOVE TOO FAST IN Y DIRECTION");
   }
   if(dis[2]>=spera)CHRw(stderr,"MOVE TOO FAST IN Z DIRECTION");
*/

 }

//	fclose(FOUT);
//        fclose(fr);

FREE(i1tran);
FREE(d1sine);
FREE(d1cosine);
FREE(d1w);
FREE(d1v);
FREE(d1u);






 }

/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/


/*! \brief explicit solver of equations
 *  \param[in]     ydc control database
 *  \param[in]     yde element database
 *  \param[in,out] ydn nodal database
 *  \param[out]    ydo output database
 *  \param[in]     ydb boundary condition database
 *  \param[in]     ydk constants (for d2csmm consistent mass matrix)
 *  \par Details:
 *  Ysd() calls Ysd3MEC() which is the mechanical solver for 3D nodes 
 *  with x,y,z d.o.f..  Adjustments are made afterward for the history 
 *  variable of total kinetic energy.
 *
 */
void Ysd(YDC ydc, YDE yde, YDN ydn, YDO ydo, YDB ydb, YDK ydk, YDP ydp,
		DBL *absorp,DBL *denf, INT index, DBL *f1dt, DBL *f1dt1, DBL *f1dt2,
        YPAR ypar, YSP ysp)
{
     INT i, ihys, iqua;
     DBL Ek, ohys, stprev, ohyp;


//-PY delete it for the rigid version
/*
     if(ydp->i1ptyp[0]==YTE3TET4ELS)
     {
	  iqua=0;
     }
     else if(ydp->i1ptyp[0]==YTE3TET10ELS)
     {
	  iqua=1;
     }
*/

     /* mechanical solver for 3D nodes with x,y,z d.o.f.  */


 //    printf("IRIGID %i \n\n\n", ydc->irigid);


 if(ydc->irigid==0)
   {
   Ysd3MEC( ydn,
	  ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2ncc[2],
	  ydn->d2nfc[0], ydn->d2nfc[1], ydn->d2nfc[2],
	  ydn->d2nft[0], ydn->d2nft[1], ydn->d2nft[2], 
	  ydn->d2nvc[0], ydn->d2nvc[1], ydn->d2nvc[2],
	  ydb->i1bnvx,   ydb->i1bnvy,   ydb->i1bnvz,
	  ydb->d1bnvx,   ydb->d1bnvy,   ydb->d1bnvz,
	  ydb->d1bnfx,   ydb->d1bnfy,   ydb->d1bnfz,
	  ydk->d2ka[0],  ydk->d2ka[1],  ydk->d2ka[2], 
	  ydk->d2csmm,
	  ydn->nnopo,
	  ydn->d1nmct,
	  ydn->i1nopr,
	  ydn->i1nobf,
	  yde->d1emct,
	  yde->nelem,
	  yde->i2elto,
	  ydc->initer,
	  ydc->dcstec,ydc->dctime,
	  ydc->dcurelx,iqua,
			   ydn->d2nfp[0],ydn->d2nfp[1],ydn->d2nfp[2], 
			 ydn->d2nfv[0],ydn->d2nfv[1],ydn->d2nfv[2],
			 ydn->d1nvct,absorp,denf,index,
			 f1dt,f1dt1,f1dt2);
   }

 else if(ydc->irigid==1)
     {Ysd3Rig(ypar,ydb,ydc,ydn,ysp);
	 nodevel(     	yde->nelemi, ydn->nnopi,ysp->nspd,
			ydn->d2ncc,  ydn->d2nvc,  ydn->d2nfc,
			yde->i1elpr, yde->i2elto, ydp->i1ptyp,
			ysp->d2pvr[0],ysp->d2pvr[1],ysp->d2pvr[2],
			ysp->d2prinn1[0],ysp->d2prinn1[1],ysp->d2prinn1[2],
			ysp->d2prinn2[0],ysp->d2prinn2[1],ysp->d2prinn2[2],
			ysp->d2prinn3[0],ysp->d2prinn3[1],ysp->d2prinn3[2],
			ysp->d2pp[0],ysp->d2pp[1],ysp->d2pp[2],
			ysp->d2pvt[0],ysp->d2pvt[1],ysp->d2pvt[2],ydn->i1nind,
			ysp->d2pa[0],ysp->d2pa[1],ysp->d2pa[2]);
  }


     for(ihys=0; ihys<ydo->nohys; ihys++) /* get history variables */
     { 
	  if(ydo->i1ohyt[ihys] == YFLEK) /*Z total kinetic energy */
	  { 
	       Ek=R0;
	       for(i=0; i<ydn->nnopo; i++)
	       {  
		  //  Ek += RP5*(SQR(ydn->d2nvc[0][i]) + SQR(ydn->d2nvc[1][i]) + SQR(ydn->d2nvc[2][i]))*ydn->d1nmct[i];
          //Ek += RP5*(SQR(ydn->d2nvc[0][i]) + SQR(ydn->d2nvc[1][i]) + SQR(ydn->d2nvc[2][i]))*0.0;
          //-PY changed it   the problem is caused by ydn->d1nmct[i]
	       }

 	       /*Z isolated macros to avoid redundant evaluation */
	       ohys = ABS(ydo->d1ohys[ihys]);
	       stprev = MAXIM(EPSILON, ohys);
	       ohyp = R1 - Ek/stprev;
	       ohyp = ABS(ohyp);

	       if(ohyp > ydo->dohyp)
	       { 
		    ydo->d1ohys[ihys] = Ek;
		    ydo->d1ohyt[ihys] = ydc->dctime;
	       } 
	  } 
     }
}


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/



