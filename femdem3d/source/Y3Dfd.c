/*! \file Y3Dfd.c
 *  \brief Y nodal forces
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


/*Z check for negative values - Z - what significance? 
*/
#if (defined(DEBUG) || defined(_DEBUG))
# define DBG_NEG(d,i) if((d)<=R0){ CHRw(stdout,"Negative Jacobian for element "); INTw(stdout,i,5); CHRwcr(stdout); }
#else
# define DBG_NEG(d,i)
#endif


/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/


/* Calculate the Jacobian matrix xj
 */
static INLINE void jacobi(DBL d2xj[NDIME][NDIME], DBL **d2ncc, DBL **d2dsh, INT *i1elto)
{
     INT i, j, k, ip;
     DBL *d1dsh;
     DBL *d1ncc, *d1xj;

     for(i=0; i<NDIME; i++)
     {
	  d1ncc = d2ncc[i];
	  d1xj  = d2xj[i];
    
	  for(j=0; j<NDIME; j++)
	  {
	       d1dsh   = d2dsh[j];
	       d1xj[j] = R0;

	       for(k=0; k<NNODE; k++)
	       {
		    ip = i1elto[k];
		    d1xj[j] += d1ncc[ip]*d1dsh[k];
	       }
	  } 
     }
}


/**********************************************************************/
/**********************************************************************/


/* Calculate the Jacobian matrices xj, xjci
 */
static INLINE void jacobi2(DBL d2xj[NDIME][NDIME], DBL **d2ncc, DBL d2xjci[NDIME][NDIME], DBL **d2nci, DBL **d2dsh, INT *i1elto)
{
     INT i, j, k, ip;
     DBL *d1dsh;
     DBL *d1ncc, *d1xj;
     DBL *d1nci, *d1xjci;

     for(i=0; i<NDIME; i++)
     {
	  d1ncc  = d2ncc[i];
	  d1nci  = d2nci[i];
	  d1xj   = d2xj[i];
	  d1xjci = d2xjci[i];

	  for(j=0; j<NDIME; j++)
	  {
	       d1dsh     = d2dsh[j];
	       d1xj[j]   = R0;
	       d1xjci[j] = R0;

	       for(k=0; k<NNODE; k++)
	       {
		    ip = i1elto[k];
		    d1xj[j]   += d1ncc[ip]*d1dsh[k];
		    d1xjci[j] += d1nci[ip]*d1dsh[k];
	       }
	  } 
     }
}


/**********************************************************************/
/**********************************************************************/


/* setup for gradient calculations
 */
static INLINE void grad(DBL d2nd[NDIME][NNODE], DBL d2xjinv[NDIME][NDIME], DBL **d2dsh)
{
     INT i;

     for(i=0; i<NNODE; i++)
     {
	  d2nd[0][i] = (d2xjinv[0][0]*d2dsh[0][i]) + (d2xjinv[1][0]*d2dsh[1][i]) + (d2xjinv[2][0]*d2dsh[2][i]);
	  d2nd[1][i] = (d2xjinv[0][1]*d2dsh[0][i]) + (d2xjinv[1][1]*d2dsh[1][i]) + (d2xjinv[2][1]*d2dsh[2][i]);
	  d2nd[2][i] = (d2xjinv[0][2]*d2dsh[0][i]) + (d2xjinv[1][2]*d2dsh[1][i]) + (d2xjinv[2][2]*d2dsh[2][i]);
     } 
}


/**********************************************************************/
/**********************************************************************/

/* calculate the Cartesian gradient d2finv
 */
static INLINE void ctgrad(DBL d2nd[NDIME][NNODE], DBL d2finv[NDIME][NDIME], DBL **d2nci, INT *i1elto)
{
     INT i, j, k, ip;
     DBL *d1nd;
     DBL *d1finv, *d1nci;

     /*Z d2nd must be complete 
      */
     for(i=0; i<NDIME; i++) 
     {	
	  d1finv = d2finv[i];
	  d1nci  = d2nci[i];

	  for(j=0; j<NDIME; j++)
	  {
	       d1nd      = d2nd[j];
	       d1finv[j] = R0;

	       for(k=0; k<NNODE; k++)
	       {
		    ip = i1elto[k];
		    d1finv[j] += d1nci[ip]*d1nd[k];
	       }
	  }
     }
}


/**********************************************************************/
/**********************************************************************/


/* calculate the Cartesian gradient d2finv, and velocity gradient L
 */
static INLINE void ctvgrad(DBL d2nd[NDIME][NNODE], DBL d2finv[NDIME][NDIME],
						   DBL **d2nci, DBL d2vgrad[NDIME][NDIME], DBL **d2nvc, INT *i1elto)
{
     INT i, j, k, ip;
     DBL *d1nd;
     DBL *d1finv,  *d1nci;
     DBL *d1vgrad, *d1nvc;

     /*Z d2nd must be complete DBL d2lvg[NDIME][NDIME], 
      */
     for(i=0; i<NDIME; i++)
     {
	  d1finv = d2finv[i];
	  d1nci  = d2nci[i];

	  d1vgrad = d2vgrad[i];
	  d1nvc   = d2nvc[i];

	  for(j=0; j<NDIME; j++)
	  {
	       d1nd       = d2nd[j];
	       d1finv[j]  = R0;
	       d1vgrad[j] = R0;

	       for(k=0; k<NNODE; k++)
	       {
		    ip = i1elto[k];
		    d1finv[j]  += d1nci[ip]*d1nd[k];
		    d1vgrad[j] += d1nvc[ip]*d1nd[k];
	       }
	  }
     }
}


/**********************************************************************/
/**********************************************************************/


/* calculate calculate tensor B, Cauchy stress T, tcs 
 */
static INLINE void cauchy(DBL T[NDIME][NDIME], DBL L[NDIME][NDIME], 
	 DBL d2ftens[NDIME][NDIME], DBL ***d3tcs, DBL dpemu, DBL dpela, DBL detf, DBL detf0, DBL dpeks, INT ielem)
{
     INT i, j;
     DBL scale1, scale2, tweak, D, B;

     scale1 = pow(detf0/(R1/detf),R2/R3);
     scale2 = dpemu/detf0;
     tweak = (dpela*log(detf0)-dpemu)/detf0;

     for(i=0; i<NDIME; i++)
     { 
	  /*Z hopefully the compiler will unroll and expand this entirely */  
	  for(j=0; j<NDIME; j++) 
	  { 
	       /* left Cauchy-Green strain tensor B */
	       B = ((d2ftens[i][0]*d2ftens[j][0]) + (d2ftens[i][1]*d2ftens[j][1]) + (d2ftens[i][2]*d2ftens[j][2])) * scale1;

	       /* deformation rate D (NOTE: [j][i] access forces L to be complete prior) */
	       D = (L[i][j] + L[j][i])*RP5;

	       /* Cauchy stress T */
	       T[i][j] = B*scale2 + D*dpeks;
	  } 

	  /* Cauchy stress along diagonal of T (Neo-Hookean's law) */
	  T[i][i] += tweak;
     }

     for(i=0;i<NDIME;i++)   
     {
	  for(j=0;j<NDIME;j++)
	  {
	       /* net Cauchy stress tcs (T must be complete)*/
	       d3tcs[ielem][i][j] += T[i][j]*WEIGHT;
	  }
     }
}


/**********************************************************************/
/**********************************************************************/


/* calculate surface traction d3st[ig] (tmp1)
 */
static INLINE void stract(DBL d2st[NDIME][NNODE], DBL T[NDIME][NDIME], DBL d2nd[NDIME][NNODE], DBL detjb)
{
     INT i;

     for(i=0; i<NNODE; i++)
     {
	  d2st[0][i] = ((d2nd[0][i]*T[0][0]) + (d2nd[1][i]*T[0][1]) + (d2nd[2][i]*T[0][2]))*detjb;
	  d2st[1][i] = ((d2nd[0][i]*T[1][0]) + (d2nd[1][i]*T[1][1]) + (d2nd[2][i]*T[1][2]))*detjb;
	  d2st[2][i] = ((d2nd[0][i]*T[2][0]) + (d2nd[1][i]*T[2][1]) + (d2nd[2][i]*T[2][2]))*detjb;
     }
}


/**********************************************************************/
/**********************************************************************/


/* calculate nodal force 
 */
static INLINE void nforce(DBL *d1nftx, DBL *d1nfty, DBL *d1nftz, DBL *d1nmct, 
	 DBL *d1lpmm, DBL d3st[NGRSH][NDIME][NNODE], DBL emct, INT *i1elto)
{
     INT i, ip;

     for(i=0; i<NNODE; i++)
     {
	  ip = i1elto[i];
 
	  /*Z nodal total current force */
	  d1nftx[ip] -= WEIGHT*(d3st[1][0][i] + d3st[2][0][i] + d3st[3][0][i] + d3st[4][0][i]); 
	  d1nfty[ip] -= WEIGHT*(d3st[1][1][i] + d3st[2][1][i] + d3st[3][1][i] + d3st[4][1][i]);
	  d1nftz[ip] -= WEIGHT*(d3st[1][2][i] + d3st[2][2][i] + d3st[3][2][i] + d3st[4][2][i]);

	  /*Z nodal mass current translation */
	  d1nmct[ip] += emct*d1lpmm[i]; 
     }
}


/**********************************************************************/
/**********************************************************************/


/* small strain elastic tetrahetra 
 */
static void Yfd3TET10ELS(  
     INT nelem,    INT iprop,
     DBL dpeks,    DBL dpela,    DBL dpemu,    DBL dpero,
     DBL **d2ncc,  DBL **d2nci,  DBL **d2nft,  DBL *d1nmct,  DBL **d2nvc,
     INT *i1elpr,  INT *i1nopr,  INT **i2elto, 
     DBL ***d3dsh, DBL ***d3tcs, DBL *d1lpmm,  DBL *d1emct
     )
{
     DBL d2nd[NDIME][NNODE], d2finv[NDIME][NDIME], d2ftens[NDIME][NDIME];
     DBL d2xj[NDIME][NDIME], d2xjci[NDIME][NDIME], d2xjinv[NDIME][NDIME];
     DBL L[NDIME][NDIME], T[NDIME][NDIME];
     DBL d3st[NGRSH][NDIME][NNODE];
     DBL voli, detci, detj, detf, detf0;
     INT ielem, ig;
     INT *i1elto;

     for(ielem=0; ielem<nelem; ielem++)
     {
	  if(i1elpr[ielem] != iprop) continue; //Z not sure if this would ever happen

	  i1elto = i2elto[ielem];

	  /*Z simpler calculations for ig=0 */
	  ig = 0; 
	  /* calculate Jacobian matrix xj */
	  jacobi(d2xj, d2ncc, d3dsh[ig], i1elto);
	  voli = R0;
	  /* calculate inverse of xj (debug: warn if determinant negative) */
	  YMATINV3(d2xj, d2xjinv, detj); DBG_NEG(detj, ielem)
					      /* calculate Cartesian gradients nd, finv */
					      grad(d2nd, d2xjinv, d3dsh[ig]);
	  ctgrad(d2nd, d2finv, d2nci, i1elto);
	  /* calculate determinant of finv */
	  YMATDET3(d2finv, detf);

	  /*Z store determinant */
	  detf0=R1/detf;

	  /*Z iterate overy Gauss points using deformation (shape) gradient */
	  for(ig=1;ig<NGRSH;ig++)
	  {
	       /*Z full calculations for ig>0 */
	       /* calculate the Jacobian matrixes xj, xjci */
	       jacobi2(d2xj, d2ncc, d2xjci, d2nci, d3dsh[ig], i1elto);
	       /* calculate determinant of xjci */
	       YMATDET3(d2xjci,detci);
	       voli += ABS(detci)*WEIGHT/R6;
	       /* calculate inverse and determinant of xj (debug: warn if negative) */
	       YMATINV3(d2xj,d2xjinv,detj); DBG_NEG(detj,ielem)
						 /* calculate Cartesian gradients nd, finv, and velocity gradient L */
						 grad(d2nd, d2xjinv, d3dsh[ig]);
	       ctvgrad(d2nd, d2finv, d2nci, L, d2nvc, i1elto);
	       /* calculate inverse and determinant of finv */
	       YMATINV3(d2finv, d2ftens, detf);

	       /* calculate left Cauchy-Green strain tensor B, deformation rate D, Cauchy stress T, net stress tcs */
	       cauchy(T, L, d2ftens, d3tcs, dpemu, dpela, detf, detf0, dpeks, ielem);

	       /* end result for ig=1,2,3,4 */
	       stract(d3st[ig], T, d2nd, detj/R6);
	  }

	  /* calculate nodal force */
	  d1emct[ielem] = dpero*voli; //Z property ro - density
	  nforce(d2nft[0], d2nft[1], d2nft[2], d1nmct, d1lpmm, d3st, d1emct[ielem], i1elto);
     }
}

static void Yfd3TET4ELS(  /* small strain elastic 4-noded tetrahedra  */
			  INT    nelem,  INT    iprop,
			  DBL dpeks,    DBL dpela,    DBL dpemu,    DBL dpero,  
			  DBL *d1nccx,  DBL *d1nccy,  DBL *d1nccz,  DBL *d1ncix,  DBL  *d1nciy,
			  DBL *d1nciz,  DBL *d1nfcx,  DBL *d1nfcy,  DBL *d1nfcz,  DBL  *d1nmct,
			  DBL *d1nvcx,  DBL *d1nvcy,  DBL *d1nvcz,  INT *i1nopr,  INT *i1elpr,  INT **i2elto,
			  DBL ***d3tcs,
			  DBL *d1bnap, DBL *d1nvct, INT **i2eltoxb, INT nelem_t,
        DBL dctime,DBL dcrmpt, INT nelemi, INT nnopi, DBL *d1emct,
        DBL *d1ntc, DBL *d1nti, DBL d1ctex, INT icoutf, INT ncstep
			  ) 
{ DBL nx,ny,nz,voli,volc;
  DBL  B[3][3]; /* left Cauchy-Green strain tensor */
  DBL  D[3][3]; /* rate of deformation (stretching) tensor */
  //DBL  E[3][3]; /* strain tensor (small strains) */
  DBL  F[3][3]; /* deformation gradient in global base delta ux/delta x */
  DBL F0[3][3]; /* initial local base */
  DBL FX[3][3]; /* current local base  also delta ux/delta X */
  DBL F0inv[3][3]; /* global base in initial local base */
  DBL FXinv[3][3]; /* global base in current local base */
  DBL  L[3][3]; /* velocity gradient in global base  delta vx/delta x    */
  DBL LX[3][3]; /* velocity gradient in current local base = delta x/delta X */
  DBL  T[3][3]; /* Cauchy stress */
  DBL  detf;
  INT ielem;
  INT i,j,k,l;
  INT *i1elto;

  DBL dredf;   /*ramping factor -added to 3D by AO  */
  INT jnopr,knopr,lnopr;
  DBL Pj,Pk,Pl;
  DBL deltaT; /* Temperature change */
	DBL dups_u; /* [ielem] poisson ratio - undrained  */
	DBL dups; /* [ielem] poisson ratio - drained  */
	DBL ALPHA;

	dups = dpela / (2 * (dpela + dpemu));
	dups_u = dups/0.6;  /* using 40% or lower approaximate at the moment - Devonian Shale; Hydrostone -   */
	
/* Biot coefficient -- input from .y in the future */
	ALPHA=3 * (dups_u - dups)/ ((1 - 2 * dups) * (1 + dups_u));


  /* staged loading */
  if((dctime<dcrmpt)&&(dcrmpt>EPSILON))
  {
 dredf=(dctime)/dcrmpt;
	    }
  else
  { 
dredf=R1;
	}


  for(ielem=0;ielem<nelem;ielem++)
    {
      if(i1elpr[ielem]==iprop)
	{
	  i1elto = i2elto[ielem];
	  deltaT =(d1ntc[i1elto[0]]+d1ntc[i1elto[1]]+d1ntc[i1elto[2]]+d1ntc[i1elto[3]]
	  - d1nti[i1elto[0]]+d1nti[i1elto[1]]+d1nti[i1elto[2]]+d1nti[i1elto[3]])/4.0;
	  for(i=1;i<4;i++)
	    {
		
	      F0[0][i-1]=(d1ncix[i1elto[i]]-d1ncix[i1elto[0]])*(1.0+d1ctex*deltaT); /* init. base */
	      F0[1][i-1]=(d1nciy[i1elto[i]]-d1nciy[i1elto[0]])*(1.0+d1ctex*deltaT);
	      F0[2][i-1]=(d1nciz[i1elto[i]]-d1nciz[i1elto[0]])*(1.0+d1ctex*deltaT);
	      FX[0][i-1]=d1nccx[i1elto[i]]-d1nccx[i1elto[0]]; /* curr. base */
	      FX[1][i-1]=d1nccy[i1elto[i]]-d1nccy[i1elto[0]];
	      FX[2][i-1]=d1nccz[i1elto[i]]-d1nccz[i1elto[0]];
	      LX[0][i-1]=d1nvcx[i1elto[i]]-d1nvcx[i1elto[0]]; /* vel. grad. */
	      LX[1][i-1]=d1nvcy[i1elto[i]]-d1nvcy[i1elto[0]];
	      LX[2][i-1]=d1nvcz[i1elto[i]]-d1nvcz[i1elto[0]];
	    }



	  YMATINV3(F0,F0inv,voli);         /* global base in initial local coordinates    */
	  YMATINV3(FX,FXinv,volc);         /* global base in current local coordinates    */  
	  for(i=0;i<3;i++)
	    {
	      for(j=0;j<3;j++)
		{
		  F[i][j]=R0;
		  L[i][j]=R0;
		  for(k=0;k<3;k++)
		    {
		      F[i][j]=F[i][j]+FX[i][k]*F0inv[k][j]; /* deformation gradient  */
		      L[i][j]=L[i][j]+LX[i][k]*FXinv[k][j]; /* velocity gradient     */
		    }
		}
	    }
	  for(i=0;i<3;i++)
	    {
	      for(j=0;j<3;j++)
		{
		  B[i][j]=R0;
		  for(k=0;k<3;k++)
		    {
		      B[i][j]=B[i][j]+F[i][k]*F[j][k];  /* left Cauchy-Green strain  */
		    }
		  D[i][j]=RP5*(L[i][j]+L[j][i]);      /* rate of deformation       */
		}
	    }


	  /*Calculate Cauchy stress*/
	  detf=volc/voli; 

	  for(i=0;i<3;i++)   
	    {
	      for(j=0;j<3;j++)
		{ 
		  T[i][j]=(dpemu/detf)*B[i][j]+dpeks*D[i][j];
		} 
	      T[i][i]=T[i][i]+(dpela*log(detf)-dpemu)/detf;   //+ALPHA*d1npore[ielem]  //-ao Asiri added for 'law of effective streses -pore pressure' FEMDEM +, theory -

	    }
	  for(i=0;i<3;i++)   
	    {
	      for(j=0;j<3;j++)
		{
		  d3tcs[ielem][i][j]=T[i][j];
		}
	    }

        d1emct[ielem]=dpero*voli/6.0;
      for(i=0;i<4;i++)      /* Nodal Forces */
      { 
		j=i+1; 
		if(j>3)j=0;
        k=j+1; 
		if(k>3)k=0;
        l=k+1; 
		if(l>3)l=0;
        nx=((d1nccy[i1elto[k]]-d1nccy[i1elto[j]])*
            (d1nccz[i1elto[l]]-d1nccz[i1elto[j]])-
            (d1nccy[i1elto[l]]-d1nccy[i1elto[j]])*
            (d1nccz[i1elto[k]]-d1nccz[i1elto[j]]))/R6;
        ny=((d1nccz[i1elto[k]]-d1nccz[i1elto[j]])*
            (d1nccx[i1elto[l]]-d1nccx[i1elto[j]])-
            (d1nccx[i1elto[k]]-d1nccx[i1elto[j]])*
            (d1nccz[i1elto[l]]-d1nccz[i1elto[j]]))/R6; 
        nz=((d1nccx[i1elto[k]]-d1nccx[i1elto[j]])*
            (d1nccy[i1elto[l]]-d1nccy[i1elto[j]])-
            (d1nccy[i1elto[k]]-d1nccy[i1elto[j]])*
            (d1nccx[i1elto[l]]-d1nccx[i1elto[j]]))/R6;

        d1nmct[i1elto[i]]=d1nmct[i1elto[i]]+dpero*voli/24;
        d1nvct[i1elto[i]]=d1nvct[i1elto[i]]+voli/24;


        if((i==0)||(i==2))
        { 
          d1nfcx[i1elto[i]]=d1nfcx[i1elto[i]]+
                                    (T[0][0]*nx+T[0][1]*ny+T[0][2]*nz);
          d1nfcy[i1elto[i]]=d1nfcy[i1elto[i]]+
                                    (T[1][0]*nx+T[1][1]*ny+T[1][2]*nz);
          d1nfcz[i1elto[i]]=d1nfcz[i1elto[i]]+
                                    (T[2][0]*nx+T[2][1]*ny+T[2][2]*nz); 
        }
        else
        { 
	  d1nfcx[i1elto[i]]=d1nfcx[i1elto[i]]-
                                    (T[0][0]*nx+T[0][1]*ny+T[0][2]*nz);
          d1nfcy[i1elto[i]]=d1nfcy[i1elto[i]]-
                                    (T[1][0]*nx+T[1][1]*ny+T[1][2]*nz);
          d1nfcz[i1elto[i]]=d1nfcz[i1elto[i]]-
                                    (T[2][0]*nx+T[2][1]*ny+T[2][2]*nz); 
        }







            if((i2eltoxb[ielem][i]>=0)&&(i2eltoxb[ielem][i]<nelem_t))
          {
            Pj=dredf*d1bnap[i2eltoxb[ielem][i]];   //added dredf* - ao   //this just needs to be from surface nodes to disc nodes
            Pk=dredf*d1bnap[i2eltoxb[ielem][i]];
            Pl=dredf*d1bnap[i2eltoxb[ielem][i]];
            if((i==0)||(i==2))
              {
                d1nfcx[i1elto[j]]=d1nfcx[i1elto[j]]-Pj*nx;
                d1nfcy[i1elto[j]]=d1nfcy[i1elto[j]]-Pj*ny;
                d1nfcz[i1elto[j]]=d1nfcz[i1elto[j]]-Pj*nz;

                d1nfcx[i1elto[k]]=d1nfcx[i1elto[k]]-Pk*nx;
                d1nfcy[i1elto[k]]=d1nfcy[i1elto[k]]-Pk*ny;
                d1nfcz[i1elto[k]]=d1nfcz[i1elto[k]]-Pk*nz;

                d1nfcx[i1elto[l]]=d1nfcx[i1elto[l]]-Pl*nx;
                d1nfcy[i1elto[l]]=d1nfcy[i1elto[l]]-Pl*ny;
                d1nfcz[i1elto[l]]=d1nfcz[i1elto[l]]-Pl*nz;
              }
            else
              {
                d1nfcx[i1elto[j]]=d1nfcx[i1elto[j]]+Pj*nx;
                d1nfcy[i1elto[j]]=d1nfcy[i1elto[j]]+Pj*ny;
                d1nfcz[i1elto[j]]=d1nfcz[i1elto[j]]+Pj*nz;

                d1nfcx[i1elto[k]]=d1nfcx[i1elto[k]]+Pk*nx;
                d1nfcy[i1elto[k]]=d1nfcy[i1elto[k]]+Pk*ny;
                d1nfcz[i1elto[k]]=d1nfcz[i1elto[k]]+Pk*nz;

                d1nfcx[i1elto[l]]=d1nfcx[i1elto[l]]+Pl*nx;
                d1nfcy[i1elto[l]]=d1nfcy[i1elto[l]]+Pl*ny;
                d1nfcz[i1elto[l]]=d1nfcz[i1elto[l]]+Pl*nz;
              }
          }

	  }
    }
 }


}

static void Yfd3TET4JOINT(		/* joint element for 4-node tetrahedra */
     INT nelem,
     INT iprop,
      DBL dpeft, DBL dpegfn, DBL dpegfs,
     DBL dpeks, DBL dpepe,
     DBL *d1nccx,DBL *d1nccy,DBL *d1nccz,
     DBL *d1nftx,DBL *d1nfty,DBL *d1nftz,
     DBL *d1nvcx,DBL *d1nvcy,DBL *d1nvcz,
     DBL *d1sdel,
     INT *i1elpr,INT **i2elto,INT *i1eljo,
//-PY changed it for 3D_fracture_coupling_with_multiphase 
     INT **i2eljp,
     INT *i0iecff, INT *i1iecn, INT *i1elcf,
     INT *i1iect, INT *n0icoup, DBL *d0iedi,
     INT nneigh, INT **i2nnei, INT **i2elfr, INT *i1elbe,
     INT *i1nei, DBL dpefr, INT *i1elcft, INT ncstep, INT nelemi, INT nnopi,
     DBL dpcoh, DBL dpicf, INT *i1elty, INT *icoutf
     )
{
     DBL dpefa=0.63;
     DBL dpefb=1.8;
     DBL dpefc=6.0;
     DBL dpefm=0.0;
     DBL small,odis,sdis,sabs,o[3],s[3],op,sp,ot,st,z,sigma,tau;
     DBL h,area,el;
     INT ielem,integ,i0,i1,i2,i3,i4,i5,nfail,nfail_n,nfail_s,nnfail;
//     INT nsoft;
     INT *i1elto;

     DBL midx[3],midy[3],midz[3];   /*coordinates of middle points; 0-05, 1-14, 2-23*/
     DBL nx,ny,nz;   /*unit normal vector of middle plane*/
//     DBL jl[3];   /*distance of each pair of nodes: 0-05, 1-14, 2-23*/
     DBL prjx[3],prjy[3],prjz[3];   /*projection of first 3 joint element nodes on middle plane*/
     DBL sx[3],sy[3],sz[3];   /*unit shear direction on middle plane*/

     DBL mpx[6],mpy[6],mpz[6];   /*middle points of the 6 edges: 0-01, 1-12, 2-20, 3-54, 4-43, 5-35*/
     DBL dpefs;
//     FILE *fp;
//     fp=fopen("Yfd3TET4JOINT_output.txt","w+");
//     FILE *fp1;
//     fp1=fopen("joint_stress.txt","a+");

//     DBL fx[6],fy[6],fz[6];
//     INT i;


//-PY changed it for 3D_fracture_coupling_with_multiphase
     INT icouco,icouta;         /* contact element and target element */
     INT ifem,jfem;             /* neighbouring finite elements of broken joint element */
     INT icoup,inei,ijno;

     INT i;
     DBL sigma_tmp;

     small=EPSILON;
//     nsoft=0;

/*     for(ielem=0;ielem<nelem;ielem++)
       {
	 fprintf(fp1,"ielem%ld %ld\n",ielem,i1eljo[ielem]);
       }
*/
     for(ielem=0;ielem<nelem;ielem++)
     {
	  if(i1elpr[ielem]==iprop)
	  {
/*	    for(i=0;i<6;i++)
	      {
		fx[i]=R0;
		fy[i]=R0;
		fz[i]=R0;
	      }
*/
               i1elto=i2elto[ielem];

	       i0=i1elto[0];
	       i1=i1elto[1];
	       i2=i1elto[2];
	       i3=i1elto[3];
	       i4=i1elto[4];
	       i5=i1elto[5];

//               fprintf(fp,"ielem=%ld\n",ielem);
//	       fprintf(fp,"0=%ld\t1=%ld\t2=%ld\t3=%ld\t4=%ld\t5=%ld\n",i0,i1,i2,i3,i4,i5);

	       midx[0]=RP5*(d1nccx[i0]+d1nccx[i5]);
	       midy[0]=RP5*(d1nccy[i0]+d1nccy[i5]);
	       midz[0]=RP5*(d1nccz[i0]+d1nccz[i5]);
	       midx[1]=RP5*(d1nccx[i1]+d1nccx[i4]);
	       midy[1]=RP5*(d1nccy[i1]+d1nccy[i4]);
	       midz[1]=RP5*(d1nccz[i1]+d1nccz[i4]);
	       midx[2]=RP5*(d1nccx[i2]+d1nccx[i3]);
	       midy[2]=RP5*(d1nccy[i2]+d1nccy[i3]);
	       midz[2]=RP5*(d1nccz[i2]+d1nccz[i3]);

	       mpx[0]=RP5*(d1nccx[i0]+d1nccx[i1]);
	       mpy[0]=RP5*(d1nccy[i0]+d1nccy[i1]);
	       mpz[0]=RP5*(d1nccz[i0]+d1nccz[i1]);
	       mpx[1]=RP5*(d1nccx[i1]+d1nccx[i2]);
	       mpy[1]=RP5*(d1nccy[i1]+d1nccy[i2]);
	       mpz[1]=RP5*(d1nccz[i1]+d1nccz[i2]);
	       mpx[2]=RP5*(d1nccx[i2]+d1nccx[i0]);
	       mpy[2]=RP5*(d1nccy[i2]+d1nccy[i0]);
	       mpz[2]=RP5*(d1nccz[i2]+d1nccz[i0]);
	       mpx[3]=RP5*(d1nccx[i5]+d1nccx[i4]);
	       mpy[3]=RP5*(d1nccy[i5]+d1nccy[i4]);
	       mpz[3]=RP5*(d1nccz[i5]+d1nccz[i4]);
	       mpx[4]=RP5*(d1nccx[i4]+d1nccx[i3]);
	       mpy[4]=RP5*(d1nccy[i4]+d1nccy[i3]);
	       mpz[4]=RP5*(d1nccz[i4]+d1nccz[i3]);
	       mpx[5]=RP5*(d1nccx[i3]+d1nccx[i5]);
	       mpy[5]=RP5*(d1nccy[i3]+d1nccy[i5]);
	       mpz[5]=RP5*(d1nccz[i3]+d1nccz[i5]);

//	       fprintf(fp,"X0[%ld]=%f\tY0[%ld]=%f\tZ0[%ld]=%f\nX1[%ld]=%f\tY1[%ld]=%f\tZ1[%ld]=%f\nX2[%ld]=%f\tY2[%ld]=%f\tZ2[%ld]=%f\nX3[%ld]=%f\tY3[%ld]=%f\tZ3[%ld]=%f\nX4[%ld]=%f\tY4[%ld]=%f\tZ4[%ld]=%f\nX5[%ld]=%f\tY5[%ld]=%f\tZ5[%ld]=%f\n",i0,d1nccx[i0],i0,d1nccy[i0],i0,d1nccz[i0],i1,d1nccx[i1],i1,d1nccy[i1],i1,d1nccz[i1],i2,d1nccx[i2],i2,d1nccy[i2],i2,d1nccz[i2],i3,d1nccx[i3],i3,d1nccy[i3],i3,d1nccz[i3],i4,d1nccx[i4],i4,d1nccy[i4],i4,d1nccz[i4],i5,d1nccx[i5],i5,d1nccy[i5],i5,d1nccz[i5]);
//	       fprintf(fp,"midx[0]=%f\tmidy[0]=%f\tmidz[0]=%f\nmidx[1]=%f\tmidy[1]=%f\tmidz[1]=%f\nmidx[2]=%f\tmidy[2]=%f\tmidz[2]=%f\n",midx[0],midy[0],midz[0],midx[1],midy[1],midz[1],midx[2],midy[2],midz[2]);

               /* 2*area of this triangle */
	       h=SQRT(((midy[1]-midy[0])*(midz[2]-midz[0])-(midz[1]-midz[0])*(midy[2]-midy[0]))
		      *((midy[1]-midy[0])*(midz[2]-midz[0])-(midz[1]-midz[0])*(midy[2]-midy[0]))
		      +((midz[1]-midz[0])*(midx[2]-midx[0])-(midx[1]-midx[0])*(midz[2]-midz[0]))
		      *((midz[1]-midz[0])*(midx[2]-midx[0])-(midx[1]-midx[0])*(midz[2]-midz[0]))
		      +((midx[1]-midx[0])*(midy[2]-midy[0])-(midy[1]-midy[0])*(midx[2]-midx[0]))
		      *((midx[1]-midx[0])*(midy[2]-midy[0])-(midy[1]-midy[0])*(midx[2]-midx[0])));

	       /* in numerical integration, the area used in the calculation of node force is 1/6 of the triangle area */
	       area=h/R12;

               /* nx,ny,nz can be positive or negative */
               if((i1eljo[ielem]==0)||(i1eljo[ielem]==2))
		 {
		   nx=((midy[1]-midy[0])*(midz[2]-midz[0])-(midz[1]-midz[0])*(midy[2]-midy[0]))/(h+small);
		   ny=((midz[1]-midz[0])*(midx[2]-midx[0])-(midx[1]-midx[0])*(midz[2]-midz[0]))/(h+small);
		   nz=((midx[1]-midx[0])*(midy[2]-midy[0])-(midy[1]-midy[0])*(midx[2]-midx[0]))/(h+small);
		 }
               else
		 {
		   nx=((midy[2]-midy[0])*(midz[1]-midz[0])-(midz[2]-midz[0])*(midy[1]-midy[0]))/(h+small);
		   ny=((midz[2]-midz[0])*(midx[1]-midx[0])-(midx[2]-midx[0])*(midz[1]-midz[0]))/(h+small);
		   nz=((midx[2]-midx[0])*(midy[1]-midy[0])-(midy[2]-midy[0])*(midx[1]-midx[0]))/(h+small);
		 }

//               fprintf(fp,"area=%f\n",area);

	       /* average length of element edge */
	       el=(SQRT((midx[0]-midx[1])*(midx[0]-midx[1])+(midy[0]-midy[1])*(midy[0]-midy[1])
		       +(midz[0]-midz[1])*(midz[0]-midz[1]))
		  +SQRT((midx[1]-midx[2])*(midx[1]-midx[2])+(midy[1]-midy[2])*(midy[1]-midy[2])
			+(midz[1]-midz[2])*(midz[1]-midz[2]))
		  +SQRT((midx[2]-midx[0])*(midx[2]-midx[0])+(midy[2]-midy[0])*(midy[2]-midy[0])
			+(midz[2]-midz[0])*(midz[2]-midz[0])))/R3;

//	       fprintf(fp,"average edge h=%f\n",h);

//	       fprintf(fp,"nx=%f\tny=%f\tnz=%f\n",nx,ny,nz);

/*	       jl[0]=SQRT((d1nccx[i5]-d1nccx[i0])*(d1nccx[i5]-d1nccx[i0])+
			(d1nccy[i5]-d1nccy[i0])*(d1nccy[i5]-d1nccy[i0])+
			(d1nccz[i5]-d1nccz[i0])*(d1nccz[i5]-d1nccz[i0]));
	       jl[1]=SQRT((d1nccx[i4]-d1nccx[i1])*(d1nccx[i4]-d1nccx[i1])+
			(d1nccy[i4]-d1nccy[i1])*(d1nccy[i4]-d1nccy[i1])+
			(d1nccz[i4]-d1nccz[i1])*(d1nccz[i4]-d1nccz[i1]));
	       jl[2]=SQRT((d1nccx[i3]-d1nccx[i2])*(d1nccx[i3]-d1nccx[i2])+
			(d1nccy[i3]-d1nccy[i2])*(d1nccy[i3]-d1nccy[i2])+
			(d1nccz[i3]-d1nccz[i2])*(d1nccz[i3]-d1nccz[i2]));

	       fprintf(fp,"jl[0]=%f\tjl[1]=%f\tjl[2]=%f\n",jl[0],jl[1],jl[2]);
*/
           o[0]=nx*(mpx[0]-mpx[3])+ny*(mpy[0]-mpy[3])+nz*(mpz[0]-mpz[3]);
           o[1]=nx*(mpx[1]-mpx[4])+ny*(mpy[1]-mpy[4])+nz*(mpz[1]-mpz[4]);
           o[2]=nx*(mpx[2]-mpx[5])+ny*(mpy[2]-mpy[5])+nz*(mpz[2]-mpz[5]);

		  prjx[0]=mpx[0]-o[0]*nx; prjy[0]=mpy[0]-o[0]*ny; prjz[0]=mpz[0]-o[0]*nz;
		  prjx[1]=mpx[1]-o[1]*nx; prjy[1]=mpy[1]-o[1]*ny; prjz[1]=mpz[1]-o[1]*nz;
		  prjx[2]=mpx[2]-o[2]*nx; prjy[2]=mpy[2]-o[2]*ny; prjz[2]=mpz[2]-o[2]*nz;

//	       fprintf(fp,"prjx[0]=%f\tprjy[0]=%f\tprjz[0]=%f\nprjx[1]=%f\tprjy[1]=%f\tprjz[1]=%f\nprjx[2]=%f\tprjy[2]=%f\tprjz[2]=%f\n",prjx[0],prjy[0],prjz[0],prjx[1],prjy[1],prjz[1],prjx[2],prjy[2],prjz[2]);

	       /* sx[i],sy[i],sz[i] can be positive or negative */
               if(DABS(prjx[0]-mpx[3])<EPSILON)
                 {
                   sx[0]=R0;
                 }
               else
                 {
                   sx[0]=(prjx[0]-mpx[3])/(SQRT((prjx[0]-mpx[3])*(prjx[0]-mpx[3])+(prjy[0]-mpy[3])*(prjy[0]-mpy[3])
                                                     +(prjz[0]-mpz[3])*(prjz[0]-mpz[3])));
                 }
               if(DABS(prjy[0]-mpy[3])<EPSILON)
                 {
                   sy[0]=R0;
                 }
               else
                 {
                   sy[0]=(prjy[0]-mpy[3])/(SQRT((prjx[0]-mpx[3])*(prjx[0]-mpx[3])+(prjy[0]-mpy[3])*(prjy[0]-mpy[3])
                                                     +(prjz[0]-mpz[3])*(prjz[0]-mpz[3])));
                 }
               if(DABS(prjz[0]-mpz[3])<EPSILON)
                 {
                   sz[0]=R0;
                 }
               else
                 {
                   sz[0]=(prjz[0]-mpz[3])/(SQRT((prjx[0]-mpx[3])*(prjx[0]-mpx[3])+(prjy[0]-mpy[3])*(prjy[0]-mpy[3])
                                                     +(prjz[0]-mpz[3])*(prjz[0]-mpz[3])));
                 }

               if(DABS(prjx[1]-mpx[4])<EPSILON)
                 {
                   sx[1]=R0;
                 }
               else
                 {
                   sx[1]=(prjx[1]-mpx[4])/(SQRT((prjx[1]-mpx[4])*(prjx[1]-mpx[4])+(prjy[1]-mpy[4])*(prjy[1]-mpy[4])
                                                     +(prjz[1]-mpz[4])*(prjz[1]-mpz[4])));
                 }
               if(DABS(prjy[1]-mpy[4])<EPSILON)
                 {
                   sy[1]=R0;
                 }
               else
                 {
                   sy[1]=(prjy[1]-mpy[4])/(SQRT((prjx[1]-mpx[4])*(prjx[1]-mpx[4])+(prjy[1]-mpy[4])*(prjy[1]-mpy[4])
                                                     +(prjz[1]-mpz[4])*(prjz[1]-mpz[4])));
                 }
               if(DABS(prjz[1]-mpz[4])<EPSILON)
                 {
                   sz[1]=R0;
                 }
               else
                 {
                   sz[1]=(prjz[1]-mpz[4])/(SQRT((prjx[1]-mpx[4])*(prjx[1]-mpx[4])+(prjy[1]-mpy[4])*(prjy[1]-mpy[4])
                                                     +(prjz[1]-mpz[4])*(prjz[1]-mpz[4])));
                 }

               if(DABS(prjx[2]-mpx[5])<EPSILON)
                 {
                   sx[2]=R0;
                 }
               else
                 {
                   sx[2]=(prjx[2]-mpx[5])/(SQRT((prjx[2]-mpx[5])*(prjx[2]-mpx[5])+(prjy[2]-mpy[5])*(prjy[2]-mpy[5])
                                                     +(prjz[2]-mpz[5])*(prjz[2]-mpz[5])));
                 }
               if(DABS(prjy[2]-mpy[5])<EPSILON)
                 {
                   sy[2]=R0;
                 }
               else
                 {
                   sy[2]=(prjy[2]-mpy[5])/(SQRT((prjx[2]-mpx[5])*(prjx[2]-mpx[5])+(prjy[2]-mpy[5])*(prjy[2]-mpy[5])
                                                     +(prjz[2]-mpz[5])*(prjz[2]-mpz[5])));
                 }
               if(DABS(prjz[2]-mpz[5])<EPSILON)
                 {
                   sz[2]=R0;
                 }
               else
                 {
                   sz[2]=(prjz[2]-mpz[5])/(SQRT((prjx[2]-mpx[5])*(prjx[2]-mpx[5])+(prjy[2]-mpy[5])*(prjy[2]-mpy[5])
                                                     +(prjz[2]-mpz[5])*(prjz[2]-mpz[5])));
                 }

	       s[0]=sx[0]*(mpx[0]-mpx[3])+sy[0]*(mpy[0]-mpy[3])+sz[0]*(mpz[0]-mpz[3]);
	       s[1]=sx[1]*(mpx[1]-mpx[4])+sy[1]*(mpy[1]-mpy[4])+sz[1]*(mpz[1]-mpz[4]);
	       s[2]=sx[2]*(mpx[2]-mpx[5])+sy[2]*(mpy[2]-mpy[5])+sz[2]*(mpz[2]-mpz[5]);
/*
               s[0]=sx[0]*DABS(mpx[3]-mpx[0])+sy[0]*DABS(mpy[3]-mpy[0])+sz[0]*DABS(mpz[3]-mpz[0]);
               s[1]=sx[1]*DABS(mpx[4]-mpx[1])+sy[1]*DABS(mpy[4]-mpy[1])+sz[1]*DABS(mpz[4]-mpz[1]);
               s[2]=sx[2]*DABS(mpx[5]-mpx[2])+sy[2]*DABS(mpy[5]-mpy[2])+sz[2]*DABS(mpz[5]-mpz[2]);

	       o[0]=nx*(mpx[0]-mpx[3])+ny*(mpy[0]-mpy[3])+nz*(mpz[0]-mpz[3]);
	       o[1]=nx*(mpx[1]-mpx[4])+ny*(mpy[1]-mpy[4])+nz*(mpz[1]-mpz[4]);
	       o[2]=nx*(mpx[2]-mpx[5])+ny*(mpy[2]-mpy[5])+nz*(mpz[2]-mpz[5]);
*/
//	       fprintf(fp,"s[0]=%f\ts[1]=%f\ts[2]=%f\n",s[0],s[1],s[2]);
//	       fprintf(fp,"o[0]=%f\to[1]=%f\to[2]=%f\n",o[0],o[1],o[2]);


//	       dpefs=2*dpeft; //-problem is that this is 0, not being assigned any value

	       op=R2*el*dpeft/dpepe;
	       //sp=R2*el*dpefs/dpepe;
	       ot=MAXIM((R2*op),(R3*dpegfn/dpeft));   /*need further investigation*/
	       //st=MAXIM((R2*sp),(R3*dpegf/dpefs));   /*need further investigation*/

	       nfail=0;
               nfail_n=0;
               nfail_s=0;
               nnfail=0;
	       //numerical integration: three integration points
	       for(integ=0;integ<3;integ++)
	       {
/*		    if(integ==0)
		    {
			 odis=RP5*(o[0]+o[1]);
			 sdis=RP5*(s[0]+s[1]);
		    }
		    else if(integ==1)
		    {
			 odis=RP5*(o[1]+o[2]);
			 sdis=RP5*(s[1]+s[2]);
		    }
		    else
		    {
			 odis=RP5*(o[2]+o[0]);
			 sdis=RP5*(s[2]+s[0]);
		    }
*/
		//    dpefs=d2elfs[ielem][integ];
               /* Mohr-Coulomb failure criterion with the tension cut-off */
                odis=o[integ];
               sigma_tmp=EPSILON;
               if(odis<R0)
               {
                   sigma_tmp=R2*odis*dpeft/op;
               }
               if(sigma_tmp>R0)
               {
                   dpefs=dpcoh;
               }
               else
               {
                   dpefs=dpcoh-dpicf*sigma_tmp;
               }
               
		    sp=R2*el*dpefs/dpepe;
		    st=MAXIM((R2*sp),(R3*dpegfs/dpefs));   /*need further investigation*/

       		    odis=o[integ];
		    sdis=s[integ];
			sabs=DABS(sdis);
   		    if((odis>op)&&(sabs>sp))
		    {
			 z=SQRT(((odis-op)/ot)*((odis-op)/ot)+((sabs-sp)/st)*((sabs-sp)/st));
                         if(z>=R1)
                         {
                         nfail_n=nfail_n+1;
                         nfail_s=nfail_s+1;
                         nnfail=nnfail+1;

                         }

		    }
		    else if(odis>op)
		    {
			 z=(odis-op)/ot;
                         if(z>=R1)
                         {
                         nfail_n=nfail_n+1;
//                         nfail_s=nfail_s+1;
                         }

		    }
		    else if(sabs>sp)
		    {
			 z=(sabs-sp)/st;
                         if(z>=R1)
                         {
//                         nfail_n=nfail_n+1;
                         nfail_s=nfail_s+1;
                         }

		    }
		    else
		    {
			 z=R0;
		    }

		    if(z>=R1)
		    {
			 nfail=nfail+1;
			 if((nfail>1)&&(i1elpr[ielem]>=0))
			 {
			      i1elpr[ielem]=iprop-YIPROPMAX;
                              if(nnfail>1)
                              {
                                 i1elty[ielem]=5;
                                 }
                               else if ((nfail_s>1)&&(nfail_n<1))
                              {
                                 i1elty[ielem]=3; //normal
                                 }
                               else if ((nfail_s<1)&&(nfail_n>1))
                              { 
                                 i1elty[ielem]=4; //shear
                                 }
                               else
                                {
                                i1elty[ielem]=5; //shear-rotation
				}

//-PY changed it for 3D_fracture_coupling_with_multiphase -ao commented this again
               if(odis<R0)
			     {
			       i1elcft[i2eljp[ielem][0]]=ncstep;
			       i1elcft[i2eljp[ielem][1]]=ncstep;
			     }

			   for(i=0;i<4;i++)
			     {
			       if(i2eljp[i2eljp[ielem][0]][i]==ielem)
				 {
				   i2eljp[i2eljp[ielem][0]][i]=-1;
				 }
			       if(i2eljp[i2eljp[ielem][1]][i]==ielem)
				 {
				   i2eljp[i2eljp[ielem][1]][i]=-1;
				 }
			     }

			   i1elbe[i2eljp[ielem][0]]=YIPROPMAX;
			   i1elbe[i2eljp[ielem][1]]=YIPROPMAX;

                              // creat new contact couple for other fem elements node-connected with the broken joint
                              for(ijno=0;ijno<3;ijno++)
                                {
                                  for(inei=0;inei<i1nei[i1elto[ijno]];inei++)
                                    {
                                      if(i2nnei[i1elto[ijno]][inei]>=0)
                                        {
                                          jfem=i2nnei[i1elto[ijno]][inei];
					  i1elbe[jfem]=YIPROPMAX;
                                        }
                                    }
                                }

                              (*d0iedi)=2000.0;
                              //(*icoutf)=1;



			 }
			 z=R1;
		    }
		    z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*
		       exp(z*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))
			 *(dpefa*(R1-z)+dpefb*pow((R1-z),dpefc));


		    if(odis<R0)
		    {
			 sigma=R2*odis*dpeft/op;
		    }
		    else if(odis>op)
		    {
			 sigma=dpeft*z;
//			 nsoft=nsoft+1;
		    }
		    else
		    {
			 sigma=(R2*odis/op-(odis/op)*(odis/op))*z*dpeft;
		    }

                    /* Mohr-Coulomb failure criterion with the tension cut-off */
/*                    if(sigma>dpeft)
                      {
                        dpefs=dpcoh-dpicf*dpeft;
                      }
                    else
                      {
                        dpefs=dpcoh-dpicf*sigma;
                      }

		    d2elfs[ielem][integ]=dpefs;
*/
//		    fprintf(fp,"integ=%ld\n",integ);
//		    fprintf(fp,"sigma=%f\ttau=%f z=%f dpefs=%f R2=%d \n",sigma,tau,z,dpefs, R2);
//		    fprintf(fp,"dpefm=%f sigma=%f sabs=%f sp=%f \n",dpefm,sigma, sabs, sp);

		    if((sigma>R0)&&(sabs>sp))
		    {

//			fprintf(fp,"1 \n\n");

			 tau=z*dpefs;
		    }
		    else if(sigma>R0)
		    {
//			fprintf(fp,"2 \n\n");

			tau=(R2*(sabs/sp)-(sabs/sp)*(sabs/sp))*z*dpefs;
		    }
		    else if(sabs>sp)
		    {
//			fprintf(fp,"3 \n\n");

			tau=z*dpefs-dpefm*sigma;
		    }
		    else
		    {
//			fprintf(fp,"4 \n\n");

			tau=(R2*(sabs/sp)-(sabs/sp)*(sabs/sp))*(z*dpefs-dpefm*sigma); //-ao getting NaNs sometimes

//			fprintf(fp, "tau=%f \n", tau);
		    }
		    		    
		    if(sdis<R0)
		    {
//			fprintf(fp,"5 \n\n");

			 tau=-tau;
		    }
		    
//		    fprintf(fp,"integ=%ld\n",integ);
//		    fprintf(fp,"sigma=%f\ttau=%f z=%f dpefs=%f R2=%d pepe=%e \n",sigma,tau,z,dpefs, R2, dpepe);
//		    fprintf(fp,"dpefm=%f sigma=%f sabs=%f sp=%f \n",dpefm,sigma, sabs, sp);

		    //if(ielem==96295)
                    //if(ielem==114775)
		    //if(ielem==110133)
                    //if((ielem==110133)&&(integ==2))
                    //if(ielem==2)
                    //if((ielem==3)&&(integ==2))
		      //{
			//fprintf(fp1,"ielem=%ld\tinteg=%ld\n",ielem,integ);
			//fprintf(fp1,"o=%e\ts=%e\n",odis,sdis);
			//fprintf(fp1,"sigma=%e\ttau=%e\n",sigma,tau);
                        //fprintf(fp1,"%e %e %e %e ",odis,sdis,sigma,tau);
		      //}

		    if(integ==0)
		    {
			 d1nftx[i0]=d1nftx[i0]-area*(tau*sx[0]+sigma*nx);
			 d1nftx[i5]=d1nftx[i5]+area*(tau*sx[0]+sigma*nx);
			 d1nfty[i0]=d1nfty[i0]-area*(tau*sy[0]+sigma*ny);
			 d1nfty[i5]=d1nfty[i5]+area*(tau*sy[0]+sigma*ny);
			 d1nftz[i0]=d1nftz[i0]-area*(tau*sz[0]+sigma*nz);
			 d1nftz[i5]=d1nftz[i5]+area*(tau*sz[0]+sigma*nz);
			 d1nftx[i1]=d1nftx[i1]-area*(tau*sx[0]+sigma*nx);
			 d1nftx[i4]=d1nftx[i4]+area*(tau*sx[0]+sigma*nx);
			 d1nfty[i1]=d1nfty[i1]-area*(tau*sy[0]+sigma*ny);
			 d1nfty[i4]=d1nfty[i4]+area*(tau*sy[0]+sigma*ny);
			 d1nftz[i1]=d1nftz[i1]-area*(tau*sz[0]+sigma*nz);
			 d1nftz[i4]=d1nftz[i4]+area*(tau*sz[0]+sigma*nz);
		    }
		    else if(integ==1)
		    {
			 d1nftx[i1]=d1nftx[i1]-area*(tau*sx[1]+sigma*nx);
			 d1nftx[i4]=d1nftx[i4]+area*(tau*sx[1]+sigma*nx);
			 d1nfty[i1]=d1nfty[i1]-area*(tau*sy[1]+sigma*ny);
			 d1nfty[i4]=d1nfty[i4]+area*(tau*sy[1]+sigma*ny);
			 d1nftz[i1]=d1nftz[i1]-area*(tau*sz[1]+sigma*nz);
			 d1nftz[i4]=d1nftz[i4]+area*(tau*sz[1]+sigma*nz);
			 d1nftx[i2]=d1nftx[i2]-area*(tau*sx[1]+sigma*nx);
			 d1nftx[i3]=d1nftx[i3]+area*(tau*sx[1]+sigma*nx);
			 d1nfty[i2]=d1nfty[i2]-area*(tau*sy[1]+sigma*ny);
			 d1nfty[i3]=d1nfty[i3]+area*(tau*sy[1]+sigma*ny);
			 d1nftz[i2]=d1nftz[i2]-area*(tau*sz[1]+sigma*nz);
			 d1nftz[i3]=d1nftz[i3]+area*(tau*sz[1]+sigma*nz);
		    }
		    else
		    {
			 d1nftx[i2]=d1nftx[i2]-area*(tau*sx[2]+sigma*nx);
			 d1nftx[i3]=d1nftx[i3]+area*(tau*sx[2]+sigma*nx);
			 d1nfty[i2]=d1nfty[i2]-area*(tau*sy[2]+sigma*ny);
			 d1nfty[i3]=d1nfty[i3]+area*(tau*sy[2]+sigma*ny);
			 d1nftz[i2]=d1nftz[i2]-area*(tau*sz[2]+sigma*nz);
			 d1nftz[i3]=d1nftz[i3]+area*(tau*sz[2]+sigma*nz);
			 d1nftx[i0]=d1nftx[i0]-area*(tau*sx[2]+sigma*nx);
			 d1nftx[i5]=d1nftx[i5]+area*(tau*sx[2]+sigma*nx);
			 d1nfty[i0]=d1nfty[i0]-area*(tau*sy[2]+sigma*ny);
			 d1nfty[i5]=d1nfty[i5]+area*(tau*sy[2]+sigma*ny);
			 d1nftz[i0]=d1nftz[i0]-area*(tau*sz[2]+sigma*nz);
			 d1nftz[i5]=d1nftz[i5]+area*(tau*sz[2]+sigma*nz);
		    }
/*
		    if(integ==0)
		    {
			 fx[0]=fx[0]-area*(tau*sx[0]+sigma*nx);
			 fx[5]=fx[5]+area*(tau*sx[0]+sigma*nx);
			 fy[0]=fy[0]-area*(tau*sy[0]+sigma*ny);
			 fy[5]=fy[5]+area*(tau*sy[0]+sigma*ny);
			 fz[0]=fz[0]-area*(tau*sz[0]+sigma*nz);
			 fz[5]=fz[5]+area*(tau*sz[0]+sigma*nz);
			 fx[1]=fx[1]-area*(tau*sx[0]+sigma*nx);
			 fx[4]=fx[4]+area*(tau*sx[0]+sigma*nx);
			 fy[1]=fy[1]-area*(tau*sy[0]+sigma*ny);
			 fy[4]=fy[4]+area*(tau*sy[0]+sigma*ny);
			 fz[1]=fz[1]-area*(tau*sz[0]+sigma*nz);
			 fz[4]=fz[4]+area*(tau*sz[0]+sigma*nz);
		    }
		    else if(integ==1)
		    {
			 fx[1]=fx[1]-area*(tau*sx[1]+sigma*nx);
			 fx[4]=fx[4]+area*(tau*sx[1]+sigma*nx);
			 fy[1]=fy[1]-area*(tau*sy[1]+sigma*ny);
			 fy[4]=fy[4]+area*(tau*sy[1]+sigma*ny);
			 fz[1]=fz[1]-area*(tau*sz[1]+sigma*nz);
			 fz[4]=fz[4]+area*(tau*sz[1]+sigma*nz);
			 fx[2]=fx[2]-area*(tau*sx[1]+sigma*nx);
			 fx[3]=fx[3]+area*(tau*sx[1]+sigma*nx);
			 fy[2]=fy[2]-area*(tau*sy[1]+sigma*ny);
			 fy[3]=fy[3]+area*(tau*sy[1]+sigma*ny);
			 fz[2]=fz[2]-area*(tau*sz[1]+sigma*nz);
			 fz[3]=fz[3]+area*(tau*sz[1]+sigma*nz);
		    }
		    else
		    {
			 fx[2]=fx[2]-area*(tau*sx[2]+sigma*nx);
			 fx[3]=fx[3]+area*(tau*sx[2]+sigma*nx);
			 fy[2]=fy[2]-area*(tau*sy[2]+sigma*ny);
			 fy[3]=fy[3]+area*(tau*sy[2]+sigma*ny);
			 fz[2]=fz[2]-area*(tau*sz[2]+sigma*nz);
			 fz[3]=fz[3]+area*(tau*sz[2]+sigma*nz);
			 fx[0]=fx[0]-area*(tau*sx[2]+sigma*nx);
			 fx[5]=fx[5]+area*(tau*sx[2]+sigma*nx);
			 fy[0]=fy[0]-area*(tau*sy[2]+sigma*ny);
			 fy[5]=fy[5]+area*(tau*sy[2]+sigma*ny);
			 fz[0]=fz[0]-area*(tau*sz[2]+sigma*nz);
			 fz[5]=fz[5]+area*(tau*sz[2]+sigma*nz);
		    }
*/
//		    fprintf(fp,"Fx[%ld]=%f\tFy[%ld]=%f\tFz[%ld]=%f\nFx[%ld]=%f\tFy[%ld]=%f\tFz[%ld]=%f\nFx[%ld]=%f\tFy[%ld]=%f\tFz[%ld]=%f\nFx[%ld]=%f\tFy[%ld]=%f\tFz[%ld]=%f\nFx[%ld]=%f\tFy[%ld]=%f\tFz[%ld]=%f\nFx[%ld]=%f\tFy[%ld]=%f\tFz[%ld]=%f\n",i0,d1nftx[i0],i0,d1nfty[i0],i0,d1nftz[i0],i1,d1nftx[i1],i1,d1nfty[i1],i1,d1nftz[i1],i2,d1nftx[i2],i2,d1nfty[i2],i2,d1nftz[i2],i3,d1nftx[i3],i3,d1nfty[i3],i3,d1nftz[i3],i4,d1nftx[i4],i4,d1nfty[i4],i4,d1nftz[i4],i5,d1nftx[i5],i5,d1nfty[i5],i5,d1nftz[i5]);
	       }

/*	       fprintf(fp1,"ielem=%ld\nFx[%ld]=%e\tFy[%ld]=%e\tFz[%ld]=%e\nFx[%ld]=%e\tFy[%ld]=%e\tFz[%ld]=%e\nFx[%ld]=%e\tFy[%ld]=%e\tFz[%ld]=%e\nFx[%ld]=%e\tFy[%ld]=%e\tFz[%ld]=%e\nFx[%ld]=%e\tFy[%ld]=%e\tFz[%ld]=%e\nFx[%ld]=%e\tFy[%ld]=%e\tFz[%ld]=%e\n",
		       ielem,i0,d1nftx[i0],i0,d1nfty[i0],i0,d1nftz[i0],i1,d1nftx[i1],i1,d1nfty[i1],i1,d1nftz[i1],i2,d1nftx[i2],i2,d1nfty[i2],i2,d1nftz[i2],
		       i3,d1nftx[i3],i3,d1nfty[i3],i3,d1nftz[i3],i4,d1nftx[i4],i4,d1nfty[i4],i4,d1nftz[i4],i5,d1nftx[i5],i5,d1nfty[i5],i5,d1nftz[i5]);

	       if(ielem==2)
		 {
		   fprintf(fp1,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
		       d1nftx[i0],d1nfty[i0],d1nftz[i0],d1nftx[i1],d1nfty[i1],d1nftz[i1],d1nftx[i2],d1nfty[i2],d1nftz[i2],
		       d1nftx[i3],d1nfty[i3],d1nftz[i3],d1nftx[i4],d1nfty[i4],d1nftz[i4],d1nftx[i5],d1nfty[i5],d1nftz[i5]);
		 }

	       if(ielem==3)
		 {
		   fprintf(fp1,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
		       fx[0],fy[0],fz[0],fx[1],fy[1],fz[1],fx[2],fy[2],fz[2],
		       fx[3],fy[3],fz[3],fx[4],fy[4],fz[4],fx[5],fy[5],fz[5]);
		 }
*/
	  }
     }

//     fclose(fp);
//     fclose(fp1);
}

/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/


/*! \brief calculate nodal forces 
 *  \param[in]     yde element database
 *  \param[in,out] ydn nodal database
 *  \param[in]     ydp property database
 *  \param[in]     ydk constants database
 *  \par Details:
 *  Yfd() ... concepts described in chapter 4, I think
 *
    \verbatim
    for each tetraheron
      for eash Gauss point
        calculate Jacobian matrix
        calculate Cartesian gradients
        calculate tensor
        calculate velocity gradients
        calculate Cauchy stress
        calculate some unnamed value for this iteration
      calculate the nodal forces
    \endverbatim
 *
 *
 */
//-PY changed it for 3D_fracture_coupling_with_multiphase     
// void Yfd(YDE yde, YDN ydn, YDP ydp, YDK ydk, YDB ydb)
 void Yfd(YDE yde, YDN ydn, YDP ydp, YDK ydk, YDB ydb, YDI ydi, YDC ydc,YDX ydx)
{ 
     INT iprop,i;
    DBL k,norm,nx,ny;

     /*Z zero nodal forces and masses */
     TzDBL1(ydn->d1nmct, ydn->mnopo);
     TzDBL1(ydn->d1nvct, ydn->mnopo);

     TzDBL2(ydn->d2nfc, ydn->mnodim, ydn->mnopo);
     TzDBL2(ydn->d2nft, ydn->mnodim, ydn->mnopo);
     TzDBL3(yde->d3tcs, yde->melem, NDIME, NDIME);
/*
     ydn->i1nopr[5]=1;
     ydn->i1nopr[6]=1;
     ydn->i1nopr[7]=1;
     ydn->i1nopr[8]=1;

     ydn->i1nopr[9]=2;
     ydn->i1nopr[10]=2;
     ydn->i1nopr[11]=2;
     ydn->i1nopr[12]=2;
*/  
     for(iprop=0; iprop<ydp->nprop; iprop++)
     { 
	  if(ydp->i1ptyp[iprop] == YTE3TET10ELS) /* small strain elastic tetrahetra  */
	  { 
	       Yfd3TET10ELS(  
		    yde->nelem, iprop, 
		    ydp->d1peks[iprop], ydp->d1pela[iprop], ydp->d1pemu[iprop], ydp->d1pero[iprop],
		    ydn->d2ncc, ydn->d2nci, ydn->d2nft, ydn->d1nmct, ydn->d2nvc,
		    yde->i1elpr, ydn->i1nopr, yde->i2elto, 
		    ydk->d3dsh, yde->d3tcs, ydk->d1lpmm, yde->d1emct
		    );
	  }
	  else if((ydp->i1ptyp[iprop])==(YTE3TET4ELS))
	  { 
	       //printf("ydn->d2nft 1 %f \n\n", ydn->d2nft[0][yde->nelem]);

	       Yfd3TET4ELS(  /* small strain elastic tetrahetra  */
		yde->nelem,iprop,
		ydp->d1peks[iprop],ydp->d1pela[iprop],ydp->d1pemu[iprop],
		ydp->d1pero[iprop],ydn->d2ncc[0],ydn->d2ncc[1],
		ydn->d2ncc[2],ydn->d2nci[0],ydn->d2nci[1],
		ydn->d2nci[2],ydn->d2nft[0],ydn->d2nft[1],
		ydn->d2nft[2],ydn->d1nmct,ydn->d2nvc[0],
		ydn->d2nvc[1],ydn->d2nvc[2],ydn->i1nopr,yde->i1elpr,
		yde->i2elto,yde->d3tcs,
		ydx->d1nap,ydn->d1nvct,ydx->i2eltoxb, ydx->nelem_t,
        ydc->dctime,ydc->dcrmpt, yde->nelemi, ydn->nnopi,yde->d1emct,
        ydn->d1ntc, ydn->d1nti, ydp->d1ctex[iprop],ydc->icoutf,ydc->ncstep);

	       //printf("ydn->d2nft 2 %f \n\n", ydn->d2nft[0][yde->nelem]);


          }
          else if((ydp->i1ptyp[iprop])==(YTE3TET4JOINT))
          {
	       Yfd3TET4JOINT(	/* joint element for 4-node tetrahedra  */
		    yde->nelem,
		    iprop,
		    ydp->d1peft[iprop],
		    ydp->d1pegfn[iprop],ydp->d1pegfs[iprop],
            ydp->d1peks[iprop],ydp->d1pepe[iprop],
		    ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2ncc[2],
		    ydn->d2nft[0],ydn->d2nft[1],ydn->d2nft[2],
		    ydn->d2nvc[0],ydn->d2nvc[1],ydn->d2nvc[2],
		    yde->d2elst[ydp->i1psde[iprop]],
		    yde->i1elpr,yde->i2elto,yde->i1eljo,
		    yde->i2eljp,
                    &(ydi->iiecff),ydi->i1iecn,yde->i1elcf,
                    ydi->i1iect,&(ydi->nicoup),&(ydi->diedi),
                    ydn->nneigh,ydn->i2nnei,yde->i2elfr,yde->i1elbe,ydn->i1nei,
		    ydp->d1pefr[iprop],yde->i1elcft,ydc->ncstep,  yde->nelemi, ydn->nnopi,
		    ydp->d1pcoh[iprop], ydp->d1picf[iprop],yde->i1elty, &(ydc->icoutf));


	  }	 
     }
    
    /* for SURE project only */
/*
    k=1.0e6;
    for(i=0;i<ydn->nnopo;i++)
    {
        if(ydn->d1nmct[i]>EPSILON)
        {
        nx=ydn->d2ncc[0][i];
        ny=ydn->d2ncc[2][i];
        V2DNor(norm,nx,ny);
        if(ydn->d2ncc[1][i]>=0.024)ydn->d2nft[1][i]=ydn->d2nft[1][i]+(0.024-ydn->d2ncc[1][i])*k-0.5*SQRT(k*ydn->d1nmct[i])*ydn->d2nvc[1][i];
        if((norm>=0.023)&&ydn->d2ncc[1][i]>0.0205)
        {
          ydn->d2nft[0][i]=ydn->d2nft[0][i]+(0.023-norm)*k*nx-0.5*SQRT(k*ydn->d1nmct[i])*ydn->d2nvc[0][i];
            ydn->d2nft[2][i]=ydn->d2nft[2][i]+(0.023-norm)*k*ny-0.5*SQRT(k*ydn->d1nmct[i])*ydn->d2nvc[2][i];
        }
    }
    }
  */  
}


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/



