/*! \file Y3Dinit.c
 *  \brief Y constants initialization
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


/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/


/* calculate consistent mass matrix
*/
static void CalConMassM(DBL **d2csmm)
{
  DBL d1shapf[10][4];
  DBL d1xi[4]   = { alpha, beta,  beta,  beta };
  DBL d1eta[4]  = { beta,  alpha, beta,  beta };
  DBL d1zeta[4] = { beta,  beta,  alpha, beta };
  INT ig,ng,i,j;

  ng=4; /*Z NOTE:  ng != NGRSH */

  for(i=0;i<NNODE;i++)
  {
    for(j=0;j<NNODE;j++)
    {
      d2csmm[i][j]=R0;
    }
  }
 
  for(ig=0;ig<ng;ig++)
  {
    d1shapf[0][ig]=(R1-R2*d1xi[ig]-R2*d1eta[ig]-R2*d1zeta[ig])*(R1-d1xi[ig]-d1eta[ig]-d1zeta[ig]);
    d1shapf[1][ig]=d1xi[ig]*(R2*d1xi[ig]-R1);
    d1shapf[2][ig]=d1eta[ig]*(R2*d1eta[ig]-R1);
    d1shapf[3][ig]=d1zeta[ig]*(R2*d1zeta[ig]-R1);
    d1shapf[4][ig]=R4*d1xi[ig]*(R1-d1xi[ig]-d1eta[ig]-d1zeta[ig]);
    d1shapf[5][ig]=R4*d1xi[ig]*d1eta[ig];
    d1shapf[6][ig]=R4*d1eta[ig]*(R1-d1xi[ig]-d1eta[ig]-d1zeta[ig]);
    d1shapf[7][ig]=R4*d1zeta[ig]*(R1-d1xi[ig]-d1eta[ig]-d1zeta[ig]);
    d1shapf[8][ig]=R4*d1xi[ig]*d1zeta[ig];
    d1shapf[9][ig]=R4*d1zeta[ig]*d1eta[ig];

    for(i=0;i<NNODE;i++)
    {
      for(j=0;j<NNODE;j++)
      {
        d2csmm[i][j] += WEIGHT*d1shapf[i][ig]*d1shapf[j][ig];
      }
    }
  }
}


/**********************************************************************/
/**********************************************************************/


/* calculate lumped mass matrix
*/
static void CalLumMassM(DBL **d2csmm, DBL *d1lpmm)
{
  DBL sum;
  INT i;

  for(sum=R0,i=0;i<NNODE;i++) sum += d2csmm[i][i];   /*Z sum of diagonal */

  for(i=0;i<NNODE;i++) d1lpmm[i] = d2csmm[i][i]/sum; /*Z normalized diagonal */
}


/**********************************************************************/
/**********************************************************************/


/* quadratic shape function (deformation gradient?)
*/
static void SHAPEGRAD(DBL ***d3dsh)
{
  DBL d1xi[NGRSH]   = { RP25, alpha, beta,  beta,  beta };
  DBL d1eta[NGRSH]  = { RP25, beta,  alpha, beta,  beta };
  DBL d1zeta[NGRSH] = { RP25, beta,  beta,  alpha, beta };
  INT ig;

  for(ig=0;ig<NGRSH;ig++)
  {    	
    d3dsh[ig][0][0]=R4*(d1xi[ig]+d1eta[ig]+d1zeta[ig])-R3;
    d3dsh[ig][0][1]=R4*d1xi[ig]-R1;
    d3dsh[ig][0][2]=R0;
    d3dsh[ig][0][3]=R0;
    d3dsh[ig][0][4]=R4*(R1-R2*d1xi[ig]-d1eta[ig]-d1zeta[ig]);
    d3dsh[ig][0][5]=R4*d1eta[ig];
    d3dsh[ig][0][6]=-R4*d1eta[ig];
    d3dsh[ig][0][7]=-R4*d1zeta[ig];
    d3dsh[ig][0][8]=R4*d1zeta[ig];
    d3dsh[ig][0][9]=R0;

    d3dsh[ig][1][0]=R4*(d1xi[ig]+d1eta[ig]+d1zeta[ig])-R3;
    d3dsh[ig][1][1]=R0;
    d3dsh[ig][1][2]=R4*d1eta[ig]-R1;
    d3dsh[ig][1][3]=R0;
    d3dsh[ig][1][4]=-R4*d1xi[ig];
    d3dsh[ig][1][5]=R4*d1xi[ig];
    d3dsh[ig][1][6]=R4*(R1-d1xi[ig]-R2*d1eta[ig]-d1zeta[ig]);
    d3dsh[ig][1][7]=-R4*d1zeta[ig];
    d3dsh[ig][1][8]=R0;
    d3dsh[ig][1][9]=R4*d1zeta[ig];

    d3dsh[ig][2][0]=R4*(d1xi[ig]+d1eta[ig]+d1zeta[ig])-R3;
    d3dsh[ig][2][1]=R0;
    d3dsh[ig][2][2]=R0;
    d3dsh[ig][2][3]=R4*d1zeta[ig]-R1;
    d3dsh[ig][2][4]=-R4*d1xi[ig];
    d3dsh[ig][2][5]=R0;
    d3dsh[ig][2][6]=-R4*d1eta[ig];
    d3dsh[ig][2][7]=R4*(R1-d1xi[ig]-d1eta[ig]-R2*d1zeta[ig]);
    d3dsh[ig][2][8]=R4*d1xi[ig];
    d3dsh[ig][2][9]=R4*d1eta[ig];
  }
}


/**********************************************************************/
/**********************************************************************/


static int testIEEE(void)
{
  double zero_flt = 0.0;
  double zero_bits = 42.0;
  memset(&zero_bits, 0, sizeof(double));
  return (0 == memcmp(&zero_flt, &zero_bits, sizeof(double))); /* (zero_flt == zero_bits) */
}


/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/


/*! \brief constants database initialization  
 *  \par Details:
 *  Ykd() initializes all global constant elements of the YKD structure...
 *
 */
void Yinit(YDK ydk)
{
  ydk->d2csmm = TalDBL2(NNODE, NNODE);
  ydk->d1lpmm = TalDBL1(NNODE);
  ydk->d3dsh  = TalDBL3(NGRSH, NDIME, NNODE);

  CalConMassM(ydk->d2csmm);

  CalLumMassM(ydk->d2csmm, ydk->d1lpmm);

  SHAPEGRAD(ydk->d3dsh);

  ydk->d2ka = DBL2NULL;

  /*Z meh */
  initcode();

  /*Z if ever this is reported returning false, TzXXXn() will need revised */
  if(!testIEEE())
  {
    CHRw(stderr,"Yinit: unexpected IEEE floating point implementation...\n");
    CHRw(stderr,"       contact authors for a software patch");
    CHRwcr(stderr); exit(1);
  }
}


/*! \brief constants database configuration for current run  
 *  \par Details:
 *  Ykd() initializes all per-run constant elements of the YKD structure...
 *
 */
void Yconfig(YDK ydk, YDC ydc, YDB ydb)   
{
  INT i;

  if(ydk->d2ka) FREE(ydk->d2ka);

  ydk->d2ka = TalDBL2(NDIME, ydb->mbcon);

  for(i=0; i<ydb->nbcon; i++)
  {
    ydk->d2ka[0][i] = ydb->d1bnax[i] + ydc->dcgrax;
    ydk->d2ka[1][i] = ydb->d1bnay[i] + ydc->dcgray;
    ydk->d2ka[2][i] = ydb->d1bnaz[i] + ydc->dcgraz;
  }
}


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/



