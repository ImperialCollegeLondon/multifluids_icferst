/*! \file Y3Dinertia.c
 *  \brief Y arrays input output and allocation
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


#include "frame.h"


/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/

void symetricToTridiagonal(double **V, double *d, double  *e) 
{

//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.
//  The routine was optimized by Xavier Garcia and implemented in C++
//  Some enhancements were performed by the later because the original 
//  algorithm had some defficiencies when two of the three eigenvalues 
//  where too similar. This is a known problem in numerical methods that 
//  can be checked in the book "Numerical Recipes in C++"
//
//  The routine was tested ok.
int i,j,k;
double f,g,scale,h,hh;
  for (j = 0; j < 3; j++)
  {d[j] = V[3-1][j];}

  // Householder reduction to tridiagonal form.
  for (i = 3-1; i > 0; i--) 
    {
    // Scale to avoid under/overflow.
    scale = 0.0;
    h = 0.0;
    for (k = 0; k < i; k++) 
	{
      scale = scale + fabs(d[k]);
    }
    if (scale == 0.0) {
      e[i] = d[i-1];
      for (j = 0; j < i; j++) {
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
        V[j][i] = 0.0;
      }
    } else {

      // Generate Householder vector.

      for (k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      f = d[i-1];
      g = sqrt(h);
      if (f > 0) {
        g = -g;
      }
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for (j = 0; j < i; j++) {
        e[j] = 0.0;
      }

      // Apply similarity transformation to remaining columns.

      for (j = 0; j < i; j++) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for (k = j+1; k <= i-1; k++) {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for (j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      hh = f / (h + h);
      for (j = 0; j < i; j++) {
        e[j] -= hh * d[j];
      }
      for (j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (k = j; k <= i-1; k++) {
          V[k][j] -= (f * e[k] + g * d[k]);
        }
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  // Accumulate transformations.

  for (i = 0; i < 3-1; i++) 
   {
    V[3-1][i] = V[i][i];
    V[i][i] = 1.0;
    h = d[i+1];
    if (h != 0.0) {
      for (k = 0; k <= i; k++) {
        d[k] = V[k][i+1]/ h;
      }
      for (j = 0; j <= i; j++) {
        g = 0.0;
        for (k = 0; k <= i; k++) {
          g += V[k][i+1] * V[k][j];
        }
        for (k = 0; k <= i; k++) {
          V[k][j] -= g * d[k];
        }
      }
    }
    for (k = 0; k <= i; k++) {
      V[k][i+1] = 0.0;
    }
  }
  for (j = 0; j < 3; j++) 
   {
    d[j] = V[3-1][j];
    V[3-1][j] = 0.0;
   }
  V[3-1][3-1] = 1.0;
  e[0] = 0.0;
} 

// Symmetric tridiagonal QL algorithm.

#define solvecubicMAX(a, b) ((a)>(b)?(a):(b))
void solveEigenTridiagonal(double **V, double  *d, double *e)
 {

//  This is derived from the Algol procedures tql2, by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.
//
//  The code in C++ was writen by Xavier Garcia and optimized 
//  for 3x3 matrices.

  int i,j,l,m,k;
  double g,p,dl1,el1,h,r,c,c2,c3,s,s2,f,tst1,eps;

  f = 0.0;
  tst1 = 0.0;
  eps = pow(2.0,-52.0);

  for (i = 1; i < 3; i++) {
    e[i-1] = e[i];
  }
  e[3-1] = 0.0;
  for (l = 0; l < 3; l++) {

    // Find small subdiagonal element

    tst1 = solvecubicMAX(tst1,fabs(d[l]) + fabs(e[l]));
    m = l;
    while (m < 3) {
      if (fabs(e[m]) <= eps*tst1) {
        break;
      }
      m++;
    }

    // If m == l, d[l] is an eigenvalue,
    // otherwise, iterate.

    if (m > l) {
      int iter = 0;
      do {
        iter = iter + 1;  // (Could check iteration count here.)

        // Compute implicit shift

        g = d[l];
        p = (d[l+1] - g) / (2.0 * e[l]);
        r = sqrt(p*p+1.0);//hypot2(p,1.0);
        if (p < 0) {
          r = -r;
        }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        dl1 = d[l+1];
        h = g - d[l];
        for (i = l+2; i <3; i++) {
          d[i] -= h;
        }
        f = f + h;

        // Implicit QL transformation.

        p = d[m];
        c = 1.0;
        c2 = c;
        c3 = c;
        el1 = e[l+1];
        s = 0.0;
        s2 = 0.0;
        for (i = m-1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          //double hyp = sqrt(p*p+e[i]*e[i]);
          r = sqrt(p*p+e[i]*e[i]);//hypot2(p,e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);

          // Accumulate transformation.

          for (k = 0; k <3; k++) {
            h = V[k][i+1];
            V[k][i+1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;

        // Check for convergence.

      } while (fabs(e[l]) > eps*tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }
  
  // Sort eigenvalues and corresponding vectors.

  for (i = 0; i < 3-1; i++) {
    k = i;
    p = d[i];
    for (j = i+1; j <3; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j <3; j++) {
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
      }
    }
  }
}


void solveSymetricEigenProblem(DBL **m, DBL  *v1, DBL *v2, 
							   DBL  *v3, DBL *a1,DBL *a2,DBL *a3) 
{
DBL *Diagonal,*nonDiagonal,s;
DBL **m2;
int i,j;
m2=TalDBL2(3,3);
Diagonal=TalDBL1(3);
nonDiagonal=TalDBL1(3);
for( i=0;i<3;i++){
Diagonal[i]=0;
nonDiagonal[i]=0;
 for(j=0;j<3;j++){
 m2[i][j]=m[i][j];
}}
 
  symetricToTridiagonal(m2, Diagonal, nonDiagonal);
  solveEigenTridiagonal(m2, Diagonal, nonDiagonal);
  *a1=Diagonal[0];*a2=Diagonal[1];*a3=Diagonal[2];
  //these are the inertia moments

  v1[0]=m2[0][0]; v1[1]=m2[1][0];v1[2]=m2[2][0];
  v2[0]=m2[0][1]; v2[1]=m2[1][1];v2[2]=m2[2][1];
  v3[0]=m2[0][2]; v3[1]=m2[1][2];v3[2]=m2[2][2];
  

 //filter the noise in the vectors
 //this is because typically there are number on the order of 10e-12 times the maximum 
 //which are really  noise
/*
 double max;
 
 double* v[3]; v[0]=v1;v[1]=v2;v[2]=v3; 

 
 for(int k=0;k<3;k++) 
  {
   max = fabs( (v[k])[0]); 
   for(int d=1;d<3;d++)
    {
     if (fabs( (v[k])[d])>max) 
     max = fabs( (v[k])[d]); 
    }

   for(int d=0;d<3;d++)
    {
     if (fabs( (v[k])[d])<1.0e-6*max) 
     (v[k])[d]=0.0; 
    }
 

  }
  */
  
 //check that they are unitary again
 //we really dont need this 
  V3DNor(s,v1[0],v1[1],v1[2]);
  V3DNor(s,v2[0],v2[1],v2[2]);
  V3DNor(s,v3[0],v3[1],v3[2]);
  FREE(m2);
  }






void subExpression(DBL w0, DBL w1, DBL w2,DBL *f1,
                   DBL *f2, DBL *f3, DBL *g0,DBL *g1, 
                   DBL *g2)
  {
	DBL temp0,temp1,temp2;
    temp0 = w0+w1; 
    *f1 = temp0+w2; 
    temp1=w0*w0; 
    temp2=temp1+w1*temp0;
    *f2 = temp2+w2*(*f1); *f3 = w0*temp1+w1*temp2+w2*(*f2);
    *g0=(*f2)+w0*((*f1)+w0);*g1= (*f2)+w1*((*f1)+w1); *g2 = (*f2)+w2*((*f1)+w2);
  }


void   mirtichRoutine(DBL **points,
                      DBL *intg, 
                      DBL dro)
  {
 // setInCenter(points,nTriangles,triangleIndexes,inertia,volume,cm);
  //const double mult[10]={1./6.,1./24.,1./24.,1./24.,1./60.,1./60.,1./60.,1./120.,1./120.,1./120.};
  DBL x1,y1,z1,x2,y2,z2,x0,y0,z0,d0,d1,d2;
//  int i0,i1,i2;
  DBL  a1,b1,c1,a2,b2,c2;
  DBL f1x,f1y,f1z,f2x,f2y,f2z,f3x,f3y,f3z,g0x,g0y,g0z,g1x,g1y,g1z,g2x,g2y,g2z;
  INT t;
  INT p1[4] = { 0, 1, 2, 3 };
  INT p2[4] = { 1, 3, 3, 1 };
  INT p3[4] = { 2, 2, 0, 0 };
  

    for(t=0;t<4;t++)
      {
      //DBL values
      x0= points[p1[t]][0];
      y0= points[p1[t]][1];
      z0= points[p1[t]][2];

      x1= points[p2[t]][0];
      y1= points[p2[t]][1];
      z1= points[p2[t]][2];
      
      x2= points[p3[t]][0];
      y2= points[p3[t]][1];
      z2= points[p3[t]][2];

      a1= x1-x0; 
      b1= y1-y0;  
      c1= z1-z0;
      
      a2= x2-x0; 
      b2= y2-y0;  
      c2= z2-z0;
      
      d0  = b1*c2-b2*c1;
      d1  = a2*c1-a1*c2;
      d2  = a1*b2-a2*b1;
 
      subExpression(x0,x1,x2,&f1x,&f2x,&f3x,&g0x,&g1x,&g2x);
      subExpression(y0,y1,y2,&f1y,&f2y,&f3y,&g0y,&g1y,&g2y);
      subExpression(z0,z1,z2,&f1z,&f2z,&f3z,&g0z,&g1z,&g2z);
     
      intg[0]+=d0*f1x*dro;
      intg[1]+=d0*f2x*dro;
      intg[2]+=d1*f2y*dro;
      intg[3]+=d2*f2z*dro;
      intg[4]+=d0*f3x*dro;
      intg[5]+=d1*f3y*dro;
      intg[6]+=d2*f3z*dro;
      intg[7]+=d0*(y0*g0x + y1*g1x + y2*g2x)*dro;
      intg[8]+=d1*(z0*g0y + z1*g1y + z2*g2y)*dro;
      intg[9]+=d2*(x0*g0z + x1*g1z + x2*g2z)*dro;
     }
/*
    for(int n=0;n<10;n++)
      intg[n]*=(REAL)(mult[n]);

     //this is an error in the original publication!!!!!!!!!!!!!!!!!
     //they dont have the three!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DBL mass = intg[0]; volume = (REAL)(mass/3.0); //density = 1.0 assumed 

     //center of mass 
    cm[0]=intg[1]/mass;
    cm[1]=intg[2]/mass;
    cm[2]=intg[3]/mass;

    
    inertia(0,0) = intg[5] + intg[6] -mass*(cm[1]*cm[1]+cm[2]*cm[2]);
    inertia(1,1) = intg[4] + intg[6] -mass*(cm[0]*cm[0]+cm[2]*cm[2]);
    inertia(2,2) = intg[4] + intg[5] -mass*(cm[1]*cm[1]+cm[0]*cm[0]);

    inertia(0,1) = -(intg[7] -mass*cm[1]*cm[0]);
    inertia(1,0)=inertia(0,1);

    inertia(1,2) = -(intg[8] -mass*cm[1]*cm[2]);
    inertia(2,1)=inertia(1,2);

    inertia(0,2) = -(intg[9] -mass*cm[2]*cm[0]);
    inertia(2,0)=inertia(0,2);

    inertia*=((REAL)(1.0/3.0));

    //this is because the vertices could have been DoubleLinkedListed in the wrong direction
    if((inertia(0,0)<0.0)||(inertia(1,1)<0.0)||(inertia(2,2)<0.0))
      {
      inertia*=(-1.0);
      volume*=(-1.0);
      mass*=(-1.0);
      }*/
  }
/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/


