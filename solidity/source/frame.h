/*! \file frame.h
 *  \brief Y macros and prototypes (Generic header for all modules)
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


#ifndef FRAME_H
#define FRAME_H


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>


/**********************************************************************/
/**********************************************************************/


/*Z code optimization tests */
//#define OPT1

/*Z timer support */
#include <time.h>
#define TIMER clock_t
#define TIME clock()
#define TIMERw(file,t,d,u) {fprintf((file), "elapsed time: %0.2f %s", (TIME - (t))/(double)((d)*CLOCKS_PER_SEC), u);}

/*Z this macro shortens code for many repeated conditionals... was put in place during testing, should keep it */
#define NXT(j,i,m) {(j)=((i)+1); if((j)>(m))(j)=0;}
/*Z notes about above macro.  There 3 obvious ways to do this:

  A)  j = i + 1;  if(j > 2) j = 0;

  B)  j = (i + 1) % 3;

  C)  int nxt[3] = { 1, 2, 0 };
      j = nxt[i];

I have used (B) for years since it is elegant, avoids a conditional, and integer math is fast...  
times have changed, so now (B) is the slowest method.

For unoptimized code, due to speculative execution, (A) is the fastest

For optimized code, due to loop unrolling, (A) and (C) are identical

I have therefore used (A), since we can gauge code improvement using either optimized or unoptimized releases, 
although I find (C) to be more elegant.  

I would expect (C) to be the fastest on older processors, if loop unrolling was not applicable 
(e.g. large modulus or large code block within loop).  This is not applicable here.

*/

/*Z these are assumed constant for Y3D? */
#define INITER 2
#define NELNO 10
#define NNODE 10
#define NNODEX 4
#define NDIME  3
#define NGRSH  5
#define WEIGHT RP25

/*Z these may need some work */
#if defined(_MSC_VER)
# define INLINE __forceinline
#elif defined(__GNUC__)
# define INLINE __attribute((always_inline))
#else
# define INLINE __inline
#endif

/*Z all reasonable platforms will have 0.0 = all bits off */
#define TzINT1(d1,m1)       memset((d1),       0, sizeof(INT)*(m1))
#define TzDPT1(d1,m1)       memset((d1),       0, sizeof(DPT)*(m1))
#define TzDMX1(d1,m1)       memset((d1),       0, sizeof(DMX)*(m1))
#define TzDBL1(d1,m1)       memset((d1),       0, sizeof(DBL)*(m1))
#define TzDBL2(d2,m2,m1)    memset((d2)[0],    0, sizeof(DBL)*((m2)*(m1)+3))     
#define TzDBL3(d3,m3,m2,m1) memset((d3)[0][0], 0, sizeof(DBL)*((m3)*(m2)*(m1)+3))

#define TalzDBL1(m1) (DBL*)MALLOCZ((m1),sizeof(DBL));
#define TalzINT1(m1) (INT*)MALLOCZ((m1),sizeof(INT));

/*Z fast alternative to TalINT2[][3] */
typedef int IPT[NDIME];
#define IPT1NULL ((IPT*)NULL)
#define TalIPT1(m1) (IPT*)MALLOC((m1)*sizeof(IPT));
#define TalzIPT1(m1) (IPT*)MALLOCZ((m1),sizeof(IPT));

/*Z fast alternative to TalDBL2[][3] */
typedef double DPT[NDIME];
#define DPT1NULL ((DPT*)NULL)
#define TalDPT1(m1) (DPT*)MALLOC((m1)*sizeof(DPT));
#define TalzDPT1(m1) (DPT*)MALLOCZ((m1),sizeof(DPT));

/*Z fast alternative to TalDBL3[][3][3] */
typedef double DMX[NDIME][NDIME];
#define DMX1NULL ((DMX*)NULL)
#define TalDMX1(m1) (DMX*)MALLOC((m1)*sizeof(DMX));
#define TalzDMX1(m1) (DMX*)MALLOCZ((m1),sizeof(DMX));

#define MALLOCZ calloc
#define SQR(a) ((a)*(a))

#define  Y_VERSION "1.0"

/*Z issue - long is inconsistent size across platforms... removed pending clarification
#define INT long
*/
#define FLT float
#define INS int
#define DBL double
#define INT long //long is needed
#define UCHR unsigned char
#define CHR char
#define UINT unsigned int
#define ULON unsigned long

#define SETLINEBUF(fcheck) setvbuf((fcheck), NULL, _IONBF, 0);
/* #define SETLINEBUF(fcheck) setlinebuf((fcheck));		SGI */

#define FILENULL ((FILE*)NULL)
#define CHR2NULL ((CHR **)NULL)
#define DBL1NULL ((DBL*)NULL)
#define DBL2NULL ((DBL**)NULL)
#define DBL3NULL ((DBL***)NULL)
#define INT1NULL ((INT*)NULL)
#define INT2NULL ((INT**)NULL)
#define INT3NULL ((INT***)NULL)

#define TRUE 1
#define FALSE 0
#define MYPI 3.1415926535897932384626
#define EPSILON 1.0e-10
#define BEPSILON 1.0e+15
#define NEPSILON -1.0e-10
#define R0 0.0
#define R1 1.0
#define R2 2.0
#define R3 3.0
#define R4 4.0
#define R5 5.0
#define R6 6.0
#define R7 7.0
#define R8 8.0
#define R9 9.0
#define R10 10.0
#define R12 12.0
#define R18 18.0
#define R19 19.0
#define R20 20.0  /*66 R20 was added by Dr J Xiang                  */
#define alpha 0.585410196624969  /*66 alpha was added by Dr J Xiang */
#define beta 0.138196601125011  /*66 beta was added by Dr J Xiang   */ 
#define R21 21.0
#define RP8 0.8
#define RP2 0.2
#define RP11 1.1
#define RP15 1.5
#define RP25 0.25
#define RP5  0.5
#define RP31  3.1
#define RP9  0.9
#define RP1  0.1
#define RP2  0.2
#define RP3  0.3
#define RP7  0.7
#define RP75  0.75
#define RP37  0.37


#define ANTEATOI atoi
#define ANTEATOF atof
#define EXP exp
#define ACOS acos
#define ASIN asin
#define ATAN atan
#define SQRT sqrt
#define COS cos
#define SIN sin

#ifdef UNIX
 #define FSQRT fsqrt
 #define FCOS fcos
 #define FSIN fsin
#endif

#ifdef LINUX
 #define FSQRT sqrt
 #define FCOS cos
 #define FSIN sin
#endif


#define ABS(x) (((x)<0)?-(x):(x))
#define DABS(x) (((x)<R0)?-(x):(x))
#define MAXIM(x,y) (((x)<(y))?(y):(x))
#define MINIM(x,y) (((x)>(y))?(y):(x))
#define FREE(x) if ((x)!=NULL) free(x)
#define MALLOCED(x) (x)
#define NOTMALLOCED(x) (!x)
#define MALLOC(x) malloc(x)


/**********************************************************************/
/* macros to perform operations on  small matrices                    */
/**********************************************************************/


/* INVERSE A SMALL MATRIX - also return determinant  
*/
#define YMATINV3(m,minv,det)                        \
{                                                   \
  det=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-    \
      m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+    \
      m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);    \
  minv[0][0]=(m[1][1]*m[2][2]-m[1][2]*m[2][1])/det; \
  minv[1][0]=(m[1][2]*m[2][0]-m[1][0]*m[2][2])/det; \
  minv[2][0]=(m[1][0]*m[2][1]-m[1][1]*m[2][0])/det; \
  minv[0][1]=(m[0][2]*m[2][1]-m[0][1]*m[2][2])/det; \
  minv[1][1]=(m[0][0]*m[2][2]-m[0][2]*m[2][0])/det; \
  minv[2][1]=(m[0][1]*m[2][0]-m[0][0]*m[2][1])/det; \
  minv[0][2]=(m[0][1]*m[1][2]-m[0][2]*m[1][1])/det; \
  minv[1][2]=(m[0][2]*m[1][0]-m[0][0]*m[1][2])/det; \
  minv[2][2]=(m[0][0]*m[1][1]-m[0][1]*m[1][0])/det; \
}

/*66 YMATDET3 was added by Dr J Xiang ---- 
*/
#define YMATDET3(m,det)                             \
{                                                   \
  det=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-    \
      m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+    \
      m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);    \
}

/*66 YMATINV2 was added by Dr J Xiang ---- 
*/
#define YMATINV2(m,minv,det)                        \
{                                                   \
  det=m[0][0]*m[1][1]-m[1][0]*m[0][1];              \
  minv[0][0]= m[1][1]/det; minv[1][0]=-m[1][0]/det; \
  minv[0][1]=-m[0][1]/det; minv[1][1]= m[0][0]/det; \
}


/**********************************************************************/
/*  macros to perform operations on 3D vectors (nice mess)            */
/**********************************************************************/

/* Assign; Null */
#define V3DAss(x1,y1,z1,x2,y2,z2) {(x1)=(x2);(y1)=(y2);z1=(z2);}
#define V3DNul(x1,y1,z1) {(x1)=R0;(y1)=R0;(z1)=R0;}
/* Vector Add; Subtract; Scale; Length square; Length; 
   Dot product; Cross product; Volume of triad */
#define V3DAdd(x1,y1,z1,x2,y2,z2,x3,y3,z3) \
              {(x1)=(x2)+(x3);(y1)=(y2)+(y3);(z1)=(z2)+(z3);}
#define V3DSub(x1,y1,z1,x2,y2,z2,x3,y3,z3) \
              {(x1)=(x2)-(x3);(y1)=(y2)-(y3);(z1)=(z2)-(z3);}
#define V3DSca(x1,y1,z1,s){(x1)=(x1)*(s);(y1)=(y1)*(s);(z1)=(z1)*(s);}
#define V3DDiv(x1,y1,z1,s){(x1)=(x1)/(s);(y1)=(y1)/(s);(z1)=(z1)/(s);}
#define V3DLe2(s,x1,y1,z1){(s)=(x1)*(x1)+(y1)*(y1)+(z1)*(z1);}
#define V3DLen(s,x1,y1,z1){(s)=SQRT((x1)*(x1)+(y1)*(y1)+(z1)*(z1));}
#define V3DDot(s,x1,y1,z1,x2,y2,z2) {(s)=((x1)*(x2))+((y1)*(y2))+((z1)*(z2));}
#define V3DCro(x1,y1,z1,x2,y2,z2,x3,y3,z3)\
              {(x1)=(y2)*(z3)-(z2)*(y3); \
               (y1)=(z2)*(x3)-(x2)*(z3); \
               (z1)=(x2)*(y3)-(y2)*(x3); }
#define V3DOmegaRot(x1,y1,z1,x2,y2,z2,x3,y3,z3)\
              {(x1)=(x3)+(y2)*(z3)-(z2)*(y3); \
               (y1)=(y3)+(z2)*(x3)-(x2)*(z3); \
               (z1)=(z3)+(x2)*(y3)-(y2)*(x3); }
#define V3DTranToGl(x1,y1,z1,a,b,c,x2,y2,z2,x3,y3,z3)\
              {(x1)=(a)*(x2)+(b)*(x3)+(c)*((y2)*(z3)-(z2)*(y3)); \
               (y1)=(a)*(y2)+(b)*(y3)+(c)*((z2)*(x3)-(x2)*(z3)); \
               (z1)=(a)*(z2)+(b)*(z3)+(c)*((x2)*(y3)-(y2)*(x3)); }
#define V3DTranToLoc(x1,y1,z1,a,b,c,x2,y2,z2,x3,y3,z3)\
              {(x1)=(a)*(x2)+(b)*(y2)+(c)*(z2);  \
               (y1)=(a)*(x3)+(b)*(y3)+(c)*(z3); \
(z1)=(a)*((y2)*(z3)-(z2)*(y3))+(b)*((x3)*(z2)-(x2)*(z3))+(c)*((x2)*(y3)-(y2)*(x3)); }
#define V3DVol(v,x2,y2,z2,x3,y3,z3,x1,y1,z1) \
              {(v)=((y2)*(z3)-(z2)*(y3))*(x1)+ \
                   ((z2)*(x3)-(x2)*(z3))*(y1)+ \
                   ((x2)*(y3)-(y2)*(x3))*(z1); }
/*  Normalize a Vector; Halfway Vectors; */ /*Z faster */
#define V3DNor(s,x1,y1,z1) \
              {(s)=SQRT((x1)*(x1)+(y1)*(y1)+(z1)*(z1)); \
               if((s)>EPSILON) {(x1)/=(s);(y1)/=(s);(z1)/=(s);} }
#define V3DMid(x1,y1,z1,x2,y2,z2,x3,y3,z3) \
              {(x1)=(RP5)*((x2)+(x3)); \
               (y1)=(RP5)*((y2)+(y3)); \
               (z1)=(RP5)*((z2)+(z3)); }
/* Vector in Between two given vectors  */
#define V3DBet(x1,y1,z1,x2,y2,z2,x3,y3,z3,s2,s3) \
              {(x1)=(x2)*(s2)+(x3)*(s3); \
               (y1)=(y2)*(s2)+(y3)*(s3); \
               (z1)=(z2)*(s2)+(z3)*(s3); }


/**********************************************************************/
/*  macros to perform operations on 2D vectors (nice mess)            */
/**********************************************************************/


/* Assign; Null */
#define V2DAss(x1,y1,x2,y2) {(x1)=(x2);(y1)=(y2)}
#define V2DNul(x1,y1) {(x1)=R0;(y1)=R0;} 
/* Vector Add; Subtract; Scale; Length square; */
/* Length; Dot product; Cross product;  */
#define V2DAdd(x1,y1,x2,y2,x3,y3) {(x1)=(x2)+(x3);(y1)=(y2)+(y3);}
#define V2DSub(x1,y1,x2,y2,x3,y3) {(x1)=(x2)-(x3);(y1)=(y2)-(y3);}
#define V2DSca(x1,y1,s) {(x1)=(x1)*s;(y1)=(y1)*s;}
#define V2DLe2(s,x1,y1,z1) {(s)=(x1)*(x1)+(y1)*(y1);}
#define V2DLen(s,x1,y1,z1) {(s)=SQRT((x1)*(x1)+(y1)*(y1);}
#define V2DDot(s,x1,y1,x2) {(s)=(x1)*(x2)+(y1)*(y2);}
#define V2DCro(s,x1,y1,x2,y2){(s)=(x1)*(y2)-(y1)*(x2);} 
/*  Normalize a Vector; */  /*Z fixed */
#define V2DNor(s,x1,y1) \
              {(s)=SQRT((x1)*(x1)+(y1)*(y1)); \
              if((s)>EPSILON) {(x1)/=(s);(y1)/=(s);} }
/*  Compute a Vector Halfway Between Two Given Vectors; */ 
#define V2DMid(x1,y1,x2,y2,x3,y3) \
              {(x1)=0.5*((x2)+(x3)); \
               (y1)=0.5*((y2)+(y3)); }
/* Vector in Between two given vectors  */
#define V2DBet(x1,y1,x2,y2,x3,y3,s2,s3) \
              {(x1)=(x2)*(s2)+(x3)*(s3); \
               (y1)=(y2)*(s2)+(y3)*(s3); }


/**********************************************************************/
/*  macros to perform operations on STREAMS                           */
/**********************************************************************/


#define FILEND(fptr) (feof(fptr))
static CHR *DBL_S[20]=
{ 
  "%le"     ,"%+1.0le" ,"%+2.0le" ,"%+3.0le" ,"%+4.0le" ,"%+5.0le" ,
  "%+6.0le" ,"%+7.0le" ,"%+8.1le" ,"%+9.2le" ,"%+10.3le",
  "%+11.4le","%+12.5le","%+13.6le","%+14.7le","%+15.8le",
  "%+16.9le","%+17.10le","%+18.11le","%+19.12le" 
};
static CHR *FLT_S[20]=
{ 
  "%le"     ,"%+1.0le" ,"%+2.0le" ,"%+3.0le" ,"%+4.0le" ,"%+5.0le" ,
  "%+6.0le" ,"%+7.0le" ,"%+8.1le" ,"%+9.2le" ,"%+10.3le",
  "%+11.4le","%+12.5le","%+13.6le","%+14.7le","%+15.8le",
  "%+16.9le","%+17.10le","%+18.11le","%+19.12le" 
}; /*Z why is is this the same as DBL_S ? */
#define DBL_SR "%le" 
static CHR *INT_S[20]=
{
  "%ld","%1ld","%2ld","%3ld","%4ld","%5ld","%6ld","%7ld","%8ld","%9ld",
  "%10ld","%11ld","%12ld","%13ld","%14ld",
  "%15ld","%16ld" ,"%17ld","%18ld","%19ld"
};
#define INT_SR   "%ld"
 
#define CHR_S "%s"
#define CHRRETURN "\n"
#define CHRSPACE " "
#define CHRTERMINATE 0
#define INTr(file,x){fscanf((file),INT_SR,((x)));}
#define DBLr(file,x){fscanf((file),DBL_SR,((x)));}
#define CHRr(file,x){fscanf((file),CHR_S,((x)));}
#define INTw(file,x,ndigit)fprintf((file),(INT_S[(ndigit)]),((x))); 
#define SINTw(file,x,ndigit)sprintf((file),(INT_S[(ndigit)]),((x))); 
#define DBLw(file,x,ndigit) fprintf((file),(DBL_S[((ndigit))]),(x));
#define FLTw(file,x,ndigit) fprintf((file),(FLT_S[((ndigit))]),(x));
 
#define SCHRw(file,x){sprintf((file),CHR_S,(x));}
#define CHRw(file,x){fprintf((file),CHR_S,(x));}
#define CHRend(c,i){c[i]='\0';}
#define CHRwcr(file){fprintf((file),CHR_S,(CHRRETURN));}
#define CHRwsp(file){fprintf((file),CHR_S,(CHRSPACE));}
#define CHRcmp(c1,c2,n)(strncmp((c1),(c2),(n)))
#define CHRcpy(c1,c2)(strcpy((c1),(c2)))
#define CHRcat(c1,c2)(strcat((c1),(c2)))


/**********************************************************************/
/*Z  function prototypes (note: should always use these)              */
/**********************************************************************/


void CHRcpynoext(  /* copy no extension */
  CHR *c1,              /* to copy to    */
  CHR *c2               /* to copy from  */
); 

DBL *TalDBL1(  /* alocate array */
  INT m1               /* size          */
); 

DBL **TalDBL2(  /* alocate array */
  INT m2,              /* m2 1D arrays size          */
  INT m1               /* of size   m1 DBL           */ 
); 
DBL ***TalDBL3(  /* alocate array */
  INT m3,              /* m3 2D arrays size m2*m1         */
  INT m2,              /* m3*m2 1D arrays                 */
  INT m1               /* of size   m1 DBL                */ 
);
INT *TalINT1(  /* alocate array */
  INT m1               /* size          */
); 

INT **TalINT2(  /* alocate array */
  INT m2,              /* m2 1D arrays size          */
  INT m1               /* of size   m1 INT           */
); 
INT ***TalINT3(  /* alocate array */
  INT m3,              /* m3 2D arrays size m2*m1         */
  INT m2,              /* m3*m2 1D arrays                 */
  INT m1               /* of size   m1 INT                */ 
);


/**********************************************************************/
/**********************************************************************/


/* return i for argv[i] that follows after "name" like -i file */
INT Getname(  
  INT argc,          /* number of argument in argv      */
  CHR **argv,        /* arguments                        */
  CHR *name          /* prefix to be looked for like -i  */
);
void TformDBL1(  /* read array */
  FILE *fptr,          /* file          */
  DBL dinit,           /* initialise to */
  INT n1,              /* size          */
  DBL **d1aray          /* array pointer */
);

void TformDBL2(  /* read array */
  FILE *fptr,          /* file          */
  DBL dinit,           /* initialise to */
  INT n1,              /* size          */
  INT n2,              /* size          */
  DBL ***d2aray         /* array pointer */
);
void TformDBL2_inv(  /* read array in transposed order */
  FILE *fptr,          /* file          */
  DBL dinit,           /* initialise to */
  INT n1,              /* size          */
  INT n2,              /* size          */
  DBL ***d2aray         /* array pointer */
);
void TformDBL3(  /* read array */
  FILE *fptr,          /* file          */
  DBL dinit,           /* initialise to */
  INT n1,              /* size          */
  INT n2,              /* size          */
  INT n3,              /* size          */
  DBL ****d3aray         /* array pointer */
);
void TformDBL3_inv(  /* read array in transposed order - sort-of */
  FILE *fptr,          /* file          */
  DBL dinit,           /* initialise to */
  INT n1,              /* size          */
  INT n2,              /* size          */
  INT n3,              /* size          */
  DBL ****d3aray         /* array pointer */
);
void TformINT1(  /* read array */
  FILE *fptr,          /* file          */
  INT iinit,           /* initialise to */
  INT n1,              /* size          */
  INT **i1aray          /* array pointer */
);

void TformINT2(  /* read array */
  FILE *fptr,          /* file          */
  INT iinit,           /* initialise to */
  INT n1,              /* size          */
  INT n2,              /* size          */
  INT ***i2aray        /* array pointer */
);
void TformINT2_inv(  /* read array in transposed order */
  FILE *fptr,          /* file          */
  INT iinit,           /* initialise to */
  INT n1,              /* size          */
  INT n2,              /* size          */
  INT ***i2aray        /* array pointer */
);
void TformINT3(  /* read array */
#if NeedFunctionPrototypes 
  FILE *fptr,          /* file          */
  INT iinit,           /* initialise to */
  INT n1,              /* size          */
  INT n2,              /* size          */
  INT n3,              /* size          */
  INT ****i3aray         /* array pointer */
#endif
);

void TreadDBL1(  /* read array */
  FILE *fptr,          /* file          */
  INT n1,              /* size          */
  DBL *d1aray          /* array pointer */
);

void TreadDBL2(  /* read array */
  FILE *fptr,          /* file          */
  INT n1,              /* size          */
  INT n2,              /* size          */
  DBL **d2aray         /* array pointer */
);

void TreadINT1(  /* read array */
  FILE *fptr,          /* file          */
  INT n1,              /* size          */
  INT *i1aray          /* array pointer */
);

void TreadINT2(  /* read array */
  FILE *fptr,          /* file          */
  INT n1,              /* size          */
  INT n2,              /* size          */
  INT **i2aray         /* array pointer */
);


/**********************************************************************/
/**********************************************************************/


FILE *TwriteCHR2(  /* write array */
  CHR *namep,          /* file name                 */
  FILE *fpt,           /* file                      */
  INT ndigit,          /*number of digits  ignored  */
  INT nperli,          /* numbers per line          */
  INT n1,              /* size                      */
  CHR **c2aray,        /* array pointer             */
  CHR *name            /* array name                */
);

FILE *TwriteDBL1(  /* write array */
  CHR *namep,          /* file name       */
  FILE *fpt,           /* file             */
  INT ndigit,          /*number of digits  */
  INT nperli,          /* numbers per line */
  INT n1,              /* size             */
  DBL *d1aray,         /* array pointer    */
  CHR *name            /* array name       */
);

FILE *TwriteDBL2(  /* write array */
  CHR  *namep,         /* file name     */
  FILE *fpt,           /* file             */
  INT ndigit,          /*number of digits  */
  INT nperli,          /* numbers per line */
  INT ihow,            /* which order      */
  INT n1,              /* size             */
  INT n2,              /* size             */
  DBL **d2aray,        /* array pointer    */
  CHR *name            /* array name       */
);
FILE *TwriteDBL3(  /* write array */
  CHR  *namep,         /* file name     */
  FILE *fpt,           /* file             */
  INT ndigit,          /*number of digits  */
  INT nperli,          /* numbers per line */
  INT ihow,            /* which order      */
  INT n1,              /* size             */
  INT n2,              /* size             */
  INT n3,              /* size             */
  DBL ***d3aray,        /* array pointer    */
  CHR *name            /* array name       */
);
 
FILE *TwriteINT1(  /* write array */
  CHR  *namep,         /* file name     */
  FILE *fpt,           /* file             */
  INT ndigit,          /*number of digits  */
  INT nperli,          /* numbers per line */
  INT n1,              /* size             */
  INT *i1aray,         /* array pointer    */
  CHR *name            /* array name       */
);

FILE *TwriteINT2(  /* write array */
  CHR  *namep,         /* file name     */
  FILE *fpt,           /* file             */
  INT ndigit,          /*number of digits  */
  INT nperli,          /* numbers per line */
  INT ihow,            /* numbers per line */
  INT n1,              /* size             */
  INT n2,              /* size             */
  INT **i2aray,        /* array pointer    */
  CHR *name            /* array name       */
);
FILE *TwriteINT3(  /* write array */
  CHR  *namep,         /* file name     */
  FILE *fpt,           /* file             */
  INT ndigit,          /*number of digits  */
  INT nperli,          /* numbers per line */
  INT ihow,            /* numbers per line */
  INT n1,              /* size             */
  INT n2,              /* size             */
  INT n3,              /* size             */
  INT ***i3aray,        /* array pointer    */
  CHR *name            /* array name       */
);


/**********************************************************************/
/* SORTING ARRAYS                                                     */
/**********************************************************************/


void TsortDBL1(  /* sort array = smallest ... largest*/
  INT n,               /* size                                  */
  DBL *d1x,            /* array pointer (array to be sorted)    */
  INT *i1y             /* array pointer (array to be rearanged) */
);
void TsortINT(   /* sort array = smallest ... largest*/
  INT n,               /* size                                  */
  INT nsort,           /* total number of arrays to be sorted   */
  INT nrear,           /* total number of arrays                */
  INT **i2             /* array pointers (sorted or rearanged)  */
); 


/**********************************************************************/
/* SPACESAVING by 90 FORMAT                                           */
/**********************************************************************/


void readcode(  /* read coded array */
  CHR *c1code,         /* coded array [0]=ndigit; [1]=nnum     */
  INT *i1num           /* INT array [0]=ndigit; [1]=nnum
                          the rest between 0 and base^ndigit    */
);
void initcode( /* initialize codes */
  void
);
void codeCHRtoINT(  /* write array as coded */
  CHR *c1code,         /* coded array [0]=ndigit; [1]=nnum     */
  INT *i1num           /* INT array [0]=ndigit; [1]=nnum
                          the rest between 0 and base^ndigit    */
);
void codeINTtoCHR(  /* write array as coded */
  CHR *c1code,         /* coded array [0]=ndigit; [1]=nnum     */
  INT *i1num           /* INT array [0]=ndigit; [1]=nnum
                          the rest between 0 and base^ndigit    */
);

void codeINTtoDBL(  /* translate INT to DBL array */
  DBL *d1num,          /* DBL array [0]=ndigit; [1]=nnum
                          all between -1 and +1                */  
  INT *i1num           /* INT array [0]=ndigit; [1]=nnum
                          the rest between 0 and base^ndigit    */
);

void codeDBLtoINT(  /* translate DBL to INT array */
  DBL *d1num,          /* DBL array [0]=ndigit; [1]=nnum
                          all between -1 and +1                */  
  INT *i1num           /* INT array [0]=ndigit; [1]=nnum
                          the rest between 0 and base^ndigit    */
);


/**********************************************************************/
/**********************************************************************/


#endif /* FRAME_H */


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/

