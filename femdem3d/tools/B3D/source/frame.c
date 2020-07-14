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
* Dr Jiansheng Xiang Prof Antonio Munjiza or Dr John-Paul Latham 
* j.xiang@imperial.ac.uk a.munjiza@qmul.ac.uk 
* or j.p.latham@imperial.ac.uk 
* ******************************************************************* */ 

/*1 File frame.c  */

#include "frame.h"
/* Arrays input output and allocation routines */ 
 
/***************PUBLIC***************************************/
void CHRcpynoext(c1,c2)  /* copy no extension */
  CHR *c1; CHR *c2;
{ INT i;
  i=0;
  while((c2[i]!='\0')&&(c2[i]!='.')&&(i<300))
  { c1[i]=c2[i];
    i=i+1;
  }
  c1[i]='\0';
}

DBL *TalDBL1(m1)
  INT m1; 
{ INT isize=m1*sizeof(DBL);
  if(isize==0)return DBL1NULL;
  return (DBL*)MALLOC(isize);
}

DBL **TalDBL2(m2,m1)
  INT m2; INT m1; 
{ INT isize,i2;
  DBL     *p1;
  DBL    **p2;
  void    *v1;

  isize=sizeof(DBL*)*(m2+3)+
        sizeof(DBL )*(m2*m1+3);
  if(isize==0)return DBL2NULL;
  v1=MALLOC(isize);
  p2=(DBL**)v1;
  p1=(DBL*)v1;
  p1=p1+((m2+1)*sizeof(DBL**))/sizeof(DBL)+2;    
  for(i2=0;i2<m2;i2++)
  { p2[i2]=p1+i2*m1; 
  }
  return p2; 
}
DBL ***TalDBL3(m3,m2,m1)
  INT m3; INT m2; INT m1; 
{ INT isize,i2,i3;
  DBL     *p1b;
  DBL    **p2b,  **p2e, **p2;
  DBL   ***p3b, ***p3e;
  void    *v1;

  isize=sizeof(DBL**)*(m3+3)+
        sizeof(DBL* )*(m3*m2+3)+
        sizeof(DBL )*(m3*m2*m1+3);
  if(isize==0)return DBL3NULL;
  v1=MALLOC(isize);
  p3b=(DBL***)v1; p3e=p3b+m3+1;
  p2b=(DBL**)p3e; p2e=p2b+m3*m2+1;
  p1b=(DBL*)p2e;
  p2=p2b;    
  for(i3=0;i3<m3;i3++)
  { p2=p2b+i3*m2;
    p3b[i3]=p2;
    for(i2=0;i2<m2;i2++)
    { p2[i2]=p1b+i3*m2*m1+i2*m1;
  } }
  return p3b; 
}
INT *TalINT1(m1)
  INT m1; 
{ INT isize=m1*sizeof(INT);
  if(isize==0)return INT1NULL;
  return (INT*)MALLOC(isize);   
} 

INT **TalINT2(m2,m1)
  INT m2; INT m1; 
{ INT isize,i2;
  INT     *p1;
  INT    **p2;
  void    *v1;

  isize=sizeof(INT*)*(m2+3)+
        sizeof(INT )*(m2*m1+3);
  if(isize==0)return INT2NULL;
  v1=MALLOC(isize);
  p2=(INT**)v1;
  p1=(INT*)v1;
  p1=p1+((m2+1)*sizeof(INT**))/sizeof(INT)+2;
  for(i2=0;i2<m2;i2++)
  { p2[i2]=p1+i2*m1; 
  }
  return p2; 
}
INT ***TalINT3(m3,m2,m1)
  INT m3; INT m2; INT m1; 
{ INT isize,i2,i3;
  INT     *p1b;
  INT    **p2b,  **p2e, **p2;
  INT   ***p3b, ***p3e;
  void    *v1;

  isize=sizeof(INT**)*(m3+3)+
        sizeof(INT* )*(m3*m2+3)+
        sizeof(INT )*(m3*m2*m1+3);
  if(isize==0)return INT3NULL;
  v1=MALLOC(isize);
  p3b=(INT***)v1; p3e=p3b+m3+1;
  p2b=(INT**)p3e; p2e=p2b+m3*m2+1;
  p1b=(INT*)p2e;
  p2=p2b;    
  for(i3=0;i3<m3;i3++)
  { p2=p2b+i3*m2;
    p3b[i3]=p2;
    for(i2=0;i2<m2;i2++)
    { p2[i2]=p1b+i3*m2*m1+i2*m1;
  } }
  return p3b; 
}
/*--------------------------------------------*/
INT  Getname(argc, argv, name)
  INT argc; char **argv; CHR *name;
{ INT i;
  for(i=0;i<argc;i++)
  { if(CHRcmp(argv[i],name,2)==0)return (i+1);
  };
  return i;
}

void TreadDBL1(fptr,n1,d1aray) 
  FILE  *fptr;  INT n1; DBL *d1aray;
{ INT i;
  DBL dnum;
  for(i=0;i<n1;i++)
  { DBLr(fptr,&dnum);
    d1aray[i]=dnum;
} } 
 
void TreadDBL2(fptr,n1,n2,d2aray) 
  FILE  *fptr;  INT n1; INT n2; DBL **d2aray;
{ INT i1,i2;
  DBL dnum;
  for(i1=0;i1<n1;i1++)
  { for(i2=0;i2<n2;i2++)
    { DBLr(fptr,&dnum);
      d2aray[i2][i1]=dnum;
} } }
 
void TreadINT1(fptr,n1,i1aray) 
  FILE  *fptr;  INT n1; INT *i1aray;
{ INT i; INT inum;
  for(i=0;i<n1;i++)
  { INTr(fptr,&inum);
    i1aray[i]=inum;
} } 
 
void TreadINT2(fptr,n1,n2,i2aray) 
  FILE  *fptr;  INT n1; INT n2; INT **i2aray;
{ INT i1,i2;
  INT inum;
  for(i1=0;i1<n1;i1++)
  { for(i2=0;i2<n2;i2++)
    { INTr(fptr,&inum);
      i2aray[i2][i1]=inum;
} } }

void TformDBL1(fptr,dinit,m1,d1aray) 
  FILE  *fptr;  INT m1;  DBL dinit; DBL **d1aray; 
{ INT i1,n1;
  DBL dnum;
  DBL *d1=*d1aray;
  if(d1==DBL1NULL)d1=TalDBL1(m1);
  *d1aray=d1;
  for(i1=0;i1<m1;i1++)
  { d1[i1]=dinit;
  }
  if(fptr!=FILENULL)
  { INTr(fptr,&n1);
    for(i1=0;i1<n1;i1++)
    { DBLr(fptr,&dnum);
      d1[i1]=dnum;
} } }
void TformDBL2(fptr,dinit,m1,m2,d2aray) 
  FILE  *fptr;  INT m1; INT m2; DBL dinit; DBL ***d2aray; 
{ INT i1,i2,n1,n2,how;
  DBL dnum;
  DBL **d2=*d2aray;
  if(d2==DBL2NULL)d2=TalDBL2(m1,m2);
  *d2aray=d2;
  for(i1=0;i1<m1;i1++)
  { for(i2=0;i2<m2;i2++)
    { d2[i1][i2]=dinit;
  } }
  if(fptr!=FILENULL)
  { INTr(fptr,&how);
    INTr(fptr,&n1);
    INTr(fptr,&n2);
    if(how==12)
    { for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { DBLr(fptr,&dnum);
          d2[i1][i2]=dnum;
    } } }
    else
    { for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { DBLr(fptr,&dnum);
          d2[i1][i2]=dnum;
} } } } }
void TformDBL3(fptr,dinit,m1,m2,m3,d3aray) 
  FILE  *fptr;  INT m1; INT m2; INT m3; DBL dinit; DBL ****d3aray; 
{ INT i1,i2,i3,n1,n2,n3,how;
  DBL dnum;
  DBL ***d3=*d3aray;
  if(d3==DBL3NULL)d3=TalDBL3(m1,m2,m3);
  *d3aray=d3;
  for(i1=0;i1<m1;i1++)
  { for(i2=0;i2<m2;i2++)
    { for(i3=0;i3<m3;i3++)
      {  d3[i1][i2][i3]=dinit;     
  } } }
  if(fptr!=FILENULL)
  { INTr(fptr,&how);
    INTr(fptr,&n1);
    INTr(fptr,&n2);
    INTr(fptr,&n3);
    if(how==123)
    { for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { for(i3=0;i3<n3;i3++)
          { DBLr(fptr,&dnum);
            d3[i1][i2][i3]=dnum;
    } } } }
    else if(how==132)
    { for(i1=0;i1<n1;i1++)
      { for(i3=0;i3<n3;i3++)
        { for(i2=0;i2<n2;i2++)
          { DBLr(fptr,&dnum);
            d3[i1][i2][i3]=dnum;
    } } } }
    else if(how==213)
    { for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { for(i3=0;i3<n3;i3++)
          { DBLr(fptr,&dnum);
            d3[i1][i2][i3]=dnum;
    } } } }
    else if(how==231)
    { for(i2=0;i2<n2;i2++)
      { for(i3=0;i3<n3;i3++)
        { for(i1=0;i1<n1;i1++)
          { DBLr(fptr,&dnum);
            d3[i1][i2][i3]=dnum;
    } } } }
    else if(how==231) 
    { for(i2=0;i2<n2;i2++)
      { for(i3=0;i3<n3;i3++)
        { for(i1=0;i1<n1;i1++)
          { DBLr(fptr,&dnum);
            d3[i1][i2][i3]=dnum;
    } } } }
    else if(how==312) 
    { for(i3=0;i3<n3;i3++)
      { for(i1=0;i1<n1;i1++)
        { for(i2=0;i2<n2;i2++)
          { DBLr(fptr,&dnum);
            d3[i1][i2][i3]=dnum;
    } } } }
    else if(how==321) 
    { for(i3=0;i3<n3;i3++)
      { for(i2=0;i2<n2;i2++)
        { for(i1=0;i1<n1;i1++)
          { DBLr(fptr,&dnum);
            d3[i1][i2][i3]=dnum;
} } } } } }

void TformINT1(fptr,iinit,m1,i1aray) 
  FILE  *fptr;  INT m1;  INT iinit; INT **i1aray; 
{ INT i1,n1;
  INT inum;
  INT *i1a=*i1aray;
  if(i1a==INT1NULL)i1a=TalINT1(m1);
  *i1aray=i1a;
  for(i1=0;i1<m1;i1++)
  { i1a[i1]=iinit;
  }
  if(fptr!=FILENULL)
  { INTr(fptr,&n1);
    for(i1=0;i1<n1;i1++)
    { INTr(fptr,&inum);
      i1a[i1]=inum;
} } }
void TformINT2(fptr,iinit,m1,m2,i2aray) 
  FILE  *fptr;  INT m1; INT m2; INT iinit; INT ***i2aray; 
{ INT i1,i2,n1,n2,how;
  INT inum;
  INT **i2a=*i2aray;
  if(i2a==INT2NULL)i2a=TalINT2(m1,m2);
  *i2aray=i2a;
  for(i1=0;i1<m1;i1++)
  { for(i2=0;i2<m2;i2++)
    { i2a[i1][i2]=iinit;
  } }
  if(fptr!=FILENULL)
  { INTr(fptr,&how);
    INTr(fptr,&n1);
    INTr(fptr,&n2);
    if(how==12)
    { for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { INTr(fptr,&inum);
          i2a[i1][i2]=inum;
    } } }
    else
    { for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { INTr(fptr,&inum);
          i2a[i1][i2]=inum;
} } } } }
void TformINT3(fptr,iinit,m1,m2,m3,i3aray) 
  FILE  *fptr;  INT m1; INT m2; INT m3; INT iinit; INT ****i3aray; 
{ INT i1,i2,i3,n1,n2,n3,how;
  INT inum;
  INT ***i3a=*i3aray;
  if(i3a==INT3NULL)i3a=TalINT3(m1,m2,m3);
  *i3aray=i3a;
  for(i1=0;i1<m1;i1++)
  { for(i2=0;i2<m2;i2++)
    { for(i3=0;i3<m3;i3++)
      { i3a[i1][i2][i3]=iinit;
  } } }
  if(fptr!=FILENULL)
  { INTr(fptr,&how);
    INTr(fptr,&n1);
    INTr(fptr,&n2);
    INTr(fptr,&n3);
    if(how==123)
    { for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { for(i3=0;i3<n3;i3++)
          { INTr(fptr,&inum);
            i3a[i1][i2][i3]=inum;
    } } } }
    else if(how==132)
    { for(i1=0;i1<n1;i1++)
      { for(i3=0;i3<n3;i3++)
        { for(i2=0;i2<n2;i2++)
          { INTr(fptr,&inum);
            i3a[i1][i2][i3]=inum;
    } } } }
    else if(how==213)
    { for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { for(i3=0;i3<n3;i3++)
          { INTr(fptr,&inum);
            i3a[i1][i2][i3]=inum;
    } } } }
    else if(how==231)
    { for(i2=0;i2<n2;i2++)
      { for(i3=0;i3<n3;i3++)
        { for(i1=0;i1<n1;i1++)
          { INTr(fptr,&inum);
            i3a[i1][i2][i3]=inum;
    } } } }
    else if(how==312)
    { for(i3=0;i3<n3;i3++)
      { for(i1=0;i1<n1;i1++)
        { for(i1=0;i1<n1;i1++)
          { INTr(fptr,&inum);
            i3a[i1][i2][i3]=inum;
    } } } }
    { for(i3=0;i3<n3;i3++)
      { for(i2=0;i2<n2;i2++)
        { for(i3=0;i3<n3;i3++)
          { INTr(fptr,&inum);
            i3a[i1][i2][i3]=inum;
} } } } } }
/*--------------------------------------------*/
FILE *TwriteCHR2(namep,fpt,ndigit,nperli,n1,c2aray,name)
  CHR   *namep; FILE *fpt; INT ndigit; INT nperli; INT  n1;
  CHR **c2aray; CHR *name;
{ INT i1;
  INT counter;
  FILE *fptr=fpt;   

  if(fptr==FILENULL)fptr=fopen(namep,"w");
  if((c2aray!=CHR2NULL)&&(fptr!=FILENULL))
  { CHRw(fptr,name) ; CHRwcr(fptr);
    counter=0;   
    for(i1=0;i1<n1;i1++)
    { CHRwsp(fptr); CHRwsp(fptr);  CHRw(fptr,c2aray[i1]);  counter++;
      if(counter>=nperli)
      {  counter=0; CHRwcr(fptr);
    } }
    if(counter>0)CHRwcr(fptr);
  }
  return fptr;
}

FILE *TwriteDBL1(namep,fpt,ndigit,nperli,n1,d1aray,name)
  CHR *namep; FILE *fpt; INT ndigit; INT nperli; INT  n1;
  DBL  *d1aray; CHR *name;
{ INT i1;
  INT counter;
  DBL dnum;
  FILE *fptr=fpt;   

  if(fptr==FILENULL)fptr=fopen(namep,"w");
  if((d1aray!=DBL1NULL)&&(fptr!=FILENULL))
  { CHRw(fptr,name) ; CHRwsp(fptr); INTw(fptr,n1,8); CHRwcr(fptr);
    counter=0;   
    for(i1=0;i1<n1;i1++)
    { dnum=d1aray[i1];
      CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
      if(counter>=nperli)
      {  counter=0; CHRwcr(fptr);
    } }
    if(counter>0)CHRwcr(fptr);
  }
  return fptr;
}

FILE *TwriteDBL2(namep,fpt,ndigit,nperli,ihow,n1,n2,d2aray,name)
  CHR *namep; FILE *fpt; INT ndigit; INT nperli; INT ihow; 
  INT  n1; INT n2; 
  DBL **d2aray; CHR *name;
{ INT i1,i2;
  INT counter;
  DBL dnum;   
  FILE *fptr=fpt;  

  if(fptr==FILENULL)fptr=fopen(namep,"w");
  if((d2aray!=DBL2NULL)&&(fptr!=FILENULL))
  { CHRw(fptr,name); INTw(fptr,ihow,8); 
    CHRwsp(fptr); INTw(fptr,n1,8);
    CHRwsp(fptr); INTw(fptr,n2,8);
    CHRwcr(fptr);
    counter=0; 
    if(ihow==12)  
    { for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { dnum=d2aray[i1][i2];
          CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
          if(counter>=nperli)
          { counter=0; CHRwcr(fptr);
    } } } }
    else  
    { for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { dnum=d2aray[i1][i2];
          CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
          if(counter>=nperli)
          { counter=0; CHRwcr(fptr);
    } } } }
    if(counter>0)CHRwcr(fptr);
  }
  return fptr;
}
FILE *TwriteDBL3(namep,fpt,ndigit,nperli,ihow,n1,n2,n3,d3aray,name)
  CHR *namep; FILE *fpt; INT ndigit; INT nperli; INT ihow;
  INT  n1; INT n2; INT n3;
  DBL ***d3aray; CHR *name;
{ INT i1,i2,i3;
  INT counter;
  DBL dnum;   
  FILE *fptr=fpt;  

  if(fptr==FILENULL)fptr=fopen(namep,"w");
  if((d3aray!=DBL3NULL)&&(fptr!=FILENULL))
  { CHRw(fptr,name); INTw(fptr,ihow,8); 
    CHRwsp(fptr); INTw(fptr,n1,8);
    CHRwsp(fptr); INTw(fptr,n2,8);
    CHRwsp(fptr); INTw(fptr,n3,8);
    CHRwcr(fptr);
    counter=0; 
    if(ihow==123)  
    { for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { for(i3=0;i3<n3;i3++)
          { dnum=d3aray[i1][i2][i3];
            CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
            if(counter>=nperli)
            { counter=0; CHRwcr(fptr);
    } } } } }
    if(ihow==132)  
    { for(i1=0;i1<n1;i1++)
      { for(i3=0;i3<n3;i3++)
        { for(i2=0;i2<n2;i2++)
          { dnum=d3aray[i1][i2][i3];
            CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
            if(counter>=nperli)
            { counter=0; CHRwcr(fptr);
    } } } } }    
    if(ihow==213)  
    { for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { for(i3=0;i3<n3;i3++)
          { dnum=d3aray[i1][i2][i3];
            CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
            if(counter>=nperli)
            { counter=0; CHRwcr(fptr);
    } } } } }    
    if(ihow==231)  
    { for(i2=0;i2<n2;i2++)
      { for(i3=0;i3<n3;i3++)
        { for(i1=0;i1<n1;i1++)
          { dnum=d3aray[i1][i2][i3];
            CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
            if(counter>=nperli)
            { counter=0; CHRwcr(fptr);
    } } } } }    
    if(ihow==312)  
    { for(i3=0;i3<n3;i3++)
      { for(i1=0;i1<n1;i1++)
        { for(i2=0;i2<n2;i2++)
          { dnum=d3aray[i1][i2][i3];
            CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
            if(counter>=nperli)
            { counter=0; CHRwcr(fptr);
    } } } } }    
    if(ihow==321)  
    { for(i3=0;i3<n3;i3++)
      { for(i2=0;i2<n2;i2++)
        { for(i1=0;i1<n1;i1++)
          { dnum=d3aray[i1][i2][i3];
            CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
            if(counter>=nperli)
            { counter=0; CHRwcr(fptr);
    } } } } }    
    if(counter>0)CHRwcr(fptr);
  }
  return fptr;
}
FILE *TwriteINT1(namep,fpt,ndigit,nperli,n1,i1aray,name)
  CHR *namep; FILE *fpt; INT ndigit; INT nperli; INT n1;
  INT  *i1aray; CHR *name;
{ INT i1;
  INT counter,istring;
  INT inum;   
  FILE *fptr=fpt;  

  if(fptr==FILENULL)fptr=fopen(namep,"w");
  if((i1aray!=INT1NULL)&&(fptr!=FILENULL))
  { CHRw(fptr,name) ; CHRwcr(fptr);
    counter=0;
    for(i1=0;i1<n1;i1++)
    { counter=MAXIM(counter,ABS(i1aray[i1]));
    }
    istring=0;
    while(counter>0)
    { counter=counter/10;
      istring++;
    }
    istring=MAXIM(istring+1,ndigit);

    counter=0;   
    for(i1=0;i1<n1;i1++)
    { inum=i1aray[i1];
      CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
      if(counter>=nperli)
      {  counter=0; CHRwcr(fptr);
    } }
    if(counter>0)CHRwcr(fptr);
  }
  return fptr;
}

FILE *TwriteINT2(namep,fpt,ndigit,nperli,ihow,n1,n2,i2aray,name)
  CHR *namep; FILE *fpt; INT ndigit; INT nperli; INT ihow; 
  INT n1; INT  n2;
  INT **i2aray; CHR *name;
{ INT i1,i2;
  INT counter,istring;
  INT inum; 
  FILE *fptr=fpt;  

  if(fptr==FILENULL)fptr=fopen(namep,"w");
  if((i2aray!=INT2NULL)&&(fptr!=FILENULL))
  { CHRw(fptr,name); INTw(fptr,ihow,8); 
    CHRwsp(fptr); INTw(fptr,n1,8);
    CHRwsp(fptr); INTw(fptr,n2,8);
    CHRwcr(fptr);

    counter=0;
    for(i1=0;i1<n1;i1++)
    { for(i2=0;i2<n2;i2++)
      { counter=MAXIM(counter,ABS(i2aray[i1][i2]));
    } }
    istring=0;
    while(counter>0)
    { counter=counter/10;
      istring++;
    }
    istring=MAXIM(istring+1,ndigit); 

    counter=0; 
    if(ihow==12)
    { counter=0;   
      for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { inum=i2aray[i1][i2];
          CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
          if(counter>=nperli)
          {  counter=0; CHRwcr(fptr);
    } } } }
    else
    { counter=0;   
      for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { inum=i2aray[i1][i2];
          CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
          if(counter>=nperli)
          {  counter=0; CHRwcr(fptr);
    } } } }
    if(counter>0)CHRwcr(fptr);
  }
  return fptr;
} 

FILE *TwriteINT3(namep,fpt,ndigit,nperli,ihow,n1,n2,n3,i3aray,name)
  CHR *namep; FILE *fpt; INT ndigit; INT nperli; INT ihow; 
  INT n1; INT  n2; INT n3;
  INT ***i3aray; CHR *name;
{ INT i1,i2,i3;
  INT counter,istring;
  INT inum; 
  FILE *fptr=fpt;  

  if(fptr==FILENULL)fptr=fopen(namep,"w");
  if((i3aray!=INT3NULL)&&(fptr!=FILENULL))
  { CHRw(fptr,name); INTw(fptr,ihow,8); 
    CHRwsp(fptr); INTw(fptr,n1,8);
    CHRwsp(fptr); INTw(fptr,n2,8);
    CHRwsp(fptr); INTw(fptr,n3,8);
    CHRwcr(fptr);

    counter=0;
    for(i1=0;i1<n1;i1++)
    { for(i2=0;i2<n2;i2++)
      { for(i3=0;i3<n3;i3++)
        { counter=MAXIM(counter,ABS(i3aray[i1][i2][i3]));
    } } }
    istring=0;
    while(counter>0)
    { counter=counter/10;
      istring++;
    }
    istring=MAXIM(istring+1,ndigit); 

    counter=0; 
    if(ihow==123)
    { counter=0;   
      for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { for(i3=0;i3<n3;i3++)
          { inum=i3aray[i1][i2][i3];
            CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
            if(counter>=nperli)
            {  counter=0; CHRwcr(fptr);
    } } } } }
    else if(ihow==132)
    { counter=0;   
      for(i1=0;i1<n1;i1++)
      { for(i3=0;i3<n3;i3++)
        { for(i2=0;i2<n2;i2++)
          { inum=i3aray[i1][i2][i3];
            CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
            if(counter>=nperli)
            {  counter=0; CHRwcr(fptr);
    } } } } }
    else if(ihow==213)
    { counter=0;   
      for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { for(i3=0;i3<n3;i3++)
          { inum=i3aray[i1][i2][i3];
            CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
            if(counter>=nperli)
            {  counter=0; CHRwcr(fptr);
    } } } } }
    else if(ihow==231)
    { counter=0;   
      for(i2=0;i2<n2;i2++)
      { for(i3=0;i3<n3;i3++)
        { for(i1=0;i1<n1;i1++)
          { inum=i3aray[i1][i2][i3];
            CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
            if(counter>=nperli)
            {  counter=0; CHRwcr(fptr);
    } } } } }
    else if(ihow==312)
    { counter=0;   
      for(i3=0;i3<n3;i3++)
      { for(i1=0;i1<n1;i1++)
        { for(i2=0;i2<n2;i2++)
          { inum=i3aray[i1][i2][i3];
            CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
            if(counter>=nperli)
            {  counter=0; CHRwcr(fptr);
    } } } } }
    else
    { counter=0;   
      for(i3=0;i3<n3;i3++)
      { for(i2=0;i2<n2;i2++)
        { for(i1=0;i1<n1;i1++)
          { inum=i3aray[i1][i2][i3];
            CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
            if(counter>=nperli)
            {  counter=0; CHRwcr(fptr);
    } } } } }    
    if(counter>0)CHRwcr(fptr);
  }
  return fptr;
}
/*4 **********SORTING ARRAYS BY SIZE - smallest ... largest number   */
/*4 Sorts double array smallest...                                   */
/*4 largest rearanges i1x in the same way                            */
void TsortDBL1(n,d1x,i1x)   
  INT n;  DBL *d1x;  INT *i1x; 
{ INT iblock,i,j,k,itmp;
  DBL tmp,m;
  INT ibig[50];
  INT iend[50];
  ibig[0]=0;
  iend[0]=n-1;
  iblock=0;
  while(iblock>=0)
  { i=ibig[iblock];
    j=iend[iblock];   
    iblock=iblock-1; 
    if(j>i)
    { tmp=d1x[j];
      m=tmp;
      for(k=i;k<j;k++)
      { tmp=MINIM(tmp,d1x[k]);
        m=MAXIM(m,d1x[k]);
      }
      if(tmp!=m)
      { m=(tmp+m)/2;
        while(i<=j)
        { while((i<=j)&&(d1x[i]<=m)) {i=i+1;}
          while((j>=i)&&(d1x[j]>m))  {j=j-1;}
          if(j>i)
          { tmp=d1x[j];
            d1x[j]=d1x[i];
            d1x[i]=tmp;
            if(i1x!=INT1NULL)
            { itmp=i1x[j];
              i1x[j]=i1x[i];
              i1x[i]=itmp;
        } } }
        ibig[iblock+2]=ibig[iblock+1];
        iend[iblock+2]=j;       
        ibig[iblock+1]=i;
        iblock=iblock+2;
} } } } 

/*4 Sorts nsort arrays     - smallest...largest    */
/*4 and rearange nrear-nsort arrays  same way      */ 
void TsortINT(n,nsort,nrear,i2)      
  INT n; INT nsort; INT nrear;  INT **i2; 
{ INT iblock,i,j,k,isort,tmp,m;
  INT *i1,*i1r;
  INT ibig[50];
  INT iend[50];
  ibig[0]=0;
  iend[0]=n-1;
  iblock=0;
  while(iblock>=0)
  { i=ibig[iblock];
    j=iend[iblock];   
    iblock=iblock-1; 
    if(j>i)
    { tmp=0; m=0; isort=0;
      while((isort<nsort)&&(m==tmp))
      { i1=i2[isort];
        isort=isort+1;
        tmp=i1[j]; 
        m=tmp;
        for(k=i;k<j;k++)
        { tmp=MINIM(tmp,i1[k]);
          m=MAXIM(m,i1[k]);
      } }
      if(tmp!=m)
      { m=(tmp+m)/2;
        while(i<=j)
        { while((i<=j)&&(i1[i]<=m)) {i=i+1;}
          while((j>=i)&&(i1[j]>m))  {j=j-1;}
          if(j>i)
          { for(isort=0;isort<nrear;isort++)
            { i1r=i2[isort];
              tmp=i1r[j];
              i1r[j]=i1r[i];
              i1r[i]=tmp;
        } } }
        ibig[iblock+2]=ibig[iblock+1];
        iend[iblock+2]=j;       
        ibig[iblock+1]=i;
        iblock=iblock+2;
} } } }
 
/*4 **************SAVING SPACE BY CODED ARRAYS**************/
#define ICODEBASE 90
static CHR *c1diga="0123456789-=qwertyuiop[]^asdfg";
static CHR *c1digb="hjkl;zxcvbnm,./~!@#$%&*()_+QWE";
static CHR *c1digc="RTYUIOP{}|ASDFGHJKL:ZXCVBNM<>?";
static CHR  c1dig[ICODEBASE+1];
static INT  i1dig[1000];
static INT  i1max[10];
static INT iffirst=0;

static void initcode()
{ INT i,j;

  iffirst=1;
  i1max[0]=1;
  for(i=1;i<10;i++)
  { i1max[i]=i1max[i-1]*ICODEBASE;
  }
  for(j=0;j<1000;j++)
  { i1dig[j]=0;
  }    
  CHRcpy(c1dig,c1diga);
  CHRcat(c1dig,c1digb);
  CHRcat(c1dig,c1digc);
  for(i=0;i<ICODEBASE;i++)
  { j=(INT)c1dig[i];
    if(j<0)j=j+500;
    if(j>1000)
    { CHRw(stderr,"Unexpected digit value");
      exit(1);
    }
    i1dig[j]=i;
} }
  
void codeCHRtoINT(c1code,i1num)
   CHR *c1code; INT *i1num;
{ INT inum,idig,ndigit,nnum,icar;

  if(iffirst==0)initcode();
  icar=(INT)c1code[0];
  if(icar<0)icar=icar+500;
  ndigit=i1dig[icar];
  icar=(INT)c1code[1];
  if(icar<0)icar=icar+500;
  nnum=i1dig[icar];
  i1num[0]=ndigit;
  i1num[1]=nnum; 
  for(inum=2;inum<nnum;inum++)
  { i1num[inum]=0;
    for(idig=ndigit-1;idig>=0;idig--)
    { icar=(INT)c1code[(inum-2)*ndigit+2+idig];
      if(icar<0)icar=icar+500;
      i1num[inum]=i1num[inum]+i1dig[icar]*i1max[idig];
} } }
 
void codeINTtoCHR(c1code,i1num)
   CHR *c1code; INT *i1num;
{ INT inum,idig,ival,ndigit,nnum;

  if(iffirst==0)initcode();
  ndigit=i1num[0];
  nnum=i1num[1];
  c1code[0]=c1dig[ndigit];
  c1code[1]=c1dig[nnum];
  c1code[(nnum-2)*ndigit+2]=CHRTERMINATE; 
  for(inum=2;inum<nnum;inum++)
  { ival=i1num[inum];
    for(idig=ndigit-1;idig>=0;idig--)
    { c1code[(inum-2)*ndigit+2+idig]=c1dig[(ival/i1max[idig])];
      ival=ival%i1max[idig];
} } }

void codeINTtoDBL(d1num,i1num)
  DBL *d1num; INT *i1num;
{ INT inum,nnum,ndigit;
  DBL dmax,dval;

  if(iffirst==0)initcode();
  ndigit=i1num[0];
  nnum=i1num[1];
  dmax=(DBL)(i1max[ndigit]-1);
  for(inum=0;inum<nnum;inum++)
  { dval=(DBL)i1num[inum];
    d1num[inum]=R2*dval/dmax-R1;
} }
 
void codeDBLtoINT(d1num,i1num)
  DBL *d1num; INT *i1num;
{ INT inum,nnum,ndigit;
  DBL dmax,dval;

  if(iffirst==0)initcode();
  ndigit=i1num[0];
  nnum=i1num[1];
  dmax=(DBL)(i1max[ndigit]-1);
  for(inum=2;inum<nnum;inum++)
  { dval=RP5*(d1num[inum]+R1)*dmax;
    i1num[inum]=(INT)dval;
} }
