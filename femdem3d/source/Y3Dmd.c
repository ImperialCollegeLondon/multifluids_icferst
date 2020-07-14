/*! \file Y3Dmd.c by LG
 *  \mesh elements
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

/*JOINT ELEMENTS*/
static void Yjd3TET4(		/* mesh 4-node tetrahedra */
     INT melem, INT mnopo, INT nelest, INT nnopst,
     INT iprop, INT ipropj,
     INT *n0elem, INT *n0nopo,
     DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
     DBL *d1ncix, DBL *d1nciy, DBL *d1nciz,
     DBL *d1nvcx, DBL *d1nvcy, DBL *d1nvcz,
     DBL *d1sdel, 
     INT *i1elpr,
     INT *i1jnef, INT *i1jnen, INT *i1nobf, INT *i1nopr,
     INT **i2elto, INT **i2eljp, INT *i1elty,
     DBL *d1ntc, DBL *d1nti)
{
     INT nelem,nnopo;
     INT i,j,k,in,jn,kn,ln,ijnew,ielem;
     INT *i1elto,*i1eltonew;

//     FILE *fp;
//     fp=fopen("Yjd3TET4_output.txt","w+");

     nelem=(*n0elem);
     nnopo=(*n0nopo);
     for(ielem=0;ielem<nelest;ielem++)
     {
          if(i1elpr[ielem]==iprop)
	  {
	       i1elto=i2elto[ielem];
	       
//	       fprintf(fp,"OLD element:%ld\n",ielem);
//	       fprintf(fp,"OLD element nodes:%ld %ld %ld %ld\n",i1elto[0],i1elto[1],i1elto[2],i1elto[3]);

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
		    in=i1elto[i];
		    jn=i1elto[j];
		    kn=i1elto[k];
		    ln=MAXIM(in,jn);
                    ln=MAXIM(ln,kn);
		    if(nnopo>=mnopo)
		    {
			 CHRw(stderr,"Yjd: MNOPO too small");
			 CHRwcr(stderr);
			 exit(1);
		    }

		    /*create new node*/
		    d1nccx[nnopo]=d1nccx[in];
		    d1nccy[nnopo]=d1nccy[in];
		    d1nccz[nnopo]=d1nccz[in];
		    d1ncix[nnopo]=d1ncix[in];
		    d1nciy[nnopo]=d1nciy[in];
		    d1nciz[nnopo]=d1nciz[in];
		    d1nvcx[nnopo]=d1nvcx[in];
		    d1nvcy[nnopo]=d1nvcy[in];
		    d1nvcz[nnopo]=d1nvcz[in];
		    i1nopr[nnopo]=i1nopr[in];
		    i1nobf[nnopo]=1;
		    d1ntc[nnopo] = d1ntc[in];
		    d1nti[nnopo] = d1nti[in];


//		    fprintf(fp,"x[%ld]=%f\ty[%ld]=%f\tz[%ld]=%f\n",in,d1nccx[in],in,d1nccy[in],in,d1nccz[in]);
//		    fprintf(fp,"X[%ld]=%f\tY[%ld]=%f\tZ[%ld]=%f\n",nnopo,d1nccx[nnopo],nnopo,d1nccy[nnopo],nnopo,d1nccz[nnopo]);

		    nnopo=nnopo+1;

		    /*check if joint element already existent*/
		    ijnew=i1jnef[ln];

//                    fprintf(fp,"in=%ld,jn=%ld,kn=%ld,ln=%ld,ijnew=i1jnef[%ld]=%ld\n",in,jn,kn,ln,ln,ijnew);

		    while(ijnew>=0)
		    {
                         i1eltonew=i2elto[ijnew];
		      
			 if((i1eltonew[3]==in)&&(i1eltonew[4]==jn)&&(i1eltonew[5]==kn))
			 {
			      i1elpr[ijnew]=MAXIM(i1elpr[ijnew],ipropj);
			      if(i<2)
			      {
				   i1eltonew[3]=nnopo-1;
				   i1eltonew[4]=nnopo;
				   i1eltonew[5]=nnopo+1;
			      }
			      else if(i==2)
			      {
				   i1eltonew[3]=nnopo-1;
				   i1eltonew[4]=nnopo;
				   i1eltonew[5]=nnopo-3;
			      }
			      else
			      {
				   i1eltonew[3]=nnopo-1;
				   i1eltonew[4]=nnopo-4;
				   i1eltonew[5]=nnopo-3;
			      }
			      ijnew=-100;

//			      fprintf(fp,"***1***,ijnew=%ld\n",ijnew);
			 }
			 else if((i1eltonew[3]==jn)&&(i1eltonew[4]==kn)&&(i1eltonew[5]==in))
			 {
			      i1elpr[ijnew]=MAXIM(i1elpr[ijnew],ipropj);
			      if(i<2)
			      {
				   i1eltonew[3]=nnopo;
				   i1eltonew[4]=nnopo+1;
				   i1eltonew[5]=nnopo-1;
			      }
			      else if(i==2)
			      {
				   i1eltonew[3]=nnopo;
				   i1eltonew[4]=nnopo-3;
				   i1eltonew[5]=nnopo-1;
			      }
			      else
			      {
				   i1eltonew[3]=nnopo-4;
				   i1eltonew[4]=nnopo-3;
				   i1eltonew[5]=nnopo-1;
			      }
			      ijnew=-100;

//			      fprintf(fp,"***2***,ijnew=%ld\n",ijnew);
			 }
			 else if((i1eltonew[3]==kn)&&(i1eltonew[4]==in)&&(i1eltonew[5]==jn))
			 {
			      i1elpr[ijnew]=MAXIM(i1elpr[ijnew],ipropj);
			      if(i<2)
			      {
				   i1eltonew[3]=nnopo+1;
				   i1eltonew[4]=nnopo-1;
				   i1eltonew[5]=nnopo;
			      }
			      else if(i==2)
			      {
				   i1eltonew[3]=nnopo-3;
				   i1eltonew[4]=nnopo-1;
				   i1eltonew[5]=nnopo;
			      }
			      else
			      {
				   i1eltonew[3]=nnopo-3;
				   i1eltonew[4]=nnopo-1;
				   i1eltonew[5]=nnopo-4;
			      }
			      ijnew=-100;

//			      fprintf(fp,"***3***,ijnew=%ld\n",ijnew);
			 }
			 else if((i1eltonew[3]==kn)&&(i1eltonew[4]==jn)&&(i1eltonew[5]==in))
			 {
			      i1elpr[ijnew]=MAXIM(i1elpr[ijnew],ipropj);
			      if(i<2)
			      {
				   i1eltonew[3]=nnopo+1;
				   i1eltonew[4]=nnopo;
				   i1eltonew[5]=nnopo-1;
			      }
			      else if(i==2)
			      {
				   i1eltonew[3]=nnopo-3;
				   i1eltonew[4]=nnopo;
				   i1eltonew[5]=nnopo-1;				   
			      }
			      else
			      {
				   i1eltonew[3]=nnopo-3;
				   i1eltonew[4]=nnopo-4;
				   i1eltonew[5]=nnopo-1;
			      }
			      ijnew=-100;

//			      fprintf(fp,"***4***,ijnew=%ld\n",ijnew);
			 }
			 else if((i1eltonew[3]==jn)&&(i1eltonew[4]==in)&&(i1eltonew[5]==kn))
			 {
			      i1elpr[ijnew]=MAXIM(i1elpr[ijnew],ipropj);
			      if(i<2)
			      {
				   i1eltonew[3]=nnopo;
				   i1eltonew[4]=nnopo-1;
				   i1eltonew[5]=nnopo+1;
			      }
			      else if(i==2)
			      {
				   i1eltonew[3]=nnopo;
				   i1eltonew[4]=nnopo-1;
				   i1eltonew[5]=nnopo-3;
			      }
			      else
			      {
				   i1eltonew[3]=nnopo-4;
				   i1eltonew[4]=nnopo-1;
				   i1eltonew[5]=nnopo-3;
			      }
			      ijnew=-100;

//			      fprintf(fp,"***5***,ijnew=%ld\n",ijnew);
			 }
			 else if((i1eltonew[3]==in)&&(i1eltonew[4]==kn)&&(i1eltonew[5]==jn))
			 {
			      i1elpr[ijnew]=MAXIM(i1elpr[ijnew],ipropj);
			      if(i<2)
			      {
				   i1eltonew[3]=nnopo-1;
				   i1eltonew[4]=nnopo+1;
				   i1eltonew[5]=nnopo;
			      }
			      else if(i==2)
			      {
				   i1eltonew[3]=nnopo-1;
				   i1eltonew[4]=nnopo-3;
				   i1eltonew[5]=nnopo;
			      }
			      else
			      {
				   i1eltonew[3]=nnopo-1;
				   i1eltonew[4]=nnopo-3;
				   i1eltonew[5]=nnopo-4;
			      }
			      ijnew=-100;

//			      fprintf(fp,"***6***,ijnew=%ld\n",ijnew);
			 }
			 else
			 {
			      ijnew=i1jnen[ijnew-nelest];

//			      fprintf(fp,"ijnew=%ld\n",ijnew);
			 }
		    }

		    /*create new mid-edge joint node*/
		    if(ijnew>(-10))
		    {
			 if(nelem>=melem)
			 {
			      CHRw(stderr,"Yjd: MELEM too small");
			      CHRwcr(stderr);
			      exit(1);
			 }

			 i1eltonew=i2elto[nelem];
			 
//			 fprintf(fp,"NEW JOINT element:%ld\n",nelem);

			 if(i<2)
			 {
			      i1eltonew[2]=nnopo+1;
			      i1eltonew[1]=nnopo;
			      i1eltonew[0]=nnopo-1;
			 }
			 else if(i==2)
			 {
			      i1eltonew[2]=nnopo-3;
			      i1eltonew[1]=nnopo;
			      i1eltonew[0]=nnopo-1;
			 }
			 else
			 {
			      i1eltonew[2]=nnopo-3;
			      i1eltonew[1]=nnopo-4;
			      i1eltonew[0]=nnopo-1;
			 }

			 i1eltonew[3]=kn;
			 i1eltonew[4]=jn;
			 i1eltonew[5]=in;
//			 i2eljp[nelem][0]=ielem;

				 
			 i1elpr[nelem]=ipropj;
			 i1jnen[nelem-nelest]=i1jnef[ln];

			 i1jnef[ln]=nelem;
			 i1elty[nelem]=0;  //-ao asiri added 021017 (all joints are i1elty =0)

			 nelem=nelem+1;
		    }
	       }

	       /*detach element*/
	       for(i=0;i<4;i++)
	       {
		    i1elto[i]=nnopo-4+i;
	       }
	       
//	       fprintf(fp,"UPDATED OLD element nodes:%ld %ld %ld %ld\n",i1elto[0],i1elto[1],i1elto[2],i1elto[3]);
	  }
     }

     for(ielem=nelest;ielem<nelem;ielem++)
     {
	  if(d1sdel!=DBL1NULL)
	  {
	       d1sdel[ielem]=R0;
	  }
	  //for(i=3;i<6;i++)
	  //{
               i1elto=i2elto[ielem];
	       
	       if(i1elto[3]<nnopst)
	       {
		    i1elto[5]=i1elto[0];
		    i1elto[4]=i1elto[1];
		    i1elto[3]=i1elto[2];
		    
		    if (i1elty[ielem]<2)
		    {
			i1elty[ielem]=1; //-ao asiri added (boundary joints are i1elty = 1)
		    }
/*		    fprintf(fp,"BOUNDARY JOINT element:%ld\n",ielem);
		    fprintf(fp,"UPDATED nodes:%ld %ld %ld %ld %ld %ld\n",i1elto[0],i1elto[1],
			   i1elto[2],i1elto[3],i1elto[4],i1elto[5]);
*/	       }
	  //}

//	       fprintf(fp,"@@@**********@@@");
     }
    
     (*n0nopo)=nnopo;
     (*n0elem)=nelem;

//     fprintf(fp,"\n******Yjd3TET4 FINISHED******\n");

//     fclose(fp);
}

/*DISTINGUISH JOINT ELEMENTS*/
static void Yjd3JOINTPR(	/* distinguish two kinds of joint elements */
			INT nelest, INT nelem,
			INT **i2elto, INT *i1eljo, INT *i1elty)
{
  INT ielem,i,inopo,jnopo,knopo;
  INT joflag;

  for(ielem=nelest;ielem<nelem;ielem++) /* joint element loop */
    {
      joflag=(-1);

      for(i=0;i<nelest;i++)	/* finite element loop */
	{
	  for(inopo=0;inopo<4;inopo++)
	    {
	      jnopo=inopo+1;
	      if(jnopo>3)
		{
		  jnopo=0;
		}
	      knopo=jnopo+1;
	      if(knopo>3)
		{
		  knopo=0;
		}
	      
	      if((i2elto[ielem][0]==i2elto[i][inopo])&&(i2elto[ielem][1]==i2elto[i][jnopo])&&(i2elto[ielem][2]==i2elto[i][knopo]))
		{
		  joflag=inopo;
		  break;
		}
	    }
	  if(joflag!=(-1))
            {
              break;
            }
	}

      i1eljo[ielem]=joflag;
    }
}



/*LINK JOINT ELEMENT TO 4-NODE FINITE ELEMENT*/
static void Yjd3ELEMENT4(
                         INT nelest, INT nelem,
                         INT **i2elto, INT **i2eljp)
{
  INT i,j,ielem;

//  FILE *fp;
//  fp=fopen("JointToFEM.txt","w+");

//  fprintf(fp,"JointToFEM\n");

  for(ielem=nelest;ielem<nelem;ielem++)
    {
      for(j=0;j<2;j++)
        {
          for(i=0;i<nelest;i++)
            {
              if(((i2elto[ielem][0+j*3]==i2elto[i][0])||
                  (i2elto[ielem][0+j*3]==i2elto[i][1])||
                  (i2elto[ielem][0+j*3]==i2elto[i][2])||
                  (i2elto[ielem][0+j*3]==i2elto[i][3]))&&
                 ((i2elto[ielem][1+j*3]==i2elto[i][0])||
                  (i2elto[ielem][1+j*3]==i2elto[i][1])||
                  (i2elto[ielem][1+j*3]==i2elto[i][2])||
                  (i2elto[ielem][1+j*3]==i2elto[i][3]))&&
                 ((i2elto[ielem][2+j*3]==i2elto[i][0])||
                  (i2elto[ielem][2+j*3]==i2elto[i][1])||
                  (i2elto[ielem][2+j*3]==i2elto[i][2])||
                  (i2elto[ielem][2+j*3]==i2elto[i][3])))
                {
                  if(i2eljp[ielem][0]>=0)
                    {
                      i2eljp[ielem][1]=i;
                    }
                  else
                    {
                      i2eljp[ielem][0]=i;
                    }

                  if(i2eljp[i][0]==ielem)
                    {
                      i2eljp[i][0]=-1;

//                      fprintf(fp,"i2eljp[%d][0]=-1\n",i);
                    }
                  else if(i2eljp[i][1]==ielem)
                    {
                      i2eljp[i][1]=-1;

//                      fprintf(fp,"i2eljp[%d][1]=-1\n",i);
                    }
                  else if(i2eljp[i][2]==ielem)
                    {
                      i2eljp[i][2]=-1;

//                      fprintf(fp,"i2eljp[%d][2]=-1\n",i);
                    }
                  else if(i2eljp[i][3]==ielem)
                    {
                      i2eljp[i][3]=-1;

//                    fprintf(fp,"i2eljp[%d][3]=-1\n",i);
                    }
                  else if(i2eljp[i][0]<0)
                    {
                      i2eljp[i][0]=ielem;
                    }
                  else if(i2eljp[i][1]<0)
                    {
                      i2eljp[i][1]=ielem;
                    }
                  else if(i2eljp[i][2]<0)
                    {
                      i2eljp[i][2]=ielem;
                    }
                  else if(i2eljp[i][3]<0)
                    {
                      i2eljp[i][3]=ielem;
                    }

                  continue;
                }
            }
        }
//        fprintf(fp,"i2eljp[%d][0]=%d\ti2eljp[%d][1]=%d\n",ielem,i2eljp[ielem][0],ielem,i2eljp[ielem][1]);
    }

/*  for(ielem=0;ielem<nelest;ielem++)
  {
       for(i=0;i<6;i++)
       {
            fprintf(fp,"FEM_eljp[%d][%d]=%d\n",ielem,i,i2eljp[ielem][i]);
       }
  }

  for(ielem=nelest;ielem<nelem;ielem++)
  {
       for(i=0;i<6;i++)
       {
            fprintf(fp,"JOINT_eljp[%d][%d]=%d\n",ielem,i,i2eljp[ielem][i]);
       }
  }
       
    fclose(fp);
*/
}

/* LINK NODE TO ITS NEIGHBOURING NODES*/
static void Yjd3NODE2NODE(INT nelest, INT nnopo, INT **i2elto, INT **i2eljp,
                          INT *i1elbe,
                          DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
                          INT *n0neigh, INT **i2nnei, INT *i1nei
                          )
{
  INT ielem,jelem,i,j,k,m,n;
  INT *i1eltoi,*i1eltoj;
  INT nmax,jcon;

  for(ielem=0;ielem<nelest;ielem++)
    {
      i1eltoi=i2elto[ielem];
      for(i=0;i<4;i++)
        {
          i2nnei[i1eltoi[i]][0]=ielem; /* i2nnei[node No.][0] - store elem No. where this node is on */
	  nmax=1;

          for(jelem=0;jelem<nelest;jelem++)
            {
              if(ielem!=jelem)
                {
                  i1eltoj=i2elto[jelem];
                  for(j=0;j<4;j++)
                    {
                      if((i1eltoi[i]!=i1eltoj[j])&&
                         ((DABS(d1nccx[i1eltoi[i]]-d1nccx[i1eltoj[j]])<EPSILON)
                         &&(DABS(d1nccy[i1eltoi[i]]-d1nccy[i1eltoj[j]])<EPSILON)
                          &&(DABS(d1nccz[i1eltoi[i]]-d1nccz[i1eltoj[j]])<EPSILON)))
                        {
                          nmax+=1;
			  for(n=1;n<100;n++)
			    {
			      if(i2nnei[i1eltoi[i]][n]<0)
				{
				  i2nnei[i1eltoi[i]][n]=jelem;
				  break;
				}
			    }
                        }
                    }
                }
            }

	  i1nei[i1eltoi[i]]=nmax; /* number of node neighbours, including itself */

          if(nmax>(*n0neigh))
            {
              (*n0neigh)=nmax;
            }
        }
    }
}


/********************* FIND PREXISTING FRACTURE ********** added by Asiri 20/10/2017 *************/
static void Ypref (INT nelemi, INT nelem, DBL *d1ncix, DBL *d1nciy,
       DBL *d1nciz, INT nnopi, INT *i1elty, INT **i2elto)
{
  INT i, j, k;
  INT ielem, jelem;
  INT in0, in1, in2, jn0, jn1, jn2;
  DBL SEP;

  SEP=1.0e-7;


  for (ielem = nelemi; ielem < nelem; ielem++) //nelemi
    {

      if (i1elty[ielem] > 0)
	{
	  // checking joint whether overlaps with another one
	  in0 = i2elto[ielem][0];
	  in1 = i2elto[ielem][1];
	  in2 = i2elto[ielem][2];

	  for (jelem = ielem + 1; jelem < nelem; jelem++) //nenlmi
	    {
	      //  if (i1elty[jelem]>0)
	      //{

	        for (i = 3; i < 6; i++)
	       {
	       j = i + 1;

	       if (j > 5)
	       j = 3;

	       k = j + 1;

	       if (k > 5)
	       k = 3;

	 /*     for (i = 0; i < 3; i++)
		{
		  j = i + 1;

		  if (j > 2)
		    j = 0;

		  k = j + 1;

		  if (k > 2)
		    k = 0;*/

		  jn0 = i2elto[jelem][i];
		  jn1 = i2elto[jelem][j];
		  jn2 = i2elto[jelem][k];

		  /* this does not work properly, need to identify better*/
		  if (DABS(d1ncix[in0]-d1ncix[jn2]) < SEP
		      && DABS(d1ncix[in1]-d1ncix[jn1]) < SEP
		      && DABS(d1ncix[in2]-d1ncix[jn0]) < SEP

		      && DABS(d1nciy[in0]-d1nciy[jn2]) < SEP
		      && DABS(d1nciy[in1]-d1nciy[jn1]) < SEP
		      && DABS(d1nciy[in2]-d1nciy[jn0]) < SEP

		      && DABS(d1nciz[in0]-d1nciz[jn2]) < SEP
		      && DABS(d1nciz[in1]-d1nciz[jn1]) < SEP
		      && DABS(d1nciz[in2]-d1nciz[jn0]) < SEP)
		    {
		      i1elty[ielem] = 2;
		      i1elty[jelem] = 2;
		      break;
		    }
		}
	    }
	}
    }

}
/**************************************************************************/
/* PUBLIC                                                                 */
/**************************************************************************/
	     
void Ymd(YDE yde, YDI ydi, YDN ydn, YDP ydp) /* mesh elements */
{
     INT nelest,nnopst;
     INT iprop;
     INT inopo;
     INT imesh,imestyp;
     INT *i1jnef;   /*joint node element first for each old node*/
     INT *i1jnen;   /*joint node element next for each new node*/
     DBL *d1sdel;
     INT ielem,integ;

     for(imesh=0;imesh<10;imesh++)
     {
	  nelest=yde->nelem;
	  nnopst=ydn->nnopo;
	  i1jnef=INT1NULL;
	  i1jnen=INT1NULL;

	  for(iprop=0;iprop<ydp->nprop;iprop++)
	  {
	       imestyp=ydp->i1pemn[iprop]%10;
	       ydp->i1pemn[iprop]=ydp->i1pemn[iprop]/10;
	       if(ydp->i1ptyp[iprop]==YTE3TET4ELS)
	       {
		    d1sdel=yde->d2elst[ydp->i1psde[iprop]];
		    if(imestyp==3)
		    {
			 if(i1jnef==INT1NULL)
			 {
			      ydi->diedi=ydi->diezon+ydi->diezon;
			      i1jnef=TalINT1(nnopst);
			      i1jnen=TalINT1(4*nelest);
			      for(inopo=0;inopo<nnopst;inopo++)
			      {
				   i1jnef[inopo]=-1;
			      }
			      Yjd3TET4(   /*create joints*/
				   yde->melem, ydn->mnopo, nelest, nnopst,
				   iprop, ydp->i1pejp[iprop],
				   &(yde->nelem), &(ydn->nnopo),
				   ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2ncc[2],
				   ydn->d2nci[0], ydn->d2nci[1], ydn->d2nci[2],
				   ydn->d2nvc[0], ydn->d2nvc[1], ydn->d2nvc[2],
				   d1sdel,
				   yde->i1elpr,
				   i1jnef, i1jnen, ydn->i1nobf, ydn->i1nopr,
				   yde->i2elto,yde->i2eljp,yde->i1elty,
				   ydn->d1ntc, ydn->d1nti);

			      Yjd3JOINTPR( /* distinguish two kinds of joint elements */
					  nelest,yde->nelem,yde->i2elto,yde->i1eljo, yde->i1elty);

                  Yjd3ELEMENT4( /* link joint element to 4-node finite element */
                               nelest,yde->nelem,yde->i2elto,yde->i2eljp);


 
 
   Yjd3NODE2NODE( /* link node to its neighbouring nodes */
		   nelest,ydn->nnopo,yde->i2elto,yde->i2eljp,
		   yde->i1elbe,
		   ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2ncc[2],
		   &(ydn->nneigh),ydn->i2nnei,ydn->i1nei);

			      /*find pre-exisitng fracture */
			      Ypref (yde->nelemi, yde->nelem, ydn->d2nci[0], ydn->d2nci[1], ydn->d2nci[2],
				     ydn->nnopi, yde->i1elty, yde->i2elto);

			 }
		    }
	       }
	  }
          /* free memory */
          FREE(i1jnef);
          FREE(i1jnen);
     }
}

/***********************************************************************************/
/* EOF                                                                             */
/***********************************************************************************/			      
			 
		    
		    
      
