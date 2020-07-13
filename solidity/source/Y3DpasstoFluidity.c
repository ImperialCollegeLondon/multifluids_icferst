/*! \file Y3Dmd.c by LG^M
 *  \mesh elements^M
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
 *
 */


/**********************************************************************/
/**********************************************************************/

#include "Yproto.h"

/* initiate surface mesh on FEM tetrahedra mesh, before inserting joint elemets */
static void YFEMsurface(
			INT ncmel, INT **i2cvto, INT *i1xbtoel,
			INT **i2xbto, INT *n0sfel, INT **i2eltoxb,
			INT nelem_t, INT **i2elto_t, INT nelno
			)
{
  INT ielem,jelem,iflag;
  INT i,j,k,m,in,jn,kn;
  INT itar,jtar,ktar,itarn,jtarn,ktarn;
  INT *i1elto,itar_tmp,jelem_tmp;

  //  FILE *fp;
  //  fp=fopen("fem_surface.txt","w+");


  for(ielem=0;ielem<nelem_t;ielem++)
    {
      i1elto=i2elto_t[ielem];


	  in=i1elto[0];
	  jn=i1elto[1];
	  kn=i1elto[2];
      i2xbto[ielem][0]=in;
      i2xbto[ielem][1]=jn;
      i2xbto[ielem][2]=kn;
        
	  iflag=0;

	  for(jelem=0;jelem<ncmel;jelem++)
	    {

		  for(itar=0;itar<4;itar++)
		    {
		      jtar=itar+1;
		      if(jtar>3)
			{
			  jtar=0;
			}
		      ktar=jtar+1;
		      if(ktar>3)
			{
			  ktar=0;
			}
		      itarn=i2cvto[jelem][itar];
		      jtarn=i2cvto[jelem][jtar];
		      ktarn=i2cvto[jelem][ktar];

		      if(((in==itarn)&&(jn==jtarn)&&(kn==ktarn))||((in==jtarn)&&(jn==ktarn)&&(kn==itarn))||((in==ktarn)&&(jn==itarn)&&(kn==jtarn))
			 ||((in==ktarn)&&(jn==jtarn)&&(kn==itarn))||((in==jtarn)&&(jn==itarn)&&(kn==ktarn))||((in==itarn)&&(jn==ktarn)&&(kn==jtarn)))
			{
			  iflag=1; /* surface triangle shared by two neighbouring tetrahedra */
                itar_tmp=itar;
                jelem_tmp=jelem;
			  break;
			}
		    }
		  if(iflag==1)
		    {
		      break;
		    }
		}

	  if(iflag==1)		//surface triangle on boundary
	    {
	      /* guarantee tetrahedron volume positive */
	      

	      i1xbtoel[ielem]=jelem_tmp;

	      i2eltoxb[jelem][itar_tmp]=ielem;
	      		  
	    }
        else
        {
            CHRw(stdout,"surface match wrong");
             CHRwcr(stdout); INTw(stdout,ielem,10); 
            exit(1);
        }
    }

  (*n0sfel)=nelem_t;
    printf("surface element number is %ld\n",*n0sfel);
  /*  
      for(i=0;i<12970;i++)
      {
      fprintf(fp,"i2xbto[%ld][0]=%ld i2xbto[%ld][1]=%ld i2xbto[%ld][2]=%ld i2xbto[%ld][3]=%ld i2xbto[%ld][4]=%ld i2xbto[%ld][5]=%ld\n",
      i,i2xbto[i][0],i,i2xbto[i][1],i,i2xbto[i][2],i,i2xbto[i][3],i,i2xbto[i][4],i,i2xbto[i][5]);
      }
  
      fclose(fp);
  */
}

/* link all continuum nodes, including both surface and volume, to joint elem nodes */
static void YlinkCONtoJOINT(
			    INT **i2xbjno, INT *i1xbjo,
			    INT nnopo, DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
			    INT ncmno, INT ncmel, INT *i1ptyp, INT *i1pejp, INT *i1elpr,
			    INT **i2elto, INT **i2cvto, INT nnopi
			    )
{
  INT inopo,jnopo,ijoint,i;
  INT ielem;
  INT *i1xbto;

  //FILE *fp1;
  //fp1=fopen("surface_node.txt","w+");


  printf("----------------------------------------------------------- ncmno=%i nnopo=%i \n\n\n\n\n\n", ncmno,nnopo);

  for(inopo=0;inopo<ncmno;inopo++)
    {
      ijoint=0;
      //i1xbjo[inopo]=ijoint;
     // for(jnopo=0;jnopo<nnopo;jnopo++)
      for(jnopo=ncmno;jnopo<nnopo;jnopo++)
	{
	  if((DABS(d1nccx[inopo]-d1nccx[jnopo])<1.0e-6)&&(DABS(d1nccy[inopo]-d1nccy[jnopo])<1.0e-6)&&(DABS(d1nccz[inopo]-d1nccz[jnopo])<1.0e-6))
	    {
	      i2xbjno[inopo][ijoint]=jnopo;
	      ijoint=ijoint+1;
	      //i1xbjo[inopo]=ijoint;
	    }
	}

      if(i2xbjno[inopo][0]<0)
	{
	  i2xbjno[inopo][0]=inopo;
	}
    }
  /*
    for(ielem=0;ielem<nelem;ielem++)
    {
    //if(i1ptyp[i1elpr[ielem]]==(YTE3TET4ELS))
    if((i1ptyp[i1elpr[ielem]]==(YTE3TET4ELS))&&(i1pejp[i1elpr[ielem]]==0))
    {
    for(i=0;i<4;i++)
    {
    inopo=i2elto[ielem][i];
    i2xbjno[inopo][0]=inopo;		
    }
    }
    }
  */
  /*
    printf("start checking nodes ...\n");
    for(ielem=0;ielem<ncmel;ielem++)
    {
    for(i=0;i<4;i++)
    {
    if((DABS(d1nccx[i2cvto[ielem][i]]-d1nccx[i2elto[ielem][i]])>=1.0e-6)||(DABS(d1nccy[i2cvto[ielem][i]]-d1nccy[i2elto[ielem][i]])>=1.0e-6)||(DABS(d1nccz[i2cvto[ielem][i]]-d1nccz[i2elto[ielem][i]])>=1.0e-6))
    {
    printf("WRONG relation: elem%ldnode%ld\n",ielem,i);
    }
    }
    }
    printf("nodes checking finished.\n");
  */
  /*  
      for(i=0;i<18700;i++)
      {
      fprintf(fp1,"i2xbjno[%ld][0]=%ld i2xbjno[%ld][1]=%ld i2xbjno[%ld][2]=%ld i2xbjno[%ld][3]=%ld i2xbjno[%ld][4]=%ld i2xbjno[%ld][5]=%ld i2xbjno[%ld][6]=%ld i2xbjno[%ld][7]=%ld i2xbjno[%ld][8]=%ld i2xbjno[%ld][9]=%ld\n",
      i,i2xbjno[i][0],i,i2xbjno[i][1],i,i2xbjno[i][2],i,i2xbjno[i][3],i,i2xbjno[i][4],
      i,i2xbjno[i][5],i,i2xbjno[i][6],i,i2xbjno[i][7],i,i2xbjno[i][8],i,i2xbjno[i][9]);
      }
  */  
  //fclose(fp1);
}

/* link surface continuum nodes, to surface joint nodes */
static void YlinkSFCONtoSFJOINT(
				INT nsfel, INT **i2xbto, INT **i2xbjno, INT **i2xbtosf,
				INT nnopo, INT nelem, INT *i1ptyp, INT *i1elpr,
				INT **i2elto
				)
{
  INT isf,isfno,ielem,i,j,k,m;
  INT jono,conno,flag;
  INT *i1sfjono;

  i1sfjono=INT1NULL;
  i1sfjono=TalINT1(nnopo);

  isfno=0;
  for(ielem=0;ielem<nelem;ielem++)
    {
      if(i1ptyp[i1elpr[ielem]]==(YTE3TET4JOINT))
	{
          //printf("ielem %ld\n",ielem);
	  for(i=0;i<6;i++)
	    {
	      jono=i2elto[ielem][i];
	      for(isf=0;isf<nsfel;isf++)
		{
		  for(j=0;j<3;j++)
		    {
		      conno=i2xbto[isf][j];
		      for(k=0;k<100;k++)
			{
			  if(i2xbjno[conno][k]==jono)
			    {
			      flag=0;
			      for(m=0;m<isfno;m++)
				{
				  if(i1sfjono[m]==jono)
				    {
				      flag=1;
				      break;
				    }
				}
			      if(flag==0)
				{
				  i1sfjono[isfno]=jono;
				  isfno=isfno+1;
				  //printf("nnopo %ld, isfno %ld\n",nnopo,isfno);
				  break;
				}
			    }
			}
		    }
		}		
	    }
	}
    }

  for(isf=0;isf<nsfel;isf++)
    {
      for(i=0;i<3;i++)
	{
	  conno=i2xbto[isf][i];
	 
	  if(i2xbjno[conno][0]==conno) 
	    {
	      i2xbtosf[isf][i]=conno;
	    }
	  else		
	    {
	      for(j=0;j<100;j++)
		{
		  for(k=0;k<isfno;k++)
		    {
		      if(i2xbjno[conno][j]==i1sfjono[k])
			{
			  i2xbtosf[isf][i]=i2xbjno[conno][j];
			  break;
			}
		    }
		  if(i2xbtosf[isf][i]>=0)
		    {
		      break;
		    }
		}
	    }
	}
    }
	  		
  FREE(i1sfjono);
  /*
    for(isf=0;isf<nsfel;isf++)
    {
    printf("i2xbtosf[%ld][] %ld %ld %ld\n",isf,i2xbtosf[isf][0],i2xbtosf[isf][1],i2xbtosf[isf][2]);
    }
  */
}

/* connect fem elem with its neighbouring fem elem */
static void YfindFEMcon(
			INT *i1ptyp, INT ncmel, INT **i2eljp, INT *i1elpr, INT **i2elcon
			)
{
  INT ielem,jono,i,j;

  for(ielem=0;ielem<ncmel;ielem++)
    {
      if(i1ptyp[i1elpr[ielem]]==(YTE3TET4ELS))
	{
	  for(i=0;i<4;i++)
	    {
	      jono=i2eljp[ielem][i]; /* joint elem number */
	      if(jono<0)
		{
		  i2elcon[ielem][i]=jono;
		}
	      else
		{
		  for(j=0;j<2;j++)
		    {
		      if(i2eljp[jono][j]!=ielem)
			{
			  i2elcon[ielem][i]=i2eljp[jono][j];
			  break;
			}
		    }
		}
	    }
	}
    }
}

/* link surface triangle element number to discontinuum node number */
static void YlinkSFtoJOINT(
			   INT nsfel, INT **i2xbto, INT **i2xbtojn, INT *i1xbtoel,
			   INT **i2xbjno,
			   INT **i2elto,
			   DBL *d1ncix, DBL *d1nciy, DBL *d1nciz
			   )
{
  INT isf,flag;
  INT i,j,k;
  INT isfel,icono,jdisno;

  for(isf=0;isf<nsfel;isf++)
    {
      isfel=i1xbtoel[isf];   //surface tri connected tet elem number
      for(i=0;i<3;i++)
	{
	  icono=i2xbto[isf][i];   //surface tri continuum node number

	  for(j=0;j<4;j++)   //four nodes of connected tet element
	    {
	      jdisno=i2elto[isfel][j];
	      if((DABS(d1ncix[icono]-d1ncix[jdisno])<1.0e-6)&&(DABS(d1nciy[icono]-d1nciy[jdisno])<1.0e-6)&&(DABS(d1nciz[icono]-d1nciz[jdisno])<1.0e-6))
		{
		  i2xbtojn[isf][i]=jdisno;

		  break;
		}
	    }
	}
    }


}

/* update neighbour continuum tetrahedron element */
static void YFEMneiUpdate(
			  INT sfcvno, INT sfcvnonew, INT femup,
			  INT *i1ptyp, INT nelem, INT *i1elpr, INT **i2elcon, INT **i2cvto,
			  INT **i2xbto, INT **i2eltoxb,
			  INT **i2xbjno, INT **i2elto
			  )
{
  INT inei,inode,i,j,flag;
  INT neielem;
  
  //FILE *fp1;
  //fp1=fopen("nei.txt","a+");

  for(inei=0;inei<4;inei++) /* 4 neighbour tetrahedra elements */
    {
      //if(sfcvno==6253)
      //{
      //printf("femorig%ldnei %ld %ld %ld %ld\n",femup,i2elcon[femup][0],i2elcon[femup][1],i2elcon[femup][2],i2elcon[femup][3]);
      //}
      if(i2elcon[femup][inei]>=0)
	{
	  neielem=i2elcon[femup][inei];
	  for(inode=0;inode<4;inode++)
	    {
	      if(i2cvto[neielem][inode]==sfcvno) /* this node on neighbouring continuum tetrahedron equals the updated node */
		{
		  //if(sfcvno==6253)
		  //{
		  //printf("femorig%ldnei%ldneielem%ldnode%ld %ld %ld %ld %ld\n",femup,inei,neielem,inode,i2cvto[neielem][0],i2cvto[neielem][1],i2cvto[neielem][2],i2cvto[neielem][3]);
		  //}

		  i2cvto[neielem][inode]=sfcvnonew; /* update this node to new number */

		  //if(sfcvno==6253)
		  //{
		  //printf("femorig%ldnei%ldneielem%ldnode%ld %ld %ld %ld %ld\n",femup,inei,neielem,inode,i2cvto[neielem][0],i2cvto[neielem][1],i2cvto[neielem][2],i2cvto[neielem][3]);
		  //}

		  /* add joint elem nodes corresponding to new continuum node */
		  flag=0;
		  for(i=0;i<100;i++) 
		    {
		      if(i2xbjno[sfcvnonew][i]==i2elto[neielem][inode])
			{
			  flag=1;
			  break;
			}
		    }
		  if(flag==0)
		    {
		      for(i=0;i<100;i++)
			{
			  if(i2xbjno[sfcvnonew][i]<0)
			    {
			      i2xbjno[sfcvnonew][i]=i2elto[neielem][inode];
			      break;
			    }
			}
		    }

		  /* delete joint elem nodes corresponding to old continuum node */
		  for(i=0;i<100;i++)
		    {
		      if(i2xbjno[sfcvno][i]==i2elto[neielem][inode])
			{
			  i2xbjno[sfcvno][i]=-1;
			  break;
			}	 
		    }

		  for(i=0;i<4;i++)	/* 4 surface triange elements */
		    {
		      if(i2eltoxb[neielem][i]>=0)
			{
			  for(j=0;j<3;j++)
			    {
			      if(i2xbto[i2eltoxb[neielem][i]][j]==sfcvno)
				{
				  i2xbto[i2eltoxb[neielem][i]][j]=sfcvnonew;
				}
			    }
			}
		    }

		  YFEMneiUpdate(
				sfcvno,sfcvnonew,neielem,
				i1ptyp,nelem,i1elpr,i2elcon,i2cvto,
                                i2xbto,i2eltoxb,
				i2xbjno,i2elto);

		  //break;
		}
	    }
	}
    }
   
  //fclose(fp1);
}

/* check connectivity of neighbour continuum tetrahedron element */
static void YFEMnei(
		    INT sfcvno, INT sfcvnonew, INT femup, INT femold, INT *femflag,
		    INT *i1ptyp, INT nelem, INT *i1elpr, INT **i2elcon, INT **i2cv0to,
		    INT **i2xbto, INT **i2eltoxb,
		    INT **i2xbjno, INT **i2elto
		    )
{
  INT inei,inode,i,j,flag;
  INT neielem;
  
  //FILE *fp1;
  //fp1=fopen("nei.txt","a+");

  //printf("1 femup %ld, femold %ld, femflag %ld\n",femup,femold,*femflag);

  for(inei=0;inei<4;inei++) /* 4 neighbour tetrahedra elements */
    {
      //if(sfcvno==6253)
      //{
      //printf("femorig%ldnei %ld %ld %ld %ld\n",femup,i2elcon[femup][0],i2elcon[femup][1],i2elcon[femup][2],i2elcon[femup][3]);
      //}
      if(i2elcon[femup][inei]>=0)
	{
	  neielem=i2elcon[femup][inei];
	  for(inode=0;inode<4;inode++)
	    {
	      if(i2cv0to[neielem][inode]==sfcvno) /* this node on neighbouring continuum tetrahedron equals the updated node */
		{
		  //if(sfcvno==6253)
		  //{
		  //printf("femorig%ldnei%ldneielem%ldnode%ld %ld %ld %ld %ld\n",femup,inei,neielem,inode,i2cvto[neielem][0],i2cvto[neielem][1],i2cvto[neielem][2],i2cvto[neielem][3]);
		  //}

		  i2cv0to[neielem][inode]=sfcvnonew; /* update this node to new number */

		  //if(sfcvno==6253)
		  //{
		  //printf("femorig%ldnei%ldneielem%ldnode%ld %ld %ld %ld %ld\n",femup,inei,neielem,inode,i2cvto[neielem][0],i2cvto[neielem][1],i2cvto[neielem][2],i2cvto[neielem][3]);
		  //}

		  /* add joint elem nodes corresponding to new continuum node */
		  /*
		    flag=0;
		    for(i=0;i<100;i++) 
		    {
		    if(i2xbjno[sfcvnonew][i]==i2elto[neielem][inode])
		    {
		    flag=1;
		    break;
		    }
		    }
		    if(flag==0)
		    {
		    for(i=0;i<100;i++)
		    {
		    if(i2xbjno[sfcvnonew][i]<0)
		    {
		    i2xbjno[sfcvnonew][i]=i2elto[neielem][inode];
		    break;
		    }
		    }
		    }
		  */
		  /* delete joint elem nodes corresponding to old continuum node */
		  /*
		    for(i=0;i<100;i++)
		    {
		    if(i2xbjno[sfcvno][i]==i2elto[neielem][inode])
		    {
		    i2xbjno[sfcvno][i]=-1;
		    break;
		    }	 
		    }

		    for(i=0;i<4;i++)	// 4 surface triange elements
		    {
		    if(i2eltoxb[neielem][i]>=0)
		    {
		    for(j=0;j<3;j++)
		    {
		    if(i2xbto[i2eltoxb[neielem][i]][j]==sfcvno)
		    {
		    i2xbto[i2eltoxb[neielem][i]][j]=sfcvnonew;
		    }
		    }
		    }
		    }
		  */
                  if(neielem==femold)
		    { 
		      (*femflag)=1;
		      //printf("2 femup %ld, femold %ld, femflag %ld, neielem %ld\n",femup,femold,*femflag,neielem);
		    }
		  YFEMnei(
			  sfcvno,sfcvnonew,neielem,femold,femflag,
			  i1ptyp,nelem,i1elpr,i2elcon,i2cv0to,
			  i2xbto,i2eltoxb,
			  i2xbjno,i2elto);

		  break;
		}
	    }
	}
    }
   
  //fclose(fp1);
}

/* update unseparated edge info */
static void UpdateEdge(
                       INT node0, INT node1, INT nelem,
                       INT **i2edge
                       )
{ 
  INT i,j,k;
  INT ielem;

  for(i=0;i<nelem;i++)
    {
      if(i2edge[i][0]>=0)
        {
          if(((i2edge[i][0]==node0)&&(i2edge[i][1]==node1))||((i2edge[i][1]==node0)&&(i2edge[i][0]==node1)))
            {
              ielem=i2edge[i][2];
              for(j=0;j<3;j++)
                {
                  i2edge[i][j]=-1;
                }
              for(k=0;k<nelem;k++)
                {
                  if(i2edge[k][2]==ielem)
                    {
                      for(j=0;j<3;j++)
                        {
                          i2edge[k][j]=-1;
                        }
                      break;
                    }
                }
            }
        }
    }
}

/* add new fracture surface and node to the continuum surface mesh */
static void YSFMeshUpdate(
			  INT nprop, INT *i1ptyp, INT nelem, INT *i1elpr, INT **i2elto,
			  INT *n0sfel, INT **i2xbto, INT **i2xbtojn, INT *i1xbtoel,
			  INT *i1xbjo, INT **i2xbjno,
			  INT **i2eljp, INT *n0cmno, INT **i2cvto, INT **i2eltoxb,
			  DBL *d1ncix, DBL *d1nciy, DBL *d1nciz,
			  INT **i2elcon, INT ncmel
			  )
{
  INT ielem,isf,i,j,k,m,n,p,q,pflag[3];
  INT jointn,nsfel,ncmno,sfn;
  INT joonsf,edonsf,femup,femold,femuppair,femoldpair,fempre,femnei,flag,temp,neiflag;
  INT femup2old,femold2up,fempairup2old,fempairold2up,femflag;
  INT sfcvnonew[4],sfcvno[4],jn[3],edno[3];
  INT internode,shno[2],shjono[2],unedg[2][2],jo2cn[3][4],jocved[3][3];
  INT **i2edge,**i2cv0to;
  DBL F0[3][3],voli;
  INT i1elto[4];
  INT oldsfn,newsfn,flag2,flag3;
  INT *i1edsf,ied0,ied1;
  INT node[2];

  //INT imesh;

//  FILE *fp;
//  fp=fopen("info.txt","a+");

  i2edge=INT2NULL;
  i2edge=TalINT2(nelem,3);

  i2cv0to=INT2NULL;
  i2cv0to=TalINT2(ncmel,4);

  i1edsf=INT1NULL;
  i1edsf=TalINT1(50);

  for(i=0;i<nelem;i++)
    {
      for(j=0;j<3;j++)
	{
	  i2edge[i][j]=-1;
	}
    }

  nsfel=(*n0sfel);
  ncmno=(*n0cmno);

  //for(imesh=0;imesh<10;imesh++)
  //{
  for(ielem=0;ielem<nelem;ielem++)  //changed from nelem to nelemi (this might need to go from nelest to nelem)
    {


/*
      if(((i2eljp[ielem][0]==43463)&&(i2eljp[ielem][1]==60914))||((i2eljp[ielem][1]==43463)&&(i2eljp[ielem][0]==60914)))
	{
	  printf("JOINT elem %ld\n",ielem);
	  printf("i2eljp[%ld][0] %ld, i2eljp[%ld][1] %ld\n",ielem,i2eljp[ielem][0],ielem,i2eljp[ielem][1]);
	}

*/
      if(((i1elpr[ielem]+YIPROPMAX)>=0)&&((i1elpr[ielem]+YIPROPMAX)<nprop))
	{
	  if(i1ptyp[i1elpr[ielem]+YIPROPMAX]==(YTE3TET4JOINT)) /* broken joint element (YIPROPMAX=1000) */
	    {
//	      fprintf(fp,"*YSFMeshUpdate* failed elem %ld\n",ielem);       

	      //edonsf=0;

	      for(i=0;i<3;i++)
		{
		  jn[i]=i2elto[ielem][i]; /* broken joint elem node number 0-2 */
		}
	      for(i=0;i<4;i++)
		{
		  sfcvno[i]=-1;
		  sfcvnonew[i]=-1;
		}

	      //printf("nsfel %ld, ncmno %ld\n",nsfel,ncmno);
	      /* find continuum node corresponding to joint elem node jn[i] */
	      for(i=0;i<3;i++)
		{
		  jo2cn[i][0]=-1;
		  for(j=0;j<ncmno;j++)
		    {
		      for(k=0;k<100;k++)
			{
			  if(jn[i]==i2xbjno[j][k])
			    {
			      jo2cn[i][0]=j; /* continuum node corresponding to joint elem node jn[i] */
			      break;
			    }
			}
		      if(jo2cn[i][0]>=0)
			{
			  break;
			}
		    }
		}
	      //printf("jo2cn[-][0] %ld %ld %ld\n",jo2cn[0][0],jo2cn[1][0],jo2cn[2][0]);

	      /* form joint elem node corresponding to continuum node */
	      m=2; 
	      for(i=0;i<4;i++)
		{
		  if(i2elto[ielem][0]==i2elto[i2eljp[ielem][0]][i])
		    {
		      m=1;
		      break;
		    }
		}
	      if(m==1)
		{
		  n=2;
		}
	      else
		{
		  n=1;
		}
	      /* jo2cn[][1] - i2eljp[][0], jo2cn[][2] - i2eljp[][1] */
	      jo2cn[0][m]=i2elto[ielem][0];
	      jo2cn[1][m]=i2elto[ielem][1];
	      jo2cn[2][m]=i2elto[ielem][2];
	      jo2cn[0][n]=i2elto[ielem][5];
	      jo2cn[1][n]=i2elto[ielem][4];
	      jo2cn[2][n]=i2elto[ielem][3];
	      for(i=0;i<3;i++)
		{
		  jo2cn[i][3]=0;
		}

	      //printf("jo2cn[-][1] %ld %ld %ld\n",jo2cn[0][1],jo2cn[1][1],jo2cn[2][1]);
	      //printf("jo2cn[-][2] %ld %ld %ld\n",jo2cn[0][2],jo2cn[1][2],jo2cn[2][2]);
	      //printf("jo2cn[-][3] %ld %ld %ld\n",jo2cn[0][3],jo2cn[1][3],jo2cn[2][3]);
	    
	      /* form joint elem edge (continuum mesh) */
	      for(i=0;i<3;i++)
		{
		  j=i+1;
		  if(j>2)
		    {
		      j=0;
		    }
		  jocved[i][0]=jo2cn[i][0];
		  jocved[i][1]=jo2cn[j][0];
		  jocved[i][2]=-1;	/* -1 - not on surface, 1 - on surface, -2 - one node on one surface triangle, one node on another surface triangle */
		}

	      //printf("jocved[-][0] %ld %ld %ld\n",jocved[0][0],jocved[1][0],jocved[2][0]);
	      //printf("jocved[-][1] %ld %ld %ld\n",jocved[0][1],jocved[1][1],jocved[2][1]);
	      //printf("jocved[-][2] %ld %ld %ld\n",jocved[0][2],jocved[1][2],jocved[2][2]);

	      /* whether joint elem edge is on surface boundary */
	      edonsf=0;		/* number of joint element edges on continuum surface */
	      for(i=0;i<3;i++)
		{
		  ied0=0;
                  ied1=0;
		  for(j=0;j<50;j++)
		    {
		      i1edsf[j]=-1;
		    }

		  for(isf=0;isf<nsfel;isf++) /* first node of two edge nodes */
		    {
		      for(k=0;k<3;k++)
			{
			  if(jocved[i][0]==i2xbto[isf][k])
			    {
			      i1edsf[ied0]=isf;
			      ied0+=1;
			      break;
			    }
			}
		    }

		  if(ied0>0)
		    {
		      for(isf=0;isf<nsfel;isf++) /* second node of two edge nodes */
			{
			  for(k=0;k<3;k++)
			    {
			      if(jocved[i][1]==i2xbto[isf][k])
				{
				  ied1+=1;
				  for(m=0;m<ied0;m++)
				    {
				      if(isf==i1edsf[m])
					{
					  jocved[i][2]=1;
					  break;
					}
				    }
				  if(jocved[i][2]==1)
				    {
				      break;
				    }
				}
			    }
			  if(jocved[i][2]==1)
			    {
			      break;
			    }
			}
		    }

		  if((jocved[i][2]==(-1))&&(ied0>0)&&(ied1>0))
		    {
		      jocved[i][2]=-2;
		    }
	      
		  if(jocved[i][2]==1)
		    {
		      edonsf+=1;
		    }
		}

	      for(i=0;i<3;i++)
		{
		  j=i+1;
		  if(j>2)
		    {
		      j=0;
		    }
		  if(jocved[i][2]==1)
		    {
		      jo2cn[i][3]=jo2cn[i][3]+1;
		      jo2cn[j][3]=jo2cn[j][3]+1;
		    }
		}
                
	      //printf("***YSFMeshUpdate***\n");
//	      fprintf(fp,"edonsf %ld\n",edonsf);

	      if(edonsf==3)	/* all three joint elem edges on continuum surface, separate all three nodes */
		{
		  //printf("edonsf %ld\n",edonsf);
	       
//		  fprintf(fp,"edonsf3_elem %ld\n",ielem);
	      
		  femup=i2eljp[ielem][1]; /* updated fem elem number, in this case always update fem elem [1], keep fem elem [0] as original */
		  femold=i2eljp[ielem][0];

		  for(i=0;i<4;i++) /* set new boundary label -1 to surface */
		    {
		      if(i2elcon[femup][i]==femold)
			{
			  i2elcon[femup][i]=-1;
			  femup2old=i;
			  break;
			}
		    }
		  for(i=0;i<4;i++)
		    {
		      if(i2elcon[femold][i]==femup)
			{
			  i2elcon[femold][i]=-1;
			  femold2up=i;
			  break;
			}
		    }

		  for(i=0;i<3;i++) /* update all three nodes */
		    {
		      sfcvno[i]=jo2cn[i][0]; /* original number of updated node on continuum mesh */
		      sfcvnonew[i]=ncmno; /* new number of updated node on continuum mesh */
//		      fprintf(fp,"oldnode %ld, newnode %ld\n",sfcvno[i],sfcvnonew[i]);
		      ncmno=ncmno+1;
		    }

                  for(i=0;i<ncmel;i++)
		    {
		      for(j=0;j<4;j++)
                        {
                          i2cv0to[i][j]=i2cvto[i][j];
                        }
		    }
                  femflag=0;
                  for(i=0;i<3;i++)
                    {
		      YFEMnei(
			      sfcvno[i],sfcvnonew[i],femup,femold,&femflag,
			      i1ptyp,nelem,i1elpr,i2elcon,i2cv0to,
			      i2xbto,i2eltoxb,
			      i2xbjno,i2elto);
		    }
                               
		  if(femflag==1)
		    {
//		      fprintf(fp,"elem%ld edonsf3 femflag %ld femup %ld femold %ld\n",ielem,femflag,femup,femold);
		      i2elcon[femup][femup2old]=femold;
		      i2elcon[femold][femold2up]=femup;
                      ncmno=ncmno-3;
		      //printf("elem%ldi2elcon reset.\n",ielem);
		      continue;
		    }

		  //printf("i2xbjno[%ld][0]=%ld, ",sfcvnonew[0],i2xbjno[sfcvnonew[0]][0]);
		  //printf("i2xbjno[%ld][1]=%ld\n",sfcvnonew[0],i2xbjno[sfcvnonew[0]][1]);
		  //printf("i2xbjno[%ld][0]=%ld, ",sfcvnonew[1],i2xbjno[sfcvnonew[1]][0]);
		  //printf("i2xbjno[%ld][1]=%ld\n",sfcvnonew[1],i2xbjno[sfcvnonew[1]][1]);
		  //printf("i2xbjno[%ld][0]=%ld, ",sfcvnonew[2],i2xbjno[sfcvnonew[2]][0]);
		  //printf("i2xbjno[%ld][1]=%ld\n",sfcvnonew[2],i2xbjno[sfcvnonew[2]][1]);

		  /* update continuum tetrahedron elem topology */
		  for(i=0;i<3;i++)
		    {
		      for(j=0;j<4;j++)
			{ 
			  if((i==0)&&(i2cvto[femup][j]!=sfcvno[0])&&(i2cvto[femup][j]!=sfcvno[1])&&(i2cvto[femup][j]!=sfcvno[2]))
			    {
			      sfcvnonew[3]=i2cvto[femup][j];
			    }
			  if(i2cvto[femup][j]==sfcvno[i])
			    {     
			      i2cvto[femup][j]=sfcvnonew[i];                                 
			    }                    
			}
		    }

		  /* update continuum tetrahedron elem connected surface triangle elem topology */
		  for(i=0;i<4;i++)
		    {
		      if(i2eltoxb[femup][i]>=0)
			{
			  for(j=0;j<3;j++)
			    {
			      for(k=0;k<3;k++)
				{
				  if(i2xbto[i2eltoxb[femup][i]][j]==sfcvno[k])
				    {
				      i2xbto[i2eltoxb[femup][i]][j]=sfcvnonew[k];
				    }
				}
			    }
			}
		    }
				  
		  /* add 1st new surface triangle element, using three new nodes */
		  for(i=0;i<3;i++) 
		    {
		      i2xbto[nsfel][i]=sfcvnonew[i];
                      i2xbtojn[nsfel][i]=jo2cn[i][2];
		    }
//		  fprintf(fp,"add surface[femup]%ld: %ld, %ld, %ld\n",nsfel,i2xbto[nsfel][0],i2xbto[nsfel][1],i2xbto[nsfel][2]);
//		  fprintf(fp,"i2eltoxb[%ld]: [0]%ld [1]%ld [2]%ld [3]%ld\n",femup,i2eltoxb[femup][0],i2eltoxb[femup][1],i2eltoxb[femup][2],i2eltoxb[femup][3]);

		  for(i=0;i<4;i++)
		    {
		      if(i2eltoxb[femup][i]<0)
			{
			  i2eltoxb[femup][i]=nsfel;
//			  fprintf(fp,"i2eltoxb[%ld][%ld]=%ld\n",femup,i,i2eltoxb[femup][i]);
			  break;
			}
		    }

                  i1xbtoel[nsfel]=femup;
			  
		  /* add 2nd new surface triangle element, using three old nodes */
		  i2xbto[nsfel+1][0]=sfcvno[2];
		  i2xbto[nsfel+1][1]=sfcvno[1];
		  i2xbto[nsfel+1][2]=sfcvno[0];

                  i2xbtojn[nsfel+1][0]=jo2cn[2][1];
                  i2xbtojn[nsfel+1][1]=jo2cn[1][1];
                  i2xbtojn[nsfel+1][2]=jo2cn[0][1];

//		  fprintf(fp,"add surface[femold]%ld: %ld, %ld, %ld\n",nsfel+1,i2xbto[nsfel+1][0],i2xbto[nsfel+1][1],i2xbto[nsfel+1][2]);
//		  fprintf(fp,"i2eltoxb[%ld]: [0]%ld [1]%ld [2]%ld [3]%ld\n",femold,i2eltoxb[femold][0],i2eltoxb[femold][1],i2eltoxb[femold][2],i2eltoxb[femold][3]);
		  for(i=0;i<4;i++)
		    {
		      if(i2eltoxb[femold][i]<0)
			{
			  i2eltoxb[femold][i]=nsfel+1;
//			  fprintf(fp,"i2eltoxb[%ld][%ld]=%ld\n",femold,i,i2eltoxb[femold][i]);
			  break;
			}
		    }

                  i1xbtoel[nsfel+1]=femold;

		  /* guarantee tetrahedron volume positive */
		  for(i=0;i<3;i++)
		    {
		      i1elto[i]=sfcvno[i];
		    }
		  i1elto[3]=sfcvnonew[3];
		  for(i=1;i<4;i++)
		    { 
		      F0[0][i-1]=d1ncix[i1elto[i]]-d1ncix[i1elto[0]]; 
		      F0[1][i-1]=d1nciy[i1elto[i]]-d1nciy[i1elto[0]];
		      F0[2][i-1]=d1nciz[i1elto[i]]-d1nciz[i1elto[0]];
		    }  
		  voli=F0[0][0]*(F0[1][1]*F0[2][2]-F0[1][2]*F0[2][1])-   
		    F0[0][1]*(F0[1][0]*F0[2][2]-F0[1][2]*F0[2][0])+    
		    F0[0][2]*(F0[1][0]*F0[2][1]-F0[1][1]*F0[2][0]);
		  if(voli<0)    /*modify by JXiang 14032018  */
		    {
		      for(i=0;i<2;i++)
			{
			  temp=i2xbto[nsfel+i][0];
			  i2xbto[nsfel+i][0]=i2xbto[nsfel+i][2];
			  i2xbto[nsfel+i][2]=temp;

                          temp=i2xbtojn[nsfel+i][0];
                          i2xbtojn[nsfel+i][0]=i2xbtojn[nsfel+i][2];
                          i2xbtojn[nsfel+i][2]=temp;
			}
		    }
	
		  nsfel=nsfel+2;

		  for(i=0;i<3;i++)
		    {		  
		      flag=0;
		      for(k=0;k<100;k++) /* add joint elem nodes corresponding to new continuum node */
			{
			  if(i2xbjno[sfcvnonew[i]][k]==jo2cn[i][2]) /* already in the list */
			    {
			      flag=1;
			      break;
			    }
			}
		      if(flag==0)
			{
			  for(k=0;k<100;k++)
			    {			  
			      if(i2xbjno[sfcvnonew[i]][k]<0)
				{
				  i2xbjno[sfcvnonew[i]][k]=jo2cn[i][2];
				  //printf("sfcvnonew[i] %ld %ld\n",sfcvnonew[i],k);
				  break;
				}
			    }
			}

		      for(k=0;k<100;k++) /* delete joint elem nodes corresponding to old continuum node */
			{
			  if(i2xbjno[sfcvno[i]][k]==jo2cn[i][2])
			    {
			      i2xbjno[sfcvno[i]][k]=-1;
			      break;
			    }	 
			}
		    }
		  
		  for(i=0;i<3;i++) /* update con fem elem topology connected with updated elem*/
		    {
		      YFEMneiUpdate(
				    sfcvno[i],sfcvnonew[i],femup,
				    i1ptyp,nelem,i1elpr,i2elcon,i2cvto,
				    i2xbto,i2eltoxb,
				    i2xbjno,i2elto);
		    }

		  i1elpr[ielem]=YIPROPMAX*(-2);
		}
	      else if(edonsf==2)		/* two joint element edges on continuum surface, separate intersection */
		{
		  //printf("edonsf %ld\n",edonsf);

//		  fprintf(fp,"edonsf2_elem %ld\n",ielem);

		  femup=i2eljp[ielem][1]; /* updated fem elem number, in this case always update fem elem [1], keep fem elem [0] as original */
		  femold=i2eljp[ielem][0];

		  for(i=0;i<4;i++) /* set new boundary label -1 to surface */
		    {
		      if(i2elcon[femup][i]==femold)
			{
			  i2elcon[femup][i]=-1;
                          femup2old=i;
			  break;
			}
		    }
		  for(i=0;i<4;i++)
		    {
		      if(i2elcon[femold][i]==femup)
			{
			  i2elcon[femold][i]=-1;
                          femold2up=i;
			  break;
			}
		    }

                  for(i=0;i<2;i++)
                    {
                      node[i]=-1;
                    }

		  for(i=0;i<3;i++) /* update intersection, duplicate the other two */
		    {
		      if(jo2cn[i][3]==2)
			{
			  sfcvno[i]=jo2cn[i][0]; /* original number of updated node on continuum mesh */
			  sfcvnonew[i]=ncmno; /* new number of updated node on continuum mesh */
//			  fprintf(fp,"oldnode %ld, newnode %ld\n",sfcvno[i],sfcvnonew[i]);
			  ncmno=ncmno+1;
			}
		      else
			{
			  sfcvno[i]=jo2cn[i][0]; /* original number of updated node on continuum mesh */
			  sfcvnonew[i]=sfcvno[i]; /* duplicate */
                          if(node[0]<0)
                            {
                              node[0]=sfcvno[i];
                            }
                          else
                            {
                              node[1]=sfcvno[i];
                            }
			}
		    }

                  for(i=0;i<ncmel;i++)
		    {
		      for(j=0;j<4;j++)
                        {
                          i2cv0to[i][j]=i2cvto[i][j];
                        }
		    }

                  femflag=0;
                  for(i=0;i<3;i++)
                    {
		      if(jo2cn[i][3]==2)
			{
			  YFEMnei(
				  sfcvno[i],sfcvnonew[i],femup,femold,&femflag,
				  i1ptyp,nelem,i1elpr,i2elcon,i2cv0to,
				  i2xbto,i2eltoxb,
				  i2xbjno,i2elto);
                            
			}
		    }

		  if(femflag==1)
		    {
//		      fprintf(fp,"elem%ld edonsf2 femflag %ld femup %ld femold %ld\n",ielem,femflag,femup,femold);
		      i2elcon[femup][femup2old]=femold;
		      i2elcon[femold][femold2up]=femup;
                      ncmno=ncmno-1;
		      //printf("elem%ldi2elcon reset.\n",ielem);
		      continue;
		    }

		  /* update continuum tetrahedron elem topology */
		  for(i=0;i<3;i++)
		    {
		      for(j=0;j<4;j++)
			{
			  if((i==0)&&(i2cvto[femup][j]!=sfcvno[0])&&(i2cvto[femup][j]!=sfcvno[1])&&(i2cvto[femup][j]!=sfcvno[2]))
			    {
			      sfcvnonew[3]=i2cvto[femup][j];
			    }
			  if(i2cvto[femup][j]==sfcvno[i])
			    {
			      i2cvto[femup][j]=sfcvnonew[i];
			    }
			}
		    }

		  /* update continuum tetrahedron elem connected surface triangle elem topology */
		  for(i=0;i<4;i++)
		    {
		      if(i2eltoxb[femup][i]>=0)
			{
			  for(j=0;j<3;j++)
			    {
			      for(k=0;k<3;k++)
				{
				  if(i2xbto[i2eltoxb[femup][i]][j]==sfcvno[k])
				    {
				      i2xbto[i2eltoxb[femup][i]][j]=sfcvnonew[k];
				    }
				}
			    }
			}
		    }

		  /* add 1st new surface triangle element, using three new nodes */
		  for(i=0;i<3;i++)
		    {
		      i2xbto[nsfel][i]=sfcvnonew[i];
                      i2xbtojn[nsfel][i]=jo2cn[i][2];
		    }
//		  fprintf(fp,"add surface[femup]%ld: %ld, %ld, %ld\n",nsfel,i2xbto[nsfel][0],i2xbto[nsfel][1],i2xbto[nsfel][2]);
//		  fprintf(fp,"i2eltoxb[%ld]: [0]%ld [1]%ld [2]%ld [3]%ld\n",femup,i2eltoxb[femup][0],i2eltoxb[femup][1],i2eltoxb[femup][2],i2eltoxb[femup][3]);
		  for(i=0;i<4;i++)
		    {
		      if(i2eltoxb[femup][i]<0)
			{
			  i2eltoxb[femup][i]=nsfel;
//			  fprintf(fp,"i2eltoxb[%ld][%ld]=%ld\n",femup,i,i2eltoxb[femup][i]);
			  break;
			}
		    }

                  i1xbtoel[nsfel]=femup;

		  /* add 2nd new surface triangle element, using three old nodes */
		  i2xbto[nsfel+1][0]=sfcvno[2];
		  i2xbto[nsfel+1][1]=sfcvno[1];
		  i2xbto[nsfel+1][2]=sfcvno[0];

                  i2xbtojn[nsfel+1][0]=jo2cn[2][1];
                  i2xbtojn[nsfel+1][1]=jo2cn[1][1];
                  i2xbtojn[nsfel+1][2]=jo2cn[0][1];

//		  fprintf(fp,"add surface[femold]%ld: %ld, %ld, %ld\n",nsfel+1,i2xbto[nsfel+1][0],i2xbto[nsfel+1][1],i2xbto[nsfel+1][2]);
//		  fprintf(fp,"i2eltoxb[%ld]: [0]%ld [1]%ld [2]%ld [3]%ld\n",femold,i2eltoxb[femold][0],i2eltoxb[femold][1],i2eltoxb[femold][2],i2eltoxb[femold][3]);
		  for(i=0;i<4;i++)
		    {
		      if(i2eltoxb[femold][i]<0)
			{
			  i2eltoxb[femold][i]=nsfel+1;
//			  fprintf(fp,"i2eltoxb[%ld][%ld]=%ld\n",femold,i,i2eltoxb[femold][i]);
			  break;
			}
		    }

                  i1xbtoel[nsfel+1]=femold;

		  /* guarantee tetrahedron volume positive */
		  for(i=0;i<3;i++)
		    {
		      i1elto[i]=sfcvno[i];
		    }
		  i1elto[3]=sfcvnonew[3];
		  for(i=1;i<4;i++)
		    {
		      F0[0][i-1]=d1ncix[i1elto[i]]-d1ncix[i1elto[0]];
		      F0[1][i-1]=d1nciy[i1elto[i]]-d1nciy[i1elto[0]];
		      F0[2][i-1]=d1nciz[i1elto[i]]-d1nciz[i1elto[0]];
		    }
		  voli=F0[0][0]*(F0[1][1]*F0[2][2]-F0[1][2]*F0[2][1])-
		    F0[0][1]*(F0[1][0]*F0[2][2]-F0[1][2]*F0[2][0])+
		    F0[0][2]*(F0[1][0]*F0[2][1]-F0[1][1]*F0[2][0]);
		  if(voli<0)    /*modify by JXiang 14032018  */
		    {
		      for(i=0;i<2;i++)
			{
			  temp=i2xbto[nsfel+i][0];
			  i2xbto[nsfel+i][0]=i2xbto[nsfel+i][2];
			  i2xbto[nsfel+i][2]=temp;

                          temp=i2xbtojn[nsfel+i][0];
                          i2xbtojn[nsfel+i][0]=i2xbtojn[nsfel+i][2];
                          i2xbtojn[nsfel+i][2]=temp;
			}
		    }

		  nsfel=nsfel+2;

		  for(i=0;i<3;i++)
		    {
		      if(jo2cn[i][3]==2)
			{
			  flag=0;
			  for(k=0;k<100;k++) /* add joint elem nodes corresponding to new continuum node */
			    {
			      if(i2xbjno[sfcvnonew[i]][k]==jo2cn[i][2]) /* already in the list */
				{
				  flag=1;
				  break;
				}
			    }
			  if(flag==0)
			    {
			      for(k=0;k<100;k++)
				{
				  if(i2xbjno[sfcvnonew[i]][k]<0)
				    {
				      i2xbjno[sfcvnonew[i]][k]=jo2cn[i][2];
				      break;
				    }
				}
			    }

			  for(k=0;k<100;k++) /* delete joint elem nodes corresponding to old continuum node */
			    {
			      if(i2xbjno[sfcvno[i]][k]==jo2cn[i][2])
				{
				  i2xbjno[sfcvno[i]][k]=-1;
				  break;
				}
			    }
			}
		    }

		  for(i=0;i<3;i++) /* update con fem elem topology connected with updated elem*/
		    {
		      if(jo2cn[i][3]==2)
			{
			  YFEMneiUpdate(
					sfcvno[i],sfcvnonew[i],femup,
					i1ptyp,nelem,i1elpr,i2elcon,i2cvto,
					i2xbto,i2eltoxb,
					i2xbjno,i2elto);
			}
		    }

                  UpdateEdge(node[0],node[1],nelem,i2edge);    

		  i1elpr[ielem]=YIPROPMAX*(-2);
		}
	      else if(edonsf==1)	/* one joint element edge on continuum surface */
		{
		  //printf("edonsf %ld\n",edonsf);

//		  fprintf(fp,"edonsf1_elem %ld\n",ielem);
		  /*
		    printf("jocved[3][3]\n");
		    printf("%ld %ld %ld\n%ld %ld %ld\n%ld %ld %ld\n",
		    jocved[0][0],jocved[0][1],jocved[0][2],
		    jocved[1][0],jocved[1][1],jocved[1][2],
		    jocved[2][0],jocved[2][1],jocved[2][2]);
		    printf("jo2cn[3][4]\n");
		    printf("%ld %ld %ld %ld\n%ld %ld %ld %ld\n%ld %ld %ld %ld\n",
		    jo2cn[0][0],jo2cn[0][1],jo2cn[0][2],jo2cn[0][3],
		    jo2cn[1][0],jo2cn[1][1],jo2cn[1][2],jo2cn[1][3],
		    jo2cn[2][0],jo2cn[2][1],jo2cn[2][2],jo2cn[2][3]);

		    for(i=0;i<nelem;i++)
		    {
		    if(i2edge[i][0]>=0)
		    {
		    printf("i2edge[%ld][] %ld %ld %ld\n",i,i2edge[i][0],i2edge[i][1],i2edge[i][2]);
		    }
		    }
		  */		  
		  /* compare unseparated edge - one end inside, one end on surface */
                  flag=0;
		  for(i=0;i<3;i++)
		    {
		      if(jocved[i][2]==(-1))
			{
			  //printf("i=%ld, ielem=%ld\n",i,ielem);
			  //flag=0;

			  for(k=0;k<3;k++)   //jocved[i][0] - surface node, jocved[i][1] - inside node
			    {
			      if((jo2cn[k][3]==1)&&(jo2cn[k][0]==jocved[i][1]))
				{
				  temp=jocved[i][0];
				  jocved[i][0]=jocved[i][1];
				  jocved[i][1]=temp;
				  break;
				}
			    }
				      
			  for(k=0;k<nelem;k++)
			    {	      
                              //if(ielem==172111)
			      //{
			      //printf("k %ld nelem %ld ",k,nelem);
			      //printf("i2edge[%ld][2] %ld ielem %ld\n",k,i2edge[k][2],ielem);
			      /*
				if(k==81)
				{
				for(m=0;m<nelem;m++)
				{
				printf("i2edge[%ld][2]=%ld ",m,i2edge[m][2]);
				}
				printf("i2edge[k][2] %ld ielem %ld\n",i2edge[k][2],ielem);
				}
			      */
			      /*
				if(k==81)
				{
				printf("i2edge[k][2] %ld\nielem %ld\ni1elpr[i2edge[k][2]]+YIPROPMAX %ld\nnprop %ld\n",
				i2edge[k][2],ielem,(i1elpr[i2edge[k][2]]+YIPROPMAX),nprop);
				printf("i1ptyp[i1elpr[i2edge[k][2]]+YIPROPMAX] %ld\n",i1ptyp[i1elpr[i2edge[k][2]]+YIPROPMAX]);
				printf("i2edge[k][0] %ld\ni2edge[k][1] %ld\n",i2edge[k][0],i2edge[k][1]);
				printf("i %ld\njocved[i][0] %ld\njocved[i][1] %ld\n",i,jocved[i][0],jocved[i][1]);
				}
			      */
			      //}
			     
			      /* two unseparated edges overlap - separate the overlap node on surface */
			      if((i2edge[k][2]>=0)&&(i2edge[k][2]<nelem)&&(i2edge[k][2]!=ielem))
				{
				  //if(ielem==172111)
				  //{
				  //printf("i1elpr[i2edge[k][2]]+YIPROPMAX %ld nprop %ld\n",(i1elpr[i2edge[k][2]]+YIPROPMAX),nprop);
				  //}

				  if(((i1elpr[i2edge[k][2]]+YIPROPMAX)>=0)&&((i1elpr[i2edge[k][2]]+YIPROPMAX)<nprop))
				    {
				      //if(ielem==172111)
				      //{
				      //printf("i1ptyp[i1elpr[i2edge[k][2]]+YIPROPMAX] %ld i %ld\n",i1ptyp[i1elpr[i2edge[k][2]]+YIPROPMAX],i);
				      //}

				      if((i1ptyp[i1elpr[i2edge[k][2]]+YIPROPMAX]==(YTE3TET4JOINT))&&
					 (i2edge[k][0]==jocved[i][0])&&(i2edge[k][1]==jocved[i][1]))
					{
//					  fprintf(fp,"ielem=%ld, broken edge-i2edge[%ld][2] %ld\n",ielem,k,i2edge[k][2]);

					  femup=i2eljp[ielem][1];
					  femold=i2eljp[ielem][0];

					  //printf("femup=%ld,femold=%ld\n",femup,femold);

					  for(m=0;m<4;m++) // set new boundary label -1 to surface
					    {
					      if(i2elcon[femup][m]==femold)
						{
						  i2elcon[femup][m]=-1;
						  femup2old=m;
						  break;
						}
					    }
					  for(m=0;m<4;m++)
					    {
					      if(i2elcon[femold][m]==femup)
						{
						  i2elcon[femold][m]=-1;
						  femold2up=m;
						  break;
						}
					    }

					  //fprintf(fp,"boundary label updated\n");
					  //printf("i2edge[k][0] %ld\n",i2edge[k][0]);
					  //printf("i2edge[k][1] %ld\n",i2edge[k][1]);

					  //printf("i2eljp[i2edge[k][2]][0] %ld\n",i2eljp[i2edge[k][2]][0]);
					  //printf("i2eljp[i2edge[k][2]][1] %ld\n",i2eljp[i2edge[k][2]][1]);
					  if(femup==i2eljp[i2edge[k][2]][0])
					    {
					      femuppair=i2eljp[i2edge[k][2]][0];
					      femoldpair=i2eljp[i2edge[k][2]][1];
					    }
					  else if(femup==i2eljp[i2edge[k][2]][1])
					    {
					      femuppair=i2eljp[i2edge[k][2]][1];
					      femoldpair=i2eljp[i2edge[k][2]][0];
					    }
					  else
					    {
					      femuppair=-1;
					      femnei=-1;
					      for(n=0;n<4;n++)
						{
						  if(i2elcon[femup][n]>=0)
						    {
						      for(p=0;p<4;p++)
							{
							  if(i2cvto[i2elcon[femup][n]][p]==i2edge[k][0])
							    {
							      for(q=0;q<4;q++)
								{
								  if((q!=p)&&(i2cvto[i2elcon[femup][n]][q]==i2edge[k][1]))
								    {
								      femnei=i2elcon[femup][n];
								      break;
								    }
								}
							      if(femnei>=0)
								{
								  break;
								}
							    }
							}
						      if(femnei>=0)
							{
							  break;
							}
						    }
						}
					      fempre=femup;
					      //printf("init: fempre %ld, femnei %ld\n",fempre,femnei);

					      for(m=0;m<100;m++)
						{
						  if(femnei>=0)
						    {
						      if(femnei==i2eljp[i2edge[k][2]][0])
							{
							  femuppair=i2eljp[i2edge[k][2]][0];
							  femoldpair=i2eljp[i2edge[k][2]][1];
							  break;
							}
						      else if(femnei==i2eljp[i2edge[k][2]][1])
							{
							  femuppair=i2eljp[i2edge[k][2]][1];
							  femoldpair=i2eljp[i2edge[k][2]][0];
							  break;
							}
						      else
							{
							  neiflag=0;
							  for(n=0;n<4;n++)
							    {
							      //printf("i2elcon[%ld][%ld]=%ld\n",femnei,n,i2elcon[femnei][n]);
							      if((i2elcon[femnei][n]>=0)&&(i2elcon[femnei][n]!=fempre))
								{
								  for(p=0;p<4;p++)
								    {
								      if(i2cvto[i2elcon[femnei][n]][p]==i2edge[k][0])
									{
									  for(q=0;q<4;q++)
									    {
									      if((q!=p)&&(i2cvto[i2elcon[femnei][n]][q]==i2edge[k][1]))
										{
										  fempre=femnei;
										  femnei=i2elcon[femnei][n];
										  neiflag=1;
										  //printf("loop: fempre %ld, femnei %ld\n",fempre,femnei);
										  break;
										}
									    }
									  if((neiflag==1)&&(femnei==i2elcon[fempre][n]))
									    {
									      break;
									    }
									}
								    }
								  if((neiflag==1)&&(femnei==i2elcon[fempre][n]))
								    {
								      break;
								    }
								}
							    }
							}
						    }
						}
					    }							      
					  //printf("femuppair=%ld,femoldpair=%ld\n",femuppair,femoldpair);

					  if(femuppair<0)
					    {
					      //flag=1;
					      //break;
					      continue;
					    }

					  for(m=0;m<4;m++) // set new boundary label -1 to surface
					    {
					      if(i2elcon[femuppair][m]==femoldpair)
						{
						  i2elcon[femuppair][m]=-1;
						  fempairup2old=m;
						  break;
						}
					    }
					  for(m=0;m<4;m++)
					    {
					      if(i2elcon[femoldpair][m]==femuppair)
						{
						  i2elcon[femoldpair][m]=-1;
						  fempairold2up=m;
						  break;
						}
					    }

					  for(m=0;m<3;m++) /* update overlap node on surface */
					    {
					      pflag[m]=YIPROPMAX;
					      if(jo2cn[m][3]==1)
						{
						  pflag[m]=-1;
						  for(n=0;n<2;n++)
						    {
						      if(jo2cn[m][0]==i2edge[k][n])
							{
							  sfcvno[m]=jo2cn[m][0]; /* original number of updated node on continuum mesh */
							  sfcvnonew[m]=ncmno; /* new number of updated node on continuum mesh */
							  //printf("sfcvnonew[%ld]=%ld\n",m,sfcvnonew[m]);
//							  fprintf(fp,"oldnode %ld, newnode %ld\n",sfcvno[m],sfcvnonew[m]);
							  ncmno=ncmno+1;
							  pflag[m]=m;
							  //printf("pflag[%ld]=%ld\n",m,pflag[m]);
							  break;
							}
						    }
						  if(pflag[m]==(-1))
						    {
						      sfcvno[m]=jo2cn[m][0]; /* original number of updated node on continuum mesh */
						      sfcvnonew[m]=sfcvno[m]; /* duplicate */

						      oldsfn=jo2cn[m][0];
						    }
						}
					      else
						{
						  sfcvno[m]=jo2cn[m][0]; /* original number of updated node on continuum mesh */
						  sfcvnonew[m]=sfcvno[m]; /* duplicate */
						}
					    }
					  /*
					    for(m=0;m<4;m++)
					    {
					    if((i2cvto[femuppair][m]!=i2cvto[femup][0])&&(i2cvto[femuppair][m]!=i2cvto[femup][1])&&(i2cvto[femuppair][m]!=i2cvto[femup][2])&&(i2cvto[femuppair][m]!=i2cvto[femup][3])
					    &&(i2cvto[femuppair][m]!=i2cvto[femold][0])&&(i2cvto[femuppair][m]!=i2cvto[femold][1])&&(i2cvto[femuppair][m]!=i2cvto[femold][2])&&(i2cvto[femuppair][m]!=i2cvto[femold][3]))
					    {
					    oldsfn=i2cvto[femuppair][m];
					    break;
					    }
					    }
					  */

					  for(m=0;m<ncmel;m++)
					    {
					      for(n=0;n<4;n++)
						{
						  i2cv0to[m][n]=i2cvto[m][n];
						}
					    }

					  femflag=0;
					  for(m=0;m<3;m++)
					    {
					      if(m==pflag[m])
						{
						  YFEMnei(
							  sfcvno[m],sfcvnonew[m],femup,femold,&femflag,
							  i1ptyp,nelem,i1elpr,i2elcon,i2cv0to,
							  i2xbto,i2eltoxb,
							  i2xbjno,i2elto);
						}
					    }
					  /*
					  if((femflag==0)&&(femold!=femoldpair))
					    {
					      for(m=0;m<ncmel;m++)
						{
						  for(n=0;n<4;n++)
						    {
						      i2cv0to[m][n]=i2cvto[m][n];
						    }
						}
					      for(m=0;m<3;m++)
						{
						  if(m==pflag[m])
						    {
						      YFEMnei(
							      sfcvno[m],sfcvnonew[m],femup,femoldpair,&femflag,
							      i1ptyp,nelem,i1elpr,i2elcon,i2cv0to,
							      i2xbto,i2eltoxb,
							      i2xbjno,i2elto);
						    }
						}
					    }
					  */

					  if(femflag==0)
					    {
//					      fprintf(fp,"elem%ld edonsf1 femflag %ld\n",ielem,femflag);
					      for(m=0;m<ncmel;m++)
						{
						  for(n=0;n<4;n++)
						    {
						      i2cv0to[m][n]=i2cvto[m][n];
						    }
						}
					      for(m=0;m<3;m++)
						{
						  if(m==pflag[m])
						    {
						      YFEMnei(
							      sfcvno[m],sfcvnonew[m],femuppair,femoldpair,&femflag,
							      i1ptyp,nelem,i1elpr,i2elcon,i2cv0to,
							      i2xbto,i2eltoxb,
							      i2xbjno,i2elto);
						    }
						}
					    }

					  if(femflag==1)
					    {
//					      fprintf(fp,"elem%ld edonsf1 femflag %ld femup %ld femold %ld femuppair %ld femoldpair %ld\n",ielem,femflag,femup,femold,femuppair,femoldpair);
					      i2elcon[femup][femup2old]=femold;
					      i2elcon[femold][femold2up]=femup;
					      i2elcon[femuppair][fempairup2old]=femoldpair;
					      i2elcon[femoldpair][fempairold2up]=femuppair;
                                              ncmno=ncmno-1;
					      //printf("elem%ldi2elcon reset.\n",ielem);
					      //printf("i2elcon[%ld][%ld]=%ld\n",femup,femup2old,i2elcon[femup][femup2old]);
					      //printf("i2elcon[%ld][%ld]=%ld\n",femold,femold2up,i2elcon[femold][femold2up]);
					      //printf("i2elcon[%ld][%ld]=%ld\n",femuppair,fempairup2old,i2elcon[femuppair][fempairup2old]);
					      //printf("i2elcon[%ld][%ld]=%ld\n",femoldpair,fempairold2up,i2elcon[femoldpair][fempairold2up]);
					      //printf("k=%ld\n",k);
					      continue;
					    }

					  //if(ielem==172111)
					  //{
					  //printf("femflag=%ld\n",femflag);
					  //}

					  /* update continuum tetrahedron elem topology */
					  for(m=0;m<3;m++)
					    {
					      for(n=0;n<4;n++)
						{
						  if((m==0)&&(i2cvto[femup][n]!=sfcvno[0])&&(i2cvto[femup][n]!=sfcvno[1])&&(i2cvto[femup][n]!=sfcvno[2]))
						    {
						      sfcvnonew[3]=i2cvto[femup][n];
						    }
						  if(i2cvto[femup][n]==sfcvno[m])
						    {
						      //if(sfcvno[m]==6253)
						      //{
						      //printf("femup%ldnode%ld %ld %ld %ld %ld\n",femup,n,i2cvto[femup][0],i2cvto[femup][1],i2cvto[femup][2],i2cvto[femup][3]);
						      //}

						      i2cvto[femup][n]=sfcvnonew[m];

						      //if(sfcvno[m]==6253)
						      //{
						      //printf("femup%ldnode%ld %ld %ld %ld %ld\n",femup,n,i2cvto[femup][0],i2cvto[femup][1],i2cvto[femup][2],i2cvto[femup][3]);
						      //}

						    }
						}
					    }

					  /* update continuum tetrahedron elem connected surface triangle elem topology */
					  for(m=0;m<4;m++)
					    {
					      if(i2eltoxb[femup][m]>=0)
						{
						  for(n=0;n<3;n++)
						    {
						      for(p=0;p<3;p++)
							{
							  if(i2xbto[i2eltoxb[femup][m]][n]==sfcvno[p])
							    {
							      i2xbto[i2eltoxb[femup][m]][n]=sfcvnonew[p];
							    }
							}
						    }
						}
					    }

					  /* add 1st new surface triangle element, using three new nodes */
					  for(m=0;m<3;m++)
					    {
					      i2xbto[nsfel][m]=sfcvnonew[m];
					      i2xbtojn[nsfel][m]=jo2cn[m][2];
					    }
//					  fprintf(fp,"add surface[femup]%ld: %ld, %ld, %ld\n",nsfel,i2xbto[nsfel][0],i2xbto[nsfel][1],i2xbto[nsfel][2]);
//					  fprintf(fp,"i2eltoxb[%ld]: [0]%ld [1]%ld [2]%ld [3]%ld\n",femup,i2eltoxb[femup][0],i2eltoxb[femup][1],i2eltoxb[femup][2],i2eltoxb[femup][3]);
					  for(m=0;m<4;m++)
					    {
					      if(i2eltoxb[femup][m]<0)
						{
						  i2eltoxb[femup][m]=nsfel;
//						  fprintf(fp,"i2eltoxb[%ld][%ld]=%ld\n",femup,m,i2eltoxb[femup][m]);
						  break;
						}
					    }

					  i1xbtoel[nsfel]=femup;

					  /* add 2nd new surface triangle element, using three old nodes */
					  i2xbto[nsfel+1][0]=sfcvno[2];
					  i2xbto[nsfel+1][1]=sfcvno[1];
					  i2xbto[nsfel+1][2]=sfcvno[0];

					  i2xbtojn[nsfel+1][0]=jo2cn[2][1];
					  i2xbtojn[nsfel+1][1]=jo2cn[1][1];
					  i2xbtojn[nsfel+1][2]=jo2cn[0][1];
					  
//					  fprintf(fp,"add surface[femold]%ld: %ld, %ld, %ld\n",nsfel+1,i2xbto[nsfel+1][0],i2xbto[nsfel+1][1],i2xbto[nsfel+1][2]);
//					  fprintf(fp,"i2eltoxb[%ld]: [0]%ld [1]%ld [2]%ld [3]%ld\n",femold,i2eltoxb[femold][0],i2eltoxb[femold][1],i2eltoxb[femold][2],i2eltoxb[femold][3]);
					  for(m=0;m<4;m++)
					    {
					      if(i2eltoxb[femold][m]<0)
						{
						  i2eltoxb[femold][m]=nsfel+1;
//						  fprintf(fp,"i2eltoxb[%ld][%ld]=%ld\n",femold,m,i2eltoxb[femold][m]);
						  break;
						}
					    }

					  i1xbtoel[nsfel+1]=femold;

					  /* guarantee tetrahedron volume positive */
					  for(m=0;m<3;m++)
					    {
					      i1elto[m]=sfcvno[m];
					    }
					  i1elto[3]=sfcvnonew[3];
					  for(m=1;m<4;m++)
					    {
					      F0[0][m-1]=d1ncix[i1elto[m]]-d1ncix[i1elto[0]];
					      F0[1][m-1]=d1nciy[i1elto[m]]-d1nciy[i1elto[0]];
					      F0[2][m-1]=d1nciz[i1elto[m]]-d1nciz[i1elto[0]];
					    }
					  voli=F0[0][0]*(F0[1][1]*F0[2][2]-F0[1][2]*F0[2][1])-
					    F0[0][1]*(F0[1][0]*F0[2][2]-F0[1][2]*F0[2][0])+
					    F0[0][2]*(F0[1][0]*F0[2][1]-F0[1][1]*F0[2][0]);
					  if(voli<0)    /*modify by JXiang 14032018  */
					    {
					      for(m=0;m<2;m++)
						{
						  temp=i2xbto[nsfel+m][0];
						  i2xbto[nsfel+m][0]=i2xbto[nsfel+m][2];
						  i2xbto[nsfel+m][2]=temp;

						  temp=i2xbtojn[nsfel+m][0];
						  i2xbtojn[nsfel+m][0]=i2xbtojn[nsfel+m][2];
						  i2xbtojn[nsfel+m][2]=temp;
						}
					    }
					  /*
					    if(nsfel==7823)
					    {
					    printf("femuppair =%ld %ld %ld %ld %ld\n",femuppair,i2cvto[femuppair][0],i2cvto[femuppair][1],i2cvto[femuppair][2],i2cvto[femuppair][3]);
					    printf("femoldpair=%ld %ld %ld %ld %ld\n",femoldpair,i2cvto[femoldpair][0],i2cvto[femoldpair][1],i2cvto[femoldpair][2],i2cvto[femoldpair][3]);
					    printf("i2edge[k][0]=%ld,i2edge[k][1]=%ld\n",i2edge[k][0],i2edge[k][0]);
					    }
					  */
			      
					  newsfn=-1;
					  for(m=0;m<4;m++)
					    {
					      for(n=0;n<4;n++)
						{
						  if((i2cvto[femuppair][m]==i2cvto[femoldpair][n])&&
						     (i2cvto[femuppair][m]!=i2edge[k][0])&&
						     (i2cvto[femuppair][m]!=i2edge[k][1]))
						    {
						      newsfn=i2cvto[femuppair][m];
						      break;
						    }
						}
					      if(newsfn>=0)
						{
						  break;
						}
					    }
					  /*
					    if(nsfel==7823)
					    {
					    printf("newsfn=%ld\n",newsfn);
					    //exit(0);
					    }
					  */
					  /* update pair continuum tetrahedron elem topology */
					  for(m=0;m<3;m++)
					    {
					      if(m==pflag[m])
						{
						  for(n=0;n<4;n++)
						    {
						      if(i2cvto[femuppair][n]==sfcvno[m])
							{
							  //if(sfcvno[m]==6253)
							  //{
							  //printf("femuppair%ldnode%ld %ld %ld %ld %ld\n",femuppair,n,i2cvto[femuppair][0],i2cvto[femuppair][1],i2cvto[femuppair][2],i2cvto[femuppair][3]);
							  //}

							  i2cvto[femuppair][n]=sfcvnonew[m];

							  //if(sfcvno[m]==6253)
							  //{
							  //printf("femuppair%ldnode%ld %ld %ld %ld %ld\n",femuppair,n,i2cvto[femuppair][0],i2cvto[femuppair][1],i2cvto[femuppair][2],i2cvto[femuppair][3]);
							  //}

							  break;
							}
						    }
						}
					    }

					  /* update pair continuum tetrahedron elem connected surface triangle elem topology */
					  for(m=0;m<4;m++)
					    {
					      if(i2eltoxb[femuppair][m]>=0)
						{
						  for(n=0;n<3;n++)
						    {
						      for(p=0;p<3;p++)
							{
							  if(i2xbto[i2eltoxb[femuppair][m]][n]==sfcvno[p])
							    {
							      i2xbto[i2eltoxb[femuppair][m]][n]=sfcvnonew[p];
							    }
							}
						    }
						}
					    }

					  /* add 3rd new surface triangle element, using three new nodes */
					  i2xbto[nsfel+2][0]=i2xbto[nsfel][2];
					  i2xbto[nsfel+2][1]=i2xbto[nsfel][1];
					  i2xbto[nsfel+2][2]=i2xbto[nsfel][0];

					  for(m=0;m<3;m++)
					    {
					      if(i2xbto[nsfel+2][m]==oldsfn)
						{
						  i2xbto[nsfel+2][m]=newsfn;
						  break;
						}
					    }
					  
//					  fprintf(fp,"add surface[femuppair]%ld: %ld, %ld, %ld\n",nsfel+2,i2xbto[nsfel+2][0],i2xbto[nsfel+2][1],i2xbto[nsfel+2][2]);
//					  fprintf(fp,"i2eltoxb[%ld]: [0]%ld [1]%ld [2]%ld [3]%ld\n",femuppair,i2eltoxb[femuppair][0],i2eltoxb[femuppair][1],i2eltoxb[femuppair][2],i2eltoxb[femuppair][3]);

					  for(m=0;m<4;m++)
					    {
					      if(i2eltoxb[femuppair][m]<0)
						{
						  i2eltoxb[femuppair][m]=nsfel+2;
//						  fprintf(fp,"i2eltoxb[%ld][%ld]=%ld\n",femuppair,m,i2eltoxb[femuppair][m]);
						  break;
						}
					    }

					  i1xbtoel[nsfel+2]=femuppair;
				  
					  /* add 4th new surface triangle element, using three old nodes */				  
					  i2xbto[nsfel+3][0]=i2xbto[nsfel+1][2];
					  i2xbto[nsfel+3][1]=i2xbto[nsfel+1][1];
					  i2xbto[nsfel+3][2]=i2xbto[nsfel+1][0];

					  for(m=0;m<3;m++)
					    {
					      if(i2xbto[nsfel+3][m]==oldsfn)
						{
						  i2xbto[nsfel+3][m]=newsfn;
						  break;
						}
					    }
//					  fprintf(fp,"add surface[femoldpair]%ld: %ld, %ld, %ld\n",nsfel+3,i2xbto[nsfel+3][0],i2xbto[nsfel+3][1],i2xbto[nsfel+3][2]);
//					  fprintf(fp,"i2eltoxb[%ld]: [0]%ld [1]%ld [2]%ld [3]%ld\n",femoldpair,i2eltoxb[femoldpair][0],i2eltoxb[femoldpair][1],i2eltoxb[femoldpair][2],i2eltoxb[femoldpair][3]);

					  for(m=0;m<4;m++)
					    {
					      if(i2eltoxb[femoldpair][m]<0)
						{
						  i2eltoxb[femoldpair][m]=nsfel+3;
//						  fprintf(fp,"i2eltoxb[%ld][%ld]=%ld\n",femoldpair,m,i2eltoxb[femoldpair][m]);
						  break;
						}
					    }

					  i1xbtoel[nsfel+3]=femoldpair;

					  /* guarantee tetrahedron volume positive */
					  /*
					    for(m=0;m<3;m++)
					    {
					    i1elto[m]=i2xbto[nsfel+3][m];
					    }
					    for(m=0;m<4;m++)
					    {
					    if((i2cvto[femoldpair][m]!=i1elto[0])&&(i2cvto[femoldpair][m]!=i1elto[1])&&(i2cvto[femoldpair][m]!=i1elto[2])&&
					    (i2cvto[femoldpair][m]!=oldsfn)&&(i2cvto[femoldpair][m]!=newsfn))
					    {
					    i1elto[3]=i2cvto[femoldpair][m];
					    break;
					    }
					    }
					    for(m=1;m<4;m++)
					    {
					    F0[0][m-1]=d1ncix[i1elto[m]]-d1ncix[i1elto[0]];
					    F0[1][m-1]=d1nciy[i1elto[m]]-d1nciy[i1elto[0]];
					    F0[2][m-1]=d1nciz[i1elto[m]]-d1nciz[i1elto[0]];
					    }
					    voli=F0[0][0]*(F0[1][1]*F0[2][2]-F0[1][2]*F0[2][1])-
					    F0[0][1]*(F0[1][0]*F0[2][2]-F0[1][2]*F0[2][0])+
					    F0[0][2]*(F0[1][0]*F0[2][1]-F0[1][1]*F0[2][0]);
					    if(voli<0)
					    {
					    for(m=2;m<4;m++)
					    {
					    temp=i2xbto[nsfel+m][0];
					    i2xbto[nsfel+m][0]=i2xbto[nsfel+m][2];
					    i2xbto[nsfel+m][2]=temp;
					    }
					    }
					  */

					  //nsfel=nsfel+4;				  

					  //nsfel=nsfel+2;

					  for(m=0;m<3;m++)
					    {
					      if(m==pflag[m])
						{
						  flag2=0;
						  for(n=0;n<100;n++) /* add joint elem nodes corresponding to new continuum node */
						    {
						      if(i2xbjno[sfcvnonew[m]][n]==jo2cn[m][2]) /* already in the list */
							{
							  flag2=1;
							  break;
							}
						    }
						  if(flag2==0)
						    {
						      for(n=0;n<100;n++)
							{
							  if(i2xbjno[sfcvnonew[m]][n]<0)
							    {
							      i2xbjno[sfcvnonew[m]][n]=jo2cn[m][2];
							      break;
							    }
							}
						    }

						  for(n=0;n<100;n++) /* delete joint elem nodes corresponding to old continuum node */
						    {
						      if(i2xbjno[sfcvno[m]][n]==jo2cn[m][2])
							{
							  i2xbjno[sfcvno[m]][n]=-1;
							  break;
							}
						    }

						  flag3=0;
						  for(p=0;p<100;p++)
						    {
						      for(q=0;q<4;q++)
							{
							  if(i2xbjno[sfcvno[m]][p]==i2elto[femuppair][q])
							    {
							      i2xbjno[sfcvno[m]][p]=-1;

							      flag2=0;
							      for(n=0;n<100;n++) /* add joint elem nodes corresponding to new continuum node */
								{
								  if(i2xbjno[sfcvnonew[m]][n]==i2elto[femuppair][q]) /* already in the list */
								    {
								      flag2=1;
								      break;
								    }
								}
							      if(flag2==0)
								{
								  for(n=0;n<100;n++)
								    {
								      if(i2xbjno[sfcvnonew[m]][n]<0)
									{
									  i2xbjno[sfcvnonew[m]][n]=i2elto[femuppair][q];
									  break;
									}
								    }
								}
							      flag3=1;
							      break;
							    }
							}
						      if(flag3==1)
							{
							  break;
							}
						    }
						}
					    }

					  for(m=0;m<3;m++)
					    {
					      for(n=0;n<100;n++)
						{
						  for(p=0;p<4;p++)
						    {
						      if(i2xbjno[i2xbto[nsfel+2][m]][n]==i2elto[femuppair][p])
							{
							  i2xbtojn[nsfel+2][m]=i2elto[femuppair][p];
							  break;
							}
						    }
						  if(i2xbtojn[nsfel+2][m]>=0)
						    {
						      break;
						    }
						}
					    }

					  for(m=0;m<3;m++)
					    {
					      for(n=0;n<100;n++)
						{
						  for(p=0;p<4;p++)
						    {
						      if(i2xbjno[i2xbto[nsfel+3][m]][n]==i2elto[femoldpair][p])
							{
							  i2xbtojn[nsfel+3][m]=i2elto[femoldpair][p];
							  break;
							}
						    }
						  if(i2xbtojn[nsfel+3][m]>=0)
						    {
						      break;
						    }
						}
					    }

					  nsfel=nsfel+4;

					  for(m=0;m<3;m++) /* update con fem elem topology connected with updated elem*/
					    {
					      if(m==pflag[m])
						{
						  YFEMneiUpdate(
								sfcvno[m],sfcvnonew[m],femup,
								i1ptyp,nelem,i1elpr,i2elcon,i2cvto,
								i2xbto,i2eltoxb,
								i2xbjno,i2elto);
						  YFEMneiUpdate(
								sfcvno[m],sfcvnonew[m],femuppair,
								i1ptyp,nelem,i1elpr,i2elcon,i2cvto,
								i2xbto,i2eltoxb,
								i2xbjno,i2elto);
						}
					    }
					  /*
					    for(m=0;m<3;m++) // update con fem elem topology connected with pair updated elem
					    {
					    if(m==pflag[m])
					    {
					    YFEMneiUpdate(
					    sfcvno[m],sfcvnonew[m],femuppair,
					    i1ptyp,nelem,i1elpr,i2elcon,i2cvto,
					    i2xbto,i2eltoxb,
					    i2xbjno,i2elto);
					    }
					    }
					  */

					  flag=1;
					  /*
					    for(m=0;m<2;m++)
					    {
					    i2edge[k][m]=-1;
					    }
					  */
					  i1elpr[ielem]=YIPROPMAX*(-2);
					  //printf("i2edge[k][2] %ld\n",i2edge[k][2]);
					  i1elpr[i2edge[k][2]]=YIPROPMAX*(-2);
					  /*
					    for(m=0;m<nelem;m++)
					    {
					    if((i2edge[m][2]==i2edge[k][2])&&(m!=k))
					    {
					    for(n=0;n<3;n++)
					    {
					    i2edge[m][n]=-1;
					    }
					    break;
					    }
					    }

					    i2edge[k][2]=-1;
					  */
                                          for(m=0;m<nelem;m++)
                                            {
                                              if((i2edge[m][2]==i2edge[k][2])&&(m!=k))
                                                {
                                                  node[0]=i2edge[m][0];
                                                  node[1]=i2edge[m][1];
                                                  UpdateEdge(node[0],node[1],nelem,i2edge);
                                                }
					    }

                                          for(m=0;m<3;m++)
                                            {
                                              if(jocved[m][2]==(-1))
                                                {
                                                  node[0]=jocved[m][0];
                                                  node[1]=jocved[m][1];
                                                  UpdateEdge(node[0],node[1],nelem,i2edge);
                                                }
                                            }

					  break;
					}
				    }
				}
			    }

                          if(flag==1)
                            {
                              break;
                            }
			  /*
			    if(flag==0)	// unseparated edge not in the list - store this edge
			    {
			    flag2=0;
			    for(k=0;k<nelem;k++)
			    {
			    if((i2edge[k][0]>=0)&&(i2edge[k][1]>=0)&&(i2edge[k][2]>=0)&&(i2edge[k][2]==ielem)&&
			    (((i2edge[k][0]==jocved[i][0])&&(i2edge[k][1]==jocved[i][1]))||
			    ((i2edge[k][0]==jocved[i][1])&&(i2edge[k][1]==jocved[i][0]))))
			    {
			    flag2=1;
			    break;
			    }
			    }

			    if(flag2==0)
			    {
			    for(k=0;k<nelem;k++)
			    {
			    if(i2edge[k][0]<0)
			    {
			    i2edge[k][0]=jocved[i][0];
			    i2edge[k][1]=jocved[i][1];
			    i2edge[k][2]=ielem;
			    //printf("added i2edge[k][2] %ld\n",i2edge[k][2]);
			    break;
			    }
			    }
			    }
			    }
			  */
			}
		    }

		  if(flag==0)   // unseparated edge not in the list - store this edge
		    {
		      flag2=0;
		      for(k=0;k<nelem;k++)
			{
			  if(i2edge[k][2]==ielem)
			    {
			      flag2=1;
			      break;
			    }
			}

		      if(flag2==0)
			{
			  for(m=0;m<3;m++)
			    {
			      if(jocved[m][2]==(-1))
				{
                                  for(k=0;k<nelem;k++)
                                    {
                                      if(i2edge[k][0]<0)
                                        {
				          i2edge[k][0]=jocved[m][0];
				          i2edge[k][1]=jocved[m][1];
				          i2edge[k][2]=ielem;
				          //printf("added i2edge[k][2] %ld\n",i2edge[k][2]);
				          break;
                                        }
                                    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  //}

  (*n0sfel)=nsfel;
  (*n0cmno)=ncmno;
  FREE(i2edge);
  FREE(i2cv0to);
  FREE(i1edsf);
  /*
    for(i=0;i<ncmno;i++)
    {
    fprintf(fp,"node %ld\n",i);
    for(j=0;j<100;j++)
    {
    fprintf(fp,"%ld ",i2xbjno[i][j]);
    }
    fprintf(fp,"\n");
    }
  */
//  fclose(fp);
}
     
/* average properties of joint element nodes to each continuum node */
static void YSFCoorUpdate(
			  INT ncmno, INT **i2xbjno,
			  DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
			  DBL *d1nvcx, DBL *d1nvcy, DBL *d1nvcz,
			  DBL *d1xnccx, DBL *d1xnccy, DBL *d1xnccz,
			  DBL *d1xnvcx, DBL *d1xnvcy, DBL *d1xnvcz, INT *i1ntoC2D)
{
  INT inopo,j,jn;
  INT jono;

  for(inopo=0;inopo<ncmno;inopo++)
    {
      //printf("continuum node %ld\n",inopo);

      d1xnccx[inopo]=R0;
      d1xnccy[inopo]=R0;
      d1xnccz[inopo]=R0;

      d1xnvcx[inopo]=R0;
      d1xnvcy[inopo]=R0;
      d1xnvcz[inopo]=R0;

      jono=0;
      for(j=0;j<100;j++)
	{
	  jn=i2xbjno[inopo][j];

	  if(jn>=0)
	    {


	      d1xnccx[inopo]=d1xnccx[inopo]+d1nccx[jn];
	      d1xnccy[inopo]=d1xnccy[inopo]+d1nccy[jn];
	      d1xnccz[inopo]=d1xnccz[inopo]+d1nccz[jn];

	      d1xnvcx[inopo]=d1xnvcx[inopo]+d1nvcx[jn];
	      d1xnvcy[inopo]=d1xnvcy[inopo]+d1nvcy[jn];
	      d1xnvcz[inopo]=d1xnvcz[inopo]+d1nvcz[jn];
	      jono+=1;
	     i1ntoC2D[jn]=inopo; //CONT to DISC topology - Asiri
	    }
	}
      if(jono>0)
	{

	  d1xnccx[inopo]=d1xnccx[inopo]/jono;
	  d1xnccy[inopo]=d1xnccy[inopo]/jono;
	  d1xnccz[inopo]=d1xnccz[inopo]/jono;
	   
	  d1xnvcx[inopo]=d1xnvcx[inopo]/jono;
	  d1xnvcy[inopo]=d1xnvcy[inopo]/jono;
	  d1xnvcz[inopo]=d1xnvcz[inopo]/jono;
	
	}
    }

}
	      	      
void Ysurface(YDE yde, YDX ydx, YDN ydn)			/* find surface boundary mesh */
{
  YFEMsurface(			/* initiate surface mesh on FEM tetrahedra mesh, before inserting joint elemets */
	      ydx->ncmel,ydx->i2cvto,ydx->i1xbtoel,
	      ydx->i2xbto,&(ydx->nsfel),ydx->i2eltoxb,
	      ydx->nelem_t, ydx->i2elto_t, ydx->nelno);
}

void YSFtoJOINT(YDE yde, YDN ydn, YDX ydx, YDXN ydxn, YDP ydp)		/* link surface nodes to joint nodes */
{

  YlinkCONtoJOINT(		/* link all continuum nodes, including both surface and volume, to joint elem nodes */
		  ydxn->i2xbjno,ydxn->i1xbjo,
		  ydn->nnopo,ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2ncc[2],
		  ydxn->ncmno,ydx->ncmel,ydp->i1ptyp,ydp->i1pejp,yde->i1elpr,yde->i2elto,ydx->i2cvto, ydn->nnopi);

  printf("*YlinkCONtoJOINT* node %ld, element %ld, surface %ld\n",ydxn->ncmno,ydx->ncmel,ydx->nsfel);
/*
  YlinkSFCONtoSFJOINT(		// link surface continuum nodes, to surface joint nodes
  		      ydx->nsfel,ydx->i2xbto,ydxn->i2xbjno,ydx->i2xbtosf,
  		      ydn->nnopo,yde->nelem,ydp->i1ptyp,yde->i1elpr,yde->i2elto
  		      );
*/
  YfindFEMcon(
	      ydp->i1ptyp,ydx->ncmel,yde->i2eljp,yde->i1elpr,ydx->i2elcon);

  printf("*YfindFEMcon* node %ld, element %ld, surface %ld\n",ydxn->ncmno,ydx->ncmel,ydx->nsfel);

  YlinkSFtoJOINT(
		 ydx->nsfel,ydx->i2xbto,ydx->i2xbtojn,ydx->i1xbtoel,
		 ydxn->i2xbjno,
		 yde->i2elto,ydn->d2nci[0],ydn->d2nci[1],ydn->d2nci[2]);
}

void YSFmesh(YDC ydc, YDX ydx, YDXN ydxn, YDN ydn, YDP ydp, YDE yde)			/* update surface mesh information */
{
  INT meshq;
  INT imesh;

  //meshq=1000;
  meshq=ydc->isave;

//  if((ydc->ncstep%meshq)==0)
    //if(ydc->ncstep==10000)
  //  {
      //printf("node %ld, element %ld, surface %ld\n",ydxn->ncmno,ydx->ncmel,ydx->nsfel);

      for(imesh=0;imesh<10;imesh++)
	{
	  YSFMeshUpdate(		/* add new fracture surface and node to the continuum surface mesh */  // -ao changed to neelemi and nnopi
			ydp->nprop,ydp->i1ptyp,yde->nelem,yde->i1elpr,yde->i2elto,
			&(ydx->nsfel),ydx->i2xbto,ydx->i2xbtojn,ydx->i1xbtoel,
			ydxn->i1xbjo,
			ydxn->i2xbjno,yde->i2eljp,&(ydxn->ncmno),
			ydx->i2cvto,ydx->i2eltoxb,ydn->d2nci[0],ydn->d2nci[1],ydn->d2nci[2],
			ydx->i2elcon,ydx->ncmel);
	}

      //printf("node %ld, element %ld, surface %ld\n",ydxn->ncmno,ydx->ncmel,ydx->nsfel);
      YSFCoorUpdate(		/* average properties of joint element nodes to each continuum node */
		    ydxn->ncmno,ydxn->i2xbjno,
		    ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2ncc[2],
		    ydn->d2nvc[0],ydn->d2nvc[1],ydn->d2nvc[2],
		    ydxn->d2xncc[0],ydxn->d2xncc[1],ydxn->d2xncc[2],
		    ydxn->d2xnvc[0],ydxn->d2xnvc[1],ydxn->d2xnvc[2], ydn->i1ntoC2D);

      //printf("node %ld, element %ld, surface %ld\n",ydxn->ncmno,ydx->ncmel,ydx->nsfel);
    }
//}
