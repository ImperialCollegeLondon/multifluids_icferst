/**********************************************************************/
/** Copyright (C) 2008,                                              **/
/** Queen Mary University of London (QMUL) & Imperial College        **/
/** of Science, Technology and Medicine (ICSTM). All rights reserved.**/
/** Implemented for you by Dr Jiansheng Xiang                        **
 
 
 
* This code is part of the Virtual Geoscience Workbench (VGW) developed
* jointly by ICSTM and QMUL through two related parallel projects at 
* ICSTM and QMUL respe
 ctively funded by EPSRC.
*
* This code is provided by copyright holders under the GNU Lesser 
* General Public License (LGPL). It is open source code; you can 
* redistribute it and/or modify it under the terms of the GNU Lesser 
*  
* This code is distributed in the hope that it will be useful, 
* but WITHOUT ANY WARRANTY; without even the implied warranty 
* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
* the GNU Lesser General Public License for more details,
* http://www.gnu.org/licenses/lgpl.txt. 
*  
* You should have received a copy of the GNU Lesser General Public 
* License along with this code; if not, write to 
* Dr Jiansheng Xiang or Dr J.P.Latham 
* j.xiang@imperial.ac.uk or j.p.latham@imperial.ac.uk 
* ******************************************************************* */ 

/*1 File   Yhpd.c */

#include "Yproto.h" 
/****************** PUBLIC ************************************/
void ring(nelem_t,nnopo,nnoder,i2elto, i2elto_t, d2ncc,d2ncN,i1r2s,i1s2r,
         thickness,d2areaw,i2elto_rt,ysp,d1nvol,i2xbtojn,nelem_cv,i2elto_cv,  d1shap)
/*** extrude a ring surface mesh  ***/
INT nelem_t; INT nnopo; INT *nnoder; INT **i2elto; 
DBL **d2ncc; DBL **d2ncN; INT *i1r2s; INT **i2elto_t; 
INT *i1s2r;  DBL thickness;DBL **d2areaw;INT **i2elto_rt;
YSP ysp; DBL *d1nvol; INT **i2xbtojn; INT nelem_cv; INT **i2elto_cv; DBL *d1shap;
{ DBL x[3],y[3],z[3];
	DBL xnc,ync,znc,norm,xc,yc,zc,nx,ny,nz;
	DBL xe,ye,ze;//,fx,fy,fz;
    DBL po11,po12,po13,po21,po22,po23,po31,po32,po33,po41,po42,po43,volume;
	INT i,j,k,n0,n1,n2,n3,index,ii,jj,e0,e1,e2,e3,e4;
    DBL x21,x32,x43,y23,y34,y12,z34,z23,z12;
	DBL xtmp,ytmp,ztmp;
  //  INT *d1nvol_num;
      // DBL *d1nvol;
	INT nelem,nelem_rt;
	//DBL **d2areaw;
	INT Max,Med,Min;
	FILE *fout=FILENULL;
	INT icount,jcount,nnodeq,*i1count,itemp;
    DBL *d1len, voli, F0[3][3];
    
    DBL ap0, ap1, ap2, t1,t2;
    DBL *d1rsize;	


    for (i=0; i<(8*nelem_cv); i++)
    {
        i1r2s[i]=-1;
    }


    
//	node list of surface mesh   */
	nnodeq=0;
	for (i=0; i<nelem_t; i++) 
	{
		for(j=0;j<3;j++)
		{
			index=0;
			n0=i2elto_t[i][j];
			for (k=0; k<nnodeq; k++)
			{
				if(n0==i1r2s[k])
				{
					index=1;
					break;
				}
			}
			if (index==0) 
			{
				i1r2s[nnodeq]=n0;
//				i1r2s[nnodeq+nnoder]=n0;
				i1s2r[n0]=nnodeq;
				nnodeq++;

			}
		}
		
	}
	
//	d2areaw=TalDBL2(4,nnode+1);
//	d2ncN=TalDBL2(3,2*nnode);
//	i2elto=TalINT2(4,ydx->nelem_t*3);
    
    *nnoder=nnodeq;

        i1count=TalzINT1(nnodeq+1);
	 d1rsize=TalzDBL1(nnodeq+1);	  
	  d1len=TalzDBL1(nnodeq+1);
    
	for (i=0; i<nnodeq; i++) 
	{

		d2areaw[0][i]=R0;
		d2areaw[1][i]=R0;
		d2areaw[2][i]=R0;
//		d1areat[i]=R0;
		
	}

	for (i=0; i<nelem_t; i++)
	{		
		n0=i2elto_t[i][0];
		n1=i2elto_t[i][1];
		n2=i2elto_t[i][2];
		x[0]=d2ncc[0][n1]-d2ncc[0][n0];
		y[0]=d2ncc[1][n1]-d2ncc[1][n0];
		z[0]=d2ncc[2][n1]-d2ncc[2][n0];
        
		x[1]=d2ncc[0][n2]-d2ncc[0][n0];
		y[1]=d2ncc[1][n2]-d2ncc[1][n0];
		z[1]=d2ncc[2][n2]-d2ncc[2][n0];

        V3DCro(xnc,ync,znc,x[0],y[0],z[0],x[1],y[1],z[1]);
       
		V3DNor(norm,xnc,ync,znc);
        
      if (d1shap[n0] > 0)
	{
	  t1 = MAXIM(d1shap[n0], d1shap[n1]);
	  t2 = MAXIM(t1, d1shap[n2]);
	  /*	    ap0=d1shap[n0];
	   ap1=d1shap[n1];
	   ap2=d1shap[n2];*/    //use max value for all nodes
	  ap0 = t2;
	  ap1 = ap0;
	  ap2 = ap0;

	}
      else
	{
	  ap0 = thickness;
	  ap1 = thickness;
	  ap2 = thickness;
	}
 
	d1rsize[i1s2r[n0]]=ap0;
        d1rsize[i1s2r[n1]]=ap0;
        d1rsize[i1s2r[n2]]=ap0;

//printf("n0,n1,n2, d1shap[n0],d1shap[n1],d1shap[n2],ap0 %d %d %d  %lf %lf %lf %lf  \n",n0,n1,n2, d1shap[n0],d1shap[n1],d1shap[n2],ap0);
	d2areaw[0][i1s2r[n0]]=d2areaw[0][i1s2r[n0]]+xnc;
	d2areaw[1][i1s2r[n0]]=d2areaw[1][i1s2r[n0]]+ync;
	d2areaw[2][i1s2r[n0]]=d2areaw[2][i1s2r[n0]]+znc;
        
        i1count[i1s2r[n0]]=i1count[i1s2r[n0]]+1;

	d2areaw[0][i1s2r[n1]]=d2areaw[0][i1s2r[n1]]+xnc;
	d2areaw[1][i1s2r[n1]]=d2areaw[1][i1s2r[n1]]+ync;
	d2areaw[2][i1s2r[n1]]=d2areaw[2][i1s2r[n1]]+znc;
        
        i1count[i1s2r[n1]]=i1count[i1s2r[n1]]+1;
        
        d2areaw[0][i1s2r[n2]]=d2areaw[0][i1s2r[n2]]+xnc;
	d2areaw[1][i1s2r[n2]]=d2areaw[1][i1s2r[n2]]+ync;
	d2areaw[2][i1s2r[n2]]=d2areaw[2][i1s2r[n2]]+znc;

        i1count[i1s2r[n2]]=i1count[i1s2r[n2]]+1;


 }
   
	/*normalise wieghted vector  */
	for (i=0; i<nnodeq; i++) 
	{
//        j=i1r2s[i];
        d2areaw[0][i]=d2areaw[0][i]/i1count[i];
        d2areaw[1][i]=d2areaw[1][i]/i1count[i];
        d2areaw[2][i]=d2areaw[2][i]/i1count[i];
	V3DNor(norm,d2areaw[0][i],d2areaw[1][i],d2areaw[2][i]);

        if(norm<=0.2)
        {

	
            d1len[i]=d1rsize[i]; //thickness; //-ao thickness

            d2areaw[0][i]=d2areaw[0][i]+R1;
            d2areaw[1][i]=d2areaw[1][i]+R1;
            d2areaw[2][i]=d2areaw[2][i]+R1;
            V3DNor(norm,d2areaw[0][i],d2areaw[1][i],d2areaw[2][i]);
            
        }
            else d1len[i]=d1rsize[i]/norm;
//printf("i, i1r2s d1len[i] d1rsize[i] %d %d  %lf %lf  \n",i,i1r2s[i], d1len[i], d1rsize[i]);
	}

    CHRw(stdout, "output-4444444 \n");
    printf("In ring, nnodeq: %d\n",nnodeq);


	/*normalise wieghted vector  */
        
       printf("In ring, nnodeq, nelem_t: %d %d\n",nnodeq,nelem_t);
	for (i=0; i<nnodeq; i++) 
	{
		d2ncN[0][i]=d2ncc[0][i1r2s[i]];
		d2ncN[1][i]=d2ncc[1][i1r2s[i]];
		d2ncN[2][i]=d2ncc[2][i1r2s[i]];
	}
    
      printf("\nIn ring2, nnodeq: %d\n",nnodeq);

	for (i=nnodeq; i<(2*nnodeq); i++) 
	{
         //   printf("In ring2, d2ncN: %d  ",i);
		ii=i-nnodeq;
		d2ncN[0][i]=d2ncc[0][i1r2s[ii]]-d1len[ii]*d2areaw[0][ii];
		d2ncN[1][i]=d2ncc[1][i1r2s[ii]]-d1len[ii]*d2areaw[1][ii];
		d2ncN[2][i]=d2ncc[2][i1r2s[ii]]-d1len[ii]*d2areaw[2][ii];

	}
    CHRw(stdout, "output-45555555 \n");
	
	nelem=0;
	printf("nelem_T: %d\n",nelem_t);
/*
	for (i=0; i<nelem_t; i++) 
	{       n0=i2elto_t[i][0];
		n1=i2elto_t[i][1];
		n2=i2elto_t[i][2];
		//printf("loop i: %d   ",i);
		if(n0>n1){
			Max=n0;Min=n1;
			if(n2>Max)
			{
				Max=n2;
				Min=n1;
				Med=n0;
			}
			else {
				Max=n0;
				if (n2>n1)
				{
					Min=n1;
					Med=n2;
				}
				else {
					Min=n2;
					Med=n1;
				}

			}

		}
        
		else {
			Max=n1;Min=n0;
			if(n2>Max)
			{
				Max=n2;
				Min=n0;
				Med=n1;
			}
			else {
				Max=n1;
				if (n2>n0) 
				{
					Min=n0;
					Med=n2;
				}
				else {
					Min=n2;
					Med=n0;
				}
				
			}			
		}

		
		i2elto[nelem][0]=i1s2r[n0];
		i2elto[nelem][1]=i1s2r[n2];
		i2elto[nelem][2]=i1s2r[n1];
		i2elto[nelem][3]=i1s2r[Min]+nnodeq;
		
         
        
		nelem++;
		
		
		i2elto[nelem][0]=i1s2r[n0]+nnodeq;
		i2elto[nelem][1]=i1s2r[n1]+nnodeq;
		i2elto[nelem][2]=i1s2r[n2]+nnodeq;
		i2elto[nelem][3]=i1s2r[Max];
		
         i2elto_rt[i+nelem_t][0]=i1s2r[n0]+nnodeq;
         i2elto_rt[i+nelem_t][1]=i1s2r[n1]+nnodeq;
         i2elto_rt[i+nelem_t][2]=i1s2r[n2]+nnodeq;

		nelem++;
		
		
		i2elto[nelem][0]=i1s2r[Med];
		i2elto[nelem][2]=i1s2r[Med]+nnodeq;
		
		if (((n1>n0)&&(n2>n1))||
			((n1>n0)&&(n0>n2))||
			((n2>n1)&&(n0>n2)))
		{
				i2elto[nelem][1]=i1s2r[Max];
				i2elto[nelem][3]=i1s2r[Min]+nnodeq;
         //    i2elto_rt[i+nelem_t][2]=i1s2r[Min]+nnodeq;
		}
		else
		{
			i2elto[nelem][1]=i1s2r[Min]+nnodeq;
			i2elto[nelem][3]=i1s2r[Max];

		}
		
		nelem++;		
		
	}
*/
    
    for (i=0; i<nelem_t; i++)
    {
        n0=i1s2r[i2elto_t[i][0]];
        n1=i1s2r[i2elto_t[i][1]];
        n2=i1s2r[i2elto_t[i][2]];

        e0=i1s2r[i2elto_t[i][0]]+nnodeq;
        e1=i1s2r[i2elto_t[i][1]]+nnodeq;
        e2=i1s2r[i2elto_t[i][2]]+nnodeq;
        
       
        i2elto[nelem][0]=n1;
        i2elto[nelem][1]=n0;
        i2elto[nelem][2]=n2;
        i2elto[nelem][3]=e0;
        
        nelem++;
        
        i2elto[nelem][0]=e1;
        i2elto[nelem][1]=e2;
        i2elto[nelem][2]=e0;
        i2elto[nelem][3]=n2;
        
        nelem++;
        
        i2elto[nelem][0]=e0;
        i2elto[nelem][1]=n2;
        i2elto[nelem][2]=n1;
        i2elto[nelem][3]=e1;
        
        nelem++;

        i2elto_rt[i][0]=n2;
        i2elto_rt[i][1]=n1;
        i2elto_rt[i][2]=n0;
        
        i2elto_rt[i+nelem_t][0]=e0;
        i2elto_rt[i+nelem_t][1]=e1;
        i2elto_rt[i+nelem_t][2]=e2;

    }
    


   CHRw(stdout, "output-5-1 \n");
  
    for(i=0;i<nelem_t;i++)
    {
        // points in a element
        n0=i2xbtojn[i][0];
        n1=i2xbtojn[i][1];
        n2=i2xbtojn[i][2];

        for(j=0;j<3;j++)
        {
        ii=i*3+j;
        e1=i2elto[ii][0];e2=i2elto[ii][1];e3=i2elto[ii][2];e4=i2elto[ii][3];

        

        x21=d2ncN[0][e2]-d2ncN[0][e1];y23=d2ncN[1][e2]-d2ncN[1][e3];z12=d2ncN[2][e1]-d2ncN[2][e2];

        x32=d2ncN[0][e3]-d2ncN[0][e2];y34=d2ncN[1][e3]-d2ncN[1][e4];z34=d2ncN[2][e3]-d2ncN[2][e4];

        x43=d2ncN[0][e4]-d2ncN[0][e3];y12=d2ncN[1][e1]-d2ncN[1][e2];z23=d2ncN[2][e2]-d2ncN[2][e3];

        

        volume=(x21*y23*z34+x32*y34*z12+x43*y12*z23-x21*y34*z23-x32*y12*z34-x43*y23*z12)/6.0;
            


        if(volume<0){
            
            for(jj=1;jj<4;jj++)
            {
                F0[0][jj-1]=d2ncN[0][i2elto[ii][jj]]-d2ncN[0][i2elto[ii][0]]; /* init. base */
                F0[1][jj-1]=d2ncN[1][i2elto[ii][jj]]-d2ncN[1][i2elto[ii][0]];
                F0[2][jj-1]=d2ncN[2][i2elto[ii][jj]]-d2ncN[2][i2elto[ii][0]];
            }
            voli=F0[0][0]*(F0[1][1]*F0[2][2]-F0[1][2]*F0[2][1])-
            F0[0][1]*(F0[1][0]*F0[2][2]-F0[1][2]*F0[2][0])+
            F0[0][2]*(F0[1][0]*F0[2][1]-F0[1][1]*F0[2][0]);

            volume=-1.0*volume;
            
            CHRw(stdout,"volume negative value  "); INTw(stdout,i,10); CHRwsp(stdout);
            INTw(stdout,j,10); CHRwsp(stdout); DBLw(stdout,voli,10); CHRwcr(stdout);
            itemp=i2elto[ii][0];
            i2elto[ii][0]=i2elto[ii][1];
            i2elto[ii][1]=itemp;
            
        }
            
        d1nvol[n0]=d1nvol[n0]+volume/R3;

       // d1nvol_num[e1]=d1nvol_num[e1]+1;

        d1nvol[n1]=d1nvol[n1]+volume/R3;

       // d1nvol_num[e2]=d1nvol_num[e2]+1;

        d1nvol[n2]=d1nvol[n2]+volume/R3;

      //  d1nvol_num[e3]=d1nvol_num[e3]+1;

//        d1nvol[e4]=d1nvol[e4]+volume/R4;

      //  d1nvol_num[e4]=d1nvol_num[e4]+1;

//      printf("j, n0, n1,n2, volume, d1nolv %d %d %d %d %lf %lf %lf %lf \n",j,n0,n1,n2,volume,d1nvol[n0],d1nvol[n1],d1nvol[n2]);
 //      DBLw(stdout,volume,10); CHRwcr(stdout);
        }
        
        
    }
    
  //  for(i=0;i<nelem_t;i++)
   // {
      //  d1nvol[i]=d1nvol[i];///d1nvol_num[i];
  //  }
    
	
  //  CHRw(stdout, "output-6666666 \n");


	fout=fopen("./ring.msh","w");
 //   CHRw(stdout, "06072017_6 \n");

/*for debug purpose only, this part export tetrahedron mesh for the ring */
	

	
	CHRw(fout,"MESH    dimension 3 ElemType Tetrahedra  nnode ");
    INTw(fout,4, 4);
      //  CHRw(stdout, "06072017_7 \n");
	CHRwsp(fout);
     //   CHRw(stdout, "06072017_8 \n");
	CHRwcr(fout);
       // CHRw(stdout, "06072017_9 \n");
	CHRw(fout,"Coordinates ");
       // CHRw(stdout, "06072017_10 \n");
	CHRwcr(fout);
    
    CHRw(stdout, "output-7777777 \n");
	
	for(icount=0;icount<(2*nnodeq);icount++)
	{
		INTw(fout,icount+1,5);
		CHRwsp(fout);
		DBLw(fout,d2ncN[0][icount],19);
		CHRwsp(fout);
		DBLw(fout,d2ncN[1][icount],19);
		CHRwsp(fout);
		DBLw(fout,d2ncN[2][icount],19);
		CHRwsp(fout);
		CHRwcr(fout);
	}
	
	CHRw(fout,"end coordinates ");
	CHRwcr(fout);
	
	CHRw(fout,"Elements ");
	CHRwcr(fout);
	
	for(icount=0;icount<nelem;icount++)
	{
		INTw(fout,icount+1,5);
		CHRwsp(fout);
		for(jcount=0;jcount<4;jcount++)
		{
			INTw(fout,i2elto[icount][jcount]+1,10);
			CHRwsp(fout);
		}

		CHRwcr(fout);
	}
	CHRw(fout,"end elements ");
	CHRwcr(fout);  
	
	fclose(fout);
    
    FREE(d1len);
	FREE(d1rsize);
    FREE(i1count);
	
	//FREE(d2areaw);
 /*
    FREE(d1nz);
    FREE(d1ny);
    FREE(d1nx);
	*/
  }


void Yring(ydx,ydxn, ydn,ysp)
/*** extrude a ring surface mesh  ***/
YDX ydx; YDXN ydxn; YDN ydn; YSP ysp;
{
    TzDBL1(ydn->d1nvol, ydn->mnopo);
    ring(ydx->nsfel,ydxn->ncmno,&ydx->nnode,ydx->i2elto_r,ydx->i2xbto,
		 ydxn->d2xncc,ydx->d2ncc,ydx->i1r2s,ydx->i1s2r,ydx->dthick,
         ydx->d2areaw,ydx->i2elto_rt,ysp,ydn->d1nvol,ydx->i2xbtojn,ydx->ncmel,ydx->i2cvto,  ydxn->d1shap);
	
}

