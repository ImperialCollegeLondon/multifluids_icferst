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

/*1 File   Y3Did.c */

#include "Yproto.h"

void V3DRot(afx,afy,afz,bfx,bfy,bfz,u,v,w,cosin,sine)
     DBL *afx; DBL *afy; DBL *afz; DBL bfx; DBL bfy; 
     DBL bfz; DBL u; DBL v; DBL w; DBL cosin; DBL sine;
{
  DBL a0,a1,a2,a3,a4,a5,a6;
  a1=u*bfx;
  a2=v*bfy;
  a3=w*bfz;
  a4=u*u;
  a5=v*v;
  a6=w*w;
  a0=a1+a2+a3;
  (*afx)=u*a0+(bfx*(a5+a6)-u*(a2+a3))*cosin+(-w*bfy+v*bfz)*sine;
  (*afy)=v*a0+(bfy*(a4+a6)-v*(a1+a3))*cosin+(w*bfx-u*bfz)*sine;
  (*afz)=w*a0+(bfz*(a5+a4)-w*(a1+a2))*cosin+(-v*bfx+u*bfy)*sine;

}

static void Yid3TET2TET(/* tetrahedra to tetrahedra */
						nelem,
            dcstec,diezon, iprop, jprop,
            d1iesl,d1nccx,d1nccy,d1nccz,
			d1nfcx,d1nfcy,d1nfcz,
            d1nvcx,d1nvcy,d1nvcz,
			d1nvx,d1nvy,d1nvz,
			d1deltat1,d1deltat2,d1deltan,
			d1t1vx,d1t1vy,d1t1vz,
			d1t2vx,d1t2vy,d1t2vz,
			d1pepe,i1elcf,i1elpr,d1pefr,
			d1pesf,d1pevf,d1pepf,d1pepsf,
						i1icff,i1iecn,i1iect,i2elto,
						i1ptyp,ncstep,i1fcstep
            )
  INT    nelem;
  DBL   dcstec; DBL   diezon; INT    iprop; INT    jprop;
  DBL  *d1iesl; DBL  *d1nccx; DBL  *d1nccy; DBL  *d1nccz; 
  DBL  *d1nfcx; DBL  *d1nfcy; DBL  *d1nfcz;
  DBL  *d1nvcx; DBL  *d1nvcy; DBL  *d1nvcz;
  DBL   *d1nvx; DBL   *d1nvy; DBL   *d1nvz;
  DBL  *d1deltat1;  DBL   *d1deltat2; DBL   *d1deltan; 
  DBL  *d1t1vx; DBL  *d1t1vy; DBL  *d1t1vz; 
  DBL  *d1t2vx; DBL  *d1t2vy; DBL  *d1t2vz; DBL *d1pefr; 
  DBL  *d1pesf; DBL  *d1pevf; DBL  *d1pepf; DBL  *d1pepsf; 
  DBL  *d1pepe; INT  *i1elcf; INT  *i1elpr;
  INT  *i1icff; INT  *i1iecn; INT  *i1iect; INT **i2elto;
  INT *i1ptyp; INT ncstep; INT *i1fcstep;
{
  DBL tmp,theigh,penetr,peneto,penetu,penetv,penalty;
  DBL force,forco,uforc,vforc,factor,fact0,facti,fact1;
  DBL xorig,yorig,zorig,xe[2],ye[2],ze[2],dct[4];
  DBL dsc[6][3],dcs[3][6],us[6],vs[6],ub[10],vb[10],anb[10],penetb[10];
  DBL xt[4],yt[4],zt[4],ut[4],vt[4],ft[4],xcent,ycent,zcent,xnt,ynt,znt; 
  DBL xc[4],yc[4],zc[4],uc[4],vc[4],fc[4],xcenc,ycenc,zcenc,xnc,ync,znc;
  DBL zone2,dmin2,factor1,ru,rv,rw,ktss;
  DBL fcx[4],fcy[4],fcz[4],ftx[4],fty[4],ftz[4],fnx,fny,fnz;
  DBL fnc[4],fnt[4],fnctot,fnttot,fn_normal,vxc,vyc,vzc,vxt,vyt,vzt;
  DBL nx,ny,nz,pnx,pny,pnz,cosin,t1x,t1y,t1z,t2x,t2y,t2z;
  DBL sine,vs1,vs2,deltat1,deltat2,fs,fs1,fs2,mu,fx,fy,fz,sfactor;
  INT kprop,icontact,ielem,jelem,icoup,jcoup,fnonzero,ihycon,index;
 
  INT i,j,k,inext,jnext,itars,icons;
  INT nspoin,ninerc,niners,nbpoin,ipt[4],ipc[4],innerc[3],inners[6];
  INT ij[10] = { -1, -1, -1, -1, 0, 1, 2, 0, 1, 2 };
  INT ik[10] = { -1, -1, -1, -1, 1, 2, 0, 3, 3, 3 };

  INT iptn[4],ipcn[4],iptn1[4],ipcn1[4],m,ii1;
  DBL mus,s_vel,mud,deltan,mudd,h_cr,vn,muss;
  INT *i1celto, *i1telto;
  INT p1[4] = { 1, 3, 3, 1 };
  INT p2[4] = { 2, 2, 0, 0 };
  INT p3[4] = { 3, 0, 1, 2 };
  DBL factor2;
  DBL peratio,pepe_n;

  zone2=(R4*diezon*diezon);
  penalty=MINIM(d1pepe[iprop],d1pepe[jprop])/100.0;

  if(d1pefr[iprop]<=d1pefr[jprop])
    {
      mud=d1pefr[iprop];
      mus=d1pesf[iprop];
      muss=d1pepf[iprop];
      s_vel=d1pevf[iprop];
      h_cr=d1pepsf[iprop]/penalty;
    }
  else
    {
      mud=d1pefr[jprop];
      mus=d1pesf[jprop];
      muss=d1pepf[jprop];
      s_vel=d1pevf[jprop];
      h_cr=d1pepsf[jprop]/penalty;
    }

  ktss=2.0/7.0*penalty;	/* modified by LGuo */
  ii1=0;
  for(ielem=0;ielem<nelem;ielem++)
    {
      if(i1elpr[ielem]==iprop)
	{
	  kprop=jprop;
	}
      else if(i1elpr[ielem]==jprop)
	{
	  kprop=iprop;
	}
      else
	{
	  kprop=-1;
	}
      if(kprop>=0)
	{
	  icoup=i1elcf[ielem];
	  jcoup=-1;
	  while(icoup>=0)
	    {
	      ihycon=1;
	      if(i1iect[icoup]<0)
		{
		  ihycon=0;
		  i1iect[icoup]=-1-i1iect[icoup];
		  d1deltat1[icoup]=R0;
		  d1deltat2[icoup]=R0;
		}

	      if((i1fcstep[icoup]>=0)&&(ncstep>=i1fcstep[icoup])&&(ncstep<=(i1fcstep[icoup]+1000)))
		{
		  peratio=1.0*(ncstep-i1fcstep[icoup])/1000.0;
                  peratio=1.0;

		}
	      else
		{
		  peratio=1.0;
		}
	      
	      pepe_n=peratio*penalty;
	    
	      jelem=i1iect[icoup];
	      icontact=-2;
 
	      index=Ycd3TETintersection(d1nccx,d1nccy,d1nccz,ielem,jelem,i2elto,R0);
	      if(index==2)
		{
		  index=Ycd3TETintersection(d1nccx,d1nccy,d1nccz,jelem,ielem,i2elto,R0);
		}

	      fnonzero=0;

	      if((i1elpr[jelem]==kprop)&&(index==1))
		{
		  icontact=-1;
		  dmin2=R0;
		  i1telto = i2elto[jelem]; /*Z */
		  i1celto = i2elto[ielem];

		  /*set centres of contactor and target object */
		  xcent=R0;
		  ycent=R0;
		  zcent=R0;
		  xcenc=R0;
		  ycenc=R0;
		  zcenc=R0;
		  for(i=0;i<4;i++)
		    {
		      xcenc=xcenc+RP25*d1nccx[i1celto[i]]; 
		      ycenc=ycenc+RP25*d1nccy[i1celto[i]];
		      zcenc=zcenc+RP25*d1nccz[i1celto[i]];
		      xcent=xcent+RP25*d1nccx[i1telto[i]];
		      ycent=ycent+RP25*d1nccy[i1telto[i]];
		      zcent=zcent+RP25*d1nccz[i1telto[i]];
		    }

		  for(i=0;i<4;i++)
		    {
		      fcx[i]=d1nfcx[i1celto[i]];
		      fcy[i]=d1nfcy[i1celto[i]];
		      fcz[i]=d1nfcz[i1celto[i]];
		      ftx[i]=d1nfcx[i1telto[i]];
		      fty[i]=d1nfcy[i1telto[i]];
		      ftz[i]=d1nfcz[i1telto[i]];
		    }

		  if(i1ptyp[i1elpr[jelem]]==YTE3TET4ELS)
		    {
		      factor2=R1;
		    }
		  else
		    {
		      factor2=R2/R5;
		    }
		  if(i1ptyp[i1elpr[ielem]]==YTE3TET4ELS)
		    {
		      factor1=R1;
		    }
		  else
		    {
		      factor1=R2/R5;
		    }

		  /*********************************************************/     
		  /*                loop over target surfaces              */
		  /*********************************************************/
		  for(itars=0;itars<4;itars++)
		    { 
		      iptn1[0] = itars;
		      ipt[0] = i1telto[iptn1[0]];

		      iptn1[1] = p1[itars];
		      ipt[1] = i1telto[iptn1[1]];

		      iptn1[2] = p2[itars];
		      ipt[2] = i1telto[iptn1[2]];

		      iptn1[3] = p3[itars];
		      ipt[3] = i1telto[iptn1[3]];

		      /*********************************************************/     
		      /*                loop over contactor surfaces              */
		      /*********************************************************/    
		      for(icons=0;icons<4;icons++)
			{ 
			  ipcn1[0] = icons;
			  ipc[0] = i1celto[ipcn1[0]];

			  ipcn1[1] = p1[icons];
			  ipc[1] = i1celto[ipcn1[1]];

			  ipcn1[2] = p2[icons];
			  ipc[2] = i1celto[ipcn1[2]];

			  ipcn1[3] = p3[icons];
			  ipc[3] = i1celto[ipcn1[3]];

			  for(m=0;m<4;m++)
			    {
			      iptn[iptn1[m]]=m;
			      ipcn[ipcn1[m]]=m;
			    }

			  /* set nodal coordinates */
			  for(i=0;i<3;i++)
			    {
			      xt[i]=d1nccx[ipt[i]];
			      yt[i]=d1nccy[ipt[i]];
			      zt[i]=d1nccz[ipt[i]];
			      xc[i]=d1nccx[ipc[i]];
			      yc[i]=d1nccy[ipc[i]];
			      zc[i]=d1nccz[ipc[i]];
			    }

			  xt[3]=xcent;
			  yt[3]=ycent;
			  zt[3]=zcent;
			  xc[3]=xcenc;
			  yc[3]=ycenc;
			  zc[3]=zcenc;
			  xorig=xc[0];
			  yorig=yc[0];
			  zorig=zc[0];

			  for(i=0;i<4;i++)
			    {
			      xt[i]=xt[i]-xorig;
			      yt[i]=yt[i]-yorig;
			      zt[i]=zt[i]-zorig;
			      xc[i]=xc[i]-xorig;
			      yc[i]=yc[i]-yorig;
			      zc[i]=zc[i]-zorig; 
			    } 
			  /* contactor normal, e-base and target points in e-base */
			  V3DCro(xnc,ync,znc,xc[1],yc[1],zc[1],xc[2],yc[2],zc[2]);
			  V3DNor(xe[0],xnc,ync,znc);
			  xe[0]=xc[1];
			  ye[0]=yc[1];
			  ze[0]=zc[1];
			  V3DNor(xe[1],xe[0],ye[0],ze[0]); 
			  V3DCro(xe[1],ye[1],ze[1],xnc,ync,znc,xe[0],ye[0],ze[0]);
			  for(i=0;i<4;i++)
			    {
			      V3DDot(dct[i],xnc,ync,znc,xt[i],yt[i],zt[i]);
			      V3DDot(ut[i],xt[i],yt[i],zt[i],xe[0],ye[0],ze[0]);
			      V3DDot(vt[i],xt[i],yt[i],zt[i],xe[1],ye[1],ze[1]);
			    }
			  if((dct[0]<=EPSILON)&&(dct[1]<=EPSILON)&&
			     (dct[2]<=EPSILON))
			    {
			      continue;  //modified by JXiang
			    }

			  /* u,v coordinates of S-points and C-points   */
			  nspoin=0;
			  for(i=0;i<3;i++)
			    {
			      for(j=0;j<2;j++)
				{
				  inext=i+1;
				  if(inext>2)
				    {
				      inext=0;
				    }
				  if(j==0)
				    {
				      inext=3;
				    }
				  if(((dct[i]>EPSILON)&&(dct[inext]<NEPSILON))||   
				     ((dct[i]<NEPSILON)&&(dct[inext]>EPSILON)))
				    ///Modified by JXiang
				    {
				      factor=ABS(dct[i]-dct[inext]);          
				      if(factor>EPSILON)
					{
					  factor=ABS(dct[i]/factor);
					  us[nspoin]=factor*ut[inext]+(R1-factor)*ut[i];
					  vs[nspoin]=factor*vt[inext]+(R1-factor)*vt[i];
					  inners[nspoin]=0;
					  nspoin=nspoin+1;
					}
				    }
				}
			    }
			  if((nspoin<3)||(nspoin>4))continue;
			  /* check odering of S-points */
			  if(((us[1]-us[0])*(vs[2]-vs[0])-
			      (vs[1]-vs[0])*(us[2]-us[0]))<R0)
			    {
			      i=0;
			      j=nspoin-1;
			      while(i<j)
				{
				  k=inners[i];
				  inners[i]=inners[j];
				  inners[j]=k;
				  tmp=us[i];
				  us[i]=us[j];
				  us[j]=tmp;
				  tmp=vs[i];
				  vs[i]=vs[j];
				  vs[j]=tmp;
				  i++; j--;
				}
			    }
			  for(i=0;i<3;i++)
			    {
			      V3DDot(uc[i],xc[i],yc[i],zc[i],xe[0],ye[0],ze[0]);
			      V3DDot(vc[i],xc[i],yc[i],zc[i],xe[1],ye[1],ze[1]);
			      innerc[i]=0;
			    }
			  /* distances of C-points from S edges */
			  niners=0;
			  ninerc=0;
			  for(i=0;i<nspoin;i++)
			    {
			      inext=i+1;
			      if(inext>=nspoin)
				{
				  inext=0;
				}
			      for(j=0;j<3;j++) 
				{
				  jnext=j+1;
				  if(jnext>2)
				    {
				      jnext=0;
				    }
				  dcs[j][i]=(uc[jnext]-uc[j])*(vs[i]-vc[j])-
				    (vc[jnext]-vc[j])*(us[i]-uc[j]);
				  dsc[i][j]=(us[inext]-us[i])*(vc[j]-vs[i])-
				    (vs[inext]-vs[i])*(uc[j]-us[i]);
				  if(dsc[i][j]>=R0)
				    {
				      innerc[j]=innerc[j]+1;
				      if(innerc[j]==nspoin)
					{
					  ninerc=ninerc+1;
					}
				    }
				  if(dcs[j][i]>=R0)
				    {
				      inners[i]=inners[i]+1;
				      if(inners[i]==3)
					{
					  niners=niners+1;
					}
				    }
				}
			    }
			  /* B-points */         
			  if(ninerc==3)           /* triangle inside poligon      */
			    {
			      nbpoin=3;
			      for(i=0;i<nbpoin;i++)
				{
				  ub[i]=uc[i];
				  vb[i]=vc[i];
				}
			    }
			  else if(niners==nspoin) /* poligon inside triangle      */
			    {
			      nbpoin=nspoin;
			      for(i=0;i<nbpoin;i++)
				{
				  ub[i]=us[i];
				  vb[i]=vs[i];
				}
			    }
			  else            /* intersection points poligon triangle */
			    {
			      nbpoin=0;
			      for(i=0;i<nspoin;i++)
				{
				  if(inners[i]==3)
				    {
				      ub[nbpoin]=us[i];
				      vb[nbpoin]=vs[i];
				      nbpoin++; 
				    }
				}
			      for(i=0;i<3;i++)  /* grab inner C-points */
				{
				  if(innerc[i]==nspoin)  
				    {
				      ub[nbpoin]=uc[i];
				      vb[nbpoin]=vc[i];
				      nbpoin++;       
				    }
				}       
			      for(i=0;i<nspoin;i++)        /* intersection points   */
				{
				  inext=i+1;
				  if(inext>=nspoin)
				    {
				      inext=0;
				    }
				  for(j=0;j<3;j++)
				    {
				      jnext=j+1;
				      if(jnext>2)
					{
					  jnext=0;
					}
				      if((((dsc[i][j]>EPSILON)&&(dsc[i][jnext]<NEPSILON))||
					  ((dsc[i][j]<NEPSILON)&&(dsc[i][jnext]>EPSILON)))&&
					 (((dcs[j][i]>EPSILON)&&(dcs[j][inext]<NEPSILON))||
					  ((dcs[j][i]<NEPSILON)&&(dcs[j][inext]>EPSILON))))
					//modified by JXiang
					{
					  factor=ABS(dsc[i][j]-dsc[i][jnext]);
					  if(factor<EPSILON)
					    {
					      factor=RP5;
					    }
					  else
					    {
					      factor=ABS(dsc[i][j]/factor);
					    }
					  ub[nbpoin]=(R1-factor)*uc[j]+factor*uc[jnext];
					  vb[nbpoin]=(R1-factor)*vc[j]+factor*vc[jnext];
					  nbpoin++;
					}
				    }
				}        
			      for(i=1;i<nbpoin;i++)
				{
				  if(vb[i]<vb[0])
				    {
				      tmp=vb[i];
				      vb[i]=vb[0];
				      vb[0]=tmp;
				      tmp=ub[i];
				      ub[i]=ub[0];
				      ub[0]=tmp;
				    }
				}
			      for(i=1;i<nbpoin;i++)
				{
				  tmp=ub[i]-ub[0];				 
				  if((tmp<R0)&&(tmp>(-EPSILON)))
				    {
				      tmp=tmp-EPSILON;
				    }
				  else if((tmp>=R0)&&(tmp<EPSILON))
				    {
				      tmp=tmp+EPSILON;
				    }
				  anb[i]=(vb[i]-vb[0]+EPSILON)/tmp;
				}
			      for(i=1;i<nbpoin;i++)  /* sort B-points */
				{
				  for(j=i+1;j<nbpoin;j++)
				    {
				      if(((anb[i]>=R0)&&(anb[j]>=R0)&&(anb[j]<anb[i]))||
					 ((anb[i]<R0)&&((anb[j]>=R0)||(anb[j]<anb[i]))))
					{
					  tmp=vb[i];
					  vb[i]=vb[j];
					  vb[j]=tmp;
					  tmp=ub[i];
					  ub[i]=ub[j];
					  ub[j]=tmp;
					  tmp=anb[i];
					  anb[i]=anb[j];
					  anb[j]=tmp;
					}
				    }
				}
			    }
			  if(nbpoin<3)
			    {
			      continue;
			    }
			  /* Target-plain normal and penetration at B-points */           
			  V3DCro(xnt,ynt,znt,xt[1]-xt[0],yt[1]-yt[0],zt[1]-zt[0],
				 xt[2]-xt[0],yt[2]-yt[0],zt[2]-zt[0]);
			  V3DDot(theigh,xt[3]-xt[0],
				 yt[3]-yt[0],zt[3]-zt[0],xnt,ynt,znt);
			  /* penetration at origin of the e-base and dp/du dp/dv; */
			  V3DDot(peneto,xc[0]-xt[0],yc[0]-yt[0],
				 zc[0]-zt[0],xnt,ynt,znt)
			    V3DDot(penetu,xe[0],ye[0],ze[0],xnt,ynt,znt);
			  V3DDot(penetv,xe[1],ye[1],ze[1],xnt,ynt,znt);
			  peneto=peneto/theigh; 
			  penetu=penetu/theigh; 
			  penetv=penetv/theigh;
			  for(i=0;i<nbpoin;i++)
			    {
			      penetb[i]=peneto+ub[i]*penetu+vb[i]*penetv;
			    }
			  /* force and center of force */
			  forco=R0;
			  uforc=R0;
			  vforc=R0;   
			  for(i=1;i<(nbpoin-1);i++)
			    {
			      penetr=penetb[0]+penetb[i]+penetb[i+1];
			      if(penetr>EPSILON)
				{
				  force=((ub[i]-ub[0])*(vb[i+1]-vb[0])-
				  (vb[i]-vb[0])*(ub[i+1]-ub[0]))*penetr*pepe_n; /* modified by LGuo */

				  fact0=(RP5*penetb[0]+
					 RP25*(penetb[i]+penetb[i+1]))/penetr;
				  facti=(RP5*penetb[i]+RP25*(penetb[0]+
							     penetb[i+1]))/penetr;
				  fact1=R1-fact0-facti;
				  if(ABS(force+forco)>EPSILON)
				    {
				      uforc=(forco*uforc+force*(fact0*ub[0]+
								facti*ub[i]+fact1*ub[i+1]))/(forco+force); 
				      vforc=(forco*vforc+force*(fact0*vb[0]+
								facti*vb[i]+fact1*vb[i+1]))/(forco+force);
				      forco=forco+force;
				    }
				}
			    } 
			  /*             resultant at C-points */
			  for(i=0;i<4;i++)
			    {
			      fc[i]=R0;
			      ft[i]=R0;
			    }
			  tmp=((uc[1]-uc[0])*(vc[2]-vc[0])-
			       (vc[1]-vc[0])*(uc[2]-uc[0]));
			  for(i=0;i<3;i++)
			    {
			      j=i+1;
			      if(j>2)
				{
				  j=0;
				}
			      k=j+1;
			      if(k>2)
				{
				  k=0;
				}
			      fc[k]=forco*(((uc[j]-uc[i])*(vforc-vc[i])-
					    (vc[j]-vc[i])*(uforc-uc[i]))/tmp); 
			    }
			  /*             resultant at T-points  */ 
			  tmp=((ut[1]-ut[0])*(vt[2]-vt[0])-
			       (vt[1]-vt[0])*(ut[2]-ut[0]));
			  inext=-1;
			  if(ABS(tmp)<RP1*theigh)
			    { inext=0;
			      tmp=ABS(ut[1]-ut[0])+ABS(vt[1]-vt[0]); 
			      for(i=0;i<3;i++)
				{
				  j=i+1;
				  if(j>2)
				    {
				      j=0;
				    }
				  if(tmp>(ABS(ut[j]-ut[i])+ABS(vt[j]-vt[i])))
				    {
				      tmp=ABS(ut[j]-ut[i])+ABS(vt[j]-vt[i]);
				      inext=i;
				    }
				}
			      j=inext+1;
			      if(j>2)j=0;
			      if(ABS(zt[j])>ABS(zt[inext]))
				{
				  inext=j;
				}
			      j=inext+1;
			      if(j>2)
				{
				  j=0;
				}
			      k=j+1;
			      if(k>2)
				{
				  k=0;
				}
			      tmp=(ut[k]-ut[j])*(vt[3]-vt[j])-
				(vt[k]-vt[j])*(ut[3]-ut[j]);
			    }
			  for(jnext=0;jnext<3;jnext++)
			    {
			      i=jnext;
			      j=i+1;
			      if(j>2)
				{
				  j=0;
				}
			      k=j+1;
			      if(k>2)
				{
				  k=0;
				}
			      if(i==inext)
				{
				  i=3;
				}
			      if(j==inext)
				{
				  j=3;
				}
			      if(k==inext)
				{
				  k=3;
				}
			      ft[k]=forco*(((ut[j]-ut[i])*(vforc-vt[i])-
					    (vt[j]-vt[i])*(uforc-ut[i]))/tmp);                 
			    }
			  ft[3]=RP25*ft[3];
			  for(i=0;i<3;i++)
			    {
			      ft[i]=ft[i]+ft[3];
			    }
			  /* add forces into global vector */ 

			  if((DABS(fc[0]+fc[1]+fc[2]+fc[3]))>EPSILON)
			    {
			      fnonzero=1;
			    }

			  for(i=0;i<4;i++)
			    {
			      d1nfcx[ipc[i]]=d1nfcx[ipc[i]]+fc[i]*xnc*factor1;
			      d1nfcy[ipc[i]]=d1nfcy[ipc[i]]+fc[i]*ync*factor1;
			      d1nfcz[ipc[i]]=d1nfcz[ipc[i]]+fc[i]*znc*factor1;
			      d1nfcx[ipt[i]]=d1nfcx[ipt[i]]-ft[i]*xnc*factor2;
			      d1nfcy[ipt[i]]=d1nfcy[ipt[i]]-ft[i]*ync*factor2;
			      d1nfcz[ipt[i]]=d1nfcz[ipt[i]]-ft[i]*znc*factor2;
			    }

 			  if(i1ptyp[i1elpr[ielem]]==YTE3TET10ELS)
			    {
             
			      for(i=4;i<10;i++)
				{
				  j = ij[i]; /*Z */
				  k = ik[i];

				  d1nfcx[i1celto[i]]=d1nfcx[i1celto[i]]+
				    (fc[ipcn[j]]+fc[ipcn[k]])*xnc*factor1/R2;
				  d1nfcy[i1celto[i]]=d1nfcy[i1celto[i]]+
				    (fc[ipcn[j]]+fc[ipcn[k]])*ync*factor1/R2;
				  d1nfcz[i1celto[i]]=d1nfcz[i1celto[i]]+
				    (fc[ipcn[j]]+fc[ipcn[k]])*znc*factor1/R2;
				}
			    }
		     
   			  if(i1ptyp[i1elpr[jelem]]==YTE3TET10ELS)
			    {
             
			      for(i=4;i<10;i++)
				{
				  j = ij[i]; /*Z */
				  k = ik[i];
				  d1nfcx[i1telto[i]]=d1nfcx[i1telto[i]]-
				    (ft[iptn[j]]+ft[iptn[k]])*xnc*factor2/R2;
				  d1nfcy[i1telto[i]]=d1nfcy[i1telto[i]]-
				    (ft[iptn[j]]+ft[iptn[k]])*ync*factor2/R2;
				  d1nfcz[i1telto[i]]=d1nfcz[i1telto[i]]-
				    (ft[iptn[j]]+ft[iptn[k]])*znc*factor2/R2;
				}
			    }
			}
		    }

		  /* Calculating sliding frictional force  */


		  /* Calculating relative velocity and direction of normal  */
		  //         fnonzero=0;
		  if(fnonzero==1)
		    {
		      fnx=R0;
		      fny=R0;
		      fnz=R0;
		      fnctot=R0;
		      fnttot=R0;
		      vxc=R0;
		      vyc=R0;
		      vzc=R0;
		      vxt=R0;
		      vyt=R0;
		      vzt=R0;
		      for(i=0;i<4;i++)
			{
			  fcx[i]=d1nfcx[i1celto[i]]-fcx[i];
			  fcy[i]=d1nfcy[i1celto[i]]-fcy[i];
			  fcz[i]=d1nfcz[i1celto[i]]-fcz[i];
			  ftx[i]=d1nfcx[i1telto[i]]-ftx[i];
			  fty[i]=d1nfcy[i1telto[i]]-fty[i];
			  ftz[i]=d1nfcz[i1telto[i]]-ftz[i];
			  fnc[i]=SQRT(fcx[i]*fcx[i]+fcy[i]*fcy[i]+fcz[i]*fcz[i]);
			  fnt[i]=SQRT(ftx[i]*ftx[i]+fty[i]*fty[i]+ftz[i]*ftz[i]);
			  fnx=fnx+fcx[i];
			  fny=fny+fcy[i];
			  fnz=fnz+fcz[i];
			  fnctot=fnctot+fnc[i];
			  fnttot=fnttot+fnt[i];
			}

		      if(fnctot>EPSILON)
			{
			  fn_normal=SQRT(fnx*fnx+fny*fny+fnz*fnz);
			  nx=fnx/fn_normal;
			  ny=fny/fn_normal;
			  nz=fnz/fn_normal;
			  for(i=0;i<4;i++)
			    {
			      vxc=vxc+fnc[i]/fnctot*d1nvcx[i1celto[i]];
			      vyc=vyc+fnc[i]/fnctot*d1nvcy[i1celto[i]];
			      vzc=vzc+fnc[i]/fnctot*d1nvcz[i1celto[i]];
			      vxt=vxt+fnt[i]/fnttot*d1nvcx[i1telto[i]];
			      vyt=vyt+fnt[i]/fnttot*d1nvcy[i1telto[i]];
			      vzt=vzt+fnt[i]/fnttot*d1nvcz[i1telto[i]];
			    }

			  if((ihycon==1)&&((DABS(d1nvx[icoup])
					    +DABS(d1nvy[icoup])+DABS(d1nvz[icoup]))>RP25))
			    {
			      if(DABS((d1nvx[icoup]*nx+d1nvy[icoup]*ny+d1nvz[icoup]*nz)-R1)
				 >EPSILON)
				{
				  V3DCro(pnx,pny,pnz,d1nvx[icoup],d1nvy[icoup],d1nvz[icoup],
					 nx,ny,nz);
				  V3DNor(tmp,pnx,pny,pnz);
				  V3DDot(cosin,nx,ny,nz,d1nvx[icoup],d1nvy[icoup],d1nvz[icoup]);
				  tmp=R1-cosin*cosin;
				  if(tmp<R0)
				    {
				      sine=R0;
				    }
				  else
				    {
				      sine=SQRT(R1-cosin*cosin);
				    }
				  V3DRot(&t1x,&t1y,&t1z,d1t1vx[icoup],d1t1vy[icoup],d1t1vz[icoup],
					 pnx,pny,pnz,cosin,sine);
				  V3DRot(&t2x,&t2y,&t2z,d1t2vx[icoup],d1t2vy[icoup],d1t2vz[icoup],
					 pnx,pny,pnz,cosin,sine);
				}
			      else
				{
				  t1x=d1t1vx[icoup];
				  t1y=d1t1vy[icoup];
				  t1z=d1t1vz[icoup];                      
				  t2x=d1t2vx[icoup];
				  t2y=d1t2vy[icoup];
				  t2z=d1t2vz[icoup];                                 
				}
			    }
			  else
			    {
			      d1deltat1[icoup]=R0;
			      d1deltat2[icoup]=R0;
			      d1deltan[icoup]=R0;
	       
			      if( (ABS(nx)>ABS(ny)) && (ABS(nx)>ABS(nz)) )
				{
				  t1y = t1z = 1; 
				  t1x = - (ny+nz)/nx;
				}
			      else if( (ABS(ny)>ABS(nx)) && (ABS(ny)>ABS(nz)) )
				{
				  t1x = t1z =1;
				  t1y = - (nx+nz) / ny;
				}
			      else
				{
				  t1x = t1y =1;
				  t1z = - (nx+ny) / nz;
				} 

			      V3DNor(tmp,t1x,t1y,t1z);
			      V3DCro(t2x,t2y,t2z,nx,ny,nz,t1x,t1y,t1z);
			    }
            
			  ru=vxt-vxc;
			  rv=vyt-vyc;
			  rw=vzt-vzc;

			  vs1=-rw*t1z-rv*t1y-ru*t1x;
			  vs2=-rw*t2z-rv*t2y-ru*t2x;
			  vn=-rw*nz-rv*ny-ru*nx;
			  deltat1=d1deltat1[icoup]+vs1*dcstec;
			  deltat2=d1deltat2[icoup]+vs2*dcstec;          
			  deltan=d1deltan[icoup]+vn*dcstec;          
			  d1deltan[icoup]=deltan;
 
			  fs1=ktss*deltat1;
			  fs2=ktss*deltat2;
			  fs=fs1*fs1+fs2*fs2;
/*
			  mudd=mud;
			  if((DABS(deltan)>h_cr)&&(DABS(muss-mud)>EPSILON))
			    {
			      mudd=mud+(muss-mud)*pow(DABS(deltan)/h_cr,-2);
			    }
			  if(DABS(mus-R0)>EPSILON)
			    {
			      mu=mudd+(mus-mudd)*exp(-s_vel*SQRT(ru*ru+rv*rv+rw*rw));
			    }
			  else
			    {
			      mu=mudd;
			    }
*/
                          mu=mud;
			  ii1=ii1+1;
			  if(SQRT(fs)>(mu*fn_normal)&&SQRT(fs)>EPSILON)
			    {
			      sfactor=mu*fn_normal/(SQRT(fs)+EPSILON);
			      fs1=fs1*sfactor;
			      fs2=fs2*sfactor;
			      deltat1=fs1/ktss;
			      deltat2=fs2/ktss;
			    }

			  fx=fs1*t1x+fs2*t2x;
			  fy=fs1*t1y+fs2*t2y;
			  fz=fs1*t1z+fs2*t2z;

			  for(i=0;i<4;i++)
			    {
			      d1nfcx[i1celto[i]]=d1nfcx[i1celto[i]]-fx*fnc[i]/fnctot*factor1;
			      d1nfcy[i1celto[i]]=d1nfcy[i1celto[i]]-fy*fnc[i]/fnctot*factor1;
			      d1nfcz[i1celto[i]]=d1nfcz[i1celto[i]]-fz*fnc[i]/fnctot*factor1;
			      d1nfcx[i1telto[i]]=d1nfcx[i1telto[i]]+fx*fnt[i]/fnttot*factor2;
			      d1nfcy[i1telto[i]]=d1nfcy[i1telto[i]]+fy*fnt[i]/fnttot*factor2;
			      d1nfcz[i1telto[i]]=d1nfcz[i1telto[i]]+fz*fnt[i]/fnttot*factor2;
			    }
			  if(i1ptyp[i1elpr[ielem]]==YTE3TET10ELS)
			    {
			      for(i=4;i<10;i++)
				{ 
				  j = ij[i]; /*Z */
				  k = ik[i];
				  d1nfcx[i1celto[i]]=d1nfcx[i1celto[i]]-
				    (fnc[j]+fnc[k])/(R2*fnctot)*fx*factor1;
				  d1nfcy[i1celto[i]]=d1nfcy[i1celto[i]]-
				    (fnc[j]+fnc[k])/(R2*fnctot)*fy*factor1;
				  d1nfcz[i1celto[i]]=d1nfcz[i1celto[i]]-
				    (fnc[j]+fnc[k])/(R2*fnctot)*fz*factor1;
				}
			    }

			  if(i1ptyp[i1elpr[jelem]]==YTE3TET10ELS)
			    {

			      for(i=4;i<10;i++)
				{ 
				  j = ij[i]; /*Z */
				  k = ik[i];
				  d1nfcx[i1telto[i]]=d1nfcx[i1telto[i]]+
				    (fnt[j]+fnt[k])/(R2*fnttot)*fx*factor1;
				  d1nfcy[i1telto[i]]=d1nfcy[i1telto[i]]+
				    (fnt[j]+fnt[k])/(R2*fnttot)*fy*factor1;
				  d1nfcz[i1telto[i]]=d1nfcz[i1telto[i]]+
				    (fnt[j]+fnt[k])/(R2*fnttot)*fz*factor1;
				}
			    }

			  d1deltat1[icoup]=deltat1;
			  d1deltat2[icoup]=deltat2;

			  d1nvx[icoup]=nx;
			  d1nvy[icoup]=ny;
			  d1nvz[icoup]=nz;
 
			  d1t1vx[icoup]=t1x;
			  d1t1vy[icoup]=t1y;
			  d1t1vz[icoup]=t1z;

			  d1t2vx[icoup]=t2x;
			  d1t2vy[icoup]=t2y;
			  d1t2vz[icoup]=t2z;

			  /* end of calculating sliding frictional force  */
			}

		      else
			{
			  d1nvx[icoup]=R0;
			  d1nvy[icoup]=R0;
			  d1nvz[icoup]=R0;
			}
		    }
		  else
		    {
		      d1nvx[icoup]=R0;
		      d1nvy[icoup]=R0;
		      d1nvz[icoup]=R0;
		    }
		}
   
	      if(fnonzero==0)
		{
		  d1deltat1[icoup]=R0;
		  d1deltat2[icoup]=R0;
		}
	      if((icontact==(-1))&&(dmin2>zone2))
		{
		  if(jcoup<0)
		    {
		      /* added by LGuo */
		      i1fcstep[icoup]=-1;

		      i1elcf[ielem]=i1iecn[icoup];
		      i1iecn[icoup]=*i1icff;
		      *i1icff=icoup;
		      icoup=i1elcf[ielem];
		    }
		  else
		    {
		      i1iecn[jcoup]=i1iecn[icoup];
		      i1iecn[icoup]=*i1icff;
		      *i1icff=icoup;           
		      icoup=i1iecn[jcoup];       
		    }
		} 
	      else 
		{
		  jcoup=icoup;
		  icoup=i1iecn[icoup];
		}		
	    }
	}
    }

  //fclose(fp);
}

/*********************PUBLIC***********************/
void Yid(  YDC ydc, YDE  yde, YDI ydi, YDN ydn, YDP ydp     /***  nodal forces  ***/
        )
//  YDC ydc; YDE yde; YDI ydi; YDN ydn; YDP ydp;
{ INT iprop,jprop;
  if(ydi->micoup>0)
  { 
    for(iprop=0;iprop<ydp->nprop;iprop++)
    { if(((ydp->i1ptyp[iprop])==(YTE3TET10ELS))||
         ((ydp->i1ptyp[iprop])==(YTE3TET4ELS)))
      { for(jprop=iprop;jprop<ydp->nprop;jprop++)
        { 
           if(((ydp->i1ptyp[jprop])==(YTE3TET10ELS))||
             ((ydp->i1ptyp[jprop])==(YTE3TET4ELS)))
		   {
			Yid3TET2TET(         /* tetrahedra to tetrahedra */
            yde->nelem ,
            ydc->dcstec,ydi->diezon,iprop,jprop,
            ydi->d1iesl,ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2ncc[2],
			ydn->d2nfc[0],ydn->d2nfc[1],ydn->d2nfc[2],
            ydn->d2nvc[0],ydn->d2nvc[1],ydn->d2nvc[2],
			ydi->d2nv[0],ydi->d2nv[1],ydi->d2nv[2],
            ydi->d1deltat1,ydi->d1deltat2,ydi->d1deltan,
			ydi->d2t1v[0],ydi->d2t1v[1],ydi->d2t1v[2],
			ydi->d2t2v[0],ydi->d2t2v[1],ydi->d2t2v[2],
			ydp->d1pepe,yde->i1elcf,yde->i1elpr,ydp->d1pefr,
			ydp->d1pesf,ydp->d1pevf,ydp->d1pepf,ydp->d1pepsf,
            &(ydi->iiecff),ydi->i1iecn,ydi->i1iect,yde->i2elto,
	    ydp->i1ptyp,ydc->ncstep,ydi->i1fcstep
            );
		   }

       } } }   
  }
}




