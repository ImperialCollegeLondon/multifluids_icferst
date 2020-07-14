#include "frame.h"

int Ycd3TETintersection( double *d1nccx, double *d1nccy, double *d1nccz, 
						INT iele, INT jele, INT **i2elto,double buf) 
{
DBL xorig,yorig,zorig,dct[4],xcent,ycent,zcent,distance[4],min;
DBL xc[4],yc[4],zc[4],xnc,ync,znc,xtmp,ytmp,ztmp,xct[4],yct[4],zct[4];
INT ipc[4],i,j,icons,k,jn,inext,jnext,index,ii,iinext;
DBL dct_tmp,xt[4],yt[4],zt[4],ratio,xline,yline,zline,xv,yv,zv,xint,yint,zint;
if(buf>EPSILON)
{
for(i=0;i<4;i++)
{j=i2elto[jele][i];
xt[i]=d1nccx[j];
yt[i]=d1nccy[j];
zt[i]=d1nccz[j];
}
xcent=0.25*(xt[0]+xt[1]+xt[2]+xt[3]);
ycent=0.25*(yt[0]+yt[1]+yt[2]+yt[3]);
zcent=0.25*(zt[0]+zt[1]+zt[2]+zt[3]);

min=BEPSILON;
for(i=0;i<4;i++)
{ 
   distance[i]=SQRT((xt[i]-xcent)*(xt[i]-xcent)+
			  (yt[i]-ycent)*(yt[i]-ycent)+
			  (zt[i]-zcent)*(zt[i]-zcent));
if(distance[i]<min)min=distance[i];
}
for(i=0;i<4;i++)
{
	ratio=buf/(distance[i]*SIN(min*MYPI/(R6*distance[i])));
    xct[i]=xt[i]+ratio*(xt[i]-xcent);
	yct[i]=yt[i]+ratio*(yt[i]-ycent);
	zct[i]=zt[i]+ratio*(zt[i]-zcent);
  }
}
else
{
for(i=0;i<4;i++)
{j=i2elto[jele][i];
xct[i]=d1nccx[j];
yct[i]=d1nccy[j];
zct[i]=d1nccz[j];
}
}
/*
for( i=0;i<4;i++)
{i2elto[i]=i;
 }
 */
for(icons=0;icons<4;icons++)
{ 
  ipc[0]=icons;
  ipc[1]=1;
  ipc[2]=2;
  if(icons>0)
  { ipc[3]=icons-1;

  }
  else
  { ipc[3]=3;

  }
  if((icons==1)||(icons==2))
  {ipc[1]=3;

  }
  if(icons>1)
  {ipc[2]=0;
  }

for(i=0;i<3;i++)
{
	xc[i]=xct[ipc[i]];
	yc[i]=yct[ipc[i]];
	zc[i]=zct[ipc[i]];
}

 xorig=xc[0]; yorig=yc[0]; zorig=zc[0];
  //xtmp=xt-xorig; ytmp=yt-yorig; ztmp=zt-zorig;
  for(i=0;i<3;i++)
  { xc[i]=xc[i]-xorig; yc[i]=yc[i]-yorig; zc[i]=zc[i]-zorig; 
  } 
  /* contactor normal, e-base and target points in e-base */
  V3DCro(xnc,ync,znc,xc[1],yc[1],zc[1],xc[2],yc[2],zc[2]);
  for(j=0;j<4;j++)
  {i=i2elto[iele][j];  
   xtmp=d1nccx[i]-xorig; ytmp=d1nccy[i]-yorig; ztmp=d1nccz[i]-zorig;
   V3DDot(dct[j],xnc,ync,znc,xtmp,ytmp,ztmp);
  }

  if((dct[0]<=R0)&&(dct[1]<=R0)&&
	  (dct[2]<=R0)&&(dct[3]<=R0))  return 0;  //modified by JXiang

  if(buf<EPSILON)
  {
			  for(i=0;i<3;i++)
			  {
				  if(i==0)jn=3;
				  if(i==1)jn=2;
				  if(i==2)jn=1;
    			  
					  for(j=0;j<jn;j++)
					{ inext=i+j+1; 
					  if(((dct[i]>R0)&&(dct[inext]<R0))||   
						 ((dct[i]<R0)&&(dct[inext]>R0)))
						  ///Modified by JXiang
					  { ratio=ABS(dct[i]-dct[inext]);          
						if(ratio>EPSILON)
						{ ratio=ABS(dct[i]/ratio);
                                                  iinext=i2elto[iele][inext];
                                                  ii=i2elto[iele][i];
						  xint=ratio*d1nccx[iinext]+(R1-ratio)*d1nccx[ii];
						  yint=ratio*d1nccy[iinext]+(R1-ratio)*d1nccy[ii];
						  zint=ratio*d1nccz[iinext]+(R1-ratio)*d1nccz[ii];
						  index=0;
						  for(k=0;k<3;k++)
						  {
							  jnext=k+1;if(k==2)jnext=0;
							  xtmp=xint-xct[ipc[k]];
							  ytmp=yint-yct[ipc[k]];
							  ztmp=zint-zct[ipc[k]];
							  xline=xc[jnext]-xc[k];
							  yline=yc[jnext]-yc[k];
							  zline=zc[jnext]-zc[k];
							  V3DCro(xv,yv,zv,xnc,ync,znc,xline,yline,zline);
							  V3DDot(dct_tmp,xv,yv,zv,xtmp,ytmp,ztmp);
							  if(dct_tmp<R0){index=1;break;}

						  }
						  if(index==0)return 1;
						  
			  } }  }
			  }

 
}


} 
if(buf<EPSILON)return 2;  
else return 1;
  
}




