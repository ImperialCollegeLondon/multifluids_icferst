
/* File   Ydrag.c */
#include "Yproto.h" 
#define V2DNorr(s,x1,y1){s=SQRT((x1)*(x1)+(y1)*(y1));  \
              (x1)=(x1)/s;(y1)=(y1)/s;                 }
void Yring_drag(ydn,ydx,d1nmu,fx,fy,fz,usl,vsl,wsl,axx,ayy,axy,azz,axz,ayz,d1vx,d1vy,d1vz)
YDN ydn; YDX ydx; DBL *d1nmu; DBL *fx; DBL *fy; DBL *fz;
DBL *usl;DBL *vsl; DBL *wsl; DBL *axx; DBL *ayy; DBL *axy; DBL *axz; DBL *ayz; DBL *azz; DBL *d1vx; DBL *d1vy;DBL *d1vz;
{
    INT nelem,i0,i1,i2,i3,i,ielem,ii;
    DBL delta_sq,ii1,jj1,kk1;
    DBL *fxs;
    DBL *fys;
    DBL *fzs;
    DBL *vxs;
    DBL *vys;
    DBL *vzs;
    DBL *nx;
    DBL *ny;
    DBL *nz;
    INT *nv;
     
    DBL norm;
    DBL min_f_ms; //fluid minimum mesh size 

    vxs=TalDBL1(ydx->nnode);
    vys=TalDBL1(ydx->nnode);
    vzs=TalDBL1(ydx->nnode);
    
    nv=TalINT1(ydx->nnode);
    
    
    fxs=TalDBL1(ydx->nnode*2);
    fys=TalDBL1(ydx->nnode*2);
    fzs=TalDBL1(ydx->nnode*2);

    nx=TalDBL1(ydx->nnode*2);
    ny=TalDBL1(ydx->nnode*2);
    nz=TalDBL1(ydx->nnode*2);



    
    DBL sumfx, sumfy, sumfz, total_vol;
    sumfx=0.0;
    sumfy=0.0;
    sumfz=0.0;
    min_f_ms=0.5;
    total_vol=0.0;
   
    // FILE *fp;
    //fp = fopen("./fxfy_viscosity_test.txt", "w");
    // Given fluid velocity, for test only.
  
    //    isuf=ydx->nelem_sur/2;
    //inner ring
    
    for(i=0;i<ydx->nnode;i++)
    {
        nv[i]=0;
        vxs[i]=R0;
        vys[i]=R0;
        vzs[i]=R0;
    }
    
    for(i=0;i<ydx->nsfel;i++)
    {
        ii=ydx->i2xbtojn[i][0];
        i1=ydx->i2elto_rt[i][0];
        nv[i1]++;
        vxs[i1]=vxs[i1]+ydn->d2nvc[0][ii];
        vys[i1]=vys[i1]+ydn->d2nvc[1][ii];
        vzs[i1]=vzs[i1]+ydn->d2nvc[2][ii];

        ii=ydx->i2xbtojn[i][1];
        i1=ydx->i2elto_rt[i][1];
        nv[i1]++;
        vxs[i1]=vxs[i1]+ydn->d2nvc[0][ii];
        vys[i1]=vys[i1]+ydn->d2nvc[1][ii];
        vzs[i1]=vzs[i1]+ydn->d2nvc[2][ii];
        
        ii=ydx->i2xbtojn[i][2];
        i1=ydx->i2elto_rt[i][2];
        nv[i1]++;
        vxs[i1]=vxs[i1]+ydn->d2nvc[0][ii];
        vys[i1]=vys[i1]+ydn->d2nvc[1][ii];
        vzs[i1]=vzs[i1]+ydn->d2nvc[2][ii];
        
    }
    for(i=0;i<ydx->nnode;i++)
    {
        vxs[i]=vxs[i]/nv[i];
        vys[i]=vys[i]/nv[i];
        vzs[i]=vzs[i]/nv[i];
    }
    
    for(i=0;i<ydx->nnode;i++)
    {
        ii=i+ydx->nnode;   
        ii1=vxs[i];
        jj1=vys[i];
        kk1=vzs[i];
        //-PY changed it 
        /*
        ii1=ydn->d2nvc[0][i];
        jj1=ydn->d2nvc[1][i];
        kk1=ydn->d2nvc[2][i];
        */   
        usl[i]=ii1-d1vx[i];
        vsl[i]=jj1-d1vy[i];
        wsl[i]=kk1-d1vz[i];
        
        usl[i]=ii1-d1vx[ii];
        vsl[i]=jj1-d1vy[ii];
        wsl[i]=kk1-d1vz[ii];

        //usl[i]=ii1-1;
        //vsl[i]=jj1-1;
        //wsl[i]=kk1-1;

     nx[i]=ydx->d2areaw[0][i];
     ny[i]=ydx->d2areaw[1][i];
     nz[i]=ydx->d2areaw[2][i];


        //fprintf(fp, "\n ii1, jj1,  kk1: %lf %lf %lf \n",ii1, jj1,  kk1);
        //fprintf(fp, "\n usl[i], vsl[i], wsl[i]: %lf %lf %lf \n",usl[i], vsl[i], wsl[i]);
        //fprintf(fp, "\n d1vx[i], d1vy[i], d1vz[i]: %lf %lf %lf \n",d1vx[i], d1vy[i], d1vz[i]);   
    }
    
    //out ring
    for(i=ydx->nnode;i<2*ydx->nnode;i++)
    {
        ii=i-ydx->nnode;
        //-PY changed it     
        ii1=vxs[ii];
        jj1=vys[ii];
        kk1=vzs[ii];
        //-PY changed it
        /*
        ii1=ydn->d2nvc[0][ii];
        jj1=ydn->d2nvc[1][ii];
        kk1=ydn->d2nvc[2][ii];
        */                                                
        usl[i]=ii1-d1vx[i];
        vsl[i]=jj1-d1vy[i];
        wsl[i]=kk1-d1vz[i];
        //    usl[i]=ii1-1;
        //    vsl[i]=jj1-1;
        //    wsl[i]=kk1-1;


       nx[i]=ydx->d2areaw[0][ii];
       ny[i]=ydx->d2areaw[1][ii];
       nz[i]=ydx->d2areaw[2][ii];


        // fprintf(fp, "\n usl[i], vsl[i], wsl[i]: %lf %lf %lf \n",usl[i], vsl[i], wsl[i]);
        // fprintf(fp, "\n d1vx[i], d1vy[i], d1vz[i]: %lf %lf %lf \n",d1vx[i], d1vy[i], d1vz[i]);
    }
 

    //delta_sq=ydfn->thickness*ydfn->thickness;
    //delta_sq=ydx->dthick*min_f_ms;
    delta_sq=ydx->dthick*ydx->dthick;
    //delta_sq=0.1*min_f_ms*ydx->dthick;
    //  fprintf(fp, "\n delta_sq: %lf \n",delta_sq);

    //     for(i=0;i<ydx->nnode;i++)
    for(i=0;i<ydx->nnode*2;i++)
    {    
        //     nx=ydfn->d1nx[i];i
        //     ny=ydfn->d1ny[i];
        //     nz=ydfn->d1nz[i];
        //nx=ydx->d2areaw[0][i];
        //ny=ydx->d2areaw[1][i];
        //nz=ydx->d2areaw[2][i];
        
        //V2DNorr(norm, nx, ny);

        //-PY changed it
        /*             
        axx[i]=d1nmu[i]*(nx*nx*(2.0-2.0/3.0)+ny*ny)/delta_sq;
        ayy[i]=d1nmu[i]*(nx*nx+ny*ny*(2.0-2.0/3.0))/delta_sq;
        axy[i]=d1nmu[i]*((-2.0/3.0)*nx*ny+ny*nx)/delta_sq;      
        */ 
        /*
        axx[i]=(nx*nx*(2.0-2.0/3.0)+ny*ny)/delta_sq;
        axy[i]=((-2.0/3.0)*nx*ny+ny*nx)/delta_sq;
        //  axz[i]=((-2.0/3.0)*nx*nz+nz*nx)/delta_sq;
        ayy[i]=(nx*nx+ny*ny*(2.0-2.0/3.0))/delta_sq;
        azz[i]=(nx*nx+ny*ny+nz*nz*(2.0-2.0/3.0))/delta_sq;
        axz[i]=((-2.0/3.0)*nx*nz+nz*nx)/delta_sq;
        ayz[i]=((-2.0/3.0)*ny*nz+nz*ny)/delta_sq;
        */

        //-PY correct it
        axx[i]=(nx[i]*nx[i]*(2.0-2.0/3.0)+ny[i]*ny[i]+nz[i]*nz[i])/delta_sq;
        axy[i]=((-2.0/3.0)*nx[i]*ny[i]+ny[i]*nx[i])/delta_sq;
        axz[i]=((-2.0/3.0)*nx[i]*nz[i]+nz[i]*nx[i])/delta_sq;
        ayy[i]=(nx[i]*nx[i]+ny[i]*ny[i]*(2.0-2.0/3.0)+nz[i]*nz[i])/delta_sq;
        ayz[i]=((-2.0/3.0)*ny[i]*nz[i]+nz[i]*ny[i])/delta_sq;
        azz[i]=(nx[i]*nx[i]+ny[i]*ny[i]+nz[i]*nz[i]*(2.0-2.0/3.0))/delta_sq;






        //-PY switch off the drag force
        /*     
        axx[i]=0.0;
        axy[i]=0.0;
        //  axz[i]=((-2.0/3.0)*nx*nz+nz*nx)/delta_sq;
        ayy[i]=0.0;
        azz[i]=0.0;
        axz[i]=0.0;
        ayz[i]=0.0;
        */
        
        //fprintf(fp, "\n nx[i], ny[i], nz[i]: %lf %lf %lf \n",nx[i], ny[i], nz[i]);
        //fprintf(fp, "\n axx[i],axy[i],axz[i],ayy[i]: %lf %lf %lf %lf\n",axx[i],axy[i],axz[i],ayy[i]);
        
        /*
        ii1=ydn->d2nvc[0][ydx->i2xbtojn[i][0]];
        jj1=ydn->d2nvc[1][ydx->i2xbtojn[i][0]];
        kk1=ydn->d2nvc[2][ydx->i2xbtojn[i][0]];
        */ 
        //fprintf(fp, "\n ii1, jj1, kk1: %lf %lf %lf \n",ii1, jj1, kk1);
        //-PY some thing wrong here
        //-PY usl need to be us

        //        fx[i]=axx[i]*ii1+axy[i]*jj1+axz[i]*kk1;
        //        fy[i]=axy[i]*ii1+ayy[i]*jj1+ayz[i]*kk1;
        //        fz[i]=axz[i]*ii1+ayz[i]*jj1+azz[i]*kk1;

        fx[i]=axx[i]*(usl[i]+d1vx[i])+axy[i]*(vsl[i]+d1vy[i])+axz[i]*(wsl[i]+d1vz[i]);
        fy[i]=axy[i]*(usl[i]+d1vx[i])+ayy[i]*(vsl[i]+d1vy[i])+ayz[i]*(wsl[i]+d1vz[i]);
        fz[i]=axz[i]*(usl[i]+d1vx[i])+ayz[i]*(vsl[i]+d1vy[i])+azz[i]*(wsl[i]+d1vz[i]);
       //-PY changed it
        /*
        DBL rho=1000.0;
        fxs[i]=rho*d1nmu[i]*(axx[i]*usl[i]+axy[i]*vsl[i]+axz[i]*wsl[i]);
        fys[i]=rho*d1nmu[i]*(axy[i]*usl[i]+ayy[i]*vsl[i]+ayz[i]*wsl[i]);
        fzs[i]=rho*d1nmu[i]*(axz[i]*usl[i]+ayz[i]*vsl[i]+azz[i]*wsl[i]);
        */

        fxs[i]=d1nmu[i]*(axx[i]*usl[i]+axy[i]*vsl[i]+axz[i]*wsl[i]);
        fys[i]=d1nmu[i]*(axy[i]*usl[i]+ayy[i]*vsl[i]+ayz[i]*wsl[i]);
        fzs[i]=d1nmu[i]*(axz[i]*usl[i]+ayz[i]*vsl[i]+azz[i]*wsl[i]);
             
        /*
        nx=ydfn->d1nx[i];
        ny=ydfn->d1ny[i];
        nz=ydfn->d1nz[i];
        
        axx[i]=d1nmu[i]*(nx*nx*(2.0-2.0/3.0)+ny*ny+nz*nz)/delta_sq;
        ayy[i]=d1nmu[i]*(nx*nx+ny*ny*(2.0-2.0/3.0)+nz*nz)/delta_sq;
        azz[i]=d1nmu[i]*(nx*nx+ny*ny+nz*nz*(2.0-2.0/3.0))/delta_sq;
        axy[i]=d1nmu[i]*((-2.0/3.0)*nx*ny+ny*nx)/delta_sq;
        axz[i]=d1nmu[i]*((-2.0/3.0)*nx*nz+nz*nx)/delta_sq;
        ayz[i]=d1nmu[i]*((-2.0/3.0)*ny*nz+nz*ny)/delta_sq;
        
        fx[i]=axx[i]*usl[i]+axy[i]*vsl[i]+axz[i]*wsl[i];
        fy[i]=axy[i]*usl[i]+ayy[i]*vsl[i]+ayz[i]*wsl[i];
        fz[i]=axz[i]*usl[i]+ayz[i]*vsl[i]+azz[i]*wsl[i];
        */
        //fprintf(fp, "\n fx[i], fy[i], fz[i],d1nmu[i]: %lf %lf %lf %lf \n",fx[i], fy[i], fz[i],d1nmu[i]);
        //fprintf(fp, "\n fxs[i], fys[i], fzs[i],d1nmu[i]: %lf %lf %lf %lf \n",fxs[i], fys[i], fzs[i],d1nmu[i]);
    }
    // average fx and axx
    //fprintf(fp, "\n ydx->nnode: %d\n",ydx->nnode);
/*
    for(i=0;i<2*ydx->nnode;i++)
    {
        if(i<ydx->nnode)
        {
            axx[i]=(axx[i]+axx[i+ydx->nnode])/2.0;
            ayy[i]=(ayy[i]+ayy[i+ydx->nnode])/2.0;
            axy[i]=(axy[i]+axy[i+ydx->nnode])/2.0;
            axz[i]=(axz[i]+axz[i+ydx->nnode])/2.0;
            ayz[i]=(ayz[i]+ayz[i+ydx->nnode])/2.0;
            azz[i]=(azz[i]+azz[i+ydx->nnode])/2.0;
            
            fx[i]=(fx[i]+fx[i+ydx->nnode])/2.0;
            fy[i]=(fy[i]+fy[i+ydx->nnode])/2.0;
            fz[i]=(fz[i]+fz[i+ydx->nnode])/2.0;

            fxs[i]=(fxs[i]+fxs[i+ydx->nnode])/2.0;
            fys[i]=(fys[i]+fys[i+ydx->nnode])/2.0;
            fzs[i]=(fzs[i]+fzs[i+ydx->nnode])/2.0;
        }
        else
        {
            axx[i]=axx[i-ydx->nnode];
            ayy[i]=ayy[i-ydx->nnode];
            axy[i]=axy[i-ydx->nnode];
            axz[i]=axz[i-ydx->nnode];
            ayz[i]=ayz[i-ydx->nnode];
            azz[i]=azz[i-ydx->nnode];
            
            fx[i]=fx[i-ydx->nnode];
            fy[i]=fy[i-ydx->nnode];
            fz[i]=fz[i-ydx->nnode];

            fxs[i]=fxs[i-ydx->nnode];
            fys[i]=fys[i-ydx->nnode];
            fzs[i]=fzs[i-ydx->nnode];
        }
    }
*/

    /*
    for(i=0;i<ydfn->nfnopo;i++)
    {
        if(i<ydfn->nfnopo/2)
        {
            
            axx[i]=(axx[i+ydfn->nfnopo/2]);
            ayy[i]=(ayy[i+ydfn->nfnopo/2]);
            axy[i]=(axy[i+ydfn->nfnopo/2]);
               
            fx[i]=(fx[i+ydfn->nfnopo/2]);
            fy[i]=(fy[i+ydfn->nfnopo/2]);

            fxs[i]=(fxs[i+ydfn->nfnopo/2]);
            fys[i]=(fys[i+ydfn->nfnopo/2]);
        }    
    }
    */

    //-PY from continuous to discontinuous
/*
    for(i=0;i<2*ydx->nnode;i++)

    {

        ii=ydx->i2xbtojn[i][0];
       

        ydn->d2nfc[0][ii]=ydn->d2nfc[0][ii]-fxs[i]*ydx->d1nvol[ii];
        ydn->d2nfc[1][ii]=ydn->d2nfc[1][ii]-fys[i]*ydx->d1nvol[ii];
        ydn->d2nfc[2][ii]=ydn->d2nfc[2][ii]-fzs[i]*ydx->d1nvol[ii];

            
        sumfx=sumfx+ydn->d2nfc[0][ii];
        sumfy=sumfy+ydn->d2nfc[1][ii];
        sumfz=sumfz+ydn->d2nfc[2][ii];
        total_vol= total_vol+ydx->d1nvol[ii];

    }
 
 */
    for(i=0;i<ydx->nsfel;i++)
    {
        ii=ydx->i2xbtojn[i][0];
        i1=ydx->i2elto_rt[i][0];
        ydn->d2nfp[0][ii]=ydn->d2nfp[0][ii]-fxs[i1]*ydn->d1nvol[ii];
        ydn->d2nfp[1][ii]=ydn->d2nfp[1][ii]-fys[i1]*ydn->d1nvol[ii];
        ydn->d2nfp[2][ii]=ydn->d2nfp[2][ii]-fzs[i1]*ydn->d1nvol[ii];

        ii=ydx->i2xbtojn[i][1];
        i1=ydx->i2elto_rt[i][1];
        ydn->d2nfp[0][ii]=ydn->d2nfp[0][ii]-fxs[i1]*ydn->d1nvol[ii];
        ydn->d2nfp[1][ii]=ydn->d2nfp[1][ii]-fys[i1]*ydn->d1nvol[ii];
        ydn->d2nfp[2][ii]=ydn->d2nfp[2][ii]-fzs[i1]*ydn->d1nvol[ii];
  
        ii=ydx->i2xbtojn[i][2];
        i1=ydx->i2elto_rt[i][2];
        ydn->d2nfp[0][ii]=ydn->d2nfp[0][ii]-fxs[i1]*ydn->d1nvol[ii];
        ydn->d2nfp[1][ii]=ydn->d2nfp[1][ii]-fys[i1]*ydn->d1nvol[ii];
        ydn->d2nfp[2][ii]=ydn->d2nfp[2][ii]-fzs[i1]*ydn->d1nvol[ii];
    }
    
    // "w" means that we are going to write on this file
    //fprintf(fp, "\n sumfx, sumfy, sumfz, total_vol: %lf %lf %lf %lf",sumfx, sumfy, sumfz, total_vol);
    //fclose(fp);
    FREE(vzs);
    FREE(vys);
    FREE(vxs);
    FREE(nv);
    FREE(fzs);
    FREE(fys);
    FREE(fxs);  
    FREE(nz);
    FREE(ny);
    FREE(nx);
}

