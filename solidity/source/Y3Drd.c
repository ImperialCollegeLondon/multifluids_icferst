/*! \file Y3Drd.c
 *  \brief read model parameters
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


static  CHR *cdig="0123456789";


/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/


/* read control data    
*/
static void Yrdc(YDC ydc, FILE *finp, CHR *name)
{ 
  CHR *namep;
  
  namep=name+8;
  if(CHRcmp(namep,"MCSTEP",6)==0) 
  { 
    INTr(finp,&(ydc->mcstep));
  }
  else if(CHRcmp(namep,"NCSTEP",6)==0) 
  { 
    INTr(finp,&(ydc->ncstep));
  }
  else if(CHRcmp(namep,"DCGRAX",6)==0) 
  { 
    DBLr(finp,&(ydc->dcgrax));
  }
  else if(CHRcmp(namep,"DCGRAY",6)==0) 
  { 
    DBLr(finp,&(ydc->dcgray));
  }
  else if(CHRcmp(namep,"DCGRAZ",6)==0) 
  { 
    DBLr(finp,&(ydc->dcgraz));
  }
  else if(CHRcmp(namep,"DCSIZC",6)==0) 
  { 
    DBLr(finp,&(ydc->dcsizc));
  }
  else if(CHRcmp(namep,"DCSIZF",6)==0) 
  { 
    DBLr(finp,&(ydc->dcsizf));
  }
  else if(CHRcmp(namep,"DCSIZS",6)==0) 
  { 
    DBLr(finp,&(ydc->dcsizs));
  }
  else if(CHRcmp(namep,"DCSIZV",6)==0) 
  { 
    DBLr(finp,&(ydc->dcsizv));
  }
  else if(CHRcmp(namep,"DCSTEC",6)==0) 
  { 
    DBLr(finp,&(ydc->dcstec));
  }
  else if(CHRcmp(namep,"DCTIME",6)==0) 
  { 
    DBLr(finp,&(ydc->dctime));
  }

  else if(CHRcmp(namep,"DCRMPT",6)==0) 
  { DBLr(finp,&(ydc->dcrmpt));
  }

  else if(CHRcmp(namep,"DCURELX",7)==0) 
  { DBLr(finp,&(ydc->dcurelx));
  }
  else if(CHRcmp(namep,"INITER",6)==0) 
  { 
    INTr(finp,&(ydc->initer));
  }
  else if(CHRcmp(namep,"ICOUTF",6)==0) 
  { 
    INTr(finp,&(ydc->icoutf));
  }
  else if(CHRcmp(namep,"ICOUTI",6)==0) 
  { 
    INTr(finp,&(ydc->icouti));
  }
  else if(CHRcmp(namep,"ICOUTP",6)==0) 
  { 
    INTr(finp,&(ydc->icoutp));
  }
  else if(CHRcmp(namep,"IWFAST",6)==0) 
  { 
    INTr(finp,&(ydc->iwfast));
  }
  else if(CHRcmp(namep,"ISAVE",5)==0) 
  { 
   INTr(finp,&(ydc->isave));
   ydc->iflag=TalINT1(3);
  }

  else
  { 
    CHRw(stderr,"Yrdc: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(103);
  } 
}


/**********************************************************************/
/**********************************************************************/


/* default values 
*/
static void Yrdd(YD yd)
{ 
  YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
  YDO ydo=&(yd->ydo);
  YDP ydp=&(yd->ydp);
  YDB ydb=&(yd->ydb);
  YPAR ypar=&(yd->ypar);
  YSP ysp=&(yd->ysp);
  YDX ydx=&(yd->ydx);
  YDXN ydxn=&(yd->ydxn);
//   YDJ ydj=&(yd->ydj);

  /* Set Control to default   */

  ydc->mcstep=0; ydc->ncstep=0;
  ydc->finp=(FILE*)NULL;   /*Z FILENULL */
  ydc->fcheck=(FILE*)NULL;
  ydc->dcgrax=R0;
  ydc->dcgray=R0;
  ydc->dcgraz=R0;
  ydc->dcsizc=R1;
  ydc->dcsizf=R1;
  ydc->dcsizv=R1;
  ydc->dcstec=R1;
  ydc->dctime=R0;
  ydc->dcrmpt=R0;
  ydc->dcurelx=R1;
  ydc->initer=2;
  ydc->icoutf=0;
  ydc->icouti=0;
  ydc->icoutp=2;
  ydc->iwfast=1;
  ydc->isave=1000;
  ydc->irigid=0; //-PY changed it for rigid tests: 1 for rigid body, 0 for elastic body
  ydc->iflag=INT1NULL;
 

  /* Set Elements to default    */
  yde->melem=0; yde->nelem=0;
  yde->melst=0; yde->nelst=0;
  yde->melno=0; yde->nelno=0;
  yde->i1elcf=INT1NULL; 
  yde->i1elpr=INT1NULL;
  yde->d2elst=DBL2NULL;
  yde->i2elto=INT2NULL;
  yde->i1elty=INT1NULL;
  yde->d3tcs=DBL3NULL;
  yde->d1emct=DBL1NULL;
  yde->i1elbe=INT1NULL;
  yde->i2elfr=INT2NULL;
  yde->nelemi=0;
  yde->i1elbe_h=INT2NULL;

  yde->d2elfs=DBL2NULL;
  yde->i1eljo=INT1NULL;

  yde->i2eljp=INT2NULL;
yde->i1elcft=INT1NULL;

  yde->i2elbnp=INT2NULL;

  yde->d2xap = DBL2NULL;
  yde->d1npore=DBL1NULL;
/* joint structure */ //  AO QL


//-PY changed it for 3D_fracture_coupling_with_multiphase
  ydx->nelem=0;
  ydx->nelno=0;
  ydx->i2elto=INT2NULL;
  ydx->d1nap=DBL1NULL;

  /* Set Interaction to default */
  ydi->micoup=0; ydi->nicoup=0;  
  ydi->iiecff=-2;
  ydi->diedi=BEPSILON;
  ydi->diezon=R0;
  ydi->d1iesl=DBL1NULL;
  ydi->i1iecn=INT1NULL;
  ydi->i1iect=INT1NULL;
  ydi->d1deltat1=DBL1NULL;
  ydi->d1deltat2=DBL1NULL;
  ydi->d1deltan=DBL1NULL;
  ydi->d2nv=DBL2NULL;
  ydi->d2t1v=DBL2NULL;
  ydi->d2t2v=DBL2NULL;

    ydi->i1fcstep=INT1NULL;

  /* Set Nodes to default  */
  ydn->mnodim=0;  ydn->nnodim=0;
  ydn->mnopo=0;  ydn->nnopo=0;
  ydn->nnopi=0;
  ydn->d1nmct=DBL1NULL;
  //ydn->d1nmcf=DBL1NULL; //Z?
  ydn->d2ncc=DBL2NULL;
  ydn->d2nci=DBL2NULL;
  ydn->d2nfc=DBL2NULL;
  ydn->d2nft=DBL2NULL;
  ydn->d2nvc=DBL2NULL;
  ydn->d1nvct=DBL1NULL;
  ydn->i1nobf=INT1NULL;
  ydn->i1nopr=INT1NULL;
  ydn->d2nfd=DBL2NULL;
  ydn->d2nfp=DBL2NULL;
  ydn->d2nfv=DBL2NULL;
  ydn->nneigh=0;
  ydn->i2nnei=INT2NULL;
  ydn->i1nei=INT1NULL;
  ydn->i1ntoC2D=INT1NULL;

  ydn->d1nti=DBL1NULL;
  ydn->d1ntc=DBL1NULL;
  ydn->d1ntp=DBL1NULL;

  /* Set Output to default  */
  ydo->mohys=0;  ydo->nohys=0;
  ydo->dohyp=0.05;  /* 5% accuracy */
  ydo->d1ohyc=DBL1NULL;
  ydo->d1ohyf=DBL1NULL;
  ydo->d1ohys=DBL1NULL;
  ydo->d1ohyt=DBL1NULL;
  ydo->d1ohyx=DBL1NULL;
  ydo->d1ohyy=DBL1NULL;
  ydo->i1ohyt=INT1NULL;

  //-PY changed it for 3D_fracture_coupling_with_multiphase
  ydn->i1nind=INT1NULL;
  ydn->d2nvt=DBL2NULL;
  ydn->d1nvol=DBL1NULL;

  /* Set Properties to default  */
  ydp->mprop=0; ydp->nprop=0;
  ydp->d1peks=DBL1NULL;
  ydp->d1pela=DBL1NULL;
  ydp->d1pemu=DBL1NULL;
  ydp->d1pepe=DBL1NULL;
  ydp->d1pero=DBL1NULL;
  ydp->d1pefr=DBL1NULL;
  ydp->i1ptyp=INT1NULL;
  ydp->d1pesf=DBL1NULL;
  ydp->d1pepsf=DBL1NULL;
  ydp->d1pevf=DBL1NULL;
  ydp->d1pepf=DBL1NULL;
  ydp->d1pefs=DBL1NULL;
  ydp->d1peft=DBL1NULL;
  ydp->d1pegfn=DBL1NULL;
  ydp->d1pegfs=DBL1NULL;

  ydp->i1pejp=INT1NULL;
  ydp->i1pemn=INT1NULL;
  ydp->i1psde=INT1NULL;

//-PY changed it for 3D_fracture_coupling_with_multiphase
  ydp->d1picf=DBL1NULL;
  ydp->d1pcoh=DBL1NULL;

  ydp->d1tcon=DBL1NULL;
  ydp->d1capa=DBL1NULL;
  ydp->d1ctex=DBL1NULL;

  /* Set Boundary Condition to default  */
  ydb->mbcon=0; ydb->nbcon=0;
  ydb->d1bnax=DBL1NULL;
  ydb->d1bnay=DBL1NULL;
  ydb->d1bnaz=DBL1NULL;
  ydb->d1bnfx=DBL1NULL;
  ydb->d1bnfy=DBL1NULL;
  ydb->d1bnfz=DBL1NULL;
  ydb->d1bnvx=DBL1NULL;
  ydb->d1bnvy=DBL1NULL;
  ydb->d1bnvz=DBL1NULL;
  ydb->i1bnvx=INT1NULL;
  ydb->i1bnvy=INT1NULL;
  ydb->i1bnvz=INT1NULL;

  ydb->d1bnap=DBL1NULL;

  ydb->i1bnhf=INT1NULL;
  ydb->d1bnhf=DBL1NULL;
  ydb->i1bntp=INT1NULL;
  ydb->d1bntp=DBL1NULL;
  ydb->d1bcvt=DBL1NULL;
  ydb->d1bcvc=DBL1NULL;
  ydb->i1bcvt=INT1NULL;

//-PY changed it for 3D_fracture_coupling_with_multiphase
  ypar->nelem=0;   /*!< maximum number of properties            */
  ypar->nnode=0;   /*!< actual number of properties             */
  ypar->i2elto=INT2NULL;
  ypar->d2xyz=DBL2NULL; /*!< [mprop] amplitude of acceleration x     */
  ypar->d2gxyz=DBL2NULL;
  ypar->d1cg=DBL1NULL;
  ypar->d1dst=DBL1NULL;
  ypar->radius=R0;
  ypar->max=R0;
  ypar->d1lv=DBL1NULL;
  ypar->nindex=0;
  ypar->i1elbe=INT1NULL;
  ypar->r=5.0;
  ypar->ir=1;
  ypar->d1vc=DBL1NULL;
  //ypar->zmax=R0;


  /* Set Surface boundary to default */
  ydx->i2xbto=INT2NULL;
  ydx->i2xbtosf=INT2NULL;
  ydx->nsfel=0;
  ydx->nsfel_old=0;

//-PY changed it for 3D_fracture_coupling_with_multiphase
   ydx->nelem_t=0;
   ydx->nelno=0;
   ydx->nnode=0;
   ydx->i2elto_t=INT2NULL;
   ydx->i2elto_rt=INT2NULL;

  ydx->nelno=3;
  ydx->ncmel=0;
  ydx->i2cvto=INT2NULL;
  ydx->i2elcon=INT2NULL;
  ydx->d2area=DBL2NULL;
  ydx->d2normal=DBL2NULL;
  ydx->i2eltoxb=INT2NULL;
  ydx->i2xbtojn=INT2NULL;
  ydx->i1xbtoel=INT1NULL;
  
  ydxn->i2xbjno=INT2NULL;
  ydxn->i1xbjo=INT1NULL;
  ydxn->ncmno=0;
  ydxn->d2xncc=DBL2NULL;
  ydxn->d2xnvc=DBL2NULL;
  ydxn->d1shap = DBL1NULL; //-ao asiri added this changed from TAL1DBL
  ydxn->d1xntc = DBL1NULL;

    ydx->i2elto_r=INT2NULL;
    ydx->i1r2s=INT1NULL;
    ydx->i1s2r=INT1NULL;
    ydx->d2ncc=DBL2NULL;
    ydx->d1nvol=DBL1NULL;
    //	ydx->d1ndrag=DBL1NULL;
    ydx->d2nfv=DBL2NULL;
    ydx->dthick=R0;
    ydx->d2areaw=DBL2NULL;


//-PY changed it for 3D_fracture_coupling_with_multiphase


  ysp->mspd=0;
  ysp->nspd=0; /*maximum (actual) number of particle*/

  ysp->d1mass=DBL1NULL;
  ysp->d2pp=DBL2NULL;    /*[nnodim][msd] sphere_position_x                     */
  ysp->d2poldp=DBL2NULL;    /*[nnodim][msd] sphere_position_x                     */
  ysp->d2pvt=DBL2NULL;     /*[nnodim][msd] sphere_velocity_x              */
//  ysp->d2pfn=DBL2NULL;    /*[nnodim][msd] sphere_normal force_x         */
//  ysp->d2pft1=DBL2NULL;   /*[nnodim][msd] sphere_tangential force_x1    */
//  ysp->d2pft2=DBL2NULL;   /*[nnodim][msd] sphere_tangential force_x2    */

  ysp->d2ptfc=DBL2NULL;    /*[nnodim][msd] sphere_total_contact force_x */
  ysp->d2ptf=DBL2NULL;    /*[nnodim][msd] sphere_total_force_x */
  ysp->d2pa=DBL2NULL;     /*[nnodim][msd] sphere_rotation angel_x      */
  ysp->d2pvr=DBL2NULL;     /*[nnodim][msd] sphere_rotation velocity_x      */
  ysp->d2pm=DBL2NULL;     /*[nnodim][msd] sphere_moment_x */
  ysp->d2pmd=DBL2NULL; 
  ysp->d2pfd=DBL2NULL; 
//  DBL    dsengy;
  ysp->d2prine=DBL2NULL;  /*[nnodim] pricinple inertia                         */
  ysp->d2prinn1=DBL2NULL;  /*[nnodim] pricinple axe                         */
  ysp->d2prinn2=DBL2NULL;  /*[nnodim] pricinple axe                         */
  ysp->d2prinn3=DBL2NULL;  /*[nnodim] pricinple axe                         */
  ysp->d2poldvr=DBL2NULL;     /*[nnodim][msd] sphere_rotation velocity_x      */
  ysp->i1sppr=INT1NULL;
  ysp->i1con=INT1NULL;






}


/**********************************************************************/
/**********************************************************************/

/* read data elements   
 */

//-PY changed it for 3D_fracture_coupling_with_multiphase
static void
Yrde (YDE yde, FILE *finp, CHR *name, INT ntelem)
{
  CHR *namep;

  INT i, j;

  namep = name + 8;
  if (CHRcmp(namep,"MELEM",5) == 0)
    {
      INTr(finp, &(yde->melem));



      yde->melem = yde->melem + ntelem;

      yde->i2elbnp = TalINT2 (yde->melem, 4);
      yde->d1npore = TalDBL1 (yde->melem);
      for (i = 0; i < yde->melem; i++)
	{
	  yde->d1npore[i] = R0;
	  for (j = 0; j < 4; j++)
	    {
	      yde->i2elbnp[i][j] = -1;
	    }
	}


  yde->i1elcft=TalINT1(yde->melem);
    for(i=0;i<yde->melem;i++)
      {
	yde->i1elcft[i]=-1;
      }
    }
  else if (CHRcmp(namep,"NELEM",5) == 0)
    {
      INTr(finp, &(yde->nelem));
      yde->nelemi = yde->nelem;

    }
  else if (CHRcmp(namep,"MELST",5) == 0)
    {
      INTr(finp, &(yde->melst));
    }
  else if (CHRcmp(namep,"NELST",5) == 0)
    {
      INTr(finp, &(yde->nelst));
    }
  else if (CHRcmp(namep,"MELNO",5) == 0)
    {
      INTr(finp, &(yde->melno));
    }
  else if (CHRcmp(namep,"NELNO",5) == 0)
    {
      INTr(finp, &(yde->nelno));
    }
  else if (CHRcmp(namep,"I1ELCF",6) == 0)
    {
      TformINT1 (finp, -1, yde->melem, &(yde->i1elcf));
    }
  else if (CHRcmp(namep,"I1ELBE",6) == 0)
    {
      TformINT1 (finp, -1, yde->melem, &(yde->i1elbe));
    }
  else if (CHRcmp(namep,"I1ELPR",6) == 0)
    {
      TformINT1 (finp, 0, yde->melem, &(yde->i1elpr));

      yde->i1eljo = TalINT1 (yde->melem);
      for (i = 0; i < yde->melem; i++)
	{
	  yde->i1eljo[i] = (-1);
	}
    }
  else if (CHRcmp(namep,"D2ELST",6) == 0)
    {
      TformDBL2 (finp, R0, yde->melst, yde->melem, &(yde->d2elst));
    }
  else if (CHRcmp(namep,"I2ELTO",6) == 0)
    {
      //Z transformed array order
      //Z TformINT2(finp,-1,yde->melno,yde->melem,&(yde->i2elto));
      TformINT2_inv (finp, -1, yde->melem, yde->melno, &(yde->i2elto));
      
       yde->d2elfs=TalDBL2(yde->melem,3);
      yde->d2xap = TalDBL2 (yde->melem, yde->nelno);
       for(i=0;i<yde->melem;i++)
       {
       for(j=0;j<3;j++)
       {
       yde->d2elfs[i][j]=R0;
	      yde->d2xap[i][j] = R0;
       }
       }       
      yde->i2eljp = TalINT2 (yde->melem, yde->melno);
      for (i = 0; i < yde->melem; i++)
	{
	  for (j = 0; j < yde->melno; j++)
	    {
	      yde->i2eljp[i][j] = -1;
	    }
	}
      yde->i1elty = TalINT1 (yde->melem);
      for (i = 0; i < yde->melem; i++)
	{
	  yde->i1elty[i] = -1;
	}
    }
  else if (CHRcmp(namep,"I2ELBNP",6) == 0)
    {
      //Z transformed array order
      //Z TformINT2(finp,-1,yde->melno,yde->melem,&(yde->i2elto));
      TformINT2_inv (finp, -1, yde->melem, yde->melno, &(yde->i2elbnp));
    }
  else if (CHRcmp(namep,"D3TCS",5) == 0)
    {
      //Z transformed array order
      //Z TformDBL3(finp,R0,3,3,yde->melem,&(yde->d3tcs));
      TformDBL3_inv (finp, R0, yde->melem, 3, 3, &(yde->d3tcs));
    }
  else if (CHRcmp(namep,"D1EMCT",6) == 0)
    {
      TformDBL1 (finp, R0, yde->melem, &(yde->d1emct));
    }
  else
    {
      CHRw(stderr, "Yrde: unknown name: ");
      CHRw(stderr, name);
      CHRwcr(stderr);
      exit (104);
    }
}


/**********************************************************************/
/**********************************************************************/


/* read data interaction 
*/
static void Yrdi(YDI ydi, FILE *finp, CHR *name)
{ 
  CHR *namep;
  INT i;
  namep=name+8;
  if(CHRcmp(namep,"MICOUP",6)==0) 
  {  
	  INTr(finp,&(ydi->micoup));

	  ydi->i1fcstep=TalINT1(ydi->micoup);
	  for(i=0;i<ydi->micoup;i++)
	    {
	      ydi->i1fcstep[i]=-1;
	    }
  }
  else if(CHRcmp(namep,"NICOUP",6)==0) 
  {  
	  INTr(finp,&(ydi->nicoup));
  }
  else if(CHRcmp(namep,"IIECFF",6)==0) 
  {  
	  INTr(finp,&(ydi->iiecff));
  }
  else if(CHRcmp(namep,"DIEDI",5)==0) 
  {  
	  DBLr(finp,&(ydi->diedi));
  }
  else if(CHRcmp(namep,"DIEZON",6)==0) 
  {  
	  DBLr(finp,&(ydi->diezon));
  }
  else if(CHRcmp(namep,"D1IESL",6)==0) 
  { 
	  TformDBL1(finp,R0,ydi->micoup,&(ydi->d1iesl)); 
  }
  else if(CHRcmp(namep,"I1IECN",6)==0) 
  { 
	  TformINT1(finp,-1,ydi->micoup,&(ydi->i1iecn)); 
  }
  else if(CHRcmp(namep,"I1IECT",6)==0) 
  { 
	  TformINT1(finp,-1,ydi->micoup,&(ydi->i1iect)); 
  }
  else if(CHRcmp(namep,"D1DELTAT1",9)==0) 
  { 
	  TformDBL1(finp,R0,ydi->micoup,&(ydi->d1deltat1)); 
  }
  else if(CHRcmp(namep,"D1DELTAT2",9)==0) 
  { 
	  TformDBL1(finp,R0,ydi->micoup,&(ydi->d1deltat2)); 
  }
  else if(CHRcmp(namep,"D1DELTAN",8)==0) 
  { 
	  TformDBL1(finp,R0,ydi->micoup,&(ydi->d1deltan)); 
  }
  else if(CHRcmp(namep,"D2NV",4)==0) 
  { 
	  TformDBL2(finp,R0,3,ydi->micoup,&(ydi->d2nv));
  }
  else if(CHRcmp(namep,"D2T1V",5)==0) 
  { 
	  TformDBL2(finp,R0,3,ydi->micoup,&(ydi->d2t1v));
  }
  else if(CHRcmp(namep,"D2T2V",5)==0) 
  { 
	  TformDBL2(finp,R0,3,ydi->micoup,&(ydi->d2t2v));
  }

  else
  { 
    CHRw(stderr,"Yrdi: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
  } 
}


/**********************************************************************/
/**********************************************************************/


/* read data nodes       
*/
static void Yrdn(YDN ydn, FILE *finp, CHR *name, INT ntnode)
{ 
  CHR *namep;

//-PY changed it for 3D_fracture_coupling_with_multiphase
  CHR name1[300];
  INT i,j,ndim,nnode;

  namep=name+8;
  if(CHRcmp(namep,"MNODIM",5)==0) 
  {  
    INTr(finp,&(ydn->mnodim));
  }
  else if(CHRcmp(namep,"NNODIM",5)==0) 
  {  
    INTr(finp,&(ydn->nnodim));
  }
  else if(CHRcmp(namep,"MNOPO",5)==0) 
  {  
    INTr(finp,&(ydn->mnopo));

    ydn->i2nnei=TalINT2(ydn->mnopo,100);
    for(i=0;i<ydn->mnopo;i++)
      {
        for(j=0;j<100;j++)
          {
            ydn->i2nnei[i][j]=-1;
          }
      }
    ydn->i1nei=TalINT1(ydn->mnopo);
    for(i=0;i<ydn->mnopo;i++)
      {
	ydn->i1nei[i]=0;
      }
//-PY changed it for 3D_fracture_coupling_with_multiphase
    ydn->mnopo=ydn->mnopo+ntnode; //-ao there could be a memory issue with this 
	ydn->i1nind=TalINT1(ydn->mnopo);

	ydn->d1nvct=TalDBL1(ydn->mnopo);
  }
  else if(CHRcmp(namep,"NNOPO",5)==0) 
  {  
    INTr(finp,&(ydn->nnopo));
    ydn->nnopi=ydn->nnopo;
  }
  else if(CHRcmp(namep,"D1NMCT",6)==0) 
  { 
    TformDBL1(finp,R0,ydn->mnopo,&(ydn->d1nmct));
  }
  else if(CHRcmp(namep,"D2NCC",5)==0) 
  { 
    TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2ncc));
	  ydn->d2nfv=TalDBL2(ydn->mnodim,ydn->mnopo);
	  ydn->d2nfp=TalDBL2(ydn->mnodim,ydn->mnopo);
	  ydn->d2nfd=TalDBL2(ydn->mnodim,ydn->mnopo);
  }
  else if(CHRcmp(namep,"D2NCI",5)==0) 
  { 
    TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nci)); 
  }
  else if(CHRcmp(namep,"D2NFC",5)==0) 
  {
    TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nfc)); 
  }  
  else if(CHRcmp(namep,"D2NFT",5)==0) 
  { 
//    TformDBL2_inv(finp,R0,ydn->mnopo,ydn->mnodim,&(ydn->d2nft));
    TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nft)); 
  }
  else if(CHRcmp(namep,"D2NVC",5)==0) 
  { 
    TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nvc)); 


//-PY changed it for 3D_fracture_coupling_with_multiphase
    ydn->d2nvt=TalDBL2(ydn->mnodim,ydn->mnopo);
    ydn->d1nvol=TalDBL1(ydn->mnopo);

      ydn->i1ntoC2D = TalINT1 (ydn->mnopo); //asiri added
    for(i=0;i<ydn->mnopo;i++)
    {
	    ydn->i1ntoC2D[i] = 0; //asiri added
        ydn->d1nvol[i]=R0;

        for(j=0;j<ydn->mnodim;j++)
        {
         ydn->d2nvt[j][i]=ydn->d2nvc[j][i];
         }
    }
  }

//-PY changed it for 3D_fracture_coupling_with_multiphase
  else if(CHRcmp(namep,"D2NVT",5)==0)
  { CHRr(finp,name1);
    INTr(finp,&(ydn->mnodim));
    INTr(finp,&(ydn->mnopo));
    for(i=0;i<ydn->mnopo;i++)
    { for(j=0;j<ydn->mnodim;j++)
        {
         DBLr(finp,&(ydn->d2nvt[j][i]));
         }

        }
  }




  else if(CHRcmp(namep,"I1NOBF",6)==0) 
  { 
    TformINT1(finp,0,ydn->mnopo,&(ydn->i1nobf));
  }
  else if(CHRcmp(namep,"I1NOPR",6)==0) 
  { 
    TformINT1(finp,1,ydn->mnopo,&(ydn->i1nopr));
  }
  else if(CHRcmp(namep,"D1NTI",6)==0) 
  	{ 
    TformDBL1(finp,R0,ydn->mnopo,&(ydn->d1nti));
  	}
  else if(CHRcmp(namep,"D1NTC",6)==0) 
  	{ 
    TformDBL1(finp,R0,ydn->mnopo,&(ydn->d1ntc));
  	}
  else if(CHRcmp(namep,"D1NTP",6)==0) 
  	{ 
    TformDBL1(finp,R0,ydn->mnopo,&(ydn->d1ntp));
  	}
  else
  { 
    CHRw(stderr,"Yrdn: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(105);
  } 
}


/**********************************************************************/
/**********************************************************************/


/* read output specification     
*/
static void Yrdo(YDO ydo, FILE *finp, CHR *name)
{ 
  INT i; CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MOHYS",5)==0) 
  {  
    INTr(finp,&(ydo->mohys))
    i=(ydo->mohys)*sizeof(FILE*);
    if(i>0)ydo->f2ohyf=(FILE**)MALLOC(i);
    for(i=0;i<(ydo->mohys);i++)
    { 
      ydo->f2ohyf[i]=FILENULL;
    }
  }
  else if(CHRcmp(namep,"NOHYS",5)==0) 
  {  
    
    INTr(finp,&(ydo->nohys));
  }
  else if(CHRcmp(namep,"DOHYP",5)==0) 
  { 
    DBLr(finp,&(ydo->dohyp));
  }
  else if(CHRcmp(namep,"D1OHYS",6)==0) 
  { 
    TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohys));
  }
  else if(CHRcmp(namep,"D1OHYC",6)==0) 
  { 
    TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyc));
  }
  else if(CHRcmp(namep,"D1OHYF",6)==0) 
  { 
    TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyf));
  }
  else if(CHRcmp(namep,"D1OHYT",6)==0) 
  { 
    TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyt));
  }
  else if(CHRcmp(namep,"D1OHYX",6)==0) 
  { 
    TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyx));
  }
  else if(CHRcmp(namep,"D1OHYY",6)==0) 
  { 
    TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyy));
  }
  else if(CHRcmp(namep,"D1OHYZ",6)==0) 
  { 
    TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyz));
  }
  else if(CHRcmp(namep,"I1OHYT",6)==0) 
  { 
    TformINT1(finp,-1,ydo->mohys,&(ydo->i1ohyt));
  }
  else
  { 
    CHRw(stderr,"Yrdo: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(105);
  } 
}


/**********************************************************************/
/**********************************************************************/


/* read data properties  
*/
static void Yrdp(YDP ydp, FILE *finp, CHR *name)
{ 
  CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MPROP",5)==0) 
  {  
    INTr(finp,&(ydp->mprop));
  }
  else if(CHRcmp(namep,"NPROP",5)==0) 
  {  
    INTr(finp,&(ydp->nprop));
  }
  else if(CHRcmp(namep,"D1PEKS",6)==0) 
  { 
    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1peks)); 
  }
  else if(CHRcmp(namep,"D1PELA",6)==0) 
  { 
    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pela)); 
  }
  else if(CHRcmp(namep,"D1PEMU",6)==0) 
  { 
    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pemu)); 
  }
  else if(CHRcmp(namep,"D1PEPE",6)==0) 
  { 
    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pepe)); 
  }
  else if(CHRcmp(namep,"D1PERO",6)==0) 
  { 
    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pero)); 
  }
  else if(CHRcmp(namep,"D1PEFR",6)==0) 
  { 
    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pefr)); 
  }
  else if(CHRcmp(namep,"D1PESF",6)==0) 
  { 
	TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pesf)); 
  }
  else if(CHRcmp(namep,"D1PEPSF",7)==0) 
  { 
	TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pepsf)); 
  }
  else if(CHRcmp(namep,"D1PEVF",6)==0) 
  { 
	TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pevf)); 
  }
  else if(CHRcmp(namep,"D1PEPF",6)==0) 
  { 
	TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pepf)); 
  }
  else if(CHRcmp(namep,"D1PEFS",6)==0)
  {
        TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pefs));
  }
  else if(CHRcmp(namep,"D1PEFT",6)==0)
  {
        TformDBL1(finp,R0,ydp->mprop,&(ydp->d1peft));

        ydp->d1pefs=TalDBL1(ydp->mprop);
  }
  else if(CHRcmp(namep,"D1PEGFN",7)==0)
  {
        TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pegfn));
  }
  else if(CHRcmp(namep,"D1PEGFS",7)==0)
  {
      TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pegfs));
  }
  else if(CHRcmp(namep,"D1PCOH",6)==0)
  {
        TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pcoh));
  }
  else if(CHRcmp(namep,"D1PICF",6)==0)
  {
        TformDBL1(finp,R0,ydp->mprop,&(ydp->d1picf));
  }

  else if(CHRcmp(namep,"I1PTYP",6)==0) 
  { 
        TformINT1(finp,-1,ydp->mprop,&(ydp->i1ptyp));
  }
  else if(CHRcmp(namep,"I1PEJP",6)==0)
  {
        TformINT1(finp,0,ydp->mprop,&(ydp->i1pejp));
  }
  else if(CHRcmp(namep,"I1PEMN",6)==0)
  {
        TformINT1(finp,0,ydp->mprop,&(ydp->i1pemn));
  }
  else if(CHRcmp(namep,"I1PSDE",6)==0)
  {
        TformINT1(finp,0,ydp->mprop,&(ydp->i1psde));
  }
  else if(CHRcmp(namep,"D1TCON",6)==0) 
  { 
	    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1tcon));
  }
  else if(CHRcmp(namep,"D1CAPA",6)==0) 
  { 
	    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1capa));
  }  
  else if(CHRcmp(namep,"D1CTEX",6)==0) 
  { 
	    TformDBL1(finp,R0,ydp->mprop,&(ydp->d1ctex));
  }
  else
  { 
    CHRw(stderr,"Yrdp: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
  } 
}


/**********************************************************************/
/**********************************************************************/


/* read data properties  
*/
static void Yrdb(YDB ydb, FILE *finp, CHR *name)
{ 
  CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MBCON",5)==0) 
  {  
    INTr(finp,&(ydb->mbcon));
  }
  else if(CHRcmp(namep,"NBCON",5)==0) 
  {  
    INTr(finp,&(ydb->nbcon));
  }
  else if(CHRcmp(namep,"D1BNAX",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnax)); 
  }
  else if(CHRcmp(namep,"D1BNAY",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnay)); 
  }
  else if(CHRcmp(namep,"D1BNAZ",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnaz)); 
  }
  else if(CHRcmp(namep,"D1BNAP",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnap)); 
  }

  else if(CHRcmp(namep,"D1BNFX",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnfx));
  }
  else if(CHRcmp(namep,"D1BNFY",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnfy));
  }
  else if(CHRcmp(namep,"D1BNFZ",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnfz));
  }
  else if(CHRcmp(namep,"I1BNVX",6)==0) 
  { 
    TformINT1(finp,0,ydb->mbcon,&(ydb->i1bnvx));
  }
  else if(CHRcmp(namep,"I1BNVY",6)==0) 
  { 
    TformINT1(finp,0,ydb->mbcon,&(ydb->i1bnvy));
  }
  else if(CHRcmp(namep,"I1BNVZ",6)==0) 
  { 
    TformINT1(finp,0,ydb->mbcon,&(ydb->i1bnvz));
  }
  else if(CHRcmp(namep,"D1BNVX",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnvx));
  }
  else if(CHRcmp(namep,"D1BNVY",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnvy));
  }
  else if(CHRcmp(namep,"D1BNVZ",6)==0) 
  { 
    TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnvz));
  }
	else if(CHRcmp(namep,"D1BNTP",6)==0)
  	{
	  TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bntp));
  	}
  else if(CHRcmp(namep,"D1BNHF",6)==0)
  	{
	  TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnhf));
  	}
  else if(CHRcmp(namep,"I1BNHF",6)==0) 
  	{ 
    TformINT1(finp,0,ydb->mbcon,&(ydb->i1bnhf));
  	}
  else if(CHRcmp(namep,"I1BNTP",6)==0) 
  	{ 
    TformINT1(finp,0,ydb->mbcon,&(ydb->i1bntp));
  	}
  else if(CHRcmp(namep,"D1BCVT",6)==0)
  	{ 
		TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bcvt));
  	}
  else if(CHRcmp(namep,"D1BCVC",6)==0)
  	{ 
		TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bcvc));
  	}
  else if(CHRcmp(namep,"I1BCVT",6)==0) 
	  { 
		TformINT1(finp,0,ydb->mbcon,&(ydb->i1bcvt));
	  }
  else
  { 
    CHRw(stderr,"Yrdp: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
  }


  //  FILE *fp1;
  //  fp1=fopen("init_boundary.txt","w+");
  //  fclose(fp1);
}

//-PY changed it for 3D_fracture_coupling_with_multiphase
static void Yrdx(ydx, yde, finp,name,nnopo)
YDX    ydx;   YDE   yde; FILE *finp; CHR *name; INT nnopo; 
{
    CHR *namep; 
    INT i,j,index,nnode,n0,k,*i1node;

    namep=name+8;
    if(CHRcmp(namep,"NELEM",5)==0)
    {  
        INTr(finp,&(ydx->nelem_t));
    }
    else if(CHRcmp(namep,"NELNO",5)==0)
    {  
        INTr(finp,&(ydx->nelno));
    }
    else if(CHRcmp(namep,"D1NAP",5)==0) 
    { 
        TformDBL1(finp,R0,yde->nelem*4+1,&(ydx->d1nap));
    }

    else if(CHRcmp(namep,"DTHICK",6)==0)
    {  
        DBLr(finp,&(ydx->dthick));
    }


//-PY changed it for 3D_fracture_coupling_with_multiphase
    else if(CHRcmp(namep,"I2ELTO",6)==0) 
    {
        //         TformINT2(finp,-1,ydx->nelno+1,ydx->nelem_t,&(ydx->i2elto_t));
	  	TformINT2_inv(finp,-1,ydx->nelem_t,ydx->nelno+1,&(ydx->i2elto_t));
     
//            ydx->d1nap=TalDBL1(ydx->nelem_t*4+1);
	    ydx->d2area=TalDBL2(ydx->nelno+1,4*yde->melem+1);
        ydx->d2normal=TalDBL2(4,4*yde->melem+1);


//Jiansheng Xiang new one
        ydx->i2elto_r=TalINT2(12*yde->nelem+1,5);
        ydx->i2elto_rt=TalINT2(8*yde->nelem+1,5);
        ydx->d2areaw=TalDBL2(ydx->nelno+1,16*yde->nelem+1);

      ydx->i2xbto=TalINT2(4*yde->melem,yde->melno);
      ydx->i2xbtosf=TalINT2(4*yde->melem,yde->melno);
      ydx->i2xbtojn=TalINT2(4*yde->melem,yde->melno);
      ydx->i2cvto=TalINT2(yde->melem,yde->melno);
      ydx->i2elcon=TalINT2(yde->melem,yde->melno);
      ydx->i2eltoxb=TalINT2(yde->melem,yde->melno);
      ydx->i1xbtoel=TalINT1(yde->melem);


  for(i=0;i<yde->melem;i++)
    {
    ydx->i1xbtoel[i]=-1;
      for(j=0;j<yde->melno;j++)
	{
	  ydx->i2xbto[i][j]=-1;
	  ydx->i2xbtosf[i][j]=-1;
	  ydx->i2xbtojn[i][j]=-1;
	  ydx->i2cvto[i][j]=yde->i2elto[i][j];
	  ydx->i2elcon[i][j]=-1;
	  ydx->i2eltoxb[i][j]=-1;
	}
    }
    ydx->ncmel=yde->nelem;

        i1node=TalINT1(6*ydx->nelem_t);
        for(i=0;i<6*ydx->nelem_t;i++)
        {
            i1node[i]=-1;
        }
        nnode=0;
        for (i=0; i<ydx->nelem_t; i++)
        {
            for(j=0;j<3;j++)
            {
                index=0;
                n0=ydx->i2elto_t[i][j];
                for (k=0; k<nnode; k++)
                {
                    if(n0==i1node[k])
                    {
                        index=1;
                        break;
                    }
                }
                if (index==0)
                {
                    i1node[nnode]=n0;
                    //				  i1n2s[n0]=nnode;
                    nnode++;
                    //				  INTw(ftmp,nnode,5);
                    //				  CHRwcr(ftmp);
                }
            }
        }
   
        ydx->d2ncc=TalDBL2(4,8*yde->nelem+1);
        ydx->d2nfv=TalDBL2(4,8*yde->nelem+1);
      
        ydx->i1r2s=TalINT1(8*yde->nelem+1);
        ydx->i1s2r=TalINT1(nnopo+1);
        
        
        ydx->nnode=nnode;
        ydx->d1nvol=TalDBL1(8*yde->nelem+1);
        
        
        
        for(i=0;i<8*yde->nelem;i++)
        {
            for(j=0;j<3;j++)
            {
                ydx->d2nfv[j][i]=R0;
                ydx->d2ncc[j][i]=R0;
            }
            ydx->i1r2s[i]=-1;
            ydx->d1nvol[i]=R0;
	}
        FREE(i1node);

        // revise jul13
    }
    else
    { 
        CHRw(stderr,"Yrdx: unknown name: ");
        CHRw(stderr,name);
        CHRwcr(stderr);
        exit(104);
    }
    CHRw(stderr,"05022017-2\n"); 
}









/* initialise surface node */
static void Yrdxn(YDXN ydxn, YDN ydn,DBL dthick)
{
  INT i,j;


  ydxn->i2xbjno=TalINT2(ydn->mnopo,100);
 
  for(i=0;i<ydn->mnopo;i++)
    {
      for(j=0;j<100;j++)
	{
	  ydxn->i2xbjno[i][j]=-1;
	}
    }



  ydxn->i1xbjo=TalINT1(ydn->mnopo);
  for(i=0;i<ydn->mnopo;i++)
    {
      ydxn->i1xbjo[i]=-1;
    }

  ydxn->ncmno=ydn->nnopo;


  ydxn->d2xncc=TalDBL2(ydn->mnodim,ydn->mnopo);
  ydxn->d2xnvc=TalDBL2(ydn->mnodim,ydn->mnopo);

  ydxn->d1shap = TalDBL1 (ydn->mnopo); //! -ao this was ydx->ncmno
  ydxn->d1xntc = TalDBL1 (ydn->mnopo);

	for(j=0;j<ydn->mnopo;j++)
	{
		for(i=0;i<ydn->mnodim;i++)
		{
	  		ydxn->d2xncc[i][j]=ydn->d2ncc[i][j];
	  		ydxn->d2xnvc[i][j]=ydn->d2nvc[i][j];
	  	}

		ydxn->d1shap[j]=dthick;
		ydxn->d1xntc[j] = ydn->d1ntc[j];
	}
}



//-PY changed it for 3D_fracture_coupling_with_multiphase
/* read data nodes       
 */
static void Yrsp(YSP ysp, FILE *finp, CHR *name, INT ntpar)
{ 
	CHR *namep;
    INT mnodim;
	mnodim=4;
	namep=name+8;
	if(CHRcmp(namep,"MSPD",4)==0) 
	{  
		INTr(finp,&(ysp->mspd));
		ysp->mspd=ysp->mspd+ntpar;
	}
	else if(CHRcmp(namep,"NSPD",4)==0) 
	{  
		INTr(finp,&(ysp->nspd));
	}
	else if(CHRcmp(namep,"D1MASS",6)==0) 
	{ 
		TformDBL1(finp,R0,ysp->mspd,&(ysp->d1mass));
	}
	else if(CHRcmp(namep,"D2PP",4)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2pp)); 
	}
        /*
	else if(CHRcmp(namep,"D2POLDP",7)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2poldp)); 
	}
        */
	else if(CHRcmp(namep,"D2PVT",5)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2pvt)); 
	}
	
	else if(CHRcmp(namep,"D2PA",4)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2pa)); 
	}
	else if(CHRcmp(namep,"D2PVR",5)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2pvr)); 
	}

	else if(CHRcmp(namep,"D2PRINE",7)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2prine)); 
	}
	else if(CHRcmp(namep,"D2PRINN1",8)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2prinn1)); 
	}
	else if(CHRcmp(namep,"D2PRINN2",8)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2prinn2)); 
	}
	else if(CHRcmp(namep,"D2PRINN3",8)==0) 
	{ 
		TformDBL2(finp,R0,mnodim,ysp->mspd,&(ysp->d2prinn3)); 
	}
	else
	{ 
		CHRw(stderr,"Yrsp: unknown name: ");
		CHRw(stderr,name); 
		CHRwcr(stderr);
		exit(105);
	} 
}



/* allocate memory for variable of super particle  */
  void Yspini(YSP ysp,INT ncindex,INT ntpar,INT icon)
  {INT i,j,mnodim;
  mnodim=4;
  if (ysp->mspd==0) 
  {

  ysp->mspd=ntpar+ncindex+1;

  printf ("ncindex:   %d\n" , ncindex );

//-PY changed it for 3D_fracture_coupling_with_multiphase
//-PY there are some thing wrong with the packing   ncindex is wrong
  //ysp->nspd=ncindex; /*maximum (actual) number of particle*/

  ysp->nspd=1;


  ysp->d1mass=TalDBL1(ysp->mspd);
  ysp->d2pp=TalDBL2(mnodim,ysp->mspd);    /*[nnodim][msd] sphere_position_x                     */
  ysp->d2poldp=TalDBL2(mnodim,ysp->mspd);    /*[nnodim][msd] sphere_position_x                     */
  ysp->d2pvt=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_velocity_x              */
  //ysp->d2pfn=TalDBL2(3,ysp->mspd);    /*[nnodim][msd] sphere_normal force_x         */
  //ysp->d2pft1=TalDBL2(3,ysp->mspd);   /*[nnodim][msd] sphere_tangential force_x1    */
  //ysp->d2pft2=TalDBL2(3,ysp->mspd);   /*[nnodim][msd] sphere_tangential force_x2    */

  ysp->d2ptfc=TalDBL2(mnodim,ysp->mspd);    /*[nnodim][msd] sphere_total_contact force_x */
  ysp->d2ptf=TalDBL2(mnodim,ysp->mspd);    /*[nnodim][msd] sphere_total_force_x */
  ysp->d2pa=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_rotation angel_x      */
  ysp->d2pvr=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_rotation velocity_x      */
  ysp->d2pm=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_moment_x */

  ysp->d2prine=TalDBL2(mnodim,ysp->mspd);  /*[nnodim] pricinple inertia                         */
  ysp->d2prinn1=TalDBL2(mnodim,ysp->mspd);  /*[nnodim] pricinple axe                         */
  ysp->d2prinn2=TalDBL2(mnodim,ysp->mspd);  /*[nnodim] pricinple axe                         */
  ysp->d2prinn3=TalDBL2(mnodim,ysp->mspd);  /*[nnodim] pricinple axe                         */
  ysp->d2poldvr=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_rotation velocity_x    */
  ysp->i1sppr=TalINT1(ysp->mspd);
  ysp->i1con=TalINT1(ysp->mspd);
  ysp->d2pfd=TalDBL2(mnodim,ysp->mspd); 
  ysp->d2pmd=TalDBL2(mnodim,ysp->mspd); 

  for(i=0;i<ysp->mspd;i++)
  {ysp->d1mass[i]=R0;
   ysp->i1sppr[i]=-1;
   if(icon==1){ysp->i1con[i]=0;}
   else {ysp->i1con[i]=1;}

	  for(j=0;j<mnodim-1;j++)
	  {ysp->d2pp[j][i]=R0;
	  ysp->d2poldp[j][i]=R0;
	  ysp->d2pvt[j][i]=R0;
	  //ysp->d2pfn
	  //ysp->d2pft1
      //ysp->d2pft2
	  ysp->d2ptfc[j][i]=R0;
	  ysp->d2ptf[j][i]=R0;
	  ysp->d2pa[j][i]=R0;
	  ysp->d2pvr[j][i]=R0;
	  ysp->d2pm[j][i]=R0;
	  ysp->d2prine[j][i]=R0;
	  ysp->d2prinn1[j][i]=R0;
	  ysp->d2prinn2[j][i]=R0;
	  ysp->d2prinn3[j][i]=R0;
          ysp->d2pfd[j][i]=R0;
          ysp->d2pmd[j][i]=R0;

	  }
  }
  }
	  else
	  {
		  
		  ysp->d2ptfc=TalDBL2(mnodim,ysp->mspd);    /*[nnodim][msd] sphere_total_contact force_x */
		  ysp->d2ptf=TalDBL2(mnodim,ysp->mspd);    /*[nnodim][msd] sphere_total_force_x */
		  ysp->d2pm=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_moment_x */
          ysp->d2poldp=TalDBL2(mnodim,ysp->mspd);    /*[nnodim][msd] sphere_position_x                     */
		  ysp->i1sppr=TalINT1(ysp->mspd);
		  ysp->i1con=TalINT1(ysp->mspd);
                  ysp->d2poldvr=TalDBL2(mnodim,ysp->mspd);     /*[nnodim][msd] sphere_rotation velocity_x    */			  
                 ysp->d2pfd=TalDBL2(mnodim,ysp->mspd);
                 ysp->d2pmd=TalDBL2(mnodim,ysp->mspd);

		  
		  for(i=0;i<ysp->mspd;i++)
		  {
			  ysp->i1sppr[i]=-1;
			  if(icon==1){ysp->i1con[i]=0;}
			  else {ysp->i1con[i]=1;}
			  
			  for(j=0;j<mnodim-1;j++)
			  {
				  ysp->d2ptfc[j][i]=R0;
				  ysp->d2ptf[j][i]=R0;
				  ysp->d2pm[j][i]=R0;
                                  ysp->d2pfd[j][i]=R0;
                                  ysp->d2pmd[j][i]=R0;

				  
			  }
		  }  
	  }
}

/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/


/*! \brief read input
 *  \param[in]     namep file name of input
 *  \param[in,out] yd Y database
 *  \par Details:
 *  Yrd() ...
 *
 */

//-PY changed it for 3D_fracture_coupling_with_multiphase
//INT Yrd(CHR *namep, YD yd)
INT Yrd( CHR *namep, CHR *name_grid,CHR *name_mesh, YD yd)
{ 
  INT icount;
  CHR name[300];
  YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
  YDO ydo=&(yd->ydo);
  YDP ydp=&(yd->ydp);
  YDB ydb=&(yd->ydb);


//-PY changed it for 3D_fracture_coupling_with_multiphase
  YPAR ypar=&(yd->ypar);
  YSP ysp=&(yd->ysp);


  /* added by LGuo */
  YDX ydx=&(yd->ydx);
  YDXN ydxn=&(yd->ydxn);
 

//-PY changed it for 3D_fracture_coupling_with_multiphase
  FILE *fgrid, *fgid;
  INT i,j,index,nelemb,**i2eltob,ntyp,icon;
  DBL dis,tmp,**points;
  INT ntelem,ntnode,*i1nobf,k,m,n,index_sp;
  const double mult[10]={1./6.,1./24.,1./24.,1./24.,1./60.,1./60.,1./60.,1./120.,1./120.,1./120.};
  DBL *v1,*v2,*v3;
  DBL **inertia,ro,cm[3],*intg;
  nelemb=0;

  if((ydc->finp)==(FILENULL))
  { 
    Yrdd(yd);

    ydc->finp=fopen(namep,"r");
    ydc->fcheck=fopen("Ytmp","w");
  }
  if(((ydc->finp)==(FILENULL))||((ydc->fcheck)==(FILENULL)))
  { 
    CHRw(stderr,"Yrd: Could not open input file - usage -i inputfile"); 
    CHRwcr(stderr);   
    return 0;
  }






/* read grid file to define particle position also include total particle number,
 initial particle number, and minimum height distance */
ypar->d1dst=TalDBL1(3);
ypar->d1lv=TalDBL1(3);
ypar->d1vc=TalDBL1(3);

v1=TalDBL1(3);
v2=TalDBL1(3);
v3=TalDBL1(3);
inertia=TalDBL2(3,3);
points=TalDBL2(4,3);
intg=TalDBL1(10);

     fgrid=fopen(name_grid,"r");
    fgid=fopen(name_mesh,"r");

  while(FILEND(fgrid)==0)
  { 


    if(CHRcmp(name,"total_number",12)==0)    
    {    
		INTr(fgrid,&ypar->ntpar); 
	}
	else if(CHRcmp(name,"particle_number_per_layer",15)==0)    
    {    
		INTr(fgrid,&ypar->nlpar); 
	}
	else if(CHRcmp(name,"grid_number",11)==0)    
    {    
		INTr(fgrid,&ypar->ngrid); 
    }
		else if(CHRcmp(name,"maximum_dimension",11)==0)    
    {    
		DBLr(fgrid,&ypar->d1dst[0]); 
		DBLr(fgrid,&ypar->d1dst[1]); 
		DBLr(fgrid,&ypar->d1dst[2]); 
    }
	else if(CHRcmp(name,"normal_direction",11)==0)    
    {    
		DBLr(fgrid,&ypar->d1lv[0]); 
		DBLr(fgrid,&ypar->d1lv[1]); 
		DBLr(fgrid,&ypar->d1lv[2]); 
    }
   else if(CHRcmp(name,"particle_velocity",11)==0)    
    {    
		DBLr(fgrid,&ypar->d1vc[0]); 
		DBLr(fgrid,&ypar->d1vc[1]); 
		DBLr(fgrid,&ypar->d1vc[2]); 
    }
	   else if(CHRcmp(name,"control_velocity",11)==0)    
    {   icon=0; 
		CHRr(fgrid,name);
		if(CHRcmp(name,"yes",3)==0)
		{
			icon=1;
		}
    }
		CHRr(fgrid,name);
  }

/*
	else if(CHRcmp(name,"distance_to_bottom",11)==0)    
    {    
 		DBLr(fgrid,&ypar->d1dxyz[0]);
		DBLr(fgrid,&ypar->d1dxyz[1]);
		DBLr(fgrid,&ypar->d1dxyz[2]);
  }*/
  fclose(fgrid);

  V3DNor(tmp,ypar->d1lv[0],ypar->d1lv[1],ypar->d1lv[2]);

  ypar->d2gxyz=TalDBL2(6,ypar->ngrid);

  fgrid=fopen(name_grid,"r");


  while(FILEND(fgrid)==0)
  { 
	if(CHRcmp(name,"grid_number",11)==0)    
    {    
		CHRr(fgrid,name); 

		for(i=0;i<ypar->ngrid;i++)
		{
			for(j=0;j<6;j++)
			{
			DBLr(fgrid,&ypar->d2gxyz[j][i]);
		}
		}
	}
	CHRr(fgrid,name);
  }



/* read mesh of particle  */

CHRr(fgid,name);
index=0;

while(FILEND(fgid)==0) 
{ if(CHRcmp(name, "Tetrahedra",10)==0)
  {
   	  CHRr(fgid,name); 
	  INTr(fgid,&ntyp);
  }
else if(CHRcmp(name, "Coordinates",11)==0)
{    
	CHRr(fgid,name); 
while(CHRcmp(name, "end",3)!=0)
{
 ypar->nnode++;
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
}

}
else if(CHRcmp(name, "Elements",8)==0)
{    if(index==0)
{index=index+1;
	CHRr(fgid,name); 
while(CHRcmp(name, "end",3)!=0)
{
 ypar->nelem++;
 if(ntyp==4){
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 }
 else if(ntyp==10)
 { CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 

 }
}
}
else 
{
	CHRr(fgid,name); 
while(CHRcmp(name, "end",3)!=0)
{
 nelemb++;
 if(ntyp==4){
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
  CHRr(fgid,name); 

  }
 else if(ntyp==10)
 { CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
 CHRr(fgid,name); 
  CHRr(fgid,name); 

 }
}
}
}
    CHRr(fgid,name); 

}
fclose(fgid);


  fgid=fopen(name_mesh,"r");

ypar->i2elto=TalINT2(ypar->nelem,ntyp);
i2eltob=TalINT2(nelemb,ntyp/2+1);
i1nobf=TalINT1(ypar->nnode);
ypar->d2xyz=TalDBL2(3,ypar->nnode);
ypar->i1elbe=TalINT1(ypar->nelem);


index=0;
while(FILEND(fgid)==0) 
{ 
	
if(CHRcmp(name, "Elements",8)==0)
{
	if(index==0)
   {
	index=index+1;
	for(i=0;i<ypar->nelem;i++)
	{    
 CHRr(fgid,name); 
 for(j=0;j<ntyp;j++)
 {
 INTr(fgid,&ypar->i2elto[i][j]);
 ypar->i2elto[i][j]=ypar->i2elto[i][j]-1;
 }
  CHRr(fgid,name); 

 }
}
else
{	for(i=0;i<nelemb;i++)
	{    
 CHRr(fgid,name); 
 for(j=0;j<(ntyp/2+1);j++)
 {
 INTr(fgid,&i2eltob[i][j]);
 i2eltob[i][j]=i2eltob[i][j]-1;
 }
  CHRr(fgid,name); 

 }
}
}
else if(CHRcmp(name, "Coordinates",11)==0)
{    	if(index==0)
   {

		for(i=0;i<ypar->nnode;i++)
		{
			  CHRr(fgid,name); 

		for(j=0;j<3;j++)
		{
			DBLr(fgid,&ypar->d2xyz[j][i]);
		}
		}
}
}

    CHRr(fgid,name); 
}



/*assign layer number as facenode  */
//if(index2==0)
//{

for(i=0;i<(ypar->nnode);i++)
{ 
	i1nobf[i]=0;   //mark 1 on boundary node 

}

for(i=0;i<(nelemb);i++)
{ for(j=0;j<(ntyp/2+1);j++)
{ 
	i1nobf[i2eltob[i][j]]=1;   //mark 1 on boundary node 
} 
}

for(i=0;i<(ypar->nelem);i++)
{
  for(j=0;j<ntyp;j++)
  {
	  if(i1nobf[ypar->i2elto[i][j]]>0)
		 ypar->i1elbe[i]=1;
  }

}

ypar->d1cg=TalDBL1(3);
    


ypar->d1cg[0]=R0;
ypar->d1cg[1]=R0;
ypar->d1cg[2]=R0;
for(i=0;i<(ypar->nnode);i++)
{
	 	for(j=0;j<3;j++)
		{
	     ypar->d1cg[j]=ypar->d1cg[j]+ypar->d2xyz[j][i];
		}
}

ypar->d1cg[0]=ypar->d1cg[0]/ypar->nnode;
ypar->d1cg[1]=ypar->d1cg[1]/ypar->nnode;
ypar->d1cg[2]=ypar->d1cg[2]/ypar->nnode;


for(i=0;i<(ypar->nnode);i++)
{dis=R0;
	 	for(j=0;j<3;j++)
		{
	     dis=dis+(ypar->d1cg[j]-ypar->d2xyz[j][i])*(ypar->d1cg[j]-ypar->d2xyz[j][i]);
		}
		dis=SQRT(dis);
		if(dis>ypar->radius)ypar->radius=dis;
}

/*generate random number for postion and rotate angle  */

  ypar->max=BEPSILON;
  //ypar->ymax=EPSILON;
  //ypar->zmax=EPSILON;
  
  for(i=0;i<ypar->ngrid;i++)
  {
	  tmp=ypar->d2gxyz[0][i]*ypar->d1lv[0]+ypar->d2gxyz[1][i]*ypar->d1lv[1]+ypar->d2gxyz[2][i]*ypar->d1lv[2];
	  if(tmp<ypar->max)ypar->max=tmp;
	  //if(ypar->d2xyz[1][i]<ypar->xmax)ypar->xmax=ypar->d2xyz[1][i];
	  //if(ypar->d2xyz[2][i]<ypar->xmax)ypar->xmax=ypar->d2xyz[2][i];
	  
  }

ntelem=ypar->ntpar*ypar->nelem;
ntnode=ypar->ntpar*ypar->nnode;

/* Read container mesh or .Y3D */


    SETLINEBUF(ydc->fcheck);  
    CHRr(ydc->finp,name);
    while(FILEND(ydc->finp)==0) 
    { 
        if(CHRcmp(name, "$YSTOP",6)==0)
        { 
            CHRw(ydc->fcheck,name); CHRwcr(ydc->fcheck);
            return 0;
        }
        else if(CHRcmp(name, "$YDOIT",6)==0)
        {
            /* added by LGuo */

/* added by LGuo */
      yde->i2elfr=TalINT2(yde->melem+1,4);
      for(i=0;i<yde->melem;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      yde->i2elfr[i][j]=-1;
	    }
	}
            //-PY changed it for 3D_fracture_coupling_with_multiphase  

            index_sp=0;


//            Yrdx(ydx, yde, ydc->finp, name, ydn->nnopo);
            Yrdxn(ydxn,ydn,ydx->dthick);
            CHRw(ydc->fcheck,name); CHRwcr(ydc->fcheck);
if(ydc->irigid==1)
{
            //-PY changed it for 3D_fracture_coupling_with_multiphase
            for(i=0;i<yde->nelem;i++)
            {
               //-PY the problem is  yde->i2elto
	        if(yde->i2elto[i][yde->nelno-1]>ypar->nindex) 
                {         
                ypar->nindex=yde->i2elto[i][yde->nelno-1];
                }

		for(j=0;j<(yde->nelno-1);j++)
		{
	            ydn->i1nind[yde->i2elto[i][j]]=yde->i2elto[i][yde->nelno-1];
		}
            }

            if(ysp->mspd==0) index_sp=1;
            Yspini(ysp,ypar->nindex,ypar->ntpar,icon);

//-PY switch it on for rigid solid
        for(i=0;i<ysp->nspd;i++)
	    {
                printf ("ysp->nspd:   %d\n" , ysp->nspd );


//                 to be modified properly
	        for(n=0;n<ydn->nnopo;n++)
    	        {   
		    ysp->i1con[i]=1;
		    if(ydn->i1nind[n]==(i+1))
	            {
	                ysp->d2pvt[0][i]=ydn->d2nvc[0][n];
			ysp->d2pvt[1][i]=ydn->d2nvc[1][n];
			ysp->d2pvt[2][i]=ydn->d2nvc[2][n];
			if(ydn->i1nopr[n]!=(-1))ysp->i1sppr[i]=ydn->i1nopr[n];
                    } 
	        }
		 
//		   end 
		  
	        if (index_sp==1)
	        {
		    for(n=0;n<10;n++) intg[n]=R0;
	     
	  	    for(k=0;k<yde->nelem;k++)
		    {	 
                        if(yde->i2elto[k][yde->nelno-1]==(i+1))
		        {
                            for(m=0;m<4;m++)
		            {
	                        points[m][0]=ydn->d2ncc[0][yde->i2elto[k][m]];
		                points[m][1]=ydn->d2ncc[1][yde->i2elto[k][m]];
		                points[m][2]=ydn->d2ncc[2][yde->i2elto[k][m]];
		            }
		            ro=ydp->d1pero[yde->i1elpr[k]];
                            mirtichRoutine(points,intg,ro);
	                }
                    }

                    for(n=0;n<10;n++)  
                         intg[n]=intg[n]*mult[n];
                    //this is an error in the original publication!!!!!!!!!!!!!!!!!
                    //they dont have the three!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ysp->d1mass[i] = intg[0];  //density = 1.0 assumed 

                    //center of mass 
                    ysp->d2pp[0][i]=intg[1]/ysp->d1mass[i];
                    ysp->d2pp[1][i]=intg[2]/ysp->d1mass[i];
                    ysp->d2pp[2][i]=intg[3]/ysp->d1mass[i];
	            cm[0]=ysp->d2pp[0][i];
	            cm[1]=ysp->d2pp[1][i];
	            cm[2]=ysp->d2pp[2][i];
    
                    inertia[0][0] = intg[5] + intg[6] -ysp->d1mass[i]*(cm[1]*cm[1]+cm[2]*cm[2]);
                    inertia[1][1] = intg[4] + intg[6] -ysp->d1mass[i]*(cm[0]*cm[0]+cm[2]*cm[2]);
                    inertia[2][2] = intg[4] + intg[5] -ysp->d1mass[i]*(cm[1]*cm[1]+cm[0]*cm[0]);

                    inertia[0][1] = -(intg[7] -ysp->d1mass[i]*cm[1]*cm[0]);
                    inertia[1][0]=inertia[0][1];

                    inertia[1][2] = -(intg[8] -ysp->d1mass[i]*cm[1]*cm[2]);
                    inertia[2][1]=inertia[1][2];

                    inertia[0][2] = -(intg[9] -ysp->d1mass[i]*cm[2]*cm[0]);
                    inertia[2][0]=inertia[0][2];

                    //    inertia*=((REAL)(1.0/3.0));


	            for(m=0;m<3;m++)
	            {
	                v1[m]=R0;
                        v2[m]=R0;
	                v3[m]=R0;


/*
                        for(n=0;n<3;n++)
	                {
                            inertia[m][n]=inertia[m][n]*(R1/R3);
	                }
*/
	            }
  
                    //this is because the vertices could have been DoubleLinkedListed in the wrong direction
                    if((inertia[0][0]<0.0)||(inertia[1][1]<0.0)||(inertia[2][2]<0.0))
                    {
                        //inertia*=(-1.0);
	            for(m=0;m<3;m++)
	            {
		        for(n=0;n<3;n++)
	                {
                            inertia[m][n]=-inertia[m][n];
	                }
	            }
                    ysp->d1mass[i]*=(-1.0);
                }
     //   CHRw(stdout, "output-14-1\n");

     printf("ysp->d2prine[0][i]: %f\n",ysp->d2prine[0][i]);
     printf("ysp->d2prine[1][i]: %f\n",ysp->d2prine[1][i]);
     printf("ysp->d2prine[2][i]: %f\n",ysp->d2prine[2][i]);


                //if(i==1){
	        solveSymetricEigenProblem(inertia, v1,v2,v3,
						&ysp->d2prine[0][i],&ysp->d2prine[1][i],&ysp->d2prine[2][i]);
                //}

   //  CHRw(stdout, "output-14-1-1\n");
		for(m=0;m<3;m++)
	        {
		    ysp->d2prinn1[m][i]=v1[m];
		    ysp->d2prinn2[m][i]=v2[m];
		    ysp->d2prinn3[m][i]=v3[m];
		}
       
    //  CHRw(stdout, "output-14-1-2\n");
				
	    }
        }
//-PY switch it on for rigid solid

//-ao 050517


	/*check condition  */
	index=0;
	for(i=0;i<ydb->nbcon;i++)
	{
            if((ydb->i1bnvx[i]==0)&&(ydb->i1bnvy[i]==0)&&(ydb->i1bnvz[i]==0))
            {
                index=i+1;
            }
        }



	if(index==0)
	{   //ydb->nbcon+=1;
            ydb->d1bnvx[ydb->nbcon]=ydb->d1bnvx[0];
            ydb->d1bnvy[ydb->nbcon]=ydb->d1bnvy[0];
            ydb->d1bnvz[ydb->nbcon]=ydb->d1bnvz[0];
            ydb->d1bnax[ydb->nbcon]=ydb->d1bnax[0];
            ydb->d1bnay[ydb->nbcon]=ydb->d1bnay[0];
            ydb->d1bnaz[ydb->nbcon]=ydb->d1bnaz[0];
            ydb->d1bnfx[ydb->nbcon]=ydb->d1bnfx[0];
            ydb->d1bnfy[ydb->nbcon]=ydb->d1bnfy[0];
            ydb->d1bnfz[ydb->nbcon]=ydb->d1bnfz[0];
            ydb->i1bnvx[ydb->nbcon]=ydb->i1bnvx[0];
            ydb->i1bnvy[ydb->nbcon]=ydb->i1bnvy[0];
            ydb->i1bnvz[ydb->nbcon]=ydb->i1bnvz[0];

            ydb->i1bnvx[0]=0;
            ydb->i1bnvy[0]=0;
            ydb->i1bnvz[0]=0;

            for(i=0;i<ydn->nnopo;i++)
            {
                if(ydn->i1nopr[i]==0)ydn->i1nopr[i]=ydb->nbcon;
            }
            for(i=0;i<ydn->mnopo;i++)
            {
                if(ydn->i1nopr[i]==(-1))ydn->i1nopr[i]=0;
            }

            for(i=0;i<ysp->mspd;i++)
            {
                if(ysp->i1sppr[i]==0) ysp->i1sppr[i]=ydb->nbcon;
            }

            for(i=0;i<ysp->mspd;i++)
            {
                if(ysp->i1sppr[i]==(-1)) ysp->i1sppr[i]=0;
            }

                ydb->nbcon+=1;

        }


        else
        {
            for(i=0;i<ydn->mnopo;i++)
            {
            if(ydn->i1nopr[i]==(-1))ydn->i1nopr[i]=index-1;
            }
            for(i=0;i<ysp->mspd;i++)
            {
            if(ysp->i1sppr[i]==(-1)) ysp->i1sppr[i]=index-1;
            }
        }
}
        FREE(i2eltob);
        FREE(i1nobf);
        FREE(v1);
	    FREE(v2);
	    FREE(v3);
        FREE(inertia);
	    FREE(points);
	    FREE(intg);


        return 1;
     }

    else if(CHRcmp(name,"/*",2)==0)    /* read and ignore comments */
    { 
      icount=0;
      do
      { 
        CHRr(ydc->finp,name); icount++;
        if(icount>100)
        { 
          CHRw(stderr,"Yrd: too long comment near - ");
          CHRw(stderr,name);
          CHRwcr(stderr);      
          return 0;
        } 
      } while((FILEND(ydc->finp)==0)&&(CHRcmp(name,"*/",2)!=0));
    }
    else if(CHRcmp(name, "/YD/",4)==0)
    { 
      CHRw(ydc->fcheck,name); CHRwcr(ydc->fcheck);
      if(CHRcmp(name,"/YD/YDC/",8)==0)        /* read control data    */
      { 
        Yrdc(ydc,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDE/",8)==0)   /* read data elements   */
      { 

//-PY changed it for 3D_fracture_coupling_with_multiphase
        Yrde(yde,ydc->finp,name,ntelem);


	//Yrdx(ydx,yde);
      }
      else if(CHRcmp(name,"/YD/YDI/",8)==0)  /* read data interaction */
      { 
        Yrdi(ydi,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDN/",8)==0)  /* read data nodes       */
      { 

//-PY changed it for 3D_fracture_coupling_with_multiphase
        Yrdn(ydn,ydc->finp,name,ntnode);



	//Yrdxn(ydxn,ydn);
      }
      else if(CHRcmp(name,"/YD/YDO/",8)==0)  /* read output spec      */
      { 
        Yrdo(ydo,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDP/",8)==0)  /* read data properties  */
      { 
        Yrdp(ydp,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDB/",8)==0)  /* read data properties  */
      { 
		    Yrdb(ydb,ydc->finp,name);
      }
//-PY changed it for 3D_fracture_coupling_with_multiphase
      else if(CHRcmp(name,"/YD/YSP/",8)==0)  /* read data properties  */
      { 
		  Yrsp(ysp,ydc->finp,name,ypar->ntpar);
      }
      else if(CHRcmp(name,"/YD/YDX/",8)==0)  /* surface mesh             */
      { 
          Yrdx(ydx, yde, ydc->finp,name,ydn->mnopo);
      }




      else
      { 
        CHRw(stderr,"Yrd: unknown name: ");
        CHRw(stderr,name); 
        CHRwcr(stderr);   
        return 0;
      } 
    }
    else
    { 
      CHRw(stderr,"Yrd: unknown name: ");
      CHRw(stderr,name); 
      CHRwcr(stderr);
      return 0;
    }
    CHRr(ydc->finp,name); 
  }


  fclose(ydc->finp);
  fclose(ydc->fcheck);
  fclose(fgrid);
  fclose(fgid);
  FREE(i2eltob);
  FREE(i1nobf);

  return 0;
}


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/


