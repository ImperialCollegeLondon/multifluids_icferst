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
#include "Yd.h"
#define N 100


/*******************   Read Input File   ******************************/
static void Yrdc(ydc,finp,name)
  YDC ydc; FILE *finp; CHR *name;
{ CHR *namep;
  
  namep=name+8;
  if(CHRcmp(namep,"MCSTEP",6)==0) 
  { INTr(finp,&(ydc->mcstep));
  }
  else if(CHRcmp(namep,"NCSTEP",6)==0) 
  { INTr(finp,&(ydc->ncstep));
  }
  else if(CHRcmp(namep,  "ISAVE",5)==0) 
  { INTr(finp,&(ydc->isave));
  }
  else if(CHRcmp(namep,"DCGRAX",6)==0) 
  { DBLr(finp,&(ydc->dcgrax));
  }
  else if(CHRcmp(namep,"DCGRAY",6)==0) 
  { DBLr(finp,&(ydc->dcgray));
  }
  else if(CHRcmp(namep,"DCGRAZ",6)==0) 
  { DBLr(finp,&(ydc->dcgraz));
  } 
//will delete DCSIZC, DCSIZF, DCSIZS, DCSIZ 
 else if(CHRcmp(namep,"DCSIZC",6)==0) 
  { DBLr(finp,&(ydc->dcsizc));

  }
  else if(CHRcmp(namep,"DCSIZF",6)==0) 
  { DBLr(finp,&(ydc->dcsizf));
  }
  else if(CHRcmp(namep,"DCSIZS",6)==0) 
  { DBLr(finp,&(ydc->dcsizs));
  }
  else if(CHRcmp(namep,"DCSIZV",6)==0) 
  { DBLr(finp,&(ydc->dcsizv));
  }
  else if(CHRcmp(namep,"DCSTEC",6)==0) 
  { DBLr(finp,&(ydc->dcstec));
  }
  else if(CHRcmp(namep,"DCTIME",6)==0) 
  { DBLr(finp,&(ydc->dctime));
  }
  else if(CHRcmp(namep,"DCRMPT",6)==0) 
  { DBLr(finp,&(ydc->dcrmpt));
  }
  else if(CHRcmp(namep,"DCURELX",7)==0) 
  { DBLr(finp,&(ydc->dcurelx));
  }
  else if(CHRcmp(namep,"INITER",6)==0) 
  { INTr(finp,&(ydc->initer));
  }
  else if(CHRcmp(namep,"ICOUTF",6)==0) 
  { INTr(finp,&(ydc->icoutf));
  }
  else if(CHRcmp(namep,"ICOUTI",6)==0) 
  { INTr(finp,&(ydc->icouti));
  }
  else if(CHRcmp(namep,"ICOUTP",6)==0) 
  { INTr(finp,&(ydc->icoutp));
  }
  else if(CHRcmp(namep,"IWFAST",6)==0) 
  { INTr(finp,&(ydc->iwfast));
  }

  else
  { CHRw(stderr,"Yrdc: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(103);
} }
 
static void Yrdd(yd) /* default values */
   YD yd;
{ YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
  YDO ydo=&(yd->ydo);
  YDP ydp=&(yd->ydp);
  YDB ydb=&(yd->ydb);
  YDX ydx=&(yd->ydx);

  /* Set Control to default   */

  ydc->mcstep=0; ydc->ncstep=0;
  ydc->finp=(FILE*)NULL;
  ydc->fcheck=(FILE*)NULL;
  ydc->dcgrax=R0;
  ydc->dcgray=R0;
  ydc->dcgraz=R0;
  ydc->dcsizc=R1;
  ydc->dcsizf=R1;
  ydc->dcsizv=R1;
  ydc->dcstec=R1;
  ydc->dctime=R0;
  ydc->dcurelx=R1;
  ydc->initer=2;
  ydc->icoutf=0;
  ydc->icouti=0;
  ydc->icoutp=2;
  ydc->iwfast=1;
  ydc->dcrmpt=R0;

/* Set Elements to default    */
  yde->melem=0; yde->nelem=0;
  yde->melst=0; yde->nelst=0;
  yde->melno=0; yde->nelno=0;
  yde->i1elcf=INT1NULL; 
  yde->i1elbe=INT1NULL;
  yde->i1elpr=INT1NULL;
  yde->d2elst=DBL2NULL;
  yde->i2elto=INT2NULL;
  yde->d3tcs=DBL3NULL;
  yde->d1emct=DBL1NULL;
 
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
 
  /* Set Nodes to default  */
  ydn->mnodim=0;  ydn->nnodim=0;
  ydn->mnopo=0;  ydn->nnopo=0;
  ydn->d1nmct=DBL1NULL;
  ydn->d1nti=DBL1NULL;
  ydn->d2ncc=DBL2NULL;
  ydn->d2nci=DBL2NULL;
  ydn->d2nfc=DBL2NULL;
  ydn->d2nft=DBL2NULL;
  ydn->d2nvc=DBL2NULL;
  ydn->i1nobf=INT1NULL;
  ydn->i1nopr=INT1NULL;
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

  /* Set Properties to default  */
  ydp->mprop=0; ydp->nprop=0;
  ydp->d1peks=DBL1NULL;
  ydp->d1pela=DBL1NULL;
  ydp->d1pemu=DBL1NULL;
  ydp->d1pepe=DBL1NULL;
  ydp->d1pero=DBL1NULL;
  ydp->d1pegfn=DBL1NULL;
  ydp->d1pegfs=DBL1NULL;
  ydp->d1peft=DBL1NULL;
  //  ydp->d1pefs=DBL1NULL;

  ydp->d1pcoh=DBL1NULL;
  ydp->d1picf=DBL1NULL;

  ydp->d1pefr=DBL1NULL;
  ydp->i1ptyp=INT1NULL;
  ydp->i1pejp=INT1NULL;
  ydp->i1pemn=INT1NULL;
  ydp->i1psde=INT1NULL;
  ydp->d1tcon=DBL1NULL;
  ydp->d1capa=DBL1NULL;
  ydp->d1ctex=DBL1NULL;
  ydp->d1pesf=DBL1NULL;
  ydp->d1pepsf=DBL1NULL;
  ydp->d1pevf=DBL1NULL;
  ydp->d1pepf=DBL1NULL;

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
  ydb->d1bntp=DBL1NULL;
  ydb->i1bntp=INT1NULL;
  ydb->d1bnhf=DBL1NULL;
  ydb->i1bnhf=INT1NULL;
  ydb->i1bnvx=INT1NULL;
  ydb->i1bnvy=INT1NULL;
  ydb->i1bnvz=INT1NULL;
  ydb->d1bcvt=DBL1NULL;
  ydb->d1bcvc=DBL1NULL;
  ydb->i1bcvt=INT1NULL;

  /* Set Boundary Mesh Condition to default  */
  ydx->nelem=0;
  ydx->melem=0;
  ydx->nelno=0;
  ydx->i2elto=INT2NULL;
}
 
static void Yrde(yde,finp,name)
  YDE    yde; FILE *finp; CHR *name; 
{ CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MELEM",5)==0) 
  {  INTr(finp,&(yde->melem));
  }
  else if(CHRcmp(namep,"NELEM",5)==0) 
  {  INTr(finp,&(yde->nelem));
  }
  else if(CHRcmp(namep,"MELST",5)==0) 
  {  INTr(finp,&(yde->melst));
  }
  else if(CHRcmp(namep,"NELST",5)==0) 
  {  INTr(finp,&(yde->nelst));
  }
  else if(CHRcmp(namep,"MELNO",5)==0) 
  {  INTr(finp,&(yde->melno));
  }
  else if(CHRcmp(namep,"NELNO",5)==0) 
  {  INTr(finp,&(yde->nelno));
  }
  else if(CHRcmp(namep,"I1ELCF",6)==0) 
  { TformINT1(finp,-1,yde->melem,&(yde->i1elcf));
  }
  else if(CHRcmp(namep,"I1ELBE",6)==0) 
  { TformINT1(finp,0,yde->melem,&(yde->i1elbe));
  }
  else if(CHRcmp(namep,"I1ELPR",6)==0) 
  { TformINT1(finp,0,yde->melem,&(yde->i1elpr));
  }
  else if(CHRcmp(namep,"D2ELST",6)==0) 
  { TformDBL2(finp,R0,yde->melst,yde->melem,&(yde->d2elst));
  }
  else if(CHRcmp(namep,"I2ELTO",6)==0) 
  {
	  TformINT2(finp,-1,yde->melno,yde->melem,&(yde->i2elto));
  }
  else if(CHRcmp(namep,"D3TCS",5)==0) 
  { TformDBL3(finp,R0,3,3,yde->melem,&(yde->d3tcs));
  }  
  else if(CHRcmp(namep,"D1EMCT",6)==0) 
  { TformDBL1(finp,R0,yde->melem,&(yde->d1emct));
  }
  else
  { CHRw(stderr,"Yrde: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
} }

static void Yrdi(ydi,finp,name)
  YDI   ydi; FILE *finp; CHR *name; 
{ CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MICOUP",6)==0) 
  {  INTr(finp,&(ydi->micoup));
  }
  else if(CHRcmp(namep,"NICOUP",6)==0) 
  {  INTr(finp,&(ydi->nicoup));
  }
  else if(CHRcmp(namep,"IIECFF",6)==0) 
  {  INTr(finp,&(ydi->iiecff));
  }
  else if(CHRcmp(namep,"DIEDI",5)==0) 
  {  DBLr(finp,&(ydi->diedi));
  }
  else if(CHRcmp(namep,"DIEZON",6)==0) 
  {  DBLr(finp,&(ydi->diezon));
  }
  else if(CHRcmp(namep,"D1IESL",6)==0) 
  { TformDBL1(finp,R0,ydi->micoup,&(ydi->d1iesl)); 
  }
  else if(CHRcmp(namep,"I1IECN",6)==0) 
  { TformINT1(finp,-1,ydi->micoup,&(ydi->i1iecn)); 
  }
  else if(CHRcmp(namep,"I1IECT",6)==0) 
  { TformINT1(finp,-1,ydi->micoup,&(ydi->i1iect)); 
  }
  else if(CHRcmp(namep,"D1DELTAT1",9)==0) 
  { TformDBL1(finp,R0,ydi->micoup,&(ydi->d1deltat1)); 
  }
  else if(CHRcmp(namep,"D1DELTAT2",9)==0) 
  { TformDBL1(finp,R0,ydi->micoup,&(ydi->d1deltat2)); 
  }
  else if(CHRcmp(namep,"D1DELTAN",8)==0) 
  { TformDBL1(finp,R0,ydi->micoup,&(ydi->d1deltan)); 
  }

  else if(CHRcmp(namep,"D2NV",4)==0) 
  { TformDBL2(finp,R0,3,ydi->micoup,&(ydi->d2nv));
  }
  else if(CHRcmp(namep,"D2T1V",5)==0) 
  { TformDBL2(finp,R0,3,ydi->micoup,&(ydi->d2t1v));
  }
  else if(CHRcmp(namep,"D2T2V",5)==0) 
  { TformDBL2(finp,R0,3,ydi->micoup,&(ydi->d2t2v));
  }
  else
  { CHRw(stderr,"Yrdi: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
} }

static void Yrdn(ydn,finp,name)
   YDN    ydn; FILE *finp; CHR *name; 
{ CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MNODIM",5)==0) 
  {  INTr(finp,&(ydn->mnodim));
  }
  else if(CHRcmp(namep,"NNODIM",5)==0) 
  {  INTr(finp,&(ydn->nnodim));
  }
  else if(CHRcmp(namep,"MNOPO",5)==0) 
  {  INTr(finp,&(ydn->mnopo));
  }
  else if(CHRcmp(namep,"NNOPO",5)==0) 
  {  INTr(finp,&(ydn->nnopo));
  }
  else if(CHRcmp(namep,"D1NMCT",6)==0) 
  { TformDBL1(finp,R0,ydn->mnopo,&(ydn->d1nmct));
  }
  else if(CHRcmp(namep,"D1NTI",5)==0) 
  { TformDBL1(finp,R0,ydn->mnopo,&(ydn->d1nti));
  }
  else if(CHRcmp(namep,"D2NCC",5)==0) 
  { TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2ncc)); 
  }
  else if(CHRcmp(namep,"D2NCI",5)==0) 
  { TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nci)); 
  }
  else if(CHRcmp(namep,"D2NFC",5)==0) 
  { TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nfc)); 
  }
  else if(CHRcmp(namep,"D2NFT",5)==0) 
  { TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nft)); 
  }
  else if(CHRcmp(namep,"D2NVC",5)==0) 
  { TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nvc)); 
  }
  else if(CHRcmp(namep,"I1NOBF",6)==0) 
  { TformINT1(finp,0,ydn->mnopo,&(ydn->i1nobf));
  }
  else if(CHRcmp(namep,"I1NOPR",6)==0) 
  { TformINT1(finp,1,ydn->mnopo,&(ydn->i1nopr));
  }
  else
  { CHRw(stderr,"Yrdn: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(105);
} }

static void Yrdo(ydo,finp,name)
   YDO    ydo; FILE *finp; CHR *name; 
{ INT i; CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MOHYS",5)==0) 
  {  INTr(finp,&(ydo->mohys))
     i=(ydo->mohys)*sizeof(FILE*);
     if(i>0)ydo->f2ohyf=(FILE**)MALLOC(i);
     for(i=0;i<(ydo->mohys);i++)
     { ydo->f2ohyf[i]=FILENULL;
     }
  }
  else if(CHRcmp(namep,"NOHYS",5)==0) 
  {  INTr(finp,&(ydo->nohys));
  }
  else if(CHRcmp(namep,"DOHYP",5)==0) 
  { DBLr(finp,&(ydo->dohyp));
  }
  else if(CHRcmp(namep,"D1OHYS",6)==0) 
  { TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohys));
  }
  else if(CHRcmp(namep,"D1OHYC",6)==0) 
  { TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyc));
  }
  else if(CHRcmp(namep,"D1OHYF",6)==0) 
  { TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyf));
  }
  else if(CHRcmp(namep,"D1OHYT",6)==0) 
  { TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyt));
  }
  else if(CHRcmp(namep,"D1OHYX",6)==0) 
  { TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyx));
  }
  else if(CHRcmp(namep,"D1OHYY",6)==0) 
  { TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyy));
  }
  else if(CHRcmp(namep,"D1OHYZ",6)==0) 
  { TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyz));
  }
  else if(CHRcmp(namep,"I1OHYT",6)==0) 
  { TformINT1(finp,-1,ydo->mohys,&(ydo->i1ohyt));
  }
  else
  { CHRw(stderr,"Yrdo: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(105);
} }

static void Yrdp(ydp,finp,name)
  YDP    ydp; FILE *finp; CHR *name; 
{ CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MPROP",5)==0) 
  {  INTr(finp,&(ydp->mprop));
  }
  else if(CHRcmp(namep,"NPROP",5)==0) 
  {  INTr(finp,&(ydp->nprop));
  }
  else if(CHRcmp(namep,"D1PEKS",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1peks)); 
  }
  else if(CHRcmp(namep,"D1PELA",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pela)); 
  }
  else if(CHRcmp(namep,"D1PEMU",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pemu)); 
  }
  else if(CHRcmp(namep,"D1PEPE",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pepe)); 
  }
  else if(CHRcmp(namep,"D1PERO",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pero)); 
  }
  else if(CHRcmp(namep,"D1PEGFS",7)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pegfs));
  }
  else if(CHRcmp(namep,"D1PEGFN",7)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pegfn));
  }
  else if(CHRcmp(namep,"D1PEFT",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1peft)); 
  }
  /*
  else if(CHRcmp(namep,"D1PEFS",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pefs)); 
  }
  */
  else if(CHRcmp(namep,"D1PCOH",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pcoh)); 
  }
  else if(CHRcmp(namep,"D1PICF",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1picf)); 
  }
  else if(CHRcmp(namep,"D1PEFR",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pefr)); 
  }
  else if(CHRcmp(namep,"D1PESF",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pesf)); 
  }
  else if(CHRcmp(namep,"D1PEPSF",7)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pepsf)); 
  }
  else if(CHRcmp(namep,"D1PEVF",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pevf)); 
  }
  else if(CHRcmp(namep,"D1PEPF",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pepf)); 
  }

  else if(CHRcmp(namep,"I1PTYP",6)==0) 
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1ptyp));
  }

  else if(CHRcmp(namep,"I1PEJP",6)==0) 
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1pejp));
  }
  else if(CHRcmp(namep,"I1PEMN",6)==0) 
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1pemn));
  }
  else if(CHRcmp(namep,"I1PSDE",6)==0) 
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1psde));
  }
  else if(CHRcmp(namep,"D1TCON",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1tcon)); 
  }
  else if(CHRcmp(namep,"D1CAPA",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1capa)); 
  }  
  else if(CHRcmp(namep,"D1CTEX",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1ctex));
  }
  else
  { CHRw(stderr,"Yrdp: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
} }


static void Yrdb(ydb,finp,name)
  YDB    ydb; FILE *finp; CHR *name; 
{ CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MBCON",5)==0) 
  {  INTr(finp,&(ydb->mbcon));
  }
  else if(CHRcmp(namep,"NBCON",5)==0) 
  {  INTr(finp,&(ydb->nbcon));
  }
  else if(CHRcmp(namep,"D1BNAX",6)==0) 
  { TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnax)); 
  }
  else if(CHRcmp(namep,"D1BNAY",6)==0) 
  { TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnay)); 
  }
  else if(CHRcmp(namep,"D1BNAZ",6)==0) 
  { TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnaz)); 
  }
  else if(CHRcmp(namep,"D1BNFX",6)==0) 
  { TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnfx));
  }
  else if(CHRcmp(namep,"D1BNFY",6)==0) 
  { TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnfy));
  }
  else if(CHRcmp(namep,"D1BNFZ",6)==0) 
  { TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnfz));
  }
  else if(CHRcmp(namep,"I1BNVX",6)==0) 
  { TformINT1(finp,0,ydb->mbcon,&(ydb->i1bnvx));
  }
  else if(CHRcmp(namep,"I1BNVY",6)==0) 
  { TformINT1(finp,0,ydb->mbcon,&(ydb->i1bnvy));
  }
  else if(CHRcmp(namep,"I1BNVZ",6)==0) 
  { TformINT1(finp,0,ydb->mbcon,&(ydb->i1bnvz));
  }
  else if(CHRcmp(namep,"D1BNVX",6)==0) 
  { TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnvx));
  }
  else if(CHRcmp(namep,"D1BNVY",6)==0) 
  { TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnvy));
  }
  else if(CHRcmp(namep,"D1BNVZ",6)==0) 
  { TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnvz));
  }
  else if(CHRcmp(namep,"D1BNTP",6)==0)
  { TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bntp));
  }
  else if(CHRcmp(namep,"D1BNHF",6)==0) 
  { TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bnhf));
  }
  else if(CHRcmp(namep,"I1BNTP",6)==0) 
  { TformINT1(finp,0,ydb->mbcon,&(ydb->i1bntp));
  }
  else if(CHRcmp(namep,"I1BNHF",6)==0) 
  { TformINT1(finp,0,ydb->mbcon,&(ydb->i1bnhf));
  }
  else if(CHRcmp(namep,"D1BCVT",6)==0)
  { TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bcvt));
  }
  else if(CHRcmp(namep,"D1BCVC",6)==0)
  { TformDBL1(finp,R0,ydb->mbcon,&(ydb->d1bcvc));
  }
  else if(CHRcmp(namep,"I1BCVT",6)==0) 
  { TformINT1(finp,0,ydb->mbcon,&(ydb->i1bcvt));
  }
  else
  { CHRw(stderr,"Yrdp: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
} }

//static void Yrdx(YDX ydx,YDP ydp, YDE yde, FILE *finp, CHR *name)
//{
//    CHR *namep; //INT i,j;
// 
//    namep=name+8;
//	if(CHRcmp(namep,"MELEM",5)==0)
//    {
//        INTr(finp,&(ydx->melem));
//        ydx->melem=ydx->melem+1;
//    }
//    else if(CHRcmp(namep,"NELEM",5)==0)
//    {
//        INTr(finp,&(ydx->nelem));
//        //ydx->nelem=ydx->nelem+ntelem;
//    }
//    else if(CHRcmp(namep,"NELNO",5)==0)
//    {
//        INTr(finp,&(ydx->nelno));
//    }
//    else if(CHRcmp(namep,"i2elto",7)==0)
//    {        
//        //Z transformed array order
//        TformINT2(finp,-1,yde->melno,yde->melem,&(yde->i2elto));
//     	TformINT2_inv(finp,-1,ydx->nelem+1,ydx->nelno+1,&(ydx->i2elto));
//     	TformINT2(finp,-1,ydx->nelno,ydx->nelem+1,&(ydx->i2elto)); 
//    }
//    else
//    { 
//        CHRw(stderr,"Yrdx: unknown name: ");
//        CHRw(stderr,name); 
//        CHRwcr(stderr);
//        exit(106);
//    } 
//}




/****************** Control ***************************************/
static void Ywdc(ydc,fout)
  YDC ydc; FILE *fout;
{ CHRw(fout,"     /*   Control     */");
  CHRwcr(fout);
  CHRw(fout,"/YD/YDC/MCSTEP");
  CHRwsp(fout);
  INTw(fout,ydc->mcstep,8);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/NCSTEP");
  CHRwsp(fout);
  INTw(fout,ydc->ncstep,8);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/ISAVE");
  CHRwsp(fout);
  INTw(fout,ydc->isave,8);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCGRAX");
  CHRwsp(fout);
  DBLw(fout,ydc->dcgrax,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCGRAY");
  CHRwsp(fout);
  DBLw(fout,ydc->dcgray,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCGRAZ");
  CHRwsp(fout);
  DBLw(fout,ydc->dcgraz,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCSIZC");
  CHRwsp(fout);
  DBLw(fout,ydc->dcsizc,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCSIZF");
  CHRwsp(fout);
  DBLw(fout,ydc->dcsizf,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCSIZS");
  CHRwsp(fout);
  DBLw(fout,ydc->dcsizs,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCSIZV");
  CHRwsp(fout);
  DBLw(fout,ydc->dcsizv,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCSTEC");
  CHRwsp(fout);
  DBLw(fout,ydc->dcstec,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCTIME");
  CHRwsp(fout);
  DBLw(fout,ydc->dctime,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCRMPT");
  CHRwsp(fout);
  DBLw(fout,ydc->dcrmpt,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCURELX");
  CHRwsp(fout);
  DBLw(fout,ydc->dcurelx,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/INITER");
  CHRwsp(fout);
  INTw(fout,ydc->initer,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/ICOUTF");
  CHRwsp(fout);
  INTw(fout,ydc->icoutf,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/ICOUTI");
  CHRwsp(fout);
  INTw(fout,ydc->icouti,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/ICOUTP");
  CHRwsp(fout);
  INTw(fout,ydc->icoutp,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/IWFAST");
  CHRwsp(fout);
  INTw(fout,ydc->iwfast,5);
  CHRwcr(fout);
  CHRwcr(fout);

}	

/****************** Elements ***************************************/

static void Ywde(yde,fout,i2elbnp)
  YDE yde; FILE *fout; INT **i2elbnp;
{ INT icount;
  INT jcount;
  
  CHRw(fout,"     /*   Elements     */");
  CHRwcr(fout);
  CHRw(fout,"/YD/YDE/MELEM");
  CHRwsp(fout);
  INTw(fout,yde->melem,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/NELEM");
  CHRwsp(fout);
  INTw(fout,yde->nelem,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/MELST");
  CHRwsp(fout);
  INTw(fout,yde->melst,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/NELST");
  CHRwsp(fout);
  INTw(fout,yde->nelst,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/MELNO");
  CHRwsp(fout);
  INTw(fout,yde->melno,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/NELNO");
  CHRwsp(fout);
  INTw(fout,yde->nelno,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/D2ELST");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,yde->nelem,5);
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/I1ELCF");
  CHRwsp(fout);
  INTw(fout,yde->nelem,5);
  CHRwcr(fout);
  for(icount=0;icount<yde->nelem;icount++)
  { CHRw(fout,"  -1");
    CHRwsp(fout);
  } 
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/I1ELBE");
  CHRwsp(fout);
  INTw(fout,yde->nelem,5);
  CHRwcr(fout);
  for(icount=0;icount<yde->nelem;icount++)
  { INTw(fout,yde->i1elbe[icount],4);
    CHRwsp(fout);
  } 
  CHRwcr(fout);


  CHRw(fout,"/YD/YDE/I1ELPR");
  CHRwsp(fout);
  INTw(fout,yde->nelem,5);
  CHRwcr(fout);

  for(icount=0;icount<yde->nelem;icount++)
  { INTw(fout,yde->i1elpr[icount],4);
    CHRwsp(fout);
  }
  CHRwcr(fout);
 
  CHRw(fout,"/YD/YDE/I2ELTO");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,yde->nelno,5);
  CHRwsp(fout);
  INTw(fout,yde->nelem,5);
  CHRwcr(fout);

  for(icount=0;icount<yde->nelem;icount++)
  { for(jcount=0;jcount<yde->nelno;jcount++)
    { INTw(fout,yde->i2elto[jcount][icount],6);
      CHRwsp(fout);
    } 
    CHRwcr(fout);
}

    CHRw(fout,"/YD/YDE/I2ELBNP");
    CHRwsp(fout);
    INTw(fout,21,5);
    CHRwsp(fout);
    INTw(fout,4,5);
    CHRwsp(fout);
    INTw(fout,yde->nelem,5);
    CHRwcr(fout);
    
    for(icount=0;icount<yde->nelem;icount++)
    { for(jcount=0;jcount<4;jcount++)
    { INTw(fout,i2elbnp[jcount][icount],6);
        CHRwsp(fout);
    }
        CHRwcr(fout);
    }

    
  CHRw(fout,"/YD/YDE/D3TCS");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/D1EMCT");
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);

 }
/****************** Interactions *************************************/
static void Ywdi(ydi,fout)
  YDI ydi; FILE *fout;
{ CHRw(fout,"     /*   Interactions     */");
  CHRwcr(fout);
  CHRw(fout,"/YD/YDI/MICOUP");
  CHRwsp(fout);
  INTw(fout,ydi->micoup,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/NICOUP");
  CHRwsp(fout);
  INTw(fout,ydi->nicoup,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/IIECFF");
  CHRwsp(fout);
  INTw(fout,ydi->iiecff,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/DIEDI");
  CHRwsp(fout);
  DBLw(fout,ydi->diedi,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/DIEZON");
  CHRwsp(fout);
  DBLw(fout,ydi->diezon,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/D1IESL");
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/I1IECN");
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/I1IECT");
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/D1DELTAT1");
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);
/*
  for(icount=0;icount<ydi->micoup;icount++)
  {
    DBLw(fout,ydi->d1deltat1[icount],19);
	CHRwsp(fout);
  }
  CHRwcr(fout);
*/
  CHRw(fout,"/YD/YDI/D1DELTAT2");
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);
/*
  for(icount=0;icount<ydi->micoup;icount++)
  {
    DBLw(fout,ydi->d1deltat2[icount],19);
	CHRwsp(fout);
  }
  CHRwcr(fout);
*/
  CHRw(fout,"/YD/YDI/D1DELTAN");
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);
/*
  for(icount=0;icount<ydi->micoup;icount++)
  {
    DBLw(fout,ydi->d1deltan[icount],19);
	CHRwsp(fout);
  }
  CHRwcr(fout);
*/
  CHRw(fout,"/YD/YDI/D2NV");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,3,5);
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwsp(fout);
  CHRwcr(fout);
/*
  for(icount=0;icount<ydi->micoup;icount++)
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwsp(fout);
  CHRwcr(fout);
/*
  for(icount=0;icount<ydi->micoup;icount++)
  {for(jcount=0;jcount<3;jcount++)
  {
    DBLw(fout,ydi->d2t2v[jcount][icount],19);
	CHRwsp(fout);
  }
  CHRwcr(fout);
  }
  CHRwcr(fout);
*/

  CHRw(fout,"/YD/YDI/D2T1V");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,3,5);
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwsp(fout);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/D2T2V");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,3,5);
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwsp(fout);
  CHRwcr(fout);
 }
 


/****************** Nodes ******************************************/
static void Ywdn(ydn,fout)
  YDN ydn; FILE *fout;
{ INT icount,jcount;
 
  CHRw(fout,"     /*   Nodes     */");
  CHRwcr(fout);
  CHRw(fout,"/YD/YDN/MNODIM");
  CHRwsp(fout);
  INTw(fout,ydn->mnodim,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/NNODIM");
  CHRwsp(fout);
  INTw(fout,ydn->nnodim,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/MNOPO");
  CHRwsp(fout);
  INTw(fout,ydn->mnopo,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/NNOPO");
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/D2NCC");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnodim,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { for(jcount=0;jcount<ydn->nnodim;jcount++)
    { DBLw(fout,ydn->d2ncc[jcount][icount],19);
      CHRwsp(fout);
    }
    CHRwcr(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/D2NCI");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnodim,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { for(jcount=0;jcount<ydn->nnodim;jcount++)
    { DBLw(fout,ydn->d2nci[jcount][icount],19);
      CHRwsp(fout);    
    }
    CHRwcr(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/D2NFC");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnodim,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { for(jcount=0;jcount<ydn->nnodim;jcount++)
    { DBLw(fout,ydn->d2nfc[jcount][icount],19);
      CHRwsp(fout);
    }
    CHRwcr(fout);
  }
  CHRwcr(fout); 

  CHRw(fout,"/YD/YDN/D2NFT");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnodim,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { for(jcount=0;jcount<ydn->nnodim;jcount++)
    { DBLw(fout,ydn->d2nfc[jcount][icount],19);
      CHRwsp(fout);
    }
    CHRwcr(fout);
  }
  CHRwcr(fout); 


  CHRw(fout,"/YD/YDN/D1NMCT");
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);
  for(icount=0;icount<0;icount++)
  { DBLw(fout,ydn->d1nmct[icount],19);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/D2NVC");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnodim,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { for(jcount=0;jcount<ydn->nnodim;jcount++)
    { DBLw(fout,ydn->d2nvc[jcount][icount],19);
      CHRwsp(fout);
    }
    CHRwcr(fout);
  }
  CHRwcr(fout); 

  CHRw(fout,"/YD/YDN/I1NOBF");
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { INTw(fout,ydn->i1nobf[icount],4);   
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/I1NOPR");
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { INTw(fout,ydn->i1nopr[icount],4);   
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/D1NTI");
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { DBLw(fout,ydn->d1nti[icount],10);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/D1NTC");
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { DBLw(fout,ydn->d1nti[icount],10);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/D1NTP");
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { DBLw(fout,ydn->d1nti[icount],10);
    CHRwsp(fout);
  }
  CHRwcr(fout);
}
/****************** Output ******************************************/
static void Ywdo(ydo,fout)
  YDO ydo; FILE *fout;
{ INT icount;

  CHRw(fout,"     /*   Output     */");
  CHRwcr(fout);
  CHRw(fout,"/YD/YDO/MOHYS");
  CHRwsp(fout);
  INTw(fout,ydo->mohys,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/NOHYS");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/DOHYP");
  CHRwsp(fout);
  DBLw(fout,ydo->dohyp,19);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/D1OHYC");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);
  for(icount=0;icount<ydo->nohys;icount++)
  { DBLw(fout,ydo->d1ohyc[icount],19);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/D1OHYF");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);
  for(icount=0;icount<ydo->nohys;icount++)
  { DBLw(fout,ydo->d1ohyf[icount],19);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/D1OHYS");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);
  for(icount=0;icount<ydo->nohys;icount++)
  { DBLw(fout,ydo->d1ohys[icount],19);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/D1OHYT");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);
  for(icount=0;icount<ydo->nohys;icount++)
  { DBLw(fout,ydo->d1ohyt[icount],19);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/D1OHYX");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);
  for(icount=0;icount<ydo->nohys;icount++)
  { DBLw(fout,ydo->d1ohyx[icount],19);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/D1OHYY");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);
  for(icount=0;icount<ydo->nohys;icount++)
  { DBLw(fout,ydo->d1ohyy[icount],19);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/I1OHYT");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);
  for(icount=0;icount<ydo->nohys;icount++)
  { INTw(fout,ydo->i1ohyt[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);
}  
/****************** Properties **************************************/
static void Ywdp(ydp,fout)
  YDP ydp; FILE *fout;
{ INT icount;

  CHRw(fout,"     /*   Properties     */");
  CHRwcr(fout);
  CHRw(fout,"/YD/YDP/MPROP");
  CHRwsp(fout);
  INTw(fout,ydp->mprop,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/NPROP");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEKS");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1peks[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEFR");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pefr[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  /*
  CHRw(fout,"/YD/YDP/D1PEFS");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pefs[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);
  */

  CHRw(fout,"/YD/YDP/D1PCOH");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pcoh[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PICF");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1picf[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEFT");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1peft[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEGFN");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pegfn[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEGFS");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pegfs[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PELA");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pela[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEMU");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pemu[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEPE");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pepe[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout); 

  CHRw(fout,"/YD/YDP/D1PERO");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pero[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PESF");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pesf[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

    CHRw(fout,"/YD/YDP/D1PEPSF");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pepsf[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEVF");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pevf[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEPF");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pepf[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);


  CHRw(fout,"/YD/YDP/I1PTYP");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { INTw(fout,ydp->i1ptyp[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/I1PEJP");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { INTw(fout,ydp->i1pejp[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/I1PEMN");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { INTw(fout,ydp->i1pemn[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/I1PSDE");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { INTw(fout,ydp->i1psde[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1TCON");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1tcon[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1CAPA");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1capa[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1CTEX");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1ctex[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);
}

/****************** Boundary Conditions **************************************/
static void Ywdb(ydb,fout)
YDB ydb; FILE *fout;
{ INT icount, jcount;

  CHRw(fout,"     /*   Boundary Conditions     */");
  CHRwcr(fout);
  CHRw(fout,"/YD/YDB/MBCON");
  CHRwsp(fout);
  INTw(fout,ydb->mbcon,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/NBCON");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/D1BNVX");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { DBLw(fout,ydb->d1bnvx[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/D1BNVY");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { DBLw(fout,ydb->d1bnvy[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/D1BNVZ");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { DBLw(fout,ydb->d1bnvz[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/I1BNVX");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { INTw(fout,ydb->i1bnvx[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/I1BNVY");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { INTw(fout,ydb->i1bnvy[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/I1BNVZ");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { INTw(fout,ydb->i1bnvz[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);


  CHRw(fout,"/YD/YDB/D1BNFX");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { DBLw(fout,ydb->d1bnfx[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/D1BNFY");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { DBLw(fout,ydb->d1bnfy[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/D1BNFZ");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { DBLw(fout,ydb->d1bnfz[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);
CHRw(fout,"/YD/YDB/D1BNAX");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { DBLw(fout,ydb->d1bnax[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/D1BNAY");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { DBLw(fout,ydb->d1bnay[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/D1BNAZ");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { DBLw(fout,ydb->d1bnaz[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/I1BNTP");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { INTw(fout,ydb->i1bntp[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/D1BNTP");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { DBLw(fout,ydb->d1bntp[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/I1BNHF");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { INTw(fout,ydb->i1bnhf[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/D1BNHF");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { DBLw(fout,ydb->d1bnhf[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/I1BCVT");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { INTw(fout,ydb->i1bcvt[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/D1BCVT");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { DBLw(fout,ydb->d1bcvt[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDB/D1BCVC");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { DBLw(fout,ydb->d1bcvc[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);
}
/****************** Boundary mesh **************************************/
static void Ywdx(ydx,fout, d1nap_sur)
  YDX ydx; FILE *fout; DBL *d1nap_sur;
{ INT icount, jcount;

  CHRw(fout,"     /*   Boundary Mesh     */");
  CHRwcr(fout);
  CHRw(fout,"/YD/YDX/NELNO");
  CHRwsp(fout);
  INTw(fout,ydx->nelno,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDX/NELEM");
  CHRwsp(fout);
  INTw(fout,ydx->nelem,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDX/DTHICK");
  CHRwsp(fout);
  DBLw(fout,ydx->dthick,10);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDX/I2ELTO");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,ydx->nelno,5);
  CHRwsp(fout);
  INTw(fout,ydx->nelem,5);
  CHRwcr(fout);

  for(icount=0;icount<ydx->nelem;icount++)
  { for(jcount=0;jcount<ydx->nelno;jcount++)
    { INTw(fout,ydx->i2elto[jcount][icount],6);
      CHRwsp(fout);
    } 
    CHRwcr(fout);
  }
    CHRwcr(fout);
    CHRw(fout,"/YD/YDX/D1NAP");
    CHRwsp(fout);
    INTw(fout,ydx->nelem,5);
    CHRwcr(fout);
    for(icount=0;icount<ydx->nelem;icount++)
    { DBLw(fout,d1nap_sur[icount],18);
        CHRwsp(fout);
    }
    CHRwcr(fout);
}

static void Yboundary(i2elto,nlist,i1nlist,i2elto_list,d1nap,i2elbnp,d1nap_sur,i2elto_t,nelem_t)
INT **i2elto; INT nlist; INT *i1nlist; INT **i2elto_list; DBL *d1nap;
INT **i2elbnp; DBL *d1nap_sur; INT **i2elto_t; INT nelem_t;
{
    INT i, j,k,m,n,ii,jj,kk,innopo,jnnopo,knnopo;
    for(n=0;n<nlist;n++)
    {
        innopo=i2elto_list[0][n];
        jnnopo=i2elto_list[1][n];
        knnopo=i2elto_list[2][n];

        for(m=0;m<nelem_t;m++)
        {
            
            ii=i2elto_t[0][m];
            jj=i2elto_t[1][m];
            kk=i2elto_t[2][m];
            if(((innopo==ii)&&(jnnopo==jj)&&(knnopo==kk))||
               ((innopo==ii)&&(jnnopo==kk)&&(knnopo==jj))||
               ((innopo==jj)&&(jnnopo==ii)&&(knnopo==kk))||
               ((innopo==jj)&&(jnnopo==kk)&&(knnopo==ii))||
               ((innopo==kk)&&(jnnopo==ii)&&(knnopo==jj))||
               ((innopo==kk)&&(jnnopo==jj)&&(knnopo==ii)))
            {
				d1nap_sur[m]=d1nap[n];
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

            
        ii=i2elto[i][i1nlist[n]-1];
        jj=i2elto[j][i1nlist[n]-1];
        kk=i2elto[k][i1nlist[n]-1];
      if(((innopo==ii)&&(jnnopo==jj)&&(knnopo==kk))||
         ((innopo==ii)&&(jnnopo==kk)&&(knnopo==jj))||
         ((innopo==jj)&&(jnnopo==ii)&&(knnopo==kk))||
         ((innopo==jj)&&(jnnopo==kk)&&(knnopo==ii))||
         ((innopo==kk)&&(jnnopo==ii)&&(knnopo==jj))||
         ((innopo==kk)&&(jnnopo==jj)&&(knnopo==ii)))i2elbnp[i+j+k-3][i1nlist[n]-1]=m;
    }
            continue;
        }
		}
	}
}

main(argc, argv)
INT argc; char **argv;
{ 
     INT icount,i,j,itmp,k,index,index1,inode,index2;
     CHR name[300];
     CHR namef[300];
     CHR c1name[100];
     CHR c1name1[100];
     CHR constraint[100];
//DBL *d1vix,*d1viy,*d1viz;
     INT *i1nobf_t,**i2elbnp,**i2elto_list;
	 INT nlist,*i1nlist;
	 INT nDBLcons = 13, nINTcons = 9;
     DBL d2cond[nDBLcons][N],d1cond_tmp[nDBLcons];
     INT Ncond,i2cond[nINTcons][N],i1cond_tmp[nINTcons];
     INT istart,ivel,idis,pindex,vindex,irota,inumber_t;
     DBL tmp1,tmp2;
     DBL *d1nap,*d1nap_sur;
     INT Nmat,m,ielem,imat_tmp,Nmat_tmp,iprop,imesh,imestyp,ptyp_tmp;
     INT i1nelem[N],i1pemnn[N],i1pemn_tmp[N],i1cond_index[N],i2mat[4][N];
     DBL **d2mat;
     INT NmatProp = 18;
     CHR strthree[100],model[100],mat_tmp[100];
   
     struct YD_struct yd;
     YDC ydc=&(yd.ydc);
     YDE yde=&(yd.yde);
     YDI ydi=&(yd.ydi);
     YDN ydn=&(yd.ydn);
     YDO ydo=&(yd.ydo);
     YDP ydp=&(yd.ydp);
     YDB ydb=&(yd.ydb);
     YDX ydx=&(yd.ydx);

     FILE *fgid=FILENULL;
     FILE *fp=FILENULL;
     FILE *ftmp=FILENULL;


     d2mat=DBL2NULL;

     CHRw(stdout,"******************************************************");
     CHRw(stdout,"****************\n");
     CHRw(stdout,"** Copyright (C) 2008,");
     CHRw(stdout,"                                              **\n");
     CHRw(stdout,"** Queen Mary University of London (QMUL) & Imperial "); 
     CHRw(stdout,"College        **\n");
     CHRw(stdout,"** of Science, Technology and Medicine (ICSTM). "); 
     CHRw(stdout,"All rights reserved.**\n");
     CHRw(stdout,"** Implemented for you by Prof Antonio Munjiza & ");
     CHRw(stdout,"Dr Jiansheng Xiang **\n");
     CHRwcr(stdout);
     CHRwcr(stdout);
     CHRwcr(stdout);
     CHRw(stdout,"* This code is part of the Virtual Geoscience ");
     CHRw(stdout,"Workbench (VGW) developed\n");
     CHRw(stdout,"* jointly by ICSTM and QMUL through two related ");
     CHRw(stdout,"parallel projects at \n");
     CHRw(stdout,"* ICSTM and QMUL respectively funded by EPSRC. \n");
     CHRw(stdout,"*  \n");
     CHRw(stdout,"* This code is provided by copyright holders under");
     CHRw(stdout," the GNU lesser \n");
     CHRw(stdout,"* General Public License (LGPL)."); 
     CHRw(stdout," It is open source code; you can \n");
     CHRw(stdout,"* redistribute it and/or modify it under the terms ");
     CHRw(stdout,"of the GNU Lesser\n");
     CHRw(stdout,"* General Public License version 3.  \n");
     CHRw(stdout,"*  \n");
     CHRw(stdout,"* This code is distributed in the hope that it will ");
     CHRw(stdout,"be useful, \n");
     CHRw(stdout,"* but WITHOUT ANY WARRANTY; without even the implied ");
     CHRw(stdout,"warranty \n");
     CHRw(stdout,"* of MERCHANTABILITY or FITNESS FOR A PARTICULAR ");
     CHRw(stdout,"PURPOSE. See \n");
     CHRw(stdout,"* the GNU Lesser General Public License");
     CHRw(stdout," for more details,\n");
     CHRw(stdout,"* http://www.gnu.org/licenses/lgpl.txt. \n");
     CHRw(stdout,"*  \n");
     CHRw(stdout,"* You should have received a copy of the GNU Lesser ");
     CHRw(stdout,"General Public\n");
     CHRw(stdout,"* License along with this code; if not, write to \n");
     CHRw(stdout,"* Dr Jiansheng Xiang Prof Antonio Munjiza ");
     CHRw(stdout,"or Dr John-Paul Latham \n");
     CHRw(stdout,"* j.xiang@imperial.ac.uk a.munjiza@qmul.ac.uk \n");
     CHRw(stdout,"* or j.p.latham@imperial.ac.uk  \n");
     CHRw(stdout,"* *************************************************");
     CHRw(stdout,"****************** *\n");

     Yrdd(&yd);
     CHRcpy(c1name,argv[1]);
     CHRcpy(c1name1,argv[2]);
     ydc->finp=fopen(c1name1,"r");
     fgid=fopen(c1name,"r");
     CHRcpynoext(namef,c1name);
     CHRcat(namef,".Y3D");
			
     fp=fopen(namef,"w");

     ydc->fcheck=fopen("Ytmp","w");

     SETLINEBUF(ydc->fcheck);  
     
     if(((ydc->finp)==(FILENULL))||((fgid)==(FILENULL)))
     { CHRw(stderr,"Could not open input file - usage -i inputfile");
	  CHRwcr(stderr);  
	  return 0;
     }
     else
     {
	  CHRr(fgid,name);

	  while(FILEND(fgid)==0) 
	  { if(CHRcmp(name, "/YD/",4)==0)
	       { CHRw(ydc->fcheck,name); CHRwcr(ydc->fcheck);
		    if(CHRcmp(name,"/YD/YDC/",8)==0)        /* read control data       */
		    { Yrdc(ydc,fgid,name);
		    }
		    else if(CHRcmp(name,"/YD/YDE/",8)==0)   /* read data elements      */
		    { Yrde(yde,fgid,name);
		    }
		    else if(CHRcmp(name,"/YD/YDI/",8)==0)   /* read data interaction   */
		    { Yrdi(ydi,fgid,name);
		    }
		    else if(CHRcmp(name,"/YD/YDN/",8)==0)  /* read data nodes          */
		    { Yrdn(ydn,fgid,name);
		    }
		    else if(CHRcmp(name,"/YD/YDO/",8)==0)  /* read data nodes          */
		    { Yrdo(ydo,fgid,name);
		    }
		    else if(CHRcmp(name,"/YD/YDP/",8)==0)  /* read data properties     */
		    { Yrdp(ydp,fgid,name);
		    }
		    else if(CHRcmp(name,"/YD/YDB/",8)==0)  /* read data properties     */
		    { Yrdb(ydb,fgid,name);
		    }
		    else if(CHRcmp(name,"/YD/YDX/",8)==0)  /* surface mesh             */
		    { 
//		    Yrdx(ydx,ydp,yde,fgid,name);
			    if(CHRcmp(name,"/YD/YDX/NELEM",13)==0)
			    {  
			    	INTr(fgid,&(ydx->nelem));
			    	d1nap_sur=TalDBL1(ydx->nelem);
			    }
			    if(CHRcmp(name,"/YD/YDX/NELNO",13)==0)
			    {  
			    	INTr(fgid,&(ydx->nelno));
			    }
			    if(CHRcmp(name,"/YD/YDX/DTHICK",14)==0)
			    {  
			    	DBLr(fgid,&(ydx->dthick));
			    }
			    if(CHRcmp(name,"/YD/YDX/I2ELTO",15)==0)
			    {   
			    	TformINT2(fgid,-1,ydx->nelno,ydx->nelem+1,&(ydx->i2elto));
			    }
		    }
		    else
		    { CHRw(stderr,"Yrd: unknown name ");
			 CHRw(stderr,name); 
			 CHRwcr(stderr);   
			 return 0;
		    } }
	       else if (CHRcmp(name,"Gen_out_only",12)==0)
	       {
		    CHRr(fgid,name);
		    if(CHRcmp(name,"Yes",3)==0)
		    {ftmp=fopen("1.tmp","w");
			 fclose(ftmp);}
	       }
           else if (CHRcmp(name,"LOAD_LIST",9)==0)
           {
               INTr(fgid,&nlist);
               i1nlist=TalINT1(nlist+1);
               i2elto_list=TalINT2(ydx->nelno+1,nlist+1);
               d1nap=TalDBL1(nlist+1);
               
               for(i=0;i<nlist;i++)
               {
                  INTr(fgid,&i1nlist[i]);
                   for(j=0;j<ydx->nelno;j++)
                   {
                       INTr(fgid,&i2elto_list[j][i]);
                       i2elto_list[j][i]--;
                   }
                   DBLr(fgid,&d1nap[i]);
               }
               
           }
	       else if (CHRcmp(name,"MATERIAL_LIST",9)==0)
	       {
		    //INTr(fgid,&ydp->nprop);
		    //ydp->mprop=ydp->nprop+1;

		    INTr(fgid,&Nmat);

		    //ydp->nprop=Nmat;
		    //ydp->mprop=ydp->nprop+1;
		    
		    //if(ydp->nprop<1)
		    if(Nmat<1)
		    {
			 CHRw(stderr,"Material properties should be setup properly"); 
			 CHRwcr(stderr);  
			 return 0;
		    }
/*		    else
		    {
			 TformDBL2(fgid,R0,15,100,&d2mat);
 		    }
*/
		    else
		    {
			 d2mat=TalDBL2(NmatProp,100);

			 for(i=0;i<Nmat;i++)
			 {
			      DBLr(fgid,&d2mat[0][i]);
			      DBLr(fgid,&d2mat[1][i]);
			      DBLr(fgid,&d2mat[2][i]);
			      DBLr(fgid,&d2mat[3][i]);
			      DBLr(fgid,&d2mat[4][i]);
			      DBLr(fgid,&d2mat[5][i]);
			      DBLr(fgid,&d2mat[6][i]);
			      DBLr(fgid,&d2mat[7][i]);
			      DBLr(fgid,&d2mat[8][i]);
			      DBLr(fgid,&d2mat[9][i]);
			      DBLr(fgid,&d2mat[10][i]);
			      INTr(fgid,&vindex);
			      DBLr(fgid,&d2mat[11][i]);
			      DBLr(fgid,&d2mat[12][i]);
			      INTr(fgid,&pindex);
			      DBLr(fgid,&d2mat[13][i]);
			      DBLr(fgid,&d2mat[14][i]);
			      DBLr(fgid,&d2mat[15][i]);
				  DBLr(fgid,&d2mat[16][i]);
			      DBLr(fgid,&d2mat[17][i]);
			      INTr(fgid,&ptyp_tmp);

			      tmp1=d2mat[1][i];
			      tmp2=d2mat[2][i];
			      d2mat[2][i]=tmp1/(R2*(R1+tmp2));
			      d2mat[1][i]=tmp1*tmp2/((R1+tmp2)*(R1-R2*tmp2));
			      if(vindex==0)
			      {
				   d2mat[11][i]=R0;
			      }
			      if(pindex==0)
			      {
				   d2mat[13][i]=d2mat[10][i];
			      }
			 }
                    }

/*		    else
		    {
			 for(i=0;i<ydp->nprop;i++)
			 {
			      DBLr(fgid,&ydp->d1pero[i]);
			      DBLr(fgid,&ydp->d1pela[i]);
			      DBLr(fgid,&ydp->d1pemu[i]);
			      DBLr(fgid,&ydp->d1peks[i]);
			      DBLr(fgid,&ydp->d1pepe[i]);
			      DBLr(fgid,&ydp->d1pegfn[i]);
						DBLr(fgid,&ydp->d1pegfs[i]);
			      DBLr(fgid,&ydp->d1peft[i]);
			      DBLr(fgid,&ydp->d1pefs[i]);
			      DBLr(fgid,&ydp->d1pefr[i]);
			      INTr(fgid,&vindex);
			      DBLr(fgid,&ydp->d1pesf[i]);
			      DBLr(fgid,&ydp->d1pevf[i]);
			      INTr(fgid,&pindex);
			      DBLr(fgid,&ydp->d1pepf[i]);
			      DBLr(fgid,&ydp->d1pepsf[i]);
			      INTr(fgid,&ydp->i1ptyp[i]);
			      tmp1=ydp->d1pela[i];
			      tmp2=ydp->d1pemu[i];
			      ydp->d1pemu[i]=tmp1/(R2*(R1+tmp2));
			      ydp->d1pela[i]=tmp1*tmp2/((R1+tmp2)*(R1-R2*tmp2));
			      if(vindex==0)
			      {
				   ydp->d1pesf[i]=R0;
			      }
			      if(pindex==0)
			      {
				   ydp->d1pepf[i]=ydp->d1pefr[i];
			      }
			 }         
		    }
*/	       }
	       
	       else if (CHRcmp(name,"Initial_Data",12)==0)
	       {
	       INTr(fgid,&icount);
		    for(i=0;i<icount;i++)
		    {   INTr(fgid,&inode);
			 inode=inode-1;
			 INTr(fgid,&ivel);
			 if(ivel==1)
			 {
			      DBLr(fgid,&ydn->d2nvc[0][inode]);
			      DBLr(fgid,&ydn->d2nvc[1][inode]);
			      DBLr(fgid,&ydn->d2nvc[2][inode]);
			 }
			 else
			 {
			      CHRr(fgid,name);
			      CHRr(fgid,name);
			      CHRr(fgid,name);
			 }
			 INTr(fgid,&idis);
			 if(idis==1)
			 {
			      DBLr(fgid,&ydn->d2ncc[0][inode]);
			      DBLr(fgid,&ydn->d2ncc[1][inode]);
			      DBLr(fgid,&ydn->d2ncc[2][inode]);
			      ydn->d2ncc[0][inode]=ydn->d2ncc[0][inode]+ydn->d2nci[0][inode];
			      ydn->d2ncc[1][inode]=ydn->d2ncc[1][inode]+ydn->d2nci[1][inode];
			      ydn->d2ncc[2][inode]=ydn->d2ncc[2][inode]+ydn->d2nci[2][inode];
			 }
			 else
			 {
			      CHRr(fgid,name);
			      CHRr(fgid,name);
			      CHRr(fgid,name);
			 }
			 INTr(fgid,&idis);
			 if(idis==1)
			 {DBLr(fgid,&ydn->d1nti[inode]);
			 }else
			 {
			 	CHRr(fgid,name);
			 }
		    }
   
	       }

	       else if (CHRcmp(name,"Element_List",10)==0)
	       {
		    for(i=0;i<yde->nelem;i++)
		    {
			 yde->i1elpr[i]=yde->i1elpr[i]-1;
		    }

		    for(i=0;i<N;i++)
		    {
			 for(j=0;j<4;j++)
			 {
			      i2mat[j][i]=0;
			 }
		    }
		    INTr(fgid,&icount);
		    for(m=0;m<icount;m++)
		    {
			 //CHRcpy(strthree,"3");
                         CHRcpy(strthree,"0");
			 INTr(fgid,&ielem);
			 ielem=ielem-1;
			 CHRr(fgid,&model);
			 CHRr(fgid,&mat_tmp);
			 itmp=0;
			 if(CHRcmp(model,"Fracture",8)==0)
			 {
			      //CHRcat(strthree,mat_tmp);
                              CHRcpy(strthree,"3");
			      itmp=3;
			 }
			 imat_tmp=atoi(strthree);
			 k=yde->i1elpr[ielem];
			 index1=0;
			 for(i=0;i<Nmat;i++)
			 {
			      index=0;
			      for(j=0;j<NmatProp;j++)
			      {
				   if(ABS(d2mat[j][i]-d2mat[j][k])>EPSILON)
				   {
					index=1;
					break;
				   }

			      }

			      if((index==0)&&(i2mat[1][i]==imat_tmp))
			      {index1=1;
				   yde->i1elpr[ielem]=i;
				   break;
			      }
			 }
			 if(index1==0)
			 {
			      Nmat=Nmat+1;
			      for(j=0;j<NmatProp;j++)
			      {
				   d2mat[j][Nmat-1]=d2mat[j][k];
			      }
			      i2mat[1][Nmat-1]=imat_tmp;
			      i2mat[3][Nmat-1]=itmp;
			      yde->i1elpr[ielem]=Nmat-1;
			 }
		    }
		    for(i=0;i<N;i++)
		    {
			 i1cond_index[i]=i;
		    }

		    i=0;
		    while(i<Nmat)
		    {
			 index=0;
			 for(j=0;j<yde->nelem;j++)
			 {
			      if(yde->i1elpr[j]==i1cond_index[i])
			      {   index=1;
				   break;
			      }
			 }
			 if(index==0)
			 {
			      for(j=0;j<NmatProp;j++)
			      {
				   d2mat[j][i]=d2mat[j][Nmat-1];

			      }
			      i2mat[1][i]=i2mat[1][Nmat-1];
			      i2mat[3][i]=i2mat[3][Nmat-1];
			      i1cond_index[Nmat-1]=i;
			      i1cond_index[i]=Nmat-1;
			      Nmat--;
			 }
			 else
			      i++;

		    }

		    for(j=0;j<yde->nelem;j++)
		    {
			 yde->i1elpr[j]=i1cond_index[yde->i1elpr[j]];

		    }
		    Nmat_tmp=Nmat;
		    for(i=0;i<Nmat;i++)
		    {i2mat[2][i]=13;

			 if(i2mat[3][i]==3)
			 {
			      Nmat_tmp++;
			      i2mat[0][i]=Nmat_tmp-1;
			      i2mat[2][Nmat_tmp-1]=15;
			      i2mat[1][Nmat_tmp-1]=0;
			      i2mat[0][Nmat_tmp-1]=0;

			      for(j=0;j<NmatProp;j++)
			      {
				   d2mat[j][Nmat_tmp-1]=d2mat[j][i];
			      }
			 }
		    }
		    Nmat=Nmat_tmp;
                    ydp->nprop=Nmat;
                    ydp->mprop=ydp->nprop+1;
              
		    for(i=0;i<Nmat;i++)
		    {
			 ydp->d1pero[i]=d2mat[0][i];
			 ydp->d1pela[i]=d2mat[1][i];
			 ydp->d1pemu[i]=d2mat[2][i];
			 ydp->d1peks[i]=d2mat[3][i];
			 ydp->d1pepe[i]=d2mat[4][i];
			 ydp->d1pegfn[i]=d2mat[5][i];
			 ydp->d1pegfs[i]=d2mat[6][i];
			 ydp->d1peft[i]=d2mat[7][i];
			 ydp->d1pcoh[i]=d2mat[8][i];
			 ydp->d1picf[i]=d2mat[9][i];
			 ydp->d1pefr[i]=d2mat[10][i];
			 ydp->d1pesf[i]=d2mat[11][i];
			 ydp->d1pevf[i]=d2mat[12][i];
			 ydp->d1pepf[i]=d2mat[13][i];
			 ydp->d1pepsf[i]=d2mat[14][i];		 
			 ydp->d1tcon[i]=d2mat[15][i];
			 ydp->d1capa[i]=d2mat[16][i];
			 ydp->d1ctex[i]=d2mat[17][i];

			 ydp->i1pemn[i]=i2mat[1][i];
			 ydp->i1ptyp[i]=i2mat[2][i];
			 ydp->i1pejp[i]=i2mat[0][i];
             ydp->i1psde[i]=0;
		    }
	       }

			else if (CHRcmp(name,"Constraints_List",16)==0)
			{
				for(i=0;i<nDBLcons;i++)
				{
					d2cond[i][0]=R0;
				}
				i2cond[0][0]=0;
				i2cond[1][0]=0;
				i2cond[2][0]=0;
				Ncond=1;
				INTr(fgid,&icount);
				for(i=0;i<icount;i++)
				{   
				     INTr(fgid,&inode);
					 inode=inode-1;
					 for(j=0;j<3;j++)
					 {
					      CHRr(fgid,&constraint);
					      if(CHRcmp(constraint,"Force",5)==0)i1cond_tmp[j]=0;
					      if(CHRcmp(constraint,"Acceleration",12)==0)i1cond_tmp[j]=1;
					      if(CHRcmp(constraint,"Velocity",8)==0)i1cond_tmp[j]=2;
					      DBLr(fgid,&d1cond_tmp[j]);
					 }

					 /* reading rotation paramenters */
					 for(j=0;j<3;j++)
					 {
					      INTr(fgid,&irota);
					      if(irota==1)
					      {
					      	i1cond_tmp[j]=3;
							DBLr(fgid,&d1cond_tmp[j]);
							DBLr(fgid,&d1cond_tmp[j+3]);
						}
					      else
					      {
					      	CHRr(fgid,&constraint);
						 	CHRr(fgid,&constraint);
						 	d1cond_tmp[j+3]=R0;
					      }
					 }

					 /* Reading temperature parameters */
					 /*Thermal constrains */
					  INTr(fgid,&irota);
					  if(irota==1)
					  {
					  	i1cond_tmp[6]=1;
					  	DBLr(fgid,&d1cond_tmp[9]);
					  }
					  else
					  { 
					  	CHRr(fgid,&constraint);
					  	i1cond_tmp[6]=0;
						d1cond_tmp[9]=R0;
					  }
					/* Flux */ 
					  INTr(fgid,&irota);
					  if(irota==1)
					  {
					  	i1cond_tmp[7]=1;
					  	DBLr(fgid,&d1cond_tmp[10]);
					  }
					  else
					  {	  
					  	CHRr(fgid,&constraint);
					  	i1cond_tmp[7]=0;
						d1cond_tmp[10]=R0;

					  }
					 /* Convection */
					  INTr(fgid,&irota);
					  if(irota==1){
						  i1cond_tmp[8]=1;
						  DBLr(fgid,&d1cond_tmp[11]);
						  DBLr(fgid,&d1cond_tmp[12]);
					  }
					  else
					  {	  CHRr(fgid,&constraint);CHRr(fgid,&constraint);
					  	  i1cond_tmp[8]=0;
						  d1cond_tmp[11]=R0;
						  d1cond_tmp[12]=R0;
					  }
					 
					 index1=0;

					/* Checking if condition already exists */ 
					for(k=0;k<Ncond;k++)
					{   
						index=0;
						for(j=0;j<nDBLcons;j++)
						{
							if(ABS(d2cond[j][k]-d1cond_tmp[j])>EPSILON)
							{
								index=1;
								break;
							}
						}

						if((index==0)
						&&(i2cond[0][k]==i1cond_tmp[0])
						&&(i2cond[1][k]==i1cond_tmp[1])
						&&(i2cond[2][k]==i1cond_tmp[2])
						&&(i2cond[3][k]==i1cond_tmp[3])
						&&(i2cond[4][k]==i1cond_tmp[4])
						&&(i2cond[5][k]==i1cond_tmp[5])
						&&(i2cond[6][k]==i1cond_tmp[6])
						&&(i2cond[7][k]==i1cond_tmp[7])
						&&(i2cond[8][k]==i1cond_tmp[8])
						) //Add new nINTcond here
						{
							ydn->i1nopr[inode]=k;
			   				index1=1;
			   				break;
						}
					}
				 /* If not existing adding new condition*/ 	 
				 if(index1==0)
				 {Ncond=Ncond+1;
				      ydn->i1nopr[inode]=Ncond-1;
				      for(j=0;j<nDBLcons;j++)
				      {
					   d2cond[j][Ncond-1]=d1cond_tmp[j];
				      }	  
				      for(j=0;j<nINTcons;j++)
				      {
					   i2cond[j][Ncond-1]=i1cond_tmp[j];
				      }
				 }
				}
				index=0;
				for(i=0;i<ydn->nnopo;i++)
				{
				 if(ydn->i1nopr[i]==0)
				 {index=1;
				      break;
				 }
				}
				if(index==1)istart=0;
				else 
				{istart=1;
				 Ncond=Ncond-1;
				}
				for(i=0;i<ydn->nnopo;i++)
				{
				 ydn->i1nopr[i]=ydn->i1nopr[i]-istart;
				}
				for(i=0;i<Ncond;i++)
				{j=i+istart;
				 if(i2cond[0][j]==0)
				 {
				      ydb->i1bnvx[i]=0;
				      ydb->d1bnvx[i]=R0;
				      ydb->d1bnax[i]=R0;
				      ydb->d1bnfx[i]=d2cond[0][j];
				 }
				 else if(i2cond[0][j]==1)
				 {
				      ydb->i1bnvx[i]=0;
				      ydb->d1bnvx[i]=R0;
				      ydb->d1bnax[i]=d2cond[0][j];
				      ydb->d1bnfx[i]=R0;
				 }
				 else if(i2cond[0][j]==2)
				 {
				      ydb->i1bnvx[i]=1;
				      ydb->d1bnvx[i]=d2cond[0][j];
				      ydb->d1bnax[i]=R0;
				      ydb->d1bnfx[i]=R0;
				 }
				 else if(i2cond[0][j]==3)
				 {
				      ydb->i1bnvx[i]=2;
				      ydb->d1bnvx[i]=d2cond[0][j];
				      ydb->d1bnfx[i]=d2cond[3][j];;
				      ydb->d1bnax[i]=R0;
				 }
				 if(i2cond[1][j]==0)
				 {
				      ydb->i1bnvy[i]=0;
				      ydb->d1bnvy[i]=R0;
				      ydb->d1bnay[i]=R0;
				      ydb->d1bnfy[i]=d2cond[1][j];
				 }
				 else if(i2cond[1][j]==1)
				 {
				      ydb->i1bnvy[i]=0;
				      ydb->d1bnvy[i]=R0;
				      ydb->d1bnay[i]=d2cond[1][j];
				      ydb->d1bnfy[i]=R0;
				 }
				 else if(i2cond[1][j]==2)
				 {
				      ydb->i1bnvy[i]=1;
				      ydb->d1bnvy[i]=d2cond[1][j];
				      ydb->d1bnay[i]=R0;
				      ydb->d1bnfy[i]=R0;
				 }
				 else if(i2cond[1][j]==3)
				 {
				      ydb->i1bnvy[i]=2;
				      ydb->d1bnvy[i]=d2cond[1][j];
				      ydb->d1bnfy[i]=d2cond[4][j];;
				      ydb->d1bnay[i]=R0;
				 }
				 if(i2cond[2][j]==0)
				 {
				      ydb->i1bnvz[i]=0;
				      ydb->d1bnvz[i]=R0;
				      ydb->d1bnaz[i]=R0;
				      ydb->d1bnfz[i]=d2cond[2][j];
				 }
				 else if(i2cond[2][j]==1)
				 {
				      ydb->i1bnvz[i]=0;
				      ydb->d1bnvz[i]=R0;
				      ydb->d1bnaz[i]=d2cond[2][j];
				      ydb->d1bnfz[i]=R0;
				 }
				 else if(i2cond[2][j]==2)
				 {
				      ydb->i1bnvz[i]=1;
				      ydb->d1bnvz[i]=d2cond[2][j];
				      ydb->d1bnaz[i]=R0;
				      ydb->d1bnfz[i]=R0;
				 }
				 else if(i2cond[2][j]==3)
				 {
				      ydb->i1bnvz[i]=2;
				      ydb->d1bnvz[i]=d2cond[2][j];
				      ydb->d1bnfz[i]=d2cond[5][j];;
				      ydb->d1bnaz[i]=R0;
				 }
				 if(i2cond[6][j]==1){
					 ydb->d1bntp[i]=d2cond[9][j];
					 ydb->i1bntp[i]=1;
				 }else{
					 ydb->d1bntp[i]=R0;
					 ydb->i1bntp[i]=0;
				 }
				 if(i2cond[7][j]==1){
				 	ydb->d1bnhf[i]=d2cond[10][j];
				 	ydb->i1bnhf[i]=1;
				 }else{
					 ydb->d1bnhf[i]=R0;
					 ydb->i1bnhf[i]=0;
				 }
				 if(i2cond[8][j]==1){
					 ydb->d1bcvt[i]=d2cond[11][j];
					 ydb->d1bcvc[i]=d2cond[12][j];
					 ydb->i1bcvt[i]=1;
				 }else{
				 	ydb->d1bcvt[i]=R0;
				 	ydb->d1bcvc[i]=R0;
				 	ydb->i1bcvt[i]=0;
				 }
				}
				ydb->nbcon=Ncond;
				ydb->mbcon=ydb->nbcon+1;
	}
	  CHRr(fgid,name);
	  }
	  
	  inumber_t=yde->nelno;
	  if((yde->nelno==11)||(yde->nelno==5))
	  {yde->nelno=yde->nelno-1;
	  }

	  for(i=0;i<N;i++)
	  {
	       i1nelem[i]=0;
	       i1pemnn[i]=0;
	  }

	  for(iprop=0;iprop<ydp->nprop;iprop++)
	  {i1pemn_tmp[iprop]=ydp->i1pemn[iprop];
	  }

	  for(imesh=0;imesh<10;imesh++)
	  {
	       for(iprop=0;iprop<ydp->nprop;iprop++)
	       {
		    imestyp=i1pemn_tmp[iprop]%10;
		    i1pemn_tmp[iprop]=i1pemn_tmp[iprop]/10;
		    if((imestyp==1)||(imestyp==2))
		    {
			 i1pemnn[iprop]++;
		    }
	       }
	  }

	  for(i=0;i<yde->nelem;i++)
	  {
	       i1nelem[i1pemnn[yde->i1elpr[i]]]++;
	  }

	  for(i=1;i<N;i++)
	  {
	       if(i1nelem[i]!=0)
	       {
		    yde->melem=yde->melem+i1nelem[i]*(pow(4,i)-1);
		    ydn->mnopo=ydn->mnopo+i1nelem[i]*pow(4,i-1);
	       }
	  }
	  yde->melem=yde->melem*10;
	  ydn->mnopo=ydn->mnopo*50;


//	  for(i=0;i<yde->nelem;i++)
//	  { yde->i1elpr[i]=yde->i1elpr[i]-1;
//	  }

	  for(i=0;i<yde->nelem;i++)
	  { for(j=0;j<yde->nelno;j++)
	       { yde->i2elto[j][i]=yde->i2elto[j][i]-1;
	       } 

	  }

	  for(i=0;i<ydx->nelem;i++)
	  { for(j=0;j<3;j++)
	       { ydx->i2elto[j][i]=ydx->i2elto[j][i]-1;
	       } 
	  }
	  i1nobf_t=TalINT1(ydn->nnopo);
	  for(i=0;i<ydn->nnopo;i++)
	  {
	       i1nobf_t[i]=0;
	  }
	  CHRr(ydc->finp,name);

	  while(FILEND(ydc->finp)==0) 
	  { if(CHRcmp(name, "$YSTOP",6)==0)
	       { CHRw(ydc->fcheck,name); CHRwcr(ydc->fcheck);
		    return 0;
	       }
	       else if(CHRcmp(name, "$YDOIT",6)==0)
	       { CHRw(ydc->fcheck,name); CHRwcr(ydc->fcheck);
		    return 1;
	       }
	       else if(CHRcmp(name,"/*",2)==0)    /* read and ignore comments */
	       { icount=0;
		    do
		    { CHRr(ydc->finp,name); icount++;
			 if(icount>100)
			 { CHRw(stderr,"Yrd: too long comment near - ");
			      CHRw(stderr,name);
			      CHRwcr(stderr);      
			      return 0;
			 } 
		    }while((FILEND(ydc->finp)==0)&&(CHRcmp(name,"*/",2)!=0));
	       }
	       else if(CHRcmp(name, "/YD/",4)==0)
	       { CHRw(ydc->fcheck,name); CHRwcr(ydc->fcheck);
		    if(CHRcmp(name,"/YD/YDC/",8)==0)        /* read control data       */
		    { Yrdc(ydc,ydc->finp,name);
		    }
		    else if(CHRcmp(name,"/YD/YDE/",8)==0)   /* read data elements      */
		    { Yrde(yde,ydc->finp,name);
		    }
		    else if(CHRcmp(name,"/YD/YDI/",8)==0)   /* read data interaction   */
		    { Yrdi(ydi,ydc->finp,name);
		    }
		    else if(CHRcmp(name,"/YD/YDN/",8)==0)  /* read data nodes          */
		    { Yrdn(ydn,ydc->finp,name);
		    }
		    else if(CHRcmp(name,"/YD/YDO/",8)==0)  /* read data nodes          */
		    { Yrdo(ydo,ydc->finp,name);
		    }
		    else if(CHRcmp(name,"/YD/YDP/",8)==0)  /* read data properties     */
		    { Yrdp(ydp,ydc->finp,name);
		    }
		    else if(CHRcmp(name,"/YD/YDB/",8)==0)  /* read data properties     */
		    { Yrdb(ydb,ydc->finp,name);
		    }
		    else if(CHRcmp(name,"/YD/YDX/",8)==0)  /* surface mesh             */
		    { 
//		    	Yrdx(ydx,ydp,yde,fgid,name);
			if(CHRcmp(name,"/YD/YDX/NELEM",13)==0)
			{  
				INTr(fgid,&(ydx->nelem));
				d1nap_sur=TalDBL1(ydx->nelem);
			}
			if(CHRcmp(name,"/YD/YDX/NELNO",13)==0)
			{
				INTr(fgid,&(ydx->nelno));
			}
			if(CHRcmp(name,"/YD/YDX/DTHICK",14)==0)
			{  
			    	DBLr(fgid,&(ydx->dthick));
			}
			if(CHRcmp(name,"/YD/YDX/I2ELTO",15)==0)
			{
				TformINT2(fgid,-1,ydx->nelno,ydx->nelem+1,&(ydx->i2elto));
				i1nobf_t=TalINT1(ydn->mnopo);
			}
		   }
		    else
		    { CHRw(stderr,"Yrd: unknown name: ");
			 CHRw(stderr,name); 
			 CHRwcr(stderr);   
			 return 0;
		    } }
	       CHRr(ydc->finp,name);
	  }

	  ydc->ncstep=0;

	  index2=0;
	  for(i=0;i<ydn->nnopo;i++)
	  {
	       if(ydn->i1nobf[i]>0)
	       { index2=1;
		    break;
	       }
	  }



/*assign layer number as facenode  */
	  if(index2==0)
	  {
	       for(i=0;i<yde->nelem;i++)
	       { for(j=0;j<yde->nelno;j++)
		    { 
			 ydn->i1nobf[yde->i2elto[j][i]]=yde->i2elto[yde->nelno][i];
		    } 
	       }
	  }
/*end */

	  for(i=0;i<ydx->nelem;i++)
	  { for(j=0;j<3;j++)
	       { if(ydn->i1nobf[ydx->i2elto[j][i]]>0)
			 i1nobf_t[ydx->i2elto[j][i]]=ydn->i1nobf[ydx->i2elto[j][i]];
	       } 
	  }
	  if(ydx->nelem!=0)
	  {
	       for(i=0;i<ydn->nnopo;i++)
	       {
		    ydn->i1nobf[i]=i1nobf_t[i];
	       }
	  }
	  for(i=0;i<yde->nelem;i++)
	  {
	       for(j=0;j<yde->nelno;j++)
	       {
		    if(ydn->i1nobf[yde->i2elto[j][i]]>0)
			 yde->i1elbe[i]=ydn->i1nobf[yde->i2elto[j][i]];
	       }

	  }

/*
  for(i=0;i<yde->nelem;i++)
  {
  for(j=0;j<yde->nelno;j++)
  {
  if(yde->i1elbe[i]>0)
  ydn->i1nobf[yde->i2elto[j][i]]=yde->i1elbe[i];
  }
  }
*/

	   yde->nelno=inumber_t;
         
         i2elbnp=TalINT2(4,yde->nelem+1);
         for(i=0;i<yde->nelem;i++)
         {
             for(j=0;j<4;j++)
             {
                 i2elbnp[j][i]=-1;
             }
         }
         for(i=0;i<ydx->nelem;i++)
         {
             d1nap_sur[i]=R0;
         }
        
        Yboundary(yde->i2elto,nlist,i1nlist,i2elto_list,d1nap,i2elbnp,d1nap_sur,ydx->i2elto,ydx->nelem);
	  Ywdc(ydc,fp);     /* Write Control variables into the input file     */
	  Ywde(yde,fp,i2elbnp);     /* Write Elements variables into the input file    */
	  Ywdi(ydi,fp);     /* Write Interaction variables into the input file */
	  Ywdn(ydn,fp);     /* Write Nodes variables into the input file       */
	  Ywdp(ydp,fp);     /* Write Properties variables into the input file  */
	  Ywdo(ydo,fp);     /* Write Output variables into the input file      */
	  Ywdb(ydb,fp);     /* Write Condition variables into the input file   */
	  Ywdx(ydx,fp,d1nap_sur);     /* Write boundary mesh data into the input file   */

	  CHRw(fp,"$YDOIT");
	  CHRwcr(fp);
	  CHRw(fp,"$YSTOP");
	  CHRwcr(fp);  

	  fclose(fp);
	  fclose(ydc->finp);
	  fclose(ydc->fcheck);

	  CHRw(stderr,"   ***** PROGRAM HAS ORDERLY FINISHED *****"); 
	  CHRwcr(stderr);
//CHRw(stderr,"Press a key to continue");
//CHRwcr(stderr);
//getchar();

//FREE(i1conindex);
	  FREE(i1nobf_t);

	  return 1;

     }

}	



