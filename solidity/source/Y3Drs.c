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

/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/



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

static void Ywde(yde,fout)
  YDE yde; FILE *fout;
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
  { INTw(fout,yde->i1elcf[icount],4);
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
    { INTw(fout,yde->i2elto[icount][jcount],6);
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
//  INTw(fout,0,5);
  INTw(fout,yde->nelem,5);
  CHRwcr(fout);
  for(icount=0;icount<yde->nelem;icount++)
  { DBLw(fout,yde->d1emct[icount],19);
    CHRwsp(fout);
  }
  CHRwcr(fout);
	


 }
/****************** Interactions *************************************/
static void Ywdi(ydi,fout)
  YDI ydi; FILE *fout;
{ INT icount,jcount;
  CHRw(fout,"     /*   Interactions     */");
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
  INTw(fout,ydi->micoup,5);
  CHRwcr(fout);

  for(icount=0;icount<ydi->micoup;icount++)
  {
    INTw(fout,ydi->i1iecn[icount],10);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/I1IECT");
  CHRwsp(fout);
  INTw(fout,ydi->nicoup,5);
  CHRwcr(fout);

  for(icount=0;icount<ydi->nicoup;icount++)
  {
    INTw(fout,ydi->i1iect[icount],10);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  
  CHRw(fout,"/YD/YDI/D1DELTAT1");
  CHRwsp(fout);
  INTw(fout,ydi->nicoup,5);
  CHRwcr(fout);

  for(icount=0;icount<ydi->nicoup;icount++)
  {
    DBLw(fout,ydi->d1deltat1[icount],19);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/D1DELTAT2");
  CHRwsp(fout);
  INTw(fout,ydi->nicoup,5);
  CHRwcr(fout);

  for(icount=0;icount<ydi->nicoup;icount++)
  {
    DBLw(fout,ydi->d1deltat2[icount],19);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/D1DELTAN");
  CHRwsp(fout);
  INTw(fout,ydi->nicoup,5);
  CHRwcr(fout);

  for(icount=0;icount<ydi->nicoup;icount++)
  {
    DBLw(fout,ydi->d1deltan[icount],19);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/D2NV");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,3,5);
  CHRwsp(fout);
  INTw(fout,ydi->nicoup,5);
  CHRwsp(fout);
  CHRwcr(fout);

  for(icount=0;icount<ydi->nicoup;icount++)
  {for(jcount=0;jcount<3;jcount++)
  {
    DBLw(fout,ydi->d2nv[jcount][icount],19);
	CHRwsp(fout);
  }
  CHRwcr(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/D2T1V");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,3,5);
  CHRwsp(fout);
  INTw(fout,ydi->nicoup,5);
  CHRwsp(fout);
  CHRwcr(fout);

  for(icount=0;icount<ydi->nicoup;icount++)
  {for(jcount=0;jcount<3;jcount++)
  {
    DBLw(fout,ydi->d2t1v[jcount][icount],19);
	CHRwsp(fout);
  }
  CHRwcr(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/D2T2V");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,3,5);
  CHRwsp(fout);
  INTw(fout,ydi->nicoup,5);
  CHRwsp(fout);
  CHRwcr(fout);

  for(icount=0;icount<ydi->nicoup;icount++)
  {for(jcount=0;jcount<3;jcount++)
  {
    DBLw(fout,ydi->d2t2v[jcount][icount],19);
	CHRwsp(fout);
  }
  CHRwcr(fout);
  }
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
  INTw(fout,0,5);
  CHRwcr(fout);
  /*for(icount=0;icount<ydn->nnopo;icount++)
  { for(jcount=0;jcount<ydn->nnodim;jcount++)
    { DBLw(fout,ydn->d2nfc[jcount][icount],19);
      CHRwsp(fout);
    }
    CHRwcr(fout);
  }
  CHRwcr(fout); 
*/
	CHRw(fout,"/YD/YDN/D2NFD");
	CHRwsp(fout);
	INTw(fout,21,5);
	CHRwsp(fout);
	INTw(fout,ydn->nnodim,5);
	CHRwsp(fout);
	INTw(fout,0,5);
	CHRwcr(fout);
	
  CHRw(fout,"/YD/YDN/D2NFT");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnodim,5);
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);
  /*for(icount=0;icount<ydn->nnopo;icount++)
  { for(jcount=0;jcount<ydn->nnodim;jcount++)
    { DBLw(fout,ydn->d2nfc[jcount][icount],19);
      CHRwsp(fout);
    }
    CHRwcr(fout);
  }
  CHRwcr(fout); 
*/

  CHRw(fout,"/YD/YDN/D1NMCT");
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { DBLw(fout,ydn->d1nmct[icount],19);
    CHRwsp(fout);
  }
  CHRwcr(fout);
	
	CHRw(fout,"/YD/YDN/D1NVCT");
	CHRwsp(fout);
	INTw(fout,ydn->nnopo,5);
	CHRwcr(fout);
	for(icount=0;icount<ydn->nnopo;icount++)
	{ DBLw(fout,ydn->d1nvct[icount],19);
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

//-PY changed it for 3D_fracture_coupling_with_multiphase
  CHRw(fout,"/YD/YDN/D2NVT");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnodim,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { for(jcount=0;jcount<ydn->nnodim;jcount++)
    { DBLw(fout,ydn->d2nvt[jcount][icount],19);
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
}

/****************** Nodes ******************************************/
static void Ywgid(yde,ydn,fout)
  YDE yde;YDN ydn; FILE *fout;
{ INT icount,jcount;
 
 
  CHRw(fout,"MESH    dimension 3 ElemType Tetrahedra  Nnode ");
  if((yde->nelno==11)||(yde->nelno==5))
  {INTw(fout,yde->nelno-1, 4);}
  else 
  {INTw(fout,yde->nelno, 4);}
  CHRwsp(fout);
  CHRwcr(fout);
  CHRw(fout,"Coordinates ");
  CHRwcr(fout);

  for(icount=0;icount<ydn->nnopo;icount++)
  {
   INTw(fout,icount+1,5);
   CHRwsp(fout);
   DBLw(fout,ydn->d2ncc[0][icount],19);
   CHRwsp(fout);
   DBLw(fout,ydn->d2ncc[1][icount],19);
   CHRwsp(fout);
   DBLw(fout,ydn->d2ncc[2][icount],19);
   CHRwsp(fout);
   CHRwcr(fout);
  }

  CHRw(fout,"end coordinates ");
  CHRwcr(fout);
  CHRw(fout,"Elements ");
  CHRwcr(fout);

  for(icount=0;icount<yde->nelem;icount++)
  {
   INTw(fout,icount+1,5);
   CHRwsp(fout);
   for(jcount=0;jcount<(yde->nelno);jcount++)
   {
   INTw(fout,yde->i2elto[icount][jcount]+1,10);
   CHRwsp(fout);
   }
//   INTw(fout,yde->i2elto[icount][yde->nelno-1],10);
//   CHRwsp(fout);
   CHRwcr(fout);
  }
  CHRw(fout,"end elements ");
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
}

/****************** Boundary Conditions **************************************/
static void Ywdb(ydb,fout)
  YDB ydb; FILE *fout;
{ INT icount;

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
//-PY changed it for 3D_fracture_coupling_with_multiphase
  CHRw(fout,"/YD/YDB/D1BNAP");
  CHRwsp(fout);
  INTw(fout,ydb->nbcon,5);
  CHRwcr(fout);
  for(icount=0;icount<ydb->nbcon;icount++)
  { DBLw(fout,ydb->d1bnap[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);
}

static void Ywdx(ydx,fout)
YDX ydx; FILE *fout;
{ INT icount;
	INT jcount;
	
	CHRw(fout,"     /*   Surface Elements     */");
	CHRwcr(fout);
	CHRw(fout,"/YD/YDX/NELEM");
	CHRwsp(fout);
	INTw(fout,ydx->nsfel,5);
	CHRwcr(fout);
	
	CHRw(fout,"/YD/YDX/NELNO");
	CHRwsp(fout);
	INTw(fout,ydx->nelno,5);
	CHRwcr(fout);
	
	
	CHRw(fout,"/YD/YDX/STRIELE");
	CHRwsp(fout);
	INTw(fout,21,5);
	CHRwsp(fout); 
	INTw(fout,ydx->nelno,5);
	CHRwsp(fout);
	INTw(fout,ydx->nsfel,5);
	CHRwcr(fout);
	
	for(icount=0;icount<ydx->nsfel;icount++)
	{ for(jcount=0;jcount<ydx->nelno;jcount++)
    { INTw(fout,ydx->i2xbto[icount][jcount],6);
		CHRwsp(fout);
    } 
		CHRwcr(fout);
	}
	
}


void Ywrs(CHR *namep, YD yd, INT cp_no)
{ 
CHR name[300];
CHR namef[300];
CHR c1name[100];
CHR c1name1[100];
	CHR cindex[50];
YDC ydc=&(yd->ydc);
YDE yde=&(yd->yde);
YDI ydi=&(yd->ydi);
YDN ydn=&(yd->ydn);
YDO ydo=&(yd->ydo);
YDP ydp=&(yd->ydp);
YDB ydb=&(yd->ydb);
YDX ydx=&(yd->ydx);
FILE *fp=FILENULL;
FILE *fgid=FILENULL;

if(((ydc->ncstep%ydc->isave)==0)||(ydc->ncstep==(ydc->mcstep-1)))   
{

	CHRcpynoext(namef,namep);
	CHRcat(namef,"cp");
	SINTw(cindex,cp_no,0);
	CHRcat(namef,cindex);
	CHRcat(namef,".Y3D");
	fp=fopen(namef,"w");
	setvbuf((fp), NULL, _IOFBF, 1024);
	CHRcpynoext(name,namep);
CHRcat(name,".msh");
fgid=fopen(name,"w");
if((fp!=FILENULL)&&(fgid!=FILENULL))
{

Ywdc(ydc,fp);     /* Write Control variables into the input file     */ 
Ywde(yde,fp);     /* Write Elements variables into the input file    */
Ywdi(ydi,fp);     /* Write Interaction variables into the input file */
Ywdn(ydn,fp);     /* Write Nodes variables into the input file       */
Ywdp(ydp,fp);     /* Write Properties variables into the input file  */
Ywdo(ydo,fp);     /* Write Output variables into the input file      */
Ywdb(ydb,fp);     /* Write Condition variables into the input file   */
Ywdx(ydx,fp);
Ywgid(yde,ydn,fgid);     /* Write GiD mesh file   */

CHRw(fp,"$YDOIT");
CHRwcr(fp);
CHRw(fp,"$YSTOP");
CHRwcr(fp);  
}
fclose(fp);
fclose(fgid);
//fclose(ydc->finp);
//fclose(ydc->fcheck);
}

//return 1;

}	

/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/


