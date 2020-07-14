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

/*1 File Y3Dd.h  Y data base description */

#include "Ytypes.h"
#ifndef FRAMEINCL
#include "frame.h"
#define FRAMEINCL
#endif
#ifndef YTYPESINCL
#include "Ytypes.h"
#define YTYPESINCL
#endif

typedef struct YDC_struct *YDC; 
/*3 YDC is control structure, the variables are as follows */
struct YDC_struct    
{ INT  mcstep;      /*2 mcstep maximum number of time steps   */
  INT  ncstep;      /*2 ncstep current number of time steps   */
  FILE *finp, *fcheck;
  
  DBL  dcgrax;  /*2 dcgray gravity x     */
  DBL  dcgray;  /*2 dcgray gravity y     */
  DBL  dcgraz;  /*2 dcgray gravity z     */

  DBL  dcsizc;  /*2 dcsizc size coord.   */
  DBL  dcsizf;  /*2 dcsizf size force    */
  DBL  dcsizs;  /*2 dcsizs size stress   */
  DBL  dcsizv;  /*2 dcsizv size velocity */
  DBL  dcstec;  /*2 dcstec current time step size                 */
  DBL  dctime;  /*2 dctime current time                           */
  DBL  dcurelx; /*2 dcurelx under relaxation for mass matrix      */

  INT  initer; 
  /*2 ininter the number of interation for multi-pass algorithm   */
  INT  icoutf;  /*2 icoutf write output frequency                 */
  INT  icouti;  /*2 icouti current write output No                */
  INT  icoutp;  /*2 icoutp output precision - digits per number   */
  INT  iwfast;  /*2 fast mode */
  INT  isave;  /*2 restart number of loop */
  DBL  dcrmpt;

}; 

typedef struct YDE_struct *YDE;  
/*3 YDE element description structure, the variables are as follows */
struct YDE_struct
{ INT melem, nelem; 
  /*2 melem, nelem maximum (actual) number of elements              */
  INT melst, nelst; 
  /*2 melst, nelst maximum (actual) number of elemen. states var.   */
  INT melno, nelno; 
  /*2 melno, nelno maximum (actual) number of elemen. nodes         */
  INT   *i1elcf;     /*2 i1elcf [melem] contacting couple first     */ 
  INT   *i1elpr;     /*2 i1elpr [melem]    element property         */
  DBL   **d2elst;    /*2 d2elst [melst][melem]   - element state    */
  INT   **i2elto;    /*2 i2elto [melno][melem]   - element topology */
  DBL   ***d3tcs;    /*2 d3tcs [ndime][ndime][melem] -Cauchy stress */
  DBL   *d1emct;     /*2 d1emct[melem]  total elemental mass        */
  INT   *i1elbe;     /*2 i11elbe [melem] boundary element           */

};

typedef struct YDI_struct *YDI; 
/*3 YDI is interaction structure, the variables are as follows      */
struct YDI_struct
{ INT micoup, nicoup;
  /*2 micoup, nicoup maximum possible number of contacting couples  */
  INT    iiecff;     
  /*2 iiecff interaction element contact. couple free first         */
  DBL    diedi;      /*2 diedi travel since last detection          */
  DBL    diezon;     /*2 diezon buffer zone size                    */
  DBL   *d1iesl;     /*2 d1iesl [mcoup] contact sliding             */

  INT   *i1iecn;     /*2 i1iecn [mcoup] couple next                 */
  INT   *i1iect;     /*2 i1iect [mcoup] couple target               */

  DBL   *d1deltat1;   /*[mcontn] tangential overlap between two elements */
  DBL   *d1deltat2;   /*[mcontn] tangential overlap between two elements */
  DBL   *d1deltan;   /*[mcontn] normal overlap between two elements */
  DBL   **d2nv;
  DBL   **d2t1v;
  DBL   **d2t2v;

}; 
               
typedef struct YDN_struct *YDN;   
 /*3 YDN node description structure, the variables are as follows   */
struct YDN_struct
{ INT mnodim, nnodim;
  /*2 mnodim, nnodim max(actual) nodal dimensions number            */
  INT mnopo, nnopo;  
  /*2 mnopo, nnopo maximum (actual) number of nodal points          */
  DBL  *d1nmct; /*2 d1nmct [mnopo] nodal mass current translation   */
  DBL  **d2ncc; /*2 d2ncc [mnodim][mnopo] nodal coordinate current  */
  DBL  **d2nci; /*2 d2nci [mnodim][mnopo] nodal coordinate initial  */
  DBL  **d2nfc; 
  /*2 d2nfc [mnodim][mnopo] nodal force current due to contact      */
  DBL  **d2nft; /*2 d2nft [mnodim][mnopo] nodal total current force */
  DBL  **d2nvc; /*2 d2nvc [mnodim][mnopo] nodal velocity current    */
  DBL *d1nti; /*!< [mnopo] nodal temperature initial                */ 
  DBL *d1ntc; /*!< [mnopo] nodal temperature current                */
  DBL *d1ntp; /*!< [mnopo] nodal temperature petsc                  */

  INT  *i1nobf; /*2 i1nobf [mnopo] nodal boundary >0 is boundary    */
  INT  *i1nopr; /*2 i1nopr [mnopo] nodal property                   */
};
  
typedef struct YDO_struct *YDO;   
/*3 YDO output description structure, the variables are as follows  */
struct YDO_struct
{ INT mohys, nohys; 
  /*2 mohys, nohys maximum (actual) number of hystory variables     */
  DBL     dohyp;  /*2 dohyp output hystory accuracy                 */ 
  DBL   *d1ohyf;  
  /*2 d1ohyf [mohys] output hystory factor to scale state           */
  DBL   *d1ohyc;  
  /*2 1ohyc [mohys] output hystory factor to scale time             */
  DBL   *d1ohys;  /*2 d1ohys [mohys] output hystory state           */
  DBL   *d1ohyt;  /*2 d1ohyt [mohys] output hystory time            */
  DBL   *d1ohyx;  
  /*2 d1ohyx [mohys] output history x coordinate of the point       */
  DBL   *d1ohyy;  
  /*2 d1ohyy [mohys] output history y coordinate of the point       */
  DBL   *d1ohyz;  
  /*2 d1ohyz [mohys] output history z coordinate of the point       */
  FILE  **f2ohyf;  /*2 f2ohyf [mohys] output history files          */
  INT   *i1ohyt;  
  /*2 i1ohyt [mohys] output hystory type, i.e. which variable       */
};

typedef struct YDP_struct *YDP;  
/*3 YDP property description structure, the variables are as follows */
struct YDP_struct
{ INT mprop, nprop; 
  /*2 mprop, nprop maximum (actual) number of properties             */
  DBL   *d1pefr;  
  /*2 d1pefr [mprop] property sliding friction coefficient           */
  DBL   *d1pesf;  
  /*2 d1pesf [mprop] property static friction coefficient            */
  DBL   *d1pepsf;  
  /*2 d1pepsf [mprop] property pressure friction coefficient         */
  DBL   *d1pevf;  
  /*2 d1pefr [mprop] property velcocity dependent friction sfactor   */
  DBL   *d1pepf;  
  /*2 d1pefr [mprop] property pressure dependent friction sfactor   */

  DBL   *d1peks;  
  /*2 d1peks [mprop] dpeks=2hbeta*sqrt(E*ro) in 2D or 3D,0<beta<1    */
  DBL   *d1pela;  
  /*2 d1pela [mprop] property lamda - Lame elastic constant          */
  DBL   *d1pemu;  
  /*2 d1pemu [mprop] property mu    - Lame elastic constant          */
  DBL   *d1pepe;  /*2 d1pepe [mprop] property penalty parameter      */
  DBL   *d1pero;  /*2 d1pero [mprop] property ro    - density        */
  INT   *i1ptyp;  /*2 i1ptyp [mprop] property type                   */

  DBL   *d1pegfs;
  DBL   *d1pegfn;

  DBL   *d1peft;  
//  DBL   *d1pefs;
  INT   *i1pejp;
  INT   *i1pemn;
  INT   *i1psde;  

  DBL   *d1picf;  /*[mprop] internal friction angle                      */
  DBL   *d1pcoh;  /*[mprop] cohesion                        */
  DBL   *d1tcon;  /*[mprop] Thermal Conductivity                    */
  DBL   *d1capa;  /*[mprop] Thermal Capacity                        */
  DBL   *d1ctex;  /*[mprop] Coefficient of thermal expansion              */
};

typedef struct YDB_struct *YDB;   
/*3 YDB boundary description structure, the variables are as follows */
struct YDB_struct
{ INT mbcon, nbcon; 
  /*2 mprop, nprop maximum (actual) number of properties             */
  DBL   *d1bnax;  /*2 d1pnax [mprop] amplitude of acceleration x     */
  DBL   *d1bnay;  /*2 d1pnay [mprop] amplitude of acceleration y     */
  DBL   *d1bnaz;  /*2 d1pnaz [mprop] amplitude of acceleration z     */

  DBL   *d1bnfx;  /*2 d1pnax [mprop] amplitude of force x            */
  DBL   *d1bnfy;  /*2 d1pnay [mprop] amplitude of force y            */
  DBL   *d1bnfz;  /*2 d1pnaz [mprop] amplitude of force z            */

  DBL   *d1bnvx;  /*2 d1pnax [mprop] amplitude of velocity x         */
  DBL   *d1bnvy;  /*2 d1pnay [mprop] amplitude of velocity y         */
  DBL   *d1bnvz;  /*2 d1pnaz [mprop] amplitude of velocity z         */

  DBL   *d1bntp;  /* d1bntp [mprop] constrained temperature         */
  INT   *i1bntp;  /* i1bntp [mprop] constrained temperature boolean  */
  DBL   *d1bnhf;  /*!< [mprop] amplitude of nodal heat flux    */
  INT   *i1bnhf;  /* i1bntp [mprop] nodal heat flux boolean  */
  DBL   *d1bcvt;  /*!< [mprop] Boundary temperature for convection condition  */
  DBL   *d1bcvc;  /*!< [mprop] Convection coefficient  */ 
  INT   *i1bcvt;  /*2 i1bntp [mprop] Convection condition boolean  */


  INT   *i1bnvx;  
  /*2 i1pnfx [mprop] fixity x direction >0 vel.       */
  INT   *i1bnvy;  
  /*2 i1pnfy [mprop] fixity y direction >0 vel.       */
  INT   *i1bnvz;  
  /*2 i1pnfz [mprop] fixity z direction >0 vel.       */  

};

/**********************************************************************/
/**********************************************************************/

typedef struct YDX_struct *YDX;

struct YDX_struct
{  INT nelem;
   INT melem;
   INT nelno;
   INT **i2elto;
   DBL dthick;
};

/**********************************************************************/
/**********************************************************************/

typedef struct YD_struct *YD;
struct YD_struct
{ struct YDC_struct ydc;   
  struct YDE_struct yde;   
  struct YDI_struct ydi;   
  struct YDN_struct ydn;   
  struct YDO_struct ydo;   
  struct YDP_struct ydp;   
  struct YDB_struct ydb;   
  struct YDX_struct ydx;   
};
   
