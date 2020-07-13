/*! \file Yproto.h
 *  \brief Y main prototypes
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


#ifndef YPROTO_H
#define YPROTO_H


#include "petscksp.h"
#include "Y3Dd.h"


/**********************************************************************/
/**********************************************************************/

void Ymd(        /* mesh elements, added by LG */
 YDE yde,        /* element database */
 YDI ydi,        /* interaction database */
 YDN ydn,        /* nodal database */
 YDP ydp         /* property database */
); 

void Ycd(        /*   contact detection                               */
 YDC ydc,         /* = control database                               */
 YDE yde,         /* = element database                               */
 YDI ydi,         /* = interaction database                           */
 YDN ydn,         /* = nodal database                                 */
 YDP ydp         /* - property database                               */
);

void Yfd(         /*   nodal forces                                   */
 YDE yde,         /* = element database                               */
 YDN ydn,         /* = nodal database                                 */
 YDP ydp,         /* - property database                              */
 YDK ydk,         /* = constants database                             */
 YDB ydb,         /* - boundary condition database                    */
 YDI ydi,         /* = interaction database                           */
 YDC ydc,          /* = control database                               */
 YDX ydx
);

void Yid(        /*   procreate                                       */
 YDC ydc,         /* = control database                               */
 YDE yde,         /* = element database                               */
 YDI ydi,         /* = interaction database                           */
 YDN ydn,         /* = nodal database                                 */
 YDP ydp          /* - property database                              */
// YSP ysp
);

void Ysd(        /*  solve equations                                  */
 YDC ydc,         /* - control database                               */
 YDE yde,         /* = element database                               */ 
 YDN ydn,         /* = nodal database                                 */
 YDO ydo,         /* = output database                                */
 YDB ydb,         /* - boundary condition database                    */
 YDK ydk,          /* = constants database                             */
 YDP ydp,          /* = property database                              */
 DBL *absorp,
 DBL *denf,
 INT index,
 DBL *f1dt,
 DBL *f1dt1,
 DBL *f1dt2,
 YPAR ypar,
 YSP ysp
);


//-PY changed it for 3D_fracture_coupling_with_multiphase
INT Yrd(         /*  read input                                       */
 CHR *namep,      /* - name of the problem i.e. input file            */
 CHR *name_grid,CHR *name_mesh,
 YD  yd           /* = database                                       */
);
 
void Yod(         /* output results in space-save format              */
 CHR *namep,      /* - name of the problem i.e. input file            */
 YD  yd           /* = database                                       */
);

void Yinit(       /* constants database initialization                */
 YDK ydk          /* = database                                       */
);

void Yconfig(     /* constants database configuration                 */
 YDK ydk,         /* = database                                       */
 YDC ydc,         /* = control database                               */
 YDB ydb          /* = boundary condition database                    */
);

void Ywrs(         /*  read input                                       */
 CHR *namep,      /* - name of the problem i.e. input file            */
 YD  yd,           /* = database                                       */
 INT cp_no         /* checkpoint no  */
);







void Ysurface(
 YDE yde,
 YDX ydx,
 YDN ydn
);

void YSFtoJOINT(
 YDE yde,
 YDN ydn,
 YDX ydx,
 YDXN ydxn,
 YDP ydp
);

void YSFmesh(
 YDC ydc,
 YDX ydx,
 YDXN ydxn,
 YDN ydn,
 YDP ydp,
 YDE yde
);



//-PY changed it for 3D_fracture_coupling_with_multiphase
void  mirtichRoutine(DBL **points,
                      DBL *intg, 
                      DBL dro);

int Ycd3TETintersection( double *d1nccx, double *d1nccy, double *d1nccz, 
						INT iele, INT jele, INT **i2elto,double buf);
void Ypacking(         /* packing              */
  YD  yd           /* = database                                       */
);

void V3DRot(
     DBL *afx, DBL *afy, DBL *afz, DBL bfx, DBL bfy,
     DBL bfz, DBL u, DBL v, DBL w, DBL cosin, DBL sine
);

void solveSymetricEigenProblem(double **m, double  *v1, double *v2, 
							   double  *v3, double *a1,double *a2,double *a3) ;

void V3DRotAxes(DBL *afx,DBL *afy,DBL *afz,DBL bfx,DBL bfy,DBL bfz,DBL u,DBL v,DBL w);





void Yhpd(
 YDX ydx,
 YDN ydn, 
 YDXN ydxn,
 DBL *d1hystp
);
void Ybnp(
          YDE yde,
          YDN ydn
          );
			 
void Yring(
           YDX ydx,
           YDXN ydxn,
           YDN ydn,
           YSP ysp
           );
void Yring_drag(
                YDN ydn,
                YDX ydx,
                DBL *d1nmu,
                DBL *fx,
                DBL *fy,
                DBL *fz,
                DBL *usl,
                DBL *vsl,
                DBL *wsl,
                DBL *axx,
                DBL *ayy,
                DBL *axy,
                DBL *axz,
                DBL *ayz,
                DBL *azz,
                DBL *d1vx,
                DBL *d1vy,
                DBL *d1vz
                );
void Yscale(
	    YDE yde,
	    YDX ydx,
	    YDN ydn,
	    YDXN ydxn
	    );

void Ythm_contact_init(
YDE yde,
YDN ydn,
YDX ydx
);

void Ythm(
 YDE yde,
 YDN ydn,
 YDP ydp,
 YDK ydk,
 YDC ydc,
 YDB ydb,
 YDX ydx,
 Mat Chi,
 DBL *d1nt_r,
 DBL *d1nt_r_lhs,
 DBL *d1nt_v,
 INT flag_thm,
 YDXN ydxn
);

void Yap(
	 YDE yde,
	 YDN ydn,
	 YDX ydx
	);

void Ypore( YDE yde,
	    YDN ydn,
	    YDX ydx,
	    DBL *d1flpore
	    );



/**********************************************************************/
/**********************************************************************/


#endif /* YPROTO_H */


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/
