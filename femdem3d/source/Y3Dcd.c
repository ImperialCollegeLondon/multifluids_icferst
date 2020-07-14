/*! \file Y3Dcd.c
 *  \brief Y contact detection
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
#include "Ycd3TETintersection.h"


/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/


/* find maximum squared velocity (3D)
 */
 static DBL Ycdvel3TET( 
      INT nelem,
      INT iprop,
      DBL *d1nvcx, DBL *d1nvcy, DBL *d1nvcz,
      INT *i1elpr,
      INT **i2elto
      ) 
 { 
      INT ielem, i, ip, *i1elto;
      DBL vx, vy, vz, v2;
      DBL mxv2 = R0;

      for(ielem=0; ielem<nelem; ielem++)
      { 
	   if(i1elpr[ielem] != iprop)
		continue;

	   i1elto = i2elto[ielem];

	   for(i=0; i<NNODEX; i++)
	   { 
		ip = i1elto[i]; 

		vx = d1nvcx[ip];
		vy = d1nvcy[ip];
		vz = d1nvcz[ip];

		v2 = SQR(vx) + SQR(vy) + SQR(vz);
		if(v2 > mxv2)
		{
		     mxv2 = v2;
		}
	   }
      }

      return mxv2;
 }


/**********************************************************************/
/**********************************************************************/


/* find element coordinates (3D)
 */
static void Ycdcor3TET(  
     INT nelem,
     INT iprop,
     DBL *d1erad,
     DPT *dp1ecc,
     DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,
     INT *i1elbe,
     INT *i1elpr,
     INT **i2elto,
     INT iwfast
     ) 
{ 
     INT ielem, *i1elto;
     INT ip0, ip1, ip2, ip3;
     DBL x0, x1, x2, x3;
     DBL y0, y1, y2, y3;
     DBL z0, z1, z2, z3;
     DBL dx0, dx1, dx2, dx3;
     DBL dy0, dy1, dy2, dy3;
     DBL dz0, dz1, dz2, dz3;
     DBL xc, yc, zc;
     DBL d2, mxd2;

     for(ielem=0; ielem<nelem; ielem++)
     { 
       if((i1elpr[ielem] != iprop) || ((i1elbe[ielem] <= 0) && (iwfast != 0)))
	 {
	   continue;
	 }

	  i1elto = i2elto[ielem];

	  ip0 = i1elto[0]; 
	  ip1 = i1elto[1]; 
	  ip2 = i1elto[2]; 
	  ip3 = i1elto[3]; 

	  x0 = d1nccx[ip0];
	  x1 = d1nccx[ip1];
	  x2 = d1nccx[ip2];
	  x3 = d1nccx[ip3];
	  y0 = d1nccy[ip0];
	  y1 = d1nccy[ip1];
	  y2 = d1nccy[ip2];
	  y3 = d1nccy[ip3];
	  z0 = d1nccz[ip0];
	  z1 = d1nccz[ip1];
	  z2 = d1nccz[ip2];
	  z3 = d1nccz[ip3];

	  xc = (x0 + x1 + x2 + x3)/R4;
	  yc = (y0 + y1 + y2 + y3)/R4;
	  zc = (z0 + z1 + z2 + z3)/R4;

	  dp1ecc[ielem][0] = xc;
	  dp1ecc[ielem][1] = yc;
	  dp1ecc[ielem][2] = zc;

	  /*Z use temp variables to avoid double-calculation within MAXIM() */
	  dx0 = xc - x0;
	  dx1 = xc - x1;
	  dx2 = xc - x2;
	  dx3 = xc - x3;
	  dy0 = yc - y0;
	  dy1 = yc - y1;
	  dy2 = yc - y2;
	  dy3 = yc - y3;
	  dz0 = zc - z0;
	  dz1 = zc - z1;
	  dz2 = zc - z2;
	  dz3 = zc - z3;

	  /*Z faster than using MAXIM() */
	  mxd2 = SQR(dx0) + SQR(dy0) + SQR(dz0);

	  d2 = SQR(dx1) + SQR(dy1) + SQR(dz1);
	  if(d2 > mxd2)
	  {
	       mxd2 = d2;
	  }

	  d2 = SQR(dx2) + SQR(dy2) + SQR(dz2);
	  if(d2 > mxd2)
	  {
	       mxd2 = d2;
	  }

	  d2 = SQR(dx3) + SQR(dy3) + SQR(dz3);
	  if(d2 > mxd2)
	  {
	       mxd2 = d2;
	  }

	  d1erad[ielem] = SQRT(mxd2);
     }
}


/**********************************************************************/
/**********************************************************************/


/* Y elements contact detection (Z - this can be split up using INLINE)
 */
static void Ycdprocess(
     INT melem,   INT micoup,  INT nelem,    INT nprop,
     DBL dcstec,  DBL diezon,  DBL *d0iedi,  INT *i0iecff,
     INT *nicoup,
     DBL *d1nccx, DBL *d1nccy, DBL *d1nccz, 
     DBL *d1nvcx, DBL *d1nvcy, DBL *d1nvcz, 
     INT *i1elcf, INT *i1elpr, INT *i1iecn,  INT *i1iect,
     INT *i1elbe, INT *i1ptyp, INT **i2elto, INT iwfast,
     INT **i2elfr, INT **i2eljp, INT ncstep, INT *i1elcft,
     INT *i1fcstep
     )
{ 
     INT ncelx,ncely,ncelz;  /* total number of x, y, z cells   */
     INT nelemd;             /* twice total number of elements  */
     DBL diam;               /* maximum diameter of element     */
     DBL dmxv2;              /* maximum displacement increment  */
     INT icoup,jcoup;        /* couple                          */
     INT icouco;             /* couple's contactor              */
     INT icouta;             /* couple's target                 */
     INT ielem;                /* element                       */
     INT ielemx;               /* element assigned to cell x    */
     INT ielemy;               /* element assigned to cell y    */
     INT ielemz;               /* element assigned to cell z    */
     INT ihx,ihy,ihz;          /* x, y, z  head of a   list     */
     INT iminx,iminy,iminz;    /* space boundaries              */
     INT imaxx,imaxy,imaxz;    /* space boundaries              */
     INT iprop;                /* element property id           */  
     INT ix,iy,iz;             /* x, y, z cell                  */
     INT jelemx;               /* element assigned to cell x    */
     DBL *d1erad;     /* element radius squared                 */
     DPT *dp1ecc;     /* element coordinate current xyz         */
     INT *i1eccx;     /* element coordinate current ix          */
     INT *i1eccy;     /* element coordinate current iy          */
     INT *i1eccz;     /* element coordinate current iz          */
     INT *i1cnx;      /* contactor next  x                      */
     INT *i1cny;      /* contactor next  y                      */
     INT *i1cnz;      /* contactor next  z                      */
     INT *i1cfz;      /* contactor first z                      */
     INT  i1heax[5];  /* heads of 5 connected lists for x cells */
     INT  i1heay[5];  /* heads of 5 connected lists for y cells */
     INT  i1heaz[2];  /* heads of 2 connected lists for z cells */
     INT *i2cfx[2];   /* contactor first x                      */
     INT *i2cfy[2];   /* contactor first y                      */
     DBL dx, dy, dz, d2, dr,tmp;
     INT index;
     INT flag,isurf,jsurf;

     nelemd = nelem*2;

     /* initialise data if needed i.e. if(i0iecff==default) */
     if(*i0iecff == -2)
     { 
	  for(ielem=0; ielem<melem; ielem++)
	  { 
	       i1elcf[ielem] = -1;
	  }
	  for(icoup=0; icoup<micoup; icoup++)
	  { 
	       i1iecn[icoup] = icoup + 1;
	  }
	  i1iecn[micoup - 1] = -1;
	  *i0iecff = 0;
     }

     /* find maximum  velocity */
     dmxv2 = R0;
     for(iprop=0; iprop<nprop; iprop++)   
     { 
	  if((i1ptyp[iprop] == YTE3TET4ELS )||
	     (i1ptyp[iprop] == YTE3TET10ELS))
	  { 
	       /*Z can dmxv2 be assigned more than once? */
	       tmp = Ycdvel3TET( 
		    nelem,
		    iprop,
		    d1nvcx, d1nvcy, d1nvcz,
		    i1elpr,
		    i2elto 
		    );
	       if(tmp>dmxv2)
	       {
		    dmxv2=tmp;
	       }
	  }
     }
     *d0iedi += R2*dcstec*SQRT(dmxv2); 

     /*Z split function here? */

     /* if not time for contact detection, return */
     if(*d0iedi < diezon)
	  return;

     /* unmark new couples */
     for(icoup=0; icoup<micoup; icoup++)
     { 
	  if(i1iect[icoup] < 0)
	  {
	       i1iect[icoup] = -1 - i1iect[icoup];
	  }
     }

     /* find element coordinates */
     d1erad = TalzDBL1(nelem);    /* element radius squared         */
     dp1ecc = TalzDPT1(nelem);    /* element coordinate current xyz */
     i1eccx = TalINT1(nelem);     /* element coordinate current ix  */
     i1eccy = TalINT1(nelem);     /* element coordinate current iy  */
     i1eccz = TalINT1(nelem);     /* element coordinate current iz  */
     *d0iedi = R0;

     //printf("******Ycdcor3TET******\n");
     for(iprop=0; iprop<nprop; iprop++)
     { 
	  if((i1ptyp[iprop] == YTE3TET4ELS )||
	     (i1ptyp[iprop] == YTE3TET10ELS))
	  { 
	       Ycdcor3TET( 
		    nelem ,
		    iprop ,
		    d1erad,
		    dp1ecc,
		    d1nccx, d1nccy, d1nccz,
		    i1elbe,
		    i1elpr,
		    i2elto,
		    iwfast
		    );
	  }	  
     }

     /* find maximum diameter */ /*Z can this not be constant for simulation? */
     diam = d1erad[0];
     for(ielem=1; ielem<nelem; ielem++)
     { 
	  if(d1erad[ielem] > diam)
	  {
	       diam = d1erad[ielem];
	  }
     }
     if(diam < EPSILON)
     {
	  diam = BEPSILON;
     }
     diam = R2*(diam + diezon);

     /* integerise element coordinates, find space boundaries */
     iminx = INT_MAX;
     iminy = INT_MAX;
     iminz = INT_MAX;
     imaxx = INT_MIN;
     imaxy = INT_MIN;
     imaxz = INT_MIN;

     for(ielem=0; ielem<nelem; ielem++)
     { 
	  i1eccx[ielem] = (INT)(dp1ecc[ielem][0]/diam);
	  i1eccy[ielem] = (INT)(dp1ecc[ielem][1]/diam);
	  i1eccz[ielem] = (INT)(dp1ecc[ielem][2]/diam); 

	  if(i1eccx[ielem] < iminx)
	  {
	       iminx = i1eccx[ielem];
	  }
	  if(i1eccx[ielem] > imaxx)
	  {
	       imaxx = i1eccx[ielem];
	  }

	  if(i1eccy[ielem] < iminy)
	  {
	       iminy = i1eccy[ielem];
	  }
	  if(i1eccy[ielem] > imaxy)
	  {
	       imaxy = i1eccy[ielem];
	  }

	  if(i1eccz[ielem] < iminz)
	  {
	       iminz = i1eccz[ielem];
	  }
	  if(i1eccz[ielem] > imaxz)
	  {
	       imaxz = i1eccz[ielem];
	  }
     }

     iminx -= 1;
     iminy -= 1;
     iminz -= 1;
     imaxx += 2;
     imaxy += 2;
     imaxz += 1;

     /* normalise coordinates */
     for(ielem=0; ielem<nelem; ielem++)
     { 
	  i1eccx[ielem] -= iminx;
	  i1eccy[ielem] -= iminy;
	  i1eccz[ielem] -= iminz;
     }

     /* calculate cells, allocate memory  */
     ncelx = imaxx - iminx;
     ncely = imaxy - iminy;
     ncelz = imaxz - iminz;
     i1cnx    = TalINT1(nelem); /* contactor next x                    */
     i1cny    = TalINT1(nelem); /* contactor next y                    */
     i1cnz    = TalINT1(nelem); /* contactor next z                    */
     i1cfz    = TalINT1(ncelz); /* contactor first z                   */
     i2cfy[0] = TalINT1(ncely); /* contactor first y  (iz-1)           */
     i2cfy[1] = TalINT1(ncely); /* contactor first y  (iz  )           */
     i2cfx[0] = TalINT1(ncelx); /* contactor first x  (iz-1,iy-1)      */
     i2cfx[1] = TalINT1(ncelx); /* contactor first x  (iz-1,iy  )      */

     /* assume no contactors at any cell */
     for(iz=0; iz<ncelz; iz++)
     { 
	  i1cfz[iz] = -1;  
     }
     for(iy=0; iy<ncely; iy++)
     { 
	  i2cfy[0][iy] = -1;
	  i2cfy[1][iy] = -1;
     }
     for(ix=0; ix<ncelx; ix++)
     { 
	  i2cfx[0][ix] = -1;
	  i2cfx[1][ix] = -1;
     } 

     /* assign all contactors to z-cells - Z can this be put in loop below? */
     for(ielem=0; ielem<nelem; ielem++) 
     { 
	  if(d1erad[ielem] > R0)
	  { 
	       i1cnz[ielem] = i1cfz[i1eccz[ielem]];
	       i1cfz[i1eccz[ielem]] = ielem;
	  } 
     }

     /* scan all loaded z cells */ /* to make parallel, i1heaz[] etc. must be local to loop */
     for(ielem=0; ielem<nelem; ielem++)
     { 
	  iz = i1eccz[ielem];
	  if(i1cfz[iz] < nelem)
	  { 
	       i1heaz[0] = i1cfz[iz];
	       i1heaz[1] = i1cfz[iz-1];
	       if(i1heaz[1] > nelem)
	       {
		    i1heaz[1] -= nelemd;
	       }
	       i1cfz[iz] += nelemd;

	       /* load elements from cells iz & iz-1 onto y cells */
	       for(ihz=0; ihz<2; ihz++)
	       { 
		    ielemz = i1heaz[ihz];
                    
		    while(ielemz >= 0)
		    {
			 i1cny[ielemz] = i2cfy[ihz][i1eccy[ielemz]];
			 i2cfy[ihz][i1eccy[ielemz]] = ielemz;
			 ielemz = i1cnz[ielemz];
		    } 
	       }

	       /* scan all loaded y cells */
	       ielemz = i1heaz[0];
	       while(ielemz >= 0)
	       { 
		    iy = i1eccy[ielemz];
		    if(i2cfy[0][iy] < nelem)
		    { 
			 i1heay[0] = i2cfy[0][iy];
			 i1heay[1] = i2cfy[0][iy-1];
			 i1heay[2] = i2cfy[1][iy+1];
			 i1heay[3] = i2cfy[1][iy];
			 i1heay[4] = i2cfy[1][iy-1];
			 if(i1heay[1] > nelem)
			 {
			      i1heay[1] -= nelemd;
			 }
			 i2cfy[0][iy] += nelemd;

			 /* load elements from y cells onto x cells */
			 ihx = 0;
			 for(ihy=0; ihy<5; ihy++)
			 { 
			      if(ihy > 0)
			      {
				   ihx = 1;
			      }
			      ielemy = i1heay[ihy];
			      while(ielemy >= 0)
			      { 
				   i1cnx[ielemy] = i2cfx[ihx][i1eccx[ielemy]];
				   i2cfx[ihx][i1eccx[ielemy]] = ielemy;
				   ielemy = i1cny[ielemy];
			      } 
			 }

			 /* scan all loaded x cells */
			 ielemy = i1heay[0];
			 while(ielemy >= 0)
			 { 
			      ix = i1eccx[ielemy];
			      if(i2cfx[0][ix] < nelem)
			      { 
				   i1heax[0] = i2cfx[0][ix];
				   i1heax[1] = i2cfx[0][ix-1];
				   i1heax[2] = i2cfx[1][ix+1];
				   i1heax[3] = i2cfx[1][ix];
				   i1heax[4] = i2cfx[1][ix-1];
				   if(i1heax[1] > nelem)
				   {
					i1heax[1] -= nelemd;
				   }
				   i2cfx[0][ix] += nelemd;

				   /* detect contacts for cell (ix,iy,iz) */
				   ielemx = i1heax[0];
				   while(ielemx >= 0)
				   { 
					for(ihx=0; ihx<5; ihx++)
					{ 
					     jelemx = i1heax[ihx];
					     while(jelemx >= 0)
					     {
					       flag=0;
					       if((i1elbe[ielemx]==YIPROPMAX)&&(i1elbe[jelemx]==YIPROPMAX))
						 {
						   flag=1;

						   for(isurf=0;isurf<4;isurf++)
						     {
						       if(i2eljp[ielemx][isurf]>0)
							 {
							   for(jsurf=0;jsurf<4;jsurf++)
							     {
							       if(i2eljp[ielemx][isurf]==i2eljp[jelemx][jsurf])
								 {
								   flag=0;
								   break;	    
								 }
							     }
							   if(flag==0)
							     {
							       break;
							     }
							 }
						     }
						 }

					       if(((ihx!=0) || (ielemx>jelemx)) && ((i1elbe[ielemx]!=i1elbe[jelemx]) || (iwfast==0) || (flag==1)))
						  { 
						       icouco = MAXIM(ielemx, jelemx);
						       icouta = MINIM(ielemx, jelemx);

						       /* exclude old couple */
						       icoup = i1elcf[icouco];
						       jcoup = icoup;
						       while((icoup>=0) && (icouta>=0))
						       { 
							    if(i1iect[icoup] == icouta) 
							    { 
								 dx = dp1ecc[icouco][0] - dp1ecc[icouta][0];
								 dy = dp1ecc[icouco][1] - dp1ecc[icouta][1];
								 dz = dp1ecc[icouco][2] - dp1ecc[icouta][2];
								 d2 = SQR(dx) + SQR(dy) + SQR(dz);
								 dr = d1erad[icouco] + d1erad[icouta] + diezon;
								 if(d2 > SQR(dr))
								 { 
								      i1iecn[jcoup] = i1iecn[icoup];
								      if(jcoup == icoup)
								      {
									   i1elcf[icouco] = i1iecn[icoup];
								      }
								      i1iecn[icoup] = *i0iecff;
								      *i0iecff = icoup;

								      /* added by LGuo */
								      i1fcstep[icoup]=-1;
								 }
								 icouta = -1;
							    }
							    jcoup = icoup;
							    icoup = i1iecn[icoup];
						       }

						       /* add new couple if close */
						       if(icouta >= 0)
						       { 
							    dx = dp1ecc[icouco][0] - dp1ecc[icouta][0];
							    dy = dp1ecc[icouco][1] - dp1ecc[icouta][1];
							    dz = dp1ecc[icouco][2] - dp1ecc[icouta][2];
							    d2 = SQR(dx) + SQR(dy) + SQR(dz);
							    dr = d1erad[icouco] + d1erad[icouta] + diezon;
							    if(d2 < SQR(dr))
							    { 
								 index=Ycd3TETintersection(d1nccx,d1nccy,d1nccz,icouco,icouta,i2elto,diezon);
								      if(index==1)
								      {
									index=Ycd3TETintersection(d1nccx,d1nccy,d1nccz,icouta,icouco,i2elto,diezon); /* swap icouco & icouta */
										}
								      if(index==1)
								      {
									   icoup = *i0iecff;
									   if(icoup < 0) /*Z ever happen? - LG Yes*/
									   { 
										CHRw(stderr,"Ycd: to small meicoc");
										CHRwcr(stderr); exit(1);
									   }
									   *i0iecff = i1iecn[icoup];
									   i1iecn[icoup] = i1elcf[icouco];
									   i1elcf[icouco] = icoup;
									   i1iect[icoup] = -1 - icouta;/* mark new */
									   if((*i0iecff)>(*nicoup))
									   {
										(*nicoup)=(*i0iecff);
									   }

									   /* added by LGuo */
									   if((i1elcft[icouco]==ncstep)&&(i1elcft[icouta]==ncstep))
									     {
									       i1fcstep[icoup]=ncstep;
									     }
								      }
							    } 
						       } 
						  }
						  jelemx = i1cnx[jelemx];
					     } 
					}
					ielemx = i1cnx[ielemx];
				   } 
			      }
			      ielemy = i1cny[ielemy];
			 }

			 /* unload elements from x cells */
			 ihx = 0;
			 for(ihy=0; ihy<5; ihy++)
			 { 
			      if(ihy > 0) ihx = 1;
			      ielemy = i1heay[ihy];
			      while(ielemy >= 0)
			      { 
				   i2cfx[ihx][i1eccx[ielemy]] = -1;
				   ielemy = i1cny[ielemy];
			      }
			 }
		    }
		    ielemz = i1cnz[ielemz];
	       }

	       /* unload elements from y cells */
	       for(ihz=0; ihz<2; ihz++)
	       { 
		    ielemz = i1heaz[ihz];
		    while(ielemz >= 0)
		    { 
			 i2cfy[ihz][i1eccy[ielemz]] = -1;
			 ielemz = i1cnz[ielemz];
		    } 
	       }
	  } 
     }

     for(icoup=0;icoup<micoup;icoup++)
       {
	 if(ncstep>(i1fcstep[icoup]+1000))
	   {
	     i1fcstep[icoup]=-1;
	   }
       }

     /*Z free memory (note: in reverse order of allocation) */
     FREE(i2cfx[1]);
     FREE(i2cfx[0]);
     FREE(i2cfy[1]);  
     FREE(i2cfy[0]);
     FREE(i1cfz);
     FREE(i1cnz);
     FREE(i1cny);
     FREE(i1cnx);
     FREE(i1eccz);
     FREE(i1eccy);
     FREE(i1eccx);
     FREE(dp1ecc);
     FREE(d1erad);
}


/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/


/*! \brief contact detection
 *  \param[in]     ydc control database
 *  \param[in]     yde element database
 *  \param[in,out] ydi interaction database
 *  \param[in]     ydn nodal database
 *  \param[in]     ydp property database
 *  \par Details:
 *  Ycd() ... chapter 3 (Munjiza NBS)
 *
 \verbatim
 initialise data if needed
 find maximum  velocity
 if not time for contact detection, return
 find element coordinates - Ycdcor3TET()
 unmark new couples
 find maximum diameter
 intigerise element coordinates
 find space boundaries
 normalise coordinates
 assume no contactors at any cell
 assign all contactors to z-cells
 scan all loaded z cells (for each element)
 load elements from cells iz & iz-1  onto y cells
 scan all loaded y cells
 load elements from y cells onto x cells
 scan all loaded x cells
 detect contacts for cell (ix,iy,iz)
 for(ihx=0;ihx<5;ihx++)
 exclude old couple
 add new couple if close
 unload elements from x cells
 unload elements from y cells
 \endverbatim
 *
 */
 void Ycd(YDC ydc, YDE yde, YDI ydi, YDN ydn, YDP ydp)
 { 
      /* Process contact detection if specified by input data */
      if(ydi->micoup > 0)
      {
	Ycdprocess( /* Y elements contact detection */
		   yde->melem,    ydi->micoup,   yde->nelem,     ydp->nprop,
		   ydc->dcstec,   ydi->diezon,  &ydi->diedi,    &ydi->iiecff,
		   &ydi->nicoup,
		   ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2ncc[2],
		   ydn->d2nvc[0], ydn->d2nvc[1], ydn->d2nvc[2],
		   yde->i1elcf,   yde->i1elpr,   ydi->i1iecn,    ydi->i1iect, 
		   yde->i1elbe,   ydp->i1ptyp,   yde->i2elto,    ydc->iwfast,
		   yde->i2elfr,   yde->i2eljp,   ydc->ncstep,    yde->i1elcft,
		   ydi->i1fcstep
		    );
      }
 }


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/


