/**********************************************************************/
/** Copyright (C) 2008,                                              **/
/** Queen Mary University of London (QMUL) & Imperial College        **/
/** of Science, Technology and Medicine (ICSTM). All rights reserved.**/
/** Implemented for you by Dr Jiansheng Xiang                        **
 
 
 
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
 * Dr J.Xiang Dr J.P.Latham
 * j.xiang@ic.ac.uk or j.p.latham@imperial.ac.uk
 * ******************************************************************* */
/* File  Yod.c */
#include "Yproto.h"
#include "petscksp.h"
static struct YD_struct yd; /* Y database                       */

//static INT *i1l2q;
//static INT *i1q2l;
static INT nnode_l;
static INT flag = 0;

/* continnum surface mesh */
static void Ywgidcontinuum (yde, ydx, ydxn, fout)
  YDE yde; YDX ydx;YDXN ydxn;FILE *fout;
{ INT icount,jcount;

  CHRw(fout, "MESH dimension 3 ElemType Tetrahedra Nnode 4");
  CHRwsp(fout);
  CHRwcr(fout);

  CHRw(fout, "Coordinates ");
  CHRwcr(fout);


  for (icount = 0; icount < ydxn->ncmno; icount++)
    {
      INTw(fout, icount + 1, 5);
      CHRwsp(fout);
      DBLw(fout, ydxn->d2xncc[0][icount], 19);
      CHRwsp(fout);
      DBLw(fout, ydxn->d2xncc[1][icount], 19);
      CHRwsp(fout);
      DBLw(fout, ydxn->d2xncc[2][icount], 19);
      CHRwsp(fout);
      CHRwcr(fout);
    }

  CHRw(fout, "End Coordinates ");
  CHRwcr(fout);
  CHRw(fout, "Elements ");
  CHRwcr(fout);

    printf("ncmel %i \n\n", ydxn->ncmno);
    printf("ncmel %i \n\n", ydx->ncmel);
  for (icount = 0; icount < ydx->ncmel; icount++)
    {
      INTw(fout, icount + 1, 5);
      CHRwsp(fout);
      for (jcount = 0; jcount < 4; jcount++)
	{
	  INTw(fout, ydx->i2cvto[icount][jcount] + 1, 10);
	//	  INTw(fout, yde->i2elto[icount][jcount] + 1, 10);

	  CHRwsp(fout);
	}
      //   INTw(fout,yde->i2elto[icount][yde->nelno-1],10);
      //   CHRwsp(fout);
      CHRwcr(fout);
    }
  CHRw(fout, "End Elements ");
  CHRwcr(fout);

}

static void Ywgidsurface (yde, ydx, ydxn, fout)
  YDE yde; YDX ydx;YDXN ydxn;FILE *fout;
{ INT icount,jcount;

  CHRw(fout, "MESH dimension 2 ElemType Triangle Nnode 3");

  CHRwsp(fout);
  CHRwcr(fout);

  CHRw(fout, "Coordinates ");
  CHRwcr(fout);

  for (icount = 0; icount < ydxn->ncmno; icount++)
    {
      INTw(fout, icount + 1, 5);
      CHRwsp(fout);
      DBLw(fout, ydxn->d2xncc[0][icount], 19);
      CHRwsp(fout);
      DBLw(fout, ydxn->d2xncc[1][icount], 19);
      CHRwsp(fout);
      DBLw(fout, ydxn->d2xncc[2][icount], 19);
      CHRwsp(fout);
      CHRwcr(fout);
    }

  CHRw(fout, "End Coordinates ");
  CHRwcr(fout);
  CHRw(fout, "Elements ");
  CHRwcr(fout);

  for (icount = 0; icount < ydx->nsfel; icount++)
    {
      INTw(fout, icount + 1, 5);
      CHRwsp(fout);
      for (jcount = 0; jcount < 3; jcount++)
	{
	  INTw(fout, ydx->i2xbto[icount][jcount] + 1, 10);
	  CHRwsp(fout);
	}
      //   INTw(fout,yde->i2elto[icount][yde->nelno-1],10);
      //   CHRwsp(fout);
      CHRwcr(fout);
    }
  CHRw(fout, "End Elements ");
  CHRwcr(fout);

}

void y3d_allocate_femdem_ (c1name, 
  nnode_r, nnelem_l_r, nnelem_t_s_r, nnode, nnelem_l, nnelem_t_s,iflag)

  CHR *c1name; //CHR *g1name;CHR *p1name;
  int *nnelem_l_r;int *nnode_r;int *nnelem_t_s_r;int *nnelem_l;int *nnode;int *nnelem_t_s; int *iflag;
{
  extern struct YD_struct yd; /* Y database                      */
  YDC ydc = &(yd.ydc); /* Y control database                  */
  YDE yde = &(yd.yde); /* Y element database                  */
  YDI ydi = &(yd.ydi); /* Y interaction database              */
  YDN ydn = &(yd.ydn); /* Y node database                     */
  YDO ydo = &(yd.ydo); /* Y output database                   */
  YDP ydp = &(yd.ydp); /* Y property database                 */
  YDB ydb = &(yd.ydb); /* Y boundary conditions database      */
  YDX ydx = &(yd.ydx);
  YDK ydk = &(yd.ydk);
  YDXN ydxn = &(yd.ydxn);

  YPAR ypar = &(yd.ypar);
  YSP ysp = &(yd.ysp); /*2 ydb is Y boundary condition database */

  INT i, j, k, index;
//  INT *i1l2q_tmp;
  extern INT nnode_l;

//-PY changed it for 3D_fracture_coupling_with_multiphase
  CHR *p1name;
  CHR *g1name;

//  extern INT *i1l2q;
  extern INT flag;

//  *ifracture=0;  

  nnode_l=0;


//  i1l2q_tmp=INT1NULL;


  printf ("in y3d_allocate \n\n");

  /* get name of the problem */

  ydc->finp = FILENULL;
  ydc->fcheck = FILENULL;
  if (flag == 0)
    {
      Yinit (ydk);

      //if(Yrd(c1name,&yd)>0)
    g1name="toe_zero.txt";
    p1name="box.msh";

      if (Yrd (c1name, g1name, p1name, &yd) > 0)

     printf ("NEW input\n ");
      /*4 calculate consistent mass matrix */
      Yconfig (ydk, ydc, ydb);
      /* added by LGuo for coupling */
      Ysurface (yde, ydx, ydn);

      printf ("yde->nelem1 -------------------- %d\n\n", yde->nelem);
      printf ("ydx->nnopo1 %d\n\n", ydn->nnopo);
      /* added by LGuo */
      Ymd (yde, ydi, ydn, ydp);  // -ao check me


     printf ("yde->nelem2 ----------------- %d\n\n", yde->nelem);
     printf ("ydx->nnopo %d\n\n", ydn->nnopo);


      flag = flag + 1;
      
      for(i=0; i < 3; i++)
       {
         ydc->iflag[i]=iflag[i];
        printf ("ydc->iflag ----------------- %d\n", ydc->iflag[i]);

       }

      /* added by LGuo for coupling */
      YSFtoJOINT(yde,ydn,ydx,ydxn,ydp);
      YSFmesh(ydc,ydx,ydxn,ydn,ydp,yde);

      ydx->nsfel_old = ydx->nsfel;

  //-ao note to JX: this is scaling part
      if(ydc->iflag[1]>0)
      {
      Yap(yde, ydn, ydx);
      Yscale(yde, ydx, ydn, ydxn);      
      }
  //-ao note to JX: this is scaling part

      Yring (ydx,ydxn, ydn, ysp); /*JXiang 3D coupling new  */
    }





  nnode_l = ydxn->ncmno;

  *nnelem_t_s = ydx->nsfel; //cont edges
  *nnelem_l = ydx->ncmel;   //cont eles
  *nnode = nnode_l;         //cont nodes

  *nnode_r = 2 * ydx->nnode; //ring nodes
  *nnelem_l_r = 3 * ydx->nsfel; //ring elements
  *nnelem_t_s_r = 2 * ydx->nsfel; //ring edges

  printf ("nnelem_t_s is %d\n\n", *nnelem_t_s);
  printf ("nnode is %d\n\n", *nnode);
  printf ("nnelem_l is %d\n\n", *nnelem_l);
  printf ("nnelem_l_r is %d %d %d\n\n", *nnode_r, *nnelem_l_r,  *nnelem_t_s_r);


  printf ("out y3d_allocate \n\n");

}

void checkpoint_femdem_ (c1name, cp_no)
  int *cp_no;CHR *c1name;
{int icp_no;
  extern struct YD_struct yd;
  icp_no = *cp_no;
  printf ("in 3D checkpoint_femdem \n\n");
  printf ("name: %s\n", c1name);
  printf ("cp_no is %d\n", icp_no);
//  Ywrs (c1name, &yd, icp_no);

}

void y3d_populate_femdem_(i1elto0,i1elto1,i1elto2,i1elto3,
                          i1elto_rt0,i1elto_rt1,i1elto_rt2,
                          d1nsvx,d1nsvy,d1nsvz,
                          pp1, pp2, pp3, pp4, pp5, pp6, pp7, pp8, pp9,
                          i1elto_v0,i1elto_v1,i1elto_v2,i1elto_v3,
                          i1elto_t_l0,i1elto_t_l1,i1elto_t_l2,
                          d1nx,d1ny,d1nz)

int *i1elto0; int *i1elto1; int *i1elto2; int *i1elto3;
int *i1elto_rt0; int *i1elto_rt1; int *i1elto_rt2;
int *i1elto_v0; int *i1elto_v1; int *i1elto_v2; int *i1elto_v3;
int *i1elto_t_l0; int *i1elto_t_l1; int *i1elto_t_l2;
  DBL *d1nsvx;DBL *d1nsvy;DBL *d1nsvz;

//DBL *d1nsvx_vol; DBL *d1nsvy_vol;DBL *d1nsvz_vol;  /* more thinking, JXiang  */
  DBL *pp1;DBL *pp2;DBL *pp3;DBL *pp4;DBL *pp5;DBL *pp6;DBL *pp7;DBL *pp8;DBL *pp9;
 DBL *d1nx;DBL *d1ny;DBL *d1nz;

{
  extern struct YD_struct yd; /* Y database                      */
  YDE yde = &(yd.yde); /* Y element database                  */
  YDN ydn = &(yd.ydn); /* Y node database                     */
  YDX ydx = &(yd.ydx);
  YDXN ydxn = &(yd.ydxn);
  YDC ydc = &(yd.ydc); /* Y node database                     */

  INT i, j, l, k, j_tmp, index, n0, n1, n2, ii;
  DBL P0,P4, P8,t1,t2, CONT;
  DBL a1,a2,a3;
  FILE *fmshv = FILENULL;
  FILE *fmshs = FILENULL;
  //extern INT *i1l2q;
  extern INT nnode_l;

  printf ("in y3d_popullate \n\n");





  printf ("in y3d_populate 2 \n\n"); //-ao ---------

  /* get name of the problem */

  for (i = 0; i < ydx->ncmel; i++)
    {

      i1elto_v0[i] = ydx->i2cvto[i][0];
      i1elto_v1[i] = ydx->i2cvto[i][1];
      i1elto_v2[i] = ydx->i2cvto[i][2];
      i1elto_v3[i] = ydx->i2cvto[i][3];

    }

  for (i = 0; i < ydx->nsfel; i++)
    {
      i1elto_t_l0[i] = ydx->i2xbto[i][0];
      i1elto_t_l1[i] = ydx->i2xbto[i][1];
      i1elto_t_l2[i] = ydx->i2xbto[i][2];
    }

  for (j = 0; j < (ydxn->ncmno); j++)
    {
      d1nx[j] = ydxn->d2xncc[0][j];
      d1ny[j] = ydxn->d2xncc[1][j];
      d1nz[j] = ydxn->d2xncc[2][j];
//-PY changed it for 3D_fracture_coupling_with_multiphase
      /*
       d1nsvx_vol[j]=ydxn->d2xnvc[0][j];
       d1nsvy_vol[j]=ydxn->d2xnvc[1][j];
       d1nsvz_vol[j]=ydxn->d2xnvc[2][j];
       */

//			printf("coordinates velocity are %d %f %f %f %f %f %f \n\n",
//				   j,d1nx[j],d1ny[j],d1nz[j],d1nsvx[j],d1nsvy[j],d1nsvz[j]);
    }

  printf("nceml ncsfel ncmno %i %i %i \n\n", ydx->ncmel, ydx->nsfel,ydxn->ncmno);

  fmshv = fopen ("continuum.msh", "w");
  Ywgidcontinuum (yde, ydx, ydxn, fmshv); /* Write GiD continuum mesh file   */
  fclose (fmshv);

  fmshs = fopen ("surface.msh", "w");
  Ywgidsurface (yde, ydx, ydxn, fmshs); /* Write GiD surface mesh file   */
  fclose (fmshs);

 // Yring(ydx, ydn, ysp);   /*JXiang 3D coupling new  */

  for (i = 0; i < (3 * ydx->nsfel); i++)
    {
      //ring element
      i1elto0[i] = ydx->i2elto_r[i][0];
      i1elto1[i] = ydx->i2elto_r[i][1];
      i1elto2[i] = ydx->i2elto_r[i][2];
      i1elto3[i] = ydx->i2elto_r[i][3];

  //    printf("femdem ring element i1elto0[i]: %i %f %f %f \n",i, i1elto0[i],i1elto1[i],i1elto3[i]);
    }

//printf("femdem mesh output size i1elto0, i1elto0[i]: %d %d\n",i,i1elto0[i-1]);

  for (i = 0; i < 2 * (ydx->nnode); i++)
    { // ring coordinates
      d1nsvx[i] = ydx->d2ncc[0][i];
      d1nsvy[i] = ydx->d2ncc[1][i];
      d1nsvz[i] = ydx->d2ncc[2][i];
//printf("femdem ring coordinates: d1nsvx,i: %f %f %f %d \n",d1nsvx[i],d1nsvy[i],d1nsvz[i],i);
    }

  //printf("femdem mesh output size d1nsvx, d1nsvx[j], d1nsvx[j-1]: %d %f %f\n",i,d1nsvx[i], d1nsvx[i-1]);

  for (i = 0; i < 2 * ydx->nsfel; i++)
    {
      // ring surface
      i1elto_rt0[i] = ydx->i2elto_rt[i][0];
      i1elto_rt1[i] = ydx->i2elto_rt[i][1];
      i1elto_rt2[i] = ydx->i2elto_rt[i][2];
      //printf("femdem ring surface i1elto_rt0[i]: %f %f %f %d\n",i1elto_rt0[i],i1elto_rt1[i],i1elto_rt2[i],i);

    }


  for (i = 0; i < ydx->nsfel; i++)

   {


      n0 = ydx->i2xbto[i][0];
      n1 = ydx->i2xbto[i][1];
      n2 = ydx->i2xbto[i][2]; 
		
      P0=R0;    /*initial perm */
     // CONT=EPSILON*2.0; /* check aperture */
      CONT=1.0e-8; /* check aperture is above the user-defined minimum */

	/*a1=ydxn->d1shap[ydx->i1r2s[n0]];
	a2=ydxn->d1shap[ydx->i1r2s[n1]];
	a3=ydxn->d1shap[ydx->i1r2s[n2]]; */

//	a1=ydxn->d1shap[ydx->i1s2r[n0]];
//	a2=ydxn->d1shap[ydx->i1s2r[n1]];
//	a3=ydxn->d1shap[ydx->i1s2r[n2]];

	a1=ydxn->d1shap[n0];
	a2=ydxn->d1shap[n1];
	a3=ydxn->d1shap[n2];
      


	t1 = MAXIM(a1,a2); 
	t2 = MAXIM(t1, a3);

      if (t2>CONT)
    	{

//	V3DLen(t2, ydxn->d1shap[n0], ydxn->d1shap[n1], ydxn->d1shap[n2]); //-ao 1)>  aperture=length of 'displacement vector' 
	P0=(t2*t2*t2)/12;
	P4=P0;
	P8=P0;

	/*  P0 = (a1*a1*a1) / 12;   //-ao 3)> apeture=local displacement length
	  P4 = (a2*a2*a2) / 12;
	  P8 = (a3*a3*a3) / 12;	*/
	}
      else
	{

	  P0=R0;
	  P4=P0;
	  P8=P0;

	}


 	for(j=0;j<3;j++)
        {
       	      ii=i*3+j;

	      pp1[ii] = P0;
	      pp2[ii] = 0.0;
	      pp3[ii] = 0.0;
	      pp4[ii] = 0.0;
	      pp5[ii] = P4;
	      pp6[ii] = 0.0;
	      pp7[ii] = 0.0;
	      pp8[ii] = 0.0;
	      pp9[ii] = P8;
	}

    }

 


  printf ("out y3d_popullate \n\n");
}

void y3dfemdem_ (c1name, dt_f, d1npres, d1vx, d1vy, d1vz, d1vx_b, d1vy_b,
		 d1vz_b, d1dvx_b, d1dvy_b, d1dvz_b, vx_b, vy_b, vz_b, d1nmu, fx,
		 fy, fz, usl, vsl, wsl, axx, axy, ayy, axz, ayz, azz, d1flpore, d1nt_r, d1nt_r_lhs, d1nt_v)

  CHR *c1name;DBL *dt_f;DBL *d1npres;
  DBL *d1vx;DBL *d1vy;DBL *d1vz;
  DBL *d1vx_b;DBL *d1vy_b;DBL *d1vz_b;
  DBL *d1dvx_b;DBL *d1dvy_b;DBL *d1dvz_b;
  DBL *vx_b;DBL *vy_b;DBL *vz_b;
  DBL *d1nmu;DBL *fx;DBL *fy;DBL *fz;DBL *usl;DBL *vsl;DBL *wsl;
  DBL *axx;DBL *axy;DBL *ayy, *axz;DBL *ayz;DBL *azz;
  DBL *d1flpore;
  DBL *d1nt_r; DBL *d1nt_r_lhs; DBL *d1nt_v;
{

  printf ("in 3D FEMDEM \n\n");

  extern struct YD_struct yd; /* Y database                          */
  YDC ydc = &(yd.ydc); /* Y control database                  */
  YDE yde = &(yd.yde); /* Y element database                  */
  YDI ydi = &(yd.ydi); /* Y interaction database              */
  YDN ydn = &(yd.ydn); /* Y node database                     */
  YDO ydo = &(yd.ydo); /* Y output database                   */
  YDP ydp = &(yd.ydp); /* Y property database                 */
  YDB ydb = &(yd.ydb); /* Y boundary conditions database      */
  YDX ydx = &(yd.ydx);
  YDXN ydxn = &(yd.ydxn);
  YDK ydk = &(yd.ydk);
  YPAR ypar = &(yd.ypar);
  YSP ysp = &(yd.ysp); /*2 ydb is Y boundary condition database */


  extern INT nnode_l;
  extern INT flag;

// cj215
	PetscErrorCode ierr;
	Mat Chi;

  /* get name of the problem */
  INT istart, itmp, jtmp, l, k, m, n, i_bulk;
  INT nstep, i, j, j_tmp, index, jj;

  INT ii1, jj1, kk1;
  INT ii;

  INT p1[4] =
    { 4, 4, 5, 7 };
  INT p2[4] =
    { 6, 5, 6, 8 };
  INT p3[4] =
    { 7, 8, 9, 9 };

  DBL *d1denf;
  DBL *d1absorp;
  DBL *d1hystp;
  DBL *d1hy_CV;
  DBL press_rmp;



  DBL *f1dt;
  DBL *f1dt1;
  DBL *f1dt2;
  f1dt = TalDBL1 (ydn->mnopo);
  f1dt1 = TalDBL1 (ydn->mnopo);
  f1dt2 = TalDBL1 (ydn->mnopo);

  printf ("in 3D FEMDEM 1 \n\n");

  d1hy_CV = TalDBL1 (ydxn->ncmno+1);
  d1denf = TalDBL1 (ydn->mnopo);
  d1absorp = TalDBL1 (ydn->mnopo);
  d1hystp = TalzDBL1 (ydxn->ncmno+1);


    FILE *fmshv = FILENULL;
    FILE *fmshs = FILENULL;

  printf ("in 3D FEMDEM 2 \n");
  printf ("name: %s\n", c1name);

  istart = ydc->ncstep;

  for (j = 0; j < ydn->nnopo; j++)
    {
      ydn->d2nfv[0][j] = R0;
      ydn->d2nfv[1][j] = R0;
      ydn->d2nfv[2][j] = R0;
      d1denf[j] = R0;
      d1absorp[j] = R0;
      f1dt[j] = R1;
      f1dt1[j] = R1;
      f1dt2[j] = R1;
    }

  /* -ao changed this so that we fill an array with size shell_node_number with the pressure from Fluid model */
    FILE *fp1;
    fp1 = fopen ("pressure.txt", "w+");
  for (j = 0; j <ydx->nnode; j++)
    {
      jj=ydx->i1r2s[j];
      d1hystp[jj]=d1npres[j];
      fprintf (fp1, "node[%i] - pressure=%e ---jj= %i -- hystp =%e \n", j, d1npres[j], jj, d1hystp[jj]);
    }
    fprintf (fp1, " time step ------------ %f ------------", ydc->icouti);
   fclose (fp1);


  index = 0;
  for (j = 0; j < (ydxn->ncmno); j++)
    {
      for (k = 0; k < 100; k++)
	{
	  jtmp = ydxn->i2xbjno[j][k];
	  if ((jtmp >= 0) && (jtmp < ydn->nnopo))
	    {
              //-PY changed it for 3D_fracture_coupling_with_multiphase
	        ydn->d2nfv[0][jtmp]=d1vx_b[j];
	       ydn->d2nfv[1][jtmp]=d1vy_b[j];
	       ydn->d2nfv[2][jtmp]=d1vz_b[j];
	      d1denf[jtmp]=1000.0; //den_f[j];
	       d1absorp[jtmp]=1000.0/(*dt_f);
	       
	    }
	}
    }
  nstep = (int) (*dt_f / ydc->dcstec);

  printf ("number of loop in 3D  FEMDEM is %d\n\n", nstep);

  index = 0;
  //   Ywrs(c1name,&yd,0);   
  for (ydc->ncstep = istart; ydc->ncstep < istart + nstep; ydc->ncstep++)
    {

      index = 0;

      /*4 calculate nodal forces         */
      Yhpd (ydx, ydn,ydxn, d1hystp);
      //	 Ywrs(c1name,&yd,1);
      Yfd (yde, ydn, ydp, ydk, ydb, ydi, ydc, ydx);

      /*4 calculate contact detection    */
      Ycd (ydc, yde, ydi, ydn, ydp);

      /*4 detect contact interaction     */
	if(ydc->iflag[2]>0)
	{
		/*cj215 Prepare thermal contact interaction matrix Chi    */
		ierr = MatCreate(PETSC_COMM_WORLD,&Chi);CHKERRQ(ierr);
		ierr = MatSetSizes(Chi,PETSC_DECIDE,PETSC_DECIDE,ydxn->ncmno,ydxn->ncmno);CHKERRQ(ierr);
		ierr = MatSetFromOptions(Chi);CHKERRQ(ierr);
		ierr = MatSetUp(Chi);CHKERRQ(ierr);
		ierr = MatZeroEntries(Chi);CHKERRQ(ierr);
		ierr = MatAssemblyBegin(Chi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(Chi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}

	Yid (ydc, yde, ydi, ydn, ydp);

	if(ydc->iflag[2]>0)
	{
		/*cj215 Solve temperature    */
		ierr = MatAssemblyBegin(Chi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(Chi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

		Ythm(yde, ydn, ydp, ydk, ydc, ydb, ydx, Chi, d1nt_r, d1nt_r_lhs, d1nt_v, ydc->iflag[2], ydxn);
		ierr = MatDestroy(&Chi);CHKERRQ(ierr);
	}
			/*4  solve equations               */
      if (ydc->ncstep == istart)
	index = 1;
      if (ydc->ncstep == (istart + nstep - 1))
	index = 2;
      //-PY switch on the drag force
      if(ydc->iflag[0]>0)Yring_drag(ydn,ydx,d1nmu,fx,fy,fz,usl,vsl,wsl,axx,ayy,axy,azz,axz,ayz,d1vx,d1vy,d1vz);
      Ysd (ydc, yde, ydn, ydo, ydb, ydk, ydp, d1absorp, d1denf, index, f1dt,
	   f1dt1, f1dt2, ypar, ysp);
      //	 Ywrs(c1name,&yd,4);
      /*4  output results                */
      Yod (c1name, &yd);
      //	 Ywrs(c1name,&yd,ydc->ncstep);
      /*4  update time                   */
      ydc->dctime = ydc->dctime + ydc->dcstec;
    }

   for (j = 0; j < ydxn->ncmno; j++)  //-PY changed it
     {
       jtmp = ydxn->i2xbjno[j][0];
       vx_b[j] = ydn->d2nvc[0][jtmp];
       vy_b[j] = ydn->d2nvc[1][jtmp];
       vz_b[j] = ydn->d2nvc[2][jtmp];

       d1dvx_b[j] = vx_b[j] - d1vx_b[j];
       d1dvy_b[j] = vy_b[j] - d1vy_b[j];
       d1dvz_b[j] = vz_b[j] - d1vz_b[j];
}


  YSFmesh (ydc, ydx, ydxn, ydn, ydp, yde);

  ydx->nsfel_old = ydx->nsfel;
  fmshv = fopen ("continuum.msh", "w");
  Ywgidcontinuum (yde, ydx, ydxn, fmshv); /* Write GiD continuum mesh file   */
  fclose (fmshv);

  fmshs = fopen ("surface.msh", "w");
  Ywgidsurface (yde, ydx, ydxn, fmshs); /* Write GiD surface mesh file   */
  fclose (fmshs);
      if(ydc->iflag[1]>0)
      {
      Yap(yde, ydn, ydx);
      Yscale(yde, ydx, ydn, ydxn);      
      }
  Yring (ydx,ydxn, ydn, ysp); /*JXiang 3D coupling new  */

  for (j = 0; j < 3; j++)
    {
      f1dt[j] = f1dt[j] / nstep;
      f1dt1[j] = f1dt1[j] / nstep;
      f1dt2[j] = f1dt2[j] / nstep;
    }

  

  FREE(d1hystp);
  FREE(d1denf);
  FREE(d1absorp);
  FREE(d1hy_CV);
  FREE(f1dt);
  FREE(f1dt1);
  FREE(f1dt2);

  printf ("out 3D FEMDEM \n\n");

}

