/**********************************************************************/
/** Copyright (C) 2008,                                              **/
/** Queen Mary University of London (QMUL) & Imperial College        **/
/** of Science, Technology and Medicine (ICSTM). All rights reserved.**/
/** Implemented for you by Dr Jiansheng Xiang                        **
 
 
 
* This code is part of the Virtual Geoscience Workbench (VGW) developed
* jointly by ICSTM and QMUL through two related parallel projects at 
* ICSTM and QMUL respe
 ctively funded by EPSRC.
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
* Dr Jiansheng Xiang or Dr J.P.Latham 
* j.xiang@imperial.ac.uk or j.p.latham@imperial.ac.uk 
* ******************************************************************* */ 

/*1 File   Y3Dthm.c */
/* The Thermal explicit, implicit, contact and coupling is the work of Clement Joulin c.joulin15@imperial.ac.uk or cljoulin@gmail.com */
#include "Yproto.h"

/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/

/* Is element has face on boundary ? */
/* Only need to run once */
/* Very heavy algorithm comment if no thermal contact interaction*/
void Yface_boundary_heat
(INT **i2elto, INT **i2elto_t, INT nelem, INT nelem_t, INT *i1elbe, 
INT **i1elbe_h, DBL **d2ncc, INT *i1ptn, DBL *surf_t, INT nnopo)
{
	INT i,i_t,j,j_t;
	INT is_face;
	INT count;
	DBL det,surf_triangle;
	DBL mat[3][3];
	count = 0;
	DBL xt[3],yt[3],zt[3];
	DBL xn,yn,zn,s;
	s=0.0;
	DBL surf_tot=0.0;

	for (i=0;i<nelem;i++)
	{
		i1elbe_h[0][i]=-1;
		i1elbe_h[1][i]=-1;
//		Link node number to solid number 
		for(i_t=0;i_t<4;i_t++)
		{	
			i1ptn[i2elto[i][i_t]]=i2elto[i][4];
		}
	}
	for (i_t=0;i_t<nelem_t;i_t++)
	{
		// Identify boundary faces & nodes manually
		xt[0]= d2ncc[0][i2elto_t[i_t][0]];
		xt[1]= d2ncc[0][i2elto_t[i_t][1]];
		xt[2]= d2ncc[0][i2elto_t[i_t][2]];

		yt[0]= d2ncc[1][i2elto_t[i_t][0]];
		yt[1]= d2ncc[1][i2elto_t[i_t][1]];
		yt[2]= d2ncc[1][i2elto_t[i_t][2]];

		zt[0]= d2ncc[2][i2elto_t[i_t][0]];
		zt[1]= d2ncc[2][i2elto_t[i_t][1]];
		zt[2]= d2ncc[2][i2elto_t[i_t][2]];

		V3DCro(xn,yn,zn,xt[1]-xt[0],yt[1]-yt[0],zt[1]-zt[0],
			xt[2]-xt[0],yt[2]-yt[0],zt[2]-zt[0]);
		V3DNor(s,xn,yn,zn);		
		// Calculate the weight of surface for each boundary node for Pressure->Force cj215
		mat[0][0]=1.0;	      mat[0][1]=1.0;          mat[0][2]=1.0;
		mat[1][0]=xt[1]-xt[0];	mat[1][1]=yt[1]-yt[0];	mat[1][2]=zt[1]-zt[0];
		mat[2][0]=xt[2]-xt[0];	mat[2][1]=yt[2]-yt[0];	mat[2][2]=zt[2]-zt[0];

		YMATDET3(mat,det);
		surf_triangle=0.5*ABS(det); // Triangle surface formula from had2know.com to ddbl check
		surf_tot+=surf_triangle;
		surf_t[i_t]=surf_triangle;
		if(ABS(ABS(zn)-1.0)>EPSILON)
		{
			for (i=0;i<nelem;i++)
			{
				is_face=0;
				if(i1elbe[i]!=0)//if the i element is on a boundary
				{
					for(j_t=0;j_t<3;j_t++) //Loop all the node of the boundary element
					{
						for(j=0;j<4;j++) //Loop on the element nodes
						{
							if(i2elto_t[i_t][j_t]==i2elto[i][j])//Are nodes of the boundary element can be found in that element? 
							//"Where is that i2elto_t boundary face is to be found in the i2elto element list?"
							{
								is_face++;
							}
						}
					}
					if (is_face > 0 && is_face>i1elbe_h[0][i])
					{
						i1elbe_h[0][i]=is_face; //Number of nodes on the boundary for that element
						i1elbe_h[1][i]=i_t; //Link to the corresponding boundary element containing the most of that nodes.
						// This assumes each element can only have one face on the boundary
						if (is_face==3) count++;
					}
				}
			}
		}
	}
}

void Ythm_contact_init(yde, ydn , ydx)
	YDE yde; YDN ydn; YDX ydx;
	{
			/* Prepare boundary element faces identification for contact heat transfer */ 
		  yde->i1elbe_h = TalDBL2(2,yde->nelem);
		  ydn->i1ptn=TalINT1(ydn->nnopo);
		  ydx->surf_t=TalDBL1(ydx->nelem_t);

		  Yface_boundary_heat(yde->i2elto,
		  ydx->i2elto_t,yde->nelem,ydx->nelem_t,
		  yde->i1elbe,yde->i1elbe_h,ydn->d2nci, ydn->i1ptn, ydx->surf_t, ydn->nnopo);
	}

void Ythm_implicit_linear(  /* small strain elastic 4-noded tetrahedra  */
	INT nelem, INT nnopo, INT iprop, INT **i2elto, DBL dcstec,
	INT *i1elpr, INT *i1nopr, INT nelem_t, INT **i2xbtojn,
	DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,  
	DBL *d1ncix, DBL *d1nciy, DBL *d1nciz, 
	DBL *d1ntc, Mat A,
	DBL dpero,   DBL d1capa,  DBL d1tcon, 
	INT *i1bntp, DBL *d1bntp, DBL *d1bnhf, INT *i1bcvt, DBL *d1bcvt, DBL *d1bcvc, 
	INT flag_thm,DBL *d1nt_r, DBL *d1nt_r_lhs, DBL *d1nt_v,
	DBL *d1nvol, DBL dthick,  INT xnelno,  INT xnnode, INT *i1r2s)
{
	DBL nx,ny,nz,voli,volc;
	DBL F[3][3];  /* deformation gradient in global base delta ux/delta x */
	DBL F0[3][3]; /* initial local base */
	DBL FX[3][3]; /* current local base  also delta ux/delta X */
	DBL F0inv[3][3]; /* global base in initial local base */
	DBL FXinv[3][3]; /* global base in current local base */
	DBL Bij[3][NNODEX]; /*Temperature diffusion term + Capacity matrix */
	DBL d2dsh[3][4];
	DBL  detf;
	INT ielem;
	INT i,j,k,l,ni,nj;
	INT *i1elto;
	INT ipropp;
	DBL Btmp[4][4];
	DBL d2csmm[4][4];
	DBL d1emct;
	DBL delta_sq,k_e;
	FILE *fp;

	INT is_lumped = TRUE; // Is the mass maxtrix lumped ?
	INT is_lifting_scheme = FALSE; // Is the lifting scheme ON for boundary conditions? (not tested for more than one dirichlet BC)

	/*PETSC DECLARATION*/
	Vec            x, b, u;      /* approx solution, RHS, exact solution */
	KSP            ksp;         /* linear solver context */
	PC             pc;           /* preconditioner context */
	PetscErrorCode ierr;
	PetscInt       rows[nnopo],jstore;
	PetscMPIInt    size;
	PetscBool      nonzeroguess = PETSC_FALSE;
	PetscViewer    view_out;
	PetscScalar MT_KTi, Mij,value_temp,Theta=0.5; 

	PetscScalar chi[1];
//	PetscOptions options;
	char opt_string[]="-ksp_gmres_restart 30";// "-ksp_type fgmres\ -snes_grid_sequence 2 \ ...

	/////////////////////// 

	/*PETSC INITIALISATION*/
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
	ierr = PetscOptionsGetInt(NULL,"-n",&nnopo,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,"-nonzero_guess",&nonzeroguess,NULL);CHKERRQ(ierr);

//	/* Compute the matrix and right-hand-side vector that define
//	the linear system, Ax = b.*/
	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
	ierr = VecSetSizes(x,PETSC_DECIDE,nnopo);CHKERRQ(ierr);
	ierr = VecSetFromOptions(x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&u);CHKERRQ(ierr);

	ierr = VecZeroEntries(b);CHKERRQ(ierr);
	ierr = VecZeroEntries(u);CHKERRQ(ierr);

	///////////////////////
	// Apply contact heat transfer interactions
	for (ni=0;ni<nnopo;ni++){ierr = VecSetValue(u,ni,-d1ntc[ni]*(1.0-Theta),ADD_VALUES);CHKERRQ(ierr);}
	MatMult(A,u,b);
	if(Theta>0){ierr = MatScale(A,Theta);}else{ierr = MatZeroEntries(A);}CHKERRQ(ierr);
	MatZeroEntries(A);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	///////////////////////

	d2csmm[0][0]=0.1;d2csmm[0][1]=0.05;d2csmm[0][2]=0.05;d2csmm[0][3]=0.05;
  	d2csmm[1][0]=0.05;d2csmm[1][1]=0.1;d2csmm[1][2]=0.05;d2csmm[1][3]=0.05;
  	d2csmm[2][0]=0.05;d2csmm[2][1]=0.05;d2csmm[2][2]=0.1;d2csmm[2][3]=0.05;
  	d2csmm[3][0]=0.05;d2csmm[3][1]=0.05;d2csmm[3][2]=0.05;d2csmm[3][3]=0.1;


	//Jacobian calculs corresponds to this shape function only do not change
	d2dsh[0][0]= -1.0;
	d2dsh[0][1]=  1.0;
	d2dsh[0][2]=  0.0;
	d2dsh[0][3]=  0.0;

	d2dsh[1][0]= -1.0;
	d2dsh[1][1]=  0.0;
	d2dsh[1][2]=  1.0;
	d2dsh[1][3]=  0.0;

	d2dsh[2][0]= -1.0;
	d2dsh[2][1]=  0.0;
	d2dsh[2][2]=  0.0;
	d2dsh[2][3]=  1.0;

	 for(ielem=0;ielem<nelem;ielem++)
	  { 
		  if(i1elpr[ielem]==iprop)
	    { 
		    i1elto = i2elto[ielem];
			for(i=1;i<4;i++)
	      {
			F0[0][i-1]=d1ncix[i1elto[i]]-d1ncix[i1elto[0]]; /* init. base */
			F0[1][i-1]=d1nciy[i1elto[i]]-d1nciy[i1elto[0]];
			F0[2][i-1]=d1nciz[i1elto[i]]-d1nciz[i1elto[0]];
			FX[0][i-1]=d1nccx[i1elto[i]]-d1nccx[i1elto[0]]; /* curr. base */
			FX[1][i-1]=d1nccy[i1elto[i]]-d1nccy[i1elto[0]];
			FX[2][i-1]=d1nccz[i1elto[i]]-d1nccz[i1elto[0]];
	      }  
	      YMATINV3(F0,F0inv,voli);         /* global base in initial local coordinates    */
	      YMATINV3(FX,FXinv,volc);         /* global base in current local coordinates    */

            for(i=1;i<4;i++) //To calculate the Jacobian matrix, not the transpose
            {
                FX[i-1][0]=d1nccx[i1elto[i]]-d1nccx[i1elto[0]]; /* curr. base */
                FX[i-1][1]=d1nccy[i1elto[i]]-d1nccy[i1elto[0]];
                FX[i-1][2]=d1nccz[i1elto[i]]-d1nccz[i1elto[0]];
            }
            YMATINV3(FX,FXinv,volc);         /* global base in current local coordinates    */

		  	for (k = 0; k < 3; k++)      /* B matrix gradient of the shape function */
			{
				for(i=0;i<NNODEX;i++)
				{
					Bij[k][i]=(FXinv[k][0]*d2dsh[0][i]+FXinv[k][1]*d2dsh[1][i]+FXinv[k][2]*d2dsh[2][i]);
				}
			}

			for(i=0;i<NNODEX;i++)
	            {
	                for (j = 0; j < NNODEX; j++)
	                {
	                    Btmp[i][j]=Bij[0][i]*Bij[0][j]+Bij[1][i]*Bij[1][j]+Bij[2][i]*Bij[2][j];
	                }
	            }

            d1emct=dpero*d1capa*voli/6.0;

      	for(i=0;i<NNODEX;i++)      /* Nodal temperatures */
        	{
			ni = i1elto[i];
			MT_KTi  = 0.0;
			for (j = 0; j < NNODEX; j++)
			{
				nj = i1elto[j];
				Mij=0.0;

				//Diffusion
				MT_KTi -= (1.0-Theta)*d1tcon*Btmp[i][j]*volc/6.0*d1ntc[nj];
				Mij += Theta*d1tcon*Btmp[i][j]*volc/6.0;

				// Mass matrix
				if(is_lumped == TRUE) 
				{
					if(i==j) Mij += RP25*dpero*d1capa*voli/6.0/dcstec; //Lumped
				}else
				{
					MT_KTi += d2csmm[i][j]*d1emct*d1ntc[nj]/dcstec; //Consistent
					Mij += d2csmm[i][j]*d1emct/dcstec; //Consistent
				}

				//Thermal coupling
				if(flag_thm==1) if(i==j) Mij += RP25*dpero*d1capa/dcstec*voli/6.0;//*(1.0-d1nt_r_lhs[ni]); //Volume temperature relaxation (method 1) removed RP25*voli/6.0
				if(flag_thm==5) if(i==j) Mij += RP25*d1tcon*voli/6.0;

				//Assign LHS
				ierr = MatSetValues(A,1,&ni,1,&nj,&Mij,ADD_VALUES);CHKERRQ(ierr);
			}
			//Mass matrix
			if(is_lumped == TRUE) MT_KTi += RP25*dpero*d1capa*voli/6.0*d1ntc[ni]/dcstec;//Lumped

			//Thermal coupling
			if(flag_thm==1) MT_KTi +=RP25*dpero*d1capa/dcstec*voli/6.0*d1nt_v[ni]; //Volume temperature relaxation (method 1)
			if(flag_thm==5) MT_KTi +=RP25*d1tcon*voli/6.0*d1nt_v[ni];

			//Assign RHS
			ierr = VecSetValue(b,ni,MT_KTi,ADD_VALUES);CHKERRQ(ierr);
      	}
	    }
	 }


	ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
//	VecView(b,PETSC_VIEWER_STDOUT_SELF);

	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//	MatView(A,PETSC_VIEWER_STDOUT_SELF);

	//BOUNDARY CONDITIONS
	i=0;
	//Find BC and lift

	for(ni=0;ni<nnopo;ni++)
	{
		ipropp = i1nopr[ni];
		if (i1bntp[ipropp]==1)
		{
			rows[i] = ni;
			i++;
			MT_KTi = d1bntp[ipropp];
			if (is_lifting_scheme == TRUE)
			{		
				jstore = ni;
				ierr = MatGetColumnVector(A,u,jstore);CHKERRQ(ierr);
				VecAXPY(b,-MT_KTi,u);
			}
		}
	}


	for(ni=0;ni<nnopo;ni++)
	{
		ipropp = i1nopr[ni];

		// DIRICHLET
		if (i1bntp[ipropp]==1)
		{
			MT_KTi = d1bntp[ipropp];
			ierr = VecSetValue(b,ni,MT_KTi,INSERT_VALUES);CHKERRQ(ierr);
		}

		//CONVECTION BOUNDARY CONDITION (NEUMANN)
		if (i1bcvt[ipropp]==1)
		{
			MT_KTi = d1bcvc[ipropp]*d1bcvt[ipropp];
			ierr = VecSetValue(b,ni,MT_KTi,INSERT_VALUES);CHKERRQ(ierr);

			Mij = d1bcvc[ipropp];
			ierr = MatSetValues(A,1,&ni,1,&ni,&Mij,ADD_VALUES);CHKERRQ(ierr);
		}
	}

	// Zero rows & columns
	if(i>0) 
	{
		if (is_lifting_scheme == FALSE) MatZeroRows(A,i,rows,1.0,0,0);//no lifting-scheme
		if (is_lifting_scheme == TRUE)  MatZeroRowsColumns(A,i,rows,1.0,0,0);//lifting-scheme
	}



//	APPLY FLUID-SOLID COUPLING BC OVER RING (method 2 & 3)
	if(flag_thm==2 || flag_thm ==3) // || flag_thm ==4)
	{
		delta_sq=dthick*dthick;
		for (i=0; i<xnnode; i++)
		{
			ni=i1r2s[i];

			MT_KTi = d1tcon/delta_sq*d1nvol[ni]*d1nt_r[i+xnnode];
//			MT_KTi = RP25*d1tcon/delta_sq*(d1nvol[i]+d1nvol[i+xnnode])*d1nt_r[i+xnnode];
			ierr = VecSetValue(b,ni,MT_KTi,ADD_VALUES);CHKERRQ(ierr);

			Mij = d1tcon/delta_sq*d1nvol[ni];
//			Mij = RP25*d1tcon/delta_sq*(d1nvol[i]+d1nvol[i+xnnode]);
			ierr = MatSetValues(A,1,&ni,1,&ni,&Mij,ADD_VALUES);CHKERRQ(ierr);
		}
	}

	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//	MatView(A,PETSC_VIEWER_STDOUT_SELF);

	ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
//	VecView(b,PETSC_VIEWER_STDOUT_SELF);

	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = PCSetType(pc,PCGAMG);CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp,1.e-100,1.e-8,PETSC_DEFAULT,200);CHKERRQ(ierr);
	ierr = KSPSetType(ksp,KSPGMRES);CHKERRQ(ierr);
	ierr = PetscOptionsInsertString(opt_string);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

	if (nonzeroguess)
	{
		PetscScalar p = .5;
		ierr = VecSet(x,p);CHKERRQ(ierr);
		ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
	}
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
//	ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//	VecView(x,PETSC_VIEWER_STDOUT_SELF);

	for (i=0; i<nnopo; i++)
	{
		ierr = VecGetValues(x,1,&i, &value_temp);CHKERRQ(ierr);
		d1ntc[i]= value_temp;
		d1nt_v[i] = value_temp;
	}

	if(flag_thm==2 || flag_thm ==3)
	{
		for (i=0; i<2*xnnode; i++)
		{
			ni = i1r2s[i];
			d1nt_r[i]=d1tcon/delta_sq*d1ntc[ni];
			d1nt_r_lhs[i] = d1tcon/delta_sq;
		}
	}

	ierr = VecDestroy(&x);CHKERRQ(ierr); 
	ierr = VecDestroy(&b);CHKERRQ(ierr);
	ierr = VecDestroy(&u);CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
}

void Ythm_explicit_linear(  /* small strain elastic 4-noded tetrahedra  */
	INT nelem, INT nnopo, INT iprop, INT **i2elto, DBL dcstec,
	INT *i1elpr, INT *i1nopr, INT nelem_t, INT **i2xbtojn,
	DBL *d1nccx, DBL *d1nccy, DBL *d1nccz,  
	DBL *d1ncix, DBL *d1nciy, DBL *d1nciz, 
	DBL *d1ntc, Mat A,
	DBL dpero,   DBL d1capa,  DBL d1tcon, 
	INT *i1bntp, DBL *d1bntp, DBL *d1bnhf, INT *i1bcvt, DBL *d1bcvt, DBL *d1bcvc, 
	INT flag_thm,DBL *d1nt_r, DBL *d1nt_r_lhs, DBL *d1nt_v,
	DBL *d1nvol, DBL dthick,  INT xnelno,  INT xnnode, INT *i1r2s)
{
	DBL nx,ny,nz,voli,volc;
	DBL F[3][3];  /* deformation gradient in global base delta ux/delta x */
	DBL F0[3][3]; /* initial local base */
	DBL FX[3][3]; /* current local base  also delta ux/delta X */
	DBL F0inv[3][3]; /* global base in initial local base */
	DBL FXinv[3][3]; /* global base in current local base */
	DBL Bij[3][NNODEX];
	DBL Btmp[4][4];
	DBL d2csmm[4][4];
	DBL d2dsh[3][4];
	DBL *Lii, *KTi;
//	DBL Lii[nnopo], KTi[nnopo]; 
	DBL  detf;
	INT ielem;
	INT i,j,k,l,ni,nj;
	INT *i1elto;
	INT ipropp;

	DBL delta_sq,k_e;

	d2csmm[0][0]=0.1;d2csmm[0][1]=0.05;d2csmm[0][2]=0.05;d2csmm[0][3]=0.05;
	d2csmm[1][0]=0.05;d2csmm[1][1]=0.1;d2csmm[1][2]=0.05;d2csmm[1][3]=0.05;
	d2csmm[2][0]=0.05;d2csmm[2][1]=0.05;d2csmm[2][2]=0.1;d2csmm[2][3]=0.05;
	d2csmm[3][0]=0.05;d2csmm[3][1]=0.05;d2csmm[3][2]=0.05;d2csmm[3][3]=0.1;


	//Jacobian calculs corresponds to this shape function only do not change
	d2dsh[0][0]= -1.0;
	d2dsh[0][1]=  1.0;
	d2dsh[0][2]=  0.0;
	d2dsh[0][3]=  0.0;

	d2dsh[1][0]= -1.0;
	d2dsh[1][1]=  0.0;
	d2dsh[1][2]=  1.0;
	d2dsh[1][3]=  0.0;

	d2dsh[2][0]= -1.0;
	d2dsh[2][1]=  0.0;
	d2dsh[2][2]=  0.0;
	d2dsh[2][3]=  1.0;

	KTi = TalDBL1(nnopo);
	Lii = TalDBL1(nnopo);

	for (ni = 0; ni < nnopo; ni++) 
	{
		Lii[ni] = 0.0;
		KTi[ni]=0.0;
	}

	for(ielem=0;ielem<nelem;ielem++)
	{ 
	//	if(i1elpr[ielem]==iprop) //When dealing with different discrete blocks, needs to be deactivated if you want to solve everything in the same matrix
	//	{ 
		i1elto = i2elto[ielem];

		for(i=1;i<4;i++)
		{ 
			F0[0][i-1]=d1ncix[i1elto[i]]-d1ncix[i1elto[0]]; /* init. base */
			F0[1][i-1]=d1nciy[i1elto[i]]-d1nciy[i1elto[0]];
			F0[2][i-1]=d1nciz[i1elto[i]]-d1nciz[i1elto[0]];
			FX[0][i-1]=d1nccx[i1elto[i]]-d1nccx[i1elto[0]]; /* curr. base */
			FX[1][i-1]=d1nccy[i1elto[i]]-d1nccy[i1elto[0]];
			FX[2][i-1]=d1nccz[i1elto[i]]-d1nccz[i1elto[0]];
		}  
		YMATINV3(F0,F0inv,voli);         /* global base in initial local coordinates    */
		YMATINV3(FX,FXinv,volc);         /* global base in current local coordinates    */

		for(i=1;i<4;i++) //To calculate the Jacobian matrix, not the transpose
		{
			FX[i-1][0]=d1nccx[i1elto[i]]-d1nccx[i1elto[0]]; /* curr. base */
			FX[i-1][1]=d1nccy[i1elto[i]]-d1nccy[i1elto[0]];
			FX[i-1][2]=d1nccz[i1elto[i]]-d1nccz[i1elto[0]];
		}
		YMATINV3(FX,FXinv,volc);         /* global base in current local coordinates    */

		for (k = 0; k < 3; k++)      /* B matrix gradient of the shape function */
		{
			for(i=0;i<NNODEX;i++)
			{
				Bij[k][i]=(FXinv[k][0]*d2dsh[0][i]+FXinv[k][1]*d2dsh[1][i]+FXinv[k][2]*d2dsh[2][i]);
			}
		}

		for(i=0;i<NNODEX;i++)
		{
			for (j = 0; j < NNODEX; j++)
			{
			    Btmp[i][j]=Bij[0][i]*Bij[0][j]+Bij[1][i]*Bij[1][j]+Bij[2][i]*Bij[2][j];
			}
		}        	        
		for(i=0;i<NNODEX;i++)      /* Nodal temperatures */
		{
			ni = i1elto[i];
			for (j = 0; j < NNODEX; j++)
			{
				nj = i1elto[j];
				KTi[ni] += d1tcon*Btmp[i][j]*volc/6.0*d1ntc[nj];
			}
			Lii[ni] += RP25*dpero*d1capa*voli/6.0;
		}
	}

	for (ni = 0; ni < nnopo; ni++) 
	{
		ipropp = i1nopr[ni];
		if (i1bntp[ipropp]==1)
		{
			d1ntc[ni]=d1ntc[ni] - dcstec*(KTi[ni] + (d1ntc[ni] - d1bntp[ipropp]))/Lii[ni];
		}else
		{
				d1ntc[ni]= d1ntc[ni] - dcstec*KTi[ni]/Lii[ni];
		}
	}

	free(KTi);
	free(Lii);
}
void Ythm(yde, ydn, ydp, ydk, ydc, ydb, ydx, Chi, d1nt_r, d1nt_r_lhs, d1nt_v, flag_thm, ydxn)
YDE yde; YDN ydn; YDP ydp; YDK ydk; YDC ydc; YDB ydb; YDX ydx; Mat Chi; DBL *d1nt_r; DBL *d1nt_r_lhs; DBL *d1nt_v; INT flag_thm; YDXN ydxn;
{

	INT iprop,i,j,k;
	PetscErrorCode ierr;
	DBL *d1xnvol;

	//Continuous ring volume
	d1xnvol = DBL1NULL;
	d1xnvol = TalDBL1 (ydn->mnopo);

	for (i = 0; i < ydn->mnopo; i++) 
	{
		d1xnvol[ydn->i1ntoC2D[i]] = 0.0;
	}
	for (i = 0; i < ydn->mnopo; i++) 
	{
		d1xnvol[ydn->i1ntoC2D[i]] += ydn->d1nvol[i];
	}

	//		  Ythm_contact_init(yde, ydn , ydx); // Not functional for discontinuous mesh, to be reviewed for the current version

	for(iprop=0; iprop<ydp->nprop; iprop++)
	{
		if((ydp->i1ptyp[iprop])==(YTE3TET4ELS))
		{
			Ythm_implicit_linear(
				ydx->ncmel, ydxn->ncmno, iprop, ydx->i2cvto, ydc->dcstec,
				yde->i1elpr, ydn->i1nopr, ydx->nelem_t, ydx->i2xbtojn,
				ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2ncc[2],
				ydn->d2nci[0], ydn->d2nci[1], ydn->d2nci[2],
				ydxn->d1xntc, Chi, 
				ydp->d1pero[iprop], ydp->d1capa[0], ydp->d1tcon[0],
				ydb->i1bntp, ydb->d1bntp, ydb->d1bnhf, ydb->i1bcvt, ydb->d1bcvt, ydb->d1bcvc,
				flag_thm, d1nt_r, d1nt_r_lhs, d1nt_v,
				d1xnvol, ydx->dthick, ydx->nelno, ydx->nnode, ydx->i1r2s);

//			Ythm_explicit_linear(
//				ydx->ncmel, ydxn->ncmno, iprop, ydx->i2cvto, ydc->dcstec,
//				yde->i1elpr, ydn->i1nopr, ydx->nelem_t, ydx->i2xbtojn,
//				ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2ncc[2],
//				ydn->d2nci[0], ydn->d2nci[1], ydn->d2nci[2],
//				ydxn->d1xntc, Chi, 
//				ydp->d1pero[iprop], ydp->d1capa[0], ydp->d1tcon[0],
//				ydb->i1bntp, ydb->d1bntp, ydb->d1bnhf, ydb->i1bcvt, ydb->d1bcvt, ydb->d1bcvc,
//				flag_thm, d1nt_r, d1nt_r_lhs, d1nt_v,
//				ydx->d1xnvol, ydx->dthick, ydx->nelno, ydx->nnode, ydx->i1r2s);

		}
	}

//  From continuous to discontinuous temperature
	for(j=0;j<100;j++)
	{
		for(i=0;i<ydxn->ncmno;i++)
		{
			k=ydxn->i2xbjno[i][j];
			/*discontinuous*/ydn->d1ntc[k] += ydxn->d1xntc[i];/*continuous*/
		}
	}

}



/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/
