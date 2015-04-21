
!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"

module solvers_module

  use fldebug
  use fields
  use Petsc_tools
  use sparse_tools_petsc
  use solvers
  use global_parameters, only: OPTION_PATH_LEN, dumping_in_sat
  use spud

  use state_module
  use halo_data_types
#ifdef HAVE_PETSC_MODULES
  use petsc 
#if PETSC_VERSION_MINOR==0
  use petscvec 
  use petscmat 
  use petscksp 
  use petscpc
#endif
#endif
  use Copy_Outof_State
  use shape_functions_Linear_Quadratic
  use shape_functions_prototype

  implicit none

#include "petscversion.h"
#ifdef HAVE_PETSC_MODULES
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscpcdef.h"
#else
#include "finclude/petsc.h"
#if PETSC_VERSION_MINOR==0
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#endif
#endif

  private

  public :: solver, PRES_DG_MULTIGRID, CMC_Agglomerator_solver, Fix_to_bad_elements,&
         BoundedSolutionCorrections, Trust_region_correction, Set_Saturation_between_bounds, Set_Saturation_to_sum_one

  interface solver
     module procedure solve_via_copy_to_petsc_csr_matrix
  end interface
  
contains

! -----------------------------------------------------------------------------

  subroutine solve_via_copy_to_petsc_csr_matrix( A, &
       x, b, findfe, colfe, option_path, block_size )

    !!< Solve a matrix Ax = b system via copying over to the
    !!< petsc_csr_matrix type and calling the femtools solver
    !!< using the spud options given by the field options path

    integer, dimension( : ), intent(in) :: findfe
    integer, dimension( : ), intent(in) :: colfe
    real, dimension( : ), intent(in) :: a, b
    real, dimension( : ), intent(inout) :: x
    character( len=* ), intent(in) :: option_path
    !This optional argument is to create a block-diagonal matrix
    !and therefore to use block-solvers
    integer, optional, intent(in) :: block_size

    ! local variables
    integer :: i, j, k, rows
    integer, dimension( : ), allocatable :: dnnz
    type(petsc_csr_matrix) :: matrix
    integer :: size_of_block

    size_of_block = 1
    if (present(block_size)) size_of_block = block_size

    rows = size( x )
    assert( size( x ) == size( b ) )
    assert( size( a ) == size( colfe ) )
    assert( size( x ) + 1 == size( findfe ) )
    ewrite(3,*) rows+1, size(findfe)

    allocate( dnnz( rows/size_of_block ) ) ; dnnz = 0
    ! find the number of non zeros per row
    do i = 1, size( dnnz )
        dnnz( i ) =(findfe( i+1 ) - findfe( i ))
    end do
    call allocate( matrix, rows/size_of_block, rows/size_of_block, dnnz, dnnz,&
    (/1, 1/), name = 'dummy', element_size=size_of_block)

    call zero( matrix )
    ! add in the entries to petsc matrix
    do i = 1, rows
        do j = findfe( i ), findfe( i+1 ) - 1
            k = colfe( j )
            call addto( matrix, blocki = 1, blockj = 1, i = i, j = k, val = a( j ) )
        end do
    end do

    call assemble( matrix )

    call petsc_solve_scalar_petsc_csr_mp( x, matrix, b, rows, trim( option_path ) )

    ! deallocate as needed
    deallocate( dnnz )
    call deallocate( matrix )

    return
  end subroutine solve_via_copy_to_petsc_csr_matrix


  subroutine petsc_solve_scalar_petsc_csr_mp( x, matrix, rhs, rows, option_path )

    real, dimension( : ), intent(inout) :: x
    type( petsc_csr_matrix ), intent(inout) :: matrix
    real, dimension( : ), intent(in) :: rhs
    integer, intent(in) :: rows
    character( len=* ), intent(in) :: option_path

    KSP :: ksp
    Vec :: y, b

    character(len=OPTION_PATH_LEN) :: solver_option_path
    integer :: ierr

    assert( size( x ) == size( rhs ) )
    assert( size( x ) == size( matrix, 2 ) )
    assert( size( rhs ) == size( matrix, 1 ) )

    solver_option_path = complete_solver_option_path( option_path )

    call SetupKSP( ksp, matrix%M, matrix%M, solver_option_path, .false., &
         matrix%column_numbering, .true. )

    b = PetscNumberingCreateVec( matrix%column_numbering )
    call VecDuplicate( b, y, ierr )


    ! copy array into PETSc vecs
    call VecSetValues( y, rows, &
         matrix%row_numbering%gnn2unn( 1:rows, 1 ), &
         x, INSERT_VALUES, ierr )
    call VecAssemblyBegin( y, ierr )
    call VecAssemblyEnd( y, ierr )

    call VecSetValues( b, rows, &
         matrix%row_numbering%gnn2unn( 1:rows, 1 ), &
         rhs, INSERT_VALUES, ierr )
    call VecAssemblyBegin( b, ierr )
    call VecAssemblyEnd( b, ierr )

    call KSPSolve( ksp, b, y, ierr )

    ! copy back the result
    call VecGetValues( y, rows, &
         matrix%row_numbering%gnn2unn( 1:rows, 1 ), &
         x, ierr )

    ! destroy all PETSc objects and the petsc_numbering
    call petsc_solve_destroy_petsc_csr( y, b, ksp, solver_option_path )

    

    return
  end subroutine petsc_solve_scalar_petsc_csr_mp

    !Clone of the same subroutine in femtools/Solvers.F90
  subroutine petsc_solve_destroy_petsc_csr( y, b, ksp, solver_option_path )

      type(Vec), intent(inout):: y
      type(Vec), intent(inout):: b
      type(KSP), intent(inout):: ksp
      character(len=*), intent(in):: solver_option_path

      type(PC) :: pc
      integer ierr

      call VecDestroy(y, ierr)
      call VecDestroy(b, ierr)
      call KSPGetPC(ksp, pc, ierr)
      call KSPDestroy(ksp, ierr)

  end subroutine petsc_solve_destroy_petsc_csr


  !DO NOT REMOVE, CURRENTLY UNUSED BUT WE PLAN TO CONTINUE!!
  SUBROUTINE PRES_DG_MULTIGRID(CMC, CMC_PRECON, IGOT_CMC_PRECON, P, RHS, &
       NCOLCMC, CV_NONODS, FINDCMC, COLCMC, MIDCMC, &
       totele, cv_nloc, x_nonods, cv_ndgln, x_ndgln )
    !
    ! Solve CMC * P = RHS for RHS.
    ! form a discontinuous pressure mesh for pressure...
    implicit none
    INTEGER, intent( in ) ::  NCOLCMC, CV_NONODS, totele, cv_nloc, x_nonods, IGOT_CMC_PRECON
! IGOT_CMC_PRECON=1 or 0 (1 if we have a preconditioning matrix)
    REAL, DIMENSION( : ), intent( in ) ::  CMC
    REAL, DIMENSION( :), intent( in ) ::  CMC_PRECON
    REAL, DIMENSION( : ), intent( inout ) ::  P
    REAL, DIMENSION( : ), intent( in ) :: RHS
    INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
    INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
    INTEGER, DIMENSION( : ), intent( in ) :: MIDCMC
    INTEGER, DIMENSION( : ), intent( in ) :: cv_ndgln, x_ndgln

    REAL ERROR, RELAX, RELAX_DIAABS, RELAX_DIA
    INTEGER N_LIN_ITS, NGL_ITS
    LOGICAL ONE_PRES_SOLVE, ONE_SOLVE_1IT

    PARAMETER(ERROR=1.E-15, RELAX=0.1, RELAX_DIAABS=0.0)
    PARAMETER(RELAX_DIA=2.0, N_LIN_ITS=2)
    PARAMETER(NGL_ITS=50) ! maybe we need to increase this...
   ! PARAMETER(NGL_ITS=50) ! maybe we need to increase this...
   ! PARAMETER(NGL_ITS=50) ! maybe we need to increase this...
    PARAMETER(ONE_PRES_SOLVE=.false., ONE_SOLVE_1IT=.false.) ! Only solve for one cty pressure - distribut to DG pressure nodes. 
!    PARAMETER(ONE_PRES_SOLVE=.true., ONE_SOLVE_1IT=.false.) ! Only solve for one cty pressure - distribut to DG pressure nodes. 

    !  PARAMETER(ERROR=1.E-15, RELAX=0.05, RELAX_DIAABS=0.0)
    !  PARAMETER(RELAX_DIA=2.0, N_LIN_ITS=2)
    !  PARAMETER(NGL_ITS=2000)

    ! RELAX: overall relaxation coeff; =1 for no relaxation. 
    ! RELAX_DIAABS: relaxation of the absolute values of the sum of the row of the matrix;
    !                      - recommend >=2 for hard problems, =0 for easy
    ! RELAX_DIA: relaxation of diagonal; =1 no relaxation (normally applied). 
    ! N_LIN_ITS = no of linear iterations
    ! ERROR= solver tolerence between 2 consecutive iterations
    ! NGL_ITS = no of global its

    integer, dimension( : ), allocatable :: findcmc_small, colcmc_small, midcmc_small, &
         MAP_DG2CTY
    real, dimension( : ), allocatable :: cmc_small, resid_dg, resid_dg2, resid_cty, &
         nods_sourou, DP_DG, DP_SMALL, USTEP
    integer :: ele,cv_iloc, dg_nod, cty_nod, jcolcmc, jcolcmc_small
    integer :: mx_ncmc_small, ncmc_small, count, count2, count3, GL_ITS, col
    real :: OPT_STEP

    character(len=OPTION_PATH_LEN) :: path = "/tmp/pressure"
    integer :: stat

    ! obtain sparcity of a new matrix 
    mx_ncmc_small = ncolcmc * 4
    allocate( FINDcmc_small(x_nonods+1) )
    allocate( colcmc_small(mx_ncmc_small) )
    allocate( midcmc_small(x_nonods) )
    allocate( MAP_DG2CTY(cv_nonods) )

    EWRITE(3,*)'BEFORE pousinmc'
    ! lump the pressure nodes to take away the discontinuity...
    DO ELE = 1, TOTELE
       DO CV_ILOC = 1, CV_NLOC
          dg_nod = (ele-1) * cv_nloc + cv_iloc
          cty_nod = x_ndgln( (ele-1) * cv_nloc + cv_iloc)
          MAP_DG2CTY(dg_nod) = cty_nod
       END DO
    END DO

    CALL GET_SPAR_CMC_SMALL(FINDCMC_SMALL, COLCMC_SMALL, MIDCMC_SMALL, &
         MX_NCMC_SMALL, NCMC_SMALL, CV_NONODS, X_NONODS, MAP_DG2CTY, &
         FINDCMC, COLCMC, NCOLCMC)

    ewrite(3,*)'FINDCMC_SMALL:', FINDCMC_SMALL
    ewrite(3,*)'COLCMC_SMALL(1:NCMC_SMALL):', COLCMC_SMALL(1:NCMC_SMALL)
    !stop 3292

    allocate( cmc_small(ncmc_small) )
    allocate( resid_dg(cv_nonods) )
    allocate( resid_cty(x_nonods) )
    allocate( dp_small(x_nonods) )
    allocate( dp_dg(cv_nonods) ) 
    allocate( USTEP(cv_nonods) )
    allocate( nods_sourou(x_nonods) )

    ewrite(3,*)'***forming cmc_small:'
    CMC_SMALL = 0.
    DO dg_nod = 1, CV_NONODS
       cty_nod = MAP_DG2CTY(dg_nod)
       ! add row dg_nod to row cty_nod of cty mesh
       DO COUNT = FINDCMC(DG_NOD), FINDCMC(DG_NOD+1) - 1
          jcolcmc = COLCMC(COUNT)
          jcolcmc_small = MAP_DG2CTY(jcolcmc)
          count2 = 0
          DO COUNT3 = FINDCMC_small(cty_NOD), FINDCMC_small(cty_NOD+1) - 1
             !ewrite(3,*)'dg_nod,cty_nod,jcolcmc_small,colcmc_small(count3):', &
             !     dg_nod,cty_nod,jcolcmc_small,colcmc_small(count3)
             if(colcmc_small(count3)==jcolcmc_small) count2=count3
          end do
          if(count2==0) then
             ewrite(3,*)'could not find coln'
             stop 3282
          end if
         ! IF((IGOT_CMC_PRECON==0).or.(ONE_PRES_SOLVE)) THEN
          IF(.true.) THEN
         ! IF(.false.) THEN
             CMC_SMALL(COUNT2) = CMC_SMALL(COUNT2) + CMC(COUNT)  
          ELSE 
          !  stop 6227
             CMC_SMALL(COUNT2) = CMC_SMALL(COUNT2) + CMC_PRECON(COUNT)  
          ENDIF            
       END DO
    END DO

! Only solve for one pressure...
    IF(ONE_PRES_SOLVE) THEN
       ! Map resid_dg to resid_cty as well as the solution:
       resid_cty=0.
       do dg_nod = 1, cv_nonods
          cty_nod = MAP_DG2CTY(dg_nod)
          resid_cty(cty_nod) = resid_cty(cty_nod)+rhs(dg_nod)
       end do

       ! Course grid solver...
       DP_SMALL = 0.
       EWRITE(3,*)'SOLVER'
       CALL SOLVER( CMC_SMALL(1:NCMC_SMALL), DP_SMALL, resid_cty, &
            FINDCMC_SMALL, COLCMC_SMALL(1:NCMC_SMALL), &
            option_path = '/material_phase[0]/scalar_field::Pressure')
       EWRITE(3,*)'OUT OF SOLVER'

       ! Map the corrections DP_SMALL to dg:
       DO dg_nod = 1, cv_nonods
          cty_nod = MAP_DG2CTY(dg_nod)
          P(DG_NOD) = DP_SMALL(CTY_NOD)
       end do


       IF(ONE_SOLVE_1IT) THEN ! A single dg pressure iteration
          if (.true.) then
             call set_solver_options(path, &
                  ksptype = "gmres", &
                  !pctype = "jacobi", & ! use this for P1DGP1DG
                  !pctype = "sor", & ! use this for P1DGP1DG
                  pctype = "none", &   ! use this for P1DGP2DG
                  rtol = 1.e-10, &
                  atol = 1.e-15, &
                  max_its =21)
          endif

       resid_dg = rhs
       ustep = 0.0
       do dg_nod = 1, cv_nonods
          DO COUNT = FINDCMC(dg_NOD), FINDCMC(dg_NOD+1) - 1
             col=COLCMC(COUNT)
             resid_dg(dg_nod) = resid_dg(dg_nod) - cmc(count) * P(Col)
          END DO
       end do


             call add_option( &
                  trim(path)//"/solver/ignore_all_solver_failures", stat)
!             CALL SOLVER( CMC, P, RHS, &
             CALL SOLVER( CMC, ustep, resid_dg, &
                  FINDCMC, COLCMC, &
                  option_path = path )

            p=p+ustep
       ENDIF

       RETURN


    END IF




          if (.true.) then
             call set_solver_options(path, &
                  ksptype = "gmres", &
                  !pctype = "jacobi", & ! use this for P1DGP1DG
                  pctype = "sor", & ! use this for P1DGP1DG
                  !pctype = "none", &   ! use this for P1DGP2DG
                  rtol = 1.e-10, &
                  atol = 1.e-15, &
                  max_its = 7)
             !     max_its = 25)
          endif


    DO GL_ITS = 1, NGL_ITS

       EWRITE(3,*)'GL_ITS=', GL_ITS
       ! SSOR smoother for the multi-grid method...
       ewrite(3,*)'before solving:',p

       if (gl_its>1) then

          if (.true.) then

             !call set_solver_options(path, &
             !     ksptype = "gmres", &
             !     pctype = "hypre", &
             !     rtol = 1.e-10, &
             !     atol = 0., &
             !     max_its = 25)
             !call add_option( &
             !     trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", stat)
             !call set_option( &
             !     trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", "boomeramg")

             ! ignore solver failures...

             ustep=p

             call add_option( &
                  trim(path)//"/solver/ignore_all_solver_failures", stat)
             CALL SOLVER( CMC, P, RHS, &
                  FINDCMC, COLCMC, &
                  option_path = path )


   if(.false.) then
! optimize smoother contribution to solution...
       DP_DG=p-ustep
       p=ustep

       resid_dg = rhs
       ustep = 0.0
       do dg_nod = 1, cv_nonods
          DO COUNT = FINDCMC(dg_NOD), FINDCMC(dg_NOD+1) - 1
             col=COLCMC(COUNT)
             resid_dg(dg_nod) = resid_dg(dg_nod) - cmc(count) * P(Col)
             ustep(dg_nod) = ustep(dg_nod) + cmc(count) * dP_DG(COL)
          END DO
       end do

       OPT_STEP=-SUM(-USTEP(:)*RESID_DG(:))/MAX(1.E-15, SUM(USTEP(:)*USTEP(:)))
! Make sure the step length is between [0,1]
       OPT_STEP=MIN(1.0,MAX(0.,OPT_STEP))
        print *,'after smoother gl_its, OPT_STEP:',gl_its, OPT_STEP

       P = P + DP_DG * OPT_STEP
    endif



          end if

          if (.false.) then
             CALL SIMPLE_SOLVER( CMC, P, RHS,  &
                  NCOLCMC, CV_NONODS, FINDCMC, COLCMC, MIDCMC,  &
                  ERROR, RELAX, RELAX_DIAABS, RELAX_DIA, N_LIN_ITS )
             ewrite(3,*)'after solving:',p
          end if

          if (.false.) then
             CALL SOLVER( CMC, P, RHS, &
                  FINDCMC, COLCMC, &
                  option_path = '/material_phase[0]/scalar_field::Pressure' )
          end if

       end if


       resid_dg = rhs
       do dg_nod = 1, cv_nonods
          DO COUNT = FINDCMC(dg_NOD), FINDCMC(dg_NOD+1) - 1
             col=COLCMC(COUNT)
             resid_dg(dg_nod) = resid_dg(dg_nod) - cmc(count) * P(Col)
          END DO
       end do
       ! Map resid_dg to resid_cty as well as the solution:
       resid_cty=0.
       nods_sourou=0.
       do dg_nod = 1, cv_nonods
          cty_nod = MAP_DG2CTY(dg_nod)
          resid_cty(cty_nod) = resid_cty(cty_nod)+resid_dg(dg_nod)
          nods_sourou(cty_nod) = nods_sourou(cty_nod)+1.
       end do
       ! We have added the rows together so no need to normalize residual. 
       ! resid_cty= resid_cty/nods_sourou

       ! Course grid solver...
       DP_SMALL = 0.
       EWRITE(3,*)'SOLVER'
       CALL SOLVER( CMC_SMALL(1:NCMC_SMALL), DP_SMALL, resid_cty, &
            FINDCMC_SMALL, COLCMC_SMALL(1:NCMC_SMALL), &
            option_path = '/material_phase[0]/scalar_field::Pressure')
       EWRITE(3,*)'OUT OF SOLVER'

       ! Map the corrections DP_SMALL to dg:
       DO dg_nod = 1, cv_nonods
          cty_nod = MAP_DG2CTY(dg_nod)
          DP_DG(DG_NOD) = DP_SMALL(CTY_NOD)
       END DO

    OPT_STEP=1.0
    if(.true.) then
! Determine optimal step length...
       USTEP=0.0
       do dg_nod = 1, cv_nonods
          DO COUNT = FINDCMC(dg_NOD), FINDCMC(dg_NOD+1) - 1
             USTEP(dg_nod) = USTEP(dg_nod) + cmc(count) * DP_DG(COLCMC(COUNT))
          END DO
       end do
       OPT_STEP=-SUM(-USTEP(:)*RESID_DG(:))/MAX(1.E-15, SUM(USTEP(:)*USTEP(:)))
! Make sure the step length is between [0,1]
       OPT_STEP=MIN(1.0,MAX(0.,OPT_STEP))
        print *,'after compact solve gl_its, OPT_STEP:',gl_its, OPT_STEP
     endif

       P = P + DP_DG * OPT_STEP

    END DO

    deallocate( FINDcmc_small, &
         colcmc_small, &
         midcmc_small, &
         MAP_DG2CTY )

    deallocate( cmc_small, &
         resid_dg, &
         resid_cty, &
         dp_small, &
         dp_dg, &
         nods_sourou )

    RETURN
  END SUBROUTINE PRES_DG_MULTIGRID

  SUBROUTINE CMC_Agglomerator_solver(state, cmc_petsc, deltap, RHS_p, &
  NCOLCMC, CV_NONODS, FINDCMC, COLCMC, MIDCMC, &
  totele, cv_nloc, x_nonods, x_ndgln,  option_path)
      !
      ! Solve CMC * P = RHS for RHS.
      ! form a discontinuous pressure mesh for pressure...
      implicit none
      INTEGER, intent( in ) ::  NCOLCMC, CV_NONODS, totele, cv_nloc, x_nonods
      ! IGOT_CMC_PRECON=1 or 0 (1 if we have a preconditioning matrix)
      type( state_type ), dimension( : ), intent( inout ) :: state
      type(petsc_csr_matrix), intent(inout)::  CMC_petsc
      type( scalar_field ), intent(inout) :: deltap
      type( scalar_field ), intent(in) :: rhs_p
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
      INTEGER, DIMENSION( : ), intent( in ) :: MIDCMC
      INTEGER, DIMENSION( : ), intent( in ) :: x_ndgln
      character(len=*), intent(in) :: option_path


      !Local variables
      logical, parameter :: DG_correction = .true.! "Post-smoothing" using the DG formulation
      integer, parameter :: multi_level_its = 1!Solve for the residual and then run a couple of iterations of the DG pressure
      !With more than 1 iteration it does not work very well...
      character(len=OPTION_PATH_LEN) :: solv_options = "/tmp/pressure"
      integer, dimension( : ), allocatable :: findcmc_small, colcmc_small, midcmc_small, &
      MAP_DG2CTY
      integer :: ele,cv_iloc, dg_nod, cty_nod, jcolcmc, jcolcmc_small
      integer :: mx_ncmc_small, ncmc_small, count, count2, count3, GL_ITS, col
      real :: opt_step, k

      integer :: ierr, i
      real, dimension(1) :: auxR
      real, dimension(CV_NONODS) :: DP_DG, USTEP
      !Variables for CMC_Small_petsc
      type(petsc_csr_matrix)::  CMC_Small_petsc
      type(mesh_type), pointer :: pmesh
      type( scalar_field ) :: deltap_small, RHS_small, Residual_DG
      integer, dimension( x_nonods ) :: dnnz
      !type(halo_type), pointer :: halo

      ! obtain sparcity of a new matrix
      mx_ncmc_small = ncolcmc * 4
      allocate( FINDcmc_small(x_nonods+1) )
      allocate( colcmc_small(mx_ncmc_small) )
      allocate( midcmc_small(x_nonods) )
      allocate( MAP_DG2CTY(cv_nonods) )

      ! lump the pressure nodes to take away the discontinuity...
      DO ELE = 1, TOTELE
          DO CV_ILOC = 1, CV_NLOC
              dg_nod = (ele-1) * cv_nloc + cv_iloc
              cty_nod = x_ndgln( (ele-1) * cv_nloc + cv_iloc)
              MAP_DG2CTY(dg_nod) = cty_nod
          END DO
      END DO

      CALL GET_SPAR_CMC_SMALL(FINDCMC_SMALL, COLCMC_SMALL, MIDCMC_SMALL, &
      MX_NCMC_SMALL, NCMC_SMALL, CV_NONODS, X_NONODS, MAP_DG2CTY, &
      FINDCMC, COLCMC, NCOLCMC)



      !Prepare variables related to CMC_Small_PETSC
      pmesh => extract_mesh(state(1), "PressureMesh_Continuous")
      !halo => pmesh%halos(1)
      ! find the number of non zeros per row
      do i = 1, size( dnnz )
          dnnz( i ) =FINDCMC_SMALL( i+1 ) - FINDCMC_SMALL( i )
      end do

      call allocate( CMC_Small_petsc, size(dnnz), size(dnnz), dnnz, dnnz,(/1, 1/)&
      ,name = 'CMC_Small_petsc')!, halo = halo)
      !Allocate P_small and rhs_small
      call allocate(deltap_small,pmesh,"deltap_small")
      call allocate(RHS_small,pmesh,"RHS_small")
      call zero(CMC_Small_petsc)!; call zero(deltap_small); call zero(RHS_small)
      !Allocate Residual_DG
      pmesh => extract_mesh(state(1), "PressureMesh_Discontinuous")
      call allocate(Residual_DG,pmesh,"Residual_DG")

      !We create CMC_Small_petsc from CMC_petsc
      DO dg_nod = 1, CV_NONODS
          cty_nod = MAP_DG2CTY(dg_nod)
          ! add row dg_nod to row cty_nod of cty mesh
          DO COUNT = FINDCMC(DG_NOD), FINDCMC(DG_NOD+1) - 1
              jcolcmc_small = MAP_DG2CTY(COLCMC(COUNT))
              count2 = 0
              DO COUNT3 = FINDCMC_small(cty_NOD), FINDCMC_small(cty_NOD+1) - 1
                  if(colcmc_small(count3)==jcolcmc_small) count2=count3
              end do
              call MatGetValues(cmc_petsc%M, 1, (/ cmc_petsc%row_numbering%gnn2unn(dg_nod,1) /), 1,&
                 (/cmc_petsc%column_numbering%gnn2unn(COLCMC(COUNT),1) /),  auxR, ierr)
              call addto( CMC_Small_petsc, blocki = 1, blockj = 1, i = cty_nod, j = colcmc_small(count2),val = auxR(1))
          END DO
      END DO
      call assemble( CMC_Small_petsc )

      !call MatView(cmc_petsc%M,PETSC_VIEWER_STDOUT_SELF)
      !call MatView(CMC_Small_petsc%M,PETSC_VIEWER_STDOUT_SELF)


      if (DG_correction) then !Configure the discontinuous solver

          call get_option( '/material_phase[0]/scalar_field::Pressure/prognostic/' // &
          'solver/max_iterations', i)
          call set_solver_options(solv_options, &
          ksptype = "gmres", &
          pctype = "none", &!the pc has to be one that does not start from zero
          rtol = 1.e-10, &
          atol = 1.e-15, &
          max_its = 3)!Since it is a smoother, only 3 iterations are needed

          !Make the GMRES method more robust than usually
          call add_option( trim(solv_options)//"/solver/iterative_method/restart", ierr)
          call set_option( trim(solv_options)//"/solver/iterative_method/restart", 200)
          !Ignore solver failures since we are not looking for its convergence
          call add_option(trim(solv_options)//"/solver/ignore_all_solver_failures", ierr)
!          !Plot the residual, for debugging purposes only
!          call add_option(trim(solv_options)//"/solver/diagnostics/monitors/preconditioned_residual", ierr)
      end if

      !Multi-level solver
      do k = 1, multi_level_its

          !Calculate the residual
          if (multi_level_its > 1) then
              do dg_nod = 1, cv_nonods
                  DO COUNT = FINDCMC(dg_NOD), FINDCMC(dg_NOD+1) - 1
                      col=COLCMC(COUNT)
                      cty_nod = MAP_DG2CTY(dg_nod)
                      call MatGetValues(cmc_petsc%M, 1, (/ cmc_petsc%row_numbering%gnn2unn(dg_nod,1) /), 1,&
                       (/ cmc_petsc%column_numbering%gnn2unn(COLCMC(COUNT),1) /),  auxR, ierr)
                      !Project residual to RHS_Small
                      Residual_DG%val(dg_nod) = rhs_p%val(dg_nod) - auxR(1) * deltap%val(Col)
                  END DO
              end do
          end if

          !Restric the residual to the CG mesh
          call zero(RHS_small)
          do dg_nod = 1, cv_nonods
              cty_nod = MAP_DG2CTY(dg_nod)
              if (k > 1) then
                  RHS_small%val(cty_nod) = RHS_small%val(cty_nod) + Residual_DG%val(dg_nod)
              else!If first itetation, or most likely only one iteration, the residual is the RHS
                  RHS_small%val(cty_nod) = RHS_small%val(cty_nod) + RHS_p%val(dg_nod)
              end if
          end do

          !We solve the system (equivalent to solve in the coarsest mesh)
          call zero(deltap_small)
          call petsc_solve(deltap_small, CMC_Small_petsc, RHS_small, trim(option_path))


          !We map back the results (equivalent to project the values from coarse to fine meshes)
          DP_DG = 0.
          do dg_nod = 1, cv_nonods
              cty_nod = MAP_DG2CTY(dg_nod)
              DP_DG(dg_nod) = DP_DG(dg_nod) + deltap_small%val(cty_nod)
!              deltap%val(dg_nod) = deltap%val(dg_nod) + deltap_small%val(cty_nod)
          end do

          opt_step=1.0
          !If we are iterating then try to get the best from the CG pressure mesh
          if(multi_level_its > 1) then
              ! Determine optimal step length...
              ustep=0.0
              do dg_nod = 1, cv_nonods
                  DO COUNT = FINDCMC(dg_NOD), FINDCMC(dg_NOD+1) - 1
                      call MatGetValues(cmc_petsc%M, 1, (/ cmc_petsc%row_numbering%gnn2unn(dg_nod,1) /), 1,&
                       (/ cmc_petsc%column_numbering%gnn2unn(COLCMC(COUNT),1) /),  auxR, ierr)
                      USTEP(dg_nod) = USTEP(dg_nod) + auxR(1) * DP_DG(COLCMC(COUNT))
                  END DO
              end do
              opt_step = - sum(-ustep(:) * Residual_DG%val(:) ) / max(1d-15, dot_product(ustep, ustep))
              ! Make sure the step length is between [0,1]
              opt_step=min(1.0,max(0.,opt_step))
              print *, OPT_STEP
          endif
          !Update the DG solution with a dumped value of CG
          deltap%val = deltap%val + DP_DG * opt_step

          !Depending on whether the CG solution is helping or not we either solve or continue
          if (opt_step < 1d-2) then
              !The CG solver cannot do anymore, solve what remains to be done
              call petsc_solve(deltap,CMC_petsc,RHS_p,trim(option_path))
              exit!Exit the loop
          else if (DG_correction) then !Some smoothings in the finest mesh to finish the loop
              call petsc_solve(deltap,CMC_petsc,RHS_p,trim(solv_options))
          end if

      end do
      call deallocate(CMC_Small_petsc)
      call deallocate(deltap_small)
      call deallocate(RHS_small)
      call deallocate(Residual_DG)
      deallocate( FINDcmc_small, colcmc_small, midcmc_small, MAP_DG2CTY )

      RETURN
  END SUBROUTINE CMC_Agglomerator_solver


  SUBROUTINE Fix_to_bad_elements(cmc_petsc, &
      NCOLCMC, FINDCMC, COLCMC, MIDCMC, &
      totele, p_nloc, P_NDGLN, Quality_list)
      !We add a term in the CMC matrix to move mass from low pressure to high pressure nodes
      !Only for linear pressure for the time being
      implicit none
      INTEGER, intent( in ) ::  NCOLCMC, totele, p_nloc
      type(petsc_csr_matrix), intent(inout)::  CMC_petsc
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC, COLCMC, MIDCMC, P_NDGLN
      type(bad_elements), DIMENSION(:), intent( in ) :: Quality_list
      !Local variables
      real, parameter :: alpha = 0.5d-4
      integer, dimension(2) :: counter
      integer :: i, i_node, j_node, ele, COUNT, P_ILOC, bad_node, k, ierr
      real :: auxR, adapted_alpha
      real, dimension(1) :: rescal
      logical :: nodeFound
      !Initialize variables
      i = 1

      !Depending on the angle of the element we consider different alphas
      !for 90 it is 0, for 135 it is 0.25 *alpha for 180 it is alpha
      adapted_alpha = alpha * ((Quality_list(i)%angle - 90.)/90.)**2.

      do while (Quality_list(i)%ele>0 .and. i < size(Quality_list))
          counter = 1
          ele = Quality_list(i)%ele
          !Bad node
          bad_node = P_NDGLN((ele-1) * p_nloc + Quality_list(i)%nodes(1))
          !We get the diagonal value to use it as a reference when adding the over-relaxation
          call MatGetValues(cmc_petsc%M, 1, (/ cmc_petsc%row_numbering%gnn2unn(bad_node ,1) /), 1, &
            (/ cmc_petsc%column_numbering%gnn2unn(bad_node ,1) /),  rescal, ierr)
          rescal(1) = adapted_alpha * rescal(1)
          DO P_ILOC = 1, P_NLOC
              i_node = P_NDGLN((ele-1) * p_nloc + P_ILOC)
              DO COUNT = FINDCMC(i_node), FINDCMC(i_node+1) - 1
                j_node = COLCMC(COUNT)

                !Check that we are touching the element we want to...
                nodefound = .false.
                do k = 1, size(Quality_list(i)%nodes)
                    if (P_NDGLN((ele-1) * p_nloc + Quality_list(i)%nodes(k)) == j_node) then
                        nodefound = .true.
                        exit
                    end if
                end do

                !...otherwise we cycle to the next node
                if (.not.nodefound) cycle
!                !Just modify the diagonals
!                !Diagonals we put a 1
!                if (i_node == j_node) then
!                  auxR = 1.0
!                else!not the diagonals
!                  auxR = 0.
!                end if

                !Add diffusion from bad node to neighbours
                if (i_node == bad_node) then!Row of the bad element
                    if ( j_node== i_node) then!diagonal of bad node
                        auxR = 1.0
                    else!not the diagonals
                        auxR = -Quality_list(i)%weights(counter(1))
                        counter(1) = counter(1) + 1
                    end if
                else
                    if (i_node == j_node) then!diagonal of the other nodes
                        auxR = 1.0
                    else if (j_node == bad_node) then!The column of the bad node
                        auxR = -Quality_list(i)%weights(counter(2))
                        counter(2) = counter(2) + 1
                    else
                        auxR = 0.0
                    end if
                end if
                !Add the new data to the matrix
                call addto( cmc_petsc, blocki = 1, blockj = 1, i = i_node, j = j_node,val = rescal(1) * auxR)
              end do
          end do
          !After setting values we have to re-assemble the matrix as we want to use MatGetValues
          call assemble( cmc_petsc )
          !Advance one element
          i = i + 1
      end do

      end subroutine Fix_to_bad_elements

  SUBROUTINE GET_SPAR_CMC_SMALL(FINDCMC_SMALL,COLCMC_SMALL,MIDCMC_SMALL, &
       MX_NCMC_SMALL,NCMC_SMALL, CV_NONODS,X_NONODS, MAP_DG2CTY, &
       FINDCMC,COLCMC,NCOLCMC)
    ! Form sparcity COLCMC_SMALL, FINDCMC_SMALL 
    ! from FINDCMC,COLCMC...
    ! It lumps the DG pressure matrix to a continuous pressure matrix...
    INTEGER, intent( in ) ::  MX_NCMC_SMALL,NCOLCMC, CV_NONODS,X_NONODS
    INTEGER, intent( inout ) ::  NCMC_SMALL
    INTEGER, DIMENSION(: ), intent( inout ) :: FINDCMC_SMALL
    INTEGER, DIMENSION( : ), intent( inout ) :: COLCMC_SMALL
    INTEGER, DIMENSION( : ), intent( inout ) :: MIDCMC_SMALL

    INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
    INTEGER, DIMENSION( : ), intent( in ) :: COLCMC

    INTEGER, DIMENSION( : ), intent( in ) :: MAP_DG2CTY
    ! Local variables
    integer, dimension( : ), allocatable :: MX_NODS_ROW_SMALL, NODS_ROW_SMALL, &
         FINDCMC_SMALL_mx
    INTEGER :: dg_nod,cty_nod,count,jcolcmc,jcolcmc_small,count2,count3

    allocate(MX_NODS_ROW_SMALL(x_nonods))
    allocate(NODS_ROW_SMALL(x_nonods))
    allocate(FINDCMC_SMALL_mx(x_nonods+1))

    MX_NODS_ROW_SMALL=0
    DO dg_nod=1,CV_NONODS
       cty_nod=MAP_DG2CTY(dg_nod)
       MX_NODS_ROW_SMALL(cty_nod)=min( MX_NODS_ROW_SMALL(cty_nod)  &
            +FINDCMC(DG_NOD+1)-FINDCMC(DG_NOD),   x_nonods)
    END DO
    FINDCMC_small_MX(1)=1
    do cty_nod=2,x_nonods+1
       FINDCMC_small_MX(cty_nod)=FINDCMC_small_MX(cty_nod-1)+mx_NODS_ROW_SMALL(cty_nod-1)
    end do
    ewrite(3,*)'MAP_DG2CTY:',MAP_DG2CTY
    ewrite(3,*)'mx_NODS_ROW_SMALL:',mx_NODS_ROW_SMALL
    ewrite(3,*)'FINDCMC_small_MX:',FINDCMC_small_MX

    NODS_ROW_SMALL=0
    COLCMC_SMALL=0
    DO dg_nod=1,CV_NONODS
       cty_nod=MAP_DG2CTY(dg_nod)
       ! add row dg_nod to row cty_nod of cty mesh
       DO COUNT=FINDCMC(DG_NOD),FINDCMC(DG_NOD+1)-1
          jcolcmc=COLCMC(COUNT)
          jcolcmc_small=MAP_DG2CTY(jcolcmc)

          count2=0
          DO COUNT3=FINDCMC_small_mx(CTY_NOD),FINDCMC_SMALL_mx(CTY_NOD)+NODS_ROW_SMALL(cty_nod)-1
             if(colcmc_small(count3)==jcolcmc_small) count2=count3
          end do
          if(count2==0) then ! then put coln in as we have not found it in row 
             NODS_ROW_SMALL(cty_nod)=NODS_ROW_SMALL(cty_nod)+1
             COLCMC_SMALL(FINDCMC_small_mx(CTY_NOD)+NODS_ROW_SMALL(cty_nod)-1)=jcolcmc_small
             if(cty_nod==1) then
                ewrite(3,*)'dg_nod,cty_nod,jcolcmc_small:',dg_nod,cty_nod,jcolcmc_small
             endif
          endif
       END DO
    END DO
    ewrite(3,*)'NODS_ROW_SMALL:',NODS_ROW_SMALL

    FINDCMC_small(1)=1
    do cty_nod=2,x_nonods+1
       FINDCMC_small(cty_nod)=FINDCMC_small(cty_nod-1)+NODS_ROW_SMALL(cty_nod-1)
    end do
    NCMC_SMALL=FINDCMC_small(x_nonods+1)-1

    ! Shrink up the pointer list COLCMC_SMALL: 
    COUNT=0
    do cty_nod=1,x_nonods
       DO COUNT2=FINDCMC_small_MX(cty_nod),FINDCMC_small_MX(cty_nod+1)-1
          IF(COLCMC_SMALL(COUNT2).NE.0) THEN
             COUNT=COUNT+1
             COLCMC_SMALL(COUNT)=COLCMC_SMALL(COUNT2)
          ENDIF
       END DO
    END DO
    ! Put in assending coln order in each row...
    do cty_nod=1,x_nonods
       ewrite(3,*)'cty_nod,FINDCMC_small(cty_nod),FINDCMC_small(cty_nod+1)-1:', &
            cty_nod,FINDCMC_small(cty_nod),FINDCMC_small(cty_nod+1)-1
       call ibubble2(COLCMC_SMALL(FINDCMC_small(cty_nod):FINDCMC_small(cty_nod+1)-1))
       ewrite(3,*)'COLCMC_SMALL(FINDCMC_small(cty_nod):FINDCMC_small(cty_nod+1)-1):', &
            COLCMC_SMALL(FINDCMC_small(cty_nod):FINDCMC_small(cty_nod+1)-1)
    end do
    ! Calculate MIDCMC_SMALL...
    do cty_nod=1,x_nonods
       DO COUNT=FINDCMC_small(cty_nod),FINDCMC_small(cty_nod+1)-1
          IF(COLCMC_SMALL(COUNT)==cty_nod) MIDCMC_SMALL(CTY_NOD)=COUNT
       END DO
    END DO

    RETURN
  END SUBROUTINE GET_SPAR_CMC_SMALL
         


    subroutine ibubble2(ivec)
      ! sort ivec in increasing order
      implicit none
      integer, dimension( : ), intent( inout ) :: ivec
      ! Local variables
      integer :: nvec, i, j, itemp

      nvec = size(ivec)

!        ewrite(3,*)'before ivec:',ivec

      do j = 1, nvec
         do i = 1, nvec - 1
            if ( ivec( i ) > ivec( i + 1 ) ) then
               itemp = ivec( i + 1 )
               ivec( i + 1 ) = ivec( i )
               ivec( i ) = itemp
            end if
         end do
      end do
!        ewrite(3,*)'after ivec:',ivec
      return
    end subroutine ibubble2


    SUBROUTINE SIMPLE_SOLVER( CMC, P, RHS,  &
         NCMC, NONODS, FINCMC, COLCMC, MIDCMC,  &
         ERROR, RELAX, RELAX_DIAABS, RELAX_DIA, N_LIN_ITS )
      !
      ! Solve CMC * P = RHS for RHS.
      ! RELAX: overall relaxation coeff; =1 for no relaxation. 
      ! RELAX_DIAABS: relaxation of the absolute values of the sum of the row of the matrix;
      !               - recommend >=2 for hard problems, =0 for easy
      ! RELAX_DIA: relaxation of diagonal; =1 no relaxation (normally applied). 
      ! N_LIN_ITS = no of linear iterations
      ! ERROR= solver tolerence between 2 consecutive iterations
      implicit none
      REAL, intent( in ) :: ERROR, RELAX, RELAX_DIAABS, RELAX_DIA
      INTEGER, intent( in ) ::  N_LIN_ITS, NCMC, NONODS
      REAL, DIMENSION( : ), intent( in ) ::  CMC
      REAL, DIMENSION( : ), intent( inout ) ::  P
      REAL, DIMENSION( : ), intent( in ) :: RHS
      INTEGER, DIMENSION( : ), intent( in ) :: FINCMC
      INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
      INTEGER, DIMENSION( : ), intent( in ) :: MIDCMC
      ! Local variables
      INTEGER :: ITS, ILOOP, ISTART, IFINI, ISTEP, NOD, COUNT
      REAL :: R, SABS_DIAG, RTOP, RBOT, POLD, MAX_ERR
      LOGICAL :: jacobi

      ewrite(3,*) 'In Solver'

      jacobi = .false. !.true.

      if(jacobi) then

      Loop_Non_Linear_Iter1: DO ITS = 1, N_LIN_ITS

            Loop_Nods1: DO NOD = 1,NONODS
               R =  CMC( MIDCMC( NOD ))*P(NOD)+ RHS( NOD )
               DO COUNT = FINCMC( NOD ), FINCMC( NOD + 1 ) - 1
                  R = R - CMC( COUNT ) * P( COLCMC( COUNT ))
               END DO
               RTOP = R 
               RBOT = CMC( MIDCMC( NOD ))
               P( NOD ) = RELAX * ( RTOP / RBOT ) + ( 1.0 - RELAX ) * P( NOD )
            END DO Loop_Nods1

         END DO Loop_Non_Linear_Iter1
     else

      Loop_Non_Linear_Iter: DO ITS = 1, N_LIN_ITS

         MAX_ERR = 0.0
         Loop_Internal: DO ILOOP = 1, 2
            IF( ILOOP == 1 ) THEN
               ISTART = 1
               IFINI = NONODS
               ISTEP = 1
            ELSE
               ISTART = NONODS
               IFINI = 1
               ISTEP = -1
            ENDIF

            Loop_Nods: DO NOD = ISTART, IFINI, ISTEP
               R = RELAX_DIA * CMC( MIDCMC( NOD )) * P( NOD ) + RHS( NOD )
               SABS_DIAG = 0.0
               DO COUNT = FINCMC( NOD ), FINCMC( NOD + 1 ) - 1
                  R = R - CMC( COUNT ) * P( COLCMC( COUNT ))
                  SABS_DIAG = SABS_DIAG + ABS( CMC( COUNT ))
               END DO
               RTOP = R + RELAX_DIAABS * SABS_DIAG * P( NOD )
               RBOT = RELAX_DIAABS * SABS_DIAG + RELAX_DIA * CMC( MIDCMC( NOD ))
               POLD = P( NOD )
               P( NOD ) = RELAX * ( RTOP / RBOT ) + ( 1.0 - RELAX ) * P( NOD )
               MAX_ERR = MAX( MAX_ERR, ABS( POLD - P( NOD )))
            END DO Loop_Nods
         END DO Loop_Internal

         IF( MAX_ERR < ERROR ) CYCLE

      END DO Loop_Non_Linear_Iter
        endif

      ewrite(3,*) 'Leaving Solver'

      RETURN
    END SUBROUTINE SIMPLE_SOLVER





    subroutine BoundedSolutionCorrections( state, packed_state, small_findrm, small_colm, StorageIndexes, cv_ele_type, &
        for_sat, IDs2CV_ndgln)

      implicit none
      ! This subroutine adjusts field_val so that it is bounded between field_min, field_max in a local way.
      ! The sparcity of the local CV connectivity is in: small_findrm, small_colm.
      ! ngl_its=max no of global iterations e.g. 100.
      ! error_tol = tolerance on the iterations.
      !
      ! nloc_its: This iteration is very good at avoiding spreading the modifications too far - however it can stagnate.
      ! nloc_its2: This iteration is very good at avoiding stagnating but does spread the modifcations far.
      ! us a single iteration because of this as default...
      ! nits_nod: iterations at a nod - this iteration is very good at avoiding spreading the modifications too far -
      ! however it can stagnate.
      integer, parameter :: nloc_its = 5, nloc_its2 = 1, nits_nod = 100, ngl_its = 500
      real, parameter :: w_relax = 0.5, error_tol = 1.0e-5

      type( state_type ), dimension( : ), intent( inout ) :: state
      type( state_type ), intent( inout ) :: packed_state
      integer, dimension( : ), intent( in ) :: small_findrm, small_colm
      integer, intent( in ) :: cv_ele_type
      integer, dimension( : ), intent( inout ) :: StorageIndexes
      logical, optional, intent(in) :: for_sat
      integer, optional, dimension(:) :: IDs2CV_ndgln
      ! local variables...
      type ( tensor_field ), pointer :: field, ufield
      ! ( ndim1, ndim2, cv_nonods )
      real, dimension( :, :, : ), allocatable :: field_dev_val, field_alt_val, field_min, field_max
      real, dimension( :, : ), allocatable :: scalar_field_dev_max, scalar_field_dev_min
      real, dimension( :, : ), allocatable :: r_min, r_max
      integer, dimension( :, : ), allocatable :: ii_min, ii_max
      real, dimension( : ), allocatable :: mass_cv, mass_cv_sur
      integer :: ndim1, ndim2, cv_nonods, i, j, knod, inod, jnod, count, ii, jj, loc_its, loc_its2, its, gl_its
      logical :: changed, changed_something
      real :: max_change, error_changed, max_max_error, scalar_field_dev, mass_off, alt_max, alt_min


      ! variables for cv_fem_shape_funs_plus_storage
      integer :: ndim, cv_ngi, cv_ngi_short, cv_nloc, cv_snloc, &
           &     u_nloc, u_snloc, scvngi, sbcvngi, nface, &
           &     totele, x_nonods, ele, iloc, jloc
      real :: mm
      real, dimension( : ), pointer :: cvweight, cvweight_short, scvfeweigh, sbcvfeweigh, &
           &                           sele_overlap_scale, detwei
      real, dimension( :, : ), pointer :: cvn, cvn_short, cvfen, cvfen_short, ufen, &
           &                              scvfen, scvfenslx, scvfensly, sufen, sufenslx, sufensly, &
           &                              sbcvn, sbcvfen, sbcvfenslx, sbcvfensly, sbufen, sbufenslx, sbufensly
      real, dimension( :, :, : ), pointer :: cvfenlx_all, cvfenlx_short_all, ufenlx_all, &
           &                                 scvfenlx_all, sufenlx_all, sbcvfenlx_all, sbufenlx_all
      logical, dimension( :, : ), allocatable :: u_on_face, ufem_on_face, &
           &                                     cv_on_face, cvfem_on_face
      integer, pointer :: ncolgpts
      integer, dimension( : ), pointer :: findgpts, colgpts, x_ndgln, cv_ndgln
      integer, dimension( :, : ), pointer :: cv_neiloc, cv_sloclist, u_sloclist
      logical :: quad_over_whole_ele, d1, d3, dcyl

      type(scalar_field) :: mass_cv_sur_halo
      type( vector_field ), pointer :: x
      real, dimension( : ), allocatable, target :: detwei2, ra2
      real, dimension( :, : , :), allocatable, target :: nx_all2
      real, target :: volume2
      real, dimension(:,:), pointer :: Immobile_fraction

      !field => extract_tensor_field( packed_state, "PackedComponentMassFraction" )


      if (present_and_true(for_sat)) then
        field => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
      else
        field => extract_tensor_field( packed_state, "PackedTemperature" )
      end if

      ndim1 = size( field%val, 1 ) ; ndim2 = size( field%val, 2 ) ; cv_nonods = size( field%val, 3 )

      ewrite(3,*) 'Bounding correction input: iphase, icomp, min, max:'
      do j = 1, ndim2
         do i = 1, ndim1
            ewrite(3,*) i, j, minval( field%val( i, j, : ) ), maxval( field%val( i, j, : ) )
         end do
      end do

      allocate( field_dev_val( ndim1, ndim2, cv_nonods ), field_alt_val( ndim1, ndim2, cv_nonods ) )
      allocate( field_min( ndim1, ndim2, cv_nonods ), field_max( ndim1, ndim2, cv_nonods ) )
      allocate( scalar_field_dev_max( ndim1, ndim2 ), scalar_field_dev_min( ndim1, ndim2 ) )
      allocate( r_min( ndim1, ndim2 ), r_max( ndim1, ndim2 ) )
      allocate( ii_min( ndim1, ndim2 ), ii_max( ndim1, ndim2 ) )
      allocate( mass_cv( cv_nonods ), mass_cv_sur( cv_nonods ) )

      call allocate(mass_cv_sur_halo,field%mesh,'mass_cv_sur_halo')

      call get_option( '/geometry/dimension', ndim )

      ufield => extract_tensor_field( packed_state, "PackedVelocity" )

      cv_nloc = ele_loc( field, 1 ) ; cv_snloc = face_loc( field, 1 )
      u_nloc = ele_loc( ufield, 1 ) ; u_snloc = face_loc( ufield, 1 )

      quad_over_whole_ele = .false.

      call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, quad_over_whole_ele )

      allocate( cv_on_face( cv_nloc, scvngi ),  cvfem_on_face( cv_nloc, scvngi ) )
      allocate( u_on_face( u_nloc, scvngi ), ufem_on_face( u_nloc, scvngi ) )

      call cv_fem_shape_funs_plus_storage( &
                                ! volume shape functions...
           ndim, cv_ele_type,  &
           cv_ngi, cv_ngi_short, cv_nloc, u_nloc, cvn, cvn_short, &
           cvweight, cvfen, cvfenlx_all, &
           cvweight_short, cvfen_short, cvfenlx_short_all, &
           ufen, ufenlx_all, &
                                ! surface of each cv shape functions...
           scvngi, cv_neiloc, cv_on_face, cvfem_on_face, &
           scvfen, scvfenslx, scvfensly, scvfeweigh, &
           scvfenlx_all,  &
           sufen, sufenslx, sufensly,  &
           sufenlx_all,  &
                                ! surface element shape funcs...
           u_on_face, ufem_on_face, nface, &
           sbcvngi, sbcvn, sbcvfen,sbcvfenslx, sbcvfensly, sbcvfeweigh, sbcvfenlx_all, &
           sbufen, sbufenslx, sbufensly, sbufenlx_all, &
           cv_sloclist, u_sloclist, cv_snloc, u_snloc, &
                                ! define the gauss points that lie on the surface of the cv...
           findgpts, colgpts, ncolgpts, &
           sele_overlap_scale, quad_over_whole_ele,&
           state, "Press_mesh" , storageindexes(1))
           !state, "bound", storageindexes( 35 ) )

      totele = ele_count( field )
      x_ndgln => get_ndglno( extract_mesh( state( 1 ), "PressureMesh_Continuous" ) )
      cv_ndgln => get_ndglno( extract_mesh( state( 1 ), "PressureMesh" ) )
      x_nonods = node_count( extract_mesh( state( 1 ), "PressureMesh_Continuous" ) )
      x => extract_vector_field( packed_state, "PressureCoordinate" )

      d1 = ( ndim == 1 ) ; d3 = ( ndim == 3 ) ; dcyl = .false.

      allocate( detwei2( cv_ngi_short ), ra2( cv_ngi_short ), nx_all2( ndim, cv_nloc, cv_ngi ) )

      mass_cv = 0.0
      do  ele = 1, totele

         call detnlxr( ele, x%val( 1, : ), x%val( 2, : ), x%val( 3, : ), &
              &        x_ndgln, totele, x_nonods, cv_nloc, cv_ngi_short, cvfen_short, &
              &        cvfenlx_short_all( 1, :, : ), cvfenlx_short_all( 2, :, : ), &
              &        cvfenlx_short_all( 3, :, : ), &
              &        cvweight_short, detwei2, ra2, volume2, d1, d3, dcyl, &
              &        nx_all2( 1, :, : ), nx_all2( 2, :, : ), nx_all2( 3, :, : ) )

         detwei => detwei2

         do iloc = 1, cv_nloc
            inod = cv_ndgln( ( ele - 1 ) * cv_nloc + iloc )
            do jloc = 1, cv_nloc
               mm = sum( cvn( iloc, : ) * cvn( jloc, : ) * detwei( : ) )
               mass_cv( inod ) = mass_cv( inod ) + mm
            end do
         end do

      end do

      mass_cv_sur = 0.0
      do inod = 1, cv_nonods
         if ( .not. node_owned( field, inod ) ) cycle
         do count = small_findrm( inod ), small_findrm( inod + 1 ) - 1
            jnod = small_colm( count )
            mass_cv_sur(inod) = mass_cv_sur(inod) + mass_cv( jnod )
         end do
      end do
! Obtain the halos of mass_cv_sur: 
      mass_cv_sur_halo%val(:)=mass_cv_sur(:)
      call halo_update(mass_cv_sur_halo)
      mass_cv_sur(:)=mass_cv_sur_halo%val(:)

      !Establish bounds
      if (present_and_true(for_sat)) then
        if (present(IDs2CV_ndgln)) then
            !Define the immobile fractions for each phase
            call get_var_from_packed_state(packed_state, immobile_fraction = immobile_fraction)
            do inod = 1, cv_nonods
                do ii = 1, ndim2!phases
                    field_min(:, ii, inod) = immobile_fraction(ii,IDs2CV_ndgln(inod))
                    field_max(:, ii, inod) = 1.0 - sum(immobile_fraction(:,IDs2CV_ndgln(inod)))&
                         + immobile_fraction(ii,IDs2CV_ndgln(inod))
                end do
            end do
        else
            ewrite(3,*) 'IDs2CV_ndgln not passed into BoundedSolutionCorrections for SAT, default bounds set'
            field_min = 0.0 ; field_max = 1.0
        end if
      else
        field_min = 0.0 ; field_max = 1.0
      end if

      do gl_its = 1, ngl_its

         ! This iteration is very good at avoiding spreading the modifications too far - however it can stagnate.
         max_change = 0.0
         do loc_its = 1, nloc_its
            changed_something = .false.
            do knod = 1, cv_nonods ! exclude the halo values of knod for parallel
               if ( .not. node_owned( field, knod ) ) cycle

               do its = 1, nits_nod

                  r_min = 0.0 ; r_max = 0.0
                  ii_min = 0 ; ii_max = 0
                  scalar_field_dev_max = 0.0 ; scalar_field_dev_min = 0.0

                  ! Find the max and min deviation from the limited values....
                  do count = small_findrm( knod ), small_findrm( knod + 1 ) - 1
                     jnod = small_colm( count )
                     if ( .not. node_owned( field, jnod ) ) cycle

                     do j = 1, ndim2
                        do i = 1, ndim1

                           if ( field%val( i, j, jnod ) > field_max( i, j, jnod ) ) then
                              scalar_field_dev = field%val( i, j, jnod ) - field_max( i, j, jnod )
                           else if ( field%val( i, j, jnod ) < field_min( i, j, jnod ) ) then
                              scalar_field_dev = field%val( i, j, jnod ) - field_min( i, j, jnod )
                           else
                              scalar_field_dev = 0.0
                           end if

                           if ( scalar_field_dev * mass_cv ( jnod ) > r_max( i, j ) ) then
                              scalar_field_dev_max( i, j ) = scalar_field_dev
                              r_max( i, j ) = scalar_field_dev * mass_cv( jnod )
                              ii_max( i, j ) = jnod
                           end if

                           if ( scalar_field_dev * mass_cv ( jnod ) < r_min( i, j ) ) then
                              scalar_field_dev_min( i, j ) = scalar_field_dev
                              r_min( i, j ) = scalar_field_dev * mass_cv( jnod )
                              ii_min( i, j ) = jnod
                           end if

                        end do
                     end do

                  end do ! do count = small_findrm( knod ), small_findrm( knod + 1 ) - 1


                  changed=.false.
                  ! Change the max and min limited deviation by sharing the deviation between them...
                  do j = 1, ndim2
                     do i = 1, ndim1

                        ii = ii_max( i, j )
                        jj = ii_min( i, j )

                        if ( ii /= 0 .and. jj /= 0 ) then
                           if ( abs( r_max( i, j ) ) > abs( r_min( i, j ) ) ) then
                              alt_max = ( r_max( i, j ) + r_min( i, j ) ) / mass_cv( ii )
                              alt_min = 0.0
                           else
                              alt_max = 0.0
                              alt_min = ( r_max( i, j ) + r_min( i, j ) ) / mass_cv( jj )
                           end if
                           max_change = max( max_change, abs( -scalar_field_dev_max( i, j ) + alt_max ) )
                           max_change = max( max_change, abs( -scalar_field_dev_min( i, j ) + alt_min ) )

                           field%val( i, j, ii ) = field%val( i, j, ii ) - scalar_field_dev_max( i, j ) + alt_max
                           field%val( i, j, jj ) = field%val( i, j, jj ) - scalar_field_dev_min( i, j ) + alt_min

                           changed = .true.
                           changed_something = .true.
                        end if

                     end do
                  end do

                  if ( .not. changed ) exit ! stop iterating and move onto next node/CV...
               end do ! do its=1,nits_nod

            end do ! do knod = 1, cv_nonods
            if ( .not. changed_something ) exit ! stop iterating and move onto next stage of iteration...
         end do ! do loc_its=1,nloc_its

         call halo_update( field )


         ! This iteration is very good at avoiding stagnating but does spread the modifcations far.
         ! use a single iteration because of this as default...
         do loc_its2 = 1, nloc_its2
            do knod = 1, cv_nonods
!               if ( .not. node_owned( field, knod ) ) cycle
               do j = 1, ndim2
                  do i = 1, ndim1
                     if ( field%val( i, j, knod ) > field_max( i, j, knod ) ) then
                        field_dev_val( i, j, knod ) = field%val( i, j, knod ) - field_max( i, j, knod )
                     else if ( field%val( i, j, knod ) < field_min( i, j, knod ) ) then
                        field_dev_val( i, j, knod ) = field%val( i, j, knod ) - field_min( i, j, knod )
                     else
                        field_dev_val( i, j, knod ) = 0.0
                     end if
                  end do
               end do
            end do

            ! matrix vector...
            field_alt_val = 0.0
            do inod = 1, cv_nonods ! exclude the halo values of knod for parallel
               if ( .not. node_owned( field, inod ) ) cycle
               do count = small_findrm( inod ), small_findrm( inod + 1 ) - 1
                  jnod = small_colm( count )
                  mass_off = mass_cv( jnod ) / mass_cv_sur( jnod )
                  field_alt_val( :, :, inod ) = field_alt_val( :, :, inod ) + mass_off * field_dev_val( :, :, jnod )
               end do ! do count = small_findrm( inod ), small_findrm( inod + 1 ) - 1
            end do ! do inod = 1, cv_nonods

            ! w_relax\in[0,1]: - This relaxation is used because we
            ! have used a mass matrix which is not diagonally dominant
            ! =0.5 is suggested.
            ! =1.0 is no relaxation.
            field_alt_val = w_relax * field_alt_val + ( 1.0 - w_relax ) * field_dev_val

            ! adjust the values...
            error_changed = maxval( abs( -field_dev_val + field_alt_val ) )
            field%val( :, :, : ) = field%val( :, :, : ) - field_dev_val( :, :, : ) + field_alt_val( :, :, : )
            call halo_update( field )

         end do ! loc_its2

! communicate the errors ( max_change, error_changed ) ...
! this could be more efficient sending a vector...
          max_max_error= max( max_change, error_changed )
          call allmax( max_max_error ) 

         if ( max_max_error < error_tol ) exit

      end do ! gl_its

      ewrite(3,*) 'Bounding correction output: iphase, icomp, min, max:'
      do j = 1, ndim2
         do i = 1, ndim1
            ewrite(3,*) i, j, minval( field%val( i, j, : ) ), maxval( field%val( i, j, : ) )
         end do
      end do

      deallocate( field_dev_val, field_alt_val )
      deallocate( field_min, field_max )
      deallocate( scalar_field_dev_max, scalar_field_dev_min )
      deallocate( r_min, r_max )
      deallocate( ii_min, ii_max )
      deallocate( mass_cv, mass_cv_sur )

      deallocate( cv_on_face, cvfem_on_face, u_on_face, ufem_on_face )
      deallocate( detwei2, ra2, nx_all2 )

      call deallocate(mass_cv_sur_halo)
      return
    end subroutine BoundedSolutionCorrections


    subroutine Trust_region_correction(packed_state, sat_bak, Dumping_from_schema, CV_NDGLN, IDs2CV_ndgln, dumping)
    !In this subroutine we applied some corrections and dumping on the saturations obtained from the saturation equation
    !this idea is based on the paper SPE-173267-MS.
    !The method ensures convergence independent on the time step.
        implicit none
        !Global variables
        type( state_type ), intent(inout) :: packed_state
        integer, dimension(:) :: CV_NDGLN, IDs2CV_ndgln
        real, dimension(:, :), intent(in) :: sat_bak
        real, intent(in) :: Dumping_from_schema
        real, optional, intent(inout) :: dumping!this is the dumping_factor used
        !Local variables
        real, parameter :: offset_inflection = 0.1!This is what we allow the jump to surpass an inflection
        real, parameter :: offset_kink = 0.05!This is what we allow the jump to surpass a kink
        real, dimension(:, :), pointer :: dSat
        real :: dumping_factor
        !Parameters for the automatic dumping
        real, save :: Oldconvergence = -1
        real, save :: OldOldconvergence = -1
        real, save :: stored_dumping = 0.20
        real, save :: OldDumping = -1
        real, save :: OldOldDumping = -1
        real, save :: convergence_tol = -1
        !First, impose physical constrains
        call Set_Saturation_to_sum_one(packed_state, IDs2CV_ndgln)

        !We re-use the phase volume fraction to reduce the ram consumption
        call get_var_from_packed_state(packed_state,PhaseVolumeFraction = dSat)
        dSat = dSat - sat_bak


        if (Dumping_from_schema < 0.0) then!Automatic method based on the history of convergence
            !Retrieve convergence factor, to make sure that if between time steps things are going great, we do not reduce the
            !dumping_parameter
            if (convergence_tol< 0) &!retrieve it just once
                call get_option( '/timestepping/nonlinear_iterations/nonlinear_iterations_automatic',&
                     convergence_tol, default = -1. )

            dumping_factor = predictedDumping(stored_dumping, &
                OldDumping, OldOldDumping, dumping_in_sat, OldConvergence, OldOldConvergence, convergence_tol)
            !Update convegence history
            OldOldconvergence = Oldconvergence
            Oldconvergence = dumping_in_sat
            !Update dumping history
            OldOldDumping = OldDumping
            OldDumping = stored_dumping
            stored_dumping = dumping_factor

        else!Use the value introduced by the user
            dumping_factor = Dumping_from_schema
        end if


        !Update, if necessary, the dumping to make sure we don't overshoot. i.e: if looping we don't update over 1.0
        if (present(dumping)) then
            if (dumping + dumping_factor > 1.0) dumping_factor = 1.0 - dumping
            dumping = dumping_factor
        end if
        !***SATURATION DUMPING SECTION***
        !Calculate the new saturation (dSat is pointing to packed_state) using a global dumping parameter as it is more conservative
        !In the paper they say that local dumping could severily improve convergence...
        dSat = sat_bak + dumping_factor * dSat

    contains
        real function predictedDumping(Dumping, OldDumping, OldOldDumping,&
             Convergence, OldConvergence, OldOldConvergence, convergence_tol)
            !This function calculates a dumping parameter based on the history of convergence
            !By getting the coefficients to get the curve y = Ax^2 + Bx + C and
            !getting the dumping that will minimize that curve
            implicit none
            real, intent(in) :: Dumping, OldDumping, OldOldDumping
            real, intent(in) :: Convergence, OldConvergence, OldOldConvergence, convergence_tol
            !Local variables
            real, dimension(3,3) :: A
            real, dimension(3) :: b

            !Specify initialy the previous value to make sure it always returns a value
            predictedDumping = Dumping

            !Set the system
            if (OldOldconvergence > 0) then
                A(1,1:3) = (/OldOldDumping**2,OldOldDumping, 1.0 /)
                A(2,1:3) = (/OldDumping**2,   OldDumping,    1.0 /)
                A(3,1:3) = (/Dumping**2,      Dumping,       1.0 /)
                b(1:3) = (/OldOldConvergence, OldConvergence,Convergence  /)
                !Solve the system to get the coefficients A, B and C
                call invert(A)
                b = matmul(A,b)
            end if

            if (  OldOldconvergence> 0 .and. .not. ISNAN(sum(b)) .and. b(1) > 0) then
                !Get the value that will give us the minimum convergence
                predictedDumping = -b(2)/(2.0*b(1))
            else if (OldDumping > 0 .and. OldConvergence > 0)then!Use a linear approach
                !Solve the system to get the coefficients B and C
                b(2)= (OldConvergence- Convergence)/(OldDumping - Dumping)
                b(3) = OldConvergence - b(2) * OldDumping
                if (b(2)> 0) then!Diverging, then minimum values are in the direction of OldConvergence
                    predictedDumping = max(max(Dumping / 1.1, 0.09),b(2) * (OldConvergence&
                             + (OldConvergence- Convergence)/abs(b(2)) ) + b(3))!The search range is scalated by the tangent
                else!Converging, then minimum values are in the direction of Convergence
                    predictedDumping = min(min(Dumping * 1.05, 0.55), b(2) * (Convergence &
                            + (Convergence - OldConvergence)/abs(b(2))) + b(3))!The search range is scalated by the tangent
                end if
            else!Simplest method, increase or decrease previous dumping parameter
                if (Convergence-Oldconvergence < 0) then!Converging
                   predictedDumping = min(Dumping * 1.05, 0.55)
                else if (Convergence < convergence_tol) then
                    continue!Convergence factor very low, do not decrease the dumping parameter
                else!Diverging, the reduce with a minimum value that will mean performing all the non-linear iterations
                   predictedDumping = max(Dumping / 1.1, 0.09)
                end if
            end if

            !Make sure it is bounded
            predictedDumping = max(min(predictedDumping, 0.8), 1d-3)
        end function predictedDumping

    end subroutine Trust_region_correction

    subroutine Set_Saturation_to_sum_one(packed_state, IDs2CV_ndgln)
        !This subroutines eliminates the oscillations in the saturation that are bigger than a
        !certain tolerance and also sets the saturation to be between bounds
        Implicit none
        !Global variables
        type( state_type ), intent(inout) :: packed_state
        integer, dimension(:), intent(in) :: IDs2CV_ndgln
        !Local variables
        integer :: iphase, jphase, nphase, ele, cv_nod
        real :: maxsat, minsat, sum_of_phases
        real, dimension(:,:), pointer :: satura
        real, dimension(:, :), pointer :: Immobile_fraction

        call get_var_from_packed_state(packed_state, PhaseVolumeFraction = satura)
        !Get Immobile_fractions
        call get_var_from_packed_state(packed_state, Immobile_fraction = Immobile_fraction)

        nphase = size(satura,1)
        !Set saturation to be between bounds
        do cv_nod = 1, size(satura,2 )
            sum_of_phases = sum(satura(:,cv_nod))
            do iphase = 1, nphase
                minsat = Immobile_fraction(iphase, IDs2CV_ndgln(cv_nod))
                maxsat = 1 - sum(Immobile_fraction(:, IDs2CV_ndgln(cv_nod))) + minsat
                !We enforce the sum to one by spreading the error to all the phases
                if (sum_of_phases /= 1.0 ) &
                    satura(iphase, cv_nod) = satura(iphase, cv_nod) + (1.0 - sum_of_phases) / nphase
                satura(iphase,cv_nod) =  min(max(minsat, satura(iphase,cv_nod)),maxsat)
            end do
        end do

    end subroutine Set_Saturation_to_sum_one


    subroutine Set_Saturation_between_bounds(packed_state, IDs2CV_ndgln)
        !This subroutines eliminates the oscillations in the saturation that are bigger than a
        !certain tolerance
        Implicit none
        !Global variables
        type( state_type ), intent(inout) :: packed_state
        integer, dimension(:), intent(in) :: IDs2CV_ndgln
        !Local variables
        integer :: iphase, jphase, nphase, ele, cv_nod
        real :: maxsat, minsat
        real, dimension(:,:), pointer :: satura, Immobile_fraction

        call get_var_from_packed_state(packed_state, PhaseVolumeFraction = satura)
        !Get corey options
        call get_var_from_packed_state(packed_state, Immobile_fraction = Immobile_fraction)

        nphase = size(satura,1)
        !Set saturation to be between bounds
        do cv_nod = 1, size(satura,2 )
            do iphase = 1, nphase
                minsat = Immobile_fraction(iphase, IDs2CV_ndgln(cv_nod))
                maxsat = 1 - sum(Immobile_fraction(:, IDs2CV_ndgln(cv_nod))) + minsat
                !We enforce sat to be between bounds
                satura(iphase,cv_nod) =  min(max(minsat, satura(iphase,cv_nod)),maxsat)
            end do
        end do

    end subroutine Set_Saturation_between_bounds

end module solvers_module


