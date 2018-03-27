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

module cv_advection

    use fldebug

    use fields
    use solvers

    use reference_counting
    use memory_diagnostics

    use solvers_module
    use spud
    use global_parameters
    use futils, only: int2str
    use adapt_state_prescribed_module
    use sparse_tools
    use sparsity_patterns
    use surface_integrals, only :integrate_over_surface_element
    use petsc_tools
    use vtk_interfaces


    use shape_functions_Linear_Quadratic
    use shape_functions_NDim
    use shape_functions_prototype
    use matrix_operations
    use Copy_Outof_State
    use boundary_conditions

    use multi_interpolation
    use multi_tools
    use multi_data_types
    use multi_pipes
    use multiphase_EOS
#ifdef HAVE_PETSC_MODULES
  use petsc
#endif



    implicit none

    interface DG_DERIVS_ALL
        module procedure DG_DERIVS_ALL1
        module procedure DG_DERIVS_ALL2
    end interface DG_DERIVS_ALL

    !sprint_to_do; get rid of public variables

    public ::  mat1, mat2



    ! Variables needed for the mesh to mesh interpolation calculations

    type(csr_matrix) :: mat1
    type(csr_matrix) :: mat2

#include "petsc_legacy.h"

    INTEGER, PARAMETER :: WIC_T_BC_DIRICHLET = 1, WIC_T_BC_ROBIN = 2, &
        WIC_T_BC_DIRI_ADV_AND_ROBIN = 3, WIC_D_BC_DIRICHLET = 1, &
        WIC_U_BC_DIRICHLET = 1, &
        WIC_U_BC_ROBIN = 2, &
        WIC_U_BC_DIRI_ADV_AND_ROBIN = 3, &
        WIC_U_BC_DIRICHLET_INOUT = 2, &
        WIC_P_BC_DIRICHLET = 1, &
        WIC_P_BC_FREE = 2

contains

    SUBROUTINE CV_ASSEMB( state, packed_state, &
        Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd, &
        tracer, velocity, density, multi_absorp, &
        DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B,&
        DEN_ALL, DENOLD_ALL, &
        TDIFFUSION, IGOT_THERM_VIS, THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL, &
        CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, &
        SUF_SIG_DIAGTEN_BC, &
        DERIV, CV_P, &
        SOURCT_ALL, ABSORBT_ALL, VOLFRA_PORE, &
        GETCV_DISC, GETCT, &
        IGOT_T2, IGOT_THETA_FLUX, GET_THETA_FLUX, USE_THETA_FLUX, &
        THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
        MEAN_PORE_CV, &
        MASS_MN_PRES, THERMAL, &
        got_free_surf,  MASS_SUF, &
        MASS_ELE_TRANSP, &
        saturation,OvRelax_param, Phase_with_Pc, Courant_number,&
        Permeability_tensor_field, calculate_mass_delta, eles_with_pipe, pipes_aux, &
        porous_heat_coef, outfluxes, solving_compositional)
        !  =====================================================================
        !     In this subroutine the advection terms in the advection-diffusion
        !     equation (in the matrix and RHS) are calculated as ACV and CV_RHS.
        !
        !     This routine uses a Control Volume (CV) formulation to compute
        !     the advection terms. The general procedure is as follows:
        !
        !        1. For each node-pair, define which node is the donor, which is
        !           the receptor and define an "upwind" value of the field being
        !           advected (and the accompanying "density"; see note below)
        !        2. Calculate the volume flux across the CV face that separates
        !           these two nodes
        !        3. Estimate the value of the advected variable at the control
        !           volume face.
        !        4. Using information from the donor, receptor and upwind nodes,
        !           limit the field face value (removes oscillations from the
        !           solution)
        !        5. Assemble the fluxes to form the matrix and rhs of the
        !           advection equation
        !
        !     This procedure is implemented by considering the CV to be made up
        !     of a number of sub-control-volumes, which represent the part of
        !     the control volume within a given element.  The assembly of terms
        !     considers each of these sub-CVs in turn, calculating (and limiting)
        !     the flux across sub-CV faces that are external to the CV...
        !
        !     NOTE: Add in note about what density is in this sub!!!
        !
        !     To define the "upwind" value of the field variable, which is
        !     necessary for the limiting scheme, either:
        !
        !        A. The upwind value of the field variable to be advected is
        !           found by interpolation and stored in a matrix (TUPWIND)
        !        B. The neighbouring nodes are searched for the local maximum
        !           and minimum
        !
        !     The subroutine has several options...
        !
        !     Discretisation option
        !     ---------------------
        !      - The estimate of the face value may be determined in one of
        !        several ways.
        !      - The face value may be centered in time by either a specified
        !        CV_THETA value, or a non-linear CV_THETA value that is determined
        !        automatically.
        !      - The face value may be limited using a univeral-limiter-type
        !        scheme, or a limited-downwind scheme that is ideal for INTERFACE
        !        TRACKING.  Alternatively no limiting can be applied.
        !
        !     These options are defined by the value of CV_DISOPT, which corresponds
        !     to the clast digit of the GEM option NDISOT for the field in question.
        !
        !     CV_DISOPT=discretisation option in space and time
        !     ------------------------------------------------------------------
        !     CV_DISOPT   Method for face-value est.    Time-stepping     Limiting
        !     ------------------------------------------------------------------
        !       =0      1st order in space          Theta=specified    UNIVERSAL
        !       =1      1st order in space          Theta=non-linear   UNIVERSAL
        !       =2      Trapezoidal rule in space   Theta=specified    UNIVERSAL
        !       =2 if isotropic limiter then FEM-quadratic & stratification adjust. Theta=non-linear
        !       =3      Trapezoidal rule in space   Theta=non-linear   UNIVERSAL
        !       =4      Finite elements in space    Theta=specified    UNIVERSAL
        !       =5      Finite elements in space    Theta=non-linear   UNIVERSAL
        !       =6      Finite elements in space    Theta=specified    NONE
        !       =7      Finite elements in space    =non-linear   NONE
        !       =8      Finite elements in space    Theta=specified    DOWNWIND+
        !       =9      Finite elements in space    Theta=non-linear   DOWNWIND+
        !
        !     CV_DG_VEL_INT_OPT=interface velocity calculation option between elements
        !
        !     Limiting scheme
        !     ---------------
        !     The limiting scheme is defined in the subroutine NVDFUNNEW;
        !     the limited values are computed in subroutine ANONVDLIM/ONVDLIM.
        !
        !     ONVDLIM is the original limiting algorithm
        !
        !     ANONVDLIM is a new anisoptropic limiting algorithm, which is
        !     called if either ALOLIM=1 (where ALOLIM is an option flag set
        !     in this subroutine), or if the interface tracking limiting option
        !     is selected (CV_DISOPT=8/9).  ***In general ALOLIM appears to be set to 1 (GSC)
        !
        !     NOTE: ANONVDLIM only works for TETS; for all other element types
        !     ONVDLIM is used.
        !
        !
        !     IMPORTANT INPUTS:
        !     ----------------
        !
        !     CSR_ACV   - Matrix for assembling the advection terms (empty on input)
        !     CV_RHS      - Right-hand side vector for advection-diffusion terms
        !     X,Y,Z    - Node co-ordinates
        !     NU       - Nodal velocity component
        !     T,TOLD   - New and old advected field values at nodes
        !     DEN,  - New and old "density" at nodes, which is actually a constant
        !     DENOLD     multiplying the advection diffusion equation for the field
        !     CV_DISOPT   - The discretisation/limiting option (see above)
        !     DT       - The time step
        !     CV_THETA    - The time-stepping discretisation parameter
        !     CV_BETA     - Conservative(1.)/non-conservative(0.) flag
        !     ELE_TYP   - Integer flag definining element type
        !
        !     IMPORTANT OUTPUTS:
        !     -----------------
        !
        !     CSR_ACV   - Matrix updated to include the advection terms
        !     CV_RHS      - Right-hand side vector updated to include advection terms
        !
        !
        !     IMPORTANT LOCAL PARAMETERS:
        !     --------------------------
        !
        !     TIMOPT    - Temporal discretisation option, derived from CV_DISOPT.
        !                (1 for non-linear theta; 0 for theta specified (THETA))
        !
        !
        !***********************************************************************
        ! Inputs/Outputs
        IMPLICIT NONE
        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_GI_dimensions), intent(in) :: CV_GIdims
        type(multi_shape_funs), intent(inout) :: CV_funs
        type(multi_sparsities), intent(in) :: Mspars
        type(multi_ndgln), intent(in) :: ndgln
        type (multi_discretization_opts) :: Mdisopt
        type (multi_matrices), intent(inout) :: Mmat
        type (porous_adv_coefs), intent(inout) :: upwnd
        type(tensor_field), intent(inout), target :: tracer
        type(tensor_field), intent(in), target :: density
        type(tensor_field), intent(in) :: velocity
        type(multi_absorption), intent(inout) :: multi_absorp
        INTEGER, intent( in ) :: CV_DISOPT, CV_DG_VEL_INT_OPT, &
            IGOT_T2, IGOT_THETA_FLUX
        ! Diagonal scaling of (distributed) pressure matrix (used to treat pressure implicitly)
        REAL, DIMENSION( :, : ), intent( inout ), allocatable :: DIAG_SCALE_PRES
        REAL, DIMENSION( :, :, : ), intent( inout ), allocatable :: DIAG_SCALE_PRES_COUP ! (Mdims%npres, Mdims%npres, Mdims%cv_nonods)
        REAL, DIMENSION( :, :, : ), intent( inout ), allocatable :: INV_B ! (Mdims%nphase, Mdims%nphase, Mdims%cv_nonods)
        REAL, DIMENSION( : ), intent( inout ) :: MASS_MN_PRES
        REAL, DIMENSION( : ), intent( inout ) :: MASS_SUF
        REAL, DIMENSION( :, : ), target, intent( inout ) :: DEN_ALL
        REAL, DIMENSION( :, : ), intent( inout ) :: DENOLD_ALL
        REAL, DIMENSION( :, : ), intent( inout ) :: THETA_GDIFF ! (Mdims%nphase,Mdims%cv_nonods)
        REAL, DIMENSION( :, : ), intent( inout ), optional :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
        REAL, DIMENSION( :, :, :, : ), intent( in ) :: TDIFFUSION
        INTEGER, intent( in ) :: IGOT_THERM_VIS
        REAL, DIMENSION(:,:,:,:), intent( in ) :: THERM_U_DIFFUSION
        REAL, DIMENSION(:,:), intent( in ) :: THERM_U_DIFFUSION_VOL
        REAL, intent( in ) :: DT, CV_THETA, CV_BETA
        REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
        REAL, DIMENSION( :, : ), intent( in ) :: DERIV ! (Mdims%nphase,Mdims%cv_nonods)
        REAL, DIMENSION( :, :, : ), intent( in ) :: CV_P ! (1,Mdims%npres,Mdims%cv_nonods)
        REAL, DIMENSION( :, : ), intent( in) :: SOURCT_ALL
        REAL, DIMENSION( :, :, : ), pointer, intent( in ) :: ABSORBT_ALL
        REAL, DIMENSION( :, : ), intent( in ) :: VOLFRA_PORE ! (Mdims%npres,Mdims%totele)
        LOGICAL, intent( in ) :: GETCV_DISC, GETCT, GET_THETA_FLUX, USE_THETA_FLUX, THERMAL, got_free_surf
        ! got_free_surf - INDICATED IF WE HAVE A FREE SURFACE - TAKEN FROM DIAMOND EVENTUALLY...
        REAL, DIMENSION( :, : ), intent( inout ) :: MEAN_PORE_CV ! (Mdims%npres,Mdims%cv_nonods)
        REAL, DIMENSION( : ), intent( inout ), OPTIONAL  :: MASS_ELE_TRANSP
        type(tensor_field), intent(in), optional, target :: saturation
        !Variables for Capillary pressure
        real, optional, dimension(:), intent(in) :: OvRelax_param
        integer, optional, intent(in) :: Phase_with_Pc
        !Variables to cache get_int_vel OLD
        real, optional, dimension(:), intent(inout) :: Courant_number
        type( tensor_field ), optional, pointer, intent(in) :: Permeability_tensor_field
        ! Calculate_mass variable
        real, dimension(:,:), optional :: calculate_mass_delta
        type(pipe_coords), dimension(:), optional, intent(in):: eles_with_pipe
        type (multi_pipe_package), intent(in) :: pipes_aux
        REAL, DIMENSION( : , : ), optional, intent(in) :: porous_heat_coef
        ! Variable to store outfluxes
        type (multi_outfluxes), optional, intent(inout) :: outfluxes
        logical, optional, intent(in) :: solving_compositional
        ! ###################Local variables############################
        REAL :: ZERO_OR_TWO_THIRDS
        ! if integrate_other_side then just integrate over a face when cv_nodj>cv_nodi
        logical, PARAMETER :: integrate_other_side= .true.
        ! if .not.correct_method_petrov_method then we can compare our results directly with previous code...
        logical, PARAMETER :: correct_method_petrov_method= .true.
        ! IF GOT_CAPDIFFUS then add a diffusion term to treat capillary pressure term implicitly
        logical, PARAMETER :: GOT_CAPDIFFUS = .true.
        ! If UPWIND_CAP_DIFFUSION then when calculating capillary pressure diffusion coefficient use the
        ! upwind value
        logical, PARAMETER :: UPWIND_CAP_DIFFUSION = .true.
        ! THETA_VEL_HAT=0.0 does not change NDOTQOLD, THETA_VEL_HAT=1.0 sets NDOTQOLD=NDOTQNEW.
        ! If THETA_VEL_HAT<0.0 then automatically choose THETA_VEL to be as close to THETA_VEL_HAT (e.g.=0) as possible.
        ! This determins how implicit velocity is in the cty eqn. (from fully implciit =1.0, to do not alter the scheme =0.)
        ! Zhi try THETA_VEL_HAT = 1.0
        real :: THETA_VEL_HAT = 1.0
        ! if APPLY_ENO then apply ENO method to T and TOLD
        LOGICAL :: APPLY_ENO
        ! Mmat%CT will not change with this option...
        LOGICAL, PARAMETER :: CT_DO_NOT_CHANGE = .FALSE.,  PIPES_1D=.TRUE.
        ! PIPES_1D uses 1D pipes along edges of elements...
        ! GRAVTY is used in the free surface method only...
        REAL :: GRAVTY
        !
        logical, parameter :: EXPLICIT_PIPES= .false.
        logical, parameter :: EXPLICIT_PIPES2= .true.
        logical, PARAMETER :: MULTB_BY_POROSITY= .false.
        ! If GET_C_IN_CV_ADVDIF_AND_CALC_C_CV then form the Mmat%C matrix in here also based on control-volume pressure.
        ! if RECAL_C_CV_RHS, calculate the RHS forr the Mmat%C_CV matrix
        logical :: GET_C_IN_CV_ADVDIF_AND_CALC_C_CV
        REAL, PARAMETER :: FEM_PIPE_CORRECTION = 0.035
        ! FEM_PIPE_CORRECTION is the FEM pipe correction factor used because the Peacement
        ! model is derived for a 7-point 3D finite difference stencil. This correction factor is obtained
        ! by correlating the IC-FERST production results with Eclipse results on a regular mesh
        ! of linear tetrahedra elements and a single well at steady state. =0.0 is no correction =0.35 recommended.
        LOGICAL, DIMENSION( : ), allocatable :: X_SHARE
        INTEGER, DIMENSION( : ), allocatable :: &
            CV_OTHER_LOC, U_OTHER_LOC, MAT_OTHER_LOC, &
            JCOUNT_KLOC, JCOUNT_KLOC2, ICOUNT_KLOC, ICOUNT_KLOC2, &
            C_JCOUNT_KLOC, C_JCOUNT_KLOC2, C_ICOUNT_KLOC, C_ICOUNT_KLOC2, &
            CV_SLOC2LOC, U_SLOC2LOC
        INTEGER, DIMENSION( : , : ), allocatable :: FACE_ELE
        REAL, DIMENSION( : ), allocatable ::  &
            MASS_CV, MASS_ELE,  &
            SUM_CV, N, RSUM_VEC
        REAL, DIMENSION( :, : ), allocatable :: CVNORMX_ALL, XC_CV_ALL
        REAL, DIMENSION( :, :, : ), allocatable :: UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL
        REAL, DIMENSION( :, : ), allocatable :: CAP_DIFFUSION
        !###Variables for shape function calculation###
        type (multi_dev_shape_funs) :: SdevFuns
        !###Pointers for Shape function calculation###
        REAL, DIMENSION( : ), allocatable :: SHAPE_CV_SNL
        REAL, DIMENSION( :, :, : ), allocatable :: DUMMY_ZERO_NDIM_NDIM_NPHASE
        REAL, DIMENSION( :, :, :, : ), allocatable :: DTX_ELE_ALL, DTOLDX_ELE_ALL
        ! Variables used to calculate CV face values:
        REAL, DIMENSION( :, : ), allocatable :: LOC_F, LOC_FEMF
        REAL, DIMENSION( :, : ), allocatable :: SLOC_F, SLOC_FEMF, SLOC2_F, SLOC2_FEMF
        REAL, DIMENSION( :, :, : ), allocatable :: LOC_UF
        INTEGER, DIMENSION( : ), allocatable :: SELE_LOC_WIC_F_BC
        REAL, DIMENSION( :, : ), allocatable :: SLOC_SUF_F_BC
        REAL, DIMENSION( : ), allocatable :: FUPWIND_IN, FUPWIND_OUT, LIMF, F_INCOME, F_NDOTQ
        ! Variables used in GET_INT_VEL_NEW:
        REAL, DIMENSION ( :, :, : ), allocatable :: LOC_U, LOC2_U
        REAL, DIMENSION ( :, :, : ), allocatable :: LOC_NU, LOC2_NU, SLOC_NU, LOC_NUOLD, LOC2_NUOLD, SLOC_NUOLD
        REAL, DIMENSION ( :, : ), allocatable :: LOC_U_HAT, LOC2_U_HAT
        INTEGER :: CV_KNOD2, U_SNODK
        REAL, DIMENSION ( :, : ), allocatable :: LOC_FEMT, LOC2_FEMT, LOC_FEMTOLD, LOC2_FEMTOLD
        REAL, DIMENSION ( :, : ), allocatable :: LOC_FEMT2, LOC2_FEMT2, LOC_FEMT2OLD, LOC2_FEMT2OLD
        ! Mdims%nphase Variables:
        REAL, DIMENSION( : ), allocatable :: CAP_DIFF_COEF_DIVDX, DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX, &
            NDOTQNEW,   INCOME_J, LIMT2OLD, LIMDTOLD, &
            NDOTQ, INCOME, LIMT2, LIMTOLD, LIMT, LIMT_HAT, &
            LIMDOLD, LIMDTT2OLD,&
            FVT, FVT2, FVD, LIMD,  &
            LIMDT, LIMDTT2, Porous_diff_coef_divdx, Porous_diff_coefold_divdx
        LOGICAL :: DISTCONTINUOUS_METHOD, QUAD_ELEMENTS, use_reflect
        !Logical to check if we using a conservative method or not, to save cpu time
        logical :: conservative_advection
        !        ===> INTEGERS <===
        INTEGER :: COUNT, ICOUNT, JCOUNT, &
            ELE, ELE2, GI, GCOUNT, SELE,   &
            CV_SILOC, U_KLOC, &
            CV_ILOC, CV_JLOC, IPHASE, JPHASE, &
            CV_NODJ, ISWITCH, &
            CV_NODI, U_NODK, TIMOPT, &
            X_NODI,  X_NODJ, &
            CV_INOD, MAT_NODI,  MAT_NODJ, FACE_ITS, NFACE_ITS
        !        ===>  REALS  <===
        REAL :: HDC, &
            RSUM, &
            THERM_FTHETA, &
            W_SUM_ONE1, W_SUM_ONE2, h, rp, Skin, cc, one_m_cv_beta, auxR
        REAL :: FTHETA(Mdims%nphase), FTHETA_T2(Mdims%nphase), ONE_M_FTHETA_T2OLD(Mdims%nphase), FTHETA_T2_J(Mdims%nphase), ONE_M_FTHETA_T2OLD_J(Mdims%nphase)
        REAL :: ROBIN1(Mdims%nphase), ROBIN2(Mdims%nphase)
        integer :: IGETCT, IANISOLIM, global_face,J
        ! Functions...
        !REAL :: R2NORM, FACE_THETA
        !        ===>  LOGICALS  <===
        LOGICAL :: GETMAT, &
            D1, D3, GOT_DIFFUS, INTEGRAT_AT_GI, &
            NORMALISE, GET_GTHETA, QUAD_OVER_WHOLE_ELE
        LOGICAL, PARAMETER :: DCYL = .FALSE.
        character( len = option_path_len ) :: option_path2
        real, dimension(:,:), allocatable :: TUPWIND_MAT_ALL, TOLDUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, &
            DENOLDUPWIND_MAT_ALL, T2UPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL
        INTEGER :: IDUM(1)
        INTEGER :: I, IDIM, U_ILOC, ELE3, k
        INTEGER :: NFIELD, CV_KLOC, CV_NODK
        INTEGER :: IFI
        INTEGER :: COUNT_IN, COUNT_OUT,CV_KLOC2,CV_NODK2,CV_SKLOC
        INTEGER :: IPT_IN, IPT_OUT
        INTEGER :: U_KLOC2,U_NODK2,U_SKLOC
        INTEGER :: IPT,ILOOP,IMID,JMID,JDIM
        LOGICAL :: STORE, integrate_other_side_and_not_boundary, GOT_VIS
        REAL :: R, NDOTQ_HAT, DeltaP
        REAL , DIMENSION( : ), ALLOCATABLE :: F_CV_NODI, F_CV_NODJ
        REAL , DIMENSION( :, : ), ALLOCATABLE :: NUGI_ALL, NU_LEV_GI, SIGMA_INV_APPROX, SIGMA_INV_APPROX_NANO, opt_vel_upwind_coefs_new_cv
        REAL , DIMENSION( :, :, :, : ), ALLOCATABLE :: VECS_STRESS, VECS_GRAD_U
        REAL , DIMENSION( :, :, : ), ALLOCATABLE :: STRESS_IJ_THERM, STRESS_IJ_THERM_J
        REAL :: BCZERO(Mdims%nphase),  T_ALL_J( Mdims%nphase ), TOLD_ALL_J( Mdims%nphase )
        INTEGER :: LOC_WIC_T_BC_ALL(Mdims%nphase)
        REAL , DIMENSION( :, : ), allocatable :: NUOLDGI_ALL
        REAL, DIMENSION( : ), allocatable :: NDOTQOLD, INCOMEOLD
        REAL, DIMENSION( :, : ), ALLOCATABLE, target :: &
            FEMDEN_ALL, FEMDENOLD_ALL, FEMT2_ALL, FEMT2OLD_ALL
        LOGICAL, DIMENSION( : ), ALLOCATABLE :: DOWNWIND_EXTRAP_INDIVIDUAL
        LOGICAL, DIMENSION( :, : ), ALLOCATABLE :: IGOT_T_PACK, IGOT_T_CONST
        REAL, DIMENSION( :, : ), ALLOCATABLE :: IGOT_T_CONST_VALUE
        REAL, DIMENSION( : ), ALLOCATABLE :: ct_rhs_phase_cv_nodi, ct_rhs_phase_cv_nodj
        !Working variables
        real, dimension(:), allocatable :: VOL_FRA_FLUID ! for solid coupling
        real, dimension(:, :), allocatable :: U_HAT_ALL ! for solid coupling
        real, dimension(:,:), allocatable, target :: T_TEMP, TOLD_TEMP
        real, dimension(:,:), pointer :: T_ALL, TOLD_ALL, T2_ALL, T2OLD_ALL, X_ALL, FEMT_ALL, FEMTOLD_ALL
        real, dimension(:, :, :), pointer :: U_ALL, NU_ALL, NUOLD_ALL
        real, dimension(Mdims%nphase, Mdims%cv_nonods) :: T_ALL_KEEP
        real, dimension(:,:), allocatable :: MASS_CV_PLUS
        real, dimension( : ), allocatable :: DIAG_SCALE_PRES_phase
        real, dimension( : ), allocatable :: R_PEACMAN
        real, dimension( Mdims%nphase ) :: ct_rhs_phase
        real, dimension( : ), allocatable :: R_PRES,R_PHASE,CV_P_PHASE_NODI,CV_P_PHASE_NODJ,MEAN_PORE_CV_PHASE, MASS_PIPE_FOR_COUP
        real, dimension( :, :, : ), allocatable :: A_GAMMA_PRES_ABS,GAMMA_PRES_ABS2, PIPE_ABS
        !! boundary_condition fields
        type(tensor_field) :: velocity_BCs,tracer_BCs, density_BCs, saturation_BCs
        type(tensor_field) :: pressure_BCs
        type(tensor_field) :: tracer_BCs_robin2, saturation_BCs_robin2
        INTEGER, DIMENSION( 1 , Mdims%nphase , surface_element_count(tracer) ) :: WIC_T_BC_ALL, WIC_D_BC_ALL, WIC_T2_BC_ALL
        INTEGER, DIMENSION( Mdims%ndim , Mdims%nphase , surface_element_count(tracer) ) :: WIC_U_BC_ALL
        type( tensor_field ), pointer ::  pressure
        INTEGER, DIMENSION ( 1,Mdims%npres,surface_element_count(tracer) ) :: WIC_P_BC_ALL
        REAL, DIMENSION( :,:,: ), pointer :: SUF_T_BC_ALL,&
            SUF_T_BC_ROB1_ALL, SUF_T_BC_ROB2_ALL
        REAL, DIMENSION(:,:,: ), pointer :: SUF_D_BC_ALL,&
            SUF_T2_BC_ALL, SUF_T2_BC_ROB1_ALL, SUF_T2_BC_ROB2_ALL
        REAL, DIMENSION(:,:,: ), pointer :: SUF_U_BC_ALL
        REAL, DIMENSION( :,:,: ), allocatable, target :: SUF_T_BC,&
            SUF_T_BC_ROB1, SUF_T_BC_ROB2
        !Working variables for subroutines that are called several times
        real, dimension( Mdims%ndim,Mdims%nphase ) :: rdum_ndim_nphase_1
        real, dimension( Mdims%nphase ) :: rdum_nphase_1, rdum_nphase_2, rdum_nphase_3
        REAL, DIMENSION( Mdims%nphase ) :: ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA, &
                wrelax, XI_LIMIT, FEMTGI_IPHA, NDOTQ_TILDE, NDOTQ_INT, DT_J, abs_tilde, NDOTQ2, DT_I, LIMT3
        REAL, DIMENSION ( Mdims%ndim,Mdims%nphase ) :: UDGI_ALL, UDGI2_ALL, UDGI_INT_ALL, ROW_SUM_INV_VI, ROW_SUM_INV_VJ, UDGI_ALL_FOR_INV
        real, dimension( Mdims%nphase ) :: LOC_CV_RHS_I, LOC_CV_RHS_J, THETA_VEL
        type( vector_field ), pointer :: MeanPoreCV
        !! femdem
        type( vector_field ), pointer :: delta_u_all, us_all
        type( scalar_field ), pointer :: solid_vol_fra
        real :: theta_cty_solid, VOL_FRA_FLUID_I, VOL_FRA_FLUID_J
        type( tensor_field_pointer ), dimension(4+2*IGOT_T2) :: psi,fempsi
        type( vector_field_pointer ), dimension(1) :: PSI_AVE,PSI_INT
        type( tensor_field ), pointer :: old_tracer, old_density, old_saturation, tfield, temp_field
        integer :: FEM_IT
        integer, dimension(:), pointer :: neighbours
        integer :: nb, i_use_volume_frac_t2
        logical :: skip
        logical :: GOT_T2, use_volume_frac_T2
        logical :: symmetric_P
        ! pipe diamter for reservior modelling
        type( scalar_field ), pointer :: pipe_Diameter, pipe_Diameter_nano, pipe_Length_nano, conductivity_pipes, well_thickness
        logical :: has_conductivity_pipes
        !Permeability
        type( tensor_field ), pointer :: perm
        !Variables for Capillary pressure
        logical :: capillary_pressure_activated, between_elements, on_domain_boundary
        !Variables for get_int_vel_porous_vel
        logical :: anisotropic_and_frontier, anisotropic_perm
        real, dimension(Mdims%nphase):: rsum_nodi, rsum_nodj
        integer :: COUNT_SUF, P_JLOC, P_JNOD, stat, ipres, jpres
        REAL :: MM_GRAVTY
        !Variable to decide, for porous media, whether to consider, locally, using high order methods or not
        logical :: use_porous_limiter = .true., local_upwinding = .false.
        !Variables to calculate flux across boundaries
        logical :: calculate_flux
        real :: reservoir_P( Mdims%npres ) ! this is the background reservoir pressure
        real, dimension( :, :, : ), pointer :: fem_p
        integer :: U_JLOC
        real :: h_nano, RP_NANO, dt_pipe_factor
        logical :: got_nano

        logical :: have_absorption
        !Variables for Flooding
        real, dimension( : , : ), pointer :: DEN_ALL_DIVID
        logical, parameter :: implicit_fs = .false.
        real, parameter :: gravity_flooding = 9.80665
        real, parameter :: K_TOP = 1.0
        real :: fs_height, K_PIPES, l_surface_pipe, q_pipes, RDUM, RDUM2, CV_PIPE_LENGTH, l_frac
        type( tensor_field ), pointer :: bathymetry
        type( scalar_field ), pointer :: depth_of_drain
        real, dimension(Mdims%nphase):: SAT_FOR_PIPE, DEN_FOR_PIPE_PHASE
        real, dimension(Mdims%nphase, Mdims%nphase):: CONT_PIPE_ABS
        real, dimension(Mdims%nphase):: PRES_FOR_PIPE_PHASE, PRES_FOR_PIPE_PHASE_FULL
        !Variables for get int_vel_porous_vel
        logical :: iv_Incomming_flow
        REAL, DIMENSION(Mdims%ndim) :: iv_SUF_SIG_DIAGTEN_BC_GI
        INTEGER :: iv_u_kloc, iv_u_skloc, iv_cv_kloc, iv_idim, iv_CV_SKLOC, iv_CV_SNODK, iv_CV_SNODK_IPHA, iv_IPHASE, iv_u_kloc2
        real, dimension(Mdims%ndim, Mdims%ndim, Mdims%nphase) :: iv_aux_tensor, iv_sigma_aver, iv_aux_tensor2
        real, dimension(Mdims%ndim, Mdims%ndim) :: iv_ones
        !Variables for outfluxes
        real, dimension(:, :,:), allocatable :: bcs_outfluxes
        ! Variables needed when doing calculate_outfluxes
        real, dimension( : , : ), pointer ::Imble_frac
        ! Additions for calculating mass conservation - the total mass entering the domain is captured by 'bcs_outfluxes'
        ! and the internal changes in mass will be captured by 'calculate_mass_internal'
        real, allocatable, dimension(:) :: calculate_mass_internal  ! used in calculate_internal_mass subroutine
        real :: tmp1, tmp2, tmp3  ! Variables for parallel mass calculations


        !If on, then upwinding is used for the parts of the domain where there is no shock-front nor rarefaction
        local_upwinding = have_option('/numerical_methods/local_upwinding') .and. .not. present(solving_compositional)
        !this is true if the user is asking for high order advection scheme
        use_porous_limiter = (Mdisopt%in_ele_upwind /= 0)

        have_absorption=.false.
        if ( associated( absorbt_all ) ) have_absorption = .true.

        if ( Mdims%npres > 1 )then
            reservoir_P( 1 ) = 0.0 !1.0e+7
            reservoir_P( 2 ) = 0.0
        else
            reservoir_P = 0.0
        end if
        dt_pipe_factor = 1.0
        if ( Mdims%npres > 1 ) then
            ! Edge approach - pipe location and radius field
            ! this should really be (Mdims%npres, Mdims%cv_nonods)
            ! we assume one extra pressure for now
            pipe_Diameter => extract_scalar_field( state(1), "DiameterPipe" )
            got_nano = .false.
            if ( got_nano ) then
                pipe_Diameter_nano => extract_scalar_field( state(1), "DiameterPipeNano1" )
                pipe_Length_nano => extract_scalar_field( state(1), "LengthPipeNano1" )
            end if
            ! factor by which to reduce the pipe eqns time step size e.g. 10^{-3}
            call get_option( "/porous_media/well_option/dt_pipe_factor", dt_pipe_factor, default = 1.0 )


        end if
        !For thermal retrieve, if present, the conductivity of the pipes to calculate the heat loss
        has_conductivity_pipes = .false.
        if (thermal .and. is_porous_media) then
            has_conductivity_pipes = have_option('/wells_and_pipes/thermal_well_properties')
            if (has_conductivity_pipes) then
                conductivity_pipes => extract_scalar_field( state(1), "Conductivity" )
                well_thickness => extract_scalar_field( state(1), "well_thickness" )
            end if
        end if
        !Check pressure matrix based on Control Volumes
        !If we do not have an index where we have stored Mmat%C_CV, then we need to calculate it
        if (Mmat%CV_pressure) then
            GET_C_IN_CV_ADVDIF_AND_CALC_C_CV = .not.Mmat%stored !.true.
        else
            GET_C_IN_CV_ADVDIF_AND_CALC_C_CV = .false.
        end if
        symmetric_P = have_option( '/material_phase[0]/scalar_field::Pressure/prognostic/symmetric_P' )

        option_path2 = trim(tracer%option_path)//"/prognostic/spatial_discretisation/control_volumes/face_value::FiniteElement/limit_face_value/limiter::ENO"
        apply_eno = have_option( option_path2 )

        !THETA_VEL_HAT has to be zero for porous media flow
        if ( is_porous_media ) then
            THETA_VEL_HAT = 0.0
        else
            THETA_VEL_HAT = 1.0
        end if
        !Check capillary pressure options
        capillary_pressure_activated = .false.
        if (GOT_CAPDIFFUS) then
            if (present(OvRelax_param) .and. present(Phase_with_Pc)) then
                capillary_pressure_activated = Phase_with_Pc >0
            end if
        end if
        call get_option( "/physical_parameters/gravity/magnitude", gravty, stat )

        !#################SET WORKING VARIABLES#################

        call get_var_from_packed_state(packed_state,PressureCoordinate = X_ALL,&
            OldNonlinearVelocity = NUOLD_ALL, NonlinearVelocity = NU_ALL, FEPressure = FEM_P)
        if (is_porous_media .and. .not. present(solving_compositional)) call get_var_from_packed_state(packed_state, Immobile_fraction = Imble_frac)
        !For every Field_selector value but 3 (saturation) we need U_ALL to be NU_ALL
        U_ALL => NU_ALL
        old_tracer=>extract_tensor_field(packed_state,GetOldName(tracer))
        old_density=>extract_tensor_field(packed_state,GetOldName(density))
        if (present(saturation)) then
            old_saturation=>extract_tensor_field(packed_state,&
                GetOldName(saturation))
        end if
        T_ALL =>tracer%val(1,:,:)
        TOLD_ALL =>old_tracer%val(1,:,:)
        if (tracer%name == "PackedPhaseVolumeFraction") call get_var_from_packed_state(packed_state,Velocity = U_ALL)
        T_ALL_KEEP = T_ALL

            !##################END OF SET VARIABLES##################


        IF( GETCT ) THEN
            ! Initialise the calculate_mass variables
            !Allocate array to pass to store mass going through the boundaries
            if (allocated( outfluxes%outlet_id )) then
                allocate(bcs_outfluxes(Mdims%nphase, Mdims%cv_nonods, 0:size(outfluxes%outlet_id))); bcs_outfluxes= 0.!position zero is to store outfluxes over all bcs
            else
                allocate(bcs_outfluxes(Mdims%nphase, Mdims%cv_nonods, 0:1)); bcs_outfluxes= 0.!position zero is to store outfluxes over all bcs
            end if

            allocate ( calculate_mass_internal(Mdims%nphase))
            calculate_mass_internal(:) = 0.0  ! calculate_internal_mass subroutine
            !Extract temperature for outfluxes if required
            if (has_temperature) then
                temp_field => extract_tensor_field( packed_state, "PackedTemperature" )
                if (present(outfluxes)) then
                    if (outfluxes%calculate_flux)outfluxes%totout(2, :,:) = -273.15
                end if
            end if
        ENDIF



        !! Get boundary conditions from field
        call get_entire_boundary_condition(tracer,&
            ['weakdirichlet','robin        '],&
            tracer_BCs,WIC_T_BC_ALL,boundary_second_value=tracer_BCs_robin2)
        call get_entire_boundary_condition(density,&
            ['weakdirichlet'],&
            density_BCs,WIC_D_BC_ALL)
        if (present(saturation))&
            call get_entire_boundary_condition(saturation,&
            ['weakdirichlet','robin        '],&
            saturation_BCs,WIC_T2_BC_ALL,&
            boundary_second_value=saturation_BCs_robin2)
        call get_entire_boundary_condition(velocity,&
            ['weakdirichlet'],&
            velocity_BCs,WIC_U_BC_ALL)
        if(got_free_surf.or.is_porous_media .or. Mmat%CV_pressure) then
            pressure => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedFEPressure" )
            call get_entire_boundary_condition(pressure,&
                ['weakdirichlet','freesurface  '],&
                pressure_BCs,WIC_P_BC_ALL)
        endif
        !! reassignments to old arrays, to be discussed
        SUF_T_BC_ALL=>tracer_BCs%val
        SUF_T_BC_ROB1_ALL=>tracer_BCs%val ! re-using memory from dirichlet bc.s for Robin bc
        SUF_T_BC_ROB2_ALL=>tracer_BCs_robin2%val
        SUF_D_BC_ALL=>density_BCs%val
        SUF_U_BC_ALL=>velocity_BCs%val
        if(present(saturation)) then
            SUF_T2_BC_ALL=>saturation_BCs%val
            SUF_T2_BC_ROB1_ALL=>saturation_BCs%val ! re-using memory from dirichlet bc.s for Robin bc
            SUF_T2_BC_ROB2_ALL=>saturation_BCs_robin2%val
        end if
         if (tracer%name == "PackedTemperature" )  then
            allocate( suf_t_bc( 1,Mdims%nphase,Mdims%cv_snloc*Mdims%stotel ), suf_t_bc_rob1( 1,Mdims%nphase,Mdims%cv_snloc*Mdims%stotel ), &
                suf_t_bc_rob2( 1,Mdims%nphase,Mdims%cv_snloc*Mdims%stotel ) )
            call update_boundary_conditions( state, Mdims%stotel, Mdims%cv_snloc, Mdims%nphase, &
                suf_t_bc, suf_t_bc_rob1, suf_t_bc_rob2, tracer)
            SUF_T_BC_ALL=>suf_t_bc
            SUF_T_BC_ROB1_ALL=>suf_t_bc_rob1
            SUF_T_BC_ROB2_ALL=>suf_t_bc_rob2
        end if
        IDUM = 0
        ewrite(3,*) 'In CV_ASSEMB'
        GOT_VIS = .FALSE.
        IF(IGOT_THERM_VIS==1) GOT_VIS = ( R2NORM( THERM_U_DIFFUSION, Mdims%mat_nonods * Mdims%ndim * Mdims%ndim * Mdims%nphase ) /= 0 ) &
            .OR. ( R2NORM( THERM_U_DIFFUSION_VOL, Mdims%mat_nonods * Mdims%nphase ) /= 0 )
        GOT_DIFFUS = ( R2NORM( TDIFFUSION, Mdims%mat_nonods * Mdims%ndim * Mdims%ndim * Mdims%nphase ) /= 0 )

        call get_option( "/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/viscosity_scheme/zero_or_two_thirds", zero_or_two_thirds, default=2./3. )
        ewrite(3,*)'CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, GOT_DIFFUS:', &
            CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, GOT_DIFFUS
        ewrite(3,*)'GETCV_DISC, GETCT', GETCV_DISC, GETCT

        !For flooding, only for the saturation equation we need to take out the density
        !from the equations, we do so by forcing it to be 1, it needs to be recalculated
        !after this before it is used somewhere else
        if (is_flooding ) then
            if (GETCV_DISC) then
                DEN_ALL(1:mdims%n_in_pres,:) = 1.0
                DENOLD_ALL(1:mdims%n_in_pres,:) = 1.0
            end if

            if(mdims%npres > 1) then
                !As the density of phase 1 is used as height we use the densities defined for the phases
                !inside the wells to define the actual density
                DEN_FOR_PIPE_PHASE(1) = density%val(1,mdims%nphase - 1, 1)
                if (mdims%nphase > 1) DEN_FOR_PIPE_PHASE(2) = density%val(1,mdims%nphase, 1)
            else
                !If no pipes, then density is hard coded here
                DEN_FOR_PIPE_PHASE(1) = 1000.
                if (mdims%nphase > 1) DEN_FOR_PIPE_PHASE(2) = 1.0
            end if
        end if

        !cv_beta == 1 means conservative, meaning that everything multiplied by one_m_cv_beta can be ignored
        one_m_cv_beta = 1.0 - cv_beta
        conservative_advection = abs(one_m_cv_beta) <= 1e-8
        QUAD_OVER_WHOLE_ELE=.FALSE.
        ! Allocate memory for the control volume surface shape functions, etc.
        IF(GETCT) THEN
            ALLOCATE( JCOUNT_KLOC( Mdims%u_nloc ))
            ALLOCATE( JCOUNT_KLOC2( Mdims%u_nloc ))
            ALLOCATE( ICOUNT_KLOC( Mdims%u_nloc ))
            ALLOCATE( ICOUNT_KLOC2( Mdims%u_nloc ))
            IF(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV) THEN
                ALLOCATE( C_JCOUNT_KLOC( Mdims%u_nloc ))
                ALLOCATE( C_JCOUNT_KLOC2( Mdims%u_nloc ))
                ALLOCATE( C_ICOUNT_KLOC( Mdims%u_nloc ))
                ALLOCATE( C_ICOUNT_KLOC2( Mdims%u_nloc ))
            ENDIF
        ENDIF
        DISTCONTINUOUS_METHOD = ( Mdims%cv_nonods == Mdims%totele * Mdims%cv_nloc )
        ! Quadratic elements
        QUAD_ELEMENTS = ( ((Mdims%ndim==2).AND.(Mdims%cv_nloc==6)).or.((Mdims%ndim==3).AND.(Mdims%cv_nloc==10)) )
        !Pointer to permeability
        if ( is_porous_media ) then
            if (present(Permeability_tensor_field)) then
                perm => Permeability_tensor_field
            else
            perm=>extract_tensor_field(packed_state,"Permeability")
            end if
            !Check if the permeability is not isotropic and the method is DG
            anisotropic_perm = .not.have_option('porous_media/scalar_field::Permeability') .and. DISTCONTINUOUS_METHOD
        end if
        !Initialize Courant number for porous media
        if (present(Courant_number) .and. is_porous_media) Courant_number = 0.
        ALLOCATE( CVNORMX_ALL( Mdims%ndim, CV_GIdims%scvngi )) ; CVNORMX_ALL=0.0
        ALLOCATE( CV_OTHER_LOC( Mdims%cv_nloc ))
        ALLOCATE( U_OTHER_LOC( Mdims%u_nloc ))
        ALLOCATE( MAT_OTHER_LOC( Mdims%mat_nloc ))
        ALLOCATE( X_SHARE( Mdims%x_nonods ))
        ALLOCATE( SHAPE_CV_SNL( Mdims%cv_nloc ))
        ALLOCATE( DUMMY_ZERO_NDIM_NDIM_NPHASE(Mdims%ndim,Mdims%ndim,Mdims%nphase))
        DUMMY_ZERO_NDIM_NDIM_NPHASE=0.0
        ALLOCATE( CV_SLOC2LOC( Mdims%cv_snloc ))
        ALLOCATE( U_SLOC2LOC( Mdims%u_snloc ))
        ALLOCATE( UGI_COEF_ELE_ALL(Mdims%ndim,Mdims%nphase,Mdims%u_nloc) )
        ALLOCATE( UGI_COEF_ELE2_ALL(Mdims%ndim,Mdims%nphase,Mdims%u_nloc) )
        ! The procity mapped to the CV nodes
        ALLOCATE( SUM_CV( Mdims%cv_nonods ))
        D1 = ( Mdims%ndim == 1 )
        D3 = ( Mdims%ndim == 3 )
        !DCYL = ( Mdims%ndim == -2 )
        GETMAT = .TRUE.
        X_SHARE = .FALSE.
        ! Determine FEMT (finite element wise) etc from T (control volume wise)
        ! Also determine the CV mass matrix MASS_CV and centre of the CV's XC_CV_ALL
        ! This is for projecting to finite element basis functions...
        ALLOCATE( MASS_CV( Mdims%cv_nonods ))
        ALLOCATE( MASS_ELE( Mdims%totele ))
        ALLOCATE( XC_CV_ALL( Mdims%ndim, Mdims%cv_nonods ))
        ALLOCATE( DTX_ELE_ALL( Mdims%ndim, Mdims%nphase, Mdims%cv_nloc, Mdims%totele ))
        ALLOCATE( DTOLDX_ELE_ALL( Mdims%ndim, Mdims%nphase, Mdims%cv_nloc, Mdims%totele ))
        IGETCT = 0
        IF ( GETCT ) IGETCT = 1
        GOT_T2=( IGOT_T2 == 1 )
        use_volume_frac_T2=( IGOT_T2 == 1 .or. thermal )
        i_use_volume_frac_t2= 0
        if (use_volume_frac_T2) i_use_volume_frac_t2= 1
        nullify(FEMT_ALL); nullify(FEMTOLD_ALL);
        ALLOCATE( FEMT_ALL( Mdims%nphase, Mdims%cv_nonods ), FEMTOLD_ALL( Mdims%nphase, Mdims%cv_nonods ) )
        ALLOCATE( FEMDEN_ALL( Mdims%nphase, Mdims%cv_nonods ), FEMDENOLD_ALL( Mdims%nphase, Mdims%cv_nonods ) )
        ALLOCATE( FEMT2_ALL( Mdims%nphase, Mdims%cv_nonods ), FEMT2OLD_ALL( Mdims%nphase, Mdims%cv_nonods ) )
        allocate(CV_P_PHASE_NODI(Mdims%nphase),CV_P_PHASE_NODJ(Mdims%nphase))
        allocate(CT_RHS_PHASE_CV_NODI(Mdims%nphase),CT_RHS_PHASE_CV_NODJ(Mdims%nphase))

        IF ( GOT_T2 .OR. THERMAL) call get_var_from_packed_state( packed_state, &
            PhaseVolumeFraction = T2_ALL, OldPhaseVolumeFraction = T2OLD_ALL )
        ! variables for get_int_tden********************
        ! Set up the fields...
        ALLOCATE( IGOT_T_PACK( Mdims%nphase, 6 ), IGOT_T_CONST( Mdims%nphase, 6 ), IGOT_T_CONST_VALUE( Mdims%nphase, 6 ) )



        ! FOR packing as well as for detemining which variables to apply interface tracking**********
        !          STORE=.TRUE.
        STORE=.FALSE.
        IGOT_T_PACK=.TRUE.
        IGOT_T_CONST=.FALSE.
        IGOT_T_CONST_VALUE=0.0
        ! If we have any bc's then assume we ave a non-uniform field...
        DO IPHASE=1,Mdims%nphase
            IF( SUM(  WIC_T_BC_ALL( :, IPHASE, : ) ) == 0)  &
                CALL IS_FIELD_CONSTANT(IGOT_T_CONST(IPHASE,1), IGOT_T_CONST_VALUE(IPHASE,1), T_ALL(IPHASE,:),Mdims%cv_nonods)
        END DO
        DO IPHASE=1,Mdims%nphase
            IF( SUM(  WIC_T_BC_ALL( :, IPHASE, : ) ) == 0)  &
                CALL IS_FIELD_CONSTANT(IGOT_T_CONST(IPHASE,2), IGOT_T_CONST_VALUE(IPHASE,2), TOLD_ALL(IPHASE,:),Mdims%cv_nonods)
        END DO
        DO IPHASE=1,Mdims%nphase
            IF( SUM(  WIC_D_BC_ALL( :, IPHASE, : ) ) == 0)  &
                CALL IS_FIELD_CONSTANT(IGOT_T_CONST(IPHASE,3), IGOT_T_CONST_VALUE(IPHASE,3), DEN_ALL(IPHASE,:),Mdims%cv_nonods)
        END DO
        DO IPHASE=1,Mdims%nphase
            IF( SUM(  WIC_D_BC_ALL( :, IPHASE, : ) ) == 0)  &
                CALL IS_FIELD_CONSTANT(IGOT_T_CONST(IPHASE,4), IGOT_T_CONST_VALUE(IPHASE,4), DENOLD_ALL(IPHASE,:),Mdims%cv_nonods)
        END DO
        DO IPHASE=1,Mdims%nphase
            IF(use_volume_frac_t2) THEN
                IF( SUM(  WIC_T2_BC_ALL(:,  IPHASE, : ) ) == 0)  &
                    CALL IS_FIELD_CONSTANT(IGOT_T_CONST(IPHASE,5), IGOT_T_CONST_VALUE(IPHASE,5), T2_ALL(IPHASE,:),Mdims%cv_nonods)
            ELSE
                IGOT_T_CONST(IPHASE,5)=.TRUE.
                IGOT_T_CONST_VALUE(IPHASE,5)=1.0
            ENDIF
        END DO
        DO IPHASE=1,Mdims%nphase
            IF(use_volume_frac_t2) THEN
                IF( SUM(  WIC_T2_BC_ALL( :, IPHASE, : ) ) == 0)  &
                    CALL IS_FIELD_CONSTANT(IGOT_T_CONST(IPHASE,6), IGOT_T_CONST_VALUE(IPHASE,6), T2OLD_ALL(IPHASE,:),Mdims%cv_nonods)
            ELSE
                IGOT_T_CONST(IPHASE,6)=.TRUE.
                IGOT_T_CONST_VALUE(IPHASE,6)=1.0
            ENDIF
        END DO
        NFIELD=0
        DO IFI=1,6
            DO IPHASE=1,Mdims%nphase
                IF(.not.IGOT_T_CONST(IPHASE,IFI)) NFIELD=NFIELD+1
            END DO
        END DO
        ALLOCATE( DOWNWIND_EXTRAP_INDIVIDUAL( NFIELD ) )
        ! Determine IGOT_T_PACK(IPHASE,:):
        IGOT_T_PACK=.FALSE.
        DO ILOOP=1,6
            DO IPHASE=1,Mdims%nphase
                IF(.NOT.IGOT_T_CONST(IPHASE,ILOOP)) THEN
                    ! here we might check to see if we have this in the local storage...
                    IGOT_T_PACK(IPHASE,ILOOP)=.TRUE.
                ENDIF
            END DO
        END DO
        ! This logical needs to be expanded...
        DOWNWIND_EXTRAP_INDIVIDUAL = .FALSE.
        IF ( CV_DISOPT>=8 ) DOWNWIND_EXTRAP_INDIVIDUAL = .TRUE.
        ! F and LOC_U:
        ALLOCATE(LOC_F(NFIELD,Mdims%cv_nloc))
        ALLOCATE(LOC_FEMF(NFIELD,Mdims%cv_nloc))
        ALLOCATE(SLOC_F(NFIELD,Mdims%cv_snloc))
        ALLOCATE(SLOC_FEMF(NFIELD,Mdims%cv_snloc))
        ALLOCATE(SLOC2_F(NFIELD,Mdims%cv_snloc))
        ALLOCATE(SLOC2_FEMF(NFIELD,Mdims%cv_snloc))
        ALLOCATE(LOC_UF(Mdims%ndim,NFIELD,Mdims%u_nloc))
        ! Local variables used in GET_INT_VEL_NEW:
        ALLOCATE( LOC_U(Mdims%ndim, Mdims%nphase, Mdims%u_nloc), LOC2_U(Mdims%ndim, Mdims%nphase, Mdims%u_nloc) )
        ALLOCATE( LOC_NU(Mdims%ndim, Mdims%nphase, Mdims%u_nloc), LOC2_NU(Mdims%ndim, Mdims%nphase, Mdims%u_nloc) )
        ALLOCATE( SLOC_NU(Mdims%ndim, Mdims%nphase, Mdims%u_snloc) )
        ALLOCATE( LOC_NUOLD(Mdims%ndim, Mdims%nphase, Mdims%u_nloc), LOC2_NUOLD(Mdims%ndim, Mdims%nphase, Mdims%u_nloc) )
        ALLOCATE( SLOC_NUOLD(Mdims%ndim, Mdims%nphase, Mdims%u_snloc) )
        ALLOCATE( LOC_FEMT(Mdims%nphase, Mdims%cv_nloc), LOC2_FEMT(Mdims%nphase, Mdims%cv_nloc) )
        ALLOCATE( LOC_FEMTOLD(Mdims%nphase, Mdims%cv_nloc), LOC2_FEMTOLD(Mdims%nphase, Mdims%cv_nloc) )
        ALLOCATE( LOC_FEMT2(Mdims%nphase, Mdims%cv_nloc), LOC2_FEMT2(Mdims%nphase, Mdims%cv_nloc) )
        ALLOCATE( LOC_FEMT2OLD(Mdims%nphase, Mdims%cv_nloc), LOC2_FEMT2OLD(Mdims%nphase, Mdims%cv_nloc) )
        ! bc's:
        ALLOCATE( SELE_LOC_WIC_F_BC( NFIELD ) )
        ALLOCATE( SLOC_SUF_F_BC(NFIELD, Mdims%cv_snloc) )
        ! limiting values...
        ALLOCATE( FUPWIND_IN( NFIELD ) )
        ALLOCATE( FUPWIND_OUT( NFIELD ) )
        ! limiting and upwinding:
        ALLOCATE( LIMF( NFIELD ) )
        ALLOCATE( F_INCOME( NFIELD ), F_NDOTQ( NFIELD ) )
        !
        ALLOCATE( F_CV_NODI( NFIELD ), F_CV_NODJ( NFIELD ) )
        !
        ALLOCATE( NUGI_ALL( Mdims%ndim, Mdims%nphase ),  NUOLDGI_ALL(Mdims%ndim, Mdims%nphase), NU_LEV_GI( Mdims%ndim, Mdims%nphase ) )
        !
        IF(THERMAL.AND.GOT_VIS) THEN
            ALLOCATE( VECS_STRESS( Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%cv_nonods ), VECS_GRAD_U( Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%cv_nonods ) ) ; VECS_STRESS=0. ; VECS_GRAD_U=0.
            ALLOCATE( STRESS_IJ_THERM( Mdims%ndim, Mdims%ndim, Mdims%nphase ), STRESS_IJ_THERM_J( Mdims%ndim, Mdims%ndim, Mdims%nphase ) ) ; STRESS_IJ_THERM=0. ; STRESS_IJ_THERM_J=0.
        ENDIF
        ! NFIELD Variables:
        ALLOCATE(CAP_DIFF_COEF_DIVDX(Mdims%nphase), DIFF_COEF_DIVDX(Mdims%nphase), DIFF_COEFOLD_DIVDX(Mdims%nphase), &
            NDOTQNEW(Mdims%nphase), LIMT2OLD(Mdims%nphase), LIMDTOLD(Mdims%nphase), INCOMEOLD(Mdims%nphase), NDOTQOLD(Mdims%nphase), &
            NDOTQ(Mdims%nphase), INCOME(Mdims%nphase), LIMT2(Mdims%nphase), LIMTOLD(Mdims%nphase), LIMT(Mdims%nphase), LIMT_HAT(Mdims%nphase), &
            LIMDOLD(Mdims%nphase), LIMDTT2OLD(Mdims%nphase),&
            FVT(Mdims%nphase), FVT2(Mdims%nphase), FVD(Mdims%nphase), LIMD(Mdims%nphase),  &
            LIMDT(Mdims%nphase), LIMDTT2(Mdims%nphase)  )
        if (thermal .and.is_porous_media) allocate(Porous_diff_coef_divdx(Mdims%nphase), Porous_diff_coefold_divdx(Mdims%nphase))
        LIMT_HAT=0.0
        ALLOCATE(INCOME_J(Mdims%nphase))
        IF ( capillary_pressure_activated) THEN
            ALLOCATE( CAP_DIFFUSION( Mdims%nphase, Mdims%mat_nonods ) )
            !Introduce the information in CAP_DIFFUSION
            CAP_DIFFUSION = 0.!Initialize to zero just in case
            do ele = 1, Mdims%totele
                do CV_ILOC = 1, Mdims%cv_nloc
                    CV_NODI = ndgln%cv(CV_ILOC + (ele-1) * Mdims%cv_nloc)
                    MAT_NODI = ndgln%mat(CV_ILOC + (ele-1) * Mdims%cv_nloc)
                    CAP_DIFFUSION(Phase_with_Pc, MAT_NODI) = &
                        - T_ALL(Phase_with_Pc, CV_NODI) * OvRelax_param(CV_NODI)
                end do
            end do
        ENDIF

        if( is_flooding ) then
            allocate(DEN_ALL_DIVID( Mdims%nphase, Mdims%cv_nonods ))
            DEN_ALL_DIVID = DEN_ALL
            DEN_ALL_DIVID( 1, : ) = 1.0
        else
            DEN_ALL_DIVID => DEN_ALL
        endif
       
        ndotq = 0. ! ;         ndotqold = 0.
        ! variables for get_int_tden********************
        !      print *,'after allocate'
        ! END OF TEMP STUFF HERE
        ALLOCATE( R_PRES( Mdims%npres ), R_PHASE( Mdims%nphase ), MEAN_PORE_CV_PHASE( Mdims%nphase ) )
        IF ( Mdims%npres > 1 ) THEN
            ALLOCATE( A_GAMMA_PRES_ABS( Mdims%nphase, Mdims%nphase, Mdims%cv_nonods ) ); A_GAMMA_PRES_ABS = 0.
            ALLOCATE( GAMMA_PRES_ABS2( Mdims%nphase, Mdims%nphase, Mdims%cv_nonods ) ); GAMMA_PRES_ABS2 = 0.
            ALLOCATE( OPT_VEL_UPWIND_COEFS_NEW_CV( Mdims%nphase, Mdims%cv_nonods ) ); OPT_VEL_UPWIND_COEFS_NEW_CV = 0.
            ALLOCATE( SIGMA_INV_APPROX( Mdims%nphase, Mdims%cv_nonods ), SIGMA_INV_APPROX_NANO( Mdims%nphase, Mdims%cv_nonods ), N( Mdims%cv_nonods ) )
            SIGMA_INV_APPROX = 0; SIGMA_INV_APPROX_NANO = 0.; N = 0.
            ALLOCATE( RSUM_VEC( Mdims%nphase ) ); RSUM_VEC = 0.
            ALLOCATE( PIPE_ABS( Mdims%nphase, Mdims%nphase, Mdims%cv_nonods ) ); PIPE_ABS = 0.
        END IF
        psi(1)%ptr=>tracer
        psi(2)%ptr=>old_tracer
        FEM_IT=2
        if (.not. is_constant(density)) then
            psi(FEM_IT+1)%ptr=>density
            FEM_IT=FEM_IT+1
        end if
        if (.not. is_constant(old_density)) then
            psi(FEM_IT+1)%ptr=>old_density
            FEM_IT=FEM_IT+1
        end if
        if (present(saturation)) then
            if (.not. is_constant(saturation)) then
                psi(FEM_IT+1)%ptr=>saturation
                FEM_IT=FEM_IT+1
            end if
            if (.not. is_constant(old_saturation)) then
                psi(FEM_IT+1)%ptr=>old_saturation
                FEM_IT=FEM_IT+1
            end if
        end if
        do i=1,FEM_IT
            if (has_tensor_field(packed_state,&
                GetFEMName(psi(i)%ptr))) then
                fempsi(i)%ptr=>extract_tensor_field(packed_state,&
                    GetFEMName(psi(i)%ptr))
            else
                allocate(fempsi(i)%ptr)
                call allocate(fempsi(i)%ptr,psi(i)%ptr%mesh,"FEMPSI"//trim(psi(i)%ptr%name),&
                    dim=psi(i)%ptr%dim)
            end if
        end do
        psi_int(1)%ptr=>extract_vector_field(packed_state,"CVIntegral")
        psi_ave(1)%ptr=>extract_vector_field(packed_state,"CVBarycentre")
        call PROJ_CV_TO_FEM(packed_state, &
            FEMPSI(1:FEM_IT),PSI(1:FEM_IT), &
            Mdims, CV_GIdims, CV_funs, Mspars, ndgln, &
            IGETCT, X_ALL, MASS_ELE, MASS_MN_PRES, &
            tracer,PSI_AVE, PSI_INT)
        XC_CV_ALL=0.0
        !sprint_to_do!use a pointer?
        XC_CV_ALL(1:Mdims%ndim,:)=psi_ave(1)%ptr%val
        MASS_CV=psi_int(1)%ptr%val(1,:)
        FEMT_ALL(:,:)=FEMPSI(1)%ptr%val(1,:,:)
        FEMTOLD_ALL(:,:)=FEMPSI(2)%ptr%val(1,:,:)
        FEM_IT=3
        if (.not. is_constant(density)) then
            FEMDEN_ALL=psi(FEM_IT)%ptr%val(1,:,:)
            tfield => extract_tensor_field( packed_state, "PackedFEDensity" )
            tfield%val = psi(FEM_IT)%ptr%val
            FEM_IT=FEM_IT+1
        else
            FEMDEN_ALL=density%val(1,:,:)
        end if
        if (.not. is_constant(old_density)) then
            FEMDENOLD_ALL=psi(FEM_IT)%ptr%val(1,:,:)
            FEM_IT=FEM_IT+1
        else
            FEMDENOLD_ALL=old_density%val(1,:,:)
        end if
        IF ( present(saturation) ) then
            if (.not. is_constant(saturation)) then
                FEMT2_ALL=psi(FEM_IT)%ptr%val(1,:,:)
                FEM_IT=FEM_IT+1
            else
                FEMT2_ALL=saturation%val(1,:,:)
            end if
            if (.not. is_constant(old_saturation)) then
                FEMt2OLD_ALL=psi(FEM_IT)%ptr%val(1,:,:)
                FEM_IT=FEM_IT+1
            else
                FEMt2OLD_ALL=old_saturation%val(1,:,:)
            end if
        end IF
        do i=1,FEM_IT-1
            if (fempsi(i)%ptr%name(1:6)=='FEMPSI') then
                call deallocate(fempsi(i)%ptr)
                deallocate(fempsi(i)%ptr)
            end if
        end do
        IF (PRESENT(MASS_ELE_TRANSP)) &
            MASS_ELE_TRANSP = MASS_ELE
        NORMALISE = .FALSE.
        IF ( NORMALISE ) THEN
            ! Make sure the FEM representation sums to unity so we don't get surprising results...
            DO CV_INOD = 1, Mdims%cv_nonods
                RSUM = SUM( FEMT_ALL( :, CV_INOD ) )
                FEMT_ALL( :, CV_INOD ) = FEMT_ALL( :, CV_INOD ) / RSUM
            END DO
        END IF
        ! Calculate MEAN_PORE_CV
        MEAN_PORE_CV = 0.0 ; SUM_CV = 0.0
        DO ELE = 1, Mdims%totele
            DO CV_ILOC = 1, Mdims%cv_nloc
                CV_INOD = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
                SUM_CV( CV_INOD ) = SUM_CV( CV_INOD ) + MASS_ELE( ELE )
                MEAN_PORE_CV( :, CV_INOD ) = MEAN_PORE_CV( :, CV_INOD ) + &
                    MASS_ELE( ELE ) * VOLFRA_PORE( :, ELE )
            END DO
        END DO
        DO IPRES = 1, Mdims%npres
            MEAN_PORE_CV(IPRES,:) = MEAN_PORE_CV(IPRES,:) / SUM_CV
        END DO

        ! Scale effectively the time step size used within the pipes...
        ! dt_pipe_factor is the factor by which to reduce the pipe eqns time step size e.g. 10^{-3}
        DO IPRES = 2, Mdims%npres
            MEAN_PORE_CV(IPRES,:) = MEAN_PORE_CV(IPRES,:) / dt_pipe_factor
        END DO
        ewrite(3,*) 'MEAN_PORE_CV MIN/MAX:', MINVAL( MEAN_PORE_CV ), MAXVAL( MEAN_PORE_CV )
        MeanPoreCV=>extract_vector_field(packed_state,"MeanPoreCV")
        MeanPoreCV%val=MEAN_PORE_CV

        IANISOLIM = 0
        IF ( CV_DISOPT >= 5 ) IANISOLIM = 1
        ALLOCATE( TUPWIND_MAT_ALL( Mdims%nphase, Mspars%small_acv%ncol ), TOLDUPWIND_MAT_ALL( Mdims%nphase, Mspars%small_acv%ncol ), &
            DENUPWIND_MAT_ALL( Mdims%nphase, Mspars%small_acv%ncol ), DENOLDUPWIND_MAT_ALL( Mdims%nphase, Mspars%small_acv%ncol ) )
        ALLOCATE( T2UPWIND_MAT_ALL( Mdims%nphase*i_use_volume_frac_t2, Mspars%small_acv%ncol* i_use_volume_frac_t2), T2OLDUPWIND_MAT_ALL( Mdims%nphase*i_use_volume_frac_t2, Mspars%small_acv%ncol*i_use_volume_frac_t2 ) )
        IF ( IANISOLIM == 0 ) THEN
            ! Isotropic limiting - calculate far field upwind maticies...
            CALL ISOTROPIC_LIMITER_ALL( &
                ! FOR SUB SURRO_CV_MINMAX:
                T_ALL, TOLD_ALL, T2_ALL, T2OLD_ALL, DEN_ALL, DENOLD_ALL, i_use_volume_frac_t2, Mdims%nphase, Mdims%cv_nonods, Mspars%small_acv%ncol, Mspars%small_acv%mid, Mspars%small_acv%fin, Mspars%small_acv%col, &
                Mdims%stotel, Mdims%cv_snloc, ndgln%suf_cv, SUF_T_BC_ALL, SUF_T2_BC_ALL, SUF_D_BC_ALL, WIC_T_BC_ALL, WIC_T2_BC_ALL, WIC_D_BC_ALL, &
                MASS_CV, &
                ! FOR SUB CALC_LIMIT_MATRIX_MAX_MIN:
                TOLDUPWIND_MAT_ALL, DENOLDUPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL, &
                TUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, T2UPWIND_MAT_ALL )
        ELSE ! endof IF ( IANISOLIM == 0 ) THEN
            use_reflect = have_option("/numerical_methods/use_reflect_method")
            CALL CALC_ANISOTROP_LIM( &
                ! Caculate the upwind values stored in matrix form...
                T_ALL,TOLD_ALL,DEN_ALL,DENOLD_ALL,T2_ALL,T2OLD_ALL, &
                FEMT_ALL,FEMTOLD_ALL,FEMDEN_ALL,FEMDENOLD_ALL,FEMT2_ALL,FEMT2OLD_ALL, (Mdims%cv_nonods.NE.Mdims%x_nonods), &
                TUPWIND_MAT_ALL, TOLDUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, DENOLDUPWIND_MAT_ALL, &
                T2UPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL, &
                ! Store the upwind element for interpolation and its weights for
                ! faster results...
                i_use_volume_frac_t2,Mdims%nphase,Mdims%cv_nonods,Mdims%cv_nloc,Mdims%totele,ndgln%cv, &
                Mspars%small_acv%fin,Mspars%small_acv%mid,Mspars%small_acv%col,Mspars%small_acv%ncol, &
                ndgln%x,Mdims%x_nonods,Mdims%ndim, &
                X_ALL, XC_CV_ALL, use_reflect)
        END IF ! endof IF ( IANISOLIM == 0 ) THEN ELSE
        ALLOCATE( FACE_ELE( CV_GIdims%nface, Mdims%totele ) ) ; FACE_ELE = 0
        CALL CALC_FACE_ELE( FACE_ELE, Mdims%totele, Mdims%stotel, CV_GIdims%nface, &
            Mspars%ELE%fin, Mspars%ELE%col, Mdims%cv_nloc, Mdims%cv_snloc, Mdims%cv_nonods, ndgln%cv, ndgln%suf_cv, &
            CV_funs%cv_sloclist, Mdims%x_nloc, ndgln%x )
        IF ( GOT_DIFFUS ) THEN
            CALL DG_DERIVS_ALL( FEMT_ALL, FEMTOLD_ALL, &
                DTX_ELE_ALL, DTOLDX_ELE_ALL, &
                Mdims%ndim, Mdims%nphase, Mdims%totele, ndgln%cv, &
                ndgln%x, Mdims%x_nloc, ndgln%x,&
                CV_GIdims%cv_ngi, Mdims%cv_nloc, CV_funs%CVWEIGHT, &
                CV_funs%CVFEN, CV_funs%CVFENLX_ALL(1,:,:), CV_funs%CVFENLX_ALL(2,:,:), CV_funs%CVFENLX_ALL(3,:,:), &
                CV_funs%CVFEN, CV_funs%CVFENLX_ALL(1,:,:), CV_funs%CVFENLX_ALL(2,:,:), CV_funs%CVFENLX_ALL(3,:,:), &
                Mdims%x_nonods, X_ALL(1,:),X_ALL(2,:),X_ALL(3,:), &
                CV_GIdims%nface, FACE_ELE, CV_funs%cv_sloclist, CV_funs%cv_sloclist, Mdims%cv_snloc, Mdims%cv_snloc, WIC_T_BC_ALL, SUF_T_BC_ALL, &
                CV_GIdims%sbcvngi, CV_funs%sbcvfen, CV_funs%sbcvfeweigh, &
                CV_funs%sbcvfen, CV_funs%sbcvfenslx, CV_funs%sbcvfensly )
        END IF
        !     =============== DEFINE THETA FOR TIME-STEPPING ===================
        ! Define the type of time integration:
        ! Timopt is 0 if CV_DISOPT is even (theta specified);
        ! Timopt is 1 if CV_DISOPT is odd (non-linear theta scheme)
        TIMOPT = MOD( CV_DISOPT, 2 )
        FTHETA(:) = CV_THETA
        IF ( GETCT ) THEN ! Obtain the CV discretised Mmat%CT eqns plus RHS
            call zero(Mmat%CT_RHS)
            Mmat%CT = 0.0
            if ( got_free_surf .and. .not.symmetric_P ) MASS_SUF=0.0
            if ( Mmat%CV_pressure ) MASS_SUF=0.0
        END IF
        IF ( GETCV_DISC ) THEN ! Obtain the CV discretised advection/diffusion eqns
            call zero(Mmat%CV_RHS)
            IF ( GETMAT ) THEN
                call zero(Mmat%petsc_ACV)
            END IF
        END IF
        GET_GTHETA = .FALSE.
        IF ( IGOT_THETA_FLUX == 1 ) THEN
            IF ( GET_THETA_FLUX ) THEN
                THETA_FLUX = 0.0
                ONE_M_THETA_FLUX = 0.0
                if(integrate_other_side) then
                    THETA_FLUX_J = 0.0
                    ONE_M_THETA_FLUX_J = 0.0
                endif
                !            IF ( IGOT_T2 == 1 ) THEN
                IF ( GOT_T2 ) THEN
                    GET_GTHETA = .TRUE.
                    THETA_GDIFF = 0.0
                END IF
            END IF
        END IF
        GLOBAL_FACE = 0

        if ( is_porous_media .and. present(calculate_mass_delta) .and. present(outfluxes) ) then
            !Initialise mass conservation check; calculation of porevolume
            call mass_conservation_check_and_outfluxes(calculate_mass_delta, outfluxes, 1)
        endif


        !Allocate derivatives of the shape functions
        call allocate_multi_dev_shape_funs(CV_funs%scvfenlx_all, CV_funs%sufenlx_all, SdevFuns)
        !###########################################
        Loop_Elements: DO ELE = 1, Mdims%totele
            if (IsParallel()) then
                if (.not. assemble_ele(tracer,ele)) then
                    skip=.true.
                    neighbours=>ele_neigh(tracer,ele)
                    do nb=1,size(neighbours)
                        if (neighbours(nb)<=0) cycle
                        if (assemble_ele(tracer,neighbours(nb))) then
                            skip=.false.
                            exit
                        end if
                    end do
                    if (skip) cycle
                end if
            end if


            ! Generate some local F variables ***************
            DO CV_KLOC = 1, Mdims%cv_nloc
                CV_NODK = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_KLOC )
                ! loc_f
                IPT=1
                CALL PACK_LOC( LOC_F(:, CV_KLOC), T_ALL( :, CV_NODK ),    Mdims%nphase, IPT, IGOT_T_PACK(:,1) )
                CALL PACK_LOC( LOC_F(:, CV_KLOC), TOLD_ALL( :, CV_NODK ), Mdims%nphase, IPT, IGOT_T_PACK(:,2) )
                CALL PACK_LOC( LOC_F(:, CV_KLOC), DEN_ALL( :, CV_NODK ),    Mdims%nphase, IPT, IGOT_T_PACK(:,3) )
                CALL PACK_LOC( LOC_F(:, CV_KLOC), DENOLD_ALL( :, CV_NODK ), Mdims%nphase, IPT, IGOT_T_PACK(:,4) )
                IF(use_volume_frac_T2) THEN
                    CALL PACK_LOC( LOC_F(:, CV_KLOC), T2_ALL( :, CV_NODK ),    Mdims%nphase, IPT, IGOT_T_PACK(:,5) )
                    CALL PACK_LOC( LOC_F(:, CV_KLOC), T2OLD_ALL( :, CV_NODK ), Mdims%nphase, IPT, IGOT_T_PACK(:,6) )
                ENDIF
                ! for FEM variables...
                IPT=1
                CALL PACK_LOC( LOC_FEMF(:, CV_KLOC), FEMT_ALL( :, CV_NODK ),    Mdims%nphase, IPT, IGOT_T_PACK(:,1) )
                CALL PACK_LOC( LOC_FEMF(:, CV_KLOC), FEMTOLD_ALL( :, CV_NODK ), Mdims%nphase, IPT, IGOT_T_PACK(:,2) )
                CALL PACK_LOC( LOC_FEMF(:, CV_KLOC), FEMDEN_ALL( :, CV_NODK ),    Mdims%nphase, IPT, IGOT_T_PACK(:,3) )
                CALL PACK_LOC( LOC_FEMF(:, CV_KLOC), FEMDENOLD_ALL( :, CV_NODK ), Mdims%nphase, IPT, IGOT_T_PACK(:,4) )
                IF(use_volume_frac_T2) THEN
                    CALL PACK_LOC( LOC_FEMF(:, CV_KLOC), FEMT2_ALL( :, CV_NODK ),    Mdims%nphase, IPT, IGOT_T_PACK(:,5) )
                    CALL PACK_LOC( LOC_FEMF(:, CV_KLOC), FEMT2OLD_ALL( :, CV_NODK ), Mdims%nphase, IPT, IGOT_T_PACK(:,6) )
                ENDIF
                ! loc_femt:
                LOC_FEMT(:, CV_KLOC) = FEMT_ALL(:, CV_NODK)
                LOC_FEMTOLD(:, CV_KLOC) = FEMTOLD_ALL(:, CV_NODK)
                IF ( use_volume_frac_t2 ) THEN
                    LOC_FEMT2(:, CV_KLOC) = FEMT2_ALL(:, CV_NODK)
                    LOC_FEMT2OLD(:, CV_KLOC) = FEMT2OLD_ALL(:, CV_NODK)
                END IF
            END DO
            DO U_KLOC = 1, Mdims%u_nloc
                U_NODK = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_KLOC )
                ! LOC_UF
                IPT=1
                IFI=1
                DO ILOOP=1,3
                    DO IPHASE=1,Mdims%nphase
                        IF(IGOT_T_PACK(IPHASE,IFI)) THEN ! T:
                            LOC_UF(:, IPT, U_KLOC) =   U_ALL( :, IPHASE, U_NODK )
                            IPT=IPT+1
                        END IF
                    END DO
                    IFI=IFI+1
                    DO IPHASE=1,Mdims%nphase
                        IF(IGOT_T_PACK(IPHASE,IFI)) THEN ! Told:
                            if(correct_method_petrov_method) then
                                LOC_UF(:, IPT, U_KLOC) =   NUOLD_ALL( :, IPHASE, U_NODK )
                            else
                                LOC_UF(:, IPT, U_KLOC) =   U_ALL( :, IPHASE, U_NODK )
                            endif
                            IPT=IPT+1
                        END IF
                    END DO
                    IFI=IFI+1
                END DO
                ! LOC_U, LOC_NU:
                LOC_U( :, :, U_KLOC)=U_ALL( :, :, U_NODK)
                LOC_NU( :, :, U_KLOC)=NU_ALL( :, :, U_NODK)
                LOC_NUOLD( :, :, U_KLOC)=NUOLD_ALL( :, :, U_NODK)
            END DO
            ! Generate some local F variables ***************...
            !
            !
            Loop_CV_ILOC: DO CV_ILOC = 1, Mdims%cv_nloc ! Loop over the nodes of the element




                ! Global node number of the local node
                CV_NODI = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
                X_NODI = ndgln%x( ( ELE - 1 ) * Mdims%x_nloc  + CV_ILOC )
                MAT_NODI = ndgln%mat( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
                IMID = Mspars%small_acv%mid(CV_NODI)
                ! Generate some local F variables ***************
                F_CV_NODI(:)= LOC_F(:, CV_ILOC)
                ! Generate some local F variables ***************
                ! Loop over quadrature (gauss) points in ELE neighbouring ILOC
                Loop_GCOUNT: DO GCOUNT = CV_funs%findgpts( CV_ILOC ), CV_funs%findgpts( CV_ILOC + 1 ) - 1
                    ! CV_funs%colgpts stores the local Gauss-point number in the ELE
                    GI = CV_funs%colgpts( GCOUNT )
                    ! Get the neighbouring node for node ILOC and Gauss point GI
                    CV_JLOC = CV_funs%cv_neiloc( CV_ILOC, GI )
                    ELE2 = 0
                    SELE = 0
                    CV_SILOC=0
                    INTEGRAT_AT_GI = .TRUE.
                    U_OTHER_LOC=0
                    CV_OTHER_LOC=0
                    Conditional_CheckingNeighbourhood: IF ( CV_JLOC == -1 ) THEN
                        ! We are on the boundary or next to another element.  Determine CV_OTHER_LOC
                        ! CV_funs%cvfem_on_face(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
                        ! Look for these nodes on the other elements.
                        CALL FIND_OTHER_SIDE( CV_OTHER_LOC, Mdims%cv_nloc, U_OTHER_LOC, Mdims%u_nloc, &
                            MAT_OTHER_LOC, INTEGRAT_AT_GI, &
                            Mdims%x_nloc, Mdims%xu_nloc, ndgln%x, ndgln%xu, &
                            Mdims%cv_snloc, CV_funs%cvfem_on_face( :, GI ), X_SHARE, ELE, ELE2,  &
                            Mspars%ELE%fin, Mspars%ELE%col, DISTCONTINUOUS_METHOD )

                        IF ( INTEGRAT_AT_GI ) THEN
                            CV_JLOC = CV_OTHER_LOC( CV_ILOC )
                            SELE = 0
                            ELE3=0
                            IF ( CV_JLOC == 0 ) THEN ! We are on the boundary of the domain or subdomain
                                CV_JLOC = CV_ILOC
                                ! Calculate SELE, CV_SILOC, U_SLOC2LOC, CV_SLOC2LOC
                                CALL CALC_SELE( ELE, ELE3, SELE, CV_SILOC, CV_ILOC, U_SLOC2LOC, CV_SLOC2LOC, &
                                    FACE_ELE, GI, CV_funs, Mdims, CV_GIdims,&
                                    ndgln%cv, ndgln%u, ndgln%suf_cv, ndgln%suf_u )
                            END IF
                            INTEGRAT_AT_GI = .NOT.( (ELE==ELE2) .AND. (SELE==0) )
                        END IF
                    END IF Conditional_CheckingNeighbourhood
                    !Avoid integrating across the middle of a CV on the boundaries of elements
                    !For disconditnuous pressure this is always true
                    Conditional_integration: IF ( INTEGRAT_AT_GI ) THEN
                        between_elements = (ELE2 /= 0) .AND. (ELE2 /= ELE)
                        on_domain_boundary = ( SELE /= 0 )
                        !Decide whether to use high order advection or not, locally.
                        !If the saturation is equal to the minimum, then no need to high order
                        if (is_porous_media .and. local_upwinding) then
                            use_porous_limiter = .true.
                            do iphase = 1, Mdims%nphase - 1
                                use_porous_limiter = use_porous_limiter &
                                    .and. abs(T_ALL(iphase, cv_inod) - Imble_frac(iphase, ELE)) > 1e-4
                            end do
                            use_porous_limiter = use_porous_limiter .or. on_domain_boundary
                        end if
                        !Check as well if the CV is neigbouring another CV with a different permeability and if any of those is anisotropic
                        anisotropic_and_frontier = anisotropic_perm .and. between_elements
                        !Check if the permeabilities are different, this may fail in the unlikely case where the contraction of
                        !the difference also sum 0 despite being different permeabilities
                        if (anisotropic_and_frontier) anisotropic_and_frontier = sum(perm%val(:,:,ele) - perm%val(:,:,ele2)) < &
                            10.**(exponent(perm%val(1,1,ele))-9)
                        IF ( between_elements ) THEN
                            CV_NODJ = ndgln%cv( ( ELE2 - 1 ) * Mdims%cv_nloc + CV_JLOC )
                            X_NODJ = ndgln%x( ( ELE2 - 1 ) * Mdims%cv_nloc + CV_JLOC )
                            MAT_NODJ = ndgln%mat( ( ELE2 - 1 ) * Mdims%cv_nloc + CV_JLOC )
                        ELSE
                            CV_NODJ = ndgln%cv( ( ELE - 1 )  * Mdims%cv_nloc + CV_JLOC )
                            X_NODJ = ndgln%x( ( ELE - 1 )  * Mdims%cv_nloc + CV_JLOC )
                            MAT_NODJ = ndgln%mat( ( ELE - 1 )  * Mdims%cv_nloc + CV_JLOC )
                        END IF
                        if((.not.integrate_other_side).or.(CV_NODJ >= CV_NODI)) then
                            ! this is for DG and boundaries of the domain
                            IF( between_elements ) THEN ! this is for DG
                                ! Calculate U_SLOC2LOC, CV_SLOC2LOC:
                                CV_SKLOC=0
                                DO CV_KLOC=1,Mdims%cv_nloc
                                    CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                                    IF(CV_KLOC2.NE.0) THEN
                                        CV_SKLOC=CV_SKLOC+1
                                        CV_SLOC2LOC(CV_SKLOC)=CV_KLOC
                                        SHAPE_CV_SNL(CV_SKLOC) = CV_funs%scvfen(CV_KLOC,GI)
                                    ENDIF
                                END DO
                                U_SKLOC=0
                                DO U_KLOC=1,Mdims%u_nloc
                                    U_KLOC2 = U_OTHER_LOC( U_KLOC )
                                    IF(U_KLOC2.NE.0) THEN
                                        U_SKLOC=U_SKLOC+1
                                        U_SLOC2LOC(U_SKLOC)=U_KLOC
                                    ENDIF
                                END DO
                            ENDIF ! ENDOF IF( between_elements ) THEN
                            integrate_other_side_and_not_boundary = integrate_other_side.and.(SELE.LE.0)
                            GLOBAL_FACE = GLOBAL_FACE + 1
                            JMID = Mspars%small_acv%mid(CV_NODJ)
                            ! Calculate the control volume normals at the Gauss pts.
                            CALL SCVDETNX_new( ELE, GI, SdevFuns%DETWEI, CVNORMX_ALL,XC_CV_ALL( 1:Mdims%ndim, CV_NODI ), X_NODI)
                            !Obtain the list of neighbouring nodes
                            IF( GETCT ) call get_neigbouring_lists(JCOUNT_KLOC, ICOUNT_KLOC, JCOUNT_KLOC2 ,ICOUNT_KLOC2,&
                                                            C_JCOUNT_KLOC, C_ICOUNT_KLOC, C_JCOUNT_KLOC2, C_ICOUNT_KLOC2 )

                            ! Compute the distance HDC between the nodes either side of the CV face
                            ! (this is needed to compute the local courant number and the non-linear theta)
                            IF ( on_domain_boundary) THEN
                                HDC = SQRT( SUM( (XC_CV_ALL(1:Mdims%ndim,CV_NODI)-X_ALL(1:Mdims%ndim,X_NODI))**2) )
                            ELSE
                                HDC = SQRT( SUM( (XC_CV_ALL(1:Mdims%ndim,CV_NODI)-XC_CV_ALL(1:Mdims%ndim,CV_NODJ))**2) )
                            END IF
                            DO COUNT = Mspars%small_acv%fin( CV_NODI ), Mspars%small_acv%fin( CV_NODI + 1 ) - 1
                                IF ( Mspars%small_acv%col( COUNT ) == CV_NODJ ) THEN
                                    count_out = COUNT
                                    EXIT
                                END IF
                            END DO
                            DO COUNT = Mspars%small_acv%fin( CV_NODJ ), Mspars%small_acv%fin( CV_NODJ + 1 ) - 1
                                IF ( Mspars%small_acv%col( COUNT ) == CV_NODI ) THEN
                                    count_in = COUNT
                                    EXIT
                                END IF
                            END DO
                            ! Generate some local F variables ***************
                            IPT=1
                            CALL PACK_LOC( F_CV_NODJ(:), T_ALL( :, CV_NODJ ),    Mdims%nphase, IPT, IGOT_T_PACK(:,1) )
                            CALL PACK_LOC( F_CV_NODJ(:), TOLD_ALL( :, CV_NODJ ), Mdims%nphase, IPT, IGOT_T_PACK(:,2) )
                            CALL PACK_LOC( F_CV_NODJ(:), DEN_ALL( :, CV_NODJ ),    Mdims%nphase, IPT, IGOT_T_PACK(:,3) )
                            CALL PACK_LOC( F_CV_NODJ(:), DENOLD_ALL( :, CV_NODJ ), Mdims%nphase, IPT, IGOT_T_PACK(:,4) )
                            IF(use_volume_frac_T2) THEN
                                CALL PACK_LOC( F_CV_NODJ(:), T2_ALL( :, CV_NODJ ),    Mdims%nphase, IPT, IGOT_T_PACK(:,5) )
                                CALL PACK_LOC( F_CV_NODJ(:), T2OLD_ALL( :, CV_NODJ ), Mdims%nphase, IPT, IGOT_T_PACK(:,6) )
                            ENDIF
                            ! Generate some local F variables ***************
                            ! local surface information***********
                            IF( between_elements .or. on_domain_boundary ) THEN
                                DO CV_SKLOC = 1, Mdims%cv_snloc
                                    CV_KLOC = CV_SLOC2LOC( CV_SKLOC )
                                    SLOC_F(:, CV_SKLOC) = LOC_F(:, CV_KLOC)
                                    SLOC_FEMF(:, CV_SKLOC) = LOC_FEMF(:, CV_KLOC)
                                    !                   IF(ELE2>0) THEN
                                    IF(between_elements) THEN
                                        CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                                        CV_NODK2 = ndgln%cv( ( ELE2 - 1 ) * Mdims%cv_nloc + CV_KLOC2 )
                                        IPT=1
                                        CALL PACK_LOC( SLOC2_F(:, CV_SKLOC), T_ALL( :, CV_NODK2 ),    Mdims%nphase, IPT, IGOT_T_PACK(:,1) )
                                        CALL PACK_LOC( SLOC2_F(:, CV_SKLOC), TOLD_ALL( :, CV_NODK2 ), Mdims%nphase, IPT, IGOT_T_PACK(:,2) )
                                        CALL PACK_LOC( SLOC2_F(:, CV_SKLOC), DEN_ALL( :, CV_NODK2 ),    Mdims%nphase, IPT, IGOT_T_PACK(:,3) )
                                        CALL PACK_LOC( SLOC2_F(:, CV_SKLOC), DENOLD_ALL( :, CV_NODK2 ), Mdims%nphase, IPT, IGOT_T_PACK(:,4) )
                                        IF(use_volume_frac_T2) THEN
                                            CALL PACK_LOC( SLOC2_F(:, CV_SKLOC), T2_ALL( :, CV_NODK2 ),    Mdims%nphase, IPT, IGOT_T_PACK(:,5) )
                                            CALL PACK_LOC( SLOC2_F(:, CV_SKLOC), T2OLD_ALL( :, CV_NODK2 ), Mdims%nphase, IPT, IGOT_T_PACK(:,6) )
                                        ENDIF
                                        ! femf:
                                        IPT=1
                                        CALL PACK_LOC( SLOC2_FEMF(:, CV_SKLOC), FEMT_ALL( :, CV_NODK2 ),    Mdims%nphase, IPT, IGOT_T_PACK(:,1) )
                                        CALL PACK_LOC( SLOC2_FEMF(:, CV_SKLOC), FEMTOLD_ALL( :, CV_NODK2 ), Mdims%nphase, IPT, IGOT_T_PACK(:,2) )
                                        CALL PACK_LOC( SLOC2_FEMF(:, CV_SKLOC), FEMDEN_ALL( :, CV_NODK2 ),    Mdims%nphase, IPT, IGOT_T_PACK(:,3) )
                                        CALL PACK_LOC( SLOC2_FEMF(:, CV_SKLOC), FEMDENOLD_ALL( :, CV_NODK2 ), Mdims%nphase, IPT, IGOT_T_PACK(:,4) )
                                        IF(use_volume_frac_T2) THEN
                                            CALL PACK_LOC( SLOC2_FEMF(:, CV_SKLOC), FEMT2_ALL( :, CV_NODK2 ),    Mdims%nphase, IPT, IGOT_T_PACK(:,5) )
                                            CALL PACK_LOC( SLOC2_FEMF(:, CV_SKLOC), FEMT2OLD_ALL( :, CV_NODK2 ), Mdims%nphase, IPT, IGOT_T_PACK(:,6) )
                                        ENDIF
                                    ELSE
                                        SLOC2_F(:, CV_SKLOC)    = SLOC_F(:, CV_SKLOC)
                                        SLOC2_FEMF(:, CV_SKLOC) = SLOC_FEMF(:, CV_SKLOC)
                                    ENDIF
                                END DO
                            ENDIF ! ENDOF IF( between_elements .or. on_domain_boundary ) THEN ...
                            IF( on_domain_boundary ) THEN
                                ! bcs:
                                ! Make allowances for no matrix stencil operating from outside the boundary.
                                !             BCZERO=1.0-INCOME
                                ! What type of b.c's -integer
                                IPT=1
                                CALL I_PACK_LOC( SELE_LOC_WIC_F_BC( : ),WIC_T_BC_ALL( : , : , SELE ),&
                                    Mdims%nphase, IPT, IGOT_T_PACK( :,1) )
                                CALL I_PACK_LOC( SELE_LOC_WIC_F_BC( : ),WIC_T_BC_ALL( : , :, SELE ),&
                                    Mdims%nphase, IPT, IGOT_T_PACK( :,2) )
                                CALL I_PACK_LOC( SELE_LOC_WIC_F_BC( : ),WIC_D_BC_ALL( :,:, SELE ),&
                                    Mdims%nphase, IPT, IGOT_T_PACK( :,3) )
                                CALL I_PACK_LOC( SELE_LOC_WIC_F_BC( : ),WIC_D_BC_ALL( :,:, SELE ),&
                                    Mdims%nphase, IPT, IGOT_T_PACK( :,4) )
                                IF(use_volume_frac_T2) THEN
                                    CALL I_PACK_LOC( SELE_LOC_WIC_F_BC( : ),WIC_T2_BC_ALL( : , :, SELE ),&
                                        Mdims%nphase, IPT, IGOT_T_PACK( :,5) )
                                    CALL I_PACK_LOC( SELE_LOC_WIC_F_BC( : ),WIC_T2_BC_ALL( : , : , SELE ),&
                                        Mdims%nphase, IPT, IGOT_T_PACK( :,6) )
                                ENDIF
                                ! The b.c's values:
                                DO CV_SKLOC=1,Mdims%cv_snloc
                                    IPT=1
                                    CALL PACK_LOC( SLOC_SUF_F_BC( :, CV_SKLOC ), SUF_T_BC_ALL( 1, :, CV_SKLOC + Mdims%cv_snloc*( SELE- 1) ),    Mdims%nphase, IPT, IGOT_T_PACK(:,1) )
                                    CALL PACK_LOC( SLOC_SUF_F_BC( :, CV_SKLOC ), SUF_T_BC_ALL( 1, :, CV_SKLOC+ Mdims%cv_snloc*( SELE- 1 ) ),    Mdims%nphase, IPT, IGOT_T_PACK(:,2) )
                                    CALL PACK_LOC( SLOC_SUF_F_BC( :, CV_SKLOC ), SUF_D_BC_ALL( 1, :, CV_SKLOC+ Mdims%cv_snloc*( SELE- 1) ),    Mdims%nphase, IPT, IGOT_T_PACK(:,3) )
                                    CALL PACK_LOC( SLOC_SUF_F_BC( :, CV_SKLOC ), SUF_D_BC_ALL( 1, :, CV_SKLOC+ Mdims%cv_snloc*( SELE- 1) ),    Mdims%nphase, IPT, IGOT_T_PACK(:,4) )
                                    IF(use_volume_frac_T2) THEN
                                        CALL PACK_LOC( SLOC_SUF_F_BC( :, CV_SKLOC ), SUF_T2_BC_ALL( 1, :, CV_SKLOC + Mdims%cv_snloc*( SELE- 1) ),    Mdims%nphase, IPT, IGOT_T_PACK(:,5) )
                                        CALL PACK_LOC( SLOC_SUF_F_BC( :, CV_SKLOC ), SUF_T2_BC_ALL( 1, :, CV_SKLOC + Mdims%cv_snloc*( SELE- 1) ),    Mdims%nphase, IPT, IGOT_T_PACK(:,6) )
                                    ENDIF
                                END DO
                            ENDIF ! IF( on_domain_boundary ) THEN
                            ! local surface information***********
                            ! limiting VALUES*************:
                            IPT_IN =1
                            IPT_OUT=1
                            CALL PACK_LOC( FUPWIND_IN( : ),  TUPWIND_MAT_ALL( :, COUNT_IN),    Mdims%nphase, IPT_IN, IGOT_T_PACK(:,1) )
                            CALL PACK_LOC( FUPWIND_OUT( : ), TUPWIND_MAT_ALL( :, COUNT_OUT),    Mdims%nphase, IPT_OUT, IGOT_T_PACK(:,1) )
                            CALL PACK_LOC( FUPWIND_IN( : ),  TOLDUPWIND_MAT_ALL( :, COUNT_IN),    Mdims%nphase, IPT_IN, IGOT_T_PACK(:,2) )
                            CALL PACK_LOC( FUPWIND_OUT( : ), TOLDUPWIND_MAT_ALL( :, COUNT_OUT),    Mdims%nphase, IPT_OUT, IGOT_T_PACK(:,2) )
                            CALL PACK_LOC( FUPWIND_IN( : ),  DENUPWIND_MAT_ALL( :, COUNT_IN),    Mdims%nphase, IPT_IN, IGOT_T_PACK(:,3) )
                            CALL PACK_LOC( FUPWIND_OUT( : ), DENUPWIND_MAT_ALL( :, COUNT_OUT),    Mdims%nphase, IPT_OUT, IGOT_T_PACK(:,3) )
                            CALL PACK_LOC( FUPWIND_IN( : ),  DENOLDUPWIND_MAT_ALL( :, COUNT_IN),    Mdims%nphase, IPT_IN, IGOT_T_PACK(:,4) )
                            CALL PACK_LOC( FUPWIND_OUT( : ), DENOLDUPWIND_MAT_ALL( :, COUNT_OUT),    Mdims%nphase, IPT_OUT, IGOT_T_PACK(:,4) )
                            IF(use_volume_frac_T2) THEN
                                CALL PACK_LOC( FUPWIND_IN( : ),  T2UPWIND_MAT_ALL( :, COUNT_IN),    Mdims%nphase, IPT_IN, IGOT_T_PACK(:,5) )
                                CALL PACK_LOC( FUPWIND_OUT( : ), T2UPWIND_MAT_ALL( :, COUNT_OUT),    Mdims%nphase, IPT_OUT, IGOT_T_PACK(:,5) )
                                CALL PACK_LOC( FUPWIND_IN( : ),  T2OLDUPWIND_MAT_ALL( :, COUNT_IN),    Mdims%nphase, IPT_IN, IGOT_T_PACK(:,6) )
                                CALL PACK_LOC( FUPWIND_OUT( : ), T2OLDUPWIND_MAT_ALL( :, COUNT_OUT),    Mdims%nphase, IPT_OUT, IGOT_T_PACK(:,6) )
                            ENDIF
                            ! limiting VALUES*************:
                            !     endif ! endof if(.false.) then
                            !
                            ! LOC2_U, LOC2_NU for GET_INT_VEL_NEW
                            !       IF (ELE2/=0) THEN
                            IF (between_elements) THEN
                                LOC2_U = 0.
                                LOC2_NU = 0.
                                LOC2_NUOLD = 0.
                                DO U_SKLOC = 1, Mdims%u_snloc
                                    U_KLOC = U_SLOC2LOC(U_SKLOC)
                                    U_KLOC2 = U_OTHER_LOC( U_KLOC )
                                    U_NODK2 = ndgln%u((ELE2-1)*Mdims%u_nloc+U_KLOC2)
                                    LOC2_U(:, :, U_KLOC) = U_ALL(:, :, U_NODK2)
                                    LOC2_NU(:, :, U_KLOC) = NU_ALL(:, :, U_NODK2)
                                    LOC2_NUOLD(:, :, U_KLOC) = NUOLD_ALL(:, :, U_NODK2)
                                END DO
                                DO CV_SKLOC=1,Mdims%cv_snloc
                                    CV_KLOC=CV_SLOC2LOC( CV_SKLOC )
                                    CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                                    CV_KNOD2 = ndgln%cv((ELE2-1)*Mdims%cv_nloc+CV_KLOC2)
                                    LOC2_FEMT(:, CV_KLOC) = FEMT_ALL(:, CV_KNOD2)
                                    LOC2_FEMTOLD(:, CV_KLOC) = FEMTOLD_ALL(:, CV_KNOD2)
                                    IF (use_volume_frac_T2) THEN
                                        LOC2_FEMT2(:, CV_KLOC) = FEMT2_ALL(:, CV_KNOD2)
                                        LOC2_FEMT2OLD(:, CV_KLOC) = FEMT2OLD_ALL(:, CV_KNOD2)
                                    END IF
                                END DO
                            END IF ! ENDOF IF (between_elements) THEN
                            IF ( on_domain_boundary ) THEN
                                DO U_SKLOC = 1, Mdims%u_snloc
                                    U_KLOC = U_SLOC2LOC( U_SKLOC )
                                    U_NODK = ndgln%u(( ELE - 1 ) * Mdims%u_nloc + U_KLOC )
                                    U_SNODK = ( SELE - 1 ) * Mdims%u_snloc + U_SKLOC
                                    SLOC_NU(:, :, U_SKLOC) = NU_ALL(:, :, U_NODK)
                                    SLOC_NUOLD(:, :, U_SKLOC) = NUOLD_ALL(:, :, U_NODK)
                                END DO
                            END IF
                            !------------------
                            If_GOT_DIFFUS2: IF ( GOT_DIFFUS ) THEN
                                ! This sub caculates the effective diffusion
                                ! coefficient DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX
                                T_ALL_J( : )   =T_ALL(:, CV_NODJ)
                                TOLD_ALL_J( : )=TOLD_ALL(:, CV_NODJ)
                                LOC_WIC_T_BC_ALL(:)=0
                                !           IF(SELE.NE.0) THEN
                                IF(on_domain_boundary) THEN
                                    DO IPHASE=1,Mdims%nphase
                                        LOC_WIC_T_BC_ALL(IPHASE)=WIC_T_BC_ALL(1, IPHASE, SELE)
                                        IF(LOC_WIC_T_BC_ALL(IPHASE)==WIC_T_BC_DIRICHLET) THEN
                                            T_ALL_J( IPHASE ) = SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + Mdims%cv_snloc*( SELE- 1) )
                                            TOLD_ALL_J( IPHASE )=SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + Mdims%cv_snloc*( SELE- 1) )
                                        ENDIF
                                    END DO
                                ENDIF
                                CALL DIFFUS_CAL_COEFF( DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX,  &
                                    Mdims%cv_nloc, Mdims%mat_nloc, Mdims%nphase, ndgln%mat, &
                                    CV_funs%scvfen, CV_funs%scvfen, GI, Mdims%ndim, TDIFFUSION, DUMMY_ZERO_NDIM_NDIM_NPHASE, &
                                    HDC, &
                                    T_ALL_J( : ), T_ALL(:, CV_NODI), &
                                    TOLD_ALL_J( : ), TOLD_ALL(:, CV_NODI), &
                                    ELE, ELE2, CVNORMX_ALL( :, GI ), &
                                    DTX_ELE_ALL(:,:,:,ELE), DTOLDX_ELE_ALL(:,:,:,ELE),  DTX_ELE_ALL(:,:,:,MAX(1,ELE2)), DTOLDX_ELE_ALL(:,:,:,MAX(ELE2,1)), &
                                    LOC_WIC_T_BC_ALL, CV_OTHER_LOC, MAT_OTHER_LOC, Mdims%cv_snloc, CV_SLOC2LOC, &
                                    on_domain_boundary, between_elements )
                            ELSE
                                DIFF_COEF_DIVDX = 0.0
                                DIFF_COEFOLD_DIVDX = 0.0
                            END IF If_GOT_DIFFUS2
                            NFACE_ITS = 1
                            FACE_ITS = 1
                            ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI.
                            !Calling the functions directly instead inside a wrapper saves a around a 5%
                            IF( GOT_T2 ) THEN
                                IF( is_porous_media ) THEN
                                    CALL GET_INT_VEL_POROUS_VEL( NDOTQNEW, NDOTQOLD, INCOMEOLD, &
                                        T2OLD_ALL(:, CV_NODI), T2OLD_ALL(:, CV_NODJ), LOC_FEMT2OLD, &
                                        LOC_NUOLD, LOC2_NUOLD, SLOC_NUOLD, &
                                        UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
                                        upwnd%adv_coef(:,:,:, MAT_NODI), upwnd%adv_coef_grad(:,:,:, MAT_NODI), &
                                        upwnd%adv_coef(:,:,:, MAT_NODJ), upwnd%adv_coef_grad(:,:,:, MAT_NODJ), &
                                        upwnd%inv_adv_coef(:,:,:,MAT_NODI), upwnd%inv_adv_coef(:,:,:,MAT_NODJ), &
                                        NUOLDGI_ALL, MASS_CV(CV_NODI), MASS_CV(CV_NODJ), &
                                        T2OLDUPWIND_MAT_ALL( :, COUNT_IN), T2OLDUPWIND_MAT_ALL( :, COUNT_OUT), &
                                        .false., anisotropic_and_frontier)
                                    CALL GET_INT_VEL_POROUS_VEL( NDOTQNEW, NDOTQ, INCOME, &
                                        T2_ALL(:, CV_NODI), T2_ALL(:, CV_NODJ), LOC_FEMT2, &
                                        LOC_NU, LOC2_NU, SLOC_NU, &
                                        UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
                                        upwnd%adv_coef(:,:,:, MAT_NODI), upwnd%adv_coef_grad(:,:,:, MAT_NODI), &
                                        upwnd%adv_coef(:,:,:, MAT_NODJ), upwnd%adv_coef_grad(:,:,:, MAT_NODJ), &
                                        upwnd%inv_adv_coef(:,:,:,MAT_NODI), upwnd%inv_adv_coef(:,:,:,MAT_NODJ), &
                                        NUGI_ALL, MASS_CV(CV_NODI), MASS_CV(CV_NODJ), &
                                        T2UPWIND_MAT_ALL( :, COUNT_IN), T2UPWIND_MAT_ALL( :, COUNT_OUT), &
                                        .true., anisotropic_and_frontier)
                                else
                                    call GET_INT_VEL_ORIG_NEW( NDOTQNEW, NDOTQOLD, INCOMEOLD, &
                                        T2OLD_ALL(:, CV_NODI), T2OLD_ALL(:, CV_NODJ), DENOLD_ALL(:, CV_NODI), DENOLD_ALL(:, CV_NODJ), &
                                        LOC_NUOLD, LOC2_NUOLD, NUOLDGI_ALL, &
                                        UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, .false. )
                                    call GET_INT_VEL_ORIG_NEW( NDOTQNEW, NDOTQ, INCOME, &
                                        T2_ALL(:, CV_NODI), T2_ALL(:, CV_NODJ), DEN_ALL(:, CV_NODI), DEN_ALL(:, CV_NODJ), &
                                        LOC_NU, LOC2_NU, NUGI_ALL, &
                                        UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, .true. )
                                end if
                            ELSE
                                IF( is_porous_media ) THEN
                                    CALL GET_INT_VEL_POROUS_VEL( NDOTQNEW, NDOTQOLD, INCOMEOLD, &
                                        TOLD_ALL(:, CV_NODI), TOLD_ALL(:, CV_NODJ), LOC_FEMTOLD, &
                                        LOC_NUOLD, LOC2_NUOLD, SLOC_NUOLD, &
                                        UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
                                        upwnd%adv_coef(:,:,:, MAT_NODI), upwnd%adv_coef_grad(:,:,:, MAT_NODI), &
                                        upwnd%adv_coef(:,:,:, MAT_NODJ), upwnd%adv_coef_grad(:,:,:, MAT_NODJ), &
                                        upwnd%inv_adv_coef(:,:,:,MAT_NODI), upwnd%inv_adv_coef(:,:,:,MAT_NODJ), &
                                        NUOLDGI_ALL, MASS_CV(CV_NODI), MASS_CV(CV_NODJ), &
                                        TOLDUPWIND_MAT_ALL( :, COUNT_IN), TOLDUPWIND_MAT_ALL( :, COUNT_OUT), &
                                        .false., anisotropic_and_frontier)
                                    CALL GET_INT_VEL_POROUS_VEL( NDOTQNEW, NDOTQ, INCOME, &
                                        T_ALL(:, CV_NODI), T_ALL(:, CV_NODJ), LOC_FEMT, &
                                        LOC_NU, LOC2_NU, SLOC_NU, &
                                        UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
                                        upwnd%adv_coef(:,:,:, MAT_NODI), upwnd%adv_coef_grad(:,:,:, MAT_NODI), &
                                        upwnd%adv_coef(:,:,:, MAT_NODJ), upwnd%adv_coef_grad(:,:,:, MAT_NODJ), &
                                        upwnd%inv_adv_coef(:,:,:,MAT_NODI), upwnd%inv_adv_coef(:,:,:,MAT_NODJ), &
                                        NUGI_ALL, MASS_CV(CV_NODI), MASS_CV(CV_NODJ), &
                                        TUPWIND_MAT_ALL( :, COUNT_IN), TUPWIND_MAT_ALL( :, COUNT_OUT), &
                                        .true., anisotropic_and_frontier)
                                else
                                    call GET_INT_VEL_ORIG_NEW( NDOTQNEW, NDOTQOLD, INCOMEOLD, &
                                        TOLD_ALL(:, CV_NODI), TOLD_ALL(:, CV_NODJ), DENOLD_ALL(:, CV_NODI), DENOLD_ALL(:, CV_NODJ), &
                                        LOC_NUOLD, LOC2_NUOLD, NUOLDGI_ALL, &
                                        UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, .false. )
                                    call GET_INT_VEL_ORIG_NEW( NDOTQNEW, NDOTQ, INCOME, &
                                        T_ALL(:, CV_NODI), T_ALL(:, CV_NODJ), DEN_ALL(:, CV_NODI), DEN_ALL(:, CV_NODJ), &
                                        LOC_NU, LOC2_NU, NUGI_ALL, UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, .true. )
                                end if
                            ENDIF
                            !Obtain income for cv_nodj
                            !When NDOTQ == 0, INCOME_J has to be 1 as well, not 0
                            WHERE ( NDOTQ <= 0. )
                                INCOME_J = 0.
                            ELSE WHERE
                                INCOME_J = 1.
                            END WHERE
                            !Calculate the courant number for porous media
                            if (present(Courant_number) .and. is_porous_media.and. .not. on_domain_boundary) then
                                !ndotq = velocity * normal                     !In the wells the flow is too fast and makes this misleading
                                Courant_number(1) = max(Courant_number(1), abs ( dt * maxval(ndotq(1:Mdims%n_in_pres)) / (VOLFRA_PORE( 1, ELE ) * hdc)))
                                !and the shock-front Courant number
                                if (shock_front_in_ele(ele, Mdims, T_ALL, ndgln, Imble_frac(:, ELE))) then
                                    !ndotq = velocity * normal
                                    Courant_number(2) = max(Courant_number(2), abs ( dt * maxval(ndotq(1:Mdims%n_in_pres)) / (VOLFRA_PORE( 1, ELE ) * hdc)))
                                end if
                            end if
                            If_GOT_CAPDIFFUS: IF ( capillary_pressure_activated ) THEN
                                IF(SELE == 0) THEN
                                    CAP_DIFF_COEF_DIVDX = 0.
                                    do iphase =1, Mdims%nphase
                                        rsum_nodi(iphase) = dot_product(CVNORMX_ALL(:, GI), matmul(upwnd%inv_adv_coef(:,:,iphase,MAT_NODI),&
                                            CVNORMX_ALL(:, GI)))
                                        rsum_nodj(iphase) = dot_product(CVNORMX_ALL(:, GI), matmul(upwnd%inv_adv_coef(:,:,iphase,MAT_NODJ),&
                                            CVNORMX_ALL(:, GI) ))
                                    end do!If we are using the non-consistent capillary pressure we want to use the central...
                                    IF( UPWIND_CAP_DIFFUSION ) THEN!...method to encourage a normal diffusion
                                        CAP_DIFF_COEF_DIVDX = (CAP_DIFFUSION( :, MAT_NODI )&
                                            * rsum_nodi*(1.-INCOME(:))  +&
                                            CAP_DIFFUSION( :, MAT_NODJ ) * rsum_nodj * INCOME) /HDC

                                    ELSE ! Central difference...
                                        CAP_DIFF_COEF_DIVDX( : ) =  -OvRelax_param(CV_NODI) * (0.5*(T_ALL( :, CV_NODI )&
                                            * rsum_nodi(:) + T_ALL( :, CV_NODJ ) * rsum_nodj(:) ) /HDC)
!                                        CAP_DIFF_COEF_DIVDX( : ) = 0.5*(CAP_DIFFUSION( :, MAT_NODI )&
!                                            * rsum_nodi(:) + CAP_DIFFUSION( :, MAT_NODJ ) * rsum_nodj(:) ) /HDC
                                    ENDIF
                                ELSE
                                    CAP_DIFF_COEF_DIVDX( : ) = 0.0
                                ENDIF
                                !Distribute the capillary coefficient over the phases to ensure mass conservation
                                !This is very important as it allows to use the over-relaxation parameter safely
                                !and reduce the cost of using capillary pressure in several orders of magnitude
                                CAP_DIFF_COEF_DIVDX(1:Mdims%n_in_pres) =  CAP_DIFF_COEF_DIVDX(phase_with_pc)/Mdims%n_in_pres
                            ELSE
                                CAP_DIFF_COEF_DIVDX( : ) = 0.0
                            END IF If_GOT_CAPDIFFUS

                            ! Pack ndotq information:
                            IPT=1
                            CALL PACK_LOC( F_INCOME(:), INCOME( : ),    Mdims%nphase, IPT, IGOT_T_PACK(:,1) ) ! t
                            CALL PACK_LOC( F_INCOME(:), INCOMEOLD( : ), Mdims%nphase, IPT, IGOT_T_PACK(:,2) ) ! TOLD
                            CALL PACK_LOC( F_INCOME(:), INCOME( : ),    Mdims%nphase, IPT, IGOT_T_PACK(:,3) )  ! d
                            CALL PACK_LOC( F_INCOME(:), INCOMEOLD( : ), Mdims%nphase, IPT, IGOT_T_PACK(:,4) )  ! DOLD
                            IF ( use_volume_frac_T2 ) THEN
                                CALL PACK_LOC( F_INCOME(:), INCOME( : ),    Mdims%nphase, IPT, IGOT_T_PACK(:,5) )  ! T2
                                CALL PACK_LOC( F_INCOME(:), INCOMEOLD( : ), Mdims%nphase, IPT, IGOT_T_PACK(:,6) )  ! T2OLD
                            ENDIF
                            IPT=1
                            CALL PACK_LOC( F_NDOTQ(:), NDOTQ( : ),    Mdims%nphase, IPT, IGOT_T_PACK(:,1) ) ! t
                            CALL PACK_LOC( F_NDOTQ(:), NDOTQOLD( : ), Mdims%nphase, IPT, IGOT_T_PACK(:,2) ) ! TOLD
                            CALL PACK_LOC( F_NDOTQ(:), NDOTQ( : ),    Mdims%nphase, IPT, IGOT_T_PACK(:,3) )  ! d
                            CALL PACK_LOC( F_NDOTQ(:), NDOTQOLD( : ), Mdims%nphase, IPT, IGOT_T_PACK(:,4) )  ! DOLD
                            IF ( use_volume_frac_T2 ) THEN
                                CALL PACK_LOC( F_NDOTQ(:), NDOTQ( : ),    Mdims%nphase, IPT, IGOT_T_PACK(:,5) )  ! T2
                                CALL PACK_LOC( F_NDOTQ(:), NDOTQOLD( : ), Mdims%nphase, IPT, IGOT_T_PACK(:,6) )  ! T2OLD
                            ENDIF
                                  !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
                                  ! Calculate T and DEN on the CV face at quadrature point GI.
                            IF(NFIELD.GT.0) THEN
                                CALL GET_INT_T_DEN_new( LIMF )
                            ENDIF

                            IF(APPLY_ENO) THEN
                                ! Calculate DETWEI, RA, NX, NY, NZ for element ELE
                                call DETNLXR_INVJAC( ELE, X_ALL, ndgln%x, CV_funs%scvfeweigh, CV_funs%scvfen, CV_funs%scvfenlx_all, SdevFuns)
                                ! Apply a simple ENO scheme to T,TOLD only which is not bounded but gets rid of most of the osillations.
                                ! Put the results in LIMF.
                                CALL APPLY_ENO_2_T(LIMF, T_ALL,TOLD_ALL, FEMT_ALL,FEMTOLD_ALL, INCOME,INCOMEOLD, IGOT_T_PACK, &
                                    CV_NODI, CV_NODJ, X_NODI, X_NODJ, CV_ILOC, CV_JLOC, &
                                    ELE, Mdims%cv_nonods, Mdims%ndim, Mdims%nphase,  &
                                    Mdims%cv_nloc,Mdims%totele, ndgln%x, ndgln%cv,  &
                                    X_ALL,FACE_ELE,CV_GIdims%nface,BETWEEN_ELEMENTS, CV_funs%scvfen, SdevFuns%NX_ALL, GI, SdevFuns%INV_JAC, &
                                    NUGI_ALL, on_domain_boundary )
                            ENDIF

                            ! it does not matter about bcs for FVT below as its zero'ed out in the eqns:
                            FVT(:)=T_ALL(:,CV_NODI)*(1.0-INCOME(:)) + T_ALL(:,CV_NODJ)*INCOME(:)
                            !FVD(:)=DEN_ALL(:,CV_NODI)*(1.0-INCOME(:)) + DEN_ALL(:,CV_NODJ)*INCOME(:)
                            ! Generate some local F variables ***************
                            ! loc_f - Unpack into the limiting variables LIMT and may be store them in the cache.
                            IPT=1
                            CALL UNPACK_LOC( LIMF(:), LIMT( : ),    Mdims%nphase, IPT, IGOT_T_PACK(:,1), IGOT_T_CONST(:,1), IGOT_T_CONST_VALUE(:,1))
                            CALL UNPACK_LOC( LIMF(:), LIMTOLD( : ), Mdims%nphase, IPT, IGOT_T_PACK(:,2), IGOT_T_CONST(:,2), IGOT_T_CONST_VALUE(:,2))
                            CALL UNPACK_LOC( LIMF(:), LIMD( : ),    Mdims%nphase, IPT, IGOT_T_PACK(:,3), IGOT_T_CONST(:,3), IGOT_T_CONST_VALUE(:,3))
                            CALL UNPACK_LOC( LIMF(:), LIMDOLD( : ), Mdims%nphase, IPT, IGOT_T_PACK(:,4), IGOT_T_CONST(:,4), IGOT_T_CONST_VALUE(:,4))
                            IF ( use_volume_frac_T2 ) THEN
                                CALL UNPACK_LOC( LIMF(:), LIMT2( : ),    Mdims%nphase, IPT, IGOT_T_PACK(:,5), IGOT_T_CONST(:,5), IGOT_T_CONST_VALUE(:,5))
                                CALL UNPACK_LOC( LIMF(:), LIMT2OLD( : ), Mdims%nphase, IPT, IGOT_T_PACK(:,6), IGOT_T_CONST(:,6), IGOT_T_CONST_VALUE(:,6))
                            else
                                LIMT2( : )=1.0; LIMT2OLD( : )=1.0
                            ENDIF
                            LIMDT=LIMD*LIMT
                            LIMDTOLD=LIMDOLD*LIMTOLD
                            LIMDTT2=LIMD*LIMT*LIMT2
                            LIMDTT2OLD=LIMDOLD*LIMTOLD*LIMT2OLD
                            ! Generate some local F variables ***************...
                            ! Make allowances for no matrix stencil operating from outside the boundary.
                            BCZERO=1.0
                            IF( on_domain_boundary ) BCZERO=1.0-INCOME
                            ! Define face value of theta
                            IF ( GOT_T2 ) THEN
                                FTHETA(:) = FACE_THETA_MANY( DT, CV_THETA, ( CV_DISOPT>=8 ), HDC, Mdims%nphase, &
                                    Mdims%n_in_pres, NDOTQ, LIMDTT2, DIFF_COEF_DIVDX, &
                                    T_ALL(:, CV_NODJ) * DEN_ALL(:, CV_NODJ) * T2_ALL(:, CV_NODJ), &
                                    T_ALL(:, CV_NODI) * DEN_ALL(:, CV_NODI) * T2_ALL(:, CV_NODI), &
                                    NDOTQOLD, LIMDTT2OLD, DIFF_COEFOLD_DIVDX, &
                                    TOLD_ALL(:, CV_NODJ) * DENOLD_ALL(:, CV_NODJ) * T2OLD_ALL(:, CV_NODJ), &
                                    TOLD_ALL(:, CV_NODI) * DENOLD_ALL(:, CV_NODI) * T2OLD_ALL(:, CV_NODI) )
                            ELSE
                                FTHETA(:) = FACE_THETA_MANY( DT, CV_THETA, ( CV_DISOPT>=8 ), HDC, Mdims%nphase, &
                                    Mdims%n_in_pres,NDOTQ, LIMDTT2, DIFF_COEF_DIVDX, &
                                    T_ALL(:, CV_NODJ) * DEN_ALL(:, CV_NODJ), &
                                    T_ALL(:, CV_NODI) * DEN_ALL(:, CV_NODI), &
                                    NDOTQOLD, LIMDTT2OLD, DIFF_COEFOLD_DIVDX, &
                                    TOLD_ALL(:, CV_NODJ) * DENOLD_ALL(:, CV_NODJ), &
                                    TOLD_ALL(:, CV_NODI) * DENOLD_ALL(:, CV_NODI)  )
                            END IF
                            ! adjust the value of FTHETA for use with velocity only so we can use FTHETA=0.0 for voln frac.
                            ! THETA_VEL_HAT=0.0 does not change NDOTQOLD, THETA_VEL_HAT=1.0 sets NDOTQOLD=NDOTQNEW.
                            ! If THETA_VEL_HAT<0.0 then automatically choose THETA_VEL to be as close to THETA_VEL_HAT (e.g.=0) as possible.
                            ! This determins how implicit velocity is in the cty eqn.
                            THETA_VEL(:)=THETA_VEL_HAT
                            IF(THETA_VEL_HAT<0.0) THEN
                                DO IPHASE=1,Mdims%nphase
                                    IF(FTHETA(IPHASE)>=0.5) THEN
                                        THETA_VEL(IPHASE) = ABS(THETA_VEL_HAT)
                                    ELSE
                                        THETA_VEL(IPHASE) = MAX( (0.5-FTHETA(IPHASE))/(1.0-FTHETA(IPHASE)),  ABS(THETA_VEL_HAT))
                                    ENDIF
                                END DO
                            ENDIF
                            NDOTQOLD = THETA_VEL*NDOTQ + (1.0-THETA_VEL)*NDOTQOLD
                            FTHETA_T2(:) = FTHETA(:) * LIMT2(:)
                            ONE_M_FTHETA_T2OLD(:) = (1.0-FTHETA(:)) * LIMT2OLD(:)
                            FTHETA_T2_J(:) = FTHETA_T2(:)!FTHETA(:) * LIMT2(:)
                            ONE_M_FTHETA_T2OLD_J(:) = ONE_M_FTHETA_T2OLD(:)!(1.0-FTHETA(:)) * LIMT2OLD(:)
                            IF(IGOT_THETA_FLUX == 1) THEN
                                IF ( GET_THETA_FLUX ) THEN
                                    THETA_FLUX( :, GLOBAL_FACE ) = FTHETA(:) * LIMDT(:) / DEN_ALL_DIVID(:, CV_NODI)
                                    ONE_M_THETA_FLUX( :, GLOBAL_FACE ) = (1.0-FTHETA(:)) * LIMDTOLD(:) / DEN_ALL_DIVID(:, CV_NODI)
                                    if(integrate_other_side) then ! for the flux on the other side of the CV face...
                                        THETA_FLUX_J( :, GLOBAL_FACE ) = FTHETA(:) * LIMDT(:) / DEN_ALL_DIVID(:, CV_NODJ)
                                        ONE_M_THETA_FLUX_J( :, GLOBAL_FACE ) = (1.0-FTHETA(:)) * LIMDTOLD(:) / DEN_ALL_DIVID(:, CV_NODJ)
                                    endif
                                END IF
                                IF ( USE_THETA_FLUX ) THEN
                                    FTHETA_T2(:) = THETA_FLUX( :, GLOBAL_FACE )
                                    ONE_M_FTHETA_T2OLD(:) = ONE_M_THETA_FLUX( :, GLOBAL_FACE )
                                    if(integrate_other_side) then ! for the flux on the other side of the CV face...
                                        FTHETA_T2_J(:) = THETA_FLUX_J( :, GLOBAL_FACE )
                                        ONE_M_FTHETA_T2OLD_J(:) = ONE_M_THETA_FLUX_J( :, GLOBAL_FACE )
                                    endif
                                END IF
                            END IF
                            !====================== ACV AND RHS ASSEMBLY ===================
                            Conditional_GETCT2: IF ( GETCT ) THEN ! Obtain the CV discretised Mmat%CT eqations plus RHS
                                IF(CT_DO_NOT_CHANGE) THEN ! Mmat%CT will not change with this option...
                                    FTHETA_T2=1.0
                                    ONE_M_FTHETA_T2OLD=0.0
                                    FTHETA_T2_J=1.0
                                    ONE_M_FTHETA_T2OLD_J=0.0
                                ENDIF
                                IF ( got_free_surf ) THEN
                                    IF ( on_domain_boundary ) THEN
                                        IF ( .not.symmetric_P ) THEN
                                            IF ( WIC_P_BC_ALL( 1,1,SELE ) == WIC_P_BC_FREE ) THEN ! on the free surface...
                                                DO P_JLOC = 1, Mdims%cv_nloc
                                                    P_JNOD = ndgln%cv( (ELE-1)*Mdims%cv_nloc + P_JLOC )
                                                    ! Use the same sparcity as the MN matrix...
                                                    COUNT_SUF = 0
                                                    DO COUNT = Mspars%CMC%fin( CV_NODI ), Mspars%CMC%fin( CV_NODI + 1 ) - 1
                                                        IF ( Mspars%CMC%col( COUNT ) == P_JNOD ) THEN
                                                            COUNT_SUF = COUNT
                                                            EXIT
                                                        END IF
                                                    END DO
                                                    MM_GRAVTY = CVNORMX_ALL( 3, GI ) * CV_funs%scvfen( P_JLOC, GI ) * SdevFuns%DETWEI(GI) / ( DT**2 * GRAVTY )
                                                    MASS_SUF( COUNT_SUF ) = MASS_SUF( COUNT_SUF ) + MM_GRAVTY
                                                END DO
                                            END IF
                                        END IF
                                    END IF
                                END IF
                                ct_rhs_phase_cv_nodi=0.0; ct_rhs_phase_cv_nodj=0.0
                                CALL PUT_IN_CT_RHS(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV, ct_rhs_phase_cv_nodi, ct_rhs_phase_cv_nodj, &
                                    Mdims, CV_funs, ndgln, Mmat, GI,  &
                                    between_elements, on_domain_boundary, ELE, ELE2, SELE, HDC, MASS_ELE, &
                                    JCOUNT_KLOC, JCOUNT_KLOC2, ICOUNT_KLOC, ICOUNT_KLOC2, C_JCOUNT_KLOC, C_JCOUNT_KLOC2, C_ICOUNT_KLOC, C_ICOUNT_KLOC2, U_OTHER_LOC,  U_SLOC2LOC, CV_SLOC2LOC,&
                                    SdevFuns%DETWEI, CVNORMX_ALL, DEN_ALL_DIVID, CV_NODI, CV_NODJ, &
                                    WIC_U_BC_ALL, WIC_P_BC_ALL, pressure_BCs%val, &
                                    UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL,  &
                                    NDOTQNEW, NDOTQOLD, LIMT, LIMDT, LIMDTOLD, LIMT_HAT, NDOTQ_HAT, &
                                    FTHETA_T2, ONE_M_FTHETA_T2OLD, FTHETA_T2_J, ONE_M_FTHETA_T2OLD_J, integrate_other_side_and_not_boundary, &
                                    theta_cty_solid, loc_u, THETA_VEL, &
                                    rdum_ndim_nphase_1, rdum_nphase_1, rdum_nphase_2, rdum_nphase_3)
                                do ipres=1,Mdims%npres
                                    call addto(Mmat%CT_RHS,ipres,cv_nodi,sum(ct_rhs_phase_cv_nodi(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres) ))
                                    if ( integrate_other_side_and_not_boundary ) then
                                        call addto(Mmat%CT_RHS,ipres,cv_nodj,sum(ct_rhs_phase_cv_nodj(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres) ))
                                    end if
                                end do

                                !Store fluxes across all the boundaries either for mass conservation check or mass outflux
                                !velocity * area * density * saturation
                                if (on_domain_boundary ) then
                                    if (surface_element_owned(old_tracer, sele)) then
                                        !Store total outflux
                                        bcs_outfluxes(1:Mdims%n_in_pres, CV_NODI, 0) =  bcs_outfluxes(1:Mdims%n_in_pres, CV_NODI,0) + &
                                            ndotqnew(1:Mdims%n_in_pres) * SdevFuns%DETWEI(gi) * LIMT(1:Mdims%n_in_pres)
                                        if (outfluxes%calculate_flux)  then
                                            do k = 1, size(outfluxes%outlet_id)!here below we just need a saturation
                                                if (integrate_over_surface_element(old_tracer, sele, (/outfluxes%outlet_id(k)/))) then
                                                    bcs_outfluxes(1:Mdims%n_in_pres, CV_NODI, k) =  bcs_outfluxes(1:Mdims%n_in_pres, CV_NODI, k) + &
                                                        ndotqnew(1:Mdims%n_in_pres) * SdevFuns%DETWEI(gi) * LIMT(1:Mdims%n_in_pres)
                                                    if (has_temperature) then
                                                        do iphase = 1, Mdims%nphase
                                                            outfluxes%totout(2, iphase, k) =  max(  temp_field%val(1,iphase,CV_NODI),&
                                                                 outfluxes%totout(2, iphase, k)   )
                                                        end do
                                                    end if
                                                end if
                                            end do
                                        end if
                                    end if
                                end if

                            ENDIF Conditional_GETCT2
                            Conditional_GETCV_DISC: IF ( GETCV_DISC ) THEN
                                ! Obtain the CV discretised advection/diffusion equations
                                ROBIN1=0.0; ROBIN2=0.0
                                IF( on_domain_boundary ) then
                                    where ( WIC_T_BC_ALL(1,:,SELE) == WIC_T_BC_ROBIN )
                                        ! this needs to be corrected (its correct but misleading)...
                                        ROBIN1(:) = SUF_T_BC_ROB1_ALL(1,:, CV_SILOC+Mdims%cv_snloc*(sele-1))
                                        ROBIN2(:) = SUF_T_BC_ROB2_ALL(1,:, CV_SILOC+Mdims%cv_snloc*(sele-1))
                                    end where
                                END IF
                                LOC_CV_RHS_I(:)=0.0
                                LOC_CV_RHS_J(:)=0.0
                                IF ( GETMAT ) THEN
                                    ! - Calculate the integration of the limited, high-order flux over a face
                                    ! Conservative discretisation. The matrix (PIVOT ON LOW ORDER SOLN)
                                    IF ( on_domain_boundary ) THEN
                                        DO IPHASE=1,Mdims%nphase
                                            IF(WIC_T_BC_ALL(1,iphase,sele) == WIC_T_BC_DIRICHLET) THEN
                                                LOC_CV_RHS_I( IPHASE ) =  LOC_CV_RHS_I( IPHASE ) &
                                                    + FTHETA(IPHASE) * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE) &
                                                    * SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + Mdims%cv_snloc*( SELE- 1))
                                                IF(GET_GTHETA) THEN
                                                    THETA_GDIFF( IPHASE, CV_NODI ) =  THETA_GDIFF( IPHASE, CV_NODI ) &
                                                        + FTHETA(IPHASE) * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE) &
                                                        * SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + Mdims%cv_snloc*( SELE- 1) )
                                                END IF
                                            END IF
                                        END DO
                                    ELSE
                                        do iphase=1,Mdims%nphase
                                            !temporary to check a possible memory problem with Valgrind!!!
                                            auxR = FTHETA_T2(iphase) * SdevFuns%DETWEI( GI ) * NDOTQNEW(iphase) * INCOME(iphase) * LIMD(iphase)
                                            call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodi,cv_nodj,auxR) ! Advection
                                            if (GOT_DIFFUS) call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodi,cv_nodj,&
                                                            - FTHETA(iphase) * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX(iphase))
                                            if (capillary_pressure_activated) call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodi,cv_nodj,&
                                                            - SdevFuns%DETWEI( GI ) * CAP_DIFF_COEF_DIVDX(iphase))
                                        end do
                                         ! integrate the other CV side contribution (the sign is changed)...
                                        if(integrate_other_side_and_not_boundary) then
                                            do iphase=1,Mdims%nphase
                                                call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodj,cv_nodi,&
                                                    - FTHETA_T2_J(IPHASE) * SdevFuns%DETWEI( GI ) * NDOTQNEW(IPHASE) * INCOME_J(IPHASE) * LIMD(IPHASE) ) ! Advection
                                            if (GOT_DIFFUS) call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodj,cv_nodi,&
                                                    - FTHETA(IPHASE) * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE))
                                            if (capillary_pressure_activated) call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodj,cv_nodi,&
                                                           - SdevFuns%DETWEI( GI ) * CAP_DIFF_COEF_DIVDX(IPHASE))
                                            end do
                                        endif
                                        IF ( GET_GTHETA ) THEN
                                            THETA_GDIFF( :, CV_NODI ) =  THETA_GDIFF( :, CV_NODI ) &
                                                + FTHETA(:) * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX(:) * T_ALL(:, CV_NODJ) ! Diffusion contribution
                                            ! integrate the other CV side contribution (the sign is changed)...
                                            if(integrate_other_side_and_not_boundary) then
                                                THETA_GDIFF( :, CV_NODJ ) =  THETA_GDIFF( :, CV_NODJ ) &
                                                    + FTHETA(:) * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX(:) * T_ALL(:, CV_NODI) ! Diffusion contribution
                                            endif
                                        END IF
                                    END IF ! endif of IF ( on_domain_boundary ) THEN ELSE
                                    do iphase=1,Mdims%nphase
                                        call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodi,cv_nodi,&
                                            +  FTHETA_T2(iphase) * SdevFuns%DETWEI( GI ) * NDOTQNEW(iphase) * ( 1. - INCOME(iphase) ) * LIMD(iphase) ) ! Advection

                                        if (GOT_DIFFUS) call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodi,cv_nodi,&
                                           +  FTHETA(iphase) * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX(iphase))
                                        if (capillary_pressure_activated) call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodi,cv_nodi,&
                                           +  SdevFuns%DETWEI( GI ) * CAP_DIFF_COEF_DIVDX(iphase))
                                        if (.not.conservative_advection) call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodi,cv_nodi,&
                                           - FTHETA_T2(iphase) * ( ONE_M_CV_BETA ) * SdevFuns%DETWEI( GI ) * NDOTQNEW(iphase) * LIMD(iphase))
                                        if (on_domain_boundary) call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodi,cv_nodi,&
                                                                SdevFuns%DETWEI( GI ) * ROBIN1(iphase))
                                    end do
                                    if(integrate_other_side_and_not_boundary) then
                                        do iphase=1,Mdims%nphase
                                            call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodj,cv_nodj,&
                                                -  FTHETA_T2_J(iphase) * SdevFuns%DETWEI( GI ) * NDOTQNEW(iphase) * ( 1. - INCOME_J(iphase) ) * LIMD(iphase) ) ! Advection

                                            if (GOT_DIFFUS) call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodj,cv_nodj,&
                                                +  FTHETA(iphase) * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX(iphase))
                                            if (capillary_pressure_activated) call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodj,cv_nodj,&
                                                +  SdevFuns%DETWEI( GI ) * CAP_DIFF_COEF_DIVDX(iphase))
                                            if (.not.conservative_advection) call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodj,cv_nodj,&
                                                + FTHETA_T2_J(iphase) * ( ONE_M_CV_BETA ) * SdevFuns%DETWEI( GI ) * NDOTQNEW(iphase) * LIMD(iphase))

                                        end do
                                    endif
                                    IF ( GET_GTHETA ) THEN
                                        THETA_GDIFF( :, CV_NODI ) =  THETA_GDIFF( :, CV_NODI ) &
                                            -  FTHETA(:) * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX(:) * T_ALL(:, CV_NODI) & ! Diffusion contribution
                                            -  SdevFuns%DETWEI( GI ) * ROBIN1(:) * T_ALL(:, CV_NODI)  ! Robin bc
                                        if(integrate_other_side_and_not_boundary) then
                                            THETA_GDIFF( :, CV_NODJ ) =  THETA_GDIFF( :, CV_NODJ ) &
                                                -  FTHETA(:) * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX(:) * T_ALL(:, CV_NODJ) ! Diffusion contribution
                                        endif
                                    END IF
                                END IF  ! ENDOF IF ( GETMAT ) THEN

                                ! Put results into the RHS vector
                                LOC_CV_RHS_I( : ) =  LOC_CV_RHS_I( : )  &
                                       ! subtract 1st order adv. soln.
                                    + FTHETA_T2(:) * NDOTQNEW(:) * SdevFuns%DETWEI( GI ) * LIMD(:) * FVT(:) * BCZERO(:) &
                                    -  SdevFuns%DETWEI( GI ) * ( FTHETA_T2(:) * NDOTQNEW(:) * LIMDT(:) &
                                    + ONE_M_FTHETA_T2OLD(:)* NDOTQOLD(:) * LIMDTOLD(:) ) ! hi order adv
                                ! Subtract out 1st order term non-conservative adv.
                                    if (GOT_DIFFUS) LOC_CV_RHS_I( : ) =  LOC_CV_RHS_I( : ) &
                                        + (1.-FTHETA(:)) * SdevFuns%DETWEI(GI) * DIFF_COEFOLD_DIVDX(:) &
                                        * ( TOLD_ALL(:, CV_NODJ) - TOLD_ALL(:, CV_NODI) )
                                    if (capillary_pressure_activated) LOC_CV_RHS_I( : ) =  LOC_CV_RHS_I( : ) &
                                        - SdevFuns%DETWEI(GI) * CAP_DIFF_COEF_DIVDX(:) &  ! capillary pressure stabilization term..
                                        * ( T_ALL(:, CV_NODJ) - T_ALL(:, CV_NODI) )
                                    if (.not.conservative_advection) LOC_CV_RHS_I( : ) =  LOC_CV_RHS_I( : ) &
                                        - FTHETA_T2(:) * ( ONE_M_CV_BETA ) * SdevFuns%DETWEI( GI ) * NDOTQNEW(:) * LIMD(:) * T_ALL(:, CV_NODI) &
                                        + ( ONE_M_CV_BETA) * SdevFuns%DETWEI( GI ) &
                                        * ( FTHETA_T2(:) * NDOTQNEW(:) * T_ALL(:, CV_NODI) * LIMD(:)  &
                                        + ONE_M_FTHETA_T2OLD(:) * NDOTQOLD(:) * LIMDOLD(:) * TOLD_ALL(:, CV_NODI) )
                                    if (on_domain_boundary) LOC_CV_RHS_I( : ) =  LOC_CV_RHS_I( : ) &
                                        + SdevFuns%DETWEI( GI ) * ROBIN2(:)

                                if(integrate_other_side_and_not_boundary) then
                                    LOC_CV_RHS_J( : ) =  LOC_CV_RHS_J( : )  &
                                           ! subtract 1st order adv. soln.
                                        - FTHETA_T2_J(:) * NDOTQNEW(:) * SdevFuns%DETWEI( GI ) * LIMD(:) * FVT(:) * BCZERO(:) &
                                        +  SdevFuns%DETWEI( GI ) * ( FTHETA_T2_J(:) * NDOTQNEW(:) * LIMDT(:) &
                                        + ONE_M_FTHETA_T2OLD_J(:) * NDOTQOLD(:) * LIMDTOLD(:) )
                                    if (GOT_DIFFUS) LOC_CV_RHS_J( : ) =  LOC_CV_RHS_J( : )  &
                                        + (1.-FTHETA(:)) * SdevFuns%DETWEI(GI) * DIFF_COEFOLD_DIVDX(:) &
                                        * ( TOLD_ALL(:, CV_NODI) - TOLD_ALL(:, CV_NODJ) )
                                    if (capillary_pressure_activated) LOC_CV_RHS_J( : ) =  LOC_CV_RHS_J( : )  &
                                        - SdevFuns%DETWEI(GI) * CAP_DIFF_COEF_DIVDX(:) & ! capillary pressure stabilization term..
                                        * ( T_ALL(:, CV_NODI) - T_ALL(:, CV_NODJ) )
                                    if (.not.conservative_advection) LOC_CV_RHS_J( : ) =  LOC_CV_RHS_J( : )  &
                                        + FTHETA_T2_J(:) * ( ONE_M_CV_BETA ) * SdevFuns%DETWEI( GI ) * NDOTQNEW(:) * LIMD(:) * T_ALL(:, CV_NODJ) &
                                        - ( ONE_M_CV_BETA) * SdevFuns%DETWEI( GI ) &
                                        * ( FTHETA_T2_J(:) * NDOTQNEW(:) * T_ALL(:, CV_NODJ) * LIMD(:)  &
                                        + ONE_M_FTHETA_T2OLD_J(:) * NDOTQOLD(:) * LIMDOLD(:) * TOLD_ALL(:, CV_NODJ) )
                                endif
                                IF ( GET_GTHETA ) THEN
                                    THETA_GDIFF( :, CV_NODI ) =  THETA_GDIFF( :, CV_NODI ) &
                                        + (1.-FTHETA(:)) * SdevFuns%DETWEI(GI) * DIFF_COEFOLD_DIVDX(:) &
                                        * ( TOLD_ALL(:, CV_NODJ) - TOLD_ALL(:, CV_NODI) ) &
                                        ! Robin bc
                                        + SdevFuns%DETWEI( GI ) * ROBIN2(:)
                                    if(integrate_other_side_and_not_boundary) then
                                        THETA_GDIFF( :, CV_NODJ ) =  THETA_GDIFF( :, CV_NODJ ) &
                                            + (1.-FTHETA(:)) * SdevFuns%DETWEI(GI) * DIFF_COEFOLD_DIVDX(:) &
                                            * ( TOLD_ALL(:, CV_NODI) - TOLD_ALL(:, CV_NODJ) )
                                    endif
                                END IF
                                ! this is for the internal energy equation source term..
                                ! This is to introduce the compressibility term due to expansion and therefore the divergence of the velocity is non-zero
                                !for wells this is not straightforward <= need to CHANGE THIS FOR COMPRESSIBILITY
                                IF ( THERMAL .and. Mdims%npres == 1) THEN
                                    THERM_FTHETA = 1.0
                                    VOL_FRA_FLUID_I = 1.0
                                    VOL_FRA_FLUID_J = 1.0
                                    DO IPRES=1,Mdims%npres
                                        CV_P_PHASE_NODI(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres)=CV_P( 1, IPRES, CV_NODI )
                                        CV_P_PHASE_NODJ(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres)=CV_P( 1, IPRES, CV_NODJ )
                                    END DO
                                    LOC_CV_RHS_I = LOC_CV_RHS_I&
                                        - CV_P_PHASE_NODI * SdevFuns%DETWEI( GI ) * ( &
                                        THERM_FTHETA * NDOTQNEW * LIMT2 &
                                        + ( 1. - THERM_FTHETA ) * NDOTQOLD * LIMT2OLD )*VOL_FRA_FLUID_I
                                    if ( integrate_other_side_and_not_boundary ) then
                                        LOC_CV_RHS_J( : ) = LOC_CV_RHS_J( : ) &
                                            + CV_P_PHASE_NODJ( : ) * SdevFuns%DETWEI( GI ) * ( &
                                            THERM_FTHETA * NDOTQNEW(:) * LIMT2(:) &
                                            + ( 1. - THERM_FTHETA ) * NDOTQOLD(:) * LIMT2OLD(:) )*VOL_FRA_FLUID_J
                                    end if
                                    IF ( GOT_VIS ) THEN
                                        ! stress form of viscosity...
                                        NU_LEV_GI(:, :) =  (1.-THERM_FTHETA) * NUOLDGI_ALL(:,:) + THERM_FTHETA * NUGI_ALL(:,:)
                                        STRESS_IJ_THERM(:,:,:) = 0.0
                                        DO IPHASE=1,Mdims%nphase
                                            CALL CALC_STRESS_TEN( STRESS_IJ_THERM(:,:,IPHASE), ZERO_OR_TWO_THIRDS, Mdims%ndim, &
                                                CVNORMX_ALL(:,GI), NU_LEV_GI(:,IPHASE) * SdevFuns%DETWEI(GI), THERM_U_DIFFUSION(:,:,IPHASE,MAT_NODI), THERM_U_DIFFUSION_VOL(IPHASE,MAT_NODI) )
                                              !UFENX_ALL(1:Mdims%ndim,U_ILOC,GI), UFENX_ALL(1:Mdims%ndim,U_JLOC,GI) * DETWEI(GI), THERM_U_DIFFUSION(:,:,IPHASE,CV_NODI), THERM_U_DIFFUSION_VOL(IPHASE,CV_NODI) )
                                            if ( integrate_other_side_and_not_boundary ) then
                                                STRESS_IJ_THERM_J(:,:,IPHASE) = 0.0
                                                CALL CALC_STRESS_TEN( STRESS_IJ_THERM_J(:,:,IPHASE), ZERO_OR_TWO_THIRDS, Mdims%ndim, &
                                                    CVNORMX_ALL(:,GI), NU_LEV_GI(:,IPHASE) * SdevFuns%DETWEI(GI), THERM_U_DIFFUSION(:,:,IPHASE,MAT_NODJ), THERM_U_DIFFUSION_VOL(IPHASE,MAT_NODJ) )
                                                 !UFENX_ALL(1:Mdims%ndim,U_ILOC,GI), UFENX_ALL(1:Mdims%ndim,U_JLOC,GI) * DETWEI(GI), THERM_U_DIFFUSION(:,:,IPHASE,CV_NODJ), THERM_U_DIFFUSION_VOL(IPHASE,CV_NODJ) )
                                            end if
                                            DO JDIM = 1, Mdims%ndim
                                                DO IDIM = 1, Mdims%ndim
                                                    VECS_STRESS(IDIM,JDIM,IPHASE,CV_NODI) = VECS_STRESS(IDIM,JDIM,IPHASE,CV_NODI) + STRESS_IJ_THERM(IDIM,JDIM,IPHASE)
                                                    VECS_GRAD_U(IDIM,JDIM,IPHASE,CV_NODI) = VECS_GRAD_U(IDIM,JDIM,IPHASE,CV_NODI) + NU_LEV_GI(IDIM,IPHASE) * CVNORMX_ALL(JDIM,GI) * SdevFuns%DETWEI(GI)
                                                    if ( integrate_other_side_and_not_boundary ) then
                                                        VECS_STRESS(IDIM,JDIM,IPHASE,CV_NODJ) = VECS_STRESS(IDIM,JDIM,IPHASE,CV_NODJ) - STRESS_IJ_THERM_J(IDIM,JDIM,IPHASE )
                                                        VECS_GRAD_U(IDIM,JDIM,IPHASE,CV_NODJ) = VECS_GRAD_U(IDIM,JDIM,IPHASE,CV_NODJ) - NU_LEV_GI(IDIM,IPHASE) * CVNORMX_ALL(JDIM,GI) * SdevFuns%DETWEI(GI)
                                                    end if
                                                END DO
                                            END DO
                                        END DO! ENDOF DO IPHASE=1,Mdims%nphase
                                    END IF ! GOT_VIS
                                END IF ! THERMAL
                                call addto(Mmat%CV_RHS,CV_NODI,LOC_CV_RHS_I)
                                call addto(Mmat%CV_RHS,CV_NODJ,LOC_CV_RHS_J)
                            ENDIF Conditional_GETCV_DISC
                        endif ! if(CV_NODJ.ge.CV_NODI) then
                    END IF Conditional_integration
                END DO Loop_GCOUNT
            END DO Loop_CV_ILOC
        END DO Loop_Elements
        IF(GET_GTHETA) THEN
            DO CV_NODI = 1, Mdims%cv_nonods
                THETA_GDIFF(:, CV_NODI) = THETA_GDIFF(:, CV_NODI) / MASS_CV(CV_NODI)
            END DO
        ENDIF
        ! Used for pipe modelling...
        ALLOCATE(MASS_CV_PLUS(Mdims%npres,Mdims%cv_nonods))
        DO CV_NODI = 1, Mdims%cv_nonods
            MASS_CV_PLUS(:,CV_NODI)= MASS_CV(CV_NODI)
        END DO

        IF ( Mdims%npres > 1 ) THEN
            OPT_VEL_UPWIND_COEFS_NEW_CV=0.0 ; N=0.0
            DO ELE = 1, Mdims%totele
!            DO k = 1, size(eles_with_pipe)
!                ELE = eles_with_pipe(k)%ele!Element with pipe
                DO CV_ILOC = 1, Mdims%cv_nloc
                    CV_NODI = ndgln%cv( CV_ILOC + (ELE-1)*Mdims%cv_nloc )
                    MAT_NODI = ndgln%mat( CV_ILOC + (ELE-1)*Mdims%cv_nloc )
                    if (is_porous_media) then
                        RSUM_VEC = 0.0
                        DO IDIM = 1, Mdims%ndim
                            RSUM_VEC = RSUM_VEC + upwnd%adv_coef( IDIM, IDIM, :, MAT_NODI ) / REAL( Mdims%ndim )
                        END DO
                    else
                        RSUM_VEC = 1.0
                    end if
                    OPT_VEL_UPWIND_COEFS_NEW_CV( :, CV_NODI ) = OPT_VEL_UPWIND_COEFS_NEW_CV( :, CV_NODI ) + &
                        RSUM_VEC * MASS_ELE( ELE )
                    N( CV_NODI ) = N( CV_NODI ) + MASS_ELE( ELE )
                END DO
            END DO
            DO CV_NODI = 1, Mdims%cv_nonods
                if (is_porous_media) then
                    SIGMA_INV_APPROX( :, CV_NODI ) = 1.0 / ( OPT_VEL_UPWIND_COEFS_NEW_CV( :, CV_NODI ) / N( CV_NODI ) )
                else if(is_flooding) then
                ! set \sigma for the pipes here      !multi_absorp%Flooding is always of memory type 1
                    SIGMA_INV_APPROX(:, cv_nodi)=1.0/multi_absorp%Flooding%val( 1, 1, :, cv_nodi )!Only has the friction inside the pipes
                else
                    SIGMA_INV_APPROX( :, CV_NODI ) = 1.0
                end if
            END DO
            IF ( PIPES_1D ) THEN
                allocate(MASS_PIPE_FOR_COUP(Mdims%cv_nonods))

                CALL MOD_1D_CT_AND_ADV( state, packed_state, Mdims, ndgln, WIC_T_BC_ALL,WIC_D_BC_ALL, WIC_U_BC_ALL, SUF_T_BC_ALL,SUF_D_BC_ALL,SUF_U_BC_ALL, &
                    getcv_disc, getct, Mmat, Mspars, DT, pipes_aux%MASS_CVFEM2PIPE, pipes_aux%MASS_PIPE2CVFEM, pipes_aux%MASS_CVFEM2PIPE_TRUE, pipes_aux%MASS_PIPE, MASS_PIPE_FOR_COUP, &
                    SIGMA_INV_APPROX, SIGMA_INV_APPROX_NANO, upwnd%adv_coef, eles_with_pipe, THERMAL, cv_beta, bcs_outfluxes, outfluxes)

                if(is_flooding) then
                    DO CV_NODI = 1, Mdims%cv_nonods
                        SIGMA_INV_APPROX(:, CV_NODI)=1.0/multi_absorp%Flooding%val( 1, 1, :, CV_NODI )!Only has the friction inside the pipes
!                        SIGMA_INV_APPROX(1:2, CV_NODI)=1.0/multi_absorp%Flooding%val( 1, 1, 1:2, CV_NODI )!Only has the friction inside the pipes
                    END DO
                endif
                ! Used for pipe modelling...
                DO CV_NODI = 1, Mdims%cv_nonods
                    MASS_CV_PLUS(2:Mdims%npres,CV_NODI) = pipes_aux%MASS_PIPE(CV_NODI)
                END DO
            END IF
            GAMMA_PRES_ABS2 = 0.0
            A_GAMMA_PRES_ABS = 0.0
            if(is_flooding ) then
                bathymetry=>extract_tensor_field(packed_state,"PackedBathymetry")
                depth_of_drain=>extract_scalar_field(state(1),"Drain_depth")
                allocate(R_PEACMAN( Mdims%nphase ) )
            endif 
            DO CV_NODI = 1, Mdims%cv_nonods

                !Only go through the nodes that have a well
                if (PIPE_DIAMETER%val(cv_nodi) < 1e-8) cycle

                ! variables used in the edge approach
                IF ( PIPES_1D ) THEN
                    h = pipes_aux%MASS_PIPE( cv_nodi )/( pi*(0.5*max(PIPE_DIAMETER%val(cv_nodi), 1.0e-10))**2 )
                ELSE
                    h = mass_cv( cv_nodi )**(1.0/Mdims%ndim)
                END IF
                h = max( h, 1.0e-10 )
                h_nano = h
                cc = 0.0
                if ( pipes_aux%MASS_PIPE( cv_nodi )>0.0 ) cc = 1.0 * (1.-FEM_PIPE_CORRECTION) ! to convert Peacement to FEM.
                h = (mass_cv( cv_nodi )/h)**(1.0/(Mdims%ndim-1)) ! This is the lengthscale normal to the wells.
                IF ( is_flooding ) h = mass_cv( cv_nodi )**(1.0/Mdims%ndim)
                rp = 0.14 * h
                rp_nano = 0.14 * h_nano
                Skin = 0.0
                if (GOT_T2) then
                    SAT_FOR_PIPE(:) = MIN( MAX( 0.0, T2_ALL( :, CV_NODI ) ), 1.0 )
                else
                    SAT_FOR_PIPE(:) = MIN( MAX( 0.0, T_ALL( :, CV_NODI ) ), 1.0 )
                end if
                DEN_FOR_PIPE_PHASE(:) =  DEN_ALL( :, CV_NODI )
                DO IPRES = 1, Mdims%npres
                   PRES_FOR_PIPE_PHASE(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres) = FEM_P( 1, IPRES, CV_NODI ) + reservoir_P( IPRES )
                END DO
                PRES_FOR_PIPE_PHASE_FULL(:) = PRES_FOR_PIPE_PHASE(:)
                IF ( is_flooding ) then ! Should really use the manhole diameter here...
                   DO IPRES = 1, Mdims%npres
                      IF(Mmat%CV_pressure) THEN
                         PRES_FOR_PIPE_PHASE(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres) = CV_P( 1, IPRES, CV_NODI ) + reservoir_P( IPRES )
                      ELSE
                         PRES_FOR_PIPE_PHASE(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres) = FEM_P( 1, IPRES, CV_NODI ) + reservoir_P( IPRES )
                      ENDIF
                   END DO
                   PRES_FOR_PIPE_PHASE_FULL(:) = PRES_FOR_PIPE_PHASE(:)
                   DEN_FOR_PIPE_PHASE(1)=DEN_FOR_PIPE_PHASE(3)
                   DEN_FOR_PIPE_PHASE(2)=DEN_FOR_PIPE_PHASE(4)
                   fs_height = DEN_ALL(1,CV_NODI )!for flooding density of the first phase is the height, to make it more readable
                                                !we set a new variable here
                   SIGMA_INV_APPROX( 2, CV_NODI )=SIGMA_INV_APPROX( 4, CV_NODI )
                   SAT_FOR_PIPE(2) = max(pipe_Diameter%val( cv_nodi )-2.*fs_height,0.0)/max( pipe_Diameter%val( cv_nodi ), 1.0e-10 )
                   SAT_FOR_PIPE(1) = 1.-SAT_FOR_PIPE(2)
                   K_PIPES = MIN(K_TOP*SAT_FOR_PIPE(1),1.0)
                   PRES_FOR_PIPE_PHASE_FULL(1) =  gravity_flooding * DEN_FOR_PIPE_PHASE(Mdims%n_in_pres + 1) *( fs_height + depth_of_drain%val( cv_nodi ))
                   PRES_FOR_PIPE_PHASE(1)      =  gravity_flooding * DEN_FOR_PIPE_PHASE(Mdims%n_in_pres + 1) *( fs_height + K_PIPES*depth_of_drain%val( cv_nodi ))

                   Q_PIPES = PRES_FOR_PIPE_PHASE(3) - gravity_flooding * DEN_FOR_PIPE_PHASE(Mdims%n_in_pres + 1) *( DEN_FOR_PIPE_PHASE(1) + depth_of_drain%val( cv_nodi ))
                   CV_PIPE_LENGTH = pipes_aux%MASS_PIPE( cv_nodi )/( pi*(0.5*max(PIPE_DIAMETER%val(cv_nodi),1.0e-10))**2)
                end if

                DO IPHASE = 1, Mdims%nphase
                    IPRES = 1 + INT( (IPHASE-1)/Mdims%n_in_pres )
                    DO JPHASE = 1, Mdims%nphase
                        JPRES = 1 + INT( (JPHASE-1)/Mdims%n_in_pres )
                        ! This is the edge approach
                        ! We do NOT divide by r**2 here because we have not multiplied by r**2 in the pipes_aux%MASS_CVFEM2PIPE matrix (in MOD_1D_CT_AND_ADV)
                        IF ( IPRES /= JPRES ) THEN
                            !Peaceman correction
                            IF ( PRES_FOR_PIPE_PHASE_FULL(IPHASE) > PRES_FOR_PIPE_PHASE_FULL(JPHASE) ) THEN
                                GAMMA_PRES_ABS2( IPHASE, JPHASE, CV_NODI ) = pipes_aux%GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) * &
                                    cc * SAT_FOR_PIPE(IPHASE) * 2.0 * SIGMA_INV_APPROX( IPHASE, CV_NODI ) &
                                    / ( 1.0*(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                                IF ( GOT_NANO ) THEN
                                    GAMMA_PRES_ABS2( IPHASE, JPHASE, CV_NODI ) = GAMMA_PRES_ABS2( IPHASE, JPHASE, CV_NODI ) +pipes_aux%GAMMA_PRES_ABS_NANO( IPHASE, JPHASE, CV_NODI ) * &
                                        cc * SAT_FOR_PIPE(IPHASE) * 2.0 * PI * SIGMA_INV_APPROX_NANO( IPHASE, CV_NODI ) * pipe_length_nano%val( cv_nodi ) &
                                        / ( 1.0*(log( rp_nano / max( 0.5*pipe_Diameter_nano%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                                END IF
                            ELSE
                                GAMMA_PRES_ABS2( IPHASE, JPHASE, CV_NODI ) = pipes_aux%GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) * &
                                    cc * SAT_FOR_PIPE(JPHASE) * 2.0 * SIGMA_INV_APPROX( JPHASE, CV_NODI ) &
                                    / ( 1.0*(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                                IF ( GOT_NANO ) THEN
                                    GAMMA_PRES_ABS2( IPHASE, JPHASE, CV_NODI ) = GAMMA_PRES_ABS2( IPHASE, JPHASE, CV_NODI ) +pipes_aux%GAMMA_PRES_ABS_NANO( IPHASE, JPHASE, CV_NODI ) * &
                                        cc * SAT_FOR_PIPE(JPHASE) * 2.0 * PI * SIGMA_INV_APPROX_NANO( JPHASE, CV_NODI ) * pipe_length_nano%val( cv_nodi ) &
                                        / ( 1.0*(log( rp_nano / max( 0.5*pipe_Diameter_nano%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                                END IF
                            END IF
                        END IF ! IF ( IPRES /= JPRES ) THEN
                    END DO
                END DO

               DO IPHASE = 1, Mdims%nphase
                   DO JPHASE = 1, Mdims%nphase
                       IPRES = 1 + INT( (IPHASE-1)/Mdims%n_in_pres )
                       JPRES = 1 + INT( (JPHASE-1)/Mdims%n_in_pres )
                       IF( IPRES /= JPRES ) THEN
                           A_GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) = - GAMMA_PRES_ABS2( IPHASE, JPHASE, CV_NODI )
                           A_GAMMA_PRES_ABS( IPHASE, IPHASE, CV_NODI ) = A_GAMMA_PRES_ABS( IPHASE, IPHASE, CV_NODI ) + GAMMA_PRES_ABS2( IPHASE, JPHASE, CV_NODI )
                       END IF
                   END DO
               END DO
               if (is_flooding .and. getct)  then
! start again by re-setting to 0.0
                   A_GAMMA_PRES_ABS( :, :, CV_NODI ) = 0.0

                            !Peaceman correction

                   R_PEACMAN=0.0
                   DO IPHASE = 1, Mdims%nphase
                            ISWITCH = MIN( max(IPHASE-2, 0) ,1) ! ISWITCH=0 (for phase 1 and 2) and ISWITCH=1 for phase 3 and 4.
                            JPHASE= (IPHASE+2)*(1-ISWITCH) + (IPHASE-2)*ISWITCH
                            IPRES = 1 + INT( (IPHASE-1)/Mdims%n_in_pres )
                            JPRES = 1 + INT( (JPHASE-1)/Mdims%n_in_pres )
                            IF ( PRES_FOR_PIPE_PHASE_FULL(IPHASE) > PRES_FOR_PIPE_PHASE_FULL(JPHASE) ) THEN
                                R_PEACMAN( IPHASE ) =  pipes_aux%GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) * &
                                    cc * SAT_FOR_PIPE(IPHASE) * 2.0 * SIGMA_INV_APPROX( IPHASE, CV_NODI ) &
                                    / ( 1.0*(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                            ELSE
                                R_PEACMAN( IPHASE ) =  pipes_aux%GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) * &
                                    cc * SAT_FOR_PIPE(JPHASE) * 2.0 * SIGMA_INV_APPROX( JPHASE, CV_NODI ) &
                                    / ( 1.0*(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                            END IF
                   END DO


                   L_surface_pipe = 0.25*pipe_Diameter%val( CV_NODI )
                   l_frac = L_surface_pipe/max(1.0e-10, CV_PIPE_LENGTH)

                   R_PEACMAN = l_frac * R_PEACMAN
!                   R_PEACMAN = R_PEACMAN*1.e+10
!                   R_PEACMAN = R_PEACMAN*1.e+4
!                   R_PEACMAN = R_PEACMAN*1.e+2
!                   R_PEACMAN=0.0


                   A_GAMMA_PRES_ABS( 1, 1, CV_NODI ) = R_PEACMAN( 1 ) * L_surface_pipe
                   A_GAMMA_PRES_ABS( 1, 2, CV_NODI ) = 0.0
                   A_GAMMA_PRES_ABS( 1, 3, CV_NODI ) = - R_PEACMAN( 1 )*L_surface_pipe*K_PIPES/DEN_FOR_PIPE_PHASE(3)
!                   A_GAMMA_PRES_ABS( 1, 3, CV_NODI ) = - R_PEACMAN( 1 )*L_surface_pipe*K_PIPES
                   A_GAMMA_PRES_ABS( 1, 4, CV_NODI ) = 0.0
                   A_GAMMA_PRES_ABS( 2, :, CV_NODI ) = 0.0
!                   A_GAMMA_PRES_ABS( 3, 1, CV_NODI ) = - l_frac*R_PEACMAN( 3  )* DEN_FOR_PIPE_PHASE(3)
!                   A_GAMMA_PRES_ABS( 3, 1, CV_NODI ) = - R_PEACMAN( 3  )* DEN_FOR_PIPE_PHASE(3)
                   A_GAMMA_PRES_ABS( 3, 1, CV_NODI ) = - R_PEACMAN( 3  )
                   A_GAMMA_PRES_ABS( 3, 2, CV_NODI ) = 0.0
!                   A_GAMMA_PRES_ABS( 3, 3, CV_NODI ) = l_frac*R_PEACMAN( 3 )* K_PIPES
                   A_GAMMA_PRES_ABS( 3, 3, CV_NODI ) = R_PEACMAN( 3 )* K_PIPES/DEN_FOR_PIPE_PHASE(3)
                   A_GAMMA_PRES_ABS( 3, 4, CV_NODI ) = 0.0
                   A_GAMMA_PRES_ABS( 4, 1:3, CV_NODI ) = 0.0
!                   A_GAMMA_PRES_ABS( 4, 4, CV_NODI ) = l_frac*R_PEACMAN( 4 )
                   A_GAMMA_PRES_ABS( 4, 4, CV_NODI ) = R_PEACMAN( 4 )/DEN_FOR_PIPE_PHASE(4)
!
! cty rhs...
                   ct_rhs_phase(1)= -R_PEACMAN( 1 ) * L_surface_pipe*gravity_flooding*(-bathymetry%val(1,1,CV_NODI ) + K_PIPES*depth_of_drain%val( CV_NODI ))
                   ct_rhs_phase(2)=0.0
                   ct_rhs_phase(3)= R_PEACMAN( 3 ) *gravity_flooding*(-bathymetry%val(1,1,CV_NODI ) + K_PIPES*depth_of_drain%val( CV_NODI ))
                   ct_rhs_phase(4)=0.0
                   ct_rhs_phase(1:2)=ct_rhs_phase(1:2)*MASS_CV( CV_NODI )
                   ct_rhs_phase(3:4)=ct_rhs_phase(3:4)*MASS_CV_PLUS( 2, CV_NODI )

! Remeber for flooding DEN_ALL_DIVID( 1, CV_NODI )=1.0
! now divid through by the free surface height to be consistent with the rest of the cty equation:
!                   ct_rhs_phase(1) = ct_rhs_phase(1)/DEN_ALL_DIVID( 1, CV_NODI )
!                   A_GAMMA_PRES_ABS( 1, :, CV_NODI ) = A_GAMMA_PRES_ABS( 1, :, CV_NODI )/DEN_ALL_DIVID( 1, CV_NODI )

                   DO IPRES = 1, Mdims%npres
                        call addto(Mmat%CT_RHS, IPRES, cv_nodi, SUM( ct_rhs_phase(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres)) )
                   END DO
               endif
            END DO  ! DO CV_NODI = 1, Mdims%cv_nonods


            if(GETCV_DISC) then
                PIPE_ABS = 0.0
                DO CV_NODI = 1, Mdims%cv_nonods
                    !Only go through the nodes that have a well
                    if (PIPE_DIAMETER%val(cv_nodi) < 1e-8) cycle

                    IF ( PIPES_1D ) THEN
                        h = pipes_aux%MASS_PIPE( cv_nodi )/( pi*(0.5*max(PIPE_DIAMETER%val(cv_nodi),1.0e-10))**2)
                    ELSE
                        h = mass_cv( cv_nodi )**(1.0/Mdims%ndim)
                    END IF
                    h = max( h, 1.0e-10 )
                    h_nano = h
                    cc = 0.0
                    if ( pipes_aux%MASS_PIPE( cv_nodi )>0.0 ) cc = 1.0 * (1.-FEM_PIPE_CORRECTION) ! to convert Peacement to FEM.
                    h = (mass_cv( cv_nodi )/h)**(1.0/(Mdims%ndim-1))  ! This is the lengthscale normal to the wells.
                    rp = 0.14 * h
                    rp_NANO = 0.14 * h_NANO
                    Skin = 0.0

                    SAT_FOR_PIPE(:) = MIN( MAX( 0.0, T_ALL( :, CV_NODI ) ), 1.0 )
                    DO IPRES = 1, Mdims%npres
                       PRES_FOR_PIPE_PHASE(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres) = FEM_P( 1, IPRES, CV_NODI ) + reservoir_P( IPRES )
                       !####For Chris: This option below changes the results, IF that is what we want we need to update the test case####
!                       PRES_FOR_PIPE_PHASE(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres) = CV_P( 1, IPRES, CV_NODI ) + reservoir_P( IPRES )
                    END DO
                    PRES_FOR_PIPE_PHASE_FULL(:) = PRES_FOR_PIPE_PHASE(:)
                    DEN_FOR_PIPE_PHASE(:) =  DEN_ALL( :, CV_NODI )
                    !If got_t2 (mainly for thermal) the equations are also multiplied by the saturation. To keep the code "tidy" we do it here
                    if (GOT_T2)  DEN_FOR_PIPE_PHASE(:) = DEN_FOR_PIPE_PHASE(:) * T2_ALL(:, CV_NODI)

                    IF ( .not. is_flooding ) then
                        DO IPHASE = 1, Mdims%nphase
                            DO JPHASE = 1, Mdims%nphase
                                IPRES = 1 + INT( (IPHASE-1)/Mdims%n_in_pres )
                                JPRES = 1 + INT( (JPHASE-1)/Mdims%n_in_pres )
                                IF ( IPRES /= JPRES ) THEN
                                    DeltaP = PRES_FOR_PIPE_PHASE(IPHASE) - PRES_FOR_PIPE_PHASE(JPHASE)
                                    !DeltaP = FEM_P( 1, IPRES, CV_NODI ) + reservoir_P( ipres ) - ( FEM_P( 1, JPRES, CV_NODI ) + reservoir_P( jpres ) )
                                    ! MEAN_PORE_CV( JPRES, CV_NODI ) is taken out of the following and will be put back only for solving for saturation...
                                    ! We do NOT divide by r**2 here because we have not multiplied by r**2 in the pipes_aux%MASS_CVFEM2PIPE matrix (in MOD_1D_CT_AND_ADV)
                                    IF ( PRES_FOR_PIPE_PHASE_FULL(IPHASE) >= PRES_FOR_PIPE_PHASE_FULL(JPHASE) ) THEN
                                        PIPE_ABS( IPHASE, IPHASE, CV_NODI ) = PIPE_ABS( IPHASE, IPHASE, CV_NODI ) + &
                                            DeltaP * pipes_aux%GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) * DEN_FOR_PIPE_PHASE( IPHASE ) * &
                                            cc * 2.0 * SIGMA_INV_APPROX( IPHASE, CV_NODI ) &
                                            / (log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin)
                                        IF ( GOT_NANO ) THEN
                                            PIPE_ABS( IPHASE, IPHASE, CV_NODI ) = PIPE_ABS( IPHASE, IPHASE, CV_NODI ) +&
                                                DeltaP * pipes_aux%GAMMA_PRES_ABS_NANO( IPHASE, JPHASE, CV_NODI ) * DEN_FOR_PIPE_PHASE( IPHASE ) * &
                                                cc * 2.0 * PI * SIGMA_INV_APPROX_NANO( IPHASE, CV_NODI ) * pipe_length_nano%val( cv_nodi ) &
                                                / ( 1.0 *(log( rp_nano / max( 0.5*pipe_Diameter_nano%val( cv_nodi ), 1.0e-10 ) ) + Skin ) )
                                        END IF
                                    ELSE
                                        PIPE_ABS( IPHASE, JPHASE, CV_NODI ) = PIPE_ABS( IPHASE, JPHASE, CV_NODI ) +&
                                            DeltaP * pipes_aux%GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) * DEN_FOR_PIPE_PHASE( JPHASE ) * &
                                            cc * 2.0 * SIGMA_INV_APPROX( JPHASE, CV_NODI ) &
                                            / ( 1.0 *(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                                        IF ( GOT_NANO ) THEN
                                            PIPE_ABS( IPHASE, JPHASE, CV_NODI ) = PIPE_ABS( IPHASE, JPHASE, CV_NODI ) + &
                                                DeltaP * pipes_aux%GAMMA_PRES_ABS_NANO( IPHASE, JPHASE, CV_NODI ) * DEN_FOR_PIPE_PHASE( JPHASE ) * &
                                                cc * 2.0 * PI * SIGMA_INV_APPROX_NANO( JPHASE, CV_NODI ) * pipe_length_nano%val( cv_nodi ) &
                                                / ( 1.0 *(log( rp_nano / max( 0.5*pipe_Diameter_nano%val( cv_nodi ), 1.0e-10 ) ) + Skin ) )
                                        END IF
                                    END IF
                                END IF
                            END DO
                        END DO
                        if (.not. conservative_advection) then!add extra terms for the non-conservative formulation
                            DO IPHASE = 1, Mdims%nphase
                                DO JPHASE = 1, Mdims%nphase
                                    IPRES = 1 + INT( (IPHASE-1)/Mdims%n_in_pres )
                                    JPRES = 1 + INT( (JPHASE-1)/Mdims%n_in_pres )
                                    IF ( IPRES /= JPRES ) THEN
                                        DeltaP = PRES_FOR_PIPE_PHASE(IPHASE) - PRES_FOR_PIPE_PHASE(JPHASE)
                                        !DeltaP = FEM_P( 1, IPRES, CV_NODI ) + reservoir_P( ipres ) - ( FEM_P( 1, JPRES, CV_NODI ) + reservoir_P( jpres ) )
                                        ! MEAN_PORE_CV( JPRES, CV_NODI ) is taken out of the following and will be put back only for solving for saturation...
                                        ! We do NOT divide by r**2 here because we have not multiplied by r**2 in the pipes_aux%MASS_CVFEM2PIPE matrix (in MOD_1D_CT_AND_ADV)
                                        IF ( PRES_FOR_PIPE_PHASE_FULL(IPHASE) >= PRES_FOR_PIPE_PHASE_FULL(JPHASE)  ) THEN!mass is leaving
                                            PIPE_ABS( IPHASE, IPHASE, CV_NODI ) = PIPE_ABS( IPHASE, IPHASE, CV_NODI ) - &!notice the negative sign
                                                DeltaP * pipes_aux%GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) * DEN_FOR_PIPE_PHASE( IPHASE ) * &
                                                cc * 2.0 * SIGMA_INV_APPROX( IPHASE, CV_NODI ) &
                                                / (log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin)
                                            IF ( GOT_NANO ) PIPE_ABS( IPHASE, IPHASE, CV_NODI ) = PIPE_ABS( IPHASE, IPHASE, CV_NODI ) -&
                                                    DeltaP * pipes_aux%GAMMA_PRES_ABS_NANO( IPHASE, JPHASE, CV_NODI ) * DEN_FOR_PIPE_PHASE( IPHASE ) * &
                                                    cc * 2.0 * PI * SIGMA_INV_APPROX_NANO( IPHASE, CV_NODI ) * pipe_length_nano%val( cv_nodi ) &
                                                    / ( 1.0 *(log( rp_nano / max( 0.5*pipe_Diameter_nano%val( cv_nodi ), 1.0e-10 ) ) + Skin ) )
                                        ELSE
                                           PIPE_ABS( IPHASE, IPHASE, CV_NODI ) = PIPE_ABS( IPHASE, IPHASE, CV_NODI ) - &!notice the negative sign
                                               DeltaP * pipes_aux%GAMMA_PRES_ABS( JPHASE, IPHASE, CV_NODI ) * DEN_FOR_PIPE_PHASE( JPHASE ) * &
                                               cc * 2.0 * SIGMA_INV_APPROX( JPHASE, CV_NODI ) &
                                               / ( 1.0 *(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                                           IF ( GOT_NANO )  PIPE_ABS( IPHASE, IPHASE, CV_NODI ) = PIPE_ABS( IPHASE, IPHASE, CV_NODI ) - &
                                                DeltaP * pipes_aux%GAMMA_PRES_ABS_NANO( IPHASE, JPHASE, CV_NODI ) * DEN_FOR_PIPE_PHASE( JPHASE ) * &
                                                cc * 2.0 * PI * SIGMA_INV_APPROX_NANO( JPHASE, CV_NODI ) * pipe_length_nano%val( cv_nodi ) &
                                                    / ( 1.0 *(log( rp_nano / max( 0.5*pipe_Diameter_nano%val( cv_nodi ), 1.0e-10 ) ) + Skin ) )
                                        END IF
                                    END IF
                                END DO
                            END DO
                        end if
                        !add the heat loss contribution due to diffusion to the nodes with pipes defined, even if it is closed
                        !this might be true only for thermal and porous media
                        if (has_conductivity_pipes) then
                            !Apply onle where wells are closed  !don't like the maxval here...
                            if (maxval(pipes_aux%GAMMA_PRES_ABS( :, :, CV_NODI ))<1d-8) then!REMOVE THIS IF LATER ON IT SHOULD DO NOT HARM...
                                IF ( PIPES_1D ) THEN
                                    h = pipes_aux%MASS_PIPE( cv_nodi )/( pi*(0.5*max(PIPE_DIAMETER%val(cv_nodi),1.0e-10))**2)
                                ELSE
                                    h = mass_cv( cv_nodi )**(1.0/Mdims%ndim)
                                END IF
                                h = max( h, 1.0e-10 )
                                h_nano = h
                                count = min(1,cv_nodi)
                                !Rp is the internal radius of the well
                                rp = max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) - well_thickness%val(count)
                                rp_NANO = rp!0.14 * h_NANO
                                DO IPHASE = 1, Mdims%nphase
                                    DO JPHASE = 1, Mdims%nphase
                                        IPRES = 1 + INT( (IPHASE-1)/Mdims%n_in_pres )
                                        JPRES = 1 + INT( (JPHASE-1)/Mdims%n_in_pres )
                                        IF ( IPRES /= JPRES ) THEN
                                            !we apply Q = (Tin-Tout) * 2 * PI * K * L/(ln(Rout/Rin))
                                            IF ( T_ALL(IPHASE, CV_NODI) >= T_ALL(JPHASE, CV_NODI) ) THEN
!                                            IF (IPRES < JPRES) THEN
                                                PIPE_ABS( IPHASE, IPHASE, CV_NODI ) = PIPE_ABS( IPHASE, IPHASE, CV_NODI ) + &
                                                    conductivity_pipes%val(count) * 2.0 * PI * h &
                                                    / log( max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) / rp )
                                                IF ( GOT_NANO ) PIPE_ABS( IPHASE, IPHASE, CV_NODI ) = PIPE_ABS( IPHASE, IPHASE, CV_NODI ) +&
                                                     conductivity_pipes%val(count) * 2.0 * PI * h_nano &
                                                    / log( max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) / rp_nano )
                                            ELSE
                                                PIPE_ABS( IPHASE, JPHASE, CV_NODI ) = PIPE_ABS( IPHASE, JPHASE, CV_NODI ) - &
                                                    conductivity_pipes%val(count) * 2.0 * PI * h &
                                                    / log( max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) / rp )
                                                IF ( GOT_NANO ) PIPE_ABS( IPHASE, JPHASE, CV_NODI ) = PIPE_ABS( IPHASE, JPHASE, CV_NODI ) - &
                                                    conductivity_pipes%val(count) * 2.0 * PI * h_nano &
                                                    / log( max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) / rp_nano )
                                            END IF
                                        END IF
                                    END DO
                                END DO
                            end if
                        end if
                    else ! flooding
                        ! Should really use the manhole diameter here...
                        DO IPRES = 1, Mdims%npres
!                           PRES_FOR_PIPE_PHASE(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres) = FEM_P( 1, IPRES, CV_NODI ) + reservoir_P( IPRES )
                           PRES_FOR_PIPE_PHASE(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres) = CV_P( 1, IPRES, CV_NODI ) + reservoir_P( IPRES )
                        END DO
                        PRES_FOR_PIPE_PHASE_FULL(:) = PRES_FOR_PIPE_PHASE(:)
                        DEN_FOR_PIPE_PHASE(1)=DEN_FOR_PIPE_PHASE(3)
                        DEN_FOR_PIPE_PHASE(2)=DEN_FOR_PIPE_PHASE(4)
                        fs_height = DEN_ALL(1,CV_NODI )!for flooding density of the first phase is the height, to make it more readable
                                                       !we set a new variable here
                        SIGMA_INV_APPROX( 2, CV_NODI )=SIGMA_INV_APPROX( 4, CV_NODI )
                        SAT_FOR_PIPE(2) = max(pipe_Diameter%val( cv_nodi )-2.*fs_height,0.0)/max( pipe_Diameter%val( cv_nodi ), 1.0e-10 )
                        SAT_FOR_PIPE(1) = 1.-SAT_FOR_PIPE(2)
                        K_PIPES = MIN(K_TOP*SAT_FOR_PIPE(1),1.0)
                        PRES_FOR_PIPE_PHASE_FULL(1) =  gravity_flooding * DEN_FOR_PIPE_PHASE(Mdims%n_in_pres + 1) *( fs_height + depth_of_drain%val( cv_nodi ))
                        PRES_FOR_PIPE_PHASE(1)      =  gravity_flooding * DEN_FOR_PIPE_PHASE(Mdims%n_in_pres + 1) *( fs_height + K_PIPES*depth_of_drain%val( cv_nodi ))

                        Q_PIPES = PRES_FOR_PIPE_PHASE(2) - gravity_flooding * DEN_FOR_PIPE_PHASE(Mdims%n_in_pres + 1) *( DEN_FOR_PIPE_PHASE(1) + depth_of_drain%val( cv_nodi ))
                        CV_PIPE_LENGTH = pipes_aux%MASS_PIPE( cv_nodi )/( pi*(0.5*max(PIPE_DIAMETER%val(cv_nodi),1.0e-10))**2)


                        CONT_PIPE_ABS=0.0
                        LOC_CV_RHS_I =0.0

                        IPHASE=1
                        JPHASE=3

                        IPRES = 1 + INT( (IPHASE-1)/Mdims%n_in_pres )
                        JPRES = 1 + INT( (JPHASE-1)/Mdims%n_in_pres )

!                        DeltaP      = PRES_FOR_PIPE_PHASE(IPRES) - PRES_FOR_PIPE_PHASE(JPRES)
                        DeltaP      = PRES_FOR_PIPE_PHASE(IPHASE) - PRES_FOR_PIPE_PHASE(JPHASE)
                        ! MEAN_PORE_CV( JPRES, CV_NODI ) is taken out of the following and will be put back only for solving for saturation...
                        ! We do NOT divide by r**2 here because we have not multiplied by r**2 in the pipes_aux%MASS_CVFEM2PIPE matrix (in MOD_1D_CT_AND_ADV)
                        IF ( PRES_FOR_PIPE_PHASE_FULL(IPHASE) - PRES_FOR_PIPE_PHASE_FULL(JPHASE) >= 0.0 ) THEN
                            CONT_PIPE_ABS( IPHASE, JPHASE ) = &
                                pipes_aux%GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) * DEN_FOR_PIPE_PHASE( IPHASE ) * &
                                cc * 2.0 * SIGMA_INV_APPROX( IPHASE, CV_NODI ) &
                                / ( 1.0 *(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                        ELSE
                            CONT_PIPE_ABS( IPHASE, JPHASE ) = &
                                pipes_aux%GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) * DEN_FOR_PIPE_PHASE( JPHASE ) * &
                                cc * 2.0 * SIGMA_INV_APPROX( JPHASE, CV_NODI ) &
                                / ( 1.0 *(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                        END IF

                        IPHASE=2
                        JPHASE=4

                        IPRES = 1 + INT( (IPHASE-1)/Mdims%n_in_pres )
                        JPRES = 1 + INT( (JPHASE-1)/Mdims%n_in_pres )

                        ! MEAN_PORE_CV( JPRES, CV_NODI ) is taken out of the following and will be put back only for solving for saturation...
                        ! We do NOT divide by r**2 here because we have not multiplied by r**2 in the pipes_aux%MASS_CVFEM2PIPE matrix (in MOD_1D_CT_AND_ADV)
                        IF ( PRES_FOR_PIPE_PHASE_FULL(IPHASE) - PRES_FOR_PIPE_PHASE_FULL(JPHASE) >= 0.0 ) THEN
                            CONT_PIPE_ABS( IPHASE, JPHASE ) = &
                                pipes_aux%GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) * DEN_FOR_PIPE_PHASE( IPHASE ) * &
                                cc * 2.0 * SIGMA_INV_APPROX( IPHASE, CV_NODI ) &
                                / ( 1.0 *(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                        ELSE
                            CONT_PIPE_ABS( IPHASE, JPHASE ) = &
                                pipes_aux%GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) * DEN_FOR_PIPE_PHASE( JPHASE ) * &
                                cc * 2.0 * SIGMA_INV_APPROX( JPHASE, CV_NODI ) &
                                / ( 1.0 *(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                        END IF

                        PIPE_ABS( :, :, CV_NODI ) = 0.0
                        L_surface_pipe = 0.25*pipe_Diameter%val( CV_NODI )
                        l_frac = L_surface_pipe/max(1.0e-10, CV_PIPE_LENGTH)
                        if(implicit_fs) then
                           RDUM = CONT_PIPE_ABS( 1, 3 ) * L_surface_pipe * gravity_flooding
                           PIPE_ABS( 1, 1, CV_NODI ) = RDUM
                        ! Subtract out this contribution so we can treat implicitly in saturation.
                           LOC_CV_RHS_I(1) = RDUM * fs_height
                        endif

                        IF(CV_P(1, 2,CV_NODI)>0.0) THEN ! Air can move to surface or disapear
                            PIPE_ABS( 4, 4, CV_NODI ) = l_frac*CONT_PIPE_ABS( 2, 4 ) * CV_P(1,2,CV_NODI)
                        ELSE
                            LOC_CV_RHS_I(4) = LOC_CV_RHS_I(4) - l_frac*CONT_PIPE_ABS( 2, 4 ) * CV_P(1, 2,CV_NODI)*SAT_FOR_PIPE(2)
                        ENDIF

                        DeltaP = PRES_FOR_PIPE_PHASE(3) - PRES_FOR_PIPE_PHASE(1)
                        IF(PRES_FOR_PIPE_PHASE_FULL(3) - PRES_FOR_PIPE_PHASE_FULL(1) >= 0.0 ) THEN
                            RDUM2 = CONT_PIPE_ABS( 1, 3 ) * DeltaP * L_surface_pipe /DEN_FOR_PIPE_PHASE( 3 )
                            if(implicit_fs) then
                                RDUM = CONT_PIPE_ABS( 1, 3 ) * gravity_flooding * DEN_FOR_PIPE_PHASE( 3 ) * l_frac ! add and subtract this term to make more implicit
                                PIPE_ABS( 3, 1, CV_NODI ) = PIPE_ABS( 3, 1, CV_NODI ) + RDUM
                                LOC_CV_RHS_I(3) = LOC_CV_RHS_I(3) + RDUM * fs_height
                                PIPE_ABS( 1, 3, CV_NODI ) = PIPE_ABS( 1, 3, CV_NODI ) + RDUM2
                            else
                                LOC_CV_RHS_I(1) = LOC_CV_RHS_I(1) - RDUM2 * SAT_FOR_PIPE(3)
                            endif
                            PIPE_ABS( 3, 3, CV_NODI ) = PIPE_ABS( 3, 3, CV_NODI ) + CONT_PIPE_ABS( 1, 3 ) * DeltaP * l_frac
                        ELSE
                            LOC_CV_RHS_I(1) = LOC_CV_RHS_I(1) + (CONT_PIPE_ABS( 1, 3 ) * DeltaP * L_surface_pipe /DEN_FOR_PIPE_PHASE( 3 )) * SAT_FOR_PIPE(1)
                            LOC_CV_RHS_I(3) = LOC_CV_RHS_I(3) - (CONT_PIPE_ABS( 1, 3 ) * DeltaP ) * l_frac * SAT_FOR_PIPE(1)
                        ENDIF


                        LOC_CV_RHS_I(:) = LOC_CV_RHS_I(:) * pipes_aux%MASS_PIPE( CV_NODI )

                        call addto(Mmat%CV_RHS, cv_nodi, LOC_CV_RHS_I)


                    end if ! IF ( .not. is_flooding ) then else

                END DO ! DO CV_NODI = 1, Mdims%cv_nonods
            endif ! if(GETCV_DISC) then

            IF ( GETCT ) THEN
                INV_B = DT * PIPE_ABS
                if (.not.is_flooding) then
                    DO IPHASE = 1, Mdims%nphase
                        INV_B( IPHASE, IPHASE, : ) = INV_B( IPHASE, IPHASE, : ) + DEN_ALL( IPHASE, : )
                      !IPRES = 1 + INT( (IPHASE-1)/Mdims%n_in_pres )
                      !MEAN_PORE_CV( IPRES, : ) * DEN_ALL( IPHASE, : )
                      !DEN_ALL( IPHASE, : )
                    END DO
                else
                    !Different from reservoir to pipes, as density mean something different, so we use the density from the pipes
                    DO IPHASE = 1, Mdims%n_in_pres!as it actually defines the density in the reservoir
                        INV_B( IPHASE, IPHASE, : ) = INV_B( IPHASE, IPHASE, : ) + DEN_ALL( Mdims%n_in_pres + IPHASE, : )
                    END DO
                    DO IPHASE = Mdims%n_in_pres + 1, Mdims%nphase
                        INV_B( IPHASE, IPHASE, : ) = INV_B( IPHASE, IPHASE, : ) + DEN_ALL( IPHASE, : )
                    END DO
                end if
                DO CV_NODI = 1, Mdims%cv_nonods

                    CALL INVERT( INV_B( :, :, CV_NODI ) )
                END DO
                IF ( MULTB_BY_POROSITY ) THEN
                    DO IPHASE = 1, Mdims%nphase
                        IPRES = 1 + INT( (IPHASE-1)/Mdims%n_in_pres )
                        INV_B( IPHASE, IPHASE, : ) = INV_B( IPHASE, IPHASE, : ) * MEAN_PORE_CV( IPRES, : )
                    END DO
                END IF
            ENDIF ! ENDOF IF ( GETCT ) THEN
        END IF ! IF ( Mdims%npres > 1 ) THEN

        Conditional_GETCV_DISC2: IF( GETCV_DISC ) THEN ! Obtain the CV discretised advection/diffusion equations
            Loop_CVNODI2: DO CV_NODI = 1, Mdims%cv_nonods ! Put onto the diagonal of the matrix
                LOC_CV_RHS_I=0.0
                DO IPRES = 1, Mdims%npres
                    R_PHASE(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres) = MEAN_PORE_CV( IPRES, CV_NODI ) * MASS_CV_PLUS( IPRES, CV_NODI ) / DT
                    CV_P_PHASE_NODI(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres) = CV_P( 1, IPRES, CV_NODI ) + reservoir_P( IPRES )
                END DO
                IF ( THERMAL .and. Mdims%npres == 1) THEN
                    IF ( GOT_VIS ) THEN
                        DO IPHASE = 1, Mdims%nphase
                            LOC_CV_RHS_I(IPHASE)=LOC_CV_RHS_I(IPHASE)  &
                                + SUM( VECS_STRESS(:,:,IPHASE,CV_NODI)*VECS_GRAD_U(:,:,IPHASE,CV_NODI)  )/MASS_CV(CV_NODI)
                        END DO
                    END IF
                    LOC_CV_RHS_I(:)=LOC_CV_RHS_I(:) &
                        - CV_P_PHASE_NODI(:) * ( MASS_CV( CV_NODI ) / DT ) * ( T2_ALL( :, CV_NODI ) - T2OLD_ALL( :, CV_NODI ) )
                END IF

                IF ( GOT_T2 ) THEN
                    LOC_CV_RHS_I(:)=LOC_CV_RHS_I(:)  &
                        + MASS_CV(CV_NODI) * SOURCT_ALL( :, CV_NODI )
                    if (thermal .and. is_porous_media) then
                        !In this case for the time-integration term the effective rho Cp is a combination of the porous media
                        ! and the fluids. Here we add the porous media contribution
                        DO IPHASE = 1,Mdims%n_in_pres
                            call addto(Mmat%petsc_ACV,iphase,iphase,&
                                cv_nodi, cv_nodi,&
                                + porous_heat_coef( IPHASE, CV_NODI ) * T2_ALL( IPHASE, CV_NODI ) &
                                * R_PHASE(IPHASE) * (1-MEAN_PORE_CV( 1, CV_NODI ))/MEAN_PORE_CV( 1, CV_NODI ))
                                !R_PHASE includes the porosity. Since in this case we are interested in what is NOT porous
                                    !we divide to remove that term and multiply by the correct term (1-porosity)
                            LOC_CV_RHS_I(iphase)=LOC_CV_RHS_I(iphase)  &
                                + (CV_BETA * porous_heat_coef( iphase, CV_NODI ) * T2OLD_ALL( iphase, CV_NODI ) &
                                + (ONE_M_CV_BETA) * porous_heat_coef( iphase, CV_NODI ) * T2_ALL( iphase, CV_NODI ) ) &
                                * R_PHASE(iphase) * TOLD_ALL( iphase, CV_NODI )* (1-MEAN_PORE_CV( 1, CV_NODI ))/MEAN_PORE_CV( 1, CV_NODI )
                        END DO
                    end if
                    DO IPHASE = 1,Mdims%nphase
                        call addto(Mmat%petsc_ACV,iphase,iphase,&
                            cv_nodi, cv_nodi,&
                            + DEN_ALL( IPHASE, CV_NODI ) * T2_ALL( IPHASE, CV_NODI ) &
                            * R_PHASE(IPHASE))
                    END DO
                    LOC_CV_RHS_I(:)=LOC_CV_RHS_I(:)  &
                        + (CV_BETA * DENOLD_ALL( :, CV_NODI ) * T2OLD_ALL( :, CV_NODI ) &
                        + (ONE_M_CV_BETA) * DEN_ALL( :, CV_NODI ) * T2_ALL( :, CV_NODI ) ) &
                        * R_PHASE(:) * TOLD_ALL( :, CV_NODI )
                ELSE
                    LOC_CV_RHS_I(:)=LOC_CV_RHS_I(:)  &
                        + MASS_CV( CV_NODI ) * SOURCT_ALL( :, CV_NODI )
                    DO IPHASE = 1,Mdims%nphase
                        call addto(Mmat%petsc_ACV,iphase,iphase,&
                            cv_nodi, cv_nodi,&
                            + DEN_ALL( IPHASE, CV_NODI )  &
                            * R_PHASE(IPHASE) )
                    END DO
                    LOC_CV_RHS_I(:)=LOC_CV_RHS_I(:)  &
                        + ( CV_BETA * DENOLD_ALL( :, CV_NODI ) &
                        + (ONE_M_CV_BETA) * DEN_ALL( :, CV_NODI ) ) &
                        * R_PHASE(:) * TOLD_ALL( :, CV_NODI )
                END IF
                IF( Mdims%npres > 1 .and. explicit_pipes ) THEN ! this is not safe...
                    DO JPHASE = 1, Mdims%nphase
                        LOC_CV_RHS_I(:)=LOC_CV_RHS_I(:)  &
                            - pipes_aux%MASS_PIPE( CV_NODI ) * A_GAMMA_PRES_ABS(:,JPHASE,CV_NODI )*CV_P_PHASE_NODI(JPHASE )
                    END DO
                END IF

                Conditional_GETMAT2: IF ( GETMAT ) THEN
                    do jphase = 1, Mdims%nphase
                        do iphase=1,Mdims%nphase
                            IF ( Mdims%npres > 1 .AND. .NOT.EXPLICIT_PIPES ) THEN
                                call addto(Mmat%petsc_ACV,iphase,jphase, &
                                    cv_nodi, cv_nodi, &
                                    MASS_PIPE_FOR_COUP( CV_NODI ) * PIPE_ABS( iphase, jphase, CV_NODI ))
                            END IF
                            if ( have_absorption ) then
                               call addto(Mmat%petsc_ACV,iphase,jphase, &
                                   cv_nodi, cv_nodi, &
                                   MASS_CV( CV_NODI ) * ABSORBT_ALL( iphase, jphase, CV_NODI ))
                            end if
                        end do
                    end do
                END IF Conditional_GETMAT2
                call addto(Mmat%CV_RHS,CV_NODI,LOC_CV_RHS_I)
            END DO Loop_CVNODI2
        END IF Conditional_GETCV_DISC2

        IF ( GETCT ) THEN
            W_SUM_ONE1 = 1.0 !If == 1.0 applies constraint to T
            W_SUM_ONE2 = 0.0 !If == 1.0 applies constraint to TOLD
            DIAG_SCALE_PRES = 0.0
            allocate(DIAG_SCALE_PRES_phase(Mdims%nphase))
            DIAG_SCALE_PRES_COUP=0.0
            DO CV_NODI = 1, Mdims%cv_nonods
                ct_rhs_phase=0.0 ; DIAG_SCALE_PRES_phase=0.0
                DO IPRES=1,Mdims%npres
                    R_PRES(IPRES) = MASS_CV_PLUS( IPRES, CV_NODI ) * MEAN_PORE_CV( IPRES, CV_NODI ) / DT
                    ! Add constraint to force sum of volume fracts to be unity...
                       ! W_SUM_ONE==1 applies the constraint
                       ! W_SUM_ONE==0 does NOT apply the constraint
                    call addto(Mmat%CT_RHS,IPRES,cv_nodi,&
                        - ( W_SUM_ONE1 - W_SUM_ONE2 ) * R_PRES(IPRES))
                    R_PHASE(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres)=R_PRES(IPRES)
                    MEAN_PORE_CV_PHASE(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres) = MEAN_PORE_CV( IPRES, CV_NODI )
                    CV_P_PHASE_NODI(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres)=CV_P( 1, IPRES, CV_NODI )
                END DO
                ct_rhs_phase(:)=ct_rhs_phase(:) &
                    - R_PHASE(:) * ( &
                    + (1.0-W_SUM_ONE1) * T_ALL( :, CV_NODI ) - (1.0-W_SUM_ONE2) * TOLD_ALL( :, CV_NODI ) &
                    + ( TOLD_ALL( :, CV_NODI ) * ( DEN_ALL( :, CV_NODI ) - DENOLD_ALL( :, CV_NODI ) ) &
                    - DERIV( :, CV_NODI ) * CV_P_PHASE_NODI( : ) * T_ALL_KEEP( :, CV_NODI ) ) / DEN_ALL_DIVID(:, CV_NODI ) )
                DIAG_SCALE_PRES_phase( : ) = DIAG_SCALE_PRES_phase( : ) &
                    +  MEAN_PORE_CV_PHASE(:) * T_ALL_KEEP( :, CV_NODI ) * DERIV( :, CV_NODI ) &
                    / ( DT * DEN_ALL_DIVID(:, CV_NODI) )
                ct_rhs_phase(:)=ct_rhs_phase(:)  &
                    + MASS_CV( CV_NODI ) * SOURCT_ALL( :, CV_NODI ) / DEN_ALL_DIVID(:, CV_NODI)
                IF ( HAVE_ABSORPTION ) THEN
                   DO JPHASE = 1, Mdims%nphase
                      ct_rhs_phase(:)=ct_rhs_phase(:)  &
                         - MASS_CV( CV_NODI ) * ABSORBT_ALL( :, JPHASE, CV_NODI ) * T_ALL( JPHASE, CV_NODI ) / DEN_ALL_DIVID(:, CV_NODI)
                   END DO
                END IF
                ! scaling coefficient...
                IF ( Mdims%npres > 1 .AND. explicit_pipes2 ) THEN
                    DO IPRES=1,Mdims%npres
                        DO JPRES=1,Mdims%npres
                            DO jphase=1+(jpres-1)*Mdims%n_in_pres, jpres*Mdims%n_in_pres
                                ! dont divid the pipe to reservoir mass exchange term by density...
                                DIAG_SCALE_PRES_COUP(IPRES,JPRES, cv_nodi) = DIAG_SCALE_PRES_COUP(IPRES,JPRES, cv_nodi)  &
                                    !+ sum( A_GAMMA_PRES_ABS( 1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres    , 1+(jpres-1)*Mdims%n_in_pres:jpres*Mdims%n_in_pres, CV_NODI ) )
                                    + sum( A_GAMMA_PRES_ABS(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres,JPHASE, CV_NODI )  )
                            !   + sum( A_GAMMA_PRES_ABS(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres,JPHASE, CV_NODI )  &
                            !   / DEN_ALL( 1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres, CV_NODI ) )
                            end do
                        END DO
                    END DO
                END IF
                IF ( Mdims%npres > 1 .AND. .NOT.explicit_pipes2 ) THEN
                    DIAG_SCALE_PRES_phase( : ) = DIAG_SCALE_PRES_phase( : ) * DEN_ALL( :, CV_NODI )
                    CT_RHS_PHASE( : ) = CT_RHS_PHASE( : ) * DEN_ALL( :, CV_NODI )
                    CT_RHS_PHASE( : ) = MATMUL( INV_B( :, :, CV_NODI) , CT_RHS_PHASE( : ) )
                END IF
                DO IPRES = 1, Mdims%npres
                    if ( Mdims%npres > 1 ) then
                        call addto(Mmat%CT_RHS, IPRES, cv_nodi, SUM( ct_rhs_phase(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres)) &
                            - MASS_CV_PLUS(IPRES,CV_NODI) * SUM(DIAG_SCALE_PRES_COUP(IPRES,:, cv_nodi) * RESERVOIR_P( : )) )
                    else
                        call addto(Mmat%CT_RHS, IPRES, cv_nodi, SUM( ct_rhs_phase(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres)) )
                    end if
                    DIAG_SCALE_PRES( IPRES,CV_NODI ) = DIAG_SCALE_PRES( IPRES,CV_NODI ) &
                        + sum( DIAG_SCALE_PRES_phase(1+(ipres-1)*Mdims%n_in_pres:ipres*Mdims%n_in_pres))
                END DO
            END DO  ! endof DO CV_NODI = 1, Mdims%cv_nonods
            deallocate(DIAG_SCALE_PRES_phase)
           !deallocate(R_PRES,R_PHASE,MEAN_PORE_CV_PHASE)
            if(is_porous_media .and. present(calculate_mass_delta) .and. present(outfluxes) ) then
                !Calculate final outfluxes and mass balance in the domain
                call mass_conservation_check_and_outfluxes(calculate_mass_delta, outfluxes, 2)
            end if
        END IF



        ! Deallocating temporary working arrays
        IF(GETCT) THEN
            DEALLOCATE( JCOUNT_KLOC )
        ENDIF
        DEALLOCATE( CV_OTHER_LOC )
        DEALLOCATE( U_OTHER_LOC )
        DEALLOCATE( MAT_OTHER_LOC )
        DEALLOCATE( X_SHARE )
        DEALLOCATE( CV_SLOC2LOC )
        DEALLOCATE( U_SLOC2LOC )
        DEALLOCATE( MASS_CV )
        DEALLOCATE( XC_CV_ALL )
        DEALLOCATE( FACE_ELE )
        DEALLOCATE(TUPWIND_MAT_ALL)
        DEALLOCATE(TOLDUPWIND_MAT_ALL)
        DEALLOCATE(DENUPWIND_MAT_ALL)
        DEALLOCATE(DENOLDUPWIND_MAT_ALL)
        DEALLOCATE(T2UPWIND_MAT_ALL)
        DEALLOCATE(T2OLDUPWIND_MAT_ALL)
        DEALLOCATE( DUMMY_ZERO_NDIM_NDIM_NPHASE )
        call deallocate_multi_dev_shape_funs(SdevFuns)
        call deallocate(tracer_BCs)
        call deallocate(tracer_BCs_robin2)
        call deallocate(density_BCs)
        call deallocate(velocity_BCs)
        if(got_free_surf .or. is_porous_media .or. Mmat%CV_pressure) call deallocate(pressure_BCs)
        if (present(saturation)) then
            call deallocate(saturation_BCs)
            call deallocate(saturation_BCs_robin2)
        end if
        deallocate(INCOMEOLD, NDOTQOLD, NUOLDGI_ALL)
        DEALLOCATE( FEMT_ALL, FEMTOLD_ALL, FEMDEN_ALL, FEMDENOLD_ALL, FEMT2_ALL, FEMT2OLD_ALL)
        nullify(FEMT_ALL); nullify(FEMTOLD_ALL);
        !if (present(T_input)) then!<==TEMPORARY
        !      deallocate( T_ALL_TARGET, TOLD_ALL_TARGET, FEMT_ALL_TARGET, FEMTOLD_ALL_TARGET)
        !end if
        !      if ( Field_selector == 1 ) then ! Temperature

         if (tracer%name == "PackedTemperature" )  then
            deallocate( suf_t_bc, suf_t_bc_rob1, suf_t_bc_rob2)
        end if
        if (capillary_pressure_activated) deallocate(CAP_DIFFUSION)
        ewrite(3,*) 'Leaving CV_ASSEMB'
        if (is_flooding) deallocate(DEN_ALL_DIVID)
        if (allocated(bcs_outfluxes)) deallocate(bcs_outfluxes)
        RETURN
    contains

    SUBROUTINE APPLY_ENO_2_T(LIMF, T_ALL,TOLD_ALL, FEMT_ALL,FEMTOLD_ALL, INCOME,INCOMEOLD, IGOT_T_PACK, &
        CV_NODI, CV_NODJ, X_NODI, X_NODJ, CV_ILOC, CV_JLOC, &
        ELE, CV_NONODS, NDIM, NPHASE,  &
        CV_NLOC,TOTELE, X_NDGLN, CV_NDGLN,  &
        X_ALL,FACE_ELE,NFACE,BETWEEN_ELEMENTS, SCVFEN, SCVFENx, GI, INV_JAC, &
        UGI, on_domain_boundary )
        ! Apply ENO method to T and TOLD_ALL and put the limited values in LIMF
        IMPLICIT NONE
        REAL INFINY
        INTEGER ngi_one
        PARAMETER(INFINY=1.E+20)
        PARAMETER(ngi_one=1)
        ! If ENO_ALL_THREE use all 3 to find TGI ENO value...
        ! else use the upwind value and current element value.
        LOGICAL, PARAMETER :: ENO_ALL_THREE =  .false.   !.TRUE.
        ! If FOR_DG_ONLY_BETWEEN_ELE use fem INSIDE element for DG
        LOGICAL, PARAMETER :: FOR_DG_ONLY_BETWEEN_ELE = .TRUE.
        ! Use ENO only where there is an oscillation as it can be a bit dissipative.
        LOGICAL, PARAMETER :: ENO_ONLY_WHERE_OSCILLATE = .TRUE.
        ! Simple ENO - somewhat dispersive
        LOGICAL, PARAMETER :: SIMPLE_ENO = .FALSE.
        ! RELAX_DOWN_WIND_2_CURRENT_ELE=1.0 (use upwind only), RELAX_DOWN_WIND_2_CURRENT_ELE=0 (use downwind in ENO)
        REAL, PARAMETER :: RELAX_DOWN_WIND_2_CURRENT_ELE = 1.0

        INTEGER, intent(in) :: ELE,CV_NONODS,CV_NLOC,TOTELE,NDIM,NPHASE,NFACE,GI
        INTEGER, intent(in) :: CV_NODI, CV_NODJ, X_NODI, X_NODJ, CV_ILOC, CV_JLOC
        REAL, dimension(:), intent(inout) :: LIMF
        REAL, dimension(:,:), intent(in) :: T_ALL,TOLD_ALL,FEMT_ALL,FEMTOLD_ALL
        REAL, dimension(:), intent(in) :: INCOME,INCOMEOLD
        INTEGER, dimension(:), intent(in) :: X_NDGLN,CV_NDGLN
        REAL, dimension(:,:), intent(in) :: X_ALL
        REAL, dimension(:,:), intent(in) :: SCVFEN, UGI
        INTEGER, dimension(:,:), intent(in) :: FACE_ELE
        LOGICAL, dimension(:,:), intent(in) :: IGOT_T_PACK
        LOGICAL, intent(in) :: BETWEEN_ELEMENTS, on_domain_boundary
        REAL, dimension(:,:,:), intent(in) :: INV_JAC, SCVFENx
        !
        !     Local variables...
        REAL, dimension(NDIM+1):: LOCCORDS
        REAL :: MINCOR,MINCORK
        !The dimension of the variables below should be NDIM, however, due to cross products
        !we need three dimensions
        REAL :: XVEC(NDIM),XPT(NDIM),XPT_GI(NDIM), UGI_normalised(ndim,nphase)
        REAL :: n(cv_nloc,ngi_one), nlx(cv_nloc,ngi_one), nly(cv_nloc,ngi_one), nlz(cv_nloc,ngi_one)
        REAL :: l1(ngi_one), l2(ngi_one), l3(ngi_one), l4(ngi_one), weight(ngi_one)
        REAL :: TGI_ELE(NPHASE*2),TGI_IN(NPHASE*2),TGI_OUT(NPHASE*2),W(NPHASE*2), TGI_DER_ELE(NDIM,NPHASE*2)
        REAL :: TGI(NPHASE*2),TUP(NPHASE*2), TGI_NEI(NPHASE*2), TGI_NEI2(NPHASE*2), INCOME_BOTH(2*NPHASE)
        REAL :: TGI_ELE_B(NPHASE*2), TGI_ELE_C(NPHASE*2)
        REAL :: ENO_ELE_MATWEI_IN(CV_NLOC),ENO_ELE_MATWEI_OUT(CV_NLOC),ENO_ELE_MATWEI(CV_NLOC,2)
        REAL :: RUP_WIN,MIN_TGI,MAX_TGI,dx, min_val, max_val
        LOGICAL :: QUADRATIC_ELEMENT,IS_CORNER_NOD_I,IS_CORNER_NOD_J,DISTCONTINUOUS_METHOD!, GOT_AN_OSC
        INTEGER :: ENO_ELE_NEI(2),LOCNODS(NDIM+1)
        INTEGER :: ELE2,SELE2,ELEWIC,ENO_ELE_NEI_IN,ENO_ELE_NEI_OUT,IPHASE2,CV_KLOC,IUP_DOWN
        INTEGER :: IFACE,NPHASE2,IPT, cv_nodk, cv_nodk_IN, cv_nodk_OUT, X_KNOD, I_OLD_NEW, IPHASE

        ! nothing to do (needed to apply bcs etc)...
        if ( on_domain_boundary ) return

        QUADRATIC_ELEMENT=( (NDIM==2).AND.(CV_NLOC==6) ) .OR. ( (NDIM==3).AND.(CV_NLOC==10) )
        NPHASE2=NPHASE*2
        DISTCONTINUOUS_METHOD=( CV_NONODS == TOTELE * CV_NLOC )


        if ( .true. ) then

            TGI_ELE=0.0
            do cv_kloc=1,cv_nloc
                cv_nodk = CV_NDGLN((ELE-1)*CV_NLOC + cv_kloc)
                TGI_ELE(1:NPHASE)=TGI_ELE(1:NPHASE) + SCVFEN(cv_Kloc,gi)*FEMT_ALL(:,cv_nodk)
                TGI_ELE(1+NPHASE:2*NPHASE)=TGI_ELE(1+NPHASE:2*NPHASE) + SCVFEN(cv_Kloc,gi)*FEMTOLD_ALL(:,cv_nodk)
            end do

        else

            TGI_ELE=0.0
            do cv_kloc=1,cv_nloc
                cv_nodk = CV_NDGLN((ELE-1)*CV_NLOC + cv_kloc)
                TGI_ELE(1:NPHASE)=TGI_ELE(1:NPHASE) + SCVFEN(cv_Kloc,gi)*T_ALL(:,cv_nodk)
                TGI_ELE(1+NPHASE:2*NPHASE)=TGI_ELE(1+NPHASE:2*NPHASE) + SCVFEN(cv_Kloc,gi)*TOLD_ALL(:,cv_nodk)
            end do

        end if


        if ( simple_eno ) then


            TGI_ELE=0.0
            do cv_kloc=1,cv_nloc
                cv_nodk = CV_NDGLN((ELE-1)*CV_NLOC + cv_kloc)
                TGI_ELE(1:NPHASE)=TGI_ELE(1:NPHASE) + SCVFEN(cv_Kloc,gi)*T_ALL(:,cv_nodk)
                TGI_ELE(1+NPHASE:2*NPHASE)=TGI_ELE(1+NPHASE:2*NPHASE) + SCVFEN(cv_Kloc,gi)*TOLD_ALL(:,cv_nodk)
            end do


            if(.true.) then

                ugi_normalised(:,1) = ugi(:,1) / max(1.e-7, sqrt(sum( ugi(:,1)**2) ))
                dx = 1.0 / max( 1.0e-7, maxval( abs(matmul( INV_JAC(:,:,GI), UGI_normalised(:,1)))  ) )

                TGI_der_ELE=0.0
                do cv_kloc=1,cv_nloc
                    cv_nodk = CV_NDGLN((ELE-1)*CV_NLOC + cv_kloc)
                    do iphase = 1, nphase
                        TGI_der_ELE(:,iphase)=TGI_der_ELE(:,iphase) + SCVFENx(:, cv_Kloc,gi)*T_ALL(iphase,cv_nodk)
                        TGI_der_ELE(:,iphase+NPHASE)=TGI_der_ELE(:,iphase+NPHASE) + SCVFENx(:,cv_Kloc,gi)*TOLD_ALL(iphase,cv_nodk)
                    end do
                end do

                do iphase = 1, nphase*2
                    TGI_ele_b( iphase ) = TGI_ele( iphase ) - 0.5 * dx * sum( ugi_normalised(:,1) * TGI_der_ELE(:,iphase) )
                    TGI_ele_c( iphase ) = TGI_ele( iphase ) + 0.5 * dx * sum( ugi_normalised(:,1) * TGI_der_ELE(:,iphase) )!opposite
                end do


                if ( .true. ) then
                    TUP(1:NPHASE)=(1.0-INCOME(1:NPHASE))*T_ALL(1:NPHASE,CV_NODI) + INCOME(1:NPHASE)*T_ALL(1:NPHASE,CV_NODJ)
                    TUP(1+NPHASE:2*NPHASE)=(1.0-INCOMEOLD(1:NPHASE))*TOLD_ALL(1:NPHASE,CV_NODI) + INCOMEOLD(1:NPHASE)*TOLD_ALL(1:NPHASE,CV_NODJ)


                    do iphase = 1, nphase*2

                        max_val=max( TGI_ele_b(iphase),TGI_ele_c(iphase), TGI_ele(iphase) )
                        min_val=min( TGI_ele_b(iphase),TGI_ele_c(iphase), TGI_ele(iphase) )


                        if ( TUP( iphase) < max_val .and. TUP( iphase) > min_val ) then
                            TGI_ele(iphase) = TUP( iphase)
                        else
                            if ( abs( TGI_ele_b( iphase ) - TUP( iphase ) )  < abs( TGI_ele( iphase ) - TUP( iphase ) ) ) then
                                if ( abs( TGI_ele_c( iphase ) - TUP( iphase ) )  < abs( TGI_ele_b( iphase ) - TUP( iphase ) ) ) then
                                    TGI_ele(iphase) = TGI_ele_C(iphase)
                                else
                                    TGI_ele(iphase) = TGI_ele_B(iphase)
                                end if
                            else
                                if ( abs( TGI_ele_c( iphase ) - TUP( iphase ) )  < abs( TGI_ele( iphase ) - TUP( iphase ) ) ) then
                                    TGI_ele(iphase) = TGI_ele_C(iphase)
                                end if
                            end if
                        end if
                    end do
                else
                    TGI_ele = TGI_ele_B
                endif

            end if


            IPT=1
            CALL PACK_LOC( LIMF(:), TGI_ele( 1:NPHASE ),    NPHASE, IPT, IGOT_T_PACK(:,1) ) ! T
            CALL PACK_LOC( LIMF(:), TGI_ele( 1+NPHASE:2*NPHASE ), NPHASE, IPT, IGOT_T_PACK(:,2) ) ! TOLD_ALL
            return
        end if


        TUP(1:NPHASE)=(1.0-INCOME(1:NPHASE))*T_ALL(1:NPHASE,CV_NODI) + INCOME(1:NPHASE)*T_ALL(1:NPHASE,CV_NODJ)
        TUP(1+NPHASE:2*NPHASE)=(1.0-INCOMEOLD(1:NPHASE))*TOLD_ALL(1:NPHASE,CV_NODI) + INCOMEOLD(1:NPHASE)*TOLD_ALL(1:NPHASE,CV_NODJ)

        !TUP(1:NPHASE)=(1.0-INCOME(1:NPHASE))*femT_ALL(1:NPHASE,CV_NODI) + INCOME(1:NPHASE)*FEMT_ALL(1:NPHASE,CV_NODJ)
        !TUP(1+NPHASE:2*NPHASE)=(1.0-INCOMEOLD(1:NPHASE))*femTOLD_ALL(1:NPHASE,CV_NODI) + INCOMEOLD(1:NPHASE)*FEMTOLD_ALL(1:NPHASE,CV_NODJ)

        W = 0.0 ! =0 means always apply ENO.
        IF(ENO_ONLY_WHERE_OSCILLATE) THEN
            ! Use ENO only where there is an oscillation as it can be a bit dissipative.
            !GOT_AN_OSC=.FALSE.
            IPT=1
            DO I_OLD_NEW=1,2
                DO IPHASE=1,NPHASE
                    IF(IGOT_T_PACK(IPHASE,1)) THEN
                        IPHASE2=IPHASE + (I_OLD_NEW-1)*NPHASE
                        W(IPHASE2)=(TUP(IPHASE2)-LIMF(IPT))/TOLFUN(TUP(IPHASE2)-TGI_ELE(IPHASE2))
                        W(IPHASE2)=MAX(0.0,MIN(1.0,W(IPHASE2) ))
                        ! W=0.0(full upwind);  W=1.0(high order no limiting)
                        !IF(W(IPHASE2)<0.9999) GOT_AN_OSC=.TRUE.
                        IPT=IPT+1
                    ENDIF
                END DO ! ENDOF DO IPHASE=1,NPHASE
            END DO ! ENDOF DO I_OLD_NEW=1,2
           ! Nothing to do if no oscillations detected...
        !        IF(.NOT.GOT_AN_OSC) RETURN
        ENDIF

        IF(DISTCONTINUOUS_METHOD.AND.(.NOT.BETWEEN_ELEMENTS).AND.FOR_DG_ONLY_BETWEEN_ELE) THEN ! Use fem INSIDE element for DG method...
            TGI = TGI_ELE
        ELSE

            xpt_GI=0.0
            do cv_kloc=1,cv_nloc
                X_knod=x_ndgln((ele-1)*cv_nloc+cv_kloc)
                xpt_GI(:)=xpt_GI(:)+SCVFEN(cv_kloc,gi)*x_all(:,X_kNOD)
            end do

            IS_CORNER_NOD_I=.TRUE.
            IS_CORNER_NOD_J=.TRUE.
            IF ( between_elements ) THEN
                XVEC(:)=0.0
            else
                IF(QUADRATIC_ELEMENT) THEN
                    XVEC(:)= - 0.75*(X_ALL(:,X_NODI)-X_ALL(:,X_NODJ)) ! Double the length scale because its a quadratic element
                    ! Is CV_JLOC a corner node...
                    IS_CORNER_NOD_I = (CV_ILOC==1).OR.(CV_ILOC==3).OR.(CV_ILOC==6).OR.(CV_ILOC==10)
                    IS_CORNER_NOD_J = (CV_JLOC==1).OR.(CV_JLOC==3).OR.(CV_JLOC==6).OR.(CV_JLOC==10)
                    !IF(IS_CORNER_NOD_J) XVEC(:)= - 0.75*(X_ALL(:,X_NODI)-X_ALL(:,X_NODJ)) ! half length scale because we are next to element boundary.
                    IF(IS_CORNER_NOD_I) XVEC(:)= + 0.75*(X_ALL(:,X_NODI)-X_ALL(:,X_NODJ)) ! half length scale because we are next to element boundary.
                ELSE ! linear element..
                    XVEC(:)= - 0.75*(X_ALL(:,X_NODI)-X_ALL(:,X_NODJ))
                 ENDIF
            endif

            DO IUP_DOWN=1,2 ! Consider both directions...

                RUP_WIN=REAL(2.-IUP_DOWN)*2.0 - 1.0
                XPT(:) = XPT_GI(:) + RUP_WIN*XVEC(:)
                IF(IUP_DOWN==2) THEN ! we might just consider the 2 neighbouring elements.
                    IF(  BETWEEN_ELEMENTS  &
                        .OR. (QUADRATIC_ELEMENT.AND.(IS_CORNER_NOD_I.or.IS_CORNER_NOD_J)) ) THEN
                        ENO_ELE_NEI(IUP_DOWN)  = ELEWIC
                        ENO_ELE_MATWEI(:,IUP_DOWN) = N(:,1)
                        EXIT ! Jump out of IUP_DOWN loop
                    ENDIF
                ENDIF

                ! Search for element this pt belongs to or nearest element excluding current element.
                MINCORK=-INFINY
                ELEWIC=0

                DO IFACE = 1, NFACE
                    ELE2 = FACE_ELE( IFACE, ELE )
                    SELE2 = MAX( 0, - ELE2 )
                    ELE2 = MAX( 0, + ELE2 )
                    IF(ELE2.NE.0) THEN

                        ! Neighbouring element...
                        IF(QUADRATIC_ELEMENT) THEN
                            LOCNODS(1)=X_NDGLN((ELE2-1)*CV_NLOC+1)
                            LOCNODS(2)=X_NDGLN((ELE2-1)*CV_NLOC+3)
                            LOCNODS(3)=X_NDGLN((ELE2-1)*CV_NLOC+6)
                            IF(NDIM==3) LOCNODS(4)=X_NDGLN((ELE2-1)*CV_NLOC+10)
                        ELSE
                            LOCNODS(1:CV_NLOC)=X_NDGLN((ELE2-1)*CV_NLOC+1:(ELE2-1)*CV_NLOC+CV_NLOC)
                        ENDIF

                        CALL TRI_tet_LOCCORDS(Xpt, LOCCORDS,  &
                             ! The 3 corners of the tri...
                             X_ALL(:,LOCNODS(:)),NDIM,NDIM+1)

                        MINCOR=MINVAL( LOCCORDS(:) )

                        IF(MINCOR.GT.MINCORK) THEN
                            MINCORK=MINCOR
                            ELEWIC=ELE2
                        ENDIF

                    ENDIF ! ENDOF IF(ELE2.NE.0) THEN

                END DO ! ENDOF DO IFACE=1,NFACE

                ! FOR ELEMENT ELEWIC FIND THE VOLN COORD OF POINT
                IF(ELEWIC==0) ELEWIC=ELE
                ELE2=ELEWIC

                ! NEIGHBOURING element ELE2...

                IF(QUADRATIC_ELEMENT) THEN
                    LOCNODS(1)=X_NDGLN((ELE2-1)*CV_NLOC+1)
                    LOCNODS(2)=X_NDGLN((ELE2-1)*CV_NLOC+3)
                    LOCNODS(3)=X_NDGLN((ELE2-1)*CV_NLOC+6)
                    IF(NDIM==3) LOCNODS(4)=X_NDGLN((ELE2-1)*CV_NLOC+10)
                ELSE
                    LOCNODS(1:CV_NLOC)=X_NDGLN((ELE2-1)*CV_NLOC+1:(ELE2-1)*CV_NLOC+CV_NLOC)
                ENDIF

                CALL TRI_tet_LOCCORDS(Xpt_GI, LOCCORDS,&
                     ! The 3 corners of the tri...
                     X_ALL(:,LOCNODS(:)),NDIM,NDIM+1)

                L1(1) = LOCCORDS(1)
                L2(1) = LOCCORDS(2)
                L3(1) = LOCCORDS(3)
                L4(1) = LOCCORDS(NDIM+1)

                ! From  the local coordinates find the shape function value...
                call shatri( l1, l2, l3, l4, weight, (NDIM==3), &
                    cv_nloc, ngi_one, &
                    n, nlx, nly, nlz )

                ! could store these as STORE_ENO_ELE_NEI(IUP_DOWN,GI), ENO_ELE_MATWEI(CV_NLOC,IUP_DOWN,GI)
                ENO_ELE_NEI(IUP_DOWN)  = ELEWIC
                ENO_ELE_MATWEI(:,IUP_DOWN) = N(:,1)

            END DO ! DO IUP_DOWN=1,2

            ! now calculate ENO values at the quadrature pt **********************************...
            ! now calculate ENO values at the quadrature pt **********************************...
            !
            ENO_ELE_NEI_IN =ENO_ELE_NEI(1)
            ENO_ELE_NEI_OUT=ENO_ELE_NEI(2)

            ENO_ELE_MATWEI_IN(:) =ENO_ELE_MATWEI(:,1)
            ENO_ELE_MATWEI_OUT(:)=ENO_ELE_MATWEI(:,2)

            TGI_IN=0.0 ; TGI_OUT=0.0

            if (.true.) then

                do cv_kloc=1,cv_nloc
                    cv_nodk_IN  = CV_NDGLN((ENO_ELE_NEI_IN -1)*CV_NLOC + cv_kloc)
                    cv_nodk_OUT = CV_NDGLN((ENO_ELE_NEI_OUT-1)*CV_NLOC + cv_kloc)

                    TGI_IN(1:NPHASE)=TGI_IN(1:NPHASE) + ENO_ELE_MATWEI_IN(CV_KLOC)*FEMT_ALL(:,cv_nodk_IN)
                    TGI_OUT(1:NPHASE)=TGI_OUT(1:NPHASE) + ENO_ELE_MATWEI_OUT(CV_KLOC)*FEMT_ALL(:,cv_nodk_OUT)

                    TGI_IN(1+NPHASE:2*NPHASE)=TGI_IN(1+NPHASE:2*NPHASE) + ENO_ELE_MATWEI_IN(CV_KLOC)*FEMTOLD_ALL(:,cv_nodk_IN)
                    TGI_OUT(1+NPHASE:2*NPHASE)=TGI_OUT(1+NPHASE:2*NPHASE) + ENO_ELE_MATWEI_OUT(CV_KLOC)*FEMTOLD_ALL(:,cv_nodk_OUT)
                end do

            else
                do cv_kloc=1,cv_nloc
                    cv_nodk_IN  = CV_NDGLN((ENO_ELE_NEI_IN -1)*CV_NLOC + cv_kloc)
                    cv_nodk_OUT = CV_NDGLN((ENO_ELE_NEI_OUT-1)*CV_NLOC + cv_kloc)

                    TGI_IN(1:NPHASE)=TGI_IN(1:NPHASE) + ENO_ELE_MATWEI_IN(CV_KLOC)*T_ALL(:,cv_nodk_IN)
                    TGI_OUT(1:NPHASE)=TGI_OUT(1:NPHASE) + ENO_ELE_MATWEI_OUT(CV_KLOC)*T_ALL(:,cv_nodk_OUT)

                    TGI_IN(1+NPHASE:2*NPHASE)=TGI_IN(1+NPHASE:2*NPHASE) + ENO_ELE_MATWEI_IN(CV_KLOC)*TOLD_ALL(:,cv_nodk_IN)
                    TGI_OUT(1+NPHASE:2*NPHASE)=TGI_OUT(1+NPHASE:2*NPHASE) + ENO_ELE_MATWEI_OUT(CV_KLOC)*TOLD_ALL(:,cv_nodk_OUT)
                end do
            end if


            IF(ENO_ALL_THREE) THEN ! Use all 3 to find TGI ENO value...

               IF( IS_CORNER_NOD_I .or. IS_CORNER_NOD_J ) THEN

                  INCOME_BOTH(1:NPHASE)         =INCOME(1:NPHASE)
                  INCOME_BOTH(1+NPHASE:2*NPHASE)=INCOMEOLD(1:NPHASE)

                  IF(IS_CORNER_NOD_I) THEN
                        TGI_NEI =INCOME_BOTH       * TGI_ELE + (1.0-INCOME_BOTH) * TGI_IN
                        TGI_NEI2=(1.0-INCOME_BOTH) * TGI_ELE + INCOME_BOTH       * TGI_IN
                    ELSE
                        TGI_NEI =INCOME_BOTH       * TGI_IN + (1.0-INCOME_BOTH) * TGI_ELE
                        TGI_NEI2=(1.0-INCOME_BOTH) * TGI_IN + INCOME_BOTH       * TGI_ELE
                    ENDIF
                    !
                    TGI_IN=TGI_NEI
                    ! reduce influence of downwind value (its needed as it oscillates without this)
                    TGI_OUT=RELAX_DOWN_WIND_2_CURRENT_ELE * TGI_ELE + (1.0-RELAX_DOWN_WIND_2_CURRENT_ELE) * TGI_NEI2
                ELSE
                    ! Upwind value...
                    TGI_NEI(1:NPHASE)         =INCOME(1:NPHASE   )*TGI_IN(1:NPHASE         ) + (1.0-INCOME(1:NPHASE))   *TGI_OUT(1:NPHASE)
                    TGI_NEI(1+NPHASE:2*NPHASE)=INCOMEOLD(1:NPHASE)*TGI_IN(1+NPHASE:2*NPHASE) + (1.0-INCOMEOLD(1:NPHASE))*TGI_OUT(1+NPHASE:2*NPHASE)
                    ! Downwind value...
                    TGI_NEI2(1:NPHASE)         =(1.0-INCOME(1:NPHASE)   )*TGI_IN(1:NPHASE)          + INCOME   (1:NPHASE)*TGI_OUT(1:NPHASE)
                    TGI_NEI2(1+NPHASE:2*NPHASE)=(1.0-INCOMEOLD(1:NPHASE))*TGI_IN(1+NPHASE:2*NPHASE) + INCOMEOLD(1:NPHASE)*TGI_OUT(1+NPHASE:2*NPHASE)
                    !
                    TGI_IN=TGI_NEI
                    ! reduce influence of downwind value (its needed as it oscillates without this)
                    TGI_OUT=RELAX_DOWN_WIND_2_CURRENT_ELE * TGI_ELE + (1.0-RELAX_DOWN_WIND_2_CURRENT_ELE) * TGI_NEI2
                ENDIF

                DO IPHASE2=1,NPHASE2
                    MIN_TGI=MIN(TGI_ELE(IPHASE2),TGI_IN(IPHASE2),TGI_OUT(IPHASE2))
                    MAX_TGI=MAX(TGI_ELE(IPHASE2),TGI_IN(IPHASE2),TGI_OUT(IPHASE2))
                    IF((TUP(IPHASE2)-MIN_TGI)*(TUP(IPHASE2)-MAX_TGI)>0) THEN
                        TGI(IPHASE2)=TUP(IPHASE2)
                    ELSE ! Choose the closest
                        IF(ABS(TUP(IPHASE2)-MIN_TGI)<ABS(TUP(IPHASE2)-MAX_TGI)) THEN
                            TGI(IPHASE2)=MIN_TGI
                        ELSE
                            TGI(IPHASE2)=MAX_TGI
                        ENDIF
                    ENDIF
                END DO
            ELSE ! Give an upwind bias (the Default)...
                TGI_NEI(1:NPHASE)=INCOME(1:NPHASE)*TGI_IN(1:NPHASE) + (1.0-INCOME(1:NPHASE))*TGI_OUT(1:NPHASE)
                TGI_NEI(1+NPHASE:2*NPHASE)=INCOMEOLD(1:NPHASE)*TGI_IN(1+NPHASE:2*NPHASE) + (1.0-INCOMEOLD(1:NPHASE))*TGI_OUT(1+NPHASE:2*NPHASE)
                ! Choose TGI that is closest to TUP...
                DO IPHASE2=1,NPHASE2

if (.true.) then

                    IF(  (TUP(IPHASE2)-TGI_ELE(IPHASE2))*(TUP(IPHASE2)-TGI_NEI(IPHASE2)) >0 ) THEN
                        TGI(IPHASE2)=TUP(IPHASE2)
                    ELSE ! Choose the closest
                        IF( ABS(TUP(IPHASE2)-TGI_ELE(IPHASE2)) < ABS(TUP(IPHASE2)-TGI_NEI(IPHASE2))  ) THEN
                            TGI(IPHASE2)=TGI_ELE(IPHASE2)
                        ELSE
                            TGI(IPHASE2)=TGI_NEI(IPHASE2)
                        ENDIF
                    ENDIF
else
                   TGI(IPHASE2) = min( TUP(IPHASE2), TGI_NEI(IPHASE2) )
end if

                END DO
            ENDIF

            ! Use ENO only where there is an oscillation as it can be a bit dissipative.
            TGI(:) = (1.0-W(:))*TGI(:) + W(:)*TGI_ELE(:)

        ENDIF ! ENDOF IF(DISTCONTINUOUS_METHOD.AND.(.NOT.BETWEEN_ELEMENTS).AND.FOR_DG_ONLY_BETWEEN_ELE) THEN else

        IPT=1
        CALL PACK_LOC( LIMF(:), TGI( 1:NPHASE ),    NPHASE, IPT, IGOT_T_PACK(:,1) ) ! T
        CALL PACK_LOC( LIMF(:), TGI( 1+NPHASE:2*NPHASE ), NPHASE, IPT, IGOT_T_PACK(:,2) ) ! TOLD_ALL

    END SUBROUTINE APPLY_ENO_2_T

        
    SUBROUTINE TRI_tet_LOCCORDS(Xpt, LOCCORDS,  &
         !     The 3 corners of the tri...
         X_CORNERS_ALL, NDIM,CV_NLOC)
      ! obtain the local coordinates LOCCORDS from a pt in or outside the tet/triangle Xpt
      ! with corner nodes X_CORNERS_ALL
      IMPLICIT NONE
      INTEGER, intent(in) :: NDIM,CV_NLOC
      REAL, dimension(NDIM), intent(in) :: Xpt
      REAL, dimension(NDIM+1), intent(inout) :: LOCCORDS
      REAL, dimension(NDIM,CV_NLOC), intent(in) :: X_CORNERS_ALL

      IF (NDIM==3) THEN
         CALL TRILOCCORDS(Xpt(1),Xpt(2),Xpt(3), &
              LOCCORDS(1),LOCCORDS(2),LOCCORDS(3),LOCCORDS(4),&
              !     The 4 corners of the tet...
              X_CORNERS_ALL(1,1),X_CORNERS_ALL(2,1),X_CORNERS_ALL(3,1),&
              X_CORNERS_ALL(1,2),X_CORNERS_ALL(2,2),X_CORNERS_ALL(3,2),&
              X_CORNERS_ALL(1,3),X_CORNERS_ALL(2,3),X_CORNERS_ALL(3,3),&
              X_CORNERS_ALL(1,4),X_CORNERS_ALL(2,4),X_CORNERS_ALL(3,4) )
      ELSE
         CALL TRILOCCORDS2D(Xpt(1),Xpt(2), &
              LOCCORDS(1),LOCCORDS(2),LOCCORDS(3),&
              !     The 3 corners of the tri...
              X_CORNERS_ALL(1,1),X_CORNERS_ALL(2,1),&
              X_CORNERS_ALL(1,2),X_CORNERS_ALL(2,2),&
              X_CORNERS_ALL(1,3),X_CORNERS_ALL(2,3) )
      END IF
      ! From  the local coordinates find the shape function value...
      RETURN
    END SUBROUTINE TRI_tet_LOCCORDS

          subroutine dump_multiphase(prefix,icp)
            character(len=*), intent(in) :: prefix
            integer, optional :: icp
            integer, save :: its=0
            integer :: ip, lcomp
            type( scalar_field ), dimension(2*Mdims%nphase) :: tcr
            type( vector_field ), dimension(2*Mdims%nphase) :: vel
            type( vector_field ), pointer :: position
            position=>extract_vector_field(state(1),"Coordinate")
            if (present(icp)) then
                lcomp=icp
            else
                lcomp=0
            end if
            do ip=1,Mdims%nphase
                tcr(ip)=wrap_scalar_field(tracer%mesh,T_ALL(ip,:),&
                    'TracerPhase'//int2str(ip))
                tcr(ip+Mdims%nphase)=wrap_scalar_field(tracer%mesh,TOLD_ALL(ip,:),&
                    'OldTracerPhase'//int2str(ip))
                vel(ip)=wrap_vector_field(velocity%mesh,NU_ALL(:,ip,:),&
                    'NLVelocityPhase'//int2str(ip))
                vel(ip+Mdims%nphase)=wrap_vector_field(velocity%mesh,NUOLD_ALL(:,ip,:),&
                    'OldNLVelocityPhase'//int2str(ip))
            end do
            call vtk_write_fields(prefix//trim(tracer%name)&
                //'Component'//int2str(lcomp),&
                its,position,tracer%mesh,sfields=tcr,vfields=vel)
            do ip=1,2*Mdims%nphase
                call deallocate(tcr(ip))
                call deallocate(vel(ip))
            end do
            its=its+1
        end subroutine dump_multiphase

        SUBROUTINE GET_INT_T_DEN_new( LIMF )
            !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
            IMPLICIT NONE
            ! Calculate T and DEN on the CV face at quadrature point GI.
            REAL, DIMENSION ( NFIELD), intent( inout ) :: LIMF
            !   REAL, DIMENSION( : ), intent( in ) :: TUPWIND_MAT, DENUPWIND_MAT
            !   REAL, DIMENSION( : ), intent( in ) :: T2UPWIND_MAT
            ! Local variables
            ! If UPWIND then use upwind flux between elements else use central.
            ! If HI_ORDER_HALF then use high order interpolation when around
            ! a volume frac of 0.5 and gradually apply limiting near 0 and 1.
            LOGICAL :: DOWNWIND_EXTRAP ! Extrapolate a downwind value for interface tracking.
            ! Scaling to reduce the downwind bias(=1downwind, =0central)
            LOGICAL, PARAMETER :: SCALE_DOWN_WIND = .true.
            ! Non-linear Petrov-Galerkin option for interface tracking...
            ! =4 is anisotropic downwind diffusion based on a projected 1D system (1st recommend)
            ! =0 is anisotropic downwind diffusion based on a velocity projection like SUPG
            ! (2nd recommend, most compressive)
            ! =2 is isotropic downwind diffusion  (3rd recommend,least compressive)
            ! =5 is isotropic downwind diffusion with magnitude of =0 option.
            ! In tests they all produce similar results.
            !      INTEGER, PARAMETER :: NON_LIN_PETROV_INTERFACE = 5
            INTEGER, PARAMETER :: NON_LIN_PETROV_INTERFACE = 3
            real, parameter :: tolerance = 1.e-10
            LOGICAL :: NOLIMI
            INTEGER :: CV_KLOC, &
                CV_SKLOC, IDIM
            INTEGER :: IFIELD
            REAL, dimension(Mdims%ndim,NFIELD) :: FXGI_ALL, UDGI_ALL, A_STAR_X_ALL,VEC_VEL2
            REAL, dimension(NFIELD) :: courant_or_minus_one_new, XI_LIMIT,&
                P_STAR, U_DOT_GRADF_GI, A_STAR_F, RESIDGI, ELE_LENGTH_SCALE,FEMFGI, RGRAY, DIFF_COEF, COEF,&
                RSCALE, COEF2, FEMFGI_CENT, FEMFGI_UP
            ! No limiting if CV_DISOPT is 6 or 7  (why not just define limt=femt and skip to assembly?)
            NOLIMI = ( INT( CV_DISOPT / 2 ) == 3 )
            ! Make a guess at the CV face value of advected field variable and density
            ! (Depends on discetisation option, CV_DISOPT)
            SELECT CASE( CV_DISOPT / 2 )
                !    CASE( 0 ) ! First-order upwinding is achived through the limiting
                !       FEMFGI(:)    = FVF(:)
                CASE( 1 ) ! Central differencing [Trapezoidal rule (2 OR 3)]
                    FEMFGI(:)    = 0.5 * ( LOC_F( :, CV_ILOC ) + LOC_F( :, CV_JLOC ) )
                CASE DEFAULT ! Finite element approximation (4 OR 5)(6 or 7)(8 or 9)
                    FEMFGI(:)    = 0.0
                    Conditional_CV_DISOPT_ELE2: IF ( on_domain_boundary ) THEN
                        ! Is on boundary of the domain
                        DO IFIELD=1,NFIELD
                            IF ( SELE_LOC_WIC_F_BC(IFIELD) /= WIC_T_BC_DIRICHLET ) THEN ! Don't apply a Dirichlet bc
                                LIMF(IFIELD)    = LOC_F( IFIELD, CV_ILOC )
                            ELSE
                                LIMF(IFIELD)    = (1.-F_INCOME(IFIELD))*LOC_F( IFIELD, CV_ILOC )   + F_INCOME(IFIELD)*  SLOC_SUF_F_BC( IFIELD,  CV_SILOC)
                            END IF
                        END DO ! END OF DO IFIELD=1,NFIELD
                    ELSE Conditional_CV_DISOPT_ELE2
                        ! Extrapolate a downwind value for interface tracking.
                        DOWNWIND_EXTRAP = ( cv_disopt>=8 )
                        DO IFIELD=1,NFIELD
                            IF( DOWNWIND_EXTRAP_INDIVIDUAL(IFIELD)  ) THEN
                                courant_or_minus_one_new(IFIELD) = abs ( dt * F_ndotq(IFIELD) / hdc )
                                XI_LIMIT(IFIELD)=MAX(1./max(tolerance,(3.*courant_or_minus_one_new(IFIELD))),2.0)
                            else
                                courant_or_minus_one_new(IFIELD) = -1.0
                                XI_LIMIT(IFIELD) = 2.0
                            end if
                        END DO
                        IF( .not. between_elements ) THEN
                            RSCALE(:) = 1.0 ! Scaling to reduce the downwind bias(=1downwind, =0central)
                            IF ( SCALE_DOWN_WIND ) THEN
                                IF ( DOWNWIND_EXTRAP  ) THEN
                                    DO IFIELD=1,MIN(2*Mdims%nphase,NFIELD)
                                        DO IDIM=1,Mdims%ndim
                                            FXGI_ALL(IDIM,IFIELD) = dot_product(SdevFuns%NX_ALL(IDIM, : , GI ) , LOC_FEMF(IFIELD,:))
                                        END DO
                                        !FEMFGI(IFIELD) = dot_product( CV_funs%scvfen( : , GI ), LOC_FEMF(IFIELD,:)  )
                                        DO IDIM=1,Mdims%ndim
                                            UDGI_ALL(IDIM,IFIELD) = dot_product(CV_funs%sufen( : , GI ) , LOC_UF(IDIM,IFIELD, :) )
                                        END DO
                                    END DO
                                    IF ( NON_LIN_PETROV_INTERFACE == 0 ) THEN ! NOT Petrov-Galerkin for interface capturing...
                                        ! no cosine rule :
                                        DO IFIELD=1,MIN(2*Mdims%nphase,NFIELD) ! It should be Mdims%nphase because interface tracking only applied to the 1st set of fields.
                                            RSCALE(IFIELD) = 1.0 / PTOLFUN( SQRT( SUM( UDGI_ALL(:,IFIELD)**2)   ) )
                                            DO IDIM = 1, Mdims%ndim
                                                VEC_VEL2(IDIM,IFIELD) = SUM( SdevFuns%INV_JAC(IDIM, 1:Mdims%ndim, GI) * UDGI_ALL(1:Mdims%ndim,IFIELD) )
                                            END DO
                                            ! normalize the velocity in here:
                                            ELE_LENGTH_SCALE(IFIELD) = 0.5 * SQRT( SUM(UDGI_ALL(:,IFIELD)**2) ) / PTOLFUN( SUM( VEC_VEL2(:,IFIELD)**2 ) )
                                            ! For discontinuous elements half the length scale...
                                            IF(DISTCONTINUOUS_METHOD) ELE_LENGTH_SCALE(IFIELD)=0.5*ELE_LENGTH_SCALE(IFIELD)
                                            ! For quadratic elements...
                                            IF( QUAD_ELEMENTS ) ELE_LENGTH_SCALE(IFIELD)=0.5*ELE_LENGTH_SCALE(IFIELD)
                                        END DO
                                    ELSE ! Interface capturing...
                                        DO IFIELD=1,MIN(2*Mdims%nphase,NFIELD) ! It should be Mdims%nphase because interface tracking only applied to the 1st set of fields.
                                            U_DOT_GRADF_GI(IFIELD) = SUM( UDGI_ALL(:,IFIELD)*FXGI_ALL(:,IFIELD)  )
                                            IF ( NON_LIN_PETROV_INTERFACE == 5 ) THEN
                                                COEF(IFIELD) = 1.0 / PTOLFUN( SQRT( SUM(FXGI_ALL(IFIELD,:)**2)   ) )
                                            ELSE
                                                COEF(IFIELD) = U_DOT_GRADF_GI(IFIELD) / PTOLFUN( SUM( FXGI_ALL(:,IFIELD)**2 )  )
                                            END IF
                                            A_STAR_F(IFIELD) = 0.0
                                            A_STAR_X_ALL(:,IFIELD) = COEF(IFIELD) * FXGI_ALL(:,IFIELD)
                                            RESIDGI(IFIELD) = SQRT ( SUM( UDGI_ALL(:,IFIELD)**2 )  ) / HDC
                                            VEC_VEL2(1:Mdims%ndim,IFIELD) = matmul( SdevFuns%INV_JAC(:,:,GI), A_STAR_X_ALL(1:Mdims%ndim,IFIELD) )
                                            ! Needs 0.25 for quadratic elements...Chris
                                            P_STAR(IFIELD) = 0.5 * HDC / PTOLFUN( SQRT( SUM( A_STAR_X_ALL(:,IFIELD)**2 )))
                                            !                                      IF( QUAD_ELEMENTS ) P_STAR(IFIELD) = 0.5 * P_STAR(IFIELD)
                                            select case (NON_LIN_PETROV_INTERFACE)
                                                case ( 1 )     ! standard approach
                                                    DIFF_COEF(IFIELD) = COEF(IFIELD) * P_STAR(IFIELD) * RESIDGI(IFIELD)
                                                case ( 2 )     ! standard approach making it +ve
                                                    DIFF_COEF(IFIELD) = MAX( 0.0, COEF(IFIELD) * P_STAR(IFIELD) * RESIDGI(IFIELD) )
                                                case ( 3 )     ! residual squared approach
                                                    DIFF_COEF(IFIELD) = P_STAR(IFIELD) * RESIDGI(IFIELD)**2 / PTOLFUN( SUM(FXGI_ALL(:,IFIELD)**2)  )
                                                case ( 4 )     ! anisotropic diffusion in the A* direction.
                                                    COEF2(IFIELD) =  SUM( CVNORMX_ALL(:,GI)*A_STAR_X_ALL(:,IFIELD) )
                                                case ( 6 )     ! accurate...
                                                    P_STAR(IFIELD)=0.5/PTOLFUN( maxval(abs(VEC_VEL2(1:Mdims%ndim,IFIELD)))  )
                                                    IF( QUAD_ELEMENTS ) P_STAR(IFIELD) = 0.5 * P_STAR(IFIELD)
                                                    RESIDGI(IFIELD)=SUM( UDGI_ALL(:,IFIELD)*SdevFuns%NX_ALL( :, CV_KLOC, GI ) )
                                                    DIFF_COEF(IFIELD) = P_STAR(IFIELD) * RESIDGI(IFIELD)**2 / PTOLFUN( SUM(FXGI_ALL(:,IFIELD)**2)  )
                                                case ( 7 )     ! accurate (could be positive or negative)...
                                                    P_STAR(IFIELD)=0.5/PTOLFUN( maxval(abs(VEC_VEL2(1:Mdims%ndim,IFIELD)))  )
                                                    IF( QUAD_ELEMENTS ) P_STAR(IFIELD) = 0.5 * P_STAR(IFIELD)
                                                    RESIDGI(IFIELD)=SUM( UDGI_ALL(:,IFIELD)*SdevFuns%NX_ALL( :, CV_KLOC, GI ) )
                                                    DIFF_COEF(IFIELD) = P_STAR(IFIELD) * RESIDGI(IFIELD)*SUM(CVNORMX_ALL(:,GI)*UDGI_ALL(:,IFIELD)) * ((LOC_F( IFIELD, CV_JLOC )-LOC_F( IFIELD, CV_ILOC ))/hdc) &
                                                        / PTOLFUN( SUM(FXGI_ALL(:,IFIELD)**2)  )
                                                case ( 8 )     ! accurate (force to be positive)...
                                                    P_STAR(IFIELD)=0.5/PTOLFUN( maxval(abs(VEC_VEL2(1:Mdims%ndim,IFIELD)))  )
                                                    IF( QUAD_ELEMENTS ) P_STAR(IFIELD) = 0.5 * P_STAR(IFIELD)
                                                    RESIDGI(IFIELD)=SUM( UDGI_ALL(:,IFIELD)*SdevFuns%NX_ALL( :, CV_KLOC, GI ) )
                                                    DIFF_COEF(IFIELD) = P_STAR(IFIELD) * RESIDGI(IFIELD)*SUM(CVNORMX_ALL(:,GI)*UDGI_ALL(:,IFIELD)) * ((LOC_F( IFIELD, CV_JLOC )-LOC_F( IFIELD, CV_ILOC ))/hdc) &
                                                        / PTOLFUN( SUM(FXGI_ALL(:,IFIELD)**2)  )
                                                    DIFF_COEF(IFIELD) = abs(   DIFF_COEF(IFIELD)  )
                                                case ( 9 )     ! accurate (simplified residual squared)...
                                                    P_STAR(IFIELD)=0.5/PTOLFUN( maxval(abs(VEC_VEL2(1:Mdims%ndim,IFIELD)))  )
                                                    IF( QUAD_ELEMENTS ) P_STAR(IFIELD) = 0.5 * P_STAR(IFIELD)
                                                    RESIDGI(IFIELD)=SUM( UDGI_ALL(:,IFIELD)*SdevFuns%NX_ALL( :, CV_KLOC, GI ) )
                                                    DIFF_COEF(IFIELD) = P_STAR(IFIELD) * (SUM(CVNORMX_ALL(:,GI)*UDGI_ALL(:,IFIELD)) * ((LOC_F( IFIELD, CV_JLOC )-LOC_F( IFIELD, CV_ILOC ))/hdc))**2 &
                                                        / PTOLFUN( SUM(FXGI_ALL(:,IFIELD)**2)  )
                                                case default   ! isotropic diffusion with u magnitide
                                                    DIFF_COEF(IFIELD) = SQRT( SUM( UDGI_ALL(:,IFIELD)**2 )) * P_STAR(IFIELD)
                                            END select
                                            ! Make the diffusion coefficient negative (compressive)
                                            DIFF_COEF(IFIELD) = -DIFF_COEF(IFIELD)
                                            RSCALE(IFIELD) = 1. / TOLFUN( SUM(CVNORMX_ALL(:,GI)*UDGI_ALL(:,IFIELD))   )
                                        END DO ! END OF DO IFIELD=1,NFIELD
                                    END IF ! Petrov-Galerkin end of IF(NON_LIN_PETROV_INTERFACE==0) THEN
                                END IF ! DOWNWIND_EXTRAP
                            END IF ! SCALE_DOWN_WIND
                            FEMFGI=0.0
                            DO IFIELD=1,NFIELD ! Only perform this loop for the 1st field which is the interface tracking field...
                                IF( DOWNWIND_EXTRAP_INDIVIDUAL(IFIELD)  ) THEN ! Extrapolate to the downwind value...
                                    DO CV_KLOC = 1, Mdims%cv_nloc
                                        IF ( NON_LIN_PETROV_INTERFACE.NE.0 ) THEN
                                            IF ( NON_LIN_PETROV_INTERFACE == 4 ) THEN ! anisotropic diffusion...
                                                RGRAY(IFIELD) = RSCALE(IFIELD) * COEF2(IFIELD) * P_STAR(IFIELD) * SUM( UDGI_ALL(:,IFIELD)*SdevFuns%NX_ALL( :, CV_KLOC, GI ) )
                                            ELSE
                                                RGRAY(IFIELD) = - DIFF_COEF(IFIELD) * RSCALE(IFIELD) * SUM( CVNORMX_ALL(:,GI)*SdevFuns%NX_ALL( :, CV_KLOC, GI )  )
                                            END IF
                                        ELSE
                                            RGRAY(IFIELD) = RSCALE(IFIELD) * ELE_LENGTH_SCALE(IFIELD) * SUM( UDGI_ALL(:,IFIELD)*SdevFuns%NX_ALL( :, CV_KLOC, GI ) )
                                        END IF
                                        FEMFGI(IFIELD)    = FEMFGI(IFIELD)     +  (CV_funs%scvfen( CV_KLOC, GI ) + RGRAY(IFIELD))   * LOC_FEMF( IFIELD, CV_KLOC)
                                    END DO ! ENDOF DO CV_KLOC = 1, Mdims%cv_nloc
                                ELSE
                                    DO CV_KLOC = 1, Mdims%cv_nloc
                                        FEMFGI(IFIELD)    = FEMFGI(IFIELD)     +  CV_funs%scvfen( CV_KLOC, GI ) * LOC_FEMF( IFIELD, CV_KLOC)
                                    END DO ! ENDOF DO CV_KLOC = 1, Mdims%cv_nloc
                                !                              FEMFGI(:)    = 0.5 * ( LOC_F( :, CV_ILOC ) + LOC_F( :, CV_JLOC ) )
                                END IF
                            END DO ! ENDOF DO IFIELD=1,Mdims%nphase
                        ELSE  ! END OF IF( .not. between_elements ) THEN  ---DG saturation across elements
                            FEMFGI_CENT(:) = 0.0
                            FEMFGI_UP(:)   = 0.0
                            DO CV_SKLOC = 1, Mdims%cv_snloc
                                ! Central for DG...
                                FEMFGI_CENT(:) = FEMFGI_CENT(:) +  SHAPE_CV_SNL( CV_SKLOC ) * 0.5 * ( SLOC_FEMF( :, CV_SKLOC ) &
                                    + SLOC2_FEMF( :, CV_SKLOC )    )
                                ! Standard DG upwinding...
                                FEMFGI_UP(:) = FEMFGI_UP(:) +  SHAPE_CV_SNL( CV_SKLOC ) * ( SLOC2_FEMF( :, CV_SKLOC)  &
                                    * F_INCOME(:) + SLOC_FEMF( :, CV_SKLOC) * ( 1. - F_INCOME(:) ) )
                            END DO
                            DO IFIELD=1,NFIELD
                                IF( DOWNWIND_EXTRAP_INDIVIDUAL(IFIELD)  ) THEN ! Extrapolate to the downwind value...
                                    FEMFGI(IFIELD) = FEMFGI_CENT(IFIELD)
                                ELSE
                                    FEMFGI(IFIELD) = FEMFGI_UP(IFIELD)
                                ENDIF
                            END DO
                        ENDIF ! ENDOF IF( ( ELE2 == 0 ) .OR. ( ELE2 == ELE ) ) THEN ELSE
                        CALL ONVDLIM_ANO_MANY( NFIELD, &
                            LIMF(:), FEMFGI(:), F_INCOME(:), &
                            F_CV_NODI(:), F_CV_NODJ(:),XI_LIMIT(:),  &
                            FUPWIND_IN(:), FUPWIND_OUT(:) )
                    ENDIF Conditional_CV_DISOPT_ELE2
            END SELECT
            RETURN
        END SUBROUTINE GET_INT_T_DEN_NEW
        SUBROUTINE GET_INT_VEL_ORIG_NEW( NDOTQNEW, NDOTQ, INCOME, &
            LOC_T_I, LOC_T_J, LOC_DEN_I, LOC_DEN_J, &
            LOC_NU, LOC2_NU, NUGI_ALL, &
            UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL,not_OLD_VEL)
            ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI.
            IMPLICIT NONE
            REAL, DIMENSION( : ), intent( inout ) :: NDOTQ, INCOME, NDOTQNEW
            REAL, DIMENSION( :, : ), intent( inout ) :: NUGI_ALL
            REAL, DIMENSION( : ), intent( in ) :: LOC_T_I, LOC_T_J, LOC_DEN_I, LOC_DEN_J
            REAL, DIMENSION( :, :, : ), intent( in ) :: LOC_NU, LOC2_NU
            REAL, DIMENSION( :, :, : ), intent( inout ) :: UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL
            logical, intent(in) :: not_OLD_VEL
            ! Local variables
            INTEGER :: U_KLOC,U_KLOC2,U_SKLOC
            ! Local variable for indirect addressing
            INTEGER :: IPHASE, IDIM
            ! coefficients for this element ELE
            UGI_COEF_ELE_ALL = 0.0
            ! coefficients for this element ELE2
            UGI_COEF_ELE2_ALL = 0.0
            !    Conditional_SELE: IF( SELE /= 0 ) THEN ! On the boundary of the domain.
            Conditional_SELE: IF( on_domain_boundary ) THEN ! On the boundary of the domain.
                DO IPHASE = 1, Mdims%nphase
                    IF( WIC_U_BC_ALL( 1, IPHASE, SELE) /= WIC_U_BC_DIRICHLET ) THEN ! velocity free boundary
                        UDGI_ALL(:, IPHASE) = 0.0
                        DO U_KLOC = 1, Mdims%u_nloc
                            UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) + CV_funs%sufen( U_KLOC, GI ) * LOC_NU( :, IPHASE, U_KLOC )
                            UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC) = 1.0
                        END DO
                    ELSE ! Specified vel bc.
                        UDGI_ALL(:, IPHASE) = 0.0
                        DO U_SKLOC = 1, Mdims%u_snloc
                            U_KLOC = U_SLOC2LOC( U_SKLOC )
                            UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) + CV_funs%sufen( U_KLOC, GI ) * SUF_U_BC_ALL(:, IPHASE, Mdims%u_snloc* (SELE-1) + U_SKLOC )
                            UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC) = 0.0
                        END DO
                    END IF
                END DO
            ELSE ! Conditional_SELE. Not on the boundary of the domain.
                UDGI_ALL = 0.0
                DO U_KLOC = 1, Mdims%u_nloc
                    UDGI_ALL = UDGI_ALL + CV_funs%sufen( U_KLOC, GI ) * LOC_NU( :, :, U_KLOC )
                    UGI_COEF_ELE_ALL(:, :, U_KLOC) = 1.0
                END DO
                !       Conditional_ELE2: IF( ELE2 /= 0 ) THEN
                Conditional_ELE2: IF( between_elements ) THEN
                    UDGI2_ALL = 0.0
                    DO U_SKLOC = 1, Mdims%u_snloc
                        U_KLOC = U_SLOC2LOC(U_SKLOC)
                        U_KLOC2 = U_OTHER_LOC( U_KLOC )
                        !          DO U_KLOC = 1, Mdims%u_nloc
                        !             U_KLOC2 = U_OTHER_LOC( U_KLOC )
                        !             IF ( U_KLOC2 /= 0 ) THEN ! MAKE SURE WE DONT NEED THIS...
                        UDGI2_ALL = UDGI2_ALL + CV_funs%sufen( U_KLOC, GI ) * LOC2_NU(:, :, U_KLOC)
                        UGI_COEF_ELE2_ALL( :, :, U_KLOC2) = 1.0
                    !             ENDIF
                    END DO
                    IF( ABS( CV_DG_VEL_INT_OPT ) == 1 ) THEN
                        DT_I=1.0
                        DT_J=1.0
                    ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 2) THEN
                        DT_I=MAX(1.E-2,LOC_T_I)
                        DT_J=MAX(1.E-2,LOC_T_J)
                    ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 3) THEN
                        DT_I=LOC_DEN_I*LOC_T_I
                        DT_J=LOC_DEN_J*LOC_T_J
                    ENDIF
                    DO IPHASE = 1, Mdims%nphase
                        UDGI_INT_ALL(:,IPHASE) = (DT_I(IPHASE) * UDGI_ALL(:,IPHASE) + DT_J(IPHASE) * UDGI2_ALL(:, IPHASE)) &
                            / (DT_I(IPHASE) + DT_J(IPHASE))
                        IF( CV_DG_VEL_INT_OPT < 0 ) THEN
                            NDOTQ_INT(IPHASE) = DOT_PRODUCT( CVNORMX_ALL(:, GI), UDGI_INT_ALL( :, IPHASE ) )
                            IF( NDOTQ_INT(IPHASE) <= 0.0 ) THEN  !Incoming
                                !   DT_I=1.0
                                DT_J(IPHASE)=DT_I(IPHASE)+DT_J(IPHASE)
                            ELSE
                                DT_I(IPHASE)=DT_I(IPHASE)+DT_J(IPHASE)
                            !   DT_J=1.0
                            ENDIF
                            UDGI_INT_ALL(:,IPHASE) = (DT_I(IPHASE) * UDGI_ALL(:,IPHASE) + &
                                DT_J(IPHASE) * UDGI2_ALL(:, IPHASE)) / (DT_I(IPHASE) + DT_J(IPHASE))
                        ENDIF
                        UDGI_ALL( :, IPHASE ) = UDGI_INT_ALL( :, IPHASE )
                        UGI_COEF_ELE_ALL(:,IPHASE,:)=DT_I(IPHASE) * UGI_COEF_ELE_ALL(:,IPHASE,:) &
                            /(DT_I(IPHASE) + DT_J(IPHASE))
                        UGI_COEF_ELE2_ALL(:,IPHASE,:)=DT_J(IPHASE) * UGI_COEF_ELE2_ALL(:,IPHASE,:) &
                            /(DT_I(IPHASE) + DT_J(IPHASE))
                    END DO
                ENDIF Conditional_ELE2
            ENDIF Conditional_SELE
            NDOTQ =  MATMUL( CVNORMX_ALL(:, GI), UDGI_ALL)
            ! Define whether flux is incoming or outgoing, depending on direction of flow
            INCOME = 0.5*( 1. + SIGN(1.0, -NDOTQ) )
            !Calculate velocity velocity
            do iphase = 1, Mdims%nphase
                NUGI_ALL(:, IPHASE) = matmul(LOC_NU( :, IPHASE, : ), CV_funs%sufen( :, GI ))
            end do
            IF( between_elements ) THEN
                ! Reduce by half and take the other half from the other side of element...
                do iphase = 1, Mdims%nphase
                    NUGI_ALL(:, IPHASE) = 0.5*NUGI_ALL(:, IPHASE) + 0.5*matmul(LOC2_NU( :, IPHASE, : ), CV_funs%sufen( :, GI ))
                end do
            end if
            ! Calculate NDOTQNEW from NDOTQ
            if (not_OLD_VEL) then
                do iphase = 1, Mdims%nphase
                    NDOTQNEW(iphase) = NDOTQ(iphase) + dot_product(matmul( CVNORMX_ALL(:, GI), UGI_COEF_ELE_ALL(:, iphase,:)*&
                        ( LOC_U(:,iphase,:)-LOC_NU(:,iphase,:))), CV_funs%sufen( :, GI ))
                end do
                IF( between_elements ) THEN
                    ! We have a discontinuity between elements so integrate along the face...
                    DO U_SKLOC = 1, Mdims%u_snloc
                        U_KLOC = U_SLOC2LOC(U_SKLOC)
                        U_KLOC2 = U_OTHER_LOC( U_KLOC )
                        DO IDIM = 1, Mdims%ndim
                            NDOTQNEW=NDOTQNEW + CV_funs%sufen( U_KLOC, GI ) * UGI_COEF_ELE2_ALL(IDIM, :,U_KLOC2) &
                                * ( LOC2_U(IDIM,:, U_KLOC ) - LOC2_NU(IDIM,:,U_KLOC ) ) * CVNORMX_ALL(IDIM, GI)
                        END DO
                    END DO
                END IF
            end if
            RETURN
        END SUBROUTINE GET_INT_VEL_ORIG_NEW

        SUBROUTINE GET_INT_VEL_POROUS_VEL(NDOTQNEW, NDOTQ, INCOME, &
            LOC_T_I, LOC_T_J, LOC_FEMT, &
            LOC_NU, LOC2_NU, SLOC_NU, &
            UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
            I_adv_coef, I_adv_coef_grad, &
            J_adv_coef, J_adv_coef_grad, &
            I_inv_adv_coef, J_inv_adv_coef, &
            UDGI_ALL,MASS_CV_I, MASS_CV_J, &
            TUPWIND_IN, TUPWIND_OUT, &
            not_OLD_VEL, anisotropic_and_frontier)
!sprint_to_do. CHANGE THIS SO WE ONLY PASS DOWN INOUT VARIABLES AND WE DON'T ALLOCATE ANYTHING INSIDE THIS
            !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===
            ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI.
            ! it assumes a compact_overlapping decomposition approach for velocity.
            IMPLICIT NONE
            REAL, DIMENSION( : ), intent( inout ) :: NDOTQNEW, NDOTQ, INCOME
            REAL, DIMENSION( : ), intent( in ) :: LOC_T_I, LOC_T_J
            REAL, DIMENSION( :, : ), intent( in ) :: LOC_FEMT
            REAL, DIMENSION( :, :, : ), intent( in ) ::  LOC_NU, LOC2_NU, SLOC_NU
            REAL, DIMENSION( :, :, : ), intent( inout ) :: UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL
            REAL, DIMENSION( :, :, : ), intent( in ) :: I_adv_coef, I_adv_coef_grad, &
                J_adv_coef, J_adv_coef_grad
            REAL, DIMENSION( :, :, : ), intent( in ) :: I_inv_adv_coef, J_inv_adv_coef
            REAL, DIMENSION( :, :  ), intent( inout ) :: UDGI_ALL
            REAL, intent( in ) :: MASS_CV_I, MASS_CV_J
            REAL, DIMENSION( : ), intent( in ) :: TUPWIND_IN, TUPWIND_OUT!(Mdims%nphase)
            logical, intent(in) :: not_OLD_VEL, anisotropic_and_frontier

            UGI_COEF_ELE_ALL=0.0 ; UGI_COEF_ELE2_ALL=0.0
            Conditional_SELE: IF( on_domain_boundary ) THEN ! On the boundary of the domain.
                !Initialize variables
                if (not_OLD_VEL) then
                    forall (iv_iphase = 1:Mdims%nphase, iv_idim = 1:Mdims%ndim)
                        ROW_SUM_INV_VI(iv_idim,iv_iphase)=SUM(I_inv_adv_coef(iv_idim,:,iv_iphase))
                    end forall
                end if
                DO iv_iphase = 1, Mdims%nphase
                    IPRES = 1 + (iv_iphase-1)/mdims%N_IN_PRES
                    IF( WIC_P_BC_ALL( 1, IPRES, SELE) == WIC_P_BC_DIRICHLET ) THEN ! Pressure boundary condition
                        !(vel * shape_functions)/sigma
                        UDGI_ALL(:, iv_iphase) = matmul(I_inv_adv_coef(:,:,iv_iphase),&
                            matmul(LOC_NU( :, iv_iphase, : ), CV_funs%sufen( :, GI )))
                        ! Here we assume that sigma_out/sigma_in is a diagonal matrix
                        ! which effectively assumes that the anisotropy just inside the domain
                        ! is the same as just outside the domain.
                        ! Multiply by a normalized sigma tensor so that we use the
                        ! sigma from just outside the boundary:
                        DO iv_cv_skloc = 1, Mdims%cv_snloc
                            iv_cv_kloc = CV_SLOC2LOC( iv_cv_skloc )
                            IF(iv_cv_kloc==CV_ILOC) THEN
                                iv_cv_snodk = ( SELE - 1 ) * Mdims%cv_snloc + iv_cv_skloc
                                iv_cv_snodk_IPHA = iv_cv_snodk + ( iv_iphase - 1 ) * Mdims%stotel*Mdims%cv_snloc
                                iv_SUF_SIG_DIAGTEN_BC_GI( 1:Mdims%ndim ) = SUF_SIG_DIAGTEN_BC( iv_cv_snodk_IPHA, 1:Mdims%ndim )
                                exit
                            ENDIF
                        END DO
                        ! Only modify boundary velocity for incoming velocity...
                        iv_incomming_flow = DOT_PRODUCT(UDGI_ALL(:, iv_iphase), CVNORMX_ALL(:, GI)) < 0.0
                        if (not_OLD_VEL) then
                            DO iv_u_kloc = 1, Mdims%u_nloc
!                               IF (.false.) THEN !<= this one for strong boundary conditions
                                IF (iv_incomming_flow) THEN ! Incomming...
                                    UGI_COEF_ELE_ALL(:, iv_iphase, iv_u_kloc)=iv_SUF_SIG_DIAGTEN_BC_GI(:)
                                ELSE
                                    UGI_COEF_ELE_ALL(:, iv_iphase, iv_u_kloc)=1.0
                                ENDIF
                                UGI_COEF_ELE_ALL(:, iv_iphase, iv_u_kloc)= matmul(I_inv_adv_coef(:,:,iv_iphase),UGI_COEF_ELE_ALL(:, iv_iphase, iv_u_kloc))
                            END DO
                        end if

                        if(iv_incomming_flow) UDGI_ALL(:, iv_iphase) = UDGI_ALL(:, iv_iphase) * iv_SUF_SIG_DIAGTEN_BC_GI(:)

                    ELSE ! Specified vel bc.
                        UDGI_ALL(:, iv_iphase) = 0.0
                        UDGI_ALL_FOR_INV(:, iv_iphase) = 0.0
                        UGI_COEF_ELE_ALL(:, iv_iphase, :) = 0.0
                        DO iv_u_skloc = 1, Mdims%u_snloc
                            iv_u_kloc = U_SLOC2LOC( iv_u_skloc )
                            IF(WIC_U_BC_ALL( 1, iv_iphase, SELE ) == 10) THEN
                                UDGI_ALL(:, iv_iphase) = UDGI_ALL(:, iv_iphase) + CV_funs%sufen( iv_u_kloc, GI ) * 0.5 * &
                                    SUF_U_BC_ALL(:, iv_iphase, Mdims%u_snloc* (SELE-1) +iv_u_skloc)
                                UDGI_ALL_FOR_INV(:, iv_iphase) = UDGI_ALL_FOR_INV(:, iv_iphase) + CV_funs%sufen( iv_u_kloc, GI ) * 0.5 * &
                                    SLOC_NU(:, iv_iphase, iv_u_skloc)
                                UGI_COEF_ELE_ALL(:, iv_iphase, iv_u_kloc)=0.5*ROW_SUM_INV_VI(:,iv_iphase)
                            ELSE
                                UDGI_ALL(:, iv_iphase) = UDGI_ALL(:, iv_iphase) + CV_funs%sufen( iv_u_kloc, GI )*SUF_U_BC_ALL(:, iv_iphase, Mdims%u_snloc* (SELE-1) +iv_u_skloc)
                            END IF

                        END DO
                        UDGI_ALL(:, iv_iphase) = UDGI_ALL(:, iv_iphase)  + matmul(I_inv_adv_coef(:,:,iv_iphase),UDGI_ALL_FOR_INV(:, iv_iphase))
                    END IF
                END DO ! PHASE LOOP
            ELSE IF( .not. between_elements) THEN!same element
                !vel(GI) = (vel * shape_functions)/sigma
                do iv_iphase = 1, Mdims%nphase
                    UDGI_ALL(:, iv_iphase) = matmul(I_inv_adv_coef(:,:,iv_iphase),&
                        matmul(LOC_NU( :, iv_iphase, : ), CV_funs%sufen( :, GI )))
                    UDGI2_ALL(:, iv_iphase) = matmul(J_inv_adv_coef(:,:,iv_iphase),&
                        matmul(LOC_NU( :, iv_iphase, : ), CV_funs%sufen( :, GI )))
                end do
                !Get the projected velocity
                NDOTQ  = MATMUL( CVNORMX_ALL(:, GI), UDGI_ALL )
                !Get the direction of the flow
                WHERE (NDOTQ < 0.0)
                    INCOME=1.0
                ELSE WHERE
                    INCOME=0.0
                END WHERE
                !Calculate velocity on the interface, either using upwinding or high order methods
!                if (Mdisopt%in_ele_upwind /= 1) then !high order
                if (use_porous_limiter) then!high order
                    !Calculate saturation at GI, necessary for the limiter
                    FEMTGI_IPHA = matmul(LOC_FEMT, CV_funs%scvfen(:,GI) )
                    ! ************NEW LIMITER**************************
                    XI_LIMIT = 2.0
                     !Call the limiter to obtain the limited saturation value at the interface
                    CALL ONVDLIM_ANO_MANY( Mdims%nphase, LIMT3(:), FEMTGI_IPHA(:), INCOME(:), &
                        LOC_T_I(:), LOC_T_J(:),XI_LIMIT(:), TUPWIND_IN(:), TUPWIND_OUT(:) )
                    !We perform: n' * sigma * n
                    DO iv_iphase = 1, Mdims%nphase
                        ABS_CV_NODI_IPHA(iv_iphase) = dot_product(CVNORMX_ALL(:, GI),matmul(I_adv_coef(:,:,iv_iphase), CVNORMX_ALL(:, GI)))
                        GRAD_ABS_CV_NODI_IPHA(iv_iphase) = dot_product(CVNORMX_ALL(:, GI),matmul(I_adv_coef_grad(:,:,iv_iphase), CVNORMX_ALL(:, GI)))
                        ABS_CV_NODJ_IPHA(iv_iphase) = dot_product(CVNORMX_ALL(:, GI),matmul(J_adv_coef(:,:,iv_iphase), CVNORMX_ALL(:, GI)))
                        GRAD_ABS_CV_NODJ_IPHA(iv_iphase) = dot_product(CVNORMX_ALL(:, GI),matmul(J_adv_coef_grad(:,:,iv_iphase), CVNORMX_ALL(:, GI)))
                    END DO
                    abs_tilde(:) = 0.5*(ABS_CV_NODI_IPHA(:) + ( LIMT3(:) - LOC_T_I(:) ) * GRAD_ABS_CV_NODI_IPHA(:)+&
                        ABS_CV_NODJ_IPHA(:) + ( LIMT3 - LOC_T_J(:) ) * GRAD_ABS_CV_NODJ_IPHA(:) )
                    !Make sure the value of sigma is between bounds
                    abs_tilde = min(max(ABS_CV_NODI_IPHA,  ABS_CV_NODJ_IPHA), &
                        max(min(ABS_CV_NODI_IPHA,  ABS_CV_NODJ_IPHA),  abs_tilde ))
                    !We need the projected velocity from the other node
                    NDOTQ2 = MATMUL( CVNORMX_ALL(:, GI), UDGI2_ALL )
                    !Calculation of the velocity at the interface using the sigma at the interface
                    NDOTQ_TILDE = 0.5*( NDOTQ*ABS_CV_NODI_IPHA + NDOTQ2*ABS_CV_NODJ_IPHA ) /abs_tilde
                    !Calculate the contribution of each side
                    INCOME = MIN(1.0, MAX(0.0, (NDOTQ_TILDE - NDOTQ)/VTOLFUN( NDOTQ2 - NDOTQ ) ))
                end if
                !Finally we calculate the velocity at the interface
                DO iv_iphase = 1, Mdims%nphase
                    !Calculate contributions from each side
                    UDGI_ALL(:, iv_iphase) = UDGI_ALL(:, iv_iphase) * (1.0-INCOME(iv_iphase)) +&
                        UDGI2_ALL(:, iv_iphase) * INCOME(iv_iphase)
                END DO ! PHASE LOOP
                if (not_OLD_VEL) then
                    forall (iv_iphase = 1:Mdims%nphase, iv_idim = 1:Mdims%ndim)
                        ROW_SUM_INV_VI(iv_idim,iv_iphase)=SUM(I_inv_adv_coef(iv_idim,:,iv_iphase))
                        ROW_SUM_INV_VJ(iv_idim,iv_iphase)=SUM(J_inv_adv_coef(iv_idim,:,iv_iphase))
                    end forall
                    DO iv_iphase = 1, Mdims%nphase
                        UGI_COEF_ELE_ALL(:, iv_iphase, :)=SPREAD(ROW_SUM_INV_VI(:,iv_iphase)* (1.0-INCOME(iv_iphase)) &
                            +ROW_SUM_INV_VJ(:,iv_iphase)* INCOME(iv_iphase), DIM=2, NCOPIES=Mdims%u_nloc)
                    END DO
                end if
            ELSE !Method initially coded for between elements. Does not use a TVD limiter but a weighting method
                !vel(GI) = (vel * shape_functions)/sigma
                do iv_iphase = 1, Mdims%nphase
                    !Velocity including sigma
                    UDGI_ALL(:, iv_iphase) = matmul(LOC_NU( :, iv_iphase, : ), CV_funs%sufen( :, GI ))
                    !Normal flow including sigma, to know direction of flow
                    NDOTQ(iv_iphase)  = dot_product( CVNORMX_ALL(:, GI),UDGI_ALL(:, iv_iphase))
                    !Actual advection velocity by removing the contribution of sigma
                    UDGI_ALL(:, iv_iphase) = matmul(I_inv_adv_coef(:,:,iv_iphase),UDGI_ALL(:, iv_iphase))
                end do
                IF( between_elements ) THEN
                    do iv_iphase = 1, Mdims%nphase
                        !Velocity including sigma
                        UDGI2_ALL(:, iv_iphase) = matmul(LOC2_NU( :, iv_iphase, : ), CV_funs%sufen( :, GI ))
                        !Normal flow including sigma, to know direction of flow
                        NDOTQ2(iv_iphase) = dot_product( CVNORMX_ALL(:, GI),UDGI2_ALL(:, iv_iphase))
                        !Actual advection velocity by removing the contribution of sigma
                        UDGI2_ALL(:, iv_iphase) = matmul(J_inv_adv_coef(:,:,iv_iphase),UDGI2_ALL(:, iv_iphase))
                    end do
                else !same element
                    do iv_iphase = 1, Mdims%nphase
                        UDGI2_ALL(:, iv_iphase) = matmul(J_inv_adv_coef(:,:,iv_iphase),&
                            matmul(LOC_NU( :, iv_iphase, : ), CV_funs%sufen( :, GI )))
                    end do
                    NDOTQ2 = NDOTQ
                end if
                !The rest of the method is different depending on the type of permeability
                if ( anisotropic_and_frontier ) then !Fully anisotropic permeability
                    !Create diagonal matrix with iv_ones in the diagonal
                    iv_ones = 0.
                    do iv_idim = 1, Mdims%ndim
                        iv_ones(iv_idim, iv_idim) = 1.0
                    end do
                    !Sigma averaged with the mass to be used as divisor
                    iv_sigma_aver = I_adv_coef*MASS_CV_I+J_adv_coef*MASS_CV_J
                    do iv_iphase = 1, Mdims%nphase
                        call invert(iv_sigma_aver(:,:, iv_iphase))
                        !Calculate the contribution of each side, considering sigma and the volume of the CVs
                        if ( ( NDOTQ(iv_iphase) + NDOTQ2(iv_iphase) ) > 0.0 ) then
                            !We redefine sigma so that it detects oscillations using first order taylor series
                            iv_aux_tensor(:,:,iv_iphase) =  I_adv_coef(:,:,iv_iphase) &
                                + 0.5*( LOC_T_J(iv_iphase) - LOC_T_I(iv_iphase) ) * I_adv_coef_grad(:,:,iv_iphase)
                            !We limit the value
                            iv_aux_tensor(:,:,iv_iphase) = min(1000.*max(I_adv_coef(:,:,iv_iphase),  J_adv_coef(:,:,iv_iphase)), &
                                max(0.001*min(I_adv_coef(:,:,iv_iphase),  J_adv_coef(:,:,iv_iphase)), iv_aux_tensor(:,:,iv_iphase) ))
                            iv_aux_tensor(:,:,iv_iphase)= min( 1.0, matmul(I_inv_adv_coef(:,:,iv_iphase), iv_aux_tensor(:,:,iv_iphase)))
                            !Calculate importance of each side
                            iv_aux_tensor2(:,:,iv_iphase) = matmul(iv_aux_tensor(:,:,iv_iphase), matmul(iv_sigma_aver(:,:,iv_iphase),I_adv_coef(:,:,iv_iphase)*MASS_CV_I ))
                            !iv_aux_tensor2 has to be calculated before since iv_aux_tensor is rewritten!
                            iv_aux_tensor(:,:,iv_iphase) = (iv_ones(:,:)-iv_aux_tensor(:,:,iv_iphase)) + matmul(iv_aux_tensor(:,:,iv_iphase),matmul(iv_sigma_aver(:,:,iv_iphase), J_adv_coef(:,:,iv_iphase)*MASS_CV_J ))
                        else
                            !We redefine sigma so that it detects oscillations using first order taylor series
                            iv_aux_tensor(:,:,iv_iphase) =  J_adv_coef(:,:,iv_iphase) &
                                + 0.5*( LOC_T_I(iv_iphase) - LOC_T_J(iv_iphase) ) * J_adv_coef_grad(:,:,iv_iphase)
                            !We limit the value
                            iv_aux_tensor(:,:,iv_iphase) = min(1000.*max(I_adv_coef(:,:,iv_iphase),  J_adv_coef(:,:,iv_iphase)), &
                                max(0.001*min(I_adv_coef(:,:,iv_iphase),  J_adv_coef(:,:,iv_iphase)), iv_aux_tensor(:,:,iv_iphase) ))
                            iv_aux_tensor(:,:,iv_iphase)= min( 1.0, matmul(J_inv_adv_coef(:,:,iv_iphase), iv_aux_tensor(:,:,iv_iphase)))
                            !Calculate importance of each side
                            iv_aux_tensor2(:,:,iv_iphase) = (iv_ones(:,:)-iv_aux_tensor(:,:,iv_iphase)) + matmul(iv_aux_tensor(:,:,iv_iphase),matmul(iv_sigma_aver(:,:,iv_iphase), I_adv_coef(:,:,iv_iphase)*MASS_CV_I ))
                            !iv_aux_tensor2 has to be calculated before since iv_aux_tensor is rewritten!
                            iv_aux_tensor(:,:,iv_iphase) = matmul(iv_aux_tensor(:,:,iv_iphase), matmul(iv_sigma_aver(:,:,iv_iphase),J_adv_coef(:,:,iv_iphase)*MASS_CV_J ))
                        end if
                        !Calculation of the velocity at the GI point
                        UDGI_ALL(:, iv_iphase) = matmul(iv_aux_tensor(:,:,iv_iphase), UDGI_ALL(:, iv_iphase)) + matmul( iv_aux_tensor2(:,:,iv_iphase),UDGI2_ALL(:, iv_iphase))
                    end do
                    !Calculation of the coefficients at the GI point
                    if (not_OLD_VEL) then
                        forall (iv_iphase = 1:Mdims%nphase, iv_idim = 1:Mdims%ndim)
                            ROW_SUM_INV_VI(iv_idim,iv_iphase)=SUM(I_inv_adv_coef(iv_idim,:,iv_iphase))
                            ROW_SUM_INV_VJ(iv_idim,iv_iphase)=SUM(J_inv_adv_coef(iv_idim,:,iv_iphase))
                        end forall
                        IF( between_elements ) then
                            DO iv_iphase = 1, Mdims%nphase
                                UGI_COEF_ELE_ALL(:, iv_iphase, :) = matmul(iv_aux_tensor(:,:,iv_iphase), SPREAD(ROW_SUM_INV_VI(:,iv_iphase), DIM=2, NCOPIES=Mdims%u_nloc))
                                UGI_COEF_ELE2_ALL(:, iv_iphase, :) = matmul(iv_aux_tensor2(:,:,iv_iphase),SPREAD(ROW_SUM_INV_VJ(:,iv_iphase), DIM=2, NCOPIES=Mdims%u_nloc))
                            END DO
                        else !same element
                            DO iv_iphase = 1, Mdims%nphase
                                UGI_COEF_ELE_ALL(:, iv_iphase, :) = matmul(iv_aux_tensor(:,:,iv_iphase), SPREAD(ROW_SUM_INV_VI(:,iv_iphase), DIM=2, NCOPIES=Mdims%u_nloc)) +&
                                    matmul(iv_aux_tensor2(:,:,iv_iphase),SPREAD(ROW_SUM_INV_VJ(:,iv_iphase), DIM=2, NCOPIES=Mdims%u_nloc))
                            END DO
                        end if
                    end if
                ELSE !DG, not a CV neighbouring another CV with different permeability
                    !CV_DG_VEL_INT_OPT <= parameter to choose different options for DG
                    !We perform: n' * sigma * n
                    DO iv_iphase = 1, Mdims%nphase
                        ABS_CV_NODI_IPHA(iv_iphase) = dot_product(CVNORMX_ALL(:, GI),matmul(I_adv_coef(:,:,iv_iphase), CVNORMX_ALL(:, GI)))
                        GRAD_ABS_CV_NODI_IPHA(iv_iphase) = dot_product(CVNORMX_ALL(:, GI),matmul(I_adv_coef_grad(:,:,iv_iphase), CVNORMX_ALL(:, GI)))
                        ABS_CV_NODJ_IPHA(iv_iphase) = dot_product(CVNORMX_ALL(:, GI),matmul(J_adv_coef(:,:,iv_iphase), CVNORMX_ALL(:, GI)))
                        GRAD_ABS_CV_NODJ_IPHA(iv_iphase) = dot_product(CVNORMX_ALL(:, GI),matmul(J_adv_coef_grad(:,:,iv_iphase), CVNORMX_ALL(:, GI)))
                        !Axuliar variable to reduce computations, LIMT3 was unused for this part of the code
                        LIMT3(iv_iphase) = ABS_CV_NODI_IPHA(iv_iphase)*MASS_CV_I+ABS_CV_NODJ_IPHA(iv_iphase)*MASS_CV_J
                        !Calculate the contribution of each side, considering sigma and the volume of the CVs
                        if (( NDOTQ(iv_iphase) + NDOTQ2(iv_iphase)) > 0.0 ) then
                            !We redefine sigma so that it detects oscillations using first order taylor series
                            abs_tilde(iv_iphase) =  ABS_CV_NODI_IPHA(iv_iphase)  + 0.5*( LOC_T_J(iv_iphase) - LOC_T_I(iv_iphase) ) * GRAD_ABS_CV_NODI_IPHA(iv_iphase)
                            !We limit the value
                            abs_tilde(iv_iphase) = min(1000.*max(ABS_CV_NODI_IPHA(iv_iphase),  ABS_CV_NODJ_IPHA(iv_iphase)), &
                                max(0.001*min(ABS_CV_NODI_IPHA(iv_iphase),  ABS_CV_NODJ_IPHA(iv_iphase)), abs_tilde(iv_iphase) ))
                            wrelax(iv_iphase)= min( 1.0, abs_tilde(iv_iphase)/ABS_CV_NODI_IPHA(iv_iphase) )
                            !Calculate importance of each side
                            DT_I(iv_iphase) = (1.-wrelax(iv_iphase)) + wrelax(iv_iphase)*ABS_CV_NODJ_IPHA(iv_iphase)*MASS_CV_J /LIMT3(iv_iphase)
                            DT_J(iv_iphase) = wrelax(iv_iphase)*ABS_CV_NODI_IPHA(iv_iphase)*MASS_CV_I /LIMT3(iv_iphase)
                        else
                            !We redefine sigma so that it detects oscillations using first order taylor series
                            abs_tilde(iv_iphase) =  ABS_CV_NODJ_IPHA(iv_iphase)  + 0.5*( LOC_T_I(iv_iphase) - LOC_T_J(iv_iphase) ) * GRAD_ABS_CV_NODJ_IPHA(iv_iphase)
                            !We limit the value
                            abs_tilde(iv_iphase) = min(1000.*max(ABS_CV_NODI_IPHA(iv_iphase),  ABS_CV_NODJ_IPHA(iv_iphase)), &
                                max(0.001*min(ABS_CV_NODI_IPHA(iv_iphase),  ABS_CV_NODJ_IPHA(iv_iphase)), abs_tilde(iv_iphase) ))
                            wrelax(iv_iphase)= min( 1.0, abs_tilde(iv_iphase)/ABS_CV_NODJ_IPHA(iv_iphase)  )
                            !Calculate importance of each side
                            DT_I(iv_iphase) = wrelax(iv_iphase)*ABS_CV_NODJ_IPHA(iv_iphase)*MASS_CV_J /LIMT3(iv_iphase)
                            DT_J(iv_iphase) = (1.-wrelax(iv_iphase)) + wrelax(iv_iphase)*ABS_CV_NODI_IPHA(iv_iphase)*MASS_CV_I /LIMT3(iv_iphase)
                        end if
                        !Calculation of the velocity at the GI point
                        UDGI_ALL(:, iv_iphase) = DT_I(iv_iphase) * UDGI_ALL(:, iv_iphase) + DT_J(iv_iphase) * UDGI2_ALL(:, iv_iphase)
                    end do
                    !Calculation of the coefficients at the GI point
                    if (not_OLD_VEL) then
                        forall (iv_iphase = 1:Mdims%nphase, iv_idim = 1:Mdims%ndim)
                            ROW_SUM_INV_VI(iv_idim,iv_iphase)=SUM(I_inv_adv_coef(iv_idim,:,iv_iphase))
                            ROW_SUM_INV_VJ(iv_idim,iv_iphase)=SUM(J_inv_adv_coef(iv_idim,:,iv_iphase))
                        end forall
                        IF( between_elements ) then
                            DO iv_iphase = 1, Mdims%nphase
                                UGI_COEF_ELE_ALL(:, iv_iphase, :) = DT_I(iv_iphase) * SPREAD(ROW_SUM_INV_VI(:,iv_iphase), DIM=2, NCOPIES=Mdims%u_nloc)
                                UGI_COEF_ELE2_ALL(:, iv_iphase, :) = DT_J(iv_iphase) * SPREAD(ROW_SUM_INV_VJ(:,iv_iphase), DIM=2, NCOPIES=Mdims%u_nloc)
                            END DO
                        else !same element
                            DO iv_iphase = 1, Mdims%nphase
                                UGI_COEF_ELE_ALL(:, iv_iphase, :) = DT_I(iv_iphase) * SPREAD(ROW_SUM_INV_VI(:,iv_iphase), DIM=2, NCOPIES=Mdims%u_nloc) +&
                                    DT_J(iv_iphase) * SPREAD(ROW_SUM_INV_VJ(:,iv_iphase), DIM=2, NCOPIES=Mdims%u_nloc)
                            END DO
                        end if
                    end if
                end if
            END IF Conditional_SELE
            ! Define whether flux is incoming or outgoing, depending on direction of flow
            NDOTQ =  MATMUL( CVNORMX_ALL(:, GI), UDGI_ALL )
            WHERE ( NDOTQ >= 0. )
                INCOME = 0.
            ELSE WHERE
                INCOME = 1.
            END WHERE
            ! Calculate NDOTQNEW from NDOTQ
            if (not_OLD_VEL) then
                do iv_iphase = 1, Mdims%nphase
                    NDOTQNEW(iv_iphase) = NDOTQ(iv_iphase) + dot_product(matmul( CVNORMX_ALL(:, GI), UGI_COEF_ELE_ALL(:, iv_iphase,:)*&
                        ( LOC_U(:,iv_iphase,:)-LOC_NU(:,iv_iphase,:))), CV_funs%sufen( :, GI ))
                end do
                IF( between_elements) THEN
                    ! We have a discontinuity between elements so integrate along the face...
                    DO iv_u_skloc = 1, Mdims%u_snloc
                        iv_u_kloc = U_SLOC2LOC(iv_u_skloc)
                        iv_u_kloc2 = U_OTHER_LOC( iv_u_kloc )
                        DO iv_idim = 1, Mdims%ndim
                            NDOTQNEW=NDOTQNEW + CV_funs%sufen( iv_u_kloc, GI ) * UGI_COEF_ELE2_ALL(iv_idim, :,iv_u_kloc2) &
                                * ( LOC2_U(iv_idim,:, iv_u_kloc ) - LOC2_NU(iv_idim,:,iv_u_kloc ) ) * CVNORMX_ALL(iv_idim, GI)
                        END DO
                    END DO
                END IF
            end if

            RETURN
        END SUBROUTINE GET_INT_VEL_POROUS_VEL
        PURE SUBROUTINE ONVDLIM_ANO_MANY( NFIELD, &
            TDLIM, TDCEN, INCOME, &
            ETDNEW_PELE, ETDNEW_PELEOT, XI_LIMIT,  &
            TUPWIN, TUPWI2 )
            implicit none
            ! This sub calculates the limited face values TDADJ(1...SNGI) from the central
            ! difference face values TDCEN(1...SNGI) using a NVD shceme.
            ! INCOME(1...SNGI)=1 for incomming to element ELE  else =0.
            ! LIBETA is the flux limiting parameter.
            ! TDMAX(PELE)=maximum of the surrounding 6 element values of element PELE.
            ! TDMIN(PELE)=minimum of the surrounding 6 element values of element PELE.
            ! PELEOT=element at other side of current face.
            ! ELEOT2=element at other side of the element ELEOTH.
            ! ELESID=element next to oposing current face.
            ! The elements are arranged in this order: ELEOT2,ELE, PELEOT, ELESID.
            ! This sub finds the neighbouring elements. Suppose that this is the face IFACE.
            !---------------------------------------------------
            !|   ELEOT2   |   ELEOTH   |   ELE     |   ELESID   |
            !---------------------------------------------------
            ! TAIN         THALF       TAOUT
            !---------------------------------------------------
            !>TEXTIN
            !TEXOUT<
            !---------------------------------------------------
            INTEGER, intent( in ) :: NFIELD
            REAL, DIMENSION( NFIELD ), intent( inout ) :: TDLIM
            REAL, DIMENSION( NFIELD ), intent( in ) :: TDCEN, INCOME, XI_LIMIT, TUPWIN, TUPWI2
            REAL, DIMENSION( NFIELD ), intent( in ) :: ETDNEW_PELE, ETDNEW_PELEOT
            ! Local variables
            REAL, PARAMETER :: TOLER=1.0E-10
            REAL :: DENOIN(NFIELD), CTILIN(NFIELD), DENOOU(NFIELD), &
                CTILOU(NFIELD), FTILIN(NFIELD), FTILOU(NFIELD)
            ! Calculate normalisation parameters for incomming velocities
            DENOIN = ( ETDNEW_PELE - TUPWIN )
            where( ABS( DENOIN ) < TOLER )
                DENOIN = SIGN( TOLER, DENOIN )
            end where
            CTILIN = ( ETDNEW_PELEOT - TUPWIN ) / DENOIN
            ! Calculate normalisation parameters for out going velocities
            DENOOU = ( ETDNEW_PELEOT - TUPWI2 )
            where( ABS( DENOOU ) < TOLER )
                DENOOU = SIGN( TOLER, DENOOU )
            end where
            CTILOU = ( ETDNEW_PELE - TUPWI2 ) / DENOOU
            FTILIN = ( TDCEN - TUPWIN ) / DENOIN
            FTILOU = ( TDCEN - TUPWI2 ) / DENOOU
            ! Velocity is going out of element
            TDLIM =        INCOME   * ( TUPWIN + NVDFUNNEW_MANY( FTILIN, CTILIN, XI_LIMIT ) * DENOIN ) &
                + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW_MANY( FTILOU, CTILOU, XI_LIMIT ) * DENOOU )
            TDLIM = MAX( TDLIM, 0.0 )
            RETURN
        END SUBROUTINE ONVDLIM_ANO_MANY

        FUNCTION FACE_THETA_MANY( DT, CV_THETA, INTERFACE_TRACK, HDC, NPHASE, n_in_pres,&
            NDOTQ, LIMDT, DIFF_COEF_DIVDX, &
            T_NODJ_IPHA, T_NODI_IPHA,  &
            NDOTQOLD, LIMDTOLD, DIFF_COEFOLD_DIVDX, TOLD_NODJ_IPHA, TOLD_NODI_IPHA )
            IMPLICIT NONE
            ! Define face value of theta
            INTEGER, intent(in) :: NPHASE, n_in_pres
            REAL, intent(in) :: DT, CV_THETA, HDC
            REAL, DIMENSION( NPHASE ), intent(in) :: NDOTQ, LIMDT, DIFF_COEF_DIVDX, T_NODJ_IPHA, T_NODI_IPHA,  &
                NDOTQOLD, LIMDTOLD, DIFF_COEFOLD_DIVDX, TOLD_NODJ_IPHA, TOLD_NODI_IPHA
            real, dimension(NPHASE) :: FACE_THETA_MANY
            LOGICAL, intent(in) :: INTERFACE_TRACK
            ! Local variables
            REAL, DIMENSION( : ), allocatable :: HF, HFOLD, GF, PINVTH, QINVTH
            INTEGER :: IPHASE

            IF( CV_THETA >= 0.0) THEN ! Specified
                FACE_THETA_MANY = CV_THETA
            ELSE ! Non-linear
                ALLOCATE( HF( NPHASE ), HFOLD( NPHASE ), GF( NPHASE ), PINVTH( NPHASE ), QINVTH( NPHASE ) )
                HF    = NDOTQ * LIMDT + DIFF_COEF_DIVDX * ( T_NODI_IPHA - T_NODJ_IPHA )
                HFOLD = NDOTQOLD * LIMDTOLD + DIFF_COEFOLD_DIVDX * ( TOLD_NODI_IPHA - TOLD_NODJ_IPHA )
                GF = TOLFUN_MANY( DT *  (HF - HFOLD) )
                PINVTH = HDC * ( T_NODI_IPHA - TOLD_NODI_IPHA ) / GF
                QINVTH = HDC * ( T_NODJ_IPHA - TOLD_NODJ_IPHA ) / GF
                ! 0.5 is the original value.
                !       FTHETA = MAX( 0.5, 1. - 0.5 * MIN( ABS( PINVTH ), ABS( QINVTH )))
                !       FTHETA = MAX( 0.5, 1. - 0.25 * MIN( ABS( PINVTH ), ABS( QINVTH )))
                IF(INTERFACE_TRACK) THEN ! For interface tracking use forward Euler as much as possible...
                    !            FTHETA = MAX( 0.0, 1. - 0.125 * MIN( ABS( PINVTH ), ABS( QINVTH )))
                    DO IPHASE=1,NPHASE
                        FACE_THETA_MANY(IPHASE) = MAX( 0.0, 1. - 0.5 * MIN( ABS( PINVTH(IPHASE) ), ABS( QINVTH(IPHASE) )))
                    END DO
                ELSE ! for Crank Nickolson time stepping base scheme...
                    DO IPHASE=1,n_in_pres!with no wells this is the same to nphase
                        FACE_THETA_MANY(IPHASE) = MAX( 0.5, 1. - 0.125 * MIN( ABS( PINVTH(IPHASE) ), ABS( QINVTH(IPHASE) )))
                    END DO
                    if (nphase /= n_in_pres) then!for wells we impose implicit euler, the Courant number is massive anyway...
                        FACE_THETA_MANY(n_in_pres + 1 : nphase) = 1.0!<=backward euler
                    end if
                ENDIF
            ENDIF
            RETURN

        END FUNCTION FACE_THETA_MANY

        SUBROUTINE SCVDETNX_new( ELE,GI,SCVDETWEI, CVNORMX_ALL,XC_ALL, X_NOD)
            !     --------------------------------------------------
            !
            !     - this subroutine calculates the control volume (CV)
            !     - CVNORMX, CVNORMY, CVNORMZ normals at the Gaussian
            !     - integration points GI. NODI = the current global
            !     - node number for the co-ordinates.
            !     - (XC,YC,ZC) is the centre of CV NODI
            !
            !     -------------------------------
            !     - date last modified : 15/03/2003
            !     -------------------------------
            IMPLICIT NONE
            INTEGER, intent( in ) :: ELE, GI, X_NOD
            REAL, DIMENSION( Mdims%ndim ), intent( in ) ::   XC_ALL
            REAL, DIMENSION( Mdims%ndim, CV_GIdims%scvngi ), intent( inout ) :: CVNORMX_ALL
            REAL, DIMENSION( : ), intent( inout ) :: SCVDETWEI
            !     - Local variables
            INTEGER :: NODJ,  JLOC
            REAL :: A, B, C
            REAL :: DETJ
            REAL :: DXDLX, DXDLY, DYDLX
            REAL :: DYDLY, DZDLX, DZDLY
            REAL :: TWOPI
            REAL, PARAMETER :: PI = 3.14159265
            REAL :: RGI, RDUM
            real, dimension(3) :: POSVGI
            !ewrite(3,*)' In SCVDETNX'
            POSVGI = 0.0
            Conditional_Dimension: IF( Mdims%ndim == 3 ) THEN

                DXDLX = 0.0;DXDLY = 0.0
                DYDLX = 0.0;DYDLY = 0.0
                DZDLX = 0.0;DZDLY = 0.0
                do  JLOC = 1, Mdims%x_nloc

                    NODJ = ndgln%x((ELE-1)*Mdims%x_nloc+JLOC)

                    DXDLX = DXDLX + CV_funs%scvfenslx(JLOC,GI)*X_ALL(1,NODJ)
                    DXDLY = DXDLY + CV_funs%scvfensly(JLOC,GI)*X_ALL(1,NODJ)
                    DYDLX = DYDLX + CV_funs%scvfenslx(JLOC,GI)*X_ALL(2,NODJ)
                    DYDLY = DYDLY + CV_funs%scvfensly(JLOC,GI)*X_ALL(2,NODJ)
                    DZDLX = DZDLX + CV_funs%scvfenslx(JLOC,GI)*X_ALL(3,NODJ)
                    DZDLY = DZDLY + CV_funs%scvfensly(JLOC,GI)*X_ALL(3,NODJ)

                    POSVGI = POSVGI + CV_funs%scvfen(JLOC,GI)*X_ALL(:,NODJ)
                end do


                !To calculate the sign of the normal an average between the center of the continuous CV and the center of mass is used
                !this is required as the center of mass has shown not to be reliable and the center of the continuous CV is a particular point that can lead
                !to failures to obtain the sign (perpendicular vectors in a flat boundary)
                POSVGI = POSVGI - (0.8*X_ALL(1:Mdims%ndim, X_NOD) + 0.2*XC_ALL(1:Mdims%ndim))

                CALL NORMGI( CVNORMX_ALL(1,GI), CVNORMX_ALL(2,GI), CVNORMX_ALL(3,GI),&
                    DXDLX,       DYDLX,       DZDLX, &
                    DXDLY,       DYDLY,       DZDLY,&
                    POSVGI(1),     POSVGI(2),     POSVGI(3) )


                A = DYDLX*DZDLY - DYDLY*DZDLX
                B = DXDLX*DZDLY - DXDLY*DZDLX
                C = DXDLX*DYDLY - DXDLY*DYDLX
                !
                !     - Calculate the determinant of the Jacobian at Gauss pnt GI.
                !
                DETJ = SQRT( A**2 + B**2 + C**2 )
                !
                !     - Calculate the determinant times the surface weight at Gauss pnt GI.
                !
                SCVDETWEI(GI) = DETJ*CV_funs%scvfeweigh(GI)
                !
                !     - Calculate the normal at the Gauss pts
                !     - TANX1 = DXDLX, TANY1 = DYDLX, TANZ1 = DZDLX,
                !     - TANX2 = DXDLY, TANY2 = DYDLY, TANZ2 = DZDLY
                !     - Perform cross-product. N = T1 x T2
                !



            ELSE IF(Mdims%ndim == 2) THEN

                TWOPI = 1.0

                RGI   = 0.0
                DXDLX = 0.0;DXDLY = 0.0
                DYDLX = 0.0;DYDLY = 0.0
                DZDLX = 0.0
                !
                !     - Note that we set the derivative wrt to y of coordinate z to 1.0
                !
                DZDLY = 1.0

                do  JLOC = 1, Mdims%x_nloc! Was loop 300

                    NODJ = ndgln%x((ELE-1)*Mdims%x_nloc+JLOC)

                    DXDLX = DXDLX + CV_funs%scvfenslx(JLOC,GI)*X_ALL(1,NODJ)
                    DYDLX = DYDLX + CV_funs%scvfenslx(JLOC,GI)*X_ALL(2,NODJ)

                    POSVGI(1:Mdims%ndim) = POSVGI(1:Mdims%ndim) + CV_funs%scvfen(JLOC,GI)*X_ALL(1:Mdims%ndim,NODJ)

                    RGI = RGI + CV_funs%scvfen(JLOC,GI)*X_ALL(2,NODJ)

                end do ! Was loop 300
                !To calculate the sign of the normal an average between the center of the COntinuous CV and the center of mass is used
                POSVGI(1:Mdims%ndim) = POSVGI(1:Mdims%ndim) - (0.8*X_ALL(1:Mdims%ndim, X_NOD) + 0.2*XC_ALL(1:Mdims%ndim))

                RGI = 1.0

                DETJ = SQRT( DXDLX**2 + DYDLX**2 )
                SCVDETWEI(GI)  = TWOPI*RGI*DETJ*CV_funs%scvfeweigh(GI)
                !
                !     - Calculate the normal at the Gauss pts
                !     - TANX1 = DXDLX, TANY1 = DYDLX, TANZ1 = DZDLX,
                !     - TANX2 = DXDLY, TANY2 = DYDLY, TANZ2 = DZDLY
                !     - Perform cross-product. N = T1 x T2
                !
                CALL NORMGI( CVNORMX_ALL(1,GI), CVNORMX_ALL(2,GI), RDUM,&
                    DXDLX,       DYDLX,       DZDLX, &
                    DXDLY,       DYDLY,       DZDLY,&
                    POSVGI(1),     POSVGI(2),     POSVGI(3) )

            ELSE
                ! For 1D...
                do  JLOC = 1, Mdims%x_nloc! Was loop 300

                    NODJ = ndgln%x((ELE-1)*Mdims%x_nloc+JLOC)

                    POSVGI(1) = POSVGI(1) + CV_funs%scvfen(JLOC,GI)*X_ALL(1,NODJ)

                end do ! Was loop 300
                !
                !     - Note that POSVGIX and POSVGIY can be considered as the components
                !     - of the Gauss pnt GI with the co-ordinate origin positioned at the
                !     - current control volume NODI.
                !
                POSVGI(1) = POSVGI(1) - XC_ALL(1)
                ! SIGN(A,B) sign of B times A.
                CVNORMX_ALL(1,GI) = SIGN( 1.0, POSVGI(1) )

                DETJ = 1.0
                SCVDETWEI(GI)  = DETJ*CV_funs%scvfeweigh(GI)


            ENDIF Conditional_Dimension

        END SUBROUTINE SCVDETNX_new

        subroutine get_neigbouring_lists(JCOUNT_KLOC, ICOUNT_KLOC, JCOUNT_KLOC2 ,ICOUNT_KLOC2,&
                                        C_JCOUNT_KLOC, C_ICOUNT_KLOC, C_JCOUNT_KLOC2, C_ICOUNT_KLOC2 )
        implicit none
            integer, dimension(:), intent(inout):: JCOUNT_KLOC, ICOUNT_KLOC, JCOUNT_KLOC2 ,ICOUNT_KLOC2
            integer, dimension(:), intent(inout):: C_JCOUNT_KLOC, C_ICOUNT_KLOC, C_JCOUNT_KLOC2, C_ICOUNT_KLOC2
            !Local variables
            integer :: U_KLOC, U_NODK, JCOUNT, COUNT, ICOUNT
          ! could retrieve JCOUNT_KLOC and ICOUNT_KLOC from storage depending on quadrature point GLOBAL_FACE
            DO U_KLOC = 1, Mdims%u_nloc
                U_NODK = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_KLOC )
                JCOUNT = 0
                DO COUNT = Mspars%CT%fin( CV_NODI ), Mspars%CT%fin( CV_NODI + 1 ) - 1
                    IF ( Mspars%CT%col( COUNT ) == U_NODK ) THEN
                        JCOUNT = COUNT
                        EXIT
                    END IF
                END DO
                JCOUNT_KLOC( U_KLOC ) = JCOUNT
                if(integrate_other_side) then
                    ! for integrating just on one side...
                    ICOUNT = 0
                    DO COUNT = Mspars%CT%fin( CV_NODJ ), Mspars%CT%fin( CV_NODJ + 1 ) - 1
                        IF ( Mspars%CT%col( COUNT ) == U_NODK ) THEN
                            ICOUNT = COUNT
                            EXIT
                        END IF
                    END DO
                    ICOUNT_KLOC( U_KLOC ) = ICOUNT
                endif
            END DO
            IF ( between_elements ) THEN
                DO U_KLOC =  1, Mdims%u_nloc
                    U_NODK = ndgln%u( ( ELE2 - 1 ) * Mdims%u_nloc + U_KLOC )
                    JCOUNT = 0
                    DO COUNT = Mspars%CT%fin( CV_NODI ), Mspars%CT%fin( CV_NODI + 1 ) - 1
                        IF ( Mspars%CT%col( COUNT ) == U_NODK ) THEN
                            JCOUNT = COUNT
                            EXIT
                        END IF
                    END DO
                    JCOUNT_KLOC2( U_KLOC ) = JCOUNT
                    if(integrate_other_side) then
                        ! for integrating just on one side...
                        ICOUNT = 0
                        DO COUNT = Mspars%CT%fin( CV_NODJ ), Mspars%CT%fin( CV_NODJ + 1 ) - 1
                            IF ( Mspars%CT%col( COUNT ) == U_NODK ) THEN
                                ICOUNT = COUNT
                                EXIT
                            END IF
                        END DO
                        ICOUNT_KLOC2( U_KLOC ) = ICOUNT
                    endif
                END DO
            END IF ! endof IF ( between_elements ) THEN
            IF(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV) THEN
                ! could retrieve JCOUNT_KLOC and ICOUNT_KLOC from storage depending on quadrature point GLOBAL_FACE
                DO U_KLOC = 1, Mdims%u_nloc
                    U_NODK = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_KLOC )
                    JCOUNT = 0
                    DO COUNT = Mspars%C%fin( U_NODK ), Mspars%C%fin( U_NODK + 1 ) - 1
                        IF ( Mspars%C%col( COUNT ) == CV_NODI ) THEN
                            JCOUNT = COUNT
                            EXIT
                        END IF
                    END DO
                    C_JCOUNT_KLOC( U_KLOC ) = JCOUNT
                    if(integrate_other_side) then
                        ! for integrating just on one side...
                        ICOUNT = 0
                        DO COUNT = Mspars%C%fin( U_NODK ), Mspars%C%fin( U_NODK + 1 ) - 1
                            IF ( Mspars%C%col( COUNT ) == CV_NODJ ) THEN
                                ICOUNT = COUNT
                                EXIT
                            END IF
                        END DO
                        C_ICOUNT_KLOC( U_KLOC ) = ICOUNT
                    endif
                END DO
                IF ( between_elements ) THEN
                    DO U_KLOC =  1, Mdims%u_nloc
                        U_NODK = ndgln%u( ( ELE2 - 1 ) * Mdims%u_nloc + U_KLOC )
                        JCOUNT = 0
                        DO COUNT = Mspars%C%fin( U_NODK ), Mspars%C%fin( U_NODK + 1 ) - 1
                            IF ( Mspars%C%col( COUNT ) == CV_NODI ) THEN
                                JCOUNT = COUNT
                                EXIT
                            END IF
                        END DO
                        C_JCOUNT_KLOC2( U_KLOC ) = JCOUNT
                        if(integrate_other_side) then
                            ! for integrating just on one side...
                            ICOUNT = 0
                            DO COUNT = Mspars%C%fin( U_NODK ), Mspars%C%fin( U_NODK + 1 ) - 1
                                IF ( Mspars%C%col( COUNT ) == CV_NODJ ) THEN
                                    ICOUNT = COUNT
                                    EXIT
                                END IF
                            END DO
                            C_ICOUNT_KLOC2( U_KLOC ) = ICOUNT
                        endif
                    END DO
                END IF ! endof IF ( between_elements ) THEN
            ENDIF ! ENDOF IF(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV) THEN
        end subroutine get_neigbouring_lists


        subroutine mass_conservation_check_and_outfluxes(calculate_mass_delta, outfluxes, flag)
            ! Subroutine to calculate the integrated flux across a boundary with the specified surface_ids given that the massflux has been already stored elsewhere
            !also used to calculate mass conservation
            implicit none
            integer, intent(in) :: flag
            type (multi_outfluxes), intent(inout) :: outfluxes
            real, dimension(:,:), intent(inout) :: calculate_mass_delta
            !Local variables
            integer :: iphase, k
            real :: tmp1, tmp2, tmp3
            type(vector_field), pointer :: Por
            ! After looping over all elements, calculate the mass change inside the domain normalised to the mass inside the domain at t=t-1
            ! Difference in Total mass

            select case (flag)
                case (1)
                    !Initialise values
                    if (first_nonlinear_time_step ) then
                        calculate_mass_delta(:,1) = 0.0 ! reinitialise
                        call calculate_internal_volume( packed_state, Mdims, Mass_ELE, &
                            calculate_mass_delta(1:Mdims%n_in_pres,1) , ndgln%cv)
                        !DISABLED AS IT DOES NOT WORK WELL AND IT DOES ACCOUNT FOR A VERY TINY FRACTION OF THE OVERALL MASS
!                        if (Mdims%npres >1)then!consider as well the pipes
!                            call calculate_internal_volume( packed_state, Mdims, pipes_aux%MASS_PIPE, &
!                                calculate_mass_delta(:,1) , ndgln%cv, eles_with_pipe)
!                        end if
                    endif
                    if (outfluxes%calculate_flux) then
                        ! Extract the Porosity
                        Por => extract_vector_field( packed_state, "Porosity" )

                        ! Calculate Pore volume
                        outfluxes%porevolume = 0.0
                        DO ELE = 1, Mdims%totele
                            if (element_owned(tracer, ele)) then
                                outfluxes%porevolume = outfluxes%porevolume + MASS_ELE(ELE) * Por%val(1,ELE)
                            end if
                        END DO

                        call allsum(outfluxes%porevolume) ! Now sum the value over all processors

                    end if
                case default!Now calculate mass conservation
                    !Calculate internal volumes of each phase
                    call calculate_internal_volume( packed_state, Mdims, Mass_ELE, &
                        calculate_mass_internal(1:Mdims%n_in_pres) , ndgln%cv)
                    !DISABLED AS IT DOES NOT WORK WELL AND IT DOES ACCOUNT FOR A VERY TINY FRACTION OF THE OVERALL MASS
!                    if (Mdims%npres >1) then!consider as well the pipes
!                        call calculate_internal_volume( packed_state, Mdims, pipes_aux%MASS_PIPE, &
!                            calculate_mass_internal(:) , ndgln%cv, eles_with_pipe)
!                    end if

                    !Loop over nphases - 1
                    calculate_mass_delta(1,2) = 0.
                    do iphase = 1, max(Mdims%n_in_pres -1, 1)
                        tmp1 = calculate_mass_internal(iphase)!<=volume phase i inside the domain
                        tmp2 = calculate_mass_delta(iphase,1)!<= volume phase i inside the domain at the beginning of the time-level
                        tmp3 = sum(bcs_outfluxes(iphase, :, 0))*dt!<= volume phase i across al the boundaries
                        !Consider also wells
                        if (Mdims%npres > 1) then
                            tmp1 = tmp1 + calculate_mass_internal(iphase+ Mdims%n_in_pres)
                            tmp2 = tmp2 + calculate_mass_delta(iphase+ Mdims%n_in_pres,1)
                            tmp3 = tmp3 + sum(bcs_outfluxes(iphase + Mdims%n_in_pres, :, 0))*dt
                        end if
                        if (isparallel()) then
                            call allsum(tmp1)
                            call allsum(tmp2)
                            call allsum(tmp3)
                        end if

                        !Calculate possible mass creation inside the code
                        if (tmp2 > 1d-8) then
                            calculate_mass_delta(1,2) = max(calculate_mass_delta(1,2), abs( tmp1 - tmp2 + tmp3 ) / tmp2)
                        else
                            calculate_mass_delta(1,2) = 0.
                        end if
                    end do

                    !If calculate outfluxes then do it now
                    if (outfluxes%calculate_flux) then
                        !Retrieve only the values of bcs_outfluxes we are interested in
                        do k = 1, size(outfluxes%outlet_id)
                            do iphase = 1, Mdims%nphase
                                outfluxes%totout(1, iphase, k) = sum(bcs_outfluxes( iphase, :, k))
                            end do
                        end do
                        ! Having finished loop over elements etc. Pass the total flux across all boundaries to the global variable totout
                        if (isParallel()) then
                            do k = 1,size(outfluxes%outlet_id)
                                ! Ensure all processors have the correct value of totout for parallel runs
                                do iphase = 1, Mdims%nphase
                                    call allsum(outfluxes%totout(1, iphase, k))
                                    if (has_temperature) call allsum(outfluxes%totout(2, iphase, k))
                                end do
                            end do
                        end if
                    end if
            end select


        end subroutine mass_conservation_check_and_outfluxes






    END SUBROUTINE CV_ASSEMB


    SUBROUTINE IS_FIELD_CONSTANT(IGOT_T_CONST, IGOT_T_CONST_VALUE, T_ALL, CV_NONODS)
        LOGICAL IGOT_T_CONST
        REAL IGOT_T_CONST_VALUE
        INTEGER CV_NONODS
        REAL T_ALL(CV_NONODS)
        real, parameter :: tolerance = 1.e-10, tolerance2 = 1.e-5
        ! Local variables...
        REAL RMAX,RMIN,PERCENT
        INTEGER I

        RMIN=+1.e+20
        RMAX=-1.e+20
        DO I=1,CV_NONODS
            RMAX=MAX(RMAX, T_ALL(I))
            RMIN=MIN(RMIN, T_ALL(I))
        END DO

        PERCENT = ABS(RMAX-RMIN)/MAX(ABS(RMIN),ABS(RMAX),tolerance)

        IGOT_T_CONST=PERCENT < tolerance2

        IGOT_T_CONST_VALUE = RMAX
        RETURN
    END SUBROUTINE IS_FIELD_CONSTANT





    SUBROUTINE PACK_LOC( LOC_F, T_ALL, NPHASE, IPT, IGOT_T_PACK )
        ! If PACK then pack T_ALL into LOC_F as long at IGOT_T==1 and STORE and not already in storage.
        IMPLICIT NONE
        INTEGER, intent( in ) :: NPHASE
        ! GLOBAL_FACE is the quadrature point which helps point into the storage memory
        INTEGER, intent( inout ) :: IPT
        LOGICAL, DIMENSION(:), intent( in ) :: IGOT_T_PACK
        REAL, DIMENSION(:), intent( in ) :: T_ALL
        REAL, DIMENSION(:), intent( inout ) :: LOC_F
        ! local variables...
        INTEGER :: IPHASE

        DO IPHASE=1,NPHASE
            IF(IGOT_T_PACK(IPHASE)) THEN ! Put into packing vector LOC_F
                LOC_F(IPT) = T_ALL(IPHASE)
                IPT=IPT+1
            ENDIF
        END DO

        RETURN
    END SUBROUTINE PACK_LOC




    SUBROUTINE I_PACK_LOC( LOC_F, T_ALL, NPHASE, IPT, IGOT_T_PACK )
        ! If PACK then pack T_ALL into LOC_F as long at IGOT_T==1 and STORE and not already in storage.
        IMPLICIT NONE
        INTEGER, intent( in ) :: NPHASE
        ! GLOBAL_FACE is the quadrature point which helps point into the storage memory
        INTEGER, intent( inout ) :: IPT
        LOGICAL, DIMENSION(NPHASE), intent( in ) :: IGOT_T_PACK
        INTEGER, DIMENSION(NPHASE), intent( in ) :: T_ALL
        INTEGER, DIMENSION(:), intent( inout ) :: LOC_F
        ! local variables...
        INTEGER :: IPHASE

        DO IPHASE=1,NPHASE
            IF(IGOT_T_PACK(IPHASE)) THEN ! Put into packing vector LOC_F
                LOC_F(IPT) = T_ALL(IPHASE)
                IPT=IPT+1
            ENDIF
        END DO

        RETURN
    END SUBROUTINE I_PACK_LOC


    SUBROUTINE UNPACK_LOC( LOC_F, T_ALL, NPHASE, IPT, IGOT_T_PACK, IGOT_T_CONST, IGOT_T_CONST_VALUE)
        ! If PACK then UNpack loc_f into T_ALL  as long at IGOT_T==1 and STORE and not already in storage.
        IMPLICIT NONE
        INTEGER, intent( in ) :: NPHASE
        !INTEGER, intent( in ) :: GLOBAL_FACE
        ! GLOBAL_FACE is the quadrature point which helps point into the storage memory
        INTEGER, intent( inout ) :: IPT
        LOGICAL, DIMENSION(NPHASE), intent( in ) :: IGOT_T_PACK, IGOT_T_CONST
        REAL, DIMENSION(NPHASE), intent( inout ) :: T_ALL
        REAL, DIMENSION(NPHASE), intent( in ) :: IGOT_T_CONST_VALUE
        REAL, DIMENSION(:), intent( in ) :: LOC_F
        ! local variables...
        INTEGER :: IPHASE


        DO IPHASE=1,NPHASE
            IF(IGOT_T_PACK(IPHASE)) THEN
                T_ALL(IPHASE) = LOC_F(IPT)
                IPT=IPT+1
                if (.not.IGOT_T_CONST(IPHASE)) then
                    !This option is not considered yet
                ENDIF
            ELSE IF(IGOT_T_CONST(IPHASE)) THEN
                T_ALL(IPHASE) = IGOT_T_CONST_VALUE(IPHASE)
            ELSE ! Set to 1 as last resort e.g. for T2, T2OLD
                T_ALL(IPHASE) = 1.0
            ENDIF
        END DO

        RETURN
    END SUBROUTINE UNPACK_LOC


    SUBROUTINE PACK_OR_UNPACK_LOC( LOC_F, T_ALL, NPHASE, NFIELD, IPT, PACK, STORE, IGOT_T )
        ! If PACK then pack T_ALL into LOC_F as long at IGOT_T==1 and STORE and not already in storage.
        LOGICAL, intent( in ) :: STORE, PACK
        INTEGER, intent( in ) :: NPHASE, IGOT_T, NFIELD
        !INTEGER, intent( in ) :: GLOBAL_FACE
        ! GLOBAL_FACE is the quadrature point which helps point into the storage memory
        INTEGER, intent( inout ) :: IPT
        REAL, DIMENSION(NPHASE), intent( inout ) :: T_ALL
        REAL, DIMENSION(NFIELD), intent( inout ) :: LOC_F
        ! local variables...
        LOGICAL :: IN_STORAGE

        IN_STORAGE = .FALSE.
        IF(STORE) THEN
            ! IF STORE then look to see if in storage
            IN_STORAGE = .FALSE.
        ENDIF

        IF(IGOT_T==1) THEN
            IF(PACK) THEN
                ! Pack solution into LOC_F
                IF(.NOT.IN_STORAGE) THEN
                    LOC_F(IPT:IPT-1+NPHASE) = T_ALL(:)
                    IPT=IPT+NPHASE
                ENDIF
            ELSE
                ! Unpack...
                IF(STORE) THEN
                    IF(.NOT.IN_STORAGE) THEN ! See if we are already storing limited value
                    ENDIF
                    ! Put LOC_F(1:NPHASE, CV_KLOC) = DEN_ALL( 1:NPHASE, CV_NODK ) into storage...
                    !                      T_ALL = LOC_F ???
                    IPT=IPT+NPHASE
                ELSE
                    ! Put LOC_F(1:NPHASE, CV_KLOC) = DEN_ALL( 1:NPHASE, CV_NODK ) into storage...
                    T_ALL(:) = LOC_F(IPT:IPT-1+NPHASE)
                    IPT=IPT+NPHASE
                ENDIF

            ENDIF ! ENDOF IF(PACK) THEN ELSE
        ENDIF ! END OF IF(IGOT_T==1) THEN
        RETURN
    END SUBROUTINE PACK_OR_UNPACK_LOC

    function CV_count_faces( Mdims, CV_ELE_TYPE, CV_GIdims) result(global_face)
      !  =====================================================================
      !     This subroutine counts then number of faces in the control volume space
      !
      ! Inputs/Outputs
      IMPLICIT NONE
      type(multi_dimensions), intent(in) :: Mdims
      INTEGER, intent( in ) :: CV_ELE_TYPE
      type( multi_GI_dimensions ), optional, intent(in) :: CV_GIdims

      ! Local variables
      type( multi_GI_dimensions ) :: GIdims
      INTEGER :: GLOBAL_FACE

        if (present(CV_GIdims)) then
            global_face=Mdims%totele * CV_GIdims%scvngi * 2
        else
            call retrieve_ngi( GIdims, Mdims, cv_ele_type, QUAD_OVER_WHOLE_ELE = .false. )
            global_face=Mdims%totele * GIdims%scvngi * 2
        end if

      return
    end function CV_COUNT_FACES





    SUBROUTINE FIND_OTHER_SIDE( CV_OTHER_LOC, CV_NLOC, U_OTHER_LOC, U_NLOC,  &
        MAT_OTHER_LOC, INTEGRAT_AT_GI, &
        X_NLOC, XU_NLOC, X_NDGLN, XU_NDGLN, &
        CV_SNLOC, CVFEM_ON_FACE, X_SHARE, ELE, ELE2,  &
        FINELE, COLELE, DISTCONTINUOUS_METHOD )
        ! We are on the boundary or next to another element. Determine CV_OTHER_LOC,
        ! U_OTHER_LOC.
        ! CVFEM_ON_FACE(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
        ! Look for these nodes on the other elements.
        ! ELE2=0 also when we are between elements but are trying to integrate across
        ! the middle of a CV.
        IMPLICIT NONE
        INTEGER, intent( in ) :: CV_NLOC, U_NLOC, X_NLOC, XU_NLOC, &
            &                   CV_SNLOC, ELE
        INTEGER, DIMENSION( : ), intent( in ) :: X_NDGLN, XU_NDGLN
        LOGICAL, DIMENSION( : ), intent( in ) :: CVFEM_ON_FACE
        INTEGER, DIMENSION( : ), intent( in ) :: FINELE
        INTEGER, DIMENSION( : ), intent( in ) :: COLELE

        INTEGER, DIMENSION( : ), intent( inout ) :: CV_OTHER_LOC, U_OTHER_LOC, MAT_OTHER_LOC
        LOGICAL, DIMENSION( : ), intent( inout ) :: X_SHARE
        INTEGER, intent( inout ) :: ELE2
        LOGICAL, intent( inout ) :: INTEGRAT_AT_GI
        LOGICAL, intent( in ) :: DISTCONTINUOUS_METHOD
        ! Local variables
        INTEGER :: X_KLOC, X_NODK, X_NODK2, COUNT, ELE3, SUF_COUNT, CV_KLOC, CV_KLOC2, &
            &     U_KLOC, U_KLOC2,  XU_NODK, XU_NODK2

        !ewrite(3,*) 'In FIND_OTHER_SIDE'

        DO X_KLOC = 1, X_NLOC
            X_NODK = X_NDGLN( ( ELE - 1) * X_NLOC + X_KLOC )
            X_SHARE( X_NODK ) = CVFEM_ON_FACE( X_KLOC )
        END DO

        ELE3 = 0
        DO COUNT = FINELE( ELE ), FINELE( ELE + 1 ) - 1, 1
            ELE2 = COLELE( COUNT )
            SUF_COUNT = 0 ! See if we share the same nodes
            IF ( ELE2 /= ELE ) THEN
                DO X_KLOC = 1, X_NLOC
                    X_NODK = X_NDGLN( ( ELE2 - 1 ) * X_NLOC + X_KLOC )
                    IF ( X_SHARE( X_NODK ) ) SUF_COUNT = SUF_COUNT + 1
                END DO
            END IF
            IF( SUF_COUNT == CV_SNLOC ) THEN
                ELE3 = ELE2
                EXIT
            ENDIF
           !ewrite(3,*)'suf_count:', ele, ele2, suf_count, cv_snloc
        END DO
        ELE2 = ELE3

        DO X_KLOC = 1, X_NLOC
            X_NODK = X_NDGLN( ( ELE - 1 ) * X_NLOC + X_KLOC )
            X_SHARE( X_NODK ) = .FALSE.
        END DO


        ! Quite because there is no work to do here...
        IF(.NOT.DISTCONTINUOUS_METHOD) THEN
            IF(ELE2.NE.0) THEN ! this is not on the boundary of the domain.
                INTEGRAT_AT_GI=.FALSE.
                RETURN
            ENDIF
        ENDIF


        IF ( ( ELE2 /= 0 ) .AND. INTEGRAT_AT_GI ) THEN ! Determine CV_OTHER_LOC(CV_KLOC)
            CV_OTHER_LOC = 0
            DO CV_KLOC = 1, CV_NLOC
                IF ( CVFEM_ON_FACE( CV_KLOC ) ) THEN ! Find opposite local node
                    X_NODK = X_NDGLN( ( ELE - 1 ) * X_NLOC + CV_KLOC )
                    DO CV_KLOC2 = 1, CV_NLOC
                        X_NODK2 = X_NDGLN( ( ELE2 - 1 ) * X_NLOC + CV_KLOC2 )
                        IF( X_NODK2 == X_NODK ) THEN
                            CV_OTHER_LOC( CV_KLOC ) = CV_KLOC2
                            EXIT
                        ENDIF
                    END DO
                END IF
            END DO

            U_OTHER_LOC = 0 ! Determine U_OTHER_LOC(U_KLOC)
            if (.not.is_P0DGP1CV) then!XU_NDGLN not defined for P0DGP1
            ! Works for non constant and constant (XU_NLOC=1) vel basis functions...
            DO U_KLOC = 1, U_NLOC ! Find opposite local node
                XU_NODK = XU_NDGLN( ( ELE - 1 ) * XU_NLOC + U_KLOC )
                DO U_KLOC2 = 1, U_NLOC
                    XU_NODK2 = XU_NDGLN(( ELE2 - 1 ) * XU_NLOC + U_KLOC2 )
                    IF ( ( XU_NODK2 == XU_NODK ) .OR. ( XU_NLOC == 1 ) ) &
                        U_OTHER_LOC( U_KLOC ) = U_KLOC2
                END DO
            END DO
            end if
            MAT_OTHER_LOC = CV_OTHER_LOC
        ELSE
            CV_OTHER_LOC = 0
            U_OTHER_LOC = 0
            MAT_OTHER_LOC = 0
        END IF

        RETURN

    END SUBROUTINE FIND_OTHER_SIDE

    SUBROUTINE PROJ_CV_TO_FEM(packed_state, &
        fempsi, psi, &
        Mdims, CV_GIdims, CV_funs, Mspars, ndgln, &
        igetct, X, mass_ele, mass_mn_pres, &
        tracer, psi_ave, psi_int)

        !------------------------------------------------
        ! Subroutine description:
        ! (1) determine FEMT (finite element wise) etc
        !     from T (control volume wise)
        ! (2) (optional) calculate psi_int (area) and
        !     psi_ave (barycentre) over each CV
        !------------------------------------------------

        implicit none

        !---------------------------------
        ! I/O variables
        !---------------------------------
        type(state_type) :: packed_state                                    ! local state data
        type(tensor_field_pointer), dimension(:), intent(inout) :: fempsi   ! finite element field data
        type(tensor_field_pointer), dimension(:), intent(inout) :: psi      ! finite volume field data
        type(multi_dimensions), intent(in) :: Mdims                         ! dimension data
        type(multi_GI_dimensions), intent(in) :: CV_GIdims                  ! gauss integer dimension data
        type(multi_shape_funs), intent(inout) :: CV_funs                    ! control volume shape function data
        type(multi_sparsities), intent(in) :: Mspars                        ! sparsity data
        type(multi_ndgln), intent(in) :: ndgln                              ! global numbering data
        integer, intent(in) :: igetct                                       ! whether to get CT matrix
        real, dimension(:,:), intent(in) :: X                               ! coordinates of the elements
        real, dimension(:), intent(inout) :: mass_ele                       ! finite element mass
        real, dimension(:), intent(inout) :: mass_mn_pres                   ! ??
        type(tensor_field), intent(in) :: tracer
        ! the following two need to be changed to optional in the future
        type(vector_field_pointer), dimension(:), intent(inout) :: psi_int ! control volume area
        type(vector_field_pointer), dimension(:), intent(inout) :: psi_ave ! control volume barycentre

        !---------------------------------
        ! local variables
        !---------------------------------
        type (multi_dev_shape_funs) :: DevFuns
        integer :: cv_nodi, cv_nodj, COUNT, max_iterations
        integer :: iele, cv_iloc, cv_jloc, cv_gi, it                         ! loop variables
        real :: nn, nm, mn, mm
        type(scalar_field) :: cv_mass
        type(tensor_field), dimension(size(psi)) :: fempsi_rhs
        type(vector_field), dimension(size(psi_ave)) :: psi_ave_temp
        type(vector_field), dimension(size(psi_int)) :: psi_int_temp
        type(tensor_field), pointer :: tfield
        type(csr_sparsity), pointer :: sparsity
        logical, parameter :: DCYL = .false.
        logical :: do_not_project, cv_test_space, is_to_update
        character(len=*), parameter :: projection_options = '/projections/control_volume_projections'
        character(len=option_path_len) :: option_path

        !---------------------------------
        ! initialisation and allocation
        !---------------------------------
        do_not_project = have_option(projection_options//'/do_not_project')
        cv_test_space = have_option(projection_options//'/test_function_space::ControlVolume')
        is_to_update = .not.associated(CV_funs%CV2FE%refcount)

        do it=1,size(fempsi)
            call zero(fempsi(it)%ptr)
            call allocate(fempsi_rhs(it),psi(it)%ptr%mesh,"RHS",dim=psi(it)%ptr%dim)
            call zero(fempsi_rhs(it))
            call halo_update(psi(it)%ptr)
        end do
        tfield=>psi(1)%ptr

        if(is_to_update) then
            call allocate(cv_mass,psi(1)%ptr%mesh,'CV_mass')
            call zero(cv_mass)

            if(tracer%mesh%continuity<0) then
                psi_ave(1)%ptr%val = X(:,ndgln%x)
            else
                call set_all(psi_ave(1)%ptr,X)
            end if
            do it = 1,size(psi_ave)
                call allocate(psi_ave_temp(it),psi_ave(it)%ptr%dim,psi_ave(it)%ptr%mesh,"PsiAveTemp")
                call set(psi_ave_temp(it),psi_ave(it)%ptr)
                call zero(psi_ave(it)%ptr)
            end do

            call set(psi_int(1)%ptr,dim=1,val=1.0)
            do it = 1,size(psi_int)
                call allocate(psi_int_temp(it),psi_int(it)%ptr%dim,psi_int(it)%ptr%mesh,"PsiIntTemp")
                call set(psi_int_temp(it),psi_int(it)%ptr)
                call zero(psi_int(it)%ptr)
            end do

            sparsity=>extract_csr_sparsity(packed_state,"PressureMassMatrixSparsity")
            call allocate(CV_funs%CV2FE,sparsity,[1,1],name="ProjectionMatrix")
            call zero(CV_funs%CV2FE)
        end if
        if(igetct/=0) mass_mn_pres=0.0

        !---------------------------------
        ! projection
        !---------------------------------
        call allocate_multi_dev_shape_funs(CV_funs, DevFuns)

        Loop_Elements: do iele = 1, Mdims%totele
            ! check parallelisation
            if(isParallel()) then
                if(.not.assemble_ele(psi_int(1)%ptr,iele)) cycle ! IS THIS A PROBLEM?
            end if

            ! calculate detwei,RA,NX,NY,NZ for the ith element
            call DETNLXR(iele, X, ndgln%x, CV_funs%cvweight, CV_funs%CVFEN, CV_funs%CVFENLX_ALL, DevFuns)

            mass_ele(iele) = DevFuns%volume

            Loop_CV_iLoc: do cv_iloc = 1, Mdims%cv_nloc
                cv_nodi = ndgln%cv((iele-1)*Mdims%cv_nloc+cv_iloc)
                if(.not.node_owned(psi(1)%ptr,cv_nodi)) cycle

                Loop_CV_jLoc: do cv_jloc = 1, Mdims%cv_nloc
                    cv_nodj = ndgln%cv((iele-1)*Mdims%cv_nloc+cv_jloc)

                    nn = 0.0; nm = 0.0; mn = 0.0; mm = 0.0
                    do cv_gi = 1, CV_GIdims%cv_ngi
                        mn = mn+CV_funs%CVN(cv_iloc,cv_gi)*CV_funs%CVFEN(cv_jloc,cv_gi)*DevFuns%detwei(cv_gi)
                        mm = mm+CV_funs%CVN(cv_iloc,cv_gi)*CV_funs%CVN(cv_jloc,cv_gi)*DevFuns%detwei(cv_gi)
                        if(.not.cv_test_space) then
                            nn = nn+CV_funs%CVFEN(cv_iloc,cv_gi)*CV_funs%CVFEN(cv_jloc,cv_gi)*DevFuns%detwei(cv_gi)
                            nm = nm+CV_funs%CVFEN(cv_iloc,cv_gi)*CV_funs%CVN(cv_jloc,cv_gi)*DevFuns%detwei(cv_gi)
                        end if
                    end do

                    if(cv_test_space) then
                        if(is_to_update) call addto(CV_funs%CV2FE,1,1,cv_nodi,cv_nodj,mn)
                        do it = 1, size(fempsi_rhs)
                            fempsi_rhs(it)%val(:,:,cv_nodi) = fempsi_rhs(it)%val(:,:,cv_nodi) &
                                +mm*psi(it)%ptr%val(:,:,cv_nodj)
                        end do
                    else
                        if(is_to_update) call addto(CV_funs%CV2FE,1,1,cv_nodi,cv_nodj,nn)
                        do it = 1, size(fempsi_rhs)
                            fempsi_rhs(it)%val(:,:,cv_nodi) = fempsi_rhs(it)%val(:,:,cv_nodi) &
                                +nm*psi(it)%ptr%val(:,:,cv_nodj)
                        end do
                    end if

                    if(igetct/=0) then
                        call PosInMat(COUNT,cv_nodi,cv_nodj,Mspars%CMC%fin,Mspars%CMC%col)
                        mass_mn_pres(COUNT) = mass_mn_pres(COUNT)+mn
                    end if

                    ! mass, barycentre and area of CVs
                    if(is_to_update) then
                        call addto(cv_mass,cv_nodi,mm)
                        do it = 1, size(psi_ave)
                            call addto(psi_ave(it)%ptr,node_number=cv_nodi,val=mn*psi_ave_temp(it)%val(:,cv_nodj))
                        end do
                        do it = 1, size(psi_int)
                            call addto(psi_int(it)%ptr,node_number=cv_nodi,val=mn*psi_int_temp(it)%val(:,cv_nodj))
                        end do
                    end if

                end do Loop_CV_jLoc
            end do Loop_CV_iLoc
        end do Loop_Elements

        !---------------------------------
        ! solving, updating & deallocating
        !---------------------------------
        ! update the halo information
        if(is_to_update) then
            call halo_update(cv_mass)
            call invert(cv_mass)
            do it = 1, size(psi_ave)
                call halo_update(psi_ave(it)%ptr)
                call scale(psi_ave(it)%ptr,cv_mass)
            end do
            do it = 1, size(psi_int) ! this seems a really bad way to do halo updates.
                call halo_update(psi_int(it)%ptr)
            end do
        end if

        ! solve the petsc matrix
        if(do_not_project) then
            do it = 1, size(fempsi)
                call set(fempsi(it)%ptr,psi(it)%ptr)
            end do
        else
            if(have_option(projection_options//'/solver')) then
                option_path=projection_options//'/solver'
            else
                call get_option(trim(psi(1)%ptr%option_path)//"/prognostic/solver/max_iterations", &
                    max_iterations,default=500)
                if(max_iterations==0) then
                    option_path="/material_phase[0]/scalar_field::Pressure/prognostic"
                else
                    option_path=trim(psi(1)%ptr%option_path)//"/prognostic"
                end if
            end if
            do it = 1, size(fempsi)
                call zero_non_owned(fempsi_rhs(it))
                call petsc_solve(fempsi(it)%ptr,CV_funs%CV2FE,fempsi_rhs(it),option_path = option_path)
            end do
        end if

        ! deallocation
        call deallocate_multi_dev_shape_funs(DevFuns)

        do it = 1,size(fempsi_rhs)
            call deallocate(fempsi_rhs(it))
        end do
        if(is_to_update) then
            call deallocate(cv_mass)
            do it = 1,size(psi_ave_temp)
                call deallocate(psi_ave_temp(it))
            end do
            do it = 1,size(psi_int_temp)
                call deallocate(psi_int_temp(it))
            end do
        end if

        return

    END SUBROUTINE PROJ_CV_TO_FEM



    SUBROUTINE DG_DERIVS_ALL1( FEMT, FEMTOLD, &
        DTX_ELE, DTOLDX_ELE, &
        NDIM, NPHASE, NCOMP, TOTELE, CV_NDGLN, & ! ncomp = ndim here
        XCV_NDGLN, X_NLOC, X_NDGLN,&
        CV_NGI, CV_NLOC, CVWEIGHT, &
        N, NLX, NLY, NLZ, &
        X_N, X_NLX, X_NLY, X_NLZ, &
        X_NONODS, X, Y, Z, &
        NFACE, FACE_ELE, CV_SLOCLIST, X_SLOCLIST, CV_SNLOC, X_SNLOC, WIC_T_BC, SUF_T_BC, &
        SBCVNGI, SBCVFEN, SBWEIGH, &
        X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY, get_gradU, state )

        ! calculates derivatives of vector fields
        IMPLICIT NONE

        INTEGER, intent( in ) :: NDIM, NPHASE, NCOMP, TOTELE, X_NLOC, CV_NGI, CV_NLOC, &
            X_NONODS, CV_SNLOC, X_SNLOC, SBCVNGI, NFACE
        REAL, DIMENSION( :, :, : ), intent( in ) :: FEMT, FEMTOLD
        REAL, DIMENSION( :, :, :, :, : ), intent( inout ) :: DTX_ELE, DTOLDX_ELE
        INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) ::  XCV_NDGLN
        INTEGER, DIMENSION( :, :, : ), intent( in ) ::  WIC_T_BC
        REAL, DIMENSION( :, :, : ), intent( in ) ::  SUF_T_BC
        INTEGER, DIMENSION( :, : ), intent( in ) ::  CV_SLOCLIST
        INTEGER, DIMENSION( :, : ), intent( in ) ::  X_SLOCLIST
        INTEGER, DIMENSION( :, : ), intent( in ) ::  FACE_ELE
        REAL, DIMENSION( : ), intent( in ) :: CVWEIGHT
        REAL, DIMENSION( :, : ), intent( in ) :: N, NLX, NLY, NLZ
        REAL, DIMENSION( :, : ), intent( in ) :: X_N, X_NLX, X_NLY, X_NLZ
        REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
        REAL, DIMENSION( :, : ), intent( in ) :: SBCVFEN
        REAL, DIMENSION( :, : ), intent( in ) :: X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY
        REAL, DIMENSION( : ), intent( in ) :: SBWEIGH
        LOGICAL, intent( in ) :: get_gradU
        TYPE( STATE_TYPE), DIMENSION( : ), intent( inout ) :: state
        ! Local variables
        REAL, DIMENSION( :, :, : ), ALLOCATABLE :: MASELE
        REAL, DIMENSION( :, :, :, :, : ), ALLOCATABLE :: VTX_ELE, VTOLDX_ELE
        LOGICAL :: D1, D3, APPLYBC( NCOMP, NPHASE )
        LOGICAL, PARAMETER :: DCYL = .FALSE.
        REAL, dimension( CV_NGI ) :: DETWEI, RA
        REAL, DIMENSION( NDIM, size(NLX,1), CV_NGI):: NX_ALL
        REAL, DIMENSION( NDIM, size(X_NLX,1),CV_NGI ) :: X_NX_ALL
        REAL :: VOLUME
        REAL, DIMENSION( CV_NLOC, CV_NLOC )  :: MASS, INV_MASS
        REAL, DIMENSION( NDIM, X_SNLOC ) :: XSL( 3, X_SNLOC ), SNORMXN( 3, SBCVNGI ), SDETWE( SBCVNGI )
        INTEGER  :: SLOC2LOC( CV_SNLOC ), X_SLOC2LOC( X_SNLOC ), ILOC_OTHER_SIDE( CV_SNLOC )
        REAL :: NN, NNX( NDIM ), NORMX( 3 ), SAREA, NRBC, RTBC, VLM_NORX( NDIM )
        INTEGER :: ELE, CV_ILOC, CV_JLOC, CV_NODI, CV_NODJ, CV_ILOC2, &
            CV_INOD, CV_INOD2, CV_JLOC2, CV_NODJ2, &
            CV_SILOC, CV_SJLOC, CV_SJLOC2, ELE2, IFACE, IPHASE, SELE2, SUF_CV_SJ2, &
            X_INOD, X_SILOC, X_ILOC, ICOMP, IDIM, STAT
        INTEGER, PARAMETER :: WIC_T_BC_DIRICHLET = 1
        TYPE( TENSOR_FIELD ), POINTER :: GRADU

        ewrite(3,*)'in DG_DERIVS'

        DTX_ELE = 0.0 ; DTOLDX_ELE = 0.0

        ALLOCATE( MASELE( CV_NLOC, CV_NLOC, TOTELE ) )
        ALLOCATE( VTX_ELE( NDIM, NCOMP, NPHASE, CV_NLOC, TOTELE ) )
        ALLOCATE( VTOLDX_ELE( NDIM, NCOMP, NPHASE, CV_NLOC, TOTELE ) )

        MASELE = 0.0
        VTX_ELE = 0.0

        VTOLDX_ELE = 0.0

        D1 = ( NDIM == 1 )
        D3 = ( NDIM == 3 )
        !DCYL = .FALSE.

        Loop_Elements1: DO ELE = 1, TOTELE

            ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
            CALL DETNLXR_PLUS_U( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
                X_NLOC, X_NLOC, CV_NGI, &
                X_N, X_NLX, X_NLY, X_NLZ, CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
                X_NX_ALL, &
                CV_NLOC, NLX, NLY, NLZ, NX_ALL)
            Loop_CV_ILOC: DO CV_ILOC = 1, CV_NLOC

                CV_NODI = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )

                Loop_CV_JLOC: DO CV_JLOC = 1, CV_NLOC

                    CV_NODJ = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_JLOC )

                    NN  = SUM( N( CV_ILOC, : ) * N(  CV_JLOC, : ) * DETWEI )
                    NNX = MATMUL( NX_ALL( :, CV_JLOC, : ), N( CV_ILOC, : )  * DETWEI )

                    MASELE( CV_ILOC, CV_JLOC, ELE) = MASELE( CV_ILOC, CV_JLOC, ELE ) + NN

                    DO IPHASE = 1, NPHASE
                        DO ICOMP = 1, NCOMP
                            VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) = &
                                VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                                + NNX (:) * FEMT( ICOMP, IPHASE, CV_NODJ )

                            VTOLDX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) = &
                                VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                                + NNX (:) * FEMTOLD( ICOMP, IPHASE, CV_NODJ )
                        END DO
                    END DO

                END DO Loop_CV_JLOC

            END DO Loop_CV_ILOC

        END DO Loop_Elements1


        Loop_Elements2: DO ELE = 1, TOTELE

            Between_Elements_And_Boundary: DO IFACE = 1, NFACE
                ELE2 = FACE_ELE( IFACE, ELE )
                SELE2 = MAX( 0, - ELE2 )
                ELE2 = MAX( 0, + ELE2 )
                !ewrite(3,*)'FACE_ELE( 1, ELE ),FACE_ELE( 2, ELE ):',FACE_ELE( 1, ELE ),FACE_ELE( 2, ELE )

                ! The surface nodes on element face IFACE.
                SLOC2LOC( : ) = CV_SLOCLIST( IFACE, : )
                X_SLOC2LOC( : ) = X_SLOCLIST( IFACE, : )

                ! Form approximate surface normal (NORMX,NORMY,NORMZ)
                CALL DGSIMPLNORM( ELE, X_SLOC2LOC, X_NLOC, X_SNLOC, X_NDGLN, &
                    X, Y, Z, NORMX( 1 ), NORMX( 2 ), NORMX( 3 ) )

                ! Recalculate the normal...
                DO X_SILOC = 1, X_SNLOC
                    X_ILOC = X_SLOC2LOC( X_SILOC )
                    X_INOD = X_NDGLN(( ELE - 1 ) * X_NLOC + X_ILOC )
                    if (NDIM >= 1) XSL( 1, X_SILOC ) = X( X_INOD )
                    if (NDIM >= 2) XSL( 2, X_SILOC ) = Y( X_INOD )
                    if (NDIM >= 3) XSL( 3, X_SILOC ) = Z( X_INOD )
                END DO

                CALL DGSDETNXLOC2(X_SNLOC, SBCVNGI, &
                    XSL( 1, : ), XSL( 2, : ), XSL( 3, : ), &
                    X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY, SBWEIGH, SDETWE, SAREA, &
                    (NDIM==1), (NDIM==3), (NDIM==-2), &
                    SNORMXN( 1, : ), SNORMXN( 2, : ), SNORMXN( 3, : ), &
                    NORMX( 1 ), NORMX( 2 ), NORMX( 3 ) )

                IF ( SELE2 == 0 ) THEN
                    ! Calculate the nodes on the other side of the face:

                    DO CV_SILOC = 1, CV_SNLOC
                        CV_ILOC = SLOC2LOC( CV_SILOC )
                        CV_INOD = XCV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )

                        DO CV_ILOC2 = 1, CV_NLOC
                            CV_INOD2 = XCV_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_ILOC2 )

                            IF( CV_INOD2 == CV_INOD ) ILOC_OTHER_SIDE( CV_SILOC ) = CV_ILOC2
                        END DO
                    END DO

                    APPLYBC = (ELE /= ELE2) .AND. (ELE2 /= 0)

                ELSE
                    APPLYBC = ( WIC_T_BC( :, :, SELE2 ) == WIC_T_BC_DIRICHLET )
                END IF

                DO CV_SILOC = 1, CV_SNLOC
                    CV_ILOC = SLOC2LOC( CV_SILOC )
                    DO CV_SJLOC = 1, CV_SNLOC
                        CV_JLOC = SLOC2LOC( CV_SJLOC )
                        CV_NODJ = CV_NDGLN( (ELE-1)*CV_NLOC + CV_JLOC )
                        IF ( SELE2 /= 0 ) THEN
                            CV_JLOC2 = CV_JLOC
                            CV_SJLOC2 = CV_SJLOC
                            CV_NODJ2 = CV_NODJ
                            SUF_CV_SJ2 = CV_SJLOC + CV_SNLOC * ( SELE2 - 1 )
                            NRBC = 0.0
                        ELSE
                            CV_JLOC2 = ILOC_OTHER_SIDE( CV_SJLOC )
                            CV_NODJ2 = CV_NDGLN( (ELE2-1)*CV_NLOC + CV_JLOC2 )
                            NRBC = 1.0
                        END IF

                        ! Have a surface integral on element boundary...
                        VLM_NORX(:) = MATMUL( SNORMXN( 1:NDIM, : ), &
                            SDETWE(:) * SBCVFEN( CV_SILOC, : ) * SBCVFEN( CV_SJLOC, : ) )

                        ! add diffusion term...
                        DO IPHASE = 1, NPHASE
                            DO ICOMP = 1, NCOMP
                                IF ( APPLYBC( ICOMP, IPHASE ) ) THEN

                                    IF ( SELE2 /= 0 ) THEN
                                        IF ( WIC_T_BC( ICOMP, IPHASE, SELE2 ) == WIC_T_BC_DIRICHLET ) THEN

                                            RTBC = SUF_T_BC( ICOMP, IPHASE, CV_SJLOC2 + CV_SNLOC * ( SELE2 - 1 ) )

                                            VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) = VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                                                - VLM_NORX(:) * 0.5 * ( FEMT( ICOMP, IPHASE, CV_NODJ ) - RTBC  )

                                            VTOLDX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) = VTOLDX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                                                - VLM_NORX(:) * 0.5 * ( FEMTOLD( ICOMP, IPHASE, CV_NODJ ) - RTBC  )

                                        END IF
                                    ELSE
                                        VTX_ELE(:, ICOMP, IPHASE, CV_ILOC, ELE ) = &
                                            VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                                            - VLM_NORX(:) * 0.5 * ( FEMT( ICOMP, IPHASE, CV_NODJ ) - FEMT( ICOMP, IPHASE, CV_NODJ2 )  )
                                        VTOLDX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) = &
                                            VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                                            - VLM_NORX(:) * 0.5 * ( FEMTOLD( ICOMP, IPHASE, CV_NODJ ) - FEMTOLD( ICOMP, IPHASE, CV_NODJ2 )  )

                                    END IF

                                END IF
                            END DO
                        END DO
                    END DO

                END DO

            END DO Between_Elements_And_Boundary


        END DO Loop_Elements2


        Loop_Elements3: DO ELE = 1, TOTELE

            MASS( :, : ) = MASELE( :, :, ELE )
            CALL MATDMATINV( MASS, INV_MASS, CV_NLOC )

            FORALL ( IDIM = 1:NDIM, ICOMP = 1:NCOMP, IPHASE = 1:NPHASE )

                DTX_ELE( ICOMP, IDIM, IPHASE, :, ELE ) = MATMUL( INV_MASS( :, : ), VTX_ELE( IDIM, ICOMP, IPHASE, :, ELE ) )
                DTOLDX_ELE( ICOMP, IDIM, IPHASE, :, ELE ) = MATMUL( INV_MASS( :, : ) , VTOLDX_ELE( IDIM, ICOMP, IPHASE, :, ELE ) )

            END FORALL

        END DO Loop_Elements3

        ! set gradU
        if ( get_gradU ) then
           do iphase = 1, nphase
              gradU => extract_tensor_field( state( iphase ), "gradU", stat )
              if ( stat == 0 ) then
                 do ele = 1, totele
                    do cv_iloc = 1, cv_nloc
                       cv_nodi = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )
                       gradU%val( :, :, cv_nodi ) = dtx_ele( :, :, iphase, cv_iloc, ele )
                    end do
                 end do
              end if
           end do
        end if


        DEALLOCATE( MASELE, VTX_ELE, VTOLDX_ELE )

        ewrite(3,*)'about to leave DG_DERIVS'

        RETURN

    END SUBROUTINE DG_DERIVS_ALL1


    SUBROUTINE DG_DERIVS_ALL2( FEMT, FEMTOLD, &
        DTX_ELE, DTOLDX_ELE, &
        NDIM, NPHASE, TOTELE, CV_NDGLN, &
        XCV_NDGLN, X_NLOC, X_NDGLN,&
        CV_NGI, CV_NLOC, CVWEIGHT, &
        N, NLX, NLY, NLZ, &
        X_N, X_NLX, X_NLY, X_NLZ, &
        X_NONODS, X, Y, Z, &
        NFACE, FACE_ELE, CV_SLOCLIST, X_SLOCLIST, CV_SNLOC, X_SNLOC, WIC_T_BC, SUF_T_BC, &
        SBCVNGI, SBCVFEN, SBWEIGH, &
        X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY)

        ! determine FEMT (finite element wise) etc from T (control volume wise)
        IMPLICIT NONE

        INTEGER, intent( in ) :: NDIM, NPHASE, TOTELE, X_NLOC, CV_NGI, CV_NLOC, &
            &                   X_NONODS, CV_SNLOC, X_SNLOC, SBCVNGI, NFACE
        REAL, DIMENSION( :, : ), intent( in ) :: FEMT, FEMTOLD
        REAL, DIMENSION( :, :, :, : ), intent( inout ) :: DTX_ELE, DTOLDX_ELE
        INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) ::  XCV_NDGLN
        INTEGER, DIMENSION( :,  :, : ), intent( in ) ::  WIC_T_BC
        REAL, DIMENSION( :, :, : ), intent( in ) ::  SUF_T_BC
        INTEGER, DIMENSION( :, : ), intent( in ) ::  CV_SLOCLIST
        INTEGER, DIMENSION( :, : ), intent( in ) ::  X_SLOCLIST
        INTEGER, DIMENSION( :, : ), intent( in ) ::  FACE_ELE
        REAL, DIMENSION( : ), intent( in ) :: CVWEIGHT
        REAL, DIMENSION( :, : ), intent( in ) :: N, NLX, NLY, NLZ
        REAL, DIMENSION( :, : ), intent( in ) :: X_N, X_NLX, X_NLY, X_NLZ
        REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
        REAL, DIMENSION( :, : ), intent( in ) :: SBCVFEN
        REAL, DIMENSION( :, : ), intent( in ) :: X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY
        REAL, DIMENSION( : ), intent( in ) :: SBWEIGH

        ! Local variables
        REAL, DIMENSION( :, :, : ), ALLOCATABLE :: MASELE
        REAL, DIMENSION( :, :, :, : ), ALLOCATABLE :: VTX_ELE, VTOLDX_ELE
        LOGICAL :: D1, D3, APPLYBC( NPHASE )
        LOGICAL, PARAMETER :: DCYL = .FALSE.
        REAL, dimension( CV_NGI ) :: DETWEI, RA
        REAL, DIMENSION( NDIM, size(NLX,1), CV_NGI):: NX_ALL
        REAL, DIMENSION( NDIM, size(X_NLX,1),CV_NGI ) :: X_NX_ALL
        REAL :: VOLUME
        REAL, DIMENSION( CV_NLOC, CV_NLOC )  :: MASS, INV_MASS
        REAL, DIMENSION( NDIM, X_SNLOC ) :: XSL( 3, X_SNLOC ), SNORMXN( NDIM, SBCVNGI ), SDETWE( SBCVNGI )
        INTEGER  :: SLOC2LOC( CV_SNLOC ), X_SLOC2LOC( X_SNLOC ), ILOC_OTHER_SIDE( CV_SNLOC )
        REAL :: NN, NNX( NDIM ), NORMX( 3 ), SAREA, NRBC, RTBC, VLM_NORX( NDIM )
        INTEGER :: ELE, CV_ILOC, CV_JLOC, CV_NODI, CV_NODJ, CV_ILOC2, &
            CV_INOD, CV_INOD2, CV_JLOC2, CV_NODJ2, &
            CV_SILOC, CV_SJLOC, CV_SJLOC2, ELE2, IFACE, IPHASE, SELE2, SUF_CV_SJ2, &
            X_INOD, X_SILOC, X_ILOC, IDIM

        ewrite(3,*)'in DG_DERIVS'

        DTX_ELE = 0.0 ; DTOLDX_ELE = 0.0

        ALLOCATE( MASELE( CV_NLOC, CV_NLOC, TOTELE ) )
        ALLOCATE( VTX_ELE( NDIM, NPHASE, CV_NLOC, TOTELE ) )
        ALLOCATE( VTOLDX_ELE( NDIM, NPHASE, CV_NLOC, TOTELE ) )

        MASELE = 0.0 ; VTX_ELE = 0.0 ; VTOLDX_ELE = 0.0

        D1 = ( NDIM == 1 )
        D3 = ( NDIM == 3 )
        !DCYL = .FALSE.
        Loop_Elements1: DO ELE = 1, TOTELE

            ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
            CALL DETNLXR_PLUS_U( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
                X_NLOC, X_NLOC, CV_NGI, &
                X_N, X_NLX, X_NLY, X_NLZ, CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
                X_NX_ALL, &
                CV_NLOC, NLX, NLY, NLZ, NX_ALL)


            Loop_CV_ILOC: DO CV_ILOC = 1, CV_NLOC

                CV_NODI = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )

                Loop_CV_JLOC: DO CV_JLOC = 1, CV_NLOC

                    CV_NODJ = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_JLOC )

                    NN  = SUM( N( CV_ILOC, : ) * N(  CV_JLOC, : ) * DETWEI )
                    NNX = MATMUL( NX_ALL( :, CV_JLOC, : ), N( CV_ILOC, : )  * DETWEI )

                    MASELE( CV_ILOC, CV_JLOC, ELE) = MASELE( CV_ILOC, CV_JLOC, ELE ) + NN

                    DO IPHASE = 1, NPHASE
                        VTX_ELE( :, IPHASE, CV_ILOC, ELE ) = &
                            VTX_ELE( :, IPHASE, CV_ILOC, ELE ) &
                            + NNX (:) * FEMT( IPHASE, CV_NODJ )

                        VTOLDX_ELE( :, IPHASE, CV_ILOC, ELE ) = &
                            VTOLDX_ELE( :, IPHASE, CV_ILOC, ELE ) &
                            + NNX (:) * FEMTOLD( IPHASE, CV_NODJ )
                    END DO

                END DO Loop_CV_JLOC

            END DO Loop_CV_ILOC

        END DO Loop_Elements1


        Loop_Elements2: DO ELE = 1, TOTELE

            Between_Elements_And_Boundary: DO IFACE = 1, NFACE
                ELE2 = FACE_ELE( IFACE, ELE )
                SELE2 = MAX( 0, - ELE2 )
                ELE2 = MAX( 0, + ELE2 )
                !ewrite(3,*)'FACE_ELE( 1, ELE ),FACE_ELE( 2, ELE ):',FACE_ELE( 1, ELE ),FACE_ELE( 2, ELE )

                ! The surface nodes on element face IFACE.
                SLOC2LOC( : ) = CV_SLOCLIST( IFACE, : )
                X_SLOC2LOC( : ) = X_SLOCLIST( IFACE, : )

                ! Form approximate surface normal (NORMX,NORMY,NORMZ)
                CALL DGSIMPLNORM( ELE, X_SLOC2LOC, X_NLOC, X_SNLOC, X_NDGLN, &
                    X, Y, Z, NORMX( 1 ), NORMX( 2 ), NORMX( 3 ) )

                ! Recalculate the normal...
                DO X_SILOC = 1, X_SNLOC
                    X_ILOC = X_SLOC2LOC( X_SILOC )
                    X_INOD = X_NDGLN(( ELE - 1 ) * X_NLOC + X_ILOC )
                    if (NDIM >= 1) XSL( 1, X_SILOC ) = X( X_INOD )
                    if (NDIM >= 2) XSL( 2, X_SILOC ) = Y( X_INOD )
                    if (NDIM >= 3) XSL( 3, X_SILOC ) = Z( X_INOD )
                END DO

                IF ( SELE2 == 0 ) THEN
                    ! Calculate the nodes on the other side of the face:

                    DO CV_SILOC = 1, CV_SNLOC
                        CV_ILOC = SLOC2LOC( CV_SILOC )
                        CV_INOD = XCV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )

                        DO CV_ILOC2 = 1, CV_NLOC
                            CV_INOD2 = XCV_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_ILOC2 )

                            IF( CV_INOD2 == CV_INOD ) ILOC_OTHER_SIDE( CV_SILOC ) = CV_ILOC2
                        END DO
                    END DO

                    APPLYBC = (ELE /= ELE2) .AND. (ELE2 /= 0)

                ELSE
                    APPLYBC = ( WIC_T_BC( 1, :, SELE2 ) == WIC_T_BC_DIRICHLET )
                END IF



                DO CV_SILOC = 1, CV_SNLOC
                    CV_ILOC = SLOC2LOC( CV_SILOC )
                    DO CV_SJLOC = 1, CV_SNLOC
                        CV_JLOC = SLOC2LOC( CV_SJLOC )
                        CV_NODJ = CV_NDGLN( (ELE-1)*CV_NLOC + CV_JLOC )
                        IF ( SELE2 /= 0 ) THEN
                            CV_JLOC2 = CV_JLOC
                            CV_SJLOC2 = CV_SJLOC
                            CV_NODJ2 = CV_NODJ
                            SUF_CV_SJ2 = CV_SJLOC + CV_SNLOC * ( SELE2 - 1 )
                            NRBC = 0.0
                        ELSE
                            CV_JLOC2 = ILOC_OTHER_SIDE( CV_SJLOC )
                            CV_NODJ2 = CV_NDGLN( (ELE2-1)*CV_NLOC + CV_JLOC2 )
                            NRBC = 1.0
                        END IF

                        ! Have a surface integral on element boundary...
                        VLM_NORX(:) = MATMUL( SNORMXN( :, : ), &
                            SDETWE(:) * SBCVFEN( CV_SILOC, : ) * SBCVFEN( CV_SJLOC, : ) )

                        ! add diffusion term...
                        DO IPHASE = 1, NPHASE
                            IF ( APPLYBC( IPHASE ) ) THEN

                                VTX_ELE( :, IPHASE, CV_ILOC, ELE ) = &
                                    VTX_ELE( :, IPHASE, CV_ILOC, ELE ) &
                                    - VLM_NORX(:) * 0.5 * ( FEMT( IPHASE, CV_NODJ ) - FEMT( IPHASE, CV_NODJ2 ) * NRBC )
                                VTOLDX_ELE( :, IPHASE, CV_ILOC, ELE ) = &
                                    VTOLDX_ELE( :, IPHASE, CV_ILOC, ELE ) &
                                    - VLM_NORX(:) * 0.5 * ( FEMTOLD( IPHASE, CV_NODJ ) - FEMTOLD( IPHASE, CV_NODJ2 ) * NRBC )

                                IF ( SELE2 /= 0 ) THEN
                                    IF ( WIC_T_BC(1,  IPHASE, SELE2 ) == WIC_T_BC_DIRICHLET ) THEN

                                        RTBC = SUF_T_BC( 1,IPHASE, CV_SJLOC2+ CV_SNLOC*( SELE2-1) )

                                        VTX_ELE( :, IPHASE, CV_ILOC, ELE ) = VTX_ELE( :, IPHASE, CV_ILOC, ELE ) &
                                            + VLM_NORX(:) * 0.5 * RTBC
                                        VTOLDX_ELE( :, IPHASE, CV_ILOC, ELE ) = VTOLDX_ELE( :, IPHASE, CV_ILOC, ELE ) &
                                            + VLM_NORX(:) * 0.5 * RTBC

                                    END IF
                                END IF
                            END IF
                        END DO
                    END DO
                END DO

            END DO Between_Elements_And_Boundary

        END DO Loop_Elements2


        Loop_Elements3: DO ELE = 1, TOTELE

            MASS( :, : ) = MASELE( :, :, ELE )
            INV_MASS=MASS
            !       CALL MATDMATINV( MASS, INV_MASS, CV_NLOC )
            CALL INVERT( INV_MASS )

            FORALL ( IDIM = 1:NDIM, IPHASE = 1:NPHASE )

                DTX_ELE( IDIM, IPHASE, :, ELE ) = MATMUL( INV_MASS( :, : ), VTX_ELE( IDIM, IPHASE, :, ELE ) )
                DTOLDX_ELE( IDIM, IPHASE, :, ELE ) = MATMUL( INV_MASS( :, : ) , VTOLDX_ELE( IDIM, IPHASE, :, ELE ) )

               !DTX_ELE( IDIM, :, IPHASE, ELE ) = MATMUL( INV_MASS( :, : ), VTX_ELE( IDIM, IPHASE, :, ELE ) )
               !DTOLDX_ELE( IDIM, :, IPHASE, ELE ) = MATMUL( INV_MASS( :, : ) , VTOLDX_ELE( IDIM, IPHASE, :, ELE ) )
            END FORALL

        END DO Loop_Elements3

        DEALLOCATE( MASELE, VTX_ELE, VTOLDX_ELE )

        ewrite(3,*)'about to leave DG_DERIVS'

    END SUBROUTINE DG_DERIVS_ALL2


    SUBROUTINE DIFFUS_CAL_COEFF(DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX,  &
        CV_NLOC, MAT_NLOC, NPHASE, MAT_NDGLN, &
        SMATFEN, SCVFEN, GI, NDIM, TDIFFUSION, DIFF_GI_ADDED, &
        HDC, &
        T_CV_NODJ, T_CV_NODI, &
        TOLD_CV_NODJ, TOLD_CV_NODI, &
        ELE, ELE2, CVNORMX_ALL,  &
        LOC_DTX_ELE_ALL, LOC_DTOLDX_ELE_ALL, LOC2_DTX_ELE_ALL, LOC2_DTOLDX_ELE_ALL, &
        LOC_WIC_T_BC, CV_OTHER_LOC, MAT_OTHER_LOC, CV_SNLOC, CV_SLOC2LOC, &
        on_domain_boundary, between_elements )
        ! This sub calculates the effective diffusion coefficientd DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX
        ! based on a non-linear method and a non-oscillating scheme.
        IMPLICIT NONE
        INTEGER, intent( in ) :: CV_NLOC, MAT_NLOC, NPHASE, &
            &                    GI, NDIM, ELE, ELE2, &
            &                    CV_SNLOC
        REAL, intent( in ) :: HDC
        LOGICAL, intent( in ) :: on_domain_boundary, between_elements
        REAL, DIMENSION( NPHASE ), intent( in ) :: T_CV_NODJ, T_CV_NODI, &
            &                                     TOLD_CV_NODJ, TOLD_CV_NODI
        REAL, DIMENSION( NPHASE ), intent( inout ) :: DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX
        INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: CV_SLOC2LOC
        INTEGER, DIMENSION( NPHASE ), intent( in ) :: LOC_WIC_T_BC
        INTEGER, DIMENSION( : ), intent( in ) :: CV_OTHER_LOC
        INTEGER, DIMENSION( : ), intent( in ) :: MAT_OTHER_LOC
        REAL, DIMENSION( :, : ), intent( in ) :: SMATFEN
        REAL, DIMENSION( :, : ), intent( in ) :: SCVFEN
        REAL, DIMENSION( :, :, :, : ), intent( in ) :: TDIFFUSION
        REAL, DIMENSION( :, :, : ), intent( in ) :: DIFF_GI_ADDED
        REAL, DIMENSION( NDIM, NPHASE, CV_NLOC ), intent( in ) :: LOC_DTX_ELE_ALL, LOC_DTOLDX_ELE_ALL, LOC2_DTX_ELE_ALL, LOC2_DTOLDX_ELE_ALL
        REAL, DIMENSION( : ), intent( in ) :: CVNORMX_ALL

        ! local variables

        ! DIFF_MIN_FRAC is the fraction of the standard diffusion coefficient to use
        ! in the non-linear diffusion scheme. DIFF_MAX_FRAC is the maximum fraction.
        REAL, PARAMETER :: DIFF_MIN_FRAC = 0.05, DIFF_MAX_FRAC = 20.0

        REAL :: COEF
        INTEGER :: CV_KLOC, CV_KLOC2, MAT_KLOC, MAT_KLOC2, MAT_NODK, MAT_NODK2, IPHASE, CV_SKLOC
        LOGICAL :: ZER_DIFF

        REAL, DIMENSION ( :, : ), allocatable :: DTDX_GI_ALL, DTOLDDX_GI_ALL, DTDX_GI2_ALL, DTOLDDX_GI2_ALL
        REAL, DIMENSION ( : ), allocatable :: N_DOT_DKDT_ALL, N_DOT_DKDTOLD_ALL, N_DOT_DKDT2_ALL, N_DOT_DKDTOLD2_ALL
        REAL, DIMENSION ( : ), allocatable :: DIFF_STAND_DIVDX_ALL, DIFF_STAND_DIVDX2_ALL
        REAL, DIMENSION ( :, :, : ), allocatable :: DIFF_GI, DIFF_GI2

        ALLOCATE( DTDX_GI_ALL( NDIM, NPHASE ), DTOLDDX_GI_ALL( NDIM, NPHASE ) )
        ALLOCATE( DTDX_GI2_ALL( NDIM, NPHASE ), DTOLDDX_GI2_ALL( NDIM, NPHASE ) )
        ALLOCATE( N_DOT_DKDT_ALL( NPHASE ), N_DOT_DKDTOLD_ALL( NPHASE ) )
        ALLOCATE( N_DOT_DKDT2_ALL( NPHASE ), N_DOT_DKDTOLD2_ALL( NPHASE ) )
        ALLOCATE( DIFF_STAND_DIVDX_ALL( NPHASE ), DIFF_STAND_DIVDX2_ALL( NPHASE ) )

        ALLOCATE( DIFF_GI( NDIM, NDIM, NPHASE ) )
        ALLOCATE( DIFF_GI2( NDIM, NDIM, NPHASE ) )

        ZER_DIFF = .FALSE.
        IF ( on_domain_boundary ) ZER_DIFF = ANY ( LOC_WIC_T_BC( : ) /= WIC_T_BC_DIRICHLET )

        Cond_ZerDiff: IF ( ZER_DIFF ) THEN

            DIFF_COEF_DIVDX = 0.0
            DIFF_COEFOLD_DIVDX = 0.0

        ELSE

            DTDX_GI_ALL = 0.0 ; DTOLDDX_GI_ALL = 0.0
            DO CV_KLOC = 1, CV_NLOC
                DTDX_GI_ALL( :, : ) = DTDX_GI_ALL( :, : ) + SCVFEN( CV_KLOC, GI ) * LOC_DTX_ELE_ALL( :, :, CV_KLOC )
                DTOLDDX_GI_ALL( :, : ) = DTOLDDX_GI_ALL( :, : ) + SCVFEN( CV_KLOC, GI ) * LOC_DTOLDX_ELE_ALL( :, :, CV_KLOC )
            END DO

            DIFF_GI = 0.0
            DO MAT_KLOC = 1, MAT_NLOC
                MAT_NODK = MAT_NDGLN( ( ELE - 1 ) * MAT_NLOC + MAT_KLOC )
                DO IPHASE = 1, NPHASE
                    DIFF_GI( :, :, IPHASE ) = DIFF_GI( :, :, IPHASE ) &
                        + SMATFEN( MAT_KLOC, GI ) * TDIFFUSION( MAT_NODK, :, :, IPHASE )
                END DO
            END DO
            DIFF_GI = DIFF_GI + DIFF_GI_ADDED
            DIFF_GI = MAX( 0.0, DIFF_GI )

            DO IPHASE = 1, NPHASE
                N_DOT_DKDT_ALL( IPHASE ) = DOT_PRODUCT( CVNORMX_ALL, MATMUL( DIFF_GI( :, :, IPHASE ), DTDX_GI_ALL( :, IPHASE ) ) )
                N_DOT_DKDTOLD_ALL( IPHASE ) = DOT_PRODUCT( CVNORMX_ALL, MATMUL( DIFF_GI( :, :, IPHASE ), DTOLDDX_GI_ALL( :, IPHASE ) ) )

                COEF = DOT_PRODUCT( CVNORMX_ALL, MATMUL( DIFF_GI( :, :, IPHASE ), CVNORMX_ALL ) )
                DIFF_STAND_DIVDX_ALL( IPHASE ) = COEF / HDC
            END DO



            !       Conditional_MAT_DISOPT_ELE2: IF ( ( ELE2 /= 0 ) .AND. ( ELE2 /= ELE ) ) THEN
            Conditional_MAT_DISOPT_ELE2: IF ( between_elements ) THEN


                DTDX_GI2_ALL = 0.0 ; DTOLDDX_GI2_ALL = 0.0
                DO CV_SKLOC=1,CV_SNLOC
                    CV_KLOC=CV_SLOC2LOC( CV_SKLOC )
                    CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                    !             IF ( CV_KLOC2 /= 0 )THEN
                    DTDX_GI2_ALL( :, : ) = DTDX_GI2_ALL( :, : ) + SCVFEN( CV_KLOC, GI ) *  LOC2_DTX_ELE_ALL( :, :, CV_KLOC )
                    DTOLDDX_GI2_ALL( :, : ) = DTOLDDX_GI2_ALL( :, : ) + SCVFEN( CV_KLOC, GI ) * LOC2_DTOLDX_ELE_ALL( :, :, CV_KLOC )
                !             END IF
                END DO

                DIFF_GI2 = 0.0
                DO MAT_KLOC = 1, MAT_NLOC
                    MAT_KLOC2 = MAT_OTHER_LOC( MAT_KLOC )
                    IF ( MAT_KLOC2 /= 0 ) THEN
                        MAT_NODK2 = MAT_NDGLN( ( ELE2 - 1 ) * MAT_NLOC + MAT_KLOC2 )
                        DO IPHASE = 1, NPHASE
                            DIFF_GI2( :, :, IPHASE ) = DIFF_GI2( :, :, IPHASE ) &
                                + SMATFEN( MAT_KLOC, GI ) * TDIFFUSION( MAT_NODK2, :, :, IPHASE )
                        END DO
                    END IF
                END DO

                DO IPHASE = 1, NPHASE
                    N_DOT_DKDT2_ALL( IPHASE ) = DOT_PRODUCT( CVNORMX_ALL, MATMUL( DIFF_GI2( :, :, IPHASE ), DTDX_GI2_ALL( :, IPHASE ) ) )
                    N_DOT_DKDTOLD2_ALL( IPHASE ) = DOT_PRODUCT( CVNORMX_ALL, MATMUL( DIFF_GI2( :, :, IPHASE ), DTOLDDX_GI2_ALL( :, IPHASE ) ) )

                    COEF = DOT_PRODUCT( CVNORMX_ALL, MATMUL( DIFF_GI2( :, :, IPHASE ), CVNORMX_ALL ) )
                    DIFF_STAND_DIVDX2_ALL( IPHASE ) = COEF  /HDC
                END DO

                N_DOT_DKDT_ALL = 0.5 * ( N_DOT_DKDT_ALL + N_DOT_DKDT2_ALL )
                N_DOT_DKDTOLD_ALL = 0.5 * ( N_DOT_DKDTOLD_ALL + N_DOT_DKDTOLD2_ALL )

                DIFF_STAND_DIVDX_ALL = MIN( DIFF_STAND_DIVDX_ALL, DIFF_STAND_DIVDX2_ALL )

            END IF Conditional_MAT_DISOPT_ELE2


            DIFF_COEF_DIVDX = MAX( DIFF_MIN_FRAC * DIFF_STAND_DIVDX_ALL, N_DOT_DKDT_ALL / TOLFUN_MANY( T_CV_NODJ - T_CV_NODI ) )
            DIFF_COEFOLD_DIVDX = MAX( DIFF_MIN_FRAC * DIFF_STAND_DIVDX_ALL, N_DOT_DKDTOLD_ALL / TOLFUN_MANY( TOLD_CV_NODJ - TOLD_CV_NODI ) )

            DIFF_COEF_DIVDX = MIN( DIFF_MAX_FRAC * DIFF_STAND_DIVDX_ALL, DIFF_COEF_DIVDX )
            DIFF_COEFOLD_DIVDX = MIN( DIFF_MAX_FRAC * DIFF_STAND_DIVDX_ALL, DIFF_COEFOLD_DIVDX )

        END IF Cond_ZerDiff

        !    IF ( SELE /= 0 ) THEN
        IF ( on_domain_boundary ) THEN
            DO IPHASE=1,NPHASE
                IF( LOC_WIC_T_BC( IPHASE ) /= WIC_T_BC_DIRICHLET ) THEN
                    DIFF_COEF_DIVDX( IPHASE ) = 0.0
                    DIFF_COEFOLD_DIVDX( IPHASE ) = 0.0
                ENDIF
            END DO
        ENDIF


        DEALLOCATE( DIFF_GI, DIFF_GI2, &
            DTDX_GI_ALL, DTOLDDX_GI_ALL, &
            DTDX_GI2_ALL, DTOLDDX_GI2_ALL, &
            N_DOT_DKDT_ALL, N_DOT_DKDTOLD_ALL, &
            N_DOT_DKDT2_ALL, N_DOT_DKDTOLD2_ALL, &
            DIFF_STAND_DIVDX_ALL, DIFF_STAND_DIVDX2_ALL )

        RETURN

    END SUBROUTINE DIFFUS_CAL_COEFF


    SUBROUTINE LINEAR_HIGH_DIFFUS_CAL_COEFF_STRESS_OR_TENSOR( STRESS_IJ_ELE_EXT,  S_INV_NNX_MAT12,  &
        STRESS_FORM, STRESS_FORM_STAB, ZERO_OR_TWO_THIRDS, &
        U_SNLOC, U_NLOC, CV_SNLOC, NPHASE,  &
        SBUFEN_REVERSED,SBCVFEN_REVERSED, SDETWEI, SBCVNGI, NDIM, SLOC_UDIFFUSION, SLOC_UDIFFUSION_VOL, SLOC2_UDIFFUSION, SLOC2_UDIFFUSION_VOL, DIFF_GI_ADDED, &
        ON_BOUNDARY, SNORMXN_ALL  )
        ! This sub calculates the effective diffusion coefficientd STRESS_IJ_ELE_EXT
        ! it only works for between element contributions.
        ! based on a high order scheme.
        ! The matrix  S_INV_NNX_MAT12 is used to calculate the rows of the matrix with STRESS_IJ_ELE_EXT.
        ! This implements the stress and tensor form of diffusion and calculates a jump conidition.
        ! which is in DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX
        ! The coefficient are in N_DOT_DKDU, N_DOT_DKDUOLD.
        ! look at the manual DG treatment of viscocity.
        IMPLICIT NONE
        LOGICAL, intent( in ) :: STRESS_FORM, STRESS_FORM_STAB, ON_BOUNDARY
        INTEGER, intent( in ) :: U_SNLOC, U_NLOC, CV_SNLOC, NPHASE,  &
            &                   SBCVNGI, NDIM
        REAL, intent( in ) :: ZERO_OR_TWO_THIRDS
        REAL, DIMENSION( NDIM, NDIM, NPHASE, U_SNLOC, 2*U_NLOC ), intent( inout ) :: STRESS_IJ_ELE_EXT
        REAL, DIMENSION( NDIM, U_SNLOC, 2*U_NLOC ), intent( inout ) :: S_INV_NNX_MAT12
        REAL, DIMENSION( SBCVNGI, CV_SNLOC ), intent( in ) :: SBCVFEN_REVERSED
        REAL, DIMENSION( SBCVNGI, U_SNLOC  ), intent( in ) :: SBUFEN_REVERSED
        REAL, DIMENSION( SBCVNGI ), intent( in ) :: SDETWEI
        REAL, DIMENSION( NDIM,NDIM,NPHASE,CV_SNLOC ), intent( in ) :: SLOC_UDIFFUSION, SLOC2_UDIFFUSION
        REAL, DIMENSION( NPHASE,CV_SNLOC ), intent( in ) :: SLOC_UDIFFUSION_VOL, SLOC2_UDIFFUSION_VOL
        ! DIFF_GI_ADDED( IDIM, :,:) is for dimension IDIM e.g IDIM=1 corresponds to U
        ! the rest is for the diffusion tensor.
        REAL, DIMENSION( NDIM, NDIM,NDIM, NPHASE, SBCVNGI), intent( in ) :: DIFF_GI_ADDED
        REAL, DIMENSION( NDIM, SBCVNGI ), intent( in ) :: SNORMXN_ALL
        ! local variables...
        ! If FAST then use the fast version that is less well tested..
        LOGICAL, PARAMETER :: FAST = .true.

        REAL, DIMENSION(NDIM,NDIM,NPHASE,SBCVNGI) :: DIFF_GI, DIFF_GI2
        REAL, DIMENSION(NPHASE,SBCVNGI) :: DIFF_VOL_GI, DIFF_VOL_GI2
        ! for rapid version of code...
        REAL, DIMENSION(NDIM,NDIM,NPHASE,U_SNLOC,U_SNLOC) :: MAT_SUFXX
        REAL, DIMENSION(NDIM,NPHASE,U_SNLOC,U_SNLOC) :: MAT_SUFVOL

        !    REAL, DIMENSION(NDIM,Ndim) :: two_dim
        !    REAL, DIMENSION(NDIM) :: one_dim

        INTEGER :: IDIM,JDIM,CV_SKLOC
        INTEGER :: SGI,IPHASE,U_JLOC12,I,U_SILOC,U_SJLOC
        REAL :: RDUM


        !          ALLOCATE( DIFF_GI(NDIM,NDIM,NPHASE,SBCVNGI) )
        !          ALLOCATE( DIFF_GI2(NDIM,NDIM,NPHASE,SBCVNGI) )

        !          ALLOCATE( DIFF_VOL_GI(NPHASE,SBCVNGI) )
        !          ALLOCATE( DIFF_VOL_GI2(NPHASE,SBCVNGI) )

        DIFF_GI = 0.0
        DIFF_VOL_GI = 0.0
        DO CV_SKLOC = 1, CV_SNLOC
            DO SGI=1,SBCVNGI
                DO IPHASE=1, NPHASE
                    DIFF_GI( 1:NDIM , 1:NDIM, IPHASE, SGI ) = DIFF_GI( 1:NDIM , 1:NDIM, IPHASE, SGI ) &
                        + SBCVFEN_REVERSED(SGI,CV_SKLOC) * SLOC_UDIFFUSION( 1:NDIM , 1:NDIM , IPHASE, CV_SKLOC )

                    DIFF_VOL_GI( IPHASE, SGI ) = DIFF_VOL_GI( IPHASE, SGI ) &
                        + SBCVFEN_REVERSED(SGI,CV_SKLOC) * SLOC_UDIFFUSION_VOL( IPHASE, CV_SKLOC )
                END DO
            END DO
        END DO
        DIFF_GI=MAX(0.0, DIFF_GI)
        DIFF_VOL_GI=MAX(0.0, DIFF_VOL_GI)

        IF(.NOT.ON_BOUNDARY) THEN
            ! neighbouring element...
            DIFF_GI2 = 0.0
            DIFF_VOL_GI2 = 0.0
            DO CV_SKLOC = 1, CV_SNLOC
                DO SGI=1,SBCVNGI
                    DO IPHASE=1, NPHASE
                        DIFF_GI2( 1:NDIM, 1:NDIM, IPHASE, SGI )= DIFF_GI2( 1:NDIM, 1:NDIM, IPHASE, SGI ) +SBCVFEN_REVERSED(SGI,CV_SKLOC) &
                            *SLOC2_UDIFFUSION(1:NDIM, 1:NDIM ,IPHASE, CV_SKLOC)

                        DIFF_VOL_GI2( IPHASE, SGI )= DIFF_VOL_GI2( IPHASE, SGI ) +SBCVFEN_REVERSED(SGI,CV_SKLOC) &
                            *SLOC2_UDIFFUSION_VOL(IPHASE, CV_SKLOC)
                    END DO
                END DO
            END DO
            DIFF_GI2=MAX(0.0, DIFF_GI2)
            DIFF_VOL_GI2=MAX(0.0, DIFF_VOL_GI2)

            DIFF_GI=0.5*(DIFF_GI+DIFF_GI2)
            DIFF_VOL_GI=0.5*(DIFF_VOL_GI+DIFF_VOL_GI2)
        ENDIF


        IF(STRESS_FORM_STAB) THEN

            DO JDIM=1,NDIM
                DO IDIM=1,NDIM
                    DIFF_GI(IDIM, JDIM, :, :) = DIFF_GI(IDIM, JDIM, :, :) &
                        + SQRT( DIFF_GI_ADDED(IDIM, 1,1, :, :) * DIFF_GI_ADDED(JDIM, 1,1, :, :) )

                END DO
            END DO

        ELSE
            DO IDIM=1,NDIM
                DIFF_GI(IDIM, IDIM, :, :) = DIFF_GI(IDIM, IDIM, :, :) &
                    +  DIFF_GI_ADDED(IDIM, 1,1, :, :)
            END DO
        ENDIF



        STRESS_IJ_ELE_EXT=0.0


        IF(STRESS_FORM) THEN
            IF(FAST) THEN ! The rapid version
                MAT_SUFXX=0.0
                MAT_SUFVOL=0.0

                DO U_SILOC=1,U_SNLOC
                    DO U_SJLOC=1,U_SNLOC
                        DO SGI=1,SBCVNGI
                            DO JDIM=1,NDIM
                                RDUM = SNORMXN_ALL(JDIM,SGI)*SBUFEN_REVERSED(SGI,U_SILOC)*SBUFEN_REVERSED(SGI,U_SJLOC)*SDETWEI( SGI )
                                DO IPHASE=1,NPHASE
                                    DO IDIM=1,NDIM
                                        ! take -ve as its a surface integral...
                                        MAT_SUFXX(IDIM,JDIM,IPHASE,U_SILOC,U_SJLOC) = MAT_SUFXX(IDIM,JDIM,IPHASE,U_SILOC,U_SJLOC) &
                                            -  RDUM*DIFF_GI( IDIM, JDIM, IPHASE, SGI )
                                    END DO
                                    MAT_SUFVOL(JDIM,IPHASE,U_SILOC,U_SJLOC) = MAT_SUFVOL(JDIM,IPHASE,U_SILOC,U_SJLOC) &
                                        -  RDUM*DIFF_VOL_GI( IPHASE, SGI )
                                END DO
                            END DO
                        END DO
                    END DO
                END DO



                DO I=1,U_SNLOC
                    DO U_JLOC12=1,U_NLOC*2
                        DO U_SILOC=1,U_SNLOC
                            DO IPHASE=1,NPHASE
                                CALL CALC_STRESS_TEN_REDUCE(STRESS_IJ_ELE_EXT( :, :, IPHASE, U_SILOC, U_JLOC12 ), ZERO_OR_TWO_THIRDS, NDIM,    &
                                    MAT_SUFXX(:,:,IPHASE,U_SILOC,I), MAT_SUFVOL(:,IPHASE,U_SILOC,I),  S_INV_NNX_MAT12( :, I, U_JLOC12 )  )

                            END DO
                        END DO

                    END DO

                END DO
            ELSE ! less rapid version...
                DO U_JLOC12=1,U_NLOC*2
                    DO U_SILOC=1,U_SNLOC
                        DO I=1,U_SNLOC
                            DO SGI=1,SBCVNGI
                                DO IPHASE=1,NPHASE
                                    ! take -ve as its a surface integral...
                                    CALL CALC_STRESS_TEN( STRESS_IJ_ELE_EXT( :, :, IPHASE, U_SILOC, U_JLOC12 ), ZERO_OR_TWO_THIRDS, NDIM, &
                                        - SNORMXN_ALL(:,SGI)*SBUFEN_REVERSED(SGI,U_SILOC)* SBUFEN_REVERSED(SGI,I)*SDETWEI( SGI ), S_INV_NNX_MAT12( 1:NDIM, I, U_JLOC12 ), DIFF_GI( :, :, IPHASE, SGI ), DIFF_VOL_GI( IPHASE, SGI) )
                                END DO
                            END DO
                        END DO

                    END DO
                END DO
            ENDIF

        ELSE ! tensor form of viscocity...
            IF(FAST) THEN ! The rapid version

                MAT_SUFXX=0.0

                DO U_SILOC=1,U_SNLOC
                    DO U_SJLOC=1,U_SNLOC
                        DO SGI=1,SBCVNGI
                            DO IDIM=1,NDIM
                                RDUM = SNORMXN_ALL(IDIM,SGI)*SBUFEN_REVERSED(SGI,U_SILOC)*SBUFEN_REVERSED(SGI,U_SJLOC)*SDETWEI( SGI )
                                DO IPHASE=1,NPHASE
                                    DO JDIM=1,NDIM
                                        ! take -ve as its a surface integral...
                                        MAT_SUFXX(IDIM,JDIM,IPHASE,U_SILOC,U_SJLOC) = MAT_SUFXX(IDIM,JDIM,IPHASE,U_SILOC,U_SJLOC) &
                                            -  RDUM*DIFF_GI( IDIM, JDIM, IPHASE, SGI )
                                    END DO
                                END DO
                            END DO
                        END DO
                    END DO
                END DO
                ! STRESS_IJ_ELE_EXT = MAT_SUFX*S_INV_NNX_MAT12
                DO U_JLOC12=1,U_NLOC*2
                    DO U_SILOC=1,U_SNLOC
                        DO I=1,U_SNLOC
                            DO IPHASE=1,NPHASE
                                DO IDIM=1,NDIM
                                    STRESS_IJ_ELE_EXT( 1,1, IPHASE, U_SILOC, U_JLOC12 ) = STRESS_IJ_ELE_EXT( 1,1, IPHASE, U_SILOC, U_JLOC12 ) &
                                        + SUM(  MAT_SUFXX(IDIM,:,IPHASE, U_SILOC, I) *  S_INV_NNX_MAT12( :, I, U_JLOC12 )  )
                                END DO
                            END DO
                        END DO

                    END DO
                END DO

            ELSE  ! The rapid version

                DO SGI=1,SBCVNGI
                    DO U_JLOC12=1,U_NLOC*2
                        DO U_SILOC=1,U_SNLOC

                            DO I=1,U_SNLOC
                                DO IPHASE=1,NPHASE
                                    DO IDIM=1,NDIM
                                        STRESS_IJ_ELE_EXT( 1,1, IPHASE, U_SILOC, U_JLOC12 ) = STRESS_IJ_ELE_EXT( 1,1, IPHASE, U_SILOC, U_JLOC12 ) &
                                            - SNORMXN_ALL(IDIM,SGI)*SBUFEN_REVERSED(SGI,U_SILOC)* SBUFEN_REVERSED(SGI,I)*SDETWEI( SGI )  &
                                            * SUM( DIFF_GI( IDIM, :, IPHASE, SGI ) *  S_INV_NNX_MAT12( :, I, U_JLOC12 ) )
                                    !                                 * DIFF_GI( IDIM, idim, IPHASE, SGI ) *  S_INV_NNX_MAT12( idim, I, U_JLOC12 )
                                    END DO
                                END DO
                            END DO
                        END DO

                    END DO
                END DO
            ENDIF ! endof if THEN ELSE OF The rapid version
            ! all other components are the same...
            DO IDIM=2,NDIM
                STRESS_IJ_ELE_EXT( IDIM,IDIM, :,:,: ) = STRESS_IJ_ELE_EXT( 1,1, :,:,: )
            END DO

        ENDIF

        !          CALL CALC_STRESS_TEN( STRESS_IJ_ELE( :, :, IPHASE, U_SILOC, U_JLOC_EXT ), ZERO_OR_TWO_THIRDS, NDIM, &
        !               UFENX_ALL( 1:NDIM, U_ILOC, GI ), UFENX_ALL( 1:NDIM, U_JLOC, GI )* DETWEI( GI ), TEN_XX( :, :, IPHASE, GI ), TEN_VOL( IPHASE, GI) )


        RETURN

    END SUBROUTINE LINEAR_HIGH_DIFFUS_CAL_COEFF_STRESS_OR_TENSOR











    SUBROUTINE CALC_STRESS_TEN(STRESS_IJ, ZERO_OR_TWO_THIRDS, NDIM,    &
        UFENX_ILOC, UFENX_JLOC,  TEN_XX, TEN_VOL )
        ! determine stress form of viscocity...
        IMPLICIT NONE
        INTEGER, intent( in )  :: NDIM
        REAL, DIMENSION( :, :  ), intent( inOUT ) :: STRESS_IJ
        REAL, DIMENSION( : ), intent( in ) :: UFENX_ILOC
        REAL, DIMENSION( :,: ), intent( in ) :: TEN_XX
        ! TEN_VOL is volumetric viscocity - mostly set to zero other than q-scheme or use with kinetic theory
        REAL, intent( in ) :: TEN_VOL
        REAL, DIMENSION( : ), intent( in ) :: UFENX_JLOC
        REAL, intent( in ) :: ZERO_OR_TWO_THIRDS
        ! Local variables...
        REAL :: FEN_TEN_XX(NDIM,NDIM),FEN_TEN_VOL(NDIM)
        INTEGER :: IDIM,JDIM


        DO IDIM=1,NDIM
            FEN_TEN_XX(IDIM,:)=UFENX_ILOC(:) * TEN_XX(IDIM,:)
        END DO
        FEN_TEN_VOL(:)=UFENX_ILOC(:) * TEN_VOL

        DO IDIM=1,NDIM
            STRESS_IJ( IDIM,IDIM ) = STRESS_IJ( IDIM,IDIM ) + SUM( FEN_TEN_XX(IDIM,:) * UFENX_JLOC(:) )
        END DO

        DO JDIM=1,NDIM
            DO IDIM=1,NDIM
                !               STRESS_IJ( IDIM,JDIM ) = STRESS_IJ( IDIM,JDIM ) + FEN_TEN_XX(IDIM,JDIM) * UFENX_JLOC(JDIM)
                STRESS_IJ( IDIM,JDIM ) = STRESS_IJ( IDIM,JDIM ) + FEN_TEN_XX(IDIM,JDIM) * UFENX_JLOC(IDIM)
            END DO
        END DO

        DO JDIM=1,NDIM
            DO IDIM=1,NDIM
                STRESS_IJ( IDIM,JDIM ) = STRESS_IJ( IDIM,JDIM ) &
                    - ZERO_OR_TWO_THIRDS * FEN_TEN_XX(IDIM,IDIM) * UFENX_JLOC(JDIM) &
                    + FEN_TEN_VOL(IDIM) * UFENX_JLOC(JDIM)
            END DO
        END DO

        RETURN

    END SUBROUTINE CALC_STRESS_TEN




    SUBROUTINE CALC_STRESS_TEN_REDUCE(STRESS_IJ, ZERO_OR_TWO_THIRDS, NDIM,    &
        FEN_TEN_XX, FEN_TEN_VOL,  UFENX_JLOC  )
        ! determine stress form of viscocity...
        IMPLICIT NONE
        INTEGER, intent( in )  :: NDIM
        REAL, DIMENSION( NDIM, NDIM  ), intent( inOUT ) :: STRESS_IJ
        REAL, DIMENSION( NDIM, NDIM ), intent( in ) :: FEN_TEN_XX
        REAL, DIMENSION( NDIM ), intent( in ) :: FEN_TEN_VOL
        ! TEN_VOL is volumetric viscocity - mostly set to zero other than q-scheme or use with kinetic theory
        REAL, DIMENSION( NDIM ), intent( in ) :: UFENX_JLOC
        REAL, intent( in ) :: ZERO_OR_TWO_THIRDS
        ! Local variables...
        INTEGER :: IDIM,JDIM

        DO IDIM=1,NDIM
            STRESS_IJ( IDIM,IDIM ) = STRESS_IJ( IDIM,IDIM ) + SUM( FEN_TEN_XX(IDIM,:) * UFENX_JLOC(:) )
        END DO

        DO JDIM=1,NDIM
            DO IDIM=1,NDIM
                !               STRESS_IJ( IDIM,JDIM ) = STRESS_IJ( IDIM,JDIM ) + FEN_TEN_XX(IDIM,JDIM) * UFENX_JLOC(JDIM)
                STRESS_IJ( IDIM,JDIM ) = STRESS_IJ( IDIM,JDIM ) + FEN_TEN_XX(IDIM,JDIM) * UFENX_JLOC(IDIM)
            END DO
        END DO

        DO JDIM=1,NDIM
            DO IDIM=1,NDIM
                STRESS_IJ( IDIM,JDIM ) = STRESS_IJ( IDIM,JDIM ) &
                    - ZERO_OR_TWO_THIRDS * FEN_TEN_XX(IDIM,IDIM) * UFENX_JLOC(JDIM) &
                    + FEN_TEN_VOL(IDIM) * UFENX_JLOC(JDIM)
            END DO
        END DO

        RETURN

    END SUBROUTINE CALC_STRESS_TEN_REDUCE




    SUBROUTINE CAL_LIM_VOL_ADJUST( TMIN_STORE, TMIN, T, TMIN_NOD, RESET_STORE, MASS_CV, &
        &                         CV_NODI_IPHA, CV_NODJ_IPHA, IPHASE, CV_NONODS, INCOME )
        ! Adjust TMIN to take into account different sized CV's.
        ! if RESET_STORE then reset TMIN to orginal values.
        implicit none
        REAL, intent( in ) :: INCOME
        INTEGER, intent( in ) :: CV_NODI_IPHA, CV_NODJ_IPHA,IPHASE, CV_NONODS
        LOGICAL, intent( in ) :: RESET_STORE
        REAL, intent( inout ) :: TMIN_STORE
        REAL, DIMENSION( : ), intent( inout ) :: TMIN
        REAL, DIMENSION( : ), intent( in ) :: T
        INTEGER, DIMENSION( : ), intent( in ) :: TMIN_NOD
        REAL, DIMENSION( : ), intent( in ) :: MASS_CV
        ! Local variables...
        REAL DX1, DX2_MIN, COEFF
        INTEGER CV_NODI, CV_NODJ
        LOGICAL, PARAMETER :: DEF_BOUNDED = .FALSE.

        CV_NODI = CV_NODI_IPHA-(IPHASE-1)*CV_NONODS
        CV_NODJ = CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS

        IF ( RESET_STORE ) THEN

            IF ( INCOME < 0.5 ) THEN
                TMIN( CV_NODI_IPHA ) = TMIN_STORE
            ELSE
                TMIN( CV_NODJ_IPHA ) = TMIN_STORE
            ENDIF

        ELSE

            DX1 = 0.5 * ( MASS_CV(CV_NODI) + MASS_CV(CV_NODJ) )

            IF ( INCOME < 0.5 ) THEN
                DX2_MIN = 0.5 * ( MASS_CV( TMIN_NOD(CV_NODI_IPHA)-(IPHASE-1)*CV_NONODS) &
                    + MASS_CV( CV_NODI) )
                TMIN_STORE = TMIN(CV_NODI_IPHA)
                IF ( DEF_BOUNDED ) THEN ! This produces strictly bounded always soln
                    COEFF = MIN( 1.0, (DX1/DX2_MIN) )
                ELSE
                    COEFF = (DX1/DX2_MIN)
                ENDIF
                TMIN( CV_NODI_IPHA ) = T(CV_NODI_IPHA) + COEFF * ( TMIN_STORE - T(CV_NODI_IPHA) )
            ELSE
                DX2_MIN = 0.5 * ( MASS_CV( TMIN_NOD(CV_NODJ_IPHA)-(IPHASE-1)*CV_NONODS) &
                    + MASS_CV( CV_NODJ))
                TMIN_STORE = TMIN(CV_NODJ_IPHA)
                IF ( DEF_BOUNDED ) THEN ! This produces strictly bounded always soln
                    COEFF = MIN( 1.0, (DX1/DX2_MIN) )
                ELSE
                    COEFF = (DX1/DX2_MIN)
                END IF
                TMIN(CV_NODJ_IPHA) = T(CV_NODJ_IPHA) + COEFF * ( TMIN_STORE - T(CV_NODJ_IPHA) )
            END IF

        END IF

        RETURN
    END SUBROUTINE CAL_LIM_VOL_ADJUST

    SUBROUTINE DGSIMPLNORM_ALL( NLOC, SNLOC, NDIM,  &
        XL_ALL, XSL_ALL, NORMX_ALL )
        ! Form approximate surface normal (NORMX_ALL(1),NORMX_ALL(2),NORMX_ALL(3))
        IMPLICIT NONE
        INTEGER, intent( in ) :: NLOC, SNLOC, NDIM
        REAL, DIMENSION( NDIM, NLOC ), intent( in ) :: XL_ALL
        REAL, DIMENSION( NDIM, SNLOC ), intent( in ) :: XSL_ALL
        REAL, DIMENSION( NDIM ), intent( inout ) :: NORMX_ALL
        ! Local variables
        REAL :: XC(NDIM), SXC(NDIM), NORM
        INTEGER :: IDIM

        DO IDIM = 1, NDIM
            XC(IDIM) = SUM( XL_ALL( IDIM, : ) )

            SXC(IDIM) = SUM( XSL_ALL( IDIM, : ) )
        END DO

        NORMX_ALL = SXC / REAL( SNLOC )  - XC /  REAL( NLOC )

        NORM = SQRT( SUM( NORMX_ALL**2 ) )

        NORMX_ALL = NORMX_ALL / NORM

        RETURN

    END SUBROUTINE DGSIMPLNORM_ALL



    !sprint_to_do where this is being called use the new one with the new memory
    !and then remove
    SUBROUTINE DGSIMPLNORM( ELE, SILOC2ILOC, NLOC, SNLOC, XONDGL, &
        X, Y, Z, NORMX, NORMY, NORMZ )
        ! Form approximate surface normal (NORMX,NORMY,NORMZ)
        IMPLICIT NONE
        INTEGER, intent( in ) :: ELE, NLOC, SNLOC
        INTEGER, DIMENSION( : ), intent( in ) ::  SILOC2ILOC
        INTEGER, DIMENSION( : ), intent( in ) :: XONDGL
        REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
        REAL, intent( inout ) :: NORMX, NORMY, NORMZ
        ! Local variables
        REAL :: XC, YC, ZC, SXC, SYC, SZC, NORM
        INTEGER :: XNODI, ILOC, SILOC

        XC = 0.0
        YC = 0.0
        ZC = 0.0
        DO ILOC = 1, NLOC
            XNODI = XONDGL(( ELE - 1 ) * NLOC + ILOC )
            XC = XC + X( XNODI ) / REAL( NLOC )
            YC = YC + Y( XNODI ) / REAL( NLOC )
            ZC = ZC + Z( XNODI ) / REAL( NLOC )
        END DO

        SXC = 0.0
        SYC = 0.0
        SZC = 0.0
        DO SILOC = 1, SNLOC
            ILOC = SILOC2ILOC( SILOC )
            XNODI = XONDGL(( ELE - 1 ) * NLOC+ ILOC )
            SXC = SXC + X( XNODI ) / REAL( SNLOC )
            SYC = SYC + Y( XNODI ) / REAL( SNLOC )
            SZC = SZC + Z( XNODI ) / REAL( SNLOC )
        END DO
        NORMX = SXC - XC
        NORMY = SYC - YC
        NORMZ = SZC - ZC

        NORM = SQRT( NORMX**2 + NORMY**2 + NORMZ**2 )

        NORMX = NORMX / NORM
        NORMY = NORMY / NORM
        NORMZ = NORMZ / NORM

        RETURN

    END SUBROUTINE DGSIMPLNORM


    !sprint_to_do!does this work???
    SUBROUTINE ISOTROPIC_LIMITER_ALL( &
        ! FOR SUB SURRO_CV_MINMAX:
        T_ALL, TOLD_ALL, T2_ALL, T2OLD_ALL, DEN_ALL, DENOLD_ALL, IGOT_T2, NPHASE, CV_NONODS, nsmall_colm, SMALL_CENTRM, SMALL_FINDRM, SMALL_COLM, &
        STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC_ALL, SUF_T2_BC_ALL, SUF_D_BC_ALL, WIC_T_BC_ALL, WIC_T2_BC_ALL, WIC_D_BC_ALL, &
        MASS_CV, &
        ! FOR SUB CALC_LIMIT_MATRIX_MAX_MIN:
        TOLDUPWIND_MAT_ALL, DENOLDUPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL, &
        TUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, T2UPWIND_MAT_ALL )
        INTEGER, intent( in ) :: IGOT_T2, NPHASE, CV_NONODS, nsmall_colm, STOTEL, CV_SNLOC
        REAL, DIMENSION( :, : ), intent( in ) :: T_ALL, TOLD_ALL, T2_ALL, T2OLD_ALL, DEN_ALL, DENOLD_ALL
        REAL, DIMENSION( : ), intent( in ) :: MASS_CV
        REAL, DIMENSION( :, : ), intent( inout ) :: TOLDUPWIND_MAT_ALL, DENOLDUPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL, &
            TUPWIND_MAT_ALL,    DENUPWIND_MAT_ALL,    T2UPWIND_MAT_ALL
        INTEGER, DIMENSION( : ), intent( in ) :: SMALL_CENTRM, SMALL_FINDRM
        INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
        INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
        REAL, DIMENSION( :, :, : ), intent( in ), pointer :: SUF_T_BC_ALL, SUF_T2_BC_ALL, SUF_D_BC_ALL
        INTEGER, DIMENSION( :,:, : ), intent( in ) :: WIC_T_BC_ALL, WIC_T2_BC_ALL, WIC_D_BC_ALL

        ! Local variables...
        REAL, DIMENSION( :, : ), allocatable :: TMIN_ALL, TMAX_ALL, TOLDMIN_ALL, TOLDMAX_ALL, &
            T2MIN_ALL, T2MAX_ALL, T2OLDMIN_ALL, T2OLDMAX_ALL, DENMIN_ALL, DENMAX_ALL,  &
            DENOLDMIN_ALL, DENOLDMAX_ALL
        INTEGER, DIMENSION( :, : ), allocatable ::TMIN_NOD_ALL, TMAX_NOD_ALL, TOLDMIN_NOD_ALL, TOLDMAX_NOD_ALL, &
            T2MIN_NOD_ALL, T2MAX_NOD_ALL, T2OLDMIN_NOD_ALL, T2OLDMAX_NOD_ALL, &
            DENMIN_NOD_ALL, DENMAX_NOD_ALL, DENOLDMIN_NOD_ALL, DENOLDMAX_NOD_ALL
        INTEGER :: CV_NODI, IMID, IPHASE

        ! Allocate memory for terms needed by GETGXYZ OR ONVDLIM
        ALLOCATE(      TMIN_ALL( NPHASE, CV_NONODS) )
        ALLOCATE(      TMAX_ALL( NPHASE, CV_NONODS) )
        ALLOCATE(   TOLDMIN_ALL( NPHASE, CV_NONODS) )
        ALLOCATE(   TOLDMAX_ALL( NPHASE, CV_NONODS) )
        ALLOCATE(      T2MIN_ALL( NPHASE, CV_NONODS* IGOT_T2) )
        ALLOCATE(      T2MAX_ALL( NPHASE, CV_NONODS* IGOT_T2) )
        ALLOCATE(   T2OLDMIN_ALL( NPHASE, CV_NONODS* IGOT_T2) )
        ALLOCATE(   T2OLDMAX_ALL( NPHASE, CV_NONODS* IGOT_T2) )
        ALLOCATE(    DENMIN_ALL( NPHASE, CV_NONODS) )
        ALLOCATE(    DENMAX_ALL( NPHASE, CV_NONODS) )
        ALLOCATE( DENOLDMIN_ALL( NPHASE, CV_NONODS) )
        ALLOCATE( DENOLDMAX_ALL( NPHASE, CV_NONODS) )

        ALLOCATE(      TMIN_NOD_ALL( NPHASE, CV_NONODS) )
        ALLOCATE(      TMAX_NOD_ALL( NPHASE, CV_NONODS) )
        ALLOCATE(   TOLDMIN_NOD_ALL( NPHASE, CV_NONODS) )
        ALLOCATE(   TOLDMAX_NOD_ALL( NPHASE, CV_NONODS) )
        ALLOCATE(      T2MIN_NOD_ALL( NPHASE, CV_NONODS* IGOT_T2) )
        ALLOCATE(      T2MAX_NOD_ALL( NPHASE, CV_NONODS* IGOT_T2) )
        ALLOCATE(   T2OLDMIN_NOD_ALL( NPHASE, CV_NONODS* IGOT_T2) )
        ALLOCATE(   T2OLDMAX_NOD_ALL( NPHASE, CV_NONODS* IGOT_T2) )
        ALLOCATE(    DENMIN_NOD_ALL( NPHASE, CV_NONODS) )
        ALLOCATE(    DENMAX_NOD_ALL( NPHASE, CV_NONODS) )
        ALLOCATE( DENOLDMIN_NOD_ALL( NPHASE, CV_NONODS) )
        ALLOCATE( DENOLDMAX_NOD_ALL( NPHASE, CV_NONODS) )

        ! For each node, find the largest and smallest value of T and
        ! DENSITY for both the current and previous timestep, out of
        ! the node value and all its surrounding nodes including Dirichlet b.c's.
        CALL SURRO_CV_MINMAX( TMAX_ALL, TMIN_ALL, TOLDMAX_ALL, TOLDMIN_ALL, DENMAX_ALL, DENMIN_ALL, DENOLDMAX_ALL, DENOLDMIN_ALL, &
            T2MAX_ALL, T2MIN_ALL, T2OLDMAX_ALL, T2OLDMIN_ALL, &
            T_ALL, TOLD_ALL, T2_ALL, T2OLD_ALL, DEN_ALL, DENOLD_ALL, IGOT_T2, NPHASE, CV_NONODS, SMALL_FINDRM, SMALL_COLM, &
            STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC_ALL, SUF_T2_BC_ALL, SUF_D_BC_ALL, WIC_T_BC_ALL, WIC_T2_BC_ALL, WIC_D_BC_ALL, &
            TMIN_NOD_ALL, TMAX_NOD_ALL, TOLDMIN_NOD_ALL, TOLDMAX_NOD_ALL, &
            T2MIN_NOD_ALL, T2MAX_NOD_ALL, T2OLDMIN_NOD_ALL, T2OLDMAX_NOD_ALL, &
            DENMIN_NOD_ALL, DENMAX_NOD_ALL, DENOLDMIN_NOD_ALL, DENOLDMAX_NOD_ALL )

        ! Populate  limiting matrix based on max and min values OLD:
        CALL CALC_LIMIT_MATRIX_MAX_MIN(TOLDMAX_ALL, TOLDMIN_ALL, DENOLDMAX_ALL, DENOLDMIN_ALL, &
            T2OLDMAX_ALL, T2OLDMIN_ALL,  &
            TOLD_ALL,  T2OLD_ALL, DENOLD_ALL, IGOT_T2, NPHASE, CV_NONODS, &
            TOLDMIN_NOD_ALL, TOLDMAX_NOD_ALL,  &
            T2OLDMIN_NOD_ALL, T2OLDMAX_NOD_ALL,  &
            DENOLDMIN_NOD_ALL, DENOLDMAX_NOD_ALL,  &
            NSMALL_COLM, SMALL_FINDRM, SMALL_COLM,   &
            TOLDUPWIND_MAT_ALL, DENOLDUPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL, MASS_CV )

        ! Populate  limiting matrix based on max and min values
        CALL CALC_LIMIT_MATRIX_MAX_MIN(TMAX_ALL, TMIN_ALL, DENMAX_ALL, DENMIN_ALL, &
            T2MAX_ALL, T2MIN_ALL,  &
            T_ALL,  T2_ALL, DEN_ALL, IGOT_T2, NPHASE, CV_NONODS, &
            TMIN_NOD_ALL, TMAX_NOD_ALL,  &
            T2MIN_NOD_ALL, T2MAX_NOD_ALL,  &
            DENMIN_NOD_ALL, DENMAX_NOD_ALL,  &
            NSMALL_COLM, SMALL_FINDRM, SMALL_COLM,   &
            TUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, T2UPWIND_MAT_ALL, MASS_CV )

        ! make sure the diagonal is equal to the value:
        DO CV_NODI=1,CV_NONODS
            IMID=SMALL_CENTRM(CV_NODI)
            DO IPHASE=1,NPHASE
                TUPWIND_MAT_ALL( IPHASE, IMID)=T_ALL( IPHASE, CV_NODI)
                TOLDUPWIND_MAT_ALL(  IPHASE, IMID)=TOLD_ALL( IPHASE, CV_NODI)

                DENUPWIND_MAT_ALL(  IPHASE, IMID)=DEN_ALL( IPHASE, CV_NODI)
                DENOLDUPWIND_MAT_ALL(  IPHASE, IMID)=DENOLD_ALL( IPHASE, CV_NODI)

                IF( IGOT_T2 == 1 ) THEN
                    T2UPWIND_MAT_ALL(IPHASE, IMID)=T2_ALL( IPHASE, CV_NODI)
                    T2OLDUPWIND_MAT_ALL(IPHASE, IMID)=T2OLD_ALL( IPHASE, CV_NODI)
                ENDIF
            END DO
        END DO

        DEALLOCATE( TMIN_ALL, TMAX_ALL, TOLDMIN_ALL, TOLDMAX_ALL, &
            T2MIN_ALL, T2MAX_ALL, T2OLDMIN_ALL, T2OLDMAX_ALL, DENMIN_ALL, DENMAX_ALL,  &
            DENOLDMIN_ALL, DENOLDMAX_ALL )
        DEALLOCATE( TMIN_NOD_ALL, TMAX_NOD_ALL, TOLDMIN_NOD_ALL, TOLDMAX_NOD_ALL, &
            T2MIN_NOD_ALL, T2MAX_NOD_ALL, T2OLDMIN_NOD_ALL, T2OLDMAX_NOD_ALL, &
            DENMIN_NOD_ALL, DENMAX_NOD_ALL, DENOLDMIN_NOD_ALL, DENOLDMAX_NOD_ALL )


        contains
            SUBROUTINE SURRO_CV_MINMAX( TMAX_ALL, TMIN_ALL, TOLDMAX_ALL, TOLDMIN_ALL, DENMAX_ALL, DENMIN_ALL, DENOLDMAX_ALL, DENOLDMIN_ALL, &
                T2MAX_ALL, T2MIN_ALL, T2OLDMAX_ALL, T2OLDMIN_ALL, &
                T_ALL, TOLD_ALL,  T2_ALL, T2OLD_ALL, DEN_ALL, DENOLD_ALL, IGOT_T2, NPHASE, CV_NONODS, FINACV, COLACV, &
                STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC_ALL,  SUF_T2_BC_ALL, SUF_D_BC_ALL, WIC_T_BC_ALL, WIC_T2_BC_ALL, WIC_D_BC_ALL, &
                TMIN_NOD_ALL, TMAX_NOD_ALL, TOLDMIN_NOD_ALL, TOLDMAX_NOD_ALL, &
                T2MIN_NOD_ALL, T2MAX_NOD_ALL, T2OLDMIN_NOD_ALL, T2OLDMAX_NOD_ALL, &
                DENMIN_NOD_ALL, DENMAX_NOD_ALL, DENOLDMIN_NOD_ALL, DENOLDMAX_NOD_ALL )
                ! For each node, find the largest and smallest value of T and
                ! DENSITY for both the current and previous timestep, out of
                ! the node value and all its surrounding nodes including Dirichlet b.c's.
                IMPLICIT NONE
                INTEGER, intent( in ) :: NPHASE,CV_NONODS, STOTEL,CV_SNLOC, &
                    IGOT_T2
                INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
                REAL, DIMENSION( :, :, : ), intent( in ) :: SUF_T_BC_ALL, SUF_D_BC_ALL
                REAL, DIMENSION( :, :, : ), intent( in ), pointer :: SUF_T2_BC_ALL
                INTEGER, DIMENSION( : , : , : ), intent( in ) :: WIC_T_BC_ALL, WIC_D_BC_ALL
                INTEGER, DIMENSION( : , : , : ), intent( in ) :: WIC_T2_BC_ALL
                INTEGER, DIMENSION( : ), intent( in ) :: FINACV
                INTEGER, DIMENSION( : ), intent( in ), target :: COLACV
                REAL, DIMENSION( :, : ), intent( inout ) :: TMAX_ALL, TMIN_ALL, TOLDMAX_ALL, TOLDMIN_ALL,  &
                    DENMAX_ALL, DENMIN_ALL, DENOLDMAX_ALL, DENOLDMIN_ALL
                REAL, DIMENSION( :, : ), intent( inout ) :: T2MAX_ALL, T2MIN_ALL, T2OLDMAX_ALL, T2OLDMIN_ALL

                REAL, DIMENSION( :, : ), intent( in ) :: T_ALL,TOLD_ALL,DEN_ALL,DENOLD_ALL
                REAL, DIMENSION( :, : ), intent( in ) :: T2_ALL,T2OLD_ALL
                INTEGER, DIMENSION( :, : ), intent( inout ) :: TMIN_NOD_ALL, TMAX_NOD_ALL, TOLDMIN_NOD_ALL, &
                    TOLDMAX_NOD_ALL, DENMIN_NOD_ALL, DENMAX_NOD_ALL, DENOLDMIN_NOD_ALL, DENOLDMAX_NOD_ALL
                INTEGER, DIMENSION( :, : ), intent( inout ) :: T2MIN_NOD_ALL, T2MAX_NOD_ALL, T2OLDMIN_NOD_ALL, &
                    T2OLDMAX_NOD_ALL
                ! Local variables
                INTEGER :: CV_NODI, IPHASE, CV_SILOC, SELE, CV_INOD
                integer, dimension(:), pointer :: cv_neigh_ptr

                Loop_CV_NODI: DO CV_NODI = 1, CV_NONODS

                    cv_neigh_ptr=>colacv(finacv(cv_nodi):finacv(cv_nodi+1)-1)

                    DO IPHASE = 1, NPHASE
                        TMAX_ALL( IPHASE, CV_NODI ) = maxval(T_ALL( IPHASE, cv_neigh_ptr ))
                        TMAX_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(maxloc(T_ALL( IPHASE, cv_neigh_ptr )))  ! COLN OF THE MAXIMUM VALUE
                        TMIN_ALL( IPHASE, CV_NODI ) = minval(T_ALL( IPHASE, cv_neigh_ptr ))
                        TMIN_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(minloc(T_ALL( IPHASE, cv_neigh_ptr )))
                        TOLDMAX_ALL( IPHASE, CV_NODI ) = maxval(TOLD_ALL( IPHASE, cv_neigh_ptr ))
                        TOLDMAX_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(maxloc(TOLD_ALL( IPHASE, cv_neigh_ptr )))
                        TOLDMIN_ALL( IPHASE, CV_NODI ) = minval(TOLD_ALL( IPHASE, cv_neigh_ptr ))
                        TOLDMIN_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(minloc(TOLD_ALL( IPHASE, cv_neigh_ptr )))
                        IF(IGOT_T2==1) THEN
                            T2MAX_ALL( IPHASE, CV_NODI ) = maxval(T2_ALL( IPHASE, cv_neigh_ptr ))
                            T2MAX_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) =cv_neigh_ptr(maxloc(T2_ALL( IPHASE, cv_neigh_ptr )))
                            T2MIN_ALL( IPHASE, CV_NODI ) = minval(T2_ALL( IPHASE, cv_neigh_ptr ))
                            T2MIN_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(minloc(T2_ALL( IPHASE, cv_neigh_ptr )))
                            T2OLDMAX_ALL( IPHASE, CV_NODI ) = maxval(T2OLD_ALL( IPHASE, cv_neigh_ptr ))
                            T2OLDMAX_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(maxloc(T2OLD_ALL( IPHASE, cv_neigh_ptr )))
                            T2OLDMIN_ALL( IPHASE, CV_NODI ) = minval(T2OLD_ALL( IPHASE, cv_neigh_ptr ))
                            T2OLDMIN_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(minloc(T2OLD_ALL( IPHASE, cv_neigh_ptr )))
                        ENDIF
                        DENMAX_ALL( IPHASE, CV_NODI ) = maxval(DEN_ALL( IPHASE, cv_neigh_ptr ))
                        DENMAX_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(maxloc(DEN_ALL( IPHASE, cv_neigh_ptr )))
                        DENMIN_ALL( IPHASE, CV_NODI ) = minval(DEN_ALL( IPHASE, cv_neigh_ptr ))
                        DENMIN_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(minloc(DEN_ALL( IPHASE, cv_neigh_ptr )))
                        DENOLDMAX_ALL( IPHASE, CV_NODI ) = maxval(DENOLD_ALL( IPHASE, cv_neigh_ptr ))
                        DENOLDMAX_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(maxloc(DENOLD_ALL( IPHASE, cv_neigh_ptr )))
                        DENOLDMIN_ALL( IPHASE, CV_NODI ) = minval(DENOLD_ALL( IPHASE, cv_neigh_ptr ))
                        DENOLDMIN_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(minloc(DENOLD_ALL( IPHASE, cv_neigh_ptr )))
                    END DO

                END DO Loop_CV_NODI

                ! Take into account the Dirichlet b.c's when working out max and min values.
                Loop_SELE: DO SELE= 1, STOTEL

                    Loop_CV_SILOC: DO CV_SILOC = 1, CV_SNLOC

                        CV_INOD=CV_SNDGLN((SELE-1)*CV_SNLOC+CV_SILOC)
                        !          SUF_CV_SI=(SELE-1)*CV_SNLOC+CV_SILOC

                        DO IPHASE=1,NPHASE
                            !             SUF_CV_SI_IPHA = SUF_CV_SI + STOTEL * CV_SNLOC * ( IPHASE - 1 )
                            !             CV_INOD_IPHA=CV_INOD + CV_NONODS*(IPHASE-1)
                            IF( (WIC_T_BC_ALL(1 ,IPHASE, SELE) == WIC_T_BC_DIRICHLET) &
                                .OR.(WIC_T_BC_ALL(1, IPHASE, SELE) == WIC_T_BC_DIRI_ADV_AND_ROBIN)) THEN
                                IF(SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC*(SELE-1) ) > TMAX_ALL( IPHASE, CV_INOD ) ) THEN
                                    TMAX_ALL( IPHASE, CV_INOD ) = SUF_T_BC_ALL( 1,  IPHASE, CV_SILOC+CV_SNLOC*(SELE-1) )
                                    TMAX_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                                ENDIF
                                IF(SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC*(SELE-1) ) < TMIN_ALL( IPHASE, CV_INOD ) ) THEN
                                    TMIN_ALL( IPHASE, CV_INOD ) = SUF_T_BC_ALL( 1, IPHASE, CV_SILOC+ CV_SNLOC*(SELE-1) )
                                    TMIN_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                                ENDIF

                                IF(SUF_T_BC_ALL( 1, IPHASE, CV_SILOC+ CV_SNLOC*(SELE-1) ) > TOLDMAX_ALL( IPHASE, CV_INOD ) ) THEN
                                    TOLDMAX_ALL( IPHASE, CV_INOD ) = SUF_T_BC_ALL(1, IPHASE, CV_SILOC + CV_SNLOC*( SELE-1 ) )
                                    TOLDMAX_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                                ENDIF
                                IF(SUF_T_BC_ALL( 1 , IPHASE, CV_SILOC + CV_SNLOC* ( SELE -1 ) ) < TOLDMIN_ALL( IPHASE, CV_INOD ) ) THEN
                                    TOLDMIN_ALL( IPHASE, CV_INOD ) = SUF_T_BC_ALL( 1 , IPHASE, CV_SILOC + CV_SNLOC* ( SELE -1) )
                                    TOLDMIN_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                                ENDIF
                            ENDIF

                            ! T2:
                            IF(IGOT_T2==1) THEN
                                IF( (WIC_T2_BC_ALL(1, IPHASE, SELE) == WIC_T_BC_DIRICHLET) &
                                    .OR.(WIC_T2_BC_ALL(1, IPHASE, SELE) == WIC_T_BC_DIRI_ADV_AND_ROBIN)) THEN
                                    IF(SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC* ( SELE-1 ) ) > T2MAX_ALL( IPHASE, CV_INOD ) ) THEN
                                        T2MAX_ALL( IPHASE, CV_INOD ) = SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) )
                                        T2MAX_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                                    ENDIF
                                    IF(SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1) ) < T2MIN_ALL( IPHASE, CV_INOD ) ) THEN
                                        T2MIN_ALL( IPHASE, CV_INOD ) = SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOc * ( SELE - 1 ) )
                                        T2MIN_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                                    ENDIF

                                    IF(SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC* ( SELE - 1 ) ) > T2OLDMAX_ALL( IPHASE, CV_INOD ) ) THEN
                                        T2OLDMAX_ALL( IPHASE, CV_INOD ) = SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) )
                                        T2OLDMAX_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                                    ENDIF
                                    IF(SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) ) < T2OLDMIN_ALL( IPHASE, CV_INOD ) ) THEN
                                        T2OLDMIN_ALL( IPHASE, CV_INOD ) = SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) )
                                        T2OLDMIN_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                                    ENDIF
                                ENDIF
                            ENDIF
                            ! DEN:
                            IF( WIC_D_BC_ALL(1 , IPHASE, SELE) == WIC_D_BC_DIRICHLET ) THEN
                                IF(SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) ) > DENMAX_ALL( IPHASE, CV_INOD ) ) THEN
                                    DENMAX_ALL( IPHASE, CV_INOD ) = SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) )
                                    DENMAX_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                                ENDIF
                                IF(SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOc * ( SELE - 1 ) ) < DENMIN_ALL( IPHASE, CV_INOD ) ) THEN
                                    DENMIN_ALL( IPHASE, CV_INOD ) = SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) )
                                    DENMIN_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                                ENDIF

                                IF(SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 )  ) > DENOLDMAX_ALL( IPHASE, CV_INOD ) ) THEN
                                    DENOLDMAX_ALL( IPHASE, CV_INOD ) = SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) )
                                    DENOLDMAX_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                                ENDIF
                                IF(SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) ) < DENOLDMIN_ALL( IPHASE, CV_INOD ) ) THEN
                                    DENOLDMIN_ALL( IPHASE, CV_INOD )= SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) )
                                    DENOLDMIN_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                                ENDIF
                            ENDIF
                        END DO

                    END DO Loop_CV_SILOC

                END DO Loop_SELE

                RETURN
            END SUBROUTINE SURRO_CV_MINMAX


            SUBROUTINE CALC_LIMIT_MATRIX_MAX_MIN(TMAX_ALL, TMIN_ALL, DENMAX_ALL, DENMIN_ALL, &
                T2MAX_ALL, T2MIN_ALL, &
                T_ALL,  T2_ALL, DEN_ALL, IGOT_T2, NPHASE, CV_NONODS, &
                TMIN_NOD_ALL, TMAX_NOD_ALL,  &
                T2MIN_NOD_ALL, T2MAX_NOD_ALL, &
                DENMIN_NOD_ALL, DENMAX_NOD_ALL, &
                NSMALL_COLM, SMALL_FINDRM, SMALL_COLM, &
                TUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, T2UPWIND_MAT_ALL, MASS_CV)
                ! Populate  limiting matrix based on max and min values
                ! For each node, find the largest and smallest value of T and
                ! DENSITY for both the current and previous timestep, out of
                ! the node value and all its surrounding nodes including Dirichlet b.c's.
                IMPLICIT NONE
                INTEGER, intent( in ) :: NPHASE, CV_NONODS, NSMALL_COLM, IGOT_T2
                INTEGER, DIMENSION( : ), intent( in ) :: SMALL_FINDRM
                INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
                REAL, DIMENSION( :, : ), intent( inout ) :: TMAX_ALL, TMIN_ALL, DENMAX_ALL, DENMIN_ALL
                REAL, DIMENSION( :, : ), intent( inout ) :: T2MAX_ALL, T2MIN_ALL
                REAL, DIMENSION( :, : ), intent( inout ) :: TUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, T2UPWIND_MAT_ALL

                REAL, DIMENSION( :, : ), intent( in ) :: T_ALL, DEN_ALL
                REAL, DIMENSION( :, :), intent( in ) :: T2_ALL
                REAL, DIMENSION( : ), intent( in ) :: MASS_CV
                INTEGER, DIMENSION( :, : ), intent( inout ) :: TMIN_NOD_ALL, TMAX_NOD_ALL, DENMIN_NOD_ALL, DENMAX_NOD_ALL
                INTEGER, DIMENSION( :, : ), intent( inout ) :: T2MIN_NOD_ALL, T2MAX_NOD_ALL
                ! Local variables
                INTEGER :: CV_NODI, CV_NODJ, IPHASE, COUNT
                INTEGER :: COUNT_OUT
                LOGICAL, PARAMETER :: LIM_VOL_ADJUST=.TRUE.
                REAL :: TMIN_STORE(NPHASE),TMAX_STORE(NPHASE),DENMIN_STORE(NPHASE),DENMAX_STORE(NPHASE)
                REAL :: T2MIN_STORE(NPHASE),T2MAX_STORE(NPHASE)
                LOGICAL :: RESET_STORE
                REAL :: INCOME


                TUPWIND_MAT_ALL(1:NPHASE, 1: NSMALL_COLM)  =0.0
                DENUPWIND_MAT_ALL(1:NPHASE, 1: NSMALL_COLM)=0.0
                IF(IGOT_T2==1) T2UPWIND_MAT_ALL(1:NPHASE, 1: NSMALL_COLM) =0.0


                DO CV_NODI=1,CV_NONODS
                    DO COUNT=SMALL_FINDRM(CV_NODI), SMALL_FINDRM(CV_NODI+1)-1

                        CV_NODJ=SMALL_COLM(COUNT)

                        ! for outgoing information to CV_NODI ...

                        INCOME=0.0

                        DO IPHASE=1,NPHASE
                            !       CV_NODI_IPHA = CV_NODI + (IPHASE-1)*CV_NONODS
                            !       CV_NODJ_IPHA = CV_NODJ + (IPHASE-1)*CV_NONODS
                            IF ( LIM_VOL_ADJUST ) THEN
                                RESET_STORE = .FALSE.
                                CALL CAL_LIM_VOL_ADJUST( TMIN_STORE(IPHASE), TMIN_ALL(IPHASE,:), T_ALL(IPHASE,:), TMIN_NOD_ALL(IPHASE,:), RESET_STORE, MASS_CV, &
                                    CV_NODI, CV_NODJ, 1, CV_NONODS, INCOME )
                                CALL CAL_LIM_VOL_ADJUST( TMAX_STORE(IPHASE), TMAX_ALL(IPHASE,:), T_ALL(IPHASE,:), TMAX_NOD_ALL(IPHASE,:), RESET_STORE, MASS_CV, &
                                    CV_NODI, CV_NODJ, 1, CV_NONODS, INCOME )
                                CALL CAL_LIM_VOL_ADJUST( DENMIN_STORE(IPHASE), DENMIN_ALL(IPHASE,:), DEN_ALL(IPHASE,:), DENMIN_NOD_ALL(IPHASE,:), RESET_STORE, MASS_CV, &
                                    CV_NODI, CV_NODJ, 1, CV_NONODS, INCOME )
                                CALL CAL_LIM_VOL_ADJUST( DENMAX_STORE(IPHASE), DENMAX_ALL(IPHASE,:), DEN_ALL(IPHASE,:), DENMAX_NOD_ALL(IPHASE,:), RESET_STORE, MASS_CV, &
                                    CV_NODI, CV_NODJ, 1, CV_NONODS, INCOME )

                                IF(IGOT_T2==1) THEN
                                    CALL CAL_LIM_VOL_ADJUST( T2MIN_STORE(IPHASE), T2MIN_ALL(IPHASE,:), T2_ALL(IPHASE,:), T2MIN_NOD_ALL(IPHASE,:), RESET_STORE, MASS_CV, &
                                        CV_NODI, CV_NODJ, 1, CV_NONODS, INCOME )
                                    CALL CAL_LIM_VOL_ADJUST( T2MAX_STORE(IPHASE), T2MAX_ALL(IPHASE,:), T2_ALL(IPHASE,:), T2MAX_NOD_ALL(IPHASE,:), RESET_STORE, MASS_CV, &
                                        CV_NODI, CV_NODJ, 1, CV_NONODS, INCOME )
                                END IF

                            END IF
                        END DO



                        ! ***PUT INTO MATRIX**************
                        ! Populate  limiting matrix based on max and min values


                        COUNT_OUT= COUNT

                        DO IPHASE=1,NPHASE

                            IF(T_ALL( IPHASE, CV_NODI ).GT.T_ALL( IPHASE, CV_NODJ )) THEN
                                TUPWIND_MAT_ALL( IPHASE, COUNT_OUT ) = TMAX_ALL( IPHASE, CV_NODI )
                                DENUPWIND_MAT_ALL( IPHASE, COUNT_OUT ) = DENMAX_ALL( IPHASE, CV_NODI )
                                IF(IGOT_T2==1) THEN
                                    T2UPWIND_MAT_ALL( IPHASE, COUNT_OUT ) = T2MAX_ALL( IPHASE, CV_NODI )
                                ENDIF
                            ELSE
                                TUPWIND_MAT_ALL( IPHASE, COUNT_OUT ) = TMIN_ALL( IPHASE, CV_NODI )
                                DENUPWIND_MAT_ALL( IPHASE, COUNT_OUT ) = DENMIN_ALL( IPHASE, CV_NODI )
                                IF(IGOT_T2==1) THEN
                                    T2UPWIND_MAT_ALL( IPHASE, COUNT_OUT ) = T2MIN_ALL( IPHASE, CV_NODI )
                                ENDIF
                            ENDIF
                        END DO




                        DO IPHASE=1,NPHASE
                            IF ( LIM_VOL_ADJUST ) THEN
                                RESET_STORE = .TRUE.
                                CALL CAL_LIM_VOL_ADJUST(TMIN_STORE(IPHASE),TMIN_ALL(IPHASE,:),T_ALL(IPHASE,:),TMIN_NOD_ALL(IPHASE,:),RESET_STORE,MASS_CV, &
                                    CV_NODI, CV_NODJ, 1, CV_NONODS, INCOME )
                                CALL CAL_LIM_VOL_ADJUST(TMAX_STORE(IPHASE),TMAX_ALL(IPHASE,:),T_ALL(IPHASE,:),TMAX_NOD_ALL(IPHASE,:),RESET_STORE,MASS_CV, &
                                    CV_NODI, CV_NODJ, 1, CV_NONODS, INCOME )

                                CALL CAL_LIM_VOL_ADJUST(DENMIN_STORE(IPHASE),DENMIN_ALL(IPHASE,:),DEN_ALL(IPHASE,:),DENMIN_NOD_ALL(IPHASE,:),RESET_STORE,MASS_CV, &
                                    CV_NODI, CV_NODJ, 1, CV_NONODS, INCOME )
                                CALL CAL_LIM_VOL_ADJUST(DENMAX_STORE(IPHASE),DENMAX_ALL(IPHASE,:),DEN_ALL(IPHASE,:),DENMAX_NOD_ALL(IPHASE,:),RESET_STORE,MASS_CV, &
                                    CV_NODI, CV_NODJ, 1, CV_NONODS, INCOME )

                                IF ( IGOT_T2 == 1 ) THEN
                                    CALL CAL_LIM_VOL_ADJUST(T2MIN_STORE(IPHASE),T2MIN_ALL(IPHASE,:),T2_ALL(IPHASE,:),T2MIN_NOD_ALL(IPHASE,:),RESET_STORE,MASS_CV, &
                                        CV_NODI, CV_NODJ, 1, CV_NONODS, INCOME )
                                    CALL CAL_LIM_VOL_ADJUST(T2MAX_STORE(IPHASE),T2MAX_ALL(IPHASE,:),T2_ALL(IPHASE,:),T2MAX_NOD_ALL(IPHASE,:),RESET_STORE,MASS_CV, &
                                        CV_NODI, CV_NODJ, 1, CV_NONODS, INCOME )
                                END IF

                            END IF
                        END DO



                    END DO
                END DO



                RETURN

            END SUBROUTINE CALC_LIMIT_MATRIX_MAX_MIN



    END SUBROUTINE ISOTROPIC_LIMITER_ALL


    SUBROUTINE CALC_SELE( ELE, ELE3, SELE, CV_SILOC, CV_ILOC, U_SLOC2LOC, CV_SLOC2LOC, &
        FACE_ELE, gi, CV_funs, Mdims, CV_GIdims, &
        CV_NDGLN, U_NDGLN, CV_SNDGLN, U_SNDGLN )
        ! Calculate SELE, CV_SILOC, U_SLOC2LOC, CV_SLOC2LOC for a face on the
        ! boundary of the domain
        IMPLICIT NONE
        INTEGER, intent( in ) :: ELE, CV_ILOC, gi
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_shape_funs), intent(in) :: CV_funs
        type(multi_GI_dimensions), intent(in) :: CV_GIdims
        INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN ! (Mdims%cv_nloc*Mdims%totele)
        INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN  ! (Mdims%u_nloc*Mdims%totele)
        INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN ! (Mdims%cv_snloc*Mdims%totele)
        INTEGER, DIMENSION( : ), intent( in ) :: U_SNDGLN  ! (Mdims%u_snloc*Mdims%totele)
        INTEGER, DIMENSION( :, : ), intent( in ) :: FACE_ELE !(Mdims%totele*IFACE)
        INTEGER, intent( inout ) :: SELE, ELE3, CV_SILOC
        INTEGER, DIMENSION( : ), intent( inout ) :: U_SLOC2LOC !(CV_SKLOC)
        INTEGER, DIMENSION( : ), intent( inout ) :: CV_SLOC2LOC
        ! local variables
        INTEGER :: IFACE, ELE2, SELE2, CV_JLOC, CV_JNOD, &
            U_JLOC, U_JNOD, CV_SKNOD, &
            U_SKLOC, U_SKNOD, CV_SKLOC, CV_SKLOC2, I
        LOGICAL :: FOUND
        INTEGER, DIMENSION( Mdims%cv_snloc ) :: LOG_ON_BOUND
        log_on_bound= -66666
        !ewrite(3,*)'In Calc_Sele'
        I = 1
        DO CV_JLOC = 1, Mdims%cv_nloc
            CV_JNOD = CV_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_JLOC )
            IF ( .NOT. cv_funs%cvfem_on_face( CV_JLOC, gi ) ) THEN
                LOG_ON_BOUND( I ) = CV_JNOD
                I = I + 1
            END IF
        END DO
        SELE = 0
        ELE3 = 0
        ! What face are we on
        DO IFACE = 1, CV_GIdims%nface
            ELE2 = FACE_ELE( IFACE, ELE )
            SELE2 = MAX( 0, - ELE2 )
            ELE2 = MAX( 0, + ELE2 )
            IF ( SELE2 /= 0 ) THEN
                FOUND = .TRUE.
                DO CV_SKLOC = 1, Mdims%cv_snloc
                    CV_SKNOD = CV_SNDGLN( ( SELE2 - 1 ) * Mdims%cv_snloc + CV_SKLOC )
                    DO CV_SKLOC2 = 1, Mdims%cv_snloc
                        IF ( CV_SKNOD == LOG_ON_BOUND( CV_SKLOC2 ) ) THEN
                            FOUND = .FALSE.
                            EXIT
                        ENDIF
                    END DO
                    IF(.NOT.FOUND) EXIT
                END DO
                IF( FOUND ) THEN
                    SELE = SELE2
                    ELE3 = ELE2
                    exit
                ENDIF
            END IF
        END DO

        Conditional_Sele: IF ( SELE /= 0 ) THEN
            DO CV_SKLOC = 1, Mdims%cv_snloc
                CV_SKNOD = CV_SNDGLN( ( SELE - 1 ) * Mdims%cv_snloc + CV_SKLOC )
                DO CV_JLOC = 1, Mdims%cv_nloc
                    CV_JNOD = CV_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_JLOC )
                    IF( CV_SKNOD == CV_JNOD ) EXIT
                END DO
                CV_SLOC2LOC( CV_SKLOC ) = CV_JLOC
                IF( CV_JLOC == CV_ILOC ) CV_SILOC = CV_SKLOC
            END DO

            ! Calculate U_SLOC2LOC
            DO U_SKLOC = 1, Mdims%u_snloc
                U_SKNOD = U_SNDGLN( ( SELE - 1 ) * Mdims%u_snloc + U_SKLOC )
                DO U_JLOC = 1, Mdims%u_nloc
                    U_JNOD = U_NDGLN( ( ELE - 1 ) * Mdims%u_nloc + U_JLOC )
                    IF( U_SKNOD == U_JNOD ) EXIT
                END DO
                U_SLOC2LOC( U_SKLOC ) = U_JLOC
            END DO
        END IF Conditional_Sele
        RETURN
    END SUBROUTINE CALC_SELE




    SUBROUTINE PUT_IN_CT_RHS( GET_C_IN_CV_ADVDIF_AND_CALC_C_CV, ct_rhs_phase_cv_nodi, ct_rhs_phase_cv_nodj, &
        Mdims, CV_funs, ndgln, Mmat, GI, between_elements, on_domain_boundary, &
        ELE, ELE2, SELE, HDC, MASS_ELE, JCOUNT_KLOC, JCOUNT_KLOC2, ICOUNT_KLOC, ICOUNT_KLOC2, &
        C_JCOUNT_KLOC, C_JCOUNT_KLOC2, C_ICOUNT_KLOC, C_ICOUNT_KLOC2, U_OTHER_LOC, &
        U_SLOC2LOC, CV_SLOC2LOC,  &
        SCVDETWEI, CVNORMX_ALL, DEN_ALL_DIVID, CV_NODI, CV_NODJ, &
        WIC_U_BC_ALL, WIC_P_BC_ALL,SUF_P_BC_ALL,&
        UGI_COEF_ELE_ALL,  &
        UGI_COEF_ELE2_ALL,  &
        NDOTQ, NDOTQOLD, LIMT, LIMDT, LIMDTOLD, LIMT_HAT, &
        NDOTQ_HAT, &
        FTHETA_T2, ONE_M_FTHETA_T2OLD, FTHETA_T2_J, ONE_M_FTHETA_T2OLD_J, integrate_other_side_and_not_boundary, &
        theta_cty_solid, &
        loc_u, THETA_VEL,&
        ! local memory sent down for speed...
        UDGI_IMP_ALL, RCON, RCON_J, NDOTQ_IMP)
        ! This subroutine caculates the discretised cty eqn acting on the velocities i.e. Mmat%CT, Mmat%CT_RHS
        IMPLICIT NONE
        ! IF more_in_ct THEN PUT AS MUCH AS POSSIBLE INTO Mmat%CT MATRIX
        !    LOGICAL, PARAMETER :: more_in_ct=.false.
        INTEGER, intent( in ) :: GI, &
            CV_NODI, CV_NODJ, ELE, ELE2, SELE
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_shape_funs), intent(in) :: CV_funs
        type(multi_ndgln), intent(in) :: ndgln
        type (multi_matrices), intent(inout) :: Mmat
        REAL, DIMENSION( :, :, : ), intent( in ) :: loc_u
        LOGICAL, intent( in ) :: integrate_other_side_and_not_boundary, between_elements, on_domain_boundary,&
            GET_C_IN_CV_ADVDIF_AND_CALC_C_CV
        INTEGER, DIMENSION( : ), intent( in ) :: JCOUNT_KLOC, JCOUNT_KLOC2, ICOUNT_KLOC, ICOUNT_KLOC2, U_OTHER_LOC
        INTEGER, DIMENSION( : ), intent( in ) :: C_JCOUNT_KLOC, C_JCOUNT_KLOC2, C_ICOUNT_KLOC, C_ICOUNT_KLOC2
        INTEGER, DIMENSION( : ), intent( in ) :: U_SLOC2LOC, CV_SLOC2LOC
        REAL, DIMENSION( :, :, : ), intent( inout ) :: SUF_P_BC_ALL
        REAL, DIMENSION( : ), intent( inout ) :: ct_rhs_phase_cv_nodi, ct_rhs_phase_cv_nodj
        REAL, DIMENSION( Mdims%ndim, Mdims%nphase, Mdims%u_nloc ), intent( in ) :: UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL
        REAL, DIMENSION( : ), intent( in ) :: SCVDETWEI, MASS_ELE
        REAL, DIMENSION( :, : ), intent( in ) :: CVNORMX_ALL
        REAL, DIMENSION( :, : ), intent( in ) :: DEN_ALL_DIVID
        REAL, DIMENSION( : ), intent( in ) :: NDOTQ, NDOTQOLD, LIMT, LIMDT, LIMDTOLD, LIMT_HAT
        REAL, intent( in ) :: NDOTQ_HAT
        REAL, DIMENSION( : ), intent( in ) :: THETA_VEL
        integer, dimension(:,:,:) :: WIC_U_BC_ALL, WIC_P_BC_ALL
        ! LIMT_HAT is the normalised voln fraction
        REAL, intent( in ) :: theta_cty_solid, HDC
        REAL,  DIMENSION( Mdims%nphase ), intent( in ) :: FTHETA_T2, ONE_M_FTHETA_T2OLD, FTHETA_T2_J, ONE_M_FTHETA_T2OLD_J
        ! local memory sent down for speed...
        REAL,  DIMENSION( Mdims%ndim, Mdims%nphase ), intent( inout ) :: UDGI_IMP_ALL
        REAL,  DIMENSION( Mdims%nphase ), intent( inout ) :: RCON, RCON_J, NDOTQ_IMP
        !Variable to account for boundary conditions if using GET_C_IN_CV_ADVDIF_AND_CALC_C_CV
        real, dimension (Mdims%ndim, Mdims%nphase, Mdims%u_nloc ) :: Bound_ele_correct
        ! Local variables...
        INTEGER :: U_KLOC, U_KLOC2, IDIM, &
            IPHASE, U_SKLOC, I, J, U_KKLOC, &
            u_iloc, u_siloc, count, COUNT_SUF
        real :: Mass_corrector
        !Local variables for CV pressure bcs
        integer :: KPHASE, CV_SNODK, CV_SNODK_IPHA, CV_SKLOC
        real, dimension(:,:), allocatable :: SUF_SIG_DIAGTEN_BC_pha_GI

        !If using Mmat%C_CV prepare Bound_ele_correct and Bound_ele2_correct to correctly apply the BCs
        if (Mmat%CV_pressure) call introduce_C_CV_boundary_conditions(Bound_ele_correct)

        DO U_KLOC = 1, Mdims%u_nloc
            RCON(:) = SCVDETWEI( GI ) * (  FTHETA_T2(:) * LIMDT(:) + ONE_M_FTHETA_T2OLD(:) * LIMDTOLD(:) * THETA_VEL(:)) &
                * CV_funs%sufen( U_KLOC, GI ) / DEN_ALL_DIVID( :, CV_NODI )
            DO IPHASE = 1, Mdims%nphase
                Mmat%CT( :, IPHASE, JCOUNT_KLOC( U_KLOC ) ) = Mmat%CT( :, IPHASE, JCOUNT_KLOC( U_KLOC ) ) &
                    + rcon(IPHASE) * UGI_COEF_ELE_ALL( :, IPHASE, U_KLOC ) * CVNORMX_ALL( :, GI )
            END DO
            IF(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV) THEN
                rcon(:) = SCVDETWEI( GI ) * CV_funs%sufen( U_KLOC, GI )
                DO IPHASE=1,Mdims%n_in_pres!Mdims%nphase
                    IF ( between_elements) THEN
                        ! bias the weighting towards bigger eles - works with 0.25 and 0.1 and not 0.01.
                        !This is to perform the average between two DG pressures (same mass => 0.5)
                        Mass_corrector = (MASS_ELE( ELE2 ) + 0.25 * MASS_ELE( ELE ))/(1.25*(MASS_ELE( ELE ) + MASS_ELE( ELE2 )))

                        Mmat%C_CV( :, IPHASE, C_JCOUNT_KLOC( U_KLOC ) ) &
                            = Mmat%C_CV( :, IPHASE, C_JCOUNT_KLOC( U_KLOC ) ) &
                            + rcon(IPHASE) * CVNORMX_ALL( :, GI ) * Mass_corrector
                    else
                        Mmat%C_CV( :, IPHASE, C_JCOUNT_KLOC( U_KLOC ) ) &
                            = Mmat%C_CV( :, IPHASE, C_JCOUNT_KLOC( U_KLOC ) ) &
                            + rcon(IPHASE) * CVNORMX_ALL( :, GI ) * Bound_ele_correct(:, IPHASE, U_KLOC)
                    endif
                END DO
            ENDIF
            ! flux from the other side (change of sign because normal is -ve)...
            if ( integrate_other_side_and_not_boundary ) then
                RCON_J(:) = SCVDETWEI( GI ) * ( FTHETA_T2_J(:)* LIMDT(:) + ONE_M_FTHETA_T2OLD_J(:) * LIMDTOLD(:) * THETA_VEL(:))  &
                    * CV_funs%sufen( U_KLOC, GI ) / DEN_ALL_DIVID( :, CV_NODJ )
                DO IPHASE = 1, Mdims%nphase
                    Mmat%CT( :, IPHASE, ICOUNT_KLOC( U_KLOC ) ) = Mmat%CT( :, IPHASE, ICOUNT_KLOC( U_KLOC ) ) &
                        - RCON_J(IPHASE) * UGI_COEF_ELE_ALL( :, IPHASE, U_KLOC ) * CVNORMX_ALL( :, GI )
                END DO
                IF(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV) THEN
                    RCON_J(:) = SCVDETWEI( GI ) * CV_funs%sufen( U_KLOC, GI )
                    DO IPHASE=1,Mdims%n_in_pres!Mdims%nphase
                        IF ( between_elements ) THEN
                            Mmat%C_CV( :, IPHASE, C_ICOUNT_KLOC( U_KLOC ) ) &
                                = Mmat%C_CV( :, IPHASE, C_ICOUNT_KLOC( U_KLOC ) ) &
                                - RCON_J(IPHASE) * CVNORMX_ALL( :, GI )* Mass_corrector!(1.- Mass_corrector)
                        else
                            Mmat%C_CV( :, IPHASE, C_ICOUNT_KLOC( U_KLOC ) ) &
                                = Mmat%C_CV( :, IPHASE, C_ICOUNT_KLOC( U_KLOC ) ) &
                                - RCON_J(IPHASE) * CVNORMX_ALL( :, GI )* Bound_ele_correct(:, IPHASE, U_KLOC)!Bound_ele_correct unnecessary here?
                        endif
                    END DO
                ENDIF
            end if  ! endof if ( integrate_other_side_and_not_boundary ) then
        END DO
        IF ( on_domain_boundary ) THEN
            UDGI_IMP_ALL=0.0
            DO U_KLOC = 1, Mdims%u_nloc
                DO IPHASE = 1, Mdims%nphase
                    UDGI_IMP_ALL(:,IPHASE) = UDGI_IMP_ALL(:,IPHASE) + CV_funs%sufen( U_KLOC, GI ) * &
                        UGI_COEF_ELE_ALL( :, IPHASE, U_KLOC ) * LOC_U( :, IPHASE, U_KLOC )
                END DO
            END DO
            DO IPHASE = 1, Mdims%nphase
                NDOTQ_IMP(IPHASE)= SUM( CVNORMX_ALL( :,GI ) * UDGI_IMP_ALL(:,IPHASE) )
            END DO
            ct_rhs_phase_cv_nodi(:)=ct_rhs_phase_cv_nodi(:) &
                - SCVDETWEI( GI ) * (  ( &
                ONE_M_FTHETA_T2OLD(:) * LIMDTOLD(:) * (NDOTQOLD(:) -NDOTQ_IMP(:)*THETA_VEL(:)) &
                + FTHETA_T2(:)  * LIMDT(:) * (NDOTQ(:)-NDOTQ_IMP(:)) &
                ) / DEN_ALL_DIVID( :, CV_NODI ) )
        ELSE
            ct_rhs_phase_cv_nodi(:)=ct_rhs_phase_cv_nodi(:) &
                - SCVDETWEI( GI ) * (  ( &
                ONE_M_FTHETA_T2OLD(:) * LIMDTOLD(:) * NDOTQOLD(:) * (1.-THETA_VEL(:)) &
                ) / DEN_ALL_DIVID( :, CV_NODI )   )
            ! flux from the other side (change of sign because normal is -ve)...
            if ( integrate_other_side_and_not_boundary ) then
                ct_rhs_phase_cv_nodj(:)=ct_rhs_phase_cv_nodj(:) &
                    + SCVDETWEI( GI ) * ( ( &
                    ONE_M_FTHETA_T2OLD_J(:) * LIMDTOLD(:) * NDOTQOLD(:) * (1.-THETA_VEL(:)) &
                    ) / DEN_ALL_DIVID( :, CV_NODJ ) )
            end if
        END IF
        IF ( between_elements ) THEN
            ! We have a discontinuity between elements so integrate along the face...
            DO U_SKLOC = 1, Mdims%u_snloc
                U_KLOC = U_SLOC2LOC(U_SKLOC)
                U_KLOC2 = U_OTHER_LOC( U_KLOC )
                RCON(:) = SCVDETWEI( GI ) * (  FTHETA_T2(:) * LIMDT(:) + ONE_M_FTHETA_T2OLD(:) * LIMDTOLD(:) * THETA_VEL(:)) &
                    * CV_funs%sufen( U_KLOC, GI ) / DEN_ALL_DIVID( :, CV_NODI  )
                DO IPHASE = 1, Mdims%nphase
                    Mmat%CT( :, IPHASE, JCOUNT_KLOC2( U_KLOC2 ) ) &
                        = Mmat%CT( :, IPHASE, JCOUNT_KLOC2( U_KLOC2 ) ) &
                        + rcon(IPHASE) * UGI_COEF_ELE2_ALL( :, IPHASE, U_KLOC2 ) * CVNORMX_ALL( :, GI )
                END DO
                IF(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV) THEN
                    RCON(:) = SCVDETWEI( GI ) * CV_funs%sufen( U_KLOC, GI )
                    DO IPHASE=1,Mdims%n_in_pres!Mdims%nphase
                        Mmat%C_CV( :, IPHASE, C_JCOUNT_KLOC2( U_KLOC2 ) ) &
                            = Mmat%C_CV( :, IPHASE, C_JCOUNT_KLOC2( U_KLOC2 ) ) &
                            + RCON(IPHASE) * CVNORMX_ALL( :, GI )* (1.- Mass_corrector)
                    END DO
                ENDIF
                ! flux from the other side (change of sign because normal is -ve)...
                if ( integrate_other_side_and_not_boundary ) then
                    RCON_J(:) = SCVDETWEI( GI ) * ( FTHETA_T2_J(:)* LIMDT(:) + ONE_M_FTHETA_T2OLD_J(:) * LIMDTOLD(:) * THETA_VEL(:)) &
                        * CV_funs%sufen( U_KLOC, GI ) / DEN_ALL_DIVID( :, CV_NODJ )
                    DO IPHASE=1,Mdims%nphase
                        Mmat%CT( :, IPHASE, ICOUNT_KLOC2( U_KLOC2 ) ) &
                            = Mmat%CT( :, IPHASE, ICOUNT_KLOC2( U_KLOC2 ) ) &
                            - RCON_J(IPHASE) * UGI_COEF_ELE2_ALL( :, IPHASE, U_KLOC2 ) * CVNORMX_ALL( :, GI )
                    END DO
                    IF(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV) THEN
                        RCON_J(:) = SCVDETWEI( GI ) * CV_funs%sufen( U_KLOC, GI )
                        DO IPHASE=1,Mdims%n_in_pres!Mdims%nphase
                            Mmat%C_CV( :, IPHASE, C_ICOUNT_KLOC2( U_KLOC2 ) ) &
                                = Mmat%C_CV( :, IPHASE, C_ICOUNT_KLOC2( U_KLOC2 ) ) &
                                - RCON_J(IPHASE) * CVNORMX_ALL( :, GI )* (1.-Mass_corrector)!Mass_corrector
                        END DO
                    ENDIF
                end if  ! endof if ( integrate_other_side_and_not_boundary ) then
            END DO
        END IF ! endof IF ( between_elements ) THEN
        RETURN
    contains

        subroutine introduce_C_CV_boundary_conditions(Bound_ele_correct)
            !This subroutine populates Bound_ele_correct and Bound_ele2_correct to properly apply the BCs when creating the
            !Mmat%C_CV matrix
            implicit none
            real, dimension(:,:,:), intent(out) :: Bound_ele_correct
            !Local variables
            integer :: U_KLOC, IPHASE, P_SJLOC, U_INOD, ipres, CV_KLOC, P_ILOC
            logical, save :: show_warn_msg = .true.
            !By default no modification is required
            Bound_ele_correct = 1.0
            IF ( on_domain_boundary ) THEN
                !By default the position must not be added to the matrix
                Bound_ele_correct = 0.!<= P in the CV == P in the BC, it is done this way

                !If Mmat%C_CV formulation, apply weak pressure boundary conditions if any
                DO IPRES = 1, 1!Mdims%npres
                    IF( WIC_P_BC_ALL( 1,IPRES,SELE ) == WIC_P_BC_DIRICHLET ) THEN
                        DO U_SILOC = 1, Mdims%u_snloc
                            U_ILOC = U_SLOC2LOC( U_SILOC )
                            U_INOD = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                            DO IPHASE =  1+(IPRES-1)*Mdims%n_in_pres, IPRES*Mdims%n_in_pres
                                !We give priority to velocity boundary conditions
                                if (WIC_U_BC_ALL( 1, IPHASE, SELE ) /= WIC_U_BC_DIRICHLET ) then
                                    !Only in the boundaries with a defined pressure it needs to be added into
                                    !the matrix and into the RHS
                                    Bound_ele_correct( :, IPHASE, U_ILOC ) = 1.
                                    Mmat%U_RHS( :, IPHASE, U_INOD ) = Mmat%U_RHS( :, IPHASE, U_INOD ) &
                                        - CVNORMX_ALL( :, GI )* CV_funs%sufen( U_ILOC, GI )*SCVDETWEI( GI )&
                                        * SUF_P_BC_ALL( 1,1,1 + Mdims%cv_snloc* ( SELE - 1 ) )
                                else
                                    if (show_warn_msg) then
                                        ewrite(0,*) "WARNING: One or more boundaries have velocity and pressure boundary conditions."
                                        show_warn_msg = .false.
                                    end if
                                end if
                            end do
                        end do
                    end if
                end do
            end if
        end subroutine introduce_C_CV_boundary_conditions
    END SUBROUTINE PUT_IN_CT_RHS


    SUBROUTINE CALC_ANISOTROP_LIM(&
        ! Caculate the upwind values stored in matrix form...
        T_ALL, TOLD_ALL, DEN_ALL, DENOLD_ALL, T2_ALL, T2OLD_ALL, &
        FEMT_ALL, FEMTOLD_ALL, FEMDEN_ALL, FEMDENOLD_ALL, FEMT2_ALL, FEMT2OLD_ALL, USE_FEMT, &
        TUPWIND_MAT_ALL, TOLDUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, DENOLDUPWIND_MAT_ALL, &
        T2UPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL, &
        IGOT_T2, NPHASE, CV_NONODS,CV_NLOC, TOTELE, CV_NDGLN, &
        SMALL_FINDRM, SMALL_CENTRM, SMALL_COLM,NSMALL_COLM, &
        X_NDGLN, X_NONODS, NDIM, &
        X_ALL, XC_CV_ALL, use_reflect)
        ! For the anisotropic limiting scheme we find the upwind values
        ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
        ! value for each node pair is stored in the matrices TUPWIND AND
        IMPLICIT NONE
        INTEGER, intent( in ) :: CV_NONODS,X_NONODS,TOTELE,CV_NLOC, &
            NSMALL_COLM, NDIM,IGOT_T2,NPHASE
        REAL, DIMENSION( :, : ), intent( in ) :: T_ALL,TOLD_ALL,DEN_ALL,DENOLD_ALL
        REAL, DIMENSION( :,:), intent( in ), pointer :: T2_ALL,T2OLD_ALL
        REAL, DIMENSION( :, :), intent( in ) :: FEMT_ALL,FEMTOLD_ALL,FEMDEN_ALL,FEMDENOLD_ALL
        REAL, DIMENSION( :, :), intent( in ) :: FEMT2_ALL,FEMT2OLD_ALL
        LOGICAL, intent( in ) :: USE_FEMT, use_reflect ! Use the FEM solns rather than CV's when interpolating soln
        REAL, DIMENSION( :, : ), intent( inout ) :: TUPWIND_MAT_ALL, TOLDUPWIND_MAT_ALL, &
            DENUPWIND_MAT_ALL, DENOLDUPWIND_MAT_ALL
        REAL, DIMENSION( :, : ), intent( inout ) :: T2UPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL
        REAL, DIMENSION( :,: ), intent( in ) :: X_ALL
        INTEGER, DIMENSION(: ), intent( in ) :: X_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
        INTEGER, DIMENSION( : ), intent( in) :: SMALL_FINDRM
        INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
        INTEGER, DIMENSION( : ), intent( in) :: SMALL_CENTRM
        REAL, DIMENSION( NDIM, CV_NONODS), intent( in ) :: XC_CV_ALL
        ! local variables...
        INTEGER :: NFIELD, IMID, NOD, IFIELD
        REAL, DIMENSION( :, : ), ALLOCATABLE :: F_ALL, FEMF_ALL, FUPWIND_MAT_ALL


        ! Pack all the variables in:
        ! Always pack in the T & TOLD & T2, T2OLD variables
        NFIELD=(4 + 2*IGOT_T2)*NPHASE

        ALLOCATE( FUPWIND_MAT_ALL(NFIELD,NSMALL_COLM) )
        ALLOCATE( F_ALL(NFIELD, CV_NONODS), FEMF_ALL(NFIELD, CV_NONODS)  )
        F_ALL(1:NPHASE,            :)=T_ALL(1:NPHASE, :)
        F_ALL(1+NPHASE:  2*NPHASE, :)=TOLD_ALL(1:NPHASE, :)
        F_ALL(1+2*NPHASE:3*NPHASE, :)=DEN_ALL(1:NPHASE, :)
        F_ALL(1+3*NPHASE:4*NPHASE, :)=DENOLD_ALL(1:NPHASE, :)

        IF(IGOT_T2.NE.0) THEN
            F_ALL(1+4*NPHASE:5*NPHASE, :)=T2_ALL(1:NPHASE, :)
            F_ALL(1+5*NPHASE:6*NPHASE, :)=T2OLD_ALL(1:NPHASE, :)
        ENDIF
        ! femf:
        FEMF_ALL(1:NPHASE,            :)=FEMT_ALL(1:NPHASE, :)
        FEMF_ALL(1+NPHASE:  2*NPHASE, :)=FEMTOLD_ALL(1:NPHASE, :)
        FEMF_ALL(1+2*NPHASE:3*NPHASE, :)=FEMDEN_ALL(1:NPHASE, :)
        FEMF_ALL(1+3*NPHASE:4*NPHASE, :)=FEMDENOLD_ALL(1:NPHASE, :)

        IF(IGOT_T2.NE.0) THEN
            FEMF_ALL(1+4*NPHASE:5*NPHASE, :)=FEMT2_ALL(1:NPHASE, :)
            FEMF_ALL(1+5*NPHASE:6*NPHASE, :)=FEMT2OLD_ALL(1:NPHASE, :)
        ENDIF


        !Find upwind field values for limiting
        CALL CALC_ANISOTROP_LIM_VALS( F_ALL, FEMF_ALL, USE_FEMT, FUPWIND_MAT_ALL,  &
            NFIELD,CV_NONODS,CV_NLOC,TOTELE,CV_NDGLN, SMALL_FINDRM,&
            SMALL_COLM,NSMALL_COLM, X_NDGLN,X_NONODS,NDIM, X_ALL, XC_CV_ALL, use_reflect)

        ! make sure the diagonal is equal to the value:
        DO NOD=1,CV_NONODS
            IMID=SMALL_CENTRM(NOD)
            DO IFIELD=1,NFIELD
                FUPWIND_MAT_ALL( IFIELD, IMID)=F_ALL( IFIELD, NOD)
            END DO
        END DO


        TUPWIND_MAT_ALL(1:NPHASE,      :)=FUPWIND_MAT_ALL(1:NPHASE, :)
        TOLDUPWIND_MAT_ALL(1:NPHASE,   :)=FUPWIND_MAT_ALL(1+NPHASE:2*NPHASE, :)
        DENUPWIND_MAT_ALL(1:NPHASE,    :)=FUPWIND_MAT_ALL(2*NPHASE+1:3*NPHASE, :)
        DENOLDUPWIND_MAT_ALL(1:NPHASE, :)=FUPWIND_MAT_ALL(3*NPHASE+1:4*NPHASE, :)

        IF(IGOT_T2.NE.0) THEN
            T2UPWIND_MAT_ALL(1:NPHASE,    :)=FUPWIND_MAT_ALL(4*NPHASE+1:5*NPHASE, :)
            T2OLDUPWIND_MAT_ALL(1:NPHASE, :)=FUPWIND_MAT_ALL(5*NPHASE+1:6*NPHASE, :)
        ENDIF

        contains

            SUBROUTINE CALC_ANISOTROP_LIM_VALS( &
                ! Caculate the upwind values stored in matrix form...
                T_ALL, &
                FEMT_ALL, USE_FEMT, &
                TUPWIND_ALL, &
                NFIELD,NONODS,CV_NLOC,TOTELE,CV_NDGLN, &
                SMALL_FINDRM,SMALL_COLM,NSMALL_COLM, &
                X_NDGLN,X_NONODS,NDIM, &
                X_ALL, XC_CV_ALL, use_reflect)
                ! For the anisotropic limiting scheme we find the upwind values
                ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
                ! value for each node pair is stored in the matrices TUPWIND AND
                IMPLICIT NONE
                INTEGER, intent(in) :: NONODS,X_NONODS,TOTELE,CV_NLOC,NSMALL_COLM,NFIELD,NDIM
                REAL, DIMENSION( :, : ), intent( in ) :: T_ALL
                REAL, DIMENSION( :, : ), intent( in ) :: FEMT_ALL
                LOGICAL, intent( in ) :: USE_FEMT, use_reflect
                REAL, DIMENSION( :, : ), intent( inout ) :: TUPWIND_ALL
                INTEGER, DIMENSION( :  ), intent( in ) :: X_NDGLN
                INTEGER, DIMENSION( :  ), intent( in ) :: CV_NDGLN
                INTEGER, DIMENSION( : ), intent( in ) :: SMALL_FINDRM
                INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
                REAL, DIMENSION( :, : ), intent( in ) :: X_ALL
                REAL, DIMENSION( NDIM, NONODS ), intent( in ) :: XC_CV_ALL
                ! the centre of each CV is: XC_CV, YC_CV, ZC_CV

                ! Allocate memory for the interpolated upwind values
                real, dimension( :, : ), allocatable :: N!, NLX, NLY, NLZ
                real, dimension (:, :, :), allocatable :: NLX_ALL
                real, dimension( : ), allocatable :: WEIGHT, L1, L2, L3, L4
                integer, dimension( : ), allocatable :: SUB_NDGLNO, SUB_XNDGLNO, ndgln_p2top1
                INTEGER :: SUB_TOTELE, NGI,NLOC, ELE, IL_LOC, IQ_LOC, &
                    LOC_ELE, SUB_ELE, SUB_LIN_TOTELE


                ! **********************Calculate linear shape functions...
                IF(NDIM==1) THEN
                    NLOC=1
                    NGI=2
                ELSE IF(NDIM==2) THEN
                    NLOC=3
                    NGI=3
                ELSE IF(NDIM==3) THEN
                    NLOC=4
                    NGI=4
                ENDIF
                ALLOCATE( N(NLOC,NGI))
                ALLOCATE( WEIGHT(NGI) )
                ALLOCATE( L1(NGI), L2(NGI), L3(NGI), L4(NGI) )
                allocate(NLX_ALL(size(X_ALL,1), NLOC, NGI))
                !
                ! Shape functions for triangles and tets...
                CALL TRIQUAold( L1, L2, L3, L4, WEIGHT, ndim==3, NGI )
                ! Work out the shape functions and there derivatives...
                call SHATRInew(L1, L2, L3, L4, WEIGHT,  NLOC,NGI,  N,NLX_ALL)
                ! ******************************************************************
                ! Calculate the sub elements for quadratic element SUB_NDGLNO ...
                IF(CV_NLOC==NLOC) THEN
                    SUB_TOTELE=TOTELE
                ELSE
                    IF(NDIM==1) THEN
                        sub_lin_totele=2
                    ELSE IF(NDIM==2) THEN
                        sub_lin_totele=4
                    ELSE IF(NDIM==3) THEN
                        sub_lin_totele=8
                    ENDIF
                    SUB_TOTELE= sub_lin_totele * totele

                    allocate( ndgln_p2top1( sub_lin_totele*nloc ) ) ; ndgln_p2top1 = 0
                    call conv_quad_to_lin_tri_tet( ndgln_p2top1, nloc, cv_nloc, sub_lin_totele )

                ENDIF

                ALLOCATE( SUB_NDGLNO( SUB_TOTELE*NLOC ) )
                ALLOCATE( SUB_XNDGLNO( SUB_TOTELE*NLOC ) )

                IF ( CV_NLOC==NLOC ) THEN
                    SUB_NDGLNO = CV_NDGLN
                    SUB_XNDGLNO = X_NDGLN
                ELSE

                    SUB_ELE=0
                    DO ELE = 1, TOTELE
                        DO LOC_ELE = 1, SUB_LIN_TOTELE

                            SUB_ELE = SUB_ELE + 1

                            DO IL_LOC = 1, NLOC
                                IQ_LOC = ndgln_p2top1( (loc_ELE-1)*NLOC + IL_LOC )
                                SUB_NDGLNO( (sub_ele-1)*nloc + il_loc ) = cv_ndgln( (ele-1)*cv_nloc + iq_loc )
                                SUB_XNDGLNO( (sub_ele-1)*nloc + il_loc ) = x_ndgln( (ele-1)*cv_nloc + iq_loc )
                            END DO

                        END DO
                    END DO
                    deallocate( ndgln_p2top1 )
                END IF

                ! Calculate the sub elements for quadratic element SUB_NDGLNO ...
                ! ******************************************************************
                CALL CALC_ANISOTROP_LIM_VALS2( &
                    ! Caculate the upwind values stored in matrix form...
                    T_ALL, &
                    FEMT_ALL, USE_FEMT, &
                    TUPWIND_ALL,  &
                    NFIELD, NONODS, NLOC, NGI, SUB_TOTELE, SUB_NDGLNO, &
                    SMALL_FINDRM,SMALL_COLM, NSMALL_COLM, &
                    SUB_XNDGLNO, X_NONODS, NDIM, &
                    X_ALL, XC_CV_ALL, &
                    N, NLX_ALL, WEIGHT, use_reflect)


                !    DEALLOCATE( N, NLX, NLY, NLZ, L1, L2, L3, L4, &
                !         WEIGHT, SUB_NDGLNO, SUB_XNDGLNO )
                DEALLOCATE( N, NLX_ALL, L1, L2, L3, L4, &
                    WEIGHT, SUB_NDGLNO, SUB_XNDGLNO )
                RETURN
            END SUBROUTINE CALC_ANISOTROP_LIM_VALS


            SUBROUTINE CALC_ANISOTROP_LIM_VALS2( &
                ! Caculate the upwind values stored in matrix form...
                T_ALL, &
                FEMT_ALL, USE_FEMT, &
                TUPWIND_ALL,  &
                NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
                FINDRM,COLM,NCOLM, &
                X_NDGLN,X_NONODS,NDIM, &
                X_ALL, XC_CV_ALL,  &
                N,NLX_ALL, WEIGHT, use_reflect)
                ! For the anisotropic limiting scheme we find the upwind values
                ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
                ! value for each node pair is stored in the matrices TUPWIND AND
                IMPLICIT NONE
                INTEGER, intent(in) :: NONODS,X_NONODS,TOTELE,NLOC,NGI,NCOLM,NFIELD,NDIM
                REAL, DIMENSION( :,: ), intent( in ) :: T_ALL
                REAL, DIMENSION(  :,: ), intent( in ) :: FEMT_ALL
                LOGICAL, intent( in ) :: USE_FEMT, use_reflect
                REAL, DIMENSION( :,:  ), intent( inout ) :: TUPWIND_ALL
                INTEGER, DIMENSION( : ), INTENT(IN) :: NDGLNO,X_NDGLN
                INTEGER, DIMENSION( : ), INTENT(IN) :: FINDRM,COLM

                REAL, DIMENSION(:,:), intent( in ) :: X_ALL
                REAL, DIMENSION( NDIM, NONODS ), intent( in ) :: XC_CV_ALL
                REAL, DIMENSION(NLOC,NGI), INTENT(IN) :: N!,NLX,NLY,NLZ
                REAL, DIMENSION(:,:,:), INTENT(IN) :: NLX_ALL!DIMENSION(NDIM, NLOC,NGI)
                REAL, DIMENSION(NGI), INTENT(IN) :: WEIGHT
                !Local variables

                INTEGER, DIMENSION( : ), ALLOCATABLE, SAVE :: ELEMATPSI
                REAL, DIMENSION( :  ), ALLOCATABLE, SAVE :: ELEMATWEI
                LOGICAL, SAVE :: STORE_ELE=.TRUE., RET_STORE_ELE=.FALSE.
                ! Allocate memory for the interpolated upwind values
                LOGICAL, PARAMETER :: BOUND  = .TRUE.! limiting options
                logical:: REFLECT ! limiting options
                INTEGER, DIMENSION( : ), allocatable :: NOD_FINDELE,NOD_COLELE, NLIST, INLIST, DUMMYINT
                REAL, DIMENSION( : ), allocatable :: DUMMYREAL
                INTEGER MXNCOLEL,NCOLEL,adapt_time_steps
                REAL current_time
                !Reflect option defined from diamond
                REFLECT = use_reflect
                ! Over-estimate the size of the COLELE array
                MXNCOLEL=20*TOTELE+500

                ALLOCATE( NOD_FINDELE(X_NONODS+1) )
                ALLOCATE( NOD_COLELE(MXNCOLEL) )
                ALLOCATE( NLIST(X_NONODS) )
                ALLOCATE( INLIST(X_NONODS) )

                ! Calculate node element list - moved from (I)FINPTS
                CALL PHILNODELE(X_NONODS,NOD_FINDELE,NOD_COLELE, &
                    NCOLEL,MXNCOLEL, &
                    TOTELE,NLOC,X_NDGLN, &
                    NLIST,INLIST)

                IF( STORE_ELE ) THEN

                    ALLOCATE( ELEMATPSI( NCOLM ) )
                    ALLOCATE( ELEMATWEI( NCOLM * NLOC ) )

                    CALL FINPTSSTORE(T_ALL,FEMT_ALL,USE_FEMT,NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
                        TUPWIND_ALL,FINDRM,COLM,NCOLM,NDIM, &
                        X_NDGLN,X_NONODS, &
                        X_ALL, XC_CV_ALL, &
                        N,NLX_ALL, WEIGHT, &
                        NOD_FINDELE,NOD_COLELE,NCOLEL, &
                        ELEMATPSI,ELEMATWEI,1, &
                        BOUND, REFLECT)

                ELSE IF( RET_STORE_ELE ) THEN

                    ! Find the weights for the interpolation
                    ! This does depend on the solns T when BOUND...
                    CALL GETSTOREELEWEI(T_ALL,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
                        TUPWIND_ALL,FINDRM,COLM,NCOLM,BOUND, &
                        ELEMATPSI,ELEMATWEI)

                ELSE

                    ! Assume we have not stored anything (elements or weights)...
                    ALLOCATE(DUMMYINT(NCOLM))
                    ALLOCATE(DUMMYREAL(NCOLM*NLOC))

                    CALL FINPTSSTORE(T_ALL,FEMT_ALL,USE_FEMT,NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
                        TUPWIND_ALL,FINDRM,COLM,NCOLM,NDIM, &
                        X_NDGLN,X_NONODS, &
                        X_ALL, XC_CV_ALL, &
                        N,NLX_ALL, WEIGHT, &
                        NOD_FINDELE,NOD_COLELE,NCOLEL, &
                        DUMMYINT,DUMMYREAL,0, &
                        BOUND, REFLECT)

                    DEALLOCATE(DUMMYINT,DUMMYREAL)

                ENDIF

                store_ele = .false. ; ret_store_ele = .true.
                if( have_option( '/mesh_adaptivity/hr_adaptivity') ) then
                    if( have_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps') ) then
                        call get_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps', &
                            adapt_time_steps )
                        if( mod( timestep, adapt_time_steps ) == 0 ) store_ele = .true.
                    else if (have_option( '/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI') ) then
                        STORE_ELE=.TRUE.; RET_STORE_ELE=.FALSE.
                    end if
                elseif( have_option( '/mesh_adaptivity/hr_adaptivity_prescribed_metric') ) then
                    if( have_option( '/mesh_adaptivity/hr_adaptivity_prescribed_metric/period_in_timesteps') ) then
                        call get_option( '/mesh_adaptivity/hr_adaptivity_prescribed_metric/period_in_timesteps', &
                            adapt_time_steps )
                    end if
                    if( mod( timestep, adapt_time_steps ) == 0 ) store_ele = .true.
                elseif( have_option( '/mesh_adaptivity/prescribed_adaptivity' ) ) then
                    call get_option( '/timestepping/current_time', current_time )
                    if( do_adapt_state_prescribed( current_time ) ) store_ele = .true.
                end if
                if ( store_ele ) then
                    ret_store_ele = .false.
                    deallocate( elematpsi, elematwei )
                end if

                DEALLOCATE( NOD_FINDELE, NOD_COLELE, NLIST, INLIST )

            END SUBROUTINE CALC_ANISOTROP_LIM_VALS2





            SUBROUTINE GETSTOREELEWEI(PSI_ALL,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
                &     MATPSI_ALL,FINDRM,COLM,NCOLM,BOUND,&
                &     ELEMATPSI,ELEMATWEI)
                ! use the stored interpolation coeffs to caclulate MATPSI.
                !     This sub finds the matrix values MATPSI for a given point on the
                !     stencil
                IMPLICIT NONE
                REAL FRALINE
                LOGICAL BOUND
                PARAMETER(FRALINE=0.001)
                INTEGER, intent(in) :: NFIELD,NONODS,NLOC,TOTELE,NDGLNO(TOTELE*NLOC)
                REAL, DIMENSION(:,:), INTENT(IN) :: PSI_ALL
                INTEGER, INTENT(IN) :: NCOLM
                INTEGER, DIMENSION(:), INTENT(IN) :: FINDRM
                INTEGER, DIMENSION(:), INTENT(IN) :: COLM
                REAL, DIMENSION(:,:), INTENT(INOUT) :: MATPSI_ALL
                INTEGER, DIMENSION(:), INTENT(IN) :: ELEMATPSI
                REAL, DIMENSION(NCOLM*NLOC),  INTENT(IN) ::  ELEMATWEI
                !  LOCAL VARIABLES...
                INTEGER NOD,COUNT,ELEWIC,ILOC,INOD,IFIELD
                REAL RMATPSI
                REAL, ALLOCATABLE, DIMENSION(:,:)::MINPSI
                REAL, ALLOCATABLE, DIMENSION(:,:)::MAXPSI

                ALLOCATE(MINPSI(NFIELD, TOTELE))
                ALLOCATE(MAXPSI(NFIELD, TOTELE))

                if ( bound ) then

                    ! find the max and min local to each element...
                    CALL MINMAXELEWIC( PSI_ALL,NONODS,NLOC,TOTELE,NDGLNO, &
                        &     FINDRM,COLM,NCOLM,&
                        &     MINPSI,MAXPSI )
                end if
                do NOD = 1, NONODS
                    do COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
                        IF(NOD.NE.COLM(COUNT)) THEN
                            ELEWIC = ELEMATPSI( COUNT )
                            DO IFIELD = 1, NFIELD
                                RMATPSI=0.0
                                DO ILOC = 1, NLOC
                                    INOD = NDGLNO( (ELEWIC-1)*NLOC + ILOC )
                                    RMATPSI = RMATPSI + ELEMATWEI( (COUNT-1)*NLOC+ILOC) * PSI_ALL(IFIELD,INOD)
                                END DO

                                RMATPSI   =PSI_ALL(IFIELD,NOD)   &
                                    +(1./FRALINE)*(RMATPSI   -PSI_ALL(IFIELD,NOD))

                                ! make locally bounded...
                                if ( bound ) then
                                    MATPSI_ALL(IFIELD, COUNT)   &
                                        =MAX(MIN(RMATPSI,   MAXPSI(IFIELD, ELEWIC)),   &
                                        &                            MINPSI(IFIELD, ELEWIC))
                                else
                                    MATPSI_ALL(IFIELD, COUNT)   =RMATPSI
                                end if
                            END DO
                        END IF
                    END DO
                END DO

                !    if ( bound ) then
                DEALLOCATE( MINPSI, MAXPSI )
                !    end if

                RETURN

            end subroutine getstoreelewei

            SUBROUTINE MINMAXELEWIC(PSI_ALL,NONODS,NLOC,TOTELE,NDGLNO, &
                &     FINDRM,COLM,NCOLM,&
                &     MINPSI,MAXPSI)
                ! This sub calculates the max and min values of PSI in local vacinity of
                ! an element.
                IMPLICIT NONE
                INTEGER, intent(in) :: NONODS,NLOC,TOTELE,NDGLNO(TOTELE*NLOC)
                REAL, DIMENSION(:,:), INTENT(IN) :: PSI_ALL
                INTEGER, INTENT(IN) :: NCOLM
                INTEGER, INTENT(IN) :: FINDRM(NONODS+1),COLM(NCOLM)
                !    REAL, INTENT(INOUT) :: MINPSI(TOTELE*NFIELD),MAXPSI(TOTELE*NFIELD)
                REAL, DIMENSION(:,:), INTENT(INOUT) :: MINPSI,MAXPSI
                !  LOCAL VARIABLES...
                INTEGER ELEWIC,ILOC
                INTEGER KNOD,COUNT2,JNOD

                MINPSI   =1.E+20
                MAXPSI   =-1.E+20
                ! find the max and min local to each element...
                DO ELEWIC=1,TOTELE! Was loop
                    DO ILOC=1,NLOC! Was loop
                        KNOD=NDGLNO((ELEWIC-1)*NLOC+ILOC)
                        ! Search around node KNOD for max and min PSI...
                        DO COUNT2 = FINDRM(KNOD), FINDRM(KNOD+1)-1
                            JNOD = COLM( COUNT2 )
                            !DO IFIELD = 1, NFIELD
                            MINPSI( :, ELEWIC )  &
                                = MIN( PSI_ALL(:, JNOD), MINPSI(:, ELEWIC) )
                            MAXPSI( :, ELEWIC )  &
                                = MAX( PSI_ALL(:, JNOD), MAXPSI(:, ELEWIC) )
                        !                     = MAX( PSI_ALL(JNOD+(IFIELD-1)*NONODS), MAXPSI(ELEWIC+(IFIELD-1)*TOTELE) )
                           !END DO
                        END DO
                    END DO
                END DO

                !ewrite(3,*) '***M-m', MAXPSI-MINPSI

                RETURN

            end subroutine minmaxelewic

            !
            !
            !
            !
            SUBROUTINE FINPTSSTORE(PSI_ALL,FEMPSI_ALL,USE_FEMPSI,NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
                MATPSI_ALL,FINDRM,COLM,NCOLM,NDIM, &
                X_NDGLN,X_NONODS, &
                X_ALL, XC_CV_ALL, &
                N,NLX_ALL, WEIGHT,&
                !     work space...
                FINDELE,COLELE,NCOLEL,&
                ELEMATPSI,ELEMATWEI,IGETSTOR,&
                BOUND, REFLECT)
                !     This sub finds the matrix values MATPSI for a given point on the
                !     stencil
                ! IF IGETSTOR=1 then get ELEMATPSI,ELEMATWEI.
                IMPLICIT NONE
                LOGICAL BOUND,REFLECT
                ! IF REFLECT then use a reflection condition at boundary to
                ! do limiting.
                INTEGER, intent(in) :: NFIELD,NONODS,NLOC,NGI,TOTELE,NDIM,X_NONODS
                INTEGER, dimension(TOTELE*NLOC),intent(in) :: NDGLNO
                REAL, dimension(:,:), intent(in) :: PSI_ALL
                REAL, dimension(:,:), intent(in) :: FEMPSI_ALL
                LOGICAL, intent(in) :: USE_FEMPSI
                INTEGER, intent(in) :: NCOLM,NCOLEL
                INTEGER, dimension(NONODS+1), intent(in) :: FINDRM
                INTEGER, dimension(NCOLM),intent(in) :: COLM
                REAL, dimension(:,:), intent(inout) :: MATPSI_ALL
                INTEGER, dimension(TOTELE*NLOC),  intent(in) :: X_NDGLN
                REAL, dimension(:,:), intent(in) :: X_ALL
                REAL, DIMENSION( NDIM, NONODS ), intent( in ) :: XC_CV_ALL
                REAL, dimension(NLOC,NGI), intent(in) :: N!,NLX,NLY,NLZ
                REAL, dimension(:, :,:), intent(in) :: NLX_ALL!dimension(NDIM, NLOC,NGI)
                REAL, dimension(:), intent(in) :: WEIGHT!dimenson(NGI)
                !     work space...
                INTEGER, dimension(X_NONODS+1),intent(in) :: FINDELE
                INTEGER, dimension(NCOLEL),intent(in) :: COLELE
                INTEGER, intent(in) :: IGETSTOR
                INTEGER, dimension(NCOLM*IGETSTOR), intent(inout) :: ELEMATPSI
                REAL, dimension(NCOLM*NLOC*IGETSTOR), intent(inout) :: ELEMATWEI
                ! ELEWIC is the element to do interpolation from
                ! LOCCORDSK contains the weights.
                !     Local variables...
                INTEGER NOD,COUNT,NODI,NODJ,ILOC,GI,ELE
                INTEGER ELEWIC,XNOD,XNODJ
                REAL LOCCORDSK(NLOC)
                REAL INVH,LENG
                !     work space...
                real, pointer :: VOLUME
                real, target :: VOLUME2
                real, dimension(:), pointer :: DETWEI, RA !dimension(NGI)
                real, dimension(:), allocatable, target :: DETWEI2, RA2
                real, dimension (size(X_ALL,1)) :: NORMX1_ALL
                real, dimension(:, :, :), pointer :: NX_ALL ! dimension (size(X_ALL,1), NLOC, NGI)
                real, dimension(:,:,:), allocatable, target :: NX_ALL2
                real, dimension (size(X_ALL,1), NONODS) :: NORMX_ALL
                REAL, ALLOCATABLE, DIMENSION(:)::MLUM
                !    REAL, ALLOCATABLE, DIMENSION(:)::MINPSI
                !    REAL, ALLOCATABLE, DIMENSION(:)::MAXPSI
                REAL, ALLOCATABLE, DIMENSION(:,:)::MINPSI
                REAL, ALLOCATABLE, DIMENSION(:,:)::MAXPSI
                INTEGER, ALLOCATABLE, DIMENSION(:)::NOD2XNOD
                real, dimension(size(X_ALL,1)) :: X1_ALL, X2_ALL

                ! Initialisation
                ELEMATPSI = 0
                ELEMATWEI = 0.0
                allocate(DETWEI2(NGI))
                allocate(RA2(NGI))
                allocate(NX_ALL2(size(X_ALL,1), NLOC, NGI))
                VOLUME=>VOLUME2; DETWEI=>DETWEI2;
                RA=>RA2; NX_ALL=>NX_ALL2

                NORMX1_ALL=0.0
                IF(REFLECT) THEN
                    !     calculate normals...********************
                    ALLOCATE(MLUM(NONODS))
                    NORMX_ALL = 0
                    MLUM(1:NONODS) = 0.0
                    DO ELE=1,TOTELE! Was loop
                        call DETNLXR( ELE, X_ALL, X_NDGLN, NLOC, NGI, &
                            N, NLX_ALL, WEIGHT, DETWEI, RA2, VOLUME2, .false., NX_ALL2)
                        DO ILOC=1,NLOC! Was loop
                            NODI=NDGLNO((ELE-1)*NLOC+ILOC)
                            DO GI=1,NGI! Was loop
                                NORMX_ALL(:,NODI) = NORMX_ALL(:,NODI) + NX_ALL(:,ILOC,GI) * DETWEI(GI)
                                MLUM(NODI) =MLUM(NODI) +N(ILOC,GI) *DETWEI(GI)
                            END DO
                        END DO
                    END DO
                    !     Renormalise
                    DO NODI=1,NONODS! Was loop
                        INVH = SUM(ABS(NORMX_ALL(:,NODI)))/MLUM(NODI)

                        !          INVH=(ABS(NORMX_ALL(1,NODI))+ABS(NORMX_ALL(2,NODI))+ABS(NORMX_ALL(3,NODI)))&
                        !               &          /MLUM(NODI)
                        IF(INVH.GT.1.E-5) THEN
                            LENG = sqrt(dot_product(NORMX_ALL(:,NODI),NORMX_ALL(:,NODI)))
                            NORMX_ALL(:,NODI) = NORMX_ALL(:,NODI) / LENG
                        ELSE
                            NORMX_ALL(:,NODI) = 0.0
                        END IF
                    END DO
                ENDIF

                ALLOCATE(MINPSI(NFIELD, TOTELE))
                ALLOCATE(MAXPSI(NFIELD, TOTELE))

                IF(BOUND) THEN
                    ! find the max and min local to each element...
                    CALL MINMAXELEWIC(PSI_ALL,NONODS,NLOC,TOTELE,NDGLNO, &
                        &     FINDRM,COLM,NCOLM,&
                        &     MINPSI,MAXPSI)
                ENDIF

                !
                !     Calculate node element list.

                ALLOCATE(NOD2XNOD(NONODS))
                DO ELE=1,TOTELE! Was loop
                    DO ILOC=1,NLOC! Was loop
                        NOD =NDGLNO((ELE-1)*NLOC+ILOC)
                        XNOD=X_NDGLN((ELE-1)*NLOC+ILOC)
                        NOD2XNOD(NOD)=XNOD
                    END DO
                END DO
                !

                MATPSI_ALL=0.
                DO NOD=1,NONODS! Was loop 10

                    XNOD=NOD2XNOD(NOD)
                    !
                    DO COUNT=FINDRM(NOD ),FINDRM(NOD+1)-1! Was loop 20

                        NODJ=COLM(COUNT)
                        XNODJ=NOD2XNOD(NODJ)
                        !
                        IF(NOD.NE.NODJ) THEN

                            IF(REFLECT) THEN
                                NORMX1_ALL = NORMX_ALL(:,NOD)
                            ENDIF
                            IF(NONODS.NE.X_NONODS) THEN ! Its a DG soln field...
                                X1_ALL = XC_CV_ALL(:,NOD)

                                X2_ALL = XC_CV_ALL(:, NODJ)
                            ELSE
                                X1_ALL = X_ALL(:,XNOD)
                                X2_ALL = X_ALL(:,XNODJ)
                            ENDIF
                            CALL MATPTSSTORE(MATPSI_ALL,COUNT,NFIELD,NOD,XNOD,&
                                PSI_ALL,FEMPSI_ALL,USE_FEMPSI,NONODS,X_NONODS,&
                                NLOC,TOTELE,X_NDGLN,NDGLNO,&
                                X1_ALL,&
                                X2_ALL,&
                                NORMX1_ALL,&
                                X_ALL,&
                                !     work space...
                                FINDELE,COLELE,NCOLEL, &
                                MINPSI,MAXPSI, &
                                ELEWIC,LOCCORDSK,BOUND,REFLECT,NDIM)
                            IF(IGETSTOR.EQ.1) THEN
                                ELEMATPSI(COUNT)=ELEWIC
                                DO ILOC=1,NLOC! Was loop
                                    ELEMATWEI((COUNT-1)*NLOC+ILOC)=LOCCORDSK(ILOC)
                                END DO
                            ENDIF
                        ENDIF

                    END DO ! Was loop 20
                END DO ! Was loop 10

                deallocate(DETWEI2)
                deallocate(RA2)
                deallocate(NX_ALL2)

                RETURN

            end subroutine finptsstore
            !
            !
            !
            SUBROUTINE MATPTSSTORE(MATPSI_ALL,COUNT,NFIELD,NOD,XNOD,&
                PSI_ALL,FEMPSI_ALL,USE_FEMPSI,NONODS,X_NONODS,&
                NLOC,TOTELE,X_NDGLN,NDGLNO,&
                X1_ALL,&
                X2_ALL,&
                NORMX1_ALL,&
                X_ALL,&
                !     work space...
                FINDELE,COLELE,NCOLEL,&
                MINPSI,MAXPSI,  &
                ELEWIC,LOCCORDSK,BOUND,REFLECT,NDIM)
                !     This sub calculates the value of PSI that would be at the
                !     other side of the stencil if we had a linear variation and within
                !     a single element.
                ! IF BOUND then make locally bounded.
                IMPLICIT NONE
                REAL INFINY,FRALINE2
                LOGICAL, intent(in) :: REFLECT
                ! IF REFLECT then use a reflection condition at boundary to
                ! do limiting.
                PARAMETER(INFINY=1.E+20,FRALINE2=0.001)
                LOGICAL, intent(in) :: BOUND
                INTEGER, intent(in) :: COUNT,NFIELD,NOD,XNOD,NONODS,X_NONODS,NLOC,TOTELE,NDIM
                REAL, dimension(:,:), intent(in) :: PSI_ALL
                REAL, dimension(:,:), intent(in) :: FEMPSI_ALL
                LOGICAL, intent(in) :: USE_FEMPSI
                REAL, dimension(:,:), intent(inout) :: MATPSI_ALL
                INTEGER, intent(in) :: X_NDGLN(NLOC*TOTELE),NDGLNO(NLOC*TOTELE)
                !      REAL, intent(in) :: X1,Y1,Z1,X2,Y2,Z2,NORMX1,NORMY1,NORMZ1
                real, dimension(:) :: X1_ALL, X2_ALL, NORMX1_ALL!dimension(NDIM)
                REAL, dimension(:,:), intent(in) :: X_ALL
                INTEGER, intent(in) :: NCOLEL
                INTEGER, intent(in) :: FINDELE(X_NONODS+1),COLELE(NCOLEL)
                !      REAL, intent(in) :: MINPSI(TOTELE*NFIELD),MAXPSI(TOTELE*NFIELD)
                REAL, DIMENSION(:, :), intent(in) :: MINPSI,MAXPSI
                INTEGER, intent(inout) :: ELEWIC
                REAL, intent(inout) :: LOCCORDSK(NLOC)
                !
                !     Local variables...
                REAL, dimension(4):: LOCCORDS
                INTEGER , dimension(4) :: LOCNODS,LOCNODSK
                INTEGER, dimension(4) :: NLOCNODS,NLOCNODSK
                INTEGER :: ELE,ILOC,IFIELD,COUNT2
                REAL :: MINCOR,MINCORK,RSUM
                REAL :: DIST12,RN,RMATPSI
                REAL :: FRALINE
                LOGICAL IS_DG
                !The dimension of the variables below should be NDIM, however, due to cross products
                !we need three dimensions
                REAL, dimension(3) ::  XC_ALL, VX_ALL, REFX_ALL, REFX2_ALL, T2X_ALL, T1X_ALL, AUXNORMX1_ALL


                IS_DG=NONODS.NE.X_NONODS
                !
                FRALINE=FRALINE2
                IF(IS_DG) FRALINE=1.0
                XC_ALL = 0.
                XC_ALL(1:NDIM) = X1_ALL - FRALINE*(X2_ALL-X1_ALL)
                !print *, "XC_ALL before reflect", XC_ALL

                IF(REFLECT) THEN
                    IF(SUM(ABS(NORMX1_ALL)).NE.0.0) THEN
                        !  if (XC,YC,ZC) is outside the domain
                        !     The rotation matrix in 3-D is R=
                        !     NX    NY    NZ
                        !     T1X   T1Y   T1Z
                        !     T2X   T2Y   T2Z
                        !
                        VX_ALL = 0.
                        VX_ALL(1:NDIM) = X1_ALL - X2_ALL
                        !
                        AUXNORMX1_ALL = 0.
                        AUXNORMX1_ALL(1:NDIM) = NORMX1_ALL
                        CALL XPROD(T2X_ALL, AUXNORMX1_ALL, VX_ALL)
                        !
                        !DIST12=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
                        DIST12 = SQRT(DOT_PRODUCT(X1_ALL-X2_ALL,X1_ALL-X2_ALL))
                        RN = SQRT(DOT_PRODUCT(T2X_ALL,T2X_ALL))
                        IF(RN.LT.(1.E-5)*DIST12) THEN
                            !     Simply have VX,VY,VZ going in the opposite direction...
                            XC_ALL(1:NDIM) = X1_ALL - VX_ALL(1:NDIM)*FRALINE
                        ELSE
                            T2X_ALL = T2X_ALL/RN
                            !     T1=Nx (-T2)
                            CALL XPROD(T1X_ALL, AUXNORMX1_ALL, -T2X_ALL)
                            !
                            REFX2_ALL(1) = SUM(NORMX1_ALL(1:NDIM)*VX_ALL(1:NDIM))
                            REFX2_ALL(2) = SUM(T1X_ALL(1:NDIM)*VX_ALL(1:NDIM))


                            !     Reflect...
                            REFX2_ALL(1) = - REFX2_ALL(1)
                            !     MAP BACK USING R^T

                            !     (REFX,REFY,REFZ) is the reflected direction...
                            REFX_ALL(1) =  NORMX1_ALL(1) * REFX2_ALL(1) + T1X_ALL(1) * REFX2_ALL(2)
                            REFX_ALL(2) =  NORMX1_ALL(2) * REFX2_ALL(1) + T1X_ALL(2) * REFX2_ALL(2)
                            IF (NDIM==3) THEN
                                !SOME MORE THINGS NEED TO BE ADDED
                                REFX2_ALL(3) = SUM(T2X_ALL(:)*VX_ALL)

                                REFX_ALL(1) = REFX_ALL(1) + T2X_ALL(1) * REFX2_ALL(3)
                                REFX_ALL(2) = REFX_ALL(2) + T2X_ALL(2) * REFX2_ALL(3)
                                REFX_ALL(3) =  NORMX1_ALL(3) * REFX2_ALL(1) + T1X_ALL(3) * REFX2_ALL(2)+ T2X_ALL(3) * REFX2_ALL(3)
                            END IF

                            XC_ALL(1:NDIM) = X1_ALL(1:NDIM) + REFX_ALL(1:NDIM)*FRALINE
                        ENDIF

                    ENDIF
                ENDIF
                !
                MINCORK=-INFINY
                !
                DO COUNT2=FINDELE(XNOD),FINDELE(XNOD+1)-1! Was loop 10
                    ELE=COLELE(COUNT2)
                    !
                    NLOCNODS(1)=NDGLNO((ELE-1)*NLOC+1)
                    NLOCNODS(2)=NDGLNO((ELE-1)*NLOC+2)
                    NLOCNODS(3)=NDGLNO((ELE-1)*NLOC+3)

                    !
                    LOCNODS(1)=X_NDGLN((ELE-1)*NLOC+1)
                    LOCNODS(2)=X_NDGLN((ELE-1)*NLOC+2)
                    LOCNODS(3)=X_NDGLN((ELE-1)*NLOC+3)
                    !
                    ! Calculate the local coord but with 4th point replaced by INOD...
                    ! Find local coords LOCCORDS of point INOD corresponding to these nodes LOCNODS...

                    IF (NDIM==3) THEN
                        !Two coordinates missing if 3D
                        NLOCNODS(4)=NDGLNO((ELE-1)*NLOC+4)
                        LOCNODS(4)=X_NDGLN((ELE-1)*NLOC+4)

                        CALL TRILOCCORDS(XC_ALL(1),XC_ALL(2),XC_ALL(3), &
                            LOCCORDS(1),LOCCORDS(2),LOCCORDS(3),LOCCORDS(4),&
                            !     The 4 corners of the tet...
                            X_ALL(1,LOCNODS(1)),X_ALL(2,LOCNODS(1)),X_ALL(3,LOCNODS(1)),&
                            X_ALL(1,LOCNODS(2)),X_ALL(2,LOCNODS(2)),X_ALL(3,LOCNODS(2)),&
                            X_ALL(1,LOCNODS(3)),X_ALL(2,LOCNODS(3)),X_ALL(3,LOCNODS(3)),&
                            X_ALL(1,LOCNODS(4)),X_ALL(2,LOCNODS(4)),X_ALL(3,LOCNODS(4)) )
                    ELSE
                        CALL TRILOCCORDS2D(XC_ALL(1),XC_ALL(2), &
                            LOCCORDS(1),LOCCORDS(2),LOCCORDS(3),&
                            !     The 3 corners of the tri...
                            X_ALL(1,LOCNODS(1)),X_ALL(2,LOCNODS(1)),&
                            X_ALL(1,LOCNODS(2)),X_ALL(2,LOCNODS(2)),&
                            X_ALL(1,LOCNODS(3)),X_ALL(2,LOCNODS(3)) )
                    END IF

                    MINCOR=MINVAL( LOCCORDS(1:NLOC) )
                    !          print *,'ele,LOCCORDS(1:NLOC):',ele,LOCCORDS(1:NLOC)

                    IF(MINCOR.GT.MINCORK) THEN
                        MINCORK=MINCOR
                        DO ILOC=1,NLOC! Was loop
                            LOCCORDSK(ILOC)=LOCCORDS(ILOC)
                            LOCNODSK(ILOC)=LOCNODS(ILOC)
                            NLOCNODSK(ILOC)=NLOCNODS(ILOC)
                        END DO
                        ELEWIC=ELE
                    ENDIF
                END DO ! Was loop 10
                !        stop 677


                !     Set all the negative basis to zero and re-normalise
                !     to put on the face of an element...
                RSUM=0.0
                DO ILOC=1,NLOC! Was loop
                    LOCCORDSK(ILOC)=MAX(0.0,LOCCORDSK(ILOC))
                    RSUM=RSUM+LOCCORDSK(ILOC)
                END DO
                IF(RSUM.LT.1.E-5) THEN ! Just in case RSUM=0.0
                    LOCCORDSK(1:NLOC)=1.0/REAL(NLOC)
                ELSE
                    DO ILOC=1,NLOC! Was loop
                        LOCCORDSK(ILOC)=LOCCORDSK(ILOC)/RSUM
                    END DO
                ENDIF
                DO IFIELD=1,NFIELD
                    RMATPSI=0.0
                    DO ILOC=1,NLOC! Was loop
                        IF(USE_FEMPSI) THEN
                            RMATPSI   =RMATPSI  +LOCCORDSK(ILOC)*FEMPSI_ALL(IFIELD,NLOCNODSK(ILOC))
                        ELSE
                            RMATPSI   =RMATPSI  +LOCCORDSK(ILOC)*PSI_ALL(IFIELD, NLOCNODSK(ILOC))
                        ENDIF
                    END DO
                    !     Exaduate difference by a factor of 1/FRALINE.
                    IF(USE_FEMPSI) THEN
                        RMATPSI   = FEMPSI_ALL(IFIELD,  NOD )  &
                            + (1./FRALINE) * ( RMATPSI - FEMPSI_ALL( IFIELD, NOD) )
                    ELSE
                        RMATPSI   = PSI_ALL( IFIELD, NOD )  &
                            + (1./FRALINE) * ( RMATPSI - PSI_ALL(IFIELD,  NOD) )
                    ENDIF

                    !     Now correct to make sure that we get a bounded soln...
                    IF(BOUND) THEN
                        RMATPSI   =MAX(MIN(RMATPSI,   MAXPSI(IFIELD, ELEWIC)),   MINPSI(IFIELD, ELEWIC))
                    ENDIF
                    MATPSI_ALL(IFIELD, COUNT)   =RMATPSI
                END DO
                !
                RETURN

            end subroutine matptsstore
            !
            !
            !
            SUBROUTINE PHILNODELE(NONODS,FINDELE,COLELE, &
                NCOLEL,MXNCOLEL, &
                TOTELE,NLOC,NDGLNO, &
                NLIST,INLIST)
                !=================================================================
                ! This sub calculates the node to element list FINDELE,COLELE
                !
                ! Note NLIST and INLIST are only used locally but are passed
                ! down from parent routine where they are dynamically allocated.
                !
                ! INPUTS:
                ! ------
                ! NDGLNO  - List of global node numbers
                !
                ! OUTPUTS:
                ! -------
                ! COLELE  - This is a list of the element numbers that each node
                !           belongs to.  So it lists all elements for node 1, then
                !           all elements for node 2, and so on...
                ! FINDELE - is the pointer to the place in COLELE that gives the
                !           first element associated with a given global node
                !
                ! Called from subroutines IFINPTS and FINPTS, which are
                ! subroutines of CONSTRUCT_ADVECTION_DIFFUSION_CV
                !
                ! Description                                   Programmer      Date
                ! ==================================================================
                ! Original version..................................CCP   2013-28-01
                !
                !================================================================
                IMPLICIT NONE
                integer, intent( in ) :: NONODS,MXNCOLEL,TOTELE,NLOC
                integer, intent( inout ) :: NCOLEL
                integer, dimension( : ), intent( inout ) :: FINDELE
                integer, dimension( : ), intent( inout ) :: COLELE
                integer, dimension( : ), intent( in ) :: NDGLNO
                integer, dimension( : ), intent( inout ) :: NLIST,INLIST
                !     Local variables...
                INTEGER NOD,ELE,ILOC,COUNT, INOD

                ! Initialisation
                NLIST = 0
                INLIST = 0
                FINDELE = 0
                COLELE = 0

                ! NLIST is the number of elements each node belongs to...
                !  print *,'NONODS,totele,MXNCOLEL:',NONODS,totele,MXNCOLEL
                do ELE=1,TOTELE! Was loop
                    do ILOC=1,NLOC! Was loop
                        !      print *,'iloc,nloc,totele,ele:', iloc,nloc,totele,ele
                        !      print *,'NDGLNO((ELE-1)*NLOC+ILOC):',NDGLNO((ELE-1)*NLOC+ILOC)
                        INOD=NDGLNO((ELE-1)*NLOC+ILOC)
                        NLIST(INOD)=NLIST(INOD)+1
                    END DO
                END DO
                !  stop 771

                ! FINDELE is a pointer to the first element
                ! associated with a given global node (NOD)
                COUNT=0
                do NOD=1,NONODS! Was loop
                    FINDELE(NOD)=COUNT+1
                    COUNT=COUNT+NLIST(NOD)
                END DO
                FINDELE(NONODS+1)=COUNT+1
                NCOLEL=COUNT

                ! COLELE is a list of the element numbers each node belongs
                ! to stored in the order of the global nodes...
                ! INLIST is the element number the node belongs to.
                DO ELE=1,TOTELE! Was loop
                    DO ILOC=1,NLOC! Was loop
                        INOD=NDGLNO((ELE-1)*NLOC+ILOC)
                        INLIST(INOD)=INLIST(INOD)+1
                        IF (FINDELE(INOD)-1+INLIST(INOD).GT.MXNCOLEL) THEN
                            STOP 'COLELE ARRAY OUT OF BOUNDS--SUB:PHILNODELE'
                        ENDIF
                        COLELE(FINDELE(INOD)-1+INLIST(INOD))=ELE
                    END DO
                END DO
                RETURN

            end subroutine philnodele


            subroutine conv_quad_to_lin_tri_tet( ndgln_p2top1, nloc_lin, cv_nloc, sub_lin_totele )
                ! convert quadratic element into a series of linear elements...
                integer, intent( in ) :: nloc_lin, cv_nloc, sub_lin_totele
                integer, intent( inout ) :: ndgln_p2top1(sub_lin_totele*nloc_lin)
                ! local variables...
                integer :: sub_ele

                if(cv_nloc==6) then ! quadratic triangle...
                    sub_ele = 1
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 1
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 2
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 4

                    sub_ele = 2
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 2
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 4
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 5

                    sub_ele = 3
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 2
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 3
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 5

                    sub_ele = 4
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 4
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 5
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 6

                else if(cv_nloc==10) then ! quadratic triangle...

                    sub_ele = 1
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 7
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 8
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 9
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 10

                    sub_ele = 2
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 1
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 2
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 4
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 7

                    sub_ele = 3
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 2
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 7
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 8
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 4

                    sub_ele = 4
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 2
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 3
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 4
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 8

                    sub_ele = 5
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 3
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 5
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 4
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 8

                    sub_ele = 6
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 4
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 5
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 9
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 8

                    sub_ele = 7
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 5
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 6
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 4
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 9

                    sub_ele = 8
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 7
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 9
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 8
                    ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 4

                else
                    ewrite(3,*) 'not a viable option for calc_sub_lin_tri_tet'
                end if

                return

            end subroutine conv_quad_to_lin_tri_tet

    END SUBROUTINE CALC_ANISOTROP_LIM

    !
    !sprint_to_do!think about this
    Subroutine TRILOCCORDS(Xp,Yp,Zp, &
        N1, N2, N3, N4, &
        X1,Y1,Z1, &
        X2,Y2,Z2, &
        X3,Y3,Z3, &
        X4,Y4,Z4  )

        IMPLICIT NONE
        Real Xp, Yp, Zp

        Real N1, N2, N3, N4

        Real X1,Y1,Z1
        Real X2,Y2,Z2
        Real X3,Y3,Z3
        Real X4,Y4,Z4

        Real Volume

        !     calculate element volume...

        Volume = TetVolume(X1, Y1, Z1, &
            X2, Y2, Z2, &
            X3, Y3, Z3, &
            X4, Y4, Z4)

        Volume = Volume /6.0


        !     vol coords...

        N1 = TetVolume(Xp, Yp, Zp, &
            X2, Y2, Z2, &
            X3, Y3, Z3, &
            X4, Y4, Z4)

        N1 = N1/(6.0*Volume)



        N2 = TetVolume(X1, Y1, Z1, &
            Xp, Yp, Zp, &
            X3, Y3, Z3, &
            X4, Y4, Z4)

        N2 = N2/(6.0*Volume)



        N3 = TetVolume(X1, Y1, Z1, &
            X2, Y2, Z2, &
            Xp, Yp, Zp, &
            X4, Y4, Z4)

        N3 = N3/(6.0*Volume)


        N4 = TetVolume(X1, Y1, Z1, &
            X2, Y2, Z2, &
            X3, Y3, Z3, &
            Xp, Yp, Zp)

        N4 = N4/(6.0*Volume)


        Return

    end subroutine triloccords

    !
    !
    !
    !sprint_to_do!think about this
    Subroutine TRILOCCORDS2D(Xp,Yp, &
        &     N1, N2, N3,  &
        X1,Y1, &
        X2,Y2, &
        X3,Y3 )

        IMPLICIT NONE
        Real Xp,Yp, &
            N1, N2, N3,  &
            X1,Y1, &
            X2,Y2, &
            X3,Y3

        Real AREA

        AREA = TRIAREAF_SIGN( X1, Y1, X2, Y2, X3, Y3)
        !     area coords...

        N1 = TRIAREAF_SIGN(Xp, Yp,  &
            &     X2, Y2,  &
            &     X3, Y3 )

        N1 = N1/AREA



        N2 = TRIAREAF_SIGN(X1, Y1, &
            &     Xp, Yp,  &
            &     X3, Y3 )

        N2 = N2/AREA



        N3 = TRIAREAF_SIGN(X1, Y1,  &
            &     X2, Y2,  &
            &     Xp, Yp )

        N3 = N3/AREA


        Return

    end subroutine triloccords2d


    logical function shock_front_in_ele(ele, Mdims, sat, ndgln, Imble_frac)
        !Detects whether the element has a shockfront or not
        implicit none
        integer :: ele
        type(multi_dimensions), intent(in) :: Mdims
        real, dimension(:,:), intent(in) :: sat !(saturations of an element)
        real, dimension(:), intent(in) ::Imble_frac
        type(multi_ndgln), intent(in) :: ndgln
        !Local variables
        integer :: iphase, cv_iloc
        real :: minival, maxival, aux
        real, parameter :: tol = 0.05!Shock fronts smaller than this are unlikely to require extra handling

        minival = 10.; maxival = 0.
        do cv_iloc = 1, Mdims%cv_nloc
            do iphase = 1, mdims%nphase - 1
                aux = sat(iphase, ndgln%cv((ELE-1)*Mdims%cv_nloc+cv_iloc)) - Imble_frac(iphase)
                minival = min(aux, minival)
                maxival = max(aux, maxival)
            end do
        end do

        shock_front_in_ele = minival < tol .and. (maxival-minival) > tol

    end function shock_front_in_ele

end module cv_advection


