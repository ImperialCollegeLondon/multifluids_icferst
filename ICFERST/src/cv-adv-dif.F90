!    Copyright (C) 2020 Imperial College London and others.
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Affero General Public License
!    as published by the Free Software Foundation,
!    version 3.0 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without seven the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"
!> This module contains all the tools to assemble and solve the equations and fields associated with the CV mesh, i.e.
!> transport equation, continuity equation, DCVFE gradient matrix and the Laplacian system for the zeta-potential
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
    !>calculates derivatives of vector fields
    interface DG_DERIVS_ALL
        module procedure DG_DERIVS_ALL1
        module procedure DG_DERIVS_ALL2
    end interface DG_DERIVS_ALL

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

    !---------------------------------------------------------------------------
    !> @author Chris Pain, Pablo Salinas
    !> @brief This subroutines generates the transport equation for a cv field. It also can generate the Continuity equation if GETCT = .true.
    !> and also generate the gradient matrix of the momentum equation for the DCVFE method if the option is activated
    !>
    !>@param  state   Linked list containing all the fields defined in diamond and considered by Fluidity
    !>@param  packed_state  Linked list containing all the fields used by IC-FERST, memory partially shared with state
    !>@param  final_phase This is the final phase to be assembled, in this way we can assemble from phase 1 to final_phase not necessarily being for all the phases
    !>@param Mdims Data type storing all the dimensions describing the mesh, fields, nodes, etc
    !>@param  CV_GIdims Gauss integration numbers for CV fields
    !>@param  FE_GIdims Gauss integration numbers for FE fields
    !>@param  CV_funs Shape functions for the CV mesh
    !>@param  FE_funs Shape functions for the FE mesh
    !>@param  Mspars Sparsity of the matrices
    !>@param  ndgln Global to local variables
    !>@param  Mdisopt Discretisation options
    !>@param  Mmat Matrices for ICFERST
    !>@param  upwnd Sigmas to compute the fluxes at the interphase for porous media
    !>@param  tracer  Tracer considered for the transport equation
    !>@param  density  Density of the field
    !>@param  velocity  Velocity of the field
    !>@param  multi_absorp  Absoprtion of associated with the transport field
    !>@param  CV_DISOPT, CV_DG_VEL_INT_OPT, IGOT_THETA_FLUX. More Discretisation options
    !>@param  IGOT_T2. True (1) if solving for a tracer, false otherwise
    !>@param  DIAG_SCALE_PRES Diagonal scaling of (distributed) pressure matrix (used to treat pressure implicitly)
    !>@param  DIAG_SCALE_PRES_COUP  Diagonal scaling of (distributed) pressure matrix (for wells)
    !>@param  INV_B   Coupling term of the wells
    !>@param  MASS_MN_PRES ??
    !>@param  MASS_SUF ?? 
    !>@param  DEN_ALL  Density of the field, different memory to the input field density, used to apply the Boussinesq approximation
    !>@param  DENOLD_ALL  Density of the field, different memory to the input field density, used to apply the Boussinesq approximation
    !>@param  THETA_GDIFF ! (nphase,Mdims%cv_nonods)
    !>@param  THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
    !>@param  DT, CV_THETA, CV_BETA: Time step size and discretisation options
    !>@param  SUF_SIG_DIAGTEN_BC. Like upwnd but for the boundary
    !>@param  DERIV   Derivative of the density against the pressure (nphase,Mdims%cv_nonods)
    !>@param  CV_P ! Control volume pressure, useless for the DCVFEM (1,Mdims%npres,Mdims%cv_nonods)
    !>@param  SOURCT_ALL  Source term of the tracer equation
    !>@param  ABSORBT_ALL Absorption to be used here
    !>@param  VOLFRA_PORE   Porosity field (Mdims%npres,Mdims%totele)
    !>@param  GETCV_DISC  obtain the transport equation
    !>@param GETCT obtain the continuity equation
    !>@param GET_THETA_FLUX, USE_THETA_FLUX
    !>@param THERMAL true if solving for heat transport 
    !>@param  MEAN_PORE_CV   Porosity defined control volume wise
    !>@param  MASS_ELE_TRANSP Mass of the elements
    !>@param  saturation PhaseVolumeFraction field
    !>@param  TDIFFUSION  Diffusion associated with the tracer field
    !>@param  Phase_with_Pc   Field that defines the capillary pressure, i.e. non-wetting phase
    !>@param  VAD_parameter  Vanishing artificial diffusion parameter
    !>@param  Courant_number  Obvious ins't it?
    !>@param  Permeability_tensor_field  Permeability field
    !>@param  calculate_mass_delta  Variable used to control the mass conservation of the system
    !>@param  eles_with_pipe  Elements that have a pipe
    !>@param  pipes_aux  Information required to define wells
    !>@param  porous_heat_coef, porous_heat_coef_old includes an average of porous and fluid heat coefficients 
    !>@param  solving_compositional Solving for composition if true
    !>@param assemble_collapsed_to_one_phase Collapses phases and solves for one single temperature. When there is thermal equilibrium
    !>@param  outfluxes  Contains all the fields required to compute the outfluxes of the model and create the outfluxes.csv file. Computed when assembling the continuity equation
    !>@param  nonlinear_iteration Current non-linear iteration
  !---------------------------------------------------------------------------
    SUBROUTINE CV_ASSEMB( state, packed_state, &
          final_phase, Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd, &
          tracer, velocity, density, multi_absorp, &
          DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B,&
          DEN_ALL, DENOLD_ALL, &
          CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, &
          SUF_SIG_DIAGTEN_BC, &
          DERIV, CV_P, &
          SOURCT_ALL, ABSORBT_ALL, VOLFRA_PORE, &
          GETCV_DISC, GETCT, &
          IGOT_T2, IGOT_THETA_FLUX, GET_THETA_FLUX, USE_THETA_FLUX, &
          THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
          MEAN_PORE_CV, &
          MASS_MN_PRES, THERMAL, MASS_SUF, MASS_ELE_TRANSP, &
          TDIFFUSION, &
          saturation, VAD_parameter, Phase_with_Pc, Courant_number,&
          Permeability_tensor_field, calculate_mass_delta, eles_with_pipe, pipes_aux, &
          porous_heat_coef,porous_heat_coef_old, outfluxes, solving_compositional, nonlinear_iteration,&
          assemble_collapsed_to_one_phase)
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
          integer, intent(in) ::  final_phase
          type(multi_dimensions), intent(in) :: Mdims
          type(multi_GI_dimensions), intent(in) :: CV_GIdims
          type(multi_GI_dimensions) :: FE_GIdims
          type(multi_shape_funs), intent(inout) :: CV_funs
          type(multi_shape_funs) :: FE_funs
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
          REAL, DIMENSION( :, : ), intent( inout ), allocatable :: DIAG_SCALE_PRES
          REAL, DIMENSION( :, :, : ), intent( inout ), allocatable :: DIAG_SCALE_PRES_COUP 
          REAL, DIMENSION( :, :, : ), intent( inout ), allocatable :: INV_B 
          REAL, DIMENSION( : ), intent( inout ) :: MASS_MN_PRES
          REAL, DIMENSION( : ), intent( inout ) :: MASS_SUF
          REAL, DIMENSION( :, : ), target, intent( inout ) :: DEN_ALL
          REAL, DIMENSION( :, : ), intent( inout ), target :: DENOLD_ALL
          REAL, DIMENSION( :, : ), intent( inout ) :: THETA_GDIFF ! (nphase,Mdims%cv_nonods)
          REAL, DIMENSION( :, : ), intent( inout ), optional :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
          REAL, intent( in ) :: DT, CV_THETA, CV_BETA
          REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
          REAL, DIMENSION( :, : ), intent( in ) :: DERIV 
          REAL, DIMENSION( :, :, : ), intent( in ) :: CV_P ! (1,Mdims%npres,Mdims%cv_nonods)
          REAL, DIMENSION( :, : ), intent( in) :: SOURCT_ALL
          REAL, DIMENSION( :, :, : ), pointer, intent( in ) :: ABSORBT_ALL
          REAL, DIMENSION( :, : ), intent( in ) :: VOLFRA_PORE 
          LOGICAL, intent( in ) :: GETCV_DISC, GETCT, GET_THETA_FLUX, USE_THETA_FLUX, THERMAL
          REAL, DIMENSION( :, : ), intent( inout ) :: MEAN_PORE_CV 
          REAL, DIMENSION( : ), intent( inout ), OPTIONAL  :: MASS_ELE_TRANSP
          type(tensor_field), intent(in), optional, target :: saturation
          REAL, DIMENSION( :, :, :, : ), intent( in ), optional :: TDIFFUSION
          !Variables for Vanishing artificial diffusion
          integer, optional, intent(in) :: Phase_with_Pc 
          real, optional, dimension(:), intent(in) :: VAD_parameter
          !Variables to cache get_int_vel OLD
          real, optional, dimension(:), intent(inout) :: Courant_number
          type( tensor_field ), optional, pointer, intent(in) :: Permeability_tensor_field
          ! Calculate_mass variable
          real, dimension(:,:), optional :: calculate_mass_delta
          type(pipe_coords), dimension(:), optional, intent(in):: eles_with_pipe
          type (multi_pipe_package), intent(in) :: pipes_aux
          REAL, DIMENSION( :), optional, intent(in) :: porous_heat_coef, porous_heat_coef_old
          logical, optional, intent(in) :: solving_compositional, assemble_collapsed_to_one_phase
          ! Variable to store outfluxes
          type (multi_outfluxes), optional, intent(inout) :: outfluxes
          !Non-linear iteration count
          integer, optional, intent(in) :: nonlinear_iteration
          ! ###################Local variables############################
          REAL :: ZERO_OR_TWO_THIRDS
          REAL, dimension(:), ALLOCATABLE :: FEMFGI, DIFF_COEF, COEF

          !        ===>  LOGICALS  <===
          LOGICAL :: GETMAT, D1, D3, GOT_DIFFUS, INTEGRAT_AT_GI, GET_GTHETA, QUAD_OVER_WHOLE_ELE
          logical :: skip, GOT_T2, use_volume_frac_T2, logical_igot_theta_flux, zero_vel_BC
          ! If GET_C_IN_CV_ADVDIF_AND_CALC_C_CV then form the Mmat%C matrix in here also based on control-volume pressure.
          logical :: GET_C_IN_CV_ADVDIF_AND_CALC_C_CV
        !   LOGICAL :: DISTCONTINUOUS_METHOD
          !Logical to check if we using a conservative method or not, to save cpu time
          logical :: conservative_advection
          !        ===> GENEREIC INTEGERS <===
          INTEGER :: COUNT, ICOUNT, JCOUNT, ELE, ELE2, GI, GCOUNT, SELE, V_SILOC, U_KLOC, CV_ILOC, CV_JLOC, IPHASE, JPHASE, &
              CV_NODJ, CV_NODI, U_NODK, TIMOPT, X_NODI,  X_NODJ, CV_INOD, MAT_NODI,  MAT_NODJ, CV_SILOC
          INTEGER :: I, IDIM, U_ILOC, ELE3, k, CV_KLOC, CV_NODK, IFI, CV_KLOC2,CV_NODK2,CV_SKLOC, iofluxes,&
              U_KLOC2,U_NODK2,U_SKLOC, JDIM, IGETCT, global_face,J, nb, i_use_volume_frac_t2
          integer, dimension(:), pointer :: neighbours
          INTEGER, dimension(final_phase) :: LOC_WIC_T_BC_ALL
          !        ===>  GENERIC REALS  <===
          REAL :: HDC, RSUM, W_SUM_ONE1, W_SUM_ONE2, one_m_cv_beta, auxR, NDOTQ_HAT
          REAL, dimension(final_phase) :: ROBIN1, ROBIN2, BCZERO
          !Local copy of tracers and densities
          real, dimension(:), pointer :: LOC_T_I, LOC_DEN_I, LOC_TOLD_I, LOC_DENOLD_I, LOC_T2_I, LOC_T2OLD_I, LOC_T_J,LOC_T2_J, AUX_T
          REAL :: GRAVTY
          LOGICAL, DIMENSION( Mdims%x_nonods ) :: X_SHARE
          integer, dimension (Mdims%cv_nloc) ::CV_OTHER_LOC, SHAPE_CV_SNL
          integer, dimension (Mdims%cv_snloc) :: CV_SLOC2LOC
          integer, dimension (Mdims%u_snloc) :: U_SLOC2LOC
          integer, dimension (Mdims%mat_nloc) ::MAT_OTHER_LOC
          !Allocate only for get_CT
          INTEGER, DIMENSION( : ), allocatable :: JCOUNT_KLOC, JCOUNT_KLOC2, ICOUNT_KLOC, ICOUNT_KLOC2, &
              C_JCOUNT_KLOC, C_JCOUNT_KLOC2, C_ICOUNT_KLOC, C_ICOUNT_KLOC2

          REAL, DIMENSION( : ), allocatable ::  N
          real, dimension (Mdims%cv_nonods) ::  SUM_CV
          real, dimension (:), pointer :: MASS_ELE
          type(vector_field), pointer :: vfield
          REAL, DIMENSION( Mdims%ndim, CV_GIdims%scvngi ) :: CVNORMX_ALL
          REAL, DIMENSION( :,: ), pointer :: XC_CV_ALL
          REAL, DIMENSION( Mdims%ndim,final_phase,Mdims%u_nloc ) :: UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL
          REAL, DIMENSION( :, : ), allocatable :: CAP_DIFFUSION
          !###Variables for shape function calculation###
          type (multi_dev_shape_funs) :: SdevFuns
          type (multi_dev_shape_funs) :: FSdevFuns
          ! Variables used in GET_INT_VEL_NEW:
          REAL, DIMENSION ( Mdims%ndim, final_phase, Mdims%u_nloc ) :: LOC_U, LOC2_U, LOC_NU, LOC2_NU
          REAL, DIMENSION ( Mdims%ndim, final_phase, Mdims%u_snloc ) :: SLOC_NU
          INTEGER :: CV_KNOD2, U_SNODK
          ! nphase Variables:
          real, dimension(final_phase)::NDOTQ, INCOME, CAP_DIFF_COEF_DIVDX, DIFF_COEF_DIVDX, NDOTQNEW, &
              LIMT, LIMT_HAT, LIMD, LIMDT, LIMT2
          REAL, DIMENSION( Mdims%ndim, final_phase, Mdims%cv_nloc, Mdims%totele ) :: DTX_ELE_ALL
          REAL , DIMENSION( Mdims%ndim, final_phase ) :: NUGI_ALL
          LOGICAL :: integrate_other_side_and_not_boundary
          !Working variables
          real, dimension(:,:), pointer :: T_ALL, TOLD_ALL, T2_ALL, T2OLD_ALL, X_ALL
          real, dimension(:, :, :), pointer :: U_ALL, NU_ALL, NUOLD_ALL

          real, dimension( final_phase ) :: DIAG_SCALE_PRES_phase
          real, dimension( final_phase ) :: ct_rhs_phase_cv_nodi, ct_rhs_phase_cv_nodj
          real, dimension(Mdims%npres) :: R_PRES
          real, dimension( final_phase ) :: R_PHASE, ct_rhs_phase

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
          real, dimension( Mdims%ndim,final_phase ) :: rdum_ndim_nphase_1
          real, dimension( final_phase ) :: rdum_nphase_1, rdum_nphase_2, rdum_nphase_3
          real, dimension( final_phase ) :: THETA_VEL
          real, dimension( final_phase) :: LOC_CV_RHS_I, LOC_CV_RHS_J, LOC_MAT_II, LOC_MAT_JJ, LOC_MAT_IJ, LOC_MAT_JI
          REAL, DIMENSION ( Mdims%ndim,final_phase ) :: UDGI_ALL, UDGI2_ALL, ROW_SUM_INV_VI, ROW_SUM_INV_VJ, UDGI_ALL_FOR_INV
          type( vector_field ), pointer :: MeanPoreCV
          !! femdem
          type( tensor_field_pointer ), dimension(4+2*IGOT_T2) :: psi,fempsi
          type( vector_field_pointer ), dimension(1) :: PSI_AVE,PSI_INT
          type( tensor_field ), pointer :: old_tracer
          type (tensor_field_pointer), allocatable, DIMENSION(:) :: outfluxes_fields

          ! variables for pipes (that are needed in cv_assemb as well), allocatable because they are big and barely used
          Real, dimension(:), pointer :: MASS_CV
          !Permeability and immobile fractions
          type( tensor_field ), pointer :: perm
          real, dimension( : , : ), pointer ::Imble_frac
          !Variables for Vanishing artificial diffusion (VAD)
          logical :: VAD_activated, between_elements, on_domain_boundary, flux_limited_vad
          !Variable to decide if we are introducing the sum of phases = 1 in Ct or elsewhere
          logical :: Solve_all_phases
          !Variables for get_int_vel_porous_vel
          real, dimension(final_phase):: rsum_nodi, rsum_nodj
          integer :: stat, ipres
          !Variables for assemble_collapsed_to_one_phase; Note that diffusion parameters etc need
          logical :: loc_assemble_collapsed_to_one_phase !to be then scaled before getting into this subroutine
          integer :: assembly_phase
          !Variables to calculate flux across boundaries
          logical :: calculate_flux
          integer :: U_JLOC

          logical :: have_absorption
          !Variables for get int_vel_porous_vel
          logical :: iv_Incomming_flow
          REAL, DIMENSION(Mdims%ndim) :: iv_SUF_SIG_DIAGTEN_BC_GI
          INTEGER :: iv_u_kloc, iv_u_skloc, iv_cv_kloc, iv_idim, iv_CV_SKLOC, iv_CV_SNODK, iv_CV_SNODK_IPHA, iv_IPHASE, iv_u_kloc2
          real, dimension(final_phase) :: iv_aux_tensor, iv_sigma_aver, iv_aux_tensor2
          ! ####Variables for outfluxes#####
          logical :: compute_outfluxes
          real, dimension(:, :,:), allocatable :: bcs_outfluxes!the total mass entering the domain is captured by 'bcs_outfluxes'
          real, allocatable, dimension(:) :: calculate_mass_internal  ! internal changes in mass will be captured by 'calculate_mass_internal'
          real :: tmp1, tmp2, tmp3  ! Variables for parallel mass calculations


          !Decide if we are solving for nphases-1
          Solve_all_phases = .not. have_option("/numerical_methods/solve_nphases_minus_one")
          !Check vanishing artificial diffusion options
          VAD_activated = .false.
          if (present(VAD_parameter) .and. present(Phase_with_Pc)) then
              VAD_activated = Phase_with_Pc >0
          end if
          flux_limited_vad = have_option("/numerical_methods/flux_limited_vad")
          logical_igot_theta_flux = IGOT_THETA_FLUX == 1


          loc_assemble_collapsed_to_one_phase = .false.
          if (present_and_true(assemble_collapsed_to_one_phase)) then
            loc_assemble_collapsed_to_one_phase = .true.
          end if
          have_absorption = associated( absorbt_all )

          !Check pressure matrix based on Control Volumes
          !If we do not have an index where we have stored Mmat%C_CV, then we need to calculate it
          if (Mmat%CV_pressure) then
              GET_C_IN_CV_ADVDIF_AND_CALC_C_CV = .not.Mmat%stored !.true.
          else
              GET_C_IN_CV_ADVDIF_AND_CALC_C_CV = .false.
          end if

          call get_option( "/physical_parameters/gravity/magnitude", gravty, stat )

          !#################SET WORKING VARIABLES#################

          call get_var_from_packed_state(packed_state,PressureCoordinate = X_ALL,&
             OldNonlinearVelocity = NUOLD_ALL, NonlinearVelocity = NU_ALL)
          if (.not. present_and_true(solving_compositional)) then
            call get_var_from_packed_state(packed_state, CV_Immobile_Fraction = Imble_frac)
          end if
          !For every Field_selector value but 3 (saturation) we need U_ALL to be NU_ALL
          U_ALL => NU_ALL

          old_tracer=>extract_tensor_field(packed_state,GetOldName(tracer))
          T_ALL =>tracer%val(1,:,:)
          TOLD_ALL =>old_tracer%val(1,:,:)
          !Initialise local t2
          LIMT2 =1.0
          if (tracer%name == "PackedPhaseVolumeFraction") call get_var_from_packed_state(packed_state,Velocity = U_ALL)
         !################## END OF SET VARIABLES ##################
          compute_outfluxes = present(calculate_mass_delta) .and. present(outfluxes) .and. GETCT
          if (compute_outfluxes) then
              ! Initialise the calculate_mass variables
              !Allocate array to pass to store mass going through the boundaries
              if (allocated( outfluxes%outlet_id )) then
                  allocate(bcs_outfluxes(Mdims%nphase, Mdims%cv_nonods, 0:size(outfluxes%outlet_id))); bcs_outfluxes= 0.!position zero is to store outfluxes over all bcs
                else
                  allocate(bcs_outfluxes(Mdims%nphase, Mdims%cv_nonods, 0:1)); bcs_outfluxes= 0.!position zero is to store outfluxes over all bcs
              end if
              allocate ( calculate_mass_internal(final_phase));calculate_mass_internal(:) = 0.0  ! calculate_internal_mass subroutine
              allocate(outfluxes_fields(size(outfluxes%field_names,2)))
              !Use this as flag to know whether we are outputting the csv file or not here
              if (allocated(outfluxes%area_outlet)) then
                !Initialize to zero the area of the outlet surfaces
                 outfluxes%area_outlet = 0.; outfluxes%totout = 0.
                !Extract fields to be used for the outfluxes section
                do k = 1, size(outfluxes%field_names,2)!We use the first phase as it is the one that must contain all the prognostic fields
                    outfluxes_fields(k)%ptr =>extract_tensor_field( packed_state, "Packed"//outfluxes%field_names(1,k) )
                    if (outfluxes%calculate_flux)outfluxes%avgout(k, :,:) = 0.0
                end do
              end if
          end if

          !! Get boundary conditions from field
          call get_entire_boundary_condition(tracer,['weakdirichlet','robin        '],tracer_BCs,WIC_T_BC_ALL,boundary_second_value=tracer_BCs_robin2)
          call get_entire_boundary_condition(density,['weakdirichlet'],density_BCs,WIC_D_BC_ALL)
          if (present(saturation))&
              call get_entire_boundary_condition(saturation,['weakdirichlet','robin        '],saturation_BCs,WIC_T2_BC_ALL,boundary_second_value=saturation_BCs_robin2)
          call get_entire_boundary_condition(velocity,['weakdirichlet'],velocity_BCs,WIC_U_BC_ALL)
          pressure => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedFEPressure" )
          call get_entire_boundary_condition(pressure,['weakdirichlet','freesurface  '],pressure_BCs,WIC_P_BC_ALL)

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
           if (tracer%name == "PackedTemperature" ) then
              !(tracer%name == "PackedEnthalpy")&
              ! .or. tracer%name == "PackedConcentration")  then !Not sure if it is required for temperature either...
              allocate( suf_t_bc( 1,size(tracer_BCs%val,2),Mdims%cv_snloc*Mdims%stotel ), suf_t_bc_rob1( 1,size(tracer_BCs%val,2),Mdims%cv_snloc*Mdims%stotel ), &
                  suf_t_bc_rob2( 1,size(tracer_BCs%val,2),Mdims%cv_snloc*Mdims%stotel ) )
              call update_boundary_conditions( state, Mdims%stotel, Mdims%cv_snloc, size(tracer_BCs%val,2), &!TEMPORARY, FIXME! sprint_to_do is this call needed?
                  suf_t_bc, suf_t_bc_rob1, suf_t_bc_rob2, tracer)                                                  !BCs are updated autoamtically
              SUF_T_BC_ALL=>suf_t_bc
              SUF_T_BC_ROB1_ALL=>suf_t_bc_rob1
              SUF_T_BC_ROB2_ALL=>suf_t_bc_rob2
          end if

          ewrite(3,*) 'In CV_ASSEMB'
          GOT_DIFFUS = .false.; 
          if (present(TDIFFUSION)) then
              GOT_DIFFUS = ( R2NORM( TDIFFUSION, size(TDIFFUSION,1) * size(TDIFFUSION,2) * size(TDIFFUSION,3) * size(TDIFFUSION,4)  ) /= 0 )!<=I hate this thing...
              call allor(GOT_DIFFUS)                                                  !it should be if present then true, but it breaks the parallel CWC P1DGP2
          end if
          call get_option( "/material_phase[0]/phase_properties/Viscosity/viscosity_scheme/zero_or_two_thirds", zero_or_two_thirds, default=2./3. )
          ewrite(3,*)'CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, GOT_DIFFUS:', &
              CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, GOT_DIFFUS
          ewrite(3,*)'GETCV_DISC, GETCT', GETCV_DISC, GETCT

          !cv_beta == 1 means conservative, meaning that everything multiplied by one_m_cv_beta can be ignored
          one_m_cv_beta = 1.0 - cv_beta
          conservative_advection = abs(one_m_cv_beta) <= 1e-8
          QUAD_OVER_WHOLE_ELE=.FALSE.
          ! Allocate memory for the control volume surface shape functions, etc.
          IF(GETCT) THEN
              ALLOCATE( JCOUNT_KLOC( Mdims%u_nloc ));ALLOCATE( JCOUNT_KLOC2( Mdims%u_nloc ))
              ALLOCATE( ICOUNT_KLOC( Mdims%u_nloc ));ALLOCATE( ICOUNT_KLOC2( Mdims%u_nloc ))
              IF(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV) THEN
                  ALLOCATE( C_JCOUNT_KLOC( Mdims%u_nloc )); ALLOCATE( C_JCOUNT_KLOC2( Mdims%u_nloc ))
                  ALLOCATE( C_ICOUNT_KLOC( Mdims%u_nloc )); ALLOCATE( C_ICOUNT_KLOC2( Mdims%u_nloc ))
              ENDIF
          ENDIF
        !   DISTCONTINUOUS_METHOD = ( Mdims%cv_nonods == Mdims%totele * Mdims%cv_nloc )
          !Pointer to permeability
        if (present(Permeability_tensor_field)) then
            perm => Permeability_tensor_field
        else
        perm=>extract_tensor_field(packed_state,"Permeability")
        end if
          !Initialize Courant number for porous media
          if (present(Courant_number)) Courant_number = 0.
          CVNORMX_ALL=0.0

          ! The procity mapped to the CV nodes
          D1 = ( Mdims%ndim == 1 )
          D3 = ( Mdims%ndim == 3 )
          GETMAT = .TRUE.

          X_SHARE = .FALSE.
          ! Determine FEMT (finite element wise) etc from T (control volume wise)
          ! Also determine the CV mass matrix MASS_CV_PLUS and centre of the CV's XC_CV_ALL
          ! This is for projecting to finite element basis functions...
          IGETCT = 0
          IF ( GETCT ) IGETCT = 1
          GOT_T2=( IGOT_T2 == 1 )
          use_volume_frac_T2=( IGOT_T2 == 1 .or. thermal )
          i_use_volume_frac_t2= 0
          if (use_volume_frac_T2) i_use_volume_frac_t2= 1


          IF ( GOT_T2 .OR. THERMAL) call get_var_from_packed_state( packed_state, &
              PhaseVolumeFraction = T2_ALL, OldPhaseVolumeFraction = T2OLD_ALL )

          !###############Conditional allocations######################
          LIMT_HAT=0.0
          IF ( VAD_activated) THEN

              ALLOCATE( CAP_DIFFUSION( final_phase, Mdims%mat_nonods ) )
              !Introduce the information in CAP_DIFFUSION
              CAP_DIFFUSION = 0.!Initialize to zero just in case
              do ele = 1, Mdims%totele
                  do CV_ILOC = 1, Mdims%cv_nloc
                      CV_NODI = ndgln%cv(CV_ILOC + (ele-1) * Mdims%cv_nloc)
                      MAT_NODI = ndgln%mat(CV_ILOC + (ele-1) * Mdims%cv_nloc)
                      CAP_DIFFUSION(Phase_with_Pc, MAT_NODI) = &
                          - T_ALL(Phase_with_Pc, CV_NODI) * VAD_parameter(CV_NODI)
                  end do
              end do
          ENDIF

          ndotq = 0.

        psi_int(1)%ptr=>extract_vector_field(packed_state,"CVIntegral")
        psi_ave(1)%ptr=>extract_vector_field(packed_state,"CVBarycentre")
        vfield => extract_vector_field(packed_state,"MASS_ELE")
        MASS_ELE => vfield%val(1,:)
        XC_CV_ALL => psi_ave(1)%ptr%val
        MASS_CV   => psi_int(1)%ptr%val(1,:)
        !This works because GETCT is the first call and therefore we will have later on masses and barycenters
        if (.not.associated(CV_funs%CV2FE%refcount) ) then!This is true after adapt and at the beginning
        psi(1)%ptr=>tracer
        psi(2)%ptr=>old_tracer
        call PROJ_CV_TO_FEM(packed_state, &!For porous media we are just pointing memory from PSI to FEMPSI
            FEMPSI(1:2),PSI(1:2), &!we need to get rid of all of this... check if for inertia is there any gain at all
            Mdims, CV_GIdims, CV_funs, Mspars, ndgln, &
            IGETCT, X_ALL, MASS_ELE, MASS_MN_PRES, &
            tracer,PSI_AVE, PSI_INT)
        end if

        !Store mass_CV in packed_state. Ideally we would do this somewhere else, but here we are...
          IF (PRESENT(MASS_ELE_TRANSP)) &
              MASS_ELE_TRANSP = MASS_ELE

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
        MeanPoreCV=>extract_vector_field(packed_state,"MeanPoreCV")
        MeanPoreCV%val=MEAN_PORE_CV
        
        !Obtain elements surrounding an element (FACE_ELE) only if it is not stored yet
        if (.not. associated(Mmat%FACE_ELE)) then 
        allocate(Mmat%FACE_ELE(CV_GIdims%nface, Mdims%totele))
        Mmat%FACE_ELE = 0.
        CALL CALC_FACE_ELE( Mmat%FACE_ELE, Mdims%totele, Mdims%stotel, CV_GIdims%nface, &
        Mspars%ELE%fin, Mspars%ELE%col, Mdims%cv_nloc, Mdims%cv_snloc, Mdims%cv_nonods, ndgln%cv, ndgln%suf_cv, &
        CV_funs%cv_sloclist, Mdims%x_nloc, ndgln%x )
        end if
        IF ( GOT_DIFFUS) THEN
            CALL DG_DERIVS_ALL( T_ALL, &
            DTX_ELE_ALL, &
            Mdims%ndim, final_phase, Mdims%totele, ndgln%cv, &
            ndgln%x, Mdims%x_nloc, ndgln%x,&
            CV_GIdims%cv_ngi, Mdims%cv_nloc, CV_funs%CVWEIGHT, &
            CV_funs%CVFEN, CV_funs%CVFENLX_ALL(1,:,:), CV_funs%CVFENLX_ALL(2,:,:), CV_funs%CVFENLX_ALL(3,:,:), &
            CV_funs%CVFEN, CV_funs%CVFENLX_ALL(1,:,:), CV_funs%CVFENLX_ALL(2,:,:), CV_funs%CVFENLX_ALL(3,:,:), &
            Mdims%x_nonods, X_ALL(1,:),X_ALL(2,:),X_ALL(3,:), &
            CV_GIdims%nface, Mmat%FACE_ELE, CV_funs%cv_sloclist, CV_funs%cv_sloclist, Mdims%cv_snloc, Mdims%cv_snloc, WIC_T_BC_ALL, SUF_T_BC_ALL, &
            CV_GIdims%sbcvngi, CV_funs%sbcvfen, CV_funs%sbcvfeweigh, &
            CV_funs%sbcvfen, CV_funs%sbcvfenslx, CV_funs%sbcvfensly )
        end if
          TIMOPT = MOD( CV_DISOPT, 2 )
          IF ( GETCT ) THEN ! Obtain the CV discretised Mmat%CT eqns plus RHS
              call zero(Mmat%CT_RHS)
              Mmat%CT = 0.0
              if ( Mmat%CV_pressure ) MASS_SUF=0.0
          END IF
          IF ( GETCV_DISC ) THEN ! Obtain the CV discretised advection/diffusion eqns
              call zero(Mmat%CV_RHS)
              IF ( GETMAT ) THEN
                  call zero(Mmat%petsc_ACV)
              END IF
          END IF
          GET_GTHETA = .FALSE.
          IF ( logical_igot_theta_flux ) THEN
              IF ( GET_THETA_FLUX ) THEN
                  THETA_FLUX = 0.0
                  ONE_M_THETA_FLUX = 0.0
                  THETA_FLUX_J = 0.0
                  ONE_M_THETA_FLUX_J = 0.0
                  !            IF ( IGOT_T2 == 1 ) THEN
                  IF ( GOT_T2 ) THEN
                      GET_GTHETA = .TRUE.
                      THETA_GDIFF = 0.0
                  END IF
              END IF
          END IF
          GLOBAL_FACE = 0

          if ( compute_outfluxes) then
              !Initialise mass conservation check; calculation of porevolume
              call mass_conservation_check_and_outfluxes(calculate_mass_delta, outfluxes, DEN_ALL, 1)
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

              !Local copies of U IS THIS NEEDED?????
              DO U_KLOC = 1, Mdims%u_nloc
                U_NODK = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_KLOC )
                LOC_U( :, :, U_KLOC)=U_ALL( :, 1:final_phase, U_NODK)
                LOC_NU( :, :, U_KLOC)=NU_ALL( :, 1:final_phase, U_NODK)
              END DO
              Loop_CV_ILOC: DO CV_ILOC = 1, Mdims%cv_nloc ! Loop over the nodes of the element

                  ! Global node number of the local node
                  CV_NODI = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
                  X_NODI = ndgln%x( ( ELE - 1 ) * Mdims%x_nloc  + CV_ILOC )
                  MAT_NODI = ndgln%mat( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )

                  ! Generate some local variables to reduce slicing (sprint_to_do this is a waste...)
                  LOC_T_I => T_ALL(1:final_phase, cv_nodi); LOC_TOLD_I => TOLD_ALL(1:final_phase, cv_nodi)
                  LOC_DEN_I => DEN_ALL(1:final_phase, cv_nodi); LOC_DENOLD_I => DENOLD_ALL(1:final_phase, cv_nodi)
                  if (use_volume_frac_T2) then
                    LOC_T2_I => T2_ALL(1:final_phase, cv_nodi); LOC_T2OLD_I => T2OLD_ALL(1:final_phase, cv_nodi)
                  end if
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
                      Conditional_CheckingNeighbourhood: IF ( CV_JLOC == -1 ) THEN
                          ! We are on the boundary or next to another element.  Determine CV_OTHER_LOC
                          ! CV_funs%cvfem_on_face(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
                          ! Look for these nodes on the other elements.
                          CALL FIND_OTHER_SIDE( CV_OTHER_LOC, Mdims%cv_nloc, Mdims%u_nloc, &
                              MAT_OTHER_LOC, INTEGRAT_AT_GI, &
                              Mdims%x_nloc, Mdims%xu_nloc, ndgln%x, ndgln%xu, &
                              Mdims%cv_snloc, CV_funs%cvfem_on_face( :, GI ), X_SHARE, ELE, ELE2,  &
                              Mspars%ELE%fin, Mspars%ELE%col )

                          IF ( INTEGRAT_AT_GI ) THEN
                              CV_JLOC = CV_OTHER_LOC( CV_ILOC )
                              SELE = 0
                              ELE3=0
                              IF ( CV_JLOC == 0 ) THEN ! We are on the boundary of the domain or subdomain
                                  CV_JLOC = CV_ILOC
                                  ! Calculate SELE, CV_SILOC, U_SLOC2LOC, CV_SLOC2LOC
                                  CALL CALC_SELE( ELE, ELE3, SELE, CV_SILOC, CV_ILOC, U_SLOC2LOC, CV_SLOC2LOC, &
                                      Mmat%FACE_ELE, GI, CV_funs, Mdims, CV_GIdims,&
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
                        !   IF ( between_elements ) THEN
                        !       CV_NODJ = ndgln%cv( ( ELE2 - 1 ) * Mdims%cv_nloc + CV_JLOC )
                        !       X_NODJ = ndgln%x( ( ELE2 - 1 ) * Mdims%cv_nloc + CV_JLOC )
                        !       MAT_NODJ = ndgln%mat( ( ELE2 - 1 ) * Mdims%cv_nloc + CV_JLOC )
                        !   ELSE
                              CV_NODJ = ndgln%cv( ( ELE - 1 )  * Mdims%cv_nloc + CV_JLOC )
                              X_NODJ = ndgln%x( ( ELE - 1 )  * Mdims%cv_nloc + CV_JLOC )
                              MAT_NODJ = ndgln%mat( ( ELE - 1 )  * Mdims%cv_nloc + CV_JLOC )
                        !   END IF

                        !   permeability_jump = .false.
                        !   !For the discontinuous formulation we want to use the discontinuous method only where there is a permeability jump
                        !   !elsewhere the normal method should work better as it is also more reliable
                        !   if (between_elements) then
                        !     permeability_jump = abs(perm%val(1,1,ele) - perm%val(1,1,ele2)/perm%val(1,1,ele)) > 1e-8
                        !   end if
                          !Create local variables to reduce slicing
                          LOC_T_J => T_ALL(1:final_phase, cv_nodj)
                          if (use_volume_frac_T2) then
                            LOC_T2_J => T2_ALL(1:final_phase, cv_nodj)
                          end if

                          if(CV_NODJ >= CV_NODI) then
                              ! this is for DG and boundaries of the domain
                            !   IF( between_elements ) THEN ! this is for DG
                            !       ! Calculate U_SLOC2LOC, CV_SLOC2LOC:
                            !       CV_SKLOC=0
                            !       DO CV_KLOC=1,Mdims%cv_nloc
                            !           CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                            !           IF(CV_KLOC2.NE.0) THEN
                            !               CV_SKLOC=CV_SKLOC+1
                            !               CV_SLOC2LOC(CV_SKLOC)=CV_KLOC
                            !               SHAPE_CV_SNL(CV_SKLOC) = CV_funs%scvfen(CV_KLOC,GI)
                            !           ENDIF
                            !       END DO
                            !       U_SKLOC=0
                            !       DO U_KLOC=1,Mdims%u_nloc
                            !           U_KLOC2 = U_OTHER_LOC( U_KLOC )
                            !           IF(U_KLOC2.NE.0) THEN
                            !               U_SKLOC=U_SKLOC+1
                            !               U_SLOC2LOC(U_SKLOC)=U_KLOC
                            !           ENDIF
                            !       END DO
                            !   ENDIF ! ENDOF IF( between_elements ) THEN
                              integrate_other_side_and_not_boundary = SELE.LE.0
                              GLOBAL_FACE = GLOBAL_FACE + 1
                              ! Calculate the control volume normals at the Gauss pts. Internal subroutine for speed
                              CALL SCVDETNX( Mdims, ndgln, X_ALL, CV_funs, CV_GIdims, on_domain_boundary, between_elements, &
                                    ELE, GI, SdevFuns%DETWEI, CVNORMX_ALL,XC_CV_ALL( :, CV_NODI ), X_NODI, X_NODJ)
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
                               ! Generate some local F variables ***************
                            !   IF (between_elements) THEN
                            !       LOC2_U = 0.
                            !       LOC2_NU = 0.
                            !       LOC2_NUOLD = 0.
                            !       DO U_SKLOC = 1, Mdims%u_snloc
                            !           U_KLOC = U_SLOC2LOC(U_SKLOC)
                            !           U_KLOC2 = U_OTHER_LOC( U_KLOC )
                            !           U_NODK2 = ndgln%u((ELE2-1)*Mdims%u_nloc+U_KLOC2)
                            !           LOC2_U(:, :, U_KLOC) = U_ALL(:, :, U_NODK2)
                            !           LOC2_NU(:, :, U_KLOC) = NU_ALL(:, :, U_NODK2)
                            !           LOC2_NUOLD(:, :, U_KLOC) = NUOLD_ALL(:, :, U_NODK2)
                            !       END DO
                            !       DO CV_SKLOC=1,Mdims%cv_snloc
                            !           CV_KLOC=CV_SLOC2LOC( CV_SKLOC )
                            !           CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                            !           CV_KNOD2 = ndgln%cv((ELE2-1)*Mdims%cv_nloc+CV_KLOC2)
                            !           LOC2_FEMT(:, CV_KLOC) = FEMT_ALL(:, CV_KNOD2)
                            !           LOC2_FEMTOLD(:, CV_KLOC) = FEMTOLD_ALL(:, CV_KNOD2)
                            !           IF (use_volume_frac_T2) THEN
                            !               LOC2_FEMT2(:, CV_KLOC) = FEMT2_ALL(:, CV_KNOD2)
                            !               LOC2_FEMT2OLD(:, CV_KLOC) = FEMT2OLD_ALL(:, CV_KNOD2)
                            !           END IF
                            !       END DO
                            !   END IF ! ENDOF IF (between_elements) THEN
                              IF ( on_domain_boundary ) THEN
                                  DO U_SKLOC = 1, Mdims%u_snloc
                                      U_KLOC = U_SLOC2LOC( U_SKLOC )
                                      U_NODK = ndgln%u(( ELE - 1 ) * Mdims%u_nloc + U_KLOC )
                                      U_SNODK = ( SELE - 1 ) * Mdims%u_snloc + U_SKLOC
                                      SLOC_NU(:, :, U_SKLOC) = NU_ALL(:, 1:final_phase, U_NODK)
                                  END DO
                              END IF
                              !------------------
                              If_GOT_DIFFUS2: IF ( GOT_DIFFUS ) THEN
                                  ! This sub caculates the effective diffusion
                                  ! coefficient DIFF_COEF_DIVDX
                                  AUX_T  =>T_ALL(1:final_phase, CV_NODJ)
                                  LOC_WIC_T_BC_ALL=0
                                  !           IF(SELE.NE.0) THEN
                                  IF(on_domain_boundary) THEN
                                      DO IPHASE=1,final_phase
                                          LOC_WIC_T_BC_ALL(IPHASE)=WIC_T_BC_ALL(1, IPHASE, SELE)
                                          IF(LOC_WIC_T_BC_ALL(IPHASE)==WIC_T_BC_DIRICHLET) THEN
                                              AUX_T( IPHASE ) = SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + Mdims%cv_snloc*( SELE- 1) )
                                          ENDIF
                                      END DO
                                  ENDIF
                                  CALL DIFFUS_CAL_COEFF( DIFF_COEF_DIVDX,  &
                                      Mdims%cv_nloc, Mdims%mat_nloc, final_phase, ndgln%mat, &
                                      CV_funs%scvfen, CV_funs%scvfen, GI, Mdims%ndim, TDIFFUSION, &
                                      HDC, &
                                      AUX_T, LOC_T_I, &
                                      ELE, ELE2, CVNORMX_ALL( :, GI ), &
                                      DTX_ELE_ALL(:,:,:,ELE),  DTX_ELE_ALL(:,:,:,MAX(1,ELE2)), &
                                      LOC_WIC_T_BC_ALL, CV_OTHER_LOC, MAT_OTHER_LOC, Mdims%cv_snloc, CV_SLOC2LOC, &
                                      on_domain_boundary, between_elements )
                                      
                              ELSE 
                                  DIFF_COEF_DIVDX = 0.0
                              END IF If_GOT_DIFFUS2
                              ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI.
                              IF( GOT_T2 ) THEN
                                CALL GET_INT_VEL_POROUS_VEL( NDOTQNEW, NDOTQ, INCOME, &
                                    LOC_T2_I, LOC_T2_J, LOC_NU, SLOC_NU, UGI_COEF_ELE_ALL, &
                                    upwnd%inv_adv_coef(1:final_phase,MAT_NODI), upwnd%inv_adv_coef(1:final_phase,MAT_NODJ), &
                                    NUGI_ALL)
                              ELSE
                                CALL GET_INT_VEL_POROUS_VEL( NDOTQNEW, NDOTQ, INCOME, &
                                    LOC_T_I, LOC_T_J, LOC_NU, SLOC_NU, UGI_COEF_ELE_ALL, &
                                    upwnd%inv_adv_coef(1:final_phase,MAT_NODI), upwnd%inv_adv_coef(1:final_phase,MAT_NODJ), &
                                    NUGI_ALL)
                              ENDIF
                              !Calculate the courant number for porous media
                              if (present(Courant_number) )then 
                                if (.not. on_domain_boundary) then
                                  !ndotq = velocity * normal                     !In the wells the flow is too fast and makes this misleading
                                  Courant_number(1) = max(Courant_number(1), abs ( dt * maxval(ndotq(1:final_phase)) / (VOLFRA_PORE( 1, ELE ) * hdc)))
                                  !and the shock-front Courant number
                                  if (shock_front_in_ele(ele, Mdims, T_ALL, ndgln, Imble_frac(:, CV_INOD))) then
                                      !ndotq = velocity * normal
                                      Courant_number(2) = max(Courant_number(2), abs ( dt * maxval(ndotq(1:final_phase)) / (VOLFRA_PORE( 1, ELE ) * hdc)))
                                  end if
                                end if
                              end if
                              If_GOT_CAPDIFFUS: IF ( VAD_activated ) THEN
                                  IF(SELE == 0) THEN
                                    CAP_DIFF_COEF_DIVDX = 0.
                                    !Project permeability at the GI point
                                    auxR = 0.
                                    if (has_anisotropic_permeability) then 
                                        forall (idim = 1:Mdims%ndim, iv_idim = 1:Mdims%ndim)
                                            auxR = auxR  + CVNORMX_ALL(idim, GI) * perm%val(idim, iv_idim,ele) * CVNORMX_ALL(iv_idim, GI)
                                        end forall                                        
                                    else
                                        do idim =1, Mdims%ndim 
                                            auxR = auxR + CVNORMX_ALL(idim, GI) * perm%val(idim,idim,ele) * CVNORMX_ALL(idim, GI)
                                        end do
                                    end if
                                    do iphase =1, final_phase
                                        rsum_nodi(iphase) = upwnd%inv_adv_coef(iphase,MAT_NODI)*auxR
                                    end do
                                    do iphase =1, final_phase
                                        rsum_nodj(iphase) = upwnd%inv_adv_coef(iphase,MAT_NODJ)*auxR
                                    end do
                                    CAP_DIFF_COEF_DIVDX = (CAP_DIFFUSION( :, MAT_NODI )&
                                        * rsum_nodi*(1.-INCOME) +&
                                        CAP_DIFFUSION( :, MAT_NODJ ) * rsum_nodj * INCOME) /HDC
                                  ELSE
                                      CAP_DIFF_COEF_DIVDX = 0.0
                                  ENDIF
                                  !Used normalised flux to disable/enable VAD for certain directions
                                  if (flux_limited_vad) CAP_DIFF_COEF_DIVDX(phase_with_pc) = &
                                    CAP_DIFF_COEF_DIVDX(phase_with_pc) * abs(dot_product(NUGI_ALL(:,phase_with_pc),&
                                      CVNORMX_ALL(:, GI))/sqrt(sum(NUGI_ALL(:,phase_with_pc)**2.) + 1e-16))
                                  !Distribute the capillary coefficient over the phases to ensure mass conservation
                                  !This is very important as it allows to use the over-relaxation parameter safely
                                  !and reduce the cost of using capillary pressure in several orders of magnitude
                                  CAP_DIFF_COEF_DIVDX(1:final_phase) =  CAP_DIFF_COEF_DIVDX(phase_with_pc)/Mdims%n_in_pres

                              ELSE
                                  CAP_DIFF_COEF_DIVDX = 0.0
                              END IF If_GOT_CAPDIFFUS
                              ! Pack ndotq information:

                              !Use upwinding to obtaing the values                        
                            IF ( on_domain_boundary ) THEN
                                !tracer
                                do iphase = 1, final_phase
                                if ( WIC_T_BC_ALL( 1,iphase, SELE ) /= WIC_T_BC_DIRICHLET )then
                                    LIMT(iphase) = LOC_T_I(iphase)
                                ELSE 
                                    LIMT(iphase) = LOC_T_I(iphase) * (1.0-INCOME(iphase)) + INCOME(iphase)* SUF_T_BC_ALL( 1, iphase, CV_SILOC + Mdims%cv_snloc*( SELE- 1) )
                                END if
                                !Density
                                if ( WIC_D_BC_ALL( 1, iphase, SELE ) /= WIC_D_BC_DIRICHLET ) then
                                    LIMD(iphase) = LOC_DEN_I(iphase)
                                ELSE 
                                    LIMD(iphase) = LOC_DEN_I(iphase) * (1.0-INCOME(iphase)) + INCOME(iphase)* SUF_D_BC_ALL( 1, iphase, CV_SILOC + Mdims%cv_snloc*( SELE- 1) )
                                END if  
                                !Saturation
                                if (use_volume_frac_T2) then 
                                    if ( WIC_T2_BC_ALL( 1, iphase, SELE ) /= WIC_T_BC_DIRICHLET ) then
                                    LIMT2(iphase) = LOC_T2_I(iphase)
                                    ELSE 
                                    LIMT2(iphase) = LOC_T2_I(iphase) * (1.0-INCOME(iphase)) + INCOME(iphase)* SUF_T2_BC_ALL( 1, iphase, CV_SILOC + Mdims%cv_snloc*( SELE- 1) )
                                    END if
                                end if  
                                end do                                                                      
                            else
                                LIMT=LOC_T_I*(1.0-INCOME) + LOC_T_J*INCOME
                                LIMD=LOC_DEN_I*(1.0-INCOME) + DEN_ALL(1:final_phase, cv_nodj)*INCOME
                                if (use_volume_frac_T2) LIMT2=LOC_T2_I*(1.0-INCOME) + LOC_T2_J*INCOME
                            end if

                              LIMDT=LIMD*LIMT
                              ! Make allowances for no matrix stencil operating from outside the boundary.
                              BCZERO=1.0
                              IF( on_domain_boundary ) BCZERO=1.0-INCOME
                              !====================== ACV AND RHS ASSEMBLY ===================
                              Conditional_GETCT2: IF ( GETCT ) THEN ! Obtain the CV discretised Mmat%CT eqations plus RHS
                                  ct_rhs_phase_cv_nodi=0.0; ct_rhs_phase_cv_nodj=0.0

                                  CALL PUT_IN_CT_RHS(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV, ct_rhs_phase_cv_nodi, ct_rhs_phase_cv_nodj, &
                                      final_phase, Mdims, CV_funs, ndgln, Mmat, GI,  &
                                      between_elements, on_domain_boundary, ELE, ELE2, SELE, HDC, MASS_ELE, &
                                      JCOUNT_KLOC, JCOUNT_KLOC2, ICOUNT_KLOC, ICOUNT_KLOC2, C_JCOUNT_KLOC, C_JCOUNT_KLOC2, &
                                      C_ICOUNT_KLOC, C_ICOUNT_KLOC2, U_SLOC2LOC, CV_SLOC2LOC,&
                                      SdevFuns%DETWEI, CVNORMX_ALL, DEN_ALL(1:final_phase,:), CV_NODI, CV_NODJ, &
                                      WIC_U_BC_ALL, WIC_P_BC_ALL, pressure_BCs%val, &
                                      UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL,  &
                                      NDOTQNEW, LIMT, LIMDT, LIMT_HAT, LIMT2, NDOTQ_HAT, &
                                      integrate_other_side_and_not_boundary, &
                                      loc_u, THETA_VEL, &
                                      rdum_ndim_nphase_1, rdum_nphase_1, rdum_nphase_2, rdum_nphase_3, X_ALL, SUF_D_BC_ALL, gravty)

                                  !Only for ipres = 1
                                  call addto(Mmat%CT_RHS,1,cv_nodi,sum(ct_rhs_phase_cv_nodi))
                                  if ( integrate_other_side_and_not_boundary ) then
                                      call addto(Mmat%CT_RHS,1,cv_nodj,sum(ct_rhs_phase_cv_nodj))
                                  end if

                              ENDIF Conditional_GETCT2

                              Conditional_GETCV_DISC: IF ( GETCV_DISC ) THEN
                                  ! Obtain the CV discretised advection/diffusion equations
                                  IF( on_domain_boundary ) then
                                      where ( WIC_T_BC_ALL(1,:,SELE) == WIC_T_BC_ROBIN )
                                          !Robin 1 contains the value of the field outside the domain times the coefficient.
                                          !This is also used for Neumann to impose a fixed flux
                                          ROBIN1 = SUF_T_BC_ROB1_ALL(1,1:final_phase, CV_SILOC+Mdims%cv_snloc*(sele-1))
                                          ROBIN2 = SUF_T_BC_ROB2_ALL(1,1:final_phase, CV_SILOC+Mdims%cv_snloc*(sele-1))
                                    else where
                                        ROBIN1=0.0; ROBIN2=0.0
                                      end where
                                  END IF
                                  LOC_CV_RHS_I=0.0; LOC_MAT_II =0.
                                  LOC_CV_RHS_J=0.0; LOC_MAT_JJ =0.
                                  LOC_MAT_IJ = 0.0; LOC_MAT_JI =0.
                                  IF ( GETMAT ) THEN
                                      ! - Calculate the integration of the limited, high-order flux over a face
                                      ! Conservative discretisation. The matrix (PIVOT ON LOW ORDER SOLN)
                                      IF ( on_domain_boundary) THEN
                                        if (GOT_DIFFUS) then
                                          DO IPHASE=1,final_phase
                                              IF(WIC_T_BC_ALL(1,iphase,sele) == WIC_T_BC_DIRICHLET) THEN
                                                  LOC_CV_RHS_I( IPHASE ) =  LOC_CV_RHS_I( IPHASE ) &
                                                      + LIMT2(IPHASE) * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE) &
                                                      * SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + Mdims%cv_snloc*( SELE- 1))
                                                  IF(GET_GTHETA) THEN
                                                      THETA_GDIFF( IPHASE, CV_NODI ) =  THETA_GDIFF( IPHASE, CV_NODI ) &
                                                          + LIMT2(IPHASE) * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE) &
                                                          * SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + Mdims%cv_snloc*( SELE- 1) )
                                                  END IF
                                              END IF
                                          END DO
                                        end if
                                      ELSE                                      
                                        !Assemble off-diagonal cv_nodi-cv_nodj
                                        LOC_MAT_IJ = LOC_MAT_IJ + LIMT2 * SdevFuns%DETWEI( GI ) * NDOTQNEW * INCOME * LIMD! Advection
                                        if (GOT_DIFFUS) LOC_MAT_IJ = LOC_MAT_IJ - LIMT2 * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX
                                        if (VAD_activated) LOC_MAT_IJ = LOC_MAT_IJ - LIMT2 * SdevFuns%DETWEI( GI ) * CAP_DIFF_COEF_DIVDX
                                        !Assemble off-diagonal cv_nodj-cv_nodi, integrate the other CV side contribution (the sign is changed)...
                                        if(integrate_other_side_and_not_boundary) then
                                          LOC_MAT_JI = LOC_MAT_JI - LIMT2 * SdevFuns%DETWEI( GI ) * NDOTQNEW * (1. - INCOME) * LIMD! Advection
                                          if (GOT_DIFFUS) LOC_MAT_JI = LOC_MAT_JI - LIMT2 * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX
                                          if (VAD_activated) LOC_MAT_JI = LOC_MAT_JI - LIMT2 * SdevFuns%DETWEI( GI ) * CAP_DIFF_COEF_DIVDX
                                        end if

                                          IF ( GET_GTHETA ) THEN
                                              THETA_GDIFF( :, CV_NODI ) =  THETA_GDIFF( :, CV_NODI ) &
                                                  + LIMT2 * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX * LOC_T_J ! Diffusion contribution
                                              ! integrate the other CV side contribution (the sign is changed)...
                                              if(integrate_other_side_and_not_boundary) then
                                                  THETA_GDIFF( :, CV_NODJ ) =  THETA_GDIFF( :, CV_NODJ ) &
                                                      + LIMT2 * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX * LOC_T_I ! Diffusion contribution
                                              endif
                                          END IF
                                      END IF ! endif of IF ( on_domain_boundary ) THEN ELSE

                                      !Assemble diagonal of the matrix of node cv_nodi
                                      LOC_MAT_II = LOC_MAT_II +  LIMT2 * SdevFuns%DETWEI( GI ) * NDOTQNEW * ( 1. - INCOME ) * LIMD! Advection
                                      if (GOT_DIFFUS) LOC_MAT_II = LOC_MAT_II + LIMT2 * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX
                                      if (VAD_activated) LOC_MAT_II = LOC_MAT_II + LIMT2 * SdevFuns%DETWEI( GI ) * CAP_DIFF_COEF_DIVDX
                                      if (.not.conservative_advection) LOC_MAT_II = LOC_MAT_II - LIMT2 * ( ONE_M_CV_BETA ) * &
                                                                                    SdevFuns%DETWEI( GI ) * NDOTQNEW * LIMD
                                      if (on_domain_boundary) LOC_MAT_II = LOC_MAT_II + SdevFuns%DETWEI( GI ) * ROBIN2

                                      !Assemble diagonal of the matrix of node cv_nodj
                                      if(integrate_other_side_and_not_boundary) then
                                        LOC_MAT_JJ = LOC_MAT_JJ -  LIMT2 * SdevFuns%DETWEI( GI ) * NDOTQNEW * INCOME * LIMD! Advection
                                        if (GOT_DIFFUS) LOC_MAT_JJ = LOC_MAT_JJ + LIMT2 * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX
                                        if (VAD_activated) LOC_MAT_JJ = LOC_MAT_JJ +  LIMT2 * SdevFuns%DETWEI( GI ) * CAP_DIFF_COEF_DIVDX
                                        if (.not.conservative_advection) LOC_MAT_JJ = LOC_MAT_JJ + LIMT2 * ( ONE_M_CV_BETA ) * SdevFuns%DETWEI( GI ) * NDOTQNEW * LIMD
                                      endif

                                      IF ( GET_GTHETA ) THEN
                                          THETA_GDIFF( :, CV_NODI ) =  THETA_GDIFF( :, CV_NODI ) &
                                              -  LIMT2 * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX * LOC_T_I & ! Diffusion contribution
                                              -  SdevFuns%DETWEI( GI ) * robin1 * LOC_T_I  ! Robin bc
                                          if(integrate_other_side_and_not_boundary) then
                                              THETA_GDIFF( :, CV_NODJ ) =  THETA_GDIFF( :, CV_NODJ ) &
                                                  -  LIMT2 * SdevFuns%DETWEI( GI ) * DIFF_COEF_DIVDX * LOC_T_J ! Diffusion contribution
                                          endif
                                      END IF
                                  END IF  ! ENDOF IF ( GETMAT ) THEN

                                  ! Put results into the RHS vector
                                  LOC_CV_RHS_I =  LOC_CV_RHS_I  &
                                         ! subtract 1st order adv. soln.
                                      + LIMT2 * NDOTQNEW * SdevFuns%DETWEI( GI ) * LIMDT * BCZERO &
                                      -  SdevFuns%DETWEI( GI ) * ( LIMT2 * NDOTQNEW * LIMDT)
                                  ! Subtract out 1st order term non-conservative adv.
                                      if (VAD_activated) LOC_CV_RHS_I =  LOC_CV_RHS_I &
                                          - LIMT2* SdevFuns%DETWEI(GI) * CAP_DIFF_COEF_DIVDX &  ! capillary pressure stabilization term..
                                          * ( LOC_T_J - LOC_T_I )
                                      if (.not.conservative_advection) LOC_CV_RHS_I =  LOC_CV_RHS_I &
                                          - LIMT2 * ( ONE_M_CV_BETA ) * SdevFuns%DETWEI( GI ) * NDOTQNEW * LIMD * LOC_T_I &
                                          + ( ONE_M_CV_BETA) * SdevFuns%DETWEI( GI ) &
                                          * ( LIMT2 * NDOTQNEW * LOC_T_I * LIMD )
                                      if (on_domain_boundary) LOC_CV_RHS_I =  LOC_CV_RHS_I &
                                          + SdevFuns%DETWEI( GI ) * ROBIN1


                                  if(integrate_other_side_and_not_boundary) then
                                      LOC_CV_RHS_J =  LOC_CV_RHS_J  &
                                             ! subtract 1st order adv. soln.
                                          - LIMT2 * NDOTQNEW * SdevFuns%DETWEI( GI ) * LIMDT * BCZERO &
                                          +  SdevFuns%DETWEI( GI ) * ( LIMT2 * NDOTQNEW * LIMDT)
                                      if (VAD_activated) LOC_CV_RHS_J =  LOC_CV_RHS_J  &
                                          - LIMT2 * SdevFuns%DETWEI(GI) * CAP_DIFF_COEF_DIVDX & ! capillary pressure stabilization term..
                                          * ( LOC_T_I - LOC_T_J )
                                      if (.not.conservative_advection) LOC_CV_RHS_J =  LOC_CV_RHS_J  &
                                          + LIMT2 * ( ONE_M_CV_BETA ) * SdevFuns%DETWEI( GI ) * NDOTQNEW * LIMD * LOC_T_J &
                                          - ( ONE_M_CV_BETA) * SdevFuns%DETWEI( GI ) &
                                          * ( LIMT2 * NDOTQNEW * LOC_T_J * LIMD )

                                  endif
                                !   IF ( GET_GTHETA ) THEN
                                !       THETA_GDIFF( :, CV_NODI ) =  THETA_GDIFF( :, CV_NODI ) &
                                !           + (1.-LIMT2) * SdevFuns%DETWEI(GI) * DIFF_COEFOLD_DIVDX &
                                !           * ( LOC_TOLD_J - LOC_TOLD_I ) &
                                !           ! Robin bc
                                !           + SdevFuns%DETWEI( GI ) * robin2
                                !       if(integrate_other_side_and_not_boundary) then
                                !           THETA_GDIFF( :, CV_NODJ ) =  THETA_GDIFF( :, CV_NODJ ) &
                                !               + (1.-LIMT2) * SdevFuns%DETWEI(GI) * DIFF_COEFOLD_DIVDX &
                                !               * ( LOC_TOLD_I - LOC_TOLD_J )
                                !       endif
                                !   END IF
                                  ! this is for the internal energy equation source term..
                                  ! This is to introduce the compressibility term due to expansion and therefore the divergence of the velocity is non-zero
                                  !for wells this is not straightforward <= need to CHANGE THIS FOR COMPRESSIBILITY
                                  IF ( THERMAL .and. Mdims%npres == 1) THEN
                                      LOC_CV_RHS_I = LOC_CV_RHS_I&
                                          - CV_P( 1, 1, CV_NODI ) * SdevFuns%DETWEI( GI ) * ( &
                                           NDOTQNEW * LIMT2 )
                                      if ( integrate_other_side_and_not_boundary ) then
                                          LOC_CV_RHS_J = LOC_CV_RHS_J &
                                              + CV_P( 1, 1, CV_NODJ ) * SdevFuns%DETWEI( GI ) * ( &
                                              NDOTQNEW * LIMT2 )
                                      end if
                                  END IF ! THERMAL

                                  do iphase = 1, final_phase
                                      assembly_phase = iphase
                                      !For the RHS collapsing to assemble into phase 1 can be done just here
                                      if (loc_assemble_collapsed_to_one_phase) assembly_phase = 1
                                      call addto(Mmat%CV_RHS,assembly_phase, CV_NODI,LOC_CV_RHS_I(iphase))
                                      call addto(Mmat%CV_RHS,assembly_phase, CV_NODJ,LOC_CV_RHS_J(iphase))
                                      !Introduce the information into the petsc_ACV matrix
                                      call addto(Mmat%petsc_ACV,assembly_phase,assembly_phase,cv_nodi,cv_nodi, LOC_MAT_II(iphase) )
                                      call addto(Mmat%petsc_ACV,assembly_phase,assembly_phase,cv_nodj,cv_nodj, LOC_MAT_JJ(iphase) )
                                      call addto(Mmat%petsc_ACV,assembly_phase,assembly_phase,cv_nodi,cv_nodj, LOC_MAT_IJ(iphase) )
                                      call addto(Mmat%petsc_ACV,assembly_phase,assembly_phase,cv_nodj,cv_nodi, LOC_MAT_JI(iphase) )
                                  end do
                              ENDIF Conditional_GETCV_DISC

                            !Finally store fluxes across all the boundaries either for mass conservation check or mass outflux
                            if (compute_outfluxes .and. on_domain_boundary) then
                              call update_outfluxes(bcs_outfluxes,outfluxes, sele, cv_nodi,  SdevFuns%DETWEI(gi), &
                                ndotqnew * SdevFuns%DETWEI(gi) * LIMT, ndotqnew * SdevFuns%DETWEI(gi) * LIMDT, & !Vol_flux and Mass_flux
                                old_tracer, outfluxes_fields, 1, final_phase )
                            end if

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
          !Add compressibility to the transport equation/add time derivative term
          Conditional_GETCV_DISC2: IF( GETCV_DISC ) THEN ! Obtain the CV discretised advection/diffusion equations
              Loop_CVNODI2: DO CV_NODI = 1, Mdims%cv_nonods ! Put onto the diagonal of the matrix

                ! Generate local variables (to avoid slicing) ***************
                LOC_T_I = T_ALL(1:final_phase, cv_nodi); LOC_TOLD_I = TOLD_ALL(1:final_phase, cv_nodi)
                LOC_DEN_I =DEN_ALL(1:final_phase, cv_nodi); LOC_DENOLD_I = DENOLD_ALL(1:final_phase, cv_nodi)
                if (use_volume_frac_T2) then
                  LOC_T2_I = T2_ALL(1:final_phase, cv_nodi); LOC_T2OLD_I = T2OLD_ALL(1:final_phase, cv_nodi)
                end if


                  LOC_CV_RHS_I=0.0; LOC_MAT_II =0.
                  R_PHASE = MEAN_PORE_CV( 1, CV_NODI ) * Mass_CV( CV_NODI ) / DT
                  IF ( THERMAL .and. Mdims%npres == 1) THEN
                      LOC_CV_RHS_I = LOC_CV_RHS_I &
                          - CV_P( 1, 1, CV_NODI ) * ( Mass_CV( CV_NODI ) / DT ) * ( LOC_T2_I - LOC_T2OLD_I)
                  END IF

                  IF ( GOT_T2 ) THEN
                    DO IPHASE = 1,final_phase!to avoid slicing sourct_all
                      LOC_CV_RHS_I(iphase) = LOC_CV_RHS_I(iphase)  + Mass_CV(CV_NODI) * SOURCT_ALL( iphase, CV_NODI )
                    end do
                      if (thermal) then
                          !In this case for the time-integration term the effective rho Cp is a combination of the porous media
                          ! and the fluids. Here we add the porous media contribution. Multiplied by the saturation so we use the same
                          !paradigm that for the phases, but in the equations it isn't, but here because we iterate over phases and collapse
                          !this is required
                          LOC_MAT_II = LOC_MAT_II + porous_heat_coef( CV_NODI ) * LOC_T2_I &
                                  * R_PHASE * (1-MEAN_PORE_CV( 1, CV_NODI ))/MEAN_PORE_CV( 1, CV_NODI )

                          !R_PHASE includes the porosity. Since in this case we are interested in what is NOT porous
                              !we divide to remove that term and multiply by the correct term (1-porosity)
                          LOC_CV_RHS_I=LOC_CV_RHS_I  &
                          + (CV_BETA * porous_heat_coef_old( CV_NODI ) * LOC_T2OLD_I &
                          + (ONE_M_CV_BETA) * porous_heat_coef( CV_NODI ) * LOC_T2_I ) &
                          * R_PHASE * LOC_TOLD_I* (1-MEAN_PORE_CV( 1, CV_NODI ))/MEAN_PORE_CV( 1, CV_NODI )

                      end if

                      LOC_MAT_II = LOC_MAT_II + LOC_DEN_I * LOC_T2_I * R_PHASE
                      LOC_CV_RHS_I=LOC_CV_RHS_I  + (CV_BETA * LOC_DENOLD_I * LOC_T2OLD_I + &
                              (ONE_M_CV_BETA) * LOC_DEN_I * LOC_T2_I ) * R_PHASE * LOC_TOLD_I

                  ELSE

                    !Diagonal term Vol/dt * rho and accompaniying rhs term
                    LOC_MAT_II = LOC_MAT_II + LOC_DEN_I * R_PHASE
                    do iphase = 1, final_phase
                      LOC_CV_RHS_I(IPHASE)=LOC_CV_RHS_I(IPHASE)  &
                          + Mass_CV( CV_NODI ) * SOURCT_ALL( IPHASE, CV_NODI )&
                          + ( CV_BETA * LOC_DENOLD_I(IPHASE) &
                          + (ONE_M_CV_BETA) * LOC_DEN_I(IPHASE) ) &
                          * R_PHASE(IPHASE) * LOC_TOLD_I(IPHASE)
                   END DO
                  END IF


                  Conditional_GETMAT2: IF ( GETMAT .and. have_absorption) THEN

                    !Absorption directly introduced into Mmat%petsc_ACV
                    DO jphase=1,final_phase
                        do iphase=1, final_phase
                             call addto(Mmat%petsc_ACV,IPHASE,JPHASE, &
                             cv_nodi, cv_nodi, Mass_CV( CV_NODI ) * ABSORBT_ALL( iphase, jphase, CV_NODI ))
                      end do
                    end do
                  END IF Conditional_GETMAT2

                  do iphase = 1, final_phase
                      assembly_phase = iphase
                      !For the RHS collapsing to assemble into phase 1 can be done just here
                      if (loc_assemble_collapsed_to_one_phase) assembly_phase = 1
                      call addto(Mmat%CV_RHS,assembly_phase, CV_NODI,LOC_CV_RHS_I(IPHASE))
                      !Introduce the information into the petsc_ACV matrix
                      call addto(Mmat%petsc_ACV,assembly_phase,assembly_phase,cv_nodi,cv_nodi, LOC_MAT_II(iphase) )
                end do
              END DO Loop_CVNODI2
          END IF Conditional_GETCV_DISC2

          IF ( GETCT ) THEN
            W_SUM_ONE1 = 0.0 !If == 1.0 applies constraint to T
            if (Solve_all_phases) W_SUM_ONE1 = 1.0
            W_SUM_ONE2 = 0.0 !If == 1.0 applies constraint to TOLD !sprint_to_do Unnecessary, should be removed
            DIAG_SCALE_PRES = 0.0
            DO CV_NODI = 1, Mdims%cv_nonods
              ! Generate local variables (to avoid slicing) ***************
              LOC_T_I = T_ALL(1:final_phase, cv_nodi); LOC_TOLD_I = TOLD_ALL(1:final_phase, cv_nodi)
              LOC_DEN_I =DEN_ALL(1:final_phase, cv_nodi); LOC_DENOLD_I = DENOLD_ALL(1:final_phase, cv_nodi)
              if (use_volume_frac_T2) then
                LOC_T2_I = T2_ALL(1:final_phase, cv_nodi)
              end if

              ct_rhs_phase=0.0 ; DIAG_SCALE_PRES_phase=0.0
              IPRES=1
              R_PRES(IPRES) = MASS_CV(CV_NODI ) * MEAN_PORE_CV( IPRES, CV_NODI ) / DT
              ! Add constraint to force sum of volume fracts to be unity...
                 ! W_SUM_ONE==1 applies the constraint
                 ! W_SUM_ONE==0 does NOT apply the constraint
              if ( Solve_all_phases) call addto(Mmat%CT_RHS,IPRES,cv_nodi,&
                  - ( W_SUM_ONE1 - W_SUM_ONE2 ) * R_PRES(IPRES))
              do iphase = 1, final_phase
                ct_rhs_phase(iphase)=ct_rhs_phase(iphase) &
                    - R_PRES(1) * ( &
                    + (1.0-W_SUM_ONE1) * LOC_T_I(iphase) - (1.0-W_SUM_ONE2) * LOC_TOLD_I(iphase) &
                    + (( LOC_TOLD_I(iphase) * ( LOC_DEN_I(iphase) - LOC_DENOLD_I(iphase) ) &
                    - DERIV( iphase, CV_NODI ) * CV_P( 1, 1, CV_NODI ) ) * T_ALL( iphase, CV_NODI ) ) / LOC_DEN_I(iphase) )
                DIAG_SCALE_PRES_phase( iphase ) = DIAG_SCALE_PRES_phase( iphase ) &
                    + MEAN_PORE_CV( 1, CV_NODI ) * LOC_T_I( iphase ) * DERIV( iphase, CV_NODI ) / ( DT * LOC_DEN_I(iphase) )
                ct_rhs_phase(iphase)=ct_rhs_phase(iphase)  &
                    + Mass_CV(CV_NODI ) * SOURCT_ALL( iphase, CV_NODI ) / LOC_DEN_I(iphase)
                IF ( HAVE_ABSORPTION ) THEN
                     DO JPHASE = 1, final_phase
                        ct_rhs_phase(iphase)=ct_rhs_phase(iphase)  &
                           - Mass_CV( CV_NODI ) * ABSORBT_ALL( iphase, JPHASE, CV_NODI ) * LOC_T_I( JPHASE ) / LOC_DEN_I(iphase)
                   END DO
                END IF
              end do

              !Introduce into the RHS
              call addto(Mmat%CT_RHS, 1, cv_nodi, SUM( ct_rhs_phase) )
              !and diagonal scaling
              DIAG_SCALE_PRES( 1,CV_NODI ) = DIAG_SCALE_PRES( 1,CV_NODI ) + sum( DIAG_SCALE_PRES_phase)
            END DO
          END IF

             !deallocate(R_PRES,R_PHASE,MEAN_PORE_CV_PHASE)
          !Assemble the part of the wells matrix and create corresponding RHS, absoprtions, etc.
          if (Mdims%npres >1) call ASSEMBLE_PIPE_TRANSPORT_AND_CTY( state, packed_state, tracer, den_all, denold_all, &
                                final_phase, &! final_phase => reservoir domain
                                Mdims, ndgln, DERIV, CV_P, SOURCT_ALL, ABSORBT_ALL, WIC_T_BC_ALL,WIC_D_BC_ALL, WIC_U_BC_ALL, &
                                SUF_T_BC_ALL,SUF_D_BC_ALL,SUF_U_BC_ALL, getcv_disc, getct, Mmat, Mspars, upwnd, GOT_T2, DT, &
                                pipes_aux, DIAG_SCALE_PRES_COUP, DIAG_SCALE_PRES,mean_pore_cv, eles_with_pipe, thermal,&
                                CV_BETA, MASS_CV, INV_B, MASS_ELE, bcs_outfluxes, outfluxes,&
                                porous_heat_coef, loc_assemble_collapsed_to_one_phase )

          if ( compute_outfluxes) then
              !Calculate final outfluxes and mass balance in the domain
              call mass_conservation_check_and_outfluxes(calculate_mass_delta, outfluxes, DEN_ALL, 2)
          endif
          ! Deallocating temporary working arrays
          call deallocate_multi_dev_shape_funs(SdevFuns)
          
        call deallocate(tracer_BCs)
        call deallocate(tracer_BCs_robin2)
        call deallocate(density_BCs)
        call deallocate(velocity_BCs)
        call deallocate(pressure_BCs)
        if (present(saturation)) then
            call deallocate(saturation_BCs)
            call deallocate(saturation_BCs_robin2)
        end if
        
          !These three variables are allocated simultaneously so only one need to be checked
          if (allocated(suf_t_bc)) deallocate( suf_t_bc, suf_t_bc_rob1, suf_t_bc_rob2)
          if (VAD_activated) deallocate(CAP_DIFFUSION)
          ewrite(3,*) 'Leaving CV_ASSEMB'
          if (allocated(bcs_outfluxes)) deallocate(bcs_outfluxes)

          RETURN
      contains



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

        !---------------------------------------------------------------------------
        !> @author Chris Pain, Pablo Salinas
        !> @brief Computes the flux between CVs for porous media. NDOTQNEW contains the fluxes for a given gauss integration point
        !---------------------------------------------------------------------------
        SUBROUTINE GET_INT_VEL_POROUS_VEL(NDOTQNEW, NDOTQ, INCOME, &
            LOC_T_I, LOC_T_J, LOC_NU, SLOC_NU, UGI_COEF_ELE_ALL,&
             I_inv_adv_coef, J_inv_adv_coef, UDGI_ALL)
            ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI.
            ! it assumes a compact_overlapping decomposition approach for velocity.
            IMPLICIT NONE
            REAL, DIMENSION( : ), intent( inout ) :: NDOTQNEW, NDOTQ, INCOME
            REAL, DIMENSION( : ), intent( in ) :: LOC_T_I, LOC_T_J
            REAL, DIMENSION( :, :, : ), intent( in ) ::  LOC_NU, SLOC_NU
            REAL, DIMENSION( :, :, : ), intent( inout ) :: UGI_COEF_ELE_ALL
            REAL, DIMENSION( : ), intent( in ) :: I_inv_adv_coef, J_inv_adv_coef
            REAL, DIMENSION( :, :  ), intent( inout ) :: UDGI_ALL

            UGI_COEF_ELE_ALL=0.0 
            Conditional_SELE: IF( on_domain_boundary ) THEN ! On the boundary of the domain.
                !Initialize variables
                forall (iv_iphase = 1:final_phase, iv_idim = 1:Mdims%ndim)
                    ROW_SUM_INV_VI(iv_idim,iv_iphase)=I_inv_adv_coef (iv_iphase) * SUM(perm%val(iv_idim,:,ele))
                end forall
                IF( WIC_P_BC_ALL( 1, 1, SELE) == WIC_P_BC_DIRICHLET ) THEN ! Pressure boundary condition
                  DO iv_iphase = 1,final_phase
                        !(vel * shape_functions)/sigma
                        ! if (is_P0DGP1) then 
                        UDGI_ALL(:, iv_iphase) = I_inv_adv_coef (iv_iphase)* matmul(perm%val(:,:,ele),&
                            LOC_NU( :, iv_iphase, 1 ) )
                        ! else
                        !     UDGI_ALL(:, iv_iphase) = I_inv_adv_coef (iv_iphase)* matmul(perm%val(:,:,ele),&
                        !         matmul(LOC_NU( :, iv_iphase, : ), CV_funs%sufen( :, GI )))
                        ! end if
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
                        DO iv_u_kloc = 1, Mdims%u_nloc
                            ! IF (.false.) THEN !<= this one for strong boundary conditions
                            IF (iv_incomming_flow) THEN ! Incomming...
                                UGI_COEF_ELE_ALL(:, iv_iphase, iv_u_kloc)=iv_SUF_SIG_DIAGTEN_BC_GI
                            ELSE
                                UGI_COEF_ELE_ALL(:, iv_iphase, iv_u_kloc)=1.0
                            ENDIF
                            UGI_COEF_ELE_ALL(:, iv_iphase, iv_u_kloc)= I_inv_adv_coef (iv_iphase) * matmul(perm%val(:,:,ele),UGI_COEF_ELE_ALL(:, iv_iphase, iv_u_kloc))
                        END DO
                        if(iv_incomming_flow) UDGI_ALL(:, iv_iphase) = UDGI_ALL(:, iv_iphase) * iv_SUF_SIG_DIAGTEN_BC_GI
                  end do
                ELSE ! Specified vel bc.
                  !Majority of the times the BCs have zero velocity, check that before doing anything to save time 
                  zero_vel_BC = maxval( abs(SUF_U_BC_ALL(:, :, &
                      1 + Mdims%u_snloc*(SELE-1) : Mdims%u_snloc*(SELE-1) + Mdims%u_snloc) ))<1e-30
                  UDGI_ALL = 0.0; UDGI_ALL_FOR_INV = 0.0; UGI_COEF_ELE_ALL = 0.0
                  DO iv_iphase = 1,final_phase
                      if (.not.zero_vel_BC) then 
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
                        UDGI_ALL(:, iv_iphase) = UDGI_ALL(:, iv_iphase)  + I_inv_adv_coef (iv_iphase) * matmul(perm%val(:,:,ele),UDGI_ALL_FOR_INV(:, iv_iphase))
                      end if
                  END DO ! PHASE LOOP
                END IF
            ELSE! IF( .not. between_elements) THEN!Only for same element/Continuous formulation !old flag: DISTCONTINUOUS_METHOD
                !vel(GI) = (vel * shape_functions)/sigma
                ! if (is_P0DGP1) then                  
                forall (iv_iphase = 1:final_phase, iv_idim = 1:Mdims%ndim)
                    UDGI_ALL(iv_idim, iv_iphase) = I_inv_adv_coef(iv_iphase)*LOC_NU( iv_idim, iv_iphase, 1 )
                    UDGI2_ALL(iv_idim, iv_iphase) = J_inv_adv_coef(iv_iphase)*LOC_NU( iv_idim, iv_iphase, 1 )
                end forall
                ! else
                !     do iv_iphase = 1,final_phase
                !         UDGI_ALL(:, iv_iphase) = I_inv_adv_coef(iv_iphase)*&
                !             matmul(LOC_NU( :, iv_iphase, : ), CV_funs%sufen( :, GI ))
                !         UDGI2_ALL(:, iv_iphase) = J_inv_adv_coef(iv_iphase)*&
                !             matmul(LOC_NU( :, iv_iphase, : ), CV_funs%sufen( :, GI ))
                !     end do
                ! end if
                !Get the projected velocity
                NDOTQ = 0.
                forall (iv_iphase = 1:final_phase, iv_idim = 1:Mdims%ndim)
                    NDOTQ(iv_iphase) = NDOTQ(iv_iphase) + CVNORMX_ALL(iv_idim, GI)* UDGI_ALL(iv_idim, iv_iphase)
                end forall                
                !Get the direction of the flow
                WHERE (NDOTQ < 0.0)
                    INCOME=1.0
                ELSE WHERE
                    INCOME=0.0
                END WHERE
                !Finally we calculate the velocity at the interface (still excluding permeability)
                DO iv_iphase = 1,final_phase
                    !Calculate contributions from each side
                    UDGI_ALL(:, iv_iphase) = UDGI_ALL(:, iv_iphase) * (1.0-INCOME(iv_iphase)) + UDGI2_ALL(:, iv_iphase) * INCOME(iv_iphase)
                END DO ! PHASE LOOP
                !Now apply the permeability to obtain finally UDGI_ALL
                if (has_anisotropic_permeability) then
                    DO iv_iphase = 1,final_phase
                        UDGI_ALL(:, iv_iphase) = matmul(UDGI_ALL(:, iv_iphase), perm%val(:,:,ele))
                    end do
                    do iv_idim = 1, Mdims%ndim
                        auxR = SUM(perm%val(iv_idim,:,ele))
                        RSUM = auxR
                        do iv_iphase = 1, final_phase
                            ROW_SUM_INV_VI(iv_idim,iv_iphase)=I_inv_adv_coef(iv_iphase) * auxR
                            ROW_SUM_INV_VJ(iv_idim,iv_iphase)=J_inv_adv_coef(iv_iphase) * RSUM
                        end do
                    end do
                else !Avoid tensor multiplication if possible 
                    forall (iv_iphase = 1:final_phase, iv_idim = 1:Mdims%ndim)
                        UDGI_ALL(iv_idim, iv_iphase) = UDGI_ALL(iv_idim, iv_iphase) * perm%val(iv_idim,iv_idim,ele)
                    end forall
                    forall (iv_iphase = 1:final_phase, iv_idim = 1:Mdims%ndim)
                        ROW_SUM_INV_VI(iv_idim,iv_iphase)=I_inv_adv_coef(iv_iphase) * perm%val(iv_idim,iv_idim,ele)
                        ROW_SUM_INV_VJ(iv_idim,iv_iphase)=J_inv_adv_coef(iv_iphase) * perm%val(iv_idim,iv_idim,ele)
                    end forall
                end if
                ! if (is_P0DGP1) then 
                forall (iv_iphase = 1:final_phase, iv_idim = 1:Mdims%ndim)
                    UGI_COEF_ELE_ALL(iv_idim, iv_iphase, 1)=ROW_SUM_INV_VI(iv_idim,iv_iphase)* (1.0-INCOME(iv_iphase)) &
                        + ROW_SUM_INV_VJ(iv_idim,iv_iphase)* INCOME(iv_iphase)
                end forall
                    ! else
                    !     DO iv_iphase = 1,final_phase
                    !         UGI_COEF_ELE_ALL(:, iv_iphase, :)=SPREAD(ROW_SUM_INV_VI(:,iv_iphase)* (1.0-INCOME(iv_iphase)) &
                    !             +ROW_SUM_INV_VJ(:,iv_iphase)* INCOME(iv_iphase), DIM=2, NCOPIES=Mdims%u_nloc)
                    !     END DO
                    ! end if
            ! ELSE !For discontinuous formulation (between elements but also can be used within). Does not use a TVD limiter but a weighting method
            !     !vel(GI) = (vel * shape_functions)/sigma
            !     do iv_iphase = 1,final_phase
            !         !Velocity including sigma
            !         UDGI_ALL(:, iv_iphase) = matmul(LOC_NU( :, iv_iphase, : ), CV_funs%sufen( :, GI ))
            !         !Normal flow including sigma, to know direction of flow
            !         NDOTQ(iv_iphase)  = dot_product( CVNORMX_ALL(:, GI),UDGI_ALL(:, iv_iphase))
            !         !Actual advection velocity by removing the contribution of sigma (still not normalised by permeability)
            !         UDGI_ALL(:, iv_iphase) = I_inv_adv_coef(iv_iphase) * UDGI_ALL(:, iv_iphase)
            !     end do
            !     do iv_iphase = 1,final_phase
            !         !Velocity including sigma
            !         UDGI2_ALL(:, iv_iphase) = matmul(LOC2_NU( :, iv_iphase, : ), CV_funs%sufen( :, GI ))
            !         !Normal flow including sigma, to know direction of flow
            !         NDOTQ2(iv_iphase) = dot_product( CVNORMX_ALL(:, GI),UDGI2_ALL(:, iv_iphase))
            !         !Actual advection velocity by removing the contribution of sigma (still not normalised by permeability)
            !         UDGI2_ALL(:, iv_iphase) = J_inv_adv_coef(iv_iphase)*UDGI2_ALL(:, iv_iphase)
            !     end do

            !     !Sigma averaged with the mass to be used as divisor
            !     iv_sigma_aver = 1./(I_adv_coef*MASS_CV_I+J_adv_coef*MASS_CV_J)
            !     do iv_iphase = 1,final_phase
            !       !Calculate the contribution of each side, considering sigma and the volume of the CVs
            !       if ( ( NDOTQ(iv_iphase) + NDOTQ2(iv_iphase) ) > 0.0 ) then
            !         !We redefine sigma so that it detects oscillations using first order taylor series
            !         iv_aux_tensor(iv_iphase) =  I_adv_coef(iv_iphase) + 0.5 * (LOC_T_J(iv_iphase) - LOC_T_I(iv_iphase)) * I_adv_coef_grad(iv_iphase)
            !         !We limit the value
            !         iv_aux_tensor(iv_iphase) = min(1000.*max(I_adv_coef(iv_iphase),  J_adv_coef(iv_iphase)), &
            !         max(0.001*min(I_adv_coef(iv_iphase),  J_adv_coef(iv_iphase)), iv_aux_tensor(iv_iphase) ))
            !         iv_aux_tensor(iv_iphase)= min( 1.0, I_inv_adv_coef(iv_iphase)* iv_aux_tensor(iv_iphase))
            !         !Calculate importance of each side
            !         iv_aux_tensor2(iv_iphase) = iv_aux_tensor(iv_iphase) * iv_sigma_aver(iv_iphase)* I_adv_coef(iv_iphase)*MASS_CV_I
            !         !iv_aux_tensor2 has to be calculated before since iv_aux_tensor is rewritten!
            !         iv_aux_tensor(iv_iphase) = (1.-iv_aux_tensor(iv_iphase)) + iv_aux_tensor(iv_iphase) * iv_sigma_aver(iv_iphase) * J_adv_coef(iv_iphase)*MASS_CV_J
            !       else
            !         !We redefine sigma so that it detects oscillations using first order taylor series
            !         iv_aux_tensor(iv_iphase) =  J_adv_coef(iv_iphase) + 0.5 * (LOC_T_I(iv_iphase) - LOC_T_J(iv_iphase)) * J_adv_coef_grad(iv_iphase)
            !         !We limit the value
            !         iv_aux_tensor(iv_iphase) = min(1000.*max(I_adv_coef(iv_iphase),  J_adv_coef(iv_iphase)), &
            !         max(0.001*min(I_adv_coef(iv_iphase),  J_adv_coef(iv_iphase)), iv_aux_tensor(iv_iphase) ))
            !         iv_aux_tensor(iv_iphase)= min( 1.0, J_inv_adv_coef(iv_iphase) * iv_aux_tensor(iv_iphase))
            !         !Calculate importance of each side
            !         iv_aux_tensor2(iv_iphase) = (1.-iv_aux_tensor(iv_iphase)) + iv_aux_tensor(iv_iphase) * iv_sigma_aver(iv_iphase) * I_adv_coef(iv_iphase)*MASS_CV_I
            !         !iv_aux_tensor2 has to be calculated before since iv_aux_tensor is rewritten!
            !         iv_aux_tensor(iv_iphase) = iv_aux_tensor(iv_iphase) * iv_sigma_aver(iv_iphase) * J_adv_coef(iv_iphase)*MASS_CV_J
            !       end if
            !       !Calculation of the velocity at the GI point
            !       UDGI_ALL(:, iv_iphase) = iv_aux_tensor(iv_iphase) * UDGI_ALL(:, iv_iphase) + iv_aux_tensor2(iv_iphase) * UDGI2_ALL(:, iv_iphase)
            !     end do
            !     !Now apply the permeability to obtain finally UDGI_ALL
            !     if (permeability_jump) then
            !       inv_harmonic_perm = 2.0*inverse((upwnd%inv_permeability(:,:,ele) + upwnd%inv_permeability(:,:,ele2)))
            !       DO iv_iphase = 1,final_phase
            !       !We use the harmonic average of the permeability. Because we have the inverse (and we want the inverse) this is simplified
            !         UDGI_ALL(:, iv_iphase) = matmul(UDGI_ALL(:, iv_iphase), inv_harmonic_perm)
            !       end do
            !     else
            !       DO iv_iphase = 1,final_phase
            !         UDGI_ALL(:, iv_iphase) = matmul(UDGI_ALL(:, iv_iphase), perm%val(:,:,ele))
            !       end do
            !     end if
            !     !Calculation of the coefficients at the GI point
            !     if (not_OLD_VEL) then
            !       do iv_idim = 1, Mdims%ndim
            !         auxR = SUM(perm%val(iv_idim,:,ele))
            !         RSUM = auxR
            !         if (between_elements) RSUM = SUM(perm%val(iv_idim,:,ele2))
            !         do iv_iphase = 1, final_phase
            !           ROW_SUM_INV_VI(iv_idim,iv_iphase)=I_inv_adv_coef(iv_iphase) * auxR
            !           ROW_SUM_INV_VJ(iv_idim,iv_iphase)=J_inv_adv_coef(iv_iphase) * RSUM
            !         end do
            !       end do
            !       IF( between_elements ) then
            !         DO iv_iphase = 1,final_phase
            !           UGI_COEF_ELE_ALL(:, iv_iphase, :) = iv_aux_tensor(iv_iphase)* SPREAD(ROW_SUM_INV_VI(:,iv_iphase), DIM=2, NCOPIES=Mdims%u_nloc)
            !           UGI_COEF_ELE2_ALL(:, iv_iphase, :) = iv_aux_tensor2(iv_iphase) *SPREAD(ROW_SUM_INV_VJ(:,iv_iphase), DIM=2, NCOPIES=Mdims%u_nloc)
            !         END DO
            !       else !same element
            !         DO iv_iphase = 1,final_phase
            !           UGI_COEF_ELE_ALL(:, iv_iphase, :) = iv_aux_tensor(iv_iphase) * SPREAD(ROW_SUM_INV_VI(:,iv_iphase), DIM=2, NCOPIES=Mdims%u_nloc) +&
            !           iv_aux_tensor2(iv_iphase) * SPREAD(ROW_SUM_INV_VJ(:,iv_iphase), DIM=2, NCOPIES=Mdims%u_nloc)
            !         END DO
            !       end if
            !   end if
            END IF Conditional_SELE

            ! ! Define whether flux is incoming or outgoing, depending on direction of flow
            NDOTQ = 0.
            forall (iv_iphase = 1:final_phase, iv_idim = 1:Mdims%ndim)
                NDOTQ(iv_iphase) = NDOTQ(iv_iphase) + CVNORMX_ALL(iv_idim, GI)* UDGI_ALL(iv_idim, iv_iphase)
            end forall 
            WHERE ( NDOTQ >= 0. )
                INCOME = 0.
            ELSE WHERE
                INCOME = 1.
            END WHERE
            ! Calculate NDOTQNEW from NDOTQ
            if (is_P0DGP1) then 
                NDOTQNEW = NDOTQ
                forall (iv_iphase = 1:final_phase, iv_idim = 1:Mdims%ndim)
                    NDOTQNEW(iv_iphase) = NDOTQNEW(iv_iphase) + CVNORMX_ALL(iv_idim, GI)* UGI_COEF_ELE_ALL(iv_idim, iv_iphase,1)*&
                        ( LOC_U(iv_idim, iv_iphase,1)-LOC_NU(iv_idim, iv_iphase,1))
                end forall 
                ! NDOTQNEW = NDOTQ + matmul(CVNORMX_ALL(:, GI), UGI_COEF_ELE_ALL(:, :,1)*&
                !     ( LOC_U(:,:,1)-LOC_NU(:,:,1)))
            else
                do iv_iphase = 1,final_phase
                    NDOTQNEW(iv_iphase) = NDOTQ(iv_iphase) + dot_product(matmul( CVNORMX_ALL(:, GI), UGI_COEF_ELE_ALL(:, iv_iphase,:)*&
                        ( LOC_U(:,iv_iphase,:)-LOC_NU(:,iv_iphase,:))), CV_funs%sufen( :, GI ))
                end do
            end if
                ! IF( between_elements) THEN
                !     ! We have a discontinuity between elements so integrate along the face...
                !     DO iv_u_skloc = 1, Mdims%u_snloc
                !         iv_u_kloc = U_SLOC2LOC(iv_u_skloc)
                !         iv_u_kloc2 = U_OTHER_LOC( iv_u_kloc )
                !         DO iv_idim = 1, Mdims%ndim
                !             NDOTQNEW=NDOTQNEW + CV_funs%sufen( iv_u_kloc, GI ) * UGI_COEF_ELE2_ALL(iv_idim, :,iv_u_kloc2) &
                !                 * ( LOC2_U(iv_idim,:, iv_u_kloc ) - LOC2_NU(iv_idim,:,iv_u_kloc ) ) * CVNORMX_ALL(iv_idim, GI)
                !         END DO
                !     END DO
                ! END IF
            RETURN
        END SUBROUTINE GET_INT_VEL_POROUS_VEL

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
                ! for integrating just on one side...
                ICOUNT = 0
                DO COUNT = Mspars%CT%fin( CV_NODJ ), Mspars%CT%fin( CV_NODJ + 1 ) - 1
                    IF ( Mspars%CT%col( COUNT ) == U_NODK ) THEN
                        ICOUNT = COUNT
                        EXIT
                    END IF
                END DO
                ICOUNT_KLOC( U_KLOC ) = ICOUNT
            END DO
            ! IF ( between_elements ) THEN
            !     DO U_KLOC =  1, Mdims%u_nloc
            !         U_NODK = ndgln%u( ( ELE2 - 1 ) * Mdims%u_nloc + U_KLOC )
            !         JCOUNT = 0
            !         DO COUNT = Mspars%CT%fin( CV_NODI ), Mspars%CT%fin( CV_NODI + 1 ) - 1
            !             IF ( Mspars%CT%col( COUNT ) == U_NODK ) THEN
            !                 JCOUNT = COUNT
            !                 EXIT
            !             END IF
            !         END DO
            !         JCOUNT_KLOC2( U_KLOC ) = JCOUNT
        !             ! for integrating just on one side...
        !             ICOUNT = 0
        !             DO COUNT = Mspars%CT%fin( CV_NODJ ), Mspars%CT%fin( CV_NODJ + 1 ) - 1
        !                 IF ( Mspars%CT%col( COUNT ) == U_NODK ) THEN
        !                     ICOUNT = COUNT
        !                     EXIT
        !                 END IF
        !             END DO
        !             ICOUNT_KLOC2( U_KLOC ) = ICOUNT
            !     END DO
            ! END IF ! endof IF ( between_elements ) THEN
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
                    ! for integrating just on one side...
                    ICOUNT = 0
                    DO COUNT = Mspars%C%fin( U_NODK ), Mspars%C%fin( U_NODK + 1 ) - 1
                        IF ( Mspars%C%col( COUNT ) == CV_NODJ ) THEN
                            ICOUNT = COUNT
                            EXIT
                        END IF
                    END DO
                    C_ICOUNT_KLOC( U_KLOC ) = ICOUNT
                END DO
                ! IF ( between_elements ) THEN
                !     DO U_KLOC =  1, Mdims%u_nloc
                !         U_NODK = ndgln%u( ( ELE2 - 1 ) * Mdims%u_nloc + U_KLOC )
                !         JCOUNT = 0
                !         DO COUNT = Mspars%C%fin( U_NODK ), Mspars%C%fin( U_NODK + 1 ) - 1
                !             IF ( Mspars%C%col( COUNT ) == CV_NODI ) THEN
                !                 JCOUNT = COUNT
                !                 EXIT
                !             END IF
                !         END DO
                !         C_JCOUNT_KLOC2( U_KLOC ) = JCOUNT
                !             ! for integrating just on one side...
                !             ICOUNT = 0
                !             DO COUNT = Mspars%C%fin( U_NODK ), Mspars%C%fin( U_NODK + 1 ) - 1
                !                 IF ( Mspars%C%col( COUNT ) == CV_NODJ ) THEN
                !                     ICOUNT = COUNT
                !                     EXIT
                !                 END IF
                !             END DO
                !             C_ICOUNT_KLOC2( U_KLOC ) = ICOUNT
                !     END DO
                ! END IF ! endof IF ( between_elements ) THEN
            ENDIF ! ENDOF IF(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV) THEN
        end subroutine get_neigbouring_lists

        subroutine mass_conservation_check_and_outfluxes(calculate_mass_delta, outfluxes, DEN_ALL, flag)
            ! Subroutine to calculate the integrated flux across a boundary with the specified surface_ids given that the massflux has been already stored elsewhere
            !also used to calculate mass conservation

            !The check is done as follows:
            !First the mass at the beginning of the time-level is calculated in the reservoir
            !Second the input and output throughtout all the boundaries is calculated
            !Third the mass at the current saturation fixed point iteration is calculated within the domain
            !Finally, due to mass conservation and given a constant total volume of the domain the mass check is done as follows:
            ! Mass_error = abs (Mass_present_resv - Mass_initial_resv + Difference_Mass_io )/Maximum_flux
            ! in this way the mass error is normalise based on the mass flux.
            implicit none
            integer, intent(in) :: flag
            type (multi_outfluxes), intent(inout) :: outfluxes
            real, dimension(:,:), intent(inout) :: calculate_mass_delta
            REAL, DIMENSION( :, : ), intent( in) :: DEN_ALL!Density including the boussinesq aprox.
            !Local variables
            integer :: iphase, k
            real :: tmp1, tmp2, tmp3, maxflux
            type(vector_field), pointer :: Por
            ! After looping over all elements, calculate the mass change inside the domain normalised to the mass inside the domain at t=t-1
            ! Difference in Total mass

            select case (flag)
                case (1)
                    !Initialise values
                    if (first_nonlinear_time_step ) then
                        calculate_mass_delta(:,1) = 0.0 ! reinitialise
                        call calculate_internal_volume( packed_state, Mdims, Mass_ELE, &
                            calculate_mass_delta(1:Mdims%n_in_pres,1) , ndgln%cv, DEN_ALL)
                        !DISABLED AS IT DOES NOT WORK WELL AND IT DOES ACCOUNT FOR A VERY TINY FRACTION OF THE OVERALL MASS
                       ! if (Mdims%npres >1)then!consider as well the pipes
                       !     call calculate_internal_volume( packed_state, Mdims, pipes_aux%MASS_PIPE, &
                       !         calculate_mass_delta(:,1) , ndgln%cv, eles_with_pipe)
                       ! end if
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
                    calculate_mass_internal = 0.
                    call calculate_internal_volume( packed_state, Mdims, Mass_ELE, &
                        calculate_mass_internal(1:Mdims%n_in_pres) , ndgln%cv, DEN_ALL)
                    !DISABLED AS IT DOES NOT WORK WELL AND IT DOES ACCOUNT FOR A VERY TINY FRACTION OF THE OVERALL MASS
                   ! if (Mdims%npres >1) then!consider as well the pipes
                   !     call calculate_internal_volume( packed_state, Mdims, pipes_aux%MASS_PIPE, &
                   !         calculate_mass_internal(:) , ndgln%cv, eles_with_pipe)
                   ! end if

                    !Loop over nphases - 1
                    calculate_mass_delta(1,2) = 0.

                    tmp1 = sum(calculate_mass_internal)!<=mass inside the domain
                    tmp2 = sum(calculate_mass_delta(:,1))!<= mass inside the domain at the beginning of the time-level
                    tmp3 = sum(bcs_outfluxes(:, :, 0))*dt!<= mass phase i across all the boundaries
                    if (isparallel()) then
                        call allsum(tmp1)
                        call allsum(tmp2)
                        call allsum(tmp3)
                    end if
                    !Obtain the maximum flux and use this to normalise
                    maxflux = 0.
                    do iphase = 1, Mdims%nphase
                        maxflux= max(maxflux, abs(sum(bcs_outfluxes(iphase, :, 0)))*dt)
                    end do

                    !Calculate possible mass creation inside the code
                    if (tmp2 > 1d-8) then
                        calculate_mass_delta(1,2) = abs( tmp1 - tmp2 + tmp3 ) / tmp2 !We normise by the maximum flux
                    else
                        calculate_mass_delta(1,2) = 0.
                    end if

                    !If calculate outfluxes then do it now
                    if (outfluxes%calculate_flux) then
                        !Retrieve only the values of bcs_outfluxes we are interested in
                        do k = 1, size(outfluxes%outlet_id)
                            do iphase = 1, Mdims%nphase
                                outfluxes%totout(iphase, k) = sum(bcs_outfluxes( iphase, :, k))
                            end do
                        end do
                        ! Having finished loop over elements etc. Pass the total flux across all boundaries to the global variable totout
                        if (isParallel()) then
                            do k = 1,size(outfluxes%outlet_id)
                                ! Ensure all processors have the correct value of totout for parallel runs
                                do iphase = 1, Mdims%nphase
                                    call allsum(outfluxes%totout(iphase, k))
                                end do
                            end do
                        end if
                    end if
            end select


        end subroutine mass_conservation_check_and_outfluxes

    END SUBROUTINE CV_ASSEMB

    !>     This subroutine counts then number of faces in the control volume space
    function CV_count_faces( Mdims, CV_ELE_TYPE, CV_GIdims) result(global_face)
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




    
    !> We are on the boundary or next to another element. Determine CV_OTHER_LOC
    !> CVFEM_ON_FACE(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
    !> Look for these nodes on the other elements.
    !> ELE2=0 also when we are between elements but are trying to integrate across
    !> the middle of a CV.
    SUBROUTINE FIND_OTHER_SIDE( CV_OTHER_LOC, CV_NLOC, U_NLOC,  &
        MAT_OTHER_LOC, INTEGRAT_AT_GI, &
        X_NLOC, XU_NLOC, X_NDGLN, XU_NDGLN, &
        CV_SNLOC, CVFEM_ON_FACE, X_SHARE, ELE, ELE2,  &
        FINELE, COLELE)!, DISTCONTINUOUS_METHOD )
        IMPLICIT NONE
        INTEGER, intent( in ) :: CV_NLOC, U_NLOC, X_NLOC, XU_NLOC, &
            &                   CV_SNLOC, ELE
        INTEGER, DIMENSION( : ), intent( in ) :: X_NDGLN, XU_NDGLN
        LOGICAL, DIMENSION( : ), intent( in ) :: CVFEM_ON_FACE
        INTEGER, DIMENSION( : ), intent( in ) :: FINELE
        INTEGER, DIMENSION( : ), intent( in ) :: COLELE

        INTEGER, DIMENSION( : ), intent( inout ) :: CV_OTHER_LOC, MAT_OTHER_LOC
        LOGICAL, DIMENSION( : ), intent( inout ) :: X_SHARE
        INTEGER, intent( inout ) :: ELE2
        LOGICAL, intent( inout ) :: INTEGRAT_AT_GI
        ! LOGICAL, intent( in ) :: DISTCONTINUOUS_METHOD
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


        ! Quit because there is no work to do here...
        ! IF(.NOT.DISTCONTINUOUS_METHOD) THEN
        IF(ELE2.NE.0) THEN ! this is not on the boundary of the domain.
            INTEGRAT_AT_GI=.FALSE.
            RETURN
        ENDIF
        ! ENDIF


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

            MAT_OTHER_LOC = CV_OTHER_LOC
        ELSE
            CV_OTHER_LOC = 0
            MAT_OTHER_LOC = 0
        END IF

        RETURN

    END SUBROUTINE FIND_OTHER_SIDE

    !---------------------------------------------------------------------------
    !> @author Chris Pain, Pablo Salinas
    !> @brief Calls to generate the transport equation for the saturation. Embeded an FPI with backtracking method is uncluded
    !> Subroutine description:
    !> (1) determine FEMT (finite element wise) etc
    !>     from T (control volume wise)
    !> (2) (optional) calculate psi_int (area) and
    !>     psi_ave (barycentre) over each CV
    !>
    !>@param packed_state  ! local state data
    !>@param fempsi   ! finite element field data
    !>@param psi      ! finite volume field data
    !>@param Mdims        ! dimension data
    !>@param CV_GIdims ! gauss integer dimension data
    !>@param CV_funs   ! control volume shape function data
    !>@param Mspars       ! sparsity data
    !>@param ndgln             ! global numbering data
    !>@param igetct     ! whether to get CT matrix
    !>@param X              ! coordinates of the elements
    !>@param mass_ele      ! finite element mass
    !>@param mass_mn_pres  ! ??
    !>@param tracer field to be projected
    !>@param psi_int ! control volume area
    !>@param psi_ave ! control volume barycentre
    SUBROUTINE PROJ_CV_TO_FEM(packed_state, &
        fempsi, psi, &
        Mdims, CV_GIdims, CV_funs, Mspars, ndgln, &
        igetct, X, mass_ele, mass_mn_pres, &
        tracer, psi_ave, psi_int)


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
        type(csr_sparsity), pointer :: sparsity
        logical :: do_not_project,  is_to_update
        character(len=option_path_len) :: option_path

        !---------------------------------
        ! initialisation and allocation
        !---------------------------------
        !Currently hard-coded. This is not used for porous_media but it is used otherwise
        do_not_project = .true.
        is_to_update = .not.associated(CV_funs%CV2FE%refcount)!I think this is only true after adapt and at the beginning
        if (.not. do_not_project) then
            do it=1,size(fempsi)
                call zero(fempsi(it)%ptr)
                call allocate(fempsi_rhs(it),psi(it)%ptr%mesh,"RHS",dim=psi(it)%ptr%dim)
                call zero(fempsi_rhs(it))
                ! call halo_update(psi(it)%ptr)!Fields should get here already updated
            end do
        end if
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
            !sprint_to_do we should not need to allocate this field unless .not. do_not_project
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
            ! calculate detwei,RA,NX,NY,NZ for the ith element
            call DETNLXR(iele, X, ndgln%x, CV_funs%cvweight, CV_funs%CVFEN, CV_funs%CVFENLX_ALL, DevFuns)

            mass_ele(iele) = DevFuns%volume

            Loop_CV_iLoc: do cv_iloc = 1, Mdims%cv_nloc
                cv_nodi = ndgln%cv((iele-1)*Mdims%cv_nloc+cv_iloc)

                Loop_CV_jLoc: do cv_jloc = 1, Mdims%cv_nloc
                    cv_nodj = ndgln%cv((iele-1)*Mdims%cv_nloc+cv_jloc)

                    nn = 0.0; nm = 0.0; mn = 0.0; mm = 0.0
                    do cv_gi = 1, CV_GIdims%cv_ngi
                        mn = mn+CV_funs%CVN(cv_iloc,cv_gi)*CV_funs%CVFEN(cv_jloc,cv_gi)*DevFuns%detwei(cv_gi)
                        if (is_to_update)  mm = mm+CV_funs%CVN(cv_iloc,cv_gi)*CV_funs%CVN(cv_jloc,cv_gi)*DevFuns%detwei(cv_gi)
                        if(is_to_update.and. .not. do_not_project) nn = nn+CV_funs%CVFEN(cv_iloc,cv_gi)*CV_funs%CVFEN(cv_jloc,cv_gi)*DevFuns%detwei(cv_gi)
                        if (.not. do_not_project) nm = nm+CV_funs%CVFEN(cv_iloc,cv_gi)*CV_funs%CVN(cv_jloc,cv_gi)*DevFuns%detwei(cv_gi)
                    end do


                        if(is_to_update.and. .not. do_not_project) call addto(CV_funs%CV2FE,1,1,cv_nodi,cv_nodj,nn)
                      if (.not. do_not_project) then
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
            call invert(cv_mass)
            do it = 1, size(psi_ave)
                call scale(psi_ave(it)%ptr,cv_mass)
            end do
        end if

        ! deallocation
        call deallocate_multi_dev_shape_funs(DevFuns)

        if (.not. do_not_project) then 
            do it = 1,size(fempsi_rhs)
                call deallocate(fempsi_rhs(it))
            end do
        end if
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


    !> calculates derivatives of vector fields
    SUBROUTINE DG_DERIVS_ALL1( FEMT, &
        DTX_ELE, &
        NDIM, NPHASE, NCOMP, TOTELE, CV_NDGLN, & ! ncomp = ndim here
        XCV_NDGLN, X_NLOC, X_NDGLN,&
        CV_NGI, CV_NLOC, CVWEIGHT, &
        N, NLX, NLY, NLZ, &
        X_N, X_NLX, X_NLY, X_NLZ, &
        X_NONODS, X, Y, Z, &
        NFACE, FACE_ELE, CV_SLOCLIST, X_SLOCLIST, CV_SNLOC, X_SNLOC, WIC_T_BC, SUF_T_BC, &
        SBCVNGI, SBCVFEN, SBWEIGH, &
        X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY, get_gradU, state, P0Mesh)

        ! calculates derivatives of vector fields
        IMPLICIT NONE

        INTEGER, intent( in ) :: NDIM, NPHASE, NCOMP, TOTELE, X_NLOC, CV_NGI, CV_NLOC, &
            X_NONODS, CV_SNLOC, X_SNLOC, SBCVNGI, NFACE
        REAL, DIMENSION( :, :, : ), intent( in ) :: FEMT
        REAL, DIMENSION( :, :, :, :, : ), intent( inout ) :: DTX_ELE
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
        logical, optional :: P0Mesh !> This is for when using the P0DGP1 element pair to substitute the method to find neighbours
        ! Local variables
        REAL, DIMENSION( CV_NLOC, CV_NLOC, TOTELE ) :: MASELE
        REAL, DIMENSION( NDIM, NCOMP, NPHASE, CV_NLOC, TOTELE ) :: VTX_ELE
        LOGICAL, dimension(NCOMP, NPHASE) :: APPLYBC
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

        DTX_ELE = 0.0 ;
        MASELE = 0.0 ;VTX_ELE = 0.0

        Loop_Elements1: DO ELE = 1, TOTELE

            ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
            CALL DETNLXR_PLUS_U( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
                X_NLOC, X_NLOC, CV_NGI, &
                X_N, X_NLX, X_NLY, X_NLZ, CVWEIGHT, DETWEI, RA, VOLUME, NDIM == 1, NDIM == 3, .false., &
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
                    if (present_and_true(P0Mesh)) then
                      ILOC_OTHER_SIDE = 1
                    else
                      DO CV_SILOC = 1, CV_SNLOC
                          CV_ILOC = SLOC2LOC( CV_SILOC )
                          CV_INOD = XCV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )

                          DO CV_ILOC2 = 1, CV_NLOC
                              CV_INOD2 = XCV_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_ILOC2 )

                              IF( CV_INOD2 == CV_INOD ) ILOC_OTHER_SIDE( CV_SILOC ) = CV_ILOC2
                          END DO
                      END DO
                    end if
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
                                        END IF
                                    ELSE
                                        VTX_ELE(:, ICOMP, IPHASE, CV_ILOC, ELE ) = &
                                            VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                                            - VLM_NORX(:) * 0.5 * ( FEMT( ICOMP, IPHASE, CV_NODJ ) - FEMT( ICOMP, IPHASE, CV_NODJ2 )  )
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


        ewrite(3,*)'about to leave DG_DERIVS'

        RETURN

    END SUBROUTINE DG_DERIVS_ALL1

    !> Computes the derivatives of vector fields
    SUBROUTINE DG_DERIVS_ALL2( FEMT, &
        DTX_ELE, &
        NDIM, NPHASE, TOTELE, CV_NDGLN, &
        XCV_NDGLN, X_NLOC, X_NDGLN,&
        CV_NGI, CV_NLOC, CVWEIGHT, &
        N, NLX, NLY, NLZ, &
        X_N, X_NLX, X_NLY, X_NLZ, &
        X_NONODS, X, Y, Z, &
        NFACE, FACE_ELE, CV_SLOCLIST, X_SLOCLIST, CV_SNLOC, X_SNLOC, WIC_T_BC, SUF_T_BC, &
        SBCVNGI, SBCVFEN, SBWEIGH, &
        X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY, P0Mesh)

        ! determine FEMT (finite element wise) etc from T (control volume wise)
        IMPLICIT NONE

        INTEGER, intent( in ) :: NDIM, NPHASE, TOTELE, X_NLOC, CV_NGI, CV_NLOC, &
            &                   X_NONODS, CV_SNLOC, X_SNLOC, SBCVNGI, NFACE
        REAL, DIMENSION( :, : ), intent( in ) :: FEMT
        REAL, DIMENSION( :, :, :, : ), intent( inout ) :: DTX_ELE
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
        logical, optional :: P0Mesh !> This is for when using the P0DGP1 element pair to substitute the method to find neighbours
        ! Local variables
        REAL, DIMENSION( CV_NLOC, CV_NLOC, TOTELE ) :: MASELE
        REAL, DIMENSION( NDIM, NPHASE, CV_NLOC, TOTELE ) :: VTX_ELE
        LOGICAL :: D1, D3, APPLYBC( NPHASE )
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

        DTX_ELE = 0.0 
        MASELE = 0.0 ; VTX_ELE = 0.0

        !DCYL = .FALSE.
        Loop_Elements1: DO ELE = 1, TOTELE

            ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
            CALL DETNLXR_PLUS_U( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
                X_NLOC, X_NLOC, CV_NGI, &
                X_N, X_NLX, X_NLY, X_NLZ, CVWEIGHT, DETWEI, RA, VOLUME, NDIM == 1, NDIM == 3, .false., &
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
                    END DO

                END DO Loop_CV_JLOC

            END DO Loop_CV_ILOC

        END DO Loop_Elements1


        Loop_Elements2: DO ELE = 1, TOTELE

            Between_Elements_And_Boundary: DO IFACE = 1, NFACE
                ELE2 = FACE_ELE( IFACE, ELE )
                SELE2 = MAX( 0, - ELE2 )
                ELE2 = MAX( 0, + ELE2 )

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

                ! Calculation of the norms and determinant
                CALL DGSDETNXLOC2(X_SNLOC, SBCVNGI, &
                    XSL( 1, : ), XSL( 2, : ), XSL( 3, : ), &
                    X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY, SBWEIGH, SDETWE, SAREA, &
                    (NDIM==1), (NDIM==3), (NDIM==-2), &
                    SNORMXN( 1, : ), SNORMXN( 2, : ), SNORMXN( 3, : ), &
                    NORMX( 1 ), NORMX( 2 ), NORMX( 3 ) )

                IF ( SELE2 == 0 ) THEN
                    ! Calculate the nodes on the other side of the face:
                    if (present_and_true(P0Mesh)) then
                      ILOC_OTHER_SIDE = 1
                    else
                      DO CV_SILOC = 1, CV_SNLOC
                          CV_ILOC = SLOC2LOC( CV_SILOC )
                          CV_INOD = XCV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )

                          DO CV_ILOC2 = 1, CV_NLOC
                              CV_INOD2 = XCV_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_ILOC2 )

                              IF( CV_INOD2 == CV_INOD ) ILOC_OTHER_SIDE( CV_SILOC ) = CV_ILOC2
                          END DO
                      END DO
                    end if
                    APPLYBC = (ELE /= ELE2) .AND. (ELE2 /= 0)

                ELSE
                    APPLYBC = ( WIC_T_BC( 1, 1:nphase, SELE2 ) == WIC_T_BC_DIRICHLET )
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

                                IF ( SELE2 /= 0 ) THEN
                                    IF ( WIC_T_BC(1,  IPHASE, SELE2 ) == WIC_T_BC_DIRICHLET ) THEN

                                        RTBC = SUF_T_BC( 1,IPHASE, CV_SJLOC2+ CV_SNLOC*( SELE2-1) )

                                        VTX_ELE( :, IPHASE, CV_ILOC, ELE ) = VTX_ELE( :, IPHASE, CV_ILOC, ELE ) &
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

            END FORALL

        END DO Loop_Elements3

        ewrite(3,*)'about to leave DG_DERIVS'

    END SUBROUTINE DG_DERIVS_ALL2


    !> This sub calculates the effective diffusion coefficientd DIFF_COEF_DIVDX
    !> based on a non-linear method and a non-oscillating scheme.
    !> It requires the derivatives of the field obtained using DG_DERIVS_ALL
    !> @ref DG_DERIVS_ALL
    SUBROUTINE DIFFUS_CAL_COEFF(DIFF_COEF_DIVDX,  &
        CV_NLOC, MAT_NLOC, NPHASE, MAT_NDGLN, &
        SMATFEN, SCVFEN, GI, NDIM, TDIFFUSION, &
        HDC, &
        T_CV_NODJ, T_CV_NODI, &
        ELE, ELE2, CVNORMX_ALL,  &
        LOC_DTX_ELE_ALL, LOC2_DTX_ELE_ALL, &
        LOC_WIC_T_BC, CV_OTHER_LOC, MAT_OTHER_LOC, CV_SNLOC, CV_SLOC2LOC, &
        on_domain_boundary, between_elements )
        IMPLICIT NONE
        INTEGER, intent( in ) :: CV_NLOC, MAT_NLOC, NPHASE, &
            &                    GI, NDIM, ELE, ELE2, &
            &                    CV_SNLOC
        REAL, intent( in ) :: HDC
        LOGICAL, intent( in ) :: on_domain_boundary, between_elements
        REAL, DIMENSION( NPHASE ), intent( in ) :: T_CV_NODJ, T_CV_NODI
        REAL, DIMENSION( NPHASE ), intent( inout ) :: DIFF_COEF_DIVDX
        INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: CV_SLOC2LOC
        INTEGER, DIMENSION( NPHASE ), intent( in ) :: LOC_WIC_T_BC
        INTEGER, DIMENSION( : ), intent( in ) :: CV_OTHER_LOC
        INTEGER, DIMENSION( : ), intent( in ) :: MAT_OTHER_LOC
        REAL, DIMENSION( :, : ), intent( in ) :: SMATFEN
        REAL, DIMENSION( :, : ), intent( in ) :: SCVFEN
        REAL, DIMENSION( :, :, :, : ), intent( in ) :: TDIFFUSION
        REAL, DIMENSION( NDIM, NPHASE, CV_NLOC ), intent( in ) :: LOC_DTX_ELE_ALL, LOC2_DTX_ELE_ALL
        REAL, DIMENSION( : ), intent( in ) :: CVNORMX_ALL

        ! local variables

        ! DIFF_MIN_FRAC is the fraction of the standard diffusion coefficient to use
        ! in the non-linear diffusion scheme. DIFF_MAX_FRAC is the maximum fraction.
        REAL, PARAMETER :: DIFF_MIN_FRAC = 0.05, DIFF_MAX_FRAC = 20.0
        INTEGER :: CV_KLOC, CV_KLOC2, MAT_KLOC, MAT_KLOC2, MAT_NODK, MAT_NODK2, IPHASE, CV_SKLOC, idim, jdim
        LOGICAL :: ZER_DIFF
        REAL :: COEF
        real, dimension(ndim) :: vCoef
        REAL, DIMENSION ( NDIM, NPHASE ) :: DTDX_GI_ALL, DTDX_GI2_ALL
        REAL, DIMENSION ( NPHASE ) :: N_DOT_DKDT_ALL, N_DOT_DKDT2_ALL
        REAL, DIMENSION ( NPHASE ) :: DIFF_STAND_DIVDX_ALL, DIFF_STAND_DIVDX2_ALL
        REAL, DIMENSION ( NDIM, NDIM, NPHASE ) :: DIFF_GI, DIFF_GI2

        ZER_DIFF = .FALSE.
        IF ( on_domain_boundary ) ZER_DIFF = ANY ( LOC_WIC_T_BC /= WIC_T_BC_DIRICHLET )

        Cond_ZerDiff: IF ( ZER_DIFF ) THEN

            DIFF_COEF_DIVDX = 0.0

        ELSE

            DTDX_GI_ALL = 0.0 
            forall (cv_kloc = 1:cv_nloc, idim = 1:ndim, iphase =1:nphase)
                DTDX_GI_ALL(idim, iphase) = DTDX_GI_ALL(idim, iphase) + SCVFEN( CV_KLOC, GI ) * LOC_DTX_ELE_ALL( idim, iphase, CV_KLOC )
            end forall

            DIFF_GI = 0.0
            DO MAT_KLOC = 1, MAT_NLOC
                MAT_NODK = MAT_NDGLN( ( ELE - 1 ) * MAT_NLOC + MAT_KLOC )
                forall (iphase = 1:nphase, idim = 1:ndim, jdim =1:ndim)
                    DIFF_GI( idim, jdim, IPHASE ) = DIFF_GI( idim, jdim, IPHASE ) &
                        + SMATFEN( MAT_KLOC, GI ) * TDIFFUSION( MAT_NODK, idim, jdim, IPHASE )
                end forall
            END DO
            DIFF_GI = MAX( 0.0, DIFF_GI )

            !Projection of the tensorial diffusion coefficient times derivatives
            N_DOT_DKDT_ALL = 0; DIFF_STAND_DIVDX_ALL = 0.
            forall (idim = 1:ndim, jdim = 1:ndim, iphase =1:nphase)
                N_DOT_DKDT_ALL( IPHASE ) = N_DOT_DKDT_ALL( IPHASE ) + CVNORMX_ALL(idim) * DIFF_GI( idim, jdim, IPHASE ) * DTDX_GI_ALL( jdim, IPHASE )
                DIFF_STAND_DIVDX_ALL( IPHASE ) = DIFF_STAND_DIVDX_ALL( IPHASE ) +  CVNORMX_ALL(idim) * DIFF_GI( idim, jdim, IPHASE ) * CVNORMX_ALL(jdim)/ HDC
            end forall

            ! Conditional_MAT_DISOPT_ELE2: IF ( between_elements ) THEN
            !     DTDX_GI2_ALL = 0.0 ;
            !     DO CV_SKLOC=1,CV_SNLOC
            !         CV_KLOC=CV_SLOC2LOC( CV_SKLOC )
            !         CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
            !         DTDX_GI2_ALL = DTDX_GI2_ALL + SCVFEN( CV_KLOC, GI ) *  LOC2_DTX_ELE_ALL( :, :, CV_KLOC )
            !     END DO

            !     DIFF_GI2 = 0.0
            !     DO MAT_KLOC = 1, MAT_NLOC
            !         MAT_KLOC2 = MAT_OTHER_LOC( MAT_KLOC )
            !         IF ( MAT_KLOC2 /= 0 ) THEN
            !             MAT_NODK2 = MAT_NDGLN( ( ELE2 - 1 ) * MAT_NLOC + MAT_KLOC2 )
            !             DO IPHASE = 1, NPHASE
            !                 DIFF_GI2( :, :, IPHASE ) = DIFF_GI2( :, :, IPHASE ) &
            !                     + SMATFEN( MAT_KLOC, GI ) * TDIFFUSION( MAT_NODK2, :, :, IPHASE )
            !             END DO
            !         END IF
            !     END DO

            !     DO IPHASE = 1, NPHASE
            !         N_DOT_DKDT2_ALL( IPHASE ) = DOT_PRODUCT( CVNORMX_ALL, MATMUL( DIFF_GI2( :, :, IPHASE ), DTDX_GI2_ALL( :, IPHASE ) ) )

            !         COEF = DOT_PRODUCT( CVNORMX_ALL, MATMUL( DIFF_GI2( :, :, IPHASE ), CVNORMX_ALL ) )
            !         DIFF_STAND_DIVDX2_ALL( IPHASE ) = COEF  /HDC
            !     END DO

            !     N_DOT_DKDT_ALL = 0.5 * ( N_DOT_DKDT_ALL + N_DOT_DKDT2_ALL )

            !     DIFF_STAND_DIVDX_ALL = MIN( DIFF_STAND_DIVDX_ALL, DIFF_STAND_DIVDX2_ALL )

            ! END IF Conditional_MAT_DISOPT_ELE2

            DIFF_COEF_DIVDX = MAX( DIFF_MIN_FRAC * DIFF_STAND_DIVDX_ALL, N_DOT_DKDT_ALL / TOLFUN_MANY( T_CV_NODJ - T_CV_NODI ) )
            DIFF_COEF_DIVDX = MIN( DIFF_MAX_FRAC * DIFF_STAND_DIVDX_ALL, DIFF_COEF_DIVDX )

        END IF Cond_ZerDiff

        !    IF ( SELE /= 0 ) THEN
        IF ( on_domain_boundary ) THEN
            DO IPHASE=1,NPHASE
                IF( LOC_WIC_T_BC( IPHASE ) /= WIC_T_BC_DIRICHLET ) THEN
                    DIFF_COEF_DIVDX( IPHASE ) = 0.0
                ENDIF
            END DO
        ENDIF



        RETURN

    END SUBROUTINE DIFFUS_CAL_COEFF

    !> Form approximate surface normal (NORMX_ALL(1),NORMX_ALL(2),NORMX_ALL(3))
    SUBROUTINE DGSIMPLNORM_ALL( NLOC, SNLOC, NDIM,  &
        XL_ALL, XSL_ALL, NORMX_ALL )
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



    !> Form approximate surface normal (NORMX,NORMY,NORMZ)
    SUBROUTINE DGSIMPLNORM( ELE, SILOC2ILOC, NLOC, SNLOC, XONDGL, &
        X, Y, Z, NORMX, NORMY, NORMZ )
        !sprint_to_do where this is being called use the new one with the new memory
        !and then remove
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



    
    !> Calculate surface element, surface control volume: SELE, CV_SILOC, U_SLOC2LOC, CV_SLOC2LOC for a face on the
    !> boundary of the domain
    SUBROUTINE CALC_SELE( ELE, ELE3, SELE, CV_SILOC, CV_ILOC, U_SLOC2LOC, CV_SLOC2LOC, &
        FACE_ELE, gi, CV_funs, Mdims, CV_GIdims, &
        CV_NDGLN, U_NDGLN, CV_SNDGLN, U_SNDGLN )
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




    !> This subroutine caculates the discretised cty eqn acting on the velocities i.e. Mmat%CT, Mmat%CT_RHS
    !> It also computes the gradient matrix using the DCVFE method
    SUBROUTINE PUT_IN_CT_RHS( GET_C_IN_CV_ADVDIF_AND_CALC_C_CV, ct_rhs_phase_cv_nodi, ct_rhs_phase_cv_nodj, &
        final_phase, Mdims, CV_funs, ndgln, Mmat, GI, between_elements, on_domain_boundary, &
        ELE, ELE2, SELE, HDC, MASS_ELE, JCOUNT_KLOC, JCOUNT_KLOC2, ICOUNT_KLOC, ICOUNT_KLOC2, &
        C_JCOUNT_KLOC, C_JCOUNT_KLOC2, C_ICOUNT_KLOC, C_ICOUNT_KLOC2, &
        U_SLOC2LOC, CV_SLOC2LOC,  &
        SCVDETWEI, CVNORMX_ALL, DEN_ALL, CV_NODI, CV_NODJ, &
        WIC_U_BC_ALL, WIC_P_BC_ALL,SUF_P_BC_ALL,&
        UGI_COEF_ELE_ALL,  &
        UGI_COEF_ELE2_ALL,  &
        NDOTQ, LIMT, LIMDT, LIMT_HAT, LIMT2, &
        NDOTQ_HAT, &
        integrate_other_side_and_not_boundary, &
        loc_u, THETA_VEL, &
        UDGI_IMP_ALL, RCON, RCON_J, NDOTQ_IMP, X_ALL, SUF_D_BC_ALL, gravty) !<= local memory sent down for speed...
        IMPLICIT NONE
        ! IF more_in_ct THEN PUT AS MUCH AS POSSIBLE INTO Mmat%CT MATRIX
        !    LOGICAL, PARAMETER :: more_in_ct=.false.
        INTEGER, intent( in ) :: GI, &
            CV_NODI, CV_NODJ, ELE, ELE2, SELE, final_phase
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_shape_funs), intent(in) :: CV_funs
        type(multi_ndgln), intent(in) :: ndgln
        type (multi_matrices), intent(inout) :: Mmat
        REAL, DIMENSION( :, :, : ), intent( in ) :: loc_u
        LOGICAL, intent( in ) :: integrate_other_side_and_not_boundary, between_elements, on_domain_boundary,&
            GET_C_IN_CV_ADVDIF_AND_CALC_C_CV
        INTEGER, DIMENSION( : ), intent( in ) :: JCOUNT_KLOC, JCOUNT_KLOC2, ICOUNT_KLOC, ICOUNT_KLOC2
        INTEGER, DIMENSION( : ), intent( in ) :: C_JCOUNT_KLOC, C_JCOUNT_KLOC2, C_ICOUNT_KLOC, C_ICOUNT_KLOC2
        INTEGER, DIMENSION( : ), intent( in ) :: U_SLOC2LOC, CV_SLOC2LOC
        REAL, DIMENSION( :, :, : ), intent( inout ) :: SUF_P_BC_ALL
        REAL, DIMENSION( : ), intent( inout ) :: ct_rhs_phase_cv_nodi, ct_rhs_phase_cv_nodj
        REAL, DIMENSION( :, :, : ), intent( in ) :: UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL
        REAL, DIMENSION( : ), intent( in ) :: SCVDETWEI, MASS_ELE
        REAL, DIMENSION( :, : ), intent( in ) :: CVNORMX_ALL
        REAL, DIMENSION( :, : ), intent( in ) :: DEN_ALL
        REAL, DIMENSION( : ), intent( in ) :: NDOTQ, LIMT, LIMDT, LIMT_HAT, limt2
        REAL, intent( in ) :: NDOTQ_HAT
        REAL, DIMENSION( : ), intent( in ) :: THETA_VEL
        integer, dimension(:,:,:) :: WIC_U_BC_ALL, WIC_P_BC_ALL
        ! LIMT_HAT is the normalised voln fraction
        REAL, intent( in ) :: HDC
        ! local memory sent down for speed...
        REAL,  DIMENSION( Mdims%ndim, final_phase ), intent( inout ) :: UDGI_IMP_ALL
        REAL,  DIMENSION( final_phase ), intent( inout ) :: RCON, RCON_J, NDOTQ_IMP
        !Variable to account for boundary conditions if using GET_C_IN_CV_ADVDIF_AND_CALC_C_CV
        real, dimension (Mdims%ndim, final_phase, Mdims%u_nloc ) :: Bound_ele_correct
        ! coordinates
        real, dimension(:,:), pointer :: X_ALL
        !density
        REAL, DIMENSION(:,:,: ), pointer :: SUF_D_BC_ALL
        !gravity
        REAL :: gravty
        ! Local variables...
        INTEGER :: U_KLOC, U_KLOC2, IDIM, &
            IPHASE, U_SKLOC, I, J, U_KKLOC, &
            u_iloc, u_siloc, count
        real :: Mass_corrector, auxR
        !Local variables for CV pressure bcs
        integer :: KPHASE, CV_SNODK, CV_SNODK_IPHA, CV_SKLOC
        real, dimension(:,:), allocatable :: SUF_SIG_DIAGTEN_BC_pha_GI

        !If using Mmat%C_CV prepare Bound_ele_correct and Bound_ele2_correct to correctly apply the BCs
        if (Mmat%CV_pressure) call introduce_C_CV_boundary_conditions(Bound_ele_correct)

        DO U_KLOC = 1, Mdims%u_nloc
            RCON = SCVDETWEI( GI ) * (  LIMT2 * LIMDT ) &
                * CV_funs%sufen( U_KLOC, GI ) / DEN_ALL( :, CV_NODI )
            DO IPHASE = 1,final_phase
                Mmat%CT( :, IPHASE, JCOUNT_KLOC( U_KLOC ) ) = Mmat%CT( :, IPHASE, JCOUNT_KLOC( U_KLOC ) ) &
                    + rcon(IPHASE) * UGI_COEF_ELE_ALL( :, IPHASE, U_KLOC ) * CVNORMX_ALL( :, GI )
            END DO
            IF(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV) THEN
                rcon = SCVDETWEI( GI ) * CV_funs%sufen( U_KLOC, GI )
                DO IPHASE=1,final_phase!Mdims%nphase
                    ! IF ( between_elements) THEN
                    !     ! bias the weighting towards bigger eles - works with 0.25 and 0.1 and not 0.01.
                    !     !This is to perform the average between two DG pressures (same mass => 0.5)
                    !     auxR = 0.25!<= original!0.5 seems similar, should try up to 2
                    !     Mass_corrector = (MASS_ELE( ELE2 ) + auxR * MASS_ELE( ELE ))/( (1.+auxR) *(MASS_ELE( ELE ) + MASS_ELE( ELE2 )))
                    !     !Mass_corrector = MASS_ELE( ELE2 )/(MASS_ELE( ELE2 ) + MASS_ELE( ELE ) )!<=this seems to work

                    !     Mmat%C_CV( :, IPHASE, C_JCOUNT_KLOC( U_KLOC ) ) &
                    !         = Mmat%C_CV( :, IPHASE, C_JCOUNT_KLOC( U_KLOC ) ) &
                    !         + rcon(IPHASE) * CVNORMX_ALL( :, GI ) * Mass_corrector
                    ! else
                        Mmat%C_CV( :, IPHASE, C_JCOUNT_KLOC( U_KLOC ) ) &
                            = Mmat%C_CV( :, IPHASE, C_JCOUNT_KLOC( U_KLOC ) ) &
                            + rcon(IPHASE) * CVNORMX_ALL( :, GI ) * Bound_ele_correct(:, IPHASE, U_KLOC)
                    ! endif
                END DO
            ENDIF
            ! flux from the other side (change of sign because normal is -ve)...
            if ( integrate_other_side_and_not_boundary ) then
                RCON_J = SCVDETWEI( GI ) * ( LIMT2* LIMDT )  &
                    * CV_funs%sufen( U_KLOC, GI ) / DEN_ALL( :, CV_NODJ )
                DO IPHASE = 1,final_phase
                    Mmat%CT( :, IPHASE, ICOUNT_KLOC( U_KLOC ) ) = Mmat%CT( :, IPHASE, ICOUNT_KLOC( U_KLOC ) ) &
                        - RCON_J(IPHASE) * UGI_COEF_ELE_ALL( :, IPHASE, U_KLOC ) * CVNORMX_ALL( :, GI )
                END DO
                IF(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV) THEN
                    RCON_J = SCVDETWEI( GI ) * CV_funs%sufen( U_KLOC, GI )
                    DO IPHASE=1,final_phase!Mdims%nphase
                        ! IF ( between_elements ) THEN
                        !     Mmat%C_CV( :, IPHASE, C_ICOUNT_KLOC( U_KLOC ) ) &
                        !         = Mmat%C_CV( :, IPHASE, C_ICOUNT_KLOC( U_KLOC ) ) &
                        !         - RCON_J(IPHASE) * CVNORMX_ALL( :, GI )* Mass_corrector!(1.- Mass_corrector)
                        ! else
                            Mmat%C_CV( :, IPHASE, C_ICOUNT_KLOC( U_KLOC ) ) &
                                = Mmat%C_CV( :, IPHASE, C_ICOUNT_KLOC( U_KLOC ) ) &
                                - RCON_J(IPHASE) * CVNORMX_ALL( :, GI )* Bound_ele_correct(:, IPHASE, U_KLOC)!Bound_ele_correct unnecessary here?
                        ! endif
                    END DO
                ENDIF
            end if  ! endof if ( integrate_other_side_and_not_boundary ) then
        END DO
        IF ( on_domain_boundary ) THEN
            UDGI_IMP_ALL=0.0
            DO U_KLOC = 1, Mdims%u_nloc
                DO IPHASE = 1,final_phase
                    UDGI_IMP_ALL(:,IPHASE) = UDGI_IMP_ALL(:,IPHASE) + CV_funs%sufen( U_KLOC, GI ) * &
                        UGI_COEF_ELE_ALL( :, IPHASE, U_KLOC ) * LOC_U( :, IPHASE, U_KLOC )
                END DO
            END DO
            DO IPHASE = 1,final_phase
                NDOTQ_IMP(IPHASE)= SUM( CVNORMX_ALL( :,GI ) * UDGI_IMP_ALL(:,IPHASE) )
            END DO
            ct_rhs_phase_cv_nodi=ct_rhs_phase_cv_nodi &
                - SCVDETWEI( GI ) * (  ( &
                + LIMT2  * LIMDT * (NDOTQ-NDOTQ_IMP) &
                ) / DEN_ALL( :, CV_NODI ) )
        END IF
        ! IF ( between_elements ) THEN
        !     ! We have a discontinuity between elements so integrate along the face...
        !     DO U_SKLOC = 1, Mdims%u_snloc
        !         U_KLOC = U_SLOC2LOC(U_SKLOC)
        !         U_KLOC2 = U_OTHER_LOC( U_KLOC )
        !         RCON = SCVDETWEI( GI ) * (  FTHETA_T2 * LIMDT + ONE_M_FTHETA_T2OLD * LIMDTOLD * THETA_VEL) &
        !             * CV_funs%sufen( U_KLOC, GI ) / DEN_ALL( :, CV_NODI  )
        !         DO IPHASE = 1,final_phase
        !             Mmat%CT( :, IPHASE, JCOUNT_KLOC2( U_KLOC2 ) ) &
        !                 = Mmat%CT( :, IPHASE, JCOUNT_KLOC2( U_KLOC2 ) ) &
        !                 + rcon(IPHASE) * UGI_COEF_ELE2_ALL( :, IPHASE, U_KLOC2 ) * CVNORMX_ALL( :, GI )
        !         END DO
        !         IF(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV) THEN
        !             RCON = SCVDETWEI( GI ) * CV_funs%sufen( U_KLOC, GI )
        !             DO IPHASE=1,final_phase!Mdims%nphase
        !                 Mmat%C_CV( :, IPHASE, C_JCOUNT_KLOC2( U_KLOC2 ) ) &
        !                     = Mmat%C_CV( :, IPHASE, C_JCOUNT_KLOC2( U_KLOC2 ) ) &
        !                     + RCON(IPHASE) * CVNORMX_ALL( :, GI )* (1.- Mass_corrector)
        !             END DO
        !         ENDIF
        !         ! flux from the other side (change of sign because normal is -ve)...
        !         if ( integrate_other_side_and_not_boundary ) then
        !             RCON_J = SCVDETWEI( GI ) * ( FTHETA_T2_J* LIMDT + ONE_M_FTHETA_T2OLD_J * LIMDTOLD * THETA_VEL) &
        !                 * CV_funs%sufen( U_KLOC, GI ) / DEN_ALL( :, CV_NODJ )
        !             DO IPHASE=1,final_phase
        !                 Mmat%CT( :, IPHASE, ICOUNT_KLOC2( U_KLOC2 ) ) &
        !                     = Mmat%CT( :, IPHASE, ICOUNT_KLOC2( U_KLOC2 ) ) &
        !                     - RCON_J(IPHASE) * UGI_COEF_ELE2_ALL( :, IPHASE, U_KLOC2 ) * CVNORMX_ALL( :, GI )
        !             END DO
        !             IF(GET_C_IN_CV_ADVDIF_AND_CALC_C_CV) THEN
        !                 RCON_J = SCVDETWEI( GI ) * CV_funs%sufen( U_KLOC, GI )
        !                 DO IPHASE=1,final_phase!Mdims%nphase
        !                     Mmat%C_CV( :, IPHASE, C_ICOUNT_KLOC2( U_KLOC2 ) ) &
        !                         = Mmat%C_CV( :, IPHASE, C_ICOUNT_KLOC2( U_KLOC2 ) ) &
        !                         - RCON_J(IPHASE) * CVNORMX_ALL( :, GI )* (1.-Mass_corrector)!Mass_corrector
        !                 END DO
        !             ENDIF
        !         end if  ! endof if ( integrate_other_side_and_not_boundary ) then
        !     END DO
        ! END IF ! endof IF ( between_elements ) THEN
        RETURN
    contains

        !>This subroutine populates Bound_ele_correct and Bound_ele2_correct to properly apply the BCs when creating the
        !>Mmat%C_CV matrix
        subroutine introduce_C_CV_boundary_conditions(Bound_ele_correct)
            implicit none
            real, dimension(:,:,:), intent(out) :: Bound_ele_correct
            !Local variables
            integer :: U_KLOC, IPHASE, P_SJLOC, U_INOD, ipres, CV_KLOC, P_ILOC, CV_ILOC
            logical, save :: show_warn_msg = .true.
            logical, save :: have_been_read = .false.
            logical, save :: hydrostatic_bc
            real, save :: top_domain
            !By default no modification is required
            Bound_ele_correct = 1.0
            if ( .not. have_been_read ) then
              hydrostatic_bc = have_option( '/material_phase[0]/scalar_field::Pressure/prognostic/hydrostatic_boundaries' )
              top_domain = maxval(X_ALL(Mdims%ndim, :))
              call allmax(top_domain)
              have_been_read = .true.
            end if
            !Get vertical coordinate of top of the domain
            IF ( on_domain_boundary ) THEN
                !By default the position must not be added to the matrix
                Bound_ele_correct = 0.!<= P in the CV == P in the BC, it is done this way
                !If Mmat%C_CV formulation, apply weak pressure boundary conditions if any
                DO IPRES = 1, 1!Mdims%npres! sprint_to_do why not npres???
                    IF( WIC_P_BC_ALL( 1,IPRES,SELE ) == WIC_P_BC_DIRICHLET ) THEN
                        DO U_SILOC = 1, Mdims%u_snloc
                            U_ILOC = U_SLOC2LOC( U_SILOC )
                            U_INOD = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                            DO IPHASE =  1+(IPRES-1)*Mdims%n_in_pres, IPRES*Mdims%n_in_pres!Will this work for Nphases-1
                                !We give priority to velocity boundary conditions
                                if (WIC_U_BC_ALL( 1, IPHASE, SELE ) /= WIC_U_BC_DIRICHLET ) then
                                    !Only in the boundaries with a defined pressure it needs to be added into
                                    !the matrix and into the RHS
                                    if (hydrostatic_bc) then!not compatible with hydrostatic presure solver
                                        Bound_ele_correct( :, IPHASE, U_ILOC ) = 1.
                                        Mmat%U_RHS( :, IPHASE, U_INOD ) = Mmat%U_RHS( :, IPHASE, U_INOD ) &
                                            - CVNORMX_ALL( :, GI )* CV_funs%sufen( U_ILOC, GI )*SCVDETWEI( GI )&
                                            * SUF_P_BC_ALL( 1,1,1 + Mdims%cv_snloc* ( SELE - 1 ) ) - (gravty*&
                                            SUF_D_BC_ALL( 1, 1, 1 + Mdims%cv_snloc* ( SELE - 1 ) )*&
                                            abs(top_domain-X_ALL(Mdims%ndim, CV_NODI))*&
                                            CVNORMX_ALL( :, GI )* CV_funs%sufen( U_ILOC, GI )*SCVDETWEI( GI ))
                                    else
                                        Bound_ele_correct( :, IPHASE, U_ILOC ) = 1.
                                        Mmat%U_RHS( :, IPHASE, U_INOD ) = Mmat%U_RHS( :, IPHASE, U_INOD ) &
                                            - CVNORMX_ALL( :, GI )* CV_funs%sufen( U_ILOC, GI )*SCVDETWEI( GI )&
                                            * SUF_P_BC_ALL( 1,1,1 + Mdims%cv_snloc* ( SELE - 1 ) )
                                    endif

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


    !>Detects whether the element has a shockfront or not
    logical function shock_front_in_ele(ele, Mdims, sat, ndgln, Imble_frac)
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

    !>@brief: In this method we assemble and solve the Laplacian system using at least P1 elements
    !> The equation solved is the following: Div sigma Grad X = - SUM (Div K Grad F) with Neuman BCs = 0
    !> where K and F are passed down as a vector. Therefore for n entries the SUM will be performed over n fields
    !> Example: F = (3, nphase, cv_nonods) would include three terms in the RHS and the same for K
    !> If harmonic average then we perform the harmonic average of sigma and K
    !> IMPORTANT: This subroutine requires the PHsparsity to be generated
    !> Note that this method solves considering FE fields. If using CV you may incur in an small error.
    !>@param Mdims Data type storing all the dimensions describing the mesh, fields, nodes, etc
    !>@param  packed_state Linked list containing all the fields used by IC-FERST, memory partially shared with state
    !>@param  ndgln Global to local variables
    !>@param  Mmat Matrices for ICFERST
    !>@param  Mspars Sparsity of the matrices
    !>@param  CV_funs Shape functions for the CV mesh
    !>@param  CV_GIdims Gauss integration numbers for CV fields
    !>@param Sigma_field Coefficient multipliying the left-hand side term
    !>@param Solution Solution field of type scalar_field, required to ensure consistency
    !>@param K_fields Coefficients multipliying the RHS terms
    !>@param F_fields RHS fierlds
    !>@param intface_val_type  0 = no interpolation; 1 Harmonic mean; !20 for SP solver, harmonic mean considering charge; negative normal mean

    subroutine generate_Laplacian_system( Mdims, packed_state, ndgln, Mmat, Mspars, CV_funs, CV_GIdims, Sigma_field, &
      Solution, K_fields, F_fields, intface_val_type)
      implicit none

      type(multi_dimensions), intent( in ) :: Mdims
      type( state_type ), intent( inout ) :: packed_state
      type(multi_ndgln), intent(in) :: ndgln
      integer, intent(in) :: intface_val_type! 0 = no interpolation; 1 Harmonic mean; !20 for SP solver, harmonic mean considering charge; negative normal mean
      real, dimension(:,:), intent(in) :: Sigma_field
      real, dimension(:,:,:), intent(in) :: K_fields, F_fields
      type( scalar_field ), intent(inout) :: Solution
      type(multi_shape_funs), intent(inout) :: CV_funs
      type(multi_sparsities), intent(in) :: Mspars
      type (multi_matrices), intent(inout) :: Mmat
      type(multi_GI_dimensions), intent(in) :: CV_GIdims
      ! local variables...
      logical :: INTEGRAT_AT_GI, skip, on_domain_boundary
      integer :: stat ,nb, ele, ele2, ele3, sele, cv_siloc, i, local_phases, CV_NODK, cv_kloc,&
      cv_iloc, cv_jloc, cv_nodi, cv_nodj, idim, iphase, GCOUNT, GI, x_nodi, x_nodj, cv_xloc
      real :: HDC, HDLi, HDLj
      type( vector_field ), pointer :: X_ALL, MASS_CV, XC_CV_ALL
      LOGICAL, DIMENSION( Mdims%x_nonods ) :: X_SHARE
      integer, dimension (Mdims%cv_nloc) ::CV_OTHER_LOC
      integer, dimension (Mdims%cv_snloc) :: CV_SLOC2LOC
      integer, dimension (Mdims%u_snloc) :: U_SLOC2LOC
      integer, dimension (Mdims%mat_nloc) ::MAT_OTHER_LOC
      type(csr_sparsity), pointer :: sparsity
      real, dimension(Mdims%ndim) :: GI_coordinate
      integer, dimension(:), pointer :: neighbours
      REAL, DIMENSION( Mdims%ndim, CV_GIdims%scvngi ) :: CVNORMX_ALL
    !   INTEGER, DIMENSION( CV_GIdims%nface, Mdims%totele ) :: FACE_ELE
      !Variables to reduce communications with PETSc when assembling the matrix
      real, dimension(size(F_fields,2)) :: LOC_CV_RHS_I, LOC_CV_RHS_J, LOC_MAT_II, LOC_MAT_JJ, LOC_MAT_IJ, LOC_MAT_JI
      !###Variables for shape function calculation###
      type (multi_dev_shape_funs) :: SdevFuns
      !Local diffusion coefficients
      real, dimension(size(F_fields,2)):: SIGMA_DIFF_COEF_DIVDX
      real, dimension(size(F_fields,1), size(F_fields,2)):: DIFF_COEF_DIVDX
      !Set the number of phases from F_fields
      local_phases = size(F_fields,2)
      !Retrieve node coordinates
      X_ALL => extract_vector_field( packed_state, "PressureCoordinate" )
      !Retrieve CV volume and CV centres
      ! MASS_CV=>extract_vector_field(packed_state,"CVIntegral")
      XC_CV_ALL=>extract_vector_field(packed_state,"CVBarycentre")
      !Allocate matrix and RHS
      sparsity=>extract_csr_sparsity(packed_state,"ACVSparsity")
      call allocate(Mmat%CV_RHS,local_phases,Solution%mesh,"RHS")
      call allocate(Mmat%petsc_ACV,sparsity,[local_phases,local_phases],"ACV_INTENERGE")
      call zero(Mmat%petsc_ACV); Mmat%CV_RHS%val = 0.0

      !Allocate derivatives of the shape functions
      call allocate_multi_dev_shape_funs(CV_funs%scvfenlx_all, CV_funs%sufenlx_all, SdevFuns)

      !Obtain elements surrounding an element (FACE_ELE) only if it is not stored yet
      if (.not. associated(Mmat%FACE_ELE)) then 
        allocate(Mmat%FACE_ELE(CV_GIdims%nface, Mdims%totele))
        Mmat%FACE_ELE = 0.
        CALL CALC_FACE_ELE( Mmat%FACE_ELE, Mdims%totele, Mdims%stotel, CV_GIdims%nface, &
        Mspars%ELE%fin, Mspars%ELE%col, Mdims%cv_nloc, Mdims%cv_snloc, Mdims%cv_nonods, ndgln%cv, ndgln%suf_cv, &
        CV_funs%cv_sloclist, Mdims%x_nloc, ndgln%x )
      end if
      Loop_Elements: do ele = 1, Mdims%totele
        if (IsParallel()) then
          if (.not. assemble_ele(solution,ele)) then
            skip=.true.
            neighbours=>ele_neigh(solution,ele)
            do nb=1,size(neighbours)
              if (neighbours(nb)<=0) cycle
              if (assemble_ele(solution,neighbours(nb))) then
                skip=.false.
                exit
              end if
            end do
            if (skip) cycle
          end if
        end if

        DO CV_KLOC = 1, Mdims%cv_nloc
          CV_NODK = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_KLOC )
          Loop_CV_ILOC: DO CV_ILOC = 1, Mdims%cv_nloc ! Loop over the nodes of the element
            ! Global node number of the local node
            CV_NODI = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
            X_NODI = ndgln%x( ( ELE - 1 ) * Mdims%x_nloc  + CV_ILOC )
            ! MAT_NODI = ndgln%mat( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
            ! Loop over quadrature (gauss) points in ELE neighbouring ILOC
            Loop_GCOUNT: DO GCOUNT = CV_funs%findgpts( CV_ILOC ), CV_funs%findgpts( CV_ILOC + 1 ) - 1
              ! CV_funs%colgpts stores the local Gauss-point number in the ELE
              GI = CV_funs%colgpts( GCOUNT )
              ! Get the neighbouring node for node ILOC and Gauss point GI
              CV_JLOC = CV_funs%cv_neiloc( CV_ILOC, GI )
              ELE2 = 0; SELE = 0; CV_SILOC=0
              INTEGRAT_AT_GI = .TRUE.
              CV_OTHER_LOC=0
              Conditional_CheckingNeighbourhood: IF ( CV_JLOC == -1 ) THEN
                ! We are on the boundary or next to another element.  Determine CV_OTHER_LOC
                ! CV_funs%cvfem_on_face(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
                ! Look for these nodes on the other elements.
                CALL FIND_OTHER_SIDE( CV_OTHER_LOC, Mdims%cv_nloc, Mdims%u_nloc, &
                MAT_OTHER_LOC, INTEGRAT_AT_GI, &
                Mdims%x_nloc, Mdims%xu_nloc, ndgln%x, ndgln%xu, &
                Mdims%cv_snloc, CV_funs%cvfem_on_face( :, GI ), X_SHARE, ELE, ELE2,  &
                Mspars%ELE%fin, Mspars%ELE%col )

                IF ( INTEGRAT_AT_GI ) THEN
                  CV_JLOC = CV_OTHER_LOC( CV_ILOC )
                  SELE = 0
                  ELE3=0
                  IF ( CV_JLOC == 0 ) THEN ! We are on the boundary of the domain or subdomain
                    CV_JLOC = CV_ILOC
                    ! Calculate SELE, CV_SILOC, U_SLOC2LOC, CV_SLOC2LOC
                    CALL CALC_SELE( ELE, ELE3, SELE, CV_SILOC, CV_ILOC, U_SLOC2LOC, CV_SLOC2LOC, &
                    mmat%FACE_ELE, GI, CV_funs, Mdims, CV_GIdims,&
                    ndgln%cv, ndgln%u, ndgln%suf_cv, ndgln%suf_u )
                  END IF
                  INTEGRAT_AT_GI = .NOT.( (ELE==ELE2) .AND. (SELE==0) )
                END IF
              END IF Conditional_CheckingNeighbourhood
              !Avoid integrating across the middle of a CV on the boundaries of elements
              !For disconditnuous pressure this is always true
              Conditional_integration: IF ( INTEGRAT_AT_GI ) THEN
                on_domain_boundary = ( SELE /= 0 )
                CV_NODJ = ndgln%cv( ( ELE - 1 )  * Mdims%cv_nloc + CV_JLOC )
                X_NODJ = ndgln%x( ( ELE - 1 )  * Mdims%cv_nloc + CV_JLOC )

                if(CV_NODJ > CV_NODI) then
                  !Compute SdevFuns%DETWEI and CVNORMX_ALL
                  CALL SCVDETNX( Mdims, ndgln, X_ALL%val, CV_funs, CV_GIdims, on_domain_boundary, .false., &!NOT FULLY DG FOR THIS METHOD
                  ELE, GI, SdevFuns%DETWEI, CVNORMX_ALL, XC_CV_ALL%val( :, CV_NODI ), X_NODI, X_NODJ)
                  ! Obtain the CV discretised advection/diffusion equations
                  IF(.not. on_domain_boundary) THEN
                    GI_coordinate = 0.
                    !Obtain the coordinate at the edge between both CVs using shape functions
                    do cv_xloc = 1, Mdims%x_nloc
                      GI_coordinate = GI_coordinate + CV_funs%scvfen( cv_xloc , GI ) * X_ALL%val(:,ndgln%x((ELE-1) * Mdims%x_nloc  + cv_xloc ))
                    end do
                    !Distance from i node to edge!
                    HDLi = SQRT( SUM( (XC_CV_ALL%val(:,CV_NODI)-GI_coordinate)**2) )
                    ! Compute the distance HDC between the nodes either side of the CV face
                    HDC = SQRT( SUM( (XC_CV_ALL%val(:,CV_NODI)-XC_CV_ALL%val(:,CV_NODJ))**2) )
                    !Distance from j node to edge
                    HDLj = SQRT( SUM( (XC_CV_ALL%val(:,CV_NODJ)-GI_coordinate)**2) )
                    do iphase = 1, local_phases
                      do i = 1, size(K_fields,1)
                        DIFF_COEF_DIVDX(i, iphase) = get_DIFF_COEF_DIVDX(2*intface_val_type, HDC, HDLi, HDLj, K_fields(i, iphase, cv_nodi), K_fields(i, iphase, cv_nodj)&
                        ,Sigma_field(iphase, cv_nodi), Sigma_field(iphase, cv_nodj))
                      end do
                      SIGMA_DIFF_COEF_DIVDX(iphase) = get_DIFF_COEF_DIVDX(intface_val_type, HDC, HDLi, HDLj, Sigma_field(iphase, cv_nodi), Sigma_field(iphase, cv_nodj))
                    end do
                  else
                    DIFF_COEF_DIVDX = 0
                    SIGMA_DIFF_COEF_DIVDX = 0.
                  ENDIF
                  LOC_CV_RHS_I=0.0; LOC_MAT_II =0.
                  LOC_CV_RHS_J=0.0; LOC_MAT_JJ =0.
                  LOC_MAT_IJ = 0.0; LOC_MAT_JI =0.
                  !Assemble
                  do iphase = 1, local_phases
                    LOC_MAT_IJ(iphase) = LOC_MAT_IJ(iphase) - SdevFuns%DETWEI( GI ) * SIGMA_DIFF_COEF_DIVDX(iphase)
                    !Assemble off-diagonal
                    LOC_MAT_JI(iphase) = LOC_MAT_JI(iphase) - SdevFuns%DETWEI( GI ) * SIGMA_DIFF_COEF_DIVDX(iphase)
                    !Assemble diagonal of the matrix of node cv_nodi
                    LOC_MAT_II(iphase) = LOC_MAT_II(iphase) + SdevFuns%DETWEI( GI ) * SIGMA_DIFF_COEF_DIVDX(iphase)
                    !Assemble diagonal of the matrix of node cv_nodj
                    LOC_MAT_JJ(iphase) = LOC_MAT_JJ(iphase) + SdevFuns%DETWEI( GI ) * SIGMA_DIFF_COEF_DIVDX(iphase)
                    ! Fill up RHS
                    do i = 1, size(K_fields,1)
                      LOC_CV_RHS_I(iphase) =  LOC_CV_RHS_I(iphase) - SdevFuns%DETWEI(GI) * DIFF_COEF_DIVDX(i, iphase) * &
                      (F_fields(i, iphase, cv_nodj) - F_fields(i, iphase, cv_nodi))
                      LOC_CV_RHS_J(iphase) =  LOC_CV_RHS_J(iphase) - SdevFuns%DETWEI(GI) * DIFF_COEF_DIVDX(i, iphase) * &
                      (F_fields(i, iphase, cv_nodi) - F_fields(i, iphase, cv_nodj))
                    end do
                  end do

                  do iphase = 1, local_phases
                    !For the RHS collapsing to assemble into phase 1 can be done just here
                    call addto(Mmat%CV_RHS,iphase, CV_NODI,LOC_CV_RHS_I(iphase))
                    call addto(Mmat%CV_RHS,iphase, CV_NODJ,LOC_CV_RHS_J(iphase))
                    !Introduce the information into the petsc_ACV matrix
                    call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodi,cv_nodi, LOC_MAT_II(iphase) )
                    call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodj,cv_nodj, LOC_MAT_JJ(iphase) )
                    call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodi,cv_nodj, LOC_MAT_IJ(iphase) )
                    call addto(Mmat%petsc_ACV,iphase,iphase,cv_nodj,cv_nodi, LOC_MAT_JI(iphase) )
                  end do
                END IF
              END IF Conditional_integration
            END DO Loop_GCOUNT
          END DO Loop_CV_ILOC
        end do
      END DO Loop_Elements
      call deallocate_multi_dev_shape_funs(SdevFuns)

    contains
      !>@brief: Computes the effective value of K or sigma. at the interface between CVs
      !> Based on the intface_type integer:
      !> 0 = no interpolation; 1 Harmonic mean; !20 for SelfPotential solver, harmonic mean considering charge; negative normal mean
      real function get_DIFF_COEF_DIVDX(intface_type, HDC,  W_i, W_j, Value_i, Value_j, sigma_i, sigma_j)
        implicit none
        integer, intent(in) :: intface_type
        real, intent(in) :: Value_i, Value_j, W_i, W_j, HDC
        real, optional, intent(in) :: sigma_i, sigma_j
        !Local variable
        logical :: div_by_zero

        div_by_zero =  abs(Value_i + Value_j) < 1e-15
        if (intface_type == 0 ) then!No mean
          get_DIFF_COEF_DIVDX = Value_i
          !Harmonic mean, also used for Rock saturated conductivity
        else if ((intface_type <= 20 .and. intface_type > 0) .and. .not. div_by_zero) then
          get_DIFF_COEF_DIVDX = Value_i * Value_j * (W_j + W_i)/(Value_i *W_j + Value_j*W_i )
        !   !40 is for the mean of the coupling terms
        else if (intface_type > 20 .and.  (abs(Sigma_i + Sigma_j) > 1e-15) ) then
          get_DIFF_COEF_DIVDX = (Value_i * Sigma_j * W_i + Sigma_i * Value_j *W_j) / &
          (W_j*Sigma_i + W_i*Sigma_j)
        else !Normal mean
          get_DIFF_COEF_DIVDX = (Value_i * W_i + Value_j *W_j) / (W_j + W_i)
        end if
        !Divide now by the distance between nodes so it is a diffusion coefficient
        get_DIFF_COEF_DIVDX = get_DIFF_COEF_DIVDX/HDC
      end function get_DIFF_COEF_DIVDX

    end subroutine generate_Laplacian_system
    !>     ---------------------------------------------------------------
    !>     - this subroutine computed the surface area at the Gi point (SCVDETWEI)
    !>     - this subroutine calculates the control volume (CV)
    !>     - CVNORMX, CVNORMY, CVNORMZ normals at the Gaussian
    !>     - integration points GI. NODI = the current global
    !>     - node number for the co-ordinates.
    !>     - (XC,YC,ZC) is the centre of CV NODI
    !>
    !>     ---------------------------------------------------------------
    !>    - date last modified : 25/08/2020
    !>     ---------------------------------------------------------------
    SUBROUTINE SCVDETNX( Mdims, ndgln, X_ALL, CV_funs, CV_GIdims, on_domain_boundary, between_elements, ELE, GI,SCVDETWEI, CVNORMX_ALL,XC_ALL, X_NOD, X_NODJ)
      IMPLICIT NONE
      type(multi_dimensions), intent( in ) :: Mdims
      type(multi_ndgln), intent(in) :: ndgln
      real, dimension(:,:) :: X_ALL
      type(multi_shape_funs), intent(in) :: CV_funs
      type(multi_GI_dimensions), intent(in) :: CV_GIdims
      INTEGER, intent( in ) :: ELE, GI, X_NOD, X_NODJ
      REAL, DIMENSION( Mdims%ndim ), intent( in ) ::   XC_ALL
      REAL, DIMENSION( Mdims%ndim, CV_GIdims%scvngi ), intent( inout ) :: CVNORMX_ALL
      REAL, DIMENSION( : ), intent( inout ) :: SCVDETWEI
      logical, intent(in) :: on_domain_boundary, between_elements
      !     - Local variables
      INTEGER :: NODJ,  JLOC
      REAL :: A, B, C
      REAL :: DETJ
      REAL :: DXDLX, DXDLY, DYDLX
      REAL :: DYDLY, DZDLX, DZDLY
      REAL :: TWOPI
      REAL :: RGI, RDUM
      real, dimension(3) :: POSVGI, DLX,DLY 
      !ewrite(3,*)' In SCVDETNX'
      POSVGI = 0.0
      Conditional_Dimension: IF( Mdims%ndim == 3 ) THEN

        !ewrite(3,*)' In SCVDETNX'
  
        POSVGI = 0.0
        DLX = 0.; DLY=0.
        do JLOC = 1, Mdims%x_nloc
          NODJ = ndgln%x((ELE-1)*Mdims%x_nloc+JLOC)
          DLX = DLX + CV_funs%scvfenslx(JLOC,GI)*X_ALL(:,NODJ)
          DLY = DLY + CV_funs%scvfensly(JLOC,GI)*X_ALL(:,NODJ)
          POSVGI = POSVGI + CV_funs%scvfen(JLOC,GI)*X_ALL(:,NODJ)
        end do
  
        !To calculate the sign of the normal an average between the center of the continuous CV and the center of mass is used
        !this is required as the center of mass has shown not to be reliable and the center of the continuous CV is a particular point that can lead
        !to failures to obtain the sign (perpendicular vectors in a flat boundary); For discontinuous and boundaries we use the old method
        IF ( on_domain_boundary .or. between_elements) then!sprint_to_do between elements use both barycentres?
          POSVGI = POSVGI - (0.8*X_ALL(1:Mdims%ndim, X_NOD) + 0.2*XC_ALL(1:Mdims%ndim))
        else !Use centres of the continuous control volumes, i.e. corners of the elements
          POSVGI = X_ALL(1:Mdims%ndim, X_NODJ) - X_ALL(1:Mdims%ndim, X_NOD)
        end if
        !Obtain normal using the rotational of the two vectors
        CVNORMX_ALL(:,GI) = cross_product(DLX,DLY)
        !     - Calculate the determinant of the Jacobian at Gauss pnt GI.
        DETJ = sum(DLX*DLX)
        !Sign and normalisation value
        A = SIGN( 1.0 / DETJ, sum(CVNORMX_ALL(:,GI) * POSVGI) )
        !Normalise
        CVNORMX_ALL(:,GI) = CVNORMX_ALL(:,GI) * A
        !Calculate the determinant times the surface weight at Gauss pnt GI.
        SCVDETWEI(GI) = DETJ*CV_funs%scvfeweigh(GI)

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

        DETJ = SQRT( DXDLX*DXDLX + DYDLX*DYDLX )
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

      ENDIF Conditional_Dimension

      END SUBROUTINE SCVDETNX


end module cv_advection
