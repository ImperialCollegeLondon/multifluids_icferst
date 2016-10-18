
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

module multiphase_1D_engine

    use elements
    use field_options
    use state_module
    use spud
    use global_parameters
    use futils, only: int2str

    use Fields_Allocates, only : allocate
    use sparse_tools_petsc

    use solvers_module
    use cv_advection
    use matrix_operations
    use shape_functions
    use shape_functions_NDim
    use shape_functions_prototype
    use matrix_operations
    use spact
    use Copy_Outof_State
    use multiphase_EOS
    use Copy_Outof_State, only: as_vector
    use fldebug
    use solvers
    use memory_diagnostics
    use reference_counting
    use multi_data_types
    use Compositional_Terms
    use multi_pipes
    use multi_surface_tension
    use multi_tools, only: CALC_FACE_ELE
    implicit none

    private :: CV_ASSEMB_FORCE_CTY, ASSEMB_FORCE_CTY

    public  :: INTENERGE_ASSEM_SOLVE, VolumeFraction_Assemble_Solve, &
    FORCE_BAL_CTY_ASSEM_SOLVE

contains

  SUBROUTINE INTENERGE_ASSEM_SOLVE( state, packed_state, &
       Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd,&
       tracer, velocity, density,  DT, &
       SUF_SIG_DIAGTEN_BC,  VOLFRA_PORE, &
       IGOT_T2, igot_theta_flux,GET_THETA_FLUX, USE_THETA_FLUX,  &
       THETA_GDIFF, IDs_ndgln, &
       option_path, &
       mass_ele_transp, &
       thermal, THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, &
       icomp, saturation, Permeability_tensor_field )
           ! Solve for internal energy using a control volume method.
           implicit none
           type( state_type ), dimension( : ), intent( inout ) :: state
           type( state_type ), intent( inout ) :: packed_state
           type(multi_dimensions), intent(in) :: Mdims
           type(multi_GI_dimensions), intent(in) :: CV_GIdims
           type(multi_shape_funs), intent(inout) :: CV_funs
           type (multi_sparsities), intent(in) :: Mspars
           type(multi_ndgln), intent(in) :: ndgln
           type (multi_discretization_opts) :: Mdisopt
           type (multi_matrices), intent(inout) :: Mmat
           type (porous_adv_coefs), intent(inout) :: upwnd
           type(tensor_field), intent(inout) :: tracer
           type(tensor_field), intent(in) :: velocity, density
           INTEGER, intent( in ) :: IGOT_T2, igot_theta_flux
           LOGICAL, intent( in ) :: GET_THETA_FLUX, USE_THETA_FLUX
           LOGICAL, intent( in ), optional ::THERMAL
           INTEGER, DIMENSION( : ), intent( in ) :: IDs_ndgln
           REAL, DIMENSION( :, : ), intent( inout ) :: THETA_GDIFF
           REAL, DIMENSION( :,: ), intent( inout ), optional :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
           REAL, intent( in ) :: DT
           REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
           REAL, DIMENSION( :, : ), intent( in ) :: VOLFRA_PORE
           character( len = * ), intent( in ), optional :: option_path
           real, dimension( : ), intent( inout ), optional :: mass_ele_transp
           type(tensor_field), intent(in), optional :: saturation
           type( tensor_field ), optional, pointer, intent(in) :: Permeability_tensor_field
           integer, optional :: icomp
           ! Local variables
           LOGICAL, PARAMETER :: GETCV_DISC = .TRUE., GETCT= .FALSE.
           integer :: nits_flux_lim, its_flux_lim
           logical :: lump_eqns
           REAL, DIMENSION( :, : ), allocatable :: DIAG_SCALE_PRES
           REAL, DIMENSION( :, :, : ), allocatable :: DIAG_SCALE_PRES_COUP, GAMMA_PRES_ABS, GAMMA_PRES_ABS_NANO, INV_B
           REAL, DIMENSION( :,:,:, : ), allocatable :: TDIFFUSION
           REAL, DIMENSION( : ), ALLOCATABLE :: MASS_PIPE, MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE
           real, dimension( size(Mspars%small_acv%col )) ::  mass_mn_pres
           REAL, DIMENSION( : , : ), allocatable :: den_all, denold_all, t_source
           REAL, DIMENSION( : ), allocatable :: CV_RHS_SUB
           type( tensor_field ), pointer :: P, Q
           INTEGER :: IPHASE
           REAL, PARAMETER :: SECOND_THETA = 1.0
           LOGICAL :: RETRIEVE_SOLID_CTY
           type( tensor_field ), pointer :: den_all2, denold_all2, a, aold, deriv, Component_Absorption
           type( vector_field ), pointer  :: MeanPoreCV
           integer :: lcomp, Field_selector, IGOT_T2_loc
           type(vector_field)  :: vtracer
           type(csr_sparsity), pointer :: sparsity
           real, dimension(:,:,:), allocatable :: Velocity_Absorption
           real, dimension(:,:,:), pointer :: T_AbsorB=>null()
           integer :: IGOT_THERM_VIS
           real, dimension(:,:), allocatable :: THERM_U_DIFFUSION_VOL
           real, dimension(:,:,:,:), allocatable :: THERM_U_DIFFUSION
           integer :: ncomp_diff_coef, comp_diffusion_opt
           real, dimension(:,:,:), allocatable :: Component_Diffusion_Operator_Coefficient
           type( tensor_field ), pointer :: perm
           integer :: cv_disopt, cv_dg_vel_int_opt
           real :: cv_theta, cv_beta
           type( scalar_field ), pointer :: sfield

            if (present(Permeability_tensor_field)) then
                perm => Permeability_tensor_field
            else
                perm=>extract_tensor_field(packed_state,"Permeability")
            end if
           IGOT_THERM_VIS = 0
           ALLOCATE( THERM_U_DIFFUSION(Mdims%ndim,Mdims%ndim,Mdims%nphase,Mdims%mat_nonods*IGOT_THERM_VIS ) )
           ALLOCATE( THERM_U_DIFFUSION_VOL(Mdims%nphase,Mdims%mat_nonods*IGOT_THERM_VIS ) )

           lcomp = 0
           if ( present( icomp ) ) lcomp = icomp

           call allocate(Mmat%CV_RHS,Mdims%nphase,tracer%mesh,"RHS")
           sparsity=>extract_csr_sparsity(packed_state,"ACVSparsity")
           call allocate(Mmat%petsc_ACV,sparsity,[Mdims%nphase,Mdims%nphase],"ACV",.false.,.false.)
           call zero( Mmat%petsc_ACV )

           allocate(den_all(Mdims%nphase,Mdims%cv_nonods),denold_all(Mdims%nphase,Mdims%cv_nonods))
           allocate( DIAG_SCALE_PRES(0,0) )
           allocate( DIAG_SCALE_PRES_COUP(0,0,0),GAMMA_PRES_ABS(0,0,0),GAMMA_PRES_ABS_NANO(0,0,0),INV_B(0,0,0) )
           allocate( MASS_PIPE(0), MASS_CVFEM2PIPE(0), MASS_PIPE2CVFEM(0), MASS_CVFEM2PIPE_TRUE(0) )
           allocate( T_SOURCE( Mdims%nphase, Mdims%cv_nonods ) ) ; T_SOURCE=0.0

           IGOT_T2_loc = 0

!!-PY changed it for k_epsilon model
           if ( thermal .or. trim( option_path ) == '/material_phase[0]/scalar_field::Temperature' &
                .or. trim( option_path ) == '/material_phase[0]/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentKineticEnergy' &
                .or. trim( option_path ) == '/material_phase[0]/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentDissipation' ) then
!            if ( thermal .or. trim( option_path ) == '/material_phase[0]/scalar_field::Temperature') then


               p => extract_tensor_field( packed_state, "PackedCVPressure" )
               den_all2 => extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )
               denold_all2 => extract_tensor_field( packed_state, "PackedOldDensityHeatCapacity" )
               den_all    = den_all2 % val ( 1, :, : )
               denold_all = denold_all2 % val ( 1, :, : )
               ! open the boiling test for two phases-gas and liquid
               if (have_option('/boiling')) then ! don't the divide int. energy equation by the volume fraction
                   a => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
                   den_all = den_all * a%val(1,:,:)
                   aold => extract_tensor_field( packed_state, "PackedOldPhaseVolumeFraction" )
                   denold_all = denold_all * aold%val(1,:,:)
               end if
               IGOT_T2_loc = 1
           else if ( lcomp > 0 ) then
               p => extract_tensor_field( packed_state, "PackedFEPressure" )
               den_all2 => extract_tensor_field( packed_state, "PackedComponentDensity" )
               denold_all2 => extract_tensor_field( packed_state, "PackedOldComponentDensity" )
               den_all = den_all2 % val ( 1, :, : )
               denold_all = denold_all2 % val ( 1,  :, : )
           else
               p => extract_tensor_field( packed_state, "PackedFEPressure" )
               den_all=1.0
               denold_all=1.0
           end if
           if( present( option_path ) ) then ! solving for Temperature or Internal Energy or k_epsilon model

!!-PY this part need work for k_epsilon model
               if( trim( option_path ) == '/material_phase[0]/scalar_field::Temperature' ) then
                   call get_option( '/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation/' // &
                       'control_volumes/number_advection_iterations', nits_flux_lim, default = 3 )
                   Field_selector = 1
                   Q => extract_tensor_field( packed_state, "PackedTemperatureSource" )
                   T_source( :, : ) = Q % val( 1, :, : )

               end if
             


               if( trim( option_path ) == '/material_phase[0]/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentKineticEnergy' ) then
                   call get_option( '/material_phase[0]/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentKineticEnergy/prognostic/temporal_discretisation/' // &
                       'control_volumes/number_advection_iterations', nits_flux_lim, default = 3 )


                   Field_selector = 1
                   Q => extract_tensor_field( packed_state, "PackedTurbulentKineticEnergySource" )
                   T_source( :, : ) = Q % val( 1, :, : )
                   
              




               else if( trim( option_path ) == '/material_phase[0]/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentDissipation' ) then
                   call get_option( '/material_phase[0]/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentDissipation/prognostic/temporal_discretisation/' // &
                       'control_volumes/number_advection_iterations', nits_flux_lim, default = 3 )


                   Field_selector = 1
                   Q => extract_tensor_field( packed_state, "PackedTurbulentDissipationSource" )
                   T_source( :, : ) = Q % val( 1, :, : )
                   
               end if


               cv_disopt = Mdisopt%t_disopt
               cv_dg_vel_int_opt = Mdisopt%t_dg_vel_int_opt
               cv_theta = Mdisopt%t_theta
               cv_beta = Mdisopt%t_beta

           else ! solving for Composition
               call get_option( '/material_phase[' // int2str( Mdims%nphase ) // ']/scalar_field::ComponentMassFractionPhase1/' // &
                   'prognostic/temporal_discretisation/control_volumes/number_advection_iterations', nits_flux_lim, default = 1 )
               Field_selector = 2
               IGOT_T2_loc = IGOT_T2
               cv_disopt = Mdisopt%v_disopt
               cv_dg_vel_int_opt = Mdisopt%v_dg_vel_int_opt
               cv_theta = Mdisopt%v_theta
               cv_beta = Mdisopt%v_beta
           end if

           lump_eqns = have_option( '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
               'spatial_discretisation/continuous_galerkin/mass_terms/lump_mass_matrix' )

           RETRIEVE_SOLID_CTY = .false.
           if ( have_option( '/blasting' ) ) RETRIEVE_SOLID_CTY = .true.

           deriv => extract_tensor_field( packed_state, "PackedDRhoDPressure" )
           allocate( TDIFFUSION( Mdims%mat_nonods, Mdims%ndim, Mdims%ndim, Mdims%nphase ) ) ; TDIFFUSION=0.0





!!-PY changed it for k_epsilon model
           if ( thermal .or. trim( option_path ) == '/material_phase[0]/scalar_field::Temperature' &
                .or. trim( option_path ) == '/material_phase[0]/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentKineticEnergy' &
                .or. trim( option_path ) == '/material_phase[0]/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentDissipation') then
 !          if ( thermal .or. trim( option_path ) == '/material_phase[0]/scalar_field::Temperature') then
 
              call calculate_diffusivity( state, Mdims, ndgln, TDIFFUSION, tracer )
           end if

           ! get diffusivity for compositional
           if ( lcomp > 0 .and. is_porous_media ) then
              ncomp_diff_coef = 0 ; comp_diffusion_opt = 0
              allocate( Component_Diffusion_Operator_Coefficient( Mdims%ncomp, ncomp_diff_coef, Mdims%nphase ) )
              Component_Diffusion_Operator_Coefficient = 0.0
              call Calculate_ComponentDiffusionTerm( packed_state, &
                 Mdims, CV_GIdims, CV_funs, &
                 ndgln%mat, ndgln%u, ndgln%x, &
                 ncomp_diff_coef, comp_diffusion_opt, &
                 Component_Diffusion_Operator_Coefficient( icomp, :, : ), &
                 TDiffusion )
              deallocate( Component_Diffusion_Operator_Coefficient )
              Component_Absorption => extract_tensor_field( packed_state, "ComponentAbsorption")
              T_ABSORB => Component_Absorption%val
           end if

           ! calculate T_ABSORB
           if (have_option('/boiling')) then
              allocate ( T_AbsorB( Mdims%nphase, Mdims%nphase, Mdims%cv_nonods ) ) ; T_AbsorB=0.0
              allocate ( Velocity_Absorption( Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase, Mdims%mat_nonods ) )
              call boiling( state, packed_state, Mdims%cv_nonods, Mdims%mat_nonods, Mdims%nphase, Mdims%ndim, &
                   velocity_absorption, T_AbsorB )
              deallocate ( Velocity_Absorption )
           end if

           if ( have_option( "/magma" ) ) then

              ! set the absorption for magma sims here
              sfield => extract_scalar_field( state(1), "TemperatureAbsorption")
              T_ABSORB(1:1,1:1,1:Mdims%cv_nonods) => sfield%val ! only phase 1
           end if


           MeanPoreCV=>extract_vector_field(packed_state,"MeanPoreCV")

           Loop_NonLinearFlux: DO ITS_FLUX_LIM = 1, NITS_FLUX_LIM
               !before the sprint in this call the small_acv sparsity was passed as cmc sparsity...
               call CV_ASSEMB( state, packed_state, &
                   Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd, &
                   tracer, velocity, density, &
                   DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, GAMMA_PRES_ABS, GAMMA_PRES_ABS_NANO, &
                   INV_B, MASS_PIPE, MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE, &
                   DEN_ALL, DENOLD_ALL, &
                   TDIFFUSION, IGOT_THERM_VIS, THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL,&
                   cv_disopt, cv_dg_vel_int_opt, DT, cv_theta, SECOND_THETA, cv_beta, &
                   SUF_SIG_DIAGTEN_BC, &
                   DERIV%val(1,:,:), P%val, &
                   T_SOURCE, T_ABSORB, VOLFRA_PORE, &
                   GETCV_DISC, GETCT, &
                   IGOT_T2_loc,IGOT_THETA_FLUX ,GET_THETA_FLUX, USE_THETA_FLUX, &
                   THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
                   MeanPoreCV%val, &
                   mass_Mn_pres, THERMAL, RETRIEVE_SOLID_CTY, &
                   .false.,  mass_Mn_pres, &
                   mass_ele_transp, IDs_ndgln, &
                   saturation=saturation, Permeability_tensor_field = perm)
               Conditional_Lumping: IF ( LUMP_EQNS ) THEN
                   ! Lump the multi-phase flow eqns together
                   ALLOCATE( CV_RHS_SUB( Mdims%cv_nonods ) )
                   CV_RHS_SUB = 0.0
                   DO IPHASE = 1, Mdims%nphase
                       CV_RHS_SUB( : ) = CV_RHS_SUB( : )&
                           + Mmat%CV_RHS%val(iphase,:)
                   END DO
               ELSE
                   IF ( IGOT_T2 == 1) THEN
                       vtracer=as_vector(tracer,dim=2)
                       call zero_non_owned(Mmat%CV_RHS)
                       call zero(vtracer)
                       call petsc_solve(vtracer,Mmat%petsc_ACV,Mmat%CV_RHS,'/material_phase::Component1/scalar_field::ComponentMassFractionPhase1/prognostic')
                   ELSE
                       vtracer=as_vector(tracer,dim=2)
                       call zero_non_owned(Mmat%CV_RHS)
                       call zero(vtracer)
                       call petsc_solve(vtracer,Mmat%petsc_ACV,Mmat%CV_RHS,trim(option_path))
                       do iphase = 1, Mdims%nphase
                           ewrite(2,*) 'T phase min_max:', iphase, &
                               minval(tracer%val(1,iphase,:)), maxval(tracer%val(1,iphase,:))
                       end do
                   END IF
               END IF Conditional_Lumping
           END DO Loop_NonLinearFlux

           if (have_option('/boiling')) deallocate(T_absorb)

           call deallocate(Mmat%petsc_ACV)
           call deallocate(Mmat%CV_RHS)
           ewrite(3,*) 'Leaving INTENERGE_ASSEM_SOLVE'
  END SUBROUTINE INTENERGE_ASSEM_SOLVE




    subroutine VolumeFraction_Assemble_Solve( state,packed_state, &
         Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd, &
         DT, SUF_SIG_DIAGTEN_BC, &
         V_SOURCE, VOLFRA_PORE, &
         igot_theta_flux, mass_ele_transp,&
         nonlinear_iteration, IDs_ndgln,&
         IDs2CV_ndgln, Courant_number,&
         option_path,&
         THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, Quality_list)
             implicit none
             type( state_type ), dimension( : ), intent( inout ) :: state
             type( state_type ) :: packed_state
             type(multi_dimensions), intent(in) :: Mdims
             type(multi_GI_dimensions), intent(in) :: CV_GIdims
             type(multi_shape_funs), intent(inout) :: CV_funs
             type (multi_sparsities), intent(in) :: Mspars
             type(multi_ndgln), intent(in) :: ndgln
             type (multi_discretization_opts) :: Mdisopt
             type (multi_matrices), intent(inout) :: Mmat
             type (porous_adv_coefs), intent(inout) :: upwnd
             INTEGER, intent( in ) :: igot_theta_flux
             INTEGER, DIMENSION( : ), intent( in ) :: IDs_ndgln
             integer, dimension(:), intent(in)  :: IDs2CV_ndgln
             REAL, intent( in ) :: DT
             REAL, DIMENSION( :, : ), intent( inout ) :: SUF_SIG_DIAGTEN_BC
             REAL, DIMENSION( :, : ), intent( in ) :: V_SOURCE
             !REAL, DIMENSION( :, :, : ), intent( in ) :: V_ABSORB
             REAL, DIMENSION( :, : ), intent( in ) :: VOLFRA_PORE
             real, dimension( : ), intent( inout ) :: mass_ele_transp
             integer, intent(in) :: nonlinear_iteration
             real, intent(inout) :: Courant_number
             character(len= * ), intent(in), optional :: option_path
             REAL, DIMENSION( :, :), intent( inout ), optional :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
             ! Local Variables
             LOGICAL, PARAMETER :: THERMAL= .false.
             integer :: igot_t2
             REAL, DIMENSION( : ), allocatable :: mass_mn_pres
             REAL, DIMENSION( :, : ), allocatable :: DIAG_SCALE_PRES
             REAL, DIMENSION( :, :, : ), allocatable :: DIAG_SCALE_PRES_COUP, GAMMA_PRES_ABS, GAMMA_PRES_ABS_NANO, INV_B
             REAL, DIMENSION( : ), ALLOCATABLE :: MASS_PIPE, MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE
             REAL, DIMENSION( :,:,:,: ), allocatable :: TDIFFUSION
             REAL, DIMENSION( :, : ), allocatable :: THETA_GDIFF
             REAL, DIMENSION( :, : ), pointer :: DEN_ALL, DENOLD_ALL
             REAL, DIMENSION( :, : ), allocatable :: T2, T2OLD
             REAL, DIMENSION( :, : ), allocatable :: MEAN_PORE_CV
             REAL, DIMENSION( :, :, :, : ), allocatable :: THERM_U_DIFFUSION
             REAL, DIMENSION( :, : ), allocatable :: THERM_U_DIFFUSION_VOL
             LOGICAL :: GET_THETA_FLUX
             REAL , PARAMETER :: SECOND_THETA = 1.0
             INTEGER :: STAT, IGOT_THERM_VIS, IPHASE, JPHASE, IPHASE_REAL, JPHASE_REAL, IPRES, JPRES
             LOGICAL, PARAMETER :: GETCV_DISC = .TRUE., GETCT= .FALSE., RETRIEVE_SOLID_CTY= .FALSE.
             type( tensor_field ), pointer :: den_all2, denold_all2
             ! Element quality fix
             type(bad_elements), dimension(:), optional :: Quality_list

             !Working pointers
             real, dimension(:,:,:), pointer :: p, V_ABSORB => null() ! this is PhaseVolumeFraction_AbsorptionTerm
             real, dimension(:, :), pointer :: satura
             type(tensor_field), pointer :: tracer, velocity, density, deriv, PorousMedia_AbsorptionTerm
             type(scalar_field), pointer :: gamma
             !Variable to assign an automatic maximum backtracking parameter based on the Courant number
             logical :: Auto_max_backtrack
             !Variables for global convergence method
             real :: backtrack_par_factor
             type(vector_field)  :: vtracer, residual
             type(csr_sparsity), pointer :: sparsity
             !Variables for capillary pressure
             real, dimension(Mdims%cv_nonods) :: OvRelax_param
             integer :: Phase_with_Pc
             !Variables to stabilize the non-linear iteration solver
             real, dimension(Mdims%nphase, Mdims%cv_nonods) :: sat_bak, backtrack_sat
             real :: Previous_convergence, updating, new_backtrack_par, aux, resold, first_res
             real, save :: res = -1
             logical :: satisfactory_convergence
             integer :: its, useful_sats


             !Extract variables from packed_state
             !call get_var_from_packed_state(packed_state,FEPressure = P)
             call get_var_from_packed_state(packed_state,CVPressure = P)
             call get_var_from_packed_state(packed_state,PhaseVolumeFraction = satura)
             !Get information for capillary pressure to be use in CV_ASSEMB
             call getOverrelaxation_parameter(packed_state, Mdims, ndgln, OvRelax_param, Phase_with_Pc, IDs2CV_ndgln)
             !Get variable for global convergence method
             if (.not. have_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration')) then
                 backtrack_par_factor = 1.1
             else !Get value with the default value of 1.
                 call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Backtracking_factor',&
                     backtrack_par_factor, default = 1.0)
             end if
             if ( backtrack_par_factor < 1.01 ) PorousMedia_AbsorptionTerm => extract_tensor_field( packed_state, "PorousMedia_AbsorptionTerm" )
             !For backtrack_par_factor == -10 we will set backtrack_par_factor based on the shock front Courant number
             Auto_max_backtrack = (backtrack_par_factor == -10)

             GET_THETA_FLUX = .FALSE.
             IGOT_T2 = 0
             deriv => extract_tensor_field( packed_state, "PackedDRhoDPressure" )
             !ALLOCATE( T2( Mdims%cv_nonods * Mdims%nphase * IGOT_T2 ))
             !ALLOCATE( T2OLD( Mdims%cv_nonods * Mdims%nphase * IGOT_T2 ))
             IF ( IGOT_T2 == 1 ) THEN
                 ALLOCATE( T2( Mdims%nphase, Mdims%cv_nonods ))
                 ALLOCATE( T2OLD( Mdims%nphase, Mdims%cv_nonods ))
             END IF
             ALLOCATE( THETA_GDIFF( Mdims%nphase * IGOT_T2, Mdims%cv_nonods * IGOT_T2 ))
             ewrite(3,*) 'In VOLFRA_ASSEM_SOLVE'
             ALLOCATE( mass_mn_pres(size(Mspars%small_acv%col)) ) ; mass_mn_pres = 0.
             ALLOCATE( Mmat%CT( 0,0,0 ) )
             ALLOCATE( DIAG_SCALE_PRES( 0,0 ) )
             ALLOCATE( DIAG_SCALE_PRES_COUP( 0,0,0 ), GAMMA_PRES_ABS( Mdims%nphase,Mdims%nphase,Mdims%cv_nonods ), GAMMA_PRES_ABS_NANO( Mdims%nphase,Mdims%nphase,Mdims%cv_nonods ), INV_B( 0,0,0 ) )
             ALLOCATE( TDIFFUSION( Mdims%mat_nonods, Mdims%ndim, Mdims%ndim, Mdims%nphase ) ) ; TDIFFUSION = 0.
             ALLOCATE( MEAN_PORE_CV( Mdims%npres, Mdims%cv_nonods ) )
             !Shouldn't these allocate bellow be only for NPRES > 1?
             allocate(MASS_PIPE(Mdims%cv_nonods), MASS_CVFEM2PIPE(Mspars%CMC%ncol), MASS_PIPE2CVFEM(Mspars%CMC%ncol), MASS_CVFEM2PIPE_TRUE(Mspars%CMC%ncol))
             gamma=>extract_scalar_field(state(1),"Gamma1",stat)
             GAMMA_PRES_ABS = 0.0
             do ipres = 1, Mdims%npres
                 do iphase = 1+(ipres-1)*Mdims%n_in_pres, ipres*Mdims%n_in_pres
                     do jpres = 1, Mdims%npres
                         if ( ipres /= jpres ) then
                             do jphase = 1+(jpres-1)*Mdims%n_in_pres, jpres*Mdims%n_in_pres
                                 iphase_real = iphase-(ipres-1)*Mdims%n_in_pres
                                 jphase_real = jphase-(jpres-1)*Mdims%n_in_pres
                                 if ( iphase_real == jphase_real ) then
                                    call assign_val(GAMMA_PRES_ABS(IPHASE,JPHASE,:),gamma%val)
                                 end if
                             end do
                         end if
                     end do
                 end do
             end do
             GAMMA_PRES_ABS_NANO = GAMMA_PRES_ABS
             IF ( IGOT_THETA_FLUX == 1 ) THEN ! We have already put density in theta...
                 ! use DEN=1 because the density is already in the theta variables
                 ALLOCATE( DEN_ALL( Mdims%nphase, Mdims%cv_nonods )); DEN_ALL = 1.
                 ALLOCATE( DENOLD_ALL( Mdims%nphase, Mdims%cv_nonods )); DENOLD_ALL = 1.
             ELSE
                 DEN_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedDensity" )
                 DENOLD_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldDensity" )
                 if (is_flooding) then
                    !For Flooding the densities not in the pipes have to be equal to unity
                     ALLOCATE( DEN_ALL( Mdims%nphase, Mdims%cv_nonods ))
                     ALLOCATE( DENOLD_ALL( Mdims%nphase, Mdims%cv_nonods ))
                     do iphase = 1, mdims%n_in_pres
                         DEN_ALL(iphase,:) = 1.0
                         DENOLD_ALL(iphase,:) = 1.0
                     end do
                     do iphase = mdims%n_in_pres+1, mdims%nphase
                         DEN_ALL(iphase,:) = DEN_ALL2%VAL( 1, iphase, : )
                         DENOLD_ALL(iphase,:) = DENOLD_ALL2%VAL( 1, iphase, : )
                     end do
                 else
                     DEN_ALL => DEN_ALL2%VAL( 1, :, : ) ; DENOLD_ALL => DENOLD_ALL2%VAL( 1, :, : )
                 end if
             END IF
             TDIFFUSION = 0.0
             Mdisopt%v_beta = 1.0
             IGOT_THERM_VIS=0
             ALLOCATE( THERM_U_DIFFUSION(Mdims%ndim,Mdims%ndim,Mdims%nphase,Mdims%mat_nonods*IGOT_THERM_VIS ) )
             ALLOCATE( THERM_U_DIFFUSION_VOL(Mdims%nphase,Mdims%mat_nonods*IGOT_THERM_VIS ) )
             tracer=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
             velocity=>extract_tensor_field(packed_state,"PackedVelocity")
             density=>extract_tensor_field(packed_state,"PackedDensity")
             sparsity=>extract_csr_sparsity(packed_state,"ACVSparsity")
             !This logical is used to loop over the saturation equation until the functional
             !explained in function get_Convergence_Functional has been reduced enough
             satisfactory_convergence = .false.
             updating = 0.0
             !We store the convergence of the previous FPI to compare with
             Previous_convergence = backtrack_or_convergence!<== deprecated?
             its = 1; useful_sats = 1;
             if (resold < 0 ) res = huge(res)!<=initialize res once
             call allocate(Mmat%CV_RHS,Mdims%nphase,tracer%mesh,"RHS")
             !Allocate residual, to compute the residual
             if (backtrack_par_factor < 1.01) call allocate(residual,Mdims%nphase,tracer%mesh,"residual")
             Loop_NonLinearFlux: do while (.not. satisfactory_convergence)
                 !If I don't re-allocate this field every iteration, PETSC complains(sometimes),
                 !it works, but it complains...
                 call allocate_global_multiphase_petsc_csr(Mmat%petsc_ACV,sparsity,tracer)
                 !Assemble the matrix and the RHS
                 !before the sprint in this call the small_acv sparsity was passed as cmc sparsity...
                 call CV_ASSEMB( state, packed_state, &
                     Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,upwnd,&
                     tracer, velocity, density, &
                     DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, GAMMA_PRES_ABS, GAMMA_PRES_ABS_NANO, &
                     INV_B, MASS_PIPE, MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE,&
                     DEN_ALL, DENOLD_ALL, &
                     TDIFFUSION, IGOT_THERM_VIS, THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL,&
                     Mdisopt%v_disopt, Mdisopt%v_dg_vel_int_opt, DT, Mdisopt%v_theta, SECOND_THETA, Mdisopt%v_beta, &
                     SUF_SIG_DIAGTEN_BC, &
                     DERIV%val(1,:,:), P, &
                     V_SOURCE, V_ABSORB, VOLFRA_PORE, &
                     GETCV_DISC, GETCT, &
                     IGOT_T2, igot_theta_flux, GET_THETA_FLUX, Mdisopt%volfra_get_theta_flux, &
                     THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
                     MEAN_PORE_CV, &
                     mass_Mn_pres, THERMAL, RETRIEVE_SOLID_CTY, &
                     .false.,  mass_Mn_pres, &
                     mass_ele_transp,IDs_ndgln, &          !Capillary variables
                     OvRelax_param = OvRelax_param, Phase_with_Pc = Phase_with_Pc,&
                     Courant_number = Courant_number)
                 !Solve the system
                 vtracer=as_vector(tracer,dim=2)
                 !If using FPI with backtracking
                 if (backtrack_par_factor < 1.01) then
                     !Backup of the saturation field, to adjust the solution
                     sat_bak = satura
                     !If using ADAPTIVE FPI with backtracking
                     if (backtrack_par_factor < 0) then
                         if (Auto_max_backtrack) then!The maximum backtracking factor depends on the Courant number
                           call auto_backtracking(backtrack_par_factor, courant_number, first_time_step, nonlinear_iteration)
                         end if

                         !Calculate the actual residual using a previous backtrack_par
                         call mult(residual, Mmat%petsc_ACV, vtracer)
                         !Calculate residual
                         residual%val = Mmat%CV_RHS%val - residual%val
                         resold = res; res = 0
                         do iphase = 1, Mdims%nphase
                             aux = sqrt(dot_product(residual%val(iphase,:),residual%val(iphase,:)))/ dble(size(residual%val,2))
                             if (aux > res) res = aux
                         end do
                         !We use the highest residual across the domain
                         if (IsParallel()) call allmax(res)
                         if (its==1) first_res = res!Variable to check total convergence of the SFPI method
                     end if
                 end if
                 call zero(vtracer)
                 call zero_non_owned(Mmat%CV_RHS)
                 call petsc_solve(vtracer,Mmat%petsc_ACV,Mmat%CV_RHS,trim(option_path))
                 !Set to zero the fields
                 call zero(Mmat%CV_RHS)
                 call deallocate(Mmat%petsc_ACV)
                 !For non-porous media make sure all the phases sum to one
                 if (.not. is_porous_media) then
                    call non_porous_ensure_sum_to_one(packed_state)
                    if (is_flooding .and. Mdims%n_in_pres > 1) then
                        !The real domain can only have water, air is automatically removed from the system
                        vtracer%val(1,:) = 1.0; vtracer%val(2,:) = 0.0
                    end if
                 end if
                 !Correct the solution obtained to make sure we are on track towards the final solution
                 if (backtrack_par_factor < 1.01) then
                     !If convergence is not good, then we calculate a new saturation using backtracking
                     if (.not. satisfactory_convergence) then
                         !Calculate a backtrack_par parameter and update saturation with that parameter, ensuring convergence
                         call FPI_backtracking(packed_state, sat_bak, backtrack_sat, backtrack_par_factor, IDs2CV_ndgln,&
                             Previous_convergence, satisfactory_convergence, new_backtrack_par, its, nonlinear_iteration,&
                             useful_sats,res, res/resold, first_res, Mdims%npres)
                         !Store the accumulated updated done
                         updating = updating + new_backtrack_par
                         !If the backtrack_par factor is not adaptive, then, just one iteration
                         if (backtrack_par_factor > 0) then
                             satisfactory_convergence = .true.
                             exit Loop_NonLinearFlux
                         end if
                         !This have to be consistent between processors
                         if (IsParallel())  call alland(satisfactory_convergence)
                         !If looping again, recalculate
                         if (.not. satisfactory_convergence) then
                             !Store old saturation to fully undo an iteration if it is very divergent
                             backtrack_sat = sat_bak
                             !Velocity is recalculated through updating the sigmas
                             call Calculate_PorousMedia_AbsorptionTerms( state, packed_state, Mdims, CV_funs, CV_GIdims, &
                                   Mspars, ndgln, upwnd, suf_sig_diagten_bc, ids_ndgln, IDs2CV_ndgln, Quality_list )
                             !Also recalculate the Over-relaxation parameter
                             call getOverrelaxation_parameter(packed_state, Mdims, ndgln, OvRelax_param, Phase_with_Pc, IDs2CV_ndgln)
                         else
                             exit Loop_NonLinearFlux
                         end if
                     end if
                 else!Just one iteration
                     exit Loop_NonLinearFlux
                 end if
                 its = its + 1
                 useful_sats = useful_sats + 1
             END DO Loop_NonLinearFlux
             !Store the final accumulated backtrack_par_factor to properly calculate the convergence functional
             if (backtrack_par_factor < 1.01) then
                 !Final effective backtrack_par to calculate properly the non linear convergence is:
                 backtrack_or_convergence = updating
             else
                 backtrack_or_convergence = 1
             end if
             !Make sure the parameter is consistent between cpus
             if (IsParallel()) call allmin(backtrack_or_convergence)
             DEALLOCATE( mass_mn_pres )
             DEALLOCATE( Mmat%CT )
             DEALLOCATE( DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, GAMMA_PRES_ABS, GAMMA_PRES_ABS_NANO )
             DEALLOCATE( TDIFFUSION )
             IF ( IGOT_T2 == 1 ) THEN
                 DEALLOCATE( T2 )
                 DEALLOCATE( T2OLD )
             END IF
             DEALLOCATE( THETA_GDIFF )
             call deallocate(Mmat%CV_RHS)
             if (backtrack_par_factor < 1.01) call deallocate(residual)
             !Deallocate pointers only if not pointing to something in packed state
             if (IGOT_THETA_FLUX == 1 .or. is_flooding) then
                 deallocate(DEN_ALL, DENOLD_ALL)
             end if
             nullify(DEN_ALL); nullify(DENOLD_ALL)
             ewrite(3,*) 'Leaving VOLFRA_ASSEM_SOLVE'

         contains

        subroutine non_porous_ensure_sum_to_one(packed_state)
            !This subroutines eliminates the oscillations in the saturation that are bigger than a
            !certain tolerance and also sets the saturation to be between bounds
            Implicit none
            !Global variables
            type( state_type ), intent(inout) :: packed_state
            !Local variables
            integer :: iphase, cv_nod, i_start, i_end, ipres
            real :: correction, sum_of_phases
            real, dimension(:,:), pointer :: satura

            call get_var_from_packed_state(packed_state, PhaseVolumeFraction = satura)
            !Impose sat to be between bounds for blocks of saturations (this is for multiple pressure, otherwise there is just one block)
            do ipres = 1, Mdims%npres
                i_start = 1 + (ipres-1) * Mdims%nphase/Mdims%npres
                i_end = ipres * Mdims%nphase/Mdims%npres
                !Set saturation to be between bounds
                do cv_nod = 1, size(satura,2 )
                    sum_of_phases = sum(satura(i_start:i_end, cv_nod))
                    correction = (1.0 - sum_of_phases)
                    !Spread the error to all the phases weighted by their presence in that CV
                    !Increase the range to look for solutions by allowing oscillations below 0.1 percent
                    if (abs(correction) > 1d-3) satura(i_start:i_end, cv_nod) = (satura(i_start:i_end, cv_nod) * (1.0 + correction/sum_of_phases))
                    !Make sure saturation is between bounds after the modification
                    do iphase = i_start, i_end
                        satura(iphase,cv_nod) =  min(max(0., satura(iphase,cv_nod)),1.0)
                    end do
                end do
            end do
        end subroutine non_porous_ensure_sum_to_one



         subroutine auto_backtracking(backtrack_par_factor, courant_number, first_time_step, nonlinear_iteration)
            !The maximum backtracking factor is calculated based on the Courant number and physical effects ocurring in the domain
            implicit none
            real, intent(inout) :: backtrack_par_factor
            real, intent(in) :: courant_number
            logical, intent(in) :: first_time_step
            integer, intent(in) :: nonlinear_iteration
            !Local variables
            real :: physics_adjustment
            logical, save :: Readed_options = .false.
            logical, save :: gravity, cap_pressure, compositional, many_phases, black_oil, ov_relaxation, one_phase


            if (.not.readed_options) then
                !We read the options just once, and then they are stored as logicals
                gravity = have_option("/physical_parameters/gravity")
                if (have_option_for_any_phase("/multiphase_properties/capillary_pressure", Mdims%nphase)) then
                    cap_pressure = .true.
                else
                    cap_pressure = .false.
                end if
                compositional = Mdims%ncomp > 0
                many_phases = Mdims%n_in_pres > 2
                black_oil = have_option( "/physical_parameters/black-oil_PVT_table")
                !Positive effects on the convergence !Need to check for shock fronts...
                if (have_option_for_any_phase("/multiphase_properties/Sat_overRelax", Mdims%nphase)) then
                    ov_relaxation = .true.
                else
                    ov_relaxation = .false.
                end if
                one_phase = (Mdims%n_in_pres == 1)
                readed_options = .true.
            end if

            !Depending on the active physics, the problem is more complex and requires more relaxation
            physics_adjustment = 1.
            !Negative effects on the convergence
            if (gravity) physics_adjustment = physics_adjustment * 1.5
            if (cap_pressure) physics_adjustment = physics_adjustment * 2.0
            if (compositional) physics_adjustment = physics_adjustment * 1.5
            if (many_phases) physics_adjustment = physics_adjustment * 1.2
            !For the first two non-linear iterations, it has to re-adjust, as the gas and oil are again mixed
            if (black_oil .and. nonlinear_iteration <= 2) physics_adjustment = physics_adjustment * 2.

            !Positive effects on the convergence !Need to check for shock fronts...
            if (ov_relaxation) physics_adjustment = physics_adjustment * 0.8
            if (one_phase) physics_adjustment = physics_adjustment * 0.5

            !For the time being, it is based on this simple table
            if (Courant_number * physics_adjustment > 40.) then
                backtrack_par_factor = -0.1
            else if (Courant_number * physics_adjustment > 25.) then
                backtrack_par_factor = -0.15
            else if (Courant_number * physics_adjustment > 15.) then
                backtrack_par_factor = -0.2
            else if (Courant_number * physics_adjustment > 8.) then
                backtrack_par_factor = -0.33
            else if (Courant_number * physics_adjustment > 5.) then
                backtrack_par_factor = -0.5
            else if (Courant_number * physics_adjustment > 1) then
                backtrack_par_factor = -0.8
            else
                backtrack_par_factor = -1.
            end if
            !For the first calculation, the Courant number is usually zero, hence we force a safe value here
            if (first_time_step .and. nonlinear_iteration == 1) backtrack_par_factor = -0.05
            !Use the most restrictive value across all the processors
            if (IsParallel()) call allmin(backtrack_par_factor)




         end subroutine auto_backtracking



    end subroutine VolumeFraction_Assemble_Solve


   SUBROUTINE FORCE_BAL_CTY_ASSEM_SOLVE( state, packed_state,  &
        Mdims, CV_GIdims, FE_GIdims, CV_funs, FE_funs, Mspars, ndgln, Mdisopt, Mmat,upwnd, &
        velocity, pressure, &
        DT, NLENMCY, &!sprint_to_do NLENMCY in Mdims?
        SUF_SIG_DIAGTEN_BC, &
        V_SOURCE, VOLFRA_PORE, &
        !THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL, &
        IGOT_THETA_FLUX, &
        THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, &
        IDs_ndgln, calculate_mass_delta )
        IMPLICIT NONE
        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_GI_dimensions), intent(in) :: CV_GIdims, FE_GIdims
        type(multi_shape_funs), intent(inout) :: CV_funs
        type(multi_shape_funs), intent(in) :: FE_funs
        type (multi_sparsities), intent(in) :: Mspars
        type(multi_ndgln), intent(in) :: ndgln
        type (multi_discretization_opts) :: Mdisopt
        type (multi_matrices), intent(inout) :: Mmat
        type (porous_adv_coefs), intent(inout) :: upwnd
        type( tensor_field ), intent(inout) :: velocity
        type( tensor_field ), intent(inout) :: pressure
        INTEGER, intent( in ) :: IGOT_THETA_FLUX, NLENMCY
        INTEGER, DIMENSION(  :  ), intent( in ) :: IDs_ndgln
        REAL, DIMENSION(  : , :  ), intent( in ) :: SUF_SIG_DIAGTEN_BC
        REAL, intent( in ) :: DT
        REAL, DIMENSION(  :, :  ), intent( in ) :: V_SOURCE
        !REAL, DIMENSION(  : ,  : ,: ), intent( in ) :: V_ABSORB
        REAL, DIMENSION(  :, :  ), intent( in ) :: VOLFRA_PORE
        REAL, DIMENSION( : ,  :  ), intent( inout ) :: &
        THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
        ! Local Variables
        LOGICAL, PARAMETER :: PIPES_1D = .TRUE. ! Switch on 1D pipe modelling
        LOGICAL, PARAMETER :: GLOBAL_SOLVE = .FALSE.
        ! If IGOT_CMC_PRECON=1 use a sym matrix as pressure preconditioner,=0 else CMC as preconditioner as well.
        INTEGER :: IGOT_CMC_PRECON
! Gidaspow model B - can use conservative from of
        LOGICAL :: SOLID_FLUID_MODEL_B = .TRUE.
! switch on solid fluid coupling (THE ONLY SWITCH THAT NEEDS TO BE SWITCHED ON FOR SOLID-FLUID COUPLING)...
        LOGICAL :: RETRIEVE_SOLID_CTY = .FALSE.
! got_free_surf - INDICATED IF WE HAVE A FREE SURFACE - TAKEN FROM DIAMOND EVENTUALLY...
        LOGICAL :: got_free_surf
        character( len = option_path_len ) :: opt, bc_type
        REAL, DIMENSION( : ), allocatable :: &
        MCY_RHS, MCY, &
        MASS_MN_PRES, MASS_SUF, MASS_CV, UP, &
        UP_VEL
        REAL, DIMENSION( :, : ), allocatable :: DIAG_SCALE_PRES
        real, dimension(:,:,:), allocatable :: velocity_absorption, U_SOURCE_CV_ALL
        real, dimension(:,:,:,:), allocatable :: UDIFFUSION_ALL

        real, dimension(:,:) :: calculate_mass_delta
        type( multi_field ) :: UDIFFUSION_VOL_ALL, U_SOURCE_ALL   ! NEED TO ALLOCATE THESE - SUBS TO DO THIS ARE MISSING... - SO SET 0.0 FOR NOW

        type( multi_field ) :: UDIFFUSION_ALL2




        REAL, DIMENSION(  :, :, :  ), allocatable :: temperature_absorption, U_ABSORBIN
        REAL, DIMENSION( :, :, : ), allocatable :: DIAG_SCALE_PRES_COUP, GAMMA_PRES_ABS, GAMMA_PRES_ABS_NANO, INV_B, CMC_PRECON
        REAL, DIMENSION( : ), ALLOCATABLE :: MASS_PIPE, MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE
        REAL, DIMENSION( :, :, : ), allocatable :: DU_VEL, U_RHS_CDP2
        INTEGER :: CV_NOD, COUNT, CV_JNOD, IPHASE, JPHASE, ndpset, i
        LOGICAL :: JUST_BL_DIAG_MAT, LINEARISE_DENSITY, diag, SUF_INT_MASS_MATRIX
        INTEGER :: stat
        !Re-scale parameter can be re-used
        real, save :: rescaleVal = -1.0
        !CMC using petsc format
        type(petsc_csr_matrix)::  CMC_petsc
        !TEMPORARY VARIABLES, ADAPT FROM OLD VARIABLES TO NEW
        INTEGER :: IPRES, JPRES, iphase_real, jphase_real
        REAL, DIMENSION( :, : ), allocatable :: UDEN_ALL, UDENOLD_ALL, UDEN3
        REAL, DIMENSION( :, : ), allocatable :: rhs_p2, sigma
        REAL, DIMENSION( :, : ), pointer :: DEN_ALL, DENOLD_ALL
        type( tensor_field ), pointer :: u_all2, uold_all2, den_all2, denold_all2, tfield, den_all3
        type( tensor_field ), pointer :: p_all, pold_all, cvp_all, deriv, PorousMedia_AbsorptionTerm, Flooding_AbsorptionTerm
        type( vector_field ), pointer :: x_all2
        type( scalar_field ), pointer ::  sf, soldf, gamma
        type( vector_field ) :: packed_vel, rhs
        type( vector_field ) :: deltap, rhs_p
        type(tensor_field) :: cdp_tensor
        type( csr_sparsity ), pointer :: sparsity
        type(halo_type), pointer :: halo
        logical :: cty_proj_after_adapt, high_order_Ph, symmetric_P, boussinesq, fem_density_buoyancy
        logical :: EXPLICIT_PIPES2
        REAL, DIMENSION ( :, :, : ), pointer :: SUF_P_BC_ALL
        INTEGER, DIMENSION ( 1, Mdims%npres, surface_element_count(pressure) ) :: WIC_P_BC_ALL
        type( tensor_field ) :: pressure_BCs
        integer :: IGOT_THERM_VIS, ierr
        real, dimension(:,:), allocatable :: THERM_U_DIFFUSION_VOL
        real, dimension(:,:,:,:), allocatable :: THERM_U_DIFFUSION
        !!$ Variables used in the diffusion-like term: capilarity and surface tension:
        type( tensor_field ), pointer :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD
        INTEGER :: IPLIKE_GRAD_SOU
        !!$ magma stuff -- to be deleted shortly
        integer :: idim, idx1, idx2, ndim
        type( scalar_field ), pointer :: beta
        !!$ end of magma stuff


        ! if q scheme allocate a field in state and use pointers..
        IGOT_THERM_VIS=0
        ALLOCATE( THERM_U_DIFFUSION(Mdims%ndim,Mdims%ndim,Mdims%nphase,Mdims%mat_nonods*IGOT_THERM_VIS ) )
        ALLOCATE( THERM_U_DIFFUSION_VOL(Mdims%nphase,Mdims%mat_nonods*IGOT_THERM_VIS ) )
        EXPLICIT_PIPES2 = .true.
        deriv => extract_tensor_field( packed_state, "PackedDRhoDPressure" )
        high_order_Ph = have_option( "/physical_parameters/gravity/hydrostatic_pressure_solver" )
        symmetric_P = have_option( "/material_phase[0]/scalar_field::Pressure/prognostic/symmetric_P" )
        cty_proj_after_adapt = have_option( "/mesh_adaptivity/hr_adaptivity/project_continuity" )
        boussinesq = have_option( "/material_phase[0]/vector_field::Velocity/prognostic/equation::Boussinesq" )
        fem_density_buoyancy = have_option( "/physical_parameters/gravity/fem_density_buoyancy" )
        got_free_surf = .false.
        do i = 1, get_boundary_condition_count(pressure)
           call get_boundary_condition( pressure, i, type=bc_type )
           if ( trim( bc_type ) == "freesurface" ) then
              got_free_surf = .true.
              exit
           end if
        end do
        IGOT_CMC_PRECON = 0
        if ( symmetric_P ) IGOT_CMC_PRECON = 1
        !sprint_to_do!this looks like a place than can be easily optimized
        ALLOCATE( UDEN_ALL( Mdims%nphase, Mdims%cv_nonods ), UDENOLD_ALL( Mdims%nphase, Mdims%cv_nonods ) )
        UDEN_ALL = 0.; UDENOLD_ALL = 0.
        ewrite(3,*) 'In FORCE_BAL_CTY_ASSEM_SOLVE'
        ALLOCATE( Mmat%CT( Mdims%ndim, Mdims%nphase, Mspars%CT%ncol )) ; Mmat%CT=0.
        call allocate(Mmat%CT_RHS,Mdims%npres,pressure%mesh,"Mmat%CT_RHS")
        ALLOCATE( Mmat%U_RHS( Mdims%ndim, Mdims%nphase, Mdims%u_nonods )) ; Mmat%U_RHS=0.
        ALLOCATE( DIAG_SCALE_PRES( Mdims%npres,Mdims%cv_nonods )) ; DIAG_SCALE_PRES=0.
        ALLOCATE(DIAG_SCALE_PRES_COUP(Mdims%npres,Mdims%npres,Mdims%cv_nonods),GAMMA_PRES_ABS(Mdims%nphase,Mdims%nphase,Mdims%cv_nonods), &
            GAMMA_PRES_ABS_NANO(Mdims%nphase,Mdims%nphase,Mdims%cv_nonods),INV_B(Mdims%nphase,Mdims%nphase,Mdims%cv_nonods))
        allocate(MASS_PIPE(Mdims%cv_nonods), MASS_CVFEM2PIPE(Mspars%CMC%ncol), MASS_PIPE2CVFEM(Mspars%CMC%ncol), MASS_CVFEM2PIPE_TRUE(Mspars%CMC%ncol))
        ALLOCATE( CMC_PRECON( Mdims%npres, Mdims%npres, Mspars%CMC%ncol*IGOT_CMC_PRECON)) ; IF(IGOT_CMC_PRECON.NE.0) CMC_PRECON=0.
        ALLOCATE( MASS_MN_PRES( Mspars%CMC%ncol )) ; MASS_MN_PRES=0.
        ALLOCATE( MASS_CV( Mdims%cv_nonods )) ; MASS_CV=0.
        gamma=>extract_scalar_field(state(1),"Gamma1",stat)
        GAMMA_PRES_ABS = 0.0
        do ipres = 1, Mdims%npres
           do iphase = 1+(ipres-1)*Mdims%n_in_pres, ipres*Mdims%n_in_pres
              do jpres = 1, Mdims%npres
                 if ( ipres /= jpres ) then
                    do jphase = 1+(jpres-1)*Mdims%n_in_pres, jpres*Mdims%n_in_pres
                       iphase_real = iphase-(ipres-1)*Mdims%n_in_pres
                       jphase_real = jphase-(jpres-1)*Mdims%n_in_pres
                       if ( iphase_real == jphase_real ) then
                          call assign_val(GAMMA_PRES_ABS(IPHASE,JPHASE,:),gamma%val)
                       end if
                    end do
                 end if
              end do
           end do
        end do
        GAMMA_PRES_ABS_NANO = GAMMA_PRES_ABS
        ! this can be optimised in the future...

        U_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedVelocity" )
        UOLD_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldVelocity" )
        X_ALL2 => EXTRACT_VECTOR_FIELD( PACKED_STATE, "PressureCoordinate" )
        P_ALL => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedFEPressure" )
        CVP_ALL => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedCVPressure" )

        linearise_density = have_option( '/material_phase[0]/linearise_density' )

        DEN_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedDensity" )
        DENOLD_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldDensity" )
        DEN_ALL(1:, 1:) => DEN_ALL2%VAL( 1, :, : )
        DENOLD_ALL(1:, 1:) => DENOLD_ALL2%VAL( 1, :, : )

        call allocate(deltaP,Mdims%npres,pressure%mesh,"DeltaP")
        call allocate(rhs_p,Mdims%npres,pressure%mesh,"PressureCorrectionRHS")
        Mmat%NO_MATRIX_STORE = ( Mspars%DGM_PHA%ncol <= 1 )
        JUST_BL_DIAG_MAT = .false.
        IF (.not. ( JUST_BL_DIAG_MAT .OR. Mmat%NO_MATRIX_STORE ) ) then
           sparsity=>extract_csr_sparsity(packed_state,"MomentumSparsity")
           Mmat%DGM_PETSC = allocate_momentum_matrix(sparsity,velocity)
        end IF

        !Calculate gravity source terms
        allocate(U_SOURCE_CV_ALL(Mdims%ndim, Mdims%nphase, Mdims%cv_nonods))
        U_SOURCE_CV_ALL=0.0
        if ( is_porous_media )then
           UDEN_ALL=0.0; UDENOLD_ALL=0.0
           call calculate_u_source_cv( state, Mdims%cv_nonods, Mdims%ndim, Mdims%nphase, DEN_ALL, U_SOURCE_CV_ALL )
        else
           if ( linearise_density ) then
              call linearise_field( DEN_ALL2, UDEN_ALL )
              call linearise_field( DENOLD_ALL2, UDENOLD_ALL )
           else
              UDEN_ALL = DEN_ALL2%VAL( 1, :, : )
              UDENOLD_ALL = DENOLD_ALL2%VAL( 1, :, : )
           end if
           allocate(  uden3( Mdims%nphase, Mdims%cv_nonods )  )
           if ( fem_density_buoyancy ) then
              DEN_ALL3 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedFEDensity" )
              UDEN3 = DEN_ALL3 % VAL( 1, :, : )
           else
              UDEN3 = uden_all
           end if
           if ( have_option( "/physical_parameters/gravity/hydrostatic_pressure_solver" ) ) UDEN3 = 0.0
           call calculate_u_source_cv( state, Mdims%cv_nonods, Mdims%ndim, Mdims%nphase, uden3, U_SOURCE_CV_ALL )
           deallocate( uden3 )
           if ( boussinesq ) then
              UDEN_ALL=1.0; UDENOLD_ALL=1.0
           end if
        end if

        if ( have_option( '/blasting' ) ) then
            RETRIEVE_SOLID_CTY = .true.
            call get_option( '/blasting/Gidaspow_model', opt )
            if ( trim( opt ) == "A" ) SOLID_FLUID_MODEL_B = .false.
        end if

        IF(RETRIEVE_SOLID_CTY) THEN
        ! if model B and solid-fluid coupling:
           sf => EXTRACT_SCALAR_FIELD( PACKED_STATE, "SolidConcentration" )
           soldf => EXTRACT_SCALAR_FIELD( PACKED_STATE, "OldSolidConcentration" )
           IF(SOLID_FLUID_MODEL_B) THEN ! Gidaspow model B - can use conservative from of momentum
              DO IPHASE=1,Mdims%nphase
                 UDEN_ALL(IPHASE,:) = UDEN_ALL(IPHASE,:) * ( 1. - sf%val )
                 UDENOLD_ALL(IPHASE,:) = UDENOLD_ALL(IPHASE,:) * ( 1. - soldf%val )
              END DO
           ENDIF
        ENDIF
        allocate(UDIFFUSION_ALL(Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%mat_nonods))
        ! calculate the viscosity for the momentum equation... (uDiffusion is initialized inside)
        call calculate_viscosity( state, Mdims, ndgln, UDIFFUSION_ALL, UDIFFUSION_ALL2 )
        !UDIFFUSION_VOL_ALL = 0.

        allocate(velocity_absorption(Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase, Mdims%mat_nonods))
        ! define velocity_absorption here...
        velocity_absorption=0.0
        ! update velocity absorption
        call update_velocity_absorption( state, Mdims%ndim, Mdims%nphase, velocity_absorption )
        call update_velocity_absorption_coriolis( state, Mdims%ndim, Mdims%nphase, velocity_absorption )


        if ( have_option( "/magma" ) ) then
           ndim = Mdims%ndim

           beta => extract_scalar_field( state( 1 ), "beta" )

           iphase=1 ; jphase=1
           do idim = 1, ndim
              idx1=idim+(iphase-1)*ndim ; idx2=idim+(jphase-1)*ndim
              velocity_absorption( idx1, idx2, : ) = beta%val
           end do

           iphase=1 ; jphase=2
           do idim = 1, ndim
              idx1=idim+(iphase-1)*ndim ; idx2=idim+(jphase-1)*ndim
              velocity_absorption( idx1, idx2, : ) = -beta%val
           end do

           iphase=2 ; jphase=1
           do idim = 1, ndim
              idx1=idim+(iphase-1)*ndim ; idx2=idim+(jphase-1)*ndim
              velocity_absorption( idx1, idx2, : ) = -beta%val
           end do

           iphase=2 ; jphase=2
           do idim = 1, ndim
              idx1=idim+(iphase-1)*ndim ; idx2=idim+(jphase-1)*ndim
              velocity_absorption( idx1, idx2, : ) = beta%val
           end do

        end if


        ! open the boiling test for two phases-gas and liquid
        if (have_option('/boiling')) then
           allocate( temperature_absorption( Mdims%nphase, Mdims%nphase, Mdims%cv_nonods ) )
           call boiling( state, packed_state, Mdims%cv_nonods, Mdims%mat_nonods, Mdims%nphase, Mdims%ndim, &
                velocity_absorption, temperature_absorption )
           deallocate( temperature_absorption )
        end if

        ! update velocity source
        call update_velocity_source( state, Mdims%ndim, Mdims%nphase, u_source_all )

        PorousMedia_AbsorptionTerm => extract_tensor_field( packed_state, "PorousMedia_AbsorptionTerm", stat )
        if ( stat == 0 ) velocity_absorption = velocity_absorption + PorousMedia_AbsorptionTerm%val

        Flooding_AbsorptionTerm => extract_tensor_field( packed_state, "Flooding_AbsorptionTerm", stat )
        if (stat == 0) velocity_absorption = velocity_absorption + Flooding_AbsorptionTerm%val

        !Check if as well the Mass matrix
        SUF_INT_MASS_MATRIX = .false.!= have_option( '/material_phase[0]/scalar_field::Pressure/prognostic/CV_P_matrix/Suf_mass_matrix' )

        !Allocation of storable matrices !SPRINT_TO_DO; move as an initialization subroutine
        if (.not.Mmat%Stored) then
             if (Mmat%CV_pressure) then
                allocate(Mmat%C_CV(Mdims%ndim, Mdims%nphase, Mspars%C%ncol)); Mmat%C_CV = 0.
            else!allocate C
                allocate(Mmat%C(Mdims%ndim, Mdims%nphase, Mspars%C%ncol)); Mmat%C = 0.
            end if
            allocate( Mmat%PIVIT_MAT( Mdims%ndim * Mdims%nphase * Mdims%u_nloc, Mdims%ndim * Mdims%nphase * Mdims%u_nloc, Mdims%totele ) ); Mmat%PIVIT_MAT=0.0
        end if

        !If it is not porous media or there are more than one pressure, PIVIT_MAT needs to be recalculated, hence we set it to zero
        if (.not.is_porous_media .or. Mdims%npres > 1) then
            !Check if it requires allocation
            if (.not.associated(Mmat%PIVIT_MAT)) allocate( Mmat%PIVIT_MAT( Mdims%ndim * Mdims%nphase * Mdims%u_nloc, Mdims%ndim * Mdims%nphase * Mdims%u_nloc, Mdims%totele ) )
            Mmat%PIVIT_MAT=0.0
        end if

        if( have_option_for_any_phase( '/multiphase_properties/capillary_pressure', Mdims%nphase ) )then
            !The first time (itime/=1 .or. its/=1) we use CVSat since FESAt is not defined yet
            call calculate_capillary_pressure(packed_state, .false., &!sprint_to_do; shouldn't this flag change after the first non-linear iteration?
                ndgln%cv, ids_ndgln, Mdims%totele, Mdims%cv_nloc)
        end if

!        IF(got_free_surf .or. Mmat%CV_pressure) THEN!sprint_to_do: pressure bcs as vel bcs, remove
        IF(got_free_surf) THEN
           ALLOCATE( MASS_SUF( Mspars%CMC%ncol )) ; MASS_SUF=0.
        ELSE
           ALLOCATE( MASS_SUF( 1 )) ; MASS_SUF=0.
        ENDIF

        ! calculate surface tension
        !!$ extended to surface tension -like term.
        IPLIKE_GRAD_SOU = 0
        if( have_option_for_any_phase( '/is_multiphase_component/surface_tension', Mdims%nphase+Mdims%ncomp ) ) then
            PLIKE_GRAD_SOU_GRAD => EXTRACT_TENSOR_FIELD( PACKED_STATE, "SurfaceTensionGrad" )
            PLIKE_GRAD_SOU_COEF => EXTRACT_TENSOR_FIELD( PACKED_STATE, "SurfaceTensionCoef" )
            CALL CALCULATE_SURFACE_TENSION_NEW( state, packed_state, Mdims, Mspars, ndgln, Mdisopt, &
                PLIKE_GRAD_SOU_COEF%val, PLIKE_GRAD_SOU_GRAD%val, IPLIKE_GRAD_SOU)
        end if

        ! solid pressure term - use the surface tension code
        if ( have_option( "/magma" ) ) IPLIKE_GRAD_SOU = 2


        IF ( GLOBAL_SOLVE ) then
            !Prepare memory specific for this
            ALLOCATE( MCY_RHS( Mdims%ndim * Mdims%nphase * Mdims%u_nonods + Mdims%cv_nonods )) ; MCY_RHS=0.
            ALLOCATE( MCY( Mspars%MCY%ncol )) ; MCY=0.
            ALLOCATE( UP( NLENMCY )) ; UP=0.
        end if
        CALL CV_ASSEMB_FORCE_CTY( state, packed_state, &
            Mdims, CV_GIdims, FE_GIdims, CV_funs, FE_funs, Mspars, ndgln, Mdisopt, Mmat,upwnd, &
             velocity, pressure, &
            X_ALL2%VAL, velocity_absorption, U_SOURCE_ALL, U_SOURCE_CV_ALL, &
            U_ALL2%VAL, UOLD_ALL2%VAL, &
            P_ALL%VAL, CVP_ALL%VAL, DEN_ALL, DENOLD_ALL, DERIV%val(1,:,:), &
            DT, &
            MASS_MN_PRES, & ! pressure matrix for projection method
            got_free_surf,  MASS_SUF, &
            SUF_SIG_DIAGTEN_BC, &
            V_SOURCE, VOLFRA_PORE, &
            MCY_RHS, DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, GAMMA_PRES_ABS, GAMMA_PRES_ABS_NANO, INV_B, MASS_PIPE, MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE, GLOBAL_SOLVE, &
            MCY, JUST_BL_DIAG_MAT, &
            UDEN_ALL, UDENOLD_ALL, UDIFFUSION_ALL,  UDIFFUSION_VOL_ALL, THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL, &
            IGOT_THETA_FLUX, &
            THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, &
            RETRIEVE_SOLID_CTY, &
            IPLIKE_GRAD_SOU,&
            symmetric_P, boussinesq, IDs_ndgln, calculate_mass_delta)

        deallocate(GAMMA_PRES_ABS, GAMMA_PRES_ABS_NANO, UDIFFUSION_ALL)
        !If pressure in CV then point the FE matrix Mmat%C to Mmat%C_CV
        if ( Mmat%CV_pressure ) Mmat%C => Mmat%C_CV
        if ( Mdims%npres > 1 ) then
           ALLOCATE( SIGMA( Mdims%nphase, Mdims%mat_nonods ) )
           DO IPHASE = 1, Mdims%nphase
              SIGMA( IPHASE, : ) = velocity_absorption( (IPHASE-1)*Mdims%ndim+1, (IPHASE-1)*Mdims%ndim+1, : )
           END DO

           call get_entire_boundary_condition( pressure,&
                ['weakdirichlet','freesurface  '],&
                pressure_BCs, WIC_P_BC_ALL )
           SUF_P_BC_ALL => pressure_BCs%val
           !Introduce well modelling
           CALL MOD_1D_FORCE_BAL_C( STATE, packed_state, Mmat%U_RHS, Mdims, Mspars, associated(Mmat%PIVIT_MAT), &
                                    Mmat%C, ndgln%cv, ndgln%u, ndgln%x, ndgln%mat, Mmat%PIVIT_MAT, &
                                    ndgln%suf_p, WIC_P_BC_ALL, SUF_P_BC_ALL, SIGMA, U_ALL2%VAL, &
                                    U_SOURCE_ALL, U_SOURCE_CV_ALL*0.0 ) ! No sources in the wells for now...
           call deallocate( pressure_BCs )
           DEALLOCATE( SIGMA )
        end if
        deallocate(velocity_absorption, U_SOURCE_CV_ALL)
        IF ( .NOT.GLOBAL_SOLVE ) THEN
            ! form pres eqn.
            if (.not.Mmat%Stored .or. (.not.is_porous_media .or. Mdims%npres > 1))  CALL PHA_BLOCK_INV(&
                                        Mmat%PIVIT_MAT, Mdims%totele, Mdims%u_nloc * Mdims%nphase * Mdims%ndim )
            sparsity=>extract_csr_sparsity(packed_state,'CMCSparsity')
            diag=.true.
            if ( Mdims%npres>1 ) diag=.false.
            call allocate(CMC_petsc,sparsity,[Mdims%npres,Mdims%npres],"CMC_petsc",diag)
            if (associated(pressure%mesh%halos)) then
               halo => pressure%mesh%halos(2)
            else
               nullify(halo)
            end if

            !Form pressure matrix
            CALL COLOR_GET_CMC_PHA( Mdims, Mspars, ndgln, Mmat,&
            DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
            CMC_petsc, CMC_PRECON, IGOT_CMC_PRECON, MASS_MN_PRES, &
            MASS_PIPE, MASS_CVFEM2PIPE, MASS_CVFEM2PIPE_TRUE, &
            got_free_surf,  MASS_SUF, symmetric_P )
        END IF

        Mmat%NO_MATRIX_STORE = ( Mspars%DGM_PHA%ncol <= 1 )
        IF ( GLOBAL_SOLVE ) THEN
            ! Global solve
            IF ( JUST_BL_DIAG_MAT ) THEN
                EWRITE(-1,*) 'OPTION NOT READY YET WITH A GLOBAL SOLVE'
                STOP 8331
            END IF
            UP = 0.0
            CALL multi_solver( MCY, UP, MCY_RHS, &
            Mspars%MCY%fin, Mspars%MCY%col, &
            option_path = '/material_phase[0]/vector_field::Velocity')
            U_ALL2 % val = reshape( UP( 1 : Mdims%u_nonods * Mdims%ndim * Mdims%nphase ), (/ Mdims%ndim, Mdims%nphase, Mdims%u_nonods /) )
            P_ALL % val( 1, 1, : ) = UP( Mdims%u_nonods * Mdims%ndim * Mdims%nphase + 1 : Mdims%u_nonods * Mdims%ndim * Mdims%nphase + Mdims%cv_nonods )
            DEALLOCATE( MCY_RHS ); DEALLOCATE( MCY ); DEALLOCATE( UP )
        ELSE ! solve using a projection method
            call allocate(cdp_tensor,velocity%mesh,"CDP",dim=velocity%dim); call zero(cdp_tensor)
            ! Put pressure in rhs of force balance eqn: CDP = Mmat%C * P
!            CALL C_MULT2( CDP_TENSOR%VAL, P_ALL%val , Mdims%cv_nonods, Mdims%u_nonods, Mdims%ndim, Mdims%nphase, Mmat%C, NCOLC, FINDC, COLC)
            DO IPRES = 1, Mdims%npres
               CALL C_MULT2( CDP_TENSOR%VAL( :, 1+(IPRES-1)*Mdims%n_in_pres : IPRES*Mdims%n_in_pres, : ), P_ALL%val( 1, IPRES, : ), &
                    Mdims%cv_nonods, Mdims%u_nonods, Mdims%ndim, Mdims%n_in_pres, Mmat%C( :, 1+(IPRES-1)*Mdims%n_in_pres : IPRES*Mdims%n_in_pres, : ), Mspars%C%ncol, Mspars%C%fin, Mspars%C%col )
            END DO
            if ( high_order_Ph ) then
               if ( .not. ( after_adapt .and. cty_proj_after_adapt ) ) then
                  allocate (U_ABSORBIN(Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase, Mdims%mat_nonods))
                  call update_velocity_absorption( state, Mdims%ndim, Mdims%nphase, U_ABSORBIN )
                  call update_velocity_absorption_coriolis( state, Mdims%ndim, Mdims%nphase, U_ABSORBIN )
                  call high_order_pressure_solve( Mdims, Mmat%u_rhs, state, packed_state, Mdims%nphase, U_ABSORBIN )
                  deallocate(U_ABSORBIN)
               end if
            end if
            ALLOCATE( UP_VEL( Mdims%ndim * Mdims%nphase * Mdims%u_nonods )) ; UP_VEL = 0.
            IF ( JUST_BL_DIAG_MAT .OR. Mmat%NO_MATRIX_STORE ) THEN
                ALLOCATE( U_RHS_CDP2( Mdims%ndim, Mdims%nphase, Mdims%u_nonods ))
                !For porous media we calculate the velocity as M^-1 * CDP, no solver is needed
                U_RHS_CDP2 = Mmat%U_RHS + CDP_tensor%val
                ! DU = BLOCK_MAT * CDP
                CALL PHA_BLOCK_MAT_VEC_old( UP_VEL, Mmat%PIVIT_MAT, U_RHS_CDP2, Mdims%ndim, Mdims%nphase, &
                Mdims%totele, Mdims%u_nloc, ndgln%u )
                U_ALL2 % VAL = RESHAPE( UP_VEL, (/ Mdims%ndim, Mdims%nphase, Mdims%u_nonods /) )
                deallocate(U_RHS_CDP2)
            ELSE
               call allocate(rhs,product(velocity%dim),velocity%mesh,"RHS")
               rhs%val=RESHAPE( Mmat%U_RHS + CDP_tensor%val, (/ Mdims%ndim * Mdims%nphase , Mdims%u_nonods /) )
               call zero_non_owned(rhs)
               if ( .not. ( after_adapt .and. cty_proj_after_adapt ) ) then
                  call zero(velocity)
                  packed_vel=as_packed_vector(velocity)
                  call petsc_solve( packed_vel, Mmat%DGM_PETSC, RHS )
#ifndef USING_GFORTRAN
                  velocity%val(:,:,:)=reshape(packed_vel%val,[size(velocity%val,1),size(velocity%val,2),size(velocity%val,3)])
#endif
                  U_ALL2 % VAL=velocity%val
               end if
               UP_VEL=[velocity%val]
               call deallocate(Mmat%DGM_PETSC)
               call deallocate(rhs)
               U_ALL2 % VAL = RESHAPE( UP_VEL, (/ Mdims%ndim, Mdims%nphase, Mdims%u_nonods /) )
            END IF
            deallocate( UP_VEL )
IF ( Mdims%npres > 1 .AND. .NOT.EXPLICIT_PIPES2 ) THEN
            if ( .not.symmetric_P ) then ! original
               ALLOCATE ( rhs_p2(Mdims%nphase,Mdims%cv_nonods) ) ; rhs_p2=0.0
               DO IPHASE = 1, Mdims%nphase
                  CALL CT_MULT2( rhs_p2(IPHASE,:), U_ALL2%VAL( :, IPHASE : IPHASE, : ), &
                       Mdims%cv_nonods, Mdims%u_nonods, Mdims%ndim, 1, Mmat%CT( :, IPHASE : IPHASE, : ), Mspars%CT%ncol, Mspars%CT%fin, Mspars%CT%col )
               END DO
            else
               DO IPHASE = 1, Mdims%nphase
                  CALL CT_MULT_WITH_C3( rhs_p2(IPHASE,:), U_ALL2%VAL( :, IPHASE : IPHASE, : ), &
                       Mdims%u_nonods, Mdims%ndim, 1, Mmat%C( :, IPHASE : IPHASE, : ), Mspars%C%ncol, Mspars%C%fin, Mspars%C%col )
               END DO
            end if
            DO CV_NOD = 1, Mdims%cv_nonods
                  rhs_p2(:,CV_NOD) = MATMUL( INV_B(:,:,CV_NOD), rhs_p2(:,CV_NOD) )
            END DO
            DO CV_NOD = 1, Mdims%cv_nonods
               DO IPRES = 1, Mdims%npres
                  rhs_p%val(IPRES,CV_NOD)= SUM( rhs_p2(1+(IPRES-1)*Mdims%n_in_pres : IPRES*Mdims%n_in_pres,CV_NOD) )
               END DO
            END DO
ELSE
            if ( .not.symmetric_P ) then ! original
               DO IPRES = 1, Mdims%npres
                  CALL CT_MULT2( rhs_p%val(IPRES,:), U_ALL2%VAL( :, 1+(IPRES-1)*Mdims%n_in_pres : IPRES*Mdims%n_in_pres, : ), &
                       Mdims%cv_nonods, Mdims%u_nonods, Mdims%ndim, Mdims%n_in_pres, &
                       Mmat%CT( :, 1+(IPRES-1)*Mdims%n_in_pres : IPRES*Mdims%n_in_pres, : ), Mspars%CT%ncol, Mspars%CT%fin, Mspars%CT%col )
               END DO
            else
               DO IPRES = 1, Mdims%npres
                  CALL CT_MULT_WITH_C3( rhs_p%val(IPRES,:), U_ALL2%VAL( :, 1+(IPRES-1)*Mdims%n_in_pres : IPRES*Mdims%n_in_pres, : ), &
                       Mdims%u_nonods, Mdims%ndim, Mdims%n_in_pres, Mmat%C( :, 1+(IPRES-1)*Mdims%n_in_pres : IPRES*Mdims%n_in_pres, : ), Mspars%C%ncol, Mspars%C%fin, Mspars%C%col )
               END DO
            end if
END IF
            rhs_p%val = -rhs_p%val + Mmat%CT_RHS%val
            if(got_free_surf) POLD_ALL => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldFEPressure" )
            ! Matrix vector involving the mass diagonal term
            DO CV_NOD = 1, Mdims%cv_nonods
               DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                  CV_JNOD = Mspars%CMC%col( COUNT )
                  DO IPRES = 1, Mdims%npres
                     IF (( Mdims%npres > 1 ).AND.PIPES_1D) THEN
                        IF(IPRES==1) THEN
                           rhs_p%val( IPRES, CV_NOD ) = rhs_p%val( IPRES, CV_NOD ) &
                           -DIAG_SCALE_PRES( IPRES, CV_NOD ) * MASS_MN_PRES( COUNT ) * P_ALL%VAL( 1, IPRES, CV_JNOD )
                        ELSE
                           rhs_p%val( IPRES, CV_NOD ) = rhs_p%val( IPRES, CV_NOD ) &
                           -DIAG_SCALE_PRES( IPRES, CV_NOD ) * MASS_CVFEM2PIPE_TRUE( COUNT ) * P_ALL%VAL( 1, IPRES, CV_JNOD )
                        ENDIF
                     ELSE
                        rhs_p%val( IPRES, CV_NOD ) = rhs_p%val( IPRES, CV_NOD ) &
                          -DIAG_SCALE_PRES( IPRES, CV_NOD ) * MASS_MN_PRES( COUNT ) * P_ALL%VAL( 1, IPRES, CV_JNOD )
                     ENDIF
                     if ( got_free_surf ) then
                        rhs_p%val( IPRES, CV_NOD ) = rhs_p%val( IPRES, CV_NOD ) &
                             -MASS_SUF( COUNT ) * ( P_ALL%VAL( 1, IPRES, CV_JNOD ) - POLD_ALL%VAL( 1, IPRES, CV_JNOD ) )
                     end if
!                     !For the pressure bcs of CV pressure formulation this is part of the projection method to make them vel bcs
!                     if ( Mmat%CV_pressure) then!sprint_to_do: pressure bcs as vel bcs, remove
!                        rhs_p%val( IPRES, CV_NOD ) = rhs_p%val( IPRES, CV_NOD ) - MASS_SUF( COUNT ) * P_ALL%VAL( 1, IPRES, CV_JNOD )
!                      end if
                     IF ( Mdims%npres > 1 ) THEN
                        DO JPRES = 1, Mdims%npres
                           IF(PIPES_1D) THEN
                              rhs_p%val( IPRES, CV_NOD ) = rhs_p%val( IPRES, CV_NOD ) &
                                -DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_NOD ) * MASS_CVFEM2PIPE( COUNT ) * P_ALL%VAL( 1, JPRES, CV_JNOD )
                           ELSE
                              rhs_p%val( IPRES, CV_NOD ) = rhs_p%val( IPRES, CV_NOD ) &
                                -DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_NOD ) * MASS_MN_PRES( COUNT ) * P_ALL%VAL( 1, JPRES, CV_JNOD )
                           ENDIF
                        END DO
                     END IF
                  END DO
               END DO
            END DO
            call zero_non_owned(rhs_p)
            call get_option( '/material_phase[0]/scalar_field::Pressure/' // &
            'prognostic/reference_node', ndpset, default = 0 )
            if ( ndpset /= 0 ) rhs_p%val( 1, ndpset ) = 0.0
            ! solve for pressure correction DP that is solve CMC*DP=P_RHS...
            ewrite(3,*)'about to solve for pressure'
            !########Solve the system#############
            !Re-scale of the matrix to allow working with small values of sigma
            !this is a hack to deal with bad preconditioners and divide by zero errors.
            if (is_porous_media) then
                !Since we save the parameter rescaleVal, we only do this one time
                if (rescaleVal < 0.) then
                    tfield => extract_tensor_field(packed_state,"Permeability")
                    rescaleVal = minval(tfield%val, MASK = tfield%val > 1d-30)
                    !If it is parallel then we want to be consistent between cpus
                    if (IsParallel()) call allmin(rescaleVal)
                end if
                call scale(cmc_petsc, 1.0/rescaleVal)
                rhs_p%val = rhs_p%val / rescaleVal
                !End of re-scaling
            end if
            call zero(deltaP)
            !Solve the system to obtain dP (difference of pressure)
            call petsc_solve(deltap,cmc_petsc,rhs_p,trim(pressure%option_path))
            P_all % val(1,:,:) = P_all % val(1,:,:) + deltap%val
            call halo_update(p_all)
            call deallocate(rhs_p)
            call deallocate(cmc_petsc)
            ewrite(3,*) 'after pressure solve DP:', minval(deltap%val), maxval(deltap%val)
            ! Use a projection method
            ! CDP = Mmat%C * DP
            !CALL C_MULT2( CDP_tensor%val, deltap%val, Mdims%cv_nonods, Mdims%u_nonods, Mdims%ndim, Mdims%nphase, Mmat%C, Mspars%C%ncol, Mspars%C%fin, Mspars%C%col )
            DO IPRES = 1, Mdims%npres
               CALL C_MULT2( CDP_tensor%val( :, 1+(ipres-1)*Mdims%n_in_pres : ipres*Mdims%n_in_pres, : ), deltap%val( IPRES, : ), &
                    Mdims%cv_nonods, Mdims%u_nonods, Mdims%ndim, Mdims%n_in_pres, Mmat%C( :, 1+(ipres-1)*Mdims%n_in_pres : ipres*Mdims%n_in_pres, : ), Mspars%C%ncol, Mspars%C%fin, Mspars%C%col )
            END DO
            call deallocate(deltaP)
            call halo_update(cdp_tensor)
            ! Correct velocity...
            ! DU = BLOCK_MAT * CDP

            ALLOCATE( DU_VEL( Mdims%ndim,  Mdims%nphase, Mdims%u_nonods )) ; DU_VEL = 0.
            CALL PHA_BLOCK_MAT_VEC2( DU_VEL, Mmat%PIVIT_MAT, CDP_tensor%val, Mdims%ndim, Mdims%nphase, &
            Mdims%totele, Mdims%u_nloc, ndgln%u )
            U_ALL2 % VAL = U_ALL2 % VAL + DU_VEL

            DEALLOCATE( DU_VEL )
            if ( after_adapt .and. cty_proj_after_adapt ) UOLD_ALL2 % VAL = U_ALL2 % VAL
            call halo_update(u_all2)
            call DEALLOCATE( CDP_tensor )
        END if
        ! Calculate control volume averaged pressure CV_P from fem pressure P
        CVP_ALL%VAL = 0.0
        IF(Mdims%npres>1.AND.PIPES_1D) THEN
           MASS_CV = 0.0
           IPRES = 1
           DO CV_NOD = 1, Mdims%cv_nonods
              if (node_owned(CVP_all,CV_NOD)) then
                 DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                    CVP_all%val( 1, IPRES, CV_NOD ) = CVP_all%val( 1, IPRES, CV_NOD ) + MASS_MN_PRES( COUNT ) * P_all%val( 1, IPRES, Mspars%CMC%col( COUNT ) )
                    MASS_CV( CV_NOD ) = MASS_CV( CV_NOD ) + MASS_MN_PRES( COUNT )
                 END DO
              else
                 Mass_CV(CV_NOD)=1.0
              end if
           END DO
           CVP_all%val(1,IPRES,:) = CVP_all%val(1,IPRES,:) / MASS_CV
           MASS_CV = 0.0
           IPRES = Mdims%npres
           DO CV_NOD = 1, Mdims%cv_nonods
              if (node_owned(CVP_all,CV_NOD)) then
                 DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                    CVP_all%val( 1, IPRES, CV_NOD ) = CVP_all%val( 1, IPRES, CV_NOD ) + MASS_CVFEM2PIPE_TRUE( COUNT ) * P_all%val( 1, IPRES, Mspars%CMC%col( COUNT ) )
                    MASS_CV( CV_NOD ) = MASS_CV( CV_NOD ) + MASS_CVFEM2PIPE_TRUE( COUNT )
                 END DO
                 MASS_CV( CV_NOD ) = max( 1.0e-15, MASS_CV( CV_NOD ) )
              else
                 Mass_CV(CV_NOD)=1.0
              end if
           END DO
           CVP_all%val(1,IPRES,:) = CVP_all%val(1,IPRES,:) / MASS_CV
        ELSE
           MASS_CV = 0.0
           DO CV_NOD = 1, Mdims%cv_nonods
              if (node_owned(CVP_all,CV_NOD)) then
                 DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                    CVP_all%val( 1, :, CV_NOD ) = CVP_all%val( 1, :, CV_NOD ) + MASS_MN_PRES( COUNT ) * P_all%val( 1, :, Mspars%CMC%col( COUNT ) )
                    MASS_CV( CV_NOD ) = MASS_CV( CV_NOD ) + MASS_MN_PRES( COUNT )
                 END DO
              else
                 Mass_CV(CV_NOD)=1.0
              end if
           END DO
           DO IPRES = 1, Mdims%npres
              CVP_all%val(1,IPRES,:) = CVP_all%val(1,IPRES,:) / MASS_CV
           END DO
        ENDIF
        call halo_update(CVP_all)
        DEALLOCATE( Mmat%CT )
        DEALLOCATE( DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B )
        DEALLOCATE( MASS_PIPE, MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE )
        DEALLOCATE( Mmat%U_RHS )
        DEALLOCATE( MASS_MN_PRES )
        call deallocate(Mmat%CT_RHS)
        if (.not.is_porous_media) then
            DEALLOCATE( Mmat%PIVIT_MAT )
            nullify(Mmat%PIVIT_MAT)
        end if
        ewrite(3,*) 'Leaving FORCE_BAL_CTY_ASSEM_SOLVE'
        return
    END SUBROUTINE FORCE_BAL_CTY_ASSEM_SOLVE


    SUBROUTINE CV_ASSEMB_FORCE_CTY( state, packed_state, &
        Mdims, CV_GIdims, FE_GIdims, CV_funs, FE_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd, &
        velocity, pressure, &
        X_ALL, velocity_absorption, U_SOURCE_ALL, U_SOURCE_CV_ALL, &
        U_ALL, UOLD_ALL, &
        P, CV_P, DEN_ALL, DENOLD_ALL, DERIV, &
        DT, &
        MASS_MN_PRES,&
        got_free_surf,  MASS_SUF, &
        SUF_SIG_DIAGTEN_BC, &
        V_SOURCE, VOLFRA_PORE, &
        MCY_RHS, DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, GAMMA_PRES_ABS, GAMMA_PRES_ABS_NANO, INV_B, MASS_PIPE, MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE, GLOBAL_SOLVE, &
        MCY, JUST_BL_DIAG_MAT, &
        UDEN_ALL, UDENOLD_ALL, UDIFFUSION_ALL, UDIFFUSION_VOL_ALL, THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL, &
        IGOT_THETA_FLUX, &
        THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, &
        RETRIEVE_SOLID_CTY, &
        IPLIKE_GRAD_SOU, &
        symmetric_P, boussinesq, IDs_ndgln, calculate_mass_delta)
        implicit none
        ! Form the global CTY and momentum eqns and combine to form one large matrix eqn.
        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_GI_dimensions), intent(in) :: CV_GIdims, FE_GIdims
        type(multi_shape_funs), intent(inout) :: CV_funs
        type(multi_shape_funs), intent(in) :: FE_funs
        type(multi_sparsities), intent(in) :: Mspars
        type(multi_ndgln), intent(in) :: ndgln
        type(multi_discretization_opts) :: Mdisopt
        type(multi_matrices), intent(inout) :: Mmat
        type (porous_adv_coefs), intent(inout) :: upwnd
        type( tensor_field ), intent(in) :: velocity
        type( tensor_field ), intent(in) :: pressure
        INTEGER, intent( in ) :: IGOT_THETA_FLUX, IPLIKE_GRAD_SOU
        LOGICAL, intent( in ) :: RETRIEVE_SOLID_CTY,got_free_surf,symmetric_P,boussinesq
        INTEGER, DIMENSION( : ), intent( in ) :: IDs_ndgln
        real, dimension(:,:), intent(in) :: X_ALL
        REAL, DIMENSION( :, :, : ), intent( in ) :: velocity_absorption
        type( multi_field ), intent( in ) :: U_SOURCE_ALL
        REAL, DIMENSION( :, :, : ), intent( in ) :: U_SOURCE_CV_ALL
        REAL, DIMENSION( :, :, : ), intent( in ) :: U_ALL, UOLD_ALL
        REAL, DIMENSION( :, :, : ), intent( in ) :: CV_P, P
        REAL, DIMENSION(  :, :  ), intent( in ), pointer :: DEN_ALL, DENOLD_ALL
        REAL, DIMENSION(  : , :  ), intent( in ) :: DERIV
        REAL, DIMENSION(  : ,  :   ), intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
        REAL, intent( in ) :: DT
        REAL, DIMENSION(  : , : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
        REAL, DIMENSION(  :, :  ), intent( in ) :: V_SOURCE
        !REAL, DIMENSION( :, :, : ), intent( in ) :: V_ABSORB
        REAL, DIMENSION( :, : ), intent( in ) :: VOLFRA_PORE
        REAL, DIMENSION( : ), intent( inout ) :: MCY_RHS
        REAL, DIMENSION( : ), intent( inout ) :: MASS_MN_PRES
        REAL, DIMENSION( : ), intent( inout ) :: MASS_SUF
        REAL, DIMENSION( :, : ), intent( inout ), allocatable :: DIAG_SCALE_PRES
        REAL, DIMENSION( :, :, : ), intent( inout ), allocatable :: DIAG_SCALE_PRES_COUP, GAMMA_PRES_ABS, GAMMA_PRES_ABS_NANO, INV_B
        REAL, DIMENSION( : ), intent( inout ) :: MASS_PIPE, MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE
        LOGICAL, intent( in ) :: GLOBAL_SOLVE
        REAL, DIMENSION( : ), intent( inout ) :: MCY
        REAL, DIMENSION( :, : ), intent( in ) :: UDEN_ALL, UDENOLD_ALL
        REAL, DIMENSION( :, :, :, : ), intent( in ) :: UDIFFUSION_ALL
        REAL, DIMENSION( :, :, :, : ), intent( inout ) :: THERM_U_DIFFUSION
        type( multi_field ), intent( in ) :: UDIFFUSION_VOL_ALL
        REAL, DIMENSION( :, : ), intent( inout ) :: THERM_U_DIFFUSION_VOL
        LOGICAL, intent( inout ) :: JUST_BL_DIAG_MAT
        ! Local variables
        REAL, PARAMETER :: v_beta = 1.0
! NEED TO CHANGE RETRIEVE_SOLID_CTY TO MAKE AN OPTION
        REAL, PARAMETER :: SECOND_THETA = 1.0
        LOGICAL, PARAMETER :: GETCV_DISC = .FALSE., GETCT= .TRUE., THERMAL= .FALSE.
        REAL, DIMENSION( : ), allocatable ::  dummy_transp
        REAL, DIMENSION( :,:,:,: ), allocatable :: TDIFFUSION
        REAL, DIMENSION( :, : ), allocatable :: THETA_GDIFF
        REAL, DIMENSION( : , : ), allocatable :: DEN_OR_ONE, DENOLD_OR_ONE
        REAL, DIMENSION( :, : ), allocatable :: MEAN_PORE_CV
        LOGICAL :: GET_THETA_FLUX
        INTEGER :: IGOT_T2, I, IGOT_THERM_VIS
        INTEGER :: ELE, U_ILOC, U_INOD, IPHASE, IDIM
        type(tensor_field), pointer :: tracer, density
        REAL, DIMENSION( : , :, : ), pointer :: V_ABSORB => null() ! this is PhaseVolumeFraction_AbsorptionTerm

        real, dimension(:,:) :: calculate_mass_delta
!        real, dimension(:) :: calculate_mass_internal_previous


        ewrite(3,*)'In CV_ASSEMB_FORCE_CTY'
        GET_THETA_FLUX = .FALSE.
        IGOT_T2 = 0
        ALLOCATE( THETA_GDIFF( Mdims%nphase * IGOT_T2, Mdims%cv_nonods * IGOT_T2 )) ; THETA_GDIFF = 0.
        ALLOCATE( TDIFFUSION( Mdims%mat_nonods, Mdims%ndim, Mdims%ndim, Mdims%nphase )) ; TDIFFUSION = 0.
        ALLOCATE( MEAN_PORE_CV( Mdims%npres, Mdims%cv_nonods )) ; MEAN_PORE_CV = 0.
        allocate( dummy_transp( Mdims%totele ) ) ; dummy_transp = 0.
        TDIFFUSION = 0.0
        IF( GLOBAL_SOLVE ) MCY = 0.0
        ! Obtain the momentum and Mmat%C matricies
        CALL ASSEMB_FORCE_CTY( state, packed_state, &
            Mdims, FE_GIdims, FE_funs, Mspars, ndgln, Mmat, &
            velocity, pressure, &
            X_ALL, velocity_absorption, U_SOURCE_ALL, U_SOURCE_CV_ALL, &
            U_ALL, UOLD_ALL, &
            U_ALL, UOLD_ALL, &    ! This is nu...
            UDEN_ALL, UDENOLD_ALL, DERIV, &
            DT, &
            JUST_BL_DIAG_MAT, &
            UDIFFUSION_ALL, UDIFFUSION_VOL_ALL, THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL, DEN_ALL, RETRIEVE_SOLID_CTY, &
            IPLIKE_GRAD_SOU, &
            P, GOT_FREE_SURF=got_free_surf, MASS_SUF=MASS_SUF, SYMMETRIC_P=symmetric_P)
        ! scale the momentum equations by the volume fraction / saturation for the matrix and rhs
        IF ( GLOBAL_SOLVE ) THEN
            ! put momentum and Mmat%C matrices into global matrix MCY...
            MCY_RHS = 0.0
            DO ELE = 1, Mdims%totele
                DO U_ILOC = 1, Mdims%u_nloc
                    U_INOD = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                    DO IPHASE = 1, Mdims%nphase
                        DO IDIM = 1, Mdims%ndim
                            I = U_INOD + (IDIM-1)*Mdims%u_nonods + (IPHASE-1)*Mdims%ndim*Mdims%u_nonods
                            MCY_RHS( I ) = Mmat%U_RHS( IDIM, IPHASE, U_INOD )
                        END DO
                    END DO
                END DO
            END DO
FLAbort('Global solve for pressure-mommentum is broken until nested matrices get impliented.')
!            CALL PUT_MOM_C_IN_GLOB_MAT( Mdims%nphase,Mdims%ndim, &
!            Mspars%DGM_PHA%ncol, Mmat%DGM_petsc, Mspars%DGM_PHA%fin, &
!            NLENMCY, Mspars%MCY%ncol, MCY, Mspars%MCY%fin, &
!            Mdims%u_nonods, Mspars%C%ncol, Mmat%C, Mspars%C%fin )
        END IF
        ALLOCATE( DEN_OR_ONE( Mdims%nphase, Mdims%cv_nonods )); DEN_OR_ONE = 1.
        ALLOCATE( DENOLD_OR_ONE( Mdims%nphase, Mdims%cv_nonods )); DENOLD_OR_ONE = 1.
        IF ( Mdisopt%volfra_use_theta_flux ) THEN ! We have already put density in theta...
           DEN_OR_ONE = 1.
           DENOLD_OR_ONE = 1.
        ELSE
           DEN_OR_ONE = DEN_ALL
           DENOLD_OR_ONE = DENOLD_ALL
        END IF
        if ( boussinesq ) then
           DEN_OR_ONE = 1.0
           DENOLD_OR_ONE = 1.0
        end if
        ! no q scheme
        IGOT_THERM_VIS = 0
        tracer=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
        density=>extract_tensor_field(packed_state,"PackedDensity")
        !For flooding ensure that the height is non-zero and positive
        if (is_flooding) density%val(1,1,:) = max(density%val(1,1,:),1e-5)
        call halo_update(density)
        call CV_ASSEMB( state, packed_state, &
            Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd, &
            tracer, velocity, density, &
            DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, GAMMA_PRES_ABS, GAMMA_PRES_ABS_NANO, &
            INV_B, MASS_PIPE, MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE, &
            DEN_OR_ONE, DENOLD_OR_ONE, &
            TDIFFUSION, IGOT_THERM_VIS, THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL, &
            Mdisopt%v_disopt, Mdisopt%v_dg_vel_int_opt, DT, Mdisopt%v_theta, SECOND_THETA, v_beta, &
            SUF_SIG_DIAGTEN_BC, &
            DERIV, CV_P, &
            V_SOURCE, V_ABSORB, VOLFRA_PORE, &
            GETCV_DISC, GETCT, &
            IGOT_T2, IGOT_THETA_FLUX, GET_THETA_FLUX, Mdisopt%volfra_use_theta_flux, &
            THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
            MEAN_PORE_CV, &
            MASS_MN_PRES, THERMAL,  RETRIEVE_SOLID_CTY,&
            got_free_surf,  MASS_SUF, &
            dummy_transp, IDs_ndgln, &
            calculate_mass_delta = calculate_mass_delta)

        !For flooding ensure that the height is non-zero and positive
        if (is_flooding) density%val(1,1,:) = max(density%val(1,1,:),1e-5)

        ewrite(3,*)'Back from cv_assemb'
        IF ( GLOBAL_SOLVE ) THEN
            ! Put Mmat%CT into global matrix MCY...
            MCY_RHS( Mdims%u_nonods * Mdims%ndim * Mdims%nphase + 1 : Mdims%u_nonods * Mdims%ndim * Mdims%nphase + Mdims%cv_nonods ) = &
            Mmat%CT_RHS%val( 1, 1 : Mdims%cv_nonods )
            CALL PUT_CT_IN_GLOB_MAT( Mmat, MCY )
        END IF
        deallocate( DEN_OR_ONE, DENOLD_OR_ONE )
        DEALLOCATE( THETA_GDIFF )
        DEALLOCATE( TDIFFUSION )
        DEALLOCATE( MEAN_PORE_CV )
        ewrite(3,*) 'Leaving CV_ASSEMB_FORCE_CTY'


        contains
            SUBROUTINE PUT_CT_IN_GLOB_MAT( Mmat, MCY)
                implicit none
                ! Put Mmat%CT into global matrix MCY
                type (multi_matrices), intent(inout) :: Mmat
                REAL, DIMENSION( : ), intent( inout ) :: MCY

                ! Local variables...
                INTEGER CV_NOD, IWID, COUNT, IPHASE, COUNT_MCY1, &
                    COUNT_MCY, COUNT_CMC, IDIM, I
                ewrite(3,*) 'In PUT_CT_IN_GLOB_MAT'
                Loop_CVNOD: DO CV_NOD = 1, Mdims%cv_nonods
                    IWID = Mspars%CT%fin( CV_NOD + 1 ) - Mspars%CT%fin( CV_NOD )
                    Loop_COUNT: DO COUNT = Mspars%CT%fin( CV_NOD ), Mspars%CT%fin( CV_NOD + 1 ) - 1
                        Loop_PHASE: DO IPHASE = 1, Mdims%nphase
                            Loop_DIM: DO IDIM = 1, Mdims%ndim
                                COUNT_MCY1 = Mspars%MCY%fin( Mdims%u_nonods * Mdims%nphase * Mdims%ndim + CV_NOD ) - 1 + (COUNT - Mspars%CT%fin( CV_NOD ) +1) &
                                    + ( IPHASE - 1 ) * IWID * Mdims%ndim &
                                    + IWID*(IDIM-1)
                                MCY( COUNT_MCY1 ) = Mmat%CT( IDIM, IPHASE, COUNT )
                            END DO Loop_DIM
                        END DO Loop_PHASE
                    END DO Loop_COUNT
                END DO Loop_CVNOD
                DO CV_NOD = 1, Mdims%cv_nonods
                    IWID = Mspars%CMC%fin( CV_NOD + 1 )- Mspars%CMC%fin( CV_NOD )
                    DO I = 1, IWID
                        COUNT_CMC = Mspars%CMC%fin( CV_NOD + 1) - I
                        COUNT_MCY = Mspars%MCY%fin( Mdims%ndim * Mdims%nphase * Mdims%u_nonods + CV_NOD + 1 ) - I
                        MCY( COUNT_MCY ) = DIAG_SCALE_PRES( 1, CV_NOD ) * MASS_MN_PRES( COUNT_CMC )
                    END DO
                END DO
                ewrite(3,*) 'Leaving PUT_CT_IN_GLOB_MAT'
                RETURN
            END SUBROUTINE PUT_CT_IN_GLOB_MAT
    END SUBROUTINE CV_ASSEMB_FORCE_CTY




    SUBROUTINE ASSEMB_FORCE_CTY( state, packed_state, &
        Mdims, FE_GIdims, FE_funs, Mspars, ndgln, Mmat, &
        velocity, pressure, &
        X_ALL, U_ABSORB, U_SOURCE, U_SOURCE_CV, &
        U_ALL, UOLD_ALL, &
        NU_ALL, NUOLD_ALL, &
        UDEN, UDENOLD, DERIV, &
        DT, JUST_BL_DIAG_MAT,  &
        UDIFFUSION, UDIFFUSION_VOL, THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL, DEN_ALL, RETRIEVE_SOLID_CTY, &
        IPLIKE_GRAD_SOU, &
        P,&
        got_free_surf, mass_suf, symmetric_P )
        implicit none
        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_GI_dimensions), intent(in) :: FE_GIdims
        type(multi_shape_funs), intent(in) :: FE_funs
        type (multi_sparsities), intent(in) :: Mspars
        type(multi_ndgln), intent(in) :: ndgln
        type (multi_matrices), intent(inout) :: Mmat
        type( tensor_field ), intent( in ) :: velocity
        type( tensor_field ), intent( in ) :: pressure
! If IGOT_VOL_X_PRESSURE=1 then have a voln fraction in the pressure term and multiply density by DevFuns%VOLUME fraction...
        INTEGER, PARAMETER :: IGOT_VOL_X_PRESSURE = 0
        INTEGER, intent( in ) :: IPLIKE_GRAD_SOU
        REAL, DIMENSION( :, : ), intent( in ) :: X_ALL
        REAL, DIMENSION( :, :, : ), intent( in ) :: U_ABSORB
        type( multi_field ), intent( in ) :: U_SOURCE
        REAL, DIMENSION( :, :, : ), intent( in ) :: U_SOURCE_CV
        REAL, DIMENSION ( :, :, : ), intent( in ) :: U_ALL, UOLD_ALL, NU_ALL, NUOLD_ALL
        REAL, DIMENSION( :, : ), intent( in ) :: UDEN, UDENOLD, DERIV
        REAL, intent( in ) :: DT
        REAL, DIMENSION( :, :, :, : ), intent( in ) :: UDIFFUSION
        type( multi_field ), intent( in ) :: UDIFFUSION_VOL
        REAL, DIMENSION( :, :, :, : ), intent( inout ) :: THERM_U_DIFFUSION
        REAL, DIMENSION( :, : ), intent( inout ) :: THERM_U_DIFFUSION_VOL
        LOGICAL, intent( inout ) :: JUST_BL_DIAG_MAT
        REAL, DIMENSION( :, :, : ), intent( in ) :: P
        REAL, DIMENSION(  :, :  ), intent( in ) :: DEN_ALL
        LOGICAL, intent( in ) :: RETRIEVE_SOLID_CTY, got_free_surf, symmetric_P
        REAL, DIMENSION( : ), intent( inout ) :: MASS_SUF
        ! Local Variables
        ! This is for decifering WIC_U_BC & WIC_P_BC
        LOGICAL, PARAMETER :: VOL_ELE_INT_PRES = .TRUE., STAB_VISC_WITH_ABS=.FALSE.
        LOGICAL :: STRESS_FORM, STRESS_FORM_STAB, THERMAL_STAB_VISC, THERMAL_LES_VISC, THERMAL_FLUID_VISC, Q_SCHEME
        ! if STAB_VISC_WITH_ABS then stabilize (in the projection mehtod) the viscosity using absorption.
        REAL, PARAMETER :: WITH_NONLIN = 1.0, TOLER = 1.E-10
        ! NON_LIN_DGFLUX = .TRUE. non-linear DG flux for momentum - if we have an oscillation use upwinding else use central scheme.
        ! UPWIND_DGFLUX=.TRUE. Upwind DG flux.. Else use central scheme. if NON_LIN_DGFLUX = .TRUE. then this option is ignored.
        LOGICAL :: NON_LIN_DGFLUX, UPWIND_DGFLUX
        ! Storage for pointers to the other side of the element.
        ! Switched off for now until this is hooked up.
        LOGICAL, PARAMETER :: STORED_OTHER_SIDE = .FALSE.
        INTEGER, PARAMETER :: ISTORED_OTHER_SIDE = 0
        ! This is for rapid access to the Mmat%C matrix...
        LOGICAL, PARAMETER :: STORED_AC_SPAR_PT=.FALSE.
        INTEGER, PARAMETER :: IDO_STORE_AC_SPAR_PT=0
        ! re-calculate Mmat%C matrix...
        LOGICAL :: got_c_matrix
        INTEGER, DIMENSION( :, : ), allocatable ::  FACE_ELE
        INTEGER, DIMENSION( : ), allocatable :: CV_SLOC2LOC, U_SLOC2LOC, &
        U_ILOC_OTHER_SIDE, U_OTHER_LOC, MAT_OTHER_LOC
        REAL, DIMENSION( : ),    ALLOCATABLE ::  &
        SDETWE, VLN,VLN_OLD, &
        MASS_ELE
        REAL, DIMENSION( :, : ),    ALLOCATABLE :: XL_ALL, XL2_ALL, XSL_ALL, SNORMXN_ALL
        REAL, DIMENSION ( : , :,  : ), allocatable :: GRAD_SOU_GI_NMX, GRAD_SOU_GI, GRAD_SOU2_GI_NMX, GRAD_SOU2_GI
        REAL, DIMENSION ( : , :,  : ), allocatable :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD
        REAL, DIMENSION( : ),    ALLOCATABLE :: NORMX_ALL
        REAL, DIMENSION ( : , :,  : ), allocatable :: SIGMAGI, SIGMAGI_STAB, SIGMAGI_STAB_SOLID_RHS, &
        DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX, FTHETA, SNDOTQ_IN, SNDOTQ_OUT, &
        SNDOTQOLD_IN, SNDOTQOLD_OUT, UD, UDOLD, UD_ND, UDOLD_ND
        REAL, DIMENSION ( : , : ), allocatable :: &
        DENGI, DENGIOLD, &
        SNDOTQ, SNDOTQOLD, SINCOME, SINCOMEOLD, SDEN, SDENOLD, &
        SDEN_KEEP, SDENOLD_KEEP, SDEN2_KEEP, SDENOLD2_KEEP, &
        SNDOTQ_KEEP, SNDOTQ2_KEEP, SNDOTQOLD_KEEP, SNDOTQOLD2_KEEP, &
        N_DOT_DU, N_DOT_DU2, N_DOT_DUOLD, N_DOT_DUOLD2, RHS_U_CV, RHS_U_CV_OLD, UDEN_VFILT, UDENOLD_VFILT
        REAL, DIMENSION ( : , :, : ), allocatable :: SUD_ALL, SUDOLD_ALL, SUD2_ALL, SUDOLD2_ALL, SUD_ALL_KEEP, &
        SUDOLD_ALL_KEEP, SUD2_ALL_KEEP, SUDOLD2_ALL_KEEP
        REAL, DIMENSION ( : ), allocatable :: vel_dot, vel_dot2, velold_dot, velold_dot2, grad_fact
        ! Nonlinear Petrov-Galerkin stuff...
        REAL, DIMENSION ( : , : ), allocatable ::LOC_MASS_INV, LOC_MASS, P_DX        
        REAL, DIMENSION ( : ), allocatable :: VLK_UVW, U_R2_COEF, U_GRAD_N_MAX2
        REAL, DIMENSION ( :, :, : ), allocatable :: &
        MAT_ELE, DIFFGI_U, RHS_DIFF_U, DIFF_VEC_U, SOUGI_X, RESID_U, U_DT, &
        DIF_STAB_U, U_GRAD_NORM2, U_GRAD_NORM, A_DOT_U, STAR_U_COEF, P_STAR_U
        REAL, DIMENSION ( :, :, :, :, : ), allocatable :: UDIFF_SUF_STAB
        !###Shape function calculation###
        type(multi_dev_shape_funs) :: Devfuns

! Local variables...
!            INTEGER, PARAMETER :: LES_DISOPT=0
            INTEGER :: LES_DISOPT
! LES_DISOPT is LES option e.g. =0 No LES
!                               =1 Anisotropic element length scale
!                               =2 Take the average length scale h
!                               =3 Take the min length scale h
!                               =4 Take the max length scale h
!                               =5 for stress form take length scale across gradient direction of the element.
!                               =6 same as 5 plus original q-scheme that is multiply by -min(0.0,divq) (5 can be switched off by using a zero coefficient for the LES)
!                               =7 same as 5 plus original q-scheme but using abs(divu) (5 can be switched off by using a zero coefficient for the LES)
!            REAL, PARAMETER :: LES_THETA=1.0
            REAL :: LES_THETA, LES_CS
! LES_THETA =1 is backward Euler for the LES viscocity.
! COEFF_SOLID_FLUID is the coeffficient that determins the magnitude of the relaxation to the solid vel...
! min_den_for_solid_fluid is the minimum density that is used in the solid-fluid coupling term.
!!-PY: COEFF_SOLID_FLUID_stab=1.0 means there is relaxation in the shell by using the immersed-shell method.
!!-PY: COEFF_SOLID_FLUID_relax=0.0 means there is no relaxtion inside the solid. If switch to COEFF_SOLID_FLUID_relax=1.0, that means there is relaxtion inside the solid.
            REAL, PARAMETER :: min_den_for_solid_fluid = 1.0, COEFF_SOLID_FLUID_stab=1.0, COEFF_SOLID_FLUID_relax=0.0
! include_viscous_solid_fluid_drag_force switches on the solid-fluid coupling viscocity boundary conditions...
!            LOGICAL, PARAMETER :: include_viscous_solid_fluid_drag_force = .FALSE.
            LOGICAL :: include_viscous_solid_fluid_drag_force
! If PIVIT_ON_VISC then place the block diaongal viscocity into the pivit matrix used in the projection method...
! PIVIT_ON_VISC is the only thing that could make highly viscouse flows stabe when using projection methods...
            LOGICAL :: PIVIT_ON_VISC !THE VALUE IS SET UP BEFORE THE get_option
! GOT_VIRTUAL_MASS ! do we have virtual mass terms for multi-phase flows...
            LOGICAL :: GOT_VIRTUAL_MASS = .false.
! If FEM_DEN then use an FEM representation of density - only used within an element (default is FEM for between elements and on boundary).
            LOGICAL, PARAMETER :: FEM_DEN = .false.
            real :: w
            real, parameter :: wv=1.0, ws=1.0 ! DevFuns%VOLUME off-diagonal and surface weights, respectively
! LINEAR_HIGHORDER_DIFFUSION is the switch for the high-order linear scheme...
!            LOGICAL, PARAMETER ::  LINEAR_HIGHORDER_DIFFUSION = .TRUE.
            LOGICAL ::  LINEAR_HIGHORDER_DIFFUSION
        !
        ! Variables used to reduce indirect addressing...
        !INTEGER, DIMENSION ( :, :, : ), allocatable :: WIC_U_BC_ALL
        REAL, DIMENSION ( :, :, : ), allocatable :: LOC_U_RHS
        REAL, DIMENSION ( :, :, : ), allocatable :: LOC_U, LOC_UOLD, LOC_US, LOC_U_ABS_STAB_SOLID_RHS
        REAL, DIMENSION ( :, :, : ), allocatable :: LOC_NU, LOC_NUOLD
        REAL, DIMENSION ( :, :, : ), allocatable :: LOC_U_ABSORB, LOC_U_ABS_STAB
        REAL, DIMENSION ( :, :, :, : ), allocatable :: LOC_UDIFFUSION, U_DX_ALL, UOLD_DX_ALL, DIFF_FOR_BETWEEN_U
        REAL, DIMENSION ( :, : ), allocatable :: LOC_UDIFFUSION_VOL
        !  tHE VIRTUAL MASS COEFFICIENT :VIRTUAL_MASS
        ! VIRTUAL_MASS_ADV_CURR =1. if advect with the velocity of current phase else =0 then use the coln phase...
        REAL, DIMENSION ( :, :, : ), allocatable :: VIRTUAL_MASS, VIRTUAL_MASS_OLD
        REAL, DIMENSION ( :, : ), allocatable :: VIRTUAL_MASS_ADV_CUR
        REAL, DIMENSION ( :, :, : ), allocatable :: LOC_VIRTUAL_MASS, LOC_VIRTUAL_MASS_OLD
        REAL, DIMENSION ( :, :, : ), allocatable :: VIRTUAL_MASS_GI, VIRTUAL_MASS_OLD_GI
        REAL, DIMENSION ( :, : ), allocatable :: VLN_CVM, VLN_OLD_CVM
        REAL, DIMENSION ( :, :, :, : ), allocatable :: CVM_SNDOTQ_IN, CVM_SNDOTQ_OUT, CVM_SNDOTQOLD_IN, CVM_SNDOTQOLD_OUT
        REAL, DIMENSION ( : ), allocatable :: CVM_NN_SNDOTQ_IN, CVM_NN_SNDOTQ_OUT, CVM_NN_SNDOTQOLD_IN, CVM_NN_SNDOTQOLD_OUT
        REAL, DIMENSION ( :, :, : ), allocatable :: SVIRTUAL_MASS_GI, SVIRTUAL_MASS_OLD_GI
        REAL, DIMENSION ( :, :, : ), allocatable :: SVIRTUAL_MASS_GI_KEEP, SVIRTUAL_MASS_GI2_KEEP
        REAL, DIMENSION ( :, :, : ), allocatable :: SVIRTUAL_MASS_OLD_GI_KEEP, SVIRTUAL_MASS_OLD_GI2_KEEP
        REAL, DIMENSION ( :, :, : ), allocatable :: SLOC_VIRTUAL_MASS, SLOC2_VIRTUAL_MASS
        REAL, DIMENSION ( :, :, : ), allocatable :: SLOC_VIRTUAL_MASS_OLD, SLOC2_VIRTUAL_MASS_OLD
        REAL :: WITH_NONLIN_CVM, CVM_BETA
        ! For q-scheme etc...
        REAL, DIMENSION ( :, : ), allocatable :: MAT_ELE_CV_LOC, INV_MAT_ELE_CV_LOC
        REAL, DIMENSION ( :, :, : ), allocatable :: DIFF_FOR_BETWEEN_CV
        REAL, DIMENSION ( :, :, : ), allocatable ::    DIFFCV
        REAL, DIMENSION ( :, :, :, : ), allocatable :: DIFFCV_TEN
        REAL, DIMENSION ( :, :, :, :, : ), allocatable :: DIFFCV_TEN_ELE
        REAL, DIMENSION ( : ), allocatable :: RCOUNT_NODS
        !REAL, DIMENSION ( :, :, :, : ), allocatable :: SUF_U_BC_ALL, SUF_MOM_BC_ALL, SUF_NU_BC_ALL, SUF_ROB1_UBC_ALL, SUF_ROB2_UBC_ALL, TEN_XX
        REAL, DIMENSION ( :, :, :, : ), allocatable :: TEN_XX
        REAL, DIMENSION ( :, : ), allocatable :: TEN_VOL
        !REAL, DIMENSION ( :, :, : ), allocatable :: SUF_P_BC_ALL
        REAL, DIMENSION ( :, : ), allocatable :: LOC_UDEN,  LOC_UDENOLD
        REAL, DIMENSION ( : ), allocatable :: LOC_P
        REAL, DIMENSION ( :, :, : ), allocatable :: LOC_PLIKE_GRAD_SOU_COEF, LOC_PLIKE_GRAD_SOU_GRAD
        REAL, DIMENSION ( :, :, : ), allocatable :: LOC_U_SOURCE, LOC_U_SOURCE_CV
        REAL, DIMENSION ( :, :, :,   :, :, :,   : ), allocatable :: DIAG_BIGM_CON, BIGM_CON
        ! memory for fast retreval of surface info...
        INTEGER, DIMENSION ( :, :, : ), allocatable :: STORED_U_ILOC_OTHER_SIDE, STORED_U_OTHER_LOC, STORED_MAT_OTHER_LOC
        INTEGER, DIMENSION ( :, :, : ), allocatable :: POSINMAT_C_STORE
        INTEGER, DIMENSION ( :, :, :, : ), allocatable :: POSINMAT_C_STORE_SUF_DG
        ! To memory access very local...
        REAL, DIMENSION ( :, :, : ), allocatable :: SLOC_U, SLOC_UOLD, SLOC2_U, SLOC2_UOLD
        REAL, DIMENSION ( :, :, : ), allocatable :: SLOC_NU, SLOC_NUOLD, SLOC2_NU, SLOC2_NUOLD
        REAL, DIMENSION ( :, :, :, : ), allocatable :: SLOC_DUX_ELE_ALL, SLOC2_DUX_ELE_ALL, SLOC_DUOLDX_ELE_ALL, SLOC2_DUOLDX_ELE_ALL
        REAL, DIMENSION ( :, : ), allocatable :: SLOC_UDEN, SLOC2_UDEN, SLOC_UDENOLD, SLOC2_UDENOLD
        REAL, DIMENSION ( :, :, :, : ), allocatable :: SLOC_UDIFFUSION, SLOC2_UDIFFUSION
        REAL, DIMENSION ( :, : ), allocatable :: SLOC_UDIFFUSION_VOL, SLOC2_UDIFFUSION_VOL
        REAL, DIMENSION ( :, :, : ), allocatable :: SLOC_DIFF_FOR_BETWEEN_U, SLOC2_DIFF_FOR_BETWEEN_U
        REAL, DIMENSION ( :, :, : ), allocatable :: U_NODI_SGI_IPHASE_ALL, U_NODJ_SGI_IPHASE_ALL, UOLD_NODI_SGI_IPHASE_ALL, UOLD_NODJ_SGI_IPHASE_ALL
        REAL, DIMENSION ( :, :, : ), allocatable :: CENT_RELAX, CENT_RELAX_OLD
        ! For derivatives...
        REAL, DIMENSION ( : ), allocatable :: NMX_ALL, VNMX_ALL
        LOGICAL :: GOT_DIFFUS, GOT_UDEN, DISC_PRES, QUAD_OVER_WHOLE_ELE
        INTEGER :: IPHASE, KPHASE, ELE, GI, IU_NOD, JCV_NOD, &
            COUNT, COUNT2, IPHA_IDIM, JPHA_JDIM, MAT_NOD, SGI, SELE, &
            U_SILOC, P_SJLOC, &
            IFACE, U_ILOC, U_JLOC, I, J, MAT_ILOC, &
            IDIM, P_ILOC, P_JLOC, ELE2, ELE3, SELE2, &
            JU_NOD_DIM_PHA, &
            U_ILOC2, U_INOD, U_INOD2, U_JLOC2, &
            IU_NOD_DIM_PHA, X_INOD,  &
            U_SJLOC, MAT_ILOC2, MAT_INOD, MAT_INOD2, MAT_SILOC, &
            CV_ILOC, CV_JLOC, CV_NOD, P_JLOC2, P_JNOD, P_JNOD2, &
            CV_SILOC, JDIM, JPHASE, &
            cv_inod, COUNT_ELE, CV_ILOC2, CV_INOD2, IDIMSF,JDIMSF, &
            IPRES, ICOMP
        REAL    :: NN, NM, SAREA, R, RR
        REAL    :: HDC, VLM, VLM_NEW,VLM_OLD, NN_SNDOTQ_IN,NN_SNDOTQ_OUT, &
            NN_SNDOTQOLD_IN,NN_SNDOTQOLD_OUT, RNN, c1(Mdims%ndim), c2(Mdims%ndim)
        REAL    :: MASSE, MASSE2
        ! Nonlinear Petrov-Galerkin stuff...
        INTEGER :: RESID_BASED_STAB_DIF
        REAL :: U_NONLIN_SHOCK_COEF,RNO_P_IN_A_DOT
        REAL :: JTT_INV
        REAL :: VLKNN, zero_or_two_thirds
        INTEGER :: P_INOD, IDIM_VEL
        logical :: mom_conserv, lump_mass, lump_mass2, lump_absorption, BETWEEN_ELE_STAB
        real :: beta
        INTEGER :: FILT_DEN, J2, JU2_NOD_DIM_PHA
        LOGICAL :: SIMPLE_DIFF_CALC
        REAL :: DIFF_MIN_FRAC, DIFF_MAX_FRAC
        REAL :: AVE_SNORMXN_ALL(Mdims%ndim)
        !Variables to improve Mmat%PIVIT_MAT creation speed
        REAL, DIMENSION ( :, :, :, : ), allocatable :: NN_SIGMAGI_ELE, NN_SIGMAGI_STAB_ELE, NN_SIGMAGI_STAB_SOLID_RHS_ELE, NN_MASS_ELE, NN_MASSOLD_ELE
        REAL, DIMENSION ( :, :, :, :, : ), allocatable :: STRESS_IJ_ELE, DUX_ELE_ALL, DUOLDX_ELE_ALL
        REAL, DIMENSION ( :, :, : ), allocatable :: VLK_ELE
        REAL, DIMENSION ( :, :, :, : ), allocatable :: UDIFFUSION_ALL, LES_UDIFFUSION
        REAL, DIMENSION ( :, : ), allocatable :: UDIFFUSION_VOL_ALL, LES_UDIFFUSION_VOL
        logical :: Porous_media_PIVIT_not_stored_yet
        !Variables to account for the lumping homogeneization
        logical :: homogenize_mass_matrix
        real :: lump_weight
        logical :: use_simple_lumped_homogenenous_mass_matrix
        ! memory for linear high order viscocity calculation...
        REAL, DIMENSION ( :, :, :, :, : ), allocatable :: STRESS_IJ_ELE_EXT
        REAL, DIMENSION ( :, :, : ), allocatable :: S_INV_NNX_MAT12
        REAL, DIMENSION ( :, :, :, : ), allocatable :: NNX_MAT_ELE
        REAL, DIMENSION ( :, :, : ), allocatable :: NN_MAT_ELE
        ! solid fluid coupling visc. bc contribution...
        REAL, DIMENSION ( :, :, :, : ), allocatable :: ABS_SOLID_FLUID_COUP
        REAL, DIMENSION ( :, :, : ), allocatable :: FOURCE_SOLID_FLUID_COUP
        type( vector_field ), pointer :: f_x
        type( tensor_field ), pointer :: a_xx
        ! for the option where we divid by voln fraction...
        REAL, DIMENSION ( :, : ), allocatable :: VOL_FRA_GI, CV_DENGI
        REAL, DIMENSION ( :, :, : ), allocatable :: VOL_FRA_GI_DX_ALL
        REAL, DIMENSION ( :, : ), allocatable :: SLOC_VOL_FRA, SLOC2_VOL_FRA,  SVOL_FRA, SVOL_FRA2, VOL_FRA_NMX_ALL
        ! revere ordering of shape functions used to get optimized code...
        REAL, DIMENSION ( :, :, : ), allocatable ::  CVFENX_ALL_REVERSED, UFENX_ALL_REVERSED
        REAL, DIMENSION ( :, : ), allocatable ::  UFEN_REVERSED, CVN_REVERSED, CVFEN_REVERSED
        REAL, DIMENSION ( :, : ), allocatable :: SBCVFEN_REVERSED, SBUFEN_REVERSED
        !Capillary pressure variables
        logical :: capillary_pressure_activated, Diffusive_cap_only
        !! femdem
        type( vector_field ), pointer :: delta_u_all, us_all
        type( scalar_field ), pointer :: sf, sfield
        real, dimension( : ), allocatable :: vol_s_gi
        !! Boundary_conditions
        INTEGER, DIMENSION ( Mdims%ndim , Mdims%nphase , surface_element_count(velocity) )  :: WIC_U_BC_ALL, WIC_U_BC_ALL_VISC, &
            WIC_U_BC_ALL_ADV, WIC_MOMU_BC_ALL, WIC_NU_BC_ALL
        INTEGER, DIMENSION ( 1, Mdims%npres, surface_element_count(pressure) ) :: WIC_P_BC_ALL
        REAL, DIMENSION ( :, :, : ), pointer  :: SUF_U_BC_ALL, SUF_U_BC_ALL_VISC, SUF_MOMU_BC_ALL, SUF_NU_BC_ALL
        REAL, DIMENSION ( :, :, : ), pointer :: SUF_U_ROB1_BC_ALL, SUF_U_ROB2_BC_ALL
        REAL, DIMENSION ( :, :, : ), pointer :: SUF_P_BC_ALL
        type(tensor_field) :: pressure_BCs
        type(tensor_field) :: velocity_BCs, velocity_BCs_visc, velocity_BCs_adv, velocity_BCs_robin2
        type(tensor_field) :: momentum_BCs
        INTEGER, DIMENSION( 4 ), PARAMETER :: ELEMENT_CORNERS=(/1,3,6,10/)
        !!$ Variables used in the diffusion-like term: capilarity and surface tension:
        type( tensor_field ), pointer :: PLIKE_GRAD_SOU_COEF_F, PLIKE_GRAD_SOU_GRAD_F, a_s
        type( scalar_field ), pointer :: P_s_hat
        ! GRAVTY is used in the free surface method only...
        REAL :: GRAVTY
        REAL :: MM_GRAVTY
        integer :: count_suf, count_p, P_SILOC, stat
        type (vector_field), pointer :: position
        integer, dimension(:), pointer :: neighbours
        integer :: nb
        logical :: skip, FEM_BUOYANCY
        !variables for linear velocity relaxation
        type(tensor_field), pointer :: fem_vol_frac_f
        real, dimension( :, : ), pointer :: fem_vol_frac
        ! If =0 then false, if =1 then true.
        ! If =1 then modify the stress term to take into
        ! account dividing through by DevFuns%VOLUME fraction.
        integer, parameter :: IDIVID_BY_VOL_FRAC = 0

        ! gradU
        logical :: get_gradU
        type(tensor_field), pointer :: gradU

        fem_vol_frac_f => extract_tensor_field( packed_state, "PackedFEPhaseVolumeFraction" )
        fem_vol_frac => fem_vol_frac_f%val( 1, :, : )
        ! open the boiling test for two phases-gas and liquid
        if (have_option('/boiling')) then
            GOT_VIRTUAL_MASS=.true.
        end if
        call get_option( "/physical_parameters/gravity/magnitude", gravty, stat )
        position=>extract_vector_field(packed_state,"PressureCoordinate")
        !Check capillary options
        capillary_pressure_activated = have_option_for_any_phase('/multiphase_properties/capillary_pressure', Mdims%nphase)
        Diffusive_cap_only = have_option_for_any_phase('/multiphase_properties/capillary_pressure/Diffusive_cap_only', Mdims%nphase)
        !We set the value of logicals
        PIVIT_ON_VISC = .false.
        call get_option( "/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/viscosity_scheme/zero_or_two_thirds", zero_or_two_thirds, default=2./3. )
        ! Stress form for the fluid viscocity
        STRESS_FORM = have_option( '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/viscosity_scheme/stress_form' )
        ! Stress form for the Petrov-Galerkin viscocity
        STRESS_FORM_STAB = have_option( '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/stabilisation/stress_form' )
        ! Use Q-scheme ...
        Q_SCHEME = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/spatial_discretisation/control_volumes/q_scheme' )
        ! Put the fluid viscocity in the thermal energy eqn...
        THERMAL_FLUID_VISC = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/spatial_discretisation/control_volumes/q_scheme/include_fluid_viscosity' )
        ! Put the Petrov-Galerkin viscocity in the thermal energy eqn...
        THERMAL_STAB_VISC = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/spatial_discretisation/control_volumes/q_scheme/include_stabilisation_viscosity' )
        ! Put the LES viscocity in the thermal energy eqn...
        THERMAL_LES_VISC = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/spatial_discretisation/control_volumes/q_scheme/include_les_viscosity' )
        ! Put the LES theta value for time stepping (LES_THETA=1 is default)...
        call get_option( '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/les_model/les_theta', LES_THETA, default=1.0 )
        ! Put the LES viscocity in the thermal energy eqn...
        call get_option( '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/les_model/model' , les_disopt, default=0 )
        call get_option( '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/les_model/smagorinsky_coefficient', les_Cs , default=0.1 )
        !
        ! bc for solid-fluid viscous drag coupling...
        include_viscous_solid_fluid_drag_force = have_option( '/blasting/include_viscous_drag_force' )
        !
        !
        ! DIFF_MIN_FRAC is the fraction of the standard diffusion coefficient to use
        ! in the non-linear diffusion scheme. DIFF_MAX_FRAC is the maximum fraction.
        ! If SIMPLE_DIFF_CALC then use a simple and fast diffusion calculation.
        !         DIFF_MIN_FRAC = 0.2
        !         DIFF_MAX_FRAC = 100.0
        !         SIMPLE_DIFF_CALC = .FALSE. ! Need switches for this
        SIMPLE_DIFF_CALC = have_option( '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/viscosity_scheme/linear_scheme' )
        LINEAR_HIGHORDER_DIFFUSION=have_option( '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/viscosity_scheme/linear_scheme/high_order' )
        if (LINEAR_HIGHORDER_DIFFUSION) SIMPLE_DIFF_CALC=.false.
        if ( .not. SIMPLE_DIFF_CALC ) then
            call get_option( '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/viscosity_scheme/nonlinear_scheme/beta_viscosity_min', DIFF_MIN_FRAC, default=0.2 )
            call get_option( '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/viscosity_scheme/nonlinear_scheme/beta_viscosity_max', DIFF_MAX_FRAC, default=100. )
        end if
        !If we have calculated already the Mmat%PIVIT_MAT and stored then we don't need to calculate it again
        Porous_media_PIVIT_not_stored_yet = (.not.Mmat%Stored .or. (.not.is_porous_media .or. Mdims%npres > 1))

        !Decide whether to use the simple mass matrix (Volume/U_NLOC) or the one using shape functions
        use_simple_lumped_homogenenous_mass_matrix = (Porous_media_PIVIT_not_stored_yet .and.is_porous_media&
                 .and. (Mmat%CV_pressure .or. have_option('/numerical_methods/simple_mass_matrix')))

        !If we do not have an index where we have stored Mmat%C, then we need to calculate it
        got_c_matrix  = Mmat%stored .or. Mmat%CV_pressure
        ewrite(3,*) 'In ASSEMB_FORCE_CTY'
        !! get boundary information
        call get_entire_boundary_condition(pressure,&
            ['weakdirichlet','freesurface  '],&
            pressure_BCs,WIC_P_BC_ALL)
        call get_entire_boundary_condition(velocity,&
            ['weakdirichlet','robin        '],&
            velocity_BCs,WIC_U_BC_ALL,boundary_second_value=velocity_BCs_robin2)
        call get_entire_boundary_condition(velocity,&
            ['weakdirichlet_viscosity'],&
            velocity_BCs_VISC,WIC_U_BC_ALL_VISC)
        call get_entire_boundary_condition(velocity,&
            ['weakdirichlet_advection'],&
            velocity_BCs_ADV,WIC_U_BC_ALL_ADV)
        call get_entire_boundary_condition(velocity,&
            ['momentum     ','momentuminout'],&
            momentum_BCs,WIC_MOMU_BC_ALL)
        WIC_NU_BC_ALL=WIC_U_BC_ALL
        SUF_P_BC_ALL=>pressure_BCs%val
        SUF_U_BC_ALL=>velocity_BCs%val
        SUF_U_BC_ALL_VISC=>velocity_BCs_VISC%val
        SUF_NU_BC_ALL=>velocity_BCs%val
        SUF_MOMU_BC_ALL=>momentum_BCs%val
        SUF_U_ROB1_BC_ALL=>velocity_BCs%val
        SUF_U_ROB2_BC_ALL=>velocity_BCs_robin2%val
        mom_conserv=.false.
        call get_option( &
            '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/conservative_advection', &
            beta )
        if (beta>=.999) mom_conserv=.true.
        ewrite(3,*) 'mom_conserv:', mom_conserv

        lump_mass = have_option( &
            '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/mass_terms/lump_mass_matrix')
        !retrieve lump_weight parameter
        call get_option( &
            '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/mass_terms/lump_mass_matrix/lump_weight', &
            lump_weight, default = -1. )
        !Act only if the parameter is above zero
        homogenize_mass_matrix = (lump_weight > 0)


        !For P1DGP1 or P1DGP2 using the new formulation this solves the problem with pressure boundary conditions
        !also, this requires to use the old way to get the Pivit Matrix
        if (Mmat%CV_pressure) then
            lump_mass = .true.
            homogenize_mass_matrix = .true.
            call get_option( &
            '/geometry/mesh::PressureMesh/from_mesh/mesh_shape/polynomial_degree', j )
            !For P1DGP1 the correct value is 100 and for P1DGP2 the correct value seems to be 10.
            lump_weight = 100.**(1./j)
        end if

        lump_absorption = .false.
        if ( have_option( &
            '/material_phase[0]/vector_field::Velocity/prognostic/vector_field::Absorption/lump_absorption') &
            ) lump_absorption = .true.
        lump_mass2 = .false.
        if ( lump_absorption ) lump_mass2 = .true.
        ! This applies a non-linear shock capturing scheme which
        ! may be used to reduce oscillations in velocity or
        ! perform implicit LES modelling of turbulence.
        ! In all residual approaches do not apply Petrov-Galerkin
        ! dissipation on the 1st non-linear iteration within a
        ! time step as there is no good guess of the (U^{n+1}-U^n)/DT.
        ! RESID_BASED_STAB_DIF decides what type of Petrov-Galerkin
        ! method to use.
        ! =1 is the residual squared approach.
        ! =2 is max(0, A . grad U * residual ).
        ! =3 is the max of 1 and 2 (the most dissipative).
        ! U_NONLIN_SHOCK_COEF \in [0,1] is the magnitude of the non-linear
        ! dissipation
        ! =0.25 is small
        ! =1.0 is large
        ! RNO_P_IN_A_DOT \in [0,1] decides if we include the pressure term in
        ! A . grad soln if
        ! =0.0 dont include pressure term.
        ! =1.0 include the pressure term.
        call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
            'spatial_discretisation/discontinuous_galerkin/stabilisation/method', &
            RESID_BASED_STAB_DIF, default=0 )
        BETWEEN_ELE_STAB = RESID_BASED_STAB_DIF/=0 ! Always switch on between element diffusion if using non-linear
        call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
            'spatial_discretisation/discontinuous_galerkin/stabilisation/nonlinear_velocity_coefficient', &
            U_NONLIN_SHOCK_COEF, default=1.)
        call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
            'spatial_discretisation/discontinuous_galerkin/stabilisation/include_pressure', &
            RNO_P_IN_A_DOT, default=1.)
        ewrite(3,*) 'RESID_BASED_STAB_DIF, U_NONLIN_SHOCK_COEF, RNO_P_IN_A_DOT:', &
            RESID_BASED_STAB_DIF, U_NONLIN_SHOCK_COEF, RNO_P_IN_A_DOT
        FEM_BUOYANCY = have_option( "/physical_parameters/gravity/fem_buoyancy" )
        GOT_DIFFUS = .FALSE.
        UPWIND_DGFLUX = .TRUE.
        if ( have_option( &
            '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/advection_scheme/central_differencing') &
            ) UPWIND_DGFLUX = .FALSE.
        NON_LIN_DGFLUX = .FALSE.
        if ( have_option( &
            '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/advection_scheme/nonlinear_flux') &
            ) NON_LIN_DGFLUX = .TRUE.
        ALLOCATE( UD( Mdims%ndim, Mdims%nphase, FE_GIdims%cv_ngi ))
        ALLOCATE( UDOLD( Mdims%ndim, Mdims%nphase, FE_GIdims%cv_ngi ))
        ALLOCATE( UD_ND( Mdims%ndim, Mdims%nphase, FE_GIdims%cv_ngi ))
        ALLOCATE( UDOLD_ND( Mdims%ndim, Mdims%nphase, FE_GIdims%cv_ngi ))
        ALLOCATE( DENGI( Mdims%nphase, FE_GIdims%cv_ngi ))
        ALLOCATE( DENGIOLD( Mdims%nphase, FE_GIdims%cv_ngi ))
        ALLOCATE( GRAD_SOU_GI( Mdims%ncomp, Mdims%nphase, FE_GIdims%cv_ngi ))
        ALLOCATE( GRAD_SOU2_GI( Mdims%ncomp, Mdims%nphase, FE_GIdims%cv_ngi ))
        ALLOCATE( RHS_U_CV( Mdims%nphase, Mdims%u_nloc ))
        ALLOCATE( RHS_U_CV_OLD( Mdims%nphase, Mdims%u_nloc ))
        ALLOCATE( UDEN_VFILT( Mdims%nphase, Mdims%u_nloc ))
        ALLOCATE( UDENOLD_VFILT( Mdims%nphase, Mdims%u_nloc ))
        ALLOCATE( SIGMAGI( Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase, FE_GIdims%cv_ngi ))
        ALLOCATE( SIGMAGI_STAB( Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase, FE_GIdims%cv_ngi ))
        ALLOCATE( XL_ALL(Mdims%ndim,Mdims%cv_nloc), XL2_ALL(Mdims%ndim,Mdims%cv_nloc), XSL_ALL(Mdims%ndim,Mdims%cv_snloc) )
        ALLOCATE( NORMX_ALL(Mdims%ndim), SNORMXN_ALL(Mdims%ndim,FE_GIdims%sbcvngi) )
        !Variables to improve Mmat%PIVIT_MAT creation speed
        !Initialization to zero is necesary for porous media, since the majority of the time they will be zero for all the elements
        ALLOCATE(NN_SIGMAGI_ELE( Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase, Mdims%u_nloc, Mdims%u_nloc )); NN_SIGMAGI_ELE = 0.
        ALLOCATE(NN_SIGMAGI_STAB_ELE( Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase, Mdims%u_nloc, Mdims%u_nloc )); NN_SIGMAGI_STAB_ELE = 0.
        ALLOCATE(NN_MASS_ELE( Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase, Mdims%u_nloc, Mdims%u_nloc )); NN_MASS_ELE = 0.
        ALLOCATE(NN_MASSOLD_ELE( Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase, Mdims%u_nloc, Mdims%u_nloc )); NN_MASSOLD_ELE = 0.
        ALLOCATE( STRESS_IJ_ELE( Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%u_nloc,Mdims%u_nloc )); STRESS_IJ_ELE = 0.
        ALLOCATE( VLK_ELE( Mdims%nphase, Mdims%u_nloc, Mdims%u_nloc )); VLK_ELE = 0.
        ALLOCATE( CENT_RELAX(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) ,CENT_RELAX_OLD(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        IF(GOT_VIRTUAL_MASS) THEN
            ! GOT_VIRTUAL_MASS ! do we have virtual mass terms for multi-phase flows...
            ! VIRTUAL_MASS_ADV_CUR DEFINES THE VELOCITY IN THE TOTAL DERIVATIVE = 1 and use the velocity that
            ! one is advecting, else =0 use the velocity of the current phase.
            ALLOCATE( VIRTUAL_MASS( Mdims%nphase, Mdims%nphase, Mdims%cv_nonods), VIRTUAL_MASS_OLD( Mdims%nphase, Mdims%nphase, Mdims%cv_nonods), VIRTUAL_MASS_ADV_CUR( Mdims%nphase, Mdims%nphase) )
            VIRTUAL_MASS=0.0
            VIRTUAL_MASS(2,1,:) = -UDEN(1,:)*0.5
            VIRTUAL_MASS(2,2,:) =  UDEN(1,:)*0.5
            VIRTUAL_MASS_OLD=VIRTUAL_MASS
            VIRTUAL_MASS_ADV_CUR=1.0
            ! For the time being VIRTUAL_MASS & VIRTUAL_MASS_OLD can be the same.
            ALLOCATE( LOC_VIRTUAL_MASS( Mdims%nphase, Mdims%nphase, Mdims%mat_nloc), LOC_VIRTUAL_MASS_OLD( Mdims%nphase, Mdims%nphase, Mdims%mat_nloc) )
            ALLOCATE( VIRTUAL_MASS_GI( Mdims%nphase, Mdims%nphase, FE_GIdims%CV_NGI), VIRTUAL_MASS_OLD_GI( Mdims%nphase, Mdims%nphase, FE_GIdims%CV_NGI))
            ALLOCATE( VLN_CVM( Mdims%nphase,Mdims%nphase ), VLN_OLD_CVM( Mdims%nphase,Mdims%nphase ) )
            ALLOCATE( CVM_SNDOTQ_IN(Mdims%ndim,Mdims%nphase,Mdims%nphase,FE_GIdims%sbcvngi), CVM_SNDOTQ_OUT(Mdims%ndim,Mdims%nphase,Mdims%nphase,FE_GIdims%sbcvngi) )
            ALLOCATE( CVM_SNDOTQOLD_IN(Mdims%ndim,Mdims%nphase,Mdims%nphase,FE_GIdims%sbcvngi), CVM_SNDOTQOLD_OUT(Mdims%ndim,Mdims%nphase,Mdims%nphase,FE_GIdims%sbcvngi) )
            ALLOCATE( CVM_NN_SNDOTQ_IN(Mdims%nphase), CVM_NN_SNDOTQ_OUT(Mdims%nphase), CVM_NN_SNDOTQOLD_IN(Mdims%nphase), CVM_NN_SNDOTQOLD_OUT(Mdims%nphase) )
            ALLOCATE( SVIRTUAL_MASS_GI(Mdims%nphase,Mdims%nphase,FE_GIdims%sbcvngi) , SVIRTUAL_MASS_OLD_GI(Mdims%nphase,Mdims%nphase,FE_GIdims%sbcvngi) )
            ALLOCATE( SVIRTUAL_MASS_GI_KEEP(Mdims%nphase,Mdims%nphase,FE_GIdims%sbcvngi), SVIRTUAL_MASS_GI2_KEEP(Mdims%nphase,Mdims%nphase,FE_GIdims%sbcvngi) )
            ALLOCATE( SVIRTUAL_MASS_OLD_GI_KEEP(Mdims%nphase,Mdims%nphase,FE_GIdims%sbcvngi), SVIRTUAL_MASS_OLD_GI2_KEEP(Mdims%nphase,Mdims%nphase,FE_GIdims%sbcvngi) )
            ALLOCATE( SLOC_VIRTUAL_MASS( Mdims%nphase, Mdims%nphase, Mdims%cv_snloc ), SLOC2_VIRTUAL_MASS( Mdims%nphase, Mdims%nphase, Mdims%cv_snloc ) )
            ALLOCATE( SLOC_VIRTUAL_MASS_OLD( Mdims%nphase, Mdims%nphase, Mdims%cv_snloc ), SLOC2_VIRTUAL_MASS_OLD( Mdims%nphase, Mdims%nphase, Mdims%cv_snloc ) )
            WITH_NONLIN_CVM=1.0
            ! CVM_BETA=0.0 (nonconservative virtual mass- standard),  =1. (conservative virtual mass)
            CVM_BETA=0.0
        ENDIF
        !sprint_to_do; make all this memory static and/or between if for the memory that it is not always required
        ALLOCATE( SDETWE( FE_GIdims%sbcvngi ))
        ALLOCATE( CV_SLOC2LOC( Mdims%cv_snloc ))
        ALLOCATE( U_SLOC2LOC( Mdims%u_snloc ))
        ALLOCATE( U_ILOC_OTHER_SIDE(Mdims%u_snloc))
        ALLOCATE( U_OTHER_LOC( Mdims%u_nloc ))
        ALLOCATE( MAT_OTHER_LOC( Mdims%mat_nloc ))
        ALLOCATE( TEN_XX( Mdims%ndim, Mdims%ndim, Mdims%nphase, FE_GIdims%cv_ngi ))
        ALLOCATE( TEN_VOL( Mdims%nphase, FE_GIdims%cv_ngi ))
        ALLOCATE( VLN( Mdims%nphase ))
        ALLOCATE( VLN_OLD( Mdims%nphase ))
        ALLOCATE( SUD_ALL(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SUDOLD_ALL(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SUD2_ALL(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SUDOLD2_ALL(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SNDOTQ(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SNDOTQOLD(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SINCOME(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SINCOMEOLD(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SDEN(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SDENOLD(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SDEN_KEEP(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SDENOLD_KEEP(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SDEN2_KEEP(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SDENOLD2_KEEP(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SUD_ALL_KEEP(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SUDOLD_ALL_KEEP(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SUD2_ALL_KEEP(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SUDOLD2_ALL_KEEP(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SNDOTQ_KEEP(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SNDOTQ2_KEEP(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SNDOTQOLD_KEEP(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SNDOTQOLD2_KEEP(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( N_DOT_DU(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( N_DOT_DU2(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( N_DOT_DUOLD(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( N_DOT_DUOLD2(Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( vel_dot(FE_GIdims%sbcvngi), vel_dot2(FE_GIdims%sbcvngi), velold_dot(FE_GIdims%sbcvngi), velold_dot2(FE_GIdims%sbcvngi), grad_fact(FE_GIdims%sbcvngi) )
        ALLOCATE( DIFF_COEF_DIVDX(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( DIFF_COEFOLD_DIVDX(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( FTHETA(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SNDOTQ_IN(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SNDOTQ_OUT(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SNDOTQOLD_IN(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( SNDOTQOLD_OUT(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( GRAD_SOU_GI_NMX( Mdims%ndim, Mdims%ncomp, Mdims%nphase ))
        ALLOCATE( GRAD_SOU2_GI_NMX( Mdims%ndim, Mdims%ncomp, Mdims%nphase ))
        GRAD_SOU_GI_NMX=0.0 ; GRAD_SOU2_GI_NMX=0.0
        ALLOCATE( MASS_ELE( Mdims%totele )) ; MASS_ELE=0.0
        ! Allocating for non-linear Petrov-Galerkin diffusion stabilization...
        ALLOCATE( LOC_MASS_INV(Mdims%u_nloc, Mdims%u_nloc) )
        ALLOCATE( LOC_MASS(Mdims%u_nloc, Mdims%u_nloc) )
        ALLOCATE( RHS_DIFF_U(Mdims%ndim,Mdims%nphase,Mdims%u_nloc) )
        ALLOCATE( DIFF_VEC_U(Mdims%ndim,Mdims%nphase,Mdims%u_nloc) )
        ALLOCATE( DIFFGI_U(Mdims%ndim,Mdims%nphase,FE_GIdims%cv_ngi) )
        ALLOCATE( U_DT(Mdims%ndim, Mdims%nphase,FE_GIdims%cv_ngi) )
        ALLOCATE( U_DX_ALL( Mdims%ndim, Mdims%ndim, Mdims%nphase, FE_GIdims%cv_ngi ) )
        ALLOCATE( UOLD_DX_ALL( Mdims%ndim, Mdims%ndim, Mdims%nphase, FE_GIdims%cv_ngi ) )
        ALLOCATE( SOUGI_X(Mdims%ndim,Mdims%nphase,FE_GIdims%cv_ngi) )
        ALLOCATE( RESID_U(Mdims%ndim,Mdims%nphase,FE_GIdims%cv_ngi) )
        ALLOCATE( P_DX(Mdims%ndim, FE_GIdims%cv_ngi)  )
        ALLOCATE( U_GRAD_NORM2(Mdims%ndim,Mdims%nphase,FE_GIdims%cv_ngi), U_GRAD_NORM(Mdims%ndim,Mdims%nphase,FE_GIdims%cv_ngi) )
        ALLOCATE( A_DOT_U(Mdims%ndim,Mdims%nphase,FE_GIdims%cv_ngi) )
        ALLOCATE( STAR_U_COEF(Mdims%ndim,Mdims%nphase,FE_GIdims%cv_ngi) )
        ALLOCATE( P_STAR_U(Mdims%ndim,Mdims%nphase,FE_GIdims%cv_ngi) )
        ALLOCATE( DIF_STAB_U(Mdims%ndim,Mdims%nphase, FE_GIdims%cv_ngi) )
        ALLOCATE( U_R2_COEF( Mdims%ndim ) )
        ALLOCATE( U_GRAD_N_MAX2( Mdims%ndim ) )
        ALLOCATE( VLK_UVW(Mdims%ndim) )
        ! Variables used to reduce indirect addressing...
        ALLOCATE( LOC_U(Mdims%ndim, Mdims%nphase, Mdims%u_nloc),  LOC_UOLD(Mdims%ndim, Mdims%nphase, Mdims%u_nloc) )
        IF(RETRIEVE_SOLID_CTY) THEN
            ALLOCATE( LOC_US(Mdims%ndim, Mdims%nphase, Mdims%u_nloc))
            ALLOCATE( LOC_U_ABS_STAB_SOLID_RHS(Mdims%ndim* Mdims%nphase, Mdims%ndim* Mdims%nphase, Mdims%mat_nloc))
        ENDIF
        ALLOCATE( LOC_NU(Mdims%ndim, Mdims%nphase, Mdims%u_nloc),  LOC_NUOLD(Mdims%ndim, Mdims%nphase, Mdims%u_nloc) )
        ALLOCATE( LOC_UDEN(Mdims%nphase, Mdims%cv_nloc),  LOC_UDENOLD(Mdims%nphase, Mdims%cv_nloc) )
        ALLOCATE( LOC_P(Mdims%p_nloc) )
        ALLOCATE( LOC_PLIKE_GRAD_SOU_COEF(Mdims%ncomp, Mdims%nphase, Mdims%cv_nloc) )
        ALLOCATE( LOC_PLIKE_GRAD_SOU_GRAD(Mdims%ncomp, Mdims%nphase, Mdims%cv_nloc) )
        ALLOCATE( LOC_U_SOURCE(Mdims%ndim, Mdims%nphase, Mdims%u_nloc) ) ; LOC_U_SOURCE=0.0
        ALLOCATE( LOC_U_SOURCE_CV(Mdims%ndim, Mdims%nphase, Mdims%cv_nloc) )
        ALLOCATE( LOC_U_ABSORB  (Mdims%ndim* Mdims%nphase, Mdims%ndim* Mdims%nphase, Mdims%mat_nloc) )
        ALLOCATE( LOC_U_ABS_STAB(Mdims%ndim* Mdims%nphase, Mdims%ndim* Mdims%nphase, Mdims%mat_nloc) )
        ALLOCATE( LOC_UDIFFUSION(Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%mat_nloc) )
        ALLOCATE( LOC_UDIFFUSION_VOL( Mdims%nphase, Mdims%mat_nloc) )
        ALLOCATE( LOC_U_RHS( Mdims%ndim, Mdims%nphase, Mdims%u_nloc ) )
        ! To memory access very local...
        ALLOCATE( SLOC_U(Mdims%ndim,Mdims%nphase,Mdims%u_snloc) )
        ALLOCATE( SLOC_UOLD(Mdims%ndim,Mdims%nphase,Mdims%u_snloc) )
        ALLOCATE( SLOC2_U(Mdims%ndim,Mdims%nphase,Mdims%u_snloc) )
        ALLOCATE( SLOC2_UOLD(Mdims%ndim,Mdims%nphase,Mdims%u_snloc) )
        ALLOCATE( SLOC_NU(Mdims%ndim,Mdims%nphase,Mdims%u_snloc) )
        ALLOCATE( SLOC_NUOLD(Mdims%ndim,Mdims%nphase,Mdims%u_snloc) )
        ALLOCATE( SLOC2_NU(Mdims%ndim,Mdims%nphase,Mdims%u_snloc) )
        ALLOCATE( SLOC2_NUOLD(Mdims%ndim,Mdims%nphase,Mdims%u_snloc) )
        ALLOCATE( SLOC_DUX_ELE_ALL( Mdims%ndim, Mdims%ndim , Mdims%nphase, Mdims%u_snloc ) )
        ALLOCATE( SLOC2_DUX_ELE_ALL( Mdims%ndim, Mdims%ndim , Mdims%nphase, Mdims%u_snloc ) )
        ALLOCATE( SLOC_DUOLDX_ELE_ALL( Mdims%ndim, Mdims%ndim , Mdims%nphase, Mdims%u_snloc ) )
        ALLOCATE( SLOC2_DUOLDX_ELE_ALL( Mdims%ndim, Mdims%ndim , Mdims%nphase, Mdims%u_snloc ) )
        ALLOCATE( SLOC_DIFF_FOR_BETWEEN_U(Mdims%ndim,Mdims%nphase,Mdims%u_snloc) )
        ALLOCATE( SLOC2_DIFF_FOR_BETWEEN_U(Mdims%ndim,Mdims%nphase,Mdims%u_snloc) )
        ALLOCATE( SLOC_UDEN(Mdims%nphase, Mdims%cv_snloc)  )
        ALLOCATE( SLOC2_UDEN(Mdims%nphase, Mdims%cv_snloc)  )
        ALLOCATE( SLOC_UDENOLD(Mdims%nphase, Mdims%cv_snloc)  )
        ALLOCATE( SLOC2_UDENOLD(Mdims%nphase, Mdims%cv_snloc) )
        ALLOCATE( SLOC_UDIFFUSION(Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%cv_snloc) )
        ALLOCATE( SLOC2_UDIFFUSION(Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%cv_snloc) )
        ALLOCATE( SLOC_UDIFFUSION_VOL( Mdims%nphase, Mdims%cv_snloc) ) ; SLOC_UDIFFUSION_VOL=0.0
        ALLOCATE( SLOC2_UDIFFUSION_VOL( Mdims%nphase, Mdims%cv_snloc) ) ; SLOC2_UDIFFUSION_VOL=0.0
        ! Derivatives...
        call allocate_multi_dev_shape_funs(FE_funs, Devfuns)
        ALLOCATE( NMX_ALL(Mdims%ndim) )
        ALLOCATE( VNMX_ALL(Mdims%ndim) )
        ALLOCATE( U_NODI_SGI_IPHASE_ALL(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( U_NODJ_SGI_IPHASE_ALL(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( UOLD_NODI_SGI_IPHASE_ALL(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        ALLOCATE( UOLD_NODJ_SGI_IPHASE_ALL(Mdims%ndim,Mdims%nphase,FE_GIdims%sbcvngi) )
        IF(IDIVID_BY_VOL_FRAC+IGOT_VOL_X_PRESSURE.GE.1) THEN
            ALLOCATE( VOL_FRA_GI(Mdims%nphase, FE_GIdims%CV_NGI), VOL_FRA_GI_DX_ALL( Mdims%ndim, Mdims%nphase, FE_GIdims%CV_NGI) )
            ALLOCATE( SLOC_VOL_FRA(Mdims%nphase, Mdims%cv_snloc), SLOC2_VOL_FRA(Mdims%nphase, Mdims%cv_snloc))
            ALLOCATE( SVOL_FRA(Mdims%nphase, FE_GIdims%sbcvngi), SVOL_FRA2(Mdims%nphase, FE_GIdims%sbcvngi) )
            ALLOCATE( VOL_FRA_NMX_ALL(Mdims%ndim,Mdims%nphase) )
        ENDIF
        GOT_DIFFUS = ( R2NORM( UDIFFUSION, Mdims%mat_nonods * Mdims%ndim * Mdims%ndim * Mdims%nphase ) /= 0.0 )  &
            .OR. ( UDIFFUSION_VOL%have_field ) .OR. BETWEEN_ELE_STAB
        IF(LES_DISOPT.NE.0) GOT_DIFFUS=.TRUE.
        IF(GOT_DIFFUS.AND.LINEAR_HIGHORDER_DIFFUSION) THEN
            ALLOCATE( STRESS_IJ_ELE_EXT( Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%u_snloc, 2*Mdims%u_nloc ) )
            ALLOCATE( S_INV_NNX_MAT12(Mdims%ndim, Mdims%u_snloc, 2*Mdims%u_nloc) )
            ALLOCATE( NNX_MAT_ELE(Mdims%ndim, Mdims%u_nloc, Mdims%u_nloc, Mdims%totele), NN_MAT_ELE(Mdims%u_nloc, Mdims%u_nloc, Mdims%totele) )
        ENDIF
        IF(RETRIEVE_SOLID_CTY) THEN
            ALLOCATE( SIGMAGI_STAB_SOLID_RHS( Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase, FE_GIdims%cv_ngi ))
            ALLOCATE(NN_SIGMAGI_STAB_SOLID_RHS_ELE( Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase, Mdims%u_nloc, Mdims%u_nloc ))
            IF( GOT_DIFFUS .AND. include_viscous_solid_fluid_drag_force ) THEN  ! include_viscous_solid_fluid_drag_force taken from diamond
                ALLOCATE(ABS_SOLID_FLUID_COUP(Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%cv_nonods))
                ALLOCATE(FOURCE_SOLID_FLUID_COUP(Mdims%ndim, Mdims%nphase, Mdims%cv_nonods))
                f_x => extract_vector_field( packed_state, "f_x" )
                do IDIMSF= 1, Mdims%ndim
                    FOURCE_SOLID_FLUID_COUP( IDIMSF, 1, : ) = 7.5*f_x%val(IDIMSF, : ) !f_x do not include Mdims%nphase
                end do
                a_xx => extract_tensor_field( packed_state, "a_xx")
                do IDIMSF=1, Mdims%ndim
                    do JDIMSF=1, Mdims%ndim
                        ABS_SOLID_FLUID_COUP(IDIMSF, JDIMSF, 1, :)=7.5*a_xx%val(IDIMSF, JDIMSF, :) ! a_xx do not include Mdims%nphase
                    end do
                end do
            ENDIF
        ENDIF
        GOT_UDEN = .FALSE.
        DO IPHASE = 1, Mdims%nphase
            GOT_UDEN = GOT_UDEN .OR. ( R2NORM( UDEN( IPHASE, : ), Mdims%cv_nonods ) /= 0.0 )
        END DO
        JUST_BL_DIAG_MAT=( ( .NOT. GOT_DIFFUS ) .AND. ( .NOT. GOT_UDEN ) )
        ALLOCATE( UDIFF_SUF_STAB( Mdims%ndim, Mdims%ndim, Mdims%ndim, Mdims%nphase, FE_GIdims%sbcvngi ) )
        UDIFF_SUF_STAB = 0.0
        IF ( BETWEEN_ELE_STAB ) THEN
            ! Calculate stabilization diffusion coefficient between elements...
            ALLOCATE( DIFF_FOR_BETWEEN_U( Mdims%ndim, Mdims%nphase, Mdims%u_nloc, Mdims%totele ) ) ; DIFF_FOR_BETWEEN_U = 0.0
            ALLOCATE( MAT_ELE( Mdims%u_nloc, Mdims%u_nloc, Mdims%totele ) ) ; MAT_ELE = 0.0
        END IF
        IF(RESID_BASED_STAB_DIF.NE.0) THEN
            IF( THERMAL_STAB_VISC ) THEN
                ALLOCATE( MAT_ELE_CV_LOC( Mdims%cv_nloc, Mdims%cv_nloc ), INV_MAT_ELE_CV_LOC( Mdims%cv_nloc, Mdims%cv_nloc ), DIFF_FOR_BETWEEN_CV( Mdims%ndim, Mdims%nphase, Mdims%cv_nloc ) )
                ALLOCATE( DIFFCV(Mdims%ndim, Mdims%nphase, Mdims%mat_nloc), DIFFCV_TEN(Mdims%ndim,Mdims%ndim, Mdims%nphase, Mdims%mat_nloc) )
                ALLOCATE( DIFFCV_TEN_ELE(Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%mat_nloc, Mdims%totele) ) ; DIFFCV_TEN_ELE=0.0
                ALLOCATE( RCOUNT_NODS(Mdims%mat_nonods) )
            END IF
        END IF

        get_gradU = .false.
        do iphase = 1, Mdims%nphase
           gradU => extract_tensor_field( state( iphase ), "gradU", stat )
           if ( stat == 0 ) then
              get_gradU = .true.
              exit
           end if
        end do

        IF ( GOT_DIFFUS .or. get_gradU ) THEN
            ALLOCATE( DUX_ELE_ALL( Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%u_nloc, Mdims%totele ) )
            ALLOCATE( DUOLDX_ELE_ALL( Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%u_nloc, Mdims%totele ) )
        ENDIF
        IF( (.NOT.JUST_BL_DIAG_MAT) .AND. (.NOT.Mmat%NO_MATRIX_STORE) ) call zero( Mmat%DGM_petsc )
        if (.not.got_c_matrix) Mmat%C = 0.0
        Mmat%U_RHS = 0.0
        IF (.NOT.Mmat%NO_MATRIX_STORE ) THEN!Only for inertia flow
            ALLOCATE( DIAG_BIGM_CON( Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%nphase, Mdims%u_nloc, Mdims%u_nloc, Mdims%totele ) )
            ALLOCATE( BIGM_CON( Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%nphase, Mdims%u_nloc, Mdims%u_nloc, Mspars%ELE%ncol ) )
            DIAG_BIGM_CON = 0.0
            BIGM_CON = 0.0
        END IF
        if ( got_free_surf ) then
            if ( symmetric_P ) then
                MASS_SUF=0.0
            end if
        end if
        quad_over_whole_ele = .true.
        ! ALLOCATE reversed ordering for computational speed****************
        ALLOCATE( CVFENX_ALL_REVERSED(Mdims%ndim,FE_GIdims%cv_ngi,Mdims%cv_nloc), UFENX_ALL_REVERSED(Mdims%ndim,FE_GIdims%cv_ngi,Mdims%u_nloc) ) ! NOT CALCULATED IN SUB cv_fem_shape_funs_plus_storage
        ALLOCATE( UFEN_REVERSED(FE_GIdims%cv_ngi,Mdims%u_nloc))
        ALLOCATE( CVN_REVERSED(FE_GIdims%cv_ngi,Mdims%cv_nloc), CVFEN_REVERSED(FE_GIdims%cv_ngi,Mdims%cv_nloc) )
        ALLOCATE( SBCVFEN_REVERSED(FE_GIdims%sbcvngi,Mdims%cv_snloc), SBUFEN_REVERSED(FE_GIdims%sbcvngi,Mdims%u_snloc) )
        DO U_ILOC=1,Mdims%u_nloc
            DO GI=1,FE_GIdims%cv_ngi
                UFEN_REVERSED(GI,U_ILOC) = FE_funs%ufen(U_ILOC,GI)
            END DO
        END DO
        DO CV_ILOC=1,Mdims%cv_nloc
            DO GI=1,FE_GIdims%cv_ngi
                CVFEN_REVERSED(GI,CV_ILOC)= FE_funs%cvfen(CV_ILOC,GI)
                !############THIS IS A FIX FOR GRAVITY######################
                !IF quad_over_whole_ele = .TRUE. NEEDS TO BE ACTIVATED
                CVN_REVERSED(GI,CV_ILOC)  = FE_funs%cvfen(CV_ILOC,GI)
                !###########################################################
                CVFEN_REVERSED(GI,CV_ILOC)= FE_funs%cvfen(CV_ILOC,GI)
            END DO
        END DO
        DO U_SILOC=1,Mdims%u_snloc
            DO SGI=1,FE_GIdims%sbcvngi
                SBUFEN_REVERSED(SGI,U_SILOC) = FE_funs%sbufen(U_SILOC,SGI)
            END DO
        END DO
        DO CV_SILOC=1,Mdims%cv_snloc
            DO SGI=1,FE_GIdims%sbcvngi
                SBCVFEN_REVERSED(SGI,CV_SILOC) = FE_funs%sbcvfen(CV_SILOC,SGI)
            END DO
        END DO
        ! ALLOCATE reversed ordering for computational speed****************
        ! Memory for rapid retreval...
        ! Storage for pointers to the other side of the element.
        ALLOCATE( STORED_U_ILOC_OTHER_SIDE( Mdims%u_snloc, FE_GIdims%nface, Mdims%totele*ISTORED_OTHER_SIDE ) )
        ALLOCATE( STORED_U_OTHER_LOC( Mdims%u_nloc, FE_GIdims%nface, Mdims%totele*ISTORED_OTHER_SIDE ) )
        ALLOCATE( STORED_MAT_OTHER_LOC( Mdims%mat_nloc, FE_GIdims%nface, Mdims%totele*ISTORED_OTHER_SIDE ) )
        ALLOCATE( POSINMAT_C_STORE( Mdims%u_nloc,Mdims%p_nloc, Mdims%totele*IDO_STORE_AC_SPAR_PT) )
        ALLOCATE( POSINMAT_C_STORE_SUF_DG( Mdims%u_snloc,Mdims%p_snloc,FE_GIdims%nface,Mdims%totele*IDO_STORE_AC_SPAR_PT ) )
        ALLOCATE( FACE_ELE( FE_GIdims%nface, Mdims%totele ) )
        ! Calculate FACE_ELE
        CALL CALC_FACE_ELE( FACE_ELE, Mdims%totele, Mdims%stotel, FE_GIdims%nface, &
            Mspars%ELE%fin, Mspars%ELE%col, Mdims%cv_nloc, Mdims%cv_snloc, Mdims%cv_nonods, ndgln%cv, ndgln%suf_cv, &
            FE_funs%cv_sloclist, Mdims%x_nloc, ndgln%x )

       IF( GOT_DIFFUS .or. get_gradU ) THEN
            CALL DG_DERIVS_ALL( U_ALL, UOLD_ALL, &
                DUX_ELE_ALL, DUOLDX_ELE_ALL, &
                Mdims%ndim, Mdims%nphase, Mdims%ndim, Mdims%totele, ndgln%u, &
                ndgln%xu, Mdims%x_nloc, ndgln%x, &
                FE_GIdims%cv_ngi, Mdims%u_nloc, FE_funs%cvweight, &
                FE_funs%ufen, FE_funs%ufenlx_all(1,:,:), FE_funs%ufenlx_all(2,:,:), FE_funs%ufenlx_all(3,:,:), &
                FE_funs%cvfen, FE_funs%cvfenlx_all(1,:,:), FE_funs%cvfenlx_all(2,:,:), FE_funs%cvfenlx_all(3,:,:), &
                Mdims%x_nonods, X_ALL(1,:), X_ALL(2,:), X_ALL(3,:), &
                FE_GIdims%nface, FACE_ELE, FE_funs%u_sloclist, FE_funs%cv_sloclist, Mdims%u_snloc, Mdims%cv_snloc, WIC_U_BC_ALL_VISC, SUF_U_BC_ALL_VISC, &
                FE_GIdims%sbcvngi, FE_funs%sbufen, FE_funs%sbcvfeweigh, &
                FE_funs%sbcvfen, FE_funs%sbcvfenslx, FE_funs%sbcvfensly, get_gradU, state )
        ENDIF
        ! LES VISCOCITY CALC.
        IF ( GOT_DIFFUS ) THEN
            ALLOCATE(UDIFFUSION_ALL(Mdims%ndim,Mdims%ndim,Mdims%nphase,Mdims%mat_nonods)) ; UDIFFUSION_ALL=0.
            ALLOCATE(UDIFFUSION_VOL_ALL(Mdims%nphase,Mdims%mat_nonods)) ; UDIFFUSION_VOL_ALL=0.
            IF ( LES_DISOPT /= 0 ) THEN
                ALLOCATE(LES_UDIFFUSION(Mdims%ndim,Mdims%ndim,Mdims%nphase,Mdims%mat_nonods)) ; LES_UDIFFUSION=0.
                ALLOCATE(LES_UDIFFUSION_VOL(Mdims%nphase,Mdims%mat_nonods)) ; LES_UDIFFUSION_VOL=0.
                CALL VISCOCITY_TENSOR_LES_CALC( LES_UDIFFUSION, LES_UDIFFUSION_VOL, LES_THETA*DUX_ELE_ALL + (1.-LES_THETA)*DUOLDX_ELE_ALL, &
                    Mdims%ndim,Mdims%nphase, Mdims%u_nloc,Mdims%x_nloc,Mdims%totele, Mdims%x_nonods, &
                    X_ALL, ndgln%x,  Mdims%mat_nonods, Mdims%mat_nloc, ndgln%mat, LES_DISOPT, LES_CS, UDEN, Mdims%cv_nonods, ndgln%cv, &
                    ndgln%u, Mdims%u_nonods, LES_THETA*U_ALL + (1.-LES_THETA)*UOLD_ALL, DERIV )
                IF ( STRESS_FORM ) THEN ! put into viscocity in stress form
                    DO IDIM=1,Mdims%ndim
                        DO JDIM=1,Mdims%ndim
                            UDIFFUSION_ALL(IDIM,JDIM,:,:) = UDIFFUSION(IDIM,JDIM,:,:) + SQRT( LES_UDIFFUSION(IDIM,IDIM,:,:) * LES_UDIFFUSION(JDIM,JDIM,:,:) )
                        END DO
                    END DO
                ELSE
                    UDIFFUSION_ALL=UDIFFUSION + LES_UDIFFUSION
                ENDIF
                !UDIFFUSION_VOL_ALL=UDIFFUSION_VOL + LES_UDIFFUSION_VOL
                if ( UDIFFUSION_VOL%have_field ) UDIFFUSION_VOL_ALL = UDIFFUSION_VOL%val(:,1,1,:)
if ( have_option( "/magma" ) ) then
   sfield => extract_scalar_field( state(1), "VolumetricViscosity" ) ! this should be on a material mesh
   UDIFFUSION_VOL_ALL(2,:) = sfield%val
end if
                UDIFFUSION_VOL_ALL = UDIFFUSION_VOL_ALL + LES_UDIFFUSION_VOL
            ELSE
                UDIFFUSION_ALL=UDIFFUSION
                if ( UDIFFUSION_VOL%have_field ) UDIFFUSION_VOL_ALL = UDIFFUSION_VOL%val(:,1,1,:)
if ( have_option( "/magma" ) ) then
   sfield => extract_scalar_field( state(1), "Ksi_s" ) ! this is the volumetric viscosity
   UDIFFUSION_VOL_ALL(2,:) = sfield%val                ! and it should be on a material mesh

end if
            ENDIF
        ENDIF
        if( RETRIEVE_SOLID_CTY ) THEN
            sf=> extract_scalar_field( packed_state, "SolidConcentration" )
            delta_u_all => extract_vector_field( packed_state, "delta_U" ) ! this is delta_u
            us_all => extract_vector_field( packed_state, "solid_U" )
            if ( .false. ) then ! Do not switch this to false.
                DO ELE = 1, Mdims%totele
                    DO CV_ILOC = 1, Mdims%cv_nloc
                        MAT_NOD = ndgln%mat( (ELE-1)*Mdims%cv_nloc + CV_ILOC )
                        CV_NOD = ndgln%cv( (ELE-1)*Mdims%cv_nloc + CV_ILOC )
                        DO IPHASE = 1, Mdims%nphase
                            UDIFFUSION_ALL( :, :, IPHASE, MAT_NOD ) = UDIFFUSION_ALL( :, :, IPHASE, MAT_NOD ) * ( 1. - sf%val( cv_nod ) )
                            UDIFFUSION_VOL_ALL( IPHASE, MAT_NOD ) = UDIFFUSION_VOL_ALL( IPHASE, MAT_NOD ) * ( 1. - sf%val( cv_nod ) )
                        END DO
                    END DO
                END DO
            end if
            allocate( vol_s_gi( FE_GIdims%cv_ngi ) )
            allocate( cv_dengi( Mdims%nphase, FE_GIdims%cv_ngi ) )
        endif

        ! surface tension-like terms
        IF ( IPLIKE_GRAD_SOU /= 0 ) THEN

            ALLOCATE( PLIKE_GRAD_SOU_GRAD(Mdims%ncomp, Mdims%nphase, Mdims%cv_nonods ) ) ; PLIKE_GRAD_SOU_GRAD=0.0
            ALLOCATE( PLIKE_GRAD_SOU_COEF(Mdims%ncomp, Mdims%nphase, Mdims%cv_nonods ) ) ; PLIKE_GRAD_SOU_COEF=0.0

            IF ( IPLIKE_GRAD_SOU == 1 ) THEN
                PLIKE_GRAD_SOU_GRAD_F => EXTRACT_TENSOR_FIELD( PACKED_STATE, "SurfaceTensionGrad" )
                PLIKE_GRAD_SOU_COEF_F => EXTRACT_TENSOR_FIELD( PACKED_STATE, "SurfaceTensionCoef" )
                PLIKE_GRAD_SOU_GRAD = PLIKE_GRAD_SOU_GRAD_F%val ; PLIKE_GRAD_SOU_COEF = PLIKE_GRAD_SOU_COEF_F%val
            ELSE IF ( IPLIKE_GRAD_SOU == 2 ) THEN
                ! this is for magma runs only at this point - need to generalise
                P_s_hat => EXTRACT_SCALAR_FIELD( STATE(1), "P_s_hat" ) ! PLIKE_GRAD_SOU_GRAD_F
                a_s => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedPhaseVolumeFraction" ) ! PLIKE_GRAD_SOU_COEF_F
                PLIKE_GRAD_SOU_GRAD(1,2,:) = P_s_hat%val ; PLIKE_GRAD_SOU_COEF(1,2,:) = a_s%val(1,1,:)
            END IF
        END IF


        Loop_Elements: DO ELE = 1, Mdims%totele ! VOLUME integral
            ! Calculate DevFuns%DETWEI,DevFuns%RA,NX,NY,NZ for element ELE
             call DETNLXR_PLUS_U(ELE, X_ALL, ndgln%x, FE_funs%cvweight, &
                FE_funs%cvfen, FE_funs%cvfenlx_all, FE_funs%ufenlx_all, Devfuns)
            DO GI = 1, FE_GIdims%CV_NGI
                CVFENX_ALL_REVERSED(:,GI,:) = DevFuns%CVFENX_ALL(:,:,GI)
                UFENX_ALL_REVERSED(:,GI,:) = DevFuns%UFENX_ALL(:,:,GI)
            END DO
            ! Adjust the DevFuns%VOLUME according to the number of levels.
            MASS_ELE( ELE ) = DevFuns%VOLUME
            if (IsParallel()) then
                if (.not. assemble_ele(pressure,ele)) then
                    skip=.true.
                    neighbours=>ele_neigh(pressure,ele)
                    do nb=1,size(neighbours)
                        if (neighbours(nb)<=0) cycle
                        if (assemble_ele(pressure,neighbours(nb))) then
                            skip=.false.
                            exit
                        end if
                    end do
                end if
            end if
            !Default mass matrix for porous media
            if (use_simple_lumped_homogenenous_mass_matrix) then
                call get_porous_Mass_matrix(ELE, Mdims, FE_funs, FE_GIdims, DevFuns, Mmat, WIC_P_BC_ALL, FACE_ELE)
            end if

            ! *********subroutine Determine local vectors...
            LOC_U_RHS = 0.0
            DO U_ILOC = 1, Mdims%u_nloc
                U_INOD = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                DO IPHASE = 1, Mdims%nphase
                    DO IDIM = 1, Mdims%ndim
                        LOC_U( IDIM, IPHASE, U_ILOC ) = U_ALL( IDIM, IPHASE, U_INOD )
                        LOC_UOLD( IDIM, IPHASE, U_ILOC ) = UOLD_ALL( IDIM, IPHASE, U_INOD )
                        if ( u_source%have_field ) LOC_U_SOURCE( IDIM, IPHASE, U_ILOC ) = U_SOURCE%val( IDIM, IPHASE, 1, U_INOD )
                        IF(RETRIEVE_SOLID_CTY) THEN
                            LOC_US( IDIM, IPHASE, U_ILOC ) = us_all%val( IDIM, U_INOD )
                        ENDIF
                    END DO
                END DO
            END DO
            DO U_ILOC = 1, Mdims%u_nloc
                U_INOD = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                DO IPHASE = 1, Mdims%nphase
                    DO IDIM = 1, Mdims%ndim
                        LOC_NU( IDIM, IPHASE, U_ILOC ) = NU_ALL( IDIM, IPHASE, U_INOD )
                        LOC_NUOLD( IDIM, IPHASE, U_ILOC ) = NUOLD_ALL( IDIM, IPHASE, U_INOD )
                    END DO
                END DO
            END DO
            DO CV_ILOC = 1, Mdims%cv_nloc
                CV_INOD = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
                IF(IGOT_VOL_X_PRESSURE==1) THEN
                    LOC_UDEN( :, CV_ILOC ) = UDEN( :, CV_INOD ) * FEM_VOL_FRAC( :, CV_INOD )
                    LOC_UDENOLD( :, CV_ILOC) = UDENOLD( :, CV_INOD ) * FEM_VOL_FRAC( :, CV_INOD )
                ELSE
                    LOC_UDEN( :, CV_ILOC ) = UDEN( :, CV_INOD )
                    LOC_UDENOLD( :, CV_ILOC) = UDENOLD( :, CV_INOD )
                ENDIF
                if (is_flooding) then
                    LOC_UDEN( :, CV_ILOC ) = 1.0
                    LOC_UDENOLD( :, CV_ILOC) = 1.0
                end if
                IF(GOT_VIRTUAL_MASS) THEN
                    LOC_VIRTUAL_MASS( :,:, CV_ILOC )         = VIRTUAL_MASS( :,:, CV_INOD )
                    LOC_VIRTUAL_MASS_OLD( :,:, CV_ILOC )     = VIRTUAL_MASS_OLD( :,:, CV_INOD )
                ENDIF
                DO IPHASE = 1, Mdims%nphase
                    IF ( IPLIKE_GRAD_SOU >= 1) THEN
                        LOC_PLIKE_GRAD_SOU_COEF( :, IPHASE, CV_ILOC ) = PLIKE_GRAD_SOU_COEF( :, IPHASE, CV_INOD )
                        LOC_PLIKE_GRAD_SOU_GRAD( :, IPHASE, CV_ILOC ) = PLIKE_GRAD_SOU_GRAD( :, IPHASE, CV_INOD )
                    END IF
                    DO IDIM = 1, Mdims%ndim
                        LOC_U_SOURCE_CV( IDIM, IPHASE, CV_ILOC ) = U_SOURCE_CV( IDIM, IPHASE, CV_INOD )
                    END DO
                END DO
            END DO
            DO P_ILOC = 1, Mdims%p_nloc
                P_INOD = ndgln%p( ( ELE - 1 ) * Mdims%p_nloc + P_ILOC )
                LOC_P( P_ILOC ) = P( 1,1,P_INOD )
            END DO
            LOC_U_ABSORB = 0.0
            IF(RETRIEVE_SOLID_CTY) LOC_U_ABS_STAB_SOLID_RHS=0.0
            DO MAT_ILOC = 1, Mdims%mat_nloc
                MAT_INOD = ndgln%mat( ( ELE - 1 ) * Mdims%mat_nloc + MAT_ILOC )
                IF(is_porous_media) THEN ! Set to the identity - NOT EFFICIENT BUT GOOD ENOUGH FOR NOW AS ITS SIMPLE...
                    DO I=1,Mdims%ndim* Mdims%nphase
                        LOC_U_ABSORB( I, I, MAT_ILOC ) = 1.0
                    END DO
                ELSE
                    LOC_U_ABSORB( :, :, MAT_ILOC ) = U_ABSORB( :, :, MAT_INOD )
                    ! Switch on for solid fluid-coupling...
                    IF(RETRIEVE_SOLID_CTY) THEN
                        CV_INOD = ndgln%cv( ( ELE - 1 ) * Mdims%mat_nloc + MAT_ILOC )
                        if(.false.) then ! delete this...
                            DO IDIM=1,Mdims%ndim
                                DO IPHASE=1,Mdims%nphase
                                    I=IDIM + (IPHASE-1)*Mdims%ndim
                                    LOC_U_ABSORB( I, I, MAT_ILOC ) = LOC_U_ABSORB( I, I, MAT_ILOC ) + &
                                        COEFF_SOLID_FLUID_stab * ( min ( 1.0, DEN_ALL( IPHASE, cv_inod )) / dt ) * sf%val( cv_inod )
                                    ! not used...
                                    LOC_U_ABS_STAB_SOLID_RHS( I, I, MAT_ILOC ) = LOC_U_ABS_STAB_SOLID_RHS( I, I, MAT_ILOC ) &
                                        + COEFF_SOLID_FLUID_stab * ( min ( 1.0, DEN_ALL( IPHASE, cv_inod ) ) / dt )
                                END DO
                            END DO
                        endif
                        ! Add in the viscocity contribution...
                        IF( GOT_DIFFUS .AND. include_viscous_solid_fluid_drag_force ) THEN
                            ! Assume visc. is isotropic (can be variable)...
                            DO IDIM=1,Mdims%ndim
                                DO IPHASE=1,Mdims%nphase
                                    I=IDIM + (IPHASE-1)*Mdims%ndim
                                    DO JDIM=1,Mdims%ndim
                                        J=JDIM + (IPHASE-1)*Mdims%ndim
                                        LOC_U_ABSORB( I, J, MAT_ILOC ) = LOC_U_ABSORB( I, J, MAT_ILOC ) &
                                            + ABS_SOLID_FLUID_COUP(IDIM, JDIM, IPHASE, CV_INOD)* UDIFFUSION_ALL( 1, 1, IPHASE, MAT_INOD )
                                    END DO
                                    LOC_U_SOURCE_CV( IDIM, IPHASE, MAT_ILOC ) = LOC_U_SOURCE_CV( IDIM, IPHASE, MAT_ILOC ) &
                                        + FOURCE_SOLID_FLUID_COUP(IDIM, IPHASE, CV_INOD)* UDIFFUSION_ALL( 1, 1, IPHASE, MAT_INOD )
                                END DO
                            END DO
                        ENDIF
                    ! ENDOF IF(RETRIEVE_SOLID_CTY) THEN...
                    ENDIF
                END IF
                LOC_U_ABS_STAB( :, :, MAT_ILOC ) = 0.
                ! Switch on for solid fluid-coupling apply stabilization term...
                IF(RETRIEVE_SOLID_CTY) THEN
                    CV_INOD = ndgln%cv( ( ELE - 1 ) * Mdims%mat_nloc + MAT_ILOC )
                    DO IDIM=1,Mdims%ndim
                        DO IPHASE=1,Mdims%nphase
                            I=IDIM + (IPHASE-1)*Mdims%ndim
                            LOC_U_ABS_STAB( I, I, MAT_ILOC ) = LOC_U_ABS_STAB( I, I, MAT_ILOC ) + &
                                COEFF_SOLID_FLUID_stab *( DEN_ALL( IPHASE, cv_inod ) / dt ) * sf%val( cv_inod )
                        END DO
                    END DO
                ENDIF
                IF ( GOT_DIFFUS ) THEN
                    LOC_UDIFFUSION( :, :, :, MAT_ILOC ) = UDIFFUSION_ALL( :, :, :, MAT_INOD )
                    LOC_UDIFFUSION_VOL( :, MAT_ILOC ) = UDIFFUSION_VOL_ALL( :, MAT_INOD )
                ELSE
                    LOC_UDIFFUSION( :, :, :, MAT_ILOC ) = 0.0
                    LOC_UDIFFUSION_VOL( :, MAT_ILOC ) = 0.0
                ENDIF
            END DO
            ! *********subroutine Determine local vectors...
            UD = 0.0 ; UDOLD = 0.0
            UD_ND = 0.0 ; UDOLD_ND = 0.0
            DO U_ILOC = 1, Mdims%u_nloc
                DO GI = 1, FE_GIdims%CV_NGI
                    UD( :, :, GI ) = UD( :, :, GI ) + UFEN_REVERSED( GI, U_ILOC ) * LOC_NU( :, :, U_ILOC )
                    UDOLD( :, :, GI ) = UDOLD( :, :, GI ) + UFEN_REVERSED( GI, U_ILOC ) * LOC_NUOLD( :, :, U_ILOC )
                END DO
            END DO
            UD_ND( 1:Mdims%ndim, :, : ) = UD
            UDOLD_ND( 1:Mdims%ndim, :, : ) = UDOLD
            IF(IDIVID_BY_VOL_FRAC+IGOT_VOL_X_PRESSURE.GE.1) THEN
                VOL_FRA_GI_DX_ALL=0.0
                VOL_FRA_GI=0.0
                DO CV_ILOC = 1, Mdims%cv_nloc
                    CV_INOD = ndgln%cv( (ELE-1)*Mdims%cv_nloc + CV_ILOC )
                    DO GI = 1, FE_GIdims%CV_NGI
                        DO IPHASE=1,Mdims%nphase
                            VOL_FRA_GI( IPHASE, GI )           = VOL_FRA_GI( IPHASE,GI )            + CVFEN_REVERSED( GI, CV_ILOC )       * FEM_VOL_FRAC( IPHASE, CV_INOD )
                            VOL_FRA_GI_DX_ALL( :, IPHASE, GI ) = VOL_FRA_GI_DX_ALL( :, IPHASE, GI ) + CVFENX_ALL_REVERSED( 1:Mdims%ndim, GI, CV_ILOC )* FEM_VOL_FRAC( IPHASE, CV_INOD )
                        END DO
                    END DO
                END DO
                VOL_FRA_GI=MAX(VOL_FRA_GI, 0.0)
            ENDIF
            DENGI = 0.0 ; DENGIOLD = 0.0
            GRAD_SOU_GI = 0.0 ; GRAD_SOU2_GI = 0.0
            IF(GOT_VIRTUAL_MASS) THEN
                VIRTUAL_MASS_GI = 0.0
                VIRTUAL_MASS_OLD_GI = 0.0
            ENDIF
            DO CV_ILOC = 1, Mdims%cv_nloc
                DO GI = 1, FE_GIdims%CV_NGI
                    IF ( FEM_DEN ) then ! FEM DEN...
                        DENGI( :, GI ) = DENGI( :, GI ) + CVFEN_REVERSED( GI, CV_ILOC ) * LOC_UDEN( :, CV_ILOC )
                        DENGIOLD( :, GI ) = DENGIOLD( :, GI ) &
                            + CVFEN_REVERSED( GI, CV_ILOC ) * LOC_UDENOLD( :, CV_ILOC )
                    ELSE ! CV DEN...
                        DENGI( :, GI ) = DENGI( :, GI ) + CVN_REVERSED( GI, CV_ILOC ) * LOC_UDEN( :, CV_ILOC )
                        DENGIOLD( :, GI ) = DENGIOLD( :, GI ) &
                            + CVN_REVERSED( GI, CV_ILOC ) * LOC_UDENOLD( :, CV_ILOC )
                    END IF
                    IF(GOT_VIRTUAL_MASS) THEN
                        IF ( FEM_DEN ) then ! FEM DEN...
                            VIRTUAL_MASS_GI( :,:, GI )         = VIRTUAL_MASS_GI( :,:, GI )         + CVFEN_REVERSED( GI, CV_ILOC ) * LOC_VIRTUAL_MASS( :,:, CV_ILOC )
                            VIRTUAL_MASS_OLD_GI( :,:, GI )     = VIRTUAL_MASS_OLD_GI( :,:, GI )     + CVFEN_REVERSED( GI, CV_ILOC ) * LOC_VIRTUAL_MASS_OLD( :,:, CV_ILOC )
                        ELSE
                            VIRTUAL_MASS_GI( :,:, GI )         = VIRTUAL_MASS_GI( :,:, GI )         + CVN_REVERSED( GI, CV_ILOC ) * LOC_VIRTUAL_MASS( :,:, CV_ILOC )
                            VIRTUAL_MASS_OLD_GI( :,:, GI )     = VIRTUAL_MASS_OLD_GI( :,:, GI )     + CVN_REVERSED( GI, CV_ILOC ) * LOC_VIRTUAL_MASS_OLD( :,:, CV_ILOC )
                        ENDIF
                    ENDIF
                    IF ( IPLIKE_GRAD_SOU >= 1 ) THEN
                        GRAD_SOU_GI( :, :, GI ) = GRAD_SOU_GI( :, :, GI ) &
                            + CVFEN_REVERSED( GI, CV_ILOC ) * LOC_PLIKE_GRAD_SOU_COEF( :, :, CV_ILOC )
                        IF ( IPLIKE_GRAD_SOU == 2 ) then
                           GRAD_SOU2_GI( :, :, GI ) = GRAD_SOU2_GI( :, :, GI ) &
                            + 2.*CVFEN_REVERSED( GI, CV_ILOC ) * LOC_PLIKE_GRAD_SOU_GRAD( :, :, CV_ILOC )
                        END IF
                    END IF
                END DO
            END DO
            ! Start filtering density
            !FILT_DEN = 1
            !FILT_DEN = 2 ! best option to use
            !Don't remove the code comment below!!!
            FILT_DEN = 0
!            IF ( FILT_DEN /= 0 ) THEN ! Filter the density...
!                DENGI = 0.0 ; DENGIOLD = 0.0
!                MASS_U = 0.0 ; MASS_U_CV = 0.0
!                DO U_ILOC = 1, Mdims%u_nloc
!                    DO U_JLOC = 1, Mdims%u_nloc
!                        NN = SUM( FE_funs%ufen( U_ILOC, : ) * FE_funs%ufen( U_JLOC, : ) * DevFuns%DETWEI(:) )
!                        IF ( FILT_DEN==2 ) THEN ! Lump the mass matrix for the filter - positive density...
!                            MASS_U( U_ILOC, U_ILOC ) = MASS_U( U_ILOC, U_ILOC ) + NN
!                        ELSE
!                            MASS_U( U_ILOC, U_JLOC ) = MASS_U( U_ILOC, U_JLOC ) + NN
!                        END IF
!                    END DO
!                END DO
!                DO U_ILOC = 1, Mdims%u_nloc
!                    DO CV_JLOC = 1, Mdims%cv_nloc
!                        NCVM = SUM( FE_funs%ufen( U_ILOC, : ) * CVN_SHORT( CV_JLOC, : ) * DevFuns%DETWEI(:) )
!                        MASS_U_CV( U_ILOC, CV_JLOC ) = MASS_U_CV( U_ILOC, CV_JLOC ) + NCVM
!                    END DO
!                END DO
!
!                STORE_MASS_U=MASS_U
!                ! Store the LU decomposition...
!                GOTDEC = .FALSE.
!
!                RHS_U_CV = 0.0 ; RHS_U_CV_OLD = 0.0
!                DO CV_JLOC = 1, Mdims%cv_nloc
!                    DO U_ILOC = 1, Mdims%u_nloc
!                        RHS_U_CV( :, U_ILOC ) = RHS_U_CV( :, U_ILOC ) + MASS_U_CV( U_ILOC, CV_JLOC ) * LOC_UDEN( :, CV_JLOC )
!                        RHS_U_CV_OLD( :, U_ILOC ) = RHS_U_CV_OLD( :, U_ILOC ) + MASS_U_CV( U_ILOC, CV_JLOC ) * LOC_UDENOLD( :, CV_JLOC )
!                    END DO
!                END DO
!
!                DO IPHASE = 1, Mdims%nphase
!                    CALL SMLINNGOT( STORE_MASS_U, UDEN_VFILT( IPHASE, : ), RHS_U_CV( IPHASE, : ), Mdims%u_nloc, IPIV, GOTDEC )
!                    GOTDEC = .TRUE.
!                    CALL SMLINNGOT( STORE_MASS_U, UDENOLD_VFILT( IPHASE, : ), RHS_U_CV_OLD( IPHASE, : ), Mdims%u_nloc, IPIV, GOTDEC )
!                END DO
!
!                DO U_ILOC = 1, Mdims%u_nloc
!                    DO GI = 1, CV_NGI
!                        DENGI( :, GI ) = DENGI( :, GI ) + FE_funs%ufen( U_ILOC, GI ) * UDEN_VFILT( :, U_ILOC )
!                        DENGIOLD( :, GI ) = DENGIOLD( :, GI ) + FE_funs%ufen( U_ILOC, GI ) * UDENOLD_VFILT( :, U_ILOC )
!                    END DO
!                END DO
!            END IF
! not good to have -ve density at quadature pt...
            DENGI = MAX( 0.0, DENGI )
            DENGIOLD = MAX( 0.0, DENGIOLD )
            SIGMAGI = 0.0 ; SIGMAGI_STAB = 0.0
            TEN_XX  = 0.0 ; TEN_VOL  = 0.0
            if (is_porous_media) then
                DO IPHA_IDIM = 1, Mdims%ndim * Mdims%nphase
                    SIGMAGI( IPHA_IDIM, IPHA_IDIM, : ) = 1.0
                end do
            else
                DO MAT_ILOC = 1, Mdims%mat_nloc
                    DO GI = 1, FE_GIdims%cv_ngi
                        DO IPHA_IDIM = 1, Mdims%ndim * Mdims%nphase
                            DO JPHA_JDIM = 1, Mdims%ndim * Mdims%nphase
                                SIGMAGI( IPHA_IDIM, JPHA_JDIM, GI ) = SIGMAGI( IPHA_IDIM, JPHA_JDIM, GI ) &
                                    + CVN_REVERSED( GI, MAT_ILOC ) * LOC_U_ABSORB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC )
                                SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM, GI ) = SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM, GI ) &
                                    + CVN_REVERSED( GI, MAT_ILOC ) * LOC_U_ABS_STAB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC )
                            END DO
                        END DO
                        TEN_XX( :, :, :, GI ) = TEN_XX( :, :, :, GI ) + CVFEN_REVERSED( GI, MAT_ILOC ) * LOC_UDIFFUSION( :, :, :, MAT_ILOC )
                        TEN_VOL( :, GI )      = TEN_VOL(  :, GI )     + CVFEN_REVERSED( GI, MAT_ILOC ) * LOC_UDIFFUSION_VOL( :, MAT_ILOC )
                    END DO
                END DO
                IF ( RETRIEVE_SOLID_CTY ) THEN
                    VOL_S_GI = 0.0
                    CV_DENGI = 0.0
                    DO CV_ILOC = 1, Mdims%cv_nloc
                        CV_INOD = ndgln%cv( (ELE-1)*Mdims%cv_nloc + CV_ILOC )
                        DO GI = 1, FE_GIdims%CV_NGI
                            !VOL_S_GI( GI ) = VOL_S_GI( GI ) + CVFEN_REVERSED( GI, CV_ILOC ) * sf%val( cv_inod )
                            VOL_S_GI( GI ) = VOL_S_GI( GI ) + CVN_REVERSED( GI, CV_ILOC ) * sf%val( cv_inod )
                            CV_DENGI(:, GI ) = CV_DENGI(:, GI ) + CVN_REVERSED( GI, CV_ILOC ) * den_all( :, cv_inod )
                        END DO
                    END DO
                    VOL_S_GI = MIN( MAX( VOL_S_GI, 0.0 ), 1.0 )
                    SIGMAGI_STAB_SOLID_RHS=0.0
                    do idim=1,Mdims%ndim
                        do iphase=1,Mdims%nphase
                            ipha_idim=idim + (iphase-1)*Mdims%ndim
                            SIGMAGI( IPHA_IDIM, IPHA_IDIM, : ) = SIGMAGI( IPHA_IDIM, IPHA_IDIM, : ) + COEFF_SOLID_FLUID_relax*max( CV_dengi( iphase, : ), min_den_for_solid_fluid ) * vol_s_gi(:) / dt
                            SIGMAGI_STAB_SOLID_RHS( IPHA_IDIM, IPHA_IDIM, : ) = SIGMAGI_STAB_SOLID_RHS( IPHA_IDIM, IPHA_IDIM, : ) + COEFF_SOLID_FLUID_relax*max( CV_dengi( iphase, : ), min_den_for_solid_fluid ) / dt
                        end do
                    end do
                end if
            end if
            RHS_DIFF_U=0.0
            if (Porous_media_PIVIT_not_stored_yet) then
                NN_SIGMAGI_ELE = 0.0
                NN_SIGMAGI_STAB_ELE = 0.0
                IF(RETRIEVE_SOLID_CTY) NN_SIGMAGI_STAB_SOLID_RHS_ELE = 0.0
                NN_MASS_ELE = 0.0
                NN_MASSOLD_ELE = 0.0
            end if
            VLK_ELE = 0.0
            STRESS_IJ_ELE = 0.0
            !Prepare data
            IF ( GOT_DIFFUS ) THEN
                DO U_JLOC = 1, Mdims%u_nloc
                    DO U_ILOC = 1, Mdims%u_nloc
                        DO GI = 1, FE_GIdims%cv_ngi
                            DO IPHASE = 1, Mdims%nphase
                                IF ( STRESS_FORM ) THEN ! stress form of viscosity...
                                    IF(IDIVID_BY_VOL_FRAC==1) THEN
                                        CALL CALC_STRESS_TEN( STRESS_IJ_ELE( :, :, IPHASE, U_ILOC, U_JLOC ), ZERO_OR_TWO_THIRDS, Mdims%ndim, &
                                            ( -UFEN_REVERSED( GI, U_ILOC )*VOL_FRA_GI_DX_ALL(1:Mdims%ndim,IPHASE,GI) + UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_ILOC )*VOL_FRA_GI(IPHASE,GI) ),  UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_JLOC )* DevFuns%DETWEI( GI ), TEN_XX( :, :, IPHASE, GI ), TEN_VOL( IPHASE, GI) )
                                    ELSE
                                        CALL CALC_STRESS_TEN( STRESS_IJ_ELE( :, :, IPHASE, U_ILOC, U_JLOC ), ZERO_OR_TWO_THIRDS, Mdims%ndim, &
                                            UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_ILOC ), UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_JLOC )* DevFuns%DETWEI( GI ), TEN_XX( :, :, IPHASE, GI ), TEN_VOL( IPHASE, GI) )
                                    ENDIF
                                ELSE
                                    DO IDIM = 1, Mdims%ndim
                                        VLK_ELE( IPHASE, U_ILOC, U_JLOC ) = VLK_ELE( IPHASE, U_ILOC, U_JLOC ) + &
                                            UFENX_ALL_REVERSED( IDIM, GI, U_ILOC ) * SUM( UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_JLOC ) * TEN_XX( IDIM, :, IPHASE, GI ) ) * DevFuns%DETWEI( GI )
                                    END DO
                                END IF
                            END DO
                        END DO
                    END DO
                END DO
            END IF
            if ((Porous_media_PIVIT_not_stored_yet .and..not. Mmat%CV_pressure).or..not.is_porous_media) then!sprint_to_do; Internal subroutine for this?
!            if (Porous_media_PIVIT_not_stored_yet) then!sprint_to_do; Internal subroutine for this?
                DO U_JLOC = 1, Mdims%u_nloc
                    DO U_ILOC = 1, Mdims%u_nloc
                        DO GI = 1, FE_GIdims%cv_ngi
                            RNN = UFEN_REVERSED( GI, U_ILOC ) * UFEN_REVERSED( GI, U_JLOC ) * DevFuns%DETWEI( GI )
                            if ( lump_absorption ) then
                                NN_SIGMAGI_ELE(:, :, U_ILOC, U_ILOC ) = &
                                    NN_SIGMAGI_ELE(:, :, U_ILOC, U_ILOC ) + RNN * LOC_U_ABSORB( :, :, U_ILOC)
                            else
                                NN_SIGMAGI_ELE(:, :, U_ILOC, U_JLOC ) = &
                                    NN_SIGMAGI_ELE(:, :, U_ILOC, U_JLOC ) + RNN *SIGMAGI( :, :, GI )
                            end if
                            NN_SIGMAGI_STAB_ELE(:, :, U_ILOC, U_JLOC ) = &
                                NN_SIGMAGI_STAB_ELE(:, :, U_ILOC, U_JLOC ) + RNN *SIGMAGI_STAB( :, :, GI )

                            ! Chris change ordering of NN_SIGMAGI_STAB_SOLID_RHS_ELE
                            IF(RETRIEVE_SOLID_CTY) NN_SIGMAGI_STAB_SOLID_RHS_ELE(:, :, U_ILOC, U_JLOC ) =&
                                NN_SIGMAGI_STAB_SOLID_RHS_ELE(:, :, U_ILOC, U_JLOC ) + RNN * SIGMAGI_STAB_SOLID_RHS( :, :, GI )
                            DO JPHASE = 1, Mdims%nphase
                                DO JDIM = 1, Mdims%ndim
                                    JPHA_JDIM = JDIM + (JPHASE-1)*Mdims%ndim
                                    if ( lump_mass2 ) then
                                        NN_MASS_ELE( JPHA_JDIM, JPHA_JDIM, U_ILOC, U_ILOC ) = NN_MASS_ELE( JPHA_JDIM, JPHA_JDIM, U_ILOC, U_ILOC ) &
                                            + DENGI(JPHASE,GI) * RNN
                                        NN_MASSOLD_ELE( JPHA_JDIM, JPHA_JDIM, U_ILOC, U_ILOC ) = NN_MASSOLD_ELE( JPHA_JDIM, JPHA_JDIM, U_ILOC, U_ILOC ) &
                                            + DENGIOLD(JPHASE, GI) * RNN
                                    else
                                        NN_MASS_ELE( JPHA_JDIM, JPHA_JDIM, U_ILOC, U_JLOC ) = NN_MASS_ELE( JPHA_JDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                            + DENGI(JPHASE,GI) * RNN
                                        NN_MASSOLD_ELE( JPHA_JDIM, JPHA_JDIM, U_ILOC, U_JLOC ) = NN_MASSOLD_ELE( JPHA_JDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                            + DENGIOLD(JPHASE, GI) * RNN
                                    end if
                                    IF(GOT_VIRTUAL_MASS) THEN
                                        DO IPHASE = 1, Mdims%nphase
                                            IPHA_IDIM = JDIM + (IPHASE-1)*Mdims%ndim
                                            ! Chris re-order NN_MASS_ELE & NN_MASSOLD_ELE...
                                            NN_MASS_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) = NN_MASS_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                                + VIRTUAL_MASS_GI(IPHASE,JPHASE,GI) * RNN
                                            NN_MASSOLD_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) = NN_MASSOLD_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                                + VIRTUAL_MASS_OLD_GI(IPHASE, JPHASE, GI) * RNN
                                        END DO
                                    ENDIF
                                END DO
                            END DO
                        END DO
                    END DO
                END DO
                                ! Stabilization for viscosity...
                ! Chris re-order NN_SIGMAGI_STAB_ELE...
                IF ( STAB_VISC_WITH_ABS ) THEN
                    DO U_JLOC = 1, Mdims%u_nloc
                        DO U_ILOC = 1, Mdims%u_nloc
                            DO JPHASE = 1, Mdims%nphase
                                DO JDIM = 1, Mdims%ndim
                                    JPHA_JDIM = JDIM + (JPHASE-1)*Mdims%ndim
                                    IF ( STRESS_FORM ) THEN
                                        NN_SIGMAGI_STAB_ELE( JPHA_JDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                            = NN_SIGMAGI_STAB_ELE( JPHA_JDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                            + MAX( 0.0, STRESS_IJ_ELE( JDIM, JDIM, JPHASE, U_ILOC, U_JLOC ) )
                                    ELSE
                                        NN_SIGMAGI_STAB_ELE( JPHA_JDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                            = NN_SIGMAGI_STAB_ELE( JPHA_JDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                            + MAX( 0.0, VLK_ELE( JPHASE, U_ILOC, U_JLOC ) )
                                    END IF
                                END DO
                            END DO
                        END DO
                    END DO
                END IF
                DO U_JLOC = 1, Mdims%u_nloc
                    DO U_ILOC = 1, Mdims%u_nloc
                        DO JPHASE = 1, Mdims%nphase
                            DO JDIM = 1, Mdims%ndim
                                JPHA_JDIM = JDIM + (JPHASE-1)*Mdims%ndim
                                J = JDIM+(JPHASE-1)*Mdims%ndim+(U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                DO IPHASE = 1, Mdims%nphase
                                    DO IDIM = 1, Mdims%ndim
                                        IPHA_IDIM = IDIM + (IPHASE-1)*Mdims%ndim
                                        I = IDIM+(IPHASE-1)*Mdims%ndim+(U_ILOC-1)*Mdims%ndim*Mdims%nphase
                                        !Assemble
                                        IF ( LUMP_MASS ) THEN
                                            Mmat%PIVIT_MAT( I, I, ELE ) =  Mmat%PIVIT_MAT( I, I, ELE ) + &
                                                NN_SIGMAGI_ELE(IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                                + NN_SIGMAGI_STAB_ELE(IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                                + NN_MASS_ELE(IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC )/DT
                                            if (homogenize_mass_matrix) then
                                                Mmat%PIVIT_MAT( I, I, ELE ) =  Mmat%PIVIT_MAT( I, I, ELE ) + &
                                                    lump_weight*NN_SIGMAGI_ELE(IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC )
                                                Mmat%PIVIT_MAT( I, J, ELE ) = Mmat%PIVIT_MAT( I, J, ELE ) - &
                                                    lump_weight*NN_SIGMAGI_ELE(IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC )
                                            end if
                                        ELSE
                                            Mmat%PIVIT_MAT( I, J, ELE ) =  &
                                                NN_SIGMAGI_ELE(IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                                + NN_SIGMAGI_STAB_ELE(IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                                + NN_MASS_ELE(IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC )/DT
                                        END IF
                                        IF ( .NOT.Mmat%NO_MATRIX_STORE ) THEN
                                            IF ( .NOT.JUST_BL_DIAG_MAT ) THEN!Only for inertia
                                                IF ( LUMP_MASS ) THEN
                                                    DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_ILOC, ELE ) =  &
                                                        DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_ILOC, ELE )  &
                                                        + NN_SIGMAGI_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                                        + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                                        + NN_MASS_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) / DT
                                                ELSE
                                                    DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) = &
                                                        DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE )  &
                                                        + NN_SIGMAGI_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                                        + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) &
                                                        + NN_MASS_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) / DT
                                                END IF
                                            END IF
                                        END IF
                                    END DO
                                END DO
                            END DO
                        END DO
                    END DO
                END DO
            end if ! endof if (Porous_media_PIVIT_not_stored_yet) then
            if (.not.is_porous_media) then!sprint_to_do; internal subroutine for this?
                !###LOOP Loop_DGNods1 IS NOT NECESSARY FOR POROUS MEDIA###
                Loop_DGNods1: DO U_ILOC = 1, Mdims%u_nloc
                    ! put CV source in...
                    IF ( LUMP_MASS .AND. ( Mdims%cv_nloc==6 .OR. (Mdims%cv_nloc==10 .AND. Mdims%ndim==3) ) ) THEN ! Quadratice
                        IF(Mdims%u_nloc.GE.Mdims%cv_nloc) STOP 28211 ! Code not ready yet for this.
                        Loop_CVNods21: DO U_JLOC = 1, Mdims%u_nloc
                            CV_JLOC = ELEMENT_CORNERS( U_JLOC )
                            ! Miss out the mid side nodes...
                            NM = SUM( UFEN_REVERSED( :, U_ILOC ) * UFEN_REVERSED( :, U_JLOC ) * DevFuns%DETWEI( : ) )
                            LOC_U_RHS( :, :, U_ILOC ) = LOC_U_RHS( :, :, U_ILOC ) + NM * LOC_U_SOURCE_CV( :, :, CV_JLOC )
                        END DO LOOP_CVNODS21
                    ELSE ! ENDOF IF ( LUMP_MASS .AND. ( Mdims%cv_nloc==6 .OR. (Mdims%cv_nloc==10 .AND. Mdims%ndim==3) ) ) THEN ! Quadratice
                        Loop_CVNods2: DO CV_JLOC = 1, Mdims%cv_nloc
                            IF ( RETRIEVE_SOLID_CTY .OR. FEM_BUOYANCY ) THEN
                                NM = SUM( UFEN_REVERSED( :, U_ILOC ) * CVFEN_REVERSED( :, CV_JLOC ) * DevFuns%DETWEI( : ) )
                            ELSE
                                NM = SUM( UFEN_REVERSED( :, U_ILOC ) * CVN_REVERSED( :, CV_JLOC ) * DevFuns%DETWEI( : ) )
                            END IF
                            IF ( LUMP_MASS ) THEN
                                CV_ILOC = U_ILOC
                                LOC_U_RHS( :, :, U_ILOC ) = LOC_U_RHS( :, :, U_ILOC ) + NM * LOC_U_SOURCE_CV( :, :, CV_ILOC )
                            ELSE
                                LOC_U_RHS( :, :, U_ILOC ) = LOC_U_RHS( :, :, U_ILOC ) + NM * LOC_U_SOURCE_CV( :, :, CV_JLOC )
                            END IF
                        END DO LOOP_CVNODS2
                    END IF ! ENDOF IF ( LUMP_MASS .AND. ( Mdims%cv_nloc==6 .OR. (Mdims%cv_nloc==10 .AND. Mdims%ndim==3) ) ) THEN ! Quadratice
                    Loop_DGNods2: DO U_JLOC = 1, Mdims%u_nloc
                        VLN = 0.0
                        VLN_OLD = 0.0

                        IF(GOT_VIRTUAL_MASS) THEN
                            VLN_CVM=0.0
                            VLN_OLD_CVM=0.0
                        ENDIF
                        Loop_Gauss2: DO GI = 1, FE_GIdims%CV_NGI
                            Loop_IPHASE: DO IPHASE = 1, Mdims%nphase ! Diffusion tensor
                                IF ( MOM_CONSERV ) THEN
                                    VLN( IPHASE ) = VLN( IPHASE ) - &
                                        DENGI( IPHASE, GI ) * SUM( UD( :, IPHASE, GI ) * UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_ILOC ) )  &
                                        * UFEN_REVERSED( GI, U_JLOC ) * DevFuns%DETWEI( GI ) * WITH_NONLIN
                                    VLN_OLD( IPHASE ) = VLN_OLD( IPHASE ) - &
                                        DENGI( IPHASE, GI ) * SUM( UDOLD( :, IPHASE, GI ) * UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_ILOC ) )  &
                                        * UFEN_REVERSED( GI, U_JLOC ) * DevFuns%DETWEI( GI ) * WITH_NONLIN
                                ELSE
                                    VLN( IPHASE ) = VLN( IPHASE ) + &
                                        FE_funs%ufen( U_ILOC, GI ) * DENGI( IPHASE, GI ) * SUM( UD( :, IPHASE, GI ) * UFENX_ALL_REVERSED(1:Mdims%ndim, GI, U_JLOC ) ) &
                                        * DevFuns%DETWEI( GI ) * WITH_NONLIN
                                    VLN_OLD( IPHASE ) = VLN_OLD( IPHASE ) + &
                                        FE_funs%ufen( U_ILOC, GI ) * DENGI( IPHASE, GI ) * SUM( UDOLD( :, IPHASE, GI ) * UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_JLOC ) ) &
                                        * DevFuns%DETWEI( GI ) * WITH_NONLIN
                                END IF
                                IF(GOT_VIRTUAL_MASS) THEN
                                    DO JPHASE = 1, Mdims%nphase
                                        VLN_CVM( IPHASE,JPHASE ) = VLN_CVM( IPHASE,JPHASE )  &
                                            ! conservative discretization
                                            - CVM_BETA*VIRTUAL_MASS_GI(IPHASE,JPHASE,GI)* SUM( (VIRTUAL_MASS_ADV_CUR(IPHASE,JPHASE)*UD( :, JPHASE, GI ) +(1.-VIRTUAL_MASS_ADV_CUR(IPHASE,JPHASE))*UD( :, IPHASE, GI ))* UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_ILOC ) )  &
                                            * UFEN_REVERSED( GI, U_JLOC ) * DevFuns%DETWEI( GI ) * WITH_NONLIN_CVM &
                                            ! non-conservative discretization
                                            + (1.-CVM_BETA)*UFEN_REVERSED( GI, U_ILOC ) *VIRTUAL_MASS_GI(IPHASE,JPHASE,GI)* SUM( (VIRTUAL_MASS_ADV_CUR(IPHASE,JPHASE)*UD( :, JPHASE, GI ) +(1.-VIRTUAL_MASS_ADV_CUR(IPHASE,JPHASE))*UD( :, IPHASE, GI ))* UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_JLOC ) )  &
                                            *  DevFuns%DETWEI( GI ) * WITH_NONLIN_CVM
                                        VLN_OLD_CVM( IPHASE,JPHASE ) = VLN_OLD_CVM( IPHASE,JPHASE )  &
                                            ! conservative discretization
                                            - CVM_BETA*VIRTUAL_MASS_OLD_GI(IPHASE,JPHASE,GI) * SUM( (VIRTUAL_MASS_ADV_CUR(IPHASE,JPHASE)*UDOLD( :, JPHASE, GI ) +(1.-VIRTUAL_MASS_ADV_CUR(IPHASE,JPHASE))*UDOLD( :, IPHASE, GI )) * UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_ILOC ) )  &
                                            * UFEN_REVERSED( GI, U_JLOC ) * DevFuns%DETWEI( GI ) * WITH_NONLIN_CVM &
                                            ! non-conservative discretization
                                            + (1.-CVM_BETA)*UFEN_REVERSED( GI, U_ILOC ) *VIRTUAL_MASS_OLD_GI(IPHASE,JPHASE,GI) * SUM( (VIRTUAL_MASS_ADV_CUR(IPHASE,JPHASE)*UDOLD( :, JPHASE, GI ) +(1.-VIRTUAL_MASS_ADV_CUR(IPHASE,JPHASE))*UDOLD( :, IPHASE, GI )) * UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_JLOC ) )  &
                                            *  DevFuns%DETWEI( GI ) * WITH_NONLIN_CVM
                                    END DO
                                ENDIF
                            END DO Loop_IPHASE
                        END DO Loop_Gauss2
                        NN = SUM( UFEN_REVERSED( :, U_ILOC ) * UFEN_REVERSED( :, U_JLOC ) * DevFuns%DETWEI( : ) )
                        LOC_U_RHS( :, :, U_ILOC ) =  LOC_U_RHS( :, :, U_ILOC ) + NN * LOC_U_SOURCE( :, :, U_JLOC  )
                        DO JPHASE = 1, Mdims%nphase
                            DO JDIM = 1, Mdims%ndim
                                JPHA_JDIM = (JPHASE-1)*Mdims%ndim + JDIM
                                DO IPHASE = 1, Mdims%nphase
                                    DO IDIM = 1, Mdims%ndim
                                        IPHA_IDIM = (IPHASE-1)*Mdims%ndim + IDIM
                                        IF ( MOM_CONSERV ) THEN
                                            IF ( LUMP_MASS ) THEN
                                                LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                    + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) * LOC_U( JDIM, JPHASE, U_JLOC )     &
                                                    + ( NN_MASSOLD_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) / DT ) * LOC_UOLD( JDIM, JPHASE, U_ILOC )
                                                IF(RETRIEVE_SOLID_CTY) LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                    + NN_SIGMAGI_STAB_SOLID_RHS_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) * LOC_US( JDIM, JPHASE, U_JLOC )
                                            ELSE
                                                LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                    + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) * LOC_U( JDIM, JPHASE, U_JLOC )  &
                                                    + ( NN_MASSOLD_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) / DT ) * LOC_UOLD( JDIM, JPHASE, U_JLOC )
                                                IF(RETRIEVE_SOLID_CTY) LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                    + NN_SIGMAGI_STAB_SOLID_RHS_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) * LOC_US( JDIM, JPHASE, U_JLOC )
                                            END IF
                                        ELSE
                                            IF ( LUMP_MASS ) THEN
                                                LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                    + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) * LOC_U( JDIM, JPHASE, U_JLOC )  &
                                                    + ( NN_MASS_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) / DT ) * LOC_UOLD( JDIM, JPHASE, U_ILOC )
                                                IF(RETRIEVE_SOLID_CTY) LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                    + NN_SIGMAGI_STAB_SOLID_RHS_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) * LOC_US( JDIM, JPHASE, U_JLOC )
                                            ELSE
                                                LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                    + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) * LOC_U( JDIM, JPHASE, U_JLOC )  &
                                                    + ( NN_MASS_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) / DT ) * LOC_UOLD( JDIM, JPHASE, U_JLOC )
                                                IF(RETRIEVE_SOLID_CTY) LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                    + NN_SIGMAGI_STAB_SOLID_RHS_ELE( IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC ) * LOC_US( JDIM, JPHASE, U_JLOC )
                                            END IF
                                        END IF
                                    END DO
                                END DO
                            END DO
                        END DO
                        IF ( .NOT.JUST_BL_DIAG_MAT ) THEN
                            IF ( STRESS_FORM ) THEN
                                DO IPHASE = 1, Mdims%nphase
                                    JPHASE = IPHASE
                                    DO JDIM = 1, Mdims%ndim
                                        DO IDIM = 1, Mdims%ndim
                                            IF ( Mmat%NO_MATRIX_STORE ) THEN
                                                LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                    - STRESS_IJ_ELE( IDIM, JDIM,  IPHASE, U_ILOC, U_JLOC ) * LOC_U( JDIM, IPHASE, U_JLOC )
                                            ELSE
                                                DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE )  &
                                                    = DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) &
                                                    + STRESS_IJ_ELE( IDIM, JDIM, IPHASE, U_ILOC, U_JLOC )
                                            END IF
                                            IF(PIVIT_ON_VISC) THEN
                                                I = IDIM+(IPHASE-1)*Mdims%ndim+(U_ILOC-1)*Mdims%ndim*Mdims%nphase
                                                J = JDIM+(JPHASE-1)*Mdims%ndim+(U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                                w=1.0
                                                !                                                if (i/=j) w = wv
                                                Mmat%PIVIT_MAT( I,J, ELE ) &
                                                    = Mmat%PIVIT_MAT( I,J, ELE ) &
                                                    +  w * STRESS_IJ_ELE( IDIM, JDIM, IPHASE, U_ILOC, U_JLOC )
                                            END IF
                                            RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) = RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) + &
                                                STRESS_IJ_ELE( IDIM, JDIM, IPHASE, U_ILOC, U_JLOC ) * LOC_U( JDIM, IPHASE, U_JLOC )
                                        END DO
                                    END DO
                                END DO
                            END IF
                            !! - Asiri to change nested loop order here
                            DO IDIM = 1, Mdims%ndim
                                DO IPHASE = 1, Mdims%nphase
                                    JDIM = IDIM
                                    JPHASE = IPHASE
                                    IF ( Mmat%NO_MATRIX_STORE ) THEN
                                        LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC )  &
                                            - VLN( IPHASE ) * LOC_U( IDIM, IPHASE, U_JLOC )
                                        IF(GOT_VIRTUAL_MASS) THEN
                                            DO KPHASE = 1, Mdims%nphase
                                                LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC )  &
                                                    - VLN_CVM( IPHASE,KPHASE ) * LOC_U( IDIM, KPHASE, U_JLOC )
                                            END DO
                                        ENDIF
                                    ELSE
                                        DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) &
                                            = DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) + VLN( IPHASE )
                                        IF(GOT_VIRTUAL_MASS) THEN
                                            DO KPHASE = 1, Mdims%nphase
                                                DIAG_BIGM_CON( IDIM, JDIM, IPHASE, KPHASE, U_ILOC, U_JLOC, ELE ) &
                                                    = DIAG_BIGM_CON( IDIM, JDIM, IPHASE, KPHASE, U_ILOC, U_JLOC, ELE ) + VLN_CVM( IPHASE, KPHASE )
                                            END DO
                                        ENDIF
                                    END IF
                                    IF ( .NOT.STRESS_FORM ) THEN
                                        IF ( Mmat%NO_MATRIX_STORE ) THEN
                                            LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                - VLK_ELE( IPHASE, U_ILOC, U_JLOC ) * LOC_U( IDIM, IPHASE, U_JLOC )
                                        ELSE
                                            DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) &
                                                = DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) + VLK_ELE( IPHASE, U_ILOC, U_JLOC )
                                        END IF
                                        IF(PIVIT_ON_VISC) THEN
                                            I = IDIM+(IPHASE-1)*Mdims%ndim+(U_ILOC-1)*Mdims%ndim*Mdims%nphase
                                            J = JDIM+(JPHASE-1)*Mdims%ndim+(U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                            w=1.0
                                            !                                            if (i/=j) w = wv
                                            Mmat%PIVIT_MAT( I,J, ELE ) &
                                                = Mmat%PIVIT_MAT( I,J, ELE ) + w * VLK_ELE( IPHASE, U_ILOC, U_JLOC )
                                        END IF
                                        RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) = RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) + &
                                            VLK_ELE( IPHASE, U_ILOC, U_JLOC ) * LOC_U( IDIM, IPHASE, U_JLOC )
                                    END IF
                                END DO
                            END DO
                        END IF ! .NOT.JUST_BL_DIAG_MAT
                    END DO Loop_DGNods2
                END DO Loop_DGNods1
            else !Adding sources to the RHS for porous media
                DO U_ILOC = 1, Mdims%u_nloc
                    DO CV_JLOC = 1, Mdims%cv_nloc
                        NM = SUM( UFEN_REVERSED( :, U_ILOC ) * CVN_REVERSED( :, CV_JLOC ) * DevFuns%DETWEI( : ) )
                        LOC_U_RHS( :, :, U_ILOC ) = LOC_U_RHS( :, :, U_ILOC ) + NM * LOC_U_SOURCE_CV( :, :, CV_JLOC )
                    end do
                    DO U_JLOC = 1, Mdims%u_nloc
                        NN = SUM( UFEN_REVERSED( :, U_ILOC ) * UFEN_REVERSED( :, U_JLOC ) * DevFuns%DETWEI( : ) )
                        LOC_U_RHS( :, :, U_ILOC ) =  LOC_U_RHS( :, :, U_ILOC ) + NN * LOC_U_SOURCE( :, :, U_JLOC  )
                    end do
                end do
            end if
                ! **********REVIEWER 1-END**********************
            IF(GOT_DIFFUS.AND.LINEAR_HIGHORDER_DIFFUSION) THEN
                NN_MAT_ELE( :, :, ELE ) = 0.0
                NNX_MAT_ELE( :, :, :, ELE ) = 0.0
                DO U_ILOC = 1, Mdims%u_nloc
                    DO U_JLOC = 1, Mdims%u_nloc
                        NN_MAT_ELE( U_ILOC, U_JLOC, ELE ) = NN_MAT_ELE( U_ILOC, U_JLOC, ELE ) + &
                            SUM( UFEN_REVERSED( :, U_ILOC ) * UFEN_REVERSED( :, U_JLOC ) * DevFuns%DETWEI( : ) )
                        DO IDIM=1,Mdims%ndim
                            NNX_MAT_ELE( IDIM, U_ILOC, U_JLOC, ELE ) = NNX_MAT_ELE( IDIM, U_ILOC, U_JLOC, ELE ) + &
                                SUM( UFEN_REVERSED( :, U_ILOC ) * UFENX_ALL_REVERSED( IDIM, :, U_JLOC ) * DevFuns%DETWEI( : ) )
                        END DO
                    END DO
                END DO
            ENDIF
            ! Add in Mmat%C matrix contribution: (DG velocities)
            Loop_U_ILOC1: DO U_ILOC = 1, Mdims%u_nloc
                if(.not.got_c_matrix) IU_NOD = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                Loop_P_JLOC1: DO P_JLOC = 1, Mdims%p_nloc
                    if(.not.got_c_matrix) JCV_NOD = ndgln%p( ( ELE - 1 ) * Mdims%p_nloc + P_JLOC )
                    !In this section we multiply the shape functions over the GI points. I.E: We perform the integration
                    !over the element of the pressure like source term.
                    if ( IPLIKE_GRAD_SOU >= 1 ) then
                        !In this section of the assembly we add the volumetric part.
                        ! Coeff * Integral(N grad FE_funs%cvfen PLIKE_GRAD_SOU dV)
                        DO ICOMP= 1, Mdims%ncomp
                            GRAD_SOU_GI_NMX( :, ICOMP, :) = matmul(CVFENX_ALL_REVERSED(:, :, P_JLOC),&
                                SPREAD(DevFuns%DETWEI( : ) *UFEN_REVERSED( :, U_ILOC ), DIM=2, NCOPIES=Mdims%nphase)*&
                                transpose(GRAD_SOU_GI(ICOMP, :, :)))
                            IF ( IPLIKE_GRAD_SOU == 2 ) THEN
                               GRAD_SOU2_GI_NMX( :, ICOMP, :) = matmul(CVFENX_ALL_REVERSED(:, :, P_JLOC),&
                                   SPREAD(DevFuns%DETWEI( : ) *UFEN_REVERSED( :, U_ILOC ), DIM=2, NCOPIES=Mdims%nphase)*&
                                   transpose(GRAD_SOU2_GI(ICOMP, :, :)))
                            END IF
                        END DO
                    end if
                    ! Coeff * Integral(N grad FE_funs%cvfen VOL_FRA dV)
                    IF(IGOT_VOL_X_PRESSURE==1) THEN
                        VOL_FRA_NMX_ALL( :, : ) = matmul(CVFENX_ALL_REVERSED(:, :, P_JLOC),&
                            SPREAD(DevFuns%DETWEI( : ) *UFEN_REVERSED( :, U_ILOC ), DIM=2, NCOPIES=Mdims%nphase)*&
                            transpose(VOL_FRA_GI(:, :)))
                    ENDIF
                    ! Put into matrix
                    IF ( .NOT.GOT_C_MATRIX ) THEN
                        CALL USE_POSINMAT_C_STORE( COUNT, IU_NOD, JCV_NOD,  &
                            Mdims%u_nonods, Mspars%C%fin, Mspars%C%col, Mspars%C%ncol, &
                            IDO_STORE_AC_SPAR_PT, STORED_AC_SPAR_PT, POSINMAT_C_STORE, ELE, U_ILOC, P_JLOC, &
                            Mdims%totele, Mdims%u_nloc, Mdims%p_nloc )
                    END IF
                    !Prepare aid variable NMX_ALL to improve the speed of the calculations
                    NMX_ALL( : ) = matmul(CVFENX_ALL_REVERSED(:,:,P_JLOC),  DevFuns%DETWEI( : ) *UFEN_REVERSED( :, U_ILOC ))
                    Loop_Phase1: DO IPHASE = 1, Mdims%nphase
                        ! Put into matrix
                        IF ( .NOT.GOT_C_MATRIX ) THEN
                            IF(IGOT_VOL_X_PRESSURE==1) THEN
                                Mmat%C( :, IPHASE, COUNT ) = Mmat%C( :, IPHASE, COUNT ) - VOL_FRA_NMX_ALL( :, IPHASE )
                            ELSE
                                DO IDIM = 1, Mdims%ndim
                                    Mmat%C( IDIM, IPHASE, COUNT ) = Mmat%C( IDIM, IPHASE, COUNT ) - NMX_ALL( IDIM )
                                END DO
                            ENDIF
                        END IF
                        IF ( IPLIKE_GRAD_SOU >= 1 ) THEN ! Pressure like term
                            DO IDIM = 1, Mdims%ndim
                                DO ICOMP = 1, Mdims%ncomp
                                    LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                        - GRAD_SOU_GI_NMX( IDIM, ICOMP, IPHASE ) * LOC_PLIKE_GRAD_SOU_GRAD( ICOMP, IPHASE, P_JLOC )
                                    IF ( IPLIKE_GRAD_SOU == 2 ) THEN
                                       LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                           - GRAD_SOU2_GI_NMX( IDIM, ICOMP, IPHASE ) * LOC_PLIKE_GRAD_SOU_COEF( ICOMP, IPHASE, P_JLOC )
                                    END IF
                                END DO
                            END DO
                        END IF
                    END DO Loop_Phase1
                END DO Loop_P_JLOC1
            END DO Loop_U_ILOC1
!            IF( (.NOT.first_time_step) .AND. (RESID_BASED_STAB_DIF/=0) ) THEN!sprint_to_do; first_time_step had to be always true
            if (RESID_BASED_STAB_DIF/=0) then
                !! *************************INNER ELEMENT STABILIZATION****************************************
                !! *************************INNER ELEMENT STABILIZATION****************************************
                DO U_JLOC = 1, Mdims%u_nloc
                    DO U_ILOC = 1, Mdims%u_nloc
                        ! Sum over quadrature pts...
                        LOC_MASS( U_ILOC, U_JLOC ) = SUM( UFEN_REVERSED( :, U_ILOC ) * UFEN_REVERSED( :, U_JLOC ) * DevFuns%DETWEI( : ) )
                    END DO
                END DO
                LOC_MASS_INV = LOC_MASS
                !CALL INVERT(LOC_MASS_INV)
                CALL MATDMATINV( LOC_MASS, LOC_MASS_INV, Mdims%u_nloc )
                DO U_ILOC = 1, Mdims%u_nloc
                    DO IPHASE = 1, Mdims%nphase
                        DO IDIM = 1, Mdims%ndim
                            ! sum cols of matrix * rows of vector...
                            DIFF_VEC_U( IDIM, IPHASE, U_ILOC ) = SUM( LOC_MASS_INV( U_ILOC, : ) * RHS_DIFF_U( IDIM, IPHASE, : ) )
                        END DO
                    END DO
                END DO
                DIFFGI_U = 0.0
                U_DX_ALL = 0.0 ; UOLD_DX_ALL = 0.0
                SOUGI_X = 0.0
                DO U_ILOC = 1, Mdims%u_nloc
                    DO GI = 1, FE_GIdims%cv_ngi
                        DO IPHASE = 1, Mdims%nphase
                            DIFFGI_U( :, IPHASE, GI ) = DIFFGI_U( :, IPHASE, GI ) + &
                                UFEN_REVERSED( GI, U_ILOC ) * DIFF_VEC_U( :, IPHASE, U_ILOC )
                            SOUGI_X( :, IPHASE, GI ) = SOUGI_X( :, IPHASE, GI ) + &
                                UFEN_REVERSED( GI, U_ILOC ) * LOC_U_SOURCE( :, IPHASE, U_ILOC )
                            DO JDIM = 1, Mdims%ndim
                                DO IDIM = 1, Mdims%ndim
                                    U_DX_ALL( IDIM, JDIM, IPHASE, GI ) = U_DX_ALL( IDIM, JDIM, IPHASE, GI ) + &
                                        LOC_U( JDIM, IPHASE, U_ILOC ) * UFENX_ALL_REVERSED( IDIM, GI, U_ILOC )
                                    UOLD_DX_ALL( IDIM, JDIM, IPHASE, GI ) = UOLD_DX_ALL( IDIM, JDIM, IPHASE, GI ) + &
                                        LOC_UOLD( JDIM, IPHASE, U_ILOC ) * UFENX_ALL_REVERSED( IDIM, GI, U_ILOC )
                                END DO
                            END DO
                        END DO
                    END DO
                END DO
                U_DT = ( UD - UDOLD ) / DT
                RESID_U = 0.0
                DO GI = 1, FE_GIdims%cv_ngi
                    DO IPHASE = 1, Mdims%nphase
                        DO IDIM = 1, Mdims%ndim
                            IPHA_IDIM = (IPHASE-1)*Mdims%ndim + IDIM
                            DO JPHASE = 1, Mdims%nphase
                                DO JDIM = 1, Mdims%ndim
                                    JPHA_JDIM = (JPHASE-1)*Mdims%ndim + JDIM
                                    RESID_U( IDIM, IPHASE, GI ) = RESID_U( IDIM, IPHASE, GI ) + &
                                        SIGMAGI( IPHA_IDIM, JPHA_JDIM, GI ) * UD( JDIM, IPHASE, GI )
                                END DO
                            END DO
                        END DO
                    END DO
                END DO
                P_DX = 0.0
                RR = 0.0
                DO P_ILOC = 1, Mdims%p_nloc
                    DO GI = 1, FE_GIdims%cv_ngi
                        P_DX( :, GI ) = P_DX( :, GI ) + CVFENX_ALL_REVERSED(1:Mdims%ndim, GI, P_ILOC ) * LOC_P( P_ILOC )
                        IF ( IPLIKE_GRAD_SOU >= 1 ) THEN ! Pressure like terms...
                            DO IPHASE = 1, Mdims%nphase
                                DO ICOMP = 1, Mdims%ncomp
                                    R = GRAD_SOU_GI( ICOMP, IPHASE, GI ) * LOC_PLIKE_GRAD_SOU_GRAD( ICOMP, IPHASE, P_ILOC )
                                    IF ( IPLIKE_GRAD_SOU == 2 ) & ! Other part of solid pressure
                                        RR = GRAD_SOU2_GI( ICOMP, IPHASE, GI ) * LOC_PLIKE_GRAD_SOU_COEF( ICOMP, IPHASE, P_ILOC )
                                    DO IDIM = 1, Mdims%ndim
                                        RESID_U( IDIM, IPHASE, GI ) = RESID_U( IDIM, IPHASE, GI ) + (R + RR) * CVFENX_ALL_REVERSED( IDIM, GI, P_ILOC )
                                    END DO
                                END DO
                            END DO
                        END IF
                    END DO
                END DO
                DO GI = 1, FE_GIdims%cv_ngi
                    DO IPHASE = 1, Mdims%nphase
                        DO IDIM = 1, Mdims%ndim
                            RESID_U( IDIM, IPHASE, GI ) = RESID_U( IDIM, IPHASE, GI ) + &
                                DENGI( IPHASE, GI ) * SUM( UD_ND( :, IPHASE, GI ) * U_DX_ALL( :, IDIM, IPHASE, GI ) ) &
                                * WITH_NONLIN &
                                + DENGI( IPHASE, GI ) * U_DT( IDIM, IPHASE, GI )   &
                                - SOUGI_X( IDIM, IPHASE, GI ) - DIFFGI_U( IDIM, IPHASE, GI ) + P_DX( IDIM, GI )
                            IF(GOT_VIRTUAL_MASS) THEN
                                DO JPHASE=1,Mdims%nphase
                                    RESID_U( IDIM, IPHASE, GI ) = RESID_U( IDIM, IPHASE, GI ) + &
                                        VIRTUAL_MASS_GI(IPHASE,JPHASE,GI) * SUM( (VIRTUAL_MASS_ADV_CUR(IPHASE,JPHASE)*UD( :, JPHASE, GI ) +(1.-VIRTUAL_MASS_ADV_CUR(IPHASE,JPHASE))*UD( :, IPHASE, GI )) * U_DX_ALL( :, IDIM, JPHASE, GI ) ) &
                                        * WITH_NONLIN_CVM
                                END DO
                            ENDIF
                            U_GRAD_NORM2( IDIM, IPHASE, GI ) = U_DT( IDIM, IPHASE, GI )**2 + SUM( U_DX_ALL( :, IDIM, IPHASE, GI )**2 )
                            U_GRAD_NORM( IDIM, IPHASE, GI ) = MAX( TOLER, SQRT( U_GRAD_NORM2( IDIM, IPHASE, GI ) ) )
                            U_GRAD_NORM2( IDIM, IPHASE, GI ) = MAX( TOLER, U_GRAD_NORM2( IDIM, IPHASE, GI ) )
                            A_DOT_U( IDIM, IPHASE, GI ) = DENGI( IPHASE, GI ) * ( SUM( UD_ND( :, IPHASE, GI ) * U_DX_ALL( :, IDIM, IPHASE, GI ) ) &
                                * WITH_NONLIN + U_DT( IDIM, IPHASE, GI ) ) + P_DX( IDIM, GI ) * RNO_P_IN_A_DOT
                            STAR_U_COEF( IDIM, IPHASE, GI ) = A_DOT_U( IDIM, IPHASE, GI ) / U_GRAD_NORM2( IDIM, IPHASE, GI )
                        END DO
                        JTT_INV = 2. / DT
                        U_GRAD_N_MAX2=0.0
                        DO U_ILOC = 1, Mdims%u_nloc
                            DO IDIM = 1, Mdims%ndim
                                U_GRAD_N_MAX2( IDIM ) = MAX( U_GRAD_N_MAX2( IDIM ), &
                                    ( JTT_INV * U_DT( IDIM, IPHASE, GI ) )**2 &
                                    + 4. * SUM( ( DevFuns%UFENX_ALL( 1:Mdims%ndim, U_ILOC, GI ) * U_DX_ALL( 1:Mdims%ndim, IDIM, IPHASE, GI ) )**2 ) )
                            END DO
                        END DO
                        DO IDIM = 1, Mdims%ndim
                            P_STAR_U( IDIM, IPHASE, GI ) = U_NONLIN_SHOCK_COEF / MAX( TOLER, SQRT( STAR_U_COEF( IDIM, IPHASE, GI )**2 * U_GRAD_N_MAX2( IDIM ) ) )
                        END DO
                        IF ( RESID_BASED_STAB_DIF==1 ) THEN
                            U_R2_COEF( : ) = RESID_U( :, IPHASE, GI )**2
                        ELSE IF ( RESID_BASED_STAB_DIF==2 ) THEN
                            U_R2_COEF( : ) = MAX( 0.0, A_DOT_U( :, IPHASE, GI ) * RESID_U( :, IPHASE, GI ) )
                        ELSE IF ( RESID_BASED_STAB_DIF==3 ) THEN ! Max of two previous methods.
                            U_R2_COEF( : ) = MAX( RESID_U( :, IPHASE, GI )**2, A_DOT_U( :, IPHASE, GI ) * RESID_U( :, IPHASE, GI ) )
                        END IF
                        DIF_STAB_U( :, IPHASE, GI ) = U_R2_COEF( : ) * P_STAR_U( :, IPHASE, GI ) / U_GRAD_NORM2( :, IPHASE, GI )
                    END DO
                END DO
                DIF_STAB_U = MAX(0.0, DIF_STAB_U)
                IF ( STRESS_FORM_STAB ) THEN! stress form of viscosity...
                    DO IDIM=1,Mdims%ndim
                        DO JDIM=1,Mdims%ndim
                            TEN_XX( IDIM, JDIM, :, : ) = SQRT( DIF_STAB_U( IDIM, :, : )  *  DIF_STAB_U( JDIM, :, : ) )
                        END DO
                    END DO
                    TEN_VOL=0.0
                    STRESS_IJ_ELE=0.0
                    DO U_JLOC = 1, Mdims%u_nloc
                        DO U_ILOC = 1, Mdims%u_nloc
                            DO GI = 1, FE_GIdims%cv_ngi
                                DO IPHASE = 1, Mdims%nphase
                                    IF(IDIVID_BY_VOL_FRAC==1) THEN
                                        CALL CALC_STRESS_TEN( STRESS_IJ_ELE( :, :, IPHASE, U_ILOC, U_JLOC ), ZERO_OR_TWO_THIRDS, Mdims%ndim, &
                                            ( -UFEN_REVERSED( GI, U_ILOC )*VOL_FRA_GI_DX_ALL(1:Mdims%ndim,IPHASE,GI) + UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_ILOC )*VOL_FRA_GI(IPHASE,GI) ),  UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_JLOC )* DevFuns%DETWEI( GI ), TEN_XX( :, :, IPHASE, GI ), TEN_VOL(IPHASE,GI) )
                                    ELSE
                                        CALL CALC_STRESS_TEN( STRESS_IJ_ELE( :, :, IPHASE, U_ILOC, U_JLOC ), ZERO_OR_TWO_THIRDS, Mdims%ndim, &
                                            UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_ILOC ), UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_JLOC )* DevFuns%DETWEI( GI ), TEN_XX( :, :, IPHASE, GI ), TEN_VOL(IPHASE,GI) )
                                    ENDIF
                                END DO
                            END DO
                            DO IPHASE = 1, Mdims%nphase
                                JPHASE = IPHASE
                                DO IDIM = 1, Mdims%ndim
                                    DO JDIM = 1, Mdims%ndim
                                        DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE )  &
                                            = DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) &
                                            + STRESS_IJ_ELE( IDIM, JDIM, IPHASE, U_ILOC, U_JLOC )
                                        IF(PIVIT_ON_VISC) THEN
                                            I = IDIM+(IPHASE-1)*Mdims%ndim+(U_ILOC-1)*Mdims%ndim*Mdims%nphase
                                            J = JDIM+(JPHASE-1)*Mdims%ndim+(U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                            !                                        w=1.0
                                            !                                        if (i/=j) w = wv
                                            !                                        Mmat%PIVIT_MAT( I,J, ELE )  &
                                            !                                        = Mmat%PIVIT_MAT( I,J, ELE ) + w * STRESS_IJ_ELE( IDIM, JDIM, IPHASE, U_ILOC, U_JLOC )
                                            Mmat%PIVIT_MAT( I,J, ELE )  &
                                                = Mmat%PIVIT_MAT( I,J, ELE ) + STRESS_IJ_ELE( IDIM, JDIM, IPHASE, U_ILOC, U_JLOC )
                                        ENDIF
                                    END DO
                                END DO
                            END DO
                        END DO
                    END DO
                ELSE ! endof IF ( STRESS_FORM_STAB ) THEN ELSE - stress form of viscosity...
                    ! Place the diffusion term into matrix...
                    DO U_JLOC = 1, Mdims%u_nloc
                        DO U_ILOC = 1, Mdims%u_nloc
                            DO IPHASE = 1, Mdims%nphase
                                JPHASE = IPHASE
                                VLK_UVW = 0.0
                                DO GI = 1, FE_GIdims%cv_ngi
                                    VLKNN = SUM( UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_ILOC ) * UFENX_ALL_REVERSED( 1:Mdims%ndim, GI, U_JLOC ) ) * DevFuns%DETWEI( GI )
                                    VLK_UVW( : ) = VLK_UVW( : ) + DIF_STAB_U( :, IPHASE, GI ) * VLKNN
                                END DO
                                DO IDIM = 1, Mdims%ndim
                                    JDIM = IDIM
                                    IF ( .NOT.JUST_BL_DIAG_MAT ) THEN
                                        IF ( Mmat%NO_MATRIX_STORE ) THEN
                                            LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                - VLK_UVW( IDIM ) * LOC_U( IDIM, IPHASE, U_JLOC )
                                        ELSE
                                            DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE )  &
                                                = DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) + VLK_UVW( IDIM )
                                        END IF
                                        IF(PIVIT_ON_VISC) THEN
                                            I = IDIM+(IPHASE-1)*Mdims%ndim+(U_ILOC-1)*Mdims%ndim*Mdims%nphase
                                            J = JDIM+(JPHASE-1)*Mdims%ndim+(U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                            !                                        w=1.0
                                            !                                        if (i/=j) w = wv
                                            !                                        Mmat%PIVIT_MAT( I,J, ELE )  &
                                            !                                        = Mmat%PIVIT_MAT( I,J, ELE ) + w * VLK_UVW( IDIM )
                                            Mmat%PIVIT_MAT( I,J, ELE ) = Mmat%PIVIT_MAT( I,J, ELE ) + VLK_UVW( IDIM )
                                        ENDIF
                                    END IF
                                END DO
                            END DO
                        END DO
                    END DO
                ENDIF  ! endof IF ( STRESS_FORM_STAB ) THEN ELSE
                ! Place the diffusion term into matrix for between element diffusion stabilization...
                IF ( BETWEEN_ELE_STAB ) THEN
                    ! we store these vectors in order to try and work out the between element
                    ! diffusion/viscocity.
                    DO U_JLOC = 1, Mdims%u_nloc
                        DO U_ILOC = 1, Mdims%u_nloc
                            MAT_ELE( U_ILOC, U_JLOC, ELE ) = MAT_ELE( U_ILOC, U_JLOC, ELE ) + &
                                SUM( UFEN_REVERSED( :, U_ILOC ) * UFEN_REVERSED( :, U_JLOC ) * DevFuns%DETWEI( : ) )
                        END DO
                    END DO
                    DO U_ILOC = 1, Mdims%u_nloc
                        DO IPHASE = 1, Mdims%nphase
                            DO GI = 1, FE_GIdims%cv_ngi
                                ! we store these vectors in order to try and work out the between element
                                ! diffusion/viscocity.
                                DIFF_FOR_BETWEEN_U( :, IPHASE, U_ILOC, ELE ) = DIFF_FOR_BETWEEN_U( :, IPHASE, U_ILOC, ELE ) &
                                    + UFEN_REVERSED( GI, U_ILOC ) * DevFuns%DETWEI( GI ) * DIF_STAB_U( :, IPHASE, GI )
                            END DO
                        END DO
                    END DO
                   ! End of IF(BETWEEN_ELE_STAB) THEN...
                END IF
                ! Place the diffusion term into matrix for between element diffusion stabilization...
                IF ( THERMAL_STAB_VISC ) THEN
                    ! we store these vectors in order to try and work out the between element
                    ! diffusion/viscocity.
                    DO CV_ILOC = 1, Mdims%cv_nloc
                        DO CV_JLOC = 1, Mdims%cv_nloc
                            MAT_ELE_CV_LOC( CV_ILOC, CV_JLOC ) = &
                                !SUM( CVFEN_REVERSED( :, CV_ILOC ) * CVFEN_REVERSED( :, CV_JLOC ) * DevFuns%DETWEI( : ) )
                                SUM( CVN_REVERSED( :, CV_ILOC ) * CVFEN_REVERSED( :, CV_JLOC ) * DevFuns%DETWEI( : ) )
                        END DO
                    END DO
                    DIFF_FOR_BETWEEN_CV = 0.0
                    DO CV_ILOC = 1, Mdims%cv_nloc
                        DO IPHASE = 1, Mdims%nphase
                            DO GI = 1, FE_GIdims%cv_ngi
                                ! we store these vectors in order to try and work out the between element
                                ! diffusion/viscocity.
                                DIFF_FOR_BETWEEN_CV( :, IPHASE, CV_ILOC ) = DIFF_FOR_BETWEEN_CV( :, IPHASE, CV_ILOC ) &
                                    !+ CVFEN_REVERSED( GI, CV_ILOC ) * DevFuns%DETWEI( GI ) * DIF_STAB_U( :, IPHASE, GI )
                                    + CVN_REVERSED( GI, CV_ILOC ) * DevFuns%DETWEI( GI ) * DIF_STAB_U( :, IPHASE, GI )
                            END DO
                        END DO
                    END DO
                    INV_MAT_ELE_CV_LOC = INVERSE( MAT_ELE_CV_LOC )
                    DO IPHASE = 1, Mdims%nphase
                        DO IDIM=1,Mdims%ndim
                            DIFFCV( IDIM,IPHASE,: ) = MATMUL( INV_MAT_ELE_CV_LOC(:,:), DIFF_FOR_BETWEEN_CV( IDIM,IPHASE,:) )
                        END DO
                    END DO
                    DIFFCV = MAX( 0.0, DIFFCV )
                    DO IDIM=1,Mdims%ndim
                        DO JDIM=1,Mdims%ndim
                            DIFFCV_TEN(IDIM,JDIM, :,:) = SQRT( DIFFCV(IDIM, :,:) * DIFFCV(JDIM, :,:) )
                        END DO
                    END DO
                    DIFFCV_TEN_ELE(:,:, :,:,ELE) = DIFFCV_TEN
                   ! End of IF( THERMAL_STAB_VISC ) THEN...
                END IF
              ! endof IF(RESID_BASED_STAB_DIF.NE.0) THEN
            END IF
            ! copy local memory
            DO U_ILOC = 1, Mdims%u_nloc
                U_INOD = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                if (.not. node_owned(velocity,u_inod)) cycle
                Mmat%U_RHS( :, :, U_INOD ) = Mmat%U_RHS( :, :, U_INOD ) + LOC_U_RHS( :, :, U_ILOC )
            END DO
        END DO Loop_Elements
        IF ( Q_SCHEME ) THEN
            THERM_U_DIFFUSION = 0.0
            THERM_U_DIFFUSION_VOL = 0.0
            IF ( THERMAL_STAB_VISC ) THEN ! Petrov-Galerkin visc...
                RCOUNT_NODS = 0.0
                DO ELE = 1, Mdims%totele
                    DO MAT_ILOC=1,Mdims%mat_nloc
                        MAT_NOD = ndgln%mat( (ELE-1)*Mdims%mat_nloc + MAT_ILOC )
                        THERM_U_DIFFUSION( :,:,:,MAT_NOD ) = THERM_U_DIFFUSION( :,:,:,MAT_NOD ) + DIFFCV_TEN_ELE( :,:,:,MAT_ILOC,ELE ) * MASS_ELE( ELE )
                        RCOUNT_NODS( MAT_NOD ) = RCOUNT_NODS( MAT_NOD ) + MASS_ELE( ELE )
                    END DO
                END DO
                DO MAT_NOD = 1, Mdims%mat_nonods
                    THERM_U_DIFFUSION( :,:,:,MAT_NOD ) = THERM_U_DIFFUSION( :,:,:,MAT_NOD ) / RCOUNT_NODS( MAT_NOD )
                    THERM_U_DIFFUSION_VOL( :,MAT_NOD ) = 0.0
                END DO
            END IF
            ! Put the fluid viscocity (also includes LES viscocity) into the Q-scheme thermal viscocity
            IF ( THERMAL_FLUID_VISC .AND. THERMAL_LES_VISC) THEN
                THERM_U_DIFFUSION = THERM_U_DIFFUSION + UDIFFUSION_ALL
                THERM_U_DIFFUSION_VOL = THERM_U_DIFFUSION_VOL + UDIFFUSION_VOL_ALL
            ELSE IF ( THERMAL_FLUID_VISC ) THEN
                THERM_U_DIFFUSION = THERM_U_DIFFUSION + UDIFFUSION
                !THERM_U_DIFFUSION_VOL = THERM_U_DIFFUSION_VOL + UDIFFUSION_VOL
                if ( UDIFFUSION_VOL%have_field ) THERM_U_DIFFUSION_VOL = THERM_U_DIFFUSION_VOL + UDIFFUSION_VOL%val(:,1,1,:)
            END IF
        END IF
        !!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!!
        !!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!!
        !! *************************loop over surfaces*********************************************
        ! at some pt we need to merge these 2 loops but there is a bug when doing that!!!!!
        !it does not work because some things, like MASS_ELE(ele2) require to have gone through all the elements first
        ! **********REVIEWER 3-START**********************
        DISC_PRES = ( Mdims%cv_nonods == Mdims%totele * Mdims%cv_nloc )
        Loop_Elements2: DO ELE = 1, Mdims%totele
            if (IsParallel()) then
                if (.not. assemble_ele(pressure,ele)) then
                    skip=.true.
                    neighbours=>ele_neigh(pressure,ele)
                    do nb=1,size(neighbours)
                        if (neighbours(nb)<=0) cycle
                        if (assemble_ele(pressure,neighbours(nb))) then
                            skip=.false.
                            exit
                        end if
                    end do
                    if (skip) cycle
                end if
            end if
            ! for copy local memory copying...
            LOC_U_RHS = 0.0
            Between_Elements_And_Boundary: DO IFACE = 1, FE_GIdims%nface
                ELE2  = FACE_ELE( IFACE, ELE )
                SELE2 = MAX( 0, - ELE2 )
                SELE  = SELE2
                ELE2  = MAX( 0, + ELE2 )
                ! Find COUNT_ELE
                IF(.NOT.Mmat%NO_MATRIX_STORE) THEN
                    DO COUNT=Mspars%ELE%fin(ELE), Mspars%ELE%fin(ELE+1)-1
                        IF(ELE2==Mspars%ELE%col(COUNT)) COUNT_ELE=COUNT
                    END DO
                ENDIF
                ! The surface nodes on element face IFACE.
                U_SLOC2LOC( : ) = FE_funs%u_sloclist( IFACE, : )
                CV_SLOC2LOC( : ) = FE_funs%cv_sloclist( IFACE, : )
                ! Create local copy of X_ALL
                DO CV_ILOC = 1, Mdims%cv_nloc
                    X_INOD = ndgln%x( (ELE-1)*Mdims%x_nloc + CV_ILOC )
                    XL_ALL(:,CV_ILOC) = X_ALL( :, X_INOD )
                END DO
                ! Create local copy of X_ALL for surface nodes
                DO CV_SILOC = 1, Mdims%cv_snloc
                    CV_ILOC = CV_SLOC2LOC( CV_SILOC )
                    X_INOD = ndgln%x( (ELE-1)*Mdims%x_nloc + CV_ILOC )
                    XSL_ALL( :, CV_SILOC ) = X_ALL( :, X_INOD )
                END DO
                CALL DGSIMPLNORM_ALL( Mdims%cv_nloc, Mdims%cv_snloc, Mdims%ndim, &
                    XL_ALL, XSL_ALL, NORMX_ALL )
                CALL DGSDETNXLOC2_ALL( Mdims%cv_snloc, FE_GIdims%sbcvngi, Mdims%ndim, &
                    XSL_ALL, &
                    FE_funs%sbcvfen, FE_funs%sbcvfenslx, FE_funs%sbcvfensly, FE_funs%sbcvfeweigh, SDETWE, SAREA, &
                    SNORMXN_ALL, &
                    NORMX_ALL )
                !Surface integral along an element
                If_ele2_notzero_1: IF(ELE2 /= 0) THEN
                    ! ***********SUBROUTINE DETERMINE_OTHER_SIDE_FACE - START************
                    If_stored: IF(STORED_OTHER_SIDE) THEN
                        U_ILOC_OTHER_SIDE( : ) = STORED_U_ILOC_OTHER_SIDE( :, IFACE, ELE )
                        U_OTHER_LOC( : )       = STORED_U_OTHER_LOC( :, IFACE, ELE )
                        MAT_OTHER_LOC( : )     = STORED_MAT_OTHER_LOC( :, IFACE, ELE )
                    ELSE If_stored
                        U_OTHER_LOC=0
                        U_ILOC_OTHER_SIDE=0
                        IF( Mdims%xu_nloc == 1 ) THEN ! For constant vel basis functions...
                            U_ILOC_OTHER_SIDE( 1 ) = 1
                            U_OTHER_LOC( 1 )= 1
                        ELSE
                            DO U_SILOC = 1, Mdims%u_snloc
                                U_ILOC = U_SLOC2LOC( U_SILOC )
                                U_INOD = ndgln%xu( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                                DO U_ILOC2 = 1, Mdims%u_nloc
                                    U_INOD2 = ndgln%xu(( ELE2 - 1 ) * Mdims%u_nloc + U_ILOC2 )
                                    IF ( U_INOD2 == U_INOD ) THEN
                                        U_ILOC_OTHER_SIDE( U_SILOC ) = U_ILOC2
                                        U_OTHER_LOC( U_ILOC )=U_ILOC2
                                        exit
                                    END IF
                                END DO
                            END DO
                        ENDIF
                        MAT_OTHER_LOC=0
                        DO MAT_SILOC = 1, Mdims%cv_snloc
                            MAT_ILOC = CV_SLOC2LOC( MAT_SILOC )
                            MAT_INOD = ndgln%x(( ELE - 1 ) * Mdims%mat_nloc + MAT_ILOC )
                            DO MAT_ILOC2 = 1, Mdims%mat_nloc
                                MAT_INOD2 = ndgln%x(( ELE2 - 1 ) * Mdims%mat_nloc + MAT_ILOC2 )
                                IF ( MAT_INOD2 == MAT_INOD ) THEN
                                    MAT_OTHER_LOC( MAT_ILOC )=MAT_ILOC2
                                    exit
                                END IF
                            END DO
                        END DO
                        IF ( ISTORED_OTHER_SIDE.NE.0 ) THEN
                            STORED_U_ILOC_OTHER_SIDE( :, IFACE, ELE ) = U_ILOC_OTHER_SIDE( : )
                            STORED_U_OTHER_LOC( :, IFACE, ELE )       = U_OTHER_LOC( : )
                            STORED_MAT_OTHER_LOC( :, IFACE, ELE )     = MAT_OTHER_LOC( : )
                        END IF
                    END IF If_stored
                   ! ***********SUBROUTINE DETERMINE_OTHER_SIDE_FACE - END************
                END IF If_ele2_notzero_1
                ! insert calculation method for matrix ready for determining high order linear method for viscocity...
                IF(GOT_DIFFUS .AND. LINEAR_HIGHORDER_DIFFUSION) THEN
                    IF(ELE2.GT.0) THEN
                        CALL DG_VISC_LIN( S_INV_NNX_MAT12, NNX_MAT_ELE(:,:,:,ELE), NNX_MAT_ELE(:,:,:,ELE2), NN_MAT_ELE(:,:,ELE), NN_MAT_ELE(:,:,ELE2),  &
                            Mdims%u_snloc, Mdims%u_nloc, SBUFEN_REVERSED, SDETWE, FE_GIdims%sbcvngi, SNORMXN_ALL, Mdims%ndim, &
                            U_SLOC2LOC, U_OTHER_LOC, 2*Mdims%u_nloc , .FALSE. )
                    ELSE
                        CALL DG_VISC_LIN( S_INV_NNX_MAT12, NNX_MAT_ELE(:,:,:,ELE), NNX_MAT_ELE(:,:,:,ELE), NN_MAT_ELE(:,:,ELE), NN_MAT_ELE(:,:,ELE),  &
                            Mdims%u_snloc, Mdims%u_nloc, SBUFEN_REVERSED, SDETWE, FE_GIdims%sbcvngi, SNORMXN_ALL, Mdims%ndim, &
                            U_SLOC2LOC, U_OTHER_LOC, 2*Mdims%u_nloc , .TRUE.)
                    ENDIF
                ENDIF
                !Calculate all the necessary stuff and introduce the CapPressure in the RHS
                if (capillary_pressure_activated.and..not. Diffusive_cap_only) call Introduce_Cap_press_term(&
                    packed_state, Mdims, FE_funs, Devfuns, X_ALL, LOC_U_RHS, ele, &
                    ndgln%cv, ndgln%x, ele2, iface,&
                    sdetwe, SNORMXN_ALL, U_SLOC2LOC, CV_SLOC2LOC, MAT_OTHER_LOC )
                ! ********Mapping to local variables****************
                ! CV variables...
                DO CV_SILOC = 1, Mdims%cv_snloc
                    CV_ILOC = CV_SLOC2LOC( CV_SILOC )
                    CV_INOD = ndgln%cv( (ELE-1)*Mdims%cv_nloc + CV_ILOC )
                    MAT_INOD = ndgln%mat( (ELE-1)*Mdims%cv_nloc + CV_ILOC )
                    IF ( ELE2 /= 0) THEN
                        CV_ILOC2 = MAT_OTHER_LOC( CV_ILOC )
                        CV_INOD2 = ndgln%cv( (ELE2-1)*Mdims%cv_nloc + CV_ILOC2 )
                        MAT_INOD2 = ndgln%mat( (ELE2-1)*Mdims%cv_nloc + CV_ILOC2 )
                    ELSE
                        CV_ILOC2 = CV_ILOC
                        CV_INOD2 = CV_INOD
                        MAT_INOD2 = MAT_INOD
                    END IF
                    IF(IGOT_VOL_X_PRESSURE==1) THEN
                        SLOC_UDEN( :, CV_SILOC )  = UDEN( :, CV_INOD ) * FEM_VOL_FRAC( :, CV_INOD )
                        SLOC2_UDEN( :, CV_SILOC ) = UDEN( :, CV_INOD2 ) * FEM_VOL_FRAC( :, CV_INOD2 )
                        SLOC_UDENOLD( :, CV_SILOC ) = UDENOLD( :, CV_INOD ) * FEM_VOL_FRAC( :, CV_INOD )
                        SLOC2_UDENOLD( :, CV_SILOC ) = UDENOLD( :, CV_INOD2 ) * FEM_VOL_FRAC( :, CV_INOD2 )
                    ELSE
                        SLOC_UDEN( :, CV_SILOC )  = UDEN( :, CV_INOD )
                        SLOC2_UDEN( :, CV_SILOC ) = UDEN( :, CV_INOD2 )
                        SLOC_UDENOLD( :, CV_SILOC ) = UDENOLD( :, CV_INOD )
                        SLOC2_UDENOLD( :, CV_SILOC ) = UDENOLD( :, CV_INOD2 )
                    ENDIF
                    if (is_flooding) then
                        SLOC_UDEN( :, CV_SILOC )  = 1.0
                        SLOC2_UDEN( :, CV_SILOC ) = 1.0
                        SLOC_UDENOLD( :, CV_SILOC ) = 1.0
                        SLOC2_UDENOLD( :, CV_SILOC ) = 1.0
                    end if
                    IF(GOT_VIRTUAL_MASS) THEN
                        SLOC_VIRTUAL_MASS( :,:, CV_SILOC )   = VIRTUAL_MASS( :,:, CV_INOD )
                        SLOC2_VIRTUAL_MASS( :,:, CV_SILOC )  = VIRTUAL_MASS( :,:, CV_INOD2 )
                        SLOC_VIRTUAL_MASS_OLD( :,:, CV_SILOC )   = VIRTUAL_MASS_OLD( :,:, CV_INOD )
                        SLOC2_VIRTUAL_MASS_OLD( :,:, CV_SILOC )  = VIRTUAL_MASS_OLD( :,:, CV_INOD2 )
                    ENDIF
                    IF(IDIVID_BY_VOL_FRAC+IGOT_VOL_X_PRESSURE.GE.1) THEN
                        SLOC_VOL_FRA( :, CV_SILOC )  = FEM_VOL_FRAC( :, CV_INOD )
                        SLOC2_VOL_FRA( :, CV_SILOC ) = FEM_VOL_FRAC( :, CV_INOD2 )
                    ENDIF
                    IF ( GOT_DIFFUS ) THEN
                        DO IPHASE = 1, Mdims%nphase
                            IF ( ELE2 /= 0) THEN ! Only put LES visc. if not on the boundary of the domain...
                                SLOC_UDIFFUSION( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, CV_SILOC ) = UDIFFUSION_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, MAT_INOD )
                                SLOC_UDIFFUSION_VOL( IPHASE, CV_SILOC ) = UDIFFUSION_VOL_ALL( IPHASE, MAT_INOD )
                                SLOC2_UDIFFUSION( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, CV_SILOC ) = UDIFFUSION_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, MAT_INOD2 )
                                SLOC2_UDIFFUSION_VOL( IPHASE, CV_SILOC ) = UDIFFUSION_VOL_ALL( IPHASE, MAT_INOD2 )
                            ELSE
                                ! set to 0.0 for free-slip
                                SLOC_UDIFFUSION( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, CV_SILOC ) = UDIFFUSION( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, MAT_INOD )
                                if ( UDIFFUSION_VOL%have_field ) SLOC_UDIFFUSION_VOL( IPHASE, CV_SILOC ) = UDIFFUSION_VOL%val( IPHASE, 1,1,MAT_INOD )
                                SLOC2_UDIFFUSION( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, CV_SILOC ) = UDIFFUSION( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, MAT_INOD2 )
                                if ( UDIFFUSION_VOL%have_field ) SLOC2_UDIFFUSION_VOL( IPHASE, CV_SILOC ) = UDIFFUSION_VOL%val( IPHASE, 1,1,MAT_INOD2 )
                            ENDIF
                        END DO
                    END IF
                END DO
                DO U_SILOC = 1, Mdims%u_snloc
                    U_ILOC = U_SLOC2LOC( U_SILOC )
                    U_INOD = ndgln%u( (ELE-1)*Mdims%u_nloc + U_ILOC )
                    IF ( ELE2 /= 0 ) THEN
                        U_ILOC2 = U_ILOC_OTHER_SIDE( U_SILOC )
                        U_INOD2 = ndgln%u((ELE2-1)*Mdims%u_nloc + U_ILOC2 )
                    ELSE
                        U_ILOC2 = U_ILOC
                        U_INOD2 = U_INOD
                    END IF
                    DO IPHASE = 1, Mdims%nphase
                        IF ( GOT_DIFFUS ) THEN
                            SLOC_DUX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_SILOC ) = DUX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_ILOC, ELE )
                            SLOC_DUOLDX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_SILOC ) = DUOLDX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_ILOC, ELE )
                            IF(ELE2 /= 0) THEN
                                SLOC2_DUX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_SILOC ) = DUX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_ILOC2, ELE2 )
                                SLOC2_DUOLDX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_SILOC ) = DUOLDX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_ILOC2, ELE2 )
                            ELSE
                                SLOC2_DUX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_SILOC ) = DUX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_ILOC, ELE )
                                SLOC2_DUOLDX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_SILOC ) = DUOLDX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_ILOC, ELE )
                            END IF
                        END IF
                    END DO
                END DO
                ! velocity variables...
                DO U_SILOC = 1, Mdims%u_snloc
                    U_ILOC = U_SLOC2LOC( U_SILOC )
                    U_INOD = ndgln%u( (ELE-1)*Mdims%u_nloc + U_ILOC )
                    IF ( ELE2 /= 0) THEN
                        U_ILOC2 = U_ILOC_OTHER_SIDE( U_SILOC )
                        U_INOD2 = ndgln%u( (ELE2-1)*Mdims%u_nloc + U_ILOC2 )
                    ELSE
                        U_ILOC2 = U_ILOC
                        U_INOD2 = U_INOD
                    END IF
                    ! for normal calc...
                    DO IPHASE = 1, Mdims%nphase
                        IF ( GOT_DIFFUS ) THEN
                            SLOC_DUX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_SILOC ) = DUX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_ILOC, ELE )
                            SLOC_DUOLDX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_SILOC ) = DUOLDX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_ILOC, ELE )
                            IF(ELE2 /= 0) THEN
                                SLOC2_DUX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_SILOC ) = DUX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_ILOC2, ELE2 )
                                SLOC2_DUOLDX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_SILOC ) = DUOLDX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_ILOC2, ELE2 )
                            ELSE
                                SLOC2_DUX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_SILOC ) = DUX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_ILOC, ELE )
                                SLOC2_DUOLDX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_SILOC ) = DUOLDX_ELE_ALL( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, U_ILOC, ELE )
                            END IF
                        END IF
                        IF ( BETWEEN_ELE_STAB ) THEN
                            ! Calculate stabilization diffusion coefficient...
                            DO IDIM_VEL = 1, Mdims%ndim
                                SLOC_DIFF_FOR_BETWEEN_U( IDIM_VEL, IPHASE, U_SILOC ) = DIFF_FOR_BETWEEN_U( IDIM_VEL, IPHASE, U_ILOC, ELE )
                                IF ( ELE2 /= 0 ) THEN
                                    SLOC2_DIFF_FOR_BETWEEN_U( IDIM_VEL, IPHASE, U_SILOC ) = DIFF_FOR_BETWEEN_U( IDIM_VEL, IPHASE, U_ILOC2, ELE2 )
                                ELSE
                                    SLOC2_DIFF_FOR_BETWEEN_U( IDIM_VEL, IPHASE, U_SILOC ) = DIFF_FOR_BETWEEN_U( IDIM_VEL, IPHASE, U_ILOC, ELE )
                                END IF
                            END DO
                        END IF
                        ! U:
                        DO IDIM = 1, Mdims%ndim
                            SLOC_U( IDIM, IPHASE, U_SILOC ) = U_ALL( IDIM, IPHASE, U_INOD )
                            SLOC_UOLD( IDIM, IPHASE, U_SILOC ) = UOLD_ALL( IDIM, IPHASE, U_INOD )
                            SLOC2_U( IDIM, IPHASE, U_SILOC ) = U_ALL( IDIM, IPHASE, U_INOD2 )
                            SLOC2_UOLD( IDIM, IPHASE, U_SILOC ) = UOLD_ALL( IDIM, IPHASE, U_INOD2 )
                        END DO
                        DO IDIM = 1, Mdims%ndim
                            SLOC_NU( IDIM, IPHASE, U_SILOC ) = NU_ALL( IDIM, IPHASE, U_INOD )
                            SLOC_NUOLD( IDIM, IPHASE, U_SILOC ) = NUOLD_ALL( IDIM, IPHASE, U_INOD )
                            SLOC2_NU( IDIM, IPHASE, U_SILOC ) = NU_ALL( IDIM, IPHASE, U_INOD2 )
                            SLOC2_NUOLD( IDIM, IPHASE, U_SILOC ) = NUOLD_ALL( IDIM, IPHASE, U_INOD2 )
                        END DO
                    END DO
                END DO
                If_diffusion_or_momentum1: IF(GOT_DIFFUS .OR. GOT_UDEN) THEN
                    SDEN=0.0
                    SDENOLD=0.0
                    SDEN_KEEP=0.0 ; SDEN2_KEEP=0.0
                    SDENOLD_KEEP=0.0 ; SDENOLD2_KEEP=0.0
                    IF(IDIVID_BY_VOL_FRAC+IGOT_VOL_X_PRESSURE.GE.1) THEN
                        SVOL_FRA =0.0
                        SVOL_FRA2=0.0
                    ENDIF
                    IF(GOT_VIRTUAL_MASS) THEN
                        SVIRTUAL_MASS_GI=0.0
                        SVIRTUAL_MASS_OLD_GI=0.0
                        SVIRTUAL_MASS_GI_KEEP=0.0
                        SVIRTUAL_MASS_GI2_KEEP=0.0
                        SVIRTUAL_MASS_OLD_GI_KEEP=0.0
                        SVIRTUAL_MASS_OLD_GI2_KEEP=0.0
                    ENDIF
                    DO CV_SILOC=1,Mdims%cv_snloc
                        DO SGI=1,FE_GIdims%sbcvngi
                            DO IPHASE=1, Mdims%nphase
                                SDEN(IPHASE,SGI)=SDEN(IPHASE,SGI) + SBCVFEN_REVERSED(SGI,CV_SILOC) &
                                    *0.5*(SLOC_UDEN(IPHASE,CV_SILOC)+SLOC2_UDEN(IPHASE,CV_SILOC)) *WITH_NONLIN
                                SDENOLD(IPHASE,SGI)=SDENOLD(IPHASE,SGI) + SBCVFEN_REVERSED(SGI,CV_SILOC) &
                                    *0.5*(SLOC_UDENOLD(IPHASE,CV_SILOC)+SLOC2_UDENOLD(IPHASE,CV_SILOC)) *WITH_NONLIN
                                SDEN_KEEP(IPHASE,SGI)=SDEN_KEEP(IPHASE,SGI) + SBCVFEN_REVERSED(SGI,CV_SILOC) &
                                    *SLOC_UDEN(IPHASE,CV_SILOC)*WITH_NONLIN
                                SDEN2_KEEP(IPHASE,SGI)=SDEN2_KEEP(IPHASE,SGI) + SBCVFEN_REVERSED(SGI,CV_SILOC) &
                                    *SLOC2_UDEN(IPHASE,CV_SILOC)*WITH_NONLIN
                                SDENOLD_KEEP(IPHASE,SGI)=SDENOLD_KEEP(IPHASE,SGI) + SBCVFEN_REVERSED(SGI,CV_SILOC) &
                                    *SLOC_UDENOLD(IPHASE,CV_SILOC)*WITH_NONLIN
                                SDENOLD2_KEEP(IPHASE,SGI)=SDENOLD2_KEEP(IPHASE,SGI) + SBCVFEN_REVERSED(SGI,CV_SILOC) &
                                    *SLOC2_UDENOLD(IPHASE,CV_SILOC)*WITH_NONLIN
                                IF(IDIVID_BY_VOL_FRAC+IGOT_VOL_X_PRESSURE.GE.1) THEN
                                    SVOL_FRA(IPHASE,SGI) =SVOL_FRA(IPHASE,SGI) + SBCVFEN_REVERSED(SGI,CV_SILOC) *SLOC_VOL_FRA(IPHASE,CV_SILOC)
                                    SVOL_FRA2(IPHASE,SGI)=SVOL_FRA2(IPHASE,SGI)+ SBCVFEN_REVERSED(SGI,CV_SILOC) *SLOC2_VOL_FRA(IPHASE,CV_SILOC)
                                ENDIF
                                IF(GOT_VIRTUAL_MASS) THEN
                                    SVIRTUAL_MASS_GI(IPHASE,:,SGI)=SVIRTUAL_MASS_GI(IPHASE,:,SGI) + SBCVFEN_REVERSED(SGI,CV_SILOC) &
                                        *0.5*(SLOC_VIRTUAL_MASS(IPHASE,:,CV_SILOC)+SLOC2_VIRTUAL_MASS(IPHASE,:,CV_SILOC)) *WITH_NONLIN_CVM
                                    SVIRTUAL_MASS_OLD_GI(IPHASE,:,SGI)=SVIRTUAL_MASS_OLD_GI(IPHASE,:,SGI) + SBCVFEN_REVERSED(SGI,CV_SILOC) &
                                        *0.5*(SLOC_VIRTUAL_MASS_OLD(IPHASE,:,CV_SILOC)+SLOC2_VIRTUAL_MASS_OLD(IPHASE,:,CV_SILOC)) *WITH_NONLIN_CVM
                                    SVIRTUAL_MASS_GI_KEEP(IPHASE,:,SGI)=SVIRTUAL_MASS_GI_KEEP(IPHASE,:,SGI) + SBCVFEN_REVERSED(SGI,CV_SILOC) &
                                        *SLOC_VIRTUAL_MASS(IPHASE,:,CV_SILOC) *WITH_NONLIN_CVM
                                    SVIRTUAL_MASS_GI2_KEEP(IPHASE,:,SGI)=SVIRTUAL_MASS_GI2_KEEP(IPHASE,:,SGI) + SBCVFEN_REVERSED(SGI,CV_SILOC) &
                                        *SLOC2_VIRTUAL_MASS(IPHASE,:,CV_SILOC) *WITH_NONLIN_CVM
                                    SVIRTUAL_MASS_OLD_GI_KEEP(IPHASE,:,SGI)=SVIRTUAL_MASS_OLD_GI_KEEP(IPHASE,:,SGI) + SBCVFEN_REVERSED(SGI,CV_SILOC) &
                                        *SLOC_VIRTUAL_MASS_OLD(IPHASE,:,CV_SILOC) *WITH_NONLIN_CVM
                                    SVIRTUAL_MASS_OLD_GI2_KEEP(IPHASE,:,SGI)=SVIRTUAL_MASS_OLD_GI2_KEEP(IPHASE,:,SGI) + SBCVFEN_REVERSED(SGI,CV_SILOC) &
                                        *SLOC2_VIRTUAL_MASS_OLD(IPHASE,:,CV_SILOC) *WITH_NONLIN_CVM
                                ENDIF
                            END DO
                        END DO
                    END DO
                    SUD_ALL=0.0
                    SUDOLD_ALL=0.0
                    DO U_SILOC=1,Mdims%u_snloc
                        DO SGI=1,FE_GIdims%sbcvngi
                            DO IPHASE=1, Mdims%nphase
                                SUD_ALL(:,IPHASE,SGI)   =SUD_ALL(:,IPHASE,SGI)    + SBUFEN_REVERSED(SGI,U_SILOC)*SLOC_NU(:,IPHASE,U_SILOC)
                                SUDOLD_ALL(:,IPHASE,SGI)=SUDOLD_ALL(:,IPHASE,SGI) + SBUFEN_REVERSED(SGI,U_SILOC)*SLOC_NUOLD(:,IPHASE,U_SILOC)
                            END DO
                        END DO
                    END DO
                    SUD_ALL_KEEP=SUD_ALL
                    SUDOLD_ALL_KEEP=SUDOLD_ALL
                ENDIF If_diffusion_or_momentum1
                If_on_boundary_domain: IF(SELE /= 0) THEN
                    ! ***********SUBROUTINE DETERMINE_SUF_PRES - START************
                    ! Put the surface integrals in for pressure b.Mmat%C.'s
                    ! that is add into Mmat%C matrix and Mmat%U_RHS. (DG velocities)
                    Loop_ILOC2: DO U_SILOC = 1, Mdims%u_snloc
                        U_ILOC = U_SLOC2LOC( U_SILOC )
                        if(.not.got_c_matrix) IU_NOD = ndgln%suf_u(( SELE - 1 ) * Mdims%u_snloc + U_SILOC )
                        Loop_JLOC2: DO P_SJLOC = 1, Mdims%p_snloc
                            P_JLOC = CV_SLOC2LOC( P_SJLOC )
                            if(.not.got_c_matrix) JCV_NOD = ndgln%suf_p(( SELE - 1 ) * Mdims%p_snloc + P_SJLOC )
                            !Calculate aid variable NMX_ALL
                            NMX_ALL = matmul(SNORMXN_ALL( :, : ), SBUFEN_REVERSED( :, U_SILOC ) * SBCVFEN_REVERSED( :, P_SJLOC ) * SDETWE( : ))

                            IF(IGOT_VOL_X_PRESSURE==1) THEN
                                DO IPHASE = 1, Mdims%nphase
                                    VOL_FRA_NMX_ALL( :, IPHASE ) = VOL_FRA_NMX_ALL( :, IPHASE ) + sum(SVOL_FRA( IPHASE, : )) * NMX_ALL( : )
                                END DO
                            END IF
                            ! Put into matrix
                            ! Find COUNT - position in matrix : Mspars%MCY%fin, Mspars%MCY%col
                            IF ( .NOT.GOT_C_MATRIX ) THEN
                                CALL USE_POSINMAT_C_STORE( COUNT, IU_NOD, JCV_NOD, &
                                    Mdims%u_nonods, Mspars%C%fin, Mspars%C%col, Mspars%C%ncol, &
                                    IDO_STORE_AC_SPAR_PT, STORED_AC_SPAR_PT, POSINMAT_C_STORE, ELE, U_ILOC, P_JLOC, &
                                    Mdims%totele, Mdims%u_nloc, Mdims%p_nloc )
                            END IF
                            if (.not.Mmat%CV_pressure) then!If Mmat%C_CV then BCs are introduced in CV_ASSEMB
                                DO IPRES = 1, Mdims%npres
                                    IF( WIC_P_BC_ALL( 1,IPRES,SELE ) == WIC_P_BC_DIRICHLET ) THEN
                                        DO IPHASE =  1+(IPRES-1)*Mdims%n_in_pres, IPRES*Mdims%n_in_pres
                                            DO IDIM = 1, Mdims%ndim
                                                IF ( IGOT_VOL_X_PRESSURE == 1 ) THEN
                                                    IF ( .NOT.GOT_C_MATRIX ) THEN
                                                        Mmat%C( IDIM, IPHASE, COUNT ) = Mmat%C( IDIM, IPHASE, COUNT ) &
                                                            + VOL_FRA_NMX_ALL( IDIM, IPHASE )
                                                    END IF
                                                    LOC_U_RHS( IDIM, IPHASE, U_ILOC) =  LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                        - VOL_FRA_NMX_ALL( IDIM, IPHASE ) * SUF_P_BC_ALL( 1,IPRES,P_SJLOC + Mdims%p_snloc * ( SELE - 1) )
                                                ELSE
                                                    IF ( .NOT.GOT_C_MATRIX ) THEN
                                                        Mmat%C( IDIM, IPHASE, COUNT ) = Mmat%C( IDIM, IPHASE, COUNT ) &
                                                            + NMX_ALL( IDIM )
                                                    END IF
                                                    LOC_U_RHS( IDIM, IPHASE, U_ILOC) =  LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                        - NMX_ALL( IDIM ) * SUF_P_BC_ALL( 1,IPRES,P_SJLOC + Mdims%p_snloc* ( SELE - 1 ) )
                                                END IF
                                            END DO
                                        END DO
                                    END IF
                                END DO
                            end if
                        END DO Loop_JLOC2
                    END DO Loop_ILOC2
                    if ( got_free_surf ) then
                        if ( symmetric_P ) then ! us this if not the other formation of MASS_SUF in cv_adv_dif
                            if ( WIC_P_BC_ALL( 1,1,SELE ) == WIC_P_BC_FREE ) then ! on the free surface...
                                DO P_SILOC = 1, Mdims%p_snloc
                                    P_ILOC = CV_SLOC2LOC( P_SILOC )
                                    P_INOD=ndgln%cv( (ELE-1)*Mdims%cv_nloc + P_ILOC )
                                    DO P_SJLOC = 1, Mdims%p_snloc
                                        P_JLOC = CV_SLOC2LOC( P_SJLOC )
                                        P_JNOD = ndgln%cv( (ELE-1)*Mdims%cv_nloc + P_JLOC )
                                        ! Use the same sparcity as the MN matrix...
                                        COUNT_SUF=0
                                        DO COUNT_P = Mspars%CMC%fin(P_INOD), Mspars%CMC%fin(P_INOD+1)-1
                                            IF ( Mspars%CMC%col(COUNT_P) == P_JNOD ) THEN
                                                COUNT_SUF = COUNT_P
                                                EXIT
                                            END IF
                                        END DO
                                        MM_GRAVTY = sum( SNORMXN_ALL( 3, : ) * SBCVFEN_REVERSED( :, P_SILOC ) * &
                                            SBCVFEN_REVERSED( :, P_SJLOC ) * SDETWE( : ) ) / ( DT**2 * GRAVTY )
                                        MASS_SUF( COUNT_SUF ) = MASS_SUF( COUNT_SUF ) + MM_GRAVTY
                                    END DO
                                END DO
                            end if
                        end if
                    end if
                   ! ***********SUBROUTINE DETERMINE_SUF_PRES - END************
                ENDIF If_on_boundary_domain
                If_ele2_notzero: IF(ELE2 /= 0) THEN
                    got_c_matrix1: if(.not.got_c_matrix) then
                        discontinuous_pres: IF(DISC_PRES) THEN
                            IF( VOL_ELE_INT_PRES ) THEN
                                ! bias the weighting towards bigger eles - works with 0.25 and 0.1 and not 0.01.
                                MASSE = MASS_ELE( ELE ) + 0.25 * MASS_ELE( ELE2 )
                                MASSE2 = MASS_ELE( ELE2 ) + 0.25 * MASS_ELE( ELE )
                            ELSE ! Simple average (works well with Mdisopt%in_ele_upwind=Mdisopt%dg_ele_upwind=2)...
                                MASSE = 1.0
                                MASSE2 = 1.0
                            END IF
                            DO P_SJLOC = 1, Mdims%cv_snloc
                                P_JLOC = CV_SLOC2LOC( P_SJLOC )
                                P_JNOD = ndgln%p(( ELE - 1 ) * Mdims%p_nloc + P_JLOC )
                                P_JLOC2 = MAT_OTHER_LOC(P_JLOC)
                                P_JNOD2 = ndgln%p(( ELE2 - 1 ) * Mdims%p_nloc + P_JLOC2 )
                                DO U_SILOC = 1, Mdims%u_snloc
                                    U_ILOC = U_SLOC2LOC( U_SILOC )
                                    U_INOD = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                                    VNMX_ALL = 0.0
                                    DO SGI = 1, FE_GIdims%sbcvngi
                                        RNN = SDETWE( SGI ) * SBUFEN_REVERSED( SGI, U_SILOC ) * SBCVFEN_REVERSED( SGI, P_SJLOC )
                                        VNMX_ALL = VNMX_ALL + SNORMXN_ALL( :, SGI ) * RNN
                                    END DO
                                    CALL USE_POSINMAT_C_STORE( COUNT, U_INOD, P_JNOD,  &
                                        Mdims%u_nonods, Mspars%C%fin, Mspars%C%col, Mspars%C%ncol, &
                                        IDO_STORE_AC_SPAR_PT, STORED_AC_SPAR_PT, POSINMAT_C_STORE, ELE, U_ILOC, P_JLOC, &
                                        Mdims%totele, Mdims%u_nloc, Mdims%p_nloc )
                                    !NOT VERY ELEGANT
                                    IF(DISC_PRES) CALL USE_POSINMAT_C_STORE_SUF_DG( COUNT2, U_INOD, P_JNOD2,  &
                                        Mdims%u_nonods, Mspars%C%fin, Mspars%C%col, Mspars%C%ncol, &
                                        IDO_STORE_AC_SPAR_PT, STORED_AC_SPAR_PT, POSINMAT_C_STORE_SUF_DG, ELE, IFACE, U_SILOC, P_SJLOC,  &
                                        Mdims%totele, FE_GIdims%nface, Mdims%u_snloc, Mdims%p_snloc )
                                    Loop_Phase5: DO IPHASE = 1, Mdims%nphase
                                        ! weight integral according to non-uniform mesh spacing otherwise it will go unstable.
                                        IF ( .NOT.GOT_C_MATRIX ) THEN
                                            DO IDIM = 1, Mdims%ndim
                                                Mmat%C( IDIM, IPHASE, COUNT ) = Mmat%C( IDIM, IPHASE, COUNT ) &
                                                    + VNMX_ALL( IDIM ) * MASSE / ( MASSE + MASSE2 )
                                                Mmat%C( IDIM, IPHASE, COUNT2 ) = Mmat%C( IDIM, IPHASE, COUNT2 ) &
                                                    - VNMX_ALL( IDIM ) *  MASSE / ( MASSE + MASSE2 )
                                            END DO
                                        END IF
                                    END DO Loop_Phase5
                                END DO
                            END DO
                           !STOP 383
                        ENDIF discontinuous_pres
                    ENDIF got_c_matrix1
                    If_diffusion_or_momentum2: IF(GOT_DIFFUS .OR. GOT_UDEN) THEN
                        ! Calculate distance between centres of elements HDC
                        DO CV_ILOC = 1, Mdims%cv_nloc
                            X_INOD = ndgln%x( (ELE2-1)*Mdims%x_nloc + CV_ILOC )
                            XL2_ALL(:,CV_ILOC) = X_ALL( :, X_INOD )
                        END DO
                        DO IDIM = 1, Mdims%ndim
                            C1 ( IDIM ) = SUM( XL_ALL( IDIM, : ) ) / REAL( Mdims%x_nloc )
                            C2 ( IDIM ) = SUM( XL2_ALL( IDIM, : ) ) / REAL( Mdims%x_nloc )
                        END DO
                        HDC = SQRT( SUM( ( C1 - C2 )**2 ) )
                        SUD2_ALL=0.0
                        SUDOLD2_ALL=0.0
                        DO U_SILOC=1,Mdims%u_snloc
                            DO SGI=1,FE_GIdims%sbcvngi
                                DO IPHASE=1, Mdims%nphase
                                    SUD2_ALL(:,IPHASE,SGI)   =SUD2_ALL(:,IPHASE,SGI)    + SBUFEN_REVERSED(SGI,U_SILOC)*SLOC_NU(:,IPHASE,U_SILOC)
                                    SUDOLD2_ALL(:,IPHASE,SGI)=SUDOLD2_ALL(:,IPHASE,SGI) + SBUFEN_REVERSED(SGI,U_SILOC)*SLOC_NUOLD(:,IPHASE,U_SILOC)
                                END DO
                            END DO
                        END DO
                        SUD2_ALL_KEEP=SUD2_ALL
                        SUDOLD2_ALL_KEEP=SUDOLD2_ALL
                        IF(MOM_CONSERV) THEN
                            SUD_ALL=0.5*(SUD_ALL+SUD2_ALL)
                            SUDOLD_ALL=0.5*(SUDOLD_ALL+SUDOLD2_ALL)
                        ENDIF
                    ENDIF If_diffusion_or_momentum2
                ELSE
                    DO IDIM = 1, Mdims%ndim
                        C1 ( IDIM ) = SUM( XL_ALL( IDIM, : ) ) / REAL( Mdims%x_nloc )
                        C2 ( IDIM ) = SUM( XSL_ALL( IDIM, : ) ) / REAL( Mdims%cv_snloc )
                        AVE_SNORMXN_ALL(IDIM) = SUM( SNORMXN_ALL( IDIM, : )) / REAL(FE_GIdims%sbcvngi)
                    END DO
                    ! use 2* because the value of HD being used is measured between the centres of neighbouring elements.
                    !                    HDC = 2.*SQRT( SUM( ( C1 - C2 )**2 ) )
                    ! Take the mean normal component of the distance vector to the surface...
                    HDC = 2.*SQRT( SUM( (( C1 - C2 )*AVE_SNORMXN_ALL)**2 ) )
                END IF If_ele2_notzero
                IF ( GOT_UDEN ) THEN
                    IF ( MOM_CONSERV ) THEN
                        IF ( SELE2 /= 0 ) THEN
                            SUD2_ALL=0.0 ; SUDOLD2_ALL=0.0
                            DO IPHASE = 1, Mdims%nphase
                                DO IDIM = 1, Mdims%ndim
                                    IF ( WIC_U_BC_ALL_ADV( IDIM, IPHASE, SELE2 ) == WIC_U_BC_DIRICHLET ) THEN
                                        DO U_SILOC = 1, Mdims%u_snloc
                                            DO SGI = 1, FE_GIdims%sbcvngi
                                                SUD2_ALL(IDIM,IPHASE,SGI) = SUD2_ALL(IDIM,IPHASE,SGI) + &
                                                    SBUFEN_REVERSED(SGI,U_SILOC) * suf_nu_bc_all( idim,iphase,u_siloc + Mdims%u_snloc * ( sele2 - 1 ) )
                                                SUDOLD2_ALL(IDIM,IPHASE,SGI) = SUDOLD2_ALL(IDIM,IPHASE,SGI) + &
                                                    SBUFEN_REVERSED(SGI,U_SILOC) * suf_nu_bc_all( idim,iphase,u_siloc + Mdims%u_snloc * ( sele2 - 1 ) )
                                            END DO
                                        END DO
                                    END IF
                                END DO
                                DO SGI = 1, FE_GIdims%sbcvngi
                                    SNDOTQ(IPHASE,SGI) = SUM( SUD_ALL(:,IPHASE,SGI) * SNORMXN_ALL(:,SGI) )
                                    SNDOTQOLD(IPHASE,SGI) = SUM( SUDOLD_ALL(:,IPHASE,SGI) * SNORMXN_ALL(:,SGI) )
                                END DO
                            END DO
                            SINCOME = 0.5 + 0.5 * SIGN( 1.0, -SNDOTQ )
                            SINCOMEOLD = 0.5 + 0.5 * SIGN( 1.0, -SNDOTQOLD )
                            DO IPHASE = 1, Mdims%nphase
                                DO IDIM = 1, Mdims%ndim
                                    IF ( WIC_U_BC_ALL( IDIM, IPHASE, SELE2 ) == WIC_U_BC_DIRICHLET ) THEN
                                        DO SGI = 1, FE_GIdims%sbcvngi
                                            SUD_ALL(IDIM,IPHASE,SGI) = SINCOME(IPHASE,SGI) * 0.5 * (SUD_ALL(IDIM,IPHASE,SGI) + SUD2_ALL(IDIM,IPHASE,SGI) )  &
                                                + (1.-SINCOME(IPHASE,SGI)) * SUD_ALL(IDIM,IPHASE,SGI)
                                            SUDOLD_ALL(IDIM,IPHASE,SGI) = SINCOMEOLD(IPHASE,SGI) * 0.5 * (SUDOLD_ALL(IDIM,IPHASE,SGI) + SUDOLD2_ALL(IDIM,IPHASE,SGI) )  &
                                                + (1.-SINCOMEOLD(IPHASE,SGI)) * SUDOLD_ALL(IDIM,IPHASE,SGI)
                                        END DO
                                    END IF
                                END DO
                            END DO
                        END IF
                    END IF
                END IF
                If_diffusion_or_momentum3: IF(GOT_DIFFUS .OR. GOT_UDEN) THEN
                    IF(BETWEEN_ELE_STAB) THEN
                        ! Calculate stabilization diffusion coefficient...
                        UDIFF_SUF_STAB=0.0
                        DO U_SILOC = 1, Mdims%u_snloc
                            DO IPHASE=1,Mdims%nphase
                                DO IDIM_VEL=1,Mdims%ndim
                                    DO IDIM=1,Mdims%ndim
                                        UDIFF_SUF_STAB(IDIM_VEL,IDIM,IDIM,IPHASE,: ) = UDIFF_SUF_STAB(IDIM_VEL,IDIM,IDIM,IPHASE,: )  &
                                            +FE_funs%sbufen(U_SILOC,:)*0.5*(  SLOC_DIFF_FOR_BETWEEN_U(IDIM_VEL, IPHASE, U_SILOC) &
                                            + SLOC2_DIFF_FOR_BETWEEN_U(IDIM_VEL, IPHASE, U_SILOC)  )
                                    END DO
                                END DO
                            END DO
                        END DO
                    END IF
                    DO IPHASE = 1, Mdims%nphase
                        DO SGI = 1, FE_GIdims%sbcvngi
                            SNDOTQ(IPHASE,SGI)    = SUM( SUD_ALL(:,IPHASE,SGI)*SNORMXN_ALL(:,SGI) )
                            SNDOTQOLD(IPHASE,SGI) = SUM( SUDOLD_ALL(:,IPHASE,SGI)*SNORMXN_ALL(:,SGI) )
                        END DO
                    END DO
                    SINCOME = 0.5 + 0.5 * SIGN( 1.0, -SNDOTQ )
                    SINCOMEOLD = 0.5 + 0.5 * SIGN( 1.0, -SNDOTQOLD )
                    SNDOTQ_IN  = 0.0
                    SNDOTQ_OUT = 0.0
                    SNDOTQOLD_IN  = 0.0
                    SNDOTQOLD_OUT = 0.0
                    IF( NON_LIN_DGFLUX ) THEN
                        DO IPHASE=1, Mdims%nphase
                            DO SGI=1,FE_GIdims%sbcvngi
                                SNDOTQ_KEEP(IPHASE,SGI)   = SUM( SUD_ALL_KEEP(:,IPHASE,SGI)*SNORMXN_ALL(:,SGI) )
                                SNDOTQ2_KEEP(IPHASE,SGI)   =SUM( SUD2_ALL_KEEP(:,IPHASE,SGI)*SNORMXN_ALL(:,SGI)  )
                                SNDOTQOLD_KEEP(IPHASE,SGI)   = SUM( SUDOLD_ALL_KEEP(:,IPHASE,SGI)*SNORMXN_ALL(:,SGI) )
                                SNDOTQOLD2_KEEP(IPHASE,SGI)   =SUM( SUDOLD2_ALL_KEEP(:,IPHASE,SGI)*SNORMXN_ALL(:,SGI)  )
                            END DO
                        END DO
                        ELE3 = ELE2
                        IF ( ELE2==0 ) ELE3 = ELE
                        N_DOT_DU=0.0
                        N_DOT_DU2=0.0
                        N_DOT_DUOLD=0.0
                        N_DOT_DUOLD2=0.0
                        DO U_SILOC = 1, Mdims%u_snloc
                            U_ILOC = U_SLOC2LOC( U_SILOC )
                            DO IPHASE = 1, Mdims%nphase
                                do sgi = 1, FE_GIdims%sbcvngi
                                    vel_dot(sgi)  =  sum( SUD_ALL(:,IPHASE,SGI) *snormxn_all(:,sgi) )
                                    vel_dot2(sgi) =  sum( SUD2_ALL(:,IPHASE,SGI)*snormxn_all(:,sgi) )
                                    velold_dot(sgi)  = sum( SUDOLD_ALL(:,IPHASE,SGI) *snormxn_all(:,sgi) )
                                    velold_dot2(sgi) = sum( SUDOLD2_ALL(:,IPHASE,SGI) *snormxn_all(:,sgi) )
                                    grad_fact(sgi) = sum( DevFuns%UFENX_ALL(1:Mdims%ndim,U_ILOC,1)*snormxn_ALL(:,SGI) )
                                end do
                                N_DOT_DU(iphase,:)  = N_DOT_DU(iphase,:)  + grad_fact(:)*vel_dot(:)
                                N_DOT_DU2(iphase,:) = N_DOT_DU2(iphase,:) + grad_fact(:)*vel_dot2(:)
                                N_DOT_DUOLD(iphase,:) = N_DOT_DUOLD(iphase,:)  + grad_fact(:)*velold_dot(:)
                                N_DOT_DUOLD2(iphase,:) = N_DOT_DUOLD2(iphase,:) + grad_fact(:)*velold_dot2(:)
                            END DO
                        END DO
                    END IF  ! endof IF( NON_LIN_DGFLUX ) THEN
                    ! Have a surface integral on element boundary...
                    ! Calculate the velocities either side of the element...
                    U_NODJ_SGI_IPHASE_ALL=0.0 ; U_NODI_SGI_IPHASE_ALL=0.0
                    UOLD_NODJ_SGI_IPHASE_ALL=0.0 ; UOLD_NODI_SGI_IPHASE_ALL=0.0
                    DO U_SILOC = 1, Mdims%u_snloc
                        DO SGI=1,FE_GIdims%sbcvngi
                            DO IPHASE=1, Mdims%nphase
                                U_NODI_SGI_IPHASE_ALL(:,IPHASE,SGI) = U_NODI_SGI_IPHASE_ALL(:,IPHASE,SGI) + SBUFEN_REVERSED(SGI,U_SILOC) * SLOC_U(:,IPHASE,U_SILOC)
                                U_NODJ_SGI_IPHASE_ALL(:,IPHASE,SGI) = U_NODJ_SGI_IPHASE_ALL(:,IPHASE,SGI) + SBUFEN_REVERSED(SGI,U_SILOC) * SLOC2_U(:,IPHASE,U_SILOC)
                                UOLD_NODI_SGI_IPHASE_ALL(:,IPHASE,SGI) = UOLD_NODI_SGI_IPHASE_ALL(:,IPHASE,SGI) + SBUFEN_REVERSED(SGI,U_SILOC) * SLOC_UOLD(:,IPHASE,U_SILOC)
                                UOLD_NODJ_SGI_IPHASE_ALL(:,IPHASE,SGI) = UOLD_NODJ_SGI_IPHASE_ALL(:,IPHASE,SGI) + SBUFEN_REVERSED(SGI,U_SILOC) * SLOC2_UOLD(:,IPHASE,U_SILOC)
                            END DO
                        END DO
                    END DO
                    IF ( SELE2 /= 0 ) THEN
                        DO IPHASE = 1, Mdims%nphase
                            DO IDIM = 1, Mdims%ndim
                                IF ( WIC_U_BC_ALL_VISC(IDIM,IPHASE,SELE2 ) == WIC_U_BC_DIRICHLET ) THEN
                                    U_NODJ_SGI_IPHASE_ALL(IDIM,IPHASE,:) = 0.0
                                    UOLD_NODJ_SGI_IPHASE_ALL(IDIM,IPHASE,:) = 0.0
                                    DO U_SILOC = 1, Mdims%u_snloc
                                        DO SGI = 1, FE_GIdims%sbcvngi
                                            U_NODJ_SGI_IPHASE_ALL(IDIM,IPHASE,SGI) = U_NODJ_SGI_IPHASE_ALL(IDIM,IPHASE,SGI) + &
                                                SBUFEN_REVERSED(SGI,U_SILOC) * SUF_U_BC_ALL_VISC( IDIM,IPHASE,U_SILOC + Mdims%u_snloc*(SELE2-1) )
                                            UOLD_NODJ_SGI_IPHASE_ALL(IDIM,IPHASE,SGI) = UOLD_NODJ_SGI_IPHASE_ALL(IDIM,IPHASE,SGI) + &
                                                SBUFEN_REVERSED(SGI,U_SILOC) * SUF_U_BC_ALL_VISC( IDIM,IPHASE,U_SILOC + Mdims%u_snloc*(SELE2-1) )
                                        END DO
                                    END DO
                                ENDIF
                            END DO
                        END DO
                    END IF
                    !                    IF(GOT_DIFFUS.and.LINEAR_HIGHORDER_DIFFUSION) STRESS_IJ_ELE_EXT=0.0
                    ! This sub should be used for stress and tensor viscocity replacing the rest...
                    If_GOT_DIFFUS2: IF(GOT_DIFFUS.and.LINEAR_HIGHORDER_DIFFUSION) THEN
                        !                    If_GOT_DIFFUS2: IF(GOT_DIFFUS.and.LINEAR_HIGHORDER_DIFFUSION.and.(ele2.ne.0)) THEN
                        ! only used between elements of the domain so no modification of b.Mmat%C's nec.
                        !                        STRESS_IJ_ELE_EXT=0.0
                        CALL LINEAR_HIGH_DIFFUS_CAL_COEFF_STRESS_OR_TENSOR( STRESS_IJ_ELE_EXT, S_INV_NNX_MAT12,  &
                            STRESS_FORM, STRESS_FORM_STAB, ZERO_OR_TWO_THIRDS, &
                            Mdims%u_snloc, Mdims%u_nloc, Mdims%cv_snloc, Mdims%nphase, &
                            SBUFEN_REVERSED,SBCVFEN_REVERSED,SDETWE, FE_GIdims%sbcvngi, Mdims%ndim, SLOC_UDIFFUSION, SLOC_UDIFFUSION_VOL, SLOC2_UDIFFUSION, SLOC2_UDIFFUSION_VOL, UDIFF_SUF_STAB, &
                            (ELE2.LE.0), SNORMXN_ALL  )
                        !                        STRESS_IJ_ELE_EXT=0.0
                        DIFF_COEF_DIVDX   =0.0
                        DIFF_COEFOLD_DIVDX=0.0
                    ELSE IF(GOT_DIFFUS) THEN If_GOT_DIFFUS2
                        CALL DIFFUS_CAL_COEFF_STRESS_OR_TENSOR( Mdims, DIFF_COEF_DIVDX, &
                            DIFF_COEFOLD_DIVDX, STRESS_FORM, STRESS_FORM_STAB, ZERO_OR_TWO_THIRDS, &
                            SBUFEN_REVERSED,SBCVFEN_REVERSED,FE_GIdims%sbcvngi, SLOC_UDIFFUSION, SLOC_UDIFFUSION_VOL, SLOC2_UDIFFUSION, SLOC2_UDIFFUSION_VOL, UDIFF_SUF_STAB, &
                            HDC, U_NODJ_SGI_IPHASE_ALL, U_NODI_SGI_IPHASE_ALL, &
                            UOLD_NODJ_SGI_IPHASE_ALL, UOLD_NODI_SGI_IPHASE_ALL, &
                            ELE, ELE2, SNORMXN_ALL,  &
                            SLOC_DUX_ELE_ALL, SLOC2_DUX_ELE_ALL,   SLOC_DUOLDX_ELE_ALL, SLOC2_DUOLDX_ELE_ALL,  &
                            SELE, WIC_U_BC_ALL_VISC, WIC_U_BC_DIRICHLET, SIMPLE_DIFF_CALC, DIFF_MIN_FRAC, DIFF_MAX_FRAC  )
                    ELSE If_GOT_DIFFUS2
                        DIFF_COEF_DIVDX   =0.0
                        DIFF_COEFOLD_DIVDX=0.0
                    END IF If_GOT_DIFFUS2
                    ! *************REVIEWER 3-END*************
                    ! *************REVIEWER 4-START*************
                    SNDOTQ_IN  = 0.0
                    SNDOTQ_OUT = 0.0
                    SNDOTQOLD_IN  = 0.0
                    SNDOTQOLD_OUT = 0.0
                    IF(GOT_VIRTUAL_MASS) THEN
                        CVM_SNDOTQ_IN = 0.0
                        CVM_SNDOTQ_OUT = 0.0
                        CVM_SNDOTQOLD_IN = 0.0
                        CVM_SNDOTQOLD_OUT = 0.0
                    ENDIF
                    DO SGI=1,FE_GIdims%sbcvngi
                        DO IPHASE=1, Mdims%nphase
                            DO IDIM=1, Mdims%ndim
                                !FTHETA( SGI,IDIM,IPHASE )=0.5 !1.0  - should be 1. as there is no theta set for the internal part of an element.
                                FTHETA( IDIM,IPHASE,SGI )=1.0 ! 0.5
                                ! CENT_RELAX=1.0 (central scheme) =0.0 (upwind scheme).
                                IF( NON_LIN_DGFLUX ) THEN
                                    ! non-linear DG flux - if we have an oscillation use upwinding else use central scheme.
                                    CENT_RELAX( IDIM,IPHASE,SGI ) = dg_oscilat_detect( SNDOTQ_KEEP(IPHASE,SGI), SNDOTQ2_KEEP(IPHASE,SGI), &
                                        N_DOT_DU(IPHASE,SGI), N_DOT_DU2(IPHASE,SGI), SINCOME(IPHASE,SGI), MASS_ELE(ELE), MASS_ELE(ELE2) )
                                    CENT_RELAX_OLD( IDIM,IPHASE,SGI )= dg_oscilat_detect( SNDOTQOLD_KEEP(IPHASE,SGI), SNDOTQOLD2_KEEP(IPHASE,SGI), &
                                        N_DOT_DUOLD(IPHASE,SGI), N_DOT_DUOLD2(IPHASE,SGI), SINCOMEOLD(IPHASE,SGI), MASS_ELE(ELE), MASS_ELE(ELE2) )
                                ELSE
                                    IF( UPWIND_DGFLUX ) THEN
                                        ! Upwind DG flux...
                                        CENT_RELAX( IDIM,IPHASE,SGI )    =0.0
                                        CENT_RELAX_OLD( IDIM,IPHASE,SGI )=0.0
                                    ELSE
                                        ! Central diff DG flux...
                                        CENT_RELAX( IDIM,IPHASE,SGI )    =1.0
                                        CENT_RELAX_OLD( IDIM,IPHASE,SGI )=1.0
                                    ENDIF
                                ENDIF
                               ! CENT_RELAX=1.0 (central scheme) =0.0 (upwind scheme).
                            END DO
                        END DO
                    END DO
                    DO SGI=1,FE_GIdims%sbcvngi
                        DO IPHASE=1, Mdims%nphase
                            DO IDIM=1, Mdims%ndim
                                ! CENT_RELAX=1.0 (central scheme) =0.0 (upwind scheme).
                                SNDOTQ_IN(IDIM,IPHASE,SGI)    =SNDOTQ_IN(IDIM,IPHASE,SGI)  &
                                    +FTHETA(IDIM,IPHASE,SGI)*SDEN(IPHASE,SGI)*SNDOTQ(IPHASE,SGI)  &
                                    * (0.5 * CENT_RELAX( IDIM,IPHASE,SGI ) + SINCOME(IPHASE,SGI)*(1.-CENT_RELAX( IDIM,IPHASE,SGI )))
                                SNDOTQ_OUT(IDIM,IPHASE,SGI)   =SNDOTQ_OUT(IDIM,IPHASE,SGI)  &
                                    +FTHETA(IDIM,IPHASE,SGI)*SDEN(IPHASE,SGI)*SNDOTQ(IPHASE,SGI) &
                                    * (0.5* CENT_RELAX( IDIM,IPHASE,SGI ) + (1.-SINCOME(IPHASE,SGI))*(1.-CENT_RELAX( IDIM,IPHASE,SGI )))
                                SNDOTQOLD_IN(IDIM,IPHASE,SGI) =SNDOTQOLD_IN(IDIM,IPHASE,SGI)  &
                                    +(1.-FTHETA(IDIM,IPHASE,SGI))*SDEN(IPHASE,SGI)*SNDOTQOLD(IPHASE,SGI)  &
                                    * (0.5* CENT_RELAX_OLD( IDIM,IPHASE,SGI ) + SINCOMEOLD(IPHASE,SGI)*(1.-CENT_RELAX_OLD( IDIM,IPHASE,SGI )))
                                SNDOTQOLD_OUT(IDIM,IPHASE,SGI)=SNDOTQOLD_OUT(IDIM,IPHASE,SGI)  &
                                    +(1.-FTHETA(IDIM,IPHASE,SGI))*SDEN(IPHASE,SGI)*SNDOTQOLD(IPHASE,SGI) &
                                    * (0.5* CENT_RELAX_OLD( IDIM,IPHASE,SGI ) + (1.-SINCOMEOLD(IPHASE,SGI))*(1.-CENT_RELAX_OLD( IDIM,IPHASE,SGI )))
                                IF(GOT_VIRTUAL_MASS) THEN
                                    DO KPHASE=1,Mdims%nphase
                                        CVM_SNDOTQ_IN(IDIM,IPHASE,KPHASE,SGI)    =CVM_SNDOTQ_IN(IDIM,IPHASE,KPHASE,SGI)  &
                                            +SVIRTUAL_MASS_GI(IPHASE,KPHASE,SGI)  &
                                            *(      FTHETA(IDIM,KPHASE,SGI)*VIRTUAL_MASS_ADV_CUR(IPHASE,KPHASE)*SNDOTQ(KPHASE,SGI)  &
                                            * ( 0.5 * CENT_RELAX( IDIM,KPHASE,SGI ) + SINCOME(KPHASE,SGI)*(1.-CENT_RELAX( IDIM,KPHASE,SGI )) ) &
                                            +FTHETA(IDIM,IPHASE,SGI)*(1.-VIRTUAL_MASS_ADV_CUR(IPHASE,KPHASE))*SNDOTQ(IPHASE,SGI)  &
                                            * ( 0.5 * CENT_RELAX( IDIM,IPHASE,SGI ) + SINCOME(IPHASE,SGI)*(1.-CENT_RELAX( IDIM,IPHASE,SGI )) )    )
                                        CVM_SNDOTQ_OUT(IDIM,IPHASE,KPHASE,SGI)   =CVM_SNDOTQ_OUT(IDIM,IPHASE,KPHASE,SGI)  &
                                            +SVIRTUAL_MASS_GI(IPHASE,KPHASE,SGI)  &
                                            *(     FTHETA(IDIM,KPHASE,SGI)*VIRTUAL_MASS_ADV_CUR(IPHASE,KPHASE)*SNDOTQ(KPHASE,SGI) &
                                            * ( 0.5* CENT_RELAX( IDIM,KPHASE,SGI ) + (1.-SINCOME(KPHASE,SGI))*(1.-CENT_RELAX( IDIM,KPHASE,SGI )) ) &
                                            +FTHETA(IDIM,IPHASE,SGI)*(1.-VIRTUAL_MASS_ADV_CUR(IPHASE,KPHASE))*SNDOTQ(IPHASE,SGI) &
                                            * ( 0.5* CENT_RELAX( IDIM,IPHASE,SGI ) + (1.-SINCOME(IPHASE,SGI))*(1.-CENT_RELAX( IDIM,IPHASE,SGI )) )    )
                                        CVM_SNDOTQOLD_IN(IDIM,IPHASE,KPHASE,SGI)    =CVM_SNDOTQOLD_IN(IDIM,IPHASE,KPHASE,SGI)  &
                                            +SVIRTUAL_MASS_OLD_GI(IPHASE,KPHASE,SGI)  &
                                            *(      (1.-FTHETA(IDIM,KPHASE,SGI))*VIRTUAL_MASS_ADV_CUR(IPHASE,KPHASE)*SNDOTQOLD(KPHASE,SGI)  &
                                            * (0.5 * CENT_RELAX_OLD( IDIM,KPHASE,SGI ) + SINCOMEOLD(KPHASE,SGI)*(1.-CENT_RELAX_OLD( IDIM,KPHASE,SGI ))) &
                                            +(1.-FTHETA(IDIM,IPHASE,SGI))*(1.-VIRTUAL_MASS_ADV_CUR(IPHASE,KPHASE))*SNDOTQOLD(IPHASE,SGI)  &
                                            * (0.5 * CENT_RELAX_OLD( IDIM,IPHASE,SGI ) + SINCOMEOLD(IPHASE,SGI)*(1.-CENT_RELAX_OLD( IDIM,IPHASE,SGI )))  )
                                        CVM_SNDOTQOLD_OUT(IDIM,IPHASE,KPHASE,SGI)   =CVM_SNDOTQOLD_OUT(IDIM,IPHASE,KPHASE,SGI)  &
                                            +SVIRTUAL_MASS_OLD_GI(IPHASE,KPHASE,SGI)  &
                                            *(     (1.-FTHETA(IDIM,IPHASE,SGI))*VIRTUAL_MASS_ADV_CUR(IPHASE,KPHASE)*SNDOTQOLD(KPHASE,SGI) &
                                            * (0.5* CENT_RELAX_OLD( IDIM,KPHASE,SGI ) + (1.-SINCOMEOLD(KPHASE,SGI))*(1.-CENT_RELAX_OLD( IDIM,KPHASE,SGI ))) &
                                            +(1.-FTHETA(IDIM,IPHASE,SGI))*(1.-VIRTUAL_MASS_ADV_CUR(IPHASE,KPHASE))*SNDOTQOLD(IPHASE,SGI) &
                                            * (0.5* CENT_RELAX_OLD( IDIM,IPHASE,SGI ) + (1.-SINCOMEOLD(IPHASE,SGI))*(1.-CENT_RELAX_OLD( IDIM,IPHASE,SGI )))    )
                                    END DO
                                ENDIF
                            END DO
                        END DO
                    END DO
                    IF(GOT_DIFFUS .AND. LINEAR_HIGHORDER_DIFFUSION) THEN
                        IF(ELE2.GT.0) THEN ! Internal to domain
                            DO U_JLOC=1,Mdims%u_nloc
                                U_JLOC2 = U_JLOC
                                DO U_SILOC=1,Mdims%u_snloc
                                    U_ILOC =U_SLOC2LOC(U_SILOC)
                                    DO IPHASE = 1, Mdims%nphase
                                        JPHASE = IPHASE
                                        DIAG_BIGM_CON(:,:,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                            =DIAG_BIGM_CON(:,:,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)       + STRESS_IJ_ELE_EXT( :,:, IPHASE, U_SILOC, U_JLOC )
                                        IF(PIVIT_ON_VISC) THEN
                                            DO IDIM=1,Mdims%ndim
                                                DO JDIM=1,Mdims%ndim
                                                    I = IDIM+(IPHASE-1)*Mdims%ndim+(U_ILOC-1)*Mdims%ndim*Mdims%nphase
                                                    J = JDIM+(JPHASE-1)*Mdims%ndim+(U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                                    Mmat%PIVIT_MAT(I,J,ELE) &
                                                        =Mmat%PIVIT_MAT(I,J,ELE) + ws * STRESS_IJ_ELE_EXT( IDIM, JDIM, IPHASE, U_SILOC, U_JLOC )
                                                END DO
                                            END DO
                                        ENDIF
                                        ! Contributions from the other element...
                                        BIGM_CON( :,:, IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)  &
                                            =BIGM_CON( :,:, IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)     + STRESS_IJ_ELE_EXT( :,:, IPHASE, U_SILOC, U_JLOC + Mdims%u_nloc )
                                    END DO
                                END DO
                            END DO
                        ELSE ! on boundary of domain
                            ! bc contribution...
                            ! ADD DIRICHLET BC'S FOR VISCOCITY...
                            DO IPHASE = 1, Mdims%nphase
                                JPHASE = IPHASE
                                DO JDIM=1,Mdims%ndim
                                    IF( WIC_U_BC_ALL_VISC( JDIM, IPHASE, SELE2 ) == WIC_U_BC_DIRICHLET   ) THEN
                                        DO U_SILOC=1,Mdims%u_snloc
                                            U_ILOC   =U_SLOC2LOC(U_SILOC)
                                            DO U_JLOC=1,Mdims%u_nloc
                                                U_JLOC2 = U_JLOC
                                                DO IDIM=1,Mdims%ndim
                                                    DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                                        =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)       + STRESS_IJ_ELE_EXT( IDIM,JDIM, IPHASE, U_SILOC, U_JLOC )
                                                    IF(PIVIT_ON_VISC) THEN
                                                        I = IDIM+(IPHASE-1)*Mdims%ndim+(U_ILOC-1)*Mdims%ndim*Mdims%nphase
                                                        J = JDIM+(JPHASE-1)*Mdims%ndim+(U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                                        Mmat%PIVIT_MAT(I,J,ELE) &
                                                            =Mmat%PIVIT_MAT(I,J,ELE) + ws * STRESS_IJ_ELE_EXT( IDIM, JDIM, IPHASE, U_SILOC, U_JLOC )
                                                    ENDIF
                                                END DO
                                            END DO
                                        END DO
                                    ENDIF
                                END DO
                            END DO
                            ! bc contribution...
                            DO IPHASE = 1, Mdims%nphase
                                JPHASE = IPHASE
                                DO JDIM=1,Mdims%ndim
                                    IF( WIC_U_BC_ALL_VISC( JDIM, IPHASE, SELE2 ) == WIC_U_BC_DIRICHLET ) THEN
                                        DO U_SILOC=1,Mdims%u_snloc
                                            U_ILOC   =U_SLOC2LOC(U_SILOC)
                                            DO U_SJLOC=1,Mdims%u_snloc
                                                U_JLOC   =U_SLOC2LOC(U_SJLOC)
                                                U_JLOC2 = U_JLOC
                                                DO IDIM=1,Mdims%ndim
                                                    ! ADD DIRICHLET BC'S FOR VISCOCITY...
                                                    LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                        - STRESS_IJ_ELE_EXT( IDIM, JDIM, IPHASE, U_SILOC, U_JLOC + Mdims%u_nloc ) * SUF_U_BC_ALL_VISC( JDIM,IPHASE,U_SJLOC + Mdims%u_snloc*(SELE2-1) )
                                                END DO
                                            END DO
                                        END DO
                                    ENDIF
                                END DO
                            END DO
                        ENDIF ! endof IF(ELE2.GT.0) THEN ELSE
                    ! END OF IF(GOT_DIFFUS .AND. LINEAR_HIGHORDER_DIFFUSION) THEN...
                    ENDIF
                    DO U_SILOC=1,Mdims%u_snloc
                        U_ILOC   =U_SLOC2LOC(U_SILOC)
                        DO U_SJLOC=1,Mdims%u_snloc
                            U_JLOC =U_SLOC2LOC(U_SJLOC)
                            IF(SELE2 /= 0) THEN
                                U_JLOC2=U_JLOC
                            ELSE
                                U_JLOC2=U_ILOC_OTHER_SIDE(U_SJLOC)
                            ENDIF
                            ! add diffusion term...
                            DO IPHASE = 1, Mdims%nphase
                                JPHASE = IPHASE
                                DO IDIM = 1, Mdims%ndim
                                    JDIM = IDIM
                                    I=IDIM + (IPHASE-1)*Mdims%ndim + (U_ILOC-1)*Mdims%ndim*Mdims%nphase
                                    J=JDIM + (JPHASE-1)*Mdims%ndim + (U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                    IU_NOD_DIM_PHA = I + (ELE-1)*Mdims%ndim*Mdims%nphase*Mdims%u_nloc
                                    JU_NOD_DIM_PHA = J + (ELE-1)*Mdims%ndim*Mdims%nphase*Mdims%u_nloc
                                    J2=JDIM + (JPHASE-1)*Mdims%ndim + (U_JLOC2-1)*Mdims%ndim*Mdims%nphase
                                    JU2_NOD_DIM_PHA = J2 + (ELE2-1)*Mdims%ndim*Mdims%nphase*Mdims%u_nloc
                                    VLM=0.0
                                    VLM_NEW=0.0
                                    VLM_OLD=0.0
                                    NN_SNDOTQ_IN    = 0.0
                                    NN_SNDOTQ_OUT   = 0.0
                                    NN_SNDOTQOLD_IN = 0.0
                                    NN_SNDOTQOLD_OUT= 0.0
                                    IF(GOT_VIRTUAL_MASS) THEN
                                        CVM_NN_SNDOTQ_IN     = 0.0
                                        CVM_NN_SNDOTQ_OUT    = 0.0
                                        CVM_NN_SNDOTQOLD_IN  = 0.0
                                        CVM_NN_SNDOTQOLD_OUT = 0.0
                                    ENDIF
                                    ! Have a surface integral on element boundary...
                                    DO SGI=1,FE_GIdims%sbcvngi
                                        RNN=SDETWE(SGI)*SBUFEN_REVERSED(SGI,U_SILOC)*SBUFEN_REVERSED(SGI,U_SJLOC)
                                        VLM=VLM+RNN
                                        IF(IDIVID_BY_VOL_FRAC==1) THEN ! We are dividing by vol fract.
                                            VLM_NEW = VLM_NEW + FTHETA( IDIM,IPHASE,SGI ) * RNN &
                                                * DIFF_COEF_DIVDX( IDIM,IPHASE,SGI ) * SVOL_FRA2(IPHASE,SGI)
                                            VLM_OLD = VLM_OLD + (1.-FTHETA( IDIM,IPHASE,SGI )) * RNN &
                                                * DIFF_COEFOLD_DIVDX( IDIM,IPHASE,SGI ) * SVOL_FRA2(IPHASE,SGI)
                                        ELSE
                                            VLM_NEW = VLM_NEW + FTHETA( IDIM,IPHASE,SGI ) * RNN &
                                                * DIFF_COEF_DIVDX( IDIM,IPHASE,SGI )
                                            VLM_OLD = VLM_OLD + (1.-FTHETA( IDIM,IPHASE,SGI )) * RNN &
                                                * DIFF_COEFOLD_DIVDX( IDIM,IPHASE,SGI )
                                        ENDIF
                                        NN_SNDOTQ_IN    = NN_SNDOTQ_IN     + SNDOTQ_IN(IDIM,IPHASE,SGI)    *RNN
                                        NN_SNDOTQ_OUT   = NN_SNDOTQ_OUT    + SNDOTQ_OUT(IDIM,IPHASE,SGI)   *RNN
                                        NN_SNDOTQOLD_IN = NN_SNDOTQOLD_IN  + SNDOTQOLD_IN(IDIM,IPHASE,SGI) *RNN
                                        NN_SNDOTQOLD_OUT= NN_SNDOTQOLD_OUT + SNDOTQOLD_OUT(IDIM,IPHASE,SGI)*RNN
                                        IF(GOT_VIRTUAL_MASS) THEN
                                            DO KPHASE=1,Mdims%nphase
                                                CVM_NN_SNDOTQ_IN(KPHASE)    = CVM_NN_SNDOTQ_IN(KPHASE)     + CVM_SNDOTQ_IN(IDIM,IPHASE,KPHASE,SGI)    *RNN
                                                CVM_NN_SNDOTQ_OUT(KPHASE)   = CVM_NN_SNDOTQ_OUT(KPHASE)    + CVM_SNDOTQ_OUT(IDIM,IPHASE,KPHASE,SGI)   *RNN
                                                CVM_NN_SNDOTQOLD_IN(KPHASE) = CVM_NN_SNDOTQOLD_IN(KPHASE)  + CVM_SNDOTQOLD_IN(IDIM,IPHASE,KPHASE,SGI) *RNN
                                                CVM_NN_SNDOTQOLD_OUT(KPHASE)= CVM_NN_SNDOTQOLD_OUT(KPHASE) + CVM_SNDOTQOLD_OUT(IDIM,IPHASE,KPHASE,SGI)*RNN
                                            END DO
                                        ENDIF
                                    END DO
                                    IF(SELE2 == 0) THEN
                                        IF(Mmat%NO_MATRIX_STORE) THEN
                                            IF(MOM_CONSERV) THEN
                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                    - NN_SNDOTQ_OUT*SLOC_U( IDIM,IPHASE,U_SJLOC ) -  NN_SNDOTQ_IN*SLOC2_U(IDIM,IPHASE,U_SJLOC)
                                            ELSE
                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                    + NN_SNDOTQ_IN*SLOC_U( IDIM,IPHASE,U_SJLOC ) -  NN_SNDOTQ_IN*SLOC2_U(IDIM,IPHASE,U_SJLOC)
                                            ENDIF
                                            ! viscosity...
                                            LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                - VLM_NEW*SLOC_U( IDIM,IPHASE,U_SJLOC )    +  VLM_NEW*SLOC2_U(IDIM,IPHASE,U_SJLOC)
                                            IF(GOT_VIRTUAL_MASS) THEN ! NO matrix free virtual mass term.
                                                STOP 1811
                                            ENDIF
                                        ELSE
                                            IF(MOM_CONSERV) THEN
                                                DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                                    =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)       + NN_SNDOTQ_OUT
                                                BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)  &
                                                    =BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)     + NN_SNDOTQ_IN
                                            ELSE
                                                DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                                    =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)     - NN_SNDOTQ_IN
                                                BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)  &
                                                    =BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)   + NN_SNDOTQ_IN
                                            ENDIF
                                            IF(GOT_VIRTUAL_MASS) THEN
                                                DO KPHASE=1,Mdims%nphase
                                                    DIAG_BIGM_CON(IDIM,JDIM,IPHASE,KPHASE,U_ILOC,U_JLOC,ELE)  &
                                                        =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,KPHASE,U_ILOC,U_JLOC,ELE)       + CVM_BETA*CVM_NN_SNDOTQ_OUT(KPHASE) &
                                                        - (1.-CVM_BETA)*CVM_NN_SNDOTQ_IN(KPHASE)
                                                    BIGM_CON(IDIM,JDIM,IPHASE,KPHASE,U_ILOC,U_JLOC2,COUNT_ELE)  &
                                                        =BIGM_CON(IDIM,JDIM,IPHASE,KPHASE,U_ILOC,U_JLOC2,COUNT_ELE)     + CVM_BETA*CVM_NN_SNDOTQ_IN(KPHASE) &
                                                        + (1.-CVM_BETA)*CVM_NN_SNDOTQ_IN(KPHASE)
                                                END DO
                                            ENDIF
                                            ! viscosity...
                                            DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                                =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)       + VLM_NEW
                                            BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)  &
                                                =BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)     - VLM_NEW
                                        ENDIF
                                        IF(PIVIT_ON_VISC) THEN
                                            I = IDIM+(IPHASE-1)*Mdims%ndim+(U_ILOC-1)*Mdims%ndim*Mdims%nphase
                                            J = JDIM+(JPHASE-1)*Mdims%ndim+(U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                            Mmat%PIVIT_MAT(I,J,ELE) &
                                                =Mmat%PIVIT_MAT(I,J,ELE) + ws * VLM_NEW
                                        ENDIF
                                        RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) = RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) &
                                            - VLM_OLD * SLOC_UOLD( IDIM, IPHASE, U_SJLOC ) + VLM_OLD * SLOC2_UOLD( IDIM, IPHASE, U_SJLOC ) &
                                            - VLM_NEW * SLOC_U( IDIM, IPHASE, U_SJLOC )    + VLM_NEW * SLOC2_U( IDIM, IPHASE, U_SJLOC )
                                        IF(MOM_CONSERV) THEN
                                            LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                -(+NN_SNDOTQOLD_OUT) * SLOC_UOLD( IDIM,IPHASE,U_SJLOC )   -(+NN_SNDOTQOLD_IN) * SLOC2_UOLD( IDIM,IPHASE,U_SJLOC )
                                        ELSE
                                            LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                -(-NN_SNDOTQOLD_IN) * SLOC_UOLD( IDIM,IPHASE,U_SJLOC )  -(+NN_SNDOTQOLD_IN) * SLOC2_UOLD( IDIM,IPHASE,U_SJLOC )
                                        ENDIF
                                        IF(GOT_VIRTUAL_MASS) THEN
                                            DO KPHASE=1,Mdims%nphase
                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                    -CVM_BETA*(+CVM_NN_SNDOTQOLD_OUT(KPHASE)) * SLOC_UOLD( IDIM,KPHASE,U_SJLOC )   -CVM_BETA*(+CVM_NN_SNDOTQOLD_IN(KPHASE)) * SLOC2_UOLD( IDIM,KPHASE,U_SJLOC )  &
                                                    ! non-conservative contribution...
                                                    -(1.-CVM_BETA)*(-CVM_NN_SNDOTQOLD_IN(KPHASE)) * SLOC_UOLD( IDIM,KPHASE,U_SJLOC )  -(1.-CVM_BETA)*(+CVM_NN_SNDOTQOLD_IN(KPHASE)) * SLOC2_UOLD( IDIM,KPHASE,U_SJLOC )
                                            END DO
                                        ENDIF
                                        ! Viscosity...
                                        LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) -VLM_OLD * SLOC_UOLD( IDIM,IPHASE,U_SJLOC )  &
                                            +VLM_OLD * SLOC2_UOLD( IDIM,IPHASE,U_SJLOC )
                                    ELSE
                                        IF( WIC_U_BC_ALL_VISC( IDIM, IPHASE, SELE2 ) == WIC_U_BC_DIRICHLET ) THEN
                                            IF(Mmat%NO_MATRIX_STORE) THEN
                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                    =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) - VLM_NEW * SLOC_U( IDIM,IPHASE,U_SJLOC )
                                            else
                                                DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                                    =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) + VLM_NEW
                                            end if
                                            IF(PIVIT_ON_VISC) THEN
                                                I = IDIM+(IPHASE-1)*Mdims%ndim+(U_ILOC-1)*Mdims%ndim*Mdims%nphase
                                                J = JDIM+(JPHASE-1)*Mdims%ndim+(U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                                Mmat%PIVIT_MAT(I,J,ELE)=Mmat%PIVIT_MAT(I,J,ELE) &
                                                    + ws * VLM_NEW
                                            ENDIF
                                            LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) -VLM_OLD * SLOC_UOLD( IDIM,IPHASE,U_SJLOC )
                                            LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) + (VLM_NEW + VLM_OLD) * SUF_U_BC_ALL_VISC( IDIM,IPHASE,U_SJLOC + Mdims%u_snloc * ( SELE2 - 1 ) )
                                            RHS_DIFF_U( IDIM,IPHASE,U_ILOC ) = RHS_DIFF_U( IDIM,IPHASE,U_ILOC ) - VLM_NEW * SLOC_U( IDIM,IPHASE,U_SJLOC ) &
                                                - VLM_OLD * SLOC_UOLD( IDIM,IPHASE,U_SJLOC ) &
                                                + (VLM_NEW + VLM_OLD) * SUF_U_BC_ALL_VISC( IDIM,IPHASE,U_SJLOC + Mdims%u_snloc * ( SELE2 - 1 ) )
                                        !   DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  + VLM_NEW
                                        ELSE IF( WIC_U_BC_ALL( IDIM, IPHASE, SELE2 ) == WIC_U_BC_ROBIN ) THEN
                                            IF(Mmat%NO_MATRIX_STORE) THEN
                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                    - VLM * SUF_U_ROB1_BC_ALL( IDIM,IPHASE,U_SJLOC + Mdims%u_snloc* ( SELE2 - 1 ) )*SLOC_U( IDIM,IPHASE,U_SJLOC )
                                            ELSE
                                                DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) &
                                                    =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)+ VLM * SUF_U_ROB1_BC_ALL( IDIM,IPHASE,U_SJLOC + Mdims%u_snloc * ( SELE2 - 1 ) )
                                            !  DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  + VLM * SUF_U_ROB1_BC_ALL( IDIM,IPHASE,U_SJLOC + Mdims%u_snloc * ( SELE2 - 1 ) )
                                            ENDIF
                                            IF(PIVIT_ON_VISC) THEN
                                                I = IDIM+(IPHASE-1)*Mdims%ndim+(U_ILOC-1)*Mdims%ndim*Mdims%nphase
                                                J = JDIM+(JPHASE-1)*Mdims%ndim+(U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                                Mmat%PIVIT_MAT(I,J,ELE)=Mmat%PIVIT_MAT(I,J,ELE) &
                                                    + ws * VLM * SUF_U_ROB1_BC_ALL( IDIM,IPHASE,U_SJLOC + Mdims%u_snloc * ( SELE2 - 1 ) )
                                            ENDIF
                                            LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) - VLM * SUF_U_ROB2_BC_ALL( IDIM,IPHASE,U_SJLOC + Mdims%u_snloc * ( SELE2 - 1 ) )
                                            RHS_DIFF_U( IDIM,IPHASE,U_ILOC ) = RHS_DIFF_U( IDIM,IPHASE,U_ILOC ) &
                                                - VLM * SUF_U_ROB1_BC_ALL( IDIM, IPHASE, U_SJLOC + Mdims%u_snloc * (  SELE2 - 1 ) ) * SLOC_U( IDIM, IPHASE, U_SJLOC ) &
                                                - VLM * SUF_U_ROB2_BC_ALL( IDIM, IPHASE, U_SJLOC + Mdims%u_snloc * (  SELE2 - 1 ) )
                                        ENDIF
                                        ! BC for incoming momentum...
                                        IF( WIC_MOMU_BC_ALL( IDIM, IPHASE, SELE2 ) == WIC_U_BC_DIRICHLET ) THEN
                                            IF(MOM_CONSERV) THEN
                                                IF(.NOT.Mmat%NO_MATRIX_STORE) THEN
                                                    DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) &
                                                        =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)+ NN_SNDOTQ_OUT
                                                ELSE
                                                    LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                        - NN_SNDOTQ_OUT * SLOC_U( IDIM,IPHASE,U_SJLOC )
                                                ENDIF
                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                    - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_MOMU_BC_ALL( IDIM,IPHASE,U_SJLOC + Mdims%u_snloc * ( SELE2 - 1 ) ) &
                                                    - NN_SNDOTQOLD_OUT * SLOC_UOLD(IDIM,IPHASE,U_SJLOC)
                                               ! ENDOF IF(MOM_CONSERV) THEN...
                                            ELSE
                                                IF(.NOT.Mmat%NO_MATRIX_STORE) THEN
                                                    DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) &
                                                        =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) - NN_SNDOTQ_IN
                                                !   DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  - NN_SNDOTQ_IN
                                                ELSE
                                                    LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                        + NN_SNDOTQ_IN * SLOC_U( IDIM,IPHASE,U_SJLOC )
                                                ENDIF
                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                    - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_MOMU_BC_ALL( IDIM,IPHASE,U_SJLOC + Mdims%u_snloc * ( SELE2 - 1 ) ) &
                                                    + NN_SNDOTQOLD_IN * SLOC_UOLD(IDIM,IPHASE,U_SJLOC)
                                               ! END OF IF(MOM_CONSERV) THEN ELSE...
                                            ENDIF
                                           ! BC for incoming and outgoing momentum (NO leaking of momentum into or out of domain for example)...
                                        ELSE IF( WIC_MOMU_BC_ALL( IDIM, IPHASE, SELE2 ) == WIC_U_BC_DIRICHLET_INOUT ) THEN
                                            IF(MOM_CONSERV) THEN
                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                    - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN + NN_SNDOTQ_OUT + NN_SNDOTQOLD_OUT)*SUF_MOMU_BC_ALL( IDIM,IPHASE,U_SJLOC + Mdims%u_snloc * ( SELE2 - 1 ) )
                                               ! ENDOF IF(MOM_CONSERV) THEN...
                                            ELSE
                                                IF(.NOT.Mmat%NO_MATRIX_STORE) THEN
                                                    DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) &
                                                        =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) - (NN_SNDOTQ_IN + NN_SNDOTQ_OUT)
                                                !   DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  - (NN_SNDOTQ_IN + NN_SNDOTQ_OUT)
                                                ELSE
                                                    LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                        + (NN_SNDOTQ_IN + NN_SNDOTQ_OUT) * SLOC_U( IDIM,IPHASE,U_SJLOC )
                                                ENDIF
                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                    - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN + NN_SNDOTQ_OUT + NN_SNDOTQOLD_OUT)*SUF_MOMU_BC_ALL( IDIM,IPHASE,U_SJLOC + Mdims%u_snloc * ( SELE2 - 1 ) ) &
                                                    + (NN_SNDOTQOLD_IN + NN_SNDOTQOLD_OUT) * SLOC_UOLD(IDIM,IPHASE,U_SJLOC)
                                               ! END OF IF(MOM_CONSERV) THEN ELSE...
                                            ENDIF
                                           ! END OF IF( WIC_MOMU_BC(SELE2+(IPHASE-1)*Mdims%stotel) == WIC_U_BC_DIRICHLET) THEN ELSE...
                                        ENDIF
                                        IF(GOT_VIRTUAL_MASS) THEN
                                            DO KPHASE=1,Mdims%nphase
                                                IF( WIC_MOMU_BC_ALL( IDIM, KPHASE, SELE2 ) == WIC_U_BC_DIRICHLET ) THEN
                                                    DIAG_BIGM_CON(IDIM,JDIM,IPHASE,KPHASE,U_ILOC,U_JLOC,ELE) &
                                                        =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,KPHASE,U_ILOC,U_JLOC,ELE)+ CVM_BETA*CVM_NN_SNDOTQ_OUT(KPHASE) &
                                                        - (1.-CVM_BETA)*CVM_NN_SNDOTQ_IN(KPHASE)
                                                    LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                        - CVM_BETA*( CVM_NN_SNDOTQ_IN(KPHASE) + CVM_NN_SNDOTQOLD_IN(KPHASE) )*SUF_MOMU_BC_ALL( IDIM,KPHASE,U_SJLOC + Mdims%u_snloc * ( SELE2 - 1 ) ) &
                                                        - CVM_BETA*CVM_NN_SNDOTQOLD_OUT(KPHASE) * SLOC_UOLD(IDIM,KPHASE,U_SJLOC) &
                                                        ! non-conservative form...
                                                        - (1.-CVM_BETA)*( CVM_NN_SNDOTQ_IN(KPHASE) + CVM_NN_SNDOTQOLD_IN(KPHASE) )*SUF_MOMU_BC_ALL( IDIM,KPHASE,U_SJLOC + Mdims%u_snloc * ( SELE2 - 1 ) ) &
                                                        + (1.-CVM_BETA)*CVM_NN_SNDOTQOLD_IN(KPHASE) * SLOC_UOLD(IDIM,KPHASE,U_SJLOC)
                                                   ! BC for incoming and outgoing momentum (NO leaking of momentum into or out of domain for example)...
                                                ELSE IF( WIC_MOMU_BC_ALL( IDIM, KPHASE, SELE2 ) == WIC_U_BC_DIRICHLET_INOUT ) THEN
                                                    LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                        - CVM_BETA*( CVM_NN_SNDOTQ_IN(KPHASE) + CVM_NN_SNDOTQOLD_IN(KPHASE) + CVM_NN_SNDOTQ_OUT(KPHASE) + CVM_NN_SNDOTQOLD_OUT(KPHASE))*SUF_MOMU_BC_ALL( IDIM,KPHASE,U_SJLOC + Mdims%u_snloc * ( SELE2 - 1 ) )  &
                                                        ! non-conservative form...
                                                        - (1.-CVM_BETA)*(CVM_NN_SNDOTQ_IN(KPHASE) + CVM_NN_SNDOTQOLD_IN(KPHASE) + CVM_NN_SNDOTQ_OUT(KPHASE) + CVM_NN_SNDOTQOLD_OUT(KPHASE))*SUF_MOMU_BC_ALL( IDIM,KPHASE,U_SJLOC + Mdims%u_snloc * ( SELE2 - 1 ) ) &
                                                        + (1.-CVM_BETA)*(CVM_NN_SNDOTQOLD_IN(KPHASE) + CVM_NN_SNDOTQOLD_OUT(KPHASE)) * SLOC_UOLD(IDIM,KPHASE,U_SJLOC)
                                                    DIAG_BIGM_CON(IDIM,JDIM,IPHASE,KPHASE,U_ILOC,U_JLOC,ELE) &
                                                        =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,KPHASE,U_ILOC,U_JLOC,ELE) - (1.-CVM_BETA)*(CVM_NN_SNDOTQ_IN(KPHASE) + CVM_NN_SNDOTQ_OUT(KPHASE))
                                                   ! END OF IF( WIC_MOMU_BC(SELE2+(KPHASE-1)*Mdims%stotel) == WIC_U_BC_DIRICHLET) THEN ELSE...
                                                ENDIF
                                            END DO ! ENDOF DO KPHASE=1,Mdims%nphase
                                        ENDIF ! ENDOF IF(GOT_VIRTUAL_MASS) THEN
                                    ENDIF
                                END DO
                            END DO
                        END DO
                    END DO
                ENDIF If_diffusion_or_momentum3
            END DO Between_Elements_And_Boundary
            !      END DO Loop_Elements2
            !! *************************end loop over surfaces*********************************************
            ! ideally insert inner element stabilization here...
            ! copy local memory
            DO U_ILOC = 1, Mdims%u_nloc
                U_INOD = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                Mmat%U_RHS( :, :, U_INOD ) = Mmat%U_RHS( :, :, U_INOD ) + LOC_U_RHS( :, :, U_ILOC )
            END DO
           !      END DO Loop_Elements
        END DO Loop_Elements2
        ! This subroutine combines the distributed and block diagonal for an element
        ! into the matrix DGM_PHA.
        IF(.NOT.Mmat%NO_MATRIX_STORE) THEN
            CALL COMB_VEL_MATRIX_DIAG_DIST(DIAG_BIGM_CON, BIGM_CON, &
                Mmat%DGM_petsc, &
                Mspars%ELE%fin, Mspars%ELE%col, Mdims%ndim, Mdims%nphase, Mdims%u_nloc, Mdims%totele, velocity, pressure)  ! Element connectivity.
            DEALLOCATE( DIAG_BIGM_CON )
            DEALLOCATE( BIGM_CON)
        ENDIF
        DEALLOCATE( UD, UD_ND )
        DEALLOCATE( UDOLD, UDOLD_ND )
        DEALLOCATE( DENGI )
        DEALLOCATE( DENGIOLD )
        DEALLOCATE( GRAD_SOU_GI )
        DEALLOCATE( SIGMAGI )
        DEALLOCATE( SIGMAGI_STAB )
        DEALLOCATE( NN_SIGMAGI_ELE )
        DEALLOCATE( NN_SIGMAGI_STAB_ELE )
        DEALLOCATE( NN_MASS_ELE )
        DEALLOCATE( NN_MASSOLD_ELE )
        DEALLOCATE( CV_SLOC2LOC )
        DEALLOCATE( U_SLOC2LOC )
        DEALLOCATE(SDETWE)
        DEALLOCATE( U_ILOC_OTHER_SIDE )
        DEALLOCATE(STORED_U_ILOC_OTHER_SIDE)
        DEALLOCATE(STORED_U_OTHER_LOC)
        DEALLOCATE(STORED_MAT_OTHER_LOC)
        DEALLOCATE( U_OTHER_LOC )
        DEALLOCATE( MAT_OTHER_LOC )
        DEALLOCATE( TEN_XX )
        DEALLOCATE( VLN )
        DEALLOCATE( VLN_OLD )
        DEALLOCATE( STRESS_IJ_ELE )
        DEALLOCATE( VLK_ELE )
        DEALLOCATE( SUD_ALL )
        DEALLOCATE( SUDOLD_ALL )
        DEALLOCATE( SUD2_ALL )
        DEALLOCATE( SUDOLD2_ALL )
        DEALLOCATE( SNDOTQ )
        DEALLOCATE( SNDOTQOLD )
        DEALLOCATE( SINCOME )
        DEALLOCATE( SINCOMEOLD )
        DEALLOCATE( SDEN )
        DEALLOCATE( SDENOLD )
        DEALLOCATE( SDEN_KEEP )
        DEALLOCATE( SDENOLD_KEEP )
        DEALLOCATE( SDEN2_KEEP )
        DEALLOCATE( SDENOLD2_KEEP )
        DEALLOCATE( DIFF_COEF_DIVDX )
        DEALLOCATE( DIFF_COEFOLD_DIVDX )
        DEALLOCATE( FTHETA )
        DEALLOCATE( SNDOTQ_IN )
        DEALLOCATE( SNDOTQ_OUT )
        DEALLOCATE( SNDOTQOLD_IN )
        DEALLOCATE( SNDOTQOLD_OUT )
        DEALLOCATE( GRAD_SOU_GI_NMX )
        DEALLOCATE( MASS_ELE )
        DEALLOCATE( FACE_ELE )
        ! Deallocating for non-linear Petrov-Galerkin diffusion stabilization...
        DEALLOCATE( LOC_MASS_INV )
        DEALLOCATE( LOC_MASS )
        DEALLOCATE( RHS_DIFF_U )
        DEALLOCATE( DIFF_VEC_U )
        DEALLOCATE( DIFFGI_U )
        DEALLOCATE( U_DT )
        DEALLOCATE( SOUGI_X )
        DEALLOCATE( RESID_U)
        DEALLOCATE( P_DX )
        DEALLOCATE( U_GRAD_NORM2, U_GRAD_NORM )
        DEALLOCATE( A_DOT_U )
        DEALLOCATE( STAR_U_COEF )
        DEALLOCATE( P_STAR_U )
        DEALLOCATE( DIF_STAB_U )
        DEALLOCATE( UDIFF_SUF_STAB )
        DEALLOCATE( VLK_UVW )
        ! reversed indicies for shape functions...
        DEALLOCATE( CVFENX_ALL_REVERSED, UFENX_ALL_REVERSED )
        DEALLOCATE( UFEN_REVERSED, CVN_REVERSED, CVFEN_REVERSED )
        DEALLOCATE( SBCVFEN_REVERSED, SBUFEN_REVERSED )
        call deallocate(velocity_BCs)
        call deallocate(velocity_BCs_visc)
        call deallocate(velocity_BCs_adv)
        call deallocate(velocity_BCs_robin2)
        call deallocate(momentum_BCs)
        call deallocate(pressure_BCs)
        call deallocate_multi_dev_shape_funs(Devfuns)
        ewrite(3,*)'Leaving assemb_force_cty'
        RETURN


        contains


        subroutine get_porous_Mass_matrix(ELE, Mdims, FE_funs, FE_GIdims, DevFuns, Mmat, WIC_P_BC_ALL, FACE_ELE)
            implicit none
            integer, intent(in) :: ELE
            type(multi_dimensions), intent(in) :: Mdims
            type(multi_GI_dimensions), intent(in) :: FE_GIdims
            type(multi_shape_funs), intent(in) :: FE_funs
            type(multi_dev_shape_funs), intent(in) :: Devfuns
            type (multi_matrices), intent(inout) :: Mmat
            INTEGER, DIMENSION ( :, :, :), intent(in) :: WIC_P_BC_ALL
            INTEGER, DIMENSION( :, : ), intent(in) ::  FACE_ELE
            !Local variables
            integer:: sele, sele2, ele2, ele3, iface, iface2, I, J, vel_degree, pres_degree,&
                    U_JLOC, U_ILOC, JPHASE, JDIM, JPHA_JDIM, IPHASE, IPHA_IDIM
            logical :: P_BC_element
            !Weight parameter, controls the strenght of the homogenization of the velocity nodes per element
            !the bigger the more P0DG it tends to be
            real :: factor, factor_default, bc_factor
            real, save :: lump_vol_factor =-1d25
            real, save :: scaling_vel_nodes = -1
            !Weights for lumping taken from Zienkiewicz vol 1 page 475
            real, parameter :: corner = 3./57.
            real, parameter :: midpoint = 16./57.
            !Obtain the scaling factor to spread the volume of the mass matrix
            if (scaling_vel_nodes<0) then
                scaling_vel_nodes = dble(Mdims%u_nloc)
                !Adjust for linear bubble functions, P1(BL)DG
                !We are adding an extra node that adds extra velocity that needs to be compensated
                if ((Mdims%ndim==2 .and. Mdims%u_nloc == 4) .or.&
                        (Mdims%ndim==3 .and. Mdims%u_nloc == 5)) then
                    scaling_vel_nodes = scaling_vel_nodes - 1.
                    lump_vol_factor = 0.!No homogenisation for bubble functions
                end if
            end if

            select case (Mdims%u_nloc)
                case (6) !Quadratic 2D
                    DO U_JLOC = 1, Mdims%u_nloc
                        DO JPHASE = 1, Mdims%nphase
                            DO JDIM = 1, Mdims%ndim
                                J = JDIM+(JPHASE-1)*Mdims%ndim+(U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                select case (U_JLOC)
                                    case (1,3,6)
                                        Mmat%PIVIT_MAT( J, J, ELE ) = DevFuns%volume * corner
                                    case default
                                        Mmat%PIVIT_MAT( J, J, ELE ) = DevFuns%volume * midpoint
                                end select
                            end do
                        end do
                    end do
                case default !Create the mass matrix normally by distributting the mass evenly between the nodes
                    do i=1,size(Mmat%PIVIT_MAT,1)
                        Mmat%PIVIT_MAT(I,I,ELE) = DevFuns%VOLUME/scaling_vel_nodes
                    END DO
            end select


!            Porous_media_PIVIT_not_stored_yet = .false.
            !If pressure boundary element, then we homogenize the velocity in the element
            if (lump_vol_factor<-1d24) then
                call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/polynomial_degree', vel_degree )
                call get_option( '/geometry/mesh::PressureMesh/from_mesh/mesh_shape/polynomial_degree', pres_degree )
                if (max(pres_degree,vel_degree)>1) then
                    factor_default = 0.
                else
                    factor_default = 100.
                end if
                !Obtain the value from diamond
                call get_option( '/numerical_methods/CV_press_homogenisation', factor, default = factor_default )
                !This value is the amount of mass used to homogenize the element
                lump_vol_factor = factor * DevFuns%VOLUME/dble(Mdims%u_nloc)
            end if

            !If CV_press_homogenisation negative or zero, then, do not apply this method
            if (lump_vol_factor <= 0.) return

            !No coupling between dimensions nor phases, only based on geometry
            DO U_JLOC = 1, Mdims%u_nloc
                DO U_ILOC = 1, Mdims%u_nloc
                    DO IPHASE = 1, Mdims%nphase
                        DO IDIM = 1, Mdims%ndim
                            J = IDIM+(IPHASE-1)*Mdims%ndim+(U_JLOC-1)*Mdims%ndim*Mdims%nphase
                            I = IDIM+(IPHASE-1)*Mdims%ndim+(U_ILOC-1)*Mdims%ndim*Mdims%nphase
                            Mmat%PIVIT_MAT(I,I,ELE) = Mmat%PIVIT_MAT(I,I,ELE) + lump_vol_factor
                            Mmat%PIVIT_MAT(I,J,ELE) = Mmat%PIVIT_MAT(I,J,ELE) - lump_vol_factor
                        end do
                    end do
                end do
            end do

        end subroutine get_porous_Mass_matrix


    END SUBROUTINE ASSEMB_FORCE_CTY


    SUBROUTINE DG_VISC_LIN( S_INV_NNX_MAT12, NNX_MAT, NNX_MAT2, NN_MAT, NN_MAT2,  &
        U_SNLOC, U_NLOC, SBUFEN_REVERSED, SDETWE, SBCVNGI, SNORMXN_ALL, NDIM, &
        U_SLOC2LOC, U_OTHER_LOC, U_NLOC_EXT, ON_BOUNDARY )

        ! This sub calculates S_INV_NNX_MAT12 which contains NDIM matricies that are used to form the
        ! inter element coupling for the viscocity discretization.
        ! NNX_MAT, NNX_MAT2 contain matricies of first derivative times basis function for current element and neightbouring element.
        ! Similarly NN_MAT, NN_MAT2, contain the element-wise mass matrices.
        ! Only call this sub if element not next to the boundary...

        INTEGER, intent( in ) :: U_SNLOC, U_NLOC, SBCVNGI, NDIM, U_NLOC_EXT
        LOGICAL, intent( in ) :: ON_BOUNDARY

        !          INTEGER, PARAMETER :: U_NLOC_EXT = U_NLOC*2 - U_SNLOC
        REAL, DIMENSION( NDIM, U_SNLOC, 2*U_NLOC ), intent( inout ) :: S_INV_NNX_MAT12
        REAL, DIMENSION( NDIM, U_NLOC, U_NLOC  ), intent( in ) :: NNX_MAT, NNX_MAT2
        REAL, DIMENSION( U_NLOC, U_NLOC  ), intent( in ) :: NN_MAT, NN_MAT2

        REAL, DIMENSION( SBCVNGI,  U_SNLOC ), intent( in ) :: SBUFEN_REVERSED
        REAL, DIMENSION( NDIM, SBCVNGI  ), intent( in ) :: SNORMXN_ALL
        REAL, DIMENSION( SBCVNGI ), intent( in ) :: SDETWE
        INTEGER, DIMENSION( U_SNLOC ), intent( in ) :: U_SLOC2LOC
        INTEGER, DIMENSION( U_NLOC ), intent( in ) :: U_OTHER_LOC
        !
        ! local variables...
        !          REAL :: NNX_MAT_SUF(NDIM, U_NLOC, U_NLOC), NNX_MAT2_SUF(NDIM, U_NLOC, U_NLOC)
        REAL :: NN_MAT12_INV(U_NLOC_EXT, U_NLOC_EXT)
        REAL :: NN_MAT12(U_NLOC_EXT, U_NLOC_EXT), NNX_MAT12(NDIM, U_NLOC_EXT, 2*U_NLOC)
        REAL :: INV_NNX_MAT12(NDIM, U_NLOC_EXT, 2*U_NLOC)
        REAL :: R_SUF_SUM
        !          INTEGER :: ELE2_LOC_GL_NODS(U_NLOC)
        INTEGER IDIM, U_SILOC, U_ILOC, U_SJLOC, U_JLOC, U_ILOC2, U_JLOC2, I
        !
        !            print *,'just inside DG_VISC_LIN ON_BOUNDARY:',ON_BOUNDARY


        NNX_MAT12  = 0.0
        NN_MAT12   = 0.0

        NNX_MAT12(1:NDIM,  1:U_NLOC, 1:U_NLOC)  = NNX_MAT(1:NDIM, 1:U_NLOC, 1:U_NLOC)
        NN_MAT12(1:U_NLOC, 1:U_NLOC)            = NN_MAT(1:U_NLOC, 1:U_NLOC)

        ! Calculate ELE2_LOC_GL_NODS:
        IF(.NOT.ON_BOUNDARY) THEN

            NNX_MAT12(1:NDIM,  1+U_NLOC:2*U_NLOC, 1+U_NLOC:2*U_NLOC)  = NNX_MAT2(1:NDIM, 1:U_NLOC, 1:U_NLOC)
            NN_MAT12(1+U_NLOC:2*U_NLOC, 1+U_NLOC:2*U_NLOC)  = NN_MAT2(1:U_NLOC, 1:U_NLOC)

            ! put surface contributions into NNX_MAT,NNX_MAT2 and call the result NNX_MAT_SUF,NNX_MAT2_SUF ...
            !
            DO U_SILOC=1,U_SNLOC
                U_ILOC=U_SLOC2LOC(U_SILOC)
                U_ILOC2=U_OTHER_LOC(U_ILOC)
                DO U_SJLOC=1,U_SNLOC
                    U_JLOC=U_SLOC2LOC(U_SJLOC)
                    U_JLOC2=U_OTHER_LOC(U_JLOC)
                    DO IDIM=1,NDIM
                        R_SUF_SUM=SUM( SNORMXN_ALL(IDIM,:)*SBUFEN_REVERSED(:,U_SILOC)*SBUFEN_REVERSED(:,U_SJLOC)*SDETWE(:) )
                        ! ELEside of the element face...
                        NNX_MAT12(IDIM,U_ILOC,U_JLOC)         = &
                            &NNX_MAT12(IDIM,U_ILOC,U_JLOC)          -  0.5* R_SUF_SUM ! we have*2 the contribution
                        NNX_MAT12(IDIM,U_ILOC,U_JLOC2+U_NLOC) = &
                            &NNX_MAT12(IDIM,U_ILOC,U_JLOC2+U_NLOC)  +  0.5* R_SUF_SUM ! to take into account both sides of shared ele face

                        NNX_MAT12(IDIM,U_ILOC2+U_NLOC,U_JLOC2+U_NLOC) = &
                            &NNX_MAT12(IDIM,U_ILOC2+U_NLOC,U_JLOC2+U_NLOC)  +  0.5* R_SUF_SUM ! we have*2 the contribution
                        NNX_MAT12(IDIM,U_ILOC2+U_NLOC,U_JLOC)         = &
                            &NNX_MAT12(IDIM,U_ILOC2+U_NLOC,U_JLOC)          -  0.5* R_SUF_SUM ! to take into account both sides of shared ele face
                    END DO
                END DO
            END DO

            NN_MAT12_INV = NN_MAT12
            CALL INVERT(NN_MAT12_INV)
            DO IDIM=1,NDIM
                INV_NNX_MAT12(IDIM,:,:) = MATMUL( NN_MAT12_INV(:,:), NNX_MAT12(IDIM,:,:) )
            END DO
            S_INV_NNX_MAT12=0.0

            DO U_SILOC=1,U_SNLOC
                U_ILOC=U_SLOC2LOC(U_SILOC)
                U_ILOC2=U_OTHER_LOC(U_ILOC)
                S_INV_NNX_MAT12( 1:NDIM, U_SILOC, : )   = 0.5*( INV_NNX_MAT12( 1:NDIM, U_ILOC, : ) + INV_NNX_MAT12( 1:NDIM, U_ILOC2+U_NLOC, : ) )
            END DO

        ELSE ! IF(.NOT.ON_BOUNDARY) THEN
            DO U_SILOC=1,U_SNLOC
                U_ILOC=U_SLOC2LOC(U_SILOC)
                DO U_SJLOC=1,U_SNLOC
                    U_JLOC=U_SLOC2LOC(U_SJLOC)
                    DO IDIM=1,NDIM
                        R_SUF_SUM=SUM( SNORMXN_ALL(IDIM,:)*SBUFEN_REVERSED(:,U_SILOC)*SBUFEN_REVERSED(:,U_SJLOC)*SDETWE(:) )
                        NNX_MAT12(IDIM,U_ILOC,U_JLOC)          = NNX_MAT12(IDIM,U_ILOC,U_JLOC)           -  0.5 * R_SUF_SUM
                        NNX_MAT12(IDIM,U_ILOC,U_JLOC+U_NLOC)   = NNX_MAT12(IDIM,U_ILOC,U_JLOC+U_NLOC)    +  0.5 * R_SUF_SUM
                    !                        NNX_MAT12(IDIM,U_ILOC,U_JLOC)          = NNX_MAT12(IDIM,U_ILOC,U_JLOC)           -  1.0 * R_SUF_SUM
                    !                        NNX_MAT12(IDIM,U_ILOC,U_JLOC+U_NLOC)   = NNX_MAT12(IDIM,U_ILOC,U_JLOC+U_NLOC)    +  1.0 * R_SUF_SUM
                    END DO
                END DO
            END DO

            DO I=1+U_NLOC,U_NLOC_EXT
                NN_MAT12(I,I)=1.0
            END DO

            NN_MAT12_INV = NN_MAT12
            CALL INVERT(NN_MAT12_INV)
            DO IDIM=1,NDIM
                INV_NNX_MAT12(IDIM,:,:) = MATMUL( NN_MAT12_INV(:,:), NNX_MAT12(IDIM,:,:) )
            END DO
            S_INV_NNX_MAT12=0.0

            DO U_SILOC=1,U_SNLOC
                U_ILOC=U_SLOC2LOC(U_SILOC)
                S_INV_NNX_MAT12( 1:NDIM, U_SILOC, : )   =INV_NNX_MAT12( 1:NDIM, U_ILOC, : )
            END DO

        ENDIF ! IF(.NOT.ON_BOUNDARY) THEN ELSE


        RETURN
    END SUBROUTINE DG_VISC_LIN







    SUBROUTINE VISCOCITY_TENSOR_LES_CALC(LES_UDIFFUSION, LES_UDIFFUSION_VOL, DUX_ELE_ALL, &
        NDIM,NPHASE, U_NLOC,X_NLOC,TOTELE, X_NONODS, &
        X_ALL, X_NDGLN,  MAT_NONODS, MAT_NLOC, MAT_NDGLN, LES_DISOPT, LES_CS, UDEN, CV_NONODS, CV_NDGLN, &
        U_NDGLN, U_NONODS, U_ALL, DERIV )
        ! This subroutine calculates a tensor of viscocity LES_UDIFFUSION, LES_UDIFFUSION_VOL
        IMPLICIT NONE
        REAL, intent( in ) :: LES_CS
        INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, TOTELE, X_NONODS, MAT_NONODS, CV_NONODS, MAT_NLOC, LES_DISOPT, U_NONODS
        INTEGER, DIMENSION( X_NLOC * TOTELE  ), intent( in ) :: X_NDGLN
        INTEGER, DIMENSION( MAT_NLOC * TOTELE  ), intent( in ) :: MAT_NDGLN, CV_NDGLN
        INTEGER, DIMENSION( U_NLOC * TOTELE  ), intent( in ) :: U_NDGLN
        REAL, DIMENSION( NDIM, X_NONODS  ), intent( in ) :: X_ALL
        REAL, DIMENSION( NDIM, NDIM, NPHASE, MAT_NONODS  ), intent( inout ) :: LES_UDIFFUSION
        REAL, DIMENSION( NPHASE, MAT_NONODS  ), intent( inout ) :: LES_UDIFFUSION_VOL
        REAL, DIMENSION( NDIM, NDIM, NPHASE, U_NLOC, TOTELE  ), intent( in ) :: DUX_ELE_ALL
        REAL, DIMENSION( NPHASE, CV_NONODS ), intent( in ) :: UDEN, DERIV
        REAL, DIMENSION( NDIM, NPHASE, U_NONODS  ), intent( in ) :: U_ALL
        ! Local variables...
        !            INTEGER, PARAMETER :: LES_DISOPT=1
        ! LES_DISOPT is LES option e.g. =1 Anisotropic element length scale
        !                               =2 Take the average length scale h
        !                               =3 Take the min length scale h
        !                               =4 Take the max length scale h
        !                               =5 Use different length scales for u,v,w across an element and map visc to stress form
        !                               =6 same as 5 plus original q-scheme that is multiply by -min(0.0,divq) (5 can be switched off by using a zero coefficient for the LES)
        !                               =7 same as 5 plus original q-scheme but using abs(divu) (5 can be switched off by using a zero coefficient for the LES)
        ! =8  same as 6, but added  - \rho C* C_L * SQRT(H2) * MIN(0.0, SIGN(1.0, DIVU) ) into \mu_vol
        ! =9  same as 7, but added  - \rho C* C_L * SQRT(H2) * MIN(0.0, SIGN(1.0, DIVU) ) into \mu_vol
        ! =10 same as 6, but added  + \rho C* C_L * SQRT(H2)  into \mu_vol
        ! =11 same as 7, but added  + \rho C* C_L * SQRT(H2)  into \mu_vol
        integer :: ele, MAT_iloc, MAT_INOD, CV_INOD, iphase
        real, dimension( :, :, :, :, : ), allocatable :: LES_U_UDIFFUSION, LES_MAT_UDIFFUSION
        real, dimension( :, :, : ), allocatable :: LES_U_UDIFFUSION_VOL, LES_MAT_UDIFFUSION_VOL,  Q_SCHEME_ABS_CONT_VOL, SOUND_SPEED
        integer, dimension( : ), allocatable :: NOD_COUNT

        ALLOCATE(LES_U_UDIFFUSION(NDIM,NDIM,NPHASE,U_NLOC,TOTELE))
        ALLOCATE(LES_MAT_UDIFFUSION(NDIM,NDIM,NPHASE,MAT_NLOC,TOTELE))

        ALLOCATE(LES_U_UDIFFUSION_VOL(NPHASE,U_NLOC,TOTELE))
        ALLOCATE(LES_MAT_UDIFFUSION_VOL(NPHASE,MAT_NLOC,TOTELE))

        ALLOCATE(Q_SCHEME_ABS_CONT_VOL(NPHASE,U_NLOC,TOTELE))
        ALLOCATE(SOUND_SPEED(NPHASE,MAT_NLOC,TOTELE))

        ALLOCATE(NOD_COUNT(MAT_NONODS))

        IF(LES_DISOPT.GE.8) THEN
            DO ELE=1,TOTELE
                DO MAT_ILOC=1,MAT_NLOC
                    CV_INOD=CV_NDGLN((ELE-1)*MAT_NLOC + MAT_ILOC)
                    DO IPHASE=1,NPHASE
                        SOUND_SPEED(IPHASE,MAT_ILOC,ELE) = 1.0/MAX( SQRT(MAX(DERIV(IPHASE,CV_INOD),0.0)), 1.0E-10 )
                    END DO
                END DO
            END DO
        ELSE
            SOUND_SPEED=0.0
        END IF

        CALL VISCOCITY_TENSOR_LES_CALC_U( LES_U_UDIFFUSION, LES_U_UDIFFUSION_VOL, Q_SCHEME_ABS_CONT_VOL, &
            DUX_ELE_ALL, NDIM,NPHASE, U_NLOC,X_NLOC,TOTELE, X_NONODS, &
            X_ALL, X_NDGLN, LES_DISOPT, LES_CS,  U_NDGLN, U_NONODS, U_ALL )

        IF(MAT_NLOC==U_NLOC) THEN
            LES_MAT_UDIFFUSION=LES_U_UDIFFUSION
            LES_MAT_UDIFFUSION_VOL=LES_U_UDIFFUSION_VOL + SOUND_SPEED * Q_SCHEME_ABS_CONT_VOL
        ELSE IF( (U_NLOC==3.AND.MAT_NLOC==6) .OR. (U_NLOC==4.AND.MAT_NLOC==10) ) THEN
            !
            LES_MAT_UDIFFUSION(:,:,:,1,:) = LES_U_UDIFFUSION(:,:,:,1,:)
            LES_MAT_UDIFFUSION_VOL(:,1,:) = LES_U_UDIFFUSION_VOL(:,1,:)  + SOUND_SPEED(:,1,:) * Q_SCHEME_ABS_CONT_VOL(:,1,:)
            LES_MAT_UDIFFUSION(:,:,:,3,:) = LES_U_UDIFFUSION(:,:,:,2,:)
            LES_MAT_UDIFFUSION_VOL(:,3,:) = LES_U_UDIFFUSION_VOL(:,2,:)  + SOUND_SPEED(:,3,:) * Q_SCHEME_ABS_CONT_VOL(:,2,:)
            LES_MAT_UDIFFUSION(:,:,:,6,:) = LES_U_UDIFFUSION(:,:,:,3,:)
            LES_MAT_UDIFFUSION_VOL(:,6,:) = LES_U_UDIFFUSION_VOL(:,3,:)  + SOUND_SPEED(:,6,:) * Q_SCHEME_ABS_CONT_VOL(:,3,:)

            LES_MAT_UDIFFUSION(:,:,:,2,:) = 0.5 * ( LES_U_UDIFFUSION(:,:,:,1,:) + LES_U_UDIFFUSION(:,:,:,2,:) )
            LES_MAT_UDIFFUSION_VOL(:,2,:) = 0.5 * ( LES_U_UDIFFUSION_VOL(:,1,:) + LES_U_UDIFFUSION_VOL(:,2,:) )  &
                + SOUND_SPEED(:,2,:) * 0.5 * ( Q_SCHEME_ABS_CONT_VOL(:,1,:)+Q_SCHEME_ABS_CONT_VOL(:,2,:) )
            LES_MAT_UDIFFUSION(:,:,:,4,:) = 0.5 * ( LES_U_UDIFFUSION(:,:,:,1,:) + LES_U_UDIFFUSION(:,:,:,3,:) )
            LES_MAT_UDIFFUSION_VOL(:,4,:) = 0.5 * ( LES_U_UDIFFUSION_VOL(:,1,:) + LES_U_UDIFFUSION_VOL(:,3,:) )  &
                + SOUND_SPEED(:,4,:) * 0.5 * ( Q_SCHEME_ABS_CONT_VOL(:,1,:)+Q_SCHEME_ABS_CONT_VOL(:,3,:) )
            LES_MAT_UDIFFUSION(:,:,:,5,:) = 0.5 * ( LES_U_UDIFFUSION(:,:,:,2,:) + LES_U_UDIFFUSION(:,:,:,3,:) )
            LES_MAT_UDIFFUSION_VOL(:,5,:) = 0.5 * ( LES_U_UDIFFUSION_VOL(:,2,:) + LES_U_UDIFFUSION_VOL(:,3,:) )  &
                + SOUND_SPEED(:,5,:) * 0.5 * ( Q_SCHEME_ABS_CONT_VOL(:,2,:)+Q_SCHEME_ABS_CONT_VOL(:,3,:) )

            if( MAT_NLOC == 10 ) then
                LES_MAT_UDIFFUSION(:,:,:,7,:) = 0.5 * ( LES_U_UDIFFUSION(:,:,:,1,:) + LES_U_UDIFFUSION(:,:,:,4,:) )
                LES_MAT_UDIFFUSION_VOL(:,7,:) = 0.5 * ( LES_U_UDIFFUSION_VOL(:,1,:) + LES_U_UDIFFUSION_VOL(:,4,:) ) &
                    + SOUND_SPEED(:,7,:) * 0.5 * ( Q_SCHEME_ABS_CONT_VOL(:,1,:)+Q_SCHEME_ABS_CONT_VOL(:,4,:) )
                LES_MAT_UDIFFUSION(:,:,:,8,:) = 0.5 * ( LES_U_UDIFFUSION(:,:,:,2,:) + LES_U_UDIFFUSION(:,:,:,4,:) )
                LES_MAT_UDIFFUSION_VOL(:,8,:) = 0.5 * ( LES_U_UDIFFUSION_VOL(:,2,:) + LES_U_UDIFFUSION_VOL(:,4,:) ) &
                    + SOUND_SPEED(:,8,:) * 0.5 * ( Q_SCHEME_ABS_CONT_VOL(:,2,:)+Q_SCHEME_ABS_CONT_VOL(:,4,:) )
                LES_MAT_UDIFFUSION(:,:,:,9,:) = 0.5 * ( LES_U_UDIFFUSION(:,:,:,3,:) + LES_U_UDIFFUSION(:,:,:,4,:) )
                LES_MAT_UDIFFUSION_VOL(:,9,:) = 0.5 * ( LES_U_UDIFFUSION_VOL(:,3,:) + LES_U_UDIFFUSION_VOL(:,4,:) ) &
                    + SOUND_SPEED(:,9,:) * 0.5 * ( Q_SCHEME_ABS_CONT_VOL(:,3,:)+Q_SCHEME_ABS_CONT_VOL(:,4,:) )

                LES_MAT_UDIFFUSION(:,:,:,10,:) = LES_U_UDIFFUSION(:,:,:,4,:)
                LES_MAT_UDIFFUSION_VOL(:,10,:) = LES_U_UDIFFUSION_VOL(:,4,:) + SOUND_SPEED(:,10,:) * Q_SCHEME_ABS_CONT_VOL(:,4,:)
            end if

        ELSE IF( (U_NLOC==6.AND.MAT_NLOC==3) .OR. (U_NLOC==10.AND.MAT_NLOC==4) ) THEN
            LES_MAT_UDIFFUSION(:,:,:,1,:) = LES_U_UDIFFUSION(:,:,:,1,:)
            LES_MAT_UDIFFUSION_VOL(:,1,:) = LES_U_UDIFFUSION_VOL(:,1,:)+ SOUND_SPEED(:,1,:) * Q_SCHEME_ABS_CONT_VOL(:,1,:)
            LES_MAT_UDIFFUSION(:,:,:,2,:) = LES_U_UDIFFUSION(:,:,:,3,:)
            LES_MAT_UDIFFUSION_VOL(:,2,:) = LES_U_UDIFFUSION_VOL(:,3,:)+ SOUND_SPEED(:,2,:) * Q_SCHEME_ABS_CONT_VOL(:,3,:)
            LES_MAT_UDIFFUSION(:,:,:,3,:) = LES_U_UDIFFUSION(:,:,:,6,:)
            LES_MAT_UDIFFUSION_VOL(:,3,:) = LES_U_UDIFFUSION_VOL(:,6,:)+ SOUND_SPEED(:,3,:) * Q_SCHEME_ABS_CONT_VOL(:,6,:)
            if( u_nloc == 10 ) then
                LES_MAT_UDIFFUSION(:,:,:,4,:) = LES_U_UDIFFUSION(:,:,:,10,:)
                LES_MAT_UDIFFUSION_VOL(:,4,:) = LES_U_UDIFFUSION_VOL(:,10,:)+ SOUND_SPEED(:,4,:) * Q_SCHEME_ABS_CONT_VOL(:,10,:)
            end if
        ELSE
            PRINT *,'not ready to onvert between these elements'
            STOP 2211
        ENDIF


        ! Now map to nodal variables from element variables...
        NOD_COUNT=0
        LES_UDIFFUSION=0.0
        LES_UDIFFUSION_VOL=0.0
        do ele = 1, totele
            do MAT_iloc = 1, MAT_nloc
                MAT_Inod = MAT_ndgln( ( ele - 1 ) * MAT_nloc + MAT_iloc )
                CV_Inod = CV_ndgln( ( ele - 1 ) * MAT_nloc + MAT_iloc )

                do iphase = 1, nphase
                    LES_UDIFFUSION(:,:,iphase,MAT_Inod) = LES_UDIFFUSION(:,:,iphase,MAT_Inod) + &
                        LES_MAT_UDIFFUSION(:,:,iphase,MAT_ILOC,ELE) * UDEN( iphase, CV_Inod )
                    LES_UDIFFUSION_VOL(iphase,MAT_Inod) = LES_UDIFFUSION_VOL(iphase,MAT_Inod) + &
                        LES_MAT_UDIFFUSION_VOL(iphase,MAT_ILOC,ELE) * UDEN( iphase, CV_Inod )
                end do
                NOD_COUNT(MAT_INOD) = NOD_COUNT(MAT_INOD) + 1
            END DO
        END DO

        DO MAT_INOD=1,MAT_NONODS
            LES_UDIFFUSION(:,:,:,MAT_INOD) = LES_UDIFFUSION(:,:,:,MAT_INOD)/REAL( NOD_COUNT(MAT_INOD) )
            LES_UDIFFUSION_VOL(:,MAT_INOD) = LES_UDIFFUSION_VOL(:,MAT_INOD)/REAL( NOD_COUNT(MAT_INOD) )
        END DO

        contains

            SUBROUTINE VISCOCITY_TENSOR_LES_CALC_U( LES_U_UDIFFUSION, LES_U_UDIFFUSION_VOL, Q_SCHEME_ABS_CONT_VOL, &
                DUX_ELE_ALL, NDIM,NPHASE, U_NLOC,X_NLOC,TOTELE, X_NONODS, &
                X_ALL, X_NDGLN, LES_DISOPT, CS,  U_NDGLN, U_NONODS, U_ALL)
                ! This subroutine calculates a tensor of viscocity LES_UDIFFUSION, LES_U_UDIFFUSION_VOL
                IMPLICIT NONE
                INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, TOTELE, X_NONODS, LES_DISOPT, U_NONODS
                REAL, intent( in ) :: CS
                ! LES_DISOPT is LES option e.g. =1 Anisotropic element length scale
                INTEGER, DIMENSION( X_NLOC * TOTELE  ), intent( in ) :: X_NDGLN
                INTEGER, DIMENSION( U_NLOC * TOTELE  ), intent( in ) :: U_NDGLN
                REAL, DIMENSION( NDIM, X_NONODS  ), intent( in ) :: X_ALL
                REAL, DIMENSION( NDIM, NDIM, NPHASE, U_NLOC, TOTELE  ), intent( inout ) :: LES_U_UDIFFUSION
                REAL, DIMENSION( NPHASE, U_NLOC, TOTELE  ), intent( inout ) :: LES_U_UDIFFUSION_VOL
                REAL, DIMENSION( NPHASE, U_NLOC, TOTELE  ), intent( inout ) :: Q_SCHEME_ABS_CONT_VOL
                REAL, DIMENSION( NDIM, NDIM, NPHASE, U_NLOC, TOTELE  ), intent( in ) :: DUX_ELE_ALL
                REAL, DIMENSION( NDIM, NPHASE, U_NONODS  ), intent( in ) :: U_ALL
                ! Local variables...
                LOGICAL, PARAMETER :: ONE_OVER_H2=.FALSE.
                !     SET to metric which has 1/h^2 in it
                ! CQ controls the amount of original q-scheme viscocity ( ~ 1.0 )
                ! C_L The liquid constant for q-like scheme, similary for C_L (stardard value 0.05)
                REAL, PARAMETER :: CQ=1.0, C_L=0.05
                !            REAL, PARAMETER :: CS=0.1
                REAL :: LOC_X_ALL(NDIM, X_NLOC), TENSXX_ALL(NDIM, NDIM, NPHASE), RSUM, FOURCS, CS2, VIS, DIVU, H2
                REAL :: MEAN_UDER_U(NDIM, NDIM, NPHASE)
                INTEGER :: ELE, X_ILOC, U_ILOC, IPHASE, X_NODI, IDIM, JDIM, U_INOD
                REAL :: RN,WEIGHT

                CS2=CS**2
                FOURCS=CS2

                !            IF(LES_DISOPT.ge.8) THEN
                !               ALLOCATE(Q_SCHEME_ABS_CONT_VOL( NPHASE, U_NLOC, TOTELE  ))
                !            ENDIF

                Q_SCHEME_ABS_CONT_VOL=0.0

                DO ELE=1,TOTELE

                    DO X_ILOC=1,X_NLOC
                        X_NODI = X_NDGLN((ELE-1)*X_NLOC+X_ILOC)
                        LOC_X_ALL(:,X_ILOC) = X_ALL(:,X_NODI)
                    END DO

                    IF(LES_DISOPT.ge.5) THEN
                        DO IDIM=1,NDIM
                            DO JDIM=1,NDIM
                                MEAN_UDER_U(IDIM, JDIM, :) =  SUM( DUX_ELE_ALL(IDIM,JDIM,:,1:U_NLOC,ELE) )/REAL(U_NLOC)
                            END DO
                            ! Normalise to be of size unity...
                            DO IPHASE=1,NPHASE
                                MEAN_UDER_U(IDIM, :,IPHASE) = MEAN_UDER_U(IDIM, :,IPHASE) / MAX(SQRT(SUM( MEAN_UDER_U(IDIM, :,IPHASE)**2 )),  1.E-10)
                            END DO
                        END DO
                    ENDIF

                    CALL ONEELETENS_ALL( LOC_X_ALL, LES_DISOPT, ONE_OVER_H2, TENSXX_ALL, X_NLOC, NDIM, MEAN_UDER_U, NPHASE )

                    DO U_ILOC=1,U_NLOC
                        U_INOD = U_NDGLN((ELE-1)*U_NLOC+U_ILOC)
                        DO IPHASE=1,NPHASE

                            RSUM=0.0
                            DO IDIM=1,NDIM
                                DO JDIM=1,NDIM
                                    RSUM=RSUM + (0.5*( DUX_ELE_ALL(IDIM,JDIM,IPHASE,U_ILOC,ELE) + DUX_ELE_ALL(JDIM,IDIM,IPHASE,U_ILOC,ELE) ))**2
                                END DO
                            END DO
                            RSUM=SQRT(RSUM)
                            VIS=RSUM

                            ! for original q-scheme:
                            DIVU=0.0
                            H2=0.0
                            IF(LES_DISOPT.GE.6) THEN
                                DO IDIM=1,NDIM
                                    DIVU=DIVU + DUX_ELE_ALL(IDIM,IDIM,IPHASE,U_ILOC,ELE)
                                END DO

                                RN=0.0
                                DO IDIM=1,NDIM
                                    WEIGHT=MAX(1.E-10, ABS(U_ALL(IDIM,IPHASE,U_INOD)) )
                                    H2=H2  + TENSXX_ALL(IDIM,IDIM,IPHASE)*WEIGHT
                                    RN=RN+WEIGHT
                                END DO
                                H2=H2/RN
                            ENDIF ! ENDOF IF(LES_DISOPT.GE.6) THEN

                            ! THEN FIND TURBULENT 'VISCOSITIES'
                            !tENSXX_ALL(:,:)=6./40.

                            ! Put a bit in here which multiplies E by FOURCS*VIS
                            LES_U_UDIFFUSION(:,:,IPHASE,U_ILOC,ELE)= FOURCS*VIS*TENSXX_ALL(:,:,IPHASE)

                            ! This is the original q-scheme...
                            IF(LES_DISOPT==6) THEN
                                LES_U_UDIFFUSION_VOL(IPHASE,U_ILOC,ELE)= -CQ*H2*MIN(0.0, DIVU)
                            ELSE IF(LES_DISOPT==7) THEN
                                LES_U_UDIFFUSION_VOL(IPHASE,U_ILOC,ELE)=  CQ*H2*ABS(DIVU)
                            ELSE IF(LES_DISOPT==8) THEN
                                LES_U_UDIFFUSION_VOL(IPHASE,U_ILOC,ELE)= -CQ*H2*MIN(0.0, DIVU)
                                Q_SCHEME_ABS_CONT_VOL(IPHASE,U_ILOC,ELE)= -C_L * SQRT(H2) * MIN(0.0, SIGN(1.0, DIVU) )
                            ELSE IF(LES_DISOPT==9) THEN
                                LES_U_UDIFFUSION_VOL(IPHASE,U_ILOC,ELE)=  CQ*H2*ABS(DIVU)
                                Q_SCHEME_ABS_CONT_VOL(IPHASE,U_ILOC,ELE)= -C_L * SQRT(H2) * MIN(0.0, SIGN(1.0, DIVU) )
                            ELSE IF(LES_DISOPT==10) THEN
                                LES_U_UDIFFUSION_VOL(IPHASE,U_ILOC,ELE)= -CQ*H2*MIN(0.0, DIVU)
                                Q_SCHEME_ABS_CONT_VOL(IPHASE,U_ILOC,ELE)= C_L * SQRT(H2)
                            ELSE IF(LES_DISOPT==11) THEN
                                LES_U_UDIFFUSION_VOL(IPHASE,U_ILOC,ELE)=  CQ*H2*ABS(DIVU)
                                Q_SCHEME_ABS_CONT_VOL(IPHASE,U_ILOC,ELE)= C_L * SQRT(H2)
                            ELSE
                                LES_U_UDIFFUSION_VOL(IPHASE,U_ILOC,ELE)=  0.0
                            ENDIF

                        END DO ! DO IPHASE=1,NPHASE
                    END DO ! DO U_ILOC=1,U_NLOC

                END DO
                RETURN
            END SUBROUTINE VISCOCITY_TENSOR_LES_CALC_U


            SUBROUTINE ONEELETENS_ALL( LOC_X_ALL, LES_DISOPT, ONE_OVER_H2, TENSXX_ALL, X_NLOC, NDIM, MEAN_UDER_U, NPHASE )
                !     This sub calculates the ELEMENT-WISE TENSOR TENS
                !     REPRESENTS THE SIZE AND SHAPE OF THE SURROUNDING ELEMENTS.
                !     LES_DISOPT=LES option.
                IMPLICIT NONE
                INTEGER, intent( in ) ::  X_NLOC, NDIM, NPHASE
                LOGICAL, intent( in ) ::  ONE_OVER_H2
                INTEGER, intent( in ) ::  LES_DISOPT

                REAL, intent( inout ) ::  TENSXX_ALL(NDIM,NDIM, NPHASE)
                REAL, intent( in ) ::  MEAN_UDER_U(NDIM,NDIM, NPHASE)
                REAL, intent( in ) ::  LOC_X_ALL(NDIM,X_NLOC)

                !     HX,HY-characteristic length scales in x,y directions.
                !     Local variables...
                ! IF ONE_OVER_H2=.TRUE. then SET to metric which has 1/h^2 in it
                REAL RN
                REAL AA(NDIM,NDIM),V(NDIM,NDIM),D(NDIM),A(NDIM,NDIM)

                REAL UDL_ALL(NDIM, X_NLOC*X_NLOC)
                REAL GAMMA(X_NLOC*X_NLOC)

                INTEGER L1,L2,ID,NID,IDIM,JDIM,KDIM,I,IPHASE

                REAL HOVERQ
                REAL RWIND, D_SCALAR
                REAL RFACT

                RWIND =1./REAL(6)
                NID=X_NLOC*X_NLOC

                TENSXX_ALL=0.0

                id=0
                do L1=1,X_NLOC
                    do L2=1,X_NLOC
                        id=id+1
                        UDL_ALL(:,ID)=LOC_X_ALL(:,L1)-LOC_X_ALL(:,L2)
                        RN=SQRT( SUM(UDL_ALL(:,ID)**2) )
                        GAMMA(ID)=RN
                    END DO
                END DO

                IF(LES_DISOPT==1) THEN ! Take the anisotropic length scales
                    !     This subroutine forms a contabution to the Right Hand Side
                    !     of Poissons pressure equation, as well as  F1 & F2.


                    !     C The first is the old filter term, the second the new one MDP getting
                    !     c different results and stabiltiy for tidal applications ????
                    !     **********calculate normalised velocitys across element...
                    ID=0
                    do L1=1,X_NLOC
                        do L2=1,X_NLOC
                            ID=ID+1
                            if(l1.eq.l2) then
                                UDL_ALL(:,ID)=0.0
                                GAMMA(ID)=0.0
                            else
                                UDL_ALL(:,ID)=LOC_X_ALL(:,L1)-LOC_X_ALL(:,L2)

                                !     Normalise
                                RN=SQRT( SUM(UDL_ALL(:,ID)**2) )
                                UDL_ALL(:,ID)=UDL_ALL(:,ID)/RN
                                !     HX,HY are the characteristic length scales in x,y directions.
                                HOVERQ=RN
                                GAMMA(ID)=RWIND*HOVERQ
                            endif
                        END DO
                    END DO
                    !     **********calculate normalised velocitys across element...


                    do  ID=1,NID

                        RFACT=GAMMA(ID)/REAL(X_NLOC)

                        DO IDIM=1,NDIM
                            DO JDIM=1,NDIM
                                TENSXX_ALL(IDIM,JDIM,1)=TENSXX_ALL(IDIM,JDIM,1) + RFACT*UDL_ALL(IDIM,ID)*UDL_ALL(JDIM,ID)
                            END DO
                        END DO

                       !     USE THE COMPONENT OF DIFLIN THE X,Y & Z-DIRECTIONS
                       !     RESPECTIVELY FOR C1T,C2T,C3T.
                    end do

                    !     nb we want 1/L^2 - at the moment we have L on the diagonal.
                    !     Make sure the eigen-values are positive...
                    AA(:,:)=TENSXX_ALL(:,:,1)

                    CALL JACDIA(AA,V,D,NDIM,A)


                    IF(ONE_OVER_H2) THEN
                        !     SET to metric which has 1/h^2 in it...
                        D(:)=1./MAX(1.E-16,D(:)**2)
                    ELSE
                        !     set to inverse of metric which is a multiple of the tensor
                        D(:)=MAX(1.E-16,D(:)**2)
                    ENDIF

                    TENSXX_ALL=0.0
                    DO IDIM=1,NDIM
                        DO JDIM=1,NDIM

                            DO KDIM=1,NDIM
                                ! TENSOR=V^T D V
                                TENSXX_ALL(IDIM,JDIM,:)=TENSXX_ALL(IDIM,JDIM,:) + V(KDIM,IDIM) * D(KDIM) * V(KDIM,JDIM)
                            END DO

                        END DO
                    END DO

                ELSE IF(LES_DISOPT==2) THEN ! Take the average length scale h
                    D_SCALAR=SUM(GAMMA(:))/REAL(X_NLOC*(X_NLOC-1))
                    DO I=1,NDIM
                        TENSXX_ALL(I,I,:)=D_SCALAR**2
                    END DO
                ELSE IF(LES_DISOPT==3) THEN ! Take the min length scale h
                    D_SCALAR=MINVAL(GAMMA(:))
                    DO I=1,NDIM
                        TENSXX_ALL(I,I,:)=D_SCALAR**2
                    END DO
                ELSE IF(LES_DISOPT==4) THEN ! Take the max length scale h
                    D_SCALAR=MAXVAL(GAMMA(:))
                    DO I=1,NDIM
                        TENSXX_ALL(I,I,:)=D_SCALAR**2
                    END DO
                ELSE  IF(LES_DISOPT.ge.5) THEN ! us option good for stress form
                    ! Take the directions that are a max across each element.
                    DO IPHASE=1,NPHASE
                        DO IDIM=1,NDIM
                            RN=0.0
                            do L1=1,X_NLOC
                                do L2=1,X_NLOC
                                    RN = MAX(RN,   SUM( (LOC_X_ALL(:,L1)-LOC_X_ALL(:,L2))*MEAN_UDER_U(IDIM,:,IPHASE))**2 )
                                END DO
                            END DO
                            TENSXX_ALL(IDIM,IDIM,IPHASE)=SQRT( RN )
                        END DO
                    END DO

                ELSE
                    !            ERROR("NOT A VALID OPTION FOR LES ASSEMBLED EQNS")
                    STOP 9331
                ENDIF

                RETURN
            END SUBROUTINE ONEELETENS_ALL

            !
            !!sprint_to_do!!!!MOVE TO FORTRAN 90
            SUBROUTINE JACDIA(AA,V,D,N, &
                ! Working arrays...
                A)
                ! This sub performs Jacobi rotations of a symmetric matrix in order to
                ! find the eigen-vectors V and the eigen values A so
                ! that AA=V^T D V & D is diagonal.
                ! It uses the algorithm of Matrix Computations 2nd edition, p196.
                IMPLICIT NONE
                REAL TOLER,CONVEG
                PARAMETER(TOLER=1.E-14,CONVEG=1.E-7)
                INTEGER N
                REAL AA(N,N),V(N,N),D(N), A(N,N)
                ! Local variables...
                REAL R,ABSA,MAXA,COSAL2,COSALF,SINAL2,SINALF,MAXEIG
                INTEGER ITS,NITS,Q,P,QQ,PP
                !
                NITS=9*(N*N-N)/2
                !
                !          CALL RCLEAR(V,N*N)
                !          CALL TOCOPY(A,AA,N*N)
                do P=1,N
                    do Q=1,N
                        V(P,Q)=0.
                        A(P,Q)=AA(P,Q)
                    !             ewrite(2,*) 'P,Q,AA:',P,Q,AA(P,Q)
                    END DO
                END DO
                !
                !
                !
                !     Check first whether matrix is diagonal
                IF(Q.EQ.0) THEN

                    do PP=1,N
                        D(PP) = A(PP,PP)
                        do QQ=1,N
                            IF(PP.EQ.QQ) THEN
                                V(PP,QQ) = 1.0
                            ELSE
                                V(PP,QQ) = 0.0
                            END IF
                        END DO
                    END DO
                    RETURN
                END IF
                !
                !
                !
                MAXEIG=0.
                do P=1,N
                    V(P,P)=1.0
                    MAXEIG=MAX(MAXEIG,ABS(A(P,P)))
                END DO
                IF(MAXEIG.LT.TOLER) THEN
                    D(1:N) = 0.0
                    GOTO 2000
                ENDIF
                !           ewrite(2,*) 'maxeig=',maxeig
                !
                do  ITS=1,NITS! Was loop 10
                    ! Find maximum on upper diagonal of matrix.
                    ! QQ is the coln; PP is the row.
                    Q=0
                    P=0
                    MAXA=0.
                    do PP=1,N-1
                        do QQ=PP+1,N
                            ABSA=ABS(A(PP,QQ))
                            IF(ABSA.GT.MAXA) THEN
                                MAXA=ABSA
                                Q=QQ
                                P=PP
                            ENDIF
                        END DO
                    END DO

                    !            IF(PRISCR) ewrite(2,*) 'MAXA,MAXEIG,its=',MAXA,MAXEIG,its
                    IF(MAXA/MAXEIG.LT.CONVEG) GOTO 2000
                    ! Rotate with (Q,P) postions.
                    R=MAX(TOLER,SQRT( (A(P,P)-A(Q,Q))**2 + 4.*A(P,Q)**2 ) )
                    IF(A(P,P).GT.A(Q,Q)) THEN
                        COSAL2=0.5+0.5*(A(P,P)-A(Q,Q))/R
                        COSALF=SQRT(COSAL2)
                        IF(ABS(COSALF).LT.TOLER) COSALF=TOLER
                        SINALF=A(Q,P)/(R*COSALF)
                    ELSE
                        SINAL2=0.5-0.5*(A(P,P)-A(Q,Q))/R
                        SINALF=SQRT(SINAL2)
                        IF(ABS(SINALF).LT.TOLER) SINALF=TOLER
                        COSALF=A(Q,P)/(R*SINALF)
                    ENDIF
                    ! Pre and Post multiply of A=R^T A R  by rotation matrix.
                    CALL JACPRE(-SINALF,COSALF,P,Q,A,N)
                    CALL JACPOS( SINALF,COSALF,P,Q,A,N)
                    ! Accumulate rotations V=R^T V
                    CALL JACPRE(-SINALF,COSALF,P,Q,V,N)
                end do ! Was loop 10
            !
2000        CONTINUE
            ! Put e-values in a vector...
            do Q=1,N
                D(Q)=A(Q,Q)
            END DO
            !
            RETURN
        END SUBROUTINE JACDIA
        !
        !
        !
        !!sprint_to_do!!!OPTIMIZE, VECTORIZE
        SUBROUTINE JACPRE(SINALF,COSALF,P,Q,A,N)
            ! This sub performs matrix-matrix multiplication A=R*A.
            ! PRE-MULTIPLY matrix A by transpose of Rotation matrix
            ! is realised by passing -SINALF down into SINALF.
            IMPLICIT NONE
            INTEGER N
            REAL SINALF,COSALF,A(N,N)
            INTEGER P,Q
            ! Local variables...
            INTEGER I
            REAL P1I
            !
            ! Premultiply by rotation matrix...
            do I=1,N
                ! Row P 1st...
                P1I   =COSALF*A(P,I)-SINALF*A(Q,I)
                ! Row 2nd put strait in A...
                A(Q,I)=SINALF*A(P,I)+COSALF*A(Q,I)
                A(P,I)=P1I
            END DO
            RETURN
        END SUBROUTINE JACPRE
        !
        !
        !
        !
        SUBROUTINE JACPOS(SINALF,COSALF,P,Q,A,N)
            ! This sub performs matrix-matrix multiplication A=A*R.
            ! POST-MULTIPLY matrix A by transpose of Rotation matrix
            ! is realised by passing -SINALF down into SINALF.
            IMPLICIT NONE
            INTEGER N
            REAL SINALF,COSALF,A(N,N)
            INTEGER P,Q
            ! Local variables...
            INTEGER I
            REAL IP1
            !
            ! Post multiply by rotation matrix...
            do I=1,N
                ! Column P 1st...
                IP1   = COSALF*A(I,P)+SINALF*A(I,Q)
                ! column 2nd put strait in A...
                A(I,Q)=-SINALF*A(I,P)+COSALF*A(I,Q)
                A(I,P)=IP1
            END DO
            !
            RETURN
        END SUBROUTINE JACPOS

    END SUBROUTINE VISCOCITY_TENSOR_LES_CALC



 SUBROUTINE COMB_VEL_MATRIX_DIAG_DIST(DIAG_BIGM_CON, BIGM_CON, &
     DGM_PETSC, &
     FINELE, COLELE,  NDIM_VEL, NPHASE, U_NLOC, TOTELE, velocity, pressure)  ! Element connectivity.
     ! This subroutine combines the distributed and block diagonal for an element
     ! into the matrix DGM_PHA.
     IMPLICIT NONE
     INTEGER, intent( in ) :: NDIM_VEL, NPHASE, U_NLOC, TOTELE
     !
     REAL, DIMENSION( :,:,:, :,:,:, : ), intent( in ) :: DIAG_BIGM_CON
     REAL, DIMENSION( :,:,:, :,:,:, : ), intent( in ) :: BIGM_CON
     type( petsc_csr_matrix ), intent( inout ) :: DGM_PETSC
     INTEGER, DIMENSION(: ), intent( in ) :: FINELE
     INTEGER, DIMENSION( : ), intent( in ) :: COLELE
     type( tensor_field ) :: velocity
     type( tensor_field ) :: pressure

     INTEGER :: ELE,ELE_ROW_START,ELE_ROW_START_NEXT,ELE_IN_ROW
     INTEGER :: U_ILOC,U_JLOC, IPHASE,JPHASE, IDIM,JDIM, I,J, GLOBI, GLOBJ
     INTEGER :: COUNT_ELE,JCOLELE
     real, dimension(:,:,:, :,:,:), allocatable :: LOC_DGM_PHA

     integer, dimension(:), pointer :: neighbours
     integer :: nb
     logical :: skip


     ALLOCATE(LOC_DGM_PHA(NDIM_VEL,NDIM_VEL,NPHASE,NPHASE,U_NLOC,U_NLOC))

     Loop_Elements20: DO ELE = 1, TOTELE

         if (IsParallel()) then
             if (.not. assemble_ele(pressure,ele)) then
                 skip=.true.
                 neighbours=>ele_neigh(pressure,ele)
                 do nb=1,size(neighbours)
                     if (neighbours(nb)<=0) cycle
                     if (assemble_ele(pressure,neighbours(nb))) then
                         skip=.false.
                         exit
                     end if
                 end do
                 if (skip) cycle
             end if
         end if

         ELE_ROW_START=FINELE(ELE)
         ELE_ROW_START_NEXT=FINELE(ELE+1)
         ELE_IN_ROW = ELE_ROW_START_NEXT - ELE_ROW_START

         ! Block diagonal and off diagonal terms...
         Between_Elements_And_Boundary20: DO COUNT_ELE=ELE_ROW_START, ELE_ROW_START_NEXT-1

             JCOLELE=COLELE(COUNT_ELE)

             IF(JCOLELE==ELE) THEN
                 ! Block diagonal terms (Assume full coupling between the phases and dimensions)...
                 LOC_DGM_PHA(:,:,:, :,:,:) = DIAG_BIGM_CON(:,:,:, :,:,:, ELE) + BIGM_CON(:,:,:, :,:,:, COUNT_ELE)
             ELSE
                 LOC_DGM_PHA(:,:,:, :,:,:) = BIGM_CON(:,:,:, :,:,:, COUNT_ELE)
             ENDIF

             DO U_JLOC=1,U_NLOC
                 DO U_ILOC=1,U_NLOC
                     DO JPHASE=1,NPHASE
                         DO IPHASE=1,NPHASE
                             DO JDIM=1,NDIM_VEL
                                 DO IDIM=1,NDIM_VEL
                                     ! New for rapid code ordering of variables...
                                     I=IDIM + (IPHASE-1)*NDIM_VEL
                                     J=JDIM + (JPHASE-1)*NDIM_VEL
                                     GLOBI=(ELE-1)*U_NLOC + U_ILOC
                                     GLOBJ=(JCOLELE-1)*U_NLOC + U_JLOC
                                     if (.not. node_owned(velocity,globi)) cycle
                                     call addto(dgm_petsc, I , J , globi , globj , &
                                         LOC_DGM_PHA(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC))
                                 END DO
                             END DO
                         END DO
                     END DO
                 END DO
             END DO


         END DO Between_Elements_And_Boundary20

     END DO Loop_Elements20


     RETURN
 END SUBROUTINE COMB_VEL_MATRIX_DIAG_DIST










 SUBROUTINE USE_POSINMAT_C_STORE(COUNT, U_INOD, P_JNOD,  &
     U_NONODS, FINDC, COLC, NCOLC, &
     IDO_STORE_AC_SPAR_PT,STORED_AC_SPAR_PT, POSINMAT_C_STORE,ELE,U_ILOC,P_JLOC, &
     TOTELE,U_NLOC,P_NLOC)
     INTEGER, intent( inout ) :: COUNT
     INTEGER, intent( in ) :: U_INOD, P_JNOD, U_NONODS,  NCOLC
     INTEGER, intent( in ) :: ELE,U_ILOC,P_JLOC,  TOTELE,U_NLOC,P_NLOC
     INTEGER, intent( in ) :: IDO_STORE_AC_SPAR_PT
     LOGICAL, intent( in ) :: STORED_AC_SPAR_PT
     INTEGER, DIMENSION( U_NLOC,P_NLOC, TOTELE*IDO_STORE_AC_SPAR_PT), intent( inout ) :: POSINMAT_C_STORE
     INTEGER, DIMENSION( U_NONODS + 1), intent( in ) :: FINDC
     INTEGER, DIMENSION( NCOLC), intent( in ) :: COLC

     ! Find COUNT - position in matrix : FINMCY, COLMCY
     IF(STORED_AC_SPAR_PT) THEN
         COUNT=POSINMAT_C_STORE(U_ILOC,P_JLOC,ELE)
     ELSE
         CALL POSINMAT( COUNT, U_INOD, P_JNOD,  &
             FINDC, COLC )
         IF(IDO_STORE_AC_SPAR_PT.NE.0) POSINMAT_C_STORE(U_ILOC,P_JLOC,ELE) = COUNT
     ENDIF
     RETURN
 END SUBROUTINE USE_POSINMAT_C_STORE




 SUBROUTINE USE_POSINMAT_C_STORE_SUF_DG(COUNT, U_INOD, P_JNOD,  &
     U_NONODS, FINDC, COLC, NCOLC, &
     IDO_STORE_AC_SPAR_PT,STORED_AC_SPAR_PT, POSINMAT_C_STORE_SUF_DG, ELE,IFACE,U_SILOC,P_SJLOC,  &
     TOTELE,NFACE,U_SNLOC,P_SNLOC)
     INTEGER, intent( inout ) :: COUNT
     INTEGER, intent( in ) :: U_INOD, P_JNOD, U_NONODS,  NCOLC
     INTEGER, intent( in ) :: ELE,IFACE,U_SILOC,P_SJLOC,  TOTELE,NFACE,U_SNLOC,P_SNLOC
     INTEGER, intent( in ) :: IDO_STORE_AC_SPAR_PT
     LOGICAL, intent( in ) :: STORED_AC_SPAR_PT
     INTEGER, DIMENSION( U_SNLOC,P_SNLOC,NFACE,TOTELE*IDO_STORE_AC_SPAR_PT ), intent( inout ) :: POSINMAT_C_STORE_SUF_DG
     INTEGER, DIMENSION( U_NONODS + 1), intent( in ) :: FINDC
     INTEGER, DIMENSION( NCOLC), intent( in ) :: COLC
     ! Find COUNT2 - position in matrix : FINMCY, COLMCY
     IF(STORED_AC_SPAR_PT) THEN
         COUNT=POSINMAT_C_STORE_SUF_DG(U_SILOC,P_SJLOC,IFACE,ELE)
     ELSE
         CALL POSINMAT( COUNT, U_INOD, P_JNOD,  &
             FINDC, COLC )
         IF(IDO_STORE_AC_SPAR_PT.NE.0) POSINMAT_C_STORE_SUF_DG(U_SILOC,P_SJLOC,IFACE,ELE)=COUNT
     ENDIF
     RETURN
 END SUBROUTINE USE_POSINMAT_C_STORE_SUF_DG







 subroutine linearise_field( field_in, field_out )
     implicit none
     type( tensor_field ), intent( in ) :: field_in
     real, dimension( :, : ), intent( inout ) :: field_out

     integer, dimension( : ), pointer :: ndglno, cv_nods
     integer :: n, totele, cv_nloc, ncomp, nphase, cv_nonods, ele, cv_iloc, cv_nod
     real, dimension( :, :, : ), allocatable :: field_tmp, field_cv_nod

     ! This sub will linearise a p2 field

     field_out = 0.0

     n = field_in%mesh%shape%degree

     if ( n==2 ) then

         ndglno => get_ndglno( field_in%mesh )

         totele = field_in%mesh%elements
         cv_nloc = field_in%mesh%shape%loc
         cv_nonods =  field_in%mesh%nodes

         ncomp = size( field_in%val, 1 )
         nphase = size( field_in%val, 2 )

         allocate( field_tmp( ncomp, nphase, cv_nonods ) ) ; field_tmp = field_in % val
         allocate( field_cv_nod( ncomp, nphase, cv_nloc ) ) ; field_cv_nod = 0.0

         do ele = 1, totele

             cv_nods => ndglno( ( ele - 1 ) * cv_nloc + 1 : ele * cv_nloc )
             field_cv_nod =  field_tmp( :, :, cv_nods )

             field_cv_nod( :, :, 2 ) = 0.5 * ( field_cv_nod( :, :, 1 ) + field_cv_nod( :, :, 3 ) )
             field_cv_nod( :, :, 4 ) = 0.5 * ( field_cv_nod( :, :, 1 ) + field_cv_nod( :, :, 6 ) )
             field_cv_nod( :, :, 5 ) = 0.5 * ( field_cv_nod( :, :, 3 ) + field_cv_nod( :, :, 6 ) )

             if ( cv_nloc == 10 ) then
                 field_cv_nod( :, :, 7 ) = 0.5 * ( field_cv_nod( :, :, 1 ) + field_cv_nod( :, :, 10 ) )
                 field_cv_nod( :, :, 8 ) = 0.5 * ( field_cv_nod( :, :, 3 ) + field_cv_nod( :, :, 10 ) )
                 field_cv_nod( :, :, 9 ) = 0.5 * ( field_cv_nod( :, :, 6 ) + field_cv_nod( :, :, 10 ) )
             end if

             do cv_iloc = 1, cv_nloc
                 cv_nod = ndglno( ( ele - 1 ) * cv_nloc + cv_iloc )
                 field_out( :, cv_nod ) = field_cv_nod( 1, :, cv_iloc )
             end do

         end do

         deallocate( field_tmp, field_cv_nod )

     end if

     return
 end subroutine linearise_field


 subroutine Introduce_Cap_press_term(packed_state, Mdims, FE_funs, Devfuns, &
     X_ALL, LOC_U_RHS, ele, cv_ndgln, x_ndgln,&
     ele2, iface, sdetwe, SNORMXN_ALL, U_SLOC2LOC, CV_SLOC2LOC, MAT_OTHER_LOC)
     !This subroutine introduces the capillary pressure term in the RHS
     Implicit none
     type( state_type ), intent( inout ) :: packed_state
     type(multi_dimensions), intent(in) :: Mdims
     type(multi_shape_funs), intent(in) :: FE_funs
     integer, intent(in) :: ele, iface, ele2
     integer, dimension(:), intent(in) :: cv_ndgln, x_ndgln
     REAL, DIMENSION ( :, :, : ), intent(inout) :: LOC_U_RHS
     real, dimension(:,:) :: X_ALL
     real, dimension(:,:), intent(in) :: SNORMXN_ALL
     integer, dimension(:), intent(in) :: U_SLOC2LOC, CV_SLOC2LOC, MAT_OTHER_LOC
     real, dimension(:), intent(in) :: sdetwe
     type(multi_dev_shape_funs), intent(inout) :: Devfuns
     !Local parameters
     !!Use a finite element projection of the CapPressure, it can only be false for PnDGPn(DG) elements
     logical, parameter :: Cap_to_FEM = .true.
     !Use integration by parts to introduce the CapPressure, otherwise it uses the integration by parts twice approach
     logical, parameter :: Int_by_part_CapPress = .false.
     !Local variables
     integer :: iphase, cv_inod, u_siloc, cv_jloc,&
         CV_SJLOC, u_iloc, cv_Xnod
     real, dimension(:,:), pointer :: CapPressure
     real, pointer, dimension(:, :) :: CV_Bound_Shape_Func
     real, pointer, dimension(:, :) :: CV_Shape_Func
     real, dimension(Mdims%NDIM) :: NMX_ALL


     call get_var_from_packed_state(packed_state, CapPressure = CapPressure)

     !Retrieve derivatives of the shape functions
     call DETNLXR_PLUS_U(ELE, X_ALL, X_NDGLN, FE_funs%cvweight, &
        FE_funs%cvfen, FE_funs%cvfenlx_all, FE_funs%ufenlx_all, Devfuns)

     !Project to FEM
     if (CAP_to_FEM) then
         !Point my pointers to the FEM shape functions
         CV_Bound_Shape_Func => FE_funs%sbcvfen
         CV_Shape_Func => FE_funs%cvfen
     else
         !Point my shape functions to the Control volume ones
         CV_Bound_Shape_Func => FE_funs%sbcvn
         CV_Shape_Func => FE_funs%cvn
     end if
     !Integration by parts
     if (Int_by_part_CapPress .or. .not. CAP_to_FEM) then
         if (iface == 1) then!The volumetric term is added just one time
             !Firstly we add the volumetric integral
             DO U_ILOC = 1, Mdims%u_nloc
                 DO CV_JLOC = 1, Mdims%cv_nloc
                     ! -Integral(FE_funs%cvn CapPressure Grad FE_funs%ufen dV)
                     CV_INOD = CV_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_JLOC )
                     DO IPHASE = 1, Mdims%nphase
                         LOC_U_RHS( :, IPHASE, U_ILOC ) = LOC_U_RHS( :, IPHASE, U_ILOC ) &
                             !(FE_funs%cvn Grad FE_funs%ufen)
                             + matmul(Devfuns%UFENX_ALL(:,U_ILOC,:),CV_Shape_Func( CV_JLOC, : ) *Devfuns%detwei )&
                             !CapPressure
                             * CapPressure(IPHASE, CV_INOD)
                     END DO
                 end do
             end do
         end if
         !Performing the surface integral, -Integral(FE_funs%cvn CapPressure ᐁFE_funs%ufen dV)
         DO U_SILOC = 1, Mdims%u_snloc
             U_ILOC = U_SLOC2LOC( U_SILOC )
             DO CV_SJLOC = 1, Mdims%cv_snloc
                 CV_JLOC = CV_SLOC2LOC( CV_SJLOC )
                 CV_INOD = CV_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_JLOC )
                 NMX_ALL = matmul(SNORMXN_ALL( :, : ), FE_funs%sbufen( U_SILOC, : ) &
                     * CV_Bound_Shape_Func( CV_SJLOC, : ) * SDETWE( : ))
                 if (ELE2 > 0) then!If neighbour then we get its value to calculate the average
                     cv_Xnod = CV_NDGLN( ( ELE2 - 1 ) * Mdims%cv_nloc + MAT_OTHER_LOC(CV_JLOC) )
                 else!If no neighbour then we use the same value.
                     cv_Xnod = CV_INOD
                 end if
                 do iphase = 1, Mdims%nphase
                     LOC_U_RHS( :, IPHASE, U_ILOC) =  LOC_U_RHS( :, IPHASE, U_ILOC ) &
                         - NMX_ALL(:) * 0.5*(CapPressure(iphase, CV_INOD)+CapPressure(iphase, cv_Xnod))
                 end do
             end do
         end do
     else!Volumetric integration only (requires the CapPressure to be in FEM)
         if (iface ==1) then!The volumetric term is added just one time
             DO U_ILOC = 1, Mdims%u_nloc
                 DO CV_JLOC = 1, Mdims%cv_nloc
                     ! Integral(Grad FE_funs%cvn CapPressure FE_funs%ufen dV)
                     CV_INOD = CV_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_JLOC )
                     DO IPHASE = 1, Mdims%nphase
                         LOC_U_RHS( :, IPHASE, U_ILOC ) = LOC_U_RHS( :, IPHASE, U_ILOC ) &
                             !(Grad FE_funs%cvn FE_funs%ufen)
                             - matmul(Devfuns%CVFENX_ALL(:,CV_JLOC,:),FE_funs%ufen( U_ILOC, : ) *Devfuns%DETWEI )&
                             !CapPressure
                             * CapPressure(IPHASE, CV_INOD)
                     END DO
                 end do
             end do
         end if
         !Get neighbouring nodes
         !Performing the surface integral, Integral(FE_funs%cvn (Average CapPressure) ᐁFE_funs%ufen dV)
         DO U_SILOC = 1, Mdims%u_snloc
             U_ILOC = U_SLOC2LOC( U_SILOC )
             DO CV_SJLOC = 1, Mdims%cv_snloc
                 CV_JLOC = CV_SLOC2LOC( CV_SJLOC )
                 CV_INOD = CV_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_JLOC )
                 NMX_ALL = matmul(SNORMXN_ALL( :, : ), FE_funs%sbufen( U_SILOC, : ) * FE_funs%sbcvfen( CV_SJLOC, : ) * SDETWE( : ))
                 if (ELE2 > 0) then!If neighbour then we get its value to calculate the average
                     cv_Xnod = CV_NDGLN( ( ELE2 - 1 ) * Mdims%cv_nloc + MAT_OTHER_LOC(CV_JLOC) )
                 else!If no neighbour then we use the same value.
                     cv_Xnod = CV_INOD
                 end if
                 do iphase = 1, Mdims%nphase
                     LOC_U_RHS( :, IPHASE, U_ILOC) =  LOC_U_RHS( :, IPHASE, U_ILOC ) &
                         + NMX_ALL(:) * 0.5* (CapPressure(iphase, CV_INOD) - CapPressure(iphase, cv_Xnod))
                 end do
             end do
         end do
     end if

 end subroutine Introduce_Cap_press_term



 subroutine getOverrelaxation_parameter(packed_state, Mdims, ndgln, Overrelaxation, Phase_with_Pc, IDs2CV_ndgln)
     !This subroutine calculates the overrelaxation parameter we introduce in the saturation equation
     !It is the derivative of the capillary pressure for each node.
     !Overrelaxation has to be alocate before calling this subroutine its size is cv_nonods
     implicit none
     type( state_type ), intent(inout) :: packed_state
     type (multi_dimensions), intent(in) :: Mdims
     type(multi_ndgln), intent(in) :: ndgln
     real, dimension(:), intent(inout) :: Overrelaxation
     integer, intent(inout) :: Phase_with_Pc
     integer, dimension(:), intent(in) :: IDs2CV_ndgln
     !Local variables
     real, save :: domain_length = -1
     integer :: iphase, nphase, cv_nodi, cv_nonods, u_inod, cv_iloc, ele, u_iloc
     real :: Pe_aux
     real, dimension(:), pointer ::Pe, Cap_exp
     logical :: Artificial_Pe, Diffusive_cap_only
     real, dimension(:,:,:), pointer :: p
     real, dimension(:,:), pointer :: satura, immobile_fraction, Cap_entry_pressure, Cap_exponent, X_ALL
     type( tensor_field ), pointer :: Velocity

     !Extract variables from packed_state
     call get_var_from_packed_state(packed_state,FEPressure = P,&
         PhaseVolumeFraction = satura, immobile_fraction = immobile_fraction, PressureCoordinate = X_ALL)

     !Initiate local variables
     nphase = size(satura,1)
     cv_nonods = size(satura,2)

     !#######Only apply this method if it has been explicitly invoked through Pe_stab or
     !non-consistent capillary pressure!######
     if (.not.(have_option_for_any_phase("/multiphase_properties/Sat_overRelax", nphase) .or. &
        have_option_for_any_phase("/multiphase_properties/capillary_pressure/Diffusive_cap_only", nphase))) then
         Overrelaxation = 0.0; Phase_with_Pc = -10
         return
     end if
     !######################################################################################

     !Check capillary pressure options
     Phase_with_Pc = -1
     do iphase = Nphase, 1, -1!Going backwards since the wetting phase should be phase 1
         !this way we try to avoid problems if someone introduces 0 capillary pressure in the second phase
         if (have_option( "/material_phase["//int2str(iphase-1)//&
             "]/multiphase_properties/capillary_pressure" ) .or.&
             have_option("/material_phase[["//int2str(iphase-1)//&
             "]/multiphase_properties/Sat_overRelax")) then
             Phase_with_Pc = iphase
         end if
     end do

     Artificial_Pe = .false.
     if (Phase_with_Pc>0) then
         !Get information for capillary pressure to be use
         if (have_option("/material_phase["//int2str(Phase_with_Pc-1)//&
             "]/multiphase_properties/capillary_pressure/type_Brooks_Corey") ) then
             call get_var_from_packed_state(packed_state, Cap_entry_pressure = Cap_entry_pressure,&
                 Cap_exponent = Cap_exponent)
         end if
         !If we want to introduce a stabilization term, this one is imposed over the capillary pressure.
         !Unless we are using the non-consistent form of the capillary pressure
         Diffusive_cap_only = have_option_for_any_phase('/multiphase_properties/capillary_pressure/Diffusive_cap_only', nphase)
         if (have_option("/material_phase["//int2str(Phase_with_Pc-1)//"]/multiphase_properties/Sat_overRelax")&
             .and..not.Diffusive_cap_only) then
             allocate(Pe(CV_NONODS), Cap_exp(CV_NONODS))
             Artificial_Pe = .true.
             call get_option("/material_phase["//int2str(Phase_with_Pc-1)//"]/multiphase_properties/Sat_overRelax", Pe_aux)
             if (Pe_aux<0) then!Automatic set up for Pe
                 !Method based on calculating an entry pressure for a given Peclet number;
                 !Peclet = V * L / Diffusivity; We consider only the entry pressure for the diffusivity
                 !Pe = Vel * L/ Peclet. At present we are using the velocity that includes the sigma. Maybe might be worth it using the Darcy velocity?
                 Velocity => extract_tensor_field( packed_state, "PackedVelocity" )
                 !Since it is an approximation, the domain length is the maximum distance, we only calculate it once
                 if (domain_length < 0) domain_length = abs(maxval(X_ALL)-minval(X_ALL))
                 Pe_aux = abs(Pe_aux)
                 !Obtain an approximation of the capillary number to obtain an entry pressure
                 do ele = 1, Mdims%totele
                     do u_iloc = 1, Mdims%u_nloc
                         u_inod = ndgln%u(( ELE - 1 ) * Mdims%u_nloc +u_iloc )
                         do cv_iloc = 1, Mdims%cv_nloc
                             cv_nodi = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                             Pe(cv_nodi) = (1./Pe_aux) * (sum(abs(Velocity%val(:,Phase_with_Pc,u_inod)))/real(Mdims%ndim) * domain_length)/real(Mdims%u_nloc)
                         end do
                     end do
                 end do
             else
                 Pe = Pe_aux
             end if
!            Cap_exp = 1.!Linear exponent
             Cap_exp = 2.!Quadratic exponent
         end if

         !Calculate the overrrelaxation parameter, the numbering might be different for Pe and real capillary
         !values, hence we calculate it differently
         if (Artificial_Pe) then
             !Calculate the Overrelaxation
             do cv_nodi = 1, size(Overrelaxation)
                 Overrelaxation(CV_NODI) =  Get_DevCapPressure(satura(Phase_with_Pc, CV_NODI),&
                     Pe(CV_NODI), Cap_Exp(CV_NODI), immobile_fraction(:,IDs2CV_ndgln(CV_NODI)), Phase_with_Pc)
             end do
         else
             !Calculate the Overrelaxation
             do cv_nodi = 1, size(Overrelaxation)
                 Overrelaxation(CV_NODI) =  Get_DevCapPressure(satura(Phase_with_Pc, CV_NODI),&
                     Cap_entry_pressure(Phase_with_Pc, IDs2CV_ndgln(CV_NODI)), &
                     Cap_exponent(Phase_with_Pc, IDs2CV_ndgln(CV_NODI)),&
                     immobile_fraction(:,IDs2CV_ndgln(CV_NODI)), Phase_with_Pc)
             end do
         end if

     else
         Overrelaxation = 0.0
     end if

     !Deallocate
     if (Artificial_Pe) then
         deallocate(Pe, Cap_exp)
     end if
     nullify(Pe, Cap_exp)

 end subroutine getOverrelaxation_parameter








subroutine high_order_pressure_solve( Mdims, u_rhs, state, packed_state, nphase, u_absorbin )

      implicit none

      type(multi_dimensions), intent( in ) :: Mdims
      real, dimension( :, :, : ), intent( inout ) :: u_rhs
      type( state_type ), dimension( : ), intent( inout ) :: state
      type( state_type ), intent( inout ) :: packed_state
      integer, intent( in ) :: nphase

      real, dimension( :, :, : ), intent( in ) :: u_absorbin

      ! local variables...
      type ( tensor_field ), pointer :: ufield
      integer :: ndim, ph_ngi, ph_nloc, ph_snloc, &
                u_nloc, u_snloc, stat, &
                totele, x_nonods, ele, x_nloc, &
                ph_ele_type, iloop, u_nonods, cv_nonods, &
                cv_iloc, cv_inod, idim, iphase, u_inod, u_iloc, cv_nloc, &
                ph_iloc, ph_inod, ph_nonods, ph_jloc, ph_jnod, tmp_cv_nloc, other_nloc
      integer, dimension( : ), pointer :: x_ndgln, cv_ndgln, ph_ndgln, u_ndgln, surface_node_list, mat_ndgln
      logical :: quad_over_whole_ele, d1, d3, dcyl
      type( vector_field ), pointer :: x
      type( mesh_type ), pointer :: phmesh

      real, dimension( :, :, : ), allocatable, target :: tmp_cvfenx_all
      real, dimension( :, :, : ), allocatable, target :: other_fenx_all
      real, dimension( : ), allocatable, target :: detwei, ra
      real :: volume

      real, dimension(:,:,:), pointer :: phfenx_all, ufenx_all

      real, dimension( :, :, : ), allocatable :: u_ph_source_vel, u_ph_source_cv
      real, dimension( :, : ), allocatable :: alpha_cv, coef_alpha_cv

      real, dimension( :, :, : ), allocatable :: u_ph_source_ph, dx_ph_gi
      real, dimension( :, : ), allocatable :: alpha_ph, coef_alpha_ph, ph

      real, dimension( :, :, : ), allocatable :: u_s_gi, dx_alpha_gi
      real, dimension( :, : ), allocatable :: coef_alpha_gi, den_gi, inv_den_gi, sigma_gi

      real, dimension( : ), pointer :: tmp_cv_weight
      real, dimension( :, : ), pointer :: tmp_cvfen
      real, dimension( :, :, : ), pointer :: tmp_cvfenlx_all

      real, dimension( :, : ), pointer :: other_fen
      real, dimension( :, :, : ), pointer :: other_fenlx_all

      real :: nxnx, gravity_magnitude, dt

      type( scalar_field ) :: rhs, ph_sol
      type( petsc_csr_matrix ) :: matrix
      type( csr_sparsity ), pointer :: sparsity

      character( len = OPTION_PATH_LEN ) :: path = "/tmp", bc_type

      type( tensor_field ), pointer :: rho, pfield
      type( scalar_field ), pointer :: printf
      type( vector_field ), pointer :: printu, gravity_direction

      logical :: boussinesq, got_free_surf
      integer :: inod, ph_jnod2, ierr, count, count2, i, j, mat_inod
      integer, dimension( : ), pointer :: findph, colph

      type( multi_GI_dimensions ) :: phGIdims
      type( multi_dimensions ) :: phdims
      type( multi_shape_funs ) :: ph_funs


      ewrite(3,*) "inside high_order_pressure_solve"


      boussinesq = have_option( "/material_phase[0]/vector_field::Velocity/prognostic/equation::Boussinesq" )

      call get_option( '/timestepping/timestep', dt )

      printu => extract_vector_field( state( 1 ), "f_x", stat )
      if ( stat == 0 ) call zero( printu  )


      call get_option( '/geometry/dimension', ndim )

      ufield => extract_tensor_field( packed_state, "PackedVelocity" )

      u_nloc = ele_loc( ufield, 1 ) ; u_snloc = face_loc( ufield, 1 )

      u_nonods = node_count( ufield )

      quad_over_whole_ele = .true.

      ! ph elements
      if ( ndim == 2 ) then
         ph_ele_type = 4
         ph_nloc = 6 ; ph_snloc = 3
      else if ( ndim == 3 ) then
         ph_ele_type = 8
         ph_nloc = 10 ; ph_snloc = 6
      else
         stop 567
      end if

      phdims%ndim = Mdims%ndim
      phdims%u_nloc = Mdims%u_nloc ; phdims%u_snloc = Mdims%u_snloc
      phdims%cv_nloc = ph_nloc ; phdims%cv_snloc = ph_snloc

      call retrieve_ngi( phGIdims, phdims, ph_ele_type, quad_over_whole_ele )

      call allocate_multi_shape_funs( ph_funs, phdims, phGIdims )
      call cv_fem_shape_funs( ph_funs, phdims, phGIdims, ph_ele_type, quad_over_whole_ele )

      ph_ngi = phGIdims%cv_ngi
      totele = ele_count( ufield )
      x_ndgln => get_ndglno( extract_mesh( state( 1 ), "PressureMesh_Continuous" ) )
      cv_ndgln => get_ndglno( extract_mesh( state( 1 ), "PressureMesh" ) )
      x_nonods = node_count( extract_mesh( state( 1 ), "PressureMesh_Continuous" ) )
      x => extract_vector_field( packed_state, "PressureCoordinate" )
      u_ndgln => get_ndglno( extract_mesh( state( 1 ), "VelocityMesh" ) )
      x_nloc = ele_loc( x, 1 )
      cv_nloc = x_nloc
      cv_nonods = node_count( extract_mesh( state( 1 ), "PressureMesh" ) )
      phmesh => extract_mesh( state( 1 ), "ph" )
      ph_ndgln => get_ndglno( phmesh )
      ph_nonods = node_count( phmesh )
      d1 = ( ndim == 1 ) ; d3 = ( ndim == 3 ) ; dcyl = .false.

      mat_ndgln => get_ndglno( extract_mesh( state( 1 ), "PressureMesh_Discontinuous" ) )

      if ( cv_nloc == u_nloc ) then

         tmp_cv_nloc = u_nloc
         tmp_cvfen => ph_funs%ufen
         tmp_cvfenlx_all => ph_funs%ufenlx_all
         tmp_cv_weight => ph_funs%cvweight

         other_nloc = ph_nloc
         other_fen => ph_funs%cvfen
         other_fenlx_all => ph_funs%cvfenlx_all

      else if ( cv_nloc == ph_nloc ) then

         tmp_cv_nloc = ph_nloc
         tmp_cvfen => ph_funs%cvfen
         tmp_cvfenlx_all => ph_funs%cvfenlx_all
         tmp_cv_weight => ph_funs%cvweight

         other_nloc = u_nloc
         other_fen => ph_funs%ufen
         other_fenlx_all => ph_funs%ufenlx_all

      else
         stop 7555
      end if

      allocate(tmp_cvfenx_all(ndim, size(tmp_cvfenlx_all,2), ph_ngi))
      allocate(other_fenx_all(ndim, size(other_fenlx_all,2) ,ph_ngi))
      allocate(detwei(ph_ngi))
      allocate(ra(ph_ngi))

      allocate( u_ph_source_vel( ndim, nphase, u_nonods ), &
           &    u_ph_source_cv( ndim, nphase, cv_nonods ), &
           &    alpha_cv( nphase, cv_nonods ), &
           &    coef_alpha_cv( nphase, cv_nonods ) )

      allocate( u_ph_source_ph( ndim, nphase, ph_nonods ), &
           &    alpha_ph( nphase, ph_nonods ), &
           &    ph( nphase, ph_nonods ), &
           &    coef_alpha_ph( nphase, ph_nonods ) )

      allocate( dx_ph_gi( ph_ngi, ndim, nphase ) )

      allocate( u_s_gi( ph_ngi, ndim, nphase ), &
           &    dx_alpha_gi( ph_ngi, ndim, nphase ), &
           &    coef_alpha_gi( ph_ngi, nphase ) )

      allocate( den_gi( ph_ngi, nphase ), inv_den_gi( ph_ngi, nphase ), sigma_gi( ph_ngi, nphase ) )

      ! initialise memory
      u_ph_source_vel = 0.0 ; alpha_cv = 0.0 ; coef_alpha_cv = 0.0
      u_ph_source_ph = 0.0 ; alpha_ph = 0.0 ; coef_alpha_ph = 0.0
      u_ph_source_cv = 0.0

      ! set the gravity term
      if ( have_option( "/physical_parameters/gravity/fem_density_buoyancy" ) ) then
         rho => extract_tensor_field( packed_state, "PackedFEDensity" )
      else
         rho => extract_tensor_field( packed_state, "PackedDensity" )
      end if


      call get_option( "/physical_parameters/gravity/magnitude", gravity_magnitude )
      gravity_direction => extract_vector_field( state( 1 ), "GravityDirection" )

      do idim = 1, ndim
         u_ph_source_cv( idim, 1, : ) = rho % val( 1, 1, : ) * &
                    gravity_magnitude * gravity_direction % val( idim, 1 )
      end do


      sparsity => extract_csr_sparsity( packed_state, "phsparsity" )
      call allocate( matrix, sparsity, [ 1, 1 ], "M", .true. )
      call zero( matrix )

      call allocate( rhs, phmesh, "rhs" )
      call zero ( rhs )
      call allocate( ph_sol, phmesh, "ph_sol" )
      call zero( ph_sol )

      do iloop = 1, 2

         ! iloop=1 form the rhs of the pressure matrix and pressure matrix and solve it.
         ! iloop=2 put the residual into the rhs of the momentum eqn.

         do  ele = 1, totele

            ! calculate detwei,ra,nx,ny,nz for element ele
            call detnlxr_plus_u( ele, x%val(1,:), x%val(2,:), x%val(3,:), &
                 x_ndgln, totele, x_nonods, x_nloc, tmp_cv_nloc, ph_ngi, &
                 tmp_cvfen, tmp_cvfenlx_all(1,:,:), tmp_cvfenlx_all(2,:,:), tmp_cvfenlx_all(3,:,:), &
                 tmp_cv_weight, detwei, ra, volume, d1, d3, dcyl, tmp_cvfenx_all, &
                 other_nloc, other_fenlx_all(1,:,:), other_fenlx_all(2,:,:), other_fenlx_all(3,:,:), &
                 other_fenx_all)

            if ( u_nloc == tmp_cv_nloc ) then
                ufenx_all => tmp_cvfenx_all
            else
                ufenx_all => other_fenx_all
            end if
            if ( ph_nloc == tmp_cv_nloc ) then
                phfenx_all => tmp_cvfenx_all
            else
                phfenx_all => other_fenx_all
            end if

            u_s_gi = 0.0 ; dx_alpha_gi = 0.0 ; coef_alpha_gi = 0.0
            dx_ph_gi = 0.0 ; den_gi = 0.0 ; sigma_gi = 0.0

            do u_iloc = 1, u_nloc
               u_inod = u_ndgln( ( ele - 1 ) * u_nloc + u_iloc )
               do iphase = 1, nphase
                  do idim = 1, ndim
                     u_s_gi( :, idim, iphase ) = u_s_gi( :, idim, iphase ) + &
                          ph_funs%ufen( u_iloc, : ) * u_ph_source_vel( idim, iphase, u_inod )
                  end do
               end do
            end do

            do cv_iloc = 1, cv_nloc
               cv_inod = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )
               mat_inod = mat_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )

               do iphase = 1, nphase
                  do idim = 1, ndim
                     dx_alpha_gi( :, idim, iphase ) = dx_alpha_gi( :, idim, iphase ) + &
                          tmp_cvfenx_all( idim, cv_iloc, : ) * alpha_cv( iphase, cv_inod )
                     u_s_gi( :, idim, iphase ) = u_s_gi( :, idim, iphase ) + &
                          tmp_cvfen( cv_iloc, : ) * u_ph_source_cv( idim, iphase, cv_inod )
                  end do
                  coef_alpha_gi( :, iphase ) = coef_alpha_gi( :, iphase ) + &
                       tmp_cvfen( cv_iloc, : ) * coef_alpha_cv( iphase, cv_inod )

                  if ( boussinesq .or. is_flooding ) then
                     den_gi( :, iphase ) = 1.0
                  else
                     den_gi( :, iphase ) = den_gi( :, iphase ) + &
                          tmp_cvfen( cv_iloc, : ) * rho % val( 1, iphase, cv_inod )
                  end if

                  sigma_gi( :, iphase ) = sigma_gi( :, iphase ) + &
                        tmp_cvfen( cv_iloc, : ) * u_absorbin(  1, iphase, mat_inod )

               end do
            end do

            !inv_den_gi = 1.0 / ( den_gi + dt * sigma_gi )
            inv_den_gi = 1.0 / den_gi


            do ph_iloc = 1, ph_nloc
               ph_inod = ph_ndgln( ( ele - 1 ) * ph_nloc + ph_iloc )
               do iphase = 1, nphase
                  do idim = 1, ndim
                     dx_alpha_gi( :, idim, iphase ) = dx_alpha_gi( :, idim, iphase ) + &
                          phfenx_all( idim, ph_iloc, : ) * alpha_ph( iphase, ph_inod )
                     dx_ph_gi( :, idim, iphase ) = dx_ph_gi( :, idim, iphase ) + &
                          phfenx_all( idim, ph_iloc, : ) * ph( iphase, ph_inod )
                  end do
                  coef_alpha_gi( :, iphase ) = coef_alpha_gi( :, iphase ) + &
                       ph_funs%cvfen( ph_iloc, : ) * coef_alpha_ph( iphase, ph_inod )
               end do
            end do


            if ( iloop == 1 ) then

               ! form the hydrostatic pressure eqn...
               do ph_iloc = 1, ph_nloc
                  ph_inod = ph_ndgln( ( ele - 1 ) * ph_nloc + ph_iloc )
                  do ph_jloc = 1, ph_nloc
                     ph_jnod = ph_ndgln( ( ele - 1 ) * ph_nloc + ph_jloc )
                     nxnx = 0.0
                     do idim = 1, ndim
                        nxnx = nxnx + sum( phfenx_all( idim, ph_iloc, : ) * &
                             phfenx_all( idim, ph_jloc, : ) * detwei * inv_den_gi( :, 1 ) )
                     end do
                     call addto( matrix, 1, 1, ph_inod, ph_jnod, nxnx )
                  end do
                  do iphase = 1, nphase
                     do idim = 1, ndim
                        call addto( rhs, ph_inod, &
                             sum( phfenx_all( idim, ph_iloc, : ) * inv_den_gi( :, iphase ) * ( &
                             u_s_gi( :, idim, iphase ) - coef_alpha_gi( :, iphase ) * &
                             dx_alpha_gi( :, idim, iphase ) ) * detwei ) )
                     end do
                  end do
               end do

            else

               ! form rhs of the momentum eqn...
               do u_iloc = 1, u_nloc
                  u_inod = u_ndgln( ( ele - 1 ) * u_nloc + u_iloc )
                  do iphase = 1, nphase
                     do idim = 1, ndim
                        u_rhs( idim, iphase, u_inod ) = u_rhs( idim, iphase, u_inod ) + &
                             sum( ph_funs%ufen( u_iloc, : ) * ( - dx_ph_gi( :, idim, iphase ) &
                             + u_s_gi( :, idim, iphase ) - coef_alpha_gi( :, iphase ) * &
                             dx_alpha_gi( :, idim, iphase ) ) * detwei )
                     end do
                  end do
               end do

            end if

         end do ! ele loop

         if ( iloop == 1 ) then

            got_free_surf = .false.
            pfield => extract_tensor_field( packed_state, "PackedFEPressure" )
            do i = 1, get_boundary_condition_count( pfield )
               call get_boundary_condition( pfield, i, type=bc_type, surface_node_list=surface_node_list )
               if ( trim( bc_type ) == "freesurface" .or.  trim( bc_type ) == "top" ) then
                  got_free_surf = .true.
                  exit
               end if
            end do

            ! if free surface apply a boundary condition
            ! else don't forget to remove the null space
            if ( got_free_surf ) then
               findph => sparsity % findrm
               colph => sparsity % colm
               do inod = 1, size( surface_node_list )
                  ph_inod = surface_node_list( inod )
                  rhs % val( ph_inod ) = 0.0
                  do count = findph( ph_inod ), findph( ph_inod + 1 ) - 1
                     ph_jnod = colph( count )
                     if ( ph_jnod /= ph_inod ) then
                        i = matrix % row_numbering % gnn2unn( ph_inod, 1 )
                        j = matrix % column_numbering % gnn2unn( ph_jnod, 1 )
                        call MatSetValue( matrix % m, i, j, 0.0, INSERT_VALUES, ierr )
                        do count2 = findph( ph_jnod ), findph( ph_jnod + 1 ) - 1
                           ph_jnod2 = colph( count2 )
                           if ( ph_jnod2 == ph_inod ) then
                              i = matrix % row_numbering % gnn2unn( ph_jnod, 1 )
                              j = matrix % column_numbering % gnn2unn( ph_jnod2, 1 )
                              call MatSetValue( matrix % m, i, j, 0.0, INSERT_VALUES, ierr )
                           end if
                        end do
                     end if
                  end do
               end do
            end if

            ! solver for pressure ph
            call set_solver_options( path, &
                 ksptype = "cg", &
                 pctype = "hypre", &
                 !ksptype = "gmres", &
                 !pctype = "sor", &
                 !ksptype = "preonly", &
                 !pctype = "lu", &
                 rtol = 1.0e-10, &
                 atol = 0.0, &
                 max_its = 10000 )
            call add_option( &
                 trim( path ) // "/solver/preconditioner[0]/hypre_type[0]/name", stat )
            call set_option( &
                 trim( path ) // "/solver/preconditioner[0]/hypre_type[0]/name", "boomeramg" )
            !call add_option( &
            !     trim( path ) // "/solver/preconditioner[0]/factorization_package[0]/name", stat )
            !call set_option( &
            !     trim( path ) // "/solver/preconditioner[0]/factorization_package/name", "petsc" )

            if ( .not.got_free_surf ) call add_option( &
                    trim( path ) // "/solver/remove_null_space", stat )
            ph_sol % option_path = path

            call petsc_solve( ph_sol, matrix, rhs )

            printf => extract_scalar_field( state( 1 ), "Ph", stat )
            if ( stat == 0 ) printf%val = ph_sol%val

            do iphase = 1, nphase
               ph( iphase, : ) = ph_sol % val ! assume 1 phase for the time being
            end do

         end if

      end do


      ! deallocate
      call deallocate( rhs )
      call deallocate( ph_sol )
      call deallocate( matrix )
      deallocate( u_ph_source_vel, u_ph_source_cv, alpha_cv, &
                 coef_alpha_cv, u_ph_source_ph, alpha_ph, &
                 ph, coef_alpha_ph, dx_ph_gi, u_s_gi, &
                 dx_alpha_gi, coef_alpha_gi, den_gi, inv_den_gi, sigma_gi, &
                 tmp_cvfenx_all, other_fenx_all, detwei, ra )
      ewrite(3,*) "leaving high_order_pressure_solve"

      return
    end subroutine high_order_pressure_solve


    REAL FUNCTION dg_oscilat_detect(SNDOTQ_KEEP, SNDOTQ2_KEEP, &
    N_DOT_DU, N_DOT_DU2, SINCOME, MASS_ELE, MASS_ELE2 )
        ! Determine if we have an oscillation in the normal direction...
        ! dg_oscilat_detect=1.0- CENTRAL SCHEME.
        ! dg_oscilat_detect=0.0- UPWIND SCHEME.
        real SNDOTQ_KEEP, SNDOTQ2_KEEP, N_DOT_DU, N_DOT_DU2, SINCOME
        REAL MASS_ELE, MASS_ELE2
        ! If cons_oscillation then apply upwinding as often as possible...
        LOGICAL, PARAMETER :: cons_oscillation = .false.
        REAL H1,H2, U1,U2,U3
        !              REAL TOLFUN

        if(cons_oscillation) then

            dg_oscilat_detect = 1.0

            if( SINCOME> 0.5 ) then
                ! velcity comming into element ELE...
                if( (SNDOTQ_KEEP - SNDOTQ2_KEEP)*N_DOT_DU2 > 0.0 ) dg_oscilat_detect = 0.0
            !                   if( (SNDOTQ_KEEP - SNDOTQ2_KEEP)*N_DOT_DU2 > 0.0 ) dg_oscilat_detect = 0.333
            !                   if( (SNDOTQ_KEEP - SNDOTQ2_KEEP)*N_DOT_DU2 > 0.0 ) dg_oscilat_detect = 0.5
            else
                ! velcity pointing out of the element ELE...
                if( (SNDOTQ2_KEEP - SNDOTQ_KEEP)*N_DOT_DU < 0.0 ) dg_oscilat_detect = 0.0
            !                   if( (SNDOTQ2_KEEP - SNDOTQ_KEEP)*N_DOT_DU < 0.0 ) dg_oscilat_detect = 0.333
            !                   if( (SNDOTQ2_KEEP - SNDOTQ_KEEP)*N_DOT_DU < 0.0 ) dg_oscilat_detect = 0.5
            end if
        else
            ! tvd in the means...

            dg_oscilat_detect = 1.0

            if( SINCOME> 0.5 ) then
                ! velcity comming into element ELE...
                if( (SNDOTQ_KEEP - SNDOTQ2_KEEP)*N_DOT_DU2 > 0.0 ) then
                    H1=MASS_ELE
                    H2=MASS_ELE2
                    U1=SNDOTQ_KEEP - N_DOT_DU* H1
                    U2=SNDOTQ2_KEEP+ N_DOT_DU2* H2
                    U3=SNDOTQ2_KEEP+ N_DOT_DU2* 3*H2
                    ! Have oscillations...
                    IF( (U1-U2)/TOLFUN(U2-U3) .LE. 0.0) dg_oscilat_detect = 0.0
                endif
            else
                ! velcity pointing out of the element ELE...
                if( (SNDOTQ2_KEEP - SNDOTQ_KEEP)*N_DOT_DU < 0.0 ) then
                    H1=MASS_ELE
                    H2=MASS_ELE2
                    U1=SNDOTQ_KEEP - N_DOT_DU* 3.*H1
                    U2=SNDOTQ_KEEP - N_DOT_DU* H1
                    U3=SNDOTQ2_KEEP+ N_DOT_DU2* H2
                    ! Have oscillations...
                    IF( (U1-U2)/TOLFUN(U2-U3) .LE. 0.0) dg_oscilat_detect = 0.0
                endif
            end if

        endif

        return
    end function dg_oscilat_detect

    SUBROUTINE DIFFUS_CAL_COEFF_STRESS_OR_TENSOR( Mdims, DIFF_COEF_DIVDX, &
        DIFF_COEFOLD_DIVDX, STRESS_FORM, STRESS_FORM_STAB, ZERO_OR_TWO_THIRDS, &
        SBUFEN_REVERSED,SBCVFEN_REVERSED,SBCVNGI, SLOC_UDIFFUSION, SLOC_UDIFFUSION_VOL, SLOC2_UDIFFUSION, SLOC2_UDIFFUSION_VOL, DIFF_GI_ADDED, &
        HDC, &
        U_CV_NODJ_IPHA_ALL, U_CV_NODI_IPHA_ALL, &
        UOLD_CV_NODJ_IPHA_ALL, UOLD_CV_NODI_IPHA_ALL, &
        ELE, ELE2, SNORMXN_ALL, &
        SLOC_DUX_ELE_ALL, SLOC2_DUX_ELE_ALL,   SLOC_DUOLDX_ELE_ALL, SLOC2_DUOLDX_ELE_ALL,  &
        SELE, WIC_U_BC, WIC_U_BC_DIRICHLET, SIMPLE_DIFF_CALC, DIFF_MIN_FRAC, DIFF_MAX_FRAC  )
        ! This sub calculates the effective diffusion coefficientd DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX
        ! based on a non-linear method and a non-oscillating scheme.
        ! This implements the stress and tensor form of diffusion and calculates a jump conidition.
        ! which is in DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX
        ! The coefficient are in N_DOT_DKDU, N_DOT_DKDUOLD.
        ! look at the manual DG treatment of viscocity.
        IMPLICIT NONE
        type(multi_dimensions), intent(in) :: Mdims
        LOGICAL, intent( in ) :: STRESS_FORM, STRESS_FORM_STAB, SIMPLE_DIFF_CALC
        INTEGER, intent( in ) :: SBCVNGI, ELE, ELE2, &
            &                   SELE, WIC_U_BC_DIRICHLET
        REAL, intent( in ) :: HDC, DIFF_MIN_FRAC, DIFF_MAX_FRAC
        REAL, DIMENSION(Mdims%ndim,Mdims%nphase,SBCVNGI), intent( in ) :: U_CV_NODJ_IPHA_ALL, U_CV_NODI_IPHA_ALL, &
            UOLD_CV_NODJ_IPHA_ALL, UOLD_CV_NODI_IPHA_ALL
        REAL, intent( in ) :: ZERO_OR_TWO_THIRDS
        REAL, DIMENSION( Mdims%ndim,Mdims%nphase,SBCVNGI ), intent( inout ) :: DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX
        INTEGER, DIMENSION( Mdims%ndim,Mdims%nphase,Mdims%stotel ), intent( in ) ::WIC_U_BC
        REAL, DIMENSION(  SBCVNGI, Mdims%cv_snloc ), intent( in ) :: SBCVFEN_REVERSED
        REAL, DIMENSION( SBCVNGI, Mdims%u_snloc ), intent( in ) :: SBUFEN_REVERSED
        REAL, DIMENSION( Mdims%ndim,Mdims%ndim,Mdims%nphase,Mdims%cv_snloc ), intent( in ) :: SLOC_UDIFFUSION, SLOC2_UDIFFUSION
        REAL, DIMENSION( Mdims%nphase,Mdims%cv_snloc ), intent( in ) :: SLOC_UDIFFUSION_VOL, SLOC2_UDIFFUSION_VOL
        ! DIFF_GI_ADDED( IDIM, :,:) is for dimension IDIM e.g IDIM=1 corresponds to U
        ! the rest is for the diffusion tensor.
        REAL, DIMENSION( Mdims%ndim, Mdims%ndim,Mdims%ndim, Mdims%nphase, SBCVNGI), intent( in ) :: DIFF_GI_ADDED
        REAL, DIMENSION( Mdims%ndim, Mdims%ndim , Mdims%nphase, Mdims%u_snloc ), intent( in ) :: SLOC_DUX_ELE_ALL, SLOC2_DUX_ELE_ALL,   SLOC_DUOLDX_ELE_ALL, SLOC2_DUOLDX_ELE_ALL
        REAL, DIMENSION( Mdims%ndim, SBCVNGI ), intent( in ) :: SNORMXN_ALL
        ! local variables
        !        ===>  REALS  <===
        ! DIFF_MIN_FRAC is the fraction of the standard diffusion coefficient to use
        ! in the non-linear diffusion scheme. DIFF_MAX_FRAC is the maximum fraction.
        ! If SIMPLE_DIFF_CALC then use a simple and fast diffusion calculation.
        !    LOGICAL, PARAMETER :: SIMPLE_DIFF_CALC2 = .false.
        !REAL, PARAMETER :: DIFF_MIN_FRAC = 0.005, DIFF_MAX_FRAC = 200.0
        !    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.1, DIFF_MAX_FRAC = 1000.0 ! works well but oscillations
        !    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.5, DIFF_MAX_FRAC = 1000000.0 ! works well no oscillations
        !    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.25, DIFF_MAX_FRAC = 100.0 ! works well no oscillations
        !    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.01, DIFF_MAX_FRAC = 100.0 ! works well no oscillations
        !    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.2, DIFF_MAX_FRAC = 100.0 ! works well no oscillations  ****recommended*****
        !    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.05, DIFF_MAX_FRAC = 200.0 ! works well no oscillations
        !    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.1, DIFF_MAX_FRAC = 200.0 ! works well no oscillations
        !    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.25, DIFF_MAX_FRAC = 10000000.0 ! works well no oscillations
        !    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.25, DIFF_MAX_FRAC = 1000.0
        REAL, DIMENSION( : , :, :, : ), allocatable :: DIFF_GI, DIFF_GI2, DIFF_GI_BOTH
        REAL, DIMENSION( : , : ), allocatable :: DIFF_VOL_GI, DIFF_VOL_GI2, DIFF_VOL_GI_BOTH
        REAL, DIMENSION( :, :, : ), allocatable :: N_DOT_DKDU, N_DOT_DKDUOLD, N_DOT_DKDU2, N_DOT_DKDUOLD2
        REAL, DIMENSION( :, :, : ), allocatable :: DIFF_STAND_DIVDX_U, DIFF_STAND_DIVDX2_U, &
            DIFF_COEF_DIVDX_U, DIFF_COEFOLD_DIVDX_U
        REAL, DIMENSION( :, : ), allocatable :: IDENT, RZER_DIFF_ALL
        REAL :: COEF
        INTEGER :: IDIM,JDIM,CV_SKLOC
        INTEGER :: SGI,IPHASE
        LOGICAL :: ZER_DIFF
        !    SIMPLE_DIFF_CALC=SIMPLE_DIFF_CALC2
        ALLOCATE( RZER_DIFF_ALL(Mdims%ndim,Mdims%nphase) )
        ZER_DIFF=.FALSE.
        RZER_DIFF_ALL=1.0
        IF(SELE /= 0) THEN
            ZER_DIFF=.TRUE.
            RZER_DIFF_ALL=0.0
            DO IPHASE=1,Mdims%nphase
                DO IDIM = 1, Mdims%ndim
                    IF(WIC_U_BC( IDIM, IPHASE, SELE) == WIC_U_BC_DIRICHLET) THEN
                        ZER_DIFF=.FALSE.
                        RZER_DIFF_ALL(IDIM,IPHASE)=1.0
                    ENDIF
                END DO
            END DO
        ENDIF
        Cond_ZerDiff: IF(ZER_DIFF) THEN
            DIFF_COEF_DIVDX    = 0.0
            DIFF_COEFOLD_DIVDX = 0.0
        ELSE
            ALLOCATE( N_DOT_DKDU( Mdims%ndim,Mdims%nphase,SBCVNGI )  )
            ALLOCATE( N_DOT_DKDUOLD( Mdims%ndim,Mdims%nphase,SBCVNGI )  )
            ALLOCATE( N_DOT_DKDU2( Mdims%ndim,Mdims%nphase,SBCVNGI )  )
            ALLOCATE( N_DOT_DKDUOLD2( Mdims%ndim,Mdims%nphase,SBCVNGI )  )
            ALLOCATE( DIFF_STAND_DIVDX_U( Mdims%ndim,Mdims%nphase,SBCVNGI )  )
            ALLOCATE( DIFF_STAND_DIVDX2_U( Mdims%ndim,Mdims%nphase,SBCVNGI )  )
            ALLOCATE( DIFF_COEF_DIVDX_U( Mdims%ndim,Mdims%nphase,SBCVNGI )  )
            ALLOCATE( DIFF_COEFOLD_DIVDX_U( Mdims%ndim,Mdims%nphase,SBCVNGI )  )
            IF(SIMPLE_DIFF_CALC) THEN ! The simplest method we can think of...
                ALLOCATE( DIFF_GI(Mdims%ndim,Mdims%ndim,Mdims%nphase,SBCVNGI) )
                ALLOCATE( DIFF_GI2(Mdims%ndim,Mdims%ndim,Mdims%nphase,SBCVNGI) )
                ALLOCATE( DIFF_GI_BOTH(Mdims%ndim,Mdims%ndim,Mdims%nphase,SBCVNGI) )
                ALLOCATE( DIFF_VOL_GI(Mdims%nphase,SBCVNGI) )
                ALLOCATE( DIFF_VOL_GI2(Mdims%nphase,SBCVNGI) )
                ALLOCATE( DIFF_VOL_GI_BOTH(Mdims%nphase,SBCVNGI) )
                ALLOCATE( IDENT(Mdims%ndim,Mdims%ndim) )
                IDENT=0.0
                DO IDIM=1,Mdims%ndim
                    IDENT(IDIM,IDIM)=1.0
                END DO
                DIFF_GI = 0.0
                DIFF_VOL_GI = 0.0
                DO CV_SKLOC = 1, Mdims%cv_snloc
                    DO SGI=1,SBCVNGI
                        DO IPHASE=1, Mdims%nphase
                            DIFF_GI( 1:Mdims%ndim , 1:Mdims%ndim, IPHASE, SGI ) = DIFF_GI( 1:Mdims%ndim , 1:Mdims%ndim, IPHASE, SGI ) &
                                + SBCVFEN_REVERSED(SGI,CV_SKLOC) * SLOC_UDIFFUSION( 1:Mdims%ndim , 1:Mdims%ndim , IPHASE, CV_SKLOC )
                            DIFF_VOL_GI( IPHASE, SGI ) = DIFF_VOL_GI( IPHASE, SGI ) &
                                + SBCVFEN_REVERSED(SGI,CV_SKLOC) * SLOC_UDIFFUSION_VOL( IPHASE, CV_SKLOC )
                        END DO
                    END DO
                END DO
                DIFF_GI=MAX(0.0, DIFF_GI)
                DIFF_VOL_GI=MAX(0.0, DIFF_VOL_GI)
                Conditional_MAT_DISOPT_ELE2_2: IF( ( ELE2 /= 0 ).AND.( ELE2 /= ELE) ) THEN
                    DIFF_GI2 = 0.0
                    DIFF_VOL_GI2 = 0.0
                    DO CV_SKLOC = 1, Mdims%cv_snloc
                        DO SGI=1,SBCVNGI
                            DO IPHASE=1, Mdims%nphase
                                DIFF_GI2( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, SGI )= DIFF_GI2( 1:Mdims%ndim, 1:Mdims%ndim, IPHASE, SGI ) +SBCVFEN_REVERSED(SGI,CV_SKLOC) &
                                    *SLOC2_UDIFFUSION(1:Mdims%ndim, 1:Mdims%ndim ,IPHASE, CV_SKLOC)
                                DIFF_VOL_GI2( IPHASE, SGI )= DIFF_VOL_GI2( IPHASE, SGI ) +SBCVFEN_REVERSED(SGI,CV_SKLOC) &
                                    *SLOC2_UDIFFUSION_VOL(IPHASE, CV_SKLOC)
                            END DO
                        END DO
                    END DO
                    DIFF_GI2=MAX(0.0, DIFF_GI2)
                    DIFF_VOL_GI2=MAX(0.0, DIFF_VOL_GI2)
                    DIFF_GI=0.5*(DIFF_GI+DIFF_GI2)
                    DIFF_VOL_GI=0.5*(DIFF_VOL_GI+DIFF_VOL_GI2)
                ENDIF Conditional_MAT_DISOPT_ELE2_2
                IF(STRESS_FORM) THEN
                    IF(STRESS_FORM_STAB) THEN
                        DIFF_GI_BOTH = DIFF_GI
                        DIFF_VOL_GI_BOTH = DIFF_VOL_GI
                        DO JDIM=1,Mdims%ndim
                            DO IDIM=1,Mdims%ndim
                                DIFF_GI_BOTH(IDIM, JDIM, :, :) = DIFF_GI_BOTH(IDIM, JDIM, :, :) &
                                    + SQRT( DIFF_GI_ADDED(IDIM, 1,1, :, :) * DIFF_GI_ADDED(JDIM, 1,1, :, :) )
                            END DO
                        END DO
                        DO SGI=1,SBCVNGI
                            DO IPHASE=1, Mdims%nphase
                                DO IDIM=1, Mdims%ndim
                                    DIFF_COEF_DIVDX(IDIM,IPHASE,SGI)=8.* SUM( (1.+IDENT(IDIM,:))*SNORMXN_ALL(:,SGI)**2*(DIFF_GI_BOTH(IDIM,:,IPHASE,SGI)+DIFF_VOL_GI_BOTH(IPHASE,SGI)) ) /HDC
                                END DO
                            END DO
                        END DO
                    ELSE
                        DO SGI=1,SBCVNGI
                            DO IPHASE=1, Mdims%nphase
                                DO IDIM=1, Mdims%ndim
                                    DIFF_COEF_DIVDX(IDIM,IPHASE,SGI)=8.*( SUM( (1.+IDENT(IDIM,:))*SNORMXN_ALL(:,SGI)**2*(DIFF_GI(IDIM,:,IPHASE,SGI)+DIFF_VOL_GI(IPHASE,SGI)) ) &
                                        +DIFF_GI_ADDED(IDIM, 1,1, IPHASE,SGI) ) /HDC
                                END DO
                            END DO
                        END DO
                    ENDIF
                ELSE
                    DO SGI=1,SBCVNGI
                        DO IPHASE=1, Mdims%nphase
                            COEF=0.0
                            DO IDIM=1,Mdims%ndim
                                COEF=COEF + SNORMXN_ALL(IDIM,SGI)*( SUM( DIFF_GI(IDIM,:,IPHASE,SGI)*SNORMXN_ALL(:,SGI) )  )
                            END DO
                            DIFF_COEF_DIVDX(:,IPHASE,SGI)=8.*( COEF + DIFF_GI_ADDED(:, 1,1, IPHASE,SGI) ) /HDC
                        END DO
                    END DO
                ENDIF
                DIFF_COEFOLD_DIVDX=DIFF_COEF_DIVDX
               ! END OF IF(SIMPLE_DIFF_CALC) THEN...
            ELSE
                ! Calculate DIFF_COEF_DIVDX, N_DOT_DKDU, N_DOT_DKDUOLD
                CALL FOR_TENS_DERIVS_NDOTS(DIFF_STAND_DIVDX_U, N_DOT_DKDU, N_DOT_DKDUOLD,  &
                    DIFF_GI_ADDED, SLOC_DUX_ELE_ALL, SLOC_DUOLDX_ELE_ALL, SLOC_UDIFFUSION, SLOC_UDIFFUSION_VOL, &
                    !  Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%u_snloc, SBCVNGI, SBCVFEN, SNORMXN_ALL, HDC, ZERO_OR_TWO_THIRDS, STRESS_FORM )
                    Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%u_snloc, Mdims%cv_snloc, SBCVNGI, SBUFEN_REVERSED, SBCVFEN_REVERSED, SNORMXN_ALL, HDC, ZERO_OR_TWO_THIRDS, &
                    STRESS_FORM, STRESS_FORM_STAB )
                Conditional_MAT_DISOPT_ELE2: IF( ( ELE2 /= 0 ).AND.( ELE2 /= ELE) ) THEN
                    ! Calculate DIFF_COEF_DIVDX, N_DOT_DKDU, N_DOT_DKDUOLD
                    CALL FOR_TENS_DERIVS_NDOTS(DIFF_STAND_DIVDX2_U, N_DOT_DKDU2, N_DOT_DKDUOLD2,  &
                        DIFF_GI_ADDED, SLOC2_DUX_ELE_ALL, SLOC2_DUOLDX_ELE_ALL, SLOC2_UDIFFUSION, SLOC2_UDIFFUSION_VOL, &
                        Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%u_snloc, Mdims%cv_snloc, SBCVNGI, SBUFEN_REVERSED, SBCVFEN_REVERSED, SNORMXN_ALL, HDC, ZERO_OR_TWO_THIRDS, &
                        STRESS_FORM, STRESS_FORM_STAB )
                    N_DOT_DKDU = 0.5*( N_DOT_DKDU + N_DOT_DKDU2 )
                    N_DOT_DKDUOLD= 0.5*( N_DOT_DKDUOLD + N_DOT_DKDUOLD2 )
                    ! This is the minimum diffusion...
                    DIFF_STAND_DIVDX_U    = 0.5*( DIFF_STAND_DIVDX_U + DIFF_STAND_DIVDX2_U )
                ENDIF Conditional_MAT_DISOPT_ELE2
                DO SGI=1,SBCVNGI
                    DO IPHASE=1, Mdims%nphase
                        DO IDIM=1,Mdims%ndim
                            DIFF_COEF_DIVDX_U(IDIM,IPHASE,SGI)    = N_DOT_DKDU(IDIM,IPHASE,SGI) / &
                                TOLFUN( U_CV_NODJ_IPHA_ALL(IDIM,IPHASE,SGI)  - U_CV_NODI_IPHA_ALL(IDIM,IPHASE,SGI) )
                            DIFF_COEFOLD_DIVDX_U(IDIM,IPHASE,SGI) = N_DOT_DKDUOLD(IDIM,IPHASE,SGI) /  &
                                TOLFUN( UOLD_CV_NODJ_IPHA_ALL(IDIM,IPHASE,SGI)  - UOLD_CV_NODI_IPHA_ALL(IDIM,IPHASE,SGI) )
                        END DO
                    END DO
                END DO
                ! Make sure the diffusion has an lower bound...
                DIFF_COEF_DIVDX_U    = MAX( DIFF_MIN_FRAC*DIFF_STAND_DIVDX_U, DIFF_COEF_DIVDX_U )
                DIFF_COEFOLD_DIVDX_U = MAX( DIFF_MIN_FRAC*DIFF_STAND_DIVDX_U, DIFF_COEFOLD_DIVDX_U )
                ! Make sure the diffusion has an upper bound...
                DIFF_COEF_DIVDX_U    = MIN( DIFF_MAX_FRAC*DIFF_STAND_DIVDX_U, DIFF_COEF_DIVDX_U )
                DIFF_COEFOLD_DIVDX_U = MIN( DIFF_MAX_FRAC*DIFF_STAND_DIVDX_U, DIFF_COEFOLD_DIVDX_U )
                ! Redfine for output...
                DIFF_COEF_DIVDX   = DIFF_COEF_DIVDX_U
                DIFF_COEFOLD_DIVDX = DIFF_COEFOLD_DIVDX_U
               ! END OF IF(SIMPLE_DIFF_CALC) THEN ELSE...
            ENDIF
        END IF Cond_ZerDiff
        !
        ! Zero if we are on boundary and applying Dirichlet b.c's
        DO IPHASE=1, Mdims%nphase
            DO IDIM=1, Mdims%ndim
                DIFF_COEF_DIVDX(IDIM,IPHASE,:)    =  RZER_DIFF_ALL(IDIM,IPHASE)*DIFF_COEF_DIVDX(IDIM,IPHASE,:)
                DIFF_COEFOLD_DIVDX(IDIM,IPHASE,:) =  RZER_DIFF_ALL(IDIM,IPHASE)*DIFF_COEFOLD_DIVDX(IDIM,IPHASE,:)
            END DO
        END DO


        contains

            SUBROUTINE FOR_TENS_DERIVS_NDOTS( DIFF_STAND_DIVDX_U, N_DOT_DKDU, N_DOT_DKDUOLD,  &
                DIFF_GI_ADDED, SLOC_DUX_ELE_ALL, SLOC_DUOLDX_ELE_ALL, SLOC_UDIFFUSION, SLOC_UDIFFUSION_VOL, &
                NDIM_VEL, NDIM, NPHASE, U_SNLOC, CV_SNLOC, SBCVNGI, SBUFEN_REVERSED, SBCVFEN_REVERSED, SNORMXN_ALL, HDC, ZERO_OR_TWO_THIRDS, &
                STRESS_FORM, STRESS_FORM_STAB )

                ! Calculate DIFF_STAND_DIVDX_U, N_DOT_DKDU, N_DOT_DKDUOLD
                ! This implements the stress and tensor form of diffusion and calculates a jump conidition.
                ! DIFF_STAND_DIVDX_U is the minimal amount of diffusion.
                ! The coefficient are in N_DOT_DKDU, N_DOT_DKDUOLD.
                ! look at the manual DG treatment of viscocity.
                IMPLICIT NONE
                INTEGER, intent( in )  :: NDIM_VEL, NDIM, NPHASE, U_SNLOC, CV_SNLOC, SBCVNGI
                REAL, intent( in )  :: HDC, ZERO_OR_TWO_THIRDS
                LOGICAL, intent( in )  :: STRESS_FORM, STRESS_FORM_STAB
                REAL, DIMENSION( NDIM,NPHASE,SBCVNGI ), intent( inout ) :: DIFF_STAND_DIVDX_U
                REAL, DIMENSION( NDIM_VEL,NPHASE,SBCVNGI ), intent( inout ) :: N_DOT_DKDU, N_DOT_DKDUOLD
                ! DIFF_GI_ADDED( IDIM, :,:) is for dimension IDIM e.g IDIM=1 corresponds to U
                ! the rest is for the diffusion tensor.
                REAL, DIMENSION( NDIM_VEL, NDIM,NDIM, NPHASE, SBCVNGI), intent( in ) :: DIFF_GI_ADDED
                REAL, DIMENSION( NDIM_VEL, NDIM , NPHASE, U_SNLOC ), intent( in ) :: SLOC_DUX_ELE_ALL, SLOC_DUOLDX_ELE_ALL
                REAL, DIMENSION( SBCVNGI, U_SNLOC ), intent( in ) :: SBUFEN_REVERSED
                REAL, DIMENSION( SBCVNGI, CV_SNLOC ), intent( in ) :: SBCVFEN_REVERSED
                REAL, DIMENSION( NDIM,NDIM,NPHASE,CV_SNLOC ), intent( in ) :: SLOC_UDIFFUSION
                REAL, DIMENSION( NPHASE,CV_SNLOC ), intent( in ) :: SLOC_UDIFFUSION_VOL
                REAL, DIMENSION( NDIM, SBCVNGI ), intent( in ) :: SNORMXN_ALL

                ! local variables
                REAL, DIMENSION( : , :, :, : ), allocatable :: DIFF_GI, STRESS_INDEX, STRESS_INDEXOLD
                REAL, DIMENSION( : , :, :, : ), allocatable :: DIFF_GI_BOTH
                REAL, DIMENSION( : , : ), allocatable :: DIFF_VOL_GI, DIFF_VOL_GI_BOTH
                REAL, DIMENSION( :, :, :, : ), allocatable :: DUDX_ALL_GI, DUOLDDX_ALL_GI
                REAL, DIMENSION( :, : ), allocatable :: IDENT
                REAL :: DIVU, DIVUOLD
                INTEGER :: IDIM,JDIM,IDIM_VEL,U_SKLOC,CV_SKLOC
                INTEGER :: SGI,IPHASE


                ALLOCATE( DIFF_GI(NDIM,NDIM,NPHASE,SBCVNGI) )
                ALLOCATE( DIFF_VOL_GI(NPHASE,SBCVNGI) )

                ALLOCATE( STRESS_INDEX(NDIM,NDIM,NPHASE,SBCVNGI) )
                ALLOCATE( STRESS_INDEXOLD(NDIM,NDIM,NPHASE,SBCVNGI) )

                ALLOCATE( DIFF_GI_BOTH(NDIM,NDIM, NPHASE,SBCVNGI) )
                ALLOCATE( DIFF_VOL_GI_BOTH(NPHASE,SBCVNGI) )

                ALLOCATE( DUDX_ALL_GI( NDIM_VEL,NDIM,NPHASE,SBCVNGI )  )
                ALLOCATE( DUOLDDX_ALL_GI( NDIM_VEL,NDIM,NPHASE,SBCVNGI )  )

                ALLOCATE( IDENT(NDIM,NDIM) )


                IDENT=0.0
                DO IDIM=1,NDIM
                    IDENT(IDIM,IDIM)=1.0
                END DO


                DUDX_ALL_GI = 0.0
                DUOLDDX_ALL_GI = 0.0

                DO U_SKLOC = 1, U_SNLOC
                    DO SGI=1,SBCVNGI
                        ! U, V & W:
                        DUDX_ALL_GI(:,:,:,SGI)    = DUDX_ALL_GI(:,:,:,SGI)    + SBUFEN_REVERSED(SGI,U_SKLOC) * SLOC_DUX_ELE_ALL(:,:,:,U_SKLOC)
                        DUOLDDX_ALL_GI(:,:,:,SGI) = DUOLDDX_ALL_GI(:,:,:,SGI) + SBUFEN_REVERSED(SGI,U_SKLOC) * SLOC_DUOLDX_ELE_ALL(:,:,:,U_SKLOC)
                    END DO
                END DO

                DIFF_GI = 0.0
                DIFF_VOL_GI = 0.0
                DO CV_SKLOC = 1, CV_SNLOC
                    DO SGI=1,SBCVNGI
                        DO IPHASE=1, NPHASE
                            DIFF_GI( 1:NDIM , 1:NDIM, IPHASE,SGI ) = DIFF_GI( 1:NDIM , 1:NDIM, IPHASE,SGI ) &
                                + SBCVFEN_REVERSED(SGI,CV_SKLOC) * SLOC_UDIFFUSION( 1:NDIM , 1:NDIM , IPHASE, CV_SKLOC )

                            DIFF_VOL_GI( IPHASE,SGI ) = DIFF_VOL_GI( IPHASE,SGI ) &
                                + SBCVFEN_REVERSED(SGI,CV_SKLOC) * SLOC_UDIFFUSION_VOL( IPHASE, CV_SKLOC )
                        END DO
                    END DO
                END DO
                DIFF_GI=MAX(0.0, DIFF_GI)
                DIFF_VOL_GI=MAX(0.0, DIFF_VOL_GI)


                IF(STRESS_FORM) THEN
                    ! FOR STRESS FORM...
                    ! BUT 1st tensor form for added diffusion from stabilization say...
                    N_DOT_DKDU=0.0
                    N_DOT_DKDUOLD=0.0
                    DIFF_STAND_DIVDX_U=0.0

                    DIFF_GI_BOTH = DIFF_GI
                    DIFF_VOL_GI_BOTH = DIFF_VOL_GI

                    IF(STRESS_FORM_STAB) THEN
                        DO JDIM=1,NDIM
                            DO IDIM=1,NDIM
                                DIFF_GI_BOTH(IDIM, JDIM, :, :) = DIFF_GI_BOTH(IDIM, JDIM, :, :) &
                                    + SQRT( DIFF_GI_ADDED(IDIM, 1,1, :, :) * DIFF_GI_ADDED(JDIM, 1,1, :, :) )
                            END DO
                        END DO
                    ELSE ! Tensor form
                        DO SGI=1,SBCVNGI
                            DO IPHASE=1, NPHASE
                                DO IDIM_VEL=1,NDIM_VEL
                                    DO IDIM=1,NDIM
                                        ! tensor form...
                                        N_DOT_DKDU(IDIM_VEL,IPHASE,SGI)   =  N_DOT_DKDU(IDIM_VEL,IPHASE,SGI)   &
                                            +  SNORMXN_ALL(IDIM,SGI)*SUM( DIFF_GI_ADDED(IDIM_VEL,IDIM,:,IPHASE,SGI) * DUDX_ALL_GI(IDIM_VEL,:,IPHASE,SGI) )
                                        ! tensor form...
                                        N_DOT_DKDUOLD(IDIM_VEL,IPHASE,SGI)= N_DOT_DKDUOLD(IDIM_VEL,IPHASE,SGI)  &
                                            +  SNORMXN_ALL(IDIM,SGI)*SUM( DIFF_GI_ADDED(IDIM_VEL,IDIM,:,IPHASE,SGI) * DUOLDDX_ALL_GI(IDIM_VEL,:,IPHASE,SGI) )
                                        ! for minimal amount of diffusion calc...
                                        DIFF_STAND_DIVDX_U(IDIM_VEL,IPHASE,SGI)   =  DIFF_STAND_DIVDX_U(IDIM_VEL,IPHASE,SGI)   &
                                            + SNORMXN_ALL(IDIM,SGI)*SUM( DIFF_GI_ADDED(IDIM_VEL,IDIM,:,IPHASE,SGI) * SNORMXN_ALL(:,SGI) )  /HDC

                                    END DO
                                END DO
                            END DO
                        END DO
                    ENDIF

                    ! stress form needs to add this...
                    DO SGI=1,SBCVNGI
                        DO IPHASE=1, NPHASE

                            DIVU=0.0
                            DIVUOLD=0.0
                            DO IDIM=1,NDIM
                                DIVU=DIVU+DUDX_ALL_GI(IDIM,IDIM,IPHASE,SGI)
                                DIVUOLD=DIVUOLD+DUOLDDX_ALL_GI(IDIM,IDIM,IPHASE,SGI)
                            END DO

                            DO IDIM_VEL=1,NDIM_VEL
                                ! Stress form...
                                N_DOT_DKDU(IDIM_VEL,IPHASE,SGI)   =  N_DOT_DKDU(IDIM_VEL,IPHASE,SGI) &
                                    + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI_BOTH(IDIM_VEL,:,IPHASE,SGI)*DUDX_ALL_GI(IDIM_VEL,:,IPHASE,SGI) )  &
                                    + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI_BOTH(IDIM_VEL,:,IPHASE,SGI)*DUDX_ALL_GI(:,IDIM_VEL,IPHASE,SGI) ) &
                                    ! stress form addition...
                                    - ZERO_OR_TWO_THIRDS*SNORMXN_ALL(IDIM_VEL,SGI)*DIFF_GI_BOTH(IDIM_VEL,IDIM_VEL,IPHASE,SGI)*DIVU &
                                    + SNORMXN_ALL(IDIM_VEL,SGI)*DIFF_VOL_GI_BOTH(IPHASE,SGI)*DIVU

                                ! Stress form...
                                N_DOT_DKDUOLD(IDIM_VEL,IPHASE,SGI)= N_DOT_DKDUOLD(IDIM_VEL,IPHASE,SGI) &
                                    + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI_BOTH(IDIM_VEL,:,IPHASE,SGI)*DUOLDDX_ALL_GI(IDIM_VEL,:,IPHASE,SGI) )  &
                                    + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI_BOTH(IDIM_VEL,:,IPHASE,SGI)*DUOLDDX_ALL_GI(:,IDIM_VEL,IPHASE,SGI) ) &
                                    ! stress form addition...
                                    - ZERO_OR_TWO_THIRDS*SNORMXN_ALL(IDIM_VEL,SGI)*DIFF_GI_BOTH(IDIM_VEL,IDIM_VEL,IPHASE,SGI)*DIVUOLD &
                                    + SNORMXN_ALL(IDIM_VEL,SGI)*DIFF_VOL_GI_BOTH(IPHASE,SGI)*DIVUOLD

                                ! This is for the minimum & max. diffusion...
                                DIFF_STAND_DIVDX_U(IDIM_VEL,IPHASE,SGI)   =  DIFF_STAND_DIVDX_U(IDIM_VEL,IPHASE,SGI) &
                                    + (   SUM( SNORMXN_ALL(:,SGI)*DIFF_GI_BOTH(IDIM_VEL,:,IPHASE,SGI)*SNORMXN_ALL(:,SGI) )  &
                                    +  SNORMXN_ALL(IDIM_VEL,SGI)*DIFF_GI_BOTH(IDIM_VEL,IDIM_VEL,IPHASE,SGI)*SNORMXN_ALL(IDIM_VEL,SGI)    &
                                    +  SNORMXN_ALL(IDIM_VEL,SGI)*DIFF_VOL_GI_BOTH(IPHASE,SGI)*SNORMXN_ALL(IDIM_VEL,SGI)     )/HDC
                            !                                                               + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI_BOTH(IDIM_VEL,:,IPHASE,SGI)*SNORMXN_ALL(IDIM_VEL,SGI) )    )/HDC

                            END DO
                        END DO
                    END DO


                ELSE  ! IF(STRESS_FORM) THEN ELSE
                    ! tensor form...
                    ! tensor form for added diffusion from stabilization as well...
                    N_DOT_DKDU=0.0
                    N_DOT_DKDUOLD=0.0
                    DIFF_STAND_DIVDX_U=0.0
                    DO SGI=1,SBCVNGI
                        DO IPHASE=1, NPHASE
                            DO IDIM=1,NDIM
                                DO IDIM_VEL=1,NDIM_VEL
                                    ! tensor form...
                                    N_DOT_DKDU(IDIM_VEL,IPHASE,SGI)   =  N_DOT_DKDU(IDIM_VEL,IPHASE,SGI)   &
                                        +  SNORMXN_ALL(IDIM,SGI)*SUM( (DIFF_GI_ADDED(IDIM_VEL,IDIM,:,IPHASE,SGI)+DIFF_GI(IDIM,:,IPHASE,SGI)) * DUDX_ALL_GI(IDIM_VEL,:,IPHASE,SGI) )
                                    ! tensor form...
                                    N_DOT_DKDUOLD(IDIM_VEL,IPHASE,SGI)= N_DOT_DKDUOLD(IDIM_VEL,IPHASE,SGI)  &
                                        +  SNORMXN_ALL(IDIM,SGI)*SUM( (DIFF_GI_ADDED(IDIM_VEL,IDIM,:,IPHASE,SGI)+DIFF_GI(IDIM,:,IPHASE,SGI)) * DUOLDDX_ALL_GI(IDIM_VEL,:,IPHASE,SGI) )
                                    ! This is for the minimum & max. diffusion...
                                    DIFF_STAND_DIVDX_U(IDIM_VEL,IPHASE,SGI)   =  DIFF_STAND_DIVDX_U(IDIM_VEL,IPHASE,SGI)   &
                                        +  SNORMXN_ALL(IDIM,SGI)*SUM( (DIFF_GI_ADDED(IDIM_VEL,IDIM,:,IPHASE,SGI)+DIFF_GI(IDIM,:,IPHASE,SGI)) * SNORMXN_ALL(:,SGI) )   /HDC

                                END DO
                            END DO
                        END DO
                    END DO


                   ! ENDOF IF(STRESS_FORM) THEN ELSE...
                ENDIF
                ! just in case...
                ! the factor of 8 is there to take into account that HD is measured between centres of elements...
                DIFF_STAND_DIVDX_U=abs( 8.*DIFF_STAND_DIVDX_U )
                !          DIFF_STAND_DIVDX_U=( 8.*DIFF_STAND_DIVDX_U )

                RETURN

            END SUBROUTINE FOR_TENS_DERIVS_NDOTS


    END SUBROUTINE DIFFUS_CAL_COEFF_STRESS_OR_TENSOR

 end module multiphase_1D_engine
