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
!> This module contain all the necessary subroutines to deal with FE-wise equations and fields. Assembling and solving of the momentum equation
!> and associated fields such as capillary pressure, hydrostatic pressure solver, and stabilisation techniques.
!> Also includes the assembly and solving of a Laplacian system (for zeta potential mainly)
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
    use multi_tools, only: CALC_FACE_ELE, tolfun
    use parallel_tools, only : allmax, allmin, isparallel
    use ieee_arithmetic
#ifdef USING_XGBOOST
    use multi_machine_learning
    use iso_c_binding
#endif

    implicit none

    private :: CV_ASSEMB_FORCE_CTY, get_diagonal_mass_matrix

    public  :: INTENERGE_ASSEM_SOLVE, VolumeFraction_Assemble_Solve, &
    generate_and_solve_Laplacian_system, Tracer_Assemble_Solve

contains
    !---------------------------------------------------------------------------
    !> @author Chris Pain, Pablo Salinas
    !> @brief Calls to generate the transport equation for the transport of energy/temperature and to solve the transport of components
    !>@param  state Linked list containing all the fields defined in diamond and considered by Fluidity
    !>@param  packed_state Linked list containing all the fields used by IC-FERST, memory partially shared with state
    !!>@param Mdims Data type storing all the dimensions describing the mesh, fields, nodes, etc
    !>@param  CV_GIdims Gauss integration numbers for CV fields
    !>@param  CV_funs Shape functions for the CV mesh
    !>@param  Mspars Sparsity of the matrices
    !>@param  ndgln Global to local variables
    !>@param  Mdisopt Discretisation options
    !>@param  Mmat Matrices for ICFERST
    !>@param  upwnd Sigmas to compute the fluxes at the interphase for porous media
    !>@param  tracer  Tracer considered for the transport equation
    !>@param  density  Density of the field
    !>@param  velocity  Velocity of the field
    !>@param  multi_absorp  Absoprtion of associated with the transport field
    !>@param  IGOT_T2 !>@param  IGOT_T2. True (1) if solving for a tracer, false otherwise
    !>@param   igot_theta_flux ????
    !>@param  GET_THETA_FLUX, USE_THETA_FLUX ?????
    !>@param THERMAL If true then we are solving for temperature
    !>@param  THETA_GDIFF ????
    !>@param  eles_with_pipe Elements that have a pipe
    !>@param  pipes_aux Information required to define wells
    !>@param  THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J ????
    !>@param  DT Time step size
    !>@param  SUF_SIG_DIAGTEN_BC Like upwnd but for the boundary
    !>@param  VOLFRA_PORE     Porosity field (Mdims%npres,Mdims%totele)
    !>@param  option_path   Option path of the tracer to be solved for 
    !>@param  mass_ele_transp Mass of the elements
    !>@param  saturation PhaseVolumeFraction field
    !>@param Permeability_tensor_field Permeability field
    !>@param  icomp Number of components, if > 0 then compositional is solved
    !>param nonlinear_iteration Current non-linear iteration
  SUBROUTINE INTENERGE_ASSEM_SOLVE( state, packed_state, &
       Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd,&
       tracer, velocity, density, multi_absorp, DT, &
       SUF_SIG_DIAGTEN_BC,  VOLFRA_PORE, &
       IGOT_T2, igot_theta_flux,GET_THETA_FLUX, USE_THETA_FLUX,  &
       THETA_GDIFF, eles_with_pipe, pipes_aux, &
       option_path, &
       mass_ele_transp, &
       thermal, THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, &
       icomp, saturation, Permeability_tensor_field, nonlinear_iteration )
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
           type(multi_absorption), intent(inout) :: multi_absorp
           INTEGER, intent( in ) :: IGOT_T2, igot_theta_flux
           LOGICAL, intent( in ) :: GET_THETA_FLUX, USE_THETA_FLUX
           LOGICAL, intent( in ), optional ::THERMAL
           REAL, DIMENSION( :, : ), intent( inout ) :: THETA_GDIFF
           REAL, DIMENSION( :,: ), intent( inout ), optional :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
           REAL, intent( in ) :: DT
           REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
           REAL, DIMENSION( :, : ), intent( in ) :: VOLFRA_PORE
           character( len = * ), intent( in ), optional :: option_path
           real, dimension( : ), intent( inout ), optional :: mass_ele_transp
           type(tensor_field), intent(in), optional :: saturation
           type( tensor_field ), optional, pointer, intent(in) :: Permeability_tensor_field
           integer, optional :: icomp, nonlinear_iteration
           type(pipe_coords), dimension(:), intent(in):: eles_with_pipe
           type (multi_pipe_package), intent(in) :: pipes_aux
           ! Local variables
           LOGICAL, PARAMETER :: GETCV_DISC = .TRUE., GETCT= .FALSE.
           integer :: nits_flux_lim, its_flux_lim
           REAL, DIMENSION( :, : ), allocatable :: DIAG_SCALE_PRES
           REAL, DIMENSION( :, :, : ), allocatable :: DIAG_SCALE_PRES_COUP, GAMMA_PRES_ABS, GAMMA_PRES_ABS_NANO, INV_B
           REAL, DIMENSION( Mdims%mat_nonods, Mdims%ndim, Mdims%ndim, Mdims%nphase ) :: TDIFFUSION
           REAL, DIMENSION( : ), ALLOCATABLE :: MASS_PIPE, MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE
           real, dimension( size(Mspars%small_acv%col )) ::  mass_mn_pres
           REAL, DIMENSION( : , : ), allocatable :: denold_all, t_source
           REAL, DIMENSION( : , : ), target, allocatable :: den_all
           REAL, DIMENSION( : ), allocatable :: CV_RHS_SUB
           type( tensor_field ), pointer :: P, Q
           INTEGER :: IPHASE, its_taken, ipres
           type( tensor_field ), pointer :: den_all2, denold_all2, a, aold, deriv, Component_Absorption
           type( vector_field ), pointer  :: MeanPoreCV, python_vfield
           integer :: lcomp, Field_selector, IGOT_T2_loc, python_stat, stat
           type(vector_field)  :: vtracer, residual
           type(csr_sparsity), pointer :: sparsity
           real, dimension(:,:,:), allocatable :: Velocity_Absorption
           real, dimension(:,:,:), pointer :: T_AbsorB=>null()
           integer :: ncomp_diff_coef, comp_diffusion_opt, nphase, n_in_pres, auxI
           real, dimension(:,:,:), allocatable :: Component_Diffusion_Operator_Coefficient
           type( tensor_field ), pointer :: perm, python_tfield
           integer :: cv_disopt, cv_dg_vel_int_opt
           real :: cv_theta, cv_beta
           type( scalar_field ), pointer :: sfield, porous_field, solid_concentration
           REAL, DIMENSION( : ), allocatable :: porous_heat_coef, porous_heat_coef_old
           character(len=option_path_len) :: solver_option_path = "/solver_options/Linear_solver"
           REAL, DIMENSION( :,:,:,: ), allocatable :: CDISPERSION
           !Variables to stabilize the non-linear iteration solver
           real, dimension(2) :: totally_min_max
           logical :: impose_min_max
           real :: aux
           real, save :: inf_tolerance = -1
           !Variables to control the PETCs solver
           integer, save :: max_allowed_its = -1
           !Variables for vanishing diffusion
           real, dimension(Mdims%cv_nonods) :: OvRelax_param
           integer :: Phase_with_Ovrel
           !temperature backup for the petsc bug
           type ( tensor_field ), pointer :: temperature
           real, dimension(Mdims%nphase, Mdims%cv_nonods) :: temp_bak
           logical :: repeat_assemb_solve, assemble_collapsed_to_one_phase
           type(vector_field) :: solution

           !Initialise with an out of range value to be able to check it hasn't been
           totally_min_max = 1e30
           lcomp = 0
           if ( present( icomp ) ) lcomp = icomp


           if (present(Permeability_tensor_field)) then
            perm => Permeability_tensor_field
           else
            perm=>extract_tensor_field(packed_state,"Permeability")
           end if

           sparsity=>extract_csr_sparsity(packed_state,"ACVSparsity")
           allocate(den_all(Mdims%nphase,Mdims%cv_nonods),denold_all(Mdims%nphase,Mdims%cv_nonods))

           allocate( T_SOURCE( Mdims%nphase, Mdims%cv_nonods ) ) ; T_SOURCE=0.0!SPRINT_TO_DO TURN THESE T_SOURCE INTO POINTERS OR DIRECTLY REMOVE THEM
           IGOT_T2_loc = 0

           assemble_collapsed_to_one_phase = .false.
           if ( thermal .or. trim( option_path ) == '/material_phase[0]/scalar_field::Temperature') then

               p => extract_tensor_field( packed_state, "PackedCVPressure", stat )
               if (stat/=0) p => extract_tensor_field( packed_state, "PackedFEPressure", stat )
                !If it is thermal and porous media we need to consider thermal equilibirum between phases and porous medium
                !in this case we need to solve then only for one temperature per region(reservoir/wells)
                assemble_collapsed_to_one_phase = .true.
                !Check that the extra parameters required for porous media thermal simulations are present
                if (.not.have_option('/porous_media/porous_properties/scalar_field::porous_density') .or. &
                    .not.have_option('/porous_media/porous_properties/scalar_field::porous_heat_capacity') .or. &
                    .not.have_option('/porous_media/porous_properties/tensor_field::porous_thermal_conductivity')) then
                    FLAbort("For thermal porous media flows the following fields are mandatory: porous_density, porous_heat_capacity and porous_thermal_conductivity ")
                end if
                !need to perform average of the effective heat capacity times density for the diffusion and time terms
                allocate(porous_heat_coef(Mdims%cv_nonods))
                allocate(porous_heat_coef_old(Mdims%cv_nonods))
                call effective_Cp_density(porous_heat_coef, porous_heat_coef_old)
               den_all2 => extract_tensor_field( packed_state, "PackedDensityHeatCapacity", stat )
               denold_all2 => extract_tensor_field( packed_state, "PackedOldDensityHeatCapacity", stat )
               if (stat /= 0) then
                 den_all2 => extract_tensor_field( packed_state, "PackedDensity", stat )
                 denold_all2 => extract_tensor_field( packed_state, "PackedOldDensity" )
               end if
               den_all    = den_all2 % val ( 1, :, : )
               denold_all = denold_all2 % val ( 1, :, : )
	       if(have_option( '/femdem_thermal/coupling/ring_and_volume') .OR. have_option( '/femdem_thermal/coupling/volume_relaxation') ) then
                   solid_concentration => extract_scalar_field( packed_state, "SolidConcentration" )
                   den_all( 1, : ) = den_all ( 1, : ) * (1.0 - solid_concentration % val)
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

           !Need to change this to use a reference density/rho_cp so for porous media the rock/fluid ratio is kept
           if (has_boussinesq_aprox) then
            do iphase = 1, Mdims%nphase
                !Retrieve CP (considered constant) to get rhoCp
                sfield => extract_scalar_field( state( iphase ), 'TemperatureHeatCapacity', stat )
                !If compositional then component Cp
                if (lcomp > 0) then
                sfield => extract_scalar_field( state( Mdims%nphase + lcomp ), 'ComponentMassFractionPhase' // int2str( iphase ) // 'HeatCapacity', stat )
                den_all((lcomp - 1 ) * Mdims%nphase + iphase,:) = sfield%val(1) * retrieve_reference_density(state, packed_state, iphase, lcomp, Mdims%nphase)
                else
                den_all(iphase,:) = sfield%val(1) * retrieve_reference_density(state, packed_state, iphase, lcomp, Mdims%nphase)
                end if
            end do
            !Copy to old to ensure no time variation
            denold_all = den_all
           end if
           if( present( option_path ) ) then ! solving for Temperature or Internal Energy or k_epsilon model

               if( trim( option_path ) == '/material_phase[0]/scalar_field::Temperature' ) then
                   call get_option( '/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation/' // &
                       'control_volumes/number_advection_iterations', nits_flux_lim, default = 3 )
                   Field_selector = 1
                   !Retrieve source term; sprint_to_do something equivalent should be done for absoprtion
                    do iphase = 1, Mdims%nphase
                      sfield => extract_scalar_field( state(iphase), "TemperatureSource", stat )
                      if (stat == 0) call assign_val(T_source( iphase, : ),sfield%val)
                    end do
               end if
               if (thermal) then
                   !We control with the infinite norm of the difference the non-linear iterations done in this sub-cycle
                   !therefore the minimum/default value of nits_flux_lim is set to 9
                   nits_flux_lim = max(nits_flux_lim, 9)!Currently overriden as we are not updating the rhs or other fields so this is not useful
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

           deriv => extract_tensor_field( packed_state, "PackedDRhoDPressure" )
           TDIFFUSION=0.0

           if ( thermal .or. trim( option_path ) == '/material_phase[0]/scalar_field::Temperature') then
                !For porous media thermal two fields are returned. Being one the diffusivity of the porous medium
                call calculate_diffusivity( state, packed_state, Mdims, ndgln, TDIFFUSION)
                !Calculates dispersion with specific longitudinal and transverse dispersivity
                if (have_option("/porous_media/Dispersion/scalar_field::Longitudinal_Dispersivity")) then
                  allocate(CDISPERSION(Mdims%mat_nonods, Mdims%ndim, Mdims%ndim, Mdims%nphase)); CDISPERSION = 0.
                  call calculate_solute_dispersity( state, packed_state, Mdims, ndgln, den_all, CDISPERSION)
                  TDIFFUSION = TDIFFUSION + CDISPERSION
                  deallocate(CDISPERSION)
                 end if
           end if

           ! get diffusivity for compositional
           if ( lcomp > 0) then
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

           ! Check for a python-set absorption field when solving for temperature/internal energy
           python_tfield => extract_tensor_field( state(1), "TAbsorB", python_stat )
           if (python_stat==0 .and. Field_selector==1) T_ABSORB = python_tfield%val

           ! Check for a python-set source field when solving for temperature/internal energy
           python_vfield => extract_vector_field( state(1), "TSourcE", python_stat )
           if (python_stat==0 .and. Field_selector==1) T_SOURCE = python_vfield%val
           !Start with the process to apply the min max principle
           call force_min_max_principle(Mdims, 1, tracer, nonlinear_iteration, totally_min_max)

           MeanPoreCV=>extract_vector_field(packed_state,"MeanPoreCV")
NITS_FLUX_LIM = 5!<= currently looping here more does not add anything as RHS and/or velocity are not updated
                !we set up 5 iterations but if it converges => we exit straigth away
temp_bak = tracer%val(1,:,:)!<= backup of the tracer field, just in case the petsc bug hits us here, we can retry

            if ( have_option( '/femdem_thermal/coupling') ) then
              Component_Absorption => extract_tensor_field( packed_state, "PackedTemperatureAbsorption")
              T_ABSORB(1:1,1:1,1:Mdims%cv_nonods)=> Component_Absorption%val (1,1,1:Mdims%cv_nonods)
            end if

           !Select solver options
           solver_option_path = "/solver_options/Linear_solver"
           IF ( IGOT_T2 == 1) THEN
             if (have_option('/solver_options/Linear_solver/Custom_solver_configuration/field::Compositional')) then
               solver_option_path = '/solver_options/Linear_solver/Custom_solver_configuration/field::Compositional'
             end if
           else
             if (have_option('/solver_options/Linear_solver/Custom_solver_configuration/field::Temperature')) then
               solver_option_path = '/solver_options/Linear_solver/Custom_solver_configuration/field::Temperature'
             end if
           end if
           if(max_allowed_its < 0)  then
               call get_option( trim(solver_option_path)//"max_iterations",&
                max_allowed_its, default = 500)
           end if

           !By default use the given values
           nphase = Mdims%nphase
           n_in_pres = Mdims%n_in_pres
           if (assemble_collapsed_to_one_phase) then
             !No need to re-scale porous diffusion or CP as they are adjusted by multipliying by the saturation
             ! which sums to one and therefore is the same
              !If collapsed solver then change nphase and n_in_pres
              nphase = Mdims%npres!One temperature per region
              n_in_pres = Mdims%n_in_pres!Need to assemble all the phases
           end if
           !Allocate the RHS
           call allocate(Mmat%CV_RHS,nphase,tracer%mesh,"RHS")
           call allocate(solution,nphase,tracer%mesh,"sol_tracer")!; call zero(solution)
           Loop_NonLinearFlux: DO ITS_FLUX_LIM = 1, NITS_FLUX_LIM

               !Get information for capillary pressure to be use in CV_ASSEMB
                !Over-relaxation options. Unless explicitly decided in diamond this will be set to zero.
               if (thermal) then
                   !Get information for capillary pressure to be use in CV_ASSEMB
                   Phase_with_Ovrel = 1
                   call getOverrelaxation_parameter(state, packed_state, Mdims, ndgln, OvRelax_param, Phase_with_Ovrel&
                                      , totally_min_max = totally_min_max, for_transport = .true.)
                   if (assemble_collapsed_to_one_phase) OvRelax_param = OvRelax_param/ dble(mdims%n_in_pres)
               else
                Phase_with_Ovrel = -1
               end if

               !Solves a PETSC warning saying that we are storing information out of range
               call allocate(Mmat%petsc_ACV,sparsity,[nphase,nphase],"ACV_INTENERGE")
               call zero(Mmat%petsc_ACV); Mmat%CV_RHS%val = 0.0

               !before the sprint in this call the small_acv sparsity was passed as cmc sparsity...
               call CV_ASSEMB( state, packed_state, &
                   n_in_pres, Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd, &
                   tracer, velocity, density, multi_absorp, &
                   DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
                   DEN_ALL, DENOLD_ALL, &
                   cv_disopt, cv_dg_vel_int_opt, DT, cv_theta, cv_beta, &
                   SUF_SIG_DIAGTEN_BC, &
                   DERIV%val(1,:,:), P%val, &
                   T_SOURCE, T_ABSORB, VOLFRA_PORE, &
                   GETCV_DISC, GETCT, &
                   IGOT_T2_loc,IGOT_THETA_FLUX ,GET_THETA_FLUX, USE_THETA_FLUX, &
                   THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
                   MeanPoreCV%val, &
                   mass_Mn_pres, THERMAL, &
                   .false.,  mass_Mn_pres, &
                   mass_ele_transp, &
                   TDIFFUSION = TDIFFUSION,&
                   saturation=saturation, Permeability_tensor_field = perm,&
                   eles_with_pipe =eles_with_pipe, pipes_aux = pipes_aux,&
                   porous_heat_coef = porous_heat_coef,porous_heat_coef_old = porous_heat_coef_old,solving_compositional = lcomp > 0, &
                   VAD_parameter = OvRelax_param, Phase_with_Pc = Phase_with_Ovrel, &
                   assemble_collapsed_to_one_phase = assemble_collapsed_to_one_phase)
                   ! vtracer=as_vector(tracer,dim=2)
                   ! call zero(vtracer)
                   call petsc_solve(solution,Mmat%petsc_ACV,Mmat%CV_RHS,trim(solver_option_path), iterations_taken = its_taken)

               !Copy solution back to tracer(not ideal...)
               do ipres =1, mdims%npres
                 do iphase = 1 , n_in_pres
                   auxI = IPHASE + (ipres-1)*n_in_pres
                   if (assemble_collapsed_to_one_phase) auxI = ipres
                   tracer%val(1,iphase+(ipres-1)*Mdims%n_in_pres,:) = solution%val(auxI,:)
                end do
               end do

                !Control how it is converging and decide
               call force_min_max_principle(Mdims, 2, tracer, nonlinear_iteration, totally_min_max)!Apply if required the min max principle

               !Just after the solvers
               call deallocate(Mmat%petsc_ACV)
              !Update halo communications
              call halo_update(tracer)

               repeat_assemb_solve = (its_taken == 0)!PETSc may fail for a bug then we want to repeat the cycle
               call allor(repeat_assemb_solve)
               !Checking solver not fully implemented
               if (repeat_assemb_solve ) then
                   solver_not_converged = .true.
                   tracer%val(1,:,:) = temp_bak!recover backup
                   cycle!repeat
               else
                   solver_not_converged = its_taken >= max_allowed_its!If failed because of too many iterations we need to continue with the non-linear loop!
                   call allor(solver_not_converged)
                   exit!good to go!
               end if

           END DO Loop_NonLinearFlux


           call deallocate(Mmat%CV_RHS); nullify(Mmat%CV_RHS%val)
           if (allocated(porous_heat_coef)) deallocate(porous_heat_coef)
           if (allocated(porous_heat_coef_old)) deallocate(porous_heat_coef_old)

           if (allocated(den_all)) deallocate(den_all)
           if (allocated(denold_all)) deallocate(denold_all)
           if (allocated(T_SOURCE)) deallocate(T_SOURCE)
           call deallocate(solution); nullify(solution%val)


           ewrite(3,*) 'Leaving INTENERGE_ASSEM_SOLVE'

      contains

      !>@brief: Checks convergence on the temperature field Calculation of the averaged heat capacity and density
      !> average = porosity * Cp_f*rho_f + (1-porosity) * CP_p*rho_p
      !> Since porous promerties is defined element-wise and fluid properties CV-wise we perform an average
      !> as it is stored cv-wise
      subroutine effective_Cp_density(porous_heat_coef, porous_heat_coef_old)
          implicit none
        REAL, DIMENSION( : ), intent(inout) :: porous_heat_coef, porous_heat_coef_old
        !Local variables
        type( scalar_field ), pointer :: porosity, density_porous, Cp_porous, density_porous_old
        integer :: ele, cv_inod, iloc, p_den, h_cap, ele_nod
        real, dimension(Mdims%cv_nonods) :: cv_counter
        real :: auxR

        density_porous => extract_scalar_field( state(1), "porous_density" )
        density_porous_old => extract_scalar_field( state(1), "porous_density_old" )
        Cp_porous => extract_scalar_field( state(1), "porous_heat_capacity" )
        porosity=>extract_scalar_field(state(1),"Porosity")
        porous_heat_coef = 0.
        porous_heat_coef_old = 0.
        cv_counter = 0
        do ele = 1, Mdims%totele
            p_den = min(size(density_porous%val), ele)
            h_cap = min(size(Cp_porous%val), ele)
            ele_nod = min(size(porosity%val), ele)
            do iloc = 1, Mdims%cv_nloc
                cv_inod = ndgln%cv((ele-1)*Mdims%cv_nloc+iloc)
                cv_counter( cv_inod ) = cv_counter( cv_inod ) + 1.0
                porous_heat_coef( cv_inod ) = porous_heat_coef( cv_inod ) + &
                density_porous%val(p_den ) * Cp_porous%val( h_cap )
                porous_heat_coef_old( cv_inod ) = porous_heat_coef_old( cv_inod ) + &
                    density_porous_old%val(p_den ) * Cp_porous%val( h_cap )
            end do
        end do
        !Since nodes are visited more than once, this performs a simple average
        !This is the order it has to be done
        auxR = 1.0
        !If thermal equilibrium for porous media then porous media is added more than once, so needs to be adjusted
        if (assemble_collapsed_to_one_phase) auxR = real(Mdims%n_in_pres)
        porous_heat_coef = porous_heat_coef/(cv_counter*auxR)!<= includes an average of porous and fluid properties
        porous_heat_coef_old = porous_heat_coef_old/(cv_counter*auxR)!<= includes an average of porous and fluid properties5
      end subroutine effective_Cp_density

  END SUBROUTINE INTENERGE_ASSEM_SOLVE

  !>@brief: To help the stability of the system,if there are no sources/sinks it is known that
  !> the temperature must fulfill the min max principle, therefore here values outside this rank are capped.
  subroutine force_min_max_principle(Mdims, entrance, tracer, nonlinear_iteration, totally_min_max)
    type(multi_dimensions), intent(in) :: Mdims
    type(tensor_field), intent(inout) :: tracer
    integer, intent(in) :: nonlinear_iteration
    integer, intent(in) :: entrance
    real, dimension(2), intent(inout) :: totally_min_max
    !Local variables
    integer, allocatable, dimension( :,:,:) :: WIC_T_BC_ALL
    logical :: apply_minmax_principle
    type(tensor_field) :: tracer_BCs
    real, parameter :: tol = 1e-30
    logical :: has_imposed_min_limit, has_imposed_max_limit, has_auto_min_limit, has_auto_max_limit
    logical, save :: WarningMsgShown = .false.

    !Check whether to apply the minmax principle
    apply_minmax_principle = have_option_for_any_phase("scalar_field::"//trim(tracer%name(7:))//"/prognostic/Impose_min_max", Mdims%ndim)
    if (apply_minmax_principle) then
      select case (entrance)
      case (1)
        !Get variable for global convergence method
        totally_min_max = (/-1d30,1d30/)
        has_imposed_min_limit = have_option("/material_phase[0]/scalar_field::"//trim(tracer%name(7:))//"/prognostic/Impose_min_max/min_limit")
        has_imposed_max_limit = have_option("/material_phase[0]/scalar_field::"//trim(tracer%name(7:))//"/prognostic/Impose_min_max/max_limit")
        if (has_imposed_min_limit) call get_option("/material_phase[0]/scalar_field::"//trim(tracer%name(7:))//"/prognostic/Impose_min_max/min_limit", totally_min_max(1))
        if (has_imposed_max_limit) call get_option("/material_phase[0]/scalar_field::"//trim(tracer%name(7:))//"/prognostic/Impose_min_max/max_limit", totally_min_max(2))
        has_auto_min_limit = have_option("/material_phase[0]/scalar_field::"//trim(tracer%name(7:))//"/prognostic/Impose_min_max/automatic_min_limit")
        has_auto_max_limit = have_option("/material_phase[0]/scalar_field::"//trim(tracer%name(7:))//"/prognostic/Impose_min_max/automatic_max_limit")
        if (.not. WarningMsgShown) then
            if (has_auto_min_limit .or. has_auto_max_limit ) then
                if (getprocno() == 1) then
                    ewrite(0,*) 'MESSAGE: Automatic MinMax activated: This should ONLY be used if there are not sources or sinks, and the BCs are Dirichtlet.'
                end if
                WarningMsgShown = .true.
            end if
        end if

        allocate (WIC_T_BC_ALL (1 , Mdims%ndim , surface_element_count(tracer) ))
        call get_entire_boundary_condition(tracer,&
        ['weakdirichlet'], tracer_BCs, WIC_T_BC_ALL)
        !Use boundaries for min/max
        if (has_auto_min_limit) totally_min_max(1)=minval(tracer_BCs%val)
        if (has_auto_max_limit) totally_min_max(2)=maxval(tracer_BCs%val)
        !Check domain
        if (has_auto_min_limit) then
          totally_min_max(1)=min(totally_min_max(1), minval(tracer%val(:,1:Mdims%n_in_pres,:))) !First the reservoir
          if (Mdims%npres > 1) & !Next the wells (this is to avoid the zero values in the well domain outside of the defined regions )
          totally_min_max(1)=min(totally_min_max(1), minval(tracer%val(:, Mdims%n_in_pres+1:,:), MASK = tracer%val(:, Mdims%n_in_pres+1:,:) > tol))
        end if
        if (has_auto_max_limit) totally_min_max(2)=max(totally_min_max(2), maxval(tracer%val))
        !For parallel
        call allmin(totally_min_max(1)); call allmax(totally_min_max(2))
        deallocate(WIC_T_BC_ALL); call deallocate(tracer_BCs)
      case (2)
        tracer%val = max(min(tracer%val,totally_min_max(2)), totally_min_max(1))
      end select
    end if

  end subroutine

    !> @author Chris Pain, Pablo Salinas
    !> @brief Calls to generate and solve the transport equation for n passive tracers defined in diamond as Passive_Tracer_N
    !> Where N is an integer which is continuous starting from 1
    !>@param  state Linked list containing all the fields defined in diamond and considered by Fluidity
    !>@param  packed_state Linked list containing all the fields used by IC-FERST, memory partially shared with state
    !!>@param Mdims Data type storing all the dimensions describing the mesh, fields, nodes, etc
    !>@param  CV_GIdims Gauss integration numbers for CV fields
    !>@param  CV_funs Shape functions for the CV mesh
    !>@param  Mspars Sparsity of the matrices
    !>@param  ndgln Global to local variables
    !>@param  Mdisopt Discretisation options
    !>@param  Mmat Matrices for ICFERST
    !>@param  upwnd Sigmas to compute the fluxes at the interphase for porous media
    !>@param  tracer  Tracer considered for the transport equation
    !>@param  velocity  Velocity of the field
    !>@param  density  Density of the field
    !>@param  multi_absorp  Absoprtion of associated with the transport field
    !>@param  DT Time step size
    !>@param  SUF_SIG_DIAGTEN_BC Like upwnd but for the boundary
    !>@param  VOLFRA_PORE     Porosity field (Mdims%npres,Mdims%totele)
    !>@param  IGOT_T2 !>@param  IGOT_T2. True (1) if solving for a tracer, false otherwise
    !>@param   igot_theta_flux ????
    !>@param  GET_THETA_FLUX, USE_THETA_FLUX ?????
    !>@param  THETA_GDIFF ????
    !>@param  eles_with_pipe Elements that have a pipe
    !>@param  pipes_aux Information required to define wells
    !>@param  mass_ele_transp Mass of the elements
    !>@param  THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J ????
    !>@param  icomp Number of components, if > 0 then compositional is solved
    !>@param  saturation PhaseVolumeFraction field
    !>@param Permeability_tensor_field Permeability field
    !>param nonlinear_iteration Current non-linear iteration
  SUBROUTINE Tracer_Assemble_Solve( Tracer_name, state, packed_state, &
       Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd,&
       tracer, velocity, density, multi_absorp, DT, &
       SUF_SIG_DIAGTEN_BC,  VOLFRA_PORE, &
       IGOT_T2, igot_theta_flux,GET_THETA_FLUX, USE_THETA_FLUX,  &
       THETA_GDIFF, eles_with_pipe, pipes_aux, &
       mass_ele_transp, &
       THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, &
       icomp, saturation, Permeability_tensor_field, nonlinear_iteration )

           implicit none
           character(len=*), intent(in) :: Tracer_name
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
           type(multi_absorption), intent(inout) :: multi_absorp
           INTEGER, intent( in ) :: IGOT_T2, igot_theta_flux
           LOGICAL, intent( in ) :: GET_THETA_FLUX, USE_THETA_FLUX
           REAL, DIMENSION( :, : ), intent( inout ) :: THETA_GDIFF
           REAL, DIMENSION( :,: ), intent( inout ), optional :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
           REAL, intent( in ) :: DT
           REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
           REAL, DIMENSION( :, : ), intent( in ) :: VOLFRA_PORE
           real, dimension( : ), intent( inout ), optional :: mass_ele_transp
           type(tensor_field), intent(in), optional :: saturation
           type( tensor_field ), optional, pointer, intent(in) :: Permeability_tensor_field
           integer, optional :: icomp, nonlinear_iteration
           type(pipe_coords), dimension(:), intent(in):: eles_with_pipe
           type (multi_pipe_package), intent(in) :: pipes_aux
           ! Local variables
           LOGICAL, PARAMETER :: GETCV_DISC = .TRUE., GETCT= .FALSE.
           integer :: nits_flux_lim, its_flux_lim
           REAL, DIMENSION( :, : ), allocatable :: DIAG_SCALE_PRES
           REAL, DIMENSION( :, :, : ), allocatable :: DIAG_SCALE_PRES_COUP, GAMMA_PRES_ABS, INV_B
           REAL, DIMENSION( Mdims%mat_nonods, Mdims%ndim, Mdims%ndim, Mdims%nphase ) :: TDIFFUSION
           REAL, DIMENSION( :,:,:,: ), allocatable :: CDISPERSION
           REAL, DIMENSION( : ), ALLOCATABLE :: MASS_PIPE, MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE
           real, dimension( size(Mspars%small_acv%col )) ::  mass_mn_pres
           REAL, DIMENSION( : , : ), allocatable :: denold_all, t_source
           REAL, DIMENSION( : , : ), target, allocatable :: den_all
           REAL, DIMENSION( : ), allocatable :: CV_RHS_SUB
           type( tensor_field ), pointer :: P, Q
           INTEGER :: IPHASE, its_taken, ipres, i, stat
           type( tensor_field ), pointer :: den_all2, denold_all2, a, aold, deriv, Component_Absorption
           type( vector_field ), pointer  :: MeanPoreCV, python_vfield
           integer :: lcomp, Field_selector, IGOT_T2_loc, python_stat
           type(vector_field)  :: residual
           type(csr_sparsity), pointer :: sparsity
           real, dimension(:,:,:), allocatable :: Velocity_Absorption
           real, dimension(:,:,:), pointer :: T_AbsorB=>null()
           integer :: ncomp_diff_coef, comp_diffusion_opt
           real, dimension(:,:,:), allocatable :: Component_Diffusion_Operator_Coefficient
           type( tensor_field ), pointer :: perm, python_tfield
           integer :: cv_disopt, cv_dg_vel_int_opt
           real :: cv_theta, cv_beta
           type( scalar_field ), pointer :: sfield, porous_field, solid_concentration
           character(len=option_path_len) :: solver_option_path = "/solver_options/Linear_solver"
           !Variables to stabilize the non-linear iteration solver
           real, dimension(2) :: totally_min_max
           logical :: impose_min_max
           real :: aux
           real, save :: inf_tolerance = -1
           !Variables to control the PETCs solver
           integer, save :: max_allowed_its = -1
           !Variables for vanishing diffusion
           real, dimension(Mdims%cv_nonods) :: OvRelax_param
           integer :: Phase_with_Ovrel
           !temperature backup for the petsc bug
           logical :: repeat_assemb_solve
           !Parameters for stabilisation and compact solving, i.e. solving only concentration for some phases
           real, parameter :: min_val = 0.
           integer :: nconc !> Number of phases with tracer, this works if the phases with concentration start from the first one and are consecutive
           integer :: nconc_in_pres
           type(vector_field) :: solution

           !Initialise with an out of range value to be able to check it hasn't been
           totally_min_max = 1e30
           !Retrieve the number of phases that have this tracer, and then if they are concecutive and start from the first one
           nconc = option_count("/material_phase/scalar_field::"//trim(Tracer_name))
           nconc_in_pres = nconc
           if (Mdims%npres > 1) nconc_in_pres = max(nconc_in_pres / 2, 1)
           do iphase = 1, nconc_in_pres
             if (.not. have_option( '/material_phase['// int2str( iphase -1 ) //']/scalar_field::'//trim(Tracer_name))) then
               FLAbort('Concentration must either be defined in all the phases or to start from the first one and consecutively from that one.')
             end if
           end do

           if (present(Permeability_tensor_field)) then
              perm => Permeability_tensor_field
           else
              perm=>extract_tensor_field(packed_state,"Permeability")
           end if

           lcomp = 0
           if ( present( icomp ) ) lcomp = icomp

           call allocate(Mmat%CV_RHS,nconc,tracer%mesh,"RHS")
           call allocate(solution,nconc,tracer%mesh,"sol_tracer")!; call zero(solution)
           sparsity=>extract_csr_sparsity(packed_state,"ACVSparsity")
           allocate(den_all(Mdims%nphase,Mdims%cv_nonods),denold_all(Mdims%nphase,Mdims%cv_nonods))

           allocate( T_SOURCE( Mdims%nphase, Mdims%cv_nonods ) ) ; T_SOURCE=0.0!SPRINT_TO_DO TURN THESE T_SOURCE INTO POINTERS OR DIRECTLY REMOVE THEM
           IGOT_T2_loc = 0

           p => extract_tensor_field( packed_state, "PackedFEPressure" )

          if (has_boussinesq_aprox) then
            !We do not consider variations of density in transport
               den_all = 1
               denold_all =1
          else
             den_all2 => extract_tensor_field( packed_state, "PackedDensity" )
             denold_all2 => extract_tensor_field( packed_state, "PackedOldDensity" )
             den_all    = den_all2 % val ( 1, :, : )
             denold_all = denold_all2 % val ( 1, :, : )
           endif

           IGOT_T2_loc = 1

           call get_option( '/material_phase[0]/scalar_field::Concentration/prognostic/temporal_discretisation/' // &
               'control_volumes/number_advection_iterations', nits_flux_lim, default = 3 )

          !Retrieve source term; sprint_to_do something equivalent should be done for absoprtion
           do iphase = 1, Mdims%nphase !IF THIS WORKS DO THE SAME FOR THE OTHER SCALAR FIELDS
             sfield => extract_scalar_field( state(iphase), trim(Tracer_name)//"Source", stat )
             if (stat == 0) call assign_val(T_source( iphase, : ),sfield%val)
           end do

           !sprint to do, just pass down the other values...
           cv_disopt = Mdisopt%t_disopt; cv_dg_vel_int_opt = Mdisopt%t_dg_vel_int_opt
           cv_theta = Mdisopt%t_theta; cv_beta = Mdisopt%t_beta
           !For passive tracers use low order, the idea is to re-use the matrix when possible
           if (is_PassiveTracer_field(Tracer_name)) then
            cv_disopt = 0!Force upwind
            cv_theta = 1!For implicit euler
           end if

           deriv => extract_tensor_field( packed_state, "PackedDRhoDPressure" )
           TDIFFUSION=0.0;
           call calculate_diffusivity( state, packed_state, Mdims, ndgln, TDIFFUSION, TracerName= trim(Tracer_name))

           !Calculates solute dispersion with specific longitudinal and transverse dispersivity
           if (have_option("/porous_media/Dispersion/scalar_field::Longitudinal_Dispersivity")) then
            allocate(CDISPERSION(Mdims%mat_nonods, Mdims%ndim, Mdims%ndim, Mdims%nphase)); CDISPERSION = 0.
            call calculate_solute_dispersity( state, packed_state, Mdims, ndgln, den_all, CDISPERSION)
            TDIFFUSION = TDIFFUSION + CDISPERSION
            deallocate(CDISPERSION)
           end if

           MeanPoreCV=>extract_vector_field(packed_state,"MeanPoreCV")

           solver_option_path = "/solver_options/Linear_solver"
           if (have_option('/solver_options/Linear_solver/Custom_solver_configuration/field::Passive_Tracers')) then
             solver_option_path = '/solver_options/Linear_solver/Custom_solver_configuration/field::Passive_Tracers'
           end if
           if(max_allowed_its < 0)  then
               call get_option( trim(solver_option_path)//"max_iterations",&
                max_allowed_its, default = 500)
           end if

           !Start with the process to apply the min max principle
           call force_min_max_principle(Mdims, 1, tracer, nonlinear_iteration, totally_min_max)

           Loop_NonLinearFlux: DO ITS_FLUX_LIM = 1, NITS_FLUX_LIM

                !Over-relaxation options. Unless explicitly decided in diamond this will be set to zero.
                !Get information for capillary pressure to be use in CV_ASSEMB
                Phase_with_Ovrel = 1
                call getOverrelaxation_parameter(state, packed_state, Mdims, ndgln, OvRelax_param, Phase_with_Ovrel,&
                                                        totally_min_max = totally_min_max, for_transport = .true.)

               !Solves a PETSC warning saying that we are storing information out of range
               call allocate(Mmat%petsc_ACV,sparsity,[nconc,nconc],"ACV_Passive_Tracer")
               call zero(Mmat%petsc_ACV); Mmat%CV_RHS%val = 0.0

               !before the sprint in this call the small_acv sparsity was passed as cmc sparsity...
               call CV_ASSEMB( state, packed_state, &
                   nconc_in_pres, Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd, &
                   tracer, velocity, density, multi_absorp, &
                   DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
                   DEN_ALL, DENOLD_ALL, &
                   cv_disopt, cv_dg_vel_int_opt, DT, cv_theta, cv_beta, &
                   SUF_SIG_DIAGTEN_BC, &
                   DERIV%val(1,:,:), P%val, &
                   T_SOURCE, T_ABSORB, VOLFRA_PORE, &
                   GETCV_DISC, GETCT, &
                   IGOT_T2_loc,IGOT_THETA_FLUX ,GET_THETA_FLUX, USE_THETA_FLUX, &
                   THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
                   MeanPoreCV%val, &
                   mass_Mn_pres, .false., &
                   .true.,  mass_Mn_pres, &
                   mass_ele_transp, &
                   TDIFFUSION = TDIFFUSION,&
                   saturation=saturation, Permeability_tensor_field = perm,&
                   eles_with_pipe =eles_with_pipe, pipes_aux = pipes_aux,&
                   solving_compositional = lcomp > 0, &
                   VAD_parameter = OvRelax_param, Phase_with_Pc = Phase_with_Ovrel)
                   ! call zero_non_owned(Mmat%CV_RHS)
                   call zero(solution)
                   call petsc_solve(solution,Mmat%petsc_ACV,Mmat%CV_RHS,trim(solver_option_path), iterations_taken = its_taken)

                   !Copy solution back to tracer(not ideal...)
                   do ipres =1, mdims%npres
                     do iphase = 1 , nconc_in_pres
                      tracer%val(1,iphase+(ipres-1)*Mdims%n_in_pres,:) = solution%val(iphase+(ipres-1)*nconc_in_pres,:)
                    end do
                   end do
                   !Apply if required the min max principle
                   call force_min_max_principle(Mdims, 2, tracer, nonlinear_iteration, totally_min_max)

                   !Just after the solvers
                   call deallocate(Mmat%petsc_ACV)!<=There is a bug, if calling Fluidity to deallocate the memory of the PETSC matrix
                   !Update halo communications
                   call halo_update(tracer)
                   !Checking solver not fully implemented
                   solver_not_converged = its_taken >= max_allowed_its!If failed because of too many iterations we need to continue with the non-linear loop!
                   call allor(solver_not_converged)
                   exit!good to go!

           END DO Loop_NonLinearFlux

           call deallocate(Mmat%CV_RHS); nullify(Mmat%CV_RHS%val)
           if (allocated(den_all)) deallocate(den_all)
           if (allocated(denold_all)) deallocate(denold_all)
           if (allocated(T_SOURCE)) deallocate(T_SOURCE)
           call deallocate(solution); nullify(solution%val)
           ewrite(3,*) 'Leaving' //trim(Tracer_name)//'_assem_solve'

  END SUBROUTINE Tracer_Assemble_Solve


    !---------------------------------------------------------------------------
    !> @author Chris Pain, Pablo Salinas
    !> @brief Calls to generate the transport equation for the saturation. Embeded an FPI with backtracking method is uncluded
    !>@param  state Linked list containing all the fields defined in diamond and considered by Fluidity
    !>@param  packed_state Linked list containing all the fields used by IC-FERST, memory partially shared with state
    !>@param  multicomponent_state Linked list containing all the fields used by compositional, just in case we are solving for components
    !!>@param Mdims Data type storing all the dimensions describing the mesh, fields, nodes, etc
    !>@param  CV_GIdims Gauss integration numbers for CV fields
    !>@param  CV_funs Shape functions for the CV mesh
    !>@param  Mspars Sparsity of the matrices
    !>@param  ndgln Global to local variables
    !>@param  Mdisopt Discretisation options
    !>@param  Mmat Matrices for ICFERST
    !>@param  multi_absorp  Absoprtion of associated with the transport field 
    !>@param  upwnd Sigmas to compute the fluxes at the interphase for porous media
    !>@param  eles_with_pipe Elements that have a pipe
    !>@param  pipes_aux Information required to define wells
    !>@param  DT Time step size
    !>@param  SUF_SIG_DIAGTEN_BC Like upwnd but for the boundary
    !>@param V_SOURCE Source term
    !>@param  VOLFRA_PORE     Porosity field (Mdims%npres,Mdims%totele)
    !>@param   igot_theta_flux ????
    !>@param  mass_ele_transp Mass of the elements
    !>@param nonlinear_iteration Current non-linear iteration
    !>@param time_step current time-step
    !>@param SFPI_its Number of saturation fixed point iterations taken 
    !>@param Courant_number Global courant number and shock front courant number
    !>@param  THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J ????
    subroutine VolumeFraction_Assemble_Solve( state,packed_state, multicomponent_state, &
         Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, multi_absorp, upwnd, &
         eles_with_pipe, pipes_aux, DT, SUF_SIG_DIAGTEN_BC, &
         V_SOURCE, VOLFRA_PORE, igot_theta_flux, mass_ele_transp,&
         nonlinear_iteration, time_step, SFPI_taken, SFPI_its, Courant_number,&
         THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J)
             implicit none
             type( state_type ), dimension( : ), intent( inout ) :: state, multicomponent_state
             type( state_type ) :: packed_state
             type(multi_dimensions), intent(in) :: Mdims
             type(multi_GI_dimensions), intent(in) :: CV_GIdims
             type(multi_shape_funs), intent(inout) :: CV_funs
             type(multi_sparsities), intent(in) :: Mspars
             type(multi_ndgln), intent(in) :: ndgln
             type (multi_discretization_opts) :: Mdisopt
             type (multi_matrices), intent(inout) :: Mmat
             type(multi_absorption), intent(inout) :: multi_absorp
             type (porous_adv_coefs), intent(inout) :: upwnd
             type(pipe_coords), dimension(:), intent(in):: eles_with_pipe
             type (multi_pipe_package), intent(in) :: pipes_aux
             INTEGER, intent( in ) :: igot_theta_flux
             REAL, intent( in ) :: DT
             REAL, DIMENSION( :, : ), intent( inout ) :: SUF_SIG_DIAGTEN_BC
             REAL, DIMENSION( :, : ), intent( in ) :: V_SOURCE
             !REAL, DIMENSION( :, :, : ), intent( in ) :: V_ABSORB
             REAL, DIMENSION( :, : ), intent( in ) :: VOLFRA_PORE
             real, dimension( : ), intent( inout ) :: mass_ele_transp
             integer, intent(in) :: nonlinear_iteration
             integer, intent(in) :: time_step
             integer, intent(inout) :: SFPI_taken
             integer, intent(inout) :: SFPI_its
             real, dimension(:), intent(inout) :: Courant_number
             REAL, DIMENSION( :, :), intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J

             ! Local Variables
             LOGICAL, PARAMETER :: THERMAL= .false.
             integer :: igot_t2, nphase, n_in_pres
             REAL, DIMENSION( : ), allocatable :: mass_mn_pres
             REAL, DIMENSION( :, : ), allocatable :: DIAG_SCALE_PRES
             REAL, DIMENSION( :, :, : ), allocatable :: DIAG_SCALE_PRES_COUP, INV_B
             REAL, DIMENSION( :, : ), allocatable :: THETA_GDIFF
             REAL, DIMENSION( :, : ), pointer :: DEN_ALL, DENOLD_ALL
             REAL, DIMENSION( :, : ), allocatable :: T2, T2OLD
             REAL, DIMENSION( :, : ), allocatable :: MEAN_PORE_CV
             LOGICAL :: GET_THETA_FLUX
             INTEGER :: STAT, IPHASE, JPHASE, IPHASE_REAL, JPHASE_REAL, IPRES, JPRES, cv_nodi
             LOGICAL, PARAMETER :: GETCV_DISC = .TRUE., GETCT= .FALSE.
             type( tensor_field ), pointer :: den_all2, denold_all2
             character(len=option_path_len) :: solver_option_path = "/solver_options/Linear_solver"
             !Working pointers
             real, dimension(:,:,:), pointer :: p, V_ABSORB => null() ! this is PhaseVolumeFraction_AbsorptionTerm
             real, dimension(:, :), pointer :: satura
             type(tensor_field), pointer :: velocity, density, deriv, sat_field
             type(scalar_field), pointer :: gamma
             type(vector_field) :: solution
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
             integer :: Max_sat_its, total_cv_nodes
             real, dimension(:,:), allocatable :: sat_bak, backtrack_sat
             real :: Previous_convergence, updating, new_backtrack_par, aux, resold, first_res
             real, save :: res = -1
             logical :: satisfactory_convergence
             integer :: its, useful_sats
             logical, save :: Solve_all_phases
             logical, save :: ML_method_activated
             !Variables to control the PETCs solver
             integer, save :: max_allowed_its = -1
             integer :: its_taken
             logical, save :: written_file = .false.
             !We check this with the global number of phases per domain
             if ( Mdims%n_in_pres == 1) return!<== No need to solve the transport of phases if there is only one phase!

             solver_option_path = "/solver_options/Linear_solver"
             if (have_option('/solver_options/Linear_solver/Custom_solver_configuration/field::PhaseVolumeFraction')) then
               solver_option_path = '/solver_options/Linear_solver/Custom_solver_configuration/field::PhaseVolumeFraction'
             end if
             if(max_allowed_its < 0)  then
                 call get_option( trim(solver_option_path)//"max_iterations",&
                  max_allowed_its, default = 500)
                  Solve_all_phases = .not. have_option("/numerical_methods/solve_nphases_minus_one")
                  ML_method_activated =  have_option("/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/ML_model_path")
             end if

             if (Solve_all_phases) then!For the number_of_normal_FPI FPIs we consider all the phases normally
               nphase = Mdims%nphase
               n_in_pres = Mdims%n_in_pres
             else !For subsequent SFPI iterations just perform nphases - 1
               !Define local nphase and n_in_pres
               nphase = (Mdims%n_in_pres - 1 ) * Mdims%npres
               n_in_pres = nphase/ Mdims%npres
             end if
             !Extract variables from packed_state
             call get_var_from_packed_state(packed_state,FEPressure = P)
             sat_field => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )

             Satura =>  sat_field%val(1,:,:)
             !Get information for capillary pressure to be use in CV_ASSEMB
             call getOverrelaxation_parameter(state, packed_state, Mdims, ndgln, OvRelax_param, Phase_with_Pc)

             !Get variable for global convergence method
             if (.not. have_option( '/solver_options/Non_Linear_Solver/Fixed_Point_Iteration')) then
                 backtrack_par_factor = 1.1
             else !Get value with the default value of 1.
                 call get_option( '/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Backtracking_factor',&
                     backtrack_par_factor, default = 1.0)
             end if

             !For backtrack_par_factor == -10 we will set backtrack_par_factor based on the shock front Courant number
             Auto_max_backtrack = (backtrack_par_factor == -10)
             !Retrieve number of saturation fixed point iterations from diamond, by default 3 if Courant_number<=1, 9 otherwise
             if (Courant_number(1) <= 1) then
                 call get_option( "/numerical_methods/max_sat_its", max_sat_its, default = 3)
             else
                 call get_option( "/numerical_methods/max_sat_its", max_sat_its, default = 9)
             end if
             ewrite(3,*) 'In VOLFRA_ASSEM_SOLVE'
             GET_THETA_FLUX = .FALSE.
             !####Create dummy variables required for_cv_assemb with no memory usage ####
             IGOT_T2 = 0
             ALLOCATE( THETA_GDIFF( 0, 0 )); ALLOCATE( Mmat%CT( 0,0,0 ) )
             ALLOCATE( DIAG_SCALE_PRES( 0,0 ) ); ALLOCATE( DIAG_SCALE_PRES_COUP( 0,0,0 ), INV_B( 0,0,0 ) )
             !#### Create dummy variables required for_cv_assemb with no memory usage ####
             deriv => extract_tensor_field( packed_state, "PackedDRhoDPressure" )
             ALLOCATE( MEAN_PORE_CV( Mdims%npres, Mdims%cv_nonods ) )
             ALLOCATE( mass_mn_pres(size(Mspars%small_acv%col)) ) ; mass_mn_pres = 0.


             Mdisopt%v_beta = 1.0
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
             total_cv_nodes = Mdims%cv_nonods
             call allsum(total_cv_nodes)!For parallel consistency when normalising the residual
             !Now total_cv_nodes includes halos, but because it is a ratio it should be fine
             if (resold < 0 ) res = huge(res)!<=initialize res once


            nullify(DEN_ALL); nullify(DENOLD_ALL)
            allocate (sat_bak(Mdims%nphase, Mdims%cv_nonods), backtrack_sat(Mdims%nphase, Mdims%cv_nonods))
            !Allocate residual, to compute the residual
            if (backtrack_par_factor < 1.01) call allocate(residual,nphase,sat_field%mesh,"residual")
            call allocate(Mmat%CV_RHS,nphase,sat_field%mesh,"RHS")
            call allocate(solution,nphase,sat_field%mesh,"Saturation")!; call zero(solution)

            IF ( IGOT_THETA_FLUX == 1 ) THEN ! We have already put density in theta...
              ! use DEN=1 because the density is already in the theta variables
              ALLOCATE( DEN_ALL( nphase, Mdims%cv_nonods )); DEN_ALL = 1.
              ALLOCATE( DENOLD_ALL( nphase, Mdims%cv_nonods )); DENOLD_ALL = 1.
            ELSE
              DEN_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedDensity" )
              DENOLD_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldDensity" )
              DEN_ALL => DEN_ALL2%VAL( 1, :, : ) ; DENOLD_ALL => DENOLD_ALL2%VAL( 1, :, : )
            END IF

             Loop_NonLinearFlux: do while (.not. satisfactory_convergence)

               !To avoid a petsc warning error we need to re-allocate the matrix always
               call allocate_global_multiphase_petsc_csr(Mmat%petsc_ACV,sparsity,sat_field, nphase)
                !Update solution field to calculate the residual
                  do ipres =1, mdims%npres
                    do iphase = 1 , n_in_pres
                     solution%val(iphase+(ipres-1)*n_in_pres,:) = sat_field%val(1,iphase+(ipres-1)*n_in_pres,:)
                   end do
                 end do

                 !Assemble the matrix and the RHS
                 !before the sprint in this call the small_acv sparsity was passed as cmc sparsity...
                 call CV_ASSEMB( state, packed_state, &
                     n_in_pres, Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,upwnd,&
                     sat_field, velocity, density, multi_absorp, &
                     DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
                     DEN_ALL, DENOLD_ALL, &
                     Mdisopt%v_disopt, Mdisopt%v_dg_vel_int_opt, DT, Mdisopt%v_theta, Mdisopt%v_beta, &
                     SUF_SIG_DIAGTEN_BC, &
                     DERIV%val(1,:,:), P, &
                     V_SOURCE, V_ABSORB, VOLFRA_PORE, &
                     GETCV_DISC, GETCT, &
                     IGOT_T2, igot_theta_flux, GET_THETA_FLUX, Mdisopt%volfra_get_theta_flux, &
                     THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
                     MEAN_PORE_CV, &
                     mass_Mn_pres, THERMAL, &
                     .false.,  mass_Mn_pres, &
                     mass_ele_transp, &          !Capillary variables
                     VAD_parameter = OvRelax_param, Phase_with_Pc = Phase_with_Pc,&
                     Courant_number = Courant_number, eles_with_pipe = eles_with_pipe, pipes_aux = pipes_aux,&
                     nonlinear_iteration = nonlinear_iteration)

                 !Make the inf norm of the Courant number across cpus
                 !Normally computed when dealing with the continuity equation but
                 !if solving for saturation it is useful to have up to date information
                 if (IsParallel()) then
                    call allmax(Courant_number(1)); call allmax(Courant_number(2))
                 end if

                 !Time to solve the system
                 !If using FPI with backtracking
                 if (backtrack_par_factor < 1.01) then
                     !Backup of the saturation field, to adjust the solution
                     sat_bak = satura

                     !If using ADAPTIVE FPI with backtracking
                     if (backtrack_par_factor < 0) then
#ifdef USING_XGBOOST
                        if (.not. ML_method_activated) then
                            if (Auto_max_backtrack) then!The maximum backtracking factor depends on the shock-front Courant number
                                call auto_backtracking(Mdims, backtrack_par_factor, courant_number, first_time_step, nonlinear_iteration)
                            end if
                        end if
#else
                        if (Auto_max_backtrack) then!The maximum backtracking factor depends on the shock-front Courant number
                            call auto_backtracking(Mdims, backtrack_par_factor, courant_number, first_time_step, nonlinear_iteration)
                        end if
#endif
                        !Calculate the actual residual using a previous backtrack_par
                        ! vtracer=as_vector(sat_field,dim=2)
                        call mult(residual, Mmat%petsc_ACV, solution)
                        !Now residual is calculated as Residual = RHS - A*X0
                        residual%val = Mmat%CV_RHS%val - residual%val
                        if (IsParallel()) call halo_update(residual)!better than zero_non_owned, important for parallel
                        resold = res; res = 0
                        do iphase = 1, nphase
                        !L2 norm of the residual; needs to be done in steps to ensure cosistency in parallel
                            aux = dot_product(residual%val(iphase,:),residual%val(iphase,:))
                            call allsum(aux)
                            aux = sqrt(aux)/ dble(total_cv_nodes)
                            if (aux > res) res = aux
                        end do
                        !We use the highest residual across the domain
                        if (IsParallel()) call allmax(res)

                        if (its == 1) then
                            first_res = res !Variable to check total convergence of the SFPI method
#ifdef USING_XGBOOST
                            if (ML_method_activated) then
                                if (Auto_max_backtrack) then!The maximum backtracking factor depends on the shock-front Courant number (auto_backtracking) or a set of dimensionless numbers (AI_backtracking_parameters)
                                    !#=================================================================================================================
                                    !# Vinicius: Added a subroutine for calculating all the dimensioless numbers required fo the ML model
                                    !#=================================================================================================================
                                    ! Calculate Darcy velocity with the most up-to-date information (necessary for AI_backtracking_parameters)
                                    call get_DarcyVelocity( Mdims, ndgln, state, packed_state, upwnd )
                                    ! Generate all the dimensionless numbers
                                    call AI_backtracking_parameters(Mdims, ndgln, packed_state, state, courant_number, backtrack_par_factor, OvRelax_param, res, resold, nonlinear_iteration)
                                    !#=================================================================================================================
                                    !# Vinicius-End: Added a subroutine for calculating all the dimensioless numbers required fo the ML model
                                    !#=================================================================================================================
                                end if
                            end if
#endif
                        end if
                     end if
                 end if

                 call zero(solution)

                 !########Solve the system#############
                 call petsc_solve(solution,Mmat%petsc_ACV,Mmat%CV_RHS,trim(solver_option_path), iterations_taken = its_taken)
                 !To avoid a petsc warning error we need to re-allocate the matrix always
                 call deallocate(Mmat%petsc_ACV)
 ! call MatView(Mmat%petsc_ACV%M,   PETSC_VIEWER_STDOUT_SELF, ipres)
                !Copy solution back to sat_field (not ideal...)
                  do ipres =1, mdims%npres
                    do iphase = 1 , n_in_pres
                     sat_field%val(1,iphase+(ipres-1)*Mdims%n_in_pres,:) = solution%val(iphase+(ipres-1)*n_in_pres,:)
                   end do
                 ! end do
                    if (backtrack_par_factor > 1.01) then
                     !Convert from nphases-1 to nphases (only if we are not using FPI_backtracking)
                     do cv_nodi = 1, Mdims%cv_nonods
                       sat_field%val(1,Mdims%n_in_pres+(ipres-1)*Mdims%n_in_pres,cv_nodi) = 1.0 - &
                            sum(sat_field%val(1,1+(ipres-1)*Mdims%n_in_pres:(Mdims%n_in_pres-1) + (ipres-1)*Mdims%n_in_pres,cv_nodi))
                     end do
                   end if
                 end do

                 !Set to zero the fields
                 call zero(Mmat%CV_RHS)
                 !Correct the solution obtained to make sure we are on track towards the final solution
                 ! if (Mdims%ncomp > 0) call update_components()
                 !Correct the solution obtained to make sure we are on track towards the final solution
                 if (backtrack_par_factor < 1.01) then
                     !If convergence is not good, then we calculate a new saturation using backtracking
                     if (.not. satisfactory_convergence) then

                         !Calculate a backtrack_par parameter and update saturation with that parameter, ensuring convergence
                         call FPI_backtracking(nphase, Mdims, ndgln, state,packed_state, sat_bak(1:nphase, :), backtrack_sat(1:nphase, :), backtrack_par_factor,&
                             Previous_convergence, satisfactory_convergence, new_backtrack_par, Max_sat_its, its, nonlinear_iteration,&
                             useful_sats,res, res/resold, first_res) !halos are updated within this subroutine

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
                             call Calculate_PorousMedia_AbsorptionTerms( nphase, state, packed_state, multi_absorp%PorousMedia, Mdims, &
                                   CV_funs, CV_GIdims, Mspars, ndgln, upwnd, suf_sig_diagten_bc )

                             !Also recalculate the Over-relaxation parameter
                            call getOverrelaxation_parameter(state, packed_state, Mdims, ndgln, OvRelax_param, Phase_with_Pc)

                         else
                             exit Loop_NonLinearFlux
                         end if

                     end if
                 else !Just one iteration
                     if (IsParallel()) call halo_update(sat_field)
                     exit Loop_NonLinearFlux
                 end if
                 its = its + 1
                 useful_sats = useful_sats + 1

             END DO Loop_NonLinearFlux

             !#=================================================================================================================
             !# Vinicius: Write the number of non-linear iterations into file
             !#=================================================================================================================
             !inquire(file="Inner_non_linear_iterations.csv", exist=file_exist)
            !  if (.not. written_file) then
            !      open(75, file="non_linear_iterations.csv", status="replace")
            !      write(75, '(8(A,",",X))') "time_step", "outer_nonlinear_iteration", "Inner_non_linear_iterations"
            !      close(75)
            !      written_file = .true.
            !  end if
            !  ! Write values
            !  open(75, file="non_linear_iterations.csv", status="unknown", position="append")
            !  write(75, '(8(I8,",",X))') time_step, nonlinear_iteration, its
            !  close(75)
            !  SFPI_its = its
             !#=================================================================================================================
             !# Vinicius-End: Write the number of non-linear iterations into file
             !#=================================================================================================================

             !Store the number of Saturation Fixed_Point iterations
             SFPI_taken = SFPI_taken + its
             !Store the final accumulated backtrack_par_factor to properly calculate the convergence functional
             if (backtrack_par_factor < 1.01) then
                 !Final effective backtrack_par to calculate properly the non linear convergence is:
                 backtrack_or_convergence = updating
             else
                 backtrack_or_convergence = 1
             end if

             !If the final saturation solve of the final non-linear FPI fails, then we ensure that the result is not accepted
             !if using adaptive time-stepping of some sort, the loop will be repeated. In all the cases a Warning message will show up
             if (its_taken >= max_allowed_its  .or. its_taken == 0 ) solver_not_converged = .true.

             !Make sure the parameter is consistent between cpus
             if (IsParallel()) call allmin(backtrack_or_convergence)

             !#### Deallocate dummy variables required for_cv_assemb with no memory usage ####
             deallocate( THETA_GDIFF, Mmat%CT, DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B )
             !#### Deallocate dummy variables required for_cv_assemb with no memory usage ####

             DEALLOCATE( mass_mn_pres, MEAN_PORE_CV)
             DEALLOCATE( sat_bak, backtrack_sat )
             call deallocate(Mmat%CV_RHS); nullify(Mmat%CV_RHS%val)
             if (backtrack_par_factor < 1.01) call deallocate(residual)
             !Deallocate pointers only if not pointing to something in packed state
             if (IGOT_THETA_FLUX == 1 ) then
                 deallocate(DEN_ALL, DENOLD_ALL)
             end if
             nullify(DEN_ALL); nullify(DENOLD_ALL)
             call deallocate(solution)
             ewrite(3,*) 'Leaving VOLFRA_ASSEM_SOLVE'

         contains

        !!!>@brief: This internal subroutine deals with the components within the Saturation Fixed Point Iteration
        !> WARNING: Still work in progress
        subroutine update_components()
          implicit none
          real, dimension(Mdims%nphase, Mdims%cv_nonods) :: comp_theta_gdiff

          !First, impose physical constrains to the saturation (important to update halos here)
          call Set_Saturation_to_sum_one(mdims, packed_state, state, do_not_update_halos = .false. )
          !Next, update compoents
          !Deallocate memory re-used for the compositional assembly solve; SPRINT_TO_DO: this can be done better!
          call deallocate(Mmat%CV_RHS); nullify(Mmat%CV_RHS%val); call deallocate(Mmat%petsc_ACV)
          call Compositional_Assemble_Solve(state, packed_state, multicomponent_state, &
               Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd,&
               multi_absorp, DT, &
               SUF_SIG_DIAGTEN_BC, &
               Mdisopt%comp_get_theta_flux, Mdisopt%comp_use_theta_flux,  &
               comp_theta_gdiff, eles_with_pipe, pipes_aux, mass_ele_transp, &
               THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J)

           !Re-allocate the fields so the saturation loop has consistency with what it expects. SPRINT_TO_DO: this can be done better!
           call allocate(Mmat%CV_RHS,nphase,sat_field%mesh,"RHS")
           call allocate_global_multiphase_petsc_csr(Mmat%petsc_ACV,sparsity,sat_field, nphase)

           !First, impose physical constrains to the saturation (important to update halos here)
           call Set_Saturation_to_sum_one(mdims, packed_state, state, do_not_update_halos = .false. )

        end subroutine update_components

    end subroutine VolumeFraction_Assemble_Solve



    !>@brief:In this subroutine the components are solved for all the phases.
    !>Systems for each component are assembled and solved by calling INTENERGE_ASSEM_SOLVE
    !>@param  state Linked list containing all the fields defined in diamond and considered by Fluidity
    !>@param  packed_state Linked list containing all the fields used by IC-FERST, memory partially shared with state
    !>@param  multicomponent_state Linked list containing all the fields used by compositional, just in case we are solving for components
    !!>@param Mdims Data type storing all the dimensions describing the mesh, fields, nodes, etc
    !>@param  CV_GIdims Gauss integration numbers for CV fields
    !>@param  CV_funs Shape functions for the CV mesh
    !>@param  Mspars Sparsity of the matrices
    !>@param  ndgln Global to local variables
    !>@param  Mdisopt Discretisation options
    !>@param  Mmat Matrices for ICFERST
    !>@param  upwnd Sigmas to compute the fluxes at the interphase for porous media
    !>@param  multi_absorp  Absoprtion of associated with the transport field
    !>@param  DT Time step size
    !>@param  SUF_SIG_DIAGTEN_BC Like upwnd but for the boundary
    !>@param  GET_THETA_FLUX, USE_THETA_FLUX ?????
    !>@param  THETA_GDIFF ????
    !>@param  eles_with_pipe Elements that have a pipe
    !>@param  pipes_aux Information required to define wells
    !>@param  mass_ele Mass of the elements
    !>@param  sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j ????
    subroutine Compositional_Assemble_Solve(state, packed_state, multicomponent_state, &
         Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd,&
         multi_absorp, DT, &
         SUF_SIG_DIAGTEN_BC, &
         GET_THETA_FLUX, USE_THETA_FLUX,  &
         THETA_GDIFF, eles_with_pipe, pipes_aux, mass_ele, &
         sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j)
         implicit none
         type( state_type ), dimension( : ), intent( inout ) :: state, multicomponent_state
         type( state_type ), intent( inout ) :: packed_state
         type(multi_dimensions), intent(in) :: Mdims
         type(multi_GI_dimensions), intent(in) :: CV_GIdims
         type(multi_shape_funs), intent(inout) :: CV_funs
         type (multi_sparsities), intent(in) :: Mspars
         type(multi_ndgln), intent(in) :: ndgln
         type (multi_discretization_opts) :: Mdisopt
         type (multi_matrices), intent(inout) :: Mmat
         type (porous_adv_coefs), intent(inout) :: upwnd
         type(multi_absorption), intent(inout) :: multi_absorp
         LOGICAL, intent( in ) :: GET_THETA_FLUX, USE_THETA_FLUX
         real, dimension(:), intent(in) :: mass_ele
         REAL, DIMENSION( :, : ), intent( inout ) :: THETA_GDIFF
         REAL, DIMENSION( :,: ), intent( inout ) :: sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j
         REAL, intent( in ) :: DT
         REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
         type(pipe_coords), dimension(:), intent(in):: eles_with_pipe
         type (multi_pipe_package), intent(in) :: pipes_aux
         !Local variables
         real, dimension(Mdims%nphase) :: RSUM
         integer :: icomp, its2, iphase, jphase, cv_nodi, cv_inod, NonLinearIteration_Components, stat, ncv_faces
         !Working pointers
         type(tensor_field), pointer :: tracer_field, velocity_field, multicomp_density_field, saturation_field, old_saturation_field   !, tracer_source
         type(tensor_field), pointer :: PhaseVolumeFractionComponentSource, PhaseDensity, ComponentDensity, OldComponentDensity
         type(tensor_field), pointer :: Component_Absorption, perm_field, ComponentMassFraction, OldComponentMassFraction
         type(vector_field), pointer :: porosity_field, MeanPoreCV
         real, dimension( :, : ), allocatable ::theta_flux, one_m_theta_flux, theta_flux_j, one_m_theta_flux_j

         !Obtain the number of faces in the control volume space
         ncv_faces=CV_count_faces( Mdims, Mdisopt%cv_ele_type, CV_GIDIMS = CV_GIdims)
         allocate(theta_flux( Mdims%nphase, ncv_faces  ), &
            one_m_theta_flux( Mdims%nphase, ncv_faces ), &
            theta_flux_j( Mdims%nphase, ncv_faces ), &
            one_m_theta_flux_j( Mdims%nphase, ncv_faces  ))

          theta_flux = 0.;one_m_theta_flux = 0.; theta_flux_j = 0; one_m_theta_flux_j = 0
         !Quick check to ensure that we need to solve for components
         if( Mdims%ncomp == 0 ) return

         call get_option( '/numerical_methods/Max_compositional_its', NonLinearIteration_Components, default = 1 )
        !Retrieve fields
        perm_field => extract_tensor_field(packed_state,"Permeability")
        velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
        saturation_field=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
        old_saturation_field=>extract_tensor_field(packed_state,"PackedOldPhaseVolumeFraction")
        MeanPoreCV=>extract_vector_field(packed_state,"MeanPoreCV")
        porosity_field=>extract_vector_field(packed_state,"Porosity")
        PhaseVolumeFractionComponentSource => extract_tensor_field(packed_state,"PackedPhaseVolumeFractionComponentSource")
        PhaseVolumeFractionComponentSource%val = 0.0

        PhaseDensity  => extract_tensor_field( packed_state, "PackedDensity" )
        ComponentDensity  => extract_tensor_field( packed_state, "PackedComponentDensity" )
        OldComponentDensity  => extract_tensor_field( packed_state, "PackedOldComponentDensity" )
        ComponentMassFraction  => extract_tensor_field( packed_state, "PackedComponentMassFraction" )
        OldComponentMassFraction  => extract_tensor_field( packed_state, "PackedOldComponentMassFraction" )

        !!$ Starting loop over components
        Loop_Components: do icomp = 1, Mdims%ncomp
            tracer_field=>extract_tensor_field(multicomponent_state(icomp),"PackedComponentMassFraction")
            multicomp_density_field=>extract_tensor_field(multicomponent_state(icomp),"PackedComponentDensity",stat)

            if( have_option( '/material_phase[' // int2str( Mdims%nstate - Mdims%ncomp ) // &
                ']/is_multiphase_component/KComp_Sigmoid' ) .and. Mdims%nphase > 1 ) then

                !!$ Computing the absorption term for the multi-components equation

                Component_Absorption => extract_tensor_field( multicomponent_state(icomp), "ComponentAbsorption")

                call Calculate_ComponentAbsorptionTerm( packed_state, icomp, ndgln%cv, &
                     Mdims, PhaseDensity%val, Porosity_field%val, mass_ele, Component_Absorption%val )

            end if
            !!$ NonLinear iteration for the components advection:
            Loop_NonLinearIteration_Components: do its2 = 1, NonLinearIteration_Components

                Mdisopt%comp_use_theta_flux = .false. ; Mdisopt%comp_get_theta_flux = .true.
                call INTENERGE_ASSEM_SOLVE( state, multicomponent_state(icomp), &
                    Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,upwnd,&
                    tracer_field,velocity_field,multicomp_density_field, multi_absorp, dt, &
                    SUF_SIG_DIAGTEN_BC, Porosity_field%val, &
                    1, 1, &!igot_t2, igot_theta_flux
                    Mdisopt%comp_get_theta_flux, Mdisopt%comp_use_theta_flux, &
                    theta_gdiff, eles_with_pipe, pipes_aux,&
                    thermal = .false.,& ! the false means that we don't add an extra source term
                    theta_flux=theta_flux, one_m_theta_flux=one_m_theta_flux, theta_flux_j=theta_flux_j, one_m_theta_flux_j=one_m_theta_flux_j,&
                    icomp=icomp, saturation=saturation_field, Permeability_tensor_field = perm_field)

                !This is to ensure boundedness of the ComponentMassFraction (OLD METHOD)
                tracer_field%val = min (max( tracer_field%val, 0.0), 1.0)

            end do Loop_NonLinearIteration_Components

            !These variables below are computed but never used...
            sum_theta_flux = sum_theta_flux + theta_flux
            sum_one_m_theta_flux = sum_one_m_theta_flux + one_m_theta_flux
            sum_theta_flux_j = sum_theta_flux_j + theta_flux_j
            sum_one_m_theta_flux_j = sum_one_m_theta_flux_j + one_m_theta_flux_j

            ! We have divided through by density
            PhaseVolumeFractionComponentSource%val(1,:,:) = PhaseVolumeFractionComponentSource%val(1,:,:) + THETA_GDIFF

        end do Loop_Components

        !New method to ensure that components sum to one, similar to the saturation method
        ! call ensure_components_sum_to_one(packed_state)

        !OLD METHOD BELOW (SPRINT_TO_DO REMOVE IT IF THE NEW ONE WORKS BETTER)
        if ( have_option( '/material_phase[' // int2str( Mdims%nstate - Mdims%ncomp ) // &
            ']/is_multiphase_component/Comp_Sum2One/Enforce_Comp_Sum2One' ) ) then
            ! Initially clip and then ensure the components sum to unity so we don't get surprising results...
            ComponentMassFraction % val = min ( max ( ComponentMassFraction % val, 0.0), 1.0)

            DO CV_INOD = 1, Mdims%cv_nonods
                DO IPHASE = 1, Mdims%nphase
                    RSUM( IPHASE ) = SUM (ComponentMassFraction % val (:, IPHASE, CV_INOD) )
                END DO
                DO IPHASE = 1, Mdims%nphase
                    ComponentMassFraction % val (:, IPHASE, CV_INOD) = ComponentMassFraction % val (:, IPHASE, CV_INOD) / RSUM( IPHASE )
                END DO
            END DO
        end if



        do icomp = 1, Mdims%ncomp

            if ( have_option( '/material_phase[' // int2str( Mdims%nstate - Mdims%ncomp ) // &
                ']/is_multiphase_component/KComp_Sigmoid' ) .and. Mdims%nphase > 1 ) then

                call Calculate_ComponentAbsorptionTerm( packed_state, icomp, ndgln%cv, &
                     Mdims, PhaseDensity%val, Porosity_field%val, mass_ele, Component_Absorption%val )

                do cv_nodi = 1, Mdims%cv_nonods
                    Loop_Phase_SourceTerm1: do iphase = 1, Mdims%nphase
                        Loop_Phase_SourceTerm2: do jphase = 1, Mdims%nphase
                            PhaseVolumeFractionComponentSource%val( 1, iphase, cv_nodi ) = &
                                PhaseVolumeFractionComponentSource%val( 1, iphase, cv_nodi ) &
                                - Component_Absorption%val( iphase, jphase, cv_nodi ) &
                                * ComponentMassFraction%val( icomp, jphase, cv_nodi ) &
                                / ComponentDensity%val( icomp, iphase, cv_nodi )
                        end do Loop_Phase_SourceTerm2
                    end do Loop_Phase_SourceTerm1
                end do

            end if

            ! For compressibility
            do iphase = 1, Mdims%nphase
                do cv_nodi = 1, Mdims%cv_nonods
                    PhaseVolumeFractionComponentSource%val( 1, iphase, cv_nodi ) =  &
                        PhaseVolumeFractionComponentSource%val( 1, iphase, cv_nodi ) &
                        + MeanPoreCV%val( 1, cv_nodi ) * OldComponentMassFraction%val( icomp, iphase, cv_nodi ) &
                        * ( OldComponentDensity%val( icomp, iphase, cv_nodi ) - ComponentDensity%val( icomp, iphase, cv_nodi ) ) &
                        * old_saturation_field%val( 1, IPHASE, Mdims%cv_nonods ) &
                        / ( ComponentDensity%val( ICOMP, IPHASE, CV_NODI ) * DT )
                end do
            end do

        end do ! icomp

        if ( have_option( '/material_phase[' // int2str( Mdims%nstate - Mdims%ncomp ) // &
            ']/is_multiphase_component/Comp_Sum2One' ) .and. ( Mdims%ncomp > 1 ) ) then
            call Cal_Comp_Sum2One_Sou( packed_state, Mdims )
        end if

        deallocate(theta_flux, one_m_theta_flux, theta_flux_j, one_m_theta_flux_j)


        !First, impose physical constrains to the saturation (important to update halos here)
        if (Mdims%n_in_pres > 1) call Set_Saturation_to_sum_one(mdims, packed_state, state, &
                                                              do_not_update_halos = .false. )

      contains

        !!>@brief: This subroutines eliminates the oscillations in the component that are bigger than a
        !certain tolerance and also sets the component to be between bounds
        !> WARNING: It is currently not working well...better to use Chris' method which uses a RHS
        subroutine ensure_components_sum_to_one(packed_state)
            Implicit none
            !Global variables
            type( state_type ), intent(inout) :: packed_state
            !Local variables
            integer :: iphase, cv_nod, i_start, i_end, ipres
            real :: correction, sum_of_components
            type(tensor_field), pointer :: ComponentMassFraction

            ComponentMassFraction  => extract_tensor_field( packed_state, "PackedComponentMassFraction" )

            !Impose sat to be between bounds for blocks of ComponentMassFraction (this is for multiple pressure, otherwise there is just one block)
            do ipres = 1, Mdims%npres
                i_start = 1 + (ipres-1) * Mdims%nphase/Mdims%npres
                i_end = ipres * Mdims%nphase/Mdims%npres
                !Set ComponentMassFraction to be between bounds
                do iphase = i_start, i_end
                  do cv_nod = 1, size(ComponentMassFraction%val, 3)
                      sum_of_components = sum(ComponentMassFraction%val(:, iphase, cv_nod))
                      correction = (1.0 - sum_of_components)
                      !Spread the error to all the components weighted by their presence in that CV
                      !Increase the range to look for solutions by allowing oscillations below 0.1 percent
                      if (abs(correction) > 1d-3) ComponentMassFraction%val(:,iphase, cv_nod) = &
                            (ComponentMassFraction%val(:,iphase, cv_nod) * (1.0 + correction/sum_of_components))
                      !Make sure ComponentMassFraction%val is between bounds after the modification
                      ComponentMassFraction%val(:,iphase,cv_nod) =  min(max(0., ComponentMassFraction%val(:, iphase,cv_nod)),1.0)
                  end do
                end do
            end do
            if (IsParallel()) call halo_update(ComponentMassFraction)

        end subroutine ensure_components_sum_to_one




    end subroutine Compositional_Assemble_Solve




    !>@brief: Form the global continuity and the momentum eqns
    !>@param  state Linked list containing all the fields defined in diamond and considered by Fluidity
    !>@param  packed_state Linked list containing all the fields used by IC-FERST, memory partially shared with state
    !!>@param Mdims Data type storing all the dimensions describing the mesh, fields, nodes, etc
    !>@param  CV_GIdims Gauss integration numbers for CV fields
    !>@param FE_GIdims Gauss integration numbers for FE fields
    !>@param  CV_funs Shape functions for the CV mesh
    !>@param  FE_funs Shape functions for the FE mesh
    !>@param  Mspars Sparsity of the matrices
    !>@param  ndgln Global to local variables
    !>@param  Mdisopt Discretisation options
    !>@param  Mmat Matrices for ICFERST
    !>@param  upwnd Sigmas to compute the fluxes at the interphase for porous media
    !>@param  velocity tensor field. Why do we pass it down instead of retrieving it from packed_state... 
    !>@param  pressure tensor fieldWhy do we pass it down instead of retrieving it from packed_state... 
    !>@param  multi_absorp  Absoprtion of associated with the transport field
    !>@param  eles_with_pipe Elements that have a pipe
    !>@param  pipes_aux Information required to define wells
    !>@param X_ALL Coordinates of the nodes
    !>@param velocity_absorption absorption associated to the momentum equation
    !>@param U_SOURCE_ALL momentum source term on the element mesh
    !>@param U_SOURCE_CV_ALL momentum source term on the CV mesh (used for gravity)
    !>@param U_ALL, UOLD_ALL velocity and velocity old (same as velocity as far as I know...)
    !>@param CV_P Pressure in using Finite volumes (same as Pressure for DCVFEM)
    !>@param DEN_ALL, DENOLD_ALL Densities to be used. This is so we can impose the boussinesq approximation
    !>@param DERIV Derivatives of the Equation of state against pressure so we can generate the M matrix for compressible flow
    !>@param  DT Time step size
    !>@param MASS_MN_PRES ???
    !>@param MASS_ELE mass of the elements
    !>@param got_free_surf ???
    !>@param   MASS_SUF area of the surface of the elements???
    !>@param  SUF_SIG_DIAGTEN_BC Like upwnd but for the boundary
    !>@param V_SOURCE Source term
    !>@param  VOLFRA_PORE     Porosity field (Mdims%npres,Mdims%totele)
    !>@param Courant_number Global courant number and shock front courant number   
    !>@param DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B Variables related with the coupling between the wells and reservoir domain
    !>@param JUST_BL_DIAG_MAT ???
    !>@param  UDEN_ALL, UDENOLD_ALL Densities to be used. This is so we can impose the boussinesq approximation??
    !>@param  UDIFFUSION_ALL, UDIFFUSION_VOL_ALL Parameters associated to the diffusion of the velocity for navier-stokes
    !>@param   IGOT_THETA_FLUX, THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J ????
    !>@param IPLIKE_GRAD_SOU ??
    !>@param FEM_continuity_equation This is to use the divergence as the transpose of the gradient matrix. Not conservative but more stable than simple CVFEM
    !>@param  calculate_mass_delta This is to compute mass conservation 
    !>@param  outfluxes variable containing the outfluxes information
    !>@param DIAG_BIGM_CON, BIGM_CON To assemble the momentum equation of the Navier-stokes equation
    SUBROUTINE CV_ASSEMB_FORCE_CTY( state, packed_state, &
        Mdims, CV_GIdims, FE_GIdims, CV_funs, FE_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd, &
        velocity, pressure, multi_absorp, eles_with_pipe, pipes_aux, &
        X_ALL, velocity_absorption, U_SOURCE_ALL, U_SOURCE_CV_ALL, &
        U_ALL, UOLD_ALL, &
        CV_P, DEN_ALL, DENOLD_ALL, DERIV, &
        DT, &
        MASS_MN_PRES, MASS_ELE,&
        got_free_surf,  MASS_SUF, &
        SUF_SIG_DIAGTEN_BC, &
        V_SOURCE, VOLFRA_PORE, Courant_number, &
        DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
        JUST_BL_DIAG_MAT, &
        UDEN_ALL, UDENOLD_ALL, UDIFFUSION_ALL, UDIFFUSION_VOL_ALL, &
        IGOT_THETA_FLUX, &
        THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, &
        IPLIKE_GRAD_SOU, &
        FEM_continuity_equation, calculate_mass_delta, outfluxes, DIAG_BIGM_CON, BIGM_CON) !-ao
        implicit none
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
        type(multi_absorption), intent(inout) :: multi_absorp
        type(pipe_coords), dimension(:), intent(in):: eles_with_pipe
        type (multi_pipe_package), intent(in) :: pipes_aux
        INTEGER, intent( in ) :: IGOT_THETA_FLUX, IPLIKE_GRAD_SOU
        LOGICAL, intent( in ) :: got_free_surf,FEM_continuity_equation
        real, dimension(:,:), intent(in) :: X_ALL
        REAL, DIMENSION( :, :, : ), intent( in ) :: velocity_absorption
        type( multi_field ), intent( in ) :: U_SOURCE_ALL
        REAL, DIMENSION( :, :, : ), intent( in ) :: U_SOURCE_CV_ALL
        REAL, DIMENSION( :, :, : ), intent( in ) :: U_ALL, UOLD_ALL
        REAL, DIMENSION( :, :, : ), intent( in ) :: CV_P
        REAL, DIMENSION(  :, :  ), intent( in ) :: DEN_ALL, DENOLD_ALL
        REAL, DIMENSION(  : , :  ), intent( in ) :: DERIV
        REAL, DIMENSION(  : ,  :   ), intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
        REAL, intent( in ) :: DT
        REAL, DIMENSION(  : , : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
        REAL, DIMENSION(  :, :  ), intent( in ) :: V_SOURCE
        REAL, DIMENSION( :, : ), intent( in ) :: VOLFRA_PORE
        REAL, DIMENSION( : ), intent( inout ) :: MASS_MN_PRES, MASS_ELE
        REAL, DIMENSION( : ), intent( inout ) :: MASS_SUF
        REAL, DIMENSION( :, : ), intent( inout ), allocatable :: DIAG_SCALE_PRES
        REAL, DIMENSION( :, :, : ), intent( inout ), allocatable :: DIAG_SCALE_PRES_COUP, INV_B
        REAL, DIMENSION( :, : ), intent( in ) :: UDEN_ALL, UDENOLD_ALL
        REAL, DIMENSION( :, :, :, : ), intent( in ) :: UDIFFUSION_ALL
        type( multi_field ), intent( in ) :: UDIFFUSION_VOL_ALL
        LOGICAL, intent( inout ) :: JUST_BL_DIAG_MAT
        type (multi_outfluxes), intent(inout) :: outfluxes
        real, dimension(:,:), intent(inout) :: calculate_mass_delta
        REAL, DIMENSION( :,:,:,:,:,:,: ), allocatable ::  DIAG_BIGM_CON
        REAL, DIMENSION( :,:,:,:,:,:,: ), allocatable ::  BIGM_CON
        real, dimension(:), intent(inout) :: Courant_number
        ! Local variables
        REAL, PARAMETER :: v_beta = 1.0
        LOGICAL, PARAMETER :: GETCV_DISC = .FALSE., GETCT= .TRUE., THERMAL= .FALSE.
        REAL, DIMENSION( : ), allocatable ::  dummy_transp
        REAL, DIMENSION( :, : ), allocatable :: THETA_GDIFF
        REAL, DIMENSION( : , : ), allocatable :: DENOLD_OR_ONE
        REAL, DIMENSION( : , : ), target, allocatable :: DEN_OR_ONE
        REAL, DIMENSION( :, : ), allocatable :: MEAN_PORE_CV
        LOGICAL :: GET_THETA_FLUX
        INTEGER :: IGOT_T2, I
        INTEGER :: ELE, U_ILOC, U_INOD, IPHASE, IDIM
        type(tensor_field), pointer :: tracer, density
        REAL, DIMENSION( : , :, : ), pointer :: V_ABSORB => null() ! this is PhaseVolumeFraction_AbsorptionTerm

        ewrite(3,*)'In CV_ASSEMB_FORCE_CTY'
        GET_THETA_FLUX = .FALSE.
        IGOT_T2 = 0
        ALLOCATE( THETA_GDIFF( Mdims%nphase * IGOT_T2, Mdims%cv_nonods * IGOT_T2 )) ; THETA_GDIFF = 0.
        ALLOCATE( MEAN_PORE_CV( Mdims%npres, Mdims%cv_nonods )) ; MEAN_PORE_CV = 0.
        allocate( dummy_transp( Mdims%totele ) ) ; dummy_transp = 0.

        !Only the Mass matrix and the RHS of the Darcy equation is assembled here
        CALL porous_assemb_force_cty( packed_state, pressure, &
        Mdims, FE_GIdims, FE_funs, Mspars, ndgln, Mmat, X_ALL, U_SOURCE_CV_ALL)

        ALLOCATE( DEN_OR_ONE( Mdims%nphase, Mdims%cv_nonods )); DEN_OR_ONE = 1.
        ALLOCATE( DENOLD_OR_ONE( Mdims%nphase, Mdims%cv_nonods )); DENOLD_OR_ONE = 1.
        IF ( Mdisopt%volfra_use_theta_flux .or. has_boussinesq_aprox) THEN ! We have already put density in theta... or
          !boussinesq so we do not consider variations of density for the continuity equation
           DEN_OR_ONE = 1.
           DENOLD_OR_ONE = 1.
        ELSE
           DEN_OR_ONE = DEN_ALL
           DENOLD_OR_ONE = DENOLD_ALL
        END IF
        ! no q scheme
        tracer=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
        density=>extract_tensor_field(packed_state,"PackedDensity")

        call CV_ASSEMB( state, packed_state, &
            Mdims%n_in_pres, Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd, &
            tracer, velocity, density, multi_absorp, &
            DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
            DEN_OR_ONE, DENOLD_OR_ONE, &
            Mdisopt%v_disopt, Mdisopt%v_dg_vel_int_opt, DT, Mdisopt%v_theta, v_beta, &
            SUF_SIG_DIAGTEN_BC, &
            DERIV, CV_P, &
            V_SOURCE, V_ABSORB, VOLFRA_PORE, &
            GETCV_DISC, GETCT, &
            IGOT_T2, IGOT_THETA_FLUX, GET_THETA_FLUX, Mdisopt%volfra_use_theta_flux, &
            THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
            MEAN_PORE_CV, &
            MASS_MN_PRES, THERMAL, &
            got_free_surf,  MASS_SUF, &
            dummy_transp, &
            eles_with_pipe = eles_with_pipe, pipes_aux = pipes_aux, &
            calculate_mass_delta = calculate_mass_delta, outfluxes = outfluxes,&
            Courant_number = Courant_number)

            !Make the inf norm of the Courant number across cpus
            if (IsParallel()) then
              call allmax(Courant_number(1)); call allmax(Courant_number(2))
           end if

        ewrite(3,*)'Back from cv_assemb'
        deallocate( DEN_OR_ONE, DENOLD_OR_ONE )
        DEALLOCATE( THETA_GDIFF )
        DEALLOCATE( MEAN_PORE_CV )
        ewrite(3,*) 'Leaving CV_ASSEMB_FORCE_CTY'

      contains

        !!>@brief: If using the CV formulation and porous media, the momemtum equation and the mass matrix are much simpler
        !> and faster to compute than for Navier-Stokes. Therefore, here the RHS and the Mass matrix are computed for this case
        SUBROUTINE porous_assemb_force_cty( packed_state, pressure,&
            Mdims, FE_GIdims, FE_funs, Mspars, ndgln, Mmat, X_ALL, U_SOURCE_CV_ALL)
            implicit none
            type( state_type ), intent( inout ) :: packed_state
            type(multi_dimensions), intent(in) :: Mdims
            type(multi_GI_dimensions), intent(in) :: FE_GIdims
            type(multi_shape_funs), intent(in) :: FE_funs
            type (multi_sparsities), intent(in) :: Mspars
            type(multi_ndgln), intent(in) :: ndgln
            type (multi_matrices), intent(inout) :: Mmat
            REAL, DIMENSION( :, : ), intent( in ) :: X_ALL
            real, dimension(:,:,:), intent(in) :: U_SOURCE_CV_ALL
            type( tensor_field ), intent(in) :: pressure
            ! Local Variables
            integer :: CV_ILOC, CV_JLOC, GI, ELE, U_ILOC, U_INOD, CV_INOD, stat
            real, dimension(FE_GIdims%cv_ngi, Mdims%u_nloc) :: UFEN_REVERSED
            real, dimension(FE_GIdims%cv_ngi, Mdims%cv_nloc) :: CVN_REVERSED
            !Variables for capillary pressure
            integer :: MAT_ILOC, MAT_ILOC2, CV_SILOC, iface, mat_inod, mat_inod2, mat_siloc, x_inod, ele2
            integer, dimension(Mdims%mat_nloc) :: MAT_OTHER_LOC
            integer, dimension(:,:), allocatable :: FACE_ELE
            real :: sarea
            real, dimension(:), allocatable :: NORMX_ALL, sdetwe
            real, dimension(:, :), allocatable ::  XL_ALL, XSL_ALL, SNORMXN_ALL!should remove all local conversions
            real, dimension(:, :, :), allocatable :: LOC_U_RHS !should remove all local conversions
            type(tensor_field), pointer :: CapPressure
            !Diamond options
            logical, save :: options_read = .false., Bubble_element_active,capillary_pressure_activated, Diffusive_cap_only, gravity_on
            real, save :: gravty = 0.0
            character(len=FIELD_NAME_LEN) :: element_type_name
            !###Shape function calculation###
            type(multi_dev_shape_funs) :: Devfuns
            !Parallel variables
            logical :: skip
            integer :: nb
            integer, dimension(:), pointer :: neighbours

            !Check if we have a bubble element to calculate the mass matrix differently
            call get_option("/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type", element_type_name, default = "lagranian")
            Bubble_element_active = trim(element_type_name)=="bubble"
            !Sprint_to_do: improvement for the bubble element pair, use only here the
            !highest quadrature set (45 GI points), or even better, only for the mass matrix calculation, elsewhere go for the 11 GI points

            !Prepare Devfuns
            call allocate_multi_dev_shape_funs(FE_funs, Devfuns)
            !Read options from diamond, only once
            if (.not.options_read) then
                gravity_on = have_option("/physical_parameters/gravity/magnitude")
                call get_option( "/physical_parameters/gravity/magnitude", gravty, default = 0. )
                !Check capillary options
                capillary_pressure_activated = have_option_for_any_phase('/multiphase_properties/capillary_pressure', Mdims%nphase)
                Diffusive_cap_only = have_option_for_any_phase('/multiphase_properties/capillary_pressure/Diffusive_cap_only', Mdims%nphase)
            end if

            CapPressure => extract_tensor_field( packed_state, "PackedCapPressure", stat )

            DO U_ILOC=1,Mdims%u_nloc
                DO GI=1,FE_GIdims%cv_ngi
                    UFEN_REVERSED(GI,U_ILOC) = FE_funs%ufen(U_ILOC,GI)
                END DO
            END DO
            DO CV_ILOC=1,Mdims%cv_nloc
                DO GI=1,FE_GIdims%cv_ngi
                    !############THIS IS A FIX FOR GRAVITY######################
                    CVN_REVERSED(GI,CV_ILOC)  = FE_funs%cvfen(CV_ILOC,GI)
                END DO
            END DO

            !Initialise RHS
            Mmat%U_RHS = 0.0

            if (gravity_on .or. .not.Mmat%Stored) then
                DO ELE = 1, Mdims%totele ! VOLUME integral
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
                            ! if (skip) cycle
                        end if
                    end if
                    ! Calculate DevFuns%DETWEI,DevFuns%RA,NX,NY,NZ for element ELE
                    call DETNLXR_PLUS_U(ELE, X_ALL, ndgln%x, FE_funs%cvweight, &
                        FE_funs%cvfen, FE_funs%cvfenlx_all, FE_funs%ufenlx_all, Devfuns)
                    !Assemble lumped mass matrix (only necessary at the beggining and after adapt)
                    !this has to go, the mass matrix should not be assembled at all as it can be done on-the-fly so long
                    !we have the mass of each element
                    if (.not.Mmat%Stored) then
                      if (.not. Bubble_element_active) then
                        !Use the method based on diagonal scaling for lagrangian elements
                        call get_diagonal_mass_matrix(ELE, Mdims, DevFuns, Mmat)
                      else !Use the row-sum method for P1DGBLP1DG(CV)
                        call get_massMatrix(ELE, Mdims, DevFuns, Mmat, X_ALL, UFEN_REVERSED)
                      end if
                    end if
                    !Introduce gravity right-hand-side
                    do U_ILOC = 1, Mdims%u_nloc
                        U_INOD = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                        do CV_JLOC = 1, Mdims%cv_nloc
                            CV_INOD = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_JLOC )
                            Mmat%U_RHS( :, :, U_INOD ) = Mmat%U_RHS( :, :, U_INOD ) + &
                                SUM( UFEN_REVERSED( :, U_ILOC ) * CVN_REVERSED( :, CV_JLOC ) * DevFuns%DETWEI( : ) )*&!shape functions
                                U_SOURCE_CV_ALL( :, :, CV_INOD )!gravity term
                        end do
                    end do
                END DO
            end if

            !! *************************loop over surfaces*********************************************
            ! at some pt we need to merge these 2 loops but there is a bug when doing that!!!!!
            !it does not work because some things, like MASS_ELE(ele2) require to have gone through all the elements first
            if (capillary_pressure_activated.and..not. Diffusive_cap_only)  then
                allocate( FACE_ELE( FE_GIdims%nface, Mdims%totele ), sdetwe(FE_GIdims%sbcvngi))
                allocate( XL_ALL(Mdims%ndim,Mdims%cv_nloc), XSL_ALL(Mdims%ndim,Mdims%cv_snloc) )
                allocate( NORMX_ALL(Mdims%ndim), SNORMXN_ALL(Mdims%ndim,FE_GIdims%sbcvngi) )
                allocate(LOC_U_RHS(Mdims%ndim, Mdims%nphase, Mdims%u_nloc))
                ! Calculate FACE_ELE
                CALL CALC_FACE_ELE( FACE_ELE, Mdims%totele, Mdims%stotel, FE_GIdims%nface, &
                    Mspars%ELE%fin, Mspars%ELE%col, Mdims%cv_nloc, Mdims%cv_snloc, Mdims%cv_nonods, ndgln%cv, ndgln%suf_cv, &
                    FE_funs%cv_sloclist, Mdims%x_nloc, ndgln%x )
                !This has to go once the new method to implement capillary pressure in the saturation equation is done
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
                        ELE2  =MAX( 0, FACE_ELE( IFACE, ELE ))
                        ! Create local copy of X_ALL
                        DO CV_ILOC = 1, Mdims%cv_nloc
                            X_INOD = ndgln%x( (ELE-1)*Mdims%x_nloc + CV_ILOC )
                            XL_ALL(:,CV_ILOC) = X_ALL( :, X_INOD )
                        END DO
                        ! Create local copy of X_ALL for surface nodes
                        DO CV_SILOC = 1, Mdims%cv_snloc
                            CV_ILOC = FE_funs%cv_sloclist( IFACE, CV_SILOC )
                            X_INOD = ndgln%x( (ELE-1)*Mdims%x_nloc + CV_ILOC )
                            XSL_ALL( :, CV_SILOC ) = X_ALL( :, X_INOD )
                        END DO
                        !Obtain normal
                        CALL DGSIMPLNORM_ALL( Mdims%cv_nloc, Mdims%cv_snloc, Mdims%ndim, &
                            XL_ALL, XSL_ALL, NORMX_ALL )
                        CALL DGSDETNXLOC2_ALL( Mdims%cv_snloc, FE_GIdims%sbcvngi, Mdims%ndim, XSL_ALL,&
                            FE_funs%sbcvfen, FE_funs%sbcvfenslx, FE_funs%sbcvfensly, &
                            FE_funs%sbcvfeweigh, SDETWE, SAREA, SNORMXN_ALL, NORMX_ALL )


                        !Surface integral along an element
                        IF(ELE2 > 0) THEN
                            ! ***********SUBROUTINE DETERMINE_OTHER_SIDE_FACE - START************
                            MAT_OTHER_LOC=0
                            DO MAT_SILOC = 1, Mdims%cv_snloc
                                MAT_ILOC = FE_funs%cv_sloclist( IFACE, MAT_SILOC )
                                MAT_INOD = ndgln%x(( ELE - 1 ) * Mdims%mat_nloc + MAT_ILOC )
                                DO MAT_ILOC2 = 1, Mdims%mat_nloc
                                    MAT_INOD2 = ndgln%x(( ELE2 - 1 ) * Mdims%mat_nloc + MAT_ILOC2 )
                                    IF ( MAT_INOD2 == MAT_INOD ) THEN
                                        MAT_OTHER_LOC( MAT_ILOC )=MAT_ILOC2
                                        exit
                                    END IF
                                END DO
                            END DO
                           ! ***********SUBROUTINE DETERMINE_OTHER_SIDE_FACE - END************
                        END IF
                        !Calculate all the necessary stuff and introduce the CapPressure in the RHS
                        call Introduce_Grad_RHS_field_term (&
                            packed_state, Mdims, Mmat, CapPressure%val, FE_funs, Devfuns, X_ALL, LOC_U_RHS, ele, &
                            ndgln%cv, ndgln%x, ele2, iface,&
                            sdetwe, SNORMXN_ALL, FE_funs%u_sloclist( IFACE, : ), FE_funs%cv_sloclist( IFACE, : ), MAT_OTHER_LOC )
                    END DO Between_Elements_And_Boundary
                    !! *************************end loop over surfaces*********************************************
                    ! copy local memory
                    DO U_ILOC = 1, Mdims%u_nloc
                        U_INOD = ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                        Mmat%U_RHS( :, :, U_INOD ) = Mmat%U_RHS( :, :, U_INOD ) + LOC_U_RHS( :, :, U_ILOC )
                    END DO
                END DO Loop_Elements2
                deallocate (FACE_ELE, XL_ALL, XSL_ALL, NORMX_ALL, SNORMXN_ALL, LOC_U_RHS)
            end if
            call deallocate_multi_dev_shape_funs(Devfuns)
        END SUBROUTINE porous_assemb_force_cty

        !!>@brief:This subroutine creates a mass matrix using various approaches
        !>Here no homogenisation can be performed.
        !> NOTE: FOR THE TIME BEING ONLY ROW_SUM IS ACTIVATED HERE, AND GET_POROUS_MASS_MATRIX IS KEPT FOR THE DIAGONAL SCALING METHOD
        subroutine get_massMatrix(ELE, Mdims, DevFuns, Mmat, X_ALL, UFEN_REVERSED)
              implicit none
              integer, intent(in) :: ELE
              type(multi_dimensions), intent(in) :: Mdims
              type(multi_dev_shape_funs), intent(in) :: Devfuns
              type (multi_matrices), intent(inout) :: Mmat
              REAL, DIMENSION( :, : ), intent( in ) :: X_ALL, UFEN_REVERSED
              !Local variables
              integer:: I, J, U_JLOC, U_ILOC, GI, JPHASE, JDIM, IPHASE, idim, JPHA_JDIM, IPHA_IDIM
              REAL, DIMENSION ( Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase, Mdims%u_nloc, Mdims%u_nloc ) :: NN_SIGMAGI_ELE ! element mass matrix
              REAL, DIMENSION ( Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase, FE_GIdims%cv_ngi ) :: SIGMAGI
              logical, parameter :: row_sum_method = .true.
              REAL, SAVE :: scaling_vel_nodes = -1.
              REAL, SAVE :: Tau = -1.

                SIGMAGI = 0.
                DO IPHA_IDIM = 1, Mdims%ndim * Mdims%nphase
                    SIGMAGI( IPHA_IDIM, IPHA_IDIM, : ) = 1.0
                end do
                  !Initialise
                  NN_SIGMAGI_ELE = 0.
                  DO U_JLOC = 1, Mdims%u_nloc
                      DO U_ILOC = 1, Mdims%u_nloc
                          DO GI = 1, FE_GIdims%cv_ngi
                                NN_SIGMAGI_ELE(:, :, U_ILOC, U_JLOC ) = NN_SIGMAGI_ELE(:, :, U_ILOC, U_JLOC ) &
                                + UFEN_REVERSED(GI, U_ILOC) * UFEN_REVERSED(GI, U_JLOC) * DevFuns%DETWEI( GI ) * SIGMAGI( :, :, GI )
                          END DO
                      END DO
                  END DO

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
                                      IF (row_sum_method) then!sprint_to_do despite hard coded I don't like this if within these loops
                                                              !decide which method to keep and just leave that one
                                        Mmat%PIVIT_MAT( I, I, ELE ) =  Mmat%PIVIT_MAT( I, I, ELE ) + &
                                        NN_SIGMAGI_ELE(IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC )
                                      else
                                        Mmat%PIVIT_MAT( I, J, ELE ) =  &
                                        NN_SIGMAGI_ELE(IPHA_IDIM, JPHA_JDIM, U_ILOC, U_JLOC )
                                      end if
                                  END DO
                              END DO
                          END DO
                      END DO
                  END DO
              END DO

            end subroutine

    END SUBROUTINE CV_ASSEMB_FORCE_CTY

    !>@brief: This subroutine performs and introduces the gradient of a RHS field (Capillary pressure for example)
    !> for the momentum equation. It has two options (hard-coded) integrations by parts (activated) or voluemtric integration, 
    !> For capillary pressure: The capillary pressure is a term introduced as a RHS which affects the effective velocity
    !>@param  packed_state Linked list containing all the fields used by IC-FERST, memory partially shared with state
    !!>@param Mdims Data type storing all the dimensions describing the mesh, fields, nodes, etc
    !>@param  Mmat Matrices for ICFERST
    !>@param RHS_field The values to introduce, for capilalry pressure these are the evaluted values of CapPressure%val
    !>@param  FE_funs Shape functions for the FE mesh
    !>@param Devfuns Derivatives of the FE_funs
    !>@param X_ALL Coordinates of the nodes
    !>@param LOC_U_RHS RHS INOUT field with the updated field
    !>@param ele Current element
    !>@param cv_ndgln Global to local variables of the CV mesh
    !>@param x_ndgln Global to local variables of the nodes
    !>@param ele2 neighbouring element
    !>@param iface current face being integrated
    !>@param sdetwe edge lenght associated to the current gauss integration point, i.e. edge lenght
    !>@param SNORMXN_ALL normal to the current surface
    !>@param U_SLOC2LOC Surface local to local for velocity nodes
    !>@param CV_SLOC2LOC Surface local to local for cv nodes
    !>@param MAT_OTHER_LOC for sub-CV the neighbouring one
 subroutine Introduce_Grad_RHS_field_term (packed_state, Mdims, Mmat, RHS_field, FE_funs, Devfuns, &
     X_ALL, LOC_U_RHS, ele, cv_ndgln, x_ndgln,&
     ele2, iface, sdetwe, SNORMXN_ALL, U_SLOC2LOC, CV_SLOC2LOC, MAT_OTHER_LOC)
     Implicit none
     type( state_type ), intent( inout ) :: packed_state
     type(multi_dimensions), intent(in) :: Mdims
     type (multi_matrices), intent(inout) :: Mmat
     type(multi_shape_funs), intent(in) :: FE_funs
     integer, intent(in) :: ele, iface, ele2
     integer, dimension(:), intent(in) :: cv_ndgln, x_ndgln
     REAL, DIMENSION ( :, :, : ), intent(inout) :: LOC_U_RHS
     real, dimension(:,:) :: X_ALL
     real, dimension(:,:), intent(in) :: SNORMXN_ALL
     integer, dimension(:), intent(in) :: U_SLOC2LOC, CV_SLOC2LOC, MAT_OTHER_LOC
     real, dimension(:), intent(in) :: sdetwe
     type(multi_dev_shape_funs), intent(inout) :: Devfuns
     real, dimension(:,:,:), intent(in) :: RHS_field
     !Local parameters
     !!Use a finite element projection of the RHS_field, it can only be false for PnDGPn(DG) elements
     logical, parameter :: Cap_to_FEM = .false.
     !Use integration by parts to introduce the RHS_field, otherwise it uses the integration by parts twice approach
     logical, parameter :: Int_by_part_CapPress = .true.
     !The combination Cap_to_FEM = .false. and Int_by_part_CapPress = .true. seems better for DCVFEM and the same for CVFEM
     !Local variables
     integer :: iphase, cv_inod, u_siloc, cv_jloc,&
         CV_SJLOC, u_iloc, cv_Xnod
     real, pointer, dimension(:, :) :: CV_Bound_Shape_Func
     real, pointer, dimension(:, :) :: CV_Shape_Func
     real, dimension(Mdims%NDIM) :: NMX_ALL
     logical :: DISC_PRES ! discontinuous pressure flag, only perform volumetric integral for the continuous pressure method, otherwise an extra surface intergal is needed
     ! following integration by parts twice which introduces the jump condition see Gomes et al 2016

     !Retrieve derivatives of the shape functions
     call DETNLXR_PLUS_U(ELE, X_ALL, X_NDGLN, FE_funs%cvweight, &
     FE_funs%cvfen, FE_funs%cvfenlx_all, FE_funs%ufenlx_all, Devfuns)

     ! discontinuous pressure flag
     DISC_PRES = ( Mdims%cv_nonods == Mdims%totele * Mdims%cv_nloc )

     !Project to FEM
     CV_Bound_Shape_Func => FE_funs%sbcvfen
     CV_Shape_Func => FE_funs%cvfen

     !Integration by parts
     if (Int_by_part_CapPress .or. .not. CAP_to_FEM) then
       if (iface == 1) then!The volumetric term is added just one time
         !Firstly we add the volumetric integral
         DO U_ILOC = 1, Mdims%u_nloc
           DO CV_JLOC = 1, Mdims%cv_nloc
             ! -Integral(FE_funs%cvn RHS_field Grad FE_funs%ufen dV)
             CV_INOD = CV_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_JLOC )
             DO IPHASE = 1, Mdims%nphase
               LOC_U_RHS( :, IPHASE, U_ILOC ) = LOC_U_RHS( :, IPHASE, U_ILOC ) &
               !(FE_funs%cvn Grad FE_funs%ufen)
               + matmul(Devfuns%UFENX_ALL(:,U_ILOC,:),CV_Shape_Func( CV_JLOC, : ) *Devfuns%detwei )&
               !RHS_field
               * RHS_field(1, IPHASE, CV_INOD)
             END DO
           end do
         end do
       end if
       !Performing the surface integral, -Integral(FE_funs%cvn RHS_field ᐁFE_funs%ufen dV)
       DO U_SILOC = 1, Mdims%u_snloc
         U_ILOC = U_SLOC2LOC( U_SILOC )
         DO CV_SJLOC = 1, Mdims%cv_snloc
           CV_JLOC = CV_SLOC2LOC( CV_SJLOC )
           CV_INOD = CV_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_JLOC )
           NMX_ALL = matmul(SNORMXN_ALL( :, : ), FE_funs%sbufen( U_SILOC, : ) &
           * CV_Bound_Shape_Func( CV_SJLOC, : ) * SDETWE( : ))
           if (ELE2 > 0) then!If neighbour then we get its value to calculate the average
             cv_Xnod = CV_NDGLN( ( ELE2 - 1 ) * Mdims%cv_nloc + MAT_OTHER_LOC(CV_JLOC) )
           else !If no neighbour then we use the same value.
             cv_Xnod = CV_INOD
           end if
           do iphase = 1, Mdims%nphase
             LOC_U_RHS( :, IPHASE, U_ILOC) =  LOC_U_RHS( :, IPHASE, U_ILOC ) &
             - NMX_ALL(:) * 0.5*(RHS_field(1, iphase, CV_INOD)+RHS_field(1, iphase, cv_Xnod))
           end do
         end do
       end do
     else !Volumetric integration only (requires the RHS_field to be in FEM)
       if (iface ==1) then!The volumetric term is added just one time
         DO U_ILOC = 1, Mdims%u_nloc
           DO CV_JLOC = 1, Mdims%cv_nloc
             ! Integral(Grad FE_funs%cvn RHS_field FE_funs%ufen dV)
             CV_INOD = CV_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_JLOC )
             DO IPHASE = 1, Mdims%nphase
               LOC_U_RHS( :, IPHASE, U_ILOC ) = LOC_U_RHS( :, IPHASE, U_ILOC ) &
               !(Grad FE_funs%cvn FE_funs%ufen)
               - matmul(Devfuns%CVFENX_ALL(:,CV_JLOC,:),FE_funs%ufen( U_ILOC, : ) *Devfuns%DETWEI )&
               !RHS_field
               * RHS_field(1, IPHASE, CV_INOD)
             END DO
           end do
         end do
       end if
       !Get neighbouring nodes!SPRINT_TO_DO Use beta instead of the 0.5 for this. Also it seems that dPc/dS grad S is more stable
       !Also if not dicsontinuous formulation do not perform this operation
       if (DISC_PRES) then
         !Get neighbouring nodes
         !Performing the surface integral, Integral(FE_funs%cvn (Average RHS_field) ᐁFE_funs%ufen dV)
         DO U_SILOC = 1, Mdims%u_snloc
           U_ILOC = U_SLOC2LOC( U_SILOC )
           DO CV_SJLOC = 1, Mdims%cv_snloc
             CV_JLOC = CV_SLOC2LOC( CV_SJLOC )
             CV_INOD = CV_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_JLOC )
             NMX_ALL = matmul(SNORMXN_ALL( :, : ), FE_funs%sbufen( U_SILOC, : ) * FE_funs%sbcvfen( CV_SJLOC, : ) * SDETWE( : ))
             if (ELE2 > 0) then!If neighbour then we get its value to calculate the average
               cv_Xnod = CV_NDGLN( ( ELE2 - 1 ) * Mdims%cv_nloc + MAT_OTHER_LOC(CV_JLOC) )
             else !If no neighbour then we use the same value.
               cv_Xnod = CV_INOD
             end if
             do iphase = 1, Mdims%nphase
               LOC_U_RHS( :, IPHASE, U_ILOC) =  LOC_U_RHS( :, IPHASE, U_ILOC ) &
               + NMX_ALL(:) * 0.5* (RHS_field(1, iphase, CV_INOD) - RHS_field(1, iphase, cv_Xnod))
             end do
           end do
         end do
       end if
     end if

 end subroutine Introduce_Grad_RHS_field_term



    !!>@brief: This subroutine calculates the overrelaxation parameter (Vanishing relaxation) we introduce in the saturation equation
    !> Overrelaxation has to be alocate before calling this subroutine its size is cv_nonods
    !> For more information read: doi.org/10.1016/j.cma.2019.07.004
    !>@param  state Linked list containing all the fields defined in diamond and considered by Fluidity
    !>@param  packed_state Linked list containing all the fields used by IC-FERST, memory partially shared with state
    !!>@param Mdims Data type storing all the dimensions describing the mesh, fields, nodes, etc
    !>@param  ndgln Global to local variables
    !>@param Overrelaxation INOUT Field containing the outcome of this subroutine
    !>@param Phase_with_Pc Phase that will contain the Overrelaxation field
    !>@param totally_min_max Min max of the field (to normalise for transport)
    !>@param for_transport True if the field is not saturation
 subroutine getOverrelaxation_parameter(state, packed_state, Mdims, ndgln, Overrelaxation, Phase_with_Pc, totally_min_max, for_transport)
     implicit none
     type( state_type ), dimension( : ), intent( inout ) :: state
     type( state_type ), intent(inout) :: packed_state
     type (multi_dimensions), intent(in) :: Mdims
     type(multi_ndgln), intent(in) :: ndgln
     real, dimension(:), intent(inout) :: Overrelaxation
     integer, intent(inout) :: Phase_with_Pc
     real, dimension(:), optional :: totally_min_max
     logical, optional, intent(in) :: for_transport
     !Local variables
     real, save :: domain_length = -1
     integer, save :: Cap_pressure_relevant = -1
     integer :: iphase, nphase, cv_nodi, cv_nonods, u_inod, cv_iloc, ele, u_iloc, idim
     real :: Pe_aux, parl_max, parl_min, Pe_max, Pe_min
     real, dimension(:), pointer ::Pe, Cap_exp
     logical :: Artificial_Pe
     real, dimension(:,:,:), pointer :: p
     real, dimension(:,:), pointer :: satura, CV_Immobile_Fraction, Cap_entry_pressure, Cap_exponent, X_ALL
     type( tensor_field ), pointer :: Velocity
     type( vector_field ), pointer :: DarcyVelocity
     type( scalar_field ), pointer :: PIPE_Diameter


     !Extract variables from packed_state
     call get_var_from_packed_state(packed_state,FEPressure = P,&
         PhaseVolumeFraction = satura, PressureCoordinate = X_ALL)

    call get_var_from_packed_state(packed_state, CV_Immobile_Fraction = CV_Immobile_Fraction)
     !Initiate local variables
     nphase = size(satura,1)
     cv_nonods = size(satura,2)

     !#######Only apply this method if it has been explicitly invoked through Pe_stab or
     !non-consistent capillary pressure!######
     Phase_with_Pc = -1
     if (.not. have_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation') ) then
         Overrelaxation = 0.0; Phase_with_Pc = -10
         return
     end if
     !If this is for transport, check if we want to apply it
     if (present_and_true(for_transport)) then
        if (.not.have_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation/Vanishing_for_transport')) then
            Overrelaxation = 0.0; Phase_with_Pc = -10
            return
        end if
     end if
     !######################################################################################
    !By default phase 1
    Phase_with_Pc = 1
     !Check capillary pressure options
     do iphase = Nphase, 1, -1!Going backwards since the wetting phase should be phase 1
         !this way we try to avoid problems if someone introduces 0 capillary pressure in the second phase
         if (have_option( "/material_phase["//int2str(iphase-1)//"]/multiphase_properties/capillary_pressure" )) then
             Phase_with_Pc = iphase
         end if
     end do
     Artificial_Pe = .false.
     Cap_exponent => null(); Cap_entry_pressure => null()!Initialize
     if (Phase_with_Pc>0) then
         !Get information for capillary pressure to be used
         if ( (have_option("/material_phase["//int2str(Phase_with_Pc-1)//&
             "]/multiphase_properties/capillary_pressure/type_Brooks_Corey") ) .or. (have_option("/material_phase["//int2str(Phase_with_Pc-1)//&
             "]/multiphase_properties/capillary_pressure/type_Power_Law") ) )then
             call get_var_from_packed_state(packed_state, Cap_entry_pressure = Cap_entry_pressure,&
                 Cap_exponent = Cap_exponent)!no need for the imbibition because we need the derivative which will be zero as it is a constant
         end if
         !If we want to introduce a stabilization term, this one is imposed over the capillary pressure.
         !Unless we are using the non-consistent form of the capillary pressure
         if ( have_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation') ) then
             allocate(Pe(CV_NONODS), Cap_exp(CV_NONODS))
             Artificial_Pe = .true.
             if (present_and_true(for_transport)) then
                call get_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation/Vanishing_for_transport', Pe_aux)
                !This method was designed for fields between 0 and 1, so for transport fields, we need to adjust Pe_aux to ensure consistency
                if (present(totally_min_max)) then
                  !Means that min_max is defined
                  if (abs(totally_min_max(2)) < 1e20) Pe_aux = Pe_aux*max(abs(totally_min_max(2) - totally_min_max(1)),1.0) !THIS means making it HIGHER
                end if
             else
                call get_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation', Pe_aux)
             end if


            !Check if the capillary pressure introduced is important enough to actually trigger the VAD for Capillary pressure
             if ( associated(Cap_entry_pressure) .and. Cap_pressure_relevant < 0) then
                Cap_pressure_relevant = 0
                if (maxval(Cap_entry_pressure)/maxval(P) > 1e-2) Cap_pressure_relevant = 1
            end if

             if (Cap_pressure_relevant > 0) then
                 Cap_exp = 2.0 !Quadratic exponent

             else
                 Cap_exp = 1.!Linear exponent
             end if

             if (Pe_aux<0) then!Automatic set up for Pe
                 !Method based on calculating an entry pressure for a given Peclet number;
                 !Peclet = V * L / Diffusivity; We consider only the entry pressure for the diffusivity
                 !Pe = Vel * L/ Peclet. At present we are using the Darcy Velocity
                 Velocity => extract_tensor_field( packed_state, "PackedVelocity" )
                ! DarcyVelocity => extract_vector_field( state(Phase_with_Pc), "DarcyVelocity" )

                 !Since it is an approximation, the domain length is the maximum distance, we only calculate it once
                 if (domain_length < 0) then
                   do idim = 1, Mdims%ndim
                    !Apples with apples! Check each dimension individually
                     parl_max = maxval(X_ALL(idim,:))
                     parl_min = minval(X_ALL(idim,:))
                     if (IsParallel()) then
                         call allmax(parl_max)
                         call allmin(parl_min)
                     end if
                     domain_length = max(domain_length, abs(parl_max-parl_min))
                   end do
                 end if
                 Pe_aux = abs(Pe_aux)
                  !Obtain an approximation of the capillary number to obtain an entry pressure
                 Pe = 0.
                 do ele = 1, Mdims%totele
                     do u_iloc = 1, Mdims%u_nloc
                         u_inod = ndgln%u(( ELE - 1 ) * Mdims%u_nloc +u_iloc )
                         do cv_iloc = 1, Mdims%cv_nloc
                             cv_nodi = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                             Pe(cv_nodi) = Pe(cv_nodi) + (1./Pe_aux) * (sum(abs(Velocity%val(:,Phase_with_Pc, u_inod)))/real(Mdims%ndim) * domain_length)/real(Mdims%u_nloc)
                         end do
                     end do
                 end do
                 Pe_max = maxval(Pe)
                 Pe_min = minval(Pe)
                 if (IsParallel()) then
                     call allmax(Pe_max)
                     call allmin(Pe_min)
                 end if

                 if (have_option('/numerical_methods/VAD_two_levels')) then
                   !Homogenise using two values only, if below the average value then use the lowest value possible (this is for stratified flows)
                   where (Pe > 0.5*Pe_max + 0.5*Pe_min)
                     Pe = Pe_max
                   elsewhere
                     Pe = Pe_min
                   end where
                 else
                   !Homogenise the value, this seems to be better to avoid problems (method used for the CMAME paper)
                   Pe = (Pe_max+Pe_min)/2.
                 end if
             else
                 Pe = Pe_aux
             end if

         end if

         !For transport we don't use the pseudo capillary pressure function
         if (present_and_true(for_transport)) then
          Overrelaxation = - Pe
          return
         end if
         !Calculate the overrrelaxation parameter, the numbering might be different for Pe and real capillary
         !values, hence we calculate it differently
         if (Artificial_Pe) then
             !Calculate the Overrelaxation
             do ele = 1, Mdims%totele
                 do cv_iloc = 1, Mdims%cv_nloc
                     cv_nodi = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                     Overrelaxation(CV_NODI) =  Get_DevCapPressure(satura(Phase_with_Pc, CV_NODI),&
                         Pe(CV_NODI), Cap_Exp(CV_NODI), CV_Immobile_Fraction(:,CV_NODI), Phase_with_Pc, nphase)
                 end do
             end do
         else
             !Calculate the Overrelaxation
             do ele = 1, Mdims%totele
                 do cv_iloc = 1, Mdims%cv_nloc
                     cv_nodi = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                     Overrelaxation(CV_NODI) =  Get_DevCapPressure(satura(Phase_with_Pc, CV_NODI),&
                         Cap_entry_pressure(Phase_with_Pc, ele), &
                         Cap_exponent(Phase_with_Pc, ele),&
                         CV_Immobile_Fraction(:,CV_NODI), Phase_with_Pc, nphase)
                 end do
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


    !>@brief:A diagonal mass matrix is obtained using the direct lump process, i.e. using the mass of the elements
    !>@param ELE current element
    !!>@param Mdims Data type storing all the dimensions describing the mesh, fields, nodes, etc
    !>@param Devfuns Derivatives of the FE_funs
    !>@param  Mmat Matrices for ICFERST
    subroutine get_diagonal_mass_matrix(ELE, Mdims, DevFuns, Mmat)
        implicit none
        integer, intent(in) :: ELE
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_dev_shape_funs), intent(in) :: Devfuns
        type (multi_matrices), intent(inout) :: Mmat
        !Local variables
        integer:: I, J, vel_degree, pres_degree,&
                U_JLOC, U_ILOC, JPHASE, JDIM, IPHASE, idim
        !Weight parameter, controls the strenght of the homogenization of the velocity nodes per element
        !the bigger the more P0DG it tends to be
        real :: factor, factor_default
        real, save :: lump_vol_factor =-1d25
        real, save :: scaling_vel_nodes = -1

        !Weights for lumping taken from Zienkiewicz vol 1 page 475
        real, parameter :: corner = 3./57.
        real, parameter :: midpoint = 16./57.


        !Weights for direct mass lumping scaling for bubble elements
         real, parameter :: Tau2D = 1.5 !Number obtained in Osman et al. 2019
         real, parameter :: Tau3D = 1.4 !Number obtained in Osman et al. 2019

         !Obtain the scaling factor to spread the volume of the mass matrix
         if (scaling_vel_nodes<0) then
             scaling_vel_nodes = dble(Mdims%u_nloc)
             !Adjust for linear bubble functions, P1(BL)DG
             !We are adding an extra node that adds extra velocity that needs to be compensated
             if (Mdims%ndim==2 .and. Mdims%u_nloc == 4) then
                 scaling_vel_nodes = scaling_vel_nodes/Tau2D
                 lump_vol_factor = 0.!No velocity homogenisation for bubble elements
             end if
             if (Mdims%ndim==3 .and. Mdims%u_nloc == 5) then
                 scaling_vel_nodes = scaling_vel_nodes/Tau3D
                 lump_vol_factor = 0.!No velocity homogenisation for bubble elements
             end if

         end if

        !No homogenisation for Pressure discontinuous formulations
        if (Mdims%mat_nonods == Mdims%p_nonods) lump_vol_factor = 0.

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
            case default !Create the mass matrix normally by distributing the mass evenly between the nodes
                do i=1,size(Mmat%PIVIT_MAT,1)
                    Mmat%PIVIT_MAT(I,I,ELE) = DevFuns%VOLUME/scaling_vel_nodes
                END DO
        end select

        !If pressure boundary element, then we homogenize the velocity in the element
        if (lump_vol_factor<-1d24) then
            call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/polynomial_degree', vel_degree )
            call get_option( '/geometry/mesh::PressureMesh/from_mesh/mesh_shape/polynomial_degree', pres_degree )
            if (pres_degree == 1 .and. vel_degree == 1) then
                factor_default = 1e4
            else
                factor_default = 0.
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
                DO IPHASE = 1, Mdims%nphase!I think this should be mdims%n_in_pres (to not affect the wells)
                    DO IDIM = 1, Mdims%ndim
                        J = IDIM+(IPHASE-1)*Mdims%ndim+(U_JLOC-1)*Mdims%ndim*Mdims%nphase
                        I = IDIM+(IPHASE-1)*Mdims%ndim+(U_ILOC-1)*Mdims%ndim*Mdims%nphase
                        Mmat%PIVIT_MAT(I,I,ELE) = Mmat%PIVIT_MAT(I,I,ELE) + lump_vol_factor
                        Mmat%PIVIT_MAT(I,J,ELE) = Mmat%PIVIT_MAT(I,J,ELE) - lump_vol_factor
                    end do
                end do
            end do
        end do

    end subroutine get_diagonal_mass_matrix


    !>@brief: In this method we assemble and solve the Laplacian system using at least P1 elements
    !> The equation solved is the following: Div sigma Grad X = - SUM (Div K Grad F) with Neuman BCs = 0
    !> where K and F are passed down as a vector. Therefore for n entries the SUM will be performed over n fields
    !> Example: F = (3, nphase, cv_nonods) would include three terms in the RHS and the same for K
    !> If harmonic average then we perform the harmonic average of sigma and K
    !> IMPORTANT: This subroutine requires the PHsparsity to be generated
    !> Note that this method solves considering FE fields. If using CV you may incur in an small error.
    !!>@param Mdims Data type storing all the dimensions describing the mesh, fields, nodes, etc
    !>@param  state Linked list containing all the fields defined in diamond and considered by Fluidity
    !>@param  packed_state Linked list containing all the fields used by IC-FERST, memory partially shared with state
    !>@param  ndgln Global to local variables
    !>@param  Mmat Matrices for ICFERST
    !>@param  Mspars Sparsity of the matrices
    !>@param  CV_funs Shape functions for the CV mesh
    !>@param  CV_GIdims Gauss integration numbers for CV fields
    !>@param Sigma_field Coefficient multipliying the left-hand side term
    !>@param field_name name of the field to be solved for
    !>@param K_fields Coefficients multipliying the RHS terms
    !>@param F_fields RHS fierlds
    !>@param intface_val_type  0 = no interpolation; 1 Harmonic mean; !20 for SP solver, harmonic mean considering charge; negative normal mean
    !>@param solver_path Path to select the type of solver
    subroutine generate_and_solve_Laplacian_system( Mdims, state, packed_state, ndgln, Mmat, Mspars, CV_funs, CV_GIdims, Sigma_field, &
                                                    field_name, K_fields, F_fields, intface_val_type, solver_path)
      implicit none

      type(multi_dimensions), intent( in ) :: Mdims
      type( state_type ), dimension(:), intent( inout ) :: state
      type( state_type ), intent( inout ) :: packed_state
      type(multi_ndgln), intent(in) :: ndgln
      integer, intent(in) :: intface_val_type
      real, dimension(:,:), intent(in) :: Sigma_field
      real, dimension(:,:,:), intent(in) :: K_fields, F_fields
      type(multi_shape_funs), intent(inout) :: CV_funs
      type(multi_sparsities), intent(in) :: Mspars
      type (multi_matrices), intent(inout) :: Mmat
      type(multi_GI_dimensions), intent(in) :: CV_GIdims
      character( len = * ), intent( in ), optional :: field_name
      character(len=option_path_len), optional, intent(in) :: solver_path
      !Local variables
      integer :: i, stat, local_phases
      character(len=option_path_len) :: solver_option_path = "/solver_options/Linear_solver"
      type(scalar_field), pointer  :: solution
      type(vector_field)  :: v_solution

      local_phases = size(F_fields,2)
      !Solver options, if specified
      solver_option_path = "/solver_options/Linear_solver"
      if (present(solver_path)) solver_option_path = solver_path

      !Retrieve the field f interest to have access to the mesh type
      do i = 1, size(state)
        solution => extract_scalar_field(state(i),trim(field_name), stat)
        if (stat == 0) exit
      end do

      !Generate system
      call generate_Laplacian_system( Mdims, packed_state, ndgln, Mmat, Mspars, CV_funs, CV_GIdims, Sigma_field, &
                                          Solution, K_fields, F_fields, intface_val_type)
      !Solve system
      call allocate(v_solution,local_phases,Solution%mesh,"Laplacian_system")
      !Add remove null_space if not bcs specified for the field since we always have natural BCs
      call add_option( trim( solver_option_path ) // "/remove_null_space", stat )
      call zero(v_solution) !; call zero_non_owned(rhs)
      call petsc_solve( v_solution, Mmat%petsc_ACV, Mmat%CV_RHS, option_path = trim(solver_option_path) )
      if (IsParallel()) call halo_update(v_solution)


! call MatView(Mmat%petsc_ACV%M,   PETSC_VIEWER_STDOUT_SELF, i)
! print *, Mmat%CV_RHS%val
      !Copy now back to the existing fields
      do i = 1, size(state)
        solution => extract_scalar_field(state(i),trim(field_name), stat)
        if (stat == 0) solution%val = v_solution%val(i,:)
      end do
      !Remove remove_null_space
      call delete_option( trim( solver_option_path ) // "/remove_null_space", stat )
      call deallocate(v_solution)
      call deallocate( Mmat%CV_RHS ); call deallocate( Mmat%petsc_ACV )
    end subroutine generate_and_solve_Laplacian_system

    !> @author Pablo Salinas, Geraldine Regnier
    !>@brief: In this subroutine, we modify the CMC matrix and the RHS for the pressure to strongly impose the
    !> pressure boundary conditions for wells ONLY. This is required only when using P0DG-P1
    !> If not imposing strong pressure BCs and using P0DG-P1, the formulation for weak pressure BCs
    !> leads to a 0 in the diagonal of the CMC matrix, therefore strong BCs are required.
    !>@param  state Linked list containing all the fields defined in diamond and considered by Fluidity
    !>@param  pipes_aux Information required to define wells
    !!>@param Mdims Data type storing all the dimensions describing the mesh, fields, nodes, etc
    !>@param  Mmat Matrices for ICFERST
    !>@param  ndgln Global to local variables
    !>@param CMC_petsc petsc_csr_matrix field containing the pressure matrix. Matrix to be modified here
    !>@param rhs_p RHS associated to the pressure matrix.
    subroutine impose_strong_bcs_wells(state, pipes_aux, Mdims, Mmat, ndgln, CMC_petsc, pressure,rhs_p)

        ! form pressure matrix CMC using a colouring approach
        type(multi_dimensions), intent(in) :: Mdims
        type (multi_matrices), intent(in) :: Mmat
        type( state_type ), dimension(:), intent( in ) :: state
        type(multi_ndgln), intent(in) :: ndgln
        type(petsc_csr_matrix), intent(inout)::  CMC_petsc
        real, dimension (:,:), intent(inout) :: rhs_p
        type( tensor_field ), intent(in) :: pressure
        type (multi_pipe_package), intent(in) :: pipes_aux
        INTEGER, DIMENSION ( 1, Mdims%npres, surface_element_count(pressure) ) :: WIC_P_BC_ALL
        INTEGER, PARAMETER :: WIC_P_BC_DIRICHLET = 1
        ! Local variables
        INTEGER :: SELE, IPRES, CV_SILOC, CV_NOD, i, ierr, k
        integer, dimension(Mdims%stotel) :: Impose_strong
        type(tensor_field) :: pressure_BCs
        type(scalar_field), pointer :: pipe_diameter

        call get_entire_boundary_condition(pressure,&
            ['weakdirichlet'],&
            pressure_BCs,WIC_P_BC_ALL)
        PIPE_Diameter => EXTRACT_SCALAR_FIELD(state(1), "DiameterPipe")
        !Initialise array to store universal numbering of strong BC positions
        Impose_strong = -1; k = 1
        !Only for wells
        CMC_petsc%is_assembled=.false.
        call assemble( CMC_petsc )
        DO SELE = 1, Mdims%stotel
          DO IPRES = 2, Mdims%npres
            !Find element where we have a pressure BC defined
            if (WIC_P_BC_ALL(1, IPRES, SELE ) == WIC_P_BC_DIRICHLET) then
              !If no flip required or positive pressure
              IF (Mmat%WIC_FLIP_P_VEL_BCS(1,IPRES,SELE) == 0 .or. Mmat%WIC_FLIP_P_VEL_BCS(1,IPRES,SELE) == 10) THEN
                DO CV_SILOC = 1, Mdims%cv_snloc
                  CV_NOD = ndgln%suf_p((SELE-1)*Mdims%cv_snloc + CV_SILOC )
                  !Check if we are in an element that may need strongly imposed BCs
                  if (.not. pipes_aux%impose_strongBCs(CV_NOD)) cycle
                  !Store the universal numbering, required by PETSc
                  Impose_strong(k) = CMC_petsc%row_numbering%gnn2unn( cv_nod, ipres )
                  k = k + 1
                  !Impose P_BC in the right hand side
                  rhs_p(ipres,cv_nod) = pressure_BCs%val(1,IPRES, (SELE-1)*Mdims%cv_snloc + CV_SILOC ) - &
                  pressure%val(1,IPRES, CV_NOD)
                  end do
                end if
            end if
          END DO
        end do

        !Now ensure that all the processors call MatZeroRows consistently
        i = 0
        do k = 1, Mdims%stotel
          if (isparallel()) call allmax(Impose_strong(k))
          if (Impose_strong(k) < 0) exit
          i = i + 1
        end do
        call MatZeroRows(CMC_petsc%M, i, Impose_strong(1:i), 1.0,PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)

        !Re-assemble just in case
        CMC_petsc%is_assembled=.false.
        call assemble( CMC_petsc )
        call deallocate( pressure_BCs )
      end subroutine

    !----------------------------------------------------------------------------------------
    !> @author Vinicius L S Silva
    !> @brief Subroutine that calculates the backtrack_par_factor based on Machine Learning.
    !> The inputs of the Machine learning model are dimensionless numbers and configurations
    !> of the system. This subroutine also generates several dimensionless numbers cv-wise.
    !----------------------------------------------------------------------------------------
#ifdef USING_XGBOOST
    subroutine AI_backtracking_parameters(Mdims, ndgln, packed_state, state, courant_number_in, backtrack_par_factor, overrelaxation, &
                                         & res, resold, outer_nonlinear_iteration, for_transport)
        implicit none
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_ndgln), intent(in) :: ndgln
        type( state_type ), intent(inout) :: packed_state
        type( state_type ), dimension( : ), intent( inout ) :: state
        real, dimension(:), intent(inout) :: courant_number_in
        real, intent(inout) :: backtrack_par_factor
        real, dimension(:), intent(in) :: overrelaxation
        real, intent(in) :: res
        real, intent(in) :: resold
        integer, intent(in) :: outer_nonlinear_iteration
        logical, optional, intent(in) :: for_transport

        !Local variables
        real, save :: backup_shockfront_Courant = 0.
        logical, save :: read_options = .false.
        logical, save :: read_perm = .false.
        !logical, save :: read_por = .false.
        logical, save :: read_sizes = .false.
        logical, save  :: written_file = .false.
        logical, save  :: loaded_file = .false.
        real, dimension(:,:), pointer :: X_ALL, saturation, Imble_frac, cap_entry_pres, end_point_relperm, exponent_relperm
        type(tensor_field), pointer :: permeability, state_viscosity, density, velocity, pressure
        type(vector_field), pointer :: porosity, gravity_direction
        type (vector_field_pointer), dimension(Mdims%nphase) :: Darcy_velocity
        real, dimension(:), allocatable :: total_mobility, counter_subcv, nDarcy_velocity_cvwise
        real, dimension(:,:), allocatable :: Darcy_velocity_cvwise
        real, dimension(Mdims%n_in_pres, Mdims%cv_nonods) :: relperm, viscosity
        real, save :: average_perm_length, average_perm_thickness, average_por, domain_length, domain_thickness
        real :: l1_max, l1_min, l2_max, l2_min, h_max, h_min, aux_shockfront_mobratio, delta_density, diffusivity, gravity_magnitude, sin_alpha
        integer :: iphase, cv_nodi, cv_iloc, ele, total_ele, u_iloc, u_inod, total_cv, Phase_with_Pc, shockfront_counter
        !integer, dimension(8) :: date_values
        !logical :: file_exist
        !Parameters
        real, dimension(:), allocatable :: longitudinal_capillary, transverse_capillary, buoyancy_number, longitudinal_buoyancy, transverse_buoyancy, invPeclet
        logical, save :: gravity, cap_pressure, black_oil, ov_relaxation, one_phase, wells
        real, save :: n_phases, n_components, aspect_ratio
        real :: courant_number, shockfront_courant_number, auxR
        real :: min_total_mobility, max_total_mobility, average_total_mobility
        real :: min_Darcy_velocity, max_Darcy_velocity, average_Darcy_velocity
        real :: min_shockfront_mobratio, max_shockfront_mobratio, average_shockfront_mobratio
        real :: average_longitudinal_capillary, average_transverse_capillary, max_longitudinal_capillary, max_transverse_capillary, min_longitudinal_capillary, min_transverse_capillary
        real :: average_buoyancy_number, average_longitudinal_buoyancy, average_transverse_buoyancy, max_buoyancy_number, max_longitudinal_buoyancy, max_transverse_buoyancy, min_buoyancy_number, min_longitudinal_buoyancy, min_transverse_buoyancy
        real :: average_overrelaxation, max_overrelaxation, min_overrelaxation
        real :: average_invPeclet, max_invPeclet, min_invPeclet
        real :: shockfront_number_ratio, btpf

        !*************************************!
        !!! ***Getting support variables*** !!!
        !*************************************!
        allocate(Darcy_velocity_cvwise(Mdims%ndim, Mdims%cv_nonods))
        allocate(total_mobility(Mdims%cv_nonods), counter_subcv(Mdims%cv_nonods), nDarcy_velocity_cvwise(Mdims%cv_nonods))
        !Retrieve Pressure, as representative of a CV field (lives on the CV mesh although it is FE)
        pressure => extract_tensor_field(packed_state,"PackedFEPressure")
        !Retrieve permeability field - element wise
        permeability => extract_tensor_field(packed_state,"Permeability")
        !Calculate the total number of control volumes
        total_cv = 0
        do cv_nodi = 1, Mdims%cv_nonods
          if (.not. node_owned(pressure, cv_nodi)) cycle
          total_cv = total_cv + 1
        end do
        if (IsParallel()) call allsum(total_cv)
        !Calculate the total number of elements
        total_ele = 0
        do ele = 1, Mdims%totele
          if (IsParallel()) then
            if (.not. element_owned(permeability, ele)) cycle
          end if
          total_ele = total_ele + 1
        end do
        if (IsParallel()) call allsum(total_ele)
        !Permeability averages
        if (.not.read_perm) then
            !Average permeability along Mdims%ndim-1 first dimensions
            average_perm_length = 0.
            if (Mdims%ndim < 3) then
                do ele = 1, Mdims%totele
                  if (IsParallel()) then
                    if (.not. element_owned(permeability, ele)) cycle
                  end if
                  average_perm_length = average_perm_length + permeability%val(1,1,ele)
                end do
                if (IsParallel()) call allsum(average_perm_length)
                average_perm_length = average_perm_length/total_ele
            else
                do ele = 1, Mdims%totele
                  if (.not. element_owned(permeability, ele)) cycle
                  average_perm_length = average_perm_length + permeability%val(1,1,ele)  + permeability%val(2,2,ele)
                end do
                if (IsParallel()) call allsum(average_perm_length)
                average_perm_length = average_perm_length/(total_ele*2)
            end if
            !Average permeability along the last dimensions
            average_perm_thickness = 0.
            do ele = 1, Mdims%totele
              if (.not. element_owned(permeability, ele)) cycle
                average_perm_thickness = average_perm_thickness + permeability%val(Mdims%ndim,Mdims%ndim,ele)
            end do
            if (IsParallel())call allsum(average_perm_thickness)
            average_perm_thickness = average_perm_thickness/total_ele
            read_perm = .true.
        end if
        !Domain length, domain thickness
        !if (.not.read_por) then
        !    !Retrieve porosity field
        !    porosity => extract_vector_field(packed_state,"Porosity")
        !    average_por = 0
        !    do ele = 1, Mdims%totele
        !        average_por = average_por + porosity%val(1, ele)
        !        total_ele = Mdims%totele
        !    if (IsParallel()) then
        !       call allsum(average_por)
        !       call allsum(total_ele)
        !    average_por = average_por/real(total_ele)
        !    read_por = .true.
        !end if

        !Domain lemgth, domain thickness
        if (.not.read_sizes) then
            !Extract variables from packed_state, X_ALL => extract_vector_field( packed_state, "PressureCoordinate" )
            call get_var_from_packed_state(packed_state, PressureCoordinate = X_ALL)
            !Domain length as the maximum distance in the (Mdims%ndim-1) dimensions, we only calculate it once
            if (Mdims%ndim < 3) then
                l1_max = maxval(X_ALL(1,:))
                l1_min = minval(X_ALL(1,:))
                if (IsParallel()) then
                    call allmax(l1_max)
                    call allmin(l1_min)
                end if
                domain_length = abs(l1_max-l1_min)
            else
                l1_max = maxval(X_ALL(1,:))
                l1_min = minval(X_ALL(1,:))
                l2_max = maxval(X_ALL(2,:))
                l2_min = minval(X_ALL(2,:))
                if (IsParallel()) then
                    call allmax(l1_max)
                    call allmin(l1_min)
                    call allmax(l2_max)
                    call allmin(l2_min)
                end if
                domain_length = max(l1_max-l1_min, l2_max-l2_min)
            end if
            !Domain thickness as the maximum distance in the last dimension
            h_max = maxval(X_ALL(Mdims%ndim,:))
            h_min = minval(X_ALL(Mdims%ndim,:))
            if (IsParallel()) then
                call allmax(h_max)
                call allmin(h_min)
            end if
            domain_thickness = abs(h_max-h_min)
            read_sizes = .true.
        end if
        !Retrieve viscosity field - control-volume wise
        do iphase = 1,  Mdims%n_in_pres
          state_viscosity => extract_tensor_field( state( iphase ), 'Viscosity' )
          viscosity(iphase, :) = state_viscosity%val(1,1,1) !Take the first term only as for porous media we consider only scalar
        end do
        !viscosity => extract_tensor_field(state(Mdims%n_in_pres),"Viscosity" )
        !Retrieve phase density - control-volume wise
        density  => extract_tensor_field( packed_state, "PackedDensity" )
        !Retrieve velocity - control-volume wise
        velocity => extract_tensor_field( packed_state, "PackedVelocity" )
        !Retrieve Immobile fraction - control-volume wise
        call get_var_from_packed_state(packed_state, CV_Immobile_Fraction = Imble_frac)
        !Retrive saturation - control-volume wise
        call get_var_from_packed_state(packed_state, PhaseVolumeFraction = saturation)
        !Retrive Capillary entry pressure - element wise
        call get_var_from_packed_state(packed_state, Cap_entry_pressure = cap_entry_pres)
        !Retrive relative permeability end points- element wise
        call get_var_from_packed_state(packed_state, EndPointRelperm = end_point_relperm)
        !Retrive relative permeability exponents - element wise
        call get_var_from_packed_state(packed_state, RelpermExponent = exponent_relperm)
        !Retrieve relative permeabilty field - element wise
        call get_relative_permeability(Mdims, saturation, Imble_frac, end_point_relperm, exponent_relperm, relperm)
        !Total mobility - control-volume wise
        total_mobility = 0.
        do cv_nodi = 1, Mdims%cv_nonods
            if (.not. node_owned(pressure, cv_nodi)) cycle
            total_mobility(cv_nodi) = sum( relperm(:, cv_nodi) / viscosity( :, cv_nodi) )
        end do
        !P0 Darcy velocity - control-volume wise
        do iphase = 1, Mdims%n_in_pres
            Darcy_velocity(iphase)%ptr => extract_vector_field(state(iphase),"DarcyVelocity")
        end do
        Darcy_velocity_cvwise = 0.
        counter_subcv = 0
        do ele = 1, Mdims%totele
            do u_iloc = 1, Mdims%u_nloc
                u_inod = ndgln%u(( ele - 1 ) * Mdims%u_nloc + u_iloc)
                do cv_iloc = 1, Mdims%cv_nloc
                    cv_nodi = ndgln%cv(( ele - 1) * Mdims%cv_nloc + cv_iloc )
                    if (.not. node_owned(pressure, cv_nodi)) cycle
                    do iphase = 1, Mdims%n_in_pres
                        Darcy_velocity_cvwise(:,cv_nodi) = Darcy_velocity_cvwise(:,cv_nodi) + Darcy_velocity(iphase)%ptr%val(:,u_inod)
                    end do
                    counter_subcv(cv_nodi) = counter_subcv(cv_nodi) + 1
                end do
            end do
        end do
        nDarcy_velocity_cvwise = 0.
        do cv_nodi = 1, Mdims%cv_nonods
          if (.not. node_owned(pressure, cv_nodi)) cycle
          Darcy_velocity_cvwise(:,cv_nodi) = Darcy_velocity_cvwise(:,cv_nodi)/counter_subcv(cv_nodi)
          nDarcy_velocity_cvwise(cv_nodi) = norm2(Darcy_velocity_cvwise(:,cv_nodi))
        end do
        !**************************************!
        !!! ***Calculating all parameters*** !!!
        !**************************************!

        if (.not.read_options) then
            !We read the options just once

            !!! Number of phases !!!
            n_phases = Mdims%n_in_pres
            !!! have wells !!!
            wells = (Mdims%npres > 1)
            !!! Number of components !!!
            n_components = Mdims%ncomp
            !!! gravity !!!
            gravity = have_option("/physical_parameters/gravity")
            !!! Capillary pressure !!!
            cap_pressure = have_option_for_any_phase("/multiphase_properties/capillary_pressure", Mdims%nphase)
            !!! Single-phase flow !!!
            one_phase = (Mdims%n_in_pres == 1)
            !!! Black Oil
            black_oil = have_option( "/physical_parameters/black-oil_PVT_table")
            !!! Over relaxation !!!
            ov_relaxation = have_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation')

            !!! Aspect ratio !!!
            ! R_L = sqrt(permz/permx)*L/H
            aspect_ratio = domain_length*SQRT(average_perm_thickness/average_perm_length)/domain_thickness

            read_options = .true.
        end if

        !!! Courant number and Shockfront courant number !!!
        ! CFL = u*delta_t/delta_x
        !Sometimes the shock-front courant number is not well calculated, then use previous value
        !if (abs(courant_number_in(2)) < 1d-8 ) courant_number_in(2) = backup_shockfront_Courant
        !Courant number
        courant_number = courant_number_in(1)
        !shockfront courant number
        shockfront_courant_number = courant_number_in(2)
        backup_shockfront_Courant = shockfront_courant_number

        !!! Total mobility  !!!
        ! calculate average, max and min values
        average_total_mobility = sum(total_mobility)
        max_total_mobility = maxval(total_mobility)
        min_total_mobility = minval(total_mobility)
        if (IsParallel()) then
            call allsum(average_total_mobility)
            call allmax(max_total_mobility)
            call allmin(min_total_mobility)
        end if
        average_total_mobility = average_total_mobility/total_cv
        !!! Darcy velocity !!!
        ! calculate average, max and min values
        average_Darcy_velocity = sum(nDarcy_velocity_cvwise)
        max_Darcy_velocity = maxval(nDarcy_velocity_cvwise)
        min_Darcy_velocity = minval(nDarcy_velocity_cvwise)
        if (IsParallel()) then
            call allsum(average_Darcy_velocity)
            call allmax(max_Darcy_velocity)
            call allmin(min_Darcy_velocity)
        end if
        average_Darcy_velocity = average_Darcy_velocity/total_cv

        !!! shockfront mobility ratio !!!
        ! M_f = lambda_T(at the front)/lambda_T(ahead of the front)
        min_shockfront_mobratio = 999.
        max_shockfront_mobratio = 0.
        average_shockfront_mobratio = 0.
        shockfront_counter = 0
        if (Mdims%n_in_pres > 1) then
            do ele = 1, Mdims%totele
              !Here we have to go over the halos to avoid problems...
              ! if (.not. element_owned(permeability, ele)) cycle
                if (shock_front_in_ele(ele, Mdims, saturation, ndgln, Imble_frac)) then
                    aux_shockfront_mobratio = shock_front_mobility_ratio(ele, Mdims, saturation, ndgln, Imble_frac, total_mobility)
                    min_shockfront_mobratio = min(min_shockfront_mobratio, aux_shockfront_mobratio)
                    max_shockfront_mobratio = max(max_shockfront_mobratio, aux_shockfront_mobratio)
                    average_shockfront_mobratio = average_shockfront_mobratio + aux_shockfront_mobratio
                    shockfront_counter = shockfront_counter + 1
                end if
            end do
            if (IsParallel()) then
                call allmax(max_shockfront_mobratio)
                call allmin(min_shockfront_mobratio)
                call allsum(average_shockfront_mobratio)
                call allsum(shockfront_counter)
            end if
            average_shockfront_mobratio = average_shockfront_mobratio/max(shockfront_counter,1)
            shockfront_number_ratio = real(shockfront_counter)/total_ele
        end if
        !Verify whether there is a shockfront or not
        if (shockfront_counter == 0) then
            min_shockfront_mobratio = 1.
            max_shockfront_mobratio = 1.
            average_shockfront_mobratio = 1.
            shockfront_number_ratio = 0.0
        end if
        !!! Capillary numbers !!!
        ! N_cv,L = (permx*k^e_rd*P^e_c)/(u*mu_d*L)
        ! N_cv,T = (permz*k^e_rd*P^e_c*L)/(u*mu_d*H**2)
        if (cap_pressure) then
            allocate( longitudinal_capillary(Mdims%cv_nonods) )
            allocate( transverse_capillary(Mdims%cv_nonods) )
            longitudinal_capillary =0;transverse_capillary=0.
            do ele = 1, Mdims%totele
              auxR = end_point_relperm(Mdims%n_in_pres,ele)*sum(cap_entry_pres(:,ele))
              do cv_iloc = 1, Mdims%cv_nloc
                cv_nodi = ndgln%cv((ele-1)*Mdims%cv_nloc+cv_iloc)
                if (.not. node_owned(pressure, cv_nodi)) cycle
                longitudinal_capillary(cv_nodi) = average_perm_length*auxR/ &
                                                & (nDarcy_velocity_cvwise(cv_nodi)*viscosity( Mdims%n_in_pres, cv_nodi)*domain_length)
                transverse_capillary(cv_nodi) = average_perm_thickness*auxR*domain_length / &
                                                & (nDarcy_velocity_cvwise(cv_nodi)*viscosity( Mdims%n_in_pres, cv_nodi)*domain_thickness**2)
              end do
            end do
            ! calculate average, max and min values
            average_longitudinal_capillary = sum(longitudinal_capillary)
            average_transverse_capillary = sum(transverse_capillary)
            max_longitudinal_capillary = maxval(longitudinal_capillary)
            max_transverse_capillary = maxval(transverse_capillary)
            min_longitudinal_capillary = minval(longitudinal_capillary)
            min_transverse_capillary = minval(transverse_capillary)
            if (IsParallel()) then
                call allsum(average_longitudinal_capillary)
                call allsum(average_transverse_capillary)
                call allmax(max_longitudinal_capillary)
                call allmax(max_transverse_capillary)
                call allmin(min_longitudinal_capillary)
                call allmin(min_transverse_capillary)
            end if
            average_longitudinal_capillary = average_longitudinal_capillary/total_cv
            average_transverse_capillary = average_transverse_capillary/total_cv
            deallocate(longitudinal_capillary)
            deallocate(transverse_capillary)
        else
            average_longitudinal_capillary = 0.
            average_transverse_capillary = 0.
            max_longitudinal_capillary = 0.
            max_transverse_capillary = 0.
            min_longitudinal_capillary = 0.
            min_transverse_capillary = 0.
        end if
        !!! gravity numbers !!!
        ! N_bv  = (permx*k^e_rd*delta_rho*g*cos(alpha))/(u*mu_d*L)
        ! N_bv,L = (permx*k^e_rd*delta_rho*g*sin(alpha))/(u*mu_d)
        ! N_bv,T = (permz*k^e_rd*delta_rho*g*cos(alpha)*L)/(u*mu_d*H)
        if (gravity) then
            allocate( buoyancy_number(Mdims%cv_nonods) )
            allocate( longitudinal_buoyancy(Mdims%cv_nonods) )
            allocate( transverse_buoyancy(Mdims%cv_nonods) )
            call get_option( "/physical_parameters/gravity/magnitude", gravity_magnitude)
            gravity_direction => extract_vector_field( state( 1 ), 'GravityDirection' )
            buoyancy_number = 0.; longitudinal_buoyancy =0.; transverse_buoyancy=0.
            do cv_nodi = 1, Mdims%cv_nonods
                if (.not. node_owned(pressure, cv_nodi)) cycle
                delta_density = (maxval(density%val(1,:,cv_nodi)) - minval(density%val(1,:,cv_nodi)))
                sin_alpha =  dot_product(Darcy_velocity_cvwise(:,cv_nodi), gravity_direction%val(:,1))/(norm2(Darcy_velocity_cvwise(:,cv_nodi))*norm2(gravity_direction%val(:,1)))
                buoyancy_number(cv_nodi) = average_perm_length*end_point_relperm(Mdims%n_in_pres,cv_nodi)*delta_density*gravity_magnitude*sqrt(1-sin_alpha**2)*domain_thickness / &
                                            & (nDarcy_velocity_cvwise(cv_nodi)*viscosity( Mdims%n_in_pres, cv_nodi)*domain_length)
                longitudinal_buoyancy(cv_nodi) = average_perm_length*end_point_relperm(Mdims%n_in_pres,cv_nodi)*delta_density*gravity_magnitude*sin_alpha / &
                                               & (nDarcy_velocity_cvwise(cv_nodi)*viscosity( Mdims%n_in_pres, cv_nodi))
                transverse_buoyancy(cv_nodi) = average_perm_thickness*end_point_relperm(Mdims%n_in_pres,cv_nodi)*delta_density*gravity_magnitude*sqrt(1-sin_alpha**2)*domain_length / &
                                                & (nDarcy_velocity_cvwise(cv_nodi)*viscosity( Mdims%n_in_pres, cv_nodi)*domain_thickness)
            end do
            ! calculate average, max and min values
            average_buoyancy_number = sum(buoyancy_number)
            average_longitudinal_buoyancy = sum(longitudinal_buoyancy)
            average_transverse_buoyancy = sum(transverse_buoyancy)
            max_buoyancy_number = maxval(buoyancy_number)
            max_longitudinal_buoyancy = maxval(longitudinal_buoyancy)
            max_transverse_buoyancy = maxval(transverse_buoyancy)
            min_buoyancy_number = minval(buoyancy_number)
            min_longitudinal_buoyancy = minval(longitudinal_buoyancy)
            min_transverse_buoyancy = minval(transverse_buoyancy)
            if (IsParallel()) then
                call allsum(average_buoyancy_number)
                call allsum(average_longitudinal_buoyancy)
                call allsum(average_transverse_buoyancy)
                call allmax(max_buoyancy_number)
                call allmax(max_longitudinal_buoyancy)
                call allmax(max_transverse_buoyancy)
                call allmin(min_buoyancy_number)
                call allmin(min_longitudinal_buoyancy)
                call allmin(min_transverse_buoyancy)
            end if
            average_buoyancy_number = average_buoyancy_number/total_cv
            average_longitudinal_buoyancy = average_longitudinal_buoyancy/total_cv
            average_transverse_buoyancy = average_transverse_buoyancy/total_cv
            deallocate(buoyancy_number)
            deallocate(longitudinal_buoyancy)
            deallocate(transverse_buoyancy)
        else
            average_buoyancy_number = 0.
            average_longitudinal_buoyancy = 0.
            average_transverse_buoyancy = 0.
            max_buoyancy_number = 0.
            max_longitudinal_buoyancy = 0.
            max_transverse_buoyancy = 0.
            min_buoyancy_number = 0.
            min_longitudinal_buoyancy = 0.
            min_transverse_buoyancy = 0.
        end if
        !!! overrelaxation parameter !!!
        if (ov_relaxation) then
            !allocate( overrelaxation(Mdims%cv_nonods) )
            !call getOverrelaxation_parameter(state, packed_state, Mdims, ndgln, overrelaxation, Phase_with_Pc)
            max_overrelaxation = maxval(overrelaxation)
            min_overrelaxation = minval(overrelaxation)
            average_overrelaxation = 0.
            do cv_nodi = 1, Mdims%cv_nonods
                if (.not. node_owned(pressure, cv_nodi)) cycle
                average_overrelaxation = average_overrelaxation +  overrelaxation(cv_nodi)
            end do
            if (IsParallel()) then
                call allmax(max_overrelaxation)
                call allmin(min_overrelaxation)
                call allsum(average_overrelaxation)
            end if
            average_overrelaxation = average_overrelaxation/total_cv
        else
            max_overrelaxation = 0.
            min_overrelaxation = 0.
            average_overrelaxation = 0.
        end if

        !!! Peclet number !!!
        ! Peclet = u * L / Diffusivity;
        ! invPeclet = 1/Peclet
        if (ov_relaxation) then
            allocate( invPeclet(Mdims%cv_nonods) )
            invPeclet = 0.
            if (present_and_true(for_transport)) then
                call get_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation/Vanishing_for_transport', diffusivity)
            else
                call get_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation', diffusivity)
            end if
            diffusivity = abs(diffusivity)
            do cv_nodi = 1, Mdims%cv_nonods
                if (.not. node_owned(pressure, cv_nodi)) cycle
                invPeclet(cv_nodi) = diffusivity/(nDarcy_velocity_cvwise(cv_nodi)*domain_length)
            end do
            max_invPeclet = maxval(invPeclet)
            min_invPeclet = minval(invPeclet)
            average_invPeclet = sum(invPeclet)
            if (IsParallel()) then
                call allmax(max_invPeclet)
                call allmin(min_invPeclet)
                call allsum(average_invPeclet)
            end if
            average_invPeclet = average_invPeclet/total_cv
            deallocate(invPeclet)
        else
            max_invPeclet = 0.
            min_invPeclet = 0.
            average_invPeclet = 0.
        end if
        deallocate(Darcy_velocity_cvwise, total_mobility, counter_subcv, nDarcy_velocity_cvwise)
        !***************************!
        !!! **Call the ML model** !!!
        !***************************!
        block

        real(c_float), dimension(17)  :: raw_input
        real(c_float), pointer        :: out_result(:)
        real, save                    :: target_nli
        backtrack_par_factor = -1. ! Assign -1 to backtrack_par_factor in all cores
        if (getprocno() == 1) then ! To let only 1 core to load and predict using the ML model
            if (.not. loaded_file) then
                call xgboost_load_model()
                call get_option("/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/ML_model_path/target_num_nonlinear", target_nli, default=1.0)
                loaded_file = .true.
            end if

            raw_input = (/aspect_ratio, &
                        & courant_number, &
                        & shockfront_courant_number, &
                        & shockfront_number_ratio, &
                        & average_total_mobility, &
                        & average_Darcy_velocity, &
                        & average_shockfront_mobratio, &
                        & average_longitudinal_capillary, &
                        & average_transverse_capillary, &
                        & average_buoyancy_number, &
                        & average_longitudinal_buoyancy, &
                        & average_transverse_buoyancy, &
                        & average_overrelaxation, &
                        & res, &
                        & resold, &
                        & res/resold, &
                        & target_nli /) !Target number of inner-nonlinear iteration
            call xgboost_predict(raw_input, out_result)
            ! write(*,*) 'XGB model prediction: ',out_result
              backtrack_par_factor = out_result(1)
            nullify(out_result)
        end if
        if (IsParallel()) call allmax(backtrack_par_factor) ! Assign the calculated values of the backtrack_par_factor
        endblock

        !*****************************!
        !!! ***Support functions*** !!!
        !*****************************!

        contains

            subroutine get_relative_permeability(Mdims, saturation, Imble_frac, end_point_relperm, exponent_relperm, relperm)
                !Calculate the relative permeabilities based on Corey correlation
                implicit none
                type(multi_dimensions), intent(in) :: Mdims
                real, dimension(:,:), intent(in) :: saturation
                real, dimension(:,:), intent(in) :: Imble_frac
                real, dimension(:,:), intent(in) :: end_point_relperm
                real, dimension(:,:), intent(in) :: exponent_relperm
                real, dimension(:,:), intent(OUT) :: relperm
                integer :: cv_nodi, iphase, ele, cv_iloc

                do ele = 1, Mdims%totele
                  do cv_iloc = 1, Mdims%cv_nloc
                    cv_nodi = ndgln%cv((ele-1)*Mdims%cv_nloc+cv_iloc)
                      do iphase = 1, Mdims%n_in_pres
                          relperm(iphase,cv_nodi) = end_point_relperm(iphase,ele)*((saturation(iphase,cv_nodi) - Imble_frac(iphase,cv_nodi))/ &
                                                  & (1-sum(Imble_frac(:Mdims%n_in_pres,cv_nodi))))**exponent_relperm(iphase,ele)
                      end do
                  end do
                end do
            end subroutine

            real function shock_front_mobility_ratio(ele, Mdims, sat, ndgln, Imble_frac, total_mobility)
                !Detects whether the element has a shockfront or not and calculates the shockfront mobility ratio
                !if there is not a shock front within the lement the fuction returns -1
                implicit none
                integer, intent(in) :: ele
                type(multi_dimensions), intent(in) :: Mdims
                real, dimension(:,:), intent(in) :: sat
                type(multi_ndgln), intent(in) :: ndgln
                real, dimension(:,:), intent(in) :: Imble_frac
                real, dimension(:), intent(in) :: total_mobility
                !Local variables
                integer :: iphase, cv_iloc, cv_nodi
                real :: minival, maxival, aux_sat, tmob_aheadfront, tmob_atfront
                real, parameter :: tol = 0.05 !Shock fronts smaller than this are unlikely to require extra handling

                !Starts the value
                shock_front_mobility_ratio = -1.0

                minival = 999.; maxival = 0.
                do iphase = 1, mdims%n_in_pres - 1
                    do cv_iloc = 1, Mdims%cv_nloc
                        cv_nodi = ndgln%cv((ele-1)*Mdims%cv_nloc+cv_iloc)
                        aux_sat = sat(iphase, cv_nodi) - Imble_frac(iphase, cv_nodi)
                        if (aux_sat < minival) then
                            minival = aux_sat
                            tmob_aheadfront = total_mobility(cv_nodi)
                        end if
                        if (aux_sat > maxival) then
                            maxival = aux_sat
                            tmob_atfront = total_mobility(cv_nodi)
                        end if
                    end do
                    if  (minival < tol .and. (maxival-minival) > tol) then
                        shock_front_mobility_ratio = tmob_atfront/tmob_aheadfront
                        exit
                    end if
                end do
            end function shock_front_mobility_ratio

            logical function shock_front_in_ele(ele, Mdims, sat, ndgln, Imble_frac)
                !Detects whether the element has a shockfront or not
                implicit none
                integer, intent(in) :: ele !***
                type(multi_dimensions), intent(in) :: Mdims
                real, dimension(:,:), intent(in) :: sat
                type(multi_ndgln), intent(in) :: ndgln
                real, dimension(:,:), intent(in) :: Imble_frac
                !Local variables
                integer :: iphase, cv_iloc, cv_inod
                real :: minival, maxival, aux
                real, parameter :: tol = 0.05 !Shock fronts smaller than this are unlikely to require extra handling

                !Starts the value
                shock_front_in_ele = .false.

                minival = 999.; maxival = 0.
                do iphase = 1, mdims%n_in_pres - 1
                    do cv_iloc = 1, Mdims%cv_nloc
                        cv_inod = ndgln%cv((ele-1)*Mdims%cv_nloc+cv_iloc)
                        aux = sat(iphase, ndgln%cv((ele-1)*Mdims%cv_nloc+cv_iloc)) - Imble_frac(iphase,cv_inod)
                        minival = min(aux, minival)
                        maxival = max(aux, maxival)
                    end do
                    if  (minival < tol .and. (maxival-minival) > tol) then
                        shock_front_in_ele = .true.
                        exit
                    end if
                end do
            end function shock_front_in_ele

    end subroutine AI_backtracking_parameters
#endif





    !> @author Chris Pain, Pablo Salinas, Asiri Obeysekara
    !> @brief Calls to generate for Porous media only the Gradient Matrix, the divergence matrix, the momentum matrix and the mass matrix
    !> Once these matrices (and corresponding RHSs) are generated the system of equations is solved using the projection method.
    !>@param  state Linked list containing all the fields defined in diamond and considered by Fluidity
    !>@param  packed_state Linked list containing all the fields used by IC-FERST, memory partially shared with state
    !!>@param Mdims Data type storing all the dimensions describing the mesh, fields, nodes, etc
    !>@param  CV_GIdims Gauss integration numbers for CV fields
    !>@param FE_GIdims Gauss integration numbers for FE fields
    !>@param  CV_funs Shape functions for the CV mesh
    !>@param  FE_funs Shape functions for the FE mesh
    !>@param  Mspars Sparsity of the matrices
    !>@param  ndgln Global to local variables
    !>@param  Mdisopt Discretisation options
    !>@param  Mmat Matrices for ICFERST
    !>@param  multi_absorp  Absoprtion of associated with the transport field
    !>@param  upwnd Sigmas to compute the fluxes at the interphase for porous media
    !>@param  eles_with_pipe Elements that have a pipe
    !>@param  pipes_aux Information required to define wells
    !>@param  velocity tensor field. Why do we pass it down instead of retrieving it from packed_state... 
    !>@param  pressure tensor fieldWhy do we pass it down instead of retrieving it from packed_state... 
    !>@param  DT Time step size
    !>@param  SUF_SIG_DIAGTEN_BC Like upwnd but for the boundary
    !>@param V_SOURCE Source term
    !>@param  VOLFRA_PORE     Porosity field (Mdims%npres,Mdims%totele)
    !>@param   IGOT_THETA_FLUX, THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J ????
    !>@param  calculate_mass_delta This is to compute mass conservation 
    !>@param  outfluxes variable containing the outfluxes information
    !>@param pres_its_taken Integer containing the number of iterations performed by the solver to solve for pressure
    !>@param nonlinear_its current non-linear iteration
    !>@param Courant_number Global courant number and shock front courant number
SUBROUTINE POROUS_FORCE_BAL_CTY_ASSEM_SOLVE( state, packed_state,  &
    Mdims, CV_GIdims, FE_GIdims, CV_funs, FE_funs, Mspars, ndgln, Mdisopt,  &
    Mmat, multi_absorp, upwnd, eles_with_pipe, pipes_aux, velocity, pressure, &
    DT, SUF_SIG_DIAGTEN_BC, V_SOURCE, VOLFRA_PORE, &
    IGOT_THETA_FLUX, THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J,&
    calculate_mass_delta, outfluxes, pres_its_taken, nonlinear_its, Courant_number)
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
    type(multi_absorption), intent(inout) :: multi_absorp
    type (porous_adv_coefs), intent(inout) :: upwnd
    type(pipe_coords), dimension(:), intent(in):: eles_with_pipe
    type (multi_pipe_package), intent(in) :: pipes_aux
    type( tensor_field ), intent(inout) :: velocity
    type( tensor_field ), intent(inout) :: pressure
    INTEGER, intent( in ) :: IGOT_THETA_FLUX, nonlinear_its
    REAL, DIMENSION(  : , :  ), intent( in ) :: SUF_SIG_DIAGTEN_BC
    REAL, intent( in ) :: DT
    REAL, DIMENSION(  :, :  ), intent( in ) :: V_SOURCE
    REAL, DIMENSION(  :, :  ), intent( in ) :: VOLFRA_PORE
    REAL, DIMENSION( : ,  :  ), intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
    type (multi_outfluxes), intent(inout) :: outfluxes
    real, dimension(:,:), intent(inout) :: calculate_mass_delta
    integer, intent(inout) :: pres_its_taken
    real, dimension(:), intent(inout) :: Courant_number
    ! Local Variables
    character(len=option_path_len) :: solver_option_pressure = "/solver_options/Linear_solver"

    ! If IGOT_CMC_PRECON=1 use a sym matrix as pressure preconditioner,=0 else CMC as preconditioner as well.
    INTEGER :: IGOT_CMC_PRECON
    REAL, DIMENSION( : ), allocatable :: MASS_MN_PRES, MASS_SUF, MASS_CV, UP
    REAL, DIMENSION( :, : ), allocatable :: DIAG_SCALE_PRES
    real, dimension(:,:,:), allocatable :: velocity_absorption, U_SOURCE_CV_ALL
    real, dimension(1,1,1,1) :: UDIFFUSION_ALL
    type( multi_field ) :: UDIFFUSION_VOL_ALL, U_SOURCE_ALL   ! NEED TO ALLOCATE THESE - SUBS TO DO THIS ARE MISSING... - SO SET 0.0 FOR NOW
    REAL, DIMENSION( :, :, : ), allocatable :: DIAG_SCALE_PRES_COUP, INV_B, CMC_PRECON
    REAL, DIMENSION( :, :, : ), allocatable :: DU_VEL
    INTEGER :: CV_NOD, COUNT, CV_JNOD, IPHASE, JPHASE, ndpset, i
    LOGICAL ::  diag, JUST_BL_DIAG_MAT = .false.
    INTEGER :: python_stat
    !Re-scale parameter can be re-used
    real, save :: rescaleVal = -1.0
    !CMC using petsc format
    type(petsc_csr_matrix)::  CMC_petsc
    !TEMPORARY VARIABLES, ADAPT FROM OLD VARIABLES TO NEW
    INTEGER :: IPRES, JPRES, iphase_real, jphase_real
    REAL, DIMENSION( :, : ), allocatable :: UDEN_ALL, UDENOLD_ALL
    REAL, DIMENSION( :, : ), allocatable :: sigma
    REAL, DIMENSION( :, : ), pointer :: DEN_ALL, DENOLD_ALL
    type( tensor_field ), pointer :: OLDvelocity, den_all2, denold_all2, tfield, den_all3!, u_all2
    type( tensor_field ), pointer :: p_all, cvp_all, deriv, python_tfield
    type( vector_field ), pointer :: x_all2, U
    type( scalar_field ), pointer :: gamma
    type( vector_field ) :: packed_vel, rhs
    type( vector_field ) :: deltap, rhs_p
    type(tensor_field) :: cdp_tensor
    type( csr_sparsity ), pointer :: sparsity
    INTEGER, DIMENSION ( 1, Mdims%npres, surface_element_count(pressure) ) :: WIC_P_BC_ALL
    type( tensor_field ) :: pressure_BCs
    integer :: idim, idx1, idx2, ndim
    !Variables to control de performance of the solvers
    integer :: its_taken
    integer, save :: max_allowed_P_its = -1, max_allowed_V_its = -1
    real, dimension(Mdims%totele) :: MASS_ELE
    REAL, DIMENSION ( :, :, :,:, :, :, :), allocatable :: DIAG_BIGM_CON
    REAL, DIMENSION ( :, :, :,:, :, :, :), allocatable :: BIGM_CON

    !Since we save the parameter rescaleVal, we only do this one time
    if (rescaleVal < 0.) then
        tfield => extract_tensor_field(packed_state,"Permeability")
        rescaleVal = minval(tfield%val, MASK = tfield%val > 1d-30)
        !If it is parallel then we want to be consistent between cpus
        if (IsParallel()) call allmin(rescaleVal)
    end if

    !Retrieve solver setting configurations
    solver_option_pressure = "/solver_options/Linear_solver"
    if (have_option('/solver_options/Linear_solver/Custom_solver_configuration/Pressure')) then
      solver_option_pressure = '/solver_options/Linear_solver/Custom_solver_configuration/Pressure'
    end if
    if(max_allowed_P_its < 0)  then
        call get_option( trim(solver_option_pressure)//'/max_iterations',&
         max_allowed_P_its, default = 10000)
    end if

    deriv => extract_tensor_field( packed_state, "PackedDRhoDPressure" )
    IGOT_CMC_PRECON = 0
    !sprint_to_do!this looks like a place than can be easily optimized
    ALLOCATE( UDEN_ALL( Mdims%nphase, Mdims%cv_nonods ), UDENOLD_ALL( Mdims%nphase, Mdims%cv_nonods ) )
    UDEN_ALL = 0.; UDENOLD_ALL = 0.
    ewrite(3,*) 'In POROUS_FORCE_BAL_CTY_ASSEM_SOLVE'
    ALLOCATE( Mmat%CT( Mdims%ndim, Mdims%nphase, Mspars%CT%ncol )) ; Mmat%CT=0.
    call allocate(Mmat%CT_RHS,Mdims%npres,pressure%mesh,"Mmat%CT_RHS")
    ALLOCATE( Mmat%U_RHS( Mdims%ndim, Mdims%nphase, Mdims%u_nonods )) ; !initialised inside the subroutines
    ALLOCATE( DIAG_SCALE_PRES( Mdims%npres,Mdims%cv_nonods )) ; DIAG_SCALE_PRES=0.
    ALLOCATE(DIAG_SCALE_PRES_COUP(Mdims%npres,Mdims%npres,Mdims%cv_nonods),INV_B(Mdims%nphase,Mdims%nphase,Mdims%cv_nonods))
    ALLOCATE( CMC_PRECON( Mdims%npres, Mdims%npres, Mspars%CMC%ncol*IGOT_CMC_PRECON)) ; IF(IGOT_CMC_PRECON.NE.0) CMC_PRECON=0.
    ALLOCATE( MASS_MN_PRES( Mspars%CMC%ncol )) ; MASS_MN_PRES=0.
    ALLOCATE( MASS_CV( Mdims%cv_nonods )) ; MASS_CV=0.

    OLDvelocity => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldVelocity" )
    X_ALL2 => EXTRACT_VECTOR_FIELD( PACKED_STATE, "PressureCoordinate" )
    P_ALL => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedFEPressure" )
    !For porous media we do not need PackedCVPressure
    CVP_ALL => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedFEPressure" )

    DEN_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedDensity" )
    DENOLD_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldDensity" )

    DEN_ALL(1:, 1:) => DEN_ALL2%VAL( 1, :, : )
    DENOLD_ALL(1:, 1:) => DENOLD_ALL2%VAL( 1, :, : )

    call allocate(deltaP,Mdims%npres,pressure%mesh,"DeltaP")
    call allocate(rhs_p,Mdims%npres,pressure%mesh,"PressureCorrectionRHS")

    !Calculate gravity source terms
    allocate(U_SOURCE_CV_ALL(Mdims%ndim, Mdims%nphase, Mdims%cv_nonods))
    U_SOURCE_CV_ALL=0.0
    call calculate_u_source_cv( Mdims, state, packed_state, DEN_ALL, U_SOURCE_CV_ALL )

    ! define velocity_absorption here...
    ! update velocity absorption
    allocate(velocity_absorption(Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase, Mdims%mat_nonods))
    velocity_absorption=0.0 
    ! Check for a python-set absorption field
    ! Assumes that python blocks are (nphase x nphase) and isotropic
    python_tfield => extract_tensor_field( state(1), "UAbsorB", python_stat )
    if (python_stat==0) then
       ewrite(3,*)"Python UAbsorB"
       velocity_absorption = 0.0
       do iphase = 1, Mdims%nphase
          do jphase = 1, Mdims%nphase
             do idim = 1, Mdims%ndim
                idx1 = idim+(iphase-1)*Mdims%ndim ; idx2 = idim+(jphase-1)*Mdims%ndim
                velocity_absorption( idx1, idx2, : ) = python_tfield%val( iphase, jphase, : )
             end do
          end do
       end do
    end if
    ! update velocity source
    if (have_option_for_any_phase('vector_field::Velocity' // &
      '/prognostic/vector_field::Source', Mdims%n_in_pres)) then
      call allocate_multi_field( Mdims, u_source_all, Mdims%u_nonods, "SourceTerm")
      call update_velocity_source( state, Mdims, u_source_all )
    end if

!Temporary conversion
if (associated(multi_absorp%PorousMedia%val))then!sprint_to_do AVOID THESE CONVERSIONS...
do cv_nod = 1, size(multi_absorp%PorousMedia%val,4)
    call add_multi_field_to_array(multi_absorp%PorousMedia, velocity_absorption(:,:,cv_nod), 1, 1, cv_nod, 1.0)
end do
end if

        !Allocation of storable matrices
        if (.not.Mmat%Stored) then
            allocate(Mmat%C_CV(Mdims%ndim, Mdims%nphase, Mspars%C%ncol)); Mmat%C_CV = 0.
            if (Mmat%compact_PIVIT_MAT) then!Use a compacted and lumped version of the mass matrix
                    !sprint_to_do for this to work with wells we need to change the sparsity, but that still needs to be done!
                allocate( Mmat%PIVIT_MAT( 1, 1, Mdims%totele ) ); Mmat%PIVIT_MAT=0.0
            else
                allocate( Mmat%PIVIT_MAT( Mdims%ndim * Mdims%nphase * Mdims%u_nloc, Mdims%ndim * Mdims%nphase * Mdims%u_nloc, Mdims%totele ) ); Mmat%PIVIT_MAT=0.0
            end if
        end if

    if( have_option_for_any_phase( '/multiphase_properties/capillary_pressure', Mdims%nphase ) )then
        call calculate_capillary_pressure(packed_state, ndgln, Mdims%totele, Mdims%cv_nloc, CV_funs)
    end if

    ALLOCATE( MASS_SUF( 1 )) ; MASS_SUF=0.

    if ( Mdims%npres > 1 ) then
       call get_entire_boundary_condition( pressure, ['weakdirichlet'],&
            pressure_BCs, WIC_P_BC_ALL )
       !Array defined to check where to apply pressure or velocity BCs (for wells)
       if (.not. associated(Mmat%WIC_FLIP_P_VEL_BCS)) then
         allocate(Mmat%WIC_FLIP_P_VEL_BCS( 1,Mdims%npres,surface_element_count(pressure)))
       end if
       !Assing values to check later in multi_pipes when applying the BCs
       call flip_p_and_v_bcs(Mdims, WIC_P_BC_ALL, pressure_BCs, Mmat%WIC_FLIP_P_VEL_BCS)
    end if

    CALL CV_ASSEMB_FORCE_CTY( state, packed_state, &
        Mdims, CV_GIdims, FE_GIdims, CV_funs, FE_funs, Mspars, ndgln, Mdisopt, Mmat,upwnd, &
         velocity, pressure, multi_absorp, eles_with_pipe, pipes_aux,&
        X_ALL2%VAL, velocity_absorption, U_SOURCE_ALL, U_SOURCE_CV_ALL, &
        velocity%VAL, OLDvelocity%VAL, &
        CVP_ALL%VAL, DEN_ALL, DENOLD_ALL, DERIV%val(1,:,:), &
        DT, MASS_MN_PRES, MASS_ELE,& ! pressure matrix for projection method
        .false.,  MASS_SUF, SUF_SIG_DIAGTEN_BC, &
        V_SOURCE, VOLFRA_PORE, Courant_number, &
        DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
        JUST_BL_DIAG_MAT, UDEN_ALL, UDENOLD_ALL, UDIFFUSION_ALL,  UDIFFUSION_VOL_ALL, &
        IGOT_THETA_FLUX, THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, &
        0,.false., calculate_mass_delta, outfluxes, DIAG_BIGM_CON, BIGM_CON ) !

    !If pressure in CV then point the FE matrix Mmat%C to Mmat%C_CV
    if ( Mmat%CV_pressure ) Mmat%C => Mmat%C_CV

    !###################For wells###########################
    if ( Mdims%npres > 1 ) then
       ALLOCATE( SIGMA( Mdims%nphase, Mdims%mat_nonods ) )
       DO IPHASE = 1, Mdims%nphase
          SIGMA( IPHASE, : ) = velocity_absorption( (IPHASE-1)*Mdims%ndim+1, (IPHASE-1)*Mdims%ndim+1, : )
       END DO

       !Introduce well modelling
       CALL MOD_1D_FORCE_BAL_C( STATE, packed_state, Mdims, Mspars, Mmat, ndgln, eles_with_pipe,&
            associated(Mmat%PIVIT_MAT) .and. .not.Mmat%Stored, WIC_P_BC_ALL, pressure_BCs%val, SIGMA,&
            velocity%VAL, U_SOURCE_ALL, U_SOURCE_CV_ALL, pipes_aux )
       call deallocate( pressure_BCs )
       DEALLOCATE( SIGMA )
    end if

    deallocate(velocity_absorption, U_SOURCE_CV_ALL)
    if (u_source_all%have_field) call deallocate_multi_field(U_SOURCE_ALL, .true.)
    !Now invert the Mass matrix
    if (.not.Mmat%Stored) then
        CALL Mass_matrix_inversion(Mmat%PIVIT_MAT, Mdims )
      end if
    ! solve using a projection method
    call allocate(cdp_tensor,velocity%mesh,"CDP",dim = velocity%dim); call zero(cdp_tensor)
    ! Put pressure in rhs of force balance eqn: CDP = Mmat%C * P
    call C_MULT2_MULTI_PRES(Mdims, Mspars, Mmat, P_ALL%val, CDP_tensor)

    !"########################UPDATE VELOCITY STEP####################################"
    !For porous media we calculate the velocity as M^-1 * CDP, no solver is needed
    CALL Mass_matrix_MATVEC( velocity % VAL, Mmat%PIVIT_MAT, Mmat%U_RHS + CDP_tensor%val, Mdims%ndim, Mdims%nphase, &
        Mdims%totele, Mdims%u_nloc, ndgln%u )
    !"########################UPDATE PRESSURE STEP####################################"
    !Form pressure matrix (Sprint_to_do move this (and the allocate!) just before the pressure solver
    sparsity=>extract_csr_sparsity(packed_state,'CMCSparsity')
    diag = Mdims%npres == 1!Make it non-diagonal to allow coupling between reservoir and pipes domains
    call allocate(CMC_petsc,sparsity,[Mdims%npres,Mdims%npres],"CMC_petsc",diag)
    CALL COLOR_GET_CMC_PHA( Mdims, Mspars, ndgln, Mmat,&
    DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
    CMC_petsc, CMC_PRECON, IGOT_CMC_PRECON, MASS_MN_PRES, &
    pipes_aux, .false.,  MASS_SUF, .false. )
! call MatView(CMC_petsc%M,   PETSC_VIEWER_STDOUT_SELF, ipres)

    !This section is to impose a pressure of zero at some node (when solving for only a gradient of pressure)
    !This is deprecated as the remove_null space method works much better (provided by PETSc)
    call get_option( '/material_phase[0]/scalar_field::Pressure/' // 'prognostic/reference_node', ndpset, default = 0 )
    if ( ndpset /= 0 ) rhs_p%val( 1, ndpset ) = 0.0
    !======================================================================================
    ! solve for pressure correction DP that is solve CMC*DP=P_RHS...
    ewrite(3,*)'about to solve for pressure'
    !Perform Div * U for the RHS of the pressure equation
    rhs_p%val = 0.
    call compute_DIV_U(Mdims, Mmat, Mspars, velocity%val, INV_B, rhs_p)
    rhs_p%val = -rhs_p%val + Mmat%CT_RHS%val
    call include_wells_and_compressibility_into_RHS(Mdims, rhs_p, DIAG_SCALE_PRES, MASS_MN_PRES, MASS_SUF, pipes_aux, DIAG_SCALE_PRES_COUP)
    !Call impose strong BCs must be here before the rescaling and after include wells and
    !compressibility - only done if using P0DG and gamma is zero at the BC see
    !subroutine initialize_pipes_package_and_gamma for more information
    if (is_P0DGP1 .and. Mdims%npres > 1) call impose_strong_bcs_wells(state, pipes_aux, Mdims, Mmat, ndgln, CMC_petsc,pressure,rhs_p%val)
    !Re-scale system so we can deal with SI units of permeability
    call scale(cmc_petsc, 1.0/rescaleVal)
    rhs_p%val = rhs_p%val / rescaleVal

    !"########################UPDATE PRESSURE STEP####################################"
    call solve_and_update_pressure(Mdims, rhs_p, P_all%val, deltap, cmc_petsc)
    call deallocate(cmc_petsc);call deallocate(rhs_p)
    if (isParallel()) call halo_update(P_all)

    !######################## CORRECTION VELOCITY STEP####################################
    !Ensure that the velocity fulfils the continuity equation before moving on
    call project_velocity_to_affine_space(Mdims, Mmat, Mspars, ndgln, velocity, deltap, cdp_tensor)
    call deallocate(deltaP)
    if (isParallel()) call halo_update(velocity)
    call DEALLOCATE( CDP_tensor )
    !######################## CORRECTION VELOCITY STEP####################################
!
    DEALLOCATE( Mmat%CT )
    DEALLOCATE( DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B )
    DEALLOCATE( Mmat%U_RHS )
    DEALLOCATE( MASS_MN_PRES )
    call deallocate(Mmat%CT_RHS)

    if (associated(UDIFFUSION_VOL_ALL%val)) call deallocate_multi_field(UDIFFUSION_VOL_ALL)

    ewrite(3,*) 'Leaving POROUS_FORCE_BAL_CTY_ASSEM_SOLVE'
    return
  contains

    !---------------------------------------------------------------------------
    !> @author Geraldine Regnier
    !> @brief Flip P and V BCs for wells.
    !>Implemented for ATES applications or other applications when  needing to
    !>switch on and off wells. For this to work, need to define boundary conditions
    !>for wells for both pressure and velocity in a time-dependent manner through
    !>Diamond using a python function. When the user wants the pressure BC to be
    !>ignored and the velocity BC to be applied, specify negative pressure.
    !>When the user wants the pressure BC to applied, specify normal positive
    !>pressure. Specify velocity BCs for both wells as a constant.
    !---------------------------------------------------------------------------
    subroutine flip_p_and_v_bcs(Mdims, WIC_P_BC_ALL, pressure_BCs, WIC_FLIP_P_VEL_BCS)
      implicit none
      type(multi_dimensions), intent(in) :: Mdims
      integer :: SELE, CV_SILOC
      INTEGER, DIMENSION( :,:,: ), INTENT( IN ) :: WIC_P_BC_ALL
      INTEGER, PARAMETER :: WIC_P_BC_DIRICHLET = 1
      type(tensor_field) :: pressure_BCs
      INTEGER, DIMENSION( 1,Mdims%npres,surface_element_count(pressure)), INTENT(OUT):: WIC_FLIP_P_VEL_BCS

      WIC_FLIP_P_VEL_BCS = 0!There is no flipping
      DO SELE = 1, Mdims%stotel
        IF ( WIC_P_BC_ALL( 1, 2, SELE ) == WIC_P_BC_DIRICHLET ) THEN
          DO CV_SILOC = 1, Mdims%cv_snloc
            if (pressure_BCs%val(1, 2, (SELE-1)*Mdims%cv_snloc + CV_SILOC) < -1) then
              WIC_FLIP_P_VEL_BCS(1, 2, SELE) = 1!Impose velocity BCs
            else !Impose Pressure BCs
              WIC_FLIP_P_VEL_BCS(1, 2, SELE) = 10
            end if
          END DO
        END IF
      END DO

    end subroutine flip_p_and_v_bcs

    !---------------------------------------------------------------------------
    !> @author Pablo Salinas
    !> @brief Compute deltaP by solving the pressure equation using the CMC matrix
    !---------------------------------------------------------------------------
    subroutine solve_and_update_pressure(Mdims, rhs_p, P_all, deltap, cmc_petsc)

      implicit none
      type(multi_dimensions), intent(in) :: Mdims
      type( vector_field ), intent(inout) :: rhs_p
      type( vector_field ), intent(inout) :: deltap
      real, dimension(Mdims%npres, Mdims%cv_nonods), intent(inout) :: P_all!Ensure dynamic conversion from three entries to two
      type(petsc_csr_matrix), intent(inout) ::  CMC_petsc
      !Local variables
      integer :: its_taken

      !Rescale RHS (it is given the the matrix has been already re-scaled)
      call petsc_solve(deltap, cmc_petsc, rhs_p, option_path = trim(solver_option_pressure), iterations_taken = its_taken)
      pres_its_taken = its_taken

      if (its_taken >= max_allowed_P_its) solver_not_converged = .true.
        !Now update the pressure
      P_all = P_all + deltap%val

    end subroutine solve_and_update_pressure

    !---------------------------------------------------------------------------
    !> @author Pablo Salinas
    !> @brief Project back the velocity from a non divergent-free space to a divergent free space
    !---------------------------------------------------------------------------
    subroutine project_velocity_to_affine_space(Mdims, Mmat, Mspars, ndgln, velocity, deltap, cdp_tensor)
      !Project back the velocity from a non divergent-free space to a divergent free space
      implicit none
      type(multi_dimensions), intent(in) :: Mdims
      type (multi_sparsities), intent(in) :: Mspars
      type (multi_matrices), intent(in) :: Mmat
      type(multi_ndgln), intent(in) :: ndgln
      type(tensor_field), intent(inout) :: velocity
      type( vector_field ), intent(inout) :: deltap
      type(tensor_field), intent(inout) :: cdp_tensor!>To reuse memory
      !Local variables
      REAL, DIMENSION( Mdims%ndim,  Mdims%nphase, Mdims%u_nonods ) :: DU_VEL
      integer :: u_inod, iphase, idim
      !Perform C * DP
      call zero(cdp_tensor)
      call C_MULT2_MULTI_PRES(Mdims, Mspars, Mmat, deltap%val, CDP_tensor)

      ! DU = BLOCK_MAT * CDP: Project back the velocity from a non divergent-free space to a divergent free space
      DU_VEL = 0.
      CALL Mass_matrix_MATVEC( DU_VEL, Mmat%PIVIT_MAT, CDP_tensor%val, Mdims%ndim, Mdims%nphase, &
      Mdims%totele, Mdims%u_nloc, ndgln%u )
      !Apply correction to velocity
      velocity%val = velocity%val + DU_VEL

    end subroutine project_velocity_to_affine_space

    !---------------------------------------------------------------------------
    !> @author Pablo Salinas
    !> @brief Calculates the divergence of the velocity by multipliying it by the Ct matrix
    !---------------------------------------------------------------------------
  subroutine compute_DIV_U(Mdims, Mmat, Mspars, velocity, INV_B, rhs_p)
    implicit none
    type(multi_dimensions), intent(in) :: Mdims
    type (multi_sparsities), intent(in) :: Mspars
    type (multi_matrices), intent(in) :: Mmat
    REAL, DIMENSION( :, :, : ), intent(in) :: INV_B, velocity
    type( vector_field ), intent(inout) :: rhs_p
    !Local variables
    integer :: iphase, CV_NOD, ipres

    DO IPRES = 1, Mdims%npres
        CALL CT_MULT2( rhs_p%val(IPRES,:), velocity( :, 1+(IPRES-1)*Mdims%n_in_pres : IPRES*Mdims%n_in_pres, : ), &
        Mdims%cv_nonods, Mdims%u_nonods, Mdims%ndim, Mdims%n_in_pres, &
        Mmat%CT( :, 1+(IPRES-1)*Mdims%n_in_pres : IPRES*Mdims%n_in_pres, : ), Mspars%CT%ncol, Mspars%CT%fin, Mspars%CT%col )
    END DO

  end subroutine compute_DIV_U

  !---------------------------------------------------------------------------
  !> @author Chris Pain, Pablo Salinas
  !> @brief Include in the pressure matrix the compressibility terms (based on taylor expansion series) to ensure that we account for this term
  !> as implicitly as possible. For wells it introduces the coupling between pressures
  !---------------------------------------------------------------------------
  subroutine include_wells_and_compressibility_into_RHS(Mdims, rhs_p, DIAG_SCALE_PRES, MASS_MN_PRES, MASS_SUF, pipes_aux, DIAG_SCALE_PRES_COUP)

    implicit none
    type(multi_dimensions), intent(in) :: Mdims
    type( vector_field ), intent(inout) :: rhs_p
    REAL, DIMENSION( : ), intent(in) :: MASS_MN_PRES, MASS_SUF
    REAL, DIMENSION( :, : ), intent(in) :: DIAG_SCALE_PRES
    REAL, DIMENSION( :, :, : ), intent(in) :: DIAG_SCALE_PRES_COUP
    type (multi_pipe_package), intent(in) :: pipes_aux
    !Local variables
    type( tensor_field ), pointer :: pold_all
    integer :: CV_NOD, COUNT, CV_JNOD, IPRES, JPRES
    logical :: incl_compress_press

    incl_compress_press = .not. (has_boussinesq_aprox .or. .not. have_option_for_any_phase('/phase_properties/Density/compressible', Mdims%nphase))
    !Check whether we have to do something here or not
    if (Mdims%npres == 1 .and. .not. incl_compress_press) return

    ! Matrix vector involving the mass diagonal term
    DO CV_NOD = 1, Mdims%cv_nonods
        DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
            CV_JNOD = Mspars%CMC%col( COUNT )
            DO IPRES = 1, Mdims%npres
                if (incl_compress_press) then
                  IF (( Mdims%npres > 1 )) THEN
                      IF(IPRES==1) THEN
                          rhs_p%val( IPRES, CV_NOD ) = rhs_p%val( IPRES, CV_NOD ) &
                              -DIAG_SCALE_PRES( IPRES, CV_NOD ) * MASS_MN_PRES( COUNT ) * P_ALL%VAL( 1, IPRES, CV_JNOD )
                      ELSE
                          rhs_p%val( IPRES, CV_NOD ) = rhs_p%val( IPRES, CV_NOD ) &
                              -DIAG_SCALE_PRES( IPRES, CV_NOD ) * pipes_aux%MASS_CVFEM2PIPE_TRUE( COUNT ) * P_ALL%VAL( 1, IPRES, CV_JNOD )
                      ENDIF
                  ELSE
                      rhs_p%val( IPRES, CV_NOD ) = rhs_p%val( IPRES, CV_NOD ) &
                          -DIAG_SCALE_PRES( IPRES, CV_NOD ) * MASS_MN_PRES( COUNT ) * P_ALL%VAL( 1, IPRES, CV_JNOD )
                  ENDIF
                end if
                IF ( Mdims%npres > 1 ) THEN
                    DO JPRES = 1, Mdims%npres
                            rhs_p%val( IPRES, CV_NOD ) = rhs_p%val( IPRES, CV_NOD ) &
                                -DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_NOD ) * pipes_aux%MASS_CVFEM2PIPE( COUNT ) * P_ALL%VAL( 1, JPRES, CV_JNOD )
                    END DO
                END IF
            END DO
        END DO
    END DO

    ! solve for pressure correction DP that is solve CMC*DP=P_RHS...
  end subroutine

 END SUBROUTINE POROUS_FORCE_BAL_CTY_ASSEM_SOLVE



 end module multiphase_1D_engine
