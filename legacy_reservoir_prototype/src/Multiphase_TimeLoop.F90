
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

   module multiphase_time_loop

    use field_options
    use write_state_module
    use diagnostic_variables
    use diagnostic_fields_wrapper
    use diagnostic_fields_new, only : &
         calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
         check_diagnostic_dependencies
    use global_parameters, only: timestep, simulation_start_time, simulation_start_cpu_time, &
                               simulation_start_wall_time, &
                               topology_mesh_name, current_time, is_porous_media
    use fldebug
    use reference_counting
    use state_module
    use fields
    use field_options
    use fields_allocates
    use spud
    use signal_vars
    use populate_state_module
    use vector_tools
    use global_parameters
    use memory_diagnostics

!!$ Modules required by adaptivity
    use qmesh_module
    use adapt_state_module
    use adapt_state_prescribed_module!, only: do_adapt_state_prescribed, adapt_state_prescribed
    use populate_sub_state_module
    use fluids_module!, only: pre_adapt_tasks, update_state_post_adapt
    use parallel_tools
    use checkpoint
    use boundary_conditions
    use boundary_conditions_from_options

!!$ Modules indigenous to the prototype Code
    use shape_functions_Linear_Quadratic
    use spact
    use Compositional_Terms
    use multiphase_EOS
    use multiphase_fractures
    use multiphase_caching, only: set_caching_level, cache_level
    use shape_functions_Linear_Quadratic
    use Compositional_Terms
    use Copy_Outof_State
    use cv_advection, only : cv_count_faces
    use multiphase_1D_engine

    use multiphase_fractures
    use boundary_conditions_from_options
    USE multiphase_rheology
    use vtk_interfaces



#ifdef HAVE_ZOLTAN
  use zoltan
#endif
    !use matrix_operations
    !use shape_functions
    !use printout

    implicit none
    private
    public :: MultiFluids_SolveTimeLoop, rheology, dump_outflux


    type(rheology_type), dimension(:), allocatable :: rheology

  contains

    subroutine MultiFluids_SolveTimeLoop( state, &
         dt, nonlinear_iterations, dump_no )
      implicit none
      type( state_type ), dimension( : ), intent( inout ) :: state
      integer, intent( inout ) :: dump_no, nonlinear_iterations
      real, intent( inout ) :: dt

!!$ additional state variables for multiphase & multicomponent

      type(state_type) :: packed_state
      type(state_type), dimension(:), pointer :: multiphase_state, multicomponent_state


!!$ Primary scalars
      integer :: nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, ph_nloc, ph_nonods

!!$ Node global numbers
      integer, dimension( : ), pointer :: x_ndgln_p1, x_ndgln, cv_ndgln, p_ndgln, &
           mat_ndgln, u_ndgln, xu_ndgln, cv_sndgln, p_sndgln, u_sndgln

!!$ Sparsity patterns
      integer :: nlenmcy, mx_nface_p1, mx_ncolacv, mxnele, mx_ncoldgm_pha, &
           mx_ncolmcy, mx_nct, mx_nc, mx_ncolcmc, mx_ncolm, mx_ncolph, &
           ncolacv, ncolmcy, ncolele, ncoldgm_pha, ncolct, ncolc, ncolcmc, ncolm, ncolph
      integer, dimension( : ), allocatable :: finacv, midacv, finmcy,  midmcy, &
           findgm_pha, middgm_pha, findct, &
           findc, findcmc, midcmc, findm, &
           midm
      integer, dimension(:), pointer :: colacv, colmcy, colct,colm,colc,colcmc,coldgm_pha
      integer, dimension(:), pointer :: finele, colele, midele, findph, colph
      integer, dimension(:), pointer :: small_finacv, small_colacv, small_midacv

!!$ Defining element-pair type and discretisation options and coefficients
      integer :: cv_ele_type, p_ele_type, u_ele_type, mat_ele_type, u_sele_type, cv_sele_type, &
           t_disopt, v_disopt, t_dg_vel_int_opt, u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
           comp_diffusion_opt, ncomp_diff_coef, in_ele_upwind, dg_ele_upwind, &
           nits_flux_lim_t, nits_flux_lim_volfra, nits_flux_lim_comp,  IDIVID_BY_VOL_FRAC
      logical :: volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux, &
           t_use_theta_flux, t_get_theta_flux, scale_momentum_by_volume_fraction, q_scheme
      real :: t_beta, v_beta, t_theta, v_theta, u_theta
      real, dimension(:,:,:,:), allocatable, target :: opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new

!!$ Defining time- and nonlinear interations-loops variables
      integer :: itime, dump_period_in_timesteps, final_timestep, &
           NonLinearIteration, NonLinearIteration_Components, dtime
      real :: acctim, finish_time

!!$ Defining problem that will be solved
      logical :: have_temperature_field, have_component_field, have_extra_DiffusionLikeTerm, &
           solve_force_balance, solve_PhaseVolumeFraction

!!$ Defining solver options
      integer :: velocity_max_iterations, PhaseVolumeFraction_max_iterations

!!$ Shape function related fields:
      integer :: cv_ngi, cv_ngi_short, scvngi_theta, sbcvngi, nface, igot_t2, igot_theta_flux, IGOT_THERM_VIS

!!$ For output:
      real, dimension( : ), allocatable :: &
           Mean_Pore_CV

!!$ Variables used in the diffusion-like term: capilarity and surface tension:
      integer :: iplike_grad_sou
      real, dimension( : ), allocatable :: plike_grad_sou_grad, plike_grad_sou_coef

!!$ Adaptivity related fields and options:
      type( tensor_field ) :: metric_tensor
      type( state_type ), dimension( : ), pointer :: sub_state => null()
      integer :: nonlinear_iterations_adapt
      logical :: do_reallocate_fields, not_to_move_det_yet = .false., initialised

!!$ Working arrays:
      real, dimension(:), pointer :: mass_ele

      real, dimension( :, :, :, : ), allocatable :: THERM_U_DIFFUSION
      real, dimension( :, : ), allocatable :: THERM_U_DIFFUSION_VOL

      real, dimension( :, : ), pointer :: THETA_GDIFF

      real, dimension( :, : ), pointer ::  DRhoDPressure, FEM_VOL_FRAC
!!$
      real, dimension( :, : ), pointer :: &
           ScalarField_Source, ScalarField_Source_Store, ScalarField_Source_Component
      real, dimension( :, :, : ), pointer :: Velocity_U_Source, Velocity_U_Source_CV
      real, dimension( :, :, : ), allocatable :: Material_Absorption, Material_Absorption_Stab, &
           Velocity_Absorption, ScalarField_Absorption, Component_Absorption, Temperature_Absorption, &
!!$
           Component_Diffusion_Operator_Coefficient
      real, dimension( :, :, :, : ), allocatable :: Momentum_Diffusion, ScalarAdvectionField_Diffusion, &
           Component_Diffusion
      real, dimension( :, : ), allocatable :: Momentum_Diffusion_Vol

      real, dimension( :, : ), allocatable ::theta_flux, one_m_theta_flux, theta_flux_j, one_m_theta_flux_j, &
           sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j

      integer :: stat, istate, iphase, jphase, icomp, its, its2, cv_nodi, adapt_time_steps, cv_inod
      real, dimension( : ), allocatable :: rsum

      real, dimension(:, :), allocatable :: SUF_SIG_DIAGTEN_BC

      type( scalar_field ), pointer :: cfl, rc_field
      real :: c, rc, minc, maxc, ic
      !Variables for adaptive time stepping based on non-linear iterations
      logical :: nonLinearAdaptTs, Repeat_time_step, ExitNonLinearLoop
      real, dimension(:,:,:), allocatable  :: reference_field

      !Variables related to the deteccion and correction of bad elements
      real, parameter :: Max_bad_angle = 95.0
      real, parameter :: Min_bad_angle = 0.0
      type(bad_elements), allocatable, dimension(:) :: Quality_list
      !Variable related to help convergence
      logical :: help_convergence


      type( tensor_field ), pointer :: D_s, DC_s, DCOLD_s
      type( tensor_field ), pointer :: MFC_s, MFCOLD_s
      !! face value storage
      integer :: ncv_faces
      real::  second_theta

      integer :: checkpoint_number

      !Variable to store where we store things. Do not oversize this array, the size has to be the last index in use
      integer, dimension (37) :: StorageIndexes
      !Distribution of the indexes of StorageIndexes:
      !cv_fem_shape_funs_plus_storage: 1 (ASSEMB_FORCE_CTY), 13 (CV_ASSEMB)
      !CALC_ANISOTROP_LIM            : 2 (DETNLXR_PLUS_U_WITH_STORAGE in the inside, maybe 14 as well?)
      !DG_DERIVS_ALL2                : 3 (DETNLXR_PLUS_U_WITH_STORAGE in the inside, maybe 14 as well?)
      !DETNLXR_INVJAC                : 4
      !UNPACK_LOC                    : 5,6,7,8,9,10 (disabled)
      !COLOR_GET_CMC_PHA             : 11 (can be optimised, now it is not using only pointers)
      !Matrix C                      : 12
      !DG_DERIVS_ALL                 : 14 (DETNLXR_PLUS_U_WITH_STORAGE in the inside)
      !DETNLXR_PLUS_U_WITH_STORAGE   : 14
      !Indexes used in SURFACE_TENSION_WRAPPER (deprecated and will be removed):[15,30]
      !PROJ_CV_TO_FEM_state          : 31 (disabled)
      !Capillary pressure            : 32 (Pe), 33 (exponent a)
      !PIVIT_MAT (inverted)          : 34
      !Bound                         : 35
      !Ph 1                          : 36
      !Ph 2                          : 37


      !Working pointers

      type( tensor_field ), pointer :: tracer_field, velocity_field, density_field, saturation_field, old_saturation_field, tracer_source
      type(scalar_field), pointer :: pressure_field, porosity_field
      type(vector_field), pointer :: positions

      logical :: write_all_stats=.true.

      ! Variables used for calculating boundary outfluxes. Logical "calculate_flux" determines if this calculation is done. Intflux is the time integrated outflux
      ! Ioutlet counts the number of boundaries over which to calculate the outflux

      integer :: ioutlet

      real, dimension(:,:),  allocatable  :: intflux

      logical :: calculate_flux


#ifdef HAVE_ZOLTAN
      real(zoltan_float) :: ver
      integer(zoltan_int) :: ierr

      ierr = Zoltan_Initialize(ver)
      assert(ierr == ZOLTAN_OK)
#endif

      !Initially we set to use Stored data and that we have a new mesh
      StorageIndexes = 0!Initialize them as zero !


      !Read info for adaptive timestep based on non_linear_iterations

      if(have_option("/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep")) then
      if(have_option("/timestepping/nonlinear_iterations/nonlinear_iterations_at_adapt")) then
         call get_option('/timestepping/nonlinear_iterations/nonlinear_iterations_at_adapt',nonlinear_iterations_adapt)
         nonlinear_iterations = nonlinear_iterations_adapt
       end if

       call adapt_state_first_timestep(state)
       call allocate_and_insert_auxilliary_fields(state)

       ! Ensure that checkpoints do not adapt at first timestep.
       call delete_option(&
            "/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep")
    end if


      if(use_sub_state()) then
         call populate_sub_state(state,sub_state)
      end if

      !! JRP changes to make a multiphasic state
      call pack_multistate(state,packed_state,multiphase_state,&
           multicomponent_state)
      call set_boundary_conditions_values(state, shift_time=.true.)

      call set_caching_level()

      !      call initialize_rheologies(state,rheology)

      IDIVID_BY_VOL_FRAC=0
      !call print_state( packed_state )
      !stop 78

      !  Access boundary conditions via a call like
      !  call get_entire_boundary_condition(extract_tensor_field(packed_state,"Packed"//name),["weakdirichlet"],tfield,bc_type_list)
      !  where tfield is type(tensor_field) and bc_type_list is integer, dimension(tfield%dim(1),tfield%dim(2),nonods)
      !  Then values are in tfield%val(1/ndim/ncomp,nphase,nonods)
      !  Type ids are in bc_type_list(1/ndim/ncomp,nphase,stotel)
      !
      !A deallocate tfield when finished!!

      Repeat_time_step = .false.!Initially has to be false
      nonLinearAdaptTs = have_option(  '/timestepping/nonlinear_iterations/nonlinear_iterations_automatic/adaptive_timestep_nonlinear')

!!$ Compute primary scalars used in most of the code
      call Get_Primary_Scalars( state, &
           nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, ph_nloc=ph_nloc, ph_nonods=ph_nonods )

!!$ Calculating Global Node Numbers
      allocate( cv_sndgln( stotel * cv_snloc ), p_sndgln( stotel * p_snloc ), &
           u_sndgln( stotel * u_snloc ) )

      call Compute_Node_Global_Numbers( state, &
           totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
           cv_snloc, p_snloc, u_snloc, &
           cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln, &
           cv_sndgln, p_sndgln, u_sndgln )

!!$
!!$ Computing Sparsity Patterns Matrices
!!$
!!$ Defining lengths and allocating space for the matrices
      call Defining_MaxLengths_for_Sparsity_Matrices( state, ndim, nphase, totele, u_nloc, cv_nloc, ph_nloc, cv_nonods, &
           ph_nonods, mx_nface_p1, mxnele, mx_nct, mx_nc, mx_ncolcmc, mx_ncoldgm_pha, mx_ncolmcy, &
           mx_ncolacv, mx_ncolm, mx_ncolph )
      nlenmcy = u_nonods * nphase * ndim + cv_nonods
      allocate( finacv( cv_nonods * nphase + 1 ), colacv( mx_ncolacv ), midacv( cv_nonods * nphase ), &
           finmcy( nlenmcy + 1 ), colmcy( mx_ncolmcy ), midmcy( nlenmcy ), &
           findgm_pha( u_nonods * nphase * ndim + 1 ), coldgm_pha( mx_ncoldgm_pha ), &
           middgm_pha( u_nonods * nphase * ndim ), &
           findct( cv_nonods + 1 ), colct( mx_nct ), &
           findc( u_nonods + 1 ), colc( mx_nc ), &
           findcmc( cv_nonods + 1 ), colcmc( 0 ), midcmc( cv_nonods ), &
           findm( cv_nonods + 1 ), colm( mx_ncolm ), midm( cv_nonods ), &
           findph( ph_nonods + 1 ), colph( mx_ncolph ) )


      colct = 0 ; findc = 0 ; colc = 0 ; findcmc = 0 ; colcmc = 0 ; midcmc = 0 ; findm = 0
      colm = 0 ; midm = 0

!!$ Defining element-pair type
      call Get_Ele_Type( x_nloc, cv_ele_type, p_ele_type, u_ele_type, &
           mat_ele_type, u_sele_type, cv_sele_type )

!!$ Sparsity Patterns Matrices
      call Get_Sparsity_Patterns( state, &
!!$ CV multi-phase eqns (e.g. vol frac, temp)
           mx_ncolacv, ncolacv, finacv, colacv, midacv, &
           small_finacv, small_colacv, small_midacv, &
!!$ Force balance plus cty multi-phase eqns
           nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, &
!!$ Element connectivity
           mxnele, ncolele, midele, finele, colele, &
!!$ Force balance sparsity
           mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, &
!!$ CT sparsity - global continuity eqn
           mx_nct, ncolct, findct, colct, &
!!$ C sparsity operating on pressure in force balance
           mx_nc, ncolc, findc, colc, &
!!$ pressure matrix for projection method
           mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, &
!!$ CV-FEM matrix
           mx_ncolm, ncolm, findm, colm, midm, &
!!$ ph matrix
           mx_ncolph, ncolph, findph, colph, &
!!$ misc
           mx_nface_p1 )


      call temp_mem_hacks()

      Q_SCHEME = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/spatial_discretisation/control_volumes/q_scheme' )
      IGOT_THERM_VIS = 0
      IF ( Q_SCHEME ) IGOT_THERM_VIS = 1

!!$ Allocating space for various arrays:
      allocate( &
!!$
           DRhoDPressure( nphase, cv_nonods ), FEM_VOL_FRAC( nphase, cv_nonods ),&
!!$
           suf_sig_diagten_bc( stotel * cv_snloc * nphase, ndim ), &
           Mean_Pore_CV( cv_nonods ), &
           mass_ele( totele ), &
!!$
           Velocity_U_Source( ndim, nphase, u_nonods ), &
           Velocity_U_Source_CV( ndim, nphase, cv_nonods ), &
!!$
           Material_Absorption( mat_nonods, ndim * nphase, ndim * nphase ), &
           Velocity_Absorption( mat_nonods, ndim * nphase, ndim * nphase ), &
           Material_Absorption_Stab( mat_nonods, ndim * nphase, ndim * nphase ), &
           ScalarField_Absorption( nphase, nphase, cv_nonods ), Component_Absorption( nphase, nphase, cv_nonods ), &
           Temperature_Absorption( nphase, nphase, cv_nonods ), &
           Momentum_Diffusion( mat_nonods, ndim, ndim, nphase ), &
           Momentum_Diffusion_Vol( mat_nonods, nphase ), &
           ScalarAdvectionField_Diffusion( mat_nonods, ndim, ndim, nphase ), &
           Component_Diffusion( mat_nonods, ndim, ndim, nphase ), &
!!$ Variables used in the diffusion-like term: capilarity and surface tension:
           plike_grad_sou_grad( cv_nonods * nphase ), &
           plike_grad_sou_coef( cv_nonods * nphase ), &
           THERM_U_DIFFUSION(NDIM,NDIM,NPHASE,MAT_NONODS*IGOT_THERM_VIS ), THERM_U_DIFFUSION_VOL(NPHASE,MAT_NONODS*IGOT_THERM_VIS ) )

      ncv_faces=CV_count_faces( packed_state, CV_ELE_TYPE, stotel, cv_sndgln, u_sndgln)

!!$
      DRhoDPressure=0.
!!$
      suf_sig_diagten_bc=0.
!!$
      Mean_Pore_CV=0.
      mass_ele=0.
!!$
      Velocity_U_Source=0.
      Velocity_U_Source_CV=0.
!!$
      Material_Absorption=0.
      Velocity_Absorption=0.
      Material_Absorption_Stab=0.
      ScalarField_Absorption=0. ; Component_Absorption=0.
      Temperature_Absorption=0.
      Momentum_Diffusion=0.
      Momentum_Diffusion_Vol=0.
      ScalarAdvectionField_Diffusion=0.
      Component_Diffusion=0.
      THERM_U_DIFFUSION=0.
      THERM_U_DIFFUSION_VOL=0.
!!$
      plike_grad_sou_grad=0.
      plike_grad_sou_coef=0.
      iplike_grad_sou=0



!!$ Extracting Mesh Dependent Fields
      initialised = .false.
      call Extracting_MeshDependentFields_From_State( state, packed_state, initialised, &
           Velocity_U_Source, Velocity_Absorption )
!!$ Calculate diagnostic fields
      call calculate_diagnostic_variables( state, exclude_nonrecalculated = .true. )
      call calculate_diagnostic_variables_new( state, exclude_nonrecalculated = .true. )
!!$
!!$ Initialising Absorption terms that do not appear in the schema
!!$
      ScalarField_Absorption = 0. ; Component_Absorption = 0. ; Temperature_Absorption = 0.

!!$ Computing shape function scalars
      igot_t2 = 0 ; igot_theta_flux = 0
      if( ncomp /= 0 )then
         igot_t2 = 1 ; igot_theta_flux = 1
      end if

      call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi_theta, sbcvngi, nface, .false. )

      allocate( theta_flux( nphase, ncv_faces * igot_theta_flux ), &
           one_m_theta_flux( nphase, ncv_faces * igot_theta_flux ), &
           theta_flux_j( nphase, ncv_faces * igot_theta_flux ), &
           one_m_theta_flux_j( nphase, ncv_faces * igot_theta_flux ), &
           sum_theta_flux( nphase, ncv_faces * igot_theta_flux ), &
           sum_one_m_theta_flux( nphase, ncv_faces * igot_theta_flux ), &
           sum_theta_flux_j( nphase, ncv_faces * igot_theta_flux ), &
           sum_one_m_theta_flux_j( nphase, ncv_faces * igot_theta_flux ), &
           theta_gdiff( nphase, cv_nonods ), ScalarField_Source( nphase, cv_nonods ), &
           ScalarField_Source_Store( nphase, cv_nonods ), &
           ScalarField_Source_Component( nphase, cv_nonods ) )

      sum_theta_flux = 1. ; sum_one_m_theta_flux = 0.
      sum_theta_flux_j = 1. ; sum_one_m_theta_flux_j = 0.
      ScalarField_Source=0. ; ScalarField_Source_Store=0. ; ScalarField_Source_Component=0.

!!$ Defining discretisation options
      call Get_Discretisation_Options( state, &
           t_disopt, v_disopt, t_beta, v_beta, t_theta, v_theta, u_theta, &
           t_dg_vel_int_opt, u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
           comp_diffusion_opt, ncomp_diff_coef, in_ele_upwind, dg_ele_upwind, &
           nits_flux_lim_t, nits_flux_lim_volfra, nits_flux_lim_comp, &
           volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux, &
           t_use_theta_flux, t_get_theta_flux, scale_momentum_by_volume_fraction )

      allocate( Component_Diffusion_Operator_Coefficient( ncomp, ncomp_diff_coef, nphase ) )
      Component_Diffusion_Operator_Coefficient = 0.

!!$ Option not currently set up in the schema and zeroed from the begining. It is used to control
!!$ the upwinding rate (in the absorption term) during advection/assembling.

      allocate(opt_vel_upwind_coefs_new(ndim, ndim, nphase, mat_nonods)); opt_vel_upwind_coefs_new =0.
      allocate(opt_vel_upwind_grad_new(ndim, ndim, nphase, mat_nonods)); opt_vel_upwind_grad_new =0.

!!$ Defining problem to be solved:
      call get_option( '/material_phase[0]/vector_field::Velocity/prognostic/solver/max_iterations', &
           velocity_max_iterations,  default =  500 )
      call get_option( '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/solver/max_iterations', &
           PhaseVolumeFraction_max_iterations,  default =  500 )

      solve_force_balance = .false. ; solve_PhaseVolumeFraction = .false.
      if( velocity_max_iterations /= 0 ) solve_force_balance = .true.
      if( PhaseVolumeFraction_max_iterations /= 0 ) solve_PhaseVolumeFraction = .true.

!!$ Setting up variables for the Time- and NonLinear Iterations-Loops:
      call get_option( '/timestepping/current_time', acctim )
      call get_option( '/timestepping/timestep', dt )
      call get_option( '/timestepping/finish_time', finish_time )
      call get_option( '/io/dump_period_in_timesteps/constant', dump_period_in_timesteps, default = 1 )
      call get_option( '/timestepping/nonlinear_iterations', NonLinearIteration, default = 3 )
      !      call get_option( '/timestepping/nonlinear_iterations/nonlinear_iterations_automatic', tolerance_between_non_linear, default = -1. )
!!$
      have_temperature_field = .false. ; have_component_field = .false. ; have_extra_DiffusionLikeTerm = .false.
      do istate = 1, nstate
         if( have_option( '/material_phase[' // int2str( istate - 1 ) // ']/scalar_field::Temperature' ) ) &
              have_temperature_field = .true.
         if( have_option( '/material_phase[' // int2str( istate - 1 ) // ']/is_multiphase_component' ) ) &
              have_component_field = .true.
!!$
         call Calculate_All_Rhos( state, packed_state, ncomp, nphase, ndim, cv_nonods, cv_nloc, totele, &
              cv_ndgln, DRhoDPressure )

         if( have_component_field ) then
            call get_option( '/material_phase[' // int2str( istate - 1 ) // 'scalar_field::' // &
                 'ComponentMassFractionPhase1/prognostic/temporal_discretisation/control_volumes' // &
                 '/number_advection_iterations', NonLinearIteration_Components, default = 3 )
         end if
      end do

      if( have_option( '/material_phase[0]/multiphase_properties/capillary_pressure' ) ) &
           have_extra_DiffusionLikeTerm = .true.

      if ( have_option( '/mesh_adaptivity/hr_adaptivity' ) ) then
         call allocate( metric_tensor, extract_mesh(state(1), topology_mesh_name), 'ErrorMetric' )
      end if


      !Look for bad elements to apply a correction on them
!      if (is_porous_media) then
!          pressure_field=>extract_scalar_field(packed_state,"FEPressure")
!          allocate(Quality_list(totele*(NDIM-1)))!this number is not very well thought...
!          if (pressure_field%mesh%shape%degree < 2) &!Does not work yet for quadratic elements
!            call CheckElementAngles(packed_state, totele, x_ndgln, X_nloc, Max_bad_angle, Min_bad_angle, Quality_list)
!      end if



!!$ Starting Time Loop
      itime = 0
      dtime = 0

      ! When outlet_id is allocated, calculate_flux is true and we want to calculate outfluxes

      calculate_flux = allocated(outlet_id)

      ! If calculating boundary fluxes, initialise to zero time integrated fluxes (intflux) and the quantity (totout) used to calculate them.

      if(calculate_flux) then

          allocate(intflux(nphase,size(outlet_id)))
          allocate(totout(nphase, size(outlet_id)))

          do ioutlet = 1, size(outlet_id)
              intflux(:, ioutlet) = 0.

              totout(:, ioutlet) = 0.
          enddo
      endif

      checkpoint_number=1
      Loop_Time: do

!Always help to converge in the first non-linear iteration if it is adaptive
help_convergence = have_option('/timestepping/nonlinear_iterations/nonlinear_iterations_automatic')

!!$

         ewrite(2,*) '    NEW DT', itime+1


         sum_theta_flux_j = 1. ; sum_one_m_theta_flux_j = 0.

         if (do_checkpoint_simulation(dtime)) then
            call checkpoint_simulation(state,cp_no=checkpoint_number,&
                 protect_simulation_name=.true.,file_type='.mpml')
            checkpoint_number=checkpoint_number+1
         end if
         dtime=dtime+1


         itime = itime + 1
         timestep = itime
         call get_option( '/timestepping/timestep', dt )


         acctim = acctim + dt
         call set_option( '/timestepping/current_time', acctim )
         new_lim = .true.

         ! Added a tolerance of 0.001dt to the condition below that stops us exiting the loop before printing the last time step.

          if(calculate_flux) then
            if ( acctim > finish_time + 0.001*dt) then
              ewrite(1,*) "Passed final time"
              exit Loop_Time
            end if

          else

            if ( acctim > finish_time) then
              ewrite(1,*) "Passed final time"
              exit Loop_Time
            end if

         endif

!         if ( acctim > finish_time) then
!            ewrite(1,*) "Passed final time"
!            exit Loop_Time
!         end if

         call get_option( '/timestepping/final_timestep', final_timestep, stat )
         if( stat == spud_no_error ) then
            if( itime > final_timestep ) then
               ewrite(1,*) "Passed final timestep"
               exit Loop_Time
            end if
         end if

         ExitNonLinearLoop = .false.
         !Store backup to be able to repeat a timestep
         if (nonLinearAdaptTs) call Adaptive_NonLinear(packed_state, reference_field, its, &
              Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,1)


         porosity_field=>extract_scalar_field(packed_state,"Porosity")

         ! evaluate prescribed fields at time = current_time+dt
         call set_prescribed_field_values( state, exclude_interpolated = .true., &
              exclude_nonreprescribed = .true., time = acctim )
         !! Update all fields from time-step 'N - 1'
         call copy_packed_new_to_old( packed_state )

         ! update velocity absorption
         call update_velocity_absorption( state, ndim, nphase, mat_nonods, velocity_absorption )

!!$ FEMDEM...
#ifdef USING_FEMDEM
        if ( have_option( '/femdem_fracture' ) ) then
            call fracking(packed_state) 
        elseif ( have_option( '/blasting') ) then 
            call blasting( packed_state, nphase )
            call update_blasting_memory( packed_state, state, timestep )  
        end if
#endif

!!$ Start non-linear loop
         Loop_NonLinearIteration: do  its = 1, NonLinearIteration

            ewrite(2,*) '  NEW ITS', its

            call calculate_rheologies(state,rheology)
            !To force the recalculation of all the stored variables uncomment the following line:
!           call Clean_Storage(state, StorageIndexes)

            !call set_nu_to_u( packed_state )
            !call boiling( state, packed_state, cv_nonods, mat_nonods, nphase, ndim, &
            !   ScalarField_Source, velocity_absorption, temperature_absorption )


            if( have_temperature_field .and. &
                 have_option( '/material_phase[0]/scalar_field::Temperature/prognostic' ) ) then
               call get_option( '/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation' // &
                    '/control_volumes/second_theta', second_theta, default=1. )
            end if

            if( have_component_field ) then
               call get_option( '/material_phase[' // int2str( nphase ) // ']/scalar_field::ComponentMassFractionPhase1/' // &
                    'prognostic/temporal_discretisation/control_volumes/second_theta', second_theta, default=1. )
            end if

            !Store the field we want to compare with to check how are the computations going
            call Adaptive_NonLinear(packed_state, reference_field, its, &
                 Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,2)

            call Calculate_All_Rhos( state, packed_state, ncomp, nphase, ndim, cv_nonods, cv_nloc, totele, &
                 cv_ndgln, DRhoDPressure )
            if( solve_force_balance ) then
               call Calculate_AbsorptionTerm( state, packed_state,&
                    cv_ndgln, mat_ndgln, &
                    opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, Material_Absorption )

               ! calculate SUF_SIG_DIAGTEN_BC this is \sigma_in^{-1} \sigma_out
               ! \sigma_in and \sigma_out have the same anisotropy so SUF_SIG_DIAGTEN_BC
               ! is diagonal
               if( is_porous_media ) then
                  call calculate_SUF_SIG_DIAGTEN_BC( packed_state, suf_sig_diagten_bc, totele, stotel, cv_nloc, &
                       cv_snloc, nphase, ndim, nface, mat_nonods, cv_nonods, x_nloc, ncolele, cv_ele_type, &
                       finele, colele, cv_ndgln, cv_sndgln, x_ndgln, mat_ndgln, material_absorption, &
                       state, x_nonods )
               end if
            end if
!!$ Solve advection of the scalar 'Temperature':



            Conditional_ScalarAdvectionField: if( have_temperature_field .and. &
                 have_option( '/material_phase[0]/scalar_field::Temperature/prognostic' ) ) then
               ewrite(3,*)'Now advecting Temperature Field'

               call set_nu_to_u( packed_state )

               call calculate_diffusivity( state, ncomp, nphase, ndim, cv_nonods, mat_nonods, &
                    mat_nloc, totele, mat_ndgln, ScalarAdvectionField_Diffusion )

               tracer_field=>extract_tensor_field(packed_state,"PackedTemperature")
               velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
               density_field=>extract_tensor_field(packed_state,"PackedDensity",stat)
               saturation_field=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")

               call INTENERGE_ASSEM_SOLVE( state, packed_state, &
                    tracer_field,velocity_field,density_field,&
                    small_FINACV, small_COLACV, small_MIDACV, &
                    NCOLCT, FINDCT, COLCT, &
                    CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
                    U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE,  &
                    NPHASE, &
                    CV_NLOC, U_NLOC, X_NLOC, &
                    CV_NDGLN, X_NDGLN, U_NDGLN, &
                    CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
!!$
                    MAT_NLOC, MAT_NDGLN, MAT_NONODS, ScalarAdvectionField_Diffusion, IGOT_THERM_VIS, THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL, &
                    t_disopt, t_dg_vel_int_opt, dt, t_theta, t_beta, &
                    suf_sig_diagten_bc,&
                    DRhoDPressure, &
                    Temperature_Absorption, Porosity_field%val, &
                    ndim, &
!!$
                    NCOLM, FINDM, COLM, MIDM, &
!!$
                    XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
!!$
                    opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
                    0, igot_theta_flux, scvngi_theta, &
                    t_get_theta_flux, t_use_theta_flux, &
                    THETA_GDIFF, &
                    in_ele_upwind, dg_ele_upwind, &
                    Mean_Pore_CV, &
                    option_path = '/material_phase[0]/scalar_field::Temperature', &
                    thermal = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/equation::InternalEnergy'),&
                    StorageIndexes=StorageIndexes, saturation=saturation_field )

               call Calculate_All_Rhos( state, packed_state, ncomp, nphase, ndim, cv_nonods, cv_nloc, totele, &
                    cv_ndgln, DRhoDPressure )

            end if Conditional_ScalarAdvectionField

            ScalarField_Source_Store = ScalarField_Source + ScalarField_Source_Component
            tracer_source => extract_tensor_field(packed_state,"PackedPhaseVolumeFractionSource")

            volfra_use_theta_flux = .true.
            if( ncomp <= 1 ) volfra_use_theta_flux = .false.

!!$ Now solving the Momentum Equation ( = Force Balance Equation )
            Conditional_ForceBalanceEquation: if ( solve_force_balance ) then

               call set_nu_to_u( packed_state )

!!$ Diffusion-like term -- here used as part of the capillary pressure for porous media. It can also be
!!$ extended to surface tension -like term.
               iplike_grad_sou = 0
               plike_grad_sou_grad = 0


               CALL CALCULATE_SURFACE_TENSION( state, packed_state, nphase, ncomp, &
                    PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, IPLIKE_GRAD_SOU, &
                    Velocity_U_Source_CV, Velocity_U_Source, &
                    NCOLACV, FINACV, COLACV, MIDACV, &
                    small_FINACV, small_COLACV, small_MIDACV, &
                    NCOLCT, FINDCT, COLCT, &
                    CV_NONODS, U_NONODS, X_NONODS, TOTELE, STOTEL, &
                    CV_ELE_TYPE, CV_SELE_TYPE, U_ELE_TYPE, &
                    CV_NLOC, U_NLOC, X_NLOC, CV_SNLOC, U_SNLOC, &
                    CV_NDGLN, CV_SNDGLN, X_NDGLN, U_NDGLN, U_SNDGLN, &
                    MAT_NLOC, MAT_NDGLN, MAT_NONODS,  &
                    NDIM,  &
                    NCOLM, FINDM, COLM, MIDM, &
                    XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
                    StorageIndexes=StorageIndexes )

               if( have_option( '/material_phase[0]/multiphase_properties/capillary_pressure' ) )then
                  call calculate_capillary_pressure( state, packed_state, .true., StorageIndexes)
               end if

               velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
               pressure_field=>extract_scalar_field(packed_state,"FEPressure")

               CALL FORCE_BAL_CTY_ASSEM_SOLVE( state, packed_state, &
                    velocity_field, pressure_field, &
                    NDIM, NPHASE, NCOMP, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
                    U_ELE_TYPE, P_ELE_TYPE, &
                    U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
                    U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN,&
                    STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
                    U_SNLOC, P_SNLOC, CV_SNLOC, &
!!$
                    Material_Absorption_Stab, Material_Absorption, Velocity_Absorption, Velocity_U_Source, Velocity_U_Source_CV, &
                    DRhoDPressure, IDIVID_BY_VOL_FRAC, FEM_VOL_FRAC, &
                    dt, &
!!$
                    NCOLC, FINDC, COLC, & ! C sparsity - global cty eqn
                    NCOLDGM_PHA, &! Force balance
                    NCOLELE, FINELE, COLELE, & ! Element connectivity.
                    NCOLCMC, FINDCMC, COLCMC, MIDCMC, & ! pressure matrix for projection method
                    size(small_colacv),small_FINACV, small_COLACV, small_MIDACV, &
                    NLENMCY, NCOLMCY, FINMCY, COLMCY, MIDMCY, & ! Force balance plus cty multi-phase eqns
                    NCOLCT, FINDCT, COLCT, & ! CT sparsity - global cty eqn.
                    CV_ELE_TYPE, &
!!$
                    v_disopt, v_dg_vel_int_opt, v_theta, &
                    SUF_SIG_DIAGTEN_BC, &
                    ScalarField_Source_Store, ScalarField_Absorption, Porosity_field%val, &
!!$
                    NCOLM, FINDM, COLM, MIDM, & ! Sparsity for the CV-FEM
                    XU_NLOC, XU_NDGLN, &
!!$
                    Momentum_Diffusion, Momentum_Diffusion_Vol, THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL, &
                    opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
                    igot_theta_flux, scvngi_theta, volfra_use_theta_flux, &
                    sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j, &
                    in_ele_upwind, dg_ele_upwind, &
                    iplike_grad_sou, plike_grad_sou_coef, plike_grad_sou_grad, &
                    scale_momentum_by_volume_fraction,&
                    StorageIndexes=StorageIndexes, Quality_list = Quality_list,&
                    nonlinear_iteration = its )
!!$ Calculate Density_Component for compositional
               if( have_component_field ) &
                    call Calculate_Component_Rho( state, packed_state, &
                    ncomp, nphase, cv_nonods )

            end if Conditional_ForceBalanceEquation

            Conditional_PhaseVolumeFraction: if ( solve_PhaseVolumeFraction ) then
               call VolumeFraction_Assemble_Solve( state, packed_state, &
                    small_FINACV, small_COLACV, small_MIDACV, &
                    NCOLCT, FINDCT, COLCT, &
                    CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
                    CV_ELE_TYPE, &
                    NPHASE, &
                    CV_NLOC, U_NLOC, X_NLOC,  &
                    CV_NDGLN, X_NDGLN, U_NDGLN, &
                    CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
!!$
                    MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
!!$
                    v_disopt, v_dg_vel_int_opt, dt, v_theta, v_beta, &
                    SUF_SIG_DIAGTEN_BC, &
                    DRhoDPressure, &
                    ScalarField_Source_Store, ScalarField_Absorption, Porosity_field%val, &
!!$
                    NDIM,nface, &
                    NCOLM, FINDM, COLM, MIDM, &
                    XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
!!$
                    opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
                    igot_theta_flux,scvngi_theta, volfra_use_theta_flux, &
                    in_ele_upwind, dg_ele_upwind, &
                    option_path = '/material_phase[0]/scalar_field::PhaseVolumeFraction', &
                    mass_ele_transp = mass_ele,&
                    theta_flux=sum_theta_flux, one_m_theta_flux=sum_one_m_theta_flux, &
                    theta_flux_j=sum_theta_flux_j, one_m_theta_flux_j=sum_one_m_theta_flux_j,&
                    StorageIndexes=StorageIndexes, Material_Absorption=Material_Absorption,&
                    nonlinear_iteration = its, help_convergence = help_convergence)

            end if Conditional_PhaseVolumeFraction


!!$ Starting loop over components
            sum_theta_flux = 0. ; sum_one_m_theta_flux = 0. ; sum_theta_flux_j = 0. ; sum_one_m_theta_flux_j = 0. ; ScalarField_Source_Component = 0.

            tracer_source%val = 0.

            velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
            saturation_field=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
            old_saturation_field=>extract_tensor_field(packed_state,"PackedOldPhaseVolumeFraction")


            Conditional_Components:if( have_component_field ) then

               D_s  => extract_tensor_field( packed_state, "PackedDensity" )

               DC_s  => extract_tensor_field( packed_state, "PackedComponentDensity" )
               DCOLD_s  => extract_tensor_field( packed_state, "PackedOldComponentDensity" )
               MFC_s  => extract_tensor_field( packed_state, "PackedComponentMassFraction" )
               MFCOLD_s  => extract_tensor_field( packed_state, "PackedOldComponentMassFraction" )


               Loop_Components: do icomp = 1, ncomp

                  tracer_field=>extract_tensor_field(multicomponent_state(icomp),"PackedComponentMassFraction")
                  density_field=>extract_tensor_field(multicomponent_state(icomp),"PackedComponentDensity",stat)

!!$ Computing the absorption term for the multi-components equation
                  call Calculate_ComponentAbsorptionTerm( state, packed_state, &
                       icomp, cv_ndgln, &
                       D_s%val, Porosity_field%val, mass_ele, &
                       Component_Absorption )

                  Conditional_SmoothAbsorption: if( have_option( '/material_phase[' // int2str( nstate - ncomp ) // &
                       ']/is_multiphase_component/KComp_Sigmoid' ) .and. nphase > 1 ) then
                     do cv_nodi = 1, cv_nonods
                        if(saturation_field%val(1,1, cv_nodi ) > 0.95 ) then
                           do iphase = 1, nphase
                              do jphase = min( iphase + 1, nphase ), nphase
                                 Component_Absorption( iphase, jphase, cv_nodi ) = &
                                      Component_Absorption( iphase, jphase, cv_nodi ) * max( 0.01, &
                                      20. * ( 1. - saturation_field%val( 1,1, cv_nodi ) ) )
                              end do
                           end do
                        end if

                     end do
                  end if Conditional_SmoothAbsorption

!!$ Computing diffusion term for the component conservative equation:
                  call Calculate_ComponentDiffusionTerm( state, packed_state, &
                       mat_ndgln, u_ndgln, x_ndgln, &
                       u_ele_type, p_ele_type, ncomp_diff_coef, comp_diffusion_opt, &
                       Component_Diffusion_Operator_Coefficient( icomp, :, : ), &
                       Component_Diffusion ,&
                       StorageIndexes=StorageIndexes )

!!$ NonLinear iteration for the components advection:
                  Loop_NonLinearIteration_Components: do its2 = 1, NonLinearIteration_Components
                     comp_use_theta_flux = .false. ; comp_get_theta_flux = .true.

                     call INTENERGE_ASSEM_SOLVE( state, multicomponent_state(icomp), &
                          tracer_field,velocity_field,density_field,&
                          SMALL_FINACV, SMALL_COLACV, small_MIDACV,&
                          NCOLCT, FINDCT, COLCT, &
                          CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
                          U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE,  &
                          NPHASE,  &
                          CV_NLOC, U_NLOC, X_NLOC,  &
                          CV_NDGLN, X_NDGLN, U_NDGLN, &
                          CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
!!$
                          MAT_NLOC, MAT_NDGLN, MAT_NONODS, Component_Diffusion, 0, THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL,&
                          v_disopt, v_dg_vel_int_opt, dt, v_theta, v_beta, &
                          SUF_SIG_DIAGTEN_BC,&
                          DRhoDPressure, &
                          Component_Absorption, Porosity_field%val, &
!!$
                          NDIM,  &
                          NCOLM, FINDM, COLM, MIDM, &
                          XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
!!$
                          opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
                          igot_t2, igot_theta_flux, scvngi_theta, &
                          comp_get_theta_flux, comp_use_theta_flux, &
                          theta_gdiff, &
                          in_ele_upwind, dg_ele_upwind, &
                          Mean_Pore_CV, &
                          thermal = .false.,& ! the false means that we don't add an extra source term
                          theta_flux=theta_flux, one_m_theta_flux=one_m_theta_flux, theta_flux_j=theta_flux_j, one_m_theta_flux_j=one_m_theta_flux_j,&
                          StorageIndexes=StorageIndexes, icomp=icomp, saturation=saturation_field )

!                      tracer_field%val = min (max( tracer_field%val, 0.0), 1.0)

                  end do Loop_NonLinearIteration_Components

                  sum_theta_flux = sum_theta_flux + theta_flux
                  sum_one_m_theta_flux = sum_one_m_theta_flux + one_m_theta_flux

                  sum_theta_flux_j = sum_theta_flux_j + theta_flux_j
                  sum_one_m_theta_flux_j = sum_one_m_theta_flux_j + one_m_theta_flux_j


                  ! We have divided through by density
                  !do cv_inod=1,cv_nonods
                  !   do iphase=1,nphase
                  !      ScalarField_Source_Component((iphase-1)*cv_nonods+cv_inod) = &
                  !               ScalarField_Source_Component((iphase-1)*cv_nonods+cv_inod) + THETA_GDIFF(iphase,cv_inod)
                  !   end do
                  !end do
                  tracer_source%val(1,:,:) = tracer_source%val(1,:,:) + THETA_GDIFF

               end do Loop_Components




               if( have_option( '/material_phase[' // int2str( nstate - ncomp ) // &
                    ']/is_multiphase_component/Comp_Sum2One/Enforce_Comp_Sum2One' ) ) then
                  ! Initially clip and then ensure the components sum to unity so we don't get surprising results...
                  MFC_s % val = min ( max ( MFC_s % val, 0.0), 1.0)

                  ALLOCATE( RSUM( NPHASE ) )
                  DO CV_INOD = 1, CV_NONODS
                     DO IPHASE = 1, NPHASE
                        RSUM( IPHASE ) = SUM (MFC_s % val (:, IPHASE, CV_INOD) )
                     END DO
                     DO IPHASE = 1, NPHASE
                        MFC_s % val (:, IPHASE, CV_INOD) = MFC_s % val (:, IPHASE, CV_INOD) / RSUM( IPHASE )
                     END DO
                  END DO
                  DEALLOCATE( RSUM )
               end if






               DO ICOMP = 1, NCOMP

                  call Calculate_ComponentAbsorptionTerm( state, packed_state,&
                       icomp, cv_ndgln, &
                       D_s%val, Porosity_field%val, mass_ele, &
                       Component_Absorption )

                  if( have_option( '/material_phase[' // int2str( nstate - ncomp ) // &
                       ']/is_multiphase_component/KComp_Sigmoid' ) .and. nphase > 1 ) then
                     do cv_nodi = 1, cv_nonods
                        if( saturation_field%val( 1, 1, cv_nodi ) > 0.95 ) then
                           do iphase = 1, nphase
                              do jphase = min( iphase + 1, nphase ), nphase
                                 Component_Absorption( iphase, jphase, cv_nodi ) = &
                                      Component_Absorption( iphase, jphase, cv_nodi ) * max( 0.01, &
                                      20. * ( 1. - saturation_field%val (1,1, cv_nodi ) ) )
                              end do
                           end do
                        end if
                     end do
                  end if

                  DO CV_NODI = 1, CV_NONODS
                     Loop_Phase_SourceTerm1: do iphase = 1, nphase
                        Loop_Phase_SourceTerm2: do jphase = 1, nphase
                           tracer_source%val(1,iphase,cv_nodi)=tracer_source%val(1,iphase,cv_nodi)- &
                                Component_Absorption( IPHASE, JPHASE, CV_NODI ) * &
                                MFC_s%val(ICOMP, JPHASE, CV_NODI) / &
                                DC_s%val( icomp, iphase, cv_nodi  )
                        end do Loop_Phase_SourceTerm2
                     end do Loop_Phase_SourceTerm1
                  END DO

                  ! For compressibility
                  DO IPHASE = 1, NPHASE
                     DO CV_NODI = 1, CV_NONODS
                        tracer_source%val(1,iphase,cv_nodi)=tracer_source%val(1,iphase,cv_nodi)&
                             + Mean_Pore_CV( CV_NODI ) * MFCOLD_s%val(ICOMP, IPHASE, CV_NODI) &
                             * ( DCOLD_s%val( ICOMP, IPHASE, CV_NODI ) - DC_s%val( ICOMP, IPHASE, CV_NODI) ) &
                             * old_saturation_field%val( 1,IPHASE, CV_NONODS ) &
                             / ( DC_s%val( ICOMP, IPHASE, CV_NODI ) * DT )
                     END DO
                  END DO

               END DO ! ICOMP


               if( have_option( '/material_phase[' // int2str( nstate - ncomp ) // &
                    ']/is_multiphase_component/Comp_Sum2One' ) .and. ( ncomp > 1 ) ) then
                  call Cal_Comp_Sum2One_Sou( packed_state, cv_nonods, nphase, ncomp, dt, its, &
                       NonLinearIteration, Mean_Pore_CV )
               end if


            end if Conditional_Components


            !Check if the results are good so far and act in consequence, only does something if requested by the user
            call Adaptive_NonLinear(packed_state, reference_field, its,&
                 Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,3, help_convergence = help_convergence)
            if (ExitNonLinearLoop) exit Loop_NonLinearIteration

         end do Loop_NonLinearIteration


         ! If calculating boundary fluxes, add up contributions to \int{totout} at each time step

         if(calculate_flux) then
             intflux = intflux + totout*dt

                 ! We will output totout normalised as a fractional flow

                 totout = totout/sum(totout)

         endif

         ! If calculating boundary fluxes, dump them to outfluxes.txt

         if(calculate_flux) then
             if(getprocno() == 1) then

                 call dump_outflux(acctim,itime,totout,intflux, dt)

             endif
         endif


         if (nonLinearAdaptTs) then
            !As the value of dt and acctim may have changed we retrieve their values
            !to make sure that everything is coherent
            call get_option( '/timestepping/current_time', acctim )
            call get_option( '/timestepping/timestep', dt)
            !If repeat timestep we don't want to adapt mesh or dump results
            if ( Repeat_time_step ) then
               itime = itime - 1
               cycle Loop_Time
            end if
         end if

         call set_option( '/timestepping/current_time', acctim )
         call set_option( '/timestepping/timestep', dt)

         current_time = acctim


!!$ Calculate diagnostic fields
         call calculate_diagnostic_variables( state, exclude_nonrecalculated = .true. )
         call calculate_diagnostic_variables_new( state, exclude_nonrecalculated = .true. )


         if (write_all_stats) call write_diagnostics( state, current_time, dt, itime )  ! Write stat file


         Conditional_TimeDump: if( ( mod( itime, dump_period_in_timesteps ) == 0 ) ) then

            dtime=dtime+1
            if (do_checkpoint_simulation(dtime)) then
               call checkpoint_simulation(state,cp_no=checkpoint_number,&
                    protect_simulation_name=.true.,file_type='.mpml')
               checkpoint_number=checkpoint_number+1
            end if

            call get_option( '/timestepping/current_time', current_time ) ! Find the current time

            if (.not. write_all_stats)call write_diagnostics( state, current_time, dt, itime/dump_period_in_timesteps )  ! Write stat file
            not_to_move_det_yet = .false. ; dump_no = itime/dump_period_in_timesteps ! Sync dump_no with itime
            call write_state( dump_no, state ) ! Now writing into the vtu files

         end if Conditional_TimeDump

!!$! ******************
!!$! *** Mesh adapt ***
!!$! ******************

         do_reallocate_fields = .false.
         Conditional_Adaptivity_ReallocatingFields: if( have_option( '/mesh_adaptivity/hr_adaptivity') ) then
            if( have_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps') ) then
               call get_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps', &
                    adapt_time_steps, default=5 )
            end if
            if( mod( itime, adapt_time_steps ) == 0 ) do_reallocate_fields = .true.
         elseif (have_option( '/mesh_adaptivity/hr_adaptivity_prescribed_metric') ) then
            if( have_option( '/mesh_adaptivity/hr_adaptivity_prescribed_metric/period_in_timesteps') ) then
               call get_option( '/mesh_adaptivity/hr_adaptivity_prescribed_metric/period_in_timesteps', &
                    adapt_time_steps, default=5 )
            end if
            if( mod( itime, adapt_time_steps ) == 0 ) do_reallocate_fields = .true.
         elseif( have_option( '/mesh_adaptivity/prescribed_adaptivity' ) ) then
            if( do_adapt_state_prescribed( current_time ) ) do_reallocate_fields = .true.
         end if Conditional_Adaptivity_ReallocatingFields

         new_mesh = do_reallocate_fields

         Conditional_ReallocatingFields: if( do_reallocate_fields ) then
            !The storaged variables must be recalculated
            call Clean_Storage(state, StorageIndexes)

!            call linearise_components()

!            tfield => extract_tensor_field( packed_state, "PackedComponentMassFraction" )
!            call linearise( tfield )

!            tfield => extract_tensor_field( packed_state, "PackedDensity" )
!            call linearise( tfield )

!            tfield => extract_tensor_field( packed_state, "PackedTemperature" )
!            call linearise( tfield )

!            tfield => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
!            call linearise( tfield )


            Conditional_Adaptivity: if( have_option( '/mesh_adaptivity/hr_adaptivity ') .or. have_option( '/mesh_adaptivity/hr_adaptivity_prescribed_metric')) then

               Conditional_Adapt_by_TimeStep: if( mod( itime, adapt_time_steps ) == 0 ) then

                  call pre_adapt_tasks( sub_state )

                  if (have_option('/mesh_adaptivity/hr_adaptivity_prescribed_metric')) then
                     positions=>extract_vector_field(state(1),"Coordinate")
                     call allocate( metric_tensor, extract_mesh(state(1), topology_mesh_name), 'MetricTensor' )
                     call initialise_field(metric_tensor,'/mesh_adaptivity/hr_adaptivity_prescribed_metric/tensor_field::MetricTensor',positions)
                     nullify(positions)
                  else
                     call qmesh( state, metric_tensor )
                  end if

                  if( have_option( '/io/stat/output_before_adapts' ) ) call write_diagnostics( state, current_time, dt, &
                       itime, not_to_move_det_yet = .true. )
!!$                  if( have_option( '/io/stat/output_before_adapts' ) ) call write_diagnostics( state, current_time, dt, &
!!$                       timestep, not_to_move_det_yet = .true. )

                  call run_diagnostics( state )

                  call adapt_state( state, metric_tensor, suppress_reference_warnings=.true. )

                  call update_state_post_adapt( state, metric_tensor, dt, sub_state, nonlinear_iterations, &
                       nonlinear_iterations_adapt )

                  if( have_option( '/io/stat/output_after_adapts' ) ) call write_diagnostics( state, current_time, dt, &
                       itime, not_to_move_det_yet = .true. )
!!$                  if( have_option( '/io/stat/output_after_adapts' ) ) call write_diagnostics( state, current_time, dt, &
!!$                       timestep, not_to_move_det_yet = .true. )

                  call run_diagnostics( state )
               end if Conditional_Adapt_by_TimeStep

            elseif( have_option( '/mesh_adaptivity/prescribed_adaptivity' ) ) then !!$ Conditional_Adaptivity:

               Conditional_Adapt_by_Time: if( do_adapt_state_prescribed( current_time ) ) then

                  call pre_adapt_tasks( sub_state )

                  if( have_option( '/io/stat/output_before_adapts' ) ) call write_diagnostics( state, current_time, dt, &
                       timestep, not_to_move_det_yet = .true. )

                  call run_diagnostics( state )

                  call adapt_state_prescribed( state, current_time )

                  call update_state_post_adapt( state, metric_tensor, dt, sub_state, nonlinear_iterations, &
                       nonlinear_iterations_adapt)

                  if(have_option( '/io/stat/output_after_adapts' ) ) call write_diagnostics( state, current_time, dt, &
                       timestep, not_to_move_det_yet = .true. )

                  call run_diagnostics( state )

               end if Conditional_Adapt_by_Time

               not_to_move_det_yet = .false.

            end if Conditional_Adaptivity

            call deallocate(packed_state)
            call deallocate(multiphase_state)
            call deallocate(multicomponent_state )
            !call unlinearise_components()
            call pack_multistate(state,packed_state,&
                 multiphase_state,multicomponent_state)
            call set_boundary_conditions_values(state, shift_time=.true.)

            if (allocated(Quality_list) ) deallocate(Quality_list)

!!$ Deallocating array variables:
            deallocate( &
!!$ Node glabal numbers
                 cv_sndgln, p_sndgln, u_sndgln, &
!!$ Sparsity patterns
                 finacv, colacv, midacv,&
                 small_finacv, small_colacv, small_midacv, &
                 finmcy, colmcy, midmcy, &
                 findgm_pha, coldgm_pha, middgm_pha, findct, &
                 colct, findc, colc, findcmc, colcmc, midcmc, findm, &
                 colm, midm, &
!!$ Defining element-pair type and discretisation options and coefficients
                 opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
!!$ For output:
                 Mean_Pore_CV, &
!!$ Variables used in the diffusion-like term: capilarity and surface tension:
                 plike_grad_sou_grad, plike_grad_sou_coef, &
!!$ Working arrays
                 DRhoDPressure, &
                 Velocity_U_Source, Velocity_U_Source_CV, &
                 suf_sig_diagten_bc, &
                 theta_gdiff, ScalarField_Source, ScalarField_Source_Store, ScalarField_Source_Component, &
                 mass_ele, &
                 Material_Absorption, Material_Absorption_Stab, &
                 Velocity_Absorption, ScalarField_Absorption, Component_Absorption, Temperature_Absorption, &
                 Component_Diffusion_Operator_Coefficient, &
                 Momentum_Diffusion, Momentum_Diffusion_Vol, ScalarAdvectionField_Diffusion, &
                 Component_Diffusion, &
                 theta_flux, one_m_theta_flux, theta_flux_j, one_m_theta_flux_j, sum_theta_flux, &
                 sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j )


!!$  Compute primary scalars used in most of the code
            call Get_Primary_Scalars( state, &
                 nphase, nstate, ncomp, totele, ndim, stotel, &
                 u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
                 x_snloc, cv_snloc, u_snloc, p_snloc, &
                 cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, ph_nloc=ph_nloc, ph_nonods=ph_nonods )
!!$ Calculating Global Node Numbers
            allocate( cv_sndgln( stotel * cv_snloc ), p_sndgln( stotel * p_snloc ), &
                 u_sndgln( stotel * u_snloc ) )

            !          x_ndgln_p1 = 0 ; x_ndgln = 0 ; cv_ndgln = 0 ; p_ndgln = 0 ; mat_ndgln = 0 ; u_ndgln = 0 ; xu_ndgln = 0 ; &
            cv_sndgln = 0 ; p_sndgln = 0 ; u_sndgln = 0

            call Compute_Node_Global_Numbers( state, &
                 totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
                 cv_snloc, p_snloc, u_snloc, &
                 cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln, &
                 cv_sndgln, p_sndgln, u_sndgln )
!!$
!!$ Computing Sparsity Patterns Matrices
!!$

!!$ Defining lengths and allocating space for the matrices
            call Defining_MaxLengths_for_Sparsity_Matrices( state, ndim, nphase, totele, u_nloc, cv_nloc, ph_nloc, cv_nonods, &
                 ph_nonods, mx_nface_p1, mxnele, mx_nct, mx_nc, mx_ncolcmc, mx_ncoldgm_pha, mx_ncolmcy, &
                 mx_ncolacv, mx_ncolm, mx_ncolph )
            nlenmcy = u_nonods * nphase * ndim + cv_nonods
            allocate( finacv( cv_nonods * nphase + 1 ), colacv( mx_ncolacv ), midacv( cv_nonods * nphase ), &
                 finmcy( nlenmcy + 1 ), colmcy( mx_ncolmcy ), midmcy( nlenmcy ), &
                 findgm_pha( u_nonods * nphase * ndim + 1 ), coldgm_pha( mx_ncoldgm_pha ), &
                 middgm_pha( u_nonods * nphase * ndim ), &
                 findct( cv_nonods + 1 ), colct( mx_nct ), &
                 findc( u_nonods + 1 ), colc( mx_nc ), &
                 findcmc( cv_nonods + 1 ), colcmc( mx_ncolcmc ), midcmc( cv_nonods ), &
                 findm( cv_nonods + 1 ), colm( mx_ncolm ), midm( cv_nonods ), &
                 findph( ph_nonods + 1 ), colph( mx_ncolph ) )

            finacv = 0 ; colacv = 0 ; midacv = 0 ; finmcy = 0 ; colmcy = 0 ; midmcy = 0 ; &
                 findgm_pha = 0 ; coldgm_pha = 0 ; middgm_pha = 0 ; findct = 0 ; &
                 colct = 0 ; findc = 0 ; colc = 0 ; findcmc = 0 ; colcmc = 0 ; midcmc = 0 ; findm = 0 ; &
                 colm = 0 ; midm = 0

!!$ Defining element-pair type
            call Get_Ele_Type( x_nloc, cv_ele_type, p_ele_type, u_ele_type, &
                 mat_ele_type, u_sele_type, cv_sele_type )

!!$ Sparsity Patterns Matrices
            call Get_Sparsity_Patterns( state, &
!!$ CV multi-phase eqns (e.g. vol frac, temp)
                 mx_ncolacv, ncolacv, finacv, colacv, midacv, &
                 small_finacv, small_colacv, small_midacv, &
!!$ Force balance plus cty multi-phase eqns
                 nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, &
!!$ Element connectivity
                 mxnele, ncolele, midele, finele, colele, &
!!$ Force balance sparsity
                 mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, &
!!$ CT sparsity - global continuity eqn
                 mx_nct, ncolct, findct, colct, &
!!$ C sparsity operating on pressure in force balance
                 mx_nc, ncolc, findc, colc, &
!!$ pressure matrix for projection method
                 mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, &
!!$ CV-FEM matrix
                 mx_ncolm, ncolm, findm, colm, midm, &
!!$ ph matrix
                 mx_ncolph, ncolph, findph, colph, &
!!$ misc
                 mx_nface_p1 )

            !Look again for bad elements
!            if (is_porous_media) then
!              pressure_field=>extract_scalar_field(packed_state,"FEPressure")
!              allocate(Quality_list(totele*(NDIM-1)))!this number is not very well thought...
!              if (pressure_field%mesh%shape%degree < 2) &!Does not work yet for quadratic elements
!                call CheckElementAngles(packed_state, totele, x_ndgln, X_nloc, Max_bad_angle, Min_bad_angle, Quality_list)
!            end if
            call temp_mem_hacks()

!!$ Allocating space for various arrays:
            allocate( &
!!$
                 DRhoDPressure( nphase, cv_nonods ), &
!!$
                 suf_sig_diagten_bc( stotel * cv_snloc * nphase, ndim ), &
                 Mean_Pore_CV( cv_nonods ), &
                 mass_ele( totele ), &
!!$
!!$
                 Velocity_U_Source( ndim, nphase, u_nonods ), &
                 Velocity_U_Source_CV( ndim, nphase, cv_nonods ), &
                 Material_Absorption( mat_nonods, ndim * nphase, ndim * nphase ), &
                 Velocity_Absorption( mat_nonods, ndim * nphase, ndim * nphase ), &
                 Material_Absorption_Stab( mat_nonods, ndim * nphase, ndim * nphase ), &
                 ScalarField_Absorption( nphase, nphase, cv_nonods ), Component_Absorption( nphase, nphase, cv_nonods ), &
                 Temperature_Absorption( nphase, nphase, cv_nonods ), &
                 Momentum_Diffusion( mat_nonods, ndim, ndim, nphase ), &
                 Momentum_Diffusion_Vol( mat_nonods, nphase ), &
                 ScalarAdvectionField_Diffusion( mat_nonods, ndim, ndim, nphase ), &
                 Component_Diffusion( mat_nonods, ndim, ndim, nphase ), &
!!$ Variables used in the diffusion-like term: capilarity and surface tension:
                 plike_grad_sou_grad( cv_nonods * nphase ), &
                 plike_grad_sou_coef( cv_nonods * nphase ) )
!!$
            Velocity_U_Source = 0. ; Velocity_Absorption = 0. ; Velocity_U_Source_CV = 0.
            Momentum_Diffusion=0.
            Momentum_Diffusion_Vol=0.
!!$
            Temperature_Absorption=0.
!!$
            Component_Diffusion=0. ; Component_Absorption=0.
!!$

!!$
            ScalarAdvectionField_Diffusion=0. ; ScalarField_Absorption=0.
!!$
            Material_Absorption=0. ; Material_Absorption_Stab=0.
!!$
            plike_grad_sou_grad=0. ; plike_grad_sou_coef=0.
!!$
            suf_sig_diagten_bc=0.
!!$


!!$ Extracting Mesh Dependent Fields
            initialised = .true.
            call Extracting_MeshDependentFields_From_State( state, packed_state, initialised, &
                 Velocity_U_Source, Velocity_Absorption )

            ncv_faces=CV_count_faces( packed_state, CV_ELE_TYPE, stotel, cv_sndgln, u_sndgln )


!!$
!!$ Initialising Absorption terms that do not appear in the schema
!!$
            ScalarField_Absorption = 0. ; Component_Absorption = 0. ; Temperature_Absorption = 0.


!!$ Computing shape function scalars
            igot_t2 = 0 ; igot_theta_flux = 0
            if( ncomp /= 0 )then
               igot_t2 = 1 ; igot_theta_flux = 1
            end if

            call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
                 cv_ngi, cv_ngi_short, scvngi_theta, sbcvngi, nface, .false. )

            allocate( theta_flux( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), &
                 one_m_theta_flux( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), &
                 theta_flux_j( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), &
                 one_m_theta_flux_j( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), &
                 sum_theta_flux( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), &
                 sum_one_m_theta_flux( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), &
                 sum_theta_flux_j( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), &
                 sum_one_m_theta_flux_j( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), &
                 theta_gdiff( nphase, cv_nonods ), ScalarField_Source( nphase, cv_nonods ), &
                 ScalarField_Source_Store( nphase, cv_nonods ), &
                 ScalarField_Source_Component( nphase, cv_nonods ) )

            sum_theta_flux = 1. ; sum_one_m_theta_flux = 0.
            sum_theta_flux_j = 1. ; sum_one_m_theta_flux_j = 0.
            ScalarField_Source=0. ; ScalarField_Source_Store=0. ; ScalarField_Source_Component=0.

            allocate( Component_Diffusion_Operator_Coefficient( ncomp, ncomp_diff_coef, nphase ) )
            allocate(opt_vel_upwind_coefs_new(ndim, ndim, nphase, mat_nonods)); opt_vel_upwind_coefs_new =0.
            allocate(opt_vel_upwind_grad_new(ndim, ndim, nphase, mat_nonods)); opt_vel_upwind_grad_new =0.



            !!call BoundedSolutionCorrections( state, packed_state, small_finacv, small_colacv, StorageIndexes, cv_ele_type )


            call Calculate_All_Rhos( state, packed_state, ncomp, nphase, ndim, cv_nonods, cv_nloc, totele, &
                 cv_ndgln, DRhoDPressure )

         end if Conditional_ReallocatingFields

!!$ Simple adaptive time stepping algorithm
         if ( have_option( '/timestepping/adaptive_timestep' ) ) then
            c = -66.6 ; minc = 0. ; maxc = 66.e6 ; ic = 66.e6
            call get_option( '/timestepping/adaptive_timestep/requested_cfl', rc )
            call get_option( '/timestepping/adaptive_timestep/minimum_timestep', minc, stat )
            call get_option( '/timestepping/adaptive_timestep/maximum_timestep', maxc, stat )
            call get_option( '/timestepping/adaptive_timestep/increase_tolerance', ic, stat )

            do iphase = 1, nphase
               ! requested cfl
               rc_field => extract_scalar_field( state( iphase ), 'RequestedCFL', stat )
               if ( stat == 0 ) rc = min( rc, minval( rc_field % val ) )
               ! max cfl
               cfl => extract_scalar_field( state( iphase ), 'CFLNumber' )
               c = max ( c, maxval( cfl % val ) )
            end do

            call get_option( '/timestepping/timestep', dt )
            dt = max( min( min( dt * rc / c, ic * dt ), maxc ), minc )
            call allmin(dt)
            call set_option( '/timestepping/timestep', dt )
         end if

         call set_boundary_conditions_values(state, shift_time=.true.)

      end do Loop_Time

       if (has_references(metric_tensor)) call deallocate(metric_tensor)

!!$ Now deallocating arrays:
      deallocate( &
!!$ Node glabal numbers
           cv_sndgln, p_sndgln, u_sndgln, &
!!$ Sparsity patterns
           finacv, colacv, midacv,&
           small_finacv, small_colacv, small_midacv, &
           finmcy, colmcy, midmcy, &
           findgm_pha, coldgm_pha, middgm_pha, findct, &
           colct, findc, colc, findcmc, colcmc, midcmc, findm, &
           colm, midm, &
!!$ Defining element-pair type and discretisation options and coefficients
           opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
!!$ For output:
           Mean_Pore_CV, &
!!$ Variables used in the diffusion-like term: capilarity and surface tension:
           plike_grad_sou_grad, plike_grad_sou_coef, &
!!$ Working arrays
           DRhoDPressure, FEM_VOL_FRAC, &
           Velocity_U_Source, Velocity_U_Source_CV, &
           theta_gdiff, ScalarField_Source, ScalarField_Source_Store, ScalarField_Source_Component, &
           mass_ele,&
           Material_Absorption, Material_Absorption_Stab, &
           Velocity_Absorption, ScalarField_Absorption, Component_Absorption, Temperature_Absorption, &
           Component_Diffusion_Operator_Coefficient, &
           Momentum_Diffusion, Momentum_Diffusion_Vol, ScalarAdvectionField_Diffusion, &
           Component_Diffusion, &
           theta_flux, one_m_theta_flux, theta_flux_j, one_m_theta_flux_j, &
           sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j )

      ! Dump at end, unless explicitly disabled
      if(.not. have_option("/io/disable_dump_at_end")) then
         call write_state(dump_no, state)
      end if

      ! If calculating boundary fluxes, dump them to outfluxes.txt at the very end to get the values at the last time step

!      if(calculate_flux) then
!         if(getprocno() == 1) then
!
!            call dump_outflux(acctim,itime,totout,intflux, dt)
!
!         endif
!      endif


      call tag_references()

      call deallocate(packed_state)
      call deallocate(multiphase_state)
      call deallocate(multicomponent_state )

      if (allocated(Quality_list)) deallocate(Quality_list)

      if (calculate_flux) deallocate(outlet_id, totout, intflux)

      return

    contains


      subroutine temp_mem_hacks()

!!! routine puts various CSR sparsities into packed_state

        use sparse_tools

        type(csr_sparsity), pointer :: sparsity
        type(scalar_field), pointer :: sfield
        type(tensor_field), pointer :: tfield

        integer ic, stat

        type(halo_type), pointer :: halo
        type(mesh_type), pointer :: ph_mesh


        allocate(sparsity)

        sparsity=wrap(finacv,midacv,colm=colacv,name='PackedAdvectionSparsity')
        call insert(packed_state,sparsity,'PackedAdvectionSparsity')
        call deallocate(sparsity)
        sparsity=wrap(findc,colm=colc,name='CMatrixSparsity')
        call insert(packed_state,sparsity,'CMatrixSparsity')
        call deallocate(sparsity)
        sparsity=wrap(findct,colm=colct,name='CTMatrixSparsity')
        call insert(packed_state,sparsity,'CTMatrixSparsity')
        call deallocate(sparsity)

        sparsity=wrap(findm,midm,colm=colm,name='CVFEMSparsity')
        call insert(packed_state,sparsity,'CVFEMSparsity')
        call deallocate(sparsity)

        sfield=>extract_scalar_field(packed_state,"FEPressure")

        if (associated( sfield%mesh%halos)) then
           sparsity=wrap(findcmc,colm=colcmc,name='CMCSparsity',&
                row_halo=sfield%mesh%halos(2),&
                column_halo=sfield%mesh%halos(2))
        else
           sparsity=wrap(findcmc,colm=colcmc,name='CMCSparsity')
        end if
        call insert(packed_state,sparsity,'CMCSparsity')
        call deallocate(sparsity)

        if (associated( sfield%mesh%halos)) then
           sparsity=wrap(small_finacv,small_midacv,colm=small_colacv,name="ACVSparsity",&
                row_halo=sfield%mesh%halos(2),&
                column_halo=sfield%mesh%halos(2))
        else
           sparsity=wrap(small_finacv,small_midacv,colm=small_colacv,name="ACVSparsity")
        end if
        call insert(packed_state,sparsity,"ACVSparsity")
        call deallocate(sparsity)


        tfield=>extract_tensor_field(packed_state,"PackedVelocity")
        if (associated(tfield%mesh%halos)) then
           halo => tfield%mesh%halos(2)
        else
           nullify(halo)
        end if
        if (associated(halo)) then
           sparsity=wrap(findgm_pha,colm=coldgm_pha,&
                name='MomentumSparsity',row_halo=halo,column_halo=halo)
        else
           sparsity=wrap(findgm_pha,colm=coldgm_pha,name="MomentumSparsity")
        end if
        call insert(packed_state,sparsity,"MomentumSparsity")
        call deallocate(sparsity)

        sparsity=make_sparsity(sfield%mesh,sfield%mesh,&
             "PressureMassMatrixSparsity")
        call insert(packed_state,sparsity,"PressureMassMatrixSparsity")
        call deallocate(sparsity)
        deallocate(sparsity)

        sparsity=> extract_csr_sparsity(packed_state,"PressureMassMatrixSparsity")
        do ic=1,size(multicomponent_state)
           call insert(multicomponent_state(ic),sparsity,"PressureMassMatrixSparsity")
        end do
        sparsity=> extract_csr_sparsity(packed_state,"ACVSparsity")
         do ic=1,size(multicomponent_state)
           call insert(multicomponent_state(ic),sparsity,"ACVSparsity")
        end do

        sparsity=> extract_csr_sparsity(state(1),"ElementConnectivity")
        call insert(packed_state,sparsity,"ElementConnectivity")

        ph_mesh => extract_mesh( state( 1 ), "ph", stat )
        if ( stat == 0 ) then
           sparsity = wrap( findph, colm = colph, name = "phsparsity" )
           call insert( packed_state, sparsity, "phsparsity" )
           call deallocate( sparsity )
        end if


      end subroutine temp_mem_hacks

      subroutine linearise_components()

        integer :: ist,ip,ele
        type ( scalar_field ), pointer :: nfield
        integer, dimension(:), pointer :: nodes
        real, allocatable, dimension(:) :: comp

        if (ncomp>1) then
           do ist=1,size(state)
              if (has_scalar_field(state(ist),"Density")) then
                 nfield=>extract_scalar_field(state(ist),"Density")
                 allocate(comp(ele_loc(nfield,1)))
                 if (nfield%mesh%shape%degree==2) then
                    select case(mesh_dim(nfield))
                    case(1)
                       do ele=1,ele_count(nfield)
                          nodes=>ele_nodes(nfield,ele)
                          comp=ele_val(nfield,ele)
                          call set(nfield,nodes(2),0.5*(comp(1)+comp(3)))
                       end do
                    case(2)
                       do ele=1,ele_count(nfield)
                          nodes=>ele_nodes(nfield,ele)
                          comp=ele_val(nfield,ele)
                          call set(nfield,nodes(2),0.5*(comp(1)+comp(3)))
                          call set(nfield,nodes(4),0.5*(comp(1)+comp(6)))
                          call set(nfield,nodes(5),0.5*(comp(3)+comp(6)))
                       end do
                    case(3)
                       do ele=1,ele_count(nfield)
                          nodes=>ele_nodes(nfield,ele)
                          comp=ele_val(nfield,ele)
                          call set(nfield,nodes(2),0.5*(comp(1)+comp(3)))
                          call set(nfield,nodes(4),0.5*(comp(1)+comp(6)))
                          call set(nfield,nodes(5),0.5*(comp(3)+comp(6)))
                          call set(nfield,nodes(7),0.5*(comp(1)+comp(10)))
                          call set(nfield,nodes(8),0.5*(comp(3)+comp(10)))
                          call set(nfield,nodes(9),0.5*(comp(6)+comp(10)))
                       end do
                    end select
                    deallocate(comp)
                 end if
              end if
           end do
         end if

        do ist=1,size(state)
           if (has_scalar_field(state(ist),"ComponentMassFractionPhase1")) then
              do ip=1,nphase
                 nfield=>extract_scalar_field(state(ist),&
                      "ComponentMassFractionPhase"//int2str(ip))
                 allocate(comp(ele_loc(nfield,1)))
                 if (nfield%mesh%shape%degree==2) then
                    select case(mesh_dim(nfield))
                    case(1)
                       do ele=1,ele_count(nfield)
                          nodes=>ele_nodes(nfield,ele)
                          comp=ele_val(nfield,ele)
                          call set(nfield,nodes(2),0.5*(comp(1)+comp(3)))
                       end do
                    case(2)
                       do ele=1,ele_count(nfield)
                          nodes=>ele_nodes(nfield,ele)
                          comp=ele_val(nfield,ele)
                          call set(nfield,nodes(2),0.5*(comp(1)+comp(3)))
                          call set(nfield,nodes(4),0.5*(comp(1)+comp(6)))
                          call set(nfield,nodes(5),0.5*(comp(3)+comp(6)))
                       end do
                    case(3)
                       do ele=1,ele_count(nfield)
                          nodes=>ele_nodes(nfield,ele)
                          comp=ele_val(nfield,ele)
                          call set(nfield,nodes(2),0.5*(comp(1)+comp(3)))
                          call set(nfield,nodes(4),0.5*(comp(1)+comp(6)))
                          call set(nfield,nodes(5),0.5*(comp(3)+comp(6)))
                          call set(nfield,nodes(7),0.5*(comp(1)+comp(10)))
                          call set(nfield,nodes(8),0.5*(comp(3)+comp(10)))
                          call set(nfield,nodes(9),0.5*(comp(6)+comp(10)))
                       end do
                    end select
                    deallocate(comp)
                 end if
              end do
           end if
        end do


      end subroutine linearise_components

      !subroutine unlinearise_components()
      !
      !  integer :: ist, ip
      !  type( scalar_field ), pointer :: cmp, cmp_p1_cg
      !  type( mesh_type ), pointer   :: cmesh, mesh
      !  type ( scalar_field ), pointer :: nfield
      !  type ( scalar_field ) :: sfield
      !
      !end subroutine unlinearise_components

    end subroutine MultiFluids_SolveTimeLoop

    subroutine Updating_Linearised_Components( totele, ndim, cv_nloc, cv_nonods, cv_ndgln, &
         component )
      implicit none
      integer, intent( in ) :: totele, ndim, cv_nloc, cv_nonods
      integer, dimension( : ), intent( in ) :: cv_ndgln
      real, dimension ( : ), intent( inout ) :: component
!!$Local variables
      integer :: ele, cv_iloc, cv_nod
      real, dimension( : ), allocatable :: density_tmp, den_cv_nod

      allocate( density_tmp( cv_nonods ), den_cv_nod( cv_nloc ) )
      density_tmp = 0. ; den_cv_nod = 0.

      Conditional_CV_Number: if( cv_nloc == 6 .or. (cv_nloc == 10 .and. ndim==3) ) then ! P2 triangle or tet
         density_tmp = component
         Loop_Elements: do ele = 1, totele
            Loop_CV: do cv_iloc = 1, cv_nloc
               cv_nod = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )
               den_cv_nod( cv_iloc ) = density_tmp( cv_nod )
            end do Loop_CV

            den_cv_nod( 2 ) = 0.5 * ( den_cv_nod( 1 ) + den_cv_nod( 3 ) )
            den_cv_nod( 4 ) = 0.5 * ( den_cv_nod( 1 ) + den_cv_nod( 6 ) )
            den_cv_nod( 5 ) = 0.5 * ( den_cv_nod( 3 ) + den_cv_nod( 6 ) )
            if( cv_nloc == 10 ) then
               den_cv_nod( 7 ) = 0.5 * ( den_cv_nod( 1 ) + den_cv_nod( 10 ) )
               den_cv_nod( 8 ) = 0.5 * ( den_cv_nod( 3 ) + den_cv_nod( 10 ) )
               den_cv_nod( 9 ) = 0.5 * ( den_cv_nod( 6 ) + den_cv_nod( 10 ) )
            end if

            Loop_CV2: do cv_iloc = 1, cv_nloc
               cv_nod = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )
               component( cv_nod ) = den_cv_nod( cv_iloc )
            end do Loop_CV2

         end do Loop_Elements

      end if Conditional_CV_Number

      deallocate( density_tmp, den_cv_nod )

      return
    end subroutine Updating_Linearised_Components

   subroutine copy_packed_new_to_old(packed_state)
     type(state_type), intent(inout) :: packed_state

     type(scalar_field), pointer :: sfield, nsfield
     type(vector_field), pointer :: vfield, nvfield
     type(tensor_field), pointer :: tfield, ntfield

     integer :: i

     do i=1,size(packed_state%scalar_fields)
        sfield=>packed_state%scalar_fields(i)%ptr
        if (sfield%name(1:9)=="PackedOld") then
           nsfield=>extract_scalar_field(packed_state,"Packed"//sfield%name(10:))
           sfield%val=nsfield%val
        end if
     end do

     do i=1,size(packed_state%vector_fields)
        vfield=>packed_state%vector_fields(i)%ptr
        if (vfield%name(1:9)=="PackedOld") then
           nvfield=>extract_vector_field(packed_state,"Packed"//vfield%name(10:))
           vfield%val=nvfield%val
        end if
     end do

     do i=1,size(packed_state%tensor_fields)
        tfield=>packed_state%tensor_fields(i)%ptr
        if (tfield%name(1:9)=="PackedOld") then
           ntfield=>extract_tensor_field(packed_state,"Packed"//tfield%name(10:))
           tfield%val=ntfield%val
        end if
     end do

     sfield=>extract_scalar_field(packed_state,"OldFEPressure")
     nsfield=>extract_scalar_field(packed_state,"FEPressure")
     sfield%val=nsfield%val

     sfield=>extract_scalar_field(packed_state,"OldCVPressure")
     nsfield=>extract_scalar_field(packed_state,"CVPressure")
     sfield%val=nsfield%val

   end subroutine copy_packed_new_to_old

   subroutine set_nu_to_u(packed_state)
     type(state_type), intent(inout) :: packed_state
     type(tensor_field), pointer :: u, uold, nu, nuold

         u  => extract_tensor_field( packed_state, "PackedVelocity" )
         uold  => extract_tensor_field( packed_state, "PackedOldVelocity" )

         nu => extract_tensor_field( packed_state, "PackedNonlinearVelocity" )
         nuold => extract_tensor_field( packed_state, "PackedOldNonlinearVelocity" )

         nu % val = u % val
         nuold % val = uold % val

   end subroutine set_nu_to_u



   subroutine linearise( field )
     implicit none
     type( tensor_field ), intent( inout ) :: field

     integer, dimension( : ), pointer :: ndglno, cv_nods
     integer :: n, totele, cv_nloc, ncomp, nphase, cv_nonods, ele, cv_iloc, cv_nod
     real, dimension( :, :, : ), allocatable :: field_tmp, field_cv_nod

     ! This sub will linearise a p2 field

     n = field%mesh%shape%degree

     if ( n==2 ) then

        ndglno => get_ndglno( field%mesh )

        totele = field%mesh%elements
        cv_nloc = field%mesh%shape%loc
        cv_nonods =  field%mesh%nodes

        ncomp = size( field%val, 1 )
        nphase = size( field%val, 2 )

        allocate( field_tmp( ncomp, nphase, cv_nonods ) ) ; field_tmp = field % val
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
              field%val( :, :, cv_nod ) = field_cv_nod( :, :, cv_iloc )
           end do

        end do

        deallocate( field_tmp, field_cv_nod )

     end if

     return
   end subroutine linearise

   subroutine dump_outflux(current_time, itime, outflux, intflux,ts)

   ! Subroutine that dumps the total flux at a given timestep across all specified boudaries to a file  called 'outfluxes.txt'. In addition, the time integrated flux
   ! up to the current timestep is also outputted to this file. Integration boundaries are specified in diamond via surface_ids.
   ! (In diamond this option can be found under "/io/dump_boundaryflux/surface_ids" and the user should specify an integer array containing the IDs of every boundary they
   !to integrate over).

   real,intent(in) :: current_time
   integer, intent(in) :: itime
   real, dimension(:,:), intent(inout) :: outflux, intflux
   real,intent(in) :: ts


   integer :: ioutlet
   type(stat_type), target :: default_stat
   character (len=1000000) :: whole_line
   character (len=1000) :: numbers


   default_stat%conv_unit=free_unit()

   open(unit=default_stat%conv_unit, file="outfluxes.txt", action="write", position="append")

   ! Write column headings to file
   !   if(itime.eq.1) then
   !       whole_line = "Current Time"
   !       do ioutlet =1, size(outflux,2)
   !           write(numbers,*), outlet_id(ioutlet)
   !           whole_line = trim(whole_line) // "Surface id:"// trim(numbers)// "Phase1 boundary flux"// &
   !           "Phase2 boundary flux"// "Phase1 integrated flux"// "Phase2 integrated flux "
   !       end do
   !        !Write a line
   !       write(default_stat%conv_unit,*), trim(whole_line)
   !   endif

   write(whole_line,*), current_time
   do ioutlet =1, size(outflux,2)
       write(numbers,*), outflux(:, ioutlet),  intflux(:, ioutlet)
       whole_line = trim(whole_line) //" "// trim(numbers)//" "
   end do
    !Write a line
   write(default_stat%conv_unit,*), trim(whole_line)


!   if(itime.eq.1) then
!       write(default_stat%conv_unit,*) "Current Time"
!       do ioutlet =1, size(outflux,2)
!           write(default_stat%conv_unit,*) "Surface id:",outlet_id(ioutlet),"""", &
!           "Phase1 boundary flux","Phase2 boundary flux","Phase1 integrated flux","Phase2 integrated flux"
!       end do
!
!   endif


!write(default_stat%conv_unit,*) current_time
!
!
!do ioutlet =1, size(outflux,2)
!  write(default_stat%conv_unit,*) outflux(:, ioutlet),  intflux(:, ioutlet)
!end do


   close (default_stat%conv_unit)


   end subroutine dump_outflux



   end module multiphase_time_loop
