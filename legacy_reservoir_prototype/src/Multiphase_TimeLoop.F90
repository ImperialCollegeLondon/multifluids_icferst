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
    use diagnostic_fields_new_multiphase, only : &
        calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
        check_diagnostic_dependencies
    use global_parameters, only: timestep, simulation_start_time, simulation_start_cpu_time, &
        simulation_start_wall_time, new_mesh, &
        topology_mesh_name, current_time, is_porous_media, after_adapt, is_multifracture, is_first_time_step, &
        OPTION_PATH_LEN, FIELD_NAME_LEN
    use fldebug
    use reference_counting
    use state_module
    use fields
    use field_options
    use fields_allocates
    use field_priority_lists
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
    use fluids_module, only: pre_adapt_tasks, update_state_post_adapt
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
    use multiphase_caching
    use Compositional_Terms
    use Copy_Outof_State
    use cv_advection, only : cv_count_faces
    use multiphase_1D_engine
    use boundary_conditions_from_options
    use multi_data_types
    use vtk_interfaces
    use multi_interpolation
    use multi_surface_tension
#ifdef HAVE_ZOLTAN
  use zoltan
#endif
    !use matrix_operations
    !use shape_functions
    implicit none
    private
    !public :: MultiFluids_SolveTimeLoop, rheology, dump_outflux
    public :: MultiFluids_SolveTimeLoop, dump_outflux
    !type(rheology_type), dimension(:), allocatable :: rheology
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
        !!$ additional state variable for storage
        type(state_type) :: storage_state
        !!Define shape functions
        type (multi_shape_funs) :: CV_funs, FE_funs
        !!$ Primary scalars
        type(multi_dimensions) :: Mdims
        type(multi_gi_dimensions) :: CV_GIdims, FE_GIdims
        !sprint_to_do !substitute all these instances by the structure Mdims
        integer :: ntsol
        !!$ Node global numbers
        type(multi_ndgln) :: ndgln
        !!$ Sparsity patterns
        type (multi_sparsities) :: Mspars
        !!$ Defining element-pair type and discretisation options and coefficients
        type (multi_discretization_opts) :: Mdisopt
        !!$ Defining the necessary matrices and corresponding RHS
        type (multi_matrices) :: Mmat
        !sprint_to_do!remove all the is store inside Mspars when Mspars is fully implemented
        integer :: nlenmcy, mx_nface_p1, mx_ncolacv, mxnele, mx_ncoldgm_pha, &
            mx_ncolmcy, mx_nct, mx_nc, mx_ncolcmc, mx_ncolm, mx_ncolph
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
        integer :: scvngi_theta, igot_t2, igot_theta_flux, IGOT_THERM_VIS
        !!$ CV-wise porosity
        real, dimension( :, : ), allocatable :: Mean_Pore_CV
        !!$ Variables used in the diffusion-like term: capilarity and surface tension:
        integer :: iplike_grad_sou
        real, dimension( :, :, : ), allocatable :: plike_grad_sou_grad, plike_grad_sou_coef
        !!$ Adaptivity related fields and options:
        type( tensor_field ) :: metric_tensor
        type( state_type ), dimension( : ), pointer :: sub_state => null()
        integer :: nonlinear_iterations_adapt
        logical :: do_reallocate_fields = .false., not_to_move_det_yet = .false.
        !!$ Working arrays:
        real, dimension(:), pointer :: mass_ele
        real, dimension( :, : ), pointer :: THETA_GDIFF
        !!$
        real, dimension( :, : ), pointer :: &
            ScalarField_Source, ScalarField_Source_Store, ScalarField_Source_Component
        real, dimension( :, :, : ), allocatable :: Material_Absorption, &
            Velocity_Absorption, ScalarField_Absorption, Component_Absorption, Temperature_Absorption
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
        type( tensor_field ), pointer :: D_s, DC_s, DCOLD_s
        type( tensor_field ), pointer :: MFC_s, MFCOLD_s
        type( mesh_type ), pointer :: ph_mesh
        !! face value storage
        integer :: ncv_faces
        !Courant number for porous media
        real :: Courant_number = -1
        !Variables for adapting the mesh within the FPI solver
        logical :: adapt_mesh_in_FPI
        !type( scalar_field ) :: Saturation_bak, ConvSats
        integer :: checkpoint_number
        !Array to map nodes to region ids
        integer, dimension(:), allocatable :: IDs_ndgln, IDs2CV_ndgln!sprint_to_do; get this into ndgln structure?
        !Variable to store where we store things. Do not oversize this array, the size has to be the last index in use
        integer, dimension (38) :: StorageIndexes
        !Distribution of the indexes of StorageIndexes:
        !cv_fem_shape_funs_plus_storage: 1 (ASSEMB_FORCE_CTY), 13 (CV_ASSEMB)   (REMOVED)
        !CALC_ANISOTROP_LIM            : 2 (REMOVED)
        !DG_DERIVS_ALL2                : 3 (REMOVED)
        !DETNLXR_INVJAC                : 4 (removed)
        !UNPACK_LOC                    : 5,6,7,8,9,10(REMOVED)
        !COLOR_GET_CMC_PHA             : 11 (REMOVED)
        !Matrix C                      : 12 (REMOVED)
        !DG_DERIVS_ALL                 : 14 (REMOVED)
        !DETNLXR_PLUS_U_WITH_STORAGE   : 14
        !Indexes used in SURFACE_TENSION_WRAPPER (deprecated and will be removed):[15,30] (REMOVED)
        !PROJ_CV_TO_FEM_state          : 31 (disabled) (removed)
        !Capillary pressure            : 32 (Pe), 33 (exponent a) (disabled) (removed?)
        !PIVIT_MAT (inverted)          : 34 (REMOVED)
        !Bound                         : 35 (removed?)
        !Ph 1                          : 36 (added back, need to remove again)
        !Ph 2                          : 37 (removed?)
        !Working pointers
        type(tensor_field), pointer :: tracer_field, velocity_field, density_field, saturation_field, old_saturation_field, tracer_source
        type(tensor_field), pointer :: pressure_field, cv_pressure, fe_pressure
        type(scalar_field), pointer :: f1, f2
        type(vector_field), pointer :: positions, porosity_field
        logical, parameter :: write_all_stats=.true.
        ! Variables used for calculating boundary outfluxes. Logical "calculate_flux" determines if this calculation is done. Intflux is the time integrated outflux
        ! Ioutlet counts the number of boundaries over which to calculate the outflux
        integer :: ioutlet
        real, dimension(:,:),  allocatable  :: intflux
        logical :: calculate_flux
        ! Variables used in the CVGalerkin interpolation calculation
        integer :: numberfields
        logical :: new_ntsol_loop
        character(len=OPTION_PATH_LEN) :: option_buffer
        character(len=FIELD_NAME_LEN) :: tmp_name
        logical :: use_advdif, multiphase_scalar
        integer :: it, it2, nphase_scalar
#ifdef HAVE_ZOLTAN
      real(zoltan_float) :: ver
      integer(zoltan_int) :: ierr
      ierr = Zoltan_Initialize(ver)
      assert(ierr == ZOLTAN_OK)
#endif
        !Initially we set to use Stored data and that we have a new mesh
        StorageIndexes = 0!Initialize them as zero !
        ! Number of pressures to solve for
        Mdims%npres = option_count("/material_phase/scalar_field::Pressure/prognostic")
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
        call pack_multistate( Mdims%npres, state, packed_state, multiphase_state, &
            multicomponent_state )
        call set_boundary_conditions_values(state, shift_time=.true.)
        !Prepare the caching
        call set_caching_level()
        call initialize_storage( packed_state, storage_state )
        !call print_state( packed_state )
        !stop 78
        !  Access boundary conditions via a call like
        !  call get_entire_boundary_condition(extract_tensor_field(packed_state,"Packed"//name),["weakdirichlet"],tfield,bc_type_list)
        !  where tfield is type(tensor_field) and bc_type_list is integer, dimension(tfield%dim(1),tfield%dim(2),nonods)
        !  Then values are in tfield%val(1/Mdims%ndim/Mdims%ncomp,Mdims%nphase,nonods)
        !  Type ids are in bc_type_list(1/Mdims%ndim/Mdims%ncomp,Mdims%nphase,Mdims%stotel)
        !
        !  A deallocate tfield when finished!!
        Repeat_time_step = .false.!Initially has to be false
        nonLinearAdaptTs = have_option(  '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear')
        !!$ Compute primary scalars used in most of the code
        call Get_Primary_Scalars_new( state, Mdims )
        !!$ Calculating Global Node Numbers
        call allocate_multi_ndgln(ndgln, Mdims)
        call Compute_Node_Global_Numbers(state, ndgln)
        !!$
        !!$ Computing Sparsity Patterns Matrices
        !!$
        !!$ Defining lengths and allocating space for the matrices


        ! move somewhere else - SPRINT TO DO -
        ph_mesh => extract_mesh( state( 1 ), "ph", stat )
        if ( stat == 0 ) then
           Mdims%ph_nonods = node_count( ph_mesh )
           Mdims%ph_nloc = ele_loc( ph_mesh, 1 )
        end if

        call Defining_MaxLengths_for_Sparsity_Matrices( Mdims%ndim, Mdims%nphase, Mdims%totele, Mdims%u_nloc, Mdims%cv_nloc, Mdims%ph_nloc, Mdims%cv_nonods, &
            mx_nface_p1, mxnele, mx_nct, mx_nc, mx_ncolcmc, mx_ncoldgm_pha, mx_ncolmcy, &
            mx_ncolacv, mx_ncolm, mx_ncolph )
        nlenmcy = Mdims%u_nonods * Mdims%nphase * Mdims%ndim + Mdims%cv_nonods
        !!$ Defining element-pair type
        call Get_Ele_Type_new( Mdims, Mdisopt )
        !Allocate and calculate the sparsity patterns matrices
        call Get_Sparsity_Patterns( state, Mdims, Mspars, ndgln, Mdisopt, mx_ncolacv, nlenmcy, mx_ncolmcy, &
                mx_ncoldgm_pha, mx_nct,mx_nc, mx_ncolcmc, mx_ncolm, mx_ncolph, mx_nface_p1 )
        call temp_mem_hacks()
        !!$ Allocating space for various arrays:
        allocate( &
            !!$
            suf_sig_diagten_bc( Mdims%stotel * Mdims%cv_snloc * Mdims%nphase, Mdims%ndim ), &
            Mean_Pore_CV( Mdims%npres, Mdims%cv_nonods ), &
            mass_ele( Mdims%totele ), &
            !!$
            Material_Absorption( Mdims%mat_nonods, Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase ), &
            ScalarField_Absorption( Mdims%nphase, Mdims%nphase, Mdims%cv_nonods ), Component_Absorption( Mdims%nphase, Mdims%nphase, Mdims%cv_nonods ), & ! fix me..move in intenerg
            !!$ Variables used in the diffusion-like term: capilarity and surface tension:
            plike_grad_sou_grad( Mdims%ncomp, Mdims%nphase, Mdims%cv_nonods ), &
            plike_grad_sou_coef( Mdims%ncomp, Mdims%nphase, Mdims%cv_nonods ), &
            )
        !!$
        suf_sig_diagten_bc=0.
        !!$
        Mean_Pore_CV=0.
        mass_ele=0.
        !!$
        Material_Absorption=0.
        ScalarField_Absorption=0. ; Component_Absorption=0.
        !!$
        plike_grad_sou_grad=0.
        plike_grad_sou_coef=0.
        iplike_grad_sou=0
        do iphase = 1, Mdims%nphase
            f1 => extract_scalar_field( state(iphase), "Temperature", stat )
            if ( stat==0 ) then
                f2 => extract_scalar_field( state(iphase),"Dummy", stat )
                if ( stat==0 ) f2%val = f1%val
            end if
        end do
        !!$ Calculate diagnostic fields
        call calculate_diagnostic_variables( state, exclude_nonrecalculated = .true. )
        call calculate_diagnostic_variables_new( state, exclude_nonrecalculated = .true. )
        !!$
        !!$ Initialising Absorption terms that do not appear in the schema
        !!$
        ScalarField_Absorption = 0. ; Component_Absorption = 0.
        !!$ Computing shape function scalars
        igot_t2 = 0 ; igot_theta_flux = 0
        if( Mdims%ncomp /= 0 )then
            igot_t2 = 1 ; igot_theta_flux = 1
        end if
        !Calculate the gauss integer numbers
        call retrieve_ngi( CV_GIdims, Mdims, Mdisopt%cv_ele_type, quad_over_whole_ele = .false. )
        call retrieve_ngi( FE_GIdims, Mdims, Mdisopt%u_ele_type, quad_over_whole_ele = .true. )
        !! Compute reference shape functions
        call allocate_multi_shape_funs( CV_funs, Mdims, CV_GIdims )
        call allocate_multi_shape_funs( FE_funs, Mdims, FE_GIdims )
        call cv_fem_shape_funs_new( CV_funs, Mdims, CV_GIdims, Mdisopt%cv_ele_type, quad_over_whole_ele = .false. )
        call cv_fem_shape_funs_new( FE_funs, Mdims, FE_GIdims, Mdisopt%u_ele_type, quad_over_whole_ele = .true. )
        !Obtain the number of faces in the control volume space
        ncv_faces=CV_count_faces( Mdims, Mdisopt%cv_ele_type, CV_GIDIMS = CV_GIdims)
        allocate( theta_flux( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            one_m_theta_flux( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            theta_flux_j( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            one_m_theta_flux_j( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            sum_theta_flux( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            sum_one_m_theta_flux( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            sum_theta_flux_j( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            sum_one_m_theta_flux_j( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            theta_gdiff( Mdims%nphase, Mdims%cv_nonods ), ScalarField_Source( Mdims%nphase, Mdims%cv_nonods ), &
            ScalarField_Source_Store( Mdims%nphase, Mdims%cv_nonods ), &
            ScalarField_Source_Component( Mdims%nphase, Mdims%cv_nonods ) )
        sum_theta_flux = 1. ; sum_one_m_theta_flux = 0.
        sum_theta_flux_j = 1. ; sum_one_m_theta_flux_j = 0.
        ScalarField_Source=0. ; ScalarField_Source_Store=0. ; ScalarField_Source_Component=0.
        !!$ Defining discretisation options
        call Get_Discretisation_Options( state, Mdims, Mdisopt )
        !!$ Option not currently set up in the schema and zeroed from the begining. It is used to control
        !!$ the upwinding rate (in the absorption term) during advection/assembling.
        allocate(opt_vel_upwind_coefs_new(Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%mat_nonods)); opt_vel_upwind_coefs_new =0.
        allocate(opt_vel_upwind_grad_new(Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%mat_nonods)); opt_vel_upwind_grad_new =0.
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
        !      call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration', tolerance_between_non_linear, default = -1. )
        !!$
        have_temperature_field = .false. ; have_component_field = .false. ; have_extra_DiffusionLikeTerm = .false.
        do istate = 1, Mdims%nstate
            if( have_option( '/material_phase[' // int2str( istate - 1 ) // ']/scalar_field::Temperature' ) ) &
                have_temperature_field = .true.
            if( have_option( '/material_phase[' // int2str( istate - 1 ) // ']/is_multiphase_component' ) ) &
                have_component_field = .true.
            !!$
            call Calculate_All_Rhos( state, packed_state, Mdims )
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
        adapt_mesh_in_FPI = have_option( '/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI')
        if (adapt_mesh_in_FPI) then !Complementary part of the code in Populate_State.F90. Line: 1290
            !Prepare state by adding two extra fields for backups (TEMPORARY USE OF DENSITY OPTION_PATH)
            do iphase = 1, Mdims%nphase!Set name and option_path (this latter to impose projection method)
                call allocate_and_insert_scalar_field('/material_phase['//int2str(iphase-1)//']/scalar_field::Density', &
                    state(iphase), field_name = "Saturation_bak", parent_mesh = "PressureMesh")
            end do
        end if

        if (is_porous_media) then
            !Get into packed state relative permeability, immobile fractions, ...
            call get_RockFluidProp(state, packed_state)
            !Convert material properties to be stored using region ids, only if porous media
            call get_regionIDs2nodes(state, packed_state, ndgln%cv, IDs_ndgln, IDs2CV_ndgln, &
                fake_IDs_ndgln = .not. is_porous_media)! .or. is_multifracture )
        end if
        !!$ Starting Time Loop
        itime = 0
        dtime = 0
        if( &
             ! if this is not a zero timestep simulation (otherwise, there would
             ! be two identical dump files)
             current_time < finish_time &! unless explicitly disabled
             .and. .not. have_option("/io/disable_dump_at_start") &
             ) then
            call write_state(dump_no, state)
        end if
        if(have_option("/io/stat/output_at_start")) then
            call write_diagnostics(state, current_time, dt,&
                timestep, not_to_move_det_yet=.true.)
        end if
        ! When outlet_id is allocated, calculate_flux is true and we want to calculate outfluxes
        calculate_flux = allocated(outlet_id)
        ! If calculating boundary fluxes, initialise to zero time integrated fluxes (intflux) and the quantity (totout) used to calculate them.
        if(calculate_flux) then
            allocate(intflux(Mdims%nphase,size(outlet_id)))
            allocate(totout(Mdims%nphase, size(outlet_id)))
            do ioutlet = 1, size(outlet_id)
                intflux(:, ioutlet) = 0.
                totout(:, ioutlet) = 0.
            enddo
        endif
        checkpoint_number=1
!!$ Time loop
        Loop_Time: do
            ewrite(2,*) '    NEW DT', itime+1
            sum_theta_flux_j = 1. ; sum_one_m_theta_flux_j = 0.
            if ( do_checkpoint_simulation( dtime ) ) then
               call checkpoint_simulation( state, cp_no=checkpoint_number, &
                    protect_simulation_name=.true., file_type='.mpml' )
               checkpoint_number=checkpoint_number+1
            end if
            dtime = dtime + 1
            ! Adapt mesh within the FPI?
            adapt_mesh_in_FPI = have_option( '/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI')
            itime = itime + 1
            timestep = itime
            call get_option( '/timestepping/timestep', dt )
            acctim = acctim + dt
            call set_option( '/timestepping/current_time', acctim )
            new_lim = .true.
            ! Added a tolerance of 0.001dt to the condition below that stops us exiting the loop before printing the last time step.
            if ( calculate_flux ) then
               if ( acctim > finish_time + 0.001*dt ) then
                  ewrite(1,*) "Passed final time"
                  exit Loop_Time
               end if
            else
               if ( acctim > finish_time ) then
                  ewrite(1,*) "Passed final time"
                  exit Loop_Time
               end if
            end if
            call get_option( '/timestepping/final_timestep', final_timestep, stat )
            if ( stat == spud_no_error ) then
                if ( itime > final_timestep ) then
                    ewrite(1,*) "Passed final timestep"
                    exit Loop_Time
                end if
            end if
            ExitNonLinearLoop = .false.
            !Store backup to be able to repeat a timestep
            if (nonLinearAdaptTs) call Adaptive_NonLinear(packed_state, reference_field, its, &
                Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,1)
            porosity_field=>extract_vector_field(packed_state,"Porosity")
            ! evaluate prescribed fields at time = current_time+dt
            call set_prescribed_field_values( state, exclude_interpolated = .true., &
                exclude_nonreprescribed = .true., time = acctim )
            !! Update all fields from time-step 'N - 1'
            call copy_packed_new_to_old( packed_state )
            !!$ FEMDEM...
#ifdef USING_FEMDEM
            if ( is_multifracture ) then
               call fracking(packed_state, state,Mdims%nphase)
            elseif ( have_option( '/blasting') ) then
               call blasting( packed_state, Mdims%nphase )
               call update_blasting_memory( packed_state, state, timestep )
            end if
#endif
            !!$ Start non-linear loop
            its = 1
            Loop_NonLinearIteration: do  while (its <= NonLinearIteration)
                ewrite(2,*) '  NEW ITS', its
                ! open the boiling test for two phases-gas and liquid
                if (have_option('/boiling') ) then
                   call set_nu_to_u( packed_state )
                   allocate ( Velocity_Absorption( Mdims%mat_nonods, Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase ), &
                              Temperature_Absorption( Mdims%nphase, Mdims%nphase, Mdims%cv_nonods ) )
                   call boiling( state, packed_state, Mdims%cv_nonods, Mdims%mat_nonods, Mdims%nphase, Mdims%ndim, &
                      ScalarField_Source, velocity_absorption, temperature_absorption )
                   deallocate ( Velocity_Absorption, temperature_absorption )
                end if
                !Store the field we want to compare with to check how are the computations going
                call Adaptive_NonLinear(packed_state, reference_field, its, &
                    Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,2)
                call Calculate_All_Rhos( state, packed_state, Mdims )
                if( solve_force_balance .and. is_porous_media ) then
                    call Calculate_PorousMedia_AbsorptionTerms( state, packed_state, Mdims, CV_GIdims, Mspars, &
                       Material_Absorption, suf_sig_diagten_bc, opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
                       ids_ndgln, IDs2CV_ndgln, ndgln%cv, ndgln%suf_cv, ndgln%mat, ndgln%x, Mdisopt%cv_ele_type )
                end if


                !!$ Solve advection of the scalar 'Temperature':
if (.true.) then
                Conditional_ScalarAdvectionField: if( have_temperature_field .and. &
                    have_option( '/material_phase[0]/scalar_field::Temperature/prognostic' ) ) then
                    ewrite(3,*)'Now advecting Temperature Field'
                    call set_nu_to_u( packed_state )
                    !call calculate_diffusivity( state, Mdims%ncomp, Mdims%nphase, Mdims%ndim, Mdims%cv_nonods, Mdims%mat_nonods, &
                    !    Mdims%mat_nloc, Mdims%totele, ndgln%mat, ScalarAdvectionField_Diffusion )
                    tracer_field=>extract_tensor_field(packed_state,"PackedTemperature")
                    velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
                    density_field=>extract_tensor_field(packed_state,"PackedDensity",stat)
                    saturation_field=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
                    call INTENERGE_ASSEM_SOLVE( state, packed_state, &
                        Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,&
                        tracer_field,velocity_field,density_field, dt, &
                        suf_sig_diagten_bc, &
                        Porosity_field%val, &
                        !!$
                        opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
                        0, igot_theta_flux, &
                        Mdisopt%t_get_theta_flux, Mdisopt%t_use_theta_flux, &
                        THETA_GDIFF, &
                        Mean_Pore_CV, &
                        option_path = '/material_phase[0]/scalar_field::Temperature', &
                        thermal = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/equation::InternalEnergy'),&
                        saturation=saturation_field, IDs_ndgln=IDs_ndgln )
                    call Calculate_All_Rhos( state, packed_state, Mdims )
                end if Conditional_ScalarAdvectionField
end if








!!$ Solve advection of the scalars.   'Temperature':

!!$ Fields...
!!-
        new_ntsol_loop = .false.

if ( new_ntsol_loop  ) then

        call get_ntsol( ntsol )
        call initialise_field_lists_from_options( state, ntsol )


        call set_nu_to_u( packed_state )
        !call calculate_diffusivity( state, Mdims%ncomp, Mdims%nphase, Mdims%ndim, Mdims%cv_nonods, Mdims%mat_nonods, &
                    !    Mdims%mat_nloc, Mdims%totele, ndgln%mat, ScalarAdvectionField_Diffusion )
        velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
        density_field=>extract_tensor_field(packed_state,"PackedDensity",stat)
        saturation_field=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")





        do it = 1, ntsol

           call get_option( trim( field_optionpath_list( it ) ) // &
                '/prognostic/equation[0]/name', &
                option_buffer, default = "UnknownEquationType" )
           select case( trim( option_buffer ) )
           case ( "AdvectionDiffusion", "InternalEnergy" )
              use_advdif = .true.
           case default
              use_advdif = .false.
           end select

           !use_advdif=.true.


           if ( use_advdif ) then

              ! figure out if scalar field is mutli-phase
              multiphase_scalar = .false.
              do it2 = it+1, ntsol
                 if ( field_name_list( it ) == field_name_list( it2 ) ) then
                    multiphase_scalar = .true.
                 end if
              end do

              tmp_name = "Packed" //field_name_list( it )
              nphase_scalar = 1
              if ( multiphase_scalar ) then
                 nphase_scalar = Mdims%nphase
                 tmp_name = "Packed" // field_name_list( it )
              end if
              tracer_field => extract_tensor_field( packed_state, trim( tmp_name ) )


              if (field_name_list( it)== 'PhaseVolumeFraction' .or.  field_name_list( it)== 'ComponentMassFractionPhase[0]') then
                    cycle
              elseif (multiphase_scalar) then

                    call INTENERGE_ASSEM_SOLVE( state, packed_state, &
                        Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,&
                        tracer_field,velocity_field,density_field, dt, &
                        suf_sig_diagten_bc, &
                        Porosity_field%val, &
                        !!$
                        opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
                        0, igot_theta_flux, &
                        Mdisopt%t_get_theta_flux, Mdisopt%t_use_theta_flux, &
                        THETA_GDIFF, &
                        Mean_Pore_CV, &
                        option_path = '/material_phase[0]/scalar_field::Temperature', &
                        thermal = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/equation::InternalEnergy'),&
                        saturation=saturation_field, IDs_ndgln=IDs_ndgln )
                    call Calculate_All_Rhos( state, packed_state, Mdims )

                    exit
              else
                    
                    call INTENERGE_ASSEM_SOLVE( state, packed_state, &
                        Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,&
                        tracer_field,velocity_field,density_field, dt, &
                        suf_sig_diagten_bc,  Porosity_field%val, &
                        opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
                        0, igot_theta_flux, &
                        Mdisopt%t_get_theta_flux, Mdisopt%t_use_theta_flux, &
                        THETA_GDIFF,  Mean_Pore_CV, &
                        option_path = '/material_phase[0]/scalar_field::Temperature', &
                        thermal = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/equation::InternalEnergy'),&
                        saturation=saturation_field, IDs_ndgln=IDs_ndgln )
                    call Calculate_All_Rhos( state, packed_state, Mdims )
               end if
                    

           end if

        end do

end if
                ScalarField_Source_Store = ScalarField_Source + ScalarField_Source_Component
                tracer_source => extract_tensor_field(packed_state,"PackedPhaseVolumeFractionSource")
                Mdisopt%volfra_use_theta_flux = .true.
                if( Mdims%ncomp <= 1 ) Mdisopt%volfra_use_theta_flux = .false.
                !!$ Now solving the Momentum Equation ( = Force Balance Equation )
                Conditional_ForceBalanceEquation: if ( solve_force_balance ) then
                    call set_nu_to_u( packed_state )
                    !!$ Diffusion-like term -- here used as part of the capillary pressure for porous media. It can also be
                    !!$ extended to surface tension -like term.
                    iplike_grad_sou = 0
                    plike_grad_sou_grad = 0
                    CALL CALCULATE_SURFACE_TENSION_NEW( state, packed_state, Mdims, Mspars, ndgln, Mdisopt, Mdims%nphase, Mdims%ncomp, &
                        PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, IPLIKE_GRAD_SOU, &
                        Mspars%ACV%ncol, Mspars%ACV%fin, Mspars%ACV%col, Mspars%ACV%mid, &
                        Mspars%small_acv%fin, Mspars%small_acv%col, Mspars%small_acv%mid, &
                        Mspars%CT%ncol, Mspars%CT%fin, Mspars%CT%col, &
                        Mdims%cv_nonods, Mdims%u_nonods, Mdims%x_nonods, Mdims%totele, Mdims%stotel, &
                        Mdisopt%cv_ele_type, Mdisopt%cv_sele_type, Mdisopt%u_ele_type, &
                        Mdims%cv_nloc, Mdims%u_nloc, Mdims%x_nloc, Mdims%cv_snloc, Mdims%u_snloc, &
                        ndgln%cv, ndgln%suf_cv, ndgln%x, ndgln%u, ndgln%suf_u, &
                        Mdims%mat_nloc, ndgln%mat, Mdims%mat_nonods,  &
                        Mdims%ndim,  &
                        Mspars%M%ncol, Mspars%M%fin, Mspars%M%col, Mspars%M%mid, &
                        Mdims%xu_nloc, ndgln%xu, Mspars%ELE%fin, Mspars%ELE%col, Mspars%ELE%ncol)
                    if( have_option_for_any_phase( '/multiphase_properties/capillary_pressure', Mdims%nphase ) )then
                                !The first time (itime/=1 .or. its/=1) we use CVSat since FESAt is not defined yet
                        call calculate_capillary_pressure(packed_state, .false., &
                            ndgln%cv, ids_ndgln, Mdims%totele, Mdims%cv_nloc)
                    end if
                    velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
                    pressure_field=>extract_tensor_field(packed_state,"PackedFEPressure")
                    CALL FORCE_BAL_CTY_ASSEM_SOLVE( state, packed_state, &
                        Mdims, CV_GIdims, FE_GIdims, CV_funs, FE_funs, Mspars, ndgln, Mdisopt, Mmat, storage_state,&
                        velocity_field, pressure_field, &
                        Material_Absorption, dt, NLENMCY, & ! Force balance plus cty multi-phase eqns
                        SUF_SIG_DIAGTEN_BC, &
                        ScalarField_Source_Store, ScalarField_Absorption, Porosity_field%val, &
                        opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
                        igot_theta_flux, &
                        sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j, &
                        iplike_grad_sou, plike_grad_sou_coef, plike_grad_sou_grad, &
                        StorageIndexes=StorageIndexes, IDs_ndgln=IDs_ndgln )
                    velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
                    !                    !Calculate actual Darcy velocity
                    !                    call get_DarcyVelocity(Mdims%totele, Mdims%cv_nloc, Mdims%u_nloc, Mdims%mat_nloc, ndgln%mat, ndgln%u, &
                    !                            ndgln%cv, state, packed_state, Material_Absorption)
                    !!$ Calculate Density_Component for compositional
                    if ( have_component_field ) call Calculate_Component_Rho( state, packed_state, Mdims )
                end if Conditional_ForceBalanceEquation
                Conditional_PhaseVolumeFraction: if ( solve_PhaseVolumeFraction ) then
                    call VolumeFraction_Assemble_Solve( state, packed_state, &
                        Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,&
                        dt, SUF_SIG_DIAGTEN_BC, &
                        ScalarField_Source_Store, ScalarField_Absorption, Porosity_field%val, &
                        opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
                        igot_theta_flux, &
                        option_path = '/material_phase[0]/scalar_field::PhaseVolumeFraction', &
                        mass_ele_transp = mass_ele,&
                        theta_flux=sum_theta_flux, one_m_theta_flux=sum_one_m_theta_flux, &
                        theta_flux_j=sum_theta_flux_j, one_m_theta_flux_j=sum_one_m_theta_flux_j,&
                        Material_Absorption=Material_Absorption,&
                        nonlinear_iteration = its, IDs_ndgln = IDs_ndgln, IDs2CV_ndgln = IDs2CV_ndgln, &
                        Courant_number = Courant_number)
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
                    Loop_Components: do icomp = 1, Mdims%ncomp
                        tracer_field=>extract_tensor_field(multicomponent_state(icomp),"PackedComponentMassFraction")
                        density_field=>extract_tensor_field(multicomponent_state(icomp),"PackedComponentDensity",stat)
                        !!$ Computing the absorption term for the multi-components equation
                        call Calculate_ComponentAbsorptionTerm( state, packed_state, &
                            icomp, ndgln%cv, Mdims, &
                            D_s%val, Porosity_field%val, mass_ele, &
                            Component_Absorption, IDs_ndgln )
                        Conditional_SmoothAbsorption: if( have_option( '/material_phase[' // int2str( Mdims%nstate - Mdims%ncomp ) // &
                            ']/is_multiphase_component/KComp_Sigmoid' ) .and. Mdims%nphase > 1 ) then
                            do cv_nodi = 1, Mdims%cv_nonods
                                if(saturation_field%val(1,1, cv_nodi ) > 0.95 ) then
                                    do iphase = 1, Mdims%nphase
                                        do jphase = min( iphase + 1, Mdims%nphase ), Mdims%nphase
                                            Component_Absorption( iphase, jphase, cv_nodi ) = &
                                                Component_Absorption( iphase, jphase, cv_nodi ) * max( 0.01, &
                                                20. * ( 1. - saturation_field%val( 1,1, cv_nodi ) ) )
                                        end do
                                    end do
                                end if
                            end do
                        end if Conditional_SmoothAbsorption
                        !!$ NonLinear iteration for the components advection:
                        Loop_NonLinearIteration_Components: do its2 = 1, NonLinearIteration_Components
                            Mdisopt%comp_use_theta_flux = .false. ; Mdisopt%comp_get_theta_flux = .true.
                            call INTENERGE_ASSEM_SOLVE( state, multicomponent_state(icomp), &
                                Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,&
                                tracer_field,velocity_field,density_field, dt, &
                                SUF_SIG_DIAGTEN_BC, Porosity_field%val, &
                                opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
                                igot_t2, igot_theta_flux, &
                                Mdisopt%comp_get_theta_flux, Mdisopt%comp_use_theta_flux, &!NEED TO REMOVE EVERYTHING OF THE NEW MDISOPT BUT THE ONES IN THIS LINE...
                                theta_gdiff, Mean_Pore_CV, &
                                thermal = .false.,& ! the false means that we don't add an extra source term
                                theta_flux=theta_flux, one_m_theta_flux=one_m_theta_flux, theta_flux_j=theta_flux_j, one_m_theta_flux_j=one_m_theta_flux_j,&
                                icomp=icomp, saturation=saturation_field, IDs_ndgln=IDs_ndgln )
                            tracer_field%val = min (max( tracer_field%val, 0.0), 1.0)
                        end do Loop_NonLinearIteration_Components
                        sum_theta_flux = sum_theta_flux + theta_flux
                        sum_one_m_theta_flux = sum_one_m_theta_flux + one_m_theta_flux
                        sum_theta_flux_j = sum_theta_flux_j + theta_flux_j
                        sum_one_m_theta_flux_j = sum_one_m_theta_flux_j + one_m_theta_flux_j
                        ! We have divided through by density
                        tracer_source%val(1,:,:) = tracer_source%val(1,:,:) + THETA_GDIFF
                    end do Loop_Components
                    if ( have_option( '/material_phase[' // int2str( Mdims%nstate - Mdims%ncomp ) // &
                         ']/is_multiphase_component/Comp_Sum2One/Enforce_Comp_Sum2One' ) ) then
                       ! Initially clip and then ensure the components sum to unity so we don't get surprising results...
                       MFC_s % val = min ( max ( MFC_s % val, 0.0), 1.0)
                       ALLOCATE( RSUM( Mdims%nphase ) )
                       DO CV_INOD = 1, Mdims%cv_nonods
                          DO IPHASE = 1, Mdims%nphase
                             RSUM( IPHASE ) = SUM (MFC_s % val (:, IPHASE, CV_INOD) )
                          END DO
                          DO IPHASE = 1, Mdims%nphase
                             MFC_s % val (:, IPHASE, CV_INOD) = MFC_s % val (:, IPHASE, CV_INOD) / RSUM( IPHASE )
                          END DO
                       END DO
                       DEALLOCATE( RSUM )
                    end if
                    DO ICOMP = 1, Mdims%ncomp
                        call Calculate_ComponentAbsorptionTerm( state, packed_state,&
                            icomp, ndgln%cv, Mdims, &
                            D_s%val, Porosity_field%val, mass_ele, &
                            Component_Absorption,IDs_ndgln )
                        if ( have_option( '/material_phase[' // int2str( Mdims%nstate - Mdims%ncomp ) // &
                            ']/is_multiphase_component/KComp_Sigmoid' ) .and. Mdims%nphase > 1 ) then
                            do cv_nodi = 1, Mdims%cv_nonods
                                if( saturation_field%val( 1, 1, cv_nodi ) > 0.95 ) then
                                    do iphase = 1, Mdims%nphase
                                        do jphase = min( iphase + 1, Mdims%nphase ), Mdims%nphase
                                            Component_Absorption( iphase, jphase, cv_nodi ) = &
                                                Component_Absorption( iphase, jphase, cv_nodi ) * max( 0.01, &
                                                20. * ( 1. - saturation_field%val (1,1, cv_nodi ) ) )
                                        end do
                                    end do
                                end if
                            end do
                            DO CV_NODI = 1, Mdims%cv_nonods
                               Loop_Phase_SourceTerm1: do iphase = 1, Mdims%nphase
                                  Loop_Phase_SourceTerm2: do jphase = 1, Mdims%nphase
                                     tracer_source%val(1,iphase,cv_nodi)=tracer_source%val(1,iphase,cv_nodi)- &
                                          Component_Absorption( IPHASE, JPHASE, CV_NODI ) * &
                                          MFC_s%val(ICOMP, JPHASE, CV_NODI) / &
                                          DC_s%val( icomp, iphase, cv_nodi  )
                                  end do Loop_Phase_SourceTerm2
                               end do Loop_Phase_SourceTerm1
                            END DO
                         end if
                        ! For compressibility
                        DO IPHASE = 1, Mdims%nphase
                            DO CV_NODI = 1, Mdims%cv_nonods
                                tracer_source%val(1, iphase, cv_nodi)=tracer_source%val(1, iphase, cv_nodi)&
                                    + Mean_Pore_CV(1, CV_NODI ) * MFCOLD_s%val(ICOMP, IPHASE, CV_NODI) &
                                    * ( DCOLD_s%val(ICOMP, IPHASE, CV_NODI) - DC_s%val(ICOMP, IPHASE, CV_NODI) ) &
                                    * old_saturation_field%val(1, IPHASE, Mdims%cv_nonods) &
                                    / ( DC_s%val(ICOMP, IPHASE, CV_NODI) * DT )
                            END DO
                        END DO
                    END DO ! ICOMP
                    if( have_option( '/material_phase[' // int2str( Mdims%nstate - Mdims%ncomp ) // &
                        ']/is_multiphase_component/Comp_Sum2One' ) .and. ( Mdims%ncomp > 1 ) ) then
                        call Cal_Comp_Sum2One_Sou( packed_state, Mdims%cv_nonods, Mdims%nphase, Mdims%ncomp, dt, its, &
                            NonLinearIteration, Mean_Pore_CV )
                    end if
                end if Conditional_Components
                !Check if the results are good so far and act in consequence, only does something if requested by the user
                if (sig_hup .or. sig_int) then
                    ewrite(1,*) "Caught signal, exiting nonlinear loop"
                    exit Loop_NonLinearIteration
                end if
                call Adaptive_NonLinear(packed_state, reference_field, its,&
                    Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,3)
                if (ExitNonLinearLoop) then
                    if (adapt_mesh_in_FPI) then
                        !Three steps
                        !1. Store OldPhaseVolumeFraction
                        do iphase = 1, Mdims%nphase
                            f1 => extract_scalar_field( state(iphase), "OldPhaseVolumeFraction" )
                            f2  => extract_scalar_field( state(iphase), "Saturation_bak" )
                            f2%val = f1%val
                        end do
                        !2.Prognostic field to adapt the mesh to, has to be a convolution of old a new saturations
                        if (have_option( '/material_phase[0]/scalar_field::SatOldNew')) then
                            f2  => extract_scalar_field( state(1), "OldPhaseVolumeFraction" )
                            f1  => extract_scalar_field( state(1), "SatOldNew" )
                            f1%val = f2%val
                            f2  => extract_scalar_field( state(1), "PhaseVolumeFraction" )
                            f1%val = abs(f1%val - f2%val)**0.8
                        end if
                        call adapt_mesh_mp()
                        !3. Copy back to OldPhaseVolumeFraction
                        do iphase = 1, Mdims%nphase
                            f1 => extract_scalar_field( state(iphase), "OldPhaseVolumeFraction" )
                            f2  => extract_scalar_field( state(iphase), "Saturation_bak" )
                            f1%val = f2%val
                        end do
                        !Point porosity again
                        porosity_field=>extract_vector_field(packed_state,"Porosity")
                        !Now we have to converge again within the same time-step
                        ExitNonLinearLoop = .false.; its = 1
                        adapt_mesh_in_FPI = .false.
                    else
                        exit Loop_NonLinearIteration
                    end if
                end if
                after_adapt=.false.
                its = its + 1
                !Flag the matrices as already calculated (only the storable ones
                Mmat%stored = .true.
            end do Loop_NonLinearIteration
            ! If calculating boundary fluxes, add up contributions to \int{totout} at each time step
            if(calculate_flux) then
                ! We will output totout normalised as f1/(f1+f2). 29/01/2016 - corrected the normalisation (was incorrect for n_outlets > 1)
                do ioutlet = 1, size(outlet_id)
                    intflux(:, ioutlet) = intflux(:, ioutlet) + totout(:, ioutlet)*dt
                    totout(:, ioutlet) = totout(:, ioutlet)/sum(totout(:, ioutlet))
                enddo
            endif
            ! If calculating boundary fluxes, dump them to outfluxes.txt
            if(calculate_flux) then
                if(getprocno() == 1) then
                    call dump_outflux(acctim,itime,totout,intflux)
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
            call Calculate_All_Rhos( state, packed_state, Mdims )
            do iphase = 1, Mdims%nphase
                f1 => extract_scalar_field( state(iphase), "Temperature", stat )
                if ( stat==0 ) then
                    f2 => extract_scalar_field( state(iphase),"Dummy", stat )
                    if ( stat==0 ) f2%val = f1%val
                end if
            end do
            !!$ Calculate diagnostic fields
            call calculate_diagnostic_variables( state, exclude_nonrecalculated = .true. )
            call calculate_diagnostic_variables_new( state, exclude_nonrecalculated = .true. )
            if (write_all_stats) call write_diagnostics( state, current_time, dt, itime )  ! Write stat file
            Conditional_TimeDump: if( ( mod( itime, dump_period_in_timesteps ) == 0 ) ) then
                dtime=dtime+1
                if (do_checkpoint_simulation(dtime)) then
                    CV_Pressure=>extract_tensor_field(packed_state,"PackedCVPressure")
                    FE_Pressure=>extract_tensor_field(packed_state,"PackedFEPressure")
                    call set(pressure_field,FE_Pressure)
                    call checkpoint_simulation(state,cp_no=checkpoint_number,&
                        protect_simulation_name=.true.,file_type='.mpml')
                    checkpoint_number=checkpoint_number+1
                    call set(pressure_field,CV_Pressure)
                end if
                call get_option( '/timestepping/current_time', current_time ) ! Find the current time
                if (.not. write_all_stats)call write_diagnostics( state, current_time, dt, itime/dump_period_in_timesteps )  ! Write stat file
                not_to_move_det_yet = .false. ; dump_no = itime/dump_period_in_timesteps ! Sync dump_no with itime
                call write_state( dump_no, state ) ! Now writing into the vtu files
            end if Conditional_TimeDump
            ! CALL INITIAL MESH TO MESH INTERPOLATION ROUTINE (Before adapting the mesh)
            numberfields=option_count('/material_phase/scalar_field/prognostic/CVgalerkin_interpolation') ! Count # instances of CVGalerkin in the input file
            if (numberfields > 0) then ! If there is at least one instance of CVgalerkin then apply the method
                if (have_option('/mesh_adaptivity')) then ! Only need to use interpolation if mesh adaptivity switched on
                    call M2MInterpolation(state, packed_state, Mdims, CV_GIdims, CV_funs, Mspars%small_acv%fin, Mspars%small_acv%col ,Mdisopt%cv_ele_type, 0)
                endif
            endif
            !!!$! ******************
            !!!$! *** Mesh adapt ***
            !!!$! ******************
            if (adapt_mesh_in_FPI .and. have_option( '/material_phase[0]/scalar_field::SatOldNew')) then
                !Point SatOldNew to PhaseVolumeFraction to remove the extra nodes used in the previous adapt within the FPI
                f1  => extract_scalar_field( state(1), "SatOldNew" )
                f2  => extract_scalar_field( state(1), "PhaseVolumeFraction" )
                f1%val = f2%val
            end if
            call adapt_mesh_mp()
            !!$ Simple adaptive time stepping algorithm
            if ( have_option( '/timestepping/adaptive_timestep' ) ) then
                c = -66.6 ; minc = 0. ; maxc = 66.e6 ; ic = 1.1!66.e6
                call get_option( '/timestepping/adaptive_timestep/requested_cfl', rc )
                call get_option( '/timestepping/adaptive_timestep/minimum_timestep', minc, stat )
                call get_option( '/timestepping/adaptive_timestep/maximum_timestep', maxc, stat )
                call get_option( '/timestepping/adaptive_timestep/increase_tolerance', ic, stat )
                !For porous media we need to use the Courant number obtained in cv_assemb
                if (is_porous_media) then
                    c = max ( c, Courant_number )
                else
                    do iphase = 1, Mdims%nphase
                        ! requested cfl
                        rc_field => extract_scalar_field( state( iphase ), 'RequestedCFL', stat )
                        if ( stat == 0 ) rc = min( rc, minval( rc_field % val ) )
                        ! max cfl
                        cfl => extract_scalar_field( state( iphase ), 'CFLNumber' )
                        c = max ( c, maxval( cfl % val ) )
                    end do
                end if
                call get_option( '/timestepping/timestep', dt )
                dt = max( min( min( dt * rc / c, ic * dt ), maxc ), minc )
                dt = max(min(dt, finish_time - current_time), 1d-15)!Make sure we finish at required time and we don't get dt = 0
                call allmin(dt)
                call set_option( '/timestepping/timestep', dt )
            end if
            if ( do_reallocate_fields ) then
                after_adapt=.true.
            else
                after_adapt=.false.
            end if
            if ( after_adapt .and. have_option( '/mesh_adaptivity/hr_adaptivity/nonlinear_iterations_after_adapt' ) ) then
                call get_option( '/mesh_adaptivity/hr_adaptivity/nonlinear_iterations_after_adapt', NonLinearIteration )
            else
                call get_option( '/timestepping/nonlinear_iterations', NonLinearIteration, default = 3 )
            end if
            call set_boundary_conditions_values(state, shift_time=.true.)
            if (sig_hup .or. sig_int) then
                ewrite(1,*) "Caught signal, exiting"
                exit Loop_Time
            end if
            if(itime==1) then
                is_first_time_step = .false.
            end if
        end do Loop_Time
        if (has_references(metric_tensor)) call deallocate(metric_tensor)
        !!$ Now deallocating arrays:
        deallocate( &
            !!$ Defining element-pair type and discretisation options and coefficients
            opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
            !!$ For output:
            Mean_Pore_CV, &
            !!$ Variables used in the diffusion-like term: capilarity and surface tension:
            plike_grad_sou_grad, plike_grad_sou_coef, &
            !!$ Working arrays
            theta_gdiff, ScalarField_Source, ScalarField_Source_Store, ScalarField_Source_Component, &
            mass_ele,&
            Material_Absorption, &
            ScalarField_Absorption, Component_Absorption, &
            theta_flux, one_m_theta_flux, theta_flux_j, one_m_theta_flux_j, &
            sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j )
        ! Dump at end, unless explicitly disabled
        if(.not. have_option("/io/disable_dump_at_end")) then
            call write_state(dump_no, state)
        end if
        call tag_references()
        call deallocate(packed_state)
        call deallocate(multiphase_state)
        deallocate(multiphase_state)
        call deallocate(multicomponent_state)
        deallocate(multicomponent_state)
        call deallocate(storage_state)
        call deallocate_multi_shape_funs(CV_funs)
        call deallocate_multi_shape_funs(FE_funs)
        call deallocate_multi_ndgln(ndgln)
        call destroy_multi_matrices(Mmat)
        !***************************************
        ! INTERPOLATION MEMORY CLEANUP
        if (numberfields > 0) then
            call MemoryCleanupInterpolation1()     ! Clean up state_old allocations here.
                                                   ! State_new cleanup happens straight after calling the interpolation with flag = 1
                                                   ! inside the adaptivity loop. Probably best to split M2Minterpolation into
                                                   ! separate flag == 0, flag == 1 subroutines and deallocate inside them
                                                   ! (future work).
        endif
        !***************************************
        if (calculate_flux) deallocate(outlet_id, totout, intflux)
        return
    contains
        !!!!!sprint_to_do!!!rename and/or move elsewhere
        subroutine temp_mem_hacks()
            !!! routine puts various CSR sparsities into packed_state
            use sparse_tools
            type(csr_sparsity), pointer :: sparsity
            type(tensor_field), pointer :: tfield
            integer ic, stat
            type(halo_type), pointer :: halo
            type(mesh_type), pointer :: ph_mesh
            ph_mesh => extract_mesh( state( 1 ), "ph", stat )
            if ( stat == 0 ) then
                allocate( sparsity )
                sparsity = wrap( Mspars%ph%fin, colm = Mspars%ph%col, name = "phsparsity" )
                call insert( packed_state, sparsity, "phsparsity" )
                call deallocate( sparsity )
                deallocate ( sparsity )
            end if
            allocate(sparsity)
            sparsity=wrap(Mspars%ACV%fin,Mspars%ACV%mid,colm=Mspars%ACV%col,name='PackedAdvectionSparsity')
            call insert(packed_state,sparsity,'PackedAdvectionSparsity')
            call deallocate(sparsity)
            sparsity=wrap(Mspars%C%fin,colm=Mspars%C%col,name='CMatrixSparsity')
            call insert(packed_state,sparsity,'CMatrixSparsity')
            call deallocate(sparsity)
            sparsity=wrap(Mspars%CT%fin,colm=Mspars%CT%col,name='CTMatrixSparsity')
            call insert(packed_state,sparsity,'CTMatrixSparsity')
            call deallocate(sparsity)
            sparsity=wrap(Mspars%M%fin,Mspars%M%mid,colm=Mspars%M%col,name='CVFEMSparsity')
            call insert(packed_state,sparsity,'CVFEMSparsity')
            call deallocate(sparsity)
            tfield=>extract_tensor_field(packed_state,"PackedFEPressure")
            if (associated( tfield%mesh%halos)) then
                sparsity=wrap(Mspars%CMC%fin,colm=Mspars%CMC%col,name='CMCSparsity',&
                    row_halo=tfield%mesh%halos(2),&
                    column_halo=tfield%mesh%halos(2))
            else
                sparsity=wrap(Mspars%CMC%fin,colm=Mspars%CMC%col,name='CMCSparsity')
            end if
            call insert(packed_state,sparsity,'CMCSparsity')
            call deallocate(sparsity)
            if (associated( tfield%mesh%halos)) then
                sparsity=wrap(Mspars%small_acv%fin,Mspars%small_acv%mid,colm=Mspars%small_acv%col,name="ACVSparsity",&
                    row_halo=tfield%mesh%halos(2),&
                    column_halo=tfield%mesh%halos(2))
            else
                sparsity=wrap(Mspars%small_acv%fin,Mspars%small_acv%mid,colm=Mspars%small_acv%col,name="ACVSparsity")
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
                sparsity=wrap(Mspars%DGM_PHA%fin,colm=Mspars%DGM_PHA%col,&
                    name='MomentumSparsity',row_halo=halo,column_halo=halo)
            else
                sparsity=wrap(Mspars%DGM_PHA%fin,colm=Mspars%DGM_PHA%col,name="MomentumSparsity")
            end if
            call insert(packed_state,sparsity,"MomentumSparsity")
            call deallocate(sparsity)
            tfield=>extract_tensor_field(packed_state,"PackedFEPressure")
            sparsity=make_sparsity(tfield%mesh,tfield%mesh,&
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
        end subroutine temp_mem_hacks
        subroutine linearise_components()
            integer :: ist,ip,ele
            type ( scalar_field ), pointer :: nfield
            integer, dimension(:), pointer :: nodes
            real, allocatable, dimension(:) :: comp
            if (Mdims%ncomp>1) then
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
                    do ip=1,Mdims%nphase
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
        subroutine adapt_mesh_mp()
            do_reallocate_fields = .false.
            Conditional_Adaptivity_ReallocatingFields: if( have_option( '/mesh_adaptivity/hr_adaptivity') ) then
                if( have_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps') ) then
                    call get_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps', &
                        adapt_time_steps, default=5 )
                else if (have_option( '/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI')) then
                    do_reallocate_fields = .true.
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
                !The stored variables must be recalculated
                Conditional_Adaptivity: if( have_option( '/mesh_adaptivity/hr_adaptivity ') .or. have_option( '/mesh_adaptivity/hr_adaptivity_prescribed_metric')) then
                    Conditional_Adapt_by_TimeStep: if( mod( itime, adapt_time_steps ) == 0 .or. have_option( '/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI')) then
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
                        call run_diagnostics( state )
                        call adapt_state( state, metric_tensor, suppress_reference_warnings = .true.)
                        call update_state_post_adapt( state, metric_tensor, dt, sub_state, nonlinear_iterations, &
                            nonlinear_iterations_adapt )
                        if( have_option( '/io/stat/output_after_adapts' ) ) call write_diagnostics( state, current_time, dt, &
                            itime, not_to_move_det_yet = .true. )
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
                deallocate(multiphase_state)
                deallocate(multicomponent_state)
                call pack_multistate(Mdims%npres,state,packed_state,&
                    multiphase_state,multicomponent_state)
                call set_boundary_conditions_values(state, shift_time=.true.)
                !!$ Deallocating array variables:
                deallocate( &
                    !!$ Defining element-pair type and discretisation options and coefficients
                    opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, &
                    !!$ For output:
                    Mean_Pore_CV, &
                    !!$ Variables used in the diffusion-like term: capilarity and surface tension:
                    plike_grad_sou_grad, plike_grad_sou_coef, &
                    !!$ Working arrays
                    suf_sig_diagten_bc, &
                    theta_gdiff, ScalarField_Source, ScalarField_Source_Store, ScalarField_Source_Component, &
                    mass_ele, &
                    Material_Absorption, &
                    ScalarField_Absorption, Component_Absorption, &
                    theta_flux, one_m_theta_flux, theta_flux_j, one_m_theta_flux_j, sum_theta_flux, &
                    sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j )
                !Deallocate sparsities
                call deallocate_multi_sparsities(Mspars)
                !Deallocate ndgln
                call deallocate_multi_ndgln(ndgln)
                !Destroy what remains of the matrices
                call destroy_multi_matrices(Mmat)
                !!$ Compute primary scalars used in most of the code
                call Get_Primary_Scalars_new( state, Mdims )
                !!$ Calculating Global Node Numbers
                call allocate_multi_ndgln(ndgln, Mdims)
                call Compute_Node_Global_Numbers(state, ndgln)
                !!$
                !!$ Computing Sparsity Patterns Matrices
                !!$
                !!$ Defining lengths and allocating space for the matrices
                call Defining_MaxLengths_for_Sparsity_Matrices( Mdims%ndim, Mdims%nphase, Mdims%totele, Mdims%u_nloc, Mdims%cv_nloc, Mdims%ph_nloc, Mdims%cv_nonods, &
                    mx_nface_p1, mxnele, mx_nct, mx_nc, mx_ncolcmc, mx_ncoldgm_pha, mx_ncolmcy, &
                    mx_ncolacv, mx_ncolm, mx_ncolph )
                nlenmcy = Mdims%u_nonods * Mdims%nphase * Mdims%ndim + Mdims%cv_nonods
                !!$ Defining element-pair type
                call Get_Ele_Type( Mdims%x_nloc, Mdisopt%cv_ele_type, Mdisopt%p_ele_type, Mdisopt%u_ele_type, &
                    Mdisopt%mat_ele_type, Mdisopt%u_sele_type, Mdisopt%cv_sele_type )
                !Allocate and calculate the sparsity patterns
                call Get_Sparsity_Patterns( state, Mdims, Mspars, ndgln, Mdisopt, mx_ncolacv, nlenmcy, mx_ncolmcy, &
                    mx_ncoldgm_pha, mx_nct,mx_nc, mx_ncolcmc, mx_ncolm, mx_ncolph, mx_nface_p1 )
                if (is_porous_media) then
                    !Re-calculate IDs_ndgln after adapting the mesh
                    call get_RockFluidProp(state, packed_state)
                    !Convert material properties to be stored using region ids, only if porous media
                    call get_regionIDs2nodes(state, packed_state, ndgln%cv, IDs_ndgln, IDs2CV_ndgln, fake_IDs_ndgln = .not. is_porous_media)
                end if
                call temp_mem_hacks()
                ! SECOND INTERPOLATION CALL - After adapting the mesh ******************************
                if (numberfields > 0) then
                    if(have_option('/mesh_adaptivity')) then ! This clause may be redundant and could be removed - think this code in only executed IF adaptivity is on
                        call M2MInterpolation(state, packed_state, Mdims, CV_GIdims, CV_funs, &
                                Mspars%small_acv%fin, Mspars%small_acv%col, Mdisopt%cv_ele_type , 1, IDs2CV_ndgln = IDs2CV_ndgln)
                        call MemoryCleanupInterpolation2()
                    endif
                endif
                ! *************************************************************
                !!$ Allocating space for various arrays:
                allocate( &
                    !!$
                    suf_sig_diagten_bc( Mdims%stotel * Mdims%cv_snloc * Mdims%nphase, Mdims%ndim ), &
                    Mean_Pore_CV( Mdims%npres, Mdims%cv_nonods ), &
                    mass_ele( Mdims%totele ), &
                    !!$
                    Material_Absorption( Mdims%mat_nonods, Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase ), &
                    ScalarField_Absorption( Mdims%nphase, Mdims%nphase, Mdims%cv_nonods ), Component_Absorption( Mdims%nphase, Mdims%nphase, Mdims%cv_nonods ), &
                    !!$ Variables used in the diffusion-like term: capilarity and surface tension:
                    plike_grad_sou_grad( Mdims%ncomp, Mdims%nphase, Mdims%cv_nonods ), &
                    plike_grad_sou_coef( Mdims%ncomp, Mdims%nphase, Mdims%cv_nonods ) )
                !!$
                Component_Absorption=0.
                !!$
                ScalarField_Absorption=0.
                !!$
                Material_Absorption=0.
                !!$
                plike_grad_sou_grad=0. ; plike_grad_sou_coef=0.
                !!$
                suf_sig_diagten_bc=0.
                !!$
                !!$
                !!$ Initialising Absorption terms that do not appear in the schema
                !!$
                ScalarField_Absorption = 0. ; Component_Absorption = 0.
                !!$ Computing shape function scalars
                igot_t2 = 0 ; igot_theta_flux = 0
                if( Mdims%ncomp /= 0 )then
                    igot_t2 = 1 ; igot_theta_flux = 1
                end if
                scvngi_theta = CV_GIdims%scvngi
                ncv_faces = CV_count_faces( Mdims, Mdisopt%cv_ele_type, CV_GIdims)
                allocate( theta_flux( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    one_m_theta_flux( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    theta_flux_j( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    one_m_theta_flux_j( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    sum_theta_flux( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    sum_one_m_theta_flux( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    sum_theta_flux_j( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    sum_one_m_theta_flux_j( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    theta_gdiff( Mdims%nphase, Mdims%cv_nonods ), ScalarField_Source( Mdims%nphase, Mdims%cv_nonods ), &
                    ScalarField_Source_Store( Mdims%nphase, Mdims%cv_nonods ), &
                    ScalarField_Source_Component( Mdims%nphase, Mdims%cv_nonods ) )
                sum_theta_flux = 1. ; sum_one_m_theta_flux = 0.
                sum_theta_flux_j = 1. ; sum_one_m_theta_flux_j = 0.
                ScalarField_Source=0. ; ScalarField_Source_Store=0. ; ScalarField_Source_Component=0.
                allocate(opt_vel_upwind_coefs_new(Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%mat_nonods)); opt_vel_upwind_coefs_new =0.
                allocate(opt_vel_upwind_grad_new(Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%mat_nonods)); opt_vel_upwind_grad_new =0.
                if( have_option( '/material_phase[' // int2str( Mdims%nstate - Mdims%ncomp ) // &
                    ']/is_multiphase_component/Comp_Sum2One/Enforce_Comp_Sum2One' ) ) then
                    ! Initially clip and then ensure the components sum to unity so we don't get surprising results...
                    MFC_s  => extract_tensor_field( packed_state, "PackedComponentMassFraction" )
                    MFC_s % val = min ( max ( MFC_s % val, 0.0), 1.0)
                    ALLOCATE( RSUM( Mdims%nphase ) )
                    DO CV_INOD = 1, Mdims%cv_nonods
                        DO IPHASE = 1, Mdims%nphase
                            RSUM( IPHASE ) = SUM (MFC_s % val (:, IPHASE, CV_INOD) )
                        END DO
                        DO IPHASE = 1, Mdims%nphase
                            MFC_s % val (:, IPHASE, CV_INOD) = MFC_s % val (:, IPHASE, CV_INOD) / RSUM( IPHASE )
                        END DO
                    END DO
                    DEALLOCATE( RSUM )
                end if
                call Calculate_All_Rhos( state, packed_state, Mdims )
            end if Conditional_ReallocatingFields
        end subroutine adapt_mesh_mp
    end subroutine MultiFluids_SolveTimeLoop
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
        tfield=>extract_tensor_field(packed_state,"PackedOldFEPressure")
        ntfield=>extract_tensor_field(packed_state,"PackedFEPressure")
        tfield%val=ntfield%val
        tfield=>extract_tensor_field(packed_state,"PackedOldCVPressure")
        ntfield=>extract_tensor_field(packed_state,"PackedCVPressure")
        tfield%val=ntfield%val
    end subroutine copy_packed_new_to_old
        !!!!!sprint_to_do!!!rename and/or move elsewhere
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
    !!!!sprint_to_do!!!!take this somewhere else
    subroutine dump_outflux(current_time, itime, outflux, intflux)
        ! Subroutine that dumps the total flux at a given timestep across all specified boudaries to a file  called 'outfluxes.txt'. In addition, the time integrated flux
        ! up to the current timestep is also outputted to this file. Integration boundaries are specified in diamond via surface_ids.
        ! (In diamond this option can be found under "/io/dump_boundaryflux/surface_ids" and the user should specify an integer array containing the IDs of every boundary they
        !wish to integrate over).
        real,intent(in) :: current_time
        integer, intent(in) :: itime
        real, dimension(:,:), intent(inout) :: outflux, intflux
        integer :: ioutlet
        integer :: counter
        type(stat_type), target :: default_stat
        character (len=1000000) :: whole_line
        character (len=1000000) :: numbers
        integer :: iphase
        ! Strictly speaking don't need character arrays for fluxstring and intfluxstring, could just overwrite each time (may change later)
        character (len = 1000000), dimension(size(outflux,1)) :: fluxstring
        character (len = 1000000), dimension(size(outflux,1)) :: intfluxstring
        default_stat%conv_unit=free_unit()
        if (itime == 1) then
            !The first time, remove file if already exists
            open(unit=default_stat%conv_unit, file="outfluxes.csv", status="replace", action="write")
        else
            open(unit=default_stat%conv_unit, file="outfluxes.csv", action="write", position="append")
        end if
        ! Write column headings to file
        counter = 0
        if(itime.eq.1) then
            write(whole_line,*) "Current Time"
            do ioutlet =1, size(outflux,2)
                write(numbers,*) "Surface_id=", outlet_id(ioutlet)
                if(counter.eq.0) then
                    whole_line = trim(numbers) //","// trim(whole_line)
                else
                    whole_line = trim(whole_line) //","// trim(numbers)
                endif
                !write(whole_line,*)trim(numbers)  //","//  "Current Time"
                do iphase = 1, size(outflux,1)
                    write(fluxstring(iphase),*) "Phase", iphase, "boundary flux"
                    whole_line = trim(whole_line) //","// trim(fluxstring(iphase))
                enddo
                do iphase = 1, size(outflux,1)
                    write(intfluxstring(iphase),*) "Phase", iphase,  "time integrated flux (volume/time)"
                    whole_line = trim(whole_line) //","// trim(intfluxstring(iphase))
                enddo
                counter = counter + 1
            end do
             ! Write out the line
            write(default_stat%conv_unit,*), trim(whole_line)
        else
            ! Write the actual numbers to the file now
            counter = 0
            write(numbers,*) current_time
            whole_line =  ","// trim(numbers)
            do ioutlet =1, size(outflux,2)
                if(counter > 0) then
                    whole_line = trim(whole_line) //","
                endif
                !write(whole_line,*) current_time
                do iphase = 1, size(outflux,1)
                    write(fluxstring(iphase),*) outflux(iphase,ioutlet)
                    whole_line = trim(whole_line) //","// trim(fluxstring(iphase))
                enddo
                do iphase = 1, size(outflux,1)
                    write(intfluxstring(iphase),*) intflux(iphase,ioutlet)
                    whole_line = trim(whole_line) //","// trim(intfluxstring(iphase))
                enddo
                counter = counter + 1
            end do
            ! Write out the line
            write(default_stat%conv_unit,*), trim(whole_line)
        endif
        close (default_stat%conv_unit)
    end subroutine dump_outflux
end module multiphase_time_loop
