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
!----------------------------------------------------------------------------------------
!> @brief Time-loop module of IC-FERST. This module contains the time-loop and the non-linear loop.
!> The time-loop consists many steps: 1) Data initialisation including shape functions, memory allocation, sparsity, porous media properties, etc.
!> 2) Initialisation of the actual time loop.
!> 3) Non-linear loop. IC-FERST uses a modified Anderson-acceleration non-linear solver, which is based in a Picard iterative non-linear solver.
!> In this way, the different equations are solved independently and coupled through the Anderson non-linear solver.
!> First the momentum and continuity equations are assembled and solved for, next the different transport equations are solved, for example: saturation,
!> temperature, concentration, etc. ActiveTracers and Species are solved within the non-linear solver while PassiveTracers are solved outside the non-linear solver.
!> 4) Once the non-linear solver has converged, we proceed to jump to the next time-level but first we check if we need to adapt the mesh and/or generate a vtu file.
!> etc.
!>
!----------------------------------------------------------------------------------------
module multiphase_time_loop

#ifdef HAVE_PETSC_MODULES
  use petsc
#endif

#ifdef USING_PHREEQC
  use multi_phreeqc
#endif

#ifdef USING_XGBOOST
    use multi_machine_learning
#endif

    use field_options
    use write_state_module
    use diagnostic_variables
    use diagnostic_fields_wrapper
    use diagnostic_fields_new, only : &
        calculate_diagnostic_variables_new => calculate_diagnostic_variables
    use global_parameters
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
    !use global_parameters ! repeated?
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
    use Compositional_Terms
    use Copy_Outof_State
    use cv_advection, only : cv_count_faces
    use multiphase_1D_engine
    use boundary_conditions_from_options
    use multi_data_types
    use vtk_interfaces
    use multi_interpolation
    use multi_pipes
    use multi_SP
    use multi_tools
    use momentum_diagnostic_fields, only: calculate_densities
    use embed_python

#ifdef HAVE_ZOLTAN
  use zoltan
#endif


    !use matrix_operations
    !use shape_functions
    implicit none

#include "petsc_legacy.h"

    private
    !public :: MultiFluids_SolveTimeLoop, rheology, dump_outflux
    public :: MultiFluids_SolveTimeLoop, dump_outflux
    !type(rheology_type), dimension(:), allocatable :: rheology



contains

    !> This is the main subroutine from which everything is called. It performs the time-loop and
    !> therefore calls all the necessary blocks to solve for the system of equations,
    !> adapt the mesh, etc.
    !> See Salinas et al. 2016 doi: 10.1002/fld.4357, for a description of the time-loop and non-linear solver
    !> @param state. Linked list type of memory where all the fields following Fluidity standard are stored
    !> @param dump_no. Specify the initial number of the vtu files, has to be zero
    !> @param nonlinear_iterations. Non-linear value. This is required just in case we adapt before the first time-step
    !> @param dt. Time-step size
    !> @retval
    !> @see
    !> @ref
    subroutine MultiFluids_SolveTimeLoop( state, &
        dt, nonlinear_iterations, dump_no )
        implicit none


        type( state_type ), dimension( : ), intent( inout ) :: state!> Linked list type of memory where all the fields following Fluidity standard are stored
        integer, intent( inout ) :: dump_no !> Specify the initial number of the vtu files, has to be zero
        integer, intent( inout ) :: nonlinear_iterations !> This value is required just in case we adapt before the first time-step
        real, intent( inout ) :: dt !> Time-step size
        !Local variables
        !!$ additional state variables for multiphase & multicomponent
        type(state_type) :: packed_state!> Linked list type of memory where all the fields following IC-FERST standard are stored, memory is shared with state
        type(state_type), dimension(:), pointer :: multiphase_state!> Not sure why we need this extra state memory...
        type(state_type), dimension(:), pointer :: multicomponent_state!> For multicomponents we have this extra state file where the compositional memory is stored
        !!Define shape functions
        type (multi_shape_funs) :: CV_funs!> Structure containing everything related with the shape functions of the CV mesh
        type (multi_shape_funs) :: FE_funs!> Structure containing everything related with the shape functions of the FE mesh
        !!$ Primary scalars
        type(multi_dimensions) :: Mdims !> Structure that contains all the dimensions that are used throughout the code, dimensions, number of nodes, ...
        type(multi_gi_dimensions) :: CV_GIdims!> Structure containing everything related with the gauss integration for the CV mesh
        type(multi_gi_dimensions) :: FE_GIdims!> Structure containing everything related with the gauss integration for the FE mesh
        !!$ Node global numbers
        type(multi_ndgln) :: ndgln !> Structure containing information about the node global numbers, conversion from local to global
        !!$ Sparsity patterns
        type (multi_sparsities) :: Mspars!> Structure containing the different sparsities for the different data required, CMC matrix, Ct, C, etc...
        !!$ Defining element-pair type and discretisation options and coefficients
        type (multi_discretization_opts) :: Mdisopt !> Defining element-pair type and discretisation options and coefficients
        !!$ Defining the necessary matrices and corresponding RHS
        type (multi_matrices) :: Mmat !> Structure containing the matrices required to solve for the system
        !!$ Defining variables to calculate the sigmas at the interface for porous media
        type (porous_adv_coefs) :: upwnd !> Structure containing the information to compute the fluxes for porous media, sigma in the papers but without the permeability
        !!$ Variable storing all the absorptions we may need
        type(multi_absorption) :: multi_absorp !> Structure storing all the absorption terms (term multipliying the velocity in the momentum equation)
        integer :: mx_nface_p1, mx_ncolacv, mxnele, mx_ncoldgm_pha, &
            mx_nct, mx_nc, mx_ncolcmc, mx_ncolm, mx_ncolph
        !!$ Defining time- and nonlinear interations-loops variables
        integer :: itime, dump_period_in_timesteps, final_timestep, &
            NonLinearIteration, NonLinearIteration_Components, itimeflag
        real :: acctim, finish_time, dump_period
        !!$ Defining problem that will be solved
        logical :: have_temperature_field, have_concentration_field, have_component_field, have_extra_DiffusionLikeTerm, &
            solve_force_balance, solve_PhaseVolumeFraction
        !!$ Shape function related fields:
        integer :: scvngi_theta, igot_t2, igot_theta_flux
        !!$ Adaptivity related fields and options:
        type( tensor_field ) :: metric_tensor

        PetscErrorCode :: ierrr
        PetscLogStage,dimension(0:9) :: stages


        type( state_type ), dimension( : ), pointer :: sub_state => null()
        integer :: nonlinear_iterations_adapt
        logical :: do_reallocate_fields = .false., not_to_move_det_yet = .false.
        !!$ Working arrays:
        real, dimension(:), pointer :: mass_ele
        real, dimension( :, : ), pointer :: THETA_GDIFF
        !!$
        real, dimension( :, : ), pointer :: &
            ScalarField_Source_Store
        real, dimension( :, :, : ), allocatable :: &
            Velocity_Absorption, Temperature_Absorption
        real, dimension( :, : ), allocatable ::sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j
        integer :: stat, python_stat, istate, iphase, jphase, icomp, its, its2, cv_nodi, adapt_time_steps, cv_inod, vtu_its_counter, SFPI_taken, pres_its_taken
        real, dimension( : ), allocatable :: rsum
        real, dimension(:, :), allocatable :: SUF_SIG_DIAGTEN_BC
        type( scalar_field ), pointer :: cfl, rc_field
        real :: c, rc, minc, maxc, ic, FPI_eq_taken
        !Variables for adaptive time stepping based on non-linear iterations
        logical :: nonLinearAdaptTs, Repeat_time_step, ExitNonLinearLoop
        real, dimension(:,:,:), allocatable  :: reference_field
        !! face value storage
        integer :: ncv_faces
        !Courant number for porous media
        real, dimension(2) :: Courant_number = -1!Stored like this[Courant_number, Shock-front Courant number]
        !Variables for adapting the mesh within the FPI solver
        logical :: adapt_mesh_in_FPI
        real :: Accum_Courant = 0., Courant_tol
        !type( scalar_field ) :: Saturation_bak, ConvSats
        integer :: checkpoint_number
        !Array to map nodes to region ids
        !Variable to store where we store things. Do not oversize this array, the size has to be the last index in use
        !Working pointers
        type(tensor_field), pointer :: tracer_field, tracer_field2, velocity_field, density_field, saturation_field, old_saturation_field   !, tracer_source
        type(tensor_field), pointer :: pressure_field, cv_pressure, fe_pressure, PhaseVolumeFractionSource, PhaseVolumeFractionComponentSource
        type(tensor_field), pointer :: perm_field
        type(vector_field), pointer :: positions, porosity_field, MeanPoreCV, PythonPhaseVolumeFractionSource
        type(scalar_field), pointer :: DensitySource, T
        !Variables that are used to define the pipe pos
        type(pipe_coords), dimension(:), allocatable:: eles_with_pipe
        type (multi_pipe_package) :: pipes_aux
        !type(scalar_field), pointer :: bathymetry
        logical :: write_all_stats=.false.
        ! Variables used for calculating boundary outfluxes. Logical "calculate_flux" determines if this calculation is done. Intflux is the time integrated outflux
        ! Ioutlet counts the number of boundaries over which to calculate the outflux
        integer :: ioutlet
        type (multi_outfluxes) :: outfluxes
        ! Variables used in the CVGalerkin interpolation calculation
        integer, save :: numberfields_CVGalerkin_interp = -1
        real :: t_adapt_threshold
        ! Calculate_mass_delta to store the change in mass calculated over the whole domain
        real, allocatable, dimension(:,:) :: calculate_mass_delta

        !!-Variable to keep track of dt reduction for meeting dump_period requirements
        real, save :: stored_dt = -1
        real :: old_acctim, nonlinear_dt

        !! Variables to initialise porous media models
        logical :: exit_initialise_porous_media = .false.
        !Generic variables
        integer :: i, k
        real :: auxR
        ! Andreas. Declare the parameters required for skipping pressure solve
        Integer:: rcp                 !Requested-cfl-for-Pressure. It is a multiple of CFLNumber
        Logical:: EnterSolve =.true., after_adapt_itime =.false.  !Flag to either enter or not the pressure solve

        integer :: SFPI_its = 0
        integer :: max_sat_its

        real :: bulk_power
        character(len = OPTION_PATH_LEN) :: func
        !Variables for passive tracers
        logical :: have_Active_Tracers = .true.
        logical :: have_Passive_Tracers = .true.
        integer :: fields
        character( len = option_path_len ) :: option_name
        integer :: phreeqc_id
        double precision, ALLOCATABLE, dimension(:,:) :: concetration_phreeqc
#ifdef HAVE_ZOLTAN
      real(zoltan_float) :: ver
      integer(zoltan_int) :: ierr
      ierr = Zoltan_Initialize(ver)
      assert(ierr == ZOLTAN_OK)
#endif

#ifdef USING_XGBOOST
#else
    if (have_option("/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/ML_model_path")) then
        ewrite(0,*) "WARNING: ML path specified for a model but ICFERST hasn't been compiled with XGboost. Classical method used instead."
    end if
#endif

        !We may not want to compute always the stats
        write_all_stats = .not.have_option("/io/Sync_stat_file_with_vtu")

        ! Check wether we are using the CV_Galerkin method
        numberfields_CVGalerkin_interp=option_count('/material_phase/scalar_field/prognostic/CVgalerkin_interpolation') ! Count # instances of CVGalerkin in the input file

        if (numberfields_CVGalerkin_interp > 0 .and. isParallel()) then
          ewrite(1,*) "WARNING: CVGalerkin projection not well tested for parallel."
        end if

        ! A SWITCH TO DELAY MESH ADAPTIVITY UNTIL SPECIFIED UNSER INPUT TIME t_adapt_threshold
        call get_option("/mesh_adaptivity/hr_adaptivity/t_adapt_delay", t_adapt_threshold, default = 0.0 )

        !Read info for adaptive timestep based on non_linear_iterations
        if(have_option("/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep")) then
            if(have_option("/solver_options/Non_Linear_Solver/nonlinear_iterations_at_adapt")) then
                call get_option('/solver_options/Non_Linear_Solver/nonlinear_iterations_at_adapt',nonlinear_iterations_adapt)
                nonlinear_iterations = nonlinear_iterations_adapt
            end if
            call adapt_state_first_timestep(state)
            call allocate_and_insert_auxilliary_fields(state)
            ! Ensure that checkpoints do not adapt at first timestep.
            call delete_option(&
                "/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep")
        end if

        !!$ Compute primary scalars used in most of the code
        call Get_Primary_Scalars_new( state, Mdims )
        !Check if the user wants to store the outfluxes
        call initialize_multi_outfluxes(outfluxes)
        if(use_sub_state()) then
            call populate_sub_state(state,sub_state)
        end if

        call pack_multistate( Mdims%npres, state, packed_state, multiphase_state, &
            multicomponent_state )
        call prepare_absorptions(state, Mdims, multi_absorp)
        call set_boundary_conditions_values(state, shift_time=.false.)

        !  Access boundary conditions via a call like
        !  call get_entire_boundary_condition(extract_tensor_field(packed_state,"Packed"//name),["weakdirichlet"],tfield,bc_type_list)
        !  where tfield is type(tensor_field) and bc_type_list is integer, dimension(tfield%dim(1),tfield%dim(2),nonods)
        !  Then values are in tfield%val(1/Mdims%ndim/Mdims%ncomp,Mdims%nphase,nonods)
        !  Type ids are in bc_type_list(1/Mdims%ndim/Mdims%ncomp,Mdims%nphase,Mdims%stotel)
        !
        !  A deallocate tfield when finished!!

        Repeat_time_step = .false.
        nonLinearAdaptTs = have_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/adaptive_timestep_nonlinear')

        !!$ Calculating Global Node Numbers
        call allocate_multi_ndgln(ndgln, Mdims)
        call Compute_Node_Global_Numbers(state, ndgln)

        call Defining_MaxLengths_for_Sparsity_Matrices( Mdims%ndim, Mdims%nphase, Mdims%totele, Mdims%u_nloc, Mdims%cv_nloc, Mdims%ph_nloc, Mdims%cv_nonods, &
            mx_nface_p1, mxnele, mx_nct, mx_nc, mx_ncolcmc, mx_ncoldgm_pha, &
            mx_ncolacv, mx_ncolm, mx_ncolph )
        !!$ Defining element-pair type
        call Get_Ele_Type_new( Mdims, Mdisopt )
        !Allocate and calculate the sparsity patterns matrices
        call Get_Sparsity_Patterns( state, Mdims, Mspars, ndgln, Mdisopt, mx_ncolacv, &
                mx_ncoldgm_pha, mx_nct,mx_nc, mx_ncolcmc, mx_ncolm, mx_ncolph, mx_nface_p1 )
        call put_CSR_spars_into_packed_state()
        !!$ Allocating space for various arrays:
        allocate( &
            !!$
            suf_sig_diagten_bc( Mdims%stotel * Mdims%cv_snloc * Mdims%nphase, Mdims%ndim ), &
            mass_ele( Mdims%totele ) )
        !!$
        suf_sig_diagten_bc=0.
        mass_ele=0.
        !!$

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
        call cv_fem_shape_funs( CV_funs, Mdims, CV_GIdims, Mdisopt%cv_ele_type, quad_over_whole_ele = .false. )
        call cv_fem_shape_funs( FE_funs, Mdims, FE_GIdims, Mdisopt%u_ele_type, quad_over_whole_ele = .true. )
        !Obtain the number of faces in the control volume space
        ncv_faces=CV_count_faces( Mdims, Mdisopt%cv_ele_type, CV_GIDIMS = CV_GIdims)
        allocate( sum_theta_flux( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            sum_one_m_theta_flux( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            sum_theta_flux_j( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            sum_one_m_theta_flux_j( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            theta_gdiff( Mdims%nphase, Mdims%cv_nonods ), &
            ScalarField_Source_Store( Mdims%nphase, Mdims%cv_nonods ) )
        sum_theta_flux = 1. ; sum_one_m_theta_flux = 0.
        sum_theta_flux_j = 1. ; sum_one_m_theta_flux_j = 0.
        ScalarField_Source_Store=0.
        !!$ Defining discretisation options
        call Get_Discretisation_Options( state, Mdims, Mdisopt )
        !Check if the pressure matrix is a CV matrix
        Mmat%CV_pressure = .not. have_option( '/geometry/Advance_options/FE_Pressure' )
        if (.not. Mmat%CV_pressure .and. ((Mdims%ndim==2 .and. Mdims%u_nloc == 4) .or. (Mdims%ndim==3 .and. Mdims%u_nloc == 5))) then
            ewrite(0, *) "WARNING: the only tested element pair using bubble shape functions is the P1DG(BL)P1DG(CV)"
        end if
        !Warning message for compressible flows. If one phase has compressibility all of them must be define as compressible
        if (have_option_for_any_phase("phase_properties/Density/compressible", Mdims%ndim) .or. &
          have_option_for_any_phase("phase_properties/Density/python_state", Mdims%ndim)) then
          do i = 1, Mdims%nphase
            if (getprocno() == 1 .and. have_option('/material_phase[' // int2str( i - 1 ) // ']/phase_properties/Density/incompressible')) then
                ewrite(0, *) "WARNING: All the phases must be defined as compressible. You can use the linear option with A ~ 0 for the incompressible phase."
                exit
            end if

          end do
        end if

!Check that we have specified Boundary conditions for density
        if (have_option_for_any_phase("phase_properties/Density/compressible", Mdims%ndim) .or. &
          have_option_for_any_phase("phase_properties/Density/python_state", Mdims%ndim) .or. &
          have_option( '/material_phase[0]/scalar_field::Pressure/prognostic/hydrostatic_boundaries' )) then
          do i = 1, Mdims%nphase
            if (getprocno() == 1 .and. &
            (.not. (have_option('/material_phase[' // int2str( i - 1 ) // ']/phase_properties/Density/boundary_conditions') .or. &
            have_option('/material_phase[' // int2str( i - 1 ) // ']/scalar_field::Density/prognostic')))) then
                ewrite(0, *) "WARNING: It is VERY important to have boundary conditions for density if your model is compressible or using hydrostatic_BCs."
            end if
          end do
        end if
        !Check if we want to use a compacted mass matrix
        if ((Mmat%CV_pressure .or. have_option('/numerical_methods/simple_mass_matrix')) &
                    .and. is_porous_media .and. Mdims%npres == 1) then
        !sprint_to_do for this to work with wells we need to change the sparsity, but that still needs to be done!
            call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/polynomial_degree', i )
            call get_option( '/geometry/mesh::PressureMesh/from_mesh/mesh_shape/polynomial_degree', k )
            Mmat%compact_PIVIT_MAT = (i == (k - 1))!This does not include the P1DG(BL)P1DG(CV) element pair...maybe we should include it as well
        end if
        !!$ Defining problem to be solved:
        solve_PhaseVolumeFraction = have_option("/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic")
        solve_force_balance = have_option("/material_phase[0]/scalar_field::Pressure/prognostic")
        !!$ Setting up variables for the Time- and NonLinear Iterations-Loops:
        call get_option( '/timestepping/current_time', acctim )
        mdims%init_time = acctim
        call get_option( '/timestepping/timestep', dt )
        call get_option( '/timestepping/finish_time', finish_time )
        if ( have_option('/io/dump_period_in_timesteps') ) then
            if ( have_option('/io/dump_period_in_timesteps/python') ) then
                call get_option("/io/dump_period_in_timesteps/python", func)
                call integer_from_python(func, acctim, dump_period_in_timesteps)
                !Ensure it is bounded
                dump_period_in_timesteps = min(max(dump_period_in_timesteps,0), 10000000)
            else
                call get_option('/io/dump_period_in_timesteps/constant', dump_period_in_timesteps, default = 1)
            end if
        elseif ( have_option('/io/dump_period') ) then
            call get_option('/io/dump_period/constant', dump_period, default = 0.01)
        end if
        call get_option( '/solver_options/Non_Linear_Solver', NonLinearIteration, default = 3 )
        !!$
        have_temperature_field = .false. ; have_component_field = .false. ; have_extra_DiffusionLikeTerm = .false.
        do istate = 1, Mdims%nstate
            if( have_option( '/material_phase[' // int2str( istate - 1 ) // ']/scalar_field::Temperature' ) ) &
                have_temperature_field = .true.
            if( have_option( '/material_phase[' // int2str( istate - 1 ) // ']/scalar_field::Concentration' ) ) &
                have_concentration_field = .true.
            if( have_option( '/material_phase[' // int2str( istate - 1 ) // ']/is_multiphase_component' ) ) &
                have_component_field = .true.
            !!$
            call Calculate_All_Rhos( state, packed_state, Mdims )
            if( have_component_field ) then
                call get_option( '/numerical_methods/Max_compositional_its', NonLinearIteration_Components, default = 1 )
            end if
        end do

        if( have_option( '/material_phase[0]/multiphase_properties/capillary_pressure' ) ) &
            have_extra_DiffusionLikeTerm = .true.
        if ( have_option( '/mesh_adaptivity/hr_adaptivity' ) ) then
            call allocate( metric_tensor, extract_mesh(state(1), topology_mesh_name), 'ErrorMetric' )
        end if

        if (is_porous_media) then
            !Check that all the elements have permeability defined
            perm_field => extract_tensor_field(packed_state,"Permeability")
            auxR = minval(perm_field%val, MASK = perm_field%val> 1e-50)
            if (abs(log(maxval(perm_field%val)/auxR) )/log(10.)> 12.) then
              if(getprocno() == 1) then
                ewrite(0,*) "WARNING: The code supports double precision and your Permeability field is trying to span more than 12 orders of magnitude."
                ewrite(0,*) "This may lead to convergence errors. Check your permeability values."
              end if
            end if
            do its = 1, size(perm_field%val,3)
                do its2 = 1, size(perm_field%val,2)
                    if (perm_field%val(its2,its2,its)< 1e-30) then
                        FLExit( "Some elements do not have permebility defined or it is zero. Use a very small value instead, currently the minimum is 1e-30." )
                    end if
                end do
            end do

            !Get into packed state relative permeability, immobile fractions, ...
            call get_RockFluidProp(state, packed_state, Mdims, ndgln, current_time)
            !Allocate the memory to obtain the sigmas at the interface between elements
            call allocate_porous_adv_coefs(Mdims, upwnd)
            !Ensure that the initial condition for the saturation sum to 1.
            call Initialise_Saturation_sums_one(Mdims, packed_state, .true.)
        end if


        !!$ Starting Time Loop
        itime = 0
        ! if this is not a zero timestep simulation (otherwise, there would
        ! be two identical dump files)! unless explicitly disabled
        if( current_time < finish_time .and. &
         .not. have_option("/io/disable_dump_at_start") ) then
!-------------------------------------------------------------------------------
! To allow checkpointing at the 0 timestep - taken from later in the subroutine (find write_state)
             if (do_checkpoint_simulation(dump_no)) then
                  checkpoint_number=0
                  call checkpoint_simulation(state,cp_no=checkpoint_number,&
                                               protect_simulation_name=.true.,file_type='.mpml')
             end if
             !not_to_move_det_yet = .false. ;
!-------------------------------------------------------------------------------
              call write_state(dump_no, state)

        end if
        !Initialise FPI_eq_taken
        FPI_eq_taken = 0
        if(have_option("/io/stat/output_at_start")) then
            call write_diagnostics(state, current_time, dt,&
                timestep, not_to_move_det_yet=.true., non_linear_iterations = FPI_eq_taken)
        end if
        ! When outlet_id is allocated, calculate_flux is true and we want to calculate outfluxes
        ! If calculating boundary fluxes, allocate and initialise to zero outfluxes variables
        if (outfluxes%calculate_flux) call allocate_multi_outfluxes(Mdims, outfluxes)
!       Allocate memory and initialise calculate_mass_global if calculate_mass_flag is switched on to store the total mass change in the domain
        allocate(calculate_mass_delta(Mdims%nphase,2))
        calculate_mass_delta(:,:) = 0.0

        checkpoint_number=1

        !Prepapre the pipes
        if (Mdims%npres > 1) then
           !Retrieve the elements with pipes and the corresponding coordinates
           call retrieve_pipes_coords(state, packed_state, Mdims, ndgln, eles_with_pipe)
        end if

        call petsc_logging(3,stages,ierrr,default=.true.)
        call petsc_logging(1,stages,ierrr,default=.true.)
        ! call petsc_logging(1,stages,ierrr,default=.false., push_no=1, stage_name="PRELIM")
        ! call petsc_logging(1,stages,ierrr,default=.false., push_no=2, stage_name="FORCE")
        ! call petsc_logging(1,stages,ierrr,default=.false., push_no=3, stage_name="SATURATION")
        ! call petsc_logging(1,stages,ierrr,default=.false., push_no=4, stage_name="TEMP")
        ! call petsc_logging(1,stages,ierrr,default=.false., push_no=5, stage_name="COMP")
        ! call petsc_logging(1,stages,ierrr,default=.false., push_no=6, stage_name="DT+VTU ")
        ! call petsc_logging(1,stages,ierrr,default=.false., push_no=7, stage_name="ADAPT")
        ! call petsc_logging(1,stages,ierrr,default=.false., push_no=8, stage_name="REST")
        call petsc_logging(2,stages,ierrr,default=.true., push_no=1)

        !!$ Time loop
        Loop_Time: do
            ewrite(2,*) '    NEW DT', itime+1

            ! initialise the porous media model if needed. Simulation will stop once gravity capillary equilibration is reached
            !Prepapre the pipes
            if (Mdims%npres > 1) call initialize_pipes_package_and_gamma(state, packed_state, pipes_aux, Mdims, Mspars, ndgln)!Re-read pipe properties such as gamma
            if (is_porous_media) then
              call initialise_porous_media(Mdims, ndgln, packed_state, state, exit_initialise_porous_media)
              if (exit_initialise_porous_media) exit Loop_Time
            end if
            !Check first time step
            sum_theta_flux_j = 1. ; sum_one_m_theta_flux_j = 0.

            ! Adapt mesh within the FPI?
            adapt_mesh_in_FPI = have_option( '/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI')
            if (adapt_mesh_in_FPI) call get_option('/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI', Courant_tol, default = -1.)


            itime = itime + 1
            timestep = itime
            call get_option( '/timestepping/timestep', dt )
            call get_option( '/timestepping/current_time', acctim )
            old_acctim = acctim
            acctim = acctim + dt
            current_time = acctim
            call set_option( '/timestepping/current_time', acctim )
            new_lim = .true.
            ! Added a tolerance of 0.001dt to the condition below that stops us exiting the loop before printing the last time step.
            if ( acctim > finish_time + 0.001*dt ) then
              ewrite(1,*) "Passed final time"
              exit Loop_Time
            end if
            call get_option( '/timestepping/final_timestep', final_timestep, stat )
            if ( stat == spud_no_error ) then
                if ( itime > final_timestep ) then
                    ewrite(1,*) "Passed final timestep"
                    exit Loop_Time
                end if
            end if
            !########DO NOT MODIFY THE ORDERING IN THIS SECTION AND TREAT IT AS A BLOCK#######

            !!$ Start non-linear loop
            first_nonlinear_time_step = .true.
            its = 1

            !Store backup to be able to repeat a timestep
            if (nonLinearAdaptTs) call Adaptive_NonLinear(mdims, packed_state, reference_field, its, &
                Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs, old_acctim, 1)

            !! Update all fields from time-step 'N - 1'
            call copy_packed_new_to_old( packed_state )
            ExitNonLinearLoop = .false.
            porosity_field=>extract_vector_field(packed_state,"Porosity")
            ! evaluate prescribed fields at time = current_time+dt
            call set_prescribed_field_values( state, exclude_interpolated = .true., &
                exclude_nonreprescribed = .true., time = acctim )

            !Initialise to zero the SFPI counter
            SFPI_taken = 0
            !########DO NOT MODIFY THE ORDERING IN THIS SECTION AND TREAT IT AS A BLOCK#######
            Loop_NonLinearIteration: do  while (its <= NonLinearIteration)
                print *, 'its', its
                !##########################Impose tunneled BCs############################################################
                !We impose it once per time-level with the exception of the first time-level where we do it after the second non-linear loop
                !This is because we need outfluxes to contain data
                if ((itime > 1 .and. its == 1 ).or. (itime == 1 .and. its== 2)) &
                                call Impose_connected_BCs(outfluxes, packed_state, Mdims, acctim )
                !##########################End tunneled BCs############################################################



                !#=================================================================================================================
                !# Vinicius: Exit simulation if it do not reach convergence
                !#=================================================================================================================
                ! call get_option( "/numerical_methods/max_sat_its", max_sat_its, default = 9)
                ! if (its == NonLinearIteration .and. SFPI_its >= max_sat_its) exit Loop_Time
                !#=================================================================================================================
                !# Vinicius-end: Exit simulation if it do not reach convergence
                !#=================================================================================================================

              !for the diagnostic field, now it seems to be working fine...
                ewrite(2,*) '  NEW ITS', its
                !if adapt_mesh_in_FPI, relax the convergence criteria, since we only want the approx position of the flow
                if (adapt_mesh_in_FPI) call adapt_mesh_within_FPI(ExitNonLinearLoop, adapt_mesh_in_FPI, its, 1)

                !Store the field we want to compare with to check how are the computations going
                call Adaptive_NonLinear(Mdims, packed_state, reference_field, its, &
                    Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs, old_acctim, 2)
                call Calculate_All_Rhos( state, packed_state, Mdims )

                if ( is_porous_media ) then
                    call Calculate_PorousMedia_AbsorptionTerms( Mdims%nphase, state, packed_state, multi_absorp%PorousMedia, Mdims, &
                        CV_funs, CV_GIdims, Mspars, ndgln, upwnd, suf_sig_diagten_bc )
                end if

                ScalarField_Source_Store = 0.0
                if ( Mdims%ncomp > 1 ) then
                   PhaseVolumeFractionComponentSource => extract_tensor_field(packed_state,"PackedPhaseVolumeFractionComponentSource")
                   ScalarField_Source_Store = PhaseVolumeFractionComponentSource%val(1,:,:)
                end if
                PhaseVolumeFractionSource => extract_tensor_field(packed_state,"PackedPhaseVolumeFractionSource", stat)
                if ( stat == 0 ) ScalarField_Source_Store = ScalarField_Source_Store + PhaseVolumeFractionSource%val(1,:,:)

                PythonPhaseVolumeFractionSource => extract_vector_field(state(1),"VSource", python_stat)
                if ( python_stat == 0 ) ScalarField_Source_Store =  PythonPhaseVolumeFractionSource%val(:,:)

                Mdisopt%volfra_use_theta_flux = Mdims%ncomp > 1

                !#=================================================================================================================
                !# Andreas. Here we find if we have asked for a requested_cfl_pressure (rcp_)
                !#          That meens that the code will skip the pressure solve for every rcp (eg rcp=3) time steps.
                !#          We check the input value has an approprate value and if not assigns the default
                !#=================================================================================================================
                if ( have_option( '/timestepping/adaptive_timestep/cfl_pressure' ) ) &
                  call EnterForceBalanceEquation(EnterSolve, its, itime, acctim, t_adapt_threshold, after_adapt, after_adapt_itime, Courant_number(1))

                !#=================================================================================================================

                call petsc_logging(3,stages,ierrr,default=.true.)
                call petsc_logging(2,stages,ierrr,default=.true., push_no=2)

                !#=================================================================================================================
                !# Andreas. I took the velocity and pressure_fields out of the Conditional_ForceBalanceEquation, to always update
                !#=================================================================================================================
                call set_nu_to_u( packed_state )

                velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
                pressure_field=>extract_tensor_field(packed_state,"PackedFEPressure")

                !#=================================================================================================================
                !# Andreas. I added a flag in the Conditional_ForceBalanceEquation to eiher enter or not.
                !#    TODO. This has to be updated with adaptivity as well.
                !#=================================================================================================================
                !$ Now solving the Momentum Equation ( = Force Balance Equation )
                Conditional_ForceBalanceEquation: if ( solve_force_balance .and. EnterSolve ) then
                    !if (getprocno() == 1 .and. its==1) print*, "Time step is:", itime
                    print *, 'calling FORCE_BAL_CTY_ASSEM_SOLVE'
                    CALL FORCE_BAL_CTY_ASSEM_SOLVE( state, packed_state, &
                        Mdims, CV_GIdims, FE_GIdims, CV_funs, FE_funs, Mspars, ndgln, Mdisopt, &
                        Mmat,multi_absorp, upwnd, eles_with_pipe, pipes_aux, velocity_field, pressure_field, &
                        dt, SUF_SIG_DIAGTEN_BC, ScalarField_Source_Store, Porosity_field%val, &
                        igot_theta_flux, sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j,&
                        calculate_mass_delta, outfluxes, pres_its_taken, its, Courant_number)
                end if Conditional_ForceBalanceEquation

                call petsc_logging(3,stages,ierrr,default=.true.)
                call petsc_logging(2,stages,ierrr,default=.true., push_no=3)
                !#=================================================================================================================
                !# End Pressure Solve -> Move to -> Saturation
                !#=================================================================================================================

                !#=================================================================================================================
                !# Initialisation of PHREEQC - must be after FORCE_BAL_CTY_ASSEM_SOLVE. Needs a field which is calculate only there.
                !#=================================================================================================================

                if (after_adapt .or. (itime == 1 .and. its == 1)) then
                  if (have_option("/porous_media/Phreeqc_coupling"))then
#ifdef USING_PHREEQC
                    call init_PHREEQC(Mdims, packed_state, phreeqc_id, concetration_phreeqc, after_adapt)
#else
                    FLAbort( "PHREEQC coupling option activated by the link to PHREEQRM is not activated." )
#endif
                  end if
                end if
                Conditional_PhaseVolumeFraction: if ( solve_PhaseVolumeFraction ) then
                    print *, 'calling solve_phaseVolumeFraction'
                    call VolumeFraction_Assemble_Solve( state, packed_state, multicomponent_state,&
                        Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, &
                        Mmat, multi_absorp, upwnd, eles_with_pipe, pipes_aux, dt, SUF_SIG_DIAGTEN_BC, &
                        ScalarField_Source_Store, Porosity_field%val, igot_theta_flux, mass_ele, its, itime, SFPI_taken, SFPI_its, Courant_number, &
                        sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j)

                end if Conditional_PhaseVolumeFraction
              call petsc_logging(3,stages,ierrr,default=.true.)
              call petsc_logging(2,stages,ierrr,default=.true., push_no=4)
                !#=================================================================================================================
                !# End Saturation -> Move to -> Velocity Update
                !#=================================================================================================================

                !!$ Calculate Darcy velocity with the most up-to-date information
                if(is_porous_media) call get_DarcyVelocity( Mdims, ndgln, state, packed_state, upwnd )
                !#=================================================================================================================
                !# End Velocity Update -> Move to ->the rest
                !#=================================================================================================================

                !!$ Solve advection of the scalar 'Temperature':
                Conditional_ScalarAdvectionField: if( have_temperature_field ) then

                    ewrite(3,*)'Now advecting Temperature Field'
                    call set_nu_to_u( packed_state )
                    !call calculate_diffusivity( state, packed_state, Mdims, ndgln, ScalarAdvectionField_Diffusion )
                    tracer_field=>extract_tensor_field(packed_state,"PackedTemperature")
                    velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
                    density_field=>extract_tensor_field(packed_state,"PackedDensity",stat)
                    saturation_field=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
                    call Calculate_All_Rhos( state, packed_state, Mdims, get_RhoCp = .true. )
                    call INTENERGE_ASSEM_SOLVE( state, packed_state, &
                        Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,upwnd,&
                        tracer_field,velocity_field,density_field, multi_absorp, dt, &
                        suf_sig_diagten_bc, Porosity_field%val, &
                        !!$
                        0, igot_theta_flux, Mdisopt%t_get_theta_flux, Mdisopt%t_use_theta_flux, &
                        THETA_GDIFF, eles_with_pipe, pipes_aux, &
                        option_path = '/material_phase[0]/scalar_field::Temperature', &
                        thermal = .true.,&
                        ! thermal = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/equation::InternalEnergy'),&
                        saturation=saturation_field, nonlinear_iteration = its)
                END IF Conditional_ScalarAdvectionField

                sum_theta_flux = 0. ; sum_one_m_theta_flux = 0. ; sum_theta_flux_j = 0. ; sum_one_m_theta_flux_j = 0.

                !#=================================================================================================================
                  call petsc_logging(3,stages,ierrr,default=.true.)
                  call petsc_logging(2,stages,ierrr,default=.true., push_no=5)
                !Solve for components here
                if (have_component_field) then
                     !!$ Calculate Density_Component for compositional
                     if ( have_component_field ) call Calculate_Component_Rho( state, packed_state, Mdims )
                     sum_theta_flux = 0. ; sum_one_m_theta_flux = 0. ; sum_theta_flux_j = 0. ; sum_one_m_theta_flux_j = 0.
                     call Compositional_Assemble_Solve(state, packed_state, multicomponent_state, &
                     Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat, upwnd,&
                     multi_absorp, DT, SUF_SIG_DIAGTEN_BC, &
                     Mdisopt%comp_get_theta_flux, Mdisopt%comp_use_theta_flux,  &
                     theta_gdiff, eles_with_pipe, pipes_aux, mass_ele, &
                     sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j)
                end if

                  call petsc_logging(3,stages,ierrr,default=.true.)
                  call petsc_logging(2,stages,ierrr,default=.true., push_no=6)

                !# End Compositional transport -> Move to -> Analysis of the non-linear convergence
                !#=================================================================================================================

                !# Solving for species fields through PHREEQC
                !#=================================================================================================================

! #ifdef USING_PHREEQC
!                 if (.not. have_option("/porous_media/Phreeqc_coupling/one_way_coupling") .or. its == 1) then
!                       call run_PHREEQC(Mdims, packed_state, id, concetration_phreeqc)
!                 end if
! #endif

                if (have_Active_Tracers) then
                  !We make sure to only enter here once if there are no passive tracers
                  have_Active_Tracers = .false.
                  velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
                  density_field=>extract_tensor_field(packed_state,"PackedDensity",stat)
                  saturation_field=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
                  call set_nu_to_u( packed_state )
                  !Passive fields taken from phase 1, automatically should detect internally if other phases need it
                  fields = option_count("/material_phase[0]/scalar_field")
                  do k = 1, fields
                    call get_option("/material_phase[0]/scalar_field["// int2str( k - 1 )//"]/name",option_name)
                    if (is_Active_Tracer_field(option_name)) then
                      have_Active_Tracers = .true.!OK there are tracers so remember for next time
                      tracer_field=>extract_tensor_field(packed_state,"Packed"//trim(option_name))
                      call Tracer_Assemble_Solve( trim(option_name), state, packed_state, &
                      Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,upwnd,&
                      tracer_field,velocity_field,density_field, multi_absorp, dt, &
                      suf_sig_diagten_bc, Porosity_field%val, &
                      !!$
                      0, igot_theta_flux, Mdisopt%t_get_theta_flux, Mdisopt%t_use_theta_flux, &
                      THETA_GDIFF, eles_with_pipe, pipes_aux, &
                      saturation=saturation_field, nonlinear_iteration = its)
                      nullify(tracer_field)
                    end if
                  end do
                end if
                !#=================================================================================================================



                !Check if the results are good so far and act in consequence, only does something if requested by the user
                if (sig_hup .or. sig_int) then
                    ewrite(1,*) "Caught signal, exiting nonlinear loop"
                    exit Loop_NonLinearIteration
                end if

                !Finally calculate if the time needs to be adapted or not
                call Adaptive_NonLinear(Mdims, packed_state, reference_field, its,&
                    Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs, old_acctim, 3, calculate_mass_delta, &
                        adapt_mesh_in_FPI, Accum_Courant, Courant_tol, Courant_number(2), first_time_step)

                !Flag the matrices as already calculated (only the storable ones
                Mmat%stored = .true.!Since the mesh can be adapted below, this has to be set to true before the adapt_mesh_in_FPI

                if (ExitNonLinearLoop) then
                    if (adapt_mesh_in_FPI) then
                      !Calculate the acumulated COurant Number
                        Accum_Courant = Accum_Courant + Courant_number(2)
                        if (Accum_Courant >= Courant_tol .or. first_time_step) then
                            call adapt_mesh_within_FPI(ExitNonLinearLoop, adapt_mesh_in_FPI, its, 2)
                            Accum_Courant = 0.
                        else
                            exit Loop_NonLinearIteration
                        end if
                    else
                        exit Loop_NonLinearIteration
                    end if
                end if
                after_adapt=.false.
                its = its + 1
                first_nonlinear_time_step = .false.
            end do Loop_NonLinearIteration


#ifdef USING_PHREEQC
            call run_PHREEQC(Mdims, packed_state, phreeqc_id, concetration_phreeqc)
#endif

            if (have_Passive_Tracers) then
              !We make sure to only enter here once if there are no passive tracers
              have_Passive_Tracers = .false.
              velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
              density_field=>extract_tensor_field(packed_state,"PackedDensity",stat)
              saturation_field=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
              call set_nu_to_u( packed_state )
              !Passive fields taken from phase 1, automatically should detect internally if other phases need it
              fields = option_count("/material_phase[0]/scalar_field")
              do k = 1, fields
                call get_option("/material_phase[0]/scalar_field["// int2str( k - 1 )//"]/name",option_name)
                if (is_PassiveTracer_field(option_name)) then
                  have_Passive_Tracers = .true.!OK there are passive tracers so remember for next time
                  tracer_field=>extract_tensor_field(packed_state,"Packed"//trim(option_name))
                  call Tracer_Assemble_Solve( trim(option_name), state, packed_state, &
                    Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,upwnd,&
                    tracer_field,velocity_field,density_field, multi_absorp, dt, &
                    suf_sig_diagten_bc, Porosity_field%val, &
                    !!$
                    0, igot_theta_flux, Mdisopt%t_get_theta_flux, Mdisopt%t_use_theta_flux, &
                    THETA_GDIFF, eles_with_pipe, pipes_aux, &
                    saturation=saturation_field, nonlinear_iteration = its)
                  nullify(tracer_field)
                end if
              end do
            end if

            !Metal dissolution happens here
            if (have_option("/porous_media/Metal_dissolution_precipitation"))call metal_dissolution_precipitation(state, packed_state, Mdims, ndgln)

            !Flash dissolution happens here JUST BEFORE get_RockFluidProp
            if (have_option("/porous_media/Gas_dissolution"))call flash_gas_dissolution(state, packed_state, Mdims, ndgln)
            !Update immobile fractions values (for hysteresis relperm models)
            !######HAS TO BE the last step of the time-loop
            if (have_option_for_any_phase("/multiphase_properties/immobile_fraction/scalar_field::Land_coefficient",&
            Mdims%n_in_pres)) call get_RockFluidProp(state, packed_state, Mdims, ndgln, update_only = .true.)

            if (have_option( '/io/Show_Convergence') .and. getprocno() == 1) then
              ewrite(1,*) "Iterations taken by the pressure linear solver:", pres_its_taken
            end if
            !Store the combination of Nonlinear iterations performed. Only account of SFPI if multiphase porous media flow
            if (.not. is_porous_media .or. mdims%n_in_pres == 1) SFPI_taken = 0
            FPI_eq_taken = dble(its) + dble(SFPI_taken)/3.!SFPI cost 1/3 roughly, this needs to be revisited when solving for nphases-1

            ! If calculating boundary fluxes, dump them to outfluxes.csv
            if(outfluxes%calculate_flux .and..not.Repeat_time_step) then
                call get_option( '/timestepping/current_time', acctim )
                call dump_outflux(acctim,itime,outfluxes)
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
            current_time = acctim
            call Calculate_All_Rhos( state, packed_state, Mdims )
            !!######################DIAGNOSTIC FIELD CALCULATION TREAT THIS LIKE A BLOCK######################
            !!$ Calculate diagnostic fields (it should be outside the timeloop as a one-way coupling field only)
            call get_option( '/timestepping/timestep', dt); !Backup of dt
            !dT may have changed for the next time level, however for diagnostics we need it as it has been done for this time level
            !therefore we compute it based on the actual difference of time
            call set_option( '/timestepping/timestep', acctim-old_acctim)
            !Now we ensure that the time-step is the correct one
            call set_option( '/timestepping/timestep', dt)
            !Time to compute the self-potential if required
            if (write_all_stats .and. have_option("/porous_media/SelfPotential")) &
                    call Assemble_and_solve_SP(Mdims, state, packed_state, ndgln, Mmat, Mspars, CV_funs, CV_GIdims)
            !Now compute diagnostics
            call calculate_diagnostic_variables( state, exclude_nonrecalculated = .true. )
            !calculate_diagnostic_variables_new <= computes other diagnostics such as python-based fields
            call calculate_diagnostic_variables_new( state, exclude_nonrecalculated = .true. )!sprint_to_do it used to zerod the pressure

            !!######################DIAGNOSTIC FIELD CALCULATION TREAT THIS LIKE A BLOCK######################

            !Now we ensure that the time-step is the correct one
            if (write_all_stats) call write_diagnostics( state, current_time, dt, itime , non_linear_iterations = FPI_eq_taken) ! Write stat file

            if (is_porous_media .and. getprocno() == 1) then
                if (have_option('/io/Courant_number')) then!printout in the terminal
                    ewrite(0,*) "Courant_number and shock-front Courant number", Courant_number
                else !printout only in the log
                    ewrite(1,*) "Courant_number and shock-front Courant number", Courant_number
                end if
            end if

            !Call to create the output vtu files, if required and also checkpoint
            call create_dump_vtu_and_checkpoints()

            call petsc_logging(3,stages,ierrr,default=.true.)
            call petsc_logging(2,stages,ierrr,default=.true., push_no=7)
            ! Call to adapt the mesh if required! If adapting within the FPI then the adaption is controlled elsewhere
            if(acctim >= t_adapt_threshold .and. .not. have_option( '/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI')) then
              call adapt_mesh_mp()
            end if
            ! ####Packing this section inside a internal subroutine breaks the code for non-debugging####
            !!$ Simple adaptive time stepping algorithm

            call petsc_logging(3,stages,ierrr,default=.true.)
            call petsc_logging(2,stages,ierrr,default=.true., push_no=8)
            if ( have_option( '/timestepping/adaptive_timestep' ) ) then
                nonlinear_dt = dt!To use if also the nonlinear adapt time-step is on
                c = -66.6 ; minc = 0. ; maxc = 66.e6 ; ic = 1.1!66.e6
                call get_option( '/timestepping/adaptive_timestep/requested_cfl', rc )
                call get_option( '/timestepping/adaptive_timestep/minimum_timestep', minc, stat, default = 0.)
                call get_option( '/timestepping/adaptive_timestep/maximum_timestep', maxc, stat, default = 66.e6)
                call get_option( '/timestepping/adaptive_timestep/increase_tolerance', ic, stat, default = 1.1)
                !For porous media we need to use the Courant number obtained in cv_assemb
                if (is_porous_media) then
                  !We use stat here as a normal integer
                  stat = 1
                  if (have_option("/timestepping/adaptive_timestep/requested_cfl/Shock_front_CFL")) stat = 2
                  c = max ( c, Courant_number(stat) )
                    ! ewrite(1,*) "maximum cfl number at", current_time, "s =", c
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
                !call get_option( '/timestepping/timestep', dt )
                !To ensure that we always create a vtu file at the desired time (if requested),
                !we control the maximum time-step size to ensure that at some point the ts changes to provide that precise time
                !Original solution slowed down simulations due to having to build up dt again after forced reduction, now fixed by using stored_dt when appropiate
                if (have_option('/io/dump_period')) then
                    maxc = max(min(maxc, abs(current_time - mdims%init_time - dump_period*dump_no)), 1d-15)
                    ! Make sure we dump at the required time and we don't get dt = 0
                    ! Storing current dt before reduction by period_dump when necessary, so we can go back to it after dump
                    if (dt>maxc) then
                        stored_dt=dt
                    end if
				            ! Checking if previous time step was reduced (dt) for meeting dump_period requirement
                    if (ic<stored_dt/dt .and. dt>0) then
                        ! If so, change increase/decrease dt tolerance (so it can catch up faster on dt-before-reduction-by-period-dump)
                        ic=stored_dt/dt
                    end if
                    ! Set stored_dt for normal case to -1 so ic is set as default/in diamond
                    if (dt<=maxc) then
                        stored_dt=-1
                    end if
                end if
                dt = min( dt * rc / c, ic * dt )
                if(have_option("/timestepping/adaptive_timestep/minimum_timestep/dump_vtu_if_reached")) then
                  if(dt < minc) then
                    ewrite(0, *) "Minimum timestep reached - outputting vtu"
                    !SIG_INT = .true.
                    call write_state_mindt(dump_no, acctim, state)
                  end if
                end if
                if(have_option("/timestepping/adaptive_timestep/minimum_timestep/dump_vtu_if_reached/terminate")) then
                  if(dt < minc) then
                    ewrite(0, *) "Minimum timestep reached - terminating"
                    SIG_INT = .true.
                    ! call write_state(dump_no, state)
                  end if
                end if
                dt = max( min( dt , maxc ), minc )
                !If we are imposing also adaptive time-step based on the stability of the non-linear solver, use the CFL as an upper bound
                if (have_option( '/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/adaptive_timestep_nonlinear/')) then
                  dt = min( dt , nonlinear_dt )
                end if
                ! dt = max( min( min( dt * rc / c, ic * dt ), maxc ), minc ) Original
                !Make sure we finish at required time and we don't get dt = 0
                dt = max(min(dt, finish_time - current_time), 1d-8); call allmin(dt)
                if (current_time+dt>=finish_time) exit Loop_Time
                call set_option( '/timestepping/timestep', dt )
            end if
            ! ####UP TO HERE####

            !Post processing if the mesh has been adapted or to recalculate fields for the new time-level
            if ( do_reallocate_fields ) then
                after_adapt=.true.
            else
                after_adapt=.false.
            end if
            if ( after_adapt .and. have_option( '/mesh_adaptivity/hr_adaptivity/nonlinear_iterations_after_adapt' ) ) then
                call get_option( '/mesh_adaptivity/hr_adaptivity/nonlinear_iterations_after_adapt', NonLinearIteration )
            else
                call get_option( '/solver_options/Non_Linear_Solver', NonLinearIteration, default = 3 )
            end if
            call set_boundary_conditions_values(state, shift_time=.false.)
            if (sig_hup .or. sig_int) then
                ewrite(1,*) "Caught signal, exiting"
                exit Loop_Time
            end if
            first_time_step = .false.

        end do Loop_Time

!#=================================================================================================================
!# Vinicius: Free XGB model
!#=================================================================================================================
#ifdef USING_XGBOOST
        if (getprocno() == 1) call xgboost_free_model()
#endif
!#=================================================================================================================
!# Vinicius-end: Free XGB model
!#=================================================================================================================

        if (has_references(metric_tensor)) call deallocate(metric_tensor)
        !!$ Now deallocating arrays:
        deallocate( &
            !!$ Defining element-pair type and discretisation options and coefficients
            !!$ Working arrays
            theta_gdiff, ScalarField_Source_Store, &
            mass_ele,&
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
        call deallocate_multi_shape_funs(CV_funs)
        call deallocate_multi_shape_funs(FE_funs)
        call deallocate_multi_ndgln(ndgln)
        call deallocate_multi_sparsities(Mspars)
        call destroy_multi_matrices(Mmat)
        call deallocate_porous_adv_coefs(upwnd)
        call deallocate_multi_absorption(multi_absorp, .true.)
        call deallocate_multi_pipe_package(pipes_aux)
#ifdef USING_PHREEQC
        call deallocate_PHREEQC(phreeqc_id)
#endif

        !***************************************
        ! INTERPOLATION MEMORY CLEANUP
        if (numberfields_CVGalerkin_interp > 0) then
            call MemoryCleanupInterpolation1()     ! Clean up state_old allocations here.
                                                   ! State_new cleanup happens straight after calling the interpolation with flag = 1
                                                   ! inside the adaptivity loop. Probably best to split M2Minterpolation into
                                                   ! separate flag == 0, flag == 1 subroutines and deallocate inside them
                                                   ! (future work).
        endif
        !***************************************
        if (outfluxes%calculate_flux) call destroy_multi_outfluxes(outfluxes)



        return

    contains
        !> routine puts various CSR sparsities into packed_state
        subroutine put_CSR_spars_into_packed_state()

            use sparse_tools
            type(csr_sparsity), pointer :: sparsity
            type(tensor_field), pointer :: tfield
            type(scalar_field), pointer :: sfield
            integer ic, stat
            type(halo_type), pointer :: halo
            type(mesh_type), pointer :: ph_mesh
            ! This creates sparsity for the Hydrostatic Pressure Solver
            ! It should work now in parallel
            ph_mesh => extract_mesh( state( 1 ), "HydrostaticPressure", stat )
            if (stat == 0) then
                allocate( sparsity )
                sfield => extract_scalar_field(state(1),"HydrostaticPressure")
                if (associated( sfield%mesh%halos)) then
                    sparsity=wrap( Mspars%HydrostaticPressure%fin, colm = Mspars%HydrostaticPressure%col, name = "phsparsity",&
                        row_halo=sfield%mesh%halos(2),&
                        column_halo=sfield%mesh%halos(2))
                else
                    sparsity = wrap( Mspars%HydrostaticPressure%fin, colm = Mspars%HydrostaticPressure%col, name = "phsparsity" )
                end if
                call insert( packed_state, sparsity, "phsparsity" )
                call deallocate( sparsity )
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
        end subroutine put_CSR_spars_into_packed_state



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

        !> @brief In this internal subrotuine the generation of output vts and checkpoint are performed
        subroutine create_dump_vtu_and_checkpoints()

            ! sprint to do: check this routine - why are we swapping pressure fields??
            !!$ Write outputs (vtu and checkpoint files)
            if (have_option('/io/dump_period_in_timesteps')) then
                ! dump based on the prescribed period of time steps
                Conditional_Dump_TimeStep: if( ( mod( itime, dump_period_in_timesteps ) == 0 ) ) then
                    if (do_checkpoint_simulation(dump_no)) then
                        call checkpoint_simulation(state,cp_no=checkpoint_number,&
                            protect_simulation_name=.true.,file_type='.mpml')
                        checkpoint_number=checkpoint_number+1
                    end if
                    call get_option( '/timestepping/current_time', current_time ) ! Find the current time
                    if (.not. write_all_stats) then
                        if (have_option("/porous_media/SelfPotential")) &
                            call Assemble_and_solve_SP(Mdims, state, packed_state, ndgln, Mmat, Mspars, CV_funs, CV_GIdims)
                        call write_diagnostics( state, current_time, dt, itime/dump_period_in_timesteps , non_linear_iterations = FPI_eq_taken)  ! Write stat file
                    end if
                    not_to_move_det_yet = .false. ;
                    call write_state( dump_no, state ) ! Now writing into the vtu files
                end if Conditional_Dump_TimeStep
            else if (have_option('/io/dump_period')) then
                ! dump based on the prescribed period of real time
                Conditional_Dump_RealTime: if( (abs(current_time-mdims%init_time - dump_period*dump_no) < 1d-12 .or. current_time-mdims%init_time >= dump_period*dump_no)&
                    .and. current_time-mdims%init_time/=finish_time) then
                    if (do_checkpoint_simulation(dump_no)) then
                        call checkpoint_simulation(state,cp_no=checkpoint_number,&
                            protect_simulation_name=.true.,file_type='.mpml')
                        checkpoint_number=checkpoint_number+1
                    end if
                    if (.not. write_all_stats) then
                        if (have_option("/porous_media/SelfPotential")) &
                            call Assemble_and_solve_SP(Mdims, state, packed_state, ndgln, Mmat, Mspars, CV_funs, CV_GIdims)
                        call write_diagnostics( state, current_time, dt, itime/dump_period_in_timesteps , non_linear_iterations = FPI_eq_taken)  ! Write stat file
                    end if
                    not_to_move_det_yet = .false. ;
                    !Time to compute the self-potential if required
                    if (have_option("/porous_media/SelfPotential")) call Assemble_and_solve_SP(Mdims, state, packed_state, ndgln, Mmat, Mspars, CV_funs, CV_GIdims)
                    call write_state( dump_no, state ) ! Now writing into the vtu files
                end if Conditional_Dump_RealTime
            end if
        end subroutine create_dump_vtu_and_checkpoints

        !>This subroutine performs all the necessary steps to adapt the mesh and create new memory
        subroutine adapt_mesh_mp()
            !local variables
            type( scalar_field ), pointer :: s_field, s_field2, s_field3
            integer                       :: idim
            real, dimension(2)            :: min_max_limits_before, concentration_min_max_limits_before
            type (tensor_field), pointer  :: tempfield, ComponentMassFraction
            type (tensor_field), pointer  :: Concentration
            character(len=PYTHON_FUNC_LEN) :: pyfunc
            real :: acctim

            if (numberfields_CVGalerkin_interp > 0) then ! If there is at least one instance of CVgalerkin then apply the method
                if (have_option('/mesh_adaptivity')) then ! Only need to use interpolation if mesh adaptivity switched on
                    call M2MInterpolation(state, packed_state, Mdims, CV_GIdims, CV_funs, Mspars%small_acv%fin, Mspars%small_acv%col , 0)
                endif
            endif

            !If has_temperature then we want to ensure than when adapting the mesh the field is between bounds
            if (has_temperature) then
                tempfield => extract_tensor_field( packed_state, "PackedTemperature" )
                min_max_limits_before(1) = minval(tempfield%val); call allmin(min_max_limits_before(1))
                min_max_limits_before(2) = maxval(tempfield%val); call allmax(min_max_limits_before(2))
            end if
            if (has_concentration) then
                concentration => extract_tensor_field( packed_state, "PackedConcentration" )
                concentration_min_max_limits_before(1) = minval(concentration%val); call allmin(concentration_min_max_limits_before(1))
                concentration_min_max_limits_before(2) = maxval(concentration%val); call allmax(concentration_min_max_limits_before(2))
                !sprint_to_do: check the boundedness of the mass fraction (and also the temperature) field
                !concentration_min_max_limits_before(1) = 0.; call allmin(concentration_min_max_limits_before(1))
                !concentration_min_max_limits_before(2) = 1.0; call allmax(concentration_min_max_limits_before(2))
            end if
            do_reallocate_fields = .false.

            call get_option( '/timestepping/current_time', acctim )

            Conditional_Adaptivity_ReallocatingFields: if( have_option( '/mesh_adaptivity/hr_adaptivity') ) then
                if( have_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps') ) then
                  if( have_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps/constant') ) then
                    call get_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps/constant', &
                        adapt_time_steps, default=5 )
                  else if ( have_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps/python') ) then
                    call get_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps/python', pyfunc)
                    call integer_from_python(pyfunc, acctim, adapt_time_steps)
                  end if
                else if (have_option( '/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI')) then
                    do_reallocate_fields = .true.
                    adapt_time_steps = 5!adapt_time_steps requires a default value
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
                            itime, not_to_move_det_yet = .true. , non_linear_iterations = FPI_eq_taken)
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
#ifdef USING_PHREEQC
                call deallocate_PHREEQC(phreeqc_id)
#endif
                call deallocate(multiphase_state)
                call deallocate(multicomponent_state )
                !call unlinearise_components()
                deallocate(multiphase_state)
                deallocate(multicomponent_state)
                call deallocate_projection_matrices(CV_funs)
                call deallocate_projection_matrices(FE_funs)
                !!$ Compute primary scalars used in most of the code
                call Get_Primary_Scalars_new( state, Mdims )
                !Check if the user wants to store the outfluxes
                call initialize_multi_outfluxes(outfluxes)
                call pack_multistate(Mdims%npres,state,packed_state,&
                    multiphase_state,multicomponent_state)
                call prepare_absorptions(state, Mdims, multi_absorp)
                call set_boundary_conditions_values(state, shift_time=.false.)
                !!$ Deallocating array variables:
                deallocate( &
                    !!$ Working arrays
                    suf_sig_diagten_bc, &
                    theta_gdiff, ScalarField_Source_Store, &
                    mass_ele, &
                    sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j )
                !Deallocate sparsities
                call deallocate_multi_sparsities(Mspars)
                !Deallocate ndgln
                call deallocate_multi_ndgln(ndgln)
                !Destroy what remains of the matrices
                call destroy_multi_matrices(Mmat)
                !!$ Calculating Global Node Numbers
                call allocate_multi_ndgln(ndgln, Mdims)
                call Compute_Node_Global_Numbers(state, ndgln)
                !!$
                !!$ Computing Sparsity Patterns Matrices
                !!$
                !!$ Defining lengths and allocating space for the matrices

                call Defining_MaxLengths_for_Sparsity_Matrices( Mdims%ndim, Mdims%nphase, Mdims%totele, Mdims%u_nloc, Mdims%cv_nloc, Mdims%ph_nloc, Mdims%cv_nonods, &
                    mx_nface_p1, mxnele, mx_nct, mx_nc, mx_ncolcmc, mx_ncoldgm_pha, &
                    mx_ncolacv, mx_ncolm, mx_ncolph )
                !!$ Defining element-pair type
                call Get_Ele_Type( Mdims%x_nloc, Mdisopt%cv_ele_type, Mdisopt%p_ele_type, Mdisopt%u_ele_type, &
                    Mdisopt%mat_ele_type, Mdisopt%u_sele_type, Mdisopt%cv_sele_type )
                !Allocate and calculate the sparsity patterns
                call Get_Sparsity_Patterns( state, Mdims, Mspars, ndgln, Mdisopt, mx_ncolacv,&
                    mx_ncoldgm_pha, mx_nct,mx_nc, mx_ncolcmc, mx_ncolm, mx_ncolph, mx_nface_p1 )
                call put_CSR_spars_into_packed_state()

                if (is_porous_media) then
                    call get_RockFluidProp(state, packed_state, Mdims, ndgln)
                    call deallocate_porous_adv_coefs(upwnd)
                    call allocate_porous_adv_coefs(Mdims, upwnd)
                    !Clean the pipes memory if required
                    if (Mdims%npres > 1) then
                        call deallocate_multi_pipe_package(pipes_aux)
                        call retrieve_pipes_coords(state, packed_state, Mdims, ndgln, eles_with_pipe)
                        call initialize_pipes_package_and_gamma(state, packed_state, pipes_aux, Mdims, Mspars, ndgln)
                    end if
                    if (.not. have_option("/numerical_methods/do_not_bound_after_adapt")) then
                      !Ensure that the saturation is physically plausible by diffusing unphysical values to neighbouring nodes
                      call BoundedSolutionCorrections(state, packed_state, Mdims, CV_funs, Mspars%small_acv%fin, Mspars%small_acv%col,"PackedPhaseVolumeFraction", for_sat=.true.)
                    end if
                    call Set_Saturation_to_sum_one(mdims, packed_state, state)!<= just in case, cap unphysical values if there are still some
                end if
                if (.not. have_option("/numerical_methods/do_not_bound_after_adapt")) then
                  if (has_temperature) call BoundedSolutionCorrections(state, packed_state, Mdims, CV_funs, Mspars%small_acv%fin, Mspars%small_acv%col, "PackedTemperature", min_max_limits = min_max_limits_before)
                  if (has_concentration) call BoundedSolutionCorrections(state, packed_state, Mdims, CV_funs, Mspars%small_acv%fin, Mspars%small_acv%col, "PackedConcentration" ,min_max_limits = concentration_min_max_limits_before)
                end if
                ! SECOND INTERPOLATION CALL - After adapting the mesh ******************************
                if (numberfields_CVGalerkin_interp > 0) then
                    if(have_option('/mesh_adaptivity')) then ! This clause may be redundant and could be removed - think this code in only executed IF adaptivity is on
                        call M2MInterpolation(state, packed_state, Mdims, CV_GIdims, CV_funs, &
                                Mspars%small_acv%fin, Mspars%small_acv%col, 1)
                        call MemoryCleanupInterpolation2()
                    endif
                endif
                ! *************************************************************
                !!$ Allocating space for various arrays:
                allocate( &
                    !!$
                    suf_sig_diagten_bc( Mdims%stotel * Mdims%cv_snloc * Mdims%nphase, Mdims%ndim ), &
                    mass_ele( Mdims%totele ) )
                !!$
                suf_sig_diagten_bc=0.
                mass_ele=0.
                !!$
                !!$ Computing shape function scalars
                igot_t2 = 0 ; igot_theta_flux = 0
                if( Mdims%ncomp /= 0 )then
                    igot_t2 = 1 ; igot_theta_flux = 1
                end if
                scvngi_theta = CV_GIdims%scvngi
                ncv_faces = CV_count_faces( Mdims, Mdisopt%cv_ele_type, CV_GIdims)
                ! allocate( sum_theta_flux( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                !     sum_one_m_theta_flux( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                !     sum_theta_flux_j( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                !     sum_one_m_theta_flux_j( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                !     theta_gdiff( Mdims%nphase, Mdims%cv_nonods ), &
                !     ScalarField_Source_Store( Mdims%nphase, Mdims%cv_nonods ) )
                    allocate( sum_theta_flux( Mdims%nphase, ncv_faces*igot_theta_flux ), &
                        sum_one_m_theta_flux( Mdims%nphase, ncv_faces*igot_theta_flux ), &
                        sum_theta_flux_j( Mdims%nphase, ncv_faces*igot_theta_flux ), &
                        sum_one_m_theta_flux_j( Mdims%nphase, ncv_faces*igot_theta_flux ), &
                        theta_gdiff( Mdims%nphase, Mdims%cv_nonods ), &
                        ScalarField_Source_Store( Mdims%nphase, Mdims%cv_nonods ) )

                sum_theta_flux = 1. ; sum_one_m_theta_flux = 0.
                sum_theta_flux_j = 1. ; sum_one_m_theta_flux_j = 0.
                ScalarField_Source_Store=0.
!                allocate(opt_vel_upwind_coefs_new(Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%mat_nonods)); opt_vel_upwind_coefs_new =0.
!                allocate(opt_vel_upwind_grad_new(Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%mat_nonods)); opt_vel_upwind_grad_new =0.

                !SPRINT_TO_DO THIS SHOULD BE DONE LIKE IT IS FOR POROUS MEDIA, DIFFUSING THE FIELD TO ENSURE MASS CONSERVATION
                if( have_option( '/material_phase[' // int2str( Mdims%nstate - Mdims%ncomp ) // &
                    ']/is_multiphase_component/Comp_Sum2One/Enforce_Comp_Sum2One' ) ) then
                    ! Initially clip and then ensure the components sum to unity so we don't get surprising results...
                    ComponentMassFraction  => extract_tensor_field( packed_state, "PackedComponentMassFraction" )
                    ComponentMassFraction % val = min ( max ( ComponentMassFraction % val, 0.0), 1.0)
                    ALLOCATE( RSUM( Mdims%nphase ) )
                    DO CV_INOD = 1, Mdims%cv_nonods
                        DO IPHASE = 1, Mdims%nphase
                            RSUM( IPHASE ) = SUM (ComponentMassFraction % val (:, IPHASE, CV_INOD) )
                        END DO
                        DO IPHASE = 1, Mdims%nphase
                            ComponentMassFraction % val (:, IPHASE, CV_INOD) = ComponentMassFraction % val (:, IPHASE, CV_INOD) / RSUM( IPHASE )
                        END DO
                    END DO
                    DEALLOCATE( RSUM )
                end if

                call Calculate_All_Rhos( state, packed_state, Mdims )
            end if Conditional_ReallocatingFields
        end subroutine adapt_mesh_mp


    !>@brief Update the memory when moving to a new time-level by copying the fields into the OLD fields
    subroutine copy_packed_new_to_old(packed_state)
        type(state_type), intent(inout) :: packed_state
        type(scalar_field), pointer :: sfield, nsfield
        type(vector_field), pointer :: vfield, nvfield
        type(tensor_field), pointer :: tfield, ntfield
        integer :: i
        do i=1,size(packed_state%scalar_fields)
            sfield=>packed_state%scalar_fields(i)%ptr
            if (sfield%name(1:4)=="BAK_") cycle!For adapt within FPI we create these fields that we do not need to update
            if (sfield%name(1:9)=="PackedOld") then
                nsfield=>extract_scalar_field(packed_state,"Packed"//sfield%name(10:))
                sfield%val=nsfield%val
            end if
        end do
        do i=1,size(packed_state%vector_fields)
            vfield=>packed_state%vector_fields(i)%ptr
            if (vfield%name(1:4)=="BAK_") cycle!For adapt within FPI we create these fields that we do not need to update
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
        if (.not. is_porous_media) then
          tfield=>extract_tensor_field(packed_state,"PackedOldCVPressure")
          ntfield=>extract_tensor_field(packed_state,"PackedCVPressure")
          tfield%val=ntfield%val
        end if
    end subroutine copy_packed_new_to_old

    !>@brief Copies the linear velocity U into the non-linear velocity field NU
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

    !>@author : Pablo Salinas
    !>@brief In this subroutine we control all the necessary steps to adapt the mesh within the non-linear loop.
    !> This is useful to ensure that the mesh is always well placed.
    !> There are two steps, relaxation of the tolerance in the first loop of the non-linear solver to get a good guess but cheaply,
    !> and secondly adaptation of the mesh based on the new and old field. For more information see: doi.org/10.1007/s10596-018-9759-z
    subroutine adapt_mesh_within_FPI(ExitNonLinearLoop, adapt_mesh_in_FPI, its, flag)
        implicit none
        logical, intent(inout) :: ExitNonLinearLoop, adapt_mesh_in_FPI
        integer, intent(inout) :: its !> Non-linear iteration
        integer, intent(in)    :: flag!> Controls the steps, either 1 or 2.
        ! Local variables
        type(scalar_field), pointer :: sat1, sat2
        type(vector_field), pointer :: vfield1, vfield2
        integer                     :: iphase, fields, k, i
        real, save                  :: Inf_tol = -1, non_linear_tol = -1
        character( len = option_path_len ) :: option_path, option_name


        if (Inf_tol<0) then
          !Store the settings selected by the user
          call get_option( '/solver_options/Non_Linear_Solver/Fixed_Point_Iteration', non_linear_tol, default = 5e-2)
          call get_option( '/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Infinite_norm_tol',Inf_tol, default = 0.01)
        end if

        select case (flag)
            case (1)
                !Relax the convergence criteria since we don't need much precision at this stage
                !Since non_linear_tol is they key convergence criterion, we use a relative value and we reduce it half-order
                call set_option( '/solver_options/Non_Linear_Solver/Fixed_Point_Iteration', min(non_linear_tol*10., 1e-1))
                if (.not.have_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Infinite_norm_tol'))then
                    !Create the option
                    call copy_option( '/timestepping/timestep/', &
                        '/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Infinite_norm_tol')
                end if
                !10% tolerance provides a good enough result
                call set_option( '/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Infinite_norm_tol',min(Inf_tol*5.,0.1))
        case default
            !Three steps
            !1. Store OldValues into BAK values
            do iphase = 1, Mdims%nphase
                !First scalar prognostic fields
                option_path = "/material_phase["// int2str( iphase - 1 )//"]/scalar_field"
                fields = option_count("/material_phase["// int2str( iphase - 1 )//"]/scalar_field")
                do k = 1, fields
                  call get_option("/material_phase["// int2str( iphase - 1 )//"]/scalar_field["// int2str( k - 1 )//"]/name",option_name)
                  option_path = "/material_phase["// int2str( iphase - 1 )//"]/scalar_field::"//trim(option_name)
                  if (option_name(1:4)=="BAK_") cycle!Not backed fields!
                  if (trim(option_name)=="Density" .or. .not.have_option(trim(option_path)//"/prognostic" )) cycle!Nothing to do for density or non prognostic fields
                  if (have_option(trim(option_path)//"/aliased" )) cycle
                  sat1 => extract_scalar_field( state(iphase),  "Old"//trim(option_name))
                  sat2  => extract_scalar_field( state(iphase), "BAK_"//trim(option_name))
                  sat2%val = sat1%val
                end do
                !Finally velocity as well
                vfield1 => extract_vector_field( state(iphase),  "OldVelocity")
                vfield2 => extract_vector_field( state(iphase),  "BAK_Velocity")
                vfield2%val = vfield1%val
            end do

            !2. Adapt the mesh (bak fields have the same adaptivity options)
            call adapt_mesh_mp()

            call copy_packed_new_to_old(packed_state)!<= if we don't do this, the results are wrong and I don't know why, after we moved to the new schema...
            !3. Copy back from BAK fields to Old fields
            do iphase = 1, Mdims%nphase
                !First scalar prognostic fields
                option_path = "/material_phase["// int2str( iphase - 1 )//"]/scalar_field"
                fields = option_count("/material_phase["// int2str( iphase - 1 )//"]/scalar_field")
                do k = 1, fields
                  call get_option("/material_phase["// int2str( iphase - 1 )//"]/scalar_field["// int2str( k - 1 )//"]/name",option_name)
                  option_path = "/material_phase["// int2str( iphase - 1 )//"]/scalar_field::"//trim(option_name)
                  if (option_name(1:4)=="BAK_") cycle!Not backed fields!
                  if (trim(option_name)=="Density" .or. .not.have_option(trim(option_path)//"/prognostic" )) cycle!Nothing to do for density or non prognostic fields
                  if (have_option(trim(option_path)//"/aliased" )) cycle
                  sat1 => extract_scalar_field( state(iphase),  "Old"//trim(option_name))
                  sat2  => extract_scalar_field( state(iphase), "BAK_"//trim(option_name))
                  sat1%val = sat2%val
                end do
                !Finally velocity as well
                vfield1 => extract_vector_field( state(iphase),  "OldVelocity")
                vfield2 => extract_vector_field( state(iphase),  "BAK_Velocity")
                vfield2%val = vfield1%val
            end do
            !Pointing to porosity again is required
            porosity_field=>extract_vector_field(packed_state,"Porosity")
            !Now we have to converge again within the same time-step
            ExitNonLinearLoop = .false.; its = 1
            adapt_mesh_in_FPI = .false.

            !Set the original convergence criteria since we are now solving the equations
            call set_option( '/solver_options/Non_Linear_Solver/Fixed_Point_Iteration', non_linear_tol)
            call set_option( '/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Infinite_norm_tol',Inf_tol)

        end select

    end subroutine adapt_mesh_within_FPI



 end subroutine MultiFluids_SolveTimeLoop



end module multiphase_time_loop
