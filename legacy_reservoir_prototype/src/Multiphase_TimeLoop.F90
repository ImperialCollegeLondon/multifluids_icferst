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
    use multiphase_fractures
    use Compositional_Terms
    use Copy_Outof_State
    use cv_advection, only : cv_count_faces
    use multiphase_1D_engine
    use boundary_conditions_from_options
    use multi_data_types
    use vtk_interfaces
    use multi_interpolation
    use multi_transport
    use multi_pipes
    use momentum_diagnostic_fields, only: calculate_densities

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
        !!Define shape functions
        type (multi_shape_funs) :: CV_funs, FE_funs
        !!$ Primary scalars
        type(multi_dimensions) :: Mdims
        type(multi_gi_dimensions) :: CV_GIdims, FE_GIdims
        !!$ Node global numbers
        type(multi_ndgln) :: ndgln
        !!$ Sparsity patterns
        type (multi_sparsities) :: Mspars
        !!$ Defining element-pair type and discretisation options and coefficients
        type (multi_discretization_opts) :: Mdisopt
        !!$ Defining the necessary matrices and corresponding RHS
        type (multi_matrices) :: Mmat
        !!$ Defining variables to calculate the sigmas at the interface for porous media
        type (porous_adv_coefs) :: upwnd
        !!$ Variable storing all the absorptions we may need
        type(multi_absorption) :: multi_absorp
        integer :: mx_nface_p1, mx_ncolacv, mxnele, mx_ncoldgm_pha, &
            mx_nct, mx_nc, mx_ncolcmc, mx_ncolm, mx_ncolph
        !!$ Defining time- and nonlinear interations-loops variables
        integer :: itime, dump_period_in_timesteps, final_timestep, &
            NonLinearIteration, NonLinearIteration_Components
        real :: acctim, finish_time, dump_period
        !!$ Defining problem that will be solved
        logical :: have_temperature_field, have_component_field, have_extra_DiffusionLikeTerm, &
            solve_force_balance, solve_PhaseVolumeFraction, simple_black_oil_model
        !!$ Defining solver options
        integer :: velocity_max_iterations, PhaseVolumeFraction_max_iterations
        !!$ Shape function related fields:
        integer :: scvngi_theta, igot_t2, igot_theta_flux
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
            ScalarField_Source_Store
        real, dimension( :, :, : ), allocatable :: &
            Velocity_Absorption, Temperature_Absorption
        real, dimension( :, : ), allocatable ::theta_flux, one_m_theta_flux, theta_flux_j, one_m_theta_flux_j, &
            sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j
        integer :: stat, python_stat, istate, iphase, jphase, icomp, its, its2, cv_nodi, adapt_time_steps, cv_inod
        real, dimension( : ), allocatable :: rsum
        real, dimension(:, :), allocatable :: SUF_SIG_DIAGTEN_BC
        type( scalar_field ), pointer :: cfl, rc_field
        real :: c, rc, minc, maxc, ic
        !Variables for adaptive time stepping based on non-linear iterations
        logical :: nonLinearAdaptTs, Repeat_time_step, ExitNonLinearLoop
        real, dimension(:,:,:), allocatable  :: reference_field
        type( tensor_field ), pointer :: D_s, DC_s, DCOLD_s
        type( tensor_field ), pointer :: MFC_s, MFCOLD_s
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
        type(tensor_field), pointer :: tracer_field, velocity_field, density_field, saturation_field, old_saturation_field   !, tracer_source
        type(tensor_field), pointer :: pressure_field, cv_pressure, fe_pressure, PhaseVolumeFractionSource, PhaseVolumeFractionComponentSource
        type(tensor_field), pointer :: Component_Absorption, perm_field
        type(vector_field), pointer :: positions, porosity_field, MeanPoreCV, PythonPhaseVolumeFractionSource
        type(scalar_field), pointer :: DensitySource, T
        !Variables that are used to define the pipe pos
        type(pipe_coords), dimension(:), allocatable:: eles_with_pipe
        type (multi_pipe_package) :: pipes_aux
        !type(scalar_field), pointer :: bathymetry
        logical, parameter :: write_all_stats=.true.
        ! Variables used for calculating boundary outfluxes. Logical "calculate_flux" determines if this calculation is done. Intflux is the time integrated outflux
        ! Ioutlet counts the number of boundaries over which to calculate the outflux
        integer :: ioutlet
        type (multi_outfluxes) :: outfluxes
        ! Variables used in the CVGalerkin interpolation calculation
        integer, save :: numberfields = -1
        real :: t_adapt_threshold
        !Variables for FPI acceleration for flooding
        ! Calculate_mass_delta to store the change in mass calculated over the whole domain
        real, allocatable, dimension(:,:) :: calculate_mass_delta

!!-Variables related to the detection and correction of bad elements
        real, save :: Max_bad_angle = -1 ! set in timeloop from diamond input
        integer :: i, k
        type(bad_elements), allocatable, dimension(:) :: Quality_list
        real, dimension(:,:), pointer:: X_ALL
        real, dimension(:), allocatable :: quality_table
        integer, dimension(:), allocatable :: diagnostics ! number of bad elements - used to generate a diagnostics table
        integer, dimension(2) :: shape
        logical :: mesh_diagnostics = .false., bad_element = .false. ! print out mesh diagnostics / change properties of bad elements to improve deltaP calculations for bad meshes (with large angles)

!!-Variable to keep track of dt reduction for meeting dump_period requirements
        real, save :: stored_dt = -1

#ifdef HAVE_ZOLTAN
      real(zoltan_float) :: ver
      integer(zoltan_int) :: ierr
      ierr = Zoltan_Initialize(ver)
      assert(ierr == ZOLTAN_OK)
#endif

        ! Check wether we are using the CV_Galerkin method
        numberfields=option_count('/material_phase/scalar_field/prognostic/CVgalerkin_interpolation') ! Count # instances of CVGalerkin in the input file

        ! A SWITCH TO DELAY MESH ADAPTIVITY UNTIL SPECIFIED UNSER INPUT TIME t_adapt_threshold
        call get_option("/mesh_adaptivity/hr_adaptivity/t_adapt_delay", t_adapt_threshold, default = 0.0 )


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
        !Since this is a hack for Flooding, we want to do this before we actually start using the density as the height
        !which depends on the pressure. However, for th initial condition we need to use the density to set up the initial Pressure
        !Therefore, we correct the initial condition for the pressure before anything is modified
        !If it is flooding we impose the initial pressure to match the equation P = gravity * (height+bathymetry)
        !The height is the initial condition of the density
!        if (is_flooding .and. have_option('/material_phase[0]/scalar_field::Temperature')) then
!            density_field => extract_tensor_field( packed_state, "PackedDensity" )!Equivalent to height
!            FE_Pressure=>extract_tensor_field(packed_state,"PackedFEPressure")
!            bathymetry => extract_scalar_field( state(1), "Temperature" )!bathymetry
!            FE_Pressure%val(1,1,:) =  9.80665 * (density_field%val(1,1,:) + bathymetry%val(1))
!        end if
        !Retrieve manning coefficient for flooding, this has to be called just after creating pack_multistate
        if (is_flooding) then
            call get_FloodingProp(state, packed_state)
            if (have_option( "/physical_parameters/gravity/magnitude")) then
                ewrite(1,*) "ERROR: For flooding DO NOT define gravity, it is already defined in the code"
                stop
            end if
        end if
        call set_boundary_conditions_values(state, shift_time=.true.)

        !  Access boundary conditions via a call like
        !  call get_entire_boundary_condition(extract_tensor_field(packed_state,"Packed"//name),["weakdirichlet"],tfield,bc_type_list)
        !  where tfield is type(tensor_field) and bc_type_list is integer, dimension(tfield%dim(1),tfield%dim(2),nonods)
        !  Then values are in tfield%val(1/Mdims%ndim/Mdims%ncomp,Mdims%nphase,nonods)
        !  Type ids are in bc_type_list(1/Mdims%ndim/Mdims%ncomp,Mdims%nphase,Mdims%stotel)
        !
        !  A deallocate tfield when finished!!

        Repeat_time_step = .false.
        nonLinearAdaptTs = have_option('/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear')

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
            mass_ele( Mdims%totele ), &
            )
        !!$
        suf_sig_diagten_bc=0.
        mass_ele=0.
        !!$

        !!$ Calculate diagnostic fields
        call calculate_diagnostic_variables( state, exclude_nonrecalculated = .true. )
        call calculate_diagnostic_variables_new( state, exclude_nonrecalculated = .true. )
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
        allocate( theta_flux( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            one_m_theta_flux( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            theta_flux_j( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            one_m_theta_flux_j( Mdims%nphase, ncv_faces * igot_theta_flux ), &
            sum_theta_flux( Mdims%nphase, ncv_faces * igot_theta_flux ), &
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
        Mmat%CV_pressure = have_option( '/material_phase[0]/scalar_field::Pressure/prognostic/CV_P_matrix' )
        if (.not. Mmat%CV_pressure .and. ((Mdims%ndim==2 .and. Mdims%u_nloc == 4) .or. (Mdims%ndim==3 .and. Mdims%u_nloc == 5))) then
            ewrite(0, *) "WARNING: the only tested element pair using bubble shape functions is the P1DG(BL)P1DG(CV)"
        end if
        !Check if we want to use a compacted mass matrix
        if ((Mmat%CV_pressure .or. have_option('/numerical_methods/simple_mass_matrix')) &
                    .and. is_porous_media .and. Mdims%npres == 1) then
        !sprint_to_do for this to work with wells we need to change the sparsity, but that still needs to be done!
            call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/polynomial_degree', i )
            call get_option( '/geometry/mesh::PressureMesh/from_mesh/mesh_shape/polynomial_degree', k )
            Mmat%compact_PIVIT_MAT = (i == (k - 1))
        end if
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
        if ( have_option('/io/dump_period_in_timesteps') ) then
            call get_option('/io/dump_period_in_timesteps/constant', dump_period_in_timesteps, default = 1)
        elseif ( have_option('/io/dump_period') ) then
            call get_option('/io/dump_period/constant', dump_period, default = 0.01)
        end if
        call get_option( '/timestepping/nonlinear_iterations', NonLinearIteration, default = 3 )
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
        simple_black_oil_model = .false.
        if (have_option( "/physical_parameters/black-oil_PVT_table" )) then
            simple_black_oil_model = is_porous_media .and..not.have_component_field .and. Mdims%nphase == 3
            if (.not. simple_black_oil_model) then
                ewrite(0,*) "WARNING: Black-oil modelling based on PVT tables requires porous media, 3 phases and no components"
            end if
            !Initialize Stock tank oil conditions. sprint_to_do. Is this necessary?
            if (simple_black_oil_model) call extended_Black_Oil(state, packed_state, Mdims, flash_flag = 10)
        end if

        if( have_option( '/material_phase[0]/multiphase_properties/capillary_pressure' ) ) &
            have_extra_DiffusionLikeTerm = .true.
        if ( have_option( '/mesh_adaptivity/hr_adaptivity' ) ) then
            call allocate( metric_tensor, extract_mesh(state(1), topology_mesh_name), 'ErrorMetric' )
        end if

        if (is_porous_media) then
            !Get into packed state relative permeability, immobile fractions, ...
            call get_RockFluidProp(state, packed_state)
            !Allocate the memory to obtain the sigmas at the interface between elements
            call allocate_porous_adv_coefs(Mdims, upwnd)
            !Ensure that the initial condition for the saturation sum to 1.
            call Ensure_initial_Saturation_to_sum_one(Mdims, ndgln, packed_state)
        end if

        !!$ Starting Time Loop
        itime = 0
        if( &
             ! if this is not a zero timestep simulation (otherwise, there would
             ! be two identical dump files)
             current_time < finish_time &! unless explicitly disabled
             .and. .not. have_option("/io/disable_dump_at_start") &
             ) then
            
!-------------------------------------------------------------------------------
! to allow checkpointing at the 0 timestep - taken from later in the subroutine (find write_state)
             if (do_checkpoint_simulation(dump_no)) then
                  checkpoint_number=0
                  CV_Pressure=>extract_tensor_field(packed_state,"PackedCVPressure")
                  FE_Pressure=>extract_tensor_field(packed_state,"PackedFEPressure")
                  pressure_field=>extract_tensor_field(packed_state,"PackedFEPressure")
                  call set(pressure_field,FE_Pressure)
                  call checkpoint_simulation(state,cp_no=checkpoint_number,&
                                               protect_simulation_name=.true.,file_type='.mpml')
                  call set(pressure_field,CV_Pressure)
             end if
             !not_to_move_det_yet = .false. ;
!-------------------------------------------------------------------------------
              call write_state(dump_no, state)

        end if
        if(have_option("/io/stat/output_at_start")) then
            call write_diagnostics(state, current_time, dt,&
                timestep, not_to_move_det_yet=.true.)
        end if
        ! When outlet_id is allocated, calculate_flux is true and we want to calculate outfluxes
        ! If calculating boundary fluxes, allocate and initialise to zero outfluxes variables
        if (outfluxes%calculate_flux) call allocate_multi_outfluxes(Mdims, outfluxes)
!       Allocate memory and initialise calculate_mass_global if calculate_mass_flag is switched on to store the total mass change in the domain
        allocate(calculate_mass_delta(Mdims%nphase,2))
        calculate_mass_delta(:,:) = 0.0

        checkpoint_number=1

        ! Reading options for bad_element test
        call BadElementTest(Quality_list, 1)


!!$ Time loop
        Loop_Time: do
            ewrite(2,*) '    NEW DT', itime+1
            ! Tests bad elements and creates table of angles for model elements
            call BadElementTest(Quality_list, 2)

            !Check first time step
            sum_theta_flux_j = 1. ; sum_one_m_theta_flux_j = 0.

            !Prepapre the pipes
            if (Mdims%npres > 1) then
               !Retrieve the elements with pipes and the corresponding coordinates
               if (after_adapt .or. first_time_step) call retrieve_pipes_coords(state, packed_state, Mdims, ndgln, eles_with_pipe)
               call initialize_pipes_package_and_gamma(state, pipes_aux, Mdims, Mspars)
            end if
            ! Adapt mesh within the FPI?
            adapt_mesh_in_FPI = have_option( '/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI')
            if (adapt_mesh_in_FPI) call get_option('/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI', Courant_tol, default = -1.)


            itime = itime + 1
            timestep = itime
            call get_option( '/timestepping/timestep', dt )
            call get_option( '/timestepping/current_time', acctim )
            acctim = acctim + dt
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
            !!$ FEMDEM...
#ifdef USING_FEMDEM
            if ( is_multifracture ) then
               call fracking(packed_state, state,Mdims%nphase)
            end if
#endif
            !########DO NOT MODIFY THE ORDERING IN THIS SECTION AND TREAT IT AS A BLOCK#######
            !!$ Start non-linear loop
            first_nonlinear_time_step = .true.
            its = 1
            !Store backup to be able to repeat a timestep
            if (nonLinearAdaptTs) call Adaptive_NonLinear(packed_state, reference_field, its, &
                Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,1)
            !! Update all fields from time-step 'N - 1'
            call copy_packed_new_to_old( packed_state )
            ExitNonLinearLoop = .false.
            porosity_field=>extract_vector_field(packed_state,"Porosity")
            ! evaluate prescribed fields at time = current_time+dt
            call set_prescribed_field_values( state, exclude_interpolated = .true., &
                exclude_nonreprescribed = .true., time = acctim )
            !Initialize gas molar fraction, this has to occur after copy_packed_new_to_old
            !since for consistency, (later it is called as well) it uses the old values of pressure,
            !however, they have to be the most updated at this point
            if (simple_black_oil_model) call extended_Black_Oil(state, packed_state, Mdims, flash_flag = 0)
            !########DO NOT MODIFY THE ORDERING IN THIS SECTION AND TREAT IT AS A BLOCK#######
            Loop_NonLinearIteration: do  while (its <= NonLinearIteration)
                ewrite(2,*) '  NEW ITS', its
                !if adapt_mesh_in_FPI, relax the convergence criteria, since we only want the approx position of the flow
                if (adapt_mesh_in_FPI) call adapt_mesh_within_FPI(ExitNonLinearLoop, adapt_mesh_in_FPI, its, 1)

                !Store the field we want to compare with to check how are the computations going
                call Adaptive_NonLinear(packed_state, reference_field, its, &
                    Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,2)
                call Calculate_All_Rhos( state, packed_state, Mdims )

                if( solve_force_balance) then
                    if ( is_porous_media ) then
                        call Calculate_PorousMedia_AbsorptionTerms( state, packed_state, multi_absorp%PorousMedia, Mdims, &
                            CV_funs, CV_GIdims, Mspars, ndgln, upwnd, suf_sig_diagten_bc, Quality_list )
                    else if (is_flooding) then
                        call Calculate_flooding_absorptionTerm(state, packed_state, multi_absorp%Flooding, Mdims, ndgln)
                    end if
                end if


!Testing multi_transport
call solve_transport()

                ScalarField_Source_Store = 0.0
                if ( Mdims%ncomp > 1 ) then
                   PhaseVolumeFractionComponentSource => extract_tensor_field(packed_state,"PackedPhaseVolumeFractionComponentSource")
                   ScalarField_Source_Store = PhaseVolumeFractionComponentSource%val(1,:,:)
                end if
                if (is_flooding) then
                    do iphase = 1, size(state)
                        DensitySource => extract_scalar_field(state(iphase),"DensitySource", stat)
                        if ( stat == 0 ) ScalarField_Source_Store(iphase,:) = ScalarField_Source_Store(iphase,:) + DensitySource%val
                    end do
                end if

                PhaseVolumeFractionSource => extract_tensor_field(packed_state,"PackedPhaseVolumeFractionSource", stat)
                if ( stat == 0 ) ScalarField_Source_Store = ScalarField_Source_Store + PhaseVolumeFractionSource%val(1,:,:)

                PythonPhaseVolumeFractionSource => extract_vector_field(state(1),"VSource", python_stat)
                if ( python_stat == 0 ) ScalarField_Source_Store = ScalarField_Source_Store + PythonPhaseVolumeFractionSource%val(:,:)

                Mdisopt%volfra_use_theta_flux = Mdims%ncomp > 1

                !!$ Now solving the Momentum Equation ( = Force Balance Equation )
                Conditional_ForceBalanceEquation: if ( solve_force_balance ) then

                    call set_nu_to_u( packed_state )

                    velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
                    pressure_field=>extract_tensor_field(packed_state,"PackedFEPressure")

                    CALL FORCE_BAL_CTY_ASSEM_SOLVE( state, packed_state, &
                        Mdims, CV_GIdims, FE_GIdims, CV_funs, FE_funs, Mspars, ndgln, Mdisopt, &
                        Mmat,multi_absorp, upwnd, eles_with_pipe, pipes_aux, velocity_field, pressure_field, &
                        dt, &
                        SUF_SIG_DIAGTEN_BC, &
                        ScalarField_Source_Store, Porosity_field%val, &
                        igot_theta_flux, &
                        sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j, &
                        calculate_mass_delta, outfluxes)

                    !!$ Calculate Darcy velocity
                    if(is_porous_media) then
                        !Do not calculate unless necessary, this is not specially efficient...
                        if(have_option('/io/output_darcy_vel').or. is_multifracture) then
                            call get_DarcyVelocity( Mdims, ndgln, packed_state, multi_absorp%PorousMedia )
                        end if
                    end if

                    !!$ Calculate Density_Component for compositional
                    if ( have_component_field ) call Calculate_Component_Rho( state, packed_state, Mdims )

                end if Conditional_ForceBalanceEquation


                Conditional_PhaseVolumeFraction: if ( solve_PhaseVolumeFraction ) then
                    call VolumeFraction_Assemble_Solve( state, packed_state, &
                        Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, &
                        Mmat, multi_absorp, upwnd, eles_with_pipe, pipes_aux, dt, SUF_SIG_DIAGTEN_BC, &
                        ScalarField_Source_Store, Porosity_field%val, &
                        igot_theta_flux, mass_ele, &
                        its, Courant_number, &
                        option_path = '/material_phase[0]/scalar_field::PhaseVolumeFraction', &
                        theta_flux=sum_theta_flux, one_m_theta_flux=sum_one_m_theta_flux, &
                        theta_flux_j=sum_theta_flux_j, one_m_theta_flux_j=sum_one_m_theta_flux_j, Quality_list=Quality_list)
                end if Conditional_PhaseVolumeFraction

                !!$ Solve advection of the scalar 'Temperature':
                Conditional_ScalarAdvectionField: if( have_temperature_field .and. &
                    have_option( '/material_phase[0]/scalar_field::Temperature/prognostic' ) ) then
                    ewrite(3,*)'Now advecting Temperature Field'
                    call set_nu_to_u( packed_state )
                    !call calculate_diffusivity( state, Mdims, ndgln, ScalarAdvectionField_Diffusion )
                    tracer_field=>extract_tensor_field(packed_state,"PackedTemperature")
                    velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
                    density_field=>extract_tensor_field(packed_state,"PackedDensity",stat)
                    saturation_field=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
                    call INTENERGE_ASSEM_SOLVE( state, packed_state, &
                        Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,upwnd,&
                        tracer_field,velocity_field,density_field, multi_absorp, dt, &
                        suf_sig_diagten_bc, &
                        Porosity_field%val, &
                        !!$
                        0, igot_theta_flux, &
                        Mdisopt%t_get_theta_flux, Mdisopt%t_use_theta_flux, &
                        THETA_GDIFF, eles_with_pipe, pipes_aux, &
                        option_path = '/material_phase[0]/scalar_field::Temperature', &
                        thermal = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/equation::InternalEnergy'),&
                        saturation=saturation_field, nonlinear_iteration = its, Courant_number = Courant_number)

                    ! Copy back memory
                    do iphase=1,Mdims%nphase
                       T=>extract_scalar_field(state(iphase),"Temperature")
                       T%val=tracer_field%val(1,iphase,:)
                    end do

                    call Calculate_All_Rhos( state, packed_state, Mdims )
                end if Conditional_ScalarAdvectionField


                sum_theta_flux = 0. ; sum_one_m_theta_flux = 0. ; sum_theta_flux_j = 0. ; sum_one_m_theta_flux_j = 0.

                if ( have_component_field ) call calc_components()

                !Check if the results are good so far and act in consequence, only does something if requested by the user
                if (sig_hup .or. sig_int) then
                    ewrite(1,*) "Caught signal, exiting nonlinear loop"
                    exit Loop_NonLinearIteration
                end if
                !Finally calculate if the time needs to be adapted or not
                call Adaptive_NonLinear(packed_state, reference_field, its,&
                    Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,3, adapt_mesh_in_FPI, calculate_mass_delta)

                !Flag the matrices as already calculated (only the storable ones
                Mmat%stored = .true.!Since the mesh can be adapted below, this has to be set to true before the adapt_mesh_in_FPI

                if (ExitNonLinearLoop) then
                    if (adapt_mesh_in_FPI) then
                        Accum_Courant = Accum_Courant + Courant_number(2)
                        if (Accum_Courant >= Courant_tol .or. first_time_step) then
                            Accum_Courant = 0.
                            call adapt_mesh_within_FPI(ExitNonLinearLoop, adapt_mesh_in_FPI, its, 2)
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

            ! If calculating boundary fluxes, dump them to outfluxes.txt
            if(outfluxes%calculate_flux .and..not.Repeat_time_step) then
                if(getprocno() == 1) call dump_outflux(acctim,itime,outfluxes)
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

            !!$ Calculate diagnostic fields
            call calculate_diagnostic_variables( state, exclude_nonrecalculated = .true. )
            call calculate_diagnostic_variables_new( state, exclude_nonrecalculated = .true. )
            if (write_all_stats) call write_diagnostics( state, current_time, dt, itime ) ! Write stat file

            if (is_porous_media) then
                if (have_option('/io/Courant_number')) then!printout in the terminal
                    ewrite(0,*) "Courant_number and shock-front Courant number", Courant_number
                else!printout only in the log
                    ewrite(1,*) "Courant_number and shock-front Courant number", Courant_number
                end if
            end if


            !Call to create the output vtu files, if required and also checkpoint
            call create_dump_vtu_and_checkpoints()

            ! Call to adapt the mesh if required!
            if(acctim >= t_adapt_threshold) call adapt_mesh_mp()

            ! ####Packing this section inside a internal subroutine breaks the code for non-debugging####
            !!$ Simple adaptive time stepping algorithm
            
            if ( have_option( '/timestepping/adaptive_timestep' ) ) then
                c = -66.6 ; minc = 0. ; maxc = 66.e6 ; ic = 1.1!66.e6
                call get_option( '/timestepping/adaptive_timestep/requested_cfl', rc )
                call get_option( '/timestepping/adaptive_timestep/minimum_timestep', minc, stat )
                call get_option( '/timestepping/adaptive_timestep/maximum_timestep', maxc, stat )
                call get_option( '/timestepping/adaptive_timestep/increase_tolerance', ic, stat )
                !For porous media we need to use the Courant number obtained in cv_assemb
                if (is_porous_media) then
                    c = max ( c, Courant_number(1) )
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
                call get_option( '/timestepping/timestep', dt )
                !To ensure that we always create a vtu file at the desired time (if requested),
                !we control the maximum time-step size to ensure that at some point the ts changes to provide that precise time
                !Original solution slowed down simulations due to having to build up dt again after forced reduction, now fixed by using stored_dt when appropiate
                if (have_option('/io/dump_period')) then
                    maxc = max(min(maxc, abs(current_time - dump_period*dump_no)), 1d-15)
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
                dt = max( min( min( dt * rc / c, ic * dt ), maxc ), minc )
                !Make sure we finish at required time and we don't get dt = 0
                dt = max(min(dt, finish_time - current_time), 1d-15)                    
                
                call allmin(dt)
                call set_option( '/timestepping/timestep', dt )
            end if
            ! ####UP TO HERE####

            !Post processing if the mesh has been adapted or to ecalculate fields for the new time-level
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
            first_time_step = .false.

        end do Loop_Time
        if (has_references(metric_tensor)) call deallocate(metric_tensor)
        !!$ Now deallocating arrays:
        deallocate( &
            !!$ Defining element-pair type and discretisation options and coefficients
            !!$ Working arrays
            theta_gdiff, ScalarField_Source_Store, &
            mass_ele,&
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
        call deallocate_multi_shape_funs(CV_funs)
        call deallocate_multi_shape_funs(FE_funs)
        call deallocate_multi_ndgln(ndgln)
        call deallocate_multi_sparsities(Mspars)
        call destroy_multi_matrices(Mmat)
        call deallocate_porous_adv_coefs(upwnd)
        call deallocate_multi_absorption(multi_absorp, .true.)
        call deallocate_multi_pipe_package(pipes_aux)
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
        if (outfluxes%calculate_flux) call destroy_multi_outfluxes(outfluxes)

        return
    contains

        subroutine put_CSR_spars_into_packed_state()
            !!! routine puts various CSR sparsities into packed_state
            use sparse_tools
            type(csr_sparsity), pointer :: sparsity
            type(tensor_field), pointer :: tfield
            type(scalar_field), pointer :: sfield
            integer ic, stat
            type(halo_type), pointer :: halo
            type(mesh_type), pointer :: ph_mesh
            ! This creates sparsity for the Hydrostatic Pressure Solver
            ! It should work now in parallel
            ph_mesh => extract_mesh( state( 1 ), "ph", stat )
            if (stat == 0) then
                allocate( sparsity )
                sfield => extract_scalar_field(state(1),"Ph")
                if (associated( sfield%mesh%halos)) then
                    sparsity=wrap( Mspars%ph%fin, colm = Mspars%ph%col, name = "phsparsity",&
                        row_halo=sfield%mesh%halos(2),&
                        column_halo=sfield%mesh%halos(2))
                else
                    sparsity = wrap( Mspars%ph%fin, colm = Mspars%ph%col, name = "phsparsity" )
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

        subroutine create_dump_vtu_and_checkpoints()
            !!$ Write outputs (vtu and checkpoint files)
            if (have_option('/io/dump_period_in_timesteps')) then
                ! dump based on the prescribed period of time steps
                Conditional_Dump_TimeStep: if( ( mod( itime, dump_period_in_timesteps ) == 0 ) ) then
                    if (do_checkpoint_simulation(dump_no)) then
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
                    not_to_move_det_yet = .false. ;

                    call write_state( dump_no, state ) ! Now writing into the vtu files
                end if Conditional_Dump_TimeStep
            else if (have_option('/io/dump_period')) then
                ! dump based on the prescribed period of real time
                Conditional_Dump_RealTime: if( (abs(current_time - dump_period*dump_no) < 1d-12 .or. current_time >= dump_period*dump_no)&
                    .and. current_time/=finish_time) then
                    if (do_checkpoint_simulation(dump_no)) then
                        CV_Pressure=>extract_tensor_field(packed_state,"PackedCVPressure")
                        FE_Pressure=>extract_tensor_field(packed_state,"PackedFEPressure")
                        call set(pressure_field,FE_Pressure)
                        call checkpoint_simulation(state,cp_no=checkpoint_number,&
                            protect_simulation_name=.true.,file_type='.mpml')
                        checkpoint_number=checkpoint_number+1
                        call set(pressure_field,CV_Pressure)
                    end if
                    if (.not. write_all_stats)call write_diagnostics( state, current_time, dt, itime/dump_period_in_timesteps )  ! Write stat file
                    not_to_move_det_yet = .false. ;
                    call write_state( dump_no, state ) ! Now writing into the vtu files
                end if Conditional_Dump_RealTime
            end if
        end subroutine create_dump_vtu_and_checkpoints

        !This subroutine performs all the necessary steps to adapt the mesh and create new memory
        subroutine adapt_mesh_mp()
            !local variables
            type( scalar_field ), pointer ::  s_field, s_field2, s_field3
            type( vector_field ), pointer ::  U_x1, U_x2
            integer :: U_x1_stat, idim



            if (numberfields > 0) then ! If there is at least one instance of CVgalerkin then apply the method
                if (have_option('/mesh_adaptivity')) then ! Only need to use interpolation if mesh adaptivity switched on
                    call M2MInterpolation(state, packed_state, Mdims, CV_GIdims, CV_funs, Mspars%small_acv%fin, Mspars%small_acv%col , 0)
                endif
            endif

            do_reallocate_fields = .false.
            Conditional_Adaptivity_ReallocatingFields: if( have_option( '/mesh_adaptivity/hr_adaptivity') ) then
                if( have_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps') ) then
                    call get_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps', &
                        adapt_time_steps, default=5 )
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
                            itime, not_to_move_det_yet = .true. )
                        call run_diagnostics( state )
                        call adapt_state( state, metric_tensor, suppress_reference_warnings = .true.)
                        ! Copy U memory
                        do iphase=1,Mdims%nphase
                           U_x1=>extract_vector_field(state(iphase),"U",U_x1_stat)
                           U_x2=>extract_vector_field(state(iphase),"Velocity")
                           if(U_x1_stat==0)then
                              do idim=1,Mdims%ndim
                                 U_x1%val(idim,:)=U_x2%val(idim,:)
                              end do
                           end if
                        end do
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
                call deallocate_projection_matrices(CV_funs)
                call deallocate_projection_matrices(FE_funs)
if (is_flooding) then
    !sprint_to_do HACK (we should solve this properly): Reconstruct pressure from density because the pressure is uncorrectly interpolated between meshes
    s_field => extract_scalar_field( state(1), "Density" )
    s_field2 => extract_scalar_field( state(1), "Bathymetry" )
    s_field3 => extract_scalar_field( state(1), "Pressure" )
    if (size(s_field2%val)/=size(s_field%val)) then
        s_field3%val(:) = 9.81 * (s_field%val(:) + s_field2%val(1))
    else
        s_field3%val(:) = 9.81 * (s_field%val(:) + s_field2%val(:))
    end if
end if
                !!$ Compute primary scalars used in most of the code
                call Get_Primary_Scalars_new( state, Mdims )
                !Check if the user wants to store the outfluxes
                call initialize_multi_outfluxes(outfluxes)
                call pack_multistate(Mdims%npres,state,packed_state,&
                    multiphase_state,multicomponent_state)
                call prepare_absorptions(state, Mdims, multi_absorp)
                !Retrieve manning coefficient for flooding, this has to be called just after creating pack_multistate
                if (is_flooding) call get_FloodingProp(state, packed_state)
                call set_boundary_conditions_values(state, shift_time=.true.)
                !!$ Deallocating array variables:
                deallocate( &
                    !!$ Working arrays
                    suf_sig_diagten_bc, &
                    theta_gdiff, ScalarField_Source_Store, &
                    mass_ele, &
                    theta_flux, one_m_theta_flux, theta_flux_j, one_m_theta_flux_j, sum_theta_flux, &
                    sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j )
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
                if (is_porous_media) then
                    call get_RockFluidProp(state, packed_state)
                    call deallocate_porous_adv_coefs(upwnd)
                    call allocate_porous_adv_coefs(Mdims, upwnd)
                    !Clean the pipes memory if required
                    if (Mdims%npres > 1) call deallocate_multi_pipe_package(pipes_aux)
                end if

                call put_CSR_spars_into_packed_state()
                ! SECOND INTERPOLATION CALL - After adapting the mesh ******************************
                if (numberfields > 0) then
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
                allocate( theta_flux( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    one_m_theta_flux( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    theta_flux_j( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    one_m_theta_flux_j( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    sum_theta_flux( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    sum_one_m_theta_flux( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    sum_theta_flux_j( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    sum_one_m_theta_flux_j( Mdims%nphase, scvngi_theta*Mdims%cv_nloc*Mdims%totele * igot_theta_flux ), &
                    theta_gdiff( Mdims%nphase, Mdims%cv_nonods ), &
                    ScalarField_Source_Store( Mdims%nphase, Mdims%cv_nonods ) )
                sum_theta_flux = 1. ; sum_one_m_theta_flux = 0.
                sum_theta_flux_j = 1. ; sum_one_m_theta_flux_j = 0.
                ScalarField_Source_Store=0.
!                allocate(opt_vel_upwind_coefs_new(Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%mat_nonods)); opt_vel_upwind_coefs_new =0.
!                allocate(opt_vel_upwind_grad_new(Mdims%ndim, Mdims%ndim, Mdims%nphase, Mdims%mat_nonods)); opt_vel_upwind_grad_new =0.
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


    subroutine calc_components()
        implicit none

        perm_field => extract_tensor_field(packed_state,"Permeability")
        PhaseVolumeFractionComponentSource%val = 0.0
        velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
        saturation_field=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
        old_saturation_field=>extract_tensor_field(packed_state,"PackedOldPhaseVolumeFraction")

        PhaseVolumeFractionComponentSource%val = 0.0
        D_s  => extract_tensor_field( packed_state, "PackedDensity" )
        DC_s  => extract_tensor_field( packed_state, "PackedComponentDensity" )
        DCOLD_s  => extract_tensor_field( packed_state, "PackedOldComponentDensity" )
        MFC_s  => extract_tensor_field( packed_state, "PackedComponentMassFraction" )
        MFCOLD_s  => extract_tensor_field( packed_state, "PackedOldComponentMassFraction" )


        !!$ Starting loop over components
        Loop_Components: do icomp = 1, Mdims%ncomp

            tracer_field=>extract_tensor_field(multicomponent_state(icomp),"PackedComponentMassFraction")
            density_field=>extract_tensor_field(multicomponent_state(icomp),"PackedComponentDensity",stat)

            if( have_option( '/material_phase[' // int2str( Mdims%nstate - Mdims%ncomp ) // &
                ']/is_multiphase_component/KComp_Sigmoid' ) .and. Mdims%nphase > 1 ) then

                !!$ Computing the absorption term for the multi-components equation

                Component_Absorption => extract_tensor_field( multicomponent_state(icomp), "ComponentAbsorption")

                call Calculate_ComponentAbsorptionTerm( packed_state, icomp, ndgln%cv, &
                     Mdims, D_s%val, Porosity_field%val, mass_ele, Component_Absorption%val )

            end if
            !!$ NonLinear iteration for the components advection:
            Loop_NonLinearIteration_Components: do its2 = 1, NonLinearIteration_Components

                Mdisopt%comp_use_theta_flux = .false. ; Mdisopt%comp_get_theta_flux = .true.
                call INTENERGE_ASSEM_SOLVE( state, multicomponent_state(icomp), &
                    Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,upwnd,&
                    tracer_field,velocity_field,density_field, multi_absorp, dt, &
                    SUF_SIG_DIAGTEN_BC, Porosity_field%val, &
                    igot_t2, igot_theta_flux, &
                    Mdisopt%comp_get_theta_flux, Mdisopt%comp_use_theta_flux, &
                    theta_gdiff, eles_with_pipe, pipes_aux,&
                    thermal = .false.,& ! the false means that we don't add an extra source term
                    theta_flux=theta_flux, one_m_theta_flux=one_m_theta_flux, theta_flux_j=theta_flux_j, one_m_theta_flux_j=one_m_theta_flux_j,&
                    icomp=icomp, saturation=saturation_field, Permeability_tensor_field = perm_field)
                tracer_field%val = min (max( tracer_field%val, 0.0), 1.0)
            end do Loop_NonLinearIteration_Components

            sum_theta_flux = sum_theta_flux + theta_flux
            sum_one_m_theta_flux = sum_one_m_theta_flux + one_m_theta_flux
            sum_theta_flux_j = sum_theta_flux_j + theta_flux_j
            sum_one_m_theta_flux_j = sum_one_m_theta_flux_j + one_m_theta_flux_j

            ! We have divided through by density
            PhaseVolumeFractionComponentSource%val(1,:,:) = PhaseVolumeFractionComponentSource%val(1,:,:) + THETA_GDIFF

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

        MeanPoreCV=>extract_vector_field(packed_state,"MeanPoreCV")

        do icomp = 1, Mdims%ncomp

            if ( have_option( '/material_phase[' // int2str( Mdims%nstate - Mdims%ncomp ) // &
                ']/is_multiphase_component/KComp_Sigmoid' ) .and. Mdims%nphase > 1 ) then

                call Calculate_ComponentAbsorptionTerm( packed_state, icomp, ndgln%cv, &
                     Mdims, D_s%val, Porosity_field%val, mass_ele, Component_Absorption%val )

                do cv_nodi = 1, Mdims%cv_nonods
                    Loop_Phase_SourceTerm1: do iphase = 1, Mdims%nphase
                        Loop_Phase_SourceTerm2: do jphase = 1, Mdims%nphase
                            PhaseVolumeFractionComponentSource%val( 1, iphase, cv_nodi ) = &
                                PhaseVolumeFractionComponentSource%val( 1, iphase, cv_nodi ) &
                                - Component_Absorption%val( iphase, jphase, cv_nodi ) &
                                * MFC_s%val( icomp, jphase, cv_nodi ) &
                                / DC_s%val( icomp, iphase, cv_nodi )
                        end do Loop_Phase_SourceTerm2
                    end do Loop_Phase_SourceTerm1
                end do

            end if

            ! For compressibility
            do iphase = 1, Mdims%nphase
                do cv_nodi = 1, Mdims%cv_nonods
                    PhaseVolumeFractionComponentSource%val( 1, iphase, cv_nodi ) =  &
                        PhaseVolumeFractionComponentSource%val( 1, iphase, cv_nodi ) &
                        + MeanPoreCV%val( 1, cv_nodi ) * MFCOLD_s%val( icomp, iphase, cv_nodi ) &
                        * ( DCOLD_s%val( icomp, iphase, cv_nodi ) - DC_s%val( icomp, iphase, cv_nodi ) ) &
                        * old_saturation_field%val( 1, IPHASE, Mdims%cv_nonods ) &
                        / ( DC_s%val( ICOMP, IPHASE, CV_NODI ) * DT )
                end do
            end do

        end do ! icomp

        if ( is_porous_media .and. have_option( '/material_phase[' // int2str( Mdims%nstate - Mdims%ncomp ) // &
            ']/is_multiphase_component/Comp_Sum2One' ) .and. ( Mdims%ncomp > 1 ) ) then
            call Cal_Comp_Sum2One_Sou( packed_state, Mdims )
        end if


    end subroutine calc_components


    subroutine adapt_mesh_within_FPI(ExitNonLinearLoop, adapt_mesh_in_FPI, its, flag)
        implicit none
        logical, intent(inout) :: ExitNonLinearLoop, adapt_mesh_in_FPI
        integer, intent(inout) :: its
        integer, intent(in) :: flag
        !Local variables
        type(scalar_field), pointer :: sat1, sat2
        integer :: iphase
        integer, save :: phaseToAdapt = -1
        real, save :: Inf_tol = -1, non_linear_tol = -1


        if (phaseToAdapt<0) then
            !Retrieve which phase has the options to adapt to
            do iphase = 1, Mdims%nphase
                if (have_option('/material_phase['//int2str(iphase-1)//']/scalar_field::PhaseVolumeFraction/prognostic/adaptivity_options')) then
                    phaseToAdapt = iphase
                    exit
                end if
            end do
            !Make sure it is not the last phase! and that the option can be used
            if (phaseToAdapt==Mdims%nphase .or. phaseToAdapt < 0 .or. .not.is_porous_media ) then
                ewrite(0,*) "WARNING: Adapt_mesh_within_FPI requires porous media flow and the last PhaseVolumeFraction NOT to be the only target for adapting the mesh"
                return
            end if
            !Store the settings selected by the user
            call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration', non_linear_tol, default = 5e-2)
            call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Infinite_norm_tol',Inf_tol, default = 0.03)
        end if


        select case (flag)
            case (1)
                !Relax the convergence criteria since we don't need much precision at this stage
                !Since non_linear_tol is they key convergence criterion, we use a relative value and we reduce it half-order
                call set_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration', min(non_linear_tol*10., 1e-1))
                if (.not.have_option('/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Infinite_norm_tol'))then
                    !Create the option
                    call copy_option( '/timestepping/timestep/', &
                        '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Infinite_norm_tol')
                end if
                !10% tolerance provides a good enough result
                call set_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Infinite_norm_tol',min(Inf_tol*5.,0.1))
        case default

            !Four steps
            !1. Store OldPhaseVolumeFraction
            do iphase = 1, Mdims%n_in_pres
                sat1 => extract_scalar_field( state(iphase), "OldPhaseVolumeFraction" )
                sat2  => extract_scalar_field( state(iphase), "Saturation_bak" )
                sat2%val = sat1%val
            end do
            !2.Prognostic field to adapt the mesh to, has to be a convolution of old a new saturations
            sat2  => extract_scalar_field( state(phaseToAdapt), "OldPhaseVolumeFraction" )
            sat1  => extract_scalar_field( state(phaseToAdapt), "PhaseVolumeFraction" )
            !It is important that the average keep the sharpness of the interfaces
            sat1%val = abs(sat1%val - sat2%val)**0.8
            call adapt_mesh_mp()
            !3.Reconstruct the Saturation of the first phase
            sat1  => extract_scalar_field( state(phaseToAdapt), "PhaseVolumeFraction" )
            sat1%val = 1.0
            do iphase = 1, Mdims%n_in_pres
                if (iphase /= phaseToAdapt) then
                    sat2  => extract_scalar_field( state(iphase), "PhaseVolumeFraction" )
                    sat1%val = sat1%val - sat2%val
                end if
            end do
            !4. Copy back to OldPhaseVolumeFraction
            do iphase = 1, Mdims%n_in_pres
                sat1 => extract_scalar_field( state(iphase), "OldPhaseVolumeFraction" )
                sat2  => extract_scalar_field( state(iphase), "Saturation_bak" )
                sat1%val = sat2%val
            end do
            !Pointing to porosity again is required
            porosity_field=>extract_vector_field(packed_state,"Porosity")
            !Now we have to converge again within the same time-step
            ExitNonLinearLoop = .false.; its = 1
            adapt_mesh_in_FPI = .false.

            !Set the original convergence criteria since we are now solving the equations
            call set_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration', non_linear_tol)
            call set_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Infinite_norm_tol',Inf_tol)
            !Do not re-adapt the mesh
            t_adapt_threshold = acctim + 1.0 !Just to ensure that we do not re-mesh again
        end select

    end subroutine adapt_mesh_within_FPI


subroutine BadElementTest(Quality_list, flag)

    implicit none
    integer, intent(in) :: flag
    type(bad_elements), allocatable, dimension(:), intent(inout) :: Quality_list
    ! local variable
    integer :: i

    select case(flag)

            ! Allocate memory for the quality_table to check the angles of each element and print out the diagnostics for the mesh
        case (1)
            
			bad_element = have_option('/numerical_methods/Bad_element_fix/')
			if (have_option('/io/Mesh_Diagnostics_Angles')) then
				mesh_diagnostics = .true.
				shape = option_shape('/io/Mesh_Diagnostics_Angles')
				allocate(quality_table(shape(1)))
				call get_option( '/io/Mesh_Diagnostics_Angles', quality_table)
				allocate(diagnostics(shape(1)))
				diagnostics(:) = -1
			end if
           

        case default
            if (bad_element .or. mesh_diagnostics)  then
                ! Check bad elements (angles larger than specified in Max_bad_angle) at first time step or after mesh adapt
                if (first_time_step .or. after_adapt) then
                    if (allocated(Quality_list)) then
                        i=1
                        do while (allocated(Quality_list(i)%rotmatrix) )
                            deallocate(Quality_list(i)%rotmatrix)
                            i = i+1
                        end do
                        deallocate(Quality_list)
                    end if
                    allocate(Quality_list(Mdims%totele+1))
                    call get_var_from_packed_state(packed_state, Coordinate = X_ALL)
                    ! create table for stats on the angles of the mesh
                    if (mesh_diagnostics) then
                        do i=1, size(quality_table)
                            call CheckElementAngles(X_ALL, Mdims%totele, ndgln%x_p1, Mdims%x_nloc_p1 ,quality_table(i), Quality_list, bad_element, diagnostics(i))
                        end do



                        if (getprocno() == 1) then
                            if (after_adapt) then
                                ewrite(0, '( 38("-") / 1X A, I6, 10X, A / 38("-") )') "Time = ", itime ,"after adapt"
                            else
                                ewrite(0, '( 38("-") / 1X A, I10 / 38("-") )') "Time = ", itime
                            end if

                            ewrite(0, '(1X A, 4X A, 3X A, T9, "|", T24, "|" / 38("-"), T9, "|", T24, "|")') "Angle", "No. elements", "Percentage"
                            do i=1, size(quality_table)
                                ewrite(0, '(1X F6.2, 1X, I10, 6X, F7.2, T9, "|", T24 "|")') quality_table(i), diagnostics(i), diagnostics(i)*100./Mdims%totele
                            end do
                            ewrite(0, '(38("-"))')
                        end if
                    end if

                    if (.not. bad_element) then
                        deallocate(Quality_list) ! don't change element properties - throw away Quality list
                    else

                        if (Max_bad_angle < 0) call get_option( '/numerical_methods/Bad_element_fix/Angle', Max_bad_angle, default = 177. )
                        call CheckElementAngles(X_ALL, Mdims%totele, ndgln%x_p1, Mdims%x_nloc_p1, Max_bad_angle, Quality_list, bad_element)
                    endif
                end if
            end if

    end select

end subroutine BadElementTest

    end subroutine MultiFluids_SolveTimeLoop





end module multiphase_time_loop
