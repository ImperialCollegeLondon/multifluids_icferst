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

subroutine multiphase_prototype_wrapper() bind(C)

    use fldebug
    use elements
    use fields
    use state_module
    use populate_state_module
    use reserve_state_module
    use diagnostic_variables
    use diagnostic_fields_wrapper
    use write_state_module
    use timeloop_utilities
    use timers
    use parallel_tools
    use reference_counting
    use global_parameters
    use diagnostic_fields_new_multiphase, only : &
        & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
        & check_diagnostic_dependencies
    use field_priority_lists
    use spud
    use checkpoint
    use boundary_conditions_from_options
    use discrete_properties_module
    use multiphase_module
    use multimaterial_module
    use field_equations_cv, only: initialise_advection_convergence
    use memory_diagnostics
    use multiphase_time_loop
    !use multiphase_rheology
    use MeshDiagnostics
    use write_gmsh
    use signals

    !use mp_prototype
    use tictoc
    implicit none



    type(state_type), dimension(:), pointer :: state

    integer :: i, dump_no = 0
    integer :: ntsol, nonlinear_iterations

    character(len = option_path_len) :: filename
    character(len = option_path_len) :: simulation_name, dump_format

    real :: finish_time, nonlinear_iteration_tolerance

    ! Establish signal handlers
    call initialise_signals()

    call get_option("/simulation_name",filename)

    call set_simulation_start_times()
    call initialise_walltime
    timestep = 0

#ifdef HAVE_MEMORY_STATS
    ! this is to make sure the option /io/log_output/memory_diagnostics is read
    call reset_memory_logs()
#endif

    call initialise_write_state


    !!Retrieve what type of simulation are we doing
    call get_simulation_type()
    !Flag the first time step
    first_time_step = .true.
    ! Read state from .mpml file
    call populate_multi_state(state)
!    call populate_state(state)

    !If desired by the user create a bin msh file
    if (have_option("/geometry/create_binary_msh")) call create_bin_msh_file(state)

    ! Check the diagnostic field dependencies for circular dependencies
    call check_diagnostic_dependencies(state)

    !allocate(rheology(size(state)))
    !call initialize_rheologies(state,rheology)
    !call calculate_rheologies(state,rheology)

    !Prepare the generic warning message for users, this message is intended to help to debug a bad input file/mesh
    call set_up_generic_warning_message()!print *, trim(multi_generic_warning)


    !--------------------------------------------------------------------------------------------------------
    ! This should be read by the Copy_Outof_Into_State subrts

    ! set the remaining timestepping options, needs to be before any diagnostics are calculated
    call get_option("/timestepping/timestep", dt)
    !  if(have_option("/timestepping/adaptive_timestep/at_first_timestep")) then
    !    call calc_cflnumber_field_based_dt(state, dt, force_calculation = .true.)
    !    call set_option("/timestepping/timestep", dt)
    !  end if

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)

    call get_option('/simulation_name',simulation_name)
    call initialise_diagnostics(trim(simulation_name),state)

    ! Calculate the number of scalar fields to solve for and their correct
    ! solve order taking into account dependencies.
    call get_ntsol(ntsol)

    call initialise_field_lists_from_options(state, ntsol)


    ! For multiphase simulations, we have to call calculate_diagnostic_phase_volume_fraction *before*
    ! copy_to_stored(state,"Old") is called below. Otherwise, OldPhaseVolumeFraction (in the phase
    ! containing the diagnostic PhaseVolumeFraction) will be zero and
    ! NonlinearPhaseVolumeFraction will be calculated incorrectly at t=0.
    if(option_count("/material_phase/vector_field::Velocity/prognostic") > 1) then
        call calculate_diagnostic_phase_volume_fraction(state)
    end if

    ! set the nonlinear timestepping options, needs to be before the adapt at first timestep
    call get_option('/timestepping/nonlinear_iterations',nonlinear_iterations,&
        & default=1)
    call get_option("/timestepping/nonlinear_iterations/tolerance", &
        & nonlinear_iteration_tolerance, default=0.0)

    ! Auxilliary fields.
    call allocate_and_insert_auxilliary_fields(state)
    call copy_to_stored_values(state,"Old")
    call copy_to_stored_values(state,"Iterated")
    call relax_to_nonlinear(state)

    call enforce_discrete_properties(state)

    call run_diagnostics(state)

    !     Determine the output format 
    dump_format = "vtk"

    ! initialise the multimaterial fields
    call initialise_diagnostic_material_properties(state)

    call calculate_diagnostic_variables(state)
    call calculate_diagnostic_variables_new(state)

    call tictoc_reset()
    call tic(TICTOC_ID_SIMULATION)

    !--------------------------------------------------------------------------------------------------------

    call initialise_convergence(filename, state)
    call initialise_steady_state(filename, state)
    call initialise_advection_convergence(state)


    ! this may already have been done in populate_state, but now
    ! we evaluate at the correct "shifted" time level:
    call set_boundary_conditions_values(state, shift_time=.true.)

    call enforce_discrete_properties(state, only_prescribed=.true., &
        exclude_interpolated=.true., &
        exclude_nonreprescribed=.true.)

    call print_tagged_references(0)

    ! Call the multiphase_prototype code
    !call multiphase_prototype(state, dt, &
    !                          nonlinear_iterations, nonlinear_iteration_tolerance, &
    !                          dump_no)
    call MultiFluids_SolveTimeLoop( state, &
        dt, nonlinear_iterations, dump_no )




    call close_diagnostic_files()

    ! Deallocate state
    do i = 1, size(state)
        call deallocate(state(i))
    end do
    deallocate(state)

    call deallocate_reserve_state()

    ! Clean up registered diagnostics
    call destroy_registered_diagnostics()

    ! Delete the transform_elements cache.
    call deallocate_transform_cache()

    ewrite(2, *) "Tagged references remaining:"
    call print_tagged_references(0)

#ifdef HAVE_MEMORY_STATS
    call print_current_memory_stats(0)
#endif

    call toc(TICTOC_ID_SIMULATION)
    call tictoc_report(2, TICTOC_ID_SIMULATION)

contains

    subroutine set_simulation_start_times()
        !!< Set the simulation start times

        call get_option("/timestepping/current_time", simulation_start_time)

        call cpu_time(simulation_start_cpu_time)
        call allmax(simulation_start_cpu_time)

        simulation_start_wall_time = wall_time()
        call allmax(simulation_start_wall_time)

    end subroutine set_simulation_start_times


    subroutine populate_multi_state(state)
        implicit none
        type(state_type), dimension(:), pointer, intent(inout) :: state
        !Local variables
        integer :: i, nphase, npres
        character( len = option_path_len ) :: option_path
        integer :: stat

        nphase = option_count("/material_phase")
        npres = option_count("/material_phase/scalar_field::Pressure/prognostic")
        !Adjust nphase to not account for the extra phases added by the wells
        nphase = nphase/npres

        !If the input file is the simplified version of the mpml
        ! file then we need to reconstruct the original mpml file
        ! the flag /is_porous_media only exists in the IC_FERST schema and triggers this
        if (have_option("/is_porous_media")) call convert_FERST_input_file(nphase, npres)

        !Prepare some specific modifications prior to populating state
        !If the extra mesh have not been created, create them here
        if (.not.is_P0DGP1CV) then!We don't need this field for P0DGP1
            if (.not. have_option("/geometry/mesh::VelocityMesh_Continuous")) then
                call copy_option("/geometry/mesh::VelocityMesh", "/geometry/mesh::VelocityMesh_Continuous")
                call set_option("/geometry/mesh::VelocityMesh_Continuous/from_mesh/mesh_continuity", "continuous")
            end if
        end if
        if (.not. have_option("/geometry/mesh::PressureMesh_Continuous")) then
            call copy_option("/geometry/mesh::PressureMesh", "/geometry/mesh::PressureMesh_Continuous")
            call set_option("/geometry/mesh::PressureMesh_Continuous/from_mesh/mesh_continuity", "continuous")
        end if

        if (.not. have_option("/geometry/mesh::PressureMesh_Discontinuous")) then
            call copy_option("/geometry/mesh::PressureMesh", "/geometry/mesh::PressureMesh_Discontinuous")
            call set_option("/geometry/mesh::PressureMesh_Discontinuous/from_mesh/mesh_continuity", "discontinuous")
        end if

        if (.not. have_option("/geometry/mesh::P0DG")) then
            call copy_option("/geometry/mesh::PressureMesh", "/geometry/mesh::P0DG")
            call set_option("/geometry/mesh::P0DG/from_mesh/mesh_shape/polynomial_degree", 0)
            call set_option("/geometry/mesh::P0DG/from_mesh/mesh_continuity", "discontinuous")
        end if


        if (have_option('/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI')) then
            !Create necessary backups for storing the saturation (in a way that it is also adapted)
            do i = 1, nphase
                option_path = "/material_phase["// int2str( i - 1 )//"]/scalar_field::Saturation_bak"
                call copy_option("/material_phase["// int2str( i - 1 )//"]/scalar_field::PhaseVolumeFraction",&
                     trim(option_path))
                !Make sure the field is not shown
                if (.not.have_option(trim(option_path)//"/prognostic/output/exclude_from_vtu")) then
                    !Copy an option that lways exists to ensure we exclude the new field from vtu
                    call copy_option("/simulation_type",trim(option_path)//"/prognostic/output/exclude_from_vtu")
                end if
                !Make sure that this field is not the objective of adaptivity
                if (have_option(trim(option_path)//"/prognostic/adaptivity_options")) then
                    call delete_option(trim(option_path)//"/prognostic/adaptivity_options")
                end if
            end do
        end if


        if (have_option('/physical_parameters/black-oil_PVT_table')) then

        !Maybe no need for a field in state? and just calculate using required conditions????

            !Create necessary memory to store the mass fraction of one of the pseudo-components (in a way that it is also adapted)
            if (nphase /= 3)then
                FLAbort('Black-Oil modelling requires three phases. Phase 1 Aqua, phase 2 liquid, phase 3 vapour')
            end if

            option_path = "/material_phase["// int2str( nphase -1 )//"]/scalar_field::VapourMassFraction"
            call copy_option("/material_phase["// int2str( nphase - 1 )//"]/scalar_field::PhaseVolumeFraction",&
                 trim(option_path))
            !Make sure the field is not shown
            if (.not.have_option(trim(option_path)//"/prognostic/output/exclude_from_vtu")) then
                !Don't know how to set exclude_from_vtu to true from the spud options, hence,
                !since Porous_media HAS to be true I copy it to obtain the same effect
                call copy_option("/simulation_type/porous_media",&
                 trim(option_path)//"/prognostic/output/exclude_from_vtu")
            end if
            !Make sure that this field is not the objective of adaptivity
            if (have_option(trim(option_path)//"/prognostic/adaptivity_options")) then
                call delete_option(trim(option_path)//"/prognostic/adaptivity_options")
            end if
        end if


        if (is_porous_media .and. have_option('/mesh_adaptivity/hr_adaptivity')) then
            !Ensure that preserve_mesh_regions is on, since otherwise it does not work
            !Don't know how to set exclude_from_vtu to true from the spud options, hence,
            !since Porous_media HAS to be true I copy it to obtain the same effect
            if (.not.have_option("/mesh_adaptivity/hr_adaptivity/preserve_mesh_regions")) then
                call copy_option("/simulation_type/porous_media",&
                 "/mesh_adaptivity/hr_adaptivity/preserve_mesh_regions")
            end if
        end if

        if (is_porous_media .and. (have_option('/io/output_darcy_vel'))) then
            !Create a copy of the velocity fields to store the DarcyVelocity in it
            do i = 1, nphase
                option_path = "/material_phase["// int2str( i - 1 )//"]/vector_field::"
                call copy_option(trim(option_path)//"Velocity", trim(option_path)//"DarcyVelocity")
                if (have_option(trim(option_path)//"DarcyVelocity"//"/prognostic/tensor_field::Viscosity")) &
                        call delete_option(trim(option_path)//"DarcyVelocity"//"/prognostic/tensor_field::Viscosity")

                if (have_option(trim(option_path)//"DarcyVelocity"//"/prognostic/vector_field::Absorption"))&
                        call delete_option(trim(option_path)//"DarcyVelocity"//"/prognostic/vector_field::Absorption")

                if (have_option(trim(option_path)//"DarcyVelocity"//"/prognostic/adaptivity_options"))&
                        call delete_option(trim(option_path)//"DarcyVelocity"//"/prognostic/adaptivity_options")
            end do
        end if

        !Call fluidity to populate state
        call populate_state(state)

    end subroutine populate_multi_state


    subroutine convert_FERST_input_file(nphase, npres)
                ! ################################
        ! In this subroutine. We populate the input file elements that are not user-selected in IC_FERST.rnc. (i.e. we set all the things that the user doesn't need to worry about).
        !IT IS VERY IMPORTANT NOT TO CHANGE THE ORDERING IN THIS SUBROUTINE!!!
        implicit none
        integer, intent(in) :: nphase, npres
        !Local variables
        integer :: stat, i, simulation_quality = 10
        real :: theta, aux
        character( len = option_path_len ) :: option_path, quality_option

        !Add simulation type like in the normal fashion
        call add_option("/simulation_type/porous_media", stat=stat)

        ! GEOMETRY OPTIONS
        !Add quadrature option (always equals to 5 I think we need this for legacy reasons)
        option_path = "/geometry/quadrature/degree"
        call add_option(trim(option_path), stat=stat)
        call set_option(trim(option_path), 5)

        !Check quality option to decide mesh type, theta and advection schemes
        option_path = "/geometry/simulation_quality"
        call get_option(trim(option_path), quality_option, stat=stat)
        if (trim(quality_option) == "fast") then
            simulation_quality = 1
        else if (trim(quality_option) == "precision") then
            simulation_quality = 100
        else if (trim(quality_option) == "discontinuous_pressure") then
            simulation_quality = 1000
        else !balanced, the recommended one
            simulation_quality = 10
        end if

        !Unless it is discontinuous_pressure (and this will change in the future) specify the CV_pressure formulation
        if (simulation_quality < 1000) then
            option_path = "/material_phase[0]/scalar_field::Pressure/prognostic/CV_P_matrix/"
            call add_option(trim(option_path), stat=stat)
        end if

        option_path = "/geometry/mesh::VelocityMesh/"
        call add_option(trim(option_path)//"from_mesh", stat=stat)
        call add_option(trim(option_path)//"from_mesh/mesh::CoordinateMesh", stat=stat)
        call add_option(trim(option_path)//"from_mesh/mesh_continuity", stat=stat)
        call set_option(trim(option_path)//"from_mesh/mesh_continuity", "discontinuous")
        call add_option(trim(option_path)//"from_mesh/mesh_shape/element_type", stat=stat)
        call set_option(trim(option_path)//"from_mesh/mesh_shape/element_type", "lagrangian")
        call add_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", stat=stat)
        if (simulation_quality < 100) then
            call set_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", 0)
        else if (simulation_quality < 1000) then
            call set_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", 1)
        else
            call set_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", 2)
        end if
        call add_option(trim(option_path)//"from_mesh/stat/exclude_from_stat", stat=stat)

        option_path = "/geometry/mesh::PressureMesh/"
        call add_option(trim(option_path)//"from_mesh", stat=stat)
        call add_option(trim(option_path)//"from_mesh/mesh::CoordinateMesh", stat=stat)
        call add_option(trim(option_path)//"from_mesh/mesh_continuity", stat=stat)
        if (simulation_quality >= 1000) then
            call set_option(trim(option_path)//"from_mesh/mesh_continuity", "discontinuous")
        else
            call set_option(trim(option_path)//"from_mesh/mesh_continuity", "continuous")
        end if
        call add_option(trim(option_path)//"from_mesh/mesh_shape/element_type", stat=stat)
        call set_option(trim(option_path)//"from_mesh/mesh_shape/element_type", "lagrangian")
        call add_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", stat=stat)
        if (simulation_quality >= 100 .and. simulation_quality < 1000) then
            call set_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", 2)
        else
            call set_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", 1)
        end if
        call add_option(trim(option_path)//"from_mesh/stat/exclude_from_stat", stat=stat)

        !Select the theta
        theta = -1.
        if (simulation_quality < 10) theta = 1.0
        !If it is balanced quality select local decision on theta
        if (simulation_quality >= 10 .and. simulation_quality < 100) call add_option("/numerical_methods/local_upwinding", stat=stat)

            ! IO STAT OPTIONS
        option_path = "/io/output_mesh[0]/name"
        call add_option(trim(option_path), stat=stat)
        call set_option(trim(option_path),"PressureMesh")

        !Set up the non-linear solver to use an automatic backtracking method
        option_path = "/timestepping/nonlinear_iterations/Fixed_Point_iteration/Bactracking_factor"
        call add_option(trim(option_path), stat=stat)
        call set_option(trim(option_path),-10.)


        do i = 1, nphase

                    ! NEW TREATMENT OF EQUATION OF STATE AND SPECIFYING THE DENSITY OF EACH PHASE
            option_path = "/material_phase["// int2str( i - 1)//"]/equation_of_state/incompressible/linear/all_equal"
            if (.not. have_option(trim(option_path))) then
                call copy_option("/material_phase["// int2str( i - 1 )//"]/Density", trim(option_path)) ! Note this is an entirely new density option.
            end if

            ! REMOVE THE DENSITY FIELD ENIRELY - IT'S CONFUSIONG - CREATE IT INTERNALLY AS A COPY OF THE SATURATIONF FIELD MINUS THE SOLVER (check how safe this step is - otherwise do it the long way)
            option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::Density"
            if (.not. have_option(trim(option_path))) then
                call copy_option("/material_phase["// int2str( i - 1 )//"]/scalar_field::PhaseVolumeFraction", trim(option_path))
                call delete_option(trim(option_path)//"/prognostic/solver")
            end if


            ! SCALAR_FIELD(PRESSURE) OPTIONS ADDED AUTOMATICALLY
            option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::Pressure/prognostic"
            if (have_option(trim(option_path))) then
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::Pressure/prognostic/spatial_discretisation/continuous_galerkin"
                call add_option(trim(option_path), stat=stat)

                            ! Stat, convergence, detectors, steady state settings
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::Pressure/prognostic/stat"
                call add_option(trim(option_path), stat=stat)
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::Pressure/prognostic/convergence/include_in_convergence"
                call add_option(trim(option_path), stat=stat)
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::Pressure/prognostic/detectors/exclude_from_detectors"
                call add_option(trim(option_path), stat=stat)
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::Pressure/prognostic/steady_state/include_in_steady_state"
                call add_option(trim(option_path), stat=stat)

                ! Solver options - Iterative Method
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::Pressure/prognostic/solver/iterative_method::gmres/restart"
                call add_option(trim(option_path), stat = stat)
                call set_option(trim(option_path), 0)
                !Preconditioner
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::Pressure/prognostic/solver/preconditioner::hypre/hypre_type[0]/name"
                call add_option(trim(option_path), stat = stat)
                call set_option(trim(option_path), "boomeramg")
            end if

                ! VECTOR_FIELD(VELOCITY) OPTIONS ADDED AUTOMATICALLY
            option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic"
            if (have_option(trim(option_path))) then

                !We don't need to repopulate much of this as the default values are good for porous media
                option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin"
                call add_option(trim(option_path), stat=stat)

                !Set this values just in case. We don't need them for porous media
                option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/temporal_discretisation/theta"
                call add_option(trim(option_path), stat=stat)
                call set_option(trim(option_path), theta)
                option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/temporal_discretisation/relaxation"
                call add_option(trim(option_path), stat=stat)
                call set_option(trim(option_path), 1.0)

                !Create viscosity from simplified version and create memory for the absoprtion
                option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/tensor_field::Viscosity/prescribed"
                call add_option(trim(option_path)//"/mesh::PressureMesh", stat=stat)
                call add_option(trim(option_path)//"/value::WholeMesh/isotropic/constant", stat=stat)!<= set up mesh
                call get_option("/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/viscosity",  aux)
                call set_option(trim(option_path)//"/value::WholeMesh/isotropic/constant", aux)!<= introduce value from the user

                !Create memory for the absoprtion
                option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/vector_field::Absorption/diagnostic"
                call add_option(trim(option_path)//"/mesh::PressureMesh_Discontinuous", stat=stat)!<= set up mesh
                call add_option(trim(option_path)//"/algorithm::Internal", stat=stat)!<= set up algorithm type
                    ! Stat, convergence, detectors, steady state settings
                    call add_option(trim(option_path)//"/stat", stat=stat)
                    call add_option(trim(option_path)//"/convergence/include_in_convergence", stat=stat)
                    call add_option(trim(option_path)//"/detectors/exclude_from_detectors", stat=stat)
                    call add_option(trim(option_path)//"/steady_state/include_in_steady_state", stat=stat)

                ! Stat, convergence, detectors, steady state settings
                option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/stat"
                call add_option(trim(option_path), stat=stat)
                option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/convergence/include_in_convergence"
                call add_option(trim(option_path), stat=stat)
                option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/detectors/exclude_from_detectors"
                call add_option(trim(option_path), stat=stat)
                option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/steady_state/include_in_steady_state"
                call add_option(trim(option_path), stat=stat)

                ! Solver options - Iterative Method (we don't actually solve for velocity, so we copy the options from the pressure just in case)
                call copy_option("/material_phase[0]/scalar_field::Pressure/prognostic/solver",trim(option_path)//"/solver" )
            end if


                ! SCALAR_FIELD(PHASE VOLUME FRACTION) OPTIONS ADDED AUTOMATICALLY
            option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::PhaseVolumeFraction/prognostic"
            if (have_option(trim(option_path))) then

                            ! Stat, convergence, detectors, steady state settings
                !We have removd this from scalar_field however we are only populating with this PhaseVolumeFraction
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::PhaseVolumeFraction/prognostic/stat"
                call add_option(trim(option_path), stat=stat)
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::PhaseVolumeFraction/prognostic/convergence/include_in_convergence"
                call add_option(trim(option_path), stat=stat)
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::PhaseVolumeFraction/prognostic/detectors/exclude_from_detectors"
                call add_option(trim(option_path), stat=stat)
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::PhaseVolumeFraction/prognostic/steady_state/include_in_steady_state"
                call add_option(trim(option_path), stat=stat)

                !Removed the selection of the equation as it was not even read by the code

                !Spatial discretisation
                if (simulation_quality < 10) then
                    option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::PhaseVolumeFraction/prognostic/spatial_discretisation/control_volumes/face_value"
                    call add_option(trim(option_path)//"::FirstOrderUpwind", stat=stat)
                else
                    option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::PhaseVolumeFraction/prognostic/spatial_discretisation/control_volumes"
                    option_path = trim(option_path) //"/face_value::FiniteElement/limit_face_value/limiter"
                    call add_option(trim(option_path)//"::Sweby", stat=stat)
                end if
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::PhaseVolumeFraction/prognostic/spatial_discretisation/conservative_advection"
                call add_option(trim(option_path), stat=stat)
                call set_option(trim(option_path), 1.)
                !Temporal discretisation
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::PhaseVolumeFraction/prognostic/temporal_discretisation/theta"
                call add_option(trim(option_path), stat=stat)
                call set_option(trim(option_path), theta)


                ! Solver options - Iterative Method
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::PhaseVolumeFraction/prognostic/solver/iterative_method::gmres/restart"
                call add_option(trim(option_path), stat = stat)
                call set_option(trim(option_path), 0)
                !Preconditioner
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::PhaseVolumeFraction/prognostic/solver/preconditioner::hypre/hypre_type[0]/name"
                call add_option(trim(option_path), stat = stat)
                call set_option(trim(option_path), "boomeramg")
            end if

        end do

        !####################################

    end subroutine convert_FERST_input_file



    subroutine set_up_generic_warning_message()
        implicit none

        character(len=1) :: aux
        multi_generic_warning = "WARNING: The set up of the numerical experiment may have some issues and a physical solution cannot be guaranteed."
        multi_generic_warning = trim(multi_generic_warning)//NEW_LINE('aux')
        multi_generic_warning = trim(multi_generic_warning)//NEW_LINE('aux')//"ONLY IF A PROBLEM APPEARS IN THE NUMERICAL SOLUTION, follow the next steps to try to find/solve the problem:"
        multi_generic_warning = trim(multi_generic_warning)//NEW_LINE('aux')
        multi_generic_warning = trim(multi_generic_warning)//NEW_LINE('aux')//"1. Reduce the time-step/use a fix time-step method"
        multi_generic_warning = trim(multi_generic_warning)//NEW_LINE('aux')//"2. If using an adaptive mesh => use a fix mesh"
        multi_generic_warning = trim(multi_generic_warning)//NEW_LINE('aux')//"3. Does the problem occur with a different mesh/element type? => Change the mesh (improve the quality, no sharp and obtuse elements) and/or element type"
        multi_generic_warning = trim(multi_generic_warning)//NEW_LINE('aux')//"4. Check the convergence criteria of the non-linear/linear solvers => make them tighter"
        multi_generic_warning = trim(multi_generic_warning)//NEW_LINE('aux')//"5. Simplify the physics of the problem until it works => remove physical effects and/or phases"
        multi_generic_warning = trim(multi_generic_warning)//NEW_LINE('aux')
        multi_generic_warning = trim(multi_generic_warning)//NEW_LINE('aux')//"May a bug be found, please file it in the corresponding repository following the standard procedure"

    end subroutine set_up_generic_warning_message

    subroutine create_bin_msh_file(state)
        implicit none
        type(state_type), dimension(:), pointer, intent(inout) :: state
        !Local variables
        type(vector_field), pointer :: output_positions
        character(len= OPTION_PATH_LEN ) :: msh_file

        if (IsParallel()) then
            ewrite(1, *) "In parallel we do not perform the conversion of the msh file."
            return
        end if

        ewrite(0, *) "Converting the input msh file into binary format..."
        call get_option('/geometry/mesh::CoordinateMesh/from_file/file_name',msh_file)
        output_positions => extract_vector_field(state(1), "Coordinate")
        call write_gmsh_file(trim(msh_file), output_positions)
        ewrite(0, *) "...Conversion done."

    end subroutine create_bin_msh_file

    subroutine get_simulation_type()
        !This subroutine selects the type of simulator to perform
        !and activates the flags from global_parameters accordingly
        implicit none
        integer :: Vdegree, Pdegree
        character( len = option_path_len ) :: quality_option
        !By default it is intertia dominated
        is_porous_media = have_option('/simulation_type/porous_media') .or. have_option('/is_porous_media')
        is_magma = have_option('/simulation_type/magma')
        is_flooding = have_option('/simulation_type/flooding')
        !Flag to set up the coupling with femdem
        is_multifracture = have_option( '/simulation_type/femdem_fracture' )
        !Flag to set up boiling
        is_boiling = have_option( '/simulation_type/boiling' )
        !Flag to set up blasting
        is_blasting = have_option( '/simulation_type/blasting' )
        !Has temperature
        has_temperature = have_option( '/material_phase[0]/scalar_field::Temperature/' )
        !Check if it is P0DGP1
        if (.not. have_option("/is_porous_media")) then
            call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/polynomial_degree', &
                Vdegree )
            call get_option( '/geometry/mesh::PressureMesh/from_mesh/mesh_shape/polynomial_degree', &
                Pdegree )
            is_P0DGP1CV = (Vdegree == 0) .and. (Pdegree == 1) .and. &
                    have_option( '/material_phase[0]/scalar_field::Pressure/prognostic/CV_P_matrix' )

            if ((Vdegree == 0) .and. (Pdegree == 1) .and. .not. is_P0DGP1CV) then
                ewrite(0, *) "P0DGP1 only works for porous media and the DCVFEM formulation."
                stop
            end if
        else
            is_P0DGP1CV = .false.
            !If using IC_FERST schema then this depends on the quality selected
            call get_option("/geometry/simulation_quality", quality_option)
            if ((trim(quality_option) == "fast") .or. (trim(quality_option) == "balanced")) is_P0DGP1CV = .true.
        end if
    end subroutine get_simulation_type

end subroutine multiphase_prototype_wrapper
