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

    !     Determine the output format.
    call get_option('/io/dump_format', dump_format)
    if(trim(dump_format) /= "vtk") then
        ewrite(-1,*) "You must specify a dump format and it must be vtk."
        FLExit("Rejig your FLML: /io/dump_format")
    end if

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
        nphase = option_count("/material_phase")
        npres = option_count("/material_phase/scalar_field::Pressure/prognostic")
        !Adjust nphase to not account for the extra phases added by the wells
        nphase = nphase/npres

        !Prepare some specific modifications prior to populating state
        !If the extra mesh have not been created, create them here
        if (.not. have_option("/geometry/mesh::VelocityMesh_Continuous")) then
            call copy_option("/geometry/mesh::VelocityMesh", "/geometry/mesh::VelocityMesh_Continuous")
            call set_option("/geometry/mesh::VelocityMesh_Continuous/from_mesh/mesh_continuity", "continuous")
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
                call copy_option("/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/Porous_media",&
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
                call copy_option("/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/Porous_media",&
                 "/mesh_adaptivity/hr_adaptivity/preserve_mesh_regions")
            end if
        end if

        !Call fluidity to populate state
        call populate_state(state)

    end subroutine populate_multi_state

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

    subroutine get_simulation_type()
        !This subroutine selects the type of simulator to perform
        !and activates the flags from global_parameters accordingly
        implicit none
        !By default it is intertia dominated
        is_porous_media = have_option('/simulation_type/porous_media')
        is_magma = have_option('/simulation_type/magma')
        is_flooding = have_option('/simulation_type/flooding')
        !Flag to set up the coupling with femdem
        is_multifracture = have_option( '/simulation_type/femdem_fracture' )
        !Flag to set up boiling
        is_boiling = have_option( '/simulation_type/boiling' )
        !Flag to set up blasting
        is_blasting = have_option( '/simulation_type/blasting' )

    end subroutine get_simulation_type

end subroutine multiphase_prototype_wrapper
