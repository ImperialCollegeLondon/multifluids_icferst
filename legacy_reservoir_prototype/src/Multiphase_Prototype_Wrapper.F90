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
    use futils
    use transform_elements
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
    use diagnostic_fields_new, only : &
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
    use multi_tools
    !use mp_prototype
    use tictoc
    implicit none

    !Local variables
    type(state_type), dimension(:), pointer :: state

    integer :: i, dump_no = 0
    integer :: ntsol, nonlinear_iterations

    character(len = option_path_len) :: filename, input_mpml
    character(len = option_path_len) :: simulation_name, dump_format

    real :: finish_time, nonlinear_iteration_tolerance, auxR, dump_period

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

    if ( have_option('/io/dump_period') ) then
      !Check, if we are not adapting the timestep that the dump_period is a multiple of the timestep
      if (.not.have_option( '/timestepping/adaptive_timestep' ) .and. &
      .not. have_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/adaptive_timestep_nonlinear')) then
        call get_option('/io/dump_period/constant', dump_period, default = 0.01)
        auxR = mod(dump_period,dt)
        if ( abs(auxR) > 1e-8 .or. dt > dump_period) then
          dt = dump_period/ceiling(dump_period/dt)
          ewrite(0, *) "WARNING: Dump period has to be a multiple of the time-step. Time-step adjusted to: ", dt
          call set_option("/timestepping/timestep", dt)
        end if
      end if
    end if
    !  if(have_option("/timestepping/adaptive_timestep/at_first_timestep")) then
    !    call calc_cflnumber_field_based_dt(state, dt, force_calculation = .true.)
    !    call set_option("/timestepping/timestep", dt)
    !  end if

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)

    call get_option('/simulation_name',simulation_name)
    call initialise_diagnostics(trim(simulation_name),state, ICFERST = .true.)

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
    call get_option('/solver_options/Non_Linear_Solver',nonlinear_iterations,&
        & default=1)
    call get_option("/solver_options/Non_Linear_Solver/tolerance", &
        & nonlinear_iteration_tolerance, default=0.0)
    ! Auxilliary fields.
    call allocate_and_insert_auxilliary_fields(state)
    call copy_to_stored_values(state,"Old")
    call copy_to_stored_values(state,"Iterated")
    call relax_to_nonlinear(state)

    call enforce_discrete_properties(state)

    call run_diagnostics(state)

    ! Determine the output format
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
    call print_tagged_references(1)

#ifdef HAVE_MEMORY_STATS
    call print_current_memory_stats(1)
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
        integer :: i, nphase, npres, aux, ncomp
        integer, dimension(2) :: shape
        character( len = option_path_len ) :: option_path, option_name, option_path_BAK
        integer :: stat, Pdegree, Vdegree, fields, k
        type (scalar_field), target :: targ_Store
        type (scalar_field), pointer :: sfield1, sfield2
        type (vector_field), pointer :: position
        integer, dimension(:), allocatable :: well_ids

        nphase = option_count("/material_phase")
        npres = option_count("/material_phase/scalar_field::Pressure/prognostic")
        ncomp = option_count("/material_phase/is_multiphase_component")
        !Adjust nphase to not account for the extra phases added by the wells
        nphase = nphase/npres - ncomp
        ncomp = ncomp/nphase
        !Read at the very beginning the property input file
        call read_fluid_and_rock_properties_from_csv()

        ! Fill in to old style schema internally
        call autocomplete_input_file(nphase, npres, ncomp)

        !Check if it is P0DGP1; do it after adapting new schema file
        call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/polynomial_degree', &
            Vdegree )
        call get_option( '/geometry/mesh::PressureMesh/from_mesh/mesh_shape/polynomial_degree', &
            Pdegree )
        is_P0DGP1CV = (Vdegree == 0) .and. (Pdegree == 1) .and. &
                .not. have_option( '/geometry/Advance_options/FE_Pressure' )

        if ((Vdegree == 0) .and. (Pdegree == 1) .and.( .not. is_P0DGP1CV &
                        .or. have_option('/inertia_dominated'))) then
            ewrite(0, *) "P0DGP1 does not work for inertia dominated simulations. If using the DCVFEM method use either one of the following options: "
            ewrite(0, *) "A. Use the P1DGP2CV formulation."
            ewrite(0, *) "B. Use the P1DGP1CV formulation with mass lumping = 100 in: /numerical_methods/lump_mass_matrix/lump_weight"
            stop
        end if

        !Prepare some specific modifications prior to populating state
        !If the extra mesh have not been created, create them here
        if (.not.is_P0DGP1CV) then!We don't need this field for P0DGP1
            if (.not. have_option("/geometry/mesh::VelocityMesh_Continuous")) then
                call copy_option("/geometry/mesh::VelocityMesh", "/geometry/mesh::VelocityMesh_Continuous")
                call set_option("/geometry/mesh::VelocityMesh_Continuous/from_mesh/mesh_continuity", "continuous")
            end if
        end if
        ewrite(1, *) "Create internally: PressureMesh_Continuous, PressureMesh_Discontinuous and P0DG mesh. Check multiphase_prototype_wrapper"
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

        !SPRINT_TO_DO generilise this to other fields
        ! if (have_option('/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI')) then
        !     ewrite(1, *) "For adapt within FPI, create necessary backups for storing the saturation. Check multiphase_prototype_wrapper"
        !     !Create necessary backups for storing the saturation (in a way that it is also adapted)
        !     do i = 1, nphase
        !         option_path = "/material_phase["// int2str( i - 1 )//"]/scalar_field::Saturation_bak"
        !         call copy_option("/material_phase["// int2str( i - 1 )//"]/scalar_field::PhaseVolumeFraction",&
        !              trim(option_path))
        !         !Make sure the field is not shown
        !         if (.not.have_option(trim(option_path)//"/prognostic/output/exclude_from_vtu")) then
        !             !Copy an option that always exists to ensure we exclude the new field from vtu
        !             call copy_option("/simulation_name",trim(option_path)//"/prognostic/output/exclude_from_vtu")
        !         end if
        !         !Make sure that this field is not the objective of adaptivity
        !         if (have_option(trim(option_path)//"/prognostic/adaptivity_options")) then
        !             call delete_option(trim(option_path)//"/prognostic/adaptivity_options")
        !         end if
        !     end do
        ! end if

        if (have_option('/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI')) then
            ewrite(1, *) "For adapt within FPI, create necessary backups for storing the old fields. Check multiphase_prototype_wrapper"
            !Create necessary backups for storing the saturation (in a way that it is also adapted)
            do i = 1, nphase * npres
                !First scalar prognostic fields
                option_path = "/material_phase["// int2str( i - 1 )//"]/scalar_field"
                fields = option_count("/material_phase["// int2str( i - 1 )//"]/scalar_field")
                do k = 1, fields
                  call get_option("/material_phase["// int2str( i - 1 )//"]/scalar_field["// int2str( k - 1 )//"]/name",option_name)
                  option_path = "/material_phase["// int2str( i - 1 )//"]/scalar_field::"//trim(option_name)
                  if (option_name(1:4)=="BAK_") cycle
                  if (trim(option_name)=="Density" .or. .not.have_option(trim(option_path)//"/prognostic" )) cycle!Nothing to do for density or non prognostic fields
                  if (have_option(trim(option_path)//"/aliased" )) cycle
                  option_path_BAK = "/material_phase["// int2str( i - 1 )//"]/scalar_field::BAK_"//trim(option_name)
                  !If a prognostic field then create backup (for incompressible flows, pressure would not be necessary)
                  call copy_option(trim(option_path),trim(option_path_BAK))!Better to copy twice and delete options because it is more robust that add_option
                  !Delete the prognostic section
                  call delete_option(trim(option_path_BAK)//"/prognostic")
                  !Now copy the interior of the prognostic as prescribed
                  call copy_option(trim(option_path)//"/prognostic",trim(option_path_BAK)//"/prescribed")
                  !Make sure the field is not shown
                  if (.not.have_option(trim(option_path_BAK)//"/prescribed/output/exclude_from_vtu")) then
                      !Copy an option that always exists to ensure we exclude the new field from vtu
                      call copy_option("/simulation_name",trim(option_path_BAK)//"/prescribed/output/exclude_from_vtu")
                  end if
                end do
                !Finally velocity as well
                option_path = "/material_phase["// int2str( i - 1 )//"]/vector_field::Velocity"
                option_path_BAK = "/material_phase["// int2str( i - 1 )//"]/vector_field::BAK_Velocity"
                call copy_option(trim(option_path),trim(option_path_BAK))
                !Delete the prognostic section
                call delete_option(trim(option_path_BAK)//"/prognostic")
                !Now copy the interior of the prognostic as prescribed
                call copy_option(trim(option_path)//"/prognostic",trim(option_path_BAK)//"/prescribed")
                !Make sure the field is not shown
                if (.not.have_option(trim(option_path_BAK)//"/prescribed/output/exclude_from_vtu")) then
                    !Copy an option that always exists to ensure we exclude the new field from vtu
                    call copy_option("/simulation_name",trim(option_path_BAK)//"/prescribed/output/exclude_from_vtu")
                end if
            end do
        end if

        if (have_option('/physical_parameters/black-oil_PVT_table')) then

        !Maybe no need for a field in state? and just calculate using required conditions????

            !Create necessary memory to store the mass fraction of one of the pseudo-components (in a way that it is also adapted)
            if (nphase /= 3)then
                FLAbort('Black-Oil modelling requires three phases. Phase 1 Aqua, phase 2 liquid, phase 3 vapour')
            end if
            if (GetProcNo() == 1) then
              print*, "WARNING: Currently for Black-Oil modelling it is recommended to disable Test_mass_consv by setting a high value. For example 1e10."
            end if


            option_path = "/material_phase["// int2str( nphase -1 )//"]/scalar_field::VapourMassFraction"
            call copy_option("/material_phase["// int2str( nphase - 1 )//"]/scalar_field::PhaseVolumeFraction",&
                 trim(option_path))
            !Make sure the field is not shown
            if (.not.have_option(trim(option_path)//"/prognostic/output/exclude_from_vtu")) then
                !Don't know how to set exclude_from_vtu to true from the spud options, hence,
                !since simulation_name HAS to exist I copy it to obtain the same effect
                call copy_option("/simulation_name",&
                 trim(option_path)//"/prognostic/output/exclude_from_vtu")
            end if
            !Make sure that this field is not the objective of adaptivity
            if (have_option(trim(option_path)//"/prognostic/adaptivity_options")) then
                call delete_option(trim(option_path)//"/prognostic/adaptivity_options")
            end if
        end if


        if (is_porous_media .and. have_option('/mesh_adaptivity/hr_adaptivity')) then
            ewrite(1, *) "Preserve regions MUST to be ON. Check multiphase_prototype_wrapper"
            !Ensure that preserve_mesh_regions is on, since otherwise it does not work
            !Don't know how to set exclude_from_vtu to true from the spud options, hence,
            !since simulation_name HAS to exist I copy it to obtain the same effect
            if (.not.have_option("/mesh_adaptivity/hr_adaptivity/preserve_mesh_regions")) then
                call copy_option("/simulation_name",&
                 "/mesh_adaptivity/hr_adaptivity/preserve_mesh_regions")
            end if
        end if
        if (is_porous_media) then
            ewrite(1, *) "For porous media we output only the Darcy Velocity, not the velocity. Check multiphase_prototype_wrapper"
            !Create a field to store the DarcyVelocity in it
            do i = 1, nphase
                option_path = "/material_phase["// int2str( i - 1 )//"]/vector_field::DarcyVelocity"
                if (.not.have_option(option_path)) then
                    call add_option(trim(option_path),  stat=stat)
                    option_path = "/material_phase["// int2str( i - 1 )//"]/vector_field::DarcyVelocity/prescribed"
                    call add_option(trim(option_path)//"/mesh::VelocityMesh",  stat=stat)
                    call add_option(trim(option_path)//"/value::WholeMesh",  stat=stat)
                    call add_option(trim(option_path)//"/value::WholeMesh/no_initial_condition",  stat=stat)
                    call add_option(trim(option_path)//"/output",  stat=stat)
                    call add_option(trim(option_path)//"/stat",  stat=stat)
                    call add_option(trim(option_path)//"/stat/include_in_stat",  stat=stat)

                    call add_option(trim(option_path)//"/detectors",  stat=stat)
                    call add_option(trim(option_path)//"/detectors/exclude_from_detectors",  stat=stat)
                    call add_option(trim(option_path)//"/do_not_recalculate",  stat=stat)
                    !Velocity is the force density which is pretty much useless so we instead show the DarcyVelocity
                    !do_not_show velocity

                    if (.not.have_option("/numerical_methods/porous_output_force_density") .and.&
                    .not.have_option("/material_phase["// int2str( i - 1 )//"]/vector_field::Velocity/prognostic/output/exclude_from_vtu"))&
                    call copy_option("simulation_name", &
                        "/material_phase["// int2str( i - 1 )//"]/vector_field::Velocity/prognostic/output/exclude_from_vtu")
                end if


!                option_path = "/material_phase["// int2str( i - 1 )//"]/vector_field::"
!                call copy_option(trim(option_path)//"Velocity", trim(option_path)//"DarcyVelocity")
!                if (have_option(trim(option_path)//"DarcyVelocity"//"/prognostic/tensor_field::Viscosity")) &
!                        call delete_option(trim(option_path)//"DarcyVelocity"//"/prognostic/tensor_field::Viscosity")
!
!                if (have_option(trim(option_path)//"DarcyVelocity"//"/prognostic/vector_field::Absorption"))&
!                        call delete_option(trim(option_path)//"DarcyVelocity"//"/prognostic/vector_field::Absorption")
!
!                if (have_option(trim(option_path)//"DarcyVelocity"//"/prognostic/adaptivity_options"))&
!                        call delete_option(trim(option_path)//"DarcyVelocity"//"/prognostic/adaptivity_options")
            end do
        end if

        !Add dummy fields to ensure that the well geometries are preserved when the mesh is adapted
        !or to show the wells in paraview
        if (npres>1 ) then
            ewrite(1, *) "For wells we define dummy variables and we enforce options in Diamond. Check multiphase_prototype_wrapper"
            if (have_option('/porous_media/wells_and_pipes/well_volume_ids')) then
                !Introduce some dummy regions to ensure that mesh adaptivity keeps the wells in place
                shape = option_shape('/porous_media/wells_and_pipes/well_volume_ids')
                assert(shape(1) >= 0)
                allocate(well_ids(shape(1)))
                call get_option( '/porous_media/wells_and_pipes/well_volume_ids', well_ids)
                !Create field by adding the fields manually
                option_path = "/porous_media/wells_and_pipes/scalar_field::Well_domains"
                call add_option(trim(option_path),  stat=stat)
                call add_option(trim(option_path)//"/prescribed",  stat=stat)
                call add_option(trim(option_path)//"/prescribed/mesh::P0DG",  stat=stat)

                !Add dummy values only for the regions with wells, the others are set to zero by default which is fine
                do i = 1, size(well_ids,1)
                    call add_option(trim(option_path)//"/prescribed/value::"//"Well_domain_"//int2str(well_ids(i)),  stat=stat)
                    call add_option(trim(option_path)//"/prescribed/value["//int2str(i-1)//"]/region_ids",  stat=stat)
                    call set_option(trim(option_path)//"/prescribed/value["//int2str(i-1)//"]/region_ids", (/well_ids(i)/),  stat=stat)
                    call add_option(trim(option_path)//"/prescribed/value["//int2str(i-1)//"]/constant",  stat=stat)
                    call set_option(trim(option_path)//"/prescribed/value["//int2str(i-1)//"]/constant", real(well_ids(i)), stat=stat)
                end do
                deallocate(well_ids)

                !It is important that the DiameterPipe is not recalculated if wells are defined using files
                if (have_option('/porous_media/wells_and_pipes/well_from_file[0]') .and. &
                    .not.have_option('/porous_media/wells_and_pipes/scalar_field::DiameterPipe/prescribed/do_not_recalculate'))&
                    call copy_option("/simulation_name", '/porous_media/wells_and_pipes/scalar_field::DiameterPipe/prescribed/do_not_recalculate')

                if (.not.have_option('/porous_media/wells_and_pipes/well_volume_ids/Show_well_volumes_ids'))&
                    call copy_option("/simulation_name", trim(option_path)//"/prescribed/output/exclude_from_vtu")

                call add_option(trim(option_path)//"/prescribed/do_not_recalculate",  stat=stat)
             else if (have_option('/porous_media/wells_and_pipes/well_from_file[0]')) then
                ewrite(0, *) "WARNING: well trajectories may not be preserved after mesh adaptivity. It is recommended to use /wells_and_pipes/well_volume_ids"
             end if
        end if
!print all the options in the diamond file and added here to the terminal
        !call print_options()
        !Call fluidity to populate state
        call populate_state(state)

    end subroutine populate_multi_state


    subroutine read_fluid_and_rock_properties_from_csv()
        !This subroutine reads the properties from the given input file, introduced in diamond
        !and create the correspoding entries in spud from inside the code.
        !The input files are as follows:
        !Property, Phase, region_id, value
        !Brooks_Corey/relperm_max,1,25,0.6
        !Permeability,0,24,1.;0.;0.;0.5
        implicit none
        !Local variables
        real :: value
        real, dimension(:,:), allocatable :: permeabilities
        integer :: Nentries, i, iphase, bar_pos, start, end, k, j, stat, Number_region_ids, ndim
        integer, dimension(10000) :: region_ids
        character( len = option_path_len ) :: path_to_file, diamond_path, property, cadena
        character(len=option_path_len), dimension(:,:),  allocatable :: csv_table_strings



        if (have_option('/io/PropertiesFromFile')) then
            call get_option('/geometry/dimension',ndim)
            allocate(permeabilities(ndim,ndim))
            !Get filepath
            call get_option('/io/PropertiesFromFile', path_to_file)

            call extract_strings_from_csv_file(csv_table_strings, path_to_file, Nentries)

            do i = 1, Nentries
                region_ids = -1
                property = trim(csv_table_strings(i,1))
                bar_pos = index(property,'-')
                !Capillary pressure is a special case
                if (index(property,'capillary_pressure') > 0 )then
                    !Unify the first two entries
                    property = property(1:bar_pos-1)//'/'//property(bar_pos+1:len_trim(property))
                    bar_pos = index(property,'-')
                end if
                !Thermal_porous data is also an special case
                if (index(property,'thermal_porous') > 0 )then
                    !Unify the first two entries
                    property = property(1:bar_pos-1)//'/'//property(bar_pos+1:len_trim(property))
                    bar_pos = index(property,'-')
                end if
                !Extract phase
                read(csv_table_strings(i,2) , *) iphase
                !Extract region ids
                cadena = csv_table_strings(i,3)
                end = min(len_trim(cadena)+1,max(0,index(cadena,'_')))
                if (end ==0) end = len_trim(cadena)+1
                k = 1
                do while (end>1)
                    read(cadena(1:end-1) , *) region_ids(k)
                    !trim cadena
                    cadena = cadena(end+1:len_trim(cadena))
                    !Calculate new end
                    end = min(len_trim(cadena)+1,max(0,index(cadena,'_')))
                    if (end ==0) end = len_trim(cadena)+1
                    !Advance index
                    k = k + 1
                end do
                Number_region_ids = k - 1
                if (bar_pos>0) then
                    !It is a fluid property
                    diamond_path = '/multiphase_properties/'//property(1:bar_pos-1)//'/scalar_field::'//property(bar_pos+1:len_trim(property))//'/prescribed/value'
                    !Get value
                    read(csv_table_strings(i,4) , *) value

                    !Proceed to first remove the present options, and next to add the new options
                    if (have_option('/material_phase['//int2str(iphase-1)//']'//trim(diamond_path)//'::WholeMesh/csv_file')) &
                    call delete_option('/material_phase['//int2str(iphase-1)//']'//trim(diamond_path)//'::WholeMesh')
                    !Next proceed to add the new values
                    if (region_ids(1) == 0) then
                        call add_option('/material_phase['//int2str(iphase-1)//']'//trim(diamond_path)//'::WholeMesh/constant',  stat=stat)
                        call set_option('/material_phase['//int2str(iphase-1)//']'//trim(diamond_path)//'::WholeMesh/constant', value,  stat=stat)
                    else                                                                            !Add name and also values
                        call add_option('/material_phase['//int2str(iphase-1)//']'//trim(diamond_path)//'::CSV_Property_'//int2str(i-1),  stat=stat)
                        call add_option('/material_phase['//int2str(iphase-1)//']'//trim(diamond_path)//'::CSV_Property_'//int2str(i-1)//"/region_ids",  stat=stat)
                        call set_option('/material_phase['//int2str(iphase-1)//']'//trim(diamond_path)//'::CSV_Property_'//int2str(i-1)//"/region_ids", region_ids(1:Number_region_ids),  stat=stat)
                        call add_option('/material_phase['//int2str(iphase-1)//']'//trim(diamond_path)//'::CSV_Property_'//int2str(i-1)//"/constant",  stat=stat)
                        call set_option('/material_phase['//int2str(iphase-1)//']'//trim(diamond_path)//'::CSV_Property_'//int2str(i-1)//"/constant", value, stat=stat)
                    end if
                    !NOTHING FOR WELLS AS IT IS NOT A POROUS MEDIUM
                else!Porous media property

                    !If permeability, then tensor field, otherwise scalar field
                    if (index(property,'Permeability') <= 0 )then
                        diamond_path = '/porous_media/scalar_field::Porosity/prescribed/value'
                        !Get value
                        read(csv_table_strings(i,4) , *) value

                        !Proceed to first remove the present options, and next to add the new options
                        if (have_option(trim(diamond_path)//'::WholeMesh/csv_file')) call delete_option(trim(diamond_path)//'::WholeMesh')
                        !Next proceed to add the new values
                        if (region_ids(1) == 0) then
                            call add_option(trim(diamond_path)//'::WholeMesh/constant',  stat=stat)
                            call set_option(trim(diamond_path)//'::WholeMesh/constant', value,  stat=stat)
                        else                                    !Add name and also values
                            call add_option(trim(diamond_path)//'::CSV_Petrophysical_property_'//int2str(i-1),  stat=stat)
                            call add_option(trim(diamond_path)//'::CSV_Petrophysical_property_'//int2str(i-1)//"/region_ids",  stat=stat)
                            call set_option(trim(diamond_path)//'::CSV_Petrophysical_property_'//int2str(i-1)//"/region_ids", region_ids(1:Number_region_ids),  stat=stat)
                            call add_option(trim(diamond_path)//'::CSV_Petrophysical_property_'//int2str(i-1)//"/constant",  stat=stat)
                            call set_option(trim(diamond_path)//'::CSV_Petrophysical_property_'//int2str(i-1)//"/constant", value, stat=stat)
                        end if
                    else!Permeability
                        !Extract permeability entries
                        cadena = csv_table_strings(i,4)
                        end = min(len_trim(cadena)+1,max(0,index(cadena,'_')))
                        if (end == 0) end = len_trim(cadena)+1
                        do k = 1, ndim
                            do j = 1, ndim

                                read(cadena(1:end-1) , *) value
                                permeabilities(k,j) = value
                                !trim cadena
                                cadena = cadena(end+1:len_trim(cadena))
                                            !Calculate new end
                                end = min(len_trim(cadena)+1,max(0,index(cadena,'_')))
                                if (end ==0) end = len_trim(cadena)+1
                            end do
                        end do
                        diamond_path = '/porous_media/tensor_field::Permeability/prescribed/value'
                        !Proceed to first remove the present options, and next to add the new options
                        if (have_option(trim(diamond_path)//'::WholeMesh/anisotropic_asymmetric/csv_file')) &
                            call delete_option(trim(diamond_path)//'::WholeMesh')
                        !Next proceed to add the new values
                        if (region_ids(1) == 0) then
                            call add_option(trim(diamond_path)//'::WholeMesh/anisotropic_asymmetric/constant',  stat=stat)
                            call set_option(trim(diamond_path)//'::WholeMesh/anisotropic_asymmetric/constant', value,  stat=stat)
                        else                                    !Add name and also values
                            call add_option(trim(diamond_path)//'::CSV_Permeability_'//int2str(i-1),  stat=stat)
                            call add_option(trim(diamond_path)//'::CSV_Permeability_'//int2str(i-1)//"/region_ids",  stat=stat)
                            call set_option(trim(diamond_path)//'::CSV_Permeability_'//int2str(i-1)//"/region_ids", region_ids(1:Number_region_ids),  stat=stat)
                            call add_option(trim(diamond_path)//'::CSV_Permeability_'//int2str(i-1)//"/anisotropic_asymmetric/constant",  stat=stat)
                            call set_option(trim(diamond_path)//'::CSV_Permeability_'//int2str(i-1)//"/anisotropic_asymmetric/constant", permeabilities, stat=stat)
                        end if
                    end if
                end if
            end do

            deallocate(permeabilities)
        else
            return
        end if
    end subroutine read_fluid_and_rock_properties_from_csv

    subroutine autocomplete_input_file(nphase, npres, ncomp)
                ! ################################
        ! In this subroutine. We populate the input file elements that are not user-selected in IC_FERST.rnc. (i.e. we set all the things that the user doesn't need to worry about).
        !IT IS VERY IMPORTANT NOT TO CHANGE THE ORDERING IN THIS SUBROUTINE!!!
        implicit none
        integer, intent(in) :: nphase, npres, ncomp
        !Local variables
        integer :: stat, i, k, simulation_quality = 10, scalarComponents
        real :: aux
        character( len = option_path_len ) :: option_path, quality_option

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

        !#########GEOMETRY AND PRECISION OPTIONS#################
        !SPRINT_TO_DO Make Fluidity to read the meshes from /geometry/Advance_options instead of having to copy them in the /geometry
        if (have_option("/geometry/Advance_options/mesh::VelocityMesh/")) Then
            !use the user input options
            call copy_option("/geometry/Advance_options/mesh::VelocityMesh", "/geometry/mesh::VelocityMesh", stat=stat)
        else
          option_path = "/geometry/mesh::VelocityMesh/"
          call add_option(trim(option_path)//"from_mesh", stat=stat)
          call add_option(trim(option_path)//"from_mesh/mesh::CoordinateMesh", stat=stat)
          call add_option(trim(option_path)//"from_mesh/mesh_continuity", stat=stat)
          call set_option(trim(option_path)//"from_mesh/mesh_continuity", "discontinuous")
          call add_option(trim(option_path)//"from_mesh/mesh_shape/element_type", stat=stat)
          call set_option(trim(option_path)//"from_mesh/mesh_shape/element_type", "lagrangian")
          call add_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", stat=stat)
          if (simulation_quality < 100) then
              if (have_option("/porous_media_simulator")) then
                call set_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", 0)
              else !Currently only for porous media P0DG works, so we use P1 for the rest
                call set_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", 1)
              end if
              call set_option(trim(option_path)//"from_mesh/mesh_shape/element_type", "lagrangian")
          else if (simulation_quality < 1000) then
              call set_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", 1)
              call set_option(trim(option_path)//"from_mesh/mesh_shape/element_type", "lagrangian")
          else
              call set_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", 1)
              call set_option(trim(option_path)//"from_mesh/mesh_shape/element_type", "bubble")
          end if
          call add_option(trim(option_path)//"from_mesh/stat/exclude_from_stat", stat=stat)
        end if

        if (have_option("/geometry/Advance_options/mesh::PressureMesh/")) Then
            !use the user input options
            call copy_option("/geometry/Advance_options/mesh::PressureMesh", "/geometry/mesh::PressureMesh", stat=stat)
        else
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
        end if

        !Extra meshes
        k = option_count("/geometry/Advance_options/mesh")
        do i =1, k
          call get_option("/geometry/Advance_options/mesh["// int2str( i - 1 )//"]/name",option_path)
          if (trim(option_path)/="VelocityMesh" .and. trim(option_path)/="PressureMesh") then
            !Put in the place where Fluidity expects the mesh to be
            call copy_option("/geometry/Advance_options/mesh::"//trim(option_path), "/geometry/mesh::"//trim(option_path), stat=stat)
          end if
        end do

        if (have_option("/physical_parameters/gravity/hydrostatic_pressure_solver")) then
          !Introduce the HydrostaticPressure mesh, quadratic and continuous
          option_path = "/geometry/mesh::HydrostaticPressure/"
          call add_option(trim(option_path)//"from_mesh", stat=stat)
          call add_option(trim(option_path)//"from_mesh/mesh::CoordinateMesh", stat=stat)
          call add_option(trim(option_path)//"from_mesh/mesh_continuity", stat=stat)
          call set_option(trim(option_path)//"from_mesh/mesh_continuity", "continuous")
          call add_option(trim(option_path)//"from_mesh/mesh_shape/element_type", stat=stat)
          call set_option(trim(option_path)//"from_mesh/mesh_shape/element_type", "lagrangian")
          call add_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", stat=stat)
          call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/polynomial_degree', k )
          if (k == 0) then !For P0DG we use a hydrostatic pressure solver of order one. sprint_to_do => this is because it is unfinished but I am not sure if it is worth it
            call set_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", 1)
          else
            call set_option(trim(option_path)//"from_mesh/mesh_shape/polynomial_degree", 2)
          end if
          call add_option(trim(option_path)//"from_mesh/stat/exclude_from_stat", stat=stat)
          !Do we need the user to create a HydrostaticPressure scalar field? maybe not
        end if

        !#########GEOMETRY AND PRECISION OPTIONS#################

!Sprint_to_do
!use allocate_and_insert_auxilliary_fields as a template to include internally the fields into state, without actually changing the input fiel using spud


        ! Linear Solver options - Iterative Method; default options for all fields
        if (.not. have_option("/solver_options/Linear_solver" )) then!The default has to exist
          if (GetProcNo() == 1) then
            print*, "MESSAGE: Using default options for the linear solver: GMRES(30) + HYPRE"
          end if
          option_path = "/solver_options/Linear_solver/iterative_method::gmres/restart"
          call add_option(trim(option_path), stat = stat)
          call set_option(trim(option_path), 30)

          !Preconditioner
          option_path = "/solver_options/Linear_solver/preconditioner::hypre/"
          call add_option(trim(option_path), stat = stat)
          option_path = "/solver_options/Linear_solver/preconditioner::hypre/hypre_type::boomeramg"
          call add_option(trim(option_path), stat = stat)

          !Convergence settings
          call add_option("/solver_options/Linear_solver/relative_error", stat = stat)
          call set_option("/solver_options/Linear_solver/relative_error", 1e-10)
          call add_option("/solver_options/Linear_solver/absolute_error", stat = stat)
          call set_option("/solver_options/Linear_solver/absolute_error", 1e-8)
          call add_option("/solver_options/Linear_solver/max_iterations", stat = stat)
          call set_option("/solver_options/Linear_solver/max_iterations", 300)
          !Copy an option that always exists to ensure ignore all solver failues
          call copy_option("/simulation_name","/solver_options/Linear_solver/ignore_all_solver_failures")

        end if

        !Non_Linear_Solver default settings
        option_path = "/solver_options/Non_Linear_Solver"
        if (.not. have_option(trim(option_path) )) then!The default has to exist
          if (GetProcNo() == 1) then
            print*, "MESSAGE: Using default options for the non-linear solver. Convergence check (Pressure) with a maximum of 20 its. ",&
             "For multiphase porous media VAD and automatic backtracking are on. The check is either Saturation, if multiphase, or the Tracer used."
          end if
          call add_option(trim(option_path), stat = stat)
          call set_option(trim(option_path), 20)
          option_path = trim(option_path)//"/Fixed_Point_Iteration"
          call add_option(trim(option_path), stat = stat)
          call set_option(trim(option_path), 5e-2)
          call add_option(trim(option_path)//"/Infinite_norm_tol", stat = stat)
          call set_option(trim(option_path)//"/Infinite_norm_tol", 0.01)
          if (have_option("/porous_media_simulator")) then
            !Add VAD options
            call add_option(trim(option_path)//"/Vanishing_relaxation", stat=stat)
            call set_option(trim(option_path)//"/Vanishing_relaxation",-1e2)
            call add_option(trim(option_path)//"/Vanishing_relaxation/Vanishing_for_transport", stat=stat)
            call set_option(trim(option_path)//"/Vanishing_relaxation/Vanishing_for_transport",-1e1)
            !If multiphase then the important thing is saturation!
            if ( nphase > 1.) then
              !Set up the non-linear solver to use an automatic backtracking method
              call add_option(trim(option_path)//"/Backtracking_factor", stat=stat)
              call set_option(trim(option_path)//"/Backtracking_factor",-10.)
            else !single phase
              call add_option(trim(option_path)//"/Infinite_norm_tol/adaptive_non_linear_iterations", stat = stat)
              if (have_option('/material_phase[0]/scalar_field::Temperature')) then
                call set_option(trim(option_path)//"/Infinite_norm_tol/adaptive_non_linear_iterations", 4)
              elseif (have_option('/material_phase[0]/scalar_field::SoluteMassFraction')) then
                call set_option(trim(option_path)//"/Infinite_norm_tol/adaptive_non_linear_iterations", 5)
              else !If nothing, then pressure
                call set_option(trim(option_path)//"/Infinite_norm_tol/adaptive_non_linear_iterations", 1)
              end if

              call add_option(trim(option_path)//"/Impose_min_max", stat = stat)
            end if

          else
              call add_option("/solver_options/Linear_solver/relative_error", stat = stat)
              call set_option("/solver_options/Linear_solver/relative_error", 1e-10)
              !Use pressure to decide when to stop the non_linear iteration process
              call add_option(trim(option_path)//"/Infinite_norm_tol/adaptive_non_linear_iterations", stat = stat)
              call set_option(trim(option_path)//"/Infinite_norm_tol/adaptive_non_linear_iterations", 1)
          end if

        end if
        ! IO STAT OPTIONS
        option_path = "/io/output_mesh[0]/name"
        call add_option(trim(option_path), stat=stat)
        call set_option(trim(option_path),"PressureMesh")

        do i = 1, nphase*npres + ncomp*nphase

          !Include that fields are part of stat, and not of detectors
          do k = 1, option_count("/material_phase[" // int2str(i-1) // "]/scalar_field")
              option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field["// int2str( k - 1)//"]/prognostic"
              if (have_option(trim(option_path))) then!Only for prognostic fields
                ! Stat, convergence, detectors, steady state settings
                !We have removd this from scalar_fields
                !SPRINT_TO_DO THESE OPTIONS TO BE REMOVED from here and the other fields, this seems unnecessary?
                call add_option(trim(option_path)//"/stat/include_in_stat", stat=stat)
                call add_option(trim(option_path)//"/convergence/include_in_convergence", stat=stat)
                call add_option(trim(option_path)//"/detectors/exclude_from_detectors", stat=stat)
              end if
          end do


          !Components are treated differently, but this should change! SPRINT_TO_DO
          if (.not.have_option("/material_phase["// int2str( i - 1 )//"]/is_multiphase_component")) then
            !Create memory for Density internally
            option_path = "/material_phase["// int2str( i - 1 )//"]/scalar_field::Density"
            if (.not.have_option(option_path)) then
                call add_option(trim(option_path),  stat=stat)
                option_path = "/material_phase["// int2str( i - 1 )//"]/scalar_field::Density/prognostic"
                call add_option(trim(option_path)//"/mesh::PressureMesh",  stat=stat)
                call add_option(trim(option_path)//"/value::WholeMesh",  stat=stat)
                call add_option(trim(option_path)//"/value::WholeMesh/no_initial_condition",  stat=stat)
                call add_option(trim(option_path)//"/output",  stat=stat)
                call add_option(trim(option_path)//"/stat",  stat=stat)
                call add_option(trim(option_path)//"/stat/include_in_stat",  stat=stat)

                call add_option(trim(option_path)//"/detectors",  stat=stat)
                call add_option(trim(option_path)//"/detectors/exclude_from_detectors",  stat=stat)
                call add_option(trim(option_path)//"/do_not_recalculate",  stat=stat)
            end if
!Easiest way to create the viscosity field is to move where it was inside velocity!SPRINT_TO_DO NEED TO CHANGE THIS!
            call copy_option("/material_phase["// int2str( i - 1 )//"]/phase_properties/Viscosity/tensor_field::Viscosity",&
            "/material_phase["// int2str( i - 1 )//"]/vector_field::Velocity/prognostic/tensor_field::Viscosity")
  !CHECK BECAUSE MAYBE THESE MEMORY IS AUTOMATICALLY ALLOCATED
!Easiest way to create the diffusivity field is to move where it was inside velocity!SPRINT_TO_DO NEED TO CHANGE THIS!
            if (have_option("/material_phase["// int2str( i - 1 )//"]/phase_properties/tensor_field::Thermal_Conductivity")) then
              if (.not. have_option ("/material_phase["// int2str( i - 1 )//"]/scalar_field::Temperature/prognostic")) then
                  FLAbort("Thermal Conductivity specified but no prognostic temperature field specified.")
              end if
                call copy_option("/material_phase["// int2str( i - 1 )//"]/phase_properties/tensor_field::Thermal_Conductivity",&
                  "/material_phase["// int2str( i - 1 )//"]/scalar_field::Temperature/prognostic/tensor_field::Diffusivity")!SPRINT_TO_DO NAME THIS THERMAL_CONDUCTIVITY
            end if
!Easiest way to create the heatcapacity field is to move where it was inside velocity!SPRINT_TO_DO NEED TO CHANGE THIS!
            if (have_option("/material_phase["// int2str( i - 1 )//"]/phase_properties/scalar_field::HeatCapacity")) then
                if (.not. have_option ("/material_phase["// int2str( i - 1 )//"]/scalar_field::Temperature/prognostic")) then
                    FLAbort("HeatCapacity specified but no prognostic temperature field specified.")
                end if
                call copy_option("/material_phase["// int2str( i - 1 )//"]/phase_properties/scalar_field::HeatCapacity",&
                  "/material_phase["// int2str( i - 1 )//"]/scalar_field::Temperature/prognostic/scalar_field::HeatCapacity")
            end if
            !Easiest way to create the diffusivity field is to move where it was inside velocity!SPRINT_TO_DO NEED TO CHANGE THIS!
            if (have_option("/material_phase["// int2str( i - 1 )//"]/phase_properties/tensor_field::Solute_Diffusivity")) then
              if (.not. have_option ("/material_phase["// int2str( i - 1 )//"]/scalar_field::SoluteMassFraction/prognostic")) then
                  FLAbort("Solute Diffusivity specified but no prognostic SoluteMassFraction field specified.")
              end if
                call copy_option("/material_phase["// int2str( i - 1 )//"]/phase_properties/tensor_field::Solute_Diffusivity",&
                  "/material_phase["// int2str( i - 1 )//"]/scalar_field::SoluteMassFraction/prognostic/tensor_field::Diffusivity")!SPRINT_TO_DO NAME THIS THERMAL_CONDUCTIVITY
            end if

            if (have_option("/physical_parameters/gravity/hydrostatic_pressure_solver") .and. i == 1) then
              !Add a prognostic field named HydrostaticPressure (do we need BCs or initial conditions for this?)
              !This is only required for the first phase
              option_path = "/material_phase["// int2str( i - 1 )//"]/scalar_field::HydrostaticPressure"
              if (.not.have_option (trim(option_path))) then
                call add_option(trim(option_path),  stat=stat)
                option_path = "/material_phase["// int2str( i - 1 )//"]/scalar_field::HydrostaticPressure/prognostic"

                call add_option(trim(option_path)//"/mesh::HydrostaticPressure",  stat=stat)
                call add_option(trim(option_path)//"/value::WholeMesh",  stat=stat)
                call add_option(trim(option_path)//"/value::WholeMesh/constant",  stat=stat)
                call set_option(trim(option_path)//"/value::WholeMesh/constant",  0.)
                call add_option(trim(option_path)//"/output",  stat=stat)
                call add_option(trim(option_path)//"/stat",  stat=stat)
                call add_option(trim(option_path)//"/stat/include_in_stat",  stat=stat)

                call add_option(trim(option_path)//"/detectors",  stat=stat)
                call add_option(trim(option_path)//"/detectors/exclude_from_detectors",  stat=stat)
                call add_option(trim(option_path)//"/do_not_recalculate",  stat=stat)

              end if

            end if


            ! SCALAR_FIELD(PRESSURE) OPTIONS ADDED AUTOMATICALLY
            option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::Pressure/prognostic"
            if (have_option(trim(option_path))) then
                !SPRINT_TO_DO THESE OPTIONS TO BE REMOVED from here and the other fields, this seems unnecessary?
                option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field::Pressure/prognostic/spatial_discretisation/continuous_galerkin"
                call add_option(trim(option_path), stat=stat)
            end if

                ! VECTOR_FIELD(VELOCITY) OPTIONS ADDED AUTOMATICALLY
            option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic"
            if (have_option(trim(option_path))) then

                !SPRINT_TO_DO THIS OPTION TO BE REMOVED
                option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin"
                call add_option(trim(option_path), stat=stat)

                ! Stat, convergence, detectors, steady state settings
                option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/stat"
                call add_option(trim(option_path), stat=stat)
                option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/convergence/include_in_convergence"
                call add_option(trim(option_path), stat=stat)
                option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/detectors/exclude_from_detectors"
                call add_option(trim(option_path), stat=stat)
                ! option_path = "/material_phase["// int2str( i - 1)//"]/vector_field::Velocity/prognostic/steady_state/include_in_steady_state"
                ! call add_option(trim(option_path), stat=stat)
            end if
          else

            scalarComponents = option_count("/material_phase[" // int2str(i-1) // "]/scalar_field")
            !Easiest way to move the fields in the root of the scalar field ComponentMassFraction. NEED TO CHANGE THIS!
            do k = 1, scalarComponents
              option_path = "/material_phase["// int2str( i - 1)//"]/scalar_field["// int2str( k - 1)//"]/prognostic"
              !Check fields and move then one level up
              if (have_option(trim(option_path)//"/phase_properties/Viscosity/tensor_field::Viscosity")) then
                !Move one level up, as it used to be, SPRINT_TO_DO, CREATE THE MEMORY BASED ON THIS NEW POSITION
                call copy_option(trim(option_path)//"/phase_properties/Viscosity/tensor_field::Viscosity",&
                trim(option_path)//"/tensor_field::Viscosity")
              end if
              if (have_option(trim(option_path)//"/phase_properties/tensor_field::Diffusivity")) then
                !Move one level up, as it used to be, SPRINT_TO_DO, CREATE THE MEMORY BASED ON THIS NEW POSITION
                call copy_option(trim(option_path)//"/phase_properties/tensor_field::Diffusivity",&
                trim(option_path)//"/tensor_field::Diffusivity")
              end if
              if (have_option(trim(option_path)//"/phase_properties/scalar_field::HeatCapacity")) then
                !Move one level up, as it used to be, SPRINT_TO_DO, CREATE THE MEMORY BASED ON THIS NEW POSITION
                call copy_option(trim(option_path)//"/phase_properties/scalar_field::HeatCapacity",&
                trim(option_path)//"/scalar_field::HeatCapacity")
              end if


            end do

          end if
        end do

        !####################################

    end subroutine autocomplete_input_file



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
        !By default it is inertia dominated
        is_porous_media = have_option('/porous_media_simulator') .or. have_option('/is_porous_media')
        is_magma = have_option('/magma_simulator')
        is_flooding = have_option('/flooding_simulator')
        is_poroelasticity = have_option('/poroelasticity')
        !Flag to set up the coupling with femdem
        is_multifracture = have_option( '/femdem_fracture' )
        !Flag to set up blasting
        is_blasting = have_option( '/blasting' )
        !Has temperature
        has_temperature = have_option( '/material_phase[0]/scalar_field::Temperature/' )
        !Arash
        has_salt = have_option( '/material_phase[0]/scalar_field::SoluteMassFraction/' )

        ! Check if Porous media model initialisation
        is_porous_initialisation =  have_option("/porous_media/FWL")

        !Tell the user which sort of simulation type he is running, just in case
        if (have_option('/inertia_dominated_simulator') .and. have_option('/porous_media')) then
          if (GetProcNo() == 1) then
            FLAbort("The simulator has been set up to inertia dominated but porous media options have been set up. ")
          end if
        end if




    end subroutine get_simulation_type

end subroutine multiphase_prototype_wrapper
