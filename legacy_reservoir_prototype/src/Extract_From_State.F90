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
!    but WITHOUT ANY WARRANTY; without seven the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module Copy_Outof_State
    !! This module enables the multiphase prototype code to interact with state by
    !! copying everything required from state to the MP-space.

    use fldebug
    use state_module
    use fields
    use field_options
    use spud
    use populate_state_module
    use diagnostic_variables
    use diagnostic_fields
    use diagnostic_fields_wrapper
    use global_parameters
    use diagnostic_fields_wrapper_new
    use element_numbering
    use shape_functions
    use fefields
    use boundary_conditions
    use futils, only: int2str
    use boundary_conditions_from_options
    use parallel_tools, only : allmax, allmin, isparallel
    use parallel_fields
    use memory_diagnostics
    use initialise_fields_module, only: initialise_field_over_regions
    use halos
    ! Need to use surface integrals since a function from here is called in the calculate_outflux() subroutine
    use surface_integrals

    use solvers
    use conservative_interpolation_module
    use multi_data_types
    use multi_tools
    implicit none

    private

    public :: Get_Primary_Scalars_new, Compute_Node_Global_Numbers, &
        Get_Ele_Type, Get_Discretisation_Options, inf_norm_scalar_normalised, &
        update_boundary_conditions, pack_multistate, finalise_multistate, get_ndglno, Adaptive_NonLinear,&
        get_var_from_packed_state, as_vector, as_packed_vector, is_constant, GetOldName, GetFEMName, PrintMatrix,&
        have_option_for_any_phase, Get_Ele_Type_new,&
        get_Convergence_Functional, get_DarcyVelocity, printCSRMatrix, dump_outflux, calculate_internal_volume, prepare_absorptions


    interface Get_SNdgln
       module procedure Get_Scalar_SNdgln, Get_Vector_SNdgln
    end interface Get_SNdgln

contains

   subroutine Get_Primary_Scalars_new( state, Mdims )
        !!$ This subroutine extracts all primary variables associated with the mesh from state,
        !!$ and associated them with the variables used in the MultiFluids model.
        implicit none
        type( state_type ), dimension( : ), intent( in ) :: state
        type (multi_dimensions) :: Mdims

        !!$ Local variables
        type( vector_field ), pointer :: positions, velocity
        type( scalar_field ), pointer :: pressure
        type( mesh_type ), pointer :: velocity_cg_mesh, pressure_cg_mesh, ph_mesh
        integer :: i, stat
        logical , save :: warning_displayed = .false.
        ewrite(3,*)' In Get_Primary_Scalars'

        !!$ Defining dimension and nstate
        call get_option( '/geometry/dimension', Mdims%ndim )
        Mdims%nstate = option_count( '/material_phase' )
        !!$ Assume there are the same number of components in each phase (will need to check this eventually)
        Mdims%ncomp = 0
        do i = 1, Mdims%nstate
            if( have_option( '/material_phase[' // int2str(i-1) // &
                ']/is_multiphase_component' ) ) then
                Mdims%ncomp = Mdims%ncomp + 1
            end if
        end do
        Mdims%nphase = Mdims%nstate - Mdims%ncomp
        assert( Mdims%nphase > 0 ) ! Check if there is more than 0 phases
        ! Number of pressures to solve for
        Mdims%npres = option_count("/material_phase/scalar_field::Pressure/prognostic")
        Mdims%n_in_pres = Mdims%nphase / Mdims%npres

        !!$ Get the vel element type.
        if (is_porous_media) then!Check that the FPI method is on
            if (.not. have_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration') .and. Mdims%n_in_pres > 1) then
                ewrite(0,*) "WARNING: The option <Fixed_Point_Iteration> is HIGHLY recommended for multiphase porous media flow."
            else!Check that the user is allowing the linear solver to fail
                if (.not. have_option( '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/'//&
                'solver/ignore_all_solver_failures') .and. .not.warning_displayed) then
                    ewrite(0,*) "WARNING: The option <PhaseVolumeFraction/prognostic/solver/ignore_all_solver_failures>"//&
                    " is HIGHLY recommended for multiphase porous media flow to allow the FPI method to find a solution."
                    warning_displayed = .true.
                end if
            end if
            !Don't use for single phase porous media flows
            if (have_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration') .and. Mdims%nphase < 2) then
                !Unless we are using dynamic control of the non-linear iterations in which case it does not matter
                if (.not.have_option('/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Infinite_norm_tol/adaptive_non_linear_iterations')) then
                    ewrite(0,*) "WARNING: The option <Fixed_Point_Iteration> SHOULD NOT be used for single phase porous media flows without wells."
                end if
            end if
        end if
        is_multifracture = have_option( '/simulation_type/femdem_fracture' ) .or. is_multifracture

        positions => extract_vector_field( state, 'Coordinate' )
        pressure_cg_mesh => extract_mesh( state, 'PressureMesh_Continuous' )

        !!$ Defining number of elements and surface elements, coordinates, locs and snlocs
        Mdims%totele = ele_count( positions )
        Mdims%stotel = surface_element_count( positions )

        !!$ Coordinates
        Mdims%x_nloc_p1 = ele_loc( positions, 1 )
        Mdims%x_nloc = ele_loc( pressure_cg_mesh, 1 )
        Mdims%x_snloc = face_loc( pressure_cg_mesh, 1 )
        Mdims%x_nonods_p1 = node_count( positions )
        Mdims%x_nonods = node_count( pressure_cg_mesh )

        !!$ Pressure, Control Volumes and Materials
        pressure => extract_scalar_field( state, 'Pressure' )
        Mdims%p_nloc = ele_loc( pressure, 1 )
        Mdims%p_snloc = face_loc( pressure, 1 )
        Mdims%p_nonods = node_count( pressure )
        Mdims%cv_nloc = Mdims%p_nloc
        Mdims%cv_snloc = Mdims%p_snloc
        Mdims%cv_nonods = node_count( pressure )
        Mdims%mat_nloc = Mdims%cv_nloc
        Mdims%mat_nonods = Mdims%mat_nloc * Mdims%totele

        !!$ Velocities and velocities (DG) associated with the continuous space (CG)
        velocity => extract_vector_field( state, 'Velocity' )
        Mdims%u_nloc = ele_loc( velocity, 1 )
        Mdims%u_snloc = face_loc( velocity, 1 )
        Mdims%u_nonods = node_count( velocity )

        !!$ Get the continuous space of the velocity field
        if (.not.is_P0DGP1CV) then
            velocity_cg_mesh => extract_mesh( state, 'VelocityMesh_Continuous' )
            Mdims%xu_nloc = ele_loc( velocity_cg_mesh, 1 )
            Mdims%xu_nonods = max(( Mdims%xu_nloc - 1 ) * Mdims%totele + 1, Mdims%totele )
        end if
        if( have_option( "/physical_parameters/gravity/hydrostatic_pressure_solver" ) ) then
            ph_mesh => extract_mesh( state( 1 ), 'ph', stat )
            if ( stat == 0 ) then
                Mdims%ph_nloc = ele_loc( ph_mesh, 1 )
                Mdims%ph_nonods = node_count( ph_mesh )
            else
                FLAbort("You need a 'ph' mesh to use the high-order hydrostatic pressure solver.")
            end if
        else
            Mdims%ph_nloc = 0
            Mdims%ph_nonods = 0
        end if

        return
    end subroutine Get_Primary_Scalars_new

    subroutine Compute_Node_Global_Numbers( state, ndgln)
        !!$ This subroutine calculates the global node numbers requested to operates in the MP-space.
        implicit none
        type( state_type ), dimension( : ), intent( in ) :: state
        type(multi_ndgln), intent(inout) :: ndgln
        !Local variables
        type( vector_field ), pointer :: velocity
        type( scalar_field ), pointer :: pressure
        !Point to state
        ndgln%x_p1=>get_ndglno(extract_mesh(state(1),"CoordinateMesh"))
        ndgln%x=>get_ndglno(extract_mesh(state(1),"PressureMesh_Continuous"))
        ndgln%cv=>get_ndglno(extract_mesh(state(1),"PressureMesh"))
        ndgln%p=>get_ndglno(extract_mesh(state(1),"PressureMesh"))
        ndgln%mat=>get_ndglno(extract_mesh(state(1),"PressureMesh_Discontinuous"))
        ndgln%u=>get_ndglno(extract_mesh(state(1),"InternalVelocityMesh"))
        if (.not.is_P0DGP1CV) ndgln%xu=>get_ndglno(extract_mesh(state(1),"VelocityMesh_Continuous"))
        !!$ Pressure, control volume and material
        pressure => extract_scalar_field( state( 1 ), 'Pressure' )
        !!$ Velocities
        velocity => extract_vector_field( state( 1 ), 'Velocity' )
        !!$ Surface-based global node numbers for control volumes and pressure
        call Get_SNdgln( ndgln%suf_cv, pressure )
        ndgln%suf_p = ndgln%suf_cv
        !!$ Velocities
        call Get_SNdgln( ndgln%suf_u, velocity )
        return
    end subroutine Compute_Node_Global_Numbers

    subroutine Get_Ele_Type( x_nloc, cv_ele_type, p_ele_type, u_ele_type, &
        mat_ele_type, u_sele_type, cv_sele_type )
        !-
        !- u_ele_type = cv_ele_type = p_ele_type will flag the dimension and
        !- type of element:
        !- = 1 or 2: 1D (linear and quadratic, respectively)
        !- = 3 or 4: triangle (linear or quadratic, respectively)
        !- = 5 or 6: quadrilateral (bi-linear or tri-linear, respectively)
        !- = 7 or 8: tetrahedron (linear or quadratic, respectively)
        !- = 9 or 10: hexahedron (bi-linear or tri-linear, respectively)
        !-
        implicit none

        integer, intent(in) :: x_nloc
        integer, intent( inout ) :: cv_ele_type, p_ele_type, u_ele_type
        integer, intent( inout ), optional :: mat_ele_type, u_sele_type, cv_sele_type
        !!$ Local variables
        integer :: ndim, degree

        call get_option( '/geometry/dimension', ndim)

        call get_option( &
            '/geometry/mesh::PressureMesh/from_mesh/mesh_shape/polynomial_degree', &
            degree )

        Select Case( ndim )
            case( 1 ) ! ndim

                Select Case( degree )

                    case( 1 ) ! degree

                        !!$ ndim=1; p=1
                        cv_ele_type = 1

                    case( 2 ) ! degree

                        ! ndim=1; p=2
                        cv_ele_type = 2

                    case default; FLAbort('Degree error')

                end Select ! degree

            case( 2 ) ! ndim

                Select Case( degree )

                    case( 1 ) ! degree

                        Select Case( x_nloc )

                            case( 3 ) ! x_nloc

                                ! ndim=2; p=1; x_nloc=3
                                ! linear triangle
                                cv_ele_type = 3

                            case( 4 ) ! x_nloc

                                ! ndim=2; p=1; x_nloc=4
                                ! bilinear quad
                                cv_ele_type = 5

                            case default; FLAbort('X_nloc error')

                        end Select ! x_nloc

                    case( 2 ) ! degree

                        Select Case( x_nloc )

                            case( 6 ) ! x_nloc

                                ! ndim=2; p=2; x_nloc=3
                                ! quadratic triangle
                                cv_ele_type = 4

                            case( 10 ) ! x_nloc

                                ! ndim=2; p=2; x_nloc=4
                                ! bi-quadratic quad
                                cv_ele_type = 6

                            case default; FLAbort('X_nloc error')

                        end Select ! x_nloc

                    case default; FLAbort('Degree error')

                end Select ! degree

            case( 3 ) ! ndim

                Select Case( degree )

                    case( 1 ) ! degree

                        Select Case( x_nloc )

                            case( 4 ) ! x_nloc

                                ! ndim=3; p=1; x_nloc=4
                                ! linear tets
                                cv_ele_type = 7

                            case( 8 ) ! x_nloc

                                ! ndim=3; p=1; x_nloc=8
                                ! tri-linear hex
                                cv_ele_type = 9

                            case default; FLAbort('X_nloc error')

                        end Select ! x_nloc

                    case( 2 ) ! degree

                        Select Case( x_nloc )

                            case( 10 ) ! x_nloc

                                ! ndim=3; p=2; x_nloc=4
                                ! quadratic tet
                                cv_ele_type = 8

                            case( 27 ) ! x_nloc

                                ! ndim=3; p=2; x_nloc=8
                                ! bilinear quad
                                cv_ele_type = 10

                            case default; FLAbort('X_nloc error')

                        end Select ! x_nloc

                    case default; FLAbort('Degree error')

                end Select ! degree

        end Select ! ndim

        p_ele_type = cv_ele_type ; u_ele_type = cv_ele_type

        !!$ The following options are hardcoded and need to be either deleted from the code tree or
        !!$ added into the schema.
        if( present( mat_ele_type ) ) mat_ele_type = 1
        if( present( u_sele_type ) ) u_sele_type = 1
        if( present( cv_sele_type ) ) cv_sele_type = 1

        return
    end subroutine Get_Ele_Type

    subroutine Get_Ele_Type_new( Mdims, Mdisopt )
        !-
        !- Mdisopt%u_ele_type = Mdisopt%cv_ele_type = Mdisopt%p_ele_type will flag the dimension and
        !- type of element:
        !- = 1 or 2: 1D (linear and quadratic, respectively)
        !- = 3 or 4: triangle (linear or quadratic, respectively)
        !- = 5 or 6: quadrilateral (bi-linear or tri-linear, respectively)
        !- = 7 or 8: tetrahedron (linear or quadratic, respectively)
        !- = 9 or 10: hexahedron (bi-linear or tri-linear, respectively)
        !-
        implicit none
        type(multi_dimensions), intent(in) :: Mdims
        type (multi_discretization_opts) :: Mdisopt
        !!$ Local variables
        integer :: ndim, degree
        call get_option( '/geometry/dimension', ndim)
        call get_option( &
            '/geometry/mesh::PressureMesh/from_mesh/mesh_shape/polynomial_degree', &
            degree )
        Select Case( ndim )
            case( 1 ) ! ndim
                Select Case( degree )
                    case( 1 ) ! degree
                        !!$ ndim=1; p=1
                        Mdisopt%cv_ele_type = 1
                    case( 2 ) ! degree
                        ! ndim=1; p=2
                        Mdisopt%cv_ele_type = 2
                    case default; FLAbort('Degree error')
                end Select ! degree
            case( 2 ) ! ndim
                Select Case( degree )
                    case( 1 ) ! degree
                        Select Case( Mdims%x_nloc )
                            case( 3 ) ! Mdims%x_nloc
                                ! ndim=2; p=1; Mdims%x_nloc=3
                                ! linear triangle
                                Mdisopt%cv_ele_type = 3
                            case( 4 ) ! Mdims%x_nloc
                                ! ndim=2; p=1; Mdims%x_nloc=4
                                ! bilinear quad
                                Mdisopt%cv_ele_type = 5
                            case default; FLAbort('Mdims%x_nloc error')
                        end Select ! Mdims%x_nloc
                    case( 2 ) ! degree
                        Select Case( Mdims%x_nloc )
                            case( 6 ) ! Mdims%x_nloc
                                ! ndim=2; p=2; Mdims%x_nloc=3
                                ! quadratic triangle
                                Mdisopt%cv_ele_type = 4
                            case( 10 ) ! Mdims%x_nloc
                                ! ndim=2; p=2; Mdims%x_nloc=4
                                ! bi-quadratic quad
                                Mdisopt%cv_ele_type = 6
                            case default; FLAbort('Mdims%x_nloc error')
                        end Select ! Mdims%x_nloc
                    case default; FLAbort('Degree error')
                end Select ! degree
            case( 3 ) ! ndim
                Select Case( degree )
                    case( 1 ) ! degree
                        Select Case( Mdims%x_nloc )
                            case( 4 ) ! Mdims%x_nloc
                                ! ndim=3; p=1; Mdims%x_nloc=4
                                ! linear tets
                                Mdisopt%cv_ele_type = 7
                            case( 8 ) ! Mdims%x_nloc
                                ! ndim=3; p=1; Mdims%x_nloc=8
                                ! tri-linear hex
                                Mdisopt%cv_ele_type = 9
                            case default; FLAbort('Mdims%x_nloc error')
                        end Select ! Mdims%x_nloc
                    case( 2 ) ! degree
                        Select Case( Mdims%x_nloc )
                            case( 10 ) ! Mdims%x_nloc
                                ! ndim=3; p=2; Mdims%x_nloc=4
                                ! quadratic tet
                                Mdisopt%cv_ele_type = 8
                            case( 27 ) ! Mdims%x_nloc
                                ! ndim=3; p=2; Mdims%x_nloc=8
                                ! bilinear quad
                                Mdisopt%cv_ele_type = 10
                            case default; FLAbort('Mdims%x_nloc error')
                        end Select ! Mdims%x_nloc
                    case default; FLAbort('Degree error')
                end Select ! degree
        end Select ! ndim
        Mdisopt%p_ele_type = Mdisopt%cv_ele_type ; Mdisopt%u_ele_type = Mdisopt%cv_ele_type
        !!$ The following options are hardcoded and need to be either deleted from the code tree or
        !!$ added into the schema.
        Mdisopt%mat_ele_type = 1
        Mdisopt%u_sele_type = 1
        Mdisopt%cv_sele_type = 1
        return
    end subroutine Get_Ele_Type_new

    subroutine Get_Discretisation_Options( state, Mdims, Mdisopt )
        !!$ This subroutine extract all discretisation options from the schema
        implicit none
        type( state_type ), dimension( : ), intent( in ) :: state
        type(multi_dimensions), intent(in) :: Mdims
        type (multi_discretization_opts) :: Mdisopt
        !!$ Local variables:
        integer :: iphase
        character( len = option_path_len ) :: option_path, option_path2, option_path3
        !!$ DISOPT Options:
        !!$ =0      1st order in space          Theta=specified    UNIVERSAL
        !!$ =1      1st order in space          Theta=non-linear   UNIVERSAL
        !!$ =2      Trapezoidal rule in space   Theta=specified    UNIVERSAL
        !!$ =2      if isotropic limiter then FEM-quadratic & stratification adjust. Theta=non-linear
        !!$ =3      Trapezoidal rule in space   Theta=non-linear   UNIVERSAL
        !!$ =4      Finite elements in space    Theta=specified    UNIVERSAL
        !!$ =5      Finite elements in space    Theta=non-linear   UNIVERSAL
        !!$ =6      Finite elements in space    Theta=specified    NONE
        !!$ =7      Finite elements in space    Theta=non-linear   NONE
        !!$ =8      Finite elements in space    Theta=specified    DOWNWIND+INTERFACE TRACKING
        !!$ =9      Finite elements in space    Theta=non-linear   DOWNWIND+INTERFACE TRACKING
        !!$ Solving Advection Field: Temperature
        option_path = '/material_phase[0]/scalar_field::Temperature'
        option_path2 = trim( option_path ) //  '/prognostic/spatial_discretisation'
        option_path3 = trim( option_path ) //  '/prognostic/temporal_discretisation/control_volumes/number_advection_iterations'
        Mdisopt%t_disopt = 1
        call get_option( trim( option_path3 ), Mdisopt%nits_flux_lim_t, default = 3 )
        Conditional_TDISOPT: if( have_option( trim( option_path2 ) ) ) then
            if( have_option( trim( option_path2 ) // '/control_volumes/face_value::FiniteElement/limit_face_value/' // &
                'limiter::CompressiveAdvection' ) ) then
                Mdisopt%t_disopt = 9
            else
                if( have_option( trim( option_path2 ) // '/control_volumes/face_value::FiniteElement/limit_face_value' ) ) &
                    Mdisopt%t_disopt = 5
            end if
        end if Conditional_TDISOPT
        call get_option( trim( option_path2 ) // '/conservative_advection', Mdisopt%t_beta, default = 0.0 )
        call get_option( '/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation/theta', &
            Mdisopt%t_theta, default = 1. )
        !!$ Solving Advection Field: Volume fraction
        option_path = '/material_phase[0]/scalar_field::PhaseVolumeFraction'
        option_path2 = trim( option_path ) // '/prognostic/spatial_discretisation/control_volumes/face_value'
        option_path3 = trim( option_path ) // '/prognostic/temporal_discretisation/control_volumes/number_advection_iterations'
        Mdisopt%v_disopt = 8
        call get_option( trim( option_path3 ), Mdisopt%nits_flux_lim_volfra, default = 3 )
        Conditional_VDISOPT: if( have_option( trim( option_path ) ) ) then
            if( have_option( trim( option_path2 ) // '::FirstOrderUpwind' ) ) Mdisopt%v_disopt = 0
            if( have_option( trim( option_path2 ) // '::Trapezoidal' ) ) Mdisopt%v_disopt = 2
            if( have_option( trim( option_path2 ) // '::FiniteElement/do_not_limit_face_value' ) ) Mdisopt%v_disopt = 6
            if( have_option( trim( option_path2 ) // '::FiniteElement/limit_face_value/limiter::Sweby' ) ) Mdisopt%v_disopt = 5
            if( have_option( trim( option_path2 ) // '::FiniteElement/limit_face_value/limiter::CompressiveAdvection' ) ) Mdisopt%v_disopt = 9
        end if Conditional_VDISOPT

        call get_option( trim( option_path ) // '/prognostic/spatial_discretisation/conservative_advection', Mdisopt%v_beta )
        call get_option( trim( option_path ) // '/prognostic/temporal_discretisation/theta', Mdisopt%v_theta )
        !!$ Solving Velocity Field
        call get_option( '/material_phase[0]/vector_field::Velocity/prognostic/temporal_discretisation/theta', Mdisopt%u_theta )
        !!$ Solving Component Field
        option_path3 = '/material_phase[' // int2str( Mdims%nphase ) // ']/scalar_field::ComponentMassFractionPhase1/' // &
            'temporal_discretisation/control_volumes/number_advection_iterations'
        call get_option( trim( option_path3 ), Mdisopt%nits_flux_lim_comp, default = 3 )
        !!$ Scaling factor for the momentum equation
        Mdisopt%scale_momentum_by_volume_fraction = .false.
        do iphase = 1, Mdims%nphase
            option_path = '/material_phase[' // int2str( iphase - 1 ) // ']/Mdisopt%scale_momentum_by_volume_fraction'
            if( have_option( trim( option_path ) ) ) Mdisopt%scale_momentum_by_volume_fraction = .true.
        end do
        !!$ Options below are hardcoded and need to be added into the schema
        Mdisopt%t_dg_vel_int_opt = 1 ; Mdisopt%u_dg_vel_int_opt = 4 ; Mdisopt%v_dg_vel_int_opt = 4 ; Mdisopt%w_dg_vel_int_opt = 0
        if(.not.is_porous_media) Mdisopt%v_dg_vel_int_opt = 1
        Mdisopt%volfra_use_theta_flux = .false. ; Mdisopt%volfra_get_theta_flux = .true.
        Mdisopt%comp_use_theta_flux = .false. ; Mdisopt%comp_get_theta_flux = .true.
        Mdisopt%t_use_theta_flux = .false. ; Mdisopt%t_get_theta_flux = .false.
        !!$ IN/Mdisopt%dg_ele_upwind are options for optimisation of upwinding across faces in the compact_overlapping
        !!$ formulation. The data structure and options for this formulation need to be added later.
        Mdisopt%in_ele_upwind = 3 ; Mdisopt%dg_ele_upwind = 3
        if (is_porous_media) Mdisopt%in_ele_upwind = Mdisopt%v_disopt



        return
    end subroutine Get_Discretisation_Options



    subroutine update_boundary_conditions( state, stotel, cv_snloc, nphase, &
        &                                 suf_t_bc, suf_t_bc_rob1, suf_t_bc_rob2, tracer )
        implicit none
        type( state_type ), dimension( : ), intent( in ) :: state
        integer, intent( in ) :: stotel, cv_snloc, nphase
        real, dimension( 1, nphase, stotel * cv_snloc ), intent( inout ) :: suf_t_bc, suf_t_bc_rob1, suf_t_bc_rob2
        !
        character( len = option_path_len ) :: option_path, option_path2, field_name, name
        integer :: shape_option(2), iphase, nobcs, kk, k, j , sele, stat
        integer, dimension( : ), allocatable :: SufID_BC
        integer, dimension( : ), pointer:: surface_element_list, face_nodes
        type(tensor_field), intent(inout), target :: tracer

        type( scalar_field ), pointer :: field, field_prot_bc, field_prot1, field_prot2
        type( scalar_field ), pointer :: field_prot_bc1, field_prot_bc2
        type( scalar_field ) :: field_prot_bc1f, field_prot_bc2f
        type( mesh_type ), pointer :: pmesh, surface_mesh

        suf_t_bc = 0. ; suf_t_bc_rob1 = 0. ; suf_t_bc_rob2 = 0.

        call set_boundary_conditions_values( state, shift_time = .true. )

        pmesh => extract_mesh( state, 'PressureMesh' )

        do iphase = 1, nphase

            field_name = 'Temperature'
            field => extract_scalar_field( state( iphase ), trim( field_name ) )

            option_path = '/material_phase['//int2str( iphase - 1 )//']/scalar_field::'//trim( field_name )

            option_path2 = trim( option_path ) // '/prognostic/boundary_conditions['

            nobcs = get_boundary_condition_count( field )

            Loop_BC: do k = 1, nobcs

                option_path = trim( option_path2 ) // int2str( k - 1 ) // ']/surface_ids'
                shape_option = option_shape( trim( option_path ) )
                allocate( SufID_BC( 1 : shape_option( 1 ) ) )
                call get_option( trim( option_path ), SufID_BC )

                option_path = trim( option_path2 ) // int2str( k - 1 ) // ']/'

                Conditional_Field_BC: if( have_option( trim( option_path ) // 'type::dirichlet' ) ) then

                    !BC_Type = 1
                    field_prot_bc => extract_surface_field( field, k, 'value' )

                    sele = 1
                    do j = 1, stotel
                        if( any ( SufID_BC == pmesh % faces % boundary_ids( j ) ) ) then
                            !wic_bc( j + ( iphase - 1 ) * stotel ) = BC_Type
                            face_nodes => ele_nodes( field_prot_bc, sele )
                            do kk = 1, cv_snloc
                                suf_t_bc( 1, iphase, ( j - 1 ) * cv_snloc + kk ) = &
                                    field_prot_bc % val( face_nodes( 1 ) )
                            end do
                            sele = sele + 1
                        end if
                    end do

                else if( have_option( trim( option_path ) // 'type::robin' ) ) then

                    !BC_Type = 2
                    if( have_option( trim( option_path ) // 'type::robin/order_zero_coefficient/from_field' ) ) then

                        call get_boundary_condition( field, k, surface_mesh = surface_mesh, &
                            surface_element_list = surface_element_list )

                        call allocate( field_prot_bc1f, surface_mesh, "Robin1" )
                        call allocate( field_prot_bc2f, surface_mesh, "Robin2" )

                        call get_option( trim( option_path ) // "type::robin/order_zero_coefficient/from_field/name", name )
                        field_prot1 => extract_scalar_field( state( iphase ), name, stat )
                        if(stat /= 0) FLExit( "Could not extract parent field 1. Check options file?" )

                        call get_option( trim( option_path ) // "type::robin/order_one_coefficient/from_field/name", name )
                        field_prot2 => extract_scalar_field( state( iphase ), name, stat )
                        if(stat /= 0) FLExit( "Could not extract parent field 2. Check options file?" )

                        call remap_field_to_surface( field_prot1, field_prot_bc1f, surface_element_list )
                        call remap_field_to_surface( field_prot2, field_prot_bc2f, surface_element_list )

                        ! copy back memory
                        sele = 1
                        do j = 1, stotel
                            if( any ( SufID_BC == pmesh % faces % boundary_ids( j ) ) ) then
                                !wic_bc( j + ( iphase - 1 ) * stotel ) = BC_Type
                                face_nodes => ele_nodes( field_prot_bc1f, sele )
                                do kk = 1, cv_snloc
                                    suf_t_bc_rob1( 1, iphase, ( j - 1 ) * cv_snloc + kk ) = &
                                        field_prot_bc1f % val( face_nodes( 1 ) )
                                    suf_t_bc_rob2( 1, iphase, ( j - 1 ) * cv_snloc + kk ) = &
                                        field_prot_bc2f % val( face_nodes( 1 ) )
                                end do
                                sele = sele + 1
                            end if
                        end do

                        call deallocate( field_prot_bc1f )
                        call deallocate( field_prot_bc2f )

                    else

                        field_prot_bc1 => extract_surface_field( field, k, 'order_zero_coefficient' )
                        field_prot_bc2 => extract_surface_field( field, k, 'order_one_coefficient' )

                        sele = 1
                        do j = 1, stotel
                            if( any ( SufID_BC == pmesh % faces % boundary_ids( j ) ) ) then
                                !wic_bc( j + ( iphase - 1 ) * stotel ) = BC_Type
                                face_nodes => ele_nodes( field_prot_bc1, sele )
                                do kk = 1, cv_snloc
                                    suf_t_bc_rob1( 1, iphase, ( j - 1 ) * cv_snloc + kk ) = &
                                        field_prot_bc1 % val( face_nodes( 1 ) )
                                    suf_t_bc_rob2( 1, iphase, ( j - 1 ) * cv_snloc + kk ) = &
                                        field_prot_bc2 % val( face_nodes( 1 ) )
                                end do
                                sele = sele + 1
                            end if
                        end do

                    end if

                end if Conditional_Field_BC

                deallocate( SufID_BC )

            end do Loop_BC

        end do

    end subroutine update_boundary_conditions

    subroutine pack_multistate(npres, state, packed_state, &
        multiphase_state, multicomponent_state, pmulti_state)

        integer, intent(in) :: npres
        type(state_type), dimension(:), intent(inout) :: state
        type(state_type), intent(inout) :: packed_state
        type(state_type), dimension(:), intent(inout), pointer :: &
            multiphase_state, multicomponent_state
        type(state_type), dimension(:,:), pointer, optional :: pmulti_state

        type(state_type), dimension(:,:), pointer :: multi_state
        type(scalar_field), pointer :: pressure, sfield
        type(vector_field), pointer :: velocity, position, vfield
        type(tensor_field), pointer :: tfield, p2, d2, drhodp
        type(vector_field) :: porosity, vec_field, porous_density, porous_heat_capacity
        type(vector_field) :: p_position, u_position, m_position
        type(tensor_field) :: permeability, ten_field, porous_thermal_conductivity
        type(mesh_type), pointer :: ovmesh, element_mesh
        type(element_type) :: element_shape
        integer, dimension( : ), pointer :: element_nodes
        logical :: has_density, has_phase_volume_fraction
        integer :: i, iphase, icomp, idim, iele, ipres
        integer :: nphase,ncomp,ndim,stat,n_in_pres

#ifdef USING_FEMDEM
        if(have_option('/simulation_type/femdem_fracture')) then
            if(have_option('/simulation_type/femdem_fracture/oneway_coupling_only')) then!This option do not exist
                sfield=>extract_scalar_field(state(1),"SolidConcentration")
                call insert(packed_state,sfield,"SolidConcentration")
                call add_new_memory(packed_state,sfield,"OldSolidConcentration")

                tfield=>extract_tensor_field(state(1),"Viscosity")
                call insert(packed_state,tfield,"Viscosity")

                sfield=>extract_scalar_field(state(1),"Dummy")
                call insert(packed_state,sfield,"Dummy")

                vfield=>extract_vector_field(state(1),"Darcy_Velocity")
                call insert(packed_state,vfield,"Darcy_Velocity")

                sfield=>extract_scalar_field(state(1),"TotalFlux")
                call insert(packed_state,sfield,"TotalFlux")
            else
                sfield=>extract_scalar_field(state(1),"SolidConcentration")
                call insert(packed_state,sfield,"SolidConcentration")
                call add_new_memory(packed_state,sfield,"OldSolidConcentration")

                tfield=>extract_tensor_field(state(1),"Viscosity")
                call insert(packed_state,tfield,"Viscosity")

                sfield=>extract_scalar_field(state(1),"Dummy")
                call insert(packed_state, sfield,"Dummy" )

                sfield=>extract_scalar_field(state(1),"TotalFlux")
                call insert(packed_state,sfield,"TotalFlux")

                vfield=>extract_vector_field(state(1),"Darcy_Velocity")
                call insert(packed_state,vfield,"Darcy_Velocity" )

                vfield=>extract_vector_field(state(1),"delta_U")
                call insert(packed_state,vfield,"delta_U")

                vfield=>extract_vector_field(state(1),"solid_U")
                call insert(packed_state,vfield,"solid_U")
            end if
        end if
#endif

        ncomp=option_count('/material_phase/is_multiphase_component')
        nphase=size(state)-ncomp
        allocate(multiphase_state(nphase))
        allocate(multicomponent_state(ncomp))
        allocate(multi_state(max(1,ncomp),nphase))

        position=>extract_vector_field(state(1),"Coordinate")
        call insert(packed_state,position,"Coordinate")
        ndim=mesh_dim(position)
        call insert(packed_state,position%mesh,"CoordinateMesh")

        if(has_scalar_field(state(1),"Porosity")) then
            sfield=>extract_scalar_field(state(1),"Porosity")
            element_mesh=>sfield%mesh
            call insert(packed_state,element_mesh,'P0DG')
        else
            element_shape=make_element_shape(position%mesh%shape,degree=0)
            allocate(element_mesh)
            element_mesh=make_mesh(position%mesh,element_shape,&
                continuity=-1,name="ElementMesh")
            call insert(packed_state,element_mesh,'P0DG')
            call deallocate(element_mesh)
            deallocate(element_mesh)
            call deallocate(element_shape)
            element_mesh=>extract_mesh(packed_state,'P0DG')
        end if

        ! pack rock-fluid properties
        ! if there is capillary pressure, we store 5 entries, otherwise just 3:
        ! (Immobile fraction, Krmax, relperm exponent, [capillary entry pressure, capillary exponent])
        if(have_option_for_any_phase('/multiphase_properties/capillary_pressure',nphase)) then
            call allocate(ten_field,element_mesh,"PackedRockFluidProp",dim=[6,nphase])
        else
            call allocate(ten_field,element_mesh,"PackedRockFluidProp",dim=[3,nphase])
        end if
        call insert(packed_state,ten_field,"PackedRockFluidProp")
        call deallocate(ten_field)

        ! For Flooding: Manning coefficient
        if(have_option('/flooding')) then
            call allocate(ten_field,element_mesh,"PackedManningcoef",dim=[1,nphase])
            call insert(packed_state,ten_field,"PackedManningcoef")
            call deallocate(ten_field)
        end if

        pressure=>extract_scalar_field(state(1),"Pressure")
        call insert(packed_state,pressure%mesh,"PressureMesh")

        ! control volume barycentres
        call allocate(vec_field,ndim,pressure%mesh,"CVBarycentre")
        call zero(vec_field)
        call insert(packed_state,vec_field,"CVBarycentre")
        do icomp = 1, ncomp
            call insert(multicomponent_state(icomp),vec_field,"CVBarycentre")
        end do
        call deallocate(vec_field)

        call allocate(vec_field,1,pressure%mesh,"CVIntegral")
        call zero(vec_field)
        call insert(packed_state,vec_field,"CVIntegral")
        do icomp = 1, ncomp
            call insert(multicomponent_state(icomp),vec_field,"CVIntegral")
        end do
        call deallocate(vec_field)

        call allocate(vec_field,npres,pressure%mesh,"MeanPoreCV")
        call zero(vec_field)
        call insert(packed_state,vec_field,"MeanPoreCV")
        do icomp = 1, ncomp
            call insert(multicomponent_state(icomp),vec_field,"MeanPoreCV")
        end do
        call deallocate(vec_field)

        call insert_sfield(packed_state,"FEPressure",1,npres)
        tfield=>extract_tensor_field(packed_state,"PackedFEPressure")
        tfield%option_path=pressure%option_path
        p2=>extract_tensor_field(packed_state,"PackedFEPressure")
        do ipres = 1, npres
            p2%val(1,ipres,:)=pressure%val
        end do
        do icomp = 1, ncomp
            call insert(multicomponent_state(icomp),p2,"PackedFEPressure")
        end do

        call insert_sfield(packed_state,"CVPressure",1,npres)
        p2=>extract_tensor_field(packed_state,"PackedCVPressure")
        do ipres = 1, npres
            p2%val(1,ipres,:)=pressure%val
        end do

        ! dummy field on the pressure mesh, used for evaluating python eos's.
        ! (this could be cleaned up in the future)
        call add_new_memory(packed_state,pressure,"Dummy")

        call insert_sfield(packed_state,"FEDensity",1,nphase)
        d2=>extract_tensor_field(packed_state,"PackedFEDensity")
        do icomp = 1, ncomp
            call insert(multicomponent_state(icomp),d2,"PackedFEDensity")
        end do

        call insert_sfield(packed_state,"Density",1,nphase,&
            add_source=.false.)
        call insert_sfield(packed_state,"DensityHeatCapacity",1,nphase)

        call insert_sfield(packed_state,"DRhoDPressure",1,nphase)
        drhodp=>extract_tensor_field(packed_state,"PackedDRhoDPressure")
        do icomp = 1, ncomp
           call insert(multicomponent_state(icomp),drhodp,"PackedDRhoDPressure")
        end do

        if (option_count("/material_phase/scalar_field::Temperature")>0) then
            call insert_sfield(packed_state,"Temperature",1,nphase,&
                add_source=.true.,add_absorption=.true.)
            call insert_sfield(packed_state,"FETemperature",1,nphase)
        end if

        if (option_count("/material_phase/scalar_field::Bathymetry")>0) then
            call insert_sfield(packed_state,"Bathymetry",1,nphase,&
                add_source=.false.,add_absorption=.false.)
            call insert_sfield(packed_state,"Bathymetry",1,nphase)
        end if

        call insert_sfield(packed_state,"PhaseVolumeFraction",1,nphase,&
            add_source=.true.)
        call insert_sfield(packed_state,"FEPhaseVolumeFraction",1,nphase)
        if (ncomp>1) then
           call insert_sfield(packed_state,"PhaseVolumeFractionComponentSource",1,nphase)
        end if

        if( have_option_for_any_phase( '/multiphase_properties/capillary_pressure', nphase ) ) then
            call allocate(ten_field,pressure%mesh,"PackedCapPressure",dim=[1,nphase])
            call insert(packed_state,ten_field,"PackedCapPressure")
            call deallocate(ten_field)
        end if

        ! pack surface tension terms
        if( have_option_for_any_phase( '/is_multiphase_component/surface_tension', nphase+ncomp ) ) then
            call allocate(ten_field,pressure%mesh,"SurfaceTensionGrad",dim=[ncomp,nphase])
            call insert(packed_state,ten_field,"SurfaceTensionGrad")
            call deallocate(ten_field)
            call allocate(ten_field,pressure%mesh,"SurfaceTensionCoef",dim=[ncomp,nphase])
            call insert(packed_state,ten_field,"SurfaceTensionCoef")
            call deallocate(ten_field)
        end if

        ! pack continuous velocity mesh
        velocity=>extract_vector_field(state(1),"Velocity")
        call insert(packed_state,velocity%mesh,"VelocityMesh")
        if (.not.is_P0DGP1CV) then
            if (.not.has_mesh(state(1),"VelocityMesh_Continuous")) then
                nullify(ovmesh)
                allocate(ovmesh)
                ovmesh=make_mesh(position%mesh,&
                    shape=velocity%mesh%shape,&
                    continuity=0,name="VelocityMesh_Continuous")
                call insert(packed_state,ovmesh,"VelocityMesh_Continuous")
                call insert(state(1),ovmesh,"VelocityMesh_Continuous")
                call deallocate(ovmesh)
                deallocate(ovmesh)
            else
                ovmesh=>extract_mesh(state(1),"VelocityMesh_Continuous")
                call insert(packed_state,ovmesh,"VelocityMesh_Continuous")
            end if
            call allocate(u_position,ndim,ovmesh,"VelocityCoordinate")
            call remap_field(position,u_position)
            call insert(packed_state,u_position,"VelocityCoordinate")
            call deallocate(u_position)
        end if

        if (.not.has_mesh(state(1),"PressureMesh_Continuous")) then
            nullify(ovmesh)
            allocate(ovmesh)
            ovmesh=make_mesh(position%mesh,&
                shape=pressure%mesh%shape,&
                continuity=0,name="PressureMesh_Continuous")
            call insert(packed_state,ovmesh,"PressureMesh_Continuous")
            call deallocate(ovmesh)
            deallocate(ovmesh)
            ovmesh=>extract_mesh(packed_state,"PressureMesh_Continuous")
            call insert(state(1),ovmesh,"PressureMesh_Continuous")
        else
            ovmesh=>extract_mesh(state(1),"PressureMesh_Continuous")
            call insert(packed_state,ovmesh,"PressureMesh_Continuous")
        end if

        call allocate(p_position,ndim,ovmesh,"PressureCoordinate")
        call remap_field(position,p_position)
        call insert(packed_state,p_position,"PressureCoordinate")
        call deallocate(p_position)

        if (.not.has_mesh(state(1),"PressureMesh_Discontinuous")) then
            nullify(ovmesh)
            allocate(ovmesh)
            ovmesh=make_mesh(position%mesh,&
                shape=pressure%mesh%shape,&
                continuity=-1,name="PressureMesh_Discontinuous")
            call insert(packed_state,ovmesh,"PressureMesh_Discontinuous")
            call deallocate(ovmesh)
            deallocate(ovmesh)
            ovmesh=>extract_mesh(packed_state,"PressureMesh_Discontinuous")
            call insert(state(1),ovmesh,"PressureMesh_Discontinuous")
        else
            ovmesh=>extract_mesh(state(1),"PressureMesh_Discontinuous")
            call insert(packed_state,ovmesh,"PressureMesh_Discontinuous")
        end if

        call allocate(m_position,ndim,ovmesh,"MaterialCoordinate")
        call remap_field(position,m_position)
        call insert(packed_state,m_position,"MaterialCoordinate")
        call deallocate(m_position)

        call insert(packed_state,velocity%mesh,"InternalVelocityMesh")
        call insert_vfield(packed_state,"Velocity",add_source=.true.)
        call insert_vfield(packed_state,"NonlinearVelocity",zerod=.true.)
        call insert(state(1),velocity%mesh,"InternalVelocityMesh")
        call unpack_multiphase(packed_state,multiphase_state)
        if (ncomp>0) then
            call insert_sfield(packed_state,"ComponentDensity",ncomp,nphase)
            call insert_sfield(packed_state,"ComponentMassFraction",ncomp,&
                nphase,add_source=.true.)
            call insert_sfield(packed_state,"FEComponentDensity",ncomp,nphase)
            call insert_sfield(packed_state,"FEComponentMassFraction",ncomp,nphase)
        end if

        if (is_porous_media) then
            ovmesh=>extract_mesh(packed_state,"PressureMesh_Discontinuous")
            if ( ncomp > 0 ) then
                ovmesh=>extract_mesh(packed_state,"PressureMesh")
                do icomp = 1, ncomp
                    ! Add component absorption (nphase, nphase, cv_nonods)
                    ! to packed_state and all multicomponent_states
                    call allocate(ten_field,ovmesh,"ComponentAbsorption",dim=[nphase,nphase])
                    call insert(multicomponent_state(icomp),ten_field,"ComponentAbsorption")
                    call deallocate(ten_field)
                end do
            end if
        end if




        call allocate(porosity,npres,element_mesh,"Porosity")
        do ipres = 1, npres
            call set(porosity,ipres,1.0)
        end do
        call insert(packed_state,porosity,"Porosity")
        call deallocate(porosity)
        if (has_scalar_field(state(1),"Porosity")) then
            sfield=>extract_scalar_field(state(1),"Porosity")
!            porosity%val(1,:)=sfield%val
            call set(porosity,1,sfield)
        end if

        ! hack to define a lateral from diamond
        if(npres>1) then
            vfield=>extract_vector_field(packed_state,"Porosity")
            sfield=>extract_scalar_field(state(1),"Pipe")
            call assign_val(vfield%val(2,:),sfield%val)
        end if

        if(has_scalar_field(state(1),"Permeability")) then
            call allocate(permeability,element_mesh,"Permeability",&
                dim=[mesh_dim(position),mesh_dim(position)])
            call zero(permeability)
            sfield=>extract_scalar_field(state(1),"Permeability")
            do idim=1,mesh_dim(position)
                call set(permeability,idim,idim,sfield)
            end do
            call insert(packed_state,permeability,"Permeability")
            call deallocate(permeability)
        else if(has_vector_field(state(1),"Permeability")) then
            call allocate(permeability,element_mesh,"Permeability",&
                dim=[mesh_dim(position),mesh_dim(position)])
            call zero(permeability)
            vfield=>extract_vector_field(state(1),"Permeability")
            call set(permeability,vfield)
            call insert(packed_state,permeability,"Permeability")
            call deallocate(permeability)
        else if(has_tensor_field(state(1),"Permeability")) then
            call allocate(permeability,element_mesh,"Permeability",&
                dim=[mesh_dim(position),mesh_dim(position)])
            call zero(permeability)
            tfield=>extract_tensor_field(state(1),trim(permeability%name))
            if(size(tfield%val,3)==1) then!constant field
                do iele=1,element_count(tfield)
                    Permeability%val(:,:,iele)=tfield%val(:,:,1)
                end do
            else ! python
                do iele=1,element_count(tfield)
                    element_nodes=>ele_nodes(tfield,iele)
                    Permeability%val(:,:,iele)=tfield%val(:,:,element_nodes(1))
                end do
            end if
            call insert(packed_state,permeability,"Permeability")
            call deallocate(permeability)
        else
            call allocate(permeability,element_mesh,"Permeability",&
                dim=[mesh_dim(position),mesh_dim(position)])
            call zero(permeability)
            call insert(packed_state,permeability,"Permeability")
            call deallocate(permeability)
        end if

        !!$ pack multi_state information
        has_density=has_scalar_field(state(1),"Density")
        has_phase_volume_fraction=has_scalar_field(state(1),"PhaseVolumeFraction")
        iphase=1; icomp=1
        do i = 1, size(state)
            if(have_option(trim(state(i)%option_path)&
                //'/is_multiphase_component')) then
                velocity=>extract_vector_field(state(i),"Velocity",stat)
                if(stat==0) velocity%wrapped=.true.
                velocity=>extract_vector_field(state(i),"OldVelocity",stat)
                if(stat==0) velocity%wrapped=.true.
                velocity=>extract_vector_field(state(i),"NonlinearVelocity",stat)
                if(stat==0) velocity%wrapped=.true.
                velocity=>extract_vector_field(state(i),"IteratedVelocity",stat)
                if(stat==0) velocity%wrapped=.true.

                call unpack_component_sfield(state(i),packed_state,"FEComponentDensity",icomp)
                call unpack_component_sfield(state(i),packed_state,"OldFEComponentDensity",icomp,prefix='Old')
                call unpack_component_sfield(state(i),packed_state,"FEComponentMassFraction",icomp)
                call unpack_component_sfield(state(i),packed_state,"OldFEComponentMassFraction",icomp,prefix='Old')
                call unpack_component_sfield(state(i),packed_state,"OldComponentDensity",icomp,prefix='Old')
                call unpack_component_sfield(state(i),packed_state,"OldComponentMassFraction",icomp,prefix='Old')
                call unpack_component_sfield(state(i),packed_state,"ComponentDensity",icomp)
                call unpack_component_sfield(state(i),packed_state,"IteratedComponentMassFraction",icomp,prefix='Iterated')
                call unpack_component_sfield(state(i),packed_state,"ComponentMassFraction",icomp)
                call insert_components_in_multi_state(multi_state(icomp,:),state(i))

                icomp=icomp+1
                cycle
            else
                pressure=>extract_scalar_field(state(i),"Pressure",stat)
                if(stat==0) pressure%wrapped=.true.
                if(has_density) then
                    call unpack_sfield(state(i),packed_state,"IteratedDensity",1,iphase,&
                        check_paired(extract_scalar_field(state(i),"Density"),&
                        extract_scalar_field(state(i),"IteratedDensity")))
                    call unpack_sfield(state(i),packed_state,"OldDensity",1,iphase,&
                        check_paired(extract_scalar_field(state(i),"Density"),&
                        extract_scalar_field(state(i),"OldDensity")))
                    call unpack_sfield(state(i),packed_state,"Density",1,iphase)
!                    call unpack_sfield(state(i),packed_state,"DensitySource",1,iphase, free=.false.)
                    call insert(multi_state(1,iphase),extract_scalar_field(state(i),"Density"),"Density")
                end if

                if(have_option(trim(state(i)%option_path)&
                    //'/scalar_field::Temperature')) then
                    call unpack_sfield(state(i),packed_state,"OldTemperature",1,iphase,&
                        check_paired(extract_scalar_field(state(i),"Temperature"),&
                        extract_scalar_field(state(i),"OldTemperature")))
                    call unpack_sfield(state(i),packed_state,"IteratedTemperature",1,iphase,&
                        check_paired(extract_scalar_field(state(i),"Temperature"),&
                        extract_scalar_field(state(i),"IteratedTemperature")))
                    call unpack_sfield(state(i),packed_state,"TemperatureSource",1,iphase)
                    call unpack_sfield(state(i),packed_state,"TemperatureAbsorption",1,iphase)
                    call unpack_sfield(state(i),packed_state,"Temperature",1,iphase)
                    call insert(multi_state(1,iphase),extract_scalar_field(state(i),"Temperature"),"Temperature")
                end if

                if(has_phase_volume_fraction) then
                    call unpack_sfield(state(i),packed_state,"IteratedPhaseVolumeFraction",1,iphase,&
                        check_paired(extract_scalar_field(state(i),"IteratedPhaseVolumeFraction"),&
                        extract_scalar_field(state(i),"PhaseVolumeFraction")))
                    call unpack_sfield(state(i),packed_state,"OldPhaseVolumeFraction",1,iphase,&
                        check_paired(extract_scalar_field(state(i),"PhaseVolumeFraction"),&
                        extract_scalar_field(state(i),"OldPhaseVolumeFraction")))
                    call unpack_sfield(state(i),packed_state,"PhaseVolumeFraction",1,iphase)
                    call unpack_sfield(state(i),packed_state,"PhaseVolumeFractionSource",1,iphase)
                    call insert(multi_state(1,iphase),extract_scalar_field(state(i),"PhaseVolumeFraction"),"PhaseVolumeFraction")
                end if

                call unpack_vfield(state(i),packed_state,"IteratedVelocity",iphase,&
                    check_vpaired(extract_vector_field(state(i),"Velocity"),&
                    extract_vector_field(state(i),"IteratedVelocity")))
                call unpack_vfield(state(i),packed_state,"OldVelocity",iphase,&
                    check_vpaired(extract_vector_field(state(i),"Velocity"),&
                    extract_vector_field(state(i),"OldVelocity")))
                call unpack_vfield(state(i),packed_state,"NonlinearVelocity",iphase,&
                    check_vpaired(extract_vector_field(state(i),"Velocity"),&
                    extract_vector_field(state(i),"NonlinearVelocity")))
                call unpack_vfield(state(i),packed_state,"Velocity",iphase)
                call insert(multi_state(1,iphase),extract_vector_field(state(i),"Velocity"),"Velocity")

                iphase=iphase+1
            end if
        end do

        n_in_pres=nphase/npres
        do ipres = 1, npres
            i=(ipres-1)*n_in_pres+1
            call unpack_sfield(state(i),packed_state,"Pressure",1,ipres)
            call insert(multi_state(1,ipres),extract_scalar_field(state(i),"Pressure"),"FEPressure")
        end do

        if (option_count("/material_phase/scalar_field::Temperature")>0) then
            call allocate_multiphase_scalar_bcs(packed_state,multi_state,"Temperature")
        end if

        call allocate_multiphase_scalar_bcs(packed_state,multi_state,"Density")
        call allocate_multiphase_scalar_bcs(packed_state,multi_state,"PhaseVolumeFraction")
        call allocate_multiphase_vector_bcs(packed_state,multi_state,"Velocity")
        call allocate_multiphase_scalar_bcs(packed_state,multi_state,"FEPressure")

        if (ncomp>0) then
            call unpack_multicomponent(packed_state,multicomponent_state)
            call allocate_multicomponent_scalar_bcs(multicomponent_state,multi_state,"ComponentMassFraction")
            call allocate_multicomponent_scalar_bcs(multicomponent_state,multi_state,"ComponentDensity")
        end if

        if (present(pmulti_state)) then
            pmulti_state=>multi_state
        else
            call deallocate(multi_state)
            deallocate(multi_state)
        end if

        !!$ memory allocation for darcy velocity
        if(is_porous_media) then
            if(have_option('/io/output_darcy_vel') .or. is_multifracture) then
                ! allocate darcy velocity[in packed_state]
                call allocate(ten_field,velocity%mesh,"PackedDarcyVelocity", dim=[ndim,nphase])
                call zero(ten_field)
                call insert(packed_state,ten_field,"PackedDarcyVelocity")
                call deallocate(ten_field)

                do iphase = 1, n_in_pres
                    call unpack_vfield(state(iphase),packed_state,"DarcyVelocity",iphase)
                end do

!                ! let velocity[in state] point to darcy velocity[packed_state]
!                tfield=>extract_tensor_field(packed_state,"PackedDarcyVelocity")
!
!                do iphase=1,size(state)
!                    vfield=>extract_vector_field(state(iphase),"Velocity")
!                    vfield%val=>tfield%val(:,iphase,:)
!                end do

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! unfinished part -- still working on the field copy for when adaptivity is switched on
                ! 3 places to be changed: here, Populate_State, Multiphase_TimeLoop
                !
                ! add velocity_int[in state] and point it to internal velocity[packed_state]
                ! tfield=>extract_tensor_field(packed_state,"PackedVelocity")
                !do iphase=1,size(state)
                    ! call allocate_and_insert_vector_field('/material_phase['//int2str(size(state)-1)//']/vector_field::Velocity', &
                    !   state(iphase), field_name = "Velocity_int", parent_mesh = "VelocityMesh", dont_assign_boundary_condition = .true.)
                    !vfield=>extract_vector_field(state(iphase),"Velocity_bak")
                    !vfield%val=>tfield%val(:,iphase,:)
                !end do
            end if
        end if


        !! // How to add density as a dependant of olddensity,  as an example
        !! // Suppose we have a previous declaration
        !! type(tensor_field), pointer :: old_density, density
        !! old_density=>extract_tensor_field(packed_state,"PackedOldDensity")
        !! density=>extract_tensor_field(packed_state,"PackedDensity")
        !! call add_dependant_field(old_density,density)
        !! // if we now make a  call
        !! call mark_as_updated(density)
        !! // then is_updated(density) is now .true. (no other changes)
        !! // If we then make a  call
        !! call mark_as_updated(old_density)
        !! // then is_updated(old_density) is now .true. and is_updated(density) is .false.

    contains

        subroutine insert_components_in_multi_state(ms,s)
            type(state_type), dimension(:) :: ms
            type(state_type) :: s

            integer :: i
            character(len=FIELD_NAME_LEN) :: full_name, short_name

            do i=1,size(ms)
                full_name="ComponentMassFractionPhase"//int2str(i)
                short_name="ComponentMassFraction"
                if (has_scalar_field(s,trim(full_name))) then
                    call insert(ms(i),extract_scalar_field(s,trim(full_name)),short_name)
                end if
            end do

            do i=1,size(ms)
                full_name="ComponentDensityPhase"//int2str(i)
                short_name="ComponentDensity"
                if (has_scalar_field(s,trim(full_name))) then
                    call insert(ms(i),extract_scalar_field(s,trim(full_name)),short_name)
                end if
            end do

        end subroutine insert_components_in_multi_state

        subroutine unpack_multicomponent(mstate,mcstate)
            type(state_type) :: mstate
            type(state_type), dimension(:) :: mcstate

            integer :: icomp

            type(tensor_field), pointer :: tfield
            type(vector_field), pointer :: vecfield
            type(tensor_field) :: vfield

            character (len=FIELD_NAME_LEN), dimension(5) ::names
            character (len=FIELD_NAME_LEN), dimension(4) ::phase_names
            character (len=FIELD_NAME_LEN), dimension(1) ::phase_vector_names
            integer :: count

            names(1)="PackedComponentMassFraction"
            names(2)="PackedComponentDensity"
            names(3)="PackedOldComponentMassFraction"
            names(4)="PackedOldComponentDensity"
            names(5)="PackedComponentMassFractionSource"

            phase_names(1)="PackedPhaseVolumeFraction"
            phase_names(2)="PackedOldPhaseVolumeFraction"
            phase_names(3)="PackedNonlinearVelocity"
            phase_names(4)="PackedOldNonlinearVelocity"

            phase_vector_names(1)="PressureCoordinate"

            do count=1,size(names)

                tfield=>extract_tensor_field(mstate,names(count))

                do icomp=1,ncomp
                    call allocate(vfield,tfield%mesh,names(count),field_type=FIELD_TYPE_DEFERRED,dim=[1,nphase])
                    vfield%option_path=tfield%option_path
                    deallocate(vfield%val)
                    deallocate(vfield%bc)
                    vfield%val=>tfield%val(icomp:icomp,:,:)
                    vfield%wrapped=.true.
                    vfield%field_type=FIELD_TYPE_NORMAL
                    call insert(mcstate(icomp),vfield,vfield%name)
                    call deallocate(vfield)
                end do

            end do

            do count=1,size(phase_names)
                tfield=>extract_tensor_field(mstate,phase_names(count))
                do icomp=1,ncomp
                    call insert(mcstate(icomp),tfield,tfield%name)
                end do
            end do

            do count=1,size(phase_vector_names)
                vecfield=>extract_vector_field(mstate,phase_vector_names(count))
                do icomp=1,ncomp
                    call insert(mcstate(icomp),vecfield,vecfield%name)
                end do

            end do

        end subroutine unpack_multicomponent

        subroutine unpack_multiphase(mstate,mpstate)
            type(state_type) :: mstate
            type(state_type), dimension(:) :: mpstate
            integer :: index, iphase, si, s1, s2

            type(tensor_field), pointer :: tfield
            type(tensor_field) :: mp_tfield

            do index=1,size(mstate%tensor_fields)
                tfield=>extract_tensor_field(mstate,index)
                si=len(trim(tfield%name))
                !s1=max(0,si) ; s2=si
                if(si<=7)then
                   ! do nothing...
                else if(tfield%name(si-7:si)=="Pressure")then
                   ! do nothing...
                else if(tfield%name(:6)=="Packed")then
                    do iphase=1,nphase
                        call allocate(mp_tfield,tfield%mesh,tfield%name(7:),field_type=FiELD_TYPE_DEFERRED,dim=[tfield%dim(1),1])
                        mp_tfield%val=>tfield%val(:,iphase:iphase,:)
                        mp_tfield%updated=>tfield%updated
                        mp_tfield%wrapped=.true.
                        call insert(mpstate(iphase),mp_tfield,mp_tfield%name)
                        call deallocate(mp_tfield)
                    end do
                end if
            end do

        end subroutine unpack_multiphase

        subroutine add_new_memory(mstate,sfield,name)
            type(state_type) :: mstate
            type(scalar_field) :: sfield
            character (len=*) :: name
            type(scalar_field) :: sfield2
            integer :: i,j

            call allocate(sfield2,sfield%mesh,name)
            sfield2%option_path=sfield%option_path
            sfield2%bc=sfield%bc
            if (associated(sfield2%bc%boundary_condition)) then
                do i=1,size(sfield2%bc%boundary_condition)
                    if (associated(sfield2%bc%boundary_condition(i)&
                        &%surface_fields)) then
                        do j=1, size(sfield2%bc%boundary_condition(i)&
                            &%surface_fields)
                            call incref(sfield2%bc%boundary_condition(i)&
                                &%surface_fields(j))
                        end do
                    end if
                    call incref(sfield2%bc%boundary_condition(i)&
                        &%surface_mesh)
                end do
            end if
            call insert(mstate,sfield2,name)
            call deallocate(sfield2)

        end subroutine add_new_memory

        subroutine insert_sfield(mstate,name,ncomp,nphase,nmesh,&
            & add_source, add_absorption)
            type(state_type), intent(inout) :: mstate
            character(len=*), intent(in) :: name
            type(mesh_type),optional, target :: nmesh
            integer, intent(in) :: ncomp,nphase
            logical, optional, intent(in) :: add_source, add_absorption

            type(scalar_field), pointer :: nfield
            type(mesh_type), pointer :: lmesh
            type(tensor_field) :: mfield
            integer :: stat
            logical :: ladd_source, ladd_absorption

            if (present(add_source)) then
                ladd_source=add_source
            else
                ladd_source=.false.
            end if

            if (present(add_absorption)) then
                ladd_absorption=add_absorption
            else
                ladd_absorption=.false.
            end if

            nfield=>extract_scalar_field(state(1),name,stat)
            if (stat/=0) then
                nfield=>extract_scalar_field(state(1),"Pressure")
            end if

            if (present(nmesh)) then
                lmesh=> nmesh
            else
                lmesh=>nfield%mesh
            end if

            call allocate(mfield,lmesh,"Packed"//name,dim=[ncomp,nphase])
            call zero(mfield)
            call insert(mstate,mfield,"Packed"//name)
            call deallocate(mfield)
            call allocate(mfield,lmesh,"PackedOld"//name,dim=[ncomp,nphase])
            call zero(mfield)
            call insert(mstate,mfield,"PackedOld"//name)
            call deallocate(mfield)
            call allocate(mfield,lmesh,"PackedIterated"//name,dim=[ncomp,nphase])
            call zero(mfield)
            call insert(mstate,mfield,"PackedIterated"//name)
            call deallocate(mfield)

            if (ladd_source) then
                call allocate(mfield,lmesh,"Packed"//trim(name)//"Source",&
                    dim=[ncomp,nphase])
                call zero(mfield)
                call insert(mstate,mfield,"Packed"//trim(name)//"Source")
                call deallocate(mfield)
             end if

            if (ladd_absorption) then
                call allocate(mfield,lmesh,"Packed"//trim(name)//"Absorption",&
                    dim=[ncomp,nphase])
                call zero(mfield)
                call insert(mstate,mfield,"Packed"//trim(name)//"Absorption")
                call deallocate(mfield)
            end if

        end subroutine insert_sfield

        subroutine insert_vfield(mstate,name,nmesh,zerod, add_source, add_absorption)
            type(state_type), intent(inout) :: mstate
            character(len=*), intent(in) :: name
            type(mesh_type), optional, target :: nmesh
            logical, optional, intent(in) :: zerod, add_source, add_absorption

            type(vector_field), pointer :: nfield
            type(mesh_type), pointer :: lmesh
            type(tensor_field)  :: mfield
            logical :: lzero, ladd_source, ladd_absorption

            nfield=>extract_vector_field(state(1),name)
            if (present(nmesh)) then
                lmesh=> nmesh
            else
                lmesh=>nfield%mesh
            end if

            if (present(zerod)) then
                lzero=zerod
            else
                lzero=.false.
            end if

            if (present(add_source)) then
                ladd_source=add_source
            else
                ladd_source=.false.
            end if

            if (present(add_absorption)) then
                ladd_absorption=add_absorption
            else
                ladd_absorption=.false.
            end if

            call allocate(mfield,lmesh,"Packed"//name,dim=[ndim,nphase]) !!,contiguous=.true.)
            if (lzero) then
                call zero(mfield)
            end if

            call insert(mstate,mfield,"Packed"//name)
            call deallocate(mfield)
            call allocate(mfield,lmesh,"PackedOld"//name,dim=[ndim,nphase])
            if (lzero) then
                call zero(mfield)
            end if
            call insert(mstate,mfield,"PackedOld"//name)
            call deallocate(mfield)
            call allocate(mfield,lmesh,"PackedIterated"//name,dim=[ndim,nphase])
            if (lzero) then
                call zero(mfield)
            end if
            call insert(mstate,mfield,"PackedIterated"//name)
            call deallocate(mfield)

            if (ladd_source) then
                call allocate(mfield,lmesh,"Packed"//trim(name)//"Source",&
                    dim=[ndim,nphase])
                call zero(mfield)
                call insert(mstate,mfield,"Packed"//trim(name)//"Source")
                call deallocate(mfield)
            end if

            if (ladd_absorption) then
                call allocate(mfield,lmesh,"Packed"//trim(name)//"Absorption",&
                    dim=[ndim,nphase])
                call zero(mfield)
                call insert(mstate,mfield,"Packed"//trim(name)//"Absorption")
                call deallocate(mfield)
            end if
        end subroutine insert_vfield

        subroutine unpack_sfield(nstate,mstate,name,icomp, iphase,free)
            type(state_type), intent(inout) :: nstate, mstate
            character(len=*) :: name
            integer :: icomp,iphase, stat
            logical, optional :: free

            type(scalar_field), pointer :: nfield
            type(tensor_field), pointer :: mfield, python_tfield
            integer python_stat
            logical lfree

            if (present(free)) then
                lfree=free
            else
                lfree=.true.
            end if
            lfree=.false.

            if (trim(name)=="Pressure") then
                mfield=>extract_tensor_field(mstate,"PackedFE"//name)
                lfree=.true.
            else
                mfield=>extract_tensor_field(mstate,"Packed"//name)
            end if

            nfield=>extract_scalar_field(nstate,name,stat)

            if (stat==0) then
                if (size(nfield%val(:))>1) then
                    mfield%val(icomp,iphase,1:size(nfield%val))=nfield%val(:)
                else
                    mfield%val(icomp,iphase,:)=nfield%val(1)
                end if
                if (icomp==1 .and. iphase == 1) then
                    mfield%option_path=nfield%option_path
                end if
                if (lfree .and. associated(nfield%val)) then
#ifdef HAVE_MEMORY_STATS
                    call register_deallocation("scalar_field", "real", &
                        size(nfield%val), nfield%name)
#endif
                    deallocate(nfield%val)
                end if

                !sprint to_do: This flag makes the python scripting to work but breakes mesh adaptivity
                python_tfield => extract_tensor_field( state(1), "UAbsorB", python_stat )
                if (python_stat==0) then
                   if (trim(name)=="Pressure") then
                      nfield%val=>mfield%val(icomp,iphase,:)
                      nfield%val_stride=ncomp*nphase
                      nfield%wrapped=.true.
                   end if
                else
                   nfield%val=>mfield%val(icomp,iphase,:)
                   nfield%val_stride=ncomp*nphase
                   nfield%wrapped=.true.
                end if
            end if

          end subroutine unpack_sfield

        subroutine unpack_vfield(nstate,mstate,name,iphase,free)
            type(state_type), intent(inout) :: nstate, mstate
            character(len=*) :: name
            integer :: iphase, stat
            logical, optional :: free

            type(vector_field), pointer :: nfield
            type(tensor_field), pointer :: mfield
            logical lfree

            if (present(free)) then
                lfree=free
            else
                lfree=.true.
            end if

            nfield=>extract_vector_field(nstate,name,stat)
            mfield=>extract_tensor_field(mstate,"Packed"//name)

            if (stat==0) then
                mfield%val(:,iphase,:)=nfield%val(:,1:node_count(mfield))
                if (icomp==1 .and. iphase == 1) then
                    mfield%option_path=nfield%option_path
                end if
                if (lfree) then
#ifdef HAVE_MEMORY_STATS
                call register_deallocation("vector_field", "real", &
                    size(nfield%val), name=nfield%name)
#endif
                deallocate(nfield%val)
            end if
                nfield%val=>mfield%val(:,iphase,:)
                nfield%wrapped=.true.
            else
                call zero(mfield)
            end if

        end subroutine unpack_vfield

        subroutine unpack_component_sfield(st,mst,name,ic,prefix)
            type(state_type) :: st, mst
            character (len=*) :: name
            integer :: ic
            integer :: ip
            character (len=*), optional :: prefix

            type(scalar_field), pointer :: nfield,pnfield
            type(tensor_field), pointer :: mfield
            logical :: free

            mfield=>extract_tensor_field(mst,"Packed"//name)

            do ip=1,nphase
                if ( has_scalar_field(st,name//"Phase"//int2str(ip)) ) then

                    nfield=>extract_scalar_field(st,name//"Phase"//int2str(ip))
                    if (present(prefix)) then
                        pnfield=>extract_scalar_field(st,name(len(prefix)+1:)//"Phase"//int2str(ip))
                        free=check_paired(pnfield,nfield)
                    else
                        free=.true.
                    end if

                    mfield%val(ic,ip,:)=nfield%val(:)
                    if (ic==1 .and. ip==1) then
                        mfield%option_path=nfield%option_path
                    end if
                    if (free) then
#ifdef HAVE_MEMORY_STATS
                         call register_deallocation("scalar_field", "real", &
                              size(nfield%val), name=nfield%name)
#endif
                        deallocate(nfield%val)
                    end if
                    nfield%val=>mfield%val(ic,ip,:)
                    nfield%val_stride=ncomp*nphase
                    nfield%wrapped=.true.
                end if
            end do

        end subroutine unpack_component_sfield

        subroutine allocate_multiphase_scalar_bcs(s,ms,name)
            type(state_type), intent(inout) :: s
            type(state_type), dimension(:,:), intent(inout) :: ms
            character( len=*) :: name

            type(tensor_field), pointer :: mfield
            type(scalar_field), pointer :: sfield
            type(tensor_boundary_condition) :: tbc
            type(scalar_boundary_condition), pointer ::  bc
            type(tensor_boundary_condition), dimension(:), pointer :: temp

            integer :: iphase,icomp,stat, n, nbc, j

            nbc=0
            mfield=>extract_tensor_field(s,"Packed"//name)

            allocate(mfield%bc)

            do iphase=1,mfield%dim(2)
                do icomp=1,mfield%dim(1)
                    sfield=>extract_scalar_field( ms(icomp,iphase),trim(name),stat)
                    if (stat/= 0 ) cycle
                    do n=1,get_boundary_condition_count(sfield)
                        bc=>sfield%bc%boundary_condition(n)
                        nullify(tbc%applies)
                        allocate(tbc%applies(mfield%dim(1),mfield%dim(2)))
                        tbc%applies= .false.
                        tbc%applies(icomp,iphase)=.true.
                        tbc%name=bc%name
                        tbc%type=bc%type
                        tbc%surface_element_list=>bc%surface_element_list
                        tbc%surface_node_list=>bc%surface_node_list
                        tbc%surface_mesh=>bc%surface_mesh
                        call incref(tbc%surface_mesh)
                        tbc%scalar_surface_fields=>bc%surface_fields
                        do j=1,size(tbc%scalar_surface_fields)
                            call incref(tbc%scalar_surface_fields(j))
                        end do

                        nbc=nbc+1
                        if (nbc>1) then
                            temp=>mfield%bc%boundary_condition
                            allocate(mfield%bc%boundary_condition(nbc))
                            mfield%bc%boundary_condition(1:nbc-1)=temp
                            deallocate(temp)
                            mfield%bc%boundary_condition(nbc)=tbc
                        else
                            allocate(mfield%bc%boundary_condition(nbc))
                            mfield%bc%boundary_condition(nbc)=tbc
                        end if
                    end do
                end do
            end do

        end subroutine allocate_multiphase_scalar_bcs

        subroutine allocate_multiphase_vector_bcs(s,ms,name)
            type(state_type), intent(inout) :: s
            type(state_type), dimension(:,:), intent(inout) :: ms
            character( len=*) :: name

            type(tensor_field), pointer :: mfield
            type(vector_field), pointer :: vfield
            type(tensor_boundary_condition) :: tbc
            type(vector_boundary_condition), pointer ::  bc
            type(tensor_boundary_condition), dimension(:), pointer :: temp

            integer :: iphase, stat, n, nbc, j

            nbc=0
            mfield=>extract_tensor_field(s,"Packed"//name)

            allocate(mfield%bc)

            do iphase=1,mfield%dim(2)
                vfield=>extract_vector_field( ms(1,iphase),trim(name),stat)
                if (stat/= 0 ) cycle

                do n=1,get_boundary_condition_count(vfield)

                    bc=>vfield%bc%boundary_condition(n)

                    tbc%name=bc%name
                    tbc%type=bc%type
                    tbc%surface_element_list=>bc%surface_element_list
                    tbc%surface_node_list=>bc%surface_node_list
                    tbc%surface_mesh=>bc%surface_mesh
                    call incref(tbc%surface_mesh)
                    tbc%vector_surface_fields=>bc%surface_fields
                    do j=1,size(tbc%vector_surface_fields)
                        call incref(tbc%vector_surface_fields(j))
                    end do

                    nullify(tbc%applies)
                    allocate(tbc%applies(mfield%dim(1),mfield%dim(2)))
                    tbc%applies=.false.
                    tbc%applies(:,iphase)=bc%applies(1:mfield%dim(1))

                    nbc=nbc+1
                    if (nbc>1) then
                        temp=>mfield%bc%boundary_condition
                        allocate(mfield%bc%boundary_condition(nbc))
                        mfield%bc%boundary_condition(1:nbc-1)=temp
                        deallocate(temp)
                        mfield%bc%boundary_condition(nbc)=tbc
                    else
                        allocate(mfield%bc%boundary_condition(nbc))
                        mfield%bc%boundary_condition(nbc)=tbc
                    end if
                end do
            end do

        end subroutine allocate_multiphase_vector_bcs

        subroutine allocate_multicomponent_scalar_bcs(s,ms,name)
            type(state_type), dimension(:), intent(inout) :: s
            type(state_type), dimension(:,:), intent(inout) :: ms
            character( len=*) :: name

            type(tensor_field), pointer :: mfield
            type(scalar_field), pointer :: sfield
            type(tensor_boundary_condition) :: tbc
            type(scalar_boundary_condition), pointer ::  bc
            type(tensor_boundary_condition), dimension(:), pointer :: temp

            integer :: iphase,icomp,stat, n, nbc, j

            nullify(tbc%applies)

            do icomp=1,size(s)
                mfield=>extract_tensor_field(s(icomp),'Packed'//name)
                allocate(mfield%bc)
                nbc=0

                do iphase=1,mfield%dim(2)
                    sfield=>extract_scalar_field( ms(icomp,iphase),trim(name),stat)
                    if (stat/= 0 ) cycle
                    do n=1,get_boundary_condition_count(sfield)
                        bc=>sfield%bc%boundary_condition(n)
                        tbc%name=bc%name
                        tbc%type=bc%type
                        tbc%surface_element_list=>bc%surface_element_list
                        tbc%surface_node_list=>bc%surface_node_list
                        tbc%surface_mesh=>bc%surface_mesh
                        call incref(tbc%surface_mesh)
                        tbc%scalar_surface_fields=>bc%surface_fields
                        do j=1,size(tbc%scalar_surface_fields)
                            call incref(tbc%scalar_surface_fields(j))
                        end do
                        allocate(tbc%applies(1,mfield%dim(2)))
                        tbc%applies= .false.
                        tbc%applies(1,iphase)=.true.
                        nbc=nbc+1
                        if (nbc>1) then
                            temp=>mfield%bc%boundary_condition
                            allocate(mfield%bc%boundary_condition(nbc))
                            mfield%bc%boundary_condition(1:nbc-1)=temp
                            deallocate(temp)
                            mfield%bc%boundary_condition(nbc)=tbc
                        else
                            allocate(mfield%bc%boundary_condition(nbc))
                            mfield%bc%boundary_condition(nbc)=tbc
                        end if
                    end do
                end do
            end do

        end subroutine allocate_multicomponent_scalar_bcs

        function check_paired(sfield1,sfield2) result(unpaired)
            type(scalar_field) :: sfield1,sfield2
            logical unpaired

            unpaired = .not. associated(sfield2%val,sfield1%val)

        end function check_paired

        function check_vpaired(vfield1,vfield2) result(unpaired)
            type(vector_field) :: vfield1,vfield2
            logical unpaired

            unpaired = .not. associated(vfield2%val,vfield1%val)

        end function check_vpaired

    end subroutine pack_multistate

    subroutine prepare_absorptions(state, Mdims, multi_absorp)
        implicit none
        type(state_type), dimension(:), intent(inout) :: state
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_absorption), intent(inout) :: multi_absorp
        !Local variables
        integer :: k
        type(mesh_type), pointer :: ovmesh

        !Make sure that the memory is cleaned before using it
        call deallocate_multi_absorption(multi_absorp, .true.)
        !Prepare array that will contain the different absorptions
        ovmesh=>extract_mesh(state(1),"PressureMesh_Discontinuous")


        if (is_porous_media) then
             call allocate_multi_field( Mdims, multi_absorp%PorousMedia, ovmesh%nodes, field_name="PorousMedia_AbsorptionTerm")
!            if ( ncomp > 0 ) !"Not ready yet"
        end if
        !Need to add this
        if (is_flooding) then
             call allocate_multi_field( Mdims, multi_absorp%Flooding, ovmesh%nodes, field_name="Flooding_AbsorptionTerm")
        end if

    end subroutine prepare_absorptions



!    function wrap_as_tensor(field) result(tfield)
!
!      type(scalar_field), intent(inout) :: field
!      type(tensor_field), pointer :: tfield
!
!      allocate(tfield)
!      call allocate(tfield,field%mesh,name=field%name,dim=[1,1])
!
!      tfield%val(1,1,:)=field%val
!      deallocate(tfield%updated)
!      tfield%updated=field%updated
!      deallocate(field%val)
!      deallocate(field%updated)
!      field%val=>tfield%val(1,1,:)
!      field%updated=>tfield%updated
!      tfield%option_path=field%option_path
!
!    end function wrap_as_tensor

function as_vector(tfield,dim,slice) result(vfield)

    type(tensor_field), intent(inout) :: tfield
    integer, intent(in) :: dim
    integer, intent(in), optional :: slice

    type(vector_field)  :: vfield
    integer :: lslice

    if (present(slice)) then
        lslice=slice
    else
        lslice=1
    end if

    vfield%name=tfield%name
    vfield%mesh=tfield%mesh
    vfield%option_path=tfield%option_path
    vfield%dim=tfield%dim(dim)
    select case(dim)
        case(1)
            vfield%val=>tfield%val(:,lslice,:)
        case(2)
            vfield%val=>tfield%val(lslice,:,:)
    end select
    vfield%wrapped=.true.

end function as_vector

function as_packed_vector(tfield) result(vfield)

    type(tensor_field), intent(inout) :: tfield

    type(vector_field) :: vfield


    vfield%name=tfield%name
    vfield%mesh=tfield%mesh
    vfield%option_path=tfield%option_path
    vfield%dim=product(tfield%dim)

#ifdef USING_GFORTRAN
      vfield%val(1:vfield%dim,1:node_count(vfield)) => tfield%val !%contiguous_val
#else
    allocate(vfield%val(1:vfield%dim,1:node_count(vfield)))
    !vfield%val=reshape(tfield%contiguous_val,[vfield%dim,&!
    vfield%val=reshape(tfield%val,[vfield%dim,&
        node_count(vfield)])
#endif

end function as_packed_vector

subroutine finalise_multistate(packed_state,multiphase_state,&
    multicomponent_state)

    type(state_type) :: packed_state
    type(state_type), dimension(:), pointer :: multiphase_state, multicomponent_state


    call deallocate(multiphase_state)
    deallocate(multiphase_state)
    call deallocate(multicomponent_state)
    deallocate(multicomponent_state)
    call deallocate(packed_state)

end subroutine finalise_multistate




subroutine Adaptive_NonLinear(packed_state, reference_field, its,&
    Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,order, adapt_mesh_in_FPI, calculate_mass_delta)
    !This subroutine either store variables before the nonlinear timeloop starts, or checks
    !how the nonlinear iterations are going and depending on that increase the timestep
    !or decreases the timestep and repeats that timestep
    Implicit none
    type(state_type), intent(inout) :: packed_state
    real, dimension(:,:,:), allocatable, intent(inout) :: reference_field
    logical, intent(inout) :: Repeat_time_step, ExitNonLinearLoop
    integer, intent(inout) :: its!not to be modified unless VERY sure
    logical, intent(in) :: nonLinearAdaptTs
    integer, intent(in) :: order
    logical, optional, intent(in) :: adapt_mesh_in_FPI
    !! 1st item holds the mass at previous Linear time step, 2nd item is the delta between mass at the current FPI and 1st item
    real, dimension(:,:), optional :: calculate_mass_delta
    !Local variables
    integer, save :: nonlinear_its=0!Needed for adapt_within_fpi to consider all the non-linear iterations together
    real, save :: stored_dt = -1
    logical, save :: adjusted_ts_to_dump = .false.
    real :: dt, auxR, dump_period
    integer :: Aim_num_FPI, auxI, incr_threshold
    integer, save :: show_FPI_conv
    real, save :: OldDt
    real, parameter :: check_sat_threshold = 1d-6
    real, dimension(:,:,:), pointer :: pressure
    real, dimension(:,:), pointer :: phasevolumefraction, temperature
    real, dimension(:,:,:), pointer :: velocity
    character (len = OPTION_PATH_LEN) :: output_message =''
    !Variables for automatic non-linear iterations
    real, save :: dt_by_user = -1
    real :: tolerance_between_non_linear, min_ts, max_ts,&
        Infinite_norm_tol, calculate_mass_tol
    !! local variable, holds the maximum mass error
    real :: max_calculate_mass_delta
    real, dimension(2) :: totally_min_max
    !Variables for PID time-step size controller
    logical :: PID_controller
    !Variables for adaptive time stepping based on non-linear iterations
    real :: increaseFactor, decreaseFactor, ts_ref_val, acctim, inf_norm_val, finish_time
    integer :: variable_selection, NonLinearIteration
    !We need an acumulative nonlinear_its if adapting within the FPI we don't want to restart the reference field neither
    !consider less iterations of the total ones if adapting time using PID
    if (.not.have_option( '/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI')) then
        nonlinear_its = its
    else
        if (.not.ExitNonLinearLoop) then
            nonlinear_its = its!Only do something different when we are suppose to exit
        else!Store when we are in theory finishing
            nonlinear_its = nonlinear_its + its
        end if
    end if
    !ewrite(0,*) "entering"
    !First of all, check if the user wants to do something
    call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration', tolerance_between_non_linear, default = -1. )
    if (tolerance_between_non_linear<0) return
    !Tolerance for the infinite norm
    call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Infinite_norm_tol',&
        Infinite_norm_tol, default = 0.03 )
    !retrieve number of Fixed Point Iterations
    call get_option( '/timestepping/nonlinear_iterations', NonLinearIteration, default = 3 )
    !Get data from diamond. Despite this is slow, as it is done in the outest loop, it should not affect the performance.
    !Variable to check how good nonlinear iterations are going 1 (Pressure), 2 (Velocity), 3 (Saturation), 4 (Temperature)
    call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Infinite_norm_tol/adaptive_non_linear_iterations', &
        variable_selection, default = 3)!by default saturation so it is effectively disabled for single phase
    if (have_option('/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear')) then
        call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear', &
            variable_selection, default = 3)!by default saturation so it is effectively disabled for single phase
    end if
    call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/increase_factor', &
        increaseFactor, default = 1.1 )
    call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/decrease_factor', &
        decreaseFactor, default = 2.0 )
    call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/max_timestep', &
        max_ts, default = huge(min_ts) )
    if (dt_by_user < 0) call get_option( '/timestepping/timestep', dt_by_user )
    call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/min_timestep', &
        min_ts, default = dt_by_user*1d-3 )
    call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/increase_threshold', &
        incr_threshold, default = int(0.25 * NonLinearIteration) )
    show_FPI_conv = .not.have_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Show_Convergence')
    call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/PID_controller/Aim_num_FPI', &
        Aim_num_FPI, default = int(0.20 * NonLinearIteration) )
    call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Test_mass_consv', &
            calculate_mass_tol, default = 5d-3)
    PID_controller = have_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/PID_controller')
    !Retrieve current time and final time
    call get_option( '/timestepping/current_time', acctim )
    call get_option( '/timestepping/finish_time', finish_time )
    !Ensure that even adapting the time, the final time is matched
    max_ts = max(min(max_ts, abs(finish_time - acctim)), 1e-8)
    if (stored_dt<0) then!for the first time only
        call get_option( '/timestepping/timestep', dt )
        stored_dt = dt
    end if
    !To ensure that we always create a vtu file at the desired time,
    !we control the maximum time-step size to ensure that at some point the ts changes to provide that precise time
    if (have_option('/io/dump_period/constant')) then
        call get_option( '/io/dump_period/constant', dump_period )
        !First get the next time for a vtu dump
        auxR = dble(ceiling(acctim/dump_period)) * dump_period
        if (abs(auxR-acctim) > 1e-12) then
            max_ts = max(min(max_ts, abs(acctim-auxR)), min_ts*1d-3)!Make sure we dump at the required time and we don't get dt = 0
        else
            max_ts = min(max_ts, dump_period)
        end if
    end if


    select case (order)
        case (1)!Store or get from backup
            !If we do not have adaptive time stepping then there is nothing to backup
            if (.not.nonLinearAdaptTs) return
            !we either store the data or we recover it if repeting a timestep
            !Procedure to repeat time-steps
            !If  Repeat_time_step then we recover values, else we store them
            if (nonlinear_its == 1) call copy_packed_new_to_iterated(packed_state, Repeat_time_step)
        case (2)!Calculate and store reference_field
            !Store variable to check afterwards
            call get_var_from_packed_state(packed_state, velocity = velocity, pressure = pressure,&
                phasevolumefraction = phasevolumefraction)

            select case (variable_selection)
                case (2)!Velocity
                    if (allocated(reference_field)) then
                        if (size(reference_field,1) /= size(velocity,1) .or. &
                            size(reference_field,2) /= size(velocity,2) .or. &
                            size(reference_field,3) /= size(velocity,3) ) then
                            deallocate(reference_field)
                            allocate (reference_field(size(velocity,1),size(velocity,2),size(velocity,3) ))
                        end if
                    else
                        allocate (reference_field(size(velocity,1),size(velocity,2),size(velocity,3) ))
                    end if
                    reference_field(:,:,:) = velocity
                case (3)!Phase volume fraction

                    if (allocated(reference_field)) then
                        if (size(reference_field,2) /= size(phasevolumefraction,1) .or. &
                            size(reference_field,3) /= size(phasevolumefraction,2) ) then
                            deallocate(reference_field)
                            allocate (reference_field(1,size(phasevolumefraction,1),size(phasevolumefraction,2) ))
                        end if
                    else
                        allocate (reference_field(1,size(phasevolumefraction,1),size(phasevolumefraction,2) ))
                    end if
                    reference_field(1,:,:) = phasevolumefraction
                case (4)!Temperature
                    call get_var_from_packed_state(packed_state, temperature = temperature)
                    if (allocated(reference_field)) then
                        if (size(reference_field,2) /= size(temperature,1) .or. &
                            size(reference_field,3) /= size(temperature,2) ) then
                            deallocate(reference_field)
                            !If temperature, also keep and eye on saturation with the other convergence criterion
                            allocate (reference_field(2,size(temperature,1),size(temperature,2) ))
                        end if
                    else
                        allocate (reference_field(2,size(temperature,1),size(temperature,2) ))
                    end if
                    reference_field(1,:,:) = temperature
                    reference_field(2,:,:) = phasevolumefraction
                case default !Default as pressure is always defined and changes more smoothly than velocity
                    if (allocated(reference_field)) then
                        if (size(reference_field,3) /= size(pressure,3) ) then
                            deallocate(reference_field)
                            allocate (reference_field(1,1,size(pressure,3) ))
                        end if
                    else
                        allocate (reference_field(1,1,size(pressure,3) ))
                    end if
                    reference_field(1,1,:) = pressure(1,1,:)
            end select

        case default!Check how is the process going on and decide
            !If Automatic_NonLinerIterations then we compare the variation of the a property from one time step to the next one
            ExitNonLinearLoop = .false.
            Repeat_time_step = .false.
            call get_var_from_packed_state(packed_state, velocity = velocity, pressure = pressure,&
                phasevolumefraction = phasevolumefraction)

            select case (variable_selection)

                case (2)!Velocity
                    inf_norm_val = maxval(abs(reference_field-velocity))
                    ts_ref_val = inf_norm_val!Use the infinite norm for the time being
                    tolerance_between_non_linear = 1d9!Only infinite norm for the time being
                case (3)!Phase volume fraction
                    !Calculate infinite norm
                    inf_norm_val = maxval(abs(reference_field(1,:,:)-phasevolumefraction))/backtrack_or_convergence

                    !Calculate value of the functional
                    ts_ref_val = get_Convergence_Functional(phasevolumefraction, reference_field(1,:,:), backtrack_or_convergence, nonlinear_its)
                    backtrack_or_convergence = get_Convergence_Functional(phasevolumefraction, reference_field(1,:,:), backtrack_or_convergence)

                case (4)!Temperature
                    call get_var_from_packed_state(packed_state, temperature = temperature)
                    !Calculate normalized infinite norm of the difference
                                                            !This Mask is important because otherwise it gets the lowest saturation value
                    totally_min_max(1)=minval(reference_field, MASK = reference_field > 1.1)!Using Kelvin it is unlikely that the temperature gets to 1 Kelvin!
                    totally_min_max(2)=maxval(reference_field)!use stored temperature
                    !For parallel
                    call allmin(totally_min_max(1)); call allmax(totally_min_max(2))
                    !Analyse the difference
                    ts_ref_val = inf_norm_scalar_normalised(temperature, reference_field(1,:,:), 1.0, totally_min_max)
                    !Calculate value of the l infinitum for the saturation as well
                    inf_norm_val = maxval(abs(reference_field(2,:,:)-phasevolumefraction))/backtrack_or_convergence

                case default!Pressure
                    !Calculate normalized infinite norm of the difference
                    totally_min_max(1)=minval(reference_field)!use stored pressure
                    totally_min_max(2)=maxval(reference_field)!use stored pressure
                    !For parallel
                    call allmin(totally_min_max(1)); call allmax(totally_min_max(2))
                    !Analyse the difference
                    inf_norm_val = inf_norm_scalar_normalised(pressure(1,:,:), reference_field(1,:,:), 1.0, totally_min_max)
                    ts_ref_val = inf_norm_val!Use the infinite norm for the time being
                    tolerance_between_non_linear = 1d9!Only infinite norm for the time being
            end select
            ! find the maximum mass error to compare with the tolerance below
            ! This is the maximum error of each indivial phase
            max_calculate_mass_delta = calculate_mass_delta(1,2)

            !If it is parallel then we want to be consistent between cpus
            if (IsParallel()) then
                call allmax(ts_ref_val)
                call allmax(max_calculate_mass_delta)
                call allmax(inf_norm_val)
            end if
            !Store output messages
            if (is_porous_media .and. variable_selection == 3) then
                write(output_message, '(a, E10.3,a,E10.3, a, i0, a, E10.3)' )"FPI convergence: ",ts_ref_val,"; L_inf:", inf_norm_val, "; Total iterations: ", its, "; Mass error:", max_calculate_mass_delta
            else if (is_porous_media .and. variable_selection == 4) then
                write(output_message, '(a, E10.3,a,E10.3, a, i0, a, E10.3)' )"Temperature (L_inf): ",ts_ref_val,"; Saturation (L_inf):", inf_norm_val, "; Total iterations: ", its, "; Mass error:", max_calculate_mass_delta
            else
                write(output_message, '(a, E10.3,a,i0)' ) "L_inf:", inf_norm_val, "; Total iterations: ", its
            end if

            !TEMPORARY, re-use of global variable backtrack_or_convergence to send
            !information about convergence to the trust_region_method
            if (is_flooding) backtrack_or_convergence = ts_ref_val
            !Automatic non-linear iteration checking
            if (is_porous_media) then
                select case (variable_selection)
                    case (4)!For temperature only infinite norms for saturation and temperature
                        ExitNonLinearLoop = ((ts_ref_val < Infinite_norm_tol .and. inf_norm_val < Infinite_norm_tol &
                            .and. max_calculate_mass_delta < calculate_mass_tol ) .or. its >= NonLinearIteration )
                    case default
                        !For very tiny time-steps ts_ref_val may not be good as is it a relative value
                        !So if the infinity norm is way better than the tolerance we consider that the convergence have been achieved
                        if (inf_norm_val * 1e1 < Infinite_norm_tol) ts_ref_val = tolerance_between_non_linear/2.
                        ExitNonLinearLoop = ((ts_ref_val < tolerance_between_non_linear .and. inf_norm_val < Infinite_norm_tol &
                            .and. max_calculate_mass_delta < calculate_mass_tol ) .or. its >= NonLinearIteration )
                end select
            else
                ExitNonLinearLoop = (inf_norm_val < Infinite_norm_tol) .or. its >= NonLinearIteration
            end if
            !At least two non-linear iterations
            ExitNonLinearLoop =  ExitNonLinearLoop .and. its >= 2

            !(Maybe unnecessary) If it is parallel then we want to be consistent between cpus
            if (IsParallel()) call alland(ExitNonLinearLoop)
            !Tell the user the number of FPI and final convergence to help improving the parameters
            if (ExitNonLinearLoop .and. getprocno() == 1) then
                ewrite(show_FPI_conv,*) trim(output_message)
            end if
            !If time adapted based on the non-linear solver then
            if (nonLinearAdaptTs .and. .not. adapt_mesh_in_FPI) then!Do not adapt time if we are adapting the mesh within the FPI and
                                                                    !this is the first guess
                !If any solver fails to converge (and the user care), we may want to repeat the time-level
                !without waiting for the last non-linear iteration
                if (solver_not_converged .and.have_option(&
                  '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/ensure_solvers_convergence')) then
                    Repeat_time_step = .true.
                    solver_not_converged = .false.
                    if (getprocno() == 1) then
                        ewrite(show_FPI_conv,*) "WARNING: A solver failed to achieve convergence in the current non-linear iteration. Repeating time-level."
                    end if
                end if
                !If maximum number of FPI reached, then repeat time-step
                if (its >= NonLinearIteration) Repeat_time_step = .true.

                !If dt was modified just to match a dump_period then we impose again the previous time-step
                if (adjusted_ts_to_dump) then
                    dt = max(min(stored_dt, max_ts), 1d-8)
                    call set_option( '/timestepping/timestep', dt )
                    if (getprocno() == 1)then
                        ewrite(show_FPI_conv,*) "Time step restored to:", dt
                    end if
                    adjusted_ts_to_dump = .false.
                    return
                end if

                !This controller is supposed to be the most effective
                if (PID_controller) then
                    !We do not follow the normal approach
                    if (ExitNonLinearLoop.and..not.Repeat_time_step) then
                        !Modify time step
                        dt = stored_dt
                        auxR = PID_time_controller()
                        if (auxR < 1.0 )then!Reduce Ts
                            dt = max(dt * max(abs(auxR), 1./(1.5*decreaseFactor)), min_ts)
                        else
                            dt = dt * min(abs(auxR), 1.5*increaseFactor)
                        end if
                        auxR = stored_dt
                        call set_option( '/timestepping/timestep', dt )
                        stored_dt = dt
                        !Ensure that period_vtus or the final time are matched, controlled by max_ts
                        dt = max(min(dt, max_ts), min_ts)
                        call set_option( '/timestepping/timestep', dt )
                        if (getprocno() == 1 .and. abs(auxR-dt)/dt > 1d-3)then
                            ewrite(show_FPI_conv,*) "Time step changed to:", dt
                        end if
                        ExitNonLinearLoop = .true.
                        return
                    end if
                end if

                !Adaptive Ts for Backtracking only based on the number of FPI
                if (ExitNonLinearLoop .and. its < incr_threshold .and..not.Repeat_time_step) then
                    !Increase time step
                    dt = stored_dt!retrieve stored_dt
                    dt = dt * increaseFactor
                    call set_option( '/timestepping/timestep', dt )
                    stored_dt = dt
                    if (getprocno() == 1) then
                        ewrite(show_FPI_conv,*) "Time step increased to:", dt
                    end if
                    ExitNonLinearLoop = .true.
                    !Ensure that period_vtus or the final time are matched, controlled by max_ts
                    dt = max(min(dt, max_ts), min_ts)
                    call set_option( '/timestepping/timestep', dt )
                    return
                end if
                if (its >= NonLinearIteration .or. Repeat_time_step) then
                    !If it has not converged when reaching the maximum number of non-linear iterations,
                    !reduce ts and repeat
                    dt = stored_dt!retrieve stored_dt
                    if ( dt - min_ts < 1d-8) then
                        !Ensure that dt = min_ts
                        dt = min_ts
                        call set_option( '/timestepping/timestep', dt )
                        stored_dt = dt
                        !Do not decrease if minimum ts is reached
                        Repeat_time_step = .false.
                        ExitNonLinearLoop = .true.
                        deallocate(reference_field)
                        !Tell the user the number of FPI and final convergence to help improving the parameters
                        if (getprocno() == 1) then
                            ewrite(show_FPI_conv,*)  "Minimum time-step(",min_ts,") reached, advancing time."
                        end if
                        !If PID_controller then update the status
                        if (PID_controller) auxR = PID_time_controller(reset=.true.)
                        return
                    end if
                    !Decrease time step, reset the time and repeat!
                    call get_option( '/timestepping/current_time', acctim )
                    acctim = acctim - dt
                    call set_option( '/timestepping/current_time', acctim )

                    if (PID_controller) then
                        auxR = PID_time_controller()
                        !Maybe the PID controller thinks is better to reduce more than just half, up to 0.25
                        dt = max(min(dt / decreaseFactor, max( auxR, 0.5*dt / decreaseFactor)), min_ts)
                        !If PID_controller then update the status
                        auxR = PID_time_controller(reset=.true.)
                    else
                        dt = max(dt / decreaseFactor,min_ts)
                    end if
                    call set_option( '/timestepping/timestep', dt )
                    stored_dt = dt
                    if (getprocno() == 1) then
                        ewrite(show_FPI_conv,*) "<<<Convergence not achieved, repeating time-level>>> Time step decreased to:", dt
                    end if
                    Repeat_time_step = .true.
                    ExitNonLinearLoop = .true.
                    return
                end if
                if (ExitNonLinearLoop.and..not.Repeat_time_step) then
                    !If adapting ts and it gets here -meaning it has not been modified-, maybe we still need to adapt ts
                    !to ensure we match the final time or a period_vtu
                    dt = stored_dt!retrieve stored_dt
                    auxR = dt!Store dt before modification to compare
                    dt = max(min(dt, max_ts), 1d-8)
                    !here we do not store dt, as its modification is not based on stability
                    call set_option( '/timestepping/timestep', dt )
                    if (abs(auxR-dt) > 1d-8) then
                        if (getprocno() == 1)then
                            ewrite(show_FPI_conv,*) "Time step modified to match final time/dump_period:", dt
                            adjusted_ts_to_dump = .true.
                        end if
                    end if
                    return
                end if
            end if
    end select

contains

    real function PID_time_controller(reset)
        !This functions calculates the multiplier to get a new time-step size based on a
        !Proportional-Integral-Derivative (PID) concept. See: SPE-182601-MS
        implicit none
        logical, optional, intent(in) :: reset
        !Local variables
        real, parameter:: Ki = 1.34, Kd = 0.01, Kp = 0.001 !Fixed values from the paper, this can be improved, see SPE-182601-MS
        real, save :: Cn1 = -1, Cn2 = -1
        real, dimension(3) :: Cn
        real :: aux
        real, parameter :: tol = 1e-8
        logical, parameter :: max_criteria = .false.!If false, use an average with different weights

        if (present_and_true(reset))then
            Cn1 = -1; Cn2 = -1
        end if
        Cn = 0.; aux = 0.
        !Calculate Cn
        if (is_porous_media.and. .not.first_time_step) then
            Cn(1) = ts_ref_val/tolerance_between_non_linear
            aux = aux + 1.0
        end if
        !Compare with infinitum norm
        Cn(2) = inf_norm_val/Infinite_norm_tol
        aux = aux + 1.0
        !Maybe consider as well aiming to a certain number of FPIs
        if (Aim_num_FPI > 0) then                     !Options for the exponent:
            Cn(3) = (dble(nonlinear_its)/dble(Aim_num_FPI))**0.9! 2.0 => too strongly enforce the number of iterations, ignores other criteria
            aux = aux + 1.0                           ! 1.0 => default value, forces the number of iterations, almost ignore other criteria
        end if                                        ! 0.6 => soft constrain, it will try but not very much, considers other criteria
        if (max_criteria) then
            Cn(1) = maxval(Cn)
        else
            Cn(1) = (sum(Cn)+maxval(Cn))/(aux+1.)
        end if
        Cn(1) = (Cn(1) + tol)! <= To avoid divisions by zero
        if (Cn2 > 0) then
            PID_time_controller = (1./Cn(1))**Ki * (Cn1/Cn(1))**Kp * (Cn1**2. / (Cn(1)*Cn2))**Kd
        else if (Cn1 > 0) then
            PID_time_controller = (1./Cn(1))**Ki * (Cn1/Cn(1))**Kp
        else!Not enough information
            PID_time_controller = (1./Cn(1))**1.0
        end if

        !Store previous values
        Cn2 = Cn1; Cn1 = Cn(1)
    end function PID_time_controller

end subroutine Adaptive_NonLinear

real function inf_norm_scalar_normalised(tracer, reference_tracer, dumping, totally_min_max)
    !Calculate the inf norm of the normalised field, so the field goes from 0 to 1
    implicit none
    real, dimension(:,:), intent(in) :: tracer, reference_tracer
    real, intent(in) :: dumping
    real, dimension(2), intent(in) :: totally_min_max
    !Local variables
    integer :: cv_inod, iphase
    !Same as normilising values but should be quicker
    inf_norm_scalar_normalised = maxval(abs(reference_tracer-tracer))/max((totally_min_max(2)-totally_min_max(1)), 1d-8)
    call allmax(inf_norm_scalar_normalised)
    !rescale with accumulated dumping, if no dumping just pass down a 1.0
    inf_norm_scalar_normalised = inf_norm_scalar_normalised/dumping

end function



real function get_Convergence_Functional(phasevolumefraction, reference_sat, dumping, its)
    !We create a potential to optimize F = sum (f**2), so the solution is when this potential
    !reaches a minimum. Typically the value to consider convergence is the sqrt(epsilon of the machine), i.e. 10^-8
    !f = (NewSat-OldSat)/Number of nodes; this is the typical approach for algebraic non linear systems
    !
    !The convergence is independent of the dumping parameter
    !and measures how the previous iteration (i.e. using the previous dumping parameter) performed
    implicit none
    real, dimension(:,:), intent(in) :: phasevolumefraction, reference_sat
    real, intent(in) :: dumping
    integer, optional, intent(in) :: its
    !Local variables
    real, save :: First_potential
    integer :: cv_inod, modified_vals, iphase
    real, parameter :: tol = 1d-5
    real :: tmp ! Variable used for parallel consistency

    modified_vals = 0
    get_Convergence_Functional = 0.0

    !(L2)**2 norm of all the elements
    do iphase = 1, size(phasevolumefraction,1)

        tmp = sum((abs(reference_sat(iphase,:)-phasevolumefraction(iphase,:))/size(phasevolumefraction,2))**2.0)
        call allsum(tmp)

        get_Convergence_Functional = max(tmp, get_Convergence_Functional)
    end do

!        !(L2)**2 norm of all the elements whose value has changed (has problems to converge at
!        !the beginning since only a shock front is happening, however it is better than simple l2 norm)
!        do cv_inod = 1, size(phasevolumefraction,2)
!            aux = maxval(abs(reference_sat(:,cv_inod)-phasevolumefraction(:,cv_inod)))
!            if (aux > tol) then
!                get_Convergence_Functional = get_Convergence_Functional + (aux)**2.0
!                modified_vals = modified_vals + 1
!            end if
!        end do
!        get_Convergence_Functional = (get_Convergence_Functional / dble(modified_vals)**2.0)

    !Rescale using the dumping in saturation to get a more efficient number to compare with
    !if the backtrack_or_convergence was 10-2 then ts_ref_val will always be small
    !To make consistent the dumping parameter with the Potential, we have to raise it to 2.0
    get_Convergence_Functional = get_Convergence_Functional / dumping**2.0

    if (present(its)) then
        if (its == 1) then
            First_potential = get_Convergence_Functional
        else
            !It could happen that the first potential is effectively zero if using pressure boundary conditions and/or small ts
            !if that is the case we allow to update the first potential up to two times more
            if (First_potential * 1d10 < get_Convergence_Functional .and. its <= 3) First_potential = 2.0*get_Convergence_Functional
            get_Convergence_Functional = get_Convergence_Functional/First_potential
        end if
    end if

end function get_Convergence_Functional

subroutine copy_packed_new_to_iterated(packed_state, viceversa)
    !Values from packed_state are stored in iterated unless viceversa is true, in that case
    !the iterated values are moved to the new values
    type(state_type), intent(inout) :: packed_state
    logical, intent(in) :: viceversa

    type(scalar_field), pointer :: sfield, nsfield
    type(vector_field), pointer :: vfield, nvfield
    type(tensor_field), pointer :: tfield, ntfield

    integer :: i

    do i=1,size(packed_state%scalar_fields)
        sfield=>packed_state%scalar_fields(i)%ptr
        if (sfield%name(1:14)=="PackedIterated") then
            nsfield=>extract_scalar_field(packed_state,"Packed"//sfield%name(15:))
            if (viceversa) then
                nsfield%val = sfield%val
            else
                sfield%val=nsfield%val
            end if
        end if
    end do

    do i=1,size(packed_state%vector_fields)
        vfield=>packed_state%vector_fields(i)%ptr
        if (vfield%name(1:14)=="PackedIterated") then
            nvfield=>extract_vector_field(packed_state,"Packed"//vfield%name(15:))
            if (viceversa) then
                nvfield%val = vfield%val
            else
                vfield%val=nvfield%val
            end if
        end if
    end do

    do i=1,size(packed_state%tensor_fields)
        tfield=>packed_state%tensor_fields(i)%ptr
        if (tfield%name(1:14)=="PackedIterated") then
            ntfield=>extract_tensor_field(packed_state,"Packed"//tfield%name(15:))
            if (viceversa) then
                ntfield%val  = tfield%val
            else
                tfield%val=ntfield%val
            end if
        end if
    end do




end subroutine copy_packed_new_to_iterated

!deprecated, do not use. Use pointers instead
subroutine get_var_from_packed_state(packed_state,FEDensity,&
    OldFEDensity,IteratedFEDensity,Density,OldDensity,IteratedDensity,PhaseVolumeFraction,&
    OldPhaseVolumeFraction,IteratedPhaseVolumeFraction, Velocity, OldVelocity, IteratedVelocity, &
    FEPhaseVolumeFraction, OldFEPhaseVolumeFraction, IteratedFEPhaseVolumeFraction,&
    NonlinearVelocity, OldNonlinearVelocity,IteratedNonlinearVelocity, ComponentDensity, &
    OldComponentDensity, IteratedComponentDensity,ComponentMassFraction, OldComponentMassFraction,&
    Temperature,OldTemperature, IteratedTemperature,FETemperature, OldFETemperature, IteratedFETemperature,&
    IteratedComponentMassFraction, FEComponentDensity, OldFEComponentDensity, IteratedFEComponentDensity,&
    FEComponentMassFraction, OldFEComponentMassFraction, IteratedFEComponentMassFraction,&
    Pressure,FEPressure, OldFEPressure, CVPressure,OldCVPressure,&
    Coordinate, VelocityCoordinate,PressureCoordinate,MaterialCoordinate, CapPressure, Immobile_fraction,&
    EndPointRelperm, RelpermExponent, Cap_entry_pressure, Cap_exponent, Imbibition_term)
    !This subroutine returns a pointer to the desired values of a variable stored in packed state
    !All the input variables (but packed_stated) are pointers following the structure of the *_ALL variables
    !and also all of them are optional, hence you can obtaine whichever you want
    !######################EXAMPLE OF USAGE OF THIS SUBROUTINE:#####################################
    !If we want to get the velocity and the phasevolumefraction one should proceed this way:
    !Define variables:
    !real, dimension(:,:,:), pointer :: Velocity_pointer
    !real, dimension(:,:), pointer :: PhaseVolumeFraction_pointer
    !Assign the pointers
    !call get_var_from_packed_state(packed_state, Velocity = Velocity_pointer, PhaseVolumeFraction = PhaseVolumeFraction_pointer)
    !
    ! In this way we only have to introduce the name of the variables we want to get from packed_state
    !########################################################################################
    implicit none
    type(state_type), intent(inout) :: packed_state
    real, optional, dimension(:,:,:), pointer :: Velocity, OldVelocity, IteratedVelocity, NonlinearVelocity, OldNonlinearVelocity,&
        IteratedNonlinearVelocity,ComponentDensity, OldComponentDensity, IteratedComponentDensity,ComponentMassFraction, OldComponentMassFraction,&
        IteratedComponentMassFraction, FEComponentDensity, OldFEComponentDensity, IteratedFEComponentDensity, FEComponentMassFraction, &
        OldFEComponentMassFraction, IteratedFEComponentMassFraction
    real, optional, dimension(:,:), pointer :: FEDensity, OldFEDensity, IteratedFEDensity, Density,&
        OldDensity,IteratedDensity,PhaseVolumeFraction,OldPhaseVolumeFraction,IteratedPhaseVolumeFraction,&
        Temperature, OldTemperature, IteratedTemperature, FETemperature, OldFETemperature, IteratedFETemperature,&
        Coordinate, VelocityCoordinate,PressureCoordinate,MaterialCoordinate, &
        FEPhaseVolumeFraction, OldFEPhaseVolumeFraction, IteratedFEPhaseVolumeFraction, CapPressure,&
        Immobile_fraction, EndPointRelperm, RelpermExponent, Cap_entry_pressure, Cap_exponent, Imbibition_term
    real, optional, dimension(:,:,:), pointer ::Pressure,FEPressure, OldFEPressure, CVPressure,OldCVPressure
    !Local variables
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield

    !Scalar stored
    if (present(Pressure)) then
        tfield => extract_tensor_field( packed_state, "PackedFEPressure" )
        Pressure => tfield%val(:,:,:)
    end if
    if (present(FEPressure)) then
        tfield => extract_tensor_field( packed_state, "PackedFEPressure" )
        FEPressure => tfield%val(:,:,:)
    end if
    if (present(OldFEPressure)) then
        tfield => extract_tensor_field( packed_state, "PackedOldFEPressure" )
        OldFEPressure => tfield%val(:,:,:)
    end if
    if (present(CVPressure)) then
        tfield => extract_tensor_field( packed_state, "PackedCVPressure" )
        CVPressure => tfield%val(:,:,:)
    end if
    if (present(OldCVPressure)) then
        tfield => extract_tensor_field( packed_state, "PackedOldCVPressure" )
        OldCVPressure => tfield%val(:,:,:)
    end if

    !Vectors stored
    if (present(Coordinate)) then
        vfield => extract_vector_field( packed_state, "Coordinate" )
        Coordinate =>  vfield%val(:,:)
    end if
    if (present(VelocityCoordinate)) then
        vfield => extract_vector_field( packed_state, "VelocityCoordinate" )
        VelocityCoordinate =>  vfield%val(:,:)
    end if
    if (present(PressureCoordinate)) then
        vfield => extract_vector_field( packed_state, "PressureCoordinate" )
        PressureCoordinate =>  vfield%val(:,:)
    end if
    if (present(MaterialCoordinate)) then
        vfield => extract_vector_field( packed_state, "MaterialCoordinate" )
        MaterialCoordinate =>  vfield%val(:,:)
    end if

    !Tensors stored
    if (present(FEDensity)) then
        tfield => extract_tensor_field( packed_state, "PackedFEDensity" )
        FEDensity =>  tfield%val(1,:,:)
    end if
    if (present(OldFEDensity)) then
        tfield => extract_tensor_field( packed_state, "PackedOldFEDensity" )
        OldFEDensity =>  tfield%val(1,:,:)
    end if
    if (present(IteratedFEDensity)) then
        tfield => extract_tensor_field( packed_state, "PackedIteratedFEDensity" )
        IteratedFEDensity => tfield%val(1,:,:)
    end if
    if (present(Density)) then
        tfield => extract_tensor_field( packed_state, "PackedDensity" )
        Density => tfield%val(1,:,:)
    end if
    if (present(OldDensity)) then
        tfield => extract_tensor_field( packed_state, "PackedOldDensity" )
        OldDensity => tfield%val(1,:,:)
    end if
    if (present(IteratedDensity)) then
        tfield => extract_tensor_field( packed_state, "PackedIteratedDensity" )
        IteratedDensity =>  tfield%val(1,:,:)
    end if
    if (present(PhaseVolumeFraction)) then
        tfield => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
        PhaseVolumeFraction =>  tfield%val(1,:,:)
    end if
    if (present(OldPhaseVolumeFraction)) then
        tfield => extract_tensor_field( packed_state, "PackedOldPhaseVolumeFraction" )
        OldPhaseVolumeFraction =>  tfield%val(1,:,:)
    end if
    if (present(IteratedPhaseVolumeFraction)) then
        tfield => extract_tensor_field( packed_state, "PackedIteratedPhaseVolumeFraction" )
        IteratedPhaseVolumeFraction =>  tfield%val(1,:,:)
    end if
    if (present(FEPhaseVolumeFraction)) then
        tfield => extract_tensor_field( packed_state, "PackedFEPhaseVolumeFraction" )
        FEPhaseVolumeFraction =>  tfield%val(1,:,:)
    end if
    if (present(OldFEPhaseVolumeFraction)) then
        tfield => extract_tensor_field( packed_state, "PackedOldFEPhaseVolumeFraction" )
        OldFEPhaseVolumeFraction =>  tfield%val(1,:,:)
    end if
    if (present(IteratedFEPhaseVolumeFraction)) then
        tfield => extract_tensor_field( packed_state, "PackedIteratedFEPhaseVolumeFraction" )
        IteratedFEPhaseVolumeFraction =>  tfield%val(1,:,:)
    end if
    if (present(CapPressure)) then
        tfield => extract_tensor_field( packed_state, "PackedCapPressure" )
        CapPressure =>  tfield%val(1,:,:)
    end if
    if (present(Temperature)) then
        tfield => extract_tensor_field( packed_state, "PackedTemperature" )
        Temperature =>  tfield%val(1,:,:)
    end if
    if (present(OldTemperature)) then
        tfield => extract_tensor_field( packed_state, "PackedOldTemperature" )
        OldTemperature =>  tfield%val(1,:,:)
    end if
    if (present(IteratedTemperature)) then
        tfield => extract_tensor_field( packed_state, "PackedIteratedTemperature" )
        IteratedTemperature =>  tfield%val(1,:,:)
    end if
    if (present(FETemperature)) then
        tfield => extract_tensor_field( packed_state, "PackedFETemperature" )
        FETemperature =>  tfield%val(1,:,:)
    end if
    if (present(OldFETemperature)) then
        tfield => extract_tensor_field( packed_state, "PackedOldFETemperature" )
        OldFETemperature =>  tfield%val(1,:,:)
    end if
    if (present(IteratedFETemperature)) then
        tfield => extract_tensor_field( packed_state, "PackedIteratedFETemperature" )
        IteratedFETemperature =>  tfield%val(1,:,:)
    end if

    if (present(Velocity)) then
        tfield => extract_tensor_field( packed_state, "PackedVelocity" )
        Velocity => tfield%val(:,:,:)
    end if
    if (present(OldVelocity)) then
        tfield => extract_tensor_field( packed_state, "PackedOldVelocity" )
        OldVelocity => tfield%val(:,:,:)
    end if
    if (present(IteratedVelocity)) then
        tfield => extract_tensor_field( packed_state, "PackedIteratedVelocity" )
        IteratedVelocity => tfield%val(:,:,:)
    end if
    if (present(NonlinearVelocity)) then
        tfield => extract_tensor_field( packed_state, "PackedNonlinearVelocity" )
        NonlinearVelocity => tfield%val(:,:,:)
    end if
    if (present(OldNonlinearVelocity)) then
        tfield => extract_tensor_field( packed_state, "PackedOldNonlinearVelocity" )
        OldNonlinearVelocity => tfield%val(:,:,:)
    end if
    if (present(IteratedNonlinearVelocity)) then
        tfield => extract_tensor_field( packed_state, "PackedIteratedNonlinearVelocity" )
        IteratedNonlinearVelocity => tfield%val(:,:,:)
    end if
    if (present(ComponentDensity)) then
        tfield => extract_tensor_field( packed_state, "PackedComponentDensity" )
        ComponentDensity => tfield%val(:,:,:)
    end if
    if (present(OldComponentDensity)) then
        tfield => extract_tensor_field( packed_state, "PackedOldComponentDensity" )
        OldComponentDensity => tfield%val(:,:,:)
    end if
    if (present(IteratedComponentDensity)) then
        tfield => extract_tensor_field( packed_state, "PackedIteratedComponentDensity" )
        IteratedComponentDensity => tfield%val(:,:,:)
    end if
    if (present(ComponentMassFraction)) then
        tfield => extract_tensor_field( packed_state, "PackedComponentMassFraction" )
        ComponentMassFraction => tfield%val(:,:,:)
    end if
    if (present(OldComponentMassFraction)) then
        tfield => extract_tensor_field( packed_state, "PackedOldComponentMassFraction" )
        OldComponentMassFraction => tfield%val(:,:,:)
    end if
    if (present(IteratedComponentMassFraction)) then
        tfield => extract_tensor_field( packed_state, "PackedIteratedComponentMassFraction" )
        IteratedComponentMassFraction => tfield%val(:,:,:)
    end if
    if (present(FEComponentDensity)) then
        tfield => extract_tensor_field( packed_state, "PackedFEComponentDensity" )
        FEComponentDensity => tfield%val(:,:,:)
    end if
    if (present(OldFEComponentDensity)) then
        tfield => extract_tensor_field( packed_state, "PackedOldFEComponentDensity" )
        OldFEComponentDensity => tfield%val(:,:,:)
    end if
    if (present(IteratedFEComponentDensity)) then
        tfield => extract_tensor_field( packed_state, "PackedIteratedFEComponentDensity" )
        IteratedFEComponentDensity => tfield%val(:,:,:)
    end if
    if (present(FEComponentMassFraction)) then
        tfield => extract_tensor_field( packed_state, "PackedFEComponentMassFraction" )
        FEComponentMassFraction => tfield%val(:,:,:)
    end if
    if (present(OldFEComponentMassFraction)) then
        tfield => extract_tensor_field( packed_state, "PackedOldFEComponentMassFraction" )
        OldFEComponentMassFraction => tfield%val(:,:,:)
    end if
    if (present(IteratedFEComponentMassFraction)) then
        tfield => extract_tensor_field( packed_state, "PackedIteratedFEComponentMassFraction" )
        IteratedFEComponentMassFraction => tfield%val(:,:,:)
    end if
    if (present(Immobile_fraction))then
        tfield => extract_tensor_field( packed_state, "PackedRockFluidProp" )
        Immobile_fraction => tfield%val(1,:,:)
    end if
    if (present(EndPointRelperm))then
        tfield => extract_tensor_field( packed_state, "PackedRockFluidProp" )
        EndPointRelperm => tfield%val(2,:,:)
    end if
    if (present(RelpermExponent))then
        tfield => extract_tensor_field( packed_state, "PackedRockFluidProp" )
        RelpermExponent => tfield%val(3,:,:)
    end if
    if (present(Cap_entry_pressure))then
        tfield => extract_tensor_field( packed_state, "PackedRockFluidProp" )
        Cap_entry_pressure => tfield%val(4,:,:)
    end if
    if (present(Cap_exponent))then
        tfield => extract_tensor_field( packed_state, "PackedRockFluidProp" )
        Cap_exponent => tfield%val(5,:,:)
    end if
    if (present(Imbibition_term))then
        tfield => extract_tensor_field( packed_state, "PackedRockFluidProp" )
        Imbibition_term => tfield%val(6,:,:)
    end if

end subroutine get_var_from_packed_state


!Subroutine to print CSR matrix by (row, column)
!Dimensions and phases are printed in different rows
!So for example Matrix(2,2,10) with two rows would be presented as
!a matrix ( 8 x 10)
subroutine printCSRMatrix(Matrix, find, col, dim_same_row)
    implicit none
    integer, intent(in), dimension(:) :: find, col
    real, intent(in), dimension(:,:,:):: Matrix
    logical, optional, intent(in) :: dim_same_row
    !Local
    Integer :: i,j, k, row, column, pos, ncols
    character (len=100000) :: cadena
    character (len=100) :: aux
    real :: val

    ncols = maxval(col)

    if (.not.present_and_true(dim_same_row)) then
        do row = 1, size(find)-1
            do i = 1, size(Matrix,2)
                do j = 1, size(Matrix,1)
                    print *,""!jump line
                    cadena = ""
                    k = 0
                    do column = 1, ncols!size(col)
                        if (col(find(row)+k) == column .and. k < find(row + 1) - find( row ) ) then
                            val = Matrix(j, i, find( row )+k)
                            k = k + 1
                        else
                            val = 0.0
                        end if
                        write(aux,*), val
                        pos = index(trim(aux),"E",.true.)
                        if (pos/=0) then
                            aux = aux(1:pos-6)//trim(aux(pos:))
                        end if
                        cadena = trim(cadena)//' '//trim(aux)//','
                    end do
                    print '(A $)', trim(cadena)!print line
                end do
            end do
        end do
    else

        do row = 1, size(find)-1
            print *,""!jump line
            cadena = ""
            k = 0
            do column = 1, ncols!size(col)
                do i = 1, size(Matrix,2)
                    do j = 1, size(Matrix,1)
                        if (col(find(row)+k) == column .and. k < find(row + 1) - find( row ) ) then
                            val = Matrix(j, i, find( row )+k)
                            if (j == size(Matrix,1) .and. i == size(Matrix,2)) k = k + 1
                        else
                            val = 0.0
                        end if
                        write(aux,*), val
                        pos = index(trim(aux),"E",.true.)
                        if (pos/=0) then
                            aux = aux(1:pos-6)//trim(aux(pos:))
                        end if
                        cadena = trim(cadena)//' '//trim(aux)//','
                    end do
                end do
            end do
            print '(A $)', trim(cadena)!print line
        end do



    end if
end subroutine printCSRMatrix

logical function is_constant(tfield)
    type(tensor_field), intent(in) :: tfield

    integer :: i
    real, dimension(product(tfield%dim)) :: tmax, tmin
    real, parameter :: tol=1.0e-8

    if (isparallel()) then

        tmax=[maxval(tfield%val,dim=3)]
        tmin=[minval(tfield%val,dim=3)]

        do i=1,size(tmax)
            call allmax(tmax(i))
            call allmin(tmin(i))
        end do

        is_constant=maxval(abs(tmax-tmin))<tol

    else
        is_constant=.true.
        do i=2,size(tfield%val,3)
            if (any(abs(tfield%val(:,:,i)-tfield%val(:,:,1))>tol)) then
                is_constant=.false.
                exit
            end if
        end do
    end if
end function is_constant

function GetOldName(tfield) result(old_name)

    type(tensor_field), intent(in) :: tfield
    character (len=FIELD_NAME_LEN) :: old_name

    if(tfield%name(1:6)=='Packed') then
        old_name='PackedOld'//tfield%name(7:)
    else
        old_name='Old'//tfield%name
    end if

end function GetOldName

function GetFEMName(tfield) result(fem_name)

    type(tensor_field), intent(in) :: tfield
    character (len=FIELD_NAME_LEN) :: fem_name

    if(tfield%name(1:6)=='Packed') then
        fem_name='PackedFE'//tfield%name(7:)
    else
        fem_name='FEM'//tfield%name
    end if

end function GetFEMName

subroutine calculate_internal_volume(packed_state, Mdims, mass_ele, calculate_mass, &
    cv_ndgln, eles_with_pipe)

    implicit none

    ! Subroutine to calculate the integrated mass inside the domain

    ! Input/output variables
    type(state_type), intent(inout) :: packed_state
    type(multi_dimensions), intent(in) :: Mdims
    real, dimension( : ), intent(in) :: mass_ele ! volume of the element, split into cv_nloc equally sized pieces (barycenter)
    real, dimension(:), intent(inout) :: calculate_mass
    integer, dimension(:), intent( in ) ::  cv_ndgln
    type(pipe_coords), dimension(:), optional, intent(in):: eles_with_pipe
    ! Local variables
    type (tensor_field), pointer :: saturation
    type (vector_field), pointer :: porosity
    integer  :: cv_knod
    integer :: cv_iloc
    integer :: ele
    integer :: i, k

    saturation => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
    ! Extract the Porosity
    porosity => extract_vector_field( packed_state, "Porosity" )

    ! Having extracted the saturation field (phase volume fraction) in cv_adv_diff at control volume nodes, need to calculate it at quadrature points gi.
    ! (Note saturation is defined on a control volume basis and so the field is stored at control volume nodes).
    ! Since the CV shape functions are 1 or 0, the value at Gauss point gi, is given by the value at the nearest CV_node. So
    ! we pass down the value of cv_iloc from cv_adv_diff and calculate cv_knod. This will be the node corresponding to a given
    ! value of gi in the gcount loop in cv_adv_diff. This then gives the value of phaseVG that we need to associate to that particular Gauss point.
    ! Similar calculation done for density.
    if (present(eles_with_pipe)) then
        !Calculate mass within pipes
        DO k = 1, size(eles_with_pipe)
            ELE = eles_with_pipe(k)%ele!Element with pipe
            if (element_owned(saturation, ELE)) then
                DO CV_ILOC =1, mdims%cv_nloc
                    cv_knod=cv_ndgln((ele-1)*mdims%cv_nloc+cv_iloc)
                    !     Porosity constant element-wise so simply extract that value associated to a given element ele
                    do i = Mdims%n_in_pres + 1, Mdims%nphase
                        calculate_mass(i) = calculate_mass(i) + (Mass_ELE(cv_knod)) *saturation%val(1, i,cv_knod)
                    enddo
                ENDDO
            end if
        end do
    else
        !Calculate mass in the reservoir
        Do ELE=1, mdims%TOTELE
            if (element_owned(saturation, ELE)) then
                DO CV_ILOC =1, mdims%cv_nloc
                    cv_knod=cv_ndgln((ele-1)*mdims%cv_nloc+cv_iloc)
                    !     Porosity constant element-wise so simply extract that value associated to a given element ele
                    do i = 1, size(calculate_mass)
                        calculate_mass(i) = calculate_mass(i) + (Mass_ELE(ele)/mdims%cv_nloc)*&
                            porosity%val(1, ele)*saturation%val(1, i,cv_knod)
                    enddo
                ENDDO
            end if
        ENDDO
    end if
    return

end subroutine calculate_internal_volume


logical function have_option_for_any_phase(path, nphase)
    !The path must be the part of the path inside the phase, i.e. /multiphase_properties/capillary_pressure
    implicit none
    character (len=*), intent(in) :: path
    integer, intent(in) :: nphase
    !Local variables
    integer :: iphase

    have_option_for_any_phase = .false.
    do iphase = 1, nphase
        have_option_for_any_phase = have_option_for_any_phase .or. have_option( '/material_phase['// int2str( iphase -1 ) //']/'&
            //trim(path) )
    end do


end function have_option_for_any_phase


!!$ This subroutine calculates the actual Darcy velocity
subroutine get_DarcyVelocity(Mdims, ndgln, packed_state, PorousMedia_absorp)

    implicit none
    type(multi_ndgln), intent(in) :: ndgln
    type(multi_dimensions), intent(in) :: Mdims
    type(state_type), intent(in) :: packed_state
    type (multi_field), intent(in) :: PorousMedia_absorp

    ! Local variables
    type(tensor_field), pointer :: darcy_velocity, velocity, saturation
    real, dimension(:,:), allocatable :: loc_absorp_matrix
    real, dimension(:), allocatable :: sat_weight_velocity
    integer :: cv_iloc, u_iloc, ele, iphase, imat, u_inod, cv_loc, idim
    ! Initialisation
    darcy_velocity => extract_tensor_field(packed_state,"PackedDarcyVelocity")
    velocity => extract_tensor_field(packed_state,"PackedVelocity")
    saturation => extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")

    call zero(darcy_velocity)
    allocate(loc_absorp_matrix(Mdims%nphase*Mdims%ndim,Mdims%nphase*Mdims%ndim))
    allocate(sat_weight_velocity(Mdims%ndim))

    ! Calculation
    do ele = 1, Mdims%totele
        do u_iloc = 1, Mdims%u_nloc
            u_inod = ndgln%u((ele-1)*Mdims%u_nloc+u_iloc)
            do cv_iloc = 1, Mdims%cv_nloc
                imat = ndgln%mat((ele-1)*Mdims%mat_nloc+cv_iloc)
                cv_loc = ndgln%cv((ele-1)*Mdims%cv_nloc+cv_iloc)
                !This is not optimal, maybe just perform when CVN(U_ILOC, CV_INOD) =/ 0
                call get_multi_field_inverse(PorousMedia_absorp, imat, loc_absorp_matrix)
                do iphase = 1, Mdims%nphase
                    !Inverse of sigma avoiding inversion
                    !loc_absorp_matrix(ndim*(iphase-1)+1:iphase*ndim,ndim*(iphase-1)+1:iphase*ndim) = matmul(permeability%val(:,:,ele),&
                    !    absorption_term%val(imat,ndim*(iphase-1)+1:iphase*ndim,ndim*(iphase-1)+1:iphase*ndim))
                    !All the elements should be equal in the diagonal now
                    !loc_absorp_matrix(ndim*(iphase-1)+1:iphase*ndim,ndim*(iphase-1)+1:iphase*ndim) = &
                    !    permeability%val(ndim*(iphase-1)+1:iphase*ndim,ndim*(iphase-1)+1:iphase*ndim,ele)/&
                    !    loc_absorp_matrix(ndim*(iphase-1)+1,ndim*(iphase-1)+1)!Inverse

                    sat_weight_velocity = matmul(loc_absorp_matrix((iphase-1)*Mdims%ndim+1:iphase*Mdims%ndim, &
                        (iphase-1)*Mdims%ndim+1:iphase*Mdims%ndim),velocity%val(:,iphase,u_inod))
                    !P0 darcy velocities per element
                    darcy_velocity%val(:,iphase,u_inod) = darcy_velocity%val(:,iphase,u_inod)+ &
                        sat_weight_velocity(:)*saturation%val(1,iphase,cv_loc)/real(Mdims%cv_nloc)
                end do
            end do
        end do
    end do
    call halo_update(darcy_velocity)
    ! Deallocation
    deallocate(loc_absorp_matrix)
    deallocate(sat_weight_velocity)

end subroutine get_DarcyVelocity

    subroutine Get_Scalar_SNdgln( sndgln, field  )
      implicit none
      type( scalar_field ), intent( in ) :: field
      integer, dimension( : ), intent( inout ) :: sndgln
      ! Local variables
      integer, dimension( : ), allocatable :: snloc
      integer :: sele, iloc

      allocate( snloc( face_loc( field, 1 ) ) )
      do sele = 1, surface_element_count( field )
         snloc = face_global_nodes( field, sele )
         do iloc = 1, face_loc( field, sele )
            sndgln( ( sele - 1 ) * face_loc( field, sele ) + iloc ) =  snloc( iloc )
         end do
      end do

      deallocate( snloc )

      return
    end subroutine Get_Scalar_SNdgln

    subroutine Get_Vector_SNdgln( sndgln, field  )
      implicit none
      type( vector_field ), intent( in ) :: field
      integer, dimension( : ), intent( inout ) :: sndgln
      ! Local variables
      integer, dimension( : ), allocatable :: snloc
      integer :: sele, iloc

      allocate( snloc( face_loc( field, 1 ) ) )
      do sele = 1, surface_element_count( field )
         snloc = face_global_nodes( field, sele )
         do iloc = 1, face_loc( field, sele )
            sndgln( ( sele - 1 ) * face_loc( field, sele ) + iloc ) =  snloc( iloc )
         end do
      end do

      deallocate( snloc )

      return
    end subroutine Get_Vector_SNdgln

    subroutine dump_outflux(current_time, itime, outfluxes)

        ! Subroutine that dumps the total flux at a given timestep across all specified boundaries to a file  called 'simulation_name_outfluxes.csv'. In addition, the time integrated flux
        ! up to the current timestep is also outputted to this file. Integration boundaries are specified in diamond via surface_ids.
        ! (In diamond this option can be found under "/io/dump_boundaryflux/surface_ids" and the user should specify an integer array containing the IDs of every boundary they
        !wish to integrate over).
        real,intent(in) :: current_time
        integer, intent(in) :: itime
        type (multi_outfluxes), intent(inout) :: outfluxes
        !Local variables
        integer :: ioutlet
        integer :: counter
        character (len=1000000) :: whole_line
        character (len=1000000) :: numbers
        integer :: iphase
        ! Strictly speaking don't need character arrays for fluxstring and intfluxstring, could just overwrite each time (may change later)
        character (len = 1000000), dimension(size(outfluxes%intflux,1)) :: fluxstring
        character (len = 1000000), dimension(size(outfluxes%intflux,1)) :: intfluxstring
        character (len = 1000000), dimension(size(outfluxes%intflux,1)) :: tempstring
        character (len = 50) :: simulation_name

        call get_option('/simulation_name', simulation_name)

        if (itime == 1) then
            !The first time, remove file if already exists
            open(unit=89, file=trim(simulation_name)//"_outfluxes.csv", status="replace", action="write")
        else
            open(unit=89, file=trim(simulation_name)//"_outfluxes.csv", action="write", position="append")
        end if


        ! If calculating boundary fluxes, add up contributions to \int{totout} at each time step
        where (outfluxes%totout /= outfluxes%totout)
            outfluxes%totout = 0.!If nan then make it zero
        end where
        outfluxes%intflux = outfluxes%intflux + outfluxes%totout(1, :, :)*dt

        ! Write column headings to file
        counter = 0
        if(itime.eq.1) then
            write(whole_line,*) "Current Time (s)" // "," // "Current Time (years)" // "," // "Pore Volume"
            whole_line = trim(whole_line)
            do ioutlet =1, size(outfluxes%intflux,2)
                do iphase = 1, size(outfluxes%intflux,1)
                    write(fluxstring(iphase),'(a, i0, a, i0, a)') "Phase", iphase, "-S", outfluxes%outlet_id(ioutlet), "- Volume rate"
                    whole_line = trim(whole_line) //","// trim(fluxstring(iphase))
                enddo
                do iphase = 1, size(outfluxes%intflux,1)
                    write(intfluxstring(iphase),'(a, i0, a, i0, a)') "Phase", iphase,  "-S", outfluxes%outlet_id(ioutlet),  "- Cumulative production"
                    whole_line = trim(whole_line) //","// trim(intfluxstring(iphase))
                enddo
                if (has_temperature) then
                    do iphase = 1, size(outfluxes%intflux,1)
                        write(tempstring(iphase),'(a, i0, a, i0, a)') "Phase", iphase,  "-S", outfluxes%outlet_id(ioutlet),  "- Maximum temperature"
                        whole_line = trim(whole_line) //","// trim(tempstring(iphase))
                    enddo
                end if
            end do
             ! Write out the line
            write(89,*), trim(whole_line)
        endif
        ! Write the actual numbers to the file now
        write(numbers,'(E17.11,a,E17.11, a, E17.11)') current_time, "," , current_time/(86400.*365.) , ",",  outfluxes%porevolume
        whole_line =  trim(numbers)
        do ioutlet =1, size(outfluxes%intflux,2)
            do iphase = 1, size(outfluxes%intflux,1)
                write(fluxstring(iphase),'(E17.11)') outfluxes%totout(1, iphase,ioutlet)
                whole_line = trim(whole_line) //","// trim(fluxstring(iphase))
            enddo
            do iphase = 1, size(outfluxes%intflux,1)
                write(intfluxstring(iphase),'(E17.11)') outfluxes%intflux(iphase,ioutlet)
                whole_line = trim(whole_line) //","// trim(intfluxstring(iphase))
            enddo
            if (has_temperature) then
                do iphase = 1, size(outfluxes%intflux,1)
                    write(tempstring(iphase),'(E17.11)') outfluxes%totout(2, iphase,ioutlet)
                    whole_line = trim(whole_line) //","// trim(tempstring(iphase))
                enddo
            end if
        end do
        ! Write out the line
        write(89,*), trim(whole_line)
        close (89)
    end subroutine dump_outflux


end module Copy_Outof_State




