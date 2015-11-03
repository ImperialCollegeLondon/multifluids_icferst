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
    use global_parameters, only: option_path_len, is_porous_media, dumping_in_sat, is_multifracture, FPI_have_converged
    use diagnostic_fields_wrapper_new
    use element_numbering
    use shape_functions
    use fefields
    use boundary_conditions
    use futils, only: int2str
    use boundary_conditions_from_options
    use parallel_tools, only : allmax, allmin, isparallel
    use memory_diagnostics
    use initialise_fields_module, only: initialise_field_over_regions
    use halos
    ! Need to use surface integrals since a function from here is called in the calculate_outflux() subroutine
    use surface_integrals

    !use printout
    !use quicksort
    use solvers
    use conservative_interpolation_module

    implicit none

    private

    public :: Get_Primary_Scalars, Compute_Node_Global_Numbers, Extracting_MeshDependentFields_From_State, &
         Extract_TensorFields_Outof_State, Extract_Position_Field, Get_Ele_Type, Get_Discretisation_Options, &
         print_from_state, update_boundary_conditions, pack_multistate, finalise_multistate, get_ndglno, Adaptive_NonLinear,&
         get_var_from_packed_state, as_vector, as_packed_vector, is_constant, GetOldName, GetFEMName, PrintMatrix, Clean_Storage,&
         CheckElementAngles, bad_elements, calculate_outflux, outlet_id, have_option_for_any_phase, get_regionIDs2nodes,&
         get_Convergence_Functional, get_DarcyVelocity, printCSRMatrix


!    interface Get_Ndgln
!       module procedure Get_Scalar_Ndgln, Get_Vector_Ndgln, Get_Mesh_Ndgln
!    end interface Get_Ndgln

    interface Get_SNdgln
       module procedure Get_Scalar_SNdgln, Get_Vector_SNdgln
    end interface Get_SNdgln

    !This structure is to store data associated with a bad element
    !We store the corresponding element, the node that is in the bad angle
    !and the weights are the values to use to diffuse the bad node to the other nodes
    type bad_elements
       integer :: ele
       integer, allocatable, dimension(:) :: nodes
       real, allocatable, dimension(:) :: weights
       real :: angle
    end type bad_elements

    ! Used in calculations of the outflux - array of integers containing the gmesh IDs of each boundary that you wish to integrate over
    integer, dimension(:), allocatable :: outlet_id

    integer, dimension(2) :: shape

  contains


    subroutine Get_Primary_Scalars( state, &
         nphase, nstate, ncomp, totele, ndim, stotel, &
         u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
         x_snloc, cv_snloc, u_snloc, p_snloc, &
         cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, ph_nloc, ph_nonods )
!!$ This subroutine extracts all primary variables associated with the mesh from state,
!!$ and associated them with the variables used in the MultiFluids model.
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      integer, intent( inout ) :: nphase, nstate, ncomp, totele, ndim, stotel
      integer, intent( inout ), optional :: u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, &
           mat_nloc, x_snloc, cv_snloc, u_snloc, p_snloc, cv_nonods, mat_nonods, u_nonods, &
           xu_nonods, x_nonods, x_nonods_p1, p_nonods, ph_nloc, ph_nonods

!!$ Local variables
      type( vector_field ), pointer :: positions, velocity
      type( scalar_field ), pointer :: pressure
      type( mesh_type ), pointer :: velocity_cg_mesh, pressure_cg_mesh, ph_mesh
      integer :: i, stat

      ewrite(3,*)' In Get_Primary_Scalars'

      ! Read in the surface IDs of the boundaries (if any) that you wish to integrate over into the (integer vector) variable outlet_id.
      ! No need to explicitly allocate outlet_id (done here internally)

      if (have_option( "/io/dump_boundaryflux/surface_ids") .and..not.(allocated(outlet_id))) then
        shape = option_shape("/io/dump_boundaryflux/surface_ids")
        assert(shape(1) >= 0)
        allocate(outlet_id(shape(1)))
        !allocate(outlet_id(1))
        call get_option( "/io/dump_boundaryflux/surface_ids", outlet_id)

      endif

!!$ Defining dimension and nstate
      call get_option( '/geometry/dimension', ndim )
      nstate = option_count( '/material_phase' )

!!$ Assume there are the same number of components in each phase (will need to check this eventually)
      ncomp = 0
      do i = 1, nstate
         if( have_option( '/material_phase[' // int2str(i-1) // &
              ']/is_multiphase_component' ) ) then
            ncomp = ncomp + 1
         end if
      end do
      nphase = nstate - ncomp
      assert( nphase > 0 ) ! Check if there is more than 0 phases

!!$ Get the vel element type.
      is_porous_media = have_option('/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/Porous_media')
      if (is_porous_media) then!Check that the FPI method is on
        if (.not. have_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration')) then
            ewrite(0,*) "WARNING: The option <Fixed_Point_Iteration> is HIGHLY recommended for multiphase porous media flow"
        end if
      end if
      is_multifracture = have_option( '/femdem_fracture' )

      positions => extract_vector_field( state, 'Coordinate' )
      pressure_cg_mesh => extract_mesh( state, 'PressureMesh_Continuous' )

!!$ Defining number of elements and surface elements, coordinates, locs and snlocs
      totele = ele_count( positions )
      stotel = surface_element_count( positions )

!!$ Coordinates
      if( present( x_nloc_p1 ) ) x_nloc_p1 = ele_loc( positions, 1 )
      if( present( x_nloc ) ) x_nloc = ele_loc( pressure_cg_mesh, 1 )
      if( present( x_snloc ) ) x_snloc = face_loc( pressure_cg_mesh, 1 )
      if( present( x_nonods_p1 ) ) x_nonods_p1 = node_count( positions )
      if( present( x_nonods ) ) x_nonods = node_count( pressure_cg_mesh )

!!$ Pressure, Control Volumes and Materials
      pressure => extract_scalar_field( state, 'Pressure' )
      if( present( p_nloc ) ) p_nloc = ele_loc( pressure, 1 )
      if( present( p_snloc ) ) p_snloc = face_loc( pressure, 1 )
      if( present( p_nonods ) ) p_nonods = node_count( pressure )
      if( present( cv_nloc ) ) cv_nloc = p_nloc
      if( present( cv_snloc ) ) cv_snloc = p_snloc
      if( present( cv_nonods ) ) cv_nonods = node_count( pressure )
      if( present( mat_nloc ) ) mat_nloc = cv_nloc
      if( present( mat_nonods ) ) mat_nonods = mat_nloc * totele

!!$ Velocities and velocities (DG) associated with the continuous space (CG)
      velocity => extract_vector_field( state, 'Velocity' )
      if( present( u_nloc ) ) u_nloc = ele_loc( velocity, 1 )
      if( present( u_snloc ) ) u_snloc = face_loc( velocity, 1 )
      if( present( u_nonods ) ) u_nonods = node_count( velocity )

!!$ Get the continuous space of the velocity field
      velocity_cg_mesh => extract_mesh( state, 'VelocityMesh_Continuous' )
      if( present( xu_nloc ) ) xu_nloc = ele_loc( velocity_cg_mesh, 1 )
      if( present( xu_nonods ) ) xu_nonods = max(( xu_nloc - 1 ) * totele + 1, totele )

      if( present( ph_nloc ) ) then
         ph_mesh => extract_mesh( state( 1 ), 'ph', stat )
         if ( stat == 0 ) then
            ph_nloc = ele_loc( ph_mesh, 1 )
            ph_nonods = node_count( ph_mesh )
         else
            ph_nloc = 0
            ph_nonods = 0
         end if
      end if

      return
    end subroutine Get_Primary_Scalars

    function get_ndglno(mesh) result(ndglno)
      type(mesh_type) :: mesh
      integer, dimension(:), pointer  ::  ndglno

      ndglno=> mesh%ndglno
    end function get_ndglno


    subroutine Compute_Node_Global_Numbers( state, &
         totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
         cv_snloc, p_snloc, u_snloc, &
         cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln, &
         cv_sndgln, p_sndgln, u_sndgln )
!!$ This subroutine calculates the global node numbers requested to operates in the MP-space.
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      type( vector_field ), pointer :: positions, velocity
      type( mesh_type ), pointer :: pressure_cg_mesh, velocity_cg_mesh
      type( scalar_field ), pointer :: pressure
      integer, intent( in ) :: totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
           cv_snloc, p_snloc, u_snloc
      integer, dimension( : ), pointer  :: cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln
      integer, dimension( : ) ::  cv_sndgln, p_sndgln, u_sndgln
!!$ Local variables

      x_ndgln_p1=>get_ndglno(extract_mesh(state(1),"CoordinateMesh"))
      x_ndgln=>get_ndglno(extract_mesh(state(1),"PressureMesh_Continuous"))
      cv_ndgln=>get_ndglno(extract_mesh(state(1),"PressureMesh"))
      p_ndgln=>get_ndglno(extract_mesh(state(1),"PressureMesh"))
      mat_ndgln=>get_ndglno(extract_mesh(state(1),"PressureMesh_Discontinuous"))
      u_ndgln=>get_ndglno(extract_mesh(state(1),"InternalVelocityMesh"))
      xu_ndgln=>get_ndglno(extract_mesh(state(1),"VelocityMesh_Continuous"))

!!$ Linear mesh coordinate
      positions => extract_vector_field( state( 1 ), 'Coordinate' )
!!$      call Get_Ndgln( x_ndgln_p1, positions )
!!$
!!$ Positions/Coordinates
      pressure_cg_mesh => extract_mesh( state( 1 ), 'PressureMesh_Continuous' )
!!$      call Get_Ndgln( x_ndgln, pressure_cg_mesh )
!!$
!!$ Pressure, control volume and material
      pressure => extract_scalar_field( state( 1 ), 'Pressure' )
!!$      call Get_Ndgln( cv_ndgln, pressure )
!!$      p_ndgln = cv_ndgln
!!$      mat_ndgln = (/ (i, i = 1, totele * cv_nloc ) /)
!!$
!!$ Velocities
      velocity => extract_vector_field( state( 1 ), 'Velocity' )
!!$      call Get_Ndgln( u_ndgln, velocity, cv_nloc )
      !$
!!$ Velocity in the continuous space
      velocity_cg_mesh => extract_mesh( state( 1 ), 'VelocityMesh_Continuous' )
!!$      call Get_Ndgln( xu_ndgln, velocity_cg_mesh )

!!$ Surface-based global node numbers for control volumes and pressure
      call Get_SNdgln( cv_sndgln, pressure )
      p_sndgln = cv_sndgln

!!$ Velocities

      call Get_SNdgln( u_sndgln, velocity )

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



    subroutine Get_Discretisation_Options( state, &
         t_disopt, v_disopt, t_beta, v_beta, t_theta, v_theta, u_theta, &
         t_dg_vel_int_opt, u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
         comp_diffusion_opt, ncomp_diff_coef, in_ele_upwind, dg_ele_upwind, &
         nits_flux_lim_t, nits_flux_lim_volfra, nits_flux_lim_comp, &
         volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux, &
         t_use_theta_flux, t_get_theta_flux, scale_momentum_by_volume_fraction )
!!$ This subroutine extract all discretisation options from the schema
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      integer, intent( inout ) :: t_disopt, v_disopt
      real, intent( inout ) :: t_beta, v_beta, t_theta, v_theta, u_theta
      integer, intent( inout ) :: t_dg_vel_int_opt, u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
           comp_diffusion_opt, ncomp_diff_coef, in_ele_upwind, dg_ele_upwind, &
           nits_flux_lim_t, nits_flux_lim_volfra, nits_flux_lim_comp
      logical, intent( inout ) :: volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, &
           comp_get_theta_flux, t_use_theta_flux, t_get_theta_flux, scale_momentum_by_volume_fraction

!!$ Local variables:
      integer :: nphase, nstate, ncomp, totele, ndim, stotel, iphase
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

!!$ Extracting primary scalars as local variables:
      call Get_Primary_Scalars( state, &
           nphase, nstate, ncomp, totele, ndim, stotel )

!!$ Solving Advection Field: Temperature
      option_path = '/material_phase[0]/scalar_field::Temperature'
      option_path2 = trim( option_path ) //  '/prognostic/spatial_discretisation'
      option_path3 = trim( option_path ) //  '/prognostic/temporal_discretisation/control_volumes/number_advection_iterations'
      t_disopt = 1

      call get_option( trim( option_path3 ), nits_flux_lim_t, default = 3 )

      Conditional_TDISOPT: if( have_option( trim( option_path2 ) ) ) then
         if( have_option( trim( option_path2 ) // '/control_volumes/face_value::FiniteElement/limit_face_value/' // &
              'limiter::CompressiveAdvection' ) ) then
            t_disopt = 9
         else
            if( have_option( trim( option_path2 ) // '/control_volumes/face_value::FiniteElement/limit_face_value' ) ) &
                 t_disopt = 5
         end if
      end if Conditional_TDISOPT

      call get_option( trim( option_path2 ) // '/conservative_advection', t_beta, default = 0.0 )
      call get_option( '/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation/theta', &
           t_theta, default = 1. )

!!$ Solving Advection Field: Volume fraction
      option_path = '/material_phase[0]/scalar_field::PhaseVolumeFraction'
      option_path2 = trim( option_path ) // '/prognostic/spatial_discretisation/control_volumes/face_value'
      option_path3 = trim( option_path ) // '/prognostic/temporal_discretisation/control_volumes/number_advection_iterations'
      v_disopt = 8

      call get_option( trim( option_path3 ), nits_flux_lim_volfra, default = 3 )

      Conditional_VDISOPT: if( have_option( trim( option_path ) ) ) then
         if( have_option( trim( option_path2 ) // '::FirstOrderUpwind' ) ) v_disopt = 0
         if( have_option( trim( option_path2 ) // '::Trapezoidal' ) ) v_disopt = 2
         if( have_option( trim( option_path2 ) // '::FiniteElement/do_not_limit_face_value' ) ) v_disopt = 6
         if( have_option( trim( option_path2 ) // '::FiniteElement/limit_face_value/limiter::Sweby' ) ) v_disopt = 5
         if( have_option( trim( option_path2 ) // '::FiniteElement/limit_face_value/limiter::CompressiveAdvection' ) ) v_disopt = 9
      end if Conditional_VDISOPT

      call get_option( trim( option_path ) // '/prognostic/spatial_discretisation/conservative_advection', v_beta )
      call get_option( trim( option_path ) // '/prognostic/temporal_discretisation/theta', v_theta )

!!$ Solving Velocity Field
      call get_option( '/material_phase[0]/vector_field::Velocity/prognostic/temporal_discretisation/theta', u_theta )

!!$ Solving Component Field
      option_path3 = '/material_phase[' // int2str( nphase ) // ']/scalar_field::ComponentMassFractionPhase1/' // &
           'temporal_discretisation/control_volumes/number_advection_iterations'
      call get_option( trim( option_path3 ), nits_flux_lim_comp, default = 3 )

!!$ Scaling factor for the momentum equation
      scale_momentum_by_volume_fraction = .false.
      do iphase = 1, nphase
         option_path = '/material_phase[' // int2str( iphase - 1 ) // ']/scale_momentum_by_volume_fraction'
         if( have_option( trim( option_path ) ) ) scale_momentum_by_volume_fraction = .true.
      end do

!!$ Options below are hardcoded and need to be added into the schema
      t_dg_vel_int_opt = 1 ; u_dg_vel_int_opt = 4 ; v_dg_vel_int_opt = 4 ; w_dg_vel_int_opt = 0
      if(is_porous_media) then
        if ( have_option( &
        '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/advection_scheme/DG_weighting') &
        ) v_dg_vel_int_opt = 10
      else
        v_dg_vel_int_opt = 1
      end if
      comp_diffusion_opt = 0 ; ncomp_diff_coef = 0
      volfra_use_theta_flux = .false. ; volfra_get_theta_flux = .true.
      comp_use_theta_flux = .false. ; comp_get_theta_flux = .true.
      t_use_theta_flux = .false. ; t_get_theta_flux = .false.

!!$ IN/DG_ELE_UPWIND are options for optimisation of upwinding across faces in the compact_overlapping
!!$ formulation. The data structure and options for this formulation need to be added later.
      in_ele_upwind = 3 ; dg_ele_upwind = 3

      return
    end subroutine Get_Discretisation_Options



    subroutine Extract_Position_Field( state, &
         xu, yu, zu )
!!$ This subroutine extracts the spatial coordinates fields from state-space and copy them into
!!$ MP-space
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      type( mesh_type ), pointer :: velocity_cg_mesh
      type( vector_field ), pointer :: positions
      type( vector_field ) :: velocity_cg
      real, dimension( : ) :: xu, yu, zu
!!$ Local variables
      integer :: ndim, totele, xu_nloc, xu_nonods, ele, iloc
      integer, dimension( : ), pointer :: element_nodes

      call get_option( '/geometry/dimension', ndim )
      positions => extract_vector_field( state, 'Coordinate' )
      totele = ele_count( positions )

!!$ Velocity in the continuous space
      velocity_cg_mesh => extract_mesh( state( 1 ), 'VelocityMesh_Continuous' )
      call allocate( velocity_cg, ndim, velocity_cg_mesh, 'Velocity_CG_Coordinates' )
      velocity_cg % val( :, : )= 0
      call project_field( positions, velocity_cg, positions )
      xu_nloc = ele_loc( velocity_cg_mesh, 1 )
      xu_nonods = max( ( xu_nloc - 1 ) * totele + 1, totele )

      Loop_Elements: do ele = 1, totele
         element_nodes => ele_nodes( velocity_cg_mesh, ele )
         Loop_Local_Nodes: do iloc = 1, xu_nloc
            xu( element_nodes( iloc ) ) = velocity_cg % val( 1, element_nodes( iloc ) )
            if( ndim > 1 )  yu( element_nodes( iloc ) ) = velocity_cg % val( 2, element_nodes( iloc ) )
            if( ndim > 2 )  zu( element_nodes( iloc ) ) = velocity_cg % val( 3, element_nodes( iloc ) )
         end do Loop_Local_Nodes
      end do Loop_Elements

      call deallocate( velocity_cg )

      return
    end subroutine Extract_Position_Field

!    subroutine xp1_2_xp2( state, &
!         x_nloc_p2, x_nloc_p1, x_nonods_p1, x_nonods_p2, &
!         x_ndgln_p1, x_ndgln_p2, &
!         x, y, z )
!      ! This subrt maps the coordinate P1 mesh into a P2 mesh.
!      implicit none
!      type( state_type ), dimension( : ), intent( in ) :: state
!      type( vector_field ), pointer :: positions
!      integer, intent( in ) :: x_nloc_p2, x_nloc_p1, x_nonods_p1, x_nonods_p2
!      integer, dimension( : ), intent( in ) :: x_ndgln_p1
!      integer, dimension( : ), intent( in ) :: x_ndgln_p2
!      real, dimension( : ), intent( inout ) :: x, y, z
!
!      ! Local variables
!      real, dimension( x_nonods_p1 ) :: x_p1, y_p1, z_p1
!      integer, dimension( x_nloc_p2 ) :: iloclist_p2
!      real, dimension( x_nloc_p2 ) :: x2, y2, z2
!      integer :: ndim, totele, ele, iloc, inod
!      real :: xnod1, xnod2, ynod1, ynod2, xtemp, ytemp
!
!      call get_option( '/geometry/dimension', ndim )
!      positions => extract_vector_field( state( 1 ), 'Coordinate' )
!      totele = ele_count( positions )
!
!        x2 = 0.;y2 = 0.;z2 = 0.
!
!      if( ndim == 2 ) then
!         iloclist_p2 = (/ 1, 4, 2, 5, 6, 3 /)
!      elseif( ndim == 3 ) then
!         iloclist_p2 = (/ 1, 5, 2, 6, 7, 3, 8, 9, 10, 4 /)
!      else
!         iloclist_p2 = (/ 1, 3, 2 /)
!      end if
!
!      Conditional_Pn: if ( x_nloc_p2 == 3 .or. x_nloc_p2 == 4 ) then
!
!         if( ( x_nloc_p2 == 3 ) .and. ( ndim == 1 ) ) then ! 1D quadratic
!            x_p1 = positions % val( 1, : )
!            do ele = 1, totele
!               x2 = 0.
!               do iloc = 1, x_nloc_p1
!                  x2( iloc ) = x_p1( x_ndgln_p1( ( ele - 1 ) * x_nloc_p1 + iloc ))
!               end do
!               do iloc = 1, x_nloc_p1 - 1
!                  xnod1 = x2( iloc )
!                  xnod2 = x2( iloc + 1 )
!                  x2( x_nloc_p1 + iloc ) = 0.5 * (  xnod1  +  xnod2  )
!               end do
!               do iloc = 1, x_nloc_p2
!                  inod = x_ndgln_p2( ( ele - 1 ) * x_nloc_p2 + iloc )
!                  x( inod ) = x2( iloclist_p2( iloc ) )
!               end do
!            end do
!
!         else
!
!            x = positions % val( 1, : )
!            if( ndim > 1 ) y = positions % val( 2, : )
!            if( ndim > 2 ) z = positions % val( 3, : )
!         end if
!
!      else if ( x_nloc_p2 == 6 .or. x_nloc_p2 == 10 ) then
!
!         x_p1 = positions % val( 1, : )
!         y_p1 = positions % val( 2, : )
!         if (ndim == 3)  z_p1 = positions % val( 3, : )
!
!         if ( ( x_nloc_p2 == 6 ) .and. ( ndim == 2 ) ) then ! 2D P2 Tri
!
!            do ele = 1, totele
!
!               x2 = 0. ; y2 = 0.
!               do iloc = 1, x_nloc_p1
!                  x2( iloc ) = x_p1( x_ndgln_p1( ( ele - 1 ) * x_nloc_p1 + iloc ))
!                  y2( iloc ) = y_p1( x_ndgln_p1( ( ele - 1 ) * x_nloc_p1 + iloc ))
!               end do
!
!               do iloc = 1, x_nloc_p1
!                  if( iloc < x_nloc_p1 ) then
!                     xnod1 = x2( iloc )      ; ynod1 = y2( iloc )
!                     xnod2 = x2( iloc + 1 ); ynod2 = y2( iloc + 1 )
!                  else
!                     xnod1 = x2( iloc ) ; ynod1 = y2( iloc )
!                     xnod2 = x2( 1 )    ; ynod2 = y2 ( 1 )
!                  end if
!                  x2( x_nloc_p1 + iloc ) = 0.5 * (  xnod1  +  xnod2  )
!                  y2( x_nloc_p1 + iloc ) = 0.5 * (  ynod1  +  ynod2  )
!               end do
!
!               xtemp = x2( 5 ) ; ytemp = y2( 5 )
!               x2( 5 ) = x2( 6 ) ; y2( 5 ) = y2( 6 )
!               x2( 6 ) = xtemp ; y2( 6 ) = ytemp
!
!               do iloc = 1, x_nloc_p2
!                  inod = x_ndgln_p2( ( ele - 1 ) * x_nloc_p2 + iloc )
!                  x( inod ) = x2( iloclist_p2( iloc ) )
!                  y( inod ) = y2( iloclist_p2( iloc ) )
!               end do
!
!            end do
!
!         else ! Quadratic Tets
!
!            do ele = 1, totele
!
!               x2 = 0. ; y2 = 0. ; z2 = 0.
!               do iloc = 1, x_nloc_p1
!                  x2( iloc ) = x_p1( x_ndgln_p1( ( ele - 1 ) * x_nloc_p1 + iloc ))
!                  y2( iloc ) = y_p1( x_ndgln_p1( ( ele - 1 ) * x_nloc_p1 + iloc ))
!                  z2( iloc ) = z_p1( x_ndgln_p1( ( ele - 1 ) * x_nloc_p1 + iloc ))
!               end do
!
!               x2( 5 ) = 0.5 * (x2(1) + x2(2) )
!               y2( 5 ) = 0.5 * (y2(1) + y2(2) )
!               z2( 5 ) = 0.5 * (z2(1) + z2(2) )
!
!               x2( 6 ) = 0.5 * (x2(1) + x2(3) )
!               y2( 6 ) = 0.5 * (y2(1) + y2(3) )
!               z2( 6 ) = 0.5 * (z2(1) + z2(3) )
!
!               x2( 7 ) = 0.5 * (x2(2) + x2(3) )
!               y2( 7 ) = 0.5 * (y2(2) + y2(3) )
!               z2( 7 ) = 0.5 * (z2(2) + z2(3) )
!
!               x2( 8 ) = 0.5 * (x2(1) + x2(4) )
!               y2( 8 ) = 0.5 * (y2(1) + y2(4) )
!               z2( 8 ) = 0.5 * (z2(1) + z2(4) )
!
!               x2( 9 ) = 0.5 * (x2(2) + x2(4) )
!               y2( 9 ) = 0.5 * (y2(2) + y2(4) )
!               z2( 9 ) = 0.5 * (z2(2) + z2(4) )
!
!               x2( 10 ) = 0.5 * (x2(3) + x2(4) )
!               y2( 10 ) = 0.5 * (y2(3) + y2(4) )
!               z2( 10 ) = 0.5 * (z2(3) + z2(4) )
!
!
!               do iloc = 1, x_nloc_p2
!                  inod = x_ndgln_p2( ( ele - 1 ) * x_nloc_p2 + iloc )
!                  x( inod ) = x2( iloclist_p2( iloc ) )
!                  y( inod ) = y2( iloclist_p2( iloc ) )
!                  z( inod ) = z2( iloclist_p2( iloc ) )
!               end do
!
!            end do
!
!         end if
!      end if Conditional_Pn
!
!
!      return
!    end subroutine xp1_2_xp2


    subroutine Extracting_MeshDependentFields_From_State( state, packed_state, initialised, &
         Velocity_U_Source, Velocity_Absorption, &
         Permeability )
      implicit none
      type( state_type ), dimension( : ), intent( inout ) :: state
      type( state_type ), intent( inout ) :: packed_state

      logical, intent( in ) :: initialised
      real, dimension( :, :, : ), intent( inout ) :: Velocity_U_Source
      real, dimension( :, :, : ), intent( inout ) :: Velocity_Absorption
      real, dimension( :, :, : ), optional, intent( inout ) :: Permeability

!!$ Local variables
      type( scalar_field ), pointer :: scalarfield
      type( vector_field ), pointer :: vectorfield, x_all
      type( tensor_field ), pointer :: tensorfield
      integer, dimension( : ), pointer :: element_nodes
      character( len = option_path_len ) :: option_path
      integer :: nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, &
           cv_ele_type, p_ele_type, u_ele_type, &
           iphase, ele, idim
      integer, dimension( : ), pointer :: cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, &
           xu_ndgln, mat_ndgln, cv_sndgln, p_sndgln, u_sndgln
      logical :: is_symmetric
      real, dimension( : ), allocatable :: dummy!, x, y, z

!!$ Extracting spatial resolution
      call Get_Primary_Scalars( state, &
           nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods )

!!$ Calculating Global Node Numbers
      allocate( cv_sndgln( stotel * cv_snloc ), p_sndgln( stotel * p_snloc ), &
           u_sndgln( stotel * u_snloc ), dummy( cv_nonods ) )

      !      x_ndgln_p1 = 0 ; x_ndgln = 0 ; cv_ndgln = 0 ; p_ndgln = 0 ; mat_ndgln = 0
      !     u_ndgln = 0 ;  xu_ndgln = 0 ;
      cv_sndgln = 0 ; p_sndgln = 0 ; u_sndgln = 0

      call Compute_Node_Global_Numbers( state, &
           totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
           cv_snloc, p_snloc, u_snloc, &
           cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln, &
           cv_sndgln, p_sndgln, u_sndgln )

      call Get_Ele_Type( x_nloc, cv_ele_type, p_ele_type, u_ele_type )


!      allocate( x( x_nonods ) , y( x_nonods ), z( x_nonods ) )
!
!      x = 0. ; y = 0. ; z = 0.
!      call xp1_2_xp2( state, &
!           x_nloc, x_nloc_p1, x_nonods_p1, x_nonods, &
!           x_ndgln_p1, x_ndgln, &
!           x, y, z )
!      do idim = 1, ndim
!         if( idim ==1 ) x_all % val( idim, : ) = x
!         if( idim ==2 ) x_all % val( idim, : ) = y
!         if( idim ==3 ) x_all % val( idim, : ) = z
!      end do

      !Get the coordinates of the nodes from the mesh
      x_all => extract_vector_field( packed_state, "PressureCoordinate" )

!!$
!!$ Extracting Pressure Field:
!!$
      scalarfield => extract_scalar_field( state( 1 ), 'Pressure' )
      call Get_ScalarFields_Outof_State( state, initialised, 1, scalarfield, &
           dummy )

!!$
!!$ Extracting Density Field:
!!$
      Loop_Density: do iphase = 1, nphase
         scalarfield => extract_scalar_field( state( iphase ), 'Density' )
         !knod = ( iphase - 1 ) * node_count( scalarfield )
         call Get_ScalarFields_Outof_State( state, initialised, iphase, scalarfield, &
              dummy)
      end do Loop_Density

!!$
!!$ Extracting Components Field:
!!$
!      Loop_Components: do icomp = nphase + 1, nphase + ncomp ! Component loop
!         Loop_Phases_Components: do iphase = 1, nphase ! Phase loop
!
!            scalarfield => extract_scalar_field( state( icomp ), 'ComponentMassFractionPhase' // int2str( iphase ) )
!
!            knod = ( icomp - ( nphase + 1 ) ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods
!            knod2 = ( icomp - ( nphase + 1 ) ) * nphase * stotel * cv_snloc + &
!                 ( iphase - 1 ) * stotel * cv_snloc
!
!            call Get_CompositionFields_Outof_State( state, initialised, nphase, icomp, iphase, scalarfield, &
!                 Component( knod + 1 : knod + cv_nonods ), &
!                 field_prot_source = Component_Source( ( iphase - 1 ) * cv_nonods + 1 : &
!                 ( iphase - 1 ) * cv_nonods + cv_nonods ) )
!
!         end do Loop_Phases_Components
!      end do Loop_Components

!!$
!!$ Extracting Velocity Field:
!!$
      Loop_Velocity: do iphase = 1, nphase
         vectorfield => extract_vector_field( state( iphase ), 'Velocity' )
         call Get_VectorFields_Outof_State( state, initialised, iphase, vectorfield, &
              field_prot_source=Velocity_U_Source, field_prot_absorption=Velocity_Absorption )
      end do Loop_Velocity

!!$
!!$ Extracting Temperature Field:
!!$
!      do iphase = 1, nphase
!         Conditional_Temperature: if( have_option( '/material_phase[' // int2str( iphase - 1 ) // &
!              ']/scalar_field::Temperature' ) ) then
!            scalarfield => extract_scalar_field( state( iphase ), 'Temperature' )
!            knod = ( iphase - 1 ) * node_count( scalarfield )
!            !call Get_ScalarFields_Outof_State( state, initialised, iphase, scalarfield, &
!            !     Temperature( knod + 1 : knod + node_count( scalarfield ) ), &
!            !     field_prot_source = Temperature_Source( knod + 1 : knod + node_count( scalarfield ) ) )
!            call Get_ScalarFields_Outof_State( state, initialised, iphase, scalarfield, &
!                 Temperature( iphase, : ), &
!                 field_prot_source = Temperature_Source( knod + 1 : knod + node_count( scalarfield ) ) )
!         end if Conditional_Temperature
!      end do

!!$
!!$ Extracting Permeability Field:
!!$
      if (present(Permeability)) then
          Permeability = 0.
          Conditional_PermeabilityField: if( have_option( '/porous_media/scalar_field::Permeability' ) ) then

             scalarfield => extract_scalar_field( state( 1 ), 'Permeability' )
             do ele = 1, element_count( scalarfield )
                element_nodes => ele_nodes( scalarfield, ele )
                forall( idim = 1 : ndim ) Permeability( ele, idim, idim ) = scalarfield % val( element_nodes( 1 ) )
             end do

          elseif( have_option( '/porous_media/tensor_field::Permeability' ) ) then

             tensorfield => extract_tensor_field( state( 1 ), 'Permeability' )
             option_path =  '/porous_media/tensor_field::Permeability'
             call Extract_TensorFields_Outof_State( state, 1, &
                  tensorfield, option_path, &
                  Permeability )

          elseif( have_option( '/porous_media/vector_field::Permeability' ) ) then
             FLAbort( 'Permeability Vector Field is not defined yet.' )

          end if Conditional_PermeabilityField
      end if
      deallocate( cv_sndgln, p_sndgln, u_sndgln, dummy )

      return
    end subroutine Extracting_MeshDependentFields_From_State




 subroutine Get_ScalarFields_Outof_State( state, initialised, iphase, field, &
         field_prot, wic_bc, suf_bc, field_prot_source, field_prot_absorption, suf_bc_rob1, suf_bc_rob2 )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      logical, intent( in ) :: initialised
      integer, intent( in ) :: iphase
      type( scalar_field ), pointer :: field, field_prot_bc
      real, dimension( : ), intent( inout ) :: field_prot
      real, dimension( : ), intent( inout ), optional :: field_prot_source,field_prot_absorption, suf_bc_rob1, suf_bc_rob2
      integer, dimension( : ), intent( inout ), optional :: wic_bc
      real, dimension( : ), intent( inout ), optional :: suf_bc

      ! Local variables
      type( mesh_type ), pointer :: pmesh, cmesh
      type(scalar_field), pointer :: pressure, field_source, field_absorption
      type(scalar_field) :: dummy
      type(vector_field), pointer :: positions
      integer, dimension(:), allocatable :: sufid_bc
      character( len = option_path_len ) :: option_path, option_path2, field_name
      integer :: stotel, nobcs, bc_type, i, j, k, kk, sele
      integer :: nstate, nphase, ncomp, snloc, stat
      integer :: shape_option(2)
      real :: initial_constant
      logical :: have_source, have_absorption
      integer, dimension(:), allocatable :: face_nodes
      character( len = 8192 ) :: func

      field_name = trim( field % name )
      have_source = .false. ; have_absorption = .false.

      Conditional_SourceField: if( present( field_prot_source ) ) then
         field_source => extract_scalar_field( state( iphase ), trim(field_name) // 'Source', stat )
         have_source = ( stat == 0 )

         if ( have_source ) then
            do j = 1, node_count( field_source )
               field_prot_source( ( iphase - 1 ) * node_count( field_source ) + j ) = &
                    field_source % val( j )
            end do
         end if
      end if Conditional_SourceField

      Conditional_AbsorptionField: if( present( field_prot_absorption ) ) then
         field_absorption => extract_scalar_field( state( iphase ), trim(field_name) // 'Absorption', stat )
         have_absorption = ( stat == 0 )
         if ( have_absorption ) then
            do j = 1, node_count( field_absorption )
               field_prot_absorption( ( iphase - 1 ) * node_count( field_absorption ) + j ) = &
                    field_absorption % val( j )
            end do
         end if
      end if Conditional_AbsorptionField

      pressure => extract_scalar_field( state( 1 ), 'Pressure' )
      pmesh => extract_mesh( state, 'PressureMesh' )
      cmesh => extract_mesh( state, 'CoordinateMesh' )
      positions => extract_vector_field( state( 1 ), 'Coordinate' )

      snloc = face_loc( pressure, 1 ) ; stotel = surface_element_count( cmesh ) ; &
           nstate = option_count( '/material_phase' )

      ncomp = 0
      do i = 1, nstate
         if( have_option( '/material_phase[' // int2str( i - 1) // ']/is_multiphase_component' ) )then
            ncomp = ncomp + 1
         end if
      end do
      nphase = nstate - ncomp
      option_path = '/material_phase['//int2str( iphase - 1 )//']/scalar_field::'//trim( field_name )
      ewrite(3,*)'option_path:', trim( option_path )


!!$ This will need to be ammended later on to take into account python functions that impose
!!$ time-dependent field changes
      Conditional_InitialisationFromFLML: if( initialised ) then ! Extracting from state after initialisation
         field_prot = field % val
!!$         field_prot( 1 : node_count( field ) ) = field % val
!!$         field_prot( ( iphase - 1 ) * node_count( field ) + 1 : iphase * node_count( field ) ) = &
!!$              field % val

      else !Initialisation before adapt
         if( have_option( trim( option_path ) // '/prognostic/initial_condition::WholeMesh/constant' ) )then
            call get_option(trim( option_path ) // '/prognostic/initial_condition::WholeMesh/constant', &
                 initial_constant )
            field_prot = initial_constant

         elseif( have_option( trim( option_path ) // '/prognostic/initial_condition::WholeMesh/python ') )then
            call get_option( trim( option_path ) // '/prognostic/initial_condition::WholeMesh/python', func )
            call allocate( dummy, field % mesh, 'dummy' )
            call get_option('/timestepping/current_time', current_time)
            call set_from_python_function(dummy, trim(func), positions, current_time)
            field_prot = dummy % val
            call deallocate( dummy )
         elseif( have_option( trim( option_path ) // '/prognostic/initial_condition/from_file')) then
            field_prot = field % val

         else if (have_option( trim( option_path ) // '/prognostic/initial_condition') )then
            call allocate( dummy, field % mesh, 'dummy' )
            call get_option('/timestepping/current_time', current_time)
            call initialise_field_over_regions(dummy, trim( option_path ) // '/prognostic/initial_condition', positions, current_time)
            field_prot = dummy%val
            call deallocate( dummy )

         else
            !ewrite(-1,*) 'No initial condition for field::', trim( field_name )
            !FLAbort( 'Check initial conditions' )
         end if
      end if Conditional_InitialisationFromFLML

!!$ Boundary conditions
      if (present( wic_bc ) ) then

         option_path2 = trim( option_path ) // '/prognostic/boundary_conditions['
         nobcs = get_boundary_condition_count( field )
         Loop_BC: do k = 1, nobcs

            option_path = trim( option_path2 ) // int2str( k - 1 ) // ']/surface_ids'
            shape_option = option_shape( trim( option_path ) )
            allocate( SufID_BC( 1 : shape_option( 1 ) ) )
            call get_option( trim( option_path ), SufID_BC )
            allocate( face_nodes( face_loc( field, 1) ) )

            option_path = trim( option_path2 ) // int2str( k - 1 ) // ']/'

            Conditional_Field_BC: if( have_option( trim( option_path ) // 'type::dirichlet' ) ) then

               BC_Type = 1
               field_prot_bc => extract_surface_field( field, k, 'value' )

               sele = 1
               do j = 1, stotel
                  if( any ( SufID_BC == pmesh % faces % boundary_ids( j ) ) ) then
                     wic_bc( j + ( iphase - 1 ) * stotel ) = BC_Type
                     face_nodes = ele_nodes( field_prot_bc, sele )
                     do kk = 1, snloc
                        suf_bc( ( iphase - 1 ) * stotel * snloc + ( j - 1 ) * snloc + kk ) = &
                             field_prot_bc % val( face_nodes( kk ) )
                     end do
                     sele = sele + 1
                  end if
               end do

            else if( have_option( trim( option_path ) // 'type::robin' ) ) then

               BC_Type = 2

               do j = 1, stotel
                  if( any ( SufID_BC == pmesh % faces % boundary_ids( j ) ) ) then
                     wic_bc( j + ( iphase - 1 ) * stotel ) = BC_Type
                  end if
               end do

               ! calculate this later on...
               suf_bc_rob1 = 0.
               suf_bc_rob2 = 0.

            end if Conditional_Field_BC

            deallocate( face_nodes, sufid_bc )

         end do Loop_BC

      end if

      return
    end subroutine Get_ScalarFields_Outof_State




!subroutine Get_ScalarFields_Outof_State2( state, initialised, iphase, field, &
!         field_prot, wic_bc, suf_bc, field_prot_source, field_prot_absorption, suf_bc_rob1, suf_bc_rob2 )
!      implicit none
!      type( state_type ), dimension( : ), intent( in ) :: state
!      logical, intent( in ) :: initialised
!      integer, intent( in ) :: iphase
!      type( scalar_field ), pointer :: field, field_prot_bc, field_prot_bc1, field_prot_bc2
!      real, dimension( : ), intent( inout ) :: field_prot
!      real, dimension( : ), intent( inout ), optional :: field_prot_absorption, suf_bc_rob1, suf_bc_rob2
!      integer, dimension( : ), intent( inout ), optional :: wic_bc
!      real, dimension( : ), intent( inout ), optional :: suf_bc
!      real, dimension( : , : ), intent( inout ), optional :: field_prot_source
!
!      ! Local variables
!      type( mesh_type ), pointer :: pmesh, cmesh
!      type(scalar_field), pointer :: pressure, field_source, field_absorption
!      type(scalar_field) :: dummy
!      type(vector_field), pointer :: positions
!      integer, dimension(:), allocatable :: sufid_bc
!      character( len = option_path_len ) :: option_path, option_path2, field_name
!      integer :: stotel, nobcs, bc_type, i, j, k, kk, sele
!      integer :: nstate, nphase, ncomp, snloc, stat
!      integer :: shape_option(2)
!      real :: initial_constant
!      logical :: have_source, have_absorption
!      integer, dimension(:), allocatable :: face_nodes
!      character( len = 8192 ) :: func
!
!      field_name = trim( field % name )
!      have_source = .false. ; have_absorption = .false.
!
!      Conditional_SourceField: if( present( field_prot_source ) ) then
!         field_source => extract_scalar_field( state( iphase ), trim(field_name) // 'Source', stat )
!         have_source = ( stat == 0 )
!
!         if ( have_source ) then
!            do j = 1, node_count( field_source )
!               !field_prot_source( ( iphase - 1 ) * node_count( field_source ) + j ) = &
!               !     field_source % val( j )
!               field_prot_source( iphase, j ) = &
!                    field_source % val( j )
!            end do
!         end if
!      end if Conditional_SourceField
!
!      Conditional_AbsorptionField: if( present( field_prot_absorption ) ) then
!         field_absorption => extract_scalar_field( state( iphase ), trim(field_name) // 'Absorption', stat )
!         have_absorption = ( stat == 0 )
!         if ( have_absorption ) then
!            do j = 1, node_count( field_absorption )
!               field_prot_absorption( ( iphase - 1 ) * node_count( field_absorption ) + j ) = &
!                    field_absorption % val( j )
!            end do
!         end if
!      end if Conditional_AbsorptionField
!
!      pressure => extract_scalar_field( state( 1 ), 'Pressure' )
!      pmesh => extract_mesh( state, 'PressureMesh' )
!      cmesh => extract_mesh( state, 'CoordinateMesh' )
!      positions => extract_vector_field( state( 1 ), 'Coordinate' )
!
!      snloc = face_loc( pressure, 1 ) ; stotel = surface_element_count( cmesh ) ; &
!           nstate = option_count( '/material_phase' )
!
!      ncomp = 0
!      do i = 1, nstate
!         if( have_option( '/material_phase[' // int2str( i - 1) // ']/is_multiphase_component' ) )then
!            ncomp = ncomp + 1
!         end if
!      end do
!      nphase = nstate - ncomp
!      option_path = '/material_phase['//int2str( iphase - 1 )//']/scalar_field::'//trim( field_name )
!      ewrite(3,*)'option_path:', trim( option_path )
!
!
!$ This will need to be ammended later on to take into account python functions that impose
!$ time-dependent field changes
!      Conditional_InitialisationFromFLML: if( initialised ) then ! Extracting from state after initialisation
!         field_prot = field % val
!$         field_prot( 1 : node_count( field ) ) = field % val
!$         field_prot( ( iphase - 1 ) * node_count( field ) + 1 : iphase * node_count( field ) ) = &
!$              field % val
!
!      else !Initialisation before adapt
!         if( have_option( trim( option_path ) // '/prognostic/initial_condition::WholeMesh/constant' ) )then
!            call get_option(trim( option_path ) // '/prognostic/initial_condition::WholeMesh/constant', &
!                 initial_constant )
!            field_prot = initial_constant
!
!         elseif( have_option( trim( option_path ) // '/prognostic/initial_condition::WholeMesh/python ') )then
!            call get_option( trim( option_path ) // '/prognostic/initial_condition::WholeMesh/python', func )
!            call allocate( dummy, field % mesh, 'dummy' )
!            call get_option('/timestepping/current_time', current_time)
!            call set_from_python_function(dummy, trim(func), positions, current_time)
!            field_prot = dummy % val
!            call deallocate( dummy )
!         elseif( have_option( trim( option_path ) // '/prognostic/initial_condition/from_file')) then
!            field_prot = field % val
!
!         else if (have_option( trim( option_path ) // '/prognostic/initial_condition') )then
!            call allocate( dummy, field % mesh, 'dummy' )
!            call get_option('/timestepping/current_time', current_time)
!            call initialise_field_over_regions(dummy, trim( option_path ) // '/prognostic/initial_condition', positions, current_time)
!            field_prot = dummy%val
!            call deallocate( dummy )
!
!         else
!            !ewrite(-1,*) 'No initial condition for field::', trim( field_name )
!            !FLAbort( 'Check initial conditions' )
!         end if
!      end if Conditional_InitialisationFromFLML
!
!$ Boundary conditions
!      if (present( wic_bc ) ) then
!
!         option_path2 = trim( option_path ) // '/prognostic/boundary_conditions['
!         nobcs = get_boundary_condition_count( field )
!         Loop_BC: do k = 1, nobcs
!
!            option_path = trim( option_path2 ) // int2str( k - 1 ) // ']/surface_ids'
!            shape_option = option_shape( trim( option_path ) )
!            allocate( SufID_BC( 1 : shape_option( 1 ) ) )
!            call get_option( trim( option_path ), SufID_BC )
!            allocate( face_nodes( face_loc( field, 1) ) )
!
!            option_path = trim( option_path2 ) // int2str( k - 1 ) // ']/'
!
!            Conditional_Field_BC: if( have_option( trim( option_path ) // 'type::dirichlet' ) ) then
!
!               BC_Type = 1
!               field_prot_bc => extract_surface_field( field, k, 'value' )
!
!               sele = 1
!               do j = 1, stotel
!                  if( any ( SufID_BC == pmesh % faces % boundary_ids( j ) ) ) then
!                     wic_bc( j + ( iphase - 1 ) * stotel ) = BC_Type
!                     face_nodes = ele_nodes( field_prot_bc, sele )
!                     do kk = 1, snloc
!                        suf_bc( ( iphase - 1 ) * stotel * snloc + ( j - 1 ) * snloc + kk ) = &
!                             field_prot_bc % val( face_nodes( kk ) )
!                     end do
!                     sele = sele + 1
!                  end if
!               end do
!
!            else if( have_option( trim( option_path ) // 'type::robin' ) ) then
!
!               BC_Type = 2
!
!               do j = 1, stotel
!                  if( any ( SufID_BC == pmesh % faces % boundary_ids( j ) ) ) then
!                     wic_bc( j + ( iphase - 1 ) * stotel ) = BC_Type
!                  end if
!               end do
!
!               ! calculate this later on...
!               suf_bc_rob1 = 0.
!               suf_bc_rob2 = 0.
!
!            end if Conditional_Field_BC
!
!            deallocate( face_nodes, sufid_bc )
!
!         end do Loop_BC
!
!      end if
!
!      return
!    end subroutine Get_ScalarFields_Outof_State2




!    subroutine Get_CompositionFields_Outof_State( state, initialised, nphase, icomp, iphase, field, &
!         field_prot, wic_bc, &
!         kprime, kprime2, &
!         suf_bc, &
!         field_prot_source, field_prot_absorption )
!      implicit none
!      type( state_type ), dimension( : ), intent( in ) :: state
!      logical, intent( in ) :: initialised
!      integer, intent( in ) :: nphase, icomp, iphase
!      type( scalar_field ), pointer :: field, field_prot_bc
!      real, dimension( : ), intent( inout ) :: field_prot
!      real, dimension( : ), intent( inout ), optional :: field_prot_source, field_prot_absorption
!      integer, dimension( : ), intent( inout ), optional :: wic_bc
!      integer, intent( in ), optional  :: kprime, kprime2
!      real, dimension(  :  ), intent( inout ), optional  :: suf_bc
!      ! Local variables
!      type( mesh_type ), pointer :: pmesh, cmesh
!      type(scalar_field), pointer :: pressure, field_source, field_absorption
!      type( scalar_field ) :: dummy
!      type( vector_field ), pointer :: positions
!      integer, dimension( : ), allocatable :: sufid_bc, face_nodes
!      integer :: shape_option( 2 )
!      character( len = option_path_len ) :: option_path, field_name
!      logical :: have_source, have_absorption
!      integer :: nstate, stotel, nobcs, bc_type, i, j, k, kk, sele, stat, snloc
!      real :: initial_constant
!      character( len = 8192 ) :: func
!
!      field_name = trim( field % name )
!      positions => extract_vector_field( state( 1 ), 'Coordinate' )
!      pressure => extract_scalar_field( state( 1 ), 'Pressure' )
!      pmesh => extract_mesh( state, 'PressureMesh' )
!      cmesh => extract_mesh( state, 'CoordinateMesh' )
!
!      field_source => extract_scalar_field( state( iphase ), field_name // 'Source', stat )
!      have_source = ( stat == 0 )
!      field_absorption => extract_scalar_field( state( iphase ), field_name // 'Absorption', stat )
!      have_absorption = ( stat == 0 )
!
!      snloc = face_loc( pressure, 1 )
!      nstate = option_count('/material_phase')
!      stotel = surface_element_count( cmesh )
!
!      option_path = '/material_phase[' // int2str( icomp - 1 ) // &
!           ']/scalar_field::ComponentMassFractionPhase' // &
!           int2str( iphase )
!
!      Conditional_InitialisedFromFLML: if( initialised ) then
!         field_prot = field % val
!      else
!         !option_path = '/material_phase[' // int2str( icomp - 1 ) // &
!         !     ']/scalar_field::ComponentMassFractionPhase' // &
!         !     int2str( iphase )
!
!         Conditional_Composition_MassFraction: if ( have_option( trim( option_path ) // &
!              '/prognostic/initial_condition::WholeMesh/constant' ) ) then
!
!            call get_option( trim( option_path ) // &
!                 '/prognostic/initial_condition::WholeMesh/constant', initial_constant )
!            field_prot = initial_constant
!
!         elseif( have_option( trim( option_path ) // &
!              '/prognostic/initial_condition::WholeMesh/python') ) then
!
!            call get_option( trim( option_path ) // &
!                 '/prognostic/initial_condition::WholeMesh/python', func )
!
!            call allocate( dummy, field % mesh, 'dummy' )
!            call get_option( '/timestepping/current_time', current_time )
!            call set_from_python_function( dummy, trim( func ), positions, current_time )
!            field_prot = dummy % val
!            call deallocate( dummy )
!         elseif( have_option( trim( option_path ) // '/prognostic/initial_condition/from_file')) then
!            field_prot = field % val
!         else
!            ewrite(-1,*) 'No initial condition for field::', trim( field_name )
!            FLAbort( ' Check initial conditions ' )
!
!         end if Conditional_Composition_MassFraction
!
!      end if Conditional_InitialisedFromFLML
!      if ( present(wic_bc)) then
!
!
!         Conditional_Composition_BC: if ( have_option( trim( option_path ) // &
!              '/prognostic/boundary_conditions[0]/type::dirichlet' )) then
!
!            BC_Type = 1
!            nobcs = get_boundary_condition_count( field )
!
!            Loop_Over_BC: do k = 1, nobcs
!               field_prot_bc => extract_surface_field( field, k, 'value' )
!               shape_option = option_shape( trim( option_path ) // &
!                    '/prognostic/boundary_conditions[' // &
!                    int2str( k - 1 ) // ']/surface_ids' )
!               allocate( sufid_bc( 1 : shape_option( 1 ) ) )
!
!               call get_option( trim( option_path ) // &
!                    '/prognostic/boundary_conditions[' // &
!                    int2str( k - 1 ) // ']/surface_ids', sufid_bc )
!
!               allocate( face_nodes( face_loc( field, 1 ) ) )
!               sele = 1
!               do j = 1, stotel
!                  if( any ( sufid_bc == pmesh % faces % boundary_ids( j ) ) ) then
!                     wic_bc( j + ( iphase - 1 ) * stotel ) = bc_type
!                     face_nodes = ele_nodes( field_prot_bc, sele )
!                     do kk = 1, snloc
!                        suf_bc( ( icomp - ( nphase + 1 ) ) * nphase * stotel * snloc + &
!                             ( iphase - 1 ) * stotel * snloc + ( j - 1 ) * snloc + kk ) = &
!                             field_prot_bc % val( face_nodes( kk ) )
!                     end do
!                     sele = sele + 1
!                  end if
!               end do
!
!               deallocate( face_nodes )
!               deallocate( sufid_bc )
!
!            end do Loop_Over_BC ! End of BC loop
!
!         end if Conditional_Composition_BC
!
!      end if
!
!
!      if ( have_source )  then
!         do j = 1, node_count( field_source )
!            field_prot_source( ( iphase - 1 ) * node_count( field_source ) + j ) = &
!                 field_source % val( j )
!         end do
!      end if
!
!      if ( have_absorption ) then
!         do j = 1, node_count( field_absorption )
!            field_prot_absorption( ( iphase - 1 ) * node_count( field_absorption ) + j ) = &
!                 field_absorption % val( j )
!         end do
!      end if
!
!      return
!    end subroutine Get_CompositionFields_Outof_State



    subroutine Get_VectorFields_Outof_State( state, initialised, iphase, field, &
         wic_bc, wic_momu_bc, suf_u_bc, suf_v_bc, suf_w_bc, &
         suf_momu_bc, suf_momv_bc, suf_momw_bc, &
         field_prot_source, field_prot_absorption )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      logical, intent( in ) :: initialised
      integer, intent( in ) :: iphase
      type( vector_field ), pointer :: field, field_prot_bc
      !     real, dimension( : ), intent( inout ) :: field_u_prot, field_v_prot, field_w_prot, &
      !          field_nu_prot, field_nv_prot, field_nw_prot
      real, dimension( :, :, : ), intent( inout ), optional :: field_prot_source
      real, dimension( : , :, : ), intent( inout ), optional :: field_prot_absorption
      integer, dimension( : ), intent( inout ), optional :: wic_bc, wic_momu_bc
      real, dimension( : ), intent( inout ), optional :: suf_u_bc, suf_v_bc, suf_w_bc
      real, dimension( : ), intent( inout ), optional :: suf_momu_bc, suf_momv_bc, suf_momw_bc

      ! Local variables
      type( mesh_type ), pointer :: pmesh, cmesh
      type(vector_field), pointer :: positions, field_absorption
      integer, dimension(:), allocatable :: sufid_bc, face_nodes
      character( len = option_path_len ) :: option_path, option_path2, field_name, bct
      integer :: ndim, stotel, snloc, snloc2, nonods, nobcs, bc_type, j, k, kk, l, &
           shape_option( 2 ), count, u_nonods, idim, stat
      logical :: have_absorption

      pmesh => extract_mesh(state, 'PressureMesh' )
      cmesh => extract_mesh(state, 'CoordinateMesh' )
      positions => extract_vector_field( state( 1 ), 'Coordinate' )

      ndim = field % dim
      stotel = surface_element_count( cmesh )
      snloc2 = face_loc( field, 1)
      snloc = snloc2
      nonods = node_count( field )
      field_name = trim( field % name )
      u_nonods = nonods


      have_absorption = .false.
      Conditional_AbsorptionField: if( present( field_prot_absorption ) ) then
         field_absorption => extract_vector_field( state( iphase ), trim(field_name) // 'Absorption', stat )
         option_path = '/material_phase[' // int2str( iphase - 1 ) // ']/vector_field::' // trim( field_name ) // &
              '/prognostic/vector_field::Absorption/diagnostic/algorithm::vector_python_diagnostic'
         have_absorption =  have_option( trim(option_path) )
         if ( have_absorption ) then
            do idim = 1, ndim
               field_prot_absorption( :, idim + (iphase-1)*ndim, idim + (iphase-1)*ndim ) =  &
                    field_absorption % val( idim, : )
            end do
         else
            do idim = 1, ndim
               field_prot_absorption( :, idim + (iphase-1)*ndim, idim + (iphase-1)*ndim ) = 0.0
            end do
         end if
      end if Conditional_AbsorptionField

      if (present(wic_bc)) then

         option_path = '/material_phase[' // int2str( iphase - 1 )// ']/vector_field::' // trim( field_name )
         option_path2 = trim( option_path ) // '/prognostic/boundary_conditions['

         nobcs = get_boundary_condition_count( field )
         Loop_BC: do k = 1, nobcs

            field_prot_bc => extract_surface_field( field, k, 'value' )

            option_path = trim( option_path2 ) // int2str( k - 1 ) // ']/surface_ids'
            shape_option = option_shape( trim( option_path ) )
            allocate( SufID_BC( 1 : shape_option( 1 ) ) )
            call get_option( trim( option_path ), SufID_BC )
            allocate( face_nodes( face_loc( field, 1) ) )

            option_path = trim( option_path2 ) // int2str( k - 1 ) // ']/'

            Conditional_Field_BC: if( have_option( trim( option_path ) // 'type::dirichlet' ) ) then

               BC_Type = 1

               face_nodes = (/ ( l, l = 1, snloc2 ) /)

               do j = 1, stotel
                  if( any ( sufid_bc == field % mesh % faces % boundary_ids( j ) ) ) then
                     wic_bc( j  + ( iphase - 1 ) * stotel ) = BC_Type
                     count = 1
                     do kk = 1, snloc
                        suf_u_bc( ( iphase - 1 ) * stotel * snloc + ( j - 1 ) * snloc + kk ) = &
                             field_prot_bc % val( 1, face_nodes( count ) )
                        if( ndim > 1 ) suf_v_bc( ( iphase - 1 ) * stotel * snloc + &
                             ( j - 1 ) * snloc + kk ) = field_prot_bc % val( 2, face_nodes( count ) )
                        if( ndim > 2 ) suf_w_bc( ( iphase - 1 ) * stotel * snloc + &
                             ( j - 1 ) * snloc + kk ) = field_prot_bc % val( 3, face_nodes( count ) )
                        count = count + 1
                        if ( mod( kk, snloc2 ) == 0. ) count = 1
                     end do
                     face_nodes = face_nodes + snloc2
                  end if
               end do

            else if( have_option( trim( option_path ) // 'type::momentum' ) ) then

               call get_option( trim( option_path ) // 'type::momentum/boundary/', bct )
               if ( trim( bct ) == "incoming") then
                  BC_Type = 1
               else if ( trim( bct ) == "open") then
                  BC_Type = 5
               else
                  FLAbort( 'Wrong momentum boundary condition in diamond.' )
               end if

               face_nodes = (/ ( l, l = 1, snloc2 ) /)
               do j = 1, stotel
                  if( any ( sufid_bc == field % mesh % faces % boundary_ids( j ) ) ) then
                     wic_momu_bc( j  + ( iphase - 1 ) * stotel ) = BC_Type
                     count = 1
                     do kk = 1, snloc
                        suf_momu_bc( ( iphase - 1 ) * stotel * snloc + ( j - 1 ) * snloc + kk ) = &
                             field_prot_bc % val( 1, face_nodes( count ) )
                        if( ndim > 1 ) suf_momv_bc( ( iphase - 1 ) * stotel * snloc + &
                             ( j - 1 ) * snloc + kk ) = field_prot_bc % val( 2, face_nodes( count ) )
                        if( ndim > 2 ) suf_momw_bc( ( iphase - 1 ) * stotel * snloc + &
                             ( j - 1 ) * snloc + kk ) = field_prot_bc % val( 3, face_nodes( count ) )
                        count = count + 1
                        if ( mod( kk, snloc2 ) == 0. ) count = 1
                     end do
                     face_nodes = face_nodes + snloc2
                  end if
               end do

            end if Conditional_Field_BC

            deallocate( face_nodes, SufID_BC )

         end do Loop_BC

      end if

      return
    end subroutine Get_VectorFields_Outof_State



    subroutine Extract_TensorFields_Outof_State( state, istate_field, &
         field, field_path, &
         field_prot_tensor, &
         GlobalNodeNumber )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      integer, intent( in ) :: istate_field
      type( tensor_field ), pointer :: field
      character( len = option_path_len ), intent( in ) :: field_path
      real, dimension( :, :, : ), intent( inout ) :: field_prot_tensor
      integer, dimension( : ), intent( in ), optional :: GlobalNodeNumber
!!$ Local variables
      type( scalar_field ), pointer :: pressure
      type( vector_field ), pointer :: positions
      type( tensor_field ), pointer :: tensorfield
      integer :: nstate, ndim, totele, istate, ele, idim, jdim, iloc, inod
      character( len = option_path_len ) :: field_name, option_path, option_permeability, option_viscosity
      logical :: compute_viscosity, compute_permeability, is_isotropic, is_diagonal, is_anisotropic
      real, dimension( :, : ), allocatable :: constant
      integer, dimension( : ), pointer :: element_nodes

      call get_option( '/geometry/dimension', ndim )
      nstate = option_count( '/material_phase' )
      positions => extract_vector_field( state, 'Coordinate' )
      pressure => extract_scalar_field( state, 'Pressure' )

      totele = ele_count( positions ) ; field_name = trim( field % name )
      compute_viscosity = .false. ; compute_permeability = .false.

!!$ Defining logicals for the field associated with the tensor. For now, it is set up for either
!!$ Permeability( totele, ndim, ndim ) and Viscosity( mat_nonods, ndim, ndim, nphase ). We may need to
!!$ extend the subroutine to also with an arbitrary tensor shape.
      Conditional_WhichField: if( trim( field_path ) == '/material_phase[' // int2str( istate_field - 1 ) // &
           ']/vector_field::Velocity/prognostic/tensor_field::' // trim( field_name ) )then
         compute_viscosity = .true.
      elseif( trim( field_path ) == '/porous_media/tensor_field::Permeability' ) then
         compute_permeability = .true.
      else
         FLAbort( 'Tensor Field was not defined in the schema - sort this out' )
      end if Conditional_WhichField


!!$ Now looping over state
      LoopOverState: do istate = 1, nstate

         if( istate /= istate_field ) cycle LoopOverState

!!$ Defining flags:
         if( compute_permeability ) then
            option_permeability = '/porous_media/tensor_field::Permeability/prescribed/value[0]'
            if( have_option( trim( option_permeability ) ) ) then
               option_path = trim( option_permeability )
            else
               FLAbort( 'Option path for the Permeability tensor is incomplete.' )
            end if
         end if

         if( compute_viscosity ) then
            option_viscosity = '/material_phase[' // int2str( istate - 1 ) // ']/vector_field::' // &
                 'Velocity/prognostic/tensor_field::' // trim( field_name ) // '/prescribed/value[0]'
            if( have_option( trim( option_viscosity ) ) ) then
               option_path = trim( option_viscosity )
            else
               FLAbort( 'Option path for the Viscosity tensor is incomplete.' )
            end if
         end if
!!$
!!$         Conditional_InitialisationFromFLML:if ( initialised ) then
!!$               tensorfield => extract_tensor_field( state( istate ), trim( field_name ) )
!!$            do ele = 1, element_count( tensorfield )
!!$               element_nodes => ele_nodes( tensorfield, ele )
!!$               do idim = 1, ndim
!!$                  do jdim = 1, ndim
!!$                     if( compute_permeability )then
!!$                        field_prot_tensor( ele, idim, jdim ) = field % val( idim, jdim, element_nodes( 1 ) )
!!$                     elseif( compute_viscosity )then
!!$                         do iloc = 1, ele_loc( pressure, 1 )
!!$                              inod = GlobalNodeNumber( ( ele - 1 ) * node_count( pressure ) + iloc )
!!$                              field_prot_tensor( inod, idim, jdim ) = &
!!$                                   tensorfield % val( idim, jdim, element_nodes( 1 ) )
!!$                         end do
!!$                     else
!!$                        FLAbort( 'Incorrect path for the tensor field' )
!!$                     end if
!!$                  end do
!!$               end do
!!$            end do
!!$         end if Conditional_InitialisationFromFLML
!!$
!!$

         is_isotropic = have_option( trim( option_path ) // '/isotropic' )
         is_diagonal = have_option( trim( option_path ) // '/diagonal' )
         is_anisotropic = have_option( trim( option_path ) // '/anisotropic_symmetric' ) .or. &
              have_option( trim( option_path ) // '/anisotropic_asymmetric' )
!!$ Endof Defining flags

!!$ Isotropic Tensor:
         Conditional_Tensor: if ( is_isotropic ) then
            option_path = trim( option_path ) // '/isotropic'

            if( have_option( trim( option_path ) // '/constant')) then
               allocate( constant( 1, 1 ) ) ; constant = 0.
               call get_option( trim( option_path ) // '/constant', constant( 1, 1 ) )
               do idim = 1, ndim
                  field_prot_tensor( : , idim, idim ) = constant( 1, 1 )
               end do
               deallocate( constant )

            elseif( have_option( trim( option_path ) // '/python' ) ) then
               tensorfield => extract_tensor_field( state( istate ), trim( field_name ) )
               do ele = 1, element_count( tensorfield )
                  element_nodes => ele_nodes( tensorfield, ele )
                  do idim = 1, ndim
                     if( compute_permeability ) then
                        field_prot_tensor( ele, idim, idim ) = &
                             tensorfield % val( idim, idim, element_nodes( 1 ) )
                     elseif( compute_viscosity ) then
                        do iloc = 1, ele_loc( pressure, 1 )
                           inod = GlobalNodeNumber( ( ele - 1 ) * node_count( pressure ) + iloc )
                           field_prot_tensor( inod , idim, idim ) =   &
                                tensorfield % val( idim, idim, element_nodes( 1 ) )
                        end do
                     else
                        FLAbort( 'Option path for the tensor field is incomplete.' )
                     end if
                  end do
               end do

            else
               FLExit( 'Incorrect initial condition for field' )
            end if
!!$ Endof Isotropic Tensor

!!$ Diagonal Tensor:
         elseif( is_diagonal )then ! If_Conditional_Tensor
            option_path = trim( option_path ) // '/diagonal'

            if( have_option( trim( option_path ) // '/constant' ) ) then
               allocate(constant( ndim, 1 ) ) ; constant = 0.
               call get_option( trim( option_path ) // '/constant', constant )
               do idim = 1, ndim
                  field_prot_tensor( : , idim, idim ) = constant( idim, 1 )
               end do
               deallocate( constant )

            elseif( have_option( trim( option_path ) // '/python' ) )then
               tensorfield => extract_tensor_field( state( istate ), trim( field_name ) )
               ! element_nodes => ele_nodes( tensorfield, ele )
               do idim = 1, ndim
                  do ele = 1, element_count( tensorfield )
                     element_nodes => ele_nodes( tensorfield, ele )
                     if( compute_permeability ) then
                        field_prot_tensor( ele, idim, idim ) = tensorfield % val( idim, idim,  ele )
                     elseif( compute_viscosity ) then
                        do iloc = 1, ele_loc( pressure, 1 )
                           inod = GlobalNodeNumber( ( ele - 1 ) * node_count( pressure ) + iloc )
                           field_prot_tensor( inod , idim, idim ) =   &
                                tensorfield % val( idim, idim, element_nodes( 1 ) )
                        end do
                     else
                        FLAbort( 'Option path for the tensor field is incomplete.' )
                     end if
                  end do
               end do

            else
               FLExit( 'Incorrect initial condition for field' )
            end if
!!$ Endof Diagonal Tensor

!!$ Anisotropic Tensor:
         elseif( is_anisotropic )then  ! If_Conditional_Tensor
            if( have_option( trim( option_path ) // '/anisotropic_symmetric' ) ) then
               option_path = trim( option_path ) // '/anisotropic_symmetric'
            elseif( have_option( trim( option_path ) // '/anisotropic_asymmetric' ) ) then
               option_path = trim( option_path ) // '/anisotropic_asymmetric'
            else
               FLAbort( 'Option path for the tensor field is incomplete.' )
            end if

            if( have_option( trim( option_path ) // '/constant' ) ) then
               allocate( constant( ndim, ndim ) ) ; constant = 0.
               call get_option( trim( option_path ) // '/constant', constant )
               do idim = 1, ndim
                  do jdim = 1, ndim
                     field_prot_tensor( : , idim, jdim ) = constant( idim, jdim )
                  end do
               end do

            elseif( have_option( trim( option_path ) // '/python' ) ) then
               tensorfield => extract_tensor_field( state( istate ), trim( field_name ) )
               do ele = 1, element_count( tensorfield )
                  element_nodes => ele_nodes( tensorfield, ele )
                  do idim = 1, ndim
                     do jdim = 1, ndim
                        if( compute_permeability ) then
                           field_prot_tensor( ele, idim, jdim ) = &
                                tensorfield % val( idim, jdim, element_nodes( 1 ) )
                        elseif( compute_viscosity ) then
                           do iloc = 1, ele_loc( pressure, 1 )
                              inod = GlobalNodeNumber( ( ele - 1 ) * node_count( pressure ) + iloc )
                              field_prot_tensor( inod, idim, jdim ) = &
                                   tensorfield % val( idim, jdim, element_nodes( 1 ) )
                           end do
                        else
                           FLAbort( 'Option path for the tensor field is incomplete.' )
                        end if
                     end do
                  end do
               end do

            else
               FLExit( 'Incorrect initial condition for field' )

            end if
!!$ Endof Anisotropic Tensor

         else
            FLExit( 'Incorrect initial condition for field' )

         end if Conditional_Tensor

      end do LoopOverState

      return
    end subroutine Extract_TensorFields_Outof_State


!!$
!!$ Module Get_Ndgln Interfaces

!    subroutine Get_Scalar_Ndgln( ndgln, field, cv_nloc )
!      implicit none
!      type( scalar_field ), intent( in ) :: field
!      integer, intent( in ), optional :: cv_nloc
!      integer, dimension( : ), intent( inout ) :: ndgln
!      ! Local variables
!      integer, dimension( : ), pointer :: nloc
!      integer :: ele, iloc
!
!      do ele = 1, ele_count( field )
!         nloc => ele_nodes( field, ele )
!         do iloc = 1, ele_loc( field, ele )
!            ndgln( ( ele - 1 ) * ele_loc( field, ele ) + iloc ) =  nloc( iloc )
!$            ewrite(3,*)'ele, iloc, ndgln:', ele, iloc, &
!$                 ndgln( ( ele - 1 ) * ele_loc( field, ele ) + iloc )
!         end do
!      end do
!
!      return
!    end subroutine Get_Scalar_Ndgln

!    subroutine Get_Vector_Ndgln( ndgln, field, cv_nloc )
!      implicit none
!      type( vector_field ), intent( in ) :: field
!      integer, intent( in ), optional :: cv_nloc
!      integer, dimension( : ), intent( inout ) :: ndgln
!      ! Local variables
!      integer, dimension( : ), pointer :: nloc
!      integer :: ele, iloc, count, cv_nloc2, ndim
!
!      call get_option( '/geometry/dimension', ndim )
!      cv_nloc2 = 1
!
!      count = 0
!      do ele = 1, ele_count( field )
!         nloc => ele_nodes( field, ele )
!         do iloc = 1, ele_loc( field, 1 ) * cv_nloc2
!            ndgln( ( ele - 1 ) * ele_loc( field, 1 ) * cv_nloc2  + iloc ) =  nloc( iloc )
!         end do
!      end do
!
!      return
!    end subroutine Get_Vector_Ndgln


!    subroutine Get_Mesh_Ndgln( ndgln, mesh, cv_nloc )
!      implicit none
!      type( mesh_type ), intent( in ) :: mesh
!      integer, intent( in ), optional :: cv_nloc
!      integer, dimension( : ), intent( inout ) :: ndgln
!      ! Local variables
!      integer, dimension( : ), pointer :: nloc
!      integer :: ele, iloc
!
!      do ele = 1, ele_count( mesh )
!         nloc => ele_nodes( mesh, ele )
!         do iloc = 1, ele_loc( mesh, ele )
!            ndgln( ( ele - 1 ) * ele_loc( mesh, ele ) + iloc ) =  nloc( iloc )
!$            ewrite(3,*)'ele, iloc, ndgln:', ele, iloc, &
!$                 ndgln( ( ele - 1 ) * ele_loc( mesh, ele ) + iloc )
!         end do
!      end do
!
!      return
!    end subroutine Get_Mesh_Ndgln

!!$
!!$ Module Get_SNdgln Interfaces

    subroutine Get_Scalar_SNdgln( sndgln, field, cv_nloc  )
      implicit none
      type( scalar_field ), intent( in ) :: field
      integer, dimension( : ), intent( inout ) :: sndgln
      integer, intent( in ), optional :: cv_nloc
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

    subroutine Get_Vector_SNdgln( sndgln, field, cv_nloc  )
      implicit none
      type( vector_field ), intent( in ) :: field
      integer, dimension( : ), intent( inout ) :: sndgln
      integer, intent( in ), optional :: cv_nloc
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




    subroutine print_from_state( state, field_prot )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      real, dimension( : ), intent( in ) :: field_prot
      !
      type( scalar_field ), pointer :: field
      integer :: nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, knod, istate
      logical :: initialised
      integer, dimension( : ), pointer :: cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln
      integer, dimension( : ), allocatable ::     cv_sndgln, p_sndgln, u_sndgln, Temperature_BC_Spatial
      real, dimension( : ), allocatable :: Temperature, Temperature_BC

      call Get_Primary_Scalars( state, &
           nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods )

!!$ Calculating Global Node Numbers
      allocate( cv_sndgln( stotel * cv_snloc ), p_sndgln( stotel * p_snloc ), &
           u_sndgln( stotel * u_snloc ) )

      !    x_ndgln_p1 = 0 ; x_ndgln = 0 ; cv_ndgln = 0 ; p_ndgln = 0 ; mat_ndgln = 0 ; u_ndgln = 0 ; xu_ndgln = 0 ; &
      cv_sndgln = 0 ; p_sndgln = 0 ; u_sndgln = 0

      allocate( temperature( nphase * cv_nonods ), temperature_bc_spatial( nphase * stotel ), &
           temperature_bc( stotel * cv_snloc * nphase ) )

      call Compute_Node_Global_Numbers( state, &
           totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
           cv_snloc, p_snloc, u_snloc, &
           cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln, &
           cv_sndgln, p_sndgln, u_sndgln )

      initialised = .true.
      do istate = 1, nstate
         Conditional_Temperature: if( have_option( '/material_phase[' // int2str( istate - 1 ) // &
              ']/scalar_field::Temperature' ) ) then
            field => extract_scalar_field( state( istate ), 'Temperature' )
            knod = ( istate - 1 ) * node_count( field )
            call Get_ScalarFields_Outof_State( state, initialised, istate, field, &
                 Temperature( knod + 1 : knod + node_count( field ) ), &
                 Temperature_BC_Spatial, Temperature_BC ) !, &
!!$                 field_prot_source = Temperature_Source( knod + 1 : knod + node_count( field ) ) )
         end if Conditional_Temperature
      end do

      ewrite(3,*)'::temperature::', norm2( temperature )
      do istate = 1, cv_nonods
         ewrite(3,*) istate, field_prot( istate ), temperature( istate ), field%val(istate)
      end do

      deallocate( cv_sndgln, p_sndgln, u_sndgln, Temperature_BC_Spatial, Temperature, Temperature_BC )

      return
    end subroutine print_from_state


    subroutine update_boundary_conditions( state, stotel, cv_snloc, nphase, &
         &                                 suf_t_bc, suf_t_bc_rob1, suf_t_bc_rob2 )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      integer, intent( in ) :: stotel, cv_snloc, nphase
      real, dimension( 1, nphase, stotel * cv_snloc ), intent( inout ) :: suf_t_bc, suf_t_bc_rob1, suf_t_bc_rob2
      !
      character( len = option_path_len ) :: option_path, option_path2, field_name, name
      integer :: shape_option(2), iphase, nobcs, kk, k, j , sele, stat
      integer, dimension( : ), allocatable :: SufID_BC
      integer, dimension( : ), pointer:: surface_element_list, face_nodes

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

    subroutine pack_multistate( npres, state, packed_state,&
         multiphase_state, multicomponent_state, pmulti_state )

      type(state_type), dimension(:), intent(inout):: state
      type(state_type), dimension(:), intent(inout), pointer :: &
           multiphase_state, multicomponent_state
      type(state_type) :: packed_state
      integer, intent(in) :: npres

      type(state_type), dimension(:,:), pointer, optional :: pmulti_state

      type(state_type), dimension(:,:), pointer :: multi_state

      integer :: i,nphase,ncomp,ndim,stat,iphase,icomp,idim,ele,ipres,n_in_pres

      type(scalar_field), pointer :: pressure, sfield
      type(vector_field), pointer :: velocity, position, vfield
      type(tensor_field), pointer :: tfield, p2

      type(vector_field) :: porosity, vec_field
      type(vector_field) :: p_position, u_position, m_position
      type(tensor_field) :: permeability, ten_field
      type(mesh_type), pointer :: ovmesh, element_mesh
      type(element_type) :: element_shape

      integer, dimension( : ), pointer :: element_nodes

      logical :: has_density, has_phase_volume_fraction

      ncomp=option_count('/material_phase/is_multiphase_component')
      nphase=size(state)-ncomp

      position=>extract_vector_field(state(1),"Coordinate")
      ndim=mesh_dim(position)

      call insert(packed_state,position%mesh,"CoordinateMesh")


#ifdef USING_FEMDEM
      if ( have_option( '/blasting' ) ) then
         sfield => extract_scalar_field( state(1), "SolidConcentration" )
         call insert( packed_state, sfield, "SolidConcentration" )
         call add_new_memory(packed_state,sfield,"OldSolidConcentration")

         vfield => extract_vector_field( state(1), "delta_U" )
         call insert( packed_state, vfield, "delta_U" )

         vfield => extract_vector_field( state(1), "solid_U" )
         call insert( packed_state, vfield, "solid_U" )


         vfield => extract_vector_field( state(1), "f_x" )
         call insert( packed_state, vfield, "f_x" )

         tfield => extract_tensor_field( state(1), "a_xx" )
         call insert( packed_state, tfield, "a_xx" )

         tfield => extract_tensor_field( state(1), "Viscosity" )
         call insert( packed_state, tfield, "Viscosity" )

	elseif ( have_option( '/femdem_fracture' ) ) then
         sfield => extract_scalar_field( state(1), "SolidConcentration" )
         call insert( packed_state, sfield, "SolidConcentration" )
         call add_new_memory(packed_state,sfield,"OldSolidConcentration")
         
         tfield => extract_tensor_field( state(1), "Viscosity" )
         call insert( packed_state, tfield, "Viscosity" )

         sfield => extract_scalar_field( state(1), "Dummy" )
         call insert( packed_state, sfield, "Dummy" )

         vfield => extract_vector_field( state(1), "Darcy_Velocity" )
         call insert( packed_state,vfield, "Darcy_Velocity" )

         vfield => extract_vector_field( state(1), "delta_U" )
         call insert( packed_state, vfield, "delta_U" )

         vfield => extract_vector_field( state(1), "solid_U" )
         call insert( packed_state, vfield, "solid_U" )

      end if
#endif

      if (has_scalar_field(state(1),"Porosity")) then
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

      !If have capillary pressure, then we store 5 entries in PackedRockFluidProp, otherwise just 3
      if( have_option_for_any_phase( '/multiphase_properties/capillary_pressure', nphase ) ) then
        call allocate(ten_field,element_mesh,"PackedRockFluidProp",dim=[5,nphase])
      else
          call allocate(ten_field,element_mesh,"PackedRockFluidProp",dim=[3,nphase])
      end if
      !Introduce the rock-fluid properties (Immobile fraction, Krmax, relperm exponent -Capillary entry pressure, capillary exponent-)
      call insert(packed_state,ten_field,"PackedRockFluidProp")
      call deallocate(ten_field)


      allocate(multiphase_state(nphase))
      allocate(multicomponent_state(ncomp))
      allocate(multi_state(max(1,ncomp),nphase))

      pressure=>extract_scalar_field(state(1),"Pressure")
      call insert(packed_state,pressure%mesh,"PressureMesh")

!      call add_new_memory(packed_state,pressure,"FEPressure")
!      call add_new_memory(packed_state,pressure,"OldFEPressure")
!      call add_new_memory(packed_state,pressure,"CVPressure")
!      call add_new_memory(packed_state,pressure,"OldCVPressure")

      call insert_sfield( packed_state,"FEPressure",1,npres )

      tfield => extract_tensor_field( packed_state, "PackedFEPressure" )
      tfield%option_path = pressure%option_path

      call insert_sfield( packed_state,"CVPressure",1,npres )




      ! dummy field on the pressure mesh, used for evaluating python eos's.
      ! this could be cleaned up in the future.
      call add_new_memory(packed_state,pressure,"Dummy")

      tfield => extract_tensor_field( state(1), "Dummy", stat )
      if ( stat==0 ) call insert( packed_state, tfield, "Dummy" )

      p2=>extract_tensor_field(packed_state,"PackedFEPressure")
      !call set( p2, pressure )
      do ipres = 1, npres
         p2%val(1,ipres,:)=pressure%val
      end do
      do icomp=1,ncomp
         call insert(multicomponent_state(icomp),p2,"PackedFEPressure")
      end do

      p2=>extract_tensor_field(packed_state,"PackedCVPressure")
      !call set( p2, pressure )
      do ipres = 1, npres
         p2%val(1,ipres,:)=pressure%val
      end do

      call insert_sfield(packed_state,"FEDensity",1,nphase)

      call insert_sfield(packed_state,"Density",1,nphase)
      call insert_sfield(packed_state,"DensityHeatCapacity",1,nphase)


      if (option_count("/material_phase/scalar_field::Temperature")>0) then
         call insert_sfield(packed_state,"Temperature",1,nphase,&
              add_source=.true.,add_absorption=.true.)
         call insert_sfield(packed_state,"FETemperature",1,nphase)
      end if
      call insert_sfield(packed_state,"PhaseVolumeFraction",1,nphase,&
           add_source=.true.)
      call insert_sfield(packed_state,"FEPhaseVolumeFraction",1,nphase)

      !If we have capillary pressure we create a field in packed_state
      if( have_option_for_any_phase( '/multiphase_properties/capillary_pressure', nphase ) ) then
         call allocate(ten_field,pressure%mesh,"PackedCapPressure",dim=[1,nphase])
         call insert(packed_state,ten_field,"PackedCapPressure")
         call deallocate(ten_field)
      end if

      velocity=>extract_vector_field(state(1),"Velocity")
      call insert(packed_state,velocity%mesh,"VelocityMesh")
      if ( .not. has_mesh(state(1),"VelocityMesh_Continuous") ) then
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

      if ( .not. has_mesh(state(1),"PressureMesh_Continuous") ) then
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

      if ( .not. has_mesh(state(1),"PressureMesh_Discontinuous") ) then
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


      call remap_field( position, u_position )
      call remap_field( position, p_position )
      call remap_field( position, m_position )

      call insert(packed_state,position,"Coordinate")
      call insert(packed_state,u_position,"VelocityCoordinate")
      call insert(packed_state,p_position,"PressureCoordinate")
      call insert(packed_state,m_position,"MaterialCoordinate")
      call deallocate(p_position)
      call deallocate(u_position)
      call deallocate(m_position)

      has_density=has_scalar_field(state(1),"Density")
      has_phase_volume_fraction=has_scalar_field(state(1),"PhaseVolumeFraction")

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

      iphase=1
      icomp=1
      do i=1,size(state)
         if(have_option(trim(state(i)%option_path)&
              //'/is_multiphase_component')) then
            velocity=>extract_vector_field(state(i),"Velocity",stat)
            if (stat==0) velocity%wrapped=.true.
            velocity=>extract_vector_field(state(i),"OldVelocity",stat)
            if (stat==0) velocity%wrapped=.true.
            velocity=>extract_vector_field(state(i),"NonlinearVelocity",stat)
            if (stat==0) velocity%wrapped=.true.
            velocity=>extract_vector_field(state(i),"IteratedVelocity",stat)
            if (stat==0) velocity%wrapped=.true.

            call unpack_component_sfield(state(i),packed_state,"FEComponentDensity",icomp)
            !call unpack_component_sfield(state(i),packed_state,"FEOldComponentDensity",icomp)
            call unpack_component_sfield(state(i),packed_state,"FEComponentMassFraction",icomp)
            !call unpack_component_sfield(state(i),packed_state,"FEOldComponentMassFraction",icomp)


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
            if (stat==0) pressure%wrapped=.true.

            if (has_density) then
               call unpack_sfield(state(i),packed_state,"IteratedDensity",1,iphase,&
                    check_paired(extract_scalar_field(state(i),"Density"),&
                    extract_scalar_field(state(i),"IteratedDensity")))
               call unpack_sfield(state(i),packed_state,"OldDensity",1,iphase,&
                    check_paired(extract_scalar_field(state(i),"Density"),&
                    extract_scalar_field(state(i),"OldDensity")))
               call unpack_sfield(state(i),packed_state,"Density",1,iphase)
               call insert(multi_state(1,iphase), extract_scalar_field(state(i),"Density"),"Density")
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
               call insert(multi_state(1,iphase), extract_scalar_field(state(i),"Temperature"),"Temperature")
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
               call insert(multi_state(1,iphase), extract_scalar_field(state(i),"PhaseVolumeFraction"),"PhaseVolumeFraction")
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
            call insert(multi_state(1,iphase), extract_vector_field(state(i),"Velocity"),"Velocity")

            iphase=iphase+1
         end if
      end do

      n_in_pres = nphase / npres
      do ipres=1,npres
         i = (ipres-1)*n_in_pres + 1
         call unpack_sfield(state(i),packed_state,"Pressure",1,ipres)
         call insert(multi_state(1,ipres), extract_scalar_field(state(i),"Pressure"),"FEPressure")
      end do


      if (option_count("/material_phase/scalar_field::Temperature")>0) call allocate_multiphase_scalar_bcs(packed_state,multi_state,"Temperature")
      call allocate_multiphase_scalar_bcs(packed_state,multi_state,"Density")
      call allocate_multiphase_scalar_bcs(packed_state,multi_state,"PhaseVolumeFraction")
      call allocate_multiphase_vector_bcs(packed_state,multi_state,"Velocity")
      call allocate_multiphase_scalar_bcs(packed_state,multi_state,"FEPressure")

      if (ncomp>0) then
         call unpack_multicomponent(packed_state,multicomponent_state)
         call allocate_multicomponent_scalar_bcs(multicomponent_state,multi_state,"ComponentMassFraction")
         call allocate_multicomponent_scalar_bcs(multicomponent_state,multi_state,"ComponentDensity")
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


      if (present(pmulti_state)) then
         pmulti_state=>multi_state
      else
         call deallocate(multi_state)
         deallocate(multi_state)
      end if

      call allocate(porosity,npres,element_mesh,"Porosity")
      do ipres = 1, npres
         call set(porosity,ipres,1.0)
      end do
      call insert(packed_state,porosity,"Porosity")
      call deallocate(porosity)
      if (has_scalar_field(state(1),"Porosity")) then
         sfield => extract_scalar_field(state(1),"Porosity")
         porosity%val(1,:) = sfield%val
      end if

      ! Hack to define a lateral from diamond
      if ( npres >  1 ) then
         vfield => extract_vector_field(packed_state,"Porosity")
         sfield => extract_scalar_field(state(1),"Pipe1")
         vfield%val(2,:) = sfield%val
      end if

      if (has_scalar_field(state(1),"Permeability")) then
         call allocate(permeability,element_mesh,"Permeability",&
              dim=[mesh_dim(position),mesh_dim(position)])
         call zero(permeability)
         sfield=>extract_scalar_field(state(1),"Permeability")
         do idim=1,mesh_dim(position)
            call set(permeability,idim,idim,sfield)
         end do
         call insert(packed_state,permeability,"Permeability")
         call deallocate(permeability)
      else if (has_vector_field(state(1),"Permeability")) then
         call allocate(permeability,element_mesh,"Permeability",&
              dim=[mesh_dim(position),mesh_dim(position)])
         call zero(permeability)
         vfield=>extract_vector_field(state(1),"Permeability")
         call set(permeability,vfield)
         call insert(packed_state,permeability,"Permeability")
         call deallocate(permeability)
      else if (has_tensor_field(state(1),"Permeability")) then
         call allocate(permeability,element_mesh,"Permeability",&
              dim=[mesh_dim(position),mesh_dim(position)])
         call zero(permeability)
         tfield => extract_tensor_field( state(1), trim( permeability % name ) )
         if( size(tfield%val,3) == 1 ) then!constant field
            do ele = 1, element_count( tfield )
               Permeability%val( :, :, ele ) = tfield % val( :, :, 1 )
            end do
         else!python
            do ele = 1, element_count( tfield )
               element_nodes => ele_nodes( tfield, ele )
               Permeability%val( :, :, ele ) = &
                    tfield % val( :, :, element_nodes( 1 ) )
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

!      !Add memory for the DarcyVelocity
!      call allocate(ten_field,velocity%mesh,"PackedDarcyVelocity")
!      call zero(ten_field)
!      call insert(packed_state,ten_field,"PackedDarcyVelocity")
!      call deallocate(ten_field)
!      !Insert into state and merge memories
!      do iphase = 1, size(state)
!          call allocate(vec_field, ndim, velocity%mesh, "DarcyVelocity")
!          call zero(vec_field)
!          call insert(state(iphase),vec_field,"DarcyVelocity")
!          call deallocate(vec_field)
!      end do
!      do iphase = 1, size(state)
!          call unpack_vfield(state(iphase),packed_state,"DarcyVelocity",iphase,&
!               check_vpaired(extract_vector_field(state(iphase),"DarcyVelocity"),&
!               extract_vector_field(state(iphase),"DarcyVelocity")))
!      end do

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
              vfield%val=>tfield%val(icomp:icomp,:,:)
              vfield%wrapped=.true.
              call insert(mcstate(icomp),vfield,vfield%name)
              call deallocate(vfield)
           END do

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
        integer :: index, iphase, si

        type(tensor_field), pointer :: tfield
        type(tensor_field) :: mp_tfield


        do index=1,size(mstate%tensor_fields)
           tfield=>extract_tensor_field(mstate,index)
           si=len(trim(tfield%name))
!!-PY changed it
           if(tfield%name(si-7:si)=="Pressure")then
!           if(tfield%name(si-3:si)=="Pressure")then
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
        type(tensor_field), pointer :: mfield
        logical lfree


        if (present(free)) then
           lfree=free
        else
           lfree=.true.
        end if

        if ( trim(name)=="Pressure" ) then
           mfield=>extract_tensor_field(mstate,"PackedFE"//name)
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
           if(lfree .and. associated(nfield%val)) then
#ifdef HAVE_MEMORY_STATS
              call register_deallocation("scalar_field", "real", &
                   size(nfield%val), nfield%name)
#endif
              deallocate(nfield%val)
           end if
           nfield%val=>mfield%val(icomp,iphase,:)
           nfield%val_stride=ncomp*nphase
           nfield%wrapped=.true.
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



!    subroutine add_dependant_fields_to_tensor_from_state(infield,state,scalar_field_names,&
!         vector_field_names,tensor_field_names)
!
!      !  Convenience subroutine to add a bunch of fields as dependants of infield in one call.
!
!      type(tensor_field) :: infield
!      type(state_type) :: state
!      character (len=*) , dimension(:), optional :: scalar_field_names, vector_field_names, &
!           tensor_field_names
!
!      integer :: i
!      type(scalar_field), pointer :: sfield
!      type(vector_field), pointer :: vfield
!      type(tensor_field), pointer :: tfield
!
!      if (present(scalar_field_names)) then
!         do i=1,size(scalar_field_names)
!            sfield=>extract_scalar_field(state,trim(scalar_field_names(i)))
!            call add_dependant_field(infield,sfield)
!         end do
!      end if
!
!      if (present(vector_field_names)) then
!         do i=1,size(vector_field_names)
!            vfield=>extract_vector_field(state,trim(vector_field_names(i)))
!            call add_dependant_field(infield,vfield)
!         end do
!      end if
!
!      if (present(tensor_field_names)) then
!         do i=1,size(tensor_field_names)
!            tfield=>extract_tensor_field(state,tensor_field_names(i))
!            call add_dependant_field(infield,tfield)
!         end do
!      end if
!
!    end subroutine add_dependant_fields_to_tensor_from_state


    subroutine Adaptive_NonLinear(packed_state, reference_field, its,&
         Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,order)
      !This subroutine either store variables before the nonlinear timeloop starts, or checks
      !how the nonlinear iterations are going and depending on that increase the timestep
      !or decreases the timestep and repeats that timestep
      Implicit none
      type(state_type), intent(inout) :: packed_state!, backup_state
      real, dimension(:,:,:), allocatable, intent(inout) :: reference_field
      logical, intent(inout) :: Repeat_time_step, ExitNonLinearLoop
      logical, intent(in) :: nonLinearAdaptTs
      integer, intent(in) :: its, order
      !Local variables
      real :: dt
      logical, save :: show_FPI_conv
      logical, save :: Time_step_decreased_with_dumping = .false.
      real, save :: OldDt
      real, save :: Accumulated_sol
      real, parameter :: check_sat_threshold = 1d-6
      real, dimension(:,:,:), pointer :: pressure
      real, dimension(:,:), pointer :: phasevolumefraction
      real, dimension(:,:,:), pointer :: velocity

      !        real, dimension(:,:,:), pointer :: tVar, tVar_it
      !        real, dimension(:,:), pointer :: vVar, vVar_it
      !        real, dimension(:,:), pointer :: sVar, sVar_it

      !Variables for automatic non-linear iterations
      real :: tolerance_between_non_linear, initial_dt, min_ts, max_ts, increase_ts_switch, decrease_ts_switch,&
        Inifinite_norm_tol
      !Variables for adaptive time stepping based on non-linear iterations
      real :: increaseFactor, decreaseFactor, ts_ref_val, acctim, inf_norm_val
      integer :: variable_selection, NonLinearIteration
!ewrite(0,*) "entering"
      !First of all, check if the user wants to do something
      call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration', tolerance_between_non_linear, default = -1. )
      if (tolerance_between_non_linear<0) return
      !Tolerance for the infinite norm
      call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Inifinite_norm_tol',&
             Inifinite_norm_tol, default = 0.03 )
      !retirve number of Fixed Point Iterations
      call get_option( '/timestepping/nonlinear_iterations', NonLinearIteration, default = 3 )
      !Get data from diamond. Despite this is slow, as it is done in the outest loop, it should not affect the performance.
      !Variable to check how good nonlinear iterations are going 1 (Pressure), 2 (Velocity), 3 (Saturation)
      call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear', &
           variable_selection, default = 3)
      call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/increase_factor', &
           increaseFactor, default = 1.05 )
      call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/decrease_factor', &
           decreaseFactor, default = 1.2 )
      call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/max_timestep', &
           max_ts, default = huge(min_ts) )
      call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/min_timestep', &
           min_ts, default = -1. )
      show_FPI_conv = have_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Show_Convergence')
      !Switches are relative to the input value unless otherwise stated
      if (have_option('/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/increase_ts_switch')) then
          call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/increase_ts_switch', &
          increase_ts_switch, default = 1d-3 )
      else
          increase_ts_switch = tolerance_between_non_linear / 10.
      end if

      if (have_option('/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/decrease_ts_switch')) then
          call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/adaptive_timestep_nonlinear/decrease_ts_switch', &
          decrease_ts_switch, default = 1d-1 )
      else
          decrease_ts_switch = min(tolerance_between_non_linear * 10.,1.0)
      end if

      !Get time step
      call get_option( '/timestepping/timestep', initial_dt )
      dt = initial_dt
      !By default the minimum time-steps is ten orders smaller than the initial timestep
      if(min_ts<0) min_ts = initial_dt * 1d-10

      select case (order)
      case (1)!Store or get from backup
         !If we do not have adaptive time stepping then there is nothing to backup
         if (.not.nonLinearAdaptTs) return
         !we either store the data or we recover it if repeting a timestep
         !Procedure to repeat time-steps
         !If  Repeat_time_step then we recover values, else we store them
         if (its == 1) call copy_packed_new_to_iterated(packed_state, Repeat_time_step)
      case (2)!Calculate and store reference_field
            !Store variable to check afterwards
            call get_var_from_packed_state(packed_state, velocity = velocity, pressure = pressure,&
                 phasevolumefraction = phasevolumefraction)

            select case (variable_selection)
            case (1)
                if (allocated(reference_field)) then
                    if (size(reference_field,3) /= size(pressure,3) ) then
                        deallocate(reference_field)
                        allocate (reference_field(1,1,size(pressure,3) ))
                    end if
                else
                    allocate (reference_field(1,1,size(pressure,3) ))
                end if
               reference_field(1,1,:) = pressure(1,1,:)
            case (2)
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
            case default

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
            end select

      case default!Check how is the process going on and decide
         !If Automatic_NonLinerIterations then we compare the variation of the a property from one time step to the next one
         ExitNonLinearLoop = .false.
         Repeat_time_step = .false.

         if (its == 1) Accumulated_sol = dumping_in_sat

         if (its > 1 ) then

            call get_var_from_packed_state(packed_state, velocity = velocity, pressure = pressure,&
                 phasevolumefraction = phasevolumefraction)

            select case (variable_selection)
            case (1)
               ts_ref_val = maxval(abs(reference_field(1,1,:)-pressure(1,1,:)))
            case (2)
               ts_ref_val = maxval(abs(reference_field-velocity))
            case default
               !Calculate infinite norm
               inf_norm_val = maxval(abs(reference_field(1,:,:)-phasevolumefraction))/dumping_in_sat

               !Calculate value of the functional
               ts_ref_val = get_Convergence_Functional(phasevolumefraction, reference_field(1,:,:), dumping_in_sat)
            end select

            !If it is parallel then we want to be consistent between cpus
            if (IsParallel()) call allmax(ts_ref_val)
            !We cannot go to the next time step until we have performed a full time step
            Accumulated_sol = Accumulated_sol + dumping_in_sat
            if (IsParallel()) call allmax(Accumulated_sol)

            !TEMPORARY, re-use of global variable dumping_in_sat to send
            !information about convergence to the trust_region_method
            dumping_in_sat = ts_ref_val

            ewrite(1,*) "FPI convergence: ", "L2^2 norm=>",ts_ref_val,"; L_inf norm =>", inf_norm_val, "; Total iterations:", its

            !If only non-linear iterations
            if (.not.nonLinearAdaptTs) then
               !Automatic non-linear iteration checking
                ExitNonLinearLoop = ((ts_ref_val < tolerance_between_non_linear .and. inf_norm_val < Inifinite_norm_tol)&
                     .or. its >= NonLinearIteration)
                   if (ExitNonLinearLoop .and. show_FPI_conv) then
                        !Tell the user the number of FPI and final convergence to help improving the parameters
                         print *, "FPI convergence: ", "L2^2 norm=>",ts_ref_val,"; L_inf norm =>", inf_norm_val, "; Total iterations:", its
                   end if
               return
            end if


            !If we have a dumping parameter we only reduce the time-step if we reach
            !the maximum number of non-linear iterations
            if (.not. have_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Dumping_factor')) then
                !Tell the user if we have not converged
                if (its == NonLinearIteration) then
                    ewrite(1,*) "Fixed point method failed to converge in ",NonLinearIteration,"iterations, final convergence is", ts_ref_val
                end if
                !Increase Ts section
                if ((ts_ref_val < increase_ts_switch .and.dt*increaseFactor<max_ts).and..not.Repeat_time_step) then
                   call get_option( '/timestepping/timestep', dt )
                   dt = dt * increaseFactor
                   call set_option( '/timestepping/timestep', dt )
                   ewrite(1,*) "Time step increased to:", dt
                   ExitNonLinearLoop = .true.
                   return
                else !Maybe it is not enough to increase the time step, but we could go to the next time step
                    ExitNonLinearLoop = (ts_ref_val < tolerance_between_non_linear)
                end if

                !Decrease Ts section only if we have done at least the 90% of the  nonLinearIterations
                if ((ts_ref_val > decrease_ts_switch.or.repeat_time_step) &
                     .and.its>=int(0.90*NonLinearIteration)) then

                   if ( dt / decreaseFactor < min_ts) then
                      !Do not decrease
                      Repeat_time_step = .false.
                      ExitNonLinearLoop = .true.
                      deallocate(reference_field)
                      return
                   end if

                   !Decrease time step, reset the time and repeat!
                   call get_option( '/timestepping/timestep', dt )
                   call get_option( '/timestepping/current_time', acctim )
                   acctim = acctim - dt
                   call set_option( '/timestepping/current_time', acctim )
                   dt = dt / decreaseFactor
                   call set_option( '/timestepping/timestep', dt )
                   ewrite(1,*) "Time step decreased to:", dt
                   Repeat_time_step = .true.
                   ExitNonLinearLoop = .true.
                end if
            else!Adaptive Ts for Dumping based on the number of FPI
                if (ts_ref_val < tolerance_between_non_linear) then
                    if (its < int(0.25 * NonLinearIteration) .and..not.Repeat_time_step) then
                       !Increase time step
                       call get_option( '/timestepping/timestep', dt )
                       dt = dt * increaseFactor
                       call set_option( '/timestepping/timestep', dt )
                       ewrite(1,*) "Time step increased to:", dt
                       ExitNonLinearLoop = .true.
                       return
                    end if
                else if (its >= NonLinearIteration) then
                !If it has not converged when reaching the maximum number of non-linear iterations,
                !reduce ts and repeat
                    !Decrease time step for next time step
                       if ( dt / decreaseFactor < min_ts) then
                          !Do not decrease
                          Repeat_time_step = .false.
                          ExitNonLinearLoop = .true.
                          deallocate(reference_field)
                          return
                       end if

                       !Decrease time step, reset the time and repeat!
                       call get_option( '/timestepping/timestep', dt )
                       call get_option( '/timestepping/current_time', acctim )
                       acctim = acctim - dt
                       call set_option( '/timestepping/current_time', acctim )
                       dt = dt / decreaseFactor
                       call set_option( '/timestepping/timestep', dt )
                       ewrite(1,*) "Time step decreased to:", dt
                       Repeat_time_step = .true.
                       ExitNonLinearLoop = .true.
                end if


                !For adaptive time stepping we need to put this again
                if (ExitNonLinearLoop .and. show_FPI_conv) then
                    !Tell the user the number of FPI and final convergence to help improving the parameters
                    ewrite(0,*) "FPI convergence:", ts_ref_val, "Total iterations:", its
                end if

            end if
         end if

      end select

    end subroutine Adaptive_NonLinear

    real function get_Convergence_Functional(phasevolumefraction, reference_sat, dumping)
       !We create the potential to optimize F = sum (f**2), so the solution is when this potential
       !reach a minimum. Typically the value to consider convergence is the sqrt(epsilon of the machine), i.e. 10^-8
       !f = (NewSat-OldSat)/Number of nodes; this is the typical approach for algebraic non linear systems
       !
       !The convergence is independent of the dumping parameter
       !and measures how the previous iteration (i.e. using the previous dumping parameter) performed
        implicit none
        real, dimension(:,:), intent(in) :: phasevolumefraction, reference_sat
        real, intent(in) :: dumping
        !Local variables
        integer :: cv_inod, modified_vals, iphase
        real :: aux
        real, parameter :: tol = 1d-5

        modified_vals = 0
        get_Convergence_Functional = 0.0

        !(L2)**2 norm of all the elements
        do iphase = 1, size(phasevolumefraction,1)
            get_Convergence_Functional = max(sum((abs(reference_sat(iphase,:)-phasevolumefraction(iphase,:))&
            /size(phasevolumefraction,2))**2.0), get_Convergence_Functional)
        end do

!        !(L2)**2 norm of all the elements whose value has changed (has problems to converge at
!        !the beginning since only a shock front is happening, however it is better than simple l2 norm)
!        do cv_inod = 1, size(phasevolumefraction,2)
!            aux = maxval(abs(reference_sat(:,cv_inod)-phasevolumefraction(:,cv_inod)))
!            if (aux > tol) then
!                get_Convergence_Functional = get_Convergence_Functional + aux**2.0
!                modified_vals = modified_vals + 1
!            end if
!        end do
!        get_Convergence_Functional = (get_Convergence_Functional / dble(modified_vals)**2.0)

        !Rescale using the dumping in saturation to get a more efficient number to compare with
        !if the dumping_in_sat was 10-2 then ts_ref_val will always be small
        !To make consistent the dumping parameter with the Potential, we have to raise it to 2.0
        get_Convergence_Functional = get_Convergence_Functional / dumping**2.0


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

      tfield=>extract_tensor_field(packed_state,"OldFEPressure")
      ntfield=>extract_tensor_field(packed_state,"FEPressure")
      tfield%val=ntfield%val

      tfield=>extract_tensor_field(packed_state,"OldCVPressure")
      ntfield=>extract_tensor_field(packed_state,"CVPressure")
      tfield%val=ntfield%val

    end subroutine copy_packed_new_to_iterated

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
         EndPointRelperm, RelpermExponent, Cap_entry_pressure, Cap_exponent)
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
           Immobile_fraction, EndPointRelperm, RelpermExponent, Cap_entry_pressure, Cap_exponent
      real, optional, dimension(:,:,:), pointer ::Pressure,FEPressure, OldFEPressure, CVPressure,OldCVPressure
      !Local variables
      type(scalar_field), pointer :: sfield
      type(vector_field), pointer :: vfield
      type(tensor_field), pointer :: tfield

      !Scalar stored
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


    end subroutine get_var_from_packed_state


    !Subroutine to print Arrays by (columns,rows)
    !Matrix = 2D Array
    subroutine printMatrix(Matrix)
      implicit none

      Integer :: length,i,j, k
      character (len=1000000) :: cadena
      character (len=100) :: aux
      real, intent(in), dimension(:,:):: Matrix
      !Local
      real, dimension(size(matrix,2),size(matrix,1)) :: auxMatrix

      auxMatrix = transpose(Matrix)

      length = size(auxMatrix,2);
      do i = 1,size(auxMatrix,1)
         print *,""
         cadena = ""
         do j = 1 , length
            write(aux,*), auxMatrix(i,j)
            k = index(trim(aux),"E",.true.)
            if (k/=0) then
               aux = aux(1:k-6)//trim(aux(k:))
            end if

            cadena = trim(cadena)//' '//trim(aux)
         end do
         print '(A $)', trim(cadena)
      end do

      print *,"";
    end subroutine PrintMatrix

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


    subroutine Clean_Storage(storage_state, StorageIndexes)
      !This subroutine removes all the storage controlled by StorageIndexes
      Implicit none
      type(state_type), intent(inout) :: storage_state
      integer, dimension(:), intent(inout) :: StorageIndexes
      !Local variables
      character (len = 100) :: StorName
      integer :: maxpos, ignored_var_pos, i

         !This loop is the most robust, so by default we still use this one
         do while (maxval(abs(StorageIndexes)) > 0)
            maxpos = maxloc(abs(StorageIndexes), dim =1)
            StorName = trim(storage_state%scalar_names(abs(StorageIndexes(maxpos))))!This lines is
            call remove_scalar_field(storage_state, trim(StorName))           !failing for Xie when using adaptive meshing
            StorageIndexes(maxpos) = 0
         end do
         !Just in case
         StorageIndexes = 0
    end subroutine Clean_Storage


    subroutine CheckElementAngles(packed_state, totele, x_ndgln, X_nloc, MaxAngle, MinAngle, Quality_list, degree)
        !This function checks the angles of an input element. If one angle is above
        !the Maxangle or below the MinAngle it will be true in the list
        Implicit none
        !Global variables
        type(state_type), intent(inout) :: packed_state
        type(bad_elements), dimension(:), intent(inout) :: Quality_list
        real, intent (in) :: MaxAngle, MinAngle
        integer, dimension(:), intent(in) :: x_ndgln
        integer, intent(in) :: x_nloc, totele, degree
        !Local variables
        integer :: ELE, i
        logical :: Bad_founded
        real :: MxAngl, MnAngl
        !Definition of Pi
        real, dimension(:,:), pointer:: X_ALL
        real, parameter :: pi = acos(0.d0) * 2d0

        !Prepare data
        do i = 1, size(Quality_list)
            Quality_list(i)%ele = -1!Initialize with negative values
            allocate(Quality_list(i)%nodes(3))!Allocate, we always consider triangles
            Quality_list(i)%nodes = -1!Initialize with negative values
            allocate(Quality_list(i)%weights(x_nloc-1))!Allocate
            Quality_list(i)%weights = 0.!Initialize with zeros
        end do
        call get_var_from_packed_state(packed_state, PressureCoordinate = X_ALL)

        !Convert input angles to radians
        MxAngl = pi/180. * MaxAngle
        MnAngl = pi/180. * MinAngle
        i = 1
        if (size(X_ALL,1)==2) then!2D triangles
            do ELE = 1, totele
                !bad_node enters as the first entry
                if (degree == 1) then
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 1, 2, 3, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                else if (degree == 2) then!quadratic
                    !Here we consider three subtriangles, the normal one, plus two formed by the midpoints plus an edge
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 4, 5, 6, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 1, 2, 4, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 2, 3, 5, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                end if
            end do
        else if(size(X_ALL,1)==3) then!3D tetrahedra
        !adjust to match the 2D case once that one works properly
            do ELE = 1, totele
                if (degree == 1) then
                    !We check the 4 triangles that form a tet
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 1, 2, 3, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 1, 2, 4, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 1, 4, 3, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 4, 2, 3, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                else if (degree == 2) then!quadratic
                    !We check the 4 triangles that form a tet, three times
                    !First face
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 2, 3, 5, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 1, 2, 4, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 4, 5, 6, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    !Second face
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 2, 3, 8, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 1, 2, 7, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 7, 8, 10, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    !Third face
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 1, 4, 7, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 4, 6, 9, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 7, 9, 10, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    !Forth face
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 10, 8, 9, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 8, 3, 5, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if
                    if (i > size(Quality_list)) exit!We cannot add more elements
                    Bad_founded = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 9, 5, 6, MxAngl, MnAngl, Quality_list(i))
                    if (Bad_founded) then
                        if (size(Quality_list)>= i) Quality_list(i)%ele = ele
                        i = i + 1
                    end if

                end if
            end do
        end if


        if (size(Quality_list) < i) then
            ewrite(1,*) 'WARNING: The number of bad elements is bigger than expected not all of them will be compensated. Reduce them or increase the size of Quality_list'
        end if

        contains

            logical function Check_element(X_ALL, x_ndgln, ele_Pos, Pos1, Pos2, Pos3, MaxAngle, MinAngle, Quality_list)
                !Checks if an angle is below a threshold and stores the information to
                !modify the matrix later
                implicit none
                real, dimension(:,:), intent(in) :: X_ALL
                integer, intent(in) :: ele_Pos, Pos1, Pos2, Pos3
                real, intent(in) :: MaxAngle, MinAngle
                integer, dimension(:), intent(in) :: x_ndgln
                type(bad_elements), intent(inout) :: Quality_list
                !Local variables
                real, dimension(size(X_ALL,1)) :: X1, X2, X3
                real, dimension(3) :: alpha, lenght
                !For 3D to project values
                !Definition of Pi
                real, parameter :: pi = acos(0.0) * 2.0

                Check_element = .false.
                !Define the vertexes
                X1 = X_ALL(:, x_ndgln(ele_Pos+Pos1))
                X2 = X_ALL(:, x_ndgln(ele_Pos+Pos2))
                X3 = X_ALL(:, x_ndgln(ele_Pos+Pos3))


                !Calculate the lenght of the edges
                lenght(1) = sqrt(dot_product(X1(:)-X3(:), X1(:)-X3(:)))
                lenght(2) = sqrt(dot_product(X1(:)-X2(:), X1(:)-X2(:)))
                lenght(3) = sqrt(dot_product(X2(:)-X3(:), X2(:)-X3(:)))

                !Alphas
                alpha(2) = acos((lenght(3)**2+lenght(2)**2-lenght(1)**2)/(2. *lenght(3)*lenght(2)))
                alpha(1) = acos((lenght(1)**2+lenght(2)**2-lenght(3)**2)/(2. *lenght(1)*lenght(2)))
                alpha(3) = pi - alpha(1)-alpha(2)

                !Check angles and if necessary calculate weights and bad node, I don't know if this work for 3D...
                !for the time being just 2D
                if (alpha(1)>=MaxAngle) then
                    !We calculate weights considering a right triangle formed by the bad node, its projection and the other
                    !corner. So the cosine of the angles time their sides sum the side in which we are projecting the bad node
                    Quality_list%weights(1) = abs(cos(alpha(2)) * lenght(2) / lenght(3))
                    Quality_list%weights(2) = abs(cos(alpha(3)) * lenght(1) / lenght(3))
                    !Store nodes, the first one is the bad node
                    Quality_list%nodes(1) = Pos1
                    Quality_list%nodes(2) = Pos2
                    Quality_list%nodes(3) = Pos3

                    !Store angle so later the over-relaxation can depend on this
                    Quality_list%angle = alpha(1) * 180 / pi

                    Check_element = .true.
                else if (alpha(2)>= MaxAngle) then
                    Quality_list%weights(1) = abs(cos(alpha(1)) * lenght(2) / lenght(1))
                    Quality_list%weights(2) = abs(cos(alpha(3)) * lenght(3) / lenght(1))
                    !Store nodes, the first one is the bad node
                    Quality_list%nodes(1) = Pos2
                    Quality_list%nodes(2) = Pos1
                    Quality_list%nodes(3) = Pos3
                    !Store angle so later the over-relaxation can depend on this
                    Quality_list%angle = alpha(2) * 180 / pi
                    Check_element = .true.
                else if (alpha(3) >= MaxAngle) then
                    Quality_list%weights(1) = abs(cos(alpha(1)) * lenght(1) / lenght(2))
                    Quality_list%weights(2) = abs(cos(alpha(2)) * lenght(3) / lenght(2))
                    !Store nodes, the first one is the bad node
                    Quality_list%nodes(1) = Pos3
                    Quality_list%nodes(2) = Pos1
                    Quality_list%nodes(3) = Pos2
                    !Store angle so later the over-relaxation can depend on this
                    Quality_list%angle = alpha(3) * 180 / pi
                    Check_element = .true.
                end if
                !Make sure it is between bounds
                Quality_list%weights = min(Quality_list%weights,1.0)

                !If we have not added elements already in that element, check for small angles
                if (.not.Check_element) then
                    if (alpha(1) < MinAngle ) then
                        Quality_list%weights(1) = abs(cos(alpha(1)) * lenght(1) / lenght(2))
                        Quality_list%weights(2) = abs(cos(alpha(2)) * lenght(3) / lenght(2))
                        !Store nodes, the first one is the bad node
                        Quality_list%nodes(1) = Pos1
                        Quality_list%nodes(2) = Pos2
                        Quality_list%nodes(3) = Pos3
                        !We fake this parameter since it is preapred for obtuse angles
                        !we consider a medium angle
                        Quality_list%angle = alpha(1) * 180 / pi
                        Check_element = .true.
                    else if (alpha(2) < MinAngle) then
                        Quality_list%weights(1) = abs(cos(alpha(2)) * lenght(2) / lenght(3))
                        Quality_list%weights(2) = abs(cos(alpha(3)) * lenght(1) / lenght(3))
                        !Store nodes, the first one is the bad node
                        Quality_list%nodes(1) = Pos2
                        Quality_list%nodes(2) = Pos3
                        Quality_list%nodes(3) = Pos1
                        !We fake this parameter since it is preapred for obtuse angles
                        !we consider a medium angle
                        Quality_list%angle = alpha(2) * 180 / pi
                        Check_element = .true.
                    else if (alpha(3) < MinAngle) then
                        Quality_list%weights(1) = abs(cos(alpha(1)) * lenght(2) / lenght(1))
                        Quality_list%weights(2) = abs(cos(alpha(3)) * lenght(3) / lenght(1))
                        !Store nodes, the first one is the bad node
                        Quality_list%nodes(1) = Pos3
                        Quality_list%nodes(2) = Pos1
                        Quality_list%nodes(3) = Pos2
                        !We fake this parameter since it is preapred for obtuse angles
                        !we consider a medium angle
                        Quality_list%angle = alpha(3) * 180 / pi
                        Check_element = .true.
                    end if
                end if

            end function Check_element

    end subroutine CheckElementAngles

    subroutine calculate_outflux(packed_state, ndotqnew, sele, surface_ids, totoutflux, ele , x_ndgln,&
         cv_ndgln, cv_nloc, SCVFEN, gi, cv_nonods, totele, nphase, detwei, IDs_ndgln, cv_snloc, cv_siloc ,SUF_T_BC_ALL)

!    subroutine calculate_outflux(packed_state, ndotqnew, sele, surface_ids, totoutflux, &
!         cv_sndgln, cv_snloc, SCVFEN, gi, cv_nonods, totele, nphase, detwei, IDs_ndgln, SUF_T_BC_ALL)

       implicit none

! Subroutine to calculate the integrated flux across a boundary with the specified surface_ids.

! Input/Output variables

       type(state_type), intent(in) :: packed_state
       real, dimension(:), intent(in) :: ndotqnew
       integer, intent(in) :: sele
       integer, dimension(1), intent(in) :: surface_ids
       real, dimension(:), intent(inout) :: totoutflux
       integer, intent(in) :: ele
       integer, dimension(:), intent( in ) ::  x_ndgln
       integer, dimension(:), intent( in ) :: IDs_ndgln
       integer, dimension(:), intent( in ) ::  cv_ndgln
       !integer, dimension(:), intent( in ) ::  cv_sndgln
       integer, intent(in) :: cv_nloc
       integer, intent(in) :: cv_snloc
       integer, intent(in) :: cv_siloc
       real, dimension( : , : ), intent(in), pointer :: SCVFEN
       integer, intent(in) :: gi
       integer, intent(in) :: cv_nonods
       integer, intent(in) :: totele
       integer, intent(in) :: nphase
       real, pointer, dimension( : ), intent(in) :: detwei
       real, dimension( :,:, : ), intent( in ) :: SUF_T_BC_ALL

! Local variables

       type(scalar_field), pointer :: sfield
       type(vector_field), pointer :: vfield
       !real, dimension(:), pointer :: CVPressure
       real, dimension(:,:,:), pointer :: CVPressure
       real, dimension(:), pointer :: Por
       logical :: test
       type(tensor_field), pointer :: tfield
       type(tensor_field), pointer :: t2field, t3field
       integer  :: x_knod
       integer  :: cv_knod
       integer  :: cv_sknod
       integer :: x_ele
       integer  :: cv_kloc
       integer  :: cv_skloc
       real, dimension( : , : ), allocatable :: phaseV
       real, dimension( : , : ), allocatable :: phaseVG
       real, dimension( : , : ), allocatable :: Dens
       real, dimension( : , : ), allocatable :: DensVG
       !real, dimension( : ), allocatable :: PorG
       real :: PorG
       integer :: surf
       integer :: i

       ! SHOULD RETHINK THESE ALLOCATIONS - only need to allocate # gauss points worth of memory
       allocate(phaseV(nphase,cv_nonods), phaseVG(nphase,cv_nonods))
       allocate(Dens(nphase,cv_nonods), DensVG(nphase,cv_nonods))
       !allocate(PorG(totele))

! Extract the pressure

      !sfield => extract_scalar_field( packed_state, "CVPressure" )
      !CVPressure =>  sfield%val(:)

      ! Modified 25/06/15 to account for the changes to the treatment of Pressure in the code (Using DPPs multiple pressures formulation)

      t3field => extract_tensor_field( packed_state, "PackedCVPressure" )
      CVPressure => t3field%val(:,:,:)

! Extract the phase volume fraction

      tfield => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
      phaseV = tfield%val(1,:,:)

 ! Extract the density

      t2field => extract_tensor_field( packed_state, "PackedDensity" )
      Dens =  t2field%val(1,:,:)

! Extract the porosity (still should confirm whether or not totoutflux needs to be divided by the porosity - matter of which velocity to use : q or v)

      vfield => extract_vector_field( packed_state, "Porosity" )
      Por =>  vfield%val(1,:)

! Having extracted the saturation field (phase volume fraction) at control volume nodes, need to calculate it at quadrature points gi.
! (Note saturation is defined on a control volume basis and so the field is stored at control volume nodes). Therefore need cv_ndgln below

      phaseVG = 0.0
      do cv_kloc=1,cv_nloc
      cv_knod=cv_ndgln((ele-1)*cv_nloc+cv_kloc)
      phaseVG(:,gi) = phaseVG(:,gi) + phaseV(:,cv_knod)*SCVFEN(cv_kloc,gi)
      end do

!      phaseVG = 0.0
!      do cv_skloc=1,cv_snloc
!      cv_sknod=cv_sndgln((sele-1)*cv_snloc+cv_skloc)
!      phaseVG(:,gi) = phaseVG(:,gi) + phaseV(:,cv_sknod)*SCVFEN(cv_skloc,gi)
!      end do

! Having extracted the density at control volume nodes, need to calculate it at quadrature points gi.
! Density is a function of pressure and therefore lives on the pressure mesh. It is therefore stored at
! control volume nodes (easiest to see in a diagram of the elements). Hence, the global variable we need in the calculation below is cv_ndgln

      DensVG = 0.0
      do cv_kloc=1,cv_nloc
      cv_knod=cv_ndgln((ele-1)*cv_nloc+cv_kloc)
      DensVG(:,gi) = DensVG(:,gi) + Dens(:,cv_knod)*SCVFEN(cv_kloc,gi)
      end do

!      DensVG = 0.0
!      do cv_skloc=1,cv_snloc
!      cv_sknod=cv_sndgln((sele-1)*cv_snloc+cv_skloc)
!      DensVG(:,gi) = DensVG(:,gi) + Dens(:,cv_sknod)*SCVFEN(cv_skloc,gi)
!      end do


! Porosity is on the P0DG mesh. The value is constant across an element. Can calculate it at any point in an element and assign this to be the value
! at the Gauss point. (Strictly this last bit is unnecessary - may tidy up in the future). x_ndgln((ele-1)*cvnloc +1) allows us to extract the
! x-coordinate of the '1st physical' triangle node for each element. We then calculate the porosity at this value of x and assign it to PorG(ele).
! [Check that porosity is constant across an element and not across a control volume - else this needs modification].

      ! x_ele calculates a coordinate on each element
      !x_ele=x_ndgln((ele-1)*cv_nloc+1)
      !PorG(ele) = Por(x_ele)

     PorG = Por(IDs_ndgln(ele))

!      PorG = Por(IDs_ndgln(sele))

! This function will return true for surfaces we should be integrating over (this entire subroutine is called in a loop over ele,(sele),gi in cv-adv-dif)
! Need the condition that sele > 0 Check why it can be zero/negative.

      if(sele > 0) then

          test = integrate_over_surface_element(t3field, sele, surface_ids)

! Need to integrate the fluxes over the boundary in question (i.e. those that test true). Totoutflux initialised to zero out of this subroutine. Ndotqnew caclulated in cv-adv-diff
! Need to add up these flow velocities multiplied by the saturation phaseVG to get the correct velocity and by the Gauss weights to get an integral. Density needed to get a mass flux
! Divide by porosity to get a physical flux (check this)

          if(test) then

              ! In the case of an inflow boundary, need to use the boundary value of saturation (not the value inside the domain)
              ! Need to pass down an array with the saturation boundary conditions to deal with these cases
              ! i.e need to pass down SUF_T_BC_ALL(1, nphase, surface_element)

              do i = 1, size(ndotqnew)
                      surf = (sele - 1 ) * cv_snloc + cv_siloc
                      if(ndotqnew(i) < 0 ) then
                          ! Inlet boundary - so use boundary phase volume fraction
                          totoutflux(i) = totoutflux(i) + ndotqnew(i)*SUF_T_BC_ALL(1, i, surf)*detwei(gi)!*DensVG(i,gi)!/PorG
                          !totoutflux(i) = totoutflux(i) + ndotqnew(i)*detwei(gi)*DensVG(i,gi)/PorG

                      else
                          ! Outlet boundary - so use internal (to the domain) phase volume fraction
                          totoutflux(i) = totoutflux(i) + ndotqnew(i)*phaseVG(i,gi)*detwei(gi)!*DensVG(i,gi)!/PorG
                      endif
              enddo

          endif

      endif

      deallocate(phaseV)
      deallocate(phaseVG)
      deallocate(dens)
      deallocate(densVG)
      !deallocate(PorG)

      return

    end subroutine calculate_outflux

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

    subroutine get_regionIDs2nodes(state, packed_state, CV_NDGLN, IDs_ndgln, IDs2CV_ndgln, fake_IDs_ndgln)
        !This subroutine creates a conversor so material variables
        !can be stored based on region ids, but accessed in a normal way (node access)
        !Also re-adapts the material properties to work in this new way.
        !IDs2CV_ndgln gives you the value regarding one node, if it happens to have many it will be the last
        !value introduced
        implicit none
        type(state_type), dimension(:), intent(inout) :: state
        type(state_type), intent( inout ) :: packed_state
        integer, dimension(:), intent(in) :: CV_NDGLN
        integer, dimension(:), allocatable, intent(inout) :: IDs_ndgln
        integer, dimension(:), allocatable, intent(inout) :: IDs2CV_ndgln
        logical, optional, intent(in) :: fake_IDs_ndgln
        !Local variables

        type (tensor_field), pointer :: t_field
        type(mesh_type), pointer :: fl_mesh
        integer :: i, j, k, number_of_ids, nphase,mtemp
        integer, dimension(:), allocatable :: region_ids
        logical :: stored, all_fields_costant
        integer, dimension(1) :: aux
        character(len=200):: path, root_path
        !Use P0DG mesh
        fl_mesh => extract_mesh( state(1), "P0DG" )

        if (.not.associated(fl_mesh%region_ids)) FLAbort("P0DG mesh not defined or if using adaptivity preserve_mesh_regions is off")

        t_field => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
        nphase = size(t_field%val,2)
        !Re-allocate if necessary
        if (allocated(IDs_ndgln)) then
            if (size(IDs_ndgln)/=element_count(fl_mesh)) then
                deallocate(IDs_ndgln)
                allocate(IDs_ndgln(element_count(fl_mesh)))
            end if
            if (size(IDs2CV_ndgln)/=size(t_field%val,3)) then
               deallocate(IDs2CV_ndgln)
               allocate(IDs2CV_ndgln(size(t_field%val,3)))
            end if
        else
            allocate(IDs2CV_ndgln(size(t_field%val,3)))
            allocate(IDs_ndgln(element_count(fl_mesh)))
        end if

        !Check if all the fields are constant whitin region ids, otherwise we cannot compact the data.
        !Despite this seems restrictive, it is unlikely to use non constant values for region ids
        !since it goes against the surface based modeling idea
        all_fields_costant = .true.
        !Check capillary
        if (have_option_for_any_phase('/multiphase_properties/capillary_pressure/', nphase)) then
            root_path = '/multiphase_properties/capillary_pressure/'&
            //'type_Brooks_Corey/scalar_field::C/prescribed/value'
            k = 0
            do i = 0, nphase-1
                k = max(k,option_count('/material_phase['// int2str( i ) //']'//trim(root_path)))
            end do
            do i = 0, k-1
                path = trim(root_path)//'['//int2str(i)//']/python'
                if (have_option_for_any_phase(trim(path), nphase))&
                    all_fields_costant = .false.
            end do
            root_path = '/multiphase_properties/capillary_pressure/'&
            //'type_Brooks_Corey/scalar_field::a/prescribed/value'
            k = 0
            do i = 0, nphase-1
                k = max(k,option_count('/material_phase['// int2str( i ) //']'//trim(root_path)))
            end do
            do i = 0, k-1
                path = trim(root_path)//'['//int2str(i)//']/python'
                if (have_option_for_any_phase(trim(path), nphase))&
                    all_fields_costant = .false.
            end do
        end if
        !Check relative permeability
        if (have_option_for_any_phase('/multiphase_properties/Relperm_Corey/', nphase)) then
            root_path = '/multiphase_properties/Relperm_Corey/relperm_max/'&
            //'scalar_field::relperm_max/prescribed/value'
            k = 0
            do i = 0, nphase-1
                k = max(k,option_count('/material_phase['// int2str( i ) //']'//trim(root_path)))
            end do
            do i = 0, k-1
                path = trim(root_path)//'['//int2str(i)//']/python'
                if (have_option_for_any_phase(trim(path), nphase))&
                    all_fields_costant = .false.
            end do

            root_path = '/multiphase_properties/Relperm_Corey/relperm_exponent/'&
            //'scalar_field::relperm_exponent/prescribed/value'
            k = 0
            do i = 0, nphase-1
                k = max(k,option_count('/material_phase['// int2str( i ) //']'//trim(root_path)))
            end do
            do i = 0, k-1
                path = trim(root_path)//'['//int2str(i)//']/python'
                if (have_option_for_any_phase(trim(path), nphase))&
                    all_fields_costant = .false.
            end do
        end if
        !Check permeability
        if (have_option('porous_media/scalar_field::Permeability')) then
            root_path = 'porous_media/scalar_field::Permeability/prescribed/value'
            k = option_count(trim(root_path))
            do i = 0, k-1
                path = trim(root_path)//'['//int2str(i)//']/python'
                if (have_option(trim(path)))&
                    all_fields_costant = .false.
            end do
        end if
        if (have_option('porous_media/tensor_field::Permeability')) then
            root_path = 'porous_media/tensor_field::Permeability/prescribed/value'
            k = option_count(trim(root_path))
            do i = 0, k-1
                path = trim(root_path)//'['//int2str(i)//']/isotropic/python'
                if (have_option(trim(path)))&
                    all_fields_costant = .false.
                path = trim(root_path)//'['//int2str(i)//']/diagonal/python'
                if (have_option(trim(path)))&
                    all_fields_costant = .false.
                path = trim(root_path)//'['//int2str(i)//']/anisotropic_symmetric/python'
                if (have_option(trim(path)))&
                    all_fields_costant = .false.
                path = trim(root_path)//'['//int2str(i)//']/anisotropic_asymmetric/python'
                if (have_option(trim(path)))&
                    all_fields_costant = .false.
            end do
        end if
        if(have_option('porous_media/vector_field::Permeability')) then
            all_fields_costant = .false.
        end if
        if(have_option('porous_media/Permeability_from_femdem')) then
            all_fields_costant = .false.
        end if
        !Check porosity
        if (have_option('porous_media/scalar_field::Porosity')) then
            root_path = 'porous_media/scalar_field::Porosity/prescribed/value'
            k = option_count(trim(root_path))
            do i = 0, k-1
                path = trim(root_path)//'['//int2str(i)//']/python'
                if (have_option(trim(path)))&
                    all_fields_costant = .false.
            end do
        end if

        !If fake_IDs_ndgln, then we are not using compacted data and
        !IDs_ndgln and IDs2CV_ndgln will point to the same position
        if (present_and_true(fake_IDs_ndgln) .or. .not. all_fields_costant) then
            do i = 1, size(IDs_ndgln)
                IDs_ndgln(i) = i
            end do

            do i = 1, size(IDs_ndgln)
                do j = 1, size(CV_ndgln)/size(IDs_ndgln)
                    k = CV_ndgln((i-1)* size(CV_ndgln)/size(IDs_ndgln) + j)
                    IDs2CV_ndgln(k) = IDs_ndgln(i)
                end do
            end do

            return
        end if

        allocate(region_ids(size(fl_mesh%region_ids)))
        region_ids = -1
        !Store all the regions ids that appear
        do i = 1, size(fl_mesh%region_ids)
            !Check if already store
            stored = .false.; j = 1
            do while (region_ids(j)>0)
                if (fl_mesh%region_ids(i) == region_ids(j)) stored = .true.
                j = j + 1
            end do
            if (.not.stored) region_ids(j) = fl_mesh%region_ids(i)
        end do
        !Return the number of region ids to properly allocate the fields with this
        number_of_ids = j - 1

        !Store the position where fl_mesh%region_ids(i) appears
        do i = 1, size(fl_mesh%region_ids)
            !The number should appear only once
            aux = MAXLOC(region_ids, MASK = region_ids == fl_mesh%region_ids(i))
            IDs_ndgln(i) = aux(1)
        end do

        !Create IDs2CV_ndgln
        mtemp = size(CV_NDGLN)/size(IDs_ndgln)
        DO i = 1, size(IDs_ndgln)
            !DO j = 1, size(CV_NDGLN)/size(IDs_ndgln)
            DO j = 1, mtemp
                !k = CV_NDGLN(( i - 1 ) * size(CV_NDGLN)/size(IDs_ndgln) + j )
                k = CV_NDGLN(( i - 1 )*mtemp +j )
                IDs2CV_ndgln(k) = IDs_ndgln(i)
            end do
        end do


        !###Compact fields###
        !Relative permeability and Immobile fractions (if cappressure, also cap parameters)
        if (has_tensor_field(packed_state,"PackedRockFluidProp")) then
            t_field=>extract_tensor_field(packed_state,"PackedRockFluidProp")
            call convert_tensor_field(t_field, IDs_ndgln )
        end if

        deallocate(region_ids)

    contains

        subroutine convert_scalar_field(s_field, IDs_ndgln )
            !This subroutine converts an scalar field to use region ids
            implicit none
            integer, dimension(:), allocatable, intent(inout) :: IDs_ndgln
            type (scalar_field), intent(inout), pointer :: s_field
            !Local variables
            real, dimension(:), allocatable :: s_field_bak
            integer :: i

            !Create backup
            allocate(s_field_bak(size(s_field%val,1))); s_field_bak = s_field%val
            !re-size the field
            deallocate(s_field%val); allocate(s_field%val(number_of_ids))
            !Re-store the data
            do i = 1, size(IDs_ndgln)
                s_field%val(IDs_ndgln(i)) = s_field_bak(i)
            end do
!To keep the registry of memory correct we have to re-adapt it as well
#ifdef HAVE_MEMORY_STATS
call register_deallocation("scalar_field", "real", &
    size(s_field_bak,1), s_field%name)
#endif
!            s_field%mesh%nodes = size(IDs_ndgln)
!            s_field%mesh%elements = size(IDs_ndgln)
#ifdef HAVE_MEMORY_STATS
call register_allocation("scalar_field", "real", &
size(s_field%val,1), s_field%name)
#endif
            deallocate(s_field_bak)
        end subroutine convert_scalar_field

        subroutine convert_tensor_field(t_field, IDs_ndgln )
            !This subroutine converts an scalar field to use region ids
            implicit none
            integer, dimension(:), allocatable, intent(inout) :: IDs_ndgln
            type (tensor_field), intent(inout), pointer :: t_field
            !Local variables
            real, dimension(:, :, :), allocatable :: t_field_bak
            integer :: i

            !Create backup
            allocate(t_field_bak(size(t_field%val,1), size(t_field%val,2)&
                , size(t_field%val,3))); t_field_bak = t_field%val
            !re-size the field
            deallocate(t_field%val); allocate(t_field%val(size(t_field%val,1)&
                , size(t_field%val,2), number_of_ids))
            !Re-store the data
            do i = 1, size(IDs_ndgln)
                t_field%val(:, :, IDs_ndgln(i)) = t_field_bak(:, :, i)
            end do

!To keep the registry of memory correct we have to re-adapt it as well
#ifdef HAVE_MEMORY_STATS
call register_deallocation("tensor_field", "real", &
    size(t_field_bak,1)*size(t_field_bak,2)*size(t_field_bak,3), t_field%name)
#endif
            t_field%mesh%nodes = size(IDs_ndgln)
            t_field%mesh%elements = size(IDs_ndgln)

#ifdef HAVE_MEMORY_STATS
call register_allocation("tensor_field", "real", &
size(t_field%val,1)*size(t_field%val,2)*size(t_field%val,3), t_field%name)
#endif
            deallocate(t_field_bak)
        end subroutine convert_tensor_field

    end subroutine get_regionIDs2nodes

    subroutine get_DarcyVelocity(totele, cv_nloc, u_nloc, mat_nloc, MAT_NDGLN, U_NDGLN, CV_NDGLN, &
            state, packed_state, Material_Absorption)
        !This subroutine calculates the actual Darcy velocity, unfinished
        implicit none
        integer, intent(in) :: totele, cv_nloc, u_nloc, mat_nloc
        integer, dimension(:) :: MAT_NDGLN, U_NDGLN, CV_NDGLN
        type(state_type) , intent(in):: packed_state
        type( state_type ), dimension( : ), intent( inout ) :: state
        real, dimension(:,:,:), intent(in):: Material_Absorption
        !Local variables
        type(tensor_field), pointer :: DarcyVelocity, Velocity, saturation, oldsaturation, perm
        type(scalar_field), pointer :: porosity
        real, dimension(:,:), allocatable :: matrix
        integer :: cv_iloc, u_iloc, ele, iphase, ndim, nphase, imat, u_inod, cv_loc, idim
        real, allocatable, dimension(:) :: aux
        real :: aux2


        DarcyVelocity => extract_tensor_field(packed_state, "PackedDarcyVelocity")
        Velocity => extract_tensor_field( packed_state, "PackedVelocity" )
        saturation => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
        oldsaturation => extract_tensor_field( packed_state, "PackedOldPhaseVolumeFraction" )
        porosity => extract_scalar_field(state(1), "Porosity")
        perm=>extract_tensor_field(packed_state,"Permeability")


        call zero(DarcyVelocity)

        allocate(matrix(size(Material_Absorption,2), size(Material_Absorption,3)))
        nphase = size(Velocity%val, 2)
        ndim  = size(Velocity%val, 1)
        allocate(aux(ndim))

        do ele = 1, totele
            do u_iloc = 1, u_nloc
                u_inod = U_NDGLN(( ELE - 1 ) * u_nloc +u_iloc )
                do cv_iloc = 1, cv_nloc
                    imat = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC +CV_ILOC )
                    cv_loc = CV_NDGLN(( ELE - 1 ) * cv_nloc +CV_ILOC )
                    !This is not optimal, maybe just perform when CVN(U_ILOC, CV_INOD) =/ 0
!                    matrix = Material_Absorption(imat,:,:)
!                    call invert(matrix)
                    do iphase = 1, nphase
                        !Inverse of sigma avoiding inversion
                        matrix(ndim*(iphase-1)+1:iphase*ndim,ndim*(iphase-1)+1:iphase*ndim) = matmul(perm%val(:,:,ele),&
                            Material_Absorption(imat,ndim*(iphase-1)+1:iphase*ndim,ndim*(iphase-1)+1:iphase*ndim))
                        !All the elements should be equal in the diagonal now
                        matrix(ndim*(iphase-1)+1:iphase*ndim,ndim*(iphase-1)+1:iphase*ndim) = perm%val(ndim*(iphase-1)+1:iphase*ndim,ndim*(iphase-1)+1:iphase*ndim,ele)/&
                            matrix(ndim*(iphase-1)+1,ndim*(iphase-1)+1)!Inverse

                        aux = matmul(matrix(ndim*(iphase-1)+1:iphase*ndim, ndim*(iphase-1)+1:iphase*ndim), &
                                    Velocity % val( :, iphase , u_inod) )
                        !P0 darcy velocities per element
                        DarcyVelocity % val ( :, iphase , u_inod) = DarcyVelocity% val ( :, iphase , u_inod)+ &
                             aux(:)*saturation%val(1,iphase, cv_loc)/ cv_nloc
                    end do
                end do
            end do
        end do

        deallocate(matrix)
        deallocate(aux)
    end subroutine get_DarcyVelocity




end module Copy_Outof_State




