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

module multiphase_fractures

    use quadrature
    use elements
    use sparse_tools
    use fields
    use fefields
    use state_module
    use copy_outof_state
    use spud
    use global_parameters, only: option_path_len, field_name_len, is_multifracture
    use futils, only: int2str
    use solvers
    use implicit_solids
    use FLDebug
    use memory_diagnostics
    use multi_tools
    implicit none

#ifdef USING_FEMDEM
    ! variable name convention:
    ! _r : ring, _v : volume, _vc : volume coarse

  interface
     subroutine y2d_allocate_femdem( string, &
          &     nodes_r, elements_r, edges_r, &
          &     nodes_v, elements_v, edges_v )
       character( len = * ), intent( in ) :: string
       integer, intent( out ) :: nodes_r, elements_r, edges_r, &
            &                    nodes_v, elements_v, edges_v
     end subroutine y2d_allocate_femdem
  end interface

  interface
     subroutine y2d_populate_femdem( ele1_r, ele2_r, ele3_r, &
          &                          face1_r, face2_r, x_r, y_r, &
          &                          p1, p2, p3, p4,  &
          &                          ele1_v, ele2_v, ele3_v, &
          &                          face1_v, face2_v, x_v, y_v )
       integer, dimension( * ), intent( out ) :: ele1_r, ele2_r, ele3_r, &
            &                                    ele1_v, ele2_v, ele3_v, &
            &                                    face1_r, face2_r, face1_v, face2_v
       real, dimension( * ), intent( out ) :: x_r, y_r, x_v, y_v, p1, p2, p3, p4
     end subroutine y2d_populate_femdem
  end interface



    ! Subroutine for coupling with viscosity (drag force, slip velocity, viscosity)
  interface
     subroutine y2dfemdem( string, dt, p, uf_r, vf_r, uf_v, vf_v, du_s, dv_s, u_s, v_s, &
         mu_f, f_x, f_y, usl, uvl, a_xx, a_xy, a_yy, p_v ) 
       character( len = * ), intent( in ) :: string
       real, intent( in ) :: dt
       real, dimension( * ), intent( in ) :: p, uf_r, vf_r, uf_v, vf_v, mu_f,p_v
       real, dimension( * ), intent( out ) :: du_s, dv_s, u_s, v_s, f_x, f_y, usl, uvl, a_xx, a_xy, a_yy
     end subroutine y2dfemdem
  end interface
#endif

    ! Soubroutine for coupling without viscosity
    !  interface
    !     subroutine y2dfemdem( string, dt, p, uf_r, vf_r, uf_v, vf_v, du_s, dv_s, u_s, v_s )
    !       character( len = * ), intent( in ) :: string
    !       real, intent( in ) :: dt
    !       real, dimension( * ), intent( in ) :: p, uf_r, vf_r, uf_v, vf_v
    !       real, dimension( * ), intent( out ) :: du_s, dv_s, u_s, v_s
    !     end subroutine y2dfemdem
    !  end interface
    !#endif



    type( vector_field ), save :: positions_r, positions_v, positions_vc
    type( tensor_field ), save :: permeability_r
    character( len = FIELD_NAME_LEN ), save :: femdem_mesh_name
    integer, save :: ndim

    private
    public :: fracking

contains


    !------------------ao---------------hydro-fracture/fracturing subroutine --------------------------------------------------
    subroutine fracking( packed_state, state, nphase )

        implicit none

        integer, intent( in ) :: nphase
        real, dimension( : ), allocatable :: p_r, muf_r, p_v
        real, dimension( :, : ), allocatable :: uf_r, uf_v, du_s, u_s, f ,u
        real, dimension( :, :, : ), allocatable :: a
        integer :: r_nonods, v_nonods
        real :: dt

        type( state_type ), intent( inout ) :: packed_state
        type( state_type ), dimension( : ), intent( inout ) :: state

        type( tensor_field), pointer :: darc_vel
        type( vector_field), pointer :: vel
        character( len = option_path_len ) :: opt

        ! read in ring and solid volume meshes
        ! and simplify the volume mesh
        call initialise_femdem

            vel => extract_vector_field( packed_state, "Darcy_Velocity" ) !!-ao darcy
            darc_vel => extract_tensor_field(packed_state,"PackedDarcyVelocity")
            vel%val=darc_vel%val(:,1,:)

        !!**************************************************************
        if ( have_option( '/simulation_type/femdem_fracture/oneway_coupling_only')) then!This option do not exist
            call get_option( '/simulation_type/femdem_fracture/oneway_coupling_only', opt )!This option do not exist
        !!**************************************************************
            if (trim( opt ) == "1way" ) then


                !calculate volume-fraction for mapping solid concentration
                call calculate_volume_fraction( packed_state)
                ! calculate porosity and permeability
                call calculate_phi_and_perm( packed_state, state )
                ! deallocate
                call deallocate_femdem

            else  !! pseudo-2-way-coupling so that fluid is solved of varying geomtery and aperture
                !!-ao two way coupling
                r_nonods = node_count( positions_r )
                v_nonods = node_count( positions_v )

                allocate( p_r( r_nonods ), uf_r( ndim, r_nonods ), muf_r( r_nonods ), &
                    f( ndim, r_nonods ),   u(ndim, r_nonods) , a( ndim, ndim, r_nonods), &
                    uf_v( ndim, v_nonods ) , du_s( ndim, v_nonods ), u_s( ndim, v_nonods ), p_v(v_nonods)  )
                p_r=0.0 ; uf_r=0.0 ; muf_r=1.0 ; f=0.0 ; a=0.0 ; uf_v=0.0 ; du_s=0.0 ; u_s=0.0; p_v=0.0;

                call get_option( "/timestepping/timestep", dt )
                print *, "pressure:", maxval(p_r), minval(p_r) !!-ao

                call y2dfemdem( trim( femdem_mesh_name ) // char( 0 ), dt, p_r, uf_r( 1, : ), uf_r( 2, : ), &
                    uf_v( 1, : ), uf_v( 2, : ), du_s( 1, : ), du_s( 2, : ), u_s( 1, : ), u_s( 2, : ), &
                    muf_r, f( 1, : ), f( 2, : ), u(1, : ), u(2, : ), a( 1, 1, : ), a( 1, 2, : ), a( 2 ,2 , : ),p_v )

                !calculate volume-fraction for mapping solid concentration
                call calculate_volume_fraction( packed_state)

                ! calculate porosity and permeability
                call calculate_phi_and_perm( packed_state, state )


                ! deallocate
                call deallocate_femdem
                deallocate( p_r, uf_r, muf_r, uf_v, du_s, u_s, f, a, u, p_v)
            end if

        else
            !!-ao two way coupling
            r_nonods = node_count( positions_r )
            v_nonods = node_count( positions_v )

            allocate( p_r( r_nonods ), uf_r( ndim, r_nonods ), muf_r( r_nonods ), &
            f( ndim, r_nonods ),   u(ndim, r_nonods) , a( ndim, ndim, r_nonods), &
            uf_v( ndim, v_nonods ) , du_s( ndim, v_nonods ), u_s( ndim, v_nonods ), p_v(v_nonods) )
            p_r=0.0 ; uf_r=0.0 ; muf_r=1.0 ; f=0.0 ; a=0.0 ; uf_v=0.0 ; du_s=0.0 ; u_s=0.0; p_v=0.0;

!           !interpolate presure, velocity and visc from fluid to solid through ring
!            call interpolate_fields_out_r( packed_state, nphase, p_r, uf_r, muf_r ) !!-ao is this causing the problem with large negative pressure
            call interpolate_fields_out_r_p( packed_state, nphase, p_r) !-ao!use for Darcy flow two-way coupling if problems arise in Y_Drag (solid)


            vel => extract_vector_field( packed_state, "Darcy_Velocity" )
            darc_vel => extract_tensor_field(packed_state, "PackedDarcyVelocity")
            vel%val=darc_vel%val(:,1,:)

            call get_option( "/timestepping/timestep", dt )


           !interpolate pressure locally in solid to calculate local stresses due to pore fluid pressure
             if ( have_option( '/simulation_type/femdem_fracture/include_pore_pressure')) then !with pore_fluid presure!This option do not exist

                call interpolate_fields_out_v_pf( packed_state, nphase, p_v)



                call y2dfemdem( trim( femdem_mesh_name ) // char( 0 ), dt, p_r, uf_r( 1, : ), uf_r( 2, : ), &
                    uf_v( 1, : ), uf_v( 2, : ), du_s( 1, : ), du_s( 2, : ), u_s( 1, : ), u_s( 2, : ), &
                    muf_r, f( 1, : ), f( 2, : ), u(1, : ), u(2, : ), a( 1, 1, : ), a( 1, 2, : ), a( 2 ,2 , : ), p_v)
                !calculate volume-fraction for mapping solid concentration
                call calculate_volume_fraction( packed_state)

                ! calculate porosity and permeability
                call calculate_phi_and_perm( packed_state, state )

                ! interpolate from solid to fluid through volume mesh
 !               call interpolate_fields_in_v( packed_state, du_s, u_s )

                ! deallocate
                call deallocate_femdem
                deallocate( p_r, uf_r, muf_r, uf_v, du_s, u_s, f, a, u, p_v)

             else !without pore_fluid pressure

                print *, "-----WARNING: PORE PRESSURE IS NOT INCLUDED!------"

                call y2dfemdem( trim( femdem_mesh_name ) // char( 0 ), dt, p_r, uf_r( 1, : ), uf_r( 2, : ), &
                uf_v( 1, : ), uf_v( 2, : ), du_s( 1, : ), du_s( 2, : ), u_s( 1, : ), u_s( 2, : ), &
                muf_r, f( 1, : ), f( 2, : ), u(1, : ), u(2, : ), a( 1, 1, : ), a( 1, 2, : ), a( 2 ,2 , : ),p_v )

                !calculate volume-fraction for mapping solid concentration
                call calculate_volume_fraction( packed_state)

                ! calculate porosity and permeability
                call calculate_phi_and_perm( packed_state, state )

                ! interpolate from solid to fluid through volume mesh
    !            call interpolate_fields_in_v( packed_state, du_s, u_s )

                ! deallocate
                call deallocate_femdem
                deallocate( p_r, uf_r, muf_r, uf_v, du_s, u_s, f, a, u, p_v)

             end if

        endif
        !!**************************************************************

        return
    end subroutine fracking

    !----------------------------------------------------------------------------------------------------------
    subroutine initialise_femdem

    implicit none

        integer :: i, loc, sloc
        integer :: nodes_r, elements_r, edges_r, &
        nodes_v, elements_v, edges_v
        integer, dimension( : ), allocatable :: ele1_r, ele2_r, ele3_r, &
        ele1_v, ele2_v, ele3_v, &
        face1_r, face2_r, face1_v, face2_v
        real, dimension( : ), allocatable :: x_r, y_r, x_v, y_v, p1, p2, p3, p4
        type( quadrature_type ) :: quad
        type( element_type ) :: shape
        integer, dimension( : ), allocatable :: sndglno_r, boundary_ids_r, &
        sndglno_v, boundary_ids_v
        integer :: quad_degree, poly_degree, continuity
        type( mesh_type ) :: mesh_r, mesh_v, mesh_r_p0

    ewrite(3,*) "inside initialise_femdem"

    if (have_option('/simulation_type/femdem_fracture') ) then
        call get_option( "/simulation_type/femdem_fracture/femdem_file/name", femdem_mesh_name )!This option do not exist
        femdem_mesh_name = trim( femdem_mesh_name ) // ".y"
    end if

    call get_option( "/geometry/quadrature/degree", quad_degree )
    call get_option( "/geometry/dimension", ndim )

    if ( ndim == 2 ) then
       loc = 3 ; sloc= 2
    else if ( ndim == 3 ) then
       loc = 4 ; sloc= 3
       FLAbort( "Fracture modelling is supported for 2D only." )
    end if

    call y2d_allocate_femdem( trim( femdem_mesh_name ) // char( 0 ), &
         nodes_r, elements_r, edges_r, nodes_v, elements_v, edges_v )

    ewrite(3,*) "nodes_r, elements_r, edges_r, nodes_v, elements_v, edges_v", &
                 nodes_r, elements_r, edges_r, nodes_v, elements_v, edges_v

    allocate( ele1_r( elements_r ), ele2_r( elements_r ), ele3_r( elements_r ) )
    allocate( face1_r( edges_r ), face2_r( edges_r ) )

    allocate( ele1_v( elements_v ), ele2_v( elements_v ), ele3_v( elements_v ) )
    allocate( face1_v( edges_v ), face2_v( edges_v ) )

    allocate( x_r( nodes_r ), y_r( nodes_r ) )
    allocate( x_v( nodes_v ), y_v( nodes_v ) )

        allocate( p1( elements_r ), p2( elements_r ) )
        allocate( p3( elements_r ), p4( elements_r ) )



        call y2d_populate_femdem( ele1_r, ele2_r, ele3_r, &
        face1_r, face2_r, x_r, y_r, p1, p2, p3, p4, &
        ele1_v, ele2_v, ele3_v, face1_v, face2_v, x_v, y_v )

print *, "passed populate here" !!-ao

    quad = make_quadrature( loc, ndim, degree = quad_degree )
    shape = make_element_shape( loc, ndim, 1, quad )

    ! create the ring mesh
    call allocate( mesh_r, nodes_r, elements_r, shape, name="CoordinateMesh" )
    call allocate( positions_r, ndim, mesh_r, name="Coordinate" )

    positions_r%val( 1, : ) = x_r
    positions_r%val( 2, : ) = y_r

    do i = 1, elements_r
       positions_r%mesh%ndglno( (i-1)*loc+1 : i*loc ) = &
            (/ ele1_r(i)+1, ele2_r(i)+1, ele3_r(i)+1 /)
    end do

    allocate( sndglno_r( edges_r * sloc ) ) ; sndglno_r = 0
    allocate( boundary_ids_r( edges_r ) ) ; boundary_ids_r = 666

    do i = 1, edges_r
       sndglno_r( (i-1)*sloc+1 : i*sloc ) = &
            (/ face1_r(i)+1, face2_r(i)+1 /)
    end do

    call add_faces( positions_r%mesh, &
         sndgln = sndglno_r, &
         boundary_ids = boundary_ids_r )

    positions_r%dim = ndim

    deallocate( boundary_ids_r, sndglno_r )
    call deallocate( mesh_r )

    ! create the volume mesh
    call allocate( mesh_v, nodes_v, elements_v, shape, name="CoordinateMesh" )
    call allocate( positions_v, ndim, mesh_v, name="Coordinate" )
    positions_v%val( 1, : ) = x_v
    positions_v%val( 2, : ) = y_v

    do i = 1, elements_v
       positions_v%mesh%ndglno( (i-1)*loc+1 : i*loc ) = &
            (/ ele1_v(i)+1, ele2_v(i)+1, ele3_v(i)+1 /)
    end do

    allocate( sndglno_v( edges_v * sloc ) ) ; sndglno_v = 0
    allocate( boundary_ids_v( edges_v ) ) ; boundary_ids_v = 666

    do i = 1, edges_v
       sndglno_v( (i-1)*sloc+1 : i*sloc ) = &
            (/ face1_v(i)+1, face2_v(i)+1 /)
    end do

    call add_faces( positions_v%mesh, &
         sndgln = sndglno_v, &
         boundary_ids = boundary_ids_v )

    positions_v%dim = ndim

    deallocate( boundary_ids_v, sndglno_v )
    call deallocate( mesh_v )

    call deallocate_element( shape )

    ! coarsen the volume mesh
    call coarsen_mesh_2d( positions_v, positions_vc )

    ! create the ring p0 mesh
    poly_degree = 0 ; continuity = -1
    shape = make_element_shape( loc, ndim, poly_degree, quad )
    mesh_r_p0 = make_mesh( positions_r%mesh, shape, continuity, "P0DG" )

    ! store ring permeability
    if (associated(permeability_r%val)) then
        call deallocate(permeability_r)
    end if
    call allocate( permeability_r, mesh_r_p0, name = "Permeability" )
    call zero( permeability_r )

        call set_all( permeability_r, 1, 1, p1 )
        call set_all( permeability_r, 1, 2, p2 )
        call set_all( permeability_r, 2, 1, p3 )
        call set_all( permeability_r, 2, 2, p4 )

        ! deallocate
        call deallocate( mesh_r_p0 )
        call deallocate_element( shape )
        call deallocate( quad )

        deallocate( ele1_r, ele2_r, ele3_r, face1_r, face2_r, &
        ele1_v, ele2_v, ele3_v, face1_v, face2_v, &
        x_r, y_r, x_v, y_v, p1, p2, p3, p4)

        ewrite(3,*) "leaving initialise_femdem"

        return
    end subroutine initialise_femdem

    !----------------------------------------------------------------------------------------------------------

    subroutine calculate_volume_fraction( packed_state )

        implicit none

        type( state_type ), intent( inout ) :: packed_state
        !Local variables
        type( state_type ) :: alg_ext, alg_fl
        type( mesh_type ), pointer :: fl_mesh, p0_fl_mesh
        type( vector_field ), pointer :: fl_positions
        type( scalar_field ), pointer :: volume_fraction

        ewrite(3,*) "inside calculate_volume_fraction"

        call insert( alg_ext, positions_vc % mesh, "Mesh" )
        call insert( alg_ext, positions_vc, "Coordinate" )

        fl_mesh => extract_mesh( packed_state, "CoordinateMesh" )
        fl_positions => extract_vector_field( packed_state, "Coordinate" )

        call insert( alg_fl, fl_mesh, "Mesh" )
        call insert( alg_fl, fl_positions, "Coordinate" )


        ! volume fraction, i.e. porosity...
        volume_fraction => extract_scalar_field( packed_state, "SolidConcentration" )
        call zero( volume_fraction )

        ewrite(3,*) "...interpolating"

        ! interpolate
        call interpolation_galerkin_femdem( alg_ext, alg_fl, field = volume_fraction )


        ! ensure vf is between 0.0 and 1.0
        call bound_volume_fraction( volume_fraction%val )

        ! deallocate
!        call deallocate(fl_mesh)
!        call deallocate( fl_positions )
        call deallocate( alg_fl )
        call deallocate( alg_ext )

        ewrite(3,*) "leaving calculate_volume_fraction"

        return
    end subroutine calculate_volume_fraction

    !----------------------------------------------------------------------------------------------------------

!  subroutine calculate_absorption( totele, cv_nloc, cv_nonods, &
!       &                           nphase, cv_ndgln, rho, vf, dt, absorption )
!
!    implicit none
!
!    integer, intent( in ) :: totele, cv_nloc, cv_nonods, nphase
!    integer, dimension( : ), intent( in ) :: cv_ndgln
!    real, intent( in ) :: dt
!    real, dimension( : ), intent( in ) :: vf
!    real, dimension( :), intent( in ) :: rho
!
!    real, dimension( :, :, : ), intent( inout ) :: absorption
!
!    integer :: ele, iloc, mi, ci, iphase, idim, idx
!    real :: sigma
!
!    ewrite(3,*) "inside calculate_absorption"
!
!    do ele = 1, totele
!       do iloc = 1, cv_nloc
!
!          mi = ( ele - 1 ) * cv_nloc + iloc
!          ci = cv_ndgln( ( ele - 1 ) * cv_nloc + iloc )
!
!          do iphase = 1, nphase
!             sigma = vf( ele ) * rho( ci + ( iphase - 1 ) * cv_nonods ) / dt
!             do idim = 1, ndim
!                idx = idim + ( iphase - 1 ) * ndim
!                absorption( mi, idx, idx ) = absorption( mi, idx, idx ) + sigma
!             end do
!          end do
!
!       end do
!    end do
!
!    ewrite(3,*) "leaving calculate_absorption"
!
!    return
!  end subroutine calculate_absorption

  !----------------------------------------------------------------------------------------------------------
    subroutine calculate_phi_and_perm( packed_state, state )

        implicit none

        type( state_type ), intent( inout ) :: packed_state
        type( state_type ), dimension( : ), intent( inout ) :: state

        !Local variables
        type( state_type ) :: alg_ext, alg_fl
        type( mesh_type ), pointer :: fl_mesh, p0_fl_mesh
        type( scalar_field ) :: rvf
        type(vector_field):: fl_positions2

        type( vector_field ), pointer :: fl_positions, dealloc, porosity
        type( tensor_field ), pointer :: permeability, perm_state
        type( scalar_field ), pointer :: vf , poro_state, dum, perm_val, perm2_val
        real, dimension( :, :, : ), allocatable :: perm


        type( scalar_field ) :: field_fl_p11, field_fl_p12, field_fl_p21, field_fl_p22, &
        field_ext_p11, field_ext_p12, field_ext_p21, field_ext_p22

        character( len = OPTION_PATH_LEN ) :: &
        path = "/tmp/galerkin_projection/continuous"
        integer :: stat, totele, ele
        real, dimension( : ), allocatable :: scale
        real, parameter :: tol = 1.0e-10
        real :: bg_poro, bg_perm

        ewrite(3,*) "inside calculate_phi_and_perm"

        call insert( alg_ext, positions_r%mesh, "Mesh" )
        call insert( alg_ext, positions_r, "Coordinate" ) 

        fl_mesh => extract_mesh( packed_state, "CoordinateMesh" )
        fl_positions => extract_vector_field( packed_state, "Coordinate" )

        call insert( alg_fl, fl_mesh, "Mesh" )
        call insert( alg_fl, fl_positions, "Coordinate" ) 


        call set_solver_options( path, &
        ksptype = "gmres", &
        pctype = "hypre", &
        rtol = 1.0e-10, &
        atol = 1.0e-15, &
        max_its = 10000 )
        call add_option( &
        trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", stat)
        call set_option( &
        trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", "boomeramg")

        path = "/tmp"

        p0_fl_mesh => extract_mesh( packed_state, "P0DG" )

        ewrite(3,*) "...generating fluids state"

        ! this is the permeability on the fluidity mesh
        call allocate( field_fl_p11, p0_fl_mesh, "Permeability11" )
        call zero( field_fl_p11 )
        call insert( alg_fl, field_fl_p11, "Permeability11" )
        field_fl_p11 % option_path = path

        call allocate( field_fl_p12, p0_fl_mesh, "Permeability12" )
        call zero( field_fl_p12 )
        call insert( alg_fl, field_fl_p12, "Permeability12" )
        field_fl_p12 % option_path = path

        call allocate( field_fl_p21, p0_fl_mesh, "Permeability21" )
        call zero( field_fl_p21 )
        call insert( alg_fl, field_fl_p21, "Permeability21" )
        field_fl_p21 % option_path = path

        call allocate( field_fl_p22, p0_fl_mesh, "Permeability22" )
        call zero( field_fl_p22 )
        call insert( alg_fl, field_fl_p22, "Permeability22" )
        field_fl_p22 % option_path = path

        ewrite(3,*) "...generating femdem/ring state"

        ! this is the permeability on the solid mesh
        field_ext_p11 = extract_scalar_field( permeability_r, 1, 1 )
        call insert( alg_ext, field_ext_p11, "Permeability11" )

        field_ext_p12 = extract_scalar_field( permeability_r, 1, 2 )
        call insert( alg_ext, field_ext_p12, "Permeability12" )

        field_ext_p21 = extract_scalar_field( permeability_r, 2, 1 )
        call insert( alg_ext, field_ext_p21, "Permeability21" )

        field_ext_p22 = extract_scalar_field( permeability_r, 2, 2 )
        call insert( alg_ext, field_ext_p22, "Permeability22" )

        ! volume fraction - this is the ring
        call allocate( rvf, p0_fl_mesh, "VolumeFraction" )
        call zero( rvf )

        ewrite(3,*) "...interpolating"

        ! interpolate through ring-mesh
        call interpolation_galerkin_femdem( alg_ext, alg_fl, field = rvf )

        ! bound ring volume fraction
        call bound_volume_fraction( rvf % val )

        ! calculate phi and perm
        ! only fractures are permeable and porous
        ! perm is non-zero only in fractures
        ! also add an isotropic, background/matrix permeability
        permeability => extract_tensor_field( packed_state, "Permeability" )
        totele=ele_count(fl_mesh)


      call get_option("/porous_media/tensor_field::Permeability/prescribed/value::WholeMesh/isotropic/constant", bg_perm)
      call get_option("/porous_media/scalar_field::Porosity/prescribed/value::WholeMesh/constant", bg_poro)



    allocate( perm(ndim, ndim, totele) ) ; perm= 0.0
        do ele=1, totele
            if ((rvf % val (ele) > 0.0) .AND. ((field_fl_p11%val(ele)>bg_perm) .OR. (field_fl_p22%val(ele)>bg_perm))) then
                  ! non-normalised and conservative permeability interpolation
                  perm( 1, 1, ele ) = field_fl_p11 % val (ele)
                  perm( 1, 2, ele ) = field_fl_p12 % val (ele)
                  perm( 2, 1, ele ) = field_fl_p21 % val (ele)
                  perm( 2, 2, ele ) = field_fl_p22 % val (ele)
            else
                perm( 1, 1, ele ) =permeability%val(1,1,ele)
                perm( 1, 2, ele ) =permeability%val(1,2,ele)
                perm( 2, 1, ele ) =permeability%val(2,1,ele)
                perm( 2, 2, ele ) =permeability%val(2,2,ele)
            end if
        end do




        ! assign permeability from interpolation to the memory
        permeability % val( 1, 1, : ) =   perm( 1, 1, : )
        permeability % val( 1, 2, : ) =   perm( 1, 2, : )
        permeability % val( 2, 1, : ) =   perm( 2, 1, : )
        permeability % val( 2, 2, : ) =   perm( 2, 2, : )

        ! extract porosity and solid concentration fields from state
        porosity => extract_vector_field( packed_state, "Porosity" )
        vf => extract_scalar_field( packed_state, "SolidConcentration" )
	
	! PRINT *, "size of porosity is ", size(porosity % val)    ====== 1


        ! for adaptivity (bound perm field)
        perm_val => extract_scalar_field( state(1), "Dummy" )
        allocate(perm_val%val(totele))
        call zero( perm_val)

!        !visualising permeability in 'totalflux'Dummy field
        perm2_val => extract_scalar_field( state(1), "TotalFlux" )
        allocate(perm2_val%val(totele))
        call zero( perm2_val)
!!-ao comment - porosity is not scaled due to problems arising in the wall
!               where porosities (rvf) can arise lower than background porosity
        do ele = 1, totele
            if (rvf % val (ele) > 0.0) then
                 porosity % val (:,ele) = bg_poro*(1-rvf % val(ele))+ rvf % val (ele) ! calcualtion of effective phi --->  Phi_bg*(1-rvf_ring)+1*(rvf_ring)
            endif

            if (rvf%val(ele)> 0.0) perm_val % val (ele)=1
            if ( maxval( permeability % val( :, :, ele ) ) <= bg_perm ) perm_val % val (ele) = 0 !     ! for adaptivity (bound porosity field)

            perm2_val % val (ele) = maxval( permeability % val( :, :, ele ) ) !        !visualising permeability in 'totalflux'Dummy field


        end do

        call bound_volume_fraction( vf%val )


        ! deallocate
        deallocate( perm)
        call deallocate( rvf )

        call deallocate( alg_ext )
        call deallocate( alg_fl )


        ewrite(3,*) "leaving calculate_phi_and_perm"

        return
    end subroutine calculate_phi_and_perm

    !----------------------------------------------------------------------------------------------------------

    subroutine interpolate_fields_out_r( packed_state, nphase, p_r, u_r, mu_r )

        implicit none

        type( state_type ), intent( in ) :: packed_state
        integer, intent( in ) :: nphase
        real, dimension( : ), intent( inout ) :: p_r, mu_r
        real, dimension( :, : ), intent( inout ) :: u_r

        !Local variables
        type( mesh_type ), pointer :: fl_mesh, u_mesh, p0_fl_mesh
        type( scalar_field ) :: field_fl_p, field_fl_u, field_fl_v, field_fl_mu, &
        field_ext_p, field_ext_u, field_ext_v, field_ext_mu, &
        u_dg, v_dg, rvf
       !! type( scalar_field ), pointer :: pressure
        type( vector_field ), pointer :: fl_positions, vel
        type( tensor_field ), pointer :: velocity, viscosity, pressure, darc_vel
        type( state_type ) :: alg_ext, alg_fl
        real, dimension( :, :, : ), allocatable :: u_tmp
        real, dimension( :, :), allocatable :: u_tmpv
        integer, dimension( : ), pointer :: fl_ele_nodes, cv_ndgln
        integer :: ele, totele, u_nloc, cv_nloc, u_nonods, &
        stat
        logical :: constant_mu
        character( len = OPTION_PATH_LEN ) :: path = "/tmp/galerkin_projection/continuous"

        p_r = 0.0 ; u_r = 0.0 ; mu_r = 0.0

        cv_ndgln => get_ndglno( extract_mesh( packed_state, "PressureMesh" ) )

        pressure => extract_tensor_field( packed_state, "PackedFEPressure" )

        cv_nloc = ele_loc( pressure, 1 )

        fl_mesh => extract_mesh( packed_state, "CoordinateMesh" )
        fl_positions => extract_vector_field( packed_state, "Coordinate" )

        viscosity => extract_tensor_field( packed_state, "Viscosity", stat )

        !have_viscosity = ( stat == 0 )

        constant_mu = .false.
        if ( is_constant( viscosity ) ) constant_mu = .true.

        totele = ele_count( fl_mesh )

        call set_solver_options( path, &
        ksptype = "gmres", &
        pctype = "hypre", &
        rtol = 1.0e-10, &
        atol = 0.0, &
        max_its = 10000 )
        call add_option( &
        trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", stat)
        call set_option( &
        trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", "boomeramg")

        path = "/tmp"

        call allocate( field_fl_p, fl_mesh, "Pressure" )
        call zero( field_fl_p )

        call allocate( field_fl_u, fl_mesh, "Velocity1" )
        call zero( field_fl_u )

        call allocate( field_fl_v, fl_mesh, "Velocity2" )
        call zero( field_fl_v )

        call allocate( field_fl_mu, fl_mesh, "Viscosity" )
        call zero( field_fl_mu )

        ! deal with pressure and viscosity
        if ( cv_nloc == 6 ) then
            ! linearise pressure for p2
            do ele = 1, totele
                fl_ele_nodes => ele_nodes( fl_mesh, ele )
                field_fl_p % val( fl_ele_nodes( 1 ) ) = pressure % val( 1,1, cv_ndgln( ( ele - 1 ) * cv_nloc + 1 ) )
                field_fl_p % val( fl_ele_nodes( 2 ) ) = pressure % val( 1,1, cv_ndgln( ( ele - 1 ) * cv_nloc + 3 ) )
                field_fl_p % val( fl_ele_nodes( 3 ) ) = pressure % val( 1,1, cv_ndgln( ( ele - 1 ) * cv_nloc + 6 ) )


                if ( constant_mu ) then
                    field_fl_mu % val(fl_ele_nodes( 1 ) ) = viscosity % val( 1, 1, 1 )
                    field_fl_mu % val(fl_ele_nodes( 2 ) ) = viscosity % val( 1, 1, 1 )
                    field_fl_mu % val(fl_ele_nodes( 3 ) ) = viscosity % val( 1, 1, 1 )
                else
                    field_fl_mu % val(fl_ele_nodes( 1 ) ) = viscosity % val( 1, 1, cv_ndgln( ( ele - 1 ) * cv_nloc + 1 ) )
                    field_fl_mu % val(fl_ele_nodes( 2 ) ) = viscosity % val( 1, 1, cv_ndgln( ( ele - 1 ) * cv_nloc + 3 ) )
                    field_fl_mu % val(fl_ele_nodes( 3 ) ) = viscosity % val( 1, 1, cv_ndgln( ( ele - 1 ) * cv_nloc + 6 ) )
                end if
            end do
        else
            ! just copy memory for p1
            field_fl_p % val = pressure % val (1,1,:)
!            field_fl_mu % val = viscosity % val( 1, 1, : )
             field_fl_mu % val = viscosity % val( 1, 1, 1 )
        end if


        if ( is_multifracture  ) then
            ! deal with velocity - this part needs optimisation...
            vel => extract_vector_field( packed_state, "Darcy_Velocity" ) !!-ao darcy
            darc_vel => extract_tensor_field(packed_state, "PackedDarcyVelocity")
            vel%val=darc_vel%val(:,1,:)
           
	    u_nonods = node_count( vel )
            u_nloc = ele_loc( vel, 1 )
            allocate( u_tmpv( ndim, u_nonods ) )
            u_tmpv = vel % val

        endif

        u_mesh => extract_mesh( packed_state, "VelocityMesh" )
        call allocate( u_dg, u_mesh, "u_dg" )
        call zero( u_dg )
        call allocate( v_dg, u_mesh, "v_dg" )
        call zero( v_dg )

        if  ( is_multifracture  ) then
            u_dg % val = u_tmpv( 1, : ) ; v_dg % val = u_tmpv( 2, : ) !!-ao darcy
        endif


        !call project_field( u_dg, field_fl_u, fl_positions )
        !call project_field( v_dg, field_fl_v, fl_positions )

        call linear2quadratic_field( u_dg, field_fl_u)
        call linear2quadratic_field( v_dg, field_fl_v)

        ! fluidity state
        call insert( alg_fl, fl_mesh, "Mesh" )
        call insert( alg_fl, fl_positions, "Coordinate" )

        call insert( alg_fl, field_fl_p, "Pressure" )
        call insert( alg_fl, field_fl_u, "Velocity1" )
        call insert( alg_fl, field_fl_v, "Velocity2" )
        call insert( alg_fl, field_fl_mu, "Viscosity" )

        ! ring state
        call insert( alg_ext, positions_r%mesh, "Mesh" )
        call insert( alg_ext, positions_r, "Coordinate" )

        call allocate( field_ext_p, positions_r%mesh, "Pressure" )
        call zero( field_ext_p )
        field_ext_p % option_path = path
        call insert( alg_ext, field_ext_p, "Pressure" )

        call allocate( field_ext_u, positions_r%mesh, "Velocity1" )
        call zero( field_ext_u )
        field_ext_u % option_path = path
        call insert( alg_ext, field_ext_u, "Velocity1" )

        call allocate( field_ext_v, positions_r%mesh, "Velocity2" )
        call zero( field_ext_v )
        field_ext_v % option_path = path
        call insert( alg_ext, field_ext_v, "Velocity2" )

        call allocate( field_ext_mu, positions_r%mesh, "Viscosity" )
        call zero( field_ext_mu )
        field_ext_mu % option_path = path
        call insert( alg_ext, field_ext_mu, "Viscosity" )

        ewrite(3,*) "...interpolating"

        ! interpolate
        call interpolation_galerkin_femdem( alg_fl, alg_ext, femdem_out = .true. ) !!-ao problem here...

        ! copy memory
        p_r = field_ext_p % val
        u_r( 1, : ) = field_ext_u % val
        u_r( 2, : ) = field_ext_v % val
        mu_r = field_ext_mu % val



        ! deallocate
        call deallocate( field_fl_p )
        call deallocate( field_fl_u )
        call deallocate( field_fl_v )
        call deallocate( field_fl_mu )


        call deallocate( field_ext_p )
        call deallocate( field_ext_u )
        call deallocate( field_ext_v )
        call deallocate( field_ext_mu )


        call deallocate( alg_fl )
        call deallocate( alg_ext )

        call deallocate( u_dg )
        call deallocate( v_dg )

        if (is_multifracture) then
            deallocate( u_tmpv )
        else
            deallocate( u_tmp )
        end if

        return

    end subroutine interpolate_fields_out_r
    !----------------------------------------------------------------------------------------------------------

 subroutine interpolate_fields_out_v( packed_state, nphase, u_v )

    implicit none

    type( state_type ), intent( in ) :: packed_state
    integer, intent( in ) :: nphase
    real, dimension( :, : ), intent( inout ) :: u_v

    !Local variables
    type( mesh_type ), pointer :: fl_mesh, u_mesh
    type( scalar_field ) :: field_fl_u, field_fl_v, &
         &                  field_ext_u, field_ext_v, &
         &                  u_dg, v_dg
    type( vector_field ), pointer :: fl_positions
    type( tensor_field ), pointer :: velocity
    type( state_type ) :: alg_ext, alg_fl
    real, dimension( :, :, : ), allocatable :: u_tmp
    integer :: stat, u_nonods
    character( len = OPTION_PATH_LEN ) :: path = "/tmp/galerkin_projection/continuous"

    u_v = 0.0

    fl_mesh => extract_mesh( packed_state, "CoordinateMesh" )
    fl_positions => extract_vector_field( packed_state, "Coordinate" )

    call set_solver_options( path, &
         ksptype = "gmres", &
         pctype = "hypre", &
         rtol = 1.0e-10, &
         atol = 0.0, &
         max_its = 10000 )
    call add_option( &
         trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", stat)
    call set_option( &
         trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", "boomeramg")

    path = "/tmp"

    call allocate( field_fl_u, fl_mesh, "Velocity1" )
    call zero( field_fl_u )

    call allocate( field_fl_v, fl_mesh, "Velocity2" )
    call zero( field_fl_v )

    ! deal with velocity - this part needs optimisation...
    velocity => extract_tensor_field( packed_state, "PackedVelocity" )
    u_nonods = node_count( velocity )

    allocate( u_tmp( ndim, nphase, u_nonods ) )
    u_tmp = velocity % val

    u_mesh => extract_mesh( packed_state, "VelocityMesh" )

    call allocate( u_dg, u_mesh, "u_dg" )
    call zero( u_dg )
    call allocate( v_dg, u_mesh, "v_dg" )
    call zero( v_dg )
    u_dg % val = u_tmp( 1, 1, : ) ; v_dg % val = u_tmp( 2, 1, : ) 


    !call project_field( u_dg, field_fl_u, fl_positions )
    !call project_field( v_dg, field_fl_v, fl_positions )

    call linear2quadratic_field( u_dg, field_fl_u)
    call linear2quadratic_field( v_dg, field_fl_v)




    ! fluidity state
    call insert( alg_fl, fl_mesh, "Mesh" )
    call insert( alg_fl, fl_positions, "Coordinate" )

    call insert( alg_fl, field_fl_u, "Velocity1" )
    call insert( alg_fl, field_fl_v, "Velocity2" )

    ! ring state
    call insert( alg_ext, positions_v%mesh, "Mesh" )
    call insert( alg_ext, positions_v, "Coordinate" )

    call allocate( field_ext_u, positions_v%mesh, "Velocity1" )
    call zero( field_ext_u )
    field_ext_u % option_path = path
    call insert( alg_ext, field_ext_u, "Velocity1" )

    call allocate( field_ext_v, positions_v%mesh, "Velocity2" )
    call zero( field_ext_v )
    field_ext_v % option_path = path
    call insert( alg_ext, field_ext_v, "Velocity2" )

    ewrite(3,*) "...interpolating"

    ! interpolate
    call interpolation_galerkin_femdem( alg_fl, alg_ext, femdem_out = .true. )

    ! copy memory
    u_v( 1, : ) = field_ext_u % val
    u_v( 2, : ) = field_ext_v % val

    ! deallocate
    call deallocate( field_fl_u )
    call deallocate( field_fl_v )

    call deallocate( field_ext_u )
    call deallocate( field_ext_v )

    call deallocate( alg_fl )
    call deallocate( alg_ext )

    call deallocate( u_dg )
    call deallocate( v_dg )

    deallocate( u_tmp )

    return
  end subroutine interpolate_fields_out_v

  !----------------------------------------------------------------------------------------------------------

  subroutine interpolate_fields_in_v( packed_state, du_s, u_s )

    implicit none

    type( state_type ), intent( in ) :: packed_state
    real, dimension( :, : ), intent( in ) :: du_s, u_s

    !Local variables
    type( mesh_type ), pointer :: p0_fl_mesh, fl_mesh, u_mesh
    type( scalar_field ), pointer :: solid, old_solid, f
    type( scalar_field ) :: field_fl_du, field_fl_dv, &
         &                  field_fl_us, field_fl_vs, &
         &                  field_ext_du, field_ext_dv, &
         &                  field_ext_us, field_ext_vs, f2, &
         &                  field_fl_solid, field_ext_solid
    type( vector_field ), pointer :: fl_positions, delta_u, solid_u
    type( state_type ) :: alg_ext, alg_fl
    integer :: stat, idim
    character( len = OPTION_PATH_LEN ) :: path = "/tmp/galerkin_projection/continuous"

    call set_solver_options( path, &
         ksptype = "gmres", &
         pctype = "hypre", &
         rtol = 1.0e-10, &
         atol = 0.0, &
         max_its = 10000 )
    call add_option( &
         trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", stat)
    call set_option( &
         trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", "boomeramg")

    path = "/tmp"

    p0_fl_mesh => extract_mesh( packed_state, "P0DG" )
    fl_mesh => extract_mesh( packed_state, "CoordinateMesh" )
    fl_positions => extract_vector_field( packed_state, "Coordinate" )

    ! fluidity state
    call insert( alg_fl, fl_mesh, "Mesh" )
    call insert( alg_fl, fl_positions, "Coordinate" )

    call allocate( field_fl_du, fl_positions%mesh, "deltaVelocity1" )
    call zero( field_fl_du )
    field_fl_du % option_path = path
    call insert( alg_fl, field_fl_du, "deltaVelocity1" )

    call allocate( field_fl_dv, fl_positions%mesh, "deltaVelocity2" )
    call zero( field_fl_dv )
    field_fl_dv % option_path = path
    call insert( alg_fl, field_fl_dv, "deltaVelocity2" )

    call allocate( field_fl_us, fl_positions%mesh, "SolidVelocity1" )
    call zero( field_fl_us )
    field_fl_us % option_path = path
    call insert( alg_fl, field_fl_us, "SolidVelocity1" )

    call allocate( field_fl_vs, fl_positions%mesh, "SolidVelocity2" )
    call zero( field_fl_vs )
    field_fl_vs % option_path = path
    call insert( alg_fl, field_fl_vs, "SolidVelocity2" )

    call allocate( field_fl_solid, fl_positions%mesh, "SolidConcentration" )
    call zero( field_fl_solid )
    field_fl_solid % option_path = path
    call insert( alg_fl, field_fl_solid, "SolidConcentration" )


    ! solid volume state
    call insert( alg_ext, positions_v%mesh, "Mesh" )
    call insert( alg_ext, positions_v, "Coordinate" )

    call allocate( field_ext_du, positions_v%mesh, "deltaVelocity1" )
    field_ext_du % val = du_s( 1, : )
    call insert( alg_ext, field_ext_du, "deltaVelocity1" )

    call allocate( field_ext_dv, positions_v%mesh, "deltaVelocity2" )
    field_ext_dv % val = du_s( 2, : )
    call insert( alg_ext, field_ext_dv, "deltaVelocity2" )

    call allocate( field_ext_us, positions_v%mesh, "SolidVelocity1" )
    field_ext_us % val = u_s( 1, : )
    call insert( alg_ext, field_ext_us, "SolidVelocity1" )

    call allocate( field_ext_vs, positions_v%mesh, "SolidVelocity2" )
    field_ext_vs % val = u_s( 2, : )
    call insert( alg_ext, field_ext_vs, "SolidVelocity2" )

    call allocate( field_ext_solid, positions_v%mesh, "SolidConcentration" )
    field_ext_solid % val = 1.0
    call insert( alg_ext, field_ext_solid, "SolidConcentration" )

    ! deal with SolidConcentration
    solid => extract_scalar_field( packed_state, "SolidConcentration" )
    old_solid => extract_scalar_field( packed_state, "OldSolidConcentration" )

    call set( old_solid, solid )
    call zero( solid )

    ! interpolate
    call interpolation_galerkin_femdem( alg_ext, alg_fl, field = solid )

    ! bound solid concentration
    call bound_volume_fraction( solid % val )

    u_mesh => extract_mesh( packed_state, "VelocityMesh" )
    call allocate( f2, u_mesh, "dummy" )

    ! to be used for the supplementary equation (Eq. 124)
    ! that is u_hat = delta_u + u_f
    delta_u => extract_vector_field( packed_state, "delta_U" )
    do idim = 1, ndim
       f => extract_scalar_field( alg_fl, "deltaVelocity" // int2str( idim ) )

       !call project_field( f, f2, fl_positions )
       call linear2quadratic_field(f, f2)

       delta_u % val( idim, : ) = f2 % val
    end do

    solid_u => extract_vector_field( packed_state, "solid_U" )
    do idim = 1, ndim
       f => extract_scalar_field( alg_fl, "SolidVelocity" // int2str( idim ) )
       !call project_field( f, f2, fl_positions )
       call linear2quadratic_field(f, f2)


       solid_u % val( idim, : ) = f2 % val
    end do


    ! deallocate
    call deallocate( field_fl_du )
    call deallocate( field_fl_dv )
    call deallocate( field_fl_us )
    call deallocate( field_fl_vs )
    call deallocate( field_fl_solid )

    call deallocate( field_ext_du )
    call deallocate( field_ext_dv )
    call deallocate( field_ext_us )
    call deallocate( field_ext_vs )
    call deallocate( field_ext_solid )

    call deallocate( alg_fl )
    call deallocate( alg_ext )

    call deallocate( f2 )

    return
  end subroutine interpolate_fields_in_v

  !----------------------------------------------------------------------------------------------------------

  subroutine interpolate_fields_in_r( packed_state, fin, ain )

  
   
    implicit none

    type( state_type ), intent( in ) :: packed_state
    real, dimension( :, : ), intent( in ) :: fin
    real, dimension( :, :, : ), intent( in ) :: ain

    !Local variables
    type( mesh_type ), pointer :: fl_mesh, p_mesh
    type( scalar_field ), pointer :: f
    type( scalar_field ) :: field_fl_f1, field_fl_f2, &
         &                  field_fl_a11, field_fl_a12, field_fl_a22, &
         &                  field_ext_f1, field_ext_f2, &
         &                  field_ext_a11, field_ext_a12, field_ext_a22, f2, dummy
    type( vector_field ), pointer :: fl_positions, f_x
    type( tensor_field ), pointer :: a_xx
    type( state_type ) :: alg_ext, alg_fl
    integer :: stat, idim
    character( len = OPTION_PATH_LEN ) :: path = "/tmp/galerkin_projection/continuous"

    call set_solver_options( path, &
         ksptype = "gmres", &
         pctype = "hypre", &
         rtol = 1.0e-10, &
         atol = 0.0, &
         max_its = 10000 )
    call add_option( &
         trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", stat)
    call set_option( &
         trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", "boomeramg")

    path = "/tmp"

    fl_mesh => extract_mesh( packed_state, "CoordinateMesh" )
    fl_positions => extract_vector_field( packed_state, "Coordinate" )

    ! fluidity state
    call insert( alg_fl, fl_mesh, "Mesh" )
    call insert( alg_fl, fl_positions, "Coordinate" )

    call allocate( field_fl_f1, fl_positions%mesh, "f1" )
    call zero( field_fl_f1 )
    field_fl_f1 % option_path = path
    call insert( alg_fl, field_fl_f1, "f1" )

    call allocate( field_fl_f2, fl_positions%mesh, "f2" )
    call zero( field_fl_f2 )
    field_fl_f2 % option_path = path
    call insert( alg_fl, field_fl_f2, "f2" )

    call allocate( field_fl_a11, fl_positions%mesh, "a11" )
    call zero( field_fl_a11 )
    field_fl_a11 % option_path = path
    call insert( alg_fl, field_fl_a11, "a11" )

    call allocate( field_fl_a12, fl_positions%mesh, "a12" )
    call zero( field_fl_a12 )
    field_fl_a12 % option_path = path
    call insert( alg_fl, field_fl_a12, "a12" )

    call allocate( field_fl_a22, fl_positions%mesh, "a22" )
    call zero( field_fl_a22 )
    field_fl_a22 % option_path = path
    call insert( alg_fl, field_fl_a22, "a22" )

    ! solid volume state
    call insert( alg_ext, positions_r%mesh, "Mesh" )
    call insert( alg_ext, positions_r, "Coordinate" )

    call allocate( field_ext_f1, positions_r%mesh, "f1" )
    field_ext_f1 % val = fin( 1, : )
    call insert( alg_ext, field_ext_f1, "f1" )

    call allocate( field_ext_f2, positions_r%mesh, "f2" )
    field_ext_f2 % val = fin( 2, : )
    call insert( alg_ext, field_ext_f2, "f2" )

    call allocate( field_ext_a11, positions_r%mesh, "a11" )
    field_ext_a11 % val = ain( 1, 1, : )
    call insert( alg_ext, field_ext_a11, "a11" )

    call allocate( field_ext_a12, positions_r%mesh, "a12" )
    field_ext_a12 % val = ain( 1, 2, : )
    call insert( alg_ext, field_ext_a12, "a12" )

    call allocate( field_ext_a22, positions_r%mesh, "a22" )
    field_ext_a22 % val = ain( 2, 2, : )
    call insert( alg_ext, field_ext_a22, "a22" )

    call allocate( dummy, fl_positions%mesh, "dummy" )
    call zero( dummy )

    ! interpolate
    call interpolation_galerkin_femdem( alg_ext, alg_fl, field = dummy )


    p_mesh => extract_mesh( packed_state, "PressureMesh" ) 
    call allocate( f2, p_mesh, "dummy" )

    f_x => extract_vector_field( packed_state, "f_x" )
    do idim = 1, ndim
       f => extract_scalar_field( alg_fl, "f" // int2str( idim ) )
       call linear2quadratic_field(f, f2)
       f_x % val( idim, : ) = f2 % val
    end do

    a_xx => extract_tensor_field( packed_state, "a_xx" )
    do idim = 1, ndim
       f => extract_scalar_field( alg_fl, "a1" // int2str( idim ) )
       call linear2quadratic_field(f, f2)
       a_xx % val( 1, idim, : ) = f2 % val
    end do
    f => extract_scalar_field( alg_fl, "a22" )
    call linear2quadratic_field(f, f2)
    a_xx % val( 2, 2, : ) = f2 % val
    a_xx % val( 2, 1, :) = a_xx % val(1, 2, :)



    ! deallocate
    call deallocate( field_fl_f1 )
    call deallocate( field_fl_f2 )
    call deallocate( field_fl_a11 )
    call deallocate( field_fl_a12 )
    call deallocate( field_fl_a22 )
    
    call deallocate( field_ext_f1 )
    call deallocate( field_ext_f2 )
    call deallocate( field_ext_a11 )
    call deallocate( field_ext_a12 )
    call deallocate( field_ext_a22 )

    call deallocate( alg_fl )
    call deallocate( alg_ext )

    call deallocate( f2 )
    call deallocate( dummy )

    return
  end subroutine interpolate_fields_in_r

    !----------------------------------------------------------------------------------------------------------

    subroutine coarsen_mesh_2d( f, c )

        implicit none

        type( vector_field ), intent( in ) :: f
        type( vector_field ), intent( out ) :: c
        !Local variables
        integer :: snloc, nloc, nele_s, nele, ele, nnodes_s, &
        ele2, siloc, sjloc, iloc, jloc, st, mv, &
        edge, nedge, bcs, nbcs, count, ne, i, j, &
        quad_degree
        integer, dimension( : ), allocatable :: ndglno, sndglno, boundary_ids, &
        ndglno_s, tmp, ele_nodes, &
        ele2_nodes, bc, n, nodes_s
        real :: area
        real, dimension( : ), allocatable :: X_s, Y_s
        logical :: on_the_wall, delete_ele2
        type( quadrature_type ) :: quad
        type( element_type ) :: shape
        type( mesh_type ) :: mesh

        logical, parameter :: coarsen_mesh = .false.

        nedge = 3 ! triangles
        nloc  = 3 ! 3 nodes per element
        snloc = 2 ! 2 nodes per edge

        nbcs = surface_element_count( f )

        if ( coarsen_mesh ) then

            nele = element_count( f )

            allocate( ndglno( nloc * nele ), sndglno( snloc * nbcs ) )
            allocate( ele_nodes( nloc ), ele2_nodes( nloc ) )
            allocate( n( snloc ), bc(snloc ) )

            ndglno = f % mesh % ndglno
            call getsndgln( f % mesh, sndglno)

            nele_s = nele ! number of elements in the simplified mesh

            do ele = 1, nele

                ele_nodes = ndglno( ( ele - 1 ) * nloc + 1 : ele * nloc )

                ! check if this element has been deleted
                if ( all( ele_nodes < 0 ) ) cycle

                do edge = 1, nedge

                    n( 1 ) = ele_nodes( edge )
                    if ( edge < nedge ) then
                        n( 2 ) = ele_nodes( edge+1 )
                    else
                        n( 2 ) = ele_nodes( 1 )
                    end if

                    on_the_wall = .false.
                    do bcs = 1, nbcs
                        bc = sndglno( (bcs-1)*snloc+1 : bcs*snloc )

                        ! on_the_wall = .true. if all nodes on the boundary
                        if ( all( n == bc ) ) then
                            on_the_wall = .true.
                            exit
                        end if
                        ! if one node is on the wall
                        ! figure out which one it is
                        st = -666
                        if ( any( n == bc ) .and. .not.on_the_wall ) then
                            do siloc = 1, snloc
                                do sjloc = 1, snloc
                                    if ( n( siloc ) == bc( sjloc ) ) st = n( siloc )
                                    if ( n( siloc ) /= bc( sjloc ) ) mv = n( siloc )
                                end do
                            end do
                        end if
                    end do

                    if ( any( sndglno == mv ) ) on_the_wall = .true.

                    if ( .not.on_the_wall ) then

                        if ( st < 0 ) then
                            st = n( 1 )
                            mv = n( 2 )
                        end if

                        do ele2 = 1, nele

                            ele2_nodes = ndglno( ( ele2 - 1 ) * nloc + 1 : ele2 * nloc )

                            ! check if this element has been deleted
                            if ( all( ele2_nodes < 0 ) ) cycle

                            do iloc = 1, nloc
                                if ( ele2_nodes( iloc ) == mv ) &
                                ndglno( ( ele2 - 1 ) * nloc + iloc ) = st
                            end do

                            ! update local memory
                            ele2_nodes = ndglno( ( ele2 - 1 ) * nloc + 1 : ele2 * nloc )

                            ! figure out if we need to delete ele2
                            delete_ele2 = .false.
                            if ( ele2 == ele ) then
                                delete_ele2 = .true.
                            else
                                do iloc = 1, nloc
                                    do jloc = iloc+1, nloc
                                        if ( ele2_nodes( iloc ) == ele2_nodes( jloc ) ) &
                                        delete_ele2 = .true.
                                    end do
                                end do
                            end if

                            if ( delete_ele2 ) then
                                ndglno( ( ele2 - 1 ) * nloc + 1 : ele2 * nloc ) = -666
                                nele_s = nele_s - 1
                            end if

                        end do ! ele2, nele
                    end if ! .not.on_the_wall
                end do ! edge, nedge
            end do ! ele, nele

            ! create simplified mesh

            allocate( ndglno_s( nele_s * nloc ), tmp( nele_s * nloc ) )

            count = 1
            do ele = 1, nele
                ele_nodes = ndglno( ( ele - 1 ) * nloc + 1 : ele * nloc )
                if ( all( ele_nodes < 0) ) cycle
                ndglno_s( ( count - 1 ) * nloc + 1 : count * nloc ) = ele_nodes
                count = count +1
            end do

            ! make sure elements aren't inside out
            do ele = 1, nele_s
                ele_nodes = ndglno_s( ( ele - 1 ) * nloc + 1 : ele * nloc )
                area = triangle_area( &
                f % val( 1, ele_nodes( 1 ) ), f % val( 2, ele_nodes( 1 ) ), &
                f % val( 1, ele_nodes( 2 ) ), f % val( 2, ele_nodes( 2 ) ), &
                f % val( 1, ele_nodes( 3 ) ), f % val( 2, ele_nodes( 3 ) ) )
                if ( area < 0.0 ) then
                    ! swap 2nd and 3rd nodes
                    ndglno_s( ( ele - 1 ) * nloc + 2 : ele * nloc ) = (/ ele_nodes( 3 ), ele_nodes( 2 ) /)
                end if
            end do

            ! make sure we've recovered all the elements
            assert( nele_s == count-1 )

            ! figure out which nodes are still in the mesh and count them
            tmp = ndglno_s
            call delete_duplicates( tmp, nnodes_s )
            allocate( nodes_s( nnodes_s ) ) ; nodes_s = tmp( 1 : nnodes_s )
            ! sort the node numbers
            call quicksort(nodes_s, 1)

            ! re-numbering

            allocate( x_s( nnodes_s ), y_s( nnodes_s ) )

            ! deal with nodes
            do i = 1, nnodes_s
                j = nodes_s( i )
                x_s( i ) = f % val( 1, j )
                y_s( i ) = f % val( 2, j )
            end do

            ! deal with ndgln
            do ele = 1, nele_s
                ele_nodes = ndglno_s( ( ele - 1 ) * nloc + 1 : ele * nloc )
                do iloc = 1, nloc
                    ! figure out new node number
                    do i = 1, nnodes_s
                        if ( nodes_s( i ) == ele_nodes( iloc ) ) then
                            ne = i
                            exit
                        end if
                    end do
                    ! amend ndglno
                    ndglno_s( ( ele - 1 ) * nloc + iloc ) = ne
                end do
            end do

            ! deal with sndgln
            do bcs = 1, nbcs
                bc = sndglno( ( bcs - 1 ) * snloc + 1 : bcs * snloc )
                do siloc = 1, snloc
                    ! figure out new node number
                    do i = 1, nnodes_s
                        if ( nodes_s( i ) == bc( siloc ) ) then
                            ne = i
                            exit
                        end if
                    end do
                    ! amend sndglno
                    sndglno( ( bcs - 1 ) * snloc + siloc ) = ne
                end do
            end do

            call get_option( "/geometry/quadrature/degree", quad_degree )

            quad = make_quadrature( nloc, ndim, degree = quad_degree )
            shape = make_element_shape( nloc, ndim, 1, quad )

            ! create the coarse mesh
            call allocate( mesh, nnodes_s, nele_s, shape, name="CoordinateMesh" )
            call allocate( c, ndim, mesh, name="Coordinate" )

            c % val( 1, : ) = x_s
            c % val( 2, : ) = y_s

            c % mesh % ndglno = ndglno_s

            allocate( boundary_ids( nbcs ) ) ; boundary_ids = 666
            call add_faces( c % mesh, &
            sndgln = sndglno, &
            boundary_ids = boundary_ids )

            c % dim = ndim

            deallocate( boundary_ids )

            call deallocate( mesh )
            call deallocate_element( shape )
            call deallocate( quad )

            deallocate( ndglno_s, tmp )
            deallocate( x_s, y_s )

            deallocate( ndglno, sndglno )
            deallocate( ele_nodes, ele2_nodes, n, bc )

        else

            call get_option( "/geometry/quadrature/degree", quad_degree )

            quad = make_quadrature( nloc, ndim, degree = quad_degree )
            shape = make_element_shape( nloc, ndim, 1, quad )

            call allocate( mesh, node_count( f ), element_count( f ), shape, name="CoordinateMesh" )
            call allocate( c, ndim, mesh, name="Coordinate" )
            c % val = f % val

            c % mesh % ndglno = f % mesh % ndglno

            allocate( sndglno( snloc*nbcs ) ) ; sndglno = 0
            call getsndgln( f % mesh, sndglno)

            allocate( boundary_ids( nbcs ) ) ; boundary_ids = 666
            call add_faces( c % mesh, &
            sndgln = sndglno, &
            boundary_ids = boundary_ids )

            c % dim = ndim

            deallocate( sndglno, boundary_ids )

            call deallocate( mesh )
            call deallocate_element( shape )
            call deallocate( quad )

        end if

        return
    end subroutine coarsen_mesh_2d

    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------

     subroutine delete_duplicates( a, count )

        implicit none

        integer, dimension( : ), intent( inout ) :: a
        integer, intent( out ) :: count
        !Local variables
        integer :: n, i, j
        integer, dimension(:), allocatable :: tmp

        n = size( a )

        allocate( tmp( n ) ) ; tmp = a

        count = 1
        tmp( 1 ) = tmp( 1 )

        outer: do i = 2, n
            do j = 1, count
                ! Found a match so start looking again
                if ( tmp( j ) == a( i ) ) cycle outer
            end do
            ! No match found so add it to the output
            count = count + 1
            tmp( count ) = a( i )

        end do outer

        a = tmp

        deallocate( tmp )

        return
    end subroutine delete_duplicates

    subroutine bound_volume_fraction ( v, v_min, v_max )

        implicit none

        real, dimension( : ), intent( inout ) :: v
        real, intent( in ), optional :: v_min, v_max
        real :: vmin, vmax

        vmin = 0.0 ; vmax = 1.0
        if ( present( v_min ) ) vmin = v_min
        if ( present( v_max ) ) vmax = v_max

        v = max( vmin, min( vmax, v ) )

        return
    end subroutine bound_volume_fraction

    subroutine deallocate_femdem

        implicit none

        call deallocate( positions_r )
        call deallocate( positions_v )
        call deallocate( positions_vc )
        call deallocate( permeability_r )

        return
    end subroutine deallocate_femdem

    real function triangle_area( x1, y1, x2, y2, x3, y3 )

        implicit none

        real :: x1, y1, x2, y2, x3, y3

        triangle_area = 0.5 * ( ( x2 * y3 - y2 * x3 ) - x1 * ( y3 - y2 ) + y1 * ( x3 - x2 ) )

        return
    end function triangle_area


    subroutine linear2quadratic_field( field_in, field_out )
        implicit none
        type( scalar_field ), intent( in ) :: field_in
        type( scalar_field ), intent( inout ) :: field_out

        integer, dimension( : ), pointer :: p1_ndglno, p2_ndglno,  p1_nods
        integer :: n1, n2,  totele, p1_nloc, p2_nloc, ele, p2_iloc, p2_nod
        real, dimension( : ), allocatable :: field_p1_loc, field_p2_loc


        n1 = field_in%mesh%shape%degree
        n2 = field_out%mesh%shape%degree

        if ( n1==1 .and. n2==2) then

            p1_ndglno => get_ndglno( field_in%mesh )
            p2_ndglno => get_ndglno( field_out%mesh )

            totele = field_in%mesh%elements
            p1_nloc = field_in%mesh%shape%loc
            p2_nloc = field_out%mesh%shape%loc

            allocate( field_p1_loc( p1_nloc ) )
            allocate( field_p2_loc( p2_nloc ) )

            do ele = 1, totele

                p1_nods => p1_ndglno( ( ele - 1 ) * p1_nloc + 1 : ele * p1_nloc )
                field_p1_loc =  field_in % val( p1_nods )


                field_p2_loc( 1 ) =field_p1_loc( 1 )
                field_p2_loc( 3 ) =field_p1_loc( 2 )
                field_p2_loc( 6 ) =field_p1_loc( 3 )

                field_p2_loc( 2 ) = 0.5 * ( field_p1_loc( 1 ) + field_p1_loc( 2 ) )
                field_p2_loc( 4 ) = 0.5 * ( field_p1_loc( 1 ) + field_p1_loc( 3 ) )
                field_p2_loc( 5 ) = 0.5 * ( field_p1_loc( 2 ) + field_p1_loc( 3 ) )

                if ( p2_nloc == 10 ) then
                    field_p2_loc( 7 ) = 0.5 * ( field_p1_loc(  1 ) + field_p1_loc(  4 ) )
                    field_p2_loc( 8 ) = 0.5 * ( field_p1_loc(  2 ) + field_p1_loc(  4 ) )
                    field_p2_loc( 9 ) = 0.5 * ( field_p1_loc(  3 ) + field_p1_loc(  4 ) )
                    field_p2_loc( 10 ) = field_p1_loc(  4 )
                end if

                do p2_iloc = 1, p2_nloc
                    p2_nod = p2_ndglno( ( ele - 1 ) * p2_nloc + p2_iloc )
                    field_out % val( p2_nod ) = field_p2_loc( p2_iloc )
                end do

            end do

            deallocate( field_p1_loc, field_p2_loc)
        else if ( n1==1 .and. n2==1) then





            p1_ndglno => get_ndglno( field_in%mesh )
            p2_ndglno => get_ndglno( field_out%mesh )

            totele = field_in%mesh%elements
            p1_nloc = field_in%mesh%shape%loc
            p2_nloc = field_out%mesh%shape%loc

            allocate( field_p1_loc( p1_nloc ) )
            allocate( field_p2_loc( p2_nloc ) )

            do ele = 1, totele

                p1_nods => p1_ndglno( ( ele - 1 ) * p1_nloc + 1 : ele * p1_nloc )
                field_p1_loc =  field_in % val( p1_nods )


                field_p2_loc( 1 ) =field_p1_loc( 1 )
                field_p2_loc( 2 ) =field_p1_loc( 2 )
                field_p2_loc( 3 ) =field_p1_loc( 3 )

                

                if ( p2_nloc == 6 ) then
                    field_p2_loc( 4 ) =field_p1_loc( 4 )
                    field_p2_loc( 5 ) =field_p1_loc( 5 )
                    field_p2_loc( 6 ) =field_p1_loc( 6 )
                end if

                do p2_iloc = 1, p2_nloc
                    p2_nod = p2_ndglno( ( ele - 1 ) * p2_nloc + p2_iloc )
                    field_out % val( p2_nod ) = field_p2_loc( p2_iloc )
                end do

            end do

            deallocate( field_p1_loc, field_p2_loc)


        else
            
            ! we have not called this routine with a p2 field
            stop 28289

        end if

        return
    end subroutine linear2quadratic_field

    !------------------------- interpolate just pressure -----------------------------------------------
    subroutine interpolate_fields_out_r_p( packed_state, nphase, p_r)

        implicit none

        type( state_type ), intent( in ) :: packed_state
        integer, intent( in ) :: nphase
        real, dimension( : ), intent( inout ) :: p_r

        !Local variables
        type( mesh_type ), pointer :: fl_mesh, p0_fl_mesh
        type( scalar_field ) :: field_fl_p,  &
        field_ext_p,  &
        rvf
        type( tensor_field), pointer :: pressure
        type( vector_field ), pointer :: fl_positions
        type( state_type ) :: alg_ext, alg_fl
        integer, dimension( : ), pointer :: fl_ele_nodes, cv_ndgln
        integer :: ele, totele, u_nloc, cv_nloc, u_nonods, i, &
		stat
        character( len = OPTION_PATH_LEN ) :: path = "/tmp/galerkin_projection/continuous"

        p_r = 0.0 ;

        cv_ndgln => get_ndglno( extract_mesh( packed_state, "PressureMesh" ) )
        pressure => extract_tensor_field( packed_state, "PackedFEPressure" )
        cv_nloc = ele_loc( pressure, 1 )

        fl_mesh => extract_mesh( packed_state, "CoordinateMesh" )
        fl_positions => extract_vector_field( packed_state, "Coordinate" )
        totele = ele_count( fl_mesh )

        call set_solver_options( path, &
        ksptype = "gmres", &
        pctype = "hypre", &
        rtol = 1.0e-10, &
        atol = 0.0, &
        max_its = 10000 )
        call add_option( &
        trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", stat)
        call set_option( &
        trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", "boomeramg")

        path = "/tmp"

        call allocate( field_fl_p, fl_mesh, "Pressure" )
        call zero( field_fl_p )

        ! deal with pressure and viscosity
        if ( cv_nloc == 6 ) then
            ! linearise pressure for p2
            do ele = 1, totele
                fl_ele_nodes => ele_nodes( fl_mesh, ele )
                field_fl_p % val( fl_ele_nodes( 1 ) ) = pressure % val(1,1, cv_ndgln( ( ele - 1 ) * cv_nloc + 1 ) )
                field_fl_p % val( fl_ele_nodes( 2 ) ) = pressure % val(1,1, cv_ndgln( ( ele - 1 ) * cv_nloc + 3 ) )
                field_fl_p % val( fl_ele_nodes( 3 ) ) = pressure % val(1,1, cv_ndgln( ( ele - 1 ) * cv_nloc + 6 ) )
            end do
        else
            ! just copy memory for p1
            field_fl_p % val = pressure % val(1,1,:)
        end if



        ! fluidity state
        call insert( alg_fl, fl_mesh, "Mesh" )
        call insert( alg_fl, fl_positions, "Coordinate" )

        call insert( alg_fl, field_fl_p, "Pressure" )

        ! ring state
        call insert( alg_ext, positions_r%mesh, "Mesh" )
        call insert( alg_ext, positions_r, "Coordinate" )

        call allocate( field_ext_p, positions_r%mesh, "Pressure" )
        call zero( field_ext_p )
        field_ext_p % option_path = path
        call insert( alg_ext, field_ext_p, "Pressure" )

        ewrite(3,*) "...interpolating"
        ! interpolate
        call interpolation_galerkin_femdem( alg_fl, alg_ext, femdem_out = .true. ) !!-ao problem here...

        ! copy memory
        p_r = field_ext_p % val

        ! deallocate
        call deallocate( field_fl_p )
        call deallocate( field_ext_p )
        call deallocate( alg_fl )
        call deallocate( alg_ext )


        return
    end subroutine interpolate_fields_out_r_p
!----------------------------------------------------------------------------------------------------------

!------------------------- pore fluid pressure -------------------------------------------------------!
subroutine interpolate_fields_out_v_pf( packed_state, nphase, p_v)

    implicit none

    type( state_type ), intent( in ) :: packed_state
    integer, intent( in ) :: nphase
    real, dimension( :), intent( inout ) :: p_v

    !Local variables
    type( mesh_type ), pointer :: fl_mesh, u_mesh
    type( scalar_field ) :: field_fl_p,  &
         &                  field_ext_p

    type( vector_field ), pointer :: fl_positions
    type( tensor_field ), pointer :: pressure

    type( state_type ) :: alg_ext, alg_fl

    integer, dimension(:), pointer:: cv_ndgln, fl_ele_nodes
    integer :: ele, totele, cv_nloc, &
    stat
    character( len = OPTION_PATH_LEN ) :: path = "/tmp/galerkin_projection/continuous"

    p_v = 0.0


    cv_ndgln => get_ndglno( extract_mesh( packed_state, "PressureMesh" ) )
    pressure => extract_tensor_field( packed_state, "PackedFEPressure" )
    cv_nloc = ele_loc( pressure, 1 )

    fl_mesh => extract_mesh( packed_state, "CoordinateMesh" )
    fl_positions => extract_vector_field( packed_state, "Coordinate" )
    totele = ele_count( fl_mesh )

    call set_solver_options( path, &
         ksptype = "gmres", &
         pctype = "hypre", &
         rtol = 1.0e-10, &
         atol = 0.0, &
         max_its = 10000 )
    call add_option( &
         trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", stat)
    call set_option( &
         trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", "boomeramg")

    path = "/tmp"


    call allocate( field_fl_p, fl_mesh, "Pressure" )
    call zero( field_fl_p )


       if ( cv_nloc == 6 ) then
            ! linearise pressure for p2
            do ele = 1, totele
                fl_ele_nodes => ele_nodes( fl_mesh, ele )
                field_fl_p % val( fl_ele_nodes( 1 ) ) = pressure % val( 1,1, cv_ndgln( ( ele - 1 ) * cv_nloc + 1 ) )
                field_fl_p % val( fl_ele_nodes( 2 ) ) = pressure % val( 1,1, cv_ndgln( ( ele - 1 ) * cv_nloc + 3 ) )
                field_fl_p % val( fl_ele_nodes( 3 ) ) = pressure % val( 1,1, cv_ndgln( ( ele - 1 ) * cv_nloc + 6 ) )
            end do
        else
            ! just copy memory for p1
            field_fl_p % val = pressure % val (1,1,:)
        end if

    ! fluidity state
    call insert( alg_fl, fl_mesh, "Mesh" )
    call insert( alg_fl, fl_positions, "Coordinate" )

    call insert( alg_fl, field_fl_p, "Pressure" )

    ! ring state
    call insert( alg_ext, positions_v%mesh, "Mesh" )
    call insert( alg_ext, positions_v, "Coordinate" )

    call allocate( field_ext_p, positions_v%mesh, "Pressure" )
    call zero( field_ext_p )

    field_ext_p % option_path = path
    call insert( alg_ext, field_ext_p, "Pressure" )

    ewrite(3,*) "...interpolating pore fluid pressure"

    ! interpolate
    call interpolation_galerkin_femdem( alg_fl, alg_ext, femdem_out = .true. )

    ! copy memory
    p_v= field_ext_p % val

    print *, "FLUIDS: entered into pore_pressure sub-------------------------------------------------------", size(p_v), size(field_ext_p%val), totele


    ! deallocate
    call deallocate( field_fl_p )
    call deallocate( field_ext_p )

    call deallocate( alg_fl )
    call deallocate( alg_ext )

    return
  end subroutine interpolate_fields_out_v_pf
! ----------------------------------------------- end pore fluid pressure ----------------!



end module multiphase_fractures
