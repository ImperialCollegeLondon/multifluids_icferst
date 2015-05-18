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


!  module multi_interpolation
!
!    use fldebug
!    use state_module
!    use fields
!    use field_options
!    use spud
!    use populate_state_module
!    use diagnostic_variables
!    use diagnostic_fields
!    use diagnostic_fields_wrapper
!    use global_parameters, only: option_path_len, is_porous_media, dumping_in_sat, is_multifracture
!    use diagnostic_fields_wrapper_new
!    use element_numbering
!    use shape_functions
!    use fefields
!    use boundary_conditions
!    use futils, only: int2str
!    use boundary_conditions_from_options
!    use parallel_tools, only : allmax, allmin, isparallel
!    use memory_diagnostics
!    use initialise_fields_module, only: initialise_field_over_regions
!    use halos
!    use fields_allocates
!
!
!    use solvers
!    use conservative_interpolation_module
!
!    use shape_functions
!    use shape_functions_Linear_Quadratic
!    use solvers_module
!    use matrix_operations
!    use Copy_Outof_State
!    use sparse_tools_petsc
!    use shape_functions_prototype
!
!    type(state_type) :: state_old, state_new
!    type(scalar_field) :: ph_sol_old, ph_sol_new
!
!
!  contains
!
!      subroutine M2MInterpolation(state, packed_state, StorageIndexes, cv_ele_type, nphase, flag )
!
!          implicit none
!
!          type( state_type ), dimension( : ), intent( inout ) :: state
!          type( state_type ), intent( inout ) :: packed_state
!          integer, intent( in ) :: cv_ele_type, nphase
!          integer, dimension( : ), intent( inout ) :: StorageIndexes
!          integer, intent(in) :: flag
!
!          ! local variables...
!          type ( tensor_field ), pointer :: ufield
!          integer :: ndim, ph_ngi, ph_ngi_short, ph_nloc, ph_snloc, &
!          u_nloc, u_snloc, sphngi, sbphngi, nface, stat, &
!          totele, x_nonods, ele, iloc, jloc, x_nloc, &
!          cv_snloc, ph_ele_type, iloop, u_nonods, cv_nonods, &
!          cv_iloc, cv_jloc, cv_inod, idim, iphase, u_inod, u_iloc, cv_nloc, &
!          ph_iloc, ph_inod, ph_nonods, ph_jloc, ph_jnod, tmp_cv_nloc, other_nloc
!          real, dimension( : ), pointer :: phweight, phweight_short, sphfeweigh, sbphfeweigh, &
!          sele_overlap_scale
!          real, dimension( :, : ), pointer :: phn, phn_short, phfen, phfen_short, ufen, &
!          sphfen, sphfenslx, sphfensly, sufen, sufenslx, sufensly, &
!          sbphn, sbphfen, sbphfenslx, sbphfensly, sbufen, sbufenslx, sbufensly
!          real, dimension( :, :, : ), pointer :: phfenlx_all, phfenlx_short_all, ufenlx_all, &
!          sphfenlx_all, sufenlx_all, sbphfenlx_all, sbufenlx_all
!          logical, dimension( :, : ), allocatable :: u_on_face, ufem_on_face, &
!          ph_on_face, phfem_on_face
!          integer, pointer :: ncolgpts
!          integer, dimension( : ), pointer :: findgpts, colgpts, x_ndgln, cv_ndgln, ph_ndgln, u_ndgln
!          integer, dimension( :, : ), pointer :: ph_neiloc, ph_sloclist, u_sloclist
!          logical :: quad_over_whole_ele, d1, d3, dcyl
!          type( vector_field ), pointer :: x
!          type( mesh_type ), pointer :: phmesh
!
!          real, dimension( : ), pointer :: detwei, ra
!          real, pointer :: volume
!          real, dimension(:,:,:), pointer :: phfenx_all, ufenx_all
!
!          real, dimension( :, :, : ), allocatable :: u_ph_source_vel, u_ph_source_cv
!          real, dimension( :, : ), allocatable :: alpha_cv, coef_alpha_cv
!
!          real, dimension( :, :, : ), allocatable :: u_ph_source_ph, dx_ph_gi
!          real, dimension( :, : ), allocatable :: alpha_ph, coef_alpha_ph, ph
!
!          real, dimension( :, :, : ), allocatable :: u_s_gi, dx_alpha_gi
!          real, dimension( :, : ), allocatable :: coef_alpha_gi, den_gi, inv_den_gi
!
!          real, dimension( : ), pointer :: tmp_cv_weight
!          real, dimension( :, : ), pointer :: tmp_cvfen
!          real, dimension( :, :, : ), pointer :: tmp_cvfenlx_all
!          real, dimension( :, :, : ), pointer :: tmp_cvfenx_all
!
!          real, dimension( :, : ), pointer :: other_fen
!          real, dimension( :, :, : ), pointer :: other_fenlx_all
!          real, dimension( :, :, : ), pointer :: other_fenx_all
!
!          real :: nxnx, nm
!
!          type( scalar_field ) :: rhs, ph_sol
!          type( petsc_csr_matrix ) :: matrix
!          type( csr_sparsity ), pointer :: sparsity
!
!          character( len = OPTION_PATH_LEN ) :: path = "/tmp"
!
!          type( tensor_field ), pointer :: rho
!          type( scalar_field ), pointer :: printf
!          type( vector_field ), pointer :: printu, x_p2
!
!
!          logical :: on_boundary
!          integer :: ph_jnod2, ierr, count, count2, i,j
!          integer, dimension(:), pointer :: findph, colph
!
!          type(petsc_csr_matrix) :: mat
!          type(scalar_field) :: fempsi_rhs
!          type( mesh_type ), pointer :: mesh_old, mesh_new
!          type( vector_field ) :: positions_old, positions_new
!          integer :: cv_nodi, cv_nodj, cv_ngi, cv_gi
!          real :: MN, MM
!          type(tensor_field), pointer :: PhaseV_old
!
!          real, dimension( :, : ), allocatable :: cvn
!          real, dimension( :, : ), allocatable :: N
!
!          ! CHECK THESE ALLOCATIONS
!
!          allocate(cvn(cv_nonods, cv_ngi))
!          allocate(N(cv_nonods, cv_ngi))
!
!          call get_option( '/geometry/dimension', ndim )
!
!          ufield => extract_tensor_field( packed_state, "PackedVelocity" )
!
!          u_nloc = ele_loc( ufield, 1 ) ; u_snloc = face_loc( ufield, 1 )
!
!          u_nonods = node_count( ufield )
!
!          quad_over_whole_ele = .true.
!
!          ! ph elements
!          if ( ndim == 2 ) then
!              ph_ele_type = 4
!              ph_nloc = 6 ; ph_snloc = 3
!          else if ( ndim == 3 ) then
!              ph_ele_type = 8
!              ph_nloc = 10 ; ph_snloc = 6
!          else
!              stop 567
!          end if
!
!          call retrieve_ngi( ndim, ph_ele_type, ph_nloc, u_nloc, &
!          ph_ngi, ph_ngi_short, sphngi, sbphngi, nface, quad_over_whole_ele )
!
!          allocate( ph_on_face( ph_nloc, sphngi ), phfem_on_face( ph_nloc, sphngi ) )
!          allocate( u_on_face( u_nloc, sphngi ), ufem_on_face( u_nloc, sphngi ) )
!
!          call cv_fem_shape_funs_plus_storage( &
!                               ! volume shape functions...
!          ndim, ph_ele_type,  &
!          ph_ngi, ph_ngi_short, ph_nloc, u_nloc, phn, phn_short, &
!          phweight, phfen, phfenlx_all, &
!          phweight_short, phfen_short, phfenlx_short_all, &
!          ufen, ufenlx_all, &
!                               ! surface of each ph shape functions...
!          sphngi, ph_neiloc, ph_on_face, phfem_on_face, &
!          sphfen, sphfenslx, sphfensly, sphfeweigh, &
!          sphfenlx_all,  &
!          sufen, sufenslx, sufensly, &
!          sufenlx_all, &
!                               ! surface element shape funcs...
!          u_on_face, ufem_on_face, nface, &
!          sbphngi, sbphn, sbphfen, sbphfenslx, sbphfensly, sbphfeweigh, sbphfenlx_all, &
!          sbufen, sbufenslx, sbufensly, sbufenlx_all, &
!          ph_sloclist, u_sloclist, ph_snloc, u_snloc, &
!                               ! define the gauss points that lie on the surface of the ph...
!          findgpts, colgpts, ncolgpts, &
!          sele_overlap_scale, quad_over_whole_ele, &
!          state, "ph_1" , storageindexes( 36 ) )
!
!          totele = ele_count( ufield )
!          x_ndgln => get_ndglno( extract_mesh( state( 1 ), "PressureMesh_Continuous" ) )
!          cv_ndgln => get_ndglno( extract_mesh( state( 1 ), "PressureMesh" ) )
!          x_nonods = node_count( extract_mesh( state( 1 ), "PressureMesh_Continuous" ) )
!          x => extract_vector_field( packed_state, "PressureCoordinate" )
!          u_ndgln => get_ndglno( extract_mesh( state( 1 ), "VelocityMesh" ) )
!          x_nloc = ele_loc( x, 1 )
!          cv_nloc = x_nloc
!          cv_nonods = node_count( extract_mesh( state( 1 ), "PressureMesh" ) )
!          phmesh => extract_mesh( state( 1 ), "ph" )
!          ph_ndgln => get_ndglno( phmesh )
!          ph_nonods = node_count( phmesh )
!          d1 = ( ndim == 1 ) ; d3 = ( ndim == 3 ) ; dcyl = .false.
!
!            ! ALLOCATIONS - Check what actually needs to be allocated
!
!          if (flag == 0) then
!
!              mesh_old => extract_mesh( packed_state, "CoordinateMesh" )
!              positions_old = extract_vector_field( packed_state, "Coordinate" )
!
!              PhaseV_old=>extract_tensor_field( packed_state, "PackedPhaseVolumeFraction"  )
!
!              sparsity=>extract_csr_sparsity(packed_state,"'CTMatrixSparsity'")
!
!                  !call allocate(petsc_acv,sparsity,[nphase,nphase],"ACV",.false.,.false.)
!
!              call allocate(mat,sparsity,[1,1],"CT",.false.,.false.)
!              call zero(mat)
!
!          else if (flag ==1) then
!
!
!              mesh_new => extract_mesh( packed_state, "CoordinateMesh" )
!              positions_new = extract_vector_field( packed_state, "Coordinate" )
!
!              !PhaseV_old=> ph_sol_new ! NEED TO POINT TO THE SOLUTION AFTER ADAPTING HERE
!
!              sparsity=>extract_csr_sparsity(packed_state,"ACVSparsity")
!                  !call allocate(petsc_acv,sparsity,[nphase,nphase],"ACV",.false.,.false.)
!
!              call allocate(mat,sparsity,[1,1],"ACV",.false.,.false.)
!              call zero(mat)
!
!          end if
!
!
!          do  ele = 1, totele
!
!              ! calculate detwei,ra,nx,ny,nz for element ele
!              call detnlxr_plus_u_with_storage( ele, x%val(1,:), x%val(2,:), x%val(3,:), &
!              x_ndgln, totele, x_nonods, x_nloc, tmp_cv_nloc, ph_ngi, &
!              tmp_cvfen, tmp_cvfenlx_all(1,:,:), tmp_cvfenlx_all(2,:,:), tmp_cvfenlx_all(3,:,:), &
!              tmp_cv_weight, detwei, ra, volume, d1, d3, dcyl, tmp_cvfenx_all, &
!              other_nloc, other_fenlx_all(1,:,:), other_fenlx_all(2,:,:), other_fenlx_all(3,:,:), &
!              other_fenx_all, state , "ph_2", StorageIndexes( 37 ) )
!
!              do cv_iloc = 1, cv_nloc
!
!                  cv_nodi = cv_ndgln(( ele - 1 ) * cv_nloc + cv_iloc )
!
!                  do cv_jloc = 1, cv_nloc
!
!                      cv_nodj = cv_ndgln(( ele - 1 ) * cv_nloc + cv_jloc )
!
!                      MN = 0.0
!                      MM = 0.0
!                      do cv_gi = 1, cv_ngi
!                          MN = MN + cvn( cv_iloc, cv_gi ) * N( cv_jloc, cv_gi )   * detwei( cv_gi )
!                          MM = MM + cvn( cv_iloc, cv_gi ) * CVN( cv_jloc, cv_gi ) * detwei( cv_gi )
!                      end do
!
!                      if(flag == 0) then
!                          call addto(mat,1, 1, cv_nodi,cv_nodj,MN)
!                          fempsi_rhs%val(cv_nodi) = fempsi_rhs%val(cv_nodi)+ MM*PhaseV_old%val(1, 1, cv_nodj)
!
!                      else if (flag == 1) then
!                          call addto(mat,1, 1, cv_nodi,cv_nodj,MM)
!                          fempsi_rhs%val(cv_nodi) = fempsi_rhs%val(cv_nodi)+ MN*PhaseV_old%val(1, 1, cv_nodj)
!
!                      end if
!
!                  enddo
!
!              enddo
!
!          enddo
!
!          ! Solve the matrix equation above and assign soluton to ph_sol_old - this is a global variable to this module
!
!          if (flag == 0) then
!
!              call set_solver_options( path, &
!              ksptype = "cg", &
!              pctype = "hypre", &
!              rtol = 1.0e-10, &
!              atol = 0.0, &
!              max_its = 10000 )
!              call add_option( &
!              trim( path ) // "/solver/preconditioner[0]/hypre_type[0]/name", stat )
!              call set_option( &
!              trim( path ) // "/solver/preconditioner[0]/hypre_type[0]/name", "boomeramg" )
!              ph_sol_old % option_path = path
!
!              call petsc_solve(ph_sol_old, mat, fempsi_rhs )
!
!              ! Insert the solutions into 'state_old'
!
!              call insert( state_old, mesh_old, "Mesh" )
!              call insert( state_old, positions_old, "Coordinate" )
!              call insert( state_old, ph_sol_old, "PhaseVolFrac" )
!
!          end if
!
!          if (flag == 1) then
!
!              call insert( state_new, mesh_new, "Mesh" )
!              call insert( state_new, positions_new, "Coordinate" )
!              call insert( state_new, ph_sol_new, "PhaseVolFrac" )
!
!              call interpolation_galerkin( state_old, state_new )
!
!              call set_solver_options( path, &
!              ksptype = "cg", &
!              pctype = "hypre", &
!              rtol = 1.0e-10, &
!              atol = 0.0, &
!              max_its = 10000 )
!              call add_option( &
!              trim( path ) // "/solver/preconditioner[0]/hypre_type[0]/name", stat )
!              call set_option( &
!              trim( path ) // "/solver/preconditioner[0]/hypre_type[0]/name", "boomeramg" )
!              ph_sol_old % option_path = path
!
!              call petsc_solve(ph_sol_new, mat, fempsi_rhs )
!
!          endif
!
!          call deallocate(mat)
!
!          deallocate(CVN, N)
!
!          ! Deallocate State types etc.
!
!      end subroutine
!
!
!  end module multi_interpolation
