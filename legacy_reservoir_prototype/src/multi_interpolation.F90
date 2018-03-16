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


module multi_interpolation

    use fldebug
    use state_module
    use fields
    use field_options
    use spud
    use populate_state_module
    use diagnostic_variables
    use diagnostic_fields
    use diagnostic_fields_wrapper
    use global_parameters, only: option_path_len, is_porous_media, is_multifracture
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
    use fields_allocates
    use futils, only: int2str


    use solvers
    use conservative_interpolation_module

    use shape_functions
    use shape_functions_Linear_Quadratic
    use solvers_module
    use matrix_operations
    use Copy_Outof_State
    use sparse_tools_petsc
    use shape_functions_prototype
    use sparsity_patterns, only: make_sparsity
    use multi_data_types

    ! Global variables to this module

    type(state_type) :: state_old, state_new

contains

  subroutine M2MInterpolation(state, packed_state, Mdims, CV_GIdims, CV_funs, small_finacv, &
            small_colacv, flag)
    implicit none
    ! IMPORTANT: flag is a switch before and after the adapt and tells us which interpolation step (1) or (3) to implement
    type( state_type ), dimension( : ), intent( inout ) :: state
    type( state_type ), intent( inout ) :: packed_state
    type(multi_dimensions), intent(in) :: Mdims
    type(multi_GI_dimensions), intent(in) :: CV_GIdims
    type(multi_shape_funs), intent(in) :: CV_funs
    integer, intent(in) :: flag
    integer, dimension(:), pointer, intent(inout) :: small_finacv, small_colacv
    !          ! local variables...checking
    type ( tensor_field ), pointer :: ufield
    integer :: ele, cv_iloc, cv_jloc, iphase, &            ! Leave iphase where it is for now (will probably need it for multiphase flow)
         tmp_cv_nloc, other_nloc
    integer, dimension( : ), pointer ::  x_ndgln, cv_ndgln, u_ndgln, dg_nodes
    logical :: quad_over_whole_ele, d1, d3, dcyl
    type( vector_field ), pointer :: x
    real, dimension( :, : ), pointer :: tmp_cvfen
    real, dimension( :, :, : ), pointer :: tmp_cvfenlx_all
    real, dimension( :, :, : ), pointer :: other_fenlx_all
    type( mesh_type ), pointer :: mesh_pres, mesh_pres_disc
    type( vector_field ) :: positions_old, positions_new
    integer :: cv_nodi, cv_nodj
    real :: MN, MM
    !Shape functions vars
    real, dimension( :, :, : ), allocatable :: tmp_cvfenx_all
    real, dimension( :, :, : ), allocatable :: other_fenx_all
    real, dimension( CV_GIdims%cv_ngi ) :: detwei, ra
    real :: volume
    ! Element by element inversion variables
    real, dimension(:,:), allocatable :: MMatrix, MNatrix
    real, dimension(:,:), allocatable :: EleRHS
    real, dimension(:,:), allocatable :: EleLHS
    integer, dimension(:), allocatable :: ipiv
    logical :: gotdec
    real, dimension(:,:), allocatable :: Long_EleRHS
    real, dimension(:), allocatable ::  mass_diag
    integer :: nfields, ifields, i, j
    character(len=10000) :: field_name
    type(scalar_field_pointer), dimension(:), allocatable :: scalar_field_list
    type(scalar_field), dimension(:), allocatable :: ph_sol_old, ph_sol_new, ph_sol_new_interm
    ! Finished variable declarations
    ! EXTRACT ALL FIELDS WHICH ARE TO HAVE CVGalerkin Interpolation APPLIED TO THEM - put in scalar_field_list
    nfields=option_count('/material_phase/scalar_field/prognostic/CVgalerkin_interpolation') ! Count # instances of CVGalerkin in the input file
    allocate(scalar_field_list(nfields))
    ifields=1
    do i = 1, size(state) ! Loop over all fields in state
       do j=1, option_count(trim(state(i)%option_path)//'/scalar_field') ! Loop over scalars
          if (have_option(trim(state(i)%option_path)//'/scalar_field['//int2str(j-1)//']/prognostic/CVgalerkin_interpolation') ) then ! Check if CVGalerkin is on
             call get_option(trim(state(i)%option_path)//'/scalar_field['//int2str(j-1)//']/name',field_name) ! If so, then grab the path of that scalar field
             scalar_field_list(ifields)%ptr=>extract_scalar_field(state(i),trim(field_name)) ! Put the scalar field into scalar_field_list
             ifields=ifields+1
          end if
       enddo
    enddo
    ufield => extract_tensor_field( packed_state, "PackedVelocity" )
    quad_over_whole_ele = .false.
    !**************************************
    !ALLOCATIONS
    allocate(EleRHS(nfields, Mdims%cv_nloc))
    allocate(EleLHS(nfields, Mdims%cv_nloc))
    allocate(MMatrix(Mdims%cv_nloc, Mdims%cv_nloc))
    allocate(MNatrix(Mdims%cv_nloc, Mdims%cv_nloc))
    allocate(ipiv(Mdims%cv_nloc))
    allocate(ph_sol_old(nfields))
    if(flag == 1) then
       allocate(ph_sol_new(nfields))
       allocate(ph_sol_new_interm(nfields))
    endif
    !**************************************
    x_ndgln => get_ndglno( extract_mesh( state( 1 ), "PressureMesh_Continuous" ) )
    cv_ndgln => get_ndglno( extract_mesh( state( 1 ), "PressureMesh" ) )
    x => extract_vector_field( packed_state, "PressureCoordinate" )
    u_ndgln => get_ndglno( extract_mesh( state( 1 ), "VelocityMesh" ) )
    d1 = ( Mdims%ndim == 1 ) ; d3 = ( Mdims%ndim == 3 ) ; dcyl = .false.
    ! Try to re-order these allocations - they need to be here so Mdims%cv_nonods is defined
    allocate(Long_EleRHS(nfields, Mdims%cv_nonods))
    allocate(mass_diag(Mdims%cv_nonods))
    ! ALLOCATE MEMORY FOR AN ARRAY OF SCALARS THAT WILL STORE OUR SOLUTION AT THE FIRST STEP
    if (flag == 0) then ! Meshes and Allocations for 1st Interpolation Calculation
       mesh_pres => extract_mesh( packed_state, "PressureMesh" )                     ! Strictly only need mesh_pres_disc in this case but keep it for parity
       mesh_pres_disc => extract_mesh( packed_state, "PressureMesh_Discontinuous" )
       positions_old = extract_vector_field( packed_state, "Coordinate" )            ! Check precisely what this "positions_old" is
       ! ALLOCATE ph_sol_old (Note state_old etc. allocated by the insert() lines later)
       do ifields = 1, nfields
          call allocate(ph_sol_old(ifields), mesh_pres_disc, "ph_sol_old" // int2str(ifields) ) ! Lives on finite element mesh ~ discontinuous control volume mesh
          call zero( ph_sol_old(ifields) )
       enddo
    endif ! end of if(flag == 0)
    ! NOTE THIS step (and all other flag == 1) steps are only called AFTER adapting the mesh - so the second time this routine is called.
    if(flag == 1) then
       ! Note: The call to the pure FE mapping through supermeshing happens here
       mesh_pres => extract_mesh( packed_state, "PressureMesh" )
       mesh_pres_disc => extract_mesh( packed_state, "PressureMesh_Discontinuous" )
       positions_new = extract_vector_field( packed_state, "Coordinate" )
       ! ALLOCATE MEMORY FOR SCALAR FIELD ARRAYS: ph_sol_new_interm, ph_sol_new
       do ifields = 1, nfields
          call allocate( ph_sol_new_interm(ifields), mesh_pres_disc, "ph_sol_interm" // int2str(ifields) ) ! Intermediate solution (i.e. after step 2 but before step 3)
          call zero( ph_sol_new_interm(ifields) )
          call allocate( ph_sol_new(ifields), mesh_pres, "ph_sol_new" // int2str(ifields) )                ! The solution ph_sol_new is on the pressure mesh (may or may not be discontinuous)
          call zero( ph_sol_new(ifields) )
       enddo
       ! SUPERMESHING I: Allocations for state_new (that will hold the results after the FEM mapping stage)
       do ifields = 1, nfields
          call insert( state_new, mesh_pres_disc, "Mesh" )
          call insert( state_new, positions_new, "Coordinate")
          call insert( state_new, ph_sol_new_interm(ifields), "Interpolant" // int2str(ifields) )
       enddo
       ! SUPERMESHING II: Call interpolation_galerkin() which via supermeshing constructs state_new from state_old (remember the latter is
       ! global to this module and is the result of the flag == 0 interpolation).
       call interpolation_galerkin(state_old, state_new)
       ! Extract all the fields from state_new having performed the supermeshing step
       do ifields = 1, nfields
          ph_sol_new_interm(ifields) = extract_scalar_field(state_new, "Interpolant" // int2str(ifields))
       enddo
    endif ! end of if(flag == 1)
    ! SETTINGS NEEDED TO CALCULATE detwei()
    if ( Mdims%cv_nloc == Mdims%u_nloc ) then
       tmp_cv_nloc = Mdims%u_nloc
       tmp_cvfen => CV_funs%ufen
       tmp_cvfenlx_all => CV_funs%ufenlx_all
       other_nloc = Mdims%cv_nloc
       other_fenlx_all => CV_funs%cvfenlx_all
    else
       tmp_cv_nloc = Mdims%cv_nloc
       tmp_cvfen => CV_funs%cvfen
       tmp_cvfenlx_all => CV_funs%cvfenlx_all
       other_nloc = Mdims%u_nloc
       other_fenlx_all => CV_funs%ufenlx_all
    end if
    ! INITIALISATIONS for the element loop
    !Allocate shape function vars
    allocate(tmp_cvfenx_all(3, size(tmp_cvfenlx_all,2), CV_GIdims%cv_ngi))
    allocate(other_fenx_all(3, size(other_fenlx_all,2), CV_GIdims%cv_ngi))
    EleLHS = 0
    if(flag ==1) then
       Long_EleRHS = 0.0
       mass_diag = 0.0
    endif

    do  ele = 1, Mdims%totele
       ! Calculate detwei related quantities
       call detnlxr_plus_u( ele, x%val(1,:), x%val(2,:), x%val(3,:), &
            x_ndgln, Mdims%totele, Mdims%x_nonods, Mdims%x_nloc, tmp_cv_nloc, CV_GIdims%cv_ngi, &
            tmp_cvfen, tmp_cvfenlx_all(1,:,:), tmp_cvfenlx_all(2,:,:), tmp_cvfenlx_all(3,:,:), &
            CV_funs%cvweight, detwei, ra, volume, d1, d3, dcyl, tmp_cvfenx_all, &
            other_nloc, other_fenlx_all(1,:,:), other_fenlx_all(2,:,:), other_fenlx_all(3,:,:), &
            other_fenx_all)
       ! LOOP to calculate the mass matrices and right hand side element by element and invert the linear problem.
       ! Problem is inverted element by element
       mesh_pres_disc => extract_mesh( packed_state, "PressureMesh_Discontinuous" )
       dg_nodes => ele_nodes(mesh_pres_disc, ele)       ! Replaces cv_nodi in DISCONTINUOUS cases. Extract discontinuous pressure mesh nodes
       EleRHS = 0.0                                     ! Can extract its components as dg_nodes(cv_iloc)
       do cv_iloc = 1, Mdims%cv_nloc
          cv_nodi = cv_ndgln(( ele - 1 ) * Mdims%cv_nloc + cv_iloc ) ! Remember this is CONTINUOUS numbering
          do cv_jloc = 1, Mdims%cv_nloc
             cv_nodj = cv_ndgln(( ele - 1 ) * Mdims%cv_nloc + cv_jloc )
             MN = sum( CV_funs%cvn( cv_iloc, : ) * CV_funs%cvfen( cv_jloc, : )   * detwei( : )  )
             MM = sum( CV_funs%cvn( cv_iloc, : ) * CV_funs%cvn( cv_jloc, : )   * detwei( : ) )
             if(flag == 0) then
                ! Matrices for element by element inversion
                MNatrix(cv_iloc,cv_jloc) = MN
                do ifields = 1, nfields
                   EleRHS(ifields, cv_iloc) = EleRHS(ifields, cv_iloc) + MM*scalar_field_list(ifields)%ptr%val(cv_nodj)
                enddo
             else if(flag ==1) then
                ! Note we cannot invert this case element by element - hence the global storage
                !(the reason being that continuous CVs span multiple elements)
                mass_diag(cv_nodi) = mass_diag(cv_nodi)+ MM
                do ifields = 1, nfields
                   Long_EleRHS(ifields, cv_nodi) = Long_EleRHS(ifields, cv_nodi) + MN*ph_sol_new_interm(ifields)%val(dg_nodes(cv_jloc))
                enddo
             endif
          enddo ! cv_jloc loop
       enddo  ! cv_iloc loop
       ! Solve matrix inversion problem element-wise (in the first case i.e. flag == 0). Second case dealt with later
       if(flag == 0) then ! Solve the element-wise matrix problem MMatrix*EleLHS = EleRHS for EleLHS
          gotdec = .false.
          do ifields = 1, nfields
             call SMLINNGOT(MNatrix, EleLHS(ifields,:), EleRHS(ifields,:), Mdims%cv_nloc, ipiv, gotdec )
             gotdec = .true.
          enddo
          ! Append this solution to the global ph_sol_old
          do cv_iloc = 1, Mdims%cv_nloc
             dg_nodes => ele_nodes(mesh_pres_disc, ele)
             do ifields = 1, nfields
                ph_sol_old(ifields)%val(dg_nodes(cv_iloc)) = EleLHS(ifields, cv_iloc)
             enddo
          enddo
       endif
    enddo ! End of loop over ele
    if (flag == 0) then  ! Insert the solutions into 'state_old' ready for projection onto the new mesh in step (2), flag == 1
       positions_old = extract_vector_field( packed_state, "Coordinate" )
       do ifields = 1, nfields
          call insert( state_old, mesh_pres_disc, "Mesh")
          call insert( state_old, positions_old, "Coordinate")
          call insert( state_old, ph_sol_old(ifields), "Interpolant" // int2str(ifields) )
       enddo
    end if
    if (flag == 1) then  ! Solve the final equation (3) to map back to a CV representation. Assign this value to ph_sol_new
       do ifields = 1, nfields
          ph_sol_new(ifields)%val(:) = Long_EleRHS(ifields, :)/mass_diag(:)
       enddo
       ! Copy the output values back into state - DONE!
       do ifields = 1, nfields
          scalar_field_list(ifields)%ptr%val = ph_sol_new(ifields)%val
       enddo
    endif
    ! BOUNDEDNESS : Our solutions are currently not bounded to be in [0,1]. The following subroutine call should fix this
    ! This section needs to be generalised to work for multi-fields (I think the boundedness subroutine may need generalisation)
    !print *, nfields
    if (have_option('/material_phase[0]/scalar_field::Temperature/prognostic/CVgalerkin_interpolation')) then
       if(flag == 1) call BoundedSolutionCorrections(state, packed_state, Mdims, CV_funs, small_finacv, small_colacv)
    else if(have_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/CVgalerkin_interpolation')) then
       if(flag == 1) call BoundedSolutionCorrections(state, packed_state, Mdims, CV_funs, small_finacv, small_colacv,.true.)
    endif
    ! DEALLOCATIONS
    deallocate(EleLHS, EleRHS, MMatrix, MNatrix, ipiv)
    deallocate(Long_EleRHS, mass_diag)
    deallocate(scalar_field_list)
    deallocate(tmp_cvfenx_all, other_fenx_all)
    if(flag == 0) then
       do ifields = 1, nfields
          call deallocate(ph_sol_old(ifields)) ! Check theoretically why it has to be call deallocate, then deallocate (in that order)
       enddo
       deallocate(ph_sol_old)
    endif
    if (flag == 1) then
       do ifields = 1, nfields
          call deallocate(ph_sol_new_interm(ifields))
          call deallocate(ph_sol_new(ifields))
       enddo
       deallocate(ph_sol_new_interm, ph_sol_new)
    endif
  end subroutine M2MInterpolation

    subroutine MemoryCleanupInterpolation1()

        call deallocate(state_old)

    end subroutine

    subroutine MemoryCleanupInterpolation2()

        call deallocate(state_new)

    end subroutine

end module multi_interpolation
