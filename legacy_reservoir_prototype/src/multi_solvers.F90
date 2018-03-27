
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

module solvers_module

    use fldebug
    use fields
    use Petsc_tools
    use sparse_tools_petsc
    use solvers
    use global_parameters, only: OPTION_PATH_LEN, FPI_have_converged
    use spud

    use state_module
    use halo_data_types
#ifdef HAVE_PETSC_MODULES
  use petsc 
#if PETSC_VERSION_MINOR==0
  use petscvec 
  use petscmat 
  use petscksp 
  use petscpc
#endif
#endif
    use Copy_Outof_State
    use shape_functions_Linear_Quadratic
    use shape_functions_prototype
    use multi_data_types

    implicit none

#include "petsc_legacy.h"

    private

    public :: multi_solver, BoundedSolutionCorrections, FPI_backtracking, Set_Saturation_to_sum_one,&
         Ensure_initial_Saturation_to_sum_one, auto_backtracking

    interface multi_solver
        module procedure solve_via_copy_to_petsc_csr_matrix
    end interface
  
contains

    ! -----------------------------------------------------------------------------
    subroutine solve_via_copy_to_petsc_csr_matrix( A, &
        x, b, findfe, colfe, option_path, block_size )

        !!< Solve a matrix Ax = b system via copying over to the
        !!< petsc_csr_matrix type and calling the femtools solver
        !!< using the spud options given by the field options path

        integer, dimension( : ), intent(in) :: findfe
        integer, dimension( : ), intent(in) :: colfe
        real, dimension( : ), intent(in) :: a, b
        real, dimension( : ), intent(inout) :: x
        character( len=* ), intent(in) :: option_path
        !This optional argument is to create a block-diagonal matrix
        !and therefore to use block-solvers
        integer, optional, intent(in) :: block_size

        ! local variables
        integer :: i, j, k, rows
        integer, dimension( : ), allocatable :: dnnz
        type(petsc_csr_matrix) :: matrix
        integer :: size_of_block

        size_of_block = 1
        if (present(block_size)) size_of_block = block_size

        rows = size( x )
        assert( size( x ) == size( b ) )
        assert( size( a ) == size( colfe ) )
        assert( size( x ) + 1 == size( findfe ) )
        ewrite(3,*) rows+1, size(findfe)

        allocate( dnnz( rows/size_of_block ) ) ; dnnz = 0
        ! find the number of non zeros per row
        do i = 1, size( dnnz )
            dnnz( i ) =(findfe( i+1 ) - findfe( i ))
        end do
        call allocate( matrix, rows/size_of_block, rows/size_of_block, dnnz, dnnz,&
            (/1, 1/), name = 'dummy', element_size=size_of_block)

        call zero( matrix )
        ! add in the entries to petsc matrix
        do i = 1, rows
            do j = findfe( i ), findfe( i+1 ) - 1
                k = colfe( j )
                call addto( matrix, blocki = 1, blockj = 1, i = i, j = k, val = a( j ) )
            end do
        end do

        call assemble( matrix )

        call petsc_solve_scalar_petsc_csr_mp( x, matrix, b, rows, trim( option_path ) )

        ! deallocate as needed
        deallocate( dnnz )
        call deallocate( matrix )

        return
    end subroutine solve_via_copy_to_petsc_csr_matrix


    subroutine petsc_solve_scalar_petsc_csr_mp( x, matrix, rhs, rows, option_path )

        real, dimension( : ), intent(inout) :: x
        type( petsc_csr_matrix ), intent(inout) :: matrix
        real, dimension( : ), intent(in) :: rhs
        integer, intent(in) :: rows
        character( len=* ), intent(in) :: option_path

        KSP :: ksp
        Vec :: y, b

        character(len=OPTION_PATH_LEN) :: solver_option_path
        integer :: ierr

        assert( size( x ) == size( rhs ) )
        assert( size( x ) == size( matrix, 2 ) )
        assert( size( rhs ) == size( matrix, 1 ) )

        solver_option_path = complete_solver_option_path( option_path )

        call SetupKSP( ksp, matrix%M, matrix%M, solver_option_path, .false., &
            matrix%column_numbering, .true. )

        b = PetscNumberingCreateVec( matrix%column_numbering )
        call VecDuplicate( b, y, ierr )


        ! copy array into PETSc vecs
        call VecSetValues( y, rows, &
            matrix%row_numbering%gnn2unn( 1:rows, 1 ), &
            x, INSERT_VALUES, ierr )
        call VecAssemblyBegin( y, ierr )
        call VecAssemblyEnd( y, ierr )

        call VecSetValues( b, rows, &
            matrix%row_numbering%gnn2unn( 1:rows, 1 ), &
            rhs, INSERT_VALUES, ierr )
        call VecAssemblyBegin( b, ierr )
        call VecAssemblyEnd( b, ierr )

        call KSPSolve( ksp, b, y, ierr )

        ! copy back the result
        call VecGetValues( y, rows, &
            matrix%row_numbering%gnn2unn( 1:rows, 1 ), &
            x, ierr )

        ! destroy all PETSc objects and the petsc_numbering
        call multi_petsc_solve_destroy_petsc_csr( y, b, ksp )


        return
        contains
        !Clone of the same subroutine in femtools/Solvers.F90
        subroutine multi_petsc_solve_destroy_petsc_csr( y, b, ksp )

            type(Vec), intent(inout):: y
            type(Vec), intent(inout):: b
            type(KSP), intent(inout):: ksp

            type(PC) :: pc
            integer ierr

            call VecDestroy(y, ierr)
            call VecDestroy(b, ierr)
            call KSPGetPC(ksp, pc, ierr)
            call KSPDestroy(ksp, ierr)

        end subroutine multi_petsc_solve_destroy_petsc_csr
    end subroutine petsc_solve_scalar_petsc_csr_mp


    subroutine BoundedSolutionCorrections( state, packed_state, &
        Mdims, CV_funs, small_findrm, small_colm, &
        for_sat)
        implicit none
        ! This subroutine adjusts field_val so that it is bounded between field_min, field_max in a local way.
        ! The sparcity of the local CV connectivity is in: small_findrm, small_colm.
        ! ngl_its=max no of global iterations e.g. 100.
        ! error_tol = tolerance on the iterations.
        !
        ! nloc_its: This iteration is very good at avoiding spreading the modifications too far - however it can stagnate.
        ! nloc_its2: This iteration is very good at avoiding stagnating but does spread the modifcations far.
        ! us a single iteration because of this as default...
        ! nits_nod: iterations at a nod - this iteration is very good at avoiding spreading the modifications too far -
        ! however it can stagnate.
        integer, parameter :: nloc_its = 5, nloc_its2 = 1, nits_nod = 100, ngl_its = 500
        real, parameter :: w_relax = 0.5, error_tol = 1.0e-5
        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_shape_funs), intent(in) :: CV_funs
        integer, dimension( : ), intent( in ) :: small_findrm, small_colm
        logical, optional, intent(in) :: for_sat
        ! local variables...
        type (multi_dev_shape_funs) :: DevFuns
        type ( tensor_field ), pointer :: field, ufield
        real, dimension( :, :, : ), allocatable :: field_dev_val, field_alt_val, field_min, field_max
        real, dimension( :, : ), allocatable :: scalar_field_dev_max, scalar_field_dev_min
        real, dimension( :, : ), allocatable :: r_min, r_max
        integer, dimension( :, : ), allocatable :: ii_min, ii_max
        real, dimension( : ), allocatable :: mass_cv, mass_cv_sur
        integer :: ndim1, ndim2,  i, j, knod, inod, jnod, count, ii, jj, loc_its, loc_its2, its, gl_its, cv_iloc
        logical :: changed, changed_something
        real :: max_change, error_changed, max_max_error, scalar_field_dev, mass_off, alt_max, alt_min
        integer :: ele, iloc, jloc
        real :: mm
        integer, dimension( : ), pointer ::  x_ndgln, cv_ndgln
        type(scalar_field) :: mass_cv_sur_halo
        type( vector_field ), pointer :: x
        real, dimension(:,:), pointer :: Immobile_fraction
        if (present_and_true(for_sat)) then
            field => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
        else
            field => extract_tensor_field( packed_state, "PackedTemperature" )
        end if
        ndim1 = size( field%val, 1 ) ; ndim2 = size( field%val, 2 ) ;
        ewrite(3,*) 'Bounding correction input: iphase, icomp, min, max:'
        do j = 1, ndim2
            do i = 1, ndim1
                ewrite(3,*) i, j, minval( field%val( i, j, : ) ), maxval( field%val( i, j, : ) )
            end do
        end do
        allocate( field_dev_val( ndim1, ndim2, Mdims%cv_nonods ), field_alt_val( ndim1, ndim2, Mdims%cv_nonods ) )
        allocate( field_min( ndim1, ndim2, Mdims%cv_nonods ), field_max( ndim1, ndim2, Mdims%cv_nonods ) )
        allocate( scalar_field_dev_max( ndim1, ndim2 ), scalar_field_dev_min( ndim1, ndim2 ) )
        allocate( r_min( ndim1, ndim2 ), r_max( ndim1, ndim2 ) )
        allocate( ii_min( ndim1, ndim2 ), ii_max( ndim1, ndim2 ) )
        allocate( mass_cv( Mdims%cv_nonods ), mass_cv_sur( Mdims%cv_nonods ) )
        call allocate(mass_cv_sur_halo,field%mesh,'mass_cv_sur_halo')
        ufield => extract_tensor_field( packed_state, "PackedVelocity" )

        x_ndgln => get_ndglno( extract_mesh( state( 1 ), "PressureMesh_Continuous" ) )
        cv_ndgln => get_ndglno( extract_mesh( state( 1 ), "PressureMesh" ) )
        x => extract_vector_field( packed_state, "PressureCoordinate" )
        mass_cv = 0.0
        call allocate_multi_dev_shape_funs(CV_funs, DevFuns)
        do  ele = 1, Mdims%totele
            !Retrieve DevFuns%detwei
            call DETNLXR(ele, X%val, x_ndgln, CV_funs%cvweight, CV_funs%CVFEN, CV_funs%CVFENLX_ALL, DevFuns)

            do iloc = 1, Mdims%cv_nloc
                inod = cv_ndgln( ( ele - 1 ) * Mdims%cv_nloc + iloc )
                do jloc = 1, Mdims%cv_nloc
                    mm = sum( CV_funs%cvn( iloc, : ) * CV_funs%cvn( jloc, : ) * DevFuns%detwei )
                    mass_cv( inod ) = mass_cv( inod ) + mm
                end do
            end do
        end do
        mass_cv_sur = 0.0
        do inod = 1, Mdims%cv_nonods
            if ( .not. node_owned( field, inod ) ) cycle
            do count = small_findrm( inod ), small_findrm( inod + 1 ) - 1
                jnod = small_colm( count )
                mass_cv_sur(inod) = mass_cv_sur(inod) + mass_cv( jnod )
            end do
        end do
        ! Obtain the halos of mass_cv_sur:
        mass_cv_sur_halo%val(:)=mass_cv_sur(:)
        call halo_update(mass_cv_sur_halo)
        mass_cv_sur(:)=mass_cv_sur_halo%val(:)
        !Establish bounds
        if (present_and_true(for_sat)) then
            !Define the immobile fractions for each phase
            call get_var_from_packed_state(packed_state, immobile_fraction = immobile_fraction)

            do ele = 1, Mdims%totele
                do cv_iloc = 1, Mdims%cv_nloc
                    inod = cv_ndgln(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                    do ii = 1, ndim2!phases
                        field_min(:, ii, inod) = immobile_fraction(ii, ele)
                        field_max(:, ii, inod) = 1.0 - sum(immobile_fraction(:,ele))&
                            + immobile_fraction(ii,ele)
                    end do
                end do
            end do
        else
            field_min = 0.0 ; field_max = 1.0
        end if
        do gl_its = 1, ngl_its
            ! This iteration is very good at avoiding spreading the modifications too far - however it can stagnate.
            max_change = 0.0
            do loc_its = 1, nloc_its
                changed_something = .false.
                do knod = 1, Mdims%cv_nonods ! exclude the halo values of knod for parallel
                    if ( .not. node_owned( field, knod ) ) cycle
                    do its = 1, nits_nod
                        r_min = 0.0 ; r_max = 0.0
                        ii_min = 0 ; ii_max = 0
                        scalar_field_dev_max = 0.0 ; scalar_field_dev_min = 0.0
                        ! Find the max and min deviation from the limited values....
                        do count = small_findrm( knod ), small_findrm( knod + 1 ) - 1
                            jnod = small_colm( count )
                            if ( .not. node_owned( field, jnod ) ) cycle
                            do j = 1, ndim2
                                do i = 1, ndim1
                                    if ( field%val( i, j, jnod ) > field_max( i, j, jnod ) ) then
                                        scalar_field_dev = field%val( i, j, jnod ) - field_max( i, j, jnod )
                                    else if ( field%val( i, j, jnod ) < field_min( i, j, jnod ) ) then
                                        scalar_field_dev = field%val( i, j, jnod ) - field_min( i, j, jnod )
                                    else
                                        scalar_field_dev = 0.0
                                    end if
                                    if ( scalar_field_dev * mass_cv ( jnod ) > r_max( i, j ) ) then
                                        scalar_field_dev_max( i, j ) = scalar_field_dev
                                        r_max( i, j ) = scalar_field_dev * mass_cv( jnod )
                                        ii_max( i, j ) = jnod
                                    end if
                                    if ( scalar_field_dev * mass_cv ( jnod ) < r_min( i, j ) ) then
                                        scalar_field_dev_min( i, j ) = scalar_field_dev
                                        r_min( i, j ) = scalar_field_dev * mass_cv( jnod )
                                        ii_min( i, j ) = jnod
                                    end if
                                end do
                            end do
                        end do ! do count = small_findrm( knod ), small_findrm( knod + 1 ) - 1
                        changed=.false.
                        ! Change the max and min limited deviation by sharing the deviation between them...
                        do j = 1, ndim2
                            do i = 1, ndim1
                                ii = ii_max( i, j )
                                jj = ii_min( i, j )
                                if ( ii /= 0 .and. jj /= 0 ) then
                                    if ( abs( r_max( i, j ) ) > abs( r_min( i, j ) ) ) then
                                        alt_max = ( r_max( i, j ) + r_min( i, j ) ) / mass_cv( ii )
                                        alt_min = 0.0
                                    else
                                        alt_max = 0.0
                                        alt_min = ( r_max( i, j ) + r_min( i, j ) ) / mass_cv( jj )
                                    end if
                                    max_change = max( max_change, abs( -scalar_field_dev_max( i, j ) + alt_max ) )
                                    max_change = max( max_change, abs( -scalar_field_dev_min( i, j ) + alt_min ) )
                                    field%val( i, j, ii ) = field%val( i, j, ii ) - scalar_field_dev_max( i, j ) + alt_max
                                    field%val( i, j, jj ) = field%val( i, j, jj ) - scalar_field_dev_min( i, j ) + alt_min
                                    changed = .true.
                                    changed_something = .true.
                                end if
                            end do
                        end do
                        if ( .not. changed ) exit ! stop iterating and move onto next node/CV...
                    end do ! do its=1,nits_nod
                end do ! do knod = 1, Mdims%cv_nonods
                if ( .not. changed_something ) exit ! stop iterating and move onto next stage of iteration...
            end do ! do loc_its=1,nloc_its
            call halo_update( field )
            ! This iteration is very good at avoiding stagnating but does spread the modifcations far.
            ! use a single iteration because of this as default...
            do loc_its2 = 1, nloc_its2
                do knod = 1, Mdims%cv_nonods
                    !               if ( .not. node_owned( field, knod ) ) cycle
                    do j = 1, ndim2
                        do i = 1, ndim1
                            if ( field%val( i, j, knod ) > field_max( i, j, knod ) ) then
                                field_dev_val( i, j, knod ) = field%val( i, j, knod ) - field_max( i, j, knod )
                            else if ( field%val( i, j, knod ) < field_min( i, j, knod ) ) then
                                field_dev_val( i, j, knod ) = field%val( i, j, knod ) - field_min( i, j, knod )
                            else
                                field_dev_val( i, j, knod ) = 0.0
                            end if
                        end do
                    end do
                end do
                ! matrix vector...
                field_alt_val = 0.0
                do inod = 1, Mdims%cv_nonods ! exclude the halo values of knod for parallel
                    if ( .not. node_owned( field, inod ) ) cycle
                    do count = small_findrm( inod ), small_findrm( inod + 1 ) - 1
                        jnod = small_colm( count )
                        mass_off = mass_cv( jnod ) / mass_cv_sur( jnod )
                        field_alt_val( :, :, inod ) = field_alt_val( :, :, inod ) + mass_off * field_dev_val( :, :, jnod )
                    end do ! do count = small_findrm( inod ), small_findrm( inod + 1 ) - 1
                end do ! do inod = 1, Mdims%cv_nonods
                ! w_relax\in[0,1]: - This relaxation is used because we
                ! have used a mass matrix which is not diagonally dominant
                ! =0.5 is suggested.
                ! =1.0 is no relaxation.
                field_alt_val = w_relax * field_alt_val + ( 1.0 - w_relax ) * field_dev_val
                ! adjust the values...
                error_changed = maxval( abs( -field_dev_val + field_alt_val ) )
                field%val( :, :, : ) = field%val( :, :, : ) - field_dev_val( :, :, : ) + field_alt_val( :, :, : )
                call halo_update( field )
            end do ! loc_its2
            ! communicate the errors ( max_change, error_changed ) ...
            ! this could be more efficient sending a vector...
            max_max_error= max( max_change, error_changed )
            call allmax( max_max_error )
            if ( max_max_error < error_tol ) exit
        end do ! gl_its
        ewrite(3,*) 'Bounding correction output: iphase, icomp, min, max:'
        do j = 1, ndim2
            do i = 1, ndim1
                ewrite(3,*) i, j, minval( field%val( i, j, : ) ), maxval( field%val( i, j, : ) )
            end do
        end do
        deallocate( field_dev_val, field_alt_val )
        deallocate( field_min, field_max )
        deallocate( scalar_field_dev_max, scalar_field_dev_min )
        deallocate( r_min, r_max )
        deallocate( ii_min, ii_max )
        deallocate( mass_cv, mass_cv_sur )
        call deallocate(mass_cv_sur_halo)
        call deallocate_multi_dev_shape_funs(DevFuns)
        return
    end subroutine BoundedSolutionCorrections

    !sprint_to_do!not use one global variable
    subroutine FPI_backtracking(Mdims, ndgln, state, packed_state, sat_bak, backtrack_sat, backtrack_par_from_schema, &
        Previous_convergence, satisfactory_convergence, new_backtrack_par, Max_sat_its, its, nonlinear_iteration, useful_sats, res, &
        res_ratio, first_res)
        !In this subroutine we applied some corrections and backtrack_par on the saturations obtained from the saturation equation
        !this idea is based on the paper SPE-173267-MS.
        !The method ensures convergence "independent" of the time step.
        implicit none
        !Global variables
        type( multi_dimensions ), intent( in ) :: Mdims
        type(multi_ndgln), intent(in) :: ndgln
        type( state_type ), dimension( : ), intent( in ) :: state
        type( state_type ), intent(inout) :: packed_state
        real, dimension(:, :), intent(in) :: sat_bak, backtrack_sat
        real, intent(in) :: backtrack_par_from_schema, res, res_ratio, first_res
        logical, intent(inout) :: satisfactory_convergence
        real, intent(inout) :: new_backtrack_par, Previous_convergence
        integer, intent(in) :: Max_sat_its, its, nonlinear_iteration
        integer, intent(inout) :: useful_sats
        !Local parameters
        real, parameter :: Conv_to_achiv = 10.0
        real, save :: anders_exp!This parameter change the importance of backtrack_sat in Anderson's acceleration (mainly for high alphas)
        !Local variables        !100 => backtrack_sat is not used; 0.3 => equally important; 0.4 => recommended; 0 => more important than sat_bak
        real, dimension(:, :), pointer :: Satura
        logical :: new_time_step, new_FPI
        real :: aux
        integer :: i
        !Parameters for the automatic backtrack_par
!        integer, parameter :: History_order = 4!<= Cubic
        integer, parameter :: History_order = 3!<= Quadratic
        real :: min_backtrack
        real, dimension(History_order+1), save :: backtrack_pars = -1
        real, dimension(History_order), save :: Convergences = -1
        real, dimension(History_order) :: Coefficients
        logical, save :: allow_undo = .true.
        real, save :: convergence_tol
        real, save :: Infinite_norm_tol
        type (tensor_field), pointer :: sat_field
        !Initialize variables
        new_backtrack_par = 1.0
        new_FPI = (its == 1); new_time_step = (nonlinear_iteration == 1)
        !First, impose physical constrains
        if (is_porous_media) then
            call Set_Saturation_to_sum_one(mdims, ndgln, state, packed_state)
            sat_field => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
            Satura =>  sat_field%val(1,:,:)
            !Stablish minimum backtracking parameter
            min_backtrack = 0.1
        else
            !Use the pressure as it is in this case the field of interest
            sat_field => extract_tensor_field( packed_state, "PackedFEPressure" )
            Satura =>  sat_field%val(1,:,:)
            !Stablish minimum backtracking parameter
            min_backtrack = 0.05
        end if

        !Automatic method based on the history of convergence
        if (backtrack_par_from_schema < 0.0) then

            !Retrieve convergence factor, to make sure that if between time steps things are going great, we do not reduce the
            !backtrack_par_parameter
            if (backtrack_pars(1) < 0) then!retrieve it just once
                call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration',&
                    convergence_tol, default = 0. )
                convergence_tol =  convergence_tol
                !Tolerance for the infinite norm
                call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Infinite_norm_tol',&
                    Infinite_norm_tol, default = 0.03 )
                backtrack_pars(1) = max(min(abs(backtrack_par_from_schema), 1.0), min_backtrack)
                !Retrieve the shape of the function to use to weight the importance of previous saturations
                call get_option( '/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Acceleration_exp',&
                    anders_exp, default = 0.4 )!0.4 the best option, based on experience
                !Should not be negative
                anders_exp = max(anders_exp, 0.)
            end if

            if (new_time_step .and.is_porous_media) then
                !Store last convergence to use it as a reference
                Previous_convergence = Convergences(1)
                !restart all the storage
                backtrack_pars(2:) = -1; Convergences = -1
                !###New time step###
                !First FPI with the parameter introduced by the user
                backtrack_pars(1) = max(min(abs(backtrack_par_from_schema), 1.0), min_backtrack)
                satisfactory_convergence = .true.
            else
                !Store convergence obtained with the previous backtrack_par parameter
                Convergences(1) = res_ratio!<=Actual residual ratio
                !Compare with the Convergence using the first backtrack_par parameter
                if (new_FPI) Previous_convergence = Convergences(1)
                !####Check convergence of the method####
                satisfactory_convergence = (its > Max_sat_its) .or. (first_res / res > Conv_to_achiv) &
                    .or. (get_Convergence_Functional(Satura, Sat_bak, backtrack_pars(2)) < convergence_tol .and.&
                    maxval(abs(Sat_bak-Satura))/backtrack_pars(2) < Infinite_norm_tol)!<= exit if final convergence is achieved
                if (IsParallel()) call alland(satisfactory_convergence)
                !If a backtrack_par parameter turns out not to be useful, then undo that iteration
                if (its > 2 .and. Convergences(2) > 0 .and. allow_undo .and. Convergences(1)>5.) then
                    Satura = backtrack_sat
                    !We do not allow two consecutive undos
                    allow_undo = .false.
                    !restart the counter of useful saturations, this is to not use Anderson acceleration
                    !with values that were rejected
                    useful_sats = 1
                    !Undo de accumulation of backtrack_par parameters
                    new_backtrack_par = - backtrack_pars(1)
                    return
                else
                    allow_undo = .true.
                end if

                !Select different backtrack_par parameter for the first Saturation iteration (SFPI)
                select case (its)
                    case (1)
                        !First, use the one introduced by the user, this is necessary since the algorithm tends to small
                        !values of alphas
                        backtrack_pars(1) = max(min(abs(backtrack_par_from_schema), 1.0), min_backtrack)
                    case default
                        !Calculate a curve that fits the historical data
                        call Cubic_fitting(backtrack_pars(2:), Convergences, Coefficients)
                        !Calculate the new optimal backtrack_par parameter
                        backtrack_pars(1) = get_optimal_backtrack_par(backtrack_pars(2:), Convergences, Coefficients)
                end select

                !Update history, 1 => newest
                do i = size(backtrack_pars), 2, -1
                    backtrack_pars(i) = backtrack_pars(i-1)
                end do
                do i = size(Convergences), 2, -1
                    Convergences(i) = Convergences(i-1)
                end do

            end if
        else!Use the value introduced by the user
            backtrack_pars(1) = backtrack_par_from_schema
            !Just one local saturation iteration
            satisfactory_convergence = .true.
        end if

        ewrite(1,*) "backtrack_par factor",backtrack_pars(1)

        !If it is parallel then we want to be consistent between cpus
        !we use the smallest value, since it is more conservative
        if (IsParallel()) then
            call allmin(backtrack_pars(1))
            call allmin(Convergences(1))
        end if
        !***Calculate new saturation***
        if (is_porous_media) then
            !Obtain new saturation using the backtracking method
            if (useful_sats < 2 .or. satisfactory_convergence) then
                !Since Anderson's acceleration is unstable, when it has converged, we use the stable form of backtracking
                Satura = sat_bak * (1.0 - backtrack_pars(1)) + backtrack_pars(1) * Satura
            else !Use Anderson acceleration, idea from "AN ACCELERATED FIXED-POINT ITERATION FOR SOLUTION OF VARIABLY SATURATED FLOW"

                !Based on making backtrack_sat small when backtrack_pars(1) is high and backtrack_sat small when backtrack_pars(1) is small
                !The highest value of backtrack_sat is displaced to low values of alpha
                aux = 1.0 - backtrack_pars(1)
                Satura = backtrack_pars(1) * Satura + aux * ( (1.-(aux**anders_exp *backtrack_pars(1)) ) * sat_bak + &
                    aux**anders_exp *backtrack_pars(1) * backtrack_sat)!<=The best option so far

            end if
            !Update halos with the new values
            if (IsParallel()) call halo_update(sat_field)
        end if
        !Inform of the new backtrack_par parameter used
        new_backtrack_par = backtrack_pars(1)

    contains

        real function get_optimal_backtrack_par(backtrack_pars, Convergences, Coefficients)
            implicit none
            real, dimension(:), intent(in) ::backtrack_pars, Convergences
            real, dimension(:), intent(inout) :: Coefficients
            !Local variables
            integer, dimension(1) :: i
            integer :: n
            logical :: Basic_method
            real :: aux
            real, dimension(2) :: X, Y2


            Basic_method = .false.
            get_optimal_backtrack_par = backtrack_pars(1)


            !Check how much data we have
            do n = 1, size(backtrack_pars)
                if (backtrack_pars(n) < 0) exit
            end do
            n = n - 1

            !Check coefficients
            if ( ISNAN(sum(Coefficients(1:n)))) Basic_method = .true.

            if (.not. Basic_method) then
                select case (n)!quadratic
                    case (3)!Quadratic
                        get_optimal_backtrack_par = abs(-Coefficients(2)/(2.0*Coefficients(1)))
                        !If it is a maximum then go to the basic method
                        if (Coefficients(1) > 0) Basic_method = .true.

                    case (4)!Cubic !Very Unstable, better just quadratic
                        !Get solution of quadratic system
                        aux = 4. * Coefficients(2)**2 - 12. * Coefficients(1) * Coefficients(3)
                        if (aux > 0) then
                            aux = sqrt(aux)
                        else
                            Basic_method = .true.
                        end if
                        !Calculate first derivative
                        X(1) = (-2. * Coefficients(2) + aux)/(6.0 * Coefficients(1))
                        X(2) = (-2. * Coefficients(2) - aux)/(6.0 * Coefficients(1))
                        Y2(:) = 6.0 * Coefficients(1) * X(:)+ 2.0 * Coefficients(2)
                        if (Y2(1) > 0 .and. Y2(2) > 0) then
                            !Get the optimal value
                            i = minloc(Y2, MASK = Y2>0)
                            i(1) = max(i(1),1)
                            get_optimal_backtrack_par = X(i(1))
                        else if (Y2(1) > 0) then
                            get_optimal_backtrack_par = X(1)
                        else if(Y2(2) > 0) then
                            get_optimal_backtrack_par = X(2)
                        else!No minimums, hence go for the basic method
                            Basic_method = .true.
                        end if
                        !This method tend to be unstable, if we get something out
                        !of bounds, go back to the simple method
                        if (get_optimal_backtrack_par < 0 .or. get_optimal_backtrack_par > 1) Basic_method = .true.
                        !Test positive value
                        if (.not. Basic_method) then
                            aux = get_optimal_backtrack_par
                            aux = Coefficients(1) * aux**3  + Coefficients(2) * aux**2 + Coefficients(3) * aux + Coefficients(4)
                            if (aux < 0) then!Just for testing purposes
                                get_optimal_backtrack_par = get_optimal_backtrack_par/2.0
                            end if
                        end if

                    case default!Linear or constant
                        Basic_method = .true.
                end select
            end if

            !If results obtained are not within range, use the basic method
            if (get_optimal_backtrack_par > 1.0 .or. get_optimal_backtrack_par < 0) Basic_method = .true.

            !If not possible to get a good value, then just a simple method based on the history
            if (basic_method) then
                if (backtrack_pars(3) < 0 .or..true.) then!SIMPLE METHOD
                    if (Convergences(1)-Convergences(2) < 0) then!Converging
                        get_optimal_backtrack_par = backtrack_pars(1) * 1.1
                    else!Diverging, the reduce with a minimum value that will mean performing all the non-linear iterations
                        get_optimal_backtrack_par = backtrack_pars(1) / 2.0!1.5
                    end if
                else
                    X = (/backtrack_pars(1), Convergences(1) /) - (/backtrack_pars(2), Convergences(2) /)
                    Y2 = (/backtrack_pars(2), Convergences(2) /) - (/backtrack_pars(3), Convergences(3) /)
                    if (X(2)/X(1) < 0) then!Negative slope => Converging
                        if (Y2(2)/Y2(1) < 0) then!It was converging already
                            if (abs(X(2)/X(1)) > abs(Y2(2)/Y2(1))) then!New slope is steeper => Optimal still to be reached
                                !New value is previous value + the vector divided by the tangent, if it is very steep then
                                !the optimal is very close

                                get_optimal_backtrack_par = backtrack_pars(1) * max(1.2, 2.0 * X(1))

                            else!Old slope was steeper => Optimal value in between the previous backtrack_par parameters
                                get_optimal_backtrack_par = min(backtrack_pars(1) / 2, 0.5 * (backtrack_pars(2) + backtrack_pars(3)))

                            end if
                        else!It started to converge now, so we encourage to get away from the bad convergence value
                            get_optimal_backtrack_par = backtrack_pars(1) * max(1.2, 2.0 * X(1))
                        end if
                    else!It is NOT converging now
                        if (Y2(2)/Y2(1) < 0) then!It was converging before
                            !                            get_optimal_backtrack_par = 0.5 * (backtrack_pars(2) + backtrack_pars(3))!So we use previous convergence factors
                            get_optimal_backtrack_par = min(backtrack_pars(1) / 2, 0.5 * (backtrack_pars(2) + backtrack_pars(3)))
                        else!It was diverging as well before
                            get_optimal_backtrack_par = backtrack_pars(1) * 0.5!We halve the backtrack_par parameter to return to convergence fast
                        end if
                    end if
                end if
            end if
            !Make sure it is bounded
            get_optimal_backtrack_par = max(min(get_optimal_backtrack_par, 1.0), min_backtrack)

            !If we are stuck with the same values force a change
            if (abs(sum(backtrack_pars)/size(backtrack_pars)-get_optimal_backtrack_par)  &
                < 1d-5 .and. get_optimal_backtrack_par < 0.2 )  get_optimal_backtrack_par = max(min(2.0 * get_optimal_backtrack_par, 1.0), 1d-1)

        !            !Avoid exactly the same backtrack_par parameter as before
        !            if (abs(backtrack_pars(1)-get_optimal_backtrack_par) < 1d-3 ) get_optimal_backtrack_par = min(get_optimal_backtrack_par*1.01,1.0)
        end function

        subroutine Cubic_fitting(backtrack_pars, Convergences, Coefficients)
            implicit none
            real, dimension(:), intent(in) ::backtrack_pars, Convergences
            real, dimension(:), intent(inout) :: Coefficients
            !Local variables
            integer :: n, i, j, m
            real, dimension(:,:), allocatable :: A, A_inv

            !Check how much data we have
            do n = 1, size(backtrack_pars)
                if (backtrack_pars(n) < 0) exit
            end do
            n = n - 1
            !Not linear so far
            if (n<=2) return


            !Maximum degree of solution => cubic
            m = 3!<= quadratic
            if (n>=4) then!Two stored points plus the two extra points added to ensure minimums
                m = 4
            end if

            allocate(A(n,m), A_inv(m,m))
            !Construct matrix
            do i = 1, m
                do j = 1, n!Fill columns
                    A(j,i) = backtrack_pars(j)**real(m-i)
                    if (j==i) A(j,i) = A(j,i) + A(j,i)*(1d-7*real(rand(j)))
                end do
            end do

            !if square system
            if (n == m) then
                !Solve system
                call invert(A)!tested, coefficients are correct
                Coefficients = matmul(A, Convergences(1:m))
            else!Calculate curve fitting
                !coefficients = (A^t*A)^-1 * A^t * Convergences
                A_inv(1:m, 1:m) = matmul(transpose(A(1:n, 1:m)),A(1:n,1:m))!(M*n)*(n*m) => (m*m)

                call invert(A_inv(1:m, 1:m))
                Coefficients(1:m) = matmul(matmul(A_inv(1:m,1:m), transpose(A(1:n,1:m))), Convergences(1:n))
            end if
            deallocate(A, A_inv)
        end subroutine Cubic_fitting

    end subroutine FPI_backtracking

    subroutine Set_Saturation_to_sum_one(mdims, ndgln, state, packed_state)
        !This subroutines eliminates the oscillations in the saturation that are bigger than a
        !certain tolerance and also sets the saturation to be between bounds
        Implicit none
        !Global variables
        type( multi_dimensions ), intent( in ) :: Mdims
        type(multi_ndgln), intent(in) :: ndgln
        type( state_type ), dimension( : ), intent( in ) :: state
        type( state_type ), intent(inout) :: packed_state
        !Local variables
        type(scalar_field), pointer :: pipe_diameter
        integer :: iphase, cv_iloc, ele, cv_nod, i_start, i_end, ipres, stat
        real :: maxsat, minsat, correction, sum_of_phases, moveable_sat
        real, dimension(:), allocatable :: Normalized_sat
        real, dimension(:,:), pointer :: satura
        real, dimension(:, :), pointer :: Immobile_fraction


        call get_var_from_packed_state(packed_state, PhaseVolumeFraction = satura)
        !Get Immobile_fractions
        call get_var_from_packed_state(packed_state, Immobile_fraction = Immobile_fraction)

        if (Mdims%npres > 1) pipe_diameter => extract_scalar_field( state(1), "DiameterPipe" , stat = stat)
        !Allocate
        allocate(Normalized_sat(Mdims%nphase))
        !Impose sat to be between bounds for blocks of saturations (this is for multiple pressure, otherwise there is just one block)
        do ipres = 1, Mdims%npres
            i_start = 1 + (ipres-1) * Mdims%nphase/Mdims%npres
            i_end = ipres * Mdims%nphase/Mdims%npres
            !Set saturation to be between bounds (FOR BLACK-OIL maybe the limits have to be based on the previous saturation to allow
            !to have saturations below the immobile fractions, and the same for BoundedSolutionCorrection )
            do ele = 1, Mdims%totele
                do cv_iloc = 1, Mdims%cv_nloc
                    cv_nod = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                    if (ipres>1 .and. stat == 0) then
                        if (pipe_diameter%val(cv_nod) <=1d-8) cycle!Do not go out of the wells domain!!!
                    end if
                    moveable_sat = 1.0 - sum(Immobile_fraction(i_start:i_end, ele))
                    !Work in normalize saturation here
                    Normalized_sat(i_start:i_end) = (satura(i_start:i_end,cv_nod) - &
                        Immobile_fraction(i_start:i_end, ele))/moveable_sat
                    sum_of_phases = sum(Normalized_sat(i_start:i_end))
                    correction = (1.0 - sum_of_phases)
                    !Spread the error to all the phases weighted by their moveable presence in that CV
                    !Increase the range to look for solutions by allowing oscillations below 0.01 percent
                    if (abs(correction) > 1d-8) satura(i_start:i_end, cv_nod) = (Normalized_sat(i_start:i_end) * (1.0 + correction/sum_of_phases))*&
                        moveable_sat + Immobile_fraction(i_start:i_end, ele)
                    !Make sure saturation is between bounds after the modification
                    do iphase = i_start, i_end
                        minsat = Immobile_fraction(iphase, ele)
                        maxsat = moveable_sat + minsat
                        satura(iphase,cv_nod) =  min(max(minsat, satura(iphase,cv_nod)),maxsat)
                    end do
                end do
            end do
        end do
        !Deallocate
        deallocate(Normalized_sat)

    end subroutine Set_Saturation_to_sum_one

    subroutine Ensure_initial_Saturation_to_sum_one(mdims, ndgln, packed_state)
        !This subroutines eliminates the oscillations in the saturation that are bigger than a
        !certain tolerance and also sets the saturation to be between bounds
        Implicit none
        !Global variables
        type( multi_dimensions ), intent( in ) :: Mdims
        type(multi_ndgln), intent(in) :: ndgln
        type( state_type ), intent(inout) :: packed_state
        !Local variables
        integer :: iphase, cv_nod, cv_iloc, ele, i_start, i_end, ipres, scapegoat_phase
        real :: maxsat, minsat, sum_of_phases, moveable_sat
        real, dimension(:), allocatable :: Normalized_sat
        real, dimension(:,:), pointer :: satura
        real, dimension(:, :), pointer :: Immobile_fraction

        call get_var_from_packed_state(packed_state, PhaseVolumeFraction = satura)
        !Get Immobile_fractions
        call get_var_from_packed_state(packed_state, Immobile_fraction = Immobile_fraction)

        !Allocate
        allocate(Normalized_sat(Mdims%nphase))
        !Impose sat to be between bounds for blocks of saturations (this is for multiple pressure, otherwise there is just one block)
        do ipres = 1, Mdims%npres
            i_start = 1 + (ipres-1) * Mdims%nphase/Mdims%npres
            i_end = ipres * Mdims%nphase/Mdims%npres
            !First find a phase with values lower than -0.5 if that is so. This is a perturbation from opal and that needs to be considered
            !otherwise the last phase will be used to normalized
            scapegoat_phase = i_end
            do iphase = i_start, i_end
                if (maxval(satura(iphase,:))<-0.5) then
                    scapegoat_phase = iphase
                    exit
                end if
            end do
            !Once the saturation is found then ensure that the saturations are between bounds
            !Set saturation to be between bounds
            !to have saturations below the immobile fractions, and the same for BoundedSolutionCorrection )
            do ele = 1, Mdims%totele
                do cv_iloc = 1, Mdims%cv_nloc
                    cv_nod = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                    moveable_sat = 1.0 - sum(Immobile_fraction(i_start:i_end, ele))
                    !Work in normalize saturation here
                    Normalized_sat(i_start:i_end) = (satura(i_start:i_end,cv_nod) - &
                        Immobile_fraction(i_start:i_end, ele))/moveable_sat
                    sum_of_phases = sum(Normalized_sat(i_start:i_end))
                    !Ensure that the phases sum to 1.
                    satura(scapegoat_phase, cv_nod) = ((1.0 - sum_of_phases) + Normalized_sat(scapegoat_phase))*&
                        moveable_sat + Immobile_fraction(scapegoat_phase, ele)
                    !Make sure saturation is between bounds after the modification
                    do iphase = i_start, i_end
                        minsat = Immobile_fraction(iphase, ele)
                        maxsat = moveable_sat + minsat
                        satura(iphase,cv_nod) =  min(max(minsat, satura(iphase,cv_nod)),maxsat)
                    end do
                end do
            end do
        end do
        !Deallocate
        deallocate(Normalized_sat)

    end subroutine Ensure_initial_Saturation_to_sum_one


    subroutine auto_backtracking(Mdims, backtrack_par_factor, courant_number_in, first_time_step, nonlinear_iteration, I_am_temperature)
        !The maximum backtracking factor is calculated based on the Courant number and physical effects ocurring in the domain
        implicit none
        type(multi_dimensions), intent(in) :: Mdims
        real, intent(inout) :: backtrack_par_factor
        real, dimension(:), intent(inout) :: courant_number_in
        logical, intent(in) :: first_time_step
        integer, intent(in) :: nonlinear_iteration
        logical, optional, intent(in) :: I_am_temperature
        !Local variables
        real, save :: backup_shockfront_Courant = 0.
        real :: physics_adjustment, courant_number
        logical, save :: Readed_options = .false.
        logical, save :: gravity, cap_pressure, compositional, many_phases, black_oil, ov_relaxation, one_phase

        !Sometimes the shock-front courant number is not well calculated, then use previous value
        if (abs(courant_number_in(2)) < 1d-8 ) courant_number_in(2) = backup_shockfront_Courant
        !Combination of the overall and the shock-front Courant number
        courant_number = 0.4 * courant_number_in(1) + 0.6 * courant_number_in(2)
        backup_shockfront_Courant = courant_number_in(2)

        if (.not.readed_options) then
            !We read the options just once, and then they are stored as logicals
            gravity = have_option("/physical_parameters/gravity")
            if (have_option_for_any_phase("/multiphase_properties/capillary_pressure", Mdims%nphase)) then
                cap_pressure = .true.
            else
                cap_pressure = .false.
            end if
            compositional = Mdims%ncomp > 0
            many_phases = Mdims%n_in_pres > 2
            black_oil = have_option( "/physical_parameters/black-oil_PVT_table")
            !Positive effects on the convergence !Need to check for shock fronts...
            ov_relaxation = have_option('/timestepping/nonlinear_iterations/Fixed_Point_Iteration/Vanishing_relaxation')

            one_phase = (Mdims%n_in_pres == 1)
            readed_options = .true.
        end if

        !Depending on the active physics, the problem is more complex and requires more relaxation
        physics_adjustment = 1.
        !Negative effects on the convergence
        if (gravity) physics_adjustment = physics_adjustment * 5.0
        if (cap_pressure) physics_adjustment = physics_adjustment * 5.0
        if (compositional) physics_adjustment = physics_adjustment * 1.5
        if (many_phases) physics_adjustment = physics_adjustment * 1.5
        !For the first two non-linear iterations, it has to re-adjust, as the gas and oil are again mixed
        if (black_oil .and. nonlinear_iteration <= 2) physics_adjustment = physics_adjustment * 2.

        !Positive effects on the convergence !Need to check for shock fronts...
        if (ov_relaxation) physics_adjustment = physics_adjustment * 0.05!huge benefits when using ov_relaxation...
        if (one_phase) physics_adjustment = physics_adjustment * 0.5
        if (present_and_true(I_am_temperature)) then
            !Much simpler for temperature?
            if (.not.ov_relaxation) physics_adjustment = physics_adjustment * 0.01
!            if (Mdims%npres>1) physics_adjustment = physics_adjustment * 5.0
        end if
        !For the time being, it is based on this simple table
        if (Courant_number * physics_adjustment > 50.) then
            backtrack_par_factor = -0.1
        else if (Courant_number * physics_adjustment > 25.) then
            backtrack_par_factor = -0.15
        else if (Courant_number * physics_adjustment > 15.) then
            backtrack_par_factor = -0.2
        else if (Courant_number * physics_adjustment > 8.) then
            backtrack_par_factor = -0.33
        else if (Courant_number * physics_adjustment > 5.) then
            backtrack_par_factor = -0.5
        else if (Courant_number * physics_adjustment > 1) then
            backtrack_par_factor = -0.8
        else
            backtrack_par_factor = -1.
        end if
        !For the first calculation, the Courant number is usually zero, hence we force a safe value here
        if (first_time_step .and. nonlinear_iteration == 1) backtrack_par_factor = -0.05
        !Use the most restrictive value across all the processors
        if (IsParallel()) call allmin(backtrack_par_factor)

    end subroutine auto_backtracking

end module solvers_module


