!    Copyright (C) 2020 Imperial College London and others.
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Affero General Public License
!    as published by the Free Software Foundation,
!    version 3.0 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without seven the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"

module solvers_module

    use fldebug
    use fields
    use Petsc_tools
    use parallel_tools
    use sparse_tools_petsc
    use solvers
    use global_parameters, only: OPTION_PATH_LEN, FPI_have_converged
    use spud
    use parallel_tools, only : allmax, allmin, isparallel, getprocno
    use parallel_fields

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

    public :: BoundedSolutionCorrections, FPI_backtracking, Set_Saturation_to_sum_one,&
         Initialise_Saturation_sums_one, auto_backtracking, get_Anderson_acceleration_new_guess, &
         non_porous_ensure_sum_to_one, duplicate_petsc_matrix, scale_PETSc_matrix, petsc_Stokes_solver


contains


  !!>@brief: This subroutine adjusts field_val so that it is bounded between field_min, field_max in a local way.
  !> The sparcity of the local CV connectivity is in: small_findrm, small_colm.
  !> ngl_its=max no of global iterations e.g. 100.
  !> error_tol = tolerance on the iterations.
  !>
  !> nloc_its: This iteration is very good at avoiding spreading the modifications too far - however it can stagnate.
  !> nloc_its2: This iteration is very good at avoiding stagnating but does spread the modifcations far.
  !> us a single iteration because of this as default...
  !> nits_nod: iterations at a nod - this iteration is very good at avoiding spreading the modifications too far -
  !> however it can stagnate.
    subroutine BoundedSolutionCorrections( state, packed_state, &
        Mdims, CV_funs, small_findrm, small_colm, Field_name, &
        for_sat, min_max_limits)
        implicit none
        integer, parameter :: nloc_its = 5, nloc_its2 = 1, nits_nod = 100, ngl_its = 500
        real, parameter :: w_relax = 0.5, error_tol = 1.0e-5
        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_shape_funs), intent(in) :: CV_funs
        integer, dimension( : ), intent( in ) :: small_findrm, small_colm
        character( len=* ), intent(in) :: Field_name
        real, optional, dimension(2) :: min_max_limits
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
        type( vector_field ), pointer :: x
        real, dimension(:,:), pointer :: Immobile_fraction
        if (present_and_true(for_sat)) then
            field => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
        else
            field => extract_tensor_field( packed_state, trim(Field_name) )

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
            !Since everything is local, we overcycle and then we can avoid halo_updates
            ! if ( .not. node_owned( field, inod ) ) cycle
            do count = small_findrm( inod ), small_findrm( inod + 1 ) - 1
                jnod = small_colm( count )
                mass_cv_sur(inod) = mass_cv_sur(inod) + mass_cv( jnod )
            end do
        end do
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
        if (present(min_max_limits)) then
            !Over-write limits with the passed down values
           field_min = min_max_limits(1)
           field_max = min_max_limits(2)
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
            !We are potentially calling halo_update 500 times...
            ! call halo_update( field )!For principles this should not even exist here in this loop

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
                !We are potentially calling halo_update 500 times...
                ! call halo_update( field )!For principles this should not even exist here in this loop
            end do ! loc_its2
            ! communicate the errors ( max_change, error_changed ) ...
            ! this could be more efficient sending a vector...
            max_max_error= max( max_change, error_changed )
            call allmax( max_max_error )
            if ( max_max_error < error_tol ) exit
        end do ! gl_its
        !After performing everything update halos only once...
        if (IsParallel()) call halo_update( field )

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
        call deallocate_multi_dev_shape_funs(DevFuns)
        return
    end subroutine BoundedSolutionCorrections

    !sprint_to_do!not use one global variable
    !!>@brief:In this subroutine we applied some corrections and backtrack_par on the saturations obtained from the saturation equation
    !>this idea is based on the paper SPE-173267-MS.
    !>The method ensures convergence "independent" of the time step.
    subroutine FPI_backtracking(nphase, Mdims, ndgln, state, packed_state, sat_bak, backtrack_sat, backtrack_par_from_schema, &
        Previous_convergence, satisfactory_convergence, new_backtrack_par, Max_sat_its, its, nonlinear_iteration, useful_sats, res, &
        res_ratio, first_res)
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
        integer, intent(in) :: Max_sat_its, its, nonlinear_iteration, nphase
        integer, intent(inout) :: useful_sats
        !Local parameters
        real, parameter :: Conv_to_achiv = 10.0
        real, save :: anders_exp!This parameter change the importance of backtrack_sat in Anderson's acceleration (mainly for high alphas)
        !Local variables        !100 => backtrack_sat is not used; 0.3 => equally important; 0.4 => recommended; 0 => more important than sat_bak
        real, dimension(:, :), pointer :: Satura
        logical :: new_time_step, new_FPI, Undo_update
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
          call Set_Saturation_to_sum_one(mdims, ndgln, packed_state, state, do_not_update_halos = .TRUE. )
        else
          call non_porous_ensure_sum_to_one(mdims, packed_state, do_not_update_halos = .TRUE.)
        end if
        sat_field => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
        Satura =>  sat_field%val(1,:,:)
        !Stablish minimum backtracking parameter
        min_backtrack = 0.1

        !Automatic method based on the history of convergence
        if (backtrack_par_from_schema < 0.0) then

            !Retrieve convergence factor, to make sure that if between time steps things are going great, we do not reduce the
            !backtrack_par_parameter
            if (backtrack_pars(1) < 0) then!retrieve it just once
                call get_option( '/solver_options/Non_Linear_Solver/Fixed_Point_Iteration',&
                    convergence_tol, default = 0. )
                convergence_tol =  convergence_tol
                !Tolerance for the infinite norm
                call get_option( '/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Infinite_norm_tol',&
                    Infinite_norm_tol, default = 0.03 )
                backtrack_pars(1) = max(min(abs(backtrack_par_from_schema), 1.0), min_backtrack)
                !Retrieve the shape of the function to use to weight the importance of previous saturations
                call get_option( '/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Acceleration_exp',&
                    anders_exp, default = 0.4 )!0.4 the best option, based on experience
                !Should not be negative
                anders_exp = max(anders_exp, 0.)
            end if

            if (new_time_step) then
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
                    maxval(abs(Sat_bak-Satura(1:nphase,:)))/backtrack_pars(2) < Infinite_norm_tol)!<= exit if final convergence is achieved
                if (IsParallel()) call alland(satisfactory_convergence)
                !If a backtrack_par parameter turns out not to be useful, then undo that iteration
                Undo_update = its > 2 .and. Convergences(2) > 0 .and. allow_undo .and. Convergences(1)>5.
                if (IsParallel()) call allor(Undo_update)!Consistently repeat an update if required
                if (Undo_update) then
                    Satura(1:nphase,:) = backtrack_sat
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
                !If it is parallel then we want to be consistent between cpus
                !we use the smallest value, since it is more conservative
                if (IsParallel()) then
                    call allmin(backtrack_pars(1))
                    call allmin(Convergences(1))
                end if
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


        !***Calculate new saturation***
        !Obtain new saturation using the backtracking method
        if (useful_sats < 2 .or. satisfactory_convergence) then
            !Since Anderson's acceleration is unstable, when it has converged, we use the stable form of backtracking
            Satura(1:nphase,:) = sat_bak * (1.0 - backtrack_pars(1)) + backtrack_pars(1) * Satura(1:nphase,:)
        else !Use Anderson acceleration, idea from "AN ACCELERATED FIXED-POINT ITERATION FOR SOLUTION OF VARIABLY SATURATED FLOW"

            !Based on making backtrack_sat small when backtrack_pars(1) is high and backtrack_sat small when backtrack_pars(1) is small
            !The highest value of backtrack_sat is displaced to low values of alpha
            aux = 1.0 - backtrack_pars(1)
            Satura(1:nphase,:) = backtrack_pars(1) * Satura(1:nphase,:) + aux * ( (1.-(aux**anders_exp *backtrack_pars(1)) ) * sat_bak + &
                aux**anders_exp *backtrack_pars(1) * backtrack_sat)!<=The best option so far

        end if
        !Re-impose physical constraints
        if (is_porous_media) then
          call Set_Saturation_to_sum_one(mdims, ndgln, packed_state, state, do_not_update_halos = .FALSE. )
        else
          call non_porous_ensure_sum_to_one(Mdims, packed_state, do_not_update_halos = .FALSE.)
        end if

        !Inform of the new backtrack_par parameter used
        new_backtrack_par = backtrack_pars(1)

    contains
        !!>@brief: Based on a history of convergence and backtracking factors an optimal bracktrack factor is computed
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
                    else !Diverging, the reduce with a minimum value that will mean performing all the non-linear iterations
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

                            else !Old slope was steeper => Optimal value in between the previous backtrack_par parameters
                                get_optimal_backtrack_par = min(backtrack_pars(1) / 2, 0.5 * (backtrack_pars(2) + backtrack_pars(3)))

                            end if
                        else !It started to converge now, so we encourage to get away from the bad convergence value
                            get_optimal_backtrack_par = backtrack_pars(1) * max(1.2, 2.0 * X(1))
                        end if
                    else !It is NOT converging now
                        if (Y2(2)/Y2(1) < 0) then!It was converging before
                            !                            get_optimal_backtrack_par = 0.5 * (backtrack_pars(2) + backtrack_pars(3))!So we use previous convergence factors
                            get_optimal_backtrack_par = min(backtrack_pars(1) / 2, 0.5 * (backtrack_pars(2) + backtrack_pars(3)))
                        else !It was diverging as well before
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

        !!>@brief: Fitting of two or three points, it solves an easy least squares system
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
            else !Calculate curve fitting
                !coefficients = (A^t*A)^-1 * A^t * Convergences
                A_inv(1:m, 1:m) = matmul(transpose(A(1:n, 1:m)),A(1:n,1:m))!(M*n)*(n*m) => (m*m)

                call invert(A_inv(1:m, 1:m))
                Coefficients(1:m) = matmul(matmul(A_inv(1:m,1:m), transpose(A(1:n,1:m))), Convergences(1:n))
            end if
            deallocate(A, A_inv)
        end subroutine Cubic_fitting

    end subroutine FPI_backtracking

    !!>@brief:This subroutines eliminates the oscillations in the saturation that are bigger than a
    !>certain tolerance and also sets the saturation to be between bounds
    subroutine Set_Saturation_to_sum_one(mdims, ndgln, packed_state, state, do_not_update_halos)
        Implicit none
        !Global variables
        type( multi_dimensions ), intent( in ) :: Mdims
        type(multi_ndgln), intent(in) :: ndgln
        type( state_type ), intent(inout) :: packed_state
        type( state_type ), dimension(:), intent(in) :: state
        logical, optional, intent(in) :: do_not_update_halos
        !Local variables
        type(scalar_field), pointer :: pipe_diameter
        type(tensor_field), pointer :: sat_field, old_saturation_field
        integer :: iphase, cv_iloc, ele, cv_nod, i_start, i_end, ipres, stat, k
        real :: maxsat, minsat, correction, sum_of_phases, moveable_sat
        real, dimension(Mdims%nphase) :: Normalized_sat
        real, dimension(:,:), pointer :: satura
        real, dimension(:, :), pointer :: Immobile_fraction
        logical, save :: Solve_all_phases = .true., first_time = .true.

        !Obtain saturation field from packed_state
        sat_field => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
        satura =>  sat_field%val(1,:,:)
        !Get Immobile_fractions
        call get_var_from_packed_state(packed_state, Immobile_fraction = Immobile_fraction)

        if ( first_time) then
          Solve_all_phases = .not. have_option("/numerical_methods/solve_nphases_minus_one")
          first_time = .false.
        end if

        if (Mdims%npres > 1) pipe_diameter => extract_scalar_field( state(1), "DiameterPipe" , stat = stat)
        !Impose sat to be between bounds for blocks of saturations (this is for multiple pressure, otherwise there is just one block)
        if (Solve_all_phases) then
          do ipres = 1, Mdims%npres
              i_start = 1 + (ipres-1) * Mdims%nphase/Mdims%npres
              i_end = ipres * Mdims%nphase/Mdims%npres
              !Set saturation to be between bounds (FOR BLACK-OIL maybe the limits have to be based on the previous saturation to allow
              !to have saturations below the immobile fractions, and the same for BoundedSolutionCorrection )
              do ele = 1, Mdims%totele
                  do cv_iloc = 1, Mdims%cv_nloc
                      cv_nod = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                      if ( .not. node_owned( sat_field, cv_nod ) ) cycle
                      !Do not go out of the wells domain!!!
                      if (ipres>1) then
                        if(pipe_diameter%val(cv_nod)<=1d-8) cycle
                      end if                      !Do not go out of the wells domain!!!
                      moveable_sat = 1.0 - sum(Immobile_fraction(i_start:i_end, ele))
                      !Work in normalized saturation here
                      Normalized_sat(i_start:i_end) = (satura(i_start:i_end,cv_nod) - &
                          Immobile_fraction(i_start:i_end, ele))/moveable_sat
                      sum_of_phases = sum(Normalized_sat(i_start:i_end))
                      correction = (1.0 - sum_of_phases)
                      !Spread the error to all the phases weighted by their moveable presence in that CV
                      !Increase the range to look for solutions by allowing oscillations below 0.01 percent
                      if (abs(correction) > 1d-8) satura(i_start:i_end, cv_nod) = (Normalized_sat(i_start:i_end) * &
                          (1.0 + correction/sum_of_phases))* moveable_sat + Immobile_fraction(i_start:i_end, ele)
                      !Make sure saturation is between bounds after the modification
                      do iphase = i_start, i_end
                          minsat = Immobile_fraction(iphase, ele)
                          maxsat = moveable_sat + minsat
                          satura(iphase,cv_nod) =  min(max(minsat, satura(iphase,cv_nod)),maxsat)
                      end do
                  end do
              end do
          end do
        else !Solving for nphases -1 requires first to limit the saturation between bounds and then impose sum phases = 1
          do ipres = 1, Mdims%npres
              i_start = 1 + (ipres-1) * Mdims%n_in_pres
              i_end = ipres * Mdims%n_in_pres
              !Set saturation to be between bounds
              !to have saturations below the immobile fractions, and the same for BoundedSolutionCorrection )
              do ele = 1, Mdims%totele
                  do cv_iloc = 1, Mdims%cv_nloc
                      cv_nod = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                      if ( .not. node_owned( sat_field, cv_nod ) ) cycle
                      !Do not go out of the wells domain!!!
                      if (ipres>1) then
                        if(pipe_diameter%val(cv_nod)<=1d-8) cycle
                      end if
                      !Make sure saturation is between bounds before the modification
                      do iphase = i_start, i_end
                          minsat = Immobile_fraction(iphase, ele)
                          maxsat = moveable_sat + minsat
                          satura(iphase,cv_nod) =  min(max(minsat, satura(iphase,cv_nod)),maxsat)
                      end do

                      moveable_sat = 1.0 - sum(Immobile_fraction(i_start:i_end, ele))
                      !Work in normalize saturation here
                      Normalized_sat(i_start:i_end) = (satura(i_start:i_end,cv_nod) - &
                          Immobile_fraction(i_start:i_end, ele))/moveable_sat
                      sum_of_phases = sum(Normalized_sat(i_start:i_end))
                      !Ensure that the phases sum to 1.
                      satura(i_end, cv_nod) = ((1.0 - sum_of_phases) + Normalized_sat(i_end))*&
                          moveable_sat + Immobile_fraction(i_end, ele)
                  end do
              end do
          end do
        end if

        if (present_and_true(do_not_update_halos)) return
        !Ensure cosistency across CPUs
        if (IsParallel())call halo_update(sat_field)

    end subroutine Set_Saturation_to_sum_one


    !!>@brief: This subroutines eliminates the oscillations in the saturation that are bigger than a
    !> certain tolerance and also sets the saturation to be between bounds
     subroutine non_porous_ensure_sum_to_one(Mdims, packed_state, do_not_update_halos)
         Implicit none
         !Global variables
         type( state_type ), intent(inout) :: packed_state
         type( multi_dimensions ), intent( in ) :: Mdims
         logical, optional, intent(in) :: do_not_update_halos
         !Local variables
         integer :: iphase, cv_nod, i_start, i_end, ipres
         real :: correction, sum_of_phases
         real, dimension(:,:), pointer :: satura
         type(tensor_field), pointer :: tfield

         tfield => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
         satura =>  tfield%val(1,:,:)

         !Impose sat to be between bounds for blocks of saturations (this is for multiple pressure, otherwise there is just one block)
         do ipres = 1, Mdims%npres
             i_start = 1 + (ipres-1) * Mdims%nphase/Mdims%npres
             i_end = ipres * Mdims%nphase/Mdims%npres
             !Set saturation to be between bounds
             do cv_nod = 1, size(satura,2 )
                 sum_of_phases = sum(satura(i_start:i_end, cv_nod))
                 correction = (1.0 - sum_of_phases)
                 !Spread the error to all the phases weighted by their presence in that CV
                 !Increase the range to look for solutions by allowing oscillations below 0.1 percent
                 if (abs(correction) > 1d-3) satura(i_start:i_end-1, cv_nod) = (satura(i_start:i_end-1, cv_nod) * (1.0 + correction/sum_of_phases))
     !if (abs(correction) > 1d-3) satura(i_start:i_end, cv_nod) = (satura(i_start:i_end, cv_nod) * (1.0 + correction/sum_of_phases))
                 !Make sure saturation is between bounds after the modification
                 do iphase = i_start, i_end
                     satura(iphase,cv_nod) =  min(max(0., satura(iphase,cv_nod)),1.0)
                 end do
             end do
         end do

         if (present_and_true(do_not_update_halos)) return
         !Ensure cosistency across CPUs
         if (IsParallel())call halo_update(tfield)

     end subroutine non_porous_ensure_sum_to_one

    !!>@brief:Ensure that the saturations at the beginning sum to one, if they do not
    !> all the error is compensated in the scapegoat_phase. Normally the last
    subroutine Initialise_Saturation_sums_one(mdims, ndgln, packed_state, find_scapegoat_phase)
        Implicit none
        !Global variables
        type( multi_dimensions ), intent( in ) :: Mdims
        type(multi_ndgln), intent(in) :: ndgln
        type( state_type ), intent(inout) :: packed_state
        logical, optional, intent(in) :: find_scapegoat_phase
        !Local variables
        integer :: iphase, cv_nod, cv_iloc, ele, i_start, i_end, ipres, scapegoat_phase
        real :: maxsat, minsat, sum_of_phases, moveable_sat
        real, dimension(:), allocatable :: Normalized_sat
        real, dimension(:,:), pointer :: satura
        type(tensor_field), pointer :: sat_field
        real, dimension(:, :), pointer :: Immobile_fraction

        !Obtain saturation field from packed_state
        sat_field => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
        satura =>  sat_field%val(1,:,:)        !Get Immobile_fractions
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
            if (present_and_true(find_scapegoat_phase)) then
              do iphase = i_start, i_end
                  if (maxval(satura(iphase,:))<-0.5) then
                      scapegoat_phase = iphase
                      exit
                  end if
              end do
            end if
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
        !Update halos
        call halo_update(sat_field)!Ensure consistency across CPUs
        !Deallocate
        deallocate(Normalized_sat)

    end subroutine Initialise_Saturation_sums_one


    !!>@brief: The maximum backtracking factor is calculated based on the Courant number and physical effects ocurring in the domain
    subroutine auto_backtracking(Mdims, backtrack_par_factor, courant_number_in, first_time_step, nonlinear_iteration, I_am_temperature)
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
            ov_relaxation = have_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation')

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
        if (IsParallel()) call allmin(backtrack_par_factor)!Should not be necessary as the Courant numbers are already parallel safe

    end subroutine auto_backtracking




    !---------------------------------------------------------------------------
    !> @author Pablo Salinas
    !> @brief This subroutine provides a new guess based on previous updates and residuals.
    !> Use this subroutine to speed up any system being solved by looping.
    !> A new guess is computed and returned
    !> Method explained in DOI.10.1137/10078356X
    !---------------------------------------------------------------------------
    subroutine get_Anderson_acceleration_new_guess(N, M, NewField, History_field, stored_residuals, max_its, prev_small_matrix,restart_now)
      implicit none
      integer, intent(in) :: N !> Size if the field of interest. Used also to turn vector/tensor fields into scalar fields internally here
      integer, intent(in) :: M !> Size of the field of iterations available - 2
      real, dimension(N), intent(out) :: NewField !> New guess obtained by the Anderson Acceleration method
      real, dimension(N, M+2), intent(inout) :: History_field !> Past results obtained by the outer solver
      real, dimension(N, M+1), intent(inout) :: stored_residuals !> Past residuals obtained as the (- guessed value + G(guess value) )
      integer, optional :: max_its
      real, dimension(:,:), allocatable, optional, intent(inout) :: prev_small_matrix!> Contains the previous A'*A, speeds up the method. On exit is updated (M,M)
                                                                        !> Prepared to be passed unallocated from M = 1 and allocated internally
      logical, intent(inout) :: restart_now
      !Local variables
      real, dimension(N,M) :: Matrix !> Matrix containing the residuals
      real, dimension(M,1) :: Small_b !>RHS containing A' * The last residual
      real, dimension(M) :: auxV !>Temporary array to perform in parallel A'*A
      real, dimension(M,M) :: Small_matrix !>The resulting from A'*A
      real, dimension(m+1) :: AA_alphas
      integer :: i, j, Q_rank, ierr, max_its2, start
      real :: auxR
      logical :: present_and_useful_bak_mat

      !Not only being present prev_small_matrix we can just use it, needs to be allocated and have some information!
      present_and_useful_bak_mat = .false.
      if (present(prev_small_matrix)) then
        if (allocated(prev_small_matrix)) then
          present_and_useful_bak_mat =  size(prev_small_matrix,1)>1 .and. size(prev_small_matrix,1) == M - 1
        end if
      end if


      max_its2 = 20
      if (present(max_its)) max_its2 = max_its

      !Need at least m to be 1
      if (m <= 0) then
        !Return the last field to proceed as a normal iteration
        NewField = History_field(:,m+2)
        return
      end if

      !Construct matrix
      do i = 1, m !Here we form the matrix by making the difference of deltaF which is define as the variation of the field between updates
        do j = 1, n
          Matrix(j,i) = (stored_residuals(j, i+1) - stored_residuals(j, i))
        end do
      end do
      !For parallel it is easier to solve the system A'*A = A'*b because the system is tiny so each processor can solve it separately
      start = 1
      !If stored matrix, we re-use that matrix
      if (present_and_useful_bak_mat) then
        Small_matrix(1:M-1, 1:M-1) = prev_small_matrix
        start = 1!M
      end if
      !Perform now A'*A
      do i = 1, M
        do j = max(i,start), M!Only the upper part, we will copy this in serial to minimise parallel communications
          Small_matrix(j,i) = dot_product(Matrix(:,j), Matrix(:,i))
        end do
        !Now add up all the column between processors
        call allsum(Small_matrix(:,i))
      end do
      !Copy now the symmetric part
      do i = 1, M
        do j = i+1, M!Not the diagonal
          Small_matrix(i,j) = Small_matrix(j,i)
        end do
      end do

      !Now the RHS
      do i = 1, M
        Small_b(i,1) = dot_product(Matrix(:,i), stored_residuals(:,m+1))
      end do
      call allsum(Small_b(:,1))

      !Now solve the least squares optimisation problem
      Q_rank = m
      if (M > 1) call Least_squares_solver(Small_matrix,Small_b, Q_rank)

      if (Q_rank < m) then
        !If the least squares matrix is not full rank, then perform normal update
        !and restart the cycle
        AA_alphas = 0; AA_alphas(size(AA_alphas)) = 1
        restart_now = .true.
      else if (M == 1) then
        AA_alphas = 0; AA_alphas(size(AA_alphas)) = 1
      else
        !Now we proceed to convert to alphas so we can compute the new guess
        !Here we calculate Gammas instead of alphas (size m-1 instead of m)
        AA_alphas(size(AA_alphas)) = 1 - Small_b(m,1)
        do i = m, 2, -1 !The last one is different, the first one is the same
          AA_alphas(i) = Small_b(i,1) - Small_b(i-1,1)
        end do
        AA_alphas(1) = Small_b(1,1)
      end if
        !Multiply previous fields by the alphas to obtain the new guess
      NewField = 0.
      do i = 1, size(AA_alphas)
        NewField = NewField + History_field(:,i+1) * AA_alphas(i)
      end do

      if (present(prev_small_matrix)) then
        !Time to update prev_small_matrix!
        if (allocated(prev_small_matrix))deallocate(prev_small_matrix)
        allocate(prev_small_matrix(M,M))
        !Backup
        prev_small_matrix = Small_matrix
      end if

    end subroutine get_Anderson_acceleration_new_guess

    !---------------------------------------------------------------------------
    !> @author Pablo Salinas
    !> @brief In this subroutine the matrix is re-scaled based on the formula
    !> D^-0.5 * A * D^-0.5 X'=  D^-0.5 b; and next X = D^-0.5 * X';
    !> IMPORTANT: the step X = D^-0.5 * X' needs to be done elsewhere store the diagonal before calling this
    !> This should allow to deal with high ranges of viscosity ratio for example
    !> A is-written
    !---------------------------------------------------------------------------
    subroutine scale_PETSc_matrix(Mat_petsc)
      implicit none
      type(petsc_csr_matrix), intent(inout)::  Mat_petsc !>  System matrix in PETSc format
      !Local variables
      integer :: ierr, m, n
      Vec, target :: scale_diag

      !Proceed to allocate memory for the diagonal of A
      call MatGetLocalSize(Mat_petsc%M,m, n, ierr)
      if (isparallel()) then
        call VecCreateMPI(MPI_COMM_FEMTOOLS, m, PETSC_DETERMINE, scale_diag, ierr)
      else
        call VecCreateSeq(MPI_COMM_SELF, m, scale_diag, ierr)
      end if
      !Need the matrix ssembled
      call MatAssemblyBegin(Mat_petsc%M, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(Mat_petsc%M, MAT_FINAL_ASSEMBLY, ierr)
      !Extract diagonal from A
      call MatGetDiagonal(Mat_petsc%M, scale_diag, ierr)
      !Compute sqrt (we do this first to reduce the span induced by high viscosity ratios)
      call VecSqrtAbs(scale_diag, ierr)
      !Compute the inverse
      call VecReciprocal(scale_diag, ierr)
      !Proceed to re-scale the matrix by doing D^-0.5 * Mat_petsc * D^-0.5
      call MatDiagonalScale(Mat_petsc%M, scale_diag, scale_diag, ierr)
      !Deallocate unnecessary memory
      call VecDestroy(scale_diag, ierr)

    end subroutine scale_PETSc_matrix

    subroutine duplicate_petsc_matrix(MAT_A,MAT_B)
      type(petsc_csr_matrix), intent(in)::MAT_A
      type(petsc_csr_matrix), intent(inout)::MAT_B
      !Local variables
      integer :: ierr

      call allocate(MAT_B, MAT_A%M, MAT_A%row_numbering, MAT_A%column_numbering, "DGM_PETSC_scaled")
      call MatDuplicate(MAT_A%M,MAT_COPY_VALUES,MAT_B%M, ierr)!Deep copy
      call assemble(MAT_B)

    end subroutine

    !---------------------------------------------------------------------------
    !> @author Pablo Salinas
    !> @brief In this subroutine the Schur complement is generated and solved using PETSc to update the pressure field
    !> Matrices need to be in petsc format and pmat is the preconditioned matrix, i.e. pmat = A11- A10(AproxA00^-1)A01
    !---------------------------------------------------------------------------
    subroutine petsc_Stokes_solver(packed_state, Mdims, Mmat, ndgln, Mspars, final_phase, pmat, P_all, deltaP, rhs_p, solver_option_path, Dmat)
        use Full_Projection
        use petsc_tools
        use solvers
        
        IMPLICIT NONE
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_matrices), intent(inout) :: Mmat
        type(petsc_csr_matrix), intent( inout )::  pmat
        type( vector_field ), intent(in) :: rhs_p
        type( vector_field ), INTENT(INOUT) :: deltap
        type( tensor_field ), INTENT(INOUT) :: P_all
        character(len=option_path_len), intent(in) :: solver_option_path
        type(multi_ndgln), intent(in) :: ndgln
        type(multi_sparsities), intent(in) :: Mspars
        INTEGER, intent( in ) :: final_phase
        type(petsc_csr_matrix), optional, intent( inout )::  Dmat
        !Local variables
        type(petsc_numbering_type) :: petsc_numbering, petsc_numbering_vel
        integer :: ierr, literations
        type (scalar_field) :: scalar_rhs_p, scalar_delta_P
        type( csr_sparsity ), pointer :: sparsity
        logical :: has_diffusion_operator
        ! type(petsc_csr_matrix)::Schur_mat
        Vec y, b
        Mat S, Schur_mat
        KSP ksp_schur

        !Convert the matrices from ICFERST to PETSC format (not ideal...)
        call Convert_C_and_CT_mat_to_PETSc_format(packed_state, Mdims, Mmat, ndgln, Mspars, final_phase) 

        !generate sparsity so we can communicate with PETSc
        sparsity=>extract_csr_sparsity(packed_state,"CMatrixSparsity")
        call allocate(petsc_numbering, node_count(deltaP), deltaP%dim, &
                            halo=sparsity%row_halo)

        !Scale Mmat%CT_PETSC%M to follow the sign convention of PETSc
        !SPRINT_TO_DO, IF EVERYTHING WORKS, THEN PROBABLY BETTER FLIP THE SIGN OF THE UPDATE ONLY??
        call MatScale(Mmat%CT_PETSC%M,real(-1.0, kind = PetscScalar_kind),ierr)

        if(present(Dmat)) then
            call MatCreateSchurComplement(Mmat%DGM_PETSC%M, Mmat%DGM_PETSC%M, Mmat%C_PETSC%M, Mmat%CT_PETSC%M, Dmat%M, Schur_mat, ierr)
        else
    ! workaround bug in petsc 3.8: missing CHKFORTRANNULLOBJECT, so PETSC_NULL_MAT isn't translated to C null
    ! this however seems to do the trick:
#if PETSC_VERSION_MINOR>=8
    S%v = 0
#else
    S = PETSC_NULL_OBJECT
#endif
            !Because we provide our own pmat = A11- A10(AproxA00^-1)A01 we just populate here with A00 the approx of A00
            call MatCreateSchurComplement(Mmat%DGM_PETSC%M, Mmat%DGM_PETSC%M, Mmat%C_PETSC%M, Mmat%CT_PETSC%M, S, Schur_mat, ierr)
        end if
        !This is to obtain using PETSc the pmat matrix, currently not working
        !MAT_SCHUR_COMPLEMENT_AINV_DIAG, MAT_SCHUR_COMPLEMENT_AINV_LUMP, or MAT_SCHUR_COMPLEMENT_AINV_BLOCK_DIAG
        ! call MatSchurComplementSetAinvType(Schur_mat, MAT_SCHUR_COMPLEMENT_AINV_DIAG, ierr)
        ! call MatSchurComplementGetPmat(Schur_mat, MAT_INITIAL_MATRIX, pmat%M, ierr)

        ! Create KSP (solver_option_path should point to pressure options or the default)
        call attach_null_space_from_options(Schur_mat, solver_option_path, pmat=pmat%M, petsc_numbering=petsc_numbering)
        call create_ksp_from_options(ksp_schur,Schur_mat,pmat%M, solver_option_path, isparallel(),petsc_numbering, .true.)

        ! Create RHS and Solution vectors in PETSc format
        b = PetscNumberingCreateVec(petsc_numbering); call VecDuplicate(b,y,ierr)
        call field2petsc(rhs_p, petsc_numbering, b)
        call field2petsc(deltaP, petsc_numbering, y)
  
        ! Solve Schur complement system
        call petsc_solve_core(y, Schur_mat, b, ksp_schur, petsc_numbering, solver_option_path, .true., &
            literations, nomatrixdump=.true., vector_x0 = deltaP, vfield = deltaP)
! print *, "Iterations taken", literations
        !Copy value back
        call petsc2field(y, petsc_numbering, deltap)
        ! Destroy all PETSc related variables
        call petsc_solve_destroy(y, Schur_mat, b, ksp_schur, petsc_numbering, solver_option_path)

        !Flip sign here if we haven't flipped the sign of the Ct matrix before
        ! deltaP%val = -1.* deltaP%val
        !Update pressure now
        P_all%val(1,:,:) = P_all%val(1,:,:) + deltaP%val
        if (isParallel()) call halo_update(P_all)

    end subroutine

    !>@brief: This subroutine converts the C (gradient) and CT (Divergence) matrices into PETSc format
    !> Mainly devoted to be used by the PETSc stokes schur solver
    SUBROUTINE Convert_C_and_CT_mat_to_PETSc_format(packed_state, Mdims, Mmat, ndgln, &
        Mspars, NPHASE)  ! Element connectivity.
        IMPLICIT NONE
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_matrices), intent(inout) :: Mmat
        type(multi_ndgln), intent(in) :: ndgln
        type(multi_sparsities), intent(in) :: Mspars
        INTEGER, intent( in ) :: NPHASE
        !Local variables
        type( csr_sparsity ), pointer :: sparsity
        type( tensor_field ), pointer :: velocity, pressure
        INTEGER :: IU_NOD, count, IPHASE, IDIM, P_JNOD, j, ele
        integer, dimension(:), pointer :: neighbours
        integer :: nb
        logical :: skip

        !Extract a vector and a scalar fields for parallel checks
        Velocity => extract_tensor_field( packed_state, "PackedVelocity" )
        Pressure => extract_tensor_field( packed_state, "PackedFEPressure" )

        !C does not need to be recomputed every time-level
        if (.not. Mmat%stored) then 
            !TEMPORARY, Mmat%C_PETSC DOES NOT NEED TO BE REDONE UNLESS THE MESH CHANGES
            if (associated(Mmat%C_PETSC%refcount)) call deallocate(Mmat%C_PETSC)
            sparsity=>extract_csr_sparsity(packed_state,"CMatrixSparsity")
            call allocate(Mmat%C_PETSC,sparsity,[Mdims%ndim,NPHASE],name="C_PETSC"); call zero(Mmat%C_PETSC)

            DO IU_NOD = 1, Mdims%U_NONODS
                DO COUNT = mspars%C%fin( IU_NOD ), mspars%C%fin( IU_NOD + 1 ) - 1
                    P_JNOD = mspars%C%col( COUNT )
                    do idim =1, Mdims%ndim
                        do iphase = 1, Nphase
        !WE SHOULD USE NPRES IN ANY CASE, NOT PHASES HERE...
                        call addto(Mmat%C_PETSC, idim, iphase, IU_NOD, P_JNOD, Mmat%C( idim, iphase, COUNT ))
                        end do 
                    end do
                end do
            end do
            call assemble( Mmat%C_PETSC )
        end if
        !Now CT matrix
        if (associated(Mmat%CT_PETSC%refcount)) call deallocate(Mmat%CT_PETSC)
        sparsity=>extract_csr_sparsity(packed_state,"CTMatrixSparsity")
        call allocate(Mmat%CT_PETSC,sparsity,[NPHASE,Mdims%ndim],name="CT_PETSC"); call zero(Mmat%CT_PETSC)
        
        DO P_JNOD = 1, Mdims%CV_NONODS
            DO COUNT = mspars%CT%fin( P_JNOD ), mspars%CT%fin( P_JNOD + 1 ) - 1
                IU_NOD = mspars%CT%col( COUNT )
                do idim =1, Mdims%ndim
                    !WE SHOULD USE NPRES IN ANY CASE, NOT PHASES HERE...
                    do iphase = 1, Nphase
                        call addto(Mmat%CT_PETSC, iphase, idim, P_JNOD, IU_NOD, Mmat%CT( idim, iphase, COUNT ))
                    end do
                end do
            END DO
        END DO
        call assemble( Mmat%CT_PETSC )

        !Now Mass matrix
        ! if (associated(Mmat%PIVIT_PETSC%refcount)) call deallocate(Mmat%PIVIT_PETSC)
        ! sparsity=>extract_csr_sparsity(packed_state,"MomentumSparsity")
        ! Mmat%PIVIT_PETSC = allocate_momentum_matrix(sparsity, velocity, nphase)
        ! !Try to make it diagonal...
        ! call allocate(Mmat%PIVIT_PETSC,sparsity,[Mdims%u_nloc*Mdims%ndim*NPHASE,1],"PIVIT_PETSC"); call zero(Mmat%PIVIT_PETSC)
        
        ! DO ELE = 1, Mdims%totele
        !     do j = 1, Mdims%u_nloc * nphase *Mdims%ndim
        !         call addto(Mmat%PIVIT_PETSC, j, 1, ELE, ELE,  Mmat%PIVIT_MAT( 1, 1, ele ))
        !     end do
        ! END DO
        ! call assemble( Mmat%PIVIT_PETSC )



        ! call MatView(Mmat%C_PETSC%M, PETSC_VIEWER_STDOUT_SELF, idim)
        ! call MatView(Mmat%CT_PETSC%M, PETSC_VIEWER_STDOUT_SELF, idim)
        ! call MatView(Mmat%PIVIT_PETSC%M, PETSC_VIEWER_STDOUT_SELF, idim)
    END SUBROUTINE Convert_C_and_CT_mat_to_PETSc_format

end module solvers_module
