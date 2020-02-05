
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

    !use multiphase_1D_engine, only: getOverrelaxation_parameter 
    !use multiphase_EOS, only:Get_DevCapPressure 

    implicit none

#include "petsc_legacy.h"

    private

    public :: BoundedSolutionCorrections, FPI_backtracking, Set_Saturation_to_sum_one,&
         Ensure_Saturation_sums_one, auto_backtracking


contains


    subroutine BoundedSolutionCorrections( state, packed_state, &
        Mdims, CV_funs, small_findrm, small_colm, Field_name, &
        for_sat, min_max_limits)
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
    subroutine FPI_backtracking(nphase, Mdims, ndgln, state, packed_state, sat_bak, backtrack_sat, backtrack_par_from_schema, &
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
        call Set_Saturation_to_sum_one(mdims, ndgln, packed_state, state, do_not_update_halos = .TRUE. )
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
        call Set_Saturation_to_sum_one(mdims, ndgln, packed_state, state, do_not_update_halos = .FALSE. )

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

    subroutine Set_Saturation_to_sum_one(mdims, ndgln, packed_state, state, do_not_update_halos)
        !This subroutines eliminates the oscillations in the saturation that are bigger than a
        !certain tolerance and also sets the saturation to be between bounds
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

    subroutine Ensure_Saturation_sums_one(mdims, ndgln, packed_state, find_scapegoat_phase)
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

    end subroutine Ensure_Saturation_sums_one


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
   

    subroutine AI_backtracking_parameters(Mdims, ndgln, packed_state, state, courant_number_in, for_transport)
        ! ...
        implicit none
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_ndgln), intent(in) :: ndgln
        type( state_type ), intent(inout) :: packed_state
        type( state_type ), dimension( : ), intent( inout ) :: state
        real, dimension(:), intent(inout) :: courant_number_in
        logical, optional, intent(in) :: for_transport
        !Local variables
        real, save :: backup_shockfront_Courant = 0.
        logical, save :: readed_options = .false.
        logical, save :: readed_perm = .false.
        !logical, save :: readed_por = .false.
        logical, save :: readed_sizes = .false.
        real, dimension(:,:), pointer :: X_ALL, saturation, Imble_frac, cap_entry_pres, end_point_relperm
        type(tensor_field), pointer :: permeability, viscosity, density, velocity
        type(vector_field), pointer :: porosity 
        type (vector_field_pointer), dimension(Mdims%nphase) :: Darcy_velocity
        real, dimension(Mdims%cv_nonods) :: total_mobility, counter_subcv, nDarcy_velocity_cvwise
        real, dimension(Mdims%ndim, Mdims%cv_nonods) :: Darcy_velocity_cvwise
        real, dimension(Mdims%n_in_pres, Mdims%cv_nonods) :: relperm
        real, dimension(:), allocatable :: longitudinal_capilary, transverse_capilary, buoyancy_number, transverse_buoyancy, overrelaxation, Peclet
        logical, save :: gravity, cap_pressure, black_oil, ov_relaxation, one_phase, wells
        real, save :: n_phases, n_components, courant_number, shockfront_courant_number, aspect_ratio
        real, save :: min_shockfront_mobratio, max_shockfront_mobratio, average_shockfront_mobratio
        real, save :: average_longitudinal_capilary, average_transverse_capilary, max_longitudinal_capilary, max_transverse_capilary, min_longitudinal_capilary, min_transverse_capilary
        real, save :: average_buoyancy_number, average_transverse_buoyancy, max_buoyancy_number, max_transverse_buoyancy, min_buoyancy_number, min_transverse_buoyancy 
        real, save :: average_overrelaxation, max_overrelaxation, min_overrelaxation
        real, save :: average_Peclet, max_Peclet, min_Peclet
        real ::  average_perm_length, average_perm_thickness, average_por, l1_max, l1_min, l2_max, l2_min, h_max, h_min, aux_shockfront_mobratio, domain_length, domain_thickness, delta_density, diffusivity
        integer :: iphase, cv_nodi, cv_iloc, ele, total_ele, u_iloc, u_inod, total_cv, counter, Phase_with_Pc
        
        !*************************************!
        !!! ***Getting support variables*** !!!
        !*************************************!

        !Permeability averages 
        if (.not.readed_perm) then 
            !Retrieve permeability field - element wise
            permeability => extract_tensor_field(packed_state,"Permeability")  
            !Average permeability along Mdims%ndim-1 first dimensions
            average_perm_length = 0.
            if (Mdims%ndim < 3) then
                do ele = 1, Mdims%totele
                    average_perm_length = average_perm_length + permeability%val(1,1,ele) 
                end do 
                total_ele = Mdims%totele
                if (IsParallel()) then
                    call allsum(average_perm_length)
                    call allsum(total_ele)
                end if
                average_perm_length = average_perm_length/total_ele 
            else 
                do ele = 1, Mdims%totele
                    average_perm_length = average_perm_length + permeability%val(1,1,ele)  + permeability%val(2,2,ele) 
                end do 
                total_ele = Mdims%totele
                if (IsParallel()) then
                    call allsum(average_perm_length)
                    call allsum(total_ele)
                end if
                average_perm_length = average_perm_length/(total_ele*2) 
            end if  
            !Average permeability along the last dimensions
            average_perm_thickness = 0.
            do ele = 1, Mdims%totele
                average_perm_thickness = average_perm_thickness + permeability%val(Mdims%ndim,Mdims%ndim,ele) 
            end do 
            total_ele = Mdims%totele
            if (IsParallel()) then
                call allsum(average_perm_thickness)
                call allsum(total_ele)
            end if
            average_perm_thickness = average_perm_thickness/total_ele     
            readed_perm = .true.
        end if    

        !Porosity average - not using
        !if (.not.readed_por) then 
        !    !Retrieve porosity field
        !    porosity => extract_vector_field(packed_state,"Porosity")
        !    average_por = 0
        !    do ele = 1, Mdims%totele    
        !        average_por = average_por + porosity%val(1, ele) 
        !        total_ele = Mdims%totele
        !    if (IsParallel()) then
        !       call allsum(average_por)
        !       call allsum(total_ele)
        !    average_por = average_por/real(total_ele)    
        !    readed_por = .true.    
        !end if
            
        !Domain lemgth, domain thickness 
        if (.not.readed_sizes) then
            !Extract variables from packed_state, X_ALL => extract_vector_field( packed_state, "PressureCoordinate" )
            call get_var_from_packed_state(packed_state, PressureCoordinate = X_ALL)    
            !Domain length as the maximum distance in the (Mdims%ndim-1) dimensions, we only calculate it once
            if (Mdims%ndim < 3) then 
                l1_max = maxval(X_ALL(1,:)) 
                l1_min = minval(X_ALL(1,:))   
                if (IsParallel()) then
                    call allmax(l1_max)
                    call allmin(l1_min)
                end if
                domain_length = abs(l1_max-l1_min)
            else 
                l1_max = maxval(X_ALL(1,:)) 
                l1_min = minval(X_ALL(1,:))
                l2_max = maxval(X_ALL(2,:)) 
                l2_min = minval(X_ALL(2,:))
                if (IsParallel()) then
                    call allmax(l1_max)
                    call allmin(l1_min)
                    call allmax(l2_max)
                    call allmin(l2_min)
                end if
                domain_length = max(l1_max-l1_min, l2_max-l2_min)
            end if
            !Domain thickness as the maximum distance in the last dimension
            h_max = maxval(X_ALL(Mdims%ndim,:)) 
            h_min = minval(X_ALL(Mdims%ndim,:))
            if (IsParallel()) then
                call allmax(h_max)
                call allmin(h_min)
            end if
            domain_thickness = abs(h_max-h_min)  
            readed_sizes = .true.    
        end if

        !Retrieve relative permeabilty field - control-volume wise
        relperm = 1.0
        !Retrieve viscosity field - control-volume wise
        viscosity => extract_tensor_field( packed_state, "Viscosity")
        !Retrieve phase density - control-volume wise
        density  => extract_tensor_field( packed_state, "PackedDensity" )
        !Retrieve velocity - control-volume wise
        velocity => extract_tensor_field( packed_state, "PackedVelocity" )
        !Retrieve Immobile fraction - element wise 
        call get_var_from_packed_state(packed_state, Immobile_fraction = Imble_frac)
        !Retrive saturation - control-volume wise
        call get_var_from_packed_state(packed_state, PhaseVolumeFraction = saturation) 
        !Retrive Capillary entry pressure - control-volume wise
        call get_var_from_packed_state(packed_state, Cap_entry_pressure = cap_entry_pres) 
        !Retrive end-point relative permeability- control-volume wise
        call get_var_from_packed_state(packed_state, EndPointRelperm = end_point_relperm) 
        !Calculate the total number of control volumes  
        total_cv = Mdims%cv_nonods 
        if (IsParallel()) then
            call allsum(total_cv)                
        end if

        !Total mobility - control-volume wise
        do cv_nodi = 1, Mdims%cv_nonods
            total_mobility(cv_nodi) = sum( relperm(:, cv_nodi) / viscosity%val( 1, :, cv_nodi) )
        end do

        !P0 Darcy velocity - control-volume wise
        do iphase = 1, Mdims%n_in_pres
            Darcy_velocity(iphase)%ptr => extract_vector_field(state(iphase),"DarcyVelocity")
        end do    
        Darcy_velocity_cvwise = 0.
        counter_subcv = 0
        do ele = 1, Mdims%totele
            do u_iloc = 1, Mdims%u_nloc
                u_inod = ndgln%u(( ele - 1 ) * Mdims%u_nloc + u_iloc)
                do cv_iloc = 1, Mdims%cv_nloc
                    cv_nodi = ndgln%cv(( ele - 1) * Mdims%cv_nloc + cv_iloc )
                    do iphase = 1, Mdims%n_in_pres
                        Darcy_velocity_cvwise(:,cv_nodi) = Darcy_velocity_cvwise(:,cv_nodi) + Darcy_velocity(iphase)%ptr%val(:,u_inod)/Mdims%u_nloc  
                    end do
                    counter_subcv(cv_nodi) = counter_subcv(cv_nodi) + 1
                end do
            end do
        end do
        do cv_nodi = 1, Mdims%cv_nonods
            Darcy_velocity_cvwise(:,cv_nodi) = Darcy_velocity_cvwise(:,cv_nodi)/counter_subcv(cv_nodi)
            nDarcy_velocity_cvwise(cv_nodi) = norm2(Darcy_velocity_cvwise(:,cv_nodi))
        end do
        

        !**************************************!
        !!! ***Calculating all parameters*** !!!
        !**************************************!

        if (.not.readed_options) then
            !We read the options just once
            !!! Number of phases !!!
            n_phases = Mdims%n_in_pres 
            !!! have wells !!! 
            wells = (Mdims%npres > 1)
            !!! Number of components !!!
            n_components = Mdims%ncomp        
            !!! gravity !!!
            gravity = have_option("/physical_parameters/gravity")        
            !!! Capillary pressure !!!
            cap_pressure = have_option_for_any_phase("/multiphase_properties/capillary_pressure", Mdims%nphase)
            !!! Single-phase flow !!!
            one_phase = (Mdims%n_in_pres == 1) 
            !!! Black Oil
            black_oil = have_option( "/physical_parameters/black-oil_PVT_table")
            !!! Over relaxation !!!
            ov_relaxation = have_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation')  
            readed_options = .true.
        end if

        !!! Courant number and Shockfront courant number !!!
        !Sometimes the shock-front courant number is not well calculated, then use previous value
        if (abs(courant_number_in(2)) < 1d-8 ) courant_number_in(2) = backup_shockfront_Courant
        !Courant number
        courant_number = courant_number_in(1) 
        !shockfront courant number
        shockfront_courant_number = courant_number_in(2)
        backup_shockfront_Courant = shockfront_courant_number

        !!! Effective aspect ratio !!!
        aspect_ratio = domain_length*SQRT(average_perm_thickness/average_perm_length)/domain_thickness 
        
        !!! shockfront mobility ratio !!! 
        min_shockfront_mobratio = 999.
        max_shockfront_mobratio = 0.
        average_shockfront_mobratio = 0.
        counter = 0
        if (Mdims%n_in_pres > 1) then
            do ele = 1, Mdims%totele
                if (shock_front_in_ele(ele, Mdims, saturation, ndgln, Imble_frac)) then 
                    aux_shockfront_mobratio = shock_front_mobility_ratio(ele, Mdims, saturation, ndgln, Imble_frac, total_mobility)
                    min_shockfront_mobratio = min(min_shockfront_mobratio, aux_shockfront_mobratio)
                    max_shockfront_mobratio = max(max_shockfront_mobratio, aux_shockfront_mobratio)
                    average_shockfront_mobratio = average_shockfront_mobratio + aux_shockfront_mobratio
                    counter = counter + 1
                end if
            end do  
            if (IsParallel()) then
                call allmax(max_shockfront_mobratio)
                call allmin(min_shockfront_mobratio)
                call allsum(average_shockfront_mobratio)
                call allsum(counter)                
            end if
            average_shockfront_mobratio = average_shockfront_mobratio/counter
        else 
            min_shockfront_mobratio = 1.
            max_shockfront_mobratio = 1.
            average_shockfront_mobratio = 1.      
        end if    
        
        !!! Capilary numbers !!!
        if (cap_pressure) then
            allocate( longitudinal_capilary(Mdims%cv_nonods) )
            allocate( transverse_capilary(Mdims%cv_nonods) )
            do cv_nodi = 1, Mdims%cv_nonods
                longitudinal_capilary(cv_nodi) = average_perm_length*end_point_relperm(Mdims%n_in_pres,cv_nodi)*sum(cap_entry_pres(:,cv_nodi)) / &
                                                & (nDarcy_velocity_cvwise(cv_nodi)*viscosity%val( 1, Mdims%n_in_pres, cv_nodi)*domain_length) 
                transverse_capilary(cv_nodi) = average_perm_thickness*end_point_relperm(Mdims%n_in_pres,cv_nodi)*sum(cap_entry_pres(:,cv_nodi))*domain_length / &
                                                & (nDarcy_velocity_cvwise(cv_nodi)*viscosity%val( 1, Mdims%n_in_pres, cv_nodi)*domain_thickness**2)               
            end do 
            ! calculate average, max and min values
            average_longitudinal_capilary = sum(longitudinal_capilary)
            average_transverse_capilary = sum(transverse_capilary)
            max_longitudinal_capilary = maxval(longitudinal_capilary)
            max_transverse_capilary = maxval(transverse_capilary)
            min_longitudinal_capilary = minval(longitudinal_capilary)
            min_transverse_capilary = minval(transverse_capilary)
            if (IsParallel()) then
                call allsum(average_longitudinal_capilary)
                call allsum(average_transverse_capilary)
                call allmax(max_longitudinal_capilary)
                call allmax(max_transverse_capilary)  
                call allmin(min_longitudinal_capilary)
                call allmin(min_transverse_capilary)               
            end if
            average_longitudinal_capilary = average_longitudinal_capilary/total_cv
            average_transverse_capilary = average_transverse_capilary/total_cv 
            deallocate(longitudinal_capilary)
            deallocate(transverse_capilary)
        else 
            average_longitudinal_capilary = 0.
            average_transverse_capilary = 0.
            max_longitudinal_capilary = 0.
            max_transverse_capilary = 0.
            min_longitudinal_capilary = 0.
            min_transverse_capilary = 0.
        end if

        !!! gravity numbers !!!
        if (gravity) then
            allocate( buoyancy_number(Mdims%cv_nonods) )
            allocate( transverse_buoyancy(Mdims%cv_nonods) )
            do cv_nodi = 1, Mdims%cv_nonods
                delta_density = (maxval(density%val(1,:,cv_nodi)) - minval(density%val(1,:,cv_nodi)))*9.81 ! ************** cos and sin alpha
                buoyancy_number(cv_nodi) = average_perm_length*end_point_relperm(Mdims%n_in_pres,cv_nodi)*delta_density*domain_thickness / &
                                            & (nDarcy_velocity_cvwise(cv_nodi)*viscosity%val( 1, Mdims%n_in_pres, cv_nodi)*domain_length) 
                !longitudinal_buoyancy(cv_nodi) = average_perm_length*end_point_relperm(Mdims%n_in_pres,cv_nodi)*delta_density / &
                !                               & (nDarcy_velocity_cvwise(cv_nodi)*viscosity%val( 1, Mdims%n_in_pres, cv_nodi)) 
                transverse_buoyancy(cv_nodi) = average_perm_thickness*end_point_relperm(Mdims%n_in_pres,cv_nodi)*delta_density*domain_length / &
                                                & (nDarcy_velocity_cvwise(cv_nodi)*viscosity%val( 1, Mdims%n_in_pres, cv_nodi)*domain_thickness)               
            end do 
            ! calculate average, max and min values
            average_buoyancy_number = sum(buoyancy_number)
            average_transverse_buoyancy = sum(transverse_buoyancy)
            max_buoyancy_number = maxval(buoyancy_number)
            max_transverse_buoyancy = maxval(transverse_buoyancy)
            min_buoyancy_number = minval(buoyancy_number)
            min_transverse_buoyancy = minval(transverse_buoyancy)
            if (IsParallel()) then
                call allsum(average_buoyancy_number)
                call allsum(average_transverse_buoyancy)
                call allmax(max_buoyancy_number)
                call allmax(max_transverse_buoyancy)  
                call allmin(min_buoyancy_number)
                call allmin(min_transverse_buoyancy)               
            end if
            average_buoyancy_number = average_buoyancy_number/total_cv
            average_transverse_buoyancy = average_transverse_buoyancy/total_cv 
            deallocate(buoyancy_number)
            deallocate(transverse_buoyancy)
        else 
            average_buoyancy_number = 0.
            average_transverse_buoyancy = 0.
            max_buoyancy_number = 0.
            max_transverse_buoyancy = 0.
            min_buoyancy_number = 0.
            min_transverse_buoyancy = 0.
        end if

        !!! overrelaxation parameter !!!
        if (ov_relaxation) then
            allocate( overrelaxation(Mdims%cv_nonods) ) 
            call getOverrelaxation_parameter(state, packed_state, Mdims, ndgln, overrelaxation, Phase_with_Pc) 
            max_overrelaxation = maxval(overrelaxation)
            min_overrelaxation = minval(overrelaxation)
            average_overrelaxation = sum(overrelaxation)
            if (IsParallel()) then
                call allmax(max_overrelaxation)
                call allmin(min_overrelaxation)
                call allsum(average_overrelaxation)              
            end if
            average_overrelaxation = average_overrelaxation/total_cv
            deallocate(overrelaxation)
        else 
            max_overrelaxation = 0.
            min_overrelaxation = 0.
            average_overrelaxation = 0.
        end if

        !!! Peclet number !!!
        !Peclet = V * L / Diffusivity;
        if (ov_relaxation) then 
            allocate( Peclet(Mdims%cv_nonods) )
            if (present_and_true(for_transport)) then
                call get_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation/Vanishing_for_transport', diffusivity)
            else
                call get_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation', diffusivity)
            end if
            diffusivity = abs(diffusivity)
            Peclet = 0.
            do cv_nodi = 1, Mdims%cv_nonods
                Peclet(cv_nodi) = nDarcy_velocity_cvwise(cv_nodi)*domain_length/diffusivity
            end do
            max_Peclet = maxval(Peclet)
            min_Peclet = minval(Peclet)
            average_Peclet = sum(Peclet)
            if (IsParallel()) then
                call allmax(max_Peclet)
                call allmin(min_Peclet)
                call allsum(average_Peclet)
            end if
            average_Peclet = average_Peclet/total_cv
            deallocate(Peclet)
        else
            max_Peclet = 0.    ! ****************** the default value is 0 or inf (if inf I may use 1_over_Peclet = 1/Peclet)
            min_Peclet = 0.
            average_Peclet = 0.
        end if

        !*****************************!
        !!! ***Support functions*** !!!
        !*****************************!

        contains 

            real function shock_front_mobility_ratio(ele, Mdims, sat, ndgln, Imble_frac, total_mobility)
                !Detects whether the element has a shockfront or not and calculates the shockfront mobility ratio
                !if there is not a shock front within the lement the fuction returns -1
                implicit none
                integer, intent(in) :: ele 
                type(multi_dimensions), intent(in) :: Mdims
                real, dimension(:,:), intent(in) :: sat 
                type(multi_ndgln), intent(in) :: ndgln
                real, dimension(:,:), intent(in) :: Imble_frac  
                real, dimension(:), intent(in) :: total_mobility           
                !Local variables
                integer :: iphase, cv_iloc, cv_nodi
                real :: minival, maxival, aux_sat, tmob_aheadfront, tmob_atfront
                real, parameter :: tol = 0.05 !Shock fronts smaller than this are unlikely to require extra handling
                
                !Starts the value                
                shock_front_mobility_ratio = -1.0

                minival = 999.; maxival = 0.
                do iphase = 1, mdims%n_in_pres - 1 
                    do cv_iloc = 1, Mdims%cv_nloc
                        cv_nodi = ndgln%cv((ele-1)*Mdims%cv_nloc+cv_iloc)
                        aux_sat = sat(iphase, cv_nodi) - Imble_frac(iphase, ele)    
                        if (aux_sat < minival) then 
                            minival = aux_sat
                            tmob_aheadfront = total_mobility(cv_nodi)    
                        end if
                        if (aux_sat > maxival) then 
                            maxival = aux_sat
                            tmob_atfront = total_mobility(cv_nodi)    
                        end if
                    end do
                    if  (minival < tol .and. (maxival-minival) > tol) then
                        shock_front_mobility_ratio = tmob_atfront/tmob_aheadfront
                        exit
                    end if
                end do
            end function shock_front_mobility_ratio


            logical function shock_front_in_ele(ele, Mdims, sat, ndgln, Imble_frac)
                !Detects whether the element has a shockfront or not
                implicit none
                integer, intent(in) :: ele !***
                type(multi_dimensions), intent(in) :: Mdims
                real, dimension(:,:), intent(in) :: sat 
                type(multi_ndgln), intent(in) :: ndgln
                real, dimension(:,:), intent(in) :: Imble_frac             
                !Local variables
                integer :: iphase, cv_iloc
                real :: minival, maxival, aux
                real, parameter :: tol = 0.05 !Shock fronts smaller than this are unlikely to require extra handling

                !Starts the value   
                shock_front_in_ele = .false.

                minival = 999.; maxival = 0.
                do iphase = 1, mdims%n_in_pres - 1 ! ***** Changed from nphase to n_in_pres / Why n_in_pres -1 ???? Only injected fluid, the displaced is the last one? 
                    do cv_iloc = 1, Mdims%cv_nloc
                        aux = sat(iphase, ndgln%cv((ele-1)*Mdims%cv_nloc+cv_iloc)) - Imble_frac(iphase,ele) !***
                        minival = min(aux, minival)
                        maxival = max(aux, maxival)
                    end do
                    if  (minival < tol .and. (maxival-minival) > tol) then
                        shock_front_in_ele = .true.
                        exit
                    end if
                end do
            end function shock_front_in_ele

            subroutine getOverrelaxation_parameter(state, packed_state, Mdims, ndgln, Overrelaxation, Phase_with_Pc, for_transport)
                !This subroutine calculates the overrelaxation parameter we introduce in the saturation equation
                !It is the derivative of the capillary pressure for each node.
                !Overrelaxation has to be alocate before calling this subroutine its size is cv_nonods
                implicit none
                type( state_type ), dimension( : ), intent( inout ) :: state
                type( state_type ), intent(inout) :: packed_state
                type (multi_dimensions), intent(in) :: Mdims
                type(multi_ndgln), intent(in) :: ndgln
                real, dimension(:), intent(inout) :: Overrelaxation
                integer, intent(inout) :: Phase_with_Pc
                logical, optional, intent(in) :: for_transport
                !Local variables
                real, save :: domain_length = -1
                integer, save :: Cap_pressure_relevant = -1
                integer :: iphase, nphase, cv_nodi, cv_nonods, u_inod, cv_iloc, ele, u_iloc
                real :: Pe_aux, parl_max, parl_min, Pe_max, Pe_min
                real, dimension(:), pointer ::Pe, Cap_exp
                logical :: Artificial_Pe
                real, dimension(:,:,:), pointer :: p
                real, dimension(:,:), pointer :: satura, immobile_fraction, Cap_entry_pressure, Cap_exponent, X_ALL
                type( tensor_field ), pointer :: Velocity
                type( scalar_field ), pointer :: PIPE_Diameter
           
           
                !Extract variables from packed_state
                call get_var_from_packed_state(packed_state,FEPressure = P,&
                    PhaseVolumeFraction = satura, immobile_fraction = immobile_fraction, PressureCoordinate = X_ALL)
           
                !Initiate local variables
                nphase = size(satura,1)
                cv_nonods = size(satura,2)
           
                !#######Only apply this method if it has been explicitly invoked through Pe_stab or
                !non-consistent capillary pressure!######
                Phase_with_Pc = -1
                if (.not. have_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation') ) then
                    Overrelaxation = 0.0; Phase_with_Pc = -10
                    return
                end if
                !If this is for transport, check if we want to apply it
                if (present_and_true(for_transport)) then
                   if (.not.have_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation/Vanishing_for_transport')) then
                       Overrelaxation = 0.0; Phase_with_Pc = -10
                       return
                   end if
                end if
                !######################################################################################
               !By default phase 1
               Phase_with_Pc = 1
                !Check capillary pressure options
                do iphase = Nphase, 1, -1!Going backwards since the wetting phase should be phase 1
                    !this way we try to avoid problems if someone introduces 0 capillary pressure in the second phase
                    if (have_option( "/material_phase["//int2str(iphase-1)//"]/multiphase_properties/capillary_pressure" )) then
                        Phase_with_Pc = iphase
                    end if
                end do
                Artificial_Pe = .false.
                Cap_exponent => null(); Cap_entry_pressure => null()!Initialize
                if (Phase_with_Pc>0) then
                    !Get information for capillary pressure to be used
                    if ( (have_option("/material_phase["//int2str(Phase_with_Pc-1)//&
                        "]/multiphase_properties/capillary_pressure/type_Brooks_Corey") ) .or. (have_option("/material_phase["//int2str(Phase_with_Pc-1)//&
                        "]/multiphase_properties/capillary_pressure/type_Power_Law") ) )then
                        call get_var_from_packed_state(packed_state, Cap_entry_pressure = Cap_entry_pressure,&
                            Cap_exponent = Cap_exponent)!no need for the imbibition because we need the derivative which will be zero as it is a constant
                    end if
                    !If we want to introduce a stabilization term, this one is imposed over the capillary pressure.
                    !Unless we are using the non-consistent form of the capillary pressure
                    if ( have_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation') ) then
                        allocate(Pe(CV_NONODS), Cap_exp(CV_NONODS))
                        Artificial_Pe = .true.
                        if (present_and_true(for_transport)) then
                           call get_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation/Vanishing_for_transport', Pe_aux)
                        else
                           call get_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/Vanishing_relaxation', Pe_aux)
                        end if
           
           
                       !Check if the capillary pressure introduced is important enough to actually trigger the VAD for Capillary pressure
                        if ( associated(Cap_entry_pressure) .and. Cap_pressure_relevant < 0) then
                           Cap_pressure_relevant = 0
                           if (maxval(Cap_entry_pressure)/maxval(P) > 1e-2) Cap_pressure_relevant = 1
                       end if
           
                        if (Cap_pressure_relevant > 0) then
                            Cap_exp = 2.0 !Quadratic exponent
           
                        else
                            Cap_exp = 1.!Linear exponent
                        end if
           
                        if (Pe_aux<0) then!Automatic set up for Pe
                            !Method based on calculating an entry pressure for a given Peclet number;
                            !Peclet = V * L / Diffusivity; We consider only the entry pressure for the diffusivity
                            !Pe = Vel * L/ Peclet. At present we are using the velocity that includes the sigma. Maybe might be worth it using the Darcy velocity?
                            Velocity => extract_tensor_field( packed_state, "PackedVelocity" )
           
                            !Since it is an approximation, the domain length is the maximum distance, we only calculate it once
                            if (domain_length < 0) then
                                parl_max = maxval(X_ALL)
                                parl_min = minval(X_ALL)
                                if (IsParallel()) then
                                    call allmax(parl_max)
                                    call allmin(parl_min)
                                end if
                                domain_length = abs(parl_max-parl_min)
                            end if
                            Pe_aux = abs(Pe_aux)
                             !Obtain an approximation of the capillary number to obtain an entry pressure
                            Pe = 0.
                            do ele = 1, Mdims%totele
                                do u_iloc = 1, Mdims%u_nloc
                                    u_inod = ndgln%u(( ELE - 1 ) * Mdims%u_nloc +u_iloc )
                                    do cv_iloc = 1, Mdims%cv_nloc
                                        cv_nodi = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                                        Pe(cv_nodi) = Pe(cv_nodi) + (1./Pe_aux) * (sum(abs(Velocity%val(:,Phase_with_Pc,u_inod)))/real(Mdims%ndim) * domain_length)/real(Mdims%u_nloc)
                                    end do
                                end do
                            end do
                            Pe_max = maxval(Pe)
                            Pe_min = minval(Pe)
                            if (IsParallel()) then
                                call allmax(Pe_max)
                                call allmin(Pe_min)
                            end if
           
                            if (have_option('/numerical_methods/VAD_two_levels')) then
                              !Homogenise using two values only, if below the average value then use the lowest value possible (this is for stratified flows)
                              where (Pe > 0.5*Pe_max + 0.5*Pe_min)
                                Pe = Pe_max
                              elsewhere
                                Pe = Pe_min
                              end where
                            else
                              !Homogenise the value, this seems to be better to avoid problems (method used for the CMAME paper)
                              Pe = (Pe_max+Pe_min)/2.
                            end if
                        else
                            Pe = Pe_aux
                        end if
           
                    end if
           
                    !Calculate the overrrelaxation parameter, the numbering might be different for Pe and real capillary
                    !values, hence we calculate it differently
                    if (Artificial_Pe) then
                        !Calculate the Overrelaxation
                        do ele = 1, Mdims%totele
                            do cv_iloc = 1, Mdims%cv_nloc
                                cv_nodi = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                                Overrelaxation(CV_NODI) =  Get_DevCapPressure(satura(Phase_with_Pc, CV_NODI),&
                                    Pe(CV_NODI), Cap_Exp(CV_NODI), immobile_fraction(:,ele), Phase_with_Pc, nphase)
                            end do
                        end do
                    else
                        !Calculate the Overrelaxation
                        do ele = 1, Mdims%totele
                            do cv_iloc = 1, Mdims%cv_nloc
                                cv_nodi = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                                Overrelaxation(CV_NODI) =  Get_DevCapPressure(satura(Phase_with_Pc, CV_NODI),&
                                    Cap_entry_pressure(Phase_with_Pc, ele), &
                                    Cap_exponent(Phase_with_Pc, ele),&
                                    immobile_fraction(:,ele), Phase_with_Pc, nphase)
                            end do
                        end do
                    end if
                else
                    Overrelaxation = 0.0
                end if
           
                !DISABLED THIS SECTION BELOW AS WITH NONDEBUGGING THE RESULTS ARE AFFECTED
                ! if (Mdims%npres >1) THEN !If we have pipes, reduce VAD in the CVs that have pipes, they do not get along...
                !   PIPE_Diameter => EXTRACT_SCALAR_FIELD(state(1), "DiameterPipe")
                !   do cv_nodi = 1, Mdims%cv_nonods
                !     IF ( PIPE_DIAMETER%VAL(CV_NODI) > 1e-8 ) THEN
                !        Overrelaxation(CV_NODI) = Overrelaxation(CV_NODI) * 1e-2!Severely reduce Overrelaxation around wells
                !     end if
                !   end do
                ! end if
           
           
                !Deallocate
                if (Artificial_Pe) then
                    deallocate(Pe, Cap_exp)
                end if
                nullify(Pe, Cap_exp)
           
            end subroutine getOverrelaxation_parameter

            real function Get_DevCapPressure(sat, Pe, a, immobile_fraction, iphase, nphase)
                !This functions returns the derivative of the capillary pressure with respect to the saturation
                Implicit none
                integer, intent(in) :: iphase, nphase
                real, intent(in) :: sat, Pe, a
                real, dimension(:), intent(in) :: immobile_fraction
                !Local
                real, parameter :: eps = 1d-3
                real :: aux
                integer :: i
                logical, save :: Cap_Brooks = .true., Cap_Power = .false.
                logical, save :: first_time = .true.
        
                aux = ( 1.0 - sum(immobile_fraction(:)) )
                ! Determine which capillary pressure model is to be used for overrelaxation. Use Brooks-Corey unless power_law Pc activated (important to allow overelax even when Pc is off).
                if (first_time) then
                        Cap_Power = have_option_for_any_phase("/multiphase_properties/capillary_pressure/type_Power_Law", nphase)
                        Cap_Brooks = .not. (Cap_Power)
                        first_time = .false.
                end if
        
                if(Cap_Power) then
                    Get_DevCapPressure = &
                        -a*Pe/(1.0 - sum(Immobile_fraction(:)) )  * ( 1.0 - ( sat - Immobile_fraction(iphase) )/( 1.0 - sum(Immobile_fraction(:)) ) ) **(a-1)
                else
                    Get_DevCapPressure = &
                        -a * Pe * aux**a * min((sat - immobile_fraction(iphase) + eps), 1.0) ** (-a-1)
                endif
    
            end function Get_DevCapPressure

    end subroutine AI_backtracking_parameters

end module solvers_module
