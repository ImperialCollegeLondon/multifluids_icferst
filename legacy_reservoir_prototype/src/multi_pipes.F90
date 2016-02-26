
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

module multi_pipes

    use fldebug

    use fields

    use reference_counting
    use memory_diagnostics

    use global_parameters, only: pi, after_adapt

    use Petsc_Tools
    use Sparse_tools
    use sparse_tools_petsc


    use shape_functions_Linear_Quadratic
    use shape_functions_NDim
    use shape_functions_prototype
    use matrix_operations
    use Copy_Outof_State
    use boundary_conditions

    use multi_tools
    use multi_data_types

    implicit none

    private

    public  :: MOD_1D_CT_AND_ADV, MOD_1D_FORCE_BAL_C

    !Parameters for boundary conditions copied from cv_adv-diff.F90
    INTEGER, PARAMETER :: WIC_T_BC_DIRICHLET = 1, WIC_T_BC_ROBIN = 2, &
        WIC_T_BC_DIRI_ADV_AND_ROBIN = 3, WIC_D_BC_DIRICHLET = 1, &
        WIC_U_BC_DIRICHLET = 1, &
        WIC_U_BC_ROBIN = 2, &
        WIC_U_BC_DIRI_ADV_AND_ROBIN = 3, &
        WIC_U_BC_DIRICHLET_INOUT = 2, &
        WIC_P_BC_DIRICHLET = 1, &
        WIC_P_BC_FREE = 2


contains


    SUBROUTINE MOD_1D_CT_AND_ADV( state, packed_state, nphase, npres, n_in_pres, ndim, u_nloc, cv_nloc, x_nloc, FINACV, COLACV, &
        U_NONODS,U_SNLOC,CV_SNLOC,STOTEL,CV_SNDGLN,U_SNDGLN, WIC_T_BC_ALL,WIC_D_BC_ALL,WIC_U_BC_ALL, SUF_T_BC_ALL,SUF_D_BC_ALL,SUF_U_BC_ALL, &
        cv_nonods, getcv_disc, getct, petsc_acv, totele, cv_ndgln, x_ndgln, u_ndgln, mat_ndgln, ct, c, findct, colct, findc, colc, CV_RHS_field, CT_RHS, &
        findcmc, colcmc, MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE, mass_pipe, MASS_PIPE_FOR_COUP, &
        INV_SIGMA, INV_SIGMA_NANO, OPT_VEL_UPWIND_COEFS_NEW )

        ! This sub modified either CT or the Advection-diffusion equation for 1D pipe modelling

        type(state_type), intent(in) :: packed_state
        type(state_type), dimension(:), intent(in) :: state

        INTEGER, intent( in ) :: nphase, npres, n_in_pres, ndim, u_nloc, cv_nloc, x_nloc, cv_nonods, totele, CV_SNLOC, U_SNLOC, STOTEL, U_NONODS
        integer, dimension(:), intent( in ), target :: FINACV, COLACV
        integer, dimension(:), intent( in ) :: cv_ndgln, x_ndgln, u_ndgln, mat_ndgln, findct, colct, findc, colc, findcmc, colcmc
        integer, dimension(:), intent( in ) :: CV_SNDGLN, U_SNDGLN
        integer, dimension(:,:,:), intent( in ) :: WIC_T_BC_ALL, WIC_D_BC_ALL, WIC_U_BC_ALL
        real, dimension(:,:,:), intent( in ) :: SUF_T_BC_ALL,SUF_D_BC_ALL,SUF_U_BC_ALL
        real, dimension(:,:,:,:), intent( in ) :: OPT_VEL_UPWIND_COEFS_NEW
        real, dimension(:,:,:),intent( inout ) :: ct, c
        real, dimension(:,:),intent( inout ) :: INV_SIGMA, INV_SIGMA_NANO
        type(vector_field), intent( inout ) :: CV_RHS_field, CT_RHS
        real, dimension(:),intent( inout ) :: MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE ! of length NCMC
        real, dimension(:),intent( inout ) :: mass_pipe, MASS_PIPE_FOR_COUP ! of length cv_nonods

        logical, intent( in ) :: getcv_disc, getct
        type(petsc_csr_matrix), intent(inout) :: petsc_acv


        ! Local variables
        INTEGER :: CV_NODI, CV_NODJ, IPHASE, COUNT, CV_SILOC, SELE, cv_iloc, cv_jloc, jphase
        INTEGER :: cv_ncorner, cv_lnloc, u_lnloc, i_indx, j_indx, ele, cv_gi, iloop, ICORNER, NPIPES, i

        integer, dimension(:), pointer :: cv_neigh_ptr
        integer, dimension(:), allocatable:: CV_LOC_CORNER, U_LOC_CORNER, CV_GL_LOC, CV_GL_GL, X_GL_GL, MAT_GL_GL, u_GL_LOC, u_GL_GL, &
            pipe_corner_nds1, pipe_corner_nds2
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: CV_MID_SIDE, U_MID_SIDE

        real, dimension(:), allocatable:: xi_limit, cvweigh, cv_nodpos, u_nodpos, lcv_b, direction, direction_NORM, &
            PIPE_DIAM_GI_vol, PIPE_DIAM_GI_suf, suf_detwei, vol_detwei, vol_detwei2, &
            TUPWIND_OUT,  DUPWIND_OUT, TUPWIND_in,  DUPWIND_in, ndotq, income, income_j, &
            FEMTGI, FEMdGI, T_CV_NODI, T_CV_NODJ, D_CV_NODI, D_CV_NODJ, limt, limd, bczero, fvt, limdt, LOC_CT_RHS_U_ILOC, INV_SIGMA_GI
        real, dimension(:), allocatable:: CVN_VOL_ADJ
        INTEGER, dimension(:), allocatable:: WIC_B_BC_ALL_NODS
        INTEGER, dimension(:,:), allocatable:: WIC_T_BC_ALL_NODS, WIC_D_BC_ALL_NODS, WIC_U_BC_ALL_NODS
        real, dimension(:,:), allocatable:: cvn, n, nlx, un, unlx, sbcvfen, sbcvfenslx, sbufen, L_CVFENX_ALL, L_UFENX_ALL,L_UFEN_REVERSED,L_UFEN,&
            UGI_ALL, ct_con, x_all_corn
        real, dimension(:,:), allocatable:: cvn_fem, cvn_femlx, cvnn_fem, cvnn_femlx
        real, dimension(:,:), allocatable:: SUF_T_BC_ALL_NODS, SUF_D_BC_ALL_NODS, RVEC_SUM_T, RVEC_SUM_D, RVEC_SUM_U, MEAN_U
        real, dimension(:,:,:), allocatable:: SUF_U_BC_ALL_NODS
        real, dimension(:,:,:), allocatable:: L_CVFENX_ALL_REVERSED
        logical :: CV_QUADRATIC, U_QUADRATIC, ndiff, diff, PIPE_INDEX_LOGICAL(ndim+1), ELE_HAS_PIPE, integrate_other_side_and_not_boundary
        logical :: UPWIND_PIPES, PIPE_MIN_DIAM, IGNORE_DIAGONAL_PIPES, SOLVE_ACTUAL_VEL, LUMP_COUPLING_RES_PIPES, CALC_SIGMA_PIPE, get_c_pipes
        real :: LOC_CV_RHS_I(NPHASE)
        real :: T1(NDIM), T2(NDIM), TT1(NDIM), TT2(NDIM), NN1(NDIM), T1TT1, T1TT2, T2TT1, T2TT2, DET_SQRT, INV_SIGMA_ND, N1NN1, INV_SIGMA_NANO_ND

        real :: cv_ldx, u_ldx, dx, ele_angle, cv_m, sigma_gi, M_CVFEM2PIPE, M_PIPE2CVFEM, rnorm_sign, suf_area, PIPE_DIAM_END, INFINY, MIN_DIAM
        real :: MIN_INV_SIG, R1(1), R2(1), RZ(1)
        real :: TMAX(NPHASE), TMIN(NPHASE), DENMAX(NPHASE), DENMIN(NPHASE)
        integer :: ierr, PIPE_NOD_COUNT, NPIPES_IN_ELE, ipipe, CV_LILOC, CV_LJLOC, U_LILOC, &
            u_iloc, x_iloc, cv_knod, idim, cv_lkloc, u_lkloc, u_knod, gi, ncorner, cv_lngi, u_lngi, cv_bngi, bgi, &
            icorner1, icorner2, icorner3, icorner4, WIC_B_BC_DIRICHLET, JCV_NOD1, JCV_NOD2, CV_NOD, JCV_NOD, JU_NOD, &
            U_NOD, U_SILOC, COUNT2, MAT_KNOD, MAT_NODI, COUNT3, IPRES

        real, dimension(:,:), allocatable:: tmax_all, tmin_all, denmax_all, denmin_all

        type(tensor_field), pointer :: t_all, den_all, u_all
        type(scalar_field), pointer :: pipe_diameter
        type(vector_field), pointer :: X

        integrate_other_side_and_not_boundary = .FALSE.
        UPWIND_PIPES=.FALSE. ! Used for testing...
        PIPE_MIN_DIAM=.TRUE. ! Use the pipe min diamter along a pipe element edge and min inv_sigma (max. drag reflcting min pipe diameter)
        WIC_B_BC_DIRICHLET = 1
        INFINY=1.0E+20
        IGNORE_DIAGONAL_PIPES = .TRUE.
        SOLVE_ACTUAL_VEL = .TRUE. ! Solve for the actual real velocity in the pipes.
        LUMP_COUPLING_RES_PIPES = .TRUE. ! Lump the coupling term which couples the pressure between the pipe and reservior.
        CALC_SIGMA_PIPE = have_option("/porous_media/well_options/calculate_sigma_pipe") ! Calculate sigma based on friction factors...


        NCORNER = NDIM + 1

        allocate( xi_limit(nphase), &
            CV_LOC_CORNER(NCORNER), CV_MID_SIDE(NCORNER, NCORNER), &
            U_LOC_CORNER(NCORNER), U_MID_SIDE(NCORNER,NCORNER) )

        ! default limiting NVD diagram...
        XI_LIMIT(:) = 2.0 ; ndiff=.false. ; diff=.true.

        CV_QUADRATIC = (CV_NLOC==6.AND.NDIM==2) .OR. (CV_NLOC==10.AND.NDIM==3)
        U_QUADRATIC = (U_NLOC==6.AND.NDIM==2) .OR. (U_NLOC==10.AND.NDIM==3)

        ! Calculate the local corner nodes...
        CALL CALC_CORNER_NODS( CV_LOC_CORNER, NDIM, CV_NLOC, CV_QUADRATIC, CV_MID_SIDE )
        CALL CALC_CORNER_NODS( U_LOC_CORNER, NDIM, U_NLOC, U_QUADRATIC, U_MID_SIDE )

        if ( CV_QUADRATIC ) then
            cv_lngi = 9 !3
            cv_lnloc = 3
            cv_bngi = 2
        else
            cv_lngi = 6 !1
            cv_lnloc = 2
            cv_bngi = 1
        end if

        if ( U_QUADRATIC ) then
            u_lnloc = 3
        else
            u_lnloc = 2
        end if

        u_lngi = cv_lngi

        allocate( cvweigh(cv_lngi), cvn(cv_lnloc, cv_lngi), n(cv_lnloc, cv_lngi), &
            nlx(cv_lnloc, cv_lngi), un(u_lnloc, u_lngi), unlx(u_lnloc, u_lngi), &
            cv_nodpos(cv_lnloc), u_nodpos(u_lnloc), lcv_b(cv_lnloc), sbcvfen(cv_lnloc, cv_bngi), sbufen(u_lnloc, cv_bngi) )

        allocate(cvn_fem(cv_lnloc, cv_lngi), cvn_femlx(cv_lnloc, cv_lngi), cvnn_fem(cv_lnloc, cv_lngi), cvnn_femlx(cv_lnloc, cv_lngi))

        ! calculate shape functions...
        call quad_1d_shape( cv_lngi, cv_lnloc, u_lnloc, cvn, cvweigh, n, nlx, un, unlx )
        call quad_1d_shape( cv_lngi, cv_lnloc, cv_lnloc, cvn, cvweigh, cvn_fem, cvn_femlx, cvnn_fem, cvnn_femlx )

        ! calculate sbcvfen.
        ! determine the CV node positions...
        cv_nodpos(1) = -1.0
        cv_ldx = 2.0 / real(cv_lnloc-1)
        do i = 2, cv_lnloc
            cv_nodpos(i) = cv_nodpos(i-1) + cv_ldx
        end do

        ! for calculating sbufen.
        ! determine the U node positions...
        u_nodpos(1) = -1.0
        u_ldx = 2.0 / real(u_lnloc-1)
        do i = 2, u_lnloc
            u_nodpos(i) = u_nodpos(i-1) + u_ldx
        end do

        do i = 2, cv_lnloc
            lcv_b(i-1) = 0.5 * ( cv_nodpos(i-1) + cv_nodpos(i) )
        end do

        do CV_Liloc = 1, cv_lnloc
            do bgi = 1, cv_bngi
                SBCVFEN( cv_liloc, bgi ) = lagran( ndiff, lcv_b(bgi), cv_liloc, cv_lnloc, cv_nodpos )
            end do
        end do

        do U_Liloc = 1, U_lnloc
            do bgi = 1, cv_bngi
                SBUFEN( u_liloc, bgi ) = lagran( ndiff, lcv_b(bgi), U_liloc, u_lnloc, u_nodpos )
            end do
        end do

        ! Adjust DX for volume integrations...
        allocate(CVN_VOL_ADJ(CV_LNLOC))
        CVN_VOL_ADJ=0.5
        IF(CV_LNLOC==3) THEN
            CVN_VOL_ADJ(1)=0.25
            CVN_VOL_ADJ(2)=0.5
            CVN_VOL_ADJ(3)=0.25
        ENDIF


        ! SET UP THE SURFACE B.C'S
        allocate(WIC_B_BC_ALL_NODS(CV_NONODS))
        allocate(WIC_T_BC_ALL_NODS(NPHASE,CV_NONODS))
        allocate(WIC_D_BC_ALL_NODS(NPHASE,CV_NONODS))
        allocate(WIC_U_BC_ALL_NODS(NPHASE,CV_NONODS))

        allocate(SUF_T_BC_ALL_NODS(NPHASE,CV_NONODS))
        allocate(SUF_D_BC_ALL_NODS(NPHASE,CV_NONODS))
        allocate(SUF_U_BC_ALL_NODS(NDIM,NPHASE,CV_NONODS))

        allocate(RVEC_SUM_T(NPHASE,CV_NONODS))
        allocate(RVEC_SUM_D(NPHASE,CV_NONODS))
        allocate(RVEC_SUM_U(NPHASE,CV_NONODS))

        allocate(MEAN_U(NDIM,NPHASE))

        WIC_B_BC_ALL_NODS=0
        WIC_T_BC_ALL_NODS=0
        WIC_D_BC_ALL_NODS=0
        WIC_U_BC_ALL_NODS=0

        RVEC_SUM_T=0.0
        RVEC_SUM_D=0.0
        RVEC_SUM_U=0.0
        SUF_T_BC_ALL_NODS=0.0
        SUF_D_BC_ALL_NODS=0.0
        SUF_U_BC_ALL_NODS=0.0

        DO SELE = 1, STOTEL
            DO CV_SILOC = 1, CV_SNLOC
                CV_NOD = CV_SNDGLN( (SELE-1)*CV_SNLOC + CV_SILOC )
                WIC_B_BC_ALL_NODS( CV_NOD ) = WIC_B_BC_DIRICHLET

                DO IPHASE = N_IN_PRES+1, NPHASE
                    IF ( WIC_T_BC_ALL( 1,IPHASE, SELE ) == WIC_T_BC_DIRICHLET ) THEN
                        WIC_T_BC_ALL_NODS( IPHASE, CV_NOD ) = WIC_T_BC_DIRICHLET
                        SUF_T_BC_ALL_NODS( IPHASE, CV_NOD ) = SUF_T_BC_ALL_NODS( IPHASE, CV_NOD ) + SUF_T_BC_ALL( 1,IPHASE, (SELE-1)*CV_SNLOC + CV_SILOC )
                        RVEC_SUM_T( IPHASE, CV_NOD ) = RVEC_SUM_T( IPHASE, CV_NOD ) + 1.0
                    END IF
                    IF ( WIC_D_BC_ALL( 1,IPHASE, SELE ) == WIC_D_BC_DIRICHLET ) THEN
                        WIC_D_BC_ALL_NODS( IPHASE, CV_NOD ) = WIC_D_BC_DIRICHLET
                        SUF_D_BC_ALL_NODS( IPHASE, CV_NOD ) = SUF_D_BC_ALL_NODS( IPHASE, CV_NOD ) + SUF_D_BC_ALL( 1,IPHASE, (SELE-1)*CV_SNLOC + CV_SILOC )
                        RVEC_SUM_D( IPHASE, CV_NOD ) = RVEC_SUM_D( IPHASE, CV_NOD ) + 1.0
                    END IF
                END DO ! ENDOF DO IPHASE=N_IN_PRES+1,NPHASE
            END DO ! DO CV_SILOC=1,CV_SNLOC

            MEAN_U = 0.0
            DO U_SILOC = 1, U_SNLOC
                DO IPHASE = N_IN_PRES+1, NPHASE
                    IF ( WIC_U_BC_ALL( 1, IPHASE, SELE ) == WIC_U_BC_DIRICHLET ) THEN
                        U_NOD = U_SNDGLN( (SELE-1)*U_SNLOC + U_SILOC )
                        MEAN_U(:,IPHASE) = MEAN_U(:,IPHASE) + SUF_U_BC_ALL( :,IPHASE, (SELE-1)*U_SNLOC + U_SILOC ) / REAL(U_SNLOC)
                    END IF
                END DO ! ENDOF DO IPHASE=N_IN_PRES+1,NPHASE
            END DO ! DO U_SILOC=1,U_SNLOC

            DO CV_SILOC = 1, CV_SNLOC
                DO IPHASE = N_IN_PRES+1, NPHASE
                    IF ( WIC_U_BC_ALL( 1, IPHASE, SELE ) == WIC_U_BC_DIRICHLET ) THEN
                        CV_NOD = CV_SNDGLN( (SELE-1)*CV_SNLOC + CV_SILOC )
                        WIC_U_BC_ALL_NODS( IPHASE, CV_NOD ) = WIC_U_BC_DIRICHLET
                        SUF_U_BC_ALL_NODS( :,IPHASE, CV_NOD ) = SUF_U_BC_ALL_NODS( :,IPHASE, CV_NOD ) + MEAN_U(:,IPHASE)
                        RVEC_SUM_U( IPHASE, CV_NOD ) = RVEC_SUM_U( IPHASE, CV_NOD ) + 1.0
                    END IF
                END DO ! ENDOF DO IPHASE=N_IN_PRES+1,NPHASE
            END DO ! DO CV_SILOC=1,CV_SNLOC
        END DO

        DO CV_NOD = 1, CV_NONODS
            DO IPHASE = N_IN_PRES+1, NPHASE
                IF ( WIC_T_BC_ALL_NODS( IPHASE, CV_NOD ) /= 0 ) SUF_T_BC_ALL_NODS( IPHASE,CV_NOD ) = SUF_T_BC_ALL_NODS( IPHASE,CV_NOD ) / RVEC_SUM_T( IPHASE,CV_NOD )
                IF ( WIC_D_BC_ALL_NODS( IPHASE, CV_NOD ) /= 0 ) SUF_D_BC_ALL_NODS( IPHASE,CV_NOD ) = SUF_D_BC_ALL_NODS( IPHASE,CV_NOD ) / RVEC_SUM_D( IPHASE,CV_NOD )
                IF ( WIC_U_BC_ALL_NODS( IPHASE, CV_NOD ) /= 0 ) SUF_U_BC_ALL_NODS( :,IPHASE,CV_NOD ) = SUF_U_BC_ALL_NODS( :,IPHASE,CV_NOD ) / RVEC_SUM_U( IPHASE,CV_NOD )
            END DO
        END DO
        ! END OF SET UP THE SURFACE B.C'S


        allocate( tmax_all(nphase, cv_nonods), tmin_all(nphase, cv_nonods), &
            denmax_all(nphase, cv_nonods), denmin_all(nphase, cv_nonods) )

        T_ALL => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
        DEN_ALL => extract_tensor_field( packed_state, "PackedDensity" )
        U_ALL => extract_tensor_field( packed_state, "PackedVelocity" )
        PIPE_DIAMETER => extract_scalar_field( state(1), "DiameterPipe1" )
        X => EXTRACT_VECTOR_FIELD( PACKED_STATE, "PressureCoordinate" )

        TMAX_ALL=-INFINY; TMIN_ALL=+INFINY; DENMAX_ALL=-INFINY; DENMIN_ALL = +INFINY
        DO CV_NODI = 1, CV_NONODS
            IF ( PIPE_DIAMETER%VAL(CV_NODI) /= 0.0 ) THEN
                do count = finacv(cv_nodi), finacv(cv_nodi+1)-1
                    cv_nodj = colacv(count)
                    IF ( PIPE_DIAMETER%VAL(CV_NODJ) /= 0.0 ) THEN
                        DO IPHASE = 1, NPHASE
                            TMAX_ALL( IPHASE, CV_NODI ) = max( TMAX_ALL( IPHASE, CV_NODI ),T_ALL%val( 1, IPHASE, cv_nodj ) )
                            TMIN_ALL( IPHASE, CV_NODI ) = min( TMIN_ALL( IPHASE, CV_NODI ),T_ALL%val( 1, IPHASE, cv_nodj ) )
                            DENMAX_ALL( IPHASE, CV_NODI ) = max( DENMAX_ALL( IPHASE, CV_NODI ),DEN_ALL%val( 1, IPHASE, cv_nodj ) )
                            DENMIN_ALL( IPHASE, CV_NODI ) = min( DENMIN_ALL( IPHASE, CV_NODI ),DEN_ALL%val( 1, IPHASE, cv_nodj ) )
                        END DO
                    END IF
                END DO
            END IF
        END DO

        IF ( GETCV_DISC ) THEN
            ! Reset the value of the matrix to 0.0 for these phases...
            do iphase = n_in_pres+1, nphase
                do cv_nodi = 1, cv_nonods
                    do count = finacv(cv_nodi), finacv(cv_nodi+1)-1
                        cv_nodj = colacv(count)
                        i_indx = petsc_acv%row_numbering%gnn2unn( cv_nodi, iphase )
                        j_indx = petsc_acv%column_numbering%gnn2unn( cv_nodj, iphase )
                        call MatSetValue( petsc_acv, i_indx, j_indx, 0.0, INSERT_VALUES, ierr )
                    end do
                end do
            end do
            CV_RHS_field%val( N_IN_PRES+1:NPHASE, : ) = 0.0
        END IF

        GET_C_PIPES = .FALSE.

        IF ( GETCT ) THEN
            CT( :, N_IN_PRES+1:NPHASE, : ) = 0.0
            CT_RHS%val( 2:NPRES, : ) = 0.0
            IF ( GET_C_PIPES ) C( :, N_IN_PRES+1:NPHASE, : ) = 0.0
        END IF


        allocate( CV_GL_LOC( cv_lnloc ), CV_GL_GL( cv_lnloc ), X_GL_GL( cv_lnloc ), MAT_GL_GL( cv_lnloc ), &
            U_GL_LOC( u_lnloc ), U_GL_GL( u_lnloc ) )

        allocate( pipe_corner_nds1(ndim), pipe_corner_nds2(ndim), direction(ndim), direction_NORM(ndim)  )
        allocate( L_CVFENX_ALL(cv_lnloc, cv_lngi), L_UFENX_ALL(u_lnloc, cv_lngi) , PIPE_DIAM_GI_vol(cv_lngi) )
        allocate( PIPE_DIAM_GI_suf(cv_bngi) )
        allocate( L_CVFENX_ALL_REVERSED(ndim, cv_lnloc, cv_lngi), L_UFEN_REVERSED(cv_lngi, u_lnloc), L_UFEN(u_lnloc, cv_lngi) )
        allocate( suf_detwei(cv_bngi), vol_detwei(cv_lngi), vol_detwei2(cv_lngi), INV_SIGMA_GI(nphase) )

        allocate( TUPWIND_OUT(nphase), DUPWIND_OUT(nphase), TUPWIND_in(nphase), DUPWIND_in(nphase) )
        allocate( UGI_ALL(ndim, nphase), ndotq(nphase), income(nphase), income_j(nphase) )

        allocate( FEMTGI(nphase), FEMdGI(nphase), T_CV_NODI(nphase), T_CV_NODJ(nphase), D_CV_NODI(nphase), D_CV_NODJ(nphase) )
        allocate( limt(nphase), limd(nphase) )
        allocate( ct_con(ndim, nphase), bczero(nphase), fvt(nphase), limdt(nphase), LOC_CT_RHS_U_ILOC(nphase) )

        allocate( x_all_corn(ndim, NCORNER) )

        mass_pipe = 0.0; MASS_PIPE_FOR_COUP = 0.0
        MASS_CVFEM2PIPE = 0.0; MASS_PIPE2CVFEM = 0.0; MASS_CVFEM2PIPE_TRUE = 0.0

        INV_SIGMA = 0.0; INV_SIGMA_NANO = 0.0


        DO ELE = 1, TOTELE

            ! Look for pipe indicator in element:
            PIPE_INDEX_LOGICAL=.FALSE.
            PIPE_NOD_COUNT=0
            DO ICORNER = 1, NCORNER ! Number of corner nodes in element
                CV_ILOC = CV_LOC_CORNER( ICORNER )
                CV_NODI = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )
                IF ( PIPE_diameter%val( CV_NODI ) > 0.0 ) THEN
                    PIPE_NOD_COUNT = PIPE_NOD_COUNT + 1
                    PIPE_INDEX_LOGICAL(ICORNER) = .TRUE.
                END IF
            END DO
            NPIPES = PIPE_NOD_COUNT-1
            ELE_HAS_PIPE = NPIPES >= 1

            IF ( ELE_HAS_PIPE ) THEN

                X_ALL_CORN(:,1:NCORNER) = x%val(:, X_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_LOC_CORNER(1:NCORNER)) )

                ! Calculate the pipes within an element...
                ! we return the pipe corner nodes for each pipe in
                ! If we have more than one pipe then choose the 2 edges with the shortest sides
                ! and have a maximum of 2 pipes per element...
                CALL CALC_PIPES_IN_ELE( X_ALL_CORN, PIPE_INDEX_LOGICAL, NDIM, &
                    pipe_corner_nds1, pipe_corner_nds2, npipes )


                DO IPIPE = 1, NPIPES

                    ! DEFINE CV_LILOC:
                    CV_LILOC = 1
                    ICORNER = pipe_corner_nds1( IPIPE )
                    ICORNER1 = ICORNER
                    CV_ILOC = CV_LOC_CORNER( ICORNER )
                    CV_GL_LOC( CV_LILOC ) = CV_ILOC
                    CV_LILOC = CV_LNLOC
                    ICORNER = pipe_corner_nds2( IPIPE )
                    ICORNER2 = ICORNER
                    CV_ILOC = CV_LOC_CORNER( ICORNER )
                    CV_GL_LOC( CV_LILOC ) = CV_ILOC

                    IF ( CV_QUADRATIC ) CV_GL_LOC( 2 ) = CV_MID_SIDE( ICORNER1, ICORNER2 )

                    DO CV_LILOC = 1, CV_LNLOC
                        CV_ILOC = CV_GL_LOC( CV_LILOC )
                        X_ILOC = CV_ILOC
                        CV_GL_GL( CV_LILOC ) = CV_NDGLN( (ELE-1)*CV_NLOC + CV_ILOC )
                        X_GL_GL( CV_LILOC ) = X_NDGLN( (ELE-1)*CV_NLOC + X_ILOC )
                        MAT_GL_GL( CV_LILOC ) = MAT_NDGLN( (ELE-1)*CV_NLOC + CV_ILOC )
                    END DO

                    ! DEFINE U_LILOC:
                    U_LILOC = 1
                    ICORNER = pipe_corner_nds1( IPIPE )
                    U_ILOC = U_LOC_CORNER( ICORNER )
                    U_GL_LOC( U_LILOC ) = U_ILOC
                    U_LILOC = U_LNLOC
                    ICORNER = pipe_corner_nds2( IPIPE )
                    U_ILOC = U_LOC_CORNER( ICORNER )
                    U_GL_LOC( U_LILOC ) = U_ILOC

                    IF ( U_QUADRATIC ) U_GL_LOC( 2 ) = U_MID_SIDE( ICORNER1, ICORNER2 )

                    DO U_LILOC = 1, U_LNLOC
                        U_ILOC = U_GL_LOC( U_LILOC )
                        U_GL_GL( U_LILOC ) = U_NDGLN( (ELE-1)*U_NLOC + U_ILOC )
                    END DO

                    DIRECTION(:) = X_ALL_CORN( :, ICORNER2 ) - X_ALL_CORN( :, ICORNER1 )
                    DX = SQRT( SUM( DIRECTION(:)**2 ) )
                    DIRECTION(:) = DIRECTION(:) / DX

                    IF ( IGNORE_DIAGONAL_PIPES ) THEN
                        IF ( ABS(DIRECTION(1))<0.99 .AND. ABS(DIRECTION(2))<0.99.AND. ABS(DIRECTION(NDIM))<0.99 ) CYCLE
                    END IF

                    ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
                    L_CVFENX_ALL(:,:) = 2.0 * cvn_femlx(:,:) / DX
                    L_UFENX_ALL(:,:) = 2.0 * UNLX(:,:) / DX

                    IF(PIPE_MIN_DIAM) THEN
                        MIN_DIAM = MINVAL( PIPE_diameter%val( CV_GL_GL( : ) ) )
                        PIPE_DIAM_GI_VOL(:) = MIN_DIAM
                        PIPE_DIAM_GI_SUF(:) = MIN_DIAM
                    ELSE
                        PIPE_DIAM_GI_VOL(:) = 0.0
                        PIPE_DIAM_GI_SUF(:) = 0.0
                        DO CV_LILOC = 1, CV_LNLOC
                            CV_KNOD = CV_GL_GL( CV_LILOC )
                            PIPE_DIAM_GI_VOL(:) = PIPE_DIAM_GI_VOL(:) + PIPE_diameter%val( CV_KNOD ) * cvn_fem( CV_LILOC, : )
                            PIPE_DIAM_GI_SUF(:) = PIPE_DIAM_GI_SUF(:) + PIPE_diameter%val( CV_KNOD ) * SBCVFEN( CV_LILOC, : )
                        END DO
                        PIPE_DIAM_GI_VOL(:) = max( 0.0, PIPE_DIAM_GI_VOL(:) )
                        PIPE_DIAM_GI_SUF(:) = max( 0.0, PIPE_DIAM_GI_SUF(:) )
                    ENDIF



                    DO IDIM = 1, NDIM
                        L_CVFENX_ALL_REVERSED(IDIM,:,:) = cvn_femlx(:,:) * DIRECTION(IDIM)
                    END DO
                    DO U_LILOC = 1, U_LNLOC
                        L_UFEN_REVERSED( :, U_LILOC ) = un( U_LILOC, : )
                    END DO

                    ! Calculate element angle sweeped out by element and pipe
                    IF ( NDIM == 2 ) THEN
                        ELE_ANGLE = PI
                    ELSE
                        ! find the nodes other than the pipe end corner nodes...
                        ICORNER3 = 0
                        DO ICORNER = 1, NCORNER
                            IF ( ICORNER /= ICORNER1 ) THEN
                                IF ( ICORNER /= ICORNER2 ) THEN
                                    IF ( ICORNER3 == 0 ) THEN
                                        ICORNER3 = ICORNER
                                    ELSE
                                        ICORNER4 = ICORNER
                                    END IF
                                END IF
                            END IF
                        END DO
                        ELE_ANGLE = CALC_ELE_ANGLE_3D( X_ALL_CORN(:, ICORNER1 ), X_ALL_CORN(:, ICORNER2) , &
                            &                         X_ALL_CORN(:, ICORNER3 ), X_ALL_CORN(:, ICORNER4) )
                    END IF

                    ! Adjust according to the volume of the pipe...
                    SUF_DETWEI( : )  = 1.0               * PI * ( ( 0.5*PIPE_DIAM_GI_suf(:) )**2 ) * ELE_ANGLE / ( 2.0 * PI )
                    VOL_DETWEI( : )  = cvweigh( : ) * DX * PI * ( ( 0.5*PIPE_DIAM_GI_vol(:) )**2 ) * ELE_ANGLE / ( 2.0 * PI )
                    ! We do NOT multiply by r**2 here because we have not divided by r**2 in the mass transfer correlation (in cv_assemb)
                    VOL_DETWEI2( : ) = cvweigh( : ) * DX * PI *                                      ELE_ANGLE / ( 2.0 * PI )


                    DO IDIM = 1, NDIM
                        L_CVFENX_ALL_REVERSED( IDIM, :, : ) = L_CVFENX_ALL( :, : ) * DIRECTION( IDIM )
                    END DO
                    DO U_LILOC = 1, U_LNLOC
                        L_UFEN_REVERSED( :, U_LILOC ) =  UN( U_LILOC, : )
                    END DO

                    ! Add contributions from the volume...
                    IF( LUMP_COUPLING_RES_PIPES ) THEN
                        DO CV_LILOC = 1, CV_LNLOC
                            CV_NODI = CV_GL_GL(CV_LILOC)
                            MASS_PIPE(CV_NODI) = MASS_PIPE(CV_NODI) + SUM( CVN(CV_LILOC,:) *  CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI( : ) )
                            MASS_PIPE_FOR_COUP(CV_NODI) = MASS_PIPE_FOR_COUP(CV_NODI) + SUM( CVN(CV_LILOC,:) *  CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI2( : ) )
                            DO COUNT2 = FINDCMC(CV_NODI), FINDCMC(CV_NODI+1)-1
                                IF ( CV_NODI==colcmc(count2)) count=count2
                            end do
                            MASS_CVFEM2PIPE(COUNT) = MASS_CVFEM2PIPE(COUNT) + SUM( CVN(CV_LILOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI2( : ) )
                            MASS_PIPE2CVFEM(COUNT) = MASS_PIPE2CVFEM(COUNT) + SUM( CVN(CV_LILOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI2( : ) ) ! this is possibly bugged...
                            DO CV_LJLOC = 1, CV_LNLOC
                                CV_NODJ = CV_GL_GL(CV_LJLOC)
                                DO COUNT = FINDCMC(CV_NODI), FINDCMC(CV_NODI+1)-1
                                    IF ( CV_NODI==CV_NODJ ) THEN
                                        MASS_CVFEM2PIPE_TRUE(COUNT) = MASS_CVFEM2PIPE_TRUE(COUNT) + SUM( CVN(CV_LILOC,:) * CVN_FEM(CV_LJLOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI( : ) )
                                    END IF
                                END DO
                            END DO
                        END DO
                    ELSE
                        DO CV_LILOC = 1, CV_LNLOC
                            CV_NODI = CV_GL_GL(CV_LILOC)
                            MASS_PIPE(CV_NODI) = MASS_PIPE(CV_NODI) + SUM( CVN(CV_LILOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI( : ) )
                            DO CV_LJLOC = 1, CV_LNLOC
                                CV_NODJ = CV_GL_GL(CV_LJLOC)
                                DO COUNT = FINDCMC(CV_NODI), FINDCMC(CV_NODI+1)-1
                                    IF ( CV_NODI==CV_NODJ ) THEN
                                        MASS_CVFEM2PIPE(COUNT) = MASS_CVFEM2PIPE(COUNT) + SUM( CVN(CV_LILOC,:) * CVN_FEM(CV_LJLOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI2( : ) )
                                        MASS_PIPE2CVFEM(COUNT) = MASS_PIPE2CVFEM(COUNT) + SUM( CVN(CV_LJLOC,:) * CVN_FEM(CV_LILOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI2( : ) )
                                        MASS_CVFEM2PIPE_TRUE(COUNT) = MASS_CVFEM2PIPE_TRUE(COUNT) + SUM( CVN(CV_LJLOC,:) * CVN_FEM(CV_LILOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI( : ) )
                                    END IF
                                END DO
                            END DO
                        END DO
                    END IF

                    ! Calculate INV_SIGMA_APPROX for the 1st N_IN_PRES phases...

                    ! Determine the tangent and bi-normal vectors from the normal NormX, NormY, NormZ:
                    IF(NDIM==2) THEN
                        RZ=0.0
                        CALL GET_TANG_BINORM( DIRECTION(1), DIRECTION(2), RZ,           T1(1), T1(2), R1,    T2(1), T2(2), R2, 1 )
                    ELSE
                        CALL GET_TANG_BINORM( DIRECTION(1), DIRECTION(2), DIRECTION(3), T1(1), T1(2), T1(3), T2(1), T2(2), T2(3), 1 )
                    END IF

                    DO CV_LILOC = 1, CV_LNLOC
                        MAT_NODI = MAT_GL_GL(CV_LILOC)
                        CV_NODI = CV_GL_GL(CV_LILOC)

                        DO IPHASE = 1, N_IN_PRES
                            TT1(:) = MATMUL( OPT_VEL_UPWIND_COEFS_NEW(:,:,IPHASE, MAT_NODI), T1(:) )
                            T1TT1 = SUM(T1(:)*TT1(:))
                            IF ( NDIM==3 ) THEN
                                TT2(:) = MATMUL( OPT_VEL_UPWIND_COEFS_NEW(:,:,IPHASE, MAT_NODI), T2(:) )

                                T1TT2 = SUM( T1(:)*TT2(:) )
                                T2TT1 = SUM( T2(:)*TT1(:) )
                                T2TT2 = SUM( T2(:)*TT2(:) )
                            ELSE
                                T1TT2 = 0.0
                                T2TT1 = 0.0
                                T2TT2 = T1TT1
                            END IF
                            DET_SQRT = SQRT( ABS( T1TT1*T2TT2 - T1TT2*T2TT1 ) )
                            INV_SIGMA_ND = 1.0 / MAX( 1.E-25, DET_SQRT)
                            INV_SIGMA(IPHASE,CV_NODI) = INV_SIGMA(IPHASE,CV_NODI) + INV_SIGMA_ND * SUM( CVN(CV_LILOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI( : ) )
                            ! For the nano laterals...
                            NN1(:) = MATMUL( OPT_VEL_UPWIND_COEFS_NEW(:,:,IPHASE, MAT_NODI), DIRECTION(:) )
                            N1NN1 = SUM( DIRECTION(:)*NN1(:) )
                            INV_SIGMA_NANO_ND = 1.0/MAX(1.E-25, N1NN1)
                            INV_SIGMA_NANO(IPHASE,CV_NODI) = INV_SIGMA_NANO(IPHASE,CV_NODI) + INV_SIGMA_NANO_ND * SUM( CVN(CV_LILOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI( : ) )
                        END DO
                    END DO



                    ! The number of CV basis functions...
                    DO BGI = 1, cv_bngi
                        DO ILOOP = 1, 2

                            IF ( ILOOP==1 ) THEN
                                CV_LILOC = BGI
                                CV_LJLOC = BGI + 1
                                CV_NODI = CV_GL_GL(CV_LILOC)
                                CV_NODJ = CV_GL_GL(CV_LJLOC)
                                direction_norm = + direction
                            ELSE
                                CV_LILOC = BGI + 1
                                CV_LJLOC = BGI
                                CV_NODI = CV_GL_GL(CV_LILOC)
                                CV_NODJ = CV_GL_GL(CV_LJLOC)
                                direction_norm = - direction
                            END IF

                            TUPWIND_OUT=0.0; DUPWIND_OUT=0.0
                            TUPWIND_IN=0.0; DUPWIND_IN=0.0
                            DO IPHASE = N_IN_PRES+1, NPHASE
                                ! CV incomming T:
                                IF ( T_ALL%val( 1, IPHASE, CV_NODI ) > T_ALL%val( 1, IPHASE, CV_NODJ ) ) THEN
                                    TUPWIND_OUT( IPHASE ) = TMAX_ALL( IPHASE, CV_NODI )
                                ELSE
                                    TUPWIND_OUT( IPHASE ) = TMIN_ALL( IPHASE, CV_NODI )
                                END IF
                                IF ( DEN_ALL%val( 1, IPHASE, CV_NODI ) > DEN_ALL%val( 1, IPHASE, CV_NODJ ) ) THEN
                                    DUPWIND_OUT( IPHASE ) = DENMAX_ALL( IPHASE, CV_NODI )
                                ELSE
                                    DUPWIND_OUT( IPHASE ) = DENMIN_ALL( IPHASE, CV_NODI )
                                END IF

                                ! CV outgoing T:
                                IF ( T_ALL%val( 1, IPHASE, CV_NODI ) < T_ALL%val( 1, IPHASE, CV_NODJ ) ) THEN
                                    TUPWIND_IN( IPHASE ) = TMAX_ALL( IPHASE, CV_NODJ )
                                ELSE
                                    TUPWIND_IN( IPHASE ) = TMIN_ALL( IPHASE, CV_NODJ )
                                END IF
                                IF ( DEN_ALL%val( 1, IPHASE, CV_NODI ) < DEN_ALL%val( 1, IPHASE, CV_NODJ ) ) THEN
                                    DUPWIND_IN( IPHASE ) = DENMAX_ALL( IPHASE, CV_NODJ )
                                ELSE
                                    DUPWIND_IN( IPHASE ) = DENMIN_ALL( IPHASE, CV_NODJ )
                                END IF

                            END DO


                            ! Value of sigma in the force balance eqn...
                            IF ( SOLVE_ACTUAL_VEL ) THEN
                                INV_SIGMA_GI = 1.0
                            ELSE
                                IF ( PIPE_MIN_DIAM ) THEN
                                    DO IPHASE = 1, NPHASE
                                        MIN_INV_SIG = MINVAL( INV_SIGMA(IPHASE,MAT_GL_GL( : ) ) )
                                        INV_SIGMA_GI(IPHASE) = MIN_INV_SIG
                                    END DO
                                ELSE
                                    INV_SIGMA_GI(:) = 0.0
                                    DO CV_LKLOC = 1, CV_LNLOC
                                        MAT_KNOD = MAT_GL_GL(CV_LKLOC)
                                        DO IPHASE = 1, NPHASE
                                            INV_SIGMA_GI(IPHASE) = INV_SIGMA_GI(IPHASE) + SBCVFEN( CV_LKLOC, BGI ) * INV_SIGMA(IPHASE,MAT_KNOD)
                                        END DO
                                    END DO
                                END IF
                                IF(CALC_SIGMA_PIPE) THEN
                                    STOP 'OPTION NOT READY YET AS WE NEED TO CALCULATE INV_SIGMA_GI WHICH IS A FUNCTION OF VELOCITY'
                                END IF
                            END IF


                            ! Velocity in the pipe
                            UGI_ALL(:,:) = 0.0
                            DO U_LKLOC = 1, U_LNLOC
                                U_KNOD = U_GL_GL(U_LKLOC)
                                UGI_ALL(:,:) = UGI_ALL(:,:) + SBUFEN( u_lkloc, bgi ) * U_ALL%val(:,:,u_kNOD)
                            END DO
                            DO IPHASE=1,NPHASE
                                UGI_ALL(:,IPHASE) = UGI_ALL(:,IPHASE) * INV_SIGMA_GI(IPHASE)
                            END DO

                            DO IPHASE = 1, NPHASE
                                NDOTQ(IPHASE) = SUM( DIRECTION_norm(:) * UGI_all(:,IPHASE) )
                            END DO

                            ! When NDOTQ == 0, INCOME_J has to be 1 as well, not 0
                            WHERE ( NDOTQ <= 0.0 )
                                INCOME_J = 0.0
                            ELSE WHERE
                                INCOME_J = 1.0
                            END WHERE
                            INCOME = 1.0 - INCOME_J


                            ! high order values...
                            FEMTGI=0.0 ; FEMDGI=0.0
                            DO CV_LKLOC = 1, CV_LNLOC
                                CV_KNOD = CV_GL_GL(CV_LKLOC)
                                FEMTGI(:) = FEMTGI(:) + CVN_FEM(CV_LJLOC,BGI) * T_ALL%val( 1, :, CV_KNOD)
                                FEMDGI(:) = FEMDGI(:) + CVN_FEM(CV_LJLOC,BGI) * DEN_ALL%val( 1, :, CV_KNOD)
                            END DO
                            FEMDGI(:) = max( 0.0,FEMDGI(:) )

                            T_CV_NODI(:) = T_ALL%val( 1, :, CV_NODI)
                            T_CV_NODJ(:) = T_ALL%val( 1, :, CV_NODJ)
                            D_CV_NODI(:) = DEN_ALL%val( 1, :, CV_NODI)
                            D_CV_NODJ(:) = DEN_ALL%val( 1, :, CV_NODJ)

                            IF ( UPWIND_PIPES ) THEN ! Used for testing...
                                LIMT(:) = T_CV_NODI(:)*(1.0-INCOME(:)) + T_CV_NODJ(:)*INCOME(:)
                                LIMD(:) = D_CV_NODI(:)*(1.0-INCOME(:)) + D_CV_NODJ(:)*INCOME(:)
                            ELSE
                                ! Call the limiter for T...
                                CALL ONVDLIM_ANO_MANY( NPHASE, &
                                    LIMT(:), FEMTGI(:), INCOME(:), &
                                    T_CV_NODI(:), T_CV_NODJ(:), XI_LIMIT(:), &
                                    TUPWIND_IN(:), TUPWIND_OUT(:) )

                                ! Call the limiter for D...
                                CALL ONVDLIM_ANO_MANY( NPHASE, &
                                    LIMD(:), FEMDGI(:), INCOME(:), &
                                    D_CV_NODI(:), D_CV_NODJ(:), XI_LIMIT(:), &
                                    DUPWIND_IN(:), DUPWIND_OUT(:) )
                            END IF

                            LIMDT(:) = LIMD(:) * LIMT(:)

                            IF ( GETCT ) THEN ! Obtain the CV discretised CT eqations plus RHS
                                DO U_LKLOC = 1, U_LNLOC
                                    U_KNOD = U_GL_GL(U_LKLOC)
                                    DO IDIM = 1, NDIM
                                        CT_CON(IDIM,:) = SBUFEN( U_LKLOC, BGI ) * LIMDT(:) * suf_DETWEI( BGI ) * DIRECTION_norm(IDIM) * INV_SIGMA_GI(:) / D_CV_NODI(:)
                                    END DO
                                    ! Put into CT matrix...
                                    DO COUNT = FINDCT(CV_NODI), FINDCT(CV_NODI+1)-1
                                        IF ( COLCT(COUNT)==U_KNOD ) CT( :, N_IN_PRES+1:NPHASE, COUNT ) = &
                                            CT( :, N_IN_PRES+1:NPHASE, COUNT ) + CT_CON( :, N_IN_PRES+1:NPHASE )
                                    END DO

                                    IF ( GET_C_PIPES ) THEN
                                        DO COUNT = FINDC(U_KNOD), FINDC(U_KNOD+1)-1
                                            IF ( COLC(COUNT)==CV_NODI )  THEN
                                                DO IDIM = 1, NDIM
                                                    C( IDIM, N_IN_PRES+1:NPHASE, COUNT ) = C( IDIM, N_IN_PRES+1:NPHASE, COUNT ) + &
                                                        SBUFEN( U_LKLOC, BGI ) * SUF_DETWEI( BGI ) * DIRECTION_NORM(IDIM)
                                                END DO
                                            END IF
                                        END DO
                                    END IF

                                END DO
                            END IF

                            IF ( GETCV_DISC ) THEN

                                FVT(:) = T_CV_NODI(:)*(1.0-INCOME(:)) + T_CV_NODJ(:)*INCOME(:)

                                ! Put results into the RHS vector
                                LOC_CV_RHS_I = 0.0
                                do iphase = n_in_pres+1, nphase
                                    LOC_CV_RHS_I( IPHASE ) =  LOC_CV_RHS_I( IPHASE ) &
                                          ! subtract 1st order adv. soln.
                                        + suf_DETWEI( bGI ) * NDOTQ(IPHASE) * LIMD(IPHASE) * FVT(IPHASE) &
                                        - suf_DETWEI( bGI ) * NDOTQ(IPHASE) * LIMDT(IPHASE) ! hi order adv
                                end do

                                ! Put into matrix...
                                do iphase = n_in_pres+1, nphase
                                    call addto( petsc_acv, iphase, iphase, cv_nodi, cv_nodi, &
                                        +suf_DETWEI( bGI ) * NDOTQ(iphase) * ( 1. - INCOME(iphase) ) * LIMD(iphase) )
                                    call addto( petsc_acv, iphase, iphase, cv_nodi, cv_nodj, &
                                        +suf_DETWEI( bGI ) * NDOTQ(iphase) * INCOME(iphase) * LIMD(iphase) )
                                end do

                                call addto( cv_rhs_field, CV_NODI, LOC_CV_RHS_I )

                            END IF

                        END DO ! DO ILOOP = 1, 2
                    END DO ! DO BGI = 1, CV_BNGI


                    ! Add velocity and saturation b.c.'s to matrix and rhs...
                    JCV_NOD1 = CV_NDGLN( (ELE-1)*CV_NLOC + CV_LOC_CORNER( ICORNER1 ) )
                    JCV_NOD2 = CV_NDGLN( (ELE-1)*CV_NLOC + CV_LOC_CORNER( ICORNER2 ) )
                    JCV_NOD=0

                    IF ( WIC_B_BC_ALL_NODS( JCV_NOD1 ) == WIC_B_BC_DIRICHLET ) THEN
                        CV_LILOC = 1
                        JCV_NOD = JCV_NOD1
                        U_LILOC = 1
                        JU_NOD = U_GL_GL( U_LILOC )
                        direction_norm = - direction ! for the b.c it must be -ve at the bottom of element
                    END IF
                    IF ( WIC_B_BC_ALL_NODS( JCV_NOD2 ) == WIC_B_BC_DIRICHLET ) THEN
                        CV_LILOC = CV_LNLOC
                        JCV_NOD = JCV_NOD2
                        U_LILOC = U_LNLOC
                        JU_NOD = U_GL_GL( U_LILOC )
                        direction_norm = + direction ! for the b.c it must be +ve at the top of element
                    END IF

                    IF ( JCV_NOD /= 0 ) THEN

                        PIPE_DIAM_END = PIPE_diameter%val( JCV_NOD )
                        NDOTQ = 0.0
                        DO IPHASE = N_IN_PRES+1, NPHASE
                            IF ( SOLVE_ACTUAL_VEL ) THEN
                                NDOTQ(IPHASE) = SUM( direction_norm(:) * U_ALL%val(:,IPHASE,JU_NOD) )
                            ELSE
                                NDOTQ(IPHASE) = SUM( direction_norm(:) * U_ALL%val(:,IPHASE,JU_NOD) * INV_SIGMA_GI(IPHASE) )
                            END IF
                        END DO

                        INCOME(:) = 0.5 * ( 1. + SIGN( 1.0, -NDOTQ(:) ) )

                        LIMT(:)=0.0
                        LIMD(:)=0.0
                        FVT(:) =0.0
                        DO IPHASE = N_IN_PRES+1, NPHASE
                            IF ( WIC_T_BC_ALL_NODS( IPHASE, JCV_NOD ) == WIC_T_BC_DIRICHLET ) THEN
                                LIMT(IPHASE)=T_ALL%val(1,IPHASE,JCV_NOD)*(1.0-INCOME(IPHASE)) + SUF_T_BC_ALL_NODS(IPHASE,JCV_NOD)*INCOME(IPHASE)
                            ELSE
                                LIMT(IPHASE)=T_ALL%val(1,IPHASE,JCV_NOD)
                            END IF
                            FVT(IPHASE) = LIMT(IPHASE)
                        END DO
                        DO IPHASE = N_IN_PRES+1, NPHASE
                            IF ( WIC_D_BC_ALL_NODS( IPHASE, JCV_NOD ) == WIC_D_BC_DIRICHLET ) THEN
                                LIMD(IPHASE)=DEN_ALL%val(1,IPHASE,JCV_NOD)*(1.0-INCOME(IPHASE)) + SUF_D_BC_ALL_NODS(IPHASE,JCV_NOD)*INCOME(IPHASE)
                            ELSE
                                LIMD(IPHASE)=DEN_ALL%val(1,IPHASE,JCV_NOD)
                            END IF
                        END DO

                        LIMDT(:) = LIMD(:) * LIMT(:)
                        ! Add in C matrix contribution: (DG velocities)
                        ! In this section we multiply the shape functions over the GI points. i.e: we perform the integration
                        ! over the element of the pressure like source term.
                        ! Put into matrix

                        ! Prepare aid variable NMX_ALL to improve the speed of the calculations
                        suf_area = PI * ( (0.5*PIPE_DIAM_END)**2 ) * ELE_ANGLE / ( 2.0 * PI )

                        IF ( GETCT ) THEN ! Obtain the CV discretised CT eqations plus RHS on the boundary...
                            DO IDIM = 1, NDIM
                                CT_CON(IDIM,:) = LIMDT(:) * suf_area * DIRECTION_NORM(IDIM) * INV_SIGMA_GI(:) / DEN_ALL%val(1,:,JCV_NOD)
                            END DO
                            ! Put into CT matrix...
                            COUNT2=0
                            DO COUNT = FINDCT(JCV_NOD), FINDCT(JCV_NOD+1)-1
                                IF ( COLCT(COUNT)==JU_NOD ) COUNT2=COUNT
                            END DO

                            IF ( GET_C_PIPES ) THEN
                                COUNT3=0
                                DO COUNT = FINDC(JU_NOD), FINDC(JU_NOD+1)-1
                                    IF ( COLC(COUNT)==JCV_NOD ) COUNT3=COUNT
                                END DO
                            END IF

                            LOC_CT_RHS_U_ILOC = 0.0
                            DO IPHASE = N_IN_PRES+1, NPHASE
                                IF ( WIC_U_BC_ALL_NODS( IPHASE, JCV_NOD ) == WIC_U_BC_DIRICHLET ) THEN
                                    DO IDIM = 1, NDIM
                                        LOC_CT_RHS_U_ILOC(IPHASE) = LOC_CT_RHS_U_ILOC(IPHASE) &
                                            - CT_CON(IDIM,IPHASE) * SUF_U_BC_ALL_NODS(IDIM,IPHASE,JCV_NOD)
                                    END DO
                                ELSE
                                    CT( :, IPHASE, COUNT2 ) = CT( :, IPHASE, COUNT2 ) + CT_CON( :, IPHASE )
                                    IF ( GET_C_PIPES ) THEN
                                        C( :, IPHASE, COUNT3 ) = C( :, IPHASE, COUNT3 ) + suf_area * DIRECTION_NORM(:)
                                    END IF
                                END IF
                            END DO

                            DO IPRES = 2, NPRES
                                call addto( ct_rhs, ipres, jcv_nod, &
                                    sum( LOC_CT_RHS_U_ILOC( 1+(ipres-1)*n_in_pres : ipres*n_in_pres ) ) )
                            END DO

                        END IF ! IF ( GETCT ) THEN

                        IF ( GETCV_DISC ) THEN ! this is on the boundary...

                            ! Put results into the RHS vector
                            LOC_CV_RHS_I = 0.0
                            do iphase = n_in_pres+1, nphase
                                LOC_CV_RHS_I( IPHASE ) =  LOC_CV_RHS_I( IPHASE ) &
                                    ! subtract 1st order adv. soln.
                                    + suf_area * NDOTQ(IPHASE) * LIMD(IPHASE) * FVT(IPHASE) &
                                    - suf_area * NDOTQ(IPHASE) * LIMDT(IPHASE) ! hi order adv
                            end do

                            ! Put into matrix...
                            do iphase = n_in_pres+1, nphase
                                call addto( petsc_acv, iphase, iphase, JCV_NOD, JCV_NOD, &
                                    + suf_area * NDOTQ(iphase) * ( 1. - INCOME(iphase) ) * LIMD(iphase) )
                            end do

                            call addto( cv_rhs_field, JCV_NOD, LOC_CV_RHS_I )

                        ENDIF ! ENDOF IF ( GETCV_DISC ) THEN

                    ENDIF ! ENDOF IF(JCV_NOD.NE.0) THEN

                END DO ! DO IPIPE2 = 1, NPIPES_IN_ELE
            END IF ! IF ( ELE_HAS_PIPE ) THEN

        END DO ! DO ELE = 1, TOTELE

        DO IPHASE = 1, N_IN_PRES
            INV_SIGMA(IPHASE,:) = INV_SIGMA(IPHASE,:) / MAX( MASS_PIPE(:), 1.E-15 )
            INV_SIGMA_NANO(IPHASE,:) = INV_SIGMA_NANO(IPHASE,:) / MAX( MASS_PIPE(:), 1.E-15 )
            ! We divide by this so that we get the right source term for the nano laterals...
            INV_SIGMA_NANO(IPHASE,:) = INV_SIGMA_NANO(IPHASE,:) / MAX( MASS_PIPE_FOR_COUP(:), 1.E-15 )
        END DO

        IF ( GETCV_DISC ) THEN
            do iphase = n_in_pres+1, nphase
                do cv_nodi = 1, cv_nonods
                    if ( pipe_diameter%val(cv_nodi)==0.0 ) then
                        cv_nodj = cv_nodi ; jphase = iphase
                        i_indx = petsc_acv%row_numbering%gnn2unn( cv_nodi, iphase )
                        j_indx = petsc_acv%column_numbering%gnn2unn( cv_nodj, jphase )
                        call MatSetValue( petsc_acv, i_indx, j_indx, 1.0, INSERT_VALUES, ierr )
                    end if
                end do
            end do
        END IF


        RETURN


    CONTAINS


        PURE SUBROUTINE ONVDLIM_ANO_MANY( NFIELD, &
            TDLIM, TDCEN, INCOME, &
            ETDNEW_PELE, ETDNEW_PELEOT, XI_LIMIT, &
            TUPWIN, TUPWI2 )
            implicit none
            ! This sub calculates the limited face values TDADJ(1...SNGI) from the central
            ! difference face values TDCEN(1...SNGI) using a NVD shceme.
            ! INCOME(1...SNGI)=1 for incomming to element ELE  else =0.
            ! LIBETA is the flux limiting parameter.
            ! TDMAX(PELE)=maximum of the surrounding 6 element values of element PELE.
            ! TDMIN(PELE)=minimum of the surrounding 6 element values of element PELE.
            ! PELEOT=element at other side of current face.
            ! ELEOT2=element at other side of the element ELEOTH.
            ! ELESID=element next to oposing current face.
            ! The elements are arranged in this order: ELEOT2,ELE, PELEOT, ELESID.
            ! This sub finds the neighbouring elements. Suppose that this is the face IFACE.
            !---------------------------------------------------
            !|   ELEOT2   |   ELEOTH   |   ELE     |   ELESID   |
            !---------------------------------------------------
            ! TAIN         THALF       TAOUT
            !---------------------------------------------------
            !>TEXTIN
            !TEXOUT<
            !---------------------------------------------------
            INTEGER, intent( in ) :: NFIELD
            REAL, DIMENSION( NFIELD ), intent( inout ) :: TDLIM
            REAL, DIMENSION( NFIELD ), intent( in ) :: TDCEN, INCOME, XI_LIMIT, TUPWIN, TUPWI2
            REAL, DIMENSION( NFIELD ), intent( in ) :: ETDNEW_PELE, ETDNEW_PELEOT
            ! Local variables
            REAL, PARAMETER :: TOLER=1.0E-10
            REAL :: DENOIN(NFIELD), CTILIN(NFIELD), DENOOU(NFIELD), &
                CTILOU(NFIELD), FTILIN(NFIELD), FTILOU(NFIELD)


            ! Calculate normalisation parameters for incomming velocities

            DENOIN = ( ETDNEW_PELE - TUPWIN )
            where( ABS( DENOIN ) < TOLER )
                DENOIN = SIGN( TOLER, DENOIN )
            end where
            CTILIN = ( ETDNEW_PELEOT - TUPWIN ) / DENOIN

            ! Calculate normalisation parameters for out going velocities

            DENOOU = ( ETDNEW_PELEOT - TUPWI2 )
            where( ABS( DENOOU ) < TOLER )
                DENOOU = SIGN( TOLER, DENOOU )
            end where
            CTILOU = ( ETDNEW_PELE - TUPWI2 ) / DENOOU

            FTILIN = ( TDCEN - TUPWIN ) / DENOIN
            FTILOU = ( TDCEN - TUPWI2 ) / DENOOU

            ! Velocity is going out of element
            TDLIM =        INCOME   * ( TUPWIN + NVDFUNNEW_MANY( FTILIN, CTILIN, XI_LIMIT ) * DENOIN ) &
                + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW_MANY( FTILOU, CTILOU, XI_LIMIT ) * DENOOU )

            TDLIM = MAX( TDLIM, 0.0 )

            RETURN
        END SUBROUTINE ONVDLIM_ANO_MANY


    END SUBROUTINE MOD_1D_CT_AND_ADV





    SUBROUTINE MOD_1D_FORCE_BAL_C( STATE, packed_state, U_RHS, Mdims, Mspars, GET_PIVIT_MAT, &
        &                         C, CV_NDGLN, U_NDGLN, X_NDGLN, MAT_NDGLN, pivit_mat, &
        &                         P_SNDGLN, WIC_P_BC_ALL,SUF_P_BC_ALL, SIGMA, NU_ALL, &
        &                         U_SOURCE, U_SOURCE_CV )
        ! This sub modifies either CT or the advection-diffusion equation for 1D pipe modelling
        IMPLICIT NONE
        TYPE(STATE_TYPE), DIMENSION( : ), INTENT( IN ) :: STATE
        TYPE(STATE_TYPE), INTENT( IN ) :: packed_STATE
        type(multi_dimensions), intent(in) :: Mdims
        type (multi_sparsities), intent(in) :: Mspars
        REAL, DIMENSION( :, :, : ), INTENT( INOUT ) :: U_RHS, C, pivit_mat
        REAL, DIMENSION( :, :, : ), INTENT( IN ) :: SUF_P_BC_ALL, U_SOURCE, U_SOURCE_CV
        REAL, DIMENSION( :, : ), INTENT( IN ) :: SIGMA
        REAL, DIMENSION( :, :, : ), INTENT( IN ) :: NU_ALL
        INTEGER, DIMENSION( : ), INTENT( IN ) :: CV_NDGLN, U_NDGLN, X_NDGLN, MAT_NDGLN, P_SNDGLN
        INTEGER, DIMENSION( :,:,: ), INTENT( IN ) :: WIC_P_BC_ALL
        LOGICAL, INTENT( IN ) :: GET_PIVIT_MAT
        LOGICAL :: CV_QUADRATIC, U_QUADRATIC, ELE_HAS_PIPE, PIPE_MIN_DIAM, IGNORE_DIAGONAL_PIPES, SOLVE_ACTUAL_VEL
        LOGICAL :: CALC_SIGMA_PIPE, DEFAULT_SIGMA_PIPE_OPTIONS, SWITCH_PIPES_ON_AND_OFF
        INTEGER :: ELE, PIPE_NOD_COUNT, ICORNER, &
            &     CV_ILOC, U_ILOC, CV_NODI, IPIPE, CV_LILOC, U_LILOC, CV_LNLOC, U_LNLOC, CV_KNOD, MAT_KNOD, IDIM, &
            &     IU_NOD, P_LJLOC, JCV_NOD, COUNT, COUNT2, IPHASE
        INTEGER, DIMENSION(:), ALLOCATABLE :: CV_LOC_CORNER, U_LOC_CORNER, CV_GL_LOC, CV_GL_GL, X_GL_GL, MAT_GL_GL, U_GL_LOC, U_GL_GL
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: CV_MID_SIDE, U_MID_SIDE, WIC_P_BC_ALL_NODS
        TYPE(SCALAR_FIELD), POINTER :: PIPE_DIAMETER, WD
        TYPE(VECTOR_FIELD), POINTER :: X
        TYPE(TENSOR_FIELD), POINTER :: WM, CV_VOL_FRAC
        REAL, DIMENSION( :, : ), ALLOCATABLE :: scvfen, scvfenslx, scvfensly, &
            &                                  scvfenlx, scvfenly, scvfenlz, &
            &                                  sufen, sufenslx, sufensly, &
            &                                  sufenlx, sufenly, sufenlz
        REAL, DIMENSION( : ), ALLOCATABLE :: scvfeweigh
        REAL, DIMENSION( :, :, : ), ALLOCATABLE :: L_CVFENX_ALL_REVERSED
        REAL, DIMENSION( :, : ), ALLOCATABLE :: L_CVFENX_ALL, L_UFENX_ALL, L_UFEN_REVERSED
        REAL, DIMENSION( :, : ), ALLOCATABLE :: X_ALL_CORN, SUF_P_BC_ALL_NODS, RVEC_SUM, LOC_U_RHS_U_ILOC
        REAL, DIMENSION( : ), ALLOCATABLE :: DETWEI, PIPE_DIAM_GI, NMX_ALL, WELL_DENSITY, WELL_VISCOSITY
        REAL, DIMENSION( :, : ), ALLOCATABLE :: SIGMA_GI, SIGMA_ON_OFF_GI
        TYPE( SCALAR_FIELD ), POINTER :: PHASE_EXCLUDE_PIPE_SAT_MIN, PHASE_EXCLUDE_PIPE_SAT_MAX, SIGMA_SWITCH_ON_OFF_PIPE
        INTEGER :: PHASE_EXCLUDE
        LOGICAL, DIMENSION( : ), ALLOCATABLE :: PIPE_INDEX_LOGICAL
        REAL :: DIRECTION( Mdims%ndim ), DIRECTION_NORM( Mdims%ndim )
        REAL :: DX, ELE_ANGLE, NN, NM, suf_area, PIPE_DIAM_END, MIN_DIAM, U_GI, E_ROUGHNESS
        REAL :: S_WATER, S_WATER_MIN, S_WATER_MAX, SIGMA_SWITCH_ON_OFF_PIPE_GI, PIPE_SWITCH
        INTEGER :: pipe_corner_nds1( Mdims%ndim ), pipe_corner_nds2( Mdims%ndim ), NPIPES, ncorner, scvngi, &
            &     i_indx, j_indx, jdim, jphase, u_ljloc, u_jloc, ICORNER1, ICORNER2, ICORNER3, ICORNER4
        INTEGER :: SELE, CV_SILOC, JCV_NOD1, JCV_NOD2, IPRES, JU_NOD, CV_NOD, CV_LOC1, CV_LOC2, IPHASE_IN_PIPE, GI, IWATER

        ncorner = Mdims%ndim + 1
        PIPE_MIN_DIAM=.TRUE. ! Take the min diamter of the pipe as the real diameter.
        IGNORE_DIAGONAL_PIPES=.TRUE.
        SOLVE_ACTUAL_VEL = .TRUE. ! Solve for the actual real velocity in the pipes.
        CALC_SIGMA_PIPE = have_option("/porous_media/well_options/calculate_sigma_pipe")
        DEFAULT_SIGMA_PIPE_OPTIONS = .FALSE. ! Use default pipe options for water and oil including density and viscocity
        call get_option("/porous_media/well_options/calculate_sigma_pipe/pipe_roughness", E_ROUGHNESS, default=1.0E-6)
        ! Add the sigma associated with the switch to switch the pipe flow on and off...
        SWITCH_PIPES_ON_AND_OFF= have_option("/porous_media/well_options/switch_wells_on_and_off")
        if ( CALC_SIGMA_PIPE ) then
            allocate( well_density(Mdims%nphase), well_viscosity(Mdims%nphase) )
            do iphase = Mdims%n_in_pres+1, Mdims%nphase
                wd => extract_scalar_field( state(iphase), "Density" )
                wm => extract_tensor_field( state(iphase), "Viscosity" )
                well_density( iphase ) = wd%val(1)
                well_viscosity( iphase ) = wm%val(1,1,1)
            end do
        end if
        if ( SWITCH_PIPES_ON_AND_OFF ) then
            ! Define PHASE_EXCLUDE, PHASE_EXCLUDE_PIPE_SAT_MIN, PHASE_EXCLUDE_PIPE_SAT_MAX, SIGMA_SWITCH_ON_OFF_PIPE
            call get_option( "/porous_media/well_options/switch_wells_on_and_off/phase_exclude", phase_exclude )
            phase_exclude_pipe_sat_min => extract_scalar_field( state(1), "phase_exclude_pipe_sat_min" )
            phase_exclude_pipe_sat_max => extract_scalar_field( state(1), "phase_exclude_pipe_sat_max" )
            sigma_switch_on_off_pipe => extract_scalar_field( state(1), "sigma_switch_on_off_pipe" )
            cv_vol_frac => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
        end if
        ! Set rhs of the force balce equation to zero just for the pipes...
        U_RHS( :, Mdims%n_in_pres+1:Mdims%nphase, : ) = 0.0
        IF ( GET_PIVIT_MAT ) THEN
            DO U_ILOC = 1, Mdims%u_nloc
                DO U_JLOC = 1, Mdims%u_nloc
                    DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                        JPHASE = IPHASE
                        DO IDIM = 1, Mdims%ndim
                            JDIM = IDIM
                            i_indx = IDIM + (IPHASE-1)*Mdims%ndim + (U_ILOC-1)*Mdims%ndim*Mdims%nphase
                            j_indx = JDIM + (JPHASE-1)*Mdims%ndim + (U_JLOC-1)*Mdims%ndim*Mdims%nphase
                            pivit_mat(i_indx, j_indx, :) = 0.0
                        END DO
                    END DO
                END DO
            END DO
        END IF
        CV_QUADRATIC = (Mdims%cv_nloc==6.AND.Mdims%ndim==2) .OR. (Mdims%cv_nloc==10.AND.Mdims%ndim==3)
        U_QUADRATIC = (Mdims%u_nloc==6.AND.Mdims%ndim==2) .OR. (Mdims%u_nloc==10.AND.Mdims%ndim==3)
        ! The following is for 1d integration...
        scvngi = 2
        if ( CV_QUADRATIC ) then
            cv_lnloc = 3
            scvngi = 3
        else
            cv_lnloc = 2
        end if
        if ( U_QUADRATIC ) then
            u_lnloc = 3
            scvngi = 3
        else
            u_lnloc = 2
        end if
        PIPE_Diameter => EXTRACT_SCALAR_FIELD(STATE(1), "DiameterPipe1")
        X => EXTRACT_VECTOR_FIELD( PACKED_STATE, "PressureCoordinate" )
        allocate( PIPE_DIAM_GI(scvngi) )
        allocate( SIGMA_GI(Mdims%nphase,scvngi), SIGMA_ON_OFF_GI(Mdims%nphase,scvngi) )
        ! Get the 1D shape functions...
        allocate( scvfeweigh(scvngi), &
            scvfen(cv_lnloc, scvngi), scvfenslx(cv_lnloc, scvngi), scvfensly(cv_lnloc, scvngi), &
            scvfenlx(cv_lnloc, scvngi), scvfenly(cv_lnloc, scvngi), scvfenlz(cv_lnloc, scvngi), &
            sufen(u_lnloc, scvngi), sufenslx(u_lnloc, scvngi), sufensly(u_lnloc, scvngi), &
            sufenlx(u_lnloc, scvngi), sufenly(u_lnloc, scvngi), sufenlz(u_lnloc, scvngi) )
        allocate( detwei(scvngi), &
            l_cvfenx_all(cv_lnloc, scvngi), l_ufenx_all(u_lnloc, scvngi), &
            l_cvfenx_all_reversed(Mdims%ndim, scvngi, cv_lnloc) , &
            l_ufen_reversed(scvngi, u_lnloc), &
            nmx_all( Mdims%ndim ), X_ALL_CORN(Mdims%ndim, ncorner) )
        allocate(PIPE_INDEX_LOGICAL(ncorner))
        allocate( cv_loc_corner(ncorner), cv_mid_side(ncorner, ncorner), &
            u_loc_corner(ncorner), u_mid_side(ncorner, ncorner) )
        call fv_1d_quad( scvngi, cv_lnloc, scvfen, scvfenslx, scvfensly, scvfeweigh, &
            scvfenlx, scvfenly, scvfenlz ) ! For scalar fields
        call fv_1d_quad( scvngi, u_lnloc, sufen, sufenslx, sufensly, scvfeweigh, &
            sufenlx, sufenly, sufenlz ) ! For U fields
        ! Calculate the local corner nodes...
        CALL CALC_CORNER_NODS(CV_LOC_CORNER, Mdims%ndim, Mdims%cv_nloc, CV_QUADRATIC, CV_MID_SIDE)
        CALL CALC_CORNER_NODS(U_LOC_CORNER, Mdims%ndim, Mdims%u_nloc, U_QUADRATIC, U_MID_SIDE)
        IF ( U_QUADRATIC ) THEN
            U_LNLOC = 3
        ELSE
            U_LNLOC = 2
        END IF
        IF ( CV_QUADRATIC ) THEN
            CV_LNLOC = 3
        ELSE
            CV_LNLOC = 2
        END IF
        allocate( CV_GL_LOC(cv_lnloc), CV_GL_GL(cv_lnloc), X_GL_GL(cv_lnloc), MAT_GL_GL(cv_lnloc) )
        allocate( U_GL_LOC(U_lnloc), U_GL_GL(U_lnloc) )
        ! set up the surface B.C.'s
        allocate( WIC_P_BC_ALL_NODS(Mdims%npres, Mdims%cv_nonods) )
        allocate( SUF_P_BC_ALL_NODS(Mdims%npres, Mdims%cv_nonods), RVEC_SUM(Mdims%npres, Mdims%cv_nonods) )
        allocate( LOC_U_RHS_U_ILOC(Mdims%ndim, Mdims%nphase) )
        WIC_P_BC_ALL_NODS=0 ; RVEC_SUM=0.0 ; SUF_P_BC_ALL_NODS=0.0
        DO SELE = 1, Mdims%stotel
            DO IPRES = 2, Mdims%npres
                IF ( WIC_P_BC_ALL( 1, IPRES, SELE ) == WIC_P_BC_DIRICHLET ) THEN
                    DO CV_SILOC = 1, Mdims%cv_snloc
                        CV_NOD = P_SNDGLN((SELE-1)*Mdims%cv_snloc + CV_SILOC )
                        WIC_P_BC_ALL_NODS( IPRES, CV_NOD ) = WIC_P_BC_DIRICHLET
                        SUF_P_BC_ALL_NODS( IPRES, CV_NOD ) = SUF_P_BC_ALL_NODS( IPRES, CV_NOD ) + SUF_P_BC_ALL( 1,IPRES, (SELE-1)*Mdims%cv_snloc + CV_SILOC )
                        RVEC_SUM( IPRES, CV_NOD ) = RVEC_SUM( IPRES, CV_NOD ) + 1.0
                    END DO
                END IF
            END DO
        END DO
        DO CV_NOD = 1, Mdims%cv_nonods
            DO IPRES = 2, Mdims%npres
                IF ( WIC_P_BC_ALL_NODS(IPRES, CV_NOD) /= 0 ) SUF_P_BC_ALL_NODS(IPRES, CV_NOD) = SUF_P_BC_ALL_NODS(IPRES, CV_NOD) / RVEC_SUM(IPRES, CV_NOD)
            END DO
        END DO
        C( :, Mdims%n_in_pres+1:Mdims%nphase, : ) = 0.0
        DO ELE = 1, Mdims%totele
            ! Look for pipe indicator in element:
            PIPE_NOD_COUNT = 0
            PIPE_INDEX_LOGICAL = .false.
            DO ICORNER = 1, NCORNER ! Number of corner nodes in element
                CV_ILOC = CV_LOC_CORNER( ICORNER )
                CV_NODI = CV_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
                IF ( PIPE_Diameter%val( CV_NODI ) > 0.0 ) THEN
                    PIPE_NOD_COUNT = PIPE_NOD_COUNT + 1
                    PIPE_INDEX_LOGICAL(ICORNER) = .true.
                END IF
            END DO
            ELE_HAS_PIPE = PIPE_NOD_COUNT > 1
            npipes = PIPE_NOD_COUNT - 1
            IF ( ELE_HAS_PIPE ) THEN
                X_ALL_CORN(:, :) = x%val(:, X_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_LOC_CORNER(1:NCORNER)) )
                ! If we have more than one pipe then choose the 2 edges with the shortest sides
                ! and have a maximum of 2 pipes per element...
                ! Calculate the pipes within an element...
                ! we return the pipe corner nodes for each pipe in
                CALL CALC_PIPES_IN_ELE( X_ALL_CORN, PIPE_INDEX_LOGICAL, Mdims%ndim, &
                    pipe_corner_nds1, pipe_corner_nds2, npipes )
                DO IPIPE = 1, NPIPES
                    ! DEFINE CV_LILOC:
                    CV_LILOC = 1
                    ICORNER = pipe_corner_nds1( IPIPE )
                    ICORNER1 = ICORNER
                    CV_ILOC = CV_LOC_CORNER( ICORNER )
                    CV_GL_LOC( CV_LILOC ) = CV_ILOC
                    CV_LILOC = CV_LNLOC
                    ICORNER = pipe_corner_nds2( IPIPE )
                    ICORNER2 = ICORNER
                    CV_ILOC = CV_LOC_CORNER( ICORNER )
                    CV_GL_LOC( CV_LILOC ) = CV_ILOC
                    IF ( CV_QUADRATIC ) CV_GL_LOC(2) = CV_MID_SIDE( ICORNER1, ICORNER2 )
                    CV_GL_GL( : ) = CV_NDGLN( (ELE-1)*Mdims%cv_nloc + CV_GL_LOC(:) )
                    MAT_GL_GL( : ) = MAT_NDGLN( (ELE-1)*Mdims%cv_nloc + CV_GL_LOC(:) )
                    ! DEFINE U_LILOC:
                    U_LILOC = 1
                    ICORNER = pipe_corner_nds1(IPIPE)
                    U_ILOC = U_LOC_CORNER( ICORNER )
                    U_GL_LOC( U_LILOC ) = U_ILOC
                    U_LILOC = U_LNLOC
                    ICORNER = pipe_corner_nds2(IPIPE)
                    U_ILOC = U_LOC_CORNER( ICORNER )
                    U_GL_LOC( U_LILOC ) = U_ILOC
                    IF ( U_QUADRATIC ) U_GL_LOC( 2 ) = U_MID_SIDE( ICORNER1, ICORNER2 )
                    U_GL_GL(:) = U_NDGLN( (ELE-1)*Mdims%u_nloc + U_GL_LOC( : ) )
                    ! Calculate element angle sweeped out by element and pipe
                    IF ( Mdims%ndim == 2 ) THEN
                        ELE_ANGLE = PI
                    ELSE
                        ! find the nodes other than the pipe end corner nodes...
                        ICORNER3 = 0
                        DO ICORNER = 1, NCORNER
                            IF ( ICORNER /= ICORNER1 ) THEN
                                IF ( ICORNER /= ICORNER2 ) THEN
                                    IF ( ICORNER3 == 0 ) THEN
                                        ICORNER3 = ICORNER
                                    ELSE
                                        ICORNER4 = ICORNER
                                    END IF
                                END IF
                            END IF
                        END DO
                        ELE_ANGLE = CALC_ELE_ANGLE_3D( X_ALL_CORN(:, ICORNER1), X_ALL_CORN(:, ICORNER2) , &
                            &                         X_ALL_CORN(:, ICORNER3), X_ALL_CORN(:, ICORNER4) )
                    END IF
                    DIRECTION(:) = X_ALL_CORN( :, ICORNER2 ) - X_ALL_CORN( :, ICORNER1 )
                    DX = SQRT( SUM( DIRECTION(:)**2 ) )
                    DIRECTION(:) = DIRECTION(:) / DX
                    IF ( IGNORE_DIAGONAL_PIPES ) THEN
                        IF ( ABS(DIRECTION(1))<0.99 .AND. ABS(DIRECTION(2))<0.99.AND. ABS(DIRECTION(Mdims%ndim))<0.99 ) CYCLE
                    END IF
                    IF ( PIPE_MIN_DIAM ) THEN
                        MIN_DIAM = MINVAL( PIPE_diameter%val( CV_GL_GL( : ) ) )
                        PIPE_DIAM_GI(:) = MIN_DIAM
                        IF ( SOLVE_ACTUAL_VEL ) THEN
                            DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                                SIGMA_GI(IPHASE,:) = MINVAL( SIGMA(IPHASE, MAT_GL_GL( : ) ) )
                            END DO
                        END IF
                    ELSE
                        PIPE_DIAM_GI(:) = 0.0
                        DO CV_LILOC = 1, CV_LNLOC
                            CV_KNOD = CV_GL_GL( CV_LILOC )
                            PIPE_DIAM_GI(:) = PIPE_DIAM_GI(:) + PIPE_DIAMETER%val( CV_KNOD ) * SCVFEN( CV_LILOC, : )
                        END DO
                        PIPE_DIAM_GI(:) = MAX(PIPE_DIAM_GI(:), 0.0)
                        IF ( SOLVE_ACTUAL_VEL ) THEN
                            SIGMA_GI(:,:)=0.0
                            DO IPHASE=Mdims%n_in_pres+1,Mdims%nphase
                                DO CV_LILOC = 1, CV_LNLOC
                                    MAT_KNOD = MAT_GL_GL( CV_LILOC )
                                    SIGMA_GI(IPHASE,:) = SIGMA_GI(IPHASE,:) + SIGMA( IPHASE, MAT_KNOD ) * SCVFEN( CV_LILOC, : )
                                END DO
                            END DO
                        END IF
                    END IF
                    ! Recalculate SIGMA if we need to...
                    IF ( CALC_SIGMA_PIPE ) THEN
                        MIN_DIAM = MINVAL( PIPE_diameter%val( CV_GL_GL( : ) ) )
                        DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                            IPHASE_IN_PIPE=IPHASE-Mdims%n_in_pres
                            DO GI = 1, scvngi
                                U_GI = 0.0
                                DO IDIM = 1, Mdims%ndim
                                    U_GI = U_GI + SUM( SUFEN( : , GI ) * NU_ALL( IDIM, IPHASE, U_GL_GL( : ) ) ) * DIRECTION( IDIM )
                                END DO
                                CALL SIGMA_PIPE_FRICTION( SIGMA_GI( IPHASE, GI ), U_GI, MIN_DIAM, WELL_DENSITY( IPHASE ), WELL_VISCOSITY( IPHASE ), E_ROUGHNESS )
                            END DO
                        END DO
                    END IF
                    ! Add the sigma associated with the switch to switch the pipe flow on and off...
                    IF ( SWITCH_PIPES_ON_AND_OFF ) THEN
                        IWATER = PHASE_EXCLUDE
                        S_WATER = MINVAL( CV_VOL_FRAC%VAL( 1, IWATER, CV_GL_GL( : ) ) )
                        S_WATER_MIN = MAXVAL( PHASE_EXCLUDE_PIPE_SAT_MIN%VAL( CV_GL_GL( : ) ) )
                        S_WATER_MAX = MAXVAL( PHASE_EXCLUDE_PIPE_SAT_MAX%VAL( CV_GL_GL( : ) ) )
                        SIGMA_SWITCH_ON_OFF_PIPE_GI = MAXVAL( SIGMA_SWITCH_ON_OFF_PIPE%VAL( CV_GL_GL( : ) ) )
                        PIPE_SWITCH = 1.0 - MIN( 1.0, MAX( 0.0, ( S_WATER_MAX - S_WATER ) / MAX( S_WATER_MAX - S_WATER_MIN, 1.E-20 ) ) )
                        SIGMA_ON_OFF_GI( Mdims%n_in_pres+1:Mdims%nphase, : ) = PIPE_SWITCH * SIGMA_SWITCH_ON_OFF_PIPE_GI
                    END IF
                    ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
                    ! Adjust according to the volume of the pipe...
                    DETWEI(:) = SCVFEWEIGH(:) * 0.5* DX * PI * ( (0.5*PIPE_DIAM_GI(:))**2 ) * ELE_ANGLE / ( 2.0 * PI )
                    L_CVFENX_ALL(:,:) = 2.0 * SCVFENLX(:,:) / DX
                    L_UFENX_ALL(:,:) = 2.0 * SUFENLX(:,:) / DX
                    DO CV_LILOC = 1, CV_LNLOC
                        DO IDIM = 1, Mdims%ndim
                            L_CVFENX_ALL_REVERSED( IDIM, :, CV_LILOC ) = L_CVFENX_ALL( CV_LILOC, : ) * DIRECTION( IDIM )
                        END DO
                    END DO
                    DO U_LILOC = 1, U_LNLOC
                        L_UFEN_REVERSED( :, U_LILOC ) = SUFEN( U_LILOC, : )
                    END DO
                    ! Add in C matrix contribution: (DG velocities)
                    DO U_LILOC = 1, U_LNLOC
                        IU_NOD = U_GL_GL( U_LILOC ) ! U_NDGLN( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                        DO P_LJLOC = 1, CV_LNLOC
                            JCV_NOD = CV_GL_GL(P_LJLOC) ! P_NDGLN( ( ELE - 1 ) * Mdims%p_nloc + P_JLOC )
                            ! In this section we multiply the shape functions over the GI points. i.e: we perform the integration
                            ! over the element of the pressure like source term.
                            ! Put into matrix
                            DO COUNT2 = Mspars%C%fin(IU_NOD), Mspars%C%fin(IU_NOD+1)-1
                                IF ( Mspars%C%col(COUNT2) == JCV_NOD ) COUNT = COUNT2
                            END DO
                            ! Prepare aid variable NMX_ALL to improve the speed of the calculations
                            NMX_ALL( : ) = matmul( L_CVFENX_ALL_REVERSED( :, :, P_LJLOC ), DETWEI( : ) * L_UFEN_REVERSED( :, U_LILOC ) )
                            NM = sum( L_UFEN_REVERSED( :, U_LILOC ) * SCVFEN( P_LJLOC, : ) * DETWEI( : ) )
                            ! Put into matrix
                            DO IDIM = 1, Mdims%ndim
                                C( IDIM, Mdims%n_in_pres+1:Mdims%nphase, COUNT ) = C( IDIM, Mdims%n_in_pres+1:Mdims%nphase, COUNT ) - NMX_ALL( IDIM )
                            END DO
                            DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                                DO IDIM = 1, Mdims%ndim
                                    ! This is ( S \dot n ) n
                                    U_RHS( IDIM, IPHASE, IU_NOD ) =  U_RHS( IDIM, IPHASE, IU_NOD ) + &
                                        NM * SUM(U_SOURCE_CV( :, IPHASE, JCV_NOD ) * DIRECTION( : ) ) * DIRECTION( IDIM )
                                END DO
                            END DO
                        END DO ! DO P_LJLOC = 1, CV_LNLOC
                        U_ILOC = U_GL_LOC( U_LILOC )
                        DO U_LJLOC = 1, U_LNLOC
                            U_JLOC = U_GL_LOC( U_LJLOC )
                            JU_NOD = U_GL_GL( U_LJLOC )
                            IF ( GET_PIVIT_MAT ) THEN
                                IF ( .NOT.SOLVE_ACTUAL_VEL ) THEN
                                    NN = SUM( L_UFEN_REVERSED( :, U_LILOC ) * L_UFEN_REVERSED( :, U_LJLOC ) * DETWEI( : ) )
                                END IF
                                DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                                    IF ( SOLVE_ACTUAL_VEL ) THEN
                                        NN = SUM( L_UFEN_REVERSED( :, U_LILOC ) * L_UFEN_REVERSED( :, U_LJLOC ) * DETWEI( : ) * SIGMA_GI( IPHASE , : ) )
                                    END IF
                                    JPHASE = IPHASE
                                    DO IDIM = 1, Mdims%ndim
                                        DO JDIM = 1, Mdims%ndim
                                            i_indx = IDIM + (IPHASE-1)*Mdims%ndim + (U_ILOC-1)*Mdims%ndim*Mdims%nphase
                                            j_indx = JDIM + (JPHASE-1)*Mdims%ndim + (U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                            pivit_mat( i_indx, j_indx, ele ) = pivit_mat( i_indx, j_indx, ele ) + NN * DIRECTION( IDIM ) * DIRECTION( JDIM )
                                           !pivit_mat( i_indx, i_indx, ele ) = pivit_mat( i_indx, i_indx, ele ) + NN * DIRECTION( IDIM ) * DIRECTION( JDIM )
                                        END DO
                                    END DO
                                END DO
                                IF ( SWITCH_PIPES_ON_AND_OFF ) THEN
                                    DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                                        IF ( SOLVE_ACTUAL_VEL ) THEN
                                            NN = SUM( L_UFEN_REVERSED( :, U_LILOC ) * L_UFEN_REVERSED( :, U_LJLOC ) * DETWEI( : ) * SIGMA_ON_OFF_GI(IPHASE,:) )
                                        END IF
                                        JPHASE = IPHASE
                                        DO IDIM = 1, Mdims%ndim
                                            DO JDIM = 1, Mdims%ndim
                                                i_indx = IDIM + (IPHASE-1)*Mdims%ndim + (U_ILOC-1)*Mdims%ndim*Mdims%nphase
                                                j_indx = JDIM + (JPHASE-1)*Mdims%ndim + (U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                                pivit_mat( i_indx, j_indx, ele ) = pivit_mat( i_indx, j_indx, ele ) + NN * DIRECTION( IDIM ) * DIRECTION( JDIM )
                                               !pivit_mat( i_indx, i_indx, ele ) = pivit_mat( i_indx, i_indx, ele ) + NN * DIRECTION( IDIM ) * DIRECTION( JDIM )
                                            END DO
                                        END DO
                                    END DO
                                END IF ! SWITCH_PIPES_ON_AND_OFF
                            END IF ! GET_PIVIT_MAT
                            NN = sum( L_UFEN_REVERSED( :, U_LILOC ) * L_UFEN_REVERSED( :, U_LJLOC ) * DETWEI( : ) )
                            DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                                DO IDIM = 1, Mdims%ndim
                                    ! This is ( S \dot n ) n
                                    U_RHS( IDIM, IPHASE, IU_NOD ) =  U_RHS( IDIM, IPHASE, IU_NOD ) + &
                                        NN * SUM(U_SOURCE( :, IPHASE, JU_NOD ) * DIRECTION( : ) ) * DIRECTION( IDIM )
                                END DO
                            END DO
                        END DO ! DO U_LJLOC = 1, U_LNLOC
                    END DO ! DO U_LILOC = 1, U_LNLOC
                    ! Add pressure b.c to matrix and rhs...
                    CV_LOC1 = CV_LOC_CORNER( ICORNER1 )
                    JCV_NOD1 = CV_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_LOC1 )
                    CV_LOC2 = CV_LOC_CORNER( ICORNER2 )
                    JCV_NOD2 = CV_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_LOC2 )
                    DO IPRES = 2, Mdims%npres
                        JCV_NOD = 0
                        IF ( WIC_P_BC_ALL_NODS( IPRES, JCV_NOD1 ) == WIC_P_BC_DIRICHLET ) THEN
                            CV_LILOC = 1
                            JCV_NOD = JCV_NOD1
                            U_LILOC = 1
                            JU_NOD = U_GL_GL( U_LILOC )
                            direction_norm = - direction ! for the b.c it must be -ve at the bottom of element
                        END IF
                        IF ( WIC_P_BC_ALL_NODS( IPRES, JCV_NOD2 ) == WIC_P_BC_DIRICHLET ) THEN
                            CV_LILOC = CV_LNLOC
                            JCV_NOD = JCV_NOD2
                            U_LILOC = U_LNLOC
                            JU_NOD = U_GL_GL( U_LILOC )
                            direction_norm = +direction ! for the b.c it must be +ve at the top of element
                        END IF
                        IF ( JCV_NOD /= 0 ) THEN
                            ! Add in C matrix contribution: (DG velocities)
                            ! In this section we multiply the shape functions over the GI points. i.e: we perform the integration
                            ! over the element of the pressure like source term.
                            ! Put into matrix
                            PIPE_DIAM_END = PIPE_diameter%val( JCV_NOD )
                            DO COUNT2 = Mspars%C%fin(JU_NOD), Mspars%C%fin(JU_NOD+1)-1
                                IF ( Mspars%C%col(COUNT2) == JCV_NOD ) COUNT = COUNT2
                            END DO
                            ! Prepare aid variable NMX_ALL to improve the speed of the calculations
                            suf_area = PI * ( (0.5*PIPE_DIAM_END)**2 ) * ELE_ANGLE/ ( 2.0 * PI )
                            NMX_ALL( : ) = direction_norm(:)* suf_area
                            LOC_U_RHS_U_ILOC = 0.0
                            DO IPHASE = 1+(IPRES-1)*Mdims%n_in_pres, IPRES*Mdims%n_in_pres
                                DO IDIM = 1, Mdims%ndim
                                    C( IDIM, IPHASE, COUNT ) = C( IDIM, IPHASE, COUNT ) + NMX_ALL( IDIM )
                                    LOC_U_RHS_U_ILOC( IDIM, IPHASE) =  LOC_U_RHS_U_ILOC( IDIM, IPHASE ) &
                                        - NMX_ALL( IDIM ) * SUF_P_BC_ALL_NODS( IPRES,JCV_NOD )
                                END DO
                            END DO
                            U_RHS( :, Mdims%n_in_pres+1:Mdims%nphase, JU_NOD ) = U_RHS( :, Mdims%n_in_pres+1:Mdims%nphase, JU_NOD ) + LOC_U_RHS_U_ILOC( :, Mdims%n_in_pres+1:Mdims%nphase)
                        END IF ! IF ( JCV_NOD /= 0 ) THEN
                    END DO ! DO IPRES = 2, Mdims%npres
                END DO ! DO IPIPE = 1, NPIPES
            END IF ! IF ( ELE_HAS_PIPE ) THEN
        END DO ! DO ELE = 1, Mdims%totele
        if ( GET_PIVIT_MAT ) then
            DO U_ILOC = 1, Mdims%u_nloc
                U_JLOC = U_ILOC
                DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                    JPHASE = IPHASE
                    DO IDIM = 1, Mdims%ndim
                        JDIM = IDIM
                        i_indx = IDIM + (IPHASE-1)*Mdims%ndim + (U_ILOC-1)*Mdims%ndim*Mdims%nphase
                        j_indx = JDIM + (JPHASE-1)*Mdims%ndim + (U_JLOC-1)*Mdims%ndim*Mdims%nphase
                        WHERE ( abs( pivit_mat(i_indx, j_indx, :) ) < 1.0e-10 )
                            pivit_mat(i_indx, j_indx, :) = 1.0
                        END WHERE
                    END DO
                END DO
            END DO
        end if
        RETURN
    END SUBROUTINE MOD_1D_FORCE_BAL_C


    SUBROUTINE SIGMA_PIPE_FRICTION( SIGMA, U, DIAM, DEN, VISC, E_ROUGHNESS )
        IMPLICIT NONE
        REAL, INTENT( IN ) :: U,DIAM,DEN,VISC,E_ROUGHNESS
        REAL, INTENT( INOUT ) :: SIGMA

        REAL, PARAMETER :: TOLER=1.E-10
        REAL :: RE, A, MAX_RE, F, F1, F2, W, MAX_U

        MAX_U = MAX( TOLER, ABS( U ) )
        RE = DEN * MAX_U * DIAM / VISC

        MAX_RE = RE
        A = LOG10( 6.9/MAX_RE + (E_ROUGHNESS/(3.7*DIAM))**(10.0/9.0) )

        F1 = 16.0 / MAX_RE
        F2 = 1.0 / MAX( -3.6 * A, TOLER )

        W = MAX( 0.0,  MIN( 1.0, (RE-2000.0)/2000.0 ) ) ! Relaxation factor between 2 friction factor expressions

        F = ( 1 - W ) * F1 + W * F2

        SIGMA = ( F / DIAM ) * 2.0 * DEN * ABS( U ) ! Based on Fanning friction factor
        !SIGMA = ( F / DIAM ) * 2.0 * DEN * 1.0 ! Based on Fanning friction factor
        !SIGMA = ( F / DIAM ) * 2.0 * DEN * MAX_U ! Based on Fanning friction factor

        RETURN
    END SUBROUTINE SIGMA_PIPE_FRICTION





    SUBROUTINE CALC_CORNER_NODS( CV_LOC_CORNER, NDIM, CV_NLOC, CV_QUADRATIC, CV_MID_SIDE )
        ! Calculate the local corner nodes...
        ! CV_MID_SIDE(ICORN,JCORN)= CV_ILOC local node number for node between these two corner nodes
        IMPLICIT NONE
        INTEGER, INTENT( IN ) :: CV_NLOC, NDIM
        INTEGER, DIMENSION( : ), INTENT( INOUT ) :: CV_LOC_CORNER
        INTEGER, DIMENSION( :, : ), INTENT( INOUT ) :: CV_MID_SIDE
        LOGICAL, INTENT( INOUT ) :: CV_QUADRATIC
        INTEGER :: CV_ILOC, ICORN, JCORN, CV_NCORNER

        CV_NCORNER = NDIM + 1

        IF( ((CV_NLOC==3).AND.(NDIM==2)) .OR. ((CV_NLOC==4).AND.(NDIM==3)) ) THEN ! assume linear...
            CV_QUADRATIC=.FALSE.
            DO CV_ILOC=1,CV_NLOC
                CV_LOC_CORNER(CV_ILOC)=CV_ILOC
            END DO
        ELSE ! assume quadratic...
            CV_QUADRATIC=.TRUE.
            CV_LOC_CORNER(1)=1
            CV_LOC_CORNER(2)=3
            CV_LOC_CORNER(3)=6

            CV_MID_SIDE=0
            CV_MID_SIDE(1,2)=2
            CV_MID_SIDE(1,3)=4
            CV_MID_SIDE(2,3)=5

            IF(NDIM==3) THEN
                CV_MID_SIDE(1,4)=7
                CV_MID_SIDE(2,4)=8
                CV_MID_SIDE(3,4)=9
                CV_LOC_CORNER(4)=10
            END IF
            DO ICORN=1,CV_NCORNER
                DO JCORN=ICORN+1,CV_NCORNER
                    CV_MID_SIDE(JCORN,ICORN)=CV_MID_SIDE(ICORN,JCORN)
                END DO
            END DO
        END IF

        RETURN
    END SUBROUTINE CALC_CORNER_NODS

    SUBROUTINE CALC_PIPES_IN_ELE( X_ALL_CORN, PIPE_INDEX_LOGICAL, NDIM, &
        pipe_corner_nds1, pipe_corner_nds2, npipes )
        ! Calculate the pipes within an element...
        ! Return the pipe corner nodes for each pipe in element.
        IMPLICIT NONE
        INTEGER, intent( in ) :: NDIM, npipes
        REAL, intent( in ) :: X_ALL_CORN(NDIM,NDIM+1)
        LOGICAL, intent( in ) :: PIPE_INDEX_LOGICAL(NDIM+1)
        integer, intent( inout ) :: pipe_corner_nds1(npipes), pipe_corner_nds2(npipes)
        ! logal variables...
        REAL :: cp(ndim)
        REAL :: area_sqr(ndim+1)
        INTEGER :: I,J,K, II(ndim+1),JJ(ndim+1),KK(ndim+1), ipipe, i_nd,j_nd, iSTORE(4), jSTORE(4), kSTORE(4)
        INTEGER :: icorn, iface, iii, ik, il
        logical :: found_i, found_j
        real :: max_area

        IF ( NPIPES==1 ) THEN

            i_nd = 0
            DO ICORN = 1, NDIM+1
                IF ( PIPE_INDEX_LOGICAL(ICORN) ) THEN
                    IF ( I_ND == 0 ) THEN
                        I_ND = ICORN
                    ELSE
                        J_ND = ICORN
                    END IF
                END IF
            END DO
            ipipe = 1
            pipe_corner_nds1(ipipe) = i_nd
            pipe_corner_nds2(ipipe) = j_nd

        ELSE  ! IF ( NPIPES==1 ) THEN

            IF ( NDIM==2 ) THEN

                i=1 ; j=2
                area_sqr(1) = sum( ( X_ALL_CORN(:,1) - X_ALL_CORN(:,2) )**2 ) ! THIS IS THE LENGTH SQUARED
                ii(1)=i ; jj(1)=j

                i=1 ; j=3
                area_sqr(2) = sum( ( X_ALL_CORN(:,1) - X_ALL_CORN(:,3) )**2 )
                ii(2)=i ; jj(2)=j

                i=2 ; j=3
                area_sqr(3) = sum( ( X_ALL_CORN(:,2) - X_ALL_CORN(:,3) )**2 )
                ii(3)=i ; jj(3)=j

                kk(:)=0

            ELSE ! ENDOF IF(NDIM==2) THEN...(THIS IS FOR 3D)

                ! NB. The cross product is directly related to the area of the triangle...
                iSTORE(1)=1 ; jSTORE(1)=2 ; kSTORE(1)=3
                iSTORE(2)=1 ; jSTORE(2)=2 ; kSTORE(2)=4
                iSTORE(3)=2 ; jSTORE(3)=3 ; kSTORE(3)=4
                iSTORE(4)=1 ; jSTORE(4)=3 ; kSTORE(4)=4

                IFACE = 0
                DO III = 1, 4
                    i=iSTORE(III) ; j=jSTORE(III) ; k=kSTORE(III)
                    IF ( PIPE_INDEX_LOGICAL(I) .AND. PIPE_INDEX_LOGICAL(J) .AND. PIPE_INDEX_LOGICAL(K) ) THEN
                        IFACE = IFACE+1
                        call CrossProduct( ndim, cp,  ( X_ALL_CORN(:,j) - X_ALL_CORN(:,i) ), ( X_ALL_CORN(:,k) - X_ALL_CORN(:,i) ) )
                        area_sqr(IFACE) =  sum(cp(:)**2 )
                        ii(IFACE)=i ; jj(IFACE)=j ; kk(IFACE)=k
                    END IF
                END DO ! ENDOF DO III=1,4

                IF(IFACE==1) THEN ! A 3d adjustment for when we have only 2 pipes...
                    i=ii(IFACE) ; j=jj(IFACE) ; k=kk(IFACE) ! nodes on the face

                    area_sqr(1) = sum( ( X_ALL_CORN(:,i) - X_ALL_CORN(:,j) )**2 ) ! THIS IS THE LENGTH SQUARED

                    area_sqr(2) = sum( ( X_ALL_CORN(:,i) - X_ALL_CORN(:,k) )**2 ) ! THIS IS THE LENGTH SQUARED

                    area_sqr(3) = sum( ( X_ALL_CORN(:,j) - X_ALL_CORN(:,k) )**2 ) ! THIS IS THE LENGTH SQUARED

                    ik = 1
                    max_area = area_sqr(ik)
                    do il = 2, 3
                        if ( max_area < area_sqr(il) ) then
                            max_area = area_sqr(il)
                            ik = il
                        end if
                    end do

                    if(ik==1) then
                        ipipe = 1
                        pipe_corner_nds1(ipipe) = i
                        pipe_corner_nds2(ipipe) = k
                        ipipe = 2
                        pipe_corner_nds1(ipipe) = j
                        pipe_corner_nds2(ipipe) = k
                    else if(ik==2) then
                        ipipe = 1
                        pipe_corner_nds1(ipipe) = i
                        pipe_corner_nds2(ipipe) = j
                        ipipe = 2
                        pipe_corner_nds1(ipipe) = j
                        pipe_corner_nds2(ipipe) = k
                    else
                        ipipe = 1
                        pipe_corner_nds1(ipipe) = i
                        pipe_corner_nds2(ipipe) = j
                        ipipe = 2
                        pipe_corner_nds1(ipipe) = i
                        pipe_corner_nds2(ipipe) = k
                    endif
                    RETURN

                ENDIF  ! IF(IFACE==1) THEN

            END IF ! ENDOF IF(NDIM==2) THEN

            ! Find the maximum surface area..
            ik = 1
            max_area = area_sqr(ik)
            do i = 2, NDIM+1
                if ( max_area < area_sqr(i) ) then
                    max_area = area_sqr(i)
                    ik = i
                end if
            end do

            ! Now work out the pipes that do not contain the edges on this maximum surface
            ! and contain
            ipipe = 0
            do i_nd = 1, ndim+1
                do j_nd = i_nd+1, ndim+1
                    ! If edge nodes i_nd with j_nd are both in the ii,jj or kk
                    IF ( PIPE_INDEX_LOGICAL(i_nd) .AND. PIPE_INDEX_LOGICAL(J_ND) ) THEN
                        ! Make sure both nodes do not belong to the face with the maximum area...
                        found_i = .false.
                        found_j = .false.

                        if ( i_nd==ii(ik) ) found_i = .true.
                        if ( i_nd==jj(ik) ) found_i = .true.
                        if ( i_nd==kk(ik) ) found_i = .true.

                        if ( j_nd==ii(ik) ) found_j = .true.
                        if ( j_nd==jj(ik) ) found_j = .true.
                        if ( j_nd==kk(ik) ) found_j = .true.

                        if ( (.not.found_i) .and. (.not.found_j) ) then
                            ipipe = ipipe + 1
                            pipe_corner_nds1(ipipe) = i_nd
                            pipe_corner_nds2(ipipe) = j_nd
                        end if

                    end if
                end do
            end do

        END IF ! IF ( NPIPES==1 ) THEN ELSE

        RETURN
    END SUBROUTINE CALC_PIPES_IN_ELE



    REAL FUNCTION CALC_ELE_ANGLE_3D( X_ALL_CORN_PIPE1, X_ALL_CORN_PIPE2, X_ALL_CORN_PIPE3, X_ALL_CORN_PIPE4 )
        ! Calculate element angle sweeped out by element and pipe
        ! X_ALL_CORN_PIPE1, X_ALL_CORN_PIPE2 are the coordinates of the ends of the pipe within an element.
        ! X_ALL_CORN_PIPE3, X_ALL_CORN_PIPE4 are the other corner 2 nodes of an element.
        IMPLICIT NONE
        REAL, intent( in ) :: X_ALL_CORN_PIPE1(3), X_ALL_CORN_PIPE2(3),  X_ALL_CORN_PIPE3(3), X_ALL_CORN_PIPE4(3)
        REAL :: X_PIPE1(3), X_PIPE2(3), X_PIPE3(3), X_PIPE4(3)
        REAL :: X_PIPE3_N(3), X_PIPE4_N(3)
        REAL :: X_PIPE3_2D(2), X_PIPE4_2D(2)
        REAL :: NORM(3), T1(3), T2(3), A(3,3)
        INTEGER :: IDIM

        ! Translate
        X_PIPE1 = 0.0
        X_PIPE2 = X_ALL_CORN_PIPE2 - X_ALL_CORN_PIPE1
        X_PIPE3 = X_ALL_CORN_PIPE3 - X_ALL_CORN_PIPE1
        X_PIPE4 = X_ALL_CORN_PIPE4 - X_ALL_CORN_PIPE1

        ! Rotate to be aligned with x-axis
        Norm = X_PIPE2 / SQRT( SUM( X_PIPE2(:)**2 ) )

        ! Determine the tangent and bi-normal vectors from the normal NormX, NormY, NormZ
        CALL GET_TANG_BINORM( Norm(1), Norm(2), Norm(3), T1(1), T1(2), T1(3), T2(1), T2(2), T2(3), 1 )

        A(1,:) = NORM(:)
        A(2,:) = T1(:)
        A(3,:) = T2(:)

        DO IDIM = 1, 3
            X_PIPE3_N(IDIM) = SUM( A(IDIM,:) * X_PIPE3(:) )
            X_PIPE4_N(IDIM) = SUM( A(IDIM,:) * X_PIPE4(:) )
        END DO

        X_PIPE3_2D(1:2) = X_PIPE3_N(2:3)
        X_PIPE4_2D(1:2) = X_PIPE4_N(2:3)

        ! Use cosine rule....
        CALC_ELE_ANGLE_3D = ACOS( SUM( X_PIPE3_2D(:) * X_PIPE4_2D(:) ) / ( SQRT( SUM( X_PIPE3_2D(:)**2 ) ) * SQRT( SUM( X_PIPE4_2D(:)**2 ) ) ) )

        RETURN
    END FUNCTION CALC_ELE_ANGLE_3D


end module multi_pipes
