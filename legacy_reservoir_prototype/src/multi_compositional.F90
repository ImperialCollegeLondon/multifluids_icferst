
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

module Compositional_Terms
    use fldebug
    use futils
    use shape_functions
    use shape_functions_prototype
    use shape_functions_linear_quadratic
    use matrix_operations
    use Copy_Outof_State
    use multi_data_types
    use futils, only: int2str
    implicit none

contains

    subroutine Calculate_ComponentAbsorptionTerm( packed_state, icomp, cv_ndgln, &
                                                  Mdims, denold, volfra_pore, mass_ele, comp_absorb )

        !!$ Calculate compositional model linkage between the phase expressed in COMP_ABSORB.
        !!$ Use values from the previous time step so its easier to converge.
        !!$ ALPHA_BETA is the scaling coeff. of the compositional model e.g. =1.0

        implicit none
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent( in ) :: Mdims
        integer, intent( in ) :: icomp
        integer, dimension( : ), intent( in ) :: cv_ndgln
        real, dimension( : ), intent( in ) :: mass_ele
        real, dimension( :, :, : ), intent( in ) :: denold
        real, dimension( :, : ), intent( in ) :: volfra_pore
        real, dimension( :, :, : ), intent( inout ) :: comp_absorb

        ! Local Variables
        integer :: &
            iphase, jphase, ele, cv_iloc, cv_nod, jcomp
        real :: dt, alpha_beta, max_k, min_k, alpha
        character( len = option_path_len ) :: option_path
        logical :: KComp_Sigmoid
        real, dimension( : ), allocatable :: sum_nod, volfra_pore_nod
        real, dimension( :, :, : ), allocatable :: k_comp
        real, dimension( :, :, :, : ), allocatable :: k_comp2
        !working pointers
        type(tensor_field), pointer :: tfield
        real, dimension(:,:), pointer :: satura
        !Initialize comp_absorb
        comp_absorb = 0.0

        if ( Mdims%nphase < 2 ) return

        tfield=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
        satura => tfield%val(1,:,:)

        allocate( sum_nod( Mdims%cv_nonods ), volfra_pore_nod( Mdims%cv_nonods ), &
            k_comp( Mdims%ncomp, Mdims%nphase, Mdims%nphase ), k_comp2( Mdims%ncomp, Mdims%cv_nonods, Mdims%nphase, Mdims%nphase ) )
        k_comp = 0.0 ; k_comp2 = 0.0

        option_path = 'material_phase[' // int2str( Mdims%nstate - Mdims%ncomp ) // &
            ']/is_multiphase_component'
        call get_option( trim( option_path ) // '/alpha_beta', alpha_beta, default = 1. )


        KComp_Sigmoid = have_option( trim( option_path ) // '/KComp_Sigmoid' )
        if( KComp_Sigmoid ) then
            do jcomp = 1, Mdims%ncomp
                call get_option( 'material_phase[' // int2str( Mdims%nphase + jcomp - 1 ) // &
                    ']/is_multiphase_component/KComp_Sigmoid/K_Comp', k_comp( jcomp, 1, 1 ) )
                k_comp( jcomp, :, : ) = k_comp( jcomp, 1, 1 )
            end do
        end if
        call get_option( '/timestepping/timestep', dt )

        !!$ Determine a node-wise representation of porosity VOLFRA_PORE_NOD.
        SUM_NOD = 0.0 ; VOLFRA_PORE_NOD = 0.0
        DO ELE = 1, Mdims%totele
            DO CV_ILOC = 1, Mdims%cv_nloc
                CV_NOD = CV_NDGLN( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
                SUM_NOD( CV_NOD ) = SUM_NOD( CV_NOD ) + mass_ele( ele )
                VOLFRA_PORE_NOD( CV_NOD ) = VOLFRA_PORE_NOD( CV_NOD ) + &
                    VOLFRA_PORE( 1, ELE ) * mass_ele( ele )
            END DO
        END DO
        VOLFRA_PORE_NOD = VOLFRA_PORE_NOD / SUM_NOD

        MIN_K = max( 1.e-1, MINVAL( K_COMP( ICOMP, : , : )))
        MAX_K = MAXVAL( K_COMP( ICOMP, : , : ) )
        CALL Calc_KComp2( Mdims%cv_nonods, Mdims%nphase, icomp, KComp_Sigmoid, &
            min( 1., max( 0., satura )), K_Comp, max_k, min_k, &
            K_Comp2 )


        DO CV_NOD = 1, Mdims%cv_nonods
            DO IPHASE = 1, Mdims%nphase
                DO JPHASE = IPHASE + 1, Mdims%nphase
                    ALPHA= ALPHA_BETA * VOLFRA_PORE_NOD( CV_NOD ) * &
                        ( max( 0.0, SATURA( IPHASE, CV_NOD ) * &
                        DENOLD( 1, IPHASE, CV_NOD ) ) / &
                        K_COMP2( ICOMP, CV_NOD, IPHASE, JPHASE ) + &
                        max( 0.0, SATURA( JPHASE,CV_NOD ) * &
                        DENOLD( 1, JPHASE, CV_NOD ) ) ) / DT
                    COMP_ABSORB( IPHASE, IPHASE, CV_NOD ) = &
                        COMP_ABSORB( IPHASE, IPHASE, CV_NOD ) + &
                        ALPHA* &
                        K_COMP2( ICOMP, CV_NOD, IPHASE, JPHASE )

                    COMP_ABSORB( IPHASE, JPHASE, CV_NOD ) = &
                        - ALPHA
                END DO
            END DO
        END DO

        DO CV_NOD = 1, Mdims%cv_nonods
            DO IPHASE = 1, Mdims%nphase
                DO JPHASE = 1, IPHASE - 1

                    ALPHA= ALPHA_BETA * VOLFRA_PORE_NOD( CV_NOD ) * &
                        ( max(0.0,SATURA (IPHASE, CV_NOD ) * &
                        DENOLD( 1, IPHASE, CV_NOD ) ) + &
                        max(0.0,SATURA (JPHASE, CV_NOD ) * &
                        DENOLD( 1, JPHASE, CV_NOD ) ) / &
                        K_COMP2( ICOMP, CV_NOD, JPHASE, IPHASE ) ) / DT

                    COMP_ABSORB( IPHASE, IPHASE, CV_NOD ) = &
                        COMP_ABSORB( IPHASE, IPHASE, CV_NOD ) + ALPHA

                    COMP_ABSORB( IPHASE, JPHASE, CV_NOD ) = - ALPHA* &
                        K_COMP2( ICOMP, CV_NOD, JPHASE, IPHASE )

                END DO
            END DO

        END DO

        do cv_nod = 1, Mdims%cv_nonods
            if( satura( 1, cv_nod ) > 0.95 ) then
              do iphase = 1, Mdims%nphase
                 do jphase = min( iphase + 1, Mdims%nphase ), Mdims%nphase
                    Comp_Absorb( iphase, jphase, cv_nod ) = &
                       Comp_Absorb( iphase, jphase, cv_nod ) * max( 0.01, &
                       20.0 * ( 1. - satura ( 1, cv_nod ) ) )
                 end do
              end do
            end if
        end do


        deallocate( sum_nod, volfra_pore_nod, k_comp, k_comp2 )

        RETURN

    end subroutine Calculate_ComponentAbsorptionTerm


    subroutine Calculate_ComponentDiffusionTerm( packed_state, &
        Mdims, CV_GIdims, CV_funs,&
        mat_ndgln, u_ndgln, x_ndgln, &
        ncomp_diff_coef, comp_diffusion_opt, &
        comp_diff_coef, &
        comp_diffusion)
        !!$ Calculate the diffusion coefficient COMP_DIFFUSION for current composition...
        !!$ based on page 136 in Reservoir-Simulation-Mathematical-Techniques-In-Oil-Recovery-(2007).pdf
        !!$ COMP_DIFFUSION_OPT, integer option defining diffusion coeff
        !!$ NCOMP_DIFF_COEF,  integer defining how many coeff's are needed to define the diffusion
        !!$ COMP_DIFF_COEF( Mdims%ncomp,  NCOMP_DIFF_COEF, Mdims%nphase  )
        implicit none
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_GI_dimensions), intent(in) :: CV_GIdims
        type(multi_shape_funs), intent(in) :: CV_funs
        integer, dimension( : ), intent( in ) :: mat_ndgln, u_ndgln, x_ndgln
        integer, intent( in ) :: ncomp_diff_coef, comp_diffusion_opt
        real, dimension( :, : ), intent( in ) :: comp_diff_coef
        real, dimension( :, :, :, : ),intent( inout ) :: comp_diffusion
        !!$ Local variables:
        integer :: ele, mat_nod, iphase, idim, u_iloc, u_inod
        real :: diff_molecular, diff_longitudinal, diff_transverse
        real, dimension( : ), allocatable :: ud, mat_u, nu, nv, nw
        type( vector_field ), pointer :: x_all
        type( tensor_field ), pointer :: nu_all

        if( comp_diffusion_opt == 0 ) then
            comp_diffusion = 0.0
            return
        endif
        x_all => extract_vector_field( packed_state, "PressureCoordinate" )

        allocate( NU(  Mdims%u_nonods ) ) ; NU = 0.0
        allocate( NV(  Mdims%u_nonods ) ) ; NV = 0.0
        allocate( NW(  Mdims%u_nonods ) ) ; NW = 0.0
        nu_all => extract_tensor_field( packed_state, "PackedNonlinearVelocity" )
        DO ELE = 1, Mdims%totele
            DO U_ILOC = 1, Mdims%u_nloc
                U_INOD = U_NDGLN( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                DO IPHASE = 1, Mdims%nphase
                    DO IDIM = 1, Mdims%ndim
                        IF ( IDIM==1 ) THEN
                            NU( U_INOD + (IPHASE-1)*Mdims%u_nonods ) = NU_ALL % VAL( IDIM, IPHASE, U_INOD )
                        ELSE IF ( IDIM==2 ) THEN
                            NV( U_INOD + (IPHASE-1)*Mdims%u_nonods ) = NU_ALL % VAL( IDIM, IPHASE, U_INOD )
                        ELSE
                            NW( U_INOD + (IPHASE-1)*Mdims%u_nonods ) = NU_ALL % VAL( IDIM, IPHASE, U_INOD )
                        END IF
                    END DO
                END DO
            END DO
        END DO
        ALLOCATE( MAT_U( Mdims%ndim * Mdims%nphase * Mdims%cv_nonods ), UD( Mdims%ndim ) ) ; mat_u = 0. ; ud = 0.
        CALL PROJ_U2MAT( COMP_DIFFUSION_OPT, &
            Mdims, CV_GIdims, CV_funs,&
            COMP_DIFFUSION, NCOMP_DIFF_COEF, COMP_DIFF_COEF, &
            X_ALL%val, NU, NV, NW, MAT_NDGLN, U_NDGLN, X_NDGLN, &
            MAT_U)
        ! Determine the diffusion coeff tensor COMP_DIFFUSION from MAT_U and COMP_DIFF_COEF
        DO MAT_NOD = 1, Mdims%mat_nonods
            DO IPHASE = 1, Mdims%nphase
                DO IDIM = 1,Mdims%ndim
                    UD( IDIM ) = MAT_U( MAT_NOD + ( IDIM - 1 ) * Mdims%mat_nonods + &
                        ( IPHASE - 1 ) * Mdims%ndim * Mdims%mat_nonods )
                END DO
                DIFF_molecular    = COMP_DIFF_COEF( 1, IPHASE )
                DIFF_longitudinal = COMP_DIFF_COEF( 2, IPHASE )
                DIFF_transverse   = COMP_DIFF_COEF( 3, IPHASE )
                CALL CALC_COMP_DIF_TEN( Mdims%ndim, UD, DIFF_molecular, DIFF_longitudinal, DIFF_transverse, &
                    COMP_DIFFUSION( MAT_NOD, : , : , IPHASE ))
            END DO
        END DO
        DEALLOCATE( MAT_U, UD )
        return
    end subroutine Calculate_ComponentDiffusionTerm

    SUBROUTINE PROJ_U2MAT( COMP_DIFFUSION_OPT, &
        Mdims, CV_GIdims, CV_funs,&
        COMP_DIFFUSION, NCOMP_DIFF_COEF, COMP_DIFF_COEF, &
        X_ALL, NU, NV, NW, MAT_NDGLN, U_NDGLN, X_NDGLN, &
        MAT_U)
        ! Determine MAT_U from NU,NV,NW which are variables mapped to material mesh.

        implicit none
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_GI_dimensions), intent(in) :: CV_GIdims
        type(multi_shape_funs), intent(in) :: CV_funs
        INTEGER, intent( in ) :: NCOMP_DIFF_COEF, COMP_DIFFUSION_OPT
        REAL, DIMENSION( :, :, :, : ), intent( in ) :: COMP_DIFFUSION
        REAL, DIMENSION( :, : ), intent( in ) :: COMP_DIFF_COEF
        REAL, DIMENSION( :, : ), intent( in ) :: X_ALL
        REAL, DIMENSION( : ), intent( in ) :: NU, NV, NW
        INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: X_NDGLN
        REAL, DIMENSION( :), intent( inout ) :: MAT_U
        ! Determine MAT_U from NU,NV,NW which are these variables mapped to material mesh.
        ! Local variables
        REAL, DIMENSION( :, : ), ALLOCATABLE :: MASS, INV_MASS, MASS2U, INV_MASS_NM
        INTEGER :: &
            ELE, MAT_ILOC, MAT_JLOC, CV_GI, U_JLOC, MAT_KLOC, MAT_NOD, &
            U_NODJ, IPHASE, U_NODJ_IP, IDIM, MAT_NOD_ID_IP, CV_GI_SHORT
        REAL :: NN, NFEMU, MASELE
        type(multi_dev_shape_funs) :: Devfuns

        ALLOCATE( MASS( Mdims%mat_nloc, Mdims%mat_nloc ))
        ALLOCATE( INV_MASS( Mdims%mat_nloc, Mdims%mat_nloc ))
        ALLOCATE( MASS2U( Mdims%mat_nloc, Mdims%u_nloc ))
        ALLOCATE( INV_MASS_NM( Mdims%mat_nloc, Mdims%u_nloc ))

        call allocate_multi_dev_shape_funs(CV_funs, Devfuns)

        Loop_Elements1: DO ELE = 1, Mdims%totele
            ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
            call DETNLXR_PLUS_U(ELE, X_ALL, X_NDGLN, CV_funs%cvweight, &
                   CV_funs%cvfen, CV_funs%cvfenlx_all, CV_funs%ufenlx_all, Devfuns)
            MASELE = 0.0
            Loop_MAT_ILOC: DO MAT_ILOC = 1, Mdims%mat_nloc
                Loop_MAT_JLOC: DO MAT_JLOC = 1, Mdims%mat_nloc
                    NN = 0.0
                    DO CV_GI_SHORT = 1, CV_GIdims%cv_ngi
                        NN = NN +  CV_funs%cvfen( MAT_ILOC, CV_GI_SHORT )  * CV_funs%cvfen(  MAT_JLOC, CV_GI_SHORT ) &
                            * DevFuns%DETWEI( CV_GI_SHORT )
                    END DO
                    MASS( MAT_ILOC,MAT_JLOC)  = MASS( MAT_ILOC,MAT_JLOC) + NN
                END DO Loop_MAT_JLOC
            END DO Loop_MAT_ILOC
            CALL MATDMATINV( MASS, INV_MASS, Mdims%mat_nloc)
            MASS2U = 0.0
            Loop_MAT_ILOC2: DO MAT_ILOC = 1, Mdims%mat_nloc
                Loop_U_JLOC: DO U_JLOC = 1, Mdims%u_nloc
                    NFEMU = 0.0
                    DO CV_GI = 1, CV_GIdims%cv_ngi
                        NFEMU = NFEMU +  CV_funs%cvfen( MAT_ILOC, CV_GI ) * CV_funs%ufen(  U_JLOC, CV_GI ) * DevFuns%DETWEI( CV_GI )
                    END DO
                    MASS2U( MAT_ILOC,U_JLOC)  = MASS2U( MAT_ILOC,U_JLOC) + NFEMU
                END DO Loop_U_JLOC
            END DO Loop_MAT_ILOC2
            INV_MASS_NM = 0.0
            DO MAT_ILOC = 1, Mdims%mat_nloc
                DO U_JLOC = 1, Mdims%u_nloc
                    DO MAT_KLOC = 1, Mdims%mat_nloc
                        INV_MASS_NM( MAT_ILOC, U_JLOC ) = INV_MASS_NM( MAT_ILOC, U_JLOC ) &
                            + INV_MASS( MAT_ILOC, MAT_KLOC ) * MASS2U( MAT_KLOC, U_JLOC )
                    END DO
                END DO
            END DO
            Loop_MAT_ILOC3: DO MAT_ILOC = 1, Mdims%mat_nloc
                MAT_NOD = MAT_NDGLN(( ELE - 1 ) * Mdims%mat_nloc + MAT_ILOC )
                Loop_U_JLOC2: DO U_JLOC = 1, Mdims%u_nloc
                    U_NODJ = U_NDGLN(( ELE - 1 ) * Mdims%u_nloc + U_JLOC )
                    Loop_IPHASE: DO IPHASE = 1, Mdims%nphase
                        U_NODJ_IP = U_NODJ + ( IPHASE - 1 ) * Mdims%u_nonods
                        IDIM = 1
                        MAT_NOD_ID_IP = MAT_NOD + ( IDIM - 1 ) * Mdims%mat_nonods  + &
                            ( IPHASE - 1 ) * Mdims%ndim * Mdims%mat_nonods
                        MAT_U( MAT_NOD_ID_IP ) = MAT_U( MAT_NOD_ID_IP ) + &
                            INV_MASS_NM( MAT_ILOC, U_JLOC ) * NU( U_NODJ_IP )
                        IF( Mdims%ndim >= 2 ) THEN
                            IDIM = 2
                            MAT_NOD_ID_IP = MAT_NOD + ( IDIM - 1 ) * Mdims%mat_nonods + &
                                ( IPHASE- 1 ) * Mdims%ndim * Mdims%mat_nonods
                            MAT_U( MAT_NOD_ID_IP ) = MAT_U( MAT_NOD_ID_IP ) + &
                                INV_MASS_NM( MAT_ILOC, U_JLOC ) * NV( U_NODJ_IP )
                        ENDIF
                        IF( Mdims%ndim >= 3 ) THEN
                            IDIM = 3
                            MAT_NOD_ID_IP = MAT_NOD + ( IDIM - 1 ) * Mdims%mat_nonods + &
                                ( IPHASE - 1 ) * Mdims%ndim * Mdims%mat_nonods
                            MAT_U( MAT_NOD_ID_IP ) = MAT_U( MAT_NOD_ID_IP ) + &
                                INV_MASS_NM( MAT_ILOC, U_JLOC ) * NW( U_NODJ_IP )
                        ENDIF
                    END DO Loop_IPHASE
                END DO Loop_U_JLOC2
            END DO Loop_MAT_ILOC3
        END DO Loop_Elements1
        ! Deallocating temporary arrays
        call deallocate_multi_dev_shape_funs(Devfuns)
        DEALLOCATE( MASS )
        DEALLOCATE( INV_MASS )
        DEALLOCATE( MASS2U )
        DEALLOCATE( INV_MASS_NM )
        RETURN
    end subroutine PROJ_U2MAT


    SUBROUTINE CALC_COMP_DIF_TEN( NDIM, UD, DIFF_molecular, DIFF_longitudinal, DIFF_transverse, &
        DIFF_TEN )
        ! Calculate the diffusion coefficient COMP_DIFFUSION for current composition...
        ! based on page 136 in Reservoir-Simulation-Mathematical-Techniques-In-Oil-Recovery-(2007).pdf
        implicit none

        INTEGER, intent( in ) :: NDIM
        REAL, intent( in ) :: DIFF_molecular, DIFF_longitudinal, DIFF_transverse
        REAL, DIMENSION( : ), intent( in ) :: UD
        REAL, DIMENSION( :, : ), intent( inout ) :: DIFF_TEN

        ! Local variables...
        REAL, PARAMETER :: TOLER = 1.0E-10
        REAL, DIMENSION( :, : ), allocatable :: E, E_ident, E_OTH
        REAL :: RN, RN2
        INTEGER :: IDIM, JDIM

        ALLOCATE( E( NDIM, NDIM ))
        ALLOCATE( E_ident( NDIM, NDIM ))
        ALLOCATE( E_OTH( NDIM, NDIM ))

        RN = 0.0
        DO IDIM = 1, NDIM
            RN = RN + UD( IDIM ) **2
        END DO

        RN2 = SQRT( RN )
        RN = MAX( RN2, TOLER )

        E_ident = 0.0
        DO IDIM = 1, NDIM

            DO JDIM = 1, NDIM
                E( IDIM, JDIM ) = UD(IDIM) * UD( JDIM ) / RN
            END DO

            E_ident( IDIM, IDIM ) = 1.0
        END DO

        E_OTH = E_ident - E

        DIFF_TEN = DIFF_molecular * E_ident  &
            + RN2 * ( DIFF_longitudinal * E + DIFF_transverse * E_OTH )


        DEALLOCATE( E )
        DEALLOCATE( E_ident )
        DEALLOCATE( E_OTH )

        RETURN

    END SUBROUTINE CALC_COMP_DIF_TEN


    subroutine Calc_KComp2( cv_nonods, nphase, icomp, KComp_Sigmoid, &
        Satura, K_Comp, max_k, min_k, &
        K_Comp2 )
        implicit none
        integer, intent( in ) :: cv_nonods, nphase, icomp
        logical, intent( in ) :: KComp_Sigmoid
        real, dimension( :, : ), intent( in ) :: Satura
        real, dimension( :, :, : ), intent( in ) :: K_Comp
        real, intent( in ) :: max_k, min_k
        real, dimension( :, :, :, : ), intent( inout ) :: K_Comp2
        ! Local variables
        integer :: iphase, jphase, cv_nod
        real, parameter :: Width = 0.1, Err = 1.e-6, Sat = 0.9
        real :: Sat0

        K_Comp2 = 0.

        Conditional_KComp_Sig: if( .not. KComp_Sigmoid ) then

            do cv_nod = 1, cv_nonods
                do iphase = 1, nphase
                    do jphase = iphase + 1, nphase, 1
                        K_Comp2( icomp, cv_nod, iphase, jphase ) = &
                            1. / K_Comp( icomp, iphase, jphase )
                    end do
                end do
            end do

        else

            Sat0 = Sat
            do cv_nod = 1, cv_nonods
                do iphase = 1, nphase
                    do jphase = 1, nphase
                        if ( jphase /= iphase ) then
                            K_Comp2( icomp, cv_nod, iphase, jphase ) = &
                                1. / sigmoid_function( satura(iphase,cv_nod ), &
                                Sat0, Width, min_k, max_k )
                        endif
                    end do
                end do
            end do

        end if Conditional_KComp_Sig

        return
    end subroutine Calc_KComp2


    real function sigmoid_function( Y, Y0, Width, LowMag, UpMag )
        implicit none
        real :: Y, Y0, Width, LowMag, UpMag
        ! Local Variables
        real :: alpha
        !
        ! Width: width of the sigmoid function.
        ! The sigmoid function, varies between ( LowMag, UpMag ).
        ! Y is the variable of the function and Y0 is the centre
        ! of the function.
        ! The function looks like:
        !             -------------
        !           /
        !          /
        ! --------

        if( Y - Y0 < - 3. * Width ) then
            sigmoid_function = UpMag
        elseif( Y - Y0 > 3. * Width ) then
            sigmoid_function = LowMag
        else
            alpha = 10. / Width
            sigmoid_function = ( UpMag - LowMag ) / &
                ( 1. + exprep( alpha * ( Y - Y0 ))) + LowMag
        end if

        return
    end function sigmoid_function

    real function exprep( M )
        implicit none
        real :: M

        if( M > 174. ) then
            exprep = 3.69e+35
        elseif( M < -180. ) then
            exprep = 0.
        else
            exprep = exp( M )
        end if

        return
    end function exprep



    SUBROUTINE CAL_COMP_SUM2ONE_SOU( packed_state, Mdims )

      ! make sure the composition sums to 1.0
      implicit none
      type( state_type ), intent( inout ) :: packed_state
      type( multi_dimensions ), intent( in ) :: Mdims

      ! the relaxing (sum2one_relax) is to help convergence.
      ! =1 is full adjustment to make sure we have sum to 1.
      ! =0 is no adjustment.
      real :: dt, sum2one_relax, comp_sum
      integer :: cv_nonods, nphase, ncomp2, iphase, cv_nodi
      logical :: ensure_positive
      !Working pointer
      real, dimension(:,:), pointer ::satura
      type( vector_field ), pointer :: MeanPoreCV
      type( tensor_field ), pointer :: MFC_s, tracer_source

      cv_nonods = Mdims%cv_nonods
      nphase = Mdims%nphase
      ncomp2 = Mdims%ncomp

      call get_option( '/timestepping/timestep', dt )

      MeanPoreCV => extract_vector_field( packed_state, "MeanPoreCV" )

      call get_var_from_packed_state(packed_state,PhaseVolumeFraction = satura)
      MFC_s  => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedComponentMassFraction" )
      tracer_source => extract_tensor_field(packed_state, "PackedPhaseVolumeFractionComponentSource")

      if( have_option( '/material_phase[' // int2str( nphase ) // &
           ']/is_multiphase_component/Comp_Sum2One/Relaxation_Coefficient' ) ) then
         call get_option( '/material_phase[' // int2str( nphase ) // &
              ']/is_multiphase_component/Comp_Sum2One/Relaxation_Coefficient', &
              sum2one_relax )

         ensure_positive = have_option( '/material_phase[' // int2str( nphase ) // &
              ']/is_multiphase_component/Comp_Sum2One/Ensure_Positive' )
      else
         FLAbort( 'Please define the relaxation coefficient for components mass conservation constraint.' )
      end if

      ewrite(3,*) 'sum2one_relax, ensure_positive', sum2one_relax, ensure_positive

      DO IPHASE = 1, NPHASE
         DO CV_NODI = 1, CV_NONODS

            COMP_SUM = SUM (MFC_s % val (:, IPHASE, CV_NODI) )

            IF ( ENSURE_POSITIVE ) THEN
               tracer_source%val(1, iphase, cv_nodi) = tracer_source%val(1, iphase, cv_nodi) &
                    - SUM2ONE_RELAX * MeanPoreCV%val( 1, CV_NODI ) * SATURA( IPHASE, CV_NODI ) * MAX( ( 1. - COMP_SUM ), 0. ) / DT
            ELSE
               tracer_source%val(1, iphase, cv_nodi) = tracer_source%val(1, iphase, cv_nodi) &
                    - SUM2ONE_RELAX * MeanPoreCV%val( 1, CV_NODI ) * SATURA( IPHASE, CV_NODI ) * ( 1. - COMP_SUM ) / DT
            END IF

         END DO
      END DO

      RETURN
    END SUBROUTINE CAL_COMP_SUM2ONE_SOU

end module Compositional_Terms
