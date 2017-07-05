
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

!!!!========================================!!!!
!!!!           MATRIX OPERATIONS            !!!!
!!!!========================================!!!!
module matrix_operations

    use fldebug
    use spud
    use Fields_Allocates, only : allocate, make_mesh
    use fields_data_types, only: mesh_type, scalar_field, vector_field, tensor_field
    use fields, only : node_count, set
    use state_module
    use sparse_tools
    use sparse_tools_petsc, only : petsc_csr_matrix, petsc_csr_matrix_pointer, &
        allocate, deallocate, &
        entries,  &
        zero, addto, addto_diag, scale, &
        extract_diagonal, assemble, incref_petsc_csr_matrix, &
        addref_petsc_csr_matrix, &
        mult_T, dump_matrix, &
        csr2petsc_csr, dump_petsc_csr_matrix
    use parallel_tools
    use parallel_fields
    use halos
    use petsc_tools
    use petsc
    use global_parameters
    use boundary_conditions
    use multi_data_types
    use multi_tools
    implicit none

contains


    SUBROUTINE MATDMATINV( DMAT, DMATINV, NLOC )
        ! calculate DMATINV
        IMPLICIT NONE
        INTEGER, intent( in ) :: NLOC
        REAL, DIMENSION( :, : ), intent( in ) :: DMAT
        REAL, DIMENSION( :, : ), intent( inout ) ::  DMATINV
        ! Local variables
        !  REAL, DIMENSION( NLOC , NLOC ) :: MAT, MAT2
        !  REAL, DIMENSION( NLOC ) :: X, B

        DMATINV = DMAT
        CALL MATINV( DMATINV, NLOC, NLOC)!, MAT, MAT2, X, B )

        RETURN
    END SUBROUTINE MATDMATINV
    !


    SUBROUTINE MATINV( A, N, NMAX)!, MAT, MAT2, X, B)
        ! This sub finds the inverse of the matrix A and puts it back in A.
        ! MAT, MAT2, X and B are working vectors.
        IMPLICIT NONE
        INTEGER, intent( in ) :: N, NMAX
        REAL, DIMENSION( :, : ), intent( inout ) ::  A
        !      REAL, DIMENSION( N, N ), intent( inout ) :: MAT, MAT2
        !      REAL, DIMENSION( N ), intent( inout ) :: X, B
        ! Local variables
        !      INTEGER , DIMENSION(N,N) :: IPIV
        INTEGER , DIMENSION(N ) :: IPIV
        REAL, DIMENSION( N, N ) :: MAT

        real, dimension(max(1,N*N)) :: WORK
        integer :: LWORK

        integer info

        interface
            subroutine dgetrf(N,M,A,NMAX,IPIV,INFO)
                integer :: N, M, NMAX, info
                real, dimension(NMAX,M) :: A
                integer, dimension(N) :: IPIV
            end subroutine dgetrf
        end interface

        interface
            subroutine dgetri(N,A,NMAX,IPIV,WORK,LWORK,INFO)
                integer :: N, NMAX, info, LWORK
                real, dimension(NMAX,N) :: A
                integer, dimension(N) :: IPIV
                real, dimension(LWORK) :: WORK
            end subroutine dgetri
        end interface

        LWORK=N*N
      
        mat=a

        call dgetrf(N,N,A,NMAX,IPIV,INFO)
        call dgetri(N,A,NMAX,IPIV,WORK,max(1,LWORK),INFO)
        if (info==0) then
            return
        else
            ewrite(2,*) MAT
            print *, 'mat', mat
            FLAbort("PIVIT Matrix block inversion failed")
        end if

    END SUBROUTINE MATINV

    SUBROUTINE MATINVold( A, N, MAT,B)
        ! This sub finds the inverse of the matrix A and puts it back in A.
        ! MAT, MAT2, X and B are working vectors.
        IMPLICIT NONE
        INTEGER, intent( in ) :: N
        REAL, DIMENSION( :, : ), intent( inout ) ::  A
        REAL, DIMENSION( :, : ), intent( inout ) :: MAT
        REAL, DIMENSION( : ), intent( inout ) :: B
        ! Local variables
        INTEGER :: ICOL
        INTEGER , DIMENSION(N) :: IPIV

        ! Solve MAT XL=BX (NB BDIAG is overwritten)
        ICOL = 1
        B( 1 : N ) = 0.
        B( ICOL ) = 1.

        MAT = A( 1:N,1:N )

        CALL SMLINNGOT( MAT, A( :, ICOL ), B, N, IPIV, .FALSE. ) ! X contains the column ICOL of inverse

        DO ICOL = 2, N ! Form column ICOL of the inverse.
            B = 0.
            B( ICOL ) = 1.0 ! Solve MAT X=B (NB MAT is overwritten).

            CALL SMLINNGOT( MAT, A( :, ICOL ), B, N, IPIV, .TRUE. ) ! X contains the column ICOL of inverse

        END DO

        RETURN
    END SUBROUTINE MATINVold

    SUBROUTINE SMLINNGOT( A, X, B, NMX, IPIV, GOTDEC )
    !Calculate the inverse using the LU decomposition
    !L can be provided, speeding up the method to O(n)
        IMPLICIT NONE
        INTEGER :: NMX
        REAL, DIMENSION( :, : ), intent( inout ) :: A
        real, DIMENSION( : ), intent( inout ) :: X ! inout as n might be < nmx
        real, DIMENSION( : ), intent( in ) ::  B
        LOGICAL, intent( in ) :: GOTDEC
        ! Local
        INTEGER , DIMENSION(:), intent(inout) :: IPIV
        INTEGER :: INFO

        real, DIMENSION( NMX,1 ) :: Bloc

        interface
            subroutine dgetrf(M,N,A, LDA, IPIV, INFO)
                implicit none
                integer :: info,lda, m,n
                integer , dimension(min(m,n)) :: IPIV
                real, dimension(LDA,N) :: A
            end subroutine dgetrf
        end interface

        interface
            subroutine dgetrs(TRANS,N,NRHS,A, LDA, IPIV, B,LDB, INFO)
                implicit none
                character(len=1) :: TRANS
                integer :: info,lda,ldb,n,nrhs
                integer , dimension(n) :: IPIV
                real, dimension(LDA,N) :: A
                real, dimension(LDA,NRHS) :: B
            end subroutine dgetrs
        end interface

         ! IF GOTDEC then assume we have already got the LU decomposition in A
         ! Form X = A^{-1} B ;  Useful subroutine for inverse
         ! This sub overwrites the matrix A.

        Cond_LUDecomp: IF( .NOT. GOTDEC ) THEN

            call dgetrf(NMX,NMX,A,NMX,IPIV,INFO)

        ENDIF Cond_LUDecomp

        Bloc(:,1)=B
        call dgetrs('N',NMX,1,A, NMX,IPIV, Bloc, NMX, INFO)
        X=Bloc(:,1)

        RETURN
    END SUBROUTINE SMLINNGOT


    SUBROUTINE COLOR_GET_CMC_PHA( Mdims, Mspars, ndgln, Mmat,&
        DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
        CMC_petsc, CMC_PRECON, IGOT_CMC_PRECON, MASS_MN_PRES, &
        pipes_aux, got_free_surf,  MASS_SUF, &
        symmetric_P )
        !use multiphase_1D_engine
        !Initialize the momentum equation (CMC) and introduces the corresponding values in it.
        implicit none
        ! form pressure matrix CMC using a colouring approach
        type(multi_dimensions), intent(in) :: Mdims
        type (multi_sparsities), intent(in) :: Mspars
        type(multi_ndgln), intent(in) :: ndgln
        type (multi_matrices), intent(inout) :: Mmat
        LOGICAL, PARAMETER :: PIPES_1D=.TRUE.
        INTEGER, intent( in ) :: IGOT_CMC_PRECON
        LOGICAL, intent( in ) :: got_free_surf, symmetric_P
        REAL, DIMENSION( :, : ), intent( in ) :: DIAG_SCALE_PRES
        REAL, DIMENSION( :, :, : ), intent( in ) :: DIAG_SCALE_PRES_COUP, INV_B
        type(petsc_csr_matrix), intent(inout)::  CMC_petsc
        REAL, DIMENSION( :, :, : ), intent( inout ) :: CMC_PRECON
        REAL, DIMENSION( : ), intent( in ) :: MASS_MN_PRES
        REAL, DIMENSION( : ), intent( in ) :: MASS_SUF
        type (multi_pipe_package), intent(in) :: pipes_aux
        ! Local variables
        REAL, PARAMETER :: INFINY = 1.0E+10
        INTEGER, DIMENSION( : ), allocatable :: ndpset
        !Initialize CMC_petsc
        call zero( CMC_petsc )
        allocate(ndpset(Mdims%npres)) ; ndpset(:)=-1
        if (ndpset(1)<0) then
            call get_option( '/material_phase[0]/scalar_field::Pressure/' // &
                'prognostic/reference_node', ndpset(1), default = -1 )
            ndpset(2:Mdims%npres)=ndpset(1)
        endif
        if (isparallel()) then
            if (GetProcNo()>1) ndpset=0
        end if

        !COLOR_GET_CMC_PHA_FAST is very memory hungry, so we let the user decide
        !or if we are using a compacted lumped mass matrix then the memory reduction compensates this extra memory usage
        IF ( Mdims%npres==1 .and.( have_option("/numerical_methods/create_P_mat_fast") .or. size(Mmat%PIVIT_MAT,1) == 1 )) THEN
            ! Fast but memory intensive... (wells not yet implemented here)
            CALL COLOR_GET_CMC_PHA_FAST( Mdims,Mspars, ndgln, Mmat,  &
                DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
                CMC_petsc, CMC_PRECON, IGOT_CMC_PRECON, MASS_MN_PRES, &
                pipes_aux%MASS_PIPE, pipes_aux%MASS_CVFEM2PIPE, pipes_aux%MASS_CVFEM2PIPE_TRUE, &
                got_free_surf,  MASS_SUF, ndpset, symmetric_P )
        ELSE
            ! Slow but memory efficient...
            CALL COLOR_GET_CMC_PHA_SLOW( Mdims,Mspars, ndgln, Mmat,&
                DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
                CMC_petsc, CMC_PRECON, IGOT_CMC_PRECON, MASS_MN_PRES, &
                pipes_aux%MASS_PIPE, pipes_aux%MASS_CVFEM2PIPE, pipes_aux%MASS_CVFEM2PIPE_TRUE, &
                got_free_surf,  MASS_SUF, ndpset, symmetric_P )
        END IF
        !Re-assemble just in case
        CMC_petsc%is_assembled=.false.
        call assemble( CMC_petsc )
    contains


       SUBROUTINE COLOR_GET_CMC_PHA_SLOW( Mdims, Mspars, ndgln, Mmat,  &
            DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
            CMC_petsc, CMC_PRECON, IGOT_CMC_PRECON, MASS_MN_PRES, &
            MASS_PIPE, MASS_CVFEM2PIPE, MASS_CVFEM2PIPE_TRUE, &
            got_free_surf,  MASS_SUF, ndpset, symmetric_P )
            !use multiphase_1D_engine
            implicit none
            ! form pressure matrix CMC using a colouring approach
            type(multi_dimensions), intent(in) :: Mdims
            type (multi_sparsities), intent(in) :: Mspars
            type(multi_ndgln), intent(in) :: ndgln
            type (multi_matrices), intent(inout) :: Mmat
            INTEGER, intent( in ) :: IGOT_CMC_PRECON
            LOGICAL, PARAMETER :: PIPES_1D=.TRUE.
            LOGICAL, intent( in ) :: got_free_surf, symmetric_P
            INTEGER, DIMENSION( : ), intent( in ) :: ndpset
            REAL, DIMENSION( :, : ), intent( in ) :: DIAG_SCALE_PRES
            REAL, DIMENSION( :, :, : ), intent( in ) :: DIAG_SCALE_PRES_COUP, INV_B
            type(petsc_csr_matrix), intent(inout)::  CMC_petsc
            REAL, DIMENSION( :, :, : ), intent( inout ) :: CMC_PRECON
            REAL, DIMENSION( : ), intent( in ) :: MASS_MN_PRES
            REAL, DIMENSION( : ), intent( in ) :: MASS_SUF
            REAL, DIMENSION( : ), intent( in ) :: MASS_PIPE, MASS_CVFEM2PIPE, MASS_CVFEM2PIPE_TRUE
            ! Local variables
            INTEGER, PARAMETER :: MX_NCOLOR = 1000
            REAL, PARAMETER :: INFINY = 1.0E+10
            LOGICAL :: UNDONE, LCOL
            logical, DIMENSION( : ), allocatable :: NEED_COLOR
            logical, DIMENSION( Mdims%cv_nonods ) :: to_color
            LOGICAL :: EXPLICIT_PIPES2
            REAL, DIMENSION( : ), allocatable :: COLOR_VEC
            REAL, DIMENSION( :, : ), allocatable :: CMC_COLOR_VEC, CMC_COLOR_VEC2, CMC_COLOR_VEC_PHASE, CMC_COLOR_VEC2_PHASE!, ld
            REAL, DIMENSION( Mdims%ndim * Mdims%nphase * Mdims%u_nonods ) :: DU_LONG
            real, dimension(3 * Mdims%nphase * Mdims%u_nonods), target :: temp_memory
            real, dimension(:), pointer :: DU, DV, DW
            REAL, DIMENSION( :, :, : ), pointer :: CDP
            INTEGER :: NCOLOR, CV_NOD, CV_JNOD, COUNT, COUNT2, IPHASE, CV_JNOD2
            INTEGER :: ierr, IV_STAR, IV_FINI, i_indx, j_indx, IPRES, JPRES
            REAL :: RSUM, RSUM_SUF
            ALLOCATE( NEED_COLOR( Mdims%cv_nonods ) )
            ALLOCATE( COLOR_VEC( Mdims%cv_nonods ) )
            ALLOCATE( CMC_COLOR_VEC( Mdims%npres, Mdims%cv_nonods ) )
            ALLOCATE( CMC_COLOR_VEC2( Mdims%npres, Mdims%cv_nonods ) )

            !CDP and DU, DV and DW can share memory as they never occur at the same time
            CDP(1:Mdims%ndim, 1:Mdims%nphase, 1:Mdims%u_nonods) => temp_memory(1:Mdims%ndim * Mdims%nphase * Mdims%u_nonods)
            DU(1:Mdims%u_nonods * Mdims%nphase) => temp_memory(1:    Mdims%u_nonods * Mdims%nphase)
            DV(1:Mdims%u_nonods * Mdims%nphase) => temp_memory(1 +   Mdims%u_nonods * Mdims%nphase:2*Mdims%u_nonods * Mdims%nphase)
            DW(1:Mdims%u_nonods * Mdims%nphase) => temp_memory(1 + 2*Mdims%u_nonods * Mdims%nphase:3*Mdims%u_nonods * Mdims%nphase)

            IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON = 0.0
            EXPLICIT_PIPES2 = .true.
            NEED_COLOR = .TRUE.
            NCOLOR = 0
            UNDONE = .TRUE.
            Loop_while: DO WHILE ( UNDONE )
                NCOLOR = NCOLOR + 1 ! Determine what nodes can be coloured with the new color
                TO_COLOR = .FALSE.
                Loop_CVNOD: DO CV_NOD = 1, Mdims%cv_nonods
                    IF ( NEED_COLOR( CV_NOD ) ) THEN
                        LCOL= .FALSE.
                        ! use a distance-2 colouring...
                        Loop_Row: DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                            CV_JNOD = Mspars%CMC%col( COUNT )
                            IF( TO_COLOR( CV_JNOD ) ) THEN
                                LCOL=.TRUE.
                                EXIT
                            END IF
                            Loop_Row2: DO COUNT2 = Mspars%CMC%fin( CV_JNOD ), Mspars%CMC%fin( CV_JNOD + 1 ) - 1
                                IF ( TO_COLOR( Mspars%CMC%col( COUNT2 ) ) )  THEN
                                    LCOL = .TRUE.
                                    EXIT
                                END IF
                            END DO Loop_Row2
                            IF ( LCOL ) EXIT
                        END DO Loop_Row
                        IF ( .NOT.LCOL ) THEN
                            TO_COLOR( CV_NOD ) = .TRUE.
                        END IF
                    END IF
                END DO Loop_CVNOD
                NEED_COLOR = NEED_COLOR .AND. .NOT.TO_COLOR
                COLOR_VEC = MERGE( 1.0, 0.0, TO_COLOR )
                CALL C_MULT2( CDP, COLOR_VEC, Mdims%cv_nonods, Mdims%u_nonods, Mdims%ndim, Mdims%nphase, &
                    Mmat%C, Mspars%C%ncol, Mspars%C%fin, Mspars%C%col )

                ! DU_LONG = BLOCK_MAT * CDP
                CALL PHA_BLOCK_MAT_VEC( DU_LONG, Mmat%PIVIT_MAT, CDP, Mdims%u_nonods, Mdims%ndim, Mdims%nphase, &
                    Mdims%totele, Mdims%u_nloc, ndgln%u )

                ! NB. P_RHS = Mmat%CT * U + CV_RHS
                ! DU_LONG = CDP
                CALL ULONG_2_UVW( DU, DV, DW, DU_LONG, Mdims%u_nonods, Mdims%ndim, Mdims%nphase )
                IF ( Mdims%npres > 1 .AND. .NOT.EXPLICIT_PIPES2 ) THEN
                    ALLOCATE( CMC_COLOR_VEC_PHASE( Mdims%nphase, Mdims%cv_nonods ) )
                    DO IPHASE = 1, Mdims%nphase
                        IV_STAR = 1+(IPHASE-1)*Mdims%u_nonods*Mdims%ndim
                        IV_FINI = IPHASE*Mdims%u_nonods*Mdims%ndim
                        CALL CT_MULT( CMC_COLOR_VEC_PHASE(IPHASE,:), DU(IV_STAR:IV_FINI), &
                            DV(IV_STAR:IV_FINI), DW(IV_STAR:IV_FINI), Mdims%cv_nonods, Mdims%u_nonods, Mdims%ndim, 1, &
                            Mmat%CT(:,IPHASE:IPHASE,:), Mspars%CT%ncol, Mspars%CT%fin, Mspars%CT%col )
                    END DO
                    DO CV_NOD = 1, Mdims%cv_nonods
                        CMC_COLOR_VEC_PHASE(:,CV_NOD) = MATMUL( INV_B(:,:,CV_NOD), CMC_COLOR_VEC_PHASE(:,CV_NOD) )
                    END DO
                    DO CV_NOD = 1, Mdims%cv_nonods
                        DO IPRES = 1, Mdims%npres
                            CMC_COLOR_VEC(IPRES,CV_NOD) = SUM(CMC_COLOR_VEC_PHASE(1+(IPRES-1)*Mdims%n_in_pres:IPRES*Mdims%n_in_pres,CV_NOD) )
                        END DO
                    END DO
                    deallocate(CMC_COLOR_VEC_PHASE)
                    IF ( IGOT_CMC_PRECON /= 0 ) THEN
                        ALLOCATE( CMC_COLOR_VEC2_PHASE( Mdims%nphase, Mdims%cv_nonods ) )
                        DO IPHASE = 1, Mdims%nphase
                            IV_STAR = 1+(IPHASE-1)*Mdims%u_nonods*Mdims%ndim
                            IV_FINI = IPHASE*Mdims%u_nonods*Mdims%ndim
                            CALL CT_MULT_WITH_C( CMC_COLOR_VEC2_PHASE(IPHASE,:), &
                                DU_LONG(IV_STAR:IV_FINI), &
                                Mdims%u_nonods, Mdims%ndim, 1,  &
                                Mmat%C(:,IPHASE:IPHASE,:), Mspars%C%ncol, Mspars%C%fin, Mspars%C%col )
                        END DO
                        DO CV_NOD = 1, Mdims%cv_nonods
                            CMC_COLOR_VEC2_PHASE(:,CV_NOD) = MATMUL( INV_B(:,:,CV_NOD), CMC_COLOR_VEC2_PHASE(:,CV_NOD) )
                        END DO
                        DO CV_NOD = 1, Mdims%cv_nonods
                            DO IPHASE = 1, Mdims%nphase
                                CMC_COLOR_VEC2(IPHASE,CV_NOD) = SUM(CMC_COLOR_VEC2_PHASE(IPHASE:IPHASE,CV_NOD) )
                            END DO
                        END DO
                        deallocate(CMC_COLOR_VEC2_PHASE)
                    END IF
                ELSE
                    DO IPRES = 1, Mdims%npres
                        IV_STAR = 1+(IPRES-1)*Mdims%n_in_pres*Mdims%u_nonods
                        IV_FINI = IPRES*Mdims%n_in_pres*Mdims%u_nonods
                        CALL CT_MULT( CMC_COLOR_VEC(IPRES,:), DU(IV_STAR:IV_FINI), &
                            DV(IV_STAR:IV_FINI), DW(IV_STAR:IV_FINI), Mdims%cv_nonods, Mdims%u_nonods, Mdims%ndim, Mdims%n_in_pres, &
                            Mmat%CT(:,1+(IPRES-1)*Mdims%n_in_pres:IPRES*Mdims%n_in_pres,:), Mspars%CT%ncol, Mspars%CT%fin, Mspars%CT%col )
                    END DO
                    IF ( IGOT_CMC_PRECON /= 0 ) THEN
                        DO IPRES = 1, Mdims%npres
                            IV_STAR = 1+(IPRES-1)*Mdims%n_in_pres*Mdims%u_nonods*Mdims%ndim
                            IV_FINI = IPRES*Mdims%n_in_pres*Mdims%u_nonods*Mdims%ndim
                            CALL CT_MULT_WITH_C( CMC_COLOR_VEC2(IPRES,:), &
                                DU_LONG(IV_STAR:IV_FINI), Mdims%u_nonods, Mdims%ndim, Mdims%n_in_pres, &
                                Mmat%C(:,1+(IPRES-1)*Mdims%n_in_pres:IPRES*Mdims%n_in_pres,:), Mspars%C%ncol, Mspars%C%fin, Mspars%C%col )
                        END DO
                    END IF
                END IF
                ! Matrix vector involving the mass diagonal term
                DO CV_NOD = 1, Mdims%cv_nonods
                    RSUM = 0.0
                    RSUM_SUF = 0.0
                    DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                        CV_JNOD = Mspars%CMC%col( COUNT )
                        DO IPRES = 1, Mdims%npres
                            IF((Mdims%npres>1).AND.PIPES_1D) THEN
                                IF(IPRES==1) THEN
                                    CMC_COLOR_VEC( IPRES, CV_NOD ) = CMC_COLOR_VEC( IPRES, CV_NOD ) &
                                        + DIAG_SCALE_PRES( IPRES, CV_NOD ) * MASS_MN_PRES( COUNT ) * COLOR_VEC( CV_JNOD )
                                ELSE
                                    CMC_COLOR_VEC( IPRES, CV_NOD ) = CMC_COLOR_VEC( IPRES, CV_NOD ) &
                                        + DIAG_SCALE_PRES( IPRES, CV_NOD ) * MASS_CVFEM2PIPE_TRUE( COUNT ) * COLOR_VEC( CV_JNOD )
                                ENDIF
                            ELSE
                                CMC_COLOR_VEC( IPRES, CV_NOD ) = CMC_COLOR_VEC( IPRES, CV_NOD ) &
                                    + DIAG_SCALE_PRES( IPRES, CV_NOD ) * MASS_MN_PRES( COUNT ) * COLOR_VEC( CV_JNOD )
                            ENDIF
                            if ( got_free_surf) then
                                CMC_COLOR_VEC( IPRES, CV_NOD ) = CMC_COLOR_VEC( IPRES, CV_NOD ) +&
                                    MASS_SUF( COUNT ) * COLOR_VEC( CV_JNOD )
                            end if
                        END DO
                        if ( got_free_surf ) then
                            RSUM_SUF = RSUM_SUF + MASS_SUF( COUNT )
                        end if
                        RSUM = RSUM + MASS_MN_PRES( COUNT )
                    END DO
                    IF ( IGOT_CMC_PRECON /= 0 ) THEN ! Use lumping of MASS_MN_PRES & MASS_SUF...
                        DO IPRES = 1, Mdims%npres
                            CMC_COLOR_VEC2( IPRES, CV_NOD ) = CMC_COLOR_VEC2( IPRES, CV_NOD ) &
                                + (DIAG_SCALE_PRES( IPRES, CV_NOD ) * RSUM + RSUM_SUF) * COLOR_VEC( CV_NOD )
                        END DO
                    END IF
                END DO
                !Put into matrix CMC
                if ( .not.symmetric_P ) then
                    ! original method
                    DO CV_NOD = 1, Mdims%cv_nonods
                        DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                            CV_JNOD = Mspars%CMC%col( COUNT )
                            DO IPRES = 1, Mdims%npres
                                JPRES = IPRES ! Add contributions to the block diagonal only.
                                call addto( CMC_petsc, blocki = IPRES, blockj = JPRES, i = cv_nod, j = CV_JNOD, &
                                    val = CMC_COLOR_VEC( IPRES, CV_NOD ) * COLOR_VEC( CV_JNOD ) )
                                !CMC( COUNT ) = CMC( COUNT ) + sum(CMC_COLOR_VEC_MANY( :, CV_NOD ) * COLOR_VEC_MANY( :, CV_JNOD ))
                                IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON( IPRES, JPRES, COUNT ) = CMC_PRECON( IPRES, JPRES, COUNT ) + &
                                    CMC_COLOR_VEC2( IPRES, CV_NOD ) * COLOR_VEC( CV_JNOD )
                            END DO
                        END DO
                    END DO
                else  ! endof if ( .not.symmetric_P ) then
                    ! symmetric P matrix
                    DO CV_NOD = 1, Mdims%cv_nonods
                        DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                            CV_JNOD = Mspars%CMC%col( COUNT )
                            DO IPRES = 1, Mdims%npres
                                JPRES = IPRES ! Add contributions to the block diagonal only.
                                call addto( CMC_petsc, blocki = IPRES, blockj = JPRES, i = cv_nod, j = CV_JNOD, &
                                    val = CMC_COLOR_VEC2( IPRES, CV_NOD ) * COLOR_VEC( CV_JNOD ) )
                                !CMC( COUNT ) = CMC( COUNT ) + sum(CMC_COLOR_VEC_MANY( :, CV_NOD ) * COLOR_VEC_MANY( :, CV_JNOD ))
                                IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON( IPRES, JPRES, COUNT ) = CMC_PRECON( IPRES, JPRES, COUNT ) + &
                                    CMC_COLOR_VEC2( IPRES, CV_NOD ) * COLOR_VEC( CV_JNOD )
                            END DO
                        END DO
                    END DO
                end if  ! endof if ( .not.symmetric_P ) then else
                UNDONE = ANY( NEED_COLOR )
                ewrite(3,*)'************ rsum,undone,NCOLOR=', rsum, undone, NCOLOR
            END DO Loop_while
            DEALLOCATE( COLOR_VEC )
            ! the matrix coupling term place immediately into matrix and not through colouring...
            ! Matrix vector involving the mass diagonal term
            IF(Mdims%npres > 1) THEN
                DO CV_NOD = 1, Mdims%cv_nonods
                    RSUM = 0.0
                    RSUM_SUF = 0.0
                    DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                        CV_JNOD = Mspars%CMC%col( COUNT )
                        DO IPRES = 1, Mdims%npres
                            DO JPRES = 1, Mdims%npres
                                IF(PIPES_1D) THEN
                                    call addto( CMC_petsc, blocki = IPRES, blockj = JPRES, i = cv_nod, j = CV_JNOD, &
                                        val = DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_NOD ) * MASS_CVFEM2PIPE( COUNT ))
                                    IF ( IGOT_CMC_PRECON /= 0 ) THEN ! Use lumping of MASS_MN_PRES & MASS_SUF...
                                        CMC_PRECON( IPRES, JPRES, COUNT ) = CMC_PRECON( IPRES, JPRES, COUNT ) &
                                            + sqrt(DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_NOD )) * &
                                            MASS_CVFEM2PIPE( COUNT )*sqrt(DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_JNOD ))
                                    END IF
                                ELSE
                                    call addto( CMC_petsc, blocki = IPRES, blockj = JPRES, i = cv_nod, j = CV_JNOD, &
                                        val = DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_NOD ) * MASS_MN_PRES( COUNT ))
                                    IF ( IGOT_CMC_PRECON /= 0 ) THEN ! Use lumping of MASS_MN_PRES & MASS_SUF...
                                        CMC_PRECON( IPRES, JPRES, COUNT ) = CMC_PRECON( IPRES, JPRES, COUNT ) &
                                            + sqrt(DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_NOD )) * &
                                            MASS_MN_PRES( COUNT )*sqrt(DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_JNOD ))
                                    END IF
                                ENDIF
                            END DO
                        END DO
                    END DO
                END DO
            END IF ! ENDOF IF(Mdims%npres > 1) THEN
            !If we have a reference node with pressure zero we impose that here.
            DO IPRES = 1, Mdims%npres
                IF ( NDPSET(IPRES) > 0 ) THEN
                    CV_NOD = NDPSET( IPRES )
                    i_indx = CMC_petsc%row_numbering%gnn2unn( cv_nod,ipres)
                    DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                        CV_JNOD = Mspars%CMC%col( COUNT )
                        IF ( CV_JNOD /= CV_NOD ) THEN
                            do jpres = 1, Mdims%npres
                                j_indx = CMC_petsc%column_numbering%gnn2unn( cv_jnod, jpres )
                                call MatSetValue(CMC_petsc%M, i_indx, j_indx, 0.0,INSERT_VALUES, ierr) ! not the diagonal
                                !CMC( COUNT ) = 0.0 ! not the diagonal
                                IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON( IPRES, JPRES, COUNT ) = 0.0
                                DO COUNT2 = Mspars%CMC%fin( CV_JNOD ), Mspars%CMC%fin( CV_JNOD + 1 ) - 1
                                    CV_JNOD2 = Mspars%CMC%col( COUNT2 )
                                    IF ( CV_JNOD2 == CV_NOD ) then
                                        i_indx = CMC_petsc%row_numbering%gnn2unn( cv_jnod, ipres )
                                        j_indx = CMC_petsc%column_numbering%gnn2unn( CV_JNOD2, jpres )
                                        call MatSetValue(CMC_petsc%M, i_indx, j_indx, 0.0, INSERT_VALUES, ierr) ! not the diagonal
                                        IF ( IGOT_CMC_PRECON/=0 ) CMC_PRECON( ipres, jpres, COUNT2 ) = 0.0
                                    END IF
                                END DO
                            end do ! do jpres=1,Mdims%npres
                        END IF
                    END DO
                END IF
            END DO

            if (Mdims%npres > 1) then
                DO CV_NOD = 1, Mdims%cv_nonods
                    if ( mass_pipe(cv_nod) == 0.0 ) then
                        CV_JNOD = CV_NOD
                        DO IPRES = 2, Mdims%npres
                            JPRES = IPRES
                            i_indx = CMC_petsc%row_numbering%gnn2unn( cv_nod, ipres )
                            j_indx = CMC_petsc%column_numbering%gnn2unn( CV_JNOD, jpres )
                            call MatSetValue(CMC_petsc%M, i_indx, j_indx, 1.0, INSERT_VALUES, ierr)
                        END DO
                    end if
                END DO
            end if
!            if ( .false. ) then
!                cv_nod = 226
!                IPRES = 2
!                i_indx = CMC_petsc%row_numbering%gnn2unn( cv_nod,ipres)
!                DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
!                    CV_JNOD = Mspars%CMC%col( COUNT )
!                    jpres = ipres
!                    j_indx = CMC_petsc%column_numbering%gnn2unn( cv_jnod,jpres)
!                    IF ( CV_JNOD /= CV_NOD ) THEN
!                        call MatSetValue(CMC_petsc%M, i_indx, j_indx, 0.0,INSERT_VALUES, ierr)
!                    ELSE
!                        call MatSetValue(CMC_petsc%M, i_indx, j_indx, 1.0,INSERT_VALUES, ierr)
!                    END IF
!                END DO
!                cv_nod = 151
!                IPRES = 2
!                i_indx = CMC_petsc%row_numbering%gnn2unn( cv_nod,ipres)
!                DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
!                    CV_JNOD = Mspars%CMC%col( COUNT )
!                    jpres = ipres
!                    j_indx = CMC_petsc%column_numbering%gnn2unn( cv_jnod,jpres)
!                    IF ( CV_JNOD /= CV_NOD ) THEN
!                        call MatSetValue(CMC_petsc%M, i_indx, j_indx, 0.0,INSERT_VALUES, ierr)
!                    ELSE
!                        call MatSetValue(CMC_petsc%M, i_indx, j_indx, 1.0,INSERT_VALUES, ierr)
!                    END IF
!                END DO
!            end if
            CMC_petsc%is_assembled = .false.
            call assemble( CMC_petsc )
            DEALLOCATE( NEED_COLOR )
            DEALLOCATE( CMC_COLOR_VEC )
            DEALLOCATE( CMC_COLOR_VEC2 )
            RETURN
        END SUBROUTINE COLOR_GET_CMC_PHA_SLOW

        SUBROUTINE COLOR_GET_CMC_PHA_FAST( Mdims, Mspars, ndgln, Mmat, &
            DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
            CMC_petsc, CMC_PRECON, IGOT_CMC_PRECON, MASS_MN_PRES, &
            MASS_PIPE, MASS_CVFEM2PIPE, MASS_CVFEM2PIPE_TRUE,  &
            got_free_surf,  MASS_SUF, ndpset, symmetric_P )
            implicit none
            ! form pressure matrix CMC using a colouring approach
            type(multi_dimensions), intent(in) :: Mdims
            type (multi_sparsities), intent(in) :: Mspars
            type(multi_ndgln), intent(in) :: ndgln
            type (multi_matrices), intent(inout) :: Mmat
            INTEGER, intent( in ) :: IGOT_CMC_PRECON
            LOGICAL, PARAMETER :: PIPES_1D=.TRUE.
            LOGICAL, intent( in ) :: got_free_surf, symmetric_P
            INTEGER, DIMENSION( : ), intent( in ) :: ndpset
            REAL, DIMENSION( :, : ), intent( in ) :: DIAG_SCALE_PRES
            REAL, DIMENSION( :, :, : ), intent( in ) :: DIAG_SCALE_PRES_COUP, INV_B
            type(petsc_csr_matrix), intent(inout)::  CMC_petsc
            REAL, DIMENSION( :, :, : ), intent( inout ) :: CMC_PRECON
            REAL, DIMENSION( : ), intent( in ) :: MASS_MN_PRES
            REAL, DIMENSION( : ), intent( in ) :: MASS_SUF
            REAL, DIMENSION( : ), intent( in ) :: MASS_PIPE, MASS_CVFEM2PIPE, MASS_CVFEM2PIPE_TRUE
            ! Local variables
            INTEGER, PARAMETER :: MX_NCOLOR = 1000
            REAL, PARAMETER :: INFINY = 1.0E+10
            LOGICAL :: EXPLICIT_PIPES2
            LOGICAL, DIMENSION( : ), allocatable :: COLOR_LOGICAL
            INTEGER, DIMENSION( : ), allocatable :: COLOR_IN_ROW
            REAL, DIMENSION( :, : ), allocatable :: COLOR_VEC_MANY
            REAL, DIMENSION( :, :, :, : ), allocatable, target :: CDP_MANY
            real, DIMENSION( :, :, :, : ), pointer :: DU_LONG_MANY
            REAL, DIMENSION( :, :, : ), allocatable :: CMC_COLOR_VEC_MANY, CMC_COLOR_VEC2_MANY, &
                CMC_COLOR_VEC_MANY_PHASE, CMC_COLOR_VEC2_MANY_PHASE
            REAL, DIMENSION( : ), allocatable :: CMC_COLOR_VEC_MANY_PHASE_SHORT
            REAL, DIMENSION( :, : ), allocatable :: RINV_B_COLNS_ZEROED, CMC_COLOR_VEC_PRES_MANY_SHORT
            INTEGER :: CV_NOD, CV_JNOD, COUNT, COUNT2, IPHASE, CV_JNOD2
            INTEGER :: MAX_COLOR_IN_ROW, I, ICAN_COLOR, MX_COLOR, NOD_COLOR
            INTEGER :: ierr, i_indx, j_indx, IPRES, JPRES
            REAL :: RSUM, RSUM_SUF
            IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON = 0.0
            EXPLICIT_PIPES2 = .true.
            MAX_COLOR_IN_ROW = 0
            DO CV_NOD = 1, Mdims%cv_nonods
                MAX_COLOR_IN_ROW = MAX( MAX_COLOR_IN_ROW, Mspars%CMC%fin( CV_NOD + 1 ) - Mspars%CMC%fin( CV_NOD ) )
            END DO

            !Calculate coloring if it is not in memory already
            IF ( .not.Mmat%Stored ) THEN
                allocate(Mmat%ICOLOR(Mdims%cv_nonods)); Mmat%ICOLOR = 0
                allocate(COLOR_LOGICAL(Mdims%cv_nonods))
                Mmat%NCOLOR = 0

                COLOR_LOGICAL = .FALSE.
                Loop_CVNOD7: DO CV_NOD = 1, Mdims%cv_nonods
                    ! Color this node CV_NOD
                    MX_COLOR=0
                    ! use a distance-2 colouring...
                    Loop_Row7: DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                        CV_JNOD = Mspars%CMC%col( COUNT )
                        Loop_Row2: DO COUNT2 = Mspars%CMC%fin( CV_JNOD ), Mspars%CMC%fin( CV_JNOD + 1 ) - 1
                            CV_JNOD2 = Mspars%CMC%col( COUNT2 )
                            IF(CV_NOD.NE.CV_JNOD2) THEN
                                NOD_COLOR=Mmat%ICOLOR(CV_JNOD2)
                                IF(NOD_COLOR.NE.0) THEN
                                    COLOR_LOGICAL( NOD_COLOR ) = .TRUE.
                                    MX_COLOR=MAX(MX_COLOR, NOD_COLOR )
                                ENDIF
                            ENDIF
                        END DO Loop_Row2
                    END DO Loop_Row7
                    ! Find the node colour to use...
                    ICAN_COLOR = 0
                    DO I=1,MX_COLOR
                        IF(.NOT.COLOR_LOGICAL(I)) THEN ! This is a colour we can use...
                            ICAN_COLOR = I
                            EXIT
                        ENDIF
                    END DO
                    IF(ICAN_COLOR == 0) THEN
                        ICAN_COLOR =MX_COLOR + 1
                    ENDIF
                    Mmat%ICOLOR(CV_NOD)=ICAN_COLOR
                    Mmat%NCOLOR=MAX(Mmat%NCOLOR,ICAN_COLOR)
                    ! Reset the colors for the next node...
                    COLOR_LOGICAL(1:max(1,MX_COLOR))=.FALSE.
                END DO Loop_CVNOD7
                IF ( Mmat%NCOLOR > MX_NCOLOR ) THEN
                    EWRITE(-1, *) 'NOT ENOUGH COLOURS STOPPING - NEED TO MAKE MX_COLOR BIGGER'
                    STOP 281
                END IF
                deallocate(COLOR_LOGICAL)
            END IF ! SAVED_CMC_COLOR



            ALLOCATE( COLOR_VEC_MANY( Mmat%NCOLOR, Mdims%cv_nonods ) )
            COLOR_VEC_MANY = 0.0
            Loop_CVNOD: DO CV_NOD = 1, Mdims%cv_nonods
                COLOR_VEC_MANY( Mmat%ICOLOR( CV_NOD ), CV_NOD ) = 1.0
            END DO Loop_CVNOD
            ! we use the same colouring for each pressure variable when Mdims%npres>1.
            ALLOCATE( CDP_MANY( Mmat%NCOLOR, Mdims%ndim, Mdims%nphase, Mdims%u_nonods ) )
            CALL C_MULT_MANY( CDP_MANY, COLOR_VEC_MANY, Mdims%cv_nonods, Mdims%u_nonods, Mdims%ndim, Mdims%nphase, Mmat%NCOLOR, &
                Mmat%C, Mspars%C%ncol, Mspars%C%fin, Mspars%C%col )

            CALL PHA_BLOCK_MAT_VEC_MANY_REUSING( Mmat%PIVIT_MAT, CDP_MANY, Mdims%ndim, Mdims%nphase, Mmat%NCOLOR, &
                Mdims%totele, Mdims%u_nloc, ndgln%u )
            DU_LONG_MANY(1:Mmat%NCOLOR, 1:Mdims%ndim, 1:Mdims%nphase, 1:Mdims%u_nonods) => CDP_MANY

                     ! NB. P_RHS = Mmat%CT * U + CV_RHS
                     ! DU_LONG = CDP
            IF ( Mdims%npres > 1 .AND. .NOT.EXPLICIT_PIPES2 ) THEN
                ALLOCATE( CMC_COLOR_VEC_MANY_PHASE( Mmat%NCOLOR, Mdims%nphase, Mdims%cv_nonods ) )
                DO IPHASE = 1, Mdims%nphase
                    CALL CT_MULT_MANY( CMC_COLOR_VEC_MANY_PHASE(:,IPHASE,:), &
                        DU_LONG_MANY(:,:,IPHASE:IPHASE,:), &
                        Mdims%cv_nonods, Mdims%u_nonods, Mdims%ndim, 1, Mmat%NCOLOR, &
                        Mmat%CT(:,IPHASE:IPHASE,:), Mspars%CT%ncol, Mspars%CT%fin, Mspars%CT%col )
                END DO
                if (IGOT_CMC_PRECON==0) deallocate(CDP_MANY)!deallocate memory, not pointer
                ALLOCATE( RINV_B_COLNS_ZEROED(Mdims%nphase,Mdims%nphase) )
                ALLOCATE( CMC_COLOR_VEC_MANY_PHASE_SHORT(Mdims%nphase) )
                ALLOCATE( CMC_COLOR_VEC_PRES_MANY_SHORT(Mmat%NCOLOR,Mdims%npres) )
                DO JPRES=1,Mdims%npres
                    DO CV_NOD = 1, Mdims%cv_nonods
                        RINV_B_COLNS_ZEROED(:,:) = 0.0
                        RINV_B_COLNS_ZEROED(:,  1+(JPRES-1)*Mdims%n_in_pres:JPRES*Mdims%n_in_pres ) = INV_B(:,  1+(JPRES-1)*Mdims%n_in_pres:JPRES*Mdims%n_in_pres,  CV_NOD)
                        DO I = 1, Mmat%NCOLOR
                            CMC_COLOR_VEC_MANY_PHASE_SHORT(:) = MATMUL( RINV_B_COLNS_ZEROED(:,:), CMC_COLOR_VEC_MANY_PHASE(I,:,CV_NOD) )
                            DO IPRES = 1, Mdims%npres
                                CMC_COLOR_VEC_PRES_MANY_SHORT(I,IPRES) = SUM( CMC_COLOR_VEC_MANY_PHASE_SHORT(1+(IPRES-1)*Mdims%n_in_pres:IPRES*Mdims%n_in_pres) )
                            END DO
                        END DO
                        DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                            CV_JNOD = Mspars%CMC%col( COUNT )
                            DO IPRES = 1, Mdims%npres
                                call addto( CMC_petsc, blocki = IPRES, blockj = JPRES, i = cv_nod, j = CV_JNOD, &
                                    val = dot_product(CMC_COLOR_VEC_PRES_MANY_SHORT( :, IPRES ), COLOR_VEC_MANY( :, CV_JNOD ) ))
                            END DO
                        END DO
                    END DO
                END DO ! ENDOF DO JPRES=1,Mdims%npres
                deallocate(CMC_COLOR_VEC_MANY_PHASE)
                IF ( IGOT_CMC_PRECON /= 0 ) THEN
                    ALLOCATE( CMC_COLOR_VEC2_MANY_PHASE( Mmat%NCOLOR, Mdims%nphase, Mdims%cv_nonods ) )
                    DO IPHASE = 1, Mdims%nphase
                        CALL CT_MULT_WITH_C_MANY( CMC_COLOR_VEC2_MANY_PHASE(:,IPHASE,:), &
                            DU_LONG_MANY(:,:,IPHASE:IPHASE,:), &
                            Mdims%u_nonods, Mdims%ndim, 1,&
                            Mmat%C(:,IPHASE:IPHASE,:), Mspars%C%fin, Mspars%C%col )
                    END DO
                    deallocate(CDP_MANY)!deallocate memory, not pointer
                    DO JPRES=1,Mdims%npres
                        DO CV_NOD = 1, Mdims%cv_nonods
                            RINV_B_COLNS_ZEROED(:,:) = 0.0
                            RINV_B_COLNS_ZEROED(:,  1+(JPRES-1)*Mdims%n_in_pres:JPRES*Mdims%n_in_pres ) = INV_B(:,  1+(JPRES-1)*Mdims%n_in_pres:JPRES*Mdims%n_in_pres,  CV_NOD)
                            DO I = 1, Mmat%NCOLOR
                                ! The following uses CMC_COLOR_VEC2_MANY_PHASE:
                                CMC_COLOR_VEC_MANY_PHASE_SHORT(:) = MATMUL( RINV_B_COLNS_ZEROED(:,:), CMC_COLOR_VEC2_MANY_PHASE(I,:,CV_NOD) )
                                DO IPRES = 1, Mdims%npres
                                    CMC_COLOR_VEC_PRES_MANY_SHORT(I,IPRES) = SUM( CMC_COLOR_VEC_MANY_PHASE_SHORT(1+(IPRES-1)*Mdims%n_in_pres:IPRES*Mdims%n_in_pres) )
                                END DO
                            END DO
                            DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                                CV_JNOD = Mspars%CMC%col( COUNT )
                                DO IPRES = 1, Mdims%npres
                                    CMC_PRECON( IPRES, JPRES, COUNT ) = CMC_PRECON( IPRES, JPRES, COUNT ) + &
                                        sum( CMC_COLOR_VEC_PRES_MANY_SHORT( :, IPRES )* COLOR_VEC_MANY( :, CV_JNOD ) )
                                END DO
                            END DO
                        END DO
                    END DO ! ENDOF DO JPRES=1,Mdims%npres
                    deallocate(CMC_COLOR_VEC2_MANY_PHASE)
                END IF
                ! Re-set CMC_COLOR_VEC_MANY & CMC_COLOR_VEC_MANY2 as we have already added the
                ! contributions up to this point into the matrix.

                deallocate(CMC_COLOR_VEC_PRES_MANY_SHORT, RINV_B_COLNS_ZEROED, CMC_COLOR_VEC_MANY_PHASE_SHORT)


            ELSE
                ALLOCATE( CMC_COLOR_VEC_MANY( Mmat%NCOLOR, Mdims%npres, Mdims%cv_nonods ) )
                DO IPRES = 1, Mdims%npres
                    CALL CT_MULT_MANY( CMC_COLOR_VEC_MANY(:,IPRES,:), &
                        DU_LONG_MANY(:,:,1+(IPRES-1)*Mdims%n_in_pres:IPRES*Mdims%n_in_pres,:), &
                        Mdims%cv_nonods, Mdims%u_nonods, Mdims%ndim, Mdims%n_in_pres, Mmat%NCOLOR, &
                        Mmat%CT(:,1+(IPRES-1)*Mdims%n_in_pres:IPRES*Mdims%n_in_pres,:), Mspars%CT%ncol, Mspars%CT%fin, Mspars%CT%col )
                END DO
                IF ( IGOT_CMC_PRECON /= 0 ) THEN
                    ALLOCATE( CMC_COLOR_VEC2_MANY( Mmat%NCOLOR, Mdims%npres, Mdims%cv_nonods ) )
                    DO IPRES = 1, Mdims%npres
                        CALL CT_MULT_WITH_C_MANY( CMC_COLOR_VEC2_MANY(:,IPRES,:), &
                            DU_LONG_MANY(:,:,1+(IPRES-1)*Mdims%n_in_pres:IPRES*Mdims%n_in_pres, :), &
                            Mdims%u_nonods, Mdims%ndim, Mdims%n_in_pres, &
                            Mmat%C(:,1+(IPRES-1)*Mdims%n_in_pres:IPRES*Mdims%n_in_pres,:), Mspars%C%fin, Mspars%C%col )
                    END DO
                END IF
                deallocate(CDP_MANY)!deallocate memory, not pointer
            END IF
            ! Matrix vector involving the mass diagonal term
            if (.not.allocated(CMC_COLOR_VEC_MANY)) ALLOCATE( CMC_COLOR_VEC_MANY( Mmat%NCOLOR, Mdims%npres, Mdims%cv_nonods ) )
            DO CV_NOD = 1, Mdims%cv_nonods
                RSUM = 0.0
                RSUM_SUF = 0.0
                DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                    CV_JNOD = Mspars%CMC%col( COUNT )
                    DO IPRES = 1, Mdims%npres
                        IF((Mdims%npres>1).AND.PIPES_1D) THEN
                            IF(IPRES==1) THEN
                                CMC_COLOR_VEC_MANY( :, IPRES, CV_NOD ) = CMC_COLOR_VEC_MANY( :, IPRES, CV_NOD ) + &
                                    DIAG_SCALE_PRES( IPRES, CV_NOD ) * MASS_MN_PRES( COUNT ) * COLOR_VEC_MANY( :, CV_JNOD )
                            ELSE
                                CMC_COLOR_VEC_MANY( :, IPRES, CV_NOD ) = CMC_COLOR_VEC_MANY( :, IPRES, CV_NOD ) + &
                                    DIAG_SCALE_PRES( IPRES, CV_NOD ) * MASS_CVFEM2PIPE_TRUE( COUNT ) * COLOR_VEC_MANY( :, CV_JNOD )
                            ENDIF
                        ELSE
                            CMC_COLOR_VEC_MANY( :, IPRES, CV_NOD ) = CMC_COLOR_VEC_MANY( :, IPRES, CV_NOD ) + &
                                DIAG_SCALE_PRES( IPRES, CV_NOD ) * MASS_MN_PRES( COUNT ) * COLOR_VEC_MANY( :, CV_JNOD )
                        ENDIF
                        CMC_COLOR_VEC_MANY( :, IPRES, CV_NOD ) = CMC_COLOR_VEC_MANY( :, IPRES, CV_NOD ) + &
                            DIAG_SCALE_PRES( IPRES, CV_NOD ) * MASS_MN_PRES( COUNT ) * COLOR_VEC_MANY( :, CV_JNOD )
                        if ( got_free_surf ) then
                            CMC_COLOR_VEC_MANY( :, IPRES, CV_NOD ) = CMC_COLOR_VEC_MANY( :, IPRES, CV_NOD ) + &
                                MASS_SUF( COUNT ) * COLOR_VEC_MANY( :, CV_JNOD )
                        end if
                    END DO
                    if ( got_free_surf ) then
                        RSUM_SUF = RSUM_SUF + MASS_SUF( COUNT )
                    end if
                    RSUM = RSUM + MASS_MN_PRES( COUNT )
                END DO
                IF ( IGOT_CMC_PRECON /= 0 ) THEN ! Use lumping of MASS_MN_PRES & MASS_SUF...
                    DO IPRES=1,Mdims%npres
                        CMC_COLOR_VEC2_MANY( :, IPRES, CV_NOD ) = CMC_COLOR_VEC2_MANY( :, IPRES, CV_NOD ) &
                            + ( DIAG_SCALE_PRES( IPRES, CV_NOD ) * RSUM + RSUM_SUF ) * COLOR_VEC_MANY( :, CV_NOD )
                    END DO
                END IF
            END DO
            !Put into matrix CMC
            if ( .not.symmetric_P ) then
                ! original method
                DO CV_NOD = 1, Mdims%cv_nonods
                    DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                        CV_JNOD = Mspars%CMC%col( COUNT )
                        DO IPRES = 1, Mdims%npres
                            JPRES = IPRES ! Add contributions to the block diagonal only.
                            call addto( CMC_petsc, blocki = IPRES, blockj = JPRES, i = cv_nod, j = CV_JNOD, &
                                val = dot_product(CMC_COLOR_VEC_MANY( :, IPRES, CV_NOD ), COLOR_VEC_MANY( :, CV_JNOD ) ))
                            IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON( IPRES, JPRES, COUNT ) = CMC_PRECON( IPRES, JPRES, COUNT ) + &
                                sum(CMC_COLOR_VEC2_MANY( :, IPRES, CV_NOD ) * COLOR_VEC_MANY( :, CV_JNOD ))
                        END DO
                    END DO
                END DO
            else
                ! symmetric P matrix
                DO CV_NOD = 1, Mdims%cv_nonods
                    DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                        CV_JNOD = Mspars%CMC%col( COUNT )
                        DO IPRES = 1, Mdims%npres
                            JPRES = IPRES ! Add contributions to the block diagonal only.
                            call addto( CMC_petsc, blocki = IPRES, blockj = JPRES, i = cv_nod, j = CV_JNOD, &
                                val = dot_product(CMC_COLOR_VEC2_MANY( :, IPRES, CV_NOD ), COLOR_VEC_MANY( :, CV_JNOD ) ))!this looks like a bug, CMC_COLOR_VEC2_MANY is not defined for not IGOT_CMC_PRECON
                            IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON( IPRES, JPRES, COUNT ) = CMC_PRECON( IPRES, JPRES, COUNT ) + &
                                sum(CMC_COLOR_VEC2_MANY( :, IPRES, CV_NOD ) * COLOR_VEC_MANY( :, CV_JNOD ))
                        END DO
                    END DO
                END DO
            end if
            deallocate(CMC_COLOR_VEC_MANY, COLOR_VEC_MANY)
            ! the matrix coupling term place immediately into matrix and not through colouring...
            ! Matrix vector involving the mass diagonal term
            IF ( Mdims%npres > 1) THEN
                DO CV_NOD = 1, Mdims%cv_nonods
                    RSUM = 0.0
                    RSUM_SUF = 0.0
                    DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                        CV_JNOD = Mspars%CMC%col( COUNT )
                        DO IPRES = 1, Mdims%npres
                            DO JPRES = 1, Mdims%npres
                                IF(PIPES_1D) THEN ! 1D PIPE MODELLING
                                    call addto( CMC_petsc, blocki = IPRES, blockj = JPRES, i = cv_nod, j = CV_JNOD,&
                                        val = DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_NOD ) * MASS_CVFEM2PIPE( COUNT ))
                                    IF ( IGOT_CMC_PRECON /= 0 ) THEN ! Use lumping of MASS_MN_PRES & MASS_SUF...
                                        CMC_PRECON( IPRES, JPRES, COUNT ) = CMC_PRECON( IPRES, JPRES, COUNT ) &
                                            + sqrt(DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_NOD ))* &
                                            MASS_CVFEM2PIPE( COUNT )*sqrt(DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_JNOD ))
                                    END IF
                                ELSE
                                    call addto( CMC_petsc, blocki = IPRES, blockj = JPRES, i = cv_nod, j = CV_JNOD,&
                                        val = DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_NOD ) * MASS_MN_PRES( COUNT ))
                                    IF ( IGOT_CMC_PRECON /= 0 ) THEN ! Use lumping of MASS_MN_PRES & MASS_SUF...
                                        CMC_PRECON( IPRES, JPRES, COUNT ) = CMC_PRECON( IPRES, JPRES, COUNT ) &
                                            + sqrt(DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_NOD ))* &
                                            MASS_MN_PRES( COUNT )*sqrt(DIAG_SCALE_PRES_COUP( IPRES, JPRES, CV_JNOD ))
                                    END IF
                                ENDIF
                            END DO
                        END DO
                    END DO
                END DO
            END IF ! ENDOF IF(Mdims%npres > 1) THEN
            !If we have a reference node with pressure zero we impose that here.
            DO IPRES = 1, Mdims%npres
                IF ( NDPSET(IPRES) > 0 ) THEN
                    CV_NOD = NDPSET(IPRES)
                    i_indx = CMC_petsc%row_numbering%gnn2unn( cv_nod, ipres )
                    DO COUNT = Mspars%CMC%fin( CV_NOD ), Mspars%CMC%fin( CV_NOD + 1 ) - 1
                        CV_JNOD = Mspars%CMC%col( COUNT )
                        IF ( CV_JNOD /= CV_NOD ) THEN
                            do jpres = 1, Mdims%npres
                                j_indx = CMC_petsc%column_numbering%gnn2unn( cv_jnod, jpres )
                                call MatSetValue(CMC_petsc%M, i_indx, j_indx, 0.0,INSERT_VALUES, ierr) ! not the diagonal
                                !CMC( COUNT ) = 0.0 ! not the diagonal
                                IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON( IPRES, JPRES, COUNT ) = 0.0
                                DO COUNT2 = Mspars%CMC%fin( CV_JNOD ), Mspars%CMC%fin( CV_JNOD + 1 ) - 1
                                    CV_JNOD2 = Mspars%CMC%col( COUNT2 )
                                    IF ( CV_JNOD2 /= CV_NOD ) then
                                        i_indx = CMC_petsc%row_numbering%gnn2unn( cv_jnod, ipres )
                                        j_indx = CMC_petsc%column_numbering%gnn2unn( CV_JNOD2, jpres )
                                        call MatSetValue(CMC_petsc%M, i_indx, j_indx, 0.0, INSERT_VALUES, ierr) ! not the diagonal
                                        IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON( ipres, jpres, COUNT2 ) = 0.0
                                    END IF
                                END DO
                            end do ! do jpres=1,Mdims%npres
                        END IF
                    END DO
                END IF
            END DO
            !Re-assemble
            CMC_petsc%is_assembled=.false.
            call assemble( CMC_petsc )
            IF ( IGOT_CMC_PRECON /= 0 ) deallocate(CMC_COLOR_VEC2_MANY)
            nullify(DU_LONG_MANY)
            RETURN
        END SUBROUTINE COLOR_GET_CMC_PHA_FAST

    END SUBROUTINE COLOR_GET_CMC_PHA





    SUBROUTINE PHA_BLOCK_INV( PIVIT_MAT, Mdims )
        implicit none
        REAL, DIMENSION( : , : , : ), intent( inout ), CONTIGUOUS :: PIVIT_MAT
        type(multi_dimensions), intent(in) :: Mdims
        ! Local variables
        logical, save :: lumped_matrix = .false., check_lumped_matrix = .true.
        real :: aux
        INTEGER :: ELE, iphase, idim, u_iloc, i, u_jloc, k
        REAL, DIMENSION( :,: ), allocatable :: MAT
        REAL, DIMENSION( : ), allocatable :: B


        if(size(PIVIT_MAT,1) == 1) then
            !If it is compacted and diagonal the inverse is straightforward
            PIVIT_MAT = 1./PIVIT_MAT
            return
        end if

        !First time only, check if the Mass matrix is lumped. If it is, the inverse can be done much quicker
        if (check_lumped_matrix) then
            aux = 0.; check_lumped_matrix = .false.
            !Check first element only
            do k = 1, size(PIVIT_MAT,1)
                aux = aux + PIVIT_MAT(k, k, 1)
            end do
            !If all the information is contained in the diagonals, then the mass matrix is lumped
            lumped_matrix = ((sum(abs(PIVIT_MAT(:, :, 1))) - aux)/sum(abs(PIVIT_MAT(:, :, 1))) < 1d-8)
        end if

        if (lumped_matrix) then
            !Only invert the diagonal
             DO ELE = 1, Mdims%TOTELE
                do k = 1, size(PIVIT_MAT,1)
                    PIVIT_MAT(k,k,ELE) = 1./PIVIT_MAT(k,k,ELE)
                end do
             end do
            return
        end if


        if (is_porous_media) then !No coupling between phases nor dimensions, inverse can be done faster
             allocate(mat(Mdims%u_nloc, Mdims%u_nloc))
             DO ELE = 1, Mdims%TOTELE
                k = 1
                !Compress into a mini matrix
                do u_iloc = 1, Mdims%u_nloc
                    do u_jloc = 1, Mdims%u_nloc
                        mat(u_iloc, u_jloc) = PIVIT_MAT( k + (u_iloc-1)*mdims%nphase * mdims%ndim, &
                                 k + (u_jloc-1)*mdims%nphase * mdims%ndim, ele )
                    end do
                end do
                !Invert
                mat = inverse(mat)
                !Populate PIVIT_MAT. Since the matrix is repeated mdims%nphase * mdims%ndim times we don't need to
                !invert it that many times
                do i = 1, mdims%nphase * mdims%ndim
                    do u_iloc = 1, Mdims%u_nloc
                        do u_jloc = 1, Mdims%u_nloc
                            PIVIT_MAT( k + (u_iloc-1)*mdims%nphase * mdims%ndim, &
                                     k + (u_jloc-1)*mdims%nphase * mdims%ndim, ele ) = mat(u_iloc, u_jloc)
                        end do
                    end do
                    k = k + 1
                end do
            END DO
        else
            allocate(MAT( Mdims%u_nloc * Mdims%nphase * Mdims%ndim , Mdims%u_nloc * Mdims%nphase * Mdims%ndim ))
            allocate(B( Mdims%u_nloc * Mdims%nphase * Mdims%ndim ))
            DO ELE = 1, Mdims%TOTELE
                CALL MATINVold( PIVIT_MAT( :, :, ele ), Mdims%u_nloc * Mdims%nphase * Mdims%ndim, MAT, B )
            END DO
            deallocate(b)
        end if
        deallocate(MAT)
        RETURN
    END SUBROUTINE PHA_BLOCK_INV

    SUBROUTINE PHA_BLOCK_MAT_VEC_old( U, BLOCK_MAT, CDP, NDIM, NPHASE, &
        TOTELE, U_NLOC, U_NDGLN )
        implicit none
        ! U = BLOCK_MAT * CDP
        INTEGER, intent( in )  :: NDIM, NPHASE, TOTELE, U_NLOC
        INTEGER, DIMENSION( : ), intent( in ), target ::  U_NDGLN
        REAL, DIMENSION( : ), intent( inout ) :: U
        REAL, DIMENSION( :, :, : ), intent( in ), target :: BLOCK_MAT
        REAL, DIMENSION( :, :, : ), intent( in ) :: CDP
        ! Local
        INTEGER :: ELE, N, U_ILOC
        INTEGER, DIMENSION(:), pointer :: U_NOD
        REAL, DIMENSION( NDIM, NPHASE, U_NLOC ) :: LU
        REAL, DIMENSION( NDIM * NPHASE * U_NLOC ) :: LCDP

        interface
            subroutine dgemv(T,M,N,alpha,MAT,NMAX,X,Xinc,beta,Y,yinc)
                implicit none
                character(len=1) :: T
                integer :: m,n,nmax,xinc,yinc
                real ::  alpha, beta
                real, dimension(nmax,n) :: MAT
                real, dimension(N) :: X
                real, dimension(M) :: Y
            end subroutine dgemv
        end interface
           

        N=U_NLOC * NDIM * NPHASE

        if (size(BLOCK_MAT,1) == 1) then!BLOCK_MAT is diagonal and compacted
            DO ELE = 1, TOTELE
                U_NOD => U_NDGLN( ( ELE - 1 ) * U_NLOC +1 : ELE * U_NLOC )
                LU = CDP( :, :, U_NOD ) * BLOCK_MAT(1,1,ELE)
                DO U_ILOC = 1, U_NLOC
                    U( 1+(U_NOD(U_ILOC)-1)*NDIM*NPHASE : U_NOD(U_ILOC)*NDIM*NPHASE ) = [LU(:,:,U_ILOC)]
                END DO
            END DO
        else
            Loop_Elements: DO ELE = 1, TOTELE

                U_NOD => U_NDGLN( ( ELE - 1 ) * U_NLOC +1 : ELE * U_NLOC )

                LCDP = RESHAPE( CDP( :, :, U_NOD ) , (/ N /) )
                CALL DGEMV( 'N', N, N, 1.0d0, BLOCK_MAT( : , : , ELE ), N, LCDP, 1, 0.0d0, LU, 1 )
                DO U_ILOC = 1, U_NLOC
                    U( 1+(U_NOD(U_ILOC)-1)*NDIM*NPHASE : U_NOD(U_ILOC)*NDIM*NPHASE ) = [LU(:,:,U_ILOC)]
                END DO

            END DO Loop_Elements
        end if
        RETURN


    END SUBROUTINE PHA_BLOCK_MAT_VEC_old


    SUBROUTINE PHA_BLOCK_MAT_VEC( U, BLOCK_MAT, CDP, U_NONODS, NDIM, NPHASE, &
        TOTELE, U_NLOC, U_NDGLN )
        implicit none
        ! U = BLOCK_MAT * CDP
        INTEGER, intent( in )  :: U_NONODS, NDIM, NPHASE, TOTELE, U_NLOC
        INTEGER, DIMENSION( : ), intent( in ), target ::  U_NDGLN
        REAL, DIMENSION( : ), intent( inout ) :: U
        REAL, DIMENSION( :, :, : ), intent( in ), target :: BLOCK_MAT
        REAL, DIMENSION( ndim, nphase, u_nonods ), intent( inout ) :: CDP
        ! Local
        INTEGER :: ELE, I, JDIM, JPHASE, J, JJ

        integer, dimension(:), pointer :: U_NOD

        real, dimension(U_NLOC*NDIM*NPHASE) :: lcdp, lu
        integer, dimension(U_NLOC*NDIM*NPHASE) :: u_nodi
        integer :: N
      
        interface
            subroutine dgemv(T,M,N,alpha,MAT,NMAX,X,Xinc,beta,Y,yinc)
                implicit none
                character(len=1) :: T
                integer :: m,n,nmax,xinc,yinc
                real ::  alpha, beta
                real, dimension(nmax,n) :: MAT
                real, dimension(N) :: X
                real, dimension(M) :: Y
            end subroutine dgemv
        end interface
           

        N=U_NLOC * NDIM * NPHASE

        if (size(BLOCK_MAT,1) == 1) then!BLOCK_MAT is diagonal and compacted
            DO ELE = 1, TOTELE
                U_NOD => U_NDGLN(( ELE - 1 ) * U_NLOC +1: ELE * U_NLOC)
                DO JPHASE = 1, NPHASE
                    DO JDIM = 1, NDIM
                        J = JDIM + (JPHASE-1)*NDIM
                        JJ = ( JDIM - 1 ) * U_NLOC + ( JPHASE - 1 ) * NDIM * U_NLOC
                        lcdp([(J+(i-1)*ndim*nphase,i=1,u_NLOC)]) = CDP( JDIM, JPHASE, U_NOD )
                        U_NODI([(J+(i-1)*ndim*nphase,i=1,u_NLOC)]) = U_NOD+(J-1)*U_NONODS
                    end do
                end do
                U( U_NODI ) = lcdp * BLOCK_MAT( 1 , 1 , ele )
            END DO
        else
            Loop_Elements: DO ELE = 1, TOTELE

                U_NOD => U_NDGLN(( ELE - 1 ) * U_NLOC +1: ELE * U_NLOC)


                Loop_PhasesJ: DO JPHASE = 1, NPHASE
                    Loop_DimensionsJ: DO JDIM = 1, NDIM

                        J = JDIM + (JPHASE-1)*NDIM
                        JJ = ( JDIM - 1 ) * U_NLOC + ( JPHASE - 1 ) * NDIM * U_NLOC

                        lcdp([(J+(i-1)*ndim*nphase,i=1,u_NLOC)]) = CDP( JDIM, JPHASE, U_NOD )
                        U_NODI([(J+(i-1)*ndim*nphase,i=1,u_NLOC)]) = U_NOD+(J-1)*U_NONODS
                    end do Loop_DimensionsJ
                end do Loop_PhasesJ

                call dgemv( 'N', N, N, 1.0d0, BLOCK_MAT( : , : , ele ), N, LCDP, 1, 0.0d0, LU, 1 )
                U( U_NODI ) = LU

            END DO Loop_Elements
        end if
        RETURN


    END SUBROUTINE PHA_BLOCK_MAT_VEC

    SUBROUTINE PHA_BLOCK_MAT_VEC2( U, BLOCK_MAT, CDP, NDIM, NPHASE, &
        TOTELE, U_NLOC, U_NDGLN )
        implicit none
        ! U = BLOCK_MAT * CDP
        INTEGER, intent( in )  :: NDIM, NPHASE, TOTELE, U_NLOC
        INTEGER, DIMENSION( : ), intent( in ), target ::  U_NDGLN
        REAL, DIMENSION( :, :, : ), intent( inout ) :: U
        REAL, DIMENSION( :, :, : ), intent( in ), target :: BLOCK_MAT
        REAL, DIMENSION( :, :, : ), intent( in ) :: CDP
        ! Local
        INTEGER :: ELE, N
        INTEGER, DIMENSION( U_NLOC ) :: U_NOD
        REAL, DIMENSION( NDIM * NPHASE * U_NLOC ) :: LCDP, LU
      
        interface
            subroutine dgemv(T,M,N,alpha,MAT,NMAX,X,Xinc,beta,Y,yinc)
                implicit none
                character(len=1) :: T
                integer :: m,n,nmax,xinc,yinc
                real ::  alpha, beta
                real, dimension(nmax,n) :: MAT
                real, dimension(N) :: X
                real, dimension(M) :: Y
            end subroutine dgemv
        end interface
           
        N = U_NLOC * NDIM * NPHASE

        if (size(BLOCK_MAT,1) == 1) then
            DO ELE = 1, TOTELE
                U_NOD = U_NDGLN( ( ELE - 1 ) * U_NLOC +1 : ELE * U_NLOC )
                U( :, :, U_NOD ) = CDP( :, :, U_NOD ) * BLOCK_MAT( 1 , 1 , ELE )
            END DO
        else
            Loop_Elements: DO ELE = 1, TOTELE
                U_NOD = U_NDGLN( ( ELE - 1 ) * U_NLOC +1 : ELE * U_NLOC )

                LCDP = RESHAPE( CDP( :, :, U_NOD ) , (/ N /) )
                CALL DGEMV( 'N', N, N, 1.0d0, BLOCK_MAT( : , : , ELE ), N, LCDP, 1, 0.0d0, LU, 1 )
                U( :, :, U_NOD ) = RESHAPE( LU, (/ NDIM, NPHASE, U_NLOC/) )

            END DO Loop_Elements
        end if

        RETURN

    END SUBROUTINE PHA_BLOCK_MAT_VEC2

    SUBROUTINE PHA_BLOCK_MAT_VEC_MANY2( U, BLOCK_MAT, CDP, U_NONODS, NDIM, NPHASE, NBLOCK, &
        TOTELE, U_NLOC, U_NDGLN )
        implicit none
        ! U = BLOCK_MAT * CDP
        INTEGER, intent( in )  :: U_NONODS, NDIM, NPHASE, TOTELE, U_NLOC, NBLOCK
        INTEGER, DIMENSION( : ), intent( in ), target ::  U_NDGLN
        REAL, DIMENSION( :, : ), intent( inout ) :: U
        REAL, DIMENSION( :, :, : ), intent( in ), target :: BLOCK_MAT
        REAL, DIMENSION( :, : ), intent( in ) :: CDP
        ! Local
        INTEGER :: ELE, I, JDIM, JPHASE, J, JJ, IVEC

        integer, dimension(:), pointer :: U_NOD

        real, dimension(U_NLOC*NDIM*NPHASE) :: lcdp, lu
        integer, dimension(U_NLOC*NDIM*NPHASE) :: u_nodi
        integer :: N
      
        interface
            subroutine dgemv(T,M,N,alpha,MAT,NMAX,X,Xinc,beta,Y,yinc)
                implicit none
                character(len=1) :: T
                integer :: m,n,nmax,xinc,yinc
                real ::  alpha, beta
                real, dimension(nmax,n) :: MAT
                real, dimension(N) :: X
                real, dimension(M) :: Y
            end subroutine dgemv
        end interface

        N = U_NLOC * NDIM * NPHASE

        Loop_Elements: DO ELE = 1, TOTELE

            U_NOD => U_NDGLN( ( ELE - 1 ) * U_NLOC +1 : ELE * U_NLOC )

            DO IVEC = 1, NBLOCK

                Loop_PhasesJ: DO JPHASE = 1, NPHASE
                    Loop_DimensionsJ: DO JDIM = 1, NDIM

                        J = JDIM + ( JPHASE - 1 ) * NDIM
                        JJ = ( JDIM - 1 ) * U_NLOC + ( JPHASE - 1 ) * NDIM * U_NLOC

                        lcdp([(J+(i-1)*ndim*nphase, i = 1, u_NLOC)]) = CDP(IVEC, U_NOD+(J-1)*U_NONODS)
                        U_NODI([(J+(i-1)*ndim*nphase, i = 1, u_NLOC)]) = U_NOD+(J-1)*U_NONODS
                    end do Loop_DimensionsJ
                end do Loop_PhasesJ

                call dgemv('N', N, N, 1.0d0, BLOCK_MAT( : , : , ELE ), N, LCDP, 1, 0.0d0, LU, 1 )
                U( IVEC, U_NODI ) = LU

            END DO

        END DO Loop_Elements

        RETURN


    END SUBROUTINE PHA_BLOCK_MAT_VEC_MANY2


 

    SUBROUTINE PHA_BLOCK_MAT_VEC_MANY( U, BLOCK_MAT, CDP, NDIM, NPHASE, NBLOCK, &
        TOTELE, U_NLOC, U_NDGLN )
        implicit none
        ! U = BLOCK_MAT * CDP
        INTEGER, intent( in )  :: NDIM, NPHASE, TOTELE, U_NLOC, NBLOCK
        INTEGER, DIMENSION( : ), intent( in ) ::  U_NDGLN
        REAL, DIMENSION( :, :, :, : ), intent( inout ), contiguous, target :: U
        REAL, DIMENSION( :, : , : ), intent( in ), contiguous :: BLOCK_MAT
        REAL, DIMENSION( :, : , :, : ), intent( in ), contiguous, target :: CDP
        ! Local
        INTEGER :: ELE

        real, dimension(:,:,:,:), pointer, contiguous :: lcdp, lu
        integer :: N
       
        interface
            subroutine dgemm(TA,TB,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC)
                implicit none
                character(len=1) :: TA,TB
                integer :: m,n,k,lda,ldb,ldc
                real ::  alpha, beta
                real, dimension(lda,*) :: A
                real, dimension(ldb,*) :: B
                real, dimension(ldc,*) :: C
            end subroutine dgemm
        end interface

        U = 0.0
        N=ndim*nphase*u_nloc

        Loop_Elements: DO ELE = 1, TOTELE

            lu=>u(:,:,:,U_NDGLN( (ELE-1)*U_NLOC + 1):U_NDGLN(ELE * U_NLOC))
            lcdp=>cdp(:,:,:,U_NDGLN( (ELE-1)*U_NLOC + 1):U_NDGLN(ELE * U_NLOC))
            call dgemm('N','T',NBLOCK,N,N,1.0,LCDP,NBLOCK,Block_mat(:,:,ele),N,1.0,LU,NBLOCK)
       
        END DO Loop_Elements

        RETURN

    END SUBROUTINE PHA_BLOCK_MAT_VEC_MANY

    SUBROUTINE PHA_BLOCK_MAT_VEC_MANY_REUSING( BLOCK_MAT, CDP, NDIM, NPHASE, NBLOCK, &
        TOTELE, U_NLOC, U_NDGLN )
        implicit none
        ! U = BLOCK_MAT * CDP
        !The difference with PHA_BLOCK_MAT_VEC_MANY is that the input CDP is overwritten with the output
        INTEGER, intent( in )  :: NDIM, NPHASE, TOTELE, U_NLOC, NBLOCK
        INTEGER, DIMENSION( : ), intent( in ) ::  U_NDGLN
        REAL, DIMENSION( :, : , : ), intent( in ), contiguous :: BLOCK_MAT
        REAL, DIMENSION( :, : , :, : ), intent( inout ), contiguous, target :: CDP
        ! Local
        INTEGER :: ELE

        real, dimension(:,:,:,:), pointer, contiguous :: lcdp, lu
        integer :: N


        interface
            subroutine dgemm(TA,TB,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC)
                implicit none
                character(len=1) :: TA,TB
                integer :: m,n,k,lda,ldb,ldc
                real ::  alpha, beta
                real, dimension(lda,*) :: A
                real, dimension(ldb,*) :: B
                real, dimension(ldc,*) :: C
            end subroutine dgemm
        end interface


        N=ndim*nphase*u_nloc
        if (size(BLOCK_MAT,1) == 1) then
            DO ELE = 1, TOTELE
                lcdp=>cdp(:,:,:,U_NDGLN( (ELE-1)*U_NLOC + 1):U_NDGLN(ELE * U_NLOC))
                lcdp = lcdp * Block_mat(1,1,ele)
            END DO
        else
            allocate(lu(size(CDP,1), size(CDP,2), size(CDP,3),U_NLOC))
            Loop_Elements: DO ELE = 1, TOTELE
                lu = 0.!Initialize memory
                lcdp=>cdp(:,:,:,U_NDGLN( (ELE-1)*U_NLOC + 1):U_NDGLN(ELE * U_NLOC))
                call dgemm('N','T',NBLOCK,N,N,1.0,LCDP,NBLOCK,Block_mat(:,:,ele),N,1.0,LU,NBLOCK)
                lcdp = lu!copy to original array
            END DO Loop_Elements
            deallocate(lu)
        end if
        RETURN

    END SUBROUTINE PHA_BLOCK_MAT_VEC_MANY_REUSING




    SUBROUTINE CT_MULT( CV_RHS, U, V, W, CV_NONODS, U_NONODS, NDIM, NPHASE, &
        CT, NCOLCT, FINDCT, COLCT )
        ! CV_RHS=CT*U
        implicit none
        INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLCT
        REAL, DIMENSION( CV_NONODS ), intent( inout) :: CV_RHS
        REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W
        INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
        INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
        REAL, DIMENSION( NDIM, NPHASE, NCOLCT ), intent( in ) :: CT

        ! Local variables
        INTEGER :: CV_INOD, COUNT, U_JNOD, IPHASE, J
        logical :: is_3d

        is_3d = (NDIM >= 3)
        CV_RHS = 0.0

        DO CV_INOD = 1, CV_NONODS
            DO COUNT = FINDCT( CV_INOD ), FINDCT( CV_INOD + 1 ) - 1, 1
                U_JNOD = COLCT( COUNT )
                DO IPHASE = 1, NPHASE
                    J = U_JNOD + ( IPHASE - 1 ) * U_NONODS
                    CV_RHS( CV_INOD ) = CV_RHS( CV_INOD ) + CT( 1, IPHASE, COUNT ) * U( J ) + CT( 2, IPHASE, COUNT ) * V( J )
                    IF( is_3d ) CV_RHS( CV_INOD ) = CV_RHS( CV_INOD ) + CT( 3, IPHASE, COUNT ) * W( J )
                END DO
            END DO
        END DO

        RETURN

    END SUBROUTINE CT_MULT


    SUBROUTINE CT_MULT2( CV_RHS, U, CV_NONODS, U_NONODS, NDIM, NPHASE, &
        CT, NCOLCT, FINDCT, COLCT )
        ! CV_RHS=CT*U
        implicit none
        INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLCT
        REAL, DIMENSION( CV_NONODS ), intent( inout ) :: CV_RHS
        REAL, DIMENSION( NDIM, NPHASE, U_NONODS ), intent( in ) :: U
        INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
        INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
        REAL, DIMENSION( NDIM, NPHASE, NCOLCT ), intent( in ) :: CT

        ! Local variables
        INTEGER :: CV_INOD, COUNT, U_JNOD
        CV_RHS = 0.0
        DO CV_INOD = 1, CV_NONODS
            DO COUNT = FINDCT( CV_INOD ), FINDCT( CV_INOD + 1 ) - 1
                U_JNOD = COLCT( COUNT )
                CV_RHS( CV_INOD ) = CV_RHS( CV_INOD ) + SUM( CT( :, :, COUNT ) * U( :, :, U_JNOD ) )
            END DO
        END DO
        RETURN

    END SUBROUTINE CT_MULT2

    SUBROUTINE CT_MULT_MANY( CV_RHS, U, CV_NONODS, U_NONODS, NDIM, NPHASE, NBLOCK, &
        CT, NCOLCT, FINDCT, COLCT )
        ! CV_RHS = CT * U
        implicit none
        INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLCT, NBLOCK
        REAL, DIMENSION( NBLOCK, CV_NONODS ), intent( inout) :: CV_RHS!NBLOCK, CV_NONODS
        REAL, DIMENSION( NBLOCK, NDIM, NPHASE, U_NONODS ), intent( in ) :: U!NBLOCK, NDIM, NPHASE, U_NONODS
        INTEGER, DIMENSION( : ), intent( in ) :: FINDCT!CV_NONODS + 1
        INTEGER, DIMENSION( : ), intent( in ) :: COLCT!NCOLCT
        REAL, DIMENSION( :, :, : ), intent( in ) :: CT!NDIM, NPHASE, NCOLCT

        ! Local variables
        INTEGER :: CV_INOD, COUNT, U_JNOD

        interface
            subroutine dgemv(T,M,N,alpha,MAT,NMAX,X,Xinc,beta,Y,yinc)
                implicit none
                character(len=1) :: T
                integer :: m,n,nmax,xinc,yinc
                real ::  alpha, beta
                real, dimension(nmax,n) :: MAT
                real, dimension(N) :: X
                real, dimension(M) :: Y
            end subroutine dgemv
        end interface

        CV_RHS = 0.0

        DO CV_INOD = 1, CV_NONODS
            DO COUNT = FINDCT( CV_INOD ), FINDCT( CV_INOD + 1 ) - 1
                U_JNOD = COLCT( COUNT )
                call dgemv('N',NBLOCK,NPHASE*NDIM,1.0,U(:,:,:,U_JNOD),NBLOCK,CT(:,:,COUNT),1,1.0,CV_RHS(:,CV_INOD),1)
            END DO
        END DO

        RETURN

    END SUBROUTINE CT_MULT_MANY

    SUBROUTINE C_MULT_MANY( CDP, DP, CV_NONODS, U_NONODS, NDIM, NPHASE, NBLOCK, &
        C, NCOLC, FINDC, COLC )
        implicit none
        ! CDP=C*DP
        INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLC, NBLOCK
        REAL, DIMENSION( :, :, :, : ), intent( inout ) :: CDP!NBLOCK, NDIM, NPHASE, U_NONODS
        REAL, DIMENSION( :, : ), intent( in )  :: DP!NBLOCK, CV_NONODS
        REAL, DIMENSION( :, :, : ), intent( in ) :: C!NDIM, NPHASE, NCOLC
        INTEGER, DIMENSION( : ), intent( in ) ::FINDC!U_NONODS + 1
        INTEGER, DIMENSION( : ), intent( in ) :: COLC!NCOLC
        ! Local variables
        INTEGER :: U_INOD, COUNT, P_JNOD, IPHASE, IDIM

        CDP = 0.0

        DO U_INOD = 1, U_NONODS

            DO COUNT = FINDC( U_INOD ), FINDC( U_INOD + 1 ) - 1

                P_JNOD = COLC( COUNT )
                FORALL (IPHASE = 1:NPHASE, IDIM = 1:NDIM)
                    CDP( :, IDIM, IPHASE, U_INOD ) = CDP( :, IDIM, IPHASE, U_INOD ) + C( IDIM, IPHASE, COUNT ) * DP( :, P_JNOD )
                END FORALL

            END DO

        END DO

        RETURN

    END SUBROUTINE C_MULT_MANY

    SUBROUTINE C_MULT2( CDP, DP, CV_NONODS, U_NONODS, NDIM, NPHASE, &
        C, NCOLC, FINDC, COLC )
        implicit none
        ! CDP=C*DP
        INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLC
        REAL, DIMENSION( NDIM, NPHASE, U_NONODS ), intent( inout ) :: CDP
        REAL, DIMENSION( CV_NONODS ), intent( in )  :: DP
        REAL, DIMENSION( NDIM, NPHASE, NCOLC ), intent( in ) :: C
        INTEGER, DIMENSION( U_NONODS + 1 ), intent( in ) ::FINDC
        INTEGER, DIMENSION( NCOLC ), intent( in ) :: COLC
        ! Local variables
        INTEGER :: U_INOD, COUNT, P_JNOD

        CDP = 0.0

        DO U_INOD = 1, U_NONODS
            DO COUNT = FINDC( U_INOD ), FINDC( U_INOD + 1 ) - 1
                P_JNOD = COLC( COUNT )
                CDP( :, :, U_INOD ) = CDP( :, :, U_INOD ) + C( :, :, COUNT ) * DP( P_JNOD )
            END DO
        END DO

        RETURN

    END SUBROUTINE C_MULT2

    SUBROUTINE CT_MULT_WITH_C( DP, U_LONG, U_NONODS, NDIM, NPHASE, &
        C, NCOLC, FINDC, COLC )
        implicit none
        ! DP = (C)^T U_LONG
        INTEGER, intent( in ) :: U_NONODS, NDIM, NPHASE, NCOLC
        REAL, DIMENSION( : ), intent( in ) :: U_LONG
        REAL, DIMENSION( : ), intent( inout )  :: DP
        REAL, DIMENSION( :, :, : ), intent( in ) :: C
        INTEGER, DIMENSION( : ), intent( in ) ::FINDC
        INTEGER, DIMENSION( : ), intent( in ) :: COLC
        ! Local variables
        INTEGER :: U_INOD, COUNT, P_JNOD, IPHASE, I1, IDIM, COUNT_DIM_PHA

        DP = 0.0

        Loop_VelNodes: DO U_INOD = 1, U_NONODS

            Loop_Crow: DO COUNT = FINDC( U_INOD ), FINDC( U_INOD + 1 ) - 1, 1
                P_JNOD = COLC( COUNT )

                Loop_Phase: DO IPHASE = 1, NPHASE
                    Loop_Dim: DO IDIM = 1, NDIM
                        COUNT_DIM_PHA = COUNT + NCOLC*(IDIM-1) + NCOLC*NDIM*(IPHASE-1)
                        I1 = U_INOD + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS
                        DP( P_JNOD ) = DP( P_JNOD ) + C( IDIM, IPHASE, COUNT ) * U_LONG( I1 )
                    END DO Loop_Dim
                END DO Loop_Phase

            END DO Loop_Crow

        END DO Loop_VelNodes

        RETURN

    END SUBROUTINE CT_MULT_WITH_C

    SUBROUTINE CT_MULT_WITH_C3( DP, U_ALL, U_NONODS, NDIM, NPHASE, &
        C, NCOLC, FINDC, COLC )
        implicit none
        ! DP = (C)^T U_ALL
        INTEGER, intent( in ) :: U_NONODS, NDIM, NPHASE, NCOLC
        REAL, DIMENSION( :, :, : ), intent( in ) :: U_ALL
        REAL, DIMENSION( : ), intent( inout )  :: DP
        REAL, DIMENSION( :, :, : ), intent( in ) :: C
        INTEGER, DIMENSION( : ), intent( in ) ::FINDC
        INTEGER, DIMENSION( : ), intent( in ) :: COLC
        ! Local variables
        INTEGER :: U_INOD, COUNT, P_JNOD, IPHASE, IDIM, COUNT_DIM_PHA

        DP = 0.0

        Loop_VelNodes: DO U_INOD = 1, U_NONODS

            Loop_Crow: DO COUNT = FINDC( U_INOD ), FINDC( U_INOD + 1 ) - 1, 1
                P_JNOD = COLC( COUNT )

                Loop_Phase: DO IPHASE = 1, NPHASE
                    Loop_Dim: DO IDIM = 1, NDIM
                        COUNT_DIM_PHA = COUNT + NCOLC*(IDIM-1) + NCOLC*NDIM*(IPHASE-1)
                        DP( P_JNOD ) = DP( P_JNOD ) + C( IDIM, IPHASE, COUNT ) * U_ALL( IDIM, IPHASE, U_INOD )
                    END DO Loop_Dim
                END DO Loop_Phase

            END DO Loop_Crow

        END DO Loop_VelNodes

        RETURN

    END SUBROUTINE CT_MULT_WITH_C3


    SUBROUTINE CT_MULT_WITH_C_MANY( DP, U_LONG, U_NONODS, NDIM, NPHASE, &
        C, FINDC, COLC )
        implicit none
        ! DP = (C)^T U_LONG
        INTEGER, intent( in ) :: U_NONODS, NDIM, NPHASE
        REAL, DIMENSION( :, :, :, : ), intent( in ) :: U_LONG
        REAL, DIMENSION( :, : ), intent( inout )  :: DP
        REAL, DIMENSION( :, :, : ), intent( in ) :: C
        INTEGER, DIMENSION( : ), intent( in ) ::FINDC
        INTEGER, DIMENSION( : ), intent( in ) :: COLC
        ! Local variables
        INTEGER :: U_INOD, COUNT, P_JNOD, IPHASE, IDIM

        DP = 0.0

        Loop_VelNodes: DO U_INOD = 1, U_NONODS

            Loop_Crow: DO COUNT = FINDC( U_INOD ), FINDC( U_INOD + 1 ) - 1
                P_JNOD = COLC( COUNT )

                Loop_Phase: DO IPHASE = 1, NPHASE
                    Loop_Dim: DO IDIM = 1, NDIM
                        DP( :, P_JNOD ) = DP( :, P_JNOD ) + C( IDIM, IPHASE, COUNT ) * U_LONG( :, IDIM, IPHASE, U_INOD )
                    END DO Loop_Dim
                END DO Loop_Phase

            END DO Loop_Crow

        END DO Loop_VelNodes

        RETURN

    END SUBROUTINE CT_MULT_WITH_C_MANY







    SUBROUTINE ULONG_2_UVW( U, V, W, UP, U_NONODS, NDIM, NPHASE)
        implicit none
        INTEGER, intent( in ) :: U_NONODS, NDIM, NPHASE
        REAL, DIMENSION( : ), intent( inout ) :: U,V,W
        REAL, DIMENSION( : ), intent( in ) :: UP
        ! local variables...
        INTEGER :: IPHASE
        DO IPHASE = 1, NPHASE
            U( 1 + ( IPHASE - 1 ) * U_NONODS : U_NONODS + ( IPHASE - 1 ) * U_NONODS ) = &
                UP( 1 + ( IPHASE - 1 ) * NDIM * U_NONODS : U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS )
            V( 1 + ( IPHASE - 1 ) * U_NONODS : U_NONODS + ( IPHASE - 1 ) * U_NONODS ) = &
                UP( 1 + U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS : 2 * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS )
            IF( NDIM >= 3 ) &
                W( 1 + ( IPHASE - 1 ) * U_NONODS : U_NONODS + ( IPHASE - 1 ) * U_NONODS ) = &
                UP( 1 + 2 * U_NONODS + ( IPHASE - 1) * NDIM * U_NONODS : 3 * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS )
        END DO
        RETURN
    END SUBROUTINE ULONG_2_UVW

    subroutine posinmat( posmat, globi, globj, &
        findrm, colm )
        ! Find position in matrix POSMAT which has column GLOBJ
        implicit none
        integer, intent( inout ) :: posmat
        integer, intent( in ) :: globi, globj
        integer, dimension( : ), intent( in ) :: findrm
        integer, dimension( : ), intent( in ) :: colm
        ! Local variables
        integer, parameter :: nmax = 1000000
        integer :: inum, lower, upper, count

        lower = findrm( globi )
        upper = findrm( globi + 1 ) - 1

        count = 1
        Loop_While: do while ( count <= nmax )
            inum = lower + ( upper - lower + 1 ) / 2

            if( globj >= colm( inum ) ) then
                lower = inum
            else
                upper = inum
            end if

            if( ( upper - lower ) <= 1 ) then
                if( globj == colm( lower ) ) then
                    posmat = lower
                else
                    posmat = upper
                end if
                return
            end if

            count = count + 1

        end do Loop_While

        return
    end subroutine posinmat

    subroutine assemble_global_multiphase_csr(global_csr,&
        block_csr,dense_block_matrix,block_to_global,global_dense_block)

        real, dimension(:), intent(out)    ::  global_csr
        real, dimension(:), intent(in)     :: block_csr
        real, dimension(:,:,:), intent(in) :: dense_block_matrix
        integer, dimension(:), intent(in)  :: block_to_global
        integer, dimension(:,:), intent(in)  :: global_dense_block

        integer :: node, jphase, node_count, nphase


        ewrite(3,*), "In  assemble_global_multiphase_csr"

        ! copy the block_csr to global using the assigned map

        ewrite(3,*) 'size(global_csr), size(block_csr)', size(global_csr), size(block_csr)
        global_csr=0.0
        global_csr(block_to_global)=block_csr(:)

        ! now for the dense block
        node_count=size(dense_block_matrix,3)
        nphase=size(dense_block_matrix,2)

        do node=1,node_count
            do jphase=1,nphase
                global_csr(global_dense_block(jphase,node):&
                    global_dense_block(jphase,node)+nphase-1)=&
                    global_csr(global_dense_block(jphase,node):&
                    global_dense_block(jphase,node)+nphase-1)+dense_block_matrix(:,jphase,node)
            end do
        end do


        ewrite(3,*), "Leaving assemble_global_multiphase_csr"


    end subroutine assemble_global_multiphase_csr


    subroutine allocate_global_multiphase_petsc_csr(global_petsc,&
        sparsity,tracer)

        type(petsc_csr_matrix)    ::  global_petsc
        type(tensor_field) :: tracer

        integer :: nphase
        type(csr_sparsity) :: sparsity

        nphase=tracer%dim(2)

        ewrite(3,*), "In assemble_global_multiphase_petsc_csr"

        !call allocate(global_petsc,sparsity,[nphase,nphase],"ACV",.true.)
        call allocate(global_petsc,sparsity,[nphase,nphase],"ACV",.false.,.false.)
        call zero(global_petsc)

        ewrite(3,*), "Leaving allocate_global_multiphase_petsc_csr"

    end subroutine allocate_global_multiphase_petsc_csr

    function allocate_momentum_matrix(sparsity,velocity) result(Mat)
        type(csr_sparsity), intent (inout) :: sparsity
        type(tensor_field), intent (inout) :: velocity
        type(halo_type), pointer:: halo
        type(petsc_csr_matrix) :: mat
        integer :: ierr


        if (associated(velocity%mesh%halos)) then
            halo => velocity%mesh%halos(2)
        else
            nullify(halo)
        end if

        mat%name="MomentumMatrix"

        if (associated(halo)) then
            allocate(mat%row_halo)
            mat%row_halo = halo
            call incref(mat%row_halo)
            allocate(mat%column_halo)
            mat%column_halo = halo
            call incref(mat%column_halo)
        else
            nullify(mat%row_halo)
            nullify(mat%column_halo)
        end if

        call allocate(mat%row_numbering,node_count(velocity),&
            product(velocity%dim),halo)
        call allocate(mat%column_numbering,node_count(velocity),&
            product(velocity%dim),halo)

        if (.not. IsParallel()) then
            mat%M=full_CreateSeqAIJ(sparsity, mat%row_numbering, &
                mat%column_numbering,.false.)
        else
            mat%M=full_CreateMPIAIJ(sparsity, mat%row_numbering, &
                mat%column_numbering,.false.)
        end if

        call MatSetOption(mat%M, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)
        call MatSetOption(mat%M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE, ierr)
        call MatSetOption(mat%M, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
        nullify(mat%refcount)
        call addref_petsc_csr_matrix(mat)


    end function allocate_momentum_matrix

end module matrix_operations


