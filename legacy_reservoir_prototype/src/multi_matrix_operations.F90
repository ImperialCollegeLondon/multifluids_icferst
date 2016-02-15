
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
    use multiphase_caching, only : get_caching_level, test_caching_level, reshape_vector2pointer
    use global_parameters, only : FIELD_NAME_LEN, is_porous_media
    use boundary_conditions
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
        INTEGER :: ICOL, IM, JM
        !      INTEGER , DIMENSION(N,N) :: IPIV
        INTEGER , DIMENSION(N ) :: IPIV
        REAL, DIMENSION( N, N ) :: MAT, MAT2
        REAL, DIMENSION( N ) :: X, B

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

    SUBROUTINE MATINVold( A, N, NMAX, MAT, MAT2, X, B)
        ! This sub finds the inverse of the matrix A and puts it back in A.
        ! MAT, MAT2, X and B are working vectors.
        IMPLICIT NONE
        INTEGER, intent( in ) :: N, NMAX
        REAL, DIMENSION( :, : ), intent( inout ) ::  A
        REAL, DIMENSION( :, : ), intent( inout ) :: MAT, MAT2
        REAL, DIMENSION( : ), intent( inout ) :: X, B
        ! Local variables
        INTEGER :: ICOL, IM, JM
        INTEGER , DIMENSION(N) :: IPIV

        ! Solve MAT XL=BX (NB BDIAG is overwritten)
        ICOL = 1
        B( 1 : N ) = 0.
        B( ICOL ) = 1.

        MAT = A( 1:N,1:N )

        CALL SMLINNGOT( MAT, X, B, N, N,IPIV, .FALSE. ) ! X contains the column ICOL of inverse

        A( :, ICOL ) = X( : )

        DO ICOL = 2, N ! Form column ICOL of the inverse.
            B = 0.
            B( ICOL ) = 1.0 ! Solve MAT X=B (NB MAT is overwritten).

            CALL SMLINNGOT( MAT, X, B, N, N,IPIV, .TRUE. ) ! X contains the column ICOL of inverse

            A( :, ICOL ) = X( : )
        END DO

        RETURN
    END SUBROUTINE MATINVold

    !sprint_to_do!add some comments
    SUBROUTINE SMLINNGOT( A, X, B, NMX, N, IPIV, GOTDEC )
        IMPLICIT NONE
        INTEGER :: NMX, N
        REAL, DIMENSION( :, : ), intent( inout ) :: A
        real, DIMENSION( : ), intent( inout ) :: X ! inout as n might be < nmx
        real, DIMENSION( : ), intent( in ) ::  B
        LOGICAL, intent( in ) :: GOTDEC
        ! Local
        REAL :: R
        INTEGER :: K, I, J
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



    SUBROUTINE COLOR_GET_CMC_PHA( CV_NONODS, U_NONODS, NDIM, NPHASE, NPRES, &
        NCOLC, FINDC, COLC, &
        INV_PIVIT_MAT,  &
        TOTELE, U_NLOC, U_NDGLN, &
        NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
        CMC_petsc, CMC_PRECON, IGOT_CMC_PRECON, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
        MASS_PIPE, MASS_CVFEM2PIPE, MASS_CVFEM2PIPE_TRUE, &
        got_free_surf,  MASS_SUF, &
        C, CT, storage_state, indx, halos, symmetric_P )
        !use multiphase_1D_engine
        !Initialize the momentum equation (CMC) and introduces the corresponding values in it.
        implicit none
        ! form pressure matrix CMC using a colouring approach
        LOGICAL, PARAMETER :: PIPES_1D=.TRUE.
        INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NPRES, NCOLC, &
            TOTELE, U_NLOC, NCOLCT, NCOLCMC, IGOT_CMC_PRECON
        LOGICAL, intent( in ) :: got_free_surf, symmetric_P
        INTEGER, DIMENSION( : ), intent( in ) ::FINDC
        INTEGER, DIMENSION( : ), intent( in ) :: COLC
        REAL, DIMENSION( :, :, : ), intent( in ) :: INV_PIVIT_MAT
        INTEGER, DIMENSION( : ), intent( in ) ::  U_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
        INTEGER, DIMENSION( : ), intent( in ) :: COLCT
        REAL, DIMENSION( :, : ), intent( in ) :: DIAG_SCALE_PRES
        REAL, DIMENSION( :, :, : ), intent( in ) :: DIAG_SCALE_PRES_COUP, INV_B
        type(petsc_csr_matrix), intent(inout)::  CMC_petsc
        REAL, DIMENSION( :, :, : ), intent( inout ) :: CMC_PRECON
        REAL, DIMENSION( : ), intent( in ) :: MASS_MN_PRES
        REAL, DIMENSION( : ), intent( in ) :: MASS_SUF
        REAL, DIMENSION( : ), intent( in ) :: MASS_PIPE, MASS_CVFEM2PIPE, MASS_CVFEM2PIPE_TRUE
        INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
        INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
        REAL, DIMENSION( :, :, : ), intent( in ) :: C
        REAL, DIMENSION( :, :, : ), intent( in ) :: CT
        type( state_type ), intent( inout ):: storage_state
        integer, intent(inout) :: indx
        type(halo_type), pointer :: halos
        ! Local variables
        REAL, PARAMETER :: INFINY = 1.0E+10
        LOGICAL :: UNDONE, LCOL
        logical, DIMENSION( : ), allocatable :: NEED_COLOR
        logical, DIMENSION( : ), allocatable :: to_color
        logical, DIMENSION( :, : ), allocatable :: COLOR_VEC_MANY_LOGICAL
        REAL, DIMENSION( : ), allocatable :: COLOR_VEC, CMC_COLOR_VEC, CMC_COLOR_VEC2
        REAL, DIMENSION( : ), allocatable :: CDP, DU, DV, DW, DU_LONG
        REAL, DIMENSION( :, : ), allocatable :: COLOR_VEC_MANY, CMC_COLOR_VECC_MANY, CMC_COLOR_VEC2C_MANY
        REAL, DIMENSION( :, : ), allocatable :: CDPC_MANY, DUC_MANY, DVC_MANY, DWC_MANY, DU_LONGC_MANY
        INTEGER :: NCOLOR, COUNT, COUNT2, IDIM, IPHASE, CV_COLJ, U_JNOD, CV_JNOD2
        INTEGER :: I, ELE,u_inod,u_nod
        INTEGER, DIMENSION( : ), allocatable :: ndpset
        !Variables for CMC_petsc
        integer, dimension( : ), allocatable :: dnnz
        integer :: cmc_rows

        !Initialize CMC_petsc

        call zero( CMC_petsc )

        allocate(ndpset(npres)) ; ndpset(:)=-1
        if (ndpset(1)<0) then
            call get_option( '/material_phase[0]/scalar_field::Pressure/' // &
                'prognostic/reference_node', ndpset(1), default = -1 )
            ndpset(2:npres)=ndpset(1)
        endif

        if (isparallel()) then
            if (GetProcNo()>1) ndpset=0
        end if

           !Currently this is false, over-ride switch for porous media since it is slower than memory hungry
        IF ( ( test_caching_level(6) .or.is_porous_media ) .and. npres==1) THEN
            ! Fast but memory intensive...
            CALL COLOR_GET_CMC_PHA_FAST( CV_NONODS, U_NONODS, NDIM, NPHASE, NPRES, &
                NCOLC, FINDC, COLC, &
                INV_PIVIT_MAT,  &
                TOTELE, U_NLOC, U_NDGLN, &
                NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
                CMC_petsc, CMC_PRECON, IGOT_CMC_PRECON, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
                MASS_PIPE, MASS_CVFEM2PIPE, MASS_CVFEM2PIPE_TRUE, &
                got_free_surf,  MASS_SUF, &
                C, CT, ndpset, storage_state, indx, symmetric_P )
        ELSE
            ! Slow but memory efficient...
            CALL COLOR_GET_CMC_PHA_SLOW( CV_NONODS, U_NONODS, NDIM, NPHASE, NPRES, &
                NCOLC, FINDC, COLC, &
                INV_PIVIT_MAT,  &
                TOTELE, U_NLOC, U_NDGLN, &
                NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
                CMC_petsc, CMC_PRECON, IGOT_CMC_PRECON, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
                MASS_PIPE, MASS_CVFEM2PIPE, MASS_CVFEM2PIPE_TRUE, &
                got_free_surf,  MASS_SUF, &
                C, CT, ndpset, symmetric_P )
        END IF


        !Re-assemble just in case
        CMC_petsc%is_assembled=.false.
        call assemble( CMC_petsc )

      !Final step to prepare the CMC_petsc
      !call assemble( CMC_petsc )

    END SUBROUTINE COLOR_GET_CMC_PHA




    !sprint_to_do!make internal subroutine
    SUBROUTINE COLOR_GET_CMC_PHA_FAST( CV_NONODS, U_NONODS, NDIM, NPHASE, NPRES, &
        NCOLC, FINDC, COLC, &
        INV_PIVIT_MAT,  &
        TOTELE, U_NLOC, U_NDGLN, &
        NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
        CMC_petsc, CMC_PRECON, IGOT_CMC_PRECON, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
        MASS_PIPE, MASS_CVFEM2PIPE, MASS_CVFEM2PIPE_TRUE,  &
        got_free_surf,  MASS_SUF, &
        C, CT, ndpset, storage_state, indx, symmetric_P )

        implicit none
        ! form pressure matrix CMC using a colouring approach
        INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NPRES, NCOLC, &
            TOTELE, U_NLOC, NCOLCT, NCOLCMC, IGOT_CMC_PRECON
        LOGICAL, PARAMETER :: PIPES_1D=.TRUE.
        LOGICAL, intent( in ) :: got_free_surf, symmetric_P
        INTEGER, DIMENSION( : ), intent( in ) :: ndpset
        INTEGER, DIMENSION( : ), intent( in ) ::FINDC
        INTEGER, DIMENSION( : ), intent( in ) :: COLC
        REAL, DIMENSION( :, :, : ), intent( in ) :: INV_PIVIT_MAT
        INTEGER, DIMENSION( : ), intent( in ) ::  U_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
        INTEGER, DIMENSION( : ), intent( in ) :: COLCT
        REAL, DIMENSION( :, : ), intent( in ) :: DIAG_SCALE_PRES
        REAL, DIMENSION( :, :, : ), intent( in ) :: DIAG_SCALE_PRES_COUP, INV_B
        type(petsc_csr_matrix), intent(inout)::  CMC_petsc
        REAL, DIMENSION( :, :, : ), intent( inout ) :: CMC_PRECON
        REAL, DIMENSION( : ), intent( in ) :: MASS_MN_PRES
        REAL, DIMENSION( : ), intent( in ) :: MASS_SUF
        REAL, DIMENSION( : ), intent( in ) :: MASS_PIPE, MASS_CVFEM2PIPE, MASS_CVFEM2PIPE_TRUE
        INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
        INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
        REAL, DIMENSION( :, :, : ), intent( in ) :: C
        REAL, DIMENSION( :, :, : ), intent( in ) :: CT
        type( state_type ), intent( inout ):: storage_state
        integer, intent(inout) :: indx
        ! Local variables
        INTEGER, PARAMETER :: MX_NCOLOR = 1000
        REAL, PARAMETER :: INFINY = 1.0E+10
        LOGICAL :: LCOL, EXPLICIT_PIPES2
        LOGICAL, DIMENSION( CV_NONODS ) :: COLOR_LOGICAL
        INTEGER, DIMENSION( CV_NONODS ) :: ICOLOR
        INTEGER, DIMENSION( : ), allocatable :: COLOR_IN_ROW, COLOR_IN_ROW2
        REAL, DIMENSION( :, : ), allocatable :: COLOR_VEC_MANY
        REAL, DIMENSION( :, :, :, : ), allocatable :: CDP_MANY, DU_LONG_MANY
        REAL, DIMENSION( :, :, : ), allocatable :: CMC_COLOR_VEC_MANY, CMC_COLOR_VEC2_MANY, &
            CMC_COLOR_VEC_MANY_PHASE, CMC_COLOR_VEC2_MANY_PHASE
        REAL, DIMENSION( : ), allocatable :: CMC_COLOR_VEC_MANY_PHASE_SHORT
        REAL, DIMENSION( :, : ), allocatable :: RINV_B_COLNS_ZEROED, CMC_COLOR_VEC_PRES_MANY_SHORT
        INTEGER :: CV_NOD, CV_JNOD, COUNT, COUNT2, COUNT3, IDIM, IPHASE, CV_COLJ, U_JNOD, CV_JNOD2
        INTEGER :: MAX_COLOR_IN_ROW, ICHOOSE, KVEC, I, ELE, U_INOD, U_NOD, ICAN_COLOR, MX_COLOR, NOD_COLOR
        INTEGER :: NCOLOR, ierr, j, N_IN_PRES, i_indx, j_indx, IPRES, JPRES
        REAL :: RSUM, RSUM_SUF
        !Variables to store things in state
        type(mesh_type), pointer :: fl_mesh
        type(mesh_type) :: Auxmesh
        type(scalar_field), target :: targ_icolor
        real, pointer, dimension(:) :: pointer_icolor

        IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON = 0.0
        N_IN_PRES = NPHASE / NPRES
        EXPLICIT_PIPES2 = .true.

        MAX_COLOR_IN_ROW = 0
        DO CV_NOD = 1, CV_NONODS
            MAX_COLOR_IN_ROW = MAX( MAX_COLOR_IN_ROW, FINDCMC( CV_NOD + 1 ) - FINDCMC( CV_NOD ) )
        END DO

        ALLOCATE( COLOR_IN_ROW( MAX_COLOR_IN_ROW**2 ) )
        ALLOCATE( COLOR_IN_ROW2( MAX_COLOR_IN_ROW**2 ) )
        IF ( indx /= 0 ) THEN!coloring stored

            pointer_icolor =>  storage_state%scalar_fields(indx)%ptr%val
            !Recover the calculated coloring
            !THIS IS NOT EFFICIENT, IF WE FINALLY DECIDE TO USE ONLY THE FAST VERSION
            !WE SHOULD CONVERT ICOLOR AND NCOLOR INTO POINTERS!
            ICOLOR = pointer_icolor(1:CV_NONODS)
            NCOLOR = pointer_icolor(CV_NONODS+1)

        ELSE

               !Prepare the variable inside state
            if (has_scalar_field(storage_state, "Col_fast")) then
                !If we are recalculating due to a mesh modification then
                !we return to the original situation
                call remove_scalar_field(storage_state, "Col_fast")
            end if
            !Get mesh file just to be able to allocate the fields we want to store
            fl_mesh => extract_mesh( storage_state, "FakeMesh" )
            Auxmesh = fl_mesh
            !The number of nodes I want does not coincide
            Auxmesh%nodes = CV_NONODS+1!in the +1 we store NCOLOR
            call allocate (targ_icolor, Auxmesh)
            !Now we insert them in state and store the index
            call insert(storage_state, targ_icolor, "Col_fast")
            call deallocate (targ_icolor)
            !call deallocate(Auxmesh)
            indx = size(storage_state%scalar_fields)
            !Get from state
            pointer_icolor =>  storage_state%scalar_fields(indx)%ptr%val

            COLOR_LOGICAL = .FALSE.

            NCOLOR=0
            ICOLOR = 0
            Loop_CVNOD7: DO CV_NOD = 1, CV_NONODS
                ! Color this node CV_NOD

                MX_COLOR=0
                ! use a distance-2 colouring...
                Loop_Row7: DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                    CV_JNOD = COLCMC( COUNT )
                    Loop_Row2: DO COUNT2 = FINDCMC( CV_JNOD ), FINDCMC( CV_JNOD + 1 ) - 1
                        CV_JNOD2 = COLCMC( COUNT2 )
                        IF(CV_NOD.NE.CV_JNOD2) THEN
                            NOD_COLOR=ICOLOR(CV_JNOD2)
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
                ICOLOR(CV_NOD)=ICAN_COLOR
                NCOLOR=MAX(NCOLOR,ICAN_COLOR)
                ! Reset the colors for the next node...
                COLOR_LOGICAL(1:max(1,MX_COLOR))=.FALSE.
            END DO Loop_CVNOD7


            IF ( NCOLOR > MX_NCOLOR ) THEN
                EWRITE(-1, *) 'NOT ENOUGH COLOURS STOPPING - NEED TO MAKE MX_COLOR BIGGER'
                STOP 281
            END IF

            !##Store the calculated coloring##
            pointer_icolor(1:CV_NONODS) = ICOLOR
            pointer_icolor(CV_NONODS+1) = NCOLOR

        END IF ! SAVED_CMC_COLOR


        ALLOCATE( COLOR_VEC_MANY( NCOLOR, CV_NONODS ) )
        COLOR_VEC_MANY = 0.0

        Loop_CVNOD: DO CV_NOD = 1, CV_NONODS
            COLOR_VEC_MANY( ICOLOR( CV_NOD ), CV_NOD ) = 1.0
        END DO Loop_CVNOD

        ! we use the same colouring for each pressure variable when npres>1.
        ALLOCATE( CDP_MANY( NCOLOR, NDIM, NPHASE, U_NONODS ) )

        CALL C_MULT_MANY( CDP_MANY, COLOR_VEC_MANY, CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLOR, &
            C, NCOLC, FINDC, COLC )
      
        ! DU_LONG = INV_PIVIT_MAT * CDP
        ALLOCATE( DU_LONG_MANY( NCOLOR, NDIM, NPHASE, U_NONODS ) )

        CALL PHA_BLOCK_MAT_VEC_MANY( DU_LONG_MANY, INV_PIVIT_MAT, CDP_MANY, U_NONODS, NDIM, NPHASE, NCOLOR, &
            TOTELE, U_NLOC, U_NDGLN )

                 ! NB. P_RHS = CT * U + CV_RHS
                 ! DU_LONG = CDP






        IF ( NPRES > 1 .AND. .NOT.EXPLICIT_PIPES2 ) THEN


            ALLOCATE( CMC_COLOR_VEC_MANY_PHASE( NCOLOR, NPHASE, CV_NONODS ) )
            ALLOCATE( CMC_COLOR_VEC_MANY( NCOLOR, NPRES, CV_NONODS ) )

            DO IPHASE = 1, NPHASE
                CALL CT_MULT_MANY( CMC_COLOR_VEC_MANY_PHASE(:,IPHASE,:), &
                    DU_LONG_MANY(:,:,IPHASE:IPHASE,:), &
                    CV_NONODS, U_NONODS, NDIM, 1, NCOLOR, &
                    CT(:,IPHASE:IPHASE,:), NCOLCT, FINDCT, COLCT )
            END DO

            ALLOCATE( RINV_B_COLNS_ZEROED(NPHASE,NPHASE) )
            ALLOCATE( CMC_COLOR_VEC_MANY_PHASE_SHORT(NPHASE) )
            ALLOCATE( CMC_COLOR_VEC_PRES_MANY_SHORT(NCOLOR,NPRES) )


            DO JPRES=1,NPRES

                DO CV_NOD = 1, CV_NONODS
                    RINV_B_COLNS_ZEROED(:,:) = 0.0
                    RINV_B_COLNS_ZEROED(:,  1+(JPRES-1)*N_IN_PRES:JPRES*N_IN_PRES ) = INV_B(:,  1+(JPRES-1)*N_IN_PRES:JPRES*N_IN_PRES,  CV_NOD)
                    DO I = 1, NCOLOR
                        !                  CMC_COLOR_VEC_MANY_PHASE(I,:,CV_NOD) = MATMUL( RINV_B_COLNS_ZEROED(:,:), CMC_COLOR_VEC_MANY_PHASE(I,:,CV_NOD) )
                        CMC_COLOR_VEC_MANY_PHASE_SHORT(:) = MATMUL( RINV_B_COLNS_ZEROED(:,:), CMC_COLOR_VEC_MANY_PHASE(I,:,CV_NOD) )
                        DO IPRES = 1, NPRES
                            CMC_COLOR_VEC_PRES_MANY_SHORT(I,IPRES) = SUM( CMC_COLOR_VEC_MANY_PHASE_SHORT(1+(IPRES-1)*N_IN_PRES:IPRES*N_IN_PRES) )
                        END DO
                    END DO

                    DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                        CV_JNOD = COLCMC( COUNT )

                        DO IPRES = 1, NPRES
                            call addto( CMC_petsc, blocki = IPRES, blockj = JPRES, i = cv_nod, j = CV_JNOD, &
                                val = dot_product(CMC_COLOR_VEC_PRES_MANY_SHORT( :, IPRES ), COLOR_VEC_MANY( :, CV_JNOD ) ))
                        END DO
                    END DO

                END DO

            END DO ! ENDOF DO JPRES=1,NPRES



            IF ( IGOT_CMC_PRECON /= 0 ) THEN
                ALLOCATE( CMC_COLOR_VEC2_MANY_PHASE( NCOLOR, NPHASE, CV_NONODS ) )
                ALLOCATE( CMC_COLOR_VEC2_MANY( NCOLOR, NPRES, CV_NONODS ) )

                DO IPHASE = 1, NPHASE
                    CALL CT_MULT_WITH_C_MANY( CMC_COLOR_VEC2_MANY_PHASE(:,IPHASE,:), &
                        DU_LONG_MANY(:,:,IPHASE:IPHASE,:), &
                        CV_NONODS, U_NONODS, NDIM, 1, NCOLOR, &
                        C(:,IPHASE:IPHASE,:), NCOLC, FINDC, COLC )
                END DO


                DO JPRES=1,NPRES

                    DO CV_NOD = 1, CV_NONODS
                        RINV_B_COLNS_ZEROED(:,:) = 0.0
                        RINV_B_COLNS_ZEROED(:,  1+(JPRES-1)*N_IN_PRES:JPRES*N_IN_PRES ) = INV_B(:,  1+(JPRES-1)*N_IN_PRES:JPRES*N_IN_PRES,  CV_NOD)
                        DO I = 1, NCOLOR
                            ! The following uses CMC_COLOR_VEC2_MANY_PHASE:
                            CMC_COLOR_VEC_MANY_PHASE_SHORT(:) = MATMUL( RINV_B_COLNS_ZEROED(:,:), CMC_COLOR_VEC2_MANY_PHASE(I,:,CV_NOD) )
                            DO IPRES = 1, NPRES
                                CMC_COLOR_VEC_PRES_MANY_SHORT(I,IPRES) = SUM( CMC_COLOR_VEC_MANY_PHASE_SHORT(1+(IPRES-1)*N_IN_PRES:IPRES*N_IN_PRES) )
                            END DO
                        END DO

                        DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                            CV_JNOD = COLCMC( COUNT )

                            DO IPRES = 1, NPRES
                                CMC_PRECON( IPRES, JPRES, COUNT ) = CMC_PRECON( IPRES, JPRES, COUNT ) + &
                                    sum( CMC_COLOR_VEC_PRES_MANY_SHORT( :, IPRES )* COLOR_VEC_MANY( :, CV_JNOD ) )

                            END DO
                        END DO

                    END DO

                END DO ! ENDOF DO JPRES=1,NPRES

                CMC_COLOR_VEC2_MANY=0.0

            END IF

            ! Re-set CMC_COLOR_VEC_MANY & CMC_COLOR_VEC_MANY2 as we have already added the
            ! contributions up to this point into the matrix.
            CMC_COLOR_VEC_MANY=0.0
        ELSE


            ALLOCATE( CMC_COLOR_VEC_MANY( NCOLOR, NPRES, CV_NONODS ) )
            DO IPRES = 1, NPRES
                CALL CT_MULT_MANY( CMC_COLOR_VEC_MANY(:,IPRES,:), &
                    DU_LONG_MANY(:,:,1+(IPRES-1)*N_IN_PRES:IPRES*N_IN_PRES,:), &
                    CV_NONODS, U_NONODS, NDIM, N_IN_PRES, NCOLOR, &
                    CT(:,1+(IPRES-1)*N_IN_PRES:IPRES*N_IN_PRES,:), NCOLCT, FINDCT, COLCT )
            END DO

            IF ( IGOT_CMC_PRECON /= 0 ) THEN
                ALLOCATE( CMC_COLOR_VEC2_MANY( NCOLOR, NPRES, CV_NONODS ) )
                DO IPRES = 1, NPRES
                    CALL CT_MULT_WITH_C_MANY( CMC_COLOR_VEC2_MANY(:,IPRES,:), &
                        DU_LONG_MANY(:,:,1+(IPRES-1)*N_IN_PRES:IPRES*N_IN_PRES, :), &
                        CV_NONODS, U_NONODS, NDIM, N_IN_PRES, NCOLOR, &
                        C(:,1+(IPRES-1)*N_IN_PRES:IPRES*N_IN_PRES,:), NCOLC, FINDC, COLC )
                END DO
            END IF


        END IF



        ! Matrix vector involving the mass diagonal term
        DO CV_NOD = 1, CV_NONODS
            RSUM = 0.0
            RSUM_SUF = 0.0
            DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                CV_JNOD = COLCMC( COUNT )

                DO IPRES = 1, NPRES

                    IF((NPRES>1).AND.PIPES_1D) THEN
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
                DO IPRES=1,NPRES
                    CMC_COLOR_VEC2_MANY( :, IPRES, CV_NOD ) = CMC_COLOR_VEC2_MANY( :, IPRES, CV_NOD ) &
                        + ( DIAG_SCALE_PRES( IPRES, CV_NOD ) * RSUM + RSUM_SUF ) * COLOR_VEC_MANY( :, CV_NOD )
                END DO
            END IF
        END DO

        !Put into matrix CMC
        if ( .not.symmetric_P ) then

            ! original method
            DO CV_NOD = 1, CV_NONODS
                DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                    CV_JNOD = COLCMC( COUNT )
                    DO IPRES = 1, NPRES
                        JPRES = IPRES ! Add contributions to the block diagonal only.
                        call addto( CMC_petsc, blocki = IPRES, blockj = JPRES, i = cv_nod, j = CV_JNOD, &
                            val = dot_product(CMC_COLOR_VEC_MANY( :, IPRES, CV_NOD ), COLOR_VEC_MANY( :, CV_JNOD ) ))
                        !CMC( COUNT ) = CMC( COUNT ) + sum(CMC_COLOR_VEC_MANY( :, CV_NOD ) * COLOR_VEC_MANY( :, CV_JNOD ))
                        IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON( IPRES, JPRES, COUNT ) = CMC_PRECON( IPRES, JPRES, COUNT ) + &
                            sum(CMC_COLOR_VEC2_MANY( :, IPRES, CV_NOD ) * COLOR_VEC_MANY( :, CV_JNOD ))
                    END DO
                END DO
            END DO

        else

            ! symmetric P matrix
            DO CV_NOD = 1, CV_NONODS
                DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                    CV_JNOD = COLCMC( COUNT )
                    DO IPRES = 1, NPRES
                        JPRES = IPRES ! Add contributions to the block diagonal only.
                        call addto( CMC_petsc, blocki = IPRES, blockj = JPRES, i = cv_nod, j = CV_JNOD, &
                            val = dot_product(CMC_COLOR_VEC2_MANY( :, IPRES, CV_NOD ), COLOR_VEC_MANY( :, CV_JNOD ) ))
                        !CMC( COUNT ) = CMC( COUNT ) + sum(CMC_COLOR_VEC_MANY( :, CV_NOD ) * COLOR_VEC_MANY( :, CV_JNOD ))
                        IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON( IPRES, JPRES, COUNT ) = CMC_PRECON( IPRES, JPRES, COUNT ) + &
                            sum(CMC_COLOR_VEC2_MANY( :, IPRES, CV_NOD ) * COLOR_VEC_MANY( :, CV_JNOD ))
                    END DO
                END DO
            END DO

        end if


        ! the matrix coupling term place immediately into matrix and not through colouring...
        ! Matrix vector involving the mass diagonal term
        IF ( NPRES > 1) THEN
            DO CV_NOD = 1, CV_NONODS
                RSUM = 0.0
                RSUM_SUF = 0.0
                DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                    CV_JNOD = COLCMC( COUNT )
                    DO IPRES = 1, NPRES
                        DO JPRES = 1, NPRES
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
        END IF ! ENDOF IF(NPRES > 1) THEN

        !If we have a reference node with pressure zero we impose that here.
        DO IPRES = 1, NPRES
            IF ( NDPSET(IPRES) > 0 ) THEN
                CV_NOD = NDPSET(IPRES)
                i_indx = CMC_petsc%row_numbering%gnn2unn( cv_nod, ipres )
                DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                    CV_JNOD = COLCMC( COUNT )
                    IF ( CV_JNOD /= CV_NOD ) THEN
                        do jpres = 1, npres
                            j_indx = CMC_petsc%column_numbering%gnn2unn( cv_jnod, jpres )
                            call MatSetValue(CMC_petsc%M, i_indx, j_indx, 0.0,INSERT_VALUES, ierr) ! not the diagonal
                            !CMC( COUNT ) = 0.0 ! not the diagonal
                            IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON( IPRES, JPRES, COUNT ) = 0.0
                            DO COUNT2 = FINDCMC( CV_JNOD ), FINDCMC( CV_JNOD + 1 ) - 1
                                CV_JNOD2 = COLCMC( COUNT2 )
                                IF ( CV_JNOD2 /= CV_NOD ) then
                                    i_indx = CMC_petsc%row_numbering%gnn2unn( cv_jnod, ipres )
                                    j_indx = CMC_petsc%column_numbering%gnn2unn( CV_JNOD2, jpres )
                                    call MatSetValue(CMC_petsc%M, i_indx, j_indx, 0.0, INSERT_VALUES, ierr) ! not the diagonal
                                    IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON( ipres, jpres, COUNT2 ) = 0.0
                                END IF
                            END DO
                        end do ! do jpres=1,npres
                    END IF
                END DO
            END IF
        END DO

        !Re-assemble
        CMC_petsc%is_assembled=.false.
        call assemble( CMC_petsc )

        IF ( IGOT_CMC_PRECON /= 0 ) deallocate(CMC_COLOR_VEC2_MANY)

        deallocate(COLOR_IN_ROW, COLOR_IN_ROW2 )
        deallocate(COLOR_VEC_MANY,CMC_COLOR_VEC_MANY, CDP_MANY, DU_LONG_MANY)

        RETURN
    END SUBROUTINE COLOR_GET_CMC_PHA_FAST


    !sprint_to_do!move to library and aid subroutines as well
    recursive  subroutine quicksort(vec,n)

        implicit none

        integer, intent(in) :: n
        integer, dimension(:), intent(inout) :: vec(n)
        integer :: ii

        if (n>20) then
            ii=partition(vec)
            call quicksort(vec(1:ii-1),ii-1)
            call quicksort(vec(ii:n),n-ii+1)
        else
            call insertion_sort(vec,n)
        end if

    end subroutine quicksort
        
    integer function partition(v)
          
        implicit none

        integer, intent(inout), dimension(:) :: v
        integer :: i,j, pivot, temp

        i=0
        j=size(v)+1

        pivot=v(j/2)
          
        do while (i < j)
            j= j-1
            do while ( v(j) > pivot)
                j=j-1
            end do
            i = i + 1
            do while  ( v(i) < pivot )
                i = i + 1
            end do
            if ( i < j ) then
                temp=v(i); v(i)=v(j); v(j)= temp
            elseif ( i == j ) then
                partition = i+1
                return
            else
                partition = i
            end if
        end do

    end function partition
           
    subroutine insertion_sort(vec,n)
          
        implicit none

        integer :: n
        integer, dimension(n) :: vec(n)

        integer :: i ,j, temp


        do i = 2, N
            j = i - 1
            temp = vec(i)

            do while ( j> 0 )
                if ( vec(j) <= temp ) exit
                vec(j+1)=vec(j)
                j=j-1
            end do

            vec(j+1) = temp
        end do

    end subroutine insertion_sort
             

    !sprint_to_do!make internal subroutine. try to delete this and only use the fast version
    !remove if we finally remove this subroutine all the subroutines that are only called in here
    SUBROUTINE COLOR_GET_CMC_PHA_SLOW( CV_NONODS, U_NONODS, NDIM, NPHASE, NPRES, &
        NCOLC, FINDC, COLC, &
        INV_PIVIT_MAT,  &
        TOTELE, U_NLOC, U_NDGLN, &
        NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, DIAG_SCALE_PRES_COUP, INV_B, &
        CMC_petsc, CMC_PRECON, IGOT_CMC_PRECON, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
        MASS_PIPE, MASS_CVFEM2PIPE, MASS_CVFEM2PIPE_TRUE, &
        got_free_surf,  MASS_SUF, &
        C, CT, ndpset, symmetric_P )
        !use multiphase_1D_engine

        implicit none
        ! form pressure matrix CMC using a colouring approach
        INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NPRES, NCOLC, &
            TOTELE, U_NLOC, NCOLCT, NCOLCMC, IGOT_CMC_PRECON
        LOGICAL, PARAMETER :: PIPES_1D=.TRUE.
        LOGICAL, intent( in ) :: got_free_surf, symmetric_P
        INTEGER, DIMENSION( : ), intent( in ) :: ndpset
        INTEGER, DIMENSION( : ), intent( in ) ::FINDC
        INTEGER, DIMENSION( : ), intent( in ) :: COLC
        REAL, DIMENSION( :, :, : ), intent( in ) :: INV_PIVIT_MAT
        INTEGER, DIMENSION( : ), intent( in ) ::  U_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
        INTEGER, DIMENSION( : ), intent( in ) :: COLCT
        REAL, DIMENSION( :, : ), intent( in ) :: DIAG_SCALE_PRES
        REAL, DIMENSION( :, :, : ), intent( in ) :: DIAG_SCALE_PRES_COUP, INV_B
        type(petsc_csr_matrix), intent(inout)::  CMC_petsc
        REAL, DIMENSION( :, :, : ), intent( inout ) :: CMC_PRECON
        REAL, DIMENSION( : ), intent( in ) :: MASS_MN_PRES
        REAL, DIMENSION( : ), intent( in ) :: MASS_SUF
        REAL, DIMENSION( : ), intent( in ) :: MASS_PIPE, MASS_CVFEM2PIPE, MASS_CVFEM2PIPE_TRUE
        INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
        INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
        REAL, DIMENSION( :, :, : ), intent( in ) :: C
        REAL, DIMENSION( :, :, : ), intent( in ) :: CT


        ! Local variables
        INTEGER, PARAMETER :: MX_NCOLOR = 1000
        REAL, PARAMETER :: INFINY = 1.0E+10
        LOGICAL :: UNDONE, LCOL
        logical, DIMENSION( : ), allocatable :: NEED_COLOR
        logical, DIMENSION( CV_NONODS ) :: to_color
        LOGICAL :: EXPLICIT_PIPES2
        REAL, DIMENSION( : ), allocatable :: COLOR_VEC
        REAL, DIMENSION( :, : ), allocatable :: CMC_COLOR_VEC, CMC_COLOR_VEC2, CMC_COLOR_VEC_PHASE, CMC_COLOR_VEC2_PHASE!, ld
        REAL, DIMENSION( : ), allocatable :: DU, DV, DW, DU_LONG
        REAL, DIMENSION( :, :, : ), allocatable :: CDP
        INTEGER :: NCOLOR, CV_NOD, CV_JNOD, COUNT, COUNT2, IDIM, IPHASE, CV_COLJ, U_JNOD, CV_JNOD2
        INTEGER :: I, ELE,u_inod,u_nod, ierr, IV_STAR, IV_FINI, N_IN_PRES, i_indx, j_indx, IPRES, JPRES
        REAL :: RSUM, RSUM_SUF, auxR(1)

        ALLOCATE( NEED_COLOR( CV_NONODS ) )
        ALLOCATE( COLOR_VEC( CV_NONODS ) )
        ALLOCATE( CDP( NDIM, NPHASE, U_NONODS ) )
        ALLOCATE( DU_LONG( U_NONODS * NDIM * NPHASE ) )
        ALLOCATE( DU( U_NONODS * NPHASE ) )
        ALLOCATE( DV( U_NONODS * NPHASE ) )
        ALLOCATE( DW( U_NONODS * NPHASE ) )
        ALLOCATE( CMC_COLOR_VEC( NPRES, CV_NONODS ) )
        ALLOCATE( CMC_COLOR_VEC2( NPRES, CV_NONODS ) )
        !      ALLOCATE( ld( NPRES, CV_NONODS ) ) ; ld = 0.0


        IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON = 0.0
        N_IN_PRES = NPHASE / NPRES
        EXPLICIT_PIPES2 = .true.

        NEED_COLOR = .TRUE.
        NCOLOR = 0
        UNDONE = .TRUE.

        Loop_while: DO WHILE ( UNDONE )

            NCOLOR = NCOLOR + 1 ! Determine what nodes can be coloured with the new color
            TO_COLOR = .FALSE.

            Loop_CVNOD: DO CV_NOD = 1, CV_NONODS
                IF ( NEED_COLOR( CV_NOD ) ) THEN

                    LCOL= .FALSE.
                    ! use a distance-2 colouring...
                    Loop_Row: DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                        CV_JNOD = COLCMC( COUNT )
                        IF( TO_COLOR( CV_JNOD ) ) THEN
                            LCOL=.TRUE.
                            EXIT
                        END IF
                        Loop_Row2: DO COUNT2 = FINDCMC( CV_JNOD ), FINDCMC( CV_JNOD + 1 ) - 1
                            IF ( TO_COLOR( COLCMC( COUNT2 ) ) )  THEN
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

            CALL C_MULT2( CDP, COLOR_VEC, CV_NONODS, U_NONODS, NDIM, NPHASE, &
                C, NCOLC, FINDC, COLC )
            ! DU_LONG = BLOCK_MAT * CDP

            CALL PHA_BLOCK_MAT_VEC( DU_LONG, INV_PIVIT_MAT, CDP, U_NONODS, NDIM, NPHASE, &
                TOTELE, U_NLOC, U_NDGLN )
            ! NB. P_RHS = CT * U + CV_RHS
            ! DU_LONG = CDP

            CALL ULONG_2_UVW( DU, DV, DW, DU_LONG, U_NONODS, NDIM, NPHASE )




            IF ( NPRES > 1 .AND. .NOT.EXPLICIT_PIPES2 ) THEN


                ALLOCATE( CMC_COLOR_VEC_PHASE( NPHASE, CV_NONODS ) )

                DO IPHASE = 1, NPHASE
                    IV_STAR = 1+(IPHASE-1)*U_NONODS*NDIM
                    IV_FINI = IPHASE*U_NONODS*NDIM
                    CALL CT_MULT( CMC_COLOR_VEC_PHASE(IPHASE,:), DU(IV_STAR:IV_FINI), &
                        DV(IV_STAR:IV_FINI), DW(IV_STAR:IV_FINI), CV_NONODS, U_NONODS, NDIM, 1, &
                        CT(:,IPHASE:IPHASE,:), NCOLCT, FINDCT, COLCT )
                END DO

                DO CV_NOD = 1, CV_NONODS
                    CMC_COLOR_VEC_PHASE(:,CV_NOD) = MATMUL( INV_B(:,:,CV_NOD), CMC_COLOR_VEC_PHASE(:,CV_NOD) )
                END DO

                DO CV_NOD = 1, CV_NONODS
                    DO IPRES = 1, NPRES
                        CMC_COLOR_VEC(IPRES,CV_NOD) = SUM(CMC_COLOR_VEC_PHASE(1+(IPRES-1)*N_IN_PRES:IPRES*N_IN_PRES,CV_NOD) )
                    END DO
                END DO


                IF ( IGOT_CMC_PRECON /= 0 ) THEN
                    ALLOCATE( CMC_COLOR_VEC2_PHASE( NPHASE, CV_NONODS ) )
                    DO IPHASE = 1, NPHASE
                        IV_STAR = 1+(IPHASE-1)*U_NONODS*NDIM
                        IV_FINI = IPHASE*U_NONODS*NDIM
                        CALL CT_MULT_WITH_C( CMC_COLOR_VEC2_PHASE(IPHASE,:), &
                            DU_LONG(IV_STAR:IV_FINI), &
                            CV_NONODS, U_NONODS, NDIM, 1,  &
                            C(:,IPHASE:IPHASE,:), NCOLC, FINDC, COLC )
                    END DO

                    DO CV_NOD = 1, CV_NONODS
                        CMC_COLOR_VEC2_PHASE(:,CV_NOD) = MATMUL( INV_B(:,:,CV_NOD), CMC_COLOR_VEC2_PHASE(:,CV_NOD) )
                    END DO

                    DO CV_NOD = 1, CV_NONODS
                        DO IPHASE = 1, NPHASE
                            CMC_COLOR_VEC2(IPHASE,CV_NOD) = SUM(CMC_COLOR_VEC2_PHASE(IPHASE:IPHASE,CV_NOD) )
                        END DO
                    END DO

                END IF

            ELSE

                DO IPRES = 1, NPRES
                    IV_STAR = 1+(IPRES-1)*N_IN_PRES*U_NONODS
                    IV_FINI = IPRES*N_IN_PRES*U_NONODS
                    CALL CT_MULT( CMC_COLOR_VEC(IPRES,:), DU(IV_STAR:IV_FINI), &
                        DV(IV_STAR:IV_FINI), DW(IV_STAR:IV_FINI), CV_NONODS, U_NONODS, NDIM, N_IN_PRES, &
                        CT(:,1+(IPRES-1)*N_IN_PRES:IPRES*N_IN_PRES,:), NCOLCT, FINDCT, COLCT )
                END DO
      
                IF ( IGOT_CMC_PRECON /= 0 ) THEN
                    DO IPRES = 1, NPRES
                        IV_STAR = 1+(IPRES-1)*N_IN_PRES*U_NONODS*NDIM
                        IV_FINI = IPRES*N_IN_PRES*U_NONODS*NDIM
                        CALL CT_MULT_WITH_C( CMC_COLOR_VEC2(IPRES,:), &
                            DU_LONG(IV_STAR:IV_FINI), CV_NONODS, U_NONODS, NDIM, N_IN_PRES, &
                            C(:,1+(IPRES-1)*N_IN_PRES:IPRES*N_IN_PRES,:), NCOLC, FINDC, COLC )
                    END DO
                END IF


            END IF




            ! Matrix vector involving the mass diagonal term
            DO CV_NOD = 1, CV_NONODS
                RSUM = 0.0
                RSUM_SUF = 0.0
                DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                    CV_JNOD = COLCMC( COUNT )
                    DO IPRES = 1, NPRES
                        IF((NPRES>1).AND.PIPES_1D) THEN
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

                        if ( got_free_surf ) then
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
                    DO IPRES = 1, NPRES
                        CMC_COLOR_VEC2( IPRES, CV_NOD ) = CMC_COLOR_VEC2( IPRES, CV_NOD ) &
                            + (DIAG_SCALE_PRES( IPRES, CV_NOD ) * RSUM + RSUM_SUF) * COLOR_VEC( CV_NOD )
                    END DO
                END IF
            END DO

            !Put into matrix CMC
            if ( .not.symmetric_P ) then

                ! original method
                DO CV_NOD = 1, CV_NONODS
                    DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                        CV_JNOD = COLCMC( COUNT )
                        DO IPRES = 1, NPRES
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
                DO CV_NOD = 1, CV_NONODS
                    DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                        CV_JNOD = COLCMC( COUNT )
                        DO IPRES = 1, NPRES
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


        ! the matrix coupling term place immediately into matrix and not through colouring...
        ! Matrix vector involving the mass diagonal term
        IF(NPRES > 1) THEN
            DO CV_NOD = 1, CV_NONODS
                RSUM = 0.0
                RSUM_SUF = 0.0
                DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                    CV_JNOD = COLCMC( COUNT )

                    DO IPRES = 1, NPRES

                        !if ( CV_NOD == CV_JNOD ) then
                         !  if ( MASS_CVFEM2PIPE( COUNT ) == 0.0 ) ld(ipres, cv_nod) = 1.0
                        ! end if

                        DO JPRES = 1, NPRES
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
        END IF ! ENDOF IF(NPRES > 1) THEN

        !If we have a reference node with pressure zero we impose that here.
        DO IPRES = 1, NPRES
            IF ( NDPSET(IPRES) > 0 ) THEN
                CV_NOD = NDPSET( IPRES )
                i_indx = CMC_petsc%row_numbering%gnn2unn( cv_nod,ipres)
                DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                    CV_JNOD = COLCMC( COUNT )
                    IF ( CV_JNOD /= CV_NOD ) THEN
                        do jpres = 1, npres
                            j_indx = CMC_petsc%column_numbering%gnn2unn( cv_jnod, jpres )
                            call MatSetValue(CMC_petsc%M, i_indx, j_indx, 0.0,INSERT_VALUES, ierr) ! not the diagonal
                            !CMC( COUNT ) = 0.0 ! not the diagonal
                            IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON( IPRES, JPRES, COUNT ) = 0.0
                            DO COUNT2 = FINDCMC( CV_JNOD ), FINDCMC( CV_JNOD + 1 ) - 1
                                CV_JNOD2 = COLCMC( COUNT2 )
                                IF ( CV_JNOD2 == CV_NOD ) then
                                    i_indx = CMC_petsc%row_numbering%gnn2unn( cv_jnod, ipres )
                                    j_indx = CMC_petsc%column_numbering%gnn2unn( CV_JNOD2, jpres )
                                    call MatSetValue(CMC_petsc%M, i_indx, j_indx, 0.0, INSERT_VALUES, ierr) ! not the diagonal
                                    IF ( IGOT_CMC_PRECON/=0 ) CMC_PRECON( ipres, jpres, COUNT2 ) = 0.0
                                END IF
                            END DO
                        end do ! do jpres=1,npres
                    END IF
                END DO
            END IF
        END DO


        DO CV_NOD = 1, CV_NONODS
            if ( mass_pipe(cv_nod) == 0.0 ) then
                CV_JNOD = CV_NOD
                DO IPRES = 2, NPRES
                    JPRES = IPRES
                    i_indx = CMC_petsc%row_numbering%gnn2unn( cv_nod, ipres )
                    j_indx = CMC_petsc%column_numbering%gnn2unn( CV_JNOD, jpres )
                    call MatSetValue(CMC_petsc%M, i_indx, j_indx, 1.0, INSERT_VALUES, ierr)
                END DO
            end if
        END DO




        if ( .false. ) then


            cv_nod = 226
            IPRES = 2
            i_indx = CMC_petsc%row_numbering%gnn2unn( cv_nod,ipres)
            DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                CV_JNOD = COLCMC( COUNT )
                jpres = ipres
                j_indx = CMC_petsc%column_numbering%gnn2unn( cv_jnod,jpres)
                IF ( CV_JNOD /= CV_NOD ) THEN
                    call MatSetValue(CMC_petsc%M, i_indx, j_indx, 0.0,INSERT_VALUES, ierr)
                ELSE
                    call MatSetValue(CMC_petsc%M, i_indx, j_indx, 1.0,INSERT_VALUES, ierr)
                END IF
            END DO




            cv_nod = 151
            IPRES = 2
            i_indx = CMC_petsc%row_numbering%gnn2unn( cv_nod,ipres)
            DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                CV_JNOD = COLCMC( COUNT )
                jpres = ipres
                j_indx = CMC_petsc%column_numbering%gnn2unn( cv_jnod,jpres)
                IF ( CV_JNOD /= CV_NOD ) THEN
                    call MatSetValue(CMC_petsc%M, i_indx, j_indx, 0.0,INSERT_VALUES, ierr)
                ELSE
                    call MatSetValue(CMC_petsc%M, i_indx, j_indx, 1.0,INSERT_VALUES, ierr)
                END IF
            END DO


        end if








        CMC_petsc%is_assembled = .false.
        call assemble( CMC_petsc )




        !DO CV_NOD = 1, CV_NONODS
        !   CV_JNOD = CV_NOD
        !   DO IPRES = 1, NPRES
        !      JPRES = IPRES
        !      i_indx = CMC_petsc%row_numbering%gnn2unn( cv_nod, ipres )
        !      j_indx = CMC_petsc%column_numbering%gnn2unn( CV_JNOD, jpres )
        !      call MatGetValues(cmc_petsc%M, 1, (/ i_indx /), 1, (/ j_indx /), auxR, ierr)
        !      print *, cv_nod, ipres, auxR, ld(ipres,cv_nod)
        !   END DO
        !END DO



        DEALLOCATE( NEED_COLOR )
        DEALLOCATE( COLOR_VEC )
        DEALLOCATE( CDP )
        DEALLOCATE( DU_LONG )
        DEALLOCATE( DU, DV, DW )
        DEALLOCATE( CMC_COLOR_VEC )
        DEALLOCATE( CMC_COLOR_VEC2 )
        !DEALLOCATE( ld )

        RETURN
    END SUBROUTINE COLOR_GET_CMC_PHA_SLOW



    !sprint_to_do!figure it out why are we not using the blas library
    SUBROUTINE PHA_BLOCK_INV( PIVIT_MAT, TOTELE, NBLOCK )
        implicit none
        INTEGER, intent( in ) :: TOTELE, NBLOCK
        REAL, DIMENSION( : , : , : ), intent( inout ), CONTIGUOUS :: PIVIT_MAT
        ! Local variables
        INTEGER :: ELE

        REAL, DIMENSION( NBLOCK , NBLOCK ) :: MAT, MAT2
        REAL, DIMENSION( NBLOCK ) :: X, B

        DO ELE = 1, TOTELE
            !         CALL MATINV( PIVIT_MAT( :, :, ele ), NBLOCK, NBLOCK )
            CALL MATINVold( PIVIT_MAT( :, :, ele ), NBLOCK, NBLOCK, MAT, MAT2, X, B )
        END DO

        RETURN
    END SUBROUTINE PHA_BLOCK_INV

    !sprint_to_do!see what we do with storage
    SUBROUTINE PHA_BLOCK_INV_plus_storage( PIVIT_MAT, TOTELE, &
        NBLOCK, Storage_state, StorName, indx)
          !Retrieves the inverse of the PIVIT_MAT fron the storage
        implicit none
        INTEGER, intent( in ) :: TOTELE, NBLOCK
        REAL, DIMENSION( : , : , : ), intent( inout ), pointer :: PIVIT_MAT
        type( state_type ), intent( inout ):: Storage_state
        character(len=*), intent(in) :: StorName
        integer, intent(inout) :: indx
        ! Local variables
        integer :: from, to
        type(mesh_type), pointer :: fl_mesh
        type(mesh_type) :: Auxmesh
        type(scalar_field), target :: targ_NX_ALL
        REAL, DIMENSION( :, :, : ), allocatable :: PIVIT_MAT2

        if (indx==0) then !The first time we need to introduce the targets in state
            if (has_scalar_field(Storage_state, trim(Storname))) then
                !If we are recalculating due to a mesh modification then
                !we return to the original situation
                call remove_scalar_field(Storage_state, trim(Storname))
            end if
            !Get mesh file just to be able to allocate the fields we want to store
            fl_mesh => extract_mesh( Storage_state, "FakeMesh" )
            Auxmesh = make_mesh(fl_mesh,name=trim(Storname))
            !The number of nodes I want does not coincide
            Auxmesh%nodes = NBLOCK * NBLOCK * TOTELE

            call allocate (Targ_NX_ALL, Auxmesh, trim(Storname))

            !Now we insert them in state and store the indexes
            call insert(Storage_state, Targ_NX_ALL, trim(Storname))
            !Store index
            indx = size(Storage_state%scalar_fields)

            call deallocate (Targ_NX_ALL)
            call deallocate (Auxmesh)

            !We have to calculate and store the inverse
            ALLOCATE( PIVIT_MAT2( NBLOCK, NBLOCK, TOTELE ))

            PIVIT_MAT2 = PIVIT_MAT!Very slow, but necessary because CX1 cannot reshape pointers...
            deallocate(PIVIT_MAT)!PIVIT_MAT comes already allocated
            nullify(PIVIT_MAT)
            call PHA_BLOCK_INV( PIVIT_MAT2, TOTELE, NBLOCK )

            !Store data
            from = 1; to = NBLOCK * NBLOCK * TOTELE
            Storage_state%scalar_fields(abs(indx))%ptr%val(from:to) =&
                reshape(PIVIT_MAT2,[NBLOCK * NBLOCK * TOTELE])
            deallocate(PIVIT_MAT2)
            indx = abs(indx)
        end if

        !Set the pointer to the  solution
        from = 1; to = NBLOCK * NBLOCK * TOTELE
        call reshape_vector2pointer(Storage_state%scalar_fields(abs(indx))%ptr%val(from:to),&
            PIVIT_MAT, NBLOCK, NBLOCK, TOTELE)


    END SUBROUTINE PHA_BLOCK_INV_plus_storage


    SUBROUTINE PHA_BLOCK_MAT_VEC_old( U, BLOCK_MAT, CDP, U_NONODS, NDIM, NPHASE, &
        TOTELE, U_NLOC, U_NDGLN )
        implicit none
        ! U = BLOCK_MAT * CDP
        INTEGER, intent( in )  :: U_NONODS, NDIM, NPHASE, TOTELE, U_NLOC
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

        Loop_Elements: DO ELE = 1, TOTELE

            U_NOD => U_NDGLN( ( ELE - 1 ) * U_NLOC +1 : ELE * U_NLOC )

            LCDP = RESHAPE( CDP( :, :, U_NOD ) , (/ N /) )
            CALL DGEMV( 'N', N, N, 1.0d0, BLOCK_MAT( : , : , ELE ), N, LCDP, 1, 0.0d0, LU, 1 )

            DO U_ILOC = 1, U_NLOC
                U( 1+(U_NOD(U_ILOC)-1)*NDIM*NPHASE : U_NOD(U_ILOC)*NDIM*NPHASE ) = [LU(:,:,U_ILOC)]
            END DO

        END DO Loop_Elements

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
        REAL, DIMENSION( :, :, : ), intent( in ) :: CDP
        ! Local
        INTEGER :: ELE, U_ILOC, U_INOD, IDIM, IPHASE, I, U_JLOC, U_JNOD, JDIM, JPHASE, J, II, JJ

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

        RETURN


    END SUBROUTINE PHA_BLOCK_MAT_VEC






    SUBROUTINE PHA_BLOCK_MAT_VEC2( U, BLOCK_MAT, CDP, U_NONODS, NDIM, NPHASE, &
        TOTELE, U_NLOC, U_NDGLN )
        implicit none
        ! U = BLOCK_MAT * CDP
        INTEGER, intent( in )  :: U_NONODS, NDIM, NPHASE, TOTELE, U_NLOC
        INTEGER, DIMENSION( : ), intent( in ), target ::  U_NDGLN
        REAL, DIMENSION( :, :, : ), intent( inout ) :: U
        REAL, DIMENSION( :, :, : ), intent( in ), target :: BLOCK_MAT
        REAL, DIMENSION( :, :, : ), intent( in ) :: CDP
        ! Local
        INTEGER :: ELE, N
        !      INTEGER, DIMENSION( : ), pointer :: U_NOD
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

        Loop_Elements: DO ELE = 1, TOTELE

            !         U_NOD => U_NDGLN( ( ELE - 1 ) * U_NLOC +1 : ELE * U_NLOC )
            U_NOD = U_NDGLN( ( ELE - 1 ) * U_NLOC +1 : ELE * U_NLOC )

            LCDP = RESHAPE( CDP( :, :, U_NOD ) , (/ N /) )
            CALL DGEMV( 'N', N, N, 1.0d0, BLOCK_MAT( : , : , ELE ), N, LCDP, 1, 0.0d0, LU, 1 )
            U( :, :, U_NOD ) = RESHAPE( LU, (/ NDIM, NPHASE, U_NLOC/) )
      
        END DO Loop_Elements

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
        INTEGER :: ELE, U_ILOC, U_INOD, IDIM, IPHASE, I, U_JLOC, U_JNOD, JDIM, JPHASE, J, II, JJ, IVEC

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


 

    SUBROUTINE PHA_BLOCK_MAT_VEC_MANY( U, BLOCK_MAT, CDP, U_NONODS, NDIM, NPHASE, NBLOCK, &
        TOTELE, U_NLOC, U_NDGLN )
        implicit none
        ! U = BLOCK_MAT * CDP
        INTEGER, intent( in )  :: U_NONODS, NDIM, NPHASE, TOTELE, U_NLOC, NBLOCK
        INTEGER, DIMENSION( : ), intent( in ) ::  U_NDGLN
        REAL, DIMENSION( :, :, :, : ), intent( inout ), contiguous, target :: U
        REAL, DIMENSION( :, : , : ), intent( in ), contiguous :: BLOCK_MAT
        REAL, DIMENSION( :, : , :, : ), intent( in ), contiguous, target :: CDP
        ! Local
        INTEGER :: ELE, U_ILOC, U_INOD, IDIM, IPHASE, I, J, U_JLOC, U_JNOD, JDIM, JPHASE, II, JJ, IORIG

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

    !sprint_to_do!have a look at this one
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

        CV_RHS = 0.0

        DO CV_INOD = 1, CV_NONODS

            DO COUNT = FINDCT( CV_INOD ), FINDCT( CV_INOD + 1 ) - 1, 1
                U_JNOD = COLCT( COUNT )

                DO IPHASE = 1, NPHASE
                    J = U_JNOD + ( IPHASE - 1 ) * U_NONODS

                    CV_RHS( CV_INOD ) = CV_RHS( CV_INOD ) + CT( 1, IPHASE, COUNT ) * U( J )
                    IF( NDIM >= 2 ) CV_RHS( CV_INOD ) = CV_RHS( CV_INOD ) + CT( 2, IPHASE, COUNT ) * V( J )
                    IF( NDIM >= 3 ) CV_RHS( CV_INOD ) = CV_RHS( CV_INOD ) + CT( 3, IPHASE, COUNT ) * W( J )
                END DO

            END DO

        END DO

        RETURN

    END SUBROUTINE CT_MULT





    !sprint_to_do !try to use the many version only
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
        REAL, DIMENSION( NBLOCK, CV_NONODS ), intent( inout) :: CV_RHS
        REAL, DIMENSION( NBLOCK, NDIM, NPHASE, U_NONODS ), intent( in ) :: U
        INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
        INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
        REAL, DIMENSION( NDIM, NPHASE, NCOLCT ), intent( in ) :: CT

        ! Local variables
        INTEGER :: CV_INOD, COUNT, U_JNOD, IPHASE, J, IVEC, IDIM

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
                !            forall ( IVEC = 1 : NBLOCK, IPHASE = 1 : NPHASE,IDIM =1:NDIM)
                !                  CV_RHS( IVEC, CV_INOD ) = CV_RHS( IVEC, CV_INOD )+ U( IVEC, IDIM, IPHASE, U_JNOD ) * CT( IDIM, IPHASE, COUNT  )
                !            end forall
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
        REAL, DIMENSION( NBLOCK, NDIM, NPHASE, U_NONODS ), intent( inout ) :: CDP
        REAL, DIMENSION( NBLOCK, CV_NONODS ), intent( in )  :: DP
        REAL, DIMENSION( NDIM, NPHASE, NCOLC ), intent( in ) :: C
        INTEGER, DIMENSION( U_NONODS + 1 ), intent( in ) ::FINDC
        INTEGER, DIMENSION( NCOLC ), intent( in ) :: COLC
        ! Local variables
        INTEGER :: U_INOD, COUNT, P_JNOD, IPHASE, I1, IDIM, COUNT_DIM_PHA

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
        INTEGER :: U_INOD, COUNT, P_JNOD, IPHASE, IDIM, COUNT_DIM_PHA, DIM_PHA

        CDP = 0.0

        DO U_INOD = 1, U_NONODS
            DO COUNT = FINDC( U_INOD ), FINDC( U_INOD + 1 ) - 1
                P_JNOD = COLC( COUNT )
                CDP( :, :, U_INOD ) = CDP( :, :, U_INOD ) + C( :, :, COUNT ) * DP( P_JNOD )
            END DO
        END DO

        RETURN

    END SUBROUTINE C_MULT2

    SUBROUTINE CT_MULT_WITH_C( DP, U_LONG, CV_NONODS, U_NONODS, NDIM, NPHASE, &
        C, NCOLC, FINDC, COLC )
        implicit none
        ! DP = (C)^T U_LONG
        INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLC
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



    !sprint_to_do!can this go?
    SUBROUTINE CT_MULT_WITH_C3( DP, U_ALL, CV_NONODS, U_NONODS, NDIM, NPHASE, &
        C, NCOLC, FINDC, COLC )
        implicit none
        ! DP = (C)^T U_ALL
        INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLC
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


    SUBROUTINE CT_MULT_WITH_C_MANY( DP, U_LONG, CV_NONODS, U_NONODS, NDIM, NPHASE, NBLOCK, &
        C, NCOLC, FINDC, COLC )
        implicit none
        ! DP = (C)^T U_LONG
        INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLC, NBLOCK
        REAL, DIMENSION( :, :, :, : ), intent( in ) :: U_LONG
        REAL, DIMENSION( :, : ), intent( inout )  :: DP
        REAL, DIMENSION( :, :, : ), intent( in ) :: C
        INTEGER, DIMENSION( : ), intent( in ) ::FINDC
        INTEGER, DIMENSION( : ), intent( in ) :: COLC
        ! Local variables
        INTEGER :: U_INOD, COUNT, P_JNOD, IPHASE, I1, IDIM, COUNT_DIM_PHA

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
            IF( NDIM >= 2 ) &
                V( 1 + ( IPHASE - 1 ) * U_NONODS : U_NONODS + ( IPHASE - 1 ) * U_NONODS ) = &
                UP( 1 + U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS : 2 * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS )
            IF( NDIM >= 3 ) &
                W( 1 + ( IPHASE - 1 ) * U_NONODS : U_NONODS + ( IPHASE - 1 ) * U_NONODS ) = &
                UP( 1 + 2 * U_NONODS + ( IPHASE - 1) * NDIM * U_NONODS : 3 * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS )
        END DO
        RETURN
    END SUBROUTINE ULONG_2_UVW



    subroutine posinmat( posmat, globi, globj, &
        nonods, findrm, colm, ncolm )
        ! Find position in matrix POSMAT which has column GLOBJ
        implicit none
        integer, intent( inout ) :: posmat
        integer, intent( in ) :: globi, globj, nonods, ncolm
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

        integer :: node, jphase, count, node_count, nphase


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


