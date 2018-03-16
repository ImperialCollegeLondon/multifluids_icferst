
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


module sparsity_1D
    use fldebug

    use spud
    use global_parameters, only: option_path_len, is_porous_media
    use multi_data_types
contains

    subroutine resize(A,n,copy)
        integer, dimension(:), pointer, intent(inout) :: A
        integer n
        logical, intent(in), optional  :: copy
        integer, dimension(:), pointer :: temp
        logical :: lcopy


        if(present(copy)) then
            lcopy=copy
        else
            lcopy=.true.
        end if

        allocate(temp(n))
        if (lcopy) temp=A(1:n)
        deallocate(a)
        a=>temp
    end subroutine resize

    SUBROUTINE DEF_SPAR( SEMI_BAND_WID, NONODS, MX_NCOLC, NCOLC, &
        CENTC, FINDC, COLC)
        ! define sparsity...
        ! SEMI_BAND_WID is the semi band width.
        IMPLICIT NONE
        INTEGER, intent( in ) :: SEMI_BAND_WID, NONODS, MX_NCOLC
        INTEGER, intent( inout) :: NCOLC
        INTEGER, DIMENSION( NONODS ), intent( inout) ::  CENTC
        INTEGER, DIMENSION( NONODS + 1 ), intent( inout) :: FINDC
        INTEGER, DIMENSION( MX_NCOLC ), intent( inout) :: COLC

        ! Local variables
        INTEGER :: NOD, II, COL, COUNT

        ewrite(3,*) 'In DEF_SPAR'

        COUNT = 0
        COL = 0
        loop_nods: DO NOD = 1, NONODS

            FINDC( NOD ) = COUNT + 1

            DO II = -SEMI_BAND_WID, SEMI_BAND_WID, 1

                COL = NOD + II

                IF( ( COL >= 1 ) .AND. ( COL <= NONODS ) ) THEN
                    COUNT = COUNT + 1
                    COLC( COUNT ) = COL
                    IF( COL == NOD) CENTC( NOD ) = COUNT
                END IF

            END DO

        END DO loop_nods

        FINDC( NONODS + 1 ) = COUNT + 1

        NCOLC = COUNT

        ewrite(3,*) 'Leaving DEF_SPAR'

        RETURN
    END SUBROUTINE DEF_SPAR

    SUBROUTINE DEF_SPAR_CT_DG( CV_NONODS, MX_NCT, NCT, FINDCT, COLCT, &
        TOTELE, CV_NLOC, U_NLOC, U_NDGLN, CV_NDGLN)
        ! define sparsity...
        ! SEMI_BAND_WID is the semi band width.
        IMPLICIT NONE
        INTEGER, intent( in ) :: CV_NONODS, MX_NCT
        INTEGER, intent( inout ) :: NCT
        INTEGER, DIMENSION( CV_NONODS + 1 ), intent (inout ) :: FINDCT
        INTEGER, DIMENSION( MX_NCT ), intent( inout ) :: COLCT
        INTEGER, intent( in ) :: TOTELE, CV_NLOC, U_NLOC
        INTEGER, DIMENSION ( U_NLOC * TOTELE ), intent( in ) :: U_NDGLN
        integer, dimension (cv_nloc * totele ), intent( in ) :: cv_ndgln
        ! Local variables...
        INTEGER :: CV_NOD, U_NOD, JLOC, COUNT, ELE, ELE1, ELE2, CV_NODI, CV_ILOC, count2, rep
        integer, dimension(cv_nonods) :: cv_ndgln_small
        logical :: repeated, finished_colct

        ewrite(3,*) 'In DEF_SPAR_CT_DG'

        COUNT = 2
        !! Get condensed form of cv_ndgln, ie without any repeats
        !! This can then be used as the index for the loop below
        cv_ndgln_small=0
        rep = 1
        cv_ndgln_small(1)=cv_ndgln(1)
        do cv_nod=2,size(cv_ndgln)
            repeated = .false.
            do cv_nodi=1,cv_nod-rep
                if (cv_ndgln(cv_nod)==cv_ndgln_small(cv_nodi)) then
                    repeated=.true.
                    rep=rep+1
                end if
            end do
            if (.not.repeated) then
                cv_ndgln_small(count) = cv_ndgln(cv_nod)
                count = count+1
            end if
        end do

        finished_colct=.false.
        COUNT = 0
        count2 = 1
        !ewrite(3,*)'cv1:', size(cv_ndgln), cv_ndgln
        !ewrite(3,*)'cv2:', size(cv_ndgln_small), cv_ndgln_small
        !ewrite(3,*)'cvnonods:', cv_nonods
        !stop 23

        IF(CV_NONODS /= CV_NLOC*TOTELE ) THEN
            ! Have a cty CV_NOD
            do while (.not.finished_colct)
                loop_cvnod: DO CV_NOD = 1, CV_NONODS

                    cv_nodi = cv_ndgln_small(cv_nod)
                    if (cv_nodi==count2) then

                        FINDCT( cv_nodi ) = COUNT + 1
                        ELE1 = 1 + ( CV_NOD - 2 ) / ( CV_NLOC - 1 )
                        ELE2 = 1 + ( CV_NOD - 1 ) / ( CV_NLOC - 1 )
                        !ewrite(3,*)'findct:', cv_nod, cv_nodi, findct( cv_nodi )

                        loop_elements: DO ELE = MAX( 1 , ELE1 ), MIN( TOTELE , ELE2 ), 1

                            DO JLOC = 1 ,U_NLOC
                                U_NOD = U_NDGLN(( ELE - 1 ) * U_NLOC + JLOC )
                                COUNT = COUNT + 1
                                !ewrite(3,*) u_nod
                                COLCT( COUNT ) = U_NOD
                               !ewrite(3,*)'colct:', ele, count, colct(count)
                            END DO

                        END DO loop_elements

                    end if

                END DO loop_cvnod
                if (count2==cv_nonods) finished_colct = .true.
                count2 = count2+1
            end do
        ELSE
            ! Have a discontinuous CV_NOD
            loop_elements2: DO ELE = 1,TOTELE
                loop_cviloc: DO CV_ILOC = 1, CV_NLOC

                    CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC

                    FINDCT( CV_NOD ) = COUNT + 1

                    ! IF(U_ELE_TYPE==2) THEN
                    IF((CV_ILOC==1).AND.(ELE /= 1)) THEN
                        ELE2=ELE-1
                        JLOC=U_NLOC
                        U_NOD = U_NDGLN(( ELE2 - 1 ) * U_NLOC + JLOC )
                        COUNT = COUNT + 1
                        COLCT( COUNT ) = U_NOD
                    ENDIF

                    DO JLOC = 1 ,U_NLOC
                        U_NOD = U_NDGLN(( ELE - 1 ) * U_NLOC + JLOC )
                        COUNT = COUNT + 1
                        COLCT( COUNT ) = U_NOD
                    END DO

                    IF((CV_ILOC==CV_NLOC).AND.(ELE /= TOTELE)) THEN
                        ELE2=ELE+1
                        JLOC=1
                        U_NOD = U_NDGLN(( ELE2 - 1 ) * U_NLOC + JLOC )
                        COUNT = COUNT + 1
                        COLCT( COUNT ) = U_NOD
                    ENDIF

                END DO loop_cviloc

            END DO loop_elements2
        ENDIF

        FINDCT( CV_NONODS + 1) = COUNT + 1
        NCT = COUNT

        IF(NCT > MX_NCT ) THEN
            EWRITE(3,*)'MX_NCT is not long enough NCT,MX_NCT:',NCT,MX_NCT
            stop 2721
        ENDIF

        !DO CV_NODI=1,CV_NONODS
        !   EWRITE(3,*)'CV_NODI=',CV_NODI
        !   EWRITE(3,*)(COLCT( COUNT ),COUNT = FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1 )
        !END DO

        !ewrite(3,*)'findct', (findct(cv_nod), cv_nod=1,cv_nonods)
        !ewrite(3,*)'colct', (colct(cv_nod), cv_nod=1,nct)

        ewrite(3,*) 'Leaving DEF_SPAR_CT_DG'

        RETURN
    END SUBROUTINE DEF_SPAR_CT_DG

end module sparsity_1D

module sparsity_ND
    use fldebug
    use sparsity_1D
    use shape_functions_prototype
    use shape_functions_Linear_Quadratic
    use cv_advection
    use multi_data_types
    use multi_tools
contains

    subroutine getfinele( totele, nloc, snloc, nonods, ndglno, mx_nface_p1, &
        mxnele, ncolele, finele, colele, midele )
        ! This sub caluculates COLELE the element connectivitiy list
        ! in order of faces.
        implicit none
        integer, intent( in ) :: totele, nloc, snloc, nonods
        integer, dimension( totele * nloc ), intent( in ) :: ndglno
        integer, intent( in ) :: mx_nface_p1, mxnele
        integer, intent( inout ) :: ncolele
        integer, dimension( mxnele ), intent( inout ) :: colele
        integer, dimension( totele + 1 ), intent( inout ) :: finele
        integer, dimension( totele ), intent( inout ) :: midele
        ! Local variables
        integer :: ele, iloc, jloc, iloc2, nod, inod, jnod, count, ele2, i, hit, &
            iface, itemp, count2
        logical :: found
        integer, allocatable, dimension( : ) :: fintran, coltran, icount

        allocate( fintran( nonods + 1 ))
        allocate( coltran( max( totele, nonods ) * mx_nface_p1 ))
        allocate( icount( max( nonods, totele) ))

        icount = 0
        do ele = 1, totele
            do iloc = 1, nloc
                nod = ndglno( ( ele - 1 ) * nloc + iloc )
                icount( nod ) = icount( nod ) + 1
            end do
        end do

        fintran = 0
        fintran( 1 ) = 1
        do nod = 1, nonods
            fintran( nod + 1 ) = fintran( nod ) + icount( nod )
        end do

        icount = 0
        coltran = 0
        do ele = 1, totele
            do iloc = 1, nloc
                nod = ndglno( ( ele - 1 ) * nloc + iloc )
                !ewrite(3,*)'nod, filtran, icount, ele:', nod, fintran( nod ), ele, totele, nonods
                coltran( fintran( nod ) + icount( nod )) = ele
                icount( nod ) = icount( nod ) + 1
            end do
        end do
        !ewrite(3,*)'coltran:', coltran( 1: max(totele, nonods ) * mx_nface_p1 )
        !ewrite(3,*)'fintran:', fintran( 1: nonods + 1 )
        !ewrite(3,*)'X_NDGLN:', ndglno( 1: totele*nloc )

        icount = 0 ; colele = 0 ; ncolele = 0
        Loop_Elements1: do ele = 1, totele

            Loop_Iloc: do iloc = 1, nloc

                nod = ndglno( ( ele - 1 ) * nloc + iloc )
                Loop_Count1: do count = fintran( nod ), fintran( nod + 1 ) - 1, 1

                    ele2 = coltran( count )
                    found = .false. ! Add ELE2 into list FINELE and COLELE
                    do i = 1, icount( ele )
                        if( colele( ( ele - 1 ) * mx_nface_p1 + i ) == ele2 ) found = .true.
                    end do

                    Conditional_Found: if ( .not. found ) then ! Do elements ELE and ELE2 share at least 3 nodes?

                        hit = 0
                        do iloc2 = 1, nloc
                            inod = ndglno( ( ele - 1 ) * nloc + iloc2 )
                            do jloc = 1, nloc
                                jnod = ndglno( ( ele2 - 1 ) * nloc + jloc )
                                if ( inod == jnod ) hit = hit + 1
                            end do
                        end do
                        if ( hit >= snloc ) then
                            icount( ele ) = icount( ele ) + 1
                            colele( ( ele - 1 ) * mx_nface_p1 + icount( ele )) = ele2
                            ncolele = ncolele + 1
                        end if

                    end if Conditional_Found

                end do Loop_Count1

            end do Loop_Iloc

        end do Loop_Elements1

        finele( 1 ) = 1
        do ele = 1, totele
            finele( ele + 1 ) = finele( ele ) + icount( ele )
        end do

        ! order elements in increasing order...
        count = 0
        Loop_Elements2: do ele = 1, totele
            ! Shorten COLELE then perform a bubble sort to get the ordering right for.
            do iface = 1, mx_nface_p1
                if ( colele( ( ele - 1 ) * mx_nface_p1 + iface ) /= 0 ) then
                    count = count + 1
                    colele( count ) = colele( ( ele - 1 ) * mx_nface_p1 + iface )
                end if
            end do
        end do Loop_Elements2

        Loop_BubbleSort: do ele = 1, totele
            do count = finele( ele ) , finele( ele + 1 ) - 2
                do count2 = finele( ele ) , finele( ele + 1 ) - 1, 1
                    if ( colele( count ) > colele( count + 1 )) then ! swop over
                        itemp = colele( count + 1 )
                        colele( count + 1 ) = colele( count )
                        colele( count ) = itemp
                    end if
                end do
            end do
        end do Loop_BubbleSort

        ! Calculate midele:
        do ele = 1, totele
            do count = finele( ele ) , finele( ele + 1 ) - 1
                if(colele(count)==ele) midele( ele ) = count
            end do
        end do

        deallocate( fintran )
        deallocate( coltran )
        deallocate( icount )

        return
    end subroutine getfinele

    subroutine exten_sparse_multi_phase_old( nonods, mxnele, finm, colm, &
        nphase, npha_nonods, ncolm_pha, &
        finm_pha, colm_pha, midm_pha )
        ! Extend the sparsity to a multiphase sparsity
        implicit none
        integer, intent( in ) :: nonods, mxnele
        integer, dimension( nonods + 1 ), intent( in ) :: finm
        integer, dimension( mxnele ), intent( in ) :: colm
        integer, intent( in ) :: nphase, npha_nonods, ncolm_pha
        integer, dimension( npha_nonods + 1 ), intent( inout ) :: finm_pha
        integer, dimension( ncolm_pha ), intent( inout ) :: colm_pha
        integer, dimension( npha_nonods ), intent( inout ) :: midm_pha
        ! Local variables
        integer :: count, count2, iphase, jphase, nod

        ewrite(3,*) 'In exten_sparse_multi_phase_old subrt.'

        count2 = 0
        Loop_Phase1: do iphase = 1, nphase
            Loop_CVNODS: do nod = 1, nonods
                Loop_Phase2: do jphase = 1, nphase
                    if( jphase == 1 ) &
                        finm_pha( ( iphase - 1 ) * nonods + nod ) = count2 + 1
                    Conditional_Phases: if( iphase == jphase ) then
                        do count = finm( nod ), finm( nod + 1 ) - 1
                            count2 = count2 + 1
                            colm_pha( count2 ) = colm( count ) + ( jphase - 1 ) * nonods
                            if( colm( count ) == nod ) &
                                midm_pha( nod + ( jphase - 1 ) * nonods ) = count2
                        end do
                    else
                        count2 = count2 + 1
                        colm_pha( count2 ) = nod + ( jphase - 1 ) * nonods
                    end if Conditional_Phases
                end do Loop_Phase2
            end do Loop_CVNODS
        end do Loop_Phase1

        finm_pha( nphase * nonods + 1 ) = count2 + 1
        !ncolm_pha = count2
        if(count2.ne.ncolm_pha) then
            ewrite(3,*) 'not correct length count2,ncolm_pha:',count2,ncolm_pha
            stop 28219
        end if

        ewrite(3,*) 'Leaving exten_sparse_multi_phase_old subrt.'
        return
    end subroutine exten_sparse_multi_phase_old

    subroutine exten_sparse_multi_phase( nonods, mxnele, finm, colm, &
        nphase, npha_nonods, ncolm_pha, &
        finm_pha, colm_pha, midm_pha )
        ! Extend the sparsity to a multiphase sparsity
        implicit none
        integer, intent( in ) :: nonods, mxnele
        integer, dimension( nonods + 1 ), intent( in ) :: finm
        integer, dimension( mxnele ), intent( in ) :: colm
        integer, intent( in ) :: nphase, npha_nonods, ncolm_pha
        integer, dimension( npha_nonods + 1 ), intent( inout ) :: finm_pha
        integer, dimension( ncolm_pha ), intent( inout ) :: colm_pha
        integer, dimension( npha_nonods ), intent( inout ) :: midm_pha
        ! Local variables
        integer :: count, count2, iphase, jphase, nod

        ewrite(3,*) 'In exten_sparse_multi_phase subrt.'

        count2 = 0
        Loop_CVNODS: do nod = 1, nonods
            Loop_Phase1: do iphase = 1, nphase
                finm_pha( ( nod - 1 ) * nphase + iphase ) = count2 + 1
                do count = finm( nod ), finm( nod + 1 ) - 1
                    if (colm(count) .ne. nod) then
                        count2 = count2 + 1
                        colm_pha( count2 ) = iphase + ( colm( count ) - 1) * nphase
                    else
                        do jphase = 1, nphase
                            count2 = count2 + 1
                            colm_pha( count2 ) = jphase + ( colm( count ) - 1) * nphase
                            if (jphase==iphase) then
                                midm_pha( iphase + (nod-1) * nphase )=count2
                            end if
                        end do
                    end if
                end do
            end do Loop_Phase1
        end do Loop_CVNODS

        finm_pha( nphase * nonods + 1 ) = count2 + 1
        !ncolm_pha = count2
        if(count2.ne.ncolm_pha) then
            ewrite(3,*) 'not correct length count2,ncolm_pha:',count2,ncolm_pha
            stop 2821
        end if

        !ewrite(3,*) 'colm_pha--',colm_pha(1:ncolm_pha)

        ewrite(3,*) 'Leaving exten_sparse_multi_phase subrt.'
        return
    end subroutine exten_sparse_multi_phase


    subroutine form_dgm_pha_sparsity( totele, nphase, u_nloc, u_pha_nonods, &
        ndim, mx_ncoldgm_pha, ncoldgm_pha, &
        coldgm_pha, findgm_pha, middgm_pha, &
        finele, colele, ncolele )
        ! Form the sparsity of the phase coupled DG discretised matrix
        ! from the element-wise multi-phase sparsity matrix.
        implicit none
        integer, intent( in ) :: totele, nphase, u_nloc, u_pha_nonods, &
            mx_ncoldgm_pha, ndim, ncolele
        integer, intent( inout ) :: ncoldgm_pha
        integer, dimension( mx_ncoldgm_pha ), intent( inout ) :: coldgm_pha
        integer, dimension( u_pha_nonods + 1 ), intent( inout ) :: findgm_pha
        integer, dimension( u_pha_nonods ), intent( inout ) :: middgm_pha
        integer, dimension( totele + 1 ), intent( in ) :: finele
        integer, dimension( ncolele ), intent( in ) :: colele

        ! Local variables
        integer :: count, count2, iloc, jloc, irow, jrow
        integer :: ele, ele2, idim, jdim, iphase, jphase, u_nonods
        logical :: new_sparcity

        new_sparcity=.true.


        ewrite(3,*) 'In form_dgm_pha_sparsity subrt.'

        u_nonods = u_pha_nonods / ( nphase * ndim )


        if ( new_sparcity ) then

            count2 = 0
            do ele = 1, totele
                do iloc = 1, u_nloc
                    do iphase = 1, nphase
                        do idim = 1, ndim
                            !irow = ( ele - 1 ) * u_nloc + iloc  + ( idim - 1 ) * u_nonods + (iphase-1)*u_nonods*ndim
                            irow = ( ele - 1 ) * u_nloc*ndim*nphase + (iloc-1)*ndim*nphase  + ( iphase - 1 ) * ndim + idim
                            !ewrite(3,*)'irow, ele, u_nloc, iloc, idim, u_nonods, iphase, count2:', &
                            !irow, ele, u_nloc, iloc, idim, u_nonods, iphase, count2+1
                            findgm_pha( irow ) = count2 + 1
                            do count = finele( ele ), finele( ele + 1 ) - 1
                                ele2 = colele( count )
                                do jloc = 1, u_nloc
                                    do jphase = 1, nphase
                                        do jdim = 1, ndim
                                            !jrow = ( ele2 - 1 ) * u_nloc + jloc  + ( jdim - 1 ) * u_nonods + (jphase-1)*u_nonods*ndim
                                            jrow = ( ele2 - 1) * u_nloc*ndim*nphase + (jloc-1)*ndim*nphase  + ( jphase - 1 ) * ndim + jdim
                                            count2 = count2 + 1
                                            coldgm_pha( count2 ) = jrow
                                            if( irow == jrow ) middgm_pha( irow ) = count2
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            findgm_pha( u_pha_nonods + 1 ) = count2 + 1
            ncoldgm_pha = count2

        else

            count2 = 0
            Loop_Phase1: do iphase = 1, nphase
                Loop_Dim1: do idim = 1, ndim
                    Loop_Element: do ele = 1, totele
                        Loop_Loc1: do iloc = 1, u_nloc
                            irow = ( ele - 1 ) * u_nloc + iloc  + ( idim - 1 ) * u_nonods + (iphase-1)*u_nonods*ndim
                            !ewrite(3,*)'irow, ele, u_nloc, iloc, idim, u_nonods, iphase, count2:', &
                            !     irow, ele, u_nloc, iloc, idim, u_nonods, iphase, count2+1
                            findgm_pha( irow ) = count2 + 1
                            Loop_Phase2: do jphase = 1, nphase
                                Loop_Dim2: do jdim = 1, ndim
                                    Loop_Count: do count = finele( ele ), finele( ele + 1 ) - 1
                                        ele2 = colele( count )
                                        Loop_Loc2: do jloc = 1, u_nloc
                                            jrow = ( ele2 - 1 ) * u_nloc + jloc  + ( jdim - 1 ) * u_nonods + (jphase-1)*u_nonods*ndim
                                            count2 = count2 + 1
                                            coldgm_pha( count2 ) = jrow
                                            if( irow == jrow ) middgm_pha( irow ) = count2
                                        end do Loop_Loc2
                                    end do Loop_Count
                                end do Loop_Dim2
                            end do Loop_Phase2
                        end do Loop_Loc1
                    end do Loop_Element
                end do Loop_Dim1
            end do Loop_Phase1
            findgm_pha( u_pha_nonods + 1 ) = count2 + 1
            ncoldgm_pha = count2

        end if


        ! perform a bubble sort to order the row in ioncreasing order
        do irow = 1, u_pha_nonods
            call quicksort( coldgm_pha( findgm_pha( irow ) : findgm_pha( irow + 1 ) - 1 ) , 1)
            do count = findgm_pha( irow ), findgm_pha( irow + 1 ) - 1
                jrow = coldgm_pha( count )
                if( irow == jrow ) middgm_pha( irow ) = count
            end do
        end do

        if( ncoldgm_pha > mx_ncoldgm_pha ) &
            FLAbort(" Incorrect number of dimension of sparsity matrix - ncoldgm_pha ")

        ewrite(3,*) 'Leaving form_dgm_pha_sparsity subrt. '

        return
    end subroutine form_dgm_pha_sparsity


    subroutine pousinmc2( totele, nloc1, nonods2, nloc2, &
        nimem, ndglno1, ndglno2, &
        lencolm, findrm, colm, centrm )
        implicit none
        integer, intent( in ) :: totele, nloc1, nonods2, nloc2, nimem
        integer, dimension( totele * nloc1 ), intent( in ) :: ndglno1
        integer, dimension( totele * nloc2 ), intent( in ) :: ndglno2
        integer, intent( inout ) :: lencolm
        integer, dimension( nonods2 + 1 ), intent( inout ) :: findrm
        integer, dimension( nimem ), intent( inout ) :: colm
        integer, dimension( nonods2 ), intent( inout ) :: centrm

        ! Local Variables
        integer :: ele, globi, globj, loci, locj, i, irow, ptr

        ! Defining derived data types for the linked list
        type node
            integer :: id                 ! id number of node
            type( node ), pointer :: next ! next node
        end type node

        type row
            type( node ), pointer :: row ! recursive data type
        end type row

        type( row ), dimension( : ), allocatable :: matrix
        type( node ), pointer :: list, current, next

        allocate( matrix( nonods2 ))
        do i = 1, nonods2
            allocate( list )
            list % id = -1
            nullify( list % next )
            matrix( i ) % row => list
            nullify( list )
        end do

        Loop_Elements1: do ele = 1, totele

            Loop_LocI: do loci = 1, nloc2

                globi = ndglno2( ( ele - 1 ) * nloc2 + loci )
                list => matrix( globi ) % row

                Loop_LocJ: do locj = 1, nloc1

                    globj = ndglno1( ( ele - 1 ) * nloc1 + locj )

                    if ( list % id == -1 ) then ! Check if the list is initalised
                        list % id = globj
                        cycle
                    end if

                    !      if ( list % id == -1 ) then ! Check if the list is initalised
                    !         list % id = globj
                    !         cycle
                    !      end if

                    Conditional1: if ( globj < list % id ) then ! Insert at start of list
                        allocate( current )
                        current % id = globj
                        current % next => list
                        matrix( globi ) % row => current
                        list => matrix( globi ) % row

                    else ! Conditional1
                        current => list
                        Loop_While1: do while ( associated( current ))

                            if ( globj == current % id ) then ! Already have this node
                                exit

                            elseif ( .not. associated( current % next )) then  ! End of list - insert this node
                                allocate( current % next )
                                nullify( current % next % next )
                                current % next % id = globj
                                exit

                            elseif ( globj < current % next % id ) then ! Insert new node here
                                allocate( next )
                                next % id = globj
                                next % next => current % next
                                current % next => next
                                exit

                            end if
                            current => current % next

                        end do Loop_While1

                    end if Conditional1

                end do Loop_LocJ

            end do Loop_LocI

        end do Loop_Elements1

        ! From matrix write COLM, FINDRM and CENTRM
        ! linked list as we go
        ptr = 1
        !ewrite(3,*),'nonods2=',nonods2

        Loop_Irow: do irow = 1, nonods2
            findrm( irow ) = ptr
            centrm( irow ) = -1

            current => matrix( irow ) % row
            Loop_While2: do while ( associated( current ))
                assert( ptr <= nimem )
                colm( ptr ) = current % id
                if ( current % id == -1 ) then
                    ewrite(0,*) "ERROR: POSINM() seriously unhappy with node", IROW
                    FLAbort( "ERROR: Mesh contains nodes that are not associated with any elements." )
                end if
                if ( current % id == irow ) then
                    centrm( irow ) = ptr
                endif
                next => current % next
                deallocate( current )
                current => next
                ptr = ptr + 1
            end do Loop_While2

        end do Loop_Irow

        lencolm = ptr - 1
        findrm( nonods2 + 1 ) = lencolm + 1

        if(nimem.lt.lencolm) then
            print *,'nimem not long enough nimem,lencolm:',nimem,lencolm
            stop 822
        endif

        deallocate( matrix )

        return
    end subroutine pousinmc2


    subroutine conv_ct2c( cv_nonods, nct, findct, colct, u_nonods, &
        mx_nc, findc, colc )
        implicit none
        integer, intent( in ) :: cv_nonods, nct
        integer, dimension( cv_nonods + 1 ), intent( in ) :: findct
        integer, dimension( nct ), intent( in ) :: colct
        integer, intent( in ) :: u_nonods, mx_nc
        integer, dimension( u_nonods + 1 ), intent( inout ) :: findc
        integer, dimension( mx_nc ), intent( inout ) :: colc
        ! Local variables
        integer :: col, u_nod, cv_nod, count, count2
        integer, dimension( : ), allocatable :: in_row_c

        ewrite(3,*) 'In conv_ct2c subrt.'
        ewrite(3,*) ' max( colct ):', maxval( colct )

        ! No. of non-zero's in row of C matrix.
        allocate( in_row_c( u_nonods ) )
        in_row_c = 0
        do count = 1, nct
            col = colct( count )
            in_row_c( col ) = in_row_c( col ) + 1
        end do

        count = 0
        do u_nod = 1, u_nonods
            findc( u_nod ) = count + 1
            count = count + in_row_c( u_nod )
        end do
        findc( u_nonods + 1 ) = count + 1
        in_row_c( 1 : u_nonods ) = 0

        Loop_CVNOD: do cv_nod = 1, cv_nonods
            Loop_Count: do count = findct( cv_nod ), findct( cv_nod + 1 ) - 1
                col = colct( count )
                in_row_c( col ) = in_row_c( col ) + 1
                count2 = findc( col ) + in_row_c( col ) - 1
                colc( count2 ) = cv_nod
            end do Loop_Count
        end do Loop_CVNOD

        deallocate( in_row_c )

        ewrite(3,*) 'Leaving conv_ct2c subrt.'

        return
    end subroutine conv_ct2c


    subroutine poscmc( totele, nonods, nimem, nct, &
        findct, colct, &
        ncmc, fincmc, colcmc, midcmc, noinod, presym )
        ! This subroutine forms the matrix operating on the pressure vector.
        ! It is found from C1T ML C1 + C2T ML C2
        ! In the first part of COLCMC contains the pressure nodes surrounding
        ! a given node.
        implicit none
        integer, intent ( in ) :: totele, nonods, nimem, nct
        integer, dimension( totele + 1 ), intent( in ) :: findct
        integer, dimension( nct ), intent( in ) :: colct
        integer, intent( inout ) :: ncmc
        integer, dimension( totele + 1 ), intent( inout ) :: fincmc
        integer, dimension(:), pointer, intent( inout ) :: colcmc
        integer, dimension( totele ), intent( inout ) :: midcmc
        integer, dimension( nonods ), intent( inout ) ::  noinod
        logical, intent( inout ) :: presym

        ! Local variables
        integer :: count, globi, globj, i, irow, inod, jrow, ptr

        ! Defining derived data types for the linked list
        type node
            integer :: id                 ! id number of node
            type( node ), pointer :: next ! next node, recursive data type
        end type node

        type row
            type( node ), pointer :: row
        end type row

        type( row ), dimension( : ), allocatable :: matrix, matrix2
        type( node ), pointer :: list, current, current2, next

        allocate( matrix( nonods ))
        do i = 1, nonods
            allocate( list )
            list % id = -1
            nullify( list % next )
            matrix( i ) % row => list
            nullify( list )
        end do

        noinod = 0

        Loop_Row1: do irow = 1, totele ! Given a vel node find the pressure nodes surrounding it
            globj = irow
            Loop_Count1: do count = findct( irow ), findct( irow + 1 ) - 1
                globi = colct( count )
                list => matrix( globi ) % row
                ! ewrite(3,*)'list%id:',globj,globi,list % id

                if ( list % id == -1 ) then ! Check if the list is initalised
                    list % id = globj
                    cycle
                end if

                Conditional1: if ( globj < list % id ) then ! Insert at start of list
                    allocate( current )
                    current % id = globj
                    current % next => list
                    matrix( globj ) % row => current
                    list => matrix( globj ) % row
                else
                    current => list
                    !ewrite(3,*)'-->', globj, current % id, current % next % id
                    Loop_While1: do while( associated( current ))
                        if ( globj == current % id ) then ! Already have this node
                            exit
                        elseif ( .not. associated( current % next )) then ! End of list - insert this node
                            allocate( current % next )
                            nullify( current % next % next )
                            current % next % id = globj
                            exit
                        elseif( globj < current % next % id ) then ! Insert new node here
                            allocate( next )
                            next % id = globj
                            next % next => current % next
                            current % next => next
                            exit
                        end if
                        current => current % next
                    end do Loop_While1
                end if Conditional1

                noinod( globi ) = noinod( globi ) + 1

            end do Loop_Count1

        end do Loop_Row1

        allocate( matrix2( totele )) ! Initalise the linked lists
        do i = 1, totele
            allocate( list )
            list % id = -1
            nullify( list % next )
            matrix2( i ) % row => list
            nullify( list )
        end do

        Loop_Row2: do irow = 1, totele
            ! Find the pressure nodes surrounding pressure node IROW
            Loop_Count2: do count = findct( irow ), findct( irow + 1 ) - 1
                inod = colct( count )

                ! Find the pressure nodes surrounding node INOD
                ! these will be connected to pressure node IROW.
                current => matrix( inod ) % row
                Loop_While2: do while( associated( current ))
                    jrow = current % id

                    Conditional2: if (( .not. presym ) .or. ( jrow >= irow )) then
                        list => matrix2( irow ) % row

                        if ( list % id == -1 ) then ! Check if the list is initialised
                            list % id = jrow
                            cycle
                        end if

                        Conditional3: if ( jrow < list % id ) then ! Insert at start of list
                            allocate( current2 )
                            current2 % id = jrow
                            current2 % next => list
                            matrix2( irow ) % row => current2
                            list => matrix2( irow ) % row

                        else
                            current2 => list
                            Loop_While3: do while( associated( current2 ))

                                Conditional4: if ( jrow == current2 % id ) then ! Already have this node
                                    exit
                                elseif( .not. associated( current2 % next )) then ! End of list - insert this node
                                    allocate( current2 % next )
                                    nullify(  current2 % next % next )
                                    current2 % next % id = jrow
                                    exit
                                elseif( jrow < current2 % next % id ) then ! Insert new node here
                                    allocate( next )
                                    next % id = jrow
                                    next % next => current2 % next
                                    current2 % next => next
                                    exit
                                end if Conditional4
                                current2 => current2 % next

                            end do Loop_While3

                        end if Conditional3

                    end if Conditional2
                    current => current % next

                end do Loop_While2

            end do Loop_Count2

        end do Loop_Row2


        do irow = 1, nonods ! Delete Matrix
            current => matrix( irow ) % row
            do while ( associated( current ))
                next => current % next
                deallocate( current )
                current => next
            end do
        end do
        deallocate( matrix )

        ptr = 1
        do irow = 1, totele
            current => matrix2( irow ) % row
            do while ( associated( current ))
                next => current % next
                current => next
                ptr = ptr + 1
            end do
        end do


        ncmc =  ptr - 1
        call resize(colcmc,ncmc,copy=.false.)

        ptr = 1



        Loop_Row3: do irow = 1, totele ! From matrix write COLCMC, FINCMC and MIDCMC
            midcmc( irow ) = -1
            fincmc( irow ) = ptr
            current => matrix2( irow ) % row ! Warning: overwriten

            Loop_While4: do while ( associated( current ))

                if ( ptr > ncmc ) then
                    ewrite( -1, * ) 'ncmc, nimem, ptr: ',ncmc, nimem, ptr
                    ewrite( -1, * ) 'totele, irow: ',totele, irow
                    FLAbort( "Integer memory too small" )
                end if

                colcmc( ptr ) = current % id
                if( current % id == irow ) then
                    midcmc( irow ) = ptr
                end if

                next => current % next
                deallocate( current )
                current => next
                ptr = ptr + 1

            end do Loop_While4
        end do Loop_Row3

        fincmc( totele + 1 ) = ncmc + 1

        return
    end subroutine poscmc

    subroutine CV_Neighboor_Sparsity( Mdims, cv_ele_type, &
         cv_ndgln, x_ndgln, &
         ncolele, finele, colele, &
         ncolm, mxnacv_loc, findm, colm, &
         ncolacv_loc, finacv_loc, colacv_loc, midacv_loc )
      implicit none
      type(multi_dimensions), intent(in) :: Mdims
      integer, intent( in ) :: cv_ele_type
      integer, dimension( Mdims%totele * Mdims%cv_nloc ), intent( in ) :: cv_ndgln
      integer, dimension( Mdims%totele * Mdims%x_nloc ), intent( in ) :: x_ndgln
      integer, intent( in ) :: ncolele
      integer, dimension( Mdims%totele + 1 ), intent( in ) :: finele
      integer, dimension( ncolele ), intent( in ) :: colele
      integer, intent( in ) :: ncolm, mxnacv_loc
      integer, dimension( Mdims%cv_nonods + 1 ), intent( in ) :: findm
      integer, dimension( ncolm ), intent( in ) :: colm
      integer, intent( inout ) :: ncolacv_loc
      integer, dimension( Mdims%cv_nonods + 1 ), intent( inout ) :: finacv_loc
      integer, dimension( mxnacv_loc ), intent( inout ) :: colacv_loc
      integer, dimension( Mdims%cv_nonods ), intent( inout ) :: midacv_loc
      ! Local variables
      type( multi_GI_dimensions ) :: GIdims
      type (multi_shape_funs) :: CV_funs
      logical, dimension( ncolm ) :: found
      integer ::  ele, ele2, cv_iloc, cv_jloc, cv_nodi, cv_nodj, &
           cv_nodj2, gcount, gcount2, gi, &
           x_nodi, x_nodi2, jcount, cv_nodi2, cv_iloc2
      ! Computing Gauss points and array containing node points on neighboors elements
      call retrieve_ngi( GIdims, Mdims, cv_ele_type, .false. )
      call allocate_multi_shape_funs( CV_funs, Mdims, GIdims )
      call cv_fem_shape_funs( CV_funs, Mdims, GIdims, cv_ele_type, quad_over_whole_ele = .false. )

      ! Allocating space
      found = .false.
      ! set the diagonal to true...
      do cv_nodi = 1, Mdims%cv_nonods
         do gcount = findm( cv_nodi ), findm( cv_nodi + 1 ) - 1
            if( colm(gcount) == cv_nodi ) found( gcount ) = .true.
         end do
      end do
      ! now the off diagonal terms...
      Loop_Elements_1: do ele = 1, Mdims%totele
         Loop_CVILOC_1: do cv_iloc = 1, Mdims%cv_nloc
            cv_nodi = cv_ndgln( ( ele - 1 ) * Mdims%cv_nloc + cv_iloc )
            Loop_GI_1: do gi = 1, GIdims%scvngi
               cv_jloc = CV_funs%cv_neiloc( cv_iloc, gi )
               if ( cv_jloc  > 0 ) then
                  cv_nodj = cv_ndgln( ( ele - 1 ) * Mdims%cv_nloc + cv_jloc )
                  Loop_GCOUNT_1: do gcount = findm( cv_nodi ), findm( cv_nodi + 1 ) - 1
                     cv_nodj2 = colm( gcount )
                     if ( cv_nodj == cv_nodj2 ) found( gcount ) = .true.
                  end do Loop_GCOUNT_1
               endif
            end do Loop_GI_1
         end do Loop_CVILOC_1
      end do Loop_Elements_1
      ! for discontinuous elements...
      if( Mdims%totele * Mdims%cv_nloc == Mdims%cv_nonods ) then
         Loop_Elements_4: do ele = 1, Mdims%totele
            Loop_FaceCount: do jcount = finele( ele ), finele( ele + 1 ) - 1
               ele2 = colele( jcount )
               if( ( ele2 > 0 ) .and. ( ele2 /= ele ) ) then
                  Loop_CVILOC_5: do cv_iloc = 1, Mdims%cv_nloc ! Loop over nodes of the elements
                     cv_nodi = cv_ndgln( ( ele - 1 ) * Mdims%cv_nloc + cv_iloc )
                     x_nodi = x_ndgln( ( ele - 1 ) * Mdims%cv_nloc + cv_iloc )
                     Loop_CVILOC_6: do cv_iloc2 = 1, Mdims%cv_nloc ! Loop over nodes of the elements
                        cv_nodi2 = cv_ndgln( ( ele2 - 1 ) * Mdims%cv_nloc + cv_iloc2 )
                        x_nodi2 = x_ndgln( ( ele2 - 1 ) * Mdims%cv_nloc + cv_iloc2 )
                        if(x_nodi==x_nodi2) then
                           do gcount = findm( cv_nodi ), findm( cv_nodi + 1 ) - 1
                              cv_nodj2 = colm( gcount )
                              if ( cv_nodi2 == cv_nodj2 ) found( gcount ) = .true.
                           end do
                        endif
                     end do Loop_CVILOC_6
                  end do Loop_CVILOC_5
               endif
            end do Loop_FaceCount
         end do Loop_Elements_4
      endif
      gcount2 = 0 ! Now reducing the size of the stencil

      do cv_nodi = 1, Mdims%cv_nonods
         finacv_loc( cv_nodi ) = gcount2 + 1
         do gcount = findm( cv_nodi ), findm( cv_nodi + 1 ) - 1
            if( found( gcount ) ) then
               gcount2 = gcount2 + 1
               colacv_loc( gcount2 ) = colm( gcount )
               if( colacv_loc( gcount2 ) == cv_nodi ) midacv_loc( cv_nodi ) = gcount2
            end if
         end do
      end do
      ncolacv_loc = gcount2
      finacv_loc( Mdims%cv_nonods + 1 ) = gcount2 + 1
      !Deallocate shape functions
      call deallocate_multi_shape_funs(CV_funs)
      return
    end subroutine CV_Neighboor_Sparsity


end module sparsity_ND



module spact

    use fldebug
    use spud
    use global_parameters, only: option_path_len
    use futils, only: int2str
    use sparsity_1D
    use sparsity_ND
    use shape_functions_prototype
    use Copy_Outof_State
    use multi_data_types
    use multi_tools
contains


    subroutine Defining_MaxLengths_for_Sparsity_Matrices( ndim, nphase, totele, u_nloc, cv_nloc, ph_nloc, cv_nonods, &
        mx_nface_p1, mxnele, mx_nct, mx_nc, mx_ncolcmc, mx_ncoldgm_pha, &
        mx_ncolacv, mx_ncolm, mx_ncolph )
        implicit none
        integer, intent( in ) :: ndim, nphase, totele, u_nloc, cv_nloc, ph_nloc, cv_nonods
        integer, intent( inout ) :: mx_nface_p1, mxnele, mx_nct, mx_nc, mx_ncolcmc, mx_ncoldgm_pha, &
            mx_ncolacv, mx_ncolm, mx_ncolph

        ewrite(3,*)'In Defining_Lengths_for_Sparsity_Matrices'

        mx_nface_p1 = 2 * ndim + 1
        mxnele = mx_nface_p1 * totele

        mx_nct = totele * u_nloc * cv_nloc * ndim * nphase
        mx_ncolm = abs(mxnele * cv_nloc * cv_nloc * nphase * nphase)!the abs() is to solve an odd bug apearing in gfortran

        if( cv_nonods == cv_nloc * totele ) then ! Discontinuous in Pressure
            !   mx_nct = mx_nct * mx_nface_p1 * 10
            !   mx_ncolm = mx_ncolm * 5
            mx_nct = mx_nct * mx_nface_p1
            mx_ncolm = mx_ncolm * 3
        end if
        mx_nc = mx_nct

        !!$ Assuming the DG representation requires the more storage space
        mx_ncolcmc = mx_nface_p1 **3 * cv_nloc * cv_nloc * totele
        if(is_porous_media) then
            mx_ncoldgm_pha = 1
        else
            mx_ncoldgm_pha = ( mxnele + totele ) * ( u_nloc * ndim * nphase) **2    ! for overlapping method =1

                ! for overlapping method =1
        endif

        mx_ncolacv = 3 * mx_nface_p1 * cv_nonods * nphase + cv_nonods * ( nphase - 1 ) * nphase

        mx_ncolph = mxnele * 4  * 6 + ph_nloc

        return
    end subroutine Defining_MaxLengths_for_Sparsity_Matrices

    subroutine Get_Sparsity_Patterns( state, Mdims, Mspars, ndgln, Mdisopt, mx_ncolacv, &
                mx_ncoldgm_pha, mx_nct,mx_nc, mx_ncolcmc, mx_ncolm, mx_ncolph, mx_nface_p1 )
        !!$ Allocate and obtain the sparsity patterns of the two types of matricies for
        !!$ (momentum + cty) and for energy
        implicit none
        type( state_type ), dimension( : ), intent( inout ) :: state
        type(multi_dimensions), intent(inout) :: Mdims
        type (multi_sparsities), intent(inout) :: Mspars
        type(multi_ndgln), intent(in) :: ndgln
        type (multi_discretization_opts) :: Mdisopt
        integer, intent( in ) :: mx_ncolacv, mx_ncoldgm_pha, &
            mx_nct, mx_nc, mx_ncolcmc, mx_ncolm, mx_ncolph, mx_nface_p1
        !!$ Local variables
        integer, dimension( : ), pointer :: ph_ndgln, colele_pha, finele_pha, midele_pha, centct, dummyvec
        integer :: mx_ncolsmall_acv, count, cv_inod, mx_ncolele_pha, nsmall_acv, nsmall_acv2, stat
        logical :: presym
        type(csr_sparsity), pointer :: sparsity
        type(mesh_type), pointer :: element_mesh, ph_mesh
        ewrite(3,*)'In Get_Sparsity_Patterns'
        !Check if sparsities have been associated (allocated), if not, allocate

        call deallocate_multi_sparsities(Mspars)
        call allocate_multi_sparsities(Mspars, Mdims, mx_ncolacv, &
             mx_ncoldgm_pha, mx_nct, mx_nc, mx_ncolm, mx_ncolph)

        !-
        !- Computing sparsity for element connectivity
        !-
        element_mesh=> extract_mesh(state(1),"P0DG")
        allocate(sparsity)
        sparsity = make_sparsity_compactdgdouble(element_mesh,&
            name="ElementConnectivity")
        call insert(state(1),sparsity,name="ElementConnectivity")
        call deallocate(sparsity)
        deallocate(sparsity)
        sparsity=> extract_csr_sparsity(state(1),name="ElementConnectivity")
        Mspars%ELE%fin => sparsity%findrm
        Mspars%ELE%mid => sparsity%centrm
        Mspars%ELE%col => sparsity%colm
        Mspars%ELE%ncol=size(Mspars%ELE%col)
        if(.not.(is_porous_media)) then
            !-
            !- Computing sparsity for force balance
            !-
            mx_ncolele_pha = Mdims%nphase * Mspars%ELE%ncol + ( Mdims%nphase - 1 ) * Mdims%nphase * Mdims%totele
            allocate( colele_pha( mx_ncolele_pha ) )
            allocate( finele_pha( Mdims%totele * Mdims%nphase + 1 ) )
            allocate( midele_pha( Mdims%totele * Mdims%nphase ) )
            colele_pha = 0 ; finele_pha = 0 ; midele_pha = 0
            call exten_sparse_multi_phase_old( Mdims%totele, Mspars%ELE%ncol, Mspars%ELE%fin, Mspars%ELE%col, &
                Mdims%nphase, Mdims%totele * Mdims%nphase, mx_ncolele_pha, &
                finele_pha, colele_pha, midele_pha )
            Mspars%DGM_PHA%fin = 0 ; Mspars%DGM_PHA%col = 0 ; Mspars%DGM_PHA%mid = 0
            call form_dgm_pha_sparsity( Mdims%totele, Mdims%nphase, Mdims%u_nloc, Mdims%nphase * Mdims%u_nonods * Mdims%ndim, &
                Mdims%ndim, mx_ncoldgm_pha, Mspars%DGM_PHA%ncol, &
                Mspars%DGM_PHA%col, Mspars%DGM_PHA%fin, Mspars%DGM_PHA%mid, &
                Mspars%ELE%fin, Mspars%ELE%col, Mspars%ELE%ncol )
            ! dealocate colele_pha...
            deallocate( colele_pha ) ; deallocate( finele_pha ) ; deallocate( midele_pha )
        else
            Mspars%DGM_PHA%ncol=0
        end if
        call resize(Mspars%DGM_PHA%col,Mspars%DGM_PHA%ncol)
        !-
        !- momentum and continuity eqns
        !-
        allocate( centct( Mdims%cv_nonods ) )
        Mspars%CT%fin = 0 ; Mspars%CT%col = 0 ; centct = 0
        Conditional_Dimensional_2: if ( ( Mdims%ndim == 1 ) .and. .false. ) then
            call def_spar_ct_dg( Mdims%cv_nonods, mx_nct, Mspars%CT%ncol, Mspars%CT%fin, Mspars%CT%col, &
                Mdims%totele, Mdims%cv_nloc, Mdims%u_nloc, ndgln%u, ndgln%cv )
        else
            if( Mdims%cv_nonods == Mdims%x_nonods ) then ! a continuous pressure mesh
                call pousinmc2( Mdims%totele, Mdims%u_nloc, Mdims%cv_nonods, Mdims%cv_nloc, &
                    mx_nct, ndgln%u, ndgln%cv, &
                    Mspars%CT%ncol, Mspars%CT%fin, Mspars%CT%col, centct )
            else ! use Mspars%ELE%fin to determine Mspars%CT%col
                call CT_DG_Sparsity( mx_nface_p1,  &
                    Mdims%totele, Mdims%cv_nloc, Mdims%u_nloc,  &
                    Mdims%cv_nonods, &
                    ndgln%cv, ndgln%u,  &
                    Mspars%ELE%ncol, Mspars%ELE%fin, Mspars%ELE%col, &
                    mx_nct, Mspars%CT%ncol, Mspars%CT%fin, Mspars%CT%col )
            endif
            call resize(Mspars%CT%col,Mspars%CT%ncol)
        end if Conditional_Dimensional_2
        Mspars%C%ncol = Mspars%CT%ncol
        !-
        !- Convert CT sparsity to C sparsity
        !-
        call conv_ct2c( Mdims%cv_nonods, Mspars%CT%ncol, Mspars%CT%fin, Mspars%CT%col, Mdims%u_nonods, &
            mx_nc, Mspars%C%fin, Mspars%C%col )
        call resize(Mspars%C%col,Mspars%C%ncol)
        !-
        !- Computing sparsity for pressure matrix of projection method
        !-
        Conditional_Dimensional_3: if ( ( Mdims%ndim == 1 ) .and. .false. ) then
            Conditional_ContinuousPressure_3: if ( Mdims%cv_nonods /= Mdims%totele * Mdims%cv_nloc ) then
                call def_spar( Mdims%cv_nloc - 1, Mdims%cv_nonods, mx_ncolcmc, Mspars%CMC%ncol, &
                    Mspars%CMC%mid, Mspars%CMC%fin, Mspars%CMC%col )
            else ! Discontinuous pressure mesh
                call def_spar( Mdims%cv_nloc + 2, Mdims%cv_nonods, mx_ncolcmc, Mspars%CMC%ncol, &
                    Mspars%CMC%mid, Mspars%CMC%fin, Mspars%CMC%col )
            end if Conditional_ContinuousPressure_3
        else
            allocate( dummyvec( Mdims%u_nonods ))
            dummyvec = 0
            presym = .false.
            call poscmc( Mdims%cv_nonods, Mdims%u_nonods, mx_ncolcmc, Mspars%CT%ncol, &
                Mspars%CT%fin, Mspars%CT%col, &
                Mspars%CMC%ncol, Mspars%CMC%fin, Mspars%CMC%col, Mspars%CMC%mid, dummyvec, presym )
            deallocate( dummyvec )
        end if Conditional_Dimensional_3
        if( Mspars%CMC%ncol<1 ) FLAbort("Incorrect number of dimension of CMC sparsity matrix")
        call resize(Mspars%CMC%col,Mspars%CMC%ncol)
        !-
        !- Computing sparsity CV-FEM
        !-
        Mspars%M%fin = 0 ; Mspars%M%col = 0 ; Mspars%M%mid = 0
        Conditional_Dimensional_4: if ( ( Mdims%ndim == 1 ) .and. .false. ) then
            call def_spar( Mdims%cv_nloc - 1, Mdims%cv_nonods, mx_ncolm, Mspars%M%ncol, &
                Mspars%M%mid, Mspars%M%fin, Mspars%M%col )
        else
            if(Mdims%cv_nonods==Mdims%x_nonods) then ! a continuous pressure mesh
                call pousinmc2( Mdims%totele, Mdims%cv_nloc, Mdims%cv_nonods, Mdims%cv_nloc, mx_ncolm, ndgln%cv, ndgln%cv, &
                    Mspars%M%ncol, Mspars%M%fin, Mspars%M%col, Mspars%M%mid )
            else ! a DG pressure field mesh
                call CT_DG_Sparsity( mx_nface_p1, &
                    Mdims%totele, Mdims%cv_nloc, Mdims%cv_nloc, &
                    Mdims%cv_nonods, &
                    ndgln%cv, ndgln%cv, &
                    Mspars%ELE%ncol, Mspars%ELE%fin, Mspars%ELE%col, &
                    mx_ncolm, Mspars%M%ncol, Mspars%M%fin, Mspars%M%col )
                ! Determining Mspars%M%mid:
                do cv_inod = 1, Mdims%cv_nonods
                    do count = Mspars%M%fin( cv_inod ), Mspars%M%fin( cv_inod + 1 ) - 1
                        if( Mspars%M%col( count ) == cv_inod ) Mspars%M%mid( cv_inod ) = count
                    end do
                end do
            end if
        end if Conditional_Dimensional_4
        call resize(Mspars%M%col,Mspars%M%ncol)
        !-
        !- Computing sparsity for CV multiphase eqns (e.g. vol frac, temp)
        !-
        mx_ncolsmall_acv=mx_ncolacv/Mdims%nphase
        allocate( Mspars%small_acv%mid( Mdims%cv_nonods ) )
        allocate( Mspars%small_acv%fin( Mdims%cv_nonods + 1 ) )
        allocate( Mspars%small_acv%col( mx_ncolsmall_acv ) )
        Mspars%small_acv%ncol = mx_ncolsmall_acv
        Mspars%small_acv%mid = 0 ; Mspars%small_acv%fin = 0 ; Mspars%small_acv%col = 0
        call CV_Neighboor_Sparsity( Mdims, Mdisopt%cv_ele_type, &
            ndgln%cv, ndgln%x, &
            Mspars%ELE%ncol, Mspars%ELE%fin, Mspars%ELE%col, &
            Mspars%M%ncol, mx_ncolsmall_acv, Mspars%M%fin, Mspars%M%col, &
            nsmall_acv, Mspars%small_acv%fin, Mspars%small_acv%col, Mspars%small_acv%mid )
        nsmall_acv2 = nsmall_acv
        call resize(Mspars%small_acv%col,nsmall_acv)
        Mspars%ACV%ncol =  Mdims%nphase * nsmall_acv + ( Mdims%nphase - 1 ) * Mdims%nphase * Mdims%cv_nonods
        nsmall_acv = Mspars%ACV%ncol
        Mspars%small_acv%ncol = Mspars%ACV%ncol!<== Not sure about this (but it is working...)
        Mspars%ACV%fin = 0 ; Mspars%ACV%col = 0 ; Mspars%ACV%mid = 0
        call exten_sparse_multi_phase( Mdims%cv_nonods, nsmall_acv2, Mspars%small_acv%fin, Mspars%small_acv%col, &
            Mdims%nphase, Mdims%nphase * Mdims%cv_nonods, Mspars%ACV%ncol, &
            Mspars%ACV%fin, Mspars%ACV%col, Mspars%ACV%mid)
        call resize(Mspars%ACV%col,Mspars%ACV%ncol)
        !-
        !- Computing sparsity for ph (hydrostatic pressure)
        !-
        ph_mesh => extract_mesh( state( 1 ), "ph", stat )
        if ( stat == 0 ) then
            Mdims%ph_nonods = node_count( ph_mesh )
            Mdims%ph_nloc = ele_loc( ph_mesh, 1 )
            ph_ndgln => get_ndglno( ph_mesh )
            Mspars%ph%fin = 0 ; Mspars%ph%col = 0
            if ( Mdims%cv_nonods == Mdims%x_nonods ) then ! a continuous pressure mesh   ! BUG HERE!!!
               call pousinmc2( Mdims%totele, Mdims%ph_nloc, Mdims%ph_nonods, Mdims%ph_nloc, &
                    mx_ncolph, ph_ndgln, ph_ndgln, Mspars%ph%ncol, Mspars%ph%fin, Mspars%ph%col, Mspars%ph%mid )
            else ! a DG pressure field mesh
                call CT_DG_Sparsity( mx_nface_p1, Mdims%totele, Mdims%ph_nloc, Mdims%ph_nloc, &
                    Mdims%ph_nonods, ph_ndgln, ph_ndgln, Mspars%ELE%ncol, Mspars%ELE%fin, Mspars%ELE%col, &
                    mx_ncolph, Mspars%ph%ncol, Mspars%ph%fin, Mspars%ph%col )
            end if
            call resize( Mspars%ph%col, Mspars%ph%ncol )
        end if
        !-
        !- Deallocating temporary arrays
        !-
        deallocate( centct )
        return
    end subroutine Get_Sparsity_Patterns

    subroutine CT_DG_Sparsity( mx_nface_p1, &
        totele, cv_nloc, u_nloc, &
        cv_nonods, &
        cv_ndgln, u_ndgln, &
        ncolele, finele, colele, &
        mx_nct, nct, findct, colct )
        implicit none
        integer, intent( in ) :: mx_nface_p1, totele, cv_nloc, &
            u_nloc, cv_nonods
        integer, dimension( totele * cv_nloc ), intent( in ) :: cv_ndgln
        integer, dimension( totele * u_nloc ), intent( in ) :: u_ndgln
        integer, intent( in ) :: ncolele
        integer, dimension( totele + 1 ), intent( in ) :: finele
        integer, dimension( ncolele ), intent( in ) :: colele
        integer, intent( in ) :: mx_nct
        integer, intent( inout ) :: nct
        integer, dimension( cv_nonods + 1 ), intent( inout ) :: findct
        integer, dimension( mx_nct ), intent( inout ) :: colct
        ! Local variables
        integer :: ele, ele2, cv_iloc, cv_nodi, count, count2, u_jloc, u_nodj, &
            gcount2, gcount, FACE_COUNT
        integer, dimension( : ), allocatable :: find_ct_temp, no_in_row

        allocate( find_ct_temp( cv_nonods + 1 ) ) ; find_ct_temp = 0
        allocate( no_in_row( cv_nonods ) ) ; no_in_row = 0

        find_ct_temp(1) = 1
        do cv_nodi = 1, cv_nonods
            find_ct_temp( cv_nodi + 1 ) = find_ct_temp( cv_nodi ) &
                + ( mx_nface_p1 + 1 ) * u_nloc
        end do

        colct = 0
        Loop_Elements_4: do ele = 1, totele
            Loop_CVILOC_5: do cv_iloc = 1, cv_nloc ! Loop over nodes of the elements
                cv_nodi = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )
                count2 = 0

                do FACE_COUNT = FINELE( ELE ), FINELE( ELE + 1 ) - 1

                    ele2 = COLELE( FACE_COUNT )
                    if( ele2 > 0 ) then
                        Loop_CVILOC_6: do u_jloc = 1, u_nloc ! Loop over nodes of the elements
                            u_nodj = u_ndgln( ( ele2 - 1 ) * u_nloc + u_jloc )
                            count2 = count2 + 1
                            count = find_ct_temp( cv_nodi ) - 1 + count2
                            colct( count ) = u_nodj
                            no_in_row( cv_nodi ) = no_in_row( cv_nodi ) + 1
                        end do Loop_CVILOC_6
                    end if
                end do

            end do Loop_CVILOC_5
        end do Loop_Elements_4

        !stop 22

        ! shrink the sparcity up a little now...
        gcount2 = 0 ! Now reducing the size of the stencil
        do cv_nodi = 1, cv_nonods
            findct( cv_nodi ) = gcount2 + 1
            do gcount = find_ct_temp( cv_nodi ), find_ct_temp( cv_nodi + 1 ) - 1
                if( colct( gcount ) /= 0 ) then
                    gcount2 = gcount2 + 1
                    colct( gcount2 ) = colct( gcount )
                end if
            end do
        end do
        nct = gcount2
        findct( cv_nonods + 1 ) = gcount2 + 1

        if(mx_nct.lt.nct) then
            print *,'mx_nct not long enough mx_nct,nct:',mx_nct,nct
            stop 221
        endif

        ! sort colct in increasing order
        do cv_nodi = 1, cv_nonods
            call quicksort( colct( findct( cv_nodi ) : findct( cv_nodi + 1 ) -1 ) , 1)
        end do

        return
    end subroutine CT_DG_Sparsity

    subroutine check_sparsity( &
        u_pha_nonods, cv_pha_nonods, &
        u_nonods, cv_nonods, totele, &
        mx_ncolacv, ncolacv, finacv, colacv, midacv, & ! CV multi-phase eqns (e.g. vol frac, temp)
        mxnele, ncolele, midele, finele, colele, & ! Element connectivity
        mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, & ! Force balance sparsity
        mx_nct, ncolct, findct, colct, & ! CT sparsity - global cty eqn
        mx_nc, ncolc, findc, colc, & ! C sparsity operating on pressure in force balance
        mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, & ! pressure matrix for projection method
        mx_ncolm, ncolm, findm, colm, midm )

        implicit none
        integer, intent( in ) :: u_pha_nonods, cv_pha_nonods, u_nonods, cv_nonods, totele
        integer, intent ( in ) :: mx_ncolacv, ncolacv
        integer, dimension( cv_pha_nonods + 1 ), intent (in ) :: finacv
        integer, dimension( mx_ncolacv ), intent (in ) :: colacv
        integer, dimension( cv_pha_nonods ), intent (in ) :: midacv
        integer, intent ( in ) :: mxnele, ncolele
        integer, dimension( totele ), intent (in ) :: midele
        integer, dimension( totele + 1 ), intent (in ) :: finele
        integer, dimension( mxnele ), intent (in ) :: colele
        integer, intent ( in ) :: mx_ncoldgm_pha, ncoldgm_pha
        integer, dimension( mx_ncoldgm_pha ), intent (in ) :: coldgm_pha
        integer, dimension( u_pha_nonods + 1 ), intent (in ) :: findgm_pha
        integer, dimension( u_pha_nonods ), intent (in ) :: middgm_pha
        integer, intent ( in ) :: mx_nct, ncolct
        integer, dimension( cv_nonods + 1 ), intent (in ) :: findct
        integer, dimension( mx_nct ), intent (in ) :: colct
        integer, intent ( in ) :: mx_nc, ncolc
        integer, dimension( u_nonods + 1 ), intent (in ) :: findc
        integer, dimension( mx_nc ), intent (in ) :: colc
        integer, intent ( in ) :: mx_ncolcmc, ncolcmc
        integer, dimension( cv_nonods + 1 ), intent (in ) :: findcmc
        integer, dimension( mx_ncolcmc ), intent (in ) :: colcmc
        integer, dimension( cv_nonods ), intent (in ) :: midcmc
        integer, intent ( in ) :: mx_ncolm, ncolm
        integer, dimension( cv_nonods + 1 ), intent (in ) :: findm
        integer, dimension( ncolm ), intent (in ) :: colm
        integer, dimension( cv_nonods ), intent (in ) :: midm

        ! Local variables
        integer, dimension( : ), allocatable :: dummy

        ewrite(3,*) 'In check_sparsity'


        write( 15, * )'########## FINACV, COLACV, MIDACV ##################'
        write(15, * )'NCOLACV:', NCOLACV
        call checksparsity( .true., 15, NCOLACV, CV_PHA_NONODS, MX_NCOLACV, FINACV, MIDACV, COLACV  )

        write( 15, * )'########## FINELE, MIDELE, COLELE  ##################'
        write(15, * )'NCOLELE:',NCOLELE
        call checksparsity( .true., 15, NCOLELE, TOTELE, MXNELE, FINELE, MIDELE, COLELE )

        allocate( dummy( CV_NONODS ))
        write( 15, * )'########## FINDCT, COLCT ##################'
        write(15, * )'NCOLCT:', NCOLCT
        call checksparsity( .false., 15, NCOLCT, CV_NONODS, MX_NCT, FINDCT, dummy, COLCT  )
        deallocate( dummy )

        allocate( dummy( U_NONODS ))
        write( 15, * )'########## FINDC, COLC ##################'
        write(15, * )'NCOLC:', NCOLC
        call checksparsity( .false., 15, NCOLC, U_NONODS, MX_NC, FINDC, dummy, COLC )
        deallocate( dummy )

        write( 15, * )'########## FINDGM_PHA, MIDDGM_PHA, COLDGM_PHA ##################'
        write(15, * )'NCOLDGM_PHA:',NCOLDGM_PHA
        call checksparsity( .true., 15, NCOLDGM_PHA, U_PHA_NONODS, MX_NCOLDGM_PHA, FINDGM_PHA, MIDDGM_PHA, COLDGM_PHA )

        write( 15, * )'########## FINDCMC, MIDCMC, COLCMC ##################'
        write(15, * )'NCOLCMC:',NCOLCMC
        call checksparsity( .true., 15, NCOLCMC, CV_NONODS, MX_NCOLCMC, FINDCMC, MIDCMC, COLCMC )

        write( 15, * )'########## FINDM, MIDM, COLM ##################'
        write(15, * )'NCOLM:',NCOLM
        call checksparsity( .true., 15, NCOLM, CV_NONODS, MX_NCOLM, FINDM, MIDM, COLM )

        close( 15 )

        ewrite(3,*) 'Leaving check_sparsity'

        return

    end subroutine check_sparsity

    subroutine checksparsity( option_mid, unit, ncol2, nonods, ncol, find, mid, col )
        implicit none
        logical, intent( in ) :: option_mid
        integer, intent( in ) :: unit, ncol2, nonods, ncol
        integer, dimension( nonods + 1 ), intent( in ) :: find
        integer, dimension( nonods ), intent( in ) :: mid
        integer, dimension( ncol ), intent( in ) :: col
        ! Local variables
        integer :: inod, icol

        ewrite(3,*) 'In checksparsity'

        write( unit, * )'find:', ( find( inod ), inod = 1, nonods + 1 )
        if( option_mid )write( unit, * )'mid:', ( mid( inod ), inod = 1, nonods )
        write( unit, * )'col:', ( col( icol ), icol = 1, ncol2  )

        ewrite(3,*) 'Leaving checksparsity'

        return
    end subroutine checksparsity

end module spact


