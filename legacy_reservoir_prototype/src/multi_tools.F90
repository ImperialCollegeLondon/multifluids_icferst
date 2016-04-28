  
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
  

module multi_tools
    use fldebug
    use fields
    use global_parameters, only: OPTION_PATH_LEN, PYTHON_FUNC_LEN, PI, is_porous_media

    implicit none

contains

    !sprint_to_do!update to use new memory, or remove it since it is unused
    SUBROUTINE SIMPNORM( NORMX, NORMY, NORMZ, D3, &
        SNDGLN, SNLOC, NONODS, ELE, &
        X, Y, Z, &
        FINDRM, COLM)
        ! Calculate an approx normal (normx,normy,normz)
        IMPLICIT NONE

        REAL, intent( inout ) :: NORMX, NORMY, NORMZ
        LOGICAL, intent( in ) :: D3
        INTEGER, intent( in ) :: SNLOC, NONODS, ELE
        INTEGER, DIMENSION( : ), intent( in ) :: SNDGLN
        REAL, DIMENSION( :), intent( in ) :: X, Y, Z
        INTEGER, DIMENSION( : ), intent( in ) :: FINDRM
        INTEGER, DIMENSION( : ), intent( in ) :: COLM
        ! Local variables
        REAL :: XCS, YCS, ZCS, XCV, YCV, ZCV, NORM
        INTEGER :: ICV, L, IGL, L2, IGL2, COUNT, COUNT2, NODJ, NODJ2
        LOGICAL :: FOUND, SURF, ALLSUF

        XCS = 0.
        YCS = 0.
        ZCS = 0.
        DO L=1,SNLOC
            IGL = SNDGLN(( ELE - 1 ) * SNLOC + L )
            XCS = XCS + X( IGL ) / REAL( SNLOC )
            YCS = YCS + Y( IGL ) / REAL( SNLOC )
            IF(D3) ZCS = ZCS + Z( IGL ) / REAL( SNLOC )
        END DO
        ewrite(3,*)'XCS,YCS,ZCS:',XCS,YCS,ZCS

        XCV = 0.
        YCV = 0.
        ZCV = 0.
        ICV = 0
        L = 1
        IGL = SNDGLN(( ELE - 1 ) * SNLOC + L )

        Loop_Row: DO COUNT = FINDRM( IGL ) , FINDRM( IGL + 1 ) - 1, 1
            NODJ = COLM( COUNT )
            Inside_matrix: IF(NODJ <= NONODS) THEN
                ! Make sure its not a surface node of surface element ELE
                SURF = .FALSE.
                DO L2 = 1, SNLOC
                    IGL2 = SNDGLN(( ELE - 1 ) * SNLOC + L2 )
                    IF( IGL2 == NODJ ) SURF = .TRUE.
                END DO

                Cond_Surf: IF( .NOT. SURF) THEN
                    ! make sure NODJ is connected to all surface nodes and is not a surface node for
                    ! surface element.
                    ALLSUF = .TRUE.
                    DO L2 = 1, SNLOC
                        IGL2 = SNDGLN(( ELE - 1 ) * SNLOC + L2 )
                        FOUND = .FALSE.
                        DO COUNT2 = FINDRM( NODJ ) , FINDRM( NODJ + 1 ) - 1, 1
                            NODJ2 = COLM( COUNT2 )
                            IF( IGL2 == NODJ2) FOUND = .TRUE.
                        END DO
                        IF( .NOT. FOUND) ALLSUF = .FALSE.
                    END DO

                    IF( ALLSUF ) THEN
                        XCV = XCV + X( NODJ )
                        YCV = YCV + Y( NODJ )
                        IF( D3 ) ZCV = ZCV + Z( NODJ )
                        ICV = ICV + 1
                    ENDIF

                ENDIF Cond_Surf
            ENDIF Inside_matrix

        END DO Loop_Row

        XCV = XCV / REAL( ICV )
        YCV = YCV / REAL( ICV )
        IF( D3 ) ZCV = ZCV / REAL( ICV )
        NORMX = XCS - XCV
        NORMY = YCS - YCV
        NORMZ = ZCS - ZCV

        NORM = NORMX **2 + NORMY **2
        IF( D3 ) NORM = NORM + NORMZ ** 2
        NORM = MAX( 1.E-15, SQRT( NORM ))

        NORMX = NORMX / NORM
        NORMY = NORMY / NORM
        IF( D3 ) NORMZ = NORMZ / NORM

        RETURN

    END SUBROUTINE SIMPNORM


    REAL FUNCTION R2NORM( VEC, NVEC )
        IMPLICIT NONE
        INTEGER :: NVEC
        REAL, DIMENSION( NVEC ) :: VEC
        ! Local variables

        R2NORM = SQRT( SUM( VEC**2 ) )

        RETURN
    END FUNCTION R2NORM

    real function ptolfun(value)
        ! This function is a tolerance function for strictly positive values used as a denominator.
        ! If the value of VALUE less than 1E-10, then it returns TOLERANCE otherwise VALUE.

        implicit none
        real, intent(in) :: value
        ! Local
        real, parameter :: tolerance = 1.e-10

        ptolfun = max( tolerance, value )

        return

    end function ptolfun

    PURE function tolfun(value) result(v_tolfun)
    !real function tolfun(value)
        ! This function is a tolerance function for a value which is used as a denominator.
        ! If the absolute value of VALUE less than 1E-10, then it returns SIGN(A,B) i.e.
        ! the absolute value of A times the sign of B where A is TOLERANCE and B is VALUE.

        implicit none
        real, intent(in) :: value
        ! Local
        real :: v_tolfun
        real, parameter :: tolerance = 1.e-10

        v_tolfun = sign( 1.0, value ) * max( tolerance, abs(value) )

        return

    end function tolfun



    PURE function tolfun_many(val) result(v_tolfun)

        implicit none
        real, dimension(:), intent(in) :: val
        real, dimension(size(val)) :: v_tolfun
        ! Local
        real, parameter :: tolerance = 1.e-10

        v_tolfun = sign( 1.0, val ) * max( tolerance, abs(val) )

        return

    end function tolfun_many

    function tetvolume(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3)
        IMPLICIT NONE

        real, intent(in) :: x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3

        real :: tetvolume

        ! tetvolume = 1.0 / 6.0 * det |three tet edge vectors|
        ! Chris' tets have a clockwise base, hence the sign change in the det
        tetvolume = &
            (  &
            - (x1 - x0) * ((y2 - y0) * (z3 - z0) - (y3 - y0) * (z2 - z0)) &
            + (y1 - y0) * ((x2 - x0) * (z3 - z0) - (x3 - x0) * (z2 - z0)) &
            - (z1 - z0) * ((x2 - x0) * (y3 - y0) - (x3 - x0) * (y2 - y0)) &
            ) / 6.0

    end function tetvolume

    PURE function vtolfun(val) result(v_tolfun)

        implicit none
        real, dimension(:), intent(in) :: val
        real, dimension(size(val)) :: v_tolfun
        ! Local
        real, parameter :: tolerance = 1.e-10

        v_tolfun = sign( 1.0, val ) * max( tolerance, abs(val) )

        return

    end function vtolfun

    subroutine nan_check(a,k)
    !Checks if a number is a Nan
        real :: a
        integer :: k

        if (a/=a) then
            print*, 'nan found! loop:', k
        end if

    end subroutine nan_check

    subroutine nan_check_arr(a,k)
    !Checks if an array is a Nan
        real, dimension(:,:) :: a
        integer :: k

        if (any(a/=a)) then
            print*, 'nan found! loop:', k
        end if

    end subroutine nan_check_arr


    PURE FUNCTION NVDFUNNEW_MANY( UF, UC, XI_LIMIT ) result(nvd_limit)
        implicit none
        ! The function computes NVDFUNNEW, the normalised value of the
        ! advected variable on the face of the control volume, based on
        ! the normalised value of the advected variable in the donor CV,
        ! UC, and the high-order estimate of the face value UF.
        ! NVDFUNNEW is limited so that it is in the non-oscillatory
        ! region of normalised variable diagram (NVD).
        !
        ! XI is the parameter in equation 38 of the Riemann paper. If XI is equal
        ! to 2 then this corresponds to a TVD condition in 1-D, a value of XI
        ! equal to 3 has been recommended elsewhere
        !
        REAL, DIMENSION( : ), intent(in)  :: UC, UF, XI_LIMIT
        real, dimension(size(uc)) :: nvd_limit
        logical, PARAMETER :: orig_limit=.false. ! original limiter is less invasive.

        ! For the region 0 < UC < 1 on the NVD, define the limiter
        if(orig_limit) then
            where( ( UC > 0.0 ) .AND. ( UC < 1.0 ) )
                nvd_limit = MIN( 1.0, XI_LIMIT * UC, MAX( 0.0, UF ) )
            !       nvd_limit = MIN( 1.0, XI_LIMIT * UC, MAX( UC, UF ) )
            !      nvd_limit= MAX(  MIN(UF, XI_LIMIT*UC, 1.0), UC)
            ELSE where ! Outside the region 0<UC<1 on the NVD, use first-order upwinding
                nvd_limit = UC
            END where
        else
            nvd_limit= MAX(  MIN(UF, XI_LIMIT*UC, 1.0), UC)
        endif

    end function nvdfunnew_many




    FUNCTION NVDFUNNEW_MANY_sqrt( UF, UC, XI_LIMIT ) result(nvd_limit)
        implicit none
        ! The function computes NVDFUNNEW, the normalised value of the
        ! advected variable on the face of the control volume, based on
        ! the normalised value of the advected variable in the donor CV,
        ! UC, and the high-order estimate of the face value UF.
        ! NVDFUNNEW is limited so that it is in the non-oscillatory
        ! region of normalised variable diagram (NVD).
        !
        ! XI is the parameter in equation 38 of the Riemann paper. If XI is equal
        ! to 2 then this corresponds to a TVD condition in 1-D, a value of XI
        ! equal to 3 has been recommended elsewhere
        !
        REAL, DIMENSION( : ), intent(in)  :: UC, UF, XI_LIMIT
        real, dimension(size(uc)) :: nvd_limit


        nvd_limit= MAX(  MIN(UF, XI_LIMIT*UC, 1.0), UC)

    end function nvdfunnew_many_sqrt



    function get_ndglno(mesh) result(ndglno)
        type(mesh_type) :: mesh
        integer, dimension(:), pointer  ::  ndglno

        ndglno=> mesh%ndglno
    end function get_ndglno

    !Sort a list in increasing order
    !Vec is the vector to sort and n is an starting point, like 1
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

    contains
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

    end subroutine quicksort


    SUBROUTINE CALC_FACE_ELE( FACE_ELE, TOTELE, STOTEL, NFACE, &
        FINELE, COLELE, CV_NLOC, CV_SNLOC, CV_NONODS, CV_NDGLN, CV_SNDGLN, &
        CV_SLOCLIST, X_NLOC, X_NDGLN)
        ! Calculate FACE_ELE - the list of elements surrounding an
        ! element and referenced with a face -ve values correspond to surface elements.
        IMPLICIT NONE
        INTEGER, intent( in ) :: TOTELE, STOTEL, NFACE, CV_NLOC, CV_SNLOC, CV_NONODS, &
            X_NLOC
        INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: X_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
        INTEGER, DIMENSION( :, : ), intent( in ) :: CV_SLOCLIST
        INTEGER, DIMENSION( :, : ), intent( inout ) :: FACE_ELE
        INTEGER, DIMENSION( : ), intent( in ) :: FINELE
        INTEGER, DIMENSION( : ), intent( in ) :: COLELE
        ! Local variables
        LOGICAL, DIMENSION( : ), allocatable :: NOD_BELONG_ELE
        INTEGER, DIMENSION( : ), allocatable :: NOD_COUNT_SELE, FIN_ND_SELE, COL_ND_SELE, ELE_ROW
        LOGICAL, DIMENSION( : ), allocatable :: NOD_ON_BOUNDARY
        INTEGER :: ELE, ELE2, SELE, CV_ILOC, CV_INOD, COUNT, IFACE, CV_SILOC, JFACE, &
            CV_ILOC2, CV_INOD2, CV_SILOC2, IFACE2, SELE2
        LOGICAL :: FOUND, SELE_FOUND, GOT_ALL

        ALLOCATE( NOD_BELONG_ELE( CV_NONODS ))
        ALLOCATE( NOD_COUNT_SELE( CV_NONODS ))
        ALLOCATE( FIN_ND_SELE( CV_NONODS + 1 ))
        ALLOCATE( NOD_ON_BOUNDARY( CV_NONODS ))
        ALLOCATE( ELE_ROW( NFACE ))

        FACE_ELE=0
        DO ELE=1,TOTELE
            IFACE=0
            DO COUNT=FINELE(ELE),FINELE(ELE+1)-1
                ELE2=COLELE(COUNT)
                IF(ELE2.NE.ELE) THEN
                    IFACE=IFACE+1
                    FACE_ELE(IFACE,ELE)=ELE2
                END IF
            END DO
        END DO

        NOD_COUNT_SELE=0
        NOD_ON_BOUNDARY=.FALSE.
        DO SELE=1,STOTEL
            DO CV_SILOC=1,CV_SNLOC
                CV_INOD=CV_SNDGLN((SELE-1)*CV_SNLOC+CV_SILOC)
                !ewrite(3,*)'sele,CV_INOD:',sele,CV_INOD
                NOD_COUNT_SELE(CV_INOD)=NOD_COUNT_SELE(CV_INOD)+1
                NOD_ON_BOUNDARY(CV_INOD)=.TRUE.
            END DO
        END DO

        FIN_ND_SELE(1)=1
        DO CV_INOD=2,CV_NONODS+1
            FIN_ND_SELE(CV_INOD)=FIN_ND_SELE(CV_INOD-1)+NOD_COUNT_SELE(CV_INOD-1)
        END DO
        ALLOCATE(COL_ND_SELE(FIN_ND_SELE(CV_NONODS+1)-1))

        NOD_COUNT_SELE=0
        DO SELE=1,STOTEL
            DO CV_SILOC=1,CV_SNLOC
                CV_INOD=CV_SNDGLN((SELE-1)*CV_SNLOC+CV_SILOC)
                NOD_COUNT_SELE(CV_INOD)=NOD_COUNT_SELE(CV_INOD)+1
                COL_ND_SELE(FIN_ND_SELE(CV_INOD)-1+NOD_COUNT_SELE(CV_INOD))=SELE
            END DO
        END DO

        ! Put surface elements into FACE_ELE
        NOD_BELONG_ELE=.FALSE.
        DO ELE=1,TOTELE
            DO CV_ILOC=1,CV_NLOC
                CV_INOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                NOD_BELONG_ELE(CV_INOD)=.TRUE.
            END DO
            DO IFACE=1,NFACE
                IF(FACE_ELE(IFACE,ELE)==0) THEN
                    DO CV_ILOC=1,CV_NLOC
                        CV_INOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                        DO COUNT=FIN_ND_SELE(CV_INOD),FIN_ND_SELE(CV_INOD+1)-1
                            SELE=COL_ND_SELE(COUNT)
                            FOUND=.FALSE.
                            DO JFACE=1,IFACE-1
                                IF(SELE == - FACE_ELE(JFACE,ELE)) THEN
                                    FOUND=.TRUE.
                                    EXIT
                                ENDIF
                            END DO
                            IF(.NOT.FOUND) THEN ! SELE is a candidate.
                                SELE_FOUND=.TRUE.
                                DO CV_SILOC=1,CV_SNLOC
                                    CV_INOD=CV_SNDGLN((SELE-1)*CV_SNLOC+CV_SILOC)
                                    IF(.NOT.NOD_BELONG_ELE(CV_INOD)) THEN
                                        SELE_FOUND=.FALSE.
                                        EXIT
                                    ENDIF
                                END DO
                                IF(SELE_FOUND) FACE_ELE(IFACE,ELE)= - SELE
                            ENDIF
                        END DO
                    END DO
                ENDIF
            END DO
            DO CV_ILOC=1,CV_NLOC
                CV_INOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                NOD_BELONG_ELE(CV_INOD)=.FALSE.
            END DO
        END DO

        !do ele=1,totele
        !   ewrite(3,*)'ele',ele,' FACE_ELE(IFACE,ELE):',(FACE_ELE(IFACE,ELE),iface=1,nface)
        !end do

        ! ***********************************************************
        ! Now re-arrange ordering of elements in FACE_ELE so they are
        ! consistent with the local nodes in CV_SLOCLIST
        DO ELE=1,TOTELE
            DO IFACE=1,NFACE
                ! Find the suitable element in row ELE of FACE_ELE for face IFACE.
                DO IFACE2=1,NFACE
                    ELE2=FACE_ELE(IFACE2,ELE)
                    SELE2 = MAX( 0, - ELE2 )
                    ELE2 = MAX( 0, + ELE2 )
                    GOT_ALL=.TRUE.
                    DO CV_SILOC=1,CV_SNLOC
                        CV_ILOC=CV_SLOCLIST(IFACE,CV_SILOC)
                        FOUND=.FALSE.
                        IF(SELE2 == 0 .and. ele2>0) THEN ! is a volume element
                            CV_INOD=X_NDGLN((ELE-1)*X_NLOC+CV_ILOC)
                            DO CV_ILOC2=1,CV_NLOC
                                CV_INOD2=X_NDGLN((ELE2-1)*X_NLOC+CV_ILOC2)
                                IF(CV_INOD == CV_INOD2) THEN
                                    FOUND=.TRUE.
                                    EXIT
                                ENDIF
                            END DO
                        ELSE if (sele2>0) then ! is a surface element
                            CV_INOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                            DO CV_SILOC2=1,CV_SNLOC
                                CV_INOD2=CV_SNDGLN((SELE2-1)*CV_SNLOC+CV_SILOC2)
                                IF(CV_INOD == CV_INOD2) THEN
                                    FOUND=.TRUE.
                                    EXIT
                                ENDIF
                            END DO
                        ENDIF
                        IF(.NOT.FOUND) THEN
                            GOT_ALL=.FALSE.
                            EXIT
                        ENDIF
                    END DO
                    IF(GOT_ALL) ELE_ROW(IFACE)=FACE_ELE(IFACE2,ELE)
                END DO
            END DO
            ! Re-order row...
            FACE_ELE(:,ELE)=ELE_ROW(:)
        END DO


        DEALLOCATE( NOD_BELONG_ELE )
        DEALLOCATE( NOD_COUNT_SELE )
        DEALLOCATE( FIN_ND_SELE )
        DEALLOCATE( NOD_ON_BOUNDARY )
        DEALLOCATE( COL_ND_SELE )
        DEALLOCATE( ELE_ROW )

        RETURN

    END SUBROUTINE CALC_FACE_ELE

    !sprint_to_do; we need to see how this behaves with many regions ids
    subroutine assign_val(outval,inval)
        !Copies the data from inval to outval safely.
        !If the sizes are different outval is populated using the first value of inval
        implicit none
        real, dimension(:), intent(inout) :: outval
        real, dimension(:), intent(in) :: inval

        if (size(outval)/=size(inval)) then
            outval = inval(1)
        else
            outval = inval
        end if
    end subroutine assign_val


    real function table_interpolation(X_points, Y_points, input_X)
        !This function returns the value at position input_X using
        !X_points, Y_points to form a linear (size == 2) or quadratic (size == 3) interpolation
        implicit none
        real, intent(in) :: input_X
        real, dimension(:), intent(in) :: X_points, Y_points
        !Local variables
        integer :: k, j, siz
        real, dimension(size(X_points),size(X_points)) :: A
        real, dimension(size(X_points)) :: b

        siz = size(X_points)

        do k = 1, siz
            do j = 1, siz!Fill columns
                A(j,k) = X_points(j)**real(siz-k)
                if (j==k) A(j,k) = A(j,k) + A(j,k)*(1d-7*real(rand(j)))!to avoid non-invertible matrices
            end do
            b(k) = Y_points(k)
        end do
        !Solve system
        call invert(A); b = matmul(A, b)
        if (siz == 2) then!Linear approximation
            table_interpolation = b(1) * input_X + b(2)
        else!quadratic approximation
            table_interpolation = b(1) * input_X**2. + b(2) * input_X + b(3)
        end if
    end function table_interpolation

    subroutine read_csv_table(data_array, path_to_table, extra_data)
        !Template of csv table
        !OPTIONAL section (header)
        !real1,real2,real3,..., size(extra_data)
        !2,3                 !<= this are the columns and rows
        !Pressure,Saturation
        !1000,0.9
        !250,0.5
        !100,0.1
        implicit none
        real, dimension(:,:), allocatable, intent(inout) :: data_array
        character( len = option_path_len ), intent(in) :: path_to_table
        real, optional, dimension(:), intent(inout) :: extra_data
        !Local variables
        integer :: i, ierr
        integer, dimension(2) :: table_size
        !Open file
        open(unit= 89, file=trim(path_to_table)//".csv", status='old', action='read')
        if (present(extra_data)) then
            !The extra data is composed of one line of headers that we ignore plus one extra line with the data
            read(89,*)!skip header
            read(89,*) extra_data
        end if

        !CSV table must start with the number of columns by rows
        read(89,*) table_size
        allocate(data_array(table_size(1), table_size(2)))
        !Skip the headers
        read(89,*);
        !Read the table
        do i = 1, table_size(2)
            read(89,*, IOSTAT=ierr) data_array(1:table_size(1), i)
            if (ierr<0) exit
        end do
        close(89)
    end subroutine read_csv_table

    !Subroutine to print Arrays by (columns,rows)
    !Matrix = 2D Array
    subroutine printMatrix(Matrix)
        implicit none

        Integer :: length,i,j, k
        character (len=1000000) :: cadena
        character (len=100) :: aux
        real, intent(in), dimension(:,:):: Matrix
        !Local
        real, dimension(size(matrix,2),size(matrix,1)) :: auxMatrix

        auxMatrix = transpose(Matrix)

        length = size(auxMatrix,2);
        do i = 1,size(auxMatrix,1)
            print *,""
            cadena = ""
            do j = 1 , length
                write(aux,*), auxMatrix(i,j)
                k = index(trim(aux),"E",.true.)
                if (k/=0) then
                    aux = aux(1:k-6)//trim(aux(k:))
                end if

                cadena = trim(cadena)//' '//trim(aux)
            end do
            print '(A $)', trim(cadena)
        end do

        print *,"";
    end subroutine PrintMatrix

end module multi_tools


