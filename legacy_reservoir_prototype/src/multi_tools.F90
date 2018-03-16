
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
    use vector_tools

    implicit none

    type bad_elements
        integer :: bad_ele
        real :: angle
        real :: perp_height ! perp height from base (assuming an isosceles triangle)
        real, allocatable, dimension(:,:) :: rotmatrix ! the rotation matrix to 'stretch' the bad element in the direction normal to the big angle
        real :: base ! length of side opposite the large angle
    end type


contains

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
        !rows,columns
        !2,3
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
        read(89,*)!skip header
        !CSV table must start with the number of columns by rows
        read(89,*) table_size
        allocate(data_array(table_size(1), table_size(2)))
        !Skip the headers
        read(89,*);
        !Read the table
        do i = 1, table_size(2)
            read(89,*, IOSTAT=ierr) data_array(1:table_size(1), i)
            if (ierr/=0) exit
        end do
        close(89)
    end subroutine read_csv_table


    subroutine extract_strings_from_csv_file(csv_table_strings, path_to_table, Nentries)
        !This subroutine reads a csv file and returns
        implicit none
        integer, intent(out) :: Nentries
        character( len = option_path_len ), intent(in) :: path_to_table
        character(len=option_path_len), dimension(:,:),  allocatable, intent(out) :: csv_table_strings
        !Local variables
        integer :: i, ierr, start
        !Allocate table
        allocate(csv_table_strings(10000,4))!This should be enough

        !Open file
        open(unit= 89, file=trim(path_to_table)//'.csv', status='old', action='read')
        !CSV table must start with the number of columns by rows
        !Read the table
        i = 1
        start = 1
        do while (.true.)
            read(89,*, IOSTAT=ierr) csv_table_strings(i,:)!cadena
            if (ierr/=0) exit
            !Extract four strings from cadena
            i = i + 1
        end do
        Nentries = i-1
        close(89)

    end subroutine extract_strings_from_csv_file

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


    subroutine CheckElementAngles(X_ALL, totele, x_ndgln, X_nloc, MaxAngle, Quality_list, bad_element_flag, number_bad)
        !This function checks the angles of an input element. If one angle is above
        !the Maxangle it will be true in the list

        Implicit none
        !Global variables
        type(bad_elements), dimension(:), intent(inout) :: Quality_list
        integer, optional, intent(out) :: number_bad
        real, intent (in) :: MaxAngle
        integer, dimension(:), intent(in) :: x_ndgln
        integer, intent(in) :: x_nloc, totele
        logical, intent(in) :: bad_element_flag
        !Local variables
        integer :: ELE, i, k
        logical :: Bad_found
        real :: MxAngl
        !Definition of Pi
        real, dimension(:,:), pointer:: X_ALL
        real, parameter :: pi = acos(0.d0) * 2d0

        !Prepare data
        do i = 1, size(Quality_list)
            Quality_list(i)%bad_ele = -1!Initialize with negative values
        end do

        !Convert input angles to radians
        MxAngl = pi/180. * MaxAngle
        k = 1 ! number of bad elements
        if (size(X_ALL,1)==2) then!2D triangles
            do ele = 1, totele
                Bad_found = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 1, 2, 3, MxAngl, Quality_list(ele))
                if (Bad_found) then
                    Quality_list(k)%bad_ele = ele
                    k = k + 1
                end if
            end do
        else if(size(X_ALL,1)==3) then!3D tetrahedra
            !adjust to match the 2D case once that one works properly
            do ele = 1, totele
                !We check the 4 triangles that form a tet
                Bad_found = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 1, 2, 3, MxAngl, Quality_list(ele), 4)
                if (Bad_found) then
                    Quality_list(k)%bad_ele = ele
                    k = k + 1
                else
                    Bad_found = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 1, 2, 4, MxAngl, Quality_list(ele), 3)
                    if (Bad_found) then
                        Quality_list(k)%bad_ele = ele
                        k = k + 1
                    else
                        Bad_found = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 1, 4, 3, MxAngl, Quality_list(ele), 2)
                        if (Bad_found) then
                            Quality_list(k)%bad_ele = ele
                            k = k + 1
                        else
                            Bad_found = Check_element(X_ALL, x_ndgln, (ele-1)*X_nloc, 4, 2, 3, MxAngl, Quality_list(ele), 1)
                            if (Bad_found) then
                                Quality_list(k)%bad_ele = ele
                                k = k + 1
                            end if
                        end if
                    end if
                end if
            end do
        end if

        if (present(number_bad)) then
            number_bad = k-1 ! pass for the diagnostic table
        else
            if (float(k-1)/totele > 0.01) then ! if bad elements found, let the user know
                ewrite(0,'(A,F6.2,A,G7.2,A,I3,A)') 'Bad Element Fix WARNING:', (float(k-1)/totele)*100., '% of the elements, (',  k  ,') are bad elements, with an angle larger than ', int(MaxAngle), ' degrees. Artificial permeability anisotropy introduced, results may not be reliable.'
            end if
        end if
            contains

logical function Check_element(X_ALL, x_ndgln, ele_Pos, Pos1, Pos2, Pos3, MaxAngle, Quality_list, Pos4)
    !Checks if an angle is below a threshold and stores the information to
    implicit none
    real, dimension(:,:), intent(in) :: X_ALL
    integer, intent(in) :: ele_Pos, Pos1, Pos2, Pos3
    integer, optional :: Pos4
    real, intent(in) :: MaxAngle
    integer, dimension(:), intent(in) :: x_ndgln
    type(bad_elements), intent(inout) :: Quality_list
    !Local variables
    real, dimension(size(X_ALL,1)) :: X1, X2, X3, X4
    real, dimension(3) :: alpha, length
    real, dimension(3) :: normal ! normal vector of edge opposite largest angle
    !Definition of Pi
    real, parameter :: pi = acos(0.0) * 2.0

    ! initialise normal
    normal(:) = 0

    Check_element = .false.
    !Define the vertexes
    X1 = X_ALL(:, x_ndgln(ele_Pos+Pos1))
    X2 = X_ALL(:, x_ndgln(ele_Pos+Pos2))
    X3 = X_ALL(:, x_ndgln(ele_Pos+Pos3))

    !Define the 4th vertex if 3D to calculate normal
    if(size(X_ALL,1)==3) then
        X4 = X_ALL(:, x_ndgln(ele_Pos+Pos4))
    end if

    !Calculate the length of the edges
    length(1) = sqrt(dot_product(X1(:)-X3(:), X1(:)-X3(:)))
    length(2) = sqrt(dot_product(X1(:)-X2(:), X1(:)-X2(:)))
    length(3) = sqrt(dot_product(X2(:)-X3(:), X2(:)-X3(:)))

    !Determine the angle between each edge
    alpha(1) = acos((length(3)**2+length(2)**2-length(1)**2)/(2. *length(3)*length(2)))
    alpha(2) = acos((length(1)**2+length(3)**2-length(2)**2)/(2. *length(1)*length(3)))
    alpha(3) = pi - alpha(1)-alpha(2)

    !Check angles
    if (alpha(1)>=MaxAngle) then
        !Store angle so later the over-relaxation can depend on this
        Quality_list%angle = alpha(1) * 180 / pi
!        perp_height = (length(1)/2) / (tan(alpha(1)/2))
!        Quality_list%perp_height = perp_height
!        Quality_list%base = length(1)

        if (bad_element_flag) then
            ! calculate normal. For 3D its the cross product, 2D will be midpoint of edge opp largest angle and angle node
            if(size(X_ALL,1)==3) then
                normal = cross_product(X3-X1,X4-X3)
            else
                normal(1:2) = X2-(X3+X1)/2
            end if
            if (.not. allocated(Quality_list%rotmatrix)) allocate(Quality_list%rotmatrix(3,3))
            call RotationMatrix(normal,Quality_list%rotmatrix)
        end if
        Check_element = .true.
    else if (alpha(2)>= MaxAngle) then
        !Store angle so later the over-relaxation can depend on this
        Quality_list%angle = alpha(2) * 180 / pi
!        perp_height = (length(2)/2) / (tan(alpha(2)/2))
!        Quality_list%perp_height = perp_height
!        Quality_list%base = length(2)

        if (bad_element_flag) then
            ! calculate normal. For 3D its the cross product, 2D will be midpoint of edge opp largest angle and angle node
            if(size(X_ALL,1)==3) then
                normal = cross_product(X2-X1,X4-X2)
            else
                normal(1:2) = X3-(X2+X1)/2
            end if
            if (.not. allocated(Quality_list%rotmatrix)) allocate(Quality_list%rotmatrix(3,3))
            call RotationMatrix(normal,Quality_list%rotmatrix)
        end if

        Check_element = .true.
    else if (alpha(3) >= MaxAngle) then
        !Store angle so later the over-relaxation can depend on this
        Quality_list%angle = alpha(3) * 180 / pi
!        perp_height = (length(3)/2) / (tan(alpha(3)/2))
!        Quality_list%perp_height = perp_height
!        Quality_list%base = length(3)

        if (bad_element_flag) then
            ! calculate normal. For 3D its the cross product, 2D will be midpoint of edge opp largest angle and angle node
            if(size(X_ALL,1)==3) then
                normal = cross_product(X2-X3,X4-X2)
            else
                normal(1:2) = X3-(X2+X3)/2
            end if
            if (.not. allocated(Quality_list%rotmatrix)) allocate(Quality_list%rotmatrix(3,3))
            call RotationMatrix(normal,Quality_list%rotmatrix)
        end if
        Check_element = .true.
    end if

end function Check_element

    end subroutine CheckElementAngles

subroutine RotationMatrix(a,R)

    real, dimension(3), intent(in)    :: a ! normal vector of length opposite the largest angle of a 'bad' element
    real, dimension(3)                :: an, bn, cx, u, v, w !an/bn normal vectors, cx cross product, u/v/w basis vectors
    real                              :: dot
    real, dimension(3,3)              :: T,G !T is the rotation matrix around z-axis !Change of basis matrix
    real, dimension(3,3), intent(out) :: R !R = G*T*inv_G
    integer :: d ! dimension of model 2 or 3

    ! Normalise vector
    an = a/NORM2(a)

    ! normal vector of our cartesian aligned element
    if (d ==3) then! 3D
        bn = [ 0, 0, 1 ]
    else ! 2D
        bn = [ 0, 1, 0]
    end if

    ! calculate dot product and cross product
    dot = dot_product(an,bn)
    cx = cross_product(an,bn)

    ! find the rotation matrix on cartesian axis around z axis
    T(1,1) = dot
    T(1,2) = -NORM2(cx)
    T(1,3) = 0
    T(2,1) = NORM2(cx)
    T(2,2) = dot
    T(2,3) = 0
    T(3,1) = 0
    T(3,2) = 0
    T(3,3) = 1


    if (all(cx(:)==0)) then ! if normal vectors are parallel then no need for new basis

        R = T

    else

        ! new basis
        u = an ! Normalised vector projection of bn onto an
        v = (bn - dot*an)/NORM2(bn - dot*an) ! Normalised vector rejection of bn onto an
        w = cx/NORM2(cx) ! Normalized cross product of an and bn creating the new orthonormal basis

        ! Basis change matrix
        G(:,1) = u
        G(:,2) = v
        G(:,3) = w

        ! Calculate rotation matrix on the element basis
        R = matmul(G,matmul(T, inverse(G)))

        !write(*,'(A)') 'R='
        !do i=1,size(R,1)
        !  write(*,'(20G12.4)') R(i,:)
        !end do
    end if

END subroutine RotationMatrix


    subroutine read_nastran_file(filepath, node, edges)
        !This subroutine reads a nastran file that contains the information defining the 1D path of a well
        !the input relative filepath should include the file format, for example: well.bdf
        implicit none
        character( len = * ), intent(in) :: filepath
        real, dimension(:,:), allocatable, intent(inout) :: node
        integer, dimension(:,:), allocatable, intent(inout) :: edges
        !Local variables
        integer, dimension(:), allocatable:: conversor
        integer :: i, k, j
        character( len = option_path_len ):: cadena
        integer :: Nnodes, Nedges
        real, dimension(4) :: edge_line

        !First we need to get the number of nodes and the number of edges to correctly allocate node and edges
        call get_nodes_edges(Nnodes, Nedges)!; Nedges = Nnodes - 1!In 1d there is always one edge less than nodes
        allocate(node(3, Nnodes), edges(2, Nedges), conversor(Nnodes))
        !Open file
        open(unit= 89, file=trim(filepath)//".bdf", status='old', action='read')
        cadena = "--"
        !skip header until reaching the nodes
        do while (cadena(1:5)/="GRID")
            read(89,'(A)') cadena
        end do
        i = 1!Read nodes
        do while (cadena(1:5)=="GRID")
            !Nastran files from trellis seem to have a bug, not leaving a space between coordinate 1 and 2
            !we extract numbers by position as these are defined
            read(cadena(25:32),*) node(1,i)
            read(cadena(33:40),*) node(2,i)
            read(cadena(41:49),*) node(3,i)
            !Extract as well the number of the node
            read(cadena(9:12),*) conversor(i)
            read(89,'(A)') cadena; i = i + 1!read line and advance the counter
        end do
        !skip until reaching the edges
        do while (cadena(1:5)/="CROD")
            read(89,'(A)') cadena
        end do
        i = 1!Read edges
        do while (cadena(1:5)=="CROD")
            read(cadena(6:len(cadena)),*) edge_line
            edges(:,i) = int(edge_line(3:4))!Only the last two columns contains the connection between nodes
            read(89,'(A)') cadena; i = i + 1!read line and advance the counter
        end do
        close(89)

        !Before leaving we normalize the edges list, making it to go from 1 to last edge instead of the numeration used
         do j = 1, size(edges,1)
             do i = 1, size(edges,2)
                conversor_loop: do k = 1, size(conversor)
                    if (edges(j,i) == conversor(k)) then
                        edges(j,i) = k
                        exit conversor_loop
                    end if
                end do conversor_loop
             end do
         end do
        deallocate(conversor)
    !Print well to gnuplot format
!print*, "WELL IN GNUPLOT"
!do j = 1, size(edges,1)
!do i = 1, size(edges,2)
!print *, node(:,edges(j,i))
!end do
!end do
!print *, "-----------------------"
!read*
    contains
        subroutine get_nodes_edges(Nnodes, Nedges)
            !Get the number of nodes and edges that conform the well
            implicit none
            integer, intent(inout)::Nnodes, Nedges

            Nnodes = 0; Nedges = 0
            !Open file
            open(unit= 89, file=trim(filepath)//".bdf", status='old', action='read')
            cadena = "--"
            do while (cadena(1:5)/="CROD")
                read(89,'(A)') cadena
                if (cadena(1:5)=="GRID")Nnodes = Nnodes + 1
            end do
            cadena = "--"
            Nedges = 1!Because the previous loop finished with CROD
            do while (cadena(1:5)/="PROD")
                read(89,'(A)') cadena
                if (cadena(1:5)=="CROD")Nedges = Nedges + 1
            end do
            close(89)

        end subroutine get_nodes_edges
    end subroutine read_nastran_file


end module multi_tools


