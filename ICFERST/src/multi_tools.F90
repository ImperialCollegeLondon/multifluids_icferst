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


module multi_tools
#include "petsc/finclude/petsc.h"
  use petsc

    use python_state
    use state_module
    use fldebug
    use futils
    use spud
    use fields
    use global_parameters, only: OPTION_PATH_LEN, PYTHON_FUNC_LEN, PI, is_porous_media
    use vector_tools

    implicit none

#include "petsc_legacy.h"


    !>@brief: Type to keep an eye on the quality of the elements
    !>@DEPRECATED
    type bad_elements
        integer :: bad_ele
        real :: angle
        real :: perp_height !> perp height from base (assuming an isosceles triangle)
        real, allocatable, dimension(:,:) :: rotmatrix !> the rotation matrix to 'stretch' the bad element in the direction normal to the big angle
        real :: base !> length of side opposite the large angle
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

    !>@brief:This function is a tolerance function for strictly positive values used as a denominator.
    !> If the value of VALUE less than 1E-10, then it returns TOLERANCE otherwise VALUE.
    real function ptolfun(value)

        implicit none
        real, intent(in) :: value
        ! Local
        real, parameter :: tolerance = 1.e-10

        ptolfun = max( tolerance, value )

        return

    end function ptolfun

    !>@brief:This function is a tolerance function for a value which is used as a denominator.
    !> If the absolute value of VALUE less than 1E-10, then it returns SIGN(A,B) i.e.
    !> the absolute value of A times the sign of B where A is TOLERANCE and B is VALUE.
    PURE function tolfun(value) result(v_tolfun)
    !real function tolfun(value)

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

    !!>@brief: VOlume of a tetrahedron defined by coordinates
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

    !>@brief:Checks if a number is a Nan
    subroutine nan_check(a,k)
        real :: a
        integer :: k

        if (a/=a) then
            print*, 'nan found! loop:', k
        end if

    end subroutine nan_check

    !>@brief:Checks if an array is a Nan
    subroutine nan_check_arr(a,k)
        real, dimension(:,:) :: a
        integer :: k

        if (any(a/=a)) then
            print*, 'nan found! loop:', k
        end if

    end subroutine nan_check_arr


    !>@brief: The function computes NVDFUNNEW, the normalised value of the
    !> advected variable on the face of the control volume, based on
    !> the normalised value of the advected variable in the donor CV,
    !> UC, and the high-order estimate of the face value UF.
    !> NVDFUNNEW is limited so that it is in the non-oscillatory
    !> region of normalised variable diagram (NVD).
    !>
    !> XI is the parameter in equation 38 of the Riemann paper. If XI is equal
    !> to 2 then this corresponds to a TVD condition in 1-D, a value of XI
    !> equal to 3 has been recommended elsewhere
    PURE FUNCTION NVDFUNNEW_MANY( UF, UC, XI_LIMIT ) result(nvd_limit)
        implicit none
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




    !>@brief: The function computes NVDFUNNEW, the normalised value of the
    !> advected variable on the face of the control volume, based on
    !> the normalised value of the advected variable in the donor CV,
    !> UC, and the high-order estimate of the face value UF.
    !> NVDFUNNEW is limited so that it is in the non-oscillatory
    !> region of normalised variable diagram (NVD).
    !>
    !> XI is the parameter in equation 38 of the Riemann paper. If XI is equal
    !> to 2 then this corresponds to a TVD condition in 1-D, a value of XI
    !> equal to 3 has been recommended elsewhere
    FUNCTION NVDFUNNEW_MANY_sqrt( UF, UC, XI_LIMIT ) result(nvd_limit)
        implicit none
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

    !>@brief:Sort a list in increasing order
    !>Vec is the vector to sort and n is an starting point, like 1
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


    !>@brief: Calculate FACE_ELE - the list of elements surrounding an
    !> element and referenced with a face -ve values correspond to surface elements.
    SUBROUTINE CALC_FACE_ELE( FACE_ELE, TOTELE, STOTEL, NFACE, &
        FINELE, COLELE, CV_NLOC, CV_SNLOC, CV_NONODS, CV_NDGLN, CV_SNDGLN, &
        CV_SLOCLIST, X_NLOC, X_NDGLN)
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
    !>@brief:Copies the data from inval to outval safely.
    !>If the sizes are different outval is populated using the first value of inval
    subroutine assign_val(outval,inval)
        implicit none
        real, dimension(:), intent(inout) :: outval
        real, dimension(:), intent(in) :: inval

        if (size(outval)/=size(inval)) then
            outval = inval(1)
        else
            outval = inval
        end if
    end subroutine assign_val


    !!>@brief:This function returns the value at position input_X using
    !>X_points, Y_points to form a linear (size == 2) or quadratic (size == 3) interpolation
    real function table_interpolation(X_points, Y_points, input_X)
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

    !>@brief:Template of csv table
    !>OPTIONAL section (header)
    !>real1,real2,real3,..., size(extra_data)
    !>rows,columns
    !>2,3
    !>Pressure,Saturation
    !>1000,0.9
    !>250,0.5
    !>100,0.1
    subroutine read_csv_table(data_array, path_to_table, extra_data)
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


    !>@brief:This subroutine reads a csv file and returns them in an array
    subroutine extract_strings_from_csv_file(csv_table_strings, path_to_table, Nentries)
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

    !>@brief:Subroutine to print Arrays by (columns,rows)
    !>Matrix = 2D Array
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

!>@brief: Roates a matrix A using the toration matrix R????
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


    !>@brief:This subroutine reads a nastran file that contains the information defining the 1D path of a well
    !>the input relative filepath should include the file format, for example: well.bdf
    subroutine read_nastran_file(filepath, node, edges)
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
            !Nastran files from trellis have assigned a certain space for the different numbers so
            !we extract numbers by position as these are defined
            read(cadena(25:32),*) node(1,i)
            read(cadena(33:40),*) node(2,i)
            read(cadena(41:49),*) node(3,i)
            !Extract as well the number of the node
            read(cadena(9:16),*) conversor(i)
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
! print*, "WELL IN GNUPLOT"
! do j = 1, size(edges,1)
! do i = 1, size(edges,2)
! print *, node(:,edges(j,i))
! end do
! end do
! print *, "-----------------------"
! read*
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


    !---------------------------------------------------------------------------
    !> @author Pablo Salinas
    !> @brief Subroutine that solves the least squares problem |Ax-b|2 using LAPACK and blas
    !> Subroutine tested and compared with Matlab (not recommended changing it since it is a pain!)
    !> Only for serial: The best option is to solve in each processor the optimisation system by performing
    !> A' * A = A' *b; so the system becomes very small as COLUMS <<< ROWS
    !---------------------------------------------------------------------------
    subroutine Least_squares_solver(A, b, rank)
      implicit none
      real, dimension(:,:), intent(inout) :: A !> Input matrix to decompose, returns the Q and R combined
      real, dimension(:,:), intent(inout) :: b !> Input RHS term, returns the X that minimise the system
      integer, intent(inout) :: rank
      !Local variables
      integer :: i, j, k, theta, m, n, nrhs, lda, ldb
      !Parameters for dgegp3 to compute the QR decomposition
      integer :: info!>If info = -i, the i-th parameter had an illegal value
      real, dimension(size(A,2)) :: tau!>Contains scalar factors of the elementary reflectors for the matrix Q.
      real, dimension(3*size(A,1)+1) :: work!>work is a workspace array, its dimension max(1, lwork).
      integer, dimension(size(A,1)) :: jpvt
      real, parameter :: tolerance_rank = 1d-12

      ! interface
      !   !> @brief QR decomposition, returned in A, Q and R mixed, no pivoting!
      !     subroutine dgeqrf(m, n, MAT, lda, tau, work, lwork, info)
      !       implicit none
      !       integer :: m!>Rows of MAT
      !       integer :: n !>Columns of MAT; Constraint: m >= n > = 0.
      !       integer :: lda !>The first dimension of MAT
      !       integer :: lwork!> The size of the work array; 0 == best performance
      !       integer :: info!>If info = -i, the i-th parameter had an illegal value
      !       real, dimension(lda,n) :: MAT!>input/output matrix
      !       real, dimension(N) :: tau!>Contains scalar factors of the elementary reflectors for the matrix Q.
      !       real, dimension(3*n+1) :: work!>work is a workspace array, its dimension max(1, lwork).
      !     end subroutine dgeqrf
      ! end interface

      interface
        !> @brief QR decomposition, returned in A, Q and R mixed, with pivoting! (PREFERRED, obviously!)
          subroutine dgeqp3(m, n, MAT, lda, jpvt, tau, work, lwork, info)
            implicit none
            integer :: m!>Rows of MAT
            integer :: n !>Columns of MAT; Constraint: m >= n > = 0.
            integer :: lda !>The first dimension of MAT
            integer :: lwork!> The size of the work array; 0 == best performance
            integer :: info!>If info = -i, the i-th parameter had an illegal value
            real, dimension(lda,n) :: MAT!>input/output matrix
            real, dimension(N) :: tau!>Contains scalar factors of the elementary reflectors for the matrix Q.
            real, dimension(3*n+1) :: work!>work is a workspace array, its dimension max(1, lwork).
            integer, dimension(n) :: jpvt!>Specifies columns that are not free to move, I guess useful if updating the QR decomposition
          end subroutine dgeqp3
      end interface

      ! interface
      !     !> @brief Interface to Lapack to show a Q matrix computed using dgeqp3
      !     subroutine dorgqr(m, n, k, mat, lda, tau, work, lwork, info)
      !       implicit none
      !       integer :: m,n !>Rows and colums respectively
      !       integer :: lda !>The first dimension of a
      !       integer :: lwork!> The size of the work array; 0 == best performance
      !       integer :: info!>If info = -i, the i-th parameter had an illegal value
      !       integer :: k !>The number of elementary reflectors whose product defines the matrix Q. Constraint 0 ≤k≤m if side='L'; 0 ≤k≤n if side='R'.
      !       real, dimension(m,n) :: MAT!>input/output matrix
      !       real, dimension(N) :: tau!>Contains scalar factors of the elementary reflectors for the matrix Q.
      !       real, dimension(3*n+1) :: work!>work is a workspace array, its dimension max(1, lwork).
      !     end subroutine dorgqr
      ! end interface

      interface
          !> @brief LAPACK subroutine to perform Q times C, Q obtained using dgeqp3
          subroutine DORMQR(side, trans, m, n, k, MAT, lda, tau, c, ldc, work, lwork, info)
            implicit none
            character(len=1) :: side, trans !> Either L or R (Left right); N or T (T == Transpose)
            integer :: m,n !>Rows and colums respectively
            integer :: lda !>The first dimension of a
            integer :: k !>The number of elementary reflectors whose product defines the matrix Q. Constraint 0 ≤k≤m if side='L'; 0 ≤k≤n if side='R'.
            integer :: ldc !>The leading dimension of c. Constraint: ldc≥ max(1, m)
            integer :: lwork!> The size of the work array;For better performance, try using lwork = n*blocksize (if side = 'L') or lwork = m*blocksize (if side = 'R') where blocksize is a machine-dependent value (typically, 16 to 64) required for optimum performance of the blocked algorithm.
            integer :: info!>If info = -i, the i-th parameter had an illegal value
            real, dimension(m,n) :: MAT!>input/output matrix
            real, dimension(ldc,m) :: C !>Overwritten by the product Q*C, QT*C, C*Q, or C*QT (as specified by side and trans).
            real, dimension(N) :: tau !>Contains scalar factors of the elementary reflectors for the matrix Q.
            real, dimension(3*n+1) :: work!>work is a workspace array, its dimension max(1, lwork).
          end subroutine DORMQR
      end interface

      interface
          !> @brief BLAS subroutine to solve the system RX = C
          subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, mat, lda, b, ldb)
            implicit none
            character(len=1) :: side, uplo, transa, diag !> Either L or R (Left right); N or T (T == Transpose)
            integer :: m,n !>Rows and colums respectively
            integer :: lda !>The first dimension of a
            integer :: info!>If info = -i, the i-th parameter had an illegal value
            real :: alpha !> alpha is zero, then a is not referenced and b need not be set before entry.
            integer :: ldb !>The first dimension of b
            real, dimension(m,n) :: MAT!>input/output matrix
            real, dimension(m,n) :: B!>input/output matrix
          end subroutine dtrsm
      end interface

      !Specify dimensions
      m = size(A,1); n = size(A,2); lda = max(1,m)
      !Compute QR decomposition
      jpvt = 0; !Could be an input
      ldb = size(b,1); nrhs = size(b,2); tau =0.
      !Obtain QR decomposition of A using pivoting
      call dgeqp3(m, n, A, lda, jpvt, tau, work, 3* n + 1, info)
      ! call dgeqrf(m, n, A, lda, tau, work, 3* n + 1, info)
      !dorgqr is to obtain Q only
      ! call dorgqr(m, n, min(n,m), A, lda, tau, work, 3* n + 1, info)

      rank = n
      do k = 1, n
        if (abs(A(k,k)) <= tolerance_rank * abs(A(1,1))) then
! print *, "entry that is affecting the result for the least squares",k
          rank = rank - 1
        end if
      end do

      !Perform Q * b
      CALL DORMQR('L','T',m, nrhs, n, A, lda, tau, B, ldb, work, 3* n + 1, info)
      !Now obtain X by solving the system RX = C
      CALL DTRSM('L','U','N','N',min(n,m), nrhs,1., A,lda,B, ldb)
      !Now permute back the results so everything is consistent
      DO J = 1, NRHS
        DO I = 1, N!ldb??
          WORK(JPVT(I)) = B(I,J)
        end do
          B(:,J) = WORK(1:ldb)
      end do
    end subroutine Least_squares_solver


    !> @brief: This subroutine uses python run string to run the python_scalar_diagnostic to read a field
    !> the only difference with the normal approach is that here the Dummy field is used and the returned field is an array.
    !> IMPORTANT: state is used here, NOT packed_state
    subroutine compute_python_scalar_field(state, option_path_python, scalar_result)
      implicit none
      type( state_type ), dimension(:), intent( inout ) :: state
      character( len = * ), intent(in) :: option_path_python
      real, dimension(:), intent(inout) :: scalar_result
      !Local variables
      type (scalar_field), pointer :: sfield
      character( len = python_func_len ) :: pycode

      if (.not.have_option("/material_phase[0]/scalar_field::Dummy")) then
          ewrite(0, *) "ERROR: Trying to compute a python scalar_field without enabling the Dummy field in the first phase."
        stop 657483
      end if


      call python_reset()
      call python_add_state( state(1) )
      sfield => extract_scalar_field(state(1), "Dummy")
      sfield%val = 0.
      call python_run_string("field = state.scalar_fields['Dummy']")
      ! call get_option("/timestepping/current_time", current_time)
      ! write(buffer,*) current_time
      ! call python_run_string("time="//trim(buffer))
      ! call get_option("/timestepping/timestep", dt)
      ! write(buffer,*) dt
      ! call python_run_string("dt="//trim(buffer))
      ! Get the code
      call get_option( trim( option_path_python ) // '/algorithm', pycode )
      ! Run the code
      call python_run_string( trim( pycode ) )
      scalar_result = sfield%val
    end subroutine


    !---------------------------------------------------------------------------
    !> @author Asiri Obeysekara
    !> @brief Subroutines that can initialise, register and start/end a petsc
    !> performance profiling routin. The defauly behaviour is initiliased for
    !> time-loop profiling
    !---------------------------------------------------------------------------
    subroutine petsc_logging(func,stage,ierr,default,push_no,stage_name)

    implicit none
    integer, intent(in) :: func
    integer, parameter  :: N = 8 !this is the default for the time_loop
    integer :: x, i
    integer, optional, intent(in)  :: push_no
    logical, optional, intent(in) :: default
    PetscErrorCode, intent(inout) :: ierr
    PetscLogStage, dimension(0:9)  :: stage

    character(len=*), dimension(1), optional :: stage_name
    character(len=*), dimension(8), parameter ::  stage_name_def &
                                = (/ &
         "PRELIM                ", &
         "FORCE                 ", &
         "SATURATION            ", &
         "TEMP                  ", &
         "COMP                  ", &
         "DT+VTU                ", &
         "ADAPT                 ", &
         "REST                  "/)

#ifdef HAVE_PETSC_DBUG
#if PETSC_VERSION_MINOR<8
print*,"***WARNING: there will be compaitbility issue with your older PETSc &
version and using & profiling, please configure WITHOUT 'petscdebug'"
#else

          !! case 1 - register; case 2 - push ; case 3 - pop
          select case(func)
            case(1)
              !this is to initialise
              if (default) then
              !default is for the main time-loop
                do x=1, N
                  call petsc_log_init(stage_name_def(x),x,stage,ierr)
                end do
              else
                call petsc_log_init(stage_name(1),push_no,stage,ierr)
              end if
            case(2) !!-PUSH
                !this is to initialise
                if (default) then
                  call petsc_log_push(push_no,stage,ierr)
                end if
            case(3) !! - POP
              !this is to initialise
              if (default) then
              !default is for the main time-loop
                  call petsc_log_pop(ierr)
              else
                  call petsc_log_pop(ierr)
              end if
          end select
#endif
#endif
    return
    contains
      !> @brief: This routine registers the stage for PETSC logging
      !> IMPORTANT:
      subroutine petsc_log_init(stage_name,no,stage, ierr)
        implicit none
        PetscErrorCode, intent(inout) :: ierr
        PetscLogStage,dimension(0:9)  :: stage
        integer :: no
        character( len = * ), intent( in ) :: stage_name


        call PetscLogStageRegister(stage_name,stage(no),ierr)

      end subroutine petsc_log_init

      !> @brief: This routine starts the current stage registered
      !> for PETSc profiling
      !> IMPORTANT:
      subroutine petsc_log_push(no,stage, ierr)
        implicit none
        PetscErrorCode, intent(inout) :: ierr
        PetscLogStage,dimension(0:9)  :: stage
        integer :: no

        call PetscLogStagePush(stage(no),ierr)
      end subroutine petsc_log_push

      !> @brief: This routine ends the current stage registered
      !> for PETSc profiling
      !> IMPORTANT:
      subroutine petsc_log_pop(ierr)
        implicit none
        PetscErrorCode, intent(inout) :: ierr

        call PetscLogStagePop(ierr)
      end subroutine petsc_log_pop

    end subroutine petsc_logging

end module multi_tools
