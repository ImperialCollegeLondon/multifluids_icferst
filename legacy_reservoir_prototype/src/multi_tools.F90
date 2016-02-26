  
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


    implicit none

contains

    !sprint_to_do!update to use new memory
    SUBROUTINE SIMPNORM( NORMX, NORMY, NORMZ, D3, &
        SNDGLN, STOTEL, SNLOC, X_NONODS, NONODS, ELE, &
        X, Y, Z, &
        NCOLM, FINDRM, COLM)
        ! Calculate an approx normal (normx,normy,normz)
        IMPLICIT NONE

        REAL, intent( inout ) :: NORMX, NORMY, NORMZ
        LOGICAL, intent( in ) :: D3
        INTEGER, intent( in ) :: STOTEL, SNLOC, X_NONODS, NONODS, ELE, NCOLM
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
        INTEGER :: I
        REAL :: RSUM

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

end module multi_tools


