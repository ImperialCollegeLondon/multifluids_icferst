
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

  module FV_Quadratures
    use fldebug
    use multi_data_types
    use NGauss

  contains

!!$ ===============
!!$ ===============

    SUBROUTINE FV_1D_QUAD( SCVNGI, CV_NLOC, SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
         SCVFENLX, SCVFENLY, SCVFENLZ )
      !
      !     - this subroutine generates the FE basis functions, weights and the
      !     - derivatives of the shape functions for a variety of elements.
      !     - The routine also generates the shape functions and derivatives
      !     - associated with the CV surfaces and also the FV basis functions.
      !     -------------------------------
      !     - date last modified : 24/05/2003
      !     -------------------------------

      IMPLICIT NONE

      INTEGER, intent( in ) :: SCVNGI, CV_NLOC
      REAL, DIMENSION( CV_NLOC, SCVNGI ), intent( inout ) :: SCVFEN, SCVFENSLX, SCVFENSLY, &
           SCVFENLX, SCVFENLY, SCVFENLZ
      REAL, DIMENSION( SCVNGI ), intent( inout ) :: SCVFEWEIGH
      ! Local variables
      INTEGER, PARAMETER :: TWO = 2, THREE = 3, FOUR = 4
      REAL, DIMENSION( : ), allocatable :: LX, CV_NODPOS, WEI
      INTEGER :: NCV_BOU, GPOI, ILOC
      LOGICAL :: DIFF, NDIFF, GETNDP

      !ewrite(3,*) 'cv_nloc:',cv_nloc

      ALLOCATE( LX( SCVNGI ))
      ALLOCATE( CV_NODPOS( CV_NLOC ))
      ALLOCATE( WEI( CV_NLOC ))

      SCVFEWEIGH = 1.0
      IF(SCVNGI==2) THEN
         NCV_BOU = 2
         LX(1) = -1.0
         LX(2) = +1.0
      ELSE IF(SCVNGI==3) THEN
         NCV_BOU = 3
         LX(1) = -1.0
         LX(2) =  0.0
         LX(3) = +1.0
      ELSE IF(SCVNGI==4) THEN
         NCV_BOU = 4
         LX(1) = -1.0
         LX(2) = -0.5
         LX(3) = +0.5
         LX(4) = +1.0
      ELSE IF(SCVNGI==5) THEN
         NCV_BOU = 5
         LX(1) = -1.0
         LX(2) = -1.0 + 1./3.
         LX(3) = +0.0
         LX(4) = +1.0 - 1./3.
         LX(5) = +1.0
      ELSE
         FLAbort(" Wrong number of computed surface quadrature points for CV ")
      ENDIF

      DIFF=.TRUE.
      NDIFF=.FALSE.

      GETNDP=.TRUE.
      !      Compute standard Gauss quadrature. weits and points
      CALL LAGROT(WEI,CV_NODPOS,CV_NLOC,GETNDP)
      !EWRITE(3,*)'CV_NODPOS:',CV_NODPOS

      Loop_P2: DO GPOI = 1, NCV_BOU

         Loop_ILX2: DO  ILOC = 1, CV_NLOC
            SCVFEN( ILOC, GPOI )  = LAGRAN( NDIFF, LX( GPOI ), ILOC, CV_NLOC, CV_NODPOS )
            SCVFENSLX( ILOC, GPOI ) = 1.0
            SCVFENSLY( ILOC, GPOI ) = 0.0
            SCVFENLX( ILOC, GPOI ) = LAGRAN( DIFF, LX( GPOI ), ILOC, CV_NLOC, CV_NODPOS )
            SCVFENLY( ILOC, GPOI ) = 0.0
            SCVFENLZ( ILOC, GPOI ) = 0.0
         END DO Loop_ILX2

      end do Loop_P2

      DEALLOCATE( LX )
      DEALLOCATE( CV_NODPOS )
      DEALLOCATE( WEI )

      return

    END SUBROUTINE FV_1D_QUAD


!!$ ===============
!!$ ===============

    SUBROUTINE FVQUAD( NGI, NLOC, SVNGI, &
         M, SVN, SVNLX,                  &
         SVWEIGH )

!!$   This routine generates the shape functions associated
!!$       with the FV's i.e. their surfaces and volume shape
!!$       functions and derivatives. The surface shape functions
!!$       are the values of the FE volume shape functions evaluated
!!$       on the surfaces of the CV's.

      IMPLICIT NONE
      INTEGER, intent(in) :: NGI, NLOC, SVNGI
      REAL, dimension(NLOC,NGI), intent(inout) :: M
      REAL, dimension(NLOC,SVNGI), intent(inout) :: SVN, SVNLX
      REAL, dimension(SVNGI), intent(inout) :: SVWEIGH
!!$ Local variables:
!!$ Note that NCOORD is the number of co-ordinates.
!!$           NFACE is the number of surfaces internal to an element
      integer, parameter ::  NCOORD = 2, NFACE = 4, ENLOC = 4, FNGI = 1, FNLOC = 2
      INTEGER :: ILOC, JLOC, IFACE, ICOORD, GI, GJ

      real, dimension(:), allocatable :: pos, dpdxi, dpdeta, xigp, xip, xi, eta
      real, dimension(:, :, :), allocatable :: corn

      allocate( POS(NCOORD), DPDXI(NCOORD), DPDETA(NCOORD), XIGP(FNGI), XIP(FNLOC), &
           XI(ENLOC), ETA(ENLOC), CORN(NFACE,NCOORD,FNLOC) )

!!$ FV basis functions for use in calculating the integral of the
!!$     shape functions over a subcell of a CV.

      M = 0.0


!!$ Note that CORN(IFACE,ICOORD,JLOC) tabulates the positions
!!$      of the corners (vertices) of the faces of a subcell in volume
!!$      co-ordinates. Note also that IFACE ranges from 1-4, ICOORD
!!$      ranges from 1-2 and JLOC ranges from 1-2. IFACE signifies
!!$      the particular subcell face, ICOORD signfies the co-ordinates
!!$      of the vertices of the face in (XI,ETA) co-ordinates
!!$      and JLOC signifies the vertex of the face of the subcell
!!$      which are straight lines.

      !     - face 1
      CORN(1,1,1) = -1.0   ;   CORN(1,2,1) =  0.0
      CORN(1,1,2) =  0.0   ;   CORN(1,2,2) =  0.0

      !     - face 2
      CORN(2,1,1) =  0.0   ;   CORN(2,2,1) =  0.0
      CORN(2,1,2) =  1.0   ;   CORN(2,2,2) =  0.0

      !     - face 3
      CORN(3,1,1) =  0.0   ;   CORN(3,2,1) = -1.0
      CORN(3,1,2) =  0.0   ;   CORN(3,2,2) =  0.0

      !     - face 4
      CORN(4,1,1) =  0.0   ;   CORN(4,2,1) =  0.0
      CORN(4,1,2) =  0.0   ;   CORN(4,2,2) =  1.0


!!$   Define the Gaussian integration points in XIP space.
      XIGP(1)  = 0.0

!!$   Define positions of vertices of line element
!!$          associated with a subcell face in XIP
      XIP(1) = -1.0   ;   XIP(2) =  1.0

!!$   Define positions of vertices of quadrilateral element in
!!$          (XI,ETA) space.
      XI(1)  = -1.0   ;   ETA(1) = -1.0
      XI(2)  =  1.0   ;   ETA(2) = -1.0
      XI(3)  = -1.0   ;   ETA(3) =  1.0
      XI(4)  =  1.0   ;   ETA(4) =  1.0

!!$   Generate values of the FE basis functions at the quadrature
!!$          points on the faces of the subcells.
      Loop_ILOC:do ILOC = 1, NLOC

         Loop_IFACE: do IFACE = 1, NFACE

            Loop_GJ: do GJ = 1, FNGI

               GI = IFACE

               Loop_ICOORD: do ICOORD = 1, NCOORD

                  POS(ICOORD) = 0.0 ; DPDXI(ICOORD) = 0.0 ; DPDETA(ICOORD) = 0.0

                  Loop_JLOC: do JLOC = 1, FNLOC
                     POS(ICOORD) = POS(ICOORD) + CORN(IFACE,ICOORD,JLOC) * &
                          0.5*(1.+XIP(JLOC)*XIGP(GJ))

                     DPDXI(ICOORD) = DPDXI(ICOORD) + CORN(IFACE,ICOORD,JLOC) * &
                          0.5*XIP(JLOC)

                  END DO Loop_JLOC

               END DO Loop_ICOORD

               SVN(ILOC,GI) = 0.25*( 1.0 + XI(ILOC)*POS(1) ) *( 1.0 + ETA(ILOC)*POS(2) )

               SVNLX(ILOC,GI) = 0.25*XI(ILOC)*DPDXI(1) *( 1.0 + ETA(ILOC)*POS(2) ) +  &
                    0.25*( 1.0 + XI(ILOC)*POS(1) ) *ETA(ILOC)*DPDXI(2)

            END DO Loop_GJ

         END DO Loop_IFACE

      END DO Loop_ILOC

!!$ set the weights for the surface integration
      SVWEIGH = 2.

      deallocate( POS, DPDXI, DPDETA, XIGP, XIP, XI, ETA, CORN )

      return
    END SUBROUTINE FVQUAD

!!$ ===============
!!$ ===============

    SUBROUTINE FVQQUAD( NGI, NLOC, SVNGI,&
         M, SVN, SVNLX, &
         SVWEIGH )

!!$   This routine generates the shape functions associated
!!$       with the FV's i.e. their surfaces and volume shape
!!$       functions and derivatives. The surface shape functions
!!$       are the values of the FE volume shape functions evaluated
!!$       on the surfaces of the CV's.

      IMPLICIT NONE
      INTEGER, intent(in) :: NGI, NLOC, SVNGI
      REAL, dimension(NLOC,NGI), intent(inout) :: M
      REAL, dimension(NLOC,SVNGI), intent(inout) :: SVN, SVNLX
      REAL, dimension(SVNGI), intent(inout) :: SVWEIGH
!!$ Local variables:
!!$ Note that NCOORD is the number of co-ordinates.
!!$           NFACE is the number of surfaces internal to an element
      integer, parameter ::  NCOORD = 2, NFACE = 12, FNGI = 2, FNLOC = 2
      INTEGER :: ILOC, JLOC, IFACE, ICOORD, COUNT, GI, GJ, I, J
      real :: xi, eta

      real, dimension(:), allocatable :: pos, dpdxi, dpdeta, xigp, xip, lijxi, &
           lijeta, dlijxidxi, dlijetadxi
      real, dimension(:, :, :), allocatable :: corn

      allocate( POS(NCOORD), DPDXI(NCOORD), DPDETA(NCOORD), XIGP(FNGI), XIP(FNLOC), &
           CORN(NFACE,NCOORD,FNLOC), LIJXI(3), LIJETA(3), DLIJXIDXI(3), DLIJETADXI(3) )

!!$ FV basis functions for use in calculating the integral of the
!!$     shape functions over a subcell of a CV.

      M = 0.0
      COUNT = 1

      do ILOC = 1, NLOC
         do GI = COUNT, COUNT + 3
            M(ILOC,GI) = 1.0
         END DO
         COUNT = COUNT + 4
      END DO

!!$ Note that CORN(IFACE,ICOORD,JLOC) tabulates the positions
!!$      of the corners (vertices) of the faces of a subcell in volume
!!$      co-ordinates. Note also that IFACE ranges from 1-4, ICOORD
!!$      ranges from 1-2 and JLOC ranges from 1-2. IFACE signifies
!!$      the particular subcell face, ICOORD signfies the co-ordinates
!!$      of the vertices of the face in (XI,ETA) co-ordinates
!!$      and JLOC signifies the vertex of the face of the subcell
!!$      which are straight lines.

      !     - face 1
      CORN(1,1,1) = -0.5   ;   CORN(1,2,1) = -1.0
      CORN(1,1,2) = -0.5   ;   CORN(1,2,2) = -0.5

      !     - face 2
      CORN(2,1,1) =  0.5   ;   CORN(2,2,1) = -1.0
      CORN(2,1,2) =  0.5   ;   CORN(2,2,2) = -0.5

      !     - face 3
      CORN(3,1,1) = -1.0   ;   CORN(3,2,1) = -0.5
      CORN(3,1,2) = -0.5   ;   CORN(3,2,2) = -0.5

      !     - face 4
      CORN(4,1,1) = -0.5   ;  CORN(4,2,1) = -0.5
      CORN(4,1,2) =  0.5   ;  CORN(4,2,2) = -0.5

      !     - face 5
      CORN(5,1,1) =  0.5   ;   CORN(5,2,1) = -0.5
      CORN(5,1,2) =  1.0   ;   CORN(5,2,2) = -0.5

      !     - face 6
      CORN(6,1,1) = -0.5   ;   CORN(6,2,1) = -0.5
      CORN(6,1,2) = -0.5   ;   CORN(6,2,2) =  0.5

      !     - face 7
      CORN(7,1,1) =  0.5   ;   CORN(7,2,1) = -0.5
      CORN(7,1,2) =  0.5   ;   CORN(7,2,2) =  0.5

      !     - face 8
      CORN(8,1,1) = -1.0   ;   CORN(8,2,1) =  0.5
      CORN(8,1,2) = -0.5   ;   CORN(8,2,2) =  0.5

      !     - face 9
      CORN(9,1,1) = -0.5   ;   CORN(9,2,1) =  0.5
      CORN(9,1,2) =  0.5   ;   CORN(9,2,2) =  0.5

      !     - face 10
      CORN(10,1,1) =  0.5   ;  CORN(10,2,1) =  0.5
      CORN(10,1,2) =  1.0   ;  CORN(10,2,2) =  0.5

      !     - face 11
      CORN(11,1,1) = -0.5   ;  CORN(11,2,1) =  0.5
      CORN(11,1,2) = -0.5   ;  CORN(11,2,2) =  1.0

      !     - face 12
      CORN(12,1,1) =  0.5   ;  CORN(12,2,1) =  0.5
      CORN(12,1,2) =  0.5   ;  CORN(12,2,2) =  1.0

!!$   Define the Gaussian integration points in XIP space.
      XIGP(1)  = -1./SQRT(3.0) ; XIGP(2)  =  1./SQRT(3.0)

!!$   Define positions of vertices of line element associated with a 
!!$          subcell face in XIP
      XIP(1) = -1.0 ; XIP(2) =  1.0

!!$   Generate values of the FE basis functions at the quadrature
!!$          points on the faces of the subcells.
      Loop_IFACE: do IFACE = 1, NFACE

         Loop_GJ: do GJ = 1, FNGI
            GI = (IFACE-1)*FNGI+GJ

            Loop_ICOORD: do ICOORD = 1, NCOORD
               POS(ICOORD) = 0.0 ; DPDXI(ICOORD) = 0.0 ; DPDETA(ICOORD) = 0.0

               Loop_JLOC: do JLOC = 1, FNLOC
                  POS(ICOORD) = POS(ICOORD) + CORN(IFACE,ICOORD,JLOC) * &
                       0.5*(1.+XIP(JLOC)*XIGP(GJ))
                  DPDXI(ICOORD) = DPDXI(ICOORD) + CORN(IFACE,ICOORD,JLOC) * &
                       0.5*XIP(JLOC)

               END DO Loop_JLOC
               !     
            END DO Loop_ICOORD

            XI  = POS(1) ; ETA = POS(2)
!!$
            LIJXI(1) = 0.5*XI*(XI-1.0) ; LIJXI(2) = 1.0-XI*XI ; LIJXI(3)   =  0.5*XI*(XI+1.0)
!!$
            LIJETA(1) = 0.5*ETA*(ETA-1.0) ; LIJETA(2) = 1.0-ETA*ETA ; LIJETA(3) = 0.5*ETA*(ETA+1.0)
!!$

            DLIJXIDXI(1) = 0.5*(2.0*XI-1.0)*DPDXI(1) ; DLIJXIDXI(2) = -2.0*XI*DPDXI(1) ; &
                 DLIJXIDXI(3) = 0.5*(2.0*XI+1.0)*DPDXI(1)
!!$

            DLIJETADXI(1) = 0.5*(2.0*ETA-1.0)*DPDXI(2) ; DLIJETADXI(2) = -2.0*ETA*DPDXI(2) ; &
                 DLIJETADXI(3) = 0.5*(2.0*ETA+1.0)*DPDXI(2)
!!$

            Loop_I: do I = 1, 3
               Loop_J: do J = 1, 3
                  ILOC = I+(J-1)*3
                  SVN(ILOC,GI) = LIJXI(I)*LIJETA(J)
                  SVNLX(ILOC,GI) = LIJXI(I)*DLIJETADXI(J) + LIJETA(J)*DLIJXIDXI(I)
               END DO Loop_J
            END DO Loop_I

         END DO Loop_GJ

      END DO Loop_IFACE

!!$   Set the weights for the surface integration
      SVWEIGH = 1.

      deallocate( POS, DPDXI, DPDETA, XIGP, XIP, CORN, LIJXI, LIJETA, DLIJXIDXI, DLIJETADXI )

      return

    END SUBROUTINE FVQQUAD

!!$ ===============
!!$ ===============



    SUBROUTINE FVHEX( NGI, NLOC, SVNGI, &
         M, SVN, SVNLX,                 &
         SVNLY, SVWEIGH )

!!$   This routine generates the shape functions associated
!!$       with the FV's i.e. their surfaces and volume shape
!!$       functions and derivatives.

      IMPLICIT NONE
      INTEGER, intent(in) :: NGI, NLOC, SVNGI
      REAL, dimension(NLOC,NGI), intent(inout) :: M
      REAL, dimension(NLOC,SVNGI), intent(inout) :: SVN, SVNLX, SVNLY
      REAL, dimension(SVNGI), intent(inout) :: SVWEIGH
!!$ Local variables:
!!$ Note that NCOORD is the number of co-ordinates.
!!$           NFACE is the number of surfaces internal to an element
      integer, parameter :: NCOORD = 3, NFACE = 12, ENLOC = 8, FNGI = 1, FNLOC = 4
      INTEGER :: ILOC, JLOC, IFACE, ICOORD, COUNT, GI, GJ

      real, dimension(:), allocatable :: pos, dpdxi, dpdeta, xigp, etagp, xip, etap, xi, eta, zeta
      real, dimension(:, :, :), allocatable :: corn

      allocate( POS(NCOORD), DPDXI(NCOORD), DPDETA(NCOORD), XIGP(FNGI), ETAGP(FNGI), XIP(FNLOC), &
           ETAP(FNLOC), XI(ENLOC), ETA(ENLOC), ZETA(ENLOC), CORN(NFACE,NCOORD,FNLOC) )

!!$ FV basis functions for use in calculating the integral of the
!!$     shape functions over a subcell of a CV.

      M = 0.0 
      COUNT = 1


!!$ Note that CORN(IFACE,ICOORD,JLOC) tabulates the positions
!!$      of the corners (vertices) of the faces of a subcell in volume
!!$      co-ordinates. Note also that IFACE ranges from 1-4, ICOORD
!!$      ranges from 1-2 and JLOC ranges from 1-2. IFACE signifies
!!$      the particular subcell face, ICOORD signfies the co-ordinates
!!$      of the vertices of the face in (XI,ETA) co-ordinates
!!$      and JLOC signifies the vertex of the face of the subcell
!!$      which are straight lines.

      !     - face 1
      CORN(1,1,1) = -1.0  ;  CORN(1,2,1) =  0.0  ;  CORN(1,3,1) = -1.0
      CORN(1,1,2) =  0.0  ;  CORN(1,2,2) =  0.0  ;  CORN(1,3,2) = -1.0
      CORN(1,1,3) = -1.0  ;  CORN(1,2,3) =  0.0  ;  CORN(1,3,3) =  0.0
      CORN(1,1,4) =  0.0  ;  CORN(1,2,4) =  0.0  ;  CORN(1,3,4) =  0.0

      !     - face 2
      CORN(2,1,1) =  0.0  ;  CORN(2,2,1) =  0.0  ;  CORN(2,3,1) = -1.0
      CORN(2,1,2) =  1.0  ;  CORN(2,2,2) =  0.0  ;  CORN(2,3,2) = -1.0
      CORN(2,1,3) =  0.0  ;  CORN(2,2,3) =  0.0  ;  CORN(2,3,3) =  0.0
      CORN(2,1,4) =  1.0  ;  CORN(2,2,4) =  0.0  ;  CORN(2,3,4) =  0.0

      !     - face 3
      CORN(3,1,1) =  0.0  ;  CORN(3,2,1) = -1.0  ;  CORN(3,3,1) = -1.0
      CORN(3,1,2) =  0.0  ;  CORN(3,2,2) =  0.0  ;  CORN(3,3,2) = -1.0
      CORN(3,1,3) =  0.0  ;  CORN(3,2,3) = -1.0  ;  CORN(3,3,3) =  0.0
      CORN(3,1,4) =  0.0  ;  CORN(3,2,4) =  0.0  ;  CORN(3,3,4) =  0.0

      !     - face 4
      CORN(4,1,1) =  0.0  ;  CORN(4,2,1) =  0.0  ;  CORN(4,3,1) = -1.0
      CORN(4,1,2) =  0.0  ;  CORN(4,2,2) =  1.0  ;  CORN(4,3,2) = -1.0
      CORN(4,1,3) =  0.0  ;  CORN(4,2,3) =  0.0  ;  CORN(4,3,3) =  0.0
      CORN(4,1,4) =  0.0  ;  CORN(4,2,4) =  1.0  ;  CORN(4,3,4) =  0.0

      !     - face 5
      CORN(5,1,1) = -1.0  ;  CORN(5,2,1) = -1.0  ;  CORN(5,3,1) =  0.0
      CORN(5,1,2) =  0.0  ;  CORN(5,2,2) = -1.0  ;  CORN(5,3,2) =  0.0
      CORN(5,1,3) = -1.0  ;  CORN(5,2,3) =  0.0  ;  CORN(5,3,3) =  0.0
      CORN(5,1,4) =  0.0  ;  CORN(5,2,4) =  0.0  ;  CORN(5,3,4) =  0.0

      !     - face 6
      CORN(6,1,1) =  0.0  ;  CORN(6,2,1) = -1.0  ;  CORN(6,3,1) =  0.0
      CORN(6,1,2) =  1.0  ;  CORN(6,2,2) = -1.0  ;  CORN(6,3,2) =  0.0
      CORN(6,1,3) =  0.0  ;  CORN(6,2,3) =  0.0  ;  CORN(6,3,3) =  0.0
      CORN(6,1,4) =  1.0  ;  CORN(6,2,4) =  0.0  ;  CORN(6,3,4) =  0.0

      !     - face 7
      CORN(7,1,1) = -1.0  ;  CORN(7,2,1) =  0.0  ;  CORN(7,3,1) =  0.0
      CORN(7,1,2) =  0.0  ;  CORN(7,2,2) =  0.0  ;  CORN(7,3,2) =  0.0
      CORN(7,1,3) = -1.0  ;  CORN(7,2,3) =  1.0  ;  CORN(7,3,3) =  0.0
      CORN(7,1,4) =  0.0  ;  CORN(7,2,4) =  1.0  ;  CORN(7,3,4) =  0.0

      !     - face 8
      CORN(8,1,1) =  0.0  ;  CORN(8,2,1) =  0.0  ;  CORN(8,3,1) =  0.0
      CORN(8,1,2) =  1.0  ;  CORN(8,2,2) =  0.0  ;  CORN(8,3,2) =  0.0
      CORN(8,1,3) =  0.0  ;  CORN(8,2,3) =  1.0  ;  CORN(8,3,3) =  0.0
      CORN(8,1,4) =  1.0  ;  CORN(8,2,4) =  1.0  ;  CORN(8,3,4) =  0.0

      !     - face 9
      CORN(9,1,1) = -1.0  ;  CORN(9,2,1) =  0.0  ;  CORN(9,3,1) =  0.0
      CORN(9,1,2) =  0.0  ;  CORN(9,2,2) =  0.0  ;  CORN(9,3,2) =  0.0
      CORN(9,1,3) = -1.0  ;  CORN(9,2,3) =  0.0  ;  CORN(9,3,3) =  1.0
      CORN(9,1,4) =  0.0  ;  CORN(9,2,4) =  0.0  ;  CORN(9,3,4) =  1.0

      !     - face 10
      CORN(10,1,1) = 0.0  ;  CORN(10,2,1) = 0.0  ;  CORN(10,3,1) = 0.0
      CORN(10,1,2) = 1.0  ;  CORN(10,2,2) = 0.0  ;  CORN(10,3,2) = 0.0
      CORN(10,1,3) = 0.0  ;  CORN(10,2,3) = 0.0  ;  CORN(10,3,3) = 1.0
      CORN(10,1,4) = 1.0  ;  CORN(10,2,4) = 0.0  ;  CORN(10,3,4) = 1.0

      !     - face 11
      CORN(11,1,1) =  0.0 ;  CORN(11,2,1) = -1.0 ;  CORN(11,3,1) = 0.0
      CORN(11,1,2) =  0.0 ;  CORN(11,2,2) =  0.0 ;  CORN(11,3,2) = 0.0
      CORN(11,1,3) =  0.0 ;  CORN(11,2,3) = -1.0 ;  CORN(11,3,3) = 1.0
      CORN(11,1,4) =  0.0 ;  CORN(11,2,4) =  0.0 ;  CORN(11,3,4) = 1.0

      !     - face 12
      CORN(12,1,1) = 0.0  ;  CORN(12,2,1) = 0.0  ;  CORN(12,3,1) = 0.0
      CORN(12,1,2) = 0.0  ;  CORN(12,2,2) = 1.0  ;  CORN(12,3,2) = 0.0
      CORN(12,1,3) = 0.0  ;  CORN(12,2,3) = 0.0  ;  CORN(12,3,3) = 1.0
      CORN(12,1,4) = 0.0  ;  CORN(12,2,4) = 1.0  ;  CORN(12,3,4) = 1.0

!!$   Define the Gaussian integration points in (XIP,ETAP) space.
      XIGP(1)  = 0.0  ;  ETAGP(1) = 0.0

!!$   Define positions of vertices of quadrilateral element associated with a 
!!$          subcell face in (XIP,ETAP) space.
      XIP(1)  = -1.0  ;  ETAP(1) = -1.0
      XIP(2)  =  1.0  ;  ETAP(2) = -1.0
      XIP(3)  = -1.0  ;  ETAP(3) =  1.0
      XIP(4)  =  1.0  ;  ETAP(4) =  1.0

!!$   Define positions of vertices of hexahedral element in 
!!$          (XI,ETA,ZETA) space.
      XI(1)   = -1.0  ;  ETA(1)  = -1.0  ;  ZETA(1) = -1.0
      XI(2)   =  1.0  ;  ETA(2)  = -1.0  ;  ZETA(2) = -1.0
      XI(3)   = -1.0  ;  ETA(3)  =  1.0  ;  ZETA(3) = -1.0
      XI(4)   =  1.0  ;  ETA(4)  =  1.0  ;  ZETA(4) = -1.0
      XI(5)   = -1.0  ;  ETA(5)  = -1.0  ;  ZETA(5) =  1.0
      XI(6)   =  1.0  ;  ETA(6)  = -1.0  ;  ZETA(6) =  1.0
      XI(7)   = -1.0  ;  ETA(7)  =  1.0  ;  ZETA(7) =  1.0
      XI(8)   =  1.0  ;  ETA(8)  =  1.0  ;  ZETA(8) =  1.0

!!$   Generate values of the FE basis functions at the quadrature
!!$          points on the faces of the subcells.
      Loop_ILOC: do ILOC = 1, NLOC

         Loop_IFACE: do IFACE = 1, NFACE

            Loop_GJ:do GJ = 1, FNGI
               GI = IFACE

               Loop_ICOORD: do ICOORD = 1, NCOORD
                  POS(ICOORD) = 0.0 ; DPDXI(ICOORD) = 0.0 ; DPDETA(ICOORD) = 0.0

                  Loop_JLOC: do JLOC = 1, FNLOC
                     POS(ICOORD) = POS(ICOORD) + CORN(IFACE,ICOORD,JLOC) * 0.25*(1.+XIP(JLOC)*XIGP(GJ)) * &
                          (1.+ETAP(JLOC)*ETAGP(GJ))
                     DPDXI(ICOORD) = DPDXI(ICOORD) + CORN(IFACE,ICOORD,JLOC) * 0.25*XIP(JLOC) * &
                          (1.+ETAP(JLOC)*ETAGP(GJ))
                     DPDETA(ICOORD) = DPDETA(ICOORD) + CORN(IFACE,ICOORD,JLOC) * 0.25*(1.+XIP(JLOC)*XIGP(GJ)) * &
                          ETAP(JLOC)
                  END DO Loop_JLOC

               END DO Loop_ICOORD

               SVN(ILOC,GI) = 0.125*( 1.0 + XI(ILOC)*POS(1) ) * ( 1.0 + ETA(ILOC)*POS(2) ) * &
                    ( 1.0 + ZETA(ILOC)*POS(3) )
               SVNLX(ILOC,GI) = 0.125*XI(ILOC)*DPDXI(1) * ( 1.0 + ETA(ILOC)*POS(2) ) * &
                    ( 1.0 + ZETA(ILOC)*POS(3) ) + 0.125*( 1.0 + XI(ILOC)*POS(1) ) * ETA(ILOC)*DPDXI(2) * &
                    ( 1.0 + ZETA(ILOC)*POS(3) ) + 0.125*( 1.0 + XI(ILOC)*POS(1) ) * ( 1.0 + ETA(ILOC)*POS(2) ) * &
                    ZETA(ILOC)*DPDXI(3)
               SVNLY(ILOC,GI) = 0.125*XI(ILOC)*DPDETA(1) * ( 1.0 + ETA(ILOC)*POS(2) ) *&
                    ( 1.0 + ZETA(ILOC)*POS(3) ) + 0.125*( 1.0 + XI(ILOC)*POS(1) ) * ETA(ILOC)*DPDETA(2) * &
                    ( 1.0 + ZETA(ILOC)*POS(3) ) + 0.125*( 1.0 + XI(ILOC)*POS(1) ) * ( 1.0 + ETA(ILOC)*POS(2) ) * &
                    ZETA(ILOC)*DPDETA(3)

            END DO Loop_GJ

         END DO Loop_IFACE

      END DO Loop_ILOC

!!$   Set the weights for the surface integration. Note that this is 2*2.0 = 4.0 for the 
!!$         weight because it is a product of two 1-D weights for integration over a 2-D quadrilateral.
      SVWEIGH = 4.

      deallocate( POS, DPDXI, DPDETA, XIGP, ETAGP, XIP, ETAP, XI, ETA, ZETA, CORN )

      return
    END SUBROUTINE FVHEX


!!$ ===============
!!$ ===============

    SUBROUTINE FVQHEX( NGI, NLOC, SVNGI,  &
         M,  SVN, SVNLX,                  &
         SVNLY, SVWEIGH )

!!$   This routine generates the shape functions associated
!!$       with the FV's i.e. their surfaces and volume shape
!!$       functions and derivatives. The surface shape functions
!!$       are the values of the FE volume shape functions evaluated
!!$       on the surfaces of the CV's.

      IMPLICIT NONE
      INTEGER, intent(in) :: NGI, NLOC, SVNGI
      REAL, dimension(NLOC,NGI), intent(inout) :: M
      REAL, dimension(NLOC,SVNGI), intent(inout) :: SVN, SVNLX, SVNLY
      REAL, dimension(SVNGI), intent(inout) :: SVWEIGH
!!$ Local variables:
!!$ Note that NCOORD is the number of co-ordinates.
!!$           NFACE is the number of surfaces internal to an element
      integer, parameter :: NCOORD = 3, NFACE = 54, FNGI = 4, FNLOC = 4
      INTEGER :: ILOC, JLOC, IFACE, ICOORD, COUNT, GI, GJ, i, j, k

      real, dimension(:), allocatable :: pos, dpdxi, dpdeta, xigp, etagp, xip, etap, lijxi, &
           lijeta, lijzeta, dlijxidxi, dlijetadxi, dlijzetadxi, dlijxideta, dlijetadeta, dlijzetadeta
      real, dimension(:, :, :), allocatable :: corn
      real :: xi, eta, zeta

      allocate( POS(NCOORD), DPDXI(NCOORD), DPDETA(NCOORD), XIGP(FNGI), ETAGP(FNGI), XIP(FNLOC),      &
           ETAP(FNLOC), LIJXI(3), LIJETA(3), LIJZETA(3), DLIJXIDXI(3), DLIJETADXI(3), DLIJZETADXI(3), &
           DLIJXIDETA(3), DLIJETADETA(3), DLIJZETADETA(3), CORN(NFACE,NCOORD,FNLOC) )

!!$ FV basis functions for use in calculating the integral of the
!!$     shape functions over a subcell of a CV.

      M = 0.0
      COUNT = 1

!!$ Note that CORN(IFACE,ICOORD,JLOC) tabulates the positions
!!$      of the corners (vertices) of the faces of a subcell in volume
!!$      co-ordinates. Note also that IFACE ranges from 1-4, ICOORD
!!$      ranges from 1-2 and JLOC ranges from 1-2. IFACE signifies
!!$      the particular subcell face, ICOORD signfies the co-ordinates
!!$      of the vertices of the face in (XI,ETA) co-ordinates
!!$      and JLOC signifies the vertex of the face of the subcell
!!$      which are straight lines.

      !     - face 1
      CORN(1,1,1) = -0.5   ;   CORN(1,2,1) = -1.0   ;   CORN(1,3,1) = -1.0
      CORN(1,1,2) = -0.5   ;   CORN(1,2,2) = -0.5   ;   CORN(1,3,2) = -1.0
      CORN(1,1,3) = -0.5   ;   CORN(1,2,3) = -1.0   ;   CORN(1,3,3) = -0.5
      CORN(1,1,4) = -0.5   ;   CORN(1,2,4) = -0.5   ;   CORN(1,3,4) = -0.5

      !     - face 2
      CORN(2,1,1) =  0.5   ;   CORN(2,2,1) = -1.0   ;   CORN(2,3,1) = -1.0
      CORN(2,1,2) =  0.5   ;   CORN(2,2,2) = -0.5   ;   CORN(2,3,2) = -1.0
      CORN(2,1,3) =  0.5   ;   CORN(2,2,3) = -1.0   ;   CORN(2,3,3) = -0.5
      CORN(2,1,4) =  0.5   ;   CORN(2,2,4) = -0.5   ;   CORN(2,3,4) = -0.5

      !     - face 3
      CORN(3,1,1) = -1.0   ;   CORN(3,2,1) = -0.5   ;   CORN(3,3,1) = -1.0
      CORN(3,1,2) = -0.5   ;   CORN(3,2,2) = -0.5   ;   CORN(3,3,2) = -1.0
      CORN(3,1,3) = -1.0   ;   CORN(3,2,3) = -0.5   ;   CORN(3,3,3) = -0.5
      CORN(3,1,4) = -0.5   ;   CORN(3,2,4) = -0.5   ;   CORN(3,3,4) = -0.5

      !     - face 4
      CORN(4,1,1) = -0.5   ;   CORN(4,2,1) = -0.5   ;   CORN(4,3,1) = -1.0
      CORN(4,1,2) =  0.5   ;   CORN(4,2,2) = -0.5   ;   CORN(4,3,2) = -1.0
      CORN(4,1,3) = -0.5   ;   CORN(4,2,3) = -0.5   ;   CORN(4,3,3) = -0.5
      CORN(4,1,4) =  0.5   ;   CORN(4,2,4) = -0.5   ;   CORN(4,3,4) = -0.5

      !     - face 5
      CORN(5,1,1) =  0.5   ;   CORN(5,2,1) = -0.5   ;   CORN(5,3,1) = -1.0
      CORN(5,1,2) =  1.0   ;   CORN(5,2,2) = -0.5   ;   CORN(5,3,2) = -1.0
      CORN(5,1,3) =  0.5   ;   CORN(5,2,3) = -0.5   ;   CORN(5,3,3) = -0.5
      CORN(5,1,4) =  1.0   ;   CORN(5,2,4) = -0.5   ;   CORN(5,3,4) = -0.5

      !     - face 6
      CORN(6,1,1) = -0.5   ;   CORN(6,2,1) = -0.5   ;   CORN(6,3,1) = -1.0
      CORN(6,1,2) = -0.5   ;   CORN(6,2,2) =  0.5   ;   CORN(6,3,2) = -1.0
      CORN(6,1,3) = -0.5   ;   CORN(6,2,3) = -0.5   ;   CORN(6,3,3) = -0.5
      CORN(6,1,4) = -0.5   ;   CORN(6,2,4) =  0.5   ;   CORN(6,3,4) = -0.5

      !     - face 7
      CORN(7,1,1) =  0.5   ;   CORN(7,2,1) = -0.5   ;   CORN(7,3,1) = -1.0
      CORN(7,1,2) =  0.5   ;   CORN(7,2,2) =  0.5   ;   CORN(7,3,2) = -1.0
      CORN(7,1,3) =  0.5   ;   CORN(7,2,3) = -0.5   ;   CORN(7,3,3) = -0.5
      CORN(7,1,4) =  0.5   ;   CORN(7,2,4) =  0.5   ;   CORN(7,3,4) = -0.5

      !     - face 8
      CORN(8,1,1) = -1.0   ;   CORN(8,2,1) =  0.5   ;   CORN(8,3,1) = -1.0
      CORN(8,1,2) = -0.5   ;   CORN(8,2,2) =  0.5   ;   CORN(8,3,2) = -1.0
      CORN(8,1,3) = -1.0   ;   CORN(8,2,3) =  0.5   ;   CORN(8,3,3) = -0.5
      CORN(8,1,4) = -0.5   ;   CORN(8,2,4) =  0.5   ;   CORN(8,3,4) = -0.5

      !     - face 9
      CORN(9,1,1) = -0.5   ;   CORN(9,2,1) =  0.5   ;   CORN(9,3,1) = -1.0
      CORN(9,1,2) =  0.5   ;   CORN(9,2,2) =  0.5   ;   CORN(9,3,2) = -1.0
      CORN(9,1,3) = -0.5   ;   CORN(9,2,3) =  0.5   ;   CORN(9,3,3) = -0.5
      CORN(9,1,4) =  0.5   ;   CORN(9,2,4) =  0.5   ;   CORN(9,3,4) = -0.5

      !     - face 10
      CORN(10,1,1) =  0.5  ;   CORN(10,2,1) =  0.5  ;   CORN(10,3,1) = -1.0
      CORN(10,1,2) =  1.0  ;   CORN(10,2,2) =  0.5  ;   CORN(10,3,2) = -1.0
      CORN(10,1,3) =  0.5  ;   CORN(10,2,3) =  0.5  ;   CORN(10,3,3) = -0.5
      CORN(10,1,4) =  1.0  ;   CORN(10,2,4) =  0.5  ;   CORN(10,3,4) = -0.5

      !     - face 11
      CORN(11,1,1) = -0.5  ;   CORN(11,2,1) =  0.5  ;   CORN(11,3,1) = -1.0
      CORN(11,1,2) = -0.5  ;   CORN(11,2,2) =  1.0  ;   CORN(11,3,2) = -1.0
      CORN(11,1,3) = -0.5  ;   CORN(11,2,3) =  0.5  ;   CORN(11,3,3) = -0.5
      CORN(11,1,4) = -0.5  ;   CORN(11,2,4) =  1.0  ;   CORN(11,3,4) = -0.5

      !     - face 12
      CORN(12,1,1) =  0.5  ;   CORN(12,2,1) =  0.5  ;   CORN(12,3,1) = -1.0
      CORN(12,1,2) =  0.5  ;   CORN(12,2,2) =  1.0  ;   CORN(12,3,2) = -1.0
      CORN(12,1,3) =  0.5  ;   CORN(12,2,3) =  0.5  ;   CORN(12,3,3) = -0.5
      CORN(12,1,4) =  0.5  ;   CORN(12,2,4) =  1.0  ;   CORN(12,3,4) = -0.5

      !     - face 13
      CORN(13,1,1) = -1.0  ;   CORN(13,2,1) = -1.0  ;   CORN(13,3,1) = -0.5
      CORN(13,1,2) = -0.5  ;   CORN(13,2,2) = -1.0  ;   CORN(13,3,2) = -0.5
      CORN(13,1,3) = -1.0  ;   CORN(13,2,3) = -0.5  ;   CORN(13,3,3) = -0.5
      CORN(13,1,4) = -0.5  ;   CORN(13,2,4) = -0.5  ;   CORN(13,3,4) = -0.5

      !     - face 14
      CORN(14,1,1) = -0.5  ;   CORN(14,2,1) = -1.0  ;   CORN(14,3,1) = -0.5
      CORN(14,1,2) =  0.5  ;   CORN(14,2,2) = -1.0  ;   CORN(14,3,2) = -0.5
      CORN(14,1,3) = -0.5  ;   CORN(14,2,3) = -0.5  ;   CORN(14,3,3) = -0.5
      CORN(14,1,4) =  0.5  ;   CORN(14,2,4) = -0.5  ;   CORN(14,3,4) = -0.5

      !     - face 15
      CORN(15,1,1) =  0.5  ;   CORN(15,2,1) = -1.0  ;   CORN(15,3,1) = -0.5
      CORN(15,1,2) =  1.0  ;   CORN(15,2,2) = -1.0  ;   CORN(15,3,2) = -0.5
      CORN(15,1,3) =  0.5  ;   CORN(15,2,3) = -0.5  ;   CORN(15,3,3) = -0.5
      CORN(15,1,4) =  1.0  ;   CORN(15,2,4) = -0.5  ;   CORN(15,3,4) = -0.5

      !     - face 16
      CORN(16,1,1) = -1.0  ;   CORN(16,2,1) = -0.5  ;   CORN(16,3,1) = -0.5
      CORN(16,1,2) = -0.5  ;   CORN(16,2,2) = -0.5  ;   CORN(16,3,2) = -0.5
      CORN(16,1,3) = -1.0  ;   CORN(16,2,3) =  0.5  ;   CORN(16,3,3) = -0.5
      CORN(16,1,4) = -0.5  ;   CORN(16,2,4) =  0.5  ;   CORN(16,3,4) = -0.5

      !     - face 17
      CORN(17,1,1) = -0.5  ;   CORN(17,2,1) = -0.5  ;   CORN(17,3,1) = -0.5
      CORN(17,1,2) =  0.5  ;   CORN(17,2,2) = -0.5  ;   CORN(17,3,2) = -0.5
      CORN(17,1,3) = -0.5  ;   CORN(17,2,3) =  0.5  ;   CORN(17,3,3) = -0.5
      CORN(17,1,4) =  0.5  ;   CORN(17,2,4) =  0.5  ;   CORN(17,3,4) = -0.5

      !     - face 18
      CORN(18,1,1) =  0.5  ;   CORN(18,2,1) = -0.5  ;   CORN(18,3,1) = -0.5
      CORN(18,1,2) =  1.0  ;   CORN(18,2,2) = -0.5  ;   CORN(18,3,2) = -0.5
      CORN(18,1,3) =  0.5  ;   CORN(18,2,3) =  0.5  ;   CORN(18,3,3) = -0.5
      CORN(18,1,4) =  1.0  ;   CORN(18,2,4) =  0.5  ;   CORN(18,3,4) = -0.5

      !     - face 19
      CORN(19,1,1) = -1.0  ;   CORN(19,2,1) =  0.5  ;   CORN(19,3,1) = -0.5
      CORN(19,1,2) = -0.5  ;   CORN(19,2,2) =  0.5  ;   CORN(19,3,2) = -0.5
      CORN(19,1,3) = -1.0  ;   CORN(19,2,3) =  1.0  ;   CORN(19,3,3) = -0.5
      CORN(19,1,4) = -0.5  ;   CORN(19,2,4) =  1.0  ;   CORN(19,3,4) = -0.5

      !     - face 20
      CORN(20,1,1) = -0.5  ;   CORN(20,2,1) =  0.5  ;   CORN(20,3,1) = -0.5
      CORN(20,1,2) =  0.5  ;   CORN(20,2,2) =  0.5  ;   CORN(20,3,2) = -0.5
      CORN(20,1,3) = -0.5  ;   CORN(20,2,3) =  1.0  ;   CORN(20,3,3) = -0.5
      CORN(20,1,4) =  0.5  ;   CORN(20,2,4) =  1.0  ;   CORN(20,3,4) = -0.5

      !     - face 21
      CORN(21,1,1) =  0.5  ;   CORN(21,2,1) =  0.5  ;   CORN(21,3,1) = -0.5
      CORN(21,1,2) =  1.0  ;   CORN(21,2,2) =  0.5  ;   CORN(21,3,2) = -0.5
      CORN(21,1,3) =  0.5  ;   CORN(21,2,3) =  1.0  ;   CORN(21,3,3) = -0.5
      CORN(21,1,4) =  1.0  ;   CORN(21,2,4) =  1.0  ;   CORN(21,3,4) = -0.5

      !     - face 22
      CORN(22,1,1) = -0.5  ;   CORN(22,2,1) = -1.0  ;   CORN(22,3,1) = -0.5
      CORN(22,1,2) = -0.5  ;   CORN(22,2,2) = -0.5  ;   CORN(22,3,2) = -0.5
      CORN(22,1,3) = -0.5  ;   CORN(22,2,3) = -1.0  ;   CORN(22,3,3) =  0.5
      CORN(22,1,4) = -0.5  ;   CORN(22,2,4) = -0.5  ;   CORN(22,3,4) =  0.5

      !     - face 23
      CORN(23,1,1) =  0.5  ;   CORN(23,2,1) = -1.0  ;   CORN(23,3,1) = -0.5
      CORN(23,1,2) =  0.5  ;   CORN(23,2,2) = -0.5  ;   CORN(23,3,2) = -0.5
      CORN(23,1,3) =  0.5  ;   CORN(23,2,3) = -1.0  ;   CORN(23,3,3) =  0.5
      CORN(23,1,4) =  0.5  ;   CORN(23,2,4) = -0.5  ;   CORN(23,3,4) =  0.5

      !     - face 24
      CORN(24,1,1) = -1.0  ;   CORN(24,2,1) = -0.5  ;   CORN(24,3,1) = -0.5
      CORN(24,1,2) = -0.5  ;   CORN(24,2,2) = -0.5  ;   CORN(24,3,2) = -0.5
      CORN(24,1,3) = -1.0  ;   CORN(24,2,3) = -0.5  ;   CORN(24,3,3) =  0.5
      CORN(24,1,4) = -0.5  ;   CORN(24,2,4) = -0.5  ;   CORN(24,3,4) =  0.5

      !     - face 25
      CORN(25,1,1) = -0.5  ;   CORN(25,2,1) = -0.5  ;   CORN(25,3,1) = -0.5
      CORN(25,1,2) =  0.5  ;   CORN(25,2,2) = -0.5  ;   CORN(25,3,2) = -0.5
      CORN(25,1,3) = -0.5  ;   CORN(25,2,3) = -0.5  ;   CORN(25,3,3) =  0.5
      CORN(25,1,4) =  0.5  ;   CORN(25,2,4) = -0.5  ;   CORN(25,3,4) =  0.5

      !     - face 26
      CORN(26,1,1) =  0.5  ;   CORN(26,2,1) = -0.5  ;   CORN(26,3,1) = -0.5
      CORN(26,1,2) =  1.0  ;   CORN(26,2,2) = -0.5  ;   CORN(26,3,2) = -0.5
      CORN(26,1,3) =  0.5  ;   CORN(26,2,3) = -0.5  ;   CORN(26,3,3) =  0.5
      CORN(26,1,4) =  1.0  ;   CORN(26,2,4) = -0.5  ;   CORN(26,3,4) =  0.5

      !     - face 27
      CORN(27,1,1) = -0.5  ;   CORN(27,2,1) = -0.5  ;   CORN(27,3,1) = -0.5
      CORN(27,1,2) = -0.5  ;   CORN(27,2,2) =  0.5  ;   CORN(27,3,2) = -0.5
      CORN(27,1,3) = -0.5  ;   CORN(27,2,3) = -0.5  ;   CORN(27,3,3) =  0.5
      CORN(27,1,4) = -0.5  ;   CORN(27,2,4) =  0.5  ;   CORN(27,3,4) =  0.5

      !     - face 28
      CORN(28,1,1) =  0.5  ;   CORN(28,2,1) = -0.5  ;   CORN(28,3,1) = -0.5
      CORN(28,1,2) =  0.5  ;   CORN(28,2,2) =  0.5  ;   CORN(28,3,2) = -0.5
      CORN(28,1,3) =  0.5  ;   CORN(28,2,3) = -0.5  ;   CORN(28,3,3) =  0.5
      CORN(28,1,4) =  0.5  ;   CORN(28,2,4) =  0.5  ;   CORN(28,3,4) =  0.5

      !     - face 29
      CORN(29,1,1) = -1.0  ;   CORN(29,2,1) =  0.5  ;   CORN(29,3,1) = -0.5
      CORN(29,1,2) = -0.5  ;   CORN(29,2,2) =  0.5  ;   CORN(29,3,2) = -0.5
      CORN(29,1,3) = -1.0  ;   CORN(29,2,3) =  0.5  ;   CORN(29,3,3) =  0.5
      CORN(29,1,4) = -0.5  ;   CORN(29,2,4) =  0.5  ;   CORN(29,3,4) =  0.5

      !     - face 30
      CORN(30,1,1) = -0.5  ;   CORN(30,2,1) =  0.5  ;   CORN(30,3,1) = -0.5
      CORN(30,1,2) =  0.5  ;   CORN(30,2,2) =  0.5  ;   CORN(30,3,2) = -0.5
      CORN(30,1,3) = -0.5  ;   CORN(30,2,3) =  0.5  ;   CORN(30,3,3) =  0.5
      CORN(30,1,4) =  0.5  ;   CORN(30,2,4) =  0.5  ;   CORN(30,3,4) =  0.5

      !     - face 31
      CORN(31,1,1) =  0.5  ;   CORN(31,2,1) =  0.5  ;   CORN(31,3,1) = -0.5
      CORN(31,1,2) =  1.0  ;   CORN(31,2,2) =  0.5  ;   CORN(31,3,2) = -0.5
      CORN(31,1,3) =  0.5  ;   CORN(31,2,3) =  0.5  ;   CORN(31,3,3) =  0.5
      CORN(31,1,4) =  1.0  ;   CORN(31,2,4) =  0.5  ;   CORN(31,3,4) =  0.5

      !     - face 32
      CORN(32,1,1) = -0.5  ;   CORN(32,2,1) =  0.5  ;  CORN(32,3,1) = -0.5
      CORN(32,1,2) = -0.5  ;   CORN(32,2,2) =  1.0  ;  CORN(32,3,2) = -0.5
      CORN(32,1,3) = -0.5  ;   CORN(32,2,3) =  0.5  ;  CORN(32,3,3) =  0.5
      CORN(32,1,4) = -0.5  ;   CORN(32,2,4) =  1.0  ;  CORN(32,3,4) =  0.5

      !     - face 33
      CORN(33,1,1) =  0.5  ;   CORN(33,2,1) =  0.5  ;  CORN(33,3,1) = -0.5
      CORN(33,1,2) =  0.5  ;   CORN(33,2,2) =  1.0  ;  CORN(33,3,2) = -0.5
      CORN(33,1,3) =  0.5  ;   CORN(33,2,3) =  0.5  ;  CORN(33,3,3) =  0.5
      CORN(33,1,4) =  0.5  ;   CORN(33,2,4) =  1.0  ;  CORN(33,3,4) =  0.5

      !     - face 34
      CORN(34,1,1) = -1.0  ;   CORN(34,2,1) = -1.0  ;  CORN(34,3,1) =  0.5
      CORN(34,1,2) = -0.5  ;   CORN(34,2,2) = -1.0  ;  CORN(34,3,2) =  0.5
      CORN(34,1,3) = -1.0  ;   CORN(34,2,3) = -0.5  ;  CORN(34,3,3) =  0.5
      CORN(34,1,4) = -0.5  ;   CORN(34,2,4) = -0.5  ;  CORN(34,3,4) =  0.5

      !     - face 35
      CORN(35,1,1) = -0.5  ;   CORN(35,2,1) = -1.0  ;  CORN(35,3,1) =  0.5
      CORN(35,1,2) =  0.5  ;   CORN(35,2,2) = -1.0  ;  CORN(35,3,2) =  0.5
      CORN(35,1,3) = -0.5  ;   CORN(35,2,3) = -0.5  ;  CORN(35,3,3) =  0.5
      CORN(35,1,4) =  0.5  ;   CORN(35,2,4) = -0.5  ;  CORN(35,3,4) =  0.5

      !     - face 36
      CORN(36,1,1) =  0.5  ;   CORN(36,2,1) = -1.0  ;  CORN(36,3,1) =  0.5
      CORN(36,1,2) =  1.0  ;   CORN(36,2,2) = -1.0  ;  CORN(36,3,2) =  0.5
      CORN(36,1,3) =  0.5  ;   CORN(36,2,3) = -0.5  ;  CORN(36,3,3) =  0.5
      CORN(36,1,4) =  1.0  ;   CORN(36,2,4) = -0.5  ;  CORN(36,3,4) =  0.5

      !     - face 37
      CORN(37,1,1) = -1.0  ;   CORN(37,2,1) = -0.5  ;   CORN(37,3,1) =  0.5
      CORN(37,1,2) = -0.5  ;   CORN(37,2,2) = -0.5  ;   CORN(37,3,2) =  0.5
      CORN(37,1,3) = -1.0  ;   CORN(37,2,3) =  0.5  ;   CORN(37,3,3) =  0.5
      CORN(37,1,4) = -0.5  ;   CORN(37,2,4) =  0.5  ;   CORN(37,3,4) =  0.5

      !     - face 38
      CORN(38,1,1) = -0.5  ;   CORN(38,2,1) = -0.5  ;   CORN(38,3,1) =  0.5
      CORN(38,1,2) =  0.5  ;   CORN(38,2,2) = -0.5  ;   CORN(38,3,2) =  0.5
      CORN(38,1,3) = -0.5  ;   CORN(38,2,3) =  0.5  ;   CORN(38,3,3) =  0.5
      CORN(38,1,4) =  0.5  ;   CORN(38,2,4) =  0.5  ;   CORN(38,3,4) =  0.5

      !     - face 39
      CORN(39,1,1) =  0.5  ;   CORN(39,2,1) = -0.5  ;   CORN(39,3,1) =  0.5
      CORN(39,1,2) =  1.0  ;   CORN(39,2,2) = -0.5  ;   CORN(39,3,2) =  0.5
      CORN(39,1,3) =  0.5  ;   CORN(39,2,3) =  0.5  ;   CORN(39,3,3) =  0.5
      CORN(39,1,4) =  1.0  ;   CORN(39,2,4) =  0.5  ;   CORN(39,3,4) =  0.5

      !     - face 40
      CORN(40,1,1) = -1.0  ;   CORN(40,2,1) =  0.5  ;   CORN(40,3,1) =  0.5
      CORN(40,1,2) = -0.5  ;   CORN(40,2,2) =  0.5  ;   CORN(40,3,2) =  0.5
      CORN(40,1,3) = -1.0  ;   CORN(40,2,3) =  1.0  ;   CORN(40,3,3) =  0.5
      CORN(40,1,4) = -0.5  ;   CORN(40,2,4) =  1.0  ;   CORN(40,3,4) =  0.5

      !     - face 41
      CORN(41,1,1) = -0.5  ;   CORN(41,2,1) =  0.5  ;   CORN(41,3,1) =  0.5
      CORN(41,1,2) =  0.5  ;   CORN(41,2,2) =  0.5  ;   CORN(41,3,2) =  0.5
      CORN(41,1,3) = -0.5  ;   CORN(41,2,3) =  1.0  ;   CORN(41,3,3) =  0.5
      CORN(41,1,4) =  0.5  ;   CORN(41,2,4) =  1.0  ;   CORN(41,3,4) =  0.5

      !     - face 42
      CORN(42,1,1) =  0.5  ;   CORN(42,2,1) =  0.5  ;   CORN(42,3,1) =  0.5
      CORN(42,1,2) =  1.0  ;   CORN(42,2,2) =  0.5  ;   CORN(42,3,2) =  0.5
      CORN(42,1,3) =  0.5  ;   CORN(42,2,3) =  1.0  ;   CORN(42,3,3) =  0.5
      CORN(42,1,4) =  1.0  ;   CORN(42,2,4) =  1.0  ;   CORN(42,3,4) =  0.5

      !     - face 43
      CORN(43,1,1) = -0.5  ;   CORN(43,2,1) = -1.0  ;   CORN(43,3,1) =  0.5
      CORN(43,1,2) = -0.5  ;   CORN(43,2,2) = -0.5  ;   CORN(43,3,2) =  0.5
      CORN(43,1,3) = -0.5  ;   CORN(43,2,3) = -1.0  ;   CORN(43,3,3) =  1.0
      CORN(43,1,4) = -0.5  ;   CORN(43,2,4) = -0.5  ;   CORN(43,3,4) =  1.0

      !     - face 44
      CORN(44,1,1) =  0.5  ;   CORN(44,2,1) = -1.0  ;   CORN(44,3,1) =  0.5
      CORN(44,1,2) =  0.5  ;   CORN(44,2,2) = -0.5  ;   CORN(44,3,2) =  0.5
      CORN(44,1,3) =  0.5  ;   CORN(44,2,3) = -1.0  ;   CORN(44,3,3) =  1.0
      CORN(44,1,4) =  0.5  ;   CORN(44,2,4) = -0.5  ;   CORN(44,3,4) =  1.0

      !     - face 45
      CORN(45,1,1) = -1.0  ;   CORN(45,2,1) = -0.5  ;   CORN(45,3,1) =  0.5
      CORN(45,1,2) = -0.5  ;   CORN(45,2,2) = -0.5  ;   CORN(45,3,2) =  0.5
      CORN(45,1,3) = -1.0  ;   CORN(45,2,3) = -0.5  ;   CORN(45,3,3) =  1.0
      CORN(45,1,4) = -0.5  ;   CORN(45,2,4) = -0.5  ;   CORN(45,3,4) =  1.0

      !     - face 46
      CORN(46,1,1) = -0.5  ;   CORN(46,2,1) = -0.5  ;   CORN(46,3,1) =  0.5
      CORN(46,1,2) =  0.5  ;   CORN(46,2,2) = -0.5  ;   CORN(46,3,2) =  0.5
      CORN(46,1,3) = -0.5  ;   CORN(46,2,3) = -0.5  ;   CORN(46,3,3) =  1.0
      CORN(46,1,4) =  0.5  ;   CORN(46,2,4) = -0.5  ;   CORN(46,3,4) =  1.0

      !     - face 47
      CORN(47,1,1) =  0.5  ;   CORN(47,2,1) = -0.5  ;   CORN(47,3,1) =  0.5
      CORN(47,1,2) =  1.0  ;   CORN(47,2,2) = -0.5  ;   CORN(47,3,2) =  0.5
      CORN(47,1,3) =  0.5  ;   CORN(47,2,3) = -0.5  ;   CORN(47,3,3) =  1.0
      CORN(47,1,4) =  1.0  ;   CORN(47,2,4) = -0.5  ;   CORN(47,3,4) =  1.0

      !     - face 48
      CORN(48,1,1) = -0.5  ;   CORN(48,2,1) = -0.5  ;   CORN(48,3,1) =  0.5
      CORN(48,1,2) = -0.5  ;   CORN(48,2,2) =  0.5  ;   CORN(48,3,2) =  0.5
      CORN(48,1,3) = -0.5  ;   CORN(48,2,3) = -0.5  ;   CORN(48,3,3) =  1.0
      CORN(48,1,4) = -0.5  ;   CORN(48,2,4) =  0.5  ;   CORN(48,3,4) =  1.0

      !     - face 49
      CORN(49,1,1) =  0.5  ;   CORN(49,2,1) = -0.5  ;   CORN(49,3,1) =  0.5
      CORN(49,1,2) =  0.5  ;   CORN(49,2,2) =  0.5  ;   CORN(49,3,2) =  0.5
      CORN(49,1,3) =  0.5  ;   CORN(49,2,3) = -0.5  ;   CORN(49,3,3) =  1.0
      CORN(49,1,4) =  0.5  ;   CORN(49,2,4) =  0.5  ;   CORN(49,3,4) =  1.0

      !     - face 50
      CORN(50,1,1) = -1.0  ;   CORN(50,2,1) =  0.5  ;   CORN(50,3,1) =  0.5
      CORN(50,1,2) = -0.5  ;   CORN(50,2,2) =  0.5  ;   CORN(50,3,2) =  0.5
      CORN(50,1,3) = -1.0  ;   CORN(50,2,3) =  0.5  ;   CORN(50,3,3) =  1.0
      CORN(50,1,4) = -0.5  ;   CORN(50,2,4) =  0.5  ;   CORN(50,3,4) =  1.0

      !     - face 51
      CORN(51,1,1) = -0.5  ;   CORN(51,2,1) =  0.5  ;   CORN(51,3,1) =  0.5
      CORN(51,1,2) =  0.5  ;   CORN(51,2,2) =  0.5  ;   CORN(51,3,2) =  0.5
      CORN(51,1,3) = -0.5  ;   CORN(51,2,3) =  0.5  ;   CORN(51,3,3) =  1.0
      CORN(51,1,4) =  0.5  ;   CORN(51,2,4) =  0.5  ;   CORN(51,3,4) =  1.0

      !     - face 52
      CORN(52,1,1) =  0.5  ;   CORN(52,2,1) =  0.5  ;   CORN(52,3,1) =  0.5
      CORN(52,1,2) =  1.0  ;   CORN(52,2,2) =  0.5  ;   CORN(52,3,2) =  0.5
      CORN(52,1,3) =  0.5  ;   CORN(52,2,3) =  0.5  ;   CORN(52,3,3) =  1.0
      CORN(52,1,4) =  1.0  ;   CORN(52,2,4) =  0.5  ;   CORN(52,3,4) =  1.0

      !     - face 53
      CORN(53,1,1) = -0.5  ;   CORN(53,2,1) =  0.5  ;   CORN(53,3,1) =  0.5
      CORN(53,1,2) = -0.5  ;   CORN(53,2,2) =  1.0  ;   CORN(53,3,2) =  0.5
      CORN(53,1,3) = -0.5  ;   CORN(53,2,3) =  0.5  ;   CORN(53,3,3) =  1.0
      CORN(53,1,4) = -0.5  ;   CORN(53,2,4) =  1.0  ;   CORN(53,3,4) =  1.0

      !     - face 54
      CORN(54,1,1) =  0.5  ;   CORN(54,2,1) =  0.5  ;   CORN(54,3,1) =  0.5
      CORN(54,1,2) =  0.5  ;   CORN(54,2,2) =  1.0  ;   CORN(54,3,2) =  0.5
      CORN(54,1,3) =  0.5  ;   CORN(54,2,3) =  0.5  ;   CORN(54,3,3) =  1.0
      CORN(54,1,4) =  0.5  ;   CORN(54,2,4) =  1.0  ;   CORN(54,3,4) =  1.0

!!$   Define the Gaussian integration points in XIP space.
      XIPGP(1)  = -1.0/SQRT(3.0)  ;   ETAPGP(1) = -1.0/SQRT(3.0)
      XIPGP(2)  =  1.0/SQRT(3.0)  ;   ETAPGP(2) = -1.0/SQRT(3.0)
      XIPGP(3)  = -1.0/SQRT(3.0)  ;   ETAPGP(3) =  1.0/SQRT(3.0)
      XIPGP(4)  =  1.0/SQRT(3.0)  ;   ETAPGP(4) =  1.0/SQRT(3.0)

!!$   Define positions of vertices of quadrilateral element associated with a 
!!$          subcell face in XIP
      XIP(1)  = -1.0  ;   ETAP(1) = -1.0
      XIP(2)  =  1.0  ;   ETAP(2) = -1.0
      XIP(3)  = -1.0  ;   ETAP(3) =  1.0
      XIP(4)  =  1.0  ;   ETAP(4) =  1.0

!!$   Generate values of the FE basis functions at the quadrature
!!$          points on the faces of the subcells.
      Loop_IFACE: do IFACE = 1, NFACE

         Loop_GJ: do GJ = 1, FNGI
            GI = (IFACE-1)*FNGI+GJ

            Loop_ICOORD: do ICOORD = 1, NCOORD
               POS(ICOORD) = 0.0 ; DPDXI(ICOORD) = 0.0 ; DPDETA(ICOORD) = 0.0

               Loop_JLOC: do JLOC = 1, FNLOC
                  POS(ICOORD) = POS(ICOORD) + CORN(IFACE,ICOORD,JLOC) * 0.25*(1.+XIP(JLOC)*XIPGP(GJ)) * &
                       (1.+ETAP(JLOC)*ETAPGP(GJ))
                  DPDXI(ICOORD) = DPDXI(ICOORD) + CORN(IFACE,ICOORD,JLOC) * 0.25*XIP(JLOC) * &
                       (1.+ETAP(JLOC)*ETAPGP(GJ))
                  DPDETA(ICOORD) = DPDETA(ICOORD) + CORN(IFACE,ICOORD,JLOC) * 0.25*(1.+XIP(JLOC)*XIPGP(GJ)) * &
                       ETAP(JLOC)

               END DO Loop_JLOC

            END DO Loop_ICOORD

            XI = POS(1) ; ETA = POS(2) ; ZETA = POS(3)
!!$
            LIJXI(1) = 0.5*XI*(XI-1.0) ; LIJXI(2) = 1.0-XI*XI ; LIJXI(3) = 0.5*XI*(XI+1.0)
!!$
            LIJETA(1) = 0.5*ETA*(ETA-1.0) ; LIJETA(2) = 1.0-ETA*ETA ; LIJETA(3) = 0.5*ETA*(ETA+1.0)
!!$
            LIJZETA(1) = 0.5*ZETA*(ZETA-1.0) ; LIJZETA(2) = 1.0-ZETA*ZETA ; LIJZETA(3) = 0.5*ZETA*(ZETA+1.0)
!!$
            DLIJXIDXI(1) = 0.5*(2.0*XI-1.0)*DPDXI(1) ; DLIJXIDXI(2) = -2.0*XI*DPDXI(1) ; &
                 DLIJXIDXI(3) = 0.5*(2.0*XI+1.0)*DPDXI(1)
!!$     
            DLIJETADXI(1) = 0.5*(2.0*ETA-1.0)*DPDXI(2) ; DLIJETADXI(2) = -2.0*ETA*DPDXI(2) ; &
                 DLIJETADXI(3)  =  0.5*(2.0*ETA+1.0)*DPDXI(2)
!!$          
            DLIJZETADXI(1) = 0.5*(2.0*ZETA-1.0)*DPDXI(3) ; DLIJZETADXI(2) = -2.0*ZETA*DPDXI(3) ; &
                 DLIJZETADXI(3) = 0.5*(2.0*ZETA+1.0)*DPDXI(3)
!!$               
            DLIJXIDETA(1) = 0.5*(2.0*XI-1.0)*DPDETA(1) ; DLIJXIDETA(2) = -2.0*XI*DPDETA(1) ; &
                 DLIJXIDETA(3) = 0.5*(2.0*XI+1.0)*DPDETA(1)
!!$               
            DLIJETADETA(1) = 0.5*(2.0*ETA-1.0)*DPDETA(2) ; DLIJETADETA(2) = -2.0*ETA*DPDETA(2) ; &
                 DLIJETADETA(3) = 0.5*(2.0*ETA+1.0)*DPDETA(2)
!!$              
            DLIJZETADETA(1) = 0.5*(2.0*ZETA-1.0)*DPDETA(3) ; DLIJZETADETA(2) = -2.0*ZETA*DPDETA(3) ; &
                 DLIJZETADETA(3) = 0.5*(2.0*ZETA+1.0)*DPDETA(3)

            Loop_I: do I = 1, 3
               Loop_J: do J = 1, 3
                  Loop_K: do K = 1, 3
                     ILOC = I+(J-1)*3+(K-1)*9
                     SVN(ILOC,GI) = LIJXI(I)*LIJETA(J)*LIJZETA(K)
                     SVNLX(ILOC,GI) = LIJXI(I)*LIJETA(J)*DLIJZETADXI(K) + LIJXI(I)*DLIJETADXI(J)*LIJZETA(K) + &
                          DLIJXIDXI(I)*LIJETA(J)*LIJZETA(K)
                     SVNLY(ILOC,GI) = LIJXI(I)*LIJETA(J)*DLIJZETADETA(K) + LIJXI(I)*DLIJETADETA(J)*LIJZETA(K) + &
                          DLIJXIDETA(I)*LIJETA(J)*LIJZETA(K)
                  END DO Loop_K
               END DO Loop_J
            END DO Loop_I

         END DO Loop_GJ

      END DO Loop_IFACE

!!$   Set the weights for the surface integration
      SVWEIGH = 1.

      deallocate( POS, DPDXI, DPDETA, XIGP, ETAGP, XIP, ETAP, LIJXI, LIJETA, LIJZETA, DLIJXIDXI, DLIJETADXI, &
           DLIJZETADXI, DLIJXIDETA, DLIJETADETA, DLIJZETADETA, CORN )

      return

    END SUBROUTINE FVQHEX


!!$ ===============
!!$ ===============


  end module FV_Quadratures
