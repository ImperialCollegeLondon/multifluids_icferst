
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
      !
      !     - face 1
      !
      CORN(1,1,1) = -1.0
      CORN(1,2,1) =  0.0
      !
      CORN(1,1,2) =  0.0
      CORN(1,2,2) =  0.0
      !
      !     - face 2
      !
      CORN(2,1,1) =  0.0
      CORN(2,2,1) =  0.0
      !
      CORN(2,1,2) =  1.0
      CORN(2,2,2) =  0.0
      !
      !     - face 3
      !
      CORN(3,1,1) =  0.0
      CORN(3,2,1) = -1.0
      !
      CORN(3,1,2) =  0.0
      CORN(3,2,2) =  0.0
      !
      !     - face 4
      !
      CORN(4,1,1) =  0.0
      CORN(4,2,1) =  0.0
      !
      CORN(4,1,2) =  0.0
      CORN(4,2,2) =  1.0
      !
      !     - define the Gaussian integration points in XIP space.
      !
      XIGP(1)  = 0.0
      !
      !     - define positions of vertices of line element
      !     - associated with a subcell face in XIP
      !
      XIP(1) = -1.0
      XIP(2) =  1.0
      !
      !     - define positions of vertices of quadrilateral element in
      !     - (XI,ETA) space.
      !
      XI(1)  = -1.0
      ETA(1) = -1.0
      !
      XI(2)  =  1.0
      ETA(2) = -1.0
      !
      XI(3)  = -1.0
      ETA(3) =  1.0
      !
      XI(4)  =  1.0
      ETA(4) =  1.0
      !
      !     - generate values of the FE basis functions at the quadrature
      !     - points on the faces of the subcells.
      !
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
      !
      !     - Note that CORN(IFACE,ICOORD,JLOC) tabulates the positions
      !     - of the corners (vertices) of the faces of a subcell in volume
      !     - co-ordinates. Note also that IFACE ranges from 1-4, ICOORD
      !     - ranges from 1-2 and JLOC ranges from 1-2. IFACE signifies
      !     - the particular subcell face, ICOORD signfies the co-ordinates
      !     - of the vertices of the face in (XI,ETA) co-ordinates
      !     - and JLOC signifies the vertex of the face of the subcell
      !     - which are straight lines.
      !
      !     - face 1
      !
      CORN(1,1,1) = -0.5
      CORN(1,2,1) = -1.0
      !
      CORN(1,1,2) = -0.5
      CORN(1,2,2) = -0.5
      !
      !     - face 2
      !
      CORN(2,1,1) =  0.5
      CORN(2,2,1) = -1.0
      !
      CORN(2,1,2) =  0.5
      CORN(2,2,2) = -0.5
      !
      !     - face 3
      !
      CORN(3,1,1) = -1.0
      CORN(3,2,1) = -0.5
      !
      CORN(3,1,2) = -0.5
      CORN(3,2,2) = -0.5
      !
      !     - face 4
      !
      CORN(4,1,1) = -0.5
      CORN(4,2,1) = -0.5
      !
      CORN(4,1,2) =  0.5
      CORN(4,2,2) = -0.5
      !
      !     - face 5
      !
      CORN(5,1,1) =  0.5
      CORN(5,2,1) = -0.5
      !
      CORN(5,1,2) =  1.0
      CORN(5,2,2) = -0.5
      !
      !     - face 6
      !
      CORN(6,1,1) = -0.5
      CORN(6,2,1) = -0.5
      !
      CORN(6,1,2) = -0.5
      CORN(6,2,2) =  0.5
      !
      !     - face 7
      !
      CORN(7,1,1) =  0.5
      CORN(7,2,1) = -0.5
      !
      CORN(7,1,2) =  0.5
      CORN(7,2,2) =  0.5
      !
      !     - face 8
      !
      CORN(8,1,1) = -1.0
      CORN(8,2,1) =  0.5
      !
      CORN(8,1,2) = -0.5
      CORN(8,2,2) =  0.5
      !
      !     - face 9
      !
      CORN(9,1,1) = -0.5
      CORN(9,2,1) =  0.5
      !
      CORN(9,1,2) =  0.5
      CORN(9,2,2) =  0.5
      !
      !     - face 10
      !
      CORN(10,1,1) =  0.5
      CORN(10,2,1) =  0.5
      !
      CORN(10,1,2) =  1.0
      CORN(10,2,2) =  0.5
      !
      !     - face 11
      !
      CORN(11,1,1) = -0.5
      CORN(11,2,1) =  0.5
      !
      CORN(11,1,2) = -0.5
      CORN(11,2,2) =  1.0
      !
      !     - face 12
      !
      CORN(12,1,1) =  0.5
      CORN(12,2,1) =  0.5
      !
      CORN(12,1,2) =  0.5
      CORN(12,2,2) =  1.0
      !
      !     - define the Gaussian integration points in XIP space.
      !
      XIGP(1)  = -1/SQRT(3.0)
      XIGP(2)  =  1/SQRT(3.0)
      !
      !     - define positions of vertices of line element
      !     - associated with a subcell face in XIP
      !
      XIP(1) = -1.0
      XIP(2) =  1.0
      !
      !     - generate values of the FE basis functions at the quadrature
      !     - points on the faces of the subcells.
      !
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
                 
            LIJXI(1)   =  0.5*XI*(XI-1.0)
            LIJXI(2)   =  1.0-XI*XI
            LIJXI(3)   =  0.5*XI*(XI+1.0)
                 
            LIJETA(1)  =  0.5*ETA*(ETA-1.0)
            LIJETA(2)  =  1.0-ETA*ETA
            LIJETA(3)  =  0.5*ETA*(ETA+1.0)
                 
            DLIJXIDXI(1) =  0.5*(2.0*XI-1.0)*DPDXI(1)
            DLIJXIDXI(2) = -2.0*XI*DPDXI(1)
            DLIJXIDXI(3) =  0.5*(2.0*XI+1.0)*DPDXI(1)
                 
            DLIJETADXI(1) =  0.5*(2.0*ETA-1.0)*DPDXI(2)
            DLIJETADXI(2) = -2.0*ETA*DPDXI(2)
            DLIJETADXI(3) =  0.5*(2.0*ETA+1.0)*DPDXI(2)
                 
            do I = 1, 3
               do J = 1, 3
                  ILOC = I+(J-1)*3
                  SVN(ILOC,GI) = LIJXI(I)*LIJETA(J)
                  SVNLX(ILOC,GI) = LIJXI(I)*DLIJETADXI(J) + LIJETA(J)*DLIJXIDXI(I)
               END DO
            END DO
            
         END DO Loop_GJ
         
      END DO Loop_IFACE

!!$ set the weights for the surface integration
      SVWEIGH = 1.
      

      return
    END SUBROUTINE FVQQUAD

!!$ ===============
!!$ ===============



  end module FV_Quadratures
