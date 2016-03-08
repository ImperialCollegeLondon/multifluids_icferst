
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

module All_ShaTri

    use fldebug


contains


    SUBROUTINE SHATRInew(L1, L2, L3, L4, WEIGHT, &
         NLOC,NGI,  N,NLX_ALL)
      ! Interface to SHATRIold using the new style variables
      IMPLICIT NONE
      INTEGER , intent(in) :: NLOC,NGI
      REAL , dimension(:), intent(in) :: L1, L2, L3, L4
      REAL , dimension(:), intent(inout) :: WEIGHT
      REAL , dimension(:, : ), intent(inout) ::N
      real, dimension (:,:,:), intent(inout) :: NLX_ALL

      call SHATRIold(L1, L2, L3, L4, WEIGHT, size(NLX_ALL,1)==3, &
           NLOC,NGI,  &
           N,NLX_ALL(1,:,:),NLX_ALL(2,:,:),NLX_ALL(3,:,:))

    end subroutine SHATRInew


!!$ ===============
!!$ ===============


    SUBROUTINE SHATRIold(L1, L2, L3, L4, WEIGHT, D3, &
         NLOC,NGI,  &
         N,NLX,NLY,NLZ) 
      ! Work out the shape functions and there derivatives...
      IMPLICIT NONE
      INTEGER , intent(in) :: NLOC,NGI
      LOGICAL , intent(in) :: D3
      REAL , dimension(:), intent(in) :: L1, L2, L3, L4
      REAL , dimension(:), intent(inout) :: WEIGHT
      REAL , dimension(:, : ), intent(inout) ::N,NLX,NLY,NLZ
      ! Local variables...
      INTEGER ::  GI
      !
      IF(.NOT.D3) THEN
         ! Assume a triangle...
         !
         IF(NLOC.EQ.1) THEN
            Loop_Gi_Nloc1: DO GI=1,NGI
               N(1,GI)=1.0
               NLX(1,GI)=0.0
               NLY(1,GI)=0.0
            end DO Loop_Gi_Nloc1
         ELSE IF((NLOC.EQ.3).OR.(NLOC.EQ.4)) THEN
            Loop_Gi_Nloc3_4: DO GI=1,NGI
               N(1,GI)=L1(GI)
               N(2,GI)=L2(GI)
               N(3,GI)=L3(GI)
               !
               NLX(1,GI)=1.0
               NLX(2,GI)=0.0
               NLX(3,GI)=-1.0
               !
               NLY(1,GI)=0.0
               NLY(2,GI)=1.0
               NLY(3,GI)=-1.0
               IF(NLOC.EQ.4) THEN
                  ! Bubble function...
                  N(4,GI)  =L1(GI)*L2(GI)*L3(GI)
                  NLX(4,GI)=L2(GI)*(1.-L2(GI))-2.*L1(GI)*L2(GI)
                  NLY(4,GI)=L1(GI)*(1.-L1(GI))-2.*L1(GI)*L2(GI)
               ENDIF
            end DO Loop_Gi_Nloc3_4
         ELSE IF((NLOC.EQ.6).OR.(NLOC.EQ.7)) THEN
            Loop_Gi_Nloc_6_7: DO GI=1,NGI
               N(1,GI)=(2.*L1(GI)-1.)*L1(GI)
               N(2,GI)=(2.*L2(GI)-1.)*L2(GI)
               N(3,GI)=(2.*L3(GI)-1.)*L3(GI)
               !
               N(4,GI)=4.*L1(GI)*L2(GI)
               N(5,GI)=4.*L2(GI)*L3(GI)
               N(6,GI)=4.*L1(GI)*L3(GI)

               !
               ! nb L1+L2+L3+L4=1
               ! x-derivative...
               NLX(1,GI)=4.*L1(GI)-1.
               NLX(2,GI)=0.
               NLX(3,GI)=-4.*(1.-L2(GI))+4.*L1(GI) + 1.
               !
               NLX(4,GI)=4.*L2(GI)
               NLX(5,GI)=-4.*L2(GI)
               NLX(6,GI)=4.*(1.-L2(GI))-8.*L1(GI)
               !
               ! y-derivative...
               NLY(1,GI)=0.
               NLY(2,GI)=4.*L2(GI)-1.0
               NLY(3,GI)=-4.*(1.-L1(GI))+4.*L2(GI) + 1.
               !
               NLY(4,GI)=4.*L1(GI)
               NLY(5,GI)=4.*(1.-L1(GI))-8.*L2(GI)
               NLY(6,GI)=-4.*L1(GI)
               IF(NLOC.EQ.7) THEN
                  ! Bubble function...
                  N(7,GI)  =L1(GI)*L2(GI)*L3(GI)
                  NLX(7,GI)=L2(GI)*(1.-L2(GI))-2.*L1(GI)*L2(GI)
                  NLY(7,GI)=L1(GI)*(1.-L1(GI))-2.*L1(GI)*L2(GI)
               ENDIF
            END DO Loop_Gi_Nloc_6_7
            ! ENDOF IF(NLOC.EQ.6) THEN...
         ELSE IF(NLOC==10) THEN ! Cubic triangle...
            ! get the shape functions for a cubic triangle...
            call shape_triangle_cubic( l1, l2, l3, l4, weight, d3, &
                 nloc, ngi, &
                 n, nlx, nly, nlz )

         ELSE ! has not found the element shape functions
            stop 811
         ENDIF
         !
         ! ENDOF IF(.NOT.D3) THEN
      ENDIF
      !
      !
      IF(D3) THEN
         ! Assume a tet...
         ! This is for 5 point quadrature.
         IF((NLOC.EQ.10).OR.(NLOC.EQ.11)) THEN
            Loop_Gi_Nloc_10_11: DO GI=1,NGI
               !ewrite(3,*)'gi,L1(GI),L2(GI),L3(GI),L4(GI):',gi,L1(GI),L2(GI),L3(GI),L4(GI)
               N(1,GI)=(2.*L1(GI)-1.)*L1(GI)
               N(3,GI)=(2.*L2(GI)-1.)*L2(GI)
               N(5,GI)=(2.*L3(GI)-1.)*L3(GI)
               N(10,GI)=(2.*L4(GI)-1.)*L4(GI)

               !if(L1(GI).gt.-1.93) ewrite(3,*)'gi,L1(GI), L2(GI), L3(GI), L4(GI),N(1,GI):', &
               !                            gi,L1(GI), L2(GI), L3(GI), L4(GI),N(1,GI)
               !
               !
               N(2,GI)=4.*L1(GI)*L2(GI)
               N(6,GI)=4.*L1(GI)*L3(GI)
               N(7,GI)=4.*L1(GI)*L4(GI)
               !
               N(4,GI) =4.*L2(GI)*L3(GI)
               N(9,GI) =4.*L3(GI)*L4(GI)
               N(8,GI)=4.*L2(GI)*L4(GI)
               ! nb L1+L2+L3+L4=1
               ! x-derivative...
               NLX(1,GI)=4.*L1(GI)-1.
               NLX(3,GI)=0.
               NLX(5,GI)=0.
               NLX(10,GI)=-4.*(1.-L2(GI)-L3(GI))+4.*L1(GI) + 1.
               !if(L1(GI).gt.-1.93) ewrite(3,*)'Nlx(1,GI):', &
               !     Nlx(1,GI)
               !
               NLX(2,GI)=4.*L2(GI)
               NLX(6,GI)=4.*L3(GI)
               NLX(7,GI)=4.*(L4(GI)-L1(GI))
               !
               NLX(4,GI) =0.
               NLX(9,GI) =-4.*L3(GI)
               NLX(8,GI)=-4.*L2(GI)
               !
               ! y-derivative...
               NLY(1,GI)=0.
               NLY(3,GI)=4.*L2(GI)-1.0
               NLY(5,GI)=0.
               NLY(10,GI)=-4.*(1.-L1(GI)-L3(GI))+4.*L2(GI) + 1.
               !
               NLY(2,GI)=4.*L1(GI)
               NLY(6,GI)=0.
               NLY(7,GI)=-4.*L1(GI)
               !
               NLY(4,GI) =4.*L3(GI)
               NLY(9,GI) =-4.*L3(GI)
               NLY(8,GI)=4.*(1-L1(GI)-L3(GI))-8.*L2(GI)
               !
               ! z-derivative...
               NLZ(1,GI)=0.
               NLZ(3,GI)=0.
               NLZ(5,GI)=4.*L3(GI)-1.
               NLZ(10,GI)=-4.*(1.-L1(GI)-L2(GI))+4.*L3(GI) + 1.
               !
               NLZ(2,GI)=0.
               NLZ(6,GI)=4.*L1(GI)
               NLZ(7,GI)=-4.*L1(GI)
               !
               NLZ(4,GI) =4.*L2(GI)
               NLZ(9,GI) =4.*(1.-L1(GI)-L2(GI))-8.*L3(GI)
               NLZ(8,GI)=-4.*L2(GI)
               IF(NLOC.EQ.11) THEN
                  ! Bubble function...
                  N(11,GI)  =L1(GI)*L2(GI)*L3(GI)*L4(GI)
                  NLX(11,GI)=L2(GI)*L3(GI)*(1.-L2(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
                  NLY(11,GI)=L1(GI)*L3(GI)*(1.-L1(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
                  NLZ(11,GI)=L1(GI)*L2(GI)*(1.-L1(GI)-L2(GI))-2.*L1(GI)*L2(GI)*L3(GI)
               ENDIF
               !
            end DO Loop_Gi_Nloc_10_11
            ! ENDOF IF(NLOC.EQ.10) THEN...
         ENDIF
         !
         IF((NLOC.EQ.4).OR.(NLOC.EQ.5)) THEN
            Loop_Gi_Nloc_4_5: DO GI=1,NGI
               N(1,GI)=L1(GI)
               N(2,GI)=L2(GI)
               N(3,GI)=L3(GI)
               N(4,GI)=L4(GI)
               !
               NLX(1,GI)=1.0
               NLX(2,GI)=0
               NLX(3,GI)=0
               NLX(4,GI)=-1.0
               !
               NLY(1,GI)=0.0
               NLY(2,GI)=1.0
               NLY(3,GI)=0.0
               NLY(4,GI)=-1.0
               !
               NLZ(1,GI)=0.0
               NLZ(2,GI)=0.0
               NLZ(3,GI)=1.0
               NLZ(4,GI)=-1.0
               IF(NLOC.EQ.5) THEN
                  ! Bubble function...
                  N(5,GI)  =L1(GI)*L2(GI)*L3(GI)*L4(GI)
                  NLX(5,GI)=L2(GI)*L3(GI)*(1.-L2(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
                  NLY(5,GI)=L1(GI)*L3(GI)*(1.-L1(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
                  NLZ(5,GI)=L1(GI)*L2(GI)*(1.-L1(GI)-L2(GI))-2.*L1(GI)*L2(GI)*L3(GI)
               ENDIF
            end DO Loop_Gi_Nloc_4_5
         ENDIF
         !
         IF(NLOC.EQ.1) THEN
            Loop_Gi_Nloc_1: DO GI=1,NGI
               N(1,GI)=1.0
               NLX(1,GI)=0.0
               NLY(1,GI)=0.0
               NLZ(1,GI)=0.0
            end DO Loop_Gi_Nloc_1
         ENDIF
         !
         ! ENDOF IF(D3) THEN...
      ENDIF
      !
      RETURN
    END SUBROUTINE SHATRIold

!!$ ===============
!!$ ===============

  subroutine shatri_hex( l1, l2, l3, l4, weight, d3, &
       nloc, ngi, &
       n, nlx, nly, nlz, &
       tri_tet )
    implicit none
    integer, intent( in ) :: nloc, ngi
    logical, intent( in ) :: tri_tet
    real, dimension( : ), intent( in ) :: l1, l2, l3, l4
    real, dimension( : ), intent( inout ) :: weight
    logical, intent( in ) :: d3
    real, dimension( :, : ), intent( inout ) :: n, nlx, nly, nlz
    ! Local variables
    logical :: lowqua
    integer :: nwicel, mloc, snloc, sngi
    real, dimension( : ), allocatable :: rdum, sweigh
    real, dimension( :, : ), allocatable ::  m, sn, snlx, snly

    !ewrite(3,*)'In shatri_hex'

    Conditional_Dimension: if( tri_tet ) then ! traingles and tets

       call shatri( l1, l2, l3, l4, weight, d3, &
            nloc, ngi, &
            n, nlx, nly, nlz )

    else
       ! Get the shape functions on lines (in 2D) and quadrilateral surfaces in 3D
       if( .not. d3 ) then ! 2D:
          if( nloc == 4 ) nwicel = 1
          if( nloc == 9 ) nwicel = 3
       else
          if( nloc == 8 ) nwicel = 1
          if( nloc == 27 ) nwicel = 3
       endif
       lowqua = .false.
       mloc = 0
       snloc = 0
       sngi = 0

       Select Case( nwicel )

       case default; FLExit( "Wrong option for NWICEL " )

       case( 1 )
          allocate( rdum( 0 ) )
          allocate( m( mloc, ngi ) )
          allocate( sweigh( sngi ) )
          allocate( sn( snloc, sngi ) )
          allocate( snlx( snloc, sngi ) )
          allocate( snly( snloc, sngi ) )
          if( .not. d3 ) then
             call re2dn4( lowqua, ngi, 0, nloc, mloc, &
                  m, weight, n, nlx, nly, &
                  sngi, snloc, sweigh, sn, snlx, &
                  rdum, rdum )
          else
             call re3dn8( lowqua, ngi, 0, nloc, mloc, &
                  m, weight, n, nlx, nly, nlz, &
                  sngi, snloc, sweigh, sn, snlx, snly, &
                  rdum, rdum, rdum )
          end if
          deallocate( rdum )
          deallocate( m )
          deallocate( sweigh )
          deallocate( sn )
          deallocate( snlx )
          deallocate( snly )

       case( 3 )
          if( .not. d3 ) then
             call re2dn9( lowqua, ngi, 0, nloc, mloc, &
                  m, weight, n, nlx, nly, &
                  rdum, rdum )
          else ! Lagrange 27 nodes 3D element - bilinear pressure
             call re3d27( lowqua, ngi, 0, nloc, mloc, &
                  m, weight, n, nlx, nly, nlz, &
                  rdum, rdum, rdum )
          end if

       end Select

    endif Conditional_Dimension

    return
  end subroutine shatri_hex

!!$ ===============
!!$ ===============


  subroutine shatri( l1, l2, l3, l4, weight, d3, &
       nloc, ngi, &
       n, nlx, nly, nlz, shatriold )
    implicit none
    integer, intent( in ) :: nloc, ngi
    real, dimension( : ), intent( in ) :: l1, l2, l3, l4, weight
    logical, intent( in ) :: d3
    real, dimension( :, : ), intent( inout ) :: n, nlx, nly, nlz
    logical, intent( in ) :: shatriold
    ! Local variables
    logical :: base_order
    integer :: gi
    real :: a,b

    base_order = ( .not. shatriold )

    Conditional_Dimensionality: if( .not. d3 ) then ! Triangle

       Conditional_NLOC: Select Case( nloc )
       case( 10 ) !!$ cubic triangle
          call shape_triangle_cubic( l1, l2, l3, l4, weight, d3, &
               nloc, ngi, &
               n, nlx, nly, nlz ) !!$ Shape functions for a cubic triangle

          if(base_order) then !!$ order so that the 1st nodes are on the base...
             call base_order_tri(n,nloc,ngi)
             call base_order_tri(nlx,nloc,ngi)
             call base_order_tri(nly,nloc,ngi)
          endif

       case( 6, 7 )
          do gi = 1, ngi
             n( 1, gi ) = ( 2. * l1( gi ) - 1. ) * l1( gi )
             n( 2, gi ) = ( 2. * l2( gi ) - 1. ) * l2( gi )
             n( 3, gi ) = ( 2. * l3( gi ) - 1. ) * l3( gi )
             n( 4, gi ) = 4. * l1( gi ) * l2( gi )
             n( 5, gi ) = 4. * l2( gi ) * l3( gi )
             n( 6, gi ) = 4. * l1( gi ) * l3( gi )
             ! x-derivative (nb. l1 + l2 + l3  = 1 )
             nlx( 1, gi ) =  4. * l1( gi ) - 1.
             nlx( 2, gi ) =  0.
             nlx( 3, gi ) = -4. * ( 1. - l2( gi ) ) + 4. * l1( gi ) + 1.
             nlx( 4, gi ) =  4. * l2( gi )
             nlx( 5, gi ) = -4. * l2( gi )
             nlx( 6, gi ) =  4. * ( 1. - l2( gi ) ) -8. * l1( gi )
             ! y derivative
             nly( 1, gi ) =  0.
             nly( 2, gi ) =  4. * l2( gi ) - 1.
             nly( 3, gi ) = -4. * ( 1. - l1( gi ) ) + 4. * l2( gi ) + 1.
             nly( 4, gi ) =  4. * l1( gi )
             nly( 5, gi ) =  4. * ( 1. - l1( gi ) ) - 8. * l2( gi )
             nly( 6, gi ) = -4 * ( l1( gi ) )
             if( nloc == 7 ) then  ! Bubble function
                n( 7, gi )   = l1( gi ) * l2( gi ) * l3( gi )
                nlx( 7, gi ) = l2( gi ) * ( 1. - l2( gi ) ) - 2. * l1( gi ) * l2( gi )
                nly( 7, gi ) = l1( gi ) * ( 1. - l1( gi ) ) - 2. * l1( gi ) * l2( gi )
             end if
          end do

          if(base_order) then
             ! order so that the 1st nodes are on the base...
             call base_order_tri(n,nloc,ngi)
             call base_order_tri(nlx,nloc,ngi)
             call base_order_tri(nly,nloc,ngi)
          endif

       case( 3, 4 )
          do gi = 1, ngi
             n( 1, gi ) = l1( gi ) ; n( 2, gi ) = l2( gi ) ; n( 3, gi ) = l3( gi )
             ! x-derivative
             nlx( 1, gi ) = 1. ; nlx( 2, gi ) = 0. ; nlx( 3, gi ) = -1.
             ! y-derivative
             nly( 1, gi ) = 0. ; nly( 2, gi ) = 1. ; nly( 3, gi ) = -1.
             if( nloc == 4 ) then ! Bubble function
                n( 4, gi )   = l1( gi ) * l2( gi ) * l3( gi )
                nlx( 4, gi ) = l2( gi ) * ( 1. - l2( gi )) - 2. * l1( gi ) * l2( gi )
                nly( 4, gi ) = l1( gi ) * ( 1. - l1( gi )) - 2. * l1( gi ) * l2( gi )
             end if
          end do

       case( 1 )
          n( 1, : ) = 0. ; nlx( 1, : ) = 0. ; nly( 1, : ) = 0.

       case default; FLExit( " Wrong number of NLOC " )
       end Select Conditional_NLOC

    else ! Assume it is a tetrahedron

       Conditional_NLOC2: Select Case( nloc )
       case( 10, 11 ) ! 5 points quadrature
          do gi = 1, ngi
             n( 1, gi )  = ( 2. * l1( gi ) - 1. ) * l1( gi )
             n( 2, gi )  = 4. * l1( gi ) * l2( gi )
             n( 3, gi )  = ( 2. * l2( gi ) - 1. ) * l2( gi )
             n( 4, gi )  = 4. * l2( gi ) * l3( gi )
             n( 5, gi )  = ( 2. * l3( gi ) - 1. ) * l3( gi )
             n( 6, gi )  = 4. * l1( gi ) * l3( gi )
             n( 7, gi )  = 4. * l1( gi ) * l4( gi )
             n( 8, gi )  = 4. * l2( gi ) * l4( gi )
             n( 9, gi )  = 4. * l3( gi ) * l4( gi )
             n( 10, gi ) = ( 2. * l4( gi ) - 1. ) * l4( gi )
!!$    x-derivative (nb. l1 + l2 + l3  = 1 )
             nlx( 1, gi )  =  4. * l1( gi ) - 1.
             nlx( 2, gi )  =  4. * l2( gi )
             nlx( 3, gi )  =  0.
             nlx( 4, gi )  =  0.
             nlx( 5, gi )  =  0.
             nlx( 6, gi )  =  4. * l3( gi )
             nlx( 7, gi )  =  4. * ( l4( gi ) - l1( gi )
             nlx( 8, gi )  = -4. * l2( gi )
             nlx( 9, gi )  = -4. * l3( gi )
             nlx( 10, gi ) = -4. * ( 1. - l2( gi ) -l3( gi ) ) + 4. * l1( gi ) + 1.
!!$    y derivative
             nly( 1, gi )  =  0.
             nly( 2, gi )  =  4. * l1( gi )
             nly( 3, gi )  =  4. * l2( gi ) - 1.
             nly( 4, gi )  =  4. * l3( gi )
             nly( 5, gi )  =  0.
             nly( 6, gi )  =  0.
             nly( 7, gi )  = -4. * l1( gi )
             nly( 8, gi ) =  4. * ( 1. - l1( gi ) - l3( gi ) ) + 8. * l2( gi )
             nly( 9, gi ) = -4. * l3( gi )
             nly( 10, gi ) = -4. * ( 1. - l1( gi ) - l3( gi ) ) + 4. * l2( gi ) + 1.
!!$    y derivative
             nlz( 1, gi )  =  0.
             nlz( 2, gi )  =  0.
             nlz( 3, gi )  =  0.
             nlz( 4, gi )  =  4. * l2( gi )
             nlz( 5, gi )  =  4. * l3( gi ) - 1.
             nlz( 6, gi )  =  4. * l1( gi )
             nlz( 7, gi )  = -4. * l1( gi )
             nlz( 8, gi )  = -4. * l2( gi )
             nlz( 9, gi )  = -4. * ( 1. - l1( gi ) - l2( gi ) ) - 8. * l3( gi )
             nlz( 10, gi ) = -4. * ( 1. - l1( gi ) - l2( gi ) ) + 4. * l3( gi ) + 1.
          end do

          if(base_order) then
             ! order so that the 1st nodes are on the base...
             call base_order_tet(n,nloc,ngi)
             call base_order_tet(nlx,nloc,ngi)
             call base_order_tet(nly,nloc,ngi)
             call base_order_tet(nlz,nloc,ngi)
          endif

       case( 4, 5 )
          do gi = 1, ngi
             n( 1, gi ) = l1( gi )
             n( 2, gi ) = l2( gi )
             n( 3, gi ) = l3( gi )
             n( 4, gi ) = l4( gi )
             ! x-derivative
             nlx( 1, gi ) = 1.
             nlx( 2, gi ) = 0.
             nlx( 3, gi ) = 0.
             nlx( 4, gi ) = -1.
             ! y-derivative
             nly( 1, gi ) = 0.
             nly( 2, gi ) = 1.
             nly( 3, gi ) = 0.
             nly( 4, gi ) = -1.
             ! z-derivative
             nlz( 1, gi ) = 0.
             nlz( 2, gi ) = 0.
             nlz( 3, gi ) = 1.
             nlz( 4, gi ) = -1.
             if( nloc == 5 ) then ! Bubble function
                n( 5, gi ) = l1( gi ) * l2( gi ) * l3( gi ) * l4( gi )
                nlx( 5, gi ) = l2( gi ) * l3( gi ) * ( 1. - l2( gi ) - l3( gi ))  &
                     -2. * ( l1( gi ) * l2( gi ) * l3( gi ) )
                nly( 5, gi ) = l1( gi ) * l3( gi ) * ( 1. - l1( gi ) - l3( gi ))  &
                     -2. * ( l1( gi ) * l2( gi ) * l3( gi ) )
                nlz( 5, gi ) = l1( gi ) * l2( gi ) * ( 1. - l2( gi ) - l2( gi ))  &
                     -2. * ( l1( gi ) * l2( gi ) * l3( gi ) )
             endif
          end do

       case( 1 )
             n( 1, : ) = 1. ; nlx( 1, : ) = 0. ; nly( 1, : ) = 0. ; nlz( 1, : ) = 0.

       case default; FLExit( "Wrong number for NLOC " )

       end Select Conditional_NLOC2

    end if Conditional_Dimensionality

    return
  end subroutine shatri

!!$ ===============
!!$ ===============



end module All_ShaTri
