
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


!!!!==============================================!!!!
!!!!        SHAPE FUNCTIONS SUBRTS FOR 3D         !!!!
!!!!         (Hexahedra and   Tetrahedra)         !!!!
!!!!==============================================!!!!


  module HexTet_ShapeFunctions

    use fldebug
    use multi_data_types
    use QuadTri_ShapeFunctions



  contains


!!$ ===============
!!$ ===============

    subroutine re3dn8( lowqua, ngi, ngi_l, nloc, mloc, &
         m, weight, n, nlx, nly, nlz, &
         sngi, snloc, sweigh, sn, snlx, snly, &
         l1, l2, l3 )
!!$    This subrt. computes the shape functions M and N and their
!!$        derivatives at the Gauss points for 3D.
!!$        If LOWQUA, then use one point quadrature else use 8 point quadrature.
!!$        !NB.: LX/YP(I) are the local X/Y coordinates of nodal point I.
      implicit none
      logical, intent( in ) :: lowqua
      integer, intent( in ) :: ngi, ngi_l, nloc, mloc
      real, dimension( :, : ), intent( inout ) :: m
      real, dimension( : ), intent( inout ) :: weight
      real, dimension( :, : ), intent( inout ) :: n, nlx, nly, nlz
      integer, intent( in ) :: sngi, snloc
      real, dimension( : ), intent( inout ) :: sweigh
      real, dimension( :, : ), intent( inout ) :: sn, snlx, snly
      real, dimension( : ), intent( in ) :: l1, l2, l3
      ! Local variables:
      integer, parameter :: nl = 4, nlp = 8, npq = 2
      real, dimension( : ), allocatable :: lx, ly, lz, lxp, lyp, lzp, weit, rdum2
      integer :: p, q, ir, corn, gpoi, ngi1d, gi
      real :: posi, rdum, rlx, rly, rlz

      ! Allocating memory
      allocate( lx( nl ), ly( nl ), lz( nl ), lxp( nlp ), lyp( nlp ), lzp( nlp ), &
           weit( nl ), rdum2( ngi_l ) )

      if( sngi >= 1 ) then ! Surface integrals
         call re2dn4( lowqua, sngi, 0, snloc, 0, &
              m, sweigh, sn, snlx, snly, &
              0, 0, sweigh, sn, snlx, &
              rdum2, rdum2 )
      end if

!!$                  Nodal Point 1
      lxp( 1 ) = -1   ;   lyp( 1 ) = -1   ;   lzp( 1 ) = -1
!!$                  Nodal Point 2
      lxp( 2 ) =  1   ;   lyp( 2 ) = -1   ;   lzp( 2 ) = -1
!!$                  Nodal Point 3
      lxp( 3 ) = -1   ;   lyp( 3 ) =  1   ;   lzp( 3 ) = -1
!!$                  Nodal Point 4
      lxp( 4 ) =  1   ;   lyp( 4 ) =  1   ;   lzp( 4 ) = -1
!!$                  Nodal Point 5
      lxp( 5 ) = -1   ;   lyp( 5 ) = -1   ;   lzp( 5 ) =  1
!!$                  Nodal Point 6
      lxp( 6 ) =  1   ;   lyp( 6 ) = -1   ;   lzp( 6 ) =  1
!!$                  Nodal Point 7
      lxp( 7 ) = -1   ;   lyp( 7 ) =  1   ;   lzp( 7 ) =  1
!!$                  Nodal Point 8
      lxp( 8 ) =  1   ;   lyp( 8 ) =  1   ;   lzp( 8 ) =  1

      Conditional_NGI: Select Case( NGI )

      case( 8 )
         posi = 1. / sqrt( 3. )
         lx( 1 ) = -posi ; ly( 1 ) = -posi ; lz( 1 ) = -posi
         lx( 2 ) =  posi ; ly( 2 ) =  posi ; lz( 2 ) =  posi

         Loop_P1: do p = 1, npq
            Loop_Q1: do q = 1, npq
               Loop_IR1: do ir = 1, npq
                  Loop_Corn1: do corn = 1, nlp
                     gpoi = ( q - 1 ) * npq + p  + ( ir - 1 ) * npq * npq
                     rlx = lx( p ) ; rly = ly( q ) ; rlz = lz( ir )
                     if( ngi_l /= 0 ) then
                        rlx = l1( gpoi ) ; rly = l2( gpoi ) ; rlz = l3( gpoi )
                     endif
                     if( mloc > 0 ) m( 1, gpoi ) = 1.

                     weight( gpoi ) = 1. ! Weight

                     n( corn, gpoi ) = 0.125 * ( 1. + lxp( corn ) * rlx ) * &
                          ( 1. + lyp( corn ) * rly ) * ( 1. + lzp( corn ) * &
                          rlz ) ! N

                     nlx( corn, gpoi ) = 0.125 * lxp( corn ) * ( 1. + &
                          lyp( corn ) * rly ) * ( 1. + lzp( corn ) * rlz ) ! x-derivative

                     nly( corn, gpoi ) = 0.125 * lyp( corn ) * ( 1. + lxp( corn ) * &
                          rlx ) * ( 1. + lzp( corn ) * rlz ) ! y-derivative

                     nlz( corn, gpoi ) = 0.125 * lzp( corn ) * ( 1. + lxp( corn ) * &
                          rlx ) * ( 1. + lyp( corn ) * rly ) ! z-derivative
                  end do Loop_Corn1
               end do Loop_IR1
            end do Loop_Q1
         end do Loop_P1

      case default

         ngi1d = int( ( ngi + 0.1 ) ** ( 1. / 3. ) + 0.1 )
         do gi = 1, ngi1d
            weit( gi ) = rgptwe( gi, ngi1d, .true. )
            lx( gi ) = rgptwe( gi, ngi1d, .false. )
            ly( gi ) = lx( gi )
            lz( gi ) = lx( gi )
         end do
         Loop_P2: do p = 1, ngi1d
            Loop_Q2: do q = 1, ngi1d
               Loop_IR2: do ir = 1, ngi1d
                  Loop_Corn2: do corn = 1, nlp
                     gpoi = ( ( p - 1 ) * ngi1d + q ) + ( ir - 1 ) * ngi1d **2
                     rlx = lx( p ) ; rly = ly( q ) ; rlz = lz( ir )
                     if( ngi_l /= 0 ) then
                        rlx = l1( gpoi ) ; rly = l2( gpoi ) ; rlz = l3( gpoi )
                     endif
                     if( mloc > 0 ) m( 1, gpoi ) = 1.

                     weight( gpoi ) = weit( p ) * weit( q ) * weit( ir ) ! Weight

                     n( corn, gpoi ) = 0.125 * ( 1. + lxp( corn ) * rlx ) * &
                          ( 1. + lyp( corn ) * rly ) * ( 1. + lzp( corn ) * rlz ) ! N

                     nlx( corn, gpoi ) = 0.125 * lxp( corn ) * ( 1. + lyp( corn ) * &
                          rly ) * ( 1. + lzp( corn ) * rlz )  ! x-derivative

                     nly( corn, gpoi ) = 0.125 * lyp( corn ) * ( 1. + lxp( corn ) * &
                          rlx ) * ( 1. + lzp( corn ) * rlz ) ! y-derivative

                     nlz( corn, gpoi ) = 0.125 * lzp( corn ) * ( 1. + lxp( corn ) * &
                          rlx ) * ( 1. + lyp( corn ) * rly ) ! z-derivative
                  end do Loop_Corn2
               end do Loop_IR2
            end do Loop_Q2
         end do Loop_P2

      end Select Conditional_NGI

      ! Deallocating
      deallocate( lx, ly, lz, lxp, lyp, lzp, weit, rdum2 )

      return
    end subroutine re3dn8




  end module HexTet_ShapeFunctions
