  
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
  
  
!!!!=======================================================!!!!
!!!!   SHAPE FUNCTIONS: SUBRTS FOR LINEAR AND QUADRATIC    !!!!
!!!!                  QUADS AND HEXAHEDRA                  !!!!
!!!!=======================================================!!!!
  
  
  module QuadsHex_ShapeFunctions

    use fldebug

  contains

!!$ ===============
!!$ ===============

    subroutine re2dn4( lowqua, ngi, ngi_l, nloc, mloc, &
         m, weight, n, nlx, nly, &
         sngi, snloc, sweigh, sn, snlx, &
         l1, l2 )
!!$      NB:  Linear quads
!!$      This subrt computes shape functions M and N and their derivatives
!!$          at the Gauss points.
!!$      NB: We may need to define surface elements for p and (u,v,w)
      implicit none
      logical, intent( in ) :: lowqua
      integer, intent( in ) :: ngi, ngi_l, nloc, mloc
      real, dimension( :, : ), intent( inout ) :: m
      real, dimension( : ), intent( inout ) :: weight
      real, dimension( :, : ), intent( inout ) :: n, nlx, nly
      integer, intent( in ) :: sngi, snloc
      real, dimension( : ), intent( inout ) :: sweigh
      real, dimension( :, : ), intent( inout ) :: sn, snlx
      real, dimension(: ), intent( in ) :: l1, l2
      ! Local variables:
      integer, parameter :: nl = 16, nlp = 4, npq = 2
      real, dimension( : ), allocatable :: lx, ly, lxp, lyp, weit
      integer :: q, p, gpoi, corn, ndgi, i, nd_quad_1d
      real :: posi, rlx, rly
      logical :: getndp

      !ewrite(3,*)' In re2dn4 subrt. '
      ! Allocating memory
      allocate( lx( nl ), ly( nl ), lxp( nlp ), lyp( nlp ) )

      ! LX/YP(I) are the local X and Y coordinates of nodal point I
      lxp( 1 ) = -1 ; lyp( 1 ) = -1
      lxp( 2 ) = 1  ; lyp( 2 ) = -1
      lxp( 3 ) = -1 ; lyp( 3 ) = 1
      lxp( 4 ) = 1  ; lyp( 4 ) = 1

      Conditional_NGI: if ( ( ngi == 4 ) .or. ( ngi == 1 ) ) then
         posi = 1. / sqrt( 3. )
         lx( 1 ) = -posi ; ly( 1 ) = -posi
         lx( 2 ) = posi  ; ly( 2 ) = posi
         nd_quad_1d = 2
         if ( ngi == 1 ) then
            lx( 1 ) = 0. ; ly( 1 ) = 0. ; nd_quad_1d = 1
         end if
         Loop_Q: do q = 1, nd_quad_1d
            Loop_P: do p = 1, nd_quad_1d
               Loop_Corn: do corn = 1, nlp
                  gpoi = ( q - 1 ) * nd_quad_1d + p
                  rlx = lx( p ) ; rly = ly( q )
                  if( ngi_l /= 0 ) then
                     rlx = l1( gpoi ) ; rly = l2( gpoi )
                  endif
                  if( mloc == 1 ) m( 1, gpoi ) = 1.
                  weight( gpoi ) = 1.
                  if( ngi == 1 ) weight( gpoi ) = 4. ! It's a convolution of two weights of 2.
                  n( corn, gpoi ) = 0.25 * ( 1. + lxp( corn ) * rlx ) * &
                       ( 1. + lyp( corn ) * rly )
                  nlx( corn, gpoi ) = 0.25 * lxp( corn ) * ( 1. + lyp( corn ) * rly )
                  nly( corn, gpoi ) = 0.25 * lyp( corn ) * ( 1. + lxp( corn ) * rlx )
               end do Loop_Corn
            end do Loop_P
         end do Loop_Q

         nd_quad_1d = 2
         if ( sngi == 1 ) then
            lx( 1 ) = 0. ; ly( 1 ) = 0. ; nd_quad_1d = 1
         end if

         if( ( sngi == 1 ) .and. ( snloc > 1 ) ) then ! Surface shape functions
            Loop_P2: do p = 1, nd_quad_1d
               Loop_Corn2: do corn = 1, 2
                  gpoi = p
                  sn( corn, gpoi ) = 0.5 * ( 1. + lxp( corn ) * lx( p ) )
                  snlx( corn, gpoi ) = 0.5 * lxp( corn )
                  sweigh( gpoi ) = 1.
                  if( sngi == 1 ) sweigh( gpoi ) = 2.
               end do Loop_Corn2
            end do Loop_P2
         end if

      else ! If ngi =/ 4
         ndgi = int( sqrt( ngi + 0.1 ) + 0.1 )
         getndp = .false.
         allocate( weit( ndgi ) )
         call lagrot( weit, lx, ndgi, getndp )
         ly( 1 : ndgi ) = lx( 1 : ndgi )
         Loop_Q3: do q = 1, ndgi
            Loop_P3: do p = 1, ndgi
               Loop_Corn3: do corn = 1, nlp
                  gpoi = ( q - 1 ) * ndgi + p
                  rlx = lx( p ) ; rly = ly( q )
                  if( ngi_l /= 0 ) then
                     rlx = l1( gpoi ) ; rly = l2( gpoi )
                  endif
                  if( mloc == 1 ) m( 1, gpoi ) = 1.
                  weight( gpoi ) = weit( p ) * weit( q )
                  n( corn, gpoi ) = 0.25 * ( 1. + lxp( corn ) * rlx * &
                       ( 1. + lyp( corn ) * rly ) )
                  nlx( corn, gpoi ) = 0.25 * lxp( corn ) * &
                       ( 1. + lyp( corn ) * rly )
                  nly( corn, gpoi ) = 0.25 * lyp( corn ) * &
                       ( 1. + lxp( corn ) * rlx )
               end do Loop_Corn3
            end do Loop_P3
         end do Loop_Q3
         deallocate( weit )

         if( sngi > 0 ) then
            getndp = .false.
            allocate( weit( sngi ) )
            call lagrot( weit, lx, sngi, getndp )
            Loop_P4: do p = 1, sngi
               Loop_Corn4: do corn = 1, npq
                  gpoi = p
                  sn( corn, gpoi ) = 0.5 * ( 1. + lxp( corn ) * lx( p ) )
                  snlx( corn, gpoi ) = 0.5 * lxp( corn )
                  sweigh( gpoi ) = weit( p )
               end do Loop_Corn4
            end do Loop_P4
            deallocate( weit )
         end if

      end if Conditional_NGI

      if( mloc == nloc ) then
         do i = 1, nlp
            do corn = 1, nlp
               m( corn, i ) = n( corn, i )
            end do
         end do
      end if

      ! Deallocating
      deallocate( lx, ly, lxp, lyp )

      return
    end subroutine re2dn4

!!$ ===============
!!$ ===============

    subroutine re3dn8( lowqua, ngi, ngi_l, nloc, mloc, &
         m, weight, n, nlx, nly, nlz, &
         sngi, snloc, sweigh, sn, snlx, snly, &
         l1, l2, l3 )
!!$      NB:  Linear Hexahedra
!!$      This subrt. computes the shape functions M and N and their
!!$          derivatives at the Gauss points for 3D.
!!$          If LOWQUA, then use one point quadrature else use 8 point quadrature.
!!$      NB.: LX/YP(I) are the local X/Y coordinates of nodal point I.
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


!!$ ===============
!!$ ===============

    subroutine re2dn9( lowqua, ngi, ngi_l, nloc, mloc, &
         m, weight, n, nlx, nly, &
         l1, l2 )
!!$      NB:  Quadratic (or bi-lienar) Quads and Hex.
!!$      Quadratic variation (2D) for velocity -- 9 node brick element.
!!$          Linear variation (2D) for pressure -- 4 node brick element.
!!$          NB.: We may need to define surface elements for p and (u,v,w).
!!$      
!!$      This is for the 2-D 27node element, which is number as follows
!!$             |  /Z
!!$            Y| /
!!$             |/
!!$             +-------X
!!$
!!$       For Z = -1 ...
!!$             7   8   9
!!$             4   5   6
!!$             1   2   3
!!$
!!$       For M the shape functions have local node numbers. For Z = -1 ...
!!$              4       3
!!$      
!!$              1       2
!!$
      implicit none
      logical, intent( in ) :: lowqua
      integer, intent( in ) :: ngi, ngi_l, nloc, mloc
      real, dimension( :, : ), intent( inout ) :: m
      real, dimension( : ), intent( inout ) :: weight
      real, dimension( :, : ), intent( inout ) :: n, nlx, nly
      real, dimension( : ), intent( in ) :: l1, l2
      ! Local variables:
      integer, parameter :: nl = 3, nlp = 9, npq = 4
      real, dimension( : ), allocatable :: lx, ly, lz, lxp, lyp, weit, &
           xn, yn, zn, dxn, dyn, dzn
      integer :: nquad, p, q, corn, gpoi, ilx, ily, nj
      real :: posi, rlx, rly

      !ewrite(3,*) 'In re2dn9'

      ! Allocating memory
      allocate( lx( nl ), ly( nl ), lz( nl ), lxp( nlp ), lyp( nlp ), weit( nl ), &
           xn( nl ), yn( nl ), zn( nl ) , dxn( nl ) , dyn( nl ), dzn( nl ) )

      Conditional_LOWQUA:if( ngi == 1 ) then
         lx( 1 ) = 0.0   ;   ly( 1 ) = 0.0   ;   lz( 1 ) = 0.0
         weit( 1 ) = 2.  ;   nquad = 1
      else if( lowqua .or. (ngi == 4) ) then
         posi = 1. / sqrt( 3. )
         lx( 1 ) = -posi   ;   ly( 1 ) = -posi   ;   lz( 1 ) = -posi
         lx( 2 ) =  posi   ;   ly( 2 ) =  posi   ;   lz( 2 ) =  posi
         weit( 1 ) =  1.   ;   weit( 2 ) =  1.   ;   weit( 3 ) =  1.
         nquad = 2
      else
         posi = 0.774596669241483
         lx( 1 ) = -posi   ;   ly( 1 ) = -posi
         lx( 2 ) = 0.      ;   ly( 2 ) = 0.
         lx( 3 ) = posi    ;   ly( 3 ) = posi
         weit( 1 ) = 0.555555555555556
         weit( 2 ) = 0.888888888888889
         weit( 3 ) = 0.555555555555556
         nquad = 3
      end if Conditional_LOWQUA

!!$                  Nodal Point 1
      lxp( 1 ) = -1   ;  lyp( 1 ) = -1
!!$                  Nodal Point 2
      lxp( 2 ) =  1   ;  lyp( 2 ) = -1
!!$                  Nodal Point 3
      lxp( 3 ) =  1   ;  lyp( 3 ) =  1
!!$                  Nodal Point 4
      lxp( 4 ) = -1   ;  lyp( 4 ) =  1

      ! Compute M
      Loop_P1: do p = 1, nquad
         Loop_Q1: do q = 1, nquad
            Loop_Corn1: do corn = 1, npq
               gpoi = ( p - 1 ) * nquad + q
               rlx = lx( p )   ;  rly = ly( q )
               if( ngi_l /= 0 ) then
                  rlx = l1( gpoi )   ;  rly = l2( gpoi )
               endif
               m( corn, gpoi ) = 0.25 * ( 1. + lxp( corn ) * rlx ) * &
                    ( 1 + lyp( corn ) * rly )
            end do Loop_Corn1
         end do Loop_Q1
      end do Loop_P1

      ! Compute N
      Loop_ILX1: do ilx = 1, nl
         Loop_ILY1: do ily = 1, nl
            nj = ilx + ( ily - 1 ) * 3
            lxp( nj ) = real( ilx - 2 )
            lyp( nj ) = real( ily - 2 )
         end do Loop_ILY1
      end do Loop_ILX1

      Loop_P2: do p = 1, nquad
         Loop_Q2: do q = 1, nquad
            gpoi = ( p - 1 ) * nquad + q
            rlx = lx( p )   ;   rly = ly( q )
            if( ngi_l /= 0 ) then
               rlx = l1( gpoi )   ;   rly = l2( gpoi )
            endif
            ! Weight
            weight( gpoi ) = weit( p ) * weit( q )
            ! XN
            xn( 1 ) = 0.5 * rlx * ( rlx - 1. )
            xn( 2 ) = 1. - rlx * rlx
            xn( 3 ) = 0.5 * rlx * ( rlx + 1. )
            ! DXDN
            dxn( 1 ) = 0.5 * ( 2. * rlx - 1. )
            dxn( 2 ) = -2. * rlx
            dxn( 3 ) = 0.5 * ( 2. * rlx + 1. )
            ! YN
            yn( 1 ) = 0.5 * rly * ( rly - 1. )
            yn( 2 ) = 1. - rly * rly
            yn( 3 ) = 0.5 * rly * ( rly + 1. )
            ! DYDN
            dyn( 1 ) = 0.5 * ( 2. * rly - 1. )
            dyn( 2 ) = -2. * rly
            dyn( 3 ) = 0.5 * ( 2. * rly + 1. )
            ! N, NLX and NLY
            Loop_ILX2: do ilx = 1, nl
               Loop_ILY2: do ily = 1, nl
                  nj = ilx + ( ily - 1 ) * 3
                  n( nj, gpoi ) = xn( ilx ) * yn( ily )
                  nlx( nj, gpoi ) = dxn( ilx ) * yn( ily )
                  nly( nj, gpoi ) = xn( ilx ) * dyn( ily )
               end do Loop_ILY2
            end do Loop_ILX2
         end do Loop_Q2
      end do Loop_P2

      ! Deallocating
      deallocate(lx, ly, lz, lxp, lyp, weit, xn, yn, zn, dxn, dyn, dzn )

      return
    end subroutine re2dn9

!!$ ===============
!!$ ===============

    subroutine re3d27( lowqua, ngi, ngi_l, nloc, mloc, &
         m, weight, n, nlx, nly, nlz, &
         l1, l2, l3 )
!!$      NB:  Quadratic (or bi-lienar) Quads and Hex.
!!$      Quadratic variation (3D) for velocity -- 27 node brick element.
!!$          Linear variation (3D) for pressure -- 8 node brick element.
!!$          NB.: We may need to define surface elements for p and (u,v,w).
!!$      
!!$          This is for the 2-D 27node element, which is number as follows
!!$               |  /Z
!!$              Y| /
!!$               |/
!!$               +-------X
!!$      
!!$          For Z = -1 ...
!!$               7   8   9
!!$               4   5   6
!!$               1   2   3
!!$      
!!$          For Z = 0 ...
!!$               16   17   18
!!$               13   14   15
!!$               10   11   12
!!$      
!!$          For Z = 1 ...
!!$               25    26   27
!!$               22    23   24
!!$               19    20   21
!!$      
!!$          For M the shape functions have local node numbers. For Z = -1 ...
!!$                4       3
!!$      
!!$                1       2
!!$      
!!$          For Z = 1 ...
!!$                8       7
!!$      
!!$                5       6
!!$      
      implicit none
      logical, intent( in ) :: lowqua
      integer, intent( in ) :: ngi, ngi_l, nloc, mloc
      real, dimension( :, : ), intent( inout ) :: m
      real, dimension( : ), intent( inout ) :: weight
      real, dimension( :, : ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( : ), intent( in ) :: l1, l2, l3
      ! Local variables:
      integer, parameter :: nl = 3, nlp = 27, npq = 4
      real, dimension( : ), allocatable :: lx, ly, lz, lxp, lyp, lzp, weit, &
           xn, yn, zn, dxn, dyn, dzn
      integer :: nquad, p, q, ir, corn, gpoi, ilx, ily, ilz, nj
      real :: posi

      !ewrite(3,*)'In re3d27'

      ! Allocating memory
      allocate( lx( nl ), ly( nl ), lz( nl ), lxp( nlp ), lyp( nlp ), lzp( nlp ), &
           weit( nl ), xn( nl ), yn( nl ), zn( nl ), dxn( nl ), dyn( nl ), dzn( nl ) )

      Conditional_LOWQUA: if((ngi == 1) ) then
         if(ngi_l /= 0) then
            lx( 1 ) = l1(1) ;   ly( 1 ) = l2(1) ;   lz( 1 ) = l3(1)
         else
            lx( 1 ) = 0.0   ;   ly( 1 ) = 0.0   ;   lz( 1 ) = 0.0
         endif
         weit( 1 ) = 2.0    ;   nquad = 1  !!$ This should be 2, but the 1pt quadrature gets the 
!!$      volume wrong so adjust it.


      else if( lowqua .or. (ngi == 8) ) then
         posi = 1. / sqrt( 3. )
         if( ngi_l /= 0 ) then
            lx( 1 ) = l1(1)  ;  ly( 1 ) = l2(1)  ;  lz( 1 ) = l3(1)
            lx( 2 ) = l1(2)  ;  ly( 2 ) = l2(2)  ;  lz( 2 ) = l3(2)
         else
            lx( 1 ) = -posi  ;  ly( 1 ) = -posi  ;  lz( 1 ) = -posi
            lx( 2 ) =  posi  ;  ly( 2 ) =  posi  ;  lz( 2 ) = posi
         end if
         weit( 1 ) = 1.  ;  weit( 2 ) = 1.  ;  weit( 3 ) = 1.
         nquad = 2
      else
         posi = 0.774596669241483
         if( ngi_l /= 0 ) then
            lx( 1 ) = l1(1)  ;  ly( 1 ) = l2(1)  ;  lz( 1 ) = l3(1)
            lx( 2 ) = l1(2)  ;  ly( 2 ) = l2(2)  ;  lz( 2 ) = l3(2)
            lx( 3 ) = l1(3)  ;  ly( 3 ) = l2(3)  ;  lz( 3 ) = l3(3)
         else
            lx( 1 ) = -posi  ;  ly( 1 ) = -posi  ;  lz( 1 ) = -posi
            lx( 2 ) = 0.     ;  ly( 2 ) = 0.     ;  lz( 2 ) = 0.
            lx( 3 ) =  posi  ;  ly( 3 ) =  posi  ;  lz( 3 ) =  posi
         end if
         weit( 1 ) = 0.555555555555556
         weit( 2 ) = 0.888888888888889
         weit( 3 ) = 0.555555555555556
         nquad = 3
      end if Conditional_LOWQUA

!!$                  Nodal Point 1
      lxp( 1 ) = -1   ;   lyp( 1 ) = -1   ;   lzp( 1 ) = -1
!!$                  Nodal Point 2
      lxp( 2 ) =  1   ;   lyp( 2 ) = -1   ;   lzp( 2 ) = -1
!!$                  Nodal Point 3
      lxp( 3 ) =  1   ;   lyp( 3 ) =  1   ;   lzp( 3 ) = -1
!!$                  Nodal Point 4
      lxp( 4 ) = -1   ;   lyp( 4 ) =  1   ;   lzp( 4 ) = -1
!!$                  Nodal Point 5
      lxp( 5 ) = -1   ;   lyp( 5 ) = -1   ;   lzp( 5 ) =  1
!!$                  Nodal Point 6
      lxp( 6 ) =  1   ;   lyp( 6 ) = -1   ;   lzp( 6 ) =  1
!!$                  Nodal Point 7
      lxp( 7 ) =  1   ;   lyp( 7 ) =  1   ;   lzp( 7 ) =  1
!!$                  Nodal Point 8
      lxp( 8 ) = -1   ;   lyp( 8 ) =  1   ;   lzp( 8 ) =  1

!!$ Compute M
      Loop_P1: do p = 1, nquad
         Loop_Q1: do q = 1, nquad
            Loop_IR1: do ir = 1, nquad
               Loop_Corn1: do corn = 1, npq
                  gpoi = ( p - 1 ) * nquad + q + ( ir - 1 ) * nquad * nquad
                  m( corn, gpoi ) = 0.125 * ( 1. + lxp( corn ) * lx( p ) ) * &
                       ( 1 + lyp( corn ) * ly( q ) * ( 1. + lzp( corn ) * lz( ir ) ) )
               end do Loop_Corn1
            end do Loop_IR1
         end do Loop_Q1
      end do Loop_P1


!!$ Compute N
      Loop_ILX1: do ilx = 1, nl
         Loop_ILY1: do ily = 1, nl
            Loop_ILZ1: do ilz = 1, nl
               nj = ilx + ( ily - 1 ) * 3 + ( ilz - 1 ) * 9
               lxp( nj ) = real( ilx - 2 )
               lyp( nj ) = real( ily - 2 )
               lzp( nj ) = real( ilz - 2 )
            end do Loop_ILZ1
         end do Loop_ILY1
      end do Loop_ILX1

      Loop_P2: do ir = 1, nquad
         Loop_Q2: do q = 1, nquad
            Loop_IR2: do p = 1, nquad
               gpoi = ( p - 1 ) * nquad + q + ( ir - 1 ) * nquad * nquad
!!$                      Weight
               weight( gpoi ) = weit( p ) * weit( q ) * weit( ir )
!!$                        XN
               xn( 1 ) = 0.5 * lx( p ) * ( lx( p ) - 1. )
               xn( 2 ) = 1. - lx( p ) * lx( p )
               xn( 3 ) = 0.5 * lx( p ) * ( lx( p ) + 1. )
!!$                       DXDN
               dxn( 1 ) = 0.5 * ( 2. * lx( p ) - 1. )
               dxn( 2 ) = -2. * lx( p )
               dxn( 3 ) = 0.5 * ( 2. * lx( p ) + 1. )
!!$                        YN
               yn( 1 ) = 0.5 * ly( q ) * ( ly( q ) - 1. )
               yn( 2 ) = 1. - ly( q ) * ly( q )
               yn( 3 ) = 0.5 * ly( q ) * ( ly( q ) + 1. )
!!$                        DYDN
               dyn( 1 ) = 0.5 * ( 2. * ly( q ) - 1. )
               dyn( 2 ) = -2. * ly( q )
               dyn( 3 ) = 0.5 * ( 2. * ly( q ) + 1. )
!!$                         ZN
               zn( 1 ) = 0.5 * lz( ir ) * ( lz( ir ) - 1. )
               zn( 2 ) = 1. - lz( ir ) * lz( ir )
               zn( 3 ) = 0.5 * lz( ir ) * ( lz( ir ) + 1. )
!!$                         DZDN
               dzn( 1 ) = 0.5 * ( 2. * lz( ir ) - 1. )
               dzn( 2 ) = -2. * lz( ir )
               dzn( 3 ) = 0.5 * ( 2. * lz( ir ) + 1. )
!!$                    N, NLX and NLY
               Loop_ILX2: do ilz = 1, nl
                  Loop_ILY2: do ily = 1, nl
                     Loop_ILZ2: do ilx = 1, nl
                        nj = ilx + ( ily - 1 ) * 3 + ( ilz - 1 ) * 9
                        n( nj, gpoi ) = xn( ilx ) * yn( ily ) * zn( ilz )
                        nlx( nj, gpoi ) = dxn( ilx ) * yn( ily ) * zn( ilz )
                        nly( nj, gpoi ) = xn( ilx ) * dyn( ily ) * zn( ilz )
                        nlz( nj, gpoi ) = xn( ilx ) * yn( ily ) * dzn( ilz )
                     end do Loop_ILZ2
                  end do Loop_ILY2
               end do Loop_ILX2

            end do Loop_IR2

         end do Loop_Q2

      end do Loop_P2


      ! Deallocating
      deallocate( lx, ly, lz, lxp, lyp, lzp, weit, xn, yn, zn, dxn, dyn, dzn )

      return
    end subroutine re3d27


!!$ ===============
!!$ ===============


  end module QuadsHex_ShapeFunctions
