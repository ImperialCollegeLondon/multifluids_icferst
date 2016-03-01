
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

  module Quadratures
    use fldebug
    use multi_data_types
    use FV_Quadratures

  contains

!!$ ===============
!!$ ===============

    subroutine quad_1d_shape( cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, n, nlx, un, unlx )
      ! For quadratic elements. Shape functions associated with volume integration
      ! using both CV basis functions CVN as well as FEM basis functions N (and
      ! its derivatives NLX, NLY, NLZ)
      implicit none
      integer, intent( in ) :: cv_ngi, cv_nloc, u_nloc
      real, dimension( :, : ), intent( inout ) :: cvn
      real, dimension( : ), intent( inout ) :: cvweigh
      real, dimension( :, : ), intent( inout ) :: n, nlx
      real, dimension( :, : ), intent( inout ) :: un, unlx
      ! Local Variables
      integer, parameter :: three = 3
      real, dimension( : ), allocatable :: lx, wei, rdummy, xi_min, xi_max, &
           cv_nodpos, u_nodpos
      integer :: gpoi, icv, p, nquad, iloc
      logical :: getndp, diff, ndiff
      real :: vol_cv, lx_tran, lxgp

      !ewrite( 3, * ) 'In QUAD_1D_SHAPE'
      !ewrite( 3, * ) 'cv_ngi, cv_nloc', cv_ngi, cv_nloc

      nquad = cv_ngi / cv_nloc

      ! Allocating memory
      allocate( lx( nquad ) )
      allocate( wei( nquad ) )
      allocate( rdummy( max( cv_ngi, u_nloc, cv_nloc ) ) )
      allocate( xi_min( nquad ) )
      allocate( xi_max( nquad ) )
      allocate( cv_nodpos( cv_nloc ) )
      allocate( u_nodpos( u_nloc ) )

      if( ( cv_ngi /= 18 ) .and. ( cv_ngi /= 15 ) .and. ( cv_ngi /= 12 ) .and. &
           ( cv_ngi /= 9 ) .and. ( cv_ngi /= 6 ) .and. ( cv_ngi /= 1 ) ) &
           FLAbort("Incorrect number of quadrature points")

      Case_CV_NLOC: Select Case( cv_nloc )
      case( 4 );
         ! Node 1
         xi_min( 1 ) = -1.
         xi_max( 1 ) = -1. + 1. / 3.
         ! Node 2
         xi_min( 2 ) = -1. + 1. / 3.
         xi_max( 2 ) =  0.
         ! Node 3
         xi_min( 3 ) =  0.
         xi_max( 3 ) =  1. - 1. / 3.
         ! Node 4
         xi_min( 4 ) =  1. - 1. / 3.
         xi_max( 4 ) =  1.
      case( 3 );
         if( .true. ) then
            ! Node 1
            xi_min( 1 ) = -1.
            xi_max( 1 ) = -0.5
            ! Node 2
            xi_min( 2 ) = -0.5
            xi_max( 2 ) =  0.5
            ! Node 3
            xi_min( 3 ) =  0.5
            xi_max( 3 ) =  1.
         else
            ! Node 1
            xi_min( 1 ) = -1.
            xi_max( 1 ) = -1. + 2. / 3.
            ! Node 2
            xi_min( 2 ) = -1. + 2. / 3.
            xi_max( 2 ) =  1. - 2. / 3.
            ! Node 3
            xi_min( 3 ) =  1. - 2. / 3.
            xi_max( 3 ) =  1.
         end if
      case( 2 );
         ! Node 1
         xi_min( 1 ) = -1.
         xi_max( 1 ) =  0.
         ! Node 2
         xi_min( 2 ) =  0.
         xi_max( 2 ) =  1.
      case( 1 )
         ! Node 1
         xi_min( 1 ) = -1.
         xi_max( 1 ) =  1.
      end Select Case_CV_NLOC

      diff = .true.
      ndiff = .false.
      ! Find the roots of the quadrature points and nodes
      ! also get the weights.
      getndp = .true.

      !     Compute standard Gauss quadrature: weights and points
      call lagrot( rdummy, cv_nodpos, cv_nloc, getndp )
      call lagrot( rdummy, u_nodpos, u_nloc, getndp )

      getndp = .false.

      !     Compute standard Gauss quadrature: weights and points
      call lagrot( wei, lx, nquad, getndp )

      !ewrite(3,*)'cv_nloc, u_nloc, nquad, cv_ngi: ', cv_nloc, u_nloc, nquad, cv_ngi
      !ewrite(3,*)'wei: ', wei
      !ewrite(3,*)'lx: ', lx
      !ewrite(3,*)'cv_nodpos: ', cv_nodpos
      !ewrite(3,*)'u_nodpos: ', u_nodpos

      gpoi = 0
      cvn = 0.
      Loop_ICV1: do icv = 1, cv_nloc
         vol_cv = xi_max( icv ) - xi_min( icv )
         Loop_P1: do p = 1, nquad
            gpoi = gpoi + 1
            cvn( icv, gpoi ) = 1. ! Mapping to a new local coordinate system
            lx_tran = 0.5 * ( xi_max( icv ) + xi_min( icv ) ) + &
                 0.5 * ( xi_max( icv ) + xi_min( icv ) ) * lx( p )
            lxgp = lx_tran
            cvweigh( gpoi ) = wei( p ) * vol_cv / 2.
            Loop_ILOC1: do iloc = 1, cv_nloc
               n( iloc, gpoi ) = lagran( ndiff, lxgp, iloc, cv_nloc, cv_nodpos )
               nlx( iloc, gpoi ) = lagran( diff, lxgp, iloc, cv_nloc, cv_nodpos )
            end do Loop_ILOC1
            Loop_ILOC2: do iloc = 1, u_nloc
               un( iloc, gpoi ) = lagran( ndiff, lxgp, iloc, u_nloc, u_nodpos )
               unlx( iloc, gpoi ) = lagran( diff, lxgp, iloc, u_nloc, u_nodpos )
            end do Loop_ILOC2
         end do Loop_P1
      end do Loop_ICV1

      deallocate( lx )
      deallocate( wei )
      deallocate( rdummy )
      deallocate( xi_min )
      deallocate( xi_max )
      deallocate( cv_nodpos )
      deallocate( u_nodpos )

      !ewrite(3,*) 'Leaving QUAD_1D_SHAPE'

      return
    end subroutine quad_1d_shape

!!$ ===============
!!$ ===============

    subroutine quad_nd_shape( ndim, cv_ele_type, cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, &
         n, nlx, nly, nlz, &
         un, unlx, unly, unlz )
      ! For quadratic elements: Shape functions associated with volume integration
      ! using both CV (CVN) and FEM (N and its derivatives NLX/Y/Z) basis functions.

      implicit none
      integer, intent( in ) :: ndim, cv_ele_type, cv_ngi, cv_nloc, u_nloc
      real, dimension( :, : ), intent( inout ) :: cvn
      real, dimension( : ), intent( inout ) :: cvweigh
      real, dimension( :, : ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( :, : ), intent( inout ) :: un, unlx, unly, unlz

      ! Local variables
      real, dimension( :, : ), allocatable :: cvn_dum, cvn_1d_dum, n_1d, &
           nlx_1d, un_1d, unlx_1d, cvn_1d
      real, dimension( : ), allocatable :: cvweigh_dum, cvweigh_1d_dum, cvweigh_1d
      integer :: cv_ngi_1d, cv_nloc_1d, u_nloc_1d

      !ewrite(3,*) 'In MD_Shapes subrt: quad_nd_shape'
      !ewrite(3,*) 'cv_ngi, cv_nloc, u_nloc::', &
      !     cv_ngi, cv_nloc, u_nloc

      Conditional_Dimensionality: Select Case( ndim )
      case( 1 )
         cv_ngi_1d = cv_ngi
         cv_nloc_1d = cv_nloc
         u_nloc_1d = u_nloc

      case( 2 )
         cv_ngi_1d = int( sqrt( 1.e-3 + cv_ngi ))
         cv_nloc_1d = int( sqrt( 1.e-3 + cv_nloc ))
         u_nloc_1d = int( sqrt( 1.e-3 + u_nloc ))

      case( 3 )
         cv_ngi_1d = int( ( 1.e-3 + cv_ngi ) ** (1./3.) )
         cv_nloc_1d = int( ( 1.e-3 + cv_nloc ) ** (1./3.) )
         u_nloc_1d = int( ( 1.e-3 + u_nloc ) ** (1./3.) )

      case default; FLExit( " Invalid integer for NDIM " )

      end Select Conditional_Dimensionality

      ! Allocating memory
      allocate( cvweigh_1d( cv_ngi_1d ) )
      allocate( n_1d( cv_nloc_1d, cv_ngi_1d ) )
      allocate( nlx_1d( cv_nloc_1d, cv_ngi_1d ) )
      allocate( un_1d( u_nloc_1d, cv_ngi_1d ) )
      allocate( unlx_1d( u_nloc_1d, cv_ngi_1d ) )
      allocate( cvn_1d( cv_nloc_1d, cv_ngi_1d ) )

      ! Computing Shape Functions for 1D:
      call quad_1d_shape( cv_ngi_1d, cv_nloc_1d, u_nloc_1d, cvn_1d, cvweigh_1d, &
           n_1d, nlx_1d, un_1d, unlx_1d )

      ! Computing Shape Functions for N:
      call quad_nd_shape_N( cv_ele_type, ndim, cv_ngi, cv_nloc, cvn, cvweigh, &
           n, nlx, nly, nlz, &
           cv_ngi_1d, cv_nloc_1d, cvn_1d, cvweigh_1d, n_1d, nlx_1d )

      ! Allocating memory for UN calculation:
      allocate( cvn_dum( cv_nloc, cv_ngi ) )
      allocate( cvweigh_dum( cv_ngi ) )
      allocate( cvn_1d_dum( cv_nloc_1d, cv_ngi_1d ) )
      allocate( cvweigh_1d_dum( cv_ngi_1d ) )

      ! Computing Shape Functions for UN:
      call quad_nd_shape_N( cv_ele_type, ndim, cv_ngi, u_nloc, cvn_dum, cvweigh_dum, &
           un, unlx, unly, unlz, &
           cv_ngi_1d, u_nloc_1d, cvn_1d_dum, cvweigh_1d_dum, un_1d, unlx_1d )

      deallocate( cvn_dum )
      deallocate( cvn_1d_dum )
      deallocate( n_1d )
      deallocate( nlx_1d )
      deallocate( un_1d )
      deallocate( unlx_1d )
      deallocate( cvweigh_dum )
      deallocate( cvweigh_1d_dum )
      deallocate( cvweigh_1d )
      deallocate( cvn_1d )

      return
    end subroutine quad_nd_shape


!!$ ===============
!!$ ===============

    subroutine quad_nd_shape_N( cv_ele_type, ndim, cv_ngi, cv_nloc, cvn, cvweigh, &
         n, nlx, nly, nlz, &
         cv_ngi_1d, cv_nloc_1d, cvn_1d, cvweigh_1d, n_1d, nlx_1d )
      ! For quadatic elements -- shape functions associated with volume
      ! integration using both CV basis functions CVN as well as FEM basis
      ! functions N (and its derivatives NLX, NLY, NLZ)
      implicit none
      integer, intent( in ) :: cv_ele_type, ndim, cv_ngi, cv_nloc
      real, dimension( :, : ), intent( inout ) :: cvn
      real, dimension( : ), intent( inout ) :: cvweigh
      real, dimension( :, : ), intent( inout ) :: n, nlx, nly, nlz
      integer, intent( in ) :: cv_ngi_1d, cv_nloc_1d
      real, dimension( :, : ), intent( in ) :: cvn_1d
      real, dimension( : ), intent( in ) :: cvweigh_1d
      real, dimension( :, : ), intent( in ) :: n_1d, nlx_1d
      ! Local variables
      integer :: cv_iloc_1d, cv_jloc_1d, cv_kloc_1d, cv_iloc, cv_jloc, cv_kloc, cv_igi, &
           cv_igi_1d, cv_jgi_1d, cv_kgi_1d

      Conditional_Dimensionality: Select Case( ndim )
      case( 1 )
         cvn = cvn_1d
         n = n_1d
         nlx = nlx_1d
         nly = 0.
         nlz = 0.
         cvweigh = cvweigh_1d

      case( 2 )
         cv_iloc = 0
         nlz = 0.
         Loop_CVILOC_1D_2D: do cv_iloc_1d = 1, cv_nloc_1d
            Loop_CVJLOC_1D_2D: do cv_jloc_1d = 1, cv_nloc_1d
               cv_iloc = cv_iloc + 1
               cv_igi = 0
               do cv_igi_1d = 1, cv_nloc_1d
                  do cv_jgi_1d = 1, cv_nloc_1d
                     cv_igi = cv_igi + 1
                     cvn( cv_iloc, cv_igi ) =  &
                          cvn_1d( cv_iloc_1d, cv_igi_1d ) * cvn_1d( cv_jloc_1d, cv_jgi_1d )
                     n( cv_iloc, cv_igi ) =  &
                          n_1d( cv_iloc_1d, cv_igi_1d ) * n_1d( cv_jloc_1d, cv_jgi_1d )
                     nlx( cv_iloc, cv_igi ) =  &
                          nlx_1d( cv_iloc_1d, cv_igi_1d ) * n_1d( cv_jloc_1d, cv_jgi_1d )
                     nly( cv_iloc, cv_igi ) =  &
                          n_1d( cv_iloc_1d, cv_igi_1d ) * nlx_1d( cv_jloc_1d, cv_jgi_1d )
                     cvweigh( cv_igi ) = cvweigh_1d( cv_igi_1d ) * cvweigh_1d( cv_jgi_1d )
                  end do
               end do
            end do Loop_CVJLOC_1D_2D
         end do Loop_CVILOC_1D_2D

      case( 3 )
         cv_iloc = 0
         Loop_CVILOC_1D_3D: do cv_iloc_1d = 1, cv_nloc_1d
            Loop_CVJLOC_1D_3D: do cv_jloc_1d = 1, cv_nloc_1d
               Loop_CVKLOC_1D_3D: do cv_kloc_1d = 1, cv_nloc_1d
                  cv_iloc = cv_iloc + 1
                  cv_igi = 0
                  do cv_igi_1d = 1, cv_nloc_1d
                     do cv_jgi_1d = 1, cv_nloc_1d
                        do cv_kgi_1d = 1, cv_nloc_1d
                           cv_igi = cv_igi + 1
                           cvn( cv_iloc, cv_igi ) =  cvn_1d( cv_iloc_1d, cv_igi_1d ) * &
                                cvn_1d( cv_jloc_1d, cv_jgi_1d ) * cvn_1d( cv_kloc_1d, cv_kgi_1d )
                           n( cv_iloc, cv_igi ) =  n_1d( cv_iloc_1d, cv_igi_1d ) * &
                                n_1d( cv_jloc_1d, cv_jgi_1d ) * n_1d( cv_kloc_1d, cv_kgi_1d )
                           nlx( cv_iloc, cv_igi ) =  nlx_1d( cv_iloc_1d, cv_igi_1d ) * &
                                n_1d( cv_jloc_1d, cv_jgi_1d ) * n_1d( cv_kloc_1d, cv_kgi_1d )
                           nly( cv_iloc, cv_igi ) =  n_1d( cv_iloc_1d, cv_igi_1d ) * &
                                nlx_1d( cv_jloc_1d, cv_jgi_1d ) * n_1d( cv_kloc_1d, cv_kgi_1d )
                           nlz( cv_iloc, cv_igi ) =  n_1d( cv_iloc_1d, cv_igi_1d ) * &
                                n_1d( cv_jloc_1d, cv_jgi_1d ) * nlx_1d( cv_kloc_1d, cv_kgi_1d )
                           cvweigh( cv_igi ) = cvweigh_1d( cv_igi_1d ) * &
                                cvweigh_1d( cv_jgi_1d ) * cvweigh_1d( cv_kgi_1d )
                        end do
                     end do
                  end do
               end do Loop_CVKLOC_1D_3D
            end do Loop_CVJLOC_1D_3D
         end do Loop_CVILOC_1D_3D

      case default; FLExit( " Invalid integer for NDIM " )

      end Select Conditional_Dimensionality

      return
    end subroutine quad_nd_shape_N


!!$ ===============
!!$ ===============

    subroutine new_pt_qua_vol_cv_tri_tet_shape( cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, &
         n, nlx, nly, nlz, &
         un, unlx, unly, unlz )
      ! new 1 or 4 pt quadrature set on each CV of a quadratic tetrahedra .....
      ! Compute shape functions N, UN etc for linear trianles. Shape functions
      ! associated with volume integration using both CV basis functions CVN, as
      ! well as FEM basis functions N (and its derivatives NLX, NLY, NLZ). Also
      ! for velocity basis functions UN, UNLX, UNLY, UNLZ.
      implicit none
      integer, intent( in ) :: cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc
      real, dimension( :, : ), intent( inout ) :: cvn
      real, dimension( : ), intent( inout ) :: cvweigh
      real, dimension( :, : ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( :, : ), intent( inout ) :: un, unlx, unly, unlz
      ! Local variables
      integer, dimension( : ), allocatable :: x_ndgln, fem_nod, x_ndgln_ideal
      real, dimension( : ), allocatable :: lx, ly, lz, x, y, z, cvweigh_dummy, &
           x_ideal, y_ideal, z_ideal
      real, dimension( : ), allocatable :: quad_l1, quad_l2, quad_l3, quad_l4, &
           rdummy
      integer, parameter :: max_totele = 10000, max_x_nonods = 10000
      logical :: d1, dcyl, d3
      integer :: ele, quad_cv_ngi, quad_cv_nloc, totele, x_nonods, &
           cv_gj, cv_gk, cv_iloc, cv_gi, totele_sub, gi_max, gi, ngi_in_a_cv
      real :: rsum, rmax

      REAL, PARAMETER :: Pi=atan(1.0)*4.0

      real, parameter :: pt1=23./288., pt2=75./288., pt3=167./288., pt4=219./288.
      real, parameter :: pt5=35./96., pt6=13./96.

      real, parameter :: w1=1./192., w2=7./288., w3=1.0/72.0

      real, parameter :: ptA=((pt2+pt3)*w1+pt5*w3)/w2
      real, parameter :: ptB=(2.0*pt1*w1+pt6*w3)/w2

      !ewrite(3,*)'In vol_cv_tri_tet_shape'

      if(cv_ngi.ne.10) then
         print *,'the wrong number of quadrature points'
      endif

      d1 = ( ndim == 1 )
      dcyl = .false.
      d3 = ( ndim == 3 )

      if( d3 ) then
         Select Case( cv_nloc )
         case( 4 ) ;  quad_cv_nloc = 8  ! Linear tetrahedron
         case( 10 ) ; quad_cv_nloc = 27 ! Quadratic tetrahedron
         case default; FLExit( "Wrong integer for CV_NLOC" )
         end Select
      else
         Select Case( cv_nloc )
         case( 3 ) ;  quad_cv_nloc = 4  ! Linear triangle
         case( 6 ) ; quad_cv_nloc = 9 ! Quadratic triangle
         case default; FLExit( "Wrong integer for CV_NLOC" )
         end Select
      endif

      ! Allocating memory
      allocate( quad_l1( cv_ngi ) ) ; quad_l1 = 0.
      allocate( quad_l2( cv_ngi ) ) ; quad_l2 = 0.
      allocate( quad_l3( cv_ngi ) ) ; quad_l3 = 0.
      allocate( quad_l4( cv_ngi ) ) ; quad_l4 = 0.
      allocate( rdummy( cv_ngi ) ) ; rdummy = 0.

      if(NEW_HIGH_ORDER_VOL_QUADRATIC_ELE_QUADRATURE) then
         !****************James insert new quadrature set here****************start
         stop 2999
         !****************James insert new quadrature set here****************end
      else
         !
         !Volumetric points
         !Node :  quadrature point  :   volume
         !A: : [ 0.07986111  0.07986111  0.07986111  0.76041667] 0.005208
         !B: : [ 0.38839286  0.11160714  0.11160714  0.38839286] 0.024306
         !C: : [ 0.76041667  0.07986111  0.07986111  0.07986111] 0.005208
         !D: : [ 0.11160714  0.38839286  0.11160714  0.38839286] 0.024306
         !E: : [ 0.38839286  0.38839286  0.11160714  0.11160714] 0.024306
         !F: : [ 0.07986111  0.76041667  0.07986111  0.07986111] 0.005208
         !G: : [ 0.11160714  0.11160714  0.38839286  0.38839286] 0.024306
         !H: : [ 0.38839286  0.11160714  0.38839286  0.11160714] 0.024306
         !I: : [ 0.11160714  0.38839286  0.38839286  0.11160714] 0.024306
         !J: : [ 0.07986111  0.07986111  0.76041667  0.07986111] 0.005208

         quad_l1(1) = pt4;         quad_l2(1) = pt1;         quad_l3(1) = pt1;         quad_l4(1) = pt1
         quad_l1(2) = ptA;         quad_l2(2) = ptA;         quad_l3(2) = ptB;         quad_l4(2) = ptB
         quad_l1(3) = pt1;         quad_l2(3) = pt4;         quad_l3(3) = pt1;         quad_l4(3) = pt1
         quad_l1(4) = ptA;         quad_l2(4) = ptB;         quad_l3(4) = ptA;         quad_l4(4) = ptB
         quad_l1(5) = ptB;         quad_l2(5) = ptA;         quad_l3(5) = ptA;         quad_l4(5) = ptB
         quad_l1(6) = pt1;         quad_l2(6) = pt1;         quad_l3(6) = pt4;         quad_l4(6) = pt1
         quad_l1(7) = ptA;         quad_l2(7) = ptB;         quad_l3(7) = ptB;         quad_l4(7) = ptA
         quad_l1(8) = ptB;         quad_l2(8) = ptA;         quad_l3(8) = ptB;         quad_l4(8) = ptA
         quad_l1(9) = ptB;         quad_l2(9) = ptB;         quad_l3(9) = ptA;         quad_l4(9) = ptA
         quad_l1(10)= pt1;         quad_l2(10)= pt1;         quad_l3(10)= pt1;         quad_l4(10)= pt4

         cvweigh(1) = w1
         cvweigh(2) = w2
         cvweigh(3) = w1
         cvweigh(4) = w2
         cvweigh(5) = w2
         cvweigh(6) = w1
         cvweigh(7) = w2
         cvweigh(8) = w2
         cvweigh(9) = w2
         cvweigh(10)= w1
         ! endof if(NEW_HIGH_ORDER_VOL_QUADRATIC_ELE_QUADRATURE) then
      endif

      ! scale taking into account we have a volume of 1./6. of the tet
      !  cvweigh(:) = cvweigh(:) * 6.0


      ! Now determine the basis functions and derivatives at the
      ! quadrature pts quad_L1, quad_L2, quad_L3, quad_L4, etc

      !  call SHATRIold(quad_L4, quad_L1,quad_L2, quad_L3, rdummy, .true. , &
      !        CV_NLOC,cv_NGI, N,NLX,NLY,NLZ)

      call shatri_hex( quad_l1, quad_l2, quad_l3, quad_l4, rdummy, .true., &
           cv_nloc, cv_ngi, n, nlx, nly, nlz, &
           .true. )

      ! calculate cvn based on maximum value of n:
      !      cvn=0.0
      !      do cv_iloc=1,cv_nloc
      !         rmax=-1.e+10
      !         gi_max=0
      !         do gi=1,cv_ngi
      !            if(n(cv_iloc,gi).gt.rmax) then
      !               rmax=n(cv_iloc,gi)
      !               gi_max=gi
      !            endif
      !         end do
      !         cvn(cv_iloc,gi_max)=1.0
      !      end do
      ngi_in_a_cv=cv_ngi/cv_nloc
      cvn=0.0
      do cv_iloc=1,cv_nloc
         cvn(cv_iloc,   1+(cv_iloc-1)*ngi_in_a_cv : cv_iloc*ngi_in_a_cv) = 1.0
      end do



      ! Now determine the basis functions and derivatives at the
      ! quadrature pts quad_L1, quad_L2, quad_L3, quad_L4, etc
      call shatri_hex( quad_l1, quad_l2, quad_l3, quad_l4, rdummy, d3, &
           u_nloc, cv_ngi, un, unlx, unly, unlz, &
           .true. )

      return
    end subroutine new_pt_qua_vol_cv_tri_tet_shape


!!$ ===============
!!$ ===============

    subroutine vol_cv_tri_tet_shape( cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, &
         n, nlx, nly, nlz, &
         un, unlx, unly, unlz )
      ! Compute shape functions N, UN etc for linear trianles. Shape functions
      ! associated with volume integration using both CV basis functions CVN, as
      ! well as FEM basis functions N (and its derivatives NLX, NLY, NLZ). Also
      ! for velocity basis functions UN, UNLX, UNLY, UNLZ.
      implicit none
      integer, intent( in ) :: cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc
      real, dimension( :, : ), intent( inout ) :: cvn
      real, dimension( : ), intent( inout ) :: cvweigh
      real, dimension( :, : ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( :, : ), intent( inout ) :: un, unlx, unly, unlz
      ! Local variables
      integer, dimension( : ), allocatable :: x_ndgln, fem_nod, x_ndgln_ideal
      real, dimension( : ), allocatable :: lx, ly, lz, x, y, z, cvweigh_dummy, &
           x_ideal, y_ideal, z_ideal
      integer, parameter :: max_totele = 10000, max_x_nonods = 10000
      logical :: d1, dcyl, d3
      integer :: ele, quad_cv_ngi, quad_cv_nloc, totele, x_nonods, &
           cv_gj, cv_gk, cv_iloc, cv_gi, totele_sub
      real :: rsum

      !ewrite(3,*)'In vol_cv_tri_tet_shape'

      d1 = ( ndim == 1 )
      dcyl = .false.
      d3 = ( ndim == 3 )

      if( d3 ) then
         Select Case( cv_nloc )
         case( 4 ) ;  quad_cv_nloc = 8  ! Linear tetrahedron
         case( 10 ) ; quad_cv_nloc = 27 ! Quadratic tetrahedron
         case default; FLExit( "Wrong integer for CV_NLOC" )
         end Select
      else
         Select Case( cv_nloc )
         case( 3 ) ;  quad_cv_nloc = 4  ! Linear triangle
         case( 6 ) ; quad_cv_nloc = 9 ! Quadratic triangle
         case default; FLExit( "Wrong integer for CV_NLOC" )
         end Select
      endif

      ! Allocating memory
      allocate( lx( max_x_nonods ) ) ; lx = 0.
      allocate( ly( max_x_nonods ) ) ; ly = 0.
      allocate( lz( max_x_nonods ) ) ; lz = 0.
      allocate( x( max_x_nonods )) ; x = 0.
      allocate( y( max_x_nonods )) ; y =0.
      allocate( z( max_x_nonods )) ; z = 0.
      allocate( fem_nod( max_x_nonods ) ) ; fem_nod = 0
      allocate( x_ndgln( max_totele * quad_cv_nloc ) ) ; x_ndgln = 0
      allocate( cvweigh_dummy( cv_ngi )) ; cvweigh_dummy = 0.
      allocate( x_ideal( max_x_nonods ) ) ; x_ideal = 0.
      allocate( y_ideal( max_x_nonods ) ) ; y_ideal = 0.
      allocate( z_ideal( max_x_nonods ) ) ; z_ideal = 0.
      allocate( x_ndgln_ideal( max_x_nonods ) ) ; x_ndgln_ideal = 0

      ! Get the x_ndgln for the nodes of the triangle or tet or hex/quad super-elements:
      x = 0. ; y = 0. ; z = 0. ; lx = 0. ; ly = 0. ; lz = 0. ; fem_nod = 0 ; x_ndgln = 0
      call Compute_XNDGLN_TriTetQuadHex( cv_ele_type, &
           max_totele, max_x_nonods, quad_cv_nloc, &
           totele, x_nonods, &
           x_ndgln, lx, ly, lz, x, y, z, fem_nod, &
           x_ideal, y_ideal, z_ideal, x_ndgln_ideal )

      ! Compute the shape functions using these quadrilaterals/hexs:
      ! For pressure:
      call shape_tri_tet( cv_ele_type, cv_nloc, &
           cv_ele_type, ndim, totele, cv_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
           n, nlx, nly, nlz, cvweigh )

      ! Compute cvn:
      call Calc_CVN_TriTetQuadHex( cv_ele_type, totele, cv_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, fem_nod, cvn )

      !ewrite(3,*)'cvweigh:',cvweigh
      !do cv_iloc=1,cv_nloc
      !   rsum=0.0
      !   do cv_gi=1,cv_ngi
      !      rsum=rsum+cvn(cv_iloc,cv_gi)*cvweigh(cv_gi)
      !   end do
      !   ewrite(3,*)'cv_iloc,rsum:',cv_iloc,rsum
      !end do
      !stop 2922
      !if(cv_nloc==10) then
      !  !ewrite(3,*)'cv_nloc=',cv_nloc
      !   totele_sub=8
      !   call test_quad_tet( cv_nloc, cv_ngi, cvn, n, nlx, nly, nlz, &
      !        cvweigh, x_ideal, y_ideal, z_ideal, cv_nloc, x_ndgln_ideal, 1)
      !endif

      ! And for velocities:
      if(u_nloc==1) then ! constant basis function throughout element...
         un=1.0
         unlx=0.0
         unly=0.0
         unlz=0.0
      else
         call shape_tri_tet( cv_ele_type, cv_nloc, &
              cv_ele_type, ndim, totele, u_nloc, cv_ngi, x_nonods, &
              quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
              un, unlx, unly, unlz, cvweigh_dummy )
      endif

      deallocate( lx )
      deallocate( ly )
      deallocate( lz )
      deallocate( x )
      deallocate( y )
      deallocate( z )
      deallocate( x_ndgln )
      deallocate( fem_nod )
      deallocate( cvweigh_dummy )

      return
    end subroutine vol_cv_tri_tet_shape


!!$ ===============
!!$ ===============

  end module Quadratures
