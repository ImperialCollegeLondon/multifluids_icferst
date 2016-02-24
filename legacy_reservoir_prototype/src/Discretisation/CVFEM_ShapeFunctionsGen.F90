
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


  module CVFEM_ShapeFunctions
    use fldebug
    use multi_data_types

  contains

    subroutine cv_fem_shape_funs_new( shape_fun, Mdims, GIdims, cv_ele_type, QUAD_OVER_WHOLE_ELE)
      ! This subrt defines the sub-control volume and FEM shape functions.
      ! Shape functions associated with volume integration using both CV basis
      ! functions CVN as well as FEM basis functions CVFEN (and its derivatives
      ! CVFENLX, CVFENLY, CVFENLZ)
      implicit none
      integer , intent(in) :: cv_ele_type
      type (multi_dimensions), intent(in) :: Mdims
      type (multi_gi_dimensions), intent(in) :: GIdims
      type(multi_shape_funs), intent(inout) :: shape_fun
      logical, intent( in ) :: QUAD_OVER_WHOLE_ELE


      ! Local variables
!!$      real, dimension( : , : ), allocatable :: cvn2 ! dimension( cv_nloc, cv_ngi )
!!$      real, dimension(  : , : ), allocatable :: cvn_short2!dimension( cv_nloc, cv_ngi_short )
!!$      real, dimension( : ), allocatable :: cvweight2!dimension( cv_ngi )
!!$      real, dimension(  : , : ), allocatable :: cvfen2!dimension( cv_nloc, cv_ngi )
!!$      real, dimension( : ), allocatable :: cvweight_short2!dimension( cv_ngi_short )
!!$      real, dimension(  : , : ), allocatable :: cvfen_short2
!!$      real, dimension(:,:,:), allocatable :: cvfenlx_all2, cvfenlx_short_all2, ufenlx_all2, scvfenlx_all2,&
!!$           sufenlx_all2, sbcvfenlx_all2, sbufenlx_all2
!!$      real, dimension(  : , : ), allocatable :: ufen2!dimension( u_nloc, cv_ngi )
!!$      integer, dimension(  : , : ), allocatable :: cv_neiloc2!dimension( cv_nloc, scvngi )
!!$      real, dimension(  : , : ), allocatable :: scvfen2, scvfenslx2, scvfensly2!dimension( cv_nloc, scvngi )
!!$      real, dimension( : ), allocatable :: scvfeweigh2!dimension( scvngi )
!!$      real, dimension(  : , : ), allocatable :: sufen2, sufenslx2, sufensly2
!!$      real, dimension(  : , : ), allocatable :: sbcvn2!dimension( cv_snloc, sbcvngi )
!!$      real, dimension(  : , : ), allocatable :: sbcvfen2, sbcvfenslx2, sbcvfensly2! dimension( cv_snloc, sbcvngi )
!!$      real, dimension( : ), allocatable :: sbcvfeweigh2!dimension( sbcvngi )
!!$      real, dimension(  : , : ), allocatable :: sbufen2, sbufenslx2, sbufensly2!dimension( u_snloc, sbcvngi )
!!$      integer, dimension(  : , : ), allocatable :: cv_sloclist2!dimension( nface, cv_snloc )
!!$      integer, dimension(  : , : ), allocatable :: u_sloclist2!dimension( nface, u_snloc )
!!$      integer, dimension( : ), allocatable :: findgpts2!dimension( cv_nloc + 1 )
!!$      integer, dimension( : ), allocatable :: colgpts2!dimension( cv_nloc * scvngi )
!!$      integer :: ncolgpts2
!!$      real, dimension( : ), allocatable :: sele_overlap_scale2!dimension( cv_nloc )


      call cv_fem_shape_funs( &
           Mdims%ndim, cv_ele_type, &
           GIdims%cv_ngi, GIdims%cv_ngi, Mdims%cv_nloc, Mdims%u_nloc, cvn2, cvn_short2, &
                                ! Volume shape functions
           cvweight2, cvfen2, cvfenlx_all2(1,:,:), cvfenlx_all2(2,:,:), cvfenlx_all2(3,:,:), &
           cvweight_short2, cvfen_short2, cvfenlx_short_all2(1,:,:), cvfenlx_short_all2(2,:,:), cvfenlx_short_all2(3,:,:), &
           ufen2, ufenlx_all2(1,:,:), ufenlx_all2(2,:,:), ufenlx_all2(3,:,:), &
                                ! Surface of each CV shape functions
           GIdims%scvngi, cv_neiloc2, shape_fun%cv_on_face, shape_fun%cvfem_on_face, &
           scvfen2, scvfenslx2, scvfensly2, scvfeweigh2, &
           scvfenlx_all2(1,:,:), scvfenlx_all2(2,:,:), scvfenlx_all2(3,:,:), &
           sufen2, sufenslx2, sufensly2, &
           sufenlx_all2(1,:,:), sufenlx_all2(2,:,:), sufenlx_all2(3,:,:), &
                                ! Surface element shape funcs
           shape_fun%u_on_face, shape_fun%ufem_on_face, GIdims%nface, &
           GIdims%sbcvngi, sbcvn2, sbcvfen2, sbcvfenslx2, sbcvfensly2, sbcvfeweigh2, sbcvfenlx_all2(1,:,:),&
           sbcvfenlx_all2(2,:,:), sbcvfenlx_all2(3,:,:), &
           sbufen2, sbufenslx2, sbufensly2, sbufenlx_all2(1,:,:), sbufenlx_all2(2,:,:), sbufenlx_all2(3,:,:), &
           cv_sloclist2, u_sloclist2, Mdims%cv_snloc, Mdims%u_snloc, &
                                ! Define the gauss points that lie on the surface of the CV
           findgpts2, colgpts2, shape_fun%ncolgpts, &
           sele_overlap_scale2, QUAD_OVER_WHOLE_ELE) 


      !Copy values into new format

      shape_fun%cvn = cvn2; shape_fun%cvweight =  cvweight2
      shape_fun%cvfen = cvfen2; shape_fun%cvfenlx_all(1:Mdims%Ndim, :, :) = cvfenlx_all2(1:Mdims%Ndim, :, :)
      shape_fun%ufen = ufen2; shape_fun%ufenlx_all(1:Mdims%Ndim, :, :) = ufenlx_all2(1:Mdims%Ndim, :, :)
      shape_fun%cv_neiloc = cv_neiloc2; shape_fun%scvfen = scvfen2
      shape_fun%scvfenslx = scvfenslx2; shape_fun%scvfensly = scvfensly2
      shape_fun%scvfeweigh = scvfeweigh2; shape_fun%scvfenlx_all(1:Mdims%Ndim, :, :) = scvfenlx_all2(1:Mdims%Ndim, :, :)
      shape_fun%sufen = sufen2; shape_fun%sufenslx = sufenslx2
      shape_fun%sufensly = sufensly2; shape_fun%sufenlx_all(1:Mdims%Ndim, :, :) = sufenlx_all2(1:Mdims%Ndim, :, :)
      shape_fun%sbcvn = sbcvn2; shape_fun%sbcvfen = sbcvfen2
      shape_fun%sbcvfenslx = sbcvfenslx2; shape_fun%sbcvfensly = sbcvfensly2
      shape_fun%sbcvfeweigh = sbcvfeweigh2; shape_fun%sbcvfenlx_all(1:Mdims%Ndim, :, :) = sbcvfenlx_all2(1:Mdims%Ndim, :, :)
      shape_fun%sbufen = sbufen2; shape_fun%sbufenslx = sbufenslx2
      shape_fun%sbufensly = sbufensly2; shape_fun%sbufenlx_all(1:Mdims%Ndim, :, :) = sbufenlx_all2(1:Mdims%Ndim, :, :)
      shape_fun%cv_sloclist = cv_sloclist2; shape_fun%u_sloclist = u_sloclist2
      shape_fun%findgpts = findgpts2; shape_fun%colgpts =   colgpts2


      ewrite(3,*) 'in  cv_fem_shape_funs subrt'
      sele_overlap_scale = 1.

      if( QUAD_OVER_WHOLE_ELE ) then ! integrate over whole element
         ewrite(3,*)'2 going into SHAPE_one_ele'
         call SHAPE_one_ele2(&
              Mdims%ndim, cv_ele_type, &
              GIdims%cv_ngi, Mdims%cv_nloc, Mdims%u_nloc, &
!!$ Volume shape functions
              shape_fun%cvweight, shape_fun%cvfen, shape_fun%cvfenlx_all(1, :, :), &
              shape_fun%cvfenlx_all(2, :, :), shape_fun%cvfenlx_all(3, :, :), &
              shape_fun%ufen, shape_fun%ufenlx_all(1, :, :), &
              shape_fun%ufenlx_all(2, :, :), shape_fun%ufenlx_all(3, :, :), &
!!$ Surface of each CV shape functions
              GIdims%sbcvngi,  &
              shape_fun%sbcvfen, shape_fun%sbcvfenslx, shape_fun%sbcvfensly, &
              shape_fun%sbcvfeweigh, &
              shape_fun%sbufen, shape_fun%sbufenslx, shape_fun%sbufensly, &
!!$ Surface element shape functions
              GIdims%nface, &
              shape_fun%cv_sloclist, shape_fun%u_sloclist, Mdims%cv_snloc, Mdims%u_snloc )

         if( GIdims%scvngi /= GIdims%sbcvngi) FLAbort("scvngi/=sbcvngi")

      else
         call shape_cv_n( Mdims%ndim, cv_ele_type, &
              GIdims%cv_ngi, Mdims%cv_nloc, Mdims%u_nloc, shape_fun%cvn, shape_fun%cvweight, &
              shape_fun%cvfen, shape_fun%cvfenlx_all(1, :, :), shape_fun%cvfenlx_all(2, :, :), &
              shape_fun%cvfenlx_all(3, :, :), &
              shape_fun%ufen, shape_fun%ufenlx_all(1, :, :), shape_fun%ufenlx_all(2, :, :), &
              shape_fun%ufenlx_all(3, :, :) )
      endif



    end subroutine cv_fem_shape_funs_new


    SUBROUTINE SHAPE_one_ele2( cv_ele_type, Mdims, GIdims, &
         shape_fun, cv_sloclist, u_sloclist )
!!$
!!$ This subrt defines the sub-control volume and FEM shape functions.
!!$    Shape functions associated with volume integration using both CV basis
!!$    functions CVN as well as FEM basis functions CVFEN (and its derivatives
!!$    CVFENLX, CVFENLY, CVFENLZ)
!!$
      implicit none
      integer , intent(in) :: cv_ele_type
      type (multi_dimensions), intent(in) :: Mdims
      type (multi_gi_dimensions), intent(in) :: GIdims
      type(multi_shape_funs), intent(inout) :: shape_fun
      integer, dimension( :, : ), intent( inout ) :: cv_sloclist
      integer, dimension( :, : ), intent( inout ) :: u_sloclist


!!$ Local variables
      integer :: MLOC, SMLOC, NWICEL
      logical :: LOWQUA
      real, dimension( :, : ), allocatable :: M, MLX, MLY, MLZ, SM, SMLX, SMLY
      real :: temp_rub(1000)

      ewrite(3,*)'just inside SHAPE_one_ele'

      LOWQUA = .false. ; MLOC = 1 ; SMLOC = 1
      ALLOCATE( M( MLOC, GIdims%cv_ngi ),  MLX( MLOC, GIdims%cv_ngi ), MLY( MLOC, GIdims%cv_ngi ), &
           MLZ( MLOC, GIdims%cv_ngi ), SM( SMLOC, GIdims%sbcvngi ), SMLX( SMLOC, GIdims%sbcvngi ), &
           SMLY( SMLOC, GIdims%sbcvngi ) )

!!$ For pressure:
      nwicel = Get_NwiCel( Mdims%ndim == 3, Mdims%cv_nloc )
      call shape( lowqua, GIdims%cv_ngi, Mdims%cv_nloc, mloc, GIdims%sbcvngi, Mdims%cv_snloc, smloc, &
           m, mlx, mly, mlz, shape_fun%cvweight, shape_fun%cvfen, shape_fun%cvfenlx_all(1, :, :), &
           shape_fun%cvfenlx_all(2, :, :), shape_fun%cvfenlx_all(3, :, :), &
           shape_fun%sbcvfeweigh, shape_fun%sbcvfen, shape_fun%sbcvfenslx, shape_fun%sbcvfensly, &
           sm, smlx, smly, &
           nwicel, Mdims%ndim == 3 )

!!$ For velocity:
      nwicel = Get_NwiCel( Mdims%ndim == 3, Mdims%u_nloc )
      call shape( lowqua, GIdims%cv_ngi, Mdims%u_nloc, mloc, GIdims%sbcvngi, Mdims%u_snloc, smloc, &
           m, mlx, mly, mlz, temp_rub, shape_fun%ufen, shape_fun%ufenlx_all(1, :, :), &
           shape_fun%ufenlx_all(2, :, :), shape_fun%ufenlx_all(3, :, :), &
           temp_rub, shape_fun%sbufen, shape_fun%sbufenslx, shape_fun%sbufensly, &
           sm, smlx, smly, &
           nwicel, Mdims%ndim == 3 )

!!$ Determine CV_SLOCLIST & U_SLOCLIST
      call determin_sloclist( cv_sloclist, Mdims%cv_nloc, Mdims%cv_snloc, GIdims%nface, &
           Mdims%ndim, cv_ele_type )

      if( Mdims%u_snloc == 1 )then
         u_sloclist( 1, 1 ) = 1 ; u_sloclist( 2, 1 ) = Mdims%u_nloc
      else
         call determin_sloclist( u_sloclist, Mdims%u_nloc, Mdims%u_snloc, GIdims%nface, &
              Mdims%ndim, cv_ele_type )
      end if

      if ( Mdims%ndim < 3 ) then
         shape_fun%cvfenlx_all(3, :, :) = 0. ; shape_fun%ufenlx_all(3, :, :) = 0. ; &
              shape_fun%sbcvfensly = 0. ; shape_fun%sbufensly = 0.
      end if
      if ( Mdims%ndim < 2 ) then
         shape_fun%cvfenlx_all(2, :, :) = 0. ; shape_fun%ufenlx_all(2, :, :) = 0. ; &
              shape_fun%sbcvfenslx = 0. ; shape_fun%sbufenslx = 0.
      end if

      deallocate( M, MLX, MLY, MLZ, SM, SMLX, SMLY )

      return
    END SUBROUTINE SHAPE_one_ele2


  subroutine shape_cv_n( cv_ele_type, Mdims, GIdims, shape_fun )
!!$
!!$ Shape functions associated with volume integration using both CV basis
!!$    functions CVN as well as FEM basis functions N (and its derivatives NLX, NLY, NLZ)
!!$    also for velocity basis functions UN, UNLX, UNLY, UNLZ
!!$
      implicit none
      integer , intent(in) :: cv_ele_type
      type (multi_dimensions), intent(in) :: Mdims
      type (multi_gi_dimensions), intent(in) :: GIdims
      type(multi_shape_funs), intent(inout) :: shape_fun



!!$       ndim, cv_ele_type, &
!!$       cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, &
!!$       n, nlx, nly, nlz, &
!!$       un, unlx, unly, unlz )

    implicit none
    integer, intent( in ) :: ndim, cv_ele_type, cv_ngi, cv_nloc, u_nloc
    real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: cvn
    real, dimension( cv_ngi ), intent( inout ) :: cvweigh
    real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: n, nlx, nly, nlz
    real, dimension( u_nloc, cv_ngi ), intent( inout ) :: un, unlx, unly, unlz

    ! new quadratic element quadrature by James and Zhi and Chris:


    ewrite(3,*) 'In SHAPE_CV_N'

    Select Case( cv_ele_type )
    case( 1, 2 ) ! 1D
       call quad_1d_shape( cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, n, nlx, un, unlx )
       nly = 0.
       nlz = 0.
       unly = 0.
       unlz = 0.

    case( 5, 6, 9, 10 ) ! Quadrilaterals and Hexahedra
       call quad_nd_shape( ndim, cv_ele_type, cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, &
            n, nlx, nly, nlz, &
            un, unlx, unly, unlz )

    case( 3, 4, 7, 8 ) ! Triangles and Tetrahedra
       if( new_quadratic_ele_quadrature .and. cv_ele_type==8) then
          call new_pt_qua_vol_cv_tri_tet_shape( cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc, cvn, &
               cvweigh, n, nlx, nly, nlz, &
               un, unlx, unly, unlz )
       else
          call vol_cv_tri_tet_shape( cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc, cvn, &
               cvweigh, n, nlx, nly, nlz, &
               un, unlx, unly, unlz )
       endif
       !stop 12
    case default; FLExit( "Wrong integer for CV_ELE_TYPE" )
    end Select

    !ewrite(3,*) 'Leaving SHAPE_CV_N'
return
  end subroutine shape_cv_n



  end module CVFEM_ShapeFunctions
