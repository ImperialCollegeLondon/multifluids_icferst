
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
    use Quadratures

  contains

    subroutine CVFEM_ShapeFunctions( shape_fun, Mdims, GIdims, cv_ele_type, QUAD_OVER_WHOLE_ELE)
!!$ This subrt defines the sub-control volume and FEM shape functions. Shape functions associated 
!!$      with volume integration using both CV basis functions CVN as well as FEM basis functions 
!!$      CVFEN (and its derivatives CVFENLX, CVFENLY, CVFENLZ).
!!$
!!$      (a) scvfen( cv_nloc, scvngi ): the shape function evaluated for each node
!!$               at each surface gauss point
!!$      (b) scvfenslx[y/z]( cv_nloc, scvngi ): the surface derivatives of the shape
!!$               function for each node at those same points, and the derivatives
!!$               of the shape
!!$      (c) scvfeweigh( scvngi ): the Gauss weights to use when integrating around
!!$               the control volume surface
!!$      (d) cv_neiloc( cv_nloc, scvngi ): neighbour node for a given node/gauss-point
!!$               pair. This also include quadature points around the element.
!!$ 
      implicit none
      integer , intent(in) :: cv_ele_type
      type (multi_dimensions), intent(in) :: Mdims
      type (multi_gi_dimensions), intent(in) :: GIdims
      type(multi_shape_funs), intent(inout) :: shape_fun
      logical, intent( in ) :: QUAD_OVER_WHOLE_ELE


      ewrite(3,*) 'in  cv_fem_shape_funs subrt'
      sele_overlap_scale = 1.

      if( QUAD_OVER_WHOLE_ELE ) then ! integrate over whole element
         call shape_one_ele2( cv_ele_type, Mdims, GIdims, shape_fun )
         if( GIdims%scvngi /= GIdims%sbcvngi) FLAbort("scvngi/=sbcvngi")

      else ! Not integrate over whole element
         call shape_cv_n( cv_ele_type, Mdims, GIdims, shape_fun )

!!$ Determine the surface element shape functions from those calculated in 
!!$   SHAPESV_FEM_PLUS and also CV_SLOCLIST( NFACE,CV_SNLOC )

         call shapesv_fem_plus( GIdims%scvngi, shape_fun%cv_neiloc, shape_fun%cv_on_face, &
              shape_fun%cvfem_on_face, &
              shape_fun%ufem_on_face, &
              cv_ele_type, Mdims%cv_nloc, shape_fun%scvfen, shape_fun%scvfenslx, shape_fun%scvfensly, &
              shape_fun%scvfeweigh, &
              shape_fun%scvfenlx_all(1, :, :), shape_fun%scvfenlx_all(2, :, :), &
              shape_fun%scvfenlx_all(3, :, :), &
              Mdims%u_nloc, shape_fun%sufen, shape_fun%sufenslx, shape_fun%sufensly, &
              shape_fun%sufenlx_all(1, :, :), shape_fun%sufenlx_all(2, :, :), &
              shape_fun%sufenlx_all(3, :, :), &
              Mdims%ndim )

         call det_suf_ele_shape( GIdims%scvngi, GIdims%nface, &
              shape_fun%cvfem_on_face, &
              Mdims%cv_nloc, shape_fun%scvfen, shape_fun%scvfenslx, shape_fun%scvfensly, &
              shape_fun%scvfeweigh, &
              shape_fun%scvfenlx_all(1, :, :), shape_fun%scvfenlx_all(2, :, :), &
              shape_fun%scvfenlx_all(3, :, :), &
              Mdims%u_nloc, shape_fun%sufen, shape_fun%sufenslx, shape_fun%sufensly, &
              shape_fun%sufenlx_all(1, :, :), shape_fun%sufenlx_all(2, :, :), &
              shape_fun%sufenlx_all(3, :, :), &
              GIdims%sbcvngi, shape_fun%sbcvfen, shape_fun%sbcvfenslx, shape_fun%sbcvfensly, &
              shape_fun%sbcvfeweigh, &
              shape_fun%sbcvfenlx_all(1, :, :), shape_fun%sbcvfenlx_all(2, :, :), &
              shape_fun%sbcvfenlx_all(3, :, :), &
              shape_fun%sbufen, shape_fun%sbufenslx, shape_fun%sbufensly, &
              shape_fun%sbufenlx_all(1, :, :), shape_fun%sbufenlx_all(2, :, :), &
              shape_fun%sbufenlx_all(3, :, :), & 
              shape_fun%cv_sloclist, shape_fun%u_sloclist, Mdims%cv_snloc, Mdims%u_snloc, &
              Mdims%ndim, cv_ele_type )

!!$ Define Gauss points that lie on the surface of the control volume surrounding 
!!$   a given local node (iloc) that is FINDGPTS, COLGPTS, NCOLGPTS
!!$ 
         call gaussiloc( shape_fun%findgpts, shape_fun%colgpts, shape_fun%ncolgpts, &
              shape_fun%cv_neiloc, Mdims%cv_nloc, GIdims%scvngi )

      endif


!!$ Set to zero anything that should be zero in case it was not pre-defined
      if( Mdims%ndim < 2 ) then
         shape_fun%cvfenlx_all(2, :, :) = 0. ; shape_fun%ufenlx_all(2, :, :) = 0.    ; &
              shape_fun%scvfenslx = 0.       ; shape_fun%scvfenlx_all(2, :, :) = 0.  ; &
              shape_fun%sufenslx = 0.        ; shape_fun%sufenlx_all(2, :, :) = 0.   ; &
              shape_fun%sbcvfenslx = 0.      ; shape_fun%sbcvfenlx_all(2, :, :) = 0. ; &
              shape_fun%sbufenslx = 0.       ; shape_fun%sbufenlx_all(2, :, :) = 0. 

      elseif( Mdims%ndim < 3 ) then
         shape_fun%cvfenlx_all(3, :, :) = 0. ; shape_fun%ufenlx_all(3, :, :) = 0.    ; &
              shape_fun%scvfensly = 0.       ; shape_fun%scvfenlx_all(3, :, :) = 0.  ; &
              shape_fun%sufensly = 0.        ; shape_fun%sufenlx_all(3, :, :) = 0.   ; &
              shape_fun%sbcvfensly = 0.      ; shape_fun%sbcvfenlx_all(3, :, :) = 0. ; &
              shape_fun%sbufensly = 0.       ; shape_fun%sbufenlx_all(3, :, :) = 0. 

      end if

!!$ Calculate sbcvn from sbcvfen - Use the max scvfen at a quadrature pt and set to 1:
      shape_fun%sbcvn = 0.0
      do sgi = 1, GIdims%sbcvngi
         cv_skloc = maxloc( shape_fun%sbcvfen( :, sgi) )
         shape_fun%sbcvn( cv_skloc(1), sgi ) = 1.
      end do


      return

    end subroutine CVFEM_ShapeFunctions


!!$ ===============
!!$ ===============

    subroutine shape_one_ele2( cv_ele_type, Mdims, GIdims, shape_fun )
!!$
!!$ This subrt defines the sub-control volume and FEM shape functions. Shape functions 
!!$    associated with volume integration using both CV basis functions CVN as well
!!$    as FEM basis functions CVFEN (and its derivatives CVFENLX, CVFENLY, CVFENLZ).
!!$    
!!$
      implicit none
      integer , intent(in) :: cv_ele_type
      type (multi_dimensions), intent(in) :: Mdims
      type (multi_gi_dimensions), intent(in) :: GIdims
      type(multi_shape_funs), intent(inout) :: shape_fun


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
      call determin_sloclist( shape_fun%cv_sloclist, Mdims%cv_nloc, Mdims%cv_snloc, GIdims%nface, &
           Mdims%ndim, cv_ele_type )

      if( Mdims%u_snloc == 1 )then
         shape_fun%u_sloclist( 1, 1 ) = 1 ; shape_fun%u_sloclist( 2, 1 ) = Mdims%u_nloc
      else
         call determin_sloclist( shape_fun%u_sloclist, Mdims%u_nloc, Mdims%u_snloc, GIdims%nface, &
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
    end subroutine shape_one_ele2


!!$ ===============
!!$ ===============

    subroutine shape_cv_n( cv_ele_type, Mdims, GIdims, shape_fun )
!!$
!!$ Shape functions associated with volume integration using both CV basis functions CVN 
!!$    as well as FEM basis functions N (and its derivatives NLX, NLY, NLZ) also for  
!!$    velocity basis functions UN, UNLX, UNLY, UNLZ.
!!$
      implicit none
      integer , intent(in) :: cv_ele_type
      type (multi_dimensions), intent(in) :: Mdims
      type (multi_gi_dimensions), intent(in) :: GIdims
      type(multi_shape_funs), intent(inout) :: shape_fun

!!$ New quadratic element quadrature by James and Zhi and Chris:

      Select Case( cv_ele_type )
      case( 1, 2 ) ! 1D
!!$      call quad_1d_shape( cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, n, nlx, un, unlx )
         call quad_1d_shape( GIdims%cv_ngi, Mdims%cv_nloc, Mdims%u_nloc, &
              shape_fun%cvn, shape_fun%cvweigh,                          &
              shape_fun%cvfen, shape_fun%cvfenlx_all(1, :, :),           &
              shape_fun%ufen, shape_fun%ufenlx_all(1, :, :) )

         shape_fun%cvfenlx_all(2, :, :) = 0. ; shape_fun%cvfenlx_all(3, :, :) = 0. ; &
              shape_fun%ufenlx_all(2, :, :) = 0. ; shape_fun%ufenlx_all(3, :, :) = 0. 

      case( 5, 6, 9, 10 ) ! Quadrilaterals and Hexahedra
!!$         call quad_nd_shape( ndim, cv_ele_type, cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, &
!!$              n, nlx, nly, nlz, &
!!$              un, unlx, unly, unlz )
         call quad_nd_shape( Mdims%ndim, cv_ele_type, GIdims%cv_ngi, Mdims%cv_nloc, Mdims%u_nloc, &
              shape_fun%cvn, shape_fun%cvweight,                                                  &
              shape_fun%cvfen, shape_fun%cvfenlx_all(1, :, :), shape_fun%cvfenlx_all(2, :, :),    &
              shape_fun%cvfenlx_all(3, :, :),                                                     &
              shape_fun%ufen, shape_fun%ufenlx_all(1, :, :), shape_fun%ufenlx_all(2, :, :),       &
              shape_fun%ufenlx_all(3, :, :) )

      case( 3, 4, 7, 8 ) ! Triangles and Tetrahedra
         if( new_quadratic_ele_quadrature .and. cv_ele_type==8) then
            call new_pt_qua_vol_cv_tri_tet_shape( cv_ele_type, Mdims%ndim, GIdims%cv_ngi,           &
                 Mdims%cv_nloc, Mdims%u_nloc, shape_fun%cvn,                                        &
                 shape_fun%cvweight, shape_fun%cvfen, shape_fun%cvfenlx_all(1, :, :),               &
                 shape_fun%cvfenlx_all(2, :, :), shape_fun%cvfenlx_all(3, :, :),                    &
                 shape_fun%ufen, shape_fun%ufenlx_all(1, :, :), shape_fun%ufenlx_all(2, :, :),      &
                 shape_fun%ufenlx_all(3, :, :) )

         else
            call vol_cv_tri_tet_shape( cv_ele_type, Mdims%ndim, GIdims%cv_ngi, Mdims%cv_nloc,       &
                 Mdims%u_nloc, shape_fun%cvn,                                                       &
                 shape_fun%cvweight, shape_fun%cvfen, shape_fun%cvfenlx_all(1, :, :),               &
                 shape_fun%cvfenlx_all(2, :, :), shape_fun%cvfenlx_all(3, :, :),                    &
                 shape_fun%ufen, shape_fun%ufenlx_all(1, :, :), shape_fun%ufenlx_all(2, :, :),      &
                 shape_fun%ufenlx_all(3, :, :) )     

         endif
         !stop 12
      case default; FLExit( "Wrong integer for CV_ELE_TYPE" )
      end Select

      return
    end subroutine shape_cv_n


!!$ ===============
!!$ ===============


    SUBROUTINE DET_SUF_ELE_SHAPE( SCVNGI, NFACE, &
         CVFEM_ON_FACE, &
         CV_NLOC, SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
         SCVFENLX, SCVFENLY, SCVFENLZ,  &
         U_NLOC,  SUFEN, SUFENSLX, SUFENSLY,  &
         SUFENLX, SUFENLY, SUFENLZ,  &
         SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
         SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
         SBUFEN, SBUFENSLX, SBUFENSLY, &
         SBUFENLX, SBUFENLY, SBUFENLZ, &
         CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC, &
         NDIM, CV_ELE_TYPE )
      !
      !     - this subroutine generates the FE basis functions, weights and the
      !     - derivatives of the shape functions for a variety of elements on the
      !     - control volume boundaries.
      !     - The routine also generates the shape functions and derivatives
      !     - associated with the CV surfaces and also the FV basis functions.
      !     -------------------------------
      !     - date last modified : 21/02/2012
      !     -------------------------------

      IMPLICIT NONE

      INTEGER, intent( in ) :: SCVNGI, CV_NLOC, U_NLOC, NFACE, &
           SBCVNGI, CV_SNLOC, U_SNLOC
      LOGICAL, DIMENSION( CV_NLOC, SCVNGI ), intent( in ) :: CVFEM_ON_FACE
      ! CV_ON_FACE(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
      REAL, DIMENSION( CV_NLOC, SCVNGI ), intent( in ) :: SCVFEN, SCVFENSLX, SCVFENSLY, &
           SCVFENLX, SCVFENLY, SCVFENLZ
      REAL, DIMENSION( SCVNGI ), intent( inout ) :: SCVFEWEIGH
      REAL, DIMENSION( U_NLOC, SCVNGI ), intent( in ) :: SUFEN, SUFENSLX, SUFENSLY, &
           SUFENLX, SUFENLY, SUFENLZ
      REAL, DIMENSION( CV_SNLOC, SBCVNGI ), intent( inout ) :: SBCVFEN, SBCVFENSLX, &
           SBCVFENSLY, SBCVFENLX, SBCVFENLY, SBCVFENLZ
      REAL, DIMENSION( SBCVNGI ), intent( inout ) :: SBCVFEWEIGH
      REAL, DIMENSION( U_SNLOC, SBCVNGI ), intent( inout ) :: SBUFEN, SBUFENSLX, SBUFENSLY, &
           SBUFENLX, SBUFENLY, SBUFENLZ
      INTEGER, DIMENSION( NFACE, CV_SNLOC ), intent( inout ) ::  CV_SLOCLIST
      INTEGER, DIMENSION( NFACE, U_SNLOC ), intent( inout ) ::  U_SLOCLIST
      INTEGER, intent( in ) :: NDIM, CV_ELE_TYPE
      ! Local variables
      INTEGER :: CV_KLOC, CV_SKLOC, U_KLOC, U_SKLOC, CV_BSNGI

      ewrite(3,*) 'In DET_SUF_ELE_SHAPE'

      ! Obtain SBCVFEN from SCVFEN:
      !ewrite(3,*)'for cv:'
      CALL SCVFEN_2_SBCVFEN( CV_NLOC, CV_SNLOC, SCVNGI, SBCVNGI, &
           CV_NLOC, CV_SNLOC, CVFEM_ON_FACE, &
           SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFENLX, SBCVFENLY, SBCVFENLZ, SBCVFEWEIGH, &
           SCVFEN, SCVFENSLX, SCVFENSLY, SCVFENLX, SCVFENLY, SCVFENLZ, SCVFEWEIGH )

      !ewrite(3,*)'U_NLOC, U_SNLOC, SCVNGI, SBCVNGI, CV_NLOC, CV_SNLOC:', &
      !         U_NLOC, U_SNLOC, SCVNGI, SBCVNGI, CV_NLOC, CV_SNLOC
      !ewrite(3,*)'for u:'
      ! Obtain SBUFEN from SUFEN:
      CALL SCVFEN_2_SBCVFEN( U_NLOC, U_SNLOC, SCVNGI, SBCVNGI, &
           CV_NLOC, CV_SNLOC, CVFEM_ON_FACE, &
           SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ, SBCVFEWEIGH, &
           SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY, SUFENLZ, SCVFEWEIGH )
      !ewrite(3,*)'SBUFEN:',SBUFEN
      !ewrite(3,*)'SUFEN:',SUFEN

      ! Determine CV_SLOCLIST & U_SLOCLIST
      CALL DETERMIN_SLOCLIST( CV_SLOCLIST, CV_NLOC, CV_SNLOC, NFACE, &
           NDIM, CV_ELE_TYPE )
      IF( U_SNLOC == 1 ) THEN
         U_SLOCLIST( 1, 1 ) = 1
         U_SLOCLIST( 2, 1 ) = U_NLOC
      ELSE
         CALL DETERMIN_SLOCLIST( U_SLOCLIST, U_NLOC, U_SNLOC, NFACE, &
              NDIM, CV_ELE_TYPE )
      ENDIF
      !ewrite(3,*)'CV_SNLOC, U_SNLOC, SCVNGI:', CV_SNLOC, U_SNLOC, SCVNGI
      !ewrite(3,*)'CV_SLOCLIST:', CV_SLOCLIST
      !ewrite(3,*)'U_SLOCLIST:', U_SLOCLIST

      RETURN
    END SUBROUTINE DET_SUF_ELE_SHAPE


!!$ ===============
!!$ ===============


    subroutine shapesv_fem_plus( scvngi, cv_neiloc, cv_on_face, cvfem_on_face, &
         ufem_on_face, &
         cv_ele_type, cv_nloc, scvfen, scvfenslx, scvfensly, scvfeweigh, &
         scvfenlx, scvfenly, scvfenlz, &
         u_nloc, sufen, sufenslx, sufensly, &
         sufenlx, sufenly, sufenlz, &
         ndim )
      implicit none
      !-
      !- This subroutine generates the FE basis functions, weights and the
      !- derivatives of the shape functions for a variety of elements on the
      !- control volume boundaries.
      !- The routine also generates the shape functions and derivatives
      !- associated with the CV surfaces and also the FV basis functions.
      !-
      !- date last modified : 29/11/2011
      !-
      integer, intent( in ) :: scvngi, cv_nloc
      integer, dimension( cv_nloc, scvngi ), intent( inout ) :: cv_neiloc
      logical, dimension( cv_nloc, scvngi ), intent( inout ) :: cv_on_face, cvfem_on_face
      logical, dimension( u_nloc, scvngi ), intent( inout ) :: ufem_on_face
      integer, intent( in ) :: cv_ele_type
      real, dimension( cv_nloc, scvngi ), intent( inout ) :: scvfen, scvfenslx, scvfensly
      real, dimension( scvngi ), intent( inout ) :: scvfeweigh
      real, dimension( cv_nloc, scvngi ), intent( inout ) :: scvfenlx, scvfenly, scvfenlz
      integer, intent( in ) :: u_nloc
      real, dimension( u_nloc, scvngi ), intent( inout ) :: sufen, sufenslx, sufensly, &
           sufenlx, sufenly, sufenlz
      integer, intent( in ) :: ndim
      ! Local variables
      integer :: iloc, gi
      logical :: tri_tet
      integer, dimension( :, : ), allocatable :: cvfem_neiloc, ufem_neiloc
      real, dimension( :, : ), allocatable :: m, mu, cvn_dummy
      real, dimension( : ), allocatable :: cvweigh_dummy

      !ewrite(3,*)' In ShapesV_Fem_Plus '

      ! Allocating space
      allocate( m( cv_nloc, scvngi ) )
      allocate( mu( cv_nloc, scvngi ) )
      allocate( cvn_dummy( cv_nloc, scvngi ) )
      allocate( cvweigh_dummy( scvngi ) )
      allocate( cvfem_neiloc( cv_nloc, scvngi ) )
      allocate( ufem_neiloc( u_nloc, scvngi ) )

      tri_tet=.false.

      Cond_ShapeType: Select Case( cv_ele_type )
      case( 1, 2 ) ! 1D
         call fv_1d_quad( scvngi, cv_nloc, scvfen, scvfenslx, scvfensly, scvfeweigh, &
              scvfenlx, scvfenly, scvfenlz ) ! For scalar fields
         call fv_1d_quad( scvngi, u_nloc, sufen, sufenslx, sufensly, scvfeweigh, &
              sufenlx, sufenly, sufenlz ) ! For U fields

      case( 3, 4, 7, 8 ) ! Triangle and Tetrahedra
         tri_tet=.true.
         call suf_cv_tri_tet_shape( cv_ele_type, ndim, scvngi, cv_nloc, u_nloc, scvfeweigh, &
              scvfen, scvfenlx, scvfenly, scvfenlz, scvfenslx, scvfensly,  &
              sufen, sufenlx, sufenly, sufenlz, sufenslx, sufensly, &
              cv_neiloc, cvfem_neiloc, ufem_neiloc )
      case( 5 ) ! Bi-linear Quadrilateral
         call fvquad( scvngi, cv_nloc, scvngi, &
              m, scvfen, scvfenslx, &
              scvfeweigh )
         call fvquad( scvngi, u_nloc, scvngi, &
              mu, sufen, sufenslx, &
              scvfeweigh )
         call quad_nd_shape( ndim, cv_ele_type, scvngi, cv_nloc, u_nloc, cvn_dummy, cvweigh_dummy, &
              scvfen, scvfenlx, scvfenly, scvfenlz, &
              sufen, sufenlx, sufenly, sufenlz )

      case( 6 ) ! Tri-linear Quadrilateral
         call fvqquad( scvngi, cv_nloc, scvngi, &
              m, scvfen, scvfenslx, &
              scvfeweigh )
         call fvqquad( scvngi, u_nloc, scvngi, &
              mu, sufen, sufenslx, &
              scvfeweigh )
         call quad_nd_shape( ndim, cv_ele_type, scvngi, cv_nloc, u_nloc, cvn_dummy, cvweigh_dummy, &
              scvfen, scvfenlx, scvfenly, scvfenlz, &
              sufen, sufenlx, sufenly, sufenlz )

      case( 9 ) ! Tri-linear Hexahedron
         call fvhex( scvngi, cv_nloc, scvngi, &
              m, scvfen, scvfenslx, &
              scvfensly, scvfeweigh )
         call fvhex( scvngi, u_nloc, scvngi, &
              mu, sufen, sufenslx, &
              sufensly, scvfeweigh )
         call quad_nd_shape( ndim, cv_ele_type, scvngi, cv_nloc, u_nloc, cvn_dummy, cvweigh_dummy, &
              scvfen, scvfenlx, scvfenly, scvfenlz, &
              sufen, sufenlx, sufenly, sufenlz )

      case( 10 ) ! Tri-linear Hexahedron
         call fvqhex( scvngi, cv_nloc, scvngi, &
              m, scvfen, scvfenslx, &
              scvfensly, scvfeweigh )
         call fvqhex( scvngi, u_nloc, scvngi, &
              mu, sufen, sufenslx, &
              sufensly, scvfeweigh )
         call quad_nd_shape( ndim, cv_ele_type, scvngi, cv_nloc, u_nloc, cvn_dummy, cvweigh_dummy, &
              scvfen, scvfenlx, scvfenly, scvfenlz, &
              sufen, sufenlx, sufenly, sufenlz )

      case default; FLExit( "Wrong integer for CV_ELE_TYPE" )

      end Select Cond_ShapeType


      if(.not.tri_tet) then
         call volnei( cv_neiloc, cvfem_neiloc, cv_nloc, scvngi, cv_ele_type )
      end if

      ewrite(3,*)'cv_ele_type:', cv_ele_type
      ewrite(3,*)'cv_nloc, scvngi:', cv_nloc, scvngi

      cv_on_face = .false. ; cvfem_on_face = .false.
      if ( ( cv_ele_type == 1 ) .or. ( cv_ele_type == 2 ) ) then ! 1D
         do iloc = 1, cv_nloc
            cv_on_face( iloc, iloc ) = .true.
            cv_on_face( iloc, iloc + 1 ) = .true.
         end do
         cv_on_face = cvfem_on_face
      else
         do iloc = 1, cv_nloc
            do gi = 1, scvngi 
               ! ewrite(3,*)'cv_neiloc, cvfem_on_face:', iloc, gi, &
               !      cv_neiloc( iloc, gi ), cvfem_neiloc( iloc, gi )
               if ( cv_neiloc( iloc, gi ) == -1 ) &
                    cv_on_face( iloc, gi ) = .true.
               if ( cvfem_neiloc( iloc, gi ) == -1 ) &
                    cvfem_on_face( iloc, gi ) = .true. 
            end do
         end do


         if(NEW_QUADRATIC_ELE_QUADRATURE.and.(cv_nloc==10).and.(ndim==3)) then
            ! Exterior faces :  1,3,6  ----James is this face the face with the 1st 6 surface quadrature points and the 1st 6 CV's.
            !                     This is only the exterior surface faces on the triangle with 1,2,3,4,5,6
            cvfem_on_face=.false.
            cvfem_on_face(1:6,1:6)=.true.
            ! Exterior faces :  1,3,10
            cvfem_on_face(1,7:12)=.true.
            cvfem_on_face(2,7:17)=.true.
            cvfem_on_face(3,7:12)=.true.
            cvfem_on_face(7,7:12)=.true.
            cvfem_on_face(8,7:12)=.true.
            cvfem_on_face(10,7:12)=.true.
            ! Exterior faces :  1,6,10
            cvfem_on_face(1,13:18)=.true.
            cvfem_on_face(4,13:18)=.true.
            cvfem_on_face(6,13:18)=.true.
            cvfem_on_face(7,13:18)=.true.
            cvfem_on_face(9,13:18)=.true.
            cvfem_on_face(10,13:18)=.true.
            ! Exterior faces :  3,6,10
            cvfem_on_face(3,19:24)=.true.
            cvfem_on_face(5,19:24)=.true.
            cvfem_on_face(6,19:24)=.true.
            cvfem_on_face(8,19:24)=.true.
            cvfem_on_face(9,19:24)=.true.
            cvfem_on_face(10,19:24)=.true.
         endif


      end if

      do iloc = 1, u_nloc
         do gi = 1, scvngi
            ufem_on_face( iloc, gi ) = ( ufem_neiloc( iloc, gi ) == -1 )
         end do
      end do

      !do iloc = 1, cv_nloc
      !   ewrite(3,*)'iloc, cv_on_face:', iloc, ( cv_on_face( iloc, gi ), gi = 1, scvngi )
      !   ewrite(3,*)'iloc, cvfem_neiloc:', iloc, ( cvfem_neiloc( iloc, gi ), gi = 1, scvngi )
      !   ewrite(3,*)'iloc, cvfem_on_face:', iloc, ( cvfem_on_face( iloc, gi ), gi = 1, scvngi )
      !end do

      deallocate( m )
      deallocate( mu )
      deallocate( cvn_dummy )
      deallocate( cvweigh_dummy )
      deallocate( cvfem_neiloc )

      return
    end subroutine shapesv_fem_plus


!!$ ===============
!!$ ===============


    SUBROUTINE GAUSSILOC( FINDGPTS, COLGPTS, NCOLGPTS,&
         NEILOC, NLOC, SVNGI )
      !     ----------------------------------------------------
      !
      ! This subroutine calculates FINDGPTS,COLGPTS,NCOLGPTS
      ! which contains given a local node ILOC the Gauss pts
      ! that are used to integrate around this local node.
      !
      !     -------------------------------
      !     - date last modified : 12/02/2002
      !     -------------------------------

      implicit none

      INTEGER, intent( in ) :: NLOC, SVNGI
      INTEGER, DIMENSION( NLOC + 1 ), intent( inout ) :: FINDGPTS
      ! We have overestimated the size of COLGPTS.
      INTEGER, DIMENSION( NLOC * SVNGI ), intent( inout ) :: COLGPTS
      INTEGER, DIMENSION( NLOC,  SVNGI ), intent( in ) :: NEILOC
      INTEGER, intent( inout ) :: NCOLGPTS
      ! Local
      INTEGER :: ILOC, COUNT, GI

      COUNT = 0

      DO ILOC = 1, NLOC

         FINDGPTS( ILOC ) = COUNT + 1

         DO GI = 1, SVNGI

            IF( NEILOC( ILOC, GI ) /= 0 ) THEN
               COUNT = COUNT + 1
               COLGPTS( COUNT ) = GI
            END IF
            !ewrite(3,*)'iloc,gi,NEILOC( ILOC, GI ):',iloc,gi,NEILOC( ILOC, GI )

         END DO

      END DO

      FINDGPTS( NLOC + 1 ) = COUNT + 1
      NCOLGPTS = COUNT
      !stop 2821

      RETURN
    END SUBROUTINE GAUSSILOC


!!$ ===============
!!$ ===============


  end module CVFEM_ShapeFunctions
