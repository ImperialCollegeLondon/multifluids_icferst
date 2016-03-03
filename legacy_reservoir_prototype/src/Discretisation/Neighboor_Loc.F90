
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


  module Neighboor_Loc
    use fldebug
    use multi_data_types


  contains

!!$ ===============
!!$ ===============

    SUBROUTINE VOLNEI( cv_ele_type, Mdims, GIdims, shape_fun, fem_neiloc )
      !NEILOC, FEM_NEILOC, NLOC, SVNGI, CV_ELE_TYPE )
!!$  This subroutine calculates NEILOC which is the array containing information 
!!$       given a local node and an integration point what is the other opposing
!!$       local node. It contains -1 if on the boundary of of the element.
      implicit none
      integer, intent( in ) :: cv_ele_type
      type (multi_dimensions), intent(in) :: Mdims
      type (multi_gi_dimensions), intent(in) :: GIdims
      type(multi_shape_funs), intent(inout) :: shape_fun
      integer, dimension( Mdims%cv_nloc, GIdims%scvngi ) :: fem_neiloc
      ! Local variables
      INTEGER :: ILOC, IEXT, IGP

      FEM_NEILOC = 0

      Conditional_Type: SELECT CASE( CV_ELE_TYPE )

      CASE( 1, 2 ) ! 1D
         DO ILOC = 1, NLOC
            SHAPE_FUN%CV_NEILOC( ILOC, ILOC ) = ILOC - 1
            SHAPE_FUN%CV_NEILOC( ILOC, ILOC + 1 ) = ILOC + 1
         END DO
         SHAPE_FUN%CV_NEILOC( 1, 1 ) = -1
         SHAPE_FUN%CV_NEILOC( NLOC, SVNGI ) = -1
         FEM_NEILOC = SHAPE_FUN%CV_NEILOC

      CASE( 3 ) ! Linear Triangle
         FLAbort( " Defined elsewhere -- it needs to be updated " )

         ILOC = 1 ! CV 1
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 3 ) = 3
         ! External
         shape_fun%cv_neiloc( iloc, 4 ) = -1
         shape_fun%cv_neiloc( iloc, 9 ) = -1

         ILOC = 2 ! CV 2
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 2 ) = 3
         ! External
         shape_fun%cv_neiloc( iloc, 5 ) = -1
         shape_fun%cv_neiloc( iloc, 6 ) = -1

         ILOC = 3 ! CV 3
         SHAPE_FUN%CV_NEILOC( ILOC, 2) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 3) = 1
         ! External
         shape_fun%cv_neiloc( iloc, 7 ) = -1
         shape_fun%cv_neiloc( iloc, 8 ) = -1

         fem_neiloc = shape_fun%cv_neiloc

      CASE( 4 ) ! Quadratic Triangle
         FLAbort( " Defined elsewhere -- it needs to be updated " )
         ILOC = 1
         ! Face 1
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 2 ) = 2
         ! Face 2
         SHAPE_FUN%CV_NEILOC( ILOC,3 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC,4 ) = 6
         ! External
         shape_fun%cv_neiloc( iloc, 19 ) = -1
         shape_fun%cv_neiloc( iloc, 20 ) = -1
         shape_fun%cv_neiloc( iloc, 35 ) = -1
         shape_fun%cv_neiloc( iloc, 36 ) = -1

         ILOC = 2
         ! Face 1
         SHAPE_FUN%CV_NEILOC( ILOC,1 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC,2 ) = 1
         ! Face 7
         SHAPE_FUN%CV_NEILOC( ILOC,13 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC,14 ) = 6
         ! Face 8
         SHAPE_FUN%CV_NEILOC( ILOC,15 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC,16 ) = 4
         ! Face 3
         SHAPE_FUN%CV_NEILOC( ILOC,5 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC,6 ) = 3
         ! External
         shape_fun%cv_neiloc( iloc, 21 ) = -1
         shape_fun%cv_neiloc( iloc, 22 ) = -1

         ILOC = 3
         ! Face 3
         SHAPE_FUN%CV_NEILOC( ILOC,5 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC,6 ) = 2
         ! Face 4
         SHAPE_FUN%CV_NEILOC( ILOC,7 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC,8 ) = 4
         ! External
         shape_fun%cv_neiloc( iloc, 23 ) = -1
         shape_fun%cv_neiloc( iloc, 24 ) = -1
         shape_fun%cv_neiloc( iloc, 25 ) = -1
         shape_fun%cv_neiloc( iloc, 26 ) = -1

         ILOC = 4
         ! Face 4
         SHAPE_FUN%CV_NEILOC( ILOC,7 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC,8 ) = 3
         ! Face 8
         SHAPE_FUN%CV_NEILOC( ILOC,15 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC,16 ) = 2
         ! Face 9
         SHAPE_FUN%CV_NEILOC( ILOC,17 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC,18 ) = 6
         ! Face 5
         SHAPE_FUN%CV_NEILOC( ILOC,9 )  = 5
         SHAPE_FUN%CV_NEILOC( ILOC,10 ) = 5
         ! External
         shape_fun%cv_neiloc( iloc, 27 ) = -1
         shape_fun%cv_neiloc( iloc, 28 ) = -1

         ILOC = 5
         ! Face 5
         SHAPE_FUN%CV_NEILOC( ILOC,9 )  = 4
         SHAPE_FUN%CV_NEILOC( ILOC,10 ) = 4
         ! Face 6
         SHAPE_FUN%CV_NEILOC( ILOC,11 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC,12 ) = 6
         ! External
         shape_fun%cv_neiloc( iloc, 29 ) = -1
         shape_fun%cv_neiloc( iloc, 30 ) = -1
         shape_fun%cv_neiloc( iloc, 31 ) = -1
         shape_fun%cv_neiloc( iloc, 32 ) = -1

         ILOC = 6
         ! Face 6
         SHAPE_FUN%CV_NEILOC( ILOC,11 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC,12 ) = 5
         ! Face 9
         SHAPE_FUN%CV_NEILOC( ILOC,17 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC,18 ) = 4
         ! Face 7
         SHAPE_FUN%CV_NEILOC( ILOC,13 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC,14 ) = 2
         ! Face 2
         SHAPE_FUN%CV_NEILOC( ILOC,3 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC,4 ) = 1
         ! External
         shape_fun%cv_neiloc( iloc, 35 ) = -1
         shape_fun%cv_neiloc( iloc, 36 ) = -1

      CASE( 5  ) ! Bi-linear Quadrilateral
         ILOC = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 3 ) = 2
         ! External
         shape_fun%cv_neiloc( iloc, 5 ) = -1
         shape_fun%cv_neiloc( iloc, 6 ) = -1

         ILOC = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 2 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 3 ) = 1
         ! External
         shape_fun%cv_neiloc( iloc, 7 ) = -1
         shape_fun%cv_neiloc( iloc, 8 ) = -1

         ILOC = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 4 ) = 4
         ! External
         shape_fun%cv_neiloc( iloc, 11 ) = -1
         shape_fun%cv_neiloc( iloc, 12 ) = -1

         ILOC = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 2 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 4 ) = 3
         ! External
         shape_fun%cv_neiloc( iloc, 9 ) = -1
         shape_fun%cv_neiloc( iloc, 10 ) = -1

      CASE( 6 ) ! Bi-quadratic Quadrilateral
         ILOC = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 2 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 5 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 6 ) = 4
         ! External
         shape_fun%cv_neiloc( iloc, 25 ) = -1
         shape_fun%cv_neiloc( iloc, 26 ) = -1
         shape_fun%cv_neiloc( iloc, 27 ) = -1
         shape_fun%cv_neiloc( iloc, 28 ) = -1

         ILOC = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 2 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 7 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 8 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 3 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 4 ) = 3
         ! External
         shape_fun%cv_neiloc( iloc, 47 ) = -1
         shape_fun%cv_neiloc( iloc, 48 ) = -1

         ILOC = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 3) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 4) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 9) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 10 ) = 6
         ! External
         shape_fun%cv_neiloc( iloc, 43 ) = -1
         shape_fun%cv_neiloc( iloc, 44 ) = -1
         shape_fun%cv_neiloc( iloc, 45 ) = -1
         shape_fun%cv_neiloc( iloc, 46 ) = -1

         ILOC = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 5) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 6) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 11 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 12 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 15 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 16 ) = 7
         ! External
         shape_fun%cv_neiloc( iloc, 29 ) = -1
         shape_fun%cv_neiloc( iloc, 30 ) = -1

         ILOC = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 7) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 8) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 13 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 14 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 18 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 17 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 12 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 11 ) = 4

         ILOC = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 9) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 10 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 13 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 14 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 19 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 20 ) = 9
         ! External
         shape_fun%cv_neiloc( iloc, 41 ) = -1
         shape_fun%cv_neiloc( iloc, 42 ) = -1

         ILOC = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 15 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 16 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 21 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 22 ) = 8
         ! External
         shape_fun%cv_neiloc( iloc, 31 ) = -1
         shape_fun%cv_neiloc( iloc, 32 ) = -1
         shape_fun%cv_neiloc( iloc, 33 ) = -1
         shape_fun%cv_neiloc( iloc, 34 ) = -1

         ILOC = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 21 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 22 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 17 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 18 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 23 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 24 ) = 9
         ! External
         shape_fun%cv_neiloc( iloc, 35 ) = -1
         shape_fun%cv_neiloc( iloc, 36 ) = -1

         ILOC = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 23 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 24 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 19 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 20 ) = 6
         ! External
         shape_fun%cv_neiloc( iloc, 37 ) = -1
         shape_fun%cv_neiloc( iloc, 38 ) = -1
         shape_fun%cv_neiloc( iloc, 39 ) = -1
         shape_fun%cv_neiloc( iloc, 40 ) = -1

      CASE( 7 ) ! Linear Tetrahedron
         iext = 6
         igp = 1
         ILOC = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 3 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 4 ) = 4
         ! External
         shape_fun%cv_neiloc( iloc, iext + 1 : iext + 3 * igp ) = -1

         ILOC = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 2 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 5 ) = 4
         ! External
         shape_fun%cv_neiloc( iloc, iext + 3 * igp + 1 : iext + 6 * igp ) = -1

         ILOC = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 2 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 3 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 6 ) = 4
         ! External
         shape_fun%cv_neiloc( iloc, iext + 6 * igp + 1 : iext + 9 * igp ) = -1

         ILOC = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 4 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 5 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 6 ) = 3
         ! External
         shape_fun%cv_neiloc( iloc, iext + 9 * igp + 1 : iext + 12 * igp ) = -1

      CASE( 8 ) ! Quadratic Tetrahedron
         iext = 96
         igp = 4
         ILOC = 1
         ! Face 1
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 2 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 3 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 4 ) = 6
         ! Face 2
         SHAPE_FUN%CV_NEILOC( ILOC, 5 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 6 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 7 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 8 ) = 2
         ! Face 3
         SHAPE_FUN%CV_NEILOC( ILOC, 9 )  = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 10 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 11 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 12 ) = 7
         ! External
         shape_fun%cv_neiloc( iloc, iext + 1 : iext + 3 * igp ) = -1

         ILOC = 2
         ! Face 2
         SHAPE_FUN%CV_NEILOC( ILOC, 5 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 6 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 7 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 8 ) = 1
         ! Face 5
         SHAPE_FUN%CV_NEILOC( ILOC, 17 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 18 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 19 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 20 ) = 3
         ! Face 20
         SHAPE_FUN%CV_NEILOC( ILOC, 77 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 78 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 79 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 80 ) = 8
         ! Face 21
         SHAPE_FUN%CV_NEILOC( ILOC, 81 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 82 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 83 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 84 ) = 7
         ! Face 23
         SHAPE_FUN%CV_NEILOC( ILOC, 89 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 80 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 91 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 92 ) = 6
         ! Face 24
         SHAPE_FUN%CV_NEILOC( ILOC, 93 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 94 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 95 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 96 ) = 4
         ! External
         shape_fun%cv_neiloc( iloc, iext + 3 * igp + 1 : iext + 5 * igp ) = -1

         ILOC = 3
         ! Face 4
         SHAPE_FUN%CV_NEILOC( ILOC, 13 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 14 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 15 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 16 ) = 4
         ! Face 5
         SHAPE_FUN%CV_NEILOC( ILOC, 17 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 18 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 19 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 20 ) = 2
         ! Face 6
         SHAPE_FUN%CV_NEILOC( ILOC, 21 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 22 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 23 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 24 ) = 8
         ! External
         shape_fun%cv_neiloc( iloc, iext + 5 * igp + 1 : iext + 8 * igp ) = -1

         ILOC = 4
         ! Face 4
         SHAPE_FUN%CV_NEILOC( ILOC, 13 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 14 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 15 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 16 ) = 3
         ! Face 8
         SHAPE_FUN%CV_NEILOC( ILOC, 29 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 30 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 31 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 32 ) = 5
         ! Face 18
         SHAPE_FUN%CV_NEILOC( ILOC, 69 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 70 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 71 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 72 ) = 9
         ! Face 19
         SHAPE_FUN%CV_NEILOC( ILOC, 73 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 74 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 75 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 76 ) = 8
         ! Face 22
         SHAPE_FUN%CV_NEILOC( ILOC, 85 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 86 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 87 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 88 ) = 6
         ! Face 24
         SHAPE_FUN%CV_NEILOC( ILOC, 93 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 94 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 95 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 96 ) = 2
         ! External
         shape_fun%cv_neiloc( iloc, iext + 8 * igp + 1 : iext + 10 * igp ) = -1

         ILOC = 5
         ! Face 7
         SHAPE_FUN%CV_NEILOC( ILOC, 25 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 26 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 27 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 28 ) = 6
         ! Face 8
         SHAPE_FUN%CV_NEILOC( ILOC, 29 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 30 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 31 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 32 ) = 4
         ! Face 9
         SHAPE_FUN%CV_NEILOC( ILOC, 33 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 34 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 35 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 36 ) = 9
         ! External
         shape_fun%cv_neiloc( iloc, iext + 10 * igp + 1 : iext + 13 * igp ) = -1

         ILOC = 6
         ! Face 1
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 2 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 3 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 4 ) = 1
         ! Face 7
         SHAPE_FUN%CV_NEILOC( ILOC, 25 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 26 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 27 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 28 ) = 5
         ! Face 16
         SHAPE_FUN%CV_NEILOC( ILOC, 61 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 62 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 63 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 64 ) = 7
         ! Face 17
         SHAPE_FUN%CV_NEILOC( ILOC, 65 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 66 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 67 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 68 ) = 9
         ! Face 22
         SHAPE_FUN%CV_NEILOC( ILOC, 85 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 86 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 87 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 88 ) = 4
         ! Face 23
         SHAPE_FUN%CV_NEILOC( ILOC, 89 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 90 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 91 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 92 ) = 2
         ! External
         shape_fun%cv_neiloc( iloc, iext + 13 * igp + 1 : iext + 15 * igp ) = -1
         !
         ILOC = 7
         ! Face 3
         SHAPE_FUN%CV_NEILOC( ILOC, 9 )  = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 10 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 11 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 12 ) = 1
         ! Face 11
         SHAPE_FUN%CV_NEILOC( ILOC, 41 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 42 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 43 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 44 ) = 10
         ! Face 13
         SHAPE_FUN%CV_NEILOC( ILOC, 49 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 50 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 51 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 52 ) = 9
         ! Face 15
         SHAPE_FUN%CV_NEILOC( ILOC, 57 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 58 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 59 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 60 ) = 8
         ! Face 16
         SHAPE_FUN%CV_NEILOC( ILOC, 61 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 62 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 63 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 64 ) = 6
         ! Face 21
         SHAPE_FUN%CV_NEILOC( ILOC, 81 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 82 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 83 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 84 ) = 2
         ! External
         shape_fun%cv_neiloc( iloc, iext + 15 * igp + 1 : iext + 17 * igp ) = -1

         ILOC = 8
         ! Face 6
         SHAPE_FUN%CV_NEILOC( ILOC, 21 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 22 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 23 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 24 ) = 3
         ! Face 12
         SHAPE_FUN%CV_NEILOC( ILOC, 45 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 46 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 47 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 48 ) = 10
         ! Face 14
         SHAPE_FUN%CV_NEILOC( ILOC, 53 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 54 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 55 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 56 ) = 9
         ! Face 15
         SHAPE_FUN%CV_NEILOC( ILOC, 57 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 58 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 59 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 60 ) = 7
         ! Face 19
         SHAPE_FUN%CV_NEILOC( ILOC, 73 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 74 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 75 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 76 ) = 4
         ! Face 20
         SHAPE_FUN%CV_NEILOC( ILOC, 77 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 78 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 79 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 80 ) = 2
         ! External
         shape_fun%cv_neiloc( iloc, iext + 17 * igp + 1 : iext + 20 * igp ) = -1

         ILOC = 9
         ! Face 9
         SHAPE_FUN%CV_NEILOC( ILOC, 33 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 34 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 35 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 36 ) = 5
         ! Face 10
         SHAPE_FUN%CV_NEILOC( ILOC, 37 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 38 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 39 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 40 ) = 10
         ! Face 13
         SHAPE_FUN%CV_NEILOC( ILOC, 49 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 50 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 51 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 52 ) = 7
         ! Face 14
         SHAPE_FUN%CV_NEILOC( ILOC, 53 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 54 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 55 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 56 ) = 8
         ! Face 17
         SHAPE_FUN%CV_NEILOC( ILOC, 65 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 66 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 67 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 68 ) = 6
         ! Face 18
         SHAPE_FUN%CV_NEILOC( ILOC, 69 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 70 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 71 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 72 ) = 4
         ! External
         shape_fun%cv_neiloc( iloc, iext + 20 * igp + 1 : iext + 22 * igp ) = -1

         ILOC = 10
         ! Face 10
         SHAPE_FUN%CV_NEILOC( ILOC, 37 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 38 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 39 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 40 ) = 9
         ! Face 11
         SHAPE_FUN%CV_NEILOC( ILOC, 41 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 42 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 43 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 44 ) = 7
         ! Face 12
         SHAPE_FUN%CV_NEILOC( ILOC, 45 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 46 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 47 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 48 ) = 8
         ! External
         shape_fun%cv_neiloc( iloc, iext + 22 * igp + 1 : iext + 24 * igp ) = -1

      CASE( 9 ) ! Tri-linear Hexahedron
         ILOC = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 3 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 5 ) = 5
         ! External
         shape_fun%cv_neiloc( iloc, 13 ) = -1
         shape_fun%cv_neiloc( iloc, 17 ) = -1
         shape_fun%cv_neiloc( iloc, 21 ) = -1

         ILOC = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 2 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 3 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 6 ) = 6
         ! External
         shape_fun%cv_neiloc( iloc, 14 ) = -1
         shape_fun%cv_neiloc( iloc, 18 ) = -1
         shape_fun%cv_neiloc( iloc, 25 ) = -1

         ILOC = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 4 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 7 ) = 7
         ! External
         shape_fun%cv_neiloc( iloc, 16 ) = -1
         shape_fun%cv_neiloc( iloc, 22 ) = -1
         shape_fun%cv_neiloc( iloc, 31 ) = -1

         ILOC = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 2 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 4 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 8 ) = 8
         ! External
         shape_fun%cv_neiloc( iloc, 15 ) = -1
         shape_fun%cv_neiloc( iloc, 26 ) = -1
         shape_fun%cv_neiloc( iloc, 30 ) = -1

         ILOC = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 5 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 9 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 11 ) = 6
         ! External
         shape_fun%cv_neiloc( iloc, 19 ) = -1
         shape_fun%cv_neiloc( iloc, 24 ) = -1
         shape_fun%cv_neiloc( iloc, 33 ) = -1

         ILOC = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 6 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 10 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 11 ) = 5
         ! External
         shape_fun%cv_neiloc( iloc, 20 ) = -1
         shape_fun%cv_neiloc( iloc, 27 ) = -1
         shape_fun%cv_neiloc( iloc, 34 ) = -1

         ILOC = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 7 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 9 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 12 ) = 8
         ! External
         shape_fun%cv_neiloc( iloc, 23 ) = -1
         shape_fun%cv_neiloc( iloc, 32 ) = -1
         shape_fun%cv_neiloc( iloc, 36 ) = -1

         ILOC = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 8 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 10 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 12 ) = 7
         ! External
         shape_fun%cv_neiloc( iloc, 28 ) = -1
         shape_fun%cv_neiloc( iloc, 29 ) = -1
         shape_fun%cv_neiloc( iloc, 35 ) = -1

      CASE( 10 ) ! Tri-quadratic Hexahedron
         iext = 216
         igp = 4
         ILOC = 1
         ! Face 1
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 2 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 3 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 4 ) = 2
         ! Face 3
         SHAPE_FUN%CV_NEILOC( ILOC, 9 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 10 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 11 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 12 ) = 4
         ! Face 13
         SHAPE_FUN%CV_NEILOC( ILOC, 49 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 50 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 51 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 52 ) = 10
         ! External
         shape_fun%cv_neiloc( iloc, iext + 1 : iext + 3 * igp ) = -1

         ILOC = 2
         ! Face 1
         SHAPE_FUN%CV_NEILOC( ILOC, 1 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 2 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 3 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 4 ) = 1
         ! Face 4
         SHAPE_FUN%CV_NEILOC( ILOC, 13 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 14 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 15 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 16 ) = 5
         ! Face 2
         SHAPE_FUN%CV_NEILOC( ILOC, 5 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 6 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 7 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 8 ) = 3
         ! Face 14
         SHAPE_FUN%CV_NEILOC( ILOC, 53 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 54 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 55 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 56 ) = 11
         ! External
         shape_fun%cv_neiloc( iloc, iext + 3 * igp + 1 : iext + 5 * igp ) = -1

         ILOC = 3
         ! Face 2
         SHAPE_FUN%CV_NEILOC( ILOC, 5 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 6 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 7 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 8 ) = 2
         ! Face 5
         SHAPE_FUN%CV_NEILOC( ILOC, 17 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 18 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 19 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 20 ) = 6
         ! Face 15
         SHAPE_FUN%CV_NEILOC( ILOC, 57 ) = 12
         SHAPE_FUN%CV_NEILOC( ILOC, 58 ) = 12
         SHAPE_FUN%CV_NEILOC( ILOC, 59 ) = 12
         SHAPE_FUN%CV_NEILOC( ILOC, 60 ) = 12
         ! External
         shape_fun%cv_neiloc( iloc, iext + 5 * igp + 1 : iext + 8 * igp ) = -1

         ILOC = 4
         ! Face 3
         SHAPE_FUN%CV_NEILOC( ILOC, 9 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 10 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 11 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 12 ) = 1
         ! Face 6
         SHAPE_FUN%CV_NEILOC( ILOC, 21 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 22 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 23 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 24 ) = 5
         ! Face 8
         SHAPE_FUN%CV_NEILOC( ILOC, 29 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 30 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 31 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 32 ) = 7
         ! Face 16
         SHAPE_FUN%CV_NEILOC( ILOC, 61 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 62 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 63 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 64 ) = 13
         ! External
         shape_fun%cv_neiloc( iloc, iext + 8 * igp + 1 : iext + 10 * igp ) = -1

         ILOC = 5
         ! Face 4
         SHAPE_FUN%CV_NEILOC( ILOC, 13 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 14 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 15 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 16 ) = 2
         ! Face 7
         SHAPE_FUN%CV_NEILOC( ILOC, 25 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 26 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 27 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 28 ) = 6
         ! Face 9
         SHAPE_FUN%CV_NEILOC( ILOC, 33 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 34 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 35 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 36 ) = 8
         ! Face 6
         SHAPE_FUN%CV_NEILOC( ILOC, 21 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 22 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 23 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 24 ) = 4
         ! Face 17
         SHAPE_FUN%CV_NEILOC( ILOC, 65 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 66 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 67 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 68 ) = 14
         ! External
         shape_fun%cv_neiloc( iloc, iext + 10 * igp + 1 : iext + 11 * igp ) = -1

         ILOC = 6
         ! Face 5
         SHAPE_FUN%CV_NEILOC( ILOC, 17 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 18 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 19 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 20 ) = 3
         ! Face 7
         SHAPE_FUN%CV_NEILOC( ILOC, 25 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 26 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 27 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 28 ) = 5
         ! Face 10
         SHAPE_FUN%CV_NEILOC( ILOC, 37 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 38 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 39 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 40 ) = 9
         ! Face 18
         SHAPE_FUN%CV_NEILOC( ILOC, 69 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 70 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 71 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 72 ) = 15
         ! External
         shape_fun%cv_neiloc( iloc, iext + 11 * igp + 1 : iext + 13 * igp ) = -1

         ILOC = 7
         ! Face 8
         SHAPE_FUN%CV_NEILOC( ILOC, 29 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 30 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 31 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 32 ) = 4
         ! Face 11
         SHAPE_FUN%CV_NEILOC( ILOC, 41 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 42 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 43 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 44 ) = 8
         ! Face 19
         SHAPE_FUN%CV_NEILOC( ILOC, 73 ) = 16
         SHAPE_FUN%CV_NEILOC( ILOC, 74 ) = 16
         SHAPE_FUN%CV_NEILOC( ILOC, 75 ) = 16
         SHAPE_FUN%CV_NEILOC( ILOC, 76 ) = 16
         ! External
         shape_fun%cv_neiloc( iloc, iext + 13 * igp + 1 : iext + 16 * igp ) = -1

         ILOC = 8
         ! Face 11
         SHAPE_FUN%CV_NEILOC( ILOC, 41 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 42 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 43 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 44 ) = 7
         ! Face 9
         SHAPE_FUN%CV_NEILOC( ILOC, 33 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 34 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 35 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 36 ) = 5
         ! Face 12
         SHAPE_FUN%CV_NEILOC( ILOC, 45 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 46 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 47 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 48 ) = 9
         ! Face 20
         SHAPE_FUN%CV_NEILOC( ILOC, 77 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 78 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 79 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 80 ) = 17
         ! External
         shape_fun%cv_neiloc( iloc, iext + 16 * igp + 1 : iext + 18 * igp ) = -1

         ILOC = 9
         ! Face 10
         SHAPE_FUN%CV_NEILOC( ILOC, 37 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 38 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 39 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 40 ) = 6
         ! Face 12
         SHAPE_FUN%CV_NEILOC( ILOC, 45 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 46 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 47 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 48 ) = 8
         ! Face 21
         SHAPE_FUN%CV_NEILOC( ILOC, 81 ) = 18
         SHAPE_FUN%CV_NEILOC( ILOC, 82 ) = 18
         SHAPE_FUN%CV_NEILOC( ILOC, 83 ) = 18
         SHAPE_FUN%CV_NEILOC( ILOC, 84 ) = 18
         ! External
         shape_fun%cv_neiloc( iloc, iext + 18 * igp + 1 : iext + 21 * igp ) = -1

         ILOC = 10
         ! Face 13
         SHAPE_FUN%CV_NEILOC( ILOC, 49 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 50 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 51 ) = 1
         SHAPE_FUN%CV_NEILOC( ILOC, 52 ) = 1
         ! Face 22
         SHAPE_FUN%CV_NEILOC( ILOC, 85 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 86 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 87 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 88 ) = 11
         ! Face 24
         SHAPE_FUN%CV_NEILOC( ILOC, 93 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 94 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 95 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 96 ) = 13
         ! Face 34
         SHAPE_FUN%CV_NEILOC( ILOC, 133 ) = 19
         SHAPE_FUN%CV_NEILOC( ILOC, 134 ) = 19
         SHAPE_FUN%CV_NEILOC( ILOC, 135 ) = 19
         SHAPE_FUN%CV_NEILOC( ILOC, 136 ) = 19
         ! External
         shape_fun%cv_neiloc( iloc, iext + 21 * igp + 1 : iext + 23 * igp ) = -1

         ILOC = 11
         ! Face 14
         SHAPE_FUN%CV_NEILOC( ILOC, 53 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 54 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 55 ) = 2
         SHAPE_FUN%CV_NEILOC( ILOC, 56 ) = 2
         ! Face 22
         SHAPE_FUN%CV_NEILOC( ILOC, 85 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 86 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 87 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 88 ) = 10
         ! Face 25
         SHAPE_FUN%CV_NEILOC( ILOC, 97 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 98 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 99 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 100 ) = 14
         ! Face 23
         SHAPE_FUN%CV_NEILOC( ILOC, 89 ) = 12
         SHAPE_FUN%CV_NEILOC( ILOC, 90 ) = 12
         SHAPE_FUN%CV_NEILOC( ILOC, 91 ) = 12
         SHAPE_FUN%CV_NEILOC( ILOC, 92 ) = 12
         ! Face 35
         SHAPE_FUN%CV_NEILOC( ILOC, 137 ) = 20
         SHAPE_FUN%CV_NEILOC( ILOC, 138 ) = 20
         SHAPE_FUN%CV_NEILOC( ILOC, 139 ) = 20
         SHAPE_FUN%CV_NEILOC( ILOC, 140 ) = 20
         ! External
         shape_fun%cv_neiloc( iloc, iext + 23 * igp + 1 : iext + 24 * igp ) = -1

         ILOC = 12
         ! Face 15
         SHAPE_FUN%CV_NEILOC( ILOC, 57 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 58 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 59 ) = 3
         SHAPE_FUN%CV_NEILOC( ILOC, 60 ) = 3
         ! Face 23
         SHAPE_FUN%CV_NEILOC( ILOC, 89 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 90 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 91 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 92 ) = 11
         ! Face 26
         SHAPE_FUN%CV_NEILOC( ILOC, 101 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 102 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 103 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 104 ) = 15
         ! Face 36
         SHAPE_FUN%CV_NEILOC( ILOC, 141 ) = 21
         SHAPE_FUN%CV_NEILOC( ILOC, 142 ) = 21
         SHAPE_FUN%CV_NEILOC( ILOC, 143 ) = 21
         SHAPE_FUN%CV_NEILOC( ILOC, 144 ) = 21
         ! External
         shape_fun%cv_neiloc( iloc, iext + 24 * igp + 1 : iext + 26 * igp ) = -1

         ILOC = 13
         ! Face 16
         SHAPE_FUN%CV_NEILOC( ILOC, 61 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 62 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 63 ) = 4
         SHAPE_FUN%CV_NEILOC( ILOC, 64 ) = 4
         ! Face 24
         SHAPE_FUN%CV_NEILOC( ILOC, 93 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 94 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 95 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 96 ) = 10
         ! Face 27
         SHAPE_FUN%CV_NEILOC( ILOC, 105 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 106 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 107 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 108 ) = 14
         ! Face 29
         SHAPE_FUN%CV_NEILOC( ILOC, 113 ) = 16
         SHAPE_FUN%CV_NEILOC( ILOC, 114 ) = 16
         SHAPE_FUN%CV_NEILOC( ILOC, 115 ) = 16
         SHAPE_FUN%CV_NEILOC( ILOC, 116 ) = 16
         ! Face 37
         SHAPE_FUN%CV_NEILOC( ILOC, 145 ) = 22
         SHAPE_FUN%CV_NEILOC( ILOC, 146 ) = 22
         SHAPE_FUN%CV_NEILOC( ILOC, 147 ) = 22
         SHAPE_FUN%CV_NEILOC( ILOC, 148 ) = 22
         ! External
         shape_fun%cv_neiloc( iloc, iext + 26 * igp + 1 : iext + 27 * igp ) = -1

         ILOC = 14
         ! Face 17
         SHAPE_FUN%CV_NEILOC( ILOC, 65 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 66 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 67 ) = 5
         SHAPE_FUN%CV_NEILOC( ILOC, 68 ) = 5
         ! Face 25
         SHAPE_FUN%CV_NEILOC( ILOC, 97 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 98 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 99 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 100 ) = 11
         ! Face 28
         SHAPE_FUN%CV_NEILOC( ILOC, 109 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 110 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 111 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 112 ) = 15
         ! Face 30
         SHAPE_FUN%CV_NEILOC( ILOC, 117 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 118 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 119 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 120 ) = 17
         ! Face 27
         SHAPE_FUN%CV_NEILOC( ILOC, 105 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 106 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 107 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 108 ) = 13
         ! Face 38
         SHAPE_FUN%CV_NEILOC( ILOC, 149 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 150 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 151 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 152 ) = 23

         ILOC = 15
         ! Face 18
         SHAPE_FUN%CV_NEILOC( ILOC, 69 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 70 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 71 ) = 6
         SHAPE_FUN%CV_NEILOC( ILOC, 72 ) = 6
         ! Face 26
         SHAPE_FUN%CV_NEILOC( ILOC, 101 ) = 12
         SHAPE_FUN%CV_NEILOC( ILOC, 102 ) = 12
         SHAPE_FUN%CV_NEILOC( ILOC, 103 ) = 12
         SHAPE_FUN%CV_NEILOC( ILOC, 104 ) = 12
         ! Face 28
         SHAPE_FUN%CV_NEILOC( ILOC, 109 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 110 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 111 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 112 ) = 14
         ! Face 31
         SHAPE_FUN%CV_NEILOC( ILOC, 121 ) = 18
         SHAPE_FUN%CV_NEILOC( ILOC, 122 ) = 18
         SHAPE_FUN%CV_NEILOC( ILOC, 123 ) = 18
         SHAPE_FUN%CV_NEILOC( ILOC, 124 ) = 18
         ! Face 39
         SHAPE_FUN%CV_NEILOC( ILOC, 153 ) = 24
         SHAPE_FUN%CV_NEILOC( ILOC, 154 ) = 24
         SHAPE_FUN%CV_NEILOC( ILOC, 155 ) = 24
         SHAPE_FUN%CV_NEILOC( ILOC, 156 ) = 24
         ! External
         shape_fun%cv_neiloc( iloc, iext + 27 * igp + 1 : iext + 28 * igp ) = -1

         ILOC = 16
         ! Face 19
         SHAPE_FUN%CV_NEILOC( ILOC, 73 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 74 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 75 ) = 7
         SHAPE_FUN%CV_NEILOC( ILOC, 76 ) = 7
         ! Face 29
         SHAPE_FUN%CV_NEILOC( ILOC, 113 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 114 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 115 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 116 ) = 13
         ! Face 32
         SHAPE_FUN%CV_NEILOC( ILOC, 125 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 126 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 127 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 128 ) = 17
         ! Face 40
         SHAPE_FUN%CV_NEILOC( ILOC, 157 ) = 25
         SHAPE_FUN%CV_NEILOC( ILOC, 158 ) = 25
         SHAPE_FUN%CV_NEILOC( ILOC, 159 ) = 25
         SHAPE_FUN%CV_NEILOC( ILOC, 160 ) = 25
         ! External
         shape_fun%cv_neiloc( iloc, iext + 28 * igp + 1 : iext + 30 * igp ) = -1

         ILOC = 17
         ! Face 20
         SHAPE_FUN%CV_NEILOC( ILOC, 77 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 78 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 79 ) = 8
         SHAPE_FUN%CV_NEILOC( ILOC, 80 ) = 8
         ! Face 32
         SHAPE_FUN%CV_NEILOC( ILOC, 125 ) = 16
         SHAPE_FUN%CV_NEILOC( ILOC, 126 ) = 16
         SHAPE_FUN%CV_NEILOC( ILOC, 127 ) = 16
         SHAPE_FUN%CV_NEILOC( ILOC, 128 ) = 16
         ! Face 30
         SHAPE_FUN%CV_NEILOC( ILOC, 117 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 118 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 119 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 120 ) = 14
         ! Face 33
         SHAPE_FUN%CV_NEILOC( ILOC, 129 ) = 18
         SHAPE_FUN%CV_NEILOC( ILOC, 130 ) = 18
         SHAPE_FUN%CV_NEILOC( ILOC, 131 ) = 18
         SHAPE_FUN%CV_NEILOC( ILOC, 132 ) = 18
         ! Face 41
         SHAPE_FUN%CV_NEILOC( ILOC, 161 ) = 26
         SHAPE_FUN%CV_NEILOC( ILOC, 162 ) = 26
         SHAPE_FUN%CV_NEILOC( ILOC, 163 ) = 26
         SHAPE_FUN%CV_NEILOC( ILOC, 164 ) = 26
         ! External
         shape_fun%cv_neiloc( iloc, iext + 30 * igp + 1 : iext + 31 * igp ) = -1

         ILOC = 18
         ! Face 21
         SHAPE_FUN%CV_NEILOC( ILOC, 81 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 82 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 83 ) = 9
         SHAPE_FUN%CV_NEILOC( ILOC, 84 ) = 9
         ! Face 31
         SHAPE_FUN%CV_NEILOC( ILOC, 121 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 122 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 123 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 124 ) = 15
         ! Face 33
         SHAPE_FUN%CV_NEILOC( ILOC, 129 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 130 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 131 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 132 ) = 17
         ! Face 42
         SHAPE_FUN%CV_NEILOC( ILOC, 165 ) = 27
         SHAPE_FUN%CV_NEILOC( ILOC, 166 ) = 27
         SHAPE_FUN%CV_NEILOC( ILOC, 167 ) = 27
         SHAPE_FUN%CV_NEILOC( ILOC, 168 ) = 27
         ! External
         shape_fun%cv_neiloc( iloc, iext + 31 * igp + 1 : iext + 33 * igp ) = -1

         ILOC = 19
         ! Face 34
         SHAPE_FUN%CV_NEILOC( ILOC, 133 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 134 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 135 ) = 10
         SHAPE_FUN%CV_NEILOC( ILOC, 136 ) = 10
         ! Face 43
         SHAPE_FUN%CV_NEILOC( ILOC, 169 ) = 20
         SHAPE_FUN%CV_NEILOC( ILOC, 170 ) = 20
         SHAPE_FUN%CV_NEILOC( ILOC, 171 ) = 20
         SHAPE_FUN%CV_NEILOC( ILOC, 172 ) = 20
         ! Face 45
         SHAPE_FUN%CV_NEILOC( ILOC, 177 ) = 22
         SHAPE_FUN%CV_NEILOC( ILOC, 178 ) = 22
         SHAPE_FUN%CV_NEILOC( ILOC, 179 ) = 22
         SHAPE_FUN%CV_NEILOC( ILOC, 180 ) = 22
         ! External
         shape_fun%cv_neiloc( iloc, iext + 33 * igp + 1 : iext + 36 * igp ) = -1

         ILOC = 20
         ! Face 35
         SHAPE_FUN%CV_NEILOC( ILOC, 137 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 138 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 139 ) = 11
         SHAPE_FUN%CV_NEILOC( ILOC, 140 ) = 11
         ! Face 43
         SHAPE_FUN%CV_NEILOC( ILOC, 169 ) = 19
         SHAPE_FUN%CV_NEILOC( ILOC, 170 ) = 19
         SHAPE_FUN%CV_NEILOC( ILOC, 171 ) = 19
         SHAPE_FUN%CV_NEILOC( ILOC, 172 ) = 19
         ! Face 46
         SHAPE_FUN%CV_NEILOC( ILOC, 181 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 182 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 183 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 184 ) = 23
         ! Face 44
         SHAPE_FUN%CV_NEILOC( ILOC, 173 ) = 21
         SHAPE_FUN%CV_NEILOC( ILOC, 174 ) = 21
         SHAPE_FUN%CV_NEILOC( ILOC, 175 ) = 21
         SHAPE_FUN%CV_NEILOC( ILOC, 176 ) = 21
         ! External
         shape_fun%cv_neiloc( iloc, iext + 36 * igp + 1 : iext + 38 * igp ) = -1

         ILOC = 21
         ! Face 36
         SHAPE_FUN%CV_NEILOC( ILOC, 141 ) = 12
         SHAPE_FUN%CV_NEILOC( ILOC, 142 ) = 12
         SHAPE_FUN%CV_NEILOC( ILOC, 143 ) = 12
         SHAPE_FUN%CV_NEILOC( ILOC, 144 ) = 12
         ! Face 44
         SHAPE_FUN%CV_NEILOC( ILOC, 173 ) = 20
         SHAPE_FUN%CV_NEILOC( ILOC, 174 ) = 20
         SHAPE_FUN%CV_NEILOC( ILOC, 175 ) = 20
         SHAPE_FUN%CV_NEILOC( ILOC, 176 ) = 20
         ! Face 47
         SHAPE_FUN%CV_NEILOC( ILOC, 185 ) = 24
         SHAPE_FUN%CV_NEILOC( ILOC, 186 ) = 24
         SHAPE_FUN%CV_NEILOC( ILOC, 187 ) = 24
         SHAPE_FUN%CV_NEILOC( ILOC, 188 ) = 24
         ! External
         shape_fun%cv_neiloc( iloc, iext + 38 * igp + 1 : iext + 41 * igp ) = -1

         ILOC = 22
         ! Face 37
         SHAPE_FUN%CV_NEILOC( ILOC, 145 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 146 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 147 ) = 13
         SHAPE_FUN%CV_NEILOC( ILOC, 148 ) = 13
         ! Face 45
         SHAPE_FUN%CV_NEILOC( ILOC, 177 ) = 19
         SHAPE_FUN%CV_NEILOC( ILOC, 178 ) = 19
         SHAPE_FUN%CV_NEILOC( ILOC, 179 ) = 19
         SHAPE_FUN%CV_NEILOC( ILOC, 180 ) = 19
         ! Face 48
         SHAPE_FUN%CV_NEILOC( ILOC, 189 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 190 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 191 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 192 ) = 23
         ! Face 50
         SHAPE_FUN%CV_NEILOC( ILOC, 197 ) = 25
         SHAPE_FUN%CV_NEILOC( ILOC, 198 ) = 25
         SHAPE_FUN%CV_NEILOC( ILOC, 199 ) = 25
         SHAPE_FUN%CV_NEILOC( ILOC, 200 ) = 25
         ! External
         shape_fun%cv_neiloc( iloc, iext + 41 * igp + 1 : iext + 43 * igp ) = -1

         ILOC = 23
         ! Face 38
         SHAPE_FUN%CV_NEILOC( ILOC, 149 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 150 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 151 ) = 14
         SHAPE_FUN%CV_NEILOC( ILOC, 152 ) = 14
         ! Face 46
         SHAPE_FUN%CV_NEILOC( ILOC, 181 ) = 20
         SHAPE_FUN%CV_NEILOC( ILOC, 182 ) = 20
         SHAPE_FUN%CV_NEILOC( ILOC, 183 ) = 20
         SHAPE_FUN%CV_NEILOC( ILOC, 184 ) = 20
         ! Face 49
         SHAPE_FUN%CV_NEILOC( ILOC, 193 ) = 24
         SHAPE_FUN%CV_NEILOC( ILOC, 194 ) = 24
         SHAPE_FUN%CV_NEILOC( ILOC, 195 ) = 24
         SHAPE_FUN%CV_NEILOC( ILOC, 196 ) = 24
         ! Face 51
         SHAPE_FUN%CV_NEILOC( ILOC, 201 ) = 26
         SHAPE_FUN%CV_NEILOC( ILOC, 202 ) = 26
         SHAPE_FUN%CV_NEILOC( ILOC, 203 ) = 26
         SHAPE_FUN%CV_NEILOC( ILOC, 204 ) = 26
         ! Face 48
         SHAPE_FUN%CV_NEILOC( ILOC, 189 ) = 22
         SHAPE_FUN%CV_NEILOC( ILOC, 190 ) = 22
         SHAPE_FUN%CV_NEILOC( ILOC, 191 ) = 22
         SHAPE_FUN%CV_NEILOC( ILOC, 192 ) = 22
         ! External
         shape_fun%cv_neiloc( iloc, iext + 43 * igp + 1 : iext + 44 * igp ) = -1

         ILOC = 24
         ! Face 39
         SHAPE_FUN%CV_NEILOC( ILOC, 153 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 154 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 155 ) = 15
         SHAPE_FUN%CV_NEILOC( ILOC, 156 ) = 15
         ! Face 47
         SHAPE_FUN%CV_NEILOC( ILOC, 185 ) = 21
         SHAPE_FUN%CV_NEILOC( ILOC, 186 ) = 21
         SHAPE_FUN%CV_NEILOC( ILOC, 187 ) = 21
         SHAPE_FUN%CV_NEILOC( ILOC, 188 ) = 21
         ! Face 49
         SHAPE_FUN%CV_NEILOC( ILOC, 193 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 194 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 195 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 196 ) = 23
         ! Face 52
         SHAPE_FUN%CV_NEILOC( ILOC, 205 ) = 27
         SHAPE_FUN%CV_NEILOC( ILOC, 206 ) = 27
         SHAPE_FUN%CV_NEILOC( ILOC, 207 ) = 27
         SHAPE_FUN%CV_NEILOC( ILOC, 208 ) = 27
         ! External
         shape_fun%cv_neiloc( iloc, iext + 44 * igp + 1 : iext + 46 * igp ) = -1

         ILOC = 25
         ! Face 40
         SHAPE_FUN%CV_NEILOC( ILOC, 157 ) = 16
         SHAPE_FUN%CV_NEILOC( ILOC, 158 ) = 16
         SHAPE_FUN%CV_NEILOC( ILOC, 159 ) = 16
         SHAPE_FUN%CV_NEILOC( ILOC, 160 ) = 16
         ! Face 50
         SHAPE_FUN%CV_NEILOC( ILOC, 197 ) = 22
         SHAPE_FUN%CV_NEILOC( ILOC, 198 ) = 22
         SHAPE_FUN%CV_NEILOC( ILOC, 199 ) = 22
         SHAPE_FUN%CV_NEILOC( ILOC, 200 ) = 22
         ! Face 53
         SHAPE_FUN%CV_NEILOC( ILOC, 209 ) = 26
         SHAPE_FUN%CV_NEILOC( ILOC, 210 ) = 26
         SHAPE_FUN%CV_NEILOC( ILOC, 211 ) = 26
         SHAPE_FUN%CV_NEILOC( ILOC, 212 ) = 26
         ! External
         shape_fun%cv_neiloc( iloc, iext + 46 * igp + 1 : iext + 49 * igp ) = -1

         ILOC = 26
         ! Face 41
         SHAPE_FUN%CV_NEILOC( ILOC, 161 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 162 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 163 ) = 17
         SHAPE_FUN%CV_NEILOC( ILOC, 164 ) = 17
         ! Face 53
         SHAPE_FUN%CV_NEILOC( ILOC, 209 ) = 25
         SHAPE_FUN%CV_NEILOC( ILOC, 210 ) = 25
         SHAPE_FUN%CV_NEILOC( ILOC, 211 ) = 25
         SHAPE_FUN%CV_NEILOC( ILOC, 212 ) = 25
         ! Face 51
         SHAPE_FUN%CV_NEILOC( ILOC, 201 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 202 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 203 ) = 23
         SHAPE_FUN%CV_NEILOC( ILOC, 204 ) = 23
         ! Face 54
         SHAPE_FUN%CV_NEILOC( ILOC, 213 ) = 27
         SHAPE_FUN%CV_NEILOC( ILOC, 214 ) = 27
         SHAPE_FUN%CV_NEILOC( ILOC, 215 ) = 27
         SHAPE_FUN%CV_NEILOC( ILOC, 216 ) = 27
         ! External
         shape_fun%cv_neiloc( iloc, iext + 49 * igp + 1 : iext + 51 * igp ) = -1

         ILOC = 27
         ! Face 42
         SHAPE_FUN%CV_NEILOC( ILOC, 165 ) = 18
         SHAPE_FUN%CV_NEILOC( ILOC, 166 ) = 18
         SHAPE_FUN%CV_NEILOC( ILOC, 167 ) = 18
         SHAPE_FUN%CV_NEILOC( ILOC, 168 ) = 18
         ! Face 52
         SHAPE_FUN%CV_NEILOC( ILOC, 205 ) = 24
         SHAPE_FUN%CV_NEILOC( ILOC, 206 ) = 24
         SHAPE_FUN%CV_NEILOC( ILOC, 207 ) = 24
         SHAPE_FUN%CV_NEILOC( ILOC, 208 ) = 24
         ! Face 54
         SHAPE_FUN%CV_NEILOC( ILOC, 213 ) = 26
         SHAPE_FUN%CV_NEILOC( ILOC, 214 ) = 26
         SHAPE_FUN%CV_NEILOC( ILOC, 215 ) = 26
         SHAPE_FUN%CV_NEILOC( ILOC, 216 ) = 26
         ! External
         shape_fun%cv_neiloc( iloc, iext + 51 * igp + 1 : iext + 54 * igp ) = -1

      CASE DEFAULT; FLExit( " Invalid integer for cv_ele_type " )

      END SELECT Conditional_Type

      RETURN
    END SUBROUTINE VOLNEI



!!$ ============================================
!!$    This subroutine seems to be deprecated.
!!$ ============================================

    subroutine U_Volnei( cv_ele_type, Mdims, GIdims, shape_fun )
      !-----------------------------------------------------------------!
      !- This subroutine calculates U_ON_FACE, a logical that works    -!
      !- in a similar way of CV_ON_FACE.                               -!
      !-----------------------------------------------------------------!
      implicit none
      integer, intent( in ) :: cv_ele_type
      type (multi_dimensions), intent(in) :: Mdims
      type (multi_gi_dimensions), intent(in) :: GIdims
      type(multi_shape_funs), intent(inout) :: shape_fun
      integer :: cv_iloc, u_iloc, gi, u_jloc

      shape_fun%u_on_face = .false.


      Conditional_ElementType: if( cv_ele_type <= 2 ) then
         u_jloc = 0
         Loop_CV1: do cv_iloc = 1, Mdims%cv_nloc
            Loop_U1: do u_iloc = 1, Mdims%u_nloc
               Loop_GI1: do gi = 1, GIdims%scvngi
                  if ( shape_fun%cv_neiloc( cv_iloc, gi ) == -1 ) then
                     u_jloc = ( cv_iloc - 1 ) * Mdims%u_nloc + u_iloc
                     shape_fun%u_on_face( u_jloc, gi ) = .true.
                     if( u_iloc == 1 ) cycle Loop_CV
                  end if
               end do Loop_GI1
            end do Loop_U1
         end do Loop_CV1

      else
         u_jloc = 0
         Loop_CV2: do cv_iloc = 1, Mdims%cv_nloc
            Loop_U2: do u_iloc = 1, Mdims%u_nloc
               Loop_GI2: do gi = 1, GIdims%scvngi
                  if ( shape_fun%cv_neiloc( cv_iloc, gi ) == -1 ) then
                     u_jloc = ( cv_iloc - 1 ) * Mdims%u_nloc + u_iloc
                     shape_fun%u_on_face( u_jloc, gi ) = .true.
                  end if
               end do Loop_GI2
            end do Loop_U2
         end do Loop_CV2

      end if Conditional_ElementType

      return
    end subroutine U_Volnei

!!$ ===============
!!$ ===============

  end module Neighboor_Loc
