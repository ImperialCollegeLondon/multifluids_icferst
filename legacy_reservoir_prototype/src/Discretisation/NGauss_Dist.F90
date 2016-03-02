
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


  module NGauss
    use fldebug
    use multi_data_types

!!$ This is a hack that need to be fixed in the near future $!!
    logical :: NEW_HIGH_ORDER_VOL_QUADRATIC_ELE_QUADRATURE = .false., &
         NEW_QUADRATIC_ELE_QUADRATURE = .false.
!!$ This is a hack that need to be fixed in the near future $!!

  contains


    subroutine retrieve_ngi( GIdims, Mdims, cv_ele_type, QUAD_OVER_WHOLE_ELE, &
         scalar_nloc, vector_nloc )
      implicit none
      type( multi_GI_dimensions ), intent( inout ) :: GIdims
      type( multi_dimensions ), intent( in ) :: Mdims
      integer, intent( in ) :: cv_ele_type
!!$ If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
      logical, intent( in ) :: QUAD_OVER_WHOLE_ELE
      integer, intent( in ), optional :: scalar_nloc, vector_nloc
!!$Local variable
      integer :: cv_nloc, u_nloc
!!$   Volume_order & Surface_order are the volume and surface order of the integration 
!!$       used within each sub-quad/hex for CV approach. Default value -ve or 0 value is 
!!$       always 1pt quadrature (=1).
      integer, PARAMETER :: volume_order = 1, surface_order = 1
!!$    integer, PARAMETER :: volume_order = 2, surface_order = 2
!!$
!!$    Whole_ele_volume_order & Whole_ele_surface_order are the volume and surface
!!$        order of the integration used within each sub-quad/hex for 
!!$        QUAD_OVER_WHOLE_ELE=.true. -ve or 0 and take on the default value.
      integer, PARAMETER :: whole_ele_volume_order = 0, whole_ele_surface_order = 0
!!$    integer, PARAMETER :: whole_ele_volume_order = 1, whole_ele_surface_order = 1
!!$    integer, PARAMETER :: whole_ele_volume_order = 2, whole_ele_surface_order = 2

      if( present( scalar_nloc ) ) then
         cv_nloc = scalar_nloc
      else
         cv_nloc = Mdims%cv_nloc
      endif

      if( present( vector_nloc ) ) then
         u_nloc = vector_nloc
      else
         u_nloc = Mdims%u_nloc
      endif

      Conditional_EleType: Select Case( cv_ele_type )

      case( 1, 2 ) ! 1D
         Conditional_CV_NLOC_1D: Select Case( cv_nloc )
         case( 1 )
            GIdims%cv_ngi = 1  ; GIdims%scvngi = 2
         case( 2 )
            GIdims%cv_ngi = 12 ; GIdims%scvngi = 3
         case( 3 )
            GIdims%cv_ngi = 12 ; GIdims%scvngi = 4
            if( ( u_nloc == 4 ) .or. ( u_nloc == 5 ) ) &
                 GIdims%cv_ngi = 18
         case default; FLExit(" Invalid integer for cv_nloc ")
         end Select Conditional_CV_NLOC_1D
         GIdims%sbcvngi = 1 ; GIdims%nface = 2
!!$
      case( 3, 4 ) ! Triangles
         Conditional_CV_NLOC_2D_Tri: Select Case( cv_nloc )
         case( 3 ) ! Linear triangle
            Conditional_LinTriangle: if( QUAD_OVER_WHOLE_ELE ) then
               GIdims%cv_ngi = 3 ; GIdims%sbcvngi = 2 ; GIdims%scvngi = 2
               if( u_nloc == 6 ) then
                  GIdims%cv_ngi = 7 ; GIdims%sbcvngi = 3 ; GIdims%scvngi = 3
               elseif( u_nloc == 10 ) then
                  GIdims%cv_ngi = 14 ; GIdims%sbcvngi = 4 ; GIdims%scvngi = 4
               end if
!!$
               Select Case( whole_ele_volume_order )
               case( 1 )
                  GIdims%cv_ngi = 1
               case( 2 )
                  GIdims%cv_ngi = 3
               end Select
!!$
               Select Case( whole_ele_surface_order )
               case( 1 )
                  GIdims%sbcvngi = 1 ; GIdims%scvngi = 1
               case( 2 )
                  GIdims%sbcvngi = 2 ; GIdims%scvngi = 2
               end Select
!!$                
            else
               Select Case( volume_order )
               case( 1 )
                  GIdims%cv_ngi = 3
               case( 2 )
                  GIdims%cv_ngi = 3*4
               end Select
!!$
               Select Case( surface_order )
               case( 1 )
                  GIdims%sbcvngi = 2 ; GIdims%scvngi = 3
               case( 2 )
                  GIdims%sbcvngi = 2*2 ; GIdims%scvngi = 3*2
               end Select
!!$             
            end if Conditional_LinTriangle
!!$ ===
         case( 6 ) ! Quadratic triangle
            Conditional_QuadTriangle: if( QUAD_OVER_WHOLE_ELE ) then
               GIdims%cv_ngi = 7 ; GIdims%sbcvngi = 3 ; GIdims%scvngi = 3
               if( u_nloc == 10 ) & ! Use a quadratic interpolation pt set
                    GIdims%cv_ngi = 14 ; GIdims%sbcvngi = 4 ; GIdims%scvngi = 4
!!$
               Select Case( whole_ele_volume_order )
               case( 1 )
                  GIdims%cv_ngi = 1
               case( 2 )
                  GIdims%cv_ngi = 3
               case( 3 )
                  GIdims%cv_ngi = 7
               end Select
!!$
               Select Case( whole_ele_surface_order )
               case( 1 )
                  GIdims%sbcvngi = 1 ; GIdims%scvngi = 1
               case( 2 )
                  GIdims%sbcvngi = 2 ; GIdims%scvngi = 2
               case( 3 )
                  GIdims%sbcvngi = 3 ; GIdims%scvngi = 3
               end Select
!!$                
            else
               Select Case( volume_order )
               case( 1 )
                  GIdims%cv_ngi = 12
               case( 2 )
                  GIdims%cv_ngi = 12*4
               end Select
!!$
               Select Case( surface_order )
               case( 1 )
                  GIdims%sbcvngi = 4 ; GIdims%scvngi = 12
               case( 2 )
                  GIdims%sbcvngi = 4*2 ; GIdims%scvngi = 12*2
               end Select
!!$             
            end if Conditional_QuadTriangle
         case default; FLExit(" Invalid integer for cv_nloc ")
         end Select Conditional_CV_NLOC_2D_Tri
         GIdims%nface = 3
!!$
      case( 5, 6 ) ! Quads       
         Conditional_CV_NLOC_2D_Quad: Select Case( cv_nloc )
         case( 4 ) ! Bi-linear Quad
            Conditional_BiLinQuad: if( QUAD_OVER_WHOLE_ELE ) then
               GIdims%cv_ngi = 4 ; GIdims%sbcvngi = 2 ; GIdims%scvngi = 2
!!$
               Select Case( whole_ele_volume_order )
               case( 1 )
                  GIdims%cv_ngi = 1
               case( 2 )
                  GIdims%cv_ngi = 4
               case( 3 )
                  GIdims%cv_ngi = 9
               end Select
!!$
               Select Case( whole_ele_surface_order )
               case( 1 )
                  GIdims%sbcvngi = 1 ; GIdims%scvngi = 1
               case( 2 )
                  GIdims%sbcvngi = 2 ; GIdims%scvngi = 2
               case( 3 )
                  GIdims%sbcvngi = 3 ; GIdims%scvngi = 3
               end Select
!!$                
            else
               Select Case( volume_order )
               case( 1 )
                  GIdims%cv_ngi = 4
               case( 2 )
                  GIdims%cv_ngi = 4*4
               end Select
!!$
               Select Case( surface_order )
               case( 1 )
                  GIdims%sbcvngi = 2 ; GIdims%scvngi = 4
               case( 2 )
                  GIdims%sbcvngi = 4 ; GIdims%scvngi = 4*2
               end Select
!!$             
            end if Conditional_BiLinQuad
!!$ ===
         case( 9 ) ! Bi-Quadratic Quad
            Conditional_BiQuadQuad: if( QUAD_OVER_WHOLE_ELE ) then
               GIdims%cv_ngi = 9 ; GIdims%sbcvngi = 3 ; GIdims%scvngi = 3
!!$
               Select Case( whole_ele_volume_order )
               case( 1 )
                  GIdims%cv_ngi = 1
               case( 2 )
                  GIdims%cv_ngi = 4
               case( 3 )
                  GIdims%cv_ngi = 9
               end Select
!!$
               Select Case( whole_ele_surface_order )
               case( 1 )
                  GIdims%sbcvngi = 1 ; GIdims%scvngi = 1
               case( 2 )
                  GIdims%sbcvngi = 2 ; GIdims%scvngi = 2
               case( 3 )
                  GIdims%sbcvngi = 3 ; GIdims%scvngi = 3
               end Select
!!$                
            else
               Select Case( volume_order )
               case( 1 )
                  GIdims%cv_ngi = 16
               case( 2 )
                  GIdims%cv_ngi = 16*4
               case( 3 )
                  GIdims%cv_ngi = 16*9
               end Select
!!$
               Select Case( surface_order )
               case( 1 )
                  GIdims%sbcvngi = 4 ; GIdims%scvngi = 16
               case( 2 )
                  GIdims%sbcvngi = 4*2 ; GIdims%scvngi = 16*2
               case( 3 )
                  GIdims%sbcvngi = 4*3 ; GIdims%scvngi = 16*3
               end Select
!!$             
            end if Conditional_BiQuadQuad
         case default; FLExit(" Invalid integer for cv_nloc ")
         end Select Conditional_CV_NLOC_2D_Quad
         GIdims%nface = 4
!!$

      case( 7, 8 ) ! Tetrahedra       
         Conditional_CV_NLOC_3D_Tets: Select Case( cv_nloc )
         case( 4 ) ! Linear
            Conditional_LinTets: if( QUAD_OVER_WHOLE_ELE ) then
               GIdims%cv_ngi = 4 ; GIdims%sbcvngi = 3 ; GIdims%scvngi = 3
               if( u_nloc == 10 ) & ! Use a quadratic interpolation pt set
                    GIdims%cv_ngi = 11 ; GIdims%sbcvngi = 7 ; GIdims%scvngi = 7
!!$
               Select Case( whole_ele_volume_order )
               case( 1 )
                  GIdims%cv_ngi = 1
               case( 2 )
                  GIdims%cv_ngi = 4
               case( 3 )
                  GIdims%cv_ngi = 11
               end Select
!!$
               Select Case( whole_ele_surface_order )
               case( 1 )
                  GIdims%sbcvngi = 1 ; GIdims%scvngi = 1
               case( 2 )
                  GIdims%sbcvngi = 3 ; GIdims%scvngi = 3
               case( 3 )
                  GIdims%sbcvngi = 7 ; GIdims%scvngi = 7
               end Select
!!$                
            else
               Select Case( volume_order )
               case( 1 )
                  GIdims%cv_ngi = 4
               case( 2 )
                  GIdims%cv_ngi = 4*8
               end Select
!!$
               Select Case( surface_order )
               case( 1 )
                  GIdims%sbcvngi = 3 ; GIdims%scvngi = 6
               case( 2 )
                  GIdims%sbcvngi = 3*4 ; GIdims%scvngi = 6*4
               end Select
!!$             
            end if Conditional_LinTets
!!$ ===
         case( 10 ) ! Quadratic
            Conditional_QuadTets: if( QUAD_OVER_WHOLE_ELE ) then
               GIdims%cv_ngi = 11 ; GIdims%sbcvngi = 7 ; GIdims%scvngi = 7
!!$
               Select Case( whole_ele_volume_order )
               case( 1 )
                  GIdims%cv_ngi = 1
               case( 2 )
                  GIdims%cv_ngi = 4
               case( 3 )
                  GIdims%cv_ngi = 11
               end Select
!!$
               Select Case( whole_ele_surface_order )
               case( 1 )
                  GIdims%sbcvngi = 1 ; GIdims%scvngi = 1
               case( 2 )
                  GIdims%sbcvngi = 3 ; GIdims%scvngi = 3
               case( 3 )
                  GIdims%sbcvngi = 7 ; GIdims%scvngi = 7
               end Select
!!$                
            else
               Select Case( volume_order )
               case( 1 )
                  GIdims%cv_ngi = 8*4*1  ! (1x1x1)
               case( 2 )
                  GIdims%cv_ngi = 8*4*8  ! (2x2x2)
               case( 3 )
                  GIdims%cv_ngi = 8*4*27 ! (3x3x3)
               end Select
!!$
               Select Case( surface_order )
               case( 1 )
                  GIdims%sbcvngi = 12 ! 1x12 (sngi x cv_faces)
                  GIdims%scvngi = 48  ! 6x8x1 (cv_faces x tets x sngi)
               case( 2 )
                  GIdims%sbcvngi = 12*4  ! 4x12 (sngi x cv_faces)
                  GIdims%scvngi = 48*4   ! 6x8x4 (cv_faces x hexs x sngi)
               end Select
!!$
               if( NEW_QUADRATIC_ELE_QUADRATURE ) then ! new 1 pt quadrature
                  GIdims%sbcvngi = 6 ; GIdims%scvngi = 60-24 
                  if( NEW_HIGH_ORDER_VOL_QUADRATIC_ELE_QUADRATURE ) then
                     GIdims%cv_ngi = 10*4 ! 1 pt quadrature put CV in an element
                  else
                     GIdims%cv_ngi = 10 ! 1 pt quadrature put CV in an element
                  endif
               endif
!!$
            endif Conditional_QuadTets
         case default; FLExit(" Invalid integer for cv_nloc ")
         end Select Conditional_CV_NLOC_3D_Tets
         GIdims%nface = 4
!!$ 

      case( 9, 10 ) ! Hexahedra      
         Conditional_CV_NLOC_3D_Hexs: Select Case( cv_nloc )
         case( 8 ) ! Tri-linear Hex
            Conditional_TriLinHex: if( QUAD_OVER_WHOLE_ELE ) then
               GIdims%cv_ngi = 8 ; GIdims%sbcvngi = 4 ; GIdims%scvngi = 4
!!$
               Select Case( whole_ele_volume_order )
               case( 1 )
                  GIdims%cv_ngi = 1
               case( 2 )
                  GIdims%cv_ngi = 8
               case( 3 )
                  GIdims%cv_ngi = 27
               end Select
!!$
               Select Case( whole_ele_surface_order )
               case( 1 )
                  GIdims%sbcvngi = 1 ; GIdims%scvngi = 1
               case( 2 )
                  GIdims%sbcvngi = 4 ; GIdims%scvngi = 4
               case( 3 )
                  GIdims%sbcvngi = 9 ; GIdims%scvngi = 9
               end Select
!!$                
            else
               Select Case( volume_order )
               case( 1 )
                  GIdims%cv_ngi = 8
               case( 2 )
                  GIdims%cv_ngi = 8*8
               end Select
!!$
               Select Case( surface_order )
               case( 1 )
                  GIdims%sbcvngi = 4 ; GIdims%scvngi = 12
               case( 2 )
                  GIdims%sbcvngi = 4*4 ; GIdims%scvngi = 12*4
               end Select
!!$             
            end if Conditional_TriLinHex
!!$ ===
         case( 27 ) ! Tri-Quad Hex
            Conditional_TriQuadHex: if( QUAD_OVER_WHOLE_ELE ) then
               GIdims%cv_ngi = 27 ; GIdims%sbcvngi = 9 ; GIdims%scvngi = 9
!!$
               Select Case( whole_ele_volume_order )
               case( 1 )
                  GIdims%cv_ngi = 1
               case( 2 )
                  GIdims%cv_ngi = 8
               case( 3 )
                  GIdims%cv_ngi = 27
               end Select
!!$
               Select Case( whole_ele_surface_order )
               case( 1 )
                  GIdims%sbcvngi = 1 ; GIdims%scvngi = 1
               case( 2 )
                  GIdims%sbcvngi = 4 ; GIdims%scvngi = 4
               case( 3 )
                  GIdims%sbcvngi = 9 ; GIdims%scvngi = 9
               end Select
!!$                
            else
               Select Case( volume_order )
               case( 1 )
                  GIdims%cv_ngi = 64 
               case( 2 )
                  GIdims%cv_ngi = 64*8 
               case( 3 )
                  GIdims%cv_ngi = 64*27
               end Select
!!$
               Select Case( surface_order )
               case( 1 )
                  GIdims%sbcvngi = 16 ; GIdims%scvngi = 12*8  
               case( 2 )
                  GIdims%sbcvngi = 16*4 ; GIdims%scvngi = 12*8*4  
               case( 3 )
                  GIdims%sbcvngi = 16*9 ; GIdims%scvngi = 12*8*9
               end Select
!!$
            endif Conditional_TriQuadHex
!!$             
         end Select Conditional_CV_NLOC_3D_Hexs
         GIdims%nface = 6

      case default; FLExit( "Invalid integer for cv_ele_type" )

      end Select Conditional_EleType

      if( .not. QUAD_OVER_WHOLE_ELE) then
         if( cv_ele_type > 2 ) &
              GIdims%scvngi = GIdims%scvngi + GIdims%nface * GIdims%sbcvngi 
      endif

      return
    end subroutine retrieve_ngi


!!$ ===============
!!$ ===============


    subroutine lagrot( weit, quadpos, ndgi, getndp )
!!$  This computes the weight and points for standard Gaussian quadrature.
!!$      If (GETNDP == T) then get the position of the nodes and neglect the weights.
      implicit none
      integer, intent( in ) :: ndgi
      real, dimension( : ), intent( inout ) :: weit, quadpos
      logical, intent( in ) :: getndp
      ! Local variables
      logical :: weight
      integer :: ig

      weit = 0. ; quadpos = 0.

      Conditional_GETNDP: if( .not. getndp ) then
         weight = .true.
         do ig = 1, ndgi
            weit( ig ) = rgptwe( ig, ndgi, weight )
         end do

         weight = .false.
         do ig = 1, ndgi
            quadpos( ig ) = rgptwe( ig, ndgi, weight )
         end do

      else
         if( ndgi == 1 ) then
            quadpos( ndgi ) = 0.
         else
            do ig = 1, ndgi
               quadpos( ig ) = -1. + 2. * real( ig - 1 ) / real( ndgi - 1 )
            end do
         end if

      end if Conditional_GETNDP

      return
    end subroutine lagrot



    REAL FUNCTION RGPTWE( IG, ND, WEIGHT )

      IMPLICIT NONE

!!$  If WEIGHT is TRUE in function RGPTWE then return the Gauss-pt weight else
!!$      return the Gauss-pt.  There are ND Gauss points -- we are looking for either
!!$      the weight or the x-coord of the IG'th Gauss point.

      INTEGER :: IG, ND
      LOGICAL :: WEIGHT


      Conditional_Weight: IF( WEIGHT ) THEN ! Gauss points weights

         Conditional_ND1: SELECT CASE( ND )

         CASE( 1 ) ; RGPTWE = 2.0

         CASE( 2 ) ; RGPTWE = 1.0

         CASE( 3 )
            SELECT CASE( IG )
            CASE( 1, 3 ) ; RGPTWE = 0.555555555555556
            CASE( 2 )    ; RGPTWE = 0.888888888888889
            END SELECT

         CASE( 4 )
            SELECT CASE( IG )
            CASE( 1, 4 ) ; RGPTWE = 0.347854845137454
            CASE( 2, 3 ) ; RGPTWE = 0.652145154862546
            END SELECT

         CASE( 5 )
            SELECT CASE( IG )
            CASE( 1, 5 ) ; RGPTWE = 0.236926885056189
            CASE( 2, 4 ) ; RGPTWE = 0.478628670499366
            CASE( 3 )    ; RGPTWE = 0.568888888888889
            END SELECT

         CASE( 6 )
            SELECT CASE( IG )
            CASE( 1, 6 ) ; RGPTWE = 0.171324492379170
            CASE( 2, 5 ) ; RGPTWE = 0.360761573048139
            CASE( 3, 4 ) ; RGPTWE = 0.467913934572691
            END SELECT

         CASE( 7 )
            SELECT CASE( IG )
            CASE( 1, 7 ) ; RGPTWE = 0.129484966168870
            CASE( 2, 6 ) ; RGPTWE = 0.279705391489277
            CASE( 3, 5 ) ; RGPTWE = 0.381830050505119
            CASE( 4 )    ; RGPTWE = 0.417959183673469
            END SELECT

         CASE( 8 )
            SELECT CASE( IG )
            CASE( 1, 8 ) ; RGPTWE = 0.101228536290376
            CASE( 2, 7 ) ; RGPTWE = 0.222381034453374
            CASE( 3, 6 ) ; RGPTWE = 0.313706645877877
            CASE( 4, 5 ) ; RGPTWE = 0.362683783378362
            END SELECT

         CASE( 9 )
            SELECT CASE( IG )
            CASE( 1, 9 ) ; RGPTWE = 0.081274388361574
            CASE( 2, 8 ) ; RGPTWE = 0.180648160694857
            CASE( 3, 7 ) ; RGPTWE = 0.260610696402935
            CASE( 4, 6 ) ; RGPTWE = 0.312347077040003
            CASE( 5 )    ; RGPTWE = 0.330239355001260
            END SELECT

         CASE( 10 )
            SELECT CASE( IG )
            CASE( 1, 10 ) ; RGPTWE = 0.066671344308688
            CASE( 2, 9 )  ; RGPTWE = 0.149451349150581
            CASE( 3, 8 )  ; RGPTWE = 0.219086362515982
            CASE( 4, 7 )  ; RGPTWE = 0.269266719309996
            CASE( 5, 6 )  ; RGPTWE = 0.295524224714753
            END SELECT

         CASE DEFAULT; FLExit( "In Lagrot subrt - wrong integer for NDGI_1 " )

         END SELECT Conditional_ND1

      ELSE
         Conditional_ND1: SELECT CASE( ND ) ! Gauss points

         CASE( 1 ) ; RGPTWE = 0.0

         CASE( 2 ) ; RGPTWE = 0.577350269189626

         CASE( 3 )
            SELECT CASE( IG )
            CASE( 1, 3 ) ; RGPTWE = 0.774596669241483
            CASE( 2 )    ; RGPTWE = 0.0
            END SELECT

         CASE( 4 )
            SELECT CASE( IG )
            CASE( 1, 4 ) ; RGPTWE = 0.861136311594953
            CASE( 2, 3 ) ; RGPTWE = 0.339981043584856
            END SELECT

         CASE( 5 )
            SELECT CASE( IG )
            CASE( 1, 5 ) ; RGPTWE = 0.906179845938664
            CASE( 2, 4 ) ; RGPTWE = 0.538469310105683
            CASE( 3 )    ; RGPTWE = 0.0
            END SELECT

         CASE( 6 )
            SELECT CASE( IG )
            CASE( 1, 6 ) ; RGPTWE = 0.932469514203152
            CASE( 2, 5 ) ; RGPTWE = 0.661209386466265
            CASE( 3, 4 ) ; RGPTWE = 0.238619186083197
            END SELECT

         CASE( 7 )
            SELECT CASE( IG )
            CASE( 1, 7 ) ; RGPTWE = 0.949107912342759
            CASE( 2, 6 ) ; RGPTWE = 0.741531185599394
            CASE( 3, 5 ) ; RGPTWE = 0.405845151377397
            CASE( 4 )    ; RGPTWE = 0.0
            END SELECT

         CASE( 8 )
            SELECT CASE( IG )
            CASE( 1, 8 ) ; RGPTWE = 0.960289856497536
            CASE( 2, 7 ) ; RGPTWE = 0.796666477413627
            CASE( 3, 6 ) ; RGPTWE = 0.525532409916329
            CASE( 4, 5 ) ; RGPTWE = 0.183434642495650
            END SELECT

         CASE( 9 )
            SELECT CASE( IG )
            CASE( 1, 9 ) ; RGPTWE = 0.968160239507626
            CASE( 2, 8 ) ; RGPTWE = 0.836031107326636
            CASE( 3, 7 ) ; RGPTWE = 0.613371432700590
            CASE( 4, 6 ) ; RGPTWE = 0.324253423403809
            CASE( 5)     ; RGPTWE = 0.0
            END SELECT

         CASE( 10 )
            SELECT CASE( IG )
            CASE( 1, 10 ) ; RGPTWE = 0.973906528517172
            CASE( 2, 9 )  ; RGPTWE = 0.865063366688985
            CASE( 3, 8 )  ; RGPTWE = 0.679409568299024
            CASE( 4, 7 )  ; RGPTWE = 0.433395394129247
            CASE( 5, 6 )  ; RGPTWE = 0.148874338981631
            END SELECT

         CASE DEFAULT; FLExit( "In Lagrot subrt - wrong integer for NDGI_2 " )

         END SELECT Conditional_ND1

         IF( IG <= INT( ( ND / 2 ) + 0.1 )) RGPTWE = -RGPTWE

      END IF Conditional_Weight

    END FUNCTION RGPTWE

!!$ ===============
!!$ ===============

    real function lagran( diff, lx, inod, ndnod, nodpos )
!!$   This return the Lagrange poly assocaited with node INOD at point LX.
!!$      If DIFF then send back the value of this poly differentiated.
      implicit none
      logical :: diff
      real :: lx
      integer :: inod, ndnod
      real, dimension(:), allocatable :: nodpos
!!$  Local variables
      integer :: n, k, i, j
      real :: denomi, over, over1


      allocate( nodpos( 0 : ndnod - 1 ) )
      n = ndnod - 1
      k = inod - 1

      denomi = 1.
      do i = 0, k - 1
         denomi = denomi * ( nodpos( k ) - nodpos( i ) )
      end do
      do i = k + 1, n
         denomi = denomi * ( nodpos( k ) - nodpos( i ) )
      end do

      Conditional_DIFF1: if( .not. diff ) then
         over = 1.
         do i = 0, k - 1
            over = over * ( lx - nodpos( i ) )
         end do
         do i = k + 1, n
            over = over * ( lx - nodpos( i ) )
         end do
         lagran = over / denomi

      else ! Conditional_DIFF1
         over = 0.
         do j = 0, n
            if( j /= k ) then
               over1 = 1.
               do i = 0, k - 1
                  if( j /= i ) over1 = over1 * ( lx - nodpos( i ) )
               end do
               do i = k + 1, n
                  if( j /= i ) over1 = over1 * ( lx - nodpos( i ) )
               end do
               over = over + over1
            end if
         end do
         lagran = over / denomi

      end if Conditional_DIFF1

      deallocate( nodpos )

      return
    end function lagran

  end module NGauss
