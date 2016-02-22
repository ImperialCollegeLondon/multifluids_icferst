
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
    GIdims%cv_ngi_short = GIdims%cv_ngi


    return
  end subroutine retrieve_ngi

end module NGauss
