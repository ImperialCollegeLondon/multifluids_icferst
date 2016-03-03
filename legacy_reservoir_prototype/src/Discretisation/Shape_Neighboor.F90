
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

  module Shape_Neighboors
    use fldebug
    use multi_data_types
    use NGauss

  contains

!!$ ===============
!!$ ===============

  integer function Get_NwiCel( d3, nloc )
    implicit none
    logical :: d3
    integer :: nloc, nwicel

    Conditional_Dimensionality: if( d3 ) then
       Select Case ( nloc )
       case( 4 ) ! Linear tets
          nwicel = 4
       case( 8 ) ! Linear hex
          nwicel = 1
       case( 10 ) ! Quadratic tets
          nwicel = 5
       case( 20 ) ! Cubic tets
          nwicel = 6 
       case( 27 ) ! Quadratic hex
          nwicel = 3
       case default
          nwicel = 2
       end Select
!!$
    else
       Select Case ( nloc )
       case( 3 ) ! Linear triangle
          nwicel = 4
       case( 4 ) ! Linear quad
          nwicel = 1
       case( 6 ) ! Quadratic triangle
          nwicel = 5
       case( 9 ) ! Quadratic quad
          nwicel = 3
       case( 10 ) ! Cubic triangle
          nwicel = 6
       case default
          nwicel = 2
       end Select
!!$
    end if Conditional_Dimensionality

    Get_NwiCel = nwicel

    return
  end function Get_NwiCel

!!$ ===============
!!$ ===============

  SUBROUTINE SHAPE(LOWQUA,NGI,NLOC,MLOC, SNGI,SNLOC,SMLOC,    &
       M,MLX,MLY,MLZ,WEIGHT,N,NLX,NLY,NLZ,                    & 
       SWEIGH,SN,SNLX,SNLY, SM,SMLX,SMLY,                     &
       NWICEL,D3)
    implicit none
    LOGICAL, INTENT(IN):: LOWQUA
    INTEGER, INTENT(IN):: NGI,NLOC,MLOC,SNGI,SNLOC,SMLOC
    REAL, dimension(:,:), INTENT(OUT):: M,MLX,MLY,MLZ
    REAL, dimension(:), INTENT(OUT):: WEIGHT
    REAL, dimension(:,:), INTENT(OUT):: N,NLX,NLY,NLZ
    REAL, dimension(:), INTENT(OUT):: SWEIGH
    REAL, dimension(:,:), INTENT(OUT):: SN,SNLX,SNLY
    REAL, dimension(:,:), INTENT(OUT):: SM,SMLX,SMLY
    INTEGER, INTENT(IN):: NWICEL
    LOGICAL, INTENT(IN):: D3
    ! Local variable
    INTEGER :: IPOLY, IQADRA, gi, gj, ggi, i, j, ii

    Conditional_NWICEL: Select Case( NWICEL )
    case( 1 )
       IF(.NOT.D3) THEN
          CALL RE2DN4(LOWQUA,NGI,0,NLOC,MLOC, &
               M,WEIGHT,N,NLX,NLY,            &
               SNGI,SNLOC,SWEIGH,SN,SNLX,     &
               m(:,1),m(:,1))
       ELSE
          CALL RE3DN8(LOWQUA,NGI,0,NLOC,MLOC,  &
               M,WEIGHT,N,NLX,NLY,NLZ,         &
               SNGI,SNLOC,SWEIGH,SN,SNLX,SNLY, &
               m(:,1),m(:,1),m(:,1))
       ENDIF
!!$
    case( 2 )
       FLAbort("Option not available yet")
!!$
    case( 3 )
       IF(.NOT.D3) THEN
          CALL RE2DN9(LOWQUA,NGI,0,NLOC,MLOC, &
               M,WEIGHT,N,NLX,NLY,            &
               m(:,1),m(:,1))
          sweigh = 0.0 ; sn = 0.0 ; snlx = 0.0
          do gi=1,3
             do gj=1,3
                ggi=(gj-1)*3+gi
                sweigh(gi)=sweigh(gi)+WEIGHT(ggi)
                do i=1,3
                   do j=1,3
                      ii=(j-1)*3+i
                      sn(i,gi)=sn(i,gi)+n(ii,ggi)
                      snlx(i,gi)=snlx(i,gi)+nlx(ii,ggi)
                   end do
                end do
             end do
          end do
       ELSE ! LAGRANGE 27 NODE 3-D ELEMENT -BILINEAR PRESSURE
          CALL RE3D27(LOWQUA,NGI,0,NLOC,MLOC, &
               M,WEIGHT,N,NLX,NLY,NLZ,        &
               m(:,1),m(:,1),m(:,1))
          CALL RE2DN9(LOWQUA,SNGI,0,SNLOC,MLOC, &
               M,SWEIGH,SN,SNLX,SNLY,           &
               m(:,1),m(:,1))
       ENDIF
!!$
    case( 4:6 ) ! works for linear or quadratic triangles or tets (also cubic triangles)...
       CALL TR2or3DQU(NGI,NLOC,MLOC,        &
            M,MLX,MLY,MLZ,                  &
            WEIGHT,N,NLX,NLY,NLZ,           &
            SNGI,SNLOC,SWEIGH,SN,SNLX,SNLY, &
            SMLOC,                          &
            SM,SMLX,SMLY,D3)
!!$
    case( 100: )! A Spectal element using Legendra, Lagrange or Chebichef polynomials.
       CALL SPECTR(NGI,NLOC,MLOC, &
            M,WEIGHT,N,NLX,NLY,NLZ,D3,.NOT.D3, IPOLY,IQADRA)
!!$
    case default
       FLAbort("Option not found")
!!$
    end Select Conditional_NWICEL

    return
  END SUBROUTINE SHAPE

!!$ ===============
!!$ ===============

!!$ ===============
!!$ ===============

end module Shape_Neighboors
