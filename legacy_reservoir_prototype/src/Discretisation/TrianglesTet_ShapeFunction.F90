
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
!!!!               TRIANGLES AND TETRAHEDRA                !!!!
!!!!=======================================================!!!!
  
  
  module TriangleTet_ShapeFunctions

    use fldebug
    use All_ShaTri

  contains


!!$ ===============
!!$ ===============


    SUBROUTINE LinearQuadratic_TrianglesTetrahedra( NGI, NLOC, MLOC, &
         M, MLX, MLY, MLZ,                                           &
         WEIGHT, N, NLX, NLY, NLZ,                                   &
         SNGI, SNLOC, SWEIGH, SN, SNLX, SNLY,                        &
         SMLOC,                                                      &
         SM, SMLX, SMLY, D3 )
!!$      This subroutine defines the shape functions M and N and their
!!$           derivatives at the Gauss points for quadratic elements. For 3-D FLOW.
      implicit none
      INTEGER , intent(in) :: SNGI, SNLOC, SMLOC, NLOC, MLOC, NGI
      real, dimension(:), intent(inout) :: SWEIGH, WEIGHT
      real, dimension(:,:), intent(inout) :: SN, SNLX, SNLY, SM, SMLX, SMLY, M, MLX, MLY, MLZ, N, NLX, NLY, NLZ
      LOGICAL , intent(in) ::D3
!!$      Local variables
      integer, parameter :: n50 = 50, n500 = 500, n25 = 25
      integer :: IQADRA, IPOLY
      real, dimension(:), allocatable :: rub, L1, L2, L3, L4
      real, dimension(:,:), allocatable :: rub2
      LOGICAL :: DD3, base_order

      allocate( RUB(n500), L1(n50), L2(n50), L3(n50), L4(n50), rub2(n25,n25) )

!!$   NB LXP(I) AND LYP(I) ARE THE LOCAL X AND Y COORDS OF NODAL POINT I

!!$   Get the quadrature positions and weights for TRIANGLES or TETRAHEDRA
      DD3 = D3
      CALL TRIQUAold( L1, L2, L3, L4, WEIGHT, DD3, NGI )

!!$   Work out the shape functions and there derivatives...
      CALL SHATRIold( L1, L2, L3, L4, WEIGHT, DD3,  &
           NLOC, NGI,                               &
           N, NLX, NLY, NLZ )

!!$   Re-arrange ordering for quadratic elements...
      Conditional_Rearranging_QuadElems: if(d3) then
         if(nloc==10) then
            base_order=.true.
            if(base_order) then !!$ order so that the 1st nodes are on the base...
               call base_order_tet(n,nloc,ngi)
               call base_order_tet(nlx,nloc,ngi)
               call base_order_tet(nly,nloc,ngi)
               call base_order_tet(nlz,nloc,ngi)
            endif
         endif
      else
         if((nloc==6).or.(nloc==10)) then
            base_order=.true.
            if(base_order) then !!$ order so that the 1st nodes are on the base...
               call base_order_tri(n,nloc,ngi)
               call base_order_tri(nlx,nloc,ngi)
               call base_order_tri(nly,nloc,ngi)
            endif
         endif
      endif Conditional_Rearranging_QuadElems

      CALL SHATRIold( L1, L2, L3, L4, WEIGHT, DD3, &
           MLOC, NGI,                              &
           M, MLX, MLY, MLZ )

      Conditional_SNGI: IF(SNGI.GT.0) THEN
         Conditional_Rearranging_SurfQuadElems: IF(D3) THEN
            DD3=.FALSE.
            CALL TRIQUAold( L1, L2, L3, L4, SWEIGH, DD3, SNGI)

!!$   Work out the shape functions and there derivatives...
            CALL SHATRIold( L1, L2, L3, L4, SWEIGH, DD3, &
                 SNLOC, SNGI,                            &
                 SN, SNLX, SNLY, RUB2 )
            if( (snloc==6) .or. (snloc==10) ) then
               base_order = .true.
               if(base_order) then !!$ order so that the 1st nodes are on the base...
                  call base_order_tri(sn,snloc,sngi)
                  call base_order_tri(snlx,snloc,sngi)
                  call base_order_tri(snly,snloc,sngi)
               endif
            endif
            CALL SHATRIold( L1, L2, L3, L4, SWEIGH, DD3, &
                 SMLOC, SNGI,                            &
                 SM, SMLX, SMLY, RUB2 )

         ELSE
            IQADRA = 1 !!$   IQADRA=1 corresponds to Gaussian quadrature.
            IPOLY = 1  !!$   IPOLY=1 is for Lagrange polynomials.
            CALL SPECTR( SNGI, SNLOC, 0,     &
                 RUB2, SWEIGH, SN, SNLX, RUB2, RUB2, .FALSE., .FALSE., IPOLY, IQADRA )

            if(.false.) CALL SPECTR( SNGI, SMLOC, 0, &
                 RUB2, SWEIGH, SM, SMLX, RUB2, RUB2, .FALSE., .FALSE., IPOLY, IQADRA )

         ENDIF Conditional_Rearranging_SurfQuadElems

      ENDIF Conditional_SNGI

!!$  The weights need to sum to 0.5 in 2D triangles and 1./6. in 3D
      IF(D3) THEN
         WEIGHT = 1. * WEIGHT ; SWEIGH = 0.5 * SWEIGH
      ELSE 
         WEIGHT = 0.5 * WEIGHT
      ENDIF

!!$   Deallocating
      deallocate( RUB, L1, L2, L3, L4, rub2 )

      return
    end SUBROUTINE LinearQuadratic_TrianglesTetrahedra

!!$ ===============
!!$ ===============

    SUBROUTINE TRIQUAold(L1, L2, L3, L4, WEIGHT, D3,NGI)
      ! This sub calculates the local corrds L1, L2, L3, L4 and
      ! weights at the quadrature points.
      ! If D3 it does this for 3Dtetrahedra elements else
      ! triangular elements.
      IMPLICIT NONE
      INTEGER , intent(in):: NGI
      LOGICAL , intent(in) :: D3
      REAL , dimension(:) , intent(inout) ::L1, L2, L3, L4, WEIGHT
      ! Local variables...
      REAL :: ALPHA,BETA
      REAL :: ALPHA1,BETA1
      REAL :: ALPHA2,BETA2
      real :: rsum
      INTEGER I
      !
      IF(D3) THEN
         ! this is for a tetrahedra element...
         ! This is for one point.
         IF(NGI.EQ.1) THEN
            ! Degree of precision is 1
            DO I=1,NGI
               L1(I)=0.25
               L2(I)=0.25
               L3(I)=0.25
               L4(I)=0.25
               WEIGHT(I)=1.0
            END DO
         ENDIF

         IF(NGI.EQ.4) THEN
            ! Degree of precision is 2
            ALPHA=0.58541020
            BETA=0.13819660
            DO I=1,NGI
               L1(I)=BETA
               L2(I)=BETA
               L3(I)=BETA
               L4(I)=BETA
               WEIGHT(I)=0.25
            END DO
            L1(1)=ALPHA
            L2(2)=ALPHA
            L3(3)=ALPHA
            L4(4)=ALPHA
         ENDIF

         IF(NGI.EQ.5) THEN
            ! Degree of precision is 3
            L1(1)=0.25
            L2(1)=0.25
            L3(1)=0.25
            L4(1)=0.25
            WEIGHT(1)=-4./5.
            !
            DO I=2,NGI
               L1(I)=1./6.
               L2(I)=1./6.
               L3(I)=1./6.
               L4(I)=1./6.
               WEIGHT(I)=9./20.
            END DO
            L1(2)=0.5
            L2(3)=0.5
            L3(4)=0.5
            L4(5)=0.5
         ENDIF
         !
         IF(NGI.EQ.11) THEN
            ! Degree of precision is 4
            ALPHA=(1.+SQRT(5./14.))/4.0
            BETA =(1.-SQRT(5./14.))/4.0
            I=1
            L1(I)=0.25
            L2(I)=0.25
            L3(I)=0.25
            WEIGHT(I)=-6.*74.0/5625.0
            DO I=2,5
               L1(I)=1./14.
               L2(I)=1./14.
               L3(I)=1./14.
               WEIGHT(I)=6.*343./45000.
            END DO
            L1(2)=11./14.
            L2(3)=11./14.
            L3(4)=11./14.
            DO I=6,11
               L1(I)=ALPHA
               L2(I)=ALPHA
               L3(I)=ALPHA
               WEIGHT(I)=6.*56.0/2250.0
            END DO
            L3(6)=BETA
            L2(7)=BETA
            L2(8)=BETA
            L3(8)=BETA
            L1(9)=BETA
            L1(10)=BETA
            L3(10)=BETA
            L1(11)=BETA
            L2(11)=BETA
            ! ENDOF IF(NGI.EQ.11) THEN...
         ENDIF
         DO I=1,NGI
            L4(I)=1.0-L1(I)-L2(I)-L3(I)
         END DO
         ! Now multiply by 1/6. to get weigts correct...
         DO I=1,NGI
            WEIGHT(I)=WEIGHT(I)/6.
         END DO
         ! ENDOF IF(D3) THEN...
      ENDIF
      !
      IF(.NOT.D3) THEN
         ! 2-D TRAINGULAR ELEMENTS...
         IF(NGI.EQ.1) THEN
            ! LINEAR
            I=1
            L1(I)=1./3.
            L2(I)=1./3.
            WEIGHT(I)=1.0
         ENDIF
         !
         IF(NGI.EQ.3) THEN
            ! QUADRASTIC
            DO I=1,NGI
               L1(I)=0.5
               L2(I)=0.5
               WEIGHT(I)=1.0/3.0
            END DO
            L1(2)=0.0
            L2(3)=0.0
         ENDIF
         !
         IF(NGI.EQ.4) THEN
            ! CUBIC
            I=1
            L1(I)=1./3.
            L2(I)=1./3.
            WEIGHT(I)=-27./48.
            DO I=2,NGI
               L1(I)=0.2
               L2(I)=0.2
               WEIGHT(I)=25./48.
            END DO
            L1(1)=0.6
            L2(2)=0.6
         ENDIF
         !
         IF(NGI.EQ.7) THEN
            ! QUNTIC
            ALPHA1=0.0597158717
            BETA1 =0.4701420641
            ALPHA2=0.7974269853
            BETA2 =0.1012865073
            I=1
            L1(I)=1./3.
            L2(I)=1./3.
            WEIGHT(I)=0.225
            DO I=2,4
               L1(I)=BETA1
               L2(I)=BETA1
               WEIGHT(I)=0.1323941527
            END DO
            L1(2)=ALPHA1
            L2(4)=ALPHA1
            DO I=5,7
               L1(I)=BETA2
               L2(I)=BETA2
               WEIGHT(I)=0.1259391805
            END DO
            L1(5)=ALPHA2
            L2(6)=ALPHA2
            ! ENDOF IF(NGI.EQ.7) THEN...
         ENDIF

         IF(NGI.EQ.14) THEN
            ! 5th order quadrature set...
            L1(1) = 6.943184420297371E-002
            L1(2) = 6.943184420297371E-002
            L1(3) = 6.943184420297371E-002
            L1(4) = 6.943184420297371E-002
            L1(5) = 6.943184420297371E-002
            L1(6) = 0.330009478207572
            L1(7) = 0.330009478207572
            L1(8) = 0.330009478207572
            L1(9) = 0.330009478207572
            L1(10) = 0.669990521792428
            L1(11) = 0.669990521792428
            L1(12) = 0.669990521792428
            L1(13) = 0.930568155797026
            L1(14) = 0.930568155797026
            ! local coord 1:
            L2(1) = 4.365302387072518E-002
            L2(2) = 0.214742881469342
            L2(3) = 0.465284077898513
            L2(4) = 0.715825274327684
            L2(5) = 0.886915131926301
            L2(6) = 4.651867752656094E-002
            L2(7) = 0.221103222500738
            L2(8) = 0.448887299291690
            L2(9) = 0.623471844265867
            L2(10) = 3.719261778493340E-002
            L2(11) = 0.165004739103786
            L2(12) = 0.292816860422638
            L2(13) = 1.467267513102734E-002
            L2(14) = 5.475916907194637E-002
            ! local coord 2:
            WEIGHT(1) = 1.917346464706755E-002
            WEIGHT(2) = 3.873334126144628E-002
            WEIGHT(3) = 4.603770904527855E-002
            WEIGHT(4) = 3.873334126144628E-002
            WEIGHT(5) = 1.917346464706755E-002
            WEIGHT(6) = 3.799714764789616E-002
            WEIGHT(7) = 7.123562049953998E-002
            WEIGHT(8) = 7.123562049953998E-002
            WEIGHT(9) = 3.799714764789616E-002
            WEIGHT(10) = 2.989084475992800E-002
            WEIGHT(11) = 4.782535161588505E-002
            WEIGHT(12) = 2.989084475992800E-002
            WEIGHT(13) = 6.038050853208200E-003
            WEIGHT(14) = 6.038050853208200E-003
            rsum=SUM(WEIGHT(1:NGI))
            WEIGHT(1:NGI)=WEIGHT(1:NGI)/RSUM
            ! ENDOF IF(NGI.EQ.14) THEN...
         ENDIF

         !
         DO I=1,NGI
            L3(I)=1.0-L1(I)-L2(I)
         END DO
         ! ENDOF IF(.NOT.D3) THEN...
      ENDIF
      !
      RETURN 
    END subroutine TRIQUAold

!!$ ===============
!!$ ===============


    SUBROUTINE SPECTR(NGI,NLOC,MLOC,&
         &      M,WEIGHT,N,NLX,NLY,NLZ,D3,D2, IPOLY,IQADRA  )
      IMPLICIT NONE
      INTEGER , intent(in) :: NGI, NLOC, MLOC
      INTEGER , intent(inout) :: IPOLY,IQADRA
      REAL, dimension(:,:), intent(inout) ::  M, N, NLX, NLY, NLZ
      real, dimension(:), intent(inout) :: WEIGHT
      !Local variables
      REAL :: RGPTWE
      REAL , dimension (30) :: WEIT,NODPOS,QUAPOS
      INTEGER :: GPOI
      LOGICAL :: DIFF,NDIFF,D3,D2
      INTEGER :: NDGI,NDNOD,NMDNOD,IGR,IGQ,IGP,KNOD,JNOD,INOD,ILOC
      REAL :: LXGP,LYGP,LZGP
      ! This subroutine defines a spectal element.
      ! IPOLY defines the element type and IQADRA the quadrature.
      ! In 2-D the spectral local node numbering is as..
      ! 7 8 9
      ! 4 5 6
      ! 1 2 3
      ! For 3-D...
      ! lz=-1
      ! 3 4
      ! 1 2
      ! and for lz=1
      ! 7 8
      ! 5 6

      !ewrite(3,*)'inside SPECTR IPOLY,IQADRA', IPOLY,IQADRA
      !
      DIFF=.TRUE.
      NDIFF=.FALSE.
      IF(D3) THEN
         NDGI  =INT((NGI**(1./3.))+0.1)
         NDNOD =INT((NLOC**(1./3.))+0.1)
         NMDNOD=INT((MLOC**(1./3.))+0.1)
         !
         ! Find the roots of the quadrature points and nodes
         ! also get the weights.
         !ewrite(3,*)'about to go into inside GTROOT IPOLY,IQADRA',IPOLY,IQADRA
         CALL GTROOT(IPOLY,IQADRA,WEIT,NODPOS,QUAPOS,NDGI,NDNOD)
         !ewrite(3,*)'outside GTROOT'
         do  IGR=1,NDGI! Was loop 101
            do  IGQ=1,NDGI! Was loop 101
               do  IGP=1,NDGI! Was loop 101
                  GPOI=IGP + (IGQ-1)*NDGI + (IGR-1)*NDGI*NDGI
                  !
                  !           WEIGHT(GPOI)
                  !     &        =RGPTWE(IGP,NDGI,.TRUE.)*RGPTWE(IGQ,NDGI,.TRUE.)
                  !     &        *RGPTWE(IGR,NDGI,.TRUE.)
                  WEIGHT(GPOI)=WEIT(IGP)*WEIT(IGQ)*WEIT(IGR)
                  !
                  LXGP=QUAPOS(IGP)
                  LYGP=QUAPOS(IGQ)
                  LZGP=QUAPOS(IGR)
                  ! NB If TRUE in function RGPTWE then return the Gauss-pt weight
                  ! else return the Gauss-pt.
                  !
                  do  KNOD=1,NDNOD! Was loop 20
                     do  JNOD=1,NDNOD! Was loop 20
                        do  INOD=1,NDNOD! Was loop 20
                           ILOC=INOD + (JNOD-1)*NDNOD + (KNOD-1)*NDNOD*NDNOD
                           !
                           N(ILOC,GPOI) = &
                                SPECFU(NDIFF,LXGP,INOD,NDNOD,IPOLY,NODPOS)&
                                *SPECFU(NDIFF,LYGP,JNOD,NDNOD,IPOLY,NODPOS)&
                                *SPECFU(NDIFF,LZGP,KNOD,NDNOD,IPOLY,NODPOS)
                           !
                           NLX(ILOC,GPOI) = &
                                SPECFU(DIFF,LXGP, INOD,NDNOD,IPOLY,NODPOS)&
                                *SPECFU(NDIFF,LYGP,JNOD,NDNOD,IPOLY,NODPOS)&
                                *SPECFU(NDIFF,LZGP,KNOD,NDNOD,IPOLY,NODPOS)
                           !
                           NLY(ILOC,GPOI) = &
                                SPECFU(NDIFF,LXGP,INOD,NDNOD,IPOLY,NODPOS)&
                                *SPECFU(DIFF, LYGP,JNOD,NDNOD,IPOLY,NODPOS)&
                                *SPECFU(NDIFF,LZGP,KNOD,NDNOD,IPOLY,NODPOS)
                           !
                           NLZ(ILOC,GPOI) = &
                                SPECFU(NDIFF,LXGP,INOD,NDNOD,IPOLY,NODPOS)&
                                *SPECFU(NDIFF,LYGP,JNOD,NDNOD,IPOLY,NODPOS)&
                                *SPECFU(DIFF,LZGP, KNOD,NDNOD,IPOLY,NODPOS)
                           !
                        end do ! Was loop 20
                     end do ! Was loop 20
                  end do ! Was loop 20
               end do ! Was loop 101
            end do ! Was loop 101
         end do ! Was loop 101
         !
         !
         ! Find the roots of the quadrature points and nodes
         ! also get the weights.
         !ewrite(3,*)'2about to go into inside GTROOT IPOLY,IQADRA',IPOLY,IQADRA
         CALL GTROOT(IPOLY,IQADRA,WEIT,NODPOS,QUAPOS,NDGI,NMDNOD)
         !ewrite(3,*)'2out of GTROOT'
         do  IGR=1,NDGI! Was loop 102
            do  IGQ=1,NDGI! Was loop 102
               do  IGP=1,NDGI! Was loop 102
                  GPOI=IGP + (IGQ-1)*NDGI + (IGR-1)*NDGI*NDGI
                  !
                  LXGP=QUAPOS(IGP)
                  LYGP=QUAPOS(IGQ)
                  LZGP=QUAPOS(IGR)

                  do  KNOD=1,NMDNOD! Was loop 30
                     do  JNOD=1,NMDNOD! Was loop 30
                        do  INOD=1,NMDNOD! Was loop 30
                           ILOC=INOD + (JNOD-1)*NMDNOD + (KNOD-1)*NMDNOD*NMDNOD
                           !
                           M(ILOC,GPOI) = &
                                SPECFU(NDIFF,LXGP,INOD,NMDNOD,IPOLY,NODPOS)  &
                                *SPECFU(NDIFF,LYGP,JNOD,NMDNOD,IPOLY,NODPOS) &
                                *SPECFU(NDIFF,LZGP,KNOD,NMDNOD,IPOLY,NODPOS)
                           !
                        end do ! Was loop 30
                     end do ! Was loop 30
                  end do ! Was loop 30
               end do ! Was loop 102
            end do ! Was loop 102
         end do ! Was loop 102
      ENDIF
      !
      IF(D2) THEN
         NDGI  =INT((NGI**(1./2.))+0.1)
         NDNOD =INT((NLOC**(1./2.))+0.1)
         NMDNOD=INT((MLOC**(1./2.))+0.1)
         !
         ! Find the roots of the quadrature points and nodes
         ! also get the weights.
         CALL GTROOT(IPOLY,IQADRA,WEIT,NODPOS,QUAPOS,NDGI,NDNOD)
         do  IGQ=1,NDGI! Was loop 10
            do  IGP=1,NDGI! Was loop 10
               GPOI=IGP + (IGQ-1)*NDGI
               !
               WEIGHT(GPOI)=WEIT(IGP)*WEIT(IGQ)
               !
               LXGP=QUAPOS(IGP)
               LYGP=QUAPOS(IGQ)
               ! NB If TRUE in function RGPTWE then return the Gauss-pt weight
               ! else return the Gauss-pt.
               !
               do  JNOD=1,NDNOD! Was loop 120
                  do  INOD=1,NDNOD! Was loop 120
                     ILOC=INOD + (JNOD-1)*NDNOD
                     !
                     N(ILOC,GPOI) = &
                          SPECFU(NDIFF,LXGP,INOD,NDNOD,IPOLY,NODPOS) &
                          *SPECFU(NDIFF,LYGP,JNOD,NDNOD,IPOLY,NODPOS)
                     !
                     NLX(ILOC,GPOI) = &
                          SPECFU(DIFF, LXGP,INOD,NDNOD,IPOLY,NODPOS) &
                          *SPECFU(NDIFF,LYGP,JNOD,NDNOD,IPOLY,NODPOS)
                     !
                     NLY(ILOC,GPOI) = &
                          SPECFU(NDIFF,LXGP,INOD,NDNOD,IPOLY,NODPOS) &
                          *SPECFU(DIFF, LYGP,JNOD,NDNOD,IPOLY,NODPOS)
                     !
                  end do ! Was loop 120
               end do ! Was loop 120
            end do ! Was loop 10
         end do ! Was loop 10
         !
         ! Find the roots of the quadrature points and nodes
         ! also get the weights.
         CALL GTROOT(IPOLY,IQADRA,WEIT,NODPOS,QUAPOS,NDGI,NMDNOD)
         do  IGQ=1,NDGI! Was loop 11
            do  IGP=1,NDGI! Was loop 11
               GPOI=IGP + (IGQ-1)*NDGI
               LXGP=QUAPOS(IGP)
               LYGP=QUAPOS(IGQ)
               do  JNOD=1,NMDNOD! Was loop 130
                  do  INOD=1,NMDNOD! Was loop 130
                     ILOC=INOD + (JNOD-1)*NMDNOD
                     !
                     M(ILOC,GPOI) = &
                          SPECFU(NDIFF,LXGP,INOD,NMDNOD,IPOLY,NODPOS) &
                          *SPECFU(NDIFF,LYGP,JNOD,NMDNOD,IPOLY,NODPOS)
                     !
                  end do ! Was loop 130
               end do ! Was loop 130
               !
            end do ! Was loop 11
         end do ! Was loop 11
      ENDIF
      !
      !
      IF((.NOT.D2).AND.(.NOT.D3)) THEN
         NDGI  =NGI
         NDNOD =NLOC
         NMDNOD=MLOC
         !
         ! Find the roots of the quadrature points and nodes
         ! also get the weights.
         CALL GTROOT(IPOLY,IQADRA,WEIT,NODPOS,QUAPOS,NDGI,NDNOD)
         !ewrite(3,*)'NDGI,NDNOD,NLOC:',NDGI,NDNOD,NLOC
         !ewrite(3,*)'WEIT(1:ndgi):',WEIT(1:ndgi)
         !ewrite(3,*)'NODPOS(1:ndnod):',NODPOS(1:ndnod)
         !ewrite(3,*)'QUAPOS(1:ndgi):',QUAPOS(1:ndgi)
         do  IGP=1,NDGI! Was loop 1000
            GPOI=IGP
            !
            WEIGHT(GPOI)=WEIT(IGP)
            !
            LXGP=QUAPOS(IGP)
            ! NB If TRUE in function RGPTWE then return the Gauss-pt weight
            ! else return the Gauss-pt.
            !
            do  INOD=1,NDNOD! Was loop 12000
               ILOC=INOD
               !
               N(ILOC,GPOI) = &
                    SPECFU(NDIFF,LXGP,INOD,NDNOD,IPOLY,NODPOS)
               !
               NLX(ILOC,GPOI) = &
                    SPECFU(DIFF, LXGP,INOD,NDNOD,IPOLY,NODPOS)
               !ewrite(3,*)'ILOC,GPOI,N(ILOC,GPOI),NLX(ILOC,GPOI):', &
               !         ILOC,GPOI,N(ILOC,GPOI),NLX(ILOC,GPOI)
               !
            end do ! Was loop 12000
         end do ! Was loop 1000
         !ewrite(3,*)'n WEIGHT:',WEIGHT
         !
         ! Find the roots of the quadrature points and nodes
         ! also get the weights.
         !ewrite(3,*)'this is for m which we dont care about:'
         CALL GTROOT(IPOLY,IQADRA,WEIT,NODPOS,QUAPOS,NDGI,NMDNOD)
         do  IGP=1,NDGI! Was loop 1100
            GPOI=IGP
            LXGP=QUAPOS(IGP)
            do  INOD=1,NMDNOD! Was loop 13000
               ILOC=INOD
               !
               M(ILOC,GPOI) = &
                    SPECFU(NDIFF,LXGP,INOD,NMDNOD,IPOLY,NODPOS)
               !
            end do ! Was loop 13000
            !
         end do ! Was loop 1100
         !ewrite(3,*)'...finished this is for m which we dont care about:'
      ENDIF
    END SUBROUTINE SPECTR


!!$ ===============
!!$ ===============

  end module TriangleTet_ShapeFunctions
