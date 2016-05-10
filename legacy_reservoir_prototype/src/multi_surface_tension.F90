  
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
  

module multi_surface_tension
    use fldebug
    use elements
    use futils, only: int2str
    use field_options
    use state_module
    use Copy_Outof_State
    use spud
    use parallel_tools
    use multi_data_types
    use shape_functions
    use shape_functions_NDim
    use shape_functions_prototype
    use fields
    use cv_advection, only : dgsimplnorm
    use matrix_operations, only : smlinngot
    use multi_tools, only: CALC_FACE_ELE

    implicit none

contains

 SUBROUTINE CALCULATE_SURFACE_TENSION_NEW( state, packed_state, Mdims, Mspars, ndgln, Mdisopt, &
     PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, IPLIKE_GRAD_SOU)

     IMPLICIT NONE

     real, dimension( :, :, : ), intent( inout ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD
     integer, intent( inout ) :: IPLIKE_GRAD_SOU

     type(state_type), dimension( : ), intent( inout ) :: state
     type(state_type), intent( inout ) :: packed_state
     type(multi_dimensions), intent(in) :: Mdims
     type(multi_sparsities), intent(in) :: Mspars
     type(multi_ndgln), intent(in) :: ndgln
     type(multi_discretization_opts), intent(in) :: Mdisopt
     
     !Local variables
     real, dimension( : ), allocatable :: X, Y, Z

     integer :: iphase, icomp
     real :: coefficient, angle
     logical :: surface_tension, use_pressure_force, use_smoothing

     type( vector_field ), pointer :: x_all
     type( tensor_field ), pointer :: MFC_s


     allocate( X(  Mdims%x_nonods ) ) ; X = 0.0
     allocate( Y(  Mdims%x_nonods ) ) ; Y = 0.0
     allocate( Z(  Mdims%x_nonods ) ) ; Z = 0.0


     x_all => extract_vector_field( packed_state, "PressureCoordinate" )
     x = x_all % val( 1, : )
     if (Mdims%ndim >=2 ) y = x_all % val( 2, : )
     if (Mdims%ndim >=3 ) z = x_all % val( 3, : )


     ! Initialise...
     IPLIKE_GRAD_SOU = 0
     PLIKE_GRAD_SOU_COEF = 0.0
     PLIKE_GRAD_SOU_GRAD = 0.0


     do icomp = 1, Mdims%ncomp

         surface_tension = have_option( '/material_phase[' // int2str( Mdims%nphase - 1 + icomp ) // &
             ']/is_multiphase_component/surface_tension' )

         if ( surface_tension ) then

             MFC_s  => extract_tensor_field( packed_state, "PackedComponentMassFraction" )

             ewrite(3,*) 'Calculating surface tension for component ', icomp

             call get_option( '/material_phase[' // int2str( Mdims%nphase - 1 + icomp ) // &
                 ']/is_multiphase_component/surface_tension/coefficient', coefficient )

             use_smoothing = have_option( '/material_phase[' // int2str( Mdims%nphase - 1 + icomp ) // &
                 ']/is_multiphase_component/surface_tension/smooth' )

             call get_option( '/material_phase[' // int2str( Mdims%nphase - 1 + icomp ) // &
                 ']/is_multiphase_component/surface_tension/angle', angle, default = -1.0 )

             USE_PRESSURE_FORCE = .TRUE.

             if ( USE_PRESSURE_FORCE ) then
                 IPLIKE_GRAD_SOU = 1
             else
                 IPLIKE_GRAD_SOU = 0
             end if

             do iphase = 1, Mdims%nphase

                 CALL SURFACE_TENSION_WRAPPER_NEW( state, packed_state, &
                     PLIKE_GRAD_SOU_COEF( icomp, iphase, :), PLIKE_GRAD_SOU_GRAD( icomp, iphase, :), &
                     !PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, &
                     COEFFICIENT, ANGLE, &
                     MFC_s%val( icomp, iphase, :), &
                     Mdims, Mspars, ndgln, Mdisopt )

             end do

             if ( .not.USE_PRESSURE_FORCE ) then

                ewrite(3,*) 'Error, should use as a pressure force'

             end if


         end if

     end do

     ewrite(3,*) 'Leaving CALCULATE_SURFACE_TENSION_NEW'

     RETURN
 END SUBROUTINE CALCULATE_SURFACE_TENSION_NEW

 SUBROUTINE SURFACE_TENSION_WRAPPER_NEW( state, packed_state, &
     PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, &
     SUF_TENSION_COEF, ANGLE, VOLUME_FRAC, &
     Mdims, Mspars, ndgln, Mdisopt )
     ! Calculate the surface tension force as a pressure force term:
     ! PLIKE_GRAD_SOU_COEF and PLIKE_GRAD_SOU_GRAD 
     ! for a given DevFuns%VOLUME fraction field VOLUME_FRAC
     ! SUF_TENSION_COEF is the surface tension coefficient.
     !use shape_functions
     !use matrix_operations
     ! Inputs/Outputs
     IMPLICIT NONE
     type(state_type), dimension( : ), intent( inout ) :: state
     type(state_type), intent( inout ) :: packed_state
     type(multi_dimensions), intent(in) :: Mdims
     type(multi_sparsities), intent(in) :: Mspars
     type(multi_ndgln), intent(in) :: ndgln
     type(multi_discretization_opts), intent(in) :: Mdisopt
     type(multi_gi_dimensions) :: CV_GIdims
     type(multi_shape_funs) :: CV_funs
     REAL, DIMENSION( :), intent( inout ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD
     REAL, intent( in ) ::  SUF_TENSION_COEF, ANGLE
     REAL, DIMENSION( : ), intent( in ) :: VOLUME_FRAC
        ! Local variables
     
        INTEGER, DIMENSION( : ), allocatable :: CV_OTHER_LOC, CV_SLOC2LOC
        INTEGER, DIMENSION( : , : ), allocatable :: FACE_ELE
        REAL, DIMENSION( : ), allocatable :: MASS_CV, MASS_ELE, CURVATURE
        !        ===> INTEGERS <===
        INTEGER :: ELE, ELE2, GI, SELE, CV_SILOC, CV_ILOC, CV_JLOC, CV_JNOD, CV_NOD
        !        ===>  REALS  <===
        REAL :: HDC, NN, RR
        REAL, PARAMETER :: TOLER=1.0E-10
        !        ===>  LOGICALS  <===
        LOGICAL :: QUAD_OVER_WHOLE_ELE, GOTDEC
        ! Approaches to calculate the curvature using Diffused Interface Approach or Distance Function Approach
        LOGICAL, PARAMETER :: Dif_Int_App = .true.
        LOGICAL, PARAMETER :: DCYL = .FALSE.
        !Derivative of the shape functions
        type(multi_dev_shape_funs) :: Devfuns
!--------------------DISTANCE FUNCTION VARIABLES-----------------------
      REAL, DIMENSION( : , : ), allocatable :: SIGN_FUN_GI, MAT_LOC
      REAL, DIMENSION( : , : , : ), allocatable :: SIGN_FUN_SGI
      REAL, DIMENSION( : ), allocatable ::  DISTANCE_FUN, DISTANCE_FUN1, DISTANCE_FUN_OLD, &
               RHS_CV_SHORT,CV_SOL, XSL, YSL, ZSL, & 
               S_INCOME, SDETWE, N_DOT_SQ,  &
               NORMALIZATION, CURV, ST, DIFF_COEF
      INTEGER :: CV_ILOC2, CV_JLOC2, CV_SJLOC, &
                 DG_NOD, DG_JNOD, DG_JNOD2, DG_INOD2, ITIME, &
                 INTERAT_FACE, IFACE, SGI, SELE2, X_INOD, DG_NONODS, X_INOD2, OUTER_ITS, INNER_ITS
      REAL :: EH, PSI_GI,PSI_SGI, SUR_DT, &
              VLN, SVLN_IN, SVLN_OUT, &
              SAREA, &
              SQRT_RR, SNN, RNN, DNN
      LOGICAL :: GTHALF, LTHALF
      LOGICAL, DIMENSION(:), allocatable :: INTERFACE_ELE, INTERFACE_ELE2
      INTEGER :: IPIV(Mdims%cv_nloc)
!------------------------------------------------------
      type(scalar_field), pointer :: ds, dx, dy, dz, dk, ck
      integer, dimension(:), pointer :: state_nodes
!-----------------------------------------------------------------------
     !Local variables
      real, dimension( : ), allocatable :: X, Y, Z, PSIGI_X, DGI_X, PSISGI_X, NORMX
      real, dimension( : , : ), allocatable :: SOL_DERIV_X, SOL_DERIV_X1, UD, SNORMXN
     
      type( vector_field ), pointer :: x_all
      QUAD_OVER_WHOLE_ELE=.true. 
      ! If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
      call retrieve_ngi( CV_GIdims, Mdims, Mdisopt%cv_ele_type, quad_over_whole_ele )
      call allocate_multi_shape_funs( CV_funs, Mdims, CV_GIdims )
      call cv_fem_shape_funs( CV_funs, Mdims, CV_GIdims, Mdisopt%cv_ele_type, quad_over_whole_ele )
      ALLOCATE( X(  Mdims%x_nonods ) ) ; X = 0.0
      ALLOCATE( Y(  Mdims%x_nonods ) ) ; Y = 0.0
      ALLOCATE( Z(  Mdims%x_nonods ) ) ; Z = 0.0
      ALLOCATE( PSIGI_X( Mdims%ndim ) ) ; PSIGI_X = 0.0
      ALLOCATE( DGI_X( Mdims%ndim ) ) ; DGI_X = 0.0
      ALLOCATE( UD( Mdims%ndim, CV_GIdims%cv_ngi )); UD = 0.0
      ALLOCATE( PSISGI_X( Mdims%ndim ) ) ; PSISGI_X = 0.0
      ALLOCATE( NORMX( 3 ) ) ; NORMX = 0.0
      ALLOCATE( SNORMXN( Mdims%ndim, CV_GIdims%sbcvngi )); SNORMXN = 0.0
      ALLOCATE( CV_OTHER_LOC( Mdims%cv_nloc ))    
      ALLOCATE( CV_SLOC2LOC( Mdims%cv_snloc ))

      x_all => extract_vector_field( packed_state, "PressureCoordinate" )
      x = x_all % val( 1, : )
      if (Mdims%ndim >=2 ) y = x_all % val( 2, : )
      if (Mdims%ndim >=3 ) z = x_all % val( 3, : )

      ds => extract_scalar_field( state(1), 'ds'  )
      dx => extract_scalar_field( state(1), 'dx'  )
      dy => extract_scalar_field( state(1), 'dy'  )
      dz => extract_scalar_field( state(1), 'dz'  )
      dk => extract_scalar_field( state(1), 'dk'  )
      ck => extract_scalar_field( state(1), 'ck'  )
 
      ALLOCATE( FACE_ELE( CV_GIdims%nface, Mdims%totele ) ) ; FACE_ELE = 0
      CALL CALC_FACE_ELE( FACE_ELE, Mdims%totele, Mdims%stotel, CV_GIdims%nface, &
           Mspars%ELE%fin, Mspars%ELE%col, Mdims%cv_nloc, Mdims%cv_snloc, Mdims%cv_nonods, ndgln%cv, ndgln%suf_cv, &
           CV_funs%cv_sloclist, Mdims%x_nloc, ndgln%x )
      ALLOCATE( MASS_ELE( Mdims%totele )); MASS_ELE=0.0
      

      call allocate_multi_dev_shape_funs(CV_funs, Devfuns)
      DO ELE=1,Mdims%totele
        ! Calculate DevFuns%DETWEI,DevFuns%RA,NX,NY,NZ for element ELE
        call DETNLXR_PLUS_U(ELE, x_all % val, ndgln%x, CV_funs%cvweight, &
              CV_funs%cvfen, CV_funs%cvfenlx_all, CV_funs%ufenlx_all, Devfuns)
         MASS_ELE( ELE ) = DevFuns%VOLUME
      END DO
      ALLOCATE( MASS_CV(  Mdims%cv_nonods ) ) ; MASS_CV=0.0
 
      DO ELE=1,Mdims%totele
         DO CV_ILOC=1,Mdims%cv_nloc
            CV_NOD=ndgln%cv((ELE-1)*Mdims%cv_nloc+CV_ILOC)
            MASS_CV(CV_NOD) = MASS_CV(CV_NOD) + MASS_ELE(ELE) 
          END DO
      END DO
!---------------- Calculate the distance function -------------------
      ALLOCATE(CURVATURE( Mdims%cv_nonods))
      DG_NONODS=Mdims%cv_nloc*Mdims%totele
      ALLOCATE( INTERFACE_ELE( Mdims%totele ))
      ALLOCATE( INTERFACE_ELE2( Mdims%totele ))
      ALLOCATE( DISTANCE_FUN( DG_NONODS ))
      ALLOCATE( DISTANCE_FUN_OLD( DG_NONODS ))
      ALLOCATE( SOL_DERIV_X( Mdims%ndim, DG_NONODS ))
      ALLOCATE( SIGN_FUN_GI( Mdims%totele, CV_GIdims%cv_ngi ))
      ALLOCATE( SIGN_FUN_SGI( Mdims%totele, CV_GIdims%nface, CV_GIdims%sbcvngi ))
      ALLOCATE( XSL(Mdims%cv_snloc) )
      ALLOCATE( YSL(Mdims%cv_snloc) )
      ALLOCATE( ZSL(Mdims%cv_snloc) )
      ALLOCATE( S_INCOME(CV_GIdims%sbcvngi) )
      ALLOCATE( SDETWE(CV_GIdims%sbcvngi) )
      ALLOCATE( N_DOT_SQ(CV_GIdims%sbcvngi) )
      ALLOCATE( CURV( DG_NONODS ))
      ALLOCATE( NORMALIZATION( DG_NONODS ))
      ALLOCATE( ST( CV_GIdims%cv_ngi ))
      ALLOCATE( SOL_DERIV_X1( Mdims%ndim,  Mdims%cv_nonods ))
      ALLOCATE( DISTANCE_FUN1(  Mdims%cv_nonods ))
      ALLOCATE( DIFF_COEF( CV_GIdims%cv_ngi ))
! Convert force to velocity space...
      ALLOCATE(MAT_LOC(Mdims%cv_nloc,Mdims%cv_nloc))
      ALLOCATE(RHS_CV_SHORT(Mdims%cv_nloc) )
      ALLOCATE(CV_SOL(Mdims%cv_nloc))
! Determine which elements to integrate over
! we just need to integrate over the elements near the boundary...
      INTERFACE_ELE=.FALSE.
      DO ELE=1,Mdims%totele
         GTHALF=.FALSE.
         LTHALF=.FALSE.
         DO CV_ILOC=1,Mdims%cv_nloc
            CV_NOD=ndgln%cv((ELE-1)*Mdims%cv_nloc+CV_ILOC) 
            IF(VOLUME_FRAC(CV_NOD).GT.0.5) THEN
               GTHALF=.TRUE.
            ELSE
               LTHALF=.TRUE.
            ENDIF
         END DO
         IF(GTHALF.AND.LTHALF) THEN 
           INTERFACE_ELE(ELE)=.TRUE.
         ELSE
           INTERFACE_ELE(ELE)=.FALSE.
         ENDIF
      END DO ! ELE
      DO INTERAT_FACE=1,3  ! have a distance 4 (in terms of elemnts) away from interface...
         INTERFACE_ELE2=INTERFACE_ELE
         DO ELE=1,Mdims%totele
            IF(INTERFACE_ELE2(ELE)) THEN
               DO IFACE = 1, CV_GIdims%nface
                  ELE2  = FACE_ELE( IFACE, ELE )
                  IF(ELE2.GT.0) INTERFACE_ELE(ELE2)=.TRUE.
               END DO
            ENDIF
         END DO ! ELE
      END DO
! Calculate the maximum lengthscale near the interface
      HDC=0.0
      DO ELE=1,Mdims%totele
         IF(INTERFACE_ELE(ELE)) THEN
         DO CV_ILOC=2,Mdims%cv_nloc
            X_INOD = ndgln%x(( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
            X_INOD2 = ndgln%x(( ELE - 1 ) * Mdims%cv_nloc + 1 )
            !HDC=MAX(HDC, SQRT( (X(X_INOD)-X(X_INOD2))**2+(Y(X_INOD)-Y(X_INOD2))**2 &
            !              +(Z(X_INOD)-Z(X_INOD2))**2  ))
            HDC = MAX(HDC, SQRT( SUM( (X_ALL%val(1:Mdims%ndim,X_INOD)-X_ALL%val(1:Mdims%ndim,X_INOD2))**2) ))
         END DO
         END IF
      END DO
      
      if (IsParallel()) call allmax(HDC)
! Initialise the distance function
      DO ELE=1,Mdims%totele
         DO CV_ILOC=1,Mdims%cv_nloc
            CV_NOD=ndgln%cv((ELE-1)*Mdims%cv_nloc+CV_ILOC)
            DG_NOD=(ELE-1)*Mdims%cv_nloc+CV_ILOC
            IF (Dif_Int_App) THEN
                DISTANCE_FUN(DG_NOD)=VOLUME_FRAC(CV_NOD)*HDC
            ELSE
                DISTANCE_FUN(DG_NOD)=2.*(VOLUME_FRAC(CV_NOD)-0.5)*HDC
            END IF
         END DO
      END DO
! Sign function S(PSI)
      IF (Dif_Int_App) THEN
          SIGN_FUN_GI = 0.0
      ELSE
         EH=0.02*HDC ! assume h=5.*DT at present
         DO ELE=1,Mdims%totele
            DO GI=1,CV_GIdims%cv_ngi
               PSI_GI=0.0
               DO CV_ILOC=1,Mdims%cv_nloc
                  DG_NOD=(ELE-1)*Mdims%cv_nloc+CV_ILOC
                  PSI_GI=PSI_GI + CV_funs%CVFEN( CV_ILOC, GI ) * DISTANCE_FUN(DG_NOD)
               END DO        
               SIGN_FUN_GI(ELE,GI)=PSI_GI/SQRT(PSI_GI**2+EH**2)
            END DO
         END DO
  
         DO ELE=1,Mdims%totele
            Between_Elements_And_Boundary0: DO IFACE = 1, CV_GIdims%nface
            ! The surface nodes on element face IFACE. 
            CV_SLOC2LOC( : ) = CV_funs%cv_sloclist( IFACE, : )
             
            DO SGI=1,CV_GIdims%sbcvngi
               PSI_SGI=0.0
               DO CV_SILOC=1,Mdims%cv_snloc
                  CV_ILOC=CV_SLOC2LOC( CV_SILOC )
                  DG_NOD=(ELE-1)*Mdims%cv_nloc+CV_ILOC
                  PSI_SGI=PSI_SGI + CV_funs%sbcvfen( CV_SILOC, SGI ) * DISTANCE_FUN(DG_NOD)
               END DO
               SIGN_FUN_SGI(ELE,IFACE,SGI)=PSI_SGI/SQRT(PSI_SGI**2+EH**2)
               !IF (ABS(PSI_SGI).GT.0.5*HDC) SIGN_FUN_SGI(ELE,IFACE,SGI)=PSI_SGI/SQRT(PSI_SGI**2)
            END DO
            END DO Between_Elements_And_Boundary0
         END DO ! ELE
      
      END IF ! Dif_Int_App
! Calculate the minimum lengthscale near the interface
      HDC=1.0E+6
      DO ELE=1,Mdims%totele
         IF(INTERFACE_ELE(ELE)) THEN
         DO CV_ILOC=2,Mdims%cv_nloc
            X_INOD = ndgln%x(( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
            X_INOD2 = ndgln%x(( ELE - 1 ) * Mdims%cv_nloc + 1 )
            !HDC=MIN(HDC, SQRT( (X(X_INOD)-X(X_INOD2))**2+(Y(X_INOD)-Y(X_INOD2))**2 &
            !              +(Z(X_INOD)-Z(X_INOD2))**2  ))
            HDC = MIN(HDC, SQRT( SUM( (X_ALL%val(1:Mdims%ndim,X_INOD)-X_ALL%val(1:Mdims%ndim,X_INOD2))**2) ))
         END DO
         END IF
      END DO
      if (IsParallel()) call allmin(HDC)
      DISTANCE_FUN_OLD=DISTANCE_FUN
      INTERFACE_ELE=.TRUE. !comment out eventually
      IF (Dif_Int_App) THEN
         SUR_DT = 1.0
      ELSE
         SUR_DT = 0.1*HDC ! assume Dtau=0.1.* H
      END IF
! Start to calculate the diffused interface DevFuns%VOLUME fraction or distance function
      DO ITIME=1,5
         DO OUTER_ITS=1,1
            SOL_DERIV_X=0.
!Calculate Laplacian of the diffused interface DevFuns%VOLUME fraction or distance function
            DO ELE=1,Mdims%totele ! ELE loop 1
               IF(INTERFACE_ELE(ELE)) THEN
                call DETNLXR_PLUS_U(ELE, x_all % val, ndgln%x, CV_funs%cvweight, &
                       CV_funs%cvfen, CV_funs%cvfenlx_all, CV_funs%ufenlx_all, Devfuns)
 
               CALL LOC_1ST_DERIV_XYZ_DG_DERIV(DISTANCE_FUN, SOL_DERIV_X(1,:), SOL_DERIV_X(2,:), SOL_DERIV_X(3,:), &
                                           Mdims%ndim,  ndgln%cv, &
                                           Mdims%x_nloc, ndgln%x, &
                                           CV_GIdims%cv_ngi, Mdims%cv_nloc, &
                                           CV_funs%CVFEN, DevFuns%CVFENX_ALL(1,:,:), DevFuns%CVFENX_ALL(2,:,:), DevFuns%CVFENX_ALL(3,:,:), &
                                           X, Y, Z, &
                                           CV_GIdims%nface, FACE_ELE, CV_funs%cv_sloclist, Mdims%cv_snloc, &
                                           CV_GIdims%sbcvngi, CV_funs%sbcvfen, CV_funs%sbcvfenslx, CV_funs%sbcvfensly, CV_funs%sbcvfeweigh, &
                                           ELE, DevFuns%DETWEI, 1, angle)
                                  
               END IF
            END DO ! ELE loop 1
            dx%val = SOL_DERIV_X(1,:)
            dy%val = SOL_DERIV_X(2,:)
            IF(Mdims%ndim.GE.3) dz%val = SOL_DERIV_X(3,:)
            call halo_update(dx)
            call halo_update(dy)
            IF(Mdims%ndim.GE.3) call halo_update(dz)
            SOL_DERIV_X(1,:) = dx%val 
            SOL_DERIV_X(2,:) = dy%val 
            IF(Mdims%ndim.GE.3) SOL_DERIV_X(3,:) = dz%val 
       
!   Distribute DG DERIV to the CV nodes...
            SOL_DERIV_X1=0.0
            DO ELE=1,Mdims%totele ! ELE loop 2
               DO CV_ILOC=1,Mdims%cv_nloc
                  CV_NOD=ndgln%cv((ELE-1)*Mdims%cv_nloc+CV_ILOC)
                  DG_NOD=(ELE-1)*Mdims%cv_nloc+CV_ILOC
                  SOL_DERIV_X1( :, CV_NOD)=SOL_DERIV_X1( :, CV_NOD)+SOL_DERIV_X( :, DG_NOD) * MASS_ELE(ELE) / MASS_CV(CV_NOD)
               END DO
            END DO ! ELE loop 2
            DO ELE=1,Mdims%totele ! ELE loop 3
               DO CV_ILOC=1,Mdims%cv_nloc
                  CV_NOD=ndgln%cv((ELE-1)*Mdims%cv_nloc+CV_ILOC)
                  DG_NOD=(ELE-1)*Mdims%cv_nloc+CV_ILOC
                  SOL_DERIV_X( :, DG_NOD)=SOL_DERIV_X1( :, CV_NOD)
               END DO
               IF ( ((Mdims%ndim==2).AND.(Mdims%cv_nloc==6)).or.((Mdims%ndim==3).AND.(Mdims%cv_nloc==10)) ) THEN
                  SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+2)=0.5*(SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+1)+SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+3))
                  SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+4)=0.5*(SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+1)+SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+6))
                  SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+5)=0.5*(SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+3)+SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+6))
                  IF ( Mdims%cv_nloc == 10 ) THEN
                     SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+7)=0.5*(SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+1)+SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+10))
                     SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+8)=0.5*(SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+3)+SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+10))
                     SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+9)=0.5*(SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+6)+SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+10))
                  END IF
               END IF
            END DO ! ELE loop 3
            DO ELE=1,Mdims%totele  ! ELE loop 4
               IF(INTERFACE_ELE(ELE)) THEN
                 ! Calculate DevFuns%DETWEI,DevFuns%RA,NX,NY,NZ for element ELE
                 call DETNLXR_PLUS_U(ELE, x_all % val, ndgln%x, CV_funs%cvweight, &
                       CV_funs%cvfen, CV_funs%cvfenlx_all, CV_funs%ufenlx_all, Devfuns)
               DO INNER_ITS=1,1
! Calculate the velocity...
                  DO GI=1,CV_GIdims%cv_ngi
                     PSIGI_X=0.0
                     DGI_X=0.0
                     DO CV_ILOC=1,Mdims%cv_nloc
                        DG_NOD=(ELE-1)*Mdims%cv_nloc+CV_ILOC
                        PSIGI_X(:)=PSIGI_X(:) + CV_funs%CVFEN( CV_ILOC, GI ) * SOL_DERIV_X(:, DG_NOD)
                        DGI_X(:)=DGI_X(:) + DevFuns%CVFENX_ALL(:, CV_ILOC, GI ) * DISTANCE_FUN(DG_NOD)
                     END DO
                     RR = SUM( PSIGI_X(:)**2 )
                     SQRT_RR=SQRT(RR)
                     UD(:, GI)=SIGN_FUN_GI(ELE,GI)*PSIGI_X(:)/MAX(TOLER, SQRT_RR)
                     IF (Dif_Int_App) THEN
                        DIFF_COEF(GI)=(1.0*HDC)**2!MIN(1.0, DIFF_COEF(GI))
                     ELSE
                        DIFF_COEF(GI)=0.0
                     END IF
                  END DO
 
                  MAT_LOC=0.0
                  RHS_CV_SHORT=0.0
                  DO CV_ILOC=1,Mdims%cv_nloc
                     DO CV_JLOC=1,Mdims%cv_nloc
                        CV_JNOD=ndgln%cv((ELE-1)*Mdims%cv_nloc+CV_JLOC)
                        DG_JNOD=(ELE-1)*Mdims%cv_nloc+CV_JLOC 
                        NN=0.0
                        VLN=0.0
                        DNN=0.0
                        SNN=0.0
                        DO GI=1,CV_GIdims%cv_ngi
                           RNN = CV_funs%CVFEN( CV_ILOC, GI ) * CV_funs%CVFEN( CV_JLOC, GI ) * DevFuns%DETWEI(GI)
                           NN = NN + RNN
                           VLN = VLN+ CV_funs%CVFEN( CV_ILOC, GI ) &
                           * SUM( UD(:,GI)*DevFuns%CVFENX_ALL(:, CV_JLOC, GI ) )* DevFuns%DETWEI(GI)
                           DNN = DNN + SUM( DevFuns%CVFENX_ALL(:, CV_ILOC, GI )*DevFuns%CVFENX_ALL(:, CV_JLOC, GI ) )* DevFuns%DETWEI(GI)*DIFF_COEF(GI)
                           SNN = SNN + SIGN_FUN_GI(ELE,GI)*RNN
                        END DO
                        ! BACK EULER TIME STEPPING
                        MAT_LOC(CV_ILOC,CV_JLOC)=MAT_LOC(CV_ILOC,CV_JLOC)+NN/SUR_DT + VLN + DNN
                        RHS_CV_SHORT(CV_ILOC)=RHS_CV_SHORT(CV_ILOC) &
                                   +(NN/SUR_DT)*DISTANCE_FUN_OLD(DG_JNOD) + SNN
                     END DO
                  END DO
 
                  Between_Elements_And_Boundary: DO IFACE = 1, CV_GIdims%nface
                  ELE2  = FACE_ELE( IFACE, ELE )
                  SELE2 = MAX( 0, - ELE2 )
                  SELE  = SELE2
                  ELE2  = MAX( 0, + ELE2 )
                  IF (ELE2.GT.0) THEN
                ! The surface nodes on element face IFACE. 
                     CV_SLOC2LOC( : ) = CV_funs%cv_sloclist( IFACE, : )
                     DO CV_SILOC = 1, Mdims%cv_snloc
                        CV_ILOC = CV_SLOC2LOC( CV_SILOC )
                        X_INOD = ndgln%x(( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
                        DO CV_ILOC2 = 1, Mdims%cv_nloc
                           X_INOD2 = ndgln%x(( ELE2 - 1 ) * Mdims%cv_nloc + CV_ILOC2 )
                           IF( X_INOD2 == X_INOD ) THEN 
                              CV_OTHER_LOC( CV_ILOC ) = CV_ILOC2 
                           ENDIF
                        END DO
                     END DO
                ! Form approximate surface normal (NORMX,NORMY,NORMZ)
                     CALL DGSIMPLNORM( ELE, CV_SLOC2LOC, Mdims%cv_nloc, Mdims%cv_snloc, ndgln%x, &
                     X, Y, Z, NORMX(1), NORMX(2), NORMX(3) )
                ! Recalculate the normal...
                     XSL=0. ; YSL=0. ; ZSL=0.
                     DO CV_SILOC=1,Mdims%cv_snloc
                        CV_ILOC=CV_SLOC2LOC(CV_SILOC)
                        X_INOD=ndgln%x((ELE-1)*Mdims%x_nloc+CV_ILOC) 
                        XSL(CV_SILOC)=X(X_INOD)
                        IF(Mdims%ndim.GE.2) YSL(CV_SILOC)=Y(X_INOD)
                        IF(Mdims%ndim.GE.3) ZSL(CV_SILOC)=Z(X_INOD)
                     END DO
                     CALL DGSDETNXLOC2(Mdims%cv_snloc,CV_GIdims%sbcvngi, &
                     XSL,YSL,ZSL, &
                     CV_funs%sbcvfen, CV_funs%sbcvfenslx, CV_funs%sbcvfensly, CV_funs%sbcvfeweigh, SDETWE,SAREA, &
                     (Mdims%ndim==1), (Mdims%ndim==3), (Mdims%ndim==-2), &
                     SNORMXN(1,:),SNORMXN(2,:),SNORMXN(3,:), &
                     NORMX(1),NORMX(2),NORMX(3))
                 
                     DO SGI=1,CV_GIdims%sbcvngi
                        PSISGI_X=0.0
                        DO CV_SILOC=1,Mdims%cv_snloc
                           CV_ILOC=CV_SLOC2LOC( CV_SILOC )
                           DG_NOD=(ELE-1)*Mdims%cv_nloc+CV_ILOC
                           CV_ILOC2=CV_OTHER_LOC( CV_ILOC )
                           DG_INOD2=(ELE2-1)*Mdims%cv_nloc+CV_ILOC2
                           PSISGI_X(:)=PSISGI_X(:) + CV_funs%sbcvfen( CV_SILOC, SGI ) * (SOL_DERIV_X(:, DG_NOD)+SOL_DERIV_X(:, DG_INOD2))*0.5
                        END DO
                        N_DOT_SQ(SGI)= SUM(PSISGI_X(:)*SNORMXN(:,SGI))*(1.0*HDC)**2
                        IF (SELE2.GT.0) N_DOT_SQ(SGI)=0.0
                        S_INCOME(SGI)=1.0
                     END DO
                     DO CV_SILOC=1,Mdims%cv_snloc
                        CV_ILOC=CV_SLOC2LOC( CV_SILOC )
                        DO CV_SJLOC=1,Mdims%cv_snloc
                           CV_JLOC=CV_SLOC2LOC( CV_SJLOC )
                           CV_JLOC2=CV_OTHER_LOC( CV_JLOC )
                           DG_JNOD2=(ELE2-1)*Mdims%cv_nloc+CV_JLOC2
                           SVLN_IN=0.0
                           SVLN_OUT=0.0
                           DO SGI=1,CV_GIdims%sbcvngi
                              RR=CV_funs%sbcvfen(CV_SILOC,SGI)*N_DOT_SQ(SGI)*CV_funs%sbcvfen(CV_SJLOC,SGI)*SDETWE(SGI)
                              SVLN_IN =SVLN_IN  + RR*S_INCOME(SGI)
                              SVLN_OUT=SVLN_OUT + RR*(1.-S_INCOME(SGI))
                           END DO
               
                           !MAT_LOC(CV_ILOC,CV_JLOC)=MAT_LOC(CV_ILOC,CV_JLOC) - SVLN_IN
                           !RHS_CV_SHORT(CV_ILOC)=RHS_CV_SHORT(CV_ILOC) &
                           !   -SVLN_IN*DISTANCE_FUN(DG_JNOD2)
                           RHS_CV_SHORT(CV_ILOC)=RHS_CV_SHORT(CV_ILOC) + SVLN_IN
                           !X_INOD=ndgln%x((ELE-1)*Mdims%cv_nloc+CV_ILOC)
                   !DG_NOD=(ELE-1)*Mdims%cv_nloc+CV_ILOC
                   !IF (X(X_INOD).EQ.4.0.AND.Y(X_INOD).EQ.5.0) THEN
                           !WRITE(*,*) 'ELE,DISTANCE,X,Y:', ELE, DISTANCE_FUN(DG_NOD), X(X_INOD), Y(X_INOD), SVLN_IN
                           !END IF
                        END DO
                     END DO
                  END IF
               END DO Between_Elements_And_Boundary 
! Invert mass matrix...
! Solve MAT_LOC *CV_SOL = RHS_CV_SHORT
! STORE_MASS is overwritten by lu decomposition which used after the 1st solve. 
              GOTDEC = .FALSE.
              CALL SMLINNGOT( MAT_LOC, CV_SOL, RHS_CV_SHORT, Mdims%cv_nloc, IPIV, GOTDEC)
! Solve mass matrix systems...
              DO CV_ILOC=1,Mdims%cv_nloc
                 DG_NOD=(ELE-1)*Mdims%cv_nloc+CV_ILOC
                 DISTANCE_FUN(DG_NOD)=CV_SOL(CV_ILOC)
              END DO
            END DO  ! INNER_ITERATION
            END IF ! INTERFACE_ELE
            END DO  ! ELE loop 4
      !DISTANCE_FUN_OLD=DISTANCE_FUN
            ds%val = DISTANCE_FUN
            call halo_update(ds)
            DISTANCE_FUN = ds%val
         END DO ! OUT_ITER
      
         DISTANCE_FUN1=0.0
!   Distribute DISTANCE_FUN to the CV nodes...
         DO ELE=1,Mdims%totele ! ELE loop 5
            DO CV_ILOC=1,Mdims%cv_nloc
               CV_NOD=ndgln%cv((ELE-1)*Mdims%cv_nloc+CV_ILOC)
               DG_NOD=(ELE-1)*Mdims%cv_nloc+CV_ILOC
               !X_INOD=ndgln%x((ELE-1)*Mdims%cv_nloc+CV_ILOC)
               !DISTANCE_FUN(DG_NOD)=2.0-SQRT((X(X_INOD)-4.0)**2+(Y(X_INOD)-4.0)**2)
               DISTANCE_FUN1(CV_NOD)=DISTANCE_FUN1(CV_NOD)+DISTANCE_FUN(DG_NOD) * MASS_ELE(ELE) / MASS_CV(CV_NOD)
            END DO
         END DO ! ELE loop 5
         DO ELE=1,Mdims%totele ! ELE loop 6
            DO CV_ILOC=1,Mdims%cv_nloc
               CV_NOD=ndgln%cv((ELE-1)*Mdims%cv_nloc+CV_ILOC)
               DG_NOD=(ELE-1)*Mdims%cv_nloc+CV_ILOC
               DISTANCE_FUN(DG_NOD)=DISTANCE_FUN1(CV_NOD)
            END DO
            IF ( ((Mdims%ndim==2).AND.(Mdims%cv_nloc==6)).or.((Mdims%ndim==3).AND.(Mdims%cv_nloc==10)) ) THEN
               DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+2)=0.5*(DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+1)+DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+3))
               DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+4)=0.5*(DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+1)+DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+6))
               DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+5)=0.5*(DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+3)+DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+6))
               IF ( Mdims%cv_nloc == 10 ) THEN
                  DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+7)=0.5*(DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+1)+DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+10))
                  DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+8)=0.5*(DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+3)+DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+10))
                  DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+9)=0.5*(DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+6)+DISTANCE_FUN((ELE-1)*Mdims%cv_nloc+10))
               END IF
            END IF
         END DO ! ELE loop 6
         DISTANCE_FUN_OLD=DISTANCE_FUN
      END DO !ITIME
      
!---------------- end of calculation of the distance function -------
    
!ds => extract_scalar_field( state(1), 'ds'  )
      ds%val = distance_fun
!       DO ELE=1,Mdims%totele
!          state_nodes => ele_nodes(ds, ele)
!     DO CV_ILOC=1,Mdims%cv_nloc
!             ds%val( state_nodes( CV_ILOC )  ) = distance_fun(  (ELE-1)*Mdims%cv_nloc+CV_ILOC  )      
!          END DO
!       END DO 
      SOL_DERIV_X=0.
!Calculate the interface normal
      DO ELE=1,Mdims%totele ! ELE loop 7
         IF (INTERFACE_ELE(ELE)) THEN
! Calculate DevFuns%DETWEI,DevFuns%RA,NX,NY,NZ for element ELE
                call DETNLXR_PLUS_U(ELE, x_all % val, ndgln%x, CV_funs%cvweight, &
                       CV_funs%cvfen, CV_funs%cvfenlx_all, CV_funs%ufenlx_all, Devfuns)
 
            CALL LOC_1ST_DERIV_XYZ_DG_DERIV(DISTANCE_FUN, SOL_DERIV_X(1,:), SOL_DERIV_X(2,:), SOL_DERIV_X(3,:), &
                                           Mdims%ndim,  ndgln%cv, &
                                           Mdims%x_nloc, ndgln%x, &
                                           CV_GIdims%cv_ngi, Mdims%cv_nloc, &
                                           CV_funs%CVFEN, DevFuns%CVFENX_ALL(1,:,:), DevFuns%CVFENX_ALL(2,:,:), DevFuns%CVFENX_ALL(3,:,:), &
                                           X, Y, Z, &
                                           CV_GIdims%nface, FACE_ELE, CV_funs%cv_sloclist, Mdims%cv_snloc, &
                                           CV_GIdims%sbcvngi, CV_funs%sbcvfen, CV_funs%sbcvfenslx, CV_funs%sbcvfensly, CV_funs%sbcvfeweigh, &
                                           ELE, DevFuns%DETWEI, 2, angle)
            
         END IF
      END DO ! ELE loop 7
             
      dx%val = SOL_DERIV_X(1,:)
      dy%val = SOL_DERIV_X(2,:)
      IF(Mdims%ndim.GE.3) dz%val = SOL_DERIV_X(3,:)
      call halo_update(dx)
      call halo_update(dy)
      IF(Mdims%ndim.GE.3) call halo_update(dz)
      SOL_DERIV_X(1,:) = dx%val 
      SOL_DERIV_X(2,:) = dy%val 
      IF(Mdims%ndim.GE.3) SOL_DERIV_X(3,:) = dz%val 
!   Distribute DG DERIV to the CV nodes...
      SOL_DERIV_X1=0.0
      DO ELE=1,Mdims%totele ! ELE loop 8
         DO CV_ILOC=1,Mdims%cv_nloc
            CV_NOD=ndgln%cv((ELE-1)*Mdims%cv_nloc+CV_ILOC)
            DG_NOD=(ELE-1)*Mdims%cv_nloc+CV_ILOC
            SOL_DERIV_X1( :, CV_NOD)=SOL_DERIV_X1( :, CV_NOD)+SOL_DERIV_X( :, DG_NOD) * MASS_ELE(ELE) / MASS_CV(CV_NOD)
         END DO
      END DO ! ELE loop 8
      DO ELE=1,Mdims%totele ! ELE loop 9
         DO CV_ILOC=1,Mdims%cv_nloc
            CV_NOD=ndgln%cv((ELE-1)*Mdims%cv_nloc+CV_ILOC)
            DG_NOD=(ELE-1)*Mdims%cv_nloc+CV_ILOC
            SOL_DERIV_X( :, DG_NOD)=SOL_DERIV_X1( :, CV_NOD)
            NORMALIZATION(DG_NOD)=SQRT(SUM(SOL_DERIV_X(:,DG_NOD)**2))
            IF (NORMALIZATION(DG_NOD).LE.1.E-4) THEN
               SOL_DERIV_X(:, DG_NOD)=0.0
            END IF
            SOL_DERIV_X( :, DG_NOD)=SOL_DERIV_X( :, DG_NOD)/MAX(TOLER,NORMALIZATION(DG_NOD))
            !X_INOD=ndgln%x((ELE-1)*Mdims%cv_nloc+CV_ILOC)
            !SOL_DERIV_X(1, DG_NOD)=-(X(X_INOD)-4.0)/MAX(SQRT((X(X_INOD)-4.0)**2+(Y(X_INOD)-4.0)**2),1.E-10)
            !SOL_DERIV_X(2, DG_NOD)=-(Y(X_INOD)-4.0)/MAX(SQRT((X(X_INOD)-4.0)**2+(Y(X_INOD)-4.0)**2),1.E-10)
         END DO
         
         IF ( ((Mdims%ndim==2).AND.(Mdims%cv_nloc==6)).or.((Mdims%ndim==3).AND.(Mdims%cv_nloc==10)) ) THEN
            SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+2)=0.5*(SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+1)+SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+3))
            SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+4)=0.5*(SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+1)+SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+6))
            SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+5)=0.5*(SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+3)+SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+6))
            IF ( Mdims%cv_nloc == 10 ) THEN
               SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+7)=0.5*(SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+1)+SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+10))
               SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+8)=0.5*(SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+3)+SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+10))
               SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+9)=0.5*(SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+6)+SOL_DERIV_X(:, (ELE-1)*Mdims%cv_nloc+10))
            END IF
         END IF
      END DO ! ELE loop 9
      DO ELE=1,Mdims%totele
         state_nodes => ele_nodes(dx, ele)
         DO CV_ILOC=1,Mdims%cv_nloc
            dx%val( state_nodes( CV_ILOC )  ) = sol_deriv_x( 1, (ELE-1)*Mdims%cv_nloc+CV_ILOC  )
            dy%val( state_nodes( CV_ILOC )  ) = sol_deriv_x( 2, (ELE-1)*Mdims%cv_nloc+CV_ILOC  )
            IF(Mdims%ndim.GE.3) dz%val( state_nodes( CV_ILOC )  ) = sol_deriv_x( 3, (ELE-1)*Mdims%cv_nloc+CV_ILOC  )
         END DO
      END DO 
      !dx%val = SOL_DERIV_X
      !dy%val = SOL_DERIV_Y
      !dz%val = SOL_DERIV_Z
      CURV=0.
      DO ELE=1,Mdims%totele ! ELE loop 10
         IF(INTERFACE_ELE(ELE)) THEN
                ! Calculate DevFuns%DETWEI,DevFuns%RA,NX,NY,NZ for element ELE
                call DETNLXR_PLUS_U(ELE, x_all % val, ndgln%x, CV_funs%cvweight, &
                       CV_funs%cvfen, CV_funs%cvfenlx_all, CV_funs%ufenlx_all, Devfuns)
! calculate curvature:
            CALL LOC_1ST_DERIV_XYZ_DG_CURV(CURV, DISTANCE_FUN, SOL_DERIV_X(1,:), SOL_DERIV_X(2,:), SOL_DERIV_X(3,:), &
                         Mdims%ndim,  ndgln%cv, &
                         Mdims%x_nloc, ndgln%x,&
                         CV_GIdims%cv_ngi, Mdims%cv_nloc, &
                         CV_funs%CVFEN, DevFuns%CVFENX_ALL(1,:,:), DevFuns%CVFENX_ALL(2,:,:), DevFuns%CVFENX_ALL(3,:,:), &
                         X, Y, Z, &
                         CV_GIdims%nface, FACE_ELE, CV_funs%cv_sloclist, Mdims%cv_snloc, &
                         CV_GIdims%sbcvngi, CV_funs%sbcvfen, CV_funs%sbcvfenslx, CV_funs%sbcvfensly, CV_funs%sbcvfeweigh, &
                         ELE, DevFuns%DETWEI)
         END IF
      END DO ! ELE loop 10
      dk%val = curv
      call halo_update(dk)
      curv = dk%val 
!   Distribute curvature to the CV nodes...
      CURVATURE=0.0
      DO ELE=1,Mdims%totele ! ELE loop 11
         DO CV_ILOC=1,Mdims%cv_nloc
            CV_NOD=ndgln%cv((ELE-1)*Mdims%cv_nloc+CV_ILOC)
            DG_NOD=(ELE-1)*Mdims%cv_nloc+CV_ILOC
            CURVATURE(CV_NOD)=CURVATURE(CV_NOD)+CURV(DG_NOD) * MASS_ELE(ELE) / MASS_CV(CV_NOD)
         END DO
      END DO ! ELE loop 11
      
      DO ELE=1,Mdims%totele ! ELE loop 12
         IF( ((Mdims%ndim==2).AND.(Mdims%cv_nloc==6)).or.((Mdims%ndim==3).AND.(Mdims%cv_nloc==10)) ) THEN
            CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+2))=0.5*(CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+1)) &
                                                       +CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+3)))
            CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+4))=0.5*(CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+1)) &
                                                       +CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+6)))
            CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+5))=0.5*(CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+3)) &
                                                       +CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+6)))
            IF ( Mdims%cv_nloc == 10 ) THEN
               CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+7))=0.5*(CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+1)) &
                                                          +CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+10)))
               CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+8))=0.5*(CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+3)) &
                                                          +CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+10)))
               CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+9))=0.5*(CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+6)) &
                                                          +CURVATURE(ndgln%cv((ELE-1)*Mdims%cv_nloc+10)))
            END IF
         END IF
      END DO ! ELE loop 12
!dk => extract_scalar_field( state(1), 'dk'  )
      ck%val = curvature
      !IF_USE_PRESSURE_FORCE: IF ( USE_PRESSURE_FORCE ) THEN
         PLIKE_GRAD_SOU_COEF = PLIKE_GRAD_SOU_COEF + SUF_TENSION_COEF * CURVATURE
         PLIKE_GRAD_SOU_GRAD = PLIKE_GRAD_SOU_GRAD + VOLUME_FRAC
         
      !END IF IF_USE_PRESSURE_FORCE
      DEALLOCATE( X, Y, Z )
      DEALLOCATE( DGI_X, UD, PSISGI_X, NORMX, SNORMXN )
      DEALLOCATE( CV_OTHER_LOC, CV_SLOC2LOC )
      DEALLOCATE( MASS_CV, FACE_ELE, MASS_ELE)
      DEALLOCATE( CURVATURE, CURV, INTERFACE_ELE, INTERFACE_ELE2, DISTANCE_FUN, DISTANCE_FUN_OLD, DISTANCE_FUN1 )
      DEALLOCATE( SOL_DERIV_X, SOL_DERIV_X1, SIGN_FUN_GI, SIGN_FUN_SGI )
      DEALLOCATE( XSL, YSL, ZSL )
      DEALLOCATE( S_INCOME, SDETWE, N_DOT_SQ, ST, NORMALIZATION, DIFF_COEF)
      DEALLOCATE( MAT_LOC, RHS_CV_SHORT, CV_SOL )
      call deallocate_multi_dev_shape_funs(Devfuns)
     ewrite(3,*) 'Leaving SURFACE_TENSION_WAPPER_NEW'
contains

      SUBROUTINE LOC_1ST_DERIV_XYZ_DG_CURV(CURV, DISTANCE_FUN, SOL_DERIV_X, SOL_DERIV_Y, SOL_DERIV_Z, &
                                           NDIM, CV_NDGLN, &
                                           X_NLOC, X_NDGLN, &
                                           CV_NGI, CV_NLOC, &
                                           CVFEN, CVFENX, CVFENY, CVFENZ, &
                                           X, Y, Z, &
                                           NFACE, FACE_ELE, CV_SLOCLIST, CV_SNLOC, &
                                           SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
                                           ELE, DETWEI) 

      IMPLICIT NONE
      INTEGER, intent( in ) :: NDIM,  X_NLOC, CV_NGI, CV_NLOC, &
           CV_SNLOC, SBCVNGI, NFACE, ELE
      REAL, DIMENSION( :), intent( inout ) :: CURV
      REAL, DIMENSION( : ), intent( in ) :: DISTANCE_FUN
      REAL, DIMENSION( :), intent( in ) :: SOL_DERIV_X, SOL_DERIV_Y, SOL_DERIV_Z
      INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( :,: ), intent( in ) ::  CV_SLOCLIST
      INTEGER, DIMENSION( :,: ), intent( in ) ::  FACE_ELE
      REAL, DIMENSION( :, : ), intent( in ) :: CVFEN, CVFENX, CVFENY, CVFENZ
      REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( :, : ), intent( in ) :: SBCVFEN, SBCVFENSLX, SBCVFENSLY
      REAL, DIMENSION( : ), intent( in ) :: SBCVFEWEIGH
      REAL, DIMENSION( : ), intent( in ) :: DETWEI
      ! Local variables
 
  
      REAL, DIMENSION( : ), allocatable :: SOL, RHS_DG, DN_GI, &
                                           SOL_X, SOL_Y, SOL_Z, &
                                       RHS_DG_X, RHS_DG_Y, RHS_DG_Z
      INTEGER, DIMENSION( : ), allocatable :: CV_SLOC2LOC, CV_OTHER_LOC
      REAL, DIMENSION( :, : ), allocatable :: MASS
      REAL  :: VLNN, VLNX, VLNY, VLNZ, VLNX2, VLNY2, VLNZ2, &
               RR, SXNN, SYNN, SZNN, RDISTX, RDISTY, RDISTZ
      INTEGER :: CV_ILOC, CV_JLOC, DG_JNOD, CV_INOD, X_INOD, &
                     CV_ILOC2, CV_JLOC2,DG_JNOD2, &  
                     CV_SILOC, CV_SJLOC, &
                     ELE2, SELE, SELE2, &
                     GI, SGI, IFACE, X_INOD2, DG_NOD
      LOGICAL :: GOTDEC
      REAL, DIMENSION( : ), allocatable :: XSL, YSL, ZSL, SNORMXN, SNORMYN, SNORMZN, SDETWE
      REAL ::  NORMX, NORMY, NORMZ, SAREA
      INTEGER :: IPIV(CV_NLOC)
      REAL :: DGI_X, DGI_Y, DGI_Z, SQRT_RR, DN!, SDGI_X, SDGI_Y, SDGI_Z, SDN
  
      ALLOCATE(SOL(CV_NLOC))
      ALLOCATE(RHS_DG(CV_NLOC))
      ALLOCATE(SOL_X(CV_NLOC),SOL_Y(CV_NLOC),SOL_Z(CV_NLOC))
      ALLOCATE(RHS_DG_X(CV_NLOC),RHS_DG_Y(CV_NLOC),RHS_DG_Z(CV_NLOC))
      ALLOCATE(MASS(CV_NLOC,CV_NLOC))
      ALLOCATE( CV_SLOC2LOC( CV_SNLOC ))
      ALLOCATE( CV_OTHER_LOC( CV_NLOC ))
      ALLOCATE( XSL( CV_SNLOC ))
      ALLOCATE( YSL( CV_SNLOC ))
      ALLOCATE( ZSL( CV_SNLOC ))
      ALLOCATE( SNORMXN( SBCVNGI ))
      ALLOCATE( SNORMYN( SBCVNGI ))
      ALLOCATE( SNORMZN( SBCVNGI ))
      ALLOCATE( SDETWE( SBCVNGI ))
      ALLOCATE( DN_GI( CV_NGI ))

  
      RHS_DG=0.0
      RHS_DG_X=0.0
      RHS_DG_Y=0.0
      RHS_DG_Z=0.0
      MASS=0.0

      DO GI=1,CV_NGI
         DGI_X=0.0
         DGI_Y=0.0
         DGI_Z=0.0
         DO CV_ILOC=1,CV_NLOC
            DG_NOD=(ELE-1)*CV_NLOC+CV_ILOC
            DGI_X=DGI_X + CVFENX( CV_ILOC, GI ) * DISTANCE_FUN(DG_NOD)
            IF(NDIM.GE.2) DGI_Y=DGI_Y + CVFENY( CV_ILOC, GI ) * DISTANCE_FUN(DG_NOD)
            IF(NDIM.GE.3) DGI_Z=DGI_Z + CVFENZ( CV_ILOC, GI ) * DISTANCE_FUN(DG_NOD)
         END DO        
         RR = DGI_X**2
         IF(NDIM.GE.2) RR = RR+ DGI_Y**2
         IF(NDIM.GE.3) RR = RR+ DGI_Z**2
         SQRT_RR=SQRT(RR)
         IF (SQRT_RR.EQ.0.0) THEN
            DN=1.0
         ELSE
            DN=1.0/SQRT_RR
         END IF
         DN_GI(GI)=DN
      END DO

      DO CV_ILOC=1,CV_NLOC
         DO CV_JLOC=1,CV_NLOC
            DG_JNOD=(ELE-1)*CV_NLOC+CV_JLOC
            
            VLNN=0.0
            VLNX=0.0
            VLNY=0.0
            VLNZ=0.0
            VLNX2=0.0
            VLNY2=0.0
            VLNZ2=0.0
            DO GI=1,CV_NGI
                   VLNN= VLNN+CVFEN( CV_ILOC, GI )*CVFEN( CV_JLOC, GI )* DETWEI(GI)

     !              VLNX= VLNX-DN_GI(GI)*CVFENX( CV_ILOC, GI )*CVFENX( CV_JLOC, GI )* DETWEI(GI)
     !IF(NDIM.GE.2) VLNY= VLNY-DN_GI(GI)*CVFENY( CV_ILOC, GI )*CVFENY( CV_JLOC, GI )* DETWEI(GI)
     !IF(NDIM.GE.3) VLNZ= VLNZ-DN_GI(GI)*CVFENZ( CV_ILOC, GI )*CVFENZ( CV_JLOC, GI )* DETWEI(GI)

     !              VLNX= VLNX-CVFENX( CV_ILOC, GI )*CVFEN( CV_JLOC, GI )* DETWEI(GI)
     !IF(NDIM.GE.2) VLNY= VLNY-CVFENY( CV_ILOC, GI )*CVFEN( CV_JLOC, GI )* DETWEI(GI)
     !IF(NDIM.GE.3) VLNZ= VLNZ-CVFENZ( CV_ILOC, GI )*CVFEN( CV_JLOC, GI )* DETWEI(GI)

                   VLNX= VLNX+CVFEN( CV_ILOC, GI )*CVFENX( CV_JLOC, GI )* DETWEI(GI)
     IF(NDIM.GE.2) VLNY= VLNY+CVFEN( CV_ILOC, GI )*CVFENY( CV_JLOC, GI )* DETWEI(GI)
     IF(NDIM.GE.3) VLNZ= VLNZ+CVFEN( CV_ILOC, GI )*CVFENZ( CV_JLOC, GI )* DETWEI(GI)

                   VLNX2= VLNX2+CVFEN( CV_ILOC, GI )*CVFENX( CV_JLOC, GI )* DETWEI(GI)
     IF(NDIM.GE.2) VLNY2= VLNY2+CVFEN( CV_ILOC, GI )*CVFENY( CV_JLOC, GI )* DETWEI(GI)
     IF(NDIM.GE.3) VLNZ2= VLNZ2+CVFEN( CV_ILOC, GI )*CVFENZ( CV_JLOC, GI )* DETWEI(GI)
            END DO
            MASS(CV_ILOC,CV_JLOC)=MASS(CV_ILOC,CV_JLOC) + VLNN
            RHS_DG(CV_ILOC)=RHS_DG(CV_ILOC)+ VLNX*SOL_DERIV_X(DG_JNOD)+ VLNY*SOL_DERIV_Y(DG_JNOD) + VLNZ*SOL_DERIV_Z(DG_JNOD)
            !RHS_DG(CV_ILOC)=RHS_DG(CV_ILOC)+ (VLNX + VLNY + VLNZ)*DISTANCE_FUN(DG_JNOD)

            RHS_DG_X(CV_ILOC)=RHS_DG_X(CV_ILOC)+ VLNX2*DISTANCE_FUN(DG_JNOD)
            RHS_DG_Y(CV_ILOC)=RHS_DG_Y(CV_ILOC)+ VLNY2*DISTANCE_FUN(DG_JNOD)
            RHS_DG_Z(CV_ILOC)=RHS_DG_Z(CV_ILOC)+ VLNZ2*DISTANCE_FUN(DG_JNOD)
         END DO
      END DO

!         SOL_X=0.
!         SOL_Y=0.
!         SOL_Z=0.
!         GOTDEC = .FALSE.
!         CALL SMLINNGOT( MASS, SOL_X, RHS_DG_X, CV_NLOC, CV_NLOC, IPIV, GOTDEC)
!         GOTDEC=.TRUE.
!         IF(NDIM.GE.2) CALL SMLINNGOT( MASS, SOL_Y, RHS_DG_Y, CV_NLOC, CV_NLOC, IPIV, GOTDEC)
!         IF(NDIM.GE.3) CALL SMLINNGOT( MASS, SOL_Z, RHS_DG_Z, CV_NLOC, CV_NLOC, IPIV, GOTDEC)

  ! Calculate DETWEI,RA,NX,NY,NZ for boundary
      Between_Elements_And_Boundary: DO IFACE = 1, NFACE
         ELE2  = FACE_ELE( IFACE, ELE )
         SELE2 = MAX( 0, - ELE2 )
         SELE  = SELE2
         ELE2  = MAX( 0, + ELE2 )

            ! The surface nodes on element face IFACE. 
         CV_SLOC2LOC( : ) = CV_SLOCLIST( IFACE, : )

         CV_OTHER_LOC = 0
         IF(SELE2 == 0) THEN
            DO CV_SILOC = 1, CV_SNLOC
               CV_ILOC = CV_SLOC2LOC( CV_SILOC )
               X_INOD = X_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )
               DO CV_ILOC2 = 1, CV_NLOC
                  X_INOD2 = X_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_ILOC2 )
                  IF( X_INOD2 == X_INOD ) THEN 
                     CV_OTHER_LOC( CV_ILOC ) = CV_ILOC2 
                  ENDIF
               END DO
            END DO
         END IF

            ! Form approximate surface normal (NORMX,NORMY,NORMZ)
         CALL DGSIMPLNORM( ELE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, X_NDGLN, &
                 X, Y, Z, NORMX, NORMY, NORMZ )

            ! Recalculate the normal...
         DO CV_SILOC=1,CV_SNLOC
            CV_ILOC=CV_SLOC2LOC(CV_SILOC)
            X_INOD=X_NDGLN((ELE-1)*X_NLOC+CV_ILOC) 
            XSL(CV_SILOC)=X(X_INOD)
            YSL(CV_SILOC)=Y(X_INOD)
            ZSL(CV_SILOC)=Z(X_INOD)
         END DO

         CALL DGSDETNXLOC2(CV_SNLOC,SBCVNGI, &
                 XSL,YSL,ZSL, &
                 SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SDETWE,SAREA, &
                 (NDIM==1), (NDIM==3), (NDIM==-2), &
                 SNORMXN,SNORMYN,SNORMZN, &
                 NORMX,NORMY,NORMZ)

!         DO SGI=1,SBCVNGI
!            SDGI_X=0.0
!            SDGI_Y=0.0
!            SDGI_Z=0.0
!            DO CV_SILOC=1,CV_SNLOC
!               CV_ILOC=CV_SLOC2LOC( CV_SILOC )
!               DG_NOD=(ELE-1)*CV_NLOC+CV_ILOC
!               SDGI_X=SDGI_X + SBCVFENSLX( CV_SILOC, SGI ) * DISTANCE_FUN(DG_NOD)
!            END DO

!         END DO

         DO CV_SILOC=1,CV_SNLOC
            CV_ILOC=CV_SLOC2LOC( CV_SILOC )
            CV_INOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
            DO CV_SJLOC=1,CV_SNLOC
               CV_JLOC=CV_SLOC2LOC( CV_SJLOC )
               DG_JNOD=(ELE-1)*CV_NLOC+CV_JLOC

               SXNN=0.0
               SYNN=0.0
               SZNN=0.0
               DO SGI=1,SBCVNGI
                  RR=SBCVFEN(CV_SILOC,SGI)*SBCVFEN(CV_SJLOC,SGI)*SDETWE(SGI)
                  SXNN=SXNN+SNORMXN(SGI)*RR
                  SYNN=SYNN+SNORMYN(SGI)*RR
                  SZNN=SZNN+SNORMZN(SGI)*RR
               END DO
               IF(SELE.GT.0) THEN
                  RDISTX=0.0!SOL_DERIV_X(DG_JNOD) 
                  RDISTY=0.0!SOL_DERIV_Y(DG_JNOD)
                  RDISTZ=0.0!SOL_DERIV_Z(DG_JNOD)
                  !RDISTX=SOL_DERIV_X(DG_JNOD) 
                  !RDISTY=SOL_DERIV_Y(DG_JNOD)
                  !RDISTZ=SOL_DERIV_Z(DG_JNOD)

!                  DGI_X=0.0
!                  DGI_Y=0.0
!                  DGI_Z=0.0
!                  DO GI=1,CV_NGI
!                     DGI_X=DGI_X + CVFEN( CV_JLOC, GI ) * SOL_X(CV_JLOC)
!                     IF(NDIM.GE.2) DGI_Y=DGI_Y + CVFEN( CV_JLOC, GI ) * SOL_Y(CV_JLOC)
!                     IF(NDIM.GE.3) DGI_Z=DGI_Z + CVFEN( CV_JLOC, GI ) * SOL_Z(CV_JLOC)  
!                  END DO
!                  RR = DGI_X**2
!                  IF(NDIM.GE.2) RR = RR+ DGI_Y**2
!                  IF(NDIM.GE.3) RR = RR+ DGI_Z**2
!                  SQRT_RR=SQRT(RR)
!                  IF (SQRT_RR.EQ.0.0) THEN
!                     DN=1.0
!                  ELSE
!                     DN=1.0/SQRT_RR
!                  END IF
                  !RDISTX=DN*DGI_X
                  !RDISTY=DN*DGI_Y
                  !RDISTY=DN*DGI_Z

               ELSE
                  CV_JLOC2=CV_OTHER_LOC( CV_JLOC )
                  DG_JNOD2=(ELE2-1)*CV_NLOC+CV_JLOC2
                  RDISTX=0.5*(-SOL_DERIV_X(DG_JNOD)+SOL_DERIV_X(DG_JNOD2))
                  RDISTY=0.5*(-SOL_DERIV_Y(DG_JNOD)+SOL_DERIV_Y(DG_JNOD2))
                  RDISTZ=0.5*(-SOL_DERIV_Z(DG_JNOD)+SOL_DERIV_Z(DG_JNOD2))
                  !RDISTX=0.5*(SOL_DERIV_X(DG_JNOD)+SOL_DERIV_X(DG_JNOD2))
                  !RDISTY=0.5*(SOL_DERIV_Y(DG_JNOD)+SOL_DERIV_Y(DG_JNOD2))
                  !RDISTZ=0.5*(SOL_DERIV_Z(DG_JNOD)+SOL_DERIV_Z(DG_JNOD2))

!                  DGI_X=0.0
!                  DGI_Y=0.0
!                  DGI_Z=0.0
!                  SDGI_X=0.0
!                  SDGI_Y=0.0
!                  SDGI_Z=0.0
!                  DO GI=1,CV_NGI
!                     DGI_X=DGI_X + CVFENX( CV_JLOC, GI ) * DISTANCE_FUN(DG_JNOD)
!                     IF(NDIM.GE.2) DGI_Y=DGI_Y + CVFENY( CV_JLOC, GI ) * DISTANCE_FUN(DG_JNOD)
!                     IF(NDIM.GE.3) DGI_Z=DGI_Z + CVFENZ( CV_JLOC, GI ) * DISTANCE_FUN(DG_JNOD)  
!                     SDGI_X=SDGI_X + CVFENX( CV_JLOC2, GI ) * DISTANCE_FUN(DG_JNOD2)
!                     IF(NDIM.GE.2) SDGI_Y=SDGI_Y + CVFENY( CV_JLOC2, GI ) * DISTANCE_FUN(DG_JNOD2)
!                     IF(NDIM.GE.3) SDGI_Z=SDGI_Z + CVFENZ( CV_JLOC2, GI ) * DISTANCE_FUN(DG_JNOD2)   
!                  END DO
!                  RR = DGI_X**2
!                  IF(NDIM.GE.2) RR = RR+ DGI_Y**2
!                  IF(NDIM.GE.3) RR = RR+ DGI_Z**2
!                  SQRT_RR=SQRT(RR)
!                  IF (SQRT_RR.EQ.0.0) THEN
!                     DN=1.0
!                  ELSE
!                     DN=1.0/SQRT_RR
!                  END IF
!                  RR = SDGI_X**2
!                  IF(NDIM.GE.2) RR = RR+ SDGI_Y**2
!                  IF(NDIM.GE.3) RR = RR+ SDGI_Z**2
!                  SQRT_RR=SQRT(RR)
!                  IF (SQRT_RR.EQ.0.0) THEN
!                     SDN=1.0
!                  ELSE
!                     SDN=1.0/SQRT_RR
!                  END IF
!                  RDISTX=0.5*(DN*DGI_X+SDN*SDGI_X)
!                  RDISTY=0.5*(DN*DGI_Y+SDN*SDGI_Y)
!                  RDISTY=0.5*(DN*DGI_Z+SDN*SDGI_Z)
               ENDIF
               RHS_DG(CV_ILOC)=RHS_DG(CV_ILOC) +  SXNN*RDISTX + SYNN*RDISTY + SZNN*RDISTZ !ZX 28062013
            END DO
         END DO

      END DO Between_Elements_And_Boundary 

      SOL=0.
      GOTDEC = .FALSE.
      CALL SMLINNGOT( MASS, SOL, RHS_DG, CV_NLOC, IPIV, GOTDEC)

      Between_Elements_And_Boundary1: DO IFACE = 1, -1!NFACE
         ELE2  = FACE_ELE( IFACE, ELE )
         SELE2 = MAX( 0, - ELE2 )
         SELE  = SELE2
         ELE2  = MAX( 0, + ELE2 )

         CV_SLOC2LOC( : ) = CV_SLOCLIST( IFACE, : )

         DO CV_SILOC=1,CV_SNLOC
            CV_ILOC=CV_SLOC2LOC( CV_SILOC )
            IF(SELE.GT.0) THEN
               SOL(CV_ILOC)=0.0
            ENDIF
         END DO

      END DO Between_Elements_And_Boundary1 

      DO CV_ILOC=1,CV_NLOC
         CURV((ELE-1)*CV_NLOC+CV_ILOC)=SOL(CV_ILOC)
      END DO
 

      DEALLOCATE(SOL)
      DEALLOCATE(RHS_DG)
      DEALLOCATE(SOL_X,SOL_Y,SOL_Z)
      DEALLOCATE(RHS_DG_X,RHS_DG_Y,RHS_DG_Z)
      DEALLOCATE(MASS)
      DEALLOCATE( CV_SLOC2LOC)
      DEALLOCATE( CV_OTHER_LOC)
      DEALLOCATE( XSL)
      DEALLOCATE( YSL)
      DEALLOCATE( ZSL)
      DEALLOCATE( SNORMXN)
      DEALLOCATE( SNORMYN)
      DEALLOCATE( SNORMZN)
      DEALLOCATE( SDETWE)
      DEALLOCATE( DN_GI)

        
      END SUBROUTINE LOC_1ST_DERIV_XYZ_DG_CURV

      SUBROUTINE LOC_1ST_DERIV_XYZ_DG_DERIV(DISTANCE_FUN, SOL_DERIV_X, SOL_DERIV_Y, SOL_DERIV_Z, &
                                           NDIM, CV_NDGLN, &
                                           X_NLOC, X_NDGLN, &
                                           CV_NGI, CV_NLOC, &
                                           CVFEN, CVFENX, CVFENY, CVFENZ, &
                                           X, Y, Z, &
                                           NFACE, FACE_ELE, CV_SLOCLIST, CV_SNLOC, &
                                           SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
                                           ELE, DETWEI, factor, angle) 

      IMPLICIT NONE
      INTEGER, intent( in ) :: NDIM,  X_NLOC, CV_NGI, CV_NLOC, &
           CV_SNLOC, SBCVNGI, NFACE, ELE, factor
      REAL, intent( in ) :: angle
      REAL, DIMENSION( : ), intent( in ) :: DISTANCE_FUN
      REAL, DIMENSION( : ), intent( inout ) :: SOL_DERIV_X, SOL_DERIV_Y, SOL_DERIV_Z
      INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( :,: ), intent( in ) ::  CV_SLOCLIST
      INTEGER, DIMENSION( :,: ), intent( in ) ::  FACE_ELE
      REAL, DIMENSION( :, :), intent( in ) :: CVFEN, CVFENX, CVFENY, CVFENZ
      REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( :, : ), intent( in ) :: SBCVFEN, SBCVFENSLX, SBCVFENSLY
      REAL, DIMENSION( : ), intent( in ) :: SBCVFEWEIGH
      REAL, DIMENSION( : ), intent( in ) :: DETWEI
      ! Local variables
 
  
      REAL, DIMENSION( : ), allocatable :: SOL_X, SOL_Y, SOL_Z, &
                                       RHS_DG_X, RHS_DG_Y, RHS_DG_Z
      INTEGER, DIMENSION( : ), allocatable :: CV_SLOC2LOC, CV_OTHER_LOC
      REAL, DIMENSION( :, : ), allocatable :: MASS, INV_MASS
      REAL  :: VLNN, VLNX, VLNY, VLNZ, &
           RR, RDIST, SXNN, SYNN, SZNN, SXTT, SYTT, SZTT
      INTEGER :: CV_ILOC, CV_JLOC, DG_JNOD, CV_INOD, X_INOD, &
                     CV_ILOC2, CV_JLOC2,DG_JNOD2, &  
                     CV_SILOC, CV_SJLOC, &
                     ELE2, SELE, SELE2, &
                     GI, SGI, IFACE, X_INOD2, DG_INOD
      LOGICAL :: GOTDEC
      REAL, DIMENSION( : ), allocatable :: XSL, YSL, ZSL, SNORMXN, SNORMYN, SNORMZN, SDETWE
      REAL ::  NORMX, NORMY, NORMZ, SAREA
      INTEGER :: IPIV(CV_NLOC)
      REAL, PARAMETER :: PI=3.14159265359
  
      ALLOCATE(SOL_X(CV_NLOC),SOL_Y(CV_NLOC),SOL_Z(CV_NLOC))
      ALLOCATE(RHS_DG_X(CV_NLOC),RHS_DG_Y(CV_NLOC),RHS_DG_Z(CV_NLOC))
      ALLOCATE(MASS(CV_NLOC,CV_NLOC),INV_MASS(CV_NLOC,CV_NLOC))
      ALLOCATE( CV_SLOC2LOC( CV_SNLOC ))
      ALLOCATE( CV_OTHER_LOC( CV_NLOC ))
      ALLOCATE( XSL( CV_SNLOC ))
      ALLOCATE( YSL( CV_SNLOC ))
      ALLOCATE( ZSL( CV_SNLOC ))
      ALLOCATE( SNORMXN( SBCVNGI ))
      ALLOCATE( SNORMYN( SBCVNGI ))
      ALLOCATE( SNORMZN( SBCVNGI ))
      ALLOCATE( SDETWE( SBCVNGI ))
  
      RHS_DG_X=0.0
      RHS_DG_Y=0.0
      RHS_DG_Z=0.0
      MASS=0.0
          

      DO CV_ILOC=1,CV_NLOC
         DO CV_JLOC=1,CV_NLOC
            DG_JNOD=(ELE-1)*CV_NLOC+CV_JLOC
            VLNN=0.0
            VLNX=0.0
            VLNY=0.0
            VLNZ=0.0
            DO GI=1,CV_NGI
                   VLNN= VLNN+CVFEN( CV_ILOC, GI )*CVFEN( CV_JLOC, GI )* DETWEI(GI)
     !              VLNX= VLNX-CVFENX( CV_ILOC, GI )*CVFEN( CV_JLOC, GI )* DETWEI(GI)
     !IF(NDIM.GE.2) VLNY= VLNY-CVFENY( CV_ILOC, GI )*CVFEN( CV_JLOC, GI )* DETWEI(GI)
     !IF(NDIM.GE.3) VLNZ= VLNZ-CVFENZ( CV_ILOC, GI )*CVFEN( CV_JLOC, GI )* DETWEI(GI)

                   VLNX= VLNX+CVFEN( CV_ILOC, GI )*CVFENX( CV_JLOC, GI )* DETWEI(GI)
     IF(NDIM.GE.2) VLNY= VLNY+CVFEN( CV_ILOC, GI )*CVFENY( CV_JLOC, GI )* DETWEI(GI)
     IF(NDIM.GE.3) VLNZ= VLNZ+CVFEN( CV_ILOC, GI )*CVFENZ( CV_JLOC, GI )* DETWEI(GI)
            END DO
            MASS(CV_ILOC,CV_JLOC)=MASS(CV_ILOC,CV_JLOC) + VLNN
            RHS_DG_X(CV_ILOC)=RHS_DG_X(CV_ILOC)+ VLNX*DISTANCE_FUN(DG_JNOD)
            RHS_DG_Y(CV_ILOC)=RHS_DG_Y(CV_ILOC)+ VLNY*DISTANCE_FUN(DG_JNOD)
            RHS_DG_Z(CV_ILOC)=RHS_DG_Z(CV_ILOC)+ VLNZ*DISTANCE_FUN(DG_JNOD)
         END DO
      END DO

  ! Calculate DETWEI,RA,NX,NY,NZ for boundary
      Between_Elements_And_Boundary: DO IFACE = 1, NFACE
         ELE2  = FACE_ELE( IFACE, ELE )
         SELE2 = MAX( 0, - ELE2 )
         SELE  = SELE2
         ELE2  = MAX( 0, + ELE2 )
!         print*, nface, iface, ele, FACE_ELE( IFACE, ELE ), ele2, sele2
            ! The surface nodes on element face IFACE. 
         CV_SLOC2LOC( : ) = CV_SLOCLIST( IFACE, : )

         CV_OTHER_LOC = 0
         IF(SELE2 == 0) THEN
            DO CV_SILOC = 1, CV_SNLOC
               CV_ILOC = CV_SLOC2LOC( CV_SILOC )
               X_INOD = X_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )
               DO CV_ILOC2 = 1, CV_NLOC
                  X_INOD2 = X_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_ILOC2 )
                  IF( X_INOD2 == X_INOD ) THEN 
                     CV_OTHER_LOC( CV_ILOC ) = CV_ILOC2 
                  ENDIF
               END DO
            END DO
         END IF

            ! Form approximate surface normal (NORMX,NORMY,NORMZ)
         CALL DGSIMPLNORM( ELE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, X_NDGLN, &
                 X, Y, Z, NORMX, NORMY, NORMZ )

            ! Recalculate the normal...
         DO CV_SILOC=1,CV_SNLOC
            CV_ILOC=CV_SLOC2LOC(CV_SILOC)
            X_INOD=X_NDGLN((ELE-1)*X_NLOC+CV_ILOC) 
            XSL(CV_SILOC)=X(X_INOD)
            YSL(CV_SILOC)=Y(X_INOD)
            ZSL(CV_SILOC)=Z(X_INOD)
         END DO

         CALL DGSDETNXLOC2(CV_SNLOC,SBCVNGI, &
                 XSL,YSL,ZSL, &
                 SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SDETWE,SAREA, &
                 (NDIM==1), (NDIM==3), (NDIM==-2), &
                 SNORMXN,SNORMYN,SNORMZN, &
                 NORMX,NORMY,NORMZ)

         DO CV_SILOC=1,CV_SNLOC
            CV_ILOC=CV_SLOC2LOC( CV_SILOC )
            CV_INOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
            DO CV_SJLOC=1,CV_SNLOC
               CV_JLOC=CV_SLOC2LOC( CV_SJLOC )
               DG_JNOD=(ELE-1)*CV_NLOC+CV_JLOC

               SXNN=0.0
               SYNN=0.0
               SZNN=0.0
               DO SGI=1,SBCVNGI
                  RR=SBCVFEN(CV_SILOC,SGI)*SBCVFEN(CV_SJLOC,SGI)*SDETWE(SGI)
                  SXNN=SXNN+SNORMXN(SGI)*RR
                  SYNN=SYNN+SNORMYN(SGI)*RR
                  SZNN=SZNN+SNORMZN(SGI)*RR
               END DO
               IF(SELE.GT.0) THEN
                  RDIST=0.0
                  !RDIST=DISTANCE_FUN(DG_JNOD) 
               ELSE
                  CV_JLOC2=CV_OTHER_LOC( CV_JLOC )
                  DG_JNOD2=(ELE2-1)*CV_NLOC+CV_JLOC2
                  RDIST=0.5*(-DISTANCE_FUN(DG_JNOD)+DISTANCE_FUN(DG_JNOD2))
                  !RDIST=0.5*(DISTANCE_FUN(DG_JNOD)+DISTANCE_FUN(DG_JNOD2))
               ENDIF
               RHS_DG_X(CV_ILOC)=RHS_DG_X(CV_ILOC) +  SXNN*RDIST !ZX28062013
               RHS_DG_Y(CV_ILOC)=RHS_DG_Y(CV_ILOC) +  SYNN*RDIST
               RHS_DG_Z(CV_ILOC)=RHS_DG_Z(CV_ILOC) +  SZNN*RDIST
            END DO
         END DO

      END DO Between_Elements_And_Boundary 

      SOL_X=0.
      SOL_Y=0.
      SOL_Z=0.
      GOTDEC = .FALSE.
      CALL SMLINNGOT( MASS, SOL_X, RHS_DG_X, CV_NLOC, IPIV, GOTDEC)
      GOTDEC=.TRUE.
      IF(NDIM.GE.2) CALL SMLINNGOT( MASS, SOL_Y, RHS_DG_Y, CV_NLOC, IPIV, GOTDEC)
      IF(NDIM.GE.3) CALL SMLINNGOT( MASS, SOL_Z, RHS_DG_Z, CV_NLOC, IPIV, GOTDEC)
!      CALL MATDMATINV( MASS, INV_MASS, CV_NLOC)

!      DO CV_ILOC=1,CV_NLOC
!         DO CV_JLOC=1,CV_NLOC
!            SOL_X(CV_ILOC)=SOL_X(CV_ILOC) +INV_MASS(CV_ILOC,CV_JLOC)*RHS_DG_X(CV_JLOC)
!            SOL_Y(CV_ILOC)=SOL_Y(CV_ILOC) +INV_MASS(CV_ILOC,CV_JLOC)*RHS_DG_Y(CV_JLOC)
!            SOL_Z(CV_ILOC)=SOL_Z(CV_ILOC) +INV_MASS(CV_ILOC,CV_JLOC)*RHS_DG_Z(CV_JLOC)
!         END DO
!       END DO
      if (factor==2 .and. angle>0) then
      Between_Elements_And_Boundary1: DO IFACE = 1, NFACE
         ELE2  = FACE_ELE( IFACE, ELE )
         SELE2 = MAX( 0, - ELE2 )
         SELE  = SELE2
         ELE2  = MAX( 0, + ELE2 )

         CV_SLOC2LOC( : ) = CV_SLOCLIST( IFACE, : )

            ! Form approximate surface normal (NORMX,NORMY,NORMZ)
         CALL DGSIMPLNORM( ELE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, X_NDGLN, &
                 X, Y, Z, NORMX, NORMY, NORMZ )

            ! Recalculate the normal...
         DO CV_SILOC=1,CV_SNLOC
            CV_ILOC=CV_SLOC2LOC(CV_SILOC)
            X_INOD=X_NDGLN((ELE-1)*X_NLOC+CV_ILOC) 
            XSL(CV_SILOC)=X(X_INOD)
            YSL(CV_SILOC)=Y(X_INOD)
            ZSL(CV_SILOC)=Z(X_INOD)
         END DO

         CALL DGSDETNXLOC2(CV_SNLOC,SBCVNGI, &
                 XSL,YSL,ZSL, &
                 SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SDETWE,SAREA, &
                 (NDIM==1), (NDIM==3), (NDIM==-2), &
                 SNORMXN,SNORMYN,SNORMZN, &
                 NORMX,NORMY,NORMZ)

         SXNN=0.0
         SYNN=0.0
         SZNN=0.0
         DO SGI=1,SBCVNGI
            SXNN=SXNN+SNORMXN(SGI)/SBCVNGI
            SYNN=SYNN+SNORMYN(SGI)/SBCVNGI
            SZNN=SZNN+SNORMZN(SGI)/SBCVNGI
         END DO

         IF (ABS(SXNN).LT.1.0E-6) SXNN=0.0
         IF (ABS(SYNN).LT.1.0E-6) SYNN=0.0
         IF (ABS(SZNN).LT.1.0E-6) SZNN=0.0

         SXTT=SXNN*cos(0.5*pi)-SYNN*sin(0.5*pi)
         SYTT=SXNN*sin(0.5*pi)+SYNN*cos(0.5*pi)
         SZTT=SZNN

         DO CV_SILOC=1,CV_SNLOC
            CV_ILOC=CV_SLOC2LOC( CV_SILOC )
            IF(SELE.GT.0) THEN
               SOL_X(CV_ILOC)=SXNN*cos(pi*angle/180.0)+sign(SXTT,SOL_X(CV_ILOC))*sin(pi*angle/180.0)
               SOL_Y(CV_ILOC)=SYNN*cos(pi*angle/180.0)+sign(SYTT,SOL_Y(CV_ILOC))*sin(pi*angle/180.0)
               SOL_Z(CV_ILOC)=SZNN*cos(pi*angle/180.0)+sign(SZTT,SOL_Z(CV_ILOC))*sin(pi*angle/180.0)
            ENDIF
         END DO

      END DO Between_Elements_And_Boundary1 
      end if

      DO CV_ILOC=1,CV_NLOC
         DG_INOD=(ELE-1)*CV_NLOC+CV_ILOC

         SOL_DERIV_X(DG_INOD)=SOL_X(CV_ILOC)
         SOL_DERIV_Y(DG_INOD)=SOL_Y(CV_ILOC)
         SOL_DERIV_Z(DG_INOD)=SOL_Z(CV_ILOC)

      END DO
 

      DEALLOCATE(SOL_X,SOL_Y,SOL_Z)
      DEALLOCATE(RHS_DG_X,RHS_DG_Y,RHS_DG_Z)
      DEALLOCATE(MASS,INV_MASS)
      DEALLOCATE( CV_SLOC2LOC)
      DEALLOCATE( CV_OTHER_LOC)
      DEALLOCATE( XSL)
      DEALLOCATE( YSL)
      DEALLOCATE( ZSL)
      DEALLOCATE( SNORMXN)
      DEALLOCATE( SNORMYN)
      DEALLOCATE( SNORMZN)
      DEALLOCATE( SDETWE)

        
      END SUBROUTINE LOC_1ST_DERIV_XYZ_DG_DERIV

      END SUBROUTINE SURFACE_TENSION_WRAPPER_NEW


end module multi_surface_tension


