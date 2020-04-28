
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

module multi_magma

    use fldebug
    use futils
    use spud
    use fields
    use global_parameters, only: OPTION_PATH_LEN, PYTHON_FUNC_LEN, PI, is_porous_media
    use vector_tools
    use state_module
    use fields
    use multi_data_types
    ! Hs- Solidus (enthalpy)
    ! Hl - Liquidus (enthalpy)
    ! Ts- Solidus (temperature)
    ! Tl - Liquidus (temperature)
    ! Lf - latent heat
    !Parameters waiting to be decided if they will become fields or inputs from diamond as scalars
    !I would rather have some more descriptive names for these parameters
    !real, parameter :: Hs = 0.,Hl = 1e5,Ts = 1123.15,Tl = 1e3, Lf = 550000.0
    ! Ae - eutetic
    ! A1,B1,C1,A2,B2,C2 - phase behaviour parameters
    !real, parameter :: A1= 50.,B1= -360,C1= 1433.15,A2= 0.,B2= 0.,C2 = 0., Ae = 1.0

    type magma_phase_diagram
      real :: A1!> Phase behaviour parameters
      real :: B1!> Phase behaviour parameters
      real :: C1!> Phase behaviour parameters
      real :: A2!> Phase behaviour parameters
      real :: B2!> Phase behaviour parameters
      real :: C2!> Phase behaviour parameters
      real :: Ae!> Eutectic point
      real :: Ts!> Solidus, uniform solidus si the liquidus at the eutectic point
      real :: Lf!> Latent heat
    end type magma_phase_diagram
    !------------------------------------------------------------------------------
    !>  The copling_term follows the Darcy permeability law C=1/a/d^2*mu*phi^(2-b)
    !>  for lower meltfraction and hindled settling C=1/d^2*mu*phi^(-5)*(1-phi).
    !>  Here coefficients a, grain size d, coefficients b and the cutting between the two demains needs to be set.
    !-----------------------------------------------------------------------------
    type coupling_term_coef
      real :: a
      real :: b
      real :: grain_size
      real :: cut_low   !>  the range below which the coupling terms follow the Darcy permeability, set this to 1 to make the entire domain Darcy like.
      real :: cut_high !>   the range above which the coupling terms follow the Hindered settling
    end type coupling_term_coef

contains
  subroutine initialize_magma_parameters(phase_coef, coupling)
    implicit none
    type(magma_phase_diagram) :: phase_coef
    type(coupling_term_coef) :: coupling

    if (has_phase_diagram) then
      call get_option('/magma_parameters/Phase_diagram_coefficients/A1' , phase_coef%A1 )
      call get_option('/magma_parameters/Phase_diagram_coefficients/B1' , phase_coef%B1 )
      call get_option('/magma_parameters/Phase_diagram_coefficients/C1' , phase_coef%C1 )
      call get_option('/magma_parameters/Phase_diagram_coefficients/ae' , phase_coef%Ae )

      if (have_option('/magma_parameters/Phase_diagram_coefficients/A2')) then
        call get_option('/magma_parameters/Phase_diagram_coefficients/A2' , phase_coef%A2 )
        call get_option('/magma_parameters/Phase_diagram_coefficients/B2' , phase_coef%B2 )
        call get_option('/magma_parameters/Phase_diagram_coefficients/C2' , phase_coef%C2 )
      end if
      call get_option('/magma_parameters/Phase_diagram_coefficients/latent_heat',phase_coef%Lf)
    end if

    call get_option('/magma_parameters/coupling_term_coefficients/coupling_power_coefficent_a' , coupling%a)
    call get_option('/magma_parameters/coupling_term_coefficients/coupling_power_coefficent_b' , coupling%b)
    call get_option('/magma_parameters/coupling_term_coefficients/grain_size' , coupling%grain_size)
    call get_option('/magma_parameters/coupling_term_coefficients/cut_a' , coupling%cut_low)
    call get_option('/magma_parameters/coupling_term_coefficients/cut_b' , coupling%cut_high)

    coupling%cut_high=max(coupling%cut_low,coupling%cut_high)

    phase_coef%Ts=get_Liquidus(phase_coef%Ae,phase_coef)  !Solidus is the liquidus at the ae



  end subroutine initialize_magma_parameters


  subroutine C_generate(series, N,   state, coupling)
    implicit none
    type( state_type ), dimension(:), intent( inout ) :: state
    !Global variables
    real, dimension(:) :: series
    integer :: N  !number of items in the series
    !Local variables
    integer :: i, stat
    real,dimension(N) :: phi ! porosity series
    real :: d !grain size
    real :: mu !liquid viscosity needs to be build later
    real :: low,high !transition points
    real :: H,s !value of the smoothing function and the smoothing factor
    logical :: Test=.true. ! set to true to have uniform Darcy-like c coefficient
    type( tensor_field ), pointer :: t_field !liquid viscosity
    type(coupling_term_coef), intent(in) :: coupling

    s= -2 !> transition coefficient of the linking function
    ! d=35e-6
    d=coupling%grain_size
    t_field => extract_tensor_field( state(2), 'Viscosity', stat )
    mu=t_field%val( 1, 1, 1) !only consider a constant mu for now
    ! mu=1e2
    low=coupling%cut_low
    high=coupling%cut_high

    do i=1, N
      phi(i) = real(i - 1) / (N - 1)
    end do

    if (Test) then
      do i=2, N
        series(i)= 1/coupling%a/d**2*mu*phi(i)**(2-coupling%b)
      end do
    else
      do i=2, N
        if (phi(i)<=low) then
          series(i)= 1/coupling%a/d**2*mu*phi(i)**(2-coupling%b)
        else if (phi(i)>=high) then
          series(i)= 1/d**2*mu*phi(i)**(-5)*(1-phi(i))
        else
          H=exp(s/((phi(i)-low)/(high-low)))/(exp(s/((phi(i)-low)/(high-low)))+exp(s/(1-(phi(i)-low)/(high-low))))
          series(i)=1/coupling%a/d**2*mu*phi(i)**(2-coupling%b)*(1-H)+1/d**2*mu*phi(i)**(-5)*(1-phi(i))*H
        end if
      end do
    end if
    series(1)=2*series(2)-series(3)
  end subroutine C_generate

  subroutine cal_bulkcomposition(state,packed_state)
    type( state_type ), dimension( : ), intent( inout ) :: state
    type( state_type ), intent( inout ) :: packed_state
    type( tensor_field), pointer :: Composition, saturation
    type( scalar_field), pointer :: BulkComposition
    Composition=>extract_tensor_field(packed_state,"PackedComposition")
    saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
    BulkComposition=> extract_scalar_field(state(1),"BulkComposition")

    BulkComposition%val=Composition%val(1,1,:)*saturation%val(1,1,:)+Composition%val(1,2,:)*saturation%val(1,2,:)
  end subroutine

  real function get_Liquidus(Bulk_comp,phase_coef)
    implicit none
    real, intent(in) :: Bulk_comp
    type(magma_phase_diagram) :: phase_coef

    if (Bulk_comp > phase_coef%Ae) then
      get_Liquidus = phase_coef%A2*(Bulk_comp**2) + phase_coef%B2*Bulk_comp + phase_coef%C2
    else
      get_Liquidus = phase_coef%A1*(Bulk_comp**2) + phase_coef%B1*Bulk_comp + phase_coef%C1
    end if

  end function

  real function get_Solidus(Bulk_comp, phase_coef)
    implicit none
    real, intent(in) :: Bulk_comp
    type(magma_phase_diagram) :: phase_coef
    !Currently we consider a flat line
    get_Solidus = phase_coef%Ts
  end function

  ! Function for the liquidus (Enthalpy)
  real function get_Enthalpy_Liquidus(Bulk_comp, Cp, rho, phase_coef)
    implicit none
    real, intent(in) :: Bulk_comp, Cp, rho
    type(magma_phase_diagram) :: phase_coef
    get_Enthalpy_Liquidus = get_Liquidus(Bulk_comp, phase_coef) * Cp  + phase_coef%Lf
  end function

  ! Function for the Solidus (Enthalpy)
  real function get_Enthalpy_Solidus(Bulk_comp, Cp, rho, phase_coef)
  implicit none
  real, intent(in) :: Bulk_comp, Cp, rho
  type(magma_phase_diagram) :: phase_coef
  !Currently we consider a flat line
    get_Enthalpy_Solidus = phase_coef%Ts * Cp !* rho

  end function

  subroutine cal_solidfluidcomposition(state, packed_state, Mdims, phase_coef)
    type( state_type ), dimension( : ), intent( inout ) :: state
    type( state_type ), intent( inout ) :: packed_state
    type( multi_dimensions), intent( in ) :: Mdims
    type( tensor_field), pointer :: Composition, temperature, enthalpy, saturation
    type( scalar_field), pointer :: BC, Cp !bulk composition and heat capacity
    type(magma_phase_diagram) :: phase_coef
    integer :: cv_nodi

    BC=> extract_scalar_field(state(1), "BulkComposition")
    Composition=>extract_tensor_field(packed_state,"PackedComposition")
    temperature =>  extract_tensor_field( packed_state, "PackedTemperature" )
    enthalpy=>extract_tensor_field(packed_state,"PackedEnthalpy")
    Cp => extract_scalar_field(state(1), 'TemperatureHeatCapacity')   !
    saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")

    do cv_nodi = 1, Mdims%cv_nonods
      ! above liquidus
        if (enthalpy%val(1,1,cv_nodi)>get_Enthalpy_Liquidus(BC%val(cv_nodi), node_val(Cp,cv_nodi) , 0., phase_coef)) then   !currently density is not used, so simply pass 0
          if(BC%val(cv_nodi)>phase_coef%Ae) then
            Composition%val(1,1,cv_nodi)=1.0 !solid
          else
            Composition%val(1,1,cv_nodi)=0.0
          end if
          Composition%Val(1,2,cv_nodi)=BC%val(cv_nodi)!fluid
            ! below solidus
        else if(enthalpy%val(1,1,cv_nodi)<get_Enthalpy_Solidus(BC%val(cv_nodi), node_val(Cp,cv_nodi) , 0., phase_coef)) then
              Composition%val(1,1,cv_nodi)=BC%val(cv_nodi)
              Composition%val(1,2,cv_nodi)=phase_coef%Ae
          ! between solidus and eutectic melting line
        else if ((enthalpy%val(1,1,cv_nodi)<get_Enthalpy_Solidus(BC%val(cv_nodi), node_val(Cp,cv_nodi) , 0., phase_coef)+ phase_coef%Lf*BC%val(cv_nodi)/phase_coef%Ae) .and. (enthalpy%val(1,1,cv_nodi)<get_Enthalpy_Solidus(BC%val(cv_nodi), node_val(Cp,cv_nodi) , 0., phase_coef)+ phase_coef%Lf*(1-BC%val(cv_nodi))/(1.0000001-phase_coef%Ae))) then
              Composition%val(1,2,cv_nodi)=phase_coef%Ae
              Composition%val(1,1,cv_nodi)=(BC%val(cv_nodi)-saturation%val(1,2,cv_nodi)*phase_coef%Ae)/saturation%val(1,1,cv_nodi)
      ! between eutectic melting and liquidus
        else
          if(BC%val(cv_nodi)>phase_coef%Ae) then
            Composition%val(1,1,cv_nodi)=1.0
            ! now that A2=B2=0, this expression is meaningless for fluid compositino and will get a error for devide A2 (=0)
            Composition%val(1,2,cv_nodi)= 1.0
            !FluidComposition%val(1,1,cv_nodi)=(-B2 + sqrt(( B1**2 ) - 4. * A2 * ( C2 - temperature%val(1,1,cv_nodi)))) / ( 2. * A2)
          else
            Composition%val(1,1,cv_nodi)=0.0
            Composition%val(1,2,cv_nodi)=(-phase_coef%B1 - sqrt(phase_coef%B1**2  - 4. * phase_coef%A1 * ( phase_coef%C1 - temperature%val(1,1,cv_nodi)))) / ( 2. * phase_coef%A1)
          end if
        end if
    end do
  end subroutine

  !Subroutine to obtain the temperature field for a given enthalpy field
  subroutine enthalpy_to_temperature(Mdims, state, packed_state, phase_coef)
  implicit none
  !Global variables
  type( state_type ), dimension( : ), intent( inout ) :: state
  type( state_type ), intent( inout ) :: packed_state
  type(multi_dimensions), intent( in ) :: Mdims
  type(magma_phase_diagram) :: phase_coef

  !Local variables
  integer :: cv_nodi, iphase
  type ( scalar_field), pointer :: Cp                       !Temporary until deciding if creating a Cp in packed_state as well
  type( tensor_field ), pointer :: enthalpy, temperature,  saturation, FluidComposition, dCp , den
  real, dimension(Mdims%cv_nonods) :: enthalpy_dim
  !real, parameter :: tol = 1e-5
    !Temporary until deciding if creating a Cp in packed_state as well
    den => extract_tensor_field( packed_state,"PackedDensity" )
    Cp => extract_scalar_field( state(1), 'TemperatureHeatCapacity')
    dCp =>extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )
    enthalpy => extract_tensor_field( packed_state,"PackedEnthalpy" )
    temperature =>  extract_tensor_field( packed_state, "PackedTemperature" )
    saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
    !saturation =>  extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )!First phase is rock presence

!IT CONSIDERS ONE SINGLE TEMPERATURE!!! WHICH MAKE SENSE BECAUSE AT THAT TIME SCALE THINGS ARE IN EQUILIBRIUM...
    ! do iphase =1, nphase
    !Calculate temperature using the generic formula T = (H - Lf*porosity)/Cp
    temperature%val(1,1,:)=(enthalpy%val(1,1,:)- phase_coef%Lf*saturation%val(1,2,:))/Cp%val(1)  ! sprint_to_do Need to fix for variable heat capacity.
    !Now add corrections if required
    where (saturation%val(1,2, :)>0.) !If porosity > 0.
      where (temperature%val(1,1,:) < phase_coef%Ts ) !If there is mixture, temperature cannot be below the solidus
        temperature%val(1,1,:) = phase_coef%Ts!It seems that the Solidus line is flat
      end where
    ! elsewhere  !not necesarry
    !   temperature%val(1,1,:) = enthalpy_dim* den%val(1,iphase,:)/Density_Cp%val(1,iphase,:)
    end where
    !Now for consistency populate the temperature of all the phase
    do iphase = 2, Mdims%ndim
      temperature%val(1,iphase,:) = temperature%val(1,1,:)
    end do
  end subroutine enthalpy_to_temperature

  !Compute porosity given an enthalpy field and a BulkComposition
  subroutine porossolve(state, packed_state, Mdims, ndgln, phase_coef)
  implicit none
  !Global variables
  type(multi_dimensions), intent( in ) :: Mdims
  type( state_type ), dimension( : ), intent( inout ) :: state
  type( state_type ), intent( inout ) :: packed_state
  type(multi_ndgln), intent(in) :: ndgln
  type(magma_phase_diagram) :: phase_coef
  !Local variables
  integer ::  iphase, k, i
  real :: test_poro, test_poro_prev        !Temporary until deciding if creating a Cp in packed_state as well
  type( tensor_field ), pointer :: enthalpy, den, saturation
  type (scalar_field), pointer :: BC, Cp
  !real, dimension(Mdims%cv_nonods) :: enthalpy_dim
  real :: fx, fdashx, Loc_Cp, rho
  !Parameters for the non_linear solvers (Maybe a newton solver here makes sense?)
  real, parameter :: tol = 1e-6 !Need to be at least 1e-5 to obtain a relative stable result
  integer, parameter :: max_its = 25
  !!
  integer :: CV_ILOC, cv_nodi
  real :: ELE, He
    !Temporary until deciding if creating a Cp in packed_state as well
    den => extract_tensor_field( packed_state,"PackedDensity" )
    Cp => extract_scalar_field( state(1), 'TemperatureHeatCapacity')
    enthalpy => extract_tensor_field( packed_state,"PackedEnthalpy" )
    saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
    BC=> extract_scalar_field(state(1), "BulkComposition")
    iphase = 1

    do cv_nodi = 1, Mdims%cv_nonods
      Loc_Cp = node_val(Cp,cv_nodi) ! Cp%val(cv_nodi)
      rho=den%val(1,iphase,cv_nodi)
      IF (enthalpy%val(1,1,cv_nodi)>get_Enthalpy_Liquidus(BC%val(cv_nodi),Loc_Cp, rho, phase_coef)) then
        saturation%val(1,2, cv_nodi)=1.
        saturation%val(1,1, cv_nodi)=0.
      ELSE IF (enthalpy%val(1,1,cv_nodi)<get_Enthalpy_Solidus(BC%val(cv_nodi),Loc_Cp, rho, phase_coef)) then
        saturation%val(1,2, cv_nodi)=0.
        saturation%val(1,1, cv_nodi)=1.
      ELSE
        if (BC%val(cv_nodi) <= phase_coef%Ae) then
          He=get_Enthalpy_Solidus(BC%val(cv_nodi),Loc_Cp, rho, phase_coef)+phase_coef%Lf*BC%val(cv_nodi)/phase_coef%Ae
          if (enthalpy%val(1,1,cv_nodi)<=He) then
            test_poro=BC%val(cv_nodi)*(enthalpy%val(1,1,cv_nodi)-get_Enthalpy_Solidus(BC%val(cv_nodi),Loc_Cp, rho, phase_coef))/(He-get_Enthalpy_Solidus(BC%val(cv_nodi),Loc_Cp, rho, phase_coef))/phase_coef%Ae
          else
            k=1
            test_poro = 1.0
            test_poro_prev=0
            do while (k<Max_its .and. abs(test_poro - test_poro_prev) > tol)
              fx = (phase_coef%Lf/Loc_Cp)*test_poro**3 + (phase_coef%C1-enthalpy%val(1,iphase,cv_nodi)/Loc_Cp)*test_poro**2 + phase_coef%B1*BC%val(cv_nodi)*test_poro + phase_coef%A1*BC%val(cv_nodi)**2
              fdashx=3*phase_coef%Lf/Loc_Cp*test_poro**2+2*(phase_coef%C1-enthalpy%val(1,iphase,cv_nodi)/Loc_Cp)*test_poro+phase_coef%B1*BC%val(cv_nodi)
              test_poro_prev = test_poro
              test_poro = test_poro - fx/fdashx
              k=k+1
            end do
          end if
        else if (BC%val(cv_nodi) > phase_coef%Ae) then
          He=get_Enthalpy_Solidus(BC%val(cv_nodi),Loc_Cp, rho, phase_coef)+phase_coef%Lf*(1-BC%val(cv_nodi))/(1-phase_coef%Ae)
          if (enthalpy%val(1,1,cv_nodi)<=He) then
            test_poro=(1-BC%val(cv_nodi))*(enthalpy%val(1,1,cv_nodi)-get_Enthalpy_Solidus(BC%val(cv_nodi),Loc_Cp, rho, phase_coef))/(He-get_Enthalpy_Solidus(BC%val(cv_nodi),Loc_Cp, rho, phase_coef))/(1-phase_coef%Ae)
          else
            k=1
            test_poro = 1.0
            test_poro_prev=0
            do while (k<Max_its .and. abs(test_poro - test_poro_prev) > tol .and. test_poro<=1)
              fx= (phase_coef%Lf/Loc_Cp)**test_poro**3+(phase_coef%C2-enthalpy%val(1,iphase,cv_nodi)/Loc_Cp+phase_coef%A2+phase_coef%B2)*test_poro**2+(2*phase_coef%A2*BC%val(cv_nodi)-2*phase_coef%A2+&
              phase_coef%B2*BC%val(cv_nodi)-phase_coef%B2)*test_poro+phase_coef%A2*(BC%val(cv_nodi)-1)**2

              fdashx = 3.*(phase_coef%Lf/Loc_Cp)*test_poro**2 + 2.*(phase_coef%C2-enthalpy%val(1,iphase, cv_nodi)/Loc_Cp + phase_coef%A2 + phase_coef%B2)*test_poro +&
              (2.*phase_coef%A2*BC%val(cv_nodi) - 2.*phase_coef%A2 + phase_coef%B2*BC%val(cv_nodi) - phase_coef%B2)

              test_poro_prev = test_poro!Store to check convergence the previous value
              test_poro = test_poro - fx/fdashx
              k=k+1
            end do
          end if
        end if
        saturation%val(1,2, cv_nodi)=test_poro
        saturation%val(1,1, cv_nodi)=1-test_poro
      END IF
    end do
  end subroutine

end module multi_magma
