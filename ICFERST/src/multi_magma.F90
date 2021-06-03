!    Copyright (C) 2020 Imperial College London and others.
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Affero General Public License
!    as published by the Free Software Foundation,
!    version 3.0 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without seven the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU General Public
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

    if (is_magma) then
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


  subroutine magma_Coupling_generate(series, state, coupling)
    ! Coefficients phi^2/C is generated and stored in coupling
    implicit none
    type( state_type ), dimension(:), intent( inout ) :: state
    !Global variables
    real, dimension(:) :: series
    !Local variables
    integer :: i, stat
    real,dimension(size(series)) :: phi ! porosity series
    real :: d, a, b !grain size, coefficient a b
    real :: mu !liquid viscosity needs to be build later
    real :: low,high !transition points
    real :: H,s !value of the smoothing function and the smoothing factor
    logical :: Test=.false. ! set to true to have uniform Darcy-like c coefficient
    type( tensor_field ), pointer :: t_field !liquid viscosity
    type(coupling_term_coef), intent(in) :: coupling
    integer :: N  !number of items in the series

    real :: scaling ! a temporal fix for the scaling difference between the viscosity in ICFERST and the models
    scaling=1.0    ! the viscosity difference between ICFERST and the model
    N = size(series)
    s= -2 !> transition coefficient of the linking function
    d=coupling%grain_size
    a=coupling%a
    b=coupling%b
    t_field => extract_tensor_field( state(2), 'Viscosity', stat )
    mu=t_field%val(1,1,1)  !currently only consider a constant liquid viscosity

    low=coupling%cut_low
    high=coupling%cut_high

    do i=1, N
      phi(i) = real(i - 1) / (N - 1)
    end do
    if (Test) then
      do i=2, N
        series(i)=d**2/a/mu*phi(i)**b !coupling%a/d**2*mu*phi(i)**(1-coupling%b)*scaling
      end do
    else
      do i=1, N-1
        if (phi(i)<=low) then
          series(i)= d**2/a/mu*phi(i)**b  !coupling%a/d**2*mu*phi(i)**(1-coupling%b)*scaling
        else if (phi(i)>=high) then
          series(i)= d**2/mu*phi(i)**7/(1-phi(i))!1/d**2*mu*phi(i)**(-6)*(1-phi(i))*scaling
        else
          H=exp(s/((phi(i)-low)/(high-low)))/(exp(s/((phi(i)-low)/(high-low)))+exp(s/(1-(phi(i)-low)/(high-low))))
          series(i)= d**2/mu*phi(i)**7/(1-phi(i))*H+d**2/a/mu*phi(i)**b*(1-H) !(coupling%a/d**2*mu*phi(i)**(1-coupling%b)*(1-H)+1/d**2*mu*phi(i)**(-6)*(1-phi(i))*H)*scaling
        end if
      end do
    end if
    series(N)=2*series(N-1)-series(N-2)
  end subroutine magma_Coupling_generate


  !>@brief: Computes the bulk Composition for two phases
  subroutine cal_bulkcomposition(state,packed_state, Mdims, BulkComposition)
    type( state_type ), dimension( : ), intent( inout ) :: state
    type( state_type ), intent( inout ) :: packed_state
    type(multi_dimensions), intent( in ) :: Mdims
    real, dimension(:), intent(inout) :: BulkComposition
    !Local variables
    type( tensor_field), pointer :: Composition, saturation
    integer :: iphase, cv_inod
    Composition=>extract_tensor_field(packed_state,"PackedConcentration")
    saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")

    BulkComposition = 0.
    do cv_inod = 1, Mdims%cv_nonods
      do iphase = 1, Mdims%nphase
        BulkComposition(cv_inod) = BulkComposition(cv_inod) + Composition%val(1,iphase,cv_inod)*saturation%val(1,iphase,cv_inod)
      end do
    end do
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
    type( scalar_field), pointer :: Cp !bulk composition and heat capacity
    real, dimension(Mdims%cv_nonods) :: BulkComposition
    type(magma_phase_diagram) :: phase_coef
    integer :: cv_nodi

    ! BulkComposition=> extract_scalar_field(state(1), "BulkComposition")
    Composition=>extract_tensor_field(packed_state,"PackedConcentration")
    temperature =>  extract_tensor_field( packed_state, "PackedTemperature" )
    enthalpy=>extract_tensor_field(packed_state,"PackedEnthalpy")
    Cp => extract_scalar_field(state(1), 'TemperatureHeatCapacity')   !
    saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")

    !Compute BulkComposition
    call cal_bulkcomposition(state,packed_state, Mdims, BulkComposition)

    do cv_nodi = 1, Mdims%cv_nonods
      ! above liquidus
        if (enthalpy%val(1,1,cv_nodi)>get_Enthalpy_Liquidus(BulkComposition(cv_nodi), node_val(Cp,cv_nodi) , 0., phase_coef)) then   !currently density is not used, so simply pass 0
          if(BulkComposition(cv_nodi)>phase_coef%Ae) then
            Composition%val(1,1,cv_nodi)=1.0 !solid
          else
            Composition%val(1,1,cv_nodi)=0.0
          end if
          Composition%Val(1,2,cv_nodi)=BulkComposition(cv_nodi)!fluid
            ! below solidus
        else if(enthalpy%val(1,1,cv_nodi)<get_Enthalpy_Solidus(BulkComposition(cv_nodi), node_val(Cp,cv_nodi) , 0., phase_coef)) then
              Composition%val(1,1,cv_nodi)=BulkComposition(cv_nodi)
              Composition%val(1,2,cv_nodi)=phase_coef%Ae
          ! between solidus and eutectic melting line
        else if ((enthalpy%val(1,1,cv_nodi)<get_Enthalpy_Solidus(BulkComposition(cv_nodi), node_val(Cp,cv_nodi) , 0., phase_coef)+&
                   phase_coef%Lf*BulkComposition(cv_nodi)/phase_coef%Ae) .and. &
                   (enthalpy%val(1,1,cv_nodi)<get_Enthalpy_Solidus(BulkComposition(cv_nodi), node_val(Cp,cv_nodi) , 0., phase_coef)+ &
                   phase_coef%Lf*(1-BulkComposition(cv_nodi))/(1.0000001-phase_coef%Ae))) then
              Composition%val(1,2,cv_nodi)=phase_coef%Ae
              Composition%val(1,1,cv_nodi)=(BulkComposition(cv_nodi)-saturation%val(1,2,cv_nodi)*phase_coef%Ae)/saturation%val(1,1,cv_nodi)
      ! between eutectic melting and liquidus
        else
          if(BulkComposition(cv_nodi)>phase_coef%Ae) then
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

  !> @brief: Subroutine to obtain the temperature field for a given enthalpy field
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
  type( tensor_field ), pointer :: enthalpy, temperature,  saturation, dCp , den
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


  !> @brief: Subroutine to obtain the enthalpy field for a given temperature field
  !> Important for the initialisation of the enthalpy, we need ENthalpyOld to start!
  subroutine temperature_to_enthalpy(Mdims, state, packed_state, phase_coef)
  implicit none
  !Global variables
  type( state_type ), dimension( : ), intent( inout ) :: state
  type( state_type ), intent( inout ) :: packed_state
  type(multi_dimensions), intent( in ) :: Mdims
  type(magma_phase_diagram) :: phase_coef

  !Local variables
  integer :: cv_nodi, iphase
  type( tensor_field ), pointer :: enthalpy, temperature,  saturation, rhoCp
  real, dimension(Mdims%cv_nonods) :: enthalpy_dim
  !real, parameter :: tol = 1e-5
    !Temporary until deciding if creating a Cp in packed_state as well
    rhoCp =>extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )
    enthalpy => extract_tensor_field( packed_state,"PackedEnthalpy" )
    temperature =>  extract_tensor_field( packed_state, "PackedTemperature" )
    saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")

    !Calculate temperature using the generic formula H = T_alpha * rho_alpha * Cp_alpha * Saturation_alpha + Lf * phi (phi = Saturation_2)
    enthalpy%val=0.0
    do cv_nodi = 1, Mdims%cv_nonods
      do iphase = 1, Mdims%nphase
        !First enthalpy stored in each phase
        enthalpy%val(1,1,cv_nodi)=enthalpy%val(1,1,cv_nodi)+ temperature%val(1,1,cv_nodi) * rhoCp%val(1,iphase, cv_nodi) * saturation%val(1,iphase,cv_nodi) ! sprint_to_do Need to fix for variable heat capacity.
      end do
      !Now consider the latent heat
      enthalpy%val(1,1,cv_nodi)= enthalpy%val(1,1,cv_nodi) + phase_coef%Lf*saturation%val(1,2,cv_nodi)
    end do
    !Now for consistency populate the temperature of all the phase
    do iphase = 2, Mdims%nphase
      enthalpy%val(1,iphase,:) = enthalpy%val(1,1,:)
    end do
  end subroutine temperature_to_enthalpy


  !> @brief: Compute porosity(Saturation) given an enthalpy field and a BulkComposition
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
  type (scalar_field), pointer :: Cp
  real, dimension(Mdims%cv_nonods) :: BulkComposition
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

    !Compute BulkComposition
    call cal_bulkcomposition(state,packed_state, Mdims, BulkComposition)

    iphase = 1

    do cv_nodi = 1, Mdims%cv_nonods
      Loc_Cp = node_val(Cp,cv_nodi) ! Cp%val(cv_nodi)
      rho=den%val(1,iphase,cv_nodi)
      IF (enthalpy%val(1,1,cv_nodi)>get_Enthalpy_Liquidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef)) then
        saturation%val(1,2, cv_nodi)=1.
        saturation%val(1,1, cv_nodi)=0.
      ELSE IF (enthalpy%val(1,1,cv_nodi)<get_Enthalpy_Solidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef)) then
        saturation%val(1,2, cv_nodi)=0.
        saturation%val(1,1, cv_nodi)=1.
      ELSE
        if (BulkComposition(cv_nodi) <= phase_coef%Ae) then
          He=get_Enthalpy_Solidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef)+phase_coef%Lf*BulkComposition(cv_nodi)/phase_coef%Ae
          if (enthalpy%val(1,1,cv_nodi)<=He) then
            test_poro=BulkComposition(cv_nodi)*(enthalpy%val(1,1,cv_nodi)-get_Enthalpy_Solidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef))/(He-get_Enthalpy_Solidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef))/phase_coef%Ae
          else
            k=1
            test_poro = 1.0
            test_poro_prev=0
            do while (k<Max_its .and. abs(test_poro - test_poro_prev) > tol)
              fx = (phase_coef%Lf/Loc_Cp)*test_poro**3 + (phase_coef%C1-enthalpy%val(1,iphase,cv_nodi)/Loc_Cp)*test_poro**2 + phase_coef%B1*BulkComposition(cv_nodi)*test_poro + phase_coef%A1*BulkComposition(cv_nodi)**2
              fdashx=3*phase_coef%Lf/Loc_Cp*test_poro**2+2*(phase_coef%C1-enthalpy%val(1,iphase,cv_nodi)/Loc_Cp)*test_poro+phase_coef%B1*BulkComposition(cv_nodi)
              test_poro_prev = test_poro
              test_poro = test_poro - fx/fdashx
              k=k+1
            end do
          end if
        else if (BulkComposition(cv_nodi) > phase_coef%Ae) then
          He=get_Enthalpy_Solidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef)+phase_coef%Lf*(1-BulkComposition(cv_nodi))/(1-phase_coef%Ae)
          if (enthalpy%val(1,1,cv_nodi)<=He) then
            test_poro=(1-BulkComposition(cv_nodi))*(enthalpy%val(1,1,cv_nodi)-get_Enthalpy_Solidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef))/(He-get_Enthalpy_Solidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef))/(1-phase_coef%Ae)
          else
            k=1
            test_poro = 1.0
            test_poro_prev=0
            do while (k<Max_its .and. abs(test_poro - test_poro_prev) > tol .and. test_poro<=1)
              fx= (phase_coef%Lf/Loc_Cp)**test_poro**3+(phase_coef%C2-enthalpy%val(1,iphase,cv_nodi)/Loc_Cp+phase_coef%A2+phase_coef%B2)*test_poro**2+(2*phase_coef%A2*BulkComposition(cv_nodi)-2*phase_coef%A2+&
              phase_coef%B2*BulkComposition(cv_nodi)-phase_coef%B2)*test_poro+phase_coef%A2*(BulkComposition(cv_nodi)-1)**2

              fdashx = 3.*(phase_coef%Lf/Loc_Cp)*test_poro**2 + 2.*(phase_coef%C2-enthalpy%val(1,iphase, cv_nodi)/Loc_Cp + phase_coef%A2 + phase_coef%B2)*test_poro +&
              (2.*phase_coef%A2*BulkComposition(cv_nodi) - 2.*phase_coef%A2 + phase_coef%B2*BulkComposition(cv_nodi) - phase_coef%B2)

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


  !>@brief:This subroutine updated the FEM-stored values of the coefficient phi/C in the field magma_absorp
  subroutine update_magma_coupling_coefficients(Mdims, state, saturation, ndgln, Magma_absorp,  c_phi_series)
    implicit none
    type( state_type ), dimension( : ), intent( inout ) :: state
    real, dimension(:,:,:), intent(in) :: saturation
    type(multi_ndgln), intent(in) :: ndgln
    real, dimension(:,:,:,:), INTENT(INOUT) :: Magma_absorp
    type( multi_dimensions ), intent( in ) :: Mdims
    type(coupling_term_coef) :: coupling
    real, dimension(:), intent(in) :: c_phi_series !generated c coefficients
    !Local variables
    integer :: mat_nod, ele, CV_ILOC, cv_inod, iphase, jphase
    real :: magma_coupling
    integer:: c_phi_size ! length of c_phi_series
    real, dimension(4):: test
    c_phi_size=size(c_phi_series)
    ! print *, 'phi2/c', phi2_over_c(saturation%val(1,2, 10))
    DO ELE = 1, Mdims%totele
      DO CV_ILOC = 1, Mdims%cv_nloc
        mat_nod = ndgln%mat( ( ELE - 1 ) * Mdims%mat_nloc + CV_ILOC )
        cv_inod = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
          Do iphase =1, Mdims%nphase-1
            Magma_absorp(1, 1, iphase, mat_nod) = phi2_over_c(saturation(1,2, cv_inod))
          end Do
      end DO
    end DO
  contains
    real function phi2_over_c(phi)
        real, intent(in) :: phi
        integer :: pos
        real:: portion

        pos= int(phi*(c_phi_size-1))+1
        if (pos==c_phi_size) then
          phi2_over_c=c_phi_series(c_phi_size)
        else
          portion=(phi-(pos-1.0)/(c_phi_size-1.0))*c_phi_size
          phi2_over_c=c_phi_series(pos)*(1-portion)+c_phi_series(pos+1)*portion
          ! c_value=c_phi_series(pos)
        end if
      end function phi2_over_c
  end subroutine update_magma_coupling_coefficients


  !>@brief:Compute the source/sink term of the phase change between the concentration living in both two phases
  !> it requires to have Compostion_temp and melt_temp stored before calling cal_solidfluidcomposition
  !> then the RHS is updated to be used in the next non-linear iteration
  subroutine compute_composition_change_source(Mdims, state, packed_state, melt_temp, Compostion_temp, dt)
    type( state_type ), dimension(:), intent( inout ) :: state
    type( state_type ), intent( inout ) :: packed_state
    type( multi_dimensions), intent( in ) :: Mdims
    real, dimension(:), intent(in) :: Compostion_temp, melt_temp
    real, intent(in) :: dt
    !Local variables
    integer :: cv_inod
    type(tensor_field), pointer :: Solute_new, saturation_new
    type(scalar_field), pointer :: Composition_source

    Composition_source => extract_scalar_field(state(1), "Magma_comp_source")
    Solute_new=>extract_tensor_field(packed_state,"PackedConcentration")
    saturation_new => extract_tensor_field(packed_state, "PackedPhaseVolumeFraction")

    do cv_inod = 1, Mdims%cv_nonods
      ! The gain of the first phase  is the loss of the second phase
      Composition_source%val(cv_inod) = -(Solute_new%val(1,2,cv_inod)*saturation_new%val(1,2,cv_inod) - Compostion_temp(cv_inod)*melt_temp(cv_inod))/DT  ! need to check the sign!
    end do
    ! Composition_source(1, :)=-Composition_source(2, :)


  end subroutine compute_composition_change_source

end module multi_magma
