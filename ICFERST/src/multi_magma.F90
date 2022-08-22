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
    use boundary_conditions, only: get_entire_boundary_condition
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
    ! Coefficients phi^2/C (C without liquid viscosity) is generated and stored in coupling
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
    integer :: index_fluid  !number of items in the series
    real,  PARAMETER :: phi_min=1e-8, cap_suspension=0.6
    real :: scaling ! a temporal fix for the scaling difference between the viscosity in ICFERST and the models
    real :: suspension_scale=20.0   !HH
    scaling=1.0    ! the viscosity difference between ICFERST and the model
    N = size(series)
    s= -2 !> transition coefficient of the linking function

    d=coupling%grain_size
    a=coupling%a
    b=coupling%b
    t_field => extract_tensor_field( state(2), 'Viscosity', stat )
    mu=1.!

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
      do i=2, N
        if (phi(i)<=phi_min) then
          series(i)= phi(i)**2/(a*mu/d**2*phi_min**(2-b)+(2-b)*a*mu/d**2*phi_min**(1-b)*(phi(i)-phi_min))
        else if (phi(i)<=low) then
          series(i)= d**2/a/mu*phi(i)**b  !coupling%a/d**2*mu*phi(i)**(1-coupling%b)*scaling
        else if (phi(i)>=high) then
          series(i)= d**2/mu*phi(i)**7/(1-phi(i))/suspension_scale!1/d**2*mu*phi(i)**(-6)*(1-phi(i))*scaling
        else
          H=exp(s/((phi(i)-low)/(high-low)))/(exp(s/((phi(i)-low)/(high-low)))+exp(s/(1-(phi(i)-low)/(high-low))))
          series(i)= d**2/mu*phi(i)**7/(1-phi(i))*H/suspension_scale+d**2/a/mu*phi(i)**b*(1-H) !(coupling%a/d**2*mu*phi(i)**(1-coupling%b)*(1-H)+1/d**2*mu*phi(i)**(-6)*(1-phi(i))*H)*scaling
        end if
      end do
    end if

    series(1)=series(2)
    series(N)=series(N-1)

    index_fluid=int(cap_suspension*N)
    if (high<1.0) series(index_fluid:N)=series(index_fluid)
  end subroutine magma_Coupling_generate


  !>@brief: Computes the bulk Composition for two phases
  subroutine cal_bulkcomposition(state,packed_state, Mdims, initilization)
    type( state_type ), dimension( : ), intent( inout ) :: state
    type( state_type ), intent( inout ) :: packed_state
    type(multi_dimensions), intent( in ) :: Mdims
    logical, optional :: initilization
    !Local variables
    type( tensor_field), pointer :: Composition, saturation
    integer :: iphase, cv_inod
    real :: init_bulk
    type( scalar_field ), pointer :: Bcomposition
    Bcomposition=>extract_scalar_field(state(1),'BulkComposition')

    Composition=>extract_tensor_field(packed_state,"PackedConcentration")
    saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")

    ! call get_option('/material_phase::phase1/scalar_field::BulkComposition/prognostic/initial_condition::WholeMesh/constant', init_bulk)
    
    ! if (present_and_true(initilization)) then 
    !   Bcomposition%val=0.2
    ! else
      Bcomposition%val = 0.
      do cv_inod = 1, Mdims%cv_nonods
        do iphase = 1, Mdims%nphase
          Bcomposition%val(cv_inod) = Bcomposition%val(cv_inod) + Composition%val(1,iphase,cv_inod)*saturation%val(1,iphase,cv_inod)
        end do
        Bcomposition%val(cv_inod)=max(Bcomposition%val(cv_inod),0.)
        Bcomposition%val(cv_inod)=min(Bcomposition%val(cv_inod),1.)
      end do

    ! end if 
  end subroutine

  !>@brief: Computes the contribution for the melt fraction and bulk composition change
  subroutine cal_contribution2(state,packed_state, Mdims, Mmat, ndgln, dt)
    type( state_type ), dimension( : ), intent( inout ) :: state
    type( state_type ), intent( inout ) :: packed_state
    type(multi_dimensions), intent( in ) :: Mdims
    type (multi_matrices), intent(inout) :: Mmat
    ! type (multi_sparsities), intent(in) :: Mspars
    type(multi_ndgln), intent(in) :: ndgln
    real, intent( in ) :: dt
    !Local variables
    type( tensor_field), pointer :: Composition, old_Composition, saturation, old_saturation
    integer :: iphase, cv_inod
    real :: init_bulk, weight
    type( scalar_field ), pointer :: Bcomposition, phidecomH, phidecomC,phidecomR,  phidecomH2, phidecomC2,phidecomR2
    ! type( vector_field ) :: rhs_p
    ! REAL, DIMENSION( :, :, : ), allocatable :: scaled_velocity
    ! integer :: ele, u_iloc, u_inod, cv_iloc, cv_loc
    Bcomposition=>extract_scalar_field(state(1),'BulkComposition')

    Composition=>extract_tensor_field(packed_state,"PackedConcentration")
    old_Composition=>extract_tensor_field(packed_state,"PackedOldConcentration")
    saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
    old_saturation=>extract_tensor_field(packed_state,"PackedOldPhaseVolumeFraction")

    

    if (have_option('/material_phase::phase1/scalar_field::PhiDecompositionC')) then
      phidecomC=>extract_scalar_field(state(1),'PhiDecompositionC')
      if (have_option('/material_phase::phase1/scalar_field::PhiDecompositionCa')) then
        phidecomC2=>extract_scalar_field(state(1),'PhiDecompositionCa')
        phidecomC%val=old_saturation%val(1,1,:)-phidecomC%val
        phidecomC2%val=phidecomC2%val+phidecomC%val
      end if

      if (have_option('/material_phase::phase1/scalar_field::PhiDecompositionH')) then 
        phidecomH=>extract_scalar_field(state(1),'PhiDecompositionH')
        weight=0.5
        phidecomH%val=1./((Composition%val(1,2,:)-Composition%val(1,1,:))*weight+(old_Composition%val(1,2,:)-old_Composition%val(1,1,:))*(1-weight))* &
        (-(Composition%val(1,1,:)-old_Composition%val(1,1,:))*(saturation%val(1,1,:)*weight+old_saturation%val(1,1,:)*(1-weight))- (Composition%val(1,2,:)-old_Composition%val(1,2,:))*(saturation%val(1,2,:)*weight+old_saturation%val(1,2,:)*(1-weight)))
      end if

      if (have_option('/material_phase::phase1/scalar_field::PhiDecompositionHa')) then
        phidecomH2=>extract_scalar_field(state(1),'PhiDecompositionHa')
        phidecomH2%val=phidecomH2%val+phidecomH%val
      end if

      if (have_option('/material_phase::phase1/scalar_field::PhiDecompositionR')) then 
        phidecomR=>extract_scalar_field(state(1),'PhiDecompositionR')
        phidecomC=>extract_scalar_field(state(1),'PhiDecompositionC')
        phidecomR%val=saturation%val(1,2,:)-old_saturation%val(1,2,:)-phidecomH%val-phidecomC%val
      end if

      if (have_option('/material_phase::phase1/scalar_field::PhiDecompositionRa')) then 
        phidecomR2=>extract_scalar_field(state(1),'PhiDecompositionRa')
        phidecomR2%val=phidecomR2%val+phidecomR%val
      end if
      print *, 'Done!', phidecomH%val(100), phidecomC%val(100), phidecomR%val(100)
    end if

    ! call get_option('/material_phase::phase1/scalar_field::BulkComposition/prognostic/initial_condition::WholeMesh/constant', init_bulk)
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

    type( scalar_field ), pointer :: Bcomposition
    Composition=>extract_tensor_field(packed_state,"PackedConcentration")
    temperature =>  extract_tensor_field( packed_state, "PackedTemperature" )
    enthalpy=>extract_tensor_field(packed_state,"PackedEnthalpy")
    Cp => extract_scalar_field(state(1), 'TemperatureHeatCapacity')   !
    saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")

   

    ! DO NOT RECALCULATE BULK COMPOSITION HERE SINCE THE MELT FRACTION HAS CHANGED! JUST LOAD THE CALCULATED ONE!!
    !Compute BulkComposition
    Bcomposition=>extract_scalar_field(state(1),'BulkComposition')
    BulkComposition=Bcomposition%val 

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
                    print *, 'eutectic melting'
              Composition%val(1,2,cv_nodi)=phase_coef%Ae
              Composition%val(1,1,cv_nodi)=(BulkComposition(cv_nodi)-saturation%val(1,2,cv_nodi)*phase_coef%Ae)/max(saturation%val(1,1,cv_nodi),1e-8)
              Composition%val(1,1,cv_nodi)=max(Composition%val(1,1,cv_nodi),0.)
              Composition%val(1,1,cv_nodi)=min(Composition%val(1,1,cv_nodi),1.)
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
  type( tensor_field ), pointer :: enthalpy, temperature,  saturation, rhoCp, Den
  real, dimension(Mdims%cv_nonods) :: enthalpy_dim
  !real, parameter :: tol = 1e-5
    !Temporary until deciding if creating a Cp in packed_state as well
    rhoCp =>extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )
    enthalpy => extract_tensor_field( packed_state,"PackedEnthalpy" )
    temperature =>  extract_tensor_field( packed_state, "PackedTemperature" )
    saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
    Den =>extract_tensor_field( packed_state, 'PackedDensity')

    !FYI: IT CONSIDERS ONE SINGLE TEMPERATURE!!! WHICH MAKE SENSE BECAUSE AT THAT TIME SCALE THINGS ARE IN EQUILIBRIUM...
    !Calculate temperature using the generic formula T = (H - Lf*porosity)/Cp
    do cv_nodi = 1, Mdims%cv_nonods
      temperature%val(1,1,cv_nodi)=(enthalpy%val(1,1,cv_nodi)- phase_coef%Lf*saturation%val(1,2,cv_nodi))/rhoCp%val(1, 1, cv_nodi)*Den%val(1,1,cv_nodi)
    end do
    !Now add corrections if required
    where (1. - saturation%val(1,1, :) > 0.) !If porosity > 0.
      where (temperature%val(1,1,:) < phase_coef%Ts ) !If there is mixture, temperature cannot be below the solidus
        temperature%val(1,1,:) = phase_coef%Ts!It seems that the Solidus line is flat
      end where
    ! elsewhere  !not necesarry
    !   temperature%val(1,1,:) = enthalpy_dim* den%val(1,iphase,:)/Density_Cp%val(1,iphase,:)
    end where
    !Now for consistency populate the temperature of all the phase
    do cv_nodi = 1, Mdims%cv_nonods
      do iphase = 2, Mdims%ndim
        temperature%val(1,iphase,cv_nodi) = temperature%val(1,1,cv_nodi)
      end do
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
  type( tensor_field ), pointer :: enthalpy, temperature,  saturation, rhoCp, Den
  ! type( scalar_field ), pointer :: Cp
  real, dimension(Mdims%cv_nonods) :: enthalpy_dim
  rhoCp =>extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )
  Den =>extract_tensor_field( packed_state, 'PackedDensity')
  enthalpy => extract_tensor_field( packed_state,"PackedEnthalpy" )
  temperature =>  extract_tensor_field( packed_state, "PackedTemperature" )
  saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")

  !Calculate temperature using the generic formula H = T_alpha * rho_alpha * Cp_alpha * Saturation_alpha + Lf * phi (phi = 1- Saturation_1)
  enthalpy%val=0.0
  do cv_nodi = 1, Mdims%cv_nonods
    do iphase = 1, Mdims%nphase
      !First enthalpy stored in each phase
      enthalpy%val(1,1,cv_nodi)=enthalpy%val(1,1,cv_nodi) + temperature%val(1,1,cv_nodi) *rhoCp%val(1, iphase, cv_nodi) /Den%val(1,iphase,cv_nodi) * saturation%val(1,iphase,cv_nodi)
    end do
    !Now consider the latent heat
    enthalpy%val(1,1,cv_nodi)= enthalpy%val(1,1,cv_nodi) + phase_coef%Lf*(1.-saturation%val(1,1,cv_nodi))
  end do
  !Now for consistency populate the temperature of all the phase
  do cv_nodi = 1, Mdims%cv_nonods
    do iphase = 2, Mdims%nphase
      enthalpy%val(1,iphase,cv_nodi) = enthalpy%val(1,1,cv_nodi)
    end do
  end do
  end subroutine temperature_to_enthalpy


  !> @brief: Compute porosity(Saturation) given an enthalpy field and a BulkComposition
  subroutine poro_component_solve(state, packed_state, Mdims, ndgln, phase_coef, initilization)
  implicit none
  !Global variables
  type(multi_dimensions), intent( in ) :: Mdims
  type( state_type ), dimension( : ), intent( inout ) :: state
  type( state_type ), intent( inout ) :: packed_state
  type(multi_ndgln), intent(in) :: ndgln
  type(magma_phase_diagram) :: phase_coef
  logical, optional :: initilization
  !Local variables
  integer ::  iphase, k, i
  real :: test_poro, test_poro_prev        !Temporary until deciding if creating a Cp in packed_state as well
  type( tensor_field ), pointer :: enthalpy, den, saturation, Composition
  type( vector_field ), pointer :: p_position
  type( scalar_field ), pointer :: Cp, Bcomposition
  real, dimension(Mdims%cv_nonods) :: BulkComposition
  !real, dimension(Mdims%cv_nonods) :: enthalpy_dim
  real :: fx, fdashx, Loc_Cp, rho
  !Parameters for the non_linear solvers (Maybe a newton solver here makes sense?)
  real, parameter :: phi_min = 1e-7 !Need to be at least 1e-5 to obtain a relative stable result
  integer, parameter :: max_its = 25
  !!
  integer :: CV_ILOC, cv_nodi
  real :: ELE, He
  !!
  logical :: simplified_phase_diagram !Using a simplified diagram for validation of the thermal cases
  real:: p, q
  simplified_phase_diagram=.true.
  if  (simplified_phase_diagram) then 
    p=0.0165
    q=-11
  end if

    !Temporary until deciding if creating a Cp in packed_state as well
    den => extract_tensor_field( packed_state,"PackedDensity" )
    Cp => extract_scalar_field( state(1), 'TemperatureHeatCapacity')
    enthalpy => extract_tensor_field( packed_state,"PackedEnthalpy" )
    saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
    Composition=>extract_tensor_field(packed_state,"PackedConcentration")

    !Compute BulkComposition
    if (present_and_true(initilization)) then 
      print *, 'Setting the initial condition by recalculating the melt fraction with the giving bulk composition.'
      call cal_bulkcomposition(state,packed_state, Mdims, initilization=initilization)
    else
      call cal_bulkcomposition(state,packed_state, Mdims)
    end if

    Bcomposition=>extract_scalar_field(state(1),'BulkComposition')
    BulkComposition=Bcomposition%val 
    iphase = 1

    p_position=>extract_vector_field(packed_state,"PressureCoordinate")
    
      
    do cv_nodi = 1, Mdims%cv_nonods
      if (.not. simplified_phase_diagram) then 
        Loc_Cp = node_val(Cp,cv_nodi) ! Cp%val(cv_nodi)
        rho=den%val(1,iphase,cv_nodi)
        IF (enthalpy%val(1,1,cv_nodi)>get_Enthalpy_Liquidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef)) then
          saturation%val(1,2, cv_nodi)=1.
          saturation%val(1,1, cv_nodi)=0.
          if(BulkComposition(cv_nodi)>phase_coef%Ae) then
            Composition%val(1,1,cv_nodi)=1.0 !solid
          else
            Composition%val(1,1,cv_nodi)=0.0
          end if
          Composition%Val(1,2,cv_nodi)=BulkComposition(cv_nodi)
        ELSE IF (enthalpy%val(1,1,cv_nodi)<get_Enthalpy_Solidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef)) then
          saturation%val(1,2, cv_nodi)=0.
          saturation%val(1,1, cv_nodi)=1.
          Composition%val(1,1,cv_nodi)=BulkComposition(cv_nodi)
          Composition%val(1,2,cv_nodi)=phase_coef%Ae
        ELSE
          if (BulkComposition(cv_nodi) <= phase_coef%Ae) then
            He=get_Enthalpy_Solidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef)+phase_coef%Lf*BulkComposition(cv_nodi)/phase_coef%Ae
            if (enthalpy%val(1,1,cv_nodi)<=He) then
              test_poro=BulkComposition(cv_nodi)*(enthalpy%val(1,1,cv_nodi)-get_Enthalpy_Solidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef))/(He-get_Enthalpy_Solidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef))/phase_coef%Ae
              Composition%val(1,2,cv_nodi)=phase_coef%Ae
              Composition%val(1,1,cv_nodi)=(BulkComposition(cv_nodi)-test_poro*phase_coef%Ae)/max(1-test_poro,1e-8)
              Composition%val(1,1,cv_nodi)=max(Composition%val(1,1,cv_nodi),0.)
              Composition%val(1,1,cv_nodi)=min(Composition%val(1,1,cv_nodi),1.)
            else
              k=1
              test_poro = 1.0
              test_poro_prev=0
              do while (k<Max_its .and. abs(test_poro - test_poro_prev) > phi_min)
                fx = (phase_coef%Lf/Loc_Cp)*test_poro**3 + (phase_coef%C1-enthalpy%val(1,iphase,cv_nodi)/Loc_Cp)*test_poro**2 + phase_coef%B1*BulkComposition(cv_nodi)*test_poro + phase_coef%A1*BulkComposition(cv_nodi)**2
                fdashx=3*phase_coef%Lf/Loc_Cp*test_poro**2+2*(phase_coef%C1-enthalpy%val(1,iphase,cv_nodi)/Loc_Cp)*test_poro+phase_coef%B1*BulkComposition(cv_nodi)
                test_poro_prev = test_poro
                test_poro = test_poro - fx/fdashx
                k=k+1
              end do
              test_poro=max(test_poro, 0.)
              test_poro=min(test_poro, 1.)
              Composition%val(1,1,cv_nodi)=0.
              Composition%val(1,2,cv_nodi)=BulkComposition(cv_nodi)/max(test_poro,1e-8)
            end if
          else if (BulkComposition(cv_nodi) > phase_coef%Ae) then
            He=get_Enthalpy_Solidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef)+phase_coef%Lf*(1-BulkComposition(cv_nodi))/(1-phase_coef%Ae)
            if (enthalpy%val(1,1,cv_nodi)<=He) then
              test_poro=(1-BulkComposition(cv_nodi))*(enthalpy%val(1,1,cv_nodi)-get_Enthalpy_Solidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef))/(He-get_Enthalpy_Solidus(BulkComposition(cv_nodi),Loc_Cp, rho, phase_coef))/(1-phase_coef%Ae)
              Composition%val(1,2,cv_nodi)=phase_coef%Ae
              Composition%val(1,1,cv_nodi)=(BulkComposition(cv_nodi)-test_poro*phase_coef%Ae)/max(1-test_poro,1e-8)
              Composition%val(1,1,cv_nodi)=max(Composition%val(1,1,cv_nodi),0.)
              Composition%val(1,1,cv_nodi)=min(Composition%val(1,1,cv_nodi),1.)
            else
              k=1
              test_poro = 1.0
              test_poro_prev=0
              do while (k<Max_its .and. abs(test_poro - test_poro_prev) > phi_min .and. test_poro<=1)
                fx= (phase_coef%Lf/Loc_Cp)**test_poro**3+(phase_coef%C2-enthalpy%val(1,iphase,cv_nodi)/Loc_Cp+phase_coef%A2+phase_coef%B2)*test_poro**2+(2*phase_coef%A2*BulkComposition(cv_nodi)-2*phase_coef%A2+&
                phase_coef%B2*BulkComposition(cv_nodi)-phase_coef%B2)*test_poro+phase_coef%A2*(BulkComposition(cv_nodi)-1)**2

                fdashx = 3.*(phase_coef%Lf/Loc_Cp)*test_poro**2 + 2.*(phase_coef%C2-enthalpy%val(1,iphase, cv_nodi)/Loc_Cp + phase_coef%A2 + phase_coef%B2)*test_poro +&
                (2.*phase_coef%A2*BulkComposition(cv_nodi) - 2.*phase_coef%A2 + phase_coef%B2*BulkComposition(cv_nodi) - phase_coef%B2)

                test_poro_prev = test_poro!Store to check convergence the previous value
                test_poro = test_poro - fx/fdashx
                k=k+1
              end do
              test_poro=max(test_poro, 0.)
              test_poro=min(test_poro, 1.)
              Composition%val(1,1,cv_nodi)=0.
              Composition%val(1,2,cv_nodi)=BulkComposition(cv_nodi)/max(test_poro,1e-8)
            end if
          end if
          saturation%val(1,2, cv_nodi)=test_poro
          saturation%val(1,1, cv_nodi)=1.-test_poro
        END IF
      else !Simplified phase diagram for thermal validations
        if (cv_nodi==5) print *,enthalpy%val(1,1,cv_nodi),BulkComposition(cv_nodi)
        if (cv_nodi==5) print *,BulkComposition(cv_nodi)*(p*enthalpy%val(1,1,cv_nodi)+q*node_val(Cp,cv_nodi)), (node_val(Cp,cv_nodi)+BulkComposition(cv_nodi)*p*phase_coef%Lf)
        saturation%val(1,2, cv_nodi)=BulkComposition(cv_nodi)*(p*enthalpy%val(1,1,cv_nodi)+q*node_val(Cp,cv_nodi))/(node_val(Cp,cv_nodi)+BulkComposition(cv_nodi)*p*phase_coef%Lf)
        saturation%val(1,1, cv_nodi)=1.-saturation%val(1,2, cv_nodi)
        Composition%val(1,2,cv_nodi)=BulkComposition(cv_nodi)/saturation%val(1,2, cv_nodi)
      end if
    end do
  end subroutine


  !>@brief:This subroutine updated the FEM-stored values of the coefficient phi/C in the field magma_absorp
  subroutine update_magma_coupling_coefficients(Mdims, state, saturation, ndgln, Magma_absorp,  c_phi_series, absorption_type, Magma_absorp_capped)
    implicit none
    type( state_type ), dimension( : ), intent( inout ) :: state
    real, dimension(:,:,:), intent(in) :: saturation
    type(multi_ndgln), intent(in) :: ndgln
    real, dimension(:,:,:,:), INTENT(INOUT) :: Magma_absorp
    type( multi_dimensions ), intent( in ) :: Mdims
    type(coupling_term_coef) :: coupling
    real, dimension(:), intent(in) :: c_phi_series !generated c coefficients
    logical, optional, intent(in) :: absorption_type

    real, optional, dimension(:,:,:,:), INTENT(INOUT) :: Magma_absorp_capped
    !Local variables
    integer :: mat_nod, ele, CV_ILOC, cv_inod, iphase, jphase
    real :: magma_coupling, phi
    integer:: c_phi_size ! length of c_phi_series
    real, dimension(4):: test
    real :: max_absorp_phi
    real :: phi_min=0.08 !the phi value that capped for the coupling coefficient

    type(scalar_field), pointer :: lcomponent_field
    real:: mu_l
    c_phi_size=size(c_phi_series)
    
    max_absorp_phi=1e-8

    if (present_and_true(absorption_type)) then 
      DO ELE = 1, Mdims%totele
        DO CV_ILOC = 1, Mdims%cv_nloc
            mat_nod = ndgln%mat( ( ELE - 1 ) * Mdims%mat_nloc + CV_ILOC )
            cv_inod = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
            DO IPHASE = 2, Mdims%nphase!Not phase 1
              magma_coupling = c_value(saturation(1,2, cv_inod))
              do jphase = 1, Mdims%nphase
                if (jphase == iphase) then
                  Magma_absorp(1, iphase, jphase, mat_nod ) = -magma_coupling
                else
                  Magma_absorp(1, iphase, jphase, mat_nod ) = magma_coupling
                end if
              end do
            end do
        END DO
      END DO
    else
      lcomponent_field=> extract_scalar_field( state(2), "Concentration" )
      DO ELE = 1, Mdims%totele
        DO CV_ILOC = 1, Mdims%cv_nloc
          mat_nod = ndgln%mat( ( ELE - 1 ) * Mdims%mat_nloc + CV_ILOC )
          cv_inod = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
          ! phi = max((1.0-saturation(1,1, cv_inod)),1e-5)
          phi =max(saturation(1,2, cv_inod),max_absorp_phi)
          mu_l=lcomponent_field%val(cv_inod)*1e4+(1-lcomponent_field%val(cv_inod))*1.   !liquid viscosity scaled with composition from 1 to 1e5
          do iphase =2, Mdims%nphase!Absorption is defined as a term mutiplying the velocity term, not the pressure
            Magma_absorp(1, 1, iphase, mat_nod) = phi/phi2_over_c(phi,mu_l)
            if (present(Magma_absorp_capped)) Magma_absorp_capped(1, 1, iphase, mat_nod) = max(phi,phi_min)/phi2_over_c(max(phi,phi_min),mu_l)
          end Do
        end DO
      end DO
    end if 
  contains
    !> TO INCLUDE INFORMATION ABOUT THE FUNCTION
    real function phi2_over_c(phi,mu_l)
        real, intent(in) :: phi, mu_l
        integer :: pos
        real:: portion, phi2
        !Ensure boundedness
        phi2 = min(max(0., phi), 1.)
        pos= int(phi2*(c_phi_size-1))+1
        if (pos==c_phi_size) then
          phi2_over_c=c_phi_series(c_phi_size)/mu_l
        else
          portion=(phi2-(pos-1.0)/(c_phi_size-1.0))*c_phi_size
          phi2_over_c=(c_phi_series(pos)*(1-portion)+c_phi_series(pos+1)*portion)/mu_l
        end if
      end function phi2_over_c

      real function c_value(phi)
        real, intent(in) :: phi
        integer :: pos
        real:: portion

        pos= int(phi*(c_phi_size-1))+1
        if (pos==c_phi_size) then
          c_value=phi**2/c_phi_series(c_phi_size)
        else
          portion=(phi-(pos-1.0)/(c_phi_size-1.0))*c_phi_size
          c_value=phi**2/c_phi_series(pos)*(1-portion)+phi**2/c_phi_series(pos+1)*portion
          ! c_value=c_phi_series(pos)
        end if
      end function c_value
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
      Composition_source%val(cv_inod) =0! -(Solute_new%val(1,2,cv_inod)*saturation_new%val(1,2,cv_inod) - Compostion_temp(cv_inod)*melt_temp(cv_inod))/DT  ! need to check the sign!
    end do
    ! Composition_source(1, :)=-Composition_source(2, :)


  end subroutine compute_composition_change_source

  !>@brief: Here we compute the absorption for magma
  !> Currently ONLY SINGLE DARCY PHASE (i.e., one solid phase and a second phase fluid). All phases are defined although for
  !> the upwnd the first phase contains ones so it works for the Stokes solid phase also
  !> and we can treat all the phases consistently when computing transport or the continuity equation
  subroutine Calculate_Magma_AbsorptionTerms( state, packed_state, Magma_absorp, Mdims, CV_funs, CV_GIdims, Mspars, ndgln, &
                                                    upwnd, suf_sig_diagten_bc , magma_c_phi_series,  Magma_absorp_capped)
     implicit none
     type( state_type ), dimension( : ), intent( inout ) :: state
     type( state_type ), intent( inout ) :: packed_state
     type (multi_field) :: Magma_absorp
     type( multi_dimensions ), intent( in ) :: Mdims
     type(multi_shape_funs), intent(inout) :: CV_funs
     type( multi_gi_dimensions ), intent( in )  :: CV_GIdims
     type (multi_sparsities), intent( in ) :: Mspars
     type(multi_ndgln), intent(in) :: ndgln
     type (porous_adv_coefs), intent(inout) :: upwnd
     real, dimension( :, : ), intent( inout ) :: suf_sig_diagten_bc
     real, dimension(:), intent(in) :: magma_c_phi_series

     type (multi_field), optional :: Magma_absorp_capped
     !Local variables
     type( tensor_field ), pointer :: state_viscosity
     real, dimension(:,:), allocatable :: viscosities
     integer :: i
     real, parameter :: pi = acos(0.0) * 2.0 ! Define pi

     !Take the first term only as for porous media we consider only scalar
     allocate(viscosities(Mdims%nphase, Mdims%cv_nonods))
     do i = 1,  Mdims%nphase
       state_viscosity => extract_tensor_field( state( i ), 'Viscosity' )
       call assign_val(viscosities(i, :),state_viscosity%val(1,1,:))!Take the first term only as for porous media we consider only scalar
     end do

     call Calculate_PorousMagma_adv_terms( Magma_absorp, Mdims, upwnd, viscosities, Magma_absorp_capped)

     call calculate_SUF_SIG_DIAGTEN_BC_magma( packed_state, suf_sig_diagten_bc, Mdims, CV_funs, CV_GIdims, &
         Mspars, ndgln, upwnd%adv_coef)

     deallocate(viscosities)

     contains
        !>@brief: Computes the absorption and its derivatives against the saturation
         subroutine Calculate_PorousMagma_adv_terms( Magma_absorp, Mdims, upwnd, viscosities, Magma_absorp_capped )

             implicit none
             type (multi_field), intent( inout ) :: Magma_absorp
             type( multi_dimensions ), intent( in ) :: Mdims
             type (porous_adv_coefs), intent(inout) :: upwnd
             real, dimension(:,:) :: viscosities
             type (multi_field), intent( inout ) :: Magma_absorp_capped
             !!$ Local variables:
             type(tensor_field), pointer :: satura, OldSatura
             integer :: cv_inod, iphase, ele, cv_iloc, mat_nod, icv
             real, dimension( :, : , : ), allocatable :: satura2
             real, dimension(:,:,:,:), allocatable :: Magma_absorp2
             real, dimension(:), allocatable :: max_sat
             real :: pert
             !Local parameters

             !retrieve saturation
             satura=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")

             !Set up advection coefficients
             do ele = 1, Mdims%totele
               do cv_iloc = 1, Mdims%cv_nloc
                 cv_inod = ndgln%cv(CV_ILOC + (ele-1) * Mdims%cv_nloc)
                 mat_nod = ndgln%mat(CV_ILOC + (ele-1) * Mdims%cv_nloc)
                 !Solid phase has a value of 1
                 upwnd%inv_adv_coef(1,1,1,mat_nod)=1.0; upwnd%adv_coef(1,1,1,mat_nod)=1.0
                 do iphase = 2, Mdims%nphase!if absorpt is the same we may want to reuse memory...
                   upwnd%adv_coef(1,1,iphase,mat_nod)        = Magma_absorp%val(1,1,iphase,mat_nod)
                   upwnd%capped_adv_coef(1,1,iphase,mat_nod) = Magma_absorp_capped%val(1,1,iphase,mat_nod)
                   !Now the inverse
                   upwnd%inv_adv_coef(1,1,iphase,mat_nod) = 1./upwnd%adv_coef(1,1,iphase,mat_nod)
                 end do
               end do
             end do

             !Introduce perturbation, positive for the increasing and negative for decreasing phase
             !Make sure that the perturbation is between bounds
             PERT = 0.00000; allocate(Max_sat(Mdims%nphase), SATURA2(1, Mdims%nphase, Mdims%cv_nonods))
             OldSatura=>extract_tensor_field(packed_state,"PackedOldPhaseVolumeFraction")
             do ele = 1, Mdims%totele
                 do cv_iloc = 1, Mdims%cv_nloc
                     cv_inod = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                     Max_sat = 1. !- sum(Immobile_fraction(:, ele)) + Immobile_fraction(:, ele)
                     DO IPHASE = 2, Mdims%nphase
                       SATURA2(1, IPHASE, cv_inod) = satura%val(1, IPHASE, cv_inod) + sign(PERT, satura%val(1, iphase, cv_inod)-Oldsatura%val(1,iphase, cv_inod))
                       !If out of bounds then we perturbate in the opposite direction
                       if (satura2(1, IPHASE, cv_inod) > Max_sat(iphase)) then
                           SATURA2(1, IPHASE, cv_inod) = SATURA2(1, IPHASE, cv_inod) - 2. * sign(PERT, satura%val(1,iphase, cv_inod)-Oldsatura%val(1,iphase, cv_inod))
                       end if
                     end do
                 end do
             end do

             !Compute absorption gven a perturbed saturation
             allocate(Magma_absorp2(1, 1, Mdims%nphase, Mdims%mat_nonods))
             call update_magma_coupling_coefficients(Mdims, state, SATURA2, ndgln, Magma_absorp2,  magma_c_phi_series)
             DO ELE = 1, Mdims%totele
               DO CV_ILOC = 1, Mdims%cv_nloc
                 mat_nod = ndgln%mat( ( ELE - 1 ) * Mdims%mat_nloc + CV_ILOC )
                 cv_inod = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
                 !Absorption for phase 1 is constant so the gradient is zero
                 upwnd%adv_coef_grad(1, 1, 1, mat_nod)=0.0;
                 DO IPHASE = 2, Mdims%nphase
                   ! This is the gradient
                   upwnd%adv_coef_grad(1, 1, iphase, mat_nod) = (Magma_absorp2( 1,1, iphase ,mat_nod) -&
                   Magma_absorp%val( 1,1, iphase ,mat_nod)) / ( SATURA2(1,iphase, cv_inod ) - satura%val(1,iphase, cv_inod))
                 END DO
               end do
             END DO
             deallocate( satura2, Max_sat, Magma_absorp2)

         end subroutine Calculate_PorousMagma_adv_terms

         !>@brief: Computes the absorption and its derivatives against the saturation on the boundary
         subroutine calculate_SUF_SIG_DIAGTEN_BC_magma( packed_state, suf_sig_diagten_bc, Mdims, CV_funs, CV_GIdims, &
             Mspars, ndgln, adv_coef)
             implicit none
             type( state_type ), intent( inout ) :: packed_state
             type(multi_dimensions), intent(in) :: Mdims
             type(multi_GI_dimensions), intent(in) :: CV_GIdims
             type(multi_shape_funs), intent(inout) :: CV_funs
             type (multi_sparsities), intent(in) :: Mspars
             type(multi_ndgln), intent(in) :: ndgln
             real, dimension(:,:,:,:) :: adv_coef
             real, dimension( Mdims%stotel * Mdims%cv_snloc * Mdims%nphase, Mdims%ndim ), intent( inout ) :: suf_sig_diagten_bc
             ! local variables
             integer :: iphase, ele, sele, cv_siloc, cv_snodi, cv_snodi_ipha, iface,  &
                 ele2, sele2, cv_iloc, idim, jdim, i, mat_nod, cv_nodi
             real :: sigma_out!, mat, mat_inv
             ! real, dimension( Mdims%ndim, Mdims%ndim ) :: mat_ones, mat, mat_inv
             integer, dimension( CV_GIdims%nface, Mdims%totele) :: face_ele
             integer, dimension( Mdims%cv_snloc ) :: cv_sloc2loc
             integer, dimension( :, :, : ),  allocatable :: wic_u_bc, wic_vol_bc
             integer, parameter :: WIC_BC_DIRICHLET = 1
             type(tensor_field), pointer :: velocity, satura, perm
             type(tensor_field) :: velocity_BCs, satura_BCs
             !Local parameters


             !Get from packed_state
             satura=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
             velocity=>extract_tensor_field(packed_state,"PackedVelocity")

             allocate(wic_u_bc(velocity%dim(1),velocity%dim(2),&
                 surface_element_count(velocity)))
             allocate(wic_vol_bc(satura%dim(1),satura%dim(2),&
                 surface_element_count(satura)))
             call get_entire_boundary_condition(velocity,&
                 ['weakdirichlet'],velocity_BCs,WIC_U_BC)
             call get_entire_boundary_condition(satura,&
                 ['weakdirichlet'],satura_BCs,WIC_vol_BC)

             !We initialise as 1 as the first phase must contain a 1
             suf_sig_diagten_bc = 1.
             face_ele = 0
             call calc_face_ele( face_ele, Mdims%totele, Mdims%stotel, CV_GIdims%nface, &
                 Mspars%ELE%fin, Mspars%ELE%col, Mdims%cv_nloc, Mdims%cv_snloc, Mdims%cv_nonods, ndgln%cv, ndgln%suf_cv, &
                 CV_funs%cv_sloclist, Mdims%x_nloc, ndgln%x )
                 do ele = 1, Mdims%totele
                     !Get properties from packed state
                     do iface = 1, CV_GIdims%nface
                         ele2  = face_ele( iface, ele )
                         sele2 = max( 0, -ele2 )
                         sele  = sele2
                         if ( sele > 0 ) then
                             do iphase = 2, Mdims%nphase
                                 if ( wic_u_bc(1,iphase,sele) /= WIC_BC_DIRICHLET .and. &
                                     wic_vol_bc(1,iphase,sele) == WIC_BC_DIRICHLET ) then
                                     cv_sloc2loc( : ) = CV_funs%cv_sloclist( iface, : )
                                     do cv_siloc = 1, Mdims%cv_snloc
                                         cv_iloc = cv_sloc2loc( cv_siloc )
                                         cv_snodi = ( sele - 1 ) * Mdims%cv_snloc + cv_siloc
                                         cv_nodi = ndgln%suf_cv(cv_snodi)
                                         cv_snodi_ipha = cv_snodi + ( iphase - 1 ) * Mdims%stotel * Mdims%cv_snloc
                                         mat_nod = ndgln%mat( (ele-1)*Mdims%cv_nloc + cv_iloc  )
                                         !For the time being use the interior absorption, this will need to be changed when multiphase like the porous media one
                                        suf_sig_diagten_bc( cv_snodi_ipha, 1 : Mdims%ndim ) = adv_coef(1,1,iphase, mat_nod)
                                     end do
                                 end if
                             end do
                         end if
                     end do
                 end do
              call deallocate(velocity_BCs)
             call deallocate(satura_BCs)
             deallocate(wic_u_bc, wic_vol_bc)
             return
         end subroutine calculate_SUF_SIG_DIAGTEN_BC_magma

  end subroutine Calculate_Magma_AbsorptionTerms

end module multi_magma
