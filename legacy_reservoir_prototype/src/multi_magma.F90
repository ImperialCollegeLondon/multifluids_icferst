
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
    !use Copy_Outof_State
    !use petsc_tools
    ! Hs- Solidus (enthalpy)
    ! Hl - Liquidus (enthalpy)
    ! Ts- Solidus (temperature)
    ! Tl - Liquidus (temperature)
    ! Lf - latent heat
    !Parameters waiting to be decided if they will become fields or inputs from diamond as scalars
    !I would rather have some more descriptive names for these parameters
    real, parameter :: Hs = 0.,Hl = 1e5,Ts = 1123.15,Tl = 1e3, Lf = 550000.0
    ! Ae - eutetic
    ! A1,B1,C1,A2,B2,C2 - phase behaviour parameters
    real, parameter :: A1= 50.,B1= -360,C1= 1433.15,A2= 0.,B2= 0.,C2 = 0., Ae = 1.0
contains

  subroutine initilize_enthalpy_from_temperature(state, packed_state, Mdims)
    type( state_type ), dimension( : ), intent( inout ) :: state
    type( state_type ), intent( inout ) :: packed_state
    type( tensor_field ), pointer :: enthalpy, temperature, Density_Cp,  den
    type( tensor_field) :: enthalpy_BCs, temperature_BCs
    type( vector_field ), pointer :: porosity
    type( multi_dimensions), intent( in ) :: Mdims
    character (50) :: field_path
    integer :: cv_nodi, i

    den => extract_tensor_field( packed_state,"PackedDensity" )
    Density_Cp =>  extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )
    enthalpy => extract_tensor_field( packed_state,"PackedEnthalpy" )
    temperature =>  extract_tensor_field( packed_state, "PackedTemperature" )
    porosity=>extract_vector_field(packed_state,"MeanPoreCV")

    do cv_nodi = 1, Mdims%cv_nonods
      enthalpy%val(1,1,cv_nodi)=temperature%val(1,1,cv_nodi)*Density_cp%val(1,1,cv_nodi)+porosity%val(1,cv_nodi)*den%val(1,1,cv_nodi)*Lf
    end do


    ! call get_entire_boundary_condition(temperature,['weakdirichlet'],tmeperature_BCs,WIC_U_BC_ALL) !TOC update boundary conditions
    !
    ! !TOC should update it w.r.t porosity
    ! enthalpy_BCs=tmeperature_BCs*Density_cp(1,1,cv_nodi)+porosity%val(1,1,cv_nodi)*den(1,cv_nodi)*Lf
    !
    ! position => extract_vector_field(states(1), "Coordinate")
    ! field_path=enthalpy%option_path
    !
    ! call set_vector_boundary_conditions_values(states(1), enthalpy, &
    !      trim(field_path)//'/prognostic/boundary_conditions', &
    !      position)

  end subroutine


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
      real function get_Liquidus(Bulk_comp)
      implicit none
      real, intent(in) :: Bulk_comp

      if (Bulk_comp > Ae) then
        get_Liquidus = A2*(Bulk_comp**2) + B2*Bulk_comp + C2
      else
        get_Liquidus = A1*(Bulk_comp**2) + B1*Bulk_comp + C1
      end if

    end function

    ! Function for the Solidus (temp)
    real function get_Solidus(Bulk_comp)
      implicit none
      real, intent(in) :: Bulk_comp
      !Currently we consider a flat line
      get_Solidus = Ts
    end function


    ! Function for the liquidus (Enthalpy)
    real function get_Enthalpy_Liquidus(Bulk_comp, Cp, rho)
      implicit none
      real, intent(in) :: Bulk_comp, Cp, rho
      get_Enthalpy_Liquidus = get_Liquidus(Bulk_comp) * Cp  + Lf
    end function

    ! Function for the Solidus (Enthalpy)
    real function get_Enthalpy_Solidus(Bulk_comp, Cp, rho)
    implicit none
    real, intent(in) :: Bulk_comp, Cp, rho
    !Currently we consider a flat line
      get_Enthalpy_Solidus = get_Solidus(Bulk_comp) * Cp!* rho

    end function

    ! Function to find solid composition
    ! real function get_Solid_Composition(Bulk_comp, Temperature)
    ! !TODO: We probably want to create a subroutine using WHERE statements and fields to do all of this at once
    ! implicit none
    ! real, intent(in) :: Bulk_comp, Temperature
    !
    !   if (Temperature > get_Liquidus(Bulk_comp)) then
    !   	if (Bulk_comp > Ae) then
    !   		get_Solid_Composition = 1.0
    !   	else
    !   		get_Solid_Composition = 0.0
    !   	end if
    !   else if (Temperature < get_Solidus(Bulk_comp)) then
    !   	get_Solid_Composition = Bulk_comp
    !   else
    !   	if (Bulk_comp > Ae) then
    !   		get_Solid_Composition = 1.0
    !   	else
    !   		get_Solid_Composition = 0.0
    !   	end if
    !   end if
    !
    ! end function

    ! subroutine for updating the solid composition and fluid composition
    subroutine cal_solidfluidcomposition(state, packed_state, Mdims)
      type( state_type ), dimension( : ), intent( inout ) :: state
      type( state_type ), intent( inout ) :: packed_state
      type( multi_dimensions), intent( in ) :: Mdims
      type( tensor_field), pointer :: Composition, temperature
      type( vector_field), pointer :: porosity
      type( scalar_field), pointer :: BC !bulk composition

      integer :: cv_nodi

      BC=> extract_scalar_field(state(1), "BulkComposition")
      Composition=>extract_tensor_field(packed_state,"PackedComposition")
      temperature =>  extract_tensor_field( packed_state, "PackedTemperature" )

      do cv_nodi = 1, Mdims%cv_nonods
          if (temperature%val(1,1,cv_nodi)>get_Liquidus(BC%val(cv_nodi))) then
            if(BC%val(cv_nodi)>Ae) then
              Composition%val(1,1,cv_nodi)=1.0 !solid
            else
              Composition%val(1,1,cv_nodi)=0.0
            end if

            Composition%Val(1,2,cv_nodi)=BC%val(cv_nodi)!fluid

          else if(temperature%val(1,1,cv_nodi)<get_Solidus(BC%val(cv_nodi))) then
            Composition%val(1,1,cv_nodi)=BC%val(cv_nodi)
            Composition%val(1,2,cv_nodi)=Ae
          else
            if(BC%val(cv_nodi)>Ae) then
              Composition%val(1,1,cv_nodi)=1.0
              ! now that A2=B2=0, this expression is meaningless for fluid compositino and will get a error for devide A2 (=0)
              Composition%val(1,2,cv_nodi)= 1.0
              !FluidComposition%val(1,1,cv_nodi)=(-B2 + sqrt(( B1**2 ) - 4. * A2 * ( C2 - temperature%val(1,1,cv_nodi)))) / ( 2. * A2)
            else
              Composition%val(1,1,cv_nodi)=0.0
              Composition%val(1,2,cv_nodi)=(-B1 - sqrt(( B1**2 ) - 4. * A1 * ( C1 - temperature%val(1,1,cv_nodi)))) / ( 2. * A1)
            end if
          end if
      end do

    end subroutine

    ! Function to find liquid composition
    ! real function get_Liquid_Composition(Bulk_comp, Temperature)
    ! !TODO: We probably want to create a subroutine using WHERE statements and fields to do all of this at once
    ! implicit none
    ! real, intent(in) :: Bulk_comp, Temperature
    !
    ! 	if (Temperature > get_Liquidus(Bulk_comp)) then
    ! 		get_Liquid_Composition = Bulk_comp
    ! 	else if (Temperature < get_Solidus(Bulk_comp)) then
    ! 		get_Liquid_Composition = Ae
    ! 	else
    ! 		if (Bulk_comp < Ae) then
    ! 			get_Liquid_Composition = (-B1 - sqrt(( B1**2 ) - 4. * A1 * ( C1 - Temperature ))) / ( 2. * A1)
    ! 		else
    ! 			get_Liquid_Composition = (-B2 + sqrt(( B1**2 ) - 4. * A2 * ( C2 - Temperature ))) / ( 2. * A2)
    ! 		end if
    ! 	end if
    !
    ! end function

    !Subroutine to obtain the temperature field for a given enthalpy field
    subroutine enthalpy_to_temperature(Mdims, packed_state)
    implicit none
    !Global variables
    type( state_type ), intent( inout ) :: packed_state
    type(multi_dimensions), intent( in ) :: Mdims

    !Local variables
    integer :: cv_nodi, iphase                          !Temporary until deciding if creating a Cp in packed_state as well
    type( tensor_field ), pointer :: enthalpy, temperature,  saturation, FluidComposition, Density_Cp , den
    real, dimension(Mdims%cv_nonods) :: enthalpy_dim
    !real, parameter :: tol = 1e-5
      !Temporary until deciding if creating a Cp in packed_state as well
      den => extract_tensor_field( packed_state,"PackedDensity" )
      Density_Cp =>  extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )
      enthalpy => extract_tensor_field( packed_state,"PackedEnthalpy" )
      temperature =>  extract_tensor_field( packed_state, "PackedTemperature" )
      saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
      !saturation =>  extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )!First phase is rock presence

  !IT CONSIDERS ONE SINGLE TEMPERATURE!!! WHICH MAKE SENSE BECAUSE AT THAT TIME SCALE THINGS ARE IN EQUILIBRIUM...
      ! do iphase =1, nphase
      iphase = 1
      !Calculate temperature using the generic formula T = (H - Lf*porosity)/Cp
      temperature%val(1,1,:) = (enthalpy%val(1,1,:) - Lf * saturation%val(1,2,:))/Density_Cp%val(1,iphase,:)

      !Now add corrections if required
      where (saturation%val(1,2, :)>0.) !If porosity > 0.
        where (temperature%val(1,1,:) < Ts ) !If there is mixture, temperature cannot be below the solidus
          temperature%val(1,1,:) = Ts!It seems that the Solidus line is flat
        end where
      ! elsewhere  !not necesarry
      !   temperature%val(1,1,:) = enthalpy_dim* den%val(1,iphase,:)/Density_Cp%val(1,iphase,:)
      end where

      !Convert to dimensionles temperature
      !call magma_field_dim_to_nondim(temperature%val(1,iphase,:), temperature%val(1,iphase,:), Ts, Tl)

      !Now for consistency populate the temperature of all the phase
      ! do iphase = 2, Mdims%ndim
      !   temperature%val(1,iphase,:) = temperature%val(1,1,:)
      ! end do

    end subroutine enthalpy_to_temperature


    !Compute porosity given an enthalpy field and a BulkComposition
    subroutine porossolve(state,packed_state, Mdims, ndgln)
    implicit none
    !Global variables
    type(multi_dimensions), intent( in ) :: Mdims
    type( state_type ), intent( inout ) :: packed_state
    type( state_type ), dimension(:), intent( inout ) :: state
    type(multi_ndgln), intent(in) :: ndgln
    !Local variables
    integer ::  iphase, k, i
                                                      !Temporary until deciding if creating a Cp in packed_state as well
    type( tensor_field ), pointer :: enthalpy, den, Density_Cp, saturation
    type (tensor_field), pointer :: FluidComposition, MatrixComposition
    type (vector_field), pointer :: porosity, FEMporosity
    type (scalar_field), pointer :: BC !BUlk compositon
    !real, dimension(Mdims%cv_nonods) :: enthalpy_dim
    real :: test_poro_prev, test_poro!Temporary porosity
    real :: fx, fdashx, Loc_Cp, rho
    !Parameters for the non_linear solvers (Maybe a newton solver here makes sense?)
    real, parameter :: tol = 1e-2
    integer, parameter :: max_its = 25
    !!
    real :: ELE,CV_ILOC, cv_nodi
      !Temporary until deciding if creating a Cp in packed_state as well
      den => extract_tensor_field( packed_state,"PackedDensity" )
      Density_Cp =>  extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )
      enthalpy => extract_tensor_field( packed_state,"PackedEnthalpy" )
      !saturation =>  extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )!First phase is rock presence
      porosity=>extract_vector_field(packed_state,"MeanPoreCV") !use mean CV porosity calculation
      FEMporosity=>extract_vector_field(packed_state,"Porosity") !update the porosity in the end
      BC=> extract_scalar_field(state(1), "BulkComposition")
      iphase = 1
      ! cv_nodi = 1
      ! print *, enthalpy%val(1,1,cv_nodi)
      ! print *, get_Enthalpy_Liquidus(BC%val(cv_nodi),Loc_Cp,
      do cv_nodi = 1, Mdims%cv_nonods
        Loc_Cp = Density_Cp%val(1,iphase,cv_nodi)/den%val(1,iphase,cv_nodi)!Temporary until we decide if we create a Cp field
        rho=den%val(1,iphase,cv_nodi)
        IF (enthalpy%val(1,1,cv_nodi)>get_Enthalpy_Liquidus(BC%val(cv_nodi),Loc_Cp, rho)) then
          porosity%val(1, cv_nodi)=1.
        ELSE IF (enthalpy%val(1,1,cv_nodi)<get_Enthalpy_Solidus(BC%val(cv_nodi),Loc_Cp, rho)) then
          porosity%val(1, cv_nodi)=0.
        ELSE
          if (BC%val(cv_nodi) < Ae) then
            k=1
            test_poro = 1.0
            do while (k<Max_its .and. abs(test_poro - test_poro_prev) < tol)
              fx = (Lf/Loc_Cp)*test_poro**3 + (C1-enthalpy%val(1,iphase,cv_nodi)/Loc_Cp)*test_poro**2 + B1*BC%val(cv_nodi)*&
              test_poro + A1*BC%val(cv_nodi)
              test_poro_prev = test_poro
              test_poro = test_poro - fx/fdashx
              k=k+1
            end do
          else if (BC%val(cv_nodi) >= Ae) then
            k=1
            test_poro = 1.0
            do while (k<Max_its .and. abs(test_poro - test_poro_prev) < tol .and. test_poro<=1)
              fx= (Lf/Loc_Cp)**test_poro**3+(C2-enthalpy%val(1,iphase,cv_nodi)/Loc_Cp+A2+B2)*test_poro**2+(2*A2*BC%val(cv_nodi)-2*A2+&
              B2*BC%val(cv_nodi)-B2)*test_poro+A2*(BC%val(cv_nodi)-1)**2

              fdashx = 3.*(Lf/Loc_Cp)*test_poro**2 + 2.*(C2-enthalpy%val(1,iphase, cv_nodi)/Loc_Cp + A2 + B2)*test_poro +&
              (2.*A2*BC%val(cv_nodi) - 2.*A2 + B2*BC%val(cv_nodi) - B2)

              test_poro_prev = test_poro!Store to check convergence the previous value
              test_poro = test_poro - fx/fdashx
              k=k+1
            end do
          end if
          porosity%val(1, cv_nodi)=test_poro
        END IF
      end do

      !HH Project CV_porosity to FEM_porosity, TOC
      DO ELE = 1, Mdims%totele
          FEMporosity%val(1,ELE)=0.
        DO CV_ILOC = 1, Mdims%cv_nloc
          CV_NODI = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
          FEMporosity%val(1,ELE)=FEMporosity%val(1,ELE)+1./Mdims%cv_nloc*Porosity%val(1,CV_NODI)
        END DO
      END DO

    end subroutine


end module multi_magma
