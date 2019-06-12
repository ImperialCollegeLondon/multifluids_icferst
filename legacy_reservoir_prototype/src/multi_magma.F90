
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
    real, parameter :: Hs = 0.,Hl = 1e5,Ts = 1123.15,Tl = 1e3, Lf = 550000.0
    ! Ae - eutetic
    ! A1,B1,C1,A2,B2,C2 - phase behaviour parameters
    real, parameter :: A1= 50.,B1= -360,C1= 1433.15,A2= 0.,B2= 0.,C2 = 0., Ae = 1.0
contains
    !========================================================
    !Subroutine to convert between Dimensional and Non-Dimensional fields (temperature or enthalpy)
    !========================================================
    subroutine magma_field_nondim_to_dim(field, field_dim, field_solidus,field_liquidus)
    implicit none
    real, dimension(:), intent(in):: field
    real, dimension(:), intent(inout):: field_dim
    real, intent(in):: field_solidus,field_liquidus
    !Local variables
        field_dim = field*(field_liquidus-field_solidus) + field_solidus

    end subroutine

    subroutine magma_field_dim_to_nondim(field, field_nondim, field_solidus,field_liquidus)
    implicit none
    real, dimension(:), intent(in):: field
    real, dimension(:), intent(inout):: field_nondim
    real, intent(in):: field_solidus,field_liquidus

        field_nondim = (field - field_solidus)/(field_liquidus-field_solidus)

    end subroutine

    ! Function for the liquidus (temp)
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
    real function get_Enthalpy_Liquidus(Bulk_comp, Cp)
    implicit none
    real, intent(in) :: Bulk_comp, Cp

      get_Enthalpy_Liquidus = get_Liquidus(Bulk_comp) * Cp + Lf

    end function

    ! Function for the Solidus (Enthalpy)
    real function get_Enthalpy_Solidus(Bulk_comp, Cp)
    implicit none
    real, intent(in) :: Bulk_comp, Cp

    !Currently we consider a flat line
      get_Enthalpy_Solidus = get_Solidus(Bulk_comp) * Cp

    end function

    ! Function to find solid composition
    real function get_Solid_Composition(Bulk_comp, Temperature)
    !TODO: We probably want to create a subroutine using WHERE statements and fields to do all of this at once
    implicit none
    real, intent(in) :: Bulk_comp, Temperature

      if (Temperature > get_Liquidus(Bulk_comp)) then
      	if (Bulk_comp > Ae) then
      		get_Solid_Composition = 1.0
      	else
      		get_Solid_Composition = 0.0
      	end if
      else if (Temperature < get_Solidus(Bulk_comp)) then
      	get_Solid_Composition = Bulk_comp
      else
      	if (Bulk_comp > Ae) then
      		get_Solid_Composition = 1.0
      	else
      		get_Solid_Composition = 0.0
      	end if
      end if

    end function

    ! Function to find liquid composition
    real function get_Liquid_Composition(Bulk_comp, Temperature)
    !TODO: We probably want to create a subroutine using WHERE statements and fields to do all of this at once
    implicit none
    real, intent(in) :: Bulk_comp, Temperature

    	if (Temperature > get_Liquidus(Bulk_comp)) then
    		get_Liquid_Composition = Bulk_comp
    	else if (Temperature < get_Solidus(Bulk_comp)) then
    		get_Liquid_Composition = Ae
    	else
    		if (Bulk_comp < Ae) then
    			get_Liquid_Composition = (-B1 - sqrt(( B1**2 ) - 4. * A1 * ( C1 - Temperature ))) / ( 2. * A1)
    		else
    			get_Liquid_Composition = (-B2 + sqrt(( B1**2 ) - 4. * A2 * ( C2 - Temperature ))) / ( 2. * A2)
    		end if
    	end if

    end function

    !Subroutine to obtain the temperature field for a given enthalpy field
    subroutine enthalpy_to_temperature(Mdims, packed_state)
    implicit none
    !Global variables
    type( state_type ), intent( inout ) :: packed_state
    type(multi_dimensions), intent( in ) :: Mdims

    !Local variables
    integer :: cv_nodi, iphase                          !Temporary until deciding if creating a Cp in packed_state as well
    type( tensor_field ), pointer :: enthalpy, temperature, den, Density_Cp, saturation
    real, dimension(Mdims%cv_nonods) :: enthalpy_dim
    real, parameter :: tol = 1e-5
      !Temporary until deciding if creating a Cp in packed_state as well
      den => extract_tensor_field( packed_state,"PackedDensity" )
      Density_Cp =>  extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )
      enthalpy => extract_tensor_field( packed_state,"PackedEnthalpy" )
      temperature =>  extract_tensor_field( packed_state, "PackedTemperature" )
      saturation =>  extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )!First phase is rock presence

  !IT CONSIDERS ONE SINGLE TEMPERATURE!!! WHICH MAKE SENSE BECAUSE AT THAT TIME SCALE THINGS ARE IN EQUILIBRIUM...
      ! do iphase =1, nphase
      iphase = 1
      call magma_field_nondim_to_dim(enthalpy%val(1,iphase,:), enthalpy_dim, Hs, Hl)

      !Calculate temperature using the generic formula T = (H - Lf*porosity)/Cp     !1 - saturation == to porosity
      temperature%val(1,1,:) = (enthalpy_dim - Lf * (1. - saturation%val(1,1,:)) * den%val(1,iphase,:))/Density_Cp%val(1,iphase,:)


      !Now add corrections if required
      where (1. > tol + saturation%val(1,1,:)) !If porosity > 0.
        where (temperature%val(1,1,:) < Ts ) !If there is mixture, temperature cannot be below the solidus
          temperature%val(1,1,:) = Ts!It seems that the Solidus line is flat
        end where
      elsewhere
        temperature%val(1,1,:) = enthalpy_dim* den%val(1,iphase,:)/Density_Cp%val(1,iphase,:)
      end where

      !Convert to dimensionles temperature
      call magma_field_dim_to_nondim(temperature%val(1,iphase,:), temperature%val(1,iphase,:), Ts, Tl)

      !Now for consistency populate the temperature of all the phase
      do iphase = 2, Mdims%ndim
        temperature%val(1,iphase,:) = temperature%val(1,1,:)
      end do

    end subroutine enthalpy_to_temperature


    ! Subroutine for porosity.
    !
    ! z1,z2 - boundary nodes
    ! ha - enthalpy at t+1
    ! pa - porosity at t+1
    ! ca - bulk composition at t+1
    ! Hs - Solidus enthalpy
    ! Hl - Liquidus enthalpy
    ! Lf - Latent heat
    ! Cp - specific heat capacity
    ! Tl - Liquidus temperature
    ! Ts - Solidus Temperature
    ! Ae - Eutetic
    ! A1, B1, C1, A2, B2, C2 - phase behaviour parameters
    !



    !========================================================
    subroutine porossolve(Mdims, packed_state, state)
    implicit none
    !Global variables
    type(multi_dimensions), intent( in ) :: Mdims
    type( state_type ), intent( inout ) :: packed_state
    type( state_type ), dimension(:), intent( inout ) :: state

    !Local variables
    integer :: cv_nodi, iphase, k
                                                      !Temporary until deciding if creating a Cp in packed_state as well
    type( tensor_field ), pointer :: enthalpy, temperature, den, Density_Cp, saturation
    type (scalar_field), pointer :: BulkComposition
    real, dimension(Mdims%cv_nonods) :: enthalpy_dim
    real, parameter :: tol = 1e-5
    real :: test_poro!Temporary porosity
    real :: fx, fdashx, Loc_Cp

      !Temporary until deciding if creating a Cp in packed_state as well
      den => extract_tensor_field( packed_state,"PackedDensity" )
      Density_Cp =>  extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )
      enthalpy => extract_tensor_field( packed_state,"PackedEnthalpy" )
      temperature =>  extract_tensor_field( packed_state, "PackedTemperature" )
      saturation =>  extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )!First phase is rock presence
      BulkComposition => extract_scalar_field(state(1), "BulkComposition")! I

      iphase = 1
      call magma_field_nondim_to_dim(enthalpy%val(1,iphase,:), enthalpy_dim, Hs, Hl)
      do cv_nodi = 1, Mdims%cv_nonods

        Loc_Cp = Density_Cp%val(1,iphase,cv_nodi)/den%val(1,iphase,cv_nodi)!Temporary until we decide if we create a Cp field

      	! H = enth_nondim_to_dim(ha(i),Hs,Hl)
        !Original method is on nondimensional space but it should be fine in dimensional space?
      	! if(enthalpy%val(1,iphase,cv_nodi) > enth_dim_to_nondim(HLiquidus(ca(i),Cp,Ae,Lf,A1,B1,C1,A2,B2,C2),Hs,Hl)) then
      	! 	test_poro = 1
        !   saturation%val(1,iphase,cv_nodi) = 1 - test_poro
        !   cycle
      	! end if
        !
      	! if(enthalpy%val(1,iphase,cv_nodi) > enth_dim_to_nondim(HSolidus(ca(i),Ts,Cp),Hs,Hl)) then
      	! 	test_poro = 0
        !   saturation%val(1,iphase,cv_nodi) = 1 - test_poro
        !   cycle
      	! end if

        !In dimensional space we save some computations, it should not matter... ask MATT or double check this
        if(enthalpy_dim(cv_nodi) > get_Enthalpy_Liquidus( BulkComposition%val(cv_nodi), Loc_Cp) ) then
      		test_poro = 1
          saturation%val(1,iphase,cv_nodi) = 1 - test_poro
          cycle
      	end if

      	if(enthalpy_dim(cv_nodi) > get_Enthalpy_Solidus(BulkComposition%val(cv_nodi), Loc_Cp)) then
      		test_poro = 0
          saturation%val(1,iphase,cv_nodi) = 1 - test_poro
          cycle
      	end if


        test_poro = 1.0
        ! There is an iteration here to update porosity pa based on composition.  Does not use velocity.
      	if (BulkComposition%val(cv_nodi) < Ae) then

          do k = 1,25!Why 25?, maybe better until test_poro does not vary more than 1%?
      			fx = (Lf/Loc_Cp)*test_poro**3 + (C1-enthalpy_dim(cv_nodi)/Loc_Cp)*test_poro**2 + (B1*BulkComposition%val(cv_nodi))*&
                                    test_poro + A1*BulkComposition%val(cv_nodi)**2

            fdashx = 3*(Lf/Loc_Cp)*test_poro**2 + 2*(C1-enthalpy_dim(cv_nodi)/Loc_Cp)*test_poro + (B1*BulkComposition%val(cv_nodi))

            test_poro = test_poro - fx/fdashx
      		end do

      	else

      		do k = 1,25!Why 25?, maybe better until test_poro does not vary more than 1%?

        		fx = (Lf/Loc_Cp)*test_poro**3 + (C2-enthalpy_dim(cv_nodi)/Loc_Cp + A2 + B2)*test_poro**2 + &
                  (2.*A2*BulkComposition%val(cv_nodi) - 2.*A2 + B2*BulkComposition%val(cv_nodi) - B2)*test_poro + &
                  A2*(BulkComposition%val(cv_nodi)**2-2.*BulkComposition%val(cv_nodi) + 1.)

        		fdashx = 3.*(Lf/Loc_Cp)*test_poro**2 + 2.*(C2-enthalpy_dim(cv_nodi)/Loc_Cp + A2 + B2)*test_poro +&
                  (2.*A2*BulkComposition%val(cv_nodi) - 2.*A2 + B2*BulkComposition%val(cv_nodi) - B2)

      			test_poro = test_poro - fx/fdashx

      			if(test_poro >= 1.0) then
      				test_poro = 1.0
              saturation%val(1,iphase,cv_nodi) = 1. - test_poro
              exit!Exit current loop and move to next cv_nodi
      			end if

      		end do

      	end if
        !Assign testing value to saturation of phase 1
        saturation%val(1,iphase,cv_nodi) = 1. - test_poro
      end do

      !Assign now the liquid saturation
      iphase = 2 !I think this current formula only accepts two phases
      saturation%val(1,iphase,:) = 1. - saturation%val(1,1,:)
    end subroutine


end module multi_magma
