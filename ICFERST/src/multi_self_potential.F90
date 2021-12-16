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

module multi_SP

    use fldebug
    use futils
    use spud
    use fields
    use global_parameters, only: OPTION_PATH_LEN, PYTHON_FUNC_LEN, PI, is_porous_media
    use vector_tools
    use state_module
    use fields
    use multi_data_types
    use multi_tools
    use Copy_Outof_State
    use multiphase_1D_engine
    use boundary_conditions_from_options
    private

    public :: Assemble_and_solve_SP

    contains


      !>@brief: This subroutine computes the saturated rock conductivity based on Tiab and Donaldson, 2004 formula.
      !> and also includes the temperature/concentration dependency detailed in Sen and Goode, 1992
      !> IMPORTANT: Water needs to be phase 1! TODO Maybe add a REMINDER when running with SP solver and to run with Kelvin, maybe this for all the temperature cases
      !> TODO include the option to project to FE using PROJ_CV_TO_FEM
      subroutine Assemble_and_solve_SP( Mdims, state, packed_state, ndgln, Mmat, Mspars, CV_funs, CV_GIdims)
        implicit none

        type(multi_dimensions), intent( in ) :: Mdims
        type( state_type ), dimension(:), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(multi_ndgln), intent(in) :: ndgln
        type(multi_shape_funs), intent(inout) :: CV_funs
        type(multi_sparsities), intent(in) :: Mspars
        type (multi_matrices), intent(inout) :: Mmat
        type(multi_GI_dimensions), intent(in) :: CV_GIdims
        !Local variables
        integer :: k, cv_inod, nfields, stat
        type( scalar_field ), pointer :: Solution
        type( tensor_field ), pointer :: Temperature, Concentration, Saturation, density, pressure
        real, dimension(:,:,:), allocatable :: F_fields, K_fields
        real, dimension(:,:), allocatable ::rock_sat_conductivity
        character(len=option_path_len) :: solver_option_path = "/solver_options/Linear_solver"
        type(scalar_field), pointer :: SelfPotential
        type( vector_field ), pointer :: X_ALL
        integer :: reference_nod
        logical :: reference_node_owned
        real :: reference_value, top_coordinate, gravity_magnitude

      !     !Retrieve fields
        X_ALL => extract_vector_field( packed_state, "PressureCoordinate" )
        Saturation=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
        Pressure=>extract_tensor_field(packed_state,"PackedFEPressure")
        density =>extract_tensor_field(packed_state,"PackedDensity")
        nfields = 1
        Concentration=>extract_tensor_field(packed_state,"PackedConcentration", stat)
        Temperature=>extract_tensor_field(packed_state,"PackedTemperature", stat)
        if (has_concentration) nfields = nfields + 1
        if (has_temperature) nfields = nfields + 1
        allocate(F_fields(nfields, 1, Mdims%cv_nonods), K_fields(nfields, 1, Mdims%cv_nonods))
        allocate(rock_sat_conductivity(1, Mdims%cv_nonods))
        !#############Fill up F_fields#######################
        !First obtain the highest point of the model (either y in 2D or Z in 3D are the vertical coordinates)
        top_coordinate = maxval(X_ALL%val(Mdims%ndim,:)); call allmax(top_coordinate)
        call get_option( "/physical_parameters/gravity/magnitude", gravity_magnitude, default = 0. )
        do cv_inod = 1, Mdims%cv_nonods
          !Water potential, i.e. = P - rho * g * h
          F_fields(1, 1, cv_inod) = Pressure%val(1,1,cv_inod) - gravity_magnitude * density%val(1,1,cv_inod) * (top_coordinate - X_ALL%val(Mdims%ndim,cv_inod))
          k = 1
          ! F_fields(1, 1, cv_inod) = Saturation%val(1, 1, cv_inod); k = 1
          if (has_concentration) then
            !Convert to Moles/Litre
            F_fields(2, 1, cv_inod) = Concentration%val(1, 1, cv_inod) * 1000 !To convert from mol/m^3 to mol/l which is how the formulae are defined
            k = 2
          end if
          !Here the temperature can be in Kelvin or Celsius as we are looking at gradients
          if (has_temperature) F_fields(k+1, 1, cv_inod) =  Temperature%val(1, 1, cv_inod)
        end do

        !Obtain the conductivity of the saturated rock
        call get_rock_sat_conductivity(state, packed_state, Mdims, ndgln, Saturation%val(1, 1, :), F_fields(2, 1, :), Temperature%val(1, 1, :), rock_sat_conductivity(1,:))
        !Compute K_fields
        do k = 1, nfields
          call get_SP_coupling_coefficients(state, packed_state, Mdims, ndgln, rock_sat_conductivity(1,:), K_fields(k,1,:), &
                      Saturation%val(1, 1, :), F_fields(2,1,:), Temperature%val(1, 1, :), flag = k )
        end do

        !Solver options
        solver_option_path = "/solver_options/Linear_solver"
        if (have_option('/solver_options/Linear_solver/Custom_solver_configuration/field::SPSolver')) then
          solver_option_path = '/solver_options/Linear_solver/Custom_solver_configuration/field::SPSolver'
        end if
        ! SP Solver elements
        SelfPotential => extract_scalar_field(state(1),"SelfPotential", stat)
        call generate_and_solve_Laplacian_system( Mdims, state, packed_state, ndgln, Mmat, Mspars, CV_funs, CV_GIdims, &
                                      rock_sat_conductivity, "SelfPotential", K_fields, F_fields, 20, solver_option_path)

        !##########Now we normalise the SP result to have the reference node with voltage = 0. We do this because is better to remove the null space###########
        !##Retrieve the coordinates of the reference position##
        reference_value = 0.
        call find_reference_node_from_coordinates(X_ALL, Saturation%mesh,"/porous_media/SelfPotential",reference_nod,reference_node_owned)
        !The processor that owns the node retrieves the value
        if (IsParallel()) then
          if (reference_node_owned) reference_value = SelfPotential%val(reference_nod)
        else
          reference_value = SelfPotential%val(reference_nod)
        end if
        !Share the value between all the processors
        call allsum(reference_value)
        !Apply the reference to ensure that the reference node is zero
        SelfPotential%val = (SelfPotential%val - reference_value)
        deallocate(rock_sat_conductivity, F_fields, K_fields)
      end subroutine Assemble_and_solve_SP


      !>@brief: This subroutine computes the saturated rock conductivity based on Tiab and Donaldson, 2004 formula.
      !> and also includes the temperature/concentration dependency detailed in Sen and Goode, 1992
      !> IMPORTANT: Water needs to be phase 1! TODO Maybe add a REMINDER when running with SP solver and to run with Kelvin, maybe this for all the temperature cases
      subroutine get_rock_sat_conductivity(state, packed_state, Mdims, ndgln, Saturation, Concentration, Temperature, rock_sat_conductivity )
        implicit none
        type( state_type ), dimension(:), intent( inout ) :: state
        type(multi_dimensions), intent( in ) :: Mdims
        type( state_type ), intent( in ) :: packed_state
        type(multi_ndgln), intent(in) :: ndgln
        real, dimension(:), intent(out) :: rock_sat_conductivity!output of the subroutine.
        real, dimension(:), intent(in) :: Concentration, Saturation, Temperature!Here Concentration needs to be in mol/litre
        !Local varibales
        integer:: cv_inod, ele, cv_iloc, stat, ele_pore
        real :: auxR, temp
        type(vector_field), pointer :: porosity
        real :: cementation_exp !This I presume should be assigned from diamond?
        real :: sat_exp !Saturation exponent, again I presume should be assigned from diamond?
        real, parameter :: Kelv_conv = 273.15
        real, dimension(Mdims%cv_nonods) :: water_conductivity, cv_counter
        logical, save :: show_msg = .true.
        real, parameter :: tol = 1e-8

        !If using python code all the problem are the users
        if (have_option("/porous_media/SelfPotential/python_Rock_sat_conductivity_code")) then
          call multi_compute_python_field(state, 1, "/porous_media/SelfPotential/python_Rock_sat_conductivity_code", rock_sat_conductivity)
        else
          !Retrieve fields from state/packed_state
          !Check the situation with the temperature field/value
          if (has_temperature) then
            if (show_msg) then
              if (any(Temperature - Kelv_conv < 0.)) then
                ewrite(0, *) "REMINDER: The S.I. units for TEMPERATURE are Kelvin not Celsius."
              end if
            end if
            show_msg = .false.
          else
            if (have_option("/porous_media/SelfPotential/Reservoir_temperature")) then
              call get_option( '/porous_media/SelfPotential/Reservoir_temperature', temp )
              if (show_msg) then
                if (GetProcNo() == 1 .and. temp - Kelv_conv < 0.) then
                  ewrite(0, *) "REMINDER: The S.I. units for TEMPERATURE are Kelvin not Celsius."
                end if
              end if
              show_msg = .false.
            else
              if (GetProcNo() == 1) then
                ewrite(0, *) "ERROR: SelfPotential requires to define a temperature either by a prognostic field or using the option: /porous_media/SelfPotential/Reservoir_temperature."
                ewrite(0, *) "SelfPotential will NOT be computed."
                return
              end if
            end if
          end if
          !Retrieve exponents
          call get_option("/porous_media/SelfPotential/Cementation_exp",cementation_exp, default = 1.8 )
          call get_option("/porous_media/SelfPotential/Sat_exponent",sat_exp, default = 2.0 )

          porosity=>extract_vector_field(packed_state,"Porosity")
          if (.not. has_concentration) then
            if (GetProcNo() == 1) then
              ewrite(0, *) "ERROR: For Self Potential calculation a concentration field is required."
              ewrite(0, *) "SelfPotential will NOT be computed."
              return
            end if
          else !Compute water conductivity based on Temperature and concentration (Sen and Goode, 1992)
            do cv_inod = 1, Mdims%cv_nonods
              if (has_temperature) temp = Temperature(cv_inod)                   !Water concentration
              water_conductivity(cv_inod) = (5.6 + 0.27 * temp - 1.5e-4 * temp**2.)*(Concentration(cv_inod)+tol) &
              - Concentration(cv_inod)**1.5 * ( 2.36 + 0.099 * temp) / (1 + 0.214 * Concentration(cv_inod)**0.5)
            end do
          end if

          cv_counter = 0.; rock_sat_conductivity = 0.
          !Now compute rock_saturated conductivity
          do  ele = 1, Mdims%totele
            do cv_iloc = 1, Mdims%cv_nloc
              cv_inod = ndgln%cv( ( ele - 1 ) * Mdims%cv_nloc + cv_iloc )                                       !Only the water phase
              rock_sat_conductivity(cv_inod) = rock_sat_conductivity(cv_inod) + porosity%val(1,ele) ** cementation_exp * water_conductivity(cv_inod) * (Saturation(cv_inod)+tol) ** sat_exp
              cv_counter( cv_inod ) = cv_counter( cv_inod ) + 1.0
            end do
          end do
          !Average of rock_sat_conductivity
          rock_sat_conductivity = rock_sat_conductivity/cv_counter
        end if

      end subroutine get_rock_sat_conductivity

      !>@brief: This subroutine computes the coupling coefficients required to compute the self potential using Jackson et al. (2012)
      !>1 => Electrokinetic; 2=> Thermal; 3=> Exclusion diffusion
      subroutine get_SP_coupling_coefficients(state, packed_state, Mdims, ndgln, rock_sat_conductivity, coupling_term, Saturation, Concentration, Temperature, flag )
        implicit none
        type( state_type ), dimension(:), intent( inout ) :: state
        type(multi_dimensions), intent( in ) :: Mdims
        type( state_type ), intent( inout ) :: packed_state
        real, dimension(:), intent(in) :: rock_sat_conductivity
        type(multi_ndgln), intent(in) :: ndgln
        real, dimension(:), intent(inout) ::  coupling_term!>output of the subroutine.
        real, dimension(:), intent(in) :: Concentration, Saturation, Temperature!Here Concentration needs to be in mol/litre
        integer, intent(in) :: flag !>1 => Electrokinetic; 2=> Thermal; 3=> Exclusion diffusion
        !Local variables
        real :: norm_water_sat, Cf, Tna, AuxR, temp
        integer :: cv_inod, ele, cv_iloc, stat, imat
        real, dimension(:, :), pointer :: Immobile_fraction
        real, dimension(Mdims%cv_nonods) :: cv_counter, coupling_coef, coupling_coef_ee, coupling_coef_ed
        real, parameter :: EK_exp = 0.6 !From Jackson et al 2012
        real, parameter :: Tol = 1e-8
        logical :: post_process

        !Retrieve fields from state/packed_state
        if (have_option("/porous_media/SelfPotential/Reservoir_temperature")) then
          call get_option( '/porous_media/SelfPotential/Reservoir_temperature', temp )
        end if

        call get_var_from_packed_state(packed_state, Immobile_fraction = Immobile_fraction)
        post_process = .false.
        coupling_coef = 0.
        cv_counter = 0.
        !Electrokinetic coupling coefficient
        if (flag == 1) then
          if (have_option("/porous_media/SelfPotential/python_ElectroKinetic_code")) then
            call multi_compute_python_field(state, 1, "/porous_media/SelfPotential/python_ElectroKinetic_code", coupling_coef)
          else
            do  ele = 1, Mdims%totele
              do cv_iloc = 1, Mdims%cv_nloc
                cv_inod = ndgln%cv( ( ele - 1 ) * Mdims%cv_nloc + cv_iloc )
                IMAT = ndgln%mat( ( ELE - 1 ) * Mdims%mat_nloc + CV_ILOC )
                !Obtain normalised saturation
                norm_water_sat = (Saturation(cv_inod) - Immobile_fraction(1, imat)) / (1.0 - sum(Immobile_fraction(1:Mdims%n_in_pres, imat)))
                coupling_coef(cv_inod) = coupling_coef(cv_inod) + (-1.36 * (Concentration(cv_inod)+tol)**-0.9123 * 1e-9 ) * norm_water_sat ** EK_exp!Not sure if exponent or times...
                ! coupling_coef(cv_inod) = coupling_coef(cv_inod) + 2.5e-9!<=I think this was used for Mutlaq et al 2019
                cv_counter( cv_inod ) = cv_counter( cv_inod ) + 1.0
              end do
            end do
            coupling_coef = coupling_coef/cv_counter
          end if
        end if
        !Exclusion diffusion coefficient
        if (flag == 2 .and. .not. have_option("/porous_media/SelfPotential/python_Electrodiffusive_code")) then
          post_process = .true.
            do cv_inod = 1, Mdims%cv_nonods
              Tna = get_Hittorf_transport_number(Concentration(cv_inod))
              Cf = Concentration(cv_inod) + Tol!To avoid divisions by zero
              if (has_temperature) temp = Temperature(cv_inod)
              coupling_coef_ed(cv_inod) = - 8.61e-2 * (2.*Tna - 1) * temp/Cf
              coupling_coef_ee(cv_inod) = - 8.61e-2 * temp/Cf
            end do
        end if
        !Thermal coupling coefficient
        if (flag == 3 .and. .not. have_option("/porous_media/SelfPotential/python_Thermoelectric_code")) then
          post_process = .true.
          do cv_inod = 1, Mdims%cv_nonods
            Tna = get_Hittorf_transport_number(Concentration(cv_inod))
            AuxR = LOG(Concentration(cv_inod) + Tol )!To avoid reaching zero
            coupling_coef_ed(cv_inod) = - 1.984e-1*(2.*Tna - 1.) * AuxR + 1.059 * Tna - 5.673e-1
            coupling_coef_ee(cv_inod) = - 1.984e-1 * AuxR + 5.953e-1
          end do
        end if
        !For thermal and diffusion-exclusion we need to combine them based on the normalised saturation
        if (post_process) then
          do ele = 1, Mdims%totele
            do cv_iloc = 1, Mdims%cv_nloc
              cv_inod = ndgln%cv( ( ele - 1 ) * Mdims%cv_nloc + cv_iloc )
              norm_water_sat = (Saturation(cv_inod) - Immobile_fraction(1, ele)) / (1.0 - sum(Immobile_fraction(1:Mdims%n_in_pres, ele)))
              coupling_coef(cv_inod) = coupling_coef(cv_inod) + (1 - norm_water_sat)**3. * (coupling_coef_ee(cv_inod) - coupling_coef_ed(cv_inod)) + coupling_coef_ed(cv_inod)
              cv_counter( cv_inod ) = cv_counter( cv_inod ) + 1.0
            end do
          end do
          !Obtain the average since we have overlooped cv nodes
          coupling_coef = coupling_coef/cv_counter *1e-3 !To convert from mV to Volts only flags 2 and 3
        end if

        !#############Using python######################
        !Thermal coupling coefficient
        if (flag == 2 .and. have_option("/porous_media/SelfPotential/python_Electrodiffusive_code")) then
          call multi_compute_python_field(state, 1, "/porous_media/SelfPotential/python_Electrodiffusive_code", coupling_coef)
        end if
        if (flag == 3 .and. have_option("/porous_media/SelfPotential/python_Thermoelectric_code")) then
          call multi_compute_python_field(state, 1, "/porous_media/SelfPotential/python_Thermoelectric_code", coupling_coef)
        end if
        !Finally obtain the coupling coefficient
        coupling_term = coupling_coef * rock_sat_conductivity
      contains
        !>@brief: Compute the macroscopic Hittorf transport number for the positive Sodium ions
        real function get_Hittorf_transport_number(Concentration_water)
          implicit none
          real, intent(in) :: Concentration_water
          if (Concentration_water < 0.09) then
            get_Hittorf_transport_number = 0.39
          else
            get_Hittorf_transport_number = 3.66e-1 - 2.12e-2 * LOG10(Concentration_water)
          end if

        end function get_Hittorf_transport_number

      end subroutine get_SP_coupling_coefficients



end module multi_SP
