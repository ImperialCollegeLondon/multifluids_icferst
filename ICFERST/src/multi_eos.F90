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

module multiphase_EOS

    use fldebug
    use unittest_tools
    use state_module
    use fields
    use global_parameters
    use spud
    use futils, only: int2str
    use vector_tools
    use python_state
    use Copy_Outof_State
    use multi_tools !!!!! WHY is the bad_elements type not being picked up by Copy_Outof_State??????
    use shape_functions_Linear_Quadratic
    use sparse_tools
    use Multiphase_module
    use sparsity_patterns_meshes, only: get_csr_sparsity_firstorder
    use boundary_conditions, only: get_entire_boundary_condition
    use Field_Options, only: get_external_coordinate_field
    use initialise_fields_module, only: initialise_field_over_regions, initialise_field
    use multi_tools, only: CALC_FACE_ELE, assign_val, table_interpolation, read_csv_table
    use checkpoint


    implicit none

    real, parameter :: flooding_hmin = 1e-5


contains

    !>@brief: Computes the density for the phases and the derivatives of the density
    subroutine Calculate_All_Rhos( state, packed_state, Mdims, get_RhoCp )

        implicit none

        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent( in ) :: Mdims
        type(multi_ndgln) :: ndgln
        logical, optional, intent(in) :: get_RhoCp !This flags computes rhoCp instead of rho
        !Local variables
        integer, dimension( : ), pointer :: cv_ndgln
        integer :: ncomp_in, nphase, ndim, cv_nonods, cv_nloc, totele
        real, dimension( : ), allocatable :: Rho, dRhodP, rho_porous, drhodp_porous, &
          Density_Bulk, DensityCp_Bulk, Density_Component, Cp, Component_l, c_cv_nod
        character( len = option_path_len ), dimension( : ), allocatable :: eos_option_path
        type( tensor_field ), pointer :: PackedDRhoDPressure ! (nphase, cv_nonods)
        type( tensor_field ), pointer :: field1, field2, field3, field4
        type( scalar_field ), pointer :: Cp_s, Density
        integer :: icomp, iphase, ncomp, sc, ec, sp, ep, ip, stat, cv_iloc, cv_nod, ele
        logical :: compute_rhoCP
        logical, parameter :: harmonic_average=.false.


        !Only obtain RhoCP if CP is defined and when solving for Temperature, else, return Rho
        compute_rhoCP = .false.
        if (present_and_true(get_RhoCp)) then
          compute_rhoCP = have_option("/material_phase[0]/phase_properties/scalar_field::HeatCapacity")
        end if
        ncomp_in = Mdims%ncomp ; nphase = Mdims%nphase ; ndim = Mdims%ndim
        cv_nonods = Mdims%cv_nonods ; cv_nloc = Mdims%cv_nloc ; totele = Mdims%totele
        cv_ndgln => get_ndglno( extract_mesh( state( 1 ), "PressureMesh" ) )

        PackedDRhoDPressure => extract_tensor_field( packed_state, "PackedDRhoDPressure" )
        PackedDRhoDPressure%val = 0.

        ncomp = ncomp_in
        if( ncomp_in == 0 ) ncomp = 1

        allocate( eos_option_path( nphase * ncomp ) )

        if( ncomp > 1 ) then
           do icomp =1, ncomp
              do iphase =1, nphase
                 eos_option_path( ( icomp - 1 ) * nphase + iphase ) = &
                      trim( '/material_phase[' // int2str( nphase + icomp - 1 ) // &
                      ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
                      '/prognostic/phase_properties/Density' )
                 call Assign_Equation_of_State( eos_option_path( ( icomp - 1 ) * nphase + iphase ) )
              end do
           end do
        else
           do iphase = 1, nphase
              eos_option_path( iphase ) = trim( '/material_phase[' // int2str( iphase - 1 ) // ']/phase_properties/Density' )
              call Assign_Equation_of_State( eos_option_path( iphase ) )
           end do
        end if

        allocate( Rho( cv_nonods ), dRhodP( cv_nonods ) )
        if (compute_rhoCP) then
          allocate( Cp( cv_nonods ) ) ; Cp = 1.0
          allocate( DensityCp_Bulk( nphase * cv_nonods ) ); DensityCp_Bulk = 0.0
        end if
        allocate( Density_Component( ncomp * nphase * cv_nonods ) )
        allocate( Density_Bulk( nphase * cv_nonods ) ); Density_Bulk = 0.0

        allocate( Component_l( cv_nonods ) ) ; Component_l = 0.
        allocate( drhodp_porous( cv_nonods ) ); drhodp_porous = 0.
        !If boussinesq then we do not consider variations in rock density either
        if (.not. has_boussinesq_aprox) then
          if (have_option('/porous_media/porous_properties/scalar_field::porous_compressibility')) &
            call Calculate_porous_Rho_dRhoP(state,packed_state,Mdims, cv_ndgln, drhodp_porous  )
        end if
        do icomp = 1, ncomp
           do iphase = 1, nphase
              sc = ( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1
              ec = ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods

              sp = ( iphase - 1 ) * cv_nonods + 1
              ep = iphase * cv_nonods

              Rho=0. ; dRhodP=0.
              call Calculate_Rho_dRhodP( state, packed_state, iphase, icomp, &
                   nphase, ncomp_in, eos_option_path( (icomp - 1 ) * nphase + iphase ), Rho, dRhodP )
              if ( ncomp > 1 ) then
                 field4 => extract_tensor_field( packed_state, "PackedComponentMassFraction" )
                 Component_l = field4 % val ( icomp, iphase, :)

                 if ( have_option( '/material_phase[0]/linearise_component' ) ) then
                    ! linearise component
                    if ( cv_nloc==6 .or. (cv_nloc==10 .and. ndim==3) ) then ! P2 triangle or tet
                       allocate( c_cv_nod( cv_nloc ) )
                       do ele = 1, totele
                          do cv_iloc = 1, cv_nloc
                             cv_nod = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )
                             c_cv_nod( cv_iloc ) = Component_l( cv_nod )
                          end do

                          c_cv_nod( 2 ) = 0.5 * ( c_cv_nod ( 1 ) + c_cv_nod( 3 ) )
                          c_cv_nod( 4 ) = 0.5 * ( c_cv_nod ( 1 ) + c_cv_nod( 6 ) )
                          c_cv_nod( 5 ) = 0.5 * ( c_cv_nod ( 3 ) + c_cv_nod( 6 ) )

                          if ( cv_nloc==10 ) then
                             c_cv_nod ( 7 ) = 0.5 * ( c_cv_nod ( 1 ) + c_cv_nod( 10 ) )
                             c_cv_nod ( 8 ) = 0.5 * ( c_cv_nod ( 3 ) + c_cv_nod( 10 ) )
                             c_cv_nod ( 9 ) = 0.5 * ( c_cv_nod ( 6 ) + c_cv_nod( 10 ) )
                          end if

                          do cv_iloc = 1, cv_nloc
                             cv_nod = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )
                             Component_l( cv_nod ) = c_cv_nod( cv_iloc )
                          end do
                       end do
                       deallocate( c_cv_nod )
                    end if
                 end if

                 if ( .not.harmonic_average ) then

                    ! rho = rho +  a_i * rho_i
                    Density_Bulk( sp : ep ) = Density_Bulk( sp : ep ) + Rho * Component_l
                    if (has_boussinesq_aprox) then !disable time-derivative terms
                      PackedDRhoDPressure%val( 1, iphase, : ) = 0.
                    else
                      PackedDRhoDPressure%val( 1, iphase, : ) = PackedDRhoDPressure%val( 1, iphase, : ) + dRhodP * Component_l / Rho +drhodp_porous
                    end if
                    Density_Component( sc : ec ) = Rho

                    Cp_s => extract_scalar_field( state( nphase + icomp ), &
                         'ComponentMassFractionPhase' // int2str( iphase ) // 'HeatCapacity', stat )
                    if ( stat == 0 .and. compute_rhoCP ) then
                      call assign_val(Cp,Cp_s % val)!Cp = Cp_s % val
                      DensityCp_Bulk( sp : ep ) = DensityCp_Bulk( sp : ep ) + Rho * Cp * Component_l
                    end if

                 else
                    Density_Bulk( sp : ep ) = Density_Bulk( sp : ep ) + Rho * Component_l
                    if (has_boussinesq_aprox) then!disable time-derivative terms
                      PackedDRhoDPressure%val( 1, iphase, : ) = 0.
                    else
                      PackedDRhoDPressure%val( 1, iphase, : ) = PackedDRhoDPressure%val( 1, iphase, : ) + dRhodP * Component_l / Rho + drhodp_porous
                    end if
                    Density_Component( sc : ec ) = Rho

                    ! harmonic average
                    ! rho = rho + 1.0 / ( a_i / rho_i )
                    Cp_s => extract_scalar_field( state( nphase + icomp ), &
                         'ComponentMassFractionPhase' // int2str( iphase ) // 'HeatCapacity', stat )
                    if ( stat == 0 .and. compute_rhoCP) then
                      call assign_val(Cp,Cp_s % val)!Cp = Cp_s % val
                      do cv_nod = 1, cv_nonods
                         ip = ( iphase - 1 ) * cv_nonods + cv_nod
                         DensityCp_Bulk( ip ) = DensityCp_Bulk( ip ) + Component_l(cv_nod) / ( Rho(cv_nod) * Cp(cv_nod) )
                      end do
                    end if
                 end if

              else
                Density_Bulk( sp : ep ) = Rho
                if (has_boussinesq_aprox) then !disable time-derivative terms
                  PackedDRhoDPressure%val( 1, iphase, : ) = 0.
                else
                 PackedDRhoDPressure%val( 1, iphase, : ) = dRhodP + drhodp_porous
                end if
                 Cp_s => extract_scalar_field( state( iphase ), 'TemperatureHeatCapacity', stat )
                 !Cp_s => extract_scalar_field( state( iphase ), 'ConcentrationHeatCapacity', stat )
                 if ( stat == 0 .and. compute_rhoCP) then
                   call assign_val(Cp,Cp_s % val)
                   DensityCp_Bulk( sp : ep ) = Rho * Cp
                 end if
              end if

           end do ! iphase
        end do ! icomp

        if ( ncomp > 1 ) then
           if ( harmonic_average .and. compute_rhoCP) DensityCp_Bulk = 1.0 / DensityCp_Bulk
           call Cap_Bulk_Rho( state, ncomp, nphase, &
                cv_nonods, Density_Component, Density_Bulk, DensityCp_Bulk, compute_rhoCP )
        end if

        field1 => extract_tensor_field( packed_state, "PackedDensity" )
        if (compute_rhoCP) field2 => extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )
        if( ncomp > 1 ) field3 => extract_tensor_field( packed_state, "PackedComponentDensity" )

        do iphase = 1, nphase
           sp = ( iphase - 1 ) * cv_nonods + 1
           ep = iphase * cv_nonods

           field1 % val ( 1, iphase, : ) = Density_Bulk( sp : ep )
           if (compute_rhoCP) field2 % val ( 1, iphase, : ) = DensityCp_Bulk( sp : ep )

           if ( ncomp > 1 ) then
              do icomp = 1, ncomp
                 sc = ( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1
                 ec = ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods
                 field3 % val ( icomp, iphase, : ) = Density_Component( sc : ec )
              end do ! icomp
           end if
        end do ! iphase

        ! if (has_boussinesq_aprox .and. .not. is_porous_media) field2 % val = 1.0
        deallocate( Rho, dRhodP, Component_l)
        if (allocated(drhodp_porous)) deallocate(drhodp_porous)
        deallocate( Density_Component, Density_Bulk )
        deallocate( eos_option_path )
        if (compute_rhoCP) deallocate(Cp, DensityCp_Bulk)
        !sprint_to_do copying meory to itself...
        do iphase = 1, nphase
           Density => extract_scalar_field( state( iphase ), "Density" )
           Density % val = field1 % val( 1, iphase, : )
        end do

    end subroutine Calculate_All_Rhos

    !>@brief: Computes the bulk density for components and the derivatives of the density
    subroutine Cap_Bulk_Rho( state, ncomp, nphase, &
        cv_nonods, Density_Component, Density, Density_Cp , compute_rhoCP)

        implicit none

        type(state_type), dimension( : ) :: state
        integer, intent( in ) :: nphase, ncomp, cv_nonods
        real, dimension( cv_nonods * nphase * ncomp ), intent( in ) :: Density_Component
        real, dimension( cv_nonods * nphase ), intent( inout ) :: Density, Density_Cp
        logical, intent(in) :: compute_rhoCP

        real, dimension( :, : ), allocatable :: Density_Component_Min, Density_Component_Max
        real, dimension( :, : ), allocatable :: Density_Cp_Component_Min, Density_Cp_Component_Max
        type( scalar_field ), pointer :: Cp_s
        real, dimension( : ), allocatable :: Cp
        integer :: sp, ep, sc, ec, iphase, icomp, stat

        allocate( Density_Component_Min( nphase, cv_nonods ) ) ; Density_Component_Min = 1.0e+15
        allocate( Density_Component_Max( nphase, cv_nonods ) ) ; Density_Component_Max = 0.0
        if (compute_rhoCP) then
          allocate( Density_Cp_Component_Min( nphase, cv_nonods ) ) ; Density_Cp_Component_Min = 1.0e+15
          allocate( Density_Cp_Component_Max( nphase, cv_nonods ) ) ; Density_Cp_Component_Max = 0.0
          allocate( Cp( cv_nonods ) ) ; Cp = 1.0
        end if
        do iphase = 1, nphase
            do icomp = 1, ncomp
                sc = ( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1
                ec = ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods

                Density_Component_Min( iphase, : ) = min( Density_Component_Min( iphase, : ), Density_Component( sc : ec ) )
                Density_Component_Max( iphase, : ) = max( Density_Component_Max( iphase, : ), Density_Component( sc : ec ) )
                if (compute_rhoCP) then
                  Cp = 1.0
                  Cp_s => extract_scalar_field( state( nphase + icomp ), &
                      'ComponentMassFractionPhase' // int2str( iphase ) // 'HeatCapacity', stat )
                  if( stat == 0 ) Cp = Cp_s % val

                  Density_Cp_Component_Min( iphase, : ) = min( Density_Cp_Component_Min( iphase, : ), Density_Component( sc : ec ) * Cp )
                  Density_Cp_Component_Max( iphase, : ) = max( Density_Cp_Component_Max( iphase, : ), Density_Component( sc : ec ) * Cp )
                end if
            end do
        end do

        do iphase = 1, nphase
            sp = ( iphase - 1 ) * cv_nonods + 1
            ep = iphase * cv_nonods

            Density( sp : ep ) = min( Density( sp : ep ), Density_Component_Max( iphase, : ) )
            Density( sp : ep ) = max( Density( sp : ep ), Density_Component_Min( iphase, : ) )
            if (compute_rhoCP) then
              Density_Cp( sp : ep ) = min( Density_Cp( sp : ep ), Density_Cp_Component_Max( iphase, : ) )
              Density_Cp( sp : ep ) = max( Density_Cp( sp : ep ), Density_Cp_Component_Min( iphase, : ) )
            end if
        end do
        if (compute_rhoCP) then
          deallocate( Cp )
          deallocate( Density_Cp_Component_Min, Density_Cp_Component_Max )
        end if
        deallocate( Density_Component_Min, Density_Component_Max )

    end subroutine Cap_Bulk_Rho

    !>@brief: Computes the density for the components and the derivatives of the density
    !> This is were the action is actually done, the other calculate rhos are wrappers
    subroutine Calculate_Component_Rho( state, packed_state, Mdims )

      implicit none

      type( state_type ), dimension( : ), intent( inout ) :: state
      type( state_type ), intent( inout ) :: packed_state
      type( multi_dimensions ), intent( in ) :: Mdims

      real, dimension( : ), allocatable :: Rho, dRhodP
      type( tensor_field ), pointer :: field
      character( len = option_path_len ) :: eos_option_path
      integer :: icomp, iphase, s, e


      allocate( Rho( Mdims%cv_nonods ), dRhodP( Mdims%cv_nonods ) )

      field => extract_tensor_field( packed_state, "PackedComponentDensity" )

      do icomp = 1, Mdims%ncomp
         do iphase = 1, Mdims%nphase
            s = ( icomp - 1 ) * Mdims%nphase * Mdims%cv_nonods + ( iphase - 1 ) * Mdims%cv_nonods + 1
            e = ( icomp - 1 ) * Mdims%nphase * Mdims%cv_nonods + iphase * Mdims%cv_nonods
            eos_option_path = trim( '/material_phase[' // int2str( Mdims%nphase + icomp - 1 ) // &
                 ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
                 '/prognostic/phase_properties/Density' )
            call Assign_Equation_of_State( eos_option_path )
            Rho=0. ; dRhodP=0.
            call Calculate_Rho_dRhodP( state, packed_state, iphase, icomp, &
                 Mdims%nphase, Mdims%ncomp, eos_option_path, Rho, dRhodP )

            field % val( icomp, iphase, : ) = Rho
         end do ! iphase
      end do ! icomp
      deallocate( Rho, dRhodP )

    end subroutine Calculate_Component_Rho

    !>@brief: Computes the density and the derivative of the density
    subroutine Calculate_Rho_dRhodP( state, packed_state, iphase, icomp, &
        nphase, ncomp, eos_option_path, rho, drhodp )

        implicit none

        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state

        integer, intent( in ) :: iphase, icomp, nphase, ncomp
        integer :: JWLn, JWLi
        character( len = option_path_len ), intent( in ) :: eos_option_path
        real, dimension( : ), intent( inout ) :: rho, drhodp
        real, dimension( : ), allocatable :: ro0

        type( tensor_field ), pointer :: pressure
        type( scalar_field ), pointer :: temperature, density, Concentration
        character( len = option_path_len ) :: option_path_comp, option_path_incomp, option_path_python, buffer, option_name
        character( len = python_func_len ) :: pycode
        logical, save :: initialised = .false.
        logical :: have_temperature_field
        logical :: have_concentration_field
        real, parameter :: toler = 1.e-10
        real, dimension( : ), allocatable, save :: reference_pressure
        real, dimension( : ), allocatable :: eos_coefs, perturbation_pressure, RhoPlus, RhoMinus
        real, dimension( : ), allocatable :: pressure_back_up, density_back_up, temperature_local
        real :: dt, current_time
        integer :: ncoef, stat, nfields, ifield
        !Variables for python function for the coefficient_B for linear density (this is for bathymetry)
        type (scalar_field) :: sfield
        type (scalar_field), pointer :: pnt_sfield
        type (vector_field), pointer :: position

        !!$ Den = c1 * ( P + c2 ) / T           :: Stiffened EOS
        !!$ Den = c1 * P + c2                   :: Linear_1 EOS
        !!$ Den = c1 * P / T + c2               :: Linear_2 EOS
        !!$ Den = Den0 * exp[ c0 * ( P - P0 ) ] :: Exponential_1 EOS
        !!$ Den = c0 * P** c1                   :: Exponential_2 EOS

        pressure => extract_tensor_field( packed_state, 'PackedCVPressure', stat )
        if (stat/=0) pressure => extract_tensor_field( packed_state, "PackedFEPressure", stat )

        temperature => extract_scalar_field( state( iphase ), 'Temperature', stat )
        have_temperature_field = ( stat == 0 )
        Concentration => extract_scalar_field( state( iphase ), 'Concentration', stat )
        have_concentration_field = ( stat == 0 )

        assert( node_count( pressure ) == size( rho ) )
        assert( node_count( pressure ) == size( drhodp ) )

        allocate( perturbation_pressure( node_count( pressure ) ) ) ; perturbation_pressure = 0.
        allocate( RhoPlus( node_count( pressure ) ) ) ; RhoPlus = 0.
        allocate( RhoMinus( node_count( pressure ) ) ) ; RhoMinus = 0.

        if ( ncomp > 0 ) then
            option_path_comp = trim( '/material_phase[' // int2str( nphase + icomp - 1 ) // &
                ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
                '/prognostic/phase_properties/Density/compressible' )
            option_path_incomp = trim( '/material_phase[' // int2str(nphase + icomp - 1 ) // &
                ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
                '/prognostic/phase_properties/Density/incompressible' )
            option_path_python = trim( '/material_phase[' // int2str( nphase + icomp - 1 ) // &
                ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
                '/prognostic/phase_properties/Density/python_state' )
        else
            option_path_comp = trim( '/material_phase[' // int2str( iphase - 1 ) // &
                ']/phase_properties/Density/compressible' )
            option_path_incomp = trim( '/material_phase[' // int2str( iphase - 1 ) // &
                ']/phase_properties/Density/incompressible' )
            option_path_python = trim( '/material_phase[' // int2str( iphase - 1 ) // &
                ']/phase_properties/Density/python_state' )
        end if

        Conditional_EOS_Option: if( trim( eos_option_path ) == trim( option_path_incomp ) ) then
            !!$ Constant representation
            allocate( temperature_local( node_count( pressure ) ) ) ; temperature_local = 0.
            if ( have_temperature_field ) temperature_local = temperature % val
            ncoef = 10 ; allocate( eos_coefs( ncoef ) ) ; eos_coefs = 0.
            call get_option( trim( eos_option_path ), eos_coefs( 1 ) )
            eos_coefs( 2 : 10 ) = 0.
            call Density_Polynomial( eos_coefs, pressure % val(1,1,:), temperature_local, &
                Rho )
            perturbation_pressure = max( toler, 1.e-3 * abs( pressure % val(1,1,:) ) )
            call Density_Polynomial( eos_coefs, pressure % val(1,1,:) + perturbation_pressure, temperature_local, &
                RhoPlus )
            call Density_Polynomial( eos_coefs, pressure % val(1,1,:) - perturbation_pressure, temperature_local, &
                RhoMinus )
            dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure
            deallocate( temperature_local, eos_coefs )

      else if( trim( eos_option_path ) == trim( option_path_comp ) // '/stiffened_gas' ) then
            !!$ Den = C0 / T * ( P - C1 )
            if( .not. have_temperature_field ) FLAbort( 'Temperature Field not defined' )
            allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
            call get_option( trim( eos_option_path) // '/eos_option1' , eos_coefs( 1 ) )
            call get_option( trim( eos_option_path )// '/eos_option2' , eos_coefs( 2 ) )
            Rho = ( pressure%val(1, 1, :) + eos_coefs( 1 ) ) * eos_coefs( 2 ) / temperature % val
            perturbation_pressure = max( toler, 1.e-3 * ( abs( pressure%val(1, 1, :) ) + eos_coefs( 1 ) ) )
            RhoPlus = ( pressure%val(1, 1, :) + perturbation_pressure + eos_coefs( 1 ) ) *  eos_coefs( 2 ) / &
                temperature % val
            RhoMinus = ( pressure%val(1 , 1, :) - perturbation_pressure + eos_coefs( 1 ) ) *  eos_coefs( 2 ) / &
                temperature % val
            if (eos_coefs(1)==0) then
              dRhodP = eos_coefs(2)/temperature%val
            else
              dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure
            end if
            deallocate( eos_coefs )


          elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/linear_in_pressure' ) then
            !!$ Den = C0 * P +C1
            allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
            !By default the pressure mesh (position 1)
            pnt_sfield => extract_scalar_field(state(1),1)
            position => get_external_coordinate_field(packed_state, pnt_sfield%mesh)
            call allocate (sfield, pnt_sfield%mesh, "Temporary_linear_Coefficient_B")
            !Retrieve coefficients
            call get_option( trim( eos_option_path ) // '/coefficient_A', eos_coefs( 1 ) )
            call initialise_field(sfield, trim( option_path_comp )//"/linear_in_pressure/coefficient_B" , position)
            Rho = eos_coefs( 1 ) * pressure % val(1,1,:) + sfield%val
            perturbation_pressure = 1.
            !RhoPlus = eos_coefs( 1 ) * ( pressure % val + perturbation_pressure ) + eos_coefs( 2 )
            !RhoMinus = eos_coefs( 1 ) * ( pressure % val - perturbation_pressure ) + eos_coefs( 2 )
            dRhodP = eos_coefs( 1 ) !0.5 * ( DensityPlus - DensityMinus ) / perturbation_pressure
            deallocate( eos_coefs )
            call deallocate(sfield)


          elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/linear_in_pressure/include_internal_energy' ) then
             !!$ Den = C0 * P/T +C1
             if( .not. have_temperature_field ) FLAbort( 'Temperature Field not defined' )
             allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
             call get_option( trim( option_path_comp ) // '/linear_in_pressure/coefficient_A', eos_coefs( 1 ) )
             call get_option( trim( option_path_comp ) // '/linear_in_pressure/coefficient_B/constant', eos_coefs( 2 ) )


             Rho = eos_coefs( 1 ) * pressure % val(1,1,:) / temperature % val + eos_coefs( 2 )
             perturbation_pressure = 1.
             RhoPlus = eos_coefs( 1 ) * ( pressure % val(1,1,:) + perturbation_pressure ) / &
                  ( max( toler, temperature % val ) ) + eos_coefs( 2 )
             RhoMinus = eos_coefs( 1 ) * ( pressure % val(1,1,:) - perturbation_pressure ) / &
                  ( max( toler, temperature % val ) ) + eos_coefs( 2 )
             !dRhodP =  eos_coefs( 1 ) / temperature % val
             dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure
             deallocate( eos_coefs )


        elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/exponential_in_pressure' ) then
            !!$ Den = C0 * ( P ^ C1 )
            allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
            call get_option( trim( eos_option_path ) // '/coefficient_A', eos_coefs( 1 ) )
            call get_option( trim( eos_option_path ) // '/coefficient_B', eos_coefs( 2 ) )
             Rho = eos_coefs( 1 ) * pressure % val(1,1,:) ** eos_coefs( 2 )
             perturbation_pressure = 1.
             RhoPlus = eos_coefs( 1 ) * ( pressure % val(1,1,:) + perturbation_pressure ) ** eos_coefs( 2 )
             RhoMinus = eos_coefs( 1 ) * ( pressure % val(1,1,:) - perturbation_pressure ) ** eos_coefs( 2 )
             dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure
            deallocate( eos_coefs )

          elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/Linear_eos' ) then
              !!$ Den = den0 * ( 1 + alpha * DeltaC - beta * DeltaT + gamma * DeltaP)

              allocate( eos_coefs( 7 ) ) ; eos_coefs = 0.
              call get_option( trim( eos_option_path ) // '/reference_density', eos_coefs( 1 ) )
              call get_option( trim( eos_option_path ) // '/C0', eos_coefs( 2 ), default = 0. )
              call get_option( trim( eos_option_path ) // '/alpha', eos_coefs( 3 ), default  = -1.)
              call get_option( trim( eos_option_path ) // '/T0', eos_coefs( 4 ), default = 298. )
              call get_option( trim( eos_option_path ) // '/beta', eos_coefs( 5 ), default = -1.)
              call get_option( trim( eos_option_path ) // '/P0', eos_coefs( 6 ), default = 1e5 )
              call get_option( trim( eos_option_path ) // '/gamma', eos_coefs( 7 ), default = -1.)
              !Now compute formula
              call linear_EOS_formula(rho)

              if (has_boussinesq_aprox .or. eos_coefs( 7 ) < 0) then
                !If boussinesq or not pressure dependency the derivative is zero
                dRhodP = 0.0
              else !Compute pressure derivatives
                ! Back up pressure and density before we start perturbing stuff...
                allocate( pressure_back_up( node_count( pressure ) ), density_back_up( node_count( pressure ) ) )
                pressure_back_up = pressure % val(1,1,:); density_back_up = density % val
                ! redefine p as p+pert and p-pert and then run python state again to get dRho / d P...
                perturbation_pressure = 1.e-5
                pressure % val(1,1,:) = pressure_back_up + perturbation_pressure
                call linear_EOS_formula(RhoPlus)
                pressure % val(1,1,:) = pressure_back_up - perturbation_pressure
                call linear_EOS_formula(RhoMinus)
                ! derivative
                dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure
                ! Restore pressure and density values in state
                pressure % val(1,1,:) = pressure_back_up; density % val = density_back_up
                deallocate(pressure_back_up, density_back_up)

              end if
              deallocate( eos_coefs )
        elseif( trim( eos_option_path ) == trim( option_path_python ) ) then

            density => extract_scalar_field( state( iphase ), 'Density', stat )
            rho = density % val
            !Obtain density from the python code
            call multi_compute_python_field(state, iphase, trim( option_path_python ), Rho)
            ! Back up pressure before we start perturbing stuff...
            density % val = Rho
            if (has_boussinesq_aprox ) then
              !If boussinesq or not pressure dependency the derivative is zero
              dRhodP = 0.0
            else
              ! Calculating dRho / dP
              allocate( pressure_back_up( node_count( pressure ) ))
              pressure_back_up = pressure % val(1,1,:)
              ! redefine p as p+pert and p-pert and then run python state again to get dRho / d P...
              perturbation_pressure = 1.e-5
              pressure % val(1,1,:) = pressure_back_up + perturbation_pressure; RhoPlus = Rho
              call multi_compute_python_field(state, iphase, trim( option_path_python ), RhoPlus)
              pressure % val(1,1,:) = pressure_back_up - perturbation_pressure; RhoMinus = Rho
              call multi_compute_python_field(state, iphase, trim( option_path_python ), RhoMinus)
              ! derivative
              dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure

              ! Restore pressure and density values in state
              pressure % val(1,1,:) = pressure_back_up
              deallocate( pressure_back_up )
            end if
        else
            FLAbort( 'No option given for choice of EOS' )
        end if Conditional_EOS_Option

        deallocate( perturbation_pressure, RhoPlus, RhoMinus )

        !No need to update halos as all the operations are local, and all the
        !input fields have to be updated already
        ! if (IsParallel()) call halo_update(density)

    contains
      subroutine linear_EOS_formula(rho_internal)
        real, dimension( : ), intent( inout ) :: rho_internal
        rho_internal = 1.0
        !Add the concentration contribution
        if (eos_coefs( 3 ) > 0 ) then
          if (.not. has_concentration) print *, "ERROR: EOS defined to use concentration but concentration field is not defined"
          rho_internal =  rho_internal + eos_coefs( 3 ) * (Concentration % val - eos_coefs( 2 ) )
        end if
        !add the temperature contribution
        if (eos_coefs( 5 ) > 0) then
          if (.not. has_temperature) print *, "ERROR: EOS defined to use temperature but temperature field is not defined"
          rho_internal = rho_internal - eos_coefs( 5 ) * (temperature % val - eos_coefs( 4 ))
        end if
        !add pressure contribution
        if (eos_coefs( 7 ) > 0. ) rho_internal =  rho_internal + eos_coefs( 7 ) * (pressure%val(1,1,:) - eos_coefs( 6 ) )
        !Now add values from scalar fields such as passive tracers
        buffer = "/material_phase["// int2str( iphase -1 ) //"]/phase_properties/Density/compressible/Linear_eos/"
        nfields = option_count(trim(buffer)//"scalar_field")
        do ifield = 1, nfields
          call get_option(trim(buffer)//"scalar_field["// int2str( ifield - 1 )//"]/name",option_name)
          pnt_sfield => extract_scalar_field( state( iphase ), trim(option_name), stat )
          if (stat /=0) then
            FLAbort("ERROR: Field defined for Boussinesq EOS does not exists. Field name: "// trim(option_name))
          end if!We now reuse coefficients
          call get_option( trim(buffer)//"scalar_field["// int2str( ifield - 1 )//"]/R0", eos_coefs( 2 ))!Reference
          call get_option( trim(buffer)//"scalar_field["// int2str( ifield - 1 )//"]/Coef", eos_coefs( 3 ))!Coefficient
          !Now include into EOS
          rho_internal =  rho_internal + eos_coefs( 3 ) * (pnt_sfield % val - eos_coefs( 2 ) )
        end do
        !Finally multiply by the reference density
        rho_internal = rho_internal * eos_coefs( 1 )
        !Ensure that the density does not vary more than 10%, in theory it should never pass 5%
        rho_internal = min(max(rho_internal, eos_coefs( 1 )/1.1), eos_coefs( 1 )*1.1)

      end subroutine
    end subroutine Calculate_Rho_dRhodP

    !--------------------------------------------
    !> @author Geraldine Regnier
    !> @brief In this subroutine we calculate and update the density of the
    !> porous media based on the compressibility given in Diamond (ele-wise). We also
    !> calculate the drho/dp term for the porous media which then goes into
    !> the DERIV term in the continuity equation (cv-wise).
    !---------------------------------------------
    subroutine Calculate_porous_Rho_dRhoP(state,packed_state,Mdims, cv_ndgln, drhodp_porous )

        implicit none
        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent(in) :: Mdims
        integer, dimension( : ), intent(in) :: cv_ndgln
        real, dimension( : ), intent(inout) :: drhodp_porous
        !Local variables
        type( tensor_field ), pointer :: pressure
        type( scalar_field ), pointer :: density_porous, density_porous_initial, density_porous_old, porous_coef
        real, dimension( : ), allocatable :: RhoPlus, RhoMinus
        real, dimension(Mdims%totele) :: cv_counter
        integer :: stat
        real, dimension( : ), allocatable :: rho_porous, rho_porous_old
        real, parameter :: perturb = 1e-5
        real :: ref_pressure

        integer :: ele, cv_inod, iloc, i, k, multiplier, multiplier2

        !!$ Den = Den_surface*exp(C0 * ( P_res-P_surf) )
        pressure => extract_tensor_field( packed_state, 'PackedCVPressure', stat )
        if (stat/=0) pressure => extract_tensor_field( packed_state, "PackedFEPressure", stat )
        density_porous => extract_scalar_field( state(1), "porous_density" )
        density_porous_initial  => extract_scalar_field( state(1), "porous_density_initial" )
        density_porous_old => extract_scalar_field(state(1), "porous_density_old")
        porous_coef => extract_scalar_field(state(1), "porous_compressibility")
        call get_option("/numerical_methods/Surface_pressure", ref_pressure, default = 1e5)


        allocate( RhoPlus( Mdims%cv_nonods ) ) ; RhoPlus = 0.
        allocate( RhoMinus( Mdims%cv_nonods ) ) ; RhoMinus = 0.
        allocate( rho_porous( size(density_porous%val) ) )
        allocate( rho_porous_old(size(density_porous%val)  ) )
        rho_porous=0.
        cv_counter = 0
        rho_porous_old = density_porous%val

        !To loop over the porous_compressibility coefficient
        multiplier = 1; if (size(porous_coef%val) == 1)  multiplier = 0
        multiplier2 = 1; if (size(density_porous%val) == 1)  multiplier2 = 0

        do ele = 1, Mdims%totele
          i = ele * multiplier + (1 - multiplier)
          k = ele * multiplier2 + (1 - multiplier2)
            do iloc = 1,Mdims%cv_nloc
                cv_inod = cv_ndgln((ele-1)*Mdims%cv_nloc+iloc)
                cv_counter( ele ) = cv_counter( ele ) + 1.0
                rho_porous(k) = rho_porous(k)+ density_porous_initial%val(k )&
                  *exp(porous_coef%val( i ) * (pressure % val(1,1,cv_inod) -ref_pressure))
                RhoPlus(cv_inod) =  density_porous_initial%val(k )*exp(porous_coef%val( i ) * (pressure % val(1,1,cv_inod)*(1.+ perturb) - ref_pressure))
                RhoMinus(cv_inod) = density_porous_initial%val(k )*exp(porous_coef%val( i ) * (pressure % val(1,1,cv_inod)*(1.- perturb) - ref_pressure))
                !Obtain derivative
                drhodp_porous(cv_inod) = 0.5 * ( RhoPlus(cv_inod) - RhoMinus(cv_inod))  / pressure % val(1,1,cv_inod)*perturb
            end do
        end do

        !make average
        rho_porous = rho_porous/cv_counter
        !Avoid instabilities
        rho_porous = min(rho_porous,density_porous_initial%val*10. )
        rho_porous = max(rho_porous,density_porous_initial%val/10. )
        if (maxval(drhodp_porous)>100.0) drhodp_porous = 1.e-5
        density_porous%val = rho_porous
        density_porous_old%val = rho_porous_old
        deallocate( RhoPlus, RhoMinus, rho_porous, rho_porous_old )



    end subroutine Calculate_porous_Rho_dRhoP

    !>@brief: Define de density as a polynomial
    subroutine Density_Polynomial( eos_coefs, pressure, temperature, &
        Density_Field )
        implicit none
        real, dimension( : ), intent( in ) :: eos_coefs, pressure, temperature
        real, dimension( : ), intent( inout ) :: Density_Field

        Density_Field = eos_coefs( 1 ) + eos_coefs( 2 ) * pressure + eos_coefs( 3 ) * temperature + &
            eos_coefs( 4 ) * pressure * temperature + eos_coefs( 5 ) * pressure **2 + &
            eos_coefs( 6 ) * temperature **2 + eos_coefs( 7 ) * ( pressure ** 2 ) * temperature + &
            eos_coefs( 8 ) * ( temperature ** 2 ) * pressure + &
            eos_coefs( 9 ) * ( temperature ** 2 ) * ( pressure ** 2 )

        return
    end subroutine Density_Polynomial

    !>@brief: Read from diamond and decide which sort of density representation do we have
    subroutine Assign_Equation_of_State( eos_option_path_out )
        implicit none
        character( len = option_path_len ), intent( inout ) :: eos_option_path_out

        Conditional_for_Compressibility: if( have_option( trim( eos_option_path_out ) // '/compressible' ) ) then
            eos_option_path_out = trim( eos_option_path_out ) // '/compressible'

            Conditional_for_Compressibility_Option: if( have_option( trim( eos_option_path_out ) // '/stiffened_gas' ) ) then
                eos_option_path_out = trim( eos_option_path_out ) // '/stiffened_gas'

            elseif( have_option( trim( eos_option_path_out ) // '/linear_in_pressure' ) ) then
                eos_option_path_out = trim( eos_option_path_out ) // '/linear_in_pressure'

                if( have_option( trim( eos_option_path_out ) // '/include_internal_energy' ) ) &
                    eos_option_path_out = trim( eos_option_path_out ) // '/include_internal_energy'

            elseif( have_option( trim( eos_option_path_out ) // '/exponential_in_pressure' ) ) then
                eos_option_path_out = trim( eos_option_path_out ) // '/exponential_in_pressure'

              elseif( have_option( trim( eos_option_path_out ) // '/Linear_eos' ) ) then
                  eos_option_path_out = trim( eos_option_path_out ) // '/Linear_eos'

            elseif( have_option( trim( eos_option_path_out ) // '/Temperature_Pressure_correlation' ) ) then
                eos_option_path_out = trim( eos_option_path_out ) // '/Temperature_Pressure_correlation'

            else
                FLAbort( 'No option given for choice of EOS - compressible fluid' )

            end if Conditional_for_Compressibility_Option

        elseif( have_option( trim( eos_option_path_out ) // '/incompressible' ) )then
            eos_option_path_out = trim( eos_option_path_out ) // '/incompressible'

        elseif( have_option( trim( eos_option_path_out ) // '/python_state' ) ) then
            eos_option_path_out = trim( eos_option_path_out ) // '/python_state'

        else

            FLAbort( 'No option given for choice of EOS' )

        end if Conditional_for_Compressibility

        return
    end subroutine Assign_Equation_of_State

    !>@brief: Here we compute the absorption for porous media
    !> This is the sigma term defined in the papers and contains The permeability, the relative permeability, the viscosity and the saturation
    subroutine Calculate_PorousMedia_AbsorptionTerms( nphase, state, packed_state, PorousMedia_absorp, Mdims, CV_funs, CV_GIdims, Mspars, ndgln, &
                                                      upwnd, suf_sig_diagten_bc )
       implicit none
       integer, intent(in) :: nphase
       type( state_type ), dimension( : ), intent( inout ) :: state
       type( state_type ), intent( inout ) :: packed_state
       type (multi_field) :: PorousMedia_absorp
       type( multi_dimensions ), intent( in ) :: Mdims
       type(multi_shape_funs), intent(inout) :: CV_funs
       type( multi_gi_dimensions ), intent( in )  :: CV_GIdims
       type (multi_sparsities), intent( in ) :: Mspars
       type(multi_ndgln), intent(in) :: ndgln
       type (porous_adv_coefs), intent(inout) :: upwnd
       real, dimension( :, : ), intent( inout ) :: suf_sig_diagten_bc
       !Local variables
       type( tensor_field ), pointer :: perm, state_viscosity
       real, dimension(:,:), allocatable :: viscosities
       integer :: i, j, ele, n_in_pres, iphase
       real :: Angle, bad_element_perm_mult, Max_aspect_ratio, height ! the height of an isosceles triangle for the top angle to be equal to the trigger angle
       real, dimension(Mdims%ndim,Mdims%ndim) :: trans_matrix, rot_trans_matrix ! for bad_element permeability transformation matrix
       real, parameter :: pi = acos(0.0) * 2.0 ! Define pi
       character( len = option_path_len ) :: option_path_python

       perm => extract_tensor_field( packed_state, "Permeability" )
       !Define n_in_pres based on the local version of nphase
       n_in_pres = nphase/Mdims%npres

       !Obtain inverse of permeability and store it
       !SPRINT_TO_DO THIS COULD BE DONE FASTER IF WE KNOW IF IT HAS OFF DIAGONALS OR NOT
       do i = 1, size(perm%val,3)
           upwnd%inv_permeability( :, :, i)=inverse(perm%val( :, :, i))
       end do

       allocate(viscosities(Mdims%nphase, Mdims%cv_nonods))
       DO IPHASE = 1, Mdims%nphase!Get viscosity for all the phases
        option_path_python = "/material_phase["// int2str( iphase - 1 )//"]/phase_properties/Viscosity/tensor_field"//&
        "::Viscosity/diagnostic/algorithm::tensor_python_diagnostic"
          if (have_option(trim(option_path_python)))then
            state_viscosity => extract_tensor_field( state( iphase ), 'Viscosity' )
            call multi_compute_python_field(state, iphase, trim( option_path_python ), tfield = state_viscosity)
            !Copy into state
            do i = 1, Mdims%cv_nonods
              viscosities(iphase,i) = state_viscosity%val(1,1,i)
            end do
          else
            call set_viscosity(nphase, Mdims, state, viscosities)
          end if
        end do
       call Calculate_PorousMedia_adv_terms( nphase, state, packed_state, PorousMedia_absorp, Mdims, ndgln, &
              upwnd, viscosities)

       ! calculate SUF_SIG_DIAGTEN_BC this is \sigma_in^{-1} \sigma_out
       ! \sigma_in and \sigma_out have the same anisotropy so SUF_SIG_DIAGTEN_BC
       ! is diagonal
       call calculate_SUF_SIG_DIAGTEN_BC( nphase, packed_state, suf_sig_diagten_bc, Mdims, CV_funs, CV_GIdims, &
                              Mspars, ndgln, PorousMedia_absorp, state, upwnd%inv_permeability, viscosities)

       deallocate(viscosities)
       contains
          !>@brief: Computes the absorption and its derivatives against the saturation
           subroutine Calculate_PorousMedia_adv_terms( nphase, state, packed_state, PorousMedia_absorp, Mdims, ndgln, &
               upwnd, viscosities )

               implicit none
               integer, intent(in) :: nphase
               type( state_type ), dimension( : ), intent( in ) :: state
               type( state_type ), intent( inout ) :: packed_state
               type (multi_field), intent( inout ) :: PorousMedia_absorp
               type( multi_dimensions ), intent( in ) :: Mdims
               type( multi_ndgln ), intent( in ) :: ndgln
               type (porous_adv_coefs), intent(inout) :: upwnd
               real, dimension(:,:), intent(in) :: viscosities
               !!$ Local variables:
               integer :: ele, imat, icv, iphase, cv_iloc, idim, jdim, ipres, loc, n_in_pres, &
                global_phase, compact_phase
               real :: Mobility, pert
               real, dimension(:), allocatable :: Max_sat
               real, dimension( :, : ), allocatable :: satura2
               type (multi_field) :: PorousMedia_absorp2
               !Working pointers
               real, dimension(:,:), pointer :: Satura, OldSatura, Immobile_fraction
               type( tensor_field ), pointer :: perm
               type( scalar_field ), pointer :: Spipe

               !Define n_in_pres based on the local version of nphase
               n_in_pres = nphase/Mdims%npres

               !Initialize variables
               upwnd%adv_coef_grad=0.0;upwnd%inv_adv_coef=0.0
               !Get from packed_state
               call get_var_from_packed_state(packed_state,PhaseVolumeFraction = Satura,&
                   OldPhaseVolumeFraction = OldSatura, Immobile_fraction = Immobile_fraction)
               perm=>extract_tensor_field(packed_state,"Permeability")

               allocate( satura2( nphase, size(SATURA,2) ) );satura2 = 0.

               CALL calculate_absorption2( nphase, packed_state, PorousMedia_absorp, Mdims, ndgln, SATURA(1:Mdims%n_in_pres,:), &
                      viscosities)

               !Introduce perturbation, positive for the increasing and negative for decreasing phase
               !Make sure that the perturbation is between bounds
               PERT = 0.0001; allocate(Max_sat(Mdims%nphase))
               do ele = 1, Mdims%totele
                   do cv_iloc = 1, Mdims%cv_nloc
                       icv = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                       Max_sat(:) = 1. - sum(Immobile_fraction(:, ele)) + Immobile_fraction(:, ele)
                       DO IPHASE = 1, n_in_pres!Not for wells
                         SATURA2(IPHASE, icv) = SATURA(IPHASE, icv) + sign(PERT, satura(iphase, icv)-OldSatura(iphase, icv))
                         !If out of bounds then we perturbate in the opposite direction
                         if (satura2(IPHASE, icv) > Max_sat(iphase) .or. &
                             satura2(IPHASE, icv) < Immobile_fraction(iphase, ele)) then
                             SATURA2(IPHASE, icv) = SATURA2(IPHASE, icv) - 2. * sign(PERT, satura(iphase, icv)-OldSatura(iphase, icv))
                         end if
                       end do
                   end do
               end do


               call allocate_multi_field( Mdims, PorousMedia_absorp2, size(PorousMedia_absorp%val,4), field_name="PorousMedia_AbsorptionTerm")
               CALL calculate_absorption2( nphase, packed_state, PorousMedia_absorp2, Mdims, ndgln, SATURA2, viscosities)

               do ipres = 2, Mdims%npres
                   Spipe => extract_scalar_field( state(1), "Sigma" )
                   do iphase = 1, n_in_pres
                     ! set \sigma for the pipes here
                     call assign_val(PorousMedia_absorp%val(1, 1, iphase + (ipres - 1)*Mdims%n_in_pres, :),Spipe%val)
                   end do
               end do

               !Temporary pointer, maybe we should unify memories
               nullify(upwnd%adv_coef)
               upwnd%adv_coef => PorousMedia_absorp%val

               DO ELE = 1, Mdims%totele
                 DO CV_ILOC = 1, Mdims%cv_nloc
                   IMAT = ndgln%mat( ( ELE - 1 ) * Mdims%mat_nloc + CV_ILOC )
                   ICV = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
                   do ipres = 1, Mdims%npres
                     DO IPHASE = 1, n_in_pres
                       global_phase = iphase + (ipres - 1)*Mdims%n_in_pres
                       compact_phase = iphase + (ipres - 1)*n_in_pres
                       if ( global_phase <= Mdims%n_in_pres ) then
                         ! This is the gradient
                         ! Assume d\sigma / dS = 0.0 for the pipes for now
                         upwnd%adv_coef_grad(1, 1, global_phase, IMAT) = (PorousMedia_absorp2%val( 1,1, global_phase ,IMAT) -&
                         PorousMedia_absorp%val( 1,1, global_phase ,IMAT)) / ( SATURA2(compact_phase, ICV ) - SATURA(compact_phase, ICV))
                       end if
                       !Obtaining the inverse the "old way" since if you obtain it directly, some problems appear
                       upwnd%inv_adv_coef(1, 1, global_phase, IMAT) = 1./upwnd%adv_coef(1, 1, global_phase, IMAT)!sprint_to_do maybe we dont need to store the inverse anymore
                     END DO
                   end do
                 END DO
               END DO

               deallocate( satura2, Max_sat)
               call deallocate_multi_field(PorousMedia_absorp2, .true.)

           end subroutine Calculate_PorousMedia_adv_terms

           !>@brief: Computes the absorption and its derivatives against the saturation on the boundary
           subroutine calculate_SUF_SIG_DIAGTEN_BC( nphase, packed_state, suf_sig_diagten_bc, Mdims, CV_funs, CV_GIdims, &
               Mspars, ndgln, PorousMedia_absorp, state, inv_perm, viscosities)
               implicit none
               integer, intent(in) :: nphase
               type( state_type ), intent( inout ) :: packed_state
               type(multi_dimensions), intent(in) :: Mdims
               type(multi_GI_dimensions), intent(in) :: CV_GIdims
               type(multi_shape_funs), intent(inout) :: CV_funs
               type (multi_sparsities), intent(in) :: Mspars
               type(multi_ndgln), intent(in) :: ndgln
               type (multi_field), intent( inout ) :: PorousMedia_absorp
               type(state_type), dimension( : ), intent(in) :: state
               real, dimension( Mdims%stotel * Mdims%cv_snloc * Mdims%nphase, Mdims%ndim ), intent( inout ) :: suf_sig_diagten_bc
               real, dimension(:, :, :), target, intent(in):: inv_perm
               real, dimension(:,:), intent(in) :: viscosities
               ! local variables
               type(tensor_field), pointer :: RockFluidProp
               real, dimension(:), pointer :: Immobile_fraction, Corey_exponent, Endpoint_relperm
               integer :: iphase, ele, sele, cv_siloc, cv_snodi, cv_snodi_ipha, iface, s, e, &
                   ele2, sele2, cv_iloc, idim, jdim, i, mat_nod, cv_nodi
               real :: sigma_out!, mat, mat_inv
               ! real, dimension( Mdims%ndim, Mdims%ndim ) :: mat_ones, mat, mat_inv
               integer, dimension( CV_GIdims%nface, Mdims%totele) :: face_ele
               integer, dimension( Mdims%cv_snloc ) :: cv_sloc2loc
               integer, dimension( :, :, : ),  allocatable :: wic_u_bc, wic_vol_bc
               integer, parameter :: WIC_BC_DIRICHLET = 1
               type(tensor_field), pointer :: velocity, volfrac, perm
               type(tensor_field) :: velocity_BCs, volfrac_BCs
               integer :: one_or_zero, visc_node, n_in_pres
               !Define n_in_pres based on the local version of nphase
               n_in_pres = nphase/Mdims%npres
               !Prepapre index for viscosity
               one_or_zero = (size(viscosities,2)==Mdims%cv_nonods)
               !Get from packed_state
               volfrac=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
               velocity=>extract_tensor_field(packed_state,"PackedVelocity")
               perm=>extract_tensor_field(packed_state,"Permeability")
               RockFluidProp=>extract_tensor_field(packed_state,"PackedRockFluidProp")

               allocate(wic_u_bc(velocity%dim(1),velocity%dim(2),&
                   surface_element_count(velocity)))
               allocate(wic_vol_bc(volfrac%dim(1),volfrac%dim(2),&
                   surface_element_count(volfrac)))
               call get_entire_boundary_condition(velocity,&
                   ['weakdirichlet'],velocity_BCs,WIC_U_BC)
               call get_entire_boundary_condition(volfrac,&
                   ['weakdirichlet'],volfrac_BCs,WIC_vol_BC)

               suf_sig_diagten_bc = 1.
               face_ele = 0
               call calc_face_ele( face_ele, Mdims%totele, Mdims%stotel, CV_GIdims%nface, &
                   Mspars%ELE%fin, Mspars%ELE%col, Mdims%cv_nloc, Mdims%cv_snloc, Mdims%cv_nonods, ndgln%cv, ndgln%suf_cv, &
                   CV_funs%cv_sloclist, Mdims%x_nloc, ndgln%x )


                   do ele = 1, Mdims%totele
                       !Get properties from packed state
                       Immobile_fraction => RockFluidProp%val(1, :, ELE)
                       Endpoint_relperm => RockFluidProp%val(2, :, ELE)
                       Corey_exponent => RockFluidProp%val(3, :, ELE)
                       do iface = 1, CV_GIdims%nface
                           ele2  = face_ele( iface, ele )
                           sele2 = max( 0, -ele2 )
                           sele  = sele2
                           if ( sele > 0 ) then
                               do iphase = 1, n_in_pres
                                   s = ( iphase - 1 ) * Mdims%ndim + 1
                                   e = iphase * Mdims%ndim
                                   if ( wic_u_bc(1,iphase,sele) /= WIC_BC_DIRICHLET .and. &
                                       wic_vol_bc(1,iphase,sele) == WIC_BC_DIRICHLET ) then
                                       cv_sloc2loc( : ) = CV_funs%cv_sloclist( iface, : )
                                       do cv_siloc = 1, Mdims%cv_snloc
                                           cv_iloc = cv_sloc2loc( cv_siloc )
                                           cv_snodi = ( sele - 1 ) * Mdims%cv_snloc + cv_siloc
                                           cv_nodi = ndgln%suf_cv(cv_snodi)
                                           visc_node = (cv_nodi-1)*one_or_zero + 1
                                           cv_snodi_ipha = cv_snodi + ( iphase - 1 ) * Mdims%stotel * Mdims%cv_snloc
                                           mat_nod = ndgln%mat( (ele-1)*Mdims%cv_nloc + cv_iloc  )
                                           call get_material_absorption(Mdims%n_in_pres, iphase, sigma_out,&
                                               ! this is the boundary condition (we pass as old sat the same BC)
                                               volfrac_BCs%val(1,:,cv_snodi), viscosities(:,visc_node),&
                                               Immobile_fraction, Corey_exponent, Endpoint_relperm)
                                           ! Adjust suf_sig_diagten_bc based on the internal absorption
                                            suf_sig_diagten_bc( cv_snodi_ipha, 1 : Mdims%ndim ) =  &
                                                  (sigma_out  +  PorousMedia_absorp%val(1,1,iphase, mat_nod)**2./sigma_out )&
                                                  /(sigma_out  +  PorousMedia_absorp%val(1,1,iphase, mat_nod))
                                       end do
                                   end if
                               end do
                           end if
                       end do
                   end do
                   call deallocate(velocity_BCs)
               call deallocate(volfrac_BCs)
               deallocate(wic_u_bc, wic_vol_bc)
               return
           end subroutine calculate_SUF_SIG_DIAGTEN_BC

        !>@brief: For porous media, just sets the viscosity as a scalar
        subroutine set_viscosity(nphase, Mdims, state, visc_phases)
            implicit none
            integer, intent(in) :: nphase
            type( multi_dimensions ), intent( in ) :: Mdims
            type( state_type ), dimension( : ), intent( in ) :: state
            real, dimension(:,:), intent(inout) :: visc_phases
            !Local variables
            integer :: iphase, ipres, global_phase, compact_phase
            real :: mobility
            type(tensor_field), pointer :: viscosity_ph


            !SPRINT_TO_DO what happens here if we have components???
            do ipres = 1, Mdims%npres
              DO IPHASE = 1, nphase/Mdims%npres!Get viscosity for all the phases
                  global_phase = iphase + (ipres - 1)*Mdims%n_in_pres
                  compact_phase = iphase + (ipres - 1)*n_in_pres
                  viscosity_ph => extract_tensor_field( state( global_phase ), 'Viscosity' )
                  visc_phases(compact_phase,:) = viscosity_ph%val( 1, 1, 1 )!So far we only consider scalar viscosity
              end do
            end do
        end subroutine set_viscosity


    end subroutine Calculate_PorousMedia_AbsorptionTerms




    !>@brief: Subroutine where the absorption for the porous media is actually computed
    SUBROUTINE calculate_absorption2( nphase, packed_state, PorousMedia_absorp, Mdims, ndgln, SATURA, &
        viscosities, inv_PorousMedia_absorp)
        ! Calculate absorption for momentum eqns
        implicit none
        integer, intent(in) :: nphase
        type( state_type ), intent( inout ) :: packed_state
        type (multi_field) :: PorousMedia_absorp
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_ndgln), intent(in) :: ndgln
        REAL, DIMENSION( :, : ), intent( in ) :: SATURA
        real, intent(in), dimension(:,:) :: viscosities
        real, dimension(:,:,:), INTENT(INOUT), optional :: inv_PorousMedia_absorp
        ! Local variable
        type (tensor_field), pointer :: RockFluidProp
        real, dimension(:), pointer :: Immobile_fraction, Corey_exponent, Endpoint_relperm
        REAL, PARAMETER :: TOLER = 1.E-10
        INTEGER :: ELE, CV_ILOC, CV_NOD, CV_PHA_NOD, MAT_NOD, JPHA_JDIM, &
            IPHA_IDIM, IDIM, JDIM, IPHASE, id_reg, n_in_pres
        integer :: one_or_zero, visc_node
        !Prepapre index for viscosity
        one_or_zero = (size(viscosities,2)==Mdims%cv_nonods)

        !Define n_in_pres based on the local version of nphase
        n_in_pres = nphase/Mdims%npres

        RockFluidProp=>extract_tensor_field(packed_state,"PackedRockFluidProp")
        ewrite(3,*) 'In calculate_absorption2'

        DO ELE = 1, Mdims%totele
            !Get properties from packed state
            Immobile_fraction => RockFluidProp%val(1, :, ELE)
            Endpoint_relperm => RockFluidProp%val(2, :, ELE)
            Corey_exponent => RockFluidProp%val(3, :, ELE)
            DO CV_ILOC = 1, Mdims%cv_nloc
                MAT_NOD = ndgln%mat(( ELE - 1 ) * Mdims%mat_nloc + CV_ILOC)
                CV_NOD = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + CV_ILOC )
                visc_node = (CV_NOD-1)*one_or_zero + 1
                DO IPHASE = 1, n_in_pres
                    CV_PHA_NOD = CV_NOD + ( IPHASE - 1 ) * Mdims%cv_nonods
                    call get_material_absorption(Mdims%n_in_pres, iphase, PorousMedia_absorp%val(1, 1, iphase, mat_nod),&
                        SATURA(:, CV_NOD), viscosities(:,visc_node),Immobile_fraction, Corey_exponent,&
                         Endpoint_relperm)
                END DO
            END DO
        END DO
        ewrite(3,*) 'Leaving calculate_absorption2'
        RETURN
    END SUBROUTINE calculate_absorption2


    !>@brief:Calculates the relative permeability for 1, 2 (Brooks-corey) or 3 (stone's model) phases
    subroutine get_material_absorption(nphase, iphase, material_absorption, sat, visc, Immobile_fraction, &
            Corey_exponent, Endpoint_relperm )
        implicit none
        real, intent(inout) :: material_absorption
        real, dimension(:), intent(in) :: sat, visc, Immobile_fraction, Corey_exponent, Endpoint_relperm
        integer, intent(in) :: iphase, nphase
        !local variables
        real :: Kr
        !Local parameters
        real, parameter :: eps = 1d-5!eps is another epsilon value, for less restrictive things

        call get_relperm(nphase, iphase, sat, Immobile_fraction, Corey_exponent, Endpoint_relperm, Kr )

        material_absorption = (visc(iphase) * max(eps, sat(iphase))) / KR !The value 1d-5 is only used if the boundaries have values of saturation of zero.

        !Not needed anymore?
        ! if (present(inv_mat_absorp)) &!This part is to ensure that the flow is stopped
        ! inv_mat_absorp = (max(0.0,Kr))/(VISC(iphase) * max(eps,sat(iphase)))
        !Otherwise, the saturation should never be zero, since immobile fraction is always bigger than zero.
      end subroutine get_material_absorption

    !>@brief:Calculates the relative permeability for 1, 2 (Brooks-corey) or 3 (stone's model) phases
    subroutine get_relperm(nphase, iphase, sat, Immobile_fraction, Corey_exponent, Endpoint_relperm, Kr)
        implicit none
        real, INTENT(INOUT) :: Kr
        real, dimension(:), intent(in) :: sat, Immobile_fraction, Corey_exponent, Endpoint_relperm
        integer, intent(in) :: iphase, nphase
        !Local parameters
        real, parameter :: eps = 1d-5!eps is another epsilon value, for less restrictive things
        real, parameter :: epsilon = 1d-8!This value should in theory never be used, the real lower limit

        select case (nphase)
            case (1)
                Kr = 1.0
            case (3)
                call relperm_stone(Kr)
            case default
                call relperm_corey_epsilon(Kr)
        end select

      contains
        !>@brief: Brooks corey model of relperm.
        !>This subroutine add a small quantity to the corey function to avoid getting a relperm=0 that may give problems
        !>when dividing it to obtain the sigma.
        SUBROUTINE relperm_corey_epsilon( Kr )
          IMPLICIT NONE
          REAL, intent( inout ) :: Kr
          ! Local variables...
          REAL :: aux
          aux = 1.0 - sum(Immobile_fraction)
          KR = Endpoint_relperm(iphase)*( max( sat(iphase) - Immobile_fraction(iphase), sat(iphase)*eps+eps) / aux ) ** Corey_exponent(iphase)
          !Make sure that the relperm is between bounds
          KR = min(max(epsilon, KR),Endpoint_relperm(iphase))!Lower value just to make sure we do not divide by zero.
        END SUBROUTINE relperm_corey_epsilon

        !>@brief:This subroutine calculates the relative permeability for three phases
        !>First phase has to be water, second oil and the third gas
        !>We use Stone's model II adapted, and for the two phases we use the Corey model
        !>Model explained in: Aziz, K. And Settari, T.:“Petroleum Reservoir Simulation” Applied Science Publishers, London, 30-38, 1979.
        subroutine relperm_stone(Kr)
          implicit none
          real, intent(inout) :: Kr
          !Local variables
          real, dimension(3) :: Norm_sat, relperm
          real :: Krow, Krog

          !Prepare data
          !We consider two models for two phase flow, water-oil and oil-gas
          if (iphase /= 3) then
            Norm_sat(1) = ( sat(1) - Immobile_fraction(1)) /( 1. - Immobile_fraction(1) - Immobile_fraction(2))!Water
            relperm(1) = Endpoint_relperm(1)* Norm_sat(1) ** Corey_exponent(1)!Water, Krw
          end if
          if (iphase /= 1) then
            Norm_sat(3) = ( sat(3) - Immobile_fraction(3)) /(1. - Immobile_fraction(2) - Immobile_fraction(1))!Gas
            !For phase 1 and 3 (water and gas respectively) we can use the Brooks Corey model
            relperm(3) = Endpoint_relperm(3)* Norm_sat(3) ** Corey_exponent(3)!Gas, Krg
          end if
          !Oil relperm is obtained as a combination
          if (iphase == 2 ) then
            Krow = Endpoint_relperm(2)* (1.0 - Norm_sat(1)) ** Corey_exponent(2)!Oil, Krow
            Krog = Endpoint_relperm(2)* (1.0 - Norm_sat(3)) ** Corey_exponent(2)!Oil, Krog
            !For the second phase, oil, we need to recalculate the real value(Stone model 2)
            relperm(2) = Endpoint_relperm(2)*( (Krow/Endpoint_relperm(2) + relperm(1))*&
            (Krog/Endpoint_relperm(2) + relperm(3)) - (relperm(1) + relperm(3)) )
          end if
          !Make sure that the relperm is between bounds
          Kr = min(max(epsilon, relperm(iphase)),Endpoint_relperm(iphase))!Lower value just to make sure we do not divide by zero.
        end subroutine relperm_stone

    end subroutine get_relperm

    !>@brief: In this subroutine the capilalry pressure is computed based on the saturation and the formula used
    SUBROUTINE calculate_capillary_pressure( packed_state, &
        NDGLN, totele, cv_nloc, CV_funs)

            ! CAPIL_PRES_OPT is the capillary pressure option for deciding what form it might take.
            ! CAPIL_PRES_COEF( NCAPIL_PRES_COEF, NPHASE, NPHASE ) are the coefficients
            ! Capillary pressure coefs have the dims CAPIL_PRES_COEF( NCAPIL_PRES_COEF, NPHASE,NPHASE )
            ! used to calculate the capillary pressure.

            IMPLICIT NONE
            type(state_type), intent(inout) :: packed_state
            type(multi_shape_funs) :: CV_funs ! shape function of a reference control volume
            integer, intent(in) :: totele, cv_nloc
            type(multi_ndgln), intent(in) :: ndgln
            ! Local Variables
            type(multi_dev_shape_funs) :: DevFuns ! derivative of the shape functions of the reference control volumes
            INTEGER :: IPHASE, JPHASE, nphase, ele, cv_iloc, cv_nod
            logical, save :: Cap_Brooks = .true., Cap_Power = .false.
            logical, save :: first_time = .true.
            !Working pointers
            real, dimension(:,:), pointer :: Satura, CapPressure, Immobile_fraction, Cap_entry_pressure, Cap_exponent, Imbibition_term, X_ALL
            real, dimension(:), allocatable :: Cont_correction
            !Get from packed_state
            call get_var_from_packed_state(packed_state,PhaseVolumeFraction = Satura)

            call get_var_from_packed_state(packed_state,CapPressure = CapPressure, &
                Immobile_fraction = Immobile_fraction, Cap_entry_pressure = Cap_entry_pressure, &
                    Cap_exponent = Cap_exponent, Imbibition_term = Imbibition_term, PressureCoordinate = X_ALL)

            nphase =size(Satura,1)
            allocate(Cont_correction(size(satura,2)))

            CapPressure = 0.
            ! Determine which capillary pressure model is to be used for overrelaxation. Use Brooks-Corey unless power_law Pc activated (important to allow overelax even when Pc is off).
            if (first_time) then
                Cap_Power = have_option_for_any_phase("/multiphase_properties/capillary_pressure/type_Power_Law", nphase)
                Cap_Brooks = .not. (Cap_Power)

                first_time = .false.
            end if

            ! here we run multi_dev_shape_funs just to calculate element volumes
            call allocate_multi_dev_shape_funs(CV_funs, DevFuns)

            DO IPHASE = 1, NPHASE

                if ( (Cap_Brooks) .or. (Cap_Power) ) then

                    !Apply Capillary model
                    do jphase = 1, nphase
                        Cont_correction = 0
                        if (jphase /= iphase) then!Don't know how this will work for more than 2 phases
                            do ele = 1, totele

                                call DETNLXR(ele, X_ALL, ndgln%x, CV_funs%cvweight, CV_funs%CVFEN, CV_funs%CVFENLX_ALL, DevFuns)

                                do cv_iloc = 1, cv_nloc
                                    cv_nod = ndgln%cv((ele-1)*cv_nloc + cv_iloc)
                                    CapPressure( jphase, cv_nod ) = CapPressure( jphase, cv_nod ) + &
                                        Get_capPressure(satura(iphase,cv_nod), Cap_entry_pressure(iphase, ele), &
                                        Cap_exponent(iphase, ele),Immobile_fraction(:,ele), &
                                        Imbibition_term(iphase, ele), iphase) * DevFuns%volume   ! volume weighted average of capillary pressure
                                    Cont_correction(cv_nod) = Cont_correction(cv_nod) + DevFuns%volume
                                end do
                            end do
                            !In continuous formulation nodes are visited more than once, hence we need to find a volume weighted average
                            CapPressure(jphase, :) = CapPressure(jphase, :) / Cont_correction(:)
                        end if
                    end do

                end if
            END DO

        deallocate(Cont_correction)
        contains
          !>@brief:This functions returns the capillary pressure for a certain input saturation
            pure real function Get_capPressure(sat, Pe, a, Immobile_fraction, Imbibition_term, iphase)
                Implicit none
                real, intent(in) :: sat, Pe, a, Imbibition_term
                real, dimension(:), intent(in) :: Immobile_fraction
                integer, intent(in) :: iphase
                !Local
                real, parameter :: eps = 1d-3 !Small values requires smaller time steps

                if(Cap_Power) then
                    ! Function is Max_Cap_Pressure * (1-S_norm) ^ a Specify Max_Cap_Pressure in C parameter and exponent a (a>0)
                    Get_capPressure = &
                        Pe * ( 1.0 - ( sat - Immobile_fraction(iphase) )/( 1.0 - sum(Immobile_fraction(:)) ) )**a
                else
                    !A*(Swn^-B) - C; entry pressure = A - C
                    Get_capPressure = &
                        Pe * min((sat - Immobile_fraction(iphase) + eps) / (1.0 - sum(Immobile_fraction(:)) ), 1.0) ** (-a) - Imbibition_term
                endif

            end function Get_capPressure
    END SUBROUTINE calculate_capillary_pressure

    !>@brief:This functions returns the derivative of the capillary pressure with respect to the saturation
    real function Get_DevCapPressure(sat, Pe, a, immobile_fraction, iphase, nphase)
        Implicit none
        integer, intent(in) :: iphase, nphase
        real, intent(in) :: sat, Pe, a
        real, dimension(:), intent(in) :: immobile_fraction
        !Local
        real, parameter :: eps = 1d-3
        real :: aux
        integer :: i
        logical, save :: Cap_Brooks = .true., Cap_Power = .false.
        logical, save :: first_time = .true.

        aux = ( 1.0 - sum(immobile_fraction(:)) )
        ! Determine which capillary pressure model is to be used for overrelaxation. Use Brooks-Corey unless power_law Pc activated (important to allow overelax even when Pc is off).
        if (first_time) then
                Cap_Power = have_option_for_any_phase("/multiphase_properties/capillary_pressure/type_Power_Law", nphase)
                Cap_Brooks = .not. (Cap_Power)
                first_time = .false.
        end if

        if(Cap_Power) then
            Get_DevCapPressure = &
                -a*Pe/(1.0 - sum(Immobile_fraction(:)) )  * ( 1.0 - ( sat - Immobile_fraction(iphase) )/( 1.0 - sum(Immobile_fraction(:)) ) ) **(a-1)
        else
            Get_DevCapPressure = &
                -a * Pe * aux**a * min((sat - immobile_fraction(iphase) + eps), 1.0) ** (-a-1)
        endif

    end function Get_DevCapPressure

    !>@brief: This subroutine computed the gravity effect, i.e. rho * g
    subroutine calculate_u_source_cv(Mdims, state, packed_state, den, u_source_cv)
        type(state_type), dimension(:), intent(in) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent(in) :: Mdims
        real, dimension(:,:), intent(in) :: den
        real, dimension(:,:,:), intent(inout) :: u_source_cv

        type(vector_field), pointer :: gravity_direction
        real, dimension(Mdims%ndim) :: g, ref_density
        logical :: have_gravity, high_order_Ph
        real :: gravity_magnitude
        integer :: idim, iphase, nod, stat, start_phase

        ref_density = 0.
        if (have_option("/physical_parameters/gravity/remove_hydrostatic_contribution")) then
          do iphase = 1, Mdims%nphase
            ref_density(iphase) = retrieve_reference_density(state, packed_state, iphase, 0, Mdims%nphase)
          end do
        end if

        call get_option( "/physical_parameters/gravity/magnitude", gravity_magnitude, stat )
        have_gravity = ( stat == 0 )

        !Initialise RHS
        u_source_cv = 0.
        high_order_Ph = have_option( "/physical_parameters/gravity/hydrostatic_pressure_solver" )
        if( have_gravity .and. .not. high_order_Ph) then
          start_phase = 1
          if (high_order_Ph .and. Mdims%npres > 1) then
            start_phase = Mdims%n_in_pres + 1 !hydrostatic_pressure_solver only for the reservoir domain, not the wells domain
          else if (high_order_Ph) then
            return
          end if
            gravity_direction => extract_vector_field( state( 1 ), 'GravityDirection' )
            u_source_cv = 0.
            do nod = 1, Mdims%cv_nonods
                g = node_val( gravity_direction, nod ) * gravity_magnitude
                do iphase = start_phase, Mdims%nphase
                    do idim = 1, Mdims%ndim
                        u_source_cv( idim, iphase, nod ) = (den( iphase, nod ) - ref_density(iphase)) * g( idim )
                    end do
                end do
            end do
        end if

    end subroutine calculate_u_source_cv

    !>@brief: Here we compute component/solute/thermal diffusion coefficient
    subroutine calculate_diffusivity(state, packed_state, Mdims, ndgln, ScalarAdvectionField_Diffusion, TracerName, divide_by_rho_CP)
      type(state_type), dimension(:), intent(in) :: state
      type( state_type ), intent( inout ) :: packed_state
      type(multi_dimensions), intent(in) :: Mdims
      type(multi_ndgln), intent(in) :: ndgln
      real, dimension(:, :, :, :), intent(inout) :: ScalarAdvectionField_Diffusion
      logical, optional, intent(in) :: divide_by_rho_CP !> If we want to normlise the equation by rho CP we can return the diffusion coefficient divided by rho Cp
      character(len=*), optional, intent(in) :: TracerName !> For PassiveTracer with diffusion we pass down the name of the tracer
      !Local variables
      type(scalar_field), pointer :: component, sfield, solid_concentration
      type(tensor_field), pointer :: diffusivity, tfield, den, saturation
      integer :: icomp, iphase, idim, stat, ele
      integer :: iloc, mat_inod, cv_inod, ele_nod, t_ele_nod
      logical, parameter :: harmonic_average=.false.
      real :: expo

      ScalarAdvectionField_Diffusion = 0.0
      if ( Mdims%ncomp > 1 ) then
        do icomp = 1, Mdims%ncomp
          do iphase = 1, Mdims%nphase
            component => extract_scalar_field( state(Mdims%nphase+icomp), 'ComponentMassFractionPhase' // int2str(iphase) )
            diffusivity => extract_tensor_field( state(Mdims%nphase+icomp), 'ComponentMassFractionPhase' // int2str(iphase) // 'Diffusivity', stat )
            if ( stat == 0 ) then
              do ele = 1, Mdims%totele
                do iloc = 1, Mdims%mat_nloc
                  mat_inod = ndgln%mat( (ele-1)*Mdims%mat_nloc + iloc )
                  cv_inod = ndgln%cv( (ele-1)*Mdims%mat_nloc + iloc )
                  if ( .not.harmonic_average ) then
                    do idim = 1, Mdims%ndim
                      ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase ) = &
                      ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase ) + &
                      node_val( component, cv_inod ) * node_val( diffusivity, idim, idim, mat_inod )
                    end do
                  else
                    do idim = 1, Mdims%ndim
                      if (  node_val( diffusivity, idim, idim, mat_inod ) > 0.0 ) then
                        ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase ) = &
                        ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase ) + &
                        node_val( component, cv_inod ) / node_val( diffusivity, idim, idim, mat_inod )
                      end if
                    end do
                  end if
                end do
              end do
            end if
          end do
        end do
      else
        !Note that for the temperature field this is actually the thermal conductivity (in S.I. watts per meter-kelvin => W/(m·K) ).
        if (is_porous_media) then
          !####DIFFUSIVITY FOR POROUS MEDIA ONLY####
          sfield=>extract_scalar_field(state(1),"Porosity")
          den => extract_tensor_field( packed_state,"PackedDensity" )
          !expo used to switch between boussinesq (density ==1) or normal
          expo = 1.; if (has_boussinesq_aprox) expo = 0.

          if (present(TracerName)) then
            do iphase = 1, Mdims%nphase
              !Check if the field is defined for that phase, if the property is defined but not the field then ignore the property
              diffusivity => extract_tensor_field( state(iphase), trim(TracerName)//'Diffusivity', stat )
              if (stat /= 0) cycle!If no field defined then cycle

              do ele = 1, Mdims%totele
                ele_nod = min(size(sfield%val), ele)
                do iloc = 1, Mdims%mat_nloc
                  mat_inod = ndgln%mat( (ele-1)*Mdims%mat_nloc + iloc )
                  cv_inod = ndgln%cv((ele-1)*Mdims%cv_nloc+iloc)
                  do idim = 1, Mdims%ndim
                    ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase ) =    &
                    ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase ) +    &
                    (sfield%val(ele_nod) *den%val(1, 1, cv_inod)**expo * node_val( diffusivity, idim, idim, mat_inod ))
                  enddo
                end do
              end do
            end do
          else
            ! Calculation of the averaged thermal diffusivity as
            ! lambda = (1-porosity) * lambda_p + SUM_of_phases Saturation * (porosity * lambda_f)
            ! Since lambda_p is defined element-wise and lambda_f CV-wise we perform an average
            ! as it is stored cv-wise
            ! NOTE: for porous media we consider thermal equilibrium and therefore unifiying lambda is a must
            ! Multiplied by the saturation so we use the same paradigm that for the phases,
            !but in the equations it isn't, but here because we iterate over phases and collapse this is required
            ! therefore: lambda = SUM_of_phases saturation * [(1-porosity) * lambda_p + porosity * lambda_f)] for classic
            !weighted average of conductivities Wiener method).
            !Default option is to use a more accurate Hashin and Shtrikman definition:
            !lambda_p+3*lambda_p*(lambda_f-lambda_p)*porosity/(3*lambda_p+(lambda_f-lambda_p)(1-porosity))
            saturation => extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
            do iphase = 1, Mdims%nphase
              diffusivity => extract_tensor_field( state(iphase), 'TemperatureDiffusivity', stat )
              if (stat /= 0) cycle!If no field defined then cycle
              tfield => extract_tensor_field( state(1), 'porous_thermal_conductivity', stat )
              do ele = 1, Mdims%totele
                ele_nod = min(size(sfield%val), ele)
                t_ele_nod = min(size(tfield%val, 3), ele)
                do iloc = 1, Mdims%mat_nloc
                  mat_inod = ndgln%mat( (ele-1)*Mdims%mat_nloc + iloc )
                  cv_inod = ndgln%cv((ele-1)*Mdims%cv_nloc+iloc)
                  if (have_option('/porous_media/porous_properties/tensor_field::porous_thermal_conductivity/Wiener_conductivity')) then
                    do idim = 1, Mdims%ndim
                      ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase ) = &
                      ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase )+ saturation%val(1, iphase, cv_inod) * &
                      (sfield%val(ele_nod) * node_val( diffusivity, idim, idim, mat_inod ) &
                      +(1.0-sfield%val(ele_nod))* tfield%val(idim, idim, t_ele_nod)) ! for classic weighted approach (Wiener approach)
                    end do
                  else
                    do idim = 1, Mdims%ndim
                      ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase ) = &
                      ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase )+ saturation%val(1, iphase, cv_inod) * &
                      (tfield%val(idim, idim, t_ele_nod)+3*tfield%val(idim, idim, t_ele_nod)* &
                      (node_val( diffusivity, idim, idim, mat_inod ) - tfield%val(idim, idim, t_ele_nod))*sfield%val(ele_nod)/ &
                      (3*tfield%val(idim, idim, t_ele_nod)+(node_val( diffusivity, idim, idim, mat_inod )-tfield%val(idim, idim, t_ele_nod))* &
                      (1-sfield%val(ele_nod))))
                    end do
                  end if
                end do
              end do
            end do
          endif
        else
          do iphase = 1, Mdims%nphase

            if (present(TracerName)) then
              diffusivity => extract_tensor_field( state(iphase), trim(TracerName)//'Diffusivity', stat )
            else
              diffusivity => extract_tensor_field( state(iphase), 'TemperatureDiffusivity', stat )
            endif
            if (stat /= 0) cycle!If no field defined then cycle
            do ele = 1, Mdims%totele
              !                     ele_nod = min(size(sfield%val), ele)
              !t_ele_nod = min(size(tfield%val, 3), ele)
              do iloc = 1, Mdims%mat_nloc
                mat_inod = ndgln%mat( (ele-1)*Mdims%mat_nloc + iloc )
                cv_inod = ndgln%cv((ele-1)*Mdims%cv_nloc+iloc)
                do idim = 1, Mdims%ndim
                  ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase ) = &
                  ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase )+&
                  node_val( diffusivity, idim, idim, mat_inod )
                end do
              end do
            end do
          end do
          !do iphase = 1, Mdims%nphase
          !    diffusivity => extract_tensor_field( state(iphase), 'TemperatureDiffusivity', stat )
          !    do idim = 1, Mdims%ndim
          !        ScalarAdvectionField_Diffusion( :, idim, idim, iphase ) = node_val( diffusivity, idim, idim, iphase )
          !    end do
          !end do
        end if
      end if
      if ( harmonic_average ) then
        ! ScalarAdvectionField_Diffusion = 1.0 / ScalarAdvectionField_Diffusion
        do iphase = 1, Mdims%nphase
          do idim = 1, Mdims%ndim
            do mat_inod = 1, Mdims%mat_nonods
              if ( ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase ) > 0.0 ) &
              ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase ) = &
              1.0 / ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase )
            end do
          end do
        end do
      end if
      do iphase = 1, Mdims%nphase
        if (present(TracerName)) then
          ewrite(3,*) trim(TracerName), 'diffusivity min_max', iphase, &
          minval( ScalarAdvectionField_Diffusion( :, 1, 1, iphase ) ), &
          maxval( ScalarAdvectionField_Diffusion( :, 1, 1, iphase ) )
        else
          ewrite(3,*) 'Thermal conductivity min_max', iphase, &
          minval( ScalarAdvectionField_Diffusion( :, 1, 1, iphase ) ), &
          maxval( ScalarAdvectionField_Diffusion( :, 1, 1, iphase ) )
        endif
      end do

      !If we want to normlise the equation by rho CP we can return the diffusion coefficient divided by rho Cp
      if (present_and_true(divide_by_rho_CP)) then
        tfield => extract_tensor_field( packed_state, "PackedDensityHeatCapacity", stat )
        do iphase = 1, Mdims%nphase
          do ele = 1, Mdims%totele
            do iloc = 1, Mdims%mat_nloc
              mat_inod = ndgln%mat( (ele-1)*Mdims%mat_nloc + iloc )
              cv_inod = ndgln%cv((ele-1)*Mdims%cv_nloc+iloc)
              do idim = 1, Mdims%ndim
                ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase ) = &
                ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase )/tfield%val(1, iphase, cv_inod)
              end do
            end do
          end do
        end do
      end if

      return
    end subroutine calculate_diffusivity

    !>@brief: Dispersion for porous media
    !> For thermal, the field density needs to be passed down, which ensures that even for boussinesq a reference density is still used
    subroutine calculate_solute_dispersity(state, packed_state, Mdims, ndgln, density, SoluteDispersion)
      type(state_type), dimension(:), intent(in) :: state
      type( state_type ), intent( inout ) :: packed_state
      type(multi_dimensions), intent(in) :: Mdims
      type(multi_ndgln), intent(in) :: ndgln
      real, dimension(:,:), target, intent(in) :: density !>If present density means that we need to use a reference density to ensure consistency with the rest of the equation
      real, dimension(:, :, :, :), intent(inout) :: SoluteDispersion
      !Local variables
      type(scalar_field), pointer :: component, sfield, ldfield, tdfield
      type (vector_field_pointer), dimension(Mdims%n_in_pres) ::darcy_velocity
      integer :: icomp, iphase, idim, stat, ele, idim1, idim2
      integer :: iloc, mat_inod, cv_inod, ele_nod, t_ele_nod, u_iloc, u_nod, u_nloc, cv_loc, cv_iloc, ele_nod_disp
      real :: vel_av
      real, dimension(3, 3) :: DispCoeffMat
      real, dimension(3) :: vel_comp, vel_comp2, DispDiaComp
      real, dimension(:), pointer :: tdisp


      SoluteDispersion = 0.
      DispCoeffMat = 0.
      DispDiaComp = 0.

      sfield=>extract_scalar_field(state(1),"Porosity")
      ldfield=>extract_scalar_field(state(1),"Longitudinal_Dispersivity")

      if (have_option("/porous_media/Dispersion/scalar_field::Transverse_Dispersivity")) then
        tdfield=>extract_scalar_field(state(1),"Transverse_Dispersivity")
        tdisp => tdfield%val
      else
        allocate(tdisp(size(ldfield%val)))
        tdisp = ldfield%val / 10.
      end if


      do iphase = 1, Mdims%n_in_pres
        ! if ( .not. have_option( "/material_phase["// int2str( iphase - 1 )//"]/scalar_field::"//trim(TracerName)//"/prognostic/tensor_field::Diffusivity")) cycle

        darcy_velocity(iphase)%ptr => extract_vector_field(state(iphase),"DarcyVelocity")

        do ele = 1, Mdims%totele
          ele_nod = min(size(sfield%val), ele)
          ele_nod_disp = min(size(ldfield%val), ele)
          do u_iloc = 1, mdims%u_nloc
            u_nod = ndgln%u(( ELE - 1) * Mdims%u_nloc + u_iloc )
            do cv_iloc = 1, Mdims%cv_nloc
              mat_inod = ndgln%mat((ele-1)*Mdims%mat_nloc+cv_iloc)
              cv_loc = ndgln%cv((ele-1)*Mdims%cv_nloc+cv_iloc)
              vel_av = 0

              do idim1 = 1, Mdims%ndim
                vel_comp2(idim1) = ((darcy_velocity(iphase)%ptr%val(idim1,u_nod))/&
                (sfield%val(ele_nod)))**2

                vel_comp(idim1) = ((darcy_velocity(iphase)%ptr%val(idim1,u_nod))/&
                (sfield%val(ele_nod)))

                if (Mdims%ndim == 2) then
                  vel_comp2(3) = 0
                  vel_comp(3) = 0
                endif

                vel_av = vel_av + vel_comp2(idim1)
              end do

              vel_av = SQRT(vel_av)

              do idim1 = 1, Mdims%ndim
                do idim2 = 1, Mdims%ndim
                  if (idim1 == idim2) then
                    DispCoeffMat(idim1, idim2) = vel_av * ldfield%val(ele_nod_disp)
                  else
                    DispCoeffMat(idim1, idim2) = vel_av * tdisp(ele_nod_disp)
                  endif
                end do
              end do

              !! Diagonal components of the dispersion tensor
              DispDiaComp = (1.0/(vel_av**2))*matmul(DispCoeffMat, vel_comp2)

              !! Off-diaginal components of the dispersion tensor
              do idim1 = 1, Mdims%ndim
                do idim2 = 1, Mdims%ndim
                  if (idim1 .NE. idim2) then
                    SoluteDispersion( mat_inod, idim1, idim2, iphase ) =&
                    sfield%val(ele_nod)*(1.0/(vel_av**2)) *&
                    ((vel_av * ldfield%val(ele_nod_disp)) - (vel_av * tdisp(ele_nod_disp))) *&
                    (ABS(vel_comp(idim1)) * ABS(vel_comp(idim2)))
                  else
                    SoluteDispersion( mat_inod, idim1, idim2, iphase ) =&
                    sfield%val(ele_nod)*DispDiaComp(idim1)
                  endif
                  SoluteDispersion( mat_inod, idim1, idim2, iphase ) =&
                  SoluteDispersion( mat_inod, idim1, idim2, iphase ) *&
                  density(iphase, cv_loc)
                end do
              end do
            end do
          end do
        end do
      end do
      if (.not.have_option("/porous_media/Dispersion/scalar_field::Transverse_Dispersivity")) then
        deallocate(tdisp)
      end if

      return
    end subroutine calculate_solute_dispersity

    !>@brief: Computes the viscosity effect as a momemtum diffusion, this is zero for porous media
    subroutine calculate_viscosity( state, Mdims, ndgln, Momentum_Diffusion, Momentum_Diffusion2 )
      implicit none
      type( multi_dimensions ), intent( in ) :: Mdims
      type( multi_ndgln ), intent( in ) :: ndgln
      type( state_type ), dimension( : ), intent( in ) :: state
      real, dimension( :, :, :, : ), intent( inout ) :: Momentum_Diffusion
      type( multi_field ), intent( inout ) :: Momentum_Diffusion2

      !Local variables
      type( tensor_field ), pointer :: t_field, tp_field, tc_field
      integer :: iphase, icomp, stat, mat_nod, cv_nod, ele
      type( scalar_field ), pointer :: component, saturation
      logical :: linearise_viscosity, cg_mesh
      real, dimension( : ), allocatable :: component_tmp
      real, dimension( :, :, : ), allocatable :: mu_tmp
      integer :: iloc, ndim1, ndim2, idim, jdim
      integer :: multiplier

      real :: exp_zeta_function

  ! DELETE Momentum_Diffusion - START USING THE NEW MEMORY ---

      if ( is_porous_media) then
         momentum_diffusion=0.0
      else
        if (have_option('/material_phase[0]/phase_properties/Viscosity/tensor_field::Viscosity/prescribed/value::WholeMesh/isotropic') &
         .AND. have_option('/material_phase[0]/phase_properties/Viscosity/viscosity_scheme/stress_form') ) then
        print *, "WARNING: PLEASE ENSURE THAT YOU USE ANISOTROPIC SYMMETRIC WHEN USING STRESS FORM OF VISCOSITY, Otherwise your results are likely to be wrong"
        end if
         momentum_diffusion=0.0
         t_field => extract_tensor_field( state( 1 ), 'Viscosity', stat )
         !Multiplier to control the index for the viscosity when the viscosity is constant
         multiplier = 1
         if (size(t_field%val,3) == 1)  multiplier = 0

         if ( stat == 0 ) then
            linearise_viscosity = have_option( '/material_phase[0]/linearise_viscosity' )
            allocate( component_tmp( Mdims%cv_nloc ), mu_tmp( t_field%dim(1), t_field%dim(2), Mdims%cv_nloc ) )
            if ( Mdims%ncomp > 1 ) then
               t_field%val=0.0
               do icomp = 1, Mdims%ncomp
                  do iphase = 1, Mdims%nphase
                     component => extract_scalar_field( state(Mdims%nphase + icomp), 'ComponentMassFractionPhase' // int2str(iphase) )
                     tc_field => extract_tensor_field( state( Mdims%nphase + icomp ), 'Viscosity' )
                     tp_field => extract_tensor_field( state( iphase ), 'Viscosity' )

                     ewrite(3,*) 'Component, Phase, Visc_min_max', icomp, iphase, minval( tc_field%val ), maxval( tc_field%val )
                     do ele = 1, ele_count( tc_field )
                        component_tmp = ele_val( component, ele )
                        mu_tmp = ele_val( tc_field, ele )
                        do iloc = 1, Mdims%cv_nloc
                           mu_tmp( :, :, iloc ) = mu_tmp( :, :, iloc ) * component_tmp( iloc )
                        end do
                        if ( linearise_viscosity ) then
                           mu_tmp( :, :, 2 ) = 0.5 * ( mu_tmp( :, :, 1 ) + mu_tmp( :, :, 3 ) )
                           mu_tmp( :, :, 4 ) = 0.5 * ( mu_tmp( :, :, 1 ) + mu_tmp( :, :, 6 ) )
                           mu_tmp( :, :, 5 ) = 0.5 * ( mu_tmp( :, :, 3 ) + mu_tmp( :, :, 6 ) )
                           if ( Mdims%cv_nloc == 10 ) then
                              mu_tmp( :, :, 7 ) = 0.5 * ( mu_tmp( :, :, 1 ) + mu_tmp( :, :, 10 ) )
                              mu_tmp( :, :, 8 ) = 0.5 * ( mu_tmp( :, :, 3 ) + mu_tmp( :, :, 10 ) )
                              mu_tmp( :, :, 9 ) = 0.5 * ( mu_tmp( :, :, 6 ) + mu_tmp( :, :, 10 ) )
                           end if
                        end if
                        do iloc = 1, Mdims%cv_nloc
                           cv_nod = ndgln%cv( (ele-1)*Mdims%cv_nloc + iloc )
                           mat_nod = ndgln%mat( (ele-1)*Mdims%cv_nloc + iloc )
                           momentum_diffusion( :, :, iphase, mat_nod ) = momentum_diffusion(  :, :, iphase, mat_nod ) + mu_tmp( 1, 1, iloc ) ! isotropic only - to be deleted...
                           cv_nod = cv_nod * multiplier + (1 - multiplier)!index has to be one if viscosity is constant
                           t_field%val( :, :, cv_nod ) = t_field%val( :, :, cv_nod ) + mu_tmp( :, :, iloc )/dble(Mdims%cv_nloc)
                        end do
                     end do
                  end do
               end do
            else
               if (is_magma)  then
                 saturation => extract_scalar_field(state(2), "PhaseVolumeFraction") !HH melt is the 2nd phase
                 call get_option('/magma_parameters/bulk_viscosity_exponential_coefficient' , exp_zeta_function)
               end if
               cg_mesh = have_option( '/material_phase[0]/phase_properties/Viscosity/tensor_field::Viscosity/diagnostic/mesh::PressureMesh')
               do iphase = 1, Mdims%nphase
                  tp_field => extract_tensor_field( state( iphase ), 'Viscosity', stat )
                  do ele = 1, ele_count( tp_field )
                     mu_tmp = ele_val( tp_field, ele )
                     if ( linearise_viscosity ) then
                        mu_tmp( :, :, 2 ) = 0.5 * ( mu_tmp( :, :, 1 ) + mu_tmp( :, :, 3 ) )
                        mu_tmp( :, :, 4 ) = 0.5 * ( mu_tmp( :, :, 1 ) + mu_tmp( :, :, 6 ) )
                        mu_tmp( :, :, 5 ) = 0.5 * ( mu_tmp( :, :, 3 ) + mu_tmp( :, :, 6 ) )
                        if ( Mdims%cv_nloc == 10 ) then
                           mu_tmp( :, :, 7 ) = 0.5 * ( mu_tmp( :, :, 1 ) + mu_tmp( :, :, 10 ) )
                           mu_tmp( :, :, 8 ) = 0.5 * ( mu_tmp( :, :, 3 ) + mu_tmp( :, :, 10 ) )
                           mu_tmp( :, :, 9 ) = 0.5 * ( mu_tmp( :, :, 6 ) + mu_tmp( :, :, 10 ) )
                        end if
                     end if
                     do iloc = 1, Mdims%cv_nloc
                        mat_nod = ndgln%mat( (ele-1)*Mdims%cv_nloc + iloc )
                        cv_nod = ndgln%cv( (ele-1)*Mdims%cv_nloc + iloc )
                        momentum_diffusion( :, :, iphase, mat_nod ) = mu_tmp( :, :, iloc )
                        !Currently only magma uses momentum_diffusion2
                        if (iphase==1 .and. is_magma) then !only the solid phase has bulk viscosity
                          momentum_diffusion2%val(iphase, 1, 1, mat_nod)  = zeta(mu_tmp( 1, 1, iloc ), exp_zeta_function, saturation%val(cv_nod))*0 !now turned off
                        end if
                        if(cg_mesh) then
                          mat_nod = cv_nod * multiplier + (1 - multiplier)! this is for CG
                        else
                          mat_nod = mat_nod * multiplier + (1 - multiplier)! this is for DG
                        end if
                        ! if ( have_option( '/blasting' ) ) then
                        !    t_field%val( :, :, 1 ) = mu_tmp( :, :, iloc )
                        ! else
                        !    t_field%val( :, :, mat_nod ) = mu_tmp( :, :, iloc )
                        ! end if

                     end do
                  end do
               end do
            end if
            deallocate( component_tmp, mu_tmp )
         end if
      end if


      !!! NEW CODE HERE !!!
      !!! deal with Momentum_Diffusion2



      return
    Contains
      !---------------------------------------------------------------------------
      !> @author Haiyang Hu
      !> @brief TO BE FILLED BY THE AUTHOR
      !---------------------------------------------------------------------------
      real function zeta(a,n,phi)
        implicit none
        ! bulk viscosity zeta=a*phi^-n
        real :: a!> TO BE FILLED
        real :: n!> TO BE FILLED
        real :: phi!> TO BE FILLED
        zeta=a*phi**(-n)
        return
      end function zeta

    end subroutine calculate_viscosity



    !sprint_to_do, re-use material_absoprtion by updating the values of the input absoprtion
    !>@brief:Computes velocity absorption from diamond information
    subroutine update_velocity_absorption( states, ndim, nphase, velocity_absorption )

        implicit none

        integer, intent( in ) :: ndim, nphase
        type( state_type ), dimension( : ), intent( in ) :: states
        real, dimension( :, :, : ), intent( inout ) :: velocity_absorption

        type( vector_field ), pointer :: absorption
        integer :: iphase, idim
        logical :: have_absorption
        character( len = option_path_len ) :: option_path

        velocity_absorption = 0.

        do iphase = 1, nphase
            have_absorption = .false.
            option_path = '/material_phase[' // int2str( iphase - 1 ) // ']/vector_field::Velocity' // &
                '/prognostic/vector_field::Absorption/diagnostic/algorithm::vector_python_diagnostic'
            have_absorption = have_option( trim(option_path) )
            if (.not.have_absorption) then!Test if it is prescribed and constant
                option_path = '/material_phase[' // int2str( iphase - 1 ) // ']/vector_field::Velocity' // &
                    '/prognostic/vector_field::Absorption/prescribed/value'
                have_absorption = have_option( trim(option_path) )
            end if

            if ( have_absorption ) then
                absorption => extract_vector_field( states( iphase ), 'VelocityAbsorption' )
                if (size(velocity_absorption,3) == size(absorption % val,2)) then
                    do idim = 1, ndim
                        velocity_absorption( idim + (iphase-1)*ndim, idim + (iphase-1)*ndim, : ) =  &
                            absorption % val( idim, : )
                    end do
                else
                    FLAbort(" The velocity absorption field has to be on the same mesh as velocity")
                    ! The code below doesn't interpolate from the absorption mesh to the velocity mesh
                    do idim = 1, ndim
                        velocity_absorption( idim + (iphase-1)*ndim, idim + (iphase-1)*ndim, : ) =  &
                            absorption % val( idim, size(absorption % val,2) )
                    end do
                end if
            else
                do idim = 1, ndim
                    velocity_absorption( idim + (iphase-1)*ndim, idim + (iphase-1)*ndim, : ) = 0.
                end do
            end if
        end do

        return
    end subroutine update_velocity_absorption

    !sprint_to_do, re-use material_absoprtion by updating the values of the input absoprtion
    !>@brief:Computes velocity absorption from diamond information
    subroutine update_velocity_absorption_coriolis( states, ndim, nphase, velocity_absorption )

      implicit none

      integer, intent( in ) :: ndim, nphase
      type( state_type ), dimension( : ), intent( in ) :: states
      real, dimension( :, :, : ), intent( inout ) :: velocity_absorption

      type( scalar_field ), pointer :: f
      integer :: iphase, stat, idx1, idx2

      do iphase = 1, nphase
         f => extract_scalar_field( states( iphase ), 'f', stat )
         if ( stat == 0 ) then
            idx1 = 1 + (iphase-1)*ndim ;  idx2 = 2 + (iphase-1)*ndim
            velocity_absorption( idx1, idx2, : ) = velocity_absorption(idx1, idx2, : ) - f % val
            velocity_absorption( idx2, idx1, : ) = velocity_absorption(idx2, idx1, : ) + f % val
         end if
      end do

      return
    end subroutine update_velocity_absorption_coriolis


    !>@brief:Computes velocity source from diamond information
    subroutine update_velocity_source( states, Mdims, u_source )

        implicit none

        type( multi_dimensions ), intent( in ) :: Mdims
        type( state_type ), dimension( : ), intent( in ) :: states
        !real, dimension( :, :, : ), intent( inout ) :: u_source
        type ( multi_field ), intent( inout ) :: u_source

        type( vector_field ), pointer :: source
        integer :: iphase, idim
        logical, dimension(Mdims%nphase) :: have_source
        character( len = option_path_len ) :: option_path

        ! Needs to be cleaned up !! SPRINT TO DO..


        !option_count("/material_phase/vector_field::Velocity/prognostic/vector_field::Source")

        !u_source = 0.

        have_source = .false.
        do iphase = 1, Mdims%nphase
            option_path = '/material_phase[' // int2str( iphase - 1 ) // ']/vector_field::Velocity' // &
                '/prognostic/vector_field::Source'
            have_source(iphase) =  have_option( trim(option_path) )
        end do

        do iphase = 1, Mdims%nphase
            if ( have_source(iphase) ) then
                source => extract_vector_field( states( iphase ), 'VelocitySource' )
                do idim = 1, Mdims%ndim
                    call assign_val(u_source%val( idim, iphase, 1, : ), source % val( idim, : ))
                end do
            end if
        end do

        return
    end subroutine update_velocity_source




    !>@brief: ????
    real function saturation_temperature( pressure )
        implicit none
        real :: pressure
        real :: p, c, pr

        p = pressure*10.0
        c=0.5
        if (p>1.0e6) c=0.3
        pr = 1.56e6 + c*(p - 1.0e6)
        saturation_temperature = (500.0*2.0/pi)*atan(0.5*pi*(pr-5.0e4)*1.0e-6)-273.15

        return
    end function saturation_temperature

    !>@brief:Gets the relperm max, the relperm exponent and the immobile fractions and stores
    !>them into packed state
    !>By index this is: 1) immobile fraction, 2) relperm max, 3)relperm exponent
    !> 4)Capillary entry pressure 5) Capillary exponent 6) Capillary imbition term
    !> The effective inmobile fraction is the min(inmobile,saturation),
    !> being the saturation the value after a succesful non-linear solver convergence!
    !> This NEEDS to be called after a succesful non-linear solver (with update_only)
    subroutine get_RockFluidProp(state, packed_state, Mdims, ndgln, current_time, update_only)
        implicit none
        type( multi_dimensions ), intent( in ) :: Mdims
        type(state_type), dimension(:), intent(inout) :: state
        type( multi_ndgln ), intent( in ) :: ndgln
        type( state_type ), intent( inout ) :: packed_state
        real, optional, intent(in) :: current_time
        logical, optional, intent(in) :: update_only
        !Local variables
        type (tensor_field), pointer :: t_field, Saturation, SaturationOld
        type (scalar_field), target :: targ_Store
        type (scalar_field), pointer :: s_field, saturation_flip
        type (vector_field), pointer :: position
        type(mesh_type), pointer :: fl_mesh
        type(mesh_type) :: Auxmesh
        real :: auxR
        integer :: iphase, nphase, ele, cv_iloc, cv_nod
        character(len=500) :: path, path2, path3

        t_field=>extract_tensor_field(packed_state,"PackedRockFluidProp")
        Saturation => extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
        SaturationOld => extract_tensor_field(packed_state,"PackedOldPhaseVolumeFraction")
        nphase = size(t_field%val,2)
        !By default the pressure mesh (position 1)
        s_field => extract_scalar_field(state(1),1)
        position => get_external_coordinate_field(packed_state, s_field%mesh)

        fl_mesh => extract_mesh( state(1), "P0DG" )
        Auxmesh = fl_mesh
        call allocate (targ_Store, Auxmesh, "Temporary_get_RockFluidProp")

        !If only updating there is no need to update the other parameters
        if (.not.present_and_true(update_only)) then
          !Now obtain relpermMax
          do iphase = 1, nphase
            path = "/material_phase["//int2str(iphase-1)//"]/multiphase_properties/Relperm_Corey/scalar_field::relperm_max/prescribed/value"
            if (have_option(trim(path))) then
              call initialise_field_over_regions(targ_Store, trim(path), position)
              t_field%val(2,iphase,:) = max(min(targ_Store%val, 1.0), 0.0)
            else !default value
              t_field%val(2,iphase,:) = 1.0
            end if
          end do

          !Retrieve relperm exponent
          do iphase = 1, nphase
              path = "/material_phase["//int2str(iphase-1)//&
                  "]/multiphase_properties/Relperm_Corey/scalar_field::relperm_exponent/prescribed/value"
              if (have_option(trim(path))) then
                  call initialise_field_over_regions(targ_Store, trim(path) , position)
                  t_field%val(3,iphase,:) = targ_Store%val
              else !default value
                  t_field%val(3,iphase,:) = 2.0
              end if
          end do
          !Initialize capillary pressure
          if (have_option_for_any_phase( '/multiphase_properties/capillary_pressure', nphase ) ) then
              !Get cap pressure constant, C
              do iphase = 1, nphase
                  path = "/material_phase["//int2str(iphase-1)//&
                      "]/multiphase_properties/capillary_pressure/type_Brooks_Corey/scalar_field::C/prescribed/value"
                  path3 = "/material_phase["//int2str(iphase-1)//&
                      "]/multiphase_properties/capillary_pressure/type_Power_Law/scalar_field::C/prescribed/value"
                  if (have_option(trim(path))) then
                      call initialise_field_over_regions(targ_Store, trim(path) , position)
                      t_field%val(4,iphase,:) = targ_Store%val
                  elseif (have_option(trim(path3))) then
                      call initialise_field_over_regions(targ_Store, trim(path3) , position)
                      t_field%val(4,iphase,:) = targ_Store%val
                  else !default value
                      t_field%val(4,iphase,:) = 0.0
                  end if
              end do

              !Get cap exponent, a
              do iphase = 1, nphase
                  path = "/material_phase["//int2str(iphase-1)//&
                      "]/multiphase_properties/capillary_pressure/type_Brooks_Corey/scalar_field::a/prescribed/value"
                  path3 = "/material_phase["//int2str(iphase-1)//&
                      "]/multiphase_properties/capillary_pressure/type_Power_Law/scalar_field::a/prescribed/value"
                  if (have_option(trim(path))) then
                      call initialise_field_over_regions(targ_Store, trim(path) , position)
                      t_field%val(5,iphase,:) = targ_Store%val
                  elseif (have_option(trim(path3))) then
                      call initialise_field_over_regions(targ_Store, trim(path3) , position)
                      t_field%val(5,iphase,:) = targ_Store%val
                  else !default value
                      t_field%val(5,iphase,:) = 1.0
                  end if
              end do
              !Get imbibition term
              do iphase = 1, nphase
                  path = "/material_phase["//int2str(iphase-1)//&
                      "]/multiphase_properties/capillary_pressure/type_Brooks_Corey/scalar_field::B/prescribed/value"
                  if (have_option(trim(path))) then
                      call initialise_field_over_regions(targ_Store, trim(path) , position)
                      t_field%val(6,iphase,:) = targ_Store%val
                  else !default value
                      t_field%val(6,iphase,:) = 0.0
                  end if
              end do
          end if
        end if
        !Retrieve Immobile fractions
        do iphase = 1, nphase
          path = "/material_phase["//int2str(iphase-1)//"]/multiphase_properties/immobile_fraction/scalar_field::value/prescribed/value"
            if (have_option(trim(path))) then
                call initialise_field_over_regions(targ_Store, trim(path) , position)
                t_field%val(1,iphase,:) = targ_Store%val
            else if (have_option("/material_phase["//int2str(iphase-1)//&
                      "]/multiphase_properties/immobile_fraction/scalar_field::Land_coefficient/prescribed/value")) then
              path = "/material_phase["//int2str(iphase-1)//&
                "]/multiphase_properties/immobile_fraction/scalar_field::Land_coefficient/prescribed/value"
              !Extract the land parameter
              call initialise_field_over_regions(targ_Store, trim(path) , position)
              !We first extract the field containing the historical point of saturation
              saturation_flip => extract_scalar_field(state(iphase), "Saturation_flipping")
              !Only for the first time ever, not for checkpointing, overwrite the saturation flipping value with the initial one
              if (present(current_time)) then
               if( current_time < 1e-8) then
                do cv_nod = 1, Mdims%cv_nonods
                  saturation_flip%val(cv_nod) = max(Saturation%val(1,iphase,cv_nod), 1e-16)!limit because we need to store signs also
                end do
               end if
              end if
              do ele = 1, Mdims%totele
                t_field%val(1,iphase,ele) = 0.
                do cv_iloc = 1, Mdims%cv_nloc
                  cv_nod = ndgln%cv((ele-1)*Mdims%cv_nloc + cv_iloc)
                  !Then the immobile fraction depends on the Land coefficient as follows (this must occur outside the non-linear loop!)
                  !Formula is: Immobile = S_flip/(1+C*S_flip). Where S_flip is the saturation
                  !when changing from imbibition to drainage or the other way round
                  call Update_saturation_flipping(saturation_flip%val(cv_nod), Saturation%val(1,iphase,cv_nod), SaturationOld%val(1,iphase,cv_nod))
                  auxR = abs(saturation_flip%val(cv_nod))
                  t_field%val(1,iphase,ele) = t_field%val(1,iphase,ele) + (auxR/(1. + targ_Store%val(ele) * auxR))/Mdims%cv_nloc
                end do
              end do
            else !default value for immiscible values
                t_field%val(1,iphase,:) = 0.0
            end if
        end do
        call deallocate(targ_Store)

    contains

    !>@brief:  This internal subroutine checks the if we are flipping from drainage to imbibition, or the other way round,
    !> and updates if required the value stored in Saturation_flipping
    !> Saturation_flipping stores both the value and the history, being positive if the phase is increasing and negative if the phase is decreasing.
    !> Therefore its minimum absolute value is non-zero
    !> NOTE: Currently the trapping can only increase, i.e. no thermal effects have been considered
    subroutine Update_saturation_flipping(sat_flip, sat, old_Sat)
      implicit none
      real, INTENT(IN) :: sat, old_Sat
      real, INTENT(INOUT) :: sat_flip
      !Local variables
      real, parameter :: tol = 1e-10
      !Check if the situation is changing and if so, store the new value with the sign
      ! Ensure that the immobile fraction does not decrease, i.e. sat_flip does not decrease
      ! this can only decrease once it has trapped a field with thermal effects, 
      ! but currently we are not considering these
      if (old_sat > abs(sat_flip) + tol ) then
         if (abs(sign(1., sat - old_sat ) - sign(1., sat_flip )) > tol ) &
          sat_flip = sign(old_sat, sat - old_sat )
      end if

    end subroutine

    end subroutine get_RockFluidProp
    !!JWL equation functions

    function JWL( A, B, w, R1, R2, E0, p,  roe, ro) result(fro)
        implicit none
        real, intent( in ) ::  A, B, w, R1, R2, E0, p,  roe, ro
        real :: fro
        real :: V
        V=roe/ro
        fro=(A*(1.0-w/(R1*V))*exp(-R1*V)+B*(1.0-w/(R2*V))*exp(-R2/V)+w*E0/V)-p
    end function JWL



    function diffJWL(A, B, w, R1, R2, E0, roe, ro)  result(difffro)
        implicit none
        real, intent( in ) ::  A, B, w, R1, R2, E0, roe, ro
        real ::  difffro

        difffro=(E0*w)/roe + (B*R2*exp(-(R2*ro)/roe)*((ro*w)/(R2*roe) - 1.0))/roe- (A*w*exp(-(R1*roe)/ro))/(R1*roe) - (B*w*exp(-(R2*ro)/roe))/(R2*roe)- (A*R1*roe*exp(-(R1*roe)/ro)*((ro*w)/(R1*roe) - 1.0))/ro**2.0

    end function diffJWL



    function JWLdensity(eos_coefs, pressure, ro0, JWLn) result(Rho)
        !      implicit none
        real, dimension( : ),   intent( in ) :: eos_coefs
        real, dimension( : ),   intent( in ) :: pressure
        real, dimension( : ),   intent( in ) :: ro0
        integer, intent( in ) :: JWLn

        real, dimension( JWLn ) :: Rho

        integer :: JWLi, JWLj

        real, dimension(JWLn) :: rozero
        rozero=ro0

        !      allocate(eos_coefs(7));
        !      allocate(pressure(JWLn));
        !      allocate(ro0(JWLn));
        !      allocate(Rho(JWLn));


        !      do JWLi=1, JWLn
        !          if (pressure%val(JWLi)<1e6) then
        !              Rho(JWLi)=1.2
        !          else
        !              Rho=JWLdensity(eos_coefs, pressure%val, ro0, JWLn)
        !          end if
        !      end do



        do JWLi=1, JWLn
            if (pressure(JWLi)<JWL(eos_coefs( 2 ), eos_coefs( 3 ), eos_coefs( 7 ), eos_coefs( 4 ), eos_coefs( 5 ), eos_coefs( 6 ), 0.0,  eos_coefs( 1 ), 1.205) ) then
                Rho(JWLi) = 1.205
            elseif(pressure(JWLi)<1.0e6) then
                Rho(JWLi)=2.5*pressure(JWLi)/210217.842

            else

                do JWLj=1, 10000
                    Rho(JWLi)=rozero(JWLi)-JWL(eos_coefs( 2 ), eos_coefs( 3 ), eos_coefs( 7 ), eos_coefs( 4 ), eos_coefs( 5 ), eos_coefs( 6 ), pressure(JWLi),  eos_coefs( 1 ), rozero(JWLi))/diffJWL( eos_coefs( 2 ), eos_coefs( 3 ), eos_coefs( 7 ), eos_coefs( 4 ), eos_coefs( 5 ), eos_coefs( 6 ), eos_coefs( 1 ), rozero(JWLi))
                    if (abs(Rho(JWLi)-rozero(JWLi))<1e-10) then
                        exit
                    end if
                    rozero(JWLi)=Rho(JWLi)
                end do

            end if
        end do



    !      deallocate(eos_coefs);
    !      deallocate(pressure);
    !      deallocate(ro0);
    !      deallocate(Rho);

    end function JWLdensity

    !>@brief: Initialising porous media models
    !!> Given a free water level (FWL) we simulate capillary gravity equilibration,
    !!> control volumes below FWL is kept at residual lighter phase saturation
    !!> This subroutine is called after each timestep and saturations overidden
    !!> below FWL with the heavier phase
    subroutine initialise_porous_media(Mdims, ndgln, packed_state, state, exit_initialise_porous_media)

        implicit none
        type(multi_dimensions), optional        :: Mdims
        type(multi_ndgln), optional             :: ndgln
        type(state_type), optional              :: packed_state
        type(state_type), dimension(:), optional:: state
        logical, intent(inout)                  :: exit_initialise_porous_media

        ! local
        real, save                      :: FWL = -1
        type(vector_field)              :: grav_vector
        integer, save                   :: grav_direc = 99, grav_sign=99, heavier_phase, lighter_phase ! down depends on gravity direction
        logical, save                   :: message = .true. ! user messages
        real, dimension(:,:), pointer   :: x_all, immobile_fraction, saturation_field, old_saturation_field, density
        integer                         :: cv_iloc, iphase, ele, xnod, cv_nod, i
        real                            :: inf_norm ! check saturation inf norm, dump checkpoint exit

        if (is_porous_initialisation) then
            if (message) then
                message = .false. ! only write at the start of simulation
                ewrite(0,*) "MODEL INITIALISATION: POROUS MEDIA EQUILIBRIATION WILL NOW BE CONDUCTED."
                ewrite(0,*) "At the end of equilibration a checkpoint will be produced"
                ewrite(0,*) "User must have the following:"
                ewrite(0,*) "- Free water level - height where capillary pressure is zero"
                ewrite(0,*) "- Capillary pressure for wetting phase"
                ewrite(0,*) "- Density of wetting phase and non-wetting phase"
                ewrite(0,*) "- Gravity magnitude and direction"

                ! check we have gravity parameters
                if (.not. have_option("/physical_parameters/gravity/magnitude") ) then
                    FLAbort(" *** GRAVITY MAGNITUDE AND DIRECTION ARE NEEDED ***")
                end if
                ! gravity can only point in one of the principle directions
                if (grav_direc == 99) then
                    grav_vector = extract_vector_field(state(1), 'GravityDirection')
                    do i = 1, size(grav_vector%val(:,1))
                        if (abs(grav_vector%val(i,1)) == 1 ) then
                            grav_direc = i
                            grav_sign = sign(1, int(grav_vector%val(i,1)))
                            if (sum(abs(grav_vector%val(:,1))) > 1.) then
                                FLAbort(" *** FOR NOW GRAVITY HAS TO POINT IN ONE OF THE PRINCIPLE DIRECTIONS ***")
                            end if
                        end if
                    end do
                end if

                ! check which phase is heavier - only works for constant density for now
                call get_var_from_packed_state(packed_state, Density = density)
                if (Mdims%nphase > 2) then
                    FLAbort(" *** ONLY WORKS FOR 2 PHASES ***")
                end if
                do i = 1, Mdims%nphase
                    if (.not. minval(density(i,:)) == maxval(density(i,:)) ) then
                        FLAbort(" *** ONLY WORKS FOR CONSTANT DENSITY ***")
                    end if
                end do

                heavier_phase = maxloc(density(:,1), DIM=1)
                lighter_phase = minloc(density(:,1), DIM=1)

                ! make sure we have the capillary pressure curve for this phase
                if (.not. have_option("/material_phase[" // int2str(heavier_phase-1) // "]/multiphase_properties/capillary_pressure")) then
                    FLAbort(" *** HEAVIER PHASE IS MISSING CAPILLARY CURVE ***")
                end if

                ! Get free water level from diamond
                if (FWL<0) then
                    call get_option("/porous_media/FWL", FWL, default = 0.0)
                end if

            else

                call get_var_from_packed_state(packed_state, Immobile_fraction = immobile_fraction, &
                    PressureCoordinate = x_all, PhaseVolumeFraction = saturation_field, OldPhaseVolumeFraction = old_saturation_field)

                ! override all saturations in control volumes below FWL with the maximum heavier phase saturation
                do ele = 1, Mdims%totele
                    do cv_iloc = 1, Mdims%cv_nloc
                        cv_nod = ndgln%cv((ele-1)*Mdims%cv_nloc + cv_iloc)
                        xnod = ndgln%x((ele-1)*Mdims%x_nloc + cv_iloc)
                        if (grav_sign*x_all(grav_direc, xnod) > grav_sign*FWL) then
                        saturation_field(heavier_phase, cv_nod) = 1 -  immobile_fraction(lighter_phase,ele) ! below FWL saturation at max water saturation
                        saturation_field(lighter_phase, cv_nod) = immobile_fraction(lighter_phase,ele) ! below FWL saturation at min oil saturation
                        end if
                    end do

                end do

                ! check if saturation has converged
                inf_norm = maxval(abs(saturation_field(heavier_phase,:) - old_saturation_field(heavier_phase,:)))
                if (inf_norm < 5.e-3) then
                    ewrite(0,*) 'initialisation porous media convergence reached'
                    ewrite(0,*) '... exiting'
                    call delete_option("/porous_media/FWL")
                    call checkpoint_simulation(state, prefix='Initialisation', protect_simulation_name=.true.,file_type='.mpml')
                    exit_initialise_porous_media = .true.
                end if

            end if
        end if

    end subroutine initialise_porous_media

    !>@brief: For boussinesq porous media we need the reference density to ensure consistency when mixing with the porous density/Cp etc.
    !!> In this subroutine we retrieve the value given a phase, component.
    real function retrieve_reference_density(state, packed_state, iphase, icomp, nphase)
      implicit none
      type( state_type ), intent( inout ) :: packed_state
      type( state_type ), dimension( : ), intent( in) :: state
      integer, intent(in) :: iphase, icomp, nphase
      !local variables
      character( len = option_path_len ) :: eos_option_path
      character( len = option_path_len ) :: option_path_comp, option_path_incomp, option_path_python
      real :: ref_rho
      type (scalar_field) :: sfield
      type (scalar_field), pointer :: pnt_sfield
      type (vector_field), pointer :: position
      !Provide input to find out EOS used
      if( icomp > 0 ) then
          eos_option_path = &
          trim( '/material_phase[' // int2str( nphase + icomp - 1 ) // &
          ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
          '/prognostic/phase_properties/Density' )
      else
          eos_option_path = trim( '/material_phase[' // int2str( iphase - 1 ) // ']/phase_properties/Density' )
      end if
      !Retrieve the equation of state path
      call Assign_Equation_of_State( eos_option_path )
      !Paths for comparison
      if ( icomp > 0 ) then
          option_path_comp = trim( '/material_phase[' // int2str( nphase + icomp - 1 ) // &
              ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
              '/prognostic/phase_properties/Density/compressible' )
          option_path_incomp = trim( '/material_phase[' // int2str(nphase + icomp - 1 ) // &
              ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
              '/prognostic/phase_properties/Density/incompressible' )
          option_path_python = trim( '/material_phase[' // int2str( nphase + icomp - 1 ) // &
              ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
              '/prognostic/phase_properties/Density/python_state' )
      else
          option_path_comp = trim( '/material_phase[' // int2str( iphase - 1 ) // &
              ']/phase_properties/Density/compressible' )
          option_path_incomp = trim( '/material_phase[' // int2str( iphase - 1 ) // &
              ']/phase_properties/Density/incompressible' )
          option_path_python = trim( '/material_phase[' // int2str( iphase - 1 ) // &
              ']/phase_properties/Density/python_state' )
      end if

      Conditional_EOS_Option: if( trim( eos_option_path ) == trim( option_path_incomp ) ) then
        !!$ Constant representation
        call get_option( trim( eos_option_path ), ref_rho )
      else if( trim( eos_option_path ) == trim( option_path_comp ) // '/stiffened_gas' ) then
        !!$ Den = C0 / T * ( P - C1 )
        call get_option( trim( eos_option_path) // '/eos_option1' , ref_rho )
      elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/linear_in_pressure' ) then
        !!$ Den = C0 * P +C1

        pnt_sfield => extract_scalar_field(state(1),1)
        position => get_external_coordinate_field(packed_state, pnt_sfield%mesh)
        call allocate (sfield, pnt_sfield%mesh, "Temporary_linear_Coefficient_B")
        !Retrieve coefficients
        call initialise_field(sfield, trim( option_path_comp )//"/linear_in_pressure/coefficient_B" , position)
        ref_rho = sfield%val(1)
        call deallocate(sfield)
      elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/linear_in_pressure/include_internal_energy' ) then
        !!$ Den = C0 * P/T +C1
        call get_option( trim( option_path_comp ) // '/linear_in_pressure/coefficient_B/constant', ref_rho )
      elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/exponential_in_pressure' ) then
        call get_option( trim( eos_option_path ) // '/coefficient_A', ref_rho )
      elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/Linear_eos' ) then
        !!$ Den = den0 * ( 1 + alpha * solute mass fraction - beta * DeltaT )
        call get_option( trim( eos_option_path ) // '/reference_density', ref_rho )
      else if( trim( eos_option_path ) == trim( option_path_comp ) // '/Temperature_Pressure_correlation' ) then
        call get_option( trim( option_path_comp ) // '/Temperature_Pressure_correlation/rho0', ref_rho)
      elseif( trim( eos_option_path ) == trim( option_path_python ) ) then
        !Here the user has to specify by hand under the boussinesq option
        if (.not.have_option(trim( option_path_python ) // '/Boussinesq_approximation/reference_density')) then
          FLAbort("Porous media with boussinesq REQUIRES a reference density for all the phases.")
        end if
        call get_option( trim( option_path_python ) // '/Boussinesq_approximation/reference_density', ref_rho)
      end if Conditional_EOS_Option
      !Copy value to send out of the function
      retrieve_reference_density = ref_rho
    end function



    !>@brief: subroutine to dissolv phase2 into phase1. Currently only for system for phase 1 = water, phase 2 = gas
    !> Dissolve instantaneously the amount introduced in diamond in mol/m3 for CO2 a reference number is 38 mol/m3.
    !> Requires the first phase to have a concentration field
    subroutine flash_gas_dissolution(packed_state, Mdims, dt)
      implicit none
      type(state_type) :: packed_state
      real, intent (in) :: dt
      type(multi_dimensions) :: Mdims
      !Local variables
      real, save  :: dissolution_parameter= -1
      real, parameter ::molar_mass= 0.0441 !in kg/mol
      type(tensor_field), pointer :: saturation_field, concentration_field, density
      type(vector_field), pointer :: MeanPoreCV, cv_volume
      real :: n_co2_diss_max, delta_n, n_co2_gas
      integer :: cv_nod, iphase, stat



      !Retrieve value from diamond
      if (dissolution_parameter < 0.) call get_option("/porous_media/Gas_dissolution", dissolution_parameter)

      saturation_field=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
      concentration_field=>extract_tensor_field(packed_state,"PackedConcentration", stat)
      if (stat /= 0) then
        FLAbort("To compute flash dissolution from second phase to the first phase, a Concentration field in the first phase is required.")
      end if
      density => extract_tensor_field(packed_state,"PackedDensity")
      MeanPoreCV=>extract_vector_field(packed_state,"MeanPoreCV")
      cv_volume=> extract_vector_field(packed_state,"CVIntegral")

      do cv_nod=1,Mdims%cv_nonods
        ! Compute maximum dissolved CO2 the CV can accomodate (in mol)
        n_co2_diss_max = MeanPoreCV%val(1,cv_nod)*cv_volume%val(1,cv_nod)*density%val(1,1,cv_nod)*dissolution_parameter*saturation_field%val(1,1,cv_nod)
        ! Compute the difference between maximum dissolved CO2 the CV can accomodate and the actual dissolved CO2 present in the CV (in mol)
        delta_n = n_co2_diss_max - concentration_field%val(1,1,cv_nod)*MeanPoreCV%val(1,cv_nod)*cv_volume%val(1,cv_nod)*saturation_field%val(1,1,cv_nod)
        ! Compute the gas CO2 present in the CV (in mol)
        n_co2_gas = saturation_field%val(1,2,cv_nod) * density%val(1,2,cv_nod)/molar_mass * MeanPoreCV%val(1,cv_nod)*cv_volume%val(1,cv_nod)
        ! Update dissolved CO2 (passive tracer, in mol/m3)
        concentration_field%val(1,1,cv_nod) = concentration_field%val(1,1,cv_nod) + min(delta_n, n_co2_gas)/(MeanPoreCV%val(1,cv_nod)*cv_volume%val(1,cv_nod))
        ! Update gas CO2 (in mol)
        n_co2_gas = n_co2_gas - min(delta_n, n_co2_gas)
        ! Update gas CO2 (change in saturation)
        saturation_field%val(1,2,cv_nod) = ( n_co2_gas*molar_mass/density%val(1,2,cv_nod) ) / (MeanPoreCV%val(1,cv_nod)*cv_volume%val(1,cv_nod))
        ! Update water saturation (change in saturation)
        saturation_field%val(1,1,cv_nod) = 1 - saturation_field%val(1,2,cv_nod)

      end do


    end subroutine flash_gas_dissolution

end module multiphase_EOS
