
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

module multiphase_EOS

    use fldebug
    use state_module
    use fields
    use state_module
    use global_parameters, only: OPTION_PATH_LEN, PYTHON_FUNC_LEN, PI, is_porous_media
    use spud
    use futils, only: int2str
    use vector_tools
    use python_state
    use Copy_Outof_State

    use shape_functions_Linear_Quadratic
    use cv_advection
    use sparse_tools
    use Multiphase_module
    use sparsity_patterns_meshes, only : get_csr_sparsity_firstorder
    use arbitrary_function

    use Field_Options, only: get_external_coordinate_field
    use initialise_fields_module, only: initialise_field_over_regions

    implicit none


contains

    subroutine Calculate_All_Rhos( state, packed_state, Mdims )

        implicit none

        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent( in ) :: Mdims

        integer, dimension( : ), pointer :: cv_ndgln
        integer :: ncomp_in, nphase, ndim, cv_nonods, cv_nloc, totele
        real, dimension( : ), allocatable :: Rho, dRhodP, Density_Bulk, DensityCp_Bulk, &
             Density_Component, Cp, Component_l, c_cv_nod
        character( len = option_path_len ), dimension( : ), allocatable :: eos_option_path
        type( tensor_field ), pointer :: PackedDRhoDPressure ! (nphase, cv_nonods)
        type( tensor_field ), pointer :: field1, field2, field3, field4
        type( scalar_field ), pointer :: Cp_s
        integer :: icomp, iphase, ncomp, sc, ec, sp, ep, stat, cv_iloc, cv_nod, ele
        logical :: boussinesq

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
                      '/prognostic/equation_of_state' )
                 call Assign_Equation_of_State( eos_option_path( ( icomp - 1 ) * nphase + iphase ) )
              end do
           end do
        else
           do iphase = 1, nphase
              eos_option_path( iphase ) = trim( '/material_phase[' // int2str( iphase - 1 ) // ']/equation_of_state' )
              call Assign_Equation_of_State( eos_option_path( iphase ) )
           end do
        end if

        allocate( Rho( cv_nonods ), dRhodP( cv_nonods ) )
        allocate( Cp( cv_nonods ) ) ; Cp = 1.
        allocate( Density_Component( ncomp * nphase * cv_nonods ) )
        allocate( Density_Bulk( nphase * cv_nonods ), DensityCp_Bulk( nphase * cv_nonods ) )
        Density_Bulk = 0. ; DensityCp_Bulk = 0.

        allocate( Component_l( cv_nonods ) ) ; Component_l = 0.

        do icomp = 1, ncomp
           do iphase = 1, nphase
              sc = ( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1
              ec = ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods

              sp = ( iphase - 1 ) * cv_nonods + 1
              ep = iphase * cv_nonods

              Rho=0. ; dRhodP=0. ; Cp=1.
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

                 Density_Bulk( sp : ep ) = Density_Bulk( sp : ep ) + Rho * Component_l
                 PackedDRhoDPressure%val( 1, iphase, : ) = PackedDRhoDPressure%val( 1, iphase, : ) + dRhodP * Component_l / Rho
                 Density_Component( sc : ec ) = Rho

                 Cp_s => extract_scalar_field( state( nphase + icomp ), &
                      'ComponentMassFractionPhase' // int2str( iphase ) // 'HeatCapacity', stat )
                 if ( stat == 0 ) Cp = Cp_s % val
                 DensityCp_Bulk( sp : ep ) = DensityCp_Bulk( sp : ep ) + Rho * Cp * Component_l

              else

                 Density_Bulk( sp : ep ) = Rho
                 PackedDRhoDPressure%val( 1, iphase, : ) = dRhodP

                 Cp_s => extract_scalar_field( state( iphase ), 'TemperatureHeatCapacity', stat )
                 if ( stat == 0 ) Cp = Cp_s % val
                 DensityCp_Bulk( sp : ep ) = Rho * Cp

              end if

           end do ! iphase
        end do ! icomp

        if ( ncomp > 1 ) then
           call Cap_Bulk_Rho( state, ncomp, nphase, &
                cv_nonods, Density_Component, Density_Bulk, DensityCp_Bulk )
        end if

        field1 => extract_tensor_field( packed_state, "PackedDensity" )
        field2 => extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )
        if( ncomp > 1 ) field3 => extract_tensor_field( packed_state, "PackedComponentDensity" )

        !sf => extract_scalar_field( packed_state, "SolidConcentration" )

        do iphase = 1, nphase
           sp = ( iphase - 1 ) * cv_nonods + 1
           ep = iphase * cv_nonods

           field1 % val ( 1, iphase, : ) = Density_Bulk( sp : ep )         !* ( 1. - sf%val)     +   1000. *  sf%val
           field2 % val ( 1, iphase, : ) = DensityCp_Bulk( sp : ep )       !* ( 1. - sf%val)     +   1000. *  sf%val

           if ( ncomp > 1 ) then
              do icomp = 1, ncomp
                 sc = ( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1
                 ec = ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods
                 field3 % val ( icomp, iphase, : ) = Density_Component( sc : ec )
              end do ! icomp
           end if
        end do ! iphase

        boussinesq = have_option( "/material_phase[0]/vector_field::Velocity/prognostic/equation::Boussinesq" )
        if ( boussinesq ) field2 % val = 1.0

        deallocate( Rho, dRhodP, Cp, Component_l)
        deallocate( Density_Component, Density_Bulk, DensityCp_Bulk )
        deallocate( eos_option_path )

    end subroutine Calculate_All_Rhos


    subroutine Cap_Bulk_Rho( state, ncomp, nphase, &
        cv_nonods, Density_Component, Density, Density_Cp )

        implicit none

        type(state_type), dimension( : ) :: state
        integer, intent( in ) :: nphase, ncomp, cv_nonods
        real, dimension( cv_nonods * nphase * ncomp ), intent( in ) :: Density_Component
        real, dimension( cv_nonods * nphase ), intent( inout ) :: Density, Density_Cp

        real, dimension( :, : ), allocatable :: Density_Component_Min, Density_Component_Max
        real, dimension( :, : ), allocatable :: Density_Cp_Component_Min, Density_Cp_Component_Max
        type( scalar_field ), pointer :: Cp_s
        real, dimension( : ), allocatable :: Cp
        integer :: sp, ep, sc, ec, iphase, icomp, stat

        allocate( Density_Component_Min( nphase, cv_nonods ) ) ; Density_Component_Min = 1.e+15
        allocate( Density_Component_Max( nphase, cv_nonods ) ) ; Density_Component_Max = 0.
        allocate( Density_Cp_Component_Min( nphase, cv_nonods ) ) ; Density_Cp_Component_Min = 1.e+15
        allocate( Density_Cp_Component_Max( nphase, cv_nonods ) ) ; Density_Cp_Component_Max = 0.
        allocate( Cp( cv_nonods ) ) ; Cp = 1.

        do iphase = 1, nphase
            do icomp = 1, ncomp
                sc = ( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1
                ec = ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods

                Density_Component_Min( iphase, : ) = min( Density_Component_Min( iphase, : ), Density_Component( sc : ec ) )
                Density_Component_Max( iphase, : ) = max( Density_Component_Max( iphase, : ), Density_Component( sc : ec ) )

                Cp = 1.
                Cp_s => extract_scalar_field( state( nphase + icomp ), &
                    'ComponentMassFractionPhase' // int2str( iphase ) // 'HeatCapacity', stat )
                if( stat == 0 ) Cp = Cp_s % val

                Density_Cp_Component_Min( iphase, : ) = min( Density_Cp_Component_Min( iphase, : ), Density_Component( sc : ec ) * Cp )
                Density_Cp_Component_Max( iphase, : ) = max( Density_Cp_Component_Max( iphase, : ), Density_Component( sc : ec ) * Cp )
            end do
        end do

        do iphase = 1, nphase
            sp = ( iphase - 1 ) * cv_nonods + 1
            ep = iphase * cv_nonods

            Density( sp : ep ) = min( Density( sp : ep ), Density_Component_Max( iphase, : ) )
            Density( sp : ep ) = max( Density( sp : ep ), Density_Component_Min( iphase, : ) )

            Density_Cp( sp : ep ) = min( Density_Cp( sp : ep ), Density_Cp_Component_Max( iphase, : ) )
            Density_Cp( sp : ep ) = max( Density_Cp( sp : ep ), Density_Cp_Component_Min( iphase, : ) )
        end do

        deallocate( Cp )
        deallocate( Density_Cp_Component_Min, Density_Cp_Component_Max )
        deallocate( Density_Component_Min, Density_Component_Max )

    end subroutine Cap_Bulk_Rho


    subroutine Calculate_Component_Rho( state, packed_state, Mdims )

      implicit none

      type( state_type ), dimension( : ), intent( inout ) :: state
      type( state_type ), intent( inout ) :: packed_state
      type( multi_dimensions ), intent( in ) :: Mdims

      integer :: ncomp, nphase, cv_nonods
      real, dimension( : ), allocatable :: Rho, dRhodP
      type( tensor_field ), pointer :: field
      character( len = option_path_len ) :: eos_option_path
      integer :: icomp, iphase, s, e

      ncomp = Mdims%ncomp ; nphase = Mdims%nphase
      cv_nonods =  Mdims%cv_nonods

      allocate( Rho( cv_nonods ), dRhodP( cv_nonods ) )

      field => extract_tensor_field( packed_state, "PackedComponentDensity" )

      do icomp = 1, ncomp

         do iphase = 1, nphase
            s = ( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1
            e = ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods

            eos_option_path = trim( '/material_phase[' // int2str( nphase + icomp - 1 ) // &
                 ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
                 '/prognostic/equation_of_state' )

            call Assign_Equation_of_State( eos_option_path )
            Rho=0. ; dRhodP=0.
            call Calculate_Rho_dRhodP( state, packed_state, iphase, icomp, &
                 nphase, ncomp, eos_option_path, Rho, dRhodP )

            field % val( icomp, iphase, : ) = Rho

         end do ! iphase
      end do ! icomp
      deallocate( Rho, dRhodP )

    end subroutine Calculate_Component_Rho


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
        type( scalar_field ), pointer :: temperature, density
        character( len = option_path_len ) :: option_path_comp, option_path_incomp, option_path_python, buffer
        character( len = python_func_len ) :: pycode
        logical, save :: initialised = .false.
        logical :: have_temperature_field
        real, parameter :: toler = 1.e-10
        real, dimension( : ), allocatable, save :: reference_pressure
        real, dimension( : ), allocatable :: eos_coefs, perturbation_pressure, RhoPlus, RhoMinus
        real, dimension( : ), allocatable :: pressure_back_up, density_back_up, temperature_local
        real :: dt, current_time
        integer :: ncoef, stat

        !!$ Den = c1 * ( P + c2 ) / T           :: Stiffened EOS
        !!$ Den = c1 * P + c2                   :: Linear_1 EOS
        !!$ Den = c1 * P / T + c2               :: Linear_2 EOS
        !!$ Den = Den0 * exp[ c0 * ( P - P0 ) ] :: Exponential_1 EOS
        !!$ Den = c0 * P** c1                   :: Exponential_2 EOS

        pressure => extract_tensor_field( packed_state, 'PackedCVPressure' )
        temperature => extract_scalar_field( state( iphase ), 'Temperature', stat )
        have_temperature_field = ( stat == 0 )

        assert( node_count( pressure ) == size( rho ) )
        assert( node_count( pressure ) == size( drhodp ) )

        allocate( perturbation_pressure( node_count( pressure ) ) ) ; perturbation_pressure = 0.
        allocate( RhoPlus( node_count( pressure ) ) ) ; RhoPlus = 0.
        allocate( RhoMinus( node_count( pressure ) ) ) ; RhoMinus = 0.

        if ( ncomp > 0 ) then
            option_path_comp = trim( '/material_phase[' // int2str( nphase + icomp - 1 ) // &
                ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
                '/prognostic/equation_of_state/compressible' )
            option_path_incomp = trim( '/material_phase[' // int2str(nphase + icomp - 1 ) // &
                ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
                '/prognostic/equation_of_state/incompressible' )
            option_path_python = trim( '/material_phase[' // int2str( nphase + icomp - 1 ) // &
                ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
                '/prognostic/equation_of_state/python_state' )
        else
            option_path_comp = trim( '/material_phase[' // int2str( iphase - 1 ) // &
                ']/equation_of_state/compressible' )
            option_path_incomp = trim( '/material_phase[' // int2str( iphase - 1 ) // &
                ']/equation_of_state/incompressible' )
            option_path_python = trim( '/material_phase[' // int2str( iphase - 1 ) // &
                ']/equation_of_state/python_state' )
        end if

        Conditional_EOS_Option: if( trim( eos_option_path ) == trim( option_path_comp ) // '/stiffened_gas' ) then
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
            dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure
            deallocate( eos_coefs )



        elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/JWL_equation' ) then
            !!P=A*(1-w/(R1*V))*exp(-R1*V)+B*(1-w/(R2*V))*exp(-R2/V)+w*E0/V;
            !!The ratio V =roe/ro is defined by using roe = density of the explosive (solid part) and ro = density of the detonation products.
            allocate( eos_coefs( 7 ) ) ; eos_coefs = 0.
            !! eos_coefs(1): density_of_explosive_roe,
            !! eos_coefs(2): A,
            !! eos_coefs(3): B,
            !! eos_coefs(4): R1,
            !! eos_coefs(5): R2,
            !! eos_coefs(6): E0,
            !! eos_coefs(7): w,

            call get_option( trim( eos_option_path ) // '/density_of_explosive_roe', eos_coefs( 1 ) )
            call get_option( trim( eos_option_path ) // '/A', eos_coefs( 2 ) )
            call get_option( trim( eos_option_path ) // '/B', eos_coefs( 3 ) )
            call get_option( trim( eos_option_path ) // '/R1', eos_coefs( 4 ) )
            call get_option( trim( eos_option_path ) // '/R2', eos_coefs( 5 ) )
            call get_option( trim( eos_option_path ) // '/E0', eos_coefs( 6 ) )
            call get_option( trim( eos_option_path ) // '/w', eos_coefs( 7 ) )
            JWLn=size(pressure%val(1, 1, :));

            allocate(ro0(JWLn))
            ro0=eos_coefs(1)

            Rho=JWLdensity(eos_coefs, pressure%val(1, 1, :), ro0, JWLn)

            perturbation_pressure = max( toler, 1.e-2 * abs( pressure%val(1, 1, :) ) )

            RhoPlus=JWLdensity(eos_coefs, pressure%val(1, 1, :) + perturbation_pressure, ro0, JWLn)
            RhoMinus=JWLdensity(eos_coefs, pressure%val(1, 1, :) - perturbation_pressure, ro0, JWLn)

            dRhodP = 0.5 * (  RhoPlus - RhoMinus ) / perturbation_pressure

            do JWLi=1, JWLn
                if (pressure%val(1, 1, JWLi)<JWL(eos_coefs( 2 ), eos_coefs( 3 ), eos_coefs( 7 ), eos_coefs( 4 ), eos_coefs( 5 ), eos_coefs( 6 ), 0.0,  eos_coefs( 1 ), 1.205)) then
                    perturbation_pressure(JWLi)=1.
                    dRhodP(JWLi)=2.5
                else
                      ! perturbation_pressure(JWLi)=1.
                      ! dRhodP(JWLi)=1000.
                end if
            end do

            temperature %val=pressure %val(1, 1, :) /(Rho*278.0)
            deallocate( eos_coefs )



        elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/linear_in_pressure' ) then
            !!$ Den = C0 * P +C1
            allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
            call get_option( trim( eos_option_path ) // '/coefficient_A', eos_coefs( 1 ) )
            call get_option( trim( eos_option_path ) // '/coefficient_B', eos_coefs( 2 ) )
            Rho = eos_coefs( 1 ) * pressure % val(1,1,:) + eos_coefs( 2 )
            perturbation_pressure = 1.
            !RhoPlus = eos_coefs( 1 ) * ( pressure % val + perturbation_pressure ) + eos_coefs( 2 )
            !RhoMinus = eos_coefs( 1 ) * ( pressure % val - perturbation_pressure ) + eos_coefs( 2 )
            dRhodP = eos_coefs( 1 ) !0.5 * ( DensityPlus - DensityMinus ) / perturbation_pressure
            deallocate( eos_coefs )
        elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/linear_in_pressure/include_internal_energy' ) then
            !!$ Den = C0 * P/T +C1
            if( .not. have_temperature_field ) FLAbort( 'Temperature Field not defined' )
            allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
            call get_option( trim( option_path_comp ) // '/linear_in_pressure/coefficient_A', eos_coefs( 1 ) )
            call get_option( trim( option_path_comp ) // '/linear_in_pressure/coefficient_B', eos_coefs( 2 ) )
            Rho = eos_coefs( 1 ) * pressure % val(1,1,:) / temperature % val + eos_coefs( 2 )
            perturbation_pressure = 1.
            !RhoPlus = eos_coefs( 1 ) * ( pressure % val + perturbation_pressure ) / &
            !     ( max( toler, temperature % val ) ) + eos_coefs( 2 )
            !RhoMinus = eos_coefs( 1 ) * ( pressure % val - perturbation_pressure ) / &
            !     ( max( toler, temperature % val ) ) + eos_coefs( 2 )
            dRhodP =  eos_coefs( 1 ) / temperature % val !0.5 * ( DensityPlus - DensityMinus ) / perturbation_pressure
            deallocate( eos_coefs )

        elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/exponential_oil_gas' ) then
            !!$ Den = Den0 * Exp[ C0 * ( P - P0 ) ]
            allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
            call get_option( trim( eos_option_path ) // '/compressibility', eos_coefs( 1 ) )   ! compressibility_factor
            call get_option( trim( eos_option_path ) // '/reference_density', eos_coefs( 2 ) ) ! reference_density
            if ( .not. initialised ) then
                allocate( reference_pressure( node_count( pressure ) ) )
                reference_pressure = pressure % val(1,1,:)
                initialised = .true.
            end if
            Rho = eos_coefs( 2 ) * exp( eos_coefs( 1 ) * ( pressure % val(1,1,:) - reference_pressure ) )
            perturbation_pressure = max( toler, 1.e-3 * ( abs( pressure % val(1,1,:) ) ) )
            RhoPlus = eos_coefs( 2 ) * exp( eos_coefs( 1 ) * ( ( pressure % val(1,1,:) + perturbation_pressure ) - &
                reference_pressure ) )
            RhoMinus = eos_coefs( 2 ) * exp( eos_coefs( 1 ) * ( ( pressure % val(1,1,:) - perturbation_pressure ) - &
                reference_pressure ) )
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

        elseif( trim( eos_option_path ) == trim( option_path_incomp ) // '/linear' ) then
            !!$ Polynomial representation
            allocate( temperature_local( node_count( pressure ) ) ) ; temperature_local = 0.
            if ( have_temperature_field ) temperature_local = temperature % val
            ncoef = 10 ; allocate( eos_coefs( ncoef ) ) ; eos_coefs = 0.
            if( have_option( trim( eos_option_path ) // '/all_equal' ) ) then
                call get_option( trim( eos_option_path ) // '/all_equal', eos_coefs( 1 ) )
                eos_coefs( 2 : 10 ) = 0.
            elseif( have_option( trim( eos_option_path ) // '/specify_all' ) ) then
                call get_option( trim( eos_option_path ) // '/specify_all', eos_coefs )
            else
                FLAbort('Unknown incompressible linear equation of state')
            end if
            call Density_Polynomial( eos_coefs, pressure % val(1,1,:), temperature_local, &
                Rho )
            perturbation_pressure = max( toler, 1.e-3 * abs( pressure % val(1,1,:) ) )
            call Density_Polynomial( eos_coefs, pressure % val(1,1,:) + perturbation_pressure, temperature_local, &
                RhoPlus )
            call Density_Polynomial( eos_coefs, pressure % val(1,1,:) - perturbation_pressure, temperature_local, &
                RhoMinus )
            dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure
            deallocate( temperature_local, eos_coefs )

        elseif( trim( eos_option_path ) == trim( option_path_python ) ) then

#ifdef HAVE_NUMPY
         ewrite(3,*) "Have both NumPy and a python eos..."
#else
            FLAbort("Python eos requires NumPy, which cannot be located.")
#endif

            density => extract_scalar_field( packed_state, 'Dummy' )
            call zero( density )

            call python_reset()
            call python_add_state( packed_state )


            call python_run_string("field = state.scalar_fields['Dummy']")
            call get_option("/timestepping/current_time", current_time)
            write(buffer,*) current_time
            call python_run_string("time="//trim(buffer))
            call get_option("/timestepping/timestep", dt)
            write(buffer,*) dt
            call python_run_string("dt="//trim(buffer))

            ! Get the code
            call get_option( trim( option_path_python ) // '/algorithm', pycode )

            ! Run the code
            call python_run_string( trim( pycode ) )

            ! Copy result to protoype memory
            Rho = density % val

            ! Back up pressure and density before we start perturbing stuff...
            allocate( pressure_back_up( node_count( pressure ) ), density_back_up( node_count( pressure ) ) )
            pressure_back_up = 0. ; density_back_up = 0.
            pressure_back_up = pressure % val(1,1,:)
            density_back_up = density % val

            call python_reset()

            ! Calculating dRho / dP
            ! redefine p as p+pert and p-pert and then run python state again to get dRho / d P...
            perturbation_pressure = 1.e-5

            pressure % val(1,1,:) = pressure % val(1,1,:) + perturbation_pressure
            call zero( density )

            call python_reset()
            call python_add_state( packed_state )

            call python_run_string("field = state.scalar_fields['Dummy']")

            call get_option("/timestepping/current_time", current_time)
            write(buffer,*) current_time
            call python_run_string("time="//trim(buffer))
            call get_option("/timestepping/timestep", dt)
            write(buffer,*) dt
            call python_run_string("dt="//trim(buffer))

            call python_run_string(trim(pycode))
            RhoPlus = density % val

            call python_reset()

            pressure % val(1,1,:) = pressure_back_up
            pressure % val(1,1,:) = pressure % val(1,1,:) - perturbation_pressure
            call zero( density )

            call python_reset()
            call python_add_state( packed_state )

            call python_run_string("field = state.scalar_fields['Dummy']")

            call get_option("/timestepping/current_time", current_time)
            write(buffer,*) current_time
            call python_run_string("time="//trim(buffer))
            call get_option("/timestepping/timestep", dt)
            write(buffer,*) dt
            call python_run_string("dt="//trim(buffer))

            call python_run_string(trim(pycode))
            RhoMinus = density % val

            call python_reset()

            ! derivative
            dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure

            ! Restore pressure and density values in state
            pressure % val(1,1,:) = pressure_back_up
            density % val = density_back_up

            deallocate( pressure_back_up, density_back_up )

        else
            FLAbort( 'No option given for choice of EOS' )
        end if Conditional_EOS_Option

        deallocate( perturbation_pressure, RhoPlus, RhoMinus )

    end subroutine Calculate_Rho_dRhodP


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


    subroutine Assign_Equation_of_State( eos_option_path_out )
        implicit none
        character( len = option_path_len ), intent( inout ) :: eos_option_path_out

        Conditional_for_Compressibility: if( have_option( trim( eos_option_path_out ) // '/compressible' ) ) then
            eos_option_path_out = trim( eos_option_path_out ) // '/compressible'

            Conditional_for_Compressibility_Option: if( have_option( trim( eos_option_path_out ) // '/stiffened_gas' ) ) then
                eos_option_path_out = trim( eos_option_path_out ) // '/stiffened_gas'


            elseif( have_option( trim( eos_option_path_out ) // '/JWL_equation' ) ) then
                eos_option_path_out = trim( eos_option_path_out ) // '/JWL_equation'


            elseif( have_option( trim( eos_option_path_out ) // '/exponential_oil_gas' ) ) then
                eos_option_path_out = trim( eos_option_path_out ) // '/exponential_oil_gas'

            elseif( have_option( trim( eos_option_path_out ) // '/linear_in_pressure' ) ) then
                eos_option_path_out = trim( eos_option_path_out ) // '/linear_in_pressure'

                if( have_option( trim( eos_option_path_out ) // '/include_internal_energy' ) ) &
                    eos_option_path_out = trim( eos_option_path_out ) // '/include_internal_energy'

            elseif( have_option( trim( eos_option_path_out ) // '/exponential_in_pressure' ) ) then
                eos_option_path_out = trim( eos_option_path_out ) // '/exponential_in_pressure'

            else
                FLAbort( 'No option given for choice of EOS - compressible fluid' )

            end if Conditional_for_Compressibility_Option

        elseif( have_option( trim( eos_option_path_out ) // '/incompressible' ) )then
            eos_option_path_out = trim( eos_option_path_out ) // '/incompressible'

            Conditional_for_Incompressibility_Option: if( have_option( trim( eos_option_path_out ) // '/linear' ) ) then
                eos_option_path_out = trim( eos_option_path_out ) // '/linear'

            else
                FLAbort( 'No option given for choice of EOS - incompressible fluid' )

            end if Conditional_for_Incompressibility_Option

        elseif( have_option( trim( eos_option_path_out ) // '/python_state' ) ) then
            eos_option_path_out = trim( eos_option_path_out ) // '/python_state'

        else

            FLAbort( 'No option given for choice of EOS' )

        end if Conditional_for_Compressibility

        return
    end subroutine Assign_Equation_of_State




    subroutine Calculate_AbsorptionTerm( state, packed_state, &
        npres, cv_ndgln, mat_ndgln, &
        opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, u_absorb, IDs_ndgln, IDs2CV_ndgln )
        ! Calculate absorption for momentum eqns
        use matrix_operations
        implicit none
        type( state_type ), dimension( : ), intent( in ) :: state
        type( state_type ), intent( inout ) :: packed_state
        integer, intent( in ) :: npres
        integer, dimension( : ), intent( in ) :: cv_ndgln, mat_ndgln, IDs_ndgln, IDs2CV_ndgln
        real, dimension( :, :, :, : ), intent( inout ) :: opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new
        real, dimension( :, :, : ), intent( inout ) :: u_absorb
        !!$ Local variables:
        type( tensor_field ), pointer :: viscosity_ph
        integer :: nphase, nstate, ncomp, totele, ndim, stotel, &
            u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, x_snloc, cv_snloc, u_snloc, &
            p_snloc, cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, &
            ele, imat, icv, iphase, cv_iloc, idim, jdim, ipres, n_in_pres, loc, loc2
        real :: Mobility, pert
        real, dimension(:), allocatable :: Max_sat
        real, dimension( :, :, : ), allocatable :: u_absorb2
        !      real, dimension( : ), allocatable :: satura2
        real, dimension( :, : ), allocatable :: satura2
        real, dimension(size(state,1)) :: visc_phases
        !Working pointers
        real, dimension(:,:), pointer :: Satura, OldSatura, Immobile_fraction
        type( tensor_field ), pointer :: perm
        type( scalar_field ), pointer :: Spipe

        opt_vel_upwind_coefs_new=0.0 ; opt_vel_upwind_grad_new=0.0

        !Get from packed_state
        call get_var_from_packed_state(packed_state,PhaseVolumeFraction = Satura,&
            OldPhaseVolumeFraction = OldSatura, Immobile_fraction = Immobile_fraction)
        perm=>extract_tensor_field(packed_state,"Permeability")

        !!$ Obtaining a few scalars that will be used in the current subroutine tree:
        call Get_Primary_Scalars( state, &
            nphase, nstate, ncomp, totele, ndim, stotel, &
            u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
            x_snloc, cv_snloc, u_snloc, p_snloc, &
            cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods )

        n_in_pres = nphase / npres

        ewrite(3,*) 'In calculate_absorption'

        if( have_option( '/physical_parameters/mobility' ) )then!This option is misleading, it should be removed
            call get_option( '/physical_parameters/mobility', mobility )
            visc_phases(1) = 1
            visc_phases(2) = mobility
        elseif( have_option( '/material_phase[1]/vector_field::Velocity/prognostic/tensor_field::Viscosity' // &
            '/prescribed/value::WholeMesh/isotropic' ) ) then

            DO IPHASE = 1, NPHASE!Get viscosity for all the phases
                viscosity_ph => extract_tensor_field( state( iphase ), 'Viscosity' )
                visc_phases(iphase) = viscosity_ph%val( 1, 1, 1 )!So far we only consider scalar viscosity
            end do
            mobility = visc_phases(2) / visc_phases(1)!For backwards compatibility only
        elseif( nphase == 1 ) then
            viscosity_ph => extract_tensor_field( state( 1 ), 'Viscosity' )
            visc_phases(1) = viscosity_ph%val( 1, 1, 1 )
            mobility = visc_phases(1)
        end if

        u_absorb = 0.0

        allocate( u_absorb2( mat_nonods, nphase * ndim, nphase * ndim ), satura2( N_IN_PRES, size(SATURA,2) ) )
        u_absorb2 = 0. ; satura2 = 0.

        CALL calculate_absorption2( packed_state, CV_NONODS, N_IN_PRES, NDIM, SATURA(1:N_IN_PRES,:), TOTELE, CV_NLOC, MAT_NLOC, &
            CV_NDGLN, MAT_NDGLN, U_ABSORB(:,1:N_IN_PRES*NDIM,1:N_IN_PRES*NDIM), PERM%val, visc_phases, IDs_ndgln)

        !Introduce perturbation, positive for the increasing and negative for decreasing phase
        !Make sure that the perturbation is between bounds
        PERT = 0.0001; allocate(Max_sat(nphase))
        do icv = 1, size(satura,2)
            Max_sat(:) = 1. - sum(Immobile_fraction(:, IDs2CV_ndgln(icv))) + Immobile_fraction(:, IDs2CV_ndgln(icv))
            do iphase = 1, n_in_pres !nphase
                SATURA2(iphase, icv) = SATURA(iphase, icv) + sign(PERT, satura(iphase, icv)-OldSatura(iphase, icv))
                !If out of bounds then we perturbate in the opposite direction
                if (satura2(iphase, icv) > Max_sat(iphase) .or. &
                    satura2(iphase, icv) < Immobile_fraction(iphase, IDs2CV_ndgln(icv))) then

                    SATURA2(iphase, icv) = SATURA2(iphase, icv) - 2. * sign(PERT, satura(iphase, icv)-OldSatura(iphase, icv))
                end if
            end do
        end do

        CALL calculate_absorption2( packed_state, CV_NONODS, N_IN_PRES, NDIM, SATURA2, TOTELE, CV_NLOC, MAT_NLOC, &
            CV_NDGLN, MAT_NDGLN, U_ABSORB2, PERM%val, visc_phases, IDs_ndgln)

        do ipres = 2, npres

            Spipe => extract_scalar_field( state(1), "Sigma1" )

            do iphase = 1, n_in_pres
                do idim = 1, ndim
                    ! set \sigma for the pipes here
                    LOC = (IPRES-1) * NDIM * N_IN_PRES + (IPHASE-1) * NDIM + IDIM
                    LOC2 = (1-1) * NDIM * N_IN_PRES + (IPHASE-1) * NDIM + IDIM

                    !U_ABSORB( :, LOC, LOC ) = 1000.0
                    U_ABSORB( :, LOC, LOC ) = Spipe%val
                   !U_ABSORB( :, LOC, LOC ) = U_ABSORB( :, LOC2, LOC2 )
                end do
            end do
        end do

        DO ELE = 1, TOTELE
            DO CV_ILOC = 1, CV_NLOC
                IMAT = MAT_NDGLN( ( ELE - 1 ) * MAT_NLOC + CV_ILOC )
                ICV = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )
                DO IPHASE = 1, NPHASE
                    DO JDIM = 1, NDIM
                        DO IDIM = 1, NDIM

                            opt_vel_upwind_coefs_new(IDIM, JDIM, IPHASE, IMAT) = &
                                U_ABSORB( IMAT, IDIM + ( IPHASE - 1 ) * NDIM, JDIM + ( IPHASE - 1 ) * NDIM )

                            if ( iphase <= n_in_pres ) then

                                ! This is the gradient
                                ! Assume d\sigma / dS = 0.0 for the pipes for now
                                opt_vel_upwind_grad_new(IDIM, JDIM, IPHASE, IMAT) = &
                                    (U_ABSORB2( IMAT, IDIM + ( IPHASE - 1 ) * NDIM, JDIM + ( IPHASE - 1 ) * NDIM ) -&
                                    U_ABSORB( IMAT, IDIM + ( IPHASE - 1 ) * NDIM, JDIM + ( IPHASE - 1 ) * NDIM )) &
                                    / ( SATURA2(IPHASE, ICV ) - SATURA(IPHASE, ICV))

                            end if

                        END DO
                    END DO
                END DO
            END DO
        END DO

        deallocate( u_absorb2, satura2, Max_sat )
        ewrite(3,*) 'Leaving calculate_absorption'

        RETURN

    END SUBROUTINE Calculate_AbsorptionTerm



    !sprint_to_do!internal subroutine
    SUBROUTINE calculate_absorption2( packed_state, CV_NONODS, NPHASE, NDIM, SATURA, TOTELE, CV_NLOC, MAT_NLOC, &
        CV_NDGLN, MAT_NDGLN, &
        U_ABSORB, PERM2, visc_phases, IDs_ndgln)
        ! Calculate absorption for momentum eqns
        use matrix_operations
        !    use cv_advection
        implicit none
        type( state_type ), intent( inout ) :: packed_state
        INTEGER, intent( in ) :: CV_NONODS, NPHASE, NDIM, TOTELE, CV_NLOC,MAT_NLOC
        REAL, DIMENSION( :, : ), intent( in ) :: SATURA
        INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN, IDs_ndgln
        REAL, DIMENSION( :, :, : ), intent( inout ) :: U_ABSORB
        REAL, DIMENSION( :, :, : ), intent( in ) :: PERM2
        real, intent(in), dimension(:) :: visc_phases
        ! Local variable
        type (tensor_field), pointer :: RockFluidProp
        real, dimension(:), pointer :: Immobile_fraction, Corey_exponent, Endpoint_relperm
        REAL, PARAMETER :: TOLER = 1.E-10
        INTEGER :: ELE, CV_ILOC, CV_NOD, CV_PHA_NOD, MAT_NOD, JPHA_JDIM, &
            IPHA_IDIM, IDIM, JDIM, IPHASE, id_reg
        REAL, DIMENSION( :, :, :), allocatable :: INV_PERM, PERM

        RockFluidProp=>extract_tensor_field(packed_state,"PackedRockFluidProp")

        ewrite(3,*) 'In calculate_absorption2'
        ALLOCATE( INV_PERM(  size(PERM2,1), size(PERM2,2), size(PERM2,3) ))
        ALLOCATE( PERM( size(PERM2,1), size(PERM2,2) , size(PERM2,3)))
        perm=perm2
        do id_reg = 1, size(perm,3)
            inv_perm( :, :, id_reg)=inverse(perm( :, :, id_reg))
        end do
        U_ABSORB = 0.0
        Loop_NPHASE: DO IPHASE = 1, NPHASE

            Loop_ELE: DO ELE = 1, TOTELE
                !Get properties from packed state
                Immobile_fraction => RockFluidProp%val(1, :, IDs_ndgln(ELE))
                Endpoint_relperm => RockFluidProp%val(2, :, IDs_ndgln(ELE))
                Corey_exponent => RockFluidProp%val(3, :, IDs_ndgln(ELE))
                Loop_CVNLOC: DO CV_ILOC = 1, CV_NLOC

                    MAT_NOD = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + CV_ILOC)
                    CV_NOD = CV_NDGLN(( ELE - 1) * CV_NLOC + CV_ILOC )

                    Loop_DimensionsI: DO IDIM = 1, NDIM

                        Loop_DimensionsJ: DO JDIM = 1, NDIM

                            CV_PHA_NOD = CV_NOD + ( IPHASE - 1 ) * CV_NONODS
                            IPHA_IDIM = ( IPHASE - 1 ) * NDIM + IDIM
                            JPHA_JDIM = ( IPHASE - 1 ) * NDIM + JDIM
                            if (is_porous_media) then
                                select case (NPHASE)
                                    case (1)!No relperm needed, we calculate directly the result
                                        U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = INV_PERM( IDIM, JDIM, ELE ) *&
                                            visc_phases(1) * min(1.0,max(1d-5,SATURA(1,CV_NOD)))
                                    case (2)
                                        CALL relperm_corey_epsilon( U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ), visc_phases, &
                                            INV_PERM( IDIM, JDIM, ELE ), SATURA(iphase,CV_NOD), IPHASE,&
                                            Immobile_fraction, Corey_exponent, Endpoint_relperm)
                                    case (3)!For three phases we use the Stone model. !With predefined order: Water, oil, gas
                                        call relperm_stone(U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ), iphase,&
                                            SATURA(:, CV_NOD), visc_phases, INV_PERM( IDIM, JDIM, ELE),&
                                            Immobile_fraction, Corey_exponent, Endpoint_relperm)
                                    case default
                                        FLAbort("No relative permeability function implemented for more than 3 phases")
                                end select
                            else
                                U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = 0.
                            end if
                        END DO Loop_DimensionsJ

                    END DO Loop_DimensionsI

                END DO Loop_CVNLOC

            END DO Loop_ELE

        END DO Loop_NPHASE

        DEALLOCATE( PERM, INV_PERM )

        ewrite(3,*) 'Leaving calculate_absorption2'

        RETURN

    END SUBROUTINE calculate_absorption2



    SUBROUTINE relperm_corey_epsilon( ABSP, visc_phase, INV_PERM, SAT, IPHASE,&
        Immobile_fraction, Corey_exponent, Endpoint_relperm )
          !This subroutine add a small quantity to the corey function to avoid getting a relperm=0 that may give problems
          !when dividing it to obtain the sigma.
        IMPLICIT NONE
        REAL, intent( inout ) :: ABSP
        REAL, intent( in ) :: SAT, INV_PERM
        INTEGER, intent( in ) :: IPHASE
        real, dimension(:), intent(in) :: visc_phase, Immobile_fraction, Corey_exponent, Endpoint_relperm
        ! Local variables...
        REAL :: KR, aux
        real, parameter :: epsilon = 1d-15!This value should in theory never be used, the real lower limit
        real, parameter :: eps = 1d-5!should be eps ** Corey_exponent
        !Kr_max should only multiply the wetting phase,
        !however as we do not know if it is phase 1 or 2, we let the decision to the user
        !and we multiply both phases by kr_max. By default kr_max= 1

        aux = 1.0 - sum(Immobile_fraction)
        KR = Endpoint_relperm(iphase)*( max( sat - Immobile_fraction(iphase), sat*eps+eps) / aux ) ** Corey_exponent(iphase)
        !Make sure that the relperm is between bounds
        KR = min(max(epsilon, KR),Endpoint_relperm(iphase))!Lower value just to make sure we do not divide by zero.
        ABSP = INV_PERM * (visc_phase(iphase) * max(1d-5, sat)) / KR !The value 1d-5 is only used if the boundaries have values of saturation of zero.
          !Otherwise, the saturation should never be zero, since immobile fraction is always bigger than zero.

        RETURN
    END SUBROUTINE relperm_corey_epsilon


    subroutine relperm_stone(ABSP, iphase, sat, visc, INV_PERM, Immobile_fraction, Corey_exponent, Endpoint_relperm )
        !This subroutine calculates the relative permeability for three phases
        !First phase has to be water, second oil and the third gas
        !We use Stone's model II adapted, and for the two phases we use the Corey model
        !Model explained in: Aziz, K. And Settari, T.:“Petroleum Reservoir Simulation” Applied Science Publishers, London, 30-38, 1979.
        implicit none
        real, intent(inout) :: ABSP
        real, intent(in) :: INV_PERM
        real, dimension(:), intent(in) :: sat, visc, Immobile_fraction, Corey_exponent, Endpoint_relperm
        integer, intent(in) :: iphase
        !Local variables
        real, dimension(3) :: Norm_sat, relperm, KR
        real :: Krow, Krog
        real, parameter :: epsilon = 1d-10

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
        KR(iphase) = min(max(epsilon, relperm(iphase)),Endpoint_relperm(iphase))!Lower value just to make sure we do not divide by zero.
        ABSP = INV_PERM * (VISC(iphase) * max(1d-5,sat(iphase))) / KR(iphase) !The value 1d-5 is only used if the boundaries have values of saturation of zero.
        !Otherwise, the saturation should never be zero, since immobile fraction is always bigger than zero.

    end subroutine relperm_stone

    SUBROUTINE calculate_capillary_pressure( packed_state, Sat_in_FEM,&
        CV_NDGLN, ids_ndgln, totele, cv_nloc)

        ! CAPIL_PRES_OPT is the capillary pressure option for deciding what form it might take.
        ! CAPIL_PRES_COEF( NCAPIL_PRES_COEF, NPHASE, NPHASE ) are the coefficients
        ! Capillary pressure coefs have the dims CAPIL_PRES_COEF( NCAPIL_PRES_COEF, NPHASE,NPHASE )
        ! used to calculate the capillary pressure.

        IMPLICIT NONE
        type(state_type), intent(inout) :: packed_state
        integer, dimension(:), intent(in) :: CV_NDGLN, ids_ndgln
        integer, intent(in) :: totele, cv_nloc
        logical, intent(in) :: Sat_in_FEM
        ! Local Variables
        INTEGER :: IPHASE, JPHASE, nphase, ele, cv_iloc, cv_nod
        !Working pointers
        real, dimension(:,:), pointer :: Satura, CapPressure, Immobile_fraction, Cap_entry_pressure, Cap_exponent
        real, dimension(:), allocatable :: Cont_correction
        !Get from packed_state
        if (Sat_in_FEM) then
            call get_var_from_packed_state(packed_state,FEPhaseVolumeFraction = Satura)
        else
            call get_var_from_packed_state(packed_state,PhaseVolumeFraction = Satura)
        end if
        call get_var_from_packed_state(packed_state,CapPressure = CapPressure, &
            Immobile_fraction = Immobile_fraction, Cap_entry_pressure = Cap_entry_pressure, Cap_exponent = Cap_exponent)
        nphase =size(Satura,1)
        allocate(Cont_correction(size(satura,2)))

        CapPressure = 0.

        DO IPHASE = 1, NPHASE


            if (have_option("/material_phase["//int2str(iphase-1)//&
                "]/multiphase_properties/capillary_pressure/type_Brooks_Corey") ) then

                !Apply Brooks-Corey model
                do jphase = 1, nphase
                    Cont_correction = 0
                    if (jphase /= iphase) then!Don't know how this will work for more than 2 phases
                        do ele = 1, totele
                            do cv_iloc = 1, cv_nloc
                                cv_nod = cv_ndgln((ele-1)*cv_nloc + cv_iloc)
                                CapPressure( jphase, cv_nod ) = CapPressure( jphase, cv_nod ) + &
                                    Get_capPressure(satura(iphase,cv_nod), Cap_entry_pressure(iphase, IDs_ndgln(ele)), &
                                    Cap_exponent(iphase, IDs_ndgln(ele)),Immobile_fraction(:,IDs_ndgln(ele)), iphase)
                                Cont_correction(cv_nod) = Cont_correction(cv_nod) + 1.0
                            end do
                        end do
                        !In continuous formulation nodes are visited more than once, hence we need to average the values added here
                        CapPressure(jphase, :) = CapPressure(jphase, :) / Cont_correction(:)
                    end if
                end do

            end if
        END DO

        deallocate(Cont_correction)
        RETURN
    END SUBROUTINE calculate_capillary_pressure


    pure real function Get_capPressure(sat, Pe, a, Immobile_fraction, iphase)
        !This functions returns the capillary pressure for a certain input saturation
        !There is another function, its derivative in cv-adv-diff called Get_DevCapPressure
        Implicit none
        real, intent(in) :: sat, Pe, a
        real, dimension(:), intent(in) :: Immobile_fraction
        integer, intent(in) :: iphase
        !Local
        real, parameter :: eps = 1d-3 !Small values requires smaller time steps

        Get_capPressure = &
            Pe * min((sat - Immobile_fraction(iphase) + eps) / (1.0 - Immobile_fraction(iphase)), 1.0) ** (-a)
    end function Get_capPressure


    PURE real function Get_DevCapPressure(sat, Pe, a, immobile_fraction, iphase)
        !This functions returns the derivative of the capillary pressure with the saturation
        Implicit none
        integer, intent(in) :: iphase
        real, intent(in) :: sat, Pe, a
        real, dimension(:), intent(in) :: immobile_fraction
        !Local
        real, parameter :: eps = 1d-3
        real :: aux

        aux = (1.0 - immobile_fraction(iphase))

        Get_DevCapPressure = &
            -a * Pe * aux**a * min((sat - immobile_fraction(iphase) + eps), 1.0) ** (-a-1)


    end function Get_DevCapPressure

    subroutine calculate_u_source_cv(state, cv_nonods, ndim, nphase, den, u_source_cv)
        type(state_type), dimension(:), intent(in) :: state
        integer, intent(in) :: cv_nonods, ndim, nphase
        real, dimension(:,:), intent(in) :: den
        real, dimension(:,:,:), intent(inout) :: u_source_cv

        type(vector_field), pointer :: gravity_direction
        real, dimension(ndim) :: g
        logical :: have_gravity, high_order_Ph
        real :: gravity_magnitude
        integer :: idim, iphase, nod, stat

        call get_option( "/physical_parameters/gravity/magnitude", gravity_magnitude, stat )
        have_gravity = ( stat == 0 )

        high_order_Ph = have_option( "/physical_parameters/gravity/hydrostatic_pressure_solver" )

        if( have_gravity .and. .not.high_order_Ph ) then
            gravity_direction => extract_vector_field( state( 1 ), 'GravityDirection' )
            g = node_val( gravity_direction, 1 ) * gravity_magnitude
            u_source_cv = 0.
            do nod = 1, cv_nonods
                do iphase = 1, nphase
                    do idim = 1, ndim
                        u_source_cv( idim, iphase, nod ) = den( iphase, nod ) * g( idim )
                    end do
                end do
            end do

        else
            u_source_cv = 0.
        end if

    end subroutine calculate_u_source_cv

    subroutine calculate_diffusivity(state, ncomp, nphase, ndim, cv_nonods, mat_nonods, &
        mat_nloc, totele, mat_ndgln, ScalarAdvectionField_Diffusion )

        type(state_type), dimension(:), intent(in) :: state
        integer, intent(in) :: ncomp, nphase, ndim, cv_nonods, mat_nonods, mat_nloc, totele
        integer, dimension(:), intent(in) :: mat_ndgln
        real, dimension(:, :, :, :), intent(inout) :: ScalarAdvectionField_Diffusion

        type(scalar_field), pointer :: component
        type(tensor_field), pointer :: diffusivity
        integer, dimension(:), pointer :: element_nodes
        integer :: icomp, iphase, idim, stat, ele
        integer :: iloc,mat_iloc

        ScalarAdvectionField_Diffusion = 0.

        if ( ncomp > 1 ) then

            do icomp = 1, ncomp
                do iphase = 1, nphase

                    component => extract_scalar_field( state(nphase+icomp), 'ComponentMassFractionPhase' // int2str(iphase) )
                    diffusivity => extract_tensor_field( state(nphase+icomp), 'ComponentMassFractionPhase' // int2str(iphase) // 'Diffusivity', stat )

                    if ( stat == 0 ) then

                        do ele = 1, totele

                            element_nodes => ele_nodes( component, ele )

                            do iloc = 1, mat_nloc
                                mat_iloc = mat_ndgln( (ele-1)*mat_nloc + iloc )

                                do idim = 1, ndim

                                    ScalarAdvectionField_Diffusion( mat_iloc, idim, idim, iphase ) = &
                                        ScalarAdvectionField_Diffusion( mat_iloc, idim, idim, iphase ) + &
                                        node_val( component, element_nodes(iloc) ) * node_val( diffusivity, idim, idim, element_nodes(iloc) )

                                end do
                            end do
                        end do
                    end if

                end do
            end do

        else

            do iphase = 1, nphase
                diffusivity => extract_tensor_field( state(iphase), 'TemperatureDiffusivity', stat )

                if ( stat == 0 ) then
                    do idim = 1, ndim
                        ScalarAdvectionField_Diffusion(:, idim, idim, iphase) = node_val( diffusivity, idim, idim, 1 )
                    end do
                end if
            end do

        end if

        return
    end subroutine calculate_diffusivity

    subroutine calculate_viscosity( state, packed_state, ncomp, nphase, ndim, mat_nonods, mat_ndgln, Momentum_Diffusion  )

        type(state_type), dimension(:), intent(in) :: state
        type(state_type), intent(in) :: packed_state

        integer, intent(in) :: ncomp, nphase, ndim, mat_nonods
        integer, dimension(:), intent(in) :: mat_ndgln

        real, dimension( :, :, :, : ), intent(inout) :: Momentum_Diffusion
        character( len = option_path_len ) :: option_path_python, buffer
        type(tensor_field), pointer :: t_field
        integer :: iphase, icomp, stat, mat_nod, cv_nloc, ele

        type(scalar_field), pointer :: component
        logical :: linearise_viscosity, python_diagnostic_field
        real, dimension( : ), allocatable :: component_tmp
        real, dimension( :, :, : ), allocatable :: mu_tmp

        character( len = python_func_len ) :: pycode
        real :: dt, current_time
        integer :: iloc



        if ( have_option( '/physical_parameters/mobility' ) .or. is_porous_media ) then

            ! if solving for porous media and mobility is calculated
            ! through the viscosity ratio this code will fail
            momentum_diffusion=0.

        else

            momentum_diffusion=0.

            t_field => extract_tensor_field( state( 1 ), 'Viscosity', stat )
            if ( stat == 0 ) then

                cv_nloc = ele_loc( t_field, ele )
                linearise_viscosity = have_option( '/material_phase[0]/linearise_viscosity' )
                allocate( component_tmp( cv_nloc ), mu_tmp( ndim, ndim, cv_nloc ) )


                if ( ncomp > 1 ) then


                    do icomp = 1, ncomp
                        do iphase = 1, nphase

                            component => extract_scalar_field( state(nphase + icomp), 'ComponentMassFractionPhase' // int2str(iphase) )

                            python_diagnostic_field = &
                                have_option( '/material_phase[' // int2str(nphase + icomp - 1 ) //  &
                                ']/scalar_field::ComponentMassFractionPhase' // int2str(iphase) // &
                                '/prognostic/tensor_field::Viscosity/diagnostic'    )

                            if ( python_diagnostic_field ) then

#ifdef HAVE_NUMPY
                        ewrite(3,*) "Have both NumPy and a python viscosity..."
#else
                                FLAbort("Python eos requires NumPy, which cannot be located.")
#endif

                                option_path_python = trim( '/material_phase[' // int2str( nphase + icomp - 1 ) // &
                                    ']/scalar_field::ComponentMassFractionPhase' // int2str( iphase ) // &
                                    '/prognostic/tensor_field::Viscosity/diagnostic' )

                                t_field => extract_tensor_field( packed_state, 'Dummy' )
                                call zero( t_field )

                                call python_reset()
                                call python_add_state( packed_state )

                                call python_run_string("field = state.tensor_fields['Dummy']")
                                call get_option("/timestepping/current_time", current_time)
                                write(buffer,*) current_time
                                call python_run_string("time="//trim(buffer))
                                call get_option("/timestepping/timestep", dt)
                                write(buffer,*) dt
                                call python_run_string("dt="//trim(buffer))

                                ! Get the code
                                call get_option( trim( option_path_python ) // '/algorithm', pycode )

                                ! Run the code
                                call python_run_string( trim( pycode ) )

                            else
                                t_field => extract_tensor_field( state( nphase + icomp ), 'Viscosity' )
                            end if

                            ewrite(3,*) 'Component, Phase, Visc_min_max', icomp, iphase, minval( t_field%val ), maxval( t_field%val )

                            do ele = 1, ele_count( t_field )

                                component_tmp = ele_val( component, ele )
                                mu_tmp = ele_val( t_field, ele )

                                do iloc = 1, cv_nloc
                                    mu_tmp( :, :, iloc ) = mu_tmp( :, :, iloc ) * component_tmp( iloc )
                                end do

                                if ( linearise_viscosity ) then
                                    mu_tmp( :, :, 2 ) = 0.5 * ( mu_tmp( :, :, 1 ) + mu_tmp( :, :, 3 ) )
                                    mu_tmp( :, :, 4 ) = 0.5 * ( mu_tmp( :, :, 1 ) + mu_tmp( :, :, 6 ) )
                                    mu_tmp( :, :, 5 ) = 0.5 * ( mu_tmp( :, :, 3 ) + mu_tmp( :, :, 6 ) )

                                    if ( cv_nloc == 10 ) then
                                        mu_tmp( :, :, 7 ) = 0.5 * ( mu_tmp( :, :, 1 ) + mu_tmp( :, :, 10 ) )
                                        mu_tmp( :, :, 8 ) = 0.5 * ( mu_tmp( :, :, 3 ) + mu_tmp( :, :, 10 ) )
                                        mu_tmp( :, :, 9 ) = 0.5 * ( mu_tmp( :, :, 6 ) + mu_tmp( :, :, 10 ) )
                                    end if
                                end if


                                do iloc = 1, cv_nloc
                                    mat_nod = mat_ndgln( (ele-1)*cv_nloc + iloc )
                                    momentum_diffusion( mat_nod, :, :, iphase ) = momentum_diffusion( mat_nod, :, :, iphase ) + mu_tmp( :, :, iloc )
                                end do
                            end do

                        end do
                    end do

                else

                    do iphase = 1, nphase


                        python_diagnostic_field = &
                            have_option( '/material_phase[' // int2str(iphase - 1 ) //  &
                            ']/vector_field::Velocity/prognostic/tensor_field::Viscosity/diagnostic' )


                        if ( python_diagnostic_field ) then

#ifdef HAVE_NUMPY
                     ewrite(3,*) "Have both NumPy and a python viscosity..."
#else
                            FLAbort("Python eos requires NumPy, which cannot be located.")
#endif

                            option_path_python = trim( '/material_phase[' // int2str( iphase - 1 ) // &
                                ']/vector_field::Velocity' // &
                                '/prognostic/tensor_field::Viscosity/diagnostic' )

                            t_field => extract_tensor_field( packed_state, 'Dummy' )
                            call zero( t_field )

                            call python_reset()
                            call python_add_state( packed_state )

                            call python_run_string("field = state.tensor_fields['Dummy']")
                            call get_option("/timestepping/current_time", current_time)
                            write(buffer,*) current_time
                            call python_run_string("time="//trim(buffer))
                            call get_option("/timestepping/timestep", dt)
                            write(buffer,*) dt
                            call python_run_string("dt="//trim(buffer))

                            ! Get the code
                            call get_option( trim( option_path_python ) // '/algorithm', pycode )

                            ! Run the code
                            call python_run_string( trim( pycode ) )

                        else

                            t_field => extract_tensor_field( state( iphase ), 'Viscosity', stat )

                        end if


                        do ele = 1, ele_count( t_field )

                            mu_tmp = ele_val( t_field, ele )

                            if ( linearise_viscosity ) then
                                mu_tmp( :, :, 2 ) = 0.5 * ( mu_tmp( :, :, 1 ) + mu_tmp( :, :, 3 ) )
                                mu_tmp( :, :, 4 ) = 0.5 * ( mu_tmp( :, :, 1 ) + mu_tmp( :, :, 6 ) )
                                mu_tmp( :, :, 5 ) = 0.5 * ( mu_tmp( :, :, 3 ) + mu_tmp( :, :, 6 ) )

                                if ( cv_nloc == 10 ) then
                                    mu_tmp( :, :, 7 ) = 0.5 * ( mu_tmp( :, :, 1 ) + mu_tmp( :, :, 10 ) )
                                    mu_tmp( :, :, 8 ) = 0.5 * ( mu_tmp( :, :, 3 ) + mu_tmp( :, :, 10 ) )
                                    mu_tmp( :, :, 9 ) = 0.5 * ( mu_tmp( :, :, 6 ) + mu_tmp( :, :, 10 ) )
                                end if
                            end if

                            do iloc = 1, cv_nloc
                                mat_nod = mat_ndgln( (ele-1)*cv_nloc + iloc )
                                momentum_diffusion( mat_nod, :, :, iphase ) = mu_tmp( :, :, iloc )
                            end do
                        end do

                    end do

                end if

                deallocate( component_tmp, mu_tmp )

            end if

        end if

        return
    end subroutine calculate_viscosity

    subroutine calculate_SUF_SIG_DIAGTEN_BC( packed_state, suf_sig_diagten_bc, totele, stotel, cv_nloc, &
        cv_snloc, nphase, nphase_total, ndim, nface, mat_nonods, cv_nonods, x_nloc, ncolele, cv_ele_type, &
        finele, colele, cv_ndgln, cv_sndgln, x_ndgln, mat_ndgln, material_absorption, &
        state, x_nonods, IDs_ndgln )

        implicit none
        type( state_type ), intent( inout ) :: packed_state
        integer, intent(in) :: totele, stotel, cv_nloc, cv_snloc, nphase, nphase_total, ndim, nface, &
            mat_nonods, cv_nonods, x_nloc, ncolele, cv_ele_type, x_nonods
        integer, dimension( : ), intent( in ) :: finele, colele, cv_ndgln, cv_sndgln, x_ndgln, mat_ndgln, IDs_ndgln
        real, dimension( :, :, : ), intent( inout ) :: material_absorption
        type(state_type), dimension( : ) :: state

        real, dimension( stotel * cv_snloc * nphase_total, ndim ), intent( inout ) :: suf_sig_diagten_bc

        ! local variables
        type(tensor_field), pointer :: viscosity_ph, RockFluidProp
        real, dimension(:), pointer :: Immobile_fraction, Corey_exponent, Endpoint_relperm
        real, dimension(nphase) :: visc_phases
        integer :: iphase, ele, sele, cv_siloc, cv_snodi, cv_snodi_ipha, iface, s, e, &
            ele2, sele2, cv_iloc, idim, jdim, i, mat_nod, cv_nodi
        real :: mobility, satura_bc
        real, dimension( ndim, ndim ) :: inv_perm, sigma_out, sigma_in, mat, mat_inv
        integer, dimension( nface, totele) :: face_ele
        integer, dimension( mat_nonods*nphase ) :: idone
        integer, dimension(nface, cv_snloc ) :: cv_sloclist
        integer, dimension( cv_snloc ) :: cv_sloc2loc

        integer, dimension( :, :, : ),  allocatable :: wic_u_bc, wic_vol_bc

        integer, parameter :: WIC_BC_DIRICHLET = 1

        !!$ for the pressure b.c. and compact_overlapping method
        !!$ make the material property change just inside the domain else on the surface only
        logical, parameter :: mat_change_inside = .false.

        !!$ if mat_perm_bc_dg use the method that is used for DG between the elements
        logical, parameter :: mat_perm_bc_dg = .true.

        type(tensor_field), pointer :: velocity, volfrac, perm
        type(tensor_field) :: velocity_BCs, volfrac_BCs

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

        if( nphase == 1 ) then
           viscosity_ph => extract_tensor_field( state( 1 ), 'Viscosity' )
           visc_phases(1) = viscosity_ph%val( 1, 1, 1 )
           mobility = visc_phases(1)
        elseif( have_option( '/physical_parameters/mobility' ) )then
           call get_option( '/physical_parameters/mobility', mobility )
           visc_phases(1) = 1
           visc_phases(2) = mobility
        elseif( have_option( '/material_phase[1]/vector_field::Velocity/prognostic/tensor_field::Viscosity' // &
             '/prescribed/value::WholeMesh/isotropic' ) ) then
           DO IPHASE = 1, NPHASE ! Get viscosity for all the phases
              viscosity_ph => extract_tensor_field( state( iphase ), 'Viscosity' )
              visc_phases(iphase) = viscosity_ph%val( 1, 1, 1 ) ! So far we only consider scalar viscosity
           end do
           mobility = visc_phases(2) / visc_phases(1)
        end if

        suf_sig_diagten_bc = 1.

        idone=0
        call determin_sloclist( cv_sloclist, cv_nloc, cv_snloc, nface,  &
            ndim, cv_ele_type )

        face_ele = 0
        call calc_face_ele( face_ele, totele, stotel, nface, &
            ncolele, finele, colele, cv_nloc, cv_snloc, cv_nonods, cv_ndgln, cv_sndgln, &
            cv_sloclist, x_nloc, x_ndgln )


        do iphase = 1, nphase

            s = ( iphase - 1 ) * ndim + 1
            e = iphase * ndim

            do ele = 1, totele

                !Get properties from packed state
                Immobile_fraction => RockFluidProp%val(1, :, IDs_ndgln(ELE))
                Endpoint_relperm => RockFluidProp%val(2, :, IDs_ndgln(ELE))
                Corey_exponent => RockFluidProp%val(3, :, IDs_ndgln(ELE))
                inv_perm = inverse( perm%val(:, :, ele) )

                do iface = 1, nface

                    ele2  = face_ele( iface, ele )
                    sele2 = max( 0, -ele2 )
                    sele  = sele2
                    if ( sele > 0 ) then
                        if ( wic_u_bc(1,iphase,sele) /= WIC_BC_DIRICHLET .and. &
                            wic_vol_bc(1,iphase,sele) == WIC_BC_DIRICHLET ) then

                            cv_sloc2loc( : ) = cv_sloclist( iface, : )
                            do cv_siloc = 1, cv_snloc

                                cv_iloc = cv_sloc2loc( cv_siloc )
                                cv_snodi = ( sele - 1 ) * cv_snloc + cv_siloc
                                cv_nodi = cv_sndgln(cv_snodi)
                                cv_snodi_ipha = cv_snodi + ( iphase - 1 ) * stotel * cv_snloc
                                mat_nod = mat_ndgln( (ele-1)*cv_nloc + cv_iloc  )
                                ! this is the boundary condition
                                satura_bc = volfrac_BCs%val(1,iphase,cv_snodi)

                                do idim = 1, ndim
                                    do jdim = 1, ndim
                                        select case ( nphase )
                                            case (1) ! No relperm needed, we calculate directly the result
                                                sigma_out( idim, jdim ) = inv_perm( idim, jdim ) *&
                                                    visc_phases(1) * min(1.0,max(1d-5,satura_bc))
                                            case (2)
                                                call relperm_corey_epsilon( sigma_out( idim, jdim ), visc_phases, &
                                                    inv_perm( idim, jdim ), satura_bc, IPHASE,&
                                                    Immobile_fraction, Corey_exponent, Endpoint_relperm)
                                            case (3) ! For three phases we use the Stone model. !With predefined order: Water, oil, gas
                                                call relperm_stone(sigma_out( idim, jdim ), iphase,&
                                                    volfrac_BCs%val(1,:,cv_snodi), visc_phases, inv_perm( idim, jdim ),&
                                                    Immobile_fraction, Corey_exponent, Endpoint_relperm)
                                            case default
                                                FLAbort("No relative permeability function implemented for more than 3 phases")
                                        end select
                                    end do
                                end do

                                if ( mat_perm_bc_dg ) then
                                    ! if mat_perm_bc_dg use the method that is used for DG between the elements.
                                    sigma_in=0.0
                                    sigma_in = material_absorption( mat_nod, s : e, s : e )
                                    mat = sigma_out  +  matmul(  sigma_in,  matmul( inverse( sigma_out ), sigma_in ) )
                                    mat_inv = matmul( inverse( sigma_in+sigma_out ), mat )
                                    suf_sig_diagten_bc( cv_snodi_ipha, 1 : ndim ) = (/ (mat_inv(i, i), i = 1, ndim) /)
                                   !suf_sig_diagten_bc( cv_snodi_ipha, 1 : ndim ) = 1.
                                else
                                    mat = matmul( sigma_out, inverse( material_absorption( mat_nod, s : e, s : e ) ) )
                                    mat_inv = inverse( mat )
                                    suf_sig_diagten_bc( cv_snodi_ipha, 1 : ndim ) = (/ (mat_inv(i, i), i = 1, ndim) /)
                                end if

                                if ( mat_change_inside ) then
                                    suf_sig_diagten_bc( cv_snodi_ipha, 1 : ndim ) = 1.

                                    if ( idone( mat_nod+(iphase-1)*mat_nonods ) == 0 ) then
                                        material_absorption( mat_nod, s : e, s : e  ) &
                                            = matmul( mat, material_absorption( mat_nod, s : e, s : e ) )
                                        idone( mat_nod+(iphase-1)*mat_nonods ) = 1
                                    end if
                                end if

                            end do

                        end if
                    end if

                end do
            end do

        end do

        call deallocate(velocity_BCs)
        call deallocate(volfrac_BCs)
        deallocate(wic_u_bc, wic_vol_bc)

        return
    end subroutine calculate_SUF_SIG_DIAGTEN_BC


    subroutine calculate_u_abs_stab( Material_Absorption_Stab, Material_Absorption, &
        opt_vel_upwind_coefs_new, nphase, ndim, totele, cv_nloc, mat_nloc, mat_nonods, mat_ndgln )

        implicit none

        real, dimension( :, :, : ), intent( inout ) :: Material_Absorption_Stab
        real, dimension( :, :, : ), intent( in ) :: Material_Absorption
        real, dimension( :, :, :, : ), intent( in ) :: opt_vel_upwind_coefs_new
        integer, intent( in ) :: nphase, ndim, totele, cv_nloc, mat_nloc, mat_nonods
        integer, dimension( : ), intent( in ) :: mat_ndgln

        logical :: use_mat_stab_stab
        integer :: apply_dim, idim, jdim, ipha_idim, iphase, ele, cv_iloc, imat
        real :: factor

        Material_Absorption_Stab = 0.

        use_mat_stab_stab = .false.

        if ( use_mat_stab_stab ) then

            apply_dim = 2

            if ( .true. ) then

                factor = 100.

                do iphase = 1, nphase
                    do idim = 1, ndim
                        if ( idim == apply_dim ) then
                            ipha_idim = ( iphase - 1 ) * ndim + idim
                            Material_Absorption_Stab( :, ipha_idim, ipha_idim ) = &
                                ( factor / 10. )**2 * Material_Absorption( :, ipha_idim, ipha_idim )
                        end if
                    end do
                end do

            else

                do ele = 1, totele
                    do cv_iloc = 1, cv_nloc
                        imat = mat_ndgln( ( ele - 1 ) * mat_nloc + cv_iloc )
                        do iphase = 1, nphase
                            do idim = 1, ndim
                                do jdim = 1, ndim
                                    if ( idim == apply_dim .and. jdim == apply_dim ) then
                                        ipha_idim = ( iphase - 1 ) * ndim + idim
                                        Material_Absorption_Stab( imat, ipha_idim, ipha_idim ) = &
                                            abs( opt_vel_upwind_coefs_new( idim, jdim, iphase, imat ) )
                                    end if
                                end do
                            end do
                        end do
                    end do
                end do

            end if
        end if ! use_mat_stab_stab

        return
    end subroutine calculate_u_abs_stab

    subroutine update_velocity_absorption( states, ndim, nphase, mat_nonods,velocity_absorption )

        implicit none

        integer, intent( in ) :: ndim, nphase, mat_nonods
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
            if ( have_absorption ) then
                absorption => extract_vector_field( states( iphase ), 'VelocityAbsorption' )
                do idim = 1, ndim
                    velocity_absorption( :, idim + (iphase-1)*ndim, idim + (iphase-1)*ndim ) =  &
                        absorption % val( idim, : )
                end do
            else
                do idim = 1, ndim
                    velocity_absorption( :, idim + (iphase-1)*ndim, idim + (iphase-1)*ndim ) = 0.0
                end do
            end if
        end do

        return
    end subroutine update_velocity_absorption


    subroutine update_velocity_absorption_coriolis( states, ndim, nphase, sigma )

      implicit none

      integer, intent( in ) :: ndim, nphase
      type( state_type ), dimension( : ), intent( in ) :: states
      real, dimension( :, :, : ), intent( inout ) :: sigma

      type( scalar_field ), pointer :: f
      integer :: iphase, stat, idx1, idx2

      do iphase = 1, nphase
         f => extract_scalar_field( states( iphase ), 'f', stat )
         if ( stat == 0 ) then
            idx1 = 1 + (iphase-1)*ndim ;  idx2 = 2 + (iphase-1)*ndim
            sigma( :, idx1, idx2 ) = sigma( :, idx1, idx2 ) - f % val
            sigma( :, idx2, idx1 ) = sigma( :, idx2, idx1 ) + f % val
         end if
      end do

      return
    end subroutine update_velocity_absorption_coriolis



    subroutine update_velocity_source( states, ndim, nphase, u_nonods, velocity_u_source )

        implicit none

        integer, intent( in ) :: ndim, nphase, u_nonods
        type( state_type ), dimension( : ), intent( in ) :: states
        real, dimension( :, :, : ), intent( inout ) :: velocity_u_source

        type( vector_field ), pointer :: source
        integer :: iphase, idim
        logical :: have_source
        character( len = option_path_len ) :: option_path

        velocity_u_source = 0.

        do iphase = 1, nphase
            have_source = .false.
            option_path = '/material_phase[' // int2str( iphase - 1 ) // ']/vector_field::Velocity' // &
                '/prognostic/vector_field::Source/prescribed'
            have_source = have_option( trim(option_path) )
            if ( have_source ) then
                source => extract_vector_field( states( iphase ), 'VelocitySource' )
                do idim = 1, ndim
                    velocity_u_source( idim, iphase, : ) =  source % val( idim, : )
                end do
            else
                do idim = 1, ndim
                    velocity_u_source( idim, iphase, : ) = 0.0
                end do
            end if
        end do

        return
    end subroutine update_velocity_source


    !sprint_to_do!do this properly
    subroutine extract_scalar_from_diamond(state, storage_state, field_values, path, StorName, indx, iphase, nphase)
        !Gets a scalar field directly from Diamond
        !Path have to end in /prescribed/value
        !Indx is for the cashing
        !If no phases, then pass iphase = nphase = 1
        !NOTE: This was initially done for capillary pressure with regions
        implicit none
        type(state_type), dimension(:), intent(inout) :: state
        type(state_type), intent(inout) :: storage_state
        real, dimension(:), pointer, intent(inout) :: field_values
        character(len=*), intent(in) :: path, StorName
        integer, intent(inout) :: indx
        integer, intent(in) :: iphase, nphase
        !Working pointers
        type (scalar_field), pointer :: Sfield
        type(vector_field), pointer :: position
        type(scalar_field), target :: targ_Store
        type(mesh_type), pointer :: fl_mesh
        type(mesh_type) :: Auxmesh
        integer :: siz
        if (indx<=0) then!Everything needs to be calculated
            if (has_scalar_field(storage_state, StorName)) then
                !If we are recalculating due to a mesh modification then
                !we return to the original situation
                call remove_scalar_field(storage_state, StorName)
            end if


            !By default I use the Pressure mesh (Number 1)
            Sfield => extract_scalar_field(state(1),1)
            position => get_external_coordinate_field(state(1), Sfield%mesh)
            fl_mesh => extract_mesh( storage_state, "FakeMesh" )
            Auxmesh = fl_mesh
            !The number of nodes I want does not coincide
            Auxmesh%nodes = size(Sfield%val,1) * nphase
            call allocate (targ_Store, Auxmesh, StorName)

            !            call allocate(targ_Store, Sfield%mesh)
            call initialise_field_over_regions(targ_Store, path, position)
            !Now we insert them in state and store the indexes
            call insert(storage_state, targ_Store, StorName)
            call deallocate (targ_Store)
            indx = size(storage_state%scalar_fields)
        end if
        !Get the data
        siz = size(storage_state%scalar_fields(abs(indx))%ptr%val(:),1)/nphase
        field_values => storage_state%scalar_fields(abs(indx))%ptr%val((iphase-1)*siz + 1: siz * iphase )


    end subroutine extract_scalar_from_diamond

    !sprint_to_do!delete before the sprint is over
    subroutine boiling( states, packed_state, cv_nonods, mat_nonods, nphase, ndim, &
        ScalarField_Source, velocity_absorption, temperature_absorption )
        implicit none

        type( state_type ), dimension(:), intent( inout ) :: states
        type( state_type ), intent( in ) :: packed_state
        integer, intent( in ) :: cv_nonods, mat_nonods, nphase, ndim
        real, dimension( :, : ), intent( inout ) :: ScalarField_Source
        real, dimension( :, :, : ), intent( inout ) :: velocity_absorption, temperature_absorption

        type( tensor_field ), pointer :: temperature, temperature_source
        real, dimension( :, : ), allocatable :: A
        real, dimension( : ), allocatable :: S_lg_l, S_lg_g, S_ls_l, S_gs_g, &
            T_sat, Svap_l, Svap_g, Gamma_l, Gamma_g, h_l, h_g, &
            St_gl, St_sl, St_sg
        integer :: iphase, jphase, idim
        real, parameter :: Le0=2375.7e3, Cp_l = 4200.0, Cp_g = 1996.0, HeatSource=17.1e+3 / (0.082 * 0.095 * 0.25)

        ewrite(3,*) 'inside boiling routine'

        ScalarField_Source=0.0
        velocity_absorption=0.0 ; temperature_absorption=0.0

        allocate( S_lg_l(mat_nonods), S_lg_g(mat_nonods), S_ls_l(mat_nonods), S_gs_g(mat_nonods), A(nphase,mat_nonods) )
        allocate( T_sat(cv_nonods), Svap_l(cv_nonods), Svap_g(cv_nonods), &
            &    Gamma_l(cv_nonods), Gamma_g(cv_nonods), h_l(cv_nonods), h_g(cv_nonods), &
            &    St_gl(cv_nonods), St_sl(cv_nonods), St_sg(cv_nonods) )


        call calculate_boiling_drag( states, packed_state, S_lg_l, S_lg_g, S_ls_l, S_gs_g, A )

        ! Momentum absorption

        ! open the boiling test for two phases-gas and liquid
        if (have_option('/boiling')) then
            S_ls_l=0.0
            S_gs_g=0.0
        end if



        iphase=1 ; jphase=1
        do idim = 1, ndim
            velocity_absorption( :, idim + (iphase-1)*ndim, idim + (jphase-1)*ndim ) = S_lg_l + S_ls_l
        end do
        iphase=1 ; jphase=2
        do idim = 1, ndim
            velocity_absorption( :, idim + (iphase-1)*ndim, idim + (jphase-1)*ndim ) = -S_lg_l
        end do
        iphase=1 ; jphase=3
        do idim = 1, ndim
            velocity_absorption( :, idim + (iphase-1)*ndim, idim + (jphase-1)*ndim ) = -S_ls_l
        end do

        iphase=2 ; jphase=1
        do idim = 1, ndim
            velocity_absorption( :, idim + (iphase-1)*ndim, idim + (jphase-1)*ndim ) = -S_lg_g
        end do
        iphase=2 ; jphase=2
        do idim = 1, ndim
            velocity_absorption( :, idim + (iphase-1)*ndim, idim + (jphase-1)*ndim ) = S_lg_g + S_gs_g
        end do
        iphase=2 ; jphase=3
        do idim = 1, ndim
            velocity_absorption( :, idim + (iphase-1)*ndim, idim + (jphase-1)*ndim ) = -S_gs_g
        end do

        iphase=3 ; jphase=1
        do idim = 1, ndim
            velocity_absorption( :, idim + (iphase-1)*ndim, idim + (jphase-1)*ndim ) = 0.0
        end do
        iphase=3 ; jphase=2
        do idim = 1, ndim
            velocity_absorption( :, idim + (iphase-1)*ndim, idim + (jphase-1)*ndim ) = 0.0
        end do
        iphase=3 ; jphase=3
        do idim = 1, ndim
            velocity_absorption( :, idim + (iphase-1)*ndim, idim + (jphase-1)*ndim ) = 1.0e+15
        end do


        call calculate_boiling_variables( states, packed_state, ndim, nphase, T_sat, Svap_l, Svap_g, &
            Gamma_l, Gamma_g, h_l, h_g, St_gl, St_sl, St_sg )


        ! Temperature absorption

        ! open the boiling test for two phases-gas and liquid
        if (have_option('/boiling')) then
            St_sl=0.0
            St_sg=0.0
        end if


        iphase=1 ; jphase=1
        temperature_absorption( iphase, jphase, : ) =  St_gl + St_sl + Svap_l + Cp_l*Gamma_l
        iphase=1 ; jphase=2
        temperature_absorption( iphase, jphase, : ) = -St_gl
        iphase=1 ; jphase=3
        temperature_absorption( iphase, jphase, : ) = -St_sl

        iphase=2 ; jphase=1
        temperature_absorption( iphase, jphase, : ) = -St_gl
        iphase=2 ; jphase=2
        temperature_absorption( iphase, jphase, : ) =  St_gl + St_sg + Svap_g + Cp_g*Gamma_g + 1.0e+6
        iphase=2 ; jphase=3
        temperature_absorption( iphase, jphase, : ) = -St_sg

        iphase=3 ; jphase=1
        temperature_absorption( iphase, jphase, : ) = -St_sl
        iphase=3 ; jphase=2
        temperature_absorption( iphase, jphase, : ) = -St_sg
        iphase=3 ; jphase=3
        temperature_absorption( iphase, jphase, : ) =  St_sl + St_sg


        ! Temperature source

        temperature => extract_tensor_field( packed_state, "PackedTemperature" )
        temperature_source => extract_tensor_field( packed_state, "PackedTemperatureSource" )

        iphase=1
        temperature_source%val( 1, iphase, : ) = Svap_l*T_sat + Gamma_l*h_l + Gamma_l*Le0 + HeatSource

        iphase=2
        temperature_source%val( 1, iphase, : ) = Svap_g*T_sat + Gamma_g*h_g + 1.0e+6 * temperature%val( 1, iphase, : )


        ! Mass source

        iphase=1
        ScalarField_Source( iphase, : ) = Gamma_l

        iphase=2
        ScalarField_Source( iphase, : ) = Gamma_g


        ! deallocate
        deallocate( S_lg_l, S_lg_g, S_ls_l, S_gs_g, A )
        deallocate( T_sat, Svap_l, Svap_g, Gamma_l, Gamma_g, h_l, h_g, St_gl, St_sl, St_sg )

        ewrite(3,*) 'leaving boiling routine'

        return
    end subroutine boiling



    !sprint_to_do!delete before the sprint is over
    subroutine calculate_boiling_drag( states, packed_state, S_lg_l, S_lg_g, S_ls_l, S_gs_g, A )
        implicit none

        type( state_type ), dimension( : ), intent( inout ) :: states
        type( state_type ), intent( in ) :: packed_state
        real, dimension( : ), intent( inout ) :: S_lg_l, S_lg_g, S_ls_l, S_gs_g
        real, dimension( :, : ), intent( inout ) :: A

        type( tensor_field ), pointer :: pressure
        type( tensor_field ), pointer :: density, velocity, volume_fraction

        integer, dimension( : ), pointer :: mat_ndgln, cv_ndgln, u_ndgln

        real :: rho_l, rho_g, u_l, u_g, u_s, u_gs, u_ls, u_gl, a_l, a_g, a_s, &
            a_sg, a_gs, a_sl, a_ls, a_lg, a_gl, d_b, Re_gl, Re_lg, CD

        real, dimension( : ), allocatable :: ul, ug, us

        integer :: ele, totele, mat_iloc, mat_nloc, u_nloc, mat_inod, cv_inod, &
            u_inod, u_inod1, u_inod2, ndim

        real, parameter :: &
            mu_l = 3.0e-4, mu_g = 1.0e-5, &
            d_p = 0.005

        pressure => extract_tensor_field( packed_state, "PackedCVPressure" )
        density => extract_tensor_field( packed_state, "PackedDensity" )
        velocity => extract_tensor_field( packed_state, "PackedNonlinearVelocity" )

        volume_fraction => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )

        mat_ndgln => get_ndglno( extract_mesh( packed_state, "PressureMesh_Discontinuous" ) )
        cv_ndgln => get_ndglno( extract_mesh( packed_state, "PressureMesh" ) )
        u_ndgln => get_ndglno( extract_mesh( packed_state, "VelocityMesh" ) )

        totele = ele_count( pressure )
        mat_nloc = ele_loc( pressure, 1 )
        u_nloc = ele_loc( velocity, 1 )
        ndim = size(velocity%val,1)

        S_lg_l=0.0 ; S_lg_g=0.0 ; S_ls_l=0.0 ; S_gs_g=0.0 ; A=0.0

        allocate( ul(ndim), ug(ndim), us(ndim) ) ; ul=0.0 ; ug=0.0 ; us=0.0

        do ele = 1, totele

            do mat_iloc = 1, mat_nloc

                mat_inod = mat_ndgln( ( ele - 1 ) * mat_nloc + mat_iloc )
                cv_inod = cv_ndgln( ( ele - 1 ) * mat_nloc + mat_iloc )

                rho_l = density%val(1,1,cv_inod) ; rho_g = density%val(1,2,cv_inod)

                if ( mat_iloc==1  ) then
                    u_inod = u_ndgln( ( ele - 1 ) * u_nloc + 1 )
                    ul = velocity%val(:,1,u_inod) ; ug = velocity%val(:,2,u_inod) ; us = velocity%val(:,3,u_inod)
                end if
                if ( mat_iloc==3 ) then
                    u_inod = u_ndgln( ( ele - 1 ) * u_nloc + 2 )
                    ul = velocity%val(:,1,u_inod) ; ug = velocity%val(:,2,u_inod) ; us = velocity%val(:,3,u_inod)
                end if
                if ( mat_iloc==6 ) then
                    u_inod = u_ndgln( ( ele - 1 ) * u_nloc + 3 )
                    ul = velocity%val(:,1,u_inod) ; ug = velocity%val(:,2,u_inod) ; us = velocity%val(:,3,u_inod)
                end if
                if ( mat_iloc==10 ) then
                    u_inod = u_ndgln( ( ele - 1 ) * u_nloc + 4 )
                    ul = velocity%val(:,1,u_inod) ; ug = velocity%val(:,2,u_inod) ; us = velocity%val(:,3,u_inod)
                end if
                u_inod1=-666 ; u_inod2=-666
                if ( mat_iloc==2 ) then
                    u_inod1 = u_ndgln( ( ele - 1 ) * u_nloc + 1 )
                    u_inod2 = u_ndgln( ( ele - 1 ) * u_nloc + 2 )
                end if
                if ( mat_iloc==4 ) then
                    u_inod1 = u_ndgln( ( ele - 1 ) * u_nloc + 1 )
                    u_inod2 = u_ndgln( ( ele - 1 ) * u_nloc + 3 )
                end if
                if ( mat_iloc==5 ) then
                    u_inod1 = u_ndgln( ( ele - 1 ) * u_nloc + 2 )
                    u_inod2 = u_ndgln( ( ele - 1 ) * u_nloc + 3 )
                end if
                if ( mat_iloc==7 ) then
                    u_inod1 = u_ndgln( ( ele - 1 ) * u_nloc + 1 )
                    u_inod2 = u_ndgln( ( ele - 1 ) * u_nloc + 4 )
                end if
                if ( mat_iloc==8 ) then
                    u_inod1 = u_ndgln( ( ele - 1 ) * u_nloc + 2 )
                    u_inod2 = u_ndgln( ( ele - 1 ) * u_nloc + 4 )
                end if
                if ( mat_iloc==9 ) then
                    u_inod1 = u_ndgln( ( ele - 1 ) * u_nloc + 3 )
                    u_inod2 = u_ndgln( ( ele - 1 ) * u_nloc + 4 )
                end if
                if ( u_inod1>0 ) then
                    ul = (velocity%val(:,1,u_inod1)+velocity%val(:,1,u_inod2))/2.0
                    ug = (velocity%val(:,2,u_inod1)+velocity%val(:,2,u_inod2))/2.0
                    us = (velocity%val(:,3,u_inod1)+velocity%val(:,3,u_inod2))/2.0
                end if

                u_l = sqrt( sum( ul**2 ) )
                u_g = sqrt( sum( ug**2 ) )
                u_s = sqrt( sum( us**2 ) )

                u_gs=abs(u_g-u_s) ; u_ls=abs(u_l-u_s) ; u_gl=abs(u_g-u_l)

                a_l = volume_fraction%val(1,1,cv_inod)
                a_g = volume_fraction%val(1,2,cv_inod)
                a_s = volume_fraction%val(1,3,cv_inod)

                A(1,mat_inod)=a_l ;  A(2,mat_inod)=a_g ; A(3,mat_inod)=a_s

                a_sg=a_s/(a_s+a_g) ; a_gs=1.0-a_sg
                a_sl=a_s/(a_s+a_l) ; a_ls=1.0-a_sl
                a_lg=a_l/(a_l+a_g) ; a_gl=1.0-a_lg

                d_b = 5.0*0.06/max((rho_l*u_gl**2),1.0e-5) ; d_b=min(0.5*d_p,max(1.0e-7,d_b))
                !d_b = 1.0*d_p

                Re_gl = max(rho_l*u_gl*d_b/mu_l,1.0e-5) ; Re_lg=Re_gl
                if(a_lg*Re_lg<1.0e3)then
                    CD=(24.0/(max(a_lg,1.0e-5)*Re_lg))*(1.0+0.15*(a_lg*Re_lg)**0.687)
                else
                    CD=0.44
                end if

                !S_gs_g(mat_inod) = 150.0 * (a_gs*mu_g) / (a_sg*d_p**2*(a_g+a_s)) + 1.75 * (rho_g*u_gs) / (d_p*(a_g+a_s))
                !S_ls_l(mat_inod) = 150.0 * (a_ls*mu_l) / (a_sl*d_p**2*(a_l+a_s)) + 1.75 * (rho_l*u_ls) / (d_p*(a_l+a_s))

                ! for boiling test: two phases-gas and liquid
                S_gs_g(mat_inod) = 0.0 !150.0 * (a_gs*mu_g) / (a_sg*d_p**2*(a_g+a_s)) + 1.75 * (rho_g*u_gs) / (d_p*(a_g+a_s))
                S_ls_l(mat_inod) = 0.0 !150.0 * (a_ls*mu_l) / (a_sl*d_p**2*(a_l+a_s)) + 1.75 * (rho_l*u_ls) / (d_p*(a_l+a_s))


                S_lg_l(mat_inod) = 0.75 * CD * ( (a_gl*rho_l*u_gl) / ( d_b*(a_l+a_g) ) ) * max(a_lg,1.0e-5)**(-2.65)
                S_lg_g(mat_inod) = 0.75 * CD * ( (a_lg*rho_l*u_gl) / ( d_b*(a_l+a_g) ) ) * max(a_lg,1.0e-5)**(-2.65)

            end do

        end do


        !dummy => extract_scalar_field( states(1), "S_gs", stat )
        !if (stat==0) dummy%val = S_gs
        !dummy => extract_scalar_field( states(1), "S_ls", stat )
        !if (stat==0) dummy%val = S_ls
        !dummy => extract_scalar_field( states(1), "S_lg", stat )
        !if (stat==0) dummy%val = S_lg


        ewrite(3,*) 'S_gs_g min/max:', minval(S_gs_g), maxval(S_gs_g)
        ewrite(3,*) 'S_ls_l min/max:', minval(S_ls_l), maxval(S_ls_l)
        ewrite(3,*) 'S_lg_g min/max:', minval(S_lg_g), maxval(S_lg_g)
        ewrite(3,*) 'S_lg_l min/max:', minval(S_lg_l), maxval(S_lg_l)


        deallocate( ul, ug, us )

        return
    end subroutine calculate_boiling_drag



    !sprint_to_do!delete before the sprint is over
    subroutine calculate_boiling_variables( states, packed_state, ndim, nphase, &
        T_sat, Svap_l, Svap_g, Gamma_l, Gamma_g, h_l, h_g, St_gl, St_sl, St_sg )
        implicit none

        type( state_type ), dimension( : ), intent( inout ) :: states
        type( state_type ), intent( in ) :: packed_state
        integer, intent( in ) :: ndim, nphase
        real, dimension( : ), intent( inout ) :: T_sat, Svap_l, Svap_g, Gamma_l, Gamma_g, h_l, h_g, St_gl, St_sl, St_sg

        type( scalar_field ), pointer :: dummy
        type( tensor_field ), pointer :: pressure
        type( tensor_field ), pointer :: density, velocity, temperature, volume_fraction

        integer, dimension( : ), pointer :: cv_ndgln, u_ndgln, xu_ndgln


        real :: p, rho_l, rho_g, u_l, u_g, u_s, t_l, t_g, Tsat, a_l, a_g, a_lg, a_gl, &
            &  a_b, Lh, d_b, dtr, rho_v, Svap_l_max, Svap_g_max, Svap_l1, Svap_l2, &
            &  C, phi, F5, Re_sl, Re_sg, Re_gl, Re_lg, Pr_l, Pr_g

        real, dimension( : ), allocatable :: ul, ug, us, cnt
        real, dimension( :, :, : ), allocatable :: uc

        integer :: ele, totele, cv_iloc, cv_nloc, u_nloc, cv_inod, &
            u_inod, u_inod1, u_inod2, u_iloc, xu_inod, &
            iphase, idim, xu_nonods, stat

        real, parameter :: &
            k_l = 0.58, k_g = 0.016, k_s = 16.2, &
            Cp_l = 4200.0, Cp_g = 1996.0, Cp_s = 500.0, &
            mu_l = 3.0e-4, mu_g = 1.0e-5, &
            d_p = 0.005, &
            Le0 = 2375.7e3, Csf = 0.006, g = 9.81

        pressure => extract_tensor_field( packed_state, "PackedCVPressure" )
        density => extract_tensor_field( packed_state, "PackedDensity" )
        velocity => extract_tensor_field( packed_state, "PackedNonlinearVelocity" )
        temperature => extract_tensor_field( packed_state, "PackedTemperature" )
        volume_fraction => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )

        cv_ndgln => get_ndglno( extract_mesh( packed_state, "PressureMesh" ) )
        u_ndgln => get_ndglno( extract_mesh( packed_state, "VelocityMesh" ) )
        xu_ndgln => get_ndglno( extract_mesh( packed_state, "VelocityMesh_Continuous" ) )


        totele = ele_count( pressure )
        cv_nloc = ele_loc( pressure, 1 )
        u_nloc = ele_loc( velocity, 1 )
        xu_nonods = node_count( extract_mesh( packed_state, "VelocityMesh_Continuous" )  )

        allocate( uc(ndim,nphase,xu_nonods), cnt(xu_nonods) ) ; uc=0.0 ; cnt=0
        do ele = 1, totele
            do u_iloc = 1, u_nloc
                xu_inod = xu_ndgln( ( ele - 1 ) * u_nloc + u_iloc )
                u_inod = u_ndgln( ( ele - 1 ) * u_nloc + u_iloc )
                do iphase = 1, nphase
                    uc( :, iphase, xu_inod ) = uc( :, iphase, xu_inod ) + velocity%val(:,iphase,u_inod)
                end do
                cnt( xu_inod ) = cnt( xu_inod ) + 1
            end do
        end do
        do iphase = 1, nphase
            do idim = 1, ndim
                uc(idim,iphase,:) = uc(idim,iphase,:) / cnt
            end do
        end do

        allocate( ul(ndim), ug(ndim), us(ndim) ) ; ul=0.0 ; ug=0.0 ; us=0.0

        do ele = 1, totele

            do cv_iloc = 1, cv_nloc

                cv_inod = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )

                Svap_l(cv_inod)=0.0 ; Svap_g(cv_inod)=0.0
                Gamma_l(cv_inod)=0.0 ; Gamma_g(cv_inod)=0.0
                h_l(cv_inod)=0.0 ; h_g(cv_inod)=0.0
                St_gl(cv_inod)=0.0 ; St_sl(cv_inod)=0.0 ; St_sg(cv_inod)=0.0

                p = pressure%val(1,1,cv_inod)

                rho_l = density%val(1,1,cv_inod) ; rho_g = density%val(1,2,cv_inod)

                if ( cv_iloc==1  ) then
                    u_inod = xu_ndgln( ( ele - 1 ) * u_nloc + 1 )
                    ul = uc(:,1,u_inod) ; ug = uc(:,2,u_inod) ; us = uc(:,3,u_inod)
                end if
                if ( cv_iloc==3 ) then
                    u_inod = xu_ndgln( ( ele - 1 ) * u_nloc + 2 )
                    ul = uc(:,1,u_inod) ; ug = uc(:,2,u_inod) ; us = uc(:,3,u_inod)
                end if
                if ( cv_iloc==6 ) then
                    u_inod = xu_ndgln( ( ele - 1 ) * u_nloc + 3 )
                    ul = uc(:,1,u_inod) ; ug = uc(:,2,u_inod) ; us = uc(:,3,u_inod)
                end if
                if ( cv_iloc==10 ) then
                    u_inod = xu_ndgln( ( ele - 1 ) * u_nloc + 4 )
                    ul = uc(:,1,u_inod) ; ug = uc(:,2,u_inod) ; us = uc(:,3,u_inod)
                end if
                u_inod1=-666 ; u_inod2=-666
                if ( cv_iloc==2 ) then
                    u_inod1 = xu_ndgln( ( ele - 1 ) * u_nloc + 1 )
                    u_inod2 = xu_ndgln( ( ele - 1 ) * u_nloc + 2 )
                end if
                if ( cv_iloc==4 ) then
                    u_inod1 = xu_ndgln( ( ele - 1 ) * u_nloc + 1 )
                    u_inod2 = xu_ndgln( ( ele - 1 ) * u_nloc + 3 )
                end if
                if ( cv_iloc==5 ) then
                    u_inod1 = xu_ndgln( ( ele - 1 ) * u_nloc + 2 )
                    u_inod2 = xu_ndgln( ( ele - 1 ) * u_nloc + 3 )
                end if
                if ( cv_iloc==7 ) then
                    u_inod1 = xu_ndgln( ( ele - 1 ) * u_nloc + 1 )
                    u_inod2 = xu_ndgln( ( ele - 1 ) * u_nloc + 4 )
                end if
                if ( cv_iloc==8 ) then
                    u_inod1 = xu_ndgln( ( ele - 1 ) * u_nloc + 2 )
                    u_inod2 = xu_ndgln( ( ele - 1 ) * u_nloc + 4 )
                end if
                if ( cv_iloc==9 ) then
                    u_inod1 = xu_ndgln( ( ele - 1 ) * u_nloc + 3 )
                    u_inod2 = xu_ndgln( ( ele - 1 ) * u_nloc + 4 )
                end if
                if ( u_inod1>0 ) then
                    ul = (uc(:,1,u_inod1)+uc(:,1,u_inod2))/2.0
                    ug = (uc(:,2,u_inod1)+uc(:,2,u_inod2))/2.0
                    us = (uc(:,3,u_inod1)+uc(:,3,u_inod2))/2.0
                end if

                u_l = sqrt( sum( ul**2 ) )
                u_g = sqrt( sum( ug**2 ) )
                u_s = sqrt( sum( us**2 ) )

                t_l = temperature%val(1,1,cv_inod) ; t_g = temperature%val(1,2,cv_inod)
                Tsat = saturation_temperature( p )
                T_sat(cv_inod) = Tsat

                a_l = volume_fraction%val(1,1,cv_inod) ; a_g = volume_fraction%val(1,2,cv_inod)
                a_lg=a_l/(a_l+a_g) ; a_gl=1.0-a_lg
                a_b = max(a_g,1.0e-5)

                if (Tsat<t_l) then
                    h_l(cv_inod) = -Le0 + Cp_l*t_l + p/rho_l
                else
                    h_l(cv_inod) = -Le0 + Cp_l*Tsat + p/rho_l
                end if

                if (Tsat<t_g) then
                    h_g(cv_inod) = Cp_g*Tsat + p/rho_g
                else
                    h_g(cv_inod) = Cp_g*t_g  + p/rho_g
                end if

                Lh = h_g(cv_inod) - h_l(cv_inod)

                !h_l(cv_inod) = -Le0 + Cp_l*Tsat + p/rho_l
                !h_g(cv_inod) = h_l(cv_inod) + Lh

                d_b = 5.0*0.06/max((rho_l*abs(u_g-u_l)**2),1.0e-5) ; d_b=min(0.5*d_p,max(1.0e-7,d_b))
                !d_b = 1.0*d_p

                Re_sl = rho_l*abs(u_s-u_l)*d_p/mu_l
                Re_sg = rho_g*abs(u_s-u_g)*d_p/mu_g
                Re_gl = rho_l*abs(u_g-u_l)*d_b/mu_l ; Re_lg=Re_gl

                Pr_l = Cp_l*mu_l/k_l ; Pr_g = Cp_g*mu_g/k_g

                St_gl(cv_inod) = (k_l/d_b)*(2.0+0.6*Re_gl**0.5*Pr_l**0.3333) + 1.0e+7 * a_lg**10
                St_sl(cv_inod) = (k_l/d_p)*(2.0+0.6*Re_sl**0.5*Pr_l**0.3333)
                St_sg(cv_inod) = (k_g/d_p)*(2.0+0.6*Re_sg**0.5*Pr_g**0.3333)

                dtr = 1.0e-2 ; rho_v = 1.0 * rho_g
                if (Tsat<t_l) then
                    Svap_l_max = ((min(rho_v,a_l*rho_l)/dtr)*Lh)/max(1.0e-10,abs(Tsat-t_l))
                else
                    Svap_l_max = (((a_g*rho_v)/dtr)*Lh)/max(1.0e-10,abs(Tsat-t_l))
                end if
                if (Tsat<t_g) then
                    Svap_g_max = ((min(rho_v,a_l*rho_l)/dtr)*Lh)/max(1.0e-10,abs(Tsat-t_g))
                else
                    Svap_g_max = (((a_g*rho_v)/dtr)*Lh)/max(1.0e-10,abs(Tsat-t_g))
                end if

                if(Tsat<t_l) then
                    Svap_l1 = (k_l/d_b)*(12.0/pi)*abs(Tsat-t_l)*(rho_l*Cp_l)/(rho_g*Lh)
                    Svap_l2 = (k_l/d_b)*(2.0+0.74*(a_l*Re_lg)**0.5)
                    Svap_l(cv_inod) = max(Svap_l1,Svap_l2)*3.6*a_b/d_b
                else
                    if(p<=1.1272*1.0e6) then
                        C=65.0-5.69e-5*(p-1.0e5)
                    else
                        C=2.5e9*p**(-1.418)
                    end if

                    if(abs(u_g-u_l)<=0.61) then
                        phi=1.0
                    else
                        phi=(1.639344*abs(u_g-u_l))**0.47
                    end if

                    if(a_g<0.25) then
                        F5=0.075+1.8*phi*C*exp(-45.0*a_b)
                    else
                        F5=0.075
                    end if

                    Svap_l(cv_inod) = F5*Lh*rho_g*rho_l*a_g/(rho_l-rho_g)
                    Svap_l(cv_inod) = min( Svap_l(cv_inod), 17539.0*max(4.724,472.4*a_g*a_l)*max(0.0,min(1.0,a_g/0.1)) )

                end if

                Svap_g(cv_inod) = 1.0e4*3.6*a_b/d_b

                Svap_l(cv_inod) = min( Svap_l(cv_inod), Svap_l_max)
                Svap_g(cv_inod) = min( Svap_g(cv_inod), Svap_g_max)

                Gamma_g(cv_inod) = (Svap_l(cv_inod)*(t_l-Tsat)+Svap_g(cv_inod)*(t_g-Tsat))/Lh
                Gamma_l(cv_inod) = -Gamma_g(cv_inod)

            end do
        end do


        dummy => extract_scalar_field( states(1), "Gamma_l", stat )
        if (stat==0) dummy%val = Gamma_l
        dummy => extract_scalar_field( states(1), "Gamma_g", stat )
        if (stat==0) dummy%val = Gamma_g

        dummy => extract_scalar_field( states(1), "Svap_l", stat )
        if (stat==0) dummy%val = Svap_l
        dummy => extract_scalar_field( states(1), "Svap_g", stat )
        if (stat==0) dummy%val = Svap_g

        dummy => extract_scalar_field( states(1), "h_l", stat )
        if (stat==0) dummy%val = h_l
        dummy => extract_scalar_field( states(1), "h_g", stat )
        if (stat==0) dummy%val = h_g

        dummy => extract_scalar_field( states(1), "St_gl", stat )
        if (stat==0) dummy%val = St_gl
        dummy => extract_scalar_field( states(1), "St_sl", stat )
        if (stat==0) dummy%val = St_sl
        dummy => extract_scalar_field( states(1), "St_sg", stat )
        if (stat==0) dummy%val = St_sg


        ewrite(3,*) 'Gamma_l min/max:', minval(Gamma_l), maxval(Gamma_l)
        ewrite(3,*) 'Gamma_g min/max:', minval(Gamma_g), maxval(Gamma_g)

        ewrite(3,*) 'Svap_l min/max:', minval(Svap_l), maxval(Svap_l)
        ewrite(3,*) 'Svap_g min/max:', minval(Svap_g), maxval(Svap_g)


        ewrite(3,*) 'h_l min/max:', minval(h_l), maxval(h_l)
        ewrite(3,*) 'h_g min/max:', minval(h_g), maxval(h_g)

        ewrite(3,*) 'T_sat min/max:', minval(T_sat), maxval(T_sat)

        ewrite(3,*) 'St_gl min/max:', minval(St_gl), maxval(St_gl)
        ewrite(3,*) 'St_sl min/max:', minval(St_sl), maxval(St_sl)
        ewrite(3,*) 'St_sg min/max:', minval(St_sg), maxval(St_sg)

        deallocate( ul, ug, us )

        return
    end subroutine calculate_boiling_variables



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

    subroutine get_RockFluidProp(state, packed_state)
        !Gets the relperm max, the relperm exponent and the immobile fractions and stores
        !them into packed state
        !The first index is the immobile fraction, the second is the relperm max
        ! and the third is the relperm exponent
        implicit none
        type(state_type), dimension(:), intent(inout) :: state
        type( state_type ), intent( inout ) :: packed_state
        !Local variables
        type (tensor_field), pointer :: t_field
        type (scalar_field), target :: targ_Store
        type (scalar_field), pointer :: s_field
        type (vector_field), pointer :: position
        type(mesh_type), pointer :: fl_mesh
        type(mesh_type) :: Auxmesh
        integer :: iphase, nphase
        character(len=500) :: path

        t_field=>extract_tensor_field(packed_state,"PackedRockFluidProp")
        nphase = size(t_field%val,2)
        !By default the pressure mesh (position 1)
        s_field => extract_scalar_field(state(1),1)
        position => get_external_coordinate_field(packed_state, s_field%mesh)

        fl_mesh => extract_mesh( state(1), "P0DG" )
        Auxmesh = fl_mesh
        call allocate (targ_Store, Auxmesh, "Temporary_get_RockFluidProp")

        !Retrieve Immobile fractions
        do iphase = 1, nphase
            path = "/material_phase["//int2str(iphase-1)//&
                "]/multiphase_properties/immobile_fraction/scalar_field::value/prescribed/value"
            if (have_option(trim(path))) then
                call initialise_field_over_regions(targ_Store, trim(path) , position)
                t_field%val(1,iphase,:) = targ_Store%val(:)
            else!default value
                t_field%val(1,iphase,:) = 0.0
            end if
        end do

        !Retrieve relperm max
        do iphase = 1, nphase
            path = "/material_phase["//int2str(iphase-1)//&
                "]/multiphase_properties/Relperm_Corey/relperm_max/scalar_field::relperm_max/prescribed/value"
            if (have_option(trim(path))) then
                call initialise_field_over_regions(targ_Store, trim(path) , position)
                t_field%val(2,iphase,:) = max(min(targ_Store%val(:), 1.0), 0.0)
            else!default value
                t_field%val(2,iphase,:) = 1.0
            end if
        end do

        !Retrieve relperm exponent
        do iphase = 1, nphase
            path = "/material_phase["//int2str(iphase-1)//&
                "]/multiphase_properties/Relperm_Corey/relperm_exponent/scalar_field::relperm_exponent/prescribed/value"
            if (have_option(trim(path))) then
                call initialise_field_over_regions(targ_Store, trim(path) , position)
                t_field%val(3,iphase,:) = targ_Store%val(:)
            else!default value
                t_field%val(3,iphase,:) = 2.0
            end if
        end do

        !Initialize capillary pressure
        if (have_option_for_any_phase( '/multiphase_properties/capillary_pressure', nphase ) ) then
            !Get cap entry pressure
            do iphase = 1, nphase
                path = "/material_phase["//int2str(iphase-1)//&
                    "]/multiphase_properties/capillary_pressure/type_Brooks_Corey/scalar_field::C/prescribed/value"
                if (have_option(trim(path))) then
                    call initialise_field_over_regions(targ_Store, trim(path) , position)
                    t_field%val(4,iphase,:) = targ_Store%val(:)
                else!default value
                    t_field%val(4,iphase,:) = 0.0
                end if
            end do

            !Get cap exponent
            do iphase = 1, nphase
                path = "/material_phase["//int2str(iphase-1)//&
                    "]/multiphase_properties/capillary_pressure/type_Brooks_Corey/scalar_field::a/prescribed/value"
                if (have_option(trim(path))) then
                    call initialise_field_over_regions(targ_Store, trim(path) , position)
                    t_field%val(5,iphase,:) = targ_Store%val(:)
                else!default value
                    t_field%val(5,iphase,:) = 1.0
                end if
            end do
        end if

        call deallocate(targ_Store)
    end subroutine get_RockFluidProp






    !!JWL eqaution functions

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

!!-JWL eqaution functions







end module multiphase_EOS
