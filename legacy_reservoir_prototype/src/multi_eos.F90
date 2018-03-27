
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
    use global_parameters
    use spud
    use futils, only: int2str
    use vector_tools
    use python_state
    use Copy_Outof_State
    use multi_tools !!!!! WHY is the bad_elements type not being picked up by Copy_Outof_State??????

    use global_parameters, only: domain_bbox
    use shape_functions_Linear_Quadratic
    use sparse_tools
    use Multiphase_module
    use sparsity_patterns_meshes, only: get_csr_sparsity_firstorder
    use arbitrary_function
    use boundary_conditions, only: get_entire_boundary_condition
    use Field_Options, only: get_external_coordinate_field
    use initialise_fields_module, only: initialise_field_over_regions, initialise_field
    use multi_tools, only: CALC_FACE_ELE, assign_val, table_interpolation, read_csv_table
    implicit none

    real, parameter :: flooding_hmin = 1e-5

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
        type( scalar_field ), pointer :: Cp_s, Density
        integer :: icomp, iphase, ncomp, sc, ec, sp, ep, ip, stat, cv_iloc, cv_nod, ele
        logical :: boussinesq
        logical, parameter :: harmonic_average=.false.


        !For simple Black-Oil modelling the density have to account for the disolved gas in it
        !and the three phases are modified simultaneously. Hence other models are overwritten
        if (have_option( "/physical_parameters/black-oil_PVT_table" ))then
            if (Mdims%nphase == 3.and. Mdims%ncomp<1) then
                call extended_Black_Oil(state, packed_state, Mdims, flash_flag = 1)
                return
            else
                ewrite(1,*) "WARNING: Black-Oil model activated but three phases are not present and/or there are components"
            end if
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
        allocate( Cp( cv_nonods ) ) ; Cp = 1.0
        allocate( Density_Component( ncomp * nphase * cv_nonods ) )
        allocate( Density_Bulk( nphase * cv_nonods ), DensityCp_Bulk( nphase * cv_nonods ) )
        Density_Bulk = 0.0 ; DensityCp_Bulk = 0.0

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

                 if ( .not.harmonic_average ) then

                    ! rho = rho +  a_i * rho_i
                    Density_Bulk( sp : ep ) = Density_Bulk( sp : ep ) + Rho * Component_l
                    PackedDRhoDPressure%val( 1, iphase, : ) = PackedDRhoDPressure%val( 1, iphase, : ) + dRhodP * Component_l / Rho

                    Density_Component( sc : ec ) = Rho

                    Cp_s => extract_scalar_field( state( nphase + icomp ), &
                         'ComponentMassFractionPhase' // int2str( iphase ) // 'HeatCapacity', stat )
                    if ( stat == 0 ) call assign_val(Cp,Cp_s % val)!Cp = Cp_s % val
                    DensityCp_Bulk( sp : ep ) = DensityCp_Bulk( sp : ep ) + Rho * Cp * Component_l

                 else

                    Density_Bulk( sp : ep ) = Density_Bulk( sp : ep ) + Rho * Component_l
                    PackedDRhoDPressure%val( 1, iphase, : ) = PackedDRhoDPressure%val( 1, iphase, : ) + dRhodP * Component_l / Rho

                    Density_Component( sc : ec ) = Rho

                    ! harmonic average
                    ! rho = rho + 1.0 / ( a_i / rho_i )
                    Cp_s => extract_scalar_field( state( nphase + icomp ), &
                         'ComponentMassFractionPhase' // int2str( iphase ) // 'HeatCapacity', stat )
                    if ( stat == 0 ) call assign_val(Cp,Cp_s % val)!Cp = Cp_s % val

                    do cv_nod = 1, cv_nonods
                       ip = ( iphase - 1 ) * cv_nonods + cv_nod
                       DensityCp_Bulk( ip ) = DensityCp_Bulk( ip ) + Component_l(cv_nod) / ( Rho(cv_nod) * Cp(cv_nod) )
                    end do

                 end if

              else

                 Density_Bulk( sp : ep ) = Rho
                 PackedDRhoDPressure%val( 1, iphase, : ) = dRhodP

                 Cp_s => extract_scalar_field( state( iphase ), 'TemperatureHeatCapacity', stat )
                 if ( stat == 0 ) call assign_val(Cp,Cp_s % val)
                 DensityCp_Bulk( sp : ep ) = Rho * Cp

              end if

           end do ! iphase
        end do ! icomp

        if ( ncomp > 1 ) then
           if ( harmonic_average ) DensityCp_Bulk = 1.0 / DensityCp_Bulk
           call Cap_Bulk_Rho( state, ncomp, nphase, &
                cv_nonods, Density_Component, Density_Bulk, DensityCp_Bulk )
        end if

        field1 => extract_tensor_field( packed_state, "PackedDensity" )
        field2 => extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )
        if( ncomp > 1 ) field3 => extract_tensor_field( packed_state, "PackedComponentDensity" )

        do iphase = 1, nphase
           sp = ( iphase - 1 ) * cv_nonods + 1
           ep = iphase * cv_nonods

           field1 % val ( 1, iphase, : ) = Density_Bulk( sp : ep )
           field2 % val ( 1, iphase, : ) = DensityCp_Bulk( sp : ep )

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

        do iphase = 1, nphase
           Density => extract_scalar_field( state( iphase ), "Density" )
           Density % val = field1 % val( 1, iphase, : )
        end do

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

        allocate( Density_Component_Min( nphase, cv_nonods ) ) ; Density_Component_Min = 1.0e+15
        allocate( Density_Component_Max( nphase, cv_nonods ) ) ; Density_Component_Max = 0.0
        allocate( Density_Cp_Component_Min( nphase, cv_nonods ) ) ; Density_Cp_Component_Min = 1.0e+15
        allocate( Density_Cp_Component_Max( nphase, cv_nonods ) ) ; Density_Cp_Component_Max = 0.0
        allocate( Cp( cv_nonods ) ) ; Cp = 1.0

        do iphase = 1, nphase
            do icomp = 1, ncomp
                sc = ( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1
                ec = ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods

                Density_Component_Min( iphase, : ) = min( Density_Component_Min( iphase, : ), Density_Component( sc : ec ) )
                Density_Component_Max( iphase, : ) = max( Density_Component_Max( iphase, : ), Density_Component( sc : ec ) )

                Cp = 1.0
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
                 '/prognostic/equation_of_state' )
            call Assign_Equation_of_State( eos_option_path )
            Rho=0. ; dRhodP=0.
            call Calculate_Rho_dRhodP( state, packed_state, iphase, icomp, &
                 Mdims%nphase, Mdims%ncomp, eos_option_path, Rho, dRhodP )

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
        !Variables for python function for the coefficient_B for linear density (this is for bathymetry)
        type (scalar_field) :: sfield
        type (scalar_field), pointer :: pnt_sfield
        type (vector_field), pointer :: position

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



        elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/linear_in_pressure/include_internal_energy' ) then
            !!$ Den = C0 * P/T +C1
            if( .not. have_temperature_field ) FLAbort( 'Temperature Field not defined' )
            allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
            call get_option( trim( option_path_comp ) // '/linear_in_pressure/coefficient_A', eos_coefs( 1 ) )
            call get_option( trim( option_path_comp ) // '/linear_in_pressure/coefficient_B/constant', eos_coefs( 2 ) )
            Rho = eos_coefs( 1 ) * pressure % val(1,1,:) / temperature % val + eos_coefs( 2 )
            perturbation_pressure = 1.
            !RhoPlus = eos_coefs( 1 ) * ( pressure % val + perturbation_pressure ) / &
            !     ( max( toler, temperature % val ) ) + eos_coefs( 2 )
            !RhoMinus = eos_coefs( 1 ) * ( pressure % val - perturbation_pressure ) / &
            !     ( max( toler, temperature % val ) ) + eos_coefs( 2 )
            dRhodP =  eos_coefs( 1 ) / temperature % val !0.5 * ( DensityPlus - DensityMinus ) / perturbation_pressure
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
            !If it is flooding, then the second coefficient is the bathymetry and it has to be negative.
            !The formula becomes: rho = P * eos_coefs( 1 ) - bathymetry
            if (is_flooding .and. minval(sfield%val, sfield%val> 0.) > 0.) sfield%val = - sfield%val
            Rho = eos_coefs( 1 ) * pressure % val(1,1,:) + sfield%val
            perturbation_pressure = 1.
            !RhoPlus = eos_coefs( 1 ) * ( pressure % val + perturbation_pressure ) + eos_coefs( 2 )
            !RhoMinus = eos_coefs( 1 ) * ( pressure % val - perturbation_pressure ) + eos_coefs( 2 )
            dRhodP = eos_coefs( 1 ) !0.5 * ( DensityPlus - DensityMinus ) / perturbation_pressure
            deallocate( eos_coefs )
            call deallocate(sfield)
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
        elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/Temperature_Pressure_correlation' ) then
            !!$ den = den0/(1+Beta(T1-T0))/(1-(P1-P0)/E)
            allocate( temperature_local( node_count( pressure ) ) ) ; temperature_local = 0.
            if ( have_temperature_field ) temperature_local = max(temperature % val,1e-8)!avoid possible oscillations introduced by unphysical values of temperature
                                                                                        !appearing while achieving convergence

            allocate( eos_coefs( 5 ) ) ; eos_coefs = 0.
            call get_option( trim( option_path_comp ) // '/Temperature_Pressure_correlation/rho0', eos_coefs( 1 ) )
            call get_option( trim( option_path_comp ) // '/Temperature_Pressure_correlation/T0/', eos_coefs( 2 ), default = 0. )
            call get_option( trim( option_path_comp ) // '/Temperature_Pressure_correlation/P0/', eos_coefs( 3 ) )
            call get_option( trim( option_path_comp ) // '/Temperature_Pressure_correlation/coefficient_Beta/', eos_coefs( 4 ), default = 0. )
            call get_option( trim( option_path_comp ) // '/Temperature_Pressure_correlation/coefficient_E/', eos_coefs( 5 ), default = 0. )

            !!$ den = den0/(1+Beta(T1-T0))
            !We use RHo as auxiliar variable here as the we do not perturbate the temperature

            Rho = eos_coefs(1)/(1 + eos_coefs(4)*(temperature_local-eos_coefs(2) )  )

            perturbation_pressure = max( toler, 1.e-3 * abs( pressure % val(1,1,:) ) )
            !we add the pressure part =>1-(P1-P0)/E
            RhoPlus = Rho /(1-((perturbation_pressure-eos_coefs(3))/eos_coefs(5)))
            RhoMinus = Rho /(1-((perturbation_pressure-eos_coefs(3))/eos_coefs(5)))

            dRhodP =  0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure
            !we add the pressure part =>1-(P1-P0)/E
            Rho = Rho /(1-( (min(max(pressure%val(1,1,:),-101325.),eos_coefs( 5 )*0.5) -eos_coefs(3))/eos_coefs(5)) ) !to avoid possible oscillations the pressure is imposed to be between the range of applicability of the formula.
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

        !For flooding ensure that the height (density of phase 1) is non-zero and positive
        if (is_flooding .and. iphase == 1) Rho = max(Rho, flooding_hmin)

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

            elseif( have_option( trim( eos_option_path_out ) // '/Temperature_Pressure_correlation' ) ) then
                eos_option_path_out = trim( eos_option_path_out ) // '/Temperature_Pressure_correlation'

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

    subroutine Calculate_flooding_absorptionTerm(state, packed_state, Flooding_absorp, Mdims, ndgln)
        implicit none
        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type (multi_field) :: Flooding_absorp
        type( multi_dimensions ), intent( in ) :: Mdims
        type(multi_ndgln), intent(in) :: ndgln
        !Local variables
        integer :: ipres, iphase, idim, loc
        type( scalar_field ), pointer :: Spipe

        if(.not.is_flooding) return !Nothing to do here, return

        !Initialise
        Flooding_absorp%val=0.
        !For the non-pipe phases, add manually the manning coefficient use equation 13 from the manual (the one that defines the b)
        call calculate_manning_coef(Flooding_absorp)

        !The following part is the absorption for the pipes
        if (Mdims%npres > 1) then
            Spipe => extract_scalar_field( state(1), "Sigma" )
            do iphase = Mdims%n_in_pres + 1, Mdims%nphase
                ! set \sigma for the pipes here
                call assign_val(Flooding_absorp%val( 1,1, iphase, : ),Spipe%val)
            end do
        end if

        contains

        subroutine calculate_manning_coef(Flooding_absorp)
            implicit none
            type (multi_field) :: Flooding_absorp
            !Local variables
            real, parameter :: hmin = max(flooding_hmin, 1d-8) * 1.1!The velocity solver is very sensitive to this parameter
            real, parameter :: u_min = 1d-2 !increase it if having problems to converge
            real, parameter :: g = 9.80665!Set default value if not specified by the user
            integer :: iphase, ele, cv_iloc, u_iloc, mat_nod, cv_nod, u_nod,  stat, i
            type( tensor_field ), pointer :: velocity, Nm, density
            type(vector_field), pointer :: gravity_direction
            real, dimension(mdims%cv_nloc) :: bathymetry, Nm_aux
            real, dimension(:), allocatable :: r_nod_count
            logical :: averaging
            real :: shallow_drag

            !Check whether to use the harmonic mean of the bathymetry
            averaging = have_option('/flooding/averaging')
            !Strenght shallow_drag
            call get_option('/flooding/shallow_drag', shallow_drag, default = 1d-1)

            Nm => extract_tensor_field( packed_state, "PackedManningcoef" )!Defined element-wise
            velocity => extract_tensor_field( packed_state, "PackedVelocity" )
            density => extract_tensor_field( packed_state, "PackedDensity" )!For flooding the first phase is the height
            iphase = 1!First phase of velocity only

            allocate(r_nod_count(size(Flooding_absorp%val,4))); r_nod_count = 0.
            do ele = 1, Mdims%totele
                do cv_iloc = 1, Mdims%cv_nloc
                    !Create bathymetry field just in case of using the mean
                    cv_nod = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                    bathymetry(cv_iloc) = max(hmin, density%val(1,1,cv_nod))
                end do                               !Normal mean                          !Harmonic mean
                if (averaging) bathymetry = (sum(bathymetry) / dble(Mdims%cv_nloc))!(sum(bathymetry**-1) / dble(Mdims%cv_nloc))**-1
                do cv_iloc = 1, Mdims%cv_nloc
                    mat_nod = ndgln%mat(( ELE - 1 ) * Mdims%mat_nloc + cv_iloc)
                    cv_nod = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                    r_nod_count(mat_nod) = r_nod_count(mat_nod) + 1
                    Nm_aux(cv_iloc) = Nm%val(1,1,ele) + max(flooding_hmin, shallow_drag*(2*hmin-density%val(1,1,cv_nod))/hmin)
                    do u_iloc = 1, mdims%u_nloc
                        u_nod = ndgln%u(( ELE - 1) * Mdims%u_nloc + u_iloc )
                        !Since Flooding_absorp is of memory_type 1 we can populate it directly
                        do i = 1, Mdims%n_in_pres!Only for the phases not in the pipes
                            Flooding_absorp%val(1,1,i, mat_nod) = Flooding_absorp%val(1,1,i, mat_nod) + Nm_aux(cv_iloc)**2. * g *&
                                max(u_min,sqrt(dot_product(velocity%val(1:2,iphase,u_nod),velocity%val(1:2,iphase,u_nod))))&!We are using only two dimensions of the velocity because
                                /(bathymetry(cv_iloc)**(4./3.)*dble(mdims%u_nloc))!This last term to get an average         !<= this is a 2D model if used in 3D, the third dimension is time!
                        end do
                    end do
                end do
            end do
            !Average considering times node has been visited
            do i = 1, Mdims%n_in_pres
                Flooding_absorp%val(1, 1, i, :) = Flooding_absorp%val(1, 1, i, :) / r_nod_count
            end do
            deallocate(r_nod_count)
        end subroutine calculate_manning_coef
    end subroutine Calculate_flooding_absorptionTerm

    subroutine Calculate_PorousMedia_AbsorptionTerms( state, packed_state, PorousMedia_absorp, Mdims, CV_funs, CV_GIdims, Mspars, ndgln, &
                                                      upwnd, suf_sig_diagten_bc, Quality_list )
       implicit none
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
       type(bad_elements), allocatable, dimension(:), optional :: Quality_list
       !Local variables
       real, save :: kv_kh_ratio = -1
       type( tensor_field ), pointer :: perm
       real, dimension(Mdims%ndim, Mdims%ndim, Mdims%totele), target:: inv_perm
       real, dimension(:,:), allocatable :: viscosities
       integer :: i, j, ele
       real :: Angle, bad_element_perm_mult, Max_aspect_ratio, height ! the height of an isosceles triangle for the top angle to be equal to the trigger angle
       real, dimension(Mdims%ndim,Mdims%ndim) :: trans_matrix, rot_trans_matrix ! for bad_element permeability transformation matrix
       logical :: kv_kh_ratio_log = .False. ! check if we use the kv_kh ratio or Aspect_ratio for the bad elements. Aspect ratio is the default
       real, parameter :: pi = acos(0.0) * 2.0 ! Define pi

       perm => extract_tensor_field( packed_state, "Permeability" )


        if (kv_kh_ratio < 0) then
            if (have_option('/numerical_methods/Bad_element_fix/') ) then
                call get_option('/numerical_methods/Bad_element_fix/KvKh_ratio', kv_kh_ratio, default = 0.01)
            else
                kv_kh_ratio = 0.
            end if
            kv_kh_ratio = abs(kv_kh_ratio)
        end if

        if ( present(Quality_list) .and. kv_kh_ratio > 1d-8 ) then
            if (allocated(Quality_list)) then
            ! create transformation matrix with Kv/kh ratio
            ! |1    0    0    |  |1    0     |
            ! |0    1    0    |  |0    kv/kh |
            ! |0    0   kv/kh |
            trans_matrix(:,:) = 0.
            do i=1,Mdims%ndim-1
                trans_matrix(i,i) = 1.
            end do
            trans_matrix(Mdims%ndim,Mdims%ndim) = kv_kh_ratio
            i=1
            do while ( Quality_list(i)%bad_ele > 0 )
                ele = Quality_list(i)%bad_ele
                rot_trans_matrix = matmul(transpose(Quality_list(ele)%rotmatrix(1:Mdims%ndim,1:Mdims%ndim)) , matmul(trans_matrix,Quality_list(ele)%rotmatrix(1:Mdims%ndim,1:Mdims%ndim))  )
                perm%val(:, :, ele) =  matmul(rot_trans_matrix, perm%val(:, :, ele))
                i = i+1
            end do

            i=1
            do while (allocated(Quality_list(i)%rotmatrix) )
                deallocate(Quality_list(i)%rotmatrix)
                i = i+1
            end do
            deallocate(Quality_list)
            end if
        end if

       if (PorousMedia_absorp%memory_type<2) then!The permeability is isotropic
           inv_perm = 0.
           do i = 1, size(perm%val,3)
               do j = 1, size(perm%val,1)
                   inv_perm( j, j, i)=1.0/perm%val( j, j, i)
               end do
           end do
       else
           do i = 1, size(perm%val,3)
               inv_perm( :, :, i)=inverse(perm%val( :, :, i))
           end do
       end if

       !For simple Black-Oil modelling the viscosity is calculated using the PVT tables
       if (have_option( "/physical_parameters/black-oil_PVT_table" ) .and. Mdims%ncomp<1)then
           allocate(viscosities(Mdims%nphase, Mdims%cv_nonods))
           call extended_Black_Oil(state, packed_state, Mdims, flash_flag = 3, viscosities = viscosities)
       else
            allocate(viscosities(Mdims%nphase, 1))
            call set_viscosity(state, Mdims, viscosities(:,1))
       end if
       call Calculate_PorousMedia_adv_terms( state, packed_state, PorousMedia_absorp, Mdims, ndgln, &
              upwnd, inv_perm, viscosities)

       ! calculate SUF_SIG_DIAGTEN_BC this is \sigma_in^{-1} \sigma_out
       ! \sigma_in and \sigma_out have the same anisotropy so SUF_SIG_DIAGTEN_BC
       ! is diagonal
       call calculate_SUF_SIG_DIAGTEN_BC( packed_state, suf_sig_diagten_bc, Mdims, CV_funs, CV_GIdims, &
                              Mspars, ndgln, PorousMedia_absorp, state, inv_perm, viscosities)

       deallocate(viscosities)
       contains

           subroutine Calculate_PorousMedia_adv_terms( state, packed_state, PorousMedia_absorp, Mdims, ndgln, &
               upwnd, inv_perm, viscosities )

               implicit none
               type( state_type ), dimension( : ), intent( in ) :: state
               type( state_type ), intent( inout ) :: packed_state
               type (multi_field), intent( inout ) :: PorousMedia_absorp
               type( multi_dimensions ), intent( in ) :: Mdims
               type( multi_ndgln ), intent( in ) :: ndgln
               type (porous_adv_coefs), intent(inout) :: upwnd
               real, dimension(:, :, :), target, intent(in):: inv_perm
               real, dimension(:,:), intent(in) :: viscosities
               !!$ Local variables:
               integer :: ele, imat, icv, iphase, cv_iloc, idim, jdim, ipres, loc
               real :: Mobility, pert
               real, dimension(:), allocatable :: Max_sat
               real, dimension( :, : ), allocatable :: satura2
               type (multi_field) :: PorousMedia_absorp2
               !Working pointers
               real, dimension(:,:), pointer :: Satura, OldSatura, Immobile_fraction
               type( tensor_field ), pointer :: perm
               type( scalar_field ), pointer :: Spipe

               !Initialize variables
               upwnd%adv_coef_grad=0.0;upwnd%inv_adv_coef=0.0
               !Get from packed_state
               call get_var_from_packed_state(packed_state,PhaseVolumeFraction = Satura,&
                   OldPhaseVolumeFraction = OldSatura, Immobile_fraction = Immobile_fraction)
               perm=>extract_tensor_field(packed_state,"Permeability")


               allocate( satura2( Mdims%n_in_pres, size(SATURA,2) ) );satura2 = 0.

               CALL calculate_absorption2( packed_state, PorousMedia_absorp, Mdims, ndgln, SATURA(1:Mdims%n_in_pres,:), &
                   PERM%val, viscosities, inv_perm1=inv_perm)

               !Introduce perturbation, positive for the increasing and negative for decreasing phase
               !Make sure that the perturbation is between bounds
               PERT = 0.0001; allocate(Max_sat(Mdims%nphase))
               do ele = 1, Mdims%totele
                   do cv_iloc = 1, Mdims%cv_nloc
                       icv = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + cv_iloc )
                       Max_sat(:) = 1. - sum(Immobile_fraction(:, ele)) + Immobile_fraction(:, ele)
                       do iphase = 1, Mdims%n_in_pres !Mdims%nphase
                           SATURA2(iphase, icv) = SATURA(iphase, icv) + sign(PERT, satura(iphase, icv)-OldSatura(iphase, icv))
                           !If out of bounds then we perturbate in the opposite direction
                           if (satura2(iphase, icv) > Max_sat(iphase) .or. &
                               satura2(iphase, icv) < Immobile_fraction(iphase, ele)) then
                               SATURA2(iphase, icv) = SATURA2(iphase, icv) - 2. * sign(PERT, satura(iphase, icv)-OldSatura(iphase, icv))
                           end if
                       end do
                   end do
               end do


               call allocate_multi_field( Mdims, PorousMedia_absorp2, size(PorousMedia_absorp%val,4), field_name="PorousMedia_AbsorptionTerm")
               CALL calculate_absorption2( packed_state, PorousMedia_absorp2, Mdims, ndgln, SATURA2, &
                   PERM%val, viscosities, inv_perm1=inv_perm)

               do ipres = 2, Mdims%npres
                   Spipe => extract_scalar_field( state(1), "Sigma" )
                   do iphase = Mdims%n_in_pres+1, Mdims%nphase
                       do idim = 1, Mdims%ndim
                           ! set \sigma for the pipes here
                           call assign_val(PorousMedia_absorp%val(idim, idim, iphase, :),Spipe%val)
                       end do
                   end do
               end do

               !Temporary pointer, maybe we should unify memories
               nullify(upwnd%adv_coef)
               upwnd%adv_coef => PorousMedia_absorp%val

               DO ELE = 1, Mdims%totele
                   DO CV_ILOC = 1, Mdims%cv_nloc
                       IMAT = ndgln%mat( ( ELE - 1 ) * Mdims%mat_nloc + CV_ILOC )
                       ICV = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
                       DO IPHASE = 1, Mdims%nphase
                           DO JDIM = 1, Mdims%ndim
                               DO IDIM = 1, Mdims%ndim
                                   if ( iphase <= Mdims%n_in_pres ) then
                                       ! This is the gradient
                                       ! Assume d\sigma / dS = 0.0 for the pipes for now
                                       upwnd%adv_coef_grad(IDIM, JDIM, IPHASE, IMAT) = (PorousMedia_absorp2%val( idim,jdim, iphase ,IMAT) -&
                                           PorousMedia_absorp%val( idim,jdim, iphase ,IMAT)) / ( SATURA2(IPHASE, ICV ) - SATURA(IPHASE, ICV))
                                   end if
                               END DO
                           !Obtaining the inverse the "old way" since if you obtain it directly, some problems appear
                           upwnd%inv_adv_coef(:, :, IPHASE, IMAT) = inverse(upwnd%adv_coef(:, :, IPHASE, IMAT))
                           END DO
                       END DO
                   END DO
               END DO

               deallocate( satura2, Max_sat)
               call deallocate_multi_field(PorousMedia_absorp2, .true.)

           end subroutine Calculate_PorousMedia_adv_terms


           subroutine calculate_SUF_SIG_DIAGTEN_BC( packed_state, suf_sig_diagten_bc, Mdims, CV_funs, CV_GIdims, &
               Mspars, ndgln, PorousMedia_absorp, state, inv_perm, viscosities)
               implicit none
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
               real, dimension( Mdims%ndim, Mdims%ndim ) :: sigma_out, sigma_in, mat, mat_inv
               real, dimension( Mdims%ndim * Mdims%nphase, Mdims%ndim * Mdims%nphase) :: mat_absorp
               integer, dimension( CV_GIdims%nface, Mdims%totele) :: face_ele
               integer, dimension( Mdims%mat_nonods*Mdims%n_in_pres ) :: idone
               integer, dimension( Mdims%cv_snloc ) :: cv_sloc2loc
               integer, dimension( :, :, : ),  allocatable :: wic_u_bc, wic_vol_bc
               integer, parameter :: WIC_BC_DIRICHLET = 1
               type(tensor_field), pointer :: velocity, volfrac, perm
               type(tensor_field) :: velocity_BCs, volfrac_BCs
               integer :: one_or_zero, visc_node
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
               idone=0; face_ele = 0
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
                               do iphase = 1, Mdims%n_in_pres
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
                                           do idim = 1, Mdims%ndim
                                               do jdim = 1, Mdims%ndim
                                                   call get_relperm(Mdims%n_in_pres, iphase, sigma_out( idim, jdim ),&
                                                       ! this is the boundary condition
                                                       volfrac_BCs%val(1,:,cv_snodi), viscosities(:,visc_node), inv_perm( idim, jdim, ele ),&
                                                       Immobile_fraction, Corey_exponent, Endpoint_relperm)
                                               end do
                                           end do
                                           ! Adjust suf_sig_diagten_bc based on the internal absorption
                                           mat = sigma_out  +  matmul(  PorousMedia_absorp%val(:,:,iphase, mat_nod),  &
                                                    matmul( inverse( sigma_out ), PorousMedia_absorp%val(:,:,iphase, mat_nod) ) )
                                           mat_inv = matmul( inverse( PorousMedia_absorp%val(:,:,iphase, mat_nod)+sigma_out ), mat )
                                           suf_sig_diagten_bc( cv_snodi_ipha, 1 : Mdims%ndim ) = (/ (mat_inv(i, i), i = 1, Mdims%ndim) /)
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


        subroutine set_viscosity(state, Mdims, visc_phases)
            implicit none
            type( state_type ), dimension( : ), intent( in ) :: state
            type(multi_dimensions), intent(in) :: Mdims
            real, dimension(:), intent(inout) :: visc_phases
            !Local variables
            integer :: iphase
            real :: mobility
            type(tensor_field), pointer :: viscosity_ph

            if( have_option( '/physical_parameters/mobility' ) )then!This option is misleading, it should be removed
                call get_option( '/physical_parameters/mobility', mobility )
                visc_phases(1) = 1
                visc_phases(2) = mobility
            elseif( have_option( '/material_phase[1]/vector_field::Velocity/prognostic/tensor_field::Viscosity' // &
                '/prescribed/value::WholeMesh/isotropic' ) ) then
                DO IPHASE = 1, Mdims%nphase!Get viscosity for all the phases
                    viscosity_ph => extract_tensor_field( state( iphase ), 'Viscosity' )
                    visc_phases(iphase) = viscosity_ph%val( 1, 1, 1 )!So far we only consider scalar viscosity
                end do
            elseif( Mdims%n_in_pres == 1 ) then
                viscosity_ph => extract_tensor_field( state( 1 ), 'Viscosity' )
                visc_phases(1) = viscosity_ph%val( 1, 1, 1 )
            end if


        end subroutine set_viscosity


    end subroutine Calculate_PorousMedia_AbsorptionTerms





    SUBROUTINE calculate_absorption2( packed_state, PorousMedia_absorp, Mdims, ndgln, SATURA, &
        PERM, viscosities, inv_mat_absorp, inv_perm1)
        ! Calculate absorption for momentum eqns
        implicit none
        type( state_type ), intent( inout ) :: packed_state
        type (multi_field) :: PorousMedia_absorp
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_ndgln), intent(in) :: ndgln
        REAL, DIMENSION( :, : ), intent( in ) :: SATURA
        REAL, DIMENSION( :, :, : ), intent( in ) :: PERM
        real, intent(in), dimension(:,:) :: viscosities
        REAL, DIMENSION( :, :, : ), optional, intent( inout ) :: inv_mat_absorp
        real, dimension(:,:,:), target, optional ::inv_perm1
        ! Local variable
        type (tensor_field), pointer :: RockFluidProp
        real, dimension(:), pointer :: Immobile_fraction, Corey_exponent, Endpoint_relperm
        REAL, PARAMETER :: TOLER = 1.E-10
        INTEGER :: ELE, CV_ILOC, CV_NOD, CV_PHA_NOD, MAT_NOD, JPHA_JDIM, &
            IPHA_IDIM, IDIM, JDIM, IPHASE, id_reg
        REAL, DIMENSION( :, :, :), pointer :: INV_PERM
        integer :: one_or_zero, visc_node
        real, dimension(Mdims%ndim*Mdims%nphase, Mdims%ndim*Mdims%nphase) :: Maux
        !Prepapre index for viscosity
        one_or_zero = (size(viscosities,2)==Mdims%cv_nonods)

        RockFluidProp=>extract_tensor_field(packed_state,"PackedRockFluidProp")
        ewrite(3,*) 'In calculate_absorption2'

        if (present(inv_perm1)) then
            INV_PERM => inv_perm1
        else
            ALLOCATE( INV_PERM(  Mdims%ndim, Mdims%ndim, Mdims%totele ))
            do id_reg = 1, size(perm,3)
                inv_perm( :, :, id_reg)=inverse(perm( :, :, id_reg))
            end do
        end if
        Maux = 0.
        DO ELE = 1, Mdims%totele
            !Get properties from packed state
            Immobile_fraction => RockFluidProp%val(1, :, ELE)
            Endpoint_relperm => RockFluidProp%val(2, :, ELE)
            Corey_exponent => RockFluidProp%val(3, :, ELE)
            DO CV_ILOC = 1, Mdims%cv_nloc
                MAT_NOD = ndgln%mat(( ELE - 1 ) * Mdims%mat_nloc + CV_ILOC)
                CV_NOD = ndgln%cv(( ELE - 1) * Mdims%cv_nloc + CV_ILOC )
                visc_node = (CV_NOD-1)*one_or_zero + 1
                DO IPHASE = 1, Mdims%n_in_pres
                    CV_PHA_NOD = CV_NOD + ( IPHASE - 1 ) * Mdims%cv_nonods
                    DO IDIM = 1, Mdims%ndim
                        IPHA_IDIM = ( IPHASE - 1 ) * Mdims%ndim + IDIM
                        DO JDIM = 1, Mdims%ndim
                            JPHA_JDIM = ( IPHASE - 1 ) * Mdims%ndim + JDIM
                            if (present(inv_mat_absorp)) then
                                call get_relperm(Mdims%n_in_pres, iphase, PorousMedia_absorp%val(idim, jdim, iphase, mat_nod),&
                                    SATURA(:, CV_NOD), viscosities(:,visc_node), INV_PERM( IDIM, JDIM, ELE),&
                                    Immobile_fraction, Corey_exponent, Endpoint_relperm, perm( IDIM, JDIM, ELE), inv_mat_absorp( IPHA_IDIM, JPHA_JDIM, MAT_NOD ))
                            else
                                call get_relperm(Mdims%n_in_pres, iphase, PorousMedia_absorp%val(idim, jdim, iphase, mat_nod),&
                                    SATURA(:, CV_NOD), viscosities(:,visc_node), INV_PERM( IDIM, JDIM, ELE),&
                                    Immobile_fraction, Corey_exponent, Endpoint_relperm)
                            end if
                        END DO
                    END DO
                END DO
            END DO
        END DO
        if (.not. present(inv_perm1)) DEALLOCATE( INV_PERM )
        ewrite(3,*) 'Leaving calculate_absorption2'
        RETURN
    END SUBROUTINE calculate_absorption2


    subroutine get_relperm(nphase, iphase, material_absorption, sat, visc, INV_PERM, Immobile_fraction, &
            Corey_exponent, Endpoint_relperm, PERM, inv_mat_absorp )
        !Calculates the relative permeability for 1, 2 (Brooks-corey) or 3 (stone's model) phases
        implicit none
        real, intent(inout) :: material_absorption
        real, intent(in) :: INV_PERM
        real, dimension(:), intent(in) :: sat, visc, Immobile_fraction, Corey_exponent, Endpoint_relperm
        integer, intent(in) :: iphase, nphase
        real, optional, intent(inout) :: inv_mat_absorp
        real, optional, intent(in) :: PERM
        !Local variables
        real, parameter :: epsilon = 1d-10!This value should in theory never be used, the real lower limit
        real, parameter :: eps = 1d-5!eps is another epsilon value, for less restrictive things

        select case (nphase)
            case (1)
                material_absorption = INV_PERM* visc(iphase) * min(1.0,max(eps,sat(iphase)))
                if (present(inv_mat_absorp).and.present(PERM)) &
                        inv_mat_absorp = PERM /(visc(iphase) * min(1.0,max(eps,sat(iphase))))
            case (3)
                call relperm_stone(material_absorption)
            case default
                call relperm_corey_epsilon(material_absorption)
        end select

        contains
            SUBROUTINE relperm_corey_epsilon( material_absorption )
                  !This subroutine add a small quantity to the corey function to avoid getting a relperm=0 that may give problems
                  !when dividing it to obtain the sigma.
                IMPLICIT NONE
                REAL, intent( inout ) :: material_absorption
                ! Local variables...
                REAL :: KR, aux
                !Kr_max should only multiply the wetting phase,
                !however as we do not know if it is phase 1 or 2, we let the decision to the user
                !and we multiply both phases by kr_max. By default kr_max= 1

                aux = 1.0 - sum(Immobile_fraction)
                if (present(inv_mat_absorp).and.present(PERM)) then
                    KR = Endpoint_relperm(iphase)*((sat(iphase) - Immobile_fraction(iphase)) / aux ) ** Corey_exponent(iphase)
                    inv_mat_absorp = (perm * max(0.0, KR))/(visc(iphase) * max(eps, sat(iphase)))
                end if
                KR = Endpoint_relperm(iphase)*( max( sat(iphase) - Immobile_fraction(iphase), sat(iphase)*eps+eps) / aux ) ** Corey_exponent(iphase)
                !Make sure that the relperm is between bounds
                KR = min(max(epsilon, KR),Endpoint_relperm(iphase))!Lower value just to make sure we do not divide by zero.
                material_absorption = INV_PERM * (visc(iphase) * max(eps, sat(iphase))) / KR !The value 1d-5 is only used if the boundaries have values of saturation of zero.
                  !Otherwise, the saturation should never be zero, since immobile fraction is always bigger than zero.
            END SUBROUTINE relperm_corey_epsilon

            subroutine relperm_stone(material_absorption)
                !This subroutine calculates the relative permeability for three phases
                !First phase has to be water, second oil and the third gas
                !We use Stone's model II adapted, and for the two phases we use the Corey model
                !Model explained in: Aziz, K. And Settari, T.:“Petroleum Reservoir Simulation” Applied Science Publishers, London, 30-38, 1979.
                implicit none
                real, intent(inout) :: material_absorption
                !Local variables
                real, dimension(3) :: Norm_sat, relperm, KR
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
                KR(iphase) = min(max(epsilon, relperm(iphase)),Endpoint_relperm(iphase))!Lower value just to make sure we do not divide by zero.
                material_absorption = INV_PERM * (VISC(iphase) * max(1d-5,sat(iphase))) / KR(iphase) !The value 1d-5 is only used if the boundaries have values of saturation of zero.
                !Otherwise, the saturation should never be zero, since immobile fraction is always bigger than zero.
                if (present(inv_mat_absorp).and.present(PERM)) &!This part is to ensure that the flow is stopped
                            inv_mat_absorp = (perm * max(0.0,relperm(iphase)))/(VISC(iphase) * max(eps,sat(iphase)))
            end subroutine relperm_stone

    end subroutine get_relperm


    SUBROUTINE calculate_capillary_pressure( packed_state, &
        NDGLN, totele, cv_nloc)

            ! CAPIL_PRES_OPT is the capillary pressure option for deciding what form it might take.
            ! CAPIL_PRES_COEF( NCAPIL_PRES_COEF, NPHASE, NPHASE ) are the coefficients
            ! Capillary pressure coefs have the dims CAPIL_PRES_COEF( NCAPIL_PRES_COEF, NPHASE,NPHASE )
            ! used to calculate the capillary pressure.

            IMPLICIT NONE
            type(state_type), intent(inout) :: packed_state
            integer, intent(in) :: totele, cv_nloc
            type(multi_ndgln), intent(in) :: ndgln
            ! Local Variables
            INTEGER :: IPHASE, JPHASE, nphase, ele, cv_iloc, cv_nod
            logical, save :: Cap_Brooks = .true., Cap_Power = .false.
            logical, save :: first_time = .true.
            !Working pointers
            real, dimension(:,:), pointer :: Satura, CapPressure, Immobile_fraction, Cap_entry_pressure, Cap_exponent, Imbibition_term
            real, dimension(:), allocatable :: Cont_correction
            !Get from packed_state
            call get_var_from_packed_state(packed_state,PhaseVolumeFraction = Satura)

            call get_var_from_packed_state(packed_state,CapPressure = CapPressure, &
                Immobile_fraction = Immobile_fraction, Cap_entry_pressure = Cap_entry_pressure, &
                    Cap_exponent = Cap_exponent, Imbibition_term = Imbibition_term)

            nphase =size(Satura,1)
            allocate(Cont_correction(size(satura,2)))

            CapPressure = 0.
            ! Determine which capillary pressure model is to be used for overrelaxation. Use Brooks-Corey unless power_law Pc activated (important to allow overelax even when Pc is off).
            if (first_time) then
                Cap_Power = have_option_for_any_phase("/multiphase_properties/capillary_pressure/type_Power_Law", nphase)
                Cap_Brooks = .not. (Cap_Power)

                first_time = .false.
            end if

            DO IPHASE = 1, NPHASE

                if ( (Cap_Brooks) .or. (Cap_Power) ) then

                    !Apply Capillary model
                    do jphase = 1, nphase
                        Cont_correction = 0
                        if (jphase /= iphase) then!Don't know how this will work for more than 2 phases
                            do ele = 1, totele
                                do cv_iloc = 1, cv_nloc
                                    cv_nod = ndgln%cv((ele-1)*cv_nloc + cv_iloc)
                                    CapPressure( jphase, cv_nod ) = CapPressure( jphase, cv_nod ) + &
                                        Get_capPressure(satura(iphase,cv_nod), Cap_entry_pressure(iphase, ele), &
                                        Cap_exponent(iphase, ele),Immobile_fraction(:,ele), &
                                        Imbibition_term(iphase, ele), iphase)
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
        contains
            pure real function Get_capPressure(sat, Pe, a, Immobile_fraction, Imbibition_term, iphase)
                !This functions returns the capillary pressure for a certain input saturation
                !There is another function, its derivative in cv-adv-diff called Get_DevCapPressure
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

    real function Get_DevCapPressure(sat, Pe, a, immobile_fraction, iphase, nphase)
        !This functions returns the derivative of the capillary pressure with respect to the saturation
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
            u_source_cv = 0.
            do nod = 1, cv_nonods
                g = node_val( gravity_direction, nod ) * gravity_magnitude
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

    subroutine calculate_diffusivity(state, Mdims, ndgln, ScalarAdvectionField_Diffusion, tracer)
      type(state_type), dimension(:), intent(in) :: state
      type(multi_dimensions), intent(in) :: Mdims
      type(multi_ndgln), intent(in) :: ndgln
      real, dimension(:, :, :, :), intent(inout) :: ScalarAdvectionField_Diffusion
      !Local variables
      type(scalar_field), pointer :: component, sfield
      type(tensor_field), pointer :: diffusivity, tfield
      integer :: icomp, iphase, idim, stat, ele
      integer :: iloc, mat_inod, cv_inod, ele_nod, t_ele_nod
      logical, parameter :: harmonic_average=.false.
      type(tensor_field), intent(inout) :: tracer

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

        diffusivity => extract_tensor_field( state(1), 'TemperatureDiffusivity', stat )
        !Note that for the temperature field this is actually the thermal conductivity (in S.I. watts per meter-kelvin => W/(m·K) ).
        if ( stat == 0 ) then

            if (is_porous_media) then
                sfield=>extract_scalar_field(state(1),"Porosity")
                tfield => extract_tensor_field( state(1), 'porous_thermal_conductivity', stat )
                ScalarAdvectionField_Diffusion = 0.
                ! Calculation of the averaged thermal diffusivity as
                ! lambda = porosity * lambda_f + (1-porosity) * lambda_p
                ! Since lambda_p is defined element-wise and lambda_f CV-wise we perform an average
                ! as it is stored cv-wise
                ! NOTE: that we are considering a unified lambda for all the phases
                do iphase = 1, Mdims%nphase
                    diffusivity => extract_tensor_field( state(iphase), 'TemperatureDiffusivity', stat )
                    do ele = 1, Mdims%totele
                        ele_nod = min(size(sfield%val), ele)
                        t_ele_nod = min(size(tfield%val, 3), ele)
                         do iloc = 1, Mdims%mat_nloc
                            mat_inod = ndgln%mat( (ele-1)*Mdims%mat_nloc + iloc )
                            cv_inod = ndgln%cv((ele-1)*Mdims%cv_nloc+iloc)
                            do idim = 1, Mdims%ndim
                                ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase ) = &
                                    ScalarAdvectionField_Diffusion( mat_inod, idim, idim, iphase )+&
                                    (sfield%val(ele_nod) * node_val( diffusivity, idim, idim, mat_inod ) &
                                    +(1.0-sfield%val(ele_nod))* tfield%val(idim, idim, t_ele_nod))
                            end do
                        end do
                    end do
                end do

            else
                do iphase = 1, Mdims%nphase
                    diffusivity => extract_tensor_field( state(iphase), 'TemperatureDiffusivity', stat )
                    do idim = 1, Mdims%ndim
                        ScalarAdvectionField_Diffusion( :, idim, idim, iphase ) = node_val( diffusivity, idim, idim, iphase )
                    end do
                end do
            end if
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
         ewrite(3,*) 'Thermal conductivity min_max', iphase, &
              minval( ScalarAdvectionField_Diffusion( :, 1, 1, iphase ) ), &
              maxval( ScalarAdvectionField_Diffusion( :, 1, 1, iphase ) )
      end do
      return
    end subroutine calculate_diffusivity

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
      type( scalar_field ), pointer :: component
      logical :: linearise_viscosity
      real, dimension( : ), allocatable :: component_tmp
      real, dimension( :, :, : ), allocatable :: mu_tmp
      integer :: iloc, ndim1, ndim2, idim, jdim
      integer :: multiplier


  ! DELETE Momentum_Diffusion - START USING THE NEW MEMORY ---

      if ( is_porous_media) then
         momentum_diffusion=0.0
      else
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
                           cv_nod = cv_nod * multiplier + (1 - multiplier)!index has to be one if viscosity is constant
                           mat_nod = ndgln%mat( (ele-1)*Mdims%cv_nloc + iloc )
                           momentum_diffusion( :, :, iphase, mat_nod ) = momentum_diffusion(  :, :, iphase, mat_nod ) + mu_tmp( 1, 1, iloc ) ! isotropic only - to be deleted...
                           t_field%val( :, :, cv_nod ) = t_field%val( :, :, cv_nod ) + mu_tmp( :, :, iloc )/dble(Mdims%cv_nloc)
                        end do
                     end do
                  end do
               end do
            else
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
                        mat_nod = mat_nod * multiplier + (1 - multiplier)!index has to be one if viscosity is constant
                        momentum_diffusion( :, :, iphase, mat_nod ) = mu_tmp( :, :, iloc )
                        t_field%val( :, :, mat_nod ) = mu_tmp( :, :, iloc )
                     end do
                  end do
               end do
            end if
            deallocate( component_tmp, mu_tmp )
         end if
      end if


      !!! NEW CODE HERE !!!
      !!! deal with Momentum_Diffusion2

      if ( is_porous_media ) then
         return
      else

         ! return here as code below untested
         return

         t_field => extract_tensor_field( state( 1 ), "Viscosity", stat ) ! need to set dimensions in diamond - Populate_State.F90:2164
         if ( stat == 0 ) then

            do iphase = 1, Mdims%nphase

               call allocate_multi_field( state( iphase ), Mdims, "Viscosity", Momentum_Diffusion2 )

               if ( Mdims%ncomp > 1 ) then

                  tp_field => extract_tensor_field( state( iphase ), "Viscosity" )
                  call zero( tp_field )
                  ndim1 = size( tp_field%val, 1 ) ; ndim2 = size( tp_field%val, 2 ) 

                  do icomp = 1, Mdims%ncomp
                     component => extract_scalar_field( state( Mdims%nphase + icomp ), "ComponentMassFractionPhase" // int2str( iphase ) )
                     tc_field => extract_tensor_field( state( Mdims%nphase + icomp ), "Viscosity" )
                     do ele = 1, Mdims%totele
                        do iloc = 1, Mdims%mat_nloc
                           mat_nod = ndgln%mat( (ele-1)*Mdims%cv_nloc + iloc )
                           cv_nod = ndgln%cv( (ele-1)*Mdims%cv_nloc + iloc )
                           call addto( tp_field, mat_nod, node_val( component, cv_nod ) * node_val( tc_field, mat_nod ) )
                        end do
                     end do
                  end do

               end if

            end do

            if ( have_option( "/material_phase[0]/linearise_viscosity" ) ) then
               call linearise_multi_field( Momentum_Diffusion2, Mdims, ndgln%mat )
            end if

         end if
      end if

      !!!!!!!!!!!!!!!!!!!!!

      return
    end subroutine calculate_viscosity



    !sprint_to_do, re-use material_absoprtion by updating the values of the input absoprtion
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

!        if (any(have_source)) then
!            u_source%have_field = .true.
!            if (associated(u_source%val))then
!                nullify(u_source%val)
!                deallocate(u_source%val)
!            end if
!            allocate(u_source%val(Mdims%ndim, Mdims%nphase, 1, Mdims%u_nonods))
!        end if

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
        character(len=500) :: path, path2, path3

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
                "]/multiphase_properties/Relperm_Corey/scalar_field::relperm_max/prescribed/value"
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
                "]/multiphase_properties/Relperm_Corey/scalar_field::relperm_exponent/prescribed/value"
            if (have_option(trim(path))) then
                call initialise_field_over_regions(targ_Store, trim(path) , position)
                t_field%val(3,iphase,:) = targ_Store%val(:)
            else!default value
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
                    t_field%val(4,iphase,:) = targ_Store%val(:)
                elseif (have_option(trim(path3))) then
                    call initialise_field_over_regions(targ_Store, trim(path3) , position)
                    t_field%val(4,iphase,:) = targ_Store%val(:)
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
                    t_field%val(5,iphase,:) = targ_Store%val(:)
                elseif (have_option(trim(path3))) then
                    call initialise_field_over_regions(targ_Store, trim(path3) , position)
                    t_field%val(5,iphase,:) = targ_Store%val(:)
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
                    t_field%val(6,iphase,:) = targ_Store%val(:)
                else !default value
                    t_field%val(6,iphase,:) = 0.0
                end if
            end do
        end if


        call deallocate(targ_Store)
    end subroutine get_RockFluidProp


    subroutine get_FloodingProp(state, packed_state)
        !Gets flooding manning coefficient
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

        t_field=>extract_tensor_field(packed_state,"PackedManningcoef")
        !By default the pressure mesh (position 1)
        s_field => extract_scalar_field(state(1),1)
        position => get_external_coordinate_field(packed_state, s_field%mesh)

        fl_mesh => extract_mesh( state(1), "P0DG" )
        Auxmesh = fl_mesh
        call allocate (targ_Store, Auxmesh, "Temporary_Manningcoef")

        !Retrieve Manning coefficent
        path = "/flooding/scalar_field::manning_coef/prescribed/value"
        if (have_option(trim(path))) then
            call initialise_field_over_regions(targ_Store, trim(path) , position)
            t_field%val(1,1,:) = targ_Store%val(:)
        end if
        call deallocate(targ_Store)


        fl_mesh => extract_mesh( state(1), "PressureMesh" )
        Auxmesh = fl_mesh
        call allocate (targ_Store, Auxmesh, "Temporary_Bathymetry")
        !Retrieve Bathymetry coefficent
        path = "/material_phase[0]/scalar_field::Bathymetry/prescribed/value"
        if (have_option(trim(path))) then
            t_field=>extract_tensor_field(packed_state,"PackedBathymetry")
            call initialise_field_over_regions(targ_Store, trim(path) , position)
            t_field%val(1,1,:) = targ_Store%val(:)
        end if

        call deallocate(targ_Store)
     end subroutine get_FloodingProp




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

    !!-JWL eqaution functions






    !###############################################################################################
    !This subroutine currently works in the oil industry standard units. However, the rest of the code works in SI units
    !Appropriate conversions are performed internally to make everything consistent.
    !The reason for not using the SI is because all the PVT tables are formulas are not in SI.
    !###############################################################################################
    subroutine extended_Black_Oil(state, packed_state, Mdims, flash_flag, viscosities, Kcomp, icomp)
        !In this subroutine, phase change, density and viscosity are calculated based on PVT tables and locally per cell.
        !This method is called extended Black-Oil. Three phases are considered. Phase 1 is aqua and is inert
        !The other two phases are liquid and vapour. Two pseudo components associated with each of this two phases are considered
        !This approach, calculates the gas equilibrium per CV and then a normal three phase simulation is run. Hence, the
        !saturation update has to occur before the call to the non-linear solver to allow the phases to reach equilibrium again.
        !According to the sources, this requires an extra Newton iteration, in our case the extra cost is also very low, around 2 more.
        !This method not only allows for black-oil but also for Low-volatile and condensates
        !Sources: SPE15156, petrowiki.org/Oil_fluid_characteristics, and DOI: 10.2516/ogst/2012001
        implicit none
        type( state_type ), dimension(:), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent(in) :: Mdims
        integer, intent(in) :: flash_flag!Selects which field to update 0 => All; 1 => Density; 2 => Saturation; 3 => viscosity
        real, optional, dimension(:,:), intent(inout) :: viscosities
        real, optional, dimension(:,:,:), intent(inout) :: Kcomp! cv_nonods, nphase,nphase
        integer, optional, intent(in) :: icomp
        !Local variables                 !4 => Saturation and viscosity;
        type( tensor_field ), pointer :: pressure, saturation, density, visc_liquid, visc_vapour,&
                    DRhoDP, visc_aqua
        type( scalar_field ), pointer :: VapourMassFraction
        real :: RhoPlus, RhoMinus
        integer :: cv_inod, iphase, k
        real, dimension(3) :: density_reference
        real, dimension(:,:), allocatable :: PVT_table
        real, dimension(:), allocatable :: perturbation_pressure
        character( len = option_path_len ) :: PVT_table_path
        character( len = option_path_len ), parameter :: option_path = "/physical_parameters/black-oil_PVT_table"
        real, pointer :: Zg
        real :: Yg,Yo,Xg,Xo, Zo, rho_stc_ratio, rho_ratio,Sg,&
            Kg, Ko, aux, aux2, fv, ro_cap, Mix_o, Mix_g, ro_min_cap

        !Default molecular weights for gas and oil respectively, just in case they are not defined elsewhere
        real :: Mg = 22.9395, Mo = 190.0
        !Black-Oil option
        logical :: Incompressible, constant_visc
        real, parameter :: Ref_press = 1e5
        !Conversion parameters
        real, parameter :: density2SI = 16.018463374!pound/feet^3 to SI (Kg/m^3)
        real, parameter :: Pressure2SI = 6894.75729!psi to SI (Pascals)
        real, parameter :: Viscosity2SI = 0.001!cp to SI (Pascals * secondws)


        if (Mdims%nphase /=3 ) then
            !Phase 1 => aqua phase; no interaction with other phases.
            !Phase 2 => liquid phase
            !Phase 3 => vapour phase
            ewrite(1,*) "WARNING: extended_Black_Oil requires exactly 3 phases defined"
            return
        end if


        !Check physics required by the user
        Incompressible = have_option(trim(option_path)//"/incompressible")
        constant_visc  = have_option(trim(option_path)//"/constant_visc")

        !Retrieve path to the .csv file with the data
        call get_option(trim(option_path), PVT_table_path)
        !Manual table for the time being
        if (trim(PVT_table_path)=="internal") then
!            call populate_with_Texas_Black_Oil(PVT_table, density_reference)
            call populate_with_Black_Oil(PVT_table, density_reference)
        else
            !READ FROM FILE
            call read_csv_table(PVT_table, PVT_table_path, density_reference)
            !If it has 8 entries, means that the volatize-oil/gas ratio is significant (it is in posiiton 8), otherwise ignored
        end if

        pressure => extract_tensor_field(packed_state,"PackedFEPressure")
        !Calculate the molecular density ratio (oil/gas)
        rho_stc_ratio = density_reference(2)/density_reference(3)* (Mg/Mo)

        if (flash_flag == 10) then
            !Initialize the MassFraction of the Gas component
            saturation => extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
            VapourMassFraction => extract_scalar_field(state(3),"VapourMassFraction")
            do cv_inod = 1, mdims%cv_nonods
                Zg => VapourMassFraction%val(cv_inod)
                if (saturation%val(1,1,cv_inod)==1.) then
                    Sg = 0.; Zg=0.
                else
                    Sg = min(max(saturation%val(1,3,cv_inod)/(1.0-saturation%val(1,1,cv_inod))  ,0.),1.)
                    !Calculate the component phase change K constants and the partial mass fractions
                    call calculate_K_comp(ko,kg, xg,xo,yg,yo, pressure%val(1,1,cv_inod))
                    !Calculate mix molecular weights
                    Mix_o = xg * Mg + xo * Mo; Mix_g = yo * Mo + yg * Mg
                    !Calculate rho_ratio_molecular(rho_gas/rho_oil)
                    rho_ratio = (get_density(pressure%val(1,1,cv_inod), PVT_table, density_reference, 3)/Mix_g)/&
                        (get_density(pressure%val(1,1,cv_inod), PVT_table, density_reference, 2)/Mix_o)
                    !Initialize mass fraction, based on the initial pressure and saturations
                    aux = 1-kg; aux2 = 1-Ko
                    Zg = gas_mass_fraction_from_gas_saturation(Sg, Kg, Ko, rho_ratio)
                end if
            end do
            return
        end if


        !Adjust densities based on the PVT table and reference densities
        if (flash_flag == 0 .or. flash_flag == 1) then
            allocate(perturbation_pressure(Mdims%cv_nonods))
            perturbation_pressure = max( 1.e-10, 1.e-3 * abs( pressure%val(1,1,:) ) )
            density => extract_tensor_field( packed_state,"PackedDensity" )
            DRhoDP => extract_tensor_field( packed_state, "PackedDRhoDPressure" )

            !Phase 1 is incompressible for the time being, when is not, remove this and change loop to account for phase 1 as well
            density%val(1,1,:) = density_reference(1) * density2SI
            DRhoDP%val(1,1,:)  = 0

            if (Incompressible) then
                do iphase = 2, Mdims%nphase
                    density%val(1,iphase,:) = density_reference(iphase)*density2SI
                    !Force DRhoDP =0 to all phases to make it incompressible
                    DRhoDP%val(1,iphase,:) = 0.
                end do
            else
                do cv_inod = 1, mdims%cv_nonods
                    !Calculate densities and derivatives of the densities
                    do iphase = 2, Mdims%nphase
                        density%val(1,iphase,cv_inod) = get_density(pressure%val(1,1,cv_inod), PVT_table, density_reference, iphase)*density2SI
                        RhoPlus = get_density(pressure%val(1,1,cv_inod) + perturbation_pressure(cv_inod), PVT_table, density_reference, iphase)*density2SI
                        RhoMinus= get_density(pressure%val(1,1,cv_inod) - perturbation_pressure(cv_inod), PVT_table, density_reference, iphase)*density2SI
                        DRhoDP%val(1,iphase,cv_inod) = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure(cv_inod)
                    end do
                end do
            end if

            deallocate(perturbation_pressure)
        end if

        if (flash_flag == 0 .or. flash_flag == 2 .or. flash_flag == 4) then
            saturation => extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
            VapourMassFraction => extract_scalar_field(state(3),"VapourMassFraction")

            do cv_inod = 1, mdims%cv_nonods
                Zg => VapourMassFraction%val(cv_inod);Zo = 1- Zg
                !Original paper, for Black-Oil only
                !Calculate the K components and the partial mass fractions
                call calculate_K_comp(ko,kg, xg,xo,yg,yo, pressure%val(1,1,cv_inod))
                !Check if gas phase exists
                if (ko*Zo + Kg*Zg /= 1.) then
                    fv = min(max(Zo / (1. - Kg)  + Zg/(1 - Ko),0.),1.)

                    !Calculate mix molecular weights
                    Mix_o = xg * Mg + xo * Mo; Mix_g = yo * Mo + yg * Mg

                    !Calculate rho_ratio_molecular(rho_gas/rho_oil)
                    rho_ratio = (get_density(pressure%val(1,1,cv_inod), PVT_table, density_reference, 3)/Mix_g)/&
                        (get_density(pressure%val(1,1,cv_inod), PVT_table, density_reference, 2)/Mix_o)
                    Sg = min(max(1.0/(1.+(1./fv - 1.)  *  rho_ratio),0.),1.)
                    saturation%val(1,3,cv_inod) = (1.0-saturation%val(1,1,cv_inod)) * Sg
                else!No gas phase
                    saturation%val(1,3,cv_inod) = 0.
                end if
                saturation%val(1,2,cv_inod) = 1. - saturation%val(1,1,cv_inod) - saturation%val(1,3,cv_inod)
            end do

        end if

        if (flash_flag == 0 .or. flash_flag == 3 .or. flash_flag == 4) then

            if (present(viscosities)) then!Modify the introduced viscosity data, not the one in packed state
                !Aqua phase set from mpml(SI units)
                visc_aqua => extract_tensor_field(state(1),"Viscosity")
                viscosities(1,:) = visc_aqua%val(1,1,1)!Only isotropic viscosity
                if(constant_visc) then
                    !Set from mpml (SI units)
                    visc_liquid => extract_tensor_field(state(2),"Viscosity")
                    visc_vapour => extract_tensor_field(state(3),"Viscosity")
                    viscosities(2,:) = visc_liquid%val(1,1,1)!Only isotropic viscosity
                    viscosities(3,:) = visc_vapour%val(1,1,1)!Only isotropic viscosity
                else
                    do cv_inod = 1, mdims%cv_nonods
                        !Obtain new viscosities(just diagonal tensor)
                        do k = 1, Mdims%ndim
                            viscosities(2,cv_inod) = eval_table(pressure%val(1,1,cv_inod), PVT_table,5)*Viscosity2SI
                            viscosities(3,cv_inod) = eval_table(pressure%val(1,1,cv_inod), PVT_table,6)*Viscosity2SI
                        end do
                    end do
                end if
            else!this part below is not in use nor tested
                visc_liquid => extract_tensor_field(state(2),"Viscosity")
                if (size(visc_liquid%val,3) > 1) then!Update only if viscosity is defined per CV, not homogenenous
                    visc_vapour => extract_tensor_field(state(3),"Viscosity")
                    if(.not.constant_visc) then
                        do cv_inod = 1, mdims%cv_nonods
                            !Have to change Calculate_PorousMedia_AbsorptionTerms to consider non-homogenenous viscosities
                            !Obtain new viscosities(just diagonal tensor)
                            do k = 1, Mdims%ndim
                                visc_liquid%val(k,k,cv_inod) = eval_table(pressure%val(1,1,cv_inod), PVT_table,5)*Viscosity2SI
                                visc_vapour%val(k,k,cv_inod) = eval_table(pressure%val(1,1,cv_inod), PVT_table,6)*Viscosity2SI
                            end do
                        end do
                    end if
                end if
            end if
        end if


        if (flash_flag == 6 .and.present(Kcomp).and.present(icomp)) then
            saturation => extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
            VapourMassFraction => extract_scalar_field(state(3),"VapourMassFraction")

            !Conversion from PVT data to compositional

            do cv_inod = 1, mdims%cv_nonods
                !Original paper, for Black-Oil only
                !Calculate the K components and the partial mass fractions
                call calculate_K_comp(ko,kg, xg,xo,yg,yo, pressure%val(1,1,cv_inod))
                !Store K component for the compositional
                select case (icomp)
                    case (1)
                        Kcomp(cv_inod,:,:) = 1.
                    case (2)
                        Kcomp(cv_inod,:,:) = Ko
                    case (3)
                        Kcomp(cv_inod,:,:) = Kg
                end select
            end do

        end if
        deallocate(PVT_table)

    contains

        real function gas_mass_fraction_from_gas_saturation(Sg, Kg, Ko, rho_ratio)
            !this is ignoring the
            implicit none
            real, intent(in) ::Sg, Kg, Ko, rho_ratio
            !Local variables
            real :: aux, aux2
            aux = 1-kg; aux2 = 1-Ko
            gas_mass_fraction_from_gas_saturation = ((Sg*rho_ratio*aux*aux2)/(1-Sg + Sg*rho_ratio) - aux2)/(aux-aux2)

        end function gas_mass_fraction_from_gas_saturation

        subroutine calculate_K_comp(ko,kg, xg,xo,yg,yo, Pres)
            !Calculates the K between the pseudo components
            !Requires rho_stc_ratio to have been defined before
            implicit none
            real, intent(inout) :: ko,kg, xg,xo,yg,yo, Pres
            !Local variables
            real ::ro_min_cap, RO_CAP

            RO_CAP = 0.178108 * eval_table(Pres, PVT_table,4) / rho_stc_ratio
            ro_min_cap = 5.61458d-6 * eval_table(Pres, PVT_table,8) * rho_stc_ratio
            !Calculate partial compositions of the phases(X is for liquid, and x for vapour)
            Xg = RO_CAP/(1+RO_CAP); Xo = 1. - Xg
            Yg = 1./(1+ro_min_cap); Yo = 1. - Yg
            !Calculate K component
            Kg = Yg/max(Xg, 1e-10)
            Ko = Yo/max(Xo, 1e-10)

        end subroutine calculate_K_comp

        real function get_density(pressure, PVT_table, density_reference, phase)
            !density_formation = (density_phase + density_other_phase * disolved_of_other_phase)*phase_volumetric_factor
            implicit none
            real, intent(in) :: pressure
            real, dimension(:,:), intent(in) :: PVT_table
            real, dimension(3), intent(in) :: density_reference
            integer, intent(in) :: phase
            select case (phase)
                case (2)
                    get_density = (density_reference(2) + density_reference(3)*&
                    eval_table(pressure, PVT_table,4)/5.61458)/eval_table(pressure, PVT_table,2)
                case (3)!still incomplete for phase 3
                    get_density = (density_reference(3) + 5.61458d-6 * density_reference(2)*&
                    eval_table(pressure, PVT_table,8))/(eval_table(pressure, PVT_table,3))!/(1+ro))
                case default
                    get_density = density_reference(phase)!/eval_table(pressure, PVT_table,9)
            end select

        end function get_density

        real function eval_table(input_pressure, PVT_table, prop)
            implicit none
            integer, intent(in) :: prop!Which property are you interested in
            real, intent(in) :: input_pressure
            real, dimension(:,:), intent(in) :: PVT_table
            !Local variables
            real, parameter :: tol = 1e-5
            integer :: i, k
            integer, parameter :: interpolation_order = 1!Either linear (1) or quadratic (2)
            real, dimension((interpolation_order + 1)) :: P_pos, prop_pos
            integer, save :: bubble_point_pos = -1
            real :: P_in_psi
            !Tables are in PSI, convert pressure to PSI
            P_in_psi = input_pressure/Pressure2SI

            !Need to change this to one that does not depend on the 7 column
            if (bubble_point_pos < 0) then
                !detect and store bubble point position
                do bubble_point_pos = 1, size(PVT_table,2)
                    if (PVT_table(bubble_point_pos,7) < 0.999999) exit
                end do
                bubble_point_pos = bubble_point_pos - 1
            end if

            !Check position of the input pressure in the table
            !Considering highest pressure in possition 1
            do i = 1, size(PVT_table,2)
                if (PVT_table(1,i) + tol <= P_in_psi) exit
            end do

            if (i == 1) then
                eval_table = PVT_table(prop,1)
            else if (i >= size(PVT_table,2)) then
                eval_table = PVT_table(prop,size(PVT_table,2))
            else
                !Detect if we are in bubble point
                if (i-1 == bubble_point_pos) then
                    !Shift the position ensuring that all the values lie after the bubble point
                    i = max(bubble_point_pos, 1)
                else if (i <= bubble_point_pos) then
                    !Shift the position ensuring that all the values lie before the bubble point
                    i = max(i - (interpolation_order + 1), 1)
                else
                    !Shift the position to make sure we use values surrounding the input pressure
                    i = max(i - interpolation_order, 1)
                end if

                do k = 0, interpolation_order
                    P_pos(k+1) = PVT_table(1,i+k)
                    prop_pos(k+1) = PVT_table(prop,i+k)
                end do
                eval_table = table_interpolation(P_pos, prop_pos, P_in_psi)
            end if

        end function eval_table

        subroutine populate_with_Black_Oil(PVT_table, density_reference)
            implicit none!Nothing in SI but the pressure
            real, dimension(:,:), allocatable, intent(inout) :: PVT_table
            real, dimension(3), intent(inout) :: density_reference

            if(allocated(PVT_table)) deallocate(PVT_table)
            allocate(PVT_table(8,9))
            !Aqua density_reference
            density_reference(1) = 62.427960576;
            !Liquid density_reference ; !Vapour density reference (RB/STB)
            density_reference(2) = 49.111; density_reference(3) = 0.06044
            !Molecular weights of gas and oil
            Mg = 22.9395; Mo = 190.0

            PVT_table(1,1)  = 9000.; PVT_table(2,1)  = 2.357
            PVT_table(1,2)  = 5000.; PVT_table(2,2)  = 1.827
            PVT_table(1,3)  = 4000.; PVT_table(2,3)  = 1.695
            PVT_table(1,4)  = 3000.; PVT_table(2,4)  = 1.565
            PVT_table(1,5)  = 2500.; PVT_table(2,5)  = 1.500
            PVT_table(1,6)  = 2000.; PVT_table(2,6)  = 1.435
            PVT_table(1,7)  = 1000. ; PVT_table(2,7)  = 1.295
            PVT_table(1,8)  = 500. ; PVT_table(2,8)  = 1.207
            PVT_table(1,9)  = 250. ; PVT_table(2,9)  = 1.150
            !#Volumetric factor vapour(RB/STB)#!#####Dissolved GOR (scf/STB)(Gas disolved in oil -Rs-)######
            PVT_table(3,1)  = 0.002167; PVT_table(4,1)  = 2984
            PVT_table(3,2)  = 0.003644; PVT_table(4,2)  = 1618
            PVT_table(3,3)  = 0.004553 ; PVT_table(4,3)  = 1270
            PVT_table(3,4)  = 0.006064 ; PVT_table(4,4)  = 930
            PVT_table(3,5)  = 0.007265 ; PVT_table(4,5)  = 775
            PVT_table(3,6)  = 0.009062 ; PVT_table(4,6)  = 636
            PVT_table(3,7)  = 0.017950 ; PVT_table(4,7)  = 371
            PVT_table(3,8)  = 0.035226 ; PVT_table(4,8)  = 180
            PVT_table(3,9)  = 0.067897 ; PVT_table(4,9)  = 90
            !NEED TO FILL THIS WITH THE INFO FROM Practical Techniques in 1\vo-Pseudocomponent Black-Oil Simulation,
            !G.D. Shank and C.R. Vestal, SPE, Marathon Oil Co
            !#####Viscosity liquid####!#####Viscosity vapour###
            PVT_table(5,1)  = 2.03d-1; PVT_table(6,1)  = 4.70d-2
            PVT_table(5,2)  = 4.49d-1; PVT_table(6,2)  = 3.09d-2
            PVT_table(5,3)  = 5.10d-1; PVT_table(6,3)  = 2.68d-2
            PVT_table(5,4)  = 5.94d-1; PVT_table(6,4)  = 2.28d-2
            PVT_table(5,5)  = 6.41d-1; PVT_table(6,5)  = 2.08d-2
            PVT_table(5,6)  = 6.95d-1; PVT_table(6,6)  = 1.89d-2
            PVT_table(5,7)  = 8.30d-1; PVT_table(6,7)  = 1.40d-2
            PVT_table(5,8)  = 9.10d-1; PVT_table(6,8)  = 1.12d-2
            PVT_table(5,9)  = 9.75d-1; PVT_table(6,9)  = 0.96d-2
            !#####Dissolved GOR (Oil disolved in gas -Rv-)######
            PVT_table(8,:)  = 0.

        end subroutine





        subroutine populate_with_Texas_Black_Oil(PVT_table, density_reference)
            implicit none!Everything in the S.I. but the Rs amd Rv, for the time being
            real, dimension(:,:), allocatable, intent(inout) :: PVT_table
            real, dimension(3), intent(inout) :: density_reference

            if(allocated(PVT_table)) deallocate(PVT_table)
            allocate(PVT_table(8,12))
            !Aqua density_reference
            density_reference(1) = 62.427960576;
            !Liquid density_reference ; !Vapour density reference (RB/STB)
            density_reference(2) = 49.111; density_reference(3) = 0.06044
            !Molecular weights of gas and oil
            Mg = 22.9395; Mo = 190.0
            !#####Pressure####        !#Volumetric factor liquid#
            PVT_table(1,1)  = 2000.; PVT_table(2,1)  = 1.467
            PVT_table(1,2)  = 1800.; PVT_table(2,2)  = 1.472
            PVT_table(1,3)  = 1700.; PVT_table(2,3)  = 1.475
            PVT_table(1,4)  = 1640.; PVT_table(2,4)  = 1.463
            PVT_table(1,5)  = 1600.; PVT_table(2,5)  = 1.453
            PVT_table(1,6)  = 1400.; PVT_table(2,6)  = 1.408
            PVT_table(1,7)  = 1200.; PVT_table(2,7)  = 1.359
            PVT_table(1,8)  = 1000.; PVT_table(2,8)  = 1.322
            PVT_table(1,9)  = 800. ; PVT_table(2,9)  = 1.278
            PVT_table(1,10) = 600. ; PVT_table(2,10) = 1.237
            PVT_table(1,11) = 400. ; PVT_table(2,11) = 1.194
            PVT_table(1,12) = 200. ; PVT_table(2,12) = 1.141
            !#Volumetric factor vapour#!#####Dissolved GOR (Gas disolved in oil -Rs-)######
            PVT_table(3,1)  = 1.920 ; PVT_table(4,1)  = 838.5
            PVT_table(3,2)  = 1.920 ; PVT_table(4,2)  = 838.5
            PVT_table(3,3)  = 1.920 ; PVT_table(4,3)  = 838.5
            PVT_table(3,4)  = 1.920 ; PVT_table(4,4)  = 816.1
            PVT_table(3,5)  = 1.977 ; PVT_table(4,5)  = 798.4
            PVT_table(3,6)  = 2.308 ; PVT_table(4,6)  = 713.4
            PVT_table(3,7)  = 2.730 ; PVT_table(4,7)  = 621.0
            PVT_table(3,8)  = 3.328 ; PVT_table(4,8)  = 548.0
            PVT_table(3,9)  = 4.163 ; PVT_table(4,9)  = 464.0
            PVT_table(3,10) = 5.471 ; PVT_table(4,10) = 383.9
            PVT_table(3,11) = 7.786 ; PVT_table(4,11) = 297.4
            PVT_table(3,12) = 13.331; PVT_table(4,12) = 190.9
            !#####Viscosity liquid####!#####Viscosity vapour###  !#####Volume liquid Fraction######
            PVT_table(5,1)  = 3.201d-1; PVT_table(6,1)  = 1.57d-2;PVT_table(7,1)  = 1.0
            PVT_table(5,2)  = 3.114d-1; PVT_table(6,2)  = 1.57d-2;PVT_table(7,2)  = 1.0
            PVT_table(5,3)  = 3.071d-1; PVT_table(6,3)  = 1.57d-2;PVT_table(7,3)  = 1.0
            PVT_table(5,4)  = 3.123d-1; PVT_table(6,4)  = 1.57d-2;PVT_table(7,4)  = 0.978
            PVT_table(5,5)  = 3.169d-1; PVT_table(6,5)  = 1.55d-2;PVT_table(7,5)  = 0.96
            PVT_table(5,6)  = 3.407d-1; PVT_table(6,6)  = 1.40d-2;PVT_table(7,6)  = 0.867
            PVT_table(5,7)  = 3.714d-1; PVT_table(6,7)  = 1.38d-2;PVT_table(7,7)  = 0.754
            PVT_table(5,8)  = 3.973d-1; PVT_table(6,8)  = 1.32d-2;PVT_table(7,8)  = 0.644
            PVT_table(5,9)  = 4.329d-1; PVT_table(6,9)  = 1.26d-2;PVT_table(7,9)  = 0.513
            PVT_table(5,10) = 4.712d-1; PVT_table(6,10) = 1.21d-2;PVT_table(7,10)  = 0.375
            PVT_table(5,11) = 5.189d-1; PVT_table(6,11) = 1.16d-2;PVT_table(7,11)  = 0.232
            PVT_table(5,12) = 5.893d-1; PVT_table(6,12) = 1.08d-2;PVT_table(7,12)  = 0.097
            !#####Dissolved GOR (Oil disolved in gas -Rv-)######
            PVT_table(8,:)  = 0.
        end subroutine populate_with_Texas_Black_Oil

    end subroutine extended_Black_Oil

    subroutine calculate_flooding_height(packed_state, flood_height, deltap, adjust_deltaP)
        !height = P * gravity - bathymetry
        !If deltap then P+deltap is used
        implicit none
        type( state_type ), intent( in ) :: packed_state
        real, dimension(:), optional, intent(inout) :: deltap
        real, dimension(:), intent( out ) :: flood_height
        logical, optional, intent(in) :: adjust_deltaP
        !Local variables
        real, parameter :: gravity_flooding = 9.80665
        type( tensor_field ), pointer :: bathymetry, Pres

        Pres => extract_tensor_field( packed_state, "PackedFEPressure" )
        bathymetry => extract_tensor_field(packed_state,"PackedBathymetry")

        if(present(deltap)) then
            flood_height = (Pres%val(1,1,:) + deltap) * gravity_flooding - bathymetry%val(1,1,:)
            if (present(adjust_deltaP)) then
                !Make sure that the flood height after delta_P is above zero
                where (flood_height < 0)
                    deltap = bathymetry%val(1,1,:)/gravity_flooding - Pres%val(1,1,:)
                    flood_height = 0.
                end where
            end if
        else
            flood_height = Pres%val(1,1,:) * gravity_flooding - bathymetry%val(1,1,:)
        end if

    end subroutine calculate_flooding_height

end module multiphase_EOS
