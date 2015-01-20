
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
    use global_parameters, only: OPTION_PATH_LEN, PYTHON_FUNC_LEN, PI, is_compact_overlapping
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

    type corey_options
       REAL :: S_GC
       real :: S_OR
       real :: s_gr
       real :: c
       real :: s_gi
       real :: cs_gi
       real :: kr1_max
       real :: kr2_max
       real :: kr3_max
       real :: kr1_exp
       real :: kr2_exp
       real :: kr3_exp
       logical :: boost_at_zero_saturation
              
    end type corey_options

    type drag_option_type
       real :: diameter
       real :: coefficient
       real :: Re_s
       integer :: type
       type(function_data) ::python_data
    end type drag_option_type

    type polydispersive_drag_option_type
       real :: lubrication_distance
    end type polydispersive_drag_option_type

    type solid_solid_interaction_phase_option_type
       integer :: type
       real :: coefficient_of_friction
       integer :: radial_distribution_function
    end type solid_solid_interaction_phase_option_type

    
    type solid_solid_interaction_option_type
       type(solid_solid_interaction_phase_option_type),&
            dimension(:), allocatable :: list
    end type solid_solid_interaction_option_type

  contains

    subroutine Calculate_All_Rhos( state, packed_state, ncomp_in, nphase, ndim, cv_nonods, cv_nloc, totele, &
         cv_ndgln, DRhoDPressure )

      implicit none

      type( state_type ), dimension( : ), intent( inout ) :: state
      type( state_type ), intent( inout ) :: packed_state
      integer, intent( in ) :: ncomp_in, nphase, ndim, cv_nonods, cv_nloc, totele
      integer, dimension( : ), intent( in ) :: cv_ndgln

      real, dimension( nphase, cv_nonods ), intent( inout ), optional :: DRhoDPressure

      real, dimension( : ), allocatable :: Rho, dRhodP, Density_Bulk, DensityCp_Bulk, &
                                               Density_Component, Cp, Component_l, c_cv_nod
      character( len = option_path_len ), dimension( : ), allocatable :: eos_option_path
      type( tensor_field ), pointer :: field1, field2, field3, field4
      type( scalar_field ), pointer :: Cp_s, sf
      integer :: icomp, iphase, ncomp, sc, ec, sp, ep, stat, cv_iloc, cv_nod, ele

      DRhoDPressure = 0.

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

            if( ncomp > 1 ) then
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
               DRhoDPressure( iphase, : ) = DRhoDPressure( iphase, : ) + dRhodP * Component_l / Rho
               Density_Component( sc : ec ) = Rho

               Cp_s => extract_scalar_field( state( nphase + icomp ), &
                    'ComponentMassFractionPhase' // int2str( iphase ) // 'HeatCapacity', stat )
               if( stat == 0 ) Cp = Cp_s % val
               DensityCp_Bulk( sp : ep ) = DensityCp_Bulk( sp : ep ) + Rho * Cp * Component_l

            else

               Density_Bulk( sp : ep ) = Rho
               DRhoDPressure( iphase, : ) = dRhodP

               Cp_s => extract_scalar_field( state( iphase ), 'TemperatureHeatCapacity', stat )
               if( stat == 0 ) Cp = Cp_s % val
               DensityCp_Bulk( sp : ep ) = Rho * Cp

            end if

         end do ! iphase
      end do ! icomp

      if( ncomp > 1 ) then
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
         
         field1 % val ( 1, iphase, :) = Density_Bulk( sp : ep )         !* ( 1. - sf%val)     +   1000. *  sf%val
         field2 % val ( 1, iphase, :) = DensityCp_Bulk( sp : ep )       !* ( 1. - sf%val)     +   1000. *  sf%val

        if( ncomp > 1 ) then
            do icomp = 1, ncomp
               sc = ( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1
               ec = ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods
               field3 % val ( icomp, iphase, :) = Density_Component( sc : ec )
            end do ! icomp
         end if
      end do ! iphase


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


    subroutine Calculate_Component_Rho( state, packed_state, &
         ncomp, nphase, cv_nonods )

      implicit none

      type( state_type ), dimension( : ), intent( inout ) :: state
      type( state_type ), intent( inout ) :: packed_state
      integer, intent( in ) :: ncomp, nphase, cv_nonods

      real, dimension( : ), allocatable :: Rho, dRhodP
      type( tensor_field ), pointer :: field
      character( len = option_path_len ) :: eos_option_path
      integer :: icomp, iphase, s, e

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
      character( len = option_path_len ), intent( in ) :: eos_option_path
      real, dimension( : ), intent( inout ) :: rho, drhodp

      type( scalar_field ), pointer :: pressure, temperature, density
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

      pressure => extract_scalar_field( packed_state, 'CVPressure' )
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
         Rho = ( pressure % val + eos_coefs( 1 ) ) * eos_coefs( 2 ) / temperature % val
         perturbation_pressure = max( toler, 1.e-3 * ( abs( pressure % val ) + eos_coefs( 1 ) ) )
         RhoPlus = ( pressure % val + perturbation_pressure + eos_coefs( 1 ) ) *  eos_coefs( 2 ) / &
              temperature % val
         RhoMinus = ( pressure % val - perturbation_pressure + eos_coefs( 1 ) ) *  eos_coefs( 2 ) / &
              temperature % val
         dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure
         deallocate( eos_coefs )

      elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/linear_in_pressure' ) then
!!$ Den = C0 * P +C1
         allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
         call get_option( trim( eos_option_path ) // '/coefficient_A', eos_coefs( 1 ) )
         call get_option( trim( eos_option_path ) // '/coefficient_B', eos_coefs( 2 ) )
         Rho = eos_coefs( 1 ) * pressure % val + eos_coefs( 2 )
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
         Rho = eos_coefs( 1 ) * pressure % val / temperature % val + eos_coefs( 2 )
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
            reference_pressure = pressure % val
            initialised = .true.
         end if
         Rho = eos_coefs( 2 ) * exp( eos_coefs( 1 ) * ( pressure % val - reference_pressure ) )
         perturbation_pressure = max( toler, 1.e-3 * ( abs( pressure % val ) ) )
         RhoPlus = eos_coefs( 2 ) * exp( eos_coefs( 1 ) * ( ( pressure % val + perturbation_pressure ) - &
              reference_pressure ) ) 
         RhoMinus = eos_coefs( 2 ) * exp( eos_coefs( 1 ) * ( ( pressure % val - perturbation_pressure ) - &
              reference_pressure ) ) 
         dRhodP = 0.5 * ( RhoPlus - RhoMinus ) / perturbation_pressure
         deallocate( eos_coefs )

      elseif( trim( eos_option_path ) == trim( option_path_comp ) // '/exponential_in_pressure' ) then 
!!$ Den = C0 * ( P ^ C1 )
         allocate( eos_coefs( 2 ) ) ; eos_coefs = 0.
         call get_option( trim( eos_option_path ) // '/coefficient_A', eos_coefs( 1 ) )
         call get_option( trim( eos_option_path ) // '/coefficient_B', eos_coefs( 2 ) )
         Rho = eos_coefs( 1 ) * pressure % val ** eos_coefs( 2 )
         perturbation_pressure = 1.
         RhoPlus = eos_coefs( 1 ) * ( pressure % val + perturbation_pressure ) ** eos_coefs( 2 )
         RhoMinus = eos_coefs( 1 ) * ( pressure % val - perturbation_pressure ) ** eos_coefs( 2 )
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
         call Density_Polynomial( eos_coefs, pressure % val, temperature_local, &
              Rho )
         perturbation_pressure = max( toler, 1.e-3 * abs( pressure % val ) )
         call Density_Polynomial( eos_coefs, pressure % val + perturbation_pressure, temperature_local, &
              RhoPlus )
         call Density_Polynomial( eos_coefs, pressure % val - perturbation_pressure, temperature_local, &
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
         pressure_back_up = pressure % val
         density_back_up = density % val

         call python_reset()

         ! Calculating dRho / dP
         ! redefine p as p+pert and p-pert and then run python state again to get dRho / d P...
         perturbation_pressure = 1.e-5

         pressure % val = pressure % val + perturbation_pressure
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

         pressure % val = pressure_back_up
         pressure % val = pressure % val - perturbation_pressure
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
         pressure % val = pressure_back_up
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




    subroutine Calculate_AbsorptionTerm( state, packed_state,&
         cv_ndgln, mat_ndgln, &
         opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new, u_absorb )
      ! Calculate absorption for momentum eqns
      use matrix_operations
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      type( state_type ), intent( inout ) :: packed_state
      integer, dimension( : ), intent( in ) :: cv_ndgln, mat_ndgln
      real, dimension( :, :, :, : ), intent( inout ) :: opt_vel_upwind_coefs_new, opt_vel_upwind_grad_new
      real, dimension( :, :, : ), intent( inout ) :: u_absorb
!!$ Local variables:
      type( tensor_field ), pointer :: viscosity_ph
      integer :: nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, x_snloc, cv_snloc, u_snloc, &
           p_snloc, cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, &
           ele, imat, icv, iphase, cv_iloc, idim, jdim, ij,i
      real :: Mobility, pert
      real, dimension( :, :, : ), allocatable :: u_absorb2
!      real, dimension( : ), allocatable :: satura2
        real, dimension( :, : ), allocatable :: satura2
      real, dimension(size(state,1)) :: visc_phases
      !Working pointers
      real, dimension(:,:), pointer :: Satura
      type( tensor_field ), pointer :: perm

    !Get from packed_state
    call get_var_from_packed_state(packed_state,PhaseVolumeFraction = Satura)

    perm=>extract_tensor_field(packed_state,"Permeability")

!!$ Obtaining a few scalars that will be used in the current subroutine tree:
      call Get_Primary_Scalars( state, &         
           nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods )


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
         mobility = 0.
      end if


      allocate( u_absorb2( mat_nonods, nphase * ndim, nphase * ndim ), satura2( size(SATURA,1), size(SATURA,2) ) )
      u_absorb2 = 0. ; satura2 = 0.

      CALL calculate_absorption2( MAT_NONODS, CV_NONODS, NPHASE, NDIM, SATURA, TOTELE, CV_NLOC, MAT_NLOC, &
           CV_NDGLN, MAT_NDGLN, U_ABSORB, PERM%val, MOBILITY, visc_phases)

      PERT = 0.0001

      SATURA2( 1,: ) = SATURA( 1,: ) + PERT
      IF ( NPHASE > 1 ) SATURA2( 2,: ) = SATURA( 2,: ) - PERT

      CALL calculate_absorption2( MAT_NONODS, CV_NONODS, NPHASE, NDIM, SATURA2, TOTELE, CV_NLOC, MAT_NLOC, &
           CV_NDGLN, MAT_NDGLN, &
           U_ABSORB2, PERM%val, MOBILITY, visc_phases)

      DO ELE = 1, TOTELE
         DO CV_ILOC = 1, CV_NLOC
            IMAT = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC +CV_ILOC )
            ICV = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )
            DO IPHASE = 1, NPHASE
              DO JDIM = 1, NDIM
                DO IDIM = 1, NDIM
                     IJ = ( IPHASE - 1 ) * MAT_NONODS * NDIM * NDIM + ( IMAT - 1 ) * NDIM * NDIM + &
                          ( IDIM - 1 ) * NDIM + JDIM
                     opt_vel_upwind_coefs_new(IDIM, JDIM, IPHASE, IMAT) = &
                        U_ABSORB( IMAT, IDIM + ( IPHASE - 1 ) * NDIM, JDIM + ( IPHASE - 1 ) * NDIM )

!                     ! This is the gradient...
                     opt_vel_upwind_grad_new(IDIM, JDIM, IPHASE, IMAT) = &
                        (U_ABSORB2( IMAT, IDIM + ( IPHASE - 1 ) * NDIM, JDIM + ( IPHASE - 1 ) * NDIM ) -&
                        U_ABSORB( IMAT, IDIM + ( IPHASE - 1 ) * NDIM, JDIM + ( IPHASE - 1 ) * NDIM ))&
                         / ( SATURA2(IPHASE, ICV ) - SATURA(IPHASE, ICV))

                  END DO
               END DO
            END DO
         END DO
      END DO

      deallocate( u_absorb2, satura2 )
      ewrite(3,*) 'Leaving calculate_absorption'

      RETURN

    END SUBROUTINE Calculate_AbsorptionTerm




    SUBROUTINE calculate_absorption2( MAT_NONODS, CV_NONODS, NPHASE, NDIM, SATURA, TOTELE, CV_NLOC, MAT_NLOC, &
         CV_NDGLN, MAT_NDGLN, &
         U_ABSORB, PERM2, MOBILITY, visc_phases)
      ! Calculate absorption for momentum eqns
      use matrix_operations
      !    use cv_advection
      implicit none
      INTEGER, intent( in ) :: MAT_NONODS, CV_NONODS, NPHASE, NDIM, TOTELE, CV_NLOC,MAT_NLOC
      REAL, DIMENSION( :, : ), intent( in ) :: SATURA
      INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
      REAL, DIMENSION( :, :, : ), intent( inout ) :: U_ABSORB
      REAL, DIMENSION( :, :, : ), intent( in ) :: PERM2
      REAL, intent( in ) :: MOBILITY
      real, intent(in), dimension(:) :: visc_phases
      ! Local variable
      REAL, PARAMETER :: TOLER = 1.E-10
      INTEGER :: ELE, CV_ILOC, CV_NOD, CV_PHA_NOD, MAT_NOD, JPHA_JDIM, &
           IPHA_IDIM, IDIM, JDIM, IPHASE, jphase
      !    integer :: ii
      REAL :: SATURATION
      !    real :: abs_sum
      REAL, DIMENSION( :, :, :), allocatable :: INV_PERM, PERM
      type(corey_options) :: options
      logical :: is_corey


      ewrite(3,*) 'In calculate_absorption2'
      ALLOCATE( INV_PERM(  NDIM, NDIM, TOTELE ))
      ALLOCATE( PERM( NDIM, NDIM , TOTELE))

      perm=perm2
      do ele = 1, totele
         inv_perm( :, :, ele)=inverse(perm( :, :, ele))
      end do
      U_ABSORB = 0.0

      Loop_NPHASE: DO IPHASE = 1, NPHASE

          is_Corey=.false.
          if (have_option("/material_phase["// int2str(iphase-1) //&
          "]/multiphase_properties/relperm_type/Corey")) then
              is_Corey=.true.
              call get_corey_options(options, nphase)
          end if
          Loop_ELE: DO ELE = 1, TOTELE

              Loop_CVNLOC: DO CV_ILOC = 1, CV_NLOC

                  MAT_NOD = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + CV_ILOC)
                  CV_NOD = CV_NDGLN(( ELE - 1) * CV_NLOC + CV_ILOC )

                  Loop_DimensionsI: DO IDIM = 1, NDIM

                      Loop_DimensionsJ: DO JDIM = 1, NDIM

                          CV_PHA_NOD = CV_NOD + ( IPHASE - 1 ) * CV_NONODS
                          SATURATION = SATURA( IPHASE,CV_NOD )

                          IPHA_IDIM = ( IPHASE - 1 ) * NDIM + IDIM
                          JPHA_JDIM = ( IPHASE - 1 ) * NDIM + JDIM
                          if (is_corey) then
                              select case (NPHASE)
                                  case (1)!No relperm needed, we calculate directly the result
                                      U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = INV_PERM( IDIM, JDIM, ELE ) *&
                                      visc_phases(1) * min(1.0,max(1d-5,SATURA(1,CV_NOD)))
                                  case (2)
                                      CALL relperm_corey_epsilon( U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ), visc_phases(1), visc_phases(2), &
                                      INV_PERM( IDIM, JDIM, ELE ), min(1.0,max(0.0,SATURA(1,CV_NOD))), IPHASE,&
                                      options)!Second phase is considered inside the subroutine
                                  case (3)!For three phases we use the Stone model. !With predefined order: Water, oil, gas
                                      call relperm_stone(U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ), iphase,&
                                      SATURA(:, CV_NOD), options, visc_phases(iphase), INV_PERM( IDIM, JDIM, ELE ))
                                  case default
                                      FLAbort("No relative permeability function implemented for more than 3 phases")
                              end select

                          else
                              U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = 0.0

                          endif
                      END DO Loop_DimensionsJ

                  END DO Loop_DimensionsI

              END DO Loop_CVNLOC

          END DO Loop_ELE

      END DO Loop_NPHASE

      DEALLOCATE( PERM, INV_PERM )

      ewrite(3,*) 'Leaving calculate_absorption2'

      RETURN

    END SUBROUTINE calculate_absorption2
    

    subroutine initialize_drag_matrix(absorption,state,nphase)
      type(block_csr_matrix), intent(inout) :: absorption
      type(state_type) :: state
      integer :: nphase

      type(vector_field), pointer :: velocity
      type(scalar_field), pointer :: volume_fraction
      type(csr_sparsity), pointer :: sparsity
      type(mesh_type), pointer :: mesh
      
      velocity=>extract_vector_field(state,"Velocity")
      volume_fraction=>extract_scalar_field(state,"PhaseVolumeFraction")
      mesh=>extract_mesh(state,"PressureMesh_Discontinuous")


      sparsity => get_csr_sparsity_firstorder(state, mesh,mesh)

      call allocate(absorption,sparsity,(/nphase*velocity%dim,nphase*velocity%dim/)&
           ,name='DragAbsorption')

      call zero(absorption)


    end subroutine initialize_drag_matrix

    subroutine add_drag_to_old_style_matrix_and_cleanup(drag_abs_matrix,U_ABSORB,state,&
              MAT_NONODS, CV_NONODS, NPHASE, NDIM, TOTELE, CV_NLOC, MAT_NLOC, &
           CV_NDGLN, MAT_NDGLN)
       type(block_csr_matrix), intent(inout) :: drag_abs_matrix
       real, dimension( :, :, : ), intent( inout ) :: u_absorb
       type(state_type) :: state

       INTEGER, intent( in ) :: MAT_NONODS, CV_NONODS, NPHASE, NDIM, TOTELE, CV_NLOC,MAT_NLOC
       INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN

       type(mesh_type), pointer :: mesh
       integer :: ELE,IPHASE,IDIM,JPHASE,JDIM,I,J,MAT_NOD,ICV_LOC
       integer, dimension(:), pointer :: nodes

       mesh=>extract_mesh(state,"PressureMesh_Discontinuous")

       DO ELE=1,TOTELE
          nodes=>ele_nodes(mesh,ele)
          do ICV_LOC=1,CV_NLOC
             MAT_NOD = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC+ICV_LOC)
             DO IPHASE=1,NPHASE
                DO IDIM=1,NDIM
                   I=(IPHASE-1)*NDIM+IDIM
                   DO JPHASE=1,NPHASE
                      DO JDIM=1,NDIM
                         J=(JPHASE-1)*NDIM+JDIM
                         U_ABSORB(MAT_NOD,I,J)=U_ABSORB(MAT_NOD,I,J)+&
                              val(drag_abs_matrix,I,J,nodes(ICV_LOC),nodes(ICV_LOC))
                      end DO
                   end Do
                end DO
             end DO
          end DO
       end do

       call deallocate(drag_abs_matrix)

    end subroutine add_drag_to_old_style_matrix_and_cleanup
      

    function radial_distribution_function_Lebowitz(&
         d_i,d_j,phi_c,phi_i,phi_j) result(g_ij)
      real, intent(in) :: d_i,d_j
      real, dimension(:), intent(in) :: phi_c,phi_i,phi_j
      real, dimension(size(phi_i))  :: g_ij

      g_ij=1.0/phi_c*(1.0+3.0*d_i*d_j/(phi_c*(d_i+d_j))&
           *(phi_i/d_i+phi_j/d_j))

    end function radial_distribution_function_Lebowitz

    subroutine add_solid_solid_interaction_term(state,absorption,&
         continuous_phase,&
         iphase,jphase,&
         iphase_drag_options,jphase_drag_options,&
         iphase_solid_solid_interaction_options,&
         jphase_solid_solid_interaction_options)
      
      type(state_type), dimension(:), intent(inout) :: state
      type(block_csr_matrix), intent(inout) :: absorption
      integer, intent(in) :: continuous_phase, iphase,jphase
      type(drag_option_type), intent(in) :: iphase_drag_options, jphase_drag_options
      type(solid_solid_interaction_option_type), intent(in) ::&
           iphase_solid_solid_interaction_options,&
           jphase_solid_solid_interaction_options


      type(scalar_field) :: ibeta,jbeta,cval,cvf_dg,ivf_dg,jvf_dg
      type(scalar_field), pointer :: cvf,ivf,jvf,rho_i,rho_j
      type(mesh_type), pointer :: mesh
      type(vector_field),pointer :: ivel,jvel
      type(vector_field) :: deltaV
      integer :: block1,block2,ele,idim,ndim,iopt,jopt,stat
      integer, dimension(:), pointer :: dg_nodes, nodes
      real :: d_i,d_j
      real, dimension(:), allocatable :: r_i,r_j,g_ij


      !! extract data from state 

      mesh=>extract_mesh(state(1),"PressureMesh_Discontinuous")
      CVF=>extract_scalar_field(state(continuous_phase),&
           "PhaseVolumeFraction",stat)
      IVF=>extract_scalar_field(state(iphase),"PhaseVolumeFraction",stat)
      JVF=>extract_scalar_field(state(jphase),"PhaseVolumeFraction",stat)
      rho_i=>extract_scalar_field(state(jphase),"Density",stat)
      rho_j=>extract_scalar_field(state(jphase),"Density",stat)
      ivel=>extract_vector_field(state(iphase),"Velocity",stat)
      jvel=>extract_vector_field(state(iphase),"Velocity",stat)


      !! Useful convenience variables
      
      d_i=iphase_drag_options%diameter
      d_j=jphase_drag_options%diameter

      call allocate(cval,mesh,"working_array")
      call allocate(CVF_DG,mesh,"IVF_DG")
      call allocate(IVF_DG,mesh,"IVF_DG")
      call allocate(JVF_DG,mesh,"JVF_DG")
      call allocate(deltaV,mesh_dim(mesh),mesh,"deltaV")
      call remap_field(CVF,CVF_DG)
      call remap_field(IVF,IVF_DG)
      call remap_field(JVF,JVF_DG)
      call remap_field(ivel,deltaV)
      call addto(deltaV,jvel, scale=-1.0)

      call zero(cval)

      allocate(r_i(ele_loc(rho_i,1)),&
           r_j(ele_loc(rho_j,1)),&
           g_ij(ele_loc(mesh,1)))
            
      do iopt=1, size(iphase_solid_solid_interaction_options%list)
         do jopt=1, size(iphase_solid_solid_interaction_options%list)
            if (iphase_solid_solid_interaction_options%list(iopt)%type/=&
                jphase_solid_solid_interaction_options%list(jopt)%type) &
                 cycle

            do ele=1,ele_count(ibeta)
               nodes=>ele_nodes(rho_i,ele)
               dg_nodes=>ele_nodes(mesh,ele)
               r_i=ele_val(rho_i,ele)
               r_j=ele_val(rho_i,ele)
               g_ij=radial_distribution_function_Lebowitz(d_i,d_j,&
                    cvf_dg%val(dg_nodes),&
                    ivf_dg%val(dg_nodes),&
                    jvf_dg%val(dg_nodes))
               cval%val(dg_nodes)= cval%val(dg_nodes)&
                    +3.0*PI/4.0*(1.0+iphase_drag_options%coefficient)&
                    *(1.0+iphase_solid_solid_interaction_options%list(iopt)%coefficient_of_friction*PI/4.0)&
                    *r_i*r_j*g_ij*(d_i+d_j)**2/&
                    (r_i*d_i**3+r_j*d_j**3)&
                 *sqrt(sum(deltaV%val*deltaV%val(:,dg_nodes),dim=1))
            end do
         end do
      end do

      ndim=mesh_dim(mesh)

      do IDIM=1,ndim
         block1=(iphase-1)*ndim+IDIM
         block2=(jphase-1)*ndim+IDIM
         
         do ele=1,ele_count(mesh)

            dg_nodes=>ele_nodes(mesh,ele)

            call addto_diag(absorption,block1,block1,dg_nodes, jvf_dg%val(dg_nodes)*cval%val(dg_nodes))
            call addto_diag(absorption,block1,block2,dg_nodes, -jvf_dg%val(dg_nodes)*cval%val(dg_nodes))
            call addto_diag(absorption,block2,block1,dg_nodes, -ivf_dg%val(dg_nodes)*cval%val(dg_nodes))
            call addto_diag(absorption,block2,block2,dg_nodes, ivf_dg%val(dg_nodes)*cval%val(dg_nodes))

         end do

      end do

      call deallocate(cval)
      call deallocate(cvf_dg)
      call deallocate(ivf_dg)
      call deallocate(jvf_dg)
      call deallocate(deltaV)
      
      deallocate(r_i,r_j,g_ij)
      


    end subroutine add_solid_solid_interaction_term

    subroutine add_polydispersive_drag_term(state,absorption,iphase,jphase,&
                       iphase_drag_options,jphase_drag_options,&
                       iphase_polydispersive_drag_options,&
                       jphase_polydispersive_drag_options)
      
      type(state_type), dimension(:), intent(inout) :: state
      type(block_csr_matrix), intent(inout) :: absorption
      integer, intent(in) :: iphase,jphase
      type(drag_option_type), intent(in) :: iphase_drag_options, jphase_drag_options
      type(polydispersive_drag_option_type), intent(in) ::&
           iphase_polydispersive_drag_options,jphase_polydispersive_drag_options

      type(scalar_field) :: ibeta,jbeta,cval,ivf_dg,jvf_dg
      type(scalar_field), pointer :: ivf,jvf
      type(mesh_type), pointer :: mesh

      real :: alpha
      integer :: block1,block2, ele,idim,ndim,stat
      integer, dimension(:), pointer :: dg_nodes

      mesh=>extract_mesh(state(1),"PressureMesh_Discontinuous")
      IVF=>extract_scalar_field(state(iphase),"PhaseVolumeFraction",stat)
      JVF=>extract_scalar_field(state(jphase),"PhaseVolumeFraction",stat)

      call allocate(ibeta,mesh,"Drag_coefficient_iphase")
      call allocate(jbeta,mesh,"Drag_coefficient_jphase")
      call allocate(cval,mesh,"working_array")
      call allocate(IVF_DG,mesh,"IVF_DG")
      call allocate(JVF_DG,mesh,"JVF_DG")
      call remap_field(IVF,IVF_DG)
      call remap_field(JVF,JVF_DG)

      alpha=1.313*log10(min(iphase_drag_options%diameter, jphase_drag_options%diameter)/&
           max(iphase_polydispersive_drag_options%lubrication_distance,&
           jphase_polydispersive_drag_options%lubrication_distance))-1.249


      cval%val=(ivf_dg%val/ibeta%val+jvf_dg%val/jbeta%val)

      call calculate_drag_as_scalar(state, ibeta,iphase,jphase,&
           iphase_drag_options,ivf_dg,jvf_dg)
      call calculate_drag_as_scalar(state, jbeta,jphase,iphase,&
           jphase_drag_options,ivf_dg,jvf_dg)


      ndim=mesh_dim(mesh)
      do IDIM=1,ndim
         block1=(iphase-1)*ndim+IDIM
         block2=(jphase-1)*ndim+IDIM

         do ele=1,ele_count(ibeta)

            dg_nodes=>ele_nodes(mesh,ele)
            call addto_diag(absorption,block1,block1,dg_nodes, 2*alpha*jvf_dg%val(dg_nodes)/cval%val(dg_nodes))
            call addto_diag(absorption,block1,block2,dg_nodes, -2*alpha*jvf_dg%val(dg_nodes)/cval%val(dg_nodes))
            call addto_diag(absorption,block2,block1,dg_nodes, -2*alpha*ivf_dg%val(dg_nodes)/cval%val(dg_nodes))
            call addto_diag(absorption,block2,block2,dg_nodes, 2*alpha*ivf_dg%val(dg_nodes)/cval%val(dg_nodes))

         end do

      end do

      call deallocate(ibeta)
      call deallocate(jbeta)
      call deallocate(cval)
      call deallocate(ivf)
      call deallocate(jvf)

    end subroutine add_polydispersive_drag_term


    subroutine calculate_drag_as_scalar(state, beta,iphase,jphase,&
         drag_options,ivf_dg,jvf_dg)
      type(state_type), dimension(:), intent(inout) :: state
      type(scalar_field) :: beta
      integer :: iphase, jphase
      type(drag_option_type) :: drag_options
      type(scalar_field) :: ivf_dg,jvf_dg
      
      type(vector_field), pointer :: ivel, jvel
      type(vector_field) :: deltaV
      
      call zero(beta)
      call allocate(deltaV,mesh_dim(beta%mesh),beta%mesh,"Veldifference")

      select case(drag_options%type)
      case(0)
         call calculate_function_at_nodes(drag_options%python_data,beta)
      case(1)
         beta%val=3.0d0*drag_options%coefficient/(4.0d0*drag_options%diameter)*ivf_dg%val
      case(-1)
         beta%val=3.0d0*drag_options%coefficient/(4.0d0*drag_options%diameter)*ivf_dg%val*jvf_dg%val
      case(2)
         call remap_field(ivel,deltaV)
         call addto(deltaV,jvel, scale=-1.0)
         beta%val=3.0d0*drag_options%coefficient/(4.0d0*drag_options%diameter)*ivf_dg%val&
              *sqrt(sum(deltaV%val(:,:)*deltaV%val(:,:),dim=1))


      case(3)

         !beta%val=merge(0.44,24.0/Re_s%val&
         !     *(1.0d0+0.15*(Re_S%val)**0.687)&
         !     ,Re_s%val>1000.0)

      end select

      call deallocate(deltaV)
      

    end subroutine calculate_drag_as_scalar
      

    subroutine add_drag_term(state,absorption,dispersed_phase,continuous_phase,drag_options)
      type(state_type), dimension(:), intent(inout) :: state
      type(block_csr_matrix), intent(inout) :: absorption
      type(vector_field), pointer :: continuous_velocity, dispersed_velocity 
      integer, intent(in) :: dispersed_phase, continuous_phase
      type(drag_option_type) :: drag_options
      type(scalar_field), pointer :: continuous_volume_fraction, dispersed_volume_fraction,rho_g
      type(tensor_field),pointer :: mu_g
      type(scalar_field) :: ratio, coeff, CVF_DG,DVF_DG,Re_s,psi
      type(vector_field) :: deltaV
      type(mesh_type), pointer :: mesh, cv_mesh

      real :: cval
      integer :: node,ndim, dim,block1,block2, IDIM, stat,ele, i,nloc
      integer, dimension(:), pointer :: nodes
      integer, dimension(:), pointer :: dg_nodes
      type(element_type) :: p_cvshape_full

      continuous_velocity=>extract_vector_field(state(continuous_phase),"Velocity",stat)

      
!!      mesh=>extract_mesh(state(1),"PressureMesh_Continuous")
    
      mesh=>extract_mesh(state(1),"PressureMesh_Discontinuous")
      
      if (stat/=0) then
         FLAbort("Attempting to add drag term, when no velocity in continuous phase")
      end if
      dispersed_velocity =>extract_vector_field(state(dispersed_phase ),"Velocity",stat)
      if (stat/=0) then
         FLAbort("Attempting to add drag term, when no velocity in dispersed phase")
      end if
      ndim=continuous_velocity%dim
      continuous_volume_fraction=>extract_scalar_field(state(continuous_phase),"PhaseVolumeFraction",stat)
      dispersed_volume_fraction=>extract_scalar_field(state(dispersed_phase),"PhaseVolumeFraction",stat)
      mu_g=>extract_tensor_field(state(continuous_phase),"Viscosity",stat)
      rho_g=>extract_scalar_field(state(continuous_phase),"Density",stat)
      nloc=ele_loc(continuous_volume_fraction,1)

      call allocate(CVF_DG,mesh,"CVF_DG")
      call allocate(DVF_DG,mesh,"CVF_DG")
      call remap_field(continuous_volume_fraction,CVF_DG)
      call remap_field(dispersed_volume_fraction,DVF_DG)
      call allocate(ratio,mesh,"RatioOfVolumeFractions")
      call zero(ratio)
      call set(ratio,CVF_DG)
      call invert(ratio)
      call scale(ratio,DVF_DG)
      call bound(ratio,0.0d0,1.0e12)

      select case(drag_options%type)
      case(0)
         call allocate(coeff,mesh,"Drag")
         call calculate_function_at_nodes(drag_options%python_data,coeff)
         do IDIM=1,ndim
            block1=(continuous_phase-1)*ndim+IDIM
            block2=(dispersed_phase-1)*ndim+IDIM
            do ele=1,ele_count(continuous_volume_fraction)
               nodes=>ele_nodes(continuous_volume_fraction,ele)
               dg_nodes=>ele_nodes(mesh,ele)
               call addto_diag(absorption,block1,block1,dg_nodes, coeff%val(nodes)*ratio%val(dg_nodes))
               call addto_diag(absorption,block1,block2,dg_nodes,-coeff%val(nodes)*ratio%val(dg_nodes))
               call addto_diag(absorption,block2,block1,dg_nodes,-coeff%val(dg_nodes))
               call addto_diag(absorption,block2,block2,dg_nodes, coeff%val(dg_nodes))
            end do
         end do
         call deallocate(coeff)
      case(1)
         cval=3.0d0*drag_options%coefficient/(4.0d0*drag_options%diameter)
         do IDIM=1,ndim
            block1=(continuous_phase-1)*ndim+IDIM
            block2=(dispersed_phase-1)*ndim+IDIM
            do ele=1,ele_count(continuous_volume_fraction)
               nodes=>ele_nodes(continuous_volume_fraction,ele)
               dg_nodes=>ele_nodes(mesh,ele)
               call addto_diag(absorption,block1,block1,dg_nodes, cval*ratio%val(dg_nodes))
               call addto_diag(absorption,block1,block2,dg_nodes,-cval*ratio%val(dg_nodes))
               call addto_diag(absorption,block2,block1,dg_nodes,spread(-cval,1,nloc))
               call addto_diag(absorption,block2,block2,dg_nodes,spread( cval,1,nloc))
            end do
         end do
      case(-1)
         cval=3.0d0*drag_options%coefficient/(4.0d0*drag_options%diameter)
         do IDIM=1,ndim
            block1=(continuous_phase-1)*ndim+IDIM
            block2=(dispersed_phase-1)*ndim+IDIM
            do ele=1,ele_count(continuous_volume_fraction)
               nodes=>ele_nodes(continuous_volume_fraction,ele)
               dg_nodes=>ele_nodes(mesh,ele)
               call addto_diag(absorption,block1,block1,dg_nodes, cval*dispersed_volume_fraction%val(nodes))
               call addto_diag(absorption,block1,block2,dg_nodes,-cval*dispersed_volume_fraction%val(nodes))
               call addto_diag(absorption,block2,block1,dg_nodes,-cval*continuous_volume_fraction%val(nodes))
               call addto_diag(absorption,block2,block2,dg_nodes, cval*continuous_volume_fraction%val(nodes))
            end do
         end do
      case(2)
         call allocate(deltaV,ndim,mesh,"DeltaV")

         call remap_field(continuous_velocity,deltaV)
         call addto(deltaV,dispersed_velocity, scale=-1.0)

         cval=3.0d0*drag_options%coefficient/(4.0d0*drag_options%diameter)

         do IDIM=1,ndim
            block1=(continuous_phase-1)*ndim+IDIM
            block2=(dispersed_phase-1)*ndim+IDIM
            do ele=1,ele_count(continuous_volume_fraction)
               nodes=>ele_nodes(continuous_volume_fraction,ele)
               dg_nodes=>ele_nodes(mesh,ele)

               call addto_diag(absorption,block1,block1,dg_nodes, cval*ratio%val(nodes)&
                    *sqrt(sum(deltaV%val(:,dg_nodes)*deltaV%val(:,dg_nodes),dim=1)))
               call addto_diag(absorption,block1,block2,dg_nodes,-cval*ratio%val(nodes)&
                    *sqrt(sum(deltaV%val(:,dg_nodes)*deltaV%val(:,dg_nodes),dim=1)))
               call addto_diag(absorption,block2,block1,dg_nodes,spread(-cval,1,nloc)&
                    *sqrt(sum(deltaV%val(:,dg_nodes)*deltaV%val(:,dg_nodes),dim=1)))
               call addto_diag(absorption,block2,block2,dg_nodes,spread( cval,1,nloc)&
                    *sqrt(sum(deltaV%val(:,dg_nodes)*deltaV%val(:,dg_nodes),dim=1)))
            end do
         end do

         call deallocate(deltaV)

         case(3)

         call allocate(deltaV,ndim,mesh,"DeltaV")
         call allocate(coeff,mesh,"Drag")
         call allocate(Re_s,mesh,"Drag")
         call remap_field(continuous_velocity,deltaV)
         call addto(deltaV,dispersed_velocity, scale=-1.0)


         do ele=1,ele_count(continuous_volume_fraction)
            nodes=>ele_nodes(mu_g,ele)
            dg_nodes=>ele_nodes(mesh,ele)
            Re_s%val(dg_nodes)=rho_g%val(nodes)/mu_g%val(1,1,1)&
                 *sqrt(sum(ele_val(deltaV,ele)**2)/nloc)&
                 *drag_options%diameter
         end do
         call scale(Re_s,CVF_DG)


         coeff%val=merge(0.44,24.0/Re_s%val&
              *(1.0d0+0.15*(Re_S%val)**0.687)&
              ,Re_s%val>1000.0)

         call bound(coeff,0.0d0,1.0e12)


  !       p_cvshape_full%n(1,:)=(/ 14.0/24.0,5.0/24.0,5.0/24.0 /)
  !       p_cvshape_full%n(2,:)=(/ 5.0/24.0,14.0/24.0,5.0/24.0 /)
  !       p_cvshape_full%n(3,:)=(/ 5.0/24.0,5.0/24.0,14.0/24.0 /)

!         p_cvshape_full%n=1.0/nloc


         do ele=1,ele_count(continuous_volume_fraction)
            nodes=>ele_nodes(mu_g,ele)
            dg_nodes=>ele_nodes(CVF_DG,ele)
            coeff%val(dg_nodes)=3.0/4.0*coeff%val(dg_nodes)&
                 *sqrt(sum(ele_val(deltaV,ele)**2)/nloc)&
              /drag_options%diameter&
              *CVF_DG%val(dg_nodes)**(-2.65)
         end do

!         print*, maxval(coeff%val), minval(coeff%val), minval(Re_s)
         
         do IDIM=1,ndim
            block1=(continuous_phase-1)*ndim+IDIM
            block2=(dispersed_phase-1)*ndim+IDIM
            do ele=1,ele_count(continuous_volume_fraction)
               nodes=>ele_nodes(continuous_volume_fraction,ele)
               dg_nodes=>ele_nodes(mesh,ele)

               call addto_diag(absorption,block1,block1,dg_nodes, rho_g%val(nodes)*coeff%val(dg_nodes)*DVF_DG%val(dg_nodes))
               call addto_diag(absorption,block1,block2,dg_nodes,-rho_g%val(nodes)*coeff%val(dg_nodes)*DVF_DG%val(dg_nodes))
               call addto_diag(absorption,block2,block1,dg_nodes,-rho_g%val(nodes)*coeff%val(dg_nodes)*CVF_DG%val(dg_nodes))
               call addto_diag(absorption,block2,block2,dg_nodes, rho_g%val(nodes)*coeff%val(dg_nodes)*CVF_DG%val(dg_nodes))
            end do
         end do

         call deallocate(deltaV)
         call deallocate(coeff)
         call deallocate(Re_s)
         call deallocate(CVF_DG)
         call deallocate(DVF_DG)
      case(4)
         call allocate(deltaV,ndim,mesh,"DeltaV")
         call allocate(coeff,mesh,"Drag")
         call allocate(psi,mesh,"Psi")
         call allocate(Re_s,mesh,"SolidReynoldsNumber")
         call remap_field(continuous_velocity,deltaV)
         call addto(deltaV,dispersed_velocity, scale=-1.0)


         do ele=1,ele_count(continuous_volume_fraction)
            nodes=>ele_nodes(mu_g,ele)
            dg_nodes=>ele_nodes(mesh,ele)
            Re_s%val(dg_nodes)=rho_g%val(nodes)/mu_g%val(1,1,1)&
                 *sqrt(sum(ele_val(deltaV,ele)**2)/nloc)&
                 *drag_options%diameter
         end do
         call scale(Re_s,CVF_DG)

         psi%val=atan(150*1.75*(CVF_DG%val-0.8))/pi+0.5

         coeff%val=merge(0.44,24.0/Re_s%val&
              *(1.0d0+0.15*(Re_S%val)**0.687)&
              ,Re_s%val>1000.0)
         call bound(coeff,0.0d0,1.0e12)


  !       p_cvshape_full%n(1,:)=(/ 14.0/24.0,5.0/24.0,5.0/24.0 /)
  !       p_cvshape_full%n(2,:)=(/ 5.0/24.0,14.0/24.0,5.0/24.0 /)
  !       p_cvshape_full%n(3,:)=(/ 5.0/24.0,5.0/24.0,14.0/24.0 /)

  !       p_cvshape_full%n=1.0/nloc

         do ele=1,ele_count(continuous_volume_fraction)
            nodes=>ele_nodes(mu_g,ele)
            dg_nodes=>ele_nodes(mesh,ele)
            coeff%val(dg_nodes)=(1.0-psi%val(dg_nodes))&
                 *(150.0*DVF_DG%val(dg_nodes)**2*mu_g%val(1,1,1)&
                 /(CVF_DG%val(dg_nodes)*drag_options%diameter)**2&
                 +1.75*rho_g%val(nodes)*DVF_DG%val(dg_nodes)&
                 *sqrt(sum(ele_val(deltaV,ele)**2)/nloc)&
                 /(CVF_DG%val(dg_nodes)*drag_options%diameter))&
                 +psi%val(dg_nodes)*3.0/4.0*coeff%val(dg_nodes)&
                 *sqrt(sum(ele_val(deltaV,ele)**2)/nloc)&
              /drag_options%diameter&
              *CVF_DG%val(dg_nodes)**(-2.65)
         end do

         print*, maxval(coeff%val), minval(coeff%val), minval(Re_s)
         
         do IDIM=1,ndim
            block1=(continuous_phase-1)*ndim+IDIM
            block2=(dispersed_phase-1)*ndim+IDIM
            do ele=1,ele_count(continuous_volume_fraction)
               nodes=>ele_nodes(continuous_volume_fraction,ele)
               dg_nodes=>ele_nodes(mesh,ele)

               call addto_diag(absorption,block1,block1,dg_nodes, rho_g%val(nodes)*coeff%val(dg_nodes)*DVF_DG%val(dg_nodes))
               call addto_diag(absorption,block1,block2,dg_nodes,-rho_g%val(nodes)*coeff%val(dg_nodes)*DVF_DG%val(dg_nodes))
               call addto_diag(absorption,block2,block1,dg_nodes,-rho_g%val(nodes)*coeff%val(dg_nodes)*CVF_DG%val(dg_nodes))
               call addto_diag(absorption,block2,block2,dg_nodes, rho_g%val(nodes)*coeff%val(dg_nodes)*CVF_DG%val(dg_nodes))
            end do
         end do

         call deallocate(deltaV)
         call deallocate(coeff)
         call deallocate(Re_s)
         call deallocate(CVF_DG)
         call deallocate(DVF_DG)
      end select

      call deallocate(ratio)

    end subroutine add_drag_term
    

    subroutine get_solid_solid_interaction_options(phase,state,options)
      integer :: phase
      type(state_type), dimension(:) :: state
      type(solid_solid_interaction_option_type) :: options

      character(len= OPTION_PATH_LEN) :: type, radial_distribution_fn

      integer :: nclosures, closure

      nclosures=option_count("/material_phase["//int2str(phase-1)//&
           "]/multiphase_properties//solid_solid_interaction/closure")

      allocate(options%list(nclosures))

      do closure=1,nclosures

         call get_option("/material_phase["//int2str(phase-1)//&
           "]/multiphase_properties//solid_solid_interaction/closure["&
         //int2str(closure-1)//"]/name",&
         type)
         select case(trim(type))
         case("Syamlal")
            options%list(closure)%type=0
            call get_option("/material_phase["//int2str(phase-1)//&
             "]/multiphase_properties//solid_solid_interaction/closure["&
             //int2str(closure-1)//"]/coefficient_of_friction",&
             options%list(closure)%coefficient_of_friction,default=0.15)
            call get_option("/material_phase["//int2str(phase-1)//&
             "]/multiphase_properties//solid_solid_interaction/closure["&
             //int2str(closure-1)//"]/radial_diatribution_function/name",&
                  radial_distribution_fn)
             select case(trim(radial_distribution_fn))
             case("Lebowitz")
                options%list(closure)%radial_distribution_function=0
             end select
          end select
       end do

    end subroutine get_solid_solid_interaction_options
  
    subroutine get_drag_options(phase,state,options)
      integer :: phase
      type(state_type), dimension(:) :: state
      type(drag_option_type) :: options
      character(len= OPTION_PATH_LEN) :: type

      call get_option("/material_phase["//int2str(phase-1)//"]/multiphase_properties/drag/name",&
           type)

      select case(trim(type))
      case("python_function")
         options%type=0
         call initialize_function_data(state,&
              "/material_phase["//int2str(phase-1)//"]/multiphase_properties/drag",&
              options%python_data)
      case("Linear")
         options%type=1
         call get_option("/material_phase["//int2str(phase-1)//"]/multiphase_properties/drag/diameter", &
              options%diameter, default=0.001)
         call get_option("/material_phase["//int2str(phase-1)//"]/multiphase_properties/drag/coefficient", &
              options%coefficient, default=0.001)
         case("Bilinear")
         options%type=-1
         call get_option("/material_phase["//int2str(phase-1)//"]/multiphase_properties/drag/diameter", &
              options%diameter, default=0.001)
         call get_option("/material_phase["//int2str(phase-1)//"]/multiphase_properties/drag/coefficient", &
              options%coefficient, default=0.001)
      case("Quadratic")
         options%type=2
         call get_option("/material_phase["//int2str(phase-1)//"]/multiphase_properties/drag/diameter", &
              options%diameter, default=0.001)
         call get_option("/material_phase["//int2str(phase-1)//"]/multiphase_properties/drag/coefficient", &
              options%coefficient, default=0.001)
      case("Wen&Yu")
         options%type=3
         call get_option("/material_phase["//int2str(phase-1)//"]/multiphase_properties/drag/diameter", &
              options%diameter, default=0.001)
         case("SmoothedGidaspow")
         options%type=4
         call get_option("/material_phase["//int2str(phase-1)//"]/multiphase_properties/drag/diameter", &
              options%diameter, default=0.001)
      case default
         FLAbort("Unknown drag type in phase "//state(phase)%name)
      end select


     end subroutine get_drag_options

    subroutine clean_drag_options(options,clean_option)
      type(drag_option_type), dimension(:) :: options
      logical, dimension(:) :: clean_option

      integer :: i

      do i=1,size(options)
         if( clean_option(i)) then
            select case(options(i)%type)
            case(0)
               call finalize_function_data(options(i)%python_data)
            end select
         end if
      end do
      
    end subroutine clean_drag_options

    subroutine get_corey_options(options, nphase)
        implicit none
        !Retrieves the values for the Corey model
        !Including the inmobile fractions
        integer, intent(in) :: nphase
      type(corey_options), intent(inout) :: options
      !Local variables
      integer :: iphase


      call get_option("/material_phase[0]/multiphase_properties/immobile_fraction", &
      options%s_gc, default=0.0)
      call get_option("/material_phase[0]/multiphase_properties/relperm_type/Corey/relperm_max", &
      options%kr1_max, default=1.0)
      call get_option("/material_phase[0]/multiphase_properties/relperm_type/Corey/relperm_exponent", &
      options%kr1_exp, default=2.0)
      if (nphase>=2) then
          call get_option("/material_phase[1]/multiphase_properties/immobile_fraction", &
          options%s_or, default=0.0)
          call get_option("/material_phase[1]/multiphase_properties/relperm_type/Corey/relperm_max", &
          options%kr2_max, default=1.0)
          call get_option("/material_phase[1]/multiphase_properties/relperm_type/Corey/relperm_exponent", &
          options%kr2_exp, default=2.0)
      end if
      if (nphase>=3) then
          call get_option("/material_phase[2]/multiphase_properties/immobile_fraction", &
          options%s_gr, default=0.0)
          call get_option("/material_phase[2]/multiphase_properties/relperm_type/Corey/relperm_max", &
          options%kr3_max, default=1.0)
          call get_option("/material_phase[2]/multiphase_properties/relperm_type/Corey/relperm_exponent", &
          options%kr3_exp, default=2.0)
      end if
    end subroutine get_corey_options

    SUBROUTINE relperm_corey_epsilon( ABSP, visc_phase1, visc_phase2, INV_PERM, SAT, IPHASE,opt )
          !This subroutine add a small quantity to the corey function to avoid getting a relperm=0 that may give problems
          !when dividing it to obtain the sigma.
        IMPLICIT NONE
        REAL, intent( inout ) :: ABSP
        REAL, intent( in ) :: visc_phase1, visc_phase2, SAT, INV_PERM
        INTEGER, intent( in ) :: IPHASE
        type(corey_options), intent(in) :: opt
        ! Local variables...
        REAL :: KR, VISC, SATURATION, Krmax
        real, parameter :: epsilon = 1d-10
        !Kr_max should only multiply the wetting phase,
        !however as we do not know if it is phase 1 or 2, we let the decision to the user
        !and we multiply both phases by kr_max. By default kr_max= 1

        SATURATION = sat
        KR = get_relperm_Brooks_Corey(SATURATION, iphase, opt, visc_phase1, visc_phase2, visc, krmax)

        !Make sure that the relperm is between bounds
        KR = min(max(epsilon, KR),Krmax)!Lower value just to make sure we do not divide by zero.
        ABSP = INV_PERM * (VISC * max(1d-5,SATURATION)) / KR !The value 1d-5 is only used if the boundaries have values of saturation of zero.
        !Otherwise, the saturation should never be zero, since immobile fraction is always bigger than zero.

      RETURN
    END SUBROUTINE relperm_corey_epsilon

    real function get_relperm_Brooks_Corey(sat, iphase, opt, visc_phase1, visc_phase2, visc, krmax)
        !Calculates the Brooks-Corey relperm. If optional oldSAT
        !is introduced, the new relperm is calculated using:
        !Kr = Kr(SwNew) + dKr/dS(SwNew-SwOld)
        implicit none
        real, intent(inout) :: sat, krmax, visc
        real, intent(in) :: visc_phase1, visc_phase2
        integer, intent(in) :: iphase
        type(corey_options), intent(in) :: opt
        !Local variables
        real :: derivative, aux


        aux = 1. - opt%s_gc - opt%s_or

        IF( IPHASE == 1 ) THEN
            krmax = opt%kr1_max
            get_relperm_Brooks_Corey = krmax*( ( sat - opt%s_gc) /&
                 ( aux )) ** opt%kr1_exp
            Visc = visc_phase1
        else
            sat = 1.0 - SAT
            krmax = opt%kr2_max
            get_relperm_Brooks_Corey = krmax * ( ( sat - opt%s_or ) /&
                 ( aux )) ** opt%kr2_exp
            VISC = visc_phase2
        end if

    end function get_relperm_Brooks_Corey


    subroutine relperm_stone(ABSP, iphase, sat, opt, visc, INV_PERM)
        !This subroutine calculates the relative permeability for three phases
        !First phase has to be water, second oil and the third gas
        !We use Stone's model II adapted, and for the two phases we use the Corey model
        !Model explained in: Aziz, K. And Settari, T.:“Petroleum Reservoir Simulation” Applied Science Publishers, London, 30-38, 1979.
        implicit none
        real, intent(inout) :: ABSP
        real, intent(in) :: INV_PERM, visc
        real, intent(in), dimension(:) :: sat
        type(corey_options), intent(in) :: opt
        integer, intent(in) :: iphase
        !Local variables
        real, dimension(3) :: satura, relperm, KR, Krmax
        real :: auxVisc, auxKrmax, aux, Krow, Krog
        real, parameter :: epsilon = 1d-10

        !Prepare data
        Krmax(1) = opt%kr1_exp; Krmax(2) = opt%kr2_exp; Krmax(3) = opt%kr3_exp
        !We consider two models for two phase flow, water-oil and oil-gas
        if (iphase /= 3) then
            satura(1) = ( sat(1) - opt%s_gc) /( 1. - opt%s_gc - opt%s_or )!Water
            relperm(1) = Krmax(1)* satura(1) ** opt%kr1_exp!Water, Krw
        end if
        if (iphase /= 1) then
            satura(3) = ( sat(3) - opt%s_gr ) /(1. - opt%s_gr - opt%s_or)!Gas
            !For phase 1 and 3 (water and gas respectively) we can use the Brooks Corey model
            relperm(3) = Krmax(3)* satura(3) ** opt%kr3_exp!Gas, Krg
        end if
        !Oil relperm is obtained as a combination
        if (iphase ==2 ) then
            Krow = ( 1.0 - satura(1)) ** opt%kr2_exp!Oil, Krow
            Krog = ( 1.0 - satura(3)) ** opt%kr2_exp!Oil, Krog
            !For the second phase, oil, we need to recalculate the real value(Stone model 2)
            relperm(2) = Krmax(2)*( (Krow/Krmax(2) + relperm(1))*(Krog/Krmax(2) + relperm(3)) - (relperm(1) + relperm(3)) )
        end if
        !Make sure that the relperm is between bounds
        KR(iphase) = min(max(epsilon, relperm(iphase)),Krmax(iphase))!Lower value just to make sure we do not divide by zero.
        ABSP = INV_PERM * (VISC * max(1d-5,sat(iphase))) / KR(iphase) !The value 1d-5 is only used if the boundaries have values of saturation of zero.
        !Otherwise, the saturation should never be zero, since immobile fraction is always bigger than zero.

    end subroutine relperm_stone


   SUBROUTINE calculate_capillary_pressure( state, packed_state, Sat_in_FEM, StorageIndexes)

      ! CAPIL_PRES_OPT is the capillary pressure option for deciding what form it might take.
      ! CAPIL_PRES_COEF( NCAPIL_PRES_COEF, NPHASE, NPHASE ) are the coefficients
      ! Capillary pressure coefs have the dims CAPIL_PRES_COEF( NCAPIL_PRES_COEF, NPHASE,NPHASE )
      ! used to calculate the capillary pressure.

      IMPLICIT NONE
      type(state_type), dimension(:), intent(inout) :: state
      type(state_type), intent(inout) :: packed_state
      integer, dimension(:), intent(inout) :: StorageIndexes
      logical, intent(in) :: Sat_in_FEM
      ! Local Variables
      INTEGER :: nstates, ncomps, nphases, IPHASE, JPHASE, i, j, k, nphase, useful_phases
      real ::  S_OR, S_GC, auxO, auxW!c, a,
      character(len=OPTION_PATH_LEN):: option_path, phase_name, cap_path
      logical :: Pc_imbibition
      !Corey options
      type(corey_options) :: options
      !Working pointers
      real, dimension(:,:), pointer :: Satura, CapPressure
      type (scalar_field), pointer :: sfield
      real, pointer, dimension(:) :: C_regions, a_regions
      type(vector_field), pointer :: position

      !Get from packed_state
      if (Sat_in_FEM) then
          call get_var_from_packed_state(packed_state,FEPhaseVolumeFraction = Satura)
      else
          call get_var_from_packed_state(packed_state,PhaseVolumeFraction = Satura)
      end if
      call get_var_from_packed_state(packed_state,CapPressure = CapPressure)
      !Get corey options
      nphase =size(Satura,1)
      call get_corey_options(options, nphase)
      s_gc=options%s_gc
      s_or=options%s_or

      CapPressure = 0.
      Pc_imbibition = .false.

      DO IPHASE = 1, NPHASE


          if (have_option("/material_phase["//int2str(iphase-1)//&
            "]/multiphase_properties/capillary_pressure/type_Brooks_Corey") ) then
              option_path = "/material_phase["//int2str(iphase-1)//&
                "]/multiphase_properties/capillary_pressure/type_Brooks_Corey"

            !Get C
            cap_path = "/material_phase["//int2str(iphase-1)//&
                "]/multiphase_properties/capillary_pressure/type_Brooks_Corey/scalar_field::C/prescribed/value"
            call extract_scalar_from_diamond(state, c_regions, cap_path, "CapPe", StorageIndexes(32), iphase, nphase)
            !Get a
            cap_path = "/material_phase["//int2str(iphase-1)//&
                "]/multiphase_properties/capillary_pressure/type_Brooks_Corey/scalar_field::a/prescribed/value"
            call extract_scalar_from_diamond(state, a_regions, cap_path, "CapA", StorageIndexes(33), iphase, nphase)

            if (have_option("/material_phase["//int2str(iphase-1)//&
               "]/multiphase_properties/capillary_pressure/type_Brooks_Corey/Pc_imbibition") ) then
               Pc_imbibition = .true.
            end if

              if (IPHASE==1) then
                  auxW = S_GC
                  auxO = S_OR
              else
                  auxW = S_OR
                  auxO = S_GC
              end if
!              call get_option(trim(option_path)//"/c", c)
!              call get_option(trim(option_path)//"/a", a)
              !Apply Brooks-Corey model
              do jphase =1, nphase
                if (jphase /= iphase) then!Don't know how this will work for more than 2 phases
                  forall (k = 1:size(CapPressure,2))
!                      CapPressure( jphase, k ) = CapPressure( jphase, k ) + &
!                      Get_capPressure(satura(iphase,k), c, a, auxW, auxO)
                      CapPressure( jphase, k ) = CapPressure( jphase, k ) + &
                      Get_capPressure(satura(iphase,k), c_regions(k), a_regions(k), auxW, auxO, Pc_imbibition)
                  end forall
                end if
              end do
          end if
      END DO

      RETURN
    END SUBROUTINE calculate_capillary_pressure


    pure real function Get_capPressure(sat, Pe, a, Own_irr, Other_irr, Pc_imbibition)
        !This functions returns the capillary pressure for a certain input saturation
        !There is another function, its derivative in cv-adv-diff called Get_DevCapPressure
        Implicit none
        real, intent(in) :: sat, Pe, a, Own_irr, Other_irr
        logical, intent(in) :: Pc_imbibition
        !Local
        real, parameter :: tol = 1d-5

        if ( Pc_imbibition ) then
           Get_capPressure = &
           Pe * (1.0 - max(min((sat - Own_irr) / (1.0 - Own_irr - Other_irr), 1.0), tol)) ** (-a)
        else
           Get_capPressure = &
           Pe * max(min((sat - Own_irr) / (1.0 - Own_irr), 1.0), tol) ** (-a)
        end if

    end function Get_capPressure


    PURE real function Get_DevCapPressure(sat, Pe, a, Own_irr, Other_irr, Pc_imbibition)
        !This functions returns the derivative of the capillary pressure with the saturation
        Implicit none
        real, intent(in) :: sat, Pe, a, Own_irr, Other_irr
        logical, intent(in) :: Pc_imbibition
        !Local
        real, parameter :: tol = 1d-2
        real :: aux

        if ( Pc_imbibition ) then
           aux = (1.0 - Own_irr - Other_irr)
        else
           aux = (1.0 - Own_irr)
        end if
        Get_DevCapPressure = &
            -a * Pe * aux**a * max(min((sat - Own_irr), 1.0), tol) ** (-a-1)
    end function Get_DevCapPressure

    subroutine calculate_u_source(state, Density_FEMT, u_source)
      !u_source has to be initialized before calling this subroutine
      type(state_type), dimension(:), intent(in) :: state
      real, dimension(:,:), intent(inout) :: Density_FEMT
      real, dimension(:,:,:), intent(inout) :: u_source

      type(vector_field), pointer :: gravity_direction
      real, dimension(:), allocatable :: g
      real, dimension(:,:), allocatable :: aux_FEM_den
      logical :: have_gravity
      real :: gravity_magnitude, aux
      integer :: idim, iphase, nod, stat

      allocate(g(size(Density_FEMT,1)))

      call get_option( "/physical_parameters/gravity/magnitude", gravity_magnitude, stat )
      have_gravity = ( stat == 0 )


      if( have_gravity ) then
          !####################TEMPORARY####################
          allocate(aux_FEM_den(size(u_source,2),size(u_source,3)))
          do iphase = 1, size(u_source,2)
              call get_option(trim( '/material_phase[' // int2str( iphase - 1 ) // &
                  ']/equation_of_state/incompressible/linear/all_equal'), aux)
            aux_FEM_den( iphase, : ) = aux
          end do
          !#################################################


         gravity_direction => extract_vector_field( state( 1 ), 'GravityDirection' )
         g = node_val(gravity_direction, 1) * gravity_magnitude

         do nod = 1, size(u_source,3)
             do iphase = 1, size(u_source,2)
                do idim = 1, size(u_source,1)
                  u_source( idim, iphase, nod) = u_source( idim, iphase, nod) + &
                       aux_FEM_den( iphase, nod ) * g( idim )!aux_FEM_den is TEMPORARY
                   end do
            end do
         end do
          deallocate(aux_FEM_den)
      end if

      deallocate(g)
    end subroutine calculate_u_source

    subroutine calculate_u_source_cv(state, cv_nonods, ndim, nphase, den, u_source_cv)
      type(state_type), dimension(:), intent(in) :: state
      integer, intent(in) :: cv_nonods, ndim, nphase
      real, dimension(:,:), intent(in) :: den
      real, dimension(:,:,:), intent(inout) :: u_source_cv

      type(vector_field), pointer :: gravity_direction
      real, dimension(ndim) :: g
      logical :: have_gravity
      real :: gravity_magnitude
      integer :: idim, iphase, nod, stat

      call get_option( "/physical_parameters/gravity/magnitude", gravity_magnitude, stat )
      have_gravity = ( stat == 0 )

      if( have_gravity ) then
         gravity_direction => extract_vector_field( state( 1 ), 'GravityDirection' )
         g = node_val(gravity_direction, 1) * gravity_magnitude

         u_source_cv = 0.
         do nod = 1, cv_nonods
            do iphase = 1, nphase
               do idim = 1, ndim
                  !u_source_cv( nod + (idim-1)*cv_nonods + ndim*cv_nonods*(iphase-1) ) = &
                  !     den( iphase, nod ) * g( idim )
                  u_source_cv( idim, iphase, nod) = den( iphase, nod ) * g( idim )
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
      character( len = option_path_len ) :: option_path, option_path_python, buffer
      type(tensor_field), pointer :: t_field
      integer :: iphase, icomp, idim, stat, cv_nod, mat_nod, cv_nloc, ele

      type(scalar_field), pointer :: component, diffusivity
      integer, dimension(:), pointer :: st_nodes, c_nodes
      logical :: linearise_viscosity, python_diagnostic_field
      real, dimension( : ), allocatable :: component_tmp
      real, dimension( :, :, : ), allocatable :: mu_tmp

      character( len = python_func_len ) :: pycode
      real :: dt, current_time
      integer :: iloc



      if ( have_option( '/physical_parameters/mobility' ) .or. is_compact_overlapping ) then

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
         cv_snloc, nphase, ndim, nface, mat_nonods, cv_nonods, x_nloc, ncolele, cv_ele_type, &
         finele, colele, cv_ndgln, cv_sndgln, x_ndgln, mat_ndgln, material_absorption, &
         state, x_nonods )

      implicit none
      type( state_type ), intent( inout ) :: packed_state
      integer, intent(in) :: totele, stotel, cv_nloc, cv_snloc, nphase, ndim, nface, &
           mat_nonods, cv_nonods, x_nloc, ncolele, cv_ele_type, x_nonods
      integer, dimension( : ), intent( in ) :: finele
      integer, dimension( : ), intent( in ) :: colele
      integer, dimension( : ), intent( in ) :: cv_ndgln
      integer, dimension( : ), intent( in ) :: cv_sndgln
      integer, dimension( : ), intent( in ) :: x_ndgln
      integer, dimension( : ), intent( in ) :: mat_ndgln
      real, dimension( :, :, : ), intent( inout ) :: material_absorption
      type(state_type), dimension( : ) :: state  

      real, dimension(stotel * cv_snloc * nphase, ndim ), intent( inout ) :: suf_sig_diagten_bc

      ! local variables
      type(tensor_field), pointer :: viscosity_ph
      real, dimension(nphase) :: visc_phases
      integer :: iphase, ele, sele, cv_siloc, cv_snodi, cv_snodi_ipha, iface, s, e, &
           ele2, sele2, cv_iloc, idim, jdim, i, mat_nod, cv_nodi
      real :: mobility, satura_bc, visc_phase1, visc_phase2
      real, dimension( ndim, ndim ) :: inv_perm, sigma_out, sigma_in, mat, mat_inv
      integer, dimension( nface, totele) :: face_ele
      integer, dimension( mat_nonods*nphase ) :: idone
      integer, dimension(nface, cv_snloc ) :: cv_sloclist
      integer, dimension( cv_snloc ) :: cv_sloc2loc

      integer, dimension( :, :, :),  allocatable :: wic_u_bc, wic_vol_bc

      integer, parameter :: WIC_BC_DIRICHLET = 1
!!$ for the pressure b.c. and compact_overlapping method
!!$ make the material property change just inside the domain else on the surface only...
!      logical, parameter :: mat_change_inside = .true.
      logical, parameter :: mat_change_inside = .false.
! if mat_perm_bc_dg use the method that is used for DG between the elements.
      logical, parameter :: mat_perm_bc_dg = .true.
      logical :: is_corey
      type(corey_options) :: options

      type(tensor_field), pointer :: velocity, volfrac, perm
      type(tensor_field) :: velocity_BCs, volfrac_BCs

    !Get from packed_state
      volfrac=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
      velocity=>extract_tensor_field(packed_state,"PackedVelocity")
      perm=>extract_tensor_field(packed_state,"Permeability")

    allocate(wic_u_bc(velocity%dim(1),velocity%dim(2),&
         surface_element_count(velocity)))
    allocate(wic_vol_bc(volfrac%dim(1),volfrac%dim(2),&
         surface_element_count(volfrac)))

    call get_entire_boundary_condition(velocity,&
           ['weakdirichlet'],velocity_BCs,WIC_U_BC)
    call get_entire_boundary_condition(volfrac,&
           ['weakdirichlet'],volfrac_BCs,WIC_vol_BC)

      if( have_option( '/physical_parameters/mobility' ) )then
         call get_option( '/physical_parameters/mobility', mobility )
         visc_phases(1) = 1
         visc_phases(2) = mobility
      elseif( have_option( '/material_phase[1]/vector_field::Velocity/prognostic/tensor_field::Viscosity' // &
           '/prescribed/value::WholeMesh/isotropic' ) ) then
         DO IPHASE = 1, NPHASE!Get viscosity for all the phases
            viscosity_ph => extract_tensor_field( state( iphase ), 'Viscosity' )
            visc_phases(iphase) = viscosity_ph%val( 1, 1, 1 )!So far we only consider scalar viscosity
         end do
         mobility = visc_phases(2) / visc_phases(1)
      elseif( nphase == 1 ) then
         mobility = 0.
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

         is_corey=.false.

         if ( have_option("/material_phase["// int2str(iphase-1) //&
              "]/multiphase_properties/relperm_type/Corey") ) then
            is_corey=.true.
            call get_corey_options(options, nphase)
         end if

         do ele = 1, totele

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
                        cv_nodi =cv_sndgln(cv_snodi)
                        cv_snodi_ipha = cv_snodi + ( iphase - 1 ) * stotel * cv_snloc
                        mat_nod = mat_ndgln( (ele-1)*cv_nloc + cv_iloc  )
                        ! this is the boundary condition
                        ! of the first phase
                        satura_bc = volfrac_BCs%val(1,1,cv_snodi)

!                        sigma_out = 0.
                        do idim = 1, ndim
                           do jdim = 1, ndim
                               if (is_corey) then
                                   select case (nphase)
                                       case (1)!No relperm needed, we calculate directly the result
                                           sigma_out( idim, jdim ) = inv_perm( idim, jdim ) *&
                                           visc_phases(1) * min(1.0,max(1d-5,satura_bc))
                                       case (2)
                                           call relperm_corey_epsilon( sigma_out( idim, jdim ), visc_phases(1), visc_phases(2), &
                                           inv_perm( idim, jdim ), satura_bc, iphase,options)
                                       case (3)!For three phases we use the Stone model. !With predefined order: Water, oil, gas
                                           call relperm_stone(sigma_out( idim, jdim ), iphase,&
                                           volfrac_BCs%val(1,:,cv_snodi), options, visc_phases(iphase), inv_perm( idim, jdim ))
                                       case default
                                           FLAbort("No relative permeability function implemented for more than 3 phases")
                                   end select
                               end if
                           end do
                        end do

                        if(mat_perm_bc_dg) then
! if mat_perm_bc_dg use the method that is used for DG between the elements.
                           sigma_in=0.0
                           sigma_in = material_absorption( mat_nod, s : e, s : e )
                           mat = sigma_out  +  matmul(  sigma_in,  matmul( inverse( sigma_out ), sigma_in ) )
                           mat_inv = matmul( inverse( sigma_in+sigma_out ), mat )
                           suf_sig_diagten_bc( cv_snodi_ipha, 1 : ndim ) = (/ (mat_inv(i, i), i = 1, ndim) /)
                     !!      suf_sig_diagten_bc( cv_snodi_ipha, 1 : ndim ) = 1.
                        else
                           mat = matmul( sigma_out, inverse( material_absorption( mat_nod, s : e, s : e ) ) )
                           mat_inv = inverse( mat )
                           suf_sig_diagten_bc( cv_snodi_ipha, 1 : ndim ) = (/ (mat_inv(i, i), i = 1, ndim) /)
                        endif

                        if( mat_change_inside ) then
                           suf_sig_diagten_bc( cv_snodi_ipha, 1 : ndim ) = 1.

                           if( idone( mat_nod+(iphase-1)*mat_nonods ) == 0 ) then
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
      integer :: apply_dim, idim, jdim, ipha_idim, iphase, ele, cv_iloc, ij, imat
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


    subroutine extract_scalar_from_diamond(state, field_values, path, StorName, indx, iphase, nphase)
    !Gets a scalar field directly from Diamond
    !Path have to end in /prescribed/value
    !Indx is for the cashing
    !If no phases, then pass iphase = nphase = 1
    !NOTE: This was initially done for capillary pressure with regions
        implicit none
        type(state_type), dimension(:), intent(inout) :: state
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
              if (has_scalar_field(state(1), StorName)) then
                  !If we are recalculating due to a mesh modification then
                  !we return to the original situation
                  call remove_scalar_field(state(1), StorName)
              end if


            !By default I use the Pressure mesh (Number 1)
            Sfield => extract_scalar_field(state(1),1)
            position => get_external_coordinate_field(state(1), Sfield%mesh)
            fl_mesh => extract_mesh( state(1), "CoordinateMesh" )
            Auxmesh = fl_mesh
            !The number of nodes I want does not coincide
            Auxmesh%nodes = size(Sfield%val,1) * nphase
            call allocate (targ_Store, Auxmesh, StorName)

!            call allocate(targ_Store, Sfield%mesh)
            call initialise_field_over_regions(targ_Store, path, position)
            !Now we insert them in state and store the indexes
            call insert(state(1), targ_Store, StorName)
            call deallocate (targ_Store)
            indx = size(state(1)%scalar_fields)
          end if
          !Get the data
          siz = size(state(1)%scalar_fields(abs(indx))%ptr%val(:),1)/nphase
          field_values => state(1)%scalar_fields(abs(indx))%ptr%val((iphase-1)*siz + 1: siz * iphase )


    end subroutine extract_scalar_from_diamond

    subroutine Set_Saturation_between_bounds(packed_state)
        !This subroutines eliminates the oscillations in the saturation that are bigger than a
        !certain tolerance
        Implicit none
        !Global variables
        type( state_type ), intent(inout) :: packed_state
        !Local variables
        integer :: iphase, jphase, nphase
        real :: top_limit
        real, dimension(:,:), pointer :: satura
        real, dimension(:), allocatable :: immobile_fractions
        real, parameter :: tolerance = 0.0

        call get_var_from_packed_state(packed_state, PhaseVolumeFraction = satura)

        nphase = size(satura,1)
        allocate(immobile_fractions(nphase))

        !Get immobile fractions for all the phases
        do iphase = 1, nphase
          call get_option("/material_phase["//int2str(iphase-1)//"]/multiphase_properties/immobile_fraction", &
               immobile_fractions(iphase), default=0.0)
        end do
        !Set saturation to be between bounds
        do iphase = 1, nphase
            top_limit = 0.
            !Get the upper limit of saturation for iphase
            do jphase = 1, nphase
                if (jphase/=iphase) top_limit = top_limit + immobile_fractions(jphase)
            end do
            top_limit = 1.0 + tolerance - top_limit
            !Limit the saturation, allowing small oscillations
            satura(iphase,:) =  min(max(immobile_fractions(iphase)-tolerance, satura(iphase,:)),top_limit)
        end do
        deallocate(immobile_fractions)
    end subroutine Set_Saturation_between_bounds




    subroutine boiling( states, packed_state, cv_nonods, mat_nonods, nphase, ndim, &
         ScalarField_Source, velocity_absorption, temperature_source, temperature_absorption )
      implicit none

         type( state_type ), dimension(:), intent( inout ) :: states
         type( state_type ), intent( in ) :: packed_state
      integer, intent( in ) :: cv_nonods, mat_nonods, nphase, ndim
      real, dimension( :, : ), intent( inout ) :: ScalarField_Source, temperature_source
      real, dimension( :, :, : ), intent( inout ) :: velocity_absorption, temperature_absorption

      type( tensor_field ), pointer :: temperature
      real, dimension( :, : ), allocatable :: A
      real, dimension( : ), allocatable :: S_lg_l, S_lg_g, S_ls_l, S_gs_g, &
           T_sat, Svap_l, Svap_g, Gamma_l, Gamma_g, h_l, h_g, &
           St_gl, St_sl, St_sg
      integer :: iphase, jphase, idim
      real, parameter :: Le0=2375.7e3, Cp_l = 4200.0, Cp_g = 1996.0

      ewrite(3,*) 'inside boiling routine'

      ScalarField_Source=0.0 ; temperature_source=0.0
      velocity_absorption=0.0 ; temperature_absorption=0.0

      allocate( S_lg_l(mat_nonods), S_lg_g(mat_nonods), S_ls_l(mat_nonods), S_gs_g(mat_nonods), A(nphase,mat_nonods) )
      allocate( T_sat(cv_nonods), Svap_l(cv_nonods), Svap_g(cv_nonods), &
           &    Gamma_l(cv_nonods), Gamma_g(cv_nonods), h_l(cv_nonods), h_g(cv_nonods), &
           &    St_gl(cv_nonods), St_sl(cv_nonods), St_sg(cv_nonods) )


      call calculate_boiling_drag( states, packed_state, S_lg_l, S_lg_g, S_ls_l, S_gs_g, A )

      ! Momentum absorption

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

      temperature => extract_tensor_field( packed_state, "PackedTemperature" )

      ! Temperature absorption

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

      iphase=1
      temperature_source( iphase, : ) = Svap_l*T_sat + Gamma_l*h_l + Gamma_l*Le0

      iphase=2
      temperature_source( iphase, : ) = Svap_g*T_sat + Gamma_g*h_g + 1.0e+6 * temperature%val(1,2,:)


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




    subroutine calculate_boiling_drag( states, packed_state, S_lg_l, S_lg_g, S_ls_l, S_gs_g, A )
      implicit none

      type( state_type ), dimension( : ), intent( inout ) :: states
      type( state_type ), intent( in ) :: packed_state
      real, dimension( : ), intent( inout ) :: S_lg_l, S_lg_g, S_ls_l, S_gs_g
      real, dimension( :, : ), intent( inout ) :: A

      type( scalar_field ), pointer :: pressure, dummy
      type( tensor_field ), pointer :: density, velocity, volume_fraction

      integer, dimension( : ), pointer :: mat_ndgln, cv_ndgln, u_ndgln

      real :: rho_l, rho_g, u_l, u_g, u_s, u_gs, u_ls, u_gl, a_l, a_g, a_s, &
           a_sg, a_gs, a_sl, a_ls, a_lg, a_gl, d_b, Re_gl, Re_lg, CD

      real, dimension( : ), allocatable :: ul, ug, us

      integer :: ele, totele, mat_iloc, mat_nloc, u_nloc, mat_inod, cv_inod, &
           u_inod, u_inod1, u_inod2, ndim, stat

      real, parameter :: &
           mu_l = 3.0e-4, mu_g = 1.0e-5, &
           d_p = 0.005

      pressure => extract_scalar_field( packed_state, "CVPressure" )
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

            S_gs_g(mat_inod) = 150.0 * (a_gs*mu_g) / (a_sg*d_p**2*(a_g+a_s)) + 1.75 * (rho_g*u_gs) / (d_p*(a_g+a_s))
            S_ls_l(mat_inod) = 150.0 * (a_ls*mu_l) / (a_sl*d_p**2*(a_l+a_s)) + 1.75 * (rho_l*u_ls) / (d_p*(a_l+a_s))

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




    subroutine calculate_boiling_variables( states, packed_state, ndim, nphase, &
         T_sat, Svap_l, Svap_g, Gamma_l, Gamma_g, h_l, h_g, St_gl, St_sl, St_sg )
      implicit none

      type( state_type ), dimension( : ), intent( inout ) :: states
      type( state_type ), intent( in ) :: packed_state
      integer, intent( in ) :: ndim, nphase
      real, dimension( : ), intent( inout ) :: T_sat, Svap_l, Svap_g, Gamma_l, Gamma_g, h_l, h_g, St_gl, St_sl, St_sg

      type( scalar_field ), pointer :: pressure, dummy
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

      pressure => extract_scalar_field( packed_state, "CVPressure" )
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

            p = pressure%val(cv_inod)

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




  end module multiphase_EOS
