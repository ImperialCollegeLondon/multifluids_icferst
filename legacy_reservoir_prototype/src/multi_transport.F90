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

module multi_transport

    use fldebug
    use multi_data_types
    use multiphase_1D_engine

    implicit none
    private

    public :: solve_transport

    interface solve_transport
        module procedure solve_transport_scalars
    end interface


contains

    subroutine solve_transport_scalars()
       implicit none

        type(multi_transport_scalar) :: test

        !!$ Solve advection of the scalars.   'Temperature':

!!$ Fields...
!!-
!!!!!!!!!DO NOT REMOVE THE CODE COMMENTED BELOW!!!!
!        new_ntsol_loop = .false.
!
!if ( new_ntsol_loop  ) then
!
!        call get_ntsol( ntsol )
!        call initialise_field_lists_from_options( state, ntsol )
!
!
!        call set_nu_to_u( packed_state )
!        !call calculate_diffusivity( state, Mdim, ndgln, ScalarAdvectionField_Diffusion )
!        velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
!        density_field=>extract_tensor_field(packed_state,"PackedDensity",stat)
!        saturation_field=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
!
!
!
!
!
!        do it = 1, ntsol
!
!           call get_option( trim( field_optionpath_list( it ) ) // &
!                '/prognostic/equation[0]/name', &
!                option_buffer, default = "UnknownEquationType" )
!           select case( trim( option_buffer ) )
!           case ( "AdvectionDiffusion", "InternalEnergy" )
!              use_advdif = .true.
!           case default
!              use_advdif = .false.
!           end select
!
!           !use_advdif=.true.
!
!
!           if ( use_advdif ) then
!
!              ! figure out if scalar field is mutli-phase
!              multiphase_scalar = .false.
!              do it2 = it+1, ntsol
!                 if ( field_name_list( it ) == field_name_list( it2 ) ) then
!                    multiphase_scalar = .true.
!                 end if
!              end do
!
!              tmp_name = "Packed" //field_name_list( it )
!              nphase_scalar = 1
!              if ( multiphase_scalar ) then
!                 nphase_scalar = Mdims%nphase
!                 tmp_name = "Packed" // field_name_list( it )
!              end if
!              tracer_field => extract_tensor_field( packed_state, trim( tmp_name ) )
!
!
!              if (field_name_list( it)== 'PhaseVolumeFraction' .or.  field_name_list( it)== 'ComponentMassFractionPhase[0]') then
!                    cycle
!              elseif (multiphase_scalar) then
!
!                    call INTENERGE_ASSEM_SOLVE( state, packed_state, &
!                        Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,upwnd,&
!                        tracer_field,velocity_field,density_field, dt, &
!                        suf_sig_diagten_bc, &
!                        Porosity_field%val, &
!                        !!$
!                        0, igot_theta_flux, &
!                        Mdisopt%t_get_theta_flux, Mdisopt%t_use_theta_flux, &
!                        THETA_GDIFF, IDs_ndgln, &
!                        option_path = '/material_phase[0]/scalar_field::Temperature', &
!                        thermal = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/equation::InternalEnergy'),&
!                        saturation=saturation_field )
!                    call Calculate_All_Rhos( state, packed_state, Mdims )
!
!                    exit
!              else
!
!                    call INTENERGE_ASSEM_SOLVE( state, packed_state, &
!                        Mdims, CV_GIdims, CV_funs, Mspars, ndgln, Mdisopt, Mmat,upwnd,&
!                        tracer_field,velocity_field,density_field, dt, &
!                        suf_sig_diagten_bc,  Porosity_field%val, &
!                        0, igot_theta_flux, &
!                        Mdisopt%t_get_theta_flux, Mdisopt%t_use_theta_flux, &
!                        THETA_GDIFF, IDs_ndgln, &
!                        option_path = '/material_phase[0]/scalar_field::Temperature', &
!                        thermal = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/equation::InternalEnergy'),&
!                        saturation=saturation_field)
!                    call Calculate_All_Rhos( state, packed_state, Mdims )
!               end if
!
!
!           end if
!
!        end do
!
!end if

    end subroutine solve_transport_scalars


end module multi_transport
