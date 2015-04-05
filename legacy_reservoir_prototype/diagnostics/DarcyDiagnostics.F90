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

module darcy_diagnostics

  use boundary_conditions
  use diagnostic_source_fields
  use field_derivatives
  use field_options
  use fields
  use fldebug
  use global_parameters, only : OPTION_PATH_LEN
  use spud
  use state_fields_module
  use state_module
  use multiphase_EOS
  
  implicit none
  
  private
  
  public :: calculate_darcy_velocity
  
contains

  subroutine calculate_darcy_velocity(state, v_field)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    
    character(len = OPTION_PATH_LEN) :: path
    integer :: stat
    real :: scale_factor
    type(scalar_field), pointer :: saturation
    type(vector_field), pointer :: force_density
    type(vector_field) :: rhs

    integer :: ele

    real, dimension(ele_ngi(v_field,1)) :: s,detwei


    call zero(v_field)

    call allocate(rhs,mesh=v_field%mesh,dim=v_field%dim,name="RHS")
    
    force_density => extract_vector_field(state, "Velocity")
    saturation => extract_scalar_field(state, "PhaseVolumeFraction")
    path = trim(complete_field_path(v_field%option_path)) // "/algorithm"

    do ele=1,element_count(v_field)

       S=ele_val_at_quad(saturation,ele)
!       call addto(rhs,ele_nodes(v_field,ele),&
!            shape_vector_rhs(v_field%mesh%shape,&
!            force_density,&
!            mobility(state,S,ele_region_id(v_field,ele))*detwei))

    end do
    
  
  end subroutine calculate_darcy_velocity

end module darcy_diagnostics
