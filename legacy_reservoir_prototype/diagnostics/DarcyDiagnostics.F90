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
  use elements
  use transform_elements
  use fields
  use fldebug
  use global_parameters, only : OPTION_PATH_LEN
  use spud
  use state_fields_module
  use state_module
  use solvers
  use initialise_fields_module
  use multiphase_EOS
  use fetools
  
  implicit none
  
  private

  public :: calculate_darcy_velocity
  
contains

  subroutine calculate_darcy_velocity(states,state_index, v_field)
    type(state_type), dimension(:), intent(inout) :: states
    integer, intent(in) :: state_index
    type(vector_field), intent(inout) :: v_field
    
    character(len = OPTION_PATH_LEN) :: path

    type(csr_matrix), pointer :: mass
    type(scalar_field), pointer :: saturation, scalar_permeability, immobile_fraction
    type(vector_field), pointer :: force_density, X, vector_permeability
    type(tensor_field), pointer :: tensor_permeability
    type(vector_field) :: rhs
    type(element_type), pointer :: v_shape

    integer :: ele, phase, dim, stat

    real :: mu, mobility, alpha,beta

    real, dimension(ele_ngi(v_field,1)) :: s,s_i,k,detwei
    real, dimension(mesh_dim(v_field),mesh_dim(v_field),1) :: k_a
    real, dimension(ele_loc(v_field,1), ele_ngi(v_field,1), mesh_dim(v_field)) :: dv_t 

    call zero(v_field)
    call allocate(rhs,mesh=v_field%mesh,dim=v_field%dim,name="RHS")
    call zero(rhs)
    
    !! THIS IS REALLY UGLY, BUT NECESSARY UNTIL THE MOBILITY OPTION IS REMOVED
    if (have_option('/physical_parameters/mobility')) then
       call get_option('/material_phase[0]/vector_field::Velocity/prognostic/tensor_field::Viscosity/prescribed/value::WholeMesh/isotropic/constant',mu)
       if (state_index==2) then
          call get_option('/physical_parameters/mobility',mobility)
          mu=mu*mobility
       end if
    else
       call get_option('/material_phase['//int2str(state_index-1)//']/vector_field::Velocity/prognostic/tensor_field::Viscosity/prescribed/value::WholeMesh/isotropic/constant',mu)
    end if

    call get_option('trim(states(state_index)%option_path)'//&
         '/multiphase_properties/relperm_max',beta,default=1.0)
    call get_option('trim(states(state_index)%option_path)'//&
         '/multiphase_properties/relperm_exponent',alpha,default=2.0)

    X => extract_vector_field(states(1), "Coordinate")

    if (.not. has_scalar_field(states(state_index), "ImmobileFraction")) then
       call allocate_and_insert(trim(states(state_index)%option_path)//&
            "/multiphase_properties/immobile_fraction/scalar_field::value",&
            states(state_index),X,field_name="ImmobileFraction")
    end if

    !! We really need to sort out the name mismatch here
    force_density => extract_vector_field(states(state_index), "Velocity")
    saturation => extract_scalar_field(states(state_index),&
         "PhaseVolumeFraction")


    !! Slightly messy code to get the permeability
    scalar_permeability => extract_scalar_field(states(1),&
         "Permeability",stat)
    if (stat/=0) then
       nullify(scalar_permeability)
       vector_permeability => extract_vector_field(states(1),&
         "Permeability",stat)
       if (stat/=0) then
          nullify(vector_permeability)
          tensor_permeability => extract_tensor_field(states(1),&
         "Permeability")
       end if
    end if
       
    immobile_fraction => extract_scalar_field(states(state_index),&
         "ImmobileFraction")

    path = trim(complete_field_path(v_field%option_path)) // "/algorithm"

    do ele=1,element_count(v_field)

       S=ele_val_at_quad(Saturation,ele)
       S_I=ele_val_at_quad(immobile_fraction,ele)

       k_a=0.0
       if (associated(scalar_permeability)) then
          do dim=1,mesh_dim(v_field)
             k_a(dim,dim,:)=ele_val(scalar_permeability,ele)
          end do
       else if (associated(vector_permeability)) then
          do dim=1,mesh_dim(v_field)
             k_a(dim,dim,:)=ele_val(vector_permeability,ele,dim)
          end do
       else
          k_a(:,:,:)=ele_val(tensor_permeability,ele)
       end if


       !! assuming Corey relative permeability here.
       !! This needs fixing for the general case 
       k=relperm(S,S_I,alpha,beta)


       !!! build the right hand side of the projection problem here
       v_shape => ele_shape(v_field, ele)

       call transform_to_physical(X, ele, v_shape, &
            dshape=dv_t, detwei = detwei)

       call addto(rhs,ele_nodes(v_field,ele),&
            shape_vector_rhs(v_shape,&
               matmul(k_a(:,:,1),ele_val_at_quad(force_density,ele)),&
               k*detwei/mu))

    end do


    !!! should probably get short circuited for the discontinuous case.
    mass => get_mass_matrix(states, v_field%mesh)
    call petsc_solve(v_field,mass,rhs,option_path=path)

    call deallocate(rhs)


  contains

    function relperm(Sat,ImSat,a,b)
      real, dimension(:) :: Sat,ImSat
      real :: a,b
      
      real, dimension(size(sat)) :: relperm
      
      relperm=b*(Sat-ImSat)**a
    end function relperm

    subroutine allocate_and_insert(option_path, state, positions, &
         parent_mesh, parent_name, field_name)
      
      character(len=*), intent(in) :: option_path
      type(state_type), intent(inout) :: state
      type(vector_field), intent(in) :: positions
      character(len=*), intent(in), optional :: parent_mesh
      character(len=*), intent(in), optional :: parent_name
      character(len=*), optional, intent(in):: field_name
      
      logical :: is_prognostic, is_prescribed, is_diagnostic
      ! paths for options and child fields
      character(len=OPTION_PATH_LEN) :: path
      ! Strings for names
      character(len=OPTION_PATH_LEN) :: lfield_name, mesh_name
      type(scalar_field) :: field
      type(mesh_type), pointer :: mesh
      logical :: is_constant

      ! Save option_path
      path=trim(option_path)

      if (present(field_name)) then
         lfield_name=field_name
      else
         call get_option(trim(path)//"/name", lfield_name)
      end if
        
      if(present(parent_name)) then
         lfield_name=trim(parent_name)//trim(lfield_name)
      end if

      ! Find out what kind of field we have
      is_prognostic=have_option(trim(path)//"/prognostic")
      is_prescribed=have_option(trim(path)//"/prescribed")
      is_diagnostic=have_option(trim(path)//"/diagnostic")

      
      if (is_prognostic) then

         path=trim(path)//"/prognostic"

      else if(is_prescribed) then
           
         path=trim(path)//"/prescribed"
           
      else if(is_diagnostic) then
       
         path=trim(path)//"/diagnostic"

      end if

      ! Get mesh
      if(present(parent_mesh).and.&
           .not.have_option(trim(path)//"/mesh[0]/name")) then
         mesh => extract_mesh(state, trim(parent_mesh))
         mesh_name=parent_mesh
      else
         call get_option(trim(path)//"/mesh[0]/name", mesh_name)
         mesh => extract_mesh(state, trim(mesh_name))
      end if

      ! Allocate field
      call allocate(field, mesh, name=trim(lfield_name))
      call initialise_field_over_regions(field,field%option_path,positions)

      ewrite(2,*) trim(lfield_name), " is on mesh ", trim(mesh%name)

      ! Set field%option_path
      field%option_path=trim(option_path)

      ! Finally! Insert field into state!
      call insert(state, field, field%name)
      call deallocate(field)

    end subroutine allocate_and_insert
  
  end subroutine calculate_darcy_velocity


end module darcy_diagnostics
