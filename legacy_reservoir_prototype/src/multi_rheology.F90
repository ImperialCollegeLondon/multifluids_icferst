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

  module multiphase_rheology

    use fldebug
    use global_parameters, only: OPTION_PATH_LEN, PYTHON_FUNC_LEN, PI
    use spud
    use fields
    use field_options
    use fefields
    use state_module
    use initialise_fields_module
    use populate_state_module

    use field_derivatives, only : strain_rate
    
    implicit none
    
    type power_law_rheology_type
       real :: yield_stress
       real :: exponent
       real :: consistency_index
       real :: tolerance
       real :: lower_bound,upper_bound
    end type power_law_rheology_type

    type rheology_type
       integer :: rheology_model
       type(power_law_rheology_type), pointer :: power_law_ptr=>null()
    end type rheology_type


    private
  
    public :: rheology_type, initialize_rheologies, finalize_rheologies, calculate_rheologies



    contains

      subroutine nullify_rheology(rheology)
        type(rheology_type), intent(inout) :: rheology

        nullify(rheology%power_law_ptr)
        
      end subroutine nullify_rheology
      
      function initialize_rheology(option_path) result(rheology)

        character (len=*), intent(in) :: option_path
        type(rheology_type) :: rheology

        call nullify_rheology(rheology)

        if (have_option(option_path//'/rheology/newtonian')) then
           rheology%rheology_model=1
        else if (have_option(option_path//'/rheology/non_newtonian::ShearDependent')) then
           rheology%rheology_model=2
           allocate(rheology%power_law_ptr)
           call get_option(option_path//&
                '/rheology/non_newtonian::ShearDependent/yield_stress',&
                rheology%power_law_ptr%yield_stress, default=0.0)
           call get_option(option_path//&
                '/rheology/non_newtonian::ShearDependent/exponent',&
                rheology%power_law_ptr%exponent, default=1.0)
           call get_option(option_path//&
                '/rheology/non_newtonian::ShearDependent/consistency_index',&
                rheology%power_law_ptr%consistency_index, default=1.0)
           call get_option(option_path//&
                '/rheology/non_newtonian::ShearDependent/tolerance',&
                rheology%power_law_ptr%tolerance, default=1.0e-10)
           call get_option(option_path//&
                '/rheology/non_newtonian::ShearDependent/lower_bound',&
                rheology%power_law_ptr%lower_bound, default=0.0)
           call get_option(option_path//&
                '/rheology/non_newtonian::ShearDependent/upper_bound',&
                rheology%power_law_ptr%upper_bound, default=1.0e10)
           
        else
           rheology%rheology_model=0
        end if

      end function initialize_rheology

      subroutine initialize_rheologies(state,rheology)
        type(state_type), dimension(:) :: state
        type(rheology_type), dimension(:) :: rheology

        integer :: i

        type(tensor_field), pointer :: viscosity
        type(vector_field), pointer :: position

        


        do i=1,size(state)
           rheology(i)=initialize_rheology(trim(state(i)%option_path))
           if (rheology(i)%rheology_model==1) then
              call allocate_and_insert_tensor_field(trim(state(i)%option_path)&
                   //'/rheology/newtonian/tensor_field',&
                   state(i))
              viscosity=>extract_tensor_field(state(i),"Viscosity")
              position => get_external_coordinate_field(state(i), viscosity%mesh)
              call initialise_field_over_regions(viscosity, &
                   trim(viscosity%option_path)//'/prescribed/value', &
                position)

           else if (rheology(i)%rheology_model==2) then
              call allocate_and_insert_tensor_field(trim(state(i)%option_path)&
                   //'/rheology/non_newtonian[0]/tensor_field',&
                   state(i))
           end if

        end do
           
      end subroutine initialize_rheologies

      subroutine finalize_rheology(rheology)
        type(rheology_type) :: rheology

        if (associated(rheology%power_law_ptr)) &
             deallocate(rheology%power_law_ptr)
        call nullify_rheology(rheology)

      endsubroutine finalize_rheology

      subroutine finalize_rheologies(rheology)
        type(rheology_type), dimension(:) :: rheology

        integer :: i

        do i=1,size(rheology)
           call finalize_rheology(rheology(i))
        end do
      end subroutine finalize_rheologies


      subroutine calculate_power_law_rheology(positions,state,viscosity,rheology,&
           continuous_velocity_mesh)
        type(tensor_field)  :: viscosity
        type(power_law_rheology_type) :: rheology
        type(state_type) :: state
        type(vector_field) :: positions
        type(mesh_type) :: continuous_velocity_mesh

        type(vector_field), pointer :: velocity
        type(tensor_field) :: strain_rate_tensor
        type(vector_field) :: continuous_velocity
        type(scalar_field) :: strain_rate_invariant
        type(scalar_field) :: scalar_viscosity, projected_scalar_viscosity

        integer :: i, j, dim


        velocity=>extract_vector_field(state,"Velocity")
        dim=mesh_dim(velocity)
        call allocate(continuous_velocity,mesh=continuous_velocity_mesh,&
             name="ContinuousVelocity",dim=dim)
        call allocate(strain_rate_tensor,velocity%mesh,name="StrainRateTensor",&
             dim=[dim,dim])
        call allocate(strain_rate_invariant,velocity%mesh,"StrainRateInvariant")
        call allocate(scalar_viscosity,velocitY%mesh,"ScalarViscosity")
        call allocate(projected_scalar_viscosity,viscosity%mesh,&
             "ProjectedScalarViscosity")
        
        call project_field(velocity, continuous_velocity, positions)
        call strain_rate(continuous_velocity,positions,strain_rate_tensor)
        call tensor_second_invariant(strain_rate_tensor, strain_rate_invariant)
        
        call bound(strain_rate_invariant,rheology%tolerance,huge(rheology%tolerance))
        scalar_viscosity%val=strain_rate_invariant%val**&
             (rheology%exponent-1.0d0)

        scalar_viscosity%val=&
             rheology%consistency_index*scalar_viscosity%val&
             +rheology%yield_stress/strain_rate_invariant%val
        
        call bound(scalar_viscosity,&
             rheology%lower_bound,&
             rheology%upper_bound)


        call project_field(scalar_viscosity,&
             projected_scalar_viscosity,positions)


        call zero(viscosity)

        do i=1,dim
           do j=1,dim
              call set(viscosity,i,j,projected_scalar_viscosity)
           end do
        end do
        
        call deallocate(strain_rate_tensor)
        call deallocate(strain_rate_invariant)
        call deallocate(scalar_viscosity)
        call deallocate(projected_scalar_viscosity)
        call deallocate(continuous_velocity)


      end subroutine calculate_power_law_rheology

      subroutine calculate_rheologies(state,rheology)
        type(state_type), dimension(:) :: state
        type(rheology_type), dimension(:) :: rheology

        integer :: i, j

        type(scalar_field), pointer :: cvf
        type(vector_field), pointer :: positions
        type(tensor_field), pointer :: viscosity, viscosity2
        
        type(mesh_type), pointer :: continuous_velocity_mesh

        positions=>extract_vector_field(state(1),"Coordinate")
        continuous_velocity_mesh=>extract_mesh(state(1),"VelocityMesh_Continuous")

        do i=1,size(state)
           if (rheology(i)%rheology_model==2) then
              viscosity=>extract_tensor_field(state(i),"Viscosity")
              call calculate_power_law_rheology(positions,state(i),&
                   viscosity,rheology(i)%power_law_ptr,continuous_velocity_mesh)
           end if
        end do

        do i=1,size(state)
           if ( have_option(trim(state(i)%option_path)//&
                '/vector_field::Velocity/prognostic/tensor_field::Viscosity/diagnostic/algorithm::rheological')) then
              viscosity=>extract_tensor_field(state(i),"Viscosity")
              call zero(viscosity)
              do j=1,size(state)
                 if (rheology(j)%rheology_model>0) then
                    cvf=>extract_scalar_field(state(j),"ComponentMassFractionPhase"//int2str(i))
                    viscosity2=>extract_tensor_field(state(j),"Viscosity")
                    call addto(viscosity,viscosity2,sscale=cvf)
                 end if
              end do
           end if
        end do
                 
           

      end subroutine calculate_rheologies

    end module multiphase_rheology
